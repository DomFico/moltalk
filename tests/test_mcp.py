import asyncio
import os
import sys
import time
import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from starlette.testclient import TestClient
from moltalk import limits
from moltalk.server import mcp, resolve_name, WIDGET_URI, WIDGET_MIME

HEADERS = {'Accept': 'application/json, text/event-stream'}
INIT = {'jsonrpc': '2.0', 'id': 1, 'method': 'initialize', 'params': {
    'protocolVersion': '2025-03-26', 'capabilities': {}, 'clientInfo': {'name': 'test', 'version': '1'}}}


def stdio_session(test):
    async def run():
        params = StdioServerParameters(command=sys.executable, args=['-m', 'moltalk.server'])
        async with stdio_client(params) as (read, write):
            async with ClientSession(read, write) as client:
                await client.initialize()
                await test(client)
    asyncio.run(run())


def test_stdio_client():
    async def test(client):
        tools = await client.list_tools()
        assert len(tools.tools) == 6
        result = await client.call_tool('analyze_molecule', {'smiles': 'CCO'})
        assert not result.isError
        assert result.structuredContent['formula'] == 'C2H6O'
        invalid = await client.call_tool('analyze_molecule', {'smiles': 'bad'})
        assert invalid.isError
        assert 'No structure was drawn or guessed' in invalid.content[0].text
    stdio_session(test)


def test_ui_resource_and_tool_metadata():
    async def test(client):
        resources = (await client.list_resources()).resources
        widget = next(r for r in resources if str(r.uri) == WIDGET_URI)
        assert widget.mimeType == WIDGET_MIME == 'text/html;profile=mcp-app'
        assert widget.meta['ui']['csp'] == {'connectDomains': [], 'resourceDomains': []}
        contents = (await client.read_resource(WIDGET_URI)).contents[0]
        assert contents.mimeType == WIDGET_MIME
        assert 'ui/initialize' in contents.text and 'ui/notifications/tool-result' in contents.text
        tools = {t.name: t for t in (await client.list_tools()).tools}
        for name in ('draw_molecule', 'enumerate_stereoisomers'):
            assert tools[name].meta['ui']['resourceUri'] == WIDGET_URI
            assert tools[name].meta['openai/outputTemplate'] == WIDGET_URI
        for name in ('analyze_molecule', 'find_substructure', 'resolve_name'):
            assert not (tools[name].meta or {}).get('ui')
        assert all(t.annotations.readOnlyHint for t in tools.values())
        assert tools['conformer_3d'].meta['ui']['visibility'] == ['app']  # viewer-only, hidden from the model
        assert tools['resolve_name'].annotations.openWorldHint
    stdio_session(test)


def test_draw_result_splits_model_and_widget_data():
    async def test(client):
        result = await client.call_tool('draw_molecule', {'smiles': 'F[C@H]1C[C@@H](C)CCC1',
                                                          'label': '(1R,3S)-1-fluoro-3-methylcyclohexane'})
        assert not result.isError
        data = result.structuredContent
        assert 'svg' not in data  # bulky SVG goes to the widget only
        assert result.meta['svg'].startswith('<?xml') and '<svg' in result.meta['svg']
        assert {c['atom_index']: c['cip'] for c in data['analysis']['stereocenters']} == {1: 'R', 3: 'S'}
        bonds = data['depicted_stereo_bonds']
        assert {(b['from_atom'], b['to_atom'], b['style'], b['stereocenter_cip']) for b in bonds} == {
            (1, 0, 'wedge', 'R'), (3, 4, 'wedge', 'S')}
        text = result.content[0].text  # with a PubChem name the centres are given as locants too
        assert 'atom 1 R, atom 3 S' in text or 'C1 (atom 1) R, C3 (atom 3) S' in text
        with_svg = await client.call_tool('draw_molecule', {'smiles': 'CCO', 'include_svg': True})
        assert '<svg' in with_svg.structuredContent['svg']
        bad = await client.call_tool('draw_molecule', {'smiles': 'C1CC'})
        assert bad.isError and 'No structure was drawn' in bad.content[0].text
    stdio_session(test)


def test_stereoisomer_tool():
    async def test(client):
        result = await client.call_tool('enumerate_stereoisomers', {'smiles': 'CC(O)C(=O)O'})
        data = result.structuredContent
        assert data['kind'] == 'stereoisomers' and len(data['isomers']) == 2
        assert data['input_stereo_summary']['status'] == 'unspecified'
        assert [i['enantiomer_index'] for i in data['isomers']] == [1, 0]
        assert len(result.meta['svgs']) == 2
    stdio_session(test)


@pytest.fixture
def http_app():
    mcp._session_manager = None  # a session manager can only run once; build a fresh one per test
    return mcp.streamable_http_app()


def test_http_protocol(http_app):
    app = http_app
    with TestClient(app, base_url='http://localhost:8000') as client:
        response = client.post('/mcp', headers=HEADERS, json=INIT)
        assert response.status_code == 200
        assert response.json()['result']['serverInfo']['name'] == 'MolTalk'
        response = client.post('/mcp', headers=HEADERS, json={'jsonrpc':'2.0', 'id':2, 'method':'tools/call', 'params':{'name':'draw_molecule','arguments':{'smiles':'C/C=C/C'}}})
        assert response.status_code == 200
        assert '<svg' in response.json()['result']['_meta']['svg']


def test_http_rejects_foreign_host_origin_and_large_bodies(http_app):
    app = http_app
    with TestClient(app, base_url='http://localhost:8000') as client:
        assert client.post('/mcp', headers={**HEADERS, 'Host': 'evil.example'}, json=INIT).status_code == 421
        assert client.post('/mcp', headers={**HEADERS, 'Origin': 'https://evil.example'}, json=INIT).status_code == 403
        big = {**INIT, 'params': {**INIT['params'], 'pad': 'x' * 70000}}
        assert client.post('/mcp', headers=HEADERS, json=big).status_code == 413


def test_worker_timeout(monkeypatch):
    monkeypatch.setattr(limits, 'TIMEOUT_S', 0.5)
    start = time.monotonic()
    with pytest.raises(ValueError, match='exceeded'):
        asyncio.run(limits.runner.run(time.sleep, 30))
    assert time.monotonic() - start < 15


def test_name_resolution():
    result = asyncio.run(resolve_name('aspirin'))
    assert result['analysis']['formula'] == 'C9H8O4'
    assert result['stereo_summary']['status'] == 'no stereo elements'
    with pytest.raises(ValueError): asyncio.run(resolve_name('an unknown name'))

def test_pubchem_lookup_mock(monkeypatch):
    import httpx
    import moltalk.server as server
    real_client = httpx.AsyncClient
    def handler(request):
        assert request.url.host == 'pubchem.ncbi.nlm.nih.gov'
        return httpx.Response(200, json={'PropertyTable': {'Properties': [
            {'CID': 702, 'SMILES': 'CCO', 'IUPACName': 'ethanol'}]}})
    monkeypatch.setattr(server.httpx, 'AsyncClient', lambda **kwargs: real_client(transport=httpx.MockTransport(handler), **kwargs))
    result = asyncio.run(resolve_name('zorbatrol', allow_network=True))  # not in the library, not systematic
    assert result['cid'] == 702
    assert result['analysis']['formula'] == 'C2H6O'


network = pytest.mark.skipif(os.getenv('MOLTALK_NETWORK_TESTS') != '1', reason='set MOLTALK_NETWORK_TESTS=1')

@pytest.mark.network
@network
def test_pubchem_live():
    racemic = asyncio.run(resolve_name('lactic acid', allow_network=True))
    assert racemic['cid'] == 612 and racemic['stereo_summary']['status'] == 'unspecified'
    assert 'unspecified' in racemic['warning']
    l_form = asyncio.run(resolve_name('L-lactic acid', allow_network=True))
    assert l_form['source'].startswith('OPSIN') and l_form['analysis']['stereocenters'][0]['cip'] == 'S'  # read offline
    glucose = asyncio.run(resolve_name('D-glucose', allow_network=True))
    assert glucose['stereo_summary']['status'] == 'partially specified'
    with pytest.raises(ValueError, match='not found in MolTalk'):
        asyncio.run(resolve_name('not_a_real_chemical_xyz', allow_network=True))
