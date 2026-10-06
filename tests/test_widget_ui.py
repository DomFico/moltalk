"""Render the MCP Apps widget in a real browser inside a sandboxed iframe, driven by the real server.

The host page imitates what an MCP Apps host (such as ChatGPT) does: it answers ui/initialize,
delivers tool-input/tool-result notifications and proxies the widget's tools/call requests.
Requires the ui-test extra and Google Chrome (or set MOLTALK_UI_BROWSER_CHANNEL=chromium).
"""
import asyncio
import json
import os
import sys
import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from moltalk.server import WIDGET_URI

playwright_api = pytest.importorskip("playwright.async_api")
SCREENSHOT_DIR = os.getenv("MOLTALK_UI_SCREENSHOTS")

HOST_PAGE = """<!doctype html><meta charset="utf-8"><body style="margin:0;background:#fff">
<iframe id="w" sandbox="allow-scripts" style="width:720px;height:900px;border:0"></iframe>
<script>
window.hostLog = [];
const frame = document.getElementById('w');
let toolCall = null;
window.startWidget = (html, args, result, theme) => {
  toolCall = { args, result, theme };
  frame.srcdoc = html;
};
window.addEventListener('message', async (event) => {
  if (event.source !== frame.contentWindow) return;
  const m = event.data; window.hostLog.push(m.method || 'response');
  const reply = (body) => frame.contentWindow.postMessage({ jsonrpc: '2.0', id: m.id, ...body }, '*');
  if (m.method === 'ui/initialize') {
    reply({ result: { protocolVersion: m.params.protocolVersion, hostInfo: { name: 'test-host', version: '1' },
                      hostCapabilities: {}, hostContext: { theme: toolCall.theme } } });
  } else if (m.method === 'ui/notifications/initialized') {
    const send = (method, params) => frame.contentWindow.postMessage({ jsonrpc: '2.0', method, params }, '*');
    send('ui/notifications/tool-input', { arguments: toolCall.args });
    send('ui/notifications/tool-result', toolCall.result);
  } else if (m.method === 'tools/call') {
    try { reply({ result: JSON.parse(await window.pyCallTool(JSON.stringify(m.params))) }); }
    catch (e) { reply({ error: { code: -32000, message: String(e) } }); }
  } else if (m.method === 'ui/notifications/size-changed') {
    window.lastHeight = m.params.height;
  }
});
</script>"""


def dump(result):
    return result.model_dump(mode="json", by_alias=True, exclude_none=True)


def test_widget_renders_and_calls_tools():
    async def run():
        params = StdioServerParameters(command=sys.executable, args=["-m", "moltalk.server"])
        async with stdio_client(params) as (read, write), ClientSession(read, write) as client:
            await client.initialize()
            html = (await client.read_resource(WIDGET_URI)).contents[0].text
            tool_calls = []

            async def py_call_tool(raw):
                request = json.loads(raw)
                tool_calls.append(request)
                return json.dumps(dump(await client.call_tool(request["name"], request["arguments"])))

            async with playwright_api.async_playwright() as pw:
                channel = os.getenv("MOLTALK_UI_BROWSER_CHANNEL", "chrome")
                browser = await pw.chromium.launch(channel=None if channel == "chromium" else channel)
                page = await browser.new_page(viewport={"width": 760, "height": 940})
                await page.expose_function("pyCallTool", py_call_tool)
                await page.set_content(HOST_PAGE)
                frame = page.frame_locator("#w")

                async def show(name, args, theme="light"):
                    result = dump(await client.call_tool(name, args))
                    await page.evaluate("([h, a, r, t]) => startWidget(h, a, r, t)", [html, args, result, theme])

                # 1. Stereo drawing: SVG, CIP chips and wedge list.
                args = {"smiles": "F[C@H]1C[C@@H](C)CCC1", "label": "(1R,3S)-1-fluoro-3-methylcyclohexane"}
                await show("draw_molecule", args)
                await frame.locator(".canvas svg path").first.wait_for()
                body = await frame.locator("main").inner_text()
                assert "(1R,3S)-1-fluoro-3-methylcyclohexane" in body
                assert ("atom 1 (C): R" in body and "atom 3 (C): S" in body) or ("C1: R" in body and "C3: S" in body)
                assert await frame.locator("ul.bonds li .wedge").count() == 2
                assert await page.evaluate("window.lastHeight") > 300
                if SCREENSHOT_DIR:
                    await page.screenshot(path=os.path.join(SCREENSHOT_DIR, "draw.png"), full_page=True)

                # 2. Follow-up tool call from inside the widget (toggle atom indices).
                svg_before = await frame.locator(".canvas").inner_html()
                await frame.get_by_role("group", name="Atom numbers").get_by_role("button", name="Off", exact=True).click()
                await frame.get_by_role("group", name="Atom numbers").get_by_role("button", name="Off", exact=True).and_(frame.locator("[aria-pressed=true]")).wait_for()
                assert tool_calls[-1] == {"name": "draw_molecule", "arguments": {**args, "atom_indices": True, "numbering": "iupac",
                                                                                     "hydrogens": False, "numbering": "none"}}
                assert await frame.locator(".canvas").inner_html() != svg_before

                # 2b. Cage: Schlegel layout with its explanatory note.
                from test_depiction import C60
                await show("draw_molecule", {"smiles": C60, "label": "Buckminsterfullerene"})
                await frame.get_by_text("Schlegel diagram.").wait_for()
                if SCREENSHOT_DIR:
                    await page.screenshot(path=os.path.join(SCREENSHOT_DIR, "c60.png"), full_page=True)

                # 3. Stereoisomer grid (dark theme).
                await show("enumerate_stereoisomers", {"smiles": "CC(O)C(=O)O"}, theme="dark")
                await frame.locator(".card svg").nth(1).wait_for()
                assert await frame.locator(".card").count() == 2
                grid_text = await frame.locator("main").inner_text()
                assert "enantiomer of #1" in grid_text and "enantiomer of #0" in grid_text
                assert await frame.locator("html").get_attribute("data-theme") == "dark"
                if SCREENSHOT_DIR:
                    await page.screenshot(path=os.path.join(SCREENSHOT_DIR, "isomers.png"), full_page=True)

                # 4. Invalid structure: error card, nothing drawn.
                await show("draw_molecule", {"smiles": "C1CC"})
                await frame.locator(".error").wait_for()
                assert "No structure drawn" in await frame.locator(".error").inner_text()
                assert await frame.locator("svg").count() == 0
                if SCREENSHOT_DIR:
                    await page.screenshot(path=os.path.join(SCREENSHOT_DIR, "error.png"), full_page=True)

                assert "ui/initialize" in await page.evaluate("window.hostLog")
                await browser.close()

    asyncio.run(run())


def test_widget_starts_from_window_openai_globals():
    """ChatGPT can inject window.openai (toolOutput/_meta) before the widget script runs."""
    async def run():
        params = StdioServerParameters(command=sys.executable, args=["-m", "moltalk.server"])
        async with stdio_client(params) as (read, write), ClientSession(read, write) as client:
            await client.initialize()
            html = (await client.read_resource(WIDGET_URI)).contents[0].text
            result = dump(await client.call_tool("draw_molecule", {"smiles": "C1C2CC3CC1CC(C2)C3", "label": "adamantane"}))
            globals_ = {"theme": "dark", "toolInput": {"smiles": "C1C2CC3CC1CC(C2)C3"},
                        "toolOutput": result["structuredContent"], "toolResponseMetadata": result["_meta"]}
            page_html = html.replace("<head>", "<head><script>window.openai = " + json.dumps(globals_) + ";</script>", 1)
            async with playwright_api.async_playwright() as pw:
                channel = os.getenv("MOLTALK_UI_BROWSER_CHANNEL", "chrome")
                browser = await pw.chromium.launch(channel=None if channel == "chromium" else channel)
                page = await browser.new_page()
                errors = []
                page.on("pageerror", lambda e: errors.append(str(e)))
                await page.set_content(page_html)
                await page.locator(".liftable svg path").first.wait_for()
                assert errors == []
                assert await page.locator("html").get_attribute("data-theme") == "dark"
                await browser.close()
    asyncio.run(run())
