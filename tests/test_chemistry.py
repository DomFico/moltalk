import pytest
from moltalk.chemistry import analyze, draw, enumerate_stereo, substructure

def test_properties_and_groups():
    result = analyze('CC(=O)Oc1ccccc1C(=O)O')
    assert result['formula'] == 'C9H8O4'
    assert set(result['functional_groups']) == {'ester', 'carboxylic acid', 'aromatic ring'}  # ester O is not an ether
    assert analyze('COc1ccccc1')['functional_groups']['ether']  # anisole still is
    assert result['properties']['molecular_weight'] == pytest.approx(180.159, abs=.01)

@pytest.mark.parametrize('smiles', ['C(C)(C)(C)(C)C', 'nonsense', ''])
def test_invalid(smiles):
    with pytest.raises(ValueError): analyze(smiles)

def test_enantiomers():
    a = analyze('N[C@@H](C)C(=O)O')
    b = analyze('N[C@H](C)C(=O)O')
    assert a['stereocenters'][0]['cip'] == 'S'
    assert b['stereocenters'][0]['cip'] == 'R'
    assert a['inchikey'] != b['inchikey']

def test_unspecified():
    assert analyze('CC(O)C(=O)O')['warnings']
    assert len(enumerate_stereo('CC(O)C(=O)O')['isomers']) == 2
    assert enumerate_stereo('CC(O)C(=O)O', 1)['truncated']

def test_double_stereo():
    assert analyze('C/C=C/C')['double_bond_stereo'][0]['configuration'] == 'E'
    assert analyze('C/C=C\\C')['double_bond_stereo'][0]['configuration'] == 'Z'

def test_drawing_and_query():
    result = draw('N[C@@H](C)C(=O)O')
    assert '<svg' in result['svg']
    assert result['depicted_stereo_bonds']
    assert substructure('CCO', '[OX2H]')['matches'] == [[2]]
    with pytest.raises(ValueError): substructure('CCO', '[')
    with pytest.raises(ValueError): draw('CCO', 1, 420)

def test_requested_cyclohexane_example():
    result = analyze('F[C@H]1C[C@@H](C)CCC1')
    assert {c['atom_index']: c['cip'] for c in result['stereocenters']} == {1: 'R', 3: 'S'}
    assert result['formula'] == 'C7H13F'

def test_invalid_errors_are_explanatory():
    with pytest.raises(ValueError, match='syntax'): analyze('C1CC(')
    with pytest.raises(ValueError, match='valence'): analyze('C(C)(C)(C)(C)C')

def test_wedges_explained_and_distinct_from_cip():
    bonds = draw('F[C@H]1C[C@@H](C)CCC1')['depicted_stereo_bonds']
    # Both bonds are wedges in this depiction although the centers are R and S.
    assert sorted((b['style'], b['stereocenter_cip']) for b in bonds) == [('wedge', 'R'), ('wedge', 'S')]
    assert all('toward the viewer' in b['explanation'] for b in bonds)

def test_meso_and_enantiomers():
    result = enumerate_stereo('OC(=O)C(O)C(O)C(=O)O')  # tartaric acid
    assert len(result['isomers']) == 3
    assert sum(i['achiral'] for i in result['isomers']) == 1
    chiral = [i for i in result['isomers'] if not i['achiral']]
    assert {chiral[0]['enantiomer_index'], chiral[1]['enantiomer_index']} == {chiral[0]['index'], chiral[1]['index']}

def test_stereo_summary():
    assert analyze('F[C@H]1C[C@@H](C)CCC1')['stereo_summary']['status'] == 'fully specified'
    assert analyze('F[C@H]1CC(C)CCC1')['stereo_summary']['status'] == 'partially specified'
