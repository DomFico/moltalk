import pytest
from rdkit import Chem
from moltalk.chemistry import draw
from moltalk.depiction import layout

# PubChem CID 123591 (buckminsterfullerene) and PCBM, copied so the test runs offline.
C60 = "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41"
PCBM = "COC(=O)CCCC1(C23C14C5=C6C7=C8C5=C9C1=C5C%10=C%11C%12=C%13C%10=C%10C1=C8C1=C%10C8=C%10C%14=C%15C%16=C%17C(=C%12C%12=C%17C%17=C%18C%16=C%16C%15=C%15C%10=C1C7=C%15C1=C%16C(=C%18C7=C2C2=C%10C(=C5C9=C42)C%11=C%12C%10=C%177)C3=C16)C%14=C%138)C1=CC=CC=C1"
DODECAHEDRANE = "C12C3C4C5C1C6C7C2C8C3C9C4C%10C5C6C%11C7C8C9C%10%11"


@pytest.mark.parametrize("smiles", [C60, PCBM, DODECAHEDRANE])
def test_cages_get_crossing_free_schlegel_layout(smiles):
    plain = Chem.MolFromSmiles(smiles)
    from rdkit.Chem import rdDepictor
    from moltalk.depiction import quality
    rdDepictor.Compute2DCoords(plain)
    assert quality(plain)["bond_crossings"] > 0  # RDKit's default layout of these cages is broken
    depiction = draw(smiles)["depiction"]
    assert depiction["method"] == "schlegel"
    assert depiction["bond_crossings"] == 0 and depiction["overlapping_atoms"] == 0
    assert "viewed through one ring" in depiction["note"] and "warning" not in depiction


def test_c60_graph_is_intact():
    mol = Chem.MolFromSmiles(C60)
    assert mol.GetNumAtoms() == 60 and mol.GetNumBonds() == 90
    assert all(a.GetDegree() == 3 for a in mol.GetAtoms())


@pytest.mark.parametrize("smiles", ["CC(=O)Oc1ccccc1C(=O)O", "F[C@H]1C[C@@H](C)CCC1", "C1COCCOCCOCCOCCOCCO1",
                                    "CC(C)CCCC(C)C1CCC2C1(CCC3C2CC=C4C3(CCC(C4)O)C)C"])
def test_ordinary_molecules_keep_rdkit_layout(smiles):
    mol = Chem.MolFromSmiles(smiles)
    result = layout(mol)
    assert result["method"] == "rdkit" and result["bond_crossings"] == 0


def test_hydroporphyrins_get_the_textbook_square():
    # F430 and chlorophyll have rings fused to the macrocycle that defeat RDKit's porphyrin template (the ring folded
    # inward and the metal bonds stretched to 2.6x); the core is now pinned to the porphine square.
    from moltalk import library
    from moltalk.chemistry import draw
    for name in ("cofactor F430", "chlorophyll a"):
        depiction = draw(library.find_name(name)[1]["smiles"], 640, 420, False)["depiction"]
        assert depiction["bond_crossings"] == 0 and depiction["overlapping_atoms"] == 0, name
        assert depiction["max_bond_length_ratio"] <= 1.6, (name, depiction)
