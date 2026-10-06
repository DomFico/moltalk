import pytest
from rdkit import Chem
from moltalk.chemistry import draw
from moltalk.depiction import layout

# PubChem CID 123591 (buckminsterfullerene) and PCBM, copied so the test runs offline.
C60 = "c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41"
PCBM = "COC(=O)CCCC1(C23C14C5=C6C7=C8C5=C9C1=C5C%10=C%11C%12=C%13C%10=C%10C1=C8C1=C%10C8=C%10C%14=C%15C%16=C%17C(=C%12C%12=C%17C%17=C%18C%16=C%16C%15=C%15C%10=C1C7=C%15C1=C%16C(=C%18C7=C2C2=C%10C(=C5C9=C42)C%11=C%12C%10=C%177)C3=C16)C%14=C%138)C1=CC=CC=C1"
DODECAHEDRANE = "C12C3C4C5C1C6C7C2C8C3C9C4C%10C5C6C%11C7C8C9C%10%11"


CUBANE, ADAMANTANE = "C12C3C4C1C5C2C3C45", "C1C2CC3CC1CC(C2)C3"


@pytest.mark.parametrize("smiles", [C60, PCBM, DODECAHEDRANE, CUBANE, ADAMANTANE])
def test_cages_are_drawn_as_a_view_of_their_3d_shape(smiles):
    plain = Chem.MolFromSmiles(smiles)
    from rdkit.Chem import rdDepictor
    from moltalk.depiction import quality
    rdDepictor.Compute2DCoords(plain)
    assert quality(plain)["bond_crossings"] or quality(plain)["stretched_bonds"]  # RDKit's flat layout is broken
    for hydrogens in (False, True):
        depiction = draw(smiles, 640, 420, True, hydrogens)["depiction"]
        assert depiction["method"] == "projection" and depiction["max_bond_length_ratio"] <= 2.0  # projection foreshortens bonds
        assert "3D shape" in depiction["note"] and "warning" not in depiction


def test_cage_hydrogens_point_outward():
    # In the old Schlegel diagram the inner atoms' hydrogens piled up in the middle; projected from 3D they point out.
    import numpy as np
    from moltalk.chemistry import depiction_mol
    mol, _ = depiction_mol(CUBANE, True)
    layout(mol)
    pos = mol.GetConformer().GetPositions()[:, :2]
    centre = pos[:8].mean(axis=0)
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1:
            carbon = atom.GetNeighbors()[0].GetIdx()
            outward = pos[carbon] - centre
            assert np.dot(pos[atom.GetIdx()] - pos[carbon], outward) >= -1e-6


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


def test_layout_does_not_depend_on_how_the_smiles_is_written():
    # The model draws resolve_name's canonical SMILES; with that atom order the porphyrin pin was rejected (a methyl
    # ester on the fused ring collided with a side chain) and chlorophyll came out round with Mg–N stretched 2.6x.
    from moltalk import library
    from moltalk.chemistry import analyze, draw
    for name in ("chlorophyll a", "cofactor F430", "heme"):
        smiles = library.find_name(name)[1]["smiles"]
        for written in (smiles, analyze(smiles)["canonical_smiles"]):
            depiction = draw(written, 640, 420, False)["depiction"]
            assert depiction["bond_crossings"] == 0 and depiction["max_bond_length_ratio"] <= 1.6, (name, written)


def test_hydrogens_never_trigger_the_cage_fallback():
    # FAD with explicit hydrogens had a few overlaps, and the Schlegel fallback (meant for C60-like cages) took over.
    from moltalk import library
    from moltalk.chemistry import analyze, draw
    smiles = analyze(library.find_name("FAD")[1]["smiles"])["canonical_smiles"]
    depiction = draw(smiles, 640, 420, False, True)["depiction"]
    assert depiction["method"] == "rdkit" and depiction["max_bond_length_ratio"] <= 1.2


# Porphyrin-type macrocycles, metal and free base, with carbon or nitrogen bridges (PubChem SMILES, offline copies).
MACROCYCLES = {
    "nickel phthalocyanine": "C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=NC2=N3)[N-]8.[Ni+2]",
    "phthalocyanine": "C1=CC=C2C(=C1)C3=NC4=C5C=CC=CC5=C(N4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=NC2=N3)N8",
    "tetraphenylporphyrin": "C1=CC=C(C=C1)C2=C3C=CC(=C(C4=NC(=C(C5=CC=C(N5)C(=C6C=CC2=N6)C7=CC=CC=C7)C8=CC=CC=C8)C=C4)C9=CC=CC=C9)N3",
    "octaethylporphyrin": "CCC1=C(C2=CC3=C(C(=C(N3)C=C4C(=C(C(=N4)C=C5C(=C(C(=N5)C=C1N2)CC)CC)CC)CC)CC)CC)CC",
}


@pytest.mark.parametrize("name", list(MACROCYCLES))
def test_porphyrin_type_macrocycles_get_the_textbook_square(name):
    # Phthalocyanines (nitrogen bridges) and free-base porphyrins were drawn as stretched or round macrocycles; the
    # nickel complex had Ni-N bonds 2.6x normal.
    import numpy as np
    from moltalk.chemistry import depiction_mol
    from moltalk.depiction import _porphyrinoid_core, _shape_rmsd
    smiles = MACROCYCLES[name]
    depiction = draw(smiles, 640, 420, False)["depiction"]
    assert depiction["bond_crossings"] == 0 and depiction["max_bond_length_ratio"] <= 1.6, depiction
    mol, _ = depiction_mol(smiles)
    layout(mol)
    query, template = _porphyrinoid_core()
    match = mol.GetSubstructMatch(query)
    pos = mol.GetConformer().GetPositions()[list(match), :2]
    assert match and _shape_rmsd(pos, template) < 0.1

MACROCYCLES.update({
    "tin phthalocyanine dichloride": "C1=CC=C2C(=C1)C3=NC4=NC(=NC5=C6C=CC=CC6=C(N5)N=C7C8=CC=CC=C8C(=N7)N=C2N3)C9=CC=CC=C94.Cl[Sn]Cl",
    "naphthalocyanine": "c1ccc2cc3c(cc2c1)-c1nc2nc(nc4[nH]c([nH]c5nc(nc-3n1)-c1cc3ccccc3cc1-5)c1cc3ccccc3cc41)-c1cc3ccccc3cc1-2",
})


@pytest.mark.parametrize("name", list(MACROCYCLES))
def test_macrocycles_keep_their_shape_with_hydrogens(name):
    import numpy as np
    from moltalk.chemistry import depiction_mol
    from moltalk.depiction import _porphyrinoid_core, _shape_rmsd
    mol, _ = depiction_mol(MACROCYCLES[name], True)
    layout(mol)
    query, template = _porphyrinoid_core()
    match = mol.GetSubstructMatch(query)
    assert match and _shape_rmsd(mol.GetConformer().GetPositions()[list(match), :2], template) < 0.1


def test_metal_with_its_own_ligands_sits_in_the_ring():
    # PubChem stores tin phthalocyanine dichloride as the free-base ring plus a separate Cl-Sn-Cl: tin floated beside
    # the ring. Now Sn bonds to the four ring nitrogens (replacing the two N-H) and keeps its chlorides.
    from moltalk.chemistry import depiction_mol
    from moltalk.depiction import _chelated_metal
    mol, note = depiction_mol(MACROCYCLES["tin phthalocyanine dichloride"])
    metal, donors = _chelated_metal(mol)
    assert mol.GetAtomWithIdx(metal).GetSymbol() == "Sn" and mol.GetAtomWithIdx(metal).GetDegree() == 6
    assert "replacing 2 N–H" in note


def test_separate_cyanide_is_bonded_to_cobalt_in_b12():
    from moltalk import library
    from moltalk.chemistry import depiction_mol
    from moltalk.depiction import _chelated_metal
    mol, note = depiction_mol(library.find_name("cyanocobalamin")[1]["smiles"])
    metal, donors = _chelated_metal(mol)
    assert {n.GetSymbol() for n in mol.GetAtomWithIdx(metal).GetNeighbors()} == {"N", "C"} and "ligand is bonded" in note
