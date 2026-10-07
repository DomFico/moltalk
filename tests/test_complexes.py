"""Lone pairs and radicals against the real 3D geometry, VSEPR shapes the force field gets wrong, and metal complexes
stored as loose pieces (assembled, counted, drawn and built in their ideal geometry)."""
import itertools

import numpy as np
import pytest
from rdkit import Chem

from moltalk import chemistry, complexes

WILKINSON = "[Cl-].[Rh].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
PDPPH34 = ("[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1."
           "c1ccc(P(c2ccccc2)c2ccccc2)cc1")
FERROCENE = "[CH]1C=CC=C1.[CH]1C=CC=C1.[Fe]"
ZIRCONOCENE = "[CH]1C=CC=C1.[CH]1C=CC=C1.[Cl][Zr][Cl]"
ZEISE = "C=C.O.[Cl-].[Cl-].[Cl-].[K+].[Pt+2]"
RUBPY = "[Cl-].[Cl-].[Ru+2].c1ccc(-c2ccccn2)nc1.c1ccc(-c2ccccn2)nc1.c1ccc(-c2ccccn2)nc1"
JACOBSEN = ("CC(C)(C)c1cc(C=N[C@@H]2CCCC[C@H]2N=Cc2cc(C(C)(C)C)cc(C(C)(C)C)c2[O-])c([O-])c(C(C)(C)C)c1.[Cl-].[Mn+3]")
FECO5 = "[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[Fe]"
MOCO6 = "[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[C-]#[O+].[Mo]"
CISPLATIN = "N.N.[Cl][Pt][Cl]"
OXALIPLATIN = "O=C([O-])C(=O)[O-].[NH-][C@@H]1CCCC[C@H]1[NH-].[Pt+4]"
GRUBBS = "Cl[Ru](Cl)(=Cc1ccccc1)(P(C1CCCCC1)(C1CCCCC1)C1CCCCC1)P(C1CCCCC1)(C1CCCCC1)C1CCCCC1"


def _model(smiles, hydrogens=True):
    model = chemistry.conformer(smiles, hydrogens)
    xyz = {a["index"]: np.array(a["xyz"]) for a in model["atoms"]}
    neighbours = {}
    for b in model["bonds"]:
        neighbours.setdefault(b["a"], []).append(b["b"])
        neighbours.setdefault(b["b"], []).append(b["a"])
    return model, xyz, neighbours


def _unit(v):
    return v / np.linalg.norm(v)


# ---- lone pairs and radicals --------------------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", [
    "[CH2-][P+](c1ccccc1)(c1ccccc1)c1ccccc1",  # Wittig ylide: the CH2- pair pointed 14 degrees from a C-H bond
    "[CH2-][S+](C)C", "[CH2-][S+](C)(C)=O", "[CH2-][N+](C)(C)C", "[CH3-]", "[CH2-]C(C)=O", "C=C(C)[O-]",
    "[cH-]1cccc1", "C=[N+]=[N-]", "CN1C=CN(C)[C]1", "Cl[C]Cl", "C=[N+](C)[O-]", "C[N+](C)(C)[O-]", "Nc1ccccc1",
    "CC(=O)N(C)C", "CS(C)=O", "CP(C)C", "O", "N", "CN=[N+]=[N-]", "CC=NO",
    "[CH3]", "C[C](C)C", "[CH2]c1ccccc1", "C=[CH]", "[CH2]", "CC1(C)CCCC(C)(C)N1[O]", "[O]c1ccccc1", "O=[N][O]",
])
def test_lone_pairs_and_radicals_never_sit_on_a_bond(smiles):
    model, xyz, neighbours = _model(smiles)
    for key, directions in model["lone_pair_dirs"].items():
        i = int(key)
        bonds = [_unit(xyz[j] - xyz[i]) for j in neighbours.get(i, [])]
        for d in directions:
            for b in bonds:
                angle = np.degrees(np.arccos(np.clip(np.dot(_unit(np.array(d)), b), -1, 1)))
                assert angle >= 80, (smiles, i, angle)


@pytest.mark.parametrize("smiles, centre, expected", [
    ("[CH3-]", 0, "pyramidal"), ("[CH2-][N+](C)(C)C", 0, "pyramidal"), ("[CH2-][S+](C)C", 0, "pyramidal"),
    ("[CH2-][P+](c1ccccc1)(c1ccccc1)c1ccccc1", 0, "planar"),  # phosphonium ylide carbon: near planar
    ("[CH2-]c1ccccc1", 0, "planar"), ("[CH3]", 0, "planar"),
])
def test_three_bond_centres_have_the_vsepr_shape(smiles, centre, expected):
    _, xyz, neighbours = _model(smiles)
    total = np.linalg.norm(sum(_unit(xyz[j] - xyz[centre]) for j in neighbours[centre]))
    assert (total > 0.3) == (expected == "pyramidal"), total


@pytest.mark.parametrize("smiles, centre, low, high", [
    ("[CH2]", 0, 125, 145),  # triplet methylene (no heteroatom neighbour), not linear
    ("Cl[C]Cl", 1, 98, 115),  # singlet dichlorocarbene
    ("C=[CH]", 1, 125, 145),  # vinyl radical: a bent sigma radical
])
def test_carbenes_and_sigma_radicals_are_bent(smiles, centre, low, high):
    _, xyz, neighbours = _model(smiles)
    a, b = (_unit(xyz[j] - xyz[centre]) for j in neighbours[centre])
    assert low <= np.degrees(np.arccos(np.dot(a, b))) <= high


def test_single_heavy_atom_models_keep_their_hydrogens_apart():
    for smiles in ("O", "N", "[CH3]", "C"):
        _, xyz, _ = _model(smiles)
        points = list(xyz.values())
        assert min(np.linalg.norm(p - q) for p, q in itertools.combinations(points, 2)) > 0.5


def test_charged_carbon_keeps_its_hydrogens_in_the_label():
    import re
    svg = chemistry.draw("[CH2-][P+](C)(C)C")["svg"]
    assert len(re.findall(r"class='atom-0'", svg)) >= 4  # C, H, 2, charge (RDKit alone wrote just 'C' and '-')


# ---- metal complexes ----------------------------------------------------------------------------------------------

@pytest.mark.parametrize("smiles, geometry, oxidation, electrons", [
    (WILKINSON, "square_planar", 1, 16), (PDPPH34, "tetrahedral", 0, 18), (FERROCENE, "linear", 2, 18),
    (ZIRCONOCENE, "tetrahedral", 4, 16), (ZEISE, "square_planar", 2, 16), (RUBPY, "octahedral", 2, 18),
    (JACOBSEN, "square_pyramidal", 3, 14), (FECO5, "trigonal_bipyramidal", 0, 18), (MOCO6, "octahedral", 0, 18),
    (CISPLATIN, "square_planar", 2, 16), (OXALIPLATIN, "square_planar", 2, 16),
])
def test_complexes_are_assembled_by_electron_counting(smiles, geometry, oxidation, electrons):
    drawn, note = chemistry.depiction_mol(smiles)
    info = complexes.info(drawn)
    assert info["geometry"] == geometry and info["oxidation_state"] == oxidation
    assert info["electron_count"] == electrons
    assert "separate pieces" in note and "source unchanged" in note
    analysis = chemistry.analyze(smiles)
    assert analysis["canonical_smiles"] == Chem.MolToSmiles(Chem.MolFromSmiles(smiles))  # the analysis is untouched


def test_counter_ions_and_solvent_stay_separate():
    drawn, _ = chemistry.depiction_mol(RUBPY)  # three bipyridines make Ru(II) 18-electron: the chlorides stay ions
    ru = next(a for a in drawn.GetAtoms() if a.GetSymbol() == "Ru")
    assert sorted(n.GetSymbol() for n in ru.GetNeighbors()) == ["N"] * 6
    drawn, _ = chemistry.depiction_mol(ZEISE)  # K+ and the water are not ligands
    pt = next(a for a in drawn.GetAtoms() if a.GetSymbol() == "Pt")
    assert all(n.GetSymbol() in ("Cl", "*") for n in pt.GetNeighbors())


def test_multi_metal_records_are_left_as_stored():
    tebbe = "[CH3-].[CH3][Al+][CH3].[Cl-].[Ti+3].c1cc[cH-]c1.c1cc[cH-]c1"
    drawn, note = chemistry.depiction_mol(tebbe)
    assert complexes.info(drawn) is None and "several metal atoms" in note


def test_wilkinson_is_drawn_joined_and_clean():
    drawn = chemistry.draw(WILKINSON)
    assert drawn["depiction"]["bond_crossings"] == 0 and drawn["depiction"]["overlapping_atoms"] == 0
    mol, _ = chemistry.depiction_mol(WILKINSON)
    assert len(Chem.GetMolFrags(mol)) == 1


@pytest.mark.parametrize("smiles", [FERROCENE, ZIRCONOCENE])
def test_metallocenes_are_drawn_with_ring_centre_bonds(smiles):
    depiction = chemistry.draw(smiles)["depiction"]
    assert depiction["bond_crossings"] == 0 and depiction["overlapping_atoms"] == 0
    mol, _ = chemistry.depiction_mol(smiles)
    assert sum(a.HasProp("_eta_atoms") for a in mol.GetAtoms()) == 2
    assert not any(a.GetNumRadicalElectrons() for a in mol.GetAtoms())  # no stray radical dots on the Cp rings


def test_bound_carbonyls_are_written_without_charges():
    mol, _ = chemistry.depiction_mol(MOCO6)
    assert not any(a.GetFormalCharge() or a.GetNumRadicalElectrons() for a in mol.GetAtoms())


_IDEAL = {"square_planar": [90, 90, 90, 90, 180, 180], "tetrahedral": [109.5] * 6, "linear": [180],
          "octahedral": [90] * 12 + [180] * 3, "trigonal_bipyramidal": [90] * 6 + [120] * 3 + [180]}


@pytest.mark.parametrize("smiles, geometry", [
    (WILKINSON, "square_planar"), (PDPPH34, "tetrahedral"), (FERROCENE, "linear"), (ZEISE, "square_planar"),
    (FECO5, "trigonal_bipyramidal"), (MOCO6, "octahedral"), (RUBPY, "octahedral"), (CISPLATIN, "square_planar"),
])
def test_3d_models_have_the_ideal_metal_geometry(smiles, geometry):
    model, xyz, neighbours = _model(smiles, hydrogens=False)
    metal = next(a["index"] for a in model["atoms"] if a["element"] in ("Rh", "Pd", "Fe", "Pt", "Mo", "Ru"))
    vectors = [_unit(xyz[j] - xyz[metal]) for j in neighbours[metal]]
    angles = sorted(np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))) for a, b in itertools.combinations(vectors, 2))
    tolerance = 12 if smiles == RUBPY else 2  # a bipyridine's bite angle is smaller than 90
    assert np.allclose(angles, _IDEAL[geometry], atol=tolerance), angles
    assert "around the metal" in model["method"]


def test_carbonyls_point_outward_in_3d():
    model, xyz, neighbours = _model(FECO5, hydrogens=False)
    fe = next(a["index"] for a in model["atoms"] if a["element"] == "Fe")
    for c in neighbours[fe]:
        o = next(j for j in neighbours[c] if j != fe)
        assert np.linalg.norm(xyz[o] - xyz[fe]) > np.linalg.norm(xyz[c] - xyz[fe])  # Fe–C≡O, not Fe–O


def test_already_bonded_complex_gets_its_geometry_too():
    model, xyz, neighbours = _model(GRUBBS, hydrogens=False)
    ru = next(a["index"] for a in model["atoms"] if a["element"] == "Ru")
    vectors = [_unit(xyz[j] - xyz[ru]) for j in neighbours[ru]]
    angles = sorted(np.degrees(np.arccos(np.clip(np.dot(a, b), -1, 1))) for a, b in itertools.combinations(vectors, 2))
    assert np.allclose(angles, [90] * 8 + [180] * 2, atol=2)  # square pyramid (d6 Ru), was 19 to 119 degrees


def test_jacobsen_keeps_its_stereocentres_in_3d():
    model = chemistry.conformer(JACOBSEN)
    mol = chemistry.parse(JACOBSEN)
    wanted = {a.GetIdx(): a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")}
    got = {a["index"]: a["cip"] for a in model["atoms"] if a.get("cip")}
    assert all(got.get(i) == cip for i, cip in wanted.items())
