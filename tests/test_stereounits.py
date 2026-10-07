"""Stereogenic units beyond R/S and E/Z: allenes, spiro compounds, helicenes and C–N axes. Positive and negative
controls, configurations through parsing, drawing, 3D and export, descriptors against IUPAC reference examples, and
the crowded-layout panel."""
import asyncio

import numpy as np
import pytest
from rdkit import Chem
from rdkit.Chem import AllChem

from moltalk import chemistry, stereounits
from moltalk.depiction import layout
from moltalk.export import export_structure
from moltalk.server import resolve_name

DICHLOROALLENE = "ClC=C=CCl"
SPIRODIONE = "O=C1CCCC12CCCC2=O"
DICHLOROSPIROHEPTANE = "ClC1CC2(C1)CC(Cl)C2"
HELICENE6 = "c1ccc2c(c1)ccc1ccc3ccc4ccc5ccccc5c4c3c12"
HELICENE5 = "c1ccc2c(c1)ccc1ccc3ccc4ccccc4c3c12"
METOLACHLOR = "CCc1cccc(C)c1N(C(C)COC)C(=O)CCl"


def _isomers(smiles):
    return chemistry.enumerate_stereo(smiles)["isomers"]


def _units(smiles, kind):
    return [u for u in chemistry.analyze(smiles)["stereo_units"] if u["type"] == kind]


# ---- sign conventions, checked against independent references -----------------------------------------------------

def test_tetrahedral_rule_matches_rdkit_on_l_alanine():
    mol = Chem.AddHs(Chem.MolFromSmiles("N[C@@H](C)C(=O)O"))  # L-alanine, (S)
    AllChem.EmbedMolecule(mol, randomSeed=3)
    p = mol.GetConformer().GetPositions()
    h = next(n.GetIdx() for n in mol.GetAtomWithIdx(1).GetNeighbors() if n.GetAtomicNum() == 1)
    assert stereounits.handedness(p[0], p[3], p[2], p[h]) == "S"  # N > COOH > CH3 > H


def test_right_handed_helix_is_p():
    t = np.arange(4.0)
    points = [np.array([np.cos(x), np.sin(x), 0.5 * x]) for x in t]  # a right-handed helix (IUPAC: P)
    assert stereounits.torsion(*points) > 0


def test_iupac_spiro_example_2r4s6r_is_2p():
    # IUPAC P-93.5.3.5, example 1: (2R,4S,6R)-2,6-dichlorospiro[3.3]heptane = (2P)-2,6-dichlorospiro[3.3]heptane.
    found = {}
    for iso in _isomers(DICHLOROSPIROHEPTANE):
        unit = _units(iso["smiles"], "spiro")[0]
        found[unit["descriptor"]] = unit["alternative_descriptor"]
    assert found == {"1R,3S,5R": "P", "1S,3R,5S": "M"}  # input indices 1, 3, 5 are IUPAC locants 2, 4, 6


@pytest.mark.parametrize("smiles", [DICHLOROALLENE, "CC=C=CC", "OC(=O)C=C=CC"])
def test_allene_rules_agree_m_is_ra(smiles):
    # IUPAC P-92.1.2.2.3 (for a chirality axis M = Ra) and P-93.4.2.2: (1M)-1,3-dichloropropa-1,2-diene = (1Ra).
    pairs = {(u["descriptor"], u["alternative_descriptor"]) for iso in _isomers(smiles)
             for u in _units(iso["smiles"], "allene")}
    assert pairs == {("M", "Ra"), ("P", "Sa")}


# ---- detection: positive and negative controls --------------------------------------------------------------------

@pytest.mark.parametrize("smiles, kind, expected", [
    (DICHLOROALLENE, "allene", 1), ("CC=C=CC", "allene", 1), ("ClC=C=C=C=CCl", "allene", 1),
    ("C=C=C", "allene", 0), ("CC(C)=C=CCl", "allene", 0), ("ClC=C=C=CCl", "allene", 0),  # odd cumulene: planar, E/Z
    (SPIRODIONE, "spiro", 1), ("C1CC12CC2", "spiro", 0), ("C1=CC12C=C2", "spiro", 0), ("C1CC2(C1)CCC2", "spiro", 0),
    ("ClC1CC2(C1)CCC2", "spiro", 0),
    (HELICENE5, "helix", 1), (HELICENE6, "helix", 1), ("c1ccc2c(c1)ccc1ccc3ccccc3c12", "helix", 0),  # [4]helicene
    ("c1ccc2c(c1)ccc1c2ccc2c3ccccc3ccc21", "helix", 0),  # picene: zigzag, not helical
    ("c1ccc2cc3cc4cc5ccccc5cc4cc3cc2c1", "helix", 0),  # pentacene
    ("c12c3c4c5c1c1c6c7c2c2c8c3c3c9c4c4c%10c5c5c1c1c6c6c%11c7c2c2c7c8c3c3c8c9c4c4c9c%10c5c5c1c1c6c6c%11c2c2c7c3c3c8c4c4c9c5c1c1c6c2c3c41",
     "helix", 0),  # C60
])
def test_units_detected_only_where_stereogenic(smiles, kind, expected):
    assert len([u for u in stereounits.units(Chem.MolFromSmiles(smiles)) if u["type"] == kind]) == expected


@pytest.mark.parametrize("smiles", ["C1CC12CC2", "C1=CC12C=C2", "C1CC2(C1)CCC2", "C=C=C", "CC(=O)Nc1ccccc1",
                                    "CC(=O)N(C)c1ccccc1", "CC(=O)N(C)c1c(C)cccc1C", "CN(C)C=O"])
def test_achiral_controls_have_no_stereo(smiles):
    analysis = chemistry.analyze(smiles)
    assert analysis["stereo_summary"]["status"] == "no stereo elements"
    assert analysis["stereocenters"] == [] and analysis["stereo_units"] == []


def test_cn_axis_of_a_tertiary_anilide():
    units = _units(METOLACHLOR, "axis")
    assert len(units) == 1 and units[0]["subtype"] == "C–N" and not units[0]["specified"]
    assert "rule of thumb" in units[0]["stability"]
    isomers = _isomers(METOLACHLOR)  # one open centre and one open axis: four stereoisomers
    assert len(isomers) == 4
    assert {(iso["stereocenters"][0]["cip"], iso["axial_stereo"][0]["cip"]) for iso in isomers} == \
        {("R", "Ra"), ("R", "Sa"), ("S", "Ra"), ("S", "Sa")}


def test_helicene_stability_is_separate_from_stereogenicity():
    five, six = _units(HELICENE5, "helix")[0], _units(HELICENE6, "helix")[0]
    assert "racemises slowly" in five["stability"] and "stable" in six["stability"]


# ---- configurations: enumeration, round trips, drawing, 3D --------------------------------------------------------

@pytest.mark.parametrize("smiles, kind, configs", [
    (DICHLOROALLENE, "allene", {"M", "P"}), (HELICENE6, "helix", {"M", "P"}),
    (SPIRODIONE, "spiro", {"configuration 1", "configuration 2"}),
])
def test_two_enantiomers_round_trip_through_cxsmiles(smiles, kind, configs):
    isomers = _isomers(smiles)
    assert len(isomers) == 2 and {i["enantiomer_index"] for i in isomers} == {0, 1}
    seen = set()
    for iso in isomers:
        assert "|(" in iso["smiles"]  # CXSMILES with 3D coordinates
        unit = _units(iso["smiles"], kind)[0]
        assert unit["specified"] and chemistry.analyze(iso["smiles"])["canonical_smiles"]
        seen.add(unit["descriptor"] or unit["configuration"])
        again = _units(chemistry.analyze(iso["smiles"])["canonical_smiles"], kind)[0]
        assert (again["descriptor"] or again["configuration"]) == (unit["descriptor"] or unit["configuration"])
    assert seen == configs


def test_spiro_xabab_has_no_guessed_descriptor():
    unit = _units(_isomers(SPIRODIONE)[0]["smiles"], "spiro")[0]
    assert unit["descriptor"] is None and unit["verification"] == "unverified" and "P-93.5.3.2" in unit["note"]


def test_open_units_are_reported_not_assumed():
    for smiles, kind in ((DICHLOROALLENE, "allene"), (HELICENE6, "helix"), (SPIRODIONE, "spiro")):
        analysis = chemistry.analyze(smiles)
        assert analysis["stereo_summary"]["unspecified_units"][0]["type"] == kind
        assert "|" not in analysis["canonical_smiles"]
        model = chemistry.conformer(smiles)
        assert "arbitrary configuration" in model["warning"]


@pytest.mark.parametrize("smiles", [DICHLOROALLENE, SPIRODIONE, HELICENE6, "C[C@H](O)C=C=CCl"])
def test_3d_model_has_the_stated_configuration(smiles):
    for iso in _isomers(smiles):
        mol = chemistry.parse(iso["smiles"])
        model = chemistry.conformer(iso["smiles"])
        conf = Chem.Conformer(mol.GetNumAtoms())
        for atom in model["atoms"][:mol.GetNumAtoms()]:
            conf.SetAtomPosition(atom["index"], atom["xyz"])
        stated = stereounits.stated(mol)
        for unit in stereounits.units(mol):
            assert stereounits.measure(mol, conf, unit)["sense"] == stated[unit["key"]]


def test_allene_stays_twisted_in_3d():
    # MMFF94 and UFF flatten allenes; the model holds the ends perpendicular.
    model = chemistry.conformer("C=C=C")
    xyz = {a["index"]: np.array(a["xyz"]) for a in model["atoms"]}
    mol = Chem.AddHs(Chem.MolFromSmiles("C=C=C"))
    assert "perpendicular" in model["method"]
    model_h = chemistry.conformer("C=C=C", True)
    pos = {a["index"]: np.array(a["xyz"]) for a in model_h["atoms"]}
    h0 = [n.GetIdx() for n in mol.GetAtomWithIdx(0).GetNeighbors() if n.GetAtomicNum() == 1]
    h2 = [n.GetIdx() for n in mol.GetAtomWithIdx(2).GetNeighbors() if n.GetAtomicNum() == 1]
    assert 70 < abs(stereounits.torsion(pos[h0[0]], pos[0], pos[2], pos[h2[0]])) < 110
    assert xyz


@pytest.mark.parametrize("smiles", [DICHLOROALLENE, SPIRODIONE])
def test_drawn_wedges_agree_with_the_configuration(smiles):
    for iso in _isomers(smiles):
        drawn = chemistry.draw(iso["smiles"])
        unit_wedges = [b for b in drawn["depicted_stereo_bonds"] if b.get("unit")]
        assert unit_wedges and "warning" not in drawn["depiction"]
        assert all("Not a stereocentre" in b["explanation"] or "spiro atom" in b["explanation"] for b in unit_wedges)
    first, second = (chemistry.draw(i["smiles"]) for i in _isomers(smiles))
    styles = lambda d: sorted((b["to_atom"], b["style"]) for b in d["depicted_stereo_bonds"] if b.get("unit"))
    assert styles(first) != styles(second)  # the two enantiomers look different


def test_helicene_drawing_has_label_not_wedges():
    for iso in _isomers(HELICENE6):
        drawn = chemistry.draw(iso["smiles"])
        assert drawn["depicted_stereo_bonds"] == [] and drawn["analysis"]["stereocenters"] == []
        assert drawn["depiction"]["method"] == "projection" and drawn["depiction"]["overlapping_atoms"] == 0
        assert "P or M" in drawn["depiction"]["note"] and "CIP_Code" in drawn["svg"]


@pytest.mark.parametrize("fmt, coordinates, keeps", [("smiles", None, True), ("sdf", "3d", True), ("mol", "3d", True),
                                                     ("mol", "2d", False), ("cdxml", "2d", False)])
def test_exports_keep_or_warn(fmt, coordinates, keeps):
    for smiles, kind in ((DICHLOROALLENE, "allene"), (HELICENE6, "helix")):
        for iso in _isomers(smiles):
            want = _units(iso["smiles"], kind)[0]["descriptor"]
            result = export_structure(iso["smiles"], fmt, coordinates)
            if not keeps:
                assert any("cannot carry" in w for w in result["warnings"])
                continue
            text = result["text"]
            mol = Chem.MolFromSmiles(text.split("\t")[0].strip()) if fmt == "smiles" else \
                Chem.RemoveHs(Chem.MolFromMolBlock(text.split("$$$$")[0], removeHs=False))
            cx = Chem.MolToCXSmiles(mol)
            assert _units(cx, kind)[0]["descriptor"] == want


def test_xaabb_spiro_survives_3d_export():
    for iso in _isomers(DICHLOROSPIROHEPTANE):
        want = _units(iso["smiles"], "spiro")[0]["descriptor"]
        text = export_structure(iso["smiles"], "sdf", "3d")["text"]
        back = Chem.RemoveHs(Chem.MolFromMolBlock(text.split("$$$$")[0], removeHs=False))
        assert _units(Chem.MolToCXSmiles(back), "spiro")[0]["descriptor"] == want


@pytest.mark.parametrize("name", ["(P)-hexahelicene", "(M)-hexahelicene", "(Ra)-1,3-dichloropropa-1,2-diene"])
def test_helical_and_axial_names_fail_closed(name):
    with pytest.raises(ValueError):
        asyncio.run(resolve_name(name))


# ---- crowded-layout panel -----------------------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", [
    "c1ccc(cc1)P(c2ccccc2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7ccccc7)c8ccccc8",  # BINAP
    "Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12",  # BINOL
    "Cc1cc2c(C(C)C)c(O)c(O)c(C=O)c2c(O)c1-c1c(C)cc2c(C(C)C)c(O)c(O)c(C=O)c2c1O",  # gossypol
    "CC1(C)c2cccc(P(c3ccccc3)c3ccccc3)c2Oc2c(P(c3ccccc3)c3ccccc3)cccc21",  # Xantphos
    "c1ccc(-c2c3ccccc3c(-c3ccccc3)c3c(-c4ccccc4)c4ccccc4c(-c4ccccc4)c23)cc1",  # rubrene
    "c1ccc(-c2c(-c3ccccc3)c(-c3ccccc3)c(-c3ccccc3)c(-c3ccccc3)c2-c2ccccc2)cc1",  # hexaphenylbenzene
    "c1ccc(C(c2ccccc2)(c2ccccc2)c2ccccc2)cc1", "c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    "[Cl-].[Rh].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    "[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1."
    "c1ccc(P(c2ccccc2)c2ccccc2)cc1",
    HELICENE5, HELICENE6,
])
def test_crowded_panel_is_drawn_cleanly(smiles):
    quality = layout(Chem.MolFromSmiles(smiles))
    assert quality["bond_crossings"] == 0 and quality["overlapping_atoms"] == 0, quality
    assert quality["max_bond_length_ratio"] <= 1.5, quality


def test_binap_open_and_specified_get_the_same_clean_layout():
    binap = "c1ccc(cc1)P(c2ccccc2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7ccccc7)c8ccccc8"
    for smiles in [binap] + [i["smiles"] for i in _isomers(binap)]:
        depiction = chemistry.draw(smiles)["depiction"]
        assert depiction["max_bond_length_ratio"] <= 1.5 and not depiction["bond_crossings"], depiction
