"""Axial chirality (atropisomers): generic detection, verified Ra/Sa, round trips, drawing, 3D, and no change to
ordinary R/S and E/Z."""
import asyncio

import pytest
from rdkit import Chem
from rdkit.Chem import rdChemDraw

from moltalk import atropisomer, chemistry
from moltalk.export import export_structure
from moltalk.server import resolve_name

BINAP = "c1ccc(cc1)P(c2ccccc2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7ccccc7)c8ccccc8"
BINOL = "Oc1ccc2ccccc2c1-c1c(O)ccc2ccccc12"
TETRA_ORTHO = "Cc1cccc(Cl)c1-c1c(C)cccc1Cl"  # no fused rings, no phosphorus: nothing BINAP-specific


def _isomers(smiles):
    return chemistry.enumerate_stereo(smiles)["isomers"]


@pytest.mark.parametrize("smiles", [BINAP, BINOL, TETRA_ORTHO, "OC(=O)c1cccc(C)c1-c1c(C)cccc1[N+](=O)[O-]",
                                    "COc1ccc2ccccc2c1-c1c(C)cccc1Br",
                                    "Cc1cc2c(C(C)C)c(O)c(O)c(C=O)c2c(O)c1-c1c(C)cc2c(C(C)C)c(O)c(O)c(C=O)c2c1O"])  # gossypol
def test_hindered_axes_are_found(smiles):
    assert len(atropisomer.candidate_axes(Chem.MolFromSmiles(smiles))) == 1
    assert chemistry.analyze(smiles)["stereo_summary"]["unspecified_axes"]


@pytest.mark.parametrize("smiles", ["c1ccccc1-c1ccccc1",  # free rotation
                                    "Cc1ccccc1-c1ccccc1C",  # two ortho groups: rotates at room temperature
                                    "Cc1cccc(C)c1-c1c(C)cccc1Cl",  # 2,6-dimethylphenyl: both sides alike, no chirality
                                    "OC(=O)c1ccccc1-c1ccccc1C(=O)O", "CC(C)c1ccccc1",
                                    "Cc1cccc(-c2nn3c(c2-c2ccnc4ccc(C(N)=O)cc24)CCC3)n1"])  # galunisertib: 5-ring, 3 ortho
def test_unhindered_or_symmetric_axes_are_not(smiles):
    assert not atropisomer.candidate_axes(Chem.MolFromSmiles(smiles))
    assert chemistry.analyze(smiles)["stereo_summary"]["unspecified_axes"] == []


@pytest.mark.parametrize("smiles", [BINAP, BINOL, TETRA_ORTHO])
def test_both_atropisomers_with_verified_descriptors(smiles):
    isomers = _isomers(smiles)
    assert len(isomers) == 2
    axial = [iso["axial_stereo"][0] for iso in isomers]
    assert all(x["verified"] for x in axial)
    assert {(x["helicity"], x["cip"]) for x in axial} == {("M", "Ra"), ("P", "Sa")}  # IUPAC: Ra = M, Sa = P
    assert isomers[0]["enantiomer_index"] == 1 and isomers[1]["enantiomer_index"] == 0
    assert not any(iso["achiral"] for iso in isomers)
    for iso in isomers:  # the CXSMILES carries the twist: reading it back gives the same descriptor
        again = chemistry.analyze(iso["smiles"])
        assert again["axial_stereo"][0]["cip"] == iso["axial_stereo"][0]["cip"]
        assert again["canonical_smiles"] == iso["smiles"]
        assert again["stereo_summary"]["status"] == "fully specified"


def test_canonical_form_does_not_depend_on_atom_order():
    ra = next(iso["smiles"] for iso in _isomers(BINOL) if iso["axial_stereo"][0]["cip"] == "Ra")
    mol = chemistry.parse(ra)
    shuffled = Chem.RenumberAtoms(mol, list(reversed(range(mol.GetNumAtoms()))))
    assert atropisomer.stereo_key(shuffled) == atropisomer.stereo_key(mol)
    assert chemistry.analyze(atropisomer.canonical_smiles(shuffled))["axial_stereo"][0]["cip"] == "Ra"


def test_stereocentre_and_axis_enumerate_together():
    isomers = _isomers("C[C@@H](O)c1cccc(C)c1-c1c(C)cccc1Cl")  # one fixed centre plus one open axis
    assert len(isomers) == 2  # diastereomers: the centre stays (R), the axis takes both twists
    assert {iso["axial_stereo"][0]["cip"] for iso in isomers} == {"Ra", "Sa"}
    assert all(iso["enantiomer_index"] is None and [c["cip"] for c in iso["stereocenters"]] == ["R"] for iso in isomers)


@pytest.mark.parametrize("cip", ["Ra", "Sa"])
def test_drawing_wedges_the_axis(cip):
    smiles = next(iso["smiles"] for iso in _isomers(BINAP) if iso["axial_stereo"][0]["cip"] == cip)
    drawn = chemistry.draw(smiles)
    axis_wedges = [b for b in drawn["depicted_stereo_bonds"] if "axis_atoms" in b]
    assert len(axis_wedges) == 1 and axis_wedges[0]["axial_cip"] == cip
    assert "chirality axis" in axis_wedges[0]["explanation"] and "not a stereocentre" in axis_wedges[0]["explanation"]
    assert "CIP_Code" in drawn["svg"]  # RDKit's own (Ra)/(Sa) annotation, toggled by the viewer's Stereo control


@pytest.mark.parametrize("smiles", [BINAP, TETRA_ORTHO])
def test_3d_model_has_the_twist(smiles):
    for iso in _isomers(smiles):
        model = chemistry.conformer(iso["smiles"])
        mol = chemistry.parse(iso["smiles"])
        conf = Chem.Conformer(mol.GetNumAtoms())
        for atom in model["atoms"][:mol.GetNumAtoms()]:
            conf.SetAtomPosition(atom["index"], atom["xyz"])
        mol.RemoveAllConformers()
        mol.AddConformer(conf)
        axis = atropisomer.specified_axes(mol)[0]
        assert atropisomer.from_geometry(mol, mol.GetConformer(), axis) == (iso["axial_stereo"][0]["helicity"],
                                                                            iso["axial_stereo"][0]["cip"])
        tagged = [b for b in model["bonds"] if b.get("cip") in ("Ra", "Sa")]
        assert [b["cip"] for b in tagged] == [iso["axial_stereo"][0]["cip"]] and not tagged[0].get("arbitrary")
        assert "warning" not in model


def test_3d_twist_turned_over_is_put_back():
    iso = _isomers(TETRA_ORTHO)[0]
    mol = chemistry.parse(iso["smiles"])
    model = atropisomer._embed(mol)
    a, c = iso["axial_stereo"][0]["atom_indices"]
    k = atropisomer._kekule(model)
    n1, f1 = atropisomer._ranked(k, a, c)[0][0], atropisomer._ranked(k, c, a)[0][0]
    from rdkit.Chem import rdMolTransforms
    conf = model.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, n1, a, c, f1, -rdMolTransforms.GetDihedralDeg(conf, n1, a, c, f1))
    assert atropisomer.from_geometry(model, conf, model.GetBondBetweenAtoms(a, c))[1] != iso["axial_stereo"][0]["cip"]
    assert atropisomer.enforce(model, conf.GetId(), iso["axial_stereo"])
    assert atropisomer.from_geometry(model, conf, model.GetBondBetweenAtoms(a, c))[1] == iso["axial_stereo"][0]["cip"]


def test_open_axis_is_arbitrary_in_3d():
    model = chemistry.conformer(BINAP)
    tagged = [b for b in model["bonds"] if b.get("cip") in ("Ra", "Sa")]
    assert len(tagged) == 1 and tagged[0]["arbitrary"]
    assert "atropisomer twist" in model["warning"]
    assert any("atropisomer twist" in w for w in export_structure(BINAP, "sdf", "3d")["warnings"])


def test_unverified_descriptor_is_withheld(monkeypatch):
    smiles = _isomers(TETRA_ORTHO)[0]["smiles"]
    monkeypatch.setattr(atropisomer, "_verified", {})
    monkeypatch.setattr(atropisomer, "from_geometry", lambda *args: ("P", "Ra"))  # rules disagree: Ra must be M
    analysis = chemistry.analyze(smiles)
    assert analysis["axial_stereo"][0]["verified"] is False and analysis["axial_stereo"][0]["cip"] is None
    mol = chemistry.parse(smiles)
    assert not any(b.HasProp("_CIPCode") for b in atropisomer.specified_axes(mol))  # RDKit's raw P/M is not shown
    monkeypatch.setattr(atropisomer, "_verified", {})


@pytest.mark.parametrize("fmt, coordinates", [("mol", "2d"), ("sdf", "2d"), ("mol", "3d"), ("sdf", "3d"),
                                              ("smiles", None), ("cdxml", "2d")])
def test_exports_keep_the_twist(fmt, coordinates):
    for iso in _isomers(BINOL):
        text = export_structure(iso["smiles"], fmt, coordinates)["text"]
        if fmt == "smiles":
            mol = Chem.MolFromSmiles(text.split("\t")[0].strip())
        elif fmt == "cdxml":
            mol = rdChemDraw.MolsFromChemDrawBlock(text)[0]
        else:
            mol = Chem.MolFromMolBlock(text.split("$$$$")[0], removeHs=False)
        assert [x["cip"] for x in atropisomer.describe(mol)] == [iso["axial_stereo"][0]["cip"]]


@pytest.mark.parametrize("name", ["(R)-BINAP", "(S)-BINAP", "(Ra)-BINAP", "(R)-BINOL", "(S)-BINOL", "(P)-BINAP"])
def test_axial_names_still_fail_closed(name):
    with pytest.raises(ValueError, match="No structure was assumed|could not be verified"):
        asyncio.run(resolve_name(name))


@pytest.mark.parametrize("name", ["BINAP", "BINOL"])
def test_plain_names_resolve_with_the_axis_open(name):
    result = asyncio.run(resolve_name(name))
    assert result["stereo_summary"]["unspecified_axes"] and result["stereo_summary"]["status"] == "unspecified"
    assert "|" not in result["canonical_smiles"]  # no twist invented


@pytest.mark.parametrize("smiles, canonical, centres, doubles", [
    ("C[C@H](O)CC", "CC[C@H](C)O", ["S"], []),
    ("N[C@@H](C)C(=O)O", "C[C@H](N)C(=O)O", ["S"], []),
    ("C/C=C/C", "C/C=C/C", [], ["E"]),
    ("C/C=C\\C", "C/C=C\\C", [], ["Z"]),
    ("CC1=CC[C@@H](CC1=O)C(C)=C", "C=C(C)[C@H]1CC=C(C)C(=O)C1", ["S"], []),  # (S)-carvone
    ("C[C@@H]1CC[C@H](C)CC1", None, ["s", "s"], []),  # pseudo-asymmetric, as before
    ("c1ccccc1-c1ccccc1", "c1ccc(-c2ccccc2)cc1", [], []),
])
def test_ordinary_stereo_is_unchanged(smiles, canonical, centres, doubles):
    analysis = chemistry.analyze(smiles)
    mol = Chem.MolFromSmiles(smiles)
    assert analysis["canonical_smiles"] == (canonical or Chem.MolToSmiles(mol))  # plain SMILES, never CXSMILES
    assert [c["cip"] for c in analysis["stereocenters"]] == centres
    assert [d["configuration"] for d in analysis["double_bond_stereo"]] == doubles
    assert analysis["axial_stereo"] == [] and analysis["stereo_summary"]["unspecified_axes"] == []


def test_ordinary_enumeration_is_unchanged():
    isomers = _isomers("CC(O)C(O)C")  # 2,3-butanediol: (R,R), (S,S) and meso
    assert len(isomers) == 3 and sum(iso["achiral"] for iso in isomers) == 1
    assert all(iso["axial_stereo"] == [] and "|" not in iso["smiles"] for iso in isomers)
    assert len(_isomers("CC=CC")) == 2
