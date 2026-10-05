"""Bundled compound library and the name-resolution chain (library -> OPSIN -> PubChem), all offline here."""
import asyncio
import shutil
import pytest
from rdkit import Chem
from moltalk import library, naming
from moltalk.server import resolve_name

needs_opsin = pytest.mark.skipif(not (shutil.which("java") and naming.OPSIN_JAR.exists()), reason="needs Java and OPSIN")


async def _run(fn, *args):
    return fn(*args)


def test_library_is_broad():
    assert library.size() >= 5000


@pytest.mark.parametrize("name, formula", [
    ("aspirin", "C9H8O4"), ("Acetylsalicylic acid", "C9H8O4"), ("caffeine", "C8H10N4O2"), ("ibuprofen", "C13H18O2"),
    ("cholesterol", "C27H46O"), ("adenosine triphosphate", "C10H16N5O13P3"), ("serotonin", "C10H12N2O"),
    ("benzene", "C6H6"), ("ethanol", "C2H6O"), ("L-phenylalanine", "C9H11NO2"), ("limonene", "C10H16"),
])
def test_common_names_resolve_offline(name, formula):
    result = asyncio.run(resolve_name(name))
    assert result["source"].startswith("MolTalk library") and result["analysis"]["formula"] == formula
    assert result["cid"] and result["pubchem_url"].endswith(str(result["cid"]))


@needs_opsin
def test_systematic_names_go_to_opsin_with_stereo():
    result = asyncio.run(resolve_name("(2R)-2-chloro-4-methylhexane"))  # unlikely to be in the library
    assert result["source"].startswith("OPSIN")
    mol = Chem.MolFromSmiles(result["canonical_smiles"])
    from rdkit.Chem import rdCIPLabeler
    rdCIPLabeler.AssignCIPLabels(mol)
    assert sorted(a.GetProp("_CIPCode") for a in mol.GetAtoms() if a.HasProp("_CIPCode")) == ["R"]
    assert result["stereo_summary"]["unspecified"] == 1 and "unspecified" in result["warning"]  # C4 is left open


def test_unknown_names_fail_closed():
    with pytest.raises(ValueError, match="No structure was assumed"):
        asyncio.run(resolve_name("florbanex quintophane"))


def test_ambiguous_names_are_reported_not_guessed():
    data = library._load()
    ambiguous = next(n for n, ids in data["tiers"][1].items()  # e.g. "lye": NaOH or KOH
                     if len({data["compounds"][i]["inchikey"] for i in ids}) > 1 and n not in data["tiers"][0])
    with pytest.raises(ValueError, match="more than one structure"):
        asyncio.run(resolve_name(ambiguous))


def test_structure_to_name_and_locants_from_the_library(monkeypatch):
    monkeypatch.setattr(naming, "OFFLINE", True)  # no PubChem: the library alone must name it
    ibuprofen = asyncio.run(naming.name_and_locants("CC(C)Cc1ccc(cc1)C(C)C(=O)O", _run))
    assert ibuprofen["iupac_name"] == "2-[4-(2-methylpropyl)phenyl]propanoic acid"
    assert ibuprofen["name_source"].startswith("MolTalk library")
    assert ibuprofen["locants"] == {"10": "2", "11": "3", "12": "1"}
    # The same structure written in a different atom order gets the same numbering on the right atoms.
    other_order = asyncio.run(naming.name_and_locants("OC(=O)C(C)c1ccc(CC(C)C)cc1", _run))
    assert other_order["locants"] == {"1": "1", "3": "2", "4": "3"}


def test_novel_structures_get_no_name(monkeypatch):
    monkeypatch.setattr(naming, "OFFLINE", True)
    novel = asyncio.run(naming.name_and_locants("CC(C)(F)C1CC(Br)C(O)C1C#N", _run))
    assert novel["iupac_name"] is None and not novel["locants"]
    assert "library" in novel["name_status"]


@pytest.mark.parametrize("name", ["alanine", "serine", "glucose", "fructose"])
def test_bare_stereo_family_names_are_not_resolved_to_one_isomer(name):
    with pytest.raises(ValueError, match="does not say which stereoisomer"):
        asyncio.run(resolve_name(name))


@pytest.mark.parametrize("name, cip", [("L-alanine", "S"), ("D-alanine", "R")])
def test_qualified_amino_acids_resolve(name, cip):
    result = asyncio.run(resolve_name(name))
    assert [c["cip"] for c in result["analysis"]["stereocenters"]] == [cip]


def test_named_natural_products_keep_their_enantiomer():
    assert asyncio.run(resolve_name("morphine"))["source"].startswith("MolTalk library")  # (+)-morphine is not "morphine"


@pytest.mark.parametrize("name, cid, title, formula", [
    ("heme", 4973, "Heme b", "C34H32FeN4O4"), ("haem", 4973, "Heme b", "C34H32FeN4O4"),
    ("hemin", 455658, "Hemin", "C34H32ClFeN4O4"), ("NAD", 5892, "NAD+", "C21H27N7O14P2"),
    ("NAD+", 5892, "NAD+", "C21H27N7O14P2"), ("NADH", 439153, "NADH", "C21H29N7O14P2"),
    ("FADH2", 446013, "FADH2", "C27H35N9O15P2"), ("CoA", 87642, "Coenzyme A", "C21H36N7O16P3S"),
])
def test_curated_cofactors(name, cid, title, formula):
    result = asyncio.run(resolve_name(name))
    assert result["cid"] == cid and result["title"].startswith(title) and result["analysis"]["formula"] == formula
    assert result["note"]


def _draw(smiles, label):
    from moltalk.server import draw_molecule
    return asyncio.run(draw_molecule(smiles, label=label)).structuredContent["label_check"]


def test_draw_label_is_checked_against_the_library(monkeypatch):
    monkeypatch.setattr(naming, "OFFLINE", True)
    assert _draw("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "ibuprofen")["status"] == "matches"
    wrong = _draw("Cn1cnc2c1c(=O)[nH]c(=O)n2C", "caffeine")  # theobromine drawn under caffeine's name
    assert wrong["status"] == "mismatch" and wrong["library_cid"] == 2519
    assert _draw("C[C@@H](C(=O)O)N", "D-alanine")["status"] == "stereo differs"  # this SMILES is L-alanine
    assert _draw("CC(C)(F)C1CC(Br)C(O)C1C#N", "my compound")["status"] == "unverified"
