"""IUPAC names and parent locants. OPSIN runs locally (Java); PubChem only in the network test."""
import asyncio
import os
import shutil
import pytest
from moltalk import naming

needs_opsin = pytest.mark.skipif(not (shutil.which("java") and naming.OPSIN_JAR.exists()), reason="needs Java and the OPSIN jar")

# (SMILES, PubChem's IUPAC name, expected {input atom index: locant} on the parent)
CASES = {
    "3-methylhexane": ("CCCC(C)CC", "3-methylhexane", {"3": "3", "6": "1", "0": "6"}),
    "ibuprofen": ("CC(C)Cc1ccc(cc1)C(C)C(=O)O", "2-[4-(2-methylpropyl)phenyl]propanoic acid", {"12": "1", "10": "2", "11": "3"}),
    "LSD": ("CCN(CC)C(=O)[C@@H]1C=C2c3cccc4[nH]cc(c34)C[C@H]2N(C)C1",
            "(6aR,9R)-N,N-diethyl-7-methyl-6,6a,8,9-tetrahydro-4H-indolo[4,3-fg]quinoline-9-carboxamide",
            {"7": "9", "20": "6a", "21": "7", "15": "4"}),
    "testosterone": ("C[C@]12CC[C@H]3[C@H]([C@@H]1CC[C@@H]2O)CCC4=CC(=O)CC[C@]34C",
                     "(8R,9S,10R,13S,14S,17S)-17-hydroxy-10,13-dimethyl-1,2,6,7,8,9,11,12,14,15,16,17-dodecahydrocyclopenta[a]phenanthren-3-one",
                     {"1": "13", "9": "17", "15": "3", "19": "10"}),
    "glucose": ("OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "(3R,4S,5S,6R)-6-(hydroxymethyl)oxane-2,3,4,5-tetrol", {"3": "1", "4": "2", "2": "6"}),
    "fluoromethylcyclohexane": ("F[C@H]1C[C@@H](C)CCC1", "cis-(1R,3S)-1-fluoro-3-methylcyclohexane", {"1": "1", "3": "3"}),
}


@needs_opsin
@pytest.mark.parametrize("case", CASES)
def test_parent_locants(case):
    smiles, name, expected = CASES[case]
    cml = asyncio.run(naming.opsin_cml(name)) or asyncio.run(naming.opsin_cml(next(naming._variants(name), name).replace("cis-", "")))
    locants, reason = naming.locants_from_cml(smiles, name, cml)
    assert locants, reason
    for atom, locant in expected.items():
        assert locants.get(atom) == locant, (case, atom, locants)
    assert "4" not in locants or case != "3-methylhexane" or locants["4"] != "3"  # the methyl carbon is not numbered


@needs_opsin
def test_refuses_when_the_parent_is_not_determinable():
    smiles, name = "CC(=O)Nc1ccc(O)cc1", "N-(4-hydroxyphenyl)acetamide"
    locants, reason = naming.locants_from_cml(smiles, name, asyncio.run(naming.opsin_cml(name)))
    assert locants is None and "could not be identified" in reason


@needs_opsin
def test_refuses_a_name_for_a_different_structure():
    cml = asyncio.run(naming.opsin_cml("2-methylbutane"))
    assert not naming.matches_structure("CCCCC", cml)


def test_offline_mode_gives_no_name(monkeypatch):
    monkeypatch.setattr(naming, "OFFLINE", True)
    result = asyncio.run(naming.name_and_locants("CCO", lambda f, *a: asyncio.sleep(0, f(*a))))
    assert result["iupac_name"] is None and "turned off" in result["name_status"] and result["locants"] == {}


@pytest.mark.network
@pytest.mark.skipif(os.getenv("MOLTALK_NETWORK_TESTS") != "1", reason="set MOLTALK_NETWORK_TESTS=1")
def test_live_pubchem_name_and_locants(monkeypatch):
    monkeypatch.setattr(naming, "OFFLINE", False)

    async def run(fn, *args):
        return fn(*args)
    lsd = asyncio.run(naming.name_and_locants(CASES["LSD"][0], run))
    assert lsd["iupac_name"].startswith("(6aR,9R)") and lsd["locants"]["20"] == "6a"
    unknown = asyncio.run(naming.name_and_locants("CC(C)(F)C1CC(Br)C(O)C1C#N", run))
    assert unknown["iupac_name"] is None and "not in PubChem" in unknown["name_status"]
    phe = asyncio.run(naming.name_and_locants("N[C@@H](Cc1ccccc1)C(=O)O", run))  # zwitterion duplicate in PubChem is ignored
    assert phe["iupac_name"] == "(2S)-2-amino-3-phenylpropanoic acid"
