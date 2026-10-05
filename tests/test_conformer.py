import numpy as np
from moltalk.chemistry import analyze, conformer, draw
from test_depiction import C60


def _bond_lengths(result):
    xyz = np.array([a["xyz"] for a in result["atoms"]])
    return [np.linalg.norm(xyz[b["a"]] - xyz[b["b"]]) for b in result["bonds"]]


def test_heavy_atoms_in_input_order_with_labels():
    result = conformer("F[C@H]1C[C@@H](C)CCC1")
    assert [a["index"] for a in result["atoms"]] == list(range(8))  # hidden hydrogens stay hidden
    assert {a["index"]: a["cip"] for a in result["atoms"] if "cip" in a} == {1: "R", 3: "S"}
    assert result["atoms"][0]["element"] == "F" and result["atoms"][1]["h"] == 1
    assert "warning" not in result and "one calculated conformer" in result["note"]


def test_alignment_agrees_with_flat_wedges():
    """Lifting must start consistent with the flat drawing: wedged substituents in front, dashed behind."""
    for smiles in ["F[C@H]1C[C@@H](C)CCC1", "C[C@@H](O)C(=O)O", "O[C@@H]1CCCC[C@H]1O", "N[C@@H](Cc1ccccc1)C(=O)O",
                   "C[C@H](N)C(=O)N[C@@H](CC(C)C)C(=O)O", "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O",
                   "CN1CC[C@]23c4c5ccc(O)c4O[C@H]2[C@@H](O)C=C[C@H]3[C@H]1C5"]:
        result = conformer(smiles)
        assert result["wedges_contradicted"] == 0, smiles
        z = {a["index"]: a["xyz"][2] for a in result["atoms"]}
        for b in draw(smiles)["depicted_stereo_bonds"]:
            if b["to_atom"] not in z:
                continue  # a hydrogen RDKit added only for the flat drawing
            in_front = z[b["to_atom"]] > z[b["from_atom"]]
            assert in_front == (b["style"] == "wedge"), (smiles, b)


def test_ring_double_bonds_know_their_ring():
    result = conformer("CC(=O)Oc1ccccc1C(=O)O")
    ring_doubles = [b for b in result["bonds"] if b["order"] == 2 and b["ring"] is not None]
    assert len(ring_doubles) == 3 and all(len(result["rings"][b["ring"]]) == 6 for b in ring_doubles)
    assert sum(b["order"] == 2 and b["ring"] is None for b in result["bonds"]) == 2  # the two C=O


def test_arbitrary_stereo_is_flagged():
    result = conformer("CC(O)C(=O)O")
    assert result["arbitrary_stereo_atoms"] == [1]
    assert result["atoms"][1]["arbitrary"] is True
    assert "arbitrary" in conformer("CC=CC")["warning"]


def test_c60_gets_a_sensible_cage():
    result = conformer(C60)
    assert result["method"].startswith("spectral cage")
    lengths = _bond_lengths(result)
    assert len(lengths) == 90 and max(lengths) / min(lengths) < 1.2


def test_adamantane_bridgeheads_are_not_stereocentres():
    assert analyze("C1C2CC3CC1CC(C2)C3")["stereo_summary"]["status"] == "no stereo elements"
    assert "warning" not in conformer("C1C2CC3CC1CC(C2)C3")
    assert analyze("CC1CCC(C)CC1")["stereo_summary"]["unspecified"] == 2  # cis/trans is real


ERYTHROMYCIN = "CC[C@@H]1[C@@]([C@@H]([C@H](C(=O)[C@@H](C[C@@]([C@@H]([C@H]([C@@H]([C@H](C(=O)O1)C)O[C@H]2C[C@@]([C@H]([C@@H](O2)C)O)(C)OC)C)O[C@H]3[C@@H]([C@H](C[C@H](O3)C)N(C)C)O)(C)O)C)C)O)(C)O"


def test_cage_carbons_are_not_stereocentres():
    """Inverting one corner of a cubane cannot be built in 3D, although the SMILES would allow it."""
    for smiles in ("CC12C3C4C1C5C2C3C45C", "ClC12C3C4C1C5C2(Cl)C3C45", "CC12C3C4C1C5C2C3C45"):
        assert analyze(smiles)["stereo_summary"]["status"] == "no stereo elements", smiles
        model = conformer(smiles)
        assert not any("cip" in a for a in model["atoms"]) and "warning" not in model
    # Real stereo in bridged rings stays: norbornan-2-ol (C1, C2, C4) and camphor (C1, C4).
    assert analyze("OC1CC2CCC1C2")["stereo_summary"]["unspecified"] == 3
    assert analyze("CC1(C)C2CCC1(C)C(=O)C2")["stereo_summary"]["unspecified"] == 2


def test_large_macrolide_gets_a_3d_model_quickly():
    import time
    start = time.monotonic()
    model = conformer(ERYTHROMYCIN)
    assert time.monotonic() - start < 10
    assert sum("cip" in a for a in model["atoms"]) == 18
