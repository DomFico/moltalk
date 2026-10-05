from rdkit import Chem
from moltalk.chemistry import analyze, conformer, draw

# PubChem's forms: a charged macrocycle plus a separate metal ion.
HEME = "CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C(=C3CCC(=O)O)C)CCC(=O)O.[Fe+2]"
ZNPC = "C1=CC=C2C(=C1)C3=NC2=NC4=C5C=CC=CC5=C([N-]4)N=C6C7=CC=CC=C7C(=N6)N=C8C9=CC=CC=C9C(=N8)[N-]3.[Zn+2]"


def test_heme_drawn_with_fe_n_bonds_but_analysis_unchanged():
    result = draw(HEME)
    note = result["depiction"]["coordination_note"]
    assert "2 covalent and 2 dative" in note and "unchanged" in note
    assert result["analysis"]["formula"] == "C34H32FeN4O4"
    assert result["analysis"]["properties"]["formal_charge"] == 0
    assert "disconnected fragments" in " ".join(result["analysis"]["warnings"])  # source structure, as given
    model = conformer(HEME)
    assert sum(b.get("dative", False) for b in model["bonds"]) == 2
    fe = next(a for a in model["atoms"] if a["element"] == "Fe")
    assert fe["charge"] == 0


def test_charges_kept_when_they_do_not_line_up():
    note = draw(ZNPC)["depiction"]["coordination_note"]
    assert "4 dative" in note


def test_ordinary_molecules_untouched():
    assert "coordination_note" not in draw("CC(=O)Oc1ccccc1C(=O)O")["depiction"]
    assert "coordination_note" not in draw("[Na+].[Cl-]")["depiction"]


def test_explicit_hydrogens_match_between_flat_and_3d():
    smiles = "F[C@H]1C[C@@H](C)CCC1"
    flat = draw(smiles, hydrogens=True)
    model = conformer(smiles, hydrogens=True)
    total = Chem.AddHs(Chem.MolFromSmiles(smiles)).GetNumAtoms()
    assert len(flat["atom_px"]) == len(model["atoms"]) == total
    assert model["hydrogens_shown"] and all(a["h"] == 0 for a in model["atoms"])
    z = {a["index"]: a["xyz"][2] for a in model["atoms"]}
    for b in flat["depicted_stereo_bonds"]:
        assert (z[b["to_atom"]] > z[b["from_atom"]]) == (b["style"] == "wedge")
    assert analyze(smiles)["formula"] == "C7H13F"


# How ChatGPT once wrote heme itself: Fe bonded to all four N, charge-separated ([Fe-2] with two [N+]).
HEME_CHARGE_SEPARATED = ("C=CC1=C(C)C2=Cc3c(C=C)c(C)c4[n]3[Fe-2]35[n]6c(c(C)c(CCC(=O)O)c6=CC6=[N+]3C(=C4)C(C)=C6CCC(=O)O)"
                         "=CC1=[N+]25")


def test_metal_is_never_a_stereocentre():
    for smiles in (HEME, HEME_CHARGE_SEPARATED):
        analysis = analyze(smiles)
        assert analysis["stereo_summary"]["status"] == "no stereo elements" and not analysis["stereocenters"]
        model = conformer(smiles)
        fe = next(a for a in model["atoms"] if a["element"] == "Fe")
        assert "cip" not in fe and "warning" not in model


def test_charge_separated_form_drawn_neutral_with_flat_ring():
    import numpy as np
    from moltalk.depiction import _cage
    note = draw(HEME_CHARGE_SEPARATED)["depiction"]["coordination_note"]
    assert "dative" in note and "unchanged" in note
    model = conformer(HEME_CHARGE_SEPARATED)
    assert all(a["charge"] == 0 for a in model["atoms"] if a["element"] in ("Fe", "N"))
    for smiles in (HEME, HEME_CHARGE_SEPARATED):
        model = conformer(smiles)
        core = _cage(Chem.MolFromSmiles(smiles))
        xyz = np.array([a["xyz"] for a in model["atoms"]])[core]
        bend = np.linalg.svd(xyz - xyz.mean(0), compute_uv=False)[2] / np.sqrt(len(core))
        assert bend < 0.2, (smiles, bend)  # the porphyrin ring is essentially flat
