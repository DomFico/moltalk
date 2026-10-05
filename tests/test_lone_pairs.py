import pytest
from rdkit import Chem
from moltalk.chemistry import depiction_mol, electrons

# (SMILES, {element+charge sign: (lone pairs, unpaired electrons)}) for every atom that should get dots.
CASES = {
    "water": ("O", {"O": (2, 0)}),
    "hydroxide": ("[OH-]", {"O-": (3, 0)}),
    "hydronium": ("[OH3+]", {"O+": (1, 0)}),
    "ammonium": ("[NH4+]", {}),
    "amide anion": ("[NH2-]", {"N-": (2, 0)}),
    "carbanion": ("[CH3-]", {"C-": (1, 0)}),
    "carbocation": ("[CH3+]", {}),
    "methyl radical": ("[CH3]", {"C": (0, 1)}),
    "chloride": ("[Na+].[Cl-]", {"Cl-": (4, 0)}),
    "acetate": ("CC(=O)[O-]", {"O": (2, 0), "O-": (3, 0)}),
    "nitromethane": ("C[N+](=O)[O-]", {"O": (2, 0), "O-": (3, 0)}),
    "protonated acetone": ("CC(=[OH+])C", {"O+": (1, 0)}),
    "iminium": ("C=[NH2+]", {}),
    "pyridine": ("c1ccncc1", {"N": (1, 0)}),
    "pyrrole": ("c1cc[nH]c1", {"N": (1, 0)}),
    "carbon monoxide": ("[C-]#[O+]", {"C-": (1, 0), "O+": (1, 0)}),
    "azide anion": ("[N-]=[N+]=[N-]", {"N-": (2, 0)}),
    "nitroxide radical": ("CC1(C)CCCC(C)(C)N1[O]", {"N": (1, 0), "O": (2, 1)}),
}


@pytest.mark.parametrize("case", CASES)
def test_lone_pairs_and_radicals(case):
    smiles, expected = CASES[case]
    mol = Chem.MolFromSmiles(smiles)
    found = {}
    for idx, info in electrons(mol).items():
        atom = mol.GetAtomWithIdx(int(idx))
        key = atom.GetSymbol() + {1: "+", -1: "-"}.get(atom.GetFormalCharge(), "")
        assert found.setdefault(key, (info["pairs"], info["radicals"])) == (info["pairs"], info["radicals"]), case
    assert found == expected, (case, found)


def test_dative_donor_has_given_its_pair_away():
    heme, _ = depiction_mol("CC1=C(C2=CC3=NC(=CC4=C(C(=C([N-]4)C=C5C(=C(C(=N5)C=C1[N-]2)C)C=C)C)C=C)C(=C3CCC(=O)O)C)CCC(=O)O.[Fe+2]")
    found = electrons(heme)
    donors = [b.GetBeginAtomIdx() for b in heme.GetBonds() if b.GetBondType() == Chem.BondType.DATIVE]
    assert len(donors) == 2 and all(str(d) not in found for d in donors)
    assert not any(heme.GetAtomWithIdx(int(i)).GetSymbol() == "Fe" for i in found)  # metals are skipped


import numpy as np
from moltalk.chemistry import conformer


def _angle(a, b):
    return float(np.degrees(np.arccos(np.clip(np.dot(a, b) / np.linalg.norm(a) / np.linalg.norm(b), -1, 1))))


def _pairs(smiles, atom):
    model = conformer(smiles)
    xyz = {a["index"]: np.array(a["xyz"]) for a in model["atoms"]}
    mol = Chem.MolFromSmiles(smiles)
    bonds = [xyz[n.GetIdx()] - xyz[atom] for n in mol.GetAtomWithIdx(atom).GetNeighbors()]
    return [np.array(d) for d in model["lone_pair_dirs"].get(str(atom), [])], bonds, model


@pytest.mark.parametrize("smiles, atom, count, to_bonds, between", [
    ("O", 0, 2, None, 109.5),            # water: tetrahedral pairs
    ("CC(C)=O", 3, 2, 120, 120),         # ketone O (sp2): pairs at 120 deg in the plane
    ("c1cc[nH]c1", 3, 1, 90, None),      # pyrrole N: the pair is in the p orbital
    ("c1ccncc1", 3, 1, 121, None),       # pyridine N: in the ring plane, pointing out
    ("CC#N", 2, 1, 180, None),           # nitrile N (sp): opposite the triple bond
    ("CF", 1, 3, 109.5, 109.5),          # fluorine on sp3 carbon: tetrahedral
])
def test_lone_pair_directions_follow_vsepr(smiles, atom, count, to_bonds, between):
    dirs, bonds, _ = _pairs(smiles, atom)
    assert len(dirs) == count
    if to_bonds:
        assert all(abs(_angle(d, b) - to_bonds) < 6 for d in dirs for b in bonds), [_angle(d, b) for d in dirs for b in bonds]
    if between:
        assert all(abs(_angle(dirs[i], dirs[j]) - between) < 6 for i in range(count) for j in range(i + 1, count))


def test_expanded_octets_get_their_vsepr_shape():
    dirs, bonds, model = _pairs("F[Xe](F)(F)F", 1)  # square planar, pairs trans
    assert sorted(round(_angle(a, b)) for i, a in enumerate(bonds) for b in bonds[i + 1:]) == [90, 90, 90, 90, 180, 180]
    assert len(dirs) == 2 and abs(_angle(dirs[0], dirs[1]) - 180) < 1 and "VSEPR" in model["note"]
    dirs, bonds, _ = _pairs("FS(F)(F)F", 1)  # seesaw, one equatorial pair
    assert sorted(round(_angle(a, b)) for i, a in enumerate(bonds) for b in bonds[i + 1:]) == [90, 90, 90, 90, 120, 180]
    assert len(dirs) == 1
    dirs, bonds, _ = _pairs("F[Xe]F", 1)  # linear, three equatorial pairs
    assert round(_angle(bonds[0], bonds[1])) == 180 and len(dirs) == 3


def test_metal_centre_geometry_is_flagged_not_trusted():
    model = conformer("N->[Pt](<-N)(Cl)Cl")
    assert "generic force field" in model["note"] and "cis/trans" in model["note"]
    assert model["lone_pair_dirs"].keys() == {"3", "4"}  # the two Cl; the ammine N donated their pairs to Pt
