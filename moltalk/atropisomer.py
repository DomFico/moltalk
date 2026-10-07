"""Axial chirality (atropisomers): a hindered single bond between two rings, as in BINAP, BINOL or a tetra-ortho-
substituted biaryl. Nothing here is specific to one compound.

Where the twist comes from: plain SMILES cannot state it. It arrives as CXSMILES with 2D coordinates and a wedge at an
axis atom ("|(...),wU:7.6|"), which RDKit reads (also with its legacy stereo perception, which MolTalk keeps for
everything else) into the axis bond's stereo, STEREOATROPCW/CCW. canonical_smiles() writes it back the same way, so
the twist survives every round trip; ETKDG embeds it; the drawing wedges a ring bond at the axis.

Descriptors are shown only when verified. RDKit's CIP labeller names the axis P or M; a 3D model of the molecule is
then measured independently: the helicity rule (sign of the torsion between the highest-ranked group on each side,
positive = P) and the CIP axial rule (seen along the axis, near groups 1 -> 2 then the far group 3; clockwise = Ra)
must agree with each other (Ra = M, Sa = P) and with RDKit, or the axis gets no descriptor at all.

Unspecified axes (BINAP written as plain SMILES) are found with a stated heuristic: a non-ring single bond between two
sp2 ring atoms, each with two ring neighbours that rank differently, and at least three of the four ortho positions
substituted or ring-fused (all four when an axis atom is in a five-membered ring, whose wider exocyclic angles hold
ortho groups further apart). Such biaryls are configurationally stable at room temperature; fewer ortho groups usually
rotate freely, so they are not reported. This is a rule of thumb, stated as such in the results, not a barrier
calculation.
"""
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdDepictor, rdMolTransforms

from .stereo import assign_cip

ATROP = (Chem.BondStereo.STEREOATROPCW, Chem.BondStereo.STEREOATROPCCW)
FLIP = {Chem.BondStereo.STEREOATROPCW: Chem.BondStereo.STEREOATROPCCW,
        Chem.BondStereo.STEREOATROPCCW: Chem.BondStereo.STEREOATROPCW}
MIN_ORTHO = 3
RANK_DEPTH = 8
_verified: dict[str, list] = {}


def specified_axes(mol) -> list:
    return [b for b in mol.GetBonds() if b.GetStereo() in ATROP]


def _kekule(mol):
    k = Chem.Mol(mol)
    with rdBase.BlockLogs():
        try:
            Chem.Kekulize(k, clearAromaticFlags=True)
        except Exception:  # noqa: BLE001 - rank on the aromatic form; aromatic bonds then count as single
            pass
    return k


def _branch(k, start: int, came_from: int, depth: int = RANK_DEPTH) -> tuple:
    """A simple CIP-style key for the branch start <- came_from: atomic numbers sphere by sphere (descending, with
    duplicate atoms for multiple bonds, which end there). Independent of RDKit's labeller; used only to cross-check
    it, so where this simplification and the full CIP rules differ, the descriptor is withheld, never guessed."""
    shells, frontier = [], [(start, came_from, None)]  # (atom, previous atom, atomic number if a duplicate)
    for _ in range(depth):
        shells.append(tuple(sorted((z if i is None else k.GetAtomWithIdx(i).GetAtomicNum() for i, _, z in frontier),
                                   reverse=True)))
        nxt = []
        for i, prev, _ in frontier:
            if i is None:
                continue
            atom = k.GetAtomWithIdx(i)
            for bond in atom.GetBonds():
                j = bond.GetOtherAtomIdx(i)
                order = 1 if bond.GetBondType() == Chem.BondType.AROMATIC else int(bond.GetBondTypeAsDouble())
                if j != prev:
                    nxt.append((j, i, None))
                nxt += [(None, i, k.GetAtomWithIdx(j).GetAtomicNum())] * (order - 1)
            nxt += [(None, i, 1)] * atom.GetTotalNumHs()
        frontier = nxt
    return tuple(shells)


def _ranked(k, atom: int, other: int) -> list[tuple[int, tuple]]:
    """The other neighbours of an axis atom, highest-ranked first, with their keys."""
    neighbours = [n.GetIdx() for n in k.GetAtomWithIdx(atom).GetNeighbors() if n.GetIdx() != other and n.GetAtomicNum() > 1]
    return sorted(((n, _branch(k, n, atom)) for n in neighbours), key=lambda item: item[1], reverse=True)


def _ortho_substituted(mol, axis_atom: int, ortho: int) -> bool:
    """An ortho atom carries a substituent or a fused ring (BINAP's C8a counts: its peri H blocks rotation)."""
    ring_atoms = set()
    for ring in mol.GetRingInfo().AtomRings():
        if axis_atom in ring and ortho in ring:
            ring_atoms |= set(ring)
    return any(n.GetAtomicNum() > 1 and n.GetIdx() != axis_atom and n.GetIdx() not in ring_atoms
               for n in mol.GetAtomWithIdx(ortho).GetNeighbors())


def candidate_axes(mol) -> list[int]:
    """Bond indices of likely stereogenic axes the input leaves unspecified (see the module docstring)."""
    k = _kekule(mol)
    found = []
    for bond in mol.GetBonds():
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing() or bond.GetStereo() in ATROP:
            continue
        ends = (bond.GetBeginAtom(), bond.GetEndAtom())
        if not all(a.IsInRing() and (a.GetIsAromatic() or a.GetHybridization() == Chem.HybridizationType.SP2) for a in ends):
            continue
        ortho_count, ok = 0, True
        for atom, other in (ends, ends[::-1]):
            ortho = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() != other.GetIdx() and n.GetAtomicNum() > 1]
            if len(ortho) != 2 or not all(mol.GetBondBetweenAtoms(atom.GetIdx(), o).IsInRing() for o in ortho):
                ok = False
                break
            ranked = _ranked(k, atom.GetIdx(), other.GetIdx())
            if ranked[0][1] == ranked[1][1]:  # e.g. 2,6-dimethylphenyl: the two sides are the same, no chirality
                ok = False
                break
            ortho_count += sum(_ortho_substituted(mol, atom.GetIdx(), o) for o in ortho)
        five = any(mol.GetRingInfo().IsAtomInRingOfSize(a.GetIdx(), 5) and not mol.GetRingInfo().IsAtomInRingOfSize(a.GetIdx(), 6)
                   for a in ends)
        if ok and ortho_count >= (4 if five else MIN_ORTHO):
            found.append(bond.GetIdx())
    return found


def from_geometry(mol, conf, bond) -> tuple[str, str] | None:
    """(helicity P/M, CIP axial Ra/Sa) measured on 3D coordinates, or None when the ranking ties or the axis is
    eclipsed. The two rules are computed separately; callers compare them."""
    k = _kekule(mol)
    a, c = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
    near, far = _ranked(k, a, c), _ranked(k, c, a)
    if len(near) < 2 or len(far) < 2 or near[0][1] == near[1][1] or far[0][1] == far[1][1]:
        return None
    n1, n2, f1 = near[0][0], near[1][0], far[0][0]
    torsion = rdMolTransforms.GetDihedralDeg(conf, n1, a, c, f1)
    if abs(torsion) < 5 or abs(torsion) > 175:
        return None
    helicity = "P" if torsion > 0 else "M"
    p = {i: np.array(conf.GetAtomPosition(i)) for i in (n1, n2, f1, a, c)}
    # Viewer at the near end looking along the axis (direction u): the turn n1 -> n2 -> f1 is clockwise when its
    # normal points away from the viewer, along u.
    turn = float(np.dot(np.cross(p[n2] - p[n1], p[f1] - p[n2]), p[c] - p[a]))
    if abs(turn) < 1e-6:
        return None
    return helicity, ("Ra" if turn > 0 else "Sa")


def _embed(mol):
    h = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 0xA7
    params.timeout = 10
    with rdBase.BlockLogs():
        if AllChem.EmbedMolecule(h, params) != 0:
            return None
    return h


def stereo_key(mol) -> str:
    """Identity including axial stereo: plain canonical SMILES drops it, so append each axis's twist, located by
    canonical atom ranks (the same molecule written in another atom order gives the same key)."""
    smiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    axes = specified_axes(mol)
    if not axes:
        return smiles
    labelled = Chem.Mol(mol)
    assign_cip(labelled, atoms=[], bonds=[b.GetIdx() for b in axes])
    ranks = list(Chem.CanonicalRankAtoms(labelled, breakTies=False))
    parts = sorted(f"{min(ranks[b.GetBeginAtomIdx()], ranks[b.GetEndAtomIdx()])}-"
                   f"{max(ranks[b.GetBeginAtomIdx()], ranks[b.GetEndAtomIdx()])}:"
                   f"{labelled.GetBondWithIdx(b.GetIdx()).GetPropsAsDict().get('_CIPCode', '?')}" for b in axes)
    return smiles + " axial " + ",".join(parts)


def canonical_smiles(mol) -> str:
    """Canonical SMILES; with a specified axis, CXSMILES carrying 2D coordinates and the axis wedge (the only way
    SMILES can state a twist). Self-checked: the written text must read back with the same descriptors."""
    if not specified_axes(mol):
        return Chem.MolToSmiles(mol, isomericSmiles=True)
    order = [int(i) for i in np.argsort(list(Chem.CanonicalRankAtoms(mol)))]
    candidates = [Chem.RenumberAtoms(mol, order), Chem.Mol(mol)]  # canonical order first: deterministic coordinates
    flags = Chem.CXSmilesFields.CX_COORDS | Chem.CXSmilesFields.CX_BOND_ATROPISOMER
    want = stereo_key(mol)
    for candidate in candidates:
        rdDepictor.Compute2DCoords(candidate)
        text = Chem.MolToCXSmiles(candidate, Chem.SmilesWriteParams(), flags)
        back = Chem.MolFromSmiles(text)
        if back is not None and stereo_key(back) == want:
            return text
    raise ValueError("The axial stereochemistry could not be written as CXSMILES; it was not dropped silently.")


def describe(mol) -> list[dict]:
    """Specified axes with verified descriptors; an axis that fails verification gets none (verified False)."""
    axes = specified_axes(mol)
    if not axes:
        return []
    key = stereo_key(mol) + "|" + ",".join(f"{b.GetBeginAtomIdx()}-{b.GetEndAtomIdx()}" for b in axes)
    if key in _verified:
        return [dict(item) for item in _verified[key]]
    labelled = Chem.Mol(mol)
    assign_cip(labelled, atoms=[], bonds=[b.GetIdx() for b in axes])
    model = _embed(mol)
    found = []
    for b in axes:
        rdkit_label = labelled.GetBondWithIdx(b.GetIdx()).GetPropsAsDict().get("_CIPCode")
        measured = from_geometry(model, model.GetConformer(), model.GetBondWithIdx(b.GetIdx())) if model else None
        entry = {"bond_index": b.GetIdx(), "atom_indices": [b.GetBeginAtomIdx(), b.GetEndAtomIdx()],
                 "helicity": None, "cip": None, "verified": False}
        if measured and rdkit_label == measured[0] and {"P": "Sa", "M": "Ra"}[measured[0]] == measured[1]:
            entry.update(helicity=measured[0], cip=measured[1], verified=True)
        else:
            entry["note"] = ("The axis is specified, but its descriptor could not be verified independently "
                             "(RDKit's label, the helicity rule and the CIP axial rule did not all agree), so none is shown.")
        found.append(entry)
    _verified[key] = found
    return [dict(item) for item in found]


def label(mol) -> list[dict]:
    """Put the verified axial descriptor (Ra/Sa) on each axis bond's _CIPCode, for the drawing's stereo annotation,
    and remove RDKit's unverified P/M. Returns describe(mol)."""
    found = describe(mol)
    for entry in found:
        bond = mol.GetBondWithIdx(entry["bond_index"])
        if entry["verified"]:
            bond.SetProp("_CIPCode", entry["cip"])
        elif bond.HasProp("_CIPCode"):
            bond.ClearProp("_CIPCode")
    return found


def with_axes(mol, assignment: dict[int, Chem.BondStereo]):
    """A copy with the given axis bonds set to a twist (for enumerating atropisomers)."""
    out = Chem.Mol(mol)
    for bond_index, stereo in assignment.items():
        out.GetBondWithIdx(bond_index).SetStereo(stereo)
    return out


def mirror_axes(mol):
    """Swap every axis's twist (part of taking the mirror image)."""
    for bond in specified_axes(mol):
        bond.SetStereo(FLIP[bond.GetStereo()])
    return mol


def enforce(mol, conf_id: int, axes: list[dict]) -> bool:
    """Make a 3D model's twist match the verified helicity, rotating one side about the axis if a later step (tail
    swing, chain straightening) turned it over. True if anything changed."""
    changed = False
    conf = mol.GetConformer(conf_id)
    k = _kekule(mol)
    for entry in axes:
        if not entry.get("verified"):
            continue
        a, c = entry["atom_indices"]
        bond = mol.GetBondBetweenAtoms(a, c)
        measured = from_geometry(mol, conf, bond)
        if measured and measured[0] != entry["helicity"]:
            n1, f1 = _ranked(k, a, c)[0][0], _ranked(k, c, a)[0][0]
            torsion = rdMolTransforms.GetDihedralDeg(conf, n1, a, c, f1)
            rdMolTransforms.SetDihedralDeg(conf, n1, a, c, f1, -torsion)
            changed = True
    return changed
