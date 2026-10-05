"""2D layout with a quality check and a Schlegel-style fallback for cages.

RDKit's fused-ring depiction can fail badly on polyhedral cages (fullerenes,
dodecahedrane): bonds stretch across the drawing and cross each other although
the graph is planar. Such layouts are scored, and a Tutte (barycentric)
embedding of the ring core is tried instead. For a 3-connected planar core,
Tutte's theorem guarantees a crossing-free drawing: a Schlegel diagram.
"""
import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor
from rdkit.Geometry import Point2D, Point3D

STRETCH = 3.0      # bonds longer than this × the median count as stretched
CLASH = 0.3        # non-bonded atoms closer than this × the median bond length overlap
ALPHAS = (1.0, 0.8, 0.65, 0.5, 0.4)  # radial expansion r**alpha to un-crowd inner rings
MAX_OUTER_FACES = 40


def _crossings(pos, bonds) -> int:
    """Number of pairs of bonds (sharing no atom) whose segments properly intersect."""
    if len(bonds) < 2:
        return 0
    b = np.asarray(bonds)
    p, q = pos[b[:, 0]], pos[b[:, 1]]
    i, j = np.triu_indices(len(b), 1)
    disjoint = ((b[i, 0] != b[j, 0]) & (b[i, 0] != b[j, 1]) & (b[i, 1] != b[j, 0]) & (b[i, 1] != b[j, 1]))
    i, j = i[disjoint], j[disjoint]

    def orient(a, c, d):
        return np.sign((c[:, 0] - a[:, 0]) * (d[:, 1] - a[:, 1]) - (c[:, 1] - a[:, 1]) * (d[:, 0] - a[:, 0]))

    hit = ((orient(p[i], q[i], p[j]) * orient(p[i], q[i], q[j]) < 0) &
           (orient(p[j], q[j], p[i]) * orient(p[j], q[j], q[i]) < 0))
    return int(hit.sum())


def quality(mol) -> dict:
    pos = mol.GetConformer().GetPositions()[:, :2]
    bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]
    if not bonds:
        return {"bond_crossings": 0, "stretched_bonds": 0, "overlapping_atoms": 0, "max_bond_length_ratio": 1.0}
    lengths = np.linalg.norm(pos[[a for a, _ in bonds]] - pos[[b for _, b in bonds]], axis=1)
    median = float(np.median(lengths)) or 1.0
    bonded = {frozenset(x) for x in bonds}
    d = np.linalg.norm(pos[:, None] - pos[None], axis=2)
    close = np.argwhere(np.triu(d < CLASH * median, 1))
    return {"bond_crossings": _crossings(pos, bonds),
            "stretched_bonds": int((lengths > STRETCH * median).sum()),
            "overlapping_atoms": sum(frozenset(map(int, x)) not in bonded for x in close),
            "max_bond_length_ratio": round(float(lengths.max() / median), 1)}


def _score(q: dict):
    return (q["bond_crossings"] + q["overlapping_atoms"], q["stretched_bonds"], q["max_bond_length_ratio"])


def _cage(mol) -> list[int]:
    """Atoms of the largest fused ring system (rings sharing atoms); the part Tutte can embed."""
    systems = []
    for ring in mol.GetRingInfo().AtomRings():
        ring = set(ring)
        for system in [x for x in systems if x & ring]:
            systems.remove(system)
            ring |= system
        systems.append(ring)
    return sorted(max(systems, key=len)) if systems else []


def _tutte(core, adjacency, outer):
    index = {a: k for k, a in enumerate(core)}
    pos = np.zeros((len(core), 2))
    for k, a in enumerate(outer):
        t = 2 * np.pi * k / len(outer)
        pos[index[a]] = (np.cos(t), np.sin(t))
    fixed = set(outer)
    free = [a for a in core if a not in fixed]
    if not free:
        return pos
    col = {a: k for k, a in enumerate(free)}
    lap = np.zeros((len(free), len(free)))
    rhs = np.zeros((len(free), 2))
    for a in free:
        neighbours = adjacency[a]
        lap[col[a], col[a]] = len(neighbours)
        for b in neighbours:
            if b in fixed:
                rhs[col[a]] += pos[index[b]]
            else:
                lap[col[a], col[b]] -= 1
    pos[[index[a] for a in free]] = np.linalg.solve(lap, rhs)
    return pos


def _schlegel_candidates(mol):
    core = _cage(mol)
    if len(core) < 4:
        return
    core_set = set(core)
    adjacency = {a: [n.GetIdx() for n in mol.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() in core_set] for a in core}
    rings = [r for r in mol.GetRingInfo().AtomRings() if set(r) <= core_set]
    faces = sorted(rings, key=len, reverse=True)[:MAX_OUTER_FACES]
    for outer in faces:
        try:
            base = _tutte(core, adjacency, list(outer))
        except np.linalg.LinAlgError:
            continue
        radius = np.linalg.norm(base, axis=1, keepdims=True)
        for alpha in ALPHAS:
            pos = base * np.where(radius > 0, np.power(radius, alpha - 1, where=radius > 0, out=np.ones_like(radius)), 0)
            yield core, pos


def _chelated_metal(mol):
    """(metal index, ligand atom indices) for a metal bonded only to four nitrogens, as in a drawn porphyrin."""
    from .coordination import METALS
    for atom in mol.GetAtoms():
        neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
        if atom.GetSymbol() in METALS and len(neighbours) == 4 and all(mol.GetAtomWithIdx(n).GetSymbol() == "N" for n in neighbours):
            return atom.GetIdx(), neighbours
    return None


def _turn_substituents_down(mol, donors):
    """Rotate by a multiple of 90° so the bulk of the side groups hangs below the macrocycle, as in the usual
    heme drawing (propionate chains at the bottom)."""
    conf = mol.GetConformer()
    pos = conf.GetPositions()
    centre = pos[donors].mean(axis=0)
    ring_atoms = set(_cage(mol))
    side = [a.GetIdx() for a in mol.GetAtoms() if a.GetIdx() not in ring_atoms and a.GetAtomicNum() > 1]
    if not side:
        return
    quarter = np.array([[0.0, -1, 0], [1, 0, 0], [0, 0, 1]])
    turns = [np.linalg.matrix_power(quarter, k) for k in range(4)]
    rot = min(turns, key=lambda r: (r @ (pos[side] - centre).T)[1].mean())  # lowest mean y = side groups at the bottom
    for i, p in enumerate(pos):
        conf.SetAtomPosition(i, Point3D(*(rot @ (p - centre) + centre)))


def layout(mol) -> dict:
    """Compute 2D coordinates in place; returns how they were made and how good they are."""
    chelate = _chelated_metal(mol)
    if chelate:
        # RDKit cannot fit the metal's four bonds into the macrocycle layout. Lay out the ligand alone (a tidy
        # macrocycle) and put the metal at the centre of its four nitrogens, as textbook drawings do.
        metal, donors = chelate
        ligand = Chem.RWMol(mol)
        for n in donors:
            ligand.RemoveBond(n, metal)
        ligand = ligand.GetMol()
        ligand.UpdatePropertyCache(strict=False)
        # RDKit's ring templates give the textbook square porphyrin (pyrroles on the four sides, N pointing in).
        rdDepictor.Compute2DCoords(ligand, useRingTemplates=True)
        _turn_substituents_down(ligand, donors)
        conf = Chem.Conformer(ligand.GetConformer())
        centre = np.mean([list(conf.GetAtomPosition(n)) for n in donors], axis=0)
        conf.SetAtomPosition(metal, Point3D(*centre))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        return {"method": "rdkit", **quality(mol)}
    rdDepictor.Compute2DCoords(mol)
    best_q = quality(mol)
    result = {"method": "rdkit", **best_q}
    if not (best_q["bond_crossings"] or best_q["stretched_bonds"] or best_q["overlapping_atoms"]):
        return result
    best_conf = Chem.Conformer(mol.GetConformer())
    core_bond_length = 1.5
    for core, pos in _schlegel_candidates(mol):
        trial = Chem.Mol(mol)
        # Scale so a typical core bond has RDKit's usual length, then let RDKit place substituents.
        lengths = [np.linalg.norm(pos[i] - pos[j]) for i, a in enumerate(core) for j, b in enumerate(core)
                   if i < j and trial.GetBondBetweenAtoms(a, b)]
        scale = core_bond_length / (float(np.median(lengths)) or 1.0)
        coord_map = {a: Point2D(*(pos[k] * scale)) for k, a in enumerate(core)}
        # Point each substituent's first atom away from the cage centre so RDKit grows it outward.
        centre = pos.mean(axis=0) * scale
        core_set = set(core)
        for k, a in enumerate(core):
            outside = [n.GetIdx() for n in trial.GetAtomWithIdx(a).GetNeighbors() if n.GetIdx() not in core_set]
            direction = pos[k] * scale - centre
            base_angle = np.arctan2(direction[1], direction[0])
            for m, b in enumerate(outside):
                if b in coord_map:
                    continue
                angle = base_angle + (m - (len(outside) - 1) / 2) * np.pi / 3
                coord_map[b] = Point2D(*(pos[k] * scale + core_bond_length * np.array([np.cos(angle), np.sin(angle)])))
        try:
            with rdBase.BlockLogs():  # rejected candidates can trip RDKit invariant messages
                rdDepictor.Compute2DCoords(trial, coordMap=coord_map)
        except Exception:
            continue
        q = quality(trial)
        if _score(q) < _score(best_q):
            best_q, best_conf = q, Chem.Conformer(trial.GetConformer())
            result = {"method": "schlegel", **q}
            if q["bond_crossings"] == 0 and q["overlapping_atoms"] == 0 and q["stretched_bonds"] == 0:
                break
    mol.RemoveAllConformers()
    mol.AddConformer(best_conf, assignId=True)
    return result
