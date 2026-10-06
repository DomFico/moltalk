"""2D layout with a quality check and a Schlegel-style fallback for cages.

RDKit's fused-ring depiction can fail badly on polyhedral cages (fullerenes,
dodecahedrane): bonds stretch across the drawing and cross each other although
the graph is planar. Such layouts are scored, and a Tutte (barycentric)
embedding of the ring core is tried instead. For a 3-connected planar core,
Tutte's theorem guarantees a crossing-free drawing: a Schlegel diagram.
"""
import re
from functools import lru_cache

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


@lru_cache(maxsize=1)
def _porphyrinoid_core():
    """(query, coordinates): the porphine skeleton (four N-rings joined by four one-carbon bridges) as a pattern that
    ignores bond orders, aromaticity and hydrogens, with its textbook square layout from RDKit's ring templates.
    Hydroporphyrins such as coenzyme F430 match it, though their extra fused rings defeat RDKit's own template."""
    porphine = Chem.MolFromSmiles("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3")
    rdDepictor.Compute2DCoords(porphine, useRingTemplates=True)
    query = Chem.MolFromSmarts(re.sub(r"[:=\-]", "~", Chem.MolToSmarts(porphine)).replace("&H1", ""))
    match = porphine.GetSubstructMatch(query)
    conf = porphine.GetConformer()
    return query, np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in match])


def _shape_rmsd(a, b) -> float:
    """RMSD after the best similarity transform (translation, rotation or reflection, uniform scale) of a onto b."""
    a, b = a - a.mean(axis=0), b - b.mean(axis=0)
    u, sv, vt = np.linalg.svd(a.T @ b)
    scale = sv.sum() / max(1e-9, (a ** 2).sum())
    return float(np.sqrt(((scale * a @ u @ vt - b) ** 2).sum(axis=1).mean()))


def _pin_porphyrinoid(ligand) -> bool:
    """If the ligand has a porphyrin-type core that RDKit laid out badly (its fused rings folded inward), lay it out
    again with the core pinned to the textbook porphine square; RDKit then places fused rings and side chains."""
    query, template = _porphyrinoid_core()
    match = ligand.GetSubstructMatch(query)
    if not match:
        return False
    conf = ligand.GetConformer()
    current = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in match])
    if _shape_rmsd(current, template) < 0.25:  # already the textbook square (heme, chlorophyll)
        return False
    coords = {a: template[k] for k, a in enumerate(match)}
    # A ring fused to the core (F430's lactam along one bond; chlorophyll's ring E along ring C and a bridge carbon)
    # is not placed well by RDKit once the core is fixed. Put its other atoms on an arc of equal bonds between the ends
    # of the shared path, on the far side from that path (for a single shared bond, from the core ring it borders).
    rings = [list(r) for r in Chem.GetSymmSSSR(ligand)]
    core = set(match)
    bond = float(np.median([np.linalg.norm(coords[b.GetBeginAtomIdx()] - coords[b.GetEndAtomIdx()])
                            for b in ligand.GetBonds() if b.GetBeginAtomIdx() in core and b.GetEndAtomIdx() in core]))
    for ring in rings:
        shared = [a for a in ring if a in core]
        free = [a for a in ring if a not in core]
        if len(shared) < 2 or not free or len(ring) > 8 or any(a in coords for a in free):
            continue
        k = next(i for i, a in enumerate(ring) if a in core and ring[(i + 1) % len(ring)] not in core)
        ring = ring[k:] + ring[:k]  # ring[0]: last shared atom before the free run; then the free atoms in order
        run = [a for a in ring[1:] if a not in core]
        if ring[1 + len(run):] != [a for a in ring if a in core][1:]:
            continue  # shared atoms are not one contiguous path
        p0, p1 = coords[ring[0]], coords[ring[len(run) + 1]]
        if len(shared) == 2:
            host = next((r for r in rings if set(shared) <= set(r) and set(r) <= core), None)
            if host is None:
                continue
            away_from = np.mean([coords[a] for a in host], axis=0)
        else:
            away_from = np.mean([coords[a] for a in shared], axis=0)
        mid, chord = (p0 + p1) / 2, p1 - p0
        span = float(np.linalg.norm(chord))
        normal = np.array([-chord[1], chord[0]]) / max(span, 1e-9)
        if np.dot(normal, mid - away_from) < 0:
            normal = -normal
        steps = len(run) + 1
        # Circle through p0 and p1 bulging along `normal`, whose arc splits into `steps` chords of length `bond`.
        lo, hi = span / 2 + 1e-6, 50 * bond
        for _ in range(60):
            r = (lo + hi) / 2
            half = np.arcsin(min(1.0, span / (2 * r)))
            theta = (2 * np.pi - 2 * half) / steps  # take the long (outer) arc
            lo, hi = (r, hi) if 2 * r * np.sin(theta / 2) < bond else (lo, r)
        half = np.arcsin(min(1.0, span / (2 * r)))
        centre = mid + normal * r * np.cos(half)  # the long arc: the centre is on the bulge side
        a0 = np.arctan2(*(p0 - centre)[::-1])
        a1 = np.arctan2(*(p1 - centre)[::-1])
        sweep = (a1 - a0) % (2 * np.pi)
        direction = -1 if sweep < np.pi else 1  # go the long way round, through the outer side
        total = sweep if direction == 1 else 2 * np.pi - sweep
        for j, a in enumerate(run, start=1):
            angle = a0 + direction * total * j / steps
            coords[a] = centre + r * np.array([np.cos(angle), np.sin(angle)])
    # Substituents on the core and its fused rings point out of the macrocycle (its middle is where the metal and its
    # bonds go): put a placed atom's one unplaced, acyclic neighbour in its widest outward-facing gap.
    middle = template.mean(axis=0)
    for a in list(coords):  # core atoms and the atoms of fused rings placed above
        nbrs = [n.GetIdx() for n in ligand.GetAtomWithIdx(a).GetNeighbors()]
        free = [n for n in nbrs if n not in coords]
        if len(free) != 1 or ligand.GetAtomWithIdx(free[0]).IsInRing():
            continue
        angles = sorted(np.arctan2(*(coords[n] - coords[a])[::-1]) for n in nbrs if n in coords)
        outward = coords[a] - middle
        best = None
        for k, start in enumerate(angles):
            end = angles[(k + 1) % len(angles)] + (2 * np.pi if k + 1 == len(angles) else 0)
            mid = (start + end) / 2
            direction = np.array([np.cos(mid), np.sin(mid)])
            key = (np.dot(direction, outward) > 0, end - start)
            if best is None or key > best[0]:
                best = (key, direction)
        coords[free[0]] = coords[a] + bond * best[1]
    trial = Chem.Mol(ligand)
    try:
        with rdBase.BlockLogs():
            rdDepictor.Compute2DCoords(trial, coordMap={a: Point2D(*xy) for a, xy in coords.items()})
    except Exception:
        return False
    # The folded layout can score well on these counts while looking wrong, so the textbook core wins unless it
    # adds crossings, overlaps or stretched bonds (a slightly longer longest bond does not count).
    if _score(quality(trial))[:2] > _score(quality(ligand))[:2]:
        return False
    ligand.RemoveAllConformers()
    ligand.AddConformer(Chem.Conformer(trial.GetConformer()), assignId=True)
    return True


def layout(mol) -> dict:
    """Compute 2D coordinates in place; returns how they were made and how good they are.
    With explicit hydrogens, the heavy-atom skeleton is laid out first (with every fix below) and the hydrogens are
    then placed around it, so the drawing with hydrogens has the same skeleton as the one without. Laying out all
    atoms at once gave a different, sometimes overlapping layout, which then wrongly triggered the cage fallback."""
    hydrogens = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    if hydrogens and len(hydrogens) < mol.GetNumAtoms():
        heavy_ids = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
        skeleton = Chem.RWMol(mol)
        for i in sorted(hydrogens, reverse=True):
            skeleton.RemoveAtom(i)
        skeleton = skeleton.GetMol()
        skeleton.UpdatePropertyCache(strict=False)
        Chem.GetSymmSSSR(skeleton)  # removing atoms clears the ring information the layout code relies on
        result = layout(skeleton)
        conf = skeleton.GetConformer()
        coord_map = {heavy_ids[k]: Point2D(conf.GetAtomPosition(k).x, conf.GetAtomPosition(k).y) for k in range(len(heavy_ids))}
        with rdBase.BlockLogs():
            rdDepictor.Compute2DCoords(mol, coordMap=coord_map)
        placed = quality(mol)
        # RDKit's own all-at-once layout is sometimes tidier for small molecules (glucose's ring hydrogens); keep it
        # only when it has strictly fewer crossings, overlaps and stretched bonds (and the skeleton needed no fix-up).
        if result["method"] == "rdkit" and not _chelated_metal(skeleton):
            whole = Chem.Mol(mol)
            with rdBase.BlockLogs():
                rdDepictor.Compute2DCoords(whole)
            alone = quality(whole)
            problems = lambda q: q["bond_crossings"] + q["overlapping_atoms"] + q["stretched_bonds"]
            if problems(alone) < problems(placed):
                mol.RemoveAllConformers()
                mol.AddConformer(Chem.Conformer(whole.GetConformer()), assignId=True)
                placed = alone
        return {**result, **placed}
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
        _pin_porphyrinoid(ligand)
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
    # Schlegel diagrams are for polyhedral cages (C60, cubane, dodecahedrane), whose atoms sit in three or more rings.
    # A merely crowded ordinary layout (FAD with hydrogens) must never be swapped for one.
    Chem.GetSymmSSSR(mol)  # make sure ring information exists (it may have been cleared by editing the molecule)
    ring_info = mol.GetRingInfo()
    if not any(ring_info.NumAtomRings(a.GetIdx()) >= 3 for a in mol.GetAtoms()):
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
