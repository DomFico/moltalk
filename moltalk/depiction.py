"""2D layout with a quality check, textbook porphyrin layouts and 3D-view drawings for cages.

RDKit's fused-ring depiction fails on polyhedral cages (fullerenes, cubane, adamantane): bonds stretch across
the drawing and cross. Such molecules are drawn as a view of their 3D shape, as textbooks do.
"""
import re
from functools import lru_cache

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import rdDepictor, rdqueries
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


# ---- Layout energy and repair --------------------------------------------------------------------------------------
# A crowded molecule gets several candidate layouts (RDKit's, RDKit's with random flips, CoordGen's). Each is repaired
# by a small search over rigid moves and scored; the lowest energy wins. Energy terms, in units of the median bond
# length L: bond-length deviation beyond 10% (a 3.8x bond costs ~730, so a stretched layout never wins), atoms closer
# than L (and coincident atoms), an atom lying on a bond it is not part of, bond crossings, uneven bond angles around
# an atom (an exocyclic bond off the ring's bisector looks like a distorted ring) and ring polygons that are not
# regular. Moves change coordinates only, never connectivity; the chosen layout's wedges are then checked to encode
# the same stereochemistry as the input (see stereo_faithful).
REPAIR_BUDGET_S = 1.5
_TURNS = np.radians([0, 30, -30, 60, -60, 90, -90, 120, -120, 180])


def _energy_context(mol):
    bonds = np.array([(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()]).reshape(-1, 2)
    pos = mol.GetConformer().GetPositions()[:, :2]
    lengths = np.linalg.norm(pos[bonds[:, 0]] - pos[bonds[:, 1]], axis=1) if len(bonds) else np.ones(1)
    ordinary = [k for k, b in enumerate(mol.GetBonds()) if not _metal_bond(b)]
    length = float(np.median(lengths[ordinary] if ordinary else lengths)) or 1.0
    n = mol.GetNumAtoms()
    far = Chem.GetDistanceMatrix(mol) >= 3  # 1-2 and 1-3 pairs are close by construction
    centres = {}
    for atom in mol.GetAtoms():
        neighbours = [nb.GetIdx() for nb in atom.GetNeighbors()]
        if len(neighbours) >= 2:
            centres[atom.GetIdx()] = (neighbours, atom.GetHybridization() == Chem.HybridizationType.SP)
    incident = np.zeros((n, len(bonds)), dtype=bool)
    incident[bonds[:, 0], np.arange(len(bonds))] = incident[bonds[:, 1], np.arange(len(bonds))] = True
    rings = [list(r) for r in mol.GetRingInfo().AtomRings() if len(r) <= 8]
    # Bond-length targets: metal–ligand bonds may be drawn up to 2.5x (room for bulky ligands, as in drawings of
    # Pd(PPh3)4); every other bond wants the median length.
    stretch = np.array([2.5 if _metal_bond(b) else 1.0 for b in mol.GetBonds()])
    return {"bonds": bonds, "L": length, "far": far, "centres": centres, "incident": incident, "rings": rings,
            "stretch": stretch}


def _metal_bond(bond) -> bool:
    from .coordination import is_metal
    return is_metal(bond.GetBeginAtom()) or is_metal(bond.GetEndAtom())


def _bond_term(rel, stretch):
    over = np.where(stretch > 1, np.maximum(0, rel - stretch), np.abs(rel - 1) - 0.1)
    return 100 * np.sum(np.maximum(0, over) ** 2) + 100 * np.sum(np.maximum(0, 0.9 - rel) ** 2) * (stretch > 1).any()


def _crossing_pairs(pos, first, second) -> int:
    """Proper intersections between bonds in `first` and bonds in `second` (arrays of atom pairs) sharing no atom."""
    if not len(first) or not len(second):
        return 0
    i, j = np.meshgrid(np.arange(len(first)), np.arange(len(second)), indexing="ij")
    i, j = i.ravel(), j.ravel()
    f, g = first[i], second[j]
    disjoint = (f[:, 0] != g[:, 0]) & (f[:, 0] != g[:, 1]) & (f[:, 1] != g[:, 0]) & (f[:, 1] != g[:, 1])
    f, g = f[disjoint], g[disjoint]
    p, q, r, t = pos[f[:, 0]], pos[f[:, 1]], pos[g[:, 0]], pos[g[:, 1]]

    def orient(a, c, d):
        return np.sign((c[:, 0] - a[:, 0]) * (d[:, 1] - a[:, 1]) - (c[:, 1] - a[:, 1]) * (d[:, 0] - a[:, 0]))

    return int(((orient(p, q, r) * orient(p, q, t) < 0) & (orient(r, t, p) * orient(r, t, q) < 0)).sum())


def _angle_term(pos, ctx, atoms) -> float:
    energy = 0.0
    for i in atoms:
        if i not in ctx["centres"]:
            continue
        neighbours, linear = ctx["centres"][i]
        w = pos[neighbours] - pos[i]
        theta = np.sort(np.arctan2(w[:, 1], w[:, 0]))
        gaps = np.diff(np.append(theta, theta[0] + 2 * np.pi))
        if len(neighbours) == 2:
            energy += 10 * (gaps.min() - (np.pi if linear else 2 * np.pi / 3)) ** 2
        else:
            energy += 10 * np.sum((gaps - 2 * np.pi / len(neighbours)) ** 2)
    return energy


def _interaction(pos, ctx, inside, outside, bonds_in, bonds_out, bonds_all) -> float:
    """Non-bonded terms between two atom sets: crowding, atoms on bonds, crossings."""
    length = ctx["L"]
    energy = 0.0
    if len(inside) and len(outside):
        d = np.linalg.norm(pos[inside][:, None] - pos[outside][None], axis=2) / length
        mask = ctx["far"][np.ix_(inside, outside)]
        d = d[mask]
        energy += 50 * np.sum(d < CLASH) + 30 * np.sum(np.maximum(0, 1 - d) ** 2)
    bonds = ctx["bonds"]
    for atoms, bond_ids in ((inside, bonds_out), (outside, bonds_in)):
        if not len(atoms) or not len(bond_ids):
            continue
        p = pos[bonds[bond_ids, 0]]
        v = pos[bonds[bond_ids, 1]] - p
        t = np.clip(((pos[atoms][:, None] - p[None]) * v[None]).sum(2) / ((v * v).sum(1)[None] + 1e-12), 0, 1)
        dist = np.linalg.norm(pos[atoms][:, None] - (p[None] + t[..., None] * v[None]), axis=2) / length
        dist = dist[~ctx["incident"][np.ix_(atoms, bond_ids)]]
        energy += 120 * np.sum(np.maximum(0, 0.5 - dist) ** 2)
    energy += 30 * _crossing_pairs(pos, bonds[bonds_in], bonds[bonds_out])
    return energy


def layout_energy(pos, ctx) -> float:
    """Total energy (see the section comment)."""
    bonds, length = ctx["bonds"], ctx["L"]
    if not len(bonds):
        return 0.0
    n = len(pos)
    rel = np.linalg.norm(pos[bonds[:, 0]] - pos[bonds[:, 1]], axis=1) / length
    energy = _bond_term(rel, ctx["stretch"])
    d = np.linalg.norm(pos[:, None] - pos[None], axis=2) / length
    upper = np.triu(ctx["far"], 1)
    energy += 50 * np.sum(d[upper] < CLASH) + 30 * np.sum(np.maximum(0, 1 - d[upper]) ** 2)
    p, v = pos[bonds[:, 0]], pos[bonds[:, 1]] - pos[bonds[:, 0]]
    t = np.clip(((pos[:, None] - p[None]) * v[None]).sum(2) / ((v * v).sum(1)[None] + 1e-12), 0, 1)
    to_bond = np.linalg.norm(pos[:, None] - (p[None] + t[..., None] * v[None]), axis=2) / length
    energy += 120 * np.sum(np.maximum(0, 0.5 - to_bond[~ctx["incident"]]) ** 2)
    energy += 30 * _crossings(pos, [tuple(b) for b in bonds])
    energy += _angle_term(pos, ctx, range(n))
    for ring in ctx["rings"]:
        q = pos[ring]
        sides = np.linalg.norm(q - np.roll(q, -1, 0), axis=1) / length
        a, b = q - np.roll(q, 1, 0), np.roll(q, -1, 0) - q
        cos = -(a * b).sum(1) / (np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1) + 1e-9)
        angles = np.arccos(np.clip(cos, -1, 1))
        energy += 50 * np.sum((sides - 1) ** 2) + 20 * np.sum((angles - np.pi * (len(ring) - 2) / len(ring)) ** 2)
    return float(energy)


def _movable_sides(mol):
    """For each acyclic bond: (bond index, anchor, joint, atoms of the smaller side beyond the joint)."""
    sides = []
    for bond in mol.GetBonds():
        if bond.IsInRing():
            continue
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        cut = Chem.RWMol(mol)
        cut.RemoveBond(a, b)
        fragments = Chem.GetMolFrags(cut)
        side = next(f for f in fragments if b in f)
        if len(side) > mol.GetNumAtoms() / 2:
            side = next(f for f in fragments if a in f)
            a, b = b, a
        sides.append((bond.GetIdx(), a, b, np.array(side)))
    return sides


def repair_layout(mol, budget_s: float = REPAIR_BUDGET_S) -> float:
    """Greedy search, in place: for each acyclic bond, give it its target length and try the smaller side turned
    about the anchor atom and mirrored across the bond; take the move that lowers the energy most, and repeat.
    A rigid move changes only the terms between the moved side and the rest (and the angles at the bond's two
    atoms), so each trial is scored by that difference alone."""
    import time
    ctx = _energy_context(mol)
    pos = mol.GetConformer().GetPositions()[:, :2].copy()
    current = layout_energy(pos, ctx)
    n = mol.GetNumAtoms()
    bonds = ctx["bonds"]
    prepared = []
    for bond_index, a, b, side in _movable_sides(mol):
        inside_mask = np.zeros(n, dtype=bool)
        inside_mask[side] = True
        outside = np.flatnonzero(~inside_mask)
        in_b = inside_mask[bonds[:, 0]] & inside_mask[bonds[:, 1]]
        out_b = ~inside_mask[bonds[:, 0]] & ~inside_mask[bonds[:, 1]]
        bonds_in = np.flatnonzero(in_b | (~in_b & ~out_b))  # the cut bond moves with the side
        bonds_out = np.flatnonzero(out_b)
        prepared.append((bond_index, a, b, side, outside, bonds_in, bonds_out))

    def partial(trial, item):
        bond_index, a, b, side, outside, bonds_in, bonds_out = item
        rel = np.linalg.norm(trial[a] - trial[b]) / ctx["L"]
        return (_bond_term(np.array([rel]), ctx["stretch"][[bond_index]]) + _angle_term(trial, ctx, (a, b))
                + _interaction(trial, ctx, side, outside, bonds_in, bonds_out, None))

    started = time.monotonic()
    for _ in range(25):
        best = None
        for item in prepared:
            if time.monotonic() - started > budget_s:
                break
            bond_index, a, b, side = item[:4]
            limit = ctx["stretch"][bond_index]
            current_rel = np.linalg.norm(pos[b] - pos[a]) / ctx["L"]
            targets = {min(max(current_rel, 1.0), limit)} | ({1.0, 1.75, limit} if limit > 1 else {1.0})
            base = partial(pos, item)
            anchor = pos[a]
            bond_len = np.linalg.norm(pos[b] - anchor) or 1.0
            u = (pos[b] - anchor) / bond_len
            normal = np.array([-u[1], u[0]])
            shapes = []
            for target in sorted(targets):
                rel = pos[side] - anchor - u * (bond_len - target * ctx["L"])
                shapes += [rel, np.outer(rel @ u, u) - np.outer(rel @ normal, normal)]
            for shape in shapes:
                for turn in _TURNS:
                    c, s_ = np.cos(turn), np.sin(turn)
                    trial = pos.copy()
                    trial[side] = anchor + shape @ np.array([[c, s_], [-s_, c]])
                    delta = partial(trial, item) - base
                    if delta < -1e-6 and (best is None or delta < best[0]):
                        best = (delta, trial)
        if best is None:
            break
        pos = best[1]
        current += best[0]
    conf = mol.GetConformer()
    for i, (x, y) in enumerate(pos):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), 0.0))
    return layout_energy(pos, ctx)


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


def _views(count=120):
    """Evenly spread viewing directions over a hemisphere (a view and its opposite give mirror-image pictures)."""
    k = np.arange(count) + 0.5
    z = k / count  # 0..1: one hemisphere
    phi = np.pi * (1 + 5 ** 0.5) * k
    r = np.sqrt(1 - z * z)
    return np.column_stack([r * np.cos(phi), r * np.sin(phi), z])


def _project_cage(mol):
    """2D coordinates for a polyhedral cage (cubane, adamantane, dodecahedrane, C60) as a view of its 3D shape, the
    way such molecules are drawn in textbooks. A flat layout cannot draw them without crossings, and a Schlegel
    diagram puts inner atoms where their hydrogens have nowhere to go. The 3D model is the same one the rotated view
    uses; of 120 viewing directions the one with fewest coinciding atoms, then fewest crossing bonds, is kept.
    Shown hydrogens are projected from the same model, so they point outward. Sets the conformer in place."""
    from .conformer import _embed, _spectral_sphere
    heavy_ids = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    skeleton = Chem.RWMol(mol)
    for i in sorted((a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1), reverse=True):
        skeleton.RemoveAtom(i)
    skeleton = skeleton.GetMol()
    skeleton.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(skeleton)
    try:
        embedded, _ = _embed(skeleton, 1)
    except ValueError:
        # A substituted fullerene (PCBM): ETKDG cannot embed it. View the cage alone (from its spectral sphere) and
        # let RDKit lay out the rest around the fixed cage.
        cage = _cage(skeleton)
        sphere = _spectral_sphere(skeleton, cage)
        flat = _best_view(sphere, [(cage.index(b.GetBeginAtomIdx()), cage.index(b.GetEndAtomIdx()))
                                   for b in skeleton.GetBonds()
                                   if b.GetBeginAtomIdx() in cage and b.GetEndAtomIdx() in cage], list(range(len(cage))))
        coord_map = {heavy_ids[a]: Point2D(*flat[k]) for k, a in enumerate(cage)}
        with rdBase.BlockLogs():
            rdDepictor.Compute2DCoords(mol, coordMap=coord_map)
        return
    pos3 = embedded.GetConformer().GetPositions()
    # 3D position of every atom of mol: heavy atoms by order, hydrogens by rank on their parent.
    where = {}
    for k, i in enumerate(heavy_ids):
        where[i] = pos3[k]
    hydrogens = {}
    for atom in embedded.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1:
            hydrogens.setdefault(atom.GetNeighbors()[0].GetIdx(), []).append(atom.GetIdx())
    rank = {}
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1:
            parent = heavy_ids.index(atom.GetNeighbors()[0].GetIdx())
            r = rank.get(parent, 0)
            rank[parent] = r + 1
            spare = hydrogens.get(parent, [])
            where[atom.GetIdx()] = pos3[spare[r]] if r < len(spare) else pos3[parent]
    xyz = np.array([where[i] for i in range(mol.GetNumAtoms())])
    xyz -= xyz[heavy_ids].mean(axis=0)
    heavy_bonds = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx()) for b in mol.GetBonds()
                   if b.GetBeginAtom().GetAtomicNum() != 1 and b.GetEndAtom().GetAtomicNum() != 1]
    flat = _best_view(xyz, heavy_bonds, heavy_ids)
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (px, py) in enumerate(flat):
        conf.SetAtomPosition(i, Point3D(float(px), float(py), 0.0))
    mol.RemoveAllConformers()
    mol.AddConformer(conf, assignId=True)


def _best_view(xyz, heavy_bonds, heavy_ids):
    """Project 3D points (all atoms; heavy_ids are the heavy ones) along the viewing direction with fewest coinciding
    heavy atoms, then fewest coinciding hydrogens, then fewest crossing bonds; longest extent horizontal; scaled to
    RDKit's 1.5 bond length."""
    xyz = xyz - xyz[heavy_ids].mean(axis=0)
    bond = float(np.median([np.linalg.norm(xyz[a] - xyz[b]) for a, b in heavy_bonds])) or 1.0
    best = None
    for view in _views():
        z = view / np.linalg.norm(view)
        x = np.cross([0.0, 1.0, 0.0] if abs(z[1]) < 0.9 else [1.0, 0.0, 0.0], z)
        x /= np.linalg.norm(x)
        y = np.cross(z, x)
        flat = xyz @ np.column_stack([x, y])
        h = flat[heavy_ids]
        d = np.linalg.norm(h[:, None] - h[None], axis=2)
        coincide = int((np.triu(d < 0.35 * bond, 1)).sum())
        everything = np.linalg.norm(flat[:, None] - flat[None], axis=2)
        crowded = int((np.triu(everything < 0.3 * bond, 1)).sum()) - coincide if len(flat) > len(h) else 0
        key = (coincide, crowded, _crossings(flat, heavy_bonds), -float(np.min(d + np.eye(len(h)) * 1e9)))
        if best is None or key < best[0]:
            best = (key, flat)
    flat = best[1]
    h = flat[heavy_ids] - flat[heavy_ids].mean(axis=0)
    _, _, vt = np.linalg.svd(h)
    return (flat - flat[heavy_ids].mean(axis=0)) @ vt.T * (1.5 / bond)


def _coordgen(mol) -> bool:
    """CoordGen's layout (in place), at its best precision: it resolves crowding that RDKit's own engine cannot
    (BINAP: RDKit put both PPh2 groups on top of the naphthalenes, 9 crossings and 10 overlaps)."""
    from rdkit.Chem import rdCoordGen
    params = rdCoordGen.CoordGenParams()
    params.minimizerPrecision = params.sketcherBestPrecision
    try:
        with rdBase.BlockLogs():
            rdCoordGen.AddCoords(mol, params)
        return True
    except Exception:  # noqa: BLE001 - CoordGen cannot lay out some structures (metals); keep the other layout
        return False


def _shorten_bridges(mol):
    """CoordGen avoids overlaps partly by stretching a connecting bond (BINAP's biaryl bond 3.8x, Xantphos's C-P
    3.4x). For each stretched bond outside rings, bring the smaller side back to a normal bond length, trying turns
    about the joint, and keep the first placement that adds no crossings or overlaps."""
    conf = mol.GetConformer()
    for _ in range(4):
        pos = conf.GetPositions()[:, :2]
        lengths = {b.GetIdx(): np.linalg.norm(pos[b.GetBeginAtomIdx()] - pos[b.GetEndAtomIdx()]) for b in mol.GetBonds()}
        median = float(np.median(list(lengths.values()))) if lengths else 1.0
        stretched = [b for b in mol.GetBonds() if not b.IsInRing() and lengths[b.GetIdx()] > STRETCH * median]
        if not stretched:
            return
        bond = max(stretched, key=lambda b: lengths[b.GetIdx()])
        a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        side = {b}
        todo = [b]
        while todo:
            for n in mol.GetAtomWithIdx(todo.pop()).GetNeighbors():
                if n.GetIdx() not in side and not (todo == [] and n.GetIdx() == a and len(side) == 1) and n.GetIdx() != a:
                    side.add(n.GetIdx()); todo.append(n.GetIdx())
        if a in side:
            return
        if len(side) > mol.GetNumAtoms() / 2:  # move the smaller side
            side = set(range(mol.GetNumAtoms())) - side
            a, b = b, a
        before = quality(mol)
        moved = sorted(side)
        rel = pos[moved] - pos[b]
        direction = (pos[b] - pos[a]) / lengths[bond.GetIdx()]
        best = None
        # Mirror image of the side as well: a 2D drawing may be flipped (wedges are recomputed from the stereo).
        along = direction
        across = np.array([-along[1], along[0]])
        mirrored = np.column_stack([rel @ along, -(rel @ across)]) @ np.vstack([along, across])
        for flip, shape in ((False, rel), (True, mirrored)):
          for turn in np.radians(np.arange(0, 360, 15)):
            c, s_ = np.cos(turn), np.sin(turn)
            rot = np.array([[c, -s_], [s_, c]])
            d = rot @ direction
            new_b = pos[a] + median * d
            # Turn the side with its bond, so it keeps its shape relative to the bond direction.
            placed = new_b + shape @ rot.T
            trial = Chem.Mol(mol)
            tconf = trial.GetConformer()
            for i, p in zip(moved, placed):
                tconf.SetAtomPosition(int(i), Point3D(float(p[0]), float(p[1]), 0.0))
            q = quality(trial)
            if q["bond_crossings"] <= before["bond_crossings"] and q["overlapping_atoms"] <= before["overlapping_atoms"] \
                    and q["stretched_bonds"] < before["stretched_bonds"]:
                key = (q["bond_crossings"] + q["overlapping_atoms"], flip, min(turn, 2 * np.pi - turn))
                if best is None or key < best[0]:
                    best = (key, placed)
        if best is None:
            return
        for i, p in zip(moved, best[1]):
            conf.SetAtomPosition(int(i), Point3D(float(p[0]), float(p[1]), 0.0))


MAX_RATIO = 1.5  # a conventional drawing's longest bond; above this the candidate search runs


def _best_candidate(mol, current_q) -> dict:
    """Crowded (BINAP, Xantphos, rubrene): candidate layouts, each repaired; the lowest-energy one that encodes the
    input's stereochemistry replaces the current layout (in place). Returns its quality."""
    candidates = [Chem.Mol(mol)]
    sampled = Chem.Mol(mol)
    with rdBase.BlockLogs():
        rdDepictor.Compute2DCoords(sampled, nFlipsPerSample=2, nSample=100, sampleSeed=7, permuteDeg4Nodes=True)
    candidates.append(sampled)
    coordgen = Chem.Mol(mol)
    if _coordgen(coordgen):
        candidates.append(coordgen)
    hub = _hub_layout(mol)
    if hub is not None:
        candidates.append(hub)
    # Each candidate as made (CoordGen with its stretched bridges pulled in, as before) and repaired: repair lowers the
    # energy, which can trade a crossing for shorter bonds; the acceptance rule below decides.
    if len(candidates) > 2:
        _shorten_bridges(candidates[2])
    trials = []
    for candidate in candidates:
        raw = Chem.Mol(candidate)
        trials.append((layout_energy(raw.GetConformer().GetPositions()[:, :2], _energy_context(raw)), raw))
        # Big molecules (cyclic peptides, 100+ atoms) get a third of the time: CoordGen alone takes seconds there.
        budget = REPAIR_BUDGET_S if mol.GetNumAtoms() <= 60 else REPAIR_BUDGET_S / 3
        trials.append((repair_layout(candidate, budget / len(candidates)), candidate))
    best = None
    for energy, candidate in trials:
        if not stereo_faithful(candidate):
            continue  # never trade stereochemistry for looks
        q = quality(candidate)
        # Never more crossings or overlaps than the current layout (a bridged ring traded for a crossing is not an
        # improvement), and something must actually improve.
        if q["bond_crossings"] > current_q["bond_crossings"] or q["overlapping_atoms"] > current_q["overlapping_atoms"] \
                or _score(q) >= _score(current_q):
            continue
        key = (q["bond_crossings"] + q["overlapping_atoms"], q["stretched_bonds"], energy)
        if best is None or key < best[0]:
            best = (key, candidate, q)
    if best is None:
        return current_q
    q = best[2]
    mol.RemoveAllConformers()
    mol.AddConformer(Chem.Conformer(best[1].GetConformer()), assignId=True)
    return q


def _hub_layout(mol):
    """A candidate for star-shaped molecules (Pd(PPh3)4, Wilkinson's catalyst, tetraphenylmethane): each branch around
    the most central branching atom is laid out on its own, then the branches are set evenly around that atom, each
    pointing outward. None when the molecule has no such hub."""
    hub, best = None, 0
    for atom in mol.GetAtoms():
        if atom.GetDegree() < 3 or atom.IsInRing():
            continue
        cut = Chem.RWMol(mol)
        cut.RemoveAtom(atom.GetIdx())
        smallest = min(len(f) for f in Chem.GetMolFrags(cut))
        if smallest > best:
            hub, best = atom.GetIdx(), smallest
    if hub is None or best < 2:
        return None
    neighbours = [n.GetIdx() for n in mol.GetAtomWithIdx(hub).GetNeighbors()]
    cut = Chem.RWMol(mol)
    cut.RemoveAtom(hub)
    branches = [[i if i < hub else i + 1 for i in frag] for frag in Chem.GetMolFrags(cut)]
    pos = np.zeros((mol.GetNumAtoms(), 2))
    order = []
    for branch in branches:
        anchors = [n for n in neighbours if n in branch]
        order += [(neighbours.index(a), a, branch) for a in anchors[:1]]
    order.sort()
    for k, (_, anchor, branch) in enumerate(order):
        # The branch with the hub atom attached, so its own layout leaves room for the bond to the hub.
        keep = set(branch) | {hub}
        piece = Chem.RWMol(mol)
        for i in sorted(set(range(mol.GetNumAtoms())) - keep, reverse=True):
            piece.RemoveAtom(i)
        piece = piece.GetMol()
        members = sorted(keep)
        try:
            piece.UpdatePropertyCache(strict=False)
            Chem.GetSymmSSSR(piece)
            with rdBase.BlockLogs():
                rdDepictor.Compute2DCoords(piece)
        except Exception:  # noqa: BLE001
            return None
        local = dict(zip(members, piece.GetConformer().GetPositions()[:, :2]))
        angle = 2 * np.pi * k / len(order) + np.pi / 2
        direction = np.array([np.cos(angle), np.sin(angle)])
        bond = local[anchor] - local[hub]
        scale = 1.5 / (np.linalg.norm(bond) or 1.0)
        turn = np.arctan2(direction[1], direction[0]) - np.arctan2(bond[1], bond[0])
        c, s_ = np.cos(turn), np.sin(turn)
        for i in branch:
            pos[i] = ((local[i] - local[hub]) * scale) @ np.array([[c, s_], [-s_, c]])
    conf = Chem.Conformer(mol.GetNumAtoms())
    for i, (x, y) in enumerate(pos):
        conf.SetAtomPosition(i, Point3D(float(x), float(y), 0.0))
    out = Chem.Mol(mol)
    out.RemoveAllConformers()
    out.AddConformer(conf, assignId=True)
    return out


def stereo_faithful(mol) -> bool:
    """Do this layout's wedges (as RDKit would draw them) read back as the input's stereochemistry? Written as a
    molfile with these coordinates and read again; compared on canonical isomeric SMILES plus axial twists."""
    from .atropisomer import stereo_key
    has_stereo = any(a.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED for a in mol.GetAtoms()) or \
        any(b.GetStereo() != Chem.BondStereo.STEREONONE for b in mol.GetBonds())
    if not has_stereo:
        return True
    try:
        with rdBase.BlockLogs():
            back = Chem.MolFromMolBlock(Chem.MolToMolBlock(mol), removeHs=False)
        return back is not None and stereo_key(back) == stereo_key(mol)
    except Exception:  # noqa: BLE001 - cannot check: do not claim it is faithful
        return False


def _is_helicene(mol) -> bool:
    from .stereounits import helices
    return bool(helices(mol))


def _has_macrocycle(mol, size=12) -> bool:
    """A ring of at least `size` atoms (porphyrin inner ring 16, corrole 15, crown ethers, macrolides)."""
    Chem.GetSymmSSSR(mol)
    return any(len(r) >= size for r in mol.GetRingInfo().AtomRings())


def _is_cage(mol) -> bool:
    """Polyhedral cages have atoms shared by three or more rings (cubane, adamantane, C60)."""
    heavy = [a for a in mol.GetAtoms() if a.GetAtomicNum() != 1]
    Chem.GetSymmSSSR(mol)  # ring information may have been cleared by editing the molecule
    ring_info = mol.GetRingInfo()
    return any(ring_info.NumAtomRings(a.GetIdx()) >= 3 for a in heavy)


def _chelated_metal(mol):
    """(metal index, four donor N indices) for a metal bonded to the four macrocycle nitrogens of a drawn porphyrin,
    phthalocyanine or corrin; it may carry further (axial) ligands, such as tin's two chlorides or B12's cyanide."""
    from .coordination import METALS, _macrocycle_nitrogens
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in METALS:
            continue
        ring_n = set(_macrocycle_nitrogens(mol))
        donors = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in ring_n]
        if len(donors) == 4:
            return atom.GetIdx(), donors
    return None


def _axial_branches(mol, metal, donors):
    """The metal's other ligands, each as the list of atoms reachable from it without passing the metal."""
    branches = []
    for start in (n.GetIdx() for n in mol.GetAtomWithIdx(metal).GetNeighbors() if n.GetIdx() not in donors):
        seen, todo = {start}, [start]
        while todo:
            for n in mol.GetAtomWithIdx(todo.pop()).GetNeighbors():
                if n.GetIdx() not in seen and n.GetIdx() != metal:
                    seen.add(n.GetIdx()); todo.append(n.GetIdx())
        if not seen & set(donors):  # a ligand tethered back into the ring system is not axial
            branches.append(sorted(seen, key=lambda i: (i != start, i)))
    return branches


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
    """(query, coordinates): the porphine skeleton (four N-rings joined by four one-atom bridges) as a pattern that
    ignores bond orders, aromaticity and hydrogens, with its textbook square layout from RDKit's ring templates.
    A bridge may be carbon (porphyrins, chlorins, hydroporphyrins such as F430) or nitrogen (phthalocyanines,
    porphyrazines); extra fused rings (F430, chlorophyll, the benzo rings of a phthalocyanine) defeat RDKit's own
    template, but not this pattern."""
    porphine = Chem.MolFromSmiles("c1cc2cc3ccc(cc4ccc(cc5ccc(cc1n2)[nH]5)n4)[nH]3")
    rdDepictor.Compute2DCoords(porphine, useRingTemplates=True)
    query = Chem.RWMol(Chem.MolFromSmarts(re.sub(r"[:=\-]", "~", Chem.MolToSmarts(porphine)).replace("&H1", "")))
    match = porphine.GetSubstructMatch(query)
    ring_info = porphine.GetRingInfo()
    for k, i in enumerate(match):
        if porphine.GetAtomWithIdx(i).GetAtomicNum() == 6 and not ring_info.IsAtomInRingOfSize(i, 5):
            bridge = rdqueries.AtomNumEqualsQueryAtom(6)  # a meso bridge: carbon or nitrogen
            bridge.ExpandQuery(rdqueries.AtomNumEqualsQueryAtom(7), Chem.CompositeQueryType.COMPOSITE_OR)
            query.ReplaceAtom(k, bridge)
    query = query.GetMol()
    conf = porphine.GetConformer()
    return query, np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in match])


def _shape_rmsd(a, b) -> float:
    """RMSD after the best similarity transform (translation, rotation or reflection, uniform scale) of a onto b."""
    a, b = a - a.mean(axis=0), b - b.mean(axis=0)
    u, sv, vt = np.linalg.svd(a.T @ b)
    scale = sv.sum() / max(1e-9, (a ** 2).sum())
    return float(np.sqrt(((scale * a @ u @ vt - b) ** 2).sum(axis=1).mean()))


def _pin_porphyrinoid(ligand, force=False) -> bool:
    """If the ligand has a porphyrin-type core that RDKit laid out badly (its fused rings folded inward), lay it out
    again with the core pinned to the textbook porphine square; RDKit then places fused rings and side chains."""
    query, template = _porphyrinoid_core()
    match = ligand.GetSubstructMatch(query)
    if not match:
        return False
    conf = ligand.GetConformer()
    current = np.array([[conf.GetAtomPosition(i).x, conf.GetAtomPosition(i).y] for i in match])
    if _shape_rmsd(current, template) < 0.25 and not force:  # already the textbook square (heme, chlorophyll)
        return False
    coords = {a: template[k] for k, a in enumerate(match)}
    # A ring fused to the core (F430's lactam along one bond; chlorophyll's ring E along ring C and a bridge carbon)
    # is not placed well by RDKit once the core is fixed. Put its other atoms on an arc of equal bonds between the ends
    # of the shared path, on the far side from that path (for a single shared bond, from the core ring it borders).
    rings = [list(r) for r in Chem.GetSymmSSSR(ligand)]
    core = set(match)
    bond = float(np.median([np.linalg.norm(coords[b.GetBeginAtomIdx()] - coords[b.GetEndAtomIdx()])
                            for b in ligand.GetBonds() if b.GetBeginAtomIdx() in core and b.GetEndAtomIdx() in core]))
    # Repeat outward, ring by ring: a ring fused to an already placed ring (the outer ring of naphthalocyanine's
    # naphtho groups) is placed the same way, instead of being left to RDKit and bent.
    progress = True
    while progress:
      progress = False
      placed = set(coords)
      for ring in rings:
          shared = [a for a in ring if a in placed]
          free = [a for a in ring if a not in placed]
          if len(shared) < 2 or not free or len(ring) > 8 or any(a in coords for a in free):
              continue
          k = next(i for i, a in enumerate(ring) if a in placed and ring[(i + 1) % len(ring)] not in placed)
          ring = ring[k:] + ring[:k]  # ring[0]: last shared atom before the free run; then the free atoms in order
          run = [a for a in ring[1:] if a not in placed]
          if ring[1 + len(run):] != [a for a in ring if a in placed][1:]:
              continue  # shared atoms are not one contiguous path
          p0, p1 = coords[ring[0]], coords[ring[len(run) + 1]]
          if len(shared) == 2:
              host = next((r for r in rings if set(shared) <= set(r) and set(r) <= placed), None)
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
          progress = True
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
        if result["method"] == "projection":
            _project_cage(mol)  # hydrogens from the same 3D view, pointing outward
            return {**result, **quality(mol)}
        conf = skeleton.GetConformer()
        coord_map = {heavy_ids[k]: Point2D(conf.GetAtomPosition(k).x, conf.GetAtomPosition(k).y) for k in range(len(heavy_ids))}
        with rdBase.BlockLogs():
            rdDepictor.Compute2DCoords(mol, coordMap=coord_map)
        placed = quality(mol)
        # RDKit's own all-at-once layout is sometimes tidier for small molecules (glucose's ring hydrogens); keep it
        # only when it has strictly fewer crossings, overlaps and stretched bonds (and the skeleton needed no fix-up).
        if result["method"] == "rdkit" and not _chelated_metal(skeleton) and not _has_macrocycle(skeleton):
            # (Never for a macrocycle: RDKit's all-at-once layout of chlorin e6 is round again.)
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
        # Axial ligands (tin's chlorides, B12's cyanide) point out of the paper; draw them along the diagonals between
        # the M–N bonds, one bond length out, keeping each ligand's own shape.
        n_dirs = [np.array(list(conf.GetAtomPosition(n)))[:2] - centre[:2] for n in donors]
        base = np.arctan2(n_dirs[0][1], n_dirs[0][0]) + np.pi / 4
        for k, branch in enumerate(_axial_branches(ligand, metal, donors)):
            angle = base + np.pi / 2 * [0, 2, 1, 3][k % 4]
            direction = np.array([np.cos(angle), np.sin(angle)])
            pts = np.array([list(conf.GetAtomPosition(i))[:2] for i in branch])
            rel = pts - pts[0]
            spread = rel[1:].mean(axis=0) if len(rel) > 1 else direction
            turn = np.arctan2(direction[1], direction[0]) - np.arctan2(spread[1], spread[0])
            c, s_ = np.cos(turn), np.sin(turn)
            rel = rel @ np.array([[c, s_], [-s_, c]])
            first = centre[:2] + 1.5 * direction
            for i, p in zip(branch, rel):
                conf.SetAtomPosition(i, Point3D(*(first + p), 0.0))
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        return {"method": "rdkit", **quality(mol)}
    rdDepictor.Compute2DCoords(mol)
    best_q = quality(mol)
    if _has_macrocycle(mol):
        # Porphyrins, corroles, corrins, phthalocyanines and other macrocycles: RDKit's ring templates give the
        # textbook shapes (a square porphyrin, not a large circle); a porphyrin-type core that the templates miss
        # (fused rings, nitrogen bridges, unusual tautomers) is pinned to the porphine square. Kept unless worse.
        trial = Chem.Mol(mol)
        rdDepictor.Compute2DCoords(trial, useRingTemplates=True)
        _pin_porphyrinoid(trial)
        candidates = [trial]
        if _score(quality(trial))[:2] > (0, 0):
            # The template's core is right but side groups collide (chlorin e6's neighbouring acids): pin anyway,
            # which points every side group outward.
            forced = Chem.Mol(trial)
            if _pin_porphyrinoid(forced, force=True):
                candidates.append(forced)
        # The textbook layout wins unless it adds bond crossings or more than two crowded atom pairs: a right core
        # with two side-group oxygens close together (chlorin e6) beats a round, distorted macrocycle.
        default_q = best_q
        chosen = None
        for candidate in candidates:
            q = quality(candidate)
            fine = (q["bond_crossings"] <= default_q["bond_crossings"] and q["stretched_bonds"] <= default_q["stretched_bonds"]
                    and q["overlapping_atoms"] <= default_q["overlapping_atoms"] + 2)
            if fine and (chosen is None or _score(q) < _score(chosen[0])):
                chosen = (q, candidate)
        if chosen:
            mol.RemoveAllConformers()
            mol.AddConformer(Chem.Conformer(chosen[1].GetConformer()), assignId=True)
            best_q = chosen[0]
    if (best_q["bond_crossings"] or best_q["overlapping_atoms"] or best_q["stretched_bonds"]
            or best_q["max_bond_length_ratio"] > MAX_RATIO):
        best_q = _best_candidate(mol, best_q)
    result = {"method": "rdkit", **best_q}
    if (best_q["bond_crossings"] or best_q["overlapping_atoms"]) and _is_helicene(mol):
        # A flat helicene ([6] and up) overlaps its end rings. A view of its 3D shape (as for cages) is the textbook
        # picture: the end rings side by side, slightly foreshortened. A projection has no depth, so it implies
        # neither P nor M; the label gives the handedness.
        trial = Chem.Mol(mol)
        try:
            _project_cage(trial)
            q = quality(trial)
            if q["bond_crossings"] + q["overlapping_atoms"] < best_q["bond_crossings"] + best_q["overlapping_atoms"]:
                mol.RemoveAllConformers()
                mol.AddConformer(Chem.Conformer(trial.GetConformer()), assignId=True)
                return {"method": "projection", "helicene": True, **q}
        except Exception:  # noqa: BLE001 - no 3D model: keep the flat layout
            pass
    if not (best_q["bond_crossings"] or best_q["stretched_bonds"] or best_q["overlapping_atoms"]):
        return result
    # A cage cannot be drawn flat without crossings: draw it as a view of its 3D shape instead. A merely crowded
    # ordinary layout keeps RDKit's.
    if not _is_cage(mol):
        return result
    trial = Chem.Mol(mol)
    try:
        _project_cage(trial)
    except Exception:  # noqa: BLE001 - no 3D model: keep RDKit's layout
        return result
    best_conf = Chem.Conformer(trial.GetConformer())
    result = {"method": "projection", **quality(trial)}
    mol.RemoveAllConformers()
    mol.AddConformer(best_conf, assignId=True)
    return result
