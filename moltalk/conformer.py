"""3D data for rotating the flat drawing.

Coordinates come from RDKit ETKDGv3 followed by MMFF94 (or UFF) cleanup. ETKDG cannot
embed fullerene-like cages, so those start from a spectral (graph Laplacian) sphere.
The conformer is rotated (never reflected) to best match the flat drawing, so a drag
can start from the drawing itself and lift it into 3D.
"""
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdBase, rdCIPLabeler
from rdkit.Geometry import Point3D
from .depiction import _cage, layout
from .stereo import stereogenic_unspecified, organic_stereo, assign_cip
from .coordination import is_metal
from .depiction import _chelated_metal

EMBED_TIMEOUT_S = 5
CANDIDATE_CONFORMERS = (8, 32)  # try more only if no conformer agrees with the flat wedges
SINGLE_TIMEOUT_S = 10  # last resort: one conformer, when several could not be made in time (big chelates on a slow CPU)
FIT_RETRY_MAX_ATOMS = 15
LARGE_ATOMS = 50
LARGE_CANDIDATES = 2
LARGE_ROTATABLE = 5
GOOD_FIT = 0.15  # weighted RMSD (bond-length units) below which a conformer already lifts smoothly from the drawing
RETRY_BUDGET_S = 1.5  # ...and only if the first round was quick


def _cage_like(mol, cage) -> bool:
    if len(cage) < 20:
        return False
    cage_set = set(cage)
    edges = sum(1 for b in mol.GetBonds() if b.GetBeginAtomIdx() in cage_set and b.GetEndAtomIdx() in cage_set)
    return edges / len(cage) >= 1.4  # fullerenes 1.5; flat polycyclic aromatics stay below ~1.3


def _spectral_sphere(mol, cage, bond_length=1.42):
    index = {a: k for k, a in enumerate(cage)}
    lap = np.zeros((len(cage), len(cage)))
    for b in mol.GetBonds():
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in index and j in index:
            i, j = index[i], index[j]
            lap[i, j] = lap[j, i] = -1
            lap[i, i] += 1
            lap[j, j] += 1
    _, vectors = np.linalg.eigh(lap)
    pos = vectors[:, 1:4]
    pos -= pos.mean(axis=0)
    pos /= np.linalg.norm(pos, axis=1, keepdims=True)  # on the unit sphere
    edges = [(index[b.GetBeginAtomIdx()], index[b.GetEndAtomIdx()]) for b in mol.GetBonds()
             if b.GetBeginAtomIdx() in index and b.GetEndAtomIdx() in index]
    mean = np.mean([np.linalg.norm(pos[i] - pos[j]) for i, j in edges])
    return pos * (bond_length / mean)


def _large_and_flexible(mol) -> bool:
    """The speed shortcut (fewer candidates, shorter force-field clean-up) is for big molecules with long flexible
    chains (chlorophyll; erythromycin has 7 rotatable bonds), which are slow. A big rigid aromatic one (naphthalocyanine) is fast anyway and needs the
    full clean-up to come out flat."""
    from rdkit.Chem import rdMolDescriptors
    if mol.GetNumHeavyAtoms() <= LARGE_ATOMS:
        return False
    heavy = Chem.RemoveHs(mol, sanitize=False)
    Chem.FastFindRings(heavy)  # removing atoms leaves no ring information
    return rdMolDescriptors.CalcNumRotatableBonds(heavy) >= LARGE_ROTATABLE


def _optimize(mol, method):
    """Force-field clean-up of every conformer on mol. Large molecules get a shorter clean-up: ETKDG geometry is
    already reasonable, and the viewer needs a picture, not an energy minimum (chlorophyll: 1.7 s of MMFF)."""
    iterations = 300 if _large_and_flexible(mol) else 1000
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=iterations)
        return method + " + MMFF94"
    if AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=iterations)
        return method + " + UFF"
    return method


def _free_base(drawn_heavy, site):
    """Embedding source for a drawn metal chelate: the metal's bonds are removed and each nitrogen that had a covalent
    (not dative) bond to it gets an N–H, i.e. the neutral free-base ring. MMFF treats charged porphyrin forms badly
    (it bends the ring out of plane); the free base comes out nearly flat, as the real ring is. The extra hydrogens
    are not drawn, and the metal is put back at the centre of the four nitrogens."""
    metal, donors = site
    rw = Chem.RWMol(drawn_heavy)
    for d in donors:
        bond = rw.GetBondBetweenAtoms(d, metal)
        covalent = bond.GetBondType() != Chem.BondType.DATIVE
        rw.RemoveBond(d, metal)
        if covalent:
            atom = rw.GetAtomWithIdx(d)
            atom.SetNumExplicitHs(atom.GetNumExplicitHs() + 1)
            atom.SetNoImplicit(True)
    rw.GetAtomWithIdx(metal).SetFormalCharge(0)
    out = rw.GetMol()
    Chem.SanitizeMol(out)
    return out


def _embed(heavy, conformers=8):
    """Return (molecule with explicit Hs and its conformers, method). Heavy-atom indices are unchanged; hydrogen
    order may differ from Chem.AddHs(heavy), so map hydrogens through _drawn_positions."""
    cage = _cage(heavy)
    with rdBase.BlockLogs():
        if _cage_like(heavy, cage) and len(cage) == heavy.GetNumAtoms():
            # ETKDG cannot embed fullerene-like cages: start from a spectral sphere and add Hs geometrically.
            sphere = _spectral_sphere(heavy, cage)
            conf = Chem.Conformer(heavy.GetNumAtoms())
            for k, a in enumerate(cage):
                conf.SetAtomPosition(a, Point3D(*sphere[k]))
            mol = Chem.Mol(heavy)
            mol.RemoveAllConformers()
            mol.AddConformer(conf, assignId=True)
            mol = Chem.AddHs(mol, addCoords=True)
            return mol, _optimize(mol, "spectral cage")
        mol = Chem.AddHs(heavy)
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        params.timeout = EMBED_TIMEOUT_S
        status = -1
        # Several conformers, so the one that best matches the flat drawing can be chosen; then random starting
        # coordinates; then, as a last resort, a single conformer with a longer budget. On timeout RDKit returns [-1]
        # and no conformers, so count the conformers actually made. A first attempt that ran out of time means the
        # second would too (Cloud Run's CPU is several times slower than a desktop), so go straight to the last.
        started = time.monotonic()
        for count, random_coords, timeout in ((conformers, False, EMBED_TIMEOUT_S), (conformers, True, EMBED_TIMEOUT_S),
                                              (1, True, SINGLE_TIMEOUT_S)):
            if count > 1 and random_coords and time.monotonic() - started > 0.5 * EMBED_TIMEOUT_S:
                continue
            params.useRandomCoords = random_coords
            params.timeout = timeout
            try:
                AllChem.EmbedMultipleConfs(mol, count, params)
                status = 0 if mol.GetNumConformers() else -1
            except RuntimeError:
                status = -1
            if status >= 0:
                break
        if status < 0:
            raise ValueError("RDKit could not generate 3D coordinates for this structure; use the flat drawing.")
        return mol, _optimize(mol, "ETKDGv3")


# Typical metal–N distances (Å) in tetrapyrrole complexes; the four N are restrained to the square these imply.
METAL_N = {"Mg": 2.05, "Fe": 2.0, "Ni": 1.95, "Co": 1.95, "Zn": 2.05, "Cu": 2.0, "Mn": 2.0, "Pd": 2.0, "Pt": 2.0}


def _square_cavity(mol, site, drawn_heavy, conf_id) -> bool:
    """Free-base embedding leaves the ring's central cavity lopsided (F430: Ni–N 1.78–2.44 Å once the metal is put at
    the centre). Re-optimise each conformer with the four donor N restrained to the square the metal needs
    (cis N···N = d·√2, trans = 2d), so the metal ends up with four equal bonds."""
    metal, donors = site
    if len(donors) != 4:
        return False
    d = METAL_N.get(drawn_heavy.GetAtomWithIdx(metal).GetSymbol(), 2.0)
    pos = mol.GetConformer(conf_id).GetPositions()
    # Order the donors round the ring (by angle about their centroid) to tell cis from trans pairs.
    centre = pos[donors].mean(axis=0)
    u, _, vt = np.linalg.svd(pos[donors] - centre)
    order = sorted(donors, key=lambda n: np.arctan2(*((pos[n] - centre) @ vt[:2].T)[::-1]))
    pairs = [(order[k], order[(k + 1) % 4], d * np.sqrt(2)) for k in range(4)] + \
            [(order[0], order[2], 2 * d), (order[1], order[3], 2 * d)]
    use_mmff = AllChem.MMFFHasAllMoleculeParams(mol)
    props = AllChem.MMFFGetMoleculeProperties(mol) if use_mmff else None
    with rdBase.BlockLogs():
        ff = (AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf_id) if use_mmff
              else AllChem.UFFGetMoleculeForceField(mol, confId=conf_id))
        if ff is None:
            return False
        for a, b, length in pairs:
            (ff.MMFFAddDistanceConstraint if use_mmff else ff.UFFAddDistanceConstraint)(
                a, b, False, length - 0.02, length + 0.02, 500.0)
        ff.Minimize(maxIts=1000)
    return True


def _extend_chains(mol, conf_id, min_chain=4):
    """A copy of mol whose conformer conf_id has every torsion along open sp3 carbon chains set to anti (180°) and
    is then force-field relaxed; None if there is no such chain of at least min_chain carbons."""
    from rdkit.Chem import rdMolTransforms
    chain = lambda a: a.GetSymbol() == "C" and not a.IsInRing() and a.GetHybridization() == Chem.HybridizationType.SP3
    carbons = {a.GetIdx() for a in mol.GetAtoms() if chain(a)}
    if len(carbons) < min_chain:
        return None
    torsions = []
    for bond in mol.GetBonds():
        b, c = bond.GetBeginAtom(), bond.GetEndAtom()
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing() or b.GetIdx() not in carbons or c.GetIdx() not in carbons:
            continue
        a = next((n.GetIdx() for n in b.GetNeighbors() if n.GetIdx() != c.GetIdx() and n.GetAtomicNum() > 1), None)
        d = next((n.GetIdx() for n in c.GetNeighbors() if n.GetIdx() != b.GetIdx() and n.GetAtomicNum() > 1), None)
        if a is not None and d is not None:
            torsions.append((a, b.GetIdx(), c.GetIdx(), d))
    if len(torsions) < min_chain - 2:
        return None
    out = Chem.Mol(mol)
    conf = out.GetConformer(conf_id)
    for t in torsions:
        try:
            rdMolTransforms.SetDihedralDeg(conf, *t, 180.0)
        except (RuntimeError, ValueError):
            continue
    with rdBase.BlockLogs():
        if AllChem.MMFFHasAllMoleculeParams(out):
            AllChem.MMFFOptimizeMolecule(out, confId=conf_id, maxIters=500)
        elif AllChem.UFFHasAllMoleculeParams(out):
            AllChem.UFFOptimizeMolecule(out, confId=conf_id, maxIters=500)
    return out


def _optional(step, *args):
    """Run a cosmetic refinement; if RDKit fails on an unusual structure, skip it rather than lose the 3D view."""
    try:
        return step(*args)
    except Exception:  # noqa: BLE001 - any failure here only means the unrefined model is used
        return None


def _swing_tail(mol, conf_id, score, heavy_count, max_bonds=4, min_tail=6):
    """A long chain hanging off a ring system (chlorophyll's phytyl ester) can leave the ring pointing back across
    it. Scan the first few torsions of its attachment (60/180/300°, staggered), keep the clash-free combination that
    best matches the flat drawing (where the chain points away), then relax. None if there is no such chain."""
    from itertools import product
    from rdkit.Chem import rdMolTransforms
    ring_atoms = {a.GetIdx() for a in mol.GetAtoms() if a.IsInRing() and a.GetIdx() < heavy_count}
    if not ring_atoms:
        return None
    dist = Chem.GetDistanceMatrix(mol)
    # The heavy atom farthest (in bonds) from every ring atom ends the tail; walk from it back to the rings.
    reach = {i: min(dist[i][r] for r in ring_atoms) for i in range(heavy_count)}
    reach = {i: d for i, d in reach.items() if d < mol.GetNumAtoms()}  # a separate metal ion is not connected
    far = max(reach, key=reach.get)
    if reach[far] < min_tail:
        return None  # no long chain: nothing to swing (and far may be a ring atom itself)
    start = min(ring_atoms, key=lambda r: dist[far][r])
    path = list(Chem.GetShortestPath(mol, start, far))
    if len(path) - 1 < min_tail:
        return None
    torsions = []
    for k in range(len(path) - 1):
        b, c = path[k], path[k + 1]
        bond = mol.GetBondBetweenAtoms(b, c)
        if bond.GetBondType() != Chem.BondType.SINGLE or bond.IsInRing():
            continue
        a = next((n.GetIdx() for n in mol.GetAtomWithIdx(b).GetNeighbors() if n.GetIdx() != c and n.GetAtomicNum() > 1), None)
        d = next((n.GetIdx() for n in mol.GetAtomWithIdx(c).GetNeighbors() if n.GetIdx() != b and n.GetAtomicNum() > 1), None)
        if a is not None and d is not None:
            torsions.append((a, b, c, d))
        if len(torsions) == max_bonds:
            break
    if not torsions:
        return None
    topological = Chem.GetDistanceMatrix(mol)
    pairs = np.argwhere(np.triu(topological[:heavy_count, :heavy_count] >= 4, 1))
    best = None
    for angles in product((60.0, 180.0, 300.0), repeat=len(torsions)):
        trial = Chem.Mol(mol)
        conf = trial.GetConformer(conf_id)
        try:
            for t, angle in zip(torsions, angles):
                rdMolTransforms.SetDihedralDeg(conf, *t, angle)
        except (RuntimeError, ValueError):
            continue
        pos = conf.GetPositions()
        if len(pairs) and np.min(np.linalg.norm(pos[pairs[:, 0]] - pos[pairs[:, 1]], axis=1)) < 2.5:
            continue  # heavy atoms four or more bonds apart closer than 2.5 Å: a clash
        key = score(trial)[0]
        if best is None or (key[0], key[2]) < best[0]:
            best = ((key[0], key[2]), trial)
    if best is None:
        return None
    out = best[1]
    with rdBase.BlockLogs():
        if AllChem.MMFFHasAllMoleculeParams(out):
            AllChem.MMFFOptimizeMolecule(out, confId=conf_id, maxIters=300)
        elif AllChem.UFFHasAllMoleculeParams(out):
            AllChem.UFFOptimizeMolecule(out, confId=conf_id, maxIters=300)
    return out


def _fit(xyz, xy, weights, scale=None):
    """Proper rotation (no reflection, so R/S is preserved), scale and offsets that best lay 3D coordinates onto the
    flat drawing's 2D layout, by weighted least squares. Returns (transform, weighted RMS misfit)."""
    base = np.asarray(weights, dtype=float)
    w = base[:, None]
    fixed = scale
    # Robust (iteratively reweighted) fit: atoms that cannot match the drawing, such as a long flexible tail folded
    # differently in 3D, are down-weighted, so the rigid part decides the rotation. Plain least squares let
    # chlorophyll's phytyl tail (15% of the weight, but residuals ~10x larger) turn the macrocycle away from the
    # drawing: 1.8-3.5 misfit, against 0.14 for the ring alone.
    for _ in range(5):
        p0, q0 = (w * xyz).sum(0) / w.sum(), (w * xy).sum(0) / w.sum()
        p = xyz - p0
        q = np.column_stack([xy - q0, np.zeros(len(xy))])
        u, _, vt = np.linalg.svd((w * p).T @ q)
        d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
        rot = vt.T @ np.diag([1, 1, d]) @ u.T
        aligned = p @ rot.T
        scale = fixed if fixed is not None else \
            float(np.sum(w * aligned[:, :2] * q[:, :2]) / max(np.sum(w * aligned[:, :2] ** 2), 1e-9))
        aligned *= scale
        residual = np.linalg.norm(aligned[:, :2] - q[:, :2], axis=1)
        c = max(0.5, 2.0 * float(np.median(residual)))
        w = (base / (1.0 + (residual / c) ** 2))[:, None]
    rmsd = float(np.sqrt((base * residual ** 2).sum() / base.sum()))
    transform = lambda pts: (pts - p0) @ rot.T * scale + np.append(q0, 0.0)
    transform.scale = scale
    return transform, rmsd, rot


def _drawn_positions(embedded, conf_id, drawn, site):
    """3D positions in the drawn molecule's atom order: heavy atoms by index, each hydrogen matched to the same-rank
    hydrogen on the same parent atom (extra free-base N–H hydrogens are simply not drawn)."""
    pos = embedded.GetConformer(conf_id).GetPositions()
    out = np.zeros((drawn.GetNumAtoms(), 3))
    hydrogens = {}
    for atom in embedded.GetAtoms():
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1:
            hydrogens.setdefault(atom.GetNeighbors()[0].GetIdx(), []).append(atom.GetIdx())
    used = {}
    for atom in drawn.GetAtoms():
        i = atom.GetIdx()
        if atom.GetAtomicNum() == 1 and i >= embedded.GetNumHeavyAtoms():
            parent = atom.GetNeighbors()[0].GetIdx()
            rank = used.get(parent, 0)
            used[parent] = rank + 1
            out[i] = pos[hydrogens[parent][rank]]
        else:
            out[i] = pos[i]
    if site:
        metal, donors = site
        out[metal] = out[donors].mean(axis=0)
        # Axial ligands (tin's chlorides, B12's cyanide) were embedded as a separate piece: put them above and below
        # the ring plane, about 2 Å from the metal, keeping each ligand's own shape.
        from .depiction import _axial_branches
        _, _, vt = np.linalg.svd(out[donors] - out[metal])
        normal = vt[2]
        for k, branch in enumerate(_axial_branches(drawn, metal, donors)):
            side = normal if k % 2 == 0 else -normal
            rel = out[branch] - out[branch[0]]
            spread = rel[1:].mean(axis=0) if len(branch) > 1 else side
            if np.linalg.norm(spread) > 1e-6:
                u = spread / np.linalg.norm(spread)
                axis = np.cross(u, side)
                sin, cos = np.linalg.norm(axis), float(np.dot(u, side))
                if sin > 1e-6:
                    axis /= sin
                    kx = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
                    rot = np.eye(3) + sin * kx + (1 - cos) * kx @ kx
                    rel = rel @ rot.T
                elif cos < 0:
                    rel = -rel
            out[branch] = out[metal] + 2.0 * side + rel
    return out


def _unit(v):
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else None


def _perpendicular(axis, hint=None):
    """A unit vector perpendicular to axis, preferably toward hint."""
    for candidate in ([hint] if hint is not None else []) + [np.array([0.0, 0.0, 1.0]), np.array([0.0, 1.0, 0.0])]:
        w = _unit(candidate - np.dot(candidate, axis) * axis)
        if w is not None:
            return w
    return np.array([1.0, 0.0, 0.0])


def _repel(bonds, count):
    """Directions for count electron domains that stay as far as possible from the bonds and from each other."""
    rng = np.random.default_rng(0)
    dirs = [_unit(v) for v in rng.normal(size=(count, 3))]
    for _ in range(400):
        moved = []
        for k, d in enumerate(dirs):
            force = np.zeros(3)
            for other in list(bonds) + [x for j, x in enumerate(dirs) if j != k]:
                diff = d - other
                force += diff / (np.linalg.norm(diff) ** 3 + 1e-6)
            moved.append(_unit(d + 0.05 * (force - np.dot(force, d) * d)))
        dirs = moved
    return dirs


# VSEPR templates for expanded octets: (ligand slots, lone-pair slots). Lone pairs take equatorial positions in a
# trigonal bipyramid and trans positions in an octahedron.
_AX, _EQ = np.array([0.0, 0.0, 1.0]), [np.array([np.cos(t), np.sin(t), 0.0]) for t in np.radians([0, 120, 240])]
_OCT = [np.array(v, dtype=float) for v in ([1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1])]


def _vsepr_slots(ligands, pairs):
    if ligands + pairs == 5:
        equatorial, axial = _EQ, [_AX, -_AX]
        lone = equatorial[:pairs]
        return axial + equatorial[pairs:], lone  # AX4E seesaw, AX3E2 T-shape, AX2E3 linear, AX5 bipyramid
    if ligands + pairs == 6:
        lone = [_OCT[5], _OCT[4]][:pairs]  # AX5E square pyramid; AX4E2 square planar (pairs trans)
        return [v for v in _OCT if not any(v is l for l in lone)], lone
    return None


def vsepr_centres(drawn, electron_map):
    """Main-group atoms with 5 or 6 electron domains whose ligands are all single atoms (SF4, XeF4, ClF3, I3-, PCl5,
    SF6): ETKDG would make them tetrahedral, so their geometry is set from the textbook VSEPR shape."""
    from .coordination import is_metal
    found = {}
    for atom in drawn.GetAtoms():
        if is_metal(atom) or atom.GetAtomicNum() == 1:
            continue
        info = electron_map.get(str(atom.GetIdx()), {"pairs": 0, "radicals": 0})
        ligands = [n.GetIdx() for n in atom.GetNeighbors()]
        slots = _vsepr_slots(len(ligands), info["pairs"])
        if slots and all(drawn.GetAtomWithIdx(n).GetDegree() == 1 for n in ligands):
            found[atom.GetIdx()] = (ligands, slots)
    return found


def _apply_vsepr(pos, centres):
    """Put each centre's ligands on its VSEPR slots (keeping their mean bond length); returns lone-pair slot vectors."""
    lone = {}
    for c, (ligands, (ligand_slots, lone_slots)) in centres.items():
        length = float(np.mean([np.linalg.norm(pos[n] - pos[c]) for n in ligands]))
        for n, slot in zip(ligands, ligand_slots):
            pos[n] = pos[c] + length * slot
        lone[c] = lone_slots
    return lone


def lone_pair_directions(atom, centre, neighbours, domains):
    """Unit vectors for an atom's lone pairs (and an unpaired electron, last), from VSEPR and its hybridisation.
    neighbours: unit vectors to bonded atoms (hidden hydrogens included); domains: pairs + unpaired electron."""
    b = [n for n in neighbours if n is not None]
    hyb = atom.GetHybridization()
    out = []
    if not b:  # a free ion (Cl-): tetrahedral
        out = [np.array(v) / np.sqrt(3) for v in ([1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1])]
    elif hyb == Chem.HybridizationType.SP3 and len(b) + domains == 4:
        if len(b) == 3:
            out = [_unit(-sum(b))]
        elif len(b) == 2:
            u, v = _unit(-(b[0] + b[1])), _unit(np.cross(b[0], b[1]))
            half = np.radians(109.47 / 2)
            out = [_unit(u * np.cos(half) + v * np.sin(half)), _unit(u * np.cos(half) - v * np.sin(half))] if v is not None else _repel(b, 2)
        else:  # one bond: three positions at 109.5° from it, staggered
            axis = b[0]
            e1 = _perpendicular(axis)
            e2 = np.cross(axis, e1)
            t = np.radians(109.47)
            out = [_unit(axis * np.cos(t) + np.sin(t) * (np.cos(phi) * e1 + np.sin(phi) * e2)) for phi in np.radians([0, 120, 240])]
    elif hyb == Chem.HybridizationType.SP2:
        normal = _unit(np.cross(b[0], b[1])) if len(b) >= 2 else None
        if normal is None:  # one bond (C=O): the plane of the neighbour's own substituents
            normal = neighbours_plane_normal = getattr(atom, "_plane_normal", None)
            if normal is None:
                normal = _perpendicular(b[0])
        in_plane = max(0, 3 - len(b))
        if len(b) == 1:
            w = _unit(np.cross(normal, b[0]))
            t = np.radians(120)
            out = [_unit(b[0] * np.cos(t) + w * np.sin(t)), _unit(b[0] * np.cos(t) - w * np.sin(t))][:in_plane]
        elif len(b) == 2 and in_plane:
            out = [_unit(-(b[0] + b[1]))]
        out += [normal, -normal][:domains - len(out)]  # the rest sit in the p orbital (pyrrole N, carboxylate O-)
    elif hyb == Chem.HybridizationType.SP and len(b) == 1:
        e1 = _perpendicular(b[0])
        out = [-b[0], e1, -e1][:domains]
    out = [d for d in out if d is not None]
    if len(out) < domains:  # anything else (hypervalent XeF4, ClF3...): spread by repulsion
        out = _repel(b, domains)
    return out[:domains]


def conformer_3d(mol_in, drawn_heavy, drawn) -> dict:
    """Data for rotating the flat drawing in 3D. mol_in is the parsed input (atom order preserved). drawn_heavy is
    the molecule as drawn without explicit hydrogens (metal–N bonds added for porphyrin-type complexes); drawn is
    what the flat drawing shows (drawn_heavy, or it with every hydrogen). The conformer is chosen and aligned from
    heavy atoms only, so showing or hiding hydrogens never changes the view."""
    heavy_count = mol_in.GetNumAtoms()
    shown = drawn.GetNumAtoms()
    flat_heavy = Chem.Mol(drawn_heavy)
    layout(flat_heavy)  # the same deterministic layout draw_molecule uses
    xy_heavy = flat_heavy.GetConformer().GetPositions()[:, :2]
    if shown > heavy_count:
        flat = Chem.Mol(drawn)
        layout(flat)
        xy = flat.GetConformer().GetPositions()[:, :2]
    else:
        xy = xy_heavy
    site = _chelated_metal(drawn_heavy)
    source = _free_base(drawn_heavy, site) if site else mol_in

    potential = organic_stereo(mol_in, Chem.FindPotentialStereo(mol_in))
    # Stereocentres and their neighbours dominate the fit: their local shape is rigid, so matching it makes the
    # flat drawing's wedges point toward the viewer in 3D too. Flexible arms swing into place as the drawing lifts.
    weights = np.ones(heavy_count)
    for st in potential:
        if str(st.type) == "Atom_Tetrahedral":
            weights[st.centeredOn] = 10
            for nb in flat_heavy.GetAtomWithIdx(st.centeredOn).GetNeighbors():
                weights[nb.GetIdx()] = 10
    if site:
        # A metal chelate's macrocycle is the rigid frame the eye follows; anchor the alignment on it (and the metal)
        # so the ring lifts in place. Otherwise five stereocentres and chlorophyll's 20-atom phytyl tail pulled the
        # ring ~3 bond lengths away from the drawing before any rotation.
        from .depiction import _porphyrinoid_core
        core = mol_in.GetSubstructMatch(_porphyrinoid_core()[0]) or tuple(_cage(mol_in))
        weights = np.where(weights > 1, 3.0, 1.0)
        for i in core:
            if i < heavy_count:
                weights[i] = 20
        weights[site[0]] = 20
    # The flat (heavy-atom) drawing's wedges: (stereocentre, substituent, +1 toward the viewer / -1 away).
    Chem.WedgeMolBonds(flat_heavy, flat_heavy.GetConformer())
    wedges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), 1 if b.GetBondDir() == Chem.BondDir.BEGINWEDGE else -1)
              for b in flat_heavy.GetBonds() if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)]
    # Porphyrin-type rings are essentially flat; among the candidates prefer the flattest ring.
    from .depiction import _porphyrinoid_core
    macrocycle = list(mol_in.GetSubstructMatch(_porphyrinoid_core()[0])) or (_cage(mol_in) if site else [])
    from .chemistry import electrons
    electron_map = electrons(drawn)
    vsepr = vsepr_centres(drawn_heavy, electron_map)
    def evaluate(candidate_mol, conf_id):
        pos = _drawn_positions(candidate_mol, conf_id, drawn, site)
        vsepr_lone = _apply_vsepr(pos, vsepr)
        transform, fit, rot = _fit(pos[:heavy_count], xy_heavy, weights)
        aligned = transform(pos)
        contradicted = sum(np.sign(aligned[end, 2] - aligned[start, 2]) != sign for start, end, sign in wedges)
        bend = 0.0
        if macrocycle:
            ring = pos[macrocycle] - pos[macrocycle].mean(axis=0)
            bend = round(float(np.linalg.svd(ring, compute_uv=False)[2] / np.sqrt(len(macrocycle))), 2)
        return (contradicted, bend, fit), conf_id, aligned, pos, rot, transform.scale, vsepr_lone

    best = mol = method = None
    started = time.monotonic()
    # Large molecules (chlorophyll) get one round of a few candidates: each embedding costs ~0.5 s on a desktop and
    # ~1.5 s on Cloud Run, and the refinements below (robust alignment, tail swing, cavity fit) do the shaping.
    rounds = (LARGE_CANDIDATES,) if _large_and_flexible(mol_in) else CANDIDATE_CONFORMERS
    for count in rounds:
        if best is not None and time.monotonic() - started > RETRY_BUDGET_S:
            break  # a big, slow molecule (erythromycin): more conformers would mostly time out
        try:
            candidate_mol, candidate_method = _embed(source, count)
        except ValueError:
            if best is not None:  # the larger retry failed (e.g. a big macrolide timing out): keep the first result
                break
            raise
        if site:
            candidate_method += ", ring embedded as the free base"
        for conf in candidate_mol.GetConformers():
            candidate = evaluate(candidate_mol, conf.GetId())
            if best is None or candidate[0] < best[0]:
                best = candidate
                mol, method = candidate_mol, candidate_method
        # Stop when the wedges agree and the shape matches the drawing. A poor fit means the random conformers missed
        # the drawing's shape (8 conformers of 2-bromobutane had no anti chain, so lifting swung the methyl ~2 bond
        # lengths, which looks like a jump); more conformers fix that, within the time budget above.
        # Only small molecules retry for shape: there one misplaced atom is a large part of the picture, while big,
        # flexible ones (heme's side chains) never fit a flat drawing closely and would just take longer.
        if candidate_mol.GetNumConformers() == 1 or (best[0][0] == 0 and (best[0][2] <= GOOD_FIT or heavy_count > FIT_RETRY_MAX_ATOMS)):
            break
    # Long open chains (chlorophyll's phytyl, fatty acids) come out of ETKDG crumpled, often folded back over the
    # rest of the molecule. Straighten them to the extended zigzag, the textbook and low-energy shape that the flat
    # drawing also shows, and keep it only if it matches the drawing better.
    straightened = _optional(_extend_chains, mol, best[1])
    if straightened is not None:
        candidate = evaluate(straightened, best[1])
        if candidate[0][0] <= best[0][0] and candidate[0][2] < best[0][2]:
            best, mol = candidate, straightened
    unfolded = _optional(_swing_tail, mol, best[1], lambda m: evaluate(m, best[1]), heavy_count)
    if unfolded is not None:
        candidate = evaluate(unfolded, best[1])
        if candidate[0][0] <= best[0][0] and candidate[0][2] < best[0][2]:
            best, mol = candidate, unfolded
    # Last, so no later unrestrained minimisation undoes it; only the chosen conformer, as it costs a minimisation.
    if site and _optional(_square_cavity, mol, site, drawn_heavy, best[1]):
        best = evaluate(mol, best[1])
        method += ", N4 cavity fitted to the metal"
    (contradicted, _, rmsd), conf_id, xyz, raw, heavy_rot, heavy_scale, vsepr_lone = best
    frame_rot = heavy_rot
    to_heavy_frame = np.eye(3)
    if shown > heavy_count:
        # The explicit-H drawing is laid out in its own orientation (often turned over relative to the heavy-atom
        # drawing), so align the same conformer to it: the lift then starts from exactly the drawing on screen.
        # to_heavy_frame lets the viewer keep the same view when hydrogens are toggled.
        transform, _, h_rot = _fit(raw[:heavy_count], xy[:heavy_count], weights, scale=heavy_scale)  # same size
        xyz = transform(raw)
        to_heavy_frame = heavy_rot @ h_rot.T
        frame_rot = h_rot
    Chem.AssignStereochemistryFrom3D(mol, confId=conf_id)

    real = stereogenic_unspecified(mol_in, potential)
    unspecified = sorted(i for kind, i in real if kind == "atom")
    # R/S is shown only at real stereocentres: specified in the input, or unspecified but genuinely stereogenic
    # (in a cubane the cage carbons are neither, although a 3D model gives every one of them a CIP label).
    stereocentres = {st.centeredOn for st in potential if str(st.type) == "Atom_Tetrahedral" and str(st.specified) == "Specified"}
    stereocentres |= set(unspecified)
    unspecified_bonds = [[mol_in.GetBondWithIdx(i).GetBeginAtomIdx(), mol_in.GetBondWithIdx(i).GetEndAtomIdx()]
                         for kind, i in real if kind == "bond"]
    # Label only what is shown: the unbounded labeller never finished on dodecahedrane's 3D model (20 tagged centres
    # in a highly symmetric cage, none of them real stereocentres).
    open_pairs = [mol.GetBondBetweenAtoms(a, b) for a, b in unspecified_bonds]
    assign_cip(mol, atoms=sorted(i for i in stereocentres if i < heavy_count),
               bonds=[bd.GetIdx() for bd in open_pairs if bd is not None])

    open_bonds = [sorted(pair) for pair in unspecified_bonds]
    kekule = Chem.Mol(drawn)
    Chem.Kekulize(kekule, clearAromaticFlags=True)
    rings = [list(r) for r in mol_in.GetRingInfo().AtomRings()]  # without metal bonds: chelate rings are not drawn rings
    bond_rings = mol_in.GetRingInfo().BondRings()
    bonds = []
    for b in kekule.GetBonds():
        containing = [k for k, r in enumerate(bond_rings) if b.GetIdx() in r] if b.GetIdx() < mol_in.GetNumBonds() else []
        entry = {"a": b.GetBeginAtomIdx(), "b": b.GetEndAtomIdx(), "order": 1 if b.GetBondType() == Chem.BondType.DATIVE
                 else int(b.GetBondTypeAsDouble()), "ring": min(containing, key=lambda k: len(rings[k])) if containing else None}
        if b.GetBondType() == Chem.BondType.DATIVE:
            entry["dative"] = True  # drawn as an arrow toward the metal
        # E/Z: from the input where it specifies the double bond; at an unspecified one, the configuration this
        # conformer happens to have (marked arbitrary), like R/S at an unspecified centre.
        pair = sorted((entry["a"], entry["b"]))
        source_bond = mol_in.GetBondBetweenAtoms(*pair) if max(pair) < heavy_count else None
        if source_bond is not None and source_bond.HasProp("_CIPCode") and source_bond.GetProp("_CIPCode") in ("E", "Z"):
            entry["cip"] = source_bond.GetProp("_CIPCode")
        elif pair in open_bonds:
            model_bond = mol.GetBondBetweenAtoms(*pair)
            if model_bond is not None and model_bond.HasProp("_CIPCode") and model_bond.GetProp("_CIPCode") in ("E", "Z"):
                entry["cip"], entry["arbitrary"] = model_bond.GetProp("_CIPCode"), True
        bonds.append(entry)
    # Lone-pair directions from the real 3D geometry (VSEPR by hybridisation), rotated into the output frame.
    embedded = mol.GetConformer(conf_id).GetPositions()
    lone_pairs = {}
    for key, info in electron_map.items():
        i = int(key)
        if i in vsepr_lone:  # expanded octet: the template's own lone-pair slots
            lone_pairs[key] = [[round(float(c), 4) for c in (d @ frame_rot.T)] for d in vsepr_lone[i]]
            continue
        atom = drawn.GetAtomWithIdx(i)
        vectors = [_unit(raw[n.GetIdx()] - raw[i]) for n in atom.GetNeighbors()]
        if shown <= heavy_count:  # hidden hydrogens: their 3D positions come from the embedded molecule
            hs = [h.GetIdx() for h in mol.GetAtomWithIdx(i).GetNeighbors() if h.GetAtomicNum() == 1][:atom.GetTotalNumHs()]
            vectors += [_unit(embedded[h] - raw[i]) for h in hs]
        if len(vectors) == 1 and atom.GetHybridization() == Chem.HybridizationType.SP2:
            other = atom.GetNeighbors()[0]  # C=O: the lone pairs lie in the plane of the carbon's substituents
            around = [_unit(raw[n.GetIdx()] - raw[other.GetIdx()]) for n in other.GetNeighbors() if n.GetIdx() != i]
            if around:
                atom._plane_normal = _unit(np.cross(vectors[0], around[0]))
        domains = info["pairs"] + (1 if info["radicals"] else 0)
        dirs = lone_pair_directions(atom, raw[i], vectors, domains)
        lone_pairs[key] = [[round(float(c), 4) for c in (d @ frame_rot.T)] for d in dirs]

    atoms = []
    for i in range(shown):
        atom = drawn.GetAtomWithIdx(i)
        entry = {"index": i, "element": atom.GetSymbol(), "charge": atom.GetFormalCharge(),
                 "h": 0 if i >= heavy_count or shown > heavy_count else atom.GetTotalNumHs(),
                 "xy": [round(float(c), 4) for c in xy[i]], "xyz": [round(float(c), 4) for c in xyz[i]]}
        if i < heavy_count and i in stereocentres and mol.GetAtomWithIdx(i).HasProp("_CIPCode") and not is_metal(atom):
            entry["cip"] = mol.GetAtomWithIdx(i).GetProp("_CIPCode")
            if i in unspecified:
                entry["arbitrary"] = True
        atoms.append(entry)
    note = ("Rotation uses one calculated conformer (" + method + "), not a measured or unique structure; "
            "flexible parts can adopt other shapes.")
    result = {"atoms": atoms, "bonds": bonds, "rings": rings, "method": method, "heavy_atom_count": heavy_count,
              "hydrogens_shown": shown > heavy_count, "fit_rmsd": round(rmsd, 3), "lone_pair_dirs": lone_pairs,
              "to_heavy_frame": [[round(float(v), 6) for v in row] for row in to_heavy_frame],
              "wedges_contradicted": int(contradicted), "note": note}
    from .coordination import is_metal as _metal
    metals = [a for a in drawn.GetAtoms() if _metal(a) and a.GetDegree() >= 2]
    if metals:
        note += (f" The arrangement of ligands around {metals[0].GetSymbol()} comes from a generic force field, which does not "
                 "know ligand-field preferences (square planar vs tetrahedral), and cis/trans placement is not specified "
                 "by the SMILES; do not read coordination geometry from this model.")
        result["note"] = note
    if vsepr:
        result["note"] = result["note"] + " Expanded-octet centres use their VSEPR shape (" + ", ".join(
            drawn.GetAtomWithIdx(c).GetSymbol() for c in vsepr) + ")."
    if unspecified or unspecified_bonds:
        parts = ([f"atom(s) {', '.join(map(str, unspecified))}"] if unspecified else []) + \
                [f"double bond {a}–{b}" for a, b in unspecified_bonds]
        result["arbitrary_stereo_atoms"] = unspecified
        result["arbitrary_stereo_bonds"] = unspecified_bonds
        result["warning"] = (f"Stereochemistry at {'; '.join(parts)} is unspecified in the input; "
                             "this model shows one arbitrary configuration there.")
    return result
