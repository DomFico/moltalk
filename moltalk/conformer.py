"""3D data for rotating the flat drawing.

Coordinates come from RDKit ETKDGv3 followed by MMFF94 (or UFF) cleanup. ETKDG cannot
embed fullerene-like cages, so those start from a spectral (graph Laplacian) sphere.
The conformer is rotated (never reflected) to best match the flat drawing, so a drag
can start from the drawing itself and lift it into 3D.
"""
import json
import time
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdBase, rdCIPLabeler, rdMolTransforms
from rdkit.Geometry import Point3D
from .depiction import _cage, layout
from .stereo import stereogenic_unspecified, organic_stereo, assign_cip
from .coordination import is_metal
from . import atropisomer, complexes, stereounits
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


def _store_energies(mol, results):
    for conf, (_, energy) in zip(mol.GetConformers(), results):
        conf.SetDoubleProp("energy", float(energy))


ENERGY_WINDOW = 8.0  # kcal/mol above the best candidate; a strained, half-collapsed shape can fit a flat drawing better


def _spiro_centres(mol):
    """(centre, ring-A neighbours, ring-B neighbours) for spiro atoms joining two small rings (3- or 4-membered)."""
    Chem.GetSymmSSSR(mol)
    rings = [set(r) for r in mol.GetRingInfo().AtomRings()]
    found = []
    for atom in mol.GetAtoms():
        i = atom.GetIdx()
        mine = [r for r in rings if i in r]
        if len(mine) != 2 or mine[0] & mine[1] != {i} or max(len(r) for r in mine) > 4:
            continue
        a = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in mine[0]]
        b = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in mine[1]]
        if len(a) == 2 and len(b) == 2:
            found.append((i, a, b))
    return found


def _square_spiro(mol, conf_id) -> bool:
    """A spiro carbon joining two small rings has them perpendicular, on a straight axis (spiropentane,
    spiropentadiene: D2d). MMFF and UFF get spiropentadiene wrong (MMFF even rates a twisted 56° shape lower in
    energy, and pulls back against restraints), so set it geometrically: rotate ring B, with everything attached to
    it, rigidly about the spiro centre so its axis points straight away from ring A and its plane is perpendicular
    to ring A's. Bond lengths and each ring's own shape are unchanged."""
    from rdkit.Geometry import Point3D
    centres = _spiro_centres(mol)
    if not centres:
        return False
    conf = mol.GetConformer(conf_id)
    pos = conf.GetPositions()
    unit = lambda v: v / np.linalg.norm(v)
    for c, a, b in centres:
        va, vb = pos[a] - pos[c], pos[b] - pos[c]
        axis_a = unit(va.sum(axis=0))
        normal_a = unit(np.cross(va[0], va[1]))
        axis_b, across_b = unit(vb.sum(axis=0)), vb[0] - vb[1]
        across_b = unit(across_b - (across_b @ axis_b) * axis_b)
        current = np.column_stack([axis_b, across_b, np.cross(axis_b, across_b)])
        # Keep ring B's atoms on the side of ring A's plane they are on: always turning them onto one fixed side
        # rotated ring B by up to 180 degrees, which inverts a chiral spiro system (2,6-dichlorospiro[3.3]heptane).
        side_sign = 1.0 if across_b @ normal_a >= 0 else -1.0
        target = np.column_stack([-axis_a, side_sign * normal_a, np.cross(-axis_a, side_sign * normal_a)])
        rot = target @ current.T
        side, todo = set(b), list(b)  # ring B and everything attached to it, not crossing the centre
        while todo:
            for n in mol.GetAtomWithIdx(todo.pop()).GetNeighbors():
                if n.GetIdx() != c and n.GetIdx() not in side:
                    side.add(n.GetIdx()); todo.append(n.GetIdx())
        if side & set(a):
            continue  # the rings are also joined elsewhere: not a simple spiro centre
        for i in side:
            conf.SetAtomPosition(i, Point3D(*(rot @ (pos[i] - pos[c]) + pos[c])))
        pos = conf.GetPositions()
    return True


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


def _allene_twists(mol) -> list[tuple[int, int, int, int]]:
    """(ligand, end, end, ligand) for each allene C=C=C, to keep its two ends perpendicular."""
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Exception:  # noqa: BLE001
        return []
    found = []
    for atom in kekule.GetAtoms():
        doubles = [b for b in atom.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
        if len(doubles) != 2 or atom.GetDegree() != 2:
            continue
        ends = [b.GetOtherAtomIdx(atom.GetIdx()) for b in doubles]
        ligands = [[n.GetIdx() for n in kekule.GetAtomWithIdx(e).GetNeighbors() if n.GetIdx() != atom.GetIdx()]
                   for e in ends]
        if all(ligands) and not any(len([b for b in kekule.GetAtomWithIdx(e).GetBonds()
                                         if b.GetBondType() == Chem.BondType.DOUBLE]) > 1 for e in ends):
            found.append((ligands[0][0], ends[0], ends[1], ligands[1][0]))
    return found


def _vsepr_angles(mol) -> tuple[list[tuple[int, int, int, float]], list[str]]:
    """Bond angles to hold at centres a force field cannot describe: MMFF94 has no parameters for radicals,
    carbenes or localised carbanions, and RDKit's hybridisation labels for them are unreliable (it calls the methyl
    radical and both carbenes sp3), so the shape comes from electron counting (VSEPR):
      * localised carbanion, three bonds and a lone pair: pyramidal, 109.5 degrees (MMFF built CH3- and the ammonium
        and sulfonium ylide carbons flat). Not when conjugated (benzyl, enolate, cyclopentadienyl: planar is right),
        and not next to P+ or As+: phosphonium ylide carbons are near planar.
      * sigma radical, two bonds and one unpaired electron (vinyl, phenyl): bent at 135 degrees (MMFF: linear).
      * carbene, two bonds and two nonbonding electrons: SMILES does not say singlet or triplet. With a heteroatom
        or halogen neighbour the singlet is the ground state (about 105 degrees: CCl2, NHCs); otherwise the triplet
        (about 136 degrees: CH2). MMFF built both linear.
    Angles inside rings of up to 6 atoms are left to the ring. Three-bond radicals (methyl, tBu) are left alone: MMFF
    builds them planar to slightly pyramidal, as they are. Returns (angles, notes for the method string)."""
    angles, notes = [], set()
    ring_info = mol.GetRingInfo()

    def in_small_ring(i, j, k):
        return any(len(r) <= 6 and {i, j, k} <= set(r) for r in ring_info.AtomRings())

    for atom in mol.GetAtoms():
        c = atom.GetIdx()
        neighbours = [n.GetIdx() for n in atom.GetNeighbors()]
        radicals, charge, z = atom.GetNumRadicalElectrons(), atom.GetFormalCharge(), atom.GetAtomicNum()
        target = None
        if z == 6 and charge == -1 and not radicals and len(neighbours) == 3 and not atom.GetIsAromatic():
            conjugated = any(b.GetBondType() != Chem.BondType.SINGLE or b.GetIsConjugated()
                             or any(nb.GetBondType() in (Chem.BondType.DOUBLE, Chem.BondType.TRIPLE, Chem.BondType.AROMATIC)
                                    for nb in b.GetOtherAtom(atom).GetBonds())
                             for b in atom.GetBonds())
            ylide_p = any(n.GetSymbol() in ("P", "As") and n.GetFormalCharge() > 0 for n in atom.GetNeighbors())
            if not conjugated and not ylide_p:
                target, note = 109.5, "carbanion pyramidal"
        elif z == 6 and not charge and len(neighbours) == 2 and radicals == 1:
            target, note = 135.0, "sigma radical bent"
        elif z == 6 and not charge and len(neighbours) == 2 and radicals == 2:
            hetero = any(n.GetAtomicNum() not in (1, 6) for n in atom.GetNeighbors())
            target, note = (105.0, "carbene drawn as singlet (heteroatom neighbour)") if hetero else \
                (136.0, "carbene drawn as triplet (no heteroatom neighbour)")
        if target is None:
            continue
        for x in range(len(neighbours)):
            for y in range(x + 1, len(neighbours)):
                i, k = neighbours[x], neighbours[y]
                if not in_small_ring(i, c, k):
                    angles.append((i, c, k, target))
                    notes.add(note)
    return angles, sorted(notes)


def _optimize(mol, method):
    """Force-field clean-up of every conformer on mol. Large molecules get a shorter clean-up: ETKDG geometry is
    already reasonable, and the viewer needs a picture, not an energy minimum (chlorophyll: 1.7 s of MMFF)."""
    iterations = 300 if _large_and_flexible(mol) else 1000
    twists = _allene_twists(mol)
    angles, notes = _vsepr_angles(mol)
    if (twists or angles) and AllChem.MMFFHasAllMoleculeParams(mol):
        # MMFF94 and UFF both flatten an allene (their torsion terms across the sp carbon are zero, so 1,4 contacts
        # win): ETKDG's 84 degrees became 180 for 1,3-dichloroallene, which also destroys its chirality. Hold each
        # allene's ends perpendicular, keeping the sense ETKDG chose. MMFF has no parameters for radicals, carbenes
        # or localised carbanions either (see _vsepr_angles): hold those centres at their VSEPR angles.
        props = AllChem.MMFFGetMoleculeProperties(mol)
        results = []
        for conf in mol.GetConformers():
            ff = AllChem.MMFFGetMoleculeForceField(mol, props, confId=conf.GetId())
            for quad in twists:
                current = rdMolTransforms.GetDihedralDeg(conf, *quad)
                target = 90.0 if current >= 0 else -90.0
                ff.MMFFAddTorsionConstraint(*quad, False, target, target, 100.0)
            for i, j, k, target in angles:
                ff.MMFFAddAngleConstraint(i, j, k, False, target - 2.0, target + 2.0, 100.0)
            ff.Minimize(maxIts=iterations)
            results.append((0, ff.CalcEnergy()))
        _store_energies(mol, results)
        held = (["allene ends held perpendicular"] if twists else []) + notes
        return method + f" + MMFF94 ({'; '.join(held)})"
    if AllChem.MMFFHasAllMoleculeParams(mol):
        results = AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=iterations)
        _store_energies(mol, results)
        return method + " + MMFF94"
    if AllChem.UFFHasAllMoleculeParams(mol):
        results = AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=iterations)
        _store_energies(mol, results)
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


def _embed(heavy, conformers=8, enforce_chirality=True):
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
        # ETKDG's "basic knowledge" keeps aromatic rings flat, which a helicene's twisted rings cannot satisfy: it ran
        # into the 10 s timeout on [6]helicene, while ETDG embeds it in milliseconds.
        from .stereounits import helices
        params = AllChem.ETDG() if helices(heavy) else AllChem.ETKDGv3()
        params.randomSeed = 0xC0FFEE
        params.enforceChirality = enforce_chirality
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
        return mol, _optimize(mol, "ETDG" if isinstance(params, type(AllChem.ETDG())) and helices(heavy) else "ETKDGv3")


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
        spread = float(np.sum(w * aligned[:, :2] ** 2))
        # One heavy atom (water, methane, a methyl radical) gives the fit nothing to size by: it returned scale 0 and
        # collapsed every shown hydrogen onto the atom. Keep the 3D model's own size then (RDKit's 2D bond length
        # of 1.5 matches a C–C bond in angstroms).
        scale = fixed if fixed is not None else \
            (float(np.sum(w * aligned[:, :2] * q[:, :2]) / spread) if spread > 1e-6 else 1.0)
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
        if atom.GetAtomicNum() == 1 and atom.GetDegree() == 1 and i >= embedded.GetNumHeavyAtoms():
            parent = atom.GetNeighbors()[0].GetIdx()
            rank = used.get(parent, 0)
            used[parent] = rank + 1
            out[i] = pos[hydrogens[parent][rank]]
        elif atom.HasProp("_eta_atoms"):
            continue  # a pi ligand's ring or alkene centre: placed below, once its atoms are
        else:
            out[i] = pos[i]
    for atom in drawn.GetAtoms():
        if atom.HasProp("_eta_atoms"):
            out[atom.GetIdx()] = out[json.loads(atom.GetProp("_eta_atoms"))].mean(axis=0)
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
    elif len(b) == 3 and domains == 1 and hyb in (Chem.HybridizationType.SP3, Chem.HybridizationType.SP2):
        # Three bonds and one lone pair (amine, carbanion, ylide carbon, pyrrole N): read the shape actually built,
        # not RDKit's hybridisation label. Pyramidal: opposite the bonds. Flat: in the p orbital, along the normal.
        # (Opposite the bonds of a flat atom is a near-zero vector: the Wittig ylide's CH2- pair pointed 14 degrees
        # from a C-H bond.)
        away = -sum(b)
        normal = _unit(np.cross(b[1] - b[0], b[2] - b[0]))
        if np.linalg.norm(away) > 0.3 or normal is None:
            out = [_unit(away)]
        else:
            out = [normal if np.dot(normal, away) >= 0 else -normal]
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
    if len(out) < domains or not _clear_of_bonds(out[:domains], b):
        # Anything else (hypervalent XeF4, ClF3...), or a template that does not fit the geometry actually built:
        # spread the domains by repulsion from the real bonds, so a lone pair never sits on a bond.
        out = _repel(b, domains)
    return out[:domains]


LONE_PAIR_MIN_ANGLE = 80.0  # degrees from any bond (a p-orbital pair is at 90; VSEPR pairs at 105 or more)


def _clear_of_bonds(dirs, bonds) -> bool:
    limit = np.cos(np.radians(LONE_PAIR_MIN_ANGLE))
    if any(np.dot(d, v) > limit for d in dirs for v in bonds):
        return False
    return all(np.dot(dirs[i], dirs[j]) < np.cos(np.radians(60)) for i in range(len(dirs)) for j in range(i + 1, len(dirs)))


def _same_centres(mol, conf, wanted: dict) -> bool:
    """Do the stereocentres of this model have the input's R/S labels? (Only those atoms: a metal with four
    neighbours would otherwise count as a new stereocentre.)"""
    probe = Chem.Mol(mol, confId=conf.GetId())
    Chem.AssignStereochemistryFrom3D(probe)
    assign_cip(probe, atoms=sorted(wanted), bonds=[])
    return all(probe.GetAtomWithIdx(i).HasProp("_CIPCode") and probe.GetAtomWithIdx(i).GetProp("_CIPCode") == cip
               for i, cip in wanted.items())


def _same_stereo(mol, conf, reference: str) -> bool:
    """Does this model's geometry give the input's isomeric SMILES?"""
    probe = Chem.Mol(mol, confId=conf.GetId())
    Chem.AssignStereochemistryFrom3D(probe)
    return Chem.MolToSmiles(Chem.RemoveHs(probe)) == reference


def _has_configuration(mol, conf, wanted) -> bool:
    for unit, sense in wanted:
        measured = stereounits.measure(mol, conf, unit)
        if measured is None or measured["sense"] != sense:
            return False
    return True


def conformer_3d(mol_in, drawn_heavy, drawn) -> dict:
    """Data for rotating the flat drawing in 3D. mol_in is the parsed input (atom order preserved). drawn_heavy is
    the molecule as drawn without explicit hydrogens (metal–N bonds added for porphyrin-type complexes); drawn is
    what the flat drawing shows (drawn_heavy, or it with every hydrogen). The conformer is chosen and aligned from
    heavy atoms only, so showing or hiding hydrogens never changes the view."""
    heavy_count = mol_in.GetNumAtoms()
    shown = drawn.GetNumAtoms()
    # The drawing without hydrogens: the input's atoms plus, for a pi complex, its ring/alkene centres.
    core_count = drawn_heavy.GetNumAtoms()
    flat_heavy = Chem.Mol(drawn_heavy)
    layout(flat_heavy)  # the same deterministic layout draw_molecule uses
    xy_core = flat_heavy.GetConformer().GetPositions()[:, :2]
    xy_heavy = xy_core[:heavy_count]
    if shown > core_count:
        flat = Chem.Mol(drawn)
        layout(flat)
        xy = flat.GetConformer().GetPositions()[:, :2]
    else:
        xy = xy_core
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
    # Allenes, helicenes, Xabab spiro atoms whose configuration the input states (from its 3D coordinates).
    stated_units = stereounits.stated(mol_in)
    wanted_units = [(u, stated_units[u["key"]]) for u in stereounits.units(mol_in) if u["key"] in stated_units] \
        if not site else []
    mirror_ok = not stereounits._other_chirality(mol_in)
    # A complex assembled from loose pieces (Wilkinson's catalyst): the pieces are embedded as stored, then set
    # around the metal in its ideal geometry (complexes.arrange).
    complex_info = complexes.info(drawn_heavy) if not site else None
    chelate_reference = None
    if complex_info is not None and complex_info.get("chelates"):
        bonded = complexes.embedding_source(drawn_heavy, heavy_count)
        if bonded is not None:
            # ETKDG could not embed Jacobsen's Mn–salen with its chirality enforced (the metal closes rings through the
            # stereocentres' neighbours); without enforcement it embeds in a fraction of a second, and conformers
            # whose stereocentres do not match the input are dropped below.
            source = bonded
            chelate_reference = {a.GetIdx(): a.GetProp("_CIPCode") for a in mol_in.GetAtoms()
                                 if a.HasProp("_CIPCode") and not is_metal(a)}
    # An Xaabb spiro system (read with the newer perception): ETKDG does not enforce its tags, so check each model's
    # stereo against the input's; if every chiral tag belongs to the spiro system, the mirror image fixes a wrong one.
    spiro_units = stereounits.xaabb(mol_in) if not site else []
    spiro_reference = Chem.MolToSmiles(mol_in) if any(u["specified"] for u in spiro_units) else None
    spiro_atoms = set().union(*(set(u["atoms"]) for u in spiro_units)) if spiro_units else set()
    spiro_mirror_ok = spiro_reference is not None and all(
        a.GetIdx() in spiro_atoms or a.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED for a in mol_in.GetAtoms())
    # Large molecules (chlorophyll) get one round of a few candidates: each embedding costs ~0.5 s on a desktop and
    # ~1.5 s on Cloud Run, and the refinements below (robust alignment, tail swing, cavity fit) do the shaping.
    rounds = (LARGE_CANDIDATES,) if _large_and_flexible(mol_in) else CANDIDATE_CONFORMERS
    for count in rounds:
        if best is not None and time.monotonic() - started > RETRY_BUDGET_S:
            break  # a big, slow molecule (erythromycin): more conformers would mostly time out
        try:
            candidate_mol, candidate_method = _embed(source, count, enforce_chirality=chelate_reference is None)
        except ValueError:
            if best is not None:  # the larger retry failed (e.g. a big macrolide timing out): keep the first result
                break
            raise
        if site:
            candidate_method += ", ring embedded as the free base"
        if complex_info is not None:
            for conf in candidate_mol.GetConformers():
                complexes.arrange(candidate_mol, conf.GetId(), complex_info)
            candidate_method += f", ligands set {complex_info['geometry'].replace('_', ' ')} around the metal"
        energies = [c.GetDoubleProp("energy") for c in candidate_mol.GetConformers() if c.HasProp("energy")]
        # Not for porphyrin-type macrocycles: MMFF rates a bowed naphthalocyanine 37 kcal/mol *below* the flat one, and
        # the flattest-ring preference already picks the right shape there.
        lowest = min(energies) if energies and not macrocycle else None
        for conf in candidate_mol.GetConformers():
            if lowest is not None and conf.HasProp("energy") and conf.GetDoubleProp("energy") > lowest + ENERGY_WINDOW:
                continue  # strained (spiropentane: one candidate half-collapsed to 45° fitted the bowtie drawing best)
            if chelate_reference and not _same_centres(candidate_mol, conf, chelate_reference):
                continue
            if spiro_reference is not None and not _same_stereo(candidate_mol, conf, spiro_reference):
                if not spiro_mirror_ok:
                    continue
                stereounits._mirror(conf)  # all the chirality is in the spiro system: its mirror image is the one asked
                if not _same_stereo(candidate_mol, conf, spiro_reference):
                    continue
            if wanted_units and not _has_configuration(candidate_mol, conf, wanted_units):
                if not mirror_ok:
                    continue  # ETKDG does not know allene or helix configurations: keep only matching conformers
                stereounits._mirror(conf)  # nothing else is chiral: the mirror image has the stated configuration
                if not _has_configuration(candidate_mol, conf, wanted_units):
                    continue
            candidate = evaluate(candidate_mol, conf.GetId())
            if best is None or candidate[0] < best[0]:
                best = candidate
                mol, method = candidate_mol, candidate_method
        # Stop when the wedges agree and the shape matches the drawing. A poor fit means the random conformers missed
        # the drawing's shape (8 conformers of 2-bromobutane had no anti chain, so lifting swung the methyl ~2 bond
        # lengths, which looks like a jump); more conformers fix that, within the time budget above.
        # Only small molecules retry for shape: there one misplaced atom is a large part of the picture, while big,
        # flexible ones (heme's side chains) never fit a flat drawing closely and would just take longer.
        if best is not None and (candidate_mol.GetNumConformers() == 1 or (
                best[0][0] == 0 and (best[0][2] <= GOOD_FIT or heavy_count > FIT_RETRY_MAX_ATOMS))):
            break
    # Long open chains (chlorophyll's phytyl, fatty acids) come out of ETKDG crumpled, often folded back over the
    # rest of the molecule. Straighten them to the extended zigzag, the textbook and low-energy shape that the flat
    # drawing also shows, and keep it only if it matches the drawing better.
    # Not for an assembled complex: these steps turn pieces about single bonds, and its loose ligands look like tails.
    straightened = None if complex_info is not None else _optional(_extend_chains, mol, best[1])
    if straightened is not None:
        candidate = evaluate(straightened, best[1])
        if candidate[0][0] <= best[0][0] and candidate[0][2] < best[0][2]:
            best, mol = candidate, straightened
    unfolded = None if complex_info is not None else \
        _optional(_swing_tail, mol, best[1], lambda m: evaluate(m, best[1]), heavy_count)
    if unfolded is not None:
        candidate = evaluate(unfolded, best[1])
        if candidate[0][0] <= best[0][0] and candidate[0][2] < best[0][2]:
            best, mol = candidate, unfolded
    if _optional(_square_spiro, mol, best[1]):
        best = evaluate(mol, best[1])
    # Last, so no later unrestrained minimisation undoes it; only the chosen conformer, as it costs a minimisation.
    if site and _optional(_square_cavity, mol, site, drawn_heavy, best[1]):
        best = evaluate(mol, best[1])
        method += ", N4 cavity fitted to the metal"
    # ETKDG embeds a specified biaryl twist, but the tail swing and chain straightening above rotate about single
    # bonds; make sure the axis still has the verified helicity (turn one side back over if not).
    axial = [entry for entry in atropisomer.describe(mol_in) if entry["verified"]]
    if axial and not site and _optional(atropisomer.enforce, mol, best[1], axial):
        best = evaluate(mol, best[1])
    if best is None:
        raise ValueError("No 3D model with the stated allene/helix/spiro configuration could be built; "
                         "use the flat drawing.")
    (contradicted, _, rmsd), conf_id, xyz, raw, heavy_rot, heavy_scale, vsepr_lone = best
    frame_rot = heavy_rot
    to_heavy_frame = np.eye(3)
    if shown > core_count:
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
    # Allene and helicene labels on their label bonds: the stated configuration, or the one this model happens to have.
    unit_labels, open_units = {}, []
    for unit in stereounits.units(mol_in):
        if unit["type"] == "spiro":
            if unit["key"] not in stated_units:
                open_units.append(unit)
            continue
        measured = stereounits.measure(mol, mol.GetConformer(conf_id), unit)
        a, b = stereounits._label_bond(unit)
        if unit["key"] in stated_units:
            unit_labels[tuple(sorted((a, b)))] = (stated_units[unit["key"]], False)
        else:
            open_units.append(unit)
            if measured and measured.get("descriptor"):
                unit_labels[tuple(sorted((a, b)))] = (measured["descriptor"], True)
    # A hindered axis the input leaves open: the model has one twist; name it (both rules must agree) as arbitrary.
    open_axes = []
    for bond_index in atropisomer.candidate_axes(mol_in):
        a, b = mol_in.GetBondWithIdx(bond_index).GetBeginAtomIdx(), mol_in.GetBondWithIdx(bond_index).GetEndAtomIdx()
        measured = atropisomer.from_geometry(mol, mol.GetConformer(conf_id), mol.GetBondBetweenAtoms(a, b))
        open_axes.append({"atom_indices": sorted([a, b]), "cip": measured[1] if measured and
                          {"P": "Sa", "M": "Ra"}[measured[0]] == measured[1] else None})
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
        axis = next((x for x in axial if sorted(x["atom_indices"]) == pair), None)
        open_axis = next((x for x in open_axes if x["atom_indices"] == pair), None)
        unit_tag = unit_labels.get(tuple(pair))
        if unit_tag is not None:  # allene M/P, helicene P/M; "arbitrary" when the input leaves it open
            entry["cip"] = unit_tag[0]
            if unit_tag[1]:
                entry["arbitrary"] = True
        elif axis is not None:
            entry["cip"] = axis["cip"]  # Ra/Sa beside the axis bond, verified against this kind of 3D geometry
        elif open_axis is not None and open_axis["cip"]:
            entry["cip"], entry["arbitrary"] = open_axis["cip"], True
        elif source_bond is not None and source_bond.HasProp("_CIPCode") and source_bond.GetProp("_CIPCode") in ("E", "Z"):
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
        if shown <= core_count:  # hidden hydrogens: their 3D positions come from the embedded molecule
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
                 "h": 0 if i >= heavy_count or shown > core_count else atom.GetTotalNumHs(),
                 "xy": [round(float(c), 4) for c in xy[i]], "xyz": [round(float(c), 4) for c in xyz[i]]}
        if atom.HasProp("_eta_atoms"):
            entry["centroid"] = True  # a pi ligand's ring or alkene centre: no label, only its bond to the metal
        if i < heavy_count and i in stereocentres and mol.GetAtomWithIdx(i).HasProp("_CIPCode") and not is_metal(atom):
            entry["cip"] = mol.GetAtomWithIdx(i).GetProp("_CIPCode")
            if i in unspecified:
                entry["arbitrary"] = True
        atoms.append(entry)
    note = ("Rotation uses one calculated conformer (" + method + "), not a measured or unique structure; "
            "flexible parts can adopt other shapes.")
    result = {"atoms": atoms, "bonds": bonds, "rings": rings, "method": method, "heavy_atom_count": heavy_count,
              "hydrogens_shown": shown > core_count, "fit_rmsd": round(rmsd, 3), "lone_pair_dirs": lone_pairs,
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
    if unspecified or unspecified_bonds or open_axes or open_units:
        parts = ([f"atom(s) {', '.join(map(str, unspecified))}"] if unspecified else []) + \
                [f"double bond {a}–{b}" for a, b in unspecified_bonds] + \
                [f"axis {x['atom_indices'][0]}–{x['atom_indices'][1]} (atropisomer twist)" for x in open_axes] + \
                [f"the {u['type']} unit at atoms {', '.join(map(str, u['atoms'][:4]))}" for u in open_units]
        if open_units:
            result["arbitrary_stereo_units"] = [{"type": u["type"], "atoms": u["atoms"]} for u in open_units]
        result["arbitrary_stereo_atoms"] = unspecified
        result["arbitrary_stereo_bonds"] = unspecified_bonds
        if open_axes:
            result["arbitrary_stereo_axes"] = [x["atom_indices"] for x in open_axes]
        result["warning"] = (f"Stereochemistry at {'; '.join(parts)} is unspecified in the input; "
                             "this model shows one arbitrary configuration there.")
    return result
