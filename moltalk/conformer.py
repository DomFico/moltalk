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
from .stereo import stereogenic_unspecified, organic_stereo
from .coordination import is_metal
from .depiction import _chelated_metal

EMBED_TIMEOUT_S = 5
CANDIDATE_CONFORMERS = (8, 32)  # try more only if no conformer agrees with the flat wedges
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


def _optimize(mol, method):
    """Force-field clean-up of every conformer on mol."""
    if AllChem.MMFFHasAllMoleculeParams(mol):
        AllChem.MMFFOptimizeMoleculeConfs(mol, maxIters=1000)
        return method + " + MMFF94"
    if AllChem.UFFHasAllMoleculeParams(mol):
        AllChem.UFFOptimizeMoleculeConfs(mol, maxIters=1000)
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
        for random_coords in (False, True):
            params.useRandomCoords = random_coords
            try:
                # Several conformers, so the one that best matches the flat drawing can be chosen. On timeout RDKit
                # returns [-1] and no conformers, so count the conformers actually made.
                AllChem.EmbedMultipleConfs(mol, conformers, params)
                status = 0 if mol.GetNumConformers() else -1
            except RuntimeError:
                status = -1
            if status >= 0:
                break
        if status < 0:
            raise ValueError("RDKit could not generate 3D coordinates for this structure; use the flat drawing.")
        return mol, _optimize(mol, "ETKDGv3")


def _fit(xyz, xy, weights, scale=None):
    """Proper rotation (no reflection, so R/S is preserved), scale and offsets that best lay 3D coordinates onto the
    flat drawing's 2D layout, by weighted least squares. Returns (transform, weighted RMS misfit)."""
    w = np.asarray(weights, dtype=float)[:, None]
    p0, q0 = (w * xyz).sum(0) / w.sum(), (w * xy).sum(0) / w.sum()
    p = xyz - p0
    q = np.column_stack([xy - q0, np.zeros(len(xy))])
    u, _, vt = np.linalg.svd((w * p).T @ q)
    d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    rot = vt.T @ np.diag([1, 1, d]) @ u.T
    aligned = p @ rot.T
    if scale is None:  # fit the scale too, unless a fixed one is given
        scale = float(np.sum(w * aligned[:, :2] * q[:, :2]) / max(np.sum(w * aligned[:, :2] ** 2), 1e-9))
    aligned *= scale
    rmsd = float(np.sqrt((w[:, 0] * np.sum((aligned[:, :2] - q[:, :2]) ** 2, axis=1)).sum() / w.sum()))
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
    # The flat (heavy-atom) drawing's wedges: (stereocentre, substituent, +1 toward the viewer / -1 away).
    Chem.WedgeMolBonds(flat_heavy, flat_heavy.GetConformer())
    wedges = [(b.GetBeginAtomIdx(), b.GetEndAtomIdx(), 1 if b.GetBondDir() == Chem.BondDir.BEGINWEDGE else -1)
              for b in flat_heavy.GetBonds() if b.GetBondDir() in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH)]
    # Porphyrin-type rings are essentially flat; among the candidates prefer the flattest ring.
    macrocycle = _cage(mol_in) if site else []
    from .chemistry import electrons
    electron_map = electrons(drawn)
    vsepr = vsepr_centres(drawn_heavy, electron_map)
    best = mol = method = None
    started = time.monotonic()
    for count in CANDIDATE_CONFORMERS:
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
            pos = _drawn_positions(candidate_mol, conf.GetId(), drawn, site)
            vsepr_lone = _apply_vsepr(pos, vsepr)
            transform, fit, rot = _fit(pos[:heavy_count], xy_heavy, weights)
            aligned = transform(pos)
            contradicted = sum(np.sign(aligned[end, 2] - aligned[start, 2]) != sign for start, end, sign in wedges)
            bend = 0.0
            if macrocycle:
                ring = pos[macrocycle] - pos[macrocycle].mean(axis=0)
                bend = round(float(np.linalg.svd(ring, compute_uv=False)[2] / np.sqrt(len(macrocycle))), 2)
            key = (contradicted, bend, fit)
            if best is None or key < best[0]:
                best = (key, conf.GetId(), aligned, pos, rot, transform.scale, vsepr_lone)
                mol, method = candidate_mol, candidate_method
        if best[0][0] == 0 or candidate_mol.GetNumConformers() == 1:
            break
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
    rdCIPLabeler.AssignCIPLabels(mol)

    real = stereogenic_unspecified(mol_in, potential)
    unspecified = sorted(i for kind, i in real if kind == "atom")
    # R/S is shown only at real stereocentres: specified in the input, or unspecified but genuinely stereogenic
    # (in a cubane the cage carbons are neither, although a 3D model gives every one of them a CIP label).
    stereocentres = {st.centeredOn for st in potential if str(st.type) == "Atom_Tetrahedral" and str(st.specified) == "Specified"}
    stereocentres |= set(unspecified)
    unspecified_bonds = [[mol_in.GetBondWithIdx(i).GetBeginAtomIdx(), mol_in.GetBondWithIdx(i).GetEndAtomIdx()]
                         for kind, i in real if kind == "bond"]

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
