"""Stereogenic units beyond what SMILES carries: allenes, spiro atoms, helicenes. One representation for every kind.

Every unit MolTalk reports (here and in chemistry.stereo_units) is a dict with:
  type                     tetrahedral | double_bond | axis (biaryl or C–N, see atropisomer.py) | allene | spiro |
                           helix
  atoms, bonds             participating atom / bond indices (input numbering)
  specified                does the input state this unit's configuration?
  configuration            the configuration the input states (descriptor form), else None
  possible_configurations  e.g. ["P", "M"]
  configuration_source     where the configuration came from (SMILES tags, CXSMILES wedge, 3D coordinates)
  verification             "verified" (with how), "unverified" or "not applicable"
  descriptor               the descriptor shown to the user, only when verified
  stability                configurational stability, kept apart from stereogenicity

Three questions are kept separate: is the unit stereogenic (graph analysis below), does the input specify it, and
has its descriptor been checked against the IUPAC definition.

Representation. RDKit's SMILES cannot state the configuration of an allene, a helicene or an 'Xabab' spiro atom
(spiro[4.4]nonane-1,6-dione). MolTalk does not invent a SMILES extension for them: their configuration travels as
standard CXSMILES with 3D coordinates (or a 3D MOL/SDF file), from which it is measured. Plain SMILES, 2D MOL and
CDXML cannot carry it, and exports say so.

'Xaabb' spiro compounds (2,6-dichlorospiro[3.3]heptane) are different: IUPAC names them with R/S at the ring atoms
and the spiro atom (P-93.5.3.5), and RDKit's newer stereo perception represents and labels exactly that. MolTalk
keeps the legacy perception for everything else (switching globally changed adamantane, chlorophyll and more), so
such molecules alone are handled with the newer perception (needs_new_perception, used by chemistry.parse).

Descriptors are computed from 3D geometry with the IUPAC definitions:
  * axial (allenes), P-92.1.2.1.2: near pair a > b, far pair c > d; Ra if a -> b -> c is clockwise looking toward d.
    This is the tetrahedral rule applied to the elongated tetrahedron; its sign convention is calibrated against
    RDKit's R/S (tests: L-alanine is S).
  * helicity, P-92.1.2.2.1: P if the torsion from the near top-ranked ligand to the far one is clockwise (positive).
    For helicenes, the handedness of the helix through the ring centres (right-handed = P).
  For allenes the two rules must agree as IUPAC states for chirality axes (P-92.1.2.2.3, M = Ra; its example
  (1M)-1,3-dichloropropa-1,2-diene = (1Ra)); a unit whose rules disagree gets no descriptor.
  'Xabab' spiro atoms need CIP's special digraph treatment (P-93.5.3.2), which neither RDKit nor MolTalk implements:
  their configuration is represented, enumerated, drawn and built in 3D, but no descriptor is given.
"""
import json
from contextlib import contextmanager

import numpy as np
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem
from rdkit.Geometry import Point3D

from .atropisomer import RANK_DEPTH, _branch, _kekule

HELIX_MIN_RINGS = 5
UNITS_PROP = "_moltalk_units3d"
_H_KEY = tuple([(1,)] + [()] * (RANK_DEPTH - 1))


@contextmanager
def new_perception(enabled: bool = True):
    """RDKit's newer stereo perception for the duration of the block (worker processes are single-threaded)."""
    before = Chem.GetUseLegacyStereoPerception()
    Chem.SetUseLegacyStereoPerception(not enabled)
    try:
        yield
    finally:
        Chem.SetUseLegacyStereoPerception(before)


# ---- geometry -------------------------------------------------------------------------------------------------------

def handedness(a, b, c, d) -> str:
    """'R' if a -> b -> c is clockwise seen with d pointing away from the viewer (the tetrahedral CIP rule), else 'S'.
    Sign calibrated in the tests against RDKit: for (S)-alanine (N, COOH, CH3, H) the determinant is positive."""
    return "R" if np.linalg.det(np.array([a - d, b - d, c - d])) < 0 else "S"


def torsion(p0, p1, p2, p3) -> float:
    b0, b1, b2 = p0 - p1, p2 - p1, p3 - p2
    b1 = b1 / (np.linalg.norm(b1) or 1.0)
    v, w = b0 - np.dot(b0, b1) * b1, b2 - np.dot(b2, b1) * b1
    return float(np.degrees(np.arctan2(np.dot(np.cross(b1, v), w), np.dot(v, w))))


def _position(mol, conf, atom: int, ligand):
    """3D position of a ligand: an atom, or an implicit hydrogen placed opposite the atom's other bonds."""
    if ligand is not None:
        return np.array(conf.GetAtomPosition(ligand))
    explicit = [n.GetIdx() for n in mol.GetAtomWithIdx(atom).GetNeighbors() if n.GetAtomicNum() == 1]
    if explicit:  # a model with explicit hydrogens
        return np.array(conf.GetAtomPosition(explicit[0]))
    centre = np.array(conf.GetAtomPosition(atom))
    others = [np.array(conf.GetAtomPosition(n.GetIdx())) - centre for n in mol.GetAtomWithIdx(atom).GetNeighbors()]
    direction = -sum(v / (np.linalg.norm(v) or 1.0) for v in others)
    return centre + direction / (np.linalg.norm(direction) or 1.0)


def _ligands(k, atom: int, exclude) -> list:
    """(key, atom index or None for an implicit H), highest-ranked first."""
    found = [(_branch(k, n.GetIdx(), atom), n.GetIdx()) for n in k.GetAtomWithIdx(atom).GetNeighbors()
             if n.GetIdx() not in exclude]
    found += [(_H_KEY, None)] * k.GetAtomWithIdx(atom).GetTotalNumHs()
    return sorted(found, key=lambda item: item[0], reverse=True)


# ---- detection ------------------------------------------------------------------------------------------------------

def allenes(mol) -> list[dict]:
    """Cumulenes with an even number of cumulated double bonds (allenes, pentatetraenes) whose two ends each carry two
    different ligands: chirality axes. (An odd number makes a planar, E/Z system, left to RDKit.)"""
    k = _kekule(mol)
    found, seen = [], set()
    for atom in k.GetAtoms():
        doubles = [b for b in atom.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE]
        if len(doubles) != 2 or atom.GetDegree() != 2 or atom.GetIdx() in seen:
            continue
        chain = [atom.GetIdx()]
        ends = []
        for bond in doubles:  # walk outward to each end of the cumulated chain
            prev, cur = atom.GetIdx(), bond.GetOtherAtomIdx(atom.GetIdx())
            while True:
                nxt = [b for b in k.GetAtomWithIdx(cur).GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE
                       and b.GetOtherAtomIdx(cur) != prev]
                if k.GetAtomWithIdx(cur).GetDegree() == 2 and nxt:
                    chain.append(cur)
                    prev, cur = cur, nxt[0].GetOtherAtomIdx(cur)
                else:
                    ends.append((cur, prev))
                    break
        seen |= set(chain)
        if len(chain) % 2 == 0:  # n inner atoms -> n + 1 double bonds; odd count of double bonds is planar
            continue
        (t0, p0), (t1, p1) = ends
        near, far = _ligands(k, t0, {p0}), _ligands(k, t1, {p1})
        if len(near) != 2 or len(far) != 2 or near[0][0] == near[1][0] or far[0][0] == far[1][0]:
            continue
        order = sorted([t0, t1])
        found.append({"type": "allene", "key": f"allene:{order[0]}-{order[1]}", "ends": (t0, t1),
                      "atoms": sorted(set(chain) | {t0, t1}),
                      "bonds": sorted(b.GetIdx() for b in k.GetBonds() if b.GetBeginAtomIdx() in set(chain) | {t0, t1}
                                      and b.GetEndAtomIdx() in set(chain) | {t0, t1}),
                      "near": (p0, near), "far": (p1, far), "possible_configurations": ["M", "P"],
                      "stability": "configurationally stable: allenes do not racemise at room temperature"})
    return found


def _spiro_atoms(mol):
    info = mol.GetRingInfo()
    for atom in mol.GetAtoms():
        if atom.GetDegree() != 4:
            continue
        rings = [set(r) for r in info.AtomRings() if atom.GetIdx() in r]
        if len(rings) != 2 or rings[0] & rings[1] != {atom.GetIdx()}:
            continue
        a = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in rings[0]]
        b = [n.GetIdx() for n in atom.GetNeighbors() if n.GetIdx() in rings[1]]
        if len(a) == 2 and len(b) == 2:
            yield atom.GetIdx(), a, b


def spiro_centres(mol) -> list[dict]:
    """'Xabab' spiro atoms: each ring's two neighbours differ (a > b), and the rings are alike (a/a', b/b'), so the
    spiro atom is a stereocentre that R/S tagging in SMILES cannot express (spiro[4.4]nonane-1,6-dione)."""
    k = _kekule(mol)
    found = []
    for s, a, b in _spiro_atoms(mol):
        ka = sorted(((_branch(k, n, s), n) for n in a), reverse=True)
        kb = sorted(((_branch(k, n, s), n) for n in b), reverse=True)
        if ka[0][0] == ka[1][0] or kb[0][0] == kb[1][0]:
            continue
        if [x[0] for x in ka] != [x[0] for x in kb]:
            continue  # Xabcd: an ordinary stereocentre, handled as tetrahedral by RDKit
        found.append({"type": "spiro", "key": f"spiro:{s}", "centre": s, "ring_a": [x[1] for x in ka],
                      "ring_b": [x[1] for x in kb], "atoms": sorted({s, *a, *b}), "bonds": [],
                      "possible_configurations": ["1", "2"],
                      "stability": "configurationally stable (a stereogenic spiro atom cannot invert)"})
    return found


def needs_new_perception(mol) -> bool:
    """An 'Xaabb' spiro system whose stereo the legacy perception drops: a spiro atom whose two neighbours are alike
    within each ring, with possible stereocentres elsewhere in both rings (2,6-dichlorospiro[3.3]heptane)."""
    k = _kekule(mol)
    potential = None
    for s, a, b in _spiro_atoms(mol):
        if _branch(k, a[0], s) != _branch(k, a[1], s) and _branch(k, b[0], s) != _branch(k, b[1], s):
            continue
        if potential is None:
            with rdBase.BlockLogs():
                potential = {st.centeredOn for st in Chem.FindPotentialStereo(Chem.Mol(mol))
                             if str(st.type) == "Atom_Tetrahedral"}
        info = mol.GetRingInfo()
        rings = [set(r) for r in info.AtomRings() if s in r]
        if all((ring - {s}) & potential for ring in rings):
            return True
    return False


def helices(mol) -> list[dict]:
    """Helicenes: at least HELIX_MIN_RINGS six-membered rings ortho-fused in a chain that always turns the same way.
    The turns are read from a flat layout of the ring system (each angular fusion turns by about 60 degrees)."""
    info = mol.GetRingInfo()
    rings = [r for r in info.AtomRings() if len(r) == 6 and all(
        mol.GetAtomWithIdx(i).GetIsAromatic() or mol.GetAtomWithIdx(i).GetHybridization() == Chem.HybridizationType.SP2
        for i in r)]
    if len(rings) < HELIX_MIN_RINGS:
        return []
    # In a helicene no ring atom belongs to three rings (the inner rim is two-ring atoms); in fullerenes, coronene and
    # other cages or sheets many do, and their ring chains are not helices.
    rings = [r for r in rings if all(info.NumAtomRings(i) <= 2 for i in r)]
    if len(rings) < HELIX_MIN_RINGS:
        return []
    shared = {(i, j) for i in range(len(rings)) for j in range(len(rings))
              if i != j and len(set(rings[i]) & set(rings[j])) == 2}
    flat = Chem.Mol(mol)
    with rdBase.BlockLogs():
        AllChem.Compute2DCoords(flat)
    xy = flat.GetConformer().GetPositions()[:, :2]
    centre = [xy[list(r)].mean(axis=0) for r in rings]

    def turn(i, j, k):
        u, v = centre[j] - centre[i], centre[k] - centre[j]
        return np.degrees(np.arctan2(u[0] * v[1] - u[1] * v[0], np.dot(u, v)))

    best = []

    def extend(path, sign):
        nonlocal best
        if len(path) > len(best):
            best = list(path)
        for nxt in range(len(rings)):
            if nxt in path or (path[-1], nxt) not in shared:
                continue
            if len(path) >= 2:
                t = turn(path[-2], path[-1], nxt)
                if abs(abs(t) - 60) > 25 or (sign and np.sign(t) != sign):
                    continue
                extend(path + [nxt], sign or np.sign(t))
            else:
                extend(path + [nxt], 0)

    for start in range(len(rings)):
        extend([start], 0)
    if len(best) < HELIX_MIN_RINGS:
        return []
    atoms = sorted(set().union(*(set(rings[i]) for i in best)))
    n = len(best)
    stability = ("configurationally stable at room temperature (racemisation barrier about 36 kcal/mol for "
                 "[6]helicene, higher for longer helicenes)" if n >= 6 else
                 "racemises slowly at room temperature ([5]helicene: barrier about 24 kcal/mol)")
    return [{"type": "helix", "key": f"helix:{atoms[0]}-{atoms[-1]}-{n}", "rings": [list(rings[i]) for i in best],
             "atoms": atoms, "bonds": [], "size": n, "possible_configurations": ["P", "M"], "stability": stability}]


def units(mol) -> list[dict]:
    """All allene, spiro (Xabab) and helix units of a molecule (graph analysis only)."""
    return allenes(mol) + spiro_centres(mol) + helices(mol)


# ---- configuration from 3D coordinates ------------------------------------------------------------------------------

def measure(mol, conf, unit) -> dict | None:
    """{'descriptor', 'helicity', 'axial', 'sense', 'verified'} for this unit in these 3D coordinates, or None."""
    if unit["type"] == "allene":
        t0, t1 = unit["ends"]
        (_, near), (_, far) = unit["near"], unit["far"]
        a, b = (_position(mol, conf, t0, x[1]) for x in near)
        c, d = (_position(mol, conf, t1, x[1]) for x in far)
        p0, p1 = np.array(conf.GetAtomPosition(t0)), np.array(conf.GetAtomPosition(t1))
        angle = torsion(a, p0, p1, c)
        if abs(angle) < 10 or abs(angle) > 170:
            return None  # flattened: no configuration to read
        helicity = "P" if angle > 0 else "M"
        axial = "Ra" if handedness(a, b, c, d) == "R" else "Sa"
        verified = (helicity == "M") == (axial == "Ra")  # IUPAC P-92.1.2.2.3
        return {"sense": helicity, "helicity": helicity, "axial": axial, "verified": verified,
                "descriptor": helicity if verified else None}
    if unit["type"] == "spiro":
        s = unit["centre"]
        points = [np.array(conf.GetAtomPosition(i)) for i in (*unit["ring_a"], *unit["ring_b"])]
        if np.linalg.norm(np.cross(points[1] - points[0], points[2] - points[0])) < 1e-3:
            return None
        sense = handedness(points[0], points[2], points[1], points[3])  # a, a', b (b' away): a fixed reference order
        return {"sense": "1" if sense == "R" else "2", "verified": False, "descriptor": None}
    if unit["type"] == "helix":
        centres = [np.mean([np.array(conf.GetAtomPosition(i)) for i in ring], axis=0) for ring in unit["rings"]]
        angles = [torsion(*centres[i:i + 4]) for i in range(len(centres) - 3)]
        if not angles or min(abs(x) for x in angles) < 2 or len({np.sign(x) for x in angles}) != 1:
            return None
        helicity = "P" if angles[0] > 0 else "M"  # right-handed helix: positive torsions (tests check a model helix)
        return {"sense": helicity, "helicity": helicity, "verified": True, "descriptor": helicity}
    return None


def read_input_configuration(mol, found=None) -> dict:
    """{unit key: sense} measured on the input's 3D coordinates (CXSMILES or MOL), stored on the molecule."""
    if mol.GetNumConformers() == 0 or not _is_3d(mol.GetConformer()):
        return {}
    found = units(mol) if found is None else found
    conf = mol.GetConformer()
    config = {}
    for unit in found:
        m = measure(mol, conf, unit)
        if m is not None:
            config[unit["key"]] = m["sense"]
    mol.SetProp(UNITS_PROP, json.dumps(config))
    return config


def stated(mol) -> dict:
    return json.loads(mol.GetProp(UNITS_PROP)) if mol.HasProp(UNITS_PROP) else {}


def _is_3d(conf) -> bool:
    z = conf.GetPositions()[:, 2]
    return bool(np.ptp(z) > 0.1)


def describe(mol) -> list[dict]:
    """The allene/spiro/helix units in the unified form (see the module docstring)."""
    config = stated(mol)
    out = []
    for unit in units(mol):
        sense = config.get(unit["key"])
        entry = {"type": unit["type"], "atoms": unit["atoms"], "bonds": unit["bonds"], "specified": sense is not None,
                 "configuration": None, "possible_configurations": unit["possible_configurations"],
                 "configuration_source": "3D coordinates in the input" if sense else None,
                 "verification": "not applicable", "descriptor": None, "stability": unit["stability"], "key": unit["key"]}
        if unit["type"] == "helix":
            entry["rings"] = unit["size"]
        if unit["type"] != "spiro" or sense is None:
            entry["label_atoms"] = list(_label_bond(unit))  # where the drawing puts its (M)/(P)/(?) label
        if sense is not None:
            if unit["type"] == "spiro":
                entry["configuration"] = f"configuration {sense}"
                entry["verification"] = "unverified"
                entry["note"] = ("This 'Xabab' spiro atom (IUPAC P-93.5.3.2) needs a CIP extension that MolTalk does "
                                 "not implement, so no R/S is given; the configuration itself is kept, drawn and built.")
            else:
                entry["configuration"] = sense
                entry["verification"] = ("verified: helicity rule and CIP axial rule agree (IUPAC P-92.1.2.2.3, M = Ra)"
                                         if unit["type"] == "allene" else
                                         "verified: handedness of the helix through the ring centres (right-handed = P)")
                entry["descriptor"] = sense
                if unit["type"] == "allene":
                    entry["alternative_descriptor"] = {"M": "Ra", "P": "Sa"}[sense]
        out.append(entry)
    return out


# ---- building molecules with a chosen configuration ----------------------------------------------------------------

def _embed(mol, count=12, seed=0xA11):
    h = Chem.AddHs(mol)
    params = AllChem.ETDG() if helices(mol) else AllChem.ETKDGv3()  # see conformer._embed
    params.randomSeed = seed
    params.timeout = 10
    with rdBase.BlockLogs():
        AllChem.EmbedMultipleConfs(h, count, params)
    return h


def _mirror(conf):
    for i in range(conf.GetNumAtoms()):
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x, p.y, -p.z))


def _other_chirality(mol) -> bool:
    """Does mirroring the coordinates change anything besides the 3D-coordinate units (R/S centres, axes)?"""
    return any(a.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
               for a in mol.GetAtoms()) or any(b.GetStereo() in (Chem.BondStereo.STEREOATROPCW,
                                                                  Chem.BondStereo.STEREOATROPCCW)
                                               for b in mol.GetBonds())


def models(mol, wanted: dict | None = None, found=None) -> dict:
    """3D models (molecules with explicit H and one conformer) for each combination of unit configurations that
    could be built: {tuple of (key, sense): model}. With `wanted`, only matching models."""
    found = units(mol) if found is None else found
    if not found:
        return {}
    embedded = _embed(mol)
    results = {}
    mirror_ok = not _other_chirality(mol)
    for conf in embedded.GetConformers():
        variants = [False, True] if mirror_ok else [False]
        for flip in variants:
            single = Chem.Mol(embedded, confId=conf.GetId())
            c = single.GetConformer()
            if flip:
                _mirror(c)
            senses = []
            for unit in found:
                m = measure(single, c, unit)
                senses.append((unit["key"], m["sense"] if m else None))
            if any(s is None for _, s in senses):
                continue
            key = tuple(senses)
            if wanted and any(wanted.get(k, s) != s for k, s in key):
                continue
            results.setdefault(key, single)
    return results


def with_coordinates(mol, model) -> "Chem.Mol":
    """The heavy-atom molecule carrying the model's 3D coordinates, so CXSMILES writes them."""
    heavy = Chem.Mol(mol)
    heavy.RemoveAllConformers()
    conf = Chem.Conformer(heavy.GetNumAtoms())
    src = model.GetConformer()
    for i in range(heavy.GetNumAtoms()):
        conf.SetAtomPosition(i, src.GetAtomPosition(i))
    conf.Set3D(True)
    heavy.AddConformer(conf, assignId=True)
    return heavy


def cxsmiles_3d(mol, model) -> str:
    """CXSMILES with the model's 3D heavy-atom coordinates, rounded; self-checked by reading it back."""
    carrier = with_coordinates(mol, model)
    params = Chem.SmilesWriteParams()
    text = Chem.MolToCXSmiles(carrier, params, Chem.CXSmilesFields.CX_COORDS)
    back = Chem.MolFromSmiles(text)
    if back is None:
        raise ValueError("Could not write this configuration as CXSMILES.")
    def senses(m):
        return sorted((key.split(":")[0], sense) for key, sense in read_input_configuration(m).items())
    if senses(back) != senses(Chem.Mol(carrier)):
        raise ValueError("The 3D configuration did not survive CXSMILES; it was not dropped silently.")
    return text


# ---- 2D depiction ---------------------------------------------------------------------------------------------------

def _pseudo_3d(drawn, lifts: dict):
    """The drawing as 3D points: everything in the paper plane, wedged ends toward the viewer (+z), hashed away."""
    conf = Chem.Conformer(drawn.GetConformer())
    for i, z in lifts.items():
        p = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, Point3D(p.x, p.y, z))
    return conf


def _bond_from(rw, start: int, end: int):
    """The bond start–end with `start` as its begin atom (a wedge's narrow end), re-adding it reversed if needed."""
    bond = rw.GetBondBetweenAtoms(start, end)
    if bond.GetBeginAtomIdx() == start:
        return bond
    order, aromatic = bond.GetBondType(), bond.GetIsAromatic()
    rw.RemoveBond(start, end)
    rw.AddBond(start, end, order)
    bond = rw.GetBondBetweenAtoms(start, end)
    bond.SetIsAromatic(aromatic)
    return bond


def _label_bond(unit) -> tuple[int, int]:
    if unit["type"] == "allene":
        return unit["ends"][0], unit["near"][0]
    if unit["type"] == "helix":
        middle = len(unit["rings"]) // 2
        shared = sorted(set(unit["rings"][middle - 1]) & set(unit["rings"][middle]))
        return shared[0], shared[1]
    return unit["centre"], unit["ring_b"][0]


def depict(drawn, mol, config: dict) -> tuple["Chem.Mol", list[dict], list[str]]:
    """Wedges and labels for allene, spiro and helix units on a prepared 2D drawing. mol is the parsed input (same
    atom numbering); config its stated configurations. Every wedge choice is checked on the pseudo-3D drawing (wedge
    toward the viewer) with the same measurement as 3D coordinates, so the drawing cannot contradict the input.
    Returns (drawing, depicted wedges, problems)."""
    found = units(mol)
    if not found:
        return drawn, [], []
    rw = Chem.RWMol(drawn)
    wedges, problems = [], []
    for unit in found:
        sense = config.get(unit["key"])
        a, b = _label_bond(unit)
        bond = rw.GetBondBetweenAtoms(a, b)
        if bond is not None:
            # Helicity of a helicene, M/P of an allene; an open unit gets "(?)". A spiro (Xabab) atom gets no label.
            if unit["type"] != "spiro" or sense is None:
                bond.SetProp("_CIPCode", sense if sense and unit["type"] != "spiro" else "?")
        if sense is None or unit["type"] == "helix":
            continue  # a helicene is drawn flat with its label: wedges would imply stereocentres it does not have
        if unit["type"] == "allene":
            centre, (_, far) = unit["ends"][1], unit["far"]
            movable = [x[1] for x in far if x[1] is not None]
        else:
            centre, movable = unit["centre"], list(unit["ring_b"])
        options = [{movable[0]: 0.8, **({movable[1]: -0.8} if len(movable) > 1 else {})},
                   {movable[0]: -0.8, **({movable[1]: 0.8} if len(movable) > 1 else {})}]
        chosen = None
        for lifts in options:
            measured = measure(mol, _pseudo_3d(rw, lifts), unit)
            if measured and measured["sense"] == sense:
                chosen = lifts
                break
        if chosen is None:
            problems.append(f"The {unit['type']} configuration could not be drawn with wedges; see the 3D view.")
            continue
        for atom, z in chosen.items():
            wedge = _bond_from(rw, centre, atom)
            wedge.SetBondDir(Chem.BondDir.BEGINWEDGE if z > 0 else Chem.BondDir.BEGINDASH)
            wedges.append({"from_atom": centre, "to_atom": atom, "style": "wedge" if z > 0 else "dash",
                           "unit": unit["type"], "configuration": sense})
    return rw.GetMol(), wedges, problems


# ---- Xaabb spiro systems (2,6-disubstituted spiro[3.3]heptanes) -----------------------------------------------------

def xaabb(mol) -> list[dict]:
    """Spiro axes of the 'Xaabb' kind, read with RDKit's newer perception (see needs_new_perception). IUPAC's
    preferred names give R/S at the spiro atom and the ring atoms that carry the substituents (RDKit labels exactly
    these), with M/P as the alternative (P-93.5.3.5). M/P is computed here only for two four-membered rings, where
    the axis runs through the spiro atom and the two opposite ring atoms: the torsion between the top-ranked
    substituents on those atoms, positive = P. Tests check IUPAC's example: (2R,4S,6R)-2,6-dichlorospiro[3.3]heptane
    is (2P)."""
    if Chem.GetUseLegacyStereoPerception():
        return []
    k = _kekule(mol)
    found = []
    for s, a, b in _spiro_atoms(mol):
        rings = [set(r) for r in mol.GetRingInfo().AtomRings() if s in r]
        if all(len(r) == 4 for r in rings):
            far = []
            for ring in rings:
                opposite = [i for i in ring if i != s and not mol.GetBondBetweenAtoms(i, s)]
                far.append(opposite[0] if opposite else None)
            if None in far:
                continue
        else:
            far = []
        atom = mol.GetAtomWithIdx(s)
        tagged = atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
        labels = {i: mol.GetAtomWithIdx(i).GetProp("_CIPCode") for i in [s, *far]
                  if mol.GetAtomWithIdx(i).HasProp("_CIPCode")}
        if not tagged and not far:
            continue
        entry = {"type": "spiro", "subtype": "Xaabb axis", "atoms": sorted({s, *far}), "bonds": [],
                 "specified": tagged and bool(labels), "configuration": None,
                 "possible_configurations": ["P", "M"] if far else ["R", "S"],
                 "configuration_source": "SMILES chirality tags (read with RDKit's newer stereo perception)" if tagged else None,
                 "verification": "not applicable", "descriptor": None,
                 "stability": "configurationally stable (a spiro axis cannot invert without breaking bonds)"}
        if entry["specified"]:
            entry["configuration"] = ",".join(f"{i}{labels[i]}" for i in sorted(labels))
            entry["descriptor"] = entry["configuration"]
            entry["verification"] = "RDKit CIP labeller (IUPAC preferred R/S form)"
            if far:
                h = Chem.AddHs(mol)
                params = AllChem.ETKDGv3()
                params.randomSeed = 0x5B
                with rdBase.BlockLogs():
                    ok = AllChem.EmbedMolecule(h, params) == 0
                if ok:
                    conf = h.GetConformer()
                    tops = []
                    for t in far:
                        ring = next(r for r in rings if t in r)
                        ligands = _ligands(_kekule(h), t, ring)
                        tops.append(_position(h, conf, t, ligands[0][1]))
                    angle = torsion(tops[0], np.array(conf.GetAtomPosition(far[0])),
                                    np.array(conf.GetAtomPosition(far[1])), tops[1])
                    if 10 < abs(angle) < 170:
                        entry["alternative_descriptor"] = "P" if angle > 0 else "M"
        found.append(entry)
    return found
