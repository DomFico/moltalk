"""Which unspecified stereo elements actually matter.

RDKit's FindPotentialStereo is deliberately generous: adamantane's bridgeheads are
listed as potential stereocentres although every assignment gives the same molecule.
An unspecified element is counted only if flipping it, in one concrete stereoisomer,
changes the canonical isomeric SMILES.
"""
from rdkit import Chem, rdBase
from rdkit.Chem import AllChem, rdCIPLabeler, rdMolDescriptors
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from .coordination import is_metal

MAX_ELEMENTS = 64
MAX_CAGE_ISOMERS = 256  # geometric check only up to 2^8 candidate stereoisomers
_cache: dict[str, set] = {}


def _canonical(mol) -> str:
    mol = Chem.Mol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    return Chem.MolToSmiles(mol, isomericSmiles=True)


CIP_MAX_ITERATIONS = 2_500_000  # about 2 s; RDKit: most structures need under 10,000


def assign_cip(mol, atoms=None, bonds=None) -> bool:
    """RDKit's CIP labeller, bounded. On a highly symmetric cage with many tagged centres (dodecahedrane's 3D model
    has 20) the unbounded labeller effectively never finishes; it is limited to the atoms and bonds actually needed
    and to CIP_MAX_ITERATIONS, and a structure that exceeds it simply gets no labels. False if it gave up."""
    if atoms is not None and bonds is not None and not atoms and not bonds:
        return True
    try:
        rdCIPLabeler.AssignCIPLabels(mol, atomsToLabel=atoms, bondsToLabel=bonds, maxRecursiveIterations=CIP_MAX_ITERATIONS)
        return True
    except Exception:  # noqa: BLE001 - MaxIterationsExceeded and friends: leave unlabelled rather than hang
        return False


def organic_stereo(mol, potential) -> list:
    """Potential stereo elements without metal centres: a four-coordinate metal (e.g. square-planar Fe in heme)
    is not a tetrahedral stereocentre, although RDKit lists it as one. Nor are the double bonds of an allene E/Z
    bonds (RDKit lists them): the allene is one chirality axis, handled in stereounits."""
    from .stereounits import allenes
    cumulated = set()
    for unit in allenes(mol):
        cumulated |= set(unit["bonds"])
    return [s for s in potential if not (str(s.type) == "Atom_Tetrahedral" and is_metal(mol.GetAtomWithIdx(s.centeredOn)))
            and not (str(s.type) == "Bond_Double" and s.centeredOn in cumulated)]


def _mirror(mol):
    mirror = Chem.Mol(mol)
    for atom in mirror.GetAtoms():
        if atom.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            atom.InvertChirality()
    return mirror


def _buildable_isomers(mol):
    """Stereoisomers (filling unspecified elements) that can actually be built in 3D, or None if not determinable.
    In a cage such as cubane, inverting one corner is impossible, although the SMILES would allow it."""
    params = AllChem.ETKDGv3()
    params.maxIterations = 5
    found = {}
    options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=MAX_CAGE_ISOMERS + 1)
    with rdBase.BlockLogs():
        candidates = list(EnumerateStereoisomers(Chem.Mol(mol), options=options))
        if len(candidates) > MAX_CAGE_ISOMERS:
            return None
        for iso in candidates:
            key = Chem.MolToSmiles(iso)
            if key in found:
                continue
            for seed in (7, 1234):  # a buildable isomer can fail a quick attempt; try a second seed
                params.randomSeed = seed
                if AllChem.EmbedMolecule(Chem.AddHs(iso), params) >= 0:
                    mirror = _mirror(iso)  # a buildable stereoisomer's mirror image is buildable too
                    found[key], found[Chem.MolToSmiles(mirror)] = iso, mirror
                    break
    return list(found.values()) or None  # nothing buildable is implausible: don't judge


def _cage_filter(mol, real):
    """Keep only elements whose configuration differs among the buildable stereoisomers."""
    isomers = _buildable_isomers(mol)
    if not isomers:
        return real
    labels = {}
    for iso in isomers:
        Chem.AssignStereochemistry(iso, cleanIt=True, force=True)
        assign_cip(iso)
        for kind, i in real:
            item = iso.GetAtomWithIdx(i) if kind == "atom" else iso.GetBondWithIdx(i) if kind == "bond" else None
            value = item.GetProp("_CIPCode") if item is not None and item.HasProp("_CIPCode") else None
            labels.setdefault((kind, i), set()).add(value)
    return {element for element in real if element[0] == "other" or len(labels.get(element, ())) > 1}


def stereogenic_unspecified(mol, potential) -> set[tuple[str, int]]:
    """Return {("atom", atom_index) | ("bond", bond_index)} for unspecified elements that are real."""
    key = Chem.MolToSmiles(mol)
    if key not in _cache:
        real = _flip_test(mol, potential)
        ring_info = mol.GetRingInfo()
        caged = rdMolDescriptors.CalcNumBridgeheadAtoms(mol) or any(ring_info.NumAtomRings(a.GetIdx()) >= 3 for a in mol.GetAtoms())
        if real and caged:  # bridged or cage-like (cubane): check which configurations can actually be built
            real = _cage_filter(mol, real)
        _cache[key] = real
    return set(_cache[key])


def _flip_test(mol, potential) -> set[tuple[str, int]]:
    """Unspecified elements whose flip changes the molecule (on paper; see _cage_filter for geometry)."""
    unspecified = [s for s in organic_stereo(mol, potential) if str(s.specified) == "Unspecified"][:MAX_ELEMENTS]
    if not unspecified:
        return set()
    options = StereoEnumerationOptions(onlyUnassigned=True, unique=False, maxIsomers=1)
    concrete = next(iter(EnumerateStereoisomers(Chem.Mol(mol), options=options)), None)
    if concrete is None:
        return {("atom" if str(s.type) == "Atom_Tetrahedral" else "bond", s.centeredOn) for s in unspecified}
    reference = _canonical(concrete)
    real = set()
    for s in unspecified:
        flipped = Chem.Mol(concrete)
        if str(s.type) == "Atom_Tetrahedral":
            flipped.GetAtomWithIdx(s.centeredOn).InvertChirality()
            key = ("atom", s.centeredOn)
        elif str(s.type) == "Bond_Double":
            bond = flipped.GetBondWithIdx(s.centeredOn)
            swap = {Chem.BondStereo.STEREOCIS: Chem.BondStereo.STEREOTRANS,
                    Chem.BondStereo.STEREOTRANS: Chem.BondStereo.STEREOCIS}
            if bond.GetStereo() not in swap:
                real.add(("bond", s.centeredOn))
                continue
            bond.SetStereo(swap[bond.GetStereo()])
            key = ("bond", s.centeredOn)
        else:  # atropisomers and other rare types: keep RDKit's verdict
            real.add(("other", s.centeredOn))
            continue
        if _canonical(flipped) != reference:
            real.add(key)
    return real
