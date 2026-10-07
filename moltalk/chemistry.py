import itertools

from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, rdCIPLabeler
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from .depiction import layout
from .conformer import conformer_3d
from .stereo import stereogenic_unspecified, organic_stereo, assign_cip
from .coordination import coordinate, normalize_coordination, is_metal
from . import atropisomer, stereounits

GROUPS = {"alcohol": "[OX2H][CX4]", "phenol": "[OX2H]c", "carboxylic acid": "[CX3](=O)[OX2H]",
          "ester": "[CX3](=O)[OX2][#6]", "amide": "[CX3](=O)[NX3]", "ketone": "[#6][CX3](=O)[#6]",
          "aldehyde": "[CX3H1](=O)[#6]", "amine": "[NX3;!$(N-C=O);!$(N-S(=O)=O)]",
          "ether": "[OD2]([#6;!$([#6]=[O,S,N])])[#6;!$([#6]=[O,S,N])]", "nitrile": "[CX2]#N", "alkene": "[CX3]=[CX3]",
          "alkyne": "[CX2]#[CX2]", "thiol": "[SX2H]", "halogen": "[F,Cl,Br,I]"}
WEDGE_NOTE = ("Wedge/dash depends on this 2D drawing and bond direction; it is not an intrinsic synonym for R/S. "
              "A wedge starts at the stereocenter (narrow end) and points toward the viewer; a dash points away.")
MAX_DRAWN_ISOMERS = 16
MAX_INPUT = 16384  # a 256-atom CXSMILES with 3D coordinates needs ~8,000 characters

def parse(smiles: str):
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > MAX_INPUT:
        raise ValueError(f"Provide a nonempty SMILES string of at most {MAX_INPUT} characters "
                         "(CXSMILES with 3D coordinates included).")
    # Legacy stereo perception, except for an 'Xaabb' spiro system (2,6-dichlorospiro[3.3]heptane), whose stereo
    # only the newer perception keeps (see stereounits). The mode is set per request's molecule: workers run one
    # request at a time, and every later step of that request (canonical SMILES, 3D, export) must use the same mode.
    Chem.SetUseLegacyStereoPerception(True)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(_parse_error(smiles))
    if mol.GetNumAtoms() > 256:
        raise ValueError("V1 supports at most 256 atoms per molecule.")
    if stereounits.needs_new_perception(mol):
        Chem.SetUseLegacyStereoPerception(False)
        mol = Chem.MolFromSmiles(smiles)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    assign_cip(mol)
    atropisomer.label(mol)  # an axis keeps a descriptor only if verified (Ra/Sa); RDKit's raw P/M is not shown
    stereounits.read_input_configuration(mol)  # allenes, helicenes, Xabab spiro atoms: from 3D coordinates if given
    return mol

def _parse_error(smiles: str) -> str:
    """Explain why RDKit rejected the input without suggesting a replacement structure."""
    suffix = " No structure was drawn or guessed; supply a corrected SMILES."
    raw = Chem.MolFromSmiles(smiles, sanitize=False)
    if raw is None:
        return f"Invalid SMILES syntax: RDKit could not parse {smiles[:80]!r}.{suffix}"
    problems = [p.Message() for p in Chem.DetectChemistryProblems(raw)]
    detail = "; ".join(problems[:3]) or "sanitization failed"
    return f"Chemically invalid structure ({detail}). Atom numbers are zero-based input indices.{suffix}"

def _stereo_summary(mol, potential) -> dict:
    real = stereogenic_unspecified(mol, potential)
    ignored = sum(str(s.specified) == 'Unspecified' for s in potential) - len(real)
    real |= {("axis", i) for i in atropisomer.candidate_axes(mol)}  # RDKit does not list unspecified axes at all
    units = stereounits.describe(mol)  # allenes, helicenes, Xabab spiro atoms
    open_units = [u for u in units if not u["specified"]]
    specified = sum(str(s.specified) == 'Specified' for s in potential) + sum(u["specified"] for u in units)
    unspecified = len(real) + len(open_units)
    if not specified and not unspecified:
        status = "no stereo elements"
    elif not unspecified:
        status = "fully specified"
    elif not specified:
        status = "unspecified"
    else:
        status = "partially specified"
    summary = {"status": status, "specified": specified, "unspecified": unspecified,
               "unspecified_atoms": sorted(i for kind, i in real if kind == "atom"),
               "unspecified_bonds": sorted([mol.GetBondWithIdx(i).GetBeginAtomIdx(), mol.GetBondWithIdx(i).GetEndAtomIdx()]
                                           for kind, i in real if kind == "bond"),
               "unspecified_axes": sorted([mol.GetBondWithIdx(i).GetBeginAtomIdx(), mol.GetBondWithIdx(i).GetEndAtomIdx()]
                                          for kind, i in real if kind == "axis"),
               "unspecified_units": [{"type": u["type"], "atoms": u["atoms"]} for u in open_units]}
    if ignored:
        summary["non_stereogenic_ignored"] = ignored  # e.g. adamantane bridgeheads: flipping them changes nothing
    return summary

_UNIT_NAMES = {"allene": "allene chirality axis", "helix": "helicene (helical chirality)",
               "spiro": "stereogenic spiro atom", "axis": "hindered axis (atropisomerism)"}


def canonical(mol) -> str:
    """Canonical SMILES; CXSMILES when a configuration needs more than SMILES: 3D coordinates for allenes, helicenes
    and Xabab spiro atoms, 2D coordinates plus a wedge for hindered axes."""
    if stereounits.stated(mol):
        return stereounits.cxsmiles_3d(mol, mol)
    return atropisomer.canonical_smiles(mol)


def stereo_units(mol, potential, summary) -> list[dict]:
    """Every stereogenic unit in one form (see stereounits): tetrahedral centres, double bonds, hindered axes,
    allenes, spiro atoms, helicenes."""
    found = []
    open_atoms, open_bonds = set(summary["unspecified_atoms"]), {tuple(b) for b in summary["unspecified_bonds"]}
    for st in potential:
        kind = str(st.type)
        if kind == "Atom_Tetrahedral":
            atom = mol.GetAtomWithIdx(st.centeredOn)
            specified = str(st.specified) == "Specified"
            if not specified and st.centeredOn not in open_atoms:
                continue  # not really stereogenic (adamantane bridgeheads)
            cip = atom.GetProp("_CIPCode") if atom.HasProp("_CIPCode") else None
            found.append({"type": "tetrahedral", "atoms": [st.centeredOn], "bonds": [], "specified": specified,
                          "configuration": cip if specified else None, "possible_configurations": ["R", "S"],
                          "configuration_source": "SMILES chirality tag" if specified else None,
                          "verification": "RDKit CIP labeller" if cip and specified else "not applicable",
                          "descriptor": cip if specified else None, "stability": "configurationally stable"})
        elif kind == "Bond_Double":
            bond = mol.GetBondWithIdx(st.centeredOn)
            pair = (bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
            specified = str(st.specified) == "Specified"
            if not specified and tuple(sorted(pair)) not in open_bonds and pair not in open_bonds:
                continue
            cip = bond.GetProp("_CIPCode") if bond.HasProp("_CIPCode") else None
            found.append({"type": "double_bond", "atoms": list(pair), "bonds": [st.centeredOn], "specified": specified,
                          "configuration": cip if specified else None, "possible_configurations": ["E", "Z"],
                          "configuration_source": "SMILES bond directions" if specified else None,
                          "verification": "RDKit CIP labeller" if cip and specified else "not applicable",
                          "descriptor": cip if specified else None, "stability": "configurationally stable"})
    for axis in atropisomer.describe(mol):
        found.append({"type": "axis", "subtype": atropisomer.axis_kind(mol, axis["bond_index"]),
                      "atoms": axis["atom_indices"], "bonds": [axis["bond_index"]], "specified": True,
                      "configuration": axis["cip"], "possible_configurations": ["Ra", "Sa"],
                      "configuration_source": "CXSMILES wedge at the axis",
                      "verification": ("verified: RDKit CIP (P/M), helicity rule and CIP axial rule agree"
                                       if axis["verified"] else "unverified"),
                      "descriptor": axis["cip"], "alternative_descriptor": axis["helicity"],
                      "stability": atropisomer.stability(mol, axis["bond_index"])})
    for a, b in summary["unspecified_axes"]:
        index = mol.GetBondBetweenAtoms(a, b).GetIdx()
        found.append({"type": "axis", "subtype": atropisomer.axis_kind(mol, index), "atoms": [a, b], "bonds": [index],
                      "specified": False, "configuration": None, "possible_configurations": ["Ra", "Sa"],
                      "configuration_source": None, "verification": "not applicable", "descriptor": None,
                      "stability": atropisomer.stability(mol, index)})
    for unit in stereounits.describe(mol):
        unit = dict(unit)
        unit.pop("key", None)
        found.append(unit)
    spiro_axes = stereounits.xaabb(mol)
    if spiro_axes:  # these atoms are reported as part of the spiro unit, not as separate centres
        covered = set().union(*(set(u["atoms"]) for u in spiro_axes))
        found = [u for u in found if not (u["type"] == "tetrahedral" and u["atoms"][0] in covered)] + spiro_axes
    return found


def analyze(smiles: str) -> dict:
    mol = parse(smiles)
    potential = organic_stereo(mol, Chem.FindPotentialStereo(mol))
    stereo = [{"type": str(s.type), "center": s.centeredOn,
               "specified": str(s.specified), "descriptor": str(s.descriptor)} for s in potential]
    centers = [{"atom_index": a.GetIdx(), "element": a.GetSymbol(), "cip": a.GetProp('_CIPCode')}
               for a in mol.GetAtoms() if a.HasProp('_CIPCode') and not is_metal(a)]
    doubles = [{"bond_index": b.GetIdx(), "atom_indices": [b.GetBeginAtomIdx(), b.GetEndAtomIdx()],
                "configuration": b.GetProp('_CIPCode') if b.HasProp('_CIPCode') else str(b.GetStereo())}
               for b in mol.GetBonds() if b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() != Chem.BondStereo.STEREONONE]
    groups = {name: [list(m) for m in mol.GetSubstructMatches(Chem.MolFromSmarts(pattern))]
              for name, pattern in GROUPS.items()}
    groups["aromatic ring"] = [list(ring) for ring in mol.GetRingInfo().AtomRings()
                               if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in ring)]
    summary = _stereo_summary(mol, potential)
    warnings = []
    if summary["unspecified"]:
        warnings.append("Some stereochemistry is unspecified; no unique stereoisomer is implied.")
    if summary["unspecified_axes"]:
        warnings.append("Hindered axis (atropisomerism) with its twist unspecified at bond(s) "
                        + ", ".join(f"{a}-{b}" for a, b in summary["unspecified_axes"])
                        + ": at least three ortho positions are substituted or fused, so the Ra and Sa atropisomers are "
                          "normally separable. Plain SMILES cannot state the twist; enumerate_stereoisomers gives both "
                          "as CXSMILES.")
    axial = atropisomer.describe(mol)
    if axial:
        warnings.append("InChI/InChIKey do not encode axial chirality: both atropisomers share this InChIKey.")
    for unit in summary["unspecified_units"]:
        warnings.append(f"{_UNIT_NAMES[unit['type']].capitalize()} (atoms {', '.join(map(str, unit['atoms'][:6]))}"
                        f"{'…' if len(unit['atoms']) > 6 else ''}) with its configuration unspecified. SMILES cannot "
                        "state it; enumerate_stereoisomers gives each configuration as CXSMILES with 3D coordinates.")
    units = stereo_units(mol, potential, summary)
    if any(u["type"] in ("allene", "helix", "spiro") and u["specified"] for u in units):
        warnings.append("The configuration of the allene/helicene/spiro unit is carried by the 3D coordinates of the "
                        "CXSMILES; plain SMILES, InChI, 2D MOL and CDXML do not carry it.")
    if len(Chem.GetMolFrags(mol)) > 1:
        warnings.append("Input has disconnected fragments; properties describe the entire input.")
    if mol.GetStereoGroups():
        warnings.append("Enhanced stereo groups are present; consult CXSMILES for relative/group semantics.")
    return {"canonical_smiles": canonical(mol),
            "cxsmiles": Chem.MolToCXSmiles(mol), "inchi": Chem.MolToInchi(mol),
            "inchikey": Chem.MolToInchiKey(mol), "formula": rdMolDescriptors.CalcMolFormula(mol),
            "atom_index_base": 0, "properties": {"molecular_weight": Descriptors.MolWt(mol),
            "exact_mass": Descriptors.ExactMolWt(mol), "logp_estimate": Descriptors.MolLogP(mol),
            "tpsa_angstrom2": rdMolDescriptors.CalcTPSA(mol), "h_bond_donors": Lipinski.NumHDonors(mol),
            "h_bond_acceptors": Lipinski.NumHAcceptors(mol), "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "formal_charge": Chem.GetFormalCharge(mol)}, "stereocenters": centers,
            "double_bond_stereo": doubles, "axial_stereo": axial, "stereo_units": units, "potential_stereo": stereo,
            "stereo_summary": summary,
            "functional_groups": {k:v for k,v in groups.items() if v},
            "functional_group_note": "SMARTS motifs can overlap; this is not an exhaustive chemical classification.",
            "warnings": warnings}

def _atom_name(i: int, locants: dict | None) -> str:
    locant = (locants or {}).get(str(i))
    return f"C{locant}, atom {i}" if locant else f"atom {i}"

def _svg(mol, width: int, height: int, atom_indices: bool, locants: dict | None = None, open_atoms=(), open_bonds=(),
         open_axes=()):
    """Depict an already-parsed molecule; returns the SVG, the drawing's wedge/dash bonds and layout quality.
    open_atoms: stereocentres the input leaves unspecified; RDKit annotates them "(?)" in the same style and with the
    same collision-avoiding placement as R/S (the viewer shows them only when asked for all stereo labels)."""
    original_atoms = mol.GetNumAtoms()
    depiction = layout(mol)
    source = mol
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    # Allenes, spiro atoms and helicenes: their own wedges and labels, checked against the configuration.
    mol, unit_wedges, unit_problems = stereounits.depict(mol, source, stereounits.stated(source))
    helix = any(u["type"] == "helix" for u in stereounits.units(source))
    for i in open_atoms:
        if i < original_atoms and not mol.GetAtomWithIdx(i).HasProp("_CIPCode"):
            mol.GetAtomWithIdx(i).SetProp("_CIPCode", "?")
    for a, b in [*open_bonds, *open_axes]:  # an unspecified double bond or axis: "(?)" where E/Z or Ra/Sa would be
        bond = mol.GetBondBetweenAtoms(a, b)
        if bond is not None and not bond.HasProp("_CIPCode"):
            bond.SetProp("_CIPCode", "?")
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addAtomIndices = atom_indices and not locants
    drawer.drawOptions().addStereoAnnotation = True
    if depiction["method"] == "projection" or helix:
        # A 3D view, or a helicene's end rings drawn flat, overlap atoms by design; no red error boxes.
        drawer.drawOptions().flagCloseContactsDist = -1
    if locants:
        # IUPAC locants as atom notes; RDKit still draws its own (R)/(S) beside them.
        for atom in mol.GetAtoms():
            locant = locants.get(str(atom.GetIdx()))
            if locant:
                atom.SetProp("atomNote", locant)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    atom_px = [[round(drawer.GetDrawCoords(i).x, 2), round(drawer.GetDrawCoords(i).y, 2)] for i in range(original_atoms)]
    bonds = []
    axes = {}
    for axis in atropisomer.specified_axes(mol):
        axes[axis.GetBeginAtomIdx()] = axes[axis.GetEndAtomIdx()] = axis
    unit_pairs = {(w["from_atom"], w["to_atom"]) for w in unit_wedges}
    for b in mol.GetBonds():
        if b.GetBondDir() not in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
            continue
        if (b.GetBeginAtomIdx(), b.GetEndAtomIdx()) in unit_pairs:
            continue  # described below
        start, end = b.GetBeginAtom(), b.GetEndAtom()
        kind = "wedge" if b.GetBondDir() == Chem.BondDir.BEGINWEDGE else "dash"
        axis = axes.get(start.GetIdx())
        if axis is not None and start.GetChiralTag() == Chem.ChiralType.CHI_UNSPECIFIED and b.GetIdx() != axis.GetIdx():
            other = axis.GetOtherAtomIdx(start.GetIdx())
            cip = axis.GetProp('_CIPCode') if axis.HasProp('_CIPCode') else None
            bonds.append({"bond_index": b.GetIdx(), "from_atom": start.GetIdx(), "to_atom": end.GetIdx(),
                          "depiction": str(b.GetBondDir()), "style": kind, "axis_atoms": [start.GetIdx(), other],
                          "axial_cip": cip, "to_element": end.GetSymbol(),
                          "explanation": f"{kind.capitalize()} at the chirality axis {_atom_name(start.GetIdx(), locants)}–"
                                         f"{_atom_name(other, locants)}{f' ({cip})' if cip else ''}: this ring bond points "
                                         f"{'toward' if kind == 'wedge' else 'away from'} the viewer, so the ring on "
                                         f"{_atom_name(start.GetIdx(), locants)} is twisted out of the drawing plane "
                                         "relative to the other ring. It marks the twist of the axis, not a stereocentre."})
            continue
        cip = start.GetProp('_CIPCode') if start.HasProp('_CIPCode') else None
        end_label = f"{end.GetSymbol()} ({_atom_name(end.GetIdx(), locants)})"
        if end.GetIdx() >= original_atoms:
            end_label = f"{end.GetSymbol()} (added only for the drawing)"
        bonds.append({"bond_index": b.GetIdx(), "from_atom": start.GetIdx(), "to_atom": end.GetIdx(),
                      "depiction": str(b.GetBondDir()), "style": kind,
                      "stereocenter_cip": cip, "to_element": end.GetSymbol(),
                      "explanation": f"{kind.capitalize()} from stereocenter "
                                     f"{_atom_name(start.GetIdx(), locants)}"
                                     f"{f' ({cip})' if cip else ''} to {end_label}: in this drawing that substituent points "
                                     f"{'toward' if kind == 'wedge' else 'away from'} the viewer."})
    for w in unit_wedges:
        toward = "toward" if w["style"] == "wedge" else "away from"
        if w["unit"] == "allene":
            text = (f"{w['style'].capitalize()} at the end of the allene, {_atom_name(w['from_atom'], locants)} to "
                    f"{_atom_name(w['to_atom'], locants)}: this substituent points {toward} the viewer, so that end's "
                    f"substituents lie across the paper while the other end's lie in it; together they show the "
                    f"allene's twist ({w['configuration']}). Not a stereocentre.")
        else:
            text = (f"{w['style'].capitalize()} in a ring at the spiro atom {_atom_name(w['from_atom'], locants)}, to "
                    f"{_atom_name(w['to_atom'], locants)}: this ring bond points {toward} the viewer, so the second "
                    "ring stands across the paper (IUPAC allows wedges within rings at spiro atoms).")
        bonds.append({"from_atom": w["from_atom"], "to_atom": w["to_atom"], "style": w["style"],
                      "depiction": "BEGINWEDGE" if w["style"] == "wedge" else "BEGINDASH", "unit": w["unit"],
                      "explanation": text})
    depiction = _depiction_note(depiction)
    if helix:
        depiction.pop("warning", None)
        depiction["note"] = ("A helicene is drawn as a view of its helical shape, as in textbooks: the end rings lie over "
                             "each other in reality, and the view sets them side by side. A flat picture shows no "
                             "handedness, so it is given as P or M (see the 3D view), not with wedges."
                             if depiction.get("helicene") else
                             "A helicene is drawn flat, as in textbooks; the molecule is a helix (see the 3D view). Its "
                             "handedness is given as P or M, not with wedges.")
    if unit_problems:
        depiction["warning"] = " ".join([depiction.get("warning", ""), *unit_problems]).strip()
    return drawer.GetDrawingText(), bonds, depiction, atom_px

def _depiction_note(depiction: dict) -> dict:
    if depiction.get("helicene"):
        return depiction  # the helicene note is set by _svg
    if depiction["method"] == "projection":
        depiction["note"] = ("This cage cannot be drawn flat without bonds crossing, so it is drawn as a view of its 3D "
                             "shape, as textbooks draw cubane or adamantane: bonds that cross pass in front of or behind "
                             "each other, and lengths and angles are foreshortened. Connectivity is exact.")
        return depiction  # crossings are expected in a 3D view; no layout warning
    if depiction["bond_crossings"] or depiction["overlapping_atoms"]:
        depiction["warning"] = (f"This 2D layout has {depiction['bond_crossings']} crossing bond pair(s) and "
                                f"{depiction['overlapping_atoms']} overlapping atom pair(s), which is common for bridged "
                                "or macrocyclic structures. Connectivity and the analysis are unaffected; do not read "
                                "geometry from the crossing region.")
    return depiction

def electrons(mol) -> dict:
    """{atom index: {"pairs": n, "radicals": r}} for drawing lone pairs. Valence electrons, minus formal charge,
    minus one electron per bond (hidden hydrogens included; aromatic rings counted in a Kekulé form), minus a pair
    for each dative bond the atom donates, minus unpaired electrons, halved. Carbon is included only when charged
    or a radical; metals are skipped; an atom whose count is odd or negative is left out rather than guessed."""
    from .coordination import is_metal
    kekule = Chem.Mol(mol)
    try:
        Chem.Kekulize(kekule, clearAromaticFlags=True)
    except Exception:
        return {}
    table = Chem.GetPeriodicTable()
    found = {}
    for atom in kekule.GetAtoms():
        z, charge, radicals = atom.GetAtomicNum(), atom.GetFormalCharge(), atom.GetNumRadicalElectrons()
        if z <= 2 or is_metal(atom) or (z == 6 and not charge and not radicals):
            continue
        bonding, donated = atom.GetTotalNumHs(), 0
        for bond in atom.GetBonds():
            if bond.GetBondType() == Chem.BondType.DATIVE:
                if bond.GetBeginAtomIdx() == atom.GetIdx():
                    donated += 1  # this atom's lone pair forms the dative bond
                continue
            bonding += int(bond.GetBondTypeAsDouble())
        free = table.GetNOuterElecs(z) - charge - bonding - 2 * donated - radicals
        if free < 0 or free % 2:
            continue
        if free or radicals:
            found[str(atom.GetIdx())] = {"pairs": free // 2, "radicals": radicals, "h": atom.GetTotalNumHs(), "charge": charge,
                                         "element": atom.GetSymbol()}
    return found

def depiction_mol(smiles: str, hydrogens: bool = False):
    """The molecule as drawn: metal–N bonds added for porphyrin-type complexes, optional explicit hydrogens.
    Atom and bond indices of the parsed input are preserved (new atoms and bonds are appended)."""
    mol = parse(smiles)
    coordinated, note = coordinate(mol)
    if coordinated is None:
        coordinated, note = normalize_coordination(mol)
    if coordinated is not None:
        mol = coordinated
    if hydrogens:
        mol = Chem.AddHs(mol)
    return mol, note

def draw(smiles: str, width: int = 640, height: int = 420, atom_indices: bool = True, hydrogens: bool = False,
         locants: dict | None = None) -> dict:
    if not 200 <= width <= 1600 or not 200 <= height <= 1600:
        raise ValueError("Drawing dimensions must be between 200 and 1600 pixels.")
    mol, coordination_note = depiction_mol(smiles, hydrogens)
    electron_map = electrons(mol)
    drawn_bonds = [[b.GetBeginAtomIdx(), b.GetEndAtomIdx()] for b in mol.GetBonds()]
    analysis = analyze(smiles)
    svg, bonds, depiction, atom_px = _svg(mol, width, height, atom_indices, locants,
                                          analysis["stereo_summary"]["unspecified_atoms"],
                                          analysis["stereo_summary"]["unspecified_bonds"],
                                          analysis["stereo_summary"]["unspecified_axes"])
    if coordination_note:
        depiction["coordination_note"] = coordination_note
    return {"svg": svg, "atom_px": atom_px, "drawn_bonds": drawn_bonds, "lone_pairs": electron_map, "depicted_stereo_bonds": bonds, "depiction": depiction, "note": WEDGE_NOTE,
            "analysis": analysis}

def conformer(smiles: str, hydrogens: bool = False) -> dict:
    heavy, _ = depiction_mol(smiles, False)
    drawn, _ = depiction_mol(smiles, hydrogens) if hydrogens else (heavy, None)
    return conformer_3d(parse(smiles), heavy, drawn)

def substructure(smiles: str, smarts: str) -> dict:
    if len(smarts) > 1024:
        raise ValueError("SMARTS is limited to 1024 characters.")
    query = Chem.MolFromSmarts(smarts)
    if query is None:
        raise ValueError("Invalid SMARTS query.")
    return {"atom_index_base": 0, "matches": [list(x) for x in parse(smiles).GetSubstructMatches(query, useChirality=True, maxMatches=100)]}

_FLIP_SENSE = {"P": "M", "M": "P", "1": "2", "2": "1"}


def _isomer_key(mol, senses=None) -> str:
    """Identity of a stereoisomer: canonical isomeric SMILES, axis twists, and allene/helix/spiro configurations."""
    senses = stereounits.stated(mol) if senses is None else senses
    return atropisomer.stereo_key(mol) + "".join(f" {k}={v}" for k, v in sorted(senses.items()))


def _mirror_key(mol) -> str:
    mirror = Chem.Mol(mol)
    for atom in mirror.GetAtoms():
        if atom.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            atom.InvertChirality()
    senses = {k: _FLIP_SENSE[v] for k, v in stereounits.stated(mol).items()}
    return _isomer_key(atropisomer.mirror_axes(mirror), senses)


def _with_units(iso, open_units) -> list:
    """The isomer with every combination of configurations of the open allene/helix/spiro units that can be built
    (each as a molecule carrying matching 3D coordinates), or the isomer itself if there are none."""
    if not open_units:
        return [iso]
    wanted = stereounits.stated(iso)
    out = []
    for combo, model in stereounits.models(iso, wanted or None).items():
        carrier = stereounits.with_coordinates(iso, model)
        stereounits.read_input_configuration(carrier)
        out.append(carrier)
    return out

def enumerate_stereo(smiles: str, limit: int = 16, with_drawings: bool = False) -> dict:
    if not 1 <= limit <= 64:
        raise ValueError("limit must be between 1 and 64.")
    options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=limit + 1)
    # Enumerated molecules keep the input atom order, so indices below are input indices. RDKit does not enumerate
    # unspecified hindered axes (it cannot detect them without a wedge), so each isomer is expanded with both twists
    # of every axis atropisomer.candidate_axes finds.
    input_mol = parse(smiles)
    axes = atropisomer.candidate_axes(input_mol)
    twists = list(itertools.product(atropisomer.FLIP, repeat=len(axes)))
    stated = stereounits.stated(input_mol)
    open_units = [u for u in stereounits.units(input_mol) if u["key"] not in stated]
    if open_units and (axes or atropisomer.specified_axes(input_mol)):
        raise ValueError("This structure has both a hindered biaryl/C–N axis and an allene, helicene or spiro unit; "
                         "enumerating both together is not supported (their configurations travel in different "
                         "CXSMILES forms).")
    found = {}
    for m in EnumerateStereoisomers(Chem.Mol(input_mol), options=options):
        for twist in twists:
            iso = atropisomer.with_axes(m, dict(zip(axes, twist)))
            if stated:
                iso.SetProp(stereounits.UNITS_PROP, input_mol.GetProp(stereounits.UNITS_PROP))
            for variant in _with_units(iso, open_units):
                found.setdefault(_isomer_key(variant), variant)
        if len(found) > limit:
            break
    values = sorted(found)
    rank = {s: i for i, s in enumerate(values[:limit])}
    isomers, drawings = [], []
    for i, key in enumerate(values[:limit]):
        mol = found[key]
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        assign_cip(mol)
        axial = atropisomer.label(mol)
        mirror = _mirror_key(mol)
        isomers.append({"index": i, "smiles": canonical(mol), "inchikey": Chem.MolToInchiKey(mol),
                        "stereocenters": [{"atom_index": a.GetIdx(), "element": a.GetSymbol(), "cip": a.GetProp('_CIPCode')}
                                          for a in mol.GetAtoms() if a.HasProp('_CIPCode')],
                        "double_bond_stereo": [{"atom_indices": [b.GetBeginAtomIdx(), b.GetEndAtomIdx()], "configuration": b.GetProp('_CIPCode')}
                                               for b in mol.GetBonds() if b.HasProp('_CIPCode') and b.GetBondType() == Chem.BondType.DOUBLE],
                        "axial_stereo": axial,
                        "stereo_units": [{k: v for k, v in u.items() if k != "key"} for u in stereounits.describe(mol)],
                        "achiral": mirror == key, "enantiomer_index": rank.get(mirror) if mirror != key else None})
        if with_drawings and i < MAX_DRAWN_ISOMERS:
            drawings.append(_svg(Chem.Mol(mol), 320, 260, True)[0])
    input_mol = parse(smiles)
    result = {"input_smiles": smiles, "input_stereo_summary": _stereo_summary(input_mol, organic_stereo(input_mol, Chem.FindPotentialStereo(input_mol))),
              "isomers": isomers, "truncated": len(values) > limit, "atom_index_base": 0,
              "note": "Enumerates unspecified stereochemistry: R/S, E/Z, the twist (Ra/Sa) of hindered biaryl and C–N "
                      "axes (as CXSMILES with a wedge), and the configuration of allenes (M/P), helicenes (P/M) and "
                      "Xabab spiro atoms (as CXSMILES with 3D coordinates); plain SMILES can state none of these. "
                      "Does not certify conformational feasibility. "
                      "'achiral' marks an isomer identical to its mirror image (e.g. meso); "
                      "isomers that are neither the same nor enantiomers are diastereomers."}
    if with_drawings:
        result["svgs"] = drawings
    return result
