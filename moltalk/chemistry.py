from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, rdMolDescriptors, rdCIPLabeler
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions
from .depiction import layout
from .conformer import conformer_3d
from .stereo import stereogenic_unspecified, organic_stereo
from .coordination import coordinate, normalize_coordination, is_metal

GROUPS = {"alcohol": "[OX2H][CX4]", "phenol": "[OX2H]c", "carboxylic acid": "[CX3](=O)[OX2H]",
          "ester": "[CX3](=O)[OX2][#6]", "amide": "[CX3](=O)[NX3]", "ketone": "[#6][CX3](=O)[#6]",
          "aldehyde": "[CX3H1](=O)[#6]", "amine": "[NX3;!$(N-C=O);!$(N-S(=O)=O)]",
          "ether": "[OD2]([#6;!$([#6]=[O,S,N])])[#6;!$([#6]=[O,S,N])]", "nitrile": "[CX2]#N", "alkene": "[CX3]=[CX3]",
          "alkyne": "[CX2]#[CX2]", "thiol": "[SX2H]", "halogen": "[F,Cl,Br,I]"}
WEDGE_NOTE = ("Wedge/dash depends on this 2D drawing and bond direction; it is not an intrinsic synonym for R/S. "
              "A wedge starts at the stereocenter (narrow end) and points toward the viewer; a dash points away.")
MAX_DRAWN_ISOMERS = 16

def parse(smiles: str):
    if not isinstance(smiles, str) or not smiles.strip() or len(smiles) > 4096:
        raise ValueError("Provide a nonempty SMILES string of at most 4096 characters.")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(_parse_error(smiles))
    if mol.GetNumAtoms() > 256:
        raise ValueError("V1 supports at most 256 atoms per molecule.")
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
    rdCIPLabeler.AssignCIPLabels(mol)
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
    specified = sum(str(s.specified) == 'Specified' for s in potential)
    unspecified = len(real)
    ignored = sum(str(s.specified) == 'Unspecified' for s in potential) - unspecified
    if not specified and not unspecified:
        status = "no stereo elements"
    elif not unspecified:
        status = "fully specified"
    elif not specified:
        status = "unspecified"
    else:
        status = "partially specified"
    summary = {"status": status, "specified": specified, "unspecified": unspecified}
    if ignored:
        summary["non_stereogenic_ignored"] = ignored  # e.g. adamantane bridgeheads: flipping them changes nothing
    return summary

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
    if len(Chem.GetMolFrags(mol)) > 1:
        warnings.append("Input has disconnected fragments; properties describe the entire input.")
    if mol.GetStereoGroups():
        warnings.append("Enhanced stereo groups are present; consult CXSMILES for relative/group semantics.")
    return {"canonical_smiles": Chem.MolToSmiles(mol, isomericSmiles=True),
            "cxsmiles": Chem.MolToCXSmiles(mol), "inchi": Chem.MolToInchi(mol),
            "inchikey": Chem.MolToInchiKey(mol), "formula": rdMolDescriptors.CalcMolFormula(mol),
            "atom_index_base": 0, "properties": {"molecular_weight": Descriptors.MolWt(mol),
            "exact_mass": Descriptors.ExactMolWt(mol), "logp_estimate": Descriptors.MolLogP(mol),
            "tpsa_angstrom2": rdMolDescriptors.CalcTPSA(mol), "h_bond_donors": Lipinski.NumHDonors(mol),
            "h_bond_acceptors": Lipinski.NumHAcceptors(mol), "rotatable_bonds": Lipinski.NumRotatableBonds(mol),
            "formal_charge": Chem.GetFormalCharge(mol)}, "stereocenters": centers,
            "double_bond_stereo": doubles, "potential_stereo": stereo,
            "stereo_summary": summary,
            "functional_groups": {k:v for k,v in groups.items() if v},
            "functional_group_note": "SMARTS motifs can overlap; this is not an exhaustive chemical classification.",
            "warnings": warnings}

def _atom_name(i: int, locants: dict | None) -> str:
    locant = (locants or {}).get(str(i))
    return f"C{locant}, atom {i}" if locant else f"atom {i}"

def _svg(mol, width: int, height: int, atom_indices: bool, locants: dict | None = None):
    """Depict an already-parsed molecule; returns the SVG, the drawing's wedge/dash bonds and layout quality."""
    original_atoms = mol.GetNumAtoms()
    depiction = layout(mol)
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
    drawer.drawOptions().addAtomIndices = atom_indices and not locants
    drawer.drawOptions().addStereoAnnotation = True
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
    for b in mol.GetBonds():
        if b.GetBondDir() not in (Chem.BondDir.BEGINWEDGE, Chem.BondDir.BEGINDASH):
            continue
        start, end = b.GetBeginAtom(), b.GetEndAtom()
        kind = "wedge" if b.GetBondDir() == Chem.BondDir.BEGINWEDGE else "dash"
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
    return drawer.GetDrawingText(), bonds, _depiction_note(depiction), atom_px

def _depiction_note(depiction: dict) -> dict:
    if depiction["method"] == "schlegel":
        depiction["note"] = ("The cage is shown as if viewed through one ring, which becomes the outer "
                             "ring; bond lengths and ring sizes are distorted by the projection, but connectivity is exact. "
                             "RDKit's standard layout was rejected because its bonds crossed or stretched.")
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
    svg, bonds, depiction, atom_px = _svg(mol, width, height, atom_indices, locants)
    if coordination_note:
        depiction["coordination_note"] = coordination_note
    return {"svg": svg, "atom_px": atom_px, "drawn_bonds": drawn_bonds, "lone_pairs": electron_map, "depicted_stereo_bonds": bonds, "depiction": depiction, "note": WEDGE_NOTE,
            "analysis": analyze(smiles)}

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

def _mirror_smiles(mol) -> str:
    mirror = Chem.Mol(mol)
    for atom in mirror.GetAtoms():
        if atom.GetChiralTag() in (Chem.ChiralType.CHI_TETRAHEDRAL_CW, Chem.ChiralType.CHI_TETRAHEDRAL_CCW):
            atom.InvertChirality()
    return Chem.MolToSmiles(mirror, isomericSmiles=True)

def enumerate_stereo(smiles: str, limit: int = 16, with_drawings: bool = False) -> dict:
    if not 1 <= limit <= 64:
        raise ValueError("limit must be between 1 and 64.")
    options = StereoEnumerationOptions(onlyUnassigned=True, unique=True, maxIsomers=limit + 1)
    # Enumerated molecules keep the input atom order, so indices below are input indices.
    found = {}
    for m in EnumerateStereoisomers(parse(smiles), options=options):
        found.setdefault(Chem.MolToSmiles(m, isomericSmiles=True), m)
    values = sorted(found)
    rank = {s: i for i, s in enumerate(values[:limit])}
    isomers, drawings = [], []
    for i, smi in enumerate(values[:limit]):
        mol = found[smi]
        Chem.AssignStereochemistry(mol, cleanIt=True, force=True)
        rdCIPLabeler.AssignCIPLabels(mol)
        mirror = _mirror_smiles(mol)
        isomers.append({"index": i, "smiles": smi, "inchikey": Chem.MolToInchiKey(mol),
                        "stereocenters": [{"atom_index": a.GetIdx(), "element": a.GetSymbol(), "cip": a.GetProp('_CIPCode')}
                                          for a in mol.GetAtoms() if a.HasProp('_CIPCode')],
                        "double_bond_stereo": [{"atom_indices": [b.GetBeginAtomIdx(), b.GetEndAtomIdx()], "configuration": b.GetProp('_CIPCode')}
                                               for b in mol.GetBonds() if b.HasProp('_CIPCode')],
                        "achiral": mirror == smi, "enantiomer_index": rank.get(mirror) if mirror != smi else None})
        if with_drawings and i < MAX_DRAWN_ISOMERS:
            drawings.append(_svg(Chem.Mol(mol), 320, 260, True)[0])
    input_mol = parse(smiles)
    result = {"input_smiles": smiles, "input_stereo_summary": _stereo_summary(input_mol, organic_stereo(input_mol, Chem.FindPotentialStereo(input_mol))),
              "isomers": isomers, "truncated": len(values) > limit, "atom_index_base": 0,
              "note": "Enumerates unspecified stereochemistry; does not certify conformational feasibility. "
                      "'achiral' marks an isomer identical to its mirror image (e.g. meso); "
                      "isomers that are neither the same nor enantiomers are diastereomers."}
    if with_drawings:
        result["svgs"] = drawings
    return result
