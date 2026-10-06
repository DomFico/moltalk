"""Draw metal–nitrogen bonds in porphyrin-type complexes (depiction only).

Databases such as PubChem often store heme, chlorophyll or metal phthalocyanines as a
charged macrocycle plus a separate metal ion ("...[N-]...[N-]... .[Fe+2]"). RDKit draws
exactly that: a floating ion. For the drawing only, the ion is bonded to the four
five-membered-ring nitrogens. When the anionic nitrogens exactly balance the metal's
charge (heme, chlorophyll) those two become covalent bonds and the other two dative,
the textbook picture with neutral atoms. Otherwise all four are dative and every formal
charge stays as the database gives it. The analysis always uses the unmodified structure.
"""
from rdkit import Chem

METALS = set("Li Be Na Mg Al K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn "
             "Cs Ba La Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi".split())


def _macrocycle_nitrogens(mol, exclude=()):
    """Nitrogens of the five-membered rings fused into a macrocycle (a ring of 12+ atoms): the four inward-facing
    donors of a porphyrin, phthalocyanine or corrin. Not other ring nitrogens (B12's benzimidazole)."""
    # Ring perception without the metal's bonds: once drawn, the M–N bonds create small chelate rings and RDKit's
    # smallest ring set no longer contains the macrocycle itself.
    bare = Chem.RWMol(mol)
    for bond in list(bare.GetBonds()):
        if bond.GetBeginAtom().GetSymbol() in METALS or bond.GetEndAtom().GetSymbol() in METALS:
            bare.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    bare = bare.GetMol()
    bare.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(bare)
    rings = [set(r) for r in bare.GetRingInfo().AtomRings()]
    macro = [r for r in rings if len(r) >= 12]
    found = set()
    for ring in rings:
        if len(ring) == 5 and any(len(ring & m) >= 2 for m in macro):
            found |= {i for i in ring if mol.GetAtomWithIdx(i).GetSymbol() == "N" and i not in exclude}
    return sorted(found)


def find_site(mol):
    """(metal index, [four donor N indices], [axial ligand atom indices]) or None.
    The metal is a separate piece: a bare ion ([Ni+2]) or a metal with its own small ligands (Cl-Sn-Cl, as PubChem
    stores tin phthalocyanine dichloride). Axial ligands stored as separate pieces (B12's cyanide, a chloride ion) are
    attached to the metal for the drawing."""
    frags = Chem.GetMolFrags(mol)
    metal_frags = [f for f in frags if len(f) <= 7 and sum(mol.GetAtomWithIdx(i).GetSymbol() in METALS for i in f) == 1]
    if len(metal_frags) != 1:
        return None
    metal = next(i for i in metal_frags[0] if mol.GetAtomWithIdx(i).GetSymbol() in METALS)
    donors = [n for n in _macrocycle_nitrogens(mol)
              if mol.GetAtomWithIdx(n).GetDegree() == 2 and mol.GetAtomWithIdx(n).GetFormalCharge() in (0, -1)]
    if len(donors) != 4:
        return None
    axial = []
    for f in frags:
        if f is metal_frags[0] or len(f) > 2:
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in f]
        if any(a.GetSymbol() in METALS for a in atoms):
            continue
        # The ligand's binding atom: the charged or radical one (carbon of cyanide).
        binding = [a for a in atoms if a.GetFormalCharge() < 0 or a.GetNumRadicalElectrons()]
        if binding:
            axial.append(min(binding, key=lambda a: (a.GetSymbol() != "C", a.GetIdx())).GetIdx())
    return metal, donors, axial


def coordinate(mol):
    """Return (molecule with metal bonds, note) or (None, None) when the pattern does not apply.
    Bonds are covalent while the metal has charge or valence to spare (anionic ligands first, then anionic or N-H
    ring nitrogens, whose H the metal replaces as in the real complex), and dative otherwise."""
    site = find_site(mol)
    if site is None:
        return None, None
    metal, donors, axial = site
    rw = Chem.RWMol(mol)
    ion = rw.GetAtomWithIdx(metal)
    symbol = ion.GetSymbol()
    valences = [v for v in Chem.GetPeriodicTable().GetValenceList(symbol) if v > 0]
    spare = max(0, ion.GetFormalCharge())
    if ion.GetFormalCharge() == 0 and valences:
        spare = max(0, max(valences) - ion.GetExplicitValence())
    # All or nothing: bonds become covalent only when the metal's spare charge (or valence) is matched exactly by
    # the anionic or radical ligands plus the anionic and N-H ring nitrogens (heme: Fe2+ and two N-; B12: Co2+, the
    # cyanide radical and one N-H; tin dichloride: two spare valences and two N-H). Otherwise every bond is dative
    # and each charge stays as the database gives it.
    anionic = [d for d in donors if rw.GetAtomWithIdx(d).GetFormalCharge() == -1]
    nh = [d for d in donors if d not in anionic and rw.GetAtomWithIdx(d).GetTotalNumHs() > 0]
    match = spare > 0 and len(axial) + len(anionic) + len(nh) == spare
    covalent, dative, removed_h = [], [], 0
    for a in axial:
        atom = rw.GetAtomWithIdx(a)
        if match:
            rw.AddBond(a, metal, Chem.BondType.SINGLE)
            if atom.GetNumRadicalElectrons():
                atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() - 1)
            elif atom.GetFormalCharge() < 0:
                atom.SetFormalCharge(atom.GetFormalCharge() + 1)
            if ion.GetFormalCharge() > 0:
                ion.SetFormalCharge(ion.GetFormalCharge() - 1)
        else:
            rw.AddBond(a, metal, Chem.BondType.DATIVE)
    for d in donors:
        atom = rw.GetAtomWithIdx(d)
        if match and (d in anionic or d in nh):
            rw.AddBond(d, metal, Chem.BondType.SINGLE)
            if d in anionic:
                atom.SetFormalCharge(0)
            else:
                atom.SetNumExplicitHs(0)
                atom.SetNoImplicit(True)
                removed_h += 1
            if ion.GetFormalCharge() > 0:
                ion.SetFormalCharge(ion.GetFormalCharge() - 1)
            covalent.append(d)
        else:
            rw.AddBond(d, metal, Chem.BondType.DATIVE)
            dative.append(d)
    out = rw.GetMol()
    try:
        # Metals are hypervalent here (tin: two chlorides and four N; RDKit counts dative bonds), so skip only the
        # valence check; everything else is still sanitised.
        Chem.SanitizeMol(out, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        out.UpdatePropertyCache(strict=False)
    except Exception:
        return None, None
    kinds = f"{len(covalent)} covalent and {len(dative)} dative" if covalent else "4 dative"
    note = (f"{symbol}–N bonds ({kinds}) are drawn to the four ring nitrogens for this picture")
    if removed_h:
        note += f", replacing {removed_h} N–H hydrogen{'s' if removed_h > 1 else ''} as in the complex"
    if axial:
        note += f", and {len(axial)} separately listed ligand{'s are' if len(axial) > 1 else ' is'} bonded to {symbol}"
    note += (f". The source structure lists {symbol} as a separate piece; the analysis, formula and charges use that "
             "structure unchanged.")
    return out, note


def is_metal(atom) -> bool:
    return atom.GetSymbol() in METALS


def normalize_coordination(mol):
    """Charge-separated metal bonds (e.g. heme written as [Fe-2] bonded to two [N+]) become dative N→M bonds with
    neutral atoms, the same picture coordinate() gives for the ionic form. Returns (molecule, note) or (None, None)."""
    rw = Chem.RWMol(mol)
    changed = []
    for metal in [a.GetIdx() for a in mol.GetAtoms() if is_metal(a) and a.GetFormalCharge() < 0]:
        for nb in list(rw.GetAtomWithIdx(metal).GetNeighbors()):
            m = rw.GetAtomWithIdx(metal)
            bond = rw.GetBondBetweenAtoms(metal, nb.GetIdx())
            if m.GetFormalCharge() >= 0 or nb.GetSymbol() != "N" or nb.GetFormalCharge() != 1 \
                    or bond.GetBondType() != Chem.BondType.SINGLE:
                continue
            rw.RemoveBond(metal, nb.GetIdx())
            rw.AddBond(nb.GetIdx(), metal, Chem.BondType.DATIVE)
            rw.GetAtomWithIdx(nb.GetIdx()).SetFormalCharge(0)
            m.SetFormalCharge(m.GetFormalCharge() + 1)
            changed.append(m.GetSymbol())
    if not changed:
        return None, None
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None, None
    symbol = changed[0]
    note = (f"{len(changed)} charge-separated {symbol}–N⁺ bond(s) in the source are drawn as dative N→{symbol} bonds with "
            "neutral atoms (the usual way to draw coordination). The analysis, formula and charges use the source unchanged.")
    return out, note
