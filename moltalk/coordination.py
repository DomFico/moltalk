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


def find_site(mol):
    """Return (metal index, [four donor N indices]) or None."""
    ions = [f[0] for f in Chem.GetMolFrags(mol) if len(f) == 1 and mol.GetAtomWithIdx(f[0]).GetSymbol() in METALS]
    if len(ions) != 1:
        return None
    ring_info = mol.GetRingInfo()
    donors = [a.GetIdx() for a in mol.GetAtoms()
              if a.GetSymbol() == "N" and a.GetDegree() == 2 and a.GetTotalNumHs() == 0
              and a.GetFormalCharge() in (0, -1) and ring_info.IsAtomInRingOfSize(a.GetIdx(), 5)]
    if len(donors) != 4:
        return None
    return ions[0], donors


def coordinate(mol):
    """Return (molecule with metal–N bonds, note) or (None, None) when the pattern does not apply."""
    site = find_site(mol)
    if site is None:
        return None, None
    metal, donors = site
    rw = Chem.RWMol(mol)
    ion = rw.GetAtomWithIdx(metal)
    anionic = [d for d in donors if rw.GetAtomWithIdx(d).GetFormalCharge() == -1]
    covalent = anionic if anionic and len(anionic) == ion.GetFormalCharge() else []
    for d in donors:
        if d in covalent:
            rw.AddBond(d, metal, Chem.BondType.SINGLE)
            rw.GetAtomWithIdx(d).SetFormalCharge(0)
            ion.SetFormalCharge(ion.GetFormalCharge() - 1)
        else:
            rw.AddBond(d, metal, Chem.BondType.DATIVE)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out)
    except Exception:
        return None, None
    symbol = ion.GetSymbol()
    kinds = f"{len(covalent)} covalent and {4 - len(covalent)} dative" if covalent else "4 dative"
    note = (f"{symbol}–N bonds ({kinds}) are drawn to the four ring nitrogens for this picture. The source structure "
            f"lists {symbol} as a separate ion; the analysis, formula and charges use that structure unchanged.")
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
