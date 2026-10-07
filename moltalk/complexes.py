"""Metal complexes stored as loose pieces, assembled for the drawing and the 3D model.

PubChem, and so MolTalk's library, stores most metal complexes as a bare metal atom plus separate ligands and
counter-ions: Wilkinson's catalyst is "[Cl-].[Rh].PPh3.PPh3.PPh3", ferrocene "[CH]1C=CC=C1.[CH]1C=CC=C1.[Fe]".
Drawn as stored, that is a free metal atom beside floating ligands; built in 3D, each piece lands anywhere. Here the
complex is put together, for the picture only (the analysis, formula and charges always use the stored structure),
with general rules rather than per-compound fixes:

1. Each piece is classified: a counter-ion (K+, Na+, H+, NH4+, PF6-, BF4-, SbF6-, ClO4-, BPh4- ...), a ligand with
   donor atoms, or neither (water of crystallisation, solvent). Donors: phosphines/arsines, CO and isocyanide carbon,
   cyanide, pyridine-type and amine nitrogen, amide anions, halides, O- of carboxylates, alkoxides, phenoxides and
   sulfonates, thiolates, hydride, carbanions; carbonyl and ether oxygen only inside a chelate (acac, salen). A piece
   with several donors binds through all of them that can close 5- or 6-membered chelate rings (bipyridine N,N;
   salen N,N,O,O; oxalate O,O; diamines). pi ligands: cyclopentadienyl (eta5) and isolated C=C (eta2: ethylene,
   COD), drawn as a bond to the ring or alkene centre.
2. Ligands bind in order of donor strength (chelates and Cp first, then neutral L donors, then anionic X donors)
   while the metal's valence electron count stays within 18 (16 for d8 Ni/Pd/Pt/Cu/Ag/Au, which are square planar),
   and six positions. Ru(bpy)3Cl2: three bipyridines make Ru(II) 18-electron, so the chlorides stay counter-ions;
   Wilkinson's RhCl(PPh3)3 is a 16-electron d8 complex with the chloride bound. The oxidation state follows the
   ionic model: the metal's charge, plus ligands already bonded to it, plus bound radical ligands (Cp, Cl as
   written by PubChem), plus any negative charge a record leaves unbalanced (PubChem writes Rh(I) in Wilkinson's
   catalyst as neutral Rh beside Cl-).
3. Bonds: neutral donors dative (arrow to the metal); anionic and radical ligands covalent when the metal's
   oxidation state accounts for them, so the textbook picture has neutral atoms; otherwise dative, charges kept.
4. Geometry for 3D: from the electron count and the number of positions (a Cp ring or an alkene is one position):
   2 linear; 3 trigonal planar; 4 square planar for d8, else tetrahedral; 5 square pyramidal with a planar
   tetradentate chelate (salen), else trigonal bipyramidal; 6 octahedral.

Not assembled: records with more than one metal (Pd2(dba)3, Tebbe's reagent, Stryker's reagent; their bonding is
not a matter of simple rules), and arene (eta6) complexes. They are drawn as stored, and the result says so.
"""
import json

from rdkit import Chem

from .coordination import METALS

PROP = "_moltalk_complex"
CENTRE_PROP = "_eta_atoms"

GROUP = {s: g for g, symbols in {3: "Sc Y La", 4: "Ti Zr Hf", 5: "V Nb Ta", 6: "Cr Mo W", 7: "Mn Tc Re",
                                 8: "Fe Ru Os", 9: "Co Rh Ir", 10: "Ni Pd Pt", 11: "Cu Ag Au", 12: "Zn Cd Hg"}.items()
         for s in symbols.split()}
COUNTER_CATIONS = {"Li", "Na", "K", "Rb", "Cs"}
SQUARE_PLANAR_D8 = {"Ni", "Pd", "Pt", "Cu", "Ag", "Au"}
CHELATE_PATH = (3, 4)  # bonds between two donors of one ligand: 5- or 6-membered chelate ring with the metal


def _is_counter_ion(mol, frag) -> bool:
    atoms = [mol.GetAtomWithIdx(i) for i in frag]
    if len(atoms) == 1 and (atoms[0].GetSymbol() in COUNTER_CATIONS or
                            (atoms[0].GetSymbol() == "H" and atoms[0].GetFormalCharge() == 1)):
        return True
    if len(atoms) == 1 and atoms[0].GetSymbol() in ("Mg", "Ca", "Sr", "Ba") and atoms[0].GetFormalCharge() > 0:
        return True
    if len(atoms) == 1 and atoms[0].GetSymbol() == "N" and atoms[0].GetFormalCharge() == 1:
        return True  # ammonium
    for a in atoms:  # non-coordinating anions: PF6-, BF4-, SbF6-, AsF6-, BPh4-, ClO4-
        if a.GetFormalCharge() < 0 and a.GetSymbol() in ("P", "B", "Sb", "As", "Al") and a.GetDegree() >= 4:
            return True
        if a.GetSymbol() == "Cl" and sum(n.GetSymbol() == "O" for n in a.GetNeighbors()) == 4:
            return True
    return False


def _donors(mol, frag) -> list[dict]:
    """Donor atoms of one piece: {'atom', 'kind' L/X/weak, 'strength', 'charge' anion/radical/neutral}."""
    found = []
    for i in frag:
        a = mol.GetAtomWithIdx(i)
        s, q, rad, deg = a.GetSymbol(), a.GetFormalCharge(), a.GetNumRadicalElectrons(), a.GetDegree()
        form = "anion" if q < 0 else "radical" if rad else "neutral"
        entry = None
        if s in ("F", "Cl", "Br", "I") and deg == 0 and (q == -1 or rad):
            entry = ("X", 3)
        elif s == "H" and deg == 0 and q == -1:
            entry = ("X", 3)
        elif s == "C" and q == -1 and any(b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds()):
            entry = ("X", 4) if any(n.GetSymbol() == "N" and n.GetFormalCharge() == 0 for n in a.GetNeighbors()) \
                else ("L", 5)  # cyanide is X; CO and isocyanides ([C-]#[O+], [C-]#[N+]R) are L
            form = "anion" if entry[0] == "X" else "neutral"
        elif s == "C" and (q == -1 or rad) and not a.GetIsAromatic() and deg <= 1 and not a.IsInRing():
            entry = ("X", 3)  # methyl and other carbanions
        elif s in ("P", "As", "Sb") and q == 0 and not rad and deg + a.GetTotalNumHs() == 3:
            entry = ("L", 5)
        elif s == "N" and q == 0 and not rad:
            amide = any(n.GetSymbol() == "C" and any(b.GetBondType() == Chem.BondType.DOUBLE and
                                                     b.GetOtherAtom(n).GetSymbol() in ("O", "S") for b in n.GetBonds())
                        for n in a.GetNeighbors())
            if a.GetIsAromatic() and a.GetTotalNumHs() == 0 and deg == 2:
                entry = ("L", 4)  # pyridine-type
            elif any(b.GetBondType() == Chem.BondType.TRIPLE for b in a.GetBonds()):
                entry = ("L", 3)  # nitrile
            elif not a.GetIsAromatic() and not amide and a.GetHybridization() == Chem.HybridizationType.SP3:
                entry = ("L", 4)  # ammonia, amines
            elif any(b.GetBondType() == Chem.BondType.DOUBLE for b in a.GetBonds()) and not amide:
                entry = ("L", 3)  # imine (salen)
        elif s == "N" and q == -1 and _disconnected_amine(a):
            entry = ("L", 4)
            form = "amine"  # see _disconnected_amine
        elif s == "N" and q == -1:
            entry = ("X", 4)
        elif s == "O" and q == -1 and deg <= 1:
            entry = ("X", 3)
        elif s == "S" and q == -1 and deg <= 1:
            entry = ("X", 3)
        elif s == "S" and q == 0 and deg == 2 and not a.GetIsAromatic():
            entry = ("L", 2)
        elif s == "O" and q == 0 and not rad and (deg == 1 and any(b.GetBondType() == Chem.BondType.DOUBLE
                                                                   for b in a.GetBonds())):
            entry = ("weak", 1)  # carbonyl: only as part of a chelate (acac)
        if entry:
            found.append({"atom": i, "kind": entry[0], "strength": entry[1], "charge": form})
    return found


def _disconnected_amine(atom) -> bool:
    """An 'amide anion' on saturated carbon ([NH-]R, [NH2-]) beside a metal ion is how PubChem's metal disconnection
    writes a coordinated amine (oxaliplatin: [Pt+4] and two [NH-]): a real free amide anion is far too basic to sit in
    a salt. Such nitrogens are rebuilt as neutral amine donors and the metal charge lowered to match."""
    return atom.GetTotalNumHs() >= 1 and all(n.GetSymbol() == "C" and n.GetHybridization() == Chem.HybridizationType.SP3
                                             and not n.GetIsAromatic() for n in atom.GetNeighbors())


def _cyclopentadienyl(mol, frag):
    """The five carbons of a Cp-type ring (anionic or radical, all sp2) in this piece, or None."""
    for ring in mol.GetRingInfo().AtomRings():
        if len(ring) != 5 or not set(ring) <= set(frag):
            continue
        atoms = [mol.GetAtomWithIdx(i) for i in ring]
        if any(a.GetSymbol() != "C" for a in atoms):
            continue
        marked = [a for a in atoms if a.GetFormalCharge() == -1 or a.GetNumRadicalElectrons() == 1]
        # The other four carbons form two C=C (or are aromatic); RDKit calls the radical carbon itself sp3.
        sp2 = all(a.GetIsAromatic() or a.GetHybridization() == Chem.HybridizationType.SP2 for a in atoms if a not in marked)
        if len(marked) == 1 and sp2:
            centre = marked[0]
            return list(ring), centre.GetIdx(), "anion" if centre.GetFormalCharge() < 0 else "radical"
    return None


def _alkenes(mol, frag):
    """Isolated, non-aromatic C=C bonds of a piece with no other donor (ethylene, COD, norbornadiene)."""
    out = []
    for b in mol.GetBonds():
        if b.GetBondType() != Chem.BondType.DOUBLE or b.GetIsAromatic():
            continue
        i, j = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
        if i in frag and j in frag and b.GetBeginAtom().GetSymbol() == b.GetEndAtom().GetSymbol() == "C" \
                and not b.GetIsConjugated():
            out.append([i, j])
    return out


def _chelate(mol, donors):
    """The donors of one piece that bind together: the largest set (up to four) linked by 5/6-ring chelate paths,
    with at least one real (non-weak) donor; else the single strongest real donor."""
    real = [d for d in donors if d["kind"] != "weak"]
    if not real:
        return []
    dist = Chem.GetDistanceMatrix(mol)
    best = [max(real, key=lambda d: d["strength"])]
    for start in real:
        group = [start]
        for d in sorted(donors, key=lambda d: -d["strength"]):
            if d in group or len(group) >= 4:
                continue
            if any(dist[d["atom"]][g["atom"]] in CHELATE_PATH for g in group) and \
                    all(dist[d["atom"]][g["atom"]] >= 3 for g in group):
                group.append(d)
        if len(group) > len(best):
            best = group
    return best


def assemble(mol):
    """(drawn molecule, note) for a complex stored as pieces, or (None, None). Atom and bond indices of the input
    are kept; ring/alkene centres for pi ligands are appended as dummy atoms (drawn without a label)."""
    frags = [list(f) for f in Chem.GetMolFrags(mol)]
    metal_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetSymbol() in METALS and a.GetSymbol() not in COUNTER_CATIONS]
    if len(frags) < 2 or not metal_atoms:
        return None, None
    if len(metal_atoms) != 1:
        return None, ("This record lists several metal atoms as separate pieces; their bonding cannot be assigned by "
                      "simple rules, so it is drawn as stored (pieces not joined).")
    metal = metal_atoms[0]
    m = mol.GetAtomWithIdx(metal)
    symbol = m.GetSymbol()
    metal_frag = next(f for f in frags if metal in f)
    prebonded = m.GetDegree()  # ligands already bonded to the metal in its own piece (Cl-Pt-Cl, Cl-Zr-Cl)
    total_charge = sum(a.GetFormalCharge() for a in mol.GetAtoms())
    deficit = max(0, -total_charge)  # PubChem writes Rh(I) as neutral Rh beside Cl-: an unbalanced record

    # Candidate ligands, strongest first.
    candidates = []
    for frag in frags:
        if frag is metal_frag or _is_counter_ion(mol, frag):
            continue
        cp = _cyclopentadienyl(mol, frag)
        if cp:
            ring, centre, form = cp
            candidates.append({"type": "eta5", "atoms": ring, "anchor": centre, "charge": form, "kind": "X",
                               "sites": 3, "electrons": 6, "positions": 1, "priority": 0})
            continue
        donors = _donors(mol, frag)
        chosen = _chelate(mol, donors) if donors else []
        if chosen:
            xs = [d for d in chosen if d["kind"] == "X"]
            candidates.append({"type": "sigma", "donors": chosen, "sites": len(chosen), "electrons": 2 * len(chosen),
                               "positions": len(chosen), "kind": "X" if xs else "L",
                               "priority": 1 if len(chosen) > 1 else (2 if chosen[0]["kind"] == "L" else 3),
                               "strength": -max(d["strength"] for d in chosen)})
            continue
        alkenes = _alkenes(mol, frag)
        if alkenes and len(frag) <= 30:
            candidates.append({"type": "eta2", "alkenes": alkenes[:2], "sites": len(alkenes[:2]),
                               "electrons": 2 * len(alkenes[:2]), "positions": len(alkenes[:2]), "kind": "L",
                               "priority": 1 if len(alkenes) > 1 else 4})
    if not candidates:
        return None, None
    candidates.sort(key=lambda c: (c["priority"], c.get("strength", 0)))

    group = GROUP.get(symbol)
    amines = sum(d["charge"] == "amine" for c in candidates for d in c.get("donors", []))
    ox = max(0, m.GetFormalCharge() - amines) + prebonded + deficit
    bound, sites, positions, electrons = [], prebonded, prebonded, 2 * prebonded
    for c in candidates:
        radicals = 1 if c["type"] == "eta5" and c["charge"] == "radical" else \
            sum(d["charge"] == "radical" for d in c.get("donors", []))
        new_ox = ox + radicals
        if group is not None:
            d_count = group - new_ox
            limit = 16 if (d_count == 8 and symbol in SQUARE_PLANAR_D8) else 18
            if d_count < 0 or d_count + electrons + c["electrons"] > limit:
                continue
        elif sites + c["sites"] > 4:
            continue
        if positions + c["positions"] > 6 or sites + c["sites"] > 9:
            continue  # six positions around the metal (a Cp ring is one), nine sites (it is three) as in Cp2ZrCl2
        bound.append(c)
        ox, sites, positions, electrons = new_ox, sites + c["sites"], positions + c["positions"], electrons + c["electrons"]
    if not bound:
        return None, None

    d_count = None if group is None else group - ox
    geometry = _geometry(positions, d_count, symbol, bound)
    rw = Chem.RWMol(mol)
    metal_atom = rw.GetAtomWithIdx(metal)
    for c in bound:
        for d in c.get("donors", []):
            if d["charge"] == "amine":  # PubChem's disconnected amine: the neutral NH2 donor again
                atom = rw.GetAtomWithIdx(d["atom"])
                atom.SetFormalCharge(0)
                atom.SetNumExplicitHs(atom.GetTotalNumHs() + 1)
                atom.SetNoImplicit(True)
                if metal_atom.GetFormalCharge() > 0:
                    metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)
                d["charge"] = "neutral"
    # Anionic ligands bind covalently when the metal's positive charge or the record's missing charge accounts for
    # each of them (all or nothing, as for heme); radical ligands always bind covalently.
    anionic = [d for c in bound for d in c.get("donors", []) if d["charge"] == "anion"] + \
        [c for c in bound if c["type"] == "eta5" and c["charge"] == "anion"]
    # As many anions bind covalently as the metal's charge (and any charge the record leaves unbalanced) accounts
    # for, strongest first; the rest bind dative and keep their charge (Zeise's salt: two Pt–Cl bonds and one
    # Cl- -> Pt, i.e. the [PtCl3(C2H4)]- anion).
    available = max(0, metal_atom.GetFormalCharge()) + deficit
    budget = {"left": available}

    def take_anion():
        if budget["left"] > 0:
            budget["left"] -= 1
            return True
        return False
    covalent_anions = False
    n_covalent = n_dative = 0
    centres = []
    for c in bound:
        if c["type"] == "sigma":
            for d in c["donors"]:
                atom = rw.GetAtomWithIdx(d["atom"])
                if d["charge"] == "radical" or (d["charge"] == "anion" and take_anion()):
                    rw.AddBond(d["atom"], metal, Chem.BondType.SINGLE)
                    if d["charge"] == "radical":
                        atom.SetNumRadicalElectrons(atom.GetNumRadicalElectrons() - 1)
                    else:
                        atom.SetFormalCharge(atom.GetFormalCharge() + 1)
                        if metal_atom.GetFormalCharge() > 0:
                            metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)
                    n_covalent += 1
                else:
                    rw.AddBond(d["atom"], metal, Chem.BondType.DATIVE)
                    n_dative += 1
                    if atom.GetSymbol() == "C" and atom.GetFormalCharge() == -1:
                        # Carbonyl and isocyanide ligands are written M–C≡O and M–C≡N–R, without the C-/O+ charges
                        # of the free molecule: the textbook convention once the carbon's lone pair is bonding.
                        partner = next((n for n in atom.GetNeighbors() if n.GetFormalCharge() == 1), None)
                        if partner is not None:
                            atom.SetFormalCharge(0)
                            atom.SetNoImplicit(True)
                            partner.SetFormalCharge(0)
        else:
            groups = [c["atoms"]] if c["type"] == "eta5" else c["alkenes"]
            if c["type"] == "eta5":
                anchor = rw.GetAtomWithIdx(c["anchor"])
                if c["charge"] == "radical":
                    anchor.SetNumRadicalElectrons(0)
                elif take_anion():
                    anchor.SetFormalCharge(0)
                    if metal_atom.GetFormalCharge() > 0:
                        metal_atom.SetFormalCharge(metal_atom.GetFormalCharge() - 1)
                if anchor.GetIsAromatic() or c["charge"] == "radical":
                    anchor.SetNoImplicit(True)
            for atoms in groups:
                dummy = rw.AddAtom(Chem.Atom(0))
                rw.GetAtomWithIdx(dummy).SetProp(CENTRE_PROP, json.dumps(atoms))
                rw.GetAtomWithIdx(dummy).SetProp("_displayLabel", "")
                rw.GetAtomWithIdx(dummy).SetProp("_displayLabelW", "")
                if c["type"] == "eta5":
                    rw.AddBond(metal, dummy, Chem.BondType.SINGLE)
                else:
                    rw.AddBond(dummy, metal, Chem.BondType.DATIVE)  # the alkene donates: the arrow points to the metal
                centres.append({"centre": dummy, "atoms": atoms, "hapto": len(atoms)})
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        out.UpdatePropertyCache(strict=False)
    except Exception:
        return None, None
    # Sanitising re-derives radicals from valence, and a Cp carbon bonded only to the ring centre looks one bond
    # short: clear the radical the ligand gave up to the metal again.
    for c in bound:
        if c["type"] == "eta5":
            out.GetAtomWithIdx(c["anchor"]).SetNumRadicalElectrons(0)
        for d in c.get("donors", []):  # likewise a carbonyl carbon whose lone pair now bonds to the metal
            atom = out.GetAtomWithIdx(d["atom"])
            if atom.GetNumRadicalElectrons() and not mol.GetAtomWithIdx(d["atom"]).GetNumRadicalElectrons():
                atom.SetNumRadicalElectrons(0)
    info = {"metal": metal, "geometry": geometry, "d_electrons": d_count, "oxidation_state": ox,
            "electron_count": None if d_count is None else d_count + electrons,
            "donors": [d["atom"] for c in bound for d in c.get("donors", [])], "centres": centres,
            "chelates": [[d["atom"] for d in c["donors"]] for c in bound if c["type"] == "sigma" and len(c["donors"]) > 1]}
    out.SetProp(PROP, json.dumps(info))
    left = [f for f in frags if f is not metal_frag and not any(
        set(f) & set(d["atom"] for d in c.get("donors", [])) or set(f) & set(c.get("atoms", []))
        or set(f) & {a for pair in c.get("alkenes", []) for a in pair} for c in bound)]
    parts = [f"{n_covalent} covalent" if n_covalent else "", f"{n_dative} dative" if n_dative else "",
             f"{len(centres)} to a ring or alkene centre (η{'/'.join(sorted({str(c['hapto']) for c in centres}))})"
             if centres else ""]
    note = (f"The source lists {symbol} and its ligands as separate pieces; for this picture they are joined: "
            f"{', '.join(p for p in parts if p)} bond(s) to {symbol}, chosen by electron counting "
            f"({symbol}({_roman(ox)})"
            + (f", d{d_count}, {d_count + electrons} electrons" if d_count is not None else "")
            + f", {geometry.replace('_', ' ')})")
    if left:
        note += f"; {len(left)} other piece(s) (counter-ions, solvent) stay separate"
    note += ". The analysis, formula and charges use the source unchanged."
    return out, note


def _geometry(positions, d_count, symbol, bound) -> str:
    if positions <= 1:
        return "terminal"
    if positions == 2:
        return "linear"
    if positions == 3:
        return "trigonal_planar"
    if positions == 4:
        return "square_planar" if d_count in (8, 9) else "tetrahedral"  # d9 Cu(II): Jahn–Teller, square planar
    if positions == 5:
        # Trigonal bipyramid for d8/d10 (Fe(CO)5) and main-group centres; square pyramid otherwise (Grubbs' d6 Ru,
        # Jacobsen's Mn with its planar salen).
        planar_tetradentate = any(c["type"] == "sigma" and len(c["donors"]) == 4 for c in bound)
        return "square_pyramidal" if planar_tetradentate or (d_count is not None and d_count not in (8, 10)) \
            else "trigonal_bipyramidal"
    return "octahedral"


def _roman(n: int) -> str:
    return {0: "0", 1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI", 7: "VII", 8: "VIII"}.get(n, str(n))


def info(mol) -> dict | None:
    return json.loads(mol.GetProp(PROP)) if mol.HasProp(PROP) else None


# ---- 3D arrangement --------------------------------------------------------------------------------------------------

import itertools  # noqa: E402

import numpy as np  # noqa: E402

_T = 1 / np.sqrt(3)
TEMPLATES = {
    "terminal": [(0, 0, 1)],
    "linear": [(0, 0, 1), (0, 0, -1)],
    "trigonal_planar": [(1, 0, 0), (-0.5, np.sqrt(3) / 2, 0), (-0.5, -np.sqrt(3) / 2, 0)],
    "square_planar": [(1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)],
    "tetrahedral": [(_T, _T, _T), (_T, -_T, -_T), (-_T, _T, -_T), (-_T, -_T, _T)],
    "trigonal_bipyramidal": [(0, 0, 1), (0, 0, -1), (1, 0, 0), (-0.5, np.sqrt(3) / 2, 0), (-0.5, -np.sqrt(3) / 2, 0)],
    "square_pyramidal": [(0, 0, 1), (1, 0, 0), (0, 1, 0), (-1, 0, 0), (0, -1, 0)],
    "octahedral": [(1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)],
}


def _kabsch(p, q):
    """Rotation r minimising |p @ r.T - q| (proper rotation)."""
    u, _, vt = np.linalg.svd(p.T @ q)
    d = np.sign(np.linalg.det(vt.T @ u.T)) or 1.0
    return vt.T @ np.diag([1, 1, d]) @ u.T


def _rotation_between(a, b):
    a, b = a / np.linalg.norm(a), b / np.linalg.norm(b)
    v, c = np.cross(a, b), float(np.dot(a, b))
    if np.linalg.norm(v) < 1e-9:
        if c > 0:
            return np.eye(3)
        axis = np.cross(a, [1, 0, 0]) if abs(a[0]) < 0.9 else np.cross(a, [0, 1, 0])
        axis /= np.linalg.norm(axis)
        return 2 * np.outer(axis, axis) - np.eye(3)
    k = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + k + k @ k * ((1 - c) / np.linalg.norm(v) ** 2)


def _distance(mol, metal, atom) -> float:
    table = Chem.GetPeriodicTable()
    return table.GetRcovalent(mol.GetAtomWithIdx(metal).GetAtomicNum()) + \
        table.GetRcovalent(mol.GetAtomWithIdx(atom).GetAtomicNum()) - 0.1


def arrange(mol, conf_id, info) -> None:
    """Put the ligands of an assembled complex (see assemble) around the metal in its ideal geometry, in place.
    mol is the embedded molecule (input atom order, explicit hydrogens appended). Each ligand keeps its own embedded
    shape and is moved rigidly: the template directions are rotated to the ligands' current directions (best
    assignment), each piece is placed with its donor(s) at a typical metal–ligand distance (a chelate by a rigid fit
    of its donors), pi ligands with their ring or alkene centre on the axis (a Cp ring perpendicular to it, an alkene
    side-on), and monodentate ligands are then spun about their bond to the metal to avoid each other."""
    metal = info["metal"]
    conf = mol.GetConformer(conf_id)
    pos = conf.GetPositions()
    centre = pos[metal].copy()
    # Ligand pieces: the molecule without the metal's own bonds, split into connected pieces.
    cut = Chem.RWMol(mol)
    for bond in list(mol.GetAtomWithIdx(metal).GetBonds()):
        cut.RemoveBond(bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
    pieces = [set(f) for f in Chem.GetMolFrags(cut) if metal not in f]
    prebonded = [n.GetIdx() for n in mol.GetAtomWithIdx(metal).GetNeighbors() if n.GetAtomicNum() > 1]
    contacts = []  # (piece, [donor atoms] or None, [pi atom groups])
    for piece in pieces:
        donors = sorted({d for d in info["donors"] + prebonded if d in piece})
        groups = [c["atoms"] for c in info["centres"] if set(c["atoms"]) <= piece]
        if donors or groups:
            contacts.append((piece, donors, groups))
    points = []  # one entry per coordination position: (contact index, donor atom or pi group)
    for k, (piece, donors, groups) in enumerate(contacts):
        points += [(k, ("donor", d)) for d in donors] + [(k, ("pi", g)) for g in groups]
    template = np.array(TEMPLATES.get(info["geometry"], TEMPLATES["octahedral"]), dtype=float)
    if len(points) != len(template):
        return  # e.g. a position count the template does not cover: leave the embedding as it is

    def where(point):
        kind, ref = point
        return pos[ref] if kind == "donor" else pos[list(ref)].mean(axis=0)

    current = np.array([where(p[1]) - centre for p in points])
    # An embedding of loose pieces can put the metal right on a ligand atom (a zero direction): give such a position
    # an arbitrary direction; the best template assignment then decides.
    fallback = np.array(TEMPLATES["octahedral"], dtype=float)
    current = np.array([v / np.linalg.norm(v) if np.linalg.norm(v) > 1e-3 else fallback[k % 6]
                        for k, v in enumerate(current)])
    best = None
    for perm in itertools.permutations(range(len(template))):
        t = template[list(perm)]
        rot = _kabsch(t, current)
        err = float(np.sum((t @ rot.T - current) ** 2))
        if best is None or err < best[0]:
            best = (err, t @ rot.T)
    directions = best[1]
    new = pos.copy()
    for k, (piece, donors, groups) in enumerate(contacts):
        mine = [(i, p[1]) for i, p in enumerate(points) if p[0] == k]
        atoms = sorted(piece)
        if len(mine) == 1:
            i, (kind, ref) = mine[0]
            u = directions[i]
            anchor = where((kind, ref))
            if kind == "donor":
                # The ligand's bulk points away from the metal, along its bond: the vector from the donor to the
                # piece's centre is turned onto the bond direction (C=O of a carbonyl outward, not the O on the
                # metal; the loose embedding's donor-to-metal direction means nothing).
                outward = pos[atoms].mean(axis=0) - anchor
                if np.linalg.norm(outward) < 1e-3:  # a single atom (Cl): nothing to orient
                    outward = u
                rot = _rotation_between(outward, u)
                target = centre + u * _distance(mol, metal, ref)
            else:
                ring = pos[list(ref)]
                radius = Chem.GetPeriodicTable().GetRcovalent(mol.GetAtomWithIdx(metal).GetAtomicNum())
                if len(ref) >= 3:  # Cp: ring normal along the metal–centroid axis (Fe 1.7 Å, Zr 2.2 Å out)
                    _, _, vt = np.linalg.svd(ring - ring.mean(axis=0))
                    normal = vt[2] if np.dot(vt[2], anchor - centre) >= 0 else -vt[2]
                    rot = _rotation_between(normal, u)
                    target = centre + u * (radius + 0.4)
                else:  # alkene side-on: C=C perpendicular to the metal–alkene axis and to the coordination plane
                    other = next((directions[j] for j in range(len(directions)) if j != i
                                  and abs(np.dot(directions[j], u)) < 0.9), np.cross(u, [1, 0, 0]))
                    want = np.cross(u, other)
                    want /= np.linalg.norm(want) or 1.0
                    rot = _rotation_between(ring[1] - ring[0], want)
                    target = centre + u * (radius + 0.65)  # Pt–ethylene midpoint about 2.0 Å
            for a in atoms:
                new[a] = target + (pos[a] - anchor) @ rot.T
        else:  # chelate: rigid fit of its donor positions onto their target positions
            src = np.array([where(p) for _, p in mine])
            dst = np.array([centre + directions[i] * (_distance(mol, metal, p[1]) if p[0] == "donor" else 2.1)
                            for i, p in mine])
            rot = _kabsch(src - src.mean(axis=0), dst - dst.mean(axis=0))
            for a in atoms:
                new[a] = dst.mean(axis=0) + (pos[a] - src.mean(axis=0)) @ rot.T
    # Spin monodentate ligands about their metal bond to keep clear of the others.
    single = [(k, mine[0]) for k in range(len(contacts))
              for mine in [[p for p in points if p[0] == k]] if len(mine) == 1 and mine[0][1][0] == "donor"]
    for k, (_, (_, donor)) in single * 3:  # a few passes: turning one ligand changes the room the others have
        atoms = sorted(contacts[k][0])
        others = np.array([new[a] for a in range(mol.GetNumAtoms()) if a not in contacts[k][0] and a != metal])
        if not len(others) or len(atoms) < 3:
            continue
        axis = new[donor] - centre
        axis /= np.linalg.norm(axis)
        rel = new[atoms] - new[donor]
        best = None
        for angle in np.radians(np.arange(0, 360, 15)):
            c, s = np.cos(angle), np.sin(angle)
            kx = np.array([[0, -axis[2], axis[1]], [axis[2], 0, -axis[0]], [-axis[1], axis[0], 0]])
            r = np.eye(3) + s * kx + (1 - c) * kx @ kx
            trial = new[donor] + rel @ r.T
            # Clearance from everything else, measured as the closest few contacts (one tight contact used to win
            # over a ring lying across a neighbouring chloride).
            d = np.sort(np.linalg.norm(trial[:, None] - others[None], axis=2).ravel())[:4]
            gap = float(d.mean())
            if best is None or gap > best[0]:
                best = (gap, trial)
        for a, p in zip(atoms, best[1]):
            new[a] = p
    # Counter-ions and solvent (Zeise's K+ and water) are not bonded: set them clear of the complex, each in the most
    # open direction, instead of wherever the embedding left them (on top of the metal).
    bonded = {metal} | set().union(*(c[0] for c in contacts)) if contacts else {metal}
    loose = [p for p in pieces if not p & bonded]
    if loose:
        placed = new[sorted(bonded)]
        reach = float(np.max(np.linalg.norm(placed - centre, axis=1))) if len(placed) else 0.0
        options = [np.array(v, dtype=float) for v in itertools.product((-1, 0, 1), repeat=3) if any(v)]
        options = [v / np.linalg.norm(v) for v in options]
        for piece in loose:
            atoms = sorted(piece)
            rel = new[atoms] - new[atoms].mean(axis=0)
            best = max(options, key=lambda v: np.min(np.linalg.norm(placed - (centre + v * (reach + 3.0)), axis=1)))
            for a, r in zip(atoms, rel):
                new[a] = centre + best * (reach + 3.0) + r
            placed = np.vstack([placed, new[atoms]])
    for i, p in enumerate(new):
        conf.SetAtomPosition(i, p.tolist())


def embedding_source(drawn_heavy, heavy_count: int):
    """The assembled complex as one molecule for ETKDG: the drawn bonds to the metal (dative ones as single bonds),
    without the pi ligands' centre atoms. Embedded bonded, a multidentate chelate (salen) wraps around the metal; as
    loose pieces it came out as a free ligand whose donors could not all reach the metal. None if not usable."""
    rw = Chem.RWMol(drawn_heavy)
    for i in range(rw.GetNumAtoms() - 1, heavy_count - 1, -1):
        rw.RemoveAtom(i)
    for k in range(rw.GetNumBonds()):
        bond = rw.GetBondWithIdx(k)
        if bond.GetBondType() == Chem.BondType.DATIVE:
            bond.SetBondType(Chem.BondType.SINGLE)
    out = rw.GetMol()
    try:
        Chem.SanitizeMol(out, Chem.SanitizeFlags.SANITIZE_ALL ^ Chem.SanitizeFlags.SANITIZE_PROPERTIES)
        out.UpdatePropertyCache(strict=False)
    except Exception:  # noqa: BLE001
        return None
    return out


def annotate_connected(mol):
    """Geometry information for a complex stored already bonded (Grubbs' catalyst, Speier's [PtCl6]2-): one metal
    with two or more bonds. Neutral P, As, N, O or S donors that keep their normal valence without the metal bond are
    counted as dative (L) ligands even when written with a plain single bond; halides, O, S, N anions and other
    single-bonded atoms as X; a metal=C double bond as a neutral carbene (L, the usual count for Grubbs, Ru(II)).
    Sets the info on the molecule (structure unchanged) and returns it, or None."""
    from .coordination import is_metal
    every = [a for a in mol.GetAtoms() if is_metal(a) and a.GetSymbol() not in COUNTER_CATIONS]
    metals = [a for a in every if a.GetDegree() >= 2]
    if len(every) != 1 or len(metals) != 1:
        return None  # several metals (Tebbe's Ti and Al): left as stored
    m = metals[0]
    table = Chem.GetPeriodicTable()
    x = 0
    donors = []
    for bond in m.GetBonds():
        n = bond.GetOtherAtom(m)
        donors.append(n.GetIdx())
        if bond.GetBondType() in (Chem.BondType.DATIVE, Chem.BondType.DOUBLE):
            continue
        normal = table.GetDefaultValence(n.GetAtomicNum())
        own = sum(b.GetBondTypeAsDouble() for b in n.GetBonds() if b.GetOtherAtomIdx(n.GetIdx()) != m.GetIdx()) \
            + n.GetTotalNumHs()
        if n.GetSymbol() in ("P", "As", "N", "O", "S") and n.GetFormalCharge() == 0 and own >= normal:
            continue  # a dative donor written as a single bond
        x += 1
    group = GROUP.get(m.GetSymbol())
    ox = max(0, m.GetFormalCharge() + x)
    d_count = None if group is None else group - ox
    positions = len(donors)
    geometry = _geometry(positions, d_count, m.GetSymbol(), [])
    info = {"metal": m.GetIdx(), "geometry": geometry, "d_electrons": d_count, "oxidation_state": ox,
            "electron_count": None, "donors": donors, "centres": [], "chelates": [], "connected": True}
    mol.SetProp(PROP, json.dumps(info))
    return info
