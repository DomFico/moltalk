"""Structure files for other programs: ChemDraw (CDXML), MDL MOL/SDF, PDB, XYZ, SMILES. All written by RDKit.

The file holds the structure as given (input atom order, charges, stereo), not drawing-only additions such as the
metal–N bonds drawn for heme. 2D coordinates are the drawing's own layout; 3D coordinates are the conformer the
viewer rotates, with explicit hydrogens. A 3D file fixes one configuration at every stereocentre, so when the input
leaves some unspecified, the result says which ones were chosen arbitrarily.
"""
import re
from io import StringIO

from rdkit import Chem
from rdkit.Chem import AllChem, rdChemDraw

from .chemistry import conformer, parse
from .depiction import layout

# format: (extension, MIME type, description, coordinates allowed)
FORMATS = {
    "cdxml": ("cdxml", "chemical/x-cdxml", "ChemDraw XML", ("2d",)),
    "mol": ("mol", "chemical/x-mdl-molfile", "MDL Molfile", ("2d", "3d")),
    "sdf": ("sdf", "chemical/x-mdl-sdfile", "MDL SD file", ("2d", "3d")),
    "pdb": ("pdb", "chemical/x-pdb", "Protein Data Bank", ("3d",)),
    "xyz": ("xyz", "chemical/x-xyz", "XYZ coordinates", ("3d",)),
    "smiles": ("smi", "chemical/x-daylight-smiles", "SMILES", ("none",)),
}
MAX_EMBED_BYTES = 200_000


def _filename(name: str | None, ext: str) -> str:
    stem = re.sub(r"[^A-Za-z0-9._()+,-]+", "_", (name or "").strip())[:60].strip("._") or "molecule"
    return f"{stem}.{ext}"


def _with_3d(smiles: str):
    """The input structure with explicit hydrogens and the viewer's conformer, plus the arbitrarily chosen
    stereo elements. Falls back to a plain ETKDG/MMFF embedding when the viewer's model does not map onto it."""
    mol = Chem.AddHs(parse(smiles))
    model = conformer(smiles, True)
    atoms = model["atoms"]
    arbitrary = {"atoms": model.get("arbitrary_stereo_atoms", []), "bonds": model.get("arbitrary_stereo_bonds", [])}
    if len(atoms) == mol.GetNumAtoms() and all(a["element"] == mol.GetAtomWithIdx(i).GetSymbol() for i, a in enumerate(atoms)):
        conf = Chem.Conformer(mol.GetNumAtoms())
        for i, a in enumerate(atoms):
            conf.SetAtomPosition(i, a["xyz"])
        conf.Set3D(True)
        mol.RemoveAllConformers()
        mol.AddConformer(conf, assignId=True)
        method = model["method"]
    else:  # e.g. a metal complex whose drawing adds bonds: embed the structure as given
        params = AllChem.ETKDGv3()
        params.randomSeed = 0xF00D
        if AllChem.EmbedMolecule(mol, params) != 0:
            raise ValueError("RDKit could not generate 3D coordinates for this structure; export it in 2D instead.")
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            method = "ETKDGv3 + MMFF94"
        else:
            AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            method = "ETKDGv3 + UFF"
    # Stereo now follows the coordinates (which matches the input wherever the input specifies it).
    Chem.AssignStereochemistryFrom3D(mol)
    return mol, method, arbitrary


def export_structure(smiles: str, fmt: str = "cdxml", coordinates: str | None = None, name: str | None = None,
                     hydrogens: bool = False) -> dict:
    fmt = (fmt or "").lower().strip().lstrip(".")
    fmt = {"chemdraw": "cdxml", "cdx": "cdx", "molfile": "mol", "sd": "sdf", "smi": "smiles"}.get(fmt, fmt)
    if fmt == "cdx":
        raise ValueError("Binary ChemDraw (.cdx) is not available; export 'cdxml', which ChemDraw opens directly.")
    if fmt not in FORMATS:
        raise ValueError(f"Unknown format {fmt!r}; choose one of: {', '.join(FORMATS)}.")
    ext, mime, description, allowed = FORMATS[fmt]
    coordinates = (coordinates or allowed[0]).lower()
    if coordinates not in allowed:
        raise ValueError(f"{description} files support {' or '.join(allowed)} coordinates, not {coordinates!r}.")
    mol = parse(smiles)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    title = (name or canonical)[:80]
    warnings, notes = [], []

    if fmt == "smiles":
        text = f"{canonical}\t{name}\n" if name else canonical + "\n"
    elif coordinates == "2d":
        if hydrogens:
            mol = Chem.AddHs(mol)
        depiction = layout(mol)  # the same deterministic layout the drawing uses
        if depiction.get("method") == "projection":
            notes.append("2D coordinates are a projection of the cage's 3D shape, as in the drawing (bonds cross).")
        mol.SetProp("_Name", title)
        if fmt == "cdxml":
            text = rdChemDraw.MolToChemDrawBlock(mol)
        elif fmt == "mol":
            text = Chem.MolToMolBlock(mol)
        else:
            text = _sdf(mol, canonical, name)
    else:
        mol, method, arbitrary = _with_3d(smiles)
        mol.SetProp("_Name", title)
        notes.append(f"3D coordinates: one calculated conformer ({method}), with explicit hydrogens; "
                     "not a measured or unique structure.")
        if arbitrary["atoms"] or arbitrary["bonds"]:
            parts = ([f"atom(s) {', '.join(map(str, arbitrary['atoms']))}"] if arbitrary["atoms"] else []) + \
                    [f"double bond {a}-{b}" for a, b in arbitrary["bonds"]]
            warnings.append(f"The input leaves stereochemistry unspecified at {'; '.join(parts)}. A 3D file must fix "
                            "a configuration there, so this file contains one arbitrary choice; it is not implied by "
                            "the input. Specify the stereochemistry, or export in 2D, if it matters.")
        if fmt == "mol":
            text = Chem.MolToMolBlock(mol)
        elif fmt == "sdf":
            text = _sdf(mol, canonical, name, warnings)
        elif fmt == "pdb":
            text = Chem.MolToPDBBlock(mol, flavor=4)  # flavor 4: CONECT records for every bond
        else:
            text = Chem.MolToXYZBlock(mol)
    if fmt in ("pdb", "xyz"):
        notes.append(f"{description} files store atoms and coordinates only; bond orders and charges are "
                     + ("not stored." if fmt == "xyz" else "only partly stored (PDB has no bond orders)."))
    data = text.encode()
    return {"kind": "file", "format": fmt, "format_name": description, "filename": _filename(name, ext),
            "mime_type": mime, "coordinates": coordinates, "explicit_hydrogens": coordinates == "3d" or hydrogens,
            "size_bytes": len(data), "canonical_smiles": canonical, "label": name, "warnings": warnings, "notes": notes,
            "text": text}


def _sdf(mol, canonical: str, name: str | None, warnings=()) -> str:
    out = StringIO()
    writer = Chem.SDWriter(out)
    mol.SetProp("SMILES", canonical)
    if name:
        mol.SetProp("NAME", name)
    if warnings:
        mol.SetProp("MOLTALK_WARNING", " ".join(warnings))
    writer.write(mol)
    writer.close()
    return out.getvalue()
