# Examples

## Things to ask

Ask in plain language; the assistant chooses MolTalk's tools. Naming MolTalk helps the first time ("Using MolTalk,
draw…").

| Ask | What you get |
|---|---|
| Draw caffeine and identify its functional groups. | Name resolved from the bundled library; inline drawing you can rotate in 3D; groups listed |
| Draw (1R,3S)-1-fluoro-3-methylcyclohexane and explain the wedges. | Both stereocentres verified (R and S); both carry wedges, showing that wedge ≠ R/S |
| Draw (2R)-2-chloro-4-methylhexane. | Read by OPSIN; C2 verified as R; C4 reported as unspecified |
| Show the stereoisomers of BINAP. | Both atropisomers with verified (Ra)/(Sa) descriptors, as enantiomers |
| Draw nitromethane with lone pairs. | Lewis structure with lone pairs and formal charges, in 2D and 3D |
| Draw ferrocene. | The sandwich complex assembled from PubChem's loose pieces, in its linear 3D geometry |
| Draw Wilkinson's catalyst. | Square-planar Rh(I), 16 electrons, explained |
| Draw cubane and replace four hydrogens with F, Cl, Br and I. | A cage drawn as a view of its 3D shape; stereo status of the new centres |
| Draw glucose. | Refused: "glucose" does not say which stereoisomer (D or L); asks which one |
| Draw (R)-BINAP. | Refused: the axial descriptor cannot be verified from the name; suggests enumerating BINAP's stereoisomers |
| Give me a ChemDraw file of that. | A CDXML download (MOL, SDF, PDB, XYZ and SMILES also available) |

## Files in this folder

All produced by MolTalk itself:

| File | What it is |
|---|---|
| `caffeine.cdxml` | ChemDraw XML export (2D) |
| `ra-binap.cxsmiles` | (Ra)-BINAP as CXSMILES: plain SMILES cannot state an axial twist, so it carries 2D coordinates and the axis wedge |
| `ra-binap-3d.sdf` | (Ra)-BINAP as a 3D SD file with explicit hydrogens |
| `ferrocene.svg` | The 2D drawing of ferrocene, assembled from PubChem's `[CH]1C=CC=C1.[CH]1C=CC=C1.[Fe]` |
| `fluoromethylcyclohexane.svg`, `.json` | The drawing and `depicted_stereo_bonds` for (1R,3S)-1-fluoro-3-methylcyclohexane |
