# Scope and limitations

MolTalk prefers saying "not determined" to showing something unverified. These are the boundaries of what it
currently does.

## General

- **Not a quantum-chemistry tool.** 3D models are single force-field conformers chosen to match the drawing, labelled
  as such; they are not energy minima, ensembles or measured structures. No reaction prediction, conformer energies,
  spectra or property prediction beyond RDKit descriptors (logP is an estimate).
- **Size:** up to 256 atoms per molecule. Very large or flexible molecules get fewer 3D candidates to stay within
  time limits.
- **Functional groups** come from a small SMARTS dictionary plus aromatic-ring detection; not an exhaustive
  classification.

## Names

- Names come from the bundled library (22,038 compounds) and OPSIN (systematic names). Trivial names outside the
  library need PubChem, which the model must explicitly allow and which may answer "server busy" from cloud servers.
- Ambiguous names (the same name for different structures) and bare D/L family names ("glucose") are refused with an
  explanation, by design.
- Structure → IUPAC name exists only for structures in the library or PubChem; a novel structure gets atom indices
  and no name. Locants are omitted when the parent structure is ambiguous.

## Stereochemistry

- Names with axial or helical descriptors ("(R)-BINAP", "(P)-hexahelicene") are not resolved; draw the plain compound
  and enumerate its stereoisomers instead.
- Hindered axes are detected by a stated rule of thumb (ortho substitution), not a rotation-barrier calculation.
- Xabab spiro atoms (spiro[4.4]nonane-1,6-dione) are represented, drawn and built, but get no descriptor (it needs a
  CIP extension neither RDKit nor MolTalk implements).
- Not handled: planar chirality (cyclophanes), heterohelicenes with five-membered rings, octahedral and other
  non-tetrahedral stereocentres, and enumerating a hindered axis together with an allene/helix/spiro unit.
- **Formats:** plain SMILES cannot state atropisomer, allene, helicene or Xabab-spiro configurations. MolTalk uses
  CXSMILES with coordinates for these. 2D MOL and CDXML cannot carry allene/helicene/spiro configurations, and PDB and
  XYZ have no bond orders; exports warn when information cannot be preserved.
- **InChI/InChIKey** do not encode axial chirality: both atropisomers of BINAP share one InChIKey.

## Metal complexes

- Mononuclear complexes are assembled by standard electron-counting rules; records with more than one metal
  (Pd₂(dba)₃, Tebbe's reagent, Stryker's reagent) are shown as stored, with a note.
- Not modelled: η⁶-arene complexes, bridging ligands (Pd(OAc)₂ is shown as a monomer, not its real trimer), clusters,
  and solid-state or aggregated structures.
- Geometry follows d-electron count and coordination number. Exceptions decided by ligand-field strength or spin
  state are not predicted (NiCl₂(PPh₃)₂ is really tetrahedral; the rules give square planar). Metal–ligand distances
  are estimated from covalent radii, and which ligands are cis or trans is not encoded in SMILES.
- Carbene spin states are not given by SMILES; the 3D shape assumes the common case (singlet with a heteroatom
  neighbour, otherwise triplet) and says so.

## Drawings

- A few bridged polycycles cannot be drawn flat without either a long bond or a crossing (bicyclic amanitin-type
  peptides); MolTalk keeps the long bond rather than add a crossing.
- Crowded metal complexes are drawn with long metal–ligand bonds so the ligands fit, as in textbooks.
- Cages (cubane, C₆₀) and helicenes are drawn as views of their 3D shape, so some bonds cross by design.

## Hosts

- Inline rendering depends on the host's MCP Apps support. Clients without UI support get the full data and can ask
  for the SVG (`include_svg=true`).
- Hosts cache the viewer template; after an update, refresh the plugin or connector and start a new chat.
- File downloads use the host's download API where available, otherwise a copy/save panel.
