# Chemistry

How MolTalk decides what a structure is, what its stereochemistry is, and how it is drawn and built in 3D. The
guiding rule throughout: **report what can be verified, and fail closed otherwise.** A name that does not identify
one structure, or a descriptor that cannot be checked, is reported as such; nothing is guessed.

- [Names and validation](#names-and-validation)
- [Stereochemistry](#stereochemistry)
- [Lone pairs, charges and radicals](#lone-pairs-charges-and-radicals)
- [2D depiction](#2d-depiction)
- [3D models](#3d-models)
- [Metal complexes](#metal-complexes)
- [IUPAC names and numbering](#iupac-names-and-numbering)

Atom indices in every result are zero-based **input-SMILES** indices, not IUPAC locants and not canonical-SMILES
order.

## Names and validation

### Name → structure

Tried in order; the first that answers wins.

1. **The bundled library** (offline). 22,038 compounds with about 200,000 names; see [DATA_SOURCES.md](DATA_SOURCES.md).
   Lookup ignores case, spacing and dash style. Names come in tiers, most trusted first: curated overrides and the
   seed list; the compound's own Wikidata label, PubChem title or IUPAC name; Wikidata aliases; PubChem synonyms. The
   first tier that knows a name decides, because lower tiers contain errors (Wikidata lists "ozone" as an alias of
   phencyclidine). A name that still points to different structures (for example "lye": NaOH or KOH) is an error
   that lists the candidates.
2. **OPSIN** (offline) reads systematic IUPAC names, including R/S and E/Z, such as `(2R)-4-chloro-2-methylheptan-3-one`.
   The result says the structure came from parsing the name and warns about stereo the name leaves open.
3. **PubChem**, only when the model passes `allow_network=true`. Best effort: on Cloud Run's shared egress addresses
   PubChem often answers "server busy".
4. Otherwise: "No structure was assumed".

**Curated overrides** cover names where database naming is messy: "heme" is heme b, Fe(II) protoporphyrin IX (PubChem
titles one copy of it "Hemin"); "hemin" is the Fe(III) chloride; also NAD⁺/NADH, NADP⁺/NADPH, FAD/FADH₂, FMN, CoA,
acetyl-CoA, ATP/ADP/AMP and cAMP, each with a note (for example on phosphate protonation).

**Stereochemistry not stated.** A bare name is refused when it names a D/L family: "alanine", "serine", "glucose" and
"fructose" ask which stereoisomer is meant, while "L-alanine" and "D-glucose" resolve. Named natural products and
drugs ("morphine") name one stereoisomer and resolve. Racemic or unspecified records ("lactic acid", "ibuprofen")
resolve with a warning that stereochemistry is unspecified.

**Stereochemistry stated in a name is verified.** If a name states R/S, E/Z or axial/helical descriptors, the
structure must encode at least as many of each, with matching CIP labels; otherwise the name is rejected. "(R)-carvone"
passes only if the centre really is R. Axial and helical descriptors in names ("(R)-BINAP", "(P)-hexahelicene") are
never resolved, because name sources give these compounds without their configuration; the plain name resolves and
`enumerate_stereoisomers` shows both forms with verified descriptors. Optical rotation, (+)/(−), is not checked.

### Label check

When the model draws a SMILES under a compound name (`label=`), MolTalk looks the label up in the library and reports
`matches`, `stereo differs`, `mismatch` (the drawing is a different compound, shown prominently) or `unverified`. A
SMILES written from memory under a known name is caught.

### Structure validation

RDKit parses and sanitises every input. Invalid syntax or impossible valences are errors with RDKit's reason; no
corrected structure is substituted.

## Stereochemistry

Every stereogenic unit is reported in one form, `analysis.stereo_units`, with `type`, `atoms`, `specified`,
`configuration`, `possible_configurations`, `configuration_source`, `verification`, `descriptor` and `stability`.
Three questions are kept apart: is the unit stereogenic (graph analysis), does the input specify it, and has its
descriptor been checked against the IUPAC definition? Stability is separate again: a stereogenic unit is not
necessarily an isolable stereoisomer.

### Tetrahedral centres and double bonds

R/S and E/Z come from RDKit's new CIP labeller. RDKit's list of *potential* stereocentres is generous (it lists
adamantane's bridgeheads); MolTalk counts an unspecified centre only if flipping it changes the molecule, and in
cages only if the alternative can actually be built in 3D. Substituted cubanes therefore report no stereocentres,
while norbornan-2-ol and camphor keep theirs.

Wedges belong to a particular drawing: a wedge starts at the stereocentre and points toward the viewer. Wedge/dash is
not a synonym for R/S; in (1R,3S)-1-fluoro-3-methylcyclohexane both centres carry wedges.

### Atropisomers (hindered axes)

- **Detection** (RDKit does not list unspecified axes): a non-ring single bond between two sp² ring atoms, each end
  with two ring neighbours that rank differently, and at least 3 of the 4 ortho positions substituted or fused (all
  4 when an axis atom is in a five-membered ring). BINAP, BINOL, gossypol and tetra-ortho-substituted biaryls
  qualify; biphenyl, 2,2′-dimethylbiphenyl and diphenic acid do not. This is a rule of thumb, not a barrier calculation.
- **C–N axes:** an aryl ring on a tertiary amide nitrogen with two different ortho substituents (metolachlor: four
  stereoisomers). Ordinary amides are never axes.
- **Representation:** plain SMILES cannot state a twist. It travels as CXSMILES with 2D coordinates and an axis wedge,
  which `canonical_smiles` returns and which MOL, SDF and CDXML preserve (round-trip tested).
- **Verified descriptors only:** RDKit names the axis P or M; a 3D model is measured independently by the helicity
  rule and the CIP axial rule. Only when all three agree (Ra = M, Sa = P) is the descriptor shown.
- **Drawing:** one ring bond at the axis is wedged or hashed and labelled (Ra)/(Sa).

### Allenes and cumulenes

A chain with an even number of cumulated double bonds whose two ends each carry two different groups (1,3-dichloroallene).
Odd cumulenes are planar and get E/Z. The descriptor is M/P (IUPAC's preferred form) with Ra/Sa as the alternative,
both measured on 3D coordinates and required to agree as IUPAC states ((1M) = (1Ra)). Drawn with one end's
substituents in the paper and the other end's wedged and hashed.

### Spiro compounds

- **Xabcd** (four different ring neighbours): an ordinary R/S centre.
- **Xaabb** (2,6-disubstituted spiro[3.3]heptanes): IUPAC uses R/S at the spiro and ring atoms, with M/P as the
  alternative; the tests reproduce IUPAC's example (2R,4S,6R) = (2P). RDKit's newer stereo perception is used for
  these molecules only.
- **Xabab** (spiro[4.4]nonane-1,6-dione): the configuration is represented, enumerated, drawn and built in 3D, but no
  descriptor is given, because it needs a CIP extension (IUPAC P-93.5.3.2) that neither RDKit nor MolTalk implements.
- Achiral spiranes (spiropentane, spiropentadiene, spiro[3.3]heptane) report nothing.

### Helicenes

Five or more six-membered rings ortho-fused in a chain that always turns the same way. P or M from the handedness of
the helix through the ring centres (right-handed = P). [6]helicene and longer are configurationally stable;
[5]helicene racemises slowly at room temperature; [4]helicene, picene and pentacene are not reported. Drawn as a view
of the helix with a (P)/(M) label, never with wedges.

### Representation and export of advanced units

SMILES cannot state the configuration of an allene, a helicene or an Xabab spiro atom. MolTalk does not invent a
SMILES extension: the configuration travels as standard CXSMILES with 3D coordinates (or a 3D MOL/SDF) and is measured
from them. 2D MOL and CDXML cannot carry it, and the export says so.

### Stereo labels in the viewer

**Stereo: Specified | All | Off** (default Off). *Specified* labels what the input fixes; *All* also marks open
elements, "(?)" in the flat drawing and "(arb. R)" in 3D (the configuration this one conformer happens to have).
MolTalk does not write R\*/S\*, which in IUPAC usage means relative configuration.

## Lone pairs, charges and radicals

- **Count:** valence electrons, minus formal charge, minus one per bond (implicit hydrogens included; aromatic rings
  in a Kekulé form), minus a pair for each dative bond donated, minus unpaired electrons, halved. Heteroatoms always;
  carbon only when charged or a radical; metals never. An odd or negative count is left out rather than guessed.
- **Flat drawing (Lewis convention):** pairs share the open gaps between bonds, all at one radius from the element
  symbol. A formal charge is written outside the pairs, conventionally at the upper right.
- **3D:** directions come from the geometry actually built: opposite the bonds at a pyramidal centre, in the p orbital
  at a flat one, tetrahedral or trigonal as the shape requires. No lone pair or radical may lie within 80° of a bond;
  anything else is re-spread by repulsion. Charges stay upright typography outside the atom's footprint.
- **Radicals:** drawn as single dots, once (RDKit's own dot is hidden when the lone-pair overlay is on).
- **Charged carbon labels** keep their hydrogens ("CH₂⁻"), which RDKit omits by default.
- **Expanded octets** (SF₄, XeF₄, PCl₅, SF₆, I₃⁻) are built in their VSEPR shapes.

## 2D depiction

- **Quality scoring:** every layout is scored for bond crossings, overlapping atoms and stretched bonds.
- **Crowded molecules** (BINAP, Xantphos, rubrene): candidate layouts are repaired by rigid moves and scored on bond
  length, crowding, crossings, bond angles and ring regularity. A candidate is accepted only if it adds no crossings
  or overlaps and its wedges read back as the input's stereochemistry. On 500 library compounds, layouts with a bond
  over 1.5× the median went from 6 to 3, with no new crossings or overlaps.
- **Cages** (cubane, adamantane, dodecahedrane, C₆₀) are drawn as a view of their 3D shape, as textbooks do.
- **Porphyrin-type macrocycles** use RDKit's ring templates or are pinned to the porphine square, with the metal at
  the centre.
- **Hydrogens:** the heavy-atom skeleton is laid out first, so drawings with and without hydrogens share it.

## 3D models

- **Geometry:** ETKDGv3 followed by MMFF94 (or UFF). Several conformers are generated and the one best matching the
  flat drawing is kept and aligned to it by a proper rotation, so wedges point toward the viewer when the drawing
  lifts. It is labelled as one calculated conformer, not a unique structure.
- **Restraints where force fields fail:**
  - allenes are held twisted (MMFF and UFF flatten them);
  - localised carbanions are pyramidal;
  - σ radicals (vinyl) and carbenes are bent. SMILES does not give a carbene's spin state: a heteroatom neighbour
    implies a singlet (about 105°), otherwise a triplet (about 136°).
- **Special shapes:**
  - spiro centres joining small rings are set to D₂d;
  - helicenes are embedded with ETDG;
  - porphyrins are embedded as the flat free base, with the metal-binding cavity fitted to the metal.
- **Large molecules:**
  - long chains are straightened to the extended zigzag;
  - long tails are swung away from the ring system;
  - strained candidates are dropped.
- **Stereochemistry is checked in 3D:** axis twists, allene and helix configurations and spiro systems are verified
  on the model and corrected or the conformer rejected.

## Metal complexes

PubChem, and so the library, stores most complexes as a bare metal atom beside loose ligands and counter-ions
(Wilkinson's catalyst is `[Cl-].[Rh].PPh3.PPh3.PPh3`). For the drawing and 3D model only, MolTalk assembles them; the
analysis, formula and charges always use the stored record.

- **Porphyrin-type chelates** (heme, chlorophyll, phthalocyanines, B12): the metal is bonded to the four ring
  nitrogens at the centre of the square ring; bonds are covalent when charges balance exactly, else dative.
- **General mononuclear complexes:**
  - **Pieces:** each one is classified as a counter-ion (K⁺, PF₆⁻, BF₄⁻ …), a ligand, or solvent.
  - **Donors and chelates:** donors are P, CO and isocyanide C, pyridine/amine/imine N, halides, O⁻, S⁻ and hydride.
    Chelates must close 5- or 6-membered rings (bipyridine, salen, oxalate).
  - **π ligands:** η⁵-Cp and η²-alkenes.
  - **Binding order:** chelates and Cp first, then neutral donors, then anions, while the electron count stays within
    18 (16 for square-planar d8 Ni/Pd/Pt/Au).
  - **Oxidation state:** the ionic model, correcting PubChem's habit of writing Rh(I) as neutral Rh beside Cl⁻ and
    coordinated amines as "[NH⁻]".
- **Bonds:**
  - neutral donors are dative arrows;
  - anions bond covalently as far as the metal's charge allows;
  - bound CO is drawn M–C≡O;
  - π ligands get a bond to the ring or alkene centre.
- **Geometry** from positions and d-electron count:
  - 2 positions: linear;
  - 3 positions: trigonal planar;
  - 4 positions: square planar (d8/d9) or tetrahedral;
  - 5 positions: trigonal bipyramidal (d8/d10) or square pyramidal;
  - 6 positions: octahedral.

  Complexes already stored bonded (Grubbs' catalyst) get the same geometry.
- **2D:** the metal is the hub, ligands sit in evenly spaced slots, chelates take adjacent slots, and crowded ligands
  are drawn on longer metal–ligand bonds, as in textbooks.
- **3D:**
  - each ligand keeps its embedded shape and is moved rigidly into the ideal geometry, its bulk pointing outward;
  - chelates are fitted rigidly, multidentate ones embedded bonded to the metal;
  - monodentate ligands are turned to avoid clashes;
  - counter-ions are set clear of the complex.

Tested on Wilkinson's, Vaska's, Crabtree's and Jacobsen's catalysts, Pd(PPh₃)₄, ferrocene, zirconocene dichloride,
Zeise's salt, the platinum drugs, metal carbonyls, Ru(bpy)₃²⁺, Grubbs I and others. See
[LIMITATIONS.md](LIMITATIONS.md#metal-complexes) for what is out of scope.

## IUPAC names and numbering

Structure → name never guesses:

1. **Lookup:** the molecule's full InChIKey (an exact structure match including stereochemistry) is looked up in the
   library, then PubChem. Names are PubChem's (OpenEye Lexichem). A structure in neither gets no name.
2. **Check:** OPSIN rebuilds each candidate name; a name is used only if its structure equals the molecule's.
3. **Numbering:** parent-structure locants come from OPSIN's numbering. They are used only when the parent is
   unambiguous and every R/S descriptor in the name matches MolTalk's own CIP label at that locant.
4. **Otherwise:** atom indices.

The viewer's **Numbers: IUPAC | Index | Off** control (default Off) switches labels.
