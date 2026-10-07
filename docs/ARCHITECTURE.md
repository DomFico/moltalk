# Architecture

MolTalk is a Python MCP server that gives AI assistants chemistry tools: name resolution, validation, 2D drawing,
interactive 3D, analysis and file export. All chemistry is computed by [RDKit](https://www.rdkit.org/) and
[OPSIN](https://github.com/dan2097/opsin), not by the language model.

```text
 ChatGPT / Claude / any MCP client
        │  MCP (streamable HTTP, stateless JSON responses)
        ▼
 Cloud Run: moltalk container ─────────────────────────────────────────────┐
 │  public.py   gateway: Host/Origin checks, per-user rate limits, logging │
 │  server.py   FastMCP tools, server instructions, UI resource            │
 │     │                                                                   │
 │     ▼  worker processes (limits.py: timeout, memory cap, bounded queue) │
 │  naming / library ── name → structure (library → OPSIN → PubChem)       │
 │  chemistry      ──── parse, analyse, stereo summary, draw (SVG)         │
 │  depiction      ──── 2D layout, quality scoring and repair              │
 │  conformer      ──── 3D model aligned to the drawing                    │
 │  stereo / atropisomer / stereounits ── stereogenic units, verification  │
 │  coordination / complexes ── porphyrins and metal complexes            │
 │  export         ──── CDXML, MOL, SDF, PDB, XYZ, SMILES                  │
 └─────────────────────────────────────────────────────────────────────────┘
        │  tool result: data for the model + SVG/3D data for the viewer
        ▼
 widget/molecule.html  (MCP Apps UI in the host's sandboxed iframe)
```

## Components

| Module | Role |
|---|---|
| `moltalk/server.py` | The MCP server (FastMCP): tools, server instructions for the model, the `ui://widget/…` resource, label checking, name-stereo verification. |
| `moltalk/public.py` | The HTTP gateway used on Cloud Run: Host/Origin validation, rate limits, `/health`, icon routes, request logging (method and tool name, never arguments). |
| `moltalk/limits.py` | Runs RDKit work in worker processes with a wall-clock timeout, a memory cap and a bounded queue, so one expensive request cannot stall the server. |
| `moltalk/library.py`, `moltalk/data/compounds.json.gz` | The bundled compound library (22,038 compounds, about 200,000 names), with tiered name lookup. |
| `moltalk/naming.py` | Name → structure (library, OPSIN, optional PubChem) and structure → IUPAC name and locants. |
| `moltalk/chemistry.py` | Parsing, analysis, the unified stereo summary, SVG drawing, stereoisomer enumeration. |
| `moltalk/depiction.py` | 2D layout: RDKit/CoordGen candidates, scoring and repair, porphyrin templates, cage projections, complex and π-ligand layouts. |
| `moltalk/conformer.py` | The 3D model: ETKDG + MMFF/UFF, alignment to the flat drawing, VSEPR and coordination geometry, lone-pair directions. |
| `moltalk/stereo.py`, `atropisomer.py`, `stereounits.py` | Stereochemistry beyond RDKit's defaults: real vs apparent stereocentres, hindered axes, allenes, spiro units, helicenes. |
| `moltalk/coordination.py`, `complexes.py` | Metal complexes stored as loose pieces: porphyrin-type chelates and general mononuclear complexes. |
| `moltalk/export.py` | Structure files written by RDKit. |
| `moltalk/widget/molecule.html` | The viewer: flat drawing, drag-to-rotate 3D, labels, lone pairs, image and file export. No network access. |

## Request flow

1. **The model calls a tool**, for example `draw_named_molecule("caffeine")`.
2. **The name is resolved** to a structure: the bundled library first, then OPSIN for systematic names, then PubChem only when the model explicitly allows network access. Ambiguous, unknown or stereo-inconsistent names fail with an explanation; no structure is guessed.
3. **The structure is parsed and validated** by RDKit. Invalid valences or syntax are errors.
4. **Analysis**: formula, identifiers, properties, CIP labels, a stereo summary that separates specified from unspecified elements, functional-group motifs, and (when the structure is known) the IUPAC name and locants.
5. **Drawing**: a 2D layout is computed and quality-checked, then drawn by RDKit as SVG with wedges and annotations.
6. **The result** goes back in two parts: chemistry data in `structuredContent` (what the model reads) and the SVG in `_meta` (only the viewer sees it).
7. **The viewer** shows the drawing. When the user drags it, the viewer calls the hidden `conformer_3d` tool once and rotates that model locally.

## Name resolution pipeline

```text
name ──► bundled library (offline, tiered: curated > titles > aliases > synonyms)
            │ not found
            ▼
         OPSIN (offline, systematic IUPAC names incl. R/S, E/Z)
            │ not readable
            ▼
         PubChem (only with allow_network=true)
            │ not found
            ▼
         error: "No structure was assumed"
```

Any stereodescriptor in the name must be present in the structure found, with matching CIP labels, or the name is
rejected. See [CHEMISTRY.md](CHEMISTRY.md#names-and-validation).

## 2D pipeline

RDKit's layout is used when it is clean. Otherwise candidate layouts (RDKit, sampled RDKit, CoordGen, a hub layout
for star-shaped and metal-centred molecules) are repaired by a small search over rigid moves and scored; the winner
must not add crossings or overlaps, and its wedges must read back as the input's stereochemistry. Special cases:
porphyrin templates, cage projections (cubane, C60), helicene projections, and π-ligand layouts (ferrocene).

## 3D pipeline

ETKDG embeds several conformers, a force field (MMFF94, or UFF) cleans them up, and the conformer that best matches
the flat drawing is kept and aligned to it by a proper rotation (never a mirror, so R/S is preserved). Corrections
on top: stereochemistry checks for axes, allenes, helicenes and spiro systems; VSEPR restraints for radicals,
carbanions and carbenes; ideal geometry for metal complexes; straightened chains; flattened porphyrins.

## Exports

`export_structure` writes CDXML, MOL, SDF, PDB, XYZ or SMILES with RDKit, from the stored structure (not drawing-only
additions). 2D files use the drawing's layout; 3D files use the viewer's conformer. The viewer itself exports SVG and
PNG of exactly what is on screen. See [MCP_AND_UI.md](MCP_AND_UI.md#exports).

## Security and limits

- RDKit runs in two worker processes with a 25 s timeout on Cloud Run (20 s default), a 2 GB memory cap and at most
  16 queued requests.
- Input limits: 16,384 characters of SMILES or CXSMILES, 256 atoms, 1,024 SMARTS characters, 100 substructure
  matches, 64 enumerated isomers (16 drawn), drawings of 200–1,600 px, names of up to 256 characters.
- HTTP mode always validates Host and Origin headers (DNS-rebinding protection); request bodies are limited to 64 KB.
- The public deployment rate-limits each user (60 requests/minute) and the service as a whole (600/minute).
- The viewer has an empty content-security allowlist: it makes no network requests.
- Privacy: names are sent to PubChem only when the model asks for network lookup; structures not in the library have
  their InChIKey looked up on PubChem for a name (disable with `MOLTALK_OFFLINE=1`).
