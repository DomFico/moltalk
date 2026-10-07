# Testing

```bash
.venv/bin/pip install -e '.[dev,ui-test]'
.venv/bin/pytest -q                              # everything offline, including the browser tests
MOLTALK_NETWORK_TESTS=1 .venv/bin/pytest -q      # also the live PubChem tests
```

The suite runs offline by default (`MOLTALK_OFFLINE=1` is set for the run); two tests that need live PubChem are
skipped unless `MOLTALK_NETWORK_TESTS=1`. The full suite takes about three minutes on a desktop.

## What is covered

| Area | Files | Examples |
|---|---|---|
| Names and validation | `test_library.py`, `test_naming.py` | Library lookups, OPSIN names with stereo, ambiguous and unknown names failing closed, D/L families, curated cofactors, label checks, axial names rejected, IUPAC numbering |
| Chemistry and analysis | `test_chemistry.py`, `test_lone_pairs.py` | Parsing errors, stereo summaries, real vs apparent stereocentres, cages, electron counting |
| Stereochemistry | `test_atropisomer.py`, `test_stereounits.py` | Atropisomers (biaryl and C–N), allenes, spiro compounds, helicenes; positive and negative controls; descriptors against IUPAC reference examples; round trips through CXSMILES, MOL, SDF, CDXML; ordinary R/S and E/Z unchanged |
| 2D depiction | `test_depiction.py`, `test_stereounits.py` | Porphyrin templates, cages, crowded-layout panel (no crossings, longest bond ≤ 1.5×), wedge faithfulness |
| 3D models | `test_conformer.py`, `test_complexes.py` | Alignment with the drawing, configurations preserved, spiro D₂d, allenes twisted, lone pairs never on a bond, VSEPR shapes for radicals and carbenes |
| Metal complexes | `test_coordination.py`, `test_complexes.py` | Porphyrins and B12; assembly of complexes stored as pieces, electron counts, counter-ions, ideal geometries in 3D |
| MCP behaviour | `test_mcp.py`, `test_public.py` | Tool schemas and results, UI resource metadata, repeated draws, exports, the HTTP gateway (Host/Origin checks, rate limits) |
| Viewer (real browser) | `test_widget_*.py` | Rendering through an MCP Apps host harness, 3D rotation, zoom, hydrogens, touch, phone-sized viewer, stereo labels, SVG/PNG export, radical dots |

The browser tests drive the real widget in Google Chrome via Playwright (set `MOLTALK_UI_BROWSER_CHANNEL=chromium` to
use Playwright's Chromium).

## Regression checks used during development

- **Stereo regression:** canonical SMILES, R/S, E/Z and stereo summaries on a random sample of 400 library compounds
  are compared against the previous release.
- **Layout regression:** crossings, overlaps and bond-length ratios on a random sample of 500 library compounds,
  before and after a depiction change.

## Release smoke test

Run against the production URL from a real host (ChatGPT web and mobile, Claude):

1. A known trivial name: "Draw caffeine and identify its functional groups."
2. A systematic stereo name: "Draw (2R)-2-chloro-4-methylhexane."
3. A custom structure (SMILES).
4. An atropisomer: "Show the stereoisomers of BINAP", then draw one.
5. Lone pairs and charges: "Draw nitromethane with lone pairs."
6. A radical or ylide: "Draw the Wittig ylide Ph₃P=CH₂ as its ylide form."
7. A cage or spiro compound: "Draw cubane."
8. A coordination complex: "Draw ferrocene."
9. An export: "Give me a ChemDraw file of that."
10. Failure: "Draw florbanex quintophane" and "Draw glucose" (unknown and ambiguous names must fail closed).
