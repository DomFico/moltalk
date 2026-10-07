# MCP tools and the viewer

## Tools

| Tool | Shows UI | Returns to the model |
|---|---|---|
| `draw_named_molecule(name, numbering, hydrogens, allow_network)` | yes | Resolves a compound name and draws it in one call; the preferred tool for "draw X". |
| `draw_molecule(smiles, label?, width, height, numbering, hydrogens, include_svg)` | yes | The full analysis, `depicted_stereo_bonds` (each wedge/dash with a plain-language explanation), layout quality, and a `label_check` when a name is given. |
| `resolve_name(name, allow_network)` | no | The structure for a name: source, PubChem CID, title, IUPAC name, `stereo_summary`, warnings. Ambiguous or unknown names are errors. |
| `analyze_molecule(smiles)` | no | Identifiers, properties, CIP labels, `stereo_units`, `stereo_summary`, functional-group motifs. |
| `enumerate_stereoisomers(smiles, limit)` | yes (grid) | Stereoisomers with descriptors, meso flags and enantiomer pairs, including atropisomers, allenes, helicenes and spiro units (as CXSMILES where SMILES cannot state them). |
| `find_substructure(smiles, smarts)` | no | Chirality-aware matches (up to 100). |
| `export_structure(smiles, format, coordinates, name, hydrogens)` | yes (download card) | A structure file: CDXML (default), MOL, SDF, PDB, XYZ or SMILES, in 2D or 3D, with warnings when a format cannot carry the stereochemistry. |
| `conformer_3d(smiles, hydrogens)` | viewer only | The 3D model the viewer rotates. Hidden from the model. |

**Server instructions** tell the model to resolve names through MolTalk, report sources and stereo status, explain
wedges only from `depicted_stereo_bonds`, reuse `canonical_smiles` in follow-ups, draw a molecule once per reply,
and never substitute a structure after an error.

**Result format.** Chemistry data is in `structuredContent` (what the model reads). The SVG is in the result's
`_meta`, which hosts pass to the viewer but not to the model. Pass `include_svg=true` for clients without UI support.

## The viewer

`moltalk/widget/molecule.html` is registered as an MCP Apps UI resource (`ui://widget/molecule-v32.html`,
`text/html;profile=mcp-app`) with an empty CSP allowlist; it makes no network requests. It uses the MCP Apps bridge
(`ui/initialize`, `ui/notifications/tool-result`, `tools/call`, `ui/notifications/size-changed`,
`ui/request-display-mode`, `ui/download-file`) and falls back to `window.openai`. Bump the version in the URI when
the widget changes, because hosts cache templates.

### Controls

- **Rotate:** drag the drawing (one finger on a phone). Over the first ~90 px it lifts off the page into 3D and keeps
  rotating as a chemical drawing: implicit carbons, element labels with their hydrogens, half-coloured bonds, and
  gaps where a bond passes in front of another. Arrow keys rotate; **2D** or Escape returns to the exact flat drawing.
- **Zoom and pan:** Ctrl/Cmd+scroll, trackpad or two-finger pinch, + and − keys; Shift+drag or a two-finger drag
  pans. Plain scrolling still scrolls the chat.
- **Show hydrogens**, **Show lone pairs**, **Numbers: IUPAC | Index | Off**, **Stereo: Specified | All | Off**.
  Changing them keeps the current view (rotation, zoom, pan).
- **Enter viewer** asks the host for its full-window view; the widget keeps clear of the host's own bars.
- **Export: SVG / PNG.**

### Fitting the host

The widget fits the host's inline height limit, so the host never adds a scroll bar. Inline, the drawing is at most
420 px tall; some hosts cap the frame without reporting a limit, and the widget detects that and shrinks the drawing.
The rotated model shrinks gradually as it lifts so it stays in the frame in every orientation.

### Repeated draws

ChatGPT sometimes re-sends the first `draw_molecule` of a chat with identical arguments. When a session's draw repeats
its previous draw within 20 s, the second viewer shows "Same drawing as above"; the model still gets the full result.

## Exports

- **Images:** **SVG** and **PNG** save exactly what is on screen (rotation, zoom, hydrogens, labels, lone pairs),
  cropped, with a transparent background. Gaps and label backings become real cut-outs, so nothing shows white on a
  coloured slide. PNG is 3× the on-screen size.
- **Structure files:** ask for a file ("give me a ChemDraw file of that") and the model calls `export_structure`. The
  result is a card with a **Download** button.
  - Formats: CDXML (ChemDraw opens it directly), MOL, SDF, PDB, XYZ, SMILES. Binary `.cdx` is not offered; CDXML
    replaces it.
  - 2D files use the drawing's layout with stereo as wedges. 3D files use the viewer's conformer with explicit
    hydrogens; when the input leaves stereo open, the file and the result say which configuration was chosen.
  - The structure as given is exported (atom order, charges, stereo), not drawing-only additions.
- **Downloads:** the standard MCP Apps `ui/download-file` when the host supports it; otherwise a browser download plus
  a panel to copy the file text or save the image, which works even where the sandbox blocks downloads. The
  structure never goes into a URL.
