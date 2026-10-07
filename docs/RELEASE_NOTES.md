# Release notes

## 1.0.0 (October 2026)

The first public release, submitted for ChatGPT plugin review. (Release candidate 1.0.0-rc.1 had the same features.)

**Policies.** [Privacy policy](PRIVACY.md), [terms of service](TERMS.md) and [support](SUPPORT.md) pages; the plugin
package carries these URLs, five positive and three negative review cases, the demo recording and these release notes.

**Deployment.** Public MCP endpoint on Google Cloud Run
(`https://moltalk-411294000488.us-central1.run.app/mcp`), stateless and scaling to zero. Production revision:
recorded in the release tag message. The earlier Secure MCP Tunnel setup is kept as a legacy option
([DEPLOYMENT.md](DEPLOYMENT.md#legacy-openai-secure-mcp-tunnel)).

**Tools.** `draw_named_molecule`, `draw_molecule`, `resolve_name`, `analyze_molecule`, `enumerate_stereoisomers`,
`find_substructure`, `export_structure`, and the viewer-only `conformer_3d`. The server reports its own version in
`serverInfo`.

**Chemistry.**
- **Names:** name resolution from a bundled library of 22,038 compounds, OPSIN and optional PubChem; it fails closed
  on unknown, ambiguous and stereo-inconsistent names, and labels are checked.
- **Stereochemistry:** R/S and E/Z; atropisomers (biaryl and C–N axes); allenes; spiro compounds; helicenes. All are
  reported as unified stereogenic units, with descriptors verified against IUPAC definitions.
- **Electrons:** lone pairs, radicals and formal charges in 2D and 3D, placed from the real geometry.
- **Metal complexes:** porphyrin-type chelates, and general mononuclear complexes assembled by electron counting and
  built in ideal geometry. η⁵-Cp and η²-alkenes are supported.
- **Depiction:** crowded molecules get repaired layouts; cages and helicenes are drawn as views of their 3D shape.

**Viewer.** Inline drawing with drag-to-rotate 3D, zoom, hydrogens, lone pairs, atom numbering and stereo labels; SVG
and PNG export; full-window viewer mode; phone support.

**Exports.** CDXML, MOL, SDF, PDB, XYZ and SMILES/CXSMILES, with warnings when a format cannot carry the
stereochemistry. CDXML no longer contains RDKit's empty `BondLength=""` attribute.

**Documentation.** A new public README, documentation under `docs/`, examples, and animated demos generated from the
real viewer (`scripts/generate_readme_demos.py`).
