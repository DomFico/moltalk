# MolTalk — draw, rotate and name organic molecules

A ChatGPT plugin (RDKit + MCP). Python package `moltalk`, command `moltalk`, background service `moltalk-tunnel`, settings in `~/.config/moltalk`.

Ask ChatGPT ordinary chemistry questions. ChatGPT calls RDKit tools on your computer, and the molecule drawings appear inline in the conversation.

```
ChatGPT (developer-mode personal plugin)
      │  MCP over OpenAI's Secure MCP Tunnel (outbound-only from your computer)
      ▼
tunnel-client  (runs on this computer, systemd user service)
      │  stdio. No HTTP port is opened, and nothing listens on the network except a loopback admin page.
      ▼
moltalk  (.venv/bin/moltalk → RDKit worker processes with time/memory limits)
```

## Status: verified vs. still to do

| Item | Status |
|---|---|
| Install in an isolated venv (`.venv`, Python 3.12, RDKit 2026.03.6, mcp 1.30.0) | **Verified** |
| Test suite: 58 tests (56 offline by default; the PubChem tests need `MOLTALK_NETWORK_TESTS=1`), including live PubChem and real-browser tests of the drawing, its 3D rotation, zoom, hydrogens and a phone touch screen | **Verified** (`pytest`) |
| Tools and widget through official MCP clients (Python SDK over stdio/HTTP, MCP Inspector CLI) | **Verified** |
| UI resource metadata (`ui://widget/molecule-v15.html`, `text/html;profile=mcp-app`, `_meta.ui.resourceUri` + `openai/outputTemplate`) | **Verified** |
| Widget rendering in a sandboxed iframe through the MCP Apps bridge, plus a follow-up `tools/call` from the widget | **Verified** in Chrome with a host harness that imitates ChatGPT (`tests/test_widget_ui.py`, `tests/test_widget_3d.py`) |
| tunnel-client → stdio server path | **Verified** with tunnel-client's local control plane (`scripts/local-tunnel-test.sh`) |
| Live PubChem resolution | **Verified** (aspirin, lactic acid, L-lactic acid, glucose, unknown name) |
| Docker image (build, read-only fs, non-root, Host/Origin rejection) | **Verified** |
| Creating the OpenAI tunnel and runtime key | **Needs you** (account action) |
| Creating and installing the plugin in ChatGPT; tool calls and inline rendering *inside ChatGPT* | **Needs you.** These were not verified, because they require your logged-in ChatGPT account |

## 1. Install and test (done once)

```bash
git clone https://github.com/DomFico/moltalk.git && cd moltalk
scripts/fetch-tools.sh                # OPSIN (IUPAC numbering) and tunnel-client, checksum-verified
python3 -m venv .venv                 # Conda's base env can stay active; the venv is isolated
.venv/bin/pip install -e '.[dev]'
.venv/bin/pytest -q                   # add MOLTALK_NETWORK_TESTS=1 to include live PubChem
# Optional browser test of the widget (uses Google Chrome):
.venv/bin/pip install -e '.[dev,ui-test]' && .venv/bin/pytest -q tests/test_widget_ui.py tests/test_widget_3d.py
```

`scripts/fetch-tools.sh` downloads `tunnel-client` v0.0.15 (linux-amd64, SHA-256 verified against the release's `SHA256SUMS.txt`) into `tools/tunnel-client/`, and OPSIN 2.9.0 (SHA-256 verified) into `vendor/`. IUPAC numbering also needs Java 11 or later on the PATH.

## 2. Connect to ChatGPT

### What you need, and cost

- **A ChatGPT Plus, Pro, Business, Enterprise or Edu account.** Developer mode is not available on Free, and plugins are created on the web.
- **An OpenAI Platform organization** (platform.openai.com) signed in with the same account. Creating one is free. OpenAI's tunnel docs do not list a price for Secure MCP Tunnel, and the tunnel does not call any models, so no token charges are expected. Check your Platform billing page if in doubt.
- **This computer must be on, awake and running the tunnel** whenever you use the plugin. If it is off, ChatGPT tool calls fail, and nothing else breaks.

> Account caveat: OpenAI documents tunnels as associated with "Platform organizations or ChatGPT workspaces". It does not say explicitly whether a *personal* Plus/Pro account can select a tunnel. If step 2.5 shows no tunnels, see [If the Tunnel option is unavailable](#if-the-tunnel-option-is-unavailable).

### 2.1 Turn on developer mode (ChatGPT, you)

ChatGPT → **Settings → Security and login → Developer mode: On** (chatgpt.com/settings/security).

### 2.2 Create a tunnel (Platform, you)

1. Open https://platform.openai.com/settings/organization/tunnels and create a tunnel, named for example `moltalk`.
2. Copy its ID (`tunnel_` + 32 hex characters). The ID is not a secret.

You need the *Tunnels Read + Manage* permission to create the tunnel and *Read + Use* to run and select it. Organization owners have both.

### 2.3 Create a runtime API key and store it locally (you)

1. Create a key at https://platform.openai.com/settings/organization/api-keys. If you can restrict its permissions, grant only Tunnels Read + Use.
2. In a terminal, run:
   ```bash
   scripts/store-tunnel-key.sh
   ```
   Paste the key at the hidden prompt. It is saved to `~/.config/moltalk/tunnel-runtime-key` with mode 600. **Never paste the key into ChatGPT, a chat with an assistant, or a command line.**

### 2.4 Configure and start the tunnel

```bash
scripts/configure-tunnel.sh tunnel_xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx   # writes the profile, runs tunnel-client doctor
scripts/start.sh                                                     # installs/starts the systemd user service
scripts/status.sh                                                    # healthz/readyz should both be OK
```

### 2.5 Create the plugin in ChatGPT (you)

1. Go to https://chatgpt.com/plugins and select **+** (shown as "Create MCP App" or similar once developer mode is on).
2. Fill in the form:
   - **Name:** `MolTalk`
   - **Description:** `Draw, rotate and name organic molecules: stereochemistry, IUPAC numbering and lone pairs (RDKit).`
   - **Connection:** **Tunnel**. Select your tunnel or enter its `tunnel_…` ID.
   - **Authentication:** **No authentication.** The server has no network listener, and only your tunnel can reach it.
3. Select **Create**. You should see 6 tools: `analyze_molecule`, `draw_molecule`, `enumerate_stereoisomers`, `find_substructure`, `resolve_name`, and `conformer_3d`. The last is viewer-only: it is marked hidden from the model and used only to rotate the drawing.

### 2.6 Use it

In a **new chat**, open **+ → Developer mode**, enable **MolTalk**, then ask your question. Naming the plugin helps at first, for example *"Using MolTalk, draw aspirin and identify its functional groups."*

Acceptance prompts and what to expect:

| Prompt | Expected |
|---|---|
| Draw aspirin and identify its functional groups. | `resolve_name` (offline dictionary) then `draw_molecule`. Inline drawing; groups: carboxylic acid, ester, aromatic ring |
| Draw (1R,3S)-1-fluoro-3-methylcyclohexane and explain the wedges/dashes in this drawing. | Drawing of `F[C@H]1C[C@@H](C)CCC1`. Atom 1 = R, atom 3 = S. **Both** bonds are wedges (1→F, 3→CH₃), toward the viewer, so wedge ≠ R/S |
| Show the stereoisomers of lactic acid when stereochemistry is unspecified. | Grid of 2: atom 1 R and atom 1 S, enantiomers of each other. If looked up via PubChem: CID 612, stereo reported as unspecified |
| Draw C1CC(C)(C)(C)C1 | Error card: "Chemically invalid structure (Explicit valence for atom # 2 C, 5 …)". No molecule is drawn or substituted |
| (follow-up) What is its molecular weight? / Is it chiral? | Reuses the earlier `canonical_smiles` |

The widget's **Hide/Show atom indices** button redraws the structure through a follow-up tool call from inside the drawing. Dragging the drawing rotates it in 3D (see below).

### 2.7 Day-to-day commands

| Action | Command |
|---|---|
| Start (background, restarts on failure) | `scripts/start.sh` |
| Run in the foreground instead (Ctrl-C stops it) | `scripts/start.sh --foreground` |
| Stop | `scripts/stop.sh` |
| Status | `scripts/status.sh` (admin UI at http://127.0.0.1:8791/ui, loopback only) |
| Logs | `journalctl --user -u moltalk-tunnel -f` |
| Start automatically at login (optional) | `systemctl --user enable moltalk-tunnel` |
| Remove the service | `scripts/uninstall-service.sh`. To also remove the key and profile: `rm -r ~/.config/moltalk` |

After changing tool metadata or the widget: restart (`scripts/start.sh`), then in ChatGPT open the plugin under **Plugins** → **Refresh**, and start a new chat. If you change the widget substantially, bump the version in its URI (`molecule-v15` → `v16`), because ChatGPT caches templates.

### If the Tunnel option is unavailable

A public HTTPS endpoint is the documented alternative. ChatGPT accepts only **OAuth** or **no auth** for MCP servers, and static API keys are not supported. This project **does not** expose its unauthenticated HTTP server publicly, so going public needs one more piece of work: an OAuth 2.1 authorization layer, using either the MCP SDK's auth support or a hosted identity provider with Dynamic Client Registration. That layer would sit in front of the Docker image on a small HTTPS host (typically about US$5–10 per month), with `MOLTALK_ALLOWED_HOSTS` set to that domain. Decide on a host and identity provider before this is built.

## 3. Tools

| Tool | Shows UI | Returns to the model |
|---|---|---|
| `analyze_molecule(smiles)` | no | identifiers, properties, CIP R/S, E/Z, `stereo_summary`, functional-group motifs |
| `draw_molecule(smiles, label?, width, height, atom_indices, include_svg)` | **yes** | full analysis + `depicted_stereo_bonds` (each wedge/dash with a plain-language explanation) |
| `enumerate_stereoisomers(smiles, limit)` | **yes** (grid) | isomers with CIP labels, `achiral` (meso) flag, `enantiomer_index` |
| `find_substructure(smiles, smarts)` | no | chirality-aware matches (≤100) |
| `conformer_3d(smiles)` | viewer-only | one calculated 3D conformer aligned to the flat drawing, for rotating it (hidden from the model) |
| `resolve_name(name, allow_network)` | no | offline dictionary, or PubChem CID, title, IUPAC name, `stereo_summary` and warnings |

**UI component.** `moltalk/widget/molecule.html` is registered as `ui://widget/molecule-v15.html` (`text/html;profile=mcp-app`) with an empty CSP allowlist, so it makes no network requests. It uses the MCP Apps bridge (`ui/initialize`, `ui/notifications/tool-result`, `tools/call`, `ui/notifications/size-changed`) and falls back to `window.openai`. The SVG is sent only in the result's `_meta`, which ChatGPT passes to the widget but not to the model; the chemistry data stays in `structuredContent` for the model. The SVG is sanitized before insertion. Pass `include_svg=true` for clients without UI support.

**Rotate the drawing in 3D.** Grab the flat drawing and drag (one finger on a phone). Over the first ~90 px it lifts off the page into 3D and keeps rotating as a chemical drawing:
- **Drawing style:** flat colour, implicit carbons, element labels with their hydrogens (OH/HO, NH₂) in RDKit's colours, half-coloured bonds, ring double bonds on the inner side, and gaps where a bond passes in front of another.
- **Zoom and pan:** **−/+** (1.2× steps that glide rather than snap), Ctrl/Cmd+scroll or a trackpad pinch (continuous), two-finger pinch on touch screens; Shift+drag or a two-finger drag pans. Zoom works on the flat drawing too. Plain scrolling still scrolls the chat.
- **Phones:** touches on the drawing belong to the drawing, so a swipe rotates it instead of scrolling the page or switching chats. To scroll past it, swipe outside the drawing.
- **Enter viewer** asks ChatGPT for its full-window viewer. There the widget fits the window exactly with no scrolling: the drawing takes the free space, and the notes and lists are hidden. ChatGPT draws its own bars over the viewer (a title bar on phones, and the "Ask ChatGPT" box at the bottom everywhere). The widget keeps clear of them, using the host's safe-area insets or minimums measured from screenshots if those are larger: 104/128 px top/bottom on phones, 88 px at the bottom on desktop. Our controls sit at the top, under the header; on touch devices the clearances are 132/140 px. Returning inline, the widget re-reports its height. On iPhone, ChatGPT shows the app in the chat as a collapsed "Open App" card and adds its own pull-up sheet; both are host behaviour. **Known iOS host bug:** closing the fullscreen view with ChatGPT's own X can make the inline widget disappear (open OpenAI community report, unresolved as of 1 Oct 2026). The widget does not call `requestClose()`; it offers **Exit viewer**, which asks the host for `inline` through `requestDisplayMode` / `ui/request-display-mode`. It declares `availableDisplayModes: [inline, fullscreen]`, reports its height only when inline and never below 280 px, and re-reports it after every mode change.
- **Reset** (or Escape) restores zoom and settles the drawing back onto the exact RDKit picture. Arrow keys rotate.
- **Inline size:** the widget fits inside ChatGPT's inline height limit (`containerDimensions.maxHeight` or `window.openai.maxHeight`), so the host never adds a scroll bar. The rotated-view note, chips, wedge list and notes share one small box that scrolls on its own. Under a tight limit the drawing shrinks rather than the controls.
- **Show hydrogens / atom indices** redraw without changing the view (the rotated picture stays on screen while the redraw loads, and the other hydrogen setting's 3D data is fetched in advance): rotation, on-screen size, zoom and pan are kept. The conformer is chosen from heavy atoms only, so it is the same shape with or without hydrogens. It is aligned to whichever drawing is on screen: RDKit lays out the explicit-H drawing in its own orientation (LSD's is turned 180°), and `to_heavy_frame` converts the remembered rotation between the two, so a toggle never flips or moves the molecule.
- **Geometry:** RDKit ETKDGv3 followed by MMFF94 (or UFF); fullerene-like cages start from a spectral sphere.
- **Matching the flat drawing:** up to 8 conformers (32 if needed) are generated, and the one best matching the flat drawing is kept. Each is placed with a proper rotation only, never a mirror image, so R/S cannot change. Stereocentres are weighted heavily in the fit, so wedges point toward the viewer and dashes away when the lift begins. Each drawing makes one 3D request.
- **Caveats shown in the view:** it is labelled as one calculated conformer, and arbitrary configurations at unspecified stereocentres are marked "?".
- **Not yet done:** wedges/dashes are not redrawn during rotation.

**IUPAC name and numbering.** Under the title the widget shows **SMILES:** and **IUPAC:**, and **Numbers: IUPAC | Index | Off** chooses the atom labels. IUPAC is the default; the stereocentre chips and the model's text use locants too (e.g. "C6a: R"). The pipeline (`moltalk/naming.py`) never guesses:
1. **Name:** looked up in PubChem by the molecule's full InChIKey, an exact structure match including stereochemistry. PubChem's names are generated by OpenEye Lexichem. A structure not in PubChem gets no name ("not guessed").
2. **Checking the name:** OPSIN 2.9.0 (open source, `vendor/`, SHA-256 matches the GitHub release, needs Java) rebuilds each candidate name. A name is used only if its structure equals ours (connectivity and charges), which drops PubChem's zwitterion and salt duplicates. If different names still fit, none is used. A redundant `cis-`/`trans-` before R/S descriptors is removed for parsing only.
3. **Numbering:** OPSIN numbers every substituent in its own scheme. The parent is the scheme that contains all top-level locants of the name and whose bonds to other carbon schemes leave from those positions; ties are broken by which positions actually carry substituents. A remaining tie means no numbering. Only parent atoms are numbered. Symmetry-equivalent numberings are equally valid; one is used only if every R/S descriptor in the name matches our own CIP label at that locant.
4. **Fallback:** when any step fails, the widget says why, keeps the name if known, and uses atom indices.

Tested offline on 25 typical course molecules (chains, rings, aromatics, steroids, sugars, LSD). All parents were numbered correctly, except acetaminophen and ethyl acetate, whose names give no positional clue to the parent (refused), and heme, where PubChem's records disagree (refused). Privacy: each new molecule's InChIKey is sent to PubChem (cached); set `MOLTALK_OFFLINE=1` to turn this off. The Docker image has no Java, so it shows names without numbering.

**Lone pairs.** **Show lone pairs** draws Lewis-structure dots on the flat and rotated views, instantly, with no server call. The setting survives redraws.
- **Count:** for each atom, take its valence electrons, subtract its formal charge, subtract one electron per bond (hidden hydrogens included; aromatic rings counted in a Kekulé form), subtract a pair for each dative bond it donates, subtract any unpaired electrons, then halve. Unpaired electrons are drawn as single dots.
- **Which atoms:** heteroatoms always; carbon only when it is charged or a radical; metals never. An atom whose count is odd or negative is left out rather than guessed.
- **Flat placement:** pairs are shared across the open gaps between bonds (and the H label) in proportion to each gap's size, then spaced evenly inside it. So a C=O oxygen gets ±120° and a terminal F 90°, 180° and 270°. Dots sit just outside the label's exact box, measured from RDKit's own label paths.
- **Rotated placement:** directions come from the 3D conformer by VSEPR and hybridisation, and turn with the molecule.
  - sp³: tetrahedral (water 109°, F on carbon 109°).
  - sp²: 120° in the plane (a ketone O, pyridine N); any extra pair sits in the p orbital (pyrrole N, carboxylate O⁻).
  - sp: opposite the bond (nitrile N).
  - Anything else is spread by repulsion.
  - A pair seen end-on falls back to the flat spacing.
- **Expanded octets** (main-group atoms with 5–6 electron domains and single-atom ligands) are built in their VSEPR shape, because ETKDG would make them tetrahedral: SF₄ seesaw, XeF₄ square planar with trans pairs, XeF₂ and I₃⁻ linear, SF₆ octahedral, PCl₅ trigonal bipyramidal. RDKit rejects some valences outright (ClF₃, BrF₅) and says so.
- **Metal complexes:** ligand lone pairs are drawn, and dative donors have none. The arrangement around the metal comes from a generic force field (no ligand-field model, and cis/trans is not in the SMILES), so the 3D view says not to read coordination geometry from it.
- **Tested** on 18 neutral, anionic, cationic and radical species, for example hydroxide 3, hydronium 1, ammonium 0, the amide anion 2, the carbanion 1, the carbocation 0, carbon monoxide C⁻ 1 and O⁺ 1, nitro O⁻ 3, and the nitroxide O with 2 pairs plus 1 unpaired electron.

**Metal complexes.** Databases often store heme, chlorophyll or metal phthalocyanines as a charged ring plus a separate metal ion. For the drawing only, the metal is bonded to the four ring nitrogens and placed at the centre of the ring (`moltalk/coordination.py`).
- **When the charges balance exactly** (heme, chlorophyll), it draws the textbook 2 covalent and 2 dative (arrow) bonds with neutral atoms. Otherwise it draws 4 dative bonds and keeps the database charges.
- **Flat layout:** RDKit's ring templates give the textbook square porphyrin (pyrroles on the four sides), turned so the side chains hang at the bottom, as in the usual heme b drawing.
- **3D shape:** MMFF bends the charged porphyrin form out of plane (about 0.6–0.8 Å RMS), so the ring is embedded as the neutral free base and the flattest candidate is kept (about 0.13 Å for heme, in about 1 s). The metal is placed in the plane at the centre of the four nitrogens.
- **The analysis is unchanged:** formula, charges and the "disconnected fragments" warning describe the source structure. The drawing says what was added.
- **Charge-separated sources** (for example heme written as [Fe-2] bonded to two [N+]) are drawn the same way, with dative N→Fe bonds and neutral atoms.
- **Metals are never treated as stereocentres.** RDKit lists a four-coordinate metal as a possible tetrahedral centre, but square-planar Fe is not one.
- **Not handled:** other metal complexes, such as ferrocene's sandwich bonding or cisplatin's ammine ligands, are drawn as stored.

**Cages.** In bridged or cage-like ring systems (bridgeheads, or atoms in three or more rings, as in cubane), an unspecified element counts only if its R/S differs among the stereoisomers that can actually be built in 3D. Each candidate is tried with a small ETKDG budget and two seeds, mirror images are added automatically, and the check is capped at 256 candidates. Substituted cubanes therefore report no stereocentres, while norbornan-2-ol (3) and camphor (2) keep theirs. The rotated view shows R/S only at real stereocentres.

**Big molecules.** If 32 conformers time out, the best of the first 8 is kept, and the 32-conformer retry is skipped when the first round took over 1.5 s. Erythromycin gets its 3D model in about 4 s instead of failing.

**Stereocentre counting.** RDKit's list of *potential* stereocentres is generous; for example, it lists adamantane's bridgeheads. An unspecified centre now counts only if flipping it actually changes the molecule (`moltalk/stereo.py`), so adamantane reports no stereo elements, while cis/trans-1,4-dimethylcyclohexane still reports two.

**Drawing quality check.** Every layout is scored for bond crossings, stretched bonds and overlapping atoms (`moltalk/depiction.py`). RDKit's standard layout fails on polyhedral cages: for C₆₀ it produced 37 crossing bond pairs and bonds up to 16× normal length. When that happens, the largest ring system is redrawn as a **Schlegel diagram** using a Tutte embedding, which is crossing-free for polyhedral cages. The server tries each ring as the outer face, places substituents outward, and keeps the result only if it scores better. The `depiction` field reports the method and quality, and an imperfect layout carries a caveat in both the widget and the model's result. Tested on C₆₀, C₇₀, PCBM, dodecahedrane and cubane; ordinary molecules keep RDKit's layout.

**Server instructions** tell the model to resolve names through `resolve_name`, check CIP labels against stereo descriptors in a name, explain wedges only from `depicted_stereo_bonds`, reuse `canonical_smiles` in follow-ups, and never substitute a structure after an error.

## 4. Security and limits

- **Tunnel deployment:** the MCP server is a stdio child of tunnel-client and opens no ports. tunnel-client makes outbound HTTPS connections only, to `api.openai.com`, authenticated with the runtime key read from a mode-600 file. Its admin UI is bound to 127.0.0.1. The systemd unit sets `NoNewPrivileges`, `MemoryMax=3G` and `TasksMax=256`.
- **RDKit work** runs in two worker processes with a 20 s per-call timeout (the worker is killed when it is exceeded), a 2 GB address-space cap, at most 16 queued requests, and stdout redirected so a worker can never corrupt the stdio stream. Configure with `MOLTALK_TIMEOUT_S`, `MOLTALK_WORKERS`, `MOLTALK_WORKER_MEMORY_MB` and `MOLTALK_MAX_QUEUED`.
- **Input limits:** 4096 SMILES characters, 256 atoms, 1024 SMARTS characters, 100 substructure matches, 64 enumerated isomers (16 drawn), drawings 200–1600 px, names ≤256 characters.
- **HTTP mode** (local testing and Docker): Host/Origin validation (DNS-rebinding protection) is **always on**, including when binding 0.0.0.0 in Docker. Only `localhost`/`127.0.0.1`/`[::1]` are accepted unless `MOLTALK_ALLOWED_HOSTS`/`MOLTALK_ALLOWED_ORIGINS` are set. Request bodies are limited to 64 KB (`MOLTALK_MAX_BODY_BYTES`). There is no authentication, so keep it on loopback.
- **PubChem:** `resolve_name` contacts PubChem only when the model sets `allow_network=true`, and the compound name is then sent to PubChem. Nothing else leaves your computer except tunnel traffic to OpenAI.

## 5. Local HTTP and Docker (testing only)

```bash
.venv/bin/moltalk --transport streamable-http          # http://127.0.0.1:8000/mcp
npx @modelcontextprotocol/inspector@latest                 # point it at the URL above

docker build -t moltalk .
docker run --rm -p 127.0.0.1:8000:8000 --read-only --tmpfs /tmp --memory 2g --cpus 2 --pids-limit 128 moltalk
```

Note on the original V1 Dockerfile: it could never start, because `python:3.12-slim` lacks `libXrender`/`libXext`/`libexpat`, which RDKit's drawing module needs. That is fixed here.

## 6. Troubleshooting

- `scripts/status.sh` shows *not ready*: check `journalctl --user -u moltalk-tunnel`. "control plane API key" errors mean the key file is missing or wrong; rerun `store-tunnel-key.sh`. Permission errors mean the key's principal lacks Tunnels Read + Use.
- The plugin shows no tools or fails to create: make sure the tunnel is running *while* you create it, and that the tunnel ID in the plugin matches `~/.config/moltalk/tunnel-profiles/moltalk.yaml`.
- The drawing doesn't appear but the answer does: refresh the plugin and start a new chat. The model still receives the full chemistry data.
- `scripts/local-tunnel-test.sh` runs the tunnel path locally without credentials and prints a local MCP URL for the MCP Inspector.

## 7. Chemistry scope

Indices are zero-based **input-SMILES** atom indices, not IUPAC locants and not canonical-SMILES order. Wedge/dash belongs to a particular 2D depiction (a wedge starts at the stereocenter and points toward the viewer). It is not a synonym for R/S. In the cyclohexane example, both the R and the S centre carry wedges.

RDKit does not parse IUPAC names. The offline dictionary covers water, ethanol, acetone, benzene, acetic acid, aspirin and caffeine; other names need PubChem, and PubChem records may leave stereochemistry unspecified. For example, "lactic acid" (CID 612) is unspecified and "glucose" (CID 5793) is partially specified. That status is reported explicitly in `stereo_summary` and the warning. logP is a calculated estimate. Functional groups come from a small SMARTS motif dictionary (ester oxygens are not reported as ethers) plus aromatic-ring detection; this is not an exhaustive classification. Not supported: reaction prediction, conformer energies, chair conformations, structure editing, general IUPAC parsing.

Sources: [OpenAI plugins: connect and test](https://developers.openai.com/plugins/deploy/connect-chatgpt), [Add UI to your MCP server](https://developers.openai.com/plugins/build/chatgpt-ui), [Plugins reference](https://developers.openai.com/plugins/reference), [Secure MCP Tunnel](https://developers.openai.com/api/docs/guides/secure-mcp-tunnels), [Developer mode](https://developers.openai.com/api/docs/guides/developer-mode), [openai/tunnel-client](https://github.com/openai/tunnel-client), [RDKit](https://www.rdkit.org/docs/), [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk).
