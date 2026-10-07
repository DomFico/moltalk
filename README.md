# MolTalk — draw, rotate and name organic molecules

<img src="assets/banner.svg" alt="MolTalk logo: dodecahedrane drawn as a view of its 3D shape" width="100%">

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
| Test suite: 106 tests (104 offline by default; 2 live PubChem tests need `MOLTALK_NETWORK_TESTS=1`), including live PubChem and real-browser tests of the drawing, its 3D rotation, zoom, hydrogens and a phone touch screen | **Verified** (`pytest`) |
| Tools and widget through official MCP clients (Python SDK over stdio/HTTP, MCP Inspector CLI) | **Verified** |
| UI resource metadata (`ui://widget/molecule-v32.html`, `text/html;profile=mcp-app`, `_meta.ui.resourceUri` + `openai/outputTemplate`) | **Verified** |
| Widget rendering in a sandboxed iframe through the MCP Apps bridge, plus a follow-up `tools/call` from the widget | **Verified** in Chrome with a host harness that imitates ChatGPT (`tests/test_widget_ui.py`, `tests/test_widget_3d.py`) |
| tunnel-client → stdio server path | **Verified** with tunnel-client's local control plane (`scripts/local-tunnel-test.sh`) |
| Name resolution: bundled library (22,038 compounds), OPSIN for systematic names, live PubChem fallback | **Verified** (`tests/test_library.py`; PubChem with `MOLTALK_NETWORK_TESTS=1`) |
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
| Draw aspirin and identify its functional groups. | `resolve_name` (MolTalk library, offline) then `draw_molecule`. Inline drawing; groups: carboxylic acid, ester, aromatic ring |
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

After changing tool metadata or the widget: restart (`scripts/start.sh`), then in ChatGPT open the plugin under **Plugins** → **Refresh**, and start a new chat. If you change the widget substantially, bump the version in its URI (`molecule-v32` → `v33`), because ChatGPT caches templates.

### If the Tunnel option is unavailable

A public HTTPS endpoint is the documented alternative. ChatGPT accepts only **OAuth** or **no auth** for MCP servers, and static API keys are not supported. This project **does not** expose its unauthenticated HTTP server publicly, so going public needs one more piece of work: an OAuth 2.1 authorization layer, using either the MCP SDK's auth support or a hosted identity provider with Dynamic Client Registration. That layer would sit in front of the Docker image on a small HTTPS host (typically about US$5–10 per month), with `MOLTALK_ALLOWED_HOSTS` set to that domain. Decide on a host and identity provider before this is built.

## 3. Tools

| Tool | Shows UI | Returns to the model |
|---|---|---|
| `analyze_molecule(smiles)` | no | identifiers, properties, CIP R/S, E/Z, `stereo_summary`, functional-group motifs |
| `draw_molecule(smiles, label?, width, height, atom_indices, include_svg)` | **yes** | full analysis + `depicted_stereo_bonds` (each wedge/dash with a plain-language explanation) |
| `enumerate_stereoisomers(smiles, limit)` | **yes** (grid) | isomers with CIP labels, verified axial Ra/Sa (atropisomers as CXSMILES), `achiral` (meso) flag, `enantiomer_index` |
| `find_substructure(smiles, smarts)` | no | chirality-aware matches (≤100) |
| `draw_named_molecule(name, numbering, hydrogens, allow_network)` | yes | resolves a compound name (as `resolve_name`) and draws it in one call; the preferred tool for "draw X" |
| `export_structure(smiles, format, coordinates, name, hydrogens)` | yes (download card) | a structure file: `cdxml` (ChemDraw, default), `mol`, `sdf`, `pdb`, `xyz` or `smiles`; 2D (the drawing's layout) or 3D (one conformer, explicit H); warnings when a 3D file fixes unspecified stereo |
| `conformer_3d(smiles)` | viewer-only | one calculated 3D conformer aligned to the flat drawing, for rotating it (hidden from the model) |
| `resolve_name(name, allow_network)` | no | structure from the bundled library, OPSIN or (opt-in) PubChem, with source, CID, title, IUPAC name, `stereo_summary` and warnings; ambiguous or unknown names are errors |

**UI component.** `moltalk/widget/molecule.html` is registered as `ui://widget/molecule-v32.html` (`text/html;profile=mcp-app`) with an empty CSP allowlist, so it makes no network requests. It uses the MCP Apps bridge (`ui/initialize`, `ui/notifications/tool-result`, `tools/call`, `ui/notifications/size-changed`) and falls back to `window.openai`. The SVG is sent only in the result's `_meta`, which ChatGPT passes to the widget but not to the model; the chemistry data stays in `structuredContent` for the model. The SVG is sanitized before insertion. Pass `include_svg=true` for clients without UI support.

**Rotate the drawing in 3D.** Grab the flat drawing and drag (one finger on a phone). Over the first ~90 px it lifts off the page into 3D and keeps rotating as a chemical drawing:
- **Drawing style:** flat colour, implicit carbons, element labels with their hydrogens (OH/HO, NH₂) in RDKit's colours, half-coloured bonds, ring double bonds on the inner side, and gaps where a bond passes in front of another.
- **Zoom and pan:** **−/+** (1.2× steps that glide rather than snap), Ctrl/Cmd+scroll or a trackpad pinch (continuous), two-finger pinch on touch screens; Shift+drag or a two-finger drag pans. Zoom works on the flat drawing too. Plain scrolling still scrolls the chat.
- **Phones:** touches on the drawing belong to the drawing, so a swipe rotates it instead of scrolling the page or switching chats. To scroll past it, swipe outside the drawing.
- **Enter viewer** asks ChatGPT for its full-window viewer. There the widget fits the window exactly with no scrolling: the drawing takes the free space, and the notes and lists are hidden. ChatGPT draws its own bars over the viewer (a title bar on phones, and the "Ask ChatGPT" box at the bottom everywhere). The widget keeps clear of them, using the host's safe-area insets or minimums measured from screenshots if those are larger: 104/128 px top/bottom on phones, 88 px at the bottom on desktop. Our controls sit at the top, under the header; on touch devices the clearances are 132/140 px. Returning inline, the widget re-reports its height. On iPhone, ChatGPT shows the app in the chat as a collapsed "Open App" card and adds its own pull-up sheet; both are host behaviour. **Known iOS host bug:** closing the fullscreen view with ChatGPT's own X can make the inline widget disappear (open OpenAI community report, unresolved as of 1 Oct 2026). The widget does not call `requestClose()`; it offers **Exit viewer**, which asks the host for `inline` through `requestDisplayMode` / `ui/request-display-mode`. It declares `availableDisplayModes: [inline, fullscreen]`, reports its height only when inline and never below 280 px, and re-reports it after every mode change.
- **Controls row:** one line under the drawing. The hint text is on the left (drag to rotate, Ctrl+scroll or pinch to zoom, Shift+drag or two fingers to pan). On the right, flush: **Export SVG / PNG**, **Enter viewer** (or **Exit viewer**) and **2D**. On a phone the hint wraps above the buttons, and the buttons stay together. There are no +/− buttons; the + and − keys still zoom.
- **2D** (or Escape) restores zoom and settles the drawing back onto the exact RDKit picture. Arrow keys rotate.
- **Fitting the rotated view:** the 3D model shrinks gradually as it lifts, so that it stays inside the canvas in every orientation; line widths shrink with it, so bonds keep the flat drawing's weight. A small molecule's flat drawing fills the canvas, and a longer 3D bond such as C–Br would otherwise leave the frame.
- **Inline size:** the widget fits inside the host's inline height limit (`containerDimensions.maxHeight` or `window.openai.maxHeight`), so the host never adds a scroll bar. Inline, the drawing is never taller than 420 px: its natural height grows with the window width, and in a wide Codex desktop window it was taller than the host would show. Some hosts cap the frame without reporting a limit. If the frame stays shorter than the height the widget asked for over two checks a second apart, the widget takes the frame height as the limit and shrinks the drawing to fit. The rotated-view note, chips, wedge list and notes share one small box that scrolls on its own. Under a tight limit the drawing shrinks rather than the controls.
- **Show hydrogens / atom indices** redraw without changing the view (the rotated picture stays on screen while the redraw loads, and the other hydrogen setting's 3D data is fetched in advance): rotation, on-screen size, zoom and pan are kept. The conformer is chosen from heavy atoms only, so it is the same shape with or without hydrogens. It is aligned to whichever drawing is on screen: RDKit lays out the explicit-H drawing in its own orientation (LSD's is turned 180°), and `to_heavy_frame` converts the remembered rotation between the two, so a toggle never flips or moves the molecule.
- **Geometry:** RDKit ETKDGv3 followed by MMFF94 (or UFF); fullerene-like cages start from a spectral sphere.
- **Matching the flat drawing:** up to 8 conformers (32 if needed) are generated, and the one best matching the flat drawing is kept. Each is placed with a proper rotation only, never a mirror image, so R/S cannot change. Stereocentres are weighted heavily in the fit, so wedges point toward the viewer and dashes away when the lift begins. Each drawing makes one 3D request.
- **Caveats shown in the view:** it is labelled as one calculated conformer, and arbitrary configurations at unspecified stereocentres are marked "?".
- **Not yet done:** wedges/dashes are not redrawn during rotation.

**IUPAC name and numbering.** Under the title the widget shows **SMILES:** and **IUPAC:**, and **Numbers: IUPAC | Index | Off** chooses the atom labels. IUPAC is the default; the stereocentre chips and the model's text use locants too (e.g. "C6a: R"). The pipeline (`moltalk/naming.py`) never guesses:
1. **Name:** looked up by the molecule's full InChIKey, an exact structure match including stereochemistry: first in MolTalk's bundled library (offline, with the numbering precomputed), then in PubChem. The names are PubChem's, generated by OpenEye Lexichem. A structure in neither gets no name ("not guessed") and keeps atom indices. OPSIN cannot generate names, so it is not a substitute here (see *Name resolution*).
2. **Checking the name:** OPSIN 2.9.0 (open source, `vendor/`, SHA-256 matches the GitHub release, needs Java) rebuilds each candidate name. A name is used only if its structure equals ours (connectivity and charges), which drops PubChem's zwitterion and salt duplicates. If different names still fit, none is used. A redundant `cis-`/`trans-` before R/S descriptors is removed for parsing only.
3. **Numbering:** OPSIN numbers every substituent in its own scheme. The parent is the scheme that contains all top-level locants of the name and whose bonds to other carbon schemes leave from those positions; ties are broken by which positions actually carry substituents. A remaining tie means no numbering. Only parent atoms are numbered. Symmetry-equivalent numberings are equally valid; one is used only if every R/S descriptor in the name matches our own CIP label at that locant.
4. **Fallback:** when any step fails, the widget says why, keeps the name if known, and uses atom indices.

Tested offline on 25 typical course molecules (chains, rings, aromatics, steroids, sugars, LSD). All parents were numbered correctly, except acetaminophen and ethyl acetate, whose names give no positional clue to the parent (refused), and heme, where PubChem's records disagree (refused). Privacy: a molecule that is not in the library has its InChIKey sent to PubChem (cached); set `MOLTALK_OFFLINE=1` to turn this off (library names and numbering keep working). The Docker image includes Java and OPSIN.

**Name resolution** (`resolve_name`, name → structure). Tried in this order; the first that answers wins:
1. **MolTalk's bundled library** (`moltalk/data/compounds.json.gz`, about 3 MB, offline). 22,038 compounds: every compound on Wikidata that has a PubChem CID and an English Wikipedia article, plus about 460 common teaching names from `scripts/library_seed_names.txt` (sugars such as glucose, heme, cofactors, terpenes, common drugs, reagents). It covers drugs, natural products, metabolites, amino acids, solvents and reagents, and holds about 200,000 names (352 of them shared by different structures, reported as ambiguous). Each entry has its CID, PubChem's isomeric SMILES and title, synonyms, PubChem's IUPAC name, whether OPSIN rebuilds exactly that structure from the name (`iupac_verified`, about 89 %), and precomputed parent locants (about 60 %). Lookup ignores case, spacing and dash style. Names come in tiers, most trusted first: the curated overrides and the seed list; the compound's own Wikidata label, PubChem title or IUPAC name; Wikidata aliases; PubChem synonyms. The first tier that knows the name decides, because the lower tiers contain errors (Wikidata lists "ozone" as an alias of phencyclidine, and both sources list "LSD" for lysergic acid). A name that still points to different structures within that tier (e.g. "lye": NaOH or KOH) is an error that lists the candidates; none is picked.
   - **Curated overrides** (`scripts/library_curated.json`) cover names where database naming is messy. "heme" is heme b, Fe(II) protoporphyrin IX (PubChem titles one copy of it "Hemin"); "hemin" is the Fe(III) chloride. Also NAD+/NADH, NADP+/NADPH, FAD/FADH2, FMN, CoA, acetyl-CoA, ATP/ADP/AMP and cAMP. Each has a corrected title and a `note` returned with the structure (for example, the phosphates' protonation state).
   - **Stereo not stated:** a bare name is refused when it names a D/L family. That applies when the record's own title is D-/L-/DL-qualified (L-Alanine, D-Fructose), or another stereoisomer's title is (D-Alanine). So "alanine", "serine", "glucose" and "fructose" ask which one is meant, while "L-alanine" and "D-glucose" resolve. Named natural products and drugs ("morphine", "epinephrine") name one stereoisomer, and resolve. Racemic or unspecified records ("lactic acid", "ibuprofen") resolve with a warning that stereo is unspecified. Synonyms do not trigger this (PubChem lists "d-morphine" for (+)-morphine).
2. **OPSIN** (offline): systematic IUPAC names, including R/S and E/Z, e.g. `(2R)-4-chloro-2-methylheptan-3-one`. It cannot read trivial names. The result says the structure came from parsing the name, and warns about stereo the name leaves open.
3. **PubChem** (best effort, only when the model passes `allow_network=true`; the name is sent to PubChem). On Cloud Run, PubChem often answers "server busy", because Cloud Run's shared outgoing addresses carry other people's traffic (see Cloud NAT below).
4. **Otherwise:** an error says the name was not found and no structure was assumed.

**Label check.** When the model draws a structure with `label=` set to a compound name, `draw_molecule` looks the label up in the library. It reports `label_check`: `matches`, `stereo differs`, `mismatch` (the drawing is a different compound; also shown in the widget and put first in the model's text) or `unverified` (the label is not a library name). So a SMILES the model writes from memory under a known name is caught. The server instructions also tell the model not to write SMILES from memory when `resolve_name` refuses a name, unless the user asks it to, and then to say it is unverified.

**Stereochemistry stated in a name is verified.** If a name states R/S, E/Z or axial descriptors, such as (Ra), (Sa), (P) or (M), the structure it resolves to must encode at least as many of each, with matching CIP labels; otherwise the name is rejected as unverified, and nothing is drawn. This applies to the library, OPSIN and PubChem alike. "(S)-BINAP" is rejected, because the records it maps to carry no axial stereo; "(R)-carvone" passes only if the structure's centre really is R. Optical rotation, (+) or (−), is not a CIP descriptor and is not checked.

**Atropisomers (axial chirality).** Supported generically (`moltalk/atropisomer.py`), not by special-casing BINAP. Everything else keeps RDKit's legacy stereo perception, so ordinary R/S and E/Z results are unchanged. On 400 library compounds, every canonical SMILES and every R/S and E/Z label matched the previous release; the only differences were the newly reported hindered axes in three polysubstituted biaryls, including gossypol.
- **Input and output:** plain SMILES cannot state a twist. It travels as CXSMILES with 2D coordinates and a wedge at an axis atom (`...|(coords),wU:20.21|`), which RDKit reads into the axis bond's stereo. When an axis is specified, `canonical_smiles` is that CXSMILES, written in canonical atom order and checked by reading it back. MOL, SDF (2D and 3D), CDXML and SMILES exports keep the twist (round-trip tested).
- **Finding open axes:** RDKit does not list unspecified axes. MolTalk uses a stated rule of thumb, not a barrier calculation. An axis counts if it is a non-ring single bond between two sp2 ring atoms, each end has two ring neighbours that rank differently, and at least 3 of the 4 ortho positions carry a substituent or a fused ring. All 4 are required when an axis atom is in a five-membered ring. Such axes are reported as `stereo_summary.unspecified_axes`, with a warning. BINAP, BINOL, gossypol and tetra-ortho-substituted biaryls qualify. Biphenyl, 2,2'-dimethylbiphenyl, diphenic acid and 2,6-dimethylphenyl (whose two sides are alike) do not.
- **Isomers:** `enumerate_stereoisomers` adds both twists of every open axis, alongside the R/S and E/Z isomers. Each atropisomer is given as CXSMILES with its enantiomer index, and taking the mirror image swaps the twist.
- **Verified descriptors only:** RDKit's CIP labeller names an axis P or M. A 3D model is then measured independently in two ways:
  - the helicity rule: the sign of the torsion between the top-ranked group on each side;
  - the CIP axial rule: seen along the axis, the turn from near group 1 to near group 2 to far group 3, where clockwise is Ra.
  
  The two must satisfy Ra = M and Sa = P, and both must agree with RDKit. Only then is `axial_stereo[]` marked `verified`, with `cip` (Ra/Sa) and `helicity` (P/M). Otherwise no descriptor is shown anywhere, and RDKit's raw P/M is removed from the drawing.
- **Drawing:** one ring bond at the axis is wedged or hashed. `depicted_stereo_bonds` explains it as the axis twist, not a stereocentre. Under the Stereo control, the axis carries (Ra)/(Sa), or (?) when it is open.
- **3D:** ETKDG embeds the specified twist. After the tail swing and chain straightening, the helicity is checked and, if a rotation turned it over, set back. The rotated view tags the axis (Ra)/(Sa). An open axis gets (arb. Ra), and the 3D export says which twist it chose.
- **Names stay strict:** an axial descriptor in a name is never resolved, so "(R)-BINAP", "(S)-BINOL" and "(Ra)-BINAP" fail closed. Plain "BINAP" resolves with the axis reported open, and `enumerate_stereoisomers` then shows both atropisomers.
- **C–N axes:** the same machinery covers an aryl ring on a tertiary amide nitrogen whose two ortho positions carry different substituents. Metolachlor is the reference case, with four stereoisomers (R/S × Ra/Sa). Ordinary amides are never treated as axes; their slowed C(O)–N rotation is not atropisomerism. Acetanilide, N-methylacetanilides, symmetric 2,6-dimethyl anilides, lidocaine, DMF and imides are all negative controls.
- **Not covered:** InChI does not encode axial chirality, so both atropisomers share an InChIKey (a warning says so).

**Stereogenic units (allenes, spiro compounds, helicenes).** Every kind of stereochemistry MolTalk handles is reported in one form, `analysis.stereo_units` (`moltalk/stereounits.py`). Each entry has:
- `type`: tetrahedral, double_bond, axis, allene, spiro or helix;
- `atoms` and `bonds`;
- `specified`, `configuration`, `possible_configurations` and `configuration_source`;
- `verification` and `descriptor`;
- `stability`.

Three questions are kept apart:
1. Is the unit stereogenic? This is a graph analysis.
2. Does the input specify its configuration?
3. Has the descriptor been checked against the IUPAC definition?

Stability is a separate field. Being stereogenic does not mean the isomers can be isolated.

- **Allenes and even cumulenes** (IUPAC P-93.4.2.2): a chain with an even number of cumulated C=C bonds whose two ends each carry two different groups. 1,3-Dichloroallene and penta-2,3-diene qualify. Propadiene, 1,1-disubstituted allenes and odd cumulenes (planar, so E/Z) do not.
  - Descriptor: M/P, the form IUPAC prefers, with Ra/Sa as the alternative.
  - How it is checked: both rules are measured on 3D coordinates. The helicity rule uses the torsion from the near top-ranked group to the far one. The CIP axial rule is the tetrahedral rule on the elongated tetrahedron; its sign is calibrated against RDKit on (S)-alanine. The two must agree as IUPAC states (M = Ra, matching its example (1M)-1,3-dichloropropa-1,2-diene = (1Ra)), or no descriptor is shown.
- **Spiro compounds:**
  - **Xabcd** (four different ring neighbours): an ordinary R/S centre.
  - **Xaabb** (2,6-dichlorospiro[3.3]heptane, Fecht's acid): IUPAC names these with R/S at the spiro and ring atoms. Only RDKit's newer stereo perception keeps them, so it is used for these molecules alone. M/P is also given for two four-membered rings, and the tests reproduce IUPAC's example (2R,4S,6R) = (2P).
  - **Xabab** (spiro[4.4]nonane-1,6-dione): IUPAC's descriptor needs a special CIP digraph (P-93.5.3.2) that neither RDKit nor MolTalk implements. The configuration is represented, enumerated, drawn and built in 3D, but no descriptor is given, and the result says why.
  - Achiral spiranes report nothing: spiropentane, spiropentadiene, spiro[3.3]heptane and 2-chlorospiro[3.3]heptane.
- **Helicenes:** five or more six-membered rings ortho-fused in a chain that always turns the same way, with no ring atom in three rings.
  - Negative controls: [4]helicene, picene, pentacene and C60.
  - Descriptor: P or M, from the handedness of the helix through the ring centres (right-handed = P), as IUPAC's (P)/(M)-hexahelicene.
  - Stability: [6] and higher are configurationally stable; [5]helicene racemises slowly at room temperature.
- **Representation:** SMILES cannot state an allene's, helicene's or Xabab spiro atom's configuration. MolTalk does not invent a SMILES extension. The configuration travels as standard CXSMILES with 3D coordinates, or as a 3D MOL/SDF file, and is measured from them. A specified configuration makes `canonical_smiles` that CXSMILES. `enumerate_stereoisomers` returns each configuration in that form, and enantiomer pairs are matched.
- **Drawing:**
  - Allenes: one end's substituents in the paper and the other end's wedged and hashed, with an (M)/(P) label.
  - Xabab spiro atoms: wedge and hash on the second ring's bonds at the spiro atom, which IUPAC's drawing rules allow.
  - Every wedge choice is checked on a pseudo-3D version of the drawing (wedge = toward the viewer) with the same measurement as real coordinates, so a drawing cannot contradict its configuration.
  - Helicenes: no wedges, which would imply stereocentres. They are drawn as a view of the helix, the textbook picture with the end rings side by side, labelled (P)/(M); a projection has no depth and implies neither hand.
- **3D:**
  - ETKDG does not know these configurations, so conformers that do not match are dropped, or mirrored when nothing else in the molecule is chiral.
  - MMFF94 and UFF both flatten allenes (1,3-dichloroallene went from 84° to 180°), so each allene's ends are held perpendicular during the clean-up. This also fixes plain allenes, which were drawn flat before.
  - Helicenes are embedded with ETDG: ETKDG's flat-aromatic-ring knowledge cannot embed [6]helicene within 10 s.
  - The spiro squaring step keeps each ring on its side, so it cannot invert a chiral spirane.
- **Export:** these configurations survive CXSMILES, MOL/SDF 3D and PDB (round-trip tested). 2D MOL/SDF and CDXML cannot carry them; the export says so instead of silently dropping them. PDB and XYZ have no bond orders, so a reader cannot re-perceive the units, and the export warns. Xaabb spiranes and biaryl/C–N axes survive 2D and 3D MOL/SDF and CDXML.
- **Names stay strict:** "(P)-hexahelicene" and "(Ra)-1,3-dichloropropa-1,2-diene" are not resolved, like "(R)-BINAP".
- **Not supported:**
  - Combining a hindered axis and an allene/helix/spiro unit in one enumeration (refused with a message).
  - Heterohelicenes with five-membered rings.
  - Planar chirality (cyclophanes).
  - Octahedral and other non-tetrahedral centres.

OPSIN replaces PubChem only for *name → structure* of systematic names. It cannot *generate* a name for a structure, so structure → IUPAC name still comes from the library or PubChem. A novel structure found in neither is shown with atom indices and no name, rather than a guessed one.

Rebuild the library with `.venv/bin/python scripts/build_library.py` (network, Java and OPSIN; about 10 minutes). RDKit must recompute PubChem's InChIKey from its SMILES, or the entry is dropped; compounds above 150 heavy atoms are also dropped. The current build kept 22,038 entries. Any seed name PubChem cannot find stops the build.

**Stereo labels.** **Stereo: Specified | All | Off** (default Off) controls only the labels in the drawing; the analysis, the chips and the model's text are unchanged. The setting survives redraws, and switching is instant (RDKit tags its labels `CIP_Code`).
- **Specified:** R/S and E/Z only where the input fixes the configuration, in both the flat and the rotated view (E/Z sits beside its double bond).
- **All:** also marks stereocentres the input leaves open. Open double bonds are marked too. In the flat drawing that is **(?)**, because there is no configuration to show. RDKit draws it as a stereo annotation, so it has the same size and collision-avoiding placement as R/S; the viewer attributes each annotation glyph to the nearest labelled atom or double bond to show or hide it. In the rotated view it is **(arb. R)**, the configuration this one conformer happens to have, with the tooltip "Configuration chosen for this displayed conformer; input stereochemistry unspecified". MolTalk does not write R\*/S\*, because in IUPAC usage that means relative configuration.
- **Off:** no stereo labels.

**Export image.** **SVG** and **PNG** under the drawing save exactly what is on screen: the rotation, zoom, hydrogens, numbers, stereo labels and lone pairs. The image is cropped to the molecule with a small margin and has a transparent background, with dark ink (the drawing is always dark-on-white, whatever the theme). PNG is 3× the on-screen size. In the rotated view, bond-crossing gaps and label backings are painted in the background colour on screen; in the export they become real cut-outs (SVG masks), so nothing shows as white on a coloured slide.

**Structure files.** Ask for a file ("give me a ChemDraw file of that") and the model calls `export_structure`. The result is a card with a **Download** button.
- **Formats:** CDXML (ChemDraw opens it directly), MOL, SDF, PDB, XYZ and SMILES, all written by RDKit; there is no Open Babel. Binary `.cdx` is refused, because RDKit's Python writer for it is broken; CDXML replaces it.
- **2D or 3D:** 2D files use the drawing's own layout, with stereo as wedges. 3D files (MOL/SDF with `coordinates="3d"`, PDB, XYZ) use the conformer the viewer rotates, with explicit hydrogens. A 3D file has to fix a configuration at every stereocentre, so when the input leaves some open, the result, the card and the SDF's `MOLTALK_WARNING` field say which ones were chosen arbitrarily.
- **What is exported:** the structure as given (atom order, charges, stereo), not drawing-only additions such as the Fe–N bonds drawn for heme.
- **Tested:** every format was round-tripped through RDKit's readers with the stereo intact. The CDXML has not been opened in ChemDraw itself.

**Downloads and privacy.** Both kinds of export use MCP Apps' standard `ui/download-file` when the host advertises `downloadFile` (the host usually asks the user to confirm). Otherwise the widget tries an ordinary browser download, and also opens a small panel: right-click or long-press the image to save it, or copy the file text. That panel works even where the sandbox blocks downloads. The structure never goes into a URL; the file content travels in the tool result (`_meta`, plus a standard embedded resource for hosts without the viewer). Which route ChatGPT and Claude take has to be checked in each app.

**Labels and lone pairs in the rotated view.** Lone pairs, charges, R/S and E/Z labels, atom numbers and double bonds are fixed in 3D and turn with the molecule, as hydrogens do, instead of keeping a fixed place on screen.
- **Lone pairs:** each pair lies along its VSEPR direction, out past the edge of its label, with its two dots split along a fixed 3D axis. A pair pointing at the viewer sits over its atom, as a hydrogen would.
- **Formal charges are typography, not geometry:** after projection, a charge is placed as an upright superscript just outside the atom's on-screen footprint (its label and lone-pair dots). It prefers the upper right, avoids bonds, and keeps its spot while that stays free, so it never orbits the atom. Charges never move between resonance-equivalent atoms: the Lewis structure as given is drawn.
- **R/S, E/Z and numbers:** they sit in the atom's most open 3D direction (for a stereocentre, usually where its hidden hydrogen points). A minimum offset keeps them off the label when that direction points at the viewer.
- **Double bonds:** the second line lies in the plane of the double bond (toward its ring or substituents). Seen edge-on, the two lines close up to 35% of the gap, but never merge. E/Z labels sit on the bond's open side.
- **Flat drawing (a Lewis drawing, not a projection):** every pair on an atom sits at the same radius from its element symbol: a circle for one-letter symbols, and for wider ones (Cl, Br) an ellipse with even clearance, centred on the symbol's own glyphs. Only the pairs' angles come from the bonds and any H text. With lone pairs shown, the charge is written outside them like an exponent: RDKit's raised sign is hidden, and the charge goes to the conventional upper right. It moves outward past a pair in the way, and only bonds or H text send it to another angle. The rotated view keeps the physical 3D directions.

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
- **Independent of SMILES order:** the model draws `resolve_name`'s canonical SMILES, whose atom order differs from the library's. Fused-ring substituents (chlorophyll's methyl ester) are now pointed outward like core substituents, so the pinned layout is accepted for either order. Previously the canonical order gave a round ring with Mg–N stretched 2.6×.
- **Any porphyrin-type macrocycle, metal or free base:** any molecule with a ring of 12 or more atoms also tries RDKit's ring templates, which give textbook porphyrins, corroles and corrins instead of a large circle. A porphyrin-type core, whose bridges may be carbon or nitrogen (phthalocyanines, porphyrazines), is then pinned to the porphine square if the templates missed it. If side groups still collide, every side group is pointed outward. The textbook layout is kept unless it adds bond crossings, stretched bonds, or more than two crowded atom pairs. Tested on 31 tetrapyrroles: Ni, Cu, Zn, Fe, Co, Mn and Sn phthalocyanines, free phthalocyanine, naphthalocyanine, porphine, tetraphenylporphyrin and its Zn complex, octaethylporphyrin and its Pt complex, protoporphyrin, coproporphyrin, temoporfin, chlorin e6, pheophytin, chlorophyll a and b, bacteriochlorophyll, heme, hemin, heme a, siroheme, Zn protoporphyrin, F430, corrole and cyanocobalamin. Also checked on 19 other macrocycles (crown ethers, macrolides, cyclic peptides, cyclodextrins): none got worse.
- **Metals with their own ligands, and separate ligands:** the metal may be stored as its own small piece with ligands (tin phthalocyanine dichloride: Cl–Sn–Cl beside a free-base ring). A separately listed anionic or radical ligand (B12's cyanide) is bonded to the metal. Bonds become covalent, with ring N–H hydrogens replaced as in the complex, only when the metal's spare charge or valence is matched exactly; otherwise all bonds are dative and the charges stay as stored. In 2D, axial ligands sit on the diagonals between the M–N bonds; in 3D, above and below the ring. Fused rings are placed outward ring by ring, so naphthalocyanine's outer rings are regular. With hydrogens shown, a macrocycle keeps its layout.
- **Hydroporphyrins with extra fused rings** (coenzyme F430's lactam and cyclohexanone, chlorophyll's ring E) defeat RDKit's template: the macrocycle folded and the metal bonds stretched to 2.6×. The porphyrin core (four N-rings joined by four bridge carbons, matched ignoring bond orders) is then pinned to the porphine square. Rings fused to it are placed on an arc of equal bonds outside the shared atoms, and lone substituents on the core point outward. This applies only when the template result is distorted, and only if the result adds no crossings, overlaps or stretched bonds.
- **3D shape:** MMFF bends the charged porphyrin form out of plane (about 0.6–0.8 Å RMS), so the ring is embedded as the neutral free base and the flattest candidate is kept (about 0.13 Å for heme, in about 1 s). The metal is placed in the plane at the centre of the four nitrogens.
- **The analysis is unchanged:** formula, charges and the "disconnected fragments" warning describe the source structure. The drawing says what was added.
- **Charge-separated sources** (for example heme written as [Fe-2] bonded to two [N+]) are drawn the same way, with dative N→Fe bonds and neutral atoms.
- **Metals are never treated as stereocentres.** RDKit lists a four-coordinate metal as a possible tetrahedral centre, but square-planar Fe is not one.
- **Not handled:** other metal complexes, such as ferrocene's sandwich bonding or cisplatin's ammine ligands, are drawn as stored.

**Cages.** In bridged or cage-like ring systems (bridgeheads, or atoms in three or more rings, as in cubane), an unspecified element counts only if its R/S differs among the stereoisomers that can actually be built in 3D. Each candidate is tried with a small ETKDG budget and two seeds, mirror images are added automatically, and the check is capped at 256 candidates. Substituted cubanes therefore report no stereocentres, while norbornan-2-ol (3) and camphor (2) keep theirs. The rotated view shows R/S only at real stereocentres.

**Alignment with the drawing.** The 3D model is rotated onto the flat drawing by a robust weighted fit: atoms that cannot match, such as a long tail folded differently in 3D, are progressively down-weighted, so the rigid part decides the rotation. For metal chelates the macrocycle and the metal anchor the fit. Before this, chlorophyll's phytyl tail turned the ring 2–3 bond lengths away from the drawing, so lifting looked scrambled; its ring now starts within 0.12. Flexible arms still swing as the drawing lifts, because the model is one real conformer.

**Tails.** A long chain hanging off a ring system can leave it pointing back across the ring (chlorophyll's phytyl ester). The first four torsions of its attachment are scanned (60/180/300°). The clash-free combination that best matches the flat drawing is kept, then relaxed. Chlorophyll a's fit to its drawing went from 4.8 to 2.6, and its tail now runs away from the ring.

**Repeated draws.** ChatGPT sometimes re-sends the first `draw_molecule` of a chat, with identical arguments, about 10 s later, although the first succeeded. When a ChatGPT session's draw repeats its immediately previous draw within 20 s, the second viewer shows one line, "Same drawing as above". The model still gets the full result, and a viewer toggling a setting back and forth is never collapsed.

**Metal–N cavity.** Embedding the ring as the free base leaves its central cavity lopsided (F430 had Ni–N of 1.78–2.44 Å once the metal was put at the centre). The chosen conformer is re-minimised with the four donor N restrained to the square the metal needs: cis N···N = d√2 and trans = 2d, with d = 2.05 Å for Mg, 2.0 Å for Fe and 1.95 Å for Ni. This is the last step, so no later minimisation undoes it. F430 now has Ni–N 1.93–1.99 Å, and chlorophyll Mg–N 2.05–2.09 Å. F430's ring stays ruffled, as in its crystal structures.

**Long chains.** ETKDG makes long open chains crumpled, often folded back over the molecule. For the chosen conformer, every torsion along open sp³ carbon chains (four or more carbons) is set to anti, the extended zigzag the flat drawing also shows, and the force field relaxes it. The result is kept only if it matches the drawing better. Palmitic acid now lifts almost exactly onto its drawing. The rotated view is sized by the 90th-percentile atom distance, so chlorophyll's tail does not shrink the ring.

**Request log.** On Cloud Run, the gateway logs each request's JSON-RPC method and tool name (never the arguments), the client's `mcp-protocol-version`, and for any 4xx response the server's reason.

**Speed.** On Cloud Run's single CPU, a 3D model costs roughly 2.5× desktop time. Molecules over 50 heavy atoms with 5 or more rotatable bonds (chlorophyll, erythromycin) get one round of 2 candidate conformers and a shorter force-field clean-up (300 instead of 1000 iterations). Rigid ones (naphthalocyanine) keep the full clean-up, or they come out bowed; chlorophyll a went from 5.3 s to 1.7 s on a desktop. The extra "poor fit" round of 32 conformers is limited to molecules of up to 15 heavy atoms; cholesterol had been paying for it every time (5.1 s → 1.0 s). The CIP labeller only labels real stereocentres and stereo double bonds, with an iteration cap: unbounded, it never finished on dodecahedrane's 3D model, which timed out at 25 s.

**Defaults.** Drawings start clean: no atom numbers (`numbering` defaults to `none`) and no stereo labels (the viewer's Stereo control starts at Off). Both can be switched on in the viewer.

**Strain and spiro centres.** Candidates more than 8 kcal/mol above the best are dropped before matching the drawing; one spiropentane candidate had collapsed to 45° and fitted the bowtie drawing best. The filter is skipped for porphyrin-type macrocycles, where MMFF rates a bowed naphthalocyanine below the flat one. A spiro carbon joining two 3- or 4-membered rings is set geometrically to the real D₂d shape (perpendicular rings on a straight axis). MMFF twists spiropentadiene to 56° and even rates that lower in energy.

**Metal complexes stored as pieces** (`moltalk/complexes.py`). PubChem, and so the library, stores most catalysts as a bare metal atom beside loose ligands and counter-ions. Wilkinson's catalyst, for example, is "[Cl-].[Rh].PPh3.PPh3.PPh3". Drawn as stored, that was a free Rh beside floating ligands, and in 3D the pieces landed anywhere (the metal sometimes on a ligand atom). For the drawing and the 3D model only, the complex is now assembled by general rules. The analysis, formula and charges still use the stored record.
- **Pieces:**
  - **Counter-ions:** K⁺, Na⁺, H⁺, NH₄⁺, PF₆⁻, BF₄⁻, SbF₆⁻, ClO₄⁻ and BPh₄⁻ are never ligands. Water and solvent also stay separate.
  - **σ donors:** P/As, CO and isocyanide carbon, cyanide, pyridine/amine/imine N, halides, O⁻ of carboxylates, alkoxides and sulfonates, thiolates, and hydride.
  - **Chelates:** several donors on one piece that can close a 5- or 6-membered chelate ring, such as bipyridine, salen, oxalate and diamines.
  - **π ligands:** η⁵-cyclopentadienyl and η² isolated alkenes, such as ethylene and COD.
- **Binding:** chelates and Cp first, then neutral donors, then anions, while the metal's valence electron count stays within 18. The cap is 16 for d8 Ni/Pd/Pt/Cu/Ag/Au, which are square planar. In Ru(bpy)₃Cl₂ the bipyridines fill Ru(II) to 18 electrons, so the chlorides stay counter-ions.
- **Oxidation state** (ionic model): the metal's charge, plus ligands already bonded to it, plus bound radical ligands, plus any negative charge the record leaves unbalanced. PubChem writes Rh(I) as neutral Rh beside Cl⁻.
- **Disconnected amines:** PubChem's metal disconnection writes coordinated amines as "[NH⁻]" beside "[Pt+4]" (oxaliplatin). These are rebuilt as neutral amine donors, so oxaliplatin is square-planar Pt(II).
- **Bonds:** neutral donors are dative (arrow to the metal). Anions bond covalently as far as the metal's charge allows, and the rest stay dative with their charge (Zeise's [PtCl₃(C₂H₄)]⁻). Bound CO is drawn M–C≡O without charges.
- **Geometry**, from the number of positions (a Cp ring or alkene counts as one) and the d-electron count:
  - 2 positions: linear;
  - 3 positions: trigonal planar;
  - 4 positions: square planar for d8/d9, else tetrahedral;
  - 5 positions: trigonal bipyramidal for d8/d10, else square pyramidal (Grubbs, Jacobsen's);
  - 6 positions: octahedral.

  Complexes stored already bonded (Grubbs' catalyst, [PtCl₆]²⁻) get the same geometry. There, a neutral P/N/O/S donor with its full valence counts as a dative ligand even when written with a single bond.
- **2D:**
  - The metal is the hub. Each ligand is laid out on its own and set in evenly spaced slots; a chelate takes as many neighbouring slots as it has donors and points along their middle.
  - Ligands are pushed outward when they would collide. Metal–ligand bonds may then be drawn long, as in textbooks; three bipyridines do not fit around Ru at normal length.
  - Small ligands (Cl, CO, NH₃) move to the most open direction.
  - π ligands get a bond from the metal to the ring or alkene centre: ferrocene is drawn as a sandwich with the Cp rings around their centres, and an alkene lies side-on. The centre bonds are not counted as crossings.
- **3D:**
  - The ligand pieces are embedded, then set in the ideal geometry. The template is rotated onto the ligands' current directions. Each ligand moves rigidly: its donor at a covalent-radius distance, and its bulk pointing away from the metal (carbonyls Fe–C≡O, not Fe–O).
  - Chelates are placed by a rigid fit. Multidentate ones (salen) are embedded bonded to the metal without enforced chirality, and any conformer whose stereocentres do not match the input is dropped.
  - Cp rings sit perpendicular to the axis (Fe–ring 1.7 Å, Zr–ring 2.15 Å), and alkenes side-on.
  - Monodentate ligands are spun about their bond, over a few passes, to keep clear of each other.
  - Counter-ions and solvent are set clear of the complex.
  - The chain-straightening and tail-swing steps are skipped for complexes.
- **Tested:**
  - Wilkinson's, Vaska's and Crabtree's catalysts; Pd(PPh₃)₄, PdCl₂(PPh₃)₂ and Pt(PPh₃)₄; ferrocene and Cp₂ZrCl₂.
  - Zeise's salt; cisplatin, carboplatin and oxaliplatin; Fe(CO)₅, Ni(CO)₄, Mo(CO)₆ and Cr(CO)₆.
  - Ru(bpy)₃Cl₂, [Co(NH₃)₆]Cl₃, Jacobsen's catalyst, Ti(OiPr)₄, Cu(acac)₂, Grubbs I and Speier's acid.
  - The 3D models have their ideal angles: 90/180° square planar, 109.5° tetrahedral, trigonal bipyramidal, square pyramidal and octahedral.
- **Not handled:**
  - Records with more than one metal (Pd₂(dba)₃, Tebbe's and Stryker's reagents) are drawn as stored, and the result says so.
  - η⁶-arene complexes.
  - Bridging ligands; Pd(OAc)₂ is drawn as a κ¹ monomer, not the real trimer.
  - Spin-state exceptions the rules cannot see: NiCl₂(PPh₃)₂ is really tetrahedral, but the rules make it square planar.
  - Real bond lengths are approximated from covalent radii.

**Lone pairs, radicals and VSEPR in 3D.** Lone pairs used to be placed from RDKit's hybridisation label, which can disagree with the shape actually built. RDKit calls a Wittig ylide's CH₂⁻ sp3 while MMFF builds it flat, so "opposite the three bonds" was a near-zero vector, and the pair pointed 14° from a C–H bond.
- **Placement:** now read from the geometry. A pyramidal centre gets its pair opposite the bonds; a flat one gets it in the p orbital, along the normal. Every result is checked: no lone pair or radical may lie within 80° of a bond, otherwise the directions are re-spread by repulsion from the real bonds. A 30-species panel, from ylides and carbenes to TEMPO, NO₂ and trityl, is tested.
- **Restraints:** MMFF has no parameters for radicals, carbenes or localised carbanions, so those centres are held at their VSEPR angles:
  - **localised carbanions:** pyramidal. CH₃⁻ and ammonium/sulfonium ylide carbons were built flat. Conjugated ones (benzyl, enolate, Cp⁻) and phosphonium ylide carbons stay planar.
  - **σ radicals** (vinyl): 135°; MMFF built them linear.
  - **carbenes:** about 105° for a singlet with a heteroatom neighbour (CCl₂), about 136° for a triplet otherwise (CH₂). SMILES does not give the spin state, and the model says which was assumed.

  Methyl and tBu radicals keep the force field's planar-to-slightly-pyramidal shape.
- **One-heavy-atom molecules:** water, methane and the methyl radical had every hydrogen collapsed onto the central atom in 3D. Aligning a single heavy atom to the drawing gave a scale of zero; it now keeps the model's own size.
- **Radical dots in 2D:** with lone pairs shown, a radical's electron was drawn twice, RDKit's dot plus the overlay's ("·CH₃·"). RDKit's dot is now hidden then, and shown alone when lone pairs are off (browser-tested).
- **Charged carbon labels:** RDKit writes a charged carbon as "C⁻" without its hydrogens. Such carbons are now labelled "CH₂⁻", facing either way.

**Crowded layouts** (`moltalk/depiction.py`). The search runs when RDKit's layout has crossings, overlapping atoms, or a bond longer than 1.5× the median. Before this change, BINAP got CoordGen's layout with its biaryl bond stretched 3.8×.
- **Candidates:** RDKit's own layout, RDKit with random ring flips, CoordGen at best precision, and a hub layout. The hub layout sets each branch around the most central branching atom on its own.
- **Repair:** each candidate is repaired by a greedy search over rigid moves. For each acyclic bond, the smaller side is set to its target length, then turned or mirrored about the bond. Only the moved part is rescored.
- **Score:** the energy terms are bond-length deviation, crowding, an atom lying on a bond, bond crossings, uneven bond angles around an atom, and ring polygons that are not regular. A 3.8× bond costs about 730, so a stretched layout never wins. Metal–ligand bonds may be drawn up to 2.5× long.
- **Acceptance:** a candidate never has more crossings or overlapping atoms than the current layout. Its wedges must read back (molfile round trip) as the input's stereochemistry. It is then picked by crossings and overlaps, stretched bonds, and energy.
- **Panel:** BINAP, BINOL, gossypol, Xantphos, rubrene, hexaphenylbenzene, tetraphenylmethane, PPh₃, Wilkinson's catalyst, Pd(PPh₃)₄ (both as stored), and [5]/[6]helicene all draw with no crossings, no overlaps, and a longest bond of 1.5× or less. The stereo-unspecified BINAP and both atropisomers get the same clean layout.
- **500 library compounds, before → after:**
  - layouts with a bond over 1.5×: 6 → 3;
  - longest bond: 2.6× (unchanged); median longest bond: 1.0×;
  - crossings: 15 → 15; overlapping atoms: 2 → 1;
  - layouts changed: 12, none with more crossings or overlaps.

  The three left are two bicyclic amanitin-type peptides and a bridged tropane. A transannular bridge cannot be drawn flat without either a long bond or a crossing, and MolTalk keeps the long bond rather than add a crossing.
- **Limitation:** covalently bonded Pd(PPh₃)₄ (four P–Pd bonds written out) cannot be drawn without crossings even with long M–P bonds, because four PPh₃ groups do not fit around one atom in 2D. Textbooks write "PPh₃" for this; the library stores it as separate ligands, which draw cleanly.
- **Cost:** layout time over the 500 compounds went from 10 to about 22 s (median unchanged). The search has a 1.5 s budget, and a third of that above 60 atoms, where CoordGen alone already takes seconds.

**Slow CPUs.** If 8 conformers cannot be embedded in time, as for F430 on Cloud Run's single slower CPU, one conformer is embedded with a 10 s budget before giving up. Molecules over 50 heavy atoms start with 4 candidates instead of 8. Small molecules (up to 30 heavy atoms) whose best conformer still fits the flat drawing poorly try 32. With only 8 random conformers, 2-bromobutane had no anti chain, so its methyl swung about 2 bond lengths on lifting.

**Big molecules.** If 32 conformers time out, the best of the first 8 is kept, and the 32-conformer retry is skipped when the first round took over 1.5 s. Erythromycin gets its 3D model in about 4 s instead of failing.

**Stereocentre counting.** RDKit's list of *potential* stereocentres is generous; for example, it lists adamantane's bridgeheads. An unspecified centre now counts only if flipping it actually changes the molecule (`moltalk/stereo.py`), so adamantane reports no stereo elements, while cis/trans-1,4-dimethylcyclohexane still reports two.

**Hydrogens.** With **Show hydrogens**, the heavy-atom skeleton is laid out first (with every fix in this section), then the hydrogens are placed around it, so both drawings share one skeleton. RDKit's all-at-once layout is used instead only when it has strictly fewer crossings, overlaps and stretched bonds (glucose, morphine). Laying out everything at once had left FAD crowded enough to trigger the cage fallback.

**Drawing quality check.** Every layout is scored for bond crossings, stretched bonds and overlapping atoms (`moltalk/depiction.py`). RDKit's flat layout fails on polyhedral cages: for C₆₀ it produced 37 crossing bond pairs and bonds up to 16× normal length. Cages (atoms shared by three or more rings) are therefore drawn the way textbooks draw them: as a view of their 3D shape. The 3D model is the same one the rotated view uses. Of 120 viewing directions, the one with the fewest coinciding atoms, then fewest crossing bonds, is kept. Cubane comes out as a cube in perspective and adamantane as its familiar cage. Shown hydrogens are projected from the same model, so they point outward; in the earlier Schlegel diagrams, inner atoms' hydrogens piled up in the middle. A substituted fullerene (PCBM) projects the cage, and RDKit places the side chain around it. Crossings are expected in such a view, so it carries a note instead of a layout warning, and RDKit's red close-contact boxes are turned off. Ordinary molecules keep RDKit's layout. Tested on cubane, adamantane, dodecahedrane, C₆₀ and PCBM, with and without hydrogens.

**Server instructions** tell the model to resolve names through `resolve_name`, check CIP labels against stereo descriptors in a name, explain wedges only from `depicted_stereo_bonds`, reuse `canonical_smiles` in follow-ups, and never substitute a structure after an error.

## 4. Security and limits

- **Tunnel deployment:** the MCP server is a stdio child of tunnel-client and opens no ports. tunnel-client makes outbound HTTPS connections only, to `api.openai.com`, authenticated with the runtime key read from a mode-600 file. Its admin UI is bound to 127.0.0.1. The systemd unit sets `NoNewPrivileges`, `MemoryMax=3G` and `TasksMax=256`.
- **RDKit work** runs in two worker processes with a 20 s per-call timeout (the worker is killed when it is exceeded), a 2 GB address-space cap, at most 16 queued requests, and stdout redirected so a worker can never corrupt the stdio stream. Configure with `MOLTALK_TIMEOUT_S`, `MOLTALK_WORKERS`, `MOLTALK_WORKER_MEMORY_MB` and `MOLTALK_MAX_QUEUED`.
- **Input limits:** 4096 SMILES characters, 256 atoms, 1024 SMARTS characters, 100 substructure matches, 64 enumerated isomers (16 drawn), drawings 200–1600 px, names ≤256 characters.
- **HTTP mode** (local testing and Docker): Host/Origin validation (DNS-rebinding protection) is **always on**, including when binding 0.0.0.0 in Docker. Only `localhost`/`127.0.0.1`/`[::1]` are accepted unless `MOLTALK_ALLOWED_HOSTS`/`MOLTALK_ALLOWED_ORIGINS` are set. Request bodies are limited to 64 KB (`MOLTALK_MAX_BODY_BYTES`). There is no authentication, so keep it on loopback.
- **PubChem:** `resolve_name` contacts PubChem only when the name is not in the library, OPSIN cannot read it, and the model sets `allow_network=true`; the compound name is then sent to PubChem. Drawings send the InChIKey of a structure that is not in the library (turn off with `MOLTALK_OFFLINE=1`). Nothing else leaves your computer except tunnel traffic to OpenAI.

## 5. Local HTTP and Docker (testing only)

```bash
.venv/bin/moltalk --transport streamable-http          # http://127.0.0.1:8000/mcp
npx @modelcontextprotocol/inspector@latest                 # point it at the URL above

docker build -t moltalk .
docker run --rm -p 127.0.0.1:8000:8000 --read-only --tmpfs /tmp --memory 2g --cpus 2 --pids-limit 128 moltalk
```

### Public deployment (Google Cloud Run)

The public instance is at `https://moltalk-411294000488.us-central1.run.app/mcp` (project `moltalk`, region `us-central1`). Deploy from this directory:

```bash
gcloud run deploy moltalk --source . --region us-central1 --allow-unauthenticated \
  --min-instances 0 --max-instances 2 --cpu 1 --memory 2Gi --concurrency 8 --timeout 60 --cpu-boost \
  --set-env-vars "^|^MOLTALK_ALLOWED_HOSTS=moltalk-411294000488.us-central1.run.app,moltalk-jo5jxgfila-uc.a.run.app,localhost:*,127.0.0.1:*"
```

- **Cost:** it scales to zero, so an idle service costs nothing, and normal use stays inside Cloud Run's free tier. `--max-instances 2` caps the bill.
- **Protection** (`moltalk/public.py`): limits of 60 requests per minute per ChatGPT user (`openai/subject`, or the client IP) and 600 per minute in total. `/health` is the health check (Cloud Run reserves `/healthz`). `/.well-known/openai-apps-challenge` serves `MOLTALK_OPENAI_CHALLENGE` for OpenAI's domain verification.
- **Names without PubChem:** the library and OPSIN run inside the container. PubChem is rate-limited to 4 requests per second, and is only a fallback.
- **Cloud NAT and a static IP (documented, not provisioned).** If real users repeatedly hit common names that are not in the library, and PubChem keeps answering "server busy" from Cloud Run's shared addresses, route outgoing traffic through a reserved IP. That needs a VPC connector or Direct VPC egress, a Cloud Router, Cloud NAT and a static address. It costs roughly US$4–5 a month while provisioned, mainly the reserved IP and the NAT gateway hours. A cheaper first step is to add the missing names to `scripts/library_seed_names.txt` and rebuild.

Note on the original V1 Dockerfile: it could never start, because `python:3.12-slim` lacks `libXrender`/`libXext`/`libexpat`, which RDKit's drawing module needs. That is fixed here.

## 6. Troubleshooting

- `scripts/status.sh` shows *not ready*: check `journalctl --user -u moltalk-tunnel`. "control plane API key" errors mean the key file is missing or wrong; rerun `store-tunnel-key.sh`. Permission errors mean the key's principal lacks Tunnels Read + Use.
- The plugin shows no tools or fails to create: make sure the tunnel is running *while* you create it, and that the tunnel ID in the plugin matches `~/.config/moltalk/tunnel-profiles/moltalk.yaml`.
- The drawing doesn't appear but the answer does: refresh the plugin and start a new chat. The model still receives the full chemistry data.
- `scripts/local-tunnel-test.sh` runs the tunnel path locally without credentials and prints a local MCP URL for the MCP Inspector.

## 7. Chemistry scope

Indices are zero-based **input-SMILES** atom indices, not IUPAC locants and not canonical-SMILES order. Wedge/dash belongs to a particular 2D depiction (a wedge starts at the stereocenter and points toward the viewer). It is not a synonym for R/S. In the cyclohexane example, both the R and the S centre carry wedges.

Names are resolved as described under *Name resolution*. Library and PubChem records may leave stereochemistry unspecified. For example, "lactic acid" (CID 612) is unspecified and "glucose" (CID 5793) is partially specified, and plain "alanine" is L-alanine, as on PubChem (the returned title says so). That status is reported explicitly in `stereo_summary` and the warning. logP is a calculated estimate. Functional groups come from a small SMARTS motif dictionary (ester oxygens are not reported as ethers) plus aromatic-ring detection; this is not an exhaustive classification. Not supported: reaction prediction, conformer energies, chair conformations, structure editing, general IUPAC parsing.

Sources: [OpenAI plugins: connect and test](https://developers.openai.com/plugins/deploy/connect-chatgpt), [Add UI to your MCP server](https://developers.openai.com/plugins/build/chatgpt-ui), [Plugins reference](https://developers.openai.com/plugins/reference), [Secure MCP Tunnel](https://developers.openai.com/api/docs/guides/secure-mcp-tunnels), [Developer mode](https://developers.openai.com/api/docs/guides/developer-mode), [openai/tunnel-client](https://github.com/openai/tunnel-client), [RDKit](https://www.rdkit.org/docs/), [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk).

## Plugin package and logo

`plugin/` is the portable plugin package: `plugin.json` (branding under `extensions.com.openai.interface`: name, descriptions, `logo` and `composerIcon`), `mcp.json` (the Cloud Run MCP URL) and `assets/moltalk-logo.png` (512×512 PNG, the logo with thicker lines on a white rounded square, so it stays legible at 48 px and on dark themes). Build the upload with `cd plugin && zip -r ../dist/moltalk-plugin.zip .`. Logo and metadata change only when the package is re-uploaded; redeploying the server does not update them.

The server also advertises the icon itself: `serverInfo.icons` (a 128 px data URI, so it works over stdio and HTTP) and, on Cloud Run, `/logo.png` and `/favicon.ico`, for hosts that use the connector's site icon. The original transparent logo is `assets/logo.svg`; the README banner is `assets/banner.svg`.

## License

MolTalk is open source under the [MIT License](LICENSE).

It builds on:
- [RDKit](https://www.rdkit.org/) (BSD 3-Clause)
- the [MCP Python SDK](https://github.com/modelcontextprotocol/python-sdk) (MIT)

`scripts/fetch-tools.sh` downloads, but the repository does not include:
- [OPSIN](https://github.com/dan2097/opsin) (MIT)
- OpenAI's [tunnel-client](https://github.com/openai/tunnel-client), under its own licence

Chemical names and structures come from [PubChem](https://pubchem.ncbi.nlm.nih.gov/) (NCBI; see their [policies](https://www.ncbi.nlm.nih.gov/home/about/policies/)) and [Wikidata](https://www.wikidata.org/) (CC0). The bundled library `moltalk/data/compounds.json.gz` is derived from both and checked with RDKit and OPSIN.

