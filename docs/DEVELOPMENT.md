# Development

## Setup

```bash
git clone https://github.com/DomFico/moltalk.git && cd moltalk
python3 -m venv .venv
.venv/bin/pip install -e '.[dev]'
scripts/fetch-tools.sh        # OPSIN 2.9.0 into vendor/ (SHA-256 verified); needs Java 11+ on the PATH
.venv/bin/pytest -q
```

Python 3.11 or later. `requirements-tested.txt` pins the versions the test suite last passed with (RDKit 2026.03,
MCP SDK 1.30).

## Running the server

```bash
.venv/bin/moltalk                                       # stdio, for local MCP clients
.venv/bin/moltalk --transport streamable-http           # http://127.0.0.1:8000/mcp
npx @modelcontextprotocol/inspector@latest              # inspect it
```

## Repository layout

```text
moltalk/            the package (server, chemistry, viewer in widget/, compound library in data/)
tests/              pytest suite, including real-browser tests of the viewer
docs/               this documentation
scripts/            library build, README demos, legacy tunnel scripts
plugin/             ChatGPT plugin package (manifest, MCP URL, logo)
assets/             logo, banner and README animations (assets/demo/)
examples/           example prompts and exported files
```

## Browser tests

The viewer tests drive the real widget in Google Chrome through Playwright, with a host harness that imitates
ChatGPT's MCP Apps bridge:

```bash
.venv/bin/pip install -e '.[dev,ui-test]'
.venv/bin/pytest -q tests/test_widget_*.py
```

`MOLTALK_UI_BROWSER_CHANNEL=chromium` uses Playwright's bundled Chromium instead of Chrome. `MOLTALK_UI_SCREENSHOTS=dir`
saves screenshots.

## README animations

The animations in `assets/demo/` are generated from the real pipeline: the script starts the MCP server, calls the
tools, shows the result in the real widget in headless Chrome, lifts it into 3D and records one seamless loop.

```bash
.venv/bin/pip install -e '.[dev,ui-test]'      # plus ffmpeg on the PATH
.venv/bin/python scripts/generate_readme_demos.py            # all
.venv/bin/python scripts/generate_readme_demos.py ferrocene  # one
```

It writes `<name>.gif` (used by the README) and `<name>.webp`; only the GIFs are committed. The widget exposes the
active viewer as `window.moltalkViewer` for this script; nothing else uses it.

## Rebuilding the compound library

```bash
.venv/bin/python scripts/build_library.py
```

Needs network access, Java and OPSIN (about 10 minutes). Sources and rules are in [DATA_SOURCES.md](DATA_SOURCES.md).
Add common names to `scripts/library_seed_names.txt`, and editorial corrections to `scripts/library_curated.json`.

## Changing the viewer

`moltalk/widget/molecule.html` is a single self-contained file. Hosts cache UI templates, so bump the version in
`WIDGET_URI` (`moltalk/server.py`) whenever it changes, then refresh the plugin or connector in the host.

## Branding assets

| File | Use |
|---|---|
| `assets/moltalk-logo.png` | Canonical 512 px app icon (identical to `plugin/assets/moltalk-logo.png` and `moltalk/static/icon-512.png`) |
| `moltalk/static/icon-128.png` | Small icon advertised in `serverInfo.icons` |
| `assets/moltalk-logo.svg` | Vector source of the app icon |
| `assets/logo.svg`, `assets/logo.png` | The original transparent logo (dodecahedrane) |
| `assets/banner.svg` | README banner |

The plugin package must keep its own copy under `plugin/assets/`, because the package is uploaded on its own.
