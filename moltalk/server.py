from importlib.resources import files
from typing import Any
import argparse
import os
import httpx
from mcp.server.fastmcp import FastMCP
from mcp.server.transport_security import TransportSecuritySettings
from mcp.types import CallToolResult, TextContent, ToolAnnotations
from .chemistry import analyze, draw, conformer, substructure, enumerate_stereo, NAMES
from .limits import runner
from .naming import name_and_locants, pubchem_slot

WIDGET_URI = "ui://widget/molecule-v15.html"
WIDGET_MIME = "text/html;profile=mcp-app"
WIDGET_HTML = files("moltalk").joinpath("widget/molecule.html").read_text(encoding="utf-8")

INSTRUCTIONS = """RDKit chemistry tools. Workflow rules:
- For a compound name, call resolve_name first (allow_network=true is needed for names outside the small offline dictionary) and use its canonical_smiles. Report the source, PubChem CID/IUPAC name and stereo_summary; say explicitly when stereochemistry is missing or partial.
- If you write SMILES yourself from a stereo-specific name, check the returned CIP labels against the name's descriptors and say whether they match.
- To show a structure, call draw_molecule; the drawing renders inline for the user, and you receive the analysis plus depicted_stereo_bonds.
- Atom indices are zero-based input-SMILES indices, not IUPAC locants. When draw_molecule returns locants (from the verified IUPAC name), refer to atoms by locant (e.g. C6a) and use indices only internally. Never invent locants or names the tools did not return.
- Explain wedges/dashes only from depicted_stereo_bonds of that drawing; wedge/dash is not a synonym for R/S.
- If draw_molecule's depiction has method 'schlegel' or a warning, tell the user what that means for the picture before describing it.
- For follow-up questions about the same molecule, reuse the canonical_smiles from the earlier result instead of re-deriving it.
- If a tool returns an error, report it; never substitute or invent a different structure."""

HOST = os.getenv("HOST", "127.0.0.1")
PORT = int(os.getenv("PORT", "8000"))
LOCAL_HOSTS = ["127.0.0.1:*", "localhost:*", "[::1]:*"]
LOCAL_ORIGINS = ["http://127.0.0.1:*", "http://localhost:*", "http://[::1]:*"]

def _env_list(name: str, default: list[str]) -> list[str]:
    value = os.getenv(name, "").strip()
    return [item.strip() for item in value.split(",") if item.strip()] if value else default

# Host/Origin validation stays on for every bind address, including 0.0.0.0 in Docker.
# A public deployment must add its own domain here explicitly (and sit behind OAuth).
SECURITY = TransportSecuritySettings(enable_dns_rebinding_protection=True,
                                     allowed_hosts=_env_list("MOLTALK_ALLOWED_HOSTS", LOCAL_HOSTS),
                                     allowed_origins=_env_list("MOLTALK_ALLOWED_ORIGINS", LOCAL_ORIGINS))

mcp = FastMCP("MolTalk", instructions=INSTRUCTIONS, host=HOST, port=PORT,
              stateless_http=True, json_response=True, transport_security=SECURITY,
              max_request_body_size=int(os.getenv("MOLTALK_MAX_BODY_BYTES", "65536")))

READ_ONLY = ToolAnnotations(readOnlyHint=True, destructiveHint=False, openWorldHint=False, idempotentHint=True)
def _ui_meta(invoking: str, invoked: str) -> dict[str, Any]:
    return {"ui": {"resourceUri": WIDGET_URI}, "openai/outputTemplate": WIDGET_URI,
            "openai/toolInvocation/invoking": invoking, "openai/toolInvocation/invoked": invoked}

@mcp.resource(WIDGET_URI, name="moltalk-viewer", title="MolTalk viewer", mime_type=WIDGET_MIME,
              description="Inline viewer for RDKit molecule drawings and stereoisomer grids.",
              meta={"ui": {"prefersBorder": True, "csp": {"connectDomains": [], "resourceDomains": []}},
                    "openai/ui": {"availableDisplayModes": ["inline", "fullscreen"]},
                    "openai/widgetDescription": "Shows the RDKit 2D drawing (with atom indices, R/S labels and wedge/dash bonds) "
                                                "or a grid of stereoisomers. The chemistry data is in the tool result."})
def molecule_widget() -> str:
    return WIDGET_HTML

@mcp.tool(annotations=READ_ONLY)
async def analyze_molecule(smiles: str) -> dict[str, Any]:
    """Validate SMILES; return identifiers, properties, CIP R/S, bond stereo, stereo_summary and SMARTS motifs. Indices are zero-based in the input molecule. Does not draw."""
    return await runner.run(analyze, smiles)

@mcp.tool(annotations=READ_ONLY, meta=_ui_meta("Drawing molecule…", "Molecule drawn"))
async def draw_molecule(smiles: str, label: str | None = None, width: int = 640, height: int = 420,
                        numbering: str = "iupac", atom_indices: bool = True, hydrogens: bool = False,
                        include_svg: bool = False) -> CallToolResult:
    """Draw a molecule inline for the user from SMILES (RDKit 2D depiction with atom indices and R/S labels). Returns the full analysis and depicted_stereo_bonds (this drawing's wedges/dashes, each with an explanation). label is an optional display caption, e.g. the compound name. numbering: "iupac" (parent-chain/ring locants from the verified IUPAC name, when they can be determined without ambiguity; otherwise atom indices), "indices" or "none". hydrogens=true draws every hydrogen explicitly. include_svg=true also returns the raw SVG text (only for clients without the inline viewer)."""
    if numbering not in ("iupac", "indices", "none"):
        raise ValueError('numbering must be "iupac", "indices" or "none".')
    if not atom_indices:
        numbering = "none"
    await runner.run(analyze, smiles)  # validate first, so an invalid structure is never sent to PubChem
    naming = await name_and_locants(smiles, runner.run)
    locants = naming["locants"] if numbering == "iupac" else None
    shown = "iupac" if locants else ("none" if numbering == "none" else "indices")
    result = await runner.run(draw, smiles, width, height, shown == "indices", hydrogens, locants)
    svg = result.pop("svg")
    atom_px = result.pop("atom_px")  # widget-only: where each atom sits in the flat drawing
    drawn_bonds = result.pop("drawn_bonds")
    a = result["analysis"]
    for centre in result["analysis"]["stereocenters"]:
        if naming["locants"].get(str(centre["atom_index"])):
            centre["locant"] = naming["locants"][str(centre["atom_index"])]
    structured = {"kind": "molecule", "label": label[:120] if label else None, "input_smiles": smiles,
                  "numbering_requested": numbering, "numbering_shown": shown, "atom_indices_shown": shown != "none",
                  "hydrogens_shown": hydrogens, **naming, **result}
    if include_svg:
        structured["svg"] = svg
    centers = ", ".join((f"C{c['locant']} (atom {c['atom_index']})" if c.get("locant") else f"atom {c['atom_index']}") + f" {c['cip']}"
                        for c in a["stereocenters"]) or "none assigned"
    wedges = " ".join(b["explanation"] for b in result["depicted_stereo_bonds"]) or "no wedge/dash bonds in this drawing."
    if naming["iupac_name"]:
        name_text = f"IUPAC name: {naming['iupac_name']} ({naming['name_source']}, exact structure match). "
        name_text += ("Drawing numbered with its parent-structure locants. " if shown == "iupac"
                      else f"Locants not shown: {naming['locant_status']}. " if naming["locant_status"] else "")
    else:
        name_text = f"No IUPAC name: {naming['name_status']}. "
    text = (f"Drew {label or a['canonical_smiles']} ({a['formula']}); canonical SMILES {a['canonical_smiles']}. " + name_text +
            f"CIP centers (zero-based input indices): {centers}. Stereo: {a['stereo_summary']['status']}. "
            f"Drawing: {wedges} Functional-group motifs: {', '.join(a['functional_groups']) or 'none matched'}.")
    depiction = result["depiction"]
    if depiction.get("note"):
        text += " " + depiction["note"]
    if depiction.get("warning"):
        text += " " + depiction["warning"]
    if depiction.get("coordination_note"):
        text += " " + depiction["coordination_note"]
    return CallToolResult(content=[TextContent(type="text", text=text)], structuredContent=structured,
                          _meta={"svg": svg, "width": width, "height": height, "atom_px": atom_px, "bonds": drawn_bonds})

@mcp.tool(annotations=READ_ONLY, meta={"ui": {"visibility": ["app"]}, "openai/widgetAccessible": True,
                                       "openai/visibility": "private"})
async def conformer_3d(smiles: str, hydrogens: bool = False) -> dict[str, Any]:
    """Viewer-only: one calculated 3D conformer aligned to the flat drawing, used to rotate the drawing in 3D."""
    return await runner.run(conformer, smiles, hydrogens)

@mcp.tool(annotations=READ_ONLY)
async def find_substructure(smiles: str, smarts: str) -> dict[str, Any]:
    """Find chirality-aware SMARTS matches; at most 100 matches."""
    return await runner.run(substructure, smiles, smarts)

@mcp.tool(annotations=READ_ONLY, meta=_ui_meta("Enumerating stereoisomers…", "Stereoisomers shown"))
async def enumerate_stereoisomers(smiles: str, limit: int = 16) -> CallToolResult:
    """Enumerate unspecified atom/bond stereochemistry (specified centers are preserved) and show the isomers inline as a grid. Each isomer has CIP labels (input atom indices), an achiral/meso flag and its enantiomer's index."""
    result = await runner.run(enumerate_stereo, smiles, limit, True)
    svgs = result.pop("svgs")
    isomers = result["isomers"]
    lines = []
    for iso in isomers:
        cip = ", ".join(f"atom {c['atom_index']} {c['cip']}" for c in iso["stereocenters"]) or "no CIP centers"
        relation = "achiral/meso" if iso["achiral"] else (f"enantiomer of #{iso['enantiomer_index']}"
                                                          if iso["enantiomer_index"] is not None else "chiral")
        lines.append(f"#{iso['index']} {iso['smiles']} ({cip}; {relation})")
    text = (f"{len(isomers)} stereoisomer(s) for {smiles}{' (truncated)' if result['truncated'] else ''}: "
            + "; ".join(lines) + ".")
    return CallToolResult(content=[TextContent(type="text", text=text)],
                          structuredContent={"kind": "stereoisomers", **result}, _meta={"svgs": svgs})

@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True, destructiveHint=False, openWorldHint=True))
async def resolve_name(name: str, allow_network: bool = False) -> dict[str, Any]:
    """Resolve a compound name to a structure. A small offline dictionary is checked first; otherwise allow_network=true sends the name to PubChem. RDKit is not a general IUPAC parser. Check the returned CID/IUPAC name and stereo_summary before relying on the identity."""
    name = name.strip()
    if not name or len(name) > 256:
        raise ValueError("Name must contain 1–256 characters.")
    if name.lower() in NAMES:
        analysis = await runner.run(analyze, NAMES[name.lower()])
        return {"source": "offline dictionary", "query": name, "canonical_smiles": analysis["canonical_smiles"],
                "stereo_summary": analysis["stereo_summary"], "analysis": analysis}
    if not allow_network:
        raise ValueError("Name is not in the offline dictionary. Supply SMILES or explicitly enable PubChem lookup.")
    from urllib.parse import quote
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + quote(name, safe='') + '/property/SMILES,ConnectivitySMILES,IUPACName,Title/JSON'
    try:
        await pubchem_slot()
        async with httpx.AsyncClient(timeout=15, follow_redirects=False) as client:
            response = await client.get(url)
    except httpx.HTTPError as exc:
        raise ValueError(f"PubChem could not be reached ({type(exc).__name__}). No structure was assumed.") from None
    if response.status_code == 404:
        raise ValueError(f"PubChem has no compound named {name!r}. No structure was assumed; supply SMILES.")
    if response.status_code != 200:
        raise ValueError(f"PubChem lookup failed (HTTP {response.status_code}).")
    records = response.json()['PropertyTable']['Properties']
    if len(records) != 1:
        raise ValueError("Name resolves to multiple structures; use an explicit SMILES.")
    row = records[0]
    smiles = row.get('SMILES') or row.get('IsomericSMILES')
    if not smiles:
        raise ValueError("PubChem did not return an isomeric SMILES.")
    analysis = await runner.run(analyze, smiles)
    stereo = analysis["stereo_summary"]
    warning = "Name lookup is a database match, not proof that the intended structure was specified."
    if stereo["unspecified"]:
        warning += (f" PubChem's record leaves {stereo['unspecified']} stereo element(s) unspecified; "
                    "the name does not identify a single stereoisomer.")
    return {"source": "PubChem", "query": name, "cid": row['CID'], "title": row.get('Title'),
            "iupac_name": row.get('IUPACName'), "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{row['CID']}",
            "canonical_smiles": analysis["canonical_smiles"], "stereo_summary": stereo,
            "analysis": analysis, "warning": warning}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--transport', choices=['stdio', 'streamable-http'], default='stdio')
    args = parser.parse_args()
    if args.transport == 'stdio':
        mcp.run(transport='stdio')
        return
    # HTTP (Cloud Run, Docker): MCP app behind the public gateway (rate limits, /healthz, domain verification).
    # Request URLs contain compound names and InChIKeys, so keep HTTP-client logging to warnings.
    import logging
    import uvicorn
    from .public import PublicGateway
    logging.getLogger("httpx").setLevel(logging.WARNING)
    uvicorn.run(PublicGateway(mcp.streamable_http_app()), host=HOST, port=PORT, proxy_headers=True,
                forwarded_allow_ips="*", log_level="warning", timeout_graceful_shutdown=10)

if __name__ == '__main__':
    main()
