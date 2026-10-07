from importlib.resources import files
from typing import Any
import argparse
import time
import json
import hashlib
import os
import httpx
from mcp.server.fastmcp import Context, FastMCP
from mcp.server.transport_security import TransportSecuritySettings
from mcp.types import Annotations, Icon, CallToolResult, EmbeddedResource, TextContent, TextResourceContents, ToolAnnotations
from .chemistry import analyze, draw, conformer, substructure, enumerate_stereo
from .export import export_structure as write_structure, MAX_EMBED_BYTES
from .limits import runner
from .naming import name_and_locants, pubchem_slot

WIDGET_URI = "ui://widget/molecule-v27.html"
WIDGET_MIME = "text/html;profile=mcp-app"
WIDGET_HTML = files("moltalk").joinpath("widget/molecule.html").read_text(encoding="utf-8")

INSTRUCTIONS = """RDKit chemistry tools. Workflow rules:
- To draw a compound by name, call draw_named_molecule once: it resolves and draws in one call. To look a name up without drawing, call resolve_name and use its canonical_smiles. It works offline for thousands of common compounds and any systematic IUPAC name; allow_network=true additionally lets it try PubChem. Report the source, PubChem CID/IUPAC name and stereo_summary; say explicitly when stereochemistry is missing or partial.
- If you write SMILES yourself from a stereo-specific name, check the returned CIP labels against the name's descriptors and say whether they match.
- To show a structure from SMILES, call draw_molecule; the drawing renders inline for the user, and you receive the analysis plus depicted_stereo_bonds. Each call opens a new viewer, so draw a molecule once per reply: the viewer itself has controls for hydrogens, atom numbers, stereo labels, lone pairs, 3D rotation and image export, so never redraw just to change those.
- Atom indices are zero-based input-SMILES indices, not IUPAC locants. When draw_molecule returns locants (from the verified IUPAC name), refer to atoms by locant (e.g. C6a) and use indices only internally. Never invent locants or names the tools did not return.
- Explain wedges/dashes only from depicted_stereo_bonds of that drawing; wedge/dash is not a synonym for R/S.
- If draw_molecule's depiction has method 'projection' (a cage drawn as a view of its 3D shape) or a warning, tell the user what that means for the picture before describing it.
- For follow-up questions about the same molecule, reuse the canonical_smiles from the earlier result instead of re-deriving it.
- If a tool returns an error, report it; never substitute or invent a different structure.
- If resolve_name cannot resolve a name (unknown, ambiguous, or stereo not stated), do not write a SMILES for it from memory. Tell the user what the tool said and ask for a structure, a SMILES or a more specific name. Only if the user then asks you to proceed from your own knowledge may you write the SMILES; say plainly that it is unverified by MolTalk.
- When the user wants a structure file (ChemDraw, MOL, SDF, PDB, XYZ, SMILES), call export_structure with the verified canonical_smiles; default format cdxml for ChemDraw. 2D files use the drawing's layout; 3D files (mol/sdf with coordinates="3d", pdb, xyz) use one calculated conformer. Relay its warnings: a 3D file fixes an arbitrary configuration at any stereocentre the input leaves unspecified. The file is offered to the user by the inline viewer's Download button; do not paste the file contents unless asked.
- When you pass label= to draw_molecule, use the compound's name only if the SMILES came from resolve_name for that name (or the user gave it). draw_molecule checks the label against MolTalk's library: if label_check.status is "mismatch", the drawing is NOT that compound; say so and do not present it under that name."""

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

# The MolTalk logo (a dodecahedrane drawn by MolTalk's own 3D view), advertised in serverInfo.icons so hosts can
# show it next to the app's name. Embedded as a data URI, so it works over stdio and HTTP alike.
_ICON_128 = "data:image/png;base64," + __import__("base64").b64encode(files("moltalk").joinpath("static/icon-128.png").read_bytes()).decode()
ICONS = [Icon(src=_ICON_128, mimeType="image/png", sizes=["128x128"])]

mcp = FastMCP("MolTalk", instructions=INSTRUCTIONS, icons=ICONS, website_url="https://github.com/DomFico/moltalk", host=HOST, port=PORT,
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
                        numbering: str = "none", atom_indices: bool = True, hydrogens: bool = False,
                        include_svg: bool = False, ctx: Context | None = None) -> CallToolResult:
    """Draw and render a molecule inline for the user from SMILES, in a single call. A successful result (rendered: true) means the drawing is already displayed to the user: do not call draw_molecule again with the same arguments to verify, inspect, retrieve or display it; read structuredContent instead. For a compound name, use draw_named_molecule. Returns the full analysis and depicted_stereo_bonds (this drawing's wedges/dashes, each with an explanation). label: display caption, e.g. the compound name. numbering (default "none": a clean drawing; the user can switch numbers on in the viewer): "iupac" (parent-chain/ring locants from the verified IUPAC name, when they can be determined without ambiguity; otherwise atom indices), "indices" or "none". Ask for numbers only when the discussion needs them. hydrogens=true draws every hydrogen explicitly. include_svg=true also returns the raw SVG text (only for clients without the inline viewer)."""
    if numbering not in ("iupac", "indices", "none"):
        raise ValueError('numbering must be "iupac", "indices" or "none".')
    if not atom_indices:
        numbering = "none"
    checked = await runner.run(analyze, smiles)  # validate first, so an invalid structure is never sent to PubChem
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
    repeated = _repeated_draw(ctx, [smiles, label, width, height, numbering, atom_indices, hydrogens, include_svg])
    structured = {"kind": "molecule", "rendered": True, "label": label[:120] if label else None, "input_smiles": smiles,
                  "numbering_requested": numbering, "numbering_shown": shown, "atom_indices_shown": shown != "none",
                  "hydrogens_shown": hydrogens, "label_check": _check_label(label, checked["inchikey"]),
                  **naming, **result}
    if include_svg:
        structured["svg"] = svg
    if repeated:
        structured["repeat_of_previous"] = True  # the viewer shows a one-line note instead of a second copy
    centers = ", ".join((f"C{c['locant']} (atom {c['atom_index']})" if c.get("locant") else f"atom {c['atom_index']}") + f" {c['cip']}"
                        for c in a["stereocenters"]) or "none assigned"
    wedges = " ".join(b["explanation"] for b in result["depicted_stereo_bonds"]) or "no wedge/dash bonds in this drawing."
    if naming["iupac_name"]:
        name_text = f"IUPAC name: {naming['iupac_name']} ({naming['name_source']}, exact structure match). "
        name_text += ("Drawing numbered with its parent-structure locants. " if shown == "iupac"
                      else f"Locants not shown: {naming['locant_status']}. " if naming["locant_status"] else "")
    else:
        name_text = f"No IUPAC name: {naming['name_status']}. "
    check = structured["label_check"]
    if check and check["status"] == "mismatch":
        name_text = f"WARNING: {check['message']} " + name_text
    elif check and check["status"] == "stereo differs":
        name_text = f"Note: {check['message']} " + name_text
    text = (f"Rendered {label or a['canonical_smiles']} inline for the user (already displayed; do not call again to show it). "
            f"Drew {label or a['canonical_smiles']} ({a['formula']}); canonical SMILES {a['canonical_smiles']}. " + name_text +
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

def _check_label(label: str | None, inchikey: str) -> dict | None:
    """Is the caption the name of the structure drawn? Checked against MolTalk's library (connectivity, then stereo),
    so a model-written SMILES shown under a compound's name is caught. None when there is no label."""
    from . import library
    if not label:
        return None
    status, found = library.find_name(label)
    candidates = [found] if status == "hit" else (found or [])
    if not candidates:
        return {"status": "unverified",
                "message": f"'{label[:120]}' is not a name in MolTalk's library, so the label was not checked against the structure."}
    same = [c for c in candidates if c["inchikey"][:14] == inchikey[:14]]
    if not same:
        c = candidates[0]
        return {"status": "mismatch", "library_cid": c["cid"], "library_smiles": c["smiles"],
                "message": f"The label '{label[:120]}' names {c['title']} (PubChem CID {c['cid']}, SMILES {c['smiles']}) "
                           "in MolTalk's library, but the structure drawn is a different compound."}
    if not any(c["inchikey"][:23] == inchikey[:23] for c in same):  # first two blocks: connectivity + stereo/isotopes
        c = same[0]
        return {"status": "stereo differs", "library_cid": c["cid"], "library_smiles": c["smiles"],
                "message": f"Same connectivity as {c['title']} (PubChem CID {c['cid']}), but different or unspecified stereochemistry."}
    return {"status": "matches", "library_cid": same[0]["cid"],
            "message": f"The structure matches {same[0]['title']} (PubChem CID {same[0]['cid']})."}

@mcp.tool(annotations=ToolAnnotations(readOnlyHint=True, destructiveHint=False, openWorldHint=True),
          meta=_ui_meta("Drawing molecule…", "Molecule drawn"))
async def draw_named_molecule(name: str, numbering: str = "none", hydrogens: bool = False, allow_network: bool = False,
                              ctx: Context | None = None) -> CallToolResult:
    """Use this when the user asks to draw a molecule by name ("draw cubane", "show me FAD"): it resolves the name and renders the drawing inline in one call. Do not call resolve_name first, and do not call draw_molecule afterwards unless the user asks for a modified structure. A successful result (rendered: true) is already displayed; never call again to verify it. Resolution is the same as resolve_name (MolTalk's library, then OPSIN for systematic names, then PubChem only if allow_network=true); an ambiguous or unknown name, or one that does not say which stereoisomer, is an error and nothing is drawn. numbering and hydrogens are as in draw_molecule."""
    resolved = await resolve_name(name, allow_network)
    result = await draw_molecule(resolved["canonical_smiles"], label=name.strip()[:120], numbering=numbering,
                                 hydrogens=hydrogens, ctx=ctx)
    identity = {k: resolved[k] for k in ("source", "cid", "title", "iupac_name", "pubchem_url", "note", "warning")
                if resolved.get(k) is not None}
    result.structuredContent["resolved"] = identity
    said = f"Resolved '{name}' with {resolved['source']}" + (f" (PubChem CID {resolved['cid']}, {resolved.get('title')})" if resolved.get("cid") else "") + ". "
    extra = " ".join(x for x in (resolved.get("note"), resolved.get("warning")) if x)
    result.content[0].text = said + (extra + " " if extra else "") + result.content[0].text
    return result

_last_draw: dict[str, tuple[str, float]] = {}


def _repeated_draw(ctx, args) -> bool:
    """True when this ChatGPT session's previous draw_molecule had identical arguments, within 20 s. ChatGPT sometimes
    re-sends the first tool call of a chat although the first one succeeded, which showed two identical viewers.
    Only the immediately previous draw counts, so toggling a viewer setting back and forth is never collapsed."""
    meta = getattr(getattr(ctx, "request_context", None), "meta", None) if ctx is not None else None
    session = ((getattr(meta, "model_extra", None) or {}).get("openai/session")) if meta is not None else None
    if not session:
        return False
    digest = hashlib.sha256(json.dumps(args, sort_keys=True, default=str).encode()).hexdigest()
    now = time.monotonic()
    previous = _last_draw.get(str(session))
    _last_draw[str(session)] = (digest, now)
    if len(_last_draw) > 5000:
        for key in [k for k, (_, t) in _last_draw.items() if now - t > 60]:
            _last_draw.pop(key, None)
    return previous is not None and previous[0] == digest and now - previous[1] < 20

@mcp.tool(annotations=READ_ONLY, meta={"ui": {"visibility": ["app"]}, "openai/widgetAccessible": True,
                                       "openai/visibility": "private"})
async def conformer_3d(smiles: str, hydrogens: bool = False) -> dict[str, Any]:
    """Viewer-only: one calculated 3D conformer aligned to the flat drawing, used to rotate the drawing in 3D."""
    return await runner.run(conformer, smiles, hydrogens)

@mcp.tool(annotations=READ_ONLY, meta=_ui_meta("Preparing structure file…", "Structure file ready"))
async def export_structure(smiles: str, format: str = "cdxml", coordinates: str | None = None, name: str | None = None,
                           hydrogens: bool = False) -> CallToolResult:
    """Write the structure as a file the user can download: format "cdxml" (ChemDraw, default), "mol", "sdf", "pdb", "xyz" or "smiles". coordinates: "2d" (the drawing's layout; default for cdxml/mol/sdf) or "3d" (one calculated conformer with explicit hydrogens; required for pdb/xyz). name: compound name for the file name and title. hydrogens=true adds explicit hydrogens to a 2D file. The file holds the structure as given (charges, stereo, atom order), not drawing-only additions. Shown inline with a Download button."""
    result = await runner.run(write_structure, smiles, format, coordinates, name, hydrogens)
    text = result.pop("text")
    content = [TextContent(type="text", text=(
        f"Prepared {result['filename']} ({result['format_name']}, {result['coordinates']} coordinates, "
        f"{result['size_bytes']} bytes) for {result['canonical_smiles']}. The user can download it from the viewer. "
        + " ".join(result["warnings"] + result["notes"])).strip())]
    if result["size_bytes"] <= MAX_EMBED_BYTES:
        # The standard MCP way to return a file, for hosts without the viewer.
        content.append(EmbeddedResource(type="resource", annotations=Annotations(audience=["user"]),
                                        resource=TextResourceContents(uri=f"file:///{result['filename']}",
                                                                      mimeType=result["mime_type"], text=text)))
    return CallToolResult(content=content, structuredContent=result, _meta={"file_text": text})

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
    """Resolve a compound name to a structure. Tried in order: (1) MolTalk's bundled library of thousands of common compounds (trivial names, synonyms, IUPAC names; offline); (2) OPSIN, which reads systematic IUPAC names offline, including stereodescriptors; (3) PubChem, only when allow_network=true (best effort; the name is sent to PubChem). A name shared by different structures is reported as ambiguous, not guessed. Check the returned identity and stereo_summary before relying on it."""
    from . import library
    from .naming import opsin_smiles
    name = name.strip()
    if not name or len(name) > 256:
        raise ValueError("Name must contain 1–256 characters.")
    status, found = library.find_name(name)
    if status == "ambiguous":
        options = "; ".join(f"{c['title']} (PubChem CID {c['cid']}, SMILES {c['smiles']})" for c in found[:6])
        raise ValueError(f"'{name}' names more than one structure in MolTalk's library: {options}. "
                         "No structure was chosen; ask which one is meant, or use a SMILES or a stereo-specific name.")
    if status == "stereo":
        options = "; ".join(f"{c['title']} (PubChem CID {c['cid']}, SMILES {c['smiles']})" for c in found[:6])
        raise ValueError(f"'{name}' does not say which stereoisomer is meant (MolTalk's library: {options}). "
                         "No stereoisomer was chosen; ask which one is meant, then resolve that name (e.g. "
                         + " or ".join(f"'{c['title']}'" for c in found[:3]) + ").")
    if status == "hit":
        analysis = await runner.run(analyze, found["smiles"])
        result = {"source": "MolTalk library (from PubChem, names from Wikidata/PubChem)", "query": name,
                  "cid": found["cid"], "title": found["title"], "iupac_name": found.get("iupac"),
                  "iupac_verified": found.get("iupac_verified", False),
                  "pubchem_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{found['cid']}",
                  "canonical_smiles": analysis["canonical_smiles"], "stereo_summary": analysis["stereo_summary"],
                  "analysis": analysis}
        if found.get("note"):
            result["note"] = found["note"]  # curated: e.g. heme b vs hemin, protonation states of cofactors
        if analysis["stereo_summary"]["unspecified"]:
            result["warning"] = (f"This record leaves {analysis['stereo_summary']['unspecified']} stereo element(s) "
                                 "unspecified; the name does not identify a single stereoisomer.")
        return result
    smiles = await opsin_smiles(name)
    if smiles:
        analysis = await runner.run(analyze, smiles)
        known = library.by_inchikey(analysis["inchikey"])
        result = {"source": "OPSIN (systematic name parsed locally)", "query": name,
                  "canonical_smiles": analysis["canonical_smiles"], "stereo_summary": analysis["stereo_summary"],
                  "analysis": analysis,
                  "warning": "Structure derived from the name's systematic nomenclature; check it matches what was meant."}
        if known:
            result.update(cid=known["cid"], title=known["title"], iupac_name=known.get("iupac"))
        if analysis["stereo_summary"]["unspecified"]:
            result["warning"] += (f" The name leaves {analysis['stereo_summary']['unspecified']} stereo element(s) "
                                  "unspecified.")
        return result
    if not allow_network:
        raise ValueError(f"'{name}' is not in MolTalk's library and is not a systematic name OPSIN can read. "
                         "No structure was assumed. Supply a SMILES, or retry with allow_network=true to ask PubChem.")
    from urllib.parse import quote
    url = 'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/' + quote(name, safe='') + '/property/SMILES,ConnectivitySMILES,IUPACName,Title/JSON'
    try:
        await pubchem_slot()
        async with httpx.AsyncClient(timeout=15, follow_redirects=False) as client:
            response = await client.get(url)
    except httpx.HTTPError as exc:
        raise ValueError(f"PubChem could not be reached ({type(exc).__name__}). No structure was assumed.") from None
    if response.status_code == 404:
        raise ValueError(f"'{name}' was not found in MolTalk's library, OPSIN or PubChem. No structure was assumed; supply SMILES.")
    if response.status_code != 200:
        from .naming import _fault
        raise ValueError(f"PubChem lookup failed (HTTP {response.status_code}: {_fault(response)[:80]}).")
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
