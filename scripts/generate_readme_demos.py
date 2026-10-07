"""Generate the README's rotating-molecule animations from MolTalk's real rendering pipeline.

Each demo starts the MolTalk MCP server (stdio), calls its tools exactly as an assistant would, shows the result in
the real widget inside a headless Chrome page (the same host harness the browser tests use), lifts the drawing into
3D and turns it one full revolution about the vertical axis in equal steps, so the loop is seamless. Frames are
cropped to the molecule's bounding box over the whole turn and encoded as GIF and animated WebP with ffmpeg.

    .venv/bin/pip install -e '.[dev,ui-test]'      # Playwright; uses Google Chrome
    .venv/bin/python scripts/generate_readme_demos.py            # all demos
    .venv/bin/python scripts/generate_readme_demos.py binap      # just one

Output: assets/demo/<name>.gif and .webp. Frames are written to a temporary directory and not kept.
Needs ffmpeg on the PATH.
"""
import asyncio
import json
import os
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from PIL import Image, ImageChops
from playwright.async_api import async_playwright

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "assets" / "demo"
FPS = 20
SECONDS = 6
WIDTH = 640  # output width cap in pixels
PAD = 18

# name: (tool, arguments, options). "isomer" picks an enumerate_stereoisomers result by its axial descriptor.
DEMOS = {
    "caffeine": ("draw_named_molecule", {"name": "caffeine"}, {"motion": "rock", "amplitude": 50}),
    "stereocentres": ("draw_molecule", {"smiles": "F[C@H]1C[C@@H](C)CCC1", "label": "(1R,3S)-1-fluoro-3-methylcyclohexane"},
                      {"motion": "rock", "amplitude": 55, "stereo": "Specified"}),
    "binap": ("draw_molecule", {"label": "(Ra)-BINAP"},
              {"isomer": ("c1ccc(cc1)P(c2ccccc2)c3ccc4ccccc4c3-c5c(ccc6ccccc56)P(c7ccccc7)c8ccccc8", "Ra"),
               "tilt": 25, "stereo": "Specified"}),
    "lone_pairs": ("draw_molecule", {"smiles": "C[N+](=O)[O-]", "label": "nitromethane"},
                   {"motion": "rock", "amplitude": 60, "lone_pairs": True}),
    "ylide": ("draw_molecule", {"smiles": "[CH2-][P+](c1ccccc1)(c1ccccc1)c1ccccc1", "label": "Wittig ylide"},
              {"tilt": 30, "lone_pairs": True}),
    "cubane": ("draw_named_molecule", {"name": "cubane"}, {"tilt": 30}),
    "ferrocene": ("draw_named_molecule", {"name": "ferrocene"}, {"tilt": 30}),
    "paclitaxel": ("draw_named_molecule", {"name": "paclitaxel"}, {"tilt": 25}),
}

HOST_PAGE = """<!doctype html><meta charset="utf-8"><body style="margin:0;background:#fff">
<iframe id="w" sandbox="allow-scripts" style="width:720px;height:900px;border:0"></iframe>
<script>
const frame = document.getElementById('w');
let toolCall = null;
window.startWidget = (html, args, result) => { toolCall = { args, result }; frame.srcdoc = html; };
window.addEventListener('message', async (event) => {
  if (event.source !== frame.contentWindow) return;
  const m = event.data;
  const reply = (body) => frame.contentWindow.postMessage({ jsonrpc: '2.0', id: m.id, ...body }, '*');
  if (m.method === 'ui/initialize') {
    reply({ result: { protocolVersion: m.params.protocolVersion, hostInfo: { name: 'demo-host', version: '1' },
                      hostCapabilities: {}, hostContext: { theme: 'light' } } });
  } else if (m.method === 'ui/notifications/initialized') {
    const send = (method, params) => frame.contentWindow.postMessage({ jsonrpc: '2.0', method, params }, '*');
    send('ui/notifications/tool-input', { arguments: toolCall.args });
    send('ui/notifications/tool-result', toolCall.result);
  } else if (m.method === 'tools/call') {
    try { reply({ result: JSON.parse(await window.pyCallTool(JSON.stringify(m.params))) }); }
    catch (e) { reply({ error: { code: -32000, message: String(e) } }); }
  }
});
</script>"""

TURN = """async ([ax, ay, d]) => {
  const v = window.moltalkViewer;
  v.rotate(ax, ay, d);
  await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
}"""

# Set the view for loop position t in [0, 1): "spin" turns once about an axis tilted toward the viewer; "rock"
# sways a planar molecule +-amplitude about the vertical axis with a slight nod, so it is never seen edge-on.
POSE = """async ([t, motion, tilt, amplitude]) => {
  const v = window.moltalkViewer;
  if (!v.base) v.base = v.R;
  const mul = (a, b) => a.map((row) => [0, 1, 2].map((j) => row[0] * b[0][j] + row[1] * b[1][j] + row[2] * b[2][j]));
  const rx = (a) => { const c = Math.cos(a), s = Math.sin(a); return [[1, 0, 0], [0, c, -s], [0, s, c]]; };
  const ry = (a) => { const c = Math.cos(a), s = Math.sin(a); return [[c, 0, s], [0, 1, 0], [-s, 0, c]]; };
  const turn = 2 * Math.PI * t;
  const R = motion === "spin" ? mul(rx(tilt), mul(ry(turn), rx(-tilt)))
                              : mul(ry(amplitude * Math.sin(turn)), rx(0.25 * amplitude * Math.sin(2 * turn)));
  v.R = mul(R, mul(rx(motion === "spin" ? tilt : 0), v.base));
  v.draw();
  await new Promise((r) => requestAnimationFrame(() => requestAnimationFrame(r)));
}"""


def dump(result):
    out = result.model_dump(by_alias=True, mode="json", exclude_none=True)
    return out


async def resolve_isomer(client, smiles, descriptor):
    result = await client.call_tool("enumerate_stereoisomers", {"smiles": smiles})
    for iso in result.structuredContent["isomers"]:
        if any(x.get("cip") == descriptor for x in iso.get("axial_stereo", [])):
            return iso["smiles"]
    raise SystemExit(f"no {descriptor} isomer for {smiles}")


async def render(client, page, html, name, tool, args, options, frames_dir):
    if "isomer" in options:
        args = {**args, "smiles": await resolve_isomer(client, *options["isomer"])}
    result = dump(await client.call_tool(tool, {**args, "numbering": "none"}))
    if result.get("isError"):
        raise SystemExit(f"{name}: {result['content'][0]['text']}")
    await page.set_content(HOST_PAGE)
    await page.evaluate("([h, a, r]) => startWidget(h, a, r)", [html, args, result])
    frame = page.frame_locator("#w")
    box = frame.locator(".liftable")
    try:
        await box.wait_for(timeout=20000)
    except Exception:
        await page.screenshot(path=str(frames_dir.parent / f"failed_{name}.png"))
        raise
    await page.wait_for_timeout(500)
    inner = page.frames[1]
    if options.get("stereo"):
        await frame.get_by_role("group", name="Stereo labels").get_by_role("button", name=options["stereo"]).click()
    if options.get("lone_pairs"):
        await frame.get_by_role("button", name="Show lone pairs").click()
    if not await inner.evaluate("() => window.moltalkViewer.ensureModel()"):
        raise SystemExit(f"{name}: no 3D model")
    # Lift fully into 3D (no turn), then pose each frame of one loop period.
    await inner.evaluate(TURN, [0, 0, 200])
    count = FPS * SECONDS
    motion = options.get("motion", "spin")
    tilt = options.get("tilt", 20) * 3.14159265 / 180
    amplitude = options.get("amplitude", 45) * 3.14159265 / 180
    paths = []
    for k in range(count):
        await inner.evaluate(POSE, [k / count, motion, tilt, amplitude])
        path = frames_dir / f"{name}_{k:03d}.png"
        await box.screenshot(path=str(path))
        paths.append(path)
    return paths


def crop_and_encode(name, paths):
    images = [Image.open(p).convert("RGB") for p in paths]
    white = Image.new("RGB", images[0].size, "white")
    box = None
    for im in images:
        b = ImageChops.difference(im, white).point(lambda v: 255 if v > 24 else 0).getbbox()
        if b:
            box = b if box is None else (min(box[0], b[0]), min(box[1], b[1]), max(box[2], b[2]), max(box[3], b[3]))
    w, h = images[0].size
    box = (max(0, box[0] - PAD), max(0, box[1] - PAD), min(w, box[2] + PAD), min(h, box[3] + PAD))
    for im, p in zip(images, paths):
        im.crop(box).save(p)
    width = min(WIDTH, box[2] - box[0])
    width -= width % 2
    pattern = str(paths[0].parent / f"{name}_%03d.png")
    OUT.mkdir(parents=True, exist_ok=True)
    scale = f"scale={width}:-2:flags=lanczos"
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pattern, "-vf",
                    f"{scale},split[a][b];[a]palettegen=max_colors=64:stats_mode=diff[p];"
                    f"[b][p]paletteuse=dither=bayer:bayer_scale=5:diff_mode=rectangle",
                    "-loop", "0", str(OUT / f"{name}.gif")], check=True)
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS), "-i", pattern, "-vf", scale,
                    "-c:v", "libwebp_anim", "-lossless", "0", "-q:v", "80", "-loop", "0",
                    str(OUT / f"{name}.webp")], check=True)
    for ext in ("gif", "webp"):
        print(f"  {name}.{ext}: {(OUT / f'{name}.{ext}').stat().st_size / 1024:.0f} KB")


async def main(names):
    env = {**os.environ, "MOLTALK_OFFLINE": "1"}
    params = StdioServerParameters(command=sys.executable, args=["-m", "moltalk.server"], env=env)
    async with stdio_client(params) as (read, write), ClientSession(read, write) as client:
        await client.initialize()
        widget = next(r for r in (await client.list_resources()).resources if str(r.uri).startswith("ui://widget/"))
        html = (await client.read_resource(widget.uri)).contents[0].text

        async def call_tool(raw):
            request = json.loads(raw)
            return json.dumps(dump(await client.call_tool(request["name"], request["arguments"])))

        async with async_playwright() as pw:
            browser = await pw.chromium.launch(channel=os.getenv("MOLTALK_UI_BROWSER_CHANNEL", "chrome"))
            for name in names:
                tool, args, options = DEMOS[name]
                print(name)
                page = await browser.new_page(viewport={"width": 760, "height": 940})  # a fresh page per demo
                await page.expose_function("pyCallTool", call_tool)
                with tempfile.TemporaryDirectory() as tmp:
                    paths = await render(client, page, html, name, tool, args, options, Path(tmp))
                    crop_and_encode(name, paths)
                await page.close()
            await browser.close()


if __name__ == "__main__":
    if not shutil.which("ffmpeg"):
        raise SystemExit("ffmpeg is required")
    chosen = sys.argv[1:] or list(DEMOS)
    unknown = [n for n in chosen if n not in DEMOS]
    if unknown:
        raise SystemExit(f"unknown demo(s) {unknown}; choose from {list(DEMOS)}")
    asyncio.run(main(chosen))
