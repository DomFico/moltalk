"""Stereo labels, image export and structure-file download in a real browser, through the test_widget_ui host
harness extended with MCP Apps' ui/download-file (and, separately, a host without it)."""
import asyncio
import base64
import io
import json
import os
import sys
import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from moltalk.server import WIDGET_URI
from test_widget_ui import HOST_PAGE, SCREENSHOT_DIR, dump

playwright_api = pytest.importorskip("playwright.async_api")

# The same host, but it advertises downloadFile (when window.canDownload) and records each download request.
EXPORT_HOST = HOST_PAGE.replace("hostCapabilities: {}", "hostCapabilities: window.canDownload ? { downloadFile: {} } : {}").replace(
    "} else if (m.method === 'ui/notifications/size-changed') {",
    "} else if (m.method === 'ui/download-file') {\n    window.downloads.push(m.params); reply({ result: {} });\n"
    "  } else if (m.method === 'ui/notifications/size-changed') {").replace(
    "window.hostLog = [];", "window.hostLog = []; window.downloads = []; window.canDownload = true;")


def _png_alpha(b64):
    from PIL import Image
    image = Image.open(io.BytesIO(base64.b64decode(b64))).convert("RGBA")
    alphas = [a for *_, a in image.getdata()]
    return image.size, min(alphas), max(alphas)


def test_stereo_labels_and_exports():
    async def run():
        params = StdioServerParameters(command=sys.executable, args=["-m", "moltalk.server"])
        async with stdio_client(params) as (read, write), ClientSession(read, write) as client:
            await client.initialize()
            html = (await client.read_resource(WIDGET_URI)).contents[0].text

            async def py_call_tool(raw):
                request = json.loads(raw)
                return json.dumps(dump(await client.call_tool(request["name"], request["arguments"])))

            async with playwright_api.async_playwright() as pw:
                channel = os.getenv("MOLTALK_UI_BROWSER_CHANNEL", "chrome")
                browser = await pw.chromium.launch(channel=None if channel == "chromium" else channel)
                page = await browser.new_page(viewport={"width": 760, "height": 1000})
                await page.expose_function("pyCallTool", py_call_tool)
                await page.set_content(EXPORT_HOST)
                frame = page.frame_locator("#w")

                async def show(name, args):
                    result = dump(await client.call_tool(name, args))
                    await page.evaluate("([h, a, r]) => startWidget(h, a, r, 'light')", [html, args, result])
                    return result

                async def downloads():
                    return await page.evaluate("window.downloads")

                # 1. Stereo labels. Atom 1 is specified (S); atom 4 is not.
                await show("draw_molecule", {"smiles": "C[C@H](O)CC(C)Cl", "label": "test"})
                box = frame.locator(".liftable")
                await box.wait_for()
                stereo = frame.get_by_role("group", name="Stereo labels")
                visible = """() => [...document.querySelectorAll('.liftable > svg:not(.lift) .CIP_Code')]
                                 .filter(n => n.style.display !== 'none').length"""
                inner = page.frames[1]
                assert await inner.evaluate(visible) == 3  # "(S)" only; the open centre's "(?)" is hidden
                await stereo.get_by_role("button", name="All").click()
                assert await inner.evaluate(visible) == 6  # plus "(?)", placed by RDKit like any stereo label
                await stereo.get_by_role("button", name="Off").click()
                assert await inner.evaluate(visible) == 0

                # Rotated view: "(S)" is shown for Specified; the open centre only under All, as "arb.".
                await stereo.get_by_role("button", name="Specified").click()
                await box.focus()
                for _ in range(4):
                    await box.press("ArrowRight")
                lift = frame.locator("svg.lift")
                await lift.wait_for()
                await page.wait_for_timeout(200)
                tags = await lift.evaluate("s => [...s.querySelectorAll('text')].map(t => t.textContent).join(' ')")
                assert "(S)" in tags and "arb." not in tags
                await stereo.get_by_role("button", name="All").click()
                tags = await lift.evaluate("s => [...s.querySelectorAll('text')].map(t => t.textContent).join(' ')")
                assert "(arb. " in tags
                await stereo.get_by_role("button", name="Off").click()
                tags = await lift.evaluate("s => [...s.querySelectorAll('text')].map(t => t.textContent).join(' ')")
                assert "(S)" not in tags and "arb." not in tags

                # 2. Image export of the rotated view through ui/download-file: transparent, no white marks.
                await frame.get_by_role("button", name="SVG").click()
                await page.wait_for_function("window.downloads.length === 1")
                svg = (await downloads())[0]["contents"][0]["resource"]
                assert svg["uri"] == "file:///test.svg" and svg["mimeType"] == "image/svg+xml"
                assert 'fill="#ffffff"' not in svg["text"].split("<mask")[0]  # no white background or label backings
                assert 'stroke="#ffffff"' not in svg["text"].replace('fill="#ffffff"', "")  # halos became masks
                await frame.get_by_role("button", name="PNG").click()
                await page.wait_for_function("window.downloads.length === 2")
                png = (await downloads())[1]["contents"][0]["resource"]
                assert png["uri"] == "file:///test.png" and png["mimeType"] == "image/png"
                (width, height), lowest, highest = _png_alpha(png["blob"])
                assert width >= 300 and lowest == 0 and highest == 255  # cropped, 3x; transparent background, opaque ink

                # 3. A crowded rotated view (crossing bonds) still exports: masks only where something is cut.
                await show("draw_molecule", {"smiles": "C1CC2CCC1CC2", "label": "bicyclo"})
                await frame.locator(".liftable").focus()
                for _ in range(6):
                    await frame.locator(".liftable").press("ArrowDown")
                await page.wait_for_timeout(200)
                await frame.get_by_role("button", name="SVG").click()
                await page.wait_for_function("window.downloads.length === 3")

                # 4. Structure file from export_structure: card with a Download button.
                await show("export_structure", {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "name": "aspirin"})
                await frame.get_by_role("button", name="Download aspirin.cdxml").click()
                await page.wait_for_function("window.downloads.length === 4")
                cdxml = (await downloads())[3]["contents"][0]["resource"]
                assert cdxml["uri"] == "file:///aspirin.cdxml" and cdxml["text"].lstrip().startswith("<?xml")
                await show("export_structure", {"smiles": "CC(O)CC(C)Cl", "format": "sdf", "coordinates": "3d"})
                assert "arbitrary" in (await frame.locator(".filecard").inner_text())
                if SCREENSHOT_DIR:
                    await page.screenshot(path=os.path.join(SCREENSHOT_DIR, "export-card.png"), full_page=True)

                # 5. A host without ui/download-file: the widget offers the file in a save/copy panel instead.
                await page.evaluate("window.canDownload = false")
                await show("export_structure", {"smiles": "CCO", "format": "mol", "name": "ethanol"})
                await frame.get_by_role("button", name="Download ethanol.mol").click()
                panel = frame.locator(".saveas")
                await panel.wait_for()
                assert "ethanol.mol" in await panel.inner_text()
                assert "M  END" in await panel.locator("textarea").input_value()
                assert len(await downloads()) == 4
                await browser.close()

    asyncio.run(run())
