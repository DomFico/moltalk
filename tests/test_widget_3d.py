"""Lifting the flat drawing into 3D, in a real browser through the test_widget_ui host harness."""
import asyncio
import json
import os
import sys
import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from moltalk.server import WIDGET_URI
from test_widget_ui import HOST_PAGE, SCREENSHOT_DIR, dump
from test_depiction import C60

playwright_api = pytest.importorskip("playwright.async_api")


def test_lift_flat_drawing_into_3d():
    async def run():
        params = StdioServerParameters(command=sys.executable, args=["-m", "moltalk.server"])
        async with stdio_client(params) as (read, write), ClientSession(read, write) as client:
            await client.initialize()
            html = (await client.read_resource(WIDGET_URI)).contents[0].text
            calls = []

            async def py_call_tool(raw):
                request = json.loads(raw)
                calls.append(request["name"])
                return json.dumps(dump(await client.call_tool(request["name"], request["arguments"])))

            async with playwright_api.async_playwright() as pw:
                channel = os.getenv("MOLTALK_UI_BROWSER_CHANNEL", "chrome")
                browser = await pw.chromium.launch(channel=None if channel == "chromium" else channel)
                page = await browser.new_page(viewport={"width": 760, "height": 1000})
                await page.expose_function("pyCallTool", py_call_tool)
                await page.set_content(HOST_PAGE)
                frame = page.frame_locator("#w")
                box = frame.locator(".liftable")
                lift = frame.locator("svg.lift")

                async def show(args):
                    result = dump(await client.call_tool("draw_molecule", args))
                    await page.evaluate("([h, a, r]) => startWidget(h, a, r, 'light')", [html, args, result])
                    await box.wait_for()
                    return result

                async def drag(dx, dy, steps=10, start=(200, 200)):
                    b = await box.bounding_box()
                    await page.mouse.move(b["x"] + start[0], b["y"] + start[1])
                    await page.mouse.down()
                    await page.mouse.move(b["x"] + start[0] + dx, b["y"] + start[1] + dy, steps=steps)
                    await page.mouse.up()
                    await page.wait_for_timeout(60)

                async def shot(name):
                    if SCREENSHOT_DIR:
                        await page.screenshot(path=os.path.join(SCREENSHOT_DIR, name), full_page=True)

                async def lines():
                    return await lift.evaluate("svg => [...svg.querySelectorAll('line')].map(l => [l.getAttribute('x1'), l.getAttribute('y1')].join()).join(' ')")

                result = await show({"smiles": "F[C@H]1C[C@@H](C)CCC1", "label": "cyclohexane test"})
                atom_px = result["_meta"]["atom_px"]
                assert await lift.is_hidden()
                await shot("lift_0_flat.png")

                # A tiny drag starts exactly on the flat drawing: the F label sits on RDKit's F position.
                await box.hover()
                await page.wait_for_timeout(300)  # prefetch on hover
                await drag(4, 0, steps=2)
                assert await lift.is_visible() and calls[-1] == "conformer_3d"
                f = await lift.evaluate("svg => { const t = [...svg.querySelectorAll('text')].find(t => t.textContent === 'F'); return [+t.getAttribute('x'), +t.getAttribute('y')]; }")
                vb = await frame.locator(".liftable svg.lift").bounding_box()
                assert abs(f[0] - atom_px[0][0]) < 6 and abs(f[1] - atom_px[0][1]) < 6, (f, atom_px[0])

                # A long drag fully lifts and rotates; halos (bond gaps) are drawn.
                before = await lines()
                await drag(160, 50)
                assert await lines() != before
                assert await lift.evaluate("svg => [...svg.querySelectorAll('line')].some(l => l.getAttribute('stroke') === '#ffffff')")
                await shot("lift_1_cyclohexane.png")
                assert await frame.get_by_text("Rotated view.").is_visible()

                # Reset animates back and restores the exact RDKit drawing.
                await frame.get_by_role("button", name="2D", exact=True).click()
                await lift.wait_for(state="hidden")
                assert await frame.locator(".liftable > svg:not(.lift)").evaluate("s => s.style.visibility") == "visible"

                # Keyboard: arrows lift and rotate; Escape flattens.
                await box.focus()
                for _ in range(4):
                    await page.keyboard.press("ArrowRight")
                await page.wait_for_timeout(60)
                assert await lift.is_visible()
                await page.keyboard.press("Escape")
                await lift.wait_for(state="hidden")

                # Aspirin: hydrogens in labels and ring double bonds drawn on the inner side.
                await show({"smiles": "CC(=O)Oc1ccccc1C(=O)O", "label": "aspirin"})
                await drag(30, 10, steps=4)
                await lift.locator("text").first.wait_for()  # a drag before the 3D data arrives is applied when it does
                texts = await lift.evaluate("svg => [...svg.querySelectorAll('text')].map(t => t.textContent)")
                assert "H" in texts and texts.count("O") == 4
                await shot("lift_2_aspirin_partial.png")
                await drag(150, 70)
                await shot("lift_3_aspirin.png")

                # Explicit-H drawings are laid out in their own orientation (LSD's is turned 180°): the lift must still
                # start exactly from the drawing on screen, with no flip.
                lsd = await show({"smiles": "CCN(CC)C(=O)[C@@H]1C=C2c3cccc4[nH]cc(c34)C[C@H]2N(C)C1", "label": "LSD",
                                  "hydrogens": True, "atom_indices": False})
                await page.wait_for_timeout(1500)
                await drag(4, 0, steps=2)
                await lift.wait_for()
                labels = await lift.evaluate("svg => [...svg.querySelectorAll('text')].filter(t => ['N', 'O'].includes(t.textContent)).map(t => [+t.getAttribute('x'), +t.getAttribute('y')])")
                px = lsd["_meta"]["atom_px"]
                hetero = [px[i] for i, ch in enumerate("CCNCCCOCCCcccccNcccCCNCC") if ch in "NO"]
                for x, y in hetero:
                    assert min(abs(x - lx) + abs(y - ly) for lx, ly in labels) < 8, (x, y, labels)

                # Lone pairs: dots on the flat drawing and the rotated one, and the setting survives a redraw.
                await show({"smiles": "CC(=O)[O-]", "label": "acetate", "numbering": "none"})
                await frame.get_by_role("button", name="Show lone pairs").click()
                assert await frame.locator(".liftable svg:not(.lift) circle.lp-dot").count() == 10  # O 2 pairs, O- 3 pairs
                await drag(120, 30)
                await lift.wait_for()
                assert await lift.locator("circle.lp-dot").count() == 10
                await frame.get_by_role("button", name="Show hydrogens").click()
                await frame.get_by_role("button", name="Hide lone pairs").wait_for()
                await page.wait_for_timeout(800)
                assert await lift.locator("circle.lp-dot").count() == 10
                await frame.get_by_role("button", name="Hide lone pairs").click()
                assert await lift.locator("circle.lp-dot").count() == 0
                await show({"smiles": "CCCC", "label": "butane"})
                assert await frame.get_by_role("button", name="Show lone pairs").is_disabled()

                # Unspecified stereo: not labelled by default; under Stereo "All" it is marked as arbitrary.
                await show({"smiles": "CC(O)C(=O)O", "label": "lactic acid"})
                await drag(120, 0)
                texts = "svg => [...svg.querySelectorAll('text')].map(t => t.textContent)"
                assert not any("arb." in t for t in await lift.evaluate(texts))
                await frame.get_by_role("group", name="Stereo labels").get_by_role("button", name="All").click()
                assert any("(arb. " in t for t in await lift.evaluate(texts))
                await frame.get_by_role("group", name="Stereo labels").get_by_role("button", name="Specified").click()
                await frame.get_by_text("Model caveat.").wait_for()

                # Zoom: the + key acts on the flat drawing too; "2D" restores the view.
                flat_svg = frame.locator(".liftable > svg:not(.lift)")
                await show({"smiles": "CC(=O)Oc1ccccc1C(=O)O", "label": "aspirin"})
                await box.focus()
                await box.press("+")
                await page.wait_for_timeout(60)
                mid = float((await flat_svg.get_attribute("viewBox")).split()[2])
                await page.wait_for_timeout(500)
                end = float((await flat_svg.get_attribute("viewBox")).split()[2])
                assert end < mid <= 640 and abs(end - 640 / 1.2) < 2  # glides smoothly to a gentle 1.2x
                await frame.get_by_role("button", name="2D", exact=True).click()
                await page.wait_for_timeout(350)
                assert await flat_svg.get_attribute("viewBox") == "0 0 640 420"

                # Hydrogens: the toggle redraws with explicit H, and the rotated view carries them too.
                await frame.get_by_role("button", name="Show hydrogens").click()
                await frame.get_by_role("button", name="Hide hydrogens").wait_for()
                await page.wait_for_timeout(400)
                await drag(120, 40)
                await lift.wait_for()
                texts = await lift.evaluate("svg => [...svg.querySelectorAll('text')].map(t => t.textContent)")
                assert texts.count("H") == 8  # aspirin's eight hydrogens as atoms
                await shot("lift_5_aspirin_h.png")

                # Heme: Fe bonded to the ring nitrogens, two of them drawn as coordinate arrows.
                from test_coordination import HEME
                await show({"smiles": HEME, "label": "heme b", "atom_indices": False})
                await frame.get_by_text("Metal bonds.").wait_for()
                await shot("lift_6_heme_flat.png")
                await drag(90, 30)
                await lift.wait_for()
                assert await lift.locator("polygon").count() == 2
                await shot("lift_6_heme.png")

                for name, smiles in [("c60", C60), ("adamantane", "C1C2CC3CC1CC(C2)C3"),
                                     ("cis_decalin", "C1CC[C@H]2CCCC[C@H]2C1")]:
                    await show({"smiles": smiles, "label": name, "atom_indices": False})
                    await drag(140, 60)  # may start before the 3D data arrives; the movement is kept
                    await lift.wait_for()
                    await shot(f"lift_4_{name}.png")
                await browser.close()

    asyncio.run(run())
