"""Phone-sized touch screen: a one-finger drag must rotate the molecule, not scroll the page; two fingers zoom."""
import asyncio
import json
import os
import sys
import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client
from moltalk.server import WIDGET_URI
from test_widget_ui import HOST_PAGE, SCREENSHOT_DIR, dump

playwright_api = pytest.importorskip("playwright.async_api")
# A host page taller than the phone screen, so a stray swipe would scroll it.
TALL_HOST = HOST_PAGE.replace('style="width:720px;height:900px;border:0"', 'style="width:100%;height:700px;border:0"') \
                     .replace("</script>", "</script><div style='height:2000px'>more chat below</div>")


def test_touch_rotates_and_pinch_zooms():
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
                context = await browser.new_context(viewport={"width": 390, "height": 844}, has_touch=True,
                                                    is_mobile=True, device_scale_factor=3)
                page = await context.new_page()
                await page.expose_function("pyCallTool", py_call_tool)
                await page.set_content(TALL_HOST)
                args = {"smiles": "F[C@H]1C[C@@H](C)CCC1", "label": "cyclohexane test"}
                result = dump(await client.call_tool("draw_molecule", args))
                await page.evaluate("([h, a, r]) => startWidget(h, a, r, 'light')", [html, args, result])
                frame = page.frame_locator("#w")
                box = frame.locator(".liftable")
                await box.wait_for()
                await page.wait_for_timeout(800)  # 3D data prefetch
                cdp = await context.new_cdp_session(page)
                b = await box.bounding_box()
                cx, cy = b["x"] + b["width"] / 2, b["y"] + b["height"] / 2

                async def touch(kind, points):
                    await cdp.send("Input.dispatchTouchEvent", {"type": kind, "touchPoints": [{"x": x, "y": y, "id": i} for i, (x, y) in enumerate(points)]})

                # One finger, mostly vertical swipe across the drawing: rotates; the page does not scroll.
                await touch("touchStart", [(cx, cy)])
                for k in range(1, 12):
                    await touch("touchMove", [(cx + 4 * k, cy + 10 * k)])
                await touch("touchEnd", [])
                await page.wait_for_timeout(100)
                assert await page.evaluate("window.scrollY") == 0
                assert await frame.locator("svg.lift").is_visible()
                assert "pinch to zoom" in await frame.locator(".lifthint").inner_text()
                if SCREENSHOT_DIR:
                    await page.screenshot(path=os.path.join(SCREENSHOT_DIR, "touch_rotated.png"))

                # Two fingers moving apart: zoom in.
                before = await frame.locator("svg.lift").get_attribute("viewBox")
                await touch("touchStart", [(cx - 20, cy), (cx + 20, cy)])
                for k in range(1, 8):
                    await touch("touchMove", [(cx - 20 - 8 * k, cy), (cx + 20 + 8 * k, cy)])
                await touch("touchEnd", [])
                await page.wait_for_timeout(100)
                after = await frame.locator("svg.lift").get_attribute("viewBox")
                assert float(after.split()[2]) < float(before.split()[2])  # narrower view = zoomed in
                assert await page.evaluate("window.scrollY") == 0
                await browser.close()

    asyncio.run(run())
