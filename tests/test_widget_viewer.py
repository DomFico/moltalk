"""Viewer (fullscreen) layout and keeping the view across redraws, through the test_widget_ui host harness."""
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


def test_viewer_fits_without_scrolling_and_view_survives_toggles():
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
                await page.set_content(HOST_PAGE)
                frame = page.frame_locator("#w")
                args = {"smiles": "F[C@H]1C[C@@H](C)CCC1", "label": "cyclohexane test"}
                result = dump(await client.call_tool("draw_molecule", args))
                await page.evaluate("([h, a, r]) => startWidget(h, a, r, 'light')", [html, args, result])
                box = frame.locator(".liftable")
                await box.wait_for()
                await page.wait_for_timeout(600)
                assert await frame.get_by_role("button", name="Enter viewer").is_visible()

                # Rotate, then toggle hydrogens and atom numbers: the view must not reset.
                b = await box.bounding_box()
                await page.mouse.move(b["x"] + 200, b["y"] + 200)
                await page.mouse.down()
                await page.mouse.move(b["x"] + 330, b["y"] + 250, steps=10)
                await page.mouse.up()
                await page.wait_for_timeout(80)
                f_label = "svg => { const t = [...svg.querySelectorAll('text')].find(t => t.textContent === 'F'); return [+t.getAttribute('x'), +t.getAttribute('y')]; }"
                before = await frame.locator("svg.lift").evaluate(f_label)
                numbers = frame.get_by_role("group", name="Atom numbers")
                for name, scope in (("Show hydrogens", frame), ("Off", numbers)):
                    await scope.get_by_role("button", name=name, exact=True).click()
                    await page.wait_for_timeout(1200)
                    assert await frame.locator("svg.lift").is_visible(), name
                    after = await frame.locator("svg.lift").evaluate(f_label)
                    assert abs(after[0] - before[0]) < 3 and abs(after[1] - before[1]) < 3, (name, before, after)

                # Inline with a host height limit: everything fits, the page never scrolls, only the details box does.
                send = "p => document.getElementById('w').contentWindow.postMessage({jsonrpc: '2.0', method: 'ui/notifications/host-context-changed', params: p}, '*')"
                await page.evaluate(send, {"displayMode": "inline", "containerDimensions": {"maxHeight": 560}})
                await page.wait_for_timeout(300)
                inline = await page.frames[1].evaluate("""() => { const d = document.querySelector('.details');
                    return {main: document.querySelector('main').getBoundingClientRect().height, pageScroll: getComputedStyle(document.documentElement).overflow,
                            detailsScrolls: d.scrollHeight > d.clientHeight, detailsOverflow: getComputedStyle(d).overflowY}; }""")
                assert inline["main"] <= 560 and inline["pageScroll"] == "hidden", inline
                assert inline["detailsScrolls"] and inline["detailsOverflow"] == "auto", inline
                assert await page.evaluate("window.lastHeight") <= 560

                # Host switches to the viewer on a phone-sized window with a 60 px title bar.
                await page.evaluate("document.getElementById('w').style.cssText = 'width:390px;height:700px;border:0'")
                await page.evaluate(send, {"displayMode": "fullscreen", "safeAreaInsets": {"top": 60, "right": 0, "bottom": 20, "left": 0}})
                await page.wait_for_timeout(300)
                inner = page.frames[1]
                fit = await inner.evaluate("""() => { const m = getComputedStyle(document.querySelector('main'));
                    return {scroll: document.documentElement.scrollHeight, view: innerHeight, top: parseFloat(m.paddingTop),
                            bottom: parseFloat(m.paddingBottom), tools: document.querySelector('.lifthint').getBoundingClientRect().top,
                            drawing: document.querySelector('.liftable').getBoundingClientRect().top}; }""")
                assert fit["scroll"] <= fit["view"], fit  # nothing to scroll
                assert fit["top"] >= 100 and fit["bottom"] >= 120, fit  # clear of ChatGPT's title bar and chat box
                assert fit["tools"] < fit["drawing"], fit  # our controls sit above the drawing, never under the chat box
                assert not await frame.get_by_role("button", name="Enter viewer").is_visible()
                assert await frame.get_by_role("button", name="Exit viewer").is_visible()
                if SCREENSHOT_DIR:
                    await page.screenshot(path=os.path.join(SCREENSHOT_DIR, "viewer_phone.png"))
                assert (await box.bounding_box())["height"] > 280  # the drawing takes the free space (the Stereo selector costs a toolbar row on phones)

                # Back inline: the widget re-reports its height (otherwise the host may collapse it).
                await page.evaluate("window.lastHeight = 0; document.getElementById('w').style.cssText = 'width:720px;height:900px;border:0'")
                await page.evaluate(send, {"displayMode": "inline"})
                await page.wait_for_timeout(500)
                assert await page.evaluate("window.lastHeight") > 300
                assert not await frame.get_by_role("button", name="Exit viewer").is_visible()
                await browser.close()

    asyncio.run(run())


def test_inline_fits_a_host_that_caps_height_silently():
    # Codex desktop: a wide window and a capped inline frame, with no maxHeight reported. The drawing is capped at
    # 420 px inline, and if the frame stays shorter than requested the widget shrinks to fit it.
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
                page = await browser.new_page(viewport={"width": 1300, "height": 1000})
                await page.expose_function("pyCallTool", py_call_tool)
                await page.set_content(HOST_PAGE)
                await page.evaluate("document.getElementById('w').style.cssText = 'width:1200px;height:900px;border:0'")
                args = {"smiles": "CC(=O)Oc1ccccc1C(=O)O", "label": "aspirin"}
                result = dump(await client.call_tool("draw_molecule", args))
                await page.evaluate("([h, a, r]) => startWidget(h, a, r, 'light')", [html, args, result])
                frame = page.frame_locator("#w")
                box = frame.locator(".liftable")
                await box.wait_for()
                await page.wait_for_timeout(300)
                assert (await box.bounding_box())["height"] <= 421  # wide window: capped, not 1200 * 420 / 640
                # The host now holds the frame at 520 px without saying so.
                await page.evaluate("document.getElementById('w').style.height = '520px'")
                await page.wait_for_timeout(2000)
                inner = page.frames[1]
                fits = await inner.evaluate("() => document.getElementById('root').getBoundingClientRect().height")
                assert fits <= 520 and await page.evaluate("window.lastHeight") <= 520
                await browser.close()

    asyncio.run(run())
