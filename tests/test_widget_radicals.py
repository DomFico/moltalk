"""A radical's unpaired electron is drawn once in the flat drawing: RDKit's own dot when lone pairs are hidden, the
viewer's electron overlay (and not RDKit's dot as well) when they are shown."""
import asyncio
import json
import os
import sys

import pytest
from mcp import ClientSession, StdioServerParameters
from mcp.client.stdio import stdio_client

from moltalk.server import WIDGET_URI
from test_widget_ui import HOST_PAGE, dump

playwright_api = pytest.importorskip("playwright.async_api")

VISIBLE_DOTS = """() => {
  const svg = document.querySelector('svg:not(.lift)');
  const rdkit = [...svg.querySelectorAll('[data-radical-glyph]')].filter(n => n.style.display !== 'none').length;
  const overlay = [...svg.querySelectorAll('g.lone-pairs circle')].length;
  return {rdkit, overlay};
}"""


def test_methyl_radical_dot_is_drawn_once():
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
                page = await browser.new_page(viewport={"width": 760, "height": 900})
                await page.expose_function("pyCallTool", py_call_tool)
                await page.set_content(HOST_PAGE)
                frame = page.frame_locator("#w")
                args = {"smiles": "[CH3]"}
                result = dump(await client.call_tool("draw_molecule", args))
                await page.evaluate("([h, a, r]) => startWidget(h, a, r, 'light')", [html, args, result])
                await frame.locator(".liftable").wait_for()
                await page.wait_for_timeout(400)
                inner = page.frames[1]
                await frame.get_by_role("button", name="Show lone pairs").click()
                await page.wait_for_timeout(300)
                shown = await inner.evaluate(VISIBLE_DOTS)
                assert shown["rdkit"] == 0 and shown["overlay"] == 1, shown  # one dot, the overlay's
                await frame.get_by_role("button", name="Hide lone pairs").click()
                await page.wait_for_timeout(300)
                hidden = await inner.evaluate(VISIBLE_DOTS)
                assert hidden["rdkit"] == 1 and hidden["overlay"] == 0, hidden  # RDKit's dot alone
                await browser.close()

    asyncio.run(run())
