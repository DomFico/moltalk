"""Render the app icon PNGs from assets/moltalk-logo.svg (headless Chrome via Playwright, transparent corners).

    .venv/bin/python scripts/render_icons.py

Writes the 512 px icon to assets/, plugin/assets/ (uploaded with the plugin package) and moltalk/static/, and the
128 px icon that the server advertises in serverInfo.icons.
"""
import asyncio
import io
from pathlib import Path

from PIL import Image
from playwright.async_api import async_playwright

ROOT = Path(__file__).resolve().parent.parent
SVG = ROOT / "assets" / "moltalk-logo.svg"
LARGE = [ROOT / "assets" / "moltalk-logo.png", ROOT / "plugin" / "assets" / "moltalk-logo.png",
         ROOT / "moltalk" / "static" / "icon-512.png"]
SMALL = ROOT / "moltalk" / "static" / "icon-128.png"


async def render(size: int) -> Image.Image:
    async with async_playwright() as pw:
        browser = await pw.chromium.launch(channel="chrome")
        page = await browser.new_page(viewport={"width": size, "height": size})
        svg = SVG.read_text().replace("<svg ", f'<svg style="display:block;width:{size}px;height:{size}px" ', 1)
        await page.set_content(f"<html><body style='margin:0;background:transparent'>{svg}</body></html>")
        png = await page.screenshot(omit_background=True)
        await browser.close()
    return Image.open(io.BytesIO(png)).convert("RGBA")


async def main():
    large = await render(512)
    for path in LARGE:
        large.save(path, optimize=True)
    # 128 px: rendered large and downsampled, which keeps the 9-unit strokes smooth.
    (await render(1024)).resize((128, 128), Image.LANCZOS).save(SMALL, optimize=True)
    print("wrote", ", ".join(str(p.relative_to(ROOT)) for p in [*LARGE, SMALL]))


if __name__ == "__main__":
    asyncio.run(main())
