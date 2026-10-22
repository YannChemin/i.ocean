"""Example scene of i.ocean: the README image, and a test of the full effect stack.

A synthetic temperate coast (fractal terrain, an offshore island, a shelf
falling to 3 km) is rendered with waves, Phillips-spectrum foam, the
latitude gradient and the maritime decorations. The decoration vector map
is drawn over it at a finer resolution than the raster, and three insets
zoom on vessels: their silhouettes stay sharp while the ocean raster shows
its 100 m cells.

Test:           PYTHONPATH=$(grass --config python_path) pytest tests/
Make the image: PYTHONPATH=$(grass --config python_path) python3 tests/example_test.py [example.png]

SPDX-License-Identifier: Unlicense
"""

import io
import os
import sys
import tempfile
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont

import grass.script as gs
from grass.tools import Tools

# Scene: UTM 30N, about 47 deg N (temperate palette), 100 km x 60 km at 100 m.
REGION = dict(w=380000, e=480000, s=5180000, n=5240000, res=100)
DECO_SCALE = 6
WIDTH, HEIGHT = 1600, 960
INSET = 520
# Hull fill per decoration category (cat = type index + 1), muted sea-chart tones.
FILL = {1: "70:100:160", 2: "140:170:190", 3: "190:210:220", 4: "90:50:30", 5: "180:145:15",
        6: "110:75:35", 7: "85:75:60", 8: "45:50:55", 9: "35:40:50", 10: "120:128:135",
        11: "45:50:58", 12: "40:45:42", 13: "80:70:55", 14: "45:52:48", 15: "45:55:70",
        16: "60:70:85", 17: "30:30:36", 18: "20:20:20", 19: "110:95:80", 20: "100:115:125",
        21: "150:80:40", 22: "140:60:40", 23: "205:210:215", 24: "70:75:70",
        25: "120:45:45", 26: "58:64:74", 27: "85:100:75"}
# Vessels zoomed in the insets, by preference: container ship, crude tanker,
# LNG carrier, attack submarine, bulk carrier, cargo ship, strike boat.
INSET_CATS = [26, 22, 23, 8, 24, 27, 10]
VESSEL_CATS = set(range(7, 18)) | set(range(22, 28))
LABEL = {"container_ship": "Container ship", "crude_tanker": "Crude oil tanker", "lng_carrier": "LNG carrier",
         "hunter_sub": "Attack submarine", "bulk_carrier": "Bulk carrier", "cargo_ship": "General cargo ship",
         "barracuda_boat": "Fast strike boat"}


def font(size):
    for f in ("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",
              "/usr/share/fonts/TTF/DejaVuSans.ttf"):
        if os.path.exists(f):
            return ImageFont.truetype(f, size)
    return ImageFont.load_default()


def scene(tools):
    """Terrain, sea mask, bathymetry and the i.ocean render (maps in the session)."""
    tools.g_region(**REGION)
    tools.r_surf_fractal(output="fractal", dimension=2.12, seed=11)
    st = tools.r_univar(map="fractal", format="json")
    lo, hi = st["min"], st["max"]
    w, e, s, n = REGION["w"], REGION["e"], REGION["s"], REGION["n"]
    # A coast running SW-NE, land to the south-east, an island offshore.
    tools.r_mapcalc(expression=(
        f"elev = 2600.0 * ((x() - {w}) / {e - w} - 0.5) - 1700.0 * ((y() - {s}) / {n - s} - 0.5)"
        f" + 900.0 * ((fractal - {lo}) / {hi - lo} - 0.5)"
        f" + 1300.0 * exp(-((x() - {w + 0.33 * (e - w)})^2 + (y() - {s + 0.62 * (n - s)})^2)"
        f" / (2.0 * 5500.0^2)) - 150.0"))
    tools.r_mapcalc(expression="sea = if(elev < 0, 1, null())")
    tools.r_mapcalc(expression="depth = if(elev < 0, -elev, null())")
    tools.r_mapcalc(expression="land = if(elev >= 0, elev, null())")
    tools.r_colors(map="land", rules=io.StringIO("0 222:214:178\n40 168:178:120\n300 118:142:92\n"
                                                 "800 146:128:98\n1500 188:180:168\n2600 245:245:245\n"))
    tools.r_relief(input="elev", output="relief_all", zscale=2)
    tools.r_mapcalc(expression="relief = if(elev >= 0, relief_all, null())")
    tools.i_ocean(input="sea", depth="depth", output="ocean", flags="wfld",
                  deco_scale=DECO_SCALE, style="auto")
    tools.r_contour(input="elev", output="coast", levels=0)
    # Area fill colours of the decorations, from the category; the deck
    # details are drawn as white boundaries and lines over the fill, so the
    # per-category colour table written by i.ocean is not used here.
    tools.v_colors(map="ocean", flags="r")
    tools.v_db_addcolumn(map="ocean", columns="rgb VARCHAR(11)")
    for cat, rgb in FILL.items():
        tools.v_db_update(map="ocean", column="rgb", value=rgb, where=f"cat = {cat}")


def draw(tools, env, png, width, height):
    """Render the current region to `png` (each d.* adds a layer to the file)."""
    Path(png).unlink(missing_ok=True)
    env.update({"GRASS_RENDER_IMMEDIATE": "cairo", "GRASS_RENDER_FILE": str(png),
                "GRASS_RENDER_WIDTH": str(width), "GRASS_RENDER_HEIGHT": str(height),
                "GRASS_RENDER_FILE_READ": "TRUE", "GRASS_RENDER_TRUECOLOR": "TRUE",
                "GRASS_RENDER_BACKGROUNDCOLOR": "000000"})
    t = Tools(env=env)
    t.d_rast(map="ocean")
    t.d_shade(shade="relief", color="land", brighten=45, flags="n")
    t.d_vect(map="coast", color="235:235:225", width=1)
    t.d_vect(map="ocean", type="area", rgb_column="rgb", color="none")
    t.d_vect(map="ocean", type="boundary", color="240:240:235", width=1)
    t.d_vect(map="ocean", type="line", color="250:250:245", width=1)


def vessel_centres(tools):
    """Centre (x, y), category and scaled length of the placed vessels."""
    out = tools.v_out_ascii(input="ocean", type="centroid", format="point", separator="comma",
                            layer=1).text
    pts = np.array([[float(v) for v in l.split(",")[:3]] for l in out.splitlines() if l.strip()])
    rows = tools.v_db_select(map="ocean", format="json")["records"]
    length = {r["cat"]: r["len_m"] for r in rows}
    seen, centres = set(), []
    for x, y, c in pts:
        c = int(c)
        if c not in VESSEL_CATS:
            continue
        near = pts[(pts[:, 2] == c) & (np.hypot(pts[:, 0] - x, pts[:, 1] - y) < length[c])]
        key = (c, round(near[:, 0].mean(), -2), round(near[:, 1].mean(), -2))
        if key[:1] + key[1:] in seen:
            continue
        seen.add(key)
        centres.append((near[:, 0].mean(), near[:, 1].mean(), c, length[c]))
    return centres


def render_example(tools, env, png):
    """Main view plus three vessel insets, composed into `png`."""
    tmp = Path(tempfile.mkdtemp())
    tools.g_region(**REGION)
    draw(tools, dict(env), tmp / "main.png", WIDTH, HEIGHT)
    main = Image.open(tmp / "main.png").convert("RGB")
    rows = tools.v_db_select(map="ocean", format="json")["records"]
    name = {r["cat"]: r["type"] for r in rows}
    centres = vessel_centres(tools)
    picks = []
    for c in INSET_CATS:
        hit = [v for v in centres if v[2] == c]
        if hit and len(picks) < 3:
            picks.append(hit[0])
    canvas = Image.new("RGB", (WIDTH, HEIGHT + INSET + 40), (12, 18, 28))
    canvas.paste(main, (0, 0))
    d = ImageDraw.Draw(main)
    dc = ImageDraw.Draw(canvas)
    f_title, f_lab = font(22), font(19)
    gap = (WIDTH - 3 * INSET) // 4
    for k, (x, y, c, L) in enumerate(picks):
        h = 0.65 * L
        tools.g_region(n=y + h, s=y - h, e=x + h, w=x - h, res=REGION["res"])
        draw(tools, dict(env), tmp / f"inset{k}.png", INSET, INSET)
        ox = gap + k * (INSET + gap)
        canvas.paste(Image.open(tmp / f"inset{k}.png").convert("RGB"), (ox, HEIGHT + 30))
        dc.rectangle([ox - 1, HEIGHT + 29, ox + INSET, HEIGHT + 30 + INSET], outline=(240, 240, 235))
        dc.rectangle([ox, HEIGHT + 31, ox + INSET - 1, HEIGHT + 64], fill=(12, 18, 28))
        dc.text((ox + 10, HEIGHT + 38), f"{LABEL.get(name[c], name[c])}  (deco_scale={DECO_SCALE})", font=f_lab,
                fill=(245, 245, 240))
        # Mark the inset on the main view.
        px = (x - REGION["w"]) / (REGION["e"] - REGION["w"]) * WIDTH
        py = (REGION["n"] - y) / (REGION["n"] - REGION["s"]) * HEIGHT
        r = h / (REGION["e"] - REGION["w"]) * WIDTH
        dc.rectangle([px - r, py - r, px + r, py + r], outline=(240, 240, 235), width=2)
        dc.text((px + r + 4, py - r), str(k + 1), font=f_lab, fill=(245, 245, 240))
        dc.text((ox + INSET - 24, HEIGHT + 38), str(k + 1), font=f_lab, fill=(245, 245, 240))
    dc.rectangle([0, 0, WIDTH, 52], fill=(12, 18, 28))
    dc.text((18, 12), "i.ocean -w -f -l -d    temperate coast, 100 km x 60 km, 100 m cells;"
            "  insets: vector silhouettes over the raster", font=f_title, fill=(245, 245, 240))
    canvas.save(png)
    tools.g_region(**REGION)
    return picks


def session_run(project_dir, png):
    project = Path(project_dir) / "ocean_example"
    gs.create_project(project, epsg="32630")
    with gs.setup.init(project, env=os.environ.copy()) as session:
        tools = Tools(session=session)
        scene(tools)
        picks = render_example(tools, session.env, png)
        return tools, picks


def test_example_scene(tmp_path):
    png = tmp_path / "example.png"
    project = tmp_path / "ocean_example"
    gs.create_project(project, epsg="32630")
    with gs.setup.init(project, env=os.environ.copy()) as session:
        tools = Tools(session=session)
        scene(tools)
        # The ocean raster covers the sea only.
        sea = tools.r_univar(map="sea", format="json")["n"]
        ocean = tools.r_univar(map="ocean", format="json")
        assert ocean["n"] == sea > 0
        # Decorations: areas with their centroids, lines, a filled attribute table.
        topo = tools.v_info(map="ocean", flags="t", format="json")
        assert topo["areas"] > 0 and topo["centroids"] > 0 and topo["lines"] > 0
        rows = tools.v_db_select(map="ocean", format="json")["records"]
        assert all(r["type"] and r["len_m"] > 0 for r in rows)
        types = {r["cat"]: r["type"] for r in rows}
        assert types[22] == "crude_tanker" and types[26] == "container_ship"
        # Merchant beams are those of real ships (scaled by deco_scale).
        wid = {r["type"]: r["wid_m"] / DECO_SCALE for r in rows}
        assert abs(wid["container_ship"] - 48.0) < 0.01 and abs(wid["crude_tanker"] - 44.0) < 0.01
        # At least two vessels were placed, and the composed image is not blank.
        assert len(vessel_centres(tools)) >= 2
        render_example(tools, session.env, png)
        img = np.asarray(Image.open(png).convert("L"), float)
        assert img.shape == (HEIGHT + INSET + 40, WIDTH)
        assert img.std() > 10


if __name__ == "__main__":
    out = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent.parent / "example.png"
    with tempfile.TemporaryDirectory() as d:
        _, picks = session_run(d, out)
    print(f"Wrote {out} (insets: {[p[2] for p in picks]})")
