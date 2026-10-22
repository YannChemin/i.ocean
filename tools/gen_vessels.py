#!/usr/bin/env python3
"""
MODULE:       i.ocean -- tools/gen_vessels.py
AUTHOR(S):    i.ocean contributors
PURPOSE:      Generate the vessel silhouettes of i.ocean (vector_deco.c) from
              cubic Bezier curves, arcs, ellipses and rounded rectangles.
SPDX-License-Identifier: Unlicense

Each vessel is a list of labelled parts in local metres: +Y towards the bow,
+X to starboard, origin amidships, as in vector_deco.c. Curves are tessellated
to about 1 m chords, so the outlines stay smooth at any map scale. Closed parts
repeat their first vertex exactly (vector_deco.c turns them into areas).

Writes ../vessel_shapes.h (C arrays, included by vector_deco.c) and
../vessel_shapes.json (the same geometry with labels, for other tools).

Usage: python3 tools/gen_vessels.py
"""

import json
import math
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent


# Geometry helpers. All return (n, 2) arrays of (x, y).

def bezier(p0, p1, p2, p3, n=16):
    """Cubic Bezier from p0 to p3 with control points p1, p2."""
    t = np.linspace(0, 1, n)[:, None]
    p0, p1, p2, p3 = (np.asarray(p, float) for p in (p0, p1, p2, p3))
    return (1 - t) ** 3 * p0 + 3 * (1 - t) ** 2 * t * p1 + 3 * (1 - t) * t ** 2 * p2 + t ** 3 * p3


def chain(*parts):
    """Join point runs, dropping the repeated joint vertices."""
    out = [np.asarray(parts[0], float)]
    for p in parts[1:]:
        p = np.asarray(p, float)
        if np.allclose(out[-1][-1], p[0]):
            p = p[1:]
        out.append(p)
    return np.vstack(out)


def close(p):
    p = np.asarray(p, float)
    if not np.allclose(p[0], p[-1]):
        p = np.vstack([p, p[:1]])
    p[-1] = p[0]
    return p


def mirror_hull(half):
    """Closed outline from a starboard half running stern centre -> bow tip."""
    port = half[::-1].copy()
    port[:, 0] *= -1
    return close(chain(half, port[1:]))


def ship_hull(L, B, stern_len, bow_len, transom=0.8, fullness=0.55, bulb=0.0):
    """Merchant or naval hull, top view.

    transom: transom half-width over half-beam; stern_len, bow_len: lengths
    of the stern run and of the bow entry (m); fullness 0 (fine) .. 1 (blunt).
    """
    ys, yb, h = -L / 2, L / 2, B / 2
    tw = transom * h
    r = min(stern_len, h - tw + 1e-6)
    half = chain(
        [[0, ys], [tw, ys]],
        bezier([tw, ys], [tw + 0.55 * (h - tw), ys], [h, ys + 0.45 * r], [h, ys + stern_len], 14),
        [[h, ys + stern_len], [h, yb - bow_len]],
        bezier([h, yb - bow_len], [h, yb - bow_len * (1 - fullness)],
               [h * (0.35 + 0.3 * fullness) + bulb, yb - 0.02 * bow_len], [0, yb], 22),
    )
    return mirror_hull(half)


def pointed_hull(L, B, stern_len, bow_len, stern_tip=0.0):
    """Hull pointed at the bow and narrow at the stern (fast boat, galley)."""
    ys, yb, h = -L / 2, L / 2, B / 2
    half = chain(
        [[0, ys], [stern_tip * h, ys]],
        bezier([stern_tip * h, ys], [h, ys + 0.1 * stern_len], [h, ys + 0.6 * stern_len],
               [h, ys + stern_len], 12),
        [[h, ys + stern_len], [h, yb - bow_len]],
        bezier([h, yb - bow_len], [h, yb - 0.45 * bow_len], [0.35 * h, yb - 0.05 * bow_len], [0, yb], 22),
    )
    return mirror_hull(half)


def myring(L, D, nose=0.18, tail=0.32, n=28):
    """Submarine hull, top view: Myring profile (elliptic nose, cubic tail)."""
    a, c = nose * L, tail * L
    y0 = -L / 2
    r = D / 2
    yn = np.linspace(0, a, n)  # Nose, from its tip.
    rn = r * (1 - ((a - yn) / a) ** 2) ** 0.5
    yt = np.linspace(0, c, n)  # Tail, from the parallel body to the shaft.
    # Convex taper to a shaft of 12 % of the radius (no wine-glass waist).
    rt = r * (0.12 + 0.88 * (1 - (yt / c) ** 2) ** 0.9)
    half = chain(
        np.c_[rt[::-1], y0 + c - yt[::-1]],
        np.c_[[r, r], [y0 + c, y0 + L - a]],
        np.c_[rn[::-1], y0 + L - yn[::-1]],
    )
    half = np.vstack([[0, y0], half[1:]]) if half[0, 0] > 1e-9 else half
    return mirror_hull(half)


def ellipse(cx, cy, rx, ry, n=40):
    t = np.linspace(0, 2 * math.pi, n + 1)
    return close(np.c_[cx + rx * np.cos(t), cy + ry * np.sin(t)])


def circle(cx, cy, r, n=None):
    return ellipse(cx, cy, r, r, n or max(16, min(48, int(2 * math.pi * r))))


def rrect(cx, cy, w, h, r=None, n=5):
    """Rounded rectangle centred at (cx, cy), width w (x) and height h (y)."""
    r = min(r if r is not None else 0.15 * min(w, h), w / 2, h / 2)
    x0, x1, y0, y1 = cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2
    pts = []
    for (ax, ay, a0) in ((x1 - r, y1 - r, 0), (x0 + r, y1 - r, 90), (x0 + r, y0 + r, 180), (x1 - r, y0 + r, 270)):
        t = np.radians(np.linspace(a0, a0 + 90, n))
        pts.append(np.c_[ax + r * np.cos(t), ay + r * np.sin(t)])
    return close(np.vstack(pts))


def line(*pts):
    return np.asarray(pts, float)


def arc(cx, cy, r, a0, a1, n=16):
    """Open arc, angles in degrees from +X (counter-clockwise)."""
    t = np.radians(np.linspace(a0, a1, n))
    return np.c_[cx + r * np.cos(t), cy + r * np.sin(t)]


def x_half_width(outline, y):
    """Half-width of a closed outline at y (largest |x| on the outline near y)."""
    p = outline
    xs = []
    for (x0, y0), (x1, y1) in zip(p[:-1], p[1:]):
        if (y0 - y) * (y1 - y) <= 0 and y0 != y1:
            xs.append(abs(x0 + (y - y0) / (y1 - y0) * (x1 - x0)))
    return max(xs) if xs else 0.0


def deck_break(hull, y, inset=1.0):
    """Transverse deck line (forecastle break, transom line) across the hull at y."""
    w = x_half_width(hull, y) - inset
    return line([-w, y], [w, y])


def superstructure(hull, B, y_aft, depth, width_frac=0.85, funnel=True):
    """Accommodation block (one outline with its bridge wings), funnel, lifeboats.

    The wings run to the ship's side at the forward face of the block; the
    funnel stands on the block; the lifeboats hang beside it. Parts are nested
    or apart, never crossing (GRASS areas).
    """
    w = width_frac * B / 2
    y1 = y_aft + depth
    ww = min(x_half_width(hull, y1 - 1.1), x_half_width(hull, y1 - 2.4)) - 0.6
    r = 1.2
    house = close(chain(
        arc(w - r, y_aft + r, r, 270, 360, 5), [[w, y1 - 2.4], [ww, y1 - 2.4]],
        arc(ww - 0.6, y1 - 0.6, 0.6, 0, 90, 4), arc(-ww + 0.6, y1 - 0.6, 0.6, 90, 180, 4),
        [[-ww, y1 - 2.4], [-w, y1 - 2.4]], arc(-w + r, y_aft + r, r, 180, 270, 5)))
    parts = [("deckhouse", house)]
    if funnel:
        parts.append(("funnel", rrect(0, y_aft + 0.28 * depth, 0.22 * B, 0.3 * depth, 2.0)))
    xb = w + 2.0
    yb = y_aft + 0.35 * depth
    if x_half_width(hull, yb - 4.5) > xb + 2.0 and x_half_width(hull, yb + 4.5) > xb + 2.0:
        parts.append(("lifeboat", rrect(xb, yb, 3.0, 8.5, 1.4)))
        parts.append(("lifeboat", rrect(-xb, yb, 3.0, 8.5, 1.4)))
    return parts


# ── vessels ───────────────────────────────────────────────────────────────

def crude_tanker(L=250.0, B=44.0, n_tanks=6):
    hull = ship_hull(L, B, 0.12 * L, 0.18 * L, transom=0.78, fullness=0.75)
    p = [("hull", hull)]
    y_fc = L / 2 - 0.09 * L
    p.append(("forecastle", deck_break(hull, y_fc)))
    ya, yb = -L / 2 + 0.2 * L, y_fc - 4
    edges = np.linspace(ya, yb, n_tanks + 1)
    for k in range(n_tanks):
        p.append(("tank", rrect(0, 0.5 * (edges[k] + edges[k + 1]), 0.7 * B, edges[k + 1] - edges[k] - 1.0, 1.0)))
    for k in range(n_tanks):
        yc = 0.5 * (edges[k] + edges[k + 1])
        p.append(("tank_hatch", circle(0.22 * B, yc, 1.3)))
        p.append(("tank_hatch", circle(-0.22 * B, yc, 1.3)))
    p.append(("pipe_rack", line([0.6, ya], [0.6, yb])))
    p.append(("pipe_rack", line([-0.6, ya], [-0.6, yb])))
    ym = 0.5 * (ya + yb)
    p.append(("manifold", line([-0.45 * B, ym], [0.45 * B, ym])))
    p.append(("manifold", line([-0.45 * B, ym + 1.5], [0.45 * B, ym + 1.5])))
    p.append(("crane_base", circle(0.3 * B, ym + 6, 1.2)))
    p.append(("crane", line([0.3 * B, ym + 6], [0.1 * B, ym + 20])))
    p.append(("crane_base", circle(-0.3 * B, ym - 6, 1.2)))
    p.append(("crane", line([-0.3 * B, ym - 6], [-0.1 * B, ym - 20])))
    p += superstructure(hull, B, -L / 2 + 0.06 * L, 0.11 * L)
    p.append(("mooring", arc(0, L / 2 - 0.13 * L, 0.25 * B, 25, 155)))
    return p


def product_tanker(L=180.0, B=32.0):
    return crude_tanker(L, B, n_tanks=6)


def lng_carrier(L=290.0, B=49.0):
    hull = ship_hull(L, B, 0.13 * L, 0.2 * L, transom=0.75, fullness=0.6)
    p = [("hull", hull)]
    y_fc = L / 2 - 0.1 * L
    p.append(("forecastle", deck_break(hull, y_fc)))
    r = 0.43 * B
    ys = np.linspace(-L / 2 + 0.25 * L + r, y_fc - r - 3, 4)
    for y in ys:
        p.append(("sphere", circle(0, y, r, 64)))
    for y in ys:
        p.append(("sphere_cap", circle(0, y, 3.0)))
    for a, b in zip(ys[:-1], ys[1:]):
        p.append(("walkway", line([0.47 * B - 1, a], [0.47 * B - 1, b])))
    p += superstructure(hull, B, -L / 2 + 0.06 * L, 0.12 * L)
    p.append(("mooring", arc(0, L / 2 - 0.13 * L, 0.25 * B, 25, 155)))
    return p


def bulk_carrier(L=300.0, B=50.0, n_holds=9):
    hull = ship_hull(L, B, 0.11 * L, 0.16 * L, transom=0.8, fullness=0.8)
    p = [("hull", hull)]
    y_fc = L / 2 - 0.08 * L
    p.append(("forecastle", deck_break(hull, y_fc)))
    ya, yb = -L / 2 + 0.2 * L, y_fc - 3
    edges = np.linspace(ya, yb, n_holds + 1)
    for k in range(n_holds):
        h = edges[k + 1] - edges[k]
        p.append(("hatch", rrect(0, 0.5 * (edges[k] + edges[k + 1]), 0.6 * B, 0.72 * h, 1.2)))
    for k in range(n_holds):
        p.append(("hatch_seam", line([0, edges[k] + 0.14 * (edges[k + 1] - edges[k])],
                                     [0, edges[k + 1] - 0.14 * (edges[k + 1] - edges[k])])))
    p += superstructure(hull, B, -L / 2 + 0.05 * L, 0.1 * L)
    p.append(("mooring", arc(0, L / 2 - 0.12 * L, 0.25 * B, 25, 155)))
    return p


def container_ship(L=300.0, B=48.0, n_bays=10):
    hull = ship_hull(L, B, 0.11 * L, 0.21 * L, transom=0.85, fullness=0.45)
    p = [("hull", hull)]
    y_fc = L / 2 - 0.08 * L
    p.append(("forecastle", deck_break(hull, y_fc)))
    y_house = -L / 2 + 0.11 * L
    ya, yb = y_house + 0.06 * L + 4, y_fc - 3
    edges = np.linspace(ya, yb, n_bays + 1)
    for k in range(n_bays):
        yc = 0.5 * (edges[k] + edges[k + 1])
        w = min(0.92 * B, 2 * x_half_width(hull, edges[k + 1]) - 2.5)
        p.append(("bay", rrect(0, yc, w, edges[k + 1] - edges[k] - 1.6, 0.6)))
    for k in range(1, n_bays):
        w = x_half_width(hull, edges[k]) - 1.0
        p.append(("lashing_bridge", line([-w, edges[k]], [w, edges[k]])))
    p += superstructure(hull, B, y_house, 0.06 * L, funnel=False)
    p.append(("funnel", rrect(0, -L / 2 + 0.065 * L, 0.25 * B, 0.03 * L, 2.0)))
    p.append(("mooring", arc(0, L / 2 - 0.12 * L, 0.2 * B, 25, 155)))
    return p


def cargo_ship(L=140.0, B=22.0):
    hull = ship_hull(L, B, 0.13 * L, 0.2 * L, transom=0.75, fullness=0.5)
    p = [("hull", hull)]
    y_fc = L / 2 - 0.1 * L
    p.append(("forecastle", deck_break(hull, y_fc)))
    ya, yb = -L / 2 + 0.24 * L, y_fc - 3
    edges = np.linspace(ya, yb, 4)
    for k in range(3):
        p.append(("hatch", rrect(0, 0.5 * (edges[k] + edges[k + 1]), 0.7 * B, edges[k + 1] - edges[k] - 6, 1.0)))
    for k, y in enumerate(edges[:3] + 1.5):
        side = 1 if k % 2 == 0 else -1
        p.append(("crane_base", circle(side * 0.3 * B, y, 1.2)))
        p.append(("crane", line([side * 0.3 * B, y], [-side * 0.15 * B, y + 0.12 * L])))
    p += superstructure(hull, B, -L / 2 + 0.06 * L, 0.13 * L)
    p.append(("mooring", arc(0, L / 2 - 0.13 * L, 0.25 * B, 25, 155)))
    return p


def hunter_sub(L=78.0, D=8.0):
    p = [("hull", myring(L, D))]
    ys = 0.14 * L
    p.append(("sail", rrect(0, ys, 0.4 * D, 0.13 * L, 0.18 * D)))
    p.append(("sail_plane", line([-0.75 * D, ys + 0.03 * L], [0.75 * D, ys + 0.03 * L])))
    yt = -L / 2 + 0.05 * L
    p.append(("stern_plane", line([-0.9 * D, yt], [0.9 * D, yt])))
    p.append(("rudder", line([0, yt - 0.04 * L], [0, yt + 0.07 * L])))
    p.append(("hatch", circle(0, 0.33 * L, 0.6)))
    p.append(("hatch", circle(0, -0.2 * L, 0.6)))
    p.append(("dome", ellipse(0, L / 2 - 0.06 * L, 0.3 * D, 0.05 * L)))
    return p


def barracuda_boat(L=50.0, B=6.0):
    p = [("hull", pointed_hull(L, B, 0.1 * L, 0.32 * L, stern_tip=0.85))]
    ys = -L / 2 + 0.2 * L
    for s in (1, -1):
        o = s * 0.15  # Clear of the hull side.
        fin = close(chain(
            bezier([s * B / 2 + o, ys + 3], [s * (B / 2 + 1.5) + o, ys + 1.5], [s * (B / 2 + 3) + o, ys - 1], [s * B + o, ys - 3], 8),
            bezier([s * B + o, ys - 3], [s * (B - 0.3) + o, ys - 4], [s * (B / 2 + 1) + o, ys - 3.5],
                   [s * B / 2 + o, ys - 2.5], 8)))
        p.append(("fin", fin))
    p.append(("superstructure", rrect(0, 0.0, 0.7 * B, 0.28 * L, 1.2)))
    p.append(("mast", circle(0, 0.06 * L, 0.6)))
    p.append(("gun", circle(0, 0.3 * L, 1.0)))
    p.append(("gun", line([0, 0.3 * L], [0, 0.3 * L + 3.5])))
    return p


def sunken_ship(L=37.0, B=7.0):
    p = [("hull", pointed_hull(L - 3, B, 0.3 * L, 0.35 * L, stern_tip=0.15))]
    yb = (L - 3) / 2
    p.append(("ram", close(np.array([[0, yb + 0.1], [0.7, yb + 1.0], [0, yb + 3.0], [-0.7, yb + 1.0], [0, yb + 0.1]]))))
    for s in (1, -1):
        for y in np.linspace(-0.3 * L, 0.3 * L, 9):
            w = x_half_width(p[0][1], y)
            p.append(("oar", line([s * w, y], [s * (w + 2.5), y - 1.0])))
    p.append(("mast", circle(0, 0.05 * L, 0.4)))
    p.append(("deck_beam", deck_break(p[0][1], -0.22 * L, 0.3)))
    return p


def _spindle(L, B, nose=0.2, tail=0.3):
    return myring(L, B, nose, tail)


def leviathan_sub():
    L, B = 95.0, 20.0
    p = [("hull", _spindle(L, B, 0.22, 0.3))]
    p += [("silo", circle(-0.25 * B, 0.12 * L, 1.6)), ("silo", circle(0.25 * B, 0.12 * L, 1.6)),
          ("silo", circle(-0.25 * B, 0.2 * L, 1.6)), ("silo", circle(0.25 * B, 0.2 * L, 1.6))]
    p.append(("reactor", rrect(0, -0.18 * L, 0.5 * B, 0.15 * L, 2.0)))
    p.append(("sail", rrect(0, 0.3 * L, 0.18 * B, 0.12 * L, 1.5)))
    return p


def kraken_carrier():
    L, B = 65.0, 18.0
    p = [("hull", _spindle(L, B, 0.2, 0.25))]
    p.append(("hangar", rrect(0, 0.08 * L, 0.85 * B, 0.32 * L, 2.0)))
    p.append(("rail", rrect(-0.32 * B, 0.08 * L, 0.1 * B, 0.28 * L, 0.6)))
    p.append(("rail", rrect(0.32 * B, 0.08 * L, 0.1 * B, 0.28 * L, 0.6)))
    p.append(("sail", rrect(0, -0.2 * L, 0.2 * B, 0.1 * L, 1.0)))
    return p


def phantom_sub():
    L, B = 60.0, 10.0
    p = [("hull", _spindle(L, B, 0.25, 0.35))]
    p.append(("core", ellipse(0, 0.0, 0.3 * B, 0.17 * L)))
    return p


def titan_barge():
    L, B = 47.0, 16.0
    # Blunt, nearly flat bow and a long rounded run to a narrow transom.
    p = [("hull", ship_hull(L, B, 0.3 * L, 0.05 * L, transom=0.3, fullness=0.95))]
    for y in (0.4 * L, 0.3 * L, 0.2 * L):
        for x in np.linspace(-0.3 * B, 0.3 * B, 4):
            p.append(("tube", circle(x, y, 0.8, 12)))
    p.append(("ammo_bay", rrect(0, -0.05 * L, 0.6 * B, 0.3 * L, 1.5)))
    return p


def wraith_sub():
    L, B = 27.0, 5.0
    return [("hull", _spindle(L, B, 0.25, 0.35)), ("hatch", circle(0, 0.05 * L, 0.9))]


def hydra_sub():
    L, B = 55.0, 10.0
    body = _spindle(0.72 * L, B, 0.15, 0.35)
    body[:, 1] -= 0.14 * L
    p = [("hull", body)]
    for s in (1, -1):
        prong = chain(bezier([s * 0.2 * B, 0.2 * L], [s * 0.3 * B, 0.3 * L], [s * 0.55 * B, 0.38 * L], [s * 0.5 * B, 0.45 * L], 12))
        p.append(("prong", prong))
        p.append(("weapon_node", circle(s * 0.5 * B, 0.45 * L, 1.0)))
    return p


def aegis_sub():
    L, B = 53.0, 14.0
    p = [("hull", _spindle(L, B, 0.2, 0.3))]
    for s in (1, -1):
        p.append(("sensor_arm", rrect(s * (B / 2 + 2.8), 0.08 * L, 5.0, 2.0, 0.8)))
    p.append(("radar_core", circle(0, 0.08 * L, 0.25 * B)))
    return p


def oblivion_sub():
    L, B = 115.0, 24.0
    p = [("hull", _spindle(L, B, 0.2, 0.3))]
    p.append(("superstructure", rrect(0, 0.05 * L, 0.5 * B, 0.26 * L, 2.5)))
    p.append(("bridge", rrect(0, 0.23 * L, 0.3 * B, 0.07 * L, 1.5)))
    p.append(("reactor", ellipse(0, -0.26 * L, 0.18 * B, 0.09 * L)))
    return p


VESSELS = {
    "sunken_ship": sunken_ship, "hunter_sub": hunter_sub, "leviathan_sub": leviathan_sub,
    "barracuda_boat": barracuda_boat, "kraken_carrier": kraken_carrier, "phantom_sub": phantom_sub,
    "titan_barge": titan_barge, "wraith_sub": wraith_sub, "hydra_sub": hydra_sub,
    "aegis_sub": aegis_sub, "oblivion_sub": oblivion_sub,
    "crude_tanker": crude_tanker, "lng_carrier": lng_carrier, "bulk_carrier": bulk_carrier,
    "product_tanker": product_tanker, "container_ship": container_ship, "cargo_ship": cargo_ship,
}


def _dist_to_edges(P, poly):
    a, b = poly[:-1], poly[1:]
    d = b - a
    t = np.clip(((P[:, None, :] - a[None]) * d[None]).sum(2) / np.maximum((d * d).sum(1), 1e-12), 0, 1)
    q = a[None] + t[..., None] * d[None]
    return np.sqrt(((P[:, None, :] - q) ** 2).sum(2)).min(1)


def _area(poly):
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * abs((x[:-1] * y[1:] - x[1:] * y[:-1]).sum())


def label_points(parts):
    """A point inside each closed part and outside the closed parts nested in it.

    GRASS puts the centroid of an area there; a centroid at the vertex average
    of a hull would fall inside a tank or hatch drawn on the deck and leave the
    deck area without one (a duplicate centroid).
    """
    from matplotlib.path import Path as MPath
    closed = [(k, p) for k, (_, p) in enumerate(parts) if len(p) >= 4 and (p[0] == p[-1]).all()]
    out = {}
    for k, p in closed:
        x0, y0 = p.min(0)
        x1, y1 = p.max(0)
        step = max(min(x1 - x0, y1 - y0) / 40, 0.05)
        gx, gy = np.meshgrid(np.arange(x0 + step / 2, x1, step), np.arange(y0 + step / 2, y1, step))
        P = np.c_[gx.ravel(), gy.ravel()]
        ok = MPath(p).contains_points(P)
        others = [q for j, q in closed if j != k and _area(q) < _area(p)]
        for q in others:
            ok &= ~MPath(q).contains_points(P)
        if not ok.any():
            out[k] = tuple(p[:-1].mean(0))
            continue
        P = P[ok]
        d = _dist_to_edges(P, p)
        for q in others:
            d = np.minimum(d, _dist_to_edges(P, q))
        out[k] = tuple(P[int(np.argmax(d))])
    return out


def crossings(parts):
    """Pairs of closed parts whose boundaries cross (GRASS needs them nested or apart)."""
    closed = [(lab, p) for lab, p in parts if len(p) >= 4 and (p[0] == p[-1]).all()]
    bad = []
    for i in range(len(closed)):
        for j in range(i + 1, len(closed)):
            a, b = closed[i][1], closed[j][1]
            if (a[:, 0].max() < b[:, 0].min() or b[:, 0].max() < a[:, 0].min() or
                    a[:, 1].max() < b[:, 1].min() or b[:, 1].max() < a[:, 1].min()):
                continue
            p1, p2 = a[:-1][:, None], a[1:][:, None]
            q1, q2 = b[:-1][None], b[1:][None]

            def orient(u, v, w):
                return np.sign((v[..., 0] - u[..., 0]) * (w[..., 1] - u[..., 1]) -
                               (v[..., 1] - u[..., 1]) * (w[..., 0] - u[..., 0]))
            hit = ((orient(p1, p2, q1) * orient(p1, p2, q2) <= 0) &
                   (orient(q1, q2, p1) * orient(q1, q2, p2) <= 0))
            if hit.any():
                bad.append((closed[i][0], closed[j][0]))
    return bad


def fmt(v):
    return f"{v:.3f}".rstrip("0").rstrip(".") if abs(v) > 5e-4 else "0"


def main():
    shapes = {}
    for name, fn in VESSELS.items():
        parts = []
        for label, pts in fn():
            pts = np.round(np.asarray(pts, float), 3)
            if len(pts) >= 4 and np.allclose(pts[0], pts[-1]):
                pts[-1] = pts[0]
            parts.append((label, pts))
        shapes[name] = parts

    js = {}
    for n, parts in shapes.items():
        lp = label_points(parts)
        js[n] = [{"label": lab, "closed": bool(len(p) >= 4 and (p[0] == p[-1]).all()),
                  "label_point": [float(v) for v in lp[k]] if k in lp else None,
                  "pts": p.tolist()} for k, (lab, p) in enumerate(parts)]
    (ROOT / "vessel_shapes.json").write_text(json.dumps(js, separators=(",", ":")))

    out = ["/* Vessel silhouettes of i.ocean, in local metres (+Y bow, +X starboard,",
           " * origin amidships). Generated by tools/gen_vessels.py: do not edit by hand.",
           " * SPDX-License-Identifier: Unlicense */", "",
           "#ifndef IOCEAN_VESSEL_SHAPES_H", "#define IOCEAN_VESSEL_SHAPES_H", "",
           "typedef struct {", "    const char *label;     /* part name */",
           "    int n;                 /* number of points */",
           "    const double (*pts)[2]; /* (x, y) in metres; closed if first == last */",
           "    double lpt[2];         /* label (centroid) point of a closed part */",
           "} VPrim;", "", "typedef struct {", "    const char *name;      /* deco_type_name */",
           "    int n_prims;", "    const VPrim *prims;", "} VShape;", ""]
    for name, parts in shapes.items():
        lp = label_points(parts)
        for k, (lab, p) in enumerate(parts):
            body = ",".join("{" + fmt(x) + "," + fmt(y) + "}" for x, y in p)
            out.append(f"static const double v_{name}_{k}[][2] = {{{body}}};")
        out.append(f"static const VPrim v_{name}_prims[] = {{")
        out += [f'    {{"{lab}", {len(p)}, v_{name}_{k}, {{{fmt(lp.get(k, (0, 0))[0])}, {fmt(lp.get(k, (0, 0))[1])}}}}},'
                for k, (lab, p) in enumerate(parts)]
        out.append("};")
        out.append(f'static const VShape vshape_{name} = {{"{name}", {len(parts)}, v_{name}_prims}};')
        out.append("")
    out.append("static const VShape *const vessel_shapes[] = {")
    out += [f"    &vshape_{n}," for n in shapes]
    out.append("};")
    out.append(f"#define N_VESSEL_SHAPES {len(shapes)}")
    out += ["", "#endif", ""]
    (ROOT / "vessel_shapes.h").write_text("\n".join(out))
    for n, parts in shapes.items():
        pts = sum(len(p) for _, p in parts)
        bad = crossings(parts)
        print(f"{n:16s} {len(parts):3d} parts {pts:5d} points" +
              (f"  CROSSING: {sorted(set(bad))}" if bad else ""))


if __name__ == "__main__":
    main()
