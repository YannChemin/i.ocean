#!/usr/bin/env python3

############################################################################
#
# MODULE:       i.ocean
# AUTHOR(S):    i.ocean contributors
# PURPOSE:      Render ocean areas of a raster map with a realistic
#               appearance, using depth, CRS scale, and latitude to drive
#               colour palettes, depth gradients, and wave textures.
# COPYRIGHT:    (C) 2026 by the GRASS Development Team
#               This program is free software under the GNU General Public
#               License (>=v2). Read the file COPYING that comes with GRASS
#               for details.
#
############################################################################

# %module
# % description: Renders ocean areas of a raster map with a realistic visual appearance
# % keyword: raster
# % keyword: ocean
# % keyword: visualization
# % keyword: color
# % keyword: depth
# %end

# %option G_OPT_R_INPUT
# % key: input
# % label: Ocean mask raster
# % description: Raster whose cells equal to ocean_value are treated as ocean; mutually optional with depth
# % required: no
# %end

# %option G_OPT_R_OUTPUT
# % key: output
# % label: Output ocean visualization raster
# % description: Floating-point raster carrying depth-index values and ocean colour rules
# % required: yes
# %end

# %option G_OPT_R_INPUT
# % key: depth
# % label: Bathymetric depth raster (positive values = metres below surface)
# % description: Drives the colour gradient; also used as the ocean source when no input mask is given
# % required: no
# %end

# %option
# % key: ocean_value
# % type: double
# % label: Cell value in the input mask that represents ocean
# % description: Cells with this value are treated as ocean; ignored when only depth is supplied
# % answer: 1
# % required: no
# %end

# %option
# % key: depth_min
# % type: double
# % label: Minimum depth value to treat as ocean (used when deriving the mask from the depth map)
# % description: Cells in the depth map with value >= depth_min become ocean; use 0 for GEBCO-style maps, a negative value for elevation/DEM maps
# % answer: 0
# % required: no
# %end

# %rules
# % required: input, depth
# %end

# %option
# % key: latitude
# % type: double
# % label: Reference latitude for colour-temperature palette (decimal degrees, –90 to 90)
# % description: Auto-detected from the map centre when omitted
# % required: no
# %end

# %option
# % key: style
# % type: string
# % label: Ocean colour style
# % description: Overrides the auto-detected latitude zone
# % answer: auto
# % options: auto,tropical,subtropical,temperate,subpolar,polar
# % descriptions: auto;Derive style from latitude;tropical;Warm turquoise (< 23.5°);subtropical;Blue-green (23.5–35°);temperate;Cool blue (35–60°);subpolar;Steel blue (60–75°);polar;Icy pale blue (> 75°)
# % required: no
# % guisection: Palette
# %end

# %flag
# % key: w
# % label: Add wave / ripple texture
# % description: Overlays sinusoidal wave patterns scaled to the current pixel size
# % guisection: Effects
# %end

# %flag
# % key: s
# % label: Simulate depth from distance to shore
# % description: When no depth map is given, estimates depth using r.grow.distance from land pixels
# % guisection: Effects
# %end

# %flag
# % key: l
# % label: Apply subtle north–south warmth gradient across the map
# % description: Mimics the visible SST gradient within the scene latitude band
# % guisection: Effects
# %end

import atexit
import math
import os
import sys

import grass.script as gs

# ── temporary map registry ──────────────────────────────────────────────────
_TMP = []


def _register(*names):
    """Register one or more temporary raster names for removal on exit."""
    _TMP.extend(names)
    return names[-1] if len(names) == 1 else names


def _cleanup():
    """Remove all registered temporary raster maps created during the run."""
    if _TMP:
        gs.run_command(
            "g.remove",
            flags="f",
            type="raster",
            name=_TMP,
            quiet=True,
            errors="ignore",
        )


atexit.register(_cleanup)


# ── CRS / scale helpers ──────────────────────────────────────────────────────

def _pixel_size_m(region, proj):
    """Return the mean pixel size in metres for the current computational region.

    For geographic CRS (degrees) the east–west pixel size is scaled by the
    cosine of the centre latitude before averaging with the north–south size.
    For projected CRS the resolution is converted from its native unit (feet,
    kilometres, or metres) to metres.
    """
    units = proj.get("units", "").lower()
    proj_name = proj.get("proj", "")

    if "degree" in units or proj_name in ("ll", "latlong", "longlat"):
        lat_rad = math.radians((region["n"] + region["s"]) / 2.0)
        ns_m = region["nsres"] * 111_320.0
        ew_m = region["ewres"] * 111_320.0 * math.cos(lat_rad)
        return (ns_m + ew_m) / 2.0

    scale = 1.0
    if "feet" in units or "foot" in units:
        scale = 0.3048
    elif "kilometer" in units or units == "km":
        scale = 1_000.0
    return ((region["nsres"] + region["ewres"]) / 2.0) * scale


def _center_latitude(region, proj):
    """Return the geographic latitude of the computational region centre in decimal degrees.

    For geographic CRS the latitude is taken directly from the region bounds.
    For projected CRS the centre point is reprojected to geographic coordinates
    using *m.proj*. If reprojection fails, 45°N is returned and a warning is
    issued.
    """
    units = proj.get("units", "").lower()
    proj_name = proj.get("proj", "")

    if "degree" in units or proj_name in ("ll", "latlong", "longlat"):
        return (region["n"] + region["s"]) / 2.0

    # Projected: reproject the centre point to geographic
    cy = (region["n"] + region["s"]) / 2.0
    cx = (region["e"] + region["w"]) / 2.0
    out = gs.read_command(
        "m.proj", coordinates=f"{cx},{cy}", flags="od"
    ).strip()
    if out:
        parts = out.split("|")
        if len(parts) >= 2:
            try:
                return float(parts[1])
            except ValueError:
                pass
    gs.warning(_("Could not detect latitude from CRS; assuming 45°N (temperate)."))
    return 45.0


# ── style / palette ──────────────────────────────────────────────────────────

def _auto_style(lat):
    """Return the colour style name that corresponds to the given latitude.

    Latitude thresholds follow standard oceanographic climate zones:
    tropical (< 23.5°), subtropical (< 35°), temperate (< 60°),
    subpolar (< 75°), polar (≥ 75°). The absolute value of the latitude
    is used so that the function is symmetric about the equator.
    """
    a = abs(lat)
    if a < 23.5:
        return "tropical"
    if a < 35.0:
        return "subtropical"
    if a < 60.0:
        return "temperate"
    if a < 75.0:
        return "subpolar"
    return "polar"


# Each palette is a list of (depth_index, R, G, B).
# depth_index runs 0 (shore) → 1000 (abyss).
_PALETTES = {
    "tropical": [
        (0,    0, 230, 255),   # bright turquoise
        (60,   0, 200, 225),   # clear cyan
        (120,  0, 165, 185),   # lagoon blue
        (220,  0, 125, 165),   # medium blue
        (420,  0,  82, 144),   # ocean blue
        (700,  0,  50, 100),   # deep blue
        (1000, 0,  20,  58),   # abyss
    ],
    "subtropical": [
        (0,    25, 210, 230),
        (60,   10, 180, 210),
        (120,   0, 145, 178),
        (220,   0, 105, 158),
        (420,   0,  68, 132),
        (700,   0,  42,  92),
        (1000,  0,  16,  56),
    ],
    "temperate": [
        (0,    65, 190, 218),
        (60,   38, 160, 198),
        (120,  18, 126, 172),
        (220,   5,  92, 152),
        (420,   0,  60, 122),
        (700,   0,  36,  86),
        (1000,  0,  14,  50),
    ],
    "subpolar": [
        (0,   125, 195, 218),
        (60,   88, 165, 198),
        (120,  58, 136, 177),
        (220,  32, 102, 152),
        (420,  14,  70, 122),
        (700,   5,  44,  86),
        (1000,  0,  20,  52),
    ],
    "polar": [
        (0,   178, 215, 232),  # icy pale blue
        (60,  140, 186, 212),
        (120,  98, 156, 192),
        (220,  62, 120, 166),
        (420,  32,  84, 136),
        (700,  14,  53,  96),
        (1000,  4,  22,  60),
    ],
}


def _build_color_rules(palette, wave_amp=0):
    """Return a GRASS *r.colors* rules string for the given depth-index palette.

    The depth index runs from 0 (shoreline) to 1000 (abyssal plain). When wave
    texture is active the actual raster values can fall slightly outside this
    range, so clamp entries at -(wave_amp + 1) and (1000 + wave_amp + 1) are
    added to prevent colour wrap-around at the extremes. Null cells are mapped
    to white (255:255:255) so that they are transparent when composited over a
    land base map.
    """
    lines = []
    lo_r, lo_g, lo_b = palette[0][1:]
    hi_r, hi_g, hi_b = palette[-1][1:]
    # Values below 0 (wave troughs near shore): clamp to shore colour
    lines.append(f"-{wave_amp + 1} {lo_r}:{lo_g}:{lo_b}")
    for v, r, g, b in palette:
        lines.append(f"{v} {r}:{g}:{b}")
    # Values above 1000 (deep + wave crests): clamp to abyss colour
    lines.append(f"{1000 + wave_amp + 1} {hi_r}:{hi_g}:{hi_b}")
    lines.append("nv 255:255:255")
    return "\n".join(lines)


# ── raster-building helpers ───────────────────────────────────────────────────

def _ocean_mask_from_input(input_map, value, pid):
    """Create a binary ocean mask by selecting cells equal to a given value.

    Cells in *input_map* whose value equals *value* are set to 1; all other
    cells become null. The result is a temporary raster used as the ocean
    extent throughout the rest of the processing chain.
    """
    name = _register(f"tmp_iocean_mask_{pid}")
    gs.mapcalc(
        f"{name} = if({input_map} == {value}, 1, null())",
        overwrite=True, quiet=True,
    )
    return name


def _ocean_mask_from_depth(depth_map, depth_min, pid):
    """Create a binary ocean mask by thresholding a depth raster.

    Cells in *depth_map* with a value greater than or equal to *depth_min* are
    set to 1; all others become null. This allows the depth map to serve as
    the sole input when no separate ocean mask is available. Use *depth_min=0*
    for GEBCO-style maps (positive values = depth below surface) or a negative
    threshold for DEM / elevation maps (negative values = below sea level).
    """
    name = _register(f"tmp_iocean_mask_{pid}")
    gs.mapcalc(
        f"{name} = if({depth_map} >= {depth_min}, 1, null())",
        overwrite=True, quiet=True,
    )
    return name


def _depth_from_map(depth_map, mask, pid):
    """Normalise a bathymetric depth map to the 0–1000 depth index scale.

    Values are taken as absolute depths (negative signs are ignored) and
    linearly rescaled so that the deepest cell in the ocean mask maps to 1000
    and the shallowest maps to 0. If the maximum depth is zero or missing
    a flat index of 500 is used and a warning is issued.
    """
    raw = _register(f"tmp_iocean_depth_raw_{pid}")
    gs.mapcalc(
        f"{raw} = if(isnull({mask}), null(), abs({depth_map}))",
        overwrite=True, quiet=True,
    )
    stats = gs.parse_command("r.univar", map=raw, flags="g", quiet=True)
    max_d = float(stats.get("max", 0) or 0)
    if max_d <= 0:
        gs.warning(_("Depth map maximum is 0 or invalid; using flat depth of 500."))
        max_d = None

    out = _register(f"tmp_iocean_depth_idx_{pid}")
    if max_d:
        gs.mapcalc(
            f"{out} = if(isnull({mask}), null(), "
            f"min(1000.0, 1000.0 * {raw} / {max_d}))",
            overwrite=True, quiet=True,
        )
    else:
        gs.mapcalc(
            f"{out} = if(isnull({mask}), null(), 500.0)",
            overwrite=True, quiet=True,
        )
    return out


def _depth_from_shore(mask, pid):
    """Estimate a depth index from each ocean pixel's distance to the nearest land pixel.

    A temporary land mask (inverse of the ocean mask) is passed to
    *r.grow.distance*, which computes the Euclidean distance in map units to
    the nearest non-null land cell. The resulting distances are linearly
    rescaled to the 0–1000 depth index range, approximating the continental
    shelf gradient without requiring external bathymetric data.
    """
    land = _register(f"tmp_iocean_land_{pid}")
    gs.mapcalc(
        f"{land} = if(isnull({mask}), 1, null())",
        overwrite=True, quiet=True,
    )
    dist = _register(f"tmp_iocean_dist_{pid}")
    gs.run_command(
        "r.grow.distance", input=land, distance=dist,
        overwrite=True, quiet=True,
    )
    stats = gs.parse_command("r.univar", map=dist, flags="g", quiet=True)
    max_dist = float(stats.get("max", 1) or 1)

    out = _register(f"tmp_iocean_depth_idx_{pid}")
    gs.mapcalc(
        f"{out} = if(isnull({mask}), null(), "
        f"min(1000.0, 1000.0 * {dist} / {max_dist}))",
        overwrite=True, quiet=True,
    )
    return out


def _flat_depth(mask, pid, value=500.0):
    """Create a constant depth index raster over all ocean pixels.

    Used as the fallback when neither a depth map nor the **-s** flag is
    provided. The default value of 500 places the colour at the mid-point of
    the palette, giving a uniform mid-ocean appearance.
    """
    out = _register(f"tmp_iocean_depth_flat_{pid}")
    gs.mapcalc(
        f"{out} = if(isnull({mask}), null(), {value})",
        overwrite=True, quiet=True,
    )
    return out


def _wave_pattern(mask, pixel_m, pid):
    """Generate a sinusoidal wave texture calibrated to the current pixel size.

    The texture is a superposition of a primary cross-wave component
    (sin × cos) and a 40 % diagonal harmonic to break the grid regularity.
    Wave period and amplitude are chosen from five scale tiers based on
    *pixel_m* (the mean pixel size in metres):

    - < 100 m: fine ripples (~30 m period, amplitude 12)
    - 100 m – 1 km: ocean swells (~600 m period, amplitude 25)
    - 1 – 10 km: large swells / mesoscale eddies (~6 km period, amplitude 35)
    - 10 – 100 km: broad surface patterns (~60 km period, amplitude 30)
    - > 100 km: large gyre-like undulations (10 px period, amplitude 20)

    Returns a tuple (raster_name, amplitude) where raster values range
    approximately from −amplitude to +amplitude.
    """
    # Choose wave parameters by scale
    if pixel_m < 100:
        # sub-100 m: fine ripples a few metres wide
        period_c = max(4, int(30 / pixel_m))
        period_r = max(4, int(20 / pixel_m))
        amp = 12
    elif pixel_m < 1_000:
        # 100 m – 1 km: visible ocean swells (hundreds of metres)
        period_c = max(5, int(600 / pixel_m))
        period_r = max(5, int(400 / pixel_m))
        amp = 25
    elif pixel_m < 10_000:
        # 1–10 km: large swells / mesoscale eddies visible
        period_c = max(5, int(6_000 / pixel_m))
        period_r = max(5, int(4_000 / pixel_m))
        amp = 35
    elif pixel_m < 100_000:
        # 10–100 km: broad surface patterns
        period_c = max(4, int(60_000 / pixel_m))
        period_r = max(4, int(40_000 / pixel_m))
        amp = 30
    else:
        # > 100 km (regional/global): large gyre-like undulations
        period_c = 10
        period_r = 7
        amp = 20

    two_pi = 2.0 * math.pi
    # Primary wave + 40 % diagonal harmonic for realism
    diag_period = int(period_c * 1.6)
    out = _register(f"tmp_iocean_wave_{pid}")
    gs.mapcalc(
        f"{out} = if(isnull({mask}), null(), "
        f"{amp:.1f} * sin({two_pi:.6f} * col() / {period_c}) "
        f"* cos({two_pi:.6f} * row() / {period_r}) "
        f"+ {amp * 0.4:.1f} * sin({two_pi:.6f} * (col() + row()) / {diag_period}))",
        overwrite=True, quiet=True,
    )
    return out, amp


def _latitude_gradient(depth_idx, lat, mask, pid):
    """Apply a north–south sea-surface-temperature gradient to the depth index.

    Equatorward pixels receive a lower depth-index value (lighter, warmer
    colour) while poleward pixels receive a higher value (darker, cooler
    colour). The total shift is ±80 depth-index units across the full height
    of the computational region, which is sufficient to introduce a visible
    warmth gradient without overpowering the depth signal.

    The sign of the correction is reversed in the southern hemisphere so that
    the gradient always points towards the equator.
    """
    out = _register(f"tmp_iocean_latgrad_{pid}")
    shift = 80.0
    # In the northern hemisphere the south edge is warmer (subtract when row grows)
    # In the southern hemisphere the north edge is warmer (add when row grows)
    if lat >= 0:
        expr = f"{depth_idx} - {shift} * row() / nrows()"
    else:
        expr = f"{depth_idx} + {shift} * row() / nrows()"
    gs.mapcalc(
        f"{out} = if(isnull({mask}), null(), {expr})",
        overwrite=True, quiet=True,
    )
    return out


# ── main ─────────────────────────────────────────────────────────────────────

def main():
    options, flags = gs.parser()

    input_map   = options["input"]
    output      = options["output"]
    depth_map   = options["depth"]
    ocean_val   = options["ocean_value"]
    depth_min   = options["depth_min"]
    lat_opt     = options["latitude"]
    style_opt   = options["style"]
    flag_waves  = flags["w"]
    flag_shore  = flags["s"]
    flag_latgrd = flags["l"]

    pid = os.getpid()

    # ── validate inputs ──────────────────────────────────────────────────────
    if input_map and not gs.find_file(input_map, element="raster")["file"]:
        gs.fatal(_(f"Raster map <{input_map}> not found."))
    if depth_map and not gs.find_file(depth_map, element="raster")["file"]:
        gs.fatal(_(f"Depth raster <{depth_map}> not found."))

    # ── region / projection ──────────────────────────────────────────────────
    gs.message(_("Analysing region and CRS…"))
    region = gs.region()
    proj   = gs.parse_command("g.proj", flags="g")

    pixel_m = _pixel_size_m(region, proj)
    gs.verbose(_(f"Mean pixel size: {pixel_m:,.1f} m"))

    if lat_opt:
        latitude = float(lat_opt)
    else:
        latitude = _center_latitude(region, proj)
    gs.verbose(_(f"Reference latitude: {latitude:.2f}°"))

    style = style_opt if style_opt != "auto" else _auto_style(latitude)
    gs.message(_(f"Ocean colour style: {style}  (latitude {latitude:.1f}°)"))

    # ── ocean mask ───────────────────────────────────────────────────────────
    gs.message(_("Extracting ocean pixels…"))
    if input_map:
        # Explicit mask: a specific cell value marks ocean
        mask = _ocean_mask_from_input(input_map, ocean_val, pid)
        stats = gs.parse_command("r.univar", map=mask, flags="g", quiet=True)
        if int(float(stats.get("n", 0))) == 0:
            gs.fatal(
                _(
                    f"No pixels with value {ocean_val} found in <{input_map}>. "
                    "Adjust ocean_value or check the input map."
                )
            )
    else:
        # Depth map is the sole source: derive mask from depth >= depth_min
        gs.message(
            _(f"No input mask given — deriving ocean extent from depth map "
              f"(depth >= {depth_min})…")
        )
        mask = _ocean_mask_from_depth(depth_map, depth_min, pid)
        stats = gs.parse_command("r.univar", map=mask, flags="g", quiet=True)
        if int(float(stats.get("n", 0))) == 0:
            gs.fatal(
                _(
                    f"No pixels with depth >= {depth_min} found in <{depth_map}>. "
                    "Adjust depth_min or check the depth map."
                )
            )

    # ── depth index (0 = shore / surface, 1000 = abyss) ─────────────────────
    if depth_map:
        gs.message(_("Building depth index from bathymetric map…"))
        depth_idx = _depth_from_map(depth_map, mask, pid)
    elif flag_shore:
        gs.message(_("Estimating depth from distance to shore…"))
        depth_idx = _depth_from_shore(mask, pid)
    else:
        depth_idx = _flat_depth(mask, pid)

    # ── latitude warmth gradient ─────────────────────────────────────────────
    if flag_latgrd:
        gs.message(_("Applying north–south temperature gradient…"))
        depth_idx = _latitude_gradient(depth_idx, latitude, mask, pid)

    # ── wave texture ─────────────────────────────────────────────────────────
    wave_amp = 0
    if flag_waves:
        gs.message(
            _(f"Generating wave texture (pixel scale {pixel_m:,.0f} m)…")
        )
        wave_rast, wave_amp = _wave_pattern(mask, pixel_m, pid)
        blended = _register(f"tmp_iocean_blended_{pid}")
        gs.mapcalc(
            f"{blended} = if(isnull({mask}), null(), {depth_idx} + {wave_rast})",
            overwrite=True, quiet=True,
        )
        final_src = blended
    else:
        final_src = depth_idx

    # ── write output ─────────────────────────────────────────────────────────
    gs.message(_("Writing output raster…"))
    gs.mapcalc(
        f"{output} = if(isnull({mask}), null(), float({final_src}))",
        overwrite=True, quiet=True,
    )

    # ── colour rules ─────────────────────────────────────────────────────────
    gs.message(_("Applying ocean colour palette…"))
    palette = _PALETTES[style]
    rules   = _build_color_rules(palette, wave_amp)
    gs.write_command(
        "r.colors", map=output, rules="-", stdin=rules, quiet=True
    )

    # ── metadata ─────────────────────────────────────────────────────────────
    effects = []
    if depth_map:
        effects.append("depth-map")
    elif flag_shore:
        effects.append("shore-distance-depth")
    if flag_waves:
        effects.append("waves")
    if flag_latgrd:
        effects.append("lat-gradient")
    effects_str = ", ".join(effects) if effects else "none"

    source_name = input_map or depth_map
    gs.run_command(
        "r.support",
        map=output,
        title=f"Ocean visualisation of {source_name}",
        description=(
            f"style={style}, lat={latitude:.1f}\u00b0, "
            f"pixel={pixel_m:.0f}m, effects=[{effects_str}]"
        ),
        quiet=True,
    )
    gs.raster_history(output)

    gs.message(_(f"Done. Ocean visualisation written to <{output}>."))
    gs.message(_(f"Display with:  d.rast map={output}"))


if __name__ == "__main__":
    sys.exit(main())
