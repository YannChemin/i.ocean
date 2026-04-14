# i.ocean

A [GRASS GIS](https://grass.osgeo.org/) addon that renders the ocean areas
of a raster map with a realistic visual appearance.

Given a raster in which one cell value marks water (default `1`), the module
produces a styled, floating-point output raster with an automatically-tuned
ocean colour table.  Non-ocean pixels become null so the result composites
cleanly over any land base-map.

## Features

- **Five latitude-driven colour palettes** — from bright tropical turquoise to
  icy polar blue, auto-detected from the map centre or set explicitly.
- **Depth-based colour gradient** — accepts a bathymetric raster; falls back to
  shore-distance estimation (`-s`) or a flat mid-depth default.
- **Scale-aware wave textures** (`-w`) — sinusoidal wave patterns whose period
  and amplitude are calibrated to the current pixel size across six scale tiers
  (sub-100 m ripples → 100 km gyre-scale undulations).
- **North–south warmth gradient** (`-l`) — mimics the visible SST gradient
  within the scene by subtly shifting colour temperature from equatorward to
  poleward.
- **Full CRS support** — geographic (degrees) and projected (metres, feet, km)
  coordinate systems; pixel sizes are converted to metres internally.

## Colour palettes

| Style | Latitude | Character |
|---|---|---|
| `tropical` | < 23.5° | Bright turquoise / cyan |
| `subtropical` | 23.5–35° | Warm blue-green |
| `temperate` | 35–60° | Clear cool blue |
| `subpolar` | 60–75° | Steel blue |
| `polar` | > 75° | Icy pale blue |

## Quick start

```sh
# 1. Create a land/sea mask (ocean = 1)
r.mapcalc "sea_mask = if(elevation < 0, 1, null())"

# 2. Basic render — auto palette, flat depth
i.ocean input=sea_mask output=ocean_view

# 3. With bathymetry, waves, and shore-depth fallback
i.ocean input=sea_mask depth=gebco_bathy output=ocean_deep -w

# 4. Shore-estimated depth + latitude warmth gradient
i.ocean input=sea_mask output=ocean_rich -s -l

# 5. Force a tropical palette
i.ocean input=sea_mask output=ocean_tropical style=tropical -w -s

# 6. Display the result
d.rast map=ocean_view
```

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input` | raster | — | Input raster map |
| `output` | raster | — | Output ocean visualisation raster |
| `depth` | raster | — | Bathymetric depth map (positive = metres below surface) |
| `ocean_value` | float | `1` | Cell value that marks ocean pixels |
| `latitude` | float | auto | Reference latitude (decimal degrees); auto-detected if omitted |
| `style` | string | `auto` | Colour style: `auto`, `tropical`, `subtropical`, `temperate`, `subpolar`, `polar` |

## Flags

| Flag | Description |
|---|---|
| `-w` | Add scale-aware wave / ripple texture |
| `-s` | Estimate depth from distance to shore (requires no depth map) |
| `-l` | Apply subtle north–south warmth gradient across the scene |

## Requirements

- GRASS GIS ≥ 8.0
- Python ≥ 3.8
- Standard GRASS modules used: `r.mapcalc`, `r.colors`, `r.grow.distance`,
  `r.univar`, `r.support`, `g.proj`, `m.proj`, `g.remove`

## Installation

See [INSTALL.md](INSTALL.md).

