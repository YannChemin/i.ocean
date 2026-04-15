# i.ocean

A [GRASS GIS](https://grass.osgeo.org/) addon that renders the ocean areas
of a raster map with a realistic visual appearance.

Written in **C with OpenMP** and **FFTW3** spectral synthesis.  Given a
raster in which one cell value marks water (default `1`), the module
produces a styled, floating-point output raster with an automatically-tuned
ocean colour table.  Non-ocean pixels become null so the result composites
cleanly over any land base-map.

## Features

- **Five latitude-driven colour palettes** — from bright tropical turquoise to
  icy polar blue, auto-detected from the map centre or set explicitly.
- **Depth-based colour gradient** — accepts a bathymetric raster; falls back to
  shore-distance estimation (`-s`) or a flat mid-depth default.
- **Automatic integer depth smoothing** — when the depth map is integer type
  (CELL), a 3×3 neighbourhood average converts it to a continuous floating-point
  (DCELL) raster before normalisation, removing staircase artefacts.
- **Scale-aware wave textures** (`-w`) — sinusoidal wave patterns calibrated to
  the current pixel size across five scale tiers (sub-100 m ripples → 100 km
  gyre-scale undulations).
- **Phillips-spectrum turbulence + coastal foam** (`-f`) — FFTW3 computes a
  physically-based ocean surface via the Phillips power spectrum; foam appears
  at wave crests using a squared-excess threshold, amplified near the coast
  with an exponential decay factor derived from a two-pass Euclidean distance
  transform (Meijster 2000).
- **Maritime decorations** (`-d`) — whale, dolphin, fish school, octopus,
  treasure chest, floating crate, and sunken ship pixel-art patterns scattered
  proportionally to ocean area with reproducible seeded placement.
- **North–south warmth gradient** (`-l`) — mimics the visible SST gradient
  within the scene by subtly shifting colour temperature equatorward to poleward.
- **Full CRS support** — geographic (degrees) and projected (metres, feet, km);
  pixel sizes are converted to metres internally.
- **Multi-core parallel processing** — OpenMP `#pragma omp parallel for` on all
  inner loops; uses `omp_get_max_threads()` automatically.
- **Progress reporting** — `G_percent()` at each named step, visible in the
  terminal and the GRASS GUI progress bar.

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
#    Integer (CELL) depth maps are smoothed automatically
i.ocean input=sea_mask depth=gebco_bathy output=ocean_deep -w

# 4. Shore-estimated depth + latitude warmth gradient
i.ocean input=sea_mask output=ocean_rich -s -l

# 5. Phillips-spectrum turbulence with coastal foam
i.ocean input=sea_mask depth=gebco_bathy output=ocean_turbulent -f

# 6. Maritime decorations (whales, ships, treasure...)
i.ocean input=sea_mask depth=gebco_bathy output=ocean_deco -f -d

# 7. Full effect stack
i.ocean input=sea_mask depth=gebco_bathy output=ocean_full \
        style=tropical -w -f -l -d

# 8. Display the result
d.rast map=ocean_full
```

## Parameters

| Parameter | Type | Default | Description |
|---|---|---|---|
| `input` | raster | — | Ocean mask raster (cells equal to `ocean_value` are ocean) |
| `output` | raster | — | Output ocean visualisation raster |
| `depth` | raster | — | Bathymetric depth map (positive = metres below surface); integer maps are smoothed automatically |
| `ocean_value` | float | `1` | Cell value that marks ocean pixels in `input` |
| `depth_min` | float | `0` | Minimum depth value treated as ocean when deriving the mask from `depth`; use a negative value for DEM/elevation maps |
| `latitude` | float | auto | Reference latitude (decimal degrees, negative = south); auto-detected from map centre if omitted |
| `style` | string | `auto` | Colour style: `auto`, `tropical`, `subtropical`, `temperate`, `subpolar`, `polar` |

## Flags

| Flag | Description |
|---|---|
| `-w` | Add scale-aware wave / ripple texture |
| `-s` | Estimate depth from distance to shore (used when no `depth` map is given) |
| `-l` | Apply subtle north–south warmth gradient across the scene |

## Requirements

- GRASS GIS ≥ 8.0
- C99 compiler with OpenMP support (`gcc -fopenmp`)
- `libfftw3-dev` — Fast Fourier Transform (Phillips spectrum)
- `libomp-dev` / `libgomp` — OpenMP runtime (usually bundled with gcc)
- Standard GRASS modules called at runtime: `r.colors`, `r.support`,
  `g.proj`, `m.proj`, `g.remove`

## Installation

See [INSTALL.md](INSTALL.md).
