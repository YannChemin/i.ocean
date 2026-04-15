## DESCRIPTION

*i.ocean* renders ocean areas of a raster map with a realistic visual
appearance. It combines depth-based colour gradients, latitude-driven colour
temperature, optional wave textures, physically-based ocean turbulence with
coastal foam, and optional maritime decorations to produce a styled,
floating-point output raster with an automatically tuned ocean colour table.

The module is written in C with OpenMP parallelism and uses FFTW3 for
spectral ocean surface synthesis. At least one of **input** or **depth**
must be provided. Non-ocean pixels are set to null in the output so the
result composites cleanly over any land base map when displayed with
*[d.rast](https://grass.osgeo.org/grass-stable/manuals/d.rast.html)*.

Processing progress is reported as a percentage through both the terminal
and the GRASS GUI progress bar. All raster computations use all available
CPU cores automatically.

The **input** parameter specifies an ocean mask raster. Cells whose value
equals **ocean_value** (default: 1) are treated as ocean; all other cells
become null in the output.

The **depth** parameter specifies a bathymetric depth raster with positive
values representing depth in metres below the surface (e.g. GEBCO). When
**depth** is supplied it drives the colour gradient: shallow areas receive
lighter, warmer colours and deep areas receive darker, cooler colours. When
**input** is omitted, **depth** additionally defines the ocean extent by
thresholding: cells with value greater than or equal to **depth_min**
(default: 0) are treated as ocean. When the **depth** map is stored as an
integer raster (CELL type), it is automatically smoothed with a 3×3
neighbourhood average before normalisation to remove staircase artefacts
from the rendered gradient (see *Integer depth smoothing* below).

The **ocean_value** option sets the cell value in the **input** mask that
identifies ocean pixels. It is ignored when only **depth** is provided.

The **depth_min** option sets the minimum value in the **depth** map that is
considered ocean. Use 0 (default) for standard bathymetric maps where
positive values indicate depth below the surface. For DEM or elevation maps
where negative values represent water, set **depth_min** to a negative
value such as -0.1 to capture cells at or below sea level.

The **latitude** option sets the reference latitude in decimal degrees
(negative = south) used to select the colour palette. When omitted the
latitude is auto-detected from the centre of the current computational
region. For projected CRS the centre point is reprojected to geographic
coordinates using
*[m.proj](https://grass.osgeo.org/grass-stable/manuals/m.proj.html)*.

The **style** option overrides the auto-detected latitude zone and forces a
specific colour palette. Accepted values are *auto*, *tropical*,
*subtropical*, *temperate*, *subpolar*, and *polar*.

The **-w** flag adds a sinusoidal wave texture to the output. The wave
period and amplitude are calibrated to the pixel size of the current
computational region so that the result looks realistic at any scale, from
sub-100 m coastal surveys to global maps.

The **-s** flag estimates a depth proxy from the Euclidean distance of each
ocean pixel to the nearest land pixel, computed with
*[r.grow.distance](https://grass.osgeo.org/grass-stable/manuals/r.grow.distance.html)*.
It is used when **depth** is not provided and provides a convincing
continental shelf gradient without external bathymetric data.

The **-l** flag applies a subtle north–south warmth gradient across the
scene. Equatorward pixels are shifted towards shallower (lighter) colours
to mimic the visible sea-surface-temperature (SST) gradient typical of
scenes that span several degrees of latitude.

The **-f** flag generates a physically-based ocean surface using the
Phillips power spectrum (Tessendorf 2001) computed by FFTW3. The spectrum
encodes multi-scale wave energy as a function of wavenumber, wind speed,
and dominant length scale. A 2-D real inverse FFT produces the spatial
height field; Gaussian noise is drawn with Box-Muller per-thread for
thread-safe reproduction. Foam appears where wave crests exceed a
threshold (a squared-excess kernel inspired by the 0AD water_high.fs foam
model), and the foam intensity is amplified near the coast using a
two-pass Euclidean distance transform.

The **-d** flag scatters maritime decorations across the ocean surface.
Decorations are stamped as small pixel-art patterns that modify the local
depth-index value. Placement density is proportional to ocean area; a
fixed-seed PRNG ensures reproducible results. Supported types: whale,
dolphin, fish school, octopus, treasure chest, floating crate, sunken ship.

## NOTES

### Ocean mask sources

The ocean extent can come from two sources, used in priority order:

1. **input** mask — an existing raster in which one cell value marks water.
   The value is specified with **ocean_value**. This is the most precise
   option and is suitable when a land/sea classification or coastline raster
   is already available.
2. **depth** map alone — when no **input** mask is given, ocean pixels are
   identified as cells in the **depth** map whose value is at or above
   **depth_min**. This is convenient when working directly with a
   bathymetric dataset.

When both **input** and **depth** are provided, **input** defines the ocean
extent and **depth** drives only the colour gradient.

### Depth index

Internally all depth information is normalised to a *depth index* ranging
from 0 (shoreline / surface) to 1000 (abyssal plain). Three strategies are
available, applied in priority order:

1. **depth** map — depths are taken as absolute values and linearly rescaled
   to 0–1000 using the maximum depth found within the ocean mask.
2. **-s** flag — shore distance is computed with *r.grow.distance* and
   rescaled to 0–1000. Pixels close to the coast are treated as shallow,
   pixels far from land as deep.
3. Default — a flat index of 500 is applied, placing all ocean pixels at the
   mid-point of the colour palette.

### Integer depth smoothing

Many widely used bathymetric datasets (ETOPO, GEBCO at coarser resolutions,
national chart data) are distributed as integer rasters (GRASS CELL type).
When the depth index is normalised from integer values, the fixed step size
between adjacent integers produces visible staircase terracing in the
rendered gradient — each integer depth level maps to an identical colour
band.

When *i.ocean* detects a CELL-type **depth** map it applies a 3×3
neighbourhood average with
*[r.neighbors](https://grass.osgeo.org/grass-stable/manuals/r.neighbors.html)*
before normalisation. The *average* method in *r.neighbors* always produces
a DCELL (double-precision floating-point) output regardless of input type
(see *r.neighbors* source, `main.c` lines 64 and 117–135): each output cell
is the mean of its 3×3 integer neighbourhood, yielding sub-pixel fractional
transitions across depth step boundaries. The resulting smooth floating-point
raster is then normalised to the 0–1000 depth index in the usual way.

FCELL and DCELL depth maps are passed through unchanged; no smoothing is
applied because the values are already continuous.

### Phillips-spectrum turbulence and foam (-f)

The ocean surface height field is generated in the spectral domain using
the Phillips power spectrum:

```
P(k) = A · exp(−1/(k·L)²) / k⁴ · |k̂ · ŵ|²
```

where *k* is the wavenumber magnitude, *L = V²/g* is the dominant length
scale derived from wind speed *V*, and *ŵ* is the unit wind direction
vector. Complex Gaussian noise (Box-Muller, thread-local `rand_r`) is
multiplied by √P(k)/2 to populate the frequency-domain array, which is
then transformed to spatial heights by a 2-D real inverse FFT
(`fftw_plan_dft_c2r_2d`). The height field is normalised so that 3σ spans
the target turbulence amplitude, which scales with pixel size.

Foam is detected where wave crests exceed a threshold fraction of the total
amplitude; a squared-excess kernel amplifies near-crest pixels (inspired by
the 0AD `water_high.fs` lines 343–356 foam accumulation formula). Coastal
foam is boosted using an exponential decay factor applied to a two-pass
Euclidean distance transform from land pixels.

### Maritime decorations (-d)

Seven decoration types are scattered proportionally to ocean area:

| Type | Pixels | Depth range | Density |
|---|---|---|---|
| Whale | 13 | 400–900 | 0.0008/km² |
| Dolphin | 5 | 100–700 | 0.002/km² |
| Fish school | 8 | 0–800 | 0.008/km² |
| Octopus | 9 | 50–600 | 0.003/km² |
| Treasure chest | 4 | 0–500 | 0.0001/km² |
| Floating crate | 6 | 0–350 | 0.001/km² |
| Sunken ship | 20 | 200–900 | 0.0003/km² |

Each decoration modifies the depth-index of its pixels. Floating objects
(dolphins, crates, treasure) use large negative deltas so they appear as
bright surface features. Sunken ships use positive deltas on the hull and
negative on the deck/mast for a shadow-and-highlight silhouette. Placement
uses a fixed seed (reproducible). Decorations are most visible at fine
resolutions (pixel size < 500 m).

### Multi-core parallelism

The module is written in C with OpenMP. All inner row loops — mask
creation, depth normalisation, 3×3 integer smoothing, wave generation,
Phillips-spectrum frequency fill, turbulence/foam, and shore distance —
use `#pragma omp parallel for` with `omp_get_max_threads()` threads.
To limit CPU usage set `OMP_NUM_THREADS` before running the module.

The two-pass Euclidean distance transform used by the **-s** and **-f**
flags is implemented directly in C (Meijster 2000 algorithm); no external
module subprocess is needed.

### Progress reporting

Processing is split into discrete named steps (region analysis, ocean mask
extraction, depth index computation, optional gradient and wave steps, output
write, colour rules, metadata). At the start of each step a percentage value
is emitted via `gs.percent()`, updating both the terminal progress indicator
and the GRASS GUI progress bar. The total step count is computed at startup
from the combination of flags supplied, so the percentage always advances
uniformly to 100 %.

### Scale-aware wave textures

When the **-w** flag is active the wave period and amplitude are chosen
according to the mean pixel size of the current computational region:

| Pixel size | Represented feature | Period | Amplitude |
|---|---|---|---|
| < 100 m | Fine ripples | ~30 m | 12 |
| 100 m – 1 km | Ocean swells | ~600 m | 25 |
| 1 – 10 km | Large swells / mesoscale eddies | ~6 km | 35 |
| 10 – 100 km | Broad surface patterns | ~60 km | 30 |
| > 100 km | Large gyre-like undulations | 10 px | 20 |

The texture is a superposition of a primary cross-wave (sin × cos) and a
40 % diagonal harmonic to break the visual regularity of the grid.

### Colour palettes

Five palettes are available, each with seven colour stops spanning the full
depth index range from shoreline turquoise / blue to abyssal dark blue:

| Style | Latitude zone | Character |
|---|---|---|
| *tropical* | < 23.5° | Bright turquoise / cyan |
| *subtropical* | 23.5 – 35° | Warm blue-green |
| *temperate* | 35 – 60° | Clear cool blue |
| *subpolar* | 60 – 75° | Steel blue |
| *polar* | > 75° | Icy pale blue |

### Output colour table

The colour table is applied directly to the output raster with *r.colors*
and is stored persistently in the map's metadata. It survives region changes
and is reproduced whenever the map is re-displayed. Non-ocean (null) cells
are mapped to white (255:255:255) so that transparency is available in
display modules that support it.

## EXAMPLES

### Basic ocean rendering from a land/sea mask

Create a simple ocean mask from an elevation raster (negative values = ocean)
and render it with the default settings. The style is auto-detected from the
map centre latitude.

```sh
r.mapcalc "sea_mask = if(elevation < 0, 1, null())"
i.ocean input=sea_mask output=ocean_view
d.rast map=ocean_view
```

### Using a bathymetric map as the sole input

When a GEBCO-style depth raster is available it can be passed directly as
**depth** without a separate mask. Cells with depth ≥ 0 become ocean. If the
raster is integer type (CELL), the stepped gradient is smoothed automatically.

```sh
i.ocean depth=gebco_bathy output=ocean_view
```

### Using a DEM as the depth source (elevation map, negative = ocean)

For a DEM where ocean pixels have negative elevation values, set
**depth_min** to a small negative threshold so only cells below sea level
are treated as ocean:

```sh
i.ocean depth=elevation output=ocean_view depth_min=-0.1
```

### Combining an explicit mask with a bathymetric depth map

The mask controls the ocean extent; the depth map drives the gradient.

```sh
i.ocean input=sea_mask depth=gebco_bathy output=ocean_deep
d.rast map=ocean_deep
```

### Adding wave texture and shore-distance depth

Estimate depth from distance to the coast and overlay scale-aware wave
patterns.

```sh
i.ocean input=sea_mask output=ocean_waves -w -s
d.rast map=ocean_waves
```

### Phillips-spectrum turbulence with coastal foam

```sh
i.ocean input=sea_mask depth=gebco_bathy output=ocean_turbulent -f
d.rast map=ocean_turbulent
```

### Full effect stack

All effects combined: waves, turbulence/foam, latitude gradient, and
maritime decorations.

```sh
i.ocean input=sea_mask depth=gebco_bathy output=ocean_full \
        style=tropical -w -f -l -d
d.rast map=ocean_full
```

### Forcing a tropical palette with all effects enabled

```sh
i.ocean input=sea_mask depth=gebco_bathy output=ocean_tropical \
        style=tropical -w -l
d.rast map=ocean_tropical
```

### Polar scene with icy palette

```sh
i.ocean input=arctic_mask output=arctic_ocean style=polar -w -s
d.rast map=arctic_ocean
```

### Overlaying ocean on a land hillshade

Render the ocean and combine it with a hillshade for a complete cartographic
basemap.

```sh
r.relief input=elevation output=hillshade
i.ocean input=sea_mask depth=gebco_bathy output=ocean_deep -w -s -l
d.shade shade=hillshade color=ocean_deep
```

## SEE ALSO

*[r.neighbors](https://grass.osgeo.org/grass-stable/manuals/r.neighbors.html)*,
*[r.colors](https://grass.osgeo.org/grass-stable/manuals/r.colors.html)*,
*[r.grow.distance](https://grass.osgeo.org/grass-stable/manuals/r.grow.distance.html)*,
*[r.relief](https://grass.osgeo.org/grass-stable/manuals/r.relief.html)*,
*[r.lake](https://grass.osgeo.org/grass-stable/manuals/r.lake.html)*,
*[d.rast](https://grass.osgeo.org/grass-stable/manuals/d.rast.html)*,
*[d.shade](https://grass.osgeo.org/grass-stable/manuals/d.shade.html)*,
*[m.proj](https://grass.osgeo.org/grass-stable/manuals/m.proj.html)*

## AUTHORS

i.ocean contributors, 2026.
