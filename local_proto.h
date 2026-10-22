/****************************************************************************
 *
 * MODULE:       i.ocean  –  local prototypes and shared types
 * AUTHORS:      i.ocean contributors
 * PURPOSE:      Shared data structures and function declarations used across
 *               the i.ocean C source files.
 *
 ****************************************************************************/

#ifndef LOCAL_PROTO_H
#define LOCAL_PROTO_H

#include <grass/gis.h>
#include <grass/raster.h>

/* ── palette entry (depth_index 0-1000, RGB) ──────────────────────────── */
typedef struct {
    int   idx;
    int   r, g, b;
} PaletteStop;

/* ── compact in-memory raster ─────────────────────────────────────────── */
typedef struct {
    DCELL  *data;   /* flat row-major array, nrows*ncols cells             */
    int     nrows;
    int     ncols;
} Rbuf;

#define RBUF_CELL(rb, r, c)  ((rb)->data[(r) * (rb)->ncols + (c)])
#define RBUF_IS_NULL(rb, r, c) \
    Rast_is_d_null_value(&RBUF_CELL(rb, r, c))

/* ── ocean processing parameters ─────────────────────────────────────── */
typedef struct {
    /* region */
    double ns_res;       /* N-S pixel size in native CRS units            */
    double ew_res;       /* E-W pixel size in native CRS units            */
    double pixel_m;      /* mean pixel size in metres                     */
    double latitude;     /* reference latitude (decimal degrees)          */
    /* flags */
    int flag_waves;
    int flag_shore;
    int flag_latgrd;
    int flag_foam;
    int flag_deco;
    /* options */
    const char *style;
    int nprocs;
} OceanParams;

/* ── decoration placement record ─────────────────────────────────────── */
/*
 * Populated by place_decorations(); consumed by create_deco_vector().
 * Stores the map-coordinate centre and heading of each placed item.
 */
typedef struct {
    int    type_idx;   /* index into DECORATIONS[] / VSHAPES[]    */
    double east;       /* map-coordinate east of centre            */
    double north;      /* map-coordinate north of centre           */
    double angle_rad;  /* heading: 0 = facing north, π/2 = east   */
} PlacedDeco;

/* ── RAM check ────────────────────────────────────────────────────────── */
typedef enum {
    RAM_OK    = 0,   /* <= 60 % available                         */
    RAM_WARN  = 1,   /* 60–80 %                                   */
    RAM_TIGHT = 2,   /* 80–95 % — caller should disable -f        */
    RAM_FATAL = 3    /* >= 95 % — G_fatal_error already called    */
} RamStatus;

/* Number of decoration types (must stay in sync with DECORATIONS[]) */
#define N_DECO_TYPES 27

/* Real-world sizes for each type (meters); defined in decorations.c */
extern const double deco_real_length[N_DECO_TYPES];
extern const double deco_real_width [N_DECO_TYPES];
extern const char  *deco_type_name  [N_DECO_TYPES];

/* ── temp-map registry ────────────────────────────────────────────────── */
void   tmp_register(const char *name);
void   tmp_cleanup(void);

/* ── raster I/O ───────────────────────────────────────────────────────── */
Rbuf  *rbuf_alloc(int nrows, int ncols);
void   rbuf_free(Rbuf *rb);
Rbuf  *rbuf_read(const char *name);
void   rbuf_write(const Rbuf *rb, const char *name, RASTER_MAP_TYPE type);

/* ── ram_check.c ──────────────────────────────────────────────────────── */
RamStatus ram_check(long nrows, long ncols,
                    int has_w, int has_f, int has_s, int has_d,
                    double *est_mb_out, double *avail_mb_out);

/* ── ocean_mask.c ─────────────────────────────────────────────────────── */
Rbuf  *make_mask_from_value(const Rbuf *input, double ocean_val,
                            int nrows, int ncols);
Rbuf  *make_mask_from_depth(const Rbuf *depth, double depth_min,
                            int nrows, int ncols);

/* ── depth_index.c ────────────────────────────────────────────────────── */
void   smooth_3x3(Rbuf *buf, const Rbuf *mask, int nprocs);
Rbuf  *depth_from_map(const char *depth_name, const Rbuf *mask,
                      int nprocs, const char *pid_str);
Rbuf  *depth_from_shore(const Rbuf *mask, double pixel_m,
                        int nprocs, const char *pid_str);
Rbuf  *depth_flat(const Rbuf *mask, double value, int nprocs);

/* ── waves.c ──────────────────────────────────────────────────────────── */
Rbuf  *wave_pattern(const Rbuf *mask, double pixel_m,
                    double *amp_out, int nprocs);

/* ── turbulence.c ─────────────────────────────────────────────────────── */

/*
 * foam_turbulence() — multi-octave Phillips-spectrum ocean turbulence
 * with coastal foam amplification.
 *
 * Technique reference (0AD water_high.fs, lines 343-356):
 *   • Dual-layer normal maps blended with a modulated time factor give
 *     multi-scale wave banding (primary swell × cross-wave envelope).
 *   • Foam appears where wave crests exceed a threshold; intensity is
 *     scaled by pow(amplitude, 2.6 - waviness/5.5).
 *   • Shore proximity is tracked via a distance field and used to boost
 *     foam intensity near land (coastal breaker zone).
 *
 * Our raster implementation uses FFTW3 (libfftw3-dev) to generate the
 * height field in frequency space via the Phillips spectrum:
 *
 *   P(k) = A · exp(−1/(k·L)²) / k⁴ · |k̂ · ŵ|²
 *
 * where  k = wavenumber,  L = V²/g  (dominant scale from wind speed V),
 * ŵ = unit wind direction vector,  A = amplitude constant.
 *
 * A 2-D real inverse FFT (fftw_plan_dft_c2r_2d) produces the height field.
 * Thread-safe Box-Muller noise is generated per-thread with rand_r().
 * OpenMP parallelises the frequency-domain fill and the spatial-domain
 * foam/coast passes.
 */
Rbuf  *foam_turbulence(const Rbuf *depth_idx, const Rbuf *mask,
                       double pixel_m, double ns_res, double ew_res,
                       double *lo_clamp_out, int nprocs,
                       const char *pid_str);

/* ── decorations.c ────────────────────────────────────────────────────── */

/*
 * place_decorations() — scatter sea-life and maritime objects onto the
 * depth-index buffer.
 *
 * Real-world sizes (deco_real_length / deco_real_width) are defined
 * in decorations.c.  Small objects (treasure chest, floating crate)
 * are bumped up to the dolphin minimum (3 m) so all items share a
 * homogeneous minimum ground size.
 *
 * Returns a G_malloc-ed array of PlacedDeco structs (*n_placed_out items).
 * The caller is responsible for calling G_free() on the returned pointer.
 * Returns NULL and sets *n_placed_out = 0 if no items were placed.
 *
 * deco_scale multiplies all vector silhouette dimensions but does NOT
 * change the raster pixel-art stamp (which is resolution-dependent).
 */
PlacedDeco *place_decorations(Rbuf *work, const Rbuf *mask,
                               const struct Cell_head *region,
                               double pixel_m, double deco_scale,
                               double lo_clamp,
                               int *n_placed_out);

/* ── vector_deco.c ────────────────────────────────────────────────────── */

/*
 * create_deco_vector() — create a GRASS vector map of maritime decoration
 * silhouettes at real-world scale.
 *
 * Each PlacedDeco entry is drawn as a set of GV_LINE primitives using the
 * 2-D silhouette shape for its type.  The shape coordinates (in metres)
 * are:
 *   • scaled by deco_scale
 *   • rotated by angle_rad (clockwise from north)
 *   • converted to map units (geographic or projected CRS)
 *   • offset from the placement centre (east, north)
 *
 * The vector map carries one attribute table with columns:
 *   cat INTEGER, type VARCHAR(32), len_m REAL, wid_m REAL, scale REAL
 *
 * Shape data is derived from 0AD structural mesh bounding boxes:
 *   - Trireme hull (hele_ship_bireme.dae): 8.8 × 49.2 m → used for sunken ship
 *   - fish_generic.dae: 8.9 × 6.1 m → scaled down for fish school
 *   - Whale/dolphin: real cetacean biology (blue whale ~28 m, dolphin ~3 m)
 *
 * proj   : GRASS projection string (e.g. "ll" for geographic, "utm" for UTM)
 * latitude: used for EW degree ↔ metre conversion in geographic CRS
 */
void create_deco_vector(const char *vmap_name,
                        const PlacedDeco *places, int n_places,
                        const struct Cell_head *region,
                        const char *proj, const char *units,
                        double latitude, double deco_scale);

/* ── color_rules.c helpers in main.c ─────────────────────────────────── */
const char *auto_style(double lat);

#endif /* LOCAL_PROTO_H */
