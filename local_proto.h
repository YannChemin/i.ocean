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

/* ── temp-map registry ────────────────────────────────────────────────── */
void   tmp_register(const char *name);
void   tmp_cleanup(void);

/* ── raster I/O ───────────────────────────────────────────────────────── */
Rbuf  *rbuf_alloc(int nrows, int ncols);
void   rbuf_free(Rbuf *rb);
Rbuf  *rbuf_read(const char *name);
void   rbuf_write(const Rbuf *rb, const char *name, RASTER_MAP_TYPE type);

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
 * Each decoration type is defined by a small pixel-pattern (offset + delta),
 * a depth-index range in which it can appear, and a density (items per km²).
 * Positions are chosen from a seeded PRNG so results are reproducible.
 *
 * Types:
 *   whale        — 13-pixel elongated body, bright surface sheen
 *   dolphin      — 5-pixel arc, very bright
 *   fish_school  — 8 scattered pixels, mottled bright
 *   octopus      — 9-pixel radial body+tentacles, alternating tones
 *   treasure     — 4-pixel bright "gold" cluster (depth_idx → 0)
 *   crate        — 6-pixel floating box (bright)
 *   ship         — 20-pixel sunken hull (dark keel, bright deck)
 */
void  place_decorations(Rbuf *work, const Rbuf *mask,
                        double pixel_m, double lo_clamp);

/* ── color_rules.c helpers in main.c ─────────────────────────────────── */
const char *auto_style(double lat);

#endif /* LOCAL_PROTO_H */
