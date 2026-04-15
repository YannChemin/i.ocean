/****************************************************************************
 *
 * MODULE:       i.ocean  —  turbulence.c
 * PURPOSE:      Multi-scale ocean surface turbulence + coastal foam using
 *               the Phillips spectrum (FFTW3).
 *
 * ALGORITHM OVERVIEW
 * ──────────────────
 * Inspired by 0AD water rendering (water_high.fs, WaterManager.cpp):
 *
 *   0AD reference points translated to raster equivalents
 *   ──────────────────────────────────────────────────────
 *   • water_high.fs lines 149-176: dual-layer normal maps blended with
 *     a time-cycled weight give multi-scale wave banding; we replace the
 *     pre-baked texture stack with a procedural IFFT over the Phillips
 *     power spectrum — same multi-frequency character, no texture storage.
 *
 *   • water_high.fs lines 343-356: foam is the squared excess above a
 *     threshold:  foam = (max(0, amplitude/A_total − thresh))^2.6  *  scale
 *     We use a slightly different exponent (2.0) for raster aesthetics.
 *
 *   • WaterManager.cpp lines 454-479: shore-proximity distance field used
 *     to amplify foam in the coastal breaker zone.  We replicate this with
 *     a two-pass Euclidean distance transform (same BFS kernel as
 *     depth_index.c:depth_from_shore) followed by an exponential decay
 *     factor:  coast_boost = 1 + 4 · exp(−dist / shore_decay_m).
 *
 * PHYSICS MODEL
 * ─────────────
 * Phillips spectrum (Tessendorf 2001 "Simulating Ocean Water"):
 *
 *   P(k) = A · exp(−1 / (k·L)²) / k⁴ · |k̂ · ŵ|²
 *
 *   k    = wavenumber magnitude (rad/m)
 *   L    = V² / g               (dominant scale from wind speed V m/s)
 *   ŵ    = unit wind direction   (default along-column, i.e. (1, 0))
 *   A    = amplitude constant    (tuned to produce ±turb_amp depth-index units)
 *
 * The height field at time t=0 is:
 *
 *   H(x) = Re{ IFFT( h₀(k) ) }
 *   h₀(k) = (1/√2) · (ξᵣ + i·ξᵢ) · √P(k)
 *
 * where ξᵣ, ξᵢ ~ N(0,1) (Box-Muller, thread-safe rand_r with per-thread seed).
 *
 * FFTW3 USAGE
 * ───────────
 * Plan: fftw_plan_dft_c2r_2d(nrows, ncols, freq, spatial, FFTW_ESTIMATE)
 *   • Input:  fftw_complex array of size nrows × (ncols/2+1)
 *     (only the non-redundant Hermitian half for a real output transform)
 *   • Output: double array of size nrows × ncols
 *   • Normalisation: divide by sqrt(nrows * ncols) after IFFT to get
 *     unit-variance output; then scale to ±turb_amp.
 *
 * OPENMP
 * ──────
 * • Frequency-domain fill: parallel over rows (independent per row).
 * • Height-to-depth-index blend: parallel over all pixels.
 * • Shore-distance DT: parallel column pass and Meijster row pass.
 * • Foam accumulation: parallel row loop.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <fftw3.h>
#include <omp.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>

#include "local_proto.h"

/* ── Phillips spectrum ───────────────────────────────────────────────── */

/*
 * phillips() — single-frequency Phillips power spectrum value.
 *
 * kx, ky   : wavenumber components (rad / pixel_unit)
 * L        : dominant length scale V^2/g
 * wx, wy   : unit wind direction vector
 * A        : amplitude constant
 *
 * Returns P(k) ≥ 0.  Returns 0 for k = 0 (DC component → zero mean).
 */
static double phillips(double kx, double ky,
                       double L, double wx, double wy, double A)
{
    double k2 = kx * kx + ky * ky;
    if (k2 < 1e-12) return 0.0;

    double k  = sqrt(k2);
    double kL = k2 * L * L;
    /* direction factor: squared cosine between k and wind */
    double kdotw = (kx * wx + ky * wy) / k;
    /* suppress waves running perpendicular to the wind */
    if (kdotw < 0.0) kdotw = 0.0;

    return A * exp(-1.0 / kL) / (k2 * k2) * kdotw * kdotw;
}

/* ── Box-Muller Gaussian sample (thread-safe via rand_r) ─────────────── */

static double gaussian_bm(unsigned int *seed)
{
    double u1, u2;
    do {
        u1 = (rand_r(seed) + 0.5) / ((double)RAND_MAX + 1.0);
    } while (u1 <= 0.0);
    u2 = (rand_r(seed) + 0.5) / ((double)RAND_MAX + 1.0);
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* ── two-pass shore distance for foam boost ─────────────────────────── */

/*
 * shore_dist_buf() — compute Euclidean distance (in metres) from each
 * ocean pixel to the nearest land (null-mask) pixel.
 *
 * Uses the same two-pass Meijster DT as depth_index.c:depth_from_shore,
 * but returns a raw distance buffer in pixel units (caller multiplies by
 * pixel_m if needed).
 */
static float *shore_dist_pixels(const Rbuf *mask, int nprocs)
{
    int   nrows = mask->nrows, ncols = mask->ncols;
    long  npx   = (long)nrows * ncols;
    float *g    = G_malloc(sizeof(float) * npx);

    omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
    for (int c = 0; c < ncols; c++) {
        float last = (float)(nrows * 2);
        for (int r = 0; r < nrows; r++) {
            last = Rast_is_d_null_value(&RBUF_CELL(mask, r, c))
                   ? 0.0f : last + 1.0f;
            g[r * ncols + c] = last;
        }
        last = (float)(nrows * 2);
        for (int r = nrows - 1; r >= 0; r--) {
            last = Rast_is_d_null_value(&RBUF_CELL(mask, r, c))
                   ? 0.0f : last + 1.0f;
            if (last < g[r * ncols + c])
                g[r * ncols + c] = last;
        }
    }

    float *dist = G_malloc(sizeof(float) * npx);

#pragma omp parallel for schedule(static)
    for (int r = 0; r < nrows; r++) {
        int   *v = G_malloc(sizeof(int)   * ncols);
        float *z = G_malloc(sizeof(float) * (ncols + 1));
        int q = 0;
        v[0] = 0; z[0] = -1e30f; z[1] = 1e30f;
        for (int c = 1; c < ncols; c++) {
            float gc = g[r * ncols + c];
            float s;
            do {
                float gv = g[r * ncols + v[q]];
                s = ((gc * gc + (float)(c * c))
                     - (gv * gv + (float)(v[q] * v[q])))
                    / (2.0f * (c - v[q]));
            } while (q > 0 && s <= z[q] && --q >= 0);
            v[++q] = c;
            z[q]   = s;
            z[q+1] = 1e30f;
        }
        q = 0;
        for (int c = 0; c < ncols; c++) {
            while (z[q+1] < (float)c) q++;
            float dv = g[r * ncols + v[q]];
            float dc = (float)(c - v[q]);
            dist[r * ncols + c] = sqrtf(dv * dv + dc * dc);
        }
        G_free(v);
        G_free(z);
    }
    G_free(g);
    return dist;   /* in pixel units */
}

/* ── scale-dependent turbulence parameters ───────────────────────────── */

typedef struct {
    double turb_amp;       /* target ±depth-index amplitude of height field */
    double foam_thresh;    /* fraction of turb_amp above which foam appears  */
    double foam_scale;     /* foam brightness multiplier                     */
    double shore_decay_px; /* shore foam boost e-folding scale in pixels     */
    double wind_speed;     /* V (m/s) for Phillips L = V^2/g                */
} TurbParams;

static TurbParams turb_params(double pixel_m)
{
    TurbParams p;
    if (pixel_m < 100.0) {
        p.turb_amp       = 18.0;
        p.foam_thresh    = 0.50;
        p.foam_scale     = 55.0;
        p.shore_decay_px = fmax(5.0,  8.0 * 100.0  / pixel_m);
        p.wind_speed     = 8.0;
    } else if (pixel_m < 1000.0) {
        p.turb_amp       = 25.0;
        p.foam_thresh    = 0.45;
        p.foam_scale     = 70.0;
        p.shore_decay_px = fmax(5.0,  6.0 * 600.0  / pixel_m);
        p.wind_speed     = 10.0;
    } else if (pixel_m < 10000.0) {
        p.turb_amp       = 30.0;
        p.foam_thresh    = 0.40;
        p.foam_scale     = 80.0;
        p.shore_decay_px = fmax(5.0,  4.0 * 4000.0 / pixel_m);
        p.wind_speed     = 14.0;
    } else {
        p.turb_amp       = 22.0;
        p.foam_thresh    = 0.45;
        p.foam_scale     = 65.0;
        p.shore_decay_px = fmax(5.0,  2.0 * 40000.0/ pixel_m);
        p.wind_speed     = 18.0;
    }
    return p;
}

/* ── foam_turbulence ─────────────────────────────────────────────────── */

Rbuf *foam_turbulence(const Rbuf *depth_idx, const Rbuf *mask,
                       double pixel_m, double ns_res, double ew_res,
                       double *lo_clamp_out, int nprocs,
                       const char *pid_str)
{
    (void)pid_str;

    int  nrows = mask->nrows, ncols = mask->ncols;
    long npx   = (long)nrows * ncols;

    TurbParams p = turb_params(pixel_m);
    G_verbose_message(
        _("Phillips turbulence: turb_amp=%.0f foam_thresh=%.2f "
          "foam_scale=%.0f wind=%.0f m/s"),
        p.turb_amp, p.foam_thresh, p.foam_scale, p.wind_speed);

    /* ── allocate FFTW arrays ────────────────────────────────────────── */
    int     ncols_h = ncols / 2 + 1;   /* half-spectrum column count      */
    fftw_complex *freq    = fftw_alloc_complex((size_t)nrows * ncols_h);
    double       *spatial = fftw_alloc_real  ((size_t)nrows * ncols);

    /* ── FFTW plan (must be created before filling the array) ────────── */
    fftw_plan plan = fftw_plan_dft_c2r_2d(nrows, ncols,
                                           freq, spatial,
                                           FFTW_ESTIMATE);

    /* ── Phillips spectrum parameters ────────────────────────────────── */
    double g_grav  = 9.81;
    double V       = p.wind_speed;
    double L       = V * V / g_grav;           /* dominant length scale (pixels) */
    /* Convert L from metres to pixels */
    L /= pixel_m;

    /* wind direction: primary swell along-column (1,0) */
    double wx = 1.0, wy = 0.0;

    /*
     * A tuning: we want the IFFT output standard deviation ≈ 1.0 so we
     * can scale it later.  Set A = 1 here; we renormalise after IFFT.
     */
    double A = 1.0;

    /* ── fill frequency domain (OpenMP, thread-local RNG) ────────────── */
    omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
    for (int i = 0; i < nrows; i++) {
        unsigned int seed = (unsigned int)(42 + i + nrows * omp_get_thread_num());
        /* row frequency: handle negative frequencies */
        int ki_signed = (i <= nrows / 2) ? i : i - nrows;
        double ky = 2.0 * M_PI * ki_signed / (nrows * ns_res / pixel_m);

        for (int j = 0; j < ncols_h; j++) {
            double kx = 2.0 * M_PI * j / (ncols * ew_res / pixel_m);

            double ph  = phillips(kx, ky, L, wx, wy, A);
            double amp = sqrt(ph * 0.5);

            double xi_r = gaussian_bm(&seed);
            double xi_i = gaussian_bm(&seed);

            long idx = (long)i * ncols_h + j;
            freq[idx][0] = amp * xi_r;
            freq[idx][1] = amp * xi_i;
        }
    }

    /* ── inverse FFT → spatial height field ─────────────────────────── */
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    fftw_free(freq);

    /* normalise FFTW output (it is not normalised by default) */
    double norm = 1.0 / sqrt((double)nrows * ncols);

    /* compute variance for rescaling */
    double sum2 = 0.0;
    long  n_ocean = 0;
    for (long i = 0; i < npx; i++) {
        spatial[i] *= norm;
        if (!Rast_is_d_null_value(&mask->data[i])) {
            sum2 += spatial[i] * spatial[i];
            n_ocean++;
        }
    }
    double stddev = (n_ocean > 0) ? sqrt(sum2 / n_ocean) : 1.0;
    if (stddev < 1e-12) stddev = 1.0;
    /* scale so 3σ ≈ turb_amp */
    double scale = p.turb_amp / (3.0 * stddev);

    /* ── shore distance (for coastal foam boost) ─────────────────────── */
    float *sdist = shore_dist_pixels(mask, nprocs);

    /* ── apply turbulence + foam to depth index ─────────────────────── */
    Rbuf *out = rbuf_alloc(nrows, ncols);

    double min_val = DBL_MAX;

#pragma omp parallel for schedule(static) reduction(min:min_val)
    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            long i = (long)r * ncols + c;
            if (Rast_is_d_null_value(&mask->data[i])) {
                Rast_set_d_null_value(&out->data[i], 1);
                continue;
            }

            double h    = spatial[i] * scale;   /* turbulence value ±turb_amp */

            /*
             * Foam detection (0AD water_high.fs lines 343-356):
             *   foam appears at CRESTS (most negative h, since h<0 lightens)
             *   excess = max(0, −h/turb_amp − thresh)
             *   foam   = excess^2 · scale · (1 + 4·exp(−dist/shore_decay))
             */
            double excess = -h / p.turb_amp - p.foam_thresh;
            double foam   = 0.0;
            if (excess > 0.0) {
                double coast = 1.0 + 4.0 * exp(-(double)sdist[i]
                                                / p.shore_decay_px);
                foam = excess * excess * p.foam_scale * coast;
            }

            double val = depth_idx->data[i] + h - foam;
            out->data[i] = val;
            if (val < min_val) min_val = val;
        }
    }

    fftw_free(spatial);
    G_free(sdist);

    *lo_clamp_out = min_val - 1.0;

    G_verbose_message(_("Turbulence height field stddev=%.3f scale=%.4f "
                        "min_val=%.1f"), stddev, scale, min_val);

    return out;
}
