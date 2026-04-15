/****************************************************************************
 *
 * MODULE:       i.ocean  —  waves.c
 * PURPOSE:      Scale-aware sinusoidal wave texture.
 *
 *               Generates a raster of wave-height values by superimposing
 *               a primary cross-wave (sin × cos) and a 40 % diagonal
 *               harmonic, following the 0AD water_high.fs two-layer normal
 *               map approach (water_high.fs lines 149–176):
 *
 *                 primary = A · sin(2π/Tc · col) · cos(2π/Tr · row)
 *                 diag    = 0.4A · sin(2π/Td · (col + row))
 *                 wave    = primary + diag
 *
 *               Wave period and amplitude are chosen from five scale tiers
 *               based on mean pixel size (metres).  All ocean pixels are
 *               computed with an OpenMP-parallel row loop.
 *
 *               Returns the wave-amplitude raster and sets *amp_out to the
 *               amplitude value (used by the colour-rule builder to extend
 *               the clamp range).
 *
 ****************************************************************************/

#include <math.h>
#include <omp.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>

#include "local_proto.h"

Rbuf *wave_pattern(const Rbuf *mask, double pixel_m,
                   double *amp_out, int nprocs)
{
    int nrows = mask->nrows, ncols = mask->ncols;

    /* ── scale-tier selection ────────────────────────────────────────── */
    int period_c, period_r, diag_period;
    double amp;

    if (pixel_m < 100.0) {
        period_c = (int)(30.0  / pixel_m); if (period_c < 4) period_c = 4;
        period_r = (int)(20.0  / pixel_m); if (period_r < 4) period_r = 4;
        amp      = 12.0;
    } else if (pixel_m < 1000.0) {
        period_c = (int)(600.0  / pixel_m); if (period_c < 5) period_c = 5;
        period_r = (int)(400.0  / pixel_m); if (period_r < 5) period_r = 5;
        amp      = 25.0;
    } else if (pixel_m < 10000.0) {
        period_c = (int)(6000.0  / pixel_m); if (period_c < 5) period_c = 5;
        period_r = (int)(4000.0  / pixel_m); if (period_r < 5) period_r = 5;
        amp      = 35.0;
    } else if (pixel_m < 100000.0) {
        period_c = (int)(60000.0 / pixel_m); if (period_c < 4) period_c = 4;
        period_r = (int)(40000.0 / pixel_m); if (period_r < 4) period_r = 4;
        amp      = 30.0;
    } else {
        period_c = 10;
        period_r =  7;
        amp      = 20.0;
    }
    diag_period = (int)(period_c * 1.6);
    if (diag_period < 4) diag_period = 4;

    *amp_out = amp;

    G_verbose_message(_("Wave texture: period_col=%d period_row=%d amplitude=%.0f"),
                      period_c, period_r, amp);

    double two_pi = 2.0 * M_PI;
    double fc  = two_pi / period_c;
    double fr  = two_pi / period_r;
    double fd  = two_pi / diag_period;
    double amp04 = amp * 0.4;

    Rbuf *out = rbuf_alloc(nrows, ncols);

    omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            if (Rast_is_d_null_value(&RBUF_CELL(mask, r, c))) {
                Rast_set_d_null_value(&RBUF_CELL(out, r, c), 1);
                continue;
            }
            /* col() and row() in r.mapcalc are 1-based; replicate that */
            double col1 = c + 1.0;
            double row1 = r + 1.0;
            double v = amp   * sin(fc * col1) * cos(fr * row1)
                     + amp04 * sin(fd * (col1 + row1));
            RBUF_CELL(out, r, c) = v;
        }
    }
    return out;
}
