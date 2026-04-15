/****************************************************************************
 *
 * MODULE:       i.ocean  —  ocean_mask.c
 * PURPOSE:      Create binary ocean masks from input rasters.
 *
 *               make_mask_from_value():
 *                 Cells equal to ocean_val → 1.0; all others → null.
 *
 *               make_mask_from_depth():
 *                 Cells >= depth_min → 1.0; all others → null.
 *
 *               Both functions operate with OpenMP-parallel row loops.
 *
 ****************************************************************************/

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>
#include <omp.h>

#include "local_proto.h"

/* ── mask from explicit ocean-value ─────────────────────────────────── */

Rbuf *make_mask_from_value(const Rbuf *input, double ocean_val,
                            int nrows, int ncols)
{
    Rbuf *mask = rbuf_alloc(nrows, ncols);
    long npx = (long)nrows * ncols;

#pragma omp parallel for schedule(static)
    for (long i = 0; i < npx; i++) {
        if (!Rast_is_d_null_value(&input->data[i]) &&
            input->data[i] == ocean_val)
            mask->data[i] = 1.0;
        else
            Rast_set_d_null_value(&mask->data[i], 1);
    }
    return mask;
}

/* ── mask from depth threshold ───────────────────────────────────────── */

Rbuf *make_mask_from_depth(const Rbuf *depth, double depth_min,
                            int nrows, int ncols)
{
    Rbuf *mask = rbuf_alloc(nrows, ncols);
    long npx = (long)nrows * ncols;

#pragma omp parallel for schedule(static)
    for (long i = 0; i < npx; i++) {
        if (!Rast_is_d_null_value(&depth->data[i]) &&
            depth->data[i] >= depth_min)
            mask->data[i] = 1.0;
        else
            Rast_set_d_null_value(&mask->data[i], 1);
    }
    return mask;
}
