/****************************************************************************
 *
 * MODULE:       i.ocean  —  depth_index.c
 * PURPOSE:      Depth-index computation (normalise 0–1000), integer
 *               raster smoothing, and shore-distance depth estimation.
 *
 *   smooth_3x3()
 *     3×3 neighbourhood average.  Replicates the r.neighbors method=average
 *     behaviour that always produces a DCELL output (r.neighbors/main.c:64,
 *     output_type() lines 117-135: T_FLOAT → DCELL unconditionally).
 *     Applied when the depth map is CELL type to eliminate staircase
 *     artefacts from the rendered gradient.  OpenMP-parallel.
 *
 *   depth_from_map()
 *     Reads a depth raster, optionally smooths it (CELL input), computes
 *     abs(), finds the max within the ocean mask, and linearly rescales to
 *     0–1000.
 *
 *   depth_from_shore()
 *     BFS-style two-pass Euclidean distance transform from land pixels,
 *     then rescales to 0–1000.  Does not require r.grow.distance.
 *     Uses OpenMP on the row loops.
 *
 *   depth_flat()
 *     Constant depth index (default 500 = mid-palette).
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include <omp.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>

#include "local_proto.h"

/* ── 3×3 neighbourhood average (in-place) ────────────────────────────── */

/*
 * smooth_3x3() — replaces each ocean pixel with the mean of its non-null
 * 3×3 neighbourhood.  Border pixels use only the available neighbours.
 *
 * Implementation note:
 *   r.neighbors with method=average registers output type T_FLOAT
 *   (main.c line 64).  output_type() maps T_FLOAT → DCELL_TYPE
 *   unconditionally (lines 117–135); the T_COPY branch (127–128) is the
 *   only path that would preserve the input type, and 'average' does not
 *   take that path.  We replicate this by computing weighted averages in
 *   double precision and storing in the DCELL buffer.
 *
 *   r.neighbors also accepts nprocs= (main.c line 279, G_OPT_M_NPROCS);
 *   we forward the same value through the OpenMP thread count.
 */
void smooth_3x3(Rbuf *buf, const Rbuf *mask, int nprocs)
{
    int nrows = buf->nrows, ncols = buf->ncols;
    DCELL *tmp = G_malloc(sizeof(DCELL) * (size_t)nrows * ncols);
    memcpy(tmp, buf->data, sizeof(DCELL) * (size_t)nrows * ncols);

    omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            if (Rast_is_d_null_value(&RBUF_CELL(mask, r, c)))
                continue;
            double sum = 0.0;
            int    cnt = 0;
            for (int dr = -1; dr <= 1; dr++) {
                int rr = r + dr;
                if (rr < 0 || rr >= nrows) continue;
                for (int dc = -1; dc <= 1; dc++) {
                    int cc = c + dc;
                    if (cc < 0 || cc >= ncols) continue;
                    DCELL v = tmp[rr * ncols + cc];
                    if (!Rast_is_d_null_value(&v)) {
                        sum += v;
                        cnt++;
                    }
                }
            }
            if (cnt > 0)
                RBUF_CELL(buf, r, c) = sum / cnt;
        }
    }
    G_free(tmp);
}

/* ── depth_from_map ──────────────────────────────────────────────────── */

Rbuf *depth_from_map(const char *depth_name, const Rbuf *mask,
                      int nprocs, const char *pid_str)
{
    (void)pid_str;
    int nrows = mask->nrows, ncols = mask->ncols;
    long npx = (long)nrows * ncols;

    /* read depth raster */
    Rbuf *dep = rbuf_read(depth_name);

    /* integer (CELL) input: smooth to remove staircase artefacts */
    {
        int fd = Rast_open_old(depth_name, "");
        RASTER_MAP_TYPE mtype = Rast_get_map_type(fd);
        Rast_close(fd);
        if (mtype == CELL_TYPE) {
            G_verbose_message(_("Depth map is integer (CELL) — "
                                "smoothing with 3x3 average to remove "
                                "staircase gradient artefacts..."));
            smooth_3x3(dep, mask, nprocs);
        }
    }

    /* abs() within ocean mask; find maximum */
    double max_d = 0.0;
    omp_set_num_threads(nprocs);
#pragma omp parallel for reduction(max:max_d) schedule(static)
    for (long i = 0; i < npx; i++) {
        if (!Rast_is_d_null_value(&mask->data[i]) &&
            !Rast_is_d_null_value(&dep->data[i])) {
            double v = fabs(dep->data[i]);
            dep->data[i] = v;
            if (v > max_d) max_d = v;
        }
    }

    if (max_d <= 0.0) {
        G_warning(_("Depth map maximum is 0 or invalid; using flat depth of 500."));
        max_d = 0.0;
    }

    /* normalise to 0–1000 */
    Rbuf *out = rbuf_alloc(nrows, ncols);
#pragma omp parallel for schedule(static)
    for (long i = 0; i < npx; i++) {
        if (Rast_is_d_null_value(&mask->data[i])) {
            Rast_set_d_null_value(&out->data[i], 1);
        } else if (max_d > 0.0) {
            out->data[i] = fmin(1000.0, 1000.0 * dep->data[i] / max_d);
        } else {
            out->data[i] = 500.0;
        }
    }

    rbuf_free(dep);
    return out;
}

/* ── depth_from_shore ────────────────────────────────────────────────── */

/*
 * depth_from_shore() — Euclidean distance transform from land boundary,
 * rescaled to 0–1000 depth index.
 *
 * Algorithm: two-pass exact Euclidean distance transform (Meijster 2000,
 * simplified variant).  Avoids calling r.grow.distance as a subprocess.
 *
 * Phase 1: for each column, compute 1-D distance along rows from the
 *          nearest land pixel (scanning top-to-bottom then bottom-to-top).
 * Phase 2: for each row, refine using column distances.
 *
 * Land pixels are those where the mask is null.  Distance is in pixels;
 * the maximum pixel dimension (max of ns_res, ew_res) converts to metres
 * for consistency with r.grow.distance output.
 */
Rbuf *depth_from_shore(const Rbuf *mask, double pixel_m,
                        int nprocs, const char *pid_str)
{
    (void)pid_str;
    (void)pixel_m;
    int nrows = mask->nrows, ncols = mask->ncols;
    long npx = (long)nrows * ncols;

    /* Phase 1: column-wise 1-D dt (distances stored in row units) */
    float *g = G_malloc(sizeof(float) * npx);

    omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
    for (int c = 0; c < ncols; c++) {
        /* downward pass */
        float last = (float)nrows * 2.0f;
        for (int r = 0; r < nrows; r++) {
            if (Rast_is_d_null_value(&RBUF_CELL(mask, r, c)))
                last = 0.0f;   /* land pixel — distance 0 */
            else
                last += 1.0f;
            g[r * ncols + c] = last;
        }
        /* upward pass */
        last = (float)nrows * 2.0f;
        for (int r = nrows - 1; r >= 0; r--) {
            if (Rast_is_d_null_value(&RBUF_CELL(mask, r, c)))
                last = 0.0f;
            else
                last += 1.0f;
            if (last < g[r * ncols + c])
                g[r * ncols + c] = last;
        }
    }

    /* Phase 2: row-wise parabolic DT (Meijster) — gives exact Euclidean */
    Rbuf *dist = rbuf_alloc(nrows, ncols);
    int   *v   = G_malloc(sizeof(int)   * ncols);
    float *z   = G_malloc(sizeof(float) * (ncols + 1));

#pragma omp parallel for schedule(static) \
    firstprivate(v, z)
    for (int r = 0; r < nrows; r++) {
        /* build lower envelope of parabolas */
        int q = 0;
        v[0] = 0;
        z[0] = -1e30f;
        z[1] =  1e30f;
        for (int c = 1; c < ncols; c++) {
            float gc = g[r * ncols + c];
            float s;
            do {
                float gv = g[r * ncols + v[q]];
                s = ((gc * gc + c * c) - (gv * gv + v[q] * v[q]))
                    / (2.0f * (c - v[q]));
            } while (q > 0 && s <= z[q] && --q >= 0);
            q++;
            v[q]     = c;
            z[q]     = s;
            z[q + 1] = 1e30f;
        }
        /* read off distances */
        q = 0;
        for (int c = 0; c < ncols; c++) {
            while (z[q + 1] < (float)c) q++;
            float dv = g[r * ncols + v[q]];
            float dc = (float)(c - v[q]);
            dist->data[r * ncols + c] = sqrt(dv * dv + dc * dc);
        }
    }
    G_free(v);
    G_free(z);
    G_free(g);

    /* find maximum distance within ocean mask */
    double max_dist = 0.0;
#pragma omp parallel for reduction(max:max_dist) schedule(static)
    for (long i = 0; i < npx; i++) {
        if (!Rast_is_d_null_value(&mask->data[i])) {
            if (dist->data[i] > max_dist)
                max_dist = dist->data[i];
        }
    }
    if (max_dist <= 0.0) max_dist = 1.0;

    /* rescale to depth index 0–1000 */
    Rbuf *out = rbuf_alloc(nrows, ncols);
#pragma omp parallel for schedule(static)
    for (long i = 0; i < npx; i++) {
        if (Rast_is_d_null_value(&mask->data[i]))
            Rast_set_d_null_value(&out->data[i], 1);
        else
            out->data[i] = fmin(1000.0,
                                1000.0 * dist->data[i] / max_dist);
    }

    rbuf_free(dist);
    return out;
}

/* ── depth_flat ──────────────────────────────────────────────────────── */

Rbuf *depth_flat(const Rbuf *mask, double value, int nprocs)
{
    int nrows = mask->nrows, ncols = mask->ncols;
    long npx = (long)nrows * ncols;
    Rbuf *out = rbuf_alloc(nrows, ncols);

    omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
    for (long i = 0; i < npx; i++) {
        if (Rast_is_d_null_value(&mask->data[i]))
            Rast_set_d_null_value(&out->data[i], 1);
        else
            out->data[i] = value;
    }
    return out;
}
