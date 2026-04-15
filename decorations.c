/****************************************************************************
 *
 * MODULE:       i.ocean  —  decorations.c
 * PURPOSE:      Scatter maritime decorations onto the depth-index buffer.
 *
 * OVERVIEW
 * ────────
 * Each decoration type is defined as a compact pixel-art pattern:
 *   • a list of (Δrow, Δcol, depth_index_delta) pixel offsets
 *   • an allowed depth-index range [min_idx, max_idx]
 *   • a placement density in decorations per km²
 *
 * The placement engine:
 *   1. Counts ocean pixels and derives the ocean area in km².
 *   2. For each type, computes how many items to place.
 *   3. Picks random ocean pixels using a seeded PRNG (reproducible).
 *   4. Stamps each pattern, clamping depth-index values to [lo_clamp, 1100].
 *
 * DECORATION TYPES
 * ─────────────────
 *
 *  WHALE  — 13 px elongated silhouette.
 *    Deep water (idx 400–900).  Body is bright (surface sheen, -180),
 *    tail flukes slightly lighter (-80).
 *    Density: 0.0008 / km²  (rare, large — one per ~1250 km²)
 *
 *  DOLPHIN  — 5 px arc.
 *    Medium depth (idx 100–700).  Very bright (-150).
 *    Density: 0.002 / km²
 *
 *  FISH_SCHOOL  — 8 scattered pixels in a 6×6 area.
 *    Any ocean depth (idx 0–800).  Mottled bright (-60 to -90).
 *    Density: 0.008 / km²
 *
 *  OCTOPUS  — 9 px radial body + tentacles.
 *    Medium depth (idx 50–600).  Alternating highlight/shadow (±55).
 *    Density: 0.003 / km²
 *
 *  TREASURE_CHEST  — 4 px bright cluster.
 *    Shallow to medium (idx 0–500).  Very bright (-550): "gold shimmer"
 *    clamped to near-zero by colour table.
 *    Density: 0.0001 / km²  (very rare)
 *
 *  FLOATING_CRATE  — 6 px rectangle.
 *    Near-shore (idx 0–350).  Bright (-250): weathered wood pale colour.
 *    Density: 0.001 / km²
 *
 *  SUNKEN_SHIP  — 20 px elongated hull.
 *    Medium-deep (idx 200–900).  Dark keel (+90), lighter deck/mast (-80):
 *    shadow-and-highlight creates the underwater silhouette effect.
 *    Density: 0.0003 / km²
 *
 * DEPTH-INDEX CONVENTION
 * ──────────────────────
 *   0    = shore colour (lightest teal / turquoise)
 *   1000 = abyss colour (darkest blue)
 *   delta < 0  → brightens (pushes toward shore colour)
 *   delta > 0  → darkens  (pushes toward abyss colour)
 *
 *   Foam/wave effects may have already pushed values below 0; the colour
 *   table clamps these to the shore colour, so very bright features
 *   (treasure, floating crate) reliably appear as the lightest ocean hue.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>

#include "local_proto.h"

/* ── pixel offset entry ──────────────────────────────────────────────── */

typedef struct { int dr, dc; double delta; } PxOff;

/* ── decoration definition ───────────────────────────────────────────── */

typedef struct {
    const char  *name;
    int          n_px;
    PxOff        px[24];
    double       min_idx;   /* depth index lower bound for placement */
    double       max_idx;   /* depth index upper bound               */
    double       density;   /* items per km²                         */
} DecoType;

/* ── pattern definitions ─────────────────────────────────────────────── */

static const DecoType DECORATIONS[] = {

    /* ── whale (13 pixels) ─────────────────────────────────────────── */
    { "whale", 13, {
        { 0,-4,-80 }, { 0,-3,-180}, { 0,-2,-180}, { 0,-1,-180},
        { 0, 0,-180}, { 0, 1,-180}, { 0, 2,-180}, { 0, 3,-80 },
        {-1,-1,-100}, { 1,-1,-100},
        { 0, 4,-50 }, {-1, 4,-80 }, { 1, 4,-80 }  /* tail flukes */
      },
      400.0, 900.0, 0.0008 },

    /* ── dolphin (5 pixels) ─────────────────────────────────────────── */
    { "dolphin", 5, {
        { 0, 0,-150}, {-1,-1,-100}, {-1, 0,-120},
        { 0, 1,-100}, { 1, 1,-80 }
      },
      100.0, 700.0, 0.002 },

    /* ── fish school (8 scattered pixels) ──────────────────────────── */
    { "fish_school", 8, {
        { 0, 0,-90 }, { 0, 2,-60 }, { 1,-1,-75 }, {-1, 1,-65 },
        { 2, 1,-85 }, {-1, 3,-70 }, { 1, 3,-55 }, { 3, 0,-80 }
      },
      0.0, 800.0, 0.008 },

    /* ── octopus (9 pixels) ─────────────────────────────────────────── */
    { "octopus", 9, {
        { 0, 0,-55 },              /* mantle centre */
        {-1, 0,-40 }, { 1, 0,-40 },{ 0,-1,-40 },{ 0, 1,-40 },  /* body */
        {-2,-1, 55 }, {-2, 1, 55 },/* tentacles — darker against bright */
        { 2,-1, 55 }, { 2, 1, 55 }
      },
      50.0, 600.0, 0.003 },

    /* ── treasure chest (4 pixels) ──────────────────────────────────── */
    { "treasure", 4, {
        { 0, 0,-550}, { 0, 1,-550},
        { 1, 0,-450}, { 1, 1,-450}
      },
      0.0, 500.0, 0.0001 },

    /* ── floating crate (6 pixels) ──────────────────────────────────── */
    { "floating_crate", 6, {
        { 0,-1,-250}, { 0, 0,-250}, { 0, 1,-250}, { 0, 2,-250},
        { 1,-1,-200}, { 1, 2,-200}                /* outline sides */
      },
      0.0, 350.0, 0.001 },

    /* ── sunken ship (20 pixels) ─────────────────────────────────────── */
    { "sunken_ship", 20, {
        /* mast  */  { 0, 0,-80 },
        /* deck  */  { 1,-3,-60 },{ 1,-2,-60 },{ 1,-1,-60 },
                     { 1, 0,-60 },{ 1, 1,-60 },{ 1, 2,-60 },{ 1, 3,-60 },
        /* hull  */  { 2,-4, 90 },{ 2,-3, 90 },{ 2,-2, 90 },{ 2,-1, 90 },
                     { 2, 0, 90 },{ 2, 1, 90 },{ 2, 2, 90 },{ 2, 3, 90 },
                     { 2, 4, 90 },
        /* keel  */  { 3,-3, 70 },{ 3, 0, 70 },{ 3, 3, 70 }
      },
      200.0, 900.0, 0.0003 }
};

#define N_DECO_TYPES ((int)(sizeof(DECORATIONS) / sizeof(DECORATIONS[0])))

/* ── stamp one decoration onto the work buffer ───────────────────────── */

static void stamp(Rbuf *work, const Rbuf *mask,
                  int center_r, int center_c,
                  const DecoType *d, double lo_clamp)
{
    int nrows = work->nrows, ncols = work->ncols;
    for (int p = 0; p < d->n_px; p++) {
        int r = center_r + d->px[p].dr;
        int c = center_c + d->px[p].dc;
        if (r < 0 || r >= nrows || c < 0 || c >= ncols) continue;
        if (Rast_is_d_null_value(&RBUF_CELL(mask, r, c)))  continue;
        double v = RBUF_CELL(work, r, c) + d->px[p].delta;
        if (v < lo_clamp) v = lo_clamp;
        RBUF_CELL(work, r, c) = v;
    }
}

/* ── place_decorations ───────────────────────────────────────────────── */

void place_decorations(Rbuf *work, const Rbuf *mask,
                        double pixel_m, double lo_clamp)
{
    int  nrows  = work->nrows, ncols = work->ncols;
    long npx    = (long)nrows * ncols;

    /* collect ocean pixel indices into an array for random sampling */
    long *ocean_idx = G_malloc(sizeof(long) * npx);
    long  n_ocean   = 0;
    for (long i = 0; i < npx; i++)
        if (!Rast_is_d_null_value(&mask->data[i]))
            ocean_idx[n_ocean++] = i;

    if (n_ocean == 0) {
        G_free(ocean_idx);
        return;
    }

    /* ocean area in km² */
    double pix_km2 = (pixel_m * pixel_m) * 1e-6;
    double ocean_km2 = n_ocean * pix_km2;
    G_verbose_message(_("Ocean area for decorations: %.0f km²"), ocean_km2);

    /* seeded PRNG (reproducible results) */
    unsigned int seed = 0xDEADBEEFu;

    int total_placed = 0;

    for (int t = 0; t < N_DECO_TYPES; t++) {
        const DecoType *d = &DECORATIONS[t];

        /* number of items proportional to ocean area */
        int count = (int)ceil(ocean_km2 * d->density);
        if (count < 1) count = 1;

        int placed = 0;
        int tries  = 0;
        int max_tries = count * 20;

        while (placed < count && tries < max_tries) {
            tries++;
            long pick = (long)((double)n_ocean * (rand_r(&seed)
                                / ((double)RAND_MAX + 1.0)));
            if (pick >= n_ocean) pick = n_ocean - 1;

            long idx   = ocean_idx[pick];
            int  r     = (int)(idx / ncols);
            int  c     = (int)(idx % ncols);
            double cur = work->data[idx];

            if (cur < d->min_idx || cur > d->max_idx)
                continue;

            stamp(work, mask, r, c, d, lo_clamp);
            placed++;
        }
        total_placed += placed;
        G_verbose_message(_("  %s: %d placed (%.0f km² at density %.4f/km²)"),
                          d->name, placed, ocean_km2, d->density);
    }

    G_message(_("Maritime decorations: %d items placed across %.0f km² of ocean."),
              total_placed, ocean_km2);
    G_free(ocean_idx);
}
