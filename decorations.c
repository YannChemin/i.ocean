/****************************************************************************
 *
 * MODULE:       i.ocean  —  decorations.c
 * PURPOSE:      Scatter maritime decorations onto the depth-index buffer
 *               and collect placement records for the vector layer.
 *
 * OVERVIEW
 * ────────
 * Each decoration type is defined by:
 *   • a list of (Δrow, Δcol, depth_index_delta) pixel offsets (raster art)
 *   • an allowed depth-index range [min_idx, max_idx]
 *   • a placement density in items per km²
 *
 * The placement engine:
 *   1. Counts ocean pixels and derives the ocean area in km².
 *   2. For each type computes how many items to place (capped at
 *      MAX_PER_TYPE to keep the vector map manageable).
 *   3. Picks random ocean pixels with a seeded PRNG (reproducible).
 *   4. Stamps each raster pattern, clamping depth-index values.
 *   5. Records the placement for the vector silhouette layer.
 *
 * REAL-WORLD SIZES AND HOMOGENEITY RULE
 * ──────────────────────────────────────
 * Real-world sizes (deco_real_length / deco_real_width) are used by the
 * vector layer in vector_deco.c.  Biologically/historically small objects
 * (treasure chest ~1 m, floating crate ~1.2 m) are bumped up to the
 * dolphin minimum (3 m) so all item types share a homogeneous minimum
 * ground size that remains visible at typical cartographic scales.
 *
 *   Whale (blue whale avg):     28 m × 7 m
 *   Dolphin (bottlenose avg):    3 m × 0.7 m
 *   Fish school:                15 m diameter
 *   Octopus (arm span):          3 m × 3 m
 *   Treasure chest (bumped):     3 m × 2 m   (real ≈ 1 m)
 *   Floating crate (bumped):     3 m × 3 m   (real ≈ 1.2 m)
 *   Sunken ship (trireme):      37 m × 7 m   (0AD hele_ship_bireme.dae)
 *
 * DECORATION TYPES (raster pixel art)
 * ─────────────────────────────────────
 *
 *  WHALE  — 13 px elongated silhouette.
 *    Deep water (idx 400–900).  Body is bright (surface sheen, −180),
 *    tail flukes slightly lighter (−80).
 *    Density: 0.0008 / km²  (rare, large — one per ~1250 km²)
 *
 *  DOLPHIN  — 5 px arc.
 *    Medium depth (idx 100–700).  Very bright (−150).
 *    Density: 0.002 / km²
 *
 *  FISH_SCHOOL  — 8 scattered pixels in a 6×6 area.
 *    Any ocean depth (idx 0–800).  Mottled bright (−60 to −90).
 *    Density: 0.008 / km²
 *
 *  OCTOPUS  — 9 px radial body + tentacles.
 *    Medium depth (idx 50–600).  Alternating highlight/shadow (±55).
 *    Density: 0.003 / km²
 *
 *  TREASURE_CHEST  — 4 px bright cluster.
 *    Shallow to medium (idx 0–500).  Very bright (−550).
 *    Density: 0.0001 / km²  (very rare)
 *
 *  FLOATING_CRATE  — 6 px rectangle.
 *    Near-shore (idx 0–350).  Bright (−250).
 *    Density: 0.001 / km²
 *
 *  SUNKEN_SHIP  — 20 px elongated hull.
 *    Medium-deep (idx 200–900).  Dark keel (+90), lighter deck/mast (−80).
 *    Density: 0.0003 / km²
 *
 * DEPTH-INDEX CONVENTION
 * ──────────────────────
 *   0    = shore colour (lightest teal / turquoise)
 *   1000 = abyss colour (darkest blue)
 *   delta < 0  → brightens (pushes toward shore colour)
 *   delta > 0  → darkens  (pushes toward abyss colour)
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

/* ── real-world decoration sizes (see homogeneity rule in header) ──────── */

const char *deco_type_name[N_DECO_TYPES] = {
    "whale", "dolphin", "fish_school", "octopus",
    "treasure", "floating_crate", "sunken_ship",
    /* war vessels */
    "hunter_sub",   "leviathan_sub", "barracuda_boat", "kraken_carrier",
    "phantom_sub",  "titan_barge",   "wraith_sub",    "hydra_sub",
    "aegis_sub",    "oblivion_sub",
    /* new biological / object types */
    "orca", "sperm_whale", "shark", "barrel_cluster",
    /* merchant vessels */
    "crude_tanker", "lng_carrier", "bulk_carrier",
    "product_tanker", "container_ship", "cargo_ship"
};

/*
 * Biologically/historically accurate lengths in metres.
 * Small objects (treasure, crate) are bumped to 3 m (dolphin minimum) so
 * all types share a homogeneous minimum ground footprint on the map.
 */
const double deco_real_length[N_DECO_TYPES] = {
   110.0,   /* whale (blue whale): y=+55 to fluke tip y=-55             */
    57.0,   /* dolphin: y=+32 (rostrum) to tail y=-25                   */
    22.0,   /* fish school: y extent ~22 m (12-position group layout)   */
     3.0,   /* octopus: arm-span (bumped to 3 m minimum)                */
     8.0,   /* treasure chest: body+lid, y=-2 to y=+6                   */
     6.0,   /* crate cluster: 8-crate stack, y=-1 to y=+5               */
    37.0,   /* sunken ship: trireme hull (0AD hele_ship_bireme.dae)      */
    /* war vessels — bow_y to stern_y */
    78.0,   /* hunter_sub:     y=+40 to y=-38                           */
    95.0,   /* leviathan_sub:  y=+50 to y=-45                           */
    50.0,   /* barracuda_boat: y=+30 to y=-20                           */
    65.0,   /* kraken_carrier: y=+35 to y=-30                           */
    60.0,   /* phantom_sub:    y=+32 to y=-28                           */
    47.0,   /* titan_barge:    y=+25 to y=-22                           */
    27.0,   /* wraith_sub:     y=+15 to y=-12                           */
    55.0,   /* hydra_sub:      prong y=+25 to y=-30                     */
    53.0,   /* aegis_sub:      y=+28 to y=-25                           */
   115.0,   /* oblivion_sub:   y=+60 to y=-55                           */
    /* new biological / object types */
    67.0,   /* orca:        y=+35 to fluke y=-32                        */
    95.0,   /* sperm_whale: y=+45 to fluke y=-50                        */
    62.0,   /* shark:       y=+30 to tail y=-32                         */
     7.0,   /* barrel_cluster: cluster y=-2 to y=+5                     */
    /* merchant vessels — full length bow-to-stern */
   250.0,   /* crude_tanker  (Aframax):  y=+125 to y=-125               */
   290.0,   /* lng_carrier   (Moss):     y=+145 to y=-145               */
   300.0,   /* bulk_carrier  (Capesize): y=+150 to y=-150               */
   180.0,   /* product_tanker:           y=+90  to y=-90                */
   300.0,   /* container_ship (Neo-Panamax): y=+150 to y=-150           */
   140.0,   /* cargo_ship    (General):  y=+70  to y=-70                */
};

const double deco_real_width[N_DECO_TYPES] = {
    30.0,   /* whale: ±15 m pectoral fin span                           */
    16.0,   /* dolphin: ±8 m tail fluke span                            */
    23.0,   /* fish school: ±11.5 m x extent (12-position layout)       */
     3.0,   /* octopus: arm span (symmetric)                            */
     8.0,   /* treasure chest: ±4 m body width                         */
    12.0,   /* crate cluster: x=-6 to x=+6 bottom row                  */
     7.0,   /* sunken ship beam (trireme ~7 m)                          */
    /* war vessels — 2 × max half-beam */
     8.0,   /* hunter_sub:     ±4 m hull                                */
    20.0,   /* leviathan_sub:  ±10 m hull                               */
    12.0,   /* barracuda_boat: ±6 m fin span                            */
    18.0,   /* kraken_carrier: ±9 m hull                                */
    10.0,   /* phantom_sub:    ±5 m hull                                */
    16.0,   /* titan_barge:    ±8 m hull                                */
     5.0,   /* wraith_sub:     ±2.5 m hull                              */
    10.0,   /* hydra_sub:      ±5 m at prongs                           */
    24.0,   /* aegis_sub:      ±12 m sensor arms                        */
    24.0,   /* oblivion_sub:   ±12 m hull                               */
    /* new biological / object types */
    20.0,   /* orca:           ±10 m fluke span                         */
    20.0,   /* sperm_whale:    ±10 m hull                               */
    24.0,   /* shark:          ±12 m pectoral fin span                  */
    12.0,   /* barrel_cluster: x=-6 to x=+6                            */
    /* merchant vessels — beam (max x-extent × 2) */
    60.0,   /* crude_tanker:   ±30 m hull beam                         */
    70.0,   /* lng_carrier:    ±35 m hull beam                         */
    64.0,   /* bulk_carrier:   ±32 m hull beam                         */
    44.0,   /* product_tanker: ±22 m hull beam                         */
    60.0,   /* container_ship: ±30 m hull beam                         */
    32.0,   /* cargo_ship:     ±16 m hull beam                         */
};

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

/* ── max vector items per type (keeps vector map manageable) ────────── */
#define MAX_PER_TYPE  2000

/* ── pattern definitions ─────────────────────────────────────────────── */

static const DecoType DECORATIONS[N_DECO_TYPES] = {

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
      200.0, 900.0, 0.0003 },

    /* ── hunter-class attack submarine (10 px) — 78 m × 8 m ─────────────── */
    { "hunter_sub", 10, {
        /* hull */   { 0,-4,-60 },{ 0,-3,-90 },{ 0,-2,-90 },{ 0,-1,-90 },
                     { 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-90 },
                     { 0, 4,-60 },
        /* tower */  {-1, 1,-110 }
      },
      300.0, 900.0, 0.0002 },

    /* ── leviathan heavy battlesub (13 px) — 95 m × 20 m ────────────────── */
    { "leviathan_sub", 13, {
        /* hull */   { 0,-5,-50 },{ 0,-4,-90 },{ 0,-3,-90 },{ 0,-2,-90 },
                     { 0,-1,-90 },{ 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },
                     { 0, 3,-90 },{ 0, 4,-90 },{ 0, 5,-50 },
        /* silos */  {-1, 0,-110 },{ 1, 0,-110 }
      },
      500.0, 900.0, 0.0001 },

    /* ── barracuda fast strike boat (9 px) — 50 m × 12 m ────────────────── */
    { "barracuda_boat", 9, {
        /* hull */   { 0,-3,-60 },{ 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },
                     { 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-60 },
        /* fins */   {-1,-1,-70 },{ 1,-1,-70 }
      },
      100.0, 700.0, 0.0003 },

    /* ── kraken drone carrier sub (13 px) — 65 m × 18 m ─────────────────── */
    { "kraken_carrier", 13, {
        /* hull */   { 0,-4,-50 },{ 0,-3,-90 },{ 0,-2,-90 },{ 0,-1,-90 },
                     { 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-90 },
                     { 0, 4,-50 },
        /* hangar */ {-2,-1,-80 },{ 2,-1,-80 },{ 2, 0,-80 },{-2, 0,-80 }
      },
      300.0, 800.0, 0.00008 },

    /* ── phantom stealth sub (9 px) — 60 m × 10 m ───────────────────────── */
    { "phantom_sub", 9, {
        /* smooth hull, no protrusions */
                     { 0,-4,-50 },{ 0,-3,-90 },{ 0,-2,-90 },{ 0,-1,-90 },
                     { 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-90 },
                     { 0, 4,-50 }
      },
       50.0, 900.0, 0.0002 },

    /* ── titan torpedo barge (12 px) — 47 m × 16 m ──────────────────────── */
    { "titan_barge", 12, {
        /* blocky front face */
                     {-1,-3,-70 },{ 0,-3,-90 },{ 1,-3,-70 },
                     {-1,-2,-70 },{ 0,-2,-90 },{ 1,-2,-70 },
                     {-1,-1,-70 },{ 0,-1,-90 },{ 1,-1,-70 },
        /* stern */  { 0, 0,-70 },{ 0, 1,-60 },{ 0, 2,-50 }
      },
      100.0, 600.0, 0.0001 },

    /* ── wraith mini sub (7 px) — 27 m × 5 m ────────────────────────────── */
    { "wraith_sub", 7, {
        /* hull */   { 0,-2,-60 },{ 0,-1,-90 },{ 0, 0,-90 },
                     { 0, 1,-90 },{ 0, 2,-60 },
        /* hatch */  {-1, 0,-110 },{ 1, 0,-110 }
      },
        0.0, 800.0, 0.0004 },

    /* ── hydra multi-head sub (12 px) — 55 m × 10 m ─────────────────────── */
    { "hydra_sub", 12, {
        /* central hull */
                     { 0,-3,-60 },{ 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },
        /* split bow */
                     {-1, 1,-80 },{ 1, 1,-80 },
                     {-2, 2,-80 },{ 2, 2,-80 },
                     {-1, 3,-90 },{ 1, 3,-90 },
        /* stern */  {-1,-3,-70 },{ 1,-3,-70 }
      },
      200.0, 700.0, 0.0001 },

    /* ── aegis defense sub (13 px) — 53 m × 24 m ────────────────────────── */
    { "aegis_sub", 13, {
        /* hull */   { 0,-3,-60 },{ 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },
                     { 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-60 },
        /* sensor arms */
                     {-3, 0,-70 },{ 3, 0,-70 },
                     {-2, 0,-80 },{ 2, 0,-80 },
        /* radar */  {-1,-1,-100 },{ 1,-1,-100 }
      },
      150.0, 650.0, 0.0002 },

    /* ── oblivion deep-sea dreadnought (17 px) — 115 m × 24 m ───────────── */
    { "oblivion_sub", 17, {
        /* hull */   { 0,-6,-50 },{ 0,-5,-80 },{ 0,-4,-80 },{ 0,-3,-80 },
                     { 0,-2,-80 },{ 0,-1,-80 },{ 0, 0,-80 },{ 0, 1,-80 },
                     { 0, 2,-80 },{ 0, 3,-80 },{ 0, 4,-80 },{ 0, 5,-80 },
                     { 0, 6,-50 },
        /* superstructure */
                     {-1, 0,-100 },{ 1, 0,-100 },
        /* bridge */ {-1,-1,-110 },{ 1,-1,-110 }
      },
      400.0, 900.0, 0.00005 },

    /* ── orca / killer whale (11 px) — 67 m × 20 m ─────────────────────── */
    { "orca", 11, {
        /* body */       { 0,-4,-70 },{ 0,-3,-100},{ 0,-2,-100},{ 0,-1,-100},
                         { 0, 0,-100},{ 0, 1,-100},{ 0, 2,-100},{ 0, 3,-70 },
        /* dorsal fin */ {-1, 1,-130},{-2, 1,-120},
        /* tail */       { 0, 4,-80 }
      },
      200.0, 800.0, 0.001 },

    /* ── sperm whale / cachalot (12 px) — 95 m × 20 m ──────────────────── */
    { "sperm_whale", 12, {
        /* block head */ {-1,-4,-80 },{ 0,-4,-90 },{ 1,-4,-80 },
                         {-1,-3,-80 },{ 0,-3,-90 },{ 1,-3,-80 },
        /* body */       { 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },
                         { 0, 1,-90 },{ 0, 2,-80 },
        /* tail */       { 0, 3,-70 }
      },
      400.0, 900.0, 0.0005 },

    /* ── shark (10 px) — 62 m × 24 m ───────────────────────────────────── */
    { "shark", 10, {
        /* body */       { 0,-3,-60 },{ 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },
                         { 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-60 },
        /* dorsal fin */ {-1, 0,-110},
        /* pec fins */   {-1,-1,-70 },{ 1,-1,-70 }
      },
      100.0, 800.0, 0.002 },

    /* ── barrel cluster (9 px) — 7 m × 12 m ────────────────────────────── */
    { "barrel_cluster", 9, {
        /* central barrel */  { 0, 0,-70 },
        /* ring of 4 */       {-2, 0,-70 },{ 2, 0,-70 },
                              { 0,-2,-70 },{ 0, 2,-70 },
        /* diagonal cluster */{-1,-1,-70 },{-1, 1,-70 },
                              { 1,-1,-70 },{ 1, 1,-70 }
      },
        0.0, 400.0, 0.002 },

    /* ── crude oil tanker / Aframax (17 px) — 250 m × 60 m ─────────────── */
    /* Surface vessel: bright hull reflection, dark cargo tank dividers.     */
    { "crude_tanker", 17, {
        /* hull bow→stern */
        { 0,-7,-60 },{ 0,-6,-90 },{ 0,-5,-90 },{ 0,-4,-90 },{ 0,-3,-90 },
        { 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },
        { 0, 3,-90 },{ 0, 4,-90 },{ 0, 5,-90 },{ 0, 6,-60 },
        /* tank bulkhead markings (darker cross-bands) */
        {-1,-2, 40 },{ 1,-2, 40 },
        /* bridge / accommodation */
        {-1, 6,-110 }
      },
       20.0, 450.0, 0.00015 },

    /* ── LNG carrier / Moss-type (18 px) — 290 m × 70 m ────────────────── */
    /* Wide hull, four prominent spherical cargo tank domes.                 */
    { "lng_carrier", 18, {
        /* hull */
        { 0,-7,-60 },{ 0,-6,-90 },{ 0,-5,-90 },{ 0,-4,-90 },{ 0,-3,-90 },
        { 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },
        { 0, 3,-90 },{ 0, 4,-90 },{ 0, 5,-90 },{ 0, 6,-60 },
        /* sphere tank dome highlights */
        {-1,-4,-120 },{-1,-1,-120 },{-1, 2,-120 },
        /* bridge */
        {-1, 6,-110 }
      },
       30.0, 500.0, 0.0001 },

    /* ── bulk carrier / Capesize (21 px) — 300 m × 64 m ────────────────── */
    /* Wide flat deck, five large hatch openings (darker), single bridge.   */
    { "bulk_carrier", 21, {
        /* hull */
        { 0,-7,-60 },{ 0,-6,-90 },{ 0,-5,-90 },{ 0,-4,-90 },{ 0,-3,-90 },
        { 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },
        { 0, 3,-90 },{ 0, 4,-90 },{ 0, 5,-90 },{ 0, 6,-60 },
        /* cargo hatch darkening (depth index raised = darker in hole) */
        {-1,-5, 40 },{-1,-3, 40 },{-1,-1, 40 },{-1, 1, 40 },{-1, 3, 40 },
        /* deck crane */
        { 1,-3,-80 },
        /* bridge */
        {-1, 6,-110 }
      },
       10.0, 450.0, 0.0002 },

    /* ── product tanker (15 px) — 180 m × 44 m ─────────────────────────── */
    /* Shorter hull, more numerous tank sub-divisions, pump room at stern.   */
    { "product_tanker", 15, {
        /* hull */
        { 0,-5,-60 },{ 0,-4,-90 },{ 0,-3,-90 },{ 0,-2,-90 },{ 0,-1,-90 },
        { 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-90 },{ 0, 4,-60 },
        /* tank dividers (darker bands at bulkheads) */
        {-1,-2, 40 },{-1, 0, 40 },{-1, 2, 40 },
        /* pump room highlight */
        { 1, 3,-70 },
        /* bridge */
        {-1, 4,-110 }
      },
       10.0, 350.0, 0.0002 },

    /* ── container ship / Neo-Panamax (20 px) — 300 m × 60 m ───────────── */
    /* Fine bow, wide stern; stacked container columns visible from above.   */
    { "container_ship", 20, {
        /* hull */
        { 0,-7,-60 },{ 0,-6,-90 },{ 0,-5,-90 },{ 0,-4,-90 },{ 0,-3,-90 },
        { 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },{ 0, 1,-90 },{ 0, 2,-90 },
        { 0, 3,-90 },{ 0, 4,-90 },{ 0, 5,-90 },{ 0, 6,-60 },
        /* container stack highlights (port side rows) */
        {-1,-5,-110 },{-1,-2,-110 },{-1, 1,-110 },
        /* container stacks (starboard) */
        { 1,-4,-90 },{ 1, 0,-90 },
        /* bridge (aft superstructure) */
        {-1, 6,-110 }
      },
       20.0, 500.0, 0.0003 },

    /* ── general cargo ship (14 px) — 140 m × 32 m ─────────────────────── */
    /* Compact multipurpose hull; cranes and cargo holds visible.           */
    { "cargo_ship", 14, {
        /* hull */
        { 0,-4,-60 },{ 0,-3,-90 },{ 0,-2,-90 },{ 0,-1,-90 },{ 0, 0,-90 },
        { 0, 1,-90 },{ 0, 2,-90 },{ 0, 3,-90 },{ 0, 4,-60 },
        /* mast */
        {-1, 0,-120 },
        /* deck cranes */
        {-1,-1,-80 },{-1, 1,-80 },
        /* cargo hold opening (darker) */
        { 1, 0,  40 },
        /* bridge */
        {-1, 4,-110 }
      },
        5.0, 400.0, 0.0004 }
};

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

PlacedDeco *place_decorations(Rbuf *work, const Rbuf *mask,
                               const struct Cell_head *region,
                               double pixel_m, double deco_scale,
                               double lo_clamp,
                               int *n_placed_out)
{
    (void)pixel_m;   /* reserved for future density scaling */
    (void)deco_scale; /* scale applied in vector_deco.c only */

    int  nrows  = work->nrows, ncols = work->ncols;
    long npx    = (long)nrows * ncols;

    /* collect ocean pixel indices for random sampling */
    long *ocean_idx = G_malloc(sizeof(long) * npx);
    long  n_ocean   = 0;
    for (long i = 0; i < npx; i++)
        if (!Rast_is_d_null_value(&mask->data[i]))
            ocean_idx[n_ocean++] = i;

    if (n_ocean == 0) {
        G_free(ocean_idx);
        *n_placed_out = 0;
        return NULL;
    }

    /* ocean area in km² */
    double pix_km2 = (region->ew_res * region->ns_res);
    /* check if geographic (degrees) and convert if needed */
    /* pixel_m is mean pixel size in metres; area = (pixel_m)^2 */
    double pix_m2 = pixel_m * pixel_m;
    pix_km2 = pix_m2 * 1e-6;
    double ocean_km2 = n_ocean * pix_km2;
    G_verbose_message(_("Ocean area for decorations: %.0f km²"), ocean_km2);

    /* allocate placement output array (upper bound) */
    int max_total = N_DECO_TYPES * MAX_PER_TYPE;
    PlacedDeco *places = G_malloc(sizeof(PlacedDeco) * (max_total + 1));
    int n_total = 0;

    /* seeded PRNG (reproducible results) */
    unsigned int seed = 0xDEADBEEFu;

    for (int t = 0; t < N_DECO_TYPES; t++) {
        const DecoType *d = &DECORATIONS[t];

        /* how many items proportional to ocean area */
        int count = (int)ceil(ocean_km2 * d->density);
        if (count < 1) count = 1;
        if (count > MAX_PER_TYPE) count = MAX_PER_TYPE;

        int placed = 0;
        int tries  = 0;
        int max_tries = count * 20;

        while (placed < count && tries < max_tries) {
            tries++;

            /* pick random ocean pixel */
            long pick = (long)((double)n_ocean *
                               (rand_r(&seed) / ((double)RAND_MAX + 1.0)));
            if (pick >= n_ocean) pick = n_ocean - 1;

            long idx = ocean_idx[pick];
            int  r   = (int)(idx / ncols);
            int  c   = (int)(idx % ncols);
            double cur = work->data[idx];

            if (cur < d->min_idx || cur > d->max_idx)
                continue;

            /* stamp raster pixel art */
            stamp(work, mask, r, c, d, lo_clamp);

            /* record placement for vector layer */
            if (n_total < max_total) {
                PlacedDeco *p = &places[n_total++];
                p->type_idx = t;
                p->east  = region->west  + (c + 0.5) * region->ew_res;
                p->north = region->north - (r + 0.5) * region->ns_res;
                /* random heading: draw from same PRNG for reproducibility */
                p->angle_rad = (rand_r(&seed) / ((double)RAND_MAX + 1.0))
                               * 2.0 * M_PI;
            }
            placed++;
        }

        G_verbose_message(
            _("  %s: %d placed (%.0f km², density %.4f/km²)"),
            d->name, placed, ocean_km2, d->density);
    }

    int total_placed = 0;
    for (int t = 0; t < N_DECO_TYPES; t++) {
        int c = 0;
        for (int i = 0; i < n_total; i++)
            if (places[i].type_idx == t) c++;
        total_placed += c;
    }
    G_message(_("Maritime decorations: %d items placed across %.0f km² of ocean."),
              total_placed, ocean_km2);

    G_free(ocean_idx);
    *n_placed_out = n_total;
    return places;
}
