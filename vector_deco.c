/****************************************************************************
 *
 * MODULE:       i.ocean  —  vector_deco.c
 * PURPOSE:      Create a GRASS GIS vector map of maritime decoration
 *               silhouettes at real-world scale.
 *
 * OVERVIEW
 * ────────
 * Each PlacedDeco (produced by place_decorations()) is rendered as a
 * set of GV_LINE features using 2-D silhouette shapes derived from:
 *
 *   Whale/dolphin — cetacean biology (blue whale avg 28 m, dolphin avg 3 m)
 *   Fish school   — typical shoal diameter 15 m
 *   Octopus       — common octopus arm span 3 m (bumped to dolphin minimum)
 *   Treasure/crate— bumped from real size (~1 m) to 3 m minimum
 *   Sunken ship   — ancient trireme hull, reference: 0AD
 *                   hele_ship_bireme.dae (beam 8.8 m, length 49.2 m);
 *                   scaled to canonical historical 37 m × 7 m
 *
 * COORDINATE SYSTEM
 * ─────────────────
 * All shapes are defined in LOCAL METRES:
 *   local-Y = along the decoration's heading (positive = bow/nose/forward)
 *   local-X = cross-track (positive = starboard / right)
 *
 * Transformation to map coordinates for a placement at (east₀, north₀)
 * with heading θ (clockwise from North):
 *
 *   Δeast  = (lx·cos θ + ly·sin θ) × m_to_ew
 *   Δnorth = (−lx·sin θ + ly·cos θ) × m_to_ns
 *
 * where  m_to_ew / m_to_ns  convert metres to map units:
 *   geographic (ll): m_to_ew = 1/(111320·cos(lat)), m_to_ns = 1/111320
 *   projected (m):   m_to_ew = m_to_ns = 1.0
 *   projected (ft):  m_to_ew = m_to_ns = 1/0.3048
 *
 * ATTRIBUTE TABLE
 * ───────────────
 * One row per decoration category (cat = type_idx + 1):
 *   cat  INTEGER  — GRASS feature category (= type index + 1)
 *   type TEXT     — decoration type name
 *   len_m REAL    — real-world length in metres (post-scale)
 *   wid_m REAL    — real-world width  in metres (post-scale)
 *   scale REAL    — deco_scale value applied
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/vector.h>
#include <grass/glocale.h>

#include "local_proto.h"

/* ── shape primitive: one GV_LINE feature in local metres ─────────────── */

typedef struct {
    int    n;          /* number of points                    */
    double pts[16][2]; /* (local_x, local_y) in metres        */
} Prim;

/* ── full decoration shape (set of primitives) ───────────────────────── */

typedef struct {
    int  n_prims;
    Prim prims[8];
} DecoShape;

/* ═══════════════════════════════════════════════════════════════════════
 * SHAPE DEFINITIONS  (local coordinates in metres; local-Y = heading)
 * ═══════════════════════════════════════════════════════════════════════ */

/* ── 0: BLUE WHALE (110 m × 30 m) ──────────────────────────────────────
 * Elongated teardrop hull + tail flukes + pectoral fins + skull/baleen.
 */
static const DecoShape shape_whale = { 6, {
    /* outer body (closed) */
    { 15, { {0,55},{6,45},{10,20},{11,0},{10,-20},{6,-40},{2,-50},
            {0,-55},{-2,-50},{-6,-40},{-10,-20},{-11,0},{-10,20},{-6,45},{0,55} } },
    /* left tail fluke */
    {  3, { {-2,-50},{-12,-55},{-8,-48} } },
    /* right tail fluke */
    {  3, { {2,-50},{12,-55},{8,-48} } },
    /* left pectoral fin */
    {  3, { {-6,5},{-15,0},{-8,-5} } },
    /* right pectoral fin */
    {  3, { {6,5},{15,0},{8,-5} } },
    /* skull / baleen zone */
    {  5, { {-4,50},{4,50},{3,35},{-3,35},{-4,50} } }
} };

/* ── 1: DOLPHIN (57 m × 16 m) ──────────────────────────────────────────
 * Sleek hull + dorsal fin + bifurcated tail + rostrum tip.
 */
static const DecoShape shape_dolphin = { 5, {
    /* outer body (closed) */
    { 13, { {0,28},{3,25},{5,15},{5,0},{4,-10},{2,-18},{0,-22},
            {-2,-18},{-4,-10},{-5,0},{-5,15},{-3,25},{0,28} } },
    /* dorsal fin */
    {  4, { {-1,2},{1,2},{0,10},{-1,2} } },
    /* left tail fluke */
    {  3, { {-1,-22},{-8,-25},{-4,-20} } },
    /* right tail fluke */
    {  3, { {1,-22},{8,-25},{4,-20} } },
    /* rostrum (snout tip triangle) */
    {  4, { {-1,28},{1,28},{0,32},{-1,28} } }
} };

/* ── 2: FISH SCHOOL (15 m × 15 m) ──────────────────────────────────────
 * 8 individual fish, each a small diamond shape (5-pt closed).
 * Centers spread over ±7.5 m from placement point.
 * Each fish is 0.6 m × 0.2 m.
 */
/* 12-position group layout from user spec (metres, pre-scaled by deco_scale) */
static double fish_centers[12][2] = {
    {-10,10},{-5,12},{ 0,10},{ 5,11},{10, 9},
    { -8, 5},{-3, 6},{ 2, 5},{ 7, 6},
    { -6, 0},{ 0, 1},{ 6, 0}
};
/* 8-pt fish body template scaled to ~0.8 m body (closed = 9 pts) */
static const double fish_tpl[9][2] = {
    { 0.00, 0.46},{ 0.17, 0.23},{ 0.17, 0.00},{ 0.11,-0.23},
    { 0.00,-0.34},{-0.11,-0.23},{-0.17, 0.00},{-0.17, 0.23},
    { 0.00, 0.46}  /* close */
};
#define N_FISH 12
#define N_FISH_TPL 9

/* Built dynamically in create_deco_vector() for this type */

/* ── 3: OCTOPUS (3 m × 3 m) ─────────────────────────────────────────────
 * Central body (8-pt circle, r=0.5 m) + 8 arms (2-pt lines to r=1.5 m)
 */
static const DecoShape shape_octopus_body = { 1, {
    { 9, { {0,0.5},{0.35,0.35},{0.5,0},{0.35,-0.35},
           {0,-0.5},{-0.35,-0.35},{-0.5,0},{-0.35,0.35},{0,0.5} } }
} };
/* arms built dynamically */
#define OCTO_BODY_R  0.5
#define OCTO_ARM_R   1.5
#define N_OCTO_ARMS  8

/* ── 4: TREASURE CHEST (8 m × 8 m) ──────────────────────────────────────
 * Rectangular body + curved lid + lock plate + internal treasure volume.
 */
static const DecoShape shape_treasure = { 4, {
    /* outer chest body */
    { 5, { {-4,-2},{4,-2},{4,2},{-4,2},{-4,-2} } },
    /* lid (curved top arc) */
    { 5, { {-4,2},{-2,5},{0,6},{2,5},{4,2} } },
    /* lock plate */
    { 5, { {-0.5,1},{0.5,1},{0.5,0},{-0.5,0},{-0.5,1} } },
    /* internal treasure volume */
    { 5, { {-3,1},{3,1},{3,-1},{-3,-1},{-3,1} } }
} };

/* ── 5: CRATE CLUSTER (6 m × 12 m) ──────────────────────────────────────
 * 8 individual 2×2 m crates stacked in a pyramid: 4 on bottom,
 * 3 on second layer, 1 on top. Each crate is a separate closed Prim.
 */
static const DecoShape shape_crate = { 8, {
    /* layer 1 — four crates centred at x=-4.5,-1.5,+1.5,+4.5, y=0 */
    { 5, { {-5.5,-1},{-3.5,-1},{-3.5,1},{-5.5,1},{-5.5,-1} } },
    { 5, { {-2.5,-1},{-0.5,-1},{-0.5,1},{-2.5,1},{-2.5,-1} } },
    { 5, { { 0.5,-1},{ 2.5,-1},{ 2.5,1},{ 0.5,1},{ 0.5,-1} } },
    { 5, { { 3.5,-1},{ 5.5,-1},{ 5.5,1},{ 3.5,1},{ 3.5,-1} } },
    /* layer 2 — three crates centred at x=-3,0,+3, y=2 */
    { 5, { {-4,1},{-2,1},{-2,3},{-4,3},{-4,1} } },
    { 5, { {-1,1},{ 1,1},{ 1,3},{-1,3},{-1,1} } },
    { 5, { { 2,1},{ 4,1},{ 4,3},{ 2,3},{ 2,1} } },
    /* layer 3 — one crate centred at x=0, y=4 */
    { 5, { {-1,3},{ 1,3},{ 1,5},{-1,5},{-1,3} } }
} };

/* ── 6: SUNKEN SHIP — trireme hull (37 m × 7 m) ─────────────────────────
 * Reference: 0AD hele_ship_bireme.dae bounding box 8.8 m × 49.2 m (beam × Y-length)
 * Scaled to historical canonical trireme (37 m hull, 7 m beam).
 * Features: hull rectangle, oar deck line, mast, ram at bow.
 */
static const DecoShape shape_ship = { 4, {
    /* hull rectangle (closed) */
    { 5, { {-3.5,-18.5},{3.5,-18.5},{3.5,18.5},{-3.5,18.5},{-3.5,-18.5} } },
    /* oar deck cross-beam */
    { 2, { {-3.5,-8},{3.5,-8} } },
    /* mast (amidships → up) */
    { 2, { {0,0},{0,8} } },
    /* bow ram (triangle) */
    { 4, { {0,18.5},{-1.5,20},{1.5,20},{0,18.5} } }
} };

/* ═══════════════════════════════════════════════════════════════════════
 * WAR VESSEL SHAPES  (indices 7–16)
 * Coordinate convention: origin at centre, +Y = bow, +X = starboard.
 * All dimensions in metres.
 * ═══════════════════════════════════════════════════════════════════════ */

/* ── 7: HUNTER-CLASS ATTACK SUBMARINE (78 m × 8 m) ──────────────────────
 * Streamlined spindle hull + conning tower + bow sonar dome.
 */
static const DecoShape shape_hunter = { 3, {
    /* outer hull */
    { 15, { {0,40},{2,35},{3.5,25},{4,10},{4,-10},{3,-25},{1.5,-35},
            {0,-38},{-1.5,-35},{-3,-25},{-4,-10},{-4,10},{-3.5,25},{-2,35},{0,40} } },
    /* conning tower */
    { 5, { {-2,8},{2,8},{2,2},{-2,2},{-2,8} } },
    /* sonar dome (bow triangle) */
    { 4, { {-1,38},{1,38},{0,36},{-1,38} } }
} };

/* ── 8: LEVIATHAN HEAVY BATTLESUB (95 m × 20 m) ─────────────────────────
 * Thick oval hull + dual missile silos + reactor compartment.
 */
static const DecoShape shape_leviathan = { 4, {
    /* outer hull */
    { 13, { {0,50},{6,45},{10,25},{10,0},{10,-25},{6,-40},{0,-45},
            {-6,-40},{-10,-25},{-10,0},{-10,25},{-6,45},{0,50} } },
    /* left missile silo */
    { 5, { {-6,12},{-4,12},{-4,8},{-6,8},{-6,12} } },
    /* right missile silo */
    { 5, { {4,12},{6,12},{6,8},{4,8},{4,12} } },
    /* reactor compartment */
    { 5, { {-5,-10},{5,-10},{5,-25},{-5,-25},{-5,-10} } }
} };

/* ── 9: BARRACUDA FAST STRIKE BOAT (50 m × 12 m) ────────────────────────
 * Dart-form hull with swept stabiliser fins at stern.
 */
static const DecoShape shape_barracuda = { 3, {
    /* outer hull */
    { 11, { {0,30},{2,20},{3,5},{2,-10},{1,-18},{0,-20},
            {-1,-18},{-2,-10},{-3,5},{-2,20},{0,30} } },
    /* left fin */
    { 4, { {-3,-5},{-6,-7},{-5,-10},{-2,-8} } },
    /* right fin */
    { 4, { {3,-5},{6,-7},{5,-10},{2,-8} } }
} };

/* ── 10: KRAKEN DRONE CARRIER SUB (65 m × 18 m) ─────────────────────────
 * Wide hull with internal hangar bay and drone launch rails.
 */
static const DecoShape shape_kraken = { 4, {
    /* outer hull */
    { 13, { {0,35},{5,30},{9,15},{9,0},{9,-20},{5,-28},{0,-30},
            {-5,-28},{-9,-20},{-9,0},{-9,15},{-5,30},{0,35} } },
    /* hangar bay */
    { 5, { {-8,15},{8,15},{8,-5},{-8,-5},{-8,15} } },
    /* left drone rail */
    { 5, { {-7,12},{-5,12},{-5,-2},{-7,-2},{-7,12} } },
    /* right drone rail */
    { 5, { {5,12},{7,12},{7,-2},{5,-2},{5,12} } }
} };

/* ── 11: PHANTOM STEALTH SUB (60 m × 10 m) ──────────────────────────────
 * Smooth blended hull — no external protrusions. Internal systems core.
 */
static const DecoShape shape_phantom = { 2, {
    /* outer hull */
    { 13, { {0,32},{3,28},{5,15},{5,0},{5,-15},{3,-25},{0,-28},
            {-3,-25},{-5,-15},{-5,0},{-5,15},{-3,28},{0,32} } },
    /* internal systems core */
    { 5, { {-3,10},{3,10},{3,-10},{-3,-10},{-3,10} } }
} };

/* ── 12: TITAN TORPEDO BARGE (47 m × 16 m) ──────────────────────────────
 * Blocky rectangular hull with flat bow + torpedo tube grid + ammo bay.
 */
static const DecoShape shape_titan = { 5, {
    /* outer hull (flat bow) */
    { 10, { {-8,25},{8,25},{8,0},{8,-15},{4,-20},{0,-22},
            {-4,-20},{-8,-15},{-8,0},{-8,25} } },
    /* torpedo tube row 1 */
    { 2, { {-6,20},{6,20} } },
    /* torpedo tube row 2 */
    { 2, { {-6,15},{6,15} } },
    /* torpedo tube row 3 */
    { 2, { {-6,10},{6,10} } },
    /* ammo storage */
    { 5, { {-5,5},{5,5},{5,-10},{-5,-10},{-5,5} } }
} };

/* ── 13: WRAITH MINI SUB (27 m × 5 m) ───────────────────────────────────
 * Compact oval hull with single diver hatch amidships.
 */
static const DecoShape shape_wraith = { 2, {
    /* outer hull */
    { 11, { {0,15},{2,12},{2.5,5},{2.5,0},{2,-8},{0,-12},
            {-2,-8},{-2.5,0},{-2.5,5},{-2,12},{0,15} } },
    /* hatch */
    { 5, { {-1,1},{1,1},{1,-1},{-1,-1},{-1,1} } }
} };

/* ── 14: HYDRA MULTI-HEAD SUB (55 m × 10 m) ─────────────────────────────
 * Central hull splitting into twin bow prongs, each with a weapon node.
 */
static const DecoShape shape_hydra = { 5, {
    /* central hull */
    { 9, { {0,10},{4,5},{5,-10},{4,-25},{0,-30},
           {-4,-25},{-5,-10},{-4,5},{0,10} } },
    /* left prong */
    { 4, { {0,10},{-3,15},{-5,25},{-2,22} } },
    /* right prong */
    { 4, { {0,10},{3,15},{5,25},{2,22} } },
    /* left weapon node (triangle) */
    { 4, { {-5,25},{-4,23},{-6,23},{-5,25} } },
    /* right weapon node (triangle) */
    { 4, { {5,25},{4,23},{6,23},{5,25} } }
} };

/* ── 15: AEGIS DEFENSE SUB (53 m × 24 m) ────────────────────────────────
 * Standard hull flanked by wide sensor arms + radar core amidships.
 */
static const DecoShape shape_aegis = { 4, {
    /* outer hull */
    { 13, { {0,28},{5,25},{7,10},{7,0},{7,-15},{3,-22},{0,-25},
            {-3,-22},{-7,-15},{-7,0},{-7,10},{-5,25},{0,28} } },
    /* left sensor arm */
    { 5, { {-7,5},{-12,5},{-12,3},{-7,3},{-7,5} } },
    /* right sensor arm */
    { 5, { {7,5},{12,5},{12,3},{7,3},{7,5} } },
    /* radar core */
    { 5, { {-3,8},{3,8},{3,0},{-3,0},{-3,8} } }
} };

/* ── 16: OBLIVION DEEP-SEA DREADNOUGHT (115 m × 24 m) ───────────────────
 * Massive elongated oval + superstructure + command bridge + reactor core.
 */
static const DecoShape shape_oblivion = { 4, {
    /* outer hull */
    { 13, { {0,60},{8,55},{12,30},{12,0},{12,-30},{8,-50},{0,-55},
            {-8,-50},{-12,-30},{-12,0},{-12,30},{-8,55},{0,60} } },
    /* superstructure */
    { 5, { {-6,20},{6,20},{6,-10},{-6,-10},{-6,20} } },
    /* command bridge */
    { 5, { {-3,25},{3,25},{3,15},{-3,15},{-3,25} } },
    /* reactor core */
    { 5, { {-4,-20},{4,-20},{4,-40},{-4,-40},{-4,-20} } }
} };

/* ── 17: ORCA / KILLER WHALE (67 m × 20 m) ──────────────────────────────
 * Robust spindle hull + tall dorsal fin + bifurcated tail + eye patches.
 */
static const DecoShape shape_orca = { 6, {
    /* outer body (closed) */
    { 15, { {0,35},{5,30},{8,15},{8,0},{7,-10},{5,-20},{2,-28},
            {0,-30},{-2,-28},{-5,-20},{-7,-10},{-8,0},{-8,15},{-5,30},{0,35} } },
    /* tall dorsal fin */
    {  4, { {-1,5},{1,5},{0,18},{-1,5} } },
    /* left tail fluke */
    {  3, { {-2,-28},{-10,-32},{-6,-26} } },
    /* right tail fluke */
    {  3, { {2,-28},{10,-32},{6,-26} } },
    /* left eye patch */
    {  4, { {-4,20},{-2,18},{-3,16},{-4,20} } },
    /* right eye patch */
    {  4, { {4,20},{2,18},{3,16},{4,20} } }
} };

/* ── 18: SPERM WHALE / CACHALOT (95 m × 20 m) ───────────────────────────
 * Block head + tapering body + spermaceti organ + tail flukes + jawline.
 */
static const DecoShape shape_sperm_whale = { 5, {
    /* outer body (closed; block head clearly visible) */
    { 14, { {-6,45},{6,45},{10,20},{10,0},{8,-20},{5,-35},{2,-45},
            {0,-48},{-2,-45},{-5,-35},{-8,-20},{-10,0},{-10,20},{-6,45} } },
    /* spermaceti head organ */
    {  5, { {-5,45},{5,45},{4,25},{-4,25},{-5,45} } },
    /* left tail fluke */
    {  3, { {-2,-45},{-10,-50},{-6,-42} } },
    /* right tail fluke */
    {  3, { {2,-45},{10,-50},{6,-42} } },
    /* lower jawline */
    {  5, { {-3,20},{3,20},{2,0},{-2,0},{-3,20} } }
} };

/* ── 19: SHARK (62 m × 24 m) ────────────────────────────────────────────
 * Fusiform body + dorsal fin + pectoral fins + heterocercal tail lobes.
 */
static const DecoShape shape_shark = { 6, {
    /* outer body (closed) */
    { 13, { {0,30},{4,25},{7,10},{7,0},{6,-10},{3,-20},{0,-25},
            {-3,-20},{-6,-10},{-7,0},{-7,10},{-4,25},{0,30} } },
    /* dorsal fin */
    {  4, { {-1,5},{1,5},{0,15},{-1,5} } },
    /* left pectoral fin */
    {  3, { {-5,0},{-12,-3},{-6,-6} } },
    /* right pectoral fin */
    {  3, { {5,0},{12,-3},{6,-6} } },
    /* upper tail lobe (longer) */
    {  3, { {0,-25},{6,-32},{2,-20} } },
    /* lower tail lobe */
    {  3, { {0,-25},{-4,-30},{-1,-20} } }
} };

/* ── 20: BARREL CLUSTER (7 m × 12 m) ────────────────────────────────────
 * Five octagonal barrels (r=2 m each) in a triangular cluster.
 * Centers: (-4,0),(0,0),(+4,0) bottom row; (-2,3),(+2,3) upper row.
 */
static const DecoShape shape_barrels = { 5, {
    /* barrel at (-4, 0) */
    { 9, { {-4,2},{-2.5,1.5},{-2,0},{-2.5,-1.5},{-4,-2},{-5.5,-1.5},{-6,0},{-5.5,1.5},{-4,2} } },
    /* barrel at (0, 0) */
    { 9, { {0,2},{1.5,1.5},{2,0},{1.5,-1.5},{0,-2},{-1.5,-1.5},{-2,0},{-1.5,1.5},{0,2} } },
    /* barrel at (+4, 0) */
    { 9, { {4,2},{5.5,1.5},{6,0},{5.5,-1.5},{4,-2},{2.5,-1.5},{2,0},{2.5,1.5},{4,2} } },
    /* barrel at (-2, 3) */
    { 9, { {-2,5},{-0.5,4.5},{0,3},{-0.5,1.5},{-2,1},{-3.5,1.5},{-4,3},{-3.5,4.5},{-2,5} } },
    /* barrel at (+2, 3) */
    { 9, { {2,5},{3.5,4.5},{4,3},{3.5,1.5},{2,1},{0.5,1.5},{0,3},{0.5,4.5},{2,5} } }
} };

/* ═══════════════════════════════════════════════════════════════════════
 * MERCHANT VESSEL SHAPES  (indices 21–26)
 * Hull vertices follow user spec with minor trimming to fit pts[16]:
 *   original specs had up to 16 unique hull vertices; closing repeat
 *   brings total to 17 — one bow/stern vertex pair is merged to 15+close.
 * All internal structures (tanks, holds, bays) are closed 4-corner rects.
 * Dimensions in metres; coordinate convention: +Y = bow, +X = starboard.
 * ═══════════════════════════════════════════════════════════════════════ */

/* ── 21: CRUDE OIL TANKER — Aframax class (250 m × 60 m) ────────────────
 * Full mid-body block coefficient ~0.85; four segmented cargo tanks;
 * pump room and bridge/accommodation at stern.
 */
static const DecoShape shape_crude_tanker = { 7, {
    /* outer hull (14 unique + close = 15) */
    { 15, { {0,125},{15,115},{25,80},{30,40},{30,0},{28,-60},{20,-110},
            {0,-125},{-20,-110},{-28,-60},{-30,0},{-30,40},{-25,80},{-15,115},{0,125} } },
    /* cargo tank 1  (fwd) */
    {  5, { {-18,80},{18,80},{18,40},{-18,40},{-18,80} } },
    /* cargo tank 2 */
    {  5, { {-18,40},{18,40},{18,0},{-18,0},{-18,40} } },
    /* cargo tank 3 */
    {  5, { {-18,0},{18,0},{18,-40},{-18,-40},{-18,0} } },
    /* cargo tank 4  (aft) */
    {  5, { {-18,-40},{18,-40},{18,-80},{-18,-80},{-18,-40} } },
    /* pump room */
    {  5, { {-8,-90},{8,-90},{8,-110},{-8,-110},{-8,-90} } },
    /* bridge / accommodation block */
    {  5, { {-10,-110},{10,-110},{10,-90},{-10,-90},{-10,-110} } }
} };

/* ── 22: LNG CARRIER — Moss-type (290 m × 70 m) ─────────────────────────
 * Wide hull carrying four spherical insulated cargo tanks; engine room
 * and bridge in the aft third.
 */
static const DecoShape shape_lng_carrier = { 7, {
    /* outer hull (14 unique + close = 15) */
    { 15, { {0,145},{18,135},{30,100},{35,50},{35,0},{32,-70},{22,-130},
            {0,-145},{-22,-130},{-32,-70},{-35,0},{-35,50},{-30,100},{-18,135},{0,145} } },
    /* spherical tank A (fwd) footprint */
    {  5, { {-12,80},{12,80},{12,50},{-12,50},{-12,80} } },
    /* spherical tank B */
    {  5, { {-12,40},{12,40},{12,10},{-12,10},{-12,40} } },
    /* spherical tank C */
    {  5, { {-12,0},{12,0},{12,-30},{-12,-30},{-12,0} } },
    /* spherical tank D (aft) */
    {  5, { {-12,-40},{12,-40},{12,-70},{-12,-70},{-12,-40} } },
    /* engine room */
    {  5, { {-10,-120},{10,-120},{10,-140},{-10,-140},{-10,-120} } },
    /* bridge */
    {  5, { {-12,-140},{12,-140},{12,-115},{-12,-115},{-12,-140} } }
} };

/* ── 23: BULK CARRIER — Capesize (300 m × 64 m) ─────────────────────────
 * Flat-deck boxy hull; five large rectangular cargo hold hatches dominate
 * the deck; single bridge/accommodation superstructure at stern.
 */
static const DecoShape shape_bulk_carrier = { 7, {
    /* outer hull — trimmed to 14 unique + close = 15
     * (original (32,0) and (-32,0) mid-points merged into neighbouring pts) */
    { 15, { {0,150},{18,140},{28,100},{32,50},{32,-50},{28,-110},{18,-140},
            {0,-150},{-18,-140},{-28,-110},{-32,-50},{-32,50},{-28,100},{-18,140},{0,150} } },
    /* cargo hold 1 */
    {  5, { {-20,100},{20,100},{20,60},{-20,60},{-20,100} } },
    /* cargo hold 2 */
    {  5, { {-20,60},{20,60},{20,20},{-20,20},{-20,60} } },
    /* cargo hold 3 */
    {  5, { {-20,20},{20,20},{20,-20},{-20,-20},{-20,20} } },
    /* cargo hold 4 */
    {  5, { {-20,-20},{20,-20},{20,-60},{-20,-60},{-20,-20} } },
    /* cargo hold 5 */
    {  5, { {-20,-60},{20,-60},{20,-100},{-20,-100},{-20,-60} } },
    /* bridge / accommodation */
    {  5, { {-12,-140},{12,-140},{12,-110},{-12,-110},{-12,-140} } }
} };

/* ── 24: PRODUCT TANKER (~180 m × 44 m) ─────────────────────────────────
 * Refined-fuel tanker; more numerous but smaller cargo tanks than crude
 * tanker; dedicated pump room aft; fine bow entry.
 */
static const DecoShape shape_product_tanker = { 7, {
    /* outer hull (14 unique + close = 15) */
    { 15, { {0,90},{10,80},{18,50},{22,20},{22,0},{20,-40},{14,-80},
            {0,-90},{-14,-80},{-20,-40},{-22,0},{-22,20},{-18,50},{-10,80},{0,90} } },
    /* product tank 1 (fwd) */
    {  5, { {-14,60},{14,60},{14,30},{-14,30},{-14,60} } },
    /* product tank 2 */
    {  5, { {-14,30},{14,30},{14,0},{-14,0},{-14,30} } },
    /* product tank 3 */
    {  5, { {-14,0},{14,0},{14,-30},{-14,-30},{-14,0} } },
    /* product tank 4 (aft) */
    {  5, { {-14,-30},{14,-30},{14,-60},{-14,-60},{-14,-30} } },
    /* pump room */
    {  5, { {-6,-70},{6,-70},{6,-85},{-6,-85},{-6,-70} } },
    /* bridge */
    {  5, { {-8,-85},{8,-85},{8,-65},{-8,-65},{-8,-85} } }
} };

/* ── 25: CONTAINER SHIP — Neo-Panamax (300 m × 60 m) ────────────────────
 * Fine bow entry, wide stern; six longitudinal container bay grids fill
 * the main deck; single aft bridge block.
 * Uses all 8 primitive slots: 1 hull + 6 bays + 1 bridge.
 */
static const DecoShape shape_container_ship = { 8, {
    /* outer hull (14 unique + close = 15) */
    { 15, { {0,150},{12,140},{22,100},{28,50},{30,0},{30,-60},{25,-120},
            {0,-150},{-25,-120},{-30,-60},{-30,0},{-28,50},{-22,100},{-12,140},{0,150} } },
    /* container bay 1 (fwd) */
    {  5, { {-20,110},{20,110},{20,80},{-20,80},{-20,110} } },
    /* container bay 2 */
    {  5, { {-20,80},{20,80},{20,50},{-20,50},{-20,80} } },
    /* container bay 3 */
    {  5, { {-20,50},{20,50},{20,20},{-20,20},{-20,50} } },
    /* container bay 4 */
    {  5, { {-20,20},{20,20},{20,-10},{-20,-10},{-20,20} } },
    /* container bay 5 */
    {  5, { {-20,-10},{20,-10},{20,-40},{-20,-40},{-20,-10} } },
    /* container bay 6 (aft) */
    {  5, { {-20,-40},{20,-40},{20,-70},{-20,-70},{-20,-40} } },
    /* bridge / aft superstructure */
    {  5, { {-12,-140},{12,-140},{12,-110},{-12,-110},{-12,-140} } }
} };

/* ── 26: GENERAL CARGO SHIP (~140 m × 32 m) ─────────────────────────────
 * Compact multipurpose hull; three cargo holds amidships;
 * twin deck cranes (open triangles); bridge at stern.
 */
static const DecoShape shape_cargo_ship = { 6, {
    /* outer hull (14 unique + close = 15) */
    { 15, { {0,70},{8,60},{14,30},{16,0},{15,-30},{10,-60},{4,-68},
            {0,-70},{-4,-68},{-10,-60},{-15,-30},{-16,0},{-14,30},{-8,60},{0,70} } },
    /* cargo hold 1 (fwd) */
    {  5, { {-10,40},{10,40},{10,10},{-10,10},{-10,40} } },
    /* cargo hold 2 (mid) */
    {  5, { {-10,10},{10,10},{10,-20},{-10,-20},{-10,10} } },
    /* cargo hold 3 (aft) */
    {  5, { {-10,-20},{10,-20},{10,-50},{-10,-50},{-10,-20} } },
    /* deck crane port (open triangle, not closed) */
    {  3, { {-12,25},{-4,5},{0,25} } },
    /* bridge */
    {  5, { {-8,-65},{8,-65},{8,-45},{-8,-45},{-8,-65} } }
} };

/* ── shape table (indexed by type_idx, matching DECORATIONS[] order) ─── */
/* NOTE: fish_school (idx 2) and octopus arms (idx 3) are built below.    */
static const DecoShape *shape_table[N_DECO_TYPES] = {
    &shape_whale,           /* 0  blue whale         */
    &shape_dolphin,         /* 1  dolphin            */
    NULL,                   /* 2  fish_school — generated dynamically */
    NULL,                   /* 3  octopus — body + dynamic arms */
    &shape_treasure,        /* 4  treasure chest     */
    &shape_crate,           /* 5  crate cluster      */
    &shape_ship,            /* 6  sunken ship        */
    &shape_hunter,          /* 7  hunter sub         */
    &shape_leviathan,       /* 8  leviathan sub      */
    &shape_barracuda,       /* 9  barracuda boat     */
    &shape_kraken,          /* 10 kraken carrier     */
    &shape_phantom,         /* 11 phantom sub        */
    &shape_titan,           /* 12 titan barge        */
    &shape_wraith,          /* 13 wraith sub         */
    &shape_hydra,           /* 14 hydra sub          */
    &shape_aegis,           /* 15 aegis sub          */
    &shape_oblivion,        /* 16 oblivion sub       */
    &shape_orca,            /* 17 orca               */
    &shape_sperm_whale,     /* 18 sperm whale        */
    &shape_shark,           /* 19 shark              */
    &shape_barrels,         /* 20 barrel cluster     */
    &shape_crude_tanker,    /* 21 crude oil tanker   */
    &shape_lng_carrier,     /* 22 LNG carrier        */
    &shape_bulk_carrier,    /* 23 bulk carrier       */
    &shape_product_tanker,  /* 24 product tanker     */
    &shape_container_ship,  /* 25 container ship     */
    &shape_cargo_ship       /* 26 general cargo ship */
};

/* ── coordinate conversion helpers ──────────────────────────────────────
 *
 * For geographic CRS (proj == "ll"):
 *   m_to_ew = 1 / (111320 * cos(lat_rad))    degrees per metre (EW)
 *   m_to_ns = 1 / 111320                      degrees per metre (NS)
 * For projected CRS in metres:
 *   m_to_ew = m_to_ns = 1.0
 * For projected CRS in feet:
 *   m_to_ew = m_to_ns = 1 / 0.3048
 */
static void crs_scale(const char *proj, const char *units,
                      double latitude,
                      double *m_to_ew, double *m_to_ns)
{
    int geographic =
        (strncmp(proj, "ll", 2) == 0) ||
        (strstr(units, "degree") != NULL);

    if (geographic) {
        double lat_rad = latitude * M_PI / 180.0;
        *m_to_ew = 1.0 / (111320.0 * cos(lat_rad));
        *m_to_ns = 1.0 / 111320.0;
    } else if (strstr(units, "feet") || strstr(units, "foot")) {
        *m_to_ew = 1.0 / 0.3048;
        *m_to_ns = 1.0 / 0.3048;
    } else {
        /* metres (default) */
        *m_to_ew = 1.0;
        *m_to_ns = 1.0;
    }
}

/* ── write one feature: GV_BOUNDARY+GV_CENTROID if closed, GV_LINE if open ─
 *
 * Closed polygon detection: n_pts >= 4 AND first point == last point
 * (coordinate equality — shapes are defined with exact repeated vertices).
 *
 * Closed path → GV_BOUNDARY (no category) + GV_CENTROID at vertex-average
 *               centroid (with category).  This registers a topological area.
 * Open path   → GV_LINE with category (unchanged behaviour).
 */

static void write_prim(struct Map_info *Map,
                        struct line_pnts *Points,
                        struct line_cats *Cats,
                        struct line_cats *Cats_empty,
                        const double (*pts_m)[2], int n_pts,
                        double cx, double cy,
                        double cos_t, double sin_t,
                        double scale, double m_ew, double m_ns,
                        int cat)
{
    if (n_pts < 2) return;

    /* closed if first == last (exact match — shapes use repeated vertex) */
    int is_closed = (n_pts >= 4 &&
                     pts_m[0][0] == pts_m[n_pts - 1][0] &&
                     pts_m[0][1] == pts_m[n_pts - 1][1]);

    /* transform all points into map coordinates */
    Vect_reset_line(Points);
    for (int i = 0; i < n_pts; i++) {
        double lx = pts_m[i][0] * scale;
        double ly = pts_m[i][1] * scale;
        double de =  lx * cos_t + ly * sin_t;
        double dn = -lx * sin_t + ly * cos_t;
        Vect_append_point(Points, cx + de * m_ew, cy + dn * m_ns, 0.0);
    }

    if (is_closed) {
        /* ── boundary (no category) ───────────────────────────────────── */
        Vect_reset_cats(Cats_empty);
        Vect_write_line(Map, GV_BOUNDARY, Points, Cats_empty);

        /* ── centroid at vertex average (skip repeated closing point) ─── */
        int n_inner = n_pts - 1;
        double sum_x = 0.0, sum_y = 0.0;
        for (int i = 0; i < n_inner; i++) {
            sum_x += Points->x[i];
            sum_y += Points->y[i];
        }
        Vect_reset_line(Points);
        Vect_append_point(Points,
                          sum_x / n_inner, sum_y / n_inner, 0.0);
        Vect_reset_cats(Cats);
        Vect_cat_set(Cats, 1, cat);
        Vect_write_line(Map, GV_CENTROID, Points, Cats);
    } else {
        /* ── open line with category ──────────────────────────────────── */
        Vect_reset_cats(Cats);
        Vect_cat_set(Cats, 1, cat);
        Vect_write_line(Map, GV_LINE, Points, Cats);
    }
}

/* ── main entry point ────────────────────────────────────────────────── */

void create_deco_vector(const char *vmap_name,
                        const PlacedDeco *places, int n_places,
                        const struct Cell_head *region,
                        const char *proj, const char *units,
                        double latitude, double deco_scale)
{
    (void)region;   /* region reserved for future per-placement CRS checks */

    if (n_places == 0) {
        G_verbose_message(_("No decoration placements; vector map not created."));
        return;
    }

    /* ── open new vector map ─────────────────────────────────────────── */
    struct Map_info Map;
    if (Vect_open_new(&Map, vmap_name, 0) < 0) {
        G_warning(_("Could not create vector map <%s>"), vmap_name);
        return;
    }

    /* set basic header info */
    Vect_set_organization(&Map, "i.ocean");
    {
        char hist[256];
        snprintf(hist, sizeof(hist),
                 "Created by i.ocean with deco_scale=%.3f", deco_scale);
        Vect_set_scale(&Map, 1);
        Vect_set_map_name(&Map, "Maritime decorations");
        Vect_hist_command(&Map);
    }

    struct line_pnts *Points     = Vect_new_line_struct();
    struct line_cats *Cats       = Vect_new_cats_struct();
    struct line_cats *Cats_empty = Vect_new_cats_struct(); /* for GV_BOUNDARY */

    /* ── CRS conversion factors ──────────────────────────────────────── */
    double m_ew, m_ns;
    crs_scale(proj, units, latitude, &m_ew, &m_ns);

    /* ── draw each placement ─────────────────────────────────────────── */
    for (int p = 0; p < n_places; p++) {
        const PlacedDeco *pd = &places[p];
        int  t    = pd->type_idx;
        int  cat  = t + 1;           /* 1-based GRASS category */

        double cx   = pd->east;
        double cy   = pd->north;
        double theta = pd->angle_rad; /* clockwise from North */
        double cos_t = cos(theta);
        double sin_t = sin(theta);
        double sc   = deco_scale;

        if (t == 2) {
            /* ── fish school: 12 fish, 8-pt body template ────────────── */
            for (int f = 0; f < N_FISH; f++) {
                double offs_x = fish_centers[f][0] * sc;
                double offs_y = fish_centers[f][1] * sc;
                double fish_pts[N_FISH_TPL][2];
                for (int j = 0; j < N_FISH_TPL; j++) {
                    fish_pts[j][0] = offs_x + fish_tpl[j][0] * sc;
                    fish_pts[j][1] = offs_y + fish_tpl[j][1] * sc;
                }
                write_prim(&Map, Points, Cats, Cats_empty,
                           (const double(*)[2])fish_pts, N_FISH_TPL,
                           cx, cy, cos_t, sin_t, 1.0, m_ew, m_ns, cat);
            }
        } else if (t == 3) {
            /* ── octopus: body + 8 arms ────────────────────────────────── */
            /* body polygon */
            for (int i = 0; i < shape_octopus_body.n_prims; i++) {
                const Prim *pr = &shape_octopus_body.prims[i];
                write_prim(&Map, Points, Cats, Cats_empty,
                           (const double(*)[2])pr->pts, pr->n,
                           cx, cy, cos_t, sin_t, sc, m_ew, m_ns, cat);
            }
            /* 8 arms */
            for (int a = 0; a < N_OCTO_ARMS; a++) {
                double alpha = a * (2.0 * M_PI / N_OCTO_ARMS);
                double ca = cos(alpha), sa = sin(alpha);
                double arm_pts[2][2] = {
                    { OCTO_BODY_R * ca,  OCTO_BODY_R * sa },
                    { OCTO_ARM_R  * ca,  OCTO_ARM_R  * sa }
                };
                write_prim(&Map, Points, Cats, Cats_empty,
                           (const double(*)[2])arm_pts, 2,
                           cx, cy, cos_t, sin_t, sc, m_ew, m_ns, cat);
            }
        } else {
            /* ── all other types: iterate shape primitives ──────────── */
            const DecoShape *sh = shape_table[t];
            if (!sh) continue;
            for (int i = 0; i < sh->n_prims; i++) {
                const Prim *pr = &sh->prims[i];
                write_prim(&Map, Points, Cats, Cats_empty,
                           (const double(*)[2])pr->pts, pr->n,
                           cx, cy, cos_t, sin_t, sc, m_ew, m_ns, cat);
            }
        }
    }

    /* ── build topology and close ────────────────────────────────────── */
    Vect_build(&Map);
    Vect_close(&Map);

    Vect_destroy_line_struct(Points);
    Vect_destroy_cats_struct(Cats);
    Vect_destroy_cats_struct(Cats_empty);

    G_message(_("Vector decoration map <%s> created (%d features)."),
              vmap_name, n_places);

    /* ── write per-type colour table ─────────────────────────────────── */
    /*
     * Realistic ocean palette:
     *   blue/grey tones for cetaceans; warm browns for treasure/crates;
     *   dark tactical greens/greys for military subs; black-white for orca;
     *   amber for barrels.
     */
    {
        static const int deco_rgb[N_DECO_TYPES][3] = {
            { 70,100,160},  /*  1 blue whale     — deep ocean blue  */
            {140,170,190},  /*  2 dolphin         — slate blue-grey */
            {190,210,220},  /*  3 fish school     — pale ice blue   */
            { 90, 50, 30},  /*  4 octopus         — rust brown      */
            {180,145, 15},  /*  5 treasure chest  — aged gold       */
            {110, 75, 35},  /*  6 crate cluster   — weathered wood  */
            { 55, 50, 45},  /*  7 sunken ship     — dark rust       */
            { 40, 55, 40},  /*  8 hunter sub      — military green  */
            { 25, 30, 45},  /*  9 leviathan sub   — deep navy       */
            { 45, 55, 55},  /* 10 barracuda boat  — teal-grey       */
            { 35, 45, 50},  /* 11 kraken carrier  — slate grey      */
            { 30, 35, 30},  /* 12 phantom sub     — shadow green    */
            { 60, 50, 35},  /* 13 titan barge     — olive-brown     */
            { 38, 45, 38},  /* 14 wraith sub      — dark sage       */
            { 35, 45, 60},  /* 15 hydra sub       — blue-grey       */
            { 45, 55, 70},  /* 16 aegis sub       — steel blue      */
            { 20, 20, 25},  /* 17 oblivion sub    — near-black      */
            {  8,  8,  8},  /* 18 orca            — black (eye patch white added via GV_LINE) */
            { 90, 78, 65},  /* 19 sperm whale     — sandy tan       */
            {100,115,125},  /* 20 shark           — blue-grey       */
            {110, 60, 30},  /* 21 barrel cluster  — rust orange     */
            /* merchant vessels */
            {120, 55, 25},  /* 22 crude tanker    — weathered rust  */
            {210,215,220},  /* 23 LNG carrier     — polished silver */
            { 45, 48, 42},  /* 24 bulk carrier    — dark charcoal   */
            { 95, 35, 35},  /* 25 product tanker  — deep maroon     */
            { 45,115,165},  /* 26 container ship  — maritime blue   */
            { 65, 78, 58}   /* 27 cargo ship      — olive grey      */
        };

        struct Colors colors;
        Rast_init_colors(&colors);
        for (int ti = 0; ti < N_DECO_TYPES; ti++) {
            CELL c = (CELL)(ti + 1);
            Rast_set_c_color(c,
                             deco_rgb[ti][0],
                             deco_rgb[ti][1],
                             deco_rgb[ti][2],
                             &colors);
        }
        Vect_write_colors(vmap_name, G_mapset(), &colors);
        G_verbose_message(_("Colour table written for vector map <%s>."), vmap_name);
        Rast_free_colors(&colors);
    }

    /* ── add attribute table via v.db.addtable ───────────────────────── */
    {
        char cmd[1024];

        /* create table with one row per type */
        snprintf(cmd, sizeof(cmd),
                 "v.db.addtable map=%s columns='type TEXT,len_m REAL,"
                 "wid_m REAL,scale REAL' --quiet 2>/dev/null",
                 vmap_name);
        if (system(cmd) != 0) {
            G_warning(_("v.db.addtable failed for <%s>"), vmap_name);
            return;
        }

        for (int t = 0; t < N_DECO_TYPES; t++) {
            double len = deco_real_length[t] * deco_scale;
            double wid = deco_real_width[t]  * deco_scale;

            /* insert row for this category */
            snprintf(cmd, sizeof(cmd),
                     "db.execute sql=\"INSERT INTO %s "
                     "(cat,type,len_m,wid_m,scale) "
                     "VALUES (%d,'%s',%.2f,%.2f,%.4f)\" --quiet 2>/dev/null",
                     vmap_name, t + 1,
                     deco_type_name[t], len, wid, deco_scale);
            if (system(cmd) != 0)
                G_verbose_message(
                    _("db.execute: could not insert record for cat %d."), t+1);
        }
    }
}
