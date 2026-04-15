/****************************************************************************
 *
 * MODULE:       i.ocean
 * AUTHOR(S):    i.ocean contributors
 * PURPOSE:      Renders ocean areas of a raster map with a realistic visual
 *               appearance, combining depth-based colour gradients,
 *               latitude-driven colour temperature, optional wave textures,
 *               FFTW3-based foam/turbulence, and maritime decorations.
 *
 *               Key implementation notes
 *               ─────────────────────────
 *               • All raster work is done on in-memory DCELL buffers (Rbuf).
 *               • OpenMP parallelism via omp_get_max_threads(); forwarded to
 *                 every row-loop that supports it.
 *               • FFTW3 (turbulence.c) generates the Phillips-spectrum ocean
 *                 surface for the -f flag.
 *               • Integer (CELL) depth maps are smoothed with a 3×3 average
 *                 before normalisation (depth_index.c:smooth_3x3).
 *               • Progress is reported with G_percent() at each named step.
 *
 * COPYRIGHT:    (C) 2026 by the GRASS Development Team
 *               This program is free software under the GNU General Public
 *               License (>=v2). Read the file COPYING that comes with GRASS
 *               for details.
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <unistd.h>

#include <omp.h>

#include <grass/gis.h>
#include <grass/raster.h>
#include <grass/glocale.h>

#include "local_proto.h"

/* ── palette data (depth_index 0-1000, RGB) ──────────────────────────── */

#define N_STOPS 7

static PaletteStop pal_tropical[N_STOPS] = {
    {0,    0,230,255}, {60,   0,200,225}, {120,  0,165,185},
    {220,  0,125,165}, {420,  0, 82,144}, {700,  0, 50,100},
    {1000, 0, 20, 58}
};
static PaletteStop pal_subtropical[N_STOPS] = {
    {0,   25,210,230}, {60,  10,180,210}, {120,  0,145,178},
    {220,  0,105,158}, {420,  0, 68,132}, {700,  0, 42, 92},
    {1000, 0, 16, 56}
};
static PaletteStop pal_temperate[N_STOPS] = {
    {0,   65,190,218}, {60,  38,160,198}, {120, 18,126,172},
    {220,  5, 92,152}, {420,  0, 60,122}, {700,  0, 36, 86},
    {1000, 0, 14, 50}
};
static PaletteStop pal_subpolar[N_STOPS] = {
    {0,  125,195,218}, {60,  88,165,198}, {120, 58,136,177},
    {220, 32,102,152}, {420, 14, 70,122}, {700,  5, 44, 86},
    {1000, 0, 20, 52}
};
static PaletteStop pal_polar[N_STOPS] = {
    {0,  178,215,232}, {60, 140,186,212}, {120, 98,156,192},
    {220, 62,120,166}, {420, 32, 84,136}, {700, 14, 53, 96},
    {1000, 4, 22, 60}
};

static const PaletteStop *get_palette(const char *style)
{
    if (strcmp(style, "tropical")    == 0) return pal_tropical;
    if (strcmp(style, "subtropical") == 0) return pal_subtropical;
    if (strcmp(style, "subpolar")    == 0) return pal_subpolar;
    if (strcmp(style, "polar")       == 0) return pal_polar;
    return pal_temperate;   /* default / temperate */
}

const char *auto_style(double lat)
{
    double a = fabs(lat);
    if (a < 23.5) return "tropical";
    if (a < 35.0) return "subtropical";
    if (a < 60.0) return "temperate";
    if (a < 75.0) return "subpolar";
    return "polar";
}

/* ── colour-rule building ─────────────────────────────────────────────── */

/*
 * build_color_rules() — write a GRASS r.colors rules file.
 *
 * The depth-index range is nominally 0–1000.  Wave texture, turbulence and
 * foam can push values outside this range:
 *   • lo_clamp  : lowest possible value (usually negative); clamped to shore
 *                 colour (palette[0]).  Accounts for wave_amp + foam depth.
 *   • hi_clamp  : highest possible value (> 1000); clamped to abyss colour
 *                 (palette[N_STOPS-1]).
 * Null cells → white (255:255:255) for transparency over land base maps.
 */
static void build_color_rules(const PaletteStop *pal, int n,
                               double lo_clamp, double hi_clamp,
                               FILE *fp)
{
    /* lower clamp — maps to shore colour */
    fprintf(fp, "%.0f %d:%d:%d\n",
            lo_clamp - 1.0,
            pal[0].r, pal[0].g, pal[0].b);
    /* palette stops */
    for (int i = 0; i < n; i++)
        fprintf(fp, "%d %d:%d:%d\n",
                pal[i].idx, pal[i].r, pal[i].g, pal[i].b);
    /* upper clamp — maps to abyss colour */
    fprintf(fp, "%.0f %d:%d:%d\n",
            hi_clamp + 1.0,
            pal[n-1].r, pal[n-1].g, pal[n-1].b);
    /* null → white */
    fprintf(fp, "nv 255:255:255\n");
}

static void apply_color_rules(const char *map, const PaletteStop *pal,
                               double lo_clamp, double hi_clamp)
{
    char rules_path[GPATH_MAX];
    FILE *fp;

    G_temp_element("i.ocean");
    snprintf(rules_path, sizeof(rules_path), "%s/i.ocean_rules_%d.txt",
             G_tempfile(), (int)getpid());

    fp = fopen(rules_path, "w");
    if (!fp)
        G_fatal_error(_("Cannot write colour rules to %s"), rules_path);
    build_color_rules(pal, N_STOPS, lo_clamp, hi_clamp, fp);
    fclose(fp);

    {
        char cmd[GPATH_MAX + 256];
        snprintf(cmd, sizeof(cmd),
                 "r.colors map=%s rules=%s --quiet", map, rules_path);
        if (system(cmd) != 0)
            G_warning(_("r.colors failed — colour table not applied"));
    }
    remove(rules_path);
}

/* ── pixel-size helpers ───────────────────────────────────────────────── */

static double pixel_size_m(const struct Cell_head *r, const char *units,
                            const char *proj)
{
    double scale = 1.0;
    int geographic = 0;

    if (strstr(units, "degree") || strcmp(proj, "ll") == 0 ||
        strcmp(proj, "latlong") == 0 || strcmp(proj, "longlat") == 0)
        geographic = 1;

    if (geographic) {
        double lat_rad = (r->north + r->south) / 2.0 * M_PI / 180.0;
        double ns_m = r->ns_res * 111320.0;
        double ew_m = r->ew_res * 111320.0 * cos(lat_rad);
        return (ns_m + ew_m) / 2.0;
    }
    if (strstr(units, "feet") || strstr(units, "foot"))
        scale = 0.3048;
    else if (strstr(units, "kilometer") || strcmp(units, "km") == 0)
        scale = 1000.0;
    return ((r->ns_res + r->ew_res) / 2.0) * scale;
}

static double center_latitude(const struct Cell_head *r,
                               const char *units, const char *proj)
{
    char cmd[512], buf[256];
    FILE *fp;

    if (strstr(units, "degree") || strcmp(proj, "ll") == 0 ||
        strcmp(proj, "latlong") == 0 || strcmp(proj, "longlat") == 0)
        return (r->north + r->south) / 2.0;

    /* projected: call m.proj to reproject the centre point */
    snprintf(cmd, sizeof(cmd),
             "m.proj coordinates=%.6f,%.6f -od --quiet",
             (r->east + r->west) / 2.0,
             (r->north + r->south) / 2.0);
    fp = popen(cmd, "r");
    if (fp) {
        if (fgets(buf, sizeof(buf), fp)) {
            double lon, lat;
            if (sscanf(buf, "%lf|%lf", &lon, &lat) == 2) {
                pclose(fp);
                return lat;
            }
        }
        pclose(fp);
    }
    G_warning(_("Could not detect latitude from CRS; assuming 45°N (temperate)."));
    return 45.0;
}

/* ── Rbuf helpers (shared across all .c files) ────────────────────────── */

Rbuf *rbuf_alloc(int nrows, int ncols)
{
    Rbuf *rb = G_malloc(sizeof(Rbuf));
    rb->nrows = nrows;
    rb->ncols = ncols;
    rb->data  = G_malloc(sizeof(DCELL) * (size_t)nrows * ncols);
    return rb;
}

void rbuf_free(Rbuf *rb)
{
    if (rb) {
        G_free(rb->data);
        G_free(rb);
    }
}

Rbuf *rbuf_read(const char *name)
{
    int fd;
    struct Cell_head region;
    Rbuf *rb;
    int r;

    G_get_set_window(&region);
    rb = rbuf_alloc(region.rows, region.cols);
    fd = Rast_open_old(name, "");

    for (r = 0; r < region.rows; r++)
        Rast_get_d_row(fd, rb->data + (size_t)r * region.cols, r);

    Rast_close(fd);
    return rb;
}

void rbuf_write(const Rbuf *rb, const char *name, RASTER_MAP_TYPE type)
{
    int fd = Rast_open_new(name, type);
    for (int r = 0; r < rb->nrows; r++)
        Rast_put_d_row(fd, rb->data + (size_t)r * rb->ncols);
    Rast_close(fd);
}

/* ── temp-map registry ────────────────────────────────────────────────── */

#define MAX_TMP 64
static char tmp_names[MAX_TMP][GNAME_MAX];
static int  tmp_count = 0;

void tmp_register(const char *name)
{
    if (tmp_count < MAX_TMP)
        G_strlcpy(tmp_names[tmp_count++], name, GNAME_MAX);
}

void tmp_cleanup(void)
{
    if (tmp_count == 0)
        return;

    /* build comma-separated name list for g.remove */
    char namelist[MAX_TMP * GNAME_MAX];
    int pos = 0;
    for (int i = 0; i < tmp_count; i++) {
        if (i > 0)
            namelist[pos++] = ',';
        int len = strlen(tmp_names[i]);
        memcpy(namelist + pos, tmp_names[i], len);
        pos += len;
    }
    namelist[pos] = '\0';

    char cmd[sizeof(namelist) + 128];
    snprintf(cmd, sizeof(cmd),
             "g.remove -f type=raster name=%s --quiet 2>/dev/null",
             namelist);
    system(cmd);
    tmp_count = 0;
}

/* ── latitude-gradient (applied in-place) ────────────────────────────── */

/*
 * latitude_gradient() — shift depth index by ±80 over the N-S extent of
 * the scene to mimic the visible SST gradient within the latitude band.
 * Equatorward → lighter (subtract); poleward → darker (add).
 * The sign is reversed in the southern hemisphere.
 *
 * Applied with an OpenMP parallel row loop.
 */
static void latitude_gradient(Rbuf *work, const Rbuf *mask,
                               double latitude, int nprocs)
{
    double shift = 80.0;
    int nrows = work->nrows, ncols = work->ncols;

    omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
    for (int r = 0; r < nrows; r++) {
        double frac = (double)r / nrows;   /* 0 = top, 1 = bottom */
        double delta = (latitude >= 0.0)
                       ? -shift * frac          /* NH: warmer southward  */
                       :  shift * frac;         /* SH: warmer northward  */
        for (int c = 0; c < ncols; c++) {
            if (Rast_is_d_null_value(&RBUF_CELL(mask, r, c)))
                continue;
            RBUF_CELL(work, r, c) += delta;
        }
    }
}

/* ── main ─────────────────────────────────────────────────────────────── */

int main(int argc, char *argv[])
{
    struct GModule *module;

    /* options */
    struct Option *opt_input, *opt_output, *opt_depth,
                  *opt_ocean_val, *opt_depth_min,
                  *opt_latitude, *opt_style;
    /* flags */
    struct Flag *flag_w, *flag_s, *flag_l, *flag_f, *flag_d;

    G_gisinit(argv[0]);

    /* ── module description ─────────────────────────────────────────── */
    module = G_define_module();
    G_add_keyword("imagery");
    G_add_keyword("raster");
    G_add_keyword("ocean");
    G_add_keyword("visualization");
    G_add_keyword("color");
    G_add_keyword("depth");
    G_add_keyword("bathymetry");
    module->description =
        _("Renders ocean areas of a raster map with a realistic visual appearance");

    /* ── options ────────────────────────────────────────────────────── */
    opt_input = G_define_standard_option(G_OPT_R_INPUT);
    opt_input->key         = "input";
    opt_input->label       = _("Ocean mask raster");
    opt_input->description = _("Raster whose cells equal to ocean_value "
                               "are treated as ocean; mutually optional with depth");
    opt_input->required    = NO;

    opt_output = G_define_standard_option(G_OPT_R_OUTPUT);
    opt_output->key         = "output";
    opt_output->label       = _("Output ocean visualization raster");
    opt_output->description = _("Floating-point raster carrying depth-index "
                                "values and ocean colour rules");

    opt_depth = G_define_standard_option(G_OPT_R_INPUT);
    opt_depth->key         = "depth";
    opt_depth->label       = _("Bathymetric depth raster (positive = metres below surface)");
    opt_depth->description = _("Drives the colour gradient; also used as ocean "
                               "source when no input mask is given; "
                               "integer (CELL) maps are smoothed automatically");
    opt_depth->required    = NO;

    opt_ocean_val = G_define_option();
    opt_ocean_val->key         = "ocean_value";
    opt_ocean_val->type        = TYPE_DOUBLE;
    opt_ocean_val->label       = _("Cell value in the input mask that represents ocean");
    opt_ocean_val->description = _("Cells with this value are treated as ocean; "
                                   "ignored when only depth is supplied");
    opt_ocean_val->answer      = "1";
    opt_ocean_val->required    = NO;

    opt_depth_min = G_define_option();
    opt_depth_min->key         = "depth_min";
    opt_depth_min->type        = TYPE_DOUBLE;
    opt_depth_min->label       = _("Minimum depth value to treat as ocean");
    opt_depth_min->description = _("Cells in the depth map with value >= depth_min "
                                   "become ocean; use 0 for GEBCO-style maps, "
                                   "a negative value for elevation/DEM maps");
    opt_depth_min->answer      = "0";
    opt_depth_min->required    = NO;

    opt_latitude = G_define_option();
    opt_latitude->key         = "latitude";
    opt_latitude->type        = TYPE_DOUBLE;
    opt_latitude->label       = _("Reference latitude (decimal degrees, -90 to 90)");
    opt_latitude->description = _("Auto-detected from the map centre when omitted");
    opt_latitude->required    = NO;

    opt_style = G_define_option();
    opt_style->key         = "style";
    opt_style->type        = TYPE_STRING;
    opt_style->label       = _("Ocean colour style");
    opt_style->description = _("Overrides the auto-detected latitude zone");
    opt_style->answer      = "auto";
    opt_style->options     = "auto,tropical,subtropical,temperate,subpolar,polar";
    opt_style->descriptions =
        "auto;Derive style from latitude;"
        "tropical;Warm turquoise (< 23.5 deg);"
        "subtropical;Blue-green (23.5-35 deg);"
        "temperate;Cool blue (35-60 deg);"
        "subpolar;Steel blue (60-75 deg);"
        "polar;Icy pale blue (> 75 deg)";
    opt_style->required    = NO;
    opt_style->guisection  = "Palette";

    /* ── flags ──────────────────────────────────────────────────────── */
    flag_w = G_define_flag();
    flag_w->key         = 'w';
    flag_w->label       = _("Add wave / ripple texture");
    flag_w->description = _("Overlays sinusoidal wave patterns scaled "
                            "to the current pixel size");
    flag_w->guisection  = "Effects";

    flag_s = G_define_flag();
    flag_s->key         = 's';
    flag_s->label       = _("Simulate depth from distance to shore");
    flag_s->description = _("When no depth map is given, estimates depth "
                            "using distance transform from land pixels");
    flag_s->guisection  = "Effects";

    flag_l = G_define_flag();
    flag_l->key         = 'l';
    flag_l->label       = _("Apply subtle north-south warmth gradient");
    flag_l->description = _("Mimics the visible SST gradient within "
                            "the scene latitude band");
    flag_l->guisection  = "Effects";

    flag_f = G_define_flag();
    flag_f->key         = 'f';
    flag_f->label       = _("Add multi-scale turbulence and coastal foam");
    flag_f->description =
        _("Generates a physically-based ocean surface using the Phillips "
          "spectrum (FFTW3) with multi-octave frequency content, "
          "coastal whitecap foam intensified near the shore, and "
          "directional swell banding — replicating the painted-water "
          "appearance of realistic cartographic ocean renderings");
    flag_f->guisection  = "Effects";

    flag_d = G_define_flag();
    flag_d->key         = 'd';
    flag_d->label       = _("Add maritime decorations");
    flag_d->description =
        _("Scatters sea wildlife (whales, dolphins, fish schools, octopus) "
          "and maritime objects (treasure chests, floating crates, "
          "sunken ships) across the ocean surface");
    flag_d->guisection  = "Effects";

    /* mutual requirement: input or depth must be given */
    G_option_required(opt_input, opt_depth, NULL);

    if (G_parser(argc, argv))
        exit(EXIT_FAILURE);

    /* ── resolve parameters ─────────────────────────────────────────── */
    const char *input_name  = opt_input->answer;
    const char *output_name = opt_output->answer;
    const char *depth_name  = opt_depth->answer;
    double ocean_val  = atof(opt_ocean_val->answer);
    double depth_min  = atof(opt_depth_min->answer);

    int nprocs = omp_get_max_threads();

    char pid_str[32];
    snprintf(pid_str, sizeof(pid_str), "%d", (int)getpid());

    /* ── progress tracking ──────────────────────────────────────────── */
    int total_steps = 6
        + (flag_l->answer ? 1 : 0)
        + (flag_w->answer ? 1 : 0)
        + (flag_f->answer ? 1 : 0)
        + (flag_d->answer ? 1 : 0);
    int step = 0;
#define PROGRESS(msg) do { \
        G_percent(++step, total_steps, 1); \
        G_message("%s", _(msg)); \
    } while (0)

    /* ── region + CRS ───────────────────────────────────────────────── */
    PROGRESS("Analysing region and CRS...");
    struct Cell_head region;
    G_get_set_window(&region);

    /* read projection info via g.proj -g */
    char proj_buf[4096] = {0};
    {
        FILE *fp = popen("g.proj -g --quiet 2>/dev/null", "r");
        if (fp) {
            fread(proj_buf, 1, sizeof(proj_buf) - 1, fp);
            pclose(fp);
        }
    }
    /* extract 'units' and 'proj' keys */
    char units_val[64] = {0}, proj_val[64] = {0};
    {
        char *p;
        if ((p = strstr(proj_buf, "units=")))
            sscanf(p + 6, "%63[^\n]", units_val);
        if ((p = strstr(proj_buf, "proj=")))
            sscanf(p + 5, "%63[^\n]", proj_val);
    }
    /* lowercase units */
    for (char *p = units_val; *p; p++)
        *p = (char)tolower((unsigned char)*p);

    double pixel_m  = pixel_size_m(&region, units_val, proj_val);
    G_verbose_message(_("Mean pixel size: %.1f m"), pixel_m);

    double latitude;
    if (opt_latitude->answer)
        latitude = atof(opt_latitude->answer);
    else
        latitude = center_latitude(&region, units_val, proj_val);
    G_verbose_message(_("Reference latitude: %.2f deg"), latitude);

    const char *style = strcmp(opt_style->answer, "auto") == 0
                        ? auto_style(latitude)
                        : opt_style->answer;
    G_message(_("Ocean colour style: %s  (latitude %.1f deg)"), style, latitude);

    /* ── ocean mask ─────────────────────────────────────────────────── */
    PROGRESS("Extracting ocean pixels...");

    Rbuf *mask = NULL;
    if (input_name) {
        Rbuf *inp = rbuf_read(input_name);
        mask = make_mask_from_value(inp, ocean_val,
                                    region.rows, region.cols);
        rbuf_free(inp);
    } else {
        Rbuf *dep = rbuf_read(depth_name);
        G_message(_("No input mask — deriving ocean extent from depth map "
                    "(depth >= %.3g)..."), depth_min);
        mask = make_mask_from_depth(dep, depth_min,
                                    region.rows, region.cols);
        rbuf_free(dep);
    }

    /* sanity: count ocean pixels */
    {
        long n_ocean = 0;
        for (int i = 0; i < region.rows * region.cols; i++)
            if (!Rast_is_d_null_value(&mask->data[i]))
                n_ocean++;
        if (n_ocean == 0)
            G_fatal_error(_("No ocean pixels found — check input/depth_min."));
        G_verbose_message(_("%ld ocean pixels"), n_ocean);
    }

    /* ── depth index ────────────────────────────────────────────────── */
    Rbuf *depth_idx = NULL;
    if (depth_name) {
        PROGRESS("Building depth index from bathymetric map...");
        depth_idx = depth_from_map(depth_name, mask, nprocs, pid_str);
    } else if (flag_s->answer) {
        PROGRESS("Estimating depth from distance to shore...");
        depth_idx = depth_from_shore(mask, pixel_m, nprocs, pid_str);
    } else {
        PROGRESS("Setting flat depth index (500)...");
        depth_idx = depth_flat(mask, 500.0, nprocs);
    }

    /* ── latitude gradient ──────────────────────────────────────────── */
    if (flag_l->answer) {
        PROGRESS("Applying north-south temperature gradient...");
        latitude_gradient(depth_idx, mask, latitude, nprocs);
    }

    /* ── wave texture ───────────────────────────────────────────────── */
    double wave_amp = 0.0;
    if (flag_w->answer) {
        PROGRESS("Generating wave texture...");
        Rbuf *wv = wave_pattern(mask, pixel_m, &wave_amp, nprocs);
        /* blend wave into depth_idx */
        long npx = (long)region.rows * region.cols;
        omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
        for (long i = 0; i < npx; i++) {
            if (!Rast_is_d_null_value(&mask->data[i]))
                depth_idx->data[i] += wv->data[i];
        }
        rbuf_free(wv);
    }

    /* ── foam / turbulence (Phillips spectrum, FFTW3) ───────────────── */
    double lo_clamp = -(wave_amp + 1.0);
    if (flag_f->answer) {
        PROGRESS("Generating Phillips-spectrum turbulence and coastal foam...");
        double foam_lo = lo_clamp;
        Rbuf *foamy = foam_turbulence(depth_idx, mask,
                                      pixel_m, region.ns_res, region.ew_res,
                                      &foam_lo, nprocs, pid_str);
        rbuf_free(depth_idx);
        depth_idx = foamy;
        if (foam_lo < lo_clamp)
            lo_clamp = foam_lo;
    }

    /* ── maritime decorations ───────────────────────────────────────── */
    if (flag_d->answer) {
        PROGRESS("Placing maritime decorations...");
        place_decorations(depth_idx, mask, pixel_m, lo_clamp);
        /* decorations can push values lower; keep clamp tight */
        if (lo_clamp > -600.0)
            lo_clamp = -600.0;
    }

    /* ── write output ───────────────────────────────────────────────── */
    PROGRESS("Writing output raster...");
    /* null-mask the result before writing */
    {
        long npx = (long)region.rows * region.cols;
        omp_set_num_threads(nprocs);
#pragma omp parallel for schedule(static)
        for (long i = 0; i < npx; i++) {
            if (Rast_is_d_null_value(&mask->data[i]))
                Rast_set_d_null_value(&depth_idx->data[i], 1);
        }
    }
    rbuf_write(depth_idx, output_name, DCELL_TYPE);

    /* ── colour rules ───────────────────────────────────────────────── */
    PROGRESS("Applying ocean colour palette...");
    {
        const PaletteStop *pal = get_palette(style);
        double hi_clamp = 1000.0 + wave_amp + 1.0;
        apply_color_rules(output_name, pal, lo_clamp, hi_clamp);
    }

    /* ── metadata ───────────────────────────────────────────────────── */
    PROGRESS("Writing metadata...");
    {
        char title[256], desc[512];
        const char *src = input_name ? input_name : depth_name;
        snprintf(title, sizeof(title),
                 "Ocean visualisation of %s", src ? src : "(none)");

        char eff[128] = "";
        if (depth_name) strncat(eff, "depth-map,", sizeof(eff)-1);
        else if (flag_s->answer) strncat(eff, "shore-dist,", sizeof(eff)-1);
        if (flag_w->answer) strncat(eff, "waves,", sizeof(eff)-1);
        if (flag_l->answer) strncat(eff, "lat-gradient,", sizeof(eff)-1);
        if (flag_f->answer) strncat(eff, "foam,", sizeof(eff)-1);
        if (flag_d->answer) strncat(eff, "decorations,", sizeof(eff)-1);
        /* strip trailing comma */
        int el = strlen(eff);
        if (el > 0 && eff[el-1] == ',')
            eff[el-1] = '\0';

        snprintf(desc, sizeof(desc),
                 "style=%s, lat=%.1fdeg, pixel=%.0fm, effects=[%s]",
                 style, latitude, pixel_m, eff[0] ? eff : "none");

        char cmd[1024];
        snprintf(cmd, sizeof(cmd),
                 "r.support map=%s title=\"%s\" description=\"%s\" --quiet",
                 output_name, title, desc);
        system(cmd);

        snprintf(cmd, sizeof(cmd),
                 "r.support map=%s history=\"Generated by i.ocean\" --quiet",
                 output_name);
        system(cmd);
    }

    /* ── done ───────────────────────────────────────────────────────── */
    G_percent(1, 1, 1);
    G_message(_("Done. Ocean visualisation written to <%s>."), output_name);
    G_message(_("Display with:  d.rast map=%s"), output_name);

    /* cleanup */
    rbuf_free(mask);
    rbuf_free(depth_idx);
    tmp_cleanup();

    return EXIT_SUCCESS;
}
