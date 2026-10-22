/****************************************************************************
 *
 * MODULE:       i.ocean  —  ram_check.c
 * PURPOSE:      Estimate peak memory usage and compare against available RAM.
 *               Prevents segfaults caused by over-allocation in large regions.
 *
 * OVERVIEW
 * ────────
 * Before allocating the main raster buffers the module estimates how many
 * bytes it will need across its lifetime, reads the Linux /proc/meminfo
 * MemAvailable counter, and returns one of four status levels:
 *
 *   RAM_OK     (<= 60 % of available)  — proceed normally
 *   RAM_WARN   (60–80 %)               — warn user, proceed
 *   RAM_TIGHT  (80–95 %)               — warn, disable foam (-f) if active
 *   RAM_FATAL  (>= 95 %)               — G_fatal_error with advice
 *
 * BUFFER ACCOUNTING (worst case, simultaneous peak)
 * ──────────────────────────────────────────────────
 *   Constant:
 *     2 × Rbuf (mask + depth_idx)              = 2 × N × 8 bytes
 *     1 × Rbuf temp (input or depth read)      = 1 × N × 8 bytes
 *   -w flag:    1 × Rbuf (wave)                = 1 × N × 8 bytes
 *   -s flag:    1 × float* shore dist (DT)     = 1 × N × 4 bytes
 *   -f flag:    fftw_complex + real output      ≈ 3 × N × 8 bytes
 *               + 1 × float* shore dist        = 1 × N × 4 bytes
 *   -d flag:    1 × long* ocean_idx            ≤ 1 × N × 8 bytes
 *
 ****************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <grass/gis.h>
#include <grass/glocale.h>

#include "local_proto.h"

/* ── read MemAvailable from /proc/meminfo (Linux) ────────────────────── */

static long read_mem_available_kb(void)
{
    FILE *fp = fopen("/proc/meminfo", "r");
    if (!fp)
        return 0;

    long avail = 0;
    char line[256];
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "MemAvailable:", 13) == 0) {
            sscanf(line + 13, " %ld", &avail);
            break;
        }
    }
    fclose(fp);
    return avail;   /* kB */
}

/* ── public entry point ───────────────────────────────────────────────── */

/*
 * ram_check() — estimate peak memory and compare against available RAM.
 *
 * Parameters:
 *   nrows, ncols  : computational region dimensions
 *   has_w/f/s/d   : which optional flags are active
 *   est_mb_out    : (out) estimated peak usage in MB
 *   avail_mb_out  : (out) available RAM in MB (0 if not readable)
 *
 * Returns one of: RAM_OK, RAM_WARN, RAM_TIGHT, RAM_FATAL.
 */
RamStatus ram_check(long nrows, long ncols,
                    int has_w, int has_f, int has_s, int has_d,
                    double *est_mb_out, double *avail_mb_out)
{
    long N = nrows * ncols;

    /* count 8-byte (DCELL / double) buffer equivalents */
    /* base: mask + depth_idx + 1 temp read buffer */
    double n8 = 3.0;
    double n4 = 0.0;   /* 4-byte float buffers */
    long   nlng = 0;   /* long integer arrays */

    if (has_w) n8 += 1.0;                    /* wave Rbuf */
    if (has_s) n4 += (double)N;              /* shore dist float* */
    if (has_f) {
        /* fftw_complex (16B each) at nrows*(ncols/2+1) ≈ N/2 */
        n8 += 1.0 * (ncols / 2 + 1.0) / ncols * 2.0;
        n8 += 1.0;                           /* spatial DCELL output */
        n4 += (double)N;                     /* shore dist float* */
    }
    if (has_d) nlng += N;                    /* ocean_idx long* */

    double est_bytes = n8 * N * 8.0
                     + n4 * 4.0
                     + (double)nlng * sizeof(long);

    double est_mb = est_bytes / (1024.0 * 1024.0);

    long avail_kb = read_mem_available_kb();
    double avail_mb = (double)avail_kb / 1024.0;

    if (est_mb_out)   *est_mb_out   = est_mb;
    if (avail_mb_out) *avail_mb_out = avail_mb;

    if (avail_mb < 1.0) {
        /* /proc/meminfo not available (non-Linux); just estimate */
        G_verbose_message(_("RAM check: estimated %.0f MB needed "
                            "(system RAM unknown)."), est_mb);
        return RAM_OK;
    }

    double ratio = est_mb / avail_mb;

    if (ratio >= 0.95) {
        G_fatal_error(
            _("Insufficient RAM: estimated peak usage is %.0f MB but "
              "only %.0f MB is available.\n"
              "Reduce the computational region (g.region) or omit "
              "memory-intensive flags (-f uses FFTW3 arrays)."),
            est_mb, avail_mb);
        /* unreachable */
        return RAM_FATAL;
    }
    if (ratio >= 0.80) {
        G_warning(
            _("Low memory: estimated %.0f MB needed vs %.0f MB available "
              "(%.0f%% of RAM).\n"
              "Disabling Phillips-spectrum foam (-f) to reduce peak usage.\n"
              "Use g.region to reduce the region, or add more RAM."),
            est_mb, avail_mb, ratio * 100.0);
        return RAM_TIGHT;
    }
    if (ratio >= 0.60) {
        G_warning(
            _("Memory usage: estimated %.0f MB (%.0f%% of %.0f MB available)."
              " Processing may be slower if swapping occurs."),
            est_mb, ratio * 100.0, avail_mb);
        return RAM_WARN;
    }

    G_verbose_message(_("RAM check: %.0f MB estimated, %.0f MB available "
                        "(%.0f%%)."),
                      est_mb, avail_mb, ratio * 100.0);
    return RAM_OK;
}
