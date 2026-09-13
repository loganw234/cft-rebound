/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * What libcft itself says when asked for a format the image lacks -
 * the second of the three refusal layers docs/BITSTREAM.md describes,
 * exercised directly rather than through this repo's own check, so the
 * record shows the library's status AND the sentence it does or does
 * not attach. Built and run by hw/cardtest-f128.sh:
 *
 *   refusal_probe <artifact> <fp32|fp64|fp128|fp256>
 *
 * Exit 0 when the format is refused everywhere it should be (elementwise
 * run, reduction, program run) with nothing written, 1 when any of them
 * ran, 2 on a setup failure. A refused run must leave the output bytes
 * exactly as they were: a "refusal" that wrote zeros would be the worst
 * possible shape of wrong answer (cft-fp256 host/src/device.c says why).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cft.h"

static const char *q(const char *s){ return s && *s ? s : "(empty)"; }

int main(int argc, char **argv){
    cft_device *dev;
    cft_status st;
    cft_format fmt;
    cft_caps caps;
    unsigned char a[4 * 32], b[4 * 32], c[4 * 32], d[4 * 32], d0[4 * 32];
    uint32_t fl = 0;
    int ran = 0, f;

    if (argc != 3){ fprintf(stderr, "usage: refusal_probe <artifact> <format>\n"); return 2; }
    if (!strcmp(argv[2], "fp32")) fmt = CFT_FP32;
    else if (!strcmp(argv[2], "fp64")) fmt = CFT_FP64;
    else if (!strcmp(argv[2], "fp128")) fmt = CFT_FP128;
    else if (!strcmp(argv[2], "fp256")) fmt = CFT_FP256;
    else { fprintf(stderr, "format fp32|fp64|fp128|fp256\n"); return 2; }

    st = cft_open(argv[1], 0, &dev);
    if (st != CFT_OK){ fprintf(stderr, "cft_open: %s (%s)\n", cft_strerror(st), q(cft_last_error())); return 2; }
    memset(&caps, 0, sizeof caps); caps.struct_size = sizeof caps;
    if (cft_get_caps(dev, &caps) != CFT_OK){ fprintf(stderr, "cft_get_caps failed\n"); return 2; }
    printf("artifact %s: CAPS[3:0] = 0x%x, formats", argv[1], (unsigned)(caps.format_mask & 0xFu));
    for (f = 0; f < 4; f++) if (caps.format_mask & (1u << f)) printf(" %s", cft_format_name((cft_format)f));
    printf("; cft_supports(FMA, %s) = %d\n", cft_format_name(fmt), cft_supports(dev, CFT_FMA, fmt));

    /* Operands that would be finite at every format; the output filled
     * with a pattern a run could never leave behind. */
    memset(a, 0, sizeof a); memset(b, 0, sizeof b); memset(c, 0, sizeof c);
    memset(d, 0xA5, sizeof d); memcpy(d0, d, sizeof d);

    st = cft_run(dev, CFT_FMA, fmt, CFT_RNE, a, b, c, d, 4, &fl, NULL);
    printf("cft_run(FMA, %s, n=4): status %d \"%s\", cft_last_error \"%s\", output %s\n",
           cft_format_name(fmt), (int)st, cft_strerror(st), q(cft_last_error()),
           memcmp(d, d0, sizeof d) ? "WRITTEN" : "untouched");
    if (st == CFT_OK) ran = 1;

    memcpy(d, d0, sizeof d);
    st = cft_reduce(dev, CFT_SUM, fmt, CFT_RNE, a, NULL, d, 4, &fl, NULL);
    printf("cft_reduce(SUM, %s, n=4): status %d \"%s\", cft_last_error \"%s\", output %s\n",
           cft_format_name(fmt), (int)st, cft_strerror(st), q(cft_last_error()),
           memcmp(d, d0, sizeof d) ? "WRITTEN" : "untouched");
    if (st == CFT_OK) ran = 1;

    cft_close(dev);
    printf("%s: %s\n", cft_format_name(fmt), ran ? "RAN" : "refused everywhere, nothing written");
    return ran ? 1 : 0;
}
