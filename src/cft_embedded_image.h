/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * A sequencer program image compiled into the binary, for the library
 * half, which loads nothing from disk. Each family's table is generated
 * into build/ by tools/embed_programs.py and includes this header; the
 * type is defined here once, so that two families can share one
 * translation unit. */
#ifndef CFT_EMBEDDED_IMAGE_H
#define CFT_EMBEDDED_IMAGE_H
#include <stddef.h>
struct cft_embedded_image {
    const char *stem;            /* the .cfta's name before -<format> */
    const char *format;          /* fp32, fp64, fp128 or fp256 */
    const unsigned char *bytes;  /* the assembled image */
    size_t length;
};
#endif
