/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The Simulationarchive field descriptors for struct cft_ias15_state,
 * as macros, so that exactly one file defines the lists and everything
 * else refers to them.
 *
 * They were in src/cft_archive.c as file-scope statics, which meant the
 * registered integrator could not name them: src/cft_ias15_fields.c is
 * linked into every drop-in target and src/cft_archive.c is not, so
 * registration installed an empty list and REBOUND could not resolve a
 * single cft_ field on read. tests/gate_real.c is the gate that caught
 * it.
 *
 * Adding a blob means adding one line to CFT_FD_BLOBS and raising
 * CFT_N_BLOBS; both halves are checked by cft_archive_selftest().
 */
#ifndef CFT_IAS15_FIELDS_H
#define CFT_IAS15_FIELDS_H

#include <stddef.h>
#include "rebound.h"
#include "binarydata.h"
#include "cft_ias15.h"

#define CFT_OFF(m)   offsetof(struct cft_ias15_state, m)

#define CFT_FD_P(NAME, MEMBER, W) \
    { "", REB_POINTER, "cft_" NAME, CFT_OFF(MEMBER), CFT_OFF(n_elem), (W), 0 },

#define CFT_FD_L(NAME, MEMBER, W) \
    CFT_FD_P(NAME "0", MEMBER[0], W) CFT_FD_P(NAME "1", MEMBER[1], W) \
    CFT_FD_P(NAME "2", MEMBER[2], W) CFT_FD_P(NAME "3", MEMBER[3], W) \
    CFT_FD_P(NAME "4", MEMBER[4], W) CFT_FD_P(NAME "5", MEMBER[5], W) \
    CFT_FD_P(NAME "6", MEMBER[6], W)

/* the wide state */
#define CFT_FD_BLOBS(W) \
    CFT_FD_P("x0",   x0,   W) \
    CFT_FD_P("v0",   v0,   W) \
    CFT_FD_P("a0",   a0,   W) \
    CFT_FD_P("csx",  csx,  W) \
    CFT_FD_P("csv",  csv,  W) \
    CFT_FD_P("csa0", csa0, W) \
    CFT_FD_L("g",   g,   W) \
    CFT_FD_L("b",   b,   W) \
    CFT_FD_L("csb", csb, W) \
    CFT_FD_L("e",   e,   W) \
    CFT_FD_L("br",  br,  W) \
    CFT_FD_L("er",  er,  W) \
    CFT_FD_P("x",   x,   W) \
    CFT_FD_P("v",   v,   W)

#define CFT_FD_S(NAME, TYPE, MEMBER) \
    { "", TYPE, "cft_" NAME, CFT_OFF(MEMBER), 0, 0, 0 },

/* Configuration and provenance. cft_abi is a char[16] inside the
 * struct, not a char*: REB_STRING dereferences the field as a pointer
 * and REB_POINTER does the same, and no simple dtype is 16 bytes wide,
 * so the two halves go out as uint64 words at their own offsets. The
 * bytes on disk are the string's bytes, in order. */
#define CFT_FD_SCALARS \
    CFT_FD_S("epsilon",          REB_DOUBLE, epsilon) \
    CFT_FD_S("min_dt",           REB_DOUBLE, min_dt) \
    CFT_FD_S("adaptive_mode",    REB_INT,    adaptive_mode) \
    CFT_FD_S("format",           REB_INT,    format) \
    CFT_FD_S("max_iter",         REB_INT,    max_iter) \
    CFT_FD_S("arith_fma",        REB_INT,    arith_fma) \
    CFT_FD_S("E",                REB_SIZE_T, E) \
    CFT_FD_S("abi_0",            REB_UINT64, cft_abi[0]) \
    CFT_FD_S("abi_1",            REB_UINT64, cft_abi[8]) \
    CFT_FD_S("constants_digest", REB_UINT64, constants_digest) \
    CFT_FD_S("n_elem",           REB_SIZE_T, n_elem)      /* LAST: see NOTE 3 */

#define CFT_FD_LIST(W) { CFT_FD_BLOBS(W) CFT_FD_SCALARS { 0 } }

#define CFT_N_BLOBS 50

/* One definition of each, in src/cft_ias15_fields.c. cft_archive.c
 * selects among them by format; the integrator registers the first,
 * which is why a load resolves names at every format - see the note in
 * that file. */
extern const struct reb_binarydata_field_descriptor cft_fd_fp64[];
extern const struct reb_binarydata_field_descriptor cft_fd_fp128[];
extern const struct reb_binarydata_field_descriptor cft_fd_fp256[];

#endif /* CFT_IAS15_FIELDS_H */
