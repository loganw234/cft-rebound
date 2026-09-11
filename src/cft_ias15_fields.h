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
 * CFT_N_BLOBS; both halves are checked by cft_archive_selftest(). The
 * same goes for CFT_FD_ALIAS and CFT_N_ALIAS_BLOBS, and for
 * CFT_FD_SCALARS and CFT_N_SCALARS. Nothing outside this file may
 * write any of those three numbers down.
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

/* --------------------------------------------------------------------
 * The alias family: a SECOND blob family, counted separately.
 *
 * These are the arrays REBOUND holds at N_allocated rather than at the
 * live 3N - the six coefficient families read at the MARK's stride,
 * which is what ias15_engine_alias_resize() keeps in its shadow, and
 * the two compensated-summation carries whose tail a regrow reads back.
 * Nothing carried either until 2026-09-11; see the state struct.
 *
 * They are not part of CFT_FD_BLOBS and not reachable through
 * cft_archive_state_blob(), for one reason: their length is
 * hiwater_n_elem and not n_elem. Everything that walks "the blobs" -
 * cft_archive_state_alloc(), the stub in tests/cft_shim_stub.c, the
 * "did it reach memory" check - means blobs of n_elem, and pretending
 * these were more of those would have sized them wrong in four places.
 *
 * ORDER. g, b, e, csb, er, br: ias15_engine_alias_resize()'s fam[6],
 * which is the order the shim publishes into state->alias[][]. The
 * order is arbitrary as long as the two agree - what matters is that a
 * level of one family is never read as a level of another - and
 * cft_archive_selftest() checks the descriptors land on alias[0][0]
 * through alias[5][6] in sequence, so a list edited out of step with
 * the member is caught rather than silently swapping two families.
 * -------------------------------------------------------------------- */
#define CFT_FD_A(NAME, MEMBER, W) \
    { "", REB_POINTER, "cft_alias_" NAME, CFT_OFF(MEMBER), CFT_OFF(hiwater_n_elem), (W), 0 },

#define CFT_FD_AL(NAME, F, W) \
    CFT_FD_A(NAME "0", alias[F][0], W) CFT_FD_A(NAME "1", alias[F][1], W) \
    CFT_FD_A(NAME "2", alias[F][2], W) CFT_FD_A(NAME "3", alias[F][3], W) \
    CFT_FD_A(NAME "4", alias[F][4], W) CFT_FD_A(NAME "5", alias[F][5], W) \
    CFT_FD_A(NAME "6", alias[F][6], W)

#define CFT_FD_ALIAS(W) \
    CFT_FD_AL("g",   0, W) \
    CFT_FD_AL("b",   1, W) \
    CFT_FD_AL("e",   2, W) \
    CFT_FD_AL("csb", 3, W) \
    CFT_FD_AL("er",  4, W) \
    CFT_FD_AL("br",  5, W) \
    CFT_FD_A("csx", alias_csx, W) \
    CFT_FD_A("csv", alias_csv, W)

/* The name every alias blob begins with, once, so that the classifier
 * in src/cft_archive.c and the names above cannot drift apart. */
#define CFT_ALIAS_PREFIX "cft_alias_"

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
    CFT_FD_S("accurate",         REB_INT,    accurate) \
    CFT_FD_S("E",                REB_SIZE_T, E) \
    CFT_FD_S("abi_0",            REB_UINT64, cft_abi[0]) \
    CFT_FD_S("abi_1",            REB_UINT64, cft_abi[8]) \
    CFT_FD_S("constants_digest", REB_UINT64, constants_digest) \
    CFT_FD_S("iterations_max_exceeded", REB_UINT64, iterations_max_exceeded) \
    CFT_FD_S("provenance",       REB_INT,    provenance) \
    CFT_FD_S("hiwater_n_elem",   REB_SIZE_T, hiwater_n_elem) \
    CFT_FD_S("n_elem",           REB_SIZE_T, n_elem)      /* LAST: see NOTE 3 */
/* cft_hiwater_n_elem sits among the scalars, i.e. AFTER the alias
 * blobs, for the reason cft_n_elem is after the others: each alias
 * blob's read leaves size_data/8 behind in hiwater_n_elem, and the
 * scalar that follows them overwrites it with the truth. Neither is
 * relied on - cft_archive_finish_load() repairs both from the file,
 * because an APPENDED snapshot is a diff and may omit an unchanged
 * scalar. NOTE 3 in src/cft_archive.c is the whole argument. */

#define CFT_FD_LIST(W) { CFT_FD_BLOBS(W) CFT_FD_ALIAS(W) CFT_FD_SCALARS { 0 } }

/* Blobs of n_elem elements - CFT_FD_BLOBS, and what
 * cft_archive_state_blob() walks. */
#define CFT_N_BLOBS 50
/* Blobs of hiwater_n_elem elements - CFT_FD_ALIAS: six families of
 * seven levels, plus csx and csv. A separate count because it is a
 * separate length; see the comment above CFT_FD_A. */
#define CFT_N_ALIAS_BLOBS 44
/* And the scalars CFT_FD_SCALARS adds after them. Here rather than
 * written out in cft_archive_selftest(), for the reason the blob
 * count is here: a number kept anywhere but beside the list it counts
 * is a number that will disagree with it. */
#define CFT_N_SCALARS 15

/* One definition of each, in src/cft_ias15_fields.c. cft_archive.c
 * selects among them by format; the integrator registers the first,
 * which is why a load resolves names at every format - see the note in
 * that file. */
extern const struct reb_binarydata_field_descriptor cft_fd_fp64[];
extern const struct reb_binarydata_field_descriptor cft_fd_fp128[];
extern const struct reb_binarydata_field_descriptor cft_fd_fp256[];

#endif /* CFT_IAS15_FIELDS_H */
