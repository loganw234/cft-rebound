/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * cft_ias15_state.h - the state struct that the REBOUND integrator shim
 * (parcel A) and the Simulationarchive support (parcel B) share.
 *
 * The definition below is ROADMAP.md's, verbatim, and is owned by the
 * integrator shim: this file exists so the archive can be written and
 * gated in parallel with it. When the parcels are integrated there must
 * be exactly one copy of this header; if the shim's copy differs from
 * this one in any member name, type or order, the archive's field
 * descriptors are wrong and the gates will say so loudly rather than
 * quietly writing the wrong bytes (see CFT_IAS15_STATE_ABI below).
 *
 * W is the wide element width in bytes: 8, 16 or 32, i.e.
 * cft_format_size(state->format).
 */
#ifndef CFT_IAS15_STATE_H
#define CFT_IAS15_STATE_H

#include <stddef.h>
#include <stdint.h>
#include "cft.h"

/* Bumped whenever a member is added, removed, renamed or retyped.
 * The archive's descriptor lists are built against this value. */
#define CFT_IAS15_STATE_ABI 1

struct cft_ias15_state {
    /* configuration, plain doubles, archived as REB_DOUBLE */
    double   epsilon;          /* as REBOUND's */
    double   min_dt;
    int      adaptive_mode;
    int      format;           /* CFT_FP64 | CFT_FP128 | CFT_FP256 */
    int      max_iter;         /* 12 in REBOUND; binary256 needs ~22 */
    int      arith_fma;        /* 0 = REBOUND's roundings, 1 = the FMA form */

    /* the wide state: byte blobs, archived as REB_POINTER with
     * element_size = W. Lengths are 3N or 3N*E elements. */
    unsigned char *x0, *v0, *a0;      /* position, velocity, acceleration */
    unsigned char *csx, *csv, *csa0;  /* compensated-summation carries */
    unsigned char *g[7], *b[7], *e[7], *br[7], *er[7], *csb[7];
    size_t   n_elem;           /* 3N, or 3NE for an ensemble */
    size_t   E;                /* ensemble members, 1 for a plain sim */

    /* provenance, so a reader knows what produced the wide bytes */
    char     cft_abi[16];      /* e.g. "0.11" */
    uint64_t constants_digest; /* the Gauss-Radau set the run used */
};

/* Bytes per element for this state's format. 0 if format is not one of
 * libcft's four. */
static inline size_t cft_ias15_state_width(const struct cft_ias15_state *s){
    return cft_format_size((cft_format)s->format);
}

#endif /* CFT_IAS15_STATE_H */
