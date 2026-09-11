/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * What this port does not do, stated once.
 *
 * It used to be stated four times: supported() in
 * src/reb_integrator_cft.c, cft_rebound_check() in
 * src/cft_rebound_run.c, case_refusals() in tools/check_dropin.c, and
 * the support table in README.md. Round 1's post-mortem (ROADMAP.md,
 * the last section) is the record of what that costs - "where a shared
 * fact appeared twice, the copies drifted", and the blob count drifted
 * in four places at once. The refusal list is the same shape of fact
 * with more copies, and round 2 has five parcels each removing one
 * refusal, so without this every parcel would edit every copy.
 *
 * Now: one row per capability. The two entry points walk the same
 * table through cft_support_first_refusal(), and the gate walks it too,
 * failing by name on any row it did not exercise - so a capability
 * cannot be added without a test, or have its test retired without
 * being noticed.
 *
 * The two paths still refuse DIFFERENT sets, deliberately, and the
 * README says which and why. That difference now lives in one place -
 * each row's `paths` mask - instead of being the emergent result of two
 * functions written months apart.
 *
 * Adding a capability: add a row, add its case to the gate. Removing a
 * refusal: delete the row, and the gate tells you its case is now
 * orphaned.
 */
#ifndef CFT_SUPPORTED_H
#define CFT_SUPPORTED_H

#include <stddef.h>
#include "rebound.h"
#include "ias15_limits.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The subprocess API's spelling of src/ias15_limits.h's cap, kept
 * because callers may have written it. One number, two names, no
 * second value to drift. */
#define CFT_REBOUND_MAX_BODIES CFT_MAX_BODIES

/* Which entry point a row applies to. They are not the same list:
 * the drop-in supports collisions and attached ODE sets that the
 * subprocess API refuses, because REBOUND's own driver runs between its
 * steps and does not run between the subprocess API's. */
#define CFT_PATH_DROPIN     0x1u   /* the registered integrator */
#define CFT_PATH_SUBPROCESS 0x2u   /* cft_rebound_check() and the API over it */
#define CFT_PATH_BOTH       (CFT_PATH_DROPIN | CFT_PATH_SUBPROCESS)

/* Everything a row is allowed to read.
 *
 * Deliberately not "the simulation plus whichever state struct this
 * path happens to have": the drop-in reads struct cft_ias15_state and
 * the subprocess path reads REBOUND's own struct
 * reb_integrator_ias15_state, and a row that had to know which would be
 * the same fact stated twice again. Each caller normalises into this,
 * and every row reads one uniform thing.
 */
struct cft_support_ctx {
    const struct reb_simulation *r;

    int    have_state;     /* 1 on the drop-in path, where the three
                            * fields below are this port's own settings.
                            * 0 on the subprocess path, which has no
                            * cft_ias15_state to read - the run is
                            * handed to the standalone program. */
    int    format;         /* CFT_FP64 | CFT_FP128 | CFT_FP256 */
    int    accurate;       /* state->accurate */
    int    max_iter;       /* after the 0 -> default fixup */
    size_t E;              /* state->E; 1 for a reb_simulation */

    int    adaptive_mode;  /* the criterion the run will actually use.
                            * 2 (PRS23) when a caller cannot tell, which
                            * is a supported value and so refuses
                            * nothing - the subprocess path can only read
                            * it when the simulation's integrator really
                            * is REBOUND's "ias15". */
};

struct cft_support_row {
    const char *name;      /* stable identifier: the gate's key, and the
                            * word to grep for. Not user-visible. */
    unsigned    paths;     /* CFT_PATH_* */
    int       (*hit)(const struct cft_support_ctx *c);   /* nonzero: refuse */
    void      (*say)(const struct cft_support_ctx *c, char *buf, size_t n);
};

/* Terminated by a row whose name is NULL. Public so the gate can walk
 * it for coverage. */
extern const struct cft_support_row cft_support_rows[];

/* The first row that applies to `path` and refuses `c`, or NULL if
 * nothing does. Both entry points are this call plus their own way of
 * reporting; the loop exists once.
 *
 * A row's message carries no prefix - the caller adds its own, because
 * "ias15_cft: " belongs to the integrator and not to the fact. */
const struct cft_support_row *
cft_support_first_refusal(const struct cft_support_ctx *c, unsigned path);

#ifdef __cplusplus
}
#endif
#endif /* CFT_SUPPORTED_H */
