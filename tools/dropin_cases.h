/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The drop-in gate's scaffolding, declared once.
 *
 * tools/check_dropin.c was one file with every case in it, one main()
 * and one case_refusals(). Round 2 has five parcels and every one adds
 * cases, so that file was five agents editing three regions - the
 * collision P0 exists to remove (PARCELS.md).
 *
 * Now a topic is a file:
 *
 *     tools/dropin_common.c   this header's definitions
 *     tools/cases_core.c      the cases that were here first
 *     tools/cases_wide.c      the binary128 smoke tests
 *     tools/check_dropin.c    main(), the refusals, the coverage check
 *
 * ADDING A TOPIC. Write tools/cases_<topic>.c, include this header,
 * keep every case function `static`, and expose ONE entry point -
 * `void cases_<topic>(void);` - declared at the bottom of this file.
 * Add the object to CASES_OBJ in the Makefile and one call in main().
 * That is the whole of what you touch outside your own file.
 *
 * A case reports by printing and by incrementing `failures`. Nothing
 * returns a verdict: main() looks at `failures` at the end, so a case
 * that forgets to report a failure is a case that passes, and that is
 * what the negative controls in PARCELS.md are for.
 */
#ifndef DROPIN_CASES_H
#define DROPIN_CASES_H

#include <stddef.h>
#include "rebound.h"
#include "cft_ias15.h"

/* ---- what every case reports through ------------------------------ */
extern int verbose;         /* -v: print every value, not only the differing ones */
extern int failures;        /* the only verdict; main() reads it once at the end */

/* ---- the problems -------------------------------------------------- */
struct body { double m, x, y, z, vx, vy, vz; };

extern const struct body kepler[];        /* 2 bodies */
extern const struct body pythagorean[];   /* 3 */
extern const struct body five[];          /* 5, and the first 4 are a system too */

/* ---- the knobs build() applies ------------------------------------- *
 * File-scope rather than parameters because they are settings of the
 * RUN, applied to both sides, and threading five of them through every
 * case signature would make each case's own arguments hard to find.
 * Set, call, set back - the existing groups in cases_core.c show the
 * shape. Leave them as you found them. */
extern int    wide_format;    /* the cft side's format; --wide makes it CFT_FP128 */
extern double soften;         /* r->softening, on BOTH sides */
extern double min_dt;         /* IAS15's step floor, on both sides */
extern int    adaptive_mode;  /* 2 PRS23, 3 AARSETH85 */
extern int    accurate;       /* state->accurate; diverges from REBOUND, so it is
                               * off for every equivalence case */

/* Both simulations from the same doubles, so the initial conditions
 * cannot be the source of a difference. use_cft picks the integrator. */
struct reb_simulation *build(const struct body *bs, size_t n, double dt,
                             double epsilon, int use_cft);

/* ---- comparison: bit patterns, never values ------------------------ */
int  bits_differ(double a, double b);
void show(const char *what, size_t i, double a, double b);
/* Every particle by nine values, plus t, dt and dt_last_done. Prints and
 * counts; returns nonzero if anything differed. */
int  compare(struct reb_simulation *ra, struct reb_simulation *rb, const char *label);

/* ---- refusals ------------------------------------------------------ *
 * `row` is the name of the row in cft_support_rows this case exercises.
 * It is not decoration: check_dropin.c walks the table afterwards and
 * fails on any drop-in row no case named, so a capability cannot be
 * added to the table without a test, or have its test quietly retired.
 * Pass the row name exactly as src/cft_supported.c spells it. */
int refused(const char *row, const char *what, void (*poison)(struct reb_simulation *));
int accepted(const char *row, const char *what, void (*prepare)(struct reb_simulation *));

/* ---- the topics ---------------------------------------------------- */
void cases_core(void);   /* tools/cases_core.c - the binary64 equivalence set */
void cases_wide(void);   /* tools/cases_wide.c - binary128, a smoke test */
void cases_forces(void); /* tools/cases_forces.c - r->additional_forces at every
                          * substage, and the timestep-modification hooks */

#endif /* DROPIN_CASES_H */
