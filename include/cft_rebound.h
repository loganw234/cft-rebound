/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * cft-rebound: REBOUND's IAS15 at binary64, binary128 or binary256.
 *
 * This is the whole public surface. A REBOUND program adds two lines to
 * its build,
 *
 *     CFLAGS  += -I$(CFT_REBOUND_PREFIX)/include
 *     LDFLAGS += -L$(CFT_REBOUND_PREFIX)/lib -lcft_rebound -lcft
 *
 * includes this header instead of (or as well as) rebound.h, and calls
 * cft_rebound_steps() where it would have called reb_simulation_steps().
 *
 *
 * WHICH FORM THIS IS
 * ------------------
 * There are two ways to reach the port, and this header describes the
 * one that exists today.
 *
 *   The subprocess form (today). cft_rebound_steps() writes the
 *   simulation's particles out as an exact binary64 problem file, runs
 *   the `ias15_cft` program on it at the format you asked for, and
 *   rounds the wide result back into r->particles. Cost is one process
 *   per call, which is nothing next to the arithmetic.
 *
 *   The wide state does not survive the call, and that is not merely an
 *   efficiency matter: IAS15 starts each step's corrector from the
 *   PREVIOUS step's b coefficients, so a call that begins with them
 *   zeroed takes a different iteration path and lands on different
 *   bits. Measured on Kepler at binary64, twenty fixed steps: one call
 *   of twenty and two calls of ten differ by one to six ulps in the
 *   final state (docs/VALIDATION.md). **Ask for a whole run in one
 *   call.** The registration form below removes the restriction,
 *   because there the state persists.
 *
 *   The registration form (ROADMAP parcel A). REBOUND accepts
 *   user-provided integrators - reb_integrator_register() plus
 *   reb_simulation_set_integrator() - so the port becomes a struct
 *   reb_integrator whose wide state lives across steps, and
 *   cft_rebound_steps() becomes a three-line wrapper around
 *   reb_simulation_steps(). The signature below is chosen so that
 *   swapping the inside changes nothing a caller wrote.
 *
 * Either way the division is the one ROADMAP.md sets out: REBOUND's
 * struct reb_particle holds double, so r->particles is the binary64
 * VIEW of the run, correctly rounded from the wide state, and every
 * REBOUND output, callback and tool keeps working on it. The wide
 * state is the integrator's, not REBOUND's.
 */
#ifndef CFT_REBOUND_H
#define CFT_REBOUND_H

#include <stddef.h>
#include "rebound.h"

#ifdef __cplusplus
extern "C" {
#endif

/* The format the integration is carried in. The particles you get back
 * are binary64 whichever you choose. */
enum cft_rebound_format {
    CFT_REBOUND_FP64  = 0,   /* bit for bit REBOUND's own IAS15 */
    CFT_REBOUND_FP128 = 1,   /* the one that is usually worth having */
    CFT_REBOUND_FP256 = 2
};

/* Zero-initialise this and set what you care about.
 *
 * `epsilon` keeps REBOUND's own meaning exactly, including that zero
 * means a fixed step of r->dt - so a zero-initialised struct asks for a
 * fixed step, not for REBOUND's default. Set epsilon = 1e-9 for that.
 * The step criterion is PRS23, REBOUND's 2024 default and the only one
 * the port implements; an ias15 state asking for another is refused. */
struct cft_rebound_options {
    enum cft_rebound_format format;
    double epsilon;        /* IAS15's tolerance: > 0 adaptive, 0 fixed */
    int    max_iter;       /* the cap on IAS15's corrector passes. 0 picks
                            * a default per format - see the note in
                            * src/cft_rebound_run.c - because the program's
                            * own default is REBOUND's 12 at every format
                            * and binary256 needs more than that. If even
                            * the cap is reached, the run says so in
                            * result.iterations_max_exceeded. */
    const char *program;   /* the ias15_cft executable. NULL searches, in
                            * order: $CFT_REBOUND_IAS15, then the install
                            * path fixed at build time if that file is
                            * there, then "ias15_cft" on PATH. */
    const char *workdir;   /* where the two temporary files go. NULL uses
                            * $TMPDIR / $TEMP / $TMP, then "." */
    int    verbose;        /* 1: echo the command line to stderr */
};

/* What the wide run reports about itself. Every double here has been
 * correctly rounded from the run's own format. */
struct cft_rebound_result {
    double   energy;        /* the total energy AS THE WIDE FORMAT COMPUTED IT.
                             * Not the same number as reb_simulation_energy(r)
                             * afterwards, which recomputes it at binary64
                             * from the rounded view - the difference between
                             * the two is the view's rounding, not an error. */
    char     energy_hex[96];/* and the same number UNROUNDED, as exact
                             * hexadecimal text in the run's own format.
                             * At binary128 and above the drift of a good
                             * run is below binary64's last bit, so
                             * `energy` alone will say a run drifted by
                             * exactly nothing; these digits are where the
                             * difference actually is. */
    double   t;             /* the run's own running sum of the steps */
    double   t_hi, t_lo;    /* and the exact elapsed time as an unevaluated
                             * sum, t_hi + t_lo, for long runs */
    double   dt_next;       /* the step the run would take next */
    double   dt_last;       /* the last step it took */
    long     steps_done;
    long     steps_rejected;
    long     iterations_max_exceeded;
    double   mean_pc_iterations;
    int      max_pc_iterations;
    double   seconds;
};

/* Is this simulation one the port can integrate?
 *
 * Returns 0 if it is. Otherwise returns non-zero and, if `why` is not
 * NULL, writes a one-line reason into it (NUL-terminated, truncated to
 * `n`). The force model here is basic pairwise Newtonian gravity and
 * nothing else; README "Scope" lists the refusals and this function is
 * what makes that list true rather than a promise.
 *
 * cft_rebound_steps() calls this first and fails the same way, so you
 * only need it directly if you want to check before setting up. */
int cft_rebound_check(const struct reb_simulation *r, char *why, size_t n);

/* Integrate `nsteps` steps of `r` through the port.
 *
 * The direct analogue of reb_simulation_steps(r, nsteps), and it leaves
 * r in the same state that would: r->particles, r->t, r->dt and
 * r->steps_done updated, everything else untouched. `opt` may be NULL
 * for all defaults; `res` may be NULL if you do not want the report.
 *
 * Returns 0 on success. On failure returns non-zero, leaves `r`
 * untouched, and prints one line to stderr saying what went wrong. */
int cft_rebound_steps(struct reb_simulation *r, long nsteps,
                      const struct cft_rebound_options *opt,
                      struct cft_rebound_result *res);

/* The name of a format, for messages: "fp64", "fp128", "fp256". */
const char *cft_rebound_format_name(enum cft_rebound_format f);

/* The ias15_cft executable this build will run, resolved by the rule in
 * struct cft_rebound_options. Never NULL. */
const char *cft_rebound_program_path(const struct cft_rebound_options *opt);

/* ------------------------------------------------------------------
 * Reserved for ROADMAP parcel A, the registration form:
 *
 *     extern const struct reb_integrator cft_ias15_integrator;
 *     void cft_rebound_register(void);
 *
 * When those land they belong in this header, beside the calls above,
 * and cft_rebound_steps() becomes a wrapper around them. Nothing
 * declared above needs to change for that to happen.
 * ------------------------------------------------------------------ */

#ifdef __cplusplus
}
#endif

#endif /* CFT_REBOUND_H */
