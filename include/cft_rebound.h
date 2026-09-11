/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * cft-rebound: REBOUND's IAS15 at binary64, binary128 or binary256.
 *
 * This is the whole public surface. A REBOUND program adds two lines to
 * its build,
 *
 *     CFLAGS  += -I$(PREFIX)/include
 *     LDLIBS  += -L$(PREFIX)/lib -lcft_rebound -lcft
 *
 * where PREFIX is what `make install` was given (the default is
 * /usr/local). examples/Makefile is those two lines in a working
 * makefile.
 *
 * includes this header instead of (or as well as) rebound.h, and calls
 * cft_rebound_steps() where it would have called reb_simulation_steps().
 *
 *
 * WHICH FORM THIS IS
 * ------------------
 * There are two ways to reach the port. Both exist; this header is the
 * subprocess one.
 *
 *   The subprocess form (this header). cft_rebound_steps() writes the
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
 *   call.** The registration form does not have this restriction,
 *   because there the state persists.
 *
 *   The registration form (src/cft_ias15.h). It LANDED - ROADMAP
 *   parcel A - and it is a different header and a different library,
 *   not a change to this one. REBOUND accepts user-provided
 *   integrators, so the port is a struct reb_integrator whose wide
 *   state lives across steps:
 *
 *       #include "cft_ias15.h"
 *       cft_ias15_register("ias15_cft");
 *       struct cft_ias15_state *s =
 *           reb_simulation_set_integrator(r, "ias15_cft");
 *       s->format = CFT_FP256;
 *       reb_simulation_steps(r, 1000);      // REBOUND's own call
 *
 *   Link -lcft_ias15 for that and -lcft_rebound for this; `make
 *   install` ships both, and examples/dropin.c is the worked example.
 *   cft_rebound_steps() was NOT rewritten as a wrapper around it: this
 *   file is still a subprocess, deliberately, for a caller that does
 *   not want the engine in its own address space.
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
 *
 * The step criterion is NOT here, deliberately: all four of REBOUND's
 * are implemented, and the one that runs is the SIMULATION's -
 * ((struct reb_integrator_ias15_state *)r->integrator.state)
 * ->adaptive_mode, when r->integrator.name really is "ias15". This call
 * forwards it to the program unconditionally. Only a number that names
 * none of the four is refused. It comes from there rather than from
 * here because it is REBOUND's own ias15 setting and a caller who has
 * one has already spelled it REBOUND's way.
 *
 * r->softening is forwarded the same way, from the simulation, and a
 * non-zero value is honoured rather than refused. */
struct cft_rebound_options {
    enum cft_rebound_format format;
    double epsilon;        /* IAS15's tolerance: > 0 adaptive, 0 fixed */
    double min_dt;         /* IAS15's floor on |dt|. 0 disables it, which
                            * is REBOUND's default: at 0 the comparison
                            * REBOUND makes is false and no floor is ever
                            * selected, so the arithmetic is unchanged.
                            * Alongside epsilon and max_iter here rather
                            * than read from the simulation, unlike
                            * adaptive_mode and r->softening - see the
                            * note above this struct for which is which. */
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
    const char *artifact;  /* the .xclbin to run on. NULL consults
                            * $CFT_REBOUND_ARTIFACT, and if that is unset
                            * or empty the software backend is used. The
                            * registered integrator resolves its artifact
                            * by the same rule, so a program that sets the
                            * variable gets the card on both paths. */
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
 * The registration form is NOT declared here.
 *
 * This header once reserved the names `cft_rebound_register()` and
 * `cft_ias15_integrator` for ROADMAP parcel A. Parcel A landed with its
 * own header instead, because the drop-in needs rebound.h's integrator
 * types and this header's callers may not want them:
 *
 *     src/cft_ias15.h    cft_ias15_register(), struct cft_ias15_state,
 *                        cft_ias15_configure(), cft_ias15_format_code()
 *     src/cft_archive.h  the Simulationarchive round trip
 *     -lcft_ias15        the library, beside -lcft_rebound
 *
 * There is no cft_rebound_register(); do not look for one.
 * ------------------------------------------------------------------ */

#ifdef __cplusplus
}
#endif

#endif /* CFT_REBOUND_H */
