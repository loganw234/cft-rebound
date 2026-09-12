/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * SEAM CASES: combinations no single parcel could have tested.
 *
 * Round 2 landed five parcels in parallel. Each gated its own work and
 * each was right, and the thing that goes untested in that arrangement
 * is the INTERSECTION - a configuration that needs two parcels' features
 * at once, which neither one's brief covered and neither one's gate
 * reaches. docs/PARCEL-ROUNDS.md calls these the lead's, because they
 * belong to no parcel; this file is where they go.
 *
 * The first one is the one the round actually produced.
 *
 * P1 ended its report with a paragraph headed "for P3": REBOUND's
 * GLOBAL step controller reads particles[].v AFTER the node loop
 * (integrator_ias15.c:632), and under a velocity-dependent
 * r->additional_forces that array holds the h[7]-PREDICTED velocities
 * (:446-449), which nothing restores. That note reached P3 hours after
 * P3 had started, which is why this repository now has a ledger.
 *
 * The port keeps its predicted velocities in FVP and leaves the wide v
 * at the step start, so the two build v2 from different numbers. A
 * documentation audit measured that the divergence is LATENT rather
 * than demonstrated - v2 feeds only the skip predicate
 * fabs(v2*dt*dt/x2) < 1e-16 and never the error estimate, and 260
 * configurations tuned to sit on that threshold all came back
 * bit-identical. What was missing was any gate on the combination at
 * all: tools/cases_forces.c pins adaptive_mode = 2, and
 * tools/cases_modes.c registers no force.
 *
 * So this asserts the usual thing - whatever REBOUND does, the port
 * does, bit for bit - for every step criterion crossed with a
 * velocity-dependent force. If the latency ever stops being latent,
 * this is what says so.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"

/* A velocity-dependent force: linear drag, large enough to matter and
 * small enough not to dominate. Both sides run this identical routine,
 * so it cannot itself be the source of a difference. */
static void drag(struct reb_simulation *r){
    const double k = 1e-3;
    for (size_t i = 0; i < r->N; i++){
        r->particles[i].ax -= k * r->particles[i].vx;
        r->particles[i].ay -= k * r->particles[i].vy;
        r->particles[i].az -= k * r->particles[i].vz;
    }
}

/* The same shape with no velocity in it, for the control that says the
 * VELOCITY dependence is what reaches the controller rather than merely
 * the presence of a force. */
static void push(struct reb_simulation *r){
    const double k = 1e-6;
    for (size_t i = 0; i < r->N; i++) r->particles[i].ax -= k;
}

static void attach(struct reb_simulation *r,
                   void (*fn)(struct reb_simulation *), int veldep){
    r->additional_forces = fn;
    r->force_is_velocity_dependent = veldep;
}

/* One criterion crossed with one force, against REBOUND's own ias15. */
static void case_mode_force(const char *label, const struct body *bs, size_t n,
                            double dt, double epsilon, size_t steps,
                            int mode, void (*fn)(struct reb_simulation *),
                            int veldep){
    static const char *const NM[4] = { "INDIVIDUAL", "GLOBAL", "PRS23", "AARSETH85" };
    struct reb_simulation *ra, *rb;

    adaptive_mode = mode;
    ra = build(bs, n, dt, epsilon, 0);
    rb = build(bs, n, dt, epsilon, 1);
    adaptive_mode = 2;

    attach(ra, fn, veldep);
    attach(rb, fn, veldep);
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);

    printf("%s (%s, %s force, N = %zu, %zu steps)\n",
           label, NM[mode], veldep ? "velocity-dependent" : "velocity-independent",
           n, light_ran(steps));
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* The controls. A gate that cannot fail is not a gate, and this one has
 * two ways to be vacuous: the force might not reach the run, and the
 * criterion might not either. Both are checked by running the port
 * against ITSELF with one thing changed - if either pair is
 * bit-identical, the case above is comparing two runs that do not
 * differ in the way it claims to test. */
/* One configuration, run to completion, its final coordinates kept.
 *
 * SEQUENTIALLY, and that is not a style choice: the engine is a single
 * global instance - every buffer in src/ias15_cft.c is a file-scope
 * static - so a second simulation on ias15_cft while the first is alive
 * is refused by name. Build, step, copy out, free, then the next. */
static void run_one(const struct body *bs, size_t n, double dt, double epsilon,
                    size_t steps, int mode,
                    void (*fn)(struct reb_simulation *), int veldep,
                    double *out, size_t *out_n){
    struct reb_simulation *r;
    adaptive_mode = mode;
    r = build(bs, n, dt, epsilon, 1);
    adaptive_mode = 2;
    attach(r, fn, veldep);
    reb_simulation_steps(r, steps);
    *out_n = r->N;
    for (size_t i = 0; i < r->N; i++){
        const struct reb_particle *p = &r->particles[i];
        out[6*i+0] = p->x;  out[6*i+1] = p->y;  out[6*i+2] = p->z;
        out[6*i+3] = p->vx; out[6*i+4] = p->vy; out[6*i+5] = p->vz;
    }
    reb_simulation_free(r);
}

static void control_differs(const char *what, const struct body *bs, size_t n,
                            double dt, double epsilon, size_t steps,
                            int mode_a, void (*fn_a)(struct reb_simulation *), int vd_a,
                            int mode_b, void (*fn_b)(struct reb_simulation *), int vd_b){
    double va[6 * 8], vb[6 * 8];
    size_t na = 0, nb = 0;
    int differ = 0;

    if (n > 8){ printf("  FAIL %s: control buffer holds 8 particles\n", what); failures++; return; }

    run_one(bs, n, dt, epsilon, steps, mode_a, fn_a, vd_a, va, &na);
    run_one(bs, n, dt, epsilon, steps, mode_b, fn_b, vd_b, vb, &nb);

    if (na != nb) differ = 1;
    else for (size_t k = 0; k < 6 * na && !differ; k++)
        if (bits_differ(va[k], vb[k])) differ = 1;

    if (differ) printf("  ok   %s: the two runs differ, so the case above is testing it\n", what);
    else{
        printf("  FAIL %s: bit-identical, so the case above compares two runs that\n"
               "       do not differ in the way it claims to test\n", what);
        failures++;
    }
}

/* ------------------------------------------------------------------ */
void cases_seam(void){
    const size_t N = 300;
    printf("\nseam cases - combinations no single parcel's gate reaches:\n");

    printf("\na step criterion crossed with a velocity-dependent force\n"
           "(REBOUND's GLOBAL controller reads the PREDICTED velocities; the\n"
           " port reads the step-start ones, and nothing gated the pair):\n");
    for (int mode = 0; mode < 4; mode++)
        case_mode_force("kepler, drag", kepler, 2, 0.05, 1e-9, N, mode, drag, 1);
    case_mode_force("pythagorean, drag", pythagorean, 3, 0.01, 1e-9, N, 1, drag, 1);
    case_mode_force("five bodies, drag", five, 5, 0.5, 1e-9, 200, 1, drag, 1);

    printf("\n  and the controls, so the case above cannot pass vacuously:\n");
    control_differs("GLOBAL with a drag against GLOBAL with none",
                    kepler, 2, 0.05, 1e-9, N,
                    1, drag, 1,
                    1, NULL, 0);
    control_differs("GLOBAL with a drag against PRS23 with the same drag",
                    kepler, 2, 0.05, 1e-9, N,
                    1, drag, 1,
                    2, drag, 1);
    control_differs("GLOBAL, drag against a velocity-INdependent push",
                    kepler, 2, 0.05, 1e-9, N,
                    1, drag, 1,
                    1, push, 0);
}
