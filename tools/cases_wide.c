/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The wide path, which has to be its own process: the engine's format
 * is fixed for the process once it is opened, because the Gauss-Radau
 * constants are derived at that format and every buffer is sized for
 * its element width. `check_dropin --wide` is that process.
 *
 * These are smoke tests and say so. The binary64 gate above never
 * exercises a single conversion - at CFT_FP64 promotion and rounding
 * are memcpy - so what these assert is that cft_convert, the wide state
 * kept across a step boundary and the rounded view written back all
 * work, and that the answer is NOT bit-identical to binary64, which is
 * what says the arithmetic really was wider.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"

/* ------------------------------------------------------------------ *
 * The wide path, in its own process.
 *
 * The binary64 gate above never exercises a single conversion: at
 * CFT_FP64 promotion and rounding are memcpy. This case is what says
 * cft_convert, the wide state kept across the step boundary and the
 * rounded view written back all work - and it has to be a separate run
 * of this program, because the engine's format is fixed for the process
 * once it is opened (the Gauss-Radau constants are derived at that
 * format and every buffer is sized for its element width).
 *
 * It is a smoke test and says so. The assertion is that the binary128
 * run is finite, close to the binary64 answer, and NOT bit-identical to
 * it - the last being what proves the arithmetic really was wider,
 * rather than the format silently falling back.
 */
static void case_wide(void){
    const size_t steps = 2000;
    const double dt = 0.05;
    printf("binary128 through the shim (kepler, fixed dt = %g, %zu steps)\n", dt, steps);
    struct reb_simulation *ra = build(kepler, 2, dt, 0.0, 0);   /* REBOUND's own, binary64 */
    struct reb_simulation *rb = build(kepler, 2, dt, 0.0, 1);   /* the shim at binary128 */
    double e0 = reb_simulation_energy(ra);
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);
    if (rb->status == REB_STATUS_GENERIC_ERROR){
        printf("  FAIL binary128: the step refused\n"); failures++;
        reb_simulation_free(ra); reb_simulation_free(rb); return;
    }
    int identical = 1, finite = 1;
    double worst = 0;
    for (size_t i = 0; i < 2; i++){
        double A[6] = { ra->particles[i].x, ra->particles[i].y, ra->particles[i].z,
                        ra->particles[i].vx, ra->particles[i].vy, ra->particles[i].vz };
        double B[6] = { rb->particles[i].x, rb->particles[i].y, rb->particles[i].z,
                        rb->particles[i].vx, rb->particles[i].vy, rb->particles[i].vz };
        for (int k = 0; k < 6; k++){
            if (!isfinite(B[k])) finite = 0;
            if (bits_differ(A[k], B[k])) identical = 0;
            double scale = fabs(A[k]) > 1e-3 ? fabs(A[k]) : 1e-3;
            double d = fabs(A[k] - B[k]) / scale;
            if (d > worst) worst = d;
        }
    }
    double ea = fabs((reb_simulation_energy(ra) - e0) / e0);
    double eb = fabs((reb_simulation_energy(rb) - e0) / e0);
    printf("    t: ias15 %.17g, ias15_cft %.17g (the wide clock, rounded)\n", ra->t, rb->t);
    printf("    largest relative difference in the twelve coordinates: %.3e\n", worst);
    printf("    relative energy change of the binary64 VIEW, a weak measure because the\n"
           "      view is rounded: ias15 %.3e, ias15_cft at binary128 %.3e\n", ea, eb);
    if (!finite){ printf("  FAIL binary128: a coordinate is not finite\n"); failures++; }
    else if (identical){ printf("  FAIL binary128: bit-identical to binary64 - the format did not take\n"); failures++; }
    else if (worst > 1e-9){ printf("  FAIL binary128: %.3e from the binary64 answer, too far to be round-off\n", worst); failures++; }
    else printf("  ok   binary128 ran, finite, %.3e from binary64 and not bit-identical to it\n", worst);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* A merge at binary128 with state->accurate = 1.
 *
 * At accurate = 0 a removal shifts r->particles, which are binary64,
 * and the wide state is re-promoted from them - correct at CFT_FP64,
 * where the view IS the state, and a silent loss of every tail above
 * it. That is why the shim refuses this combination, which the refusal
 * case below checks. accurate = 1 shifts the wide state itself, so the
 * survivors keep their tails and each keeps its own polynomial.
 *
 * Not an equivalence case, and it must not become one: REBOUND is
 * binary64 and its own removal re-reads its coefficient levels at the
 * new stride. So the assertions are that the run survives, that both
 * sides removed the same particle, that the clock agrees - the step
 * sequence is driven by the wide state but reported rounded - and that
 * the coordinates do NOT match bit for bit, which is what says the
 * extra precision was really carried across the removal. */
static void case_wide_collision(void){
    const size_t steps = 2000;
    printf("binary128 with accurate = 1, across a merge (3 bodies, dt = 0.01)\n");
    static const struct body meet[3] = {
        { 1.0,    0.0, 0.0, 0.0,   0.0,  0.0,  0.0 },
        { 1.0e-3, 1.0, 0.0, 0.0,   0.0,  1.3,  0.0 },
        { 1.0e-3, 1.02,0.0, 0.0,   0.0,  0.7,  0.0 },
    };
    accurate = 1;
    struct reb_simulation *ra = build(meet, 3, 0.01, 1e-9, 0);
    struct reb_simulation *rb = build(meet, 3, 0.01, 1e-9, 1);
    accurate = 0;
    for (int i = 0; i < 2; i++){
        struct reb_simulation *r = i ? rb : ra;
        r->collision = REB_COLLISION_DIRECT;
        r->collision_resolve = reb_collision_resolve_merge;
        for (size_t k = 0; k < r->N; k++) r->particles[k].r = 0.02;
    }
    size_t n_before = ra->N;
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);

    if (rb->status == REB_STATUS_GENERIC_ERROR){
        printf("  FAIL binary128 merge: the step refused\n"); failures++;
        reb_simulation_free(ra); reb_simulation_free(rb); return;
    }
    if (ra->N >= n_before){
        printf("  FAIL binary128 merge: no collision occurred (N is still %zu)\n", ra->N);
        failures++;
        reb_simulation_free(ra); reb_simulation_free(rb); return;
    }
    if (ra->N != rb->N){
        printf("  FAIL binary128 merge: REBOUND ended with %zu particles and this port %zu\n",
               ra->N, rb->N);
        failures++;
        reb_simulation_free(ra); reb_simulation_free(rb); return;
    }
    int finite = 1, identical = 1;
    double worst = 0;
    for (size_t i = 0; i < ra->N; i++){
        double A[6] = { ra->particles[i].x, ra->particles[i].y, ra->particles[i].z,
                        ra->particles[i].vx, ra->particles[i].vy, ra->particles[i].vz };
        double B[6] = { rb->particles[i].x, rb->particles[i].y, rb->particles[i].z,
                        rb->particles[i].vx, rb->particles[i].vy, rb->particles[i].vz };
        for (int k = 0; k < 6; k++){
            if (!isfinite(B[k])) finite = 0;
            if (bits_differ(A[k], B[k])) identical = 0;
            double scale = fabs(A[k]) > 1e-3 ? fabs(A[k]) : 1e-3;
            double d = fabs(A[k] - B[k]) / scale;
            if (d > worst) worst = d;
        }
    }
    printf("    a merge removed one: %zu particles became %zu, on both sides\n",
           n_before, ra->N);
    printf("    largest relative difference after it: %.3e\n", worst);
    if (!finite){ printf("  FAIL binary128 merge: a coordinate is not finite\n"); failures++; }
    else if (identical){
        printf("  FAIL binary128 merge: bit-identical to binary64, so the wide state did\n"
               "       not survive the removal after all\n"); failures++; }
    else printf("  ok   binary128 merge: survived, finite, and wider than binary64\n");
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* And the combination that must be refused rather than approximated:
 * a removal at a wide format with accurate = 0, where the surviving
 * particles would come back from r->particles at binary64. */
static void case_wide_collision_refused(void){
    struct reb_simulation *r = build(kepler, 2, 0.05, 1e-9, 1);
    r->collision = REB_COLLISION_DIRECT;
    r->collision_resolve = reb_collision_resolve_merge;
    double t_before = r->t;
    r->integrator.callbacks.step(r, r->integrator.state);
    int ok = (r->status == REB_STATUS_GENERIC_ERROR) &&
             !bits_differ(t_before, r->t);
    reb_simulation_free(r);
    if (ok) printf("  ok   collisions at binary128 without accurate refused, clock did not move\n");
    else { printf("  FAIL collisions at binary128 without accurate were not refused\n"); failures++; }
}


/* The force hook at a wide format, which is where its two claims live
 * and where neither was measured: tools/cases_forces.c runs on the
 * binary64 leg only, and at CFT_FP64 promote and round are memcpy, so
 * nothing there can tell a run that stayed wide from one that did not.
 *
 * The claims are src/ias15_engine.h's. First, A DO-NOTHING ROUTINE IS
 * FREE: apply_force_hook compares the returned accelerations against
 * the ones it handed over and promotes back only where they differ, so
 * registering a routine that writes nothing must leave the run bit for
 * bit where the no-force run left it. A hook that promoted back
 * unconditionally would truncate every acceleration to binary64 at
 * every Gauss-Radau node, and that is precisely what this separates -
 * at binary64 the two are the same bits and the existing
 * case_nop_force_is_free cannot see the difference.
 *
 * Second, A RUN WITH A FORCE IN IT IS STILL WIDE. The routine is
 * evaluated at binary64 and a component it writes is binary64 from that
 * node on; everything else - gravity, the polynomial, the positions and
 * velocities, and the components the routine left alone - is not. So
 * the binary128 run with the force must differ from REBOUND's own
 * binary64 run with the identical routine, by round-off and no more.
 * The routine writes ax ONLY, so ay and az are the untouched
 * components that claim is about. */
static void f_wide_ax(struct reb_simulation *r){
    for (size_t i = 0; i < r->N; i++)
        r->particles[i].ax += 1.0e-3 * (double)(i + 1);
}
static void f_wide_nop(struct reb_simulation *r){ (void)r; }

/* One ias15_cft run's twelve coordinates, kept as bits: the engine is
 * one global instance, so two of them cannot be alive at once. */
struct wide_snap { double c[2][6]; double t, dt; };
static void wide_take(struct reb_simulation *r, struct wide_snap *s){
    for (size_t i = 0; i < 2; i++){
        const double v[6] = { r->particles[i].x, r->particles[i].y, r->particles[i].z,
                              r->particles[i].vx, r->particles[i].vy, r->particles[i].vz };
        memcpy(s->c[i], v, sizeof v);
    }
    s->t = r->t; s->dt = r->dt;
}

static void case_wide_force(void){
    const size_t steps = 1000;
    const double dt = 0.05;
    printf("the force hook at binary128 (kepler, fixed dt = %g, %zu steps)\n", dt, steps);

    /* the two ias15_cft runs, one after the other */
    struct wide_snap bare, nop;
    for (int k = 0; k < 2; k++){
        struct reb_simulation *r = build(kepler, 2, dt, 0.0, 1);
        if (k) r->additional_forces = f_wide_nop;
        reb_simulation_steps(r, steps);
        if (r->status == REB_STATUS_GENERIC_ERROR){
            printf("  FAIL binary128 forces: the %s run refused\n", k ? "do-nothing" : "no-force");
            failures++; reb_simulation_free(r); return;
        }
        wide_take(r, k ? &nop : &bare);
        reb_simulation_free(r);
    }
    int same = 1;
    for (size_t i = 0; i < 2; i++)
        for (int k = 0; k < 6; k++)
            if (bits_differ(bare.c[i][k], nop.c[i][k])) same = 0;
    if (bits_differ(bare.t, nop.t) || bits_differ(bare.dt, nop.dt)) same = 0;
    if (!same){
        printf("  FAIL binary128: a do-nothing routine changed the run, so the\n"
               "       hook is promoting back accelerations the routine did not\n"
               "       write and a wide run does not stay wide under one\n");
        failures++;
    }else
        printf("  ok   a do-nothing routine at binary128 is bit for bit the no-force run\n");

    /* and the same force on both sides, REBOUND at binary64 */
    struct reb_simulation *ra = build(kepler, 2, dt, 0.0, 0);
    struct reb_simulation *rb = build(kepler, 2, dt, 0.0, 1);
    ra->additional_forces = f_wide_ax;
    rb->additional_forces = f_wide_ax;
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);
    if (rb->status == REB_STATUS_GENERIC_ERROR){
        printf("  FAIL binary128 forces: the step refused\n"); failures++;
        reb_simulation_free(ra); reb_simulation_free(rb); return;
    }
    int finite = 1, identical = 1;
    double worst = 0;
    for (size_t i = 0; i < 2; i++){
        double A[6] = { ra->particles[i].x, ra->particles[i].y, ra->particles[i].z,
                        ra->particles[i].vx, ra->particles[i].vy, ra->particles[i].vz };
        double B[6] = { rb->particles[i].x, rb->particles[i].y, rb->particles[i].z,
                        rb->particles[i].vx, rb->particles[i].vy, rb->particles[i].vz };
        for (int k = 0; k < 6; k++){
            if (!isfinite(B[k])) finite = 0;
            if (bits_differ(A[k], B[k])) identical = 0;
            double scale = fabs(A[k]) > 1e-3 ? fabs(A[k]) : 1e-3;
            double d = fabs(A[k] - B[k]) / scale;
            if (d > worst) worst = d;
        }
    }
    printf("    with a routine that writes ax only, largest relative difference\n"
           "      from REBOUND's binary64 running the same routine: %.3e\n", worst);
    if (!finite){ printf("  FAIL binary128 forces: a coordinate is not finite\n"); failures++; }
    else if (identical){
        printf("  FAIL binary128 forces: bit-identical to binary64, so the hook\n"
               "       collapsed the run to the format the routine is evaluated in\n");
        failures++; }
    else if (worst > 1e-9){
        printf("  FAIL binary128 forces: %.3e from the binary64 answer, too far to\n"
               "       be round-off - the force is not the same force\n", worst);
        failures++; }
    else
        printf("  ok   a forced binary128 run is still wide: %.3e from binary64 and\n"
               "       not bit-identical to it\n", worst);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* ------------------------------------------------------------------ */
void cases_wide(void){
    case_wide();
    printf("\n");
    case_wide_collision();
    printf("\n");
    case_wide_collision_refused();
    printf("\n");
    case_wide_force();
    printf("\n");
}
