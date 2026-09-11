/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The cases that were in tools/check_dropin.c before it was split: the
 * binary64 equivalence set. Every one of these asserts the same thing -
 * REBOUND's own ias15 and the registered ias15_cft agree BIT FOR BIT.
 *
 * cases_core() at the bottom is the entry point and the running order.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"

/* ------------------------------------------------------------------ */
/* Cases                                                                */
/* ------------------------------------------------------------------ */
static void case_steps(const char *label, const struct body *bs, size_t n,
                       double dt, double epsilon, size_t steps){
    printf("%s (N = %zu, dt = %g, epsilon = %g, %zu steps)\n", label, n, dt, epsilon, steps);
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);
    if (verbose) printf("    t = %.17g, dt = %.17g\n", ra->t, ra->dt);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}


/* A collision that really happens, at binary64.
 *
 * REBOUND's driver does the search and the resolution after the step
 * callback returns, so what is under test is whether this integrator
 * survives the removal the same way REBOUND's does - which, since
 * accurate = 0, means leaving the polynomial exactly as stale as
 * REBOUND leaves it. */
static void case_collision(const char *label, double dt, double epsilon, size_t steps){
    printf("%s (3 bodies, dt = %g, epsilon = %g, %zu steps, merge on contact)\n",
           label, dt, epsilon, steps);
    /* A star and two small bodies that meet, so that a merge still
     * leaves a system with a pair in it.
     *
     * Both tangential velocities straddle the circular speed at
     * r = 1 rather than zero, and that is the whole design. With
     * +0.3 and -0.3 the remnant carries almost no angular momentum,
     * falls into the star, and the run ends as one body - where
     * IAS15 has no pairs, divides 0 by 0 in the step controller, and
     * the two sides differ only in the SIGN of the quiet NaN, which
     * 754 does not specify and which says nothing about a removal.
     * With 1.3 and 0.7 the remnant stays in orbit: measured, the
     * merge lands at step 1291 and the two survivors are 1.0079
     * apart at step 2000.
     *
     * Late is also the point. Bodies that overlap at t = 0 merge on
     * step 1, where the polynomial being re-read is still all zeros
     * and re-reading it proves nothing. */
    static const struct body meet[3] = {
        { 1.0,    0.0, 0.0, 0.0,   0.0,  0.0,  0.0 },
        { 1.0e-3, 1.0, 0.0, 0.0,   0.0,  1.3,  0.0 },
        { 1.0e-3, 1.02,0.0, 0.0,   0.0,  0.7,  0.0 },
    };
    struct reb_simulation *ra = build(meet, 3, dt, epsilon, 0);
    struct reb_simulation *rb = build(meet, 3, dt, epsilon, 1);
    for (int i = 0; i < 2; i++){
        struct reb_simulation *r = i ? rb : ra;
        r->collision = REB_COLLISION_DIRECT;
        r->collision_resolve = reb_collision_resolve_merge;
        for (size_t k = 0; k < r->N; k++) r->particles[k].r = 0.02;
    }
    size_t n_before = ra->N;
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);

    /* Did one actually happen? Without this the case passes whether or
     * not the collision search ever fired, which would make it a test
     * of nothing. */
    if (ra->N >= n_before){
        printf("  FAIL %s: no collision occurred (N is still %zu), so this case\n"
               "       proves nothing; the bodies need to actually meet\n", label, ra->N);
        failures++;
        reb_simulation_free(ra); reb_simulation_free(rb);
        return;
    }
    if (ra->N != rb->N){
        printf("  FAIL %s: REBOUND ended with %zu particles and this port with %zu\n",
               label, ra->N, rb->N);
        failures++;
        reb_simulation_free(ra); reb_simulation_free(rb);
        return;
    }
    /* And did it leave a system behind? One body has no pairs, so
     * IAS15's controller divides 0 by 0 and the case would go on to
     * compare two NaNs whose sign 754 leaves to the implementation. */
    if (ra->N < 2){
        printf("  FAIL %s: the run collapsed to %zu particle(s); a one-body\n"
               "       IAS15 run has no pairs and compares 0/0 against 0/0\n",
               label, ra->N);
        failures++;
        reb_simulation_free(ra); reb_simulation_free(rb);
        return;
    }
    printf("  (a merge removed one: %zu particles became %zu)\n", n_before, ra->N);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

static void case_integrate(const char *label, const struct body *bs, size_t n,
                           double dt, double epsilon, double tmax){
    printf("%s (N = %zu, dt = %g, epsilon = %g, integrate to t = %g, exact_finish_time = 1)\n",
           label, n, dt, epsilon, tmax);
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    reb_simulation_integrate(ra, tmax);
    reb_simulation_integrate(rb, tmax);
    if (verbose) printf("    t = %.17g, dt = %.17g, dt_last_done = %.17g\n", ra->t, ra->dt, ra->dt_last_done);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* Two integrate() calls in a row: the second resets r->dt_last_done to
 * 0, which is what tells step_try not to predict e and b if its first
 * attempt is rejected. The shim has to carry that reset into the wide
 * state or the two diverge at the first rejection of the second call. */
static void case_two_calls(const char *label, const struct body *bs, size_t n,
                           double dt, double epsilon, double t1, double t2){
    printf("%s (two integrate() calls, to %g then %g)\n", label, t1, t2);
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    reb_simulation_integrate(ra, t1); reb_simulation_integrate(ra, t2);
    reb_simulation_integrate(rb, t1); reb_simulation_integrate(rb, t2);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* A user editing a coordinate between steps. The shim keeps the wide
 * state as the truth and only re-reads r->particles where they differ
 * from the view it wrote; this is the case that says the "differ" branch
 * works and is REBOUND's own behaviour. */
static void case_user_edit(const char *label, const struct body *bs, size_t n,
                           double dt, double epsilon){
    printf("%s (200 steps, a coordinate nudged, 200 more)\n", label);
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    reb_simulation_steps(ra, 200); reb_simulation_steps(rb, 200);
    ra->particles[1].vy += 1e-6; rb->particles[1].vy += 1e-6;
    reb_simulation_steps(ra, 200); reb_simulation_steps(rb, 200);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* A particle added mid-run, past the high-water mark. REBOUND's
 * ias15_alloc reallocates and realloc_dp7 zeroes the whole array; the
 * shim's did_add_particle invalidates and the next step does the same.
 * An equivalence case, like every other one here - including, since
 * the port performs REBOUND's re-slice rather than approximating it,
 * the two below that change the count without passing the mark. */
static void case_add_particle(const char *label, double dt, double epsilon){
    printf("%s (150 steps of 2 bodies, a third added, 150 more)\n", label);
    struct reb_simulation *ra = build(kepler, 2, dt, epsilon, 0);
    struct reb_simulation *rb = build(kepler, 2, dt, epsilon, 1);
    reb_simulation_steps(ra, 150); reb_simulation_steps(rb, 150);
    struct reb_particle p = {0};
    p.m = 1e-4; p.x = -2.5; p.y = 0.3; p.z = 0.05; p.vy = -0.62; p.vz = 0.01;
    reb_simulation_add(ra, p);
    reb_simulation_add(rb, p);
    reb_simulation_steps(ra, 150); reb_simulation_steps(rb, 150);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}


/* Down and back up again, both times below the high-water mark.
 *
 * REBOUND reallocates nothing in either direction, so both are a
 * re-reading: its seven coefficient levels share one buffer that
 * dpcast() slices at the CURRENT 3N. Going down leaves a tail above the
 * live region that no stride reaches. Coming back up reads that tail
 * straight back, at offsets the shrunk configuration never wrote.
 *
 * That second half is the point of this case. Nothing else in the gate
 * reaches it, and a port that re-sliced on a shrink but did not keep
 * what the shrink stranded would pass every other case and read zeros
 * here. Four bodies so that the mark is 4, a removal leaves 3, and the
 * one added back lands under it. */
static void case_remove_then_add(const char *label, double dt, double epsilon){
    printf("%s (4 bodies, 200 steps, one removed, 200, one added back under the mark, 200)\n",
           label);
    struct reb_simulation *ra = build(five, 4, dt, epsilon, 0);
    struct reb_simulation *rb = build(five, 4, dt, epsilon, 1);
    reb_simulation_steps(ra, 200); reb_simulation_steps(rb, 200);
    if (reb_simulation_remove_particle(ra, 2) != 0 ||   /* 0 is success */
        reb_simulation_remove_particle(rb, 2) != 0){
        printf("  FAIL %s: REBOUND refused the removal\n", label); failures++;
        reb_simulation_free(ra); reb_simulation_free(rb); return;
    }
    reb_simulation_steps(ra, 200); reb_simulation_steps(rb, 200);
    struct reb_particle p = {0};
    p.m = 2e-4; p.x = -3.1; p.y = 0.7; p.z = -0.02; p.vy = -0.55; p.vz = 0.004;
    reb_simulation_add(ra, p);
    reb_simulation_add(rb, p);
    if (ra->N != 4){
        printf("  FAIL %s: N is %zu, so the add passed the mark and this case\n"
               "       tests the zeroing path instead of the re-reading one\n", label, ra->N);
        failures++;
        reb_simulation_free(ra); reb_simulation_free(rb); return;
    }
    reb_simulation_steps(ra, 200); reb_simulation_steps(rb, 200);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}


/* ------------------------------------------------------------------ */
/* The binary64 equivalence set, in the order it is run                 */
/* ------------------------------------------------------------------ */
void cases_core(void){
    case_steps("kepler, fixed step",      kepler, 2,      0.05,  0.0,   400);
    case_steps("kepler, adaptive",        kepler, 2,      0.05,  1e-9,  400);
    case_steps("kepler, tight tolerance", kepler, 2,      0.05,  1e-12, 200);
    case_steps("pythagorean, adaptive",   pythagorean, 3, 0.01,  1e-9,  400);
    case_steps("five bodies, adaptive",   five, 5,        0.5,   1e-9,  300);
    case_steps("five bodies, fixed step", five, 5,        0.5,   0.0,   200);
    /* dt = 8 on a Kepler orbit of period 2 pi is rejected repeatedly at
     * the start: the rejection branch, its er/br restore and its
     * dt/dt_last_done ratio all run here */
    case_steps("kepler, rejected steps",  kepler, 2,      8.0,   1e-9,  200);
    case_steps("pythagorean, rejected",   pythagorean, 3, 5.0,   1e-9,  200);

    case_integrate("kepler, integrate to tmax", kepler, 2, 0.05, 1e-9, 37.0);
    case_integrate("five bodies, to tmax",      five, 5,   0.5,  1e-9, 211.0);
    case_two_calls("kepler, two integrate calls", kepler, 2, 8.0, 1e-9, 20.0, 53.0);
    case_user_edit("kepler, a coordinate edited", kepler, 2, 0.05, 1e-9);
    case_add_particle("kepler, a particle added", 0.05, 1e-9);
    case_collision("a merge removes a particle", 0.01, 1e-9, 2000);
    case_remove_then_add("one removed, one added back", 0.5, 1e-9);

    printf("\n");
    /* Softening: not "is it accepted" but "is it still REBOUND's IAS15
     * with one set", which is the only question worth asking. Both
     * sides get the same value, and 0.01 is the order of the kepler
     * problem's closest approach - a softening far below that would
     * change nothing and the case would pass while exercising
     * nothing. */
    printf("\nwith r->softening = 0.01, which REBOUND adds inside the "
           "square root:\n");
    soften = 0.01;
    case_steps("kepler, softened, fixed step",   kepler, 2,      0.05, 0.0,  400);
    case_steps("kepler, softened, adaptive",     kepler, 2,      0.05, 1e-9, 400);
    case_steps("pythagorean, softened",          pythagorean, 3, 0.01, 1e-9, 400);
    soften = 0.0;

    /* min_dt. The trap here is a floor the control never reaches: the
     * case would pass and prove nothing. kepler at dt = 8 with
     * epsilon = 1e-9 is the rejected-steps configuration, where the
     * control must shrink hard, so a floor of 0.5 is well above where
     * it wants to go and is selected. The run without the floor is
     * there to prove exactly that: if the two agreed, the floor never
     * engaged. */
    printf("\nwith IAS15's min_dt, REBOUND's copysign(min_dt, dt_new) floor:\n");
    min_dt = 0.5;
    case_steps("kepler, floored at 0.5",      kepler, 2,      8.0,  1e-9, 200);
    case_steps("pythagorean, floored at 0.5", pythagorean, 3, 5.0,  1e-9, 200);
    min_dt = 0.0;

    /* AARSETH85. It shares every sum with PRS23 and differs in one
     * expression, so the risk is not that it is wrong but that it is
     * never reached - a mode the port ignored would leave these cases
     * passing as PRS23 against PRS23. The program-side check beside
     * this file measures that the two criteria choose different steps;
     * here the assertion is the usual one, that whichever REBOUND
     * picks, the port picks the same bits. */
    printf("\nwith IAS15's AARSETH85 step criterion (adaptive_mode 3):\n");
    adaptive_mode = 3;
    case_steps("kepler, A85",        kepler, 2,      0.05, 1e-9, 400);
    case_steps("pythagorean, A85",   pythagorean, 3, 0.01, 1e-9, 400);
    case_steps("five bodies, A85",   five, 5,        0.5,  1e-9, 300);
    adaptive_mode = 2;
}
