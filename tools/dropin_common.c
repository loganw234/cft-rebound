/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The drop-in gate's scaffolding: the problems every case is built
 * from, the two simulations, and the comparison.
 *
 * Split out of tools/check_dropin.c so that a topic can be a file and
 * five parcels can add cases without editing each other - see
 * tools/dropin_cases.h, which declares all of it, and PARCELS.md for
 * why. Nothing here changed in the move except `static`.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"

#include <stdint.h>


int verbose;
int failures;
int light;      /* --light: every step count scaled down */

/* A tenth of the steps, never fewer than twenty - below that a case has
 * not got past its own transient. The floor is deliberately NOT a
 * safety net for a case that must REACH something: a merge that fires
 * at step 1291 does not fire at 200, and such a case fails in light
 * mode rather than passing vacuously. Which cases those are is the
 * thing to measure, not to guess. */
long light_steps(long n){
    if (!light || n <= 20) return n;
    return n / 10 < 20 ? 20 : n / 10;
}

/* ------------------------------------------------------------------ */
/* Problems. Both simulations are built from the same doubles, so the   */
/* initial conditions cannot be the source of a difference.             */
/* ------------------------------------------------------------------ */
/* struct body, the problem arrays and the knobs are declared in
 * tools/dropin_cases.h; this file defines them. */

/* A Kepler pair at pericentre, a = 1, e = 1/2, barycentric and with
 * zero total momentum, so the pair stays where it is put and a round-off
 * difference is not buried under a linear drift. v_rel is
 * sqrt(G(M+m)(1+e)/(a(1-e))) = 1.7329166165744962, against an escape
 * velocity of 2.0009997501249219; the period is 6.2800460687587085. */
const struct body kepler[] = {
    { 1.0,  -0.00049950049950049961, 0.0, 0.0,  0.0, -0.0017311854311433533, 0.0 },
    { 1e-3,  0.49950049950049957,    0.0, 0.0,  0.0,  1.7311854311433532,    0.0 },
};

/* Burrau's Pythagorean problem: three bodies at rest, masses 3, 4, 5 */
const struct body pythagorean[] = {
    { 3.0,  1.0,  3.0, 0.0,  0.0, 0.0, 0.0 },
    { 4.0, -2.0, -1.0, 0.0,  0.0, 0.0, 0.0 },
    { 5.0,  1.0, -1.0, 0.0,  0.0, 0.0, 0.0 },
};

/* five bodies out of the plane, so that z is exercised too */
const struct body five[] = {
    { 1.0,      0.0,  0.0,   0.0,     0.0,   0.0,      0.0     },
    { 9.55e-4,  5.2,  0.0,   0.0625,  0.0,   0.438,    0.00125 },
    { 2.86e-4, -9.5,  0.0,  -0.125,   0.0,  -0.324,    0.0025  },
    { 4.37e-5,  0.0, 19.2,   0.25,   -0.228, 0.0,     -0.00375 },
    { 5.18e-5,  0.0,-30.1,  -0.5,     0.182, 0.0,      0.005   },
};

int wide_format = CFT_FP64;   /* --wide runs the cft side at binary128 */
double soften = 0.0;          /* applied by build() to BOTH sides, so
                                      * the softening case compares like with
                                      * like; set and reset around it */
double min_dt = 0.0;          /* the same, for IAS15's step floor */
int accurate = 0;             /* state->accurate: shift the wide state on a
                                      * removal instead of re-promoting it from
                                      * r->particles. Diverges from REBOUND, so it is
                                      * off for every equivalence case. */
int adaptive_mode = 2;        /* and for the step criterion: 2 PRS23,
                                      * 3 AARSETH85 */

struct reb_simulation *build(const struct body *bs, size_t n, double dt, double epsilon,
                                    int use_cft){
    struct reb_simulation *r = reb_simulation_create();
    r->G = 1.0;
    r->dt = dt;
    r->exact_finish_time = 1;
    r->softening = soften;
    if (use_cft){
        struct cft_ias15_state *s = reb_simulation_set_integrator(r, "ias15_cft");
        if (!s){ fprintf(stderr, "check_dropin: ias15_cft not registered\n"); exit(2); }
        s->epsilon = epsilon;
        s->format = wide_format;
        s->min_dt = min_dt;
        s->adaptive_mode = adaptive_mode;
        s->accurate = accurate;
    }else{
        struct reb_integrator_ias15_state *s = reb_simulation_set_integrator(r, "ias15");
        if (!s){ fprintf(stderr, "check_dropin: ias15 missing\n"); exit(2); }
        s->epsilon = epsilon;
        s->min_dt = min_dt;
        s->adaptive_mode = adaptive_mode;
    }
    for (size_t i = 0; i < n; i++){
        struct reb_particle p = {0};
        p.m = bs[i].m;
        p.x = bs[i].x; p.y = bs[i].y; p.z = bs[i].z;
        p.vx = bs[i].vx; p.vy = bs[i].vy; p.vz = bs[i].vz;
        reb_simulation_add(r, p);
    }
    return r;
}

/* ------------------------------------------------------------------ */
/* Comparison: bit patterns, never values                              */
/* ------------------------------------------------------------------ */
int bits_differ(double a, double b){ return memcmp(&a, &b, sizeof a) != 0; }

void show(const char *what, size_t i, double a, double b){
    uint64_t ua, ub;
    memcpy(&ua, &a, 8); memcpy(&ub, &b, 8);
    printf("      %-4s[%zu]  ias15 %.17g (%016llx)\n"
           "                 cft   %.17g (%016llx)\n",
           what, i, a, (unsigned long long)ua, b, (unsigned long long)ub);
}

int compare(struct reb_simulation *ra, struct reb_simulation *rb, const char *label){
    int bad = 0;
    size_t nshown = 0;
    if (ra->N != rb->N){ printf("  FAIL %s: N %zu vs %zu\n", label, ra->N, rb->N); failures++; return 1; }
    const char *names[9] = { "x", "y", "z", "vx", "vy", "vz", "ax", "ay", "az" };
    for (size_t i = 0; i < ra->N; i++){
        double A[9] = { ra->particles[i].x, ra->particles[i].y, ra->particles[i].z,
                        ra->particles[i].vx, ra->particles[i].vy, ra->particles[i].vz,
                        ra->particles[i].ax, ra->particles[i].ay, ra->particles[i].az };
        double B[9] = { rb->particles[i].x, rb->particles[i].y, rb->particles[i].z,
                        rb->particles[i].vx, rb->particles[i].vy, rb->particles[i].vz,
                        rb->particles[i].ax, rb->particles[i].ay, rb->particles[i].az };
        for (int k = 0; k < 9; k++)
            if (bits_differ(A[k], B[k])){
                bad++;
                if (nshown < 4){ printf("    %s differs for particle %zu\n", names[k], i); show(names[k], i, A[k], B[k]); nshown++; }
            }
    }
    struct { const char *n; double a, b; } clock[3] = {
        { "t",            ra->t,            rb->t            },
        { "dt",           ra->dt,           rb->dt           },
        { "dt_last_done", ra->dt_last_done, rb->dt_last_done },
    };
    for (int k = 0; k < 3; k++)
        if (bits_differ(clock[k].a, clock[k].b)){
            bad++;
            printf("    %s differs\n", clock[k].n);
            show(clock[k].n, 0, clock[k].a, clock[k].b);
        }
    if (bad){ printf("  FAIL %s: %d values differ\n", label, bad); failures++; return 1; }
    printf("  ok   %s: %zu particles x 9 values and the clock, bit for bit\n", label, ra->N);
    return 0;
}

