/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * FORCES AT EVERY SUBSTAGE: r->additional_forces, velocity-dependent
 * forces, and the two timestep-modification hooks.
 *
 * REBOUND's IAS15 does not call gravity. It calls
 * reb_simulation_update_acceleration(), which is gravity followed by
 * r->additional_forces(r), at the top of every step attempt and again
 * at each of the seven Gauss-Radau nodes, and it reads
 * particles[mk].ax straight back into at[]. So "the port supports
 * additional forces" is the claim that the port calls the SAME routine
 * the SAME number of times, with the same particles, the same
 * velocities and the same r->t, and folds the answer back the same
 * way. Every case here is that claim, asserted bit for bit against
 * REBOUND's own ias15 running the identical routine.
 *
 * Four things would each pass a weaker gate and are each tested
 * directly:
 *
 *   1. r->t at the nodes. REBOUND sets r->t = t_beginning + r->dt*h[n]
 *      before every node call. A port that left r->t alone, or that
 *      rounded its own wide clock instead of recomputing REBOUND's
 *      binary64 expression, would hand a time-dependent force a
 *      different time - and most of REBOUNDx is time-dependent.
 *      case_call_stream compares the whole sequence of (t, dt) the
 *      routine is called with, entry by entry, bit for bit.
 *   2. The node velocities. REBOUND predicts them from the b
 *      polynomial ONLY under `r->calculate_megno || (r->additional_forces
 *      && r->force_is_velocity_dependent)`; otherwise the routine sees
 *      the step-start velocities. Both halves of that are tested:
 *      case_call_stream runs once with the flag clear and once with it
 *      set, and compares the velocity it was handed each time.
 *   3. The number of calls. The stream comparison fails on a length
 *      mismatch, so a port that called the routine six times per
 *      iteration, or once per step, is caught by name rather than by a
 *      coordinate drifting.
 *   4. That the force does anything at all. case_force_changes_run
 *      requires the constant-force run to DIFFER from the no-force run;
 *      without it every case below would be a pass against itself, which
 *      is the trap this repository has walked into twice (ROADMAP.md,
 *      and PARCELS.md's rule).
 *
 * And its twin, case_nop_force_is_free: a routine that does nothing must
 * leave the run bit-identical to no routine at all. That is what says
 * the mechanism is inert when it should be, and it is the reason the
 * engine compares the returned accelerations bit for bit before
 * promoting any of them back.
 *
 * WHAT THIS CANNOT ASSERT. The user's routine is evaluated at binary64
 * even in a binary256 run: r->particles are binary64 and they are the
 * only interface REBOUND offers a callback. Gravity stays wide; a
 * component the routine writes is binary64 from that node on. There is
 * no wider REBOUND to compare against, so that ceiling is documented
 * (src/ias15_engine.h, and the integrator's own .documentation) rather
 * than tested here.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"

/* ------------------------------------------------------------------ */
/* The routines under test                                             */
/* ------------------------------------------------------------------ */
/* Index-dependent rather than uniform, deliberately: a force that is
 * the same on every body is a translation of the whole system and
 * leaves the internal dynamics - and so the step controller, and so
 * most of what is being compared - untouched. */
static void f_constant(struct reb_simulation *r){
    for (size_t i = 0; i < r->N; i++){
        const double f = 1.0e-3 * (double)(i + 1);
        r->particles[i].ax += f;
        r->particles[i].ay += -0.25 * f;
        r->particles[i].az += 0.125 * f;
    }
}

/* Reads r->t, which is the whole point. Both simulations call the same
 * sin() on the same bits in the same process, so the transcendental is
 * not a source of difference - only the argument is, which is what is
 * under test. */
static void f_time(struct reb_simulation *r){
    const double s = sin(0.7 * r->t), c = cos(0.7 * r->t);
    for (size_t i = 0; i < r->N; i++){
        const double f = 1.0e-3 * (double)(i + 1);
        r->particles[i].ax += f * s;
        r->particles[i].ay += f * c;
        r->particles[i].az += 0.5 * f * s;
    }
}

/* Stokes drag. Reads the velocities, so it is the case that fails if
 * the node velocities are not predicted from the b polynomial the way
 * integrator_ias15.c:434 predicts them. Set with
 * r->force_is_velocity_dependent = 1 on both sides. */
static void f_drag(struct reb_simulation *r){
    const double k = 1.0e-3;
    for (size_t i = 0; i < r->N; i++){
        r->particles[i].ax -= k * r->particles[i].vx;
        r->particles[i].ay -= k * r->particles[i].vy;
        r->particles[i].az -= k * r->particles[i].vz;
    }
}

static void f_nop(struct reb_simulation *r){ (void)r; }

/* ------------------------------------------------------------------ */
/* Comparing two runs that cannot be alive at the same time            */
/* ------------------------------------------------------------------ */
/* The engine is one global instance and the shim refuses a second
 * simulation while the first owns it, so a case that compares two
 * ias15_cft runs has to run them one after the other and keep the
 * first as bits. compare() in dropin_common.c takes two live
 * simulations and is no use here. */
#define SNAP_MAX 8
struct snap { size_t N; double p[SNAP_MAX][9]; double t, dt, dtl; };

static int take(struct reb_simulation *r, struct snap *s, const char *label){
    if (r->N > SNAP_MAX){
        printf("  FAIL %s: %zu particles and SNAP_MAX is %d\n", label, r->N, SNAP_MAX);
        failures++; return 1;
    }
    memset(s, 0, sizeof *s);
    s->N = r->N;
    for (size_t i = 0; i < r->N; i++){
        const struct reb_particle *q = &r->particles[i];
        const double v[9] = { q->x, q->y, q->z, q->vx, q->vy, q->vz, q->ax, q->ay, q->az };
        memcpy(s->p[i], v, sizeof v);
    }
    s->t = r->t; s->dt = r->dt; s->dtl = r->dt_last_done;
    return 0;
}

/* How many of the 9N + 3 values differ, by bit pattern. */
static int snap_diff_count(const struct snap *a, const struct snap *b){
    int bad = 0;
    if (a->N != b->N) return 1 << 20;
    for (size_t i = 0; i < a->N; i++)
        for (int k = 0; k < 9; k++)
            if (bits_differ(a->p[i][k], b->p[i][k])) bad++;
    if (bits_differ(a->t, b->t)) bad++;
    if (bits_differ(a->dt, b->dt)) bad++;
    if (bits_differ(a->dtl, b->dtl)) bad++;
    return bad;
}

/* One ias15_cft run, snapshotted and freed, so the next may have the
 * engine. force may be NULL. */
static int run_cft(struct snap *out, const char *label,
                   const struct body *bs, size_t n, double dt, double epsilon,
                   size_t steps, void (*force)(struct reb_simulation *), int veldep){
    struct reb_simulation *r = build(bs, n, dt, epsilon, 1);
    r->additional_forces = force;
    r->force_is_velocity_dependent = (unsigned)veldep;
    reb_simulation_steps(r, steps);
    int bad = take(r, out, label);
    reb_simulation_free(r);
    return bad;
}

/* ------------------------------------------------------------------ */
/* The equivalence cases                                               */
/* ------------------------------------------------------------------ */
static void case_force(const char *label, const struct body *bs, size_t n,
                       double dt, double epsilon, size_t steps,
                       void (*force)(struct reb_simulation *), int veldep){
    printf("%s (N = %zu, dt = %g, epsilon = %g, %zu steps%s)\n",
           label, n, dt, epsilon, steps, veldep ? ", velocity-dependent" : "");
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    for (int i = 0; i < 2; i++){
        struct reb_simulation *r = i ? rb : ra;
        r->additional_forces = force;
        r->force_is_velocity_dependent = (unsigned)veldep;
    }
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);
    if (verbose) printf("    t = %.17g, dt = %.17g\n", ra->t, ra->dt);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* ------------------------------------------------------------------ */
/* The call stream: every argument the routine is handed, in order     */
/* ------------------------------------------------------------------ */
/* This is the case that pins r->t at the nodes, the predicted
 * velocities and the number of calls, none of which the coordinate
 * comparison localises. The routine records what it was given and then
 * applies a real force, so the two runs are the same perturbed problem
 * rather than two unperturbed ones that agree for the wrong reason. */
#define REC_MAX 30000
struct rec_entry { double t, dt, x, v, a; };
static struct rec_entry rec[2][REC_MAX];
static size_t rec_n[2];
static int    rec_side;        /* 0 REBOUND's ias15, 1 ias15_cft */
static int    rec_overflow;

static void f_record(struct reb_simulation *r){
    if (rec_n[rec_side] < REC_MAX){
        struct rec_entry *e = &rec[rec_side][rec_n[rec_side]++];
        e->t = r->t;
        e->dt = r->dt;
        /* particle 1 stands for all of them: the position it was
         * predicted to, the velocity it was handed, and the gravity
         * computed there before this routine touched it. */
        e->x = r->particles[1].x;
        e->v = r->particles[1].vx;
        e->a = r->particles[1].ax;
    }else rec_overflow = 1;
    for (size_t i = 0; i < r->N; i++){
        const double f = 1.0e-4 * (double)(i + 1);
        r->particles[i].ax += f * sin(0.7 * r->t);
        r->particles[i].ay += f;
        r->particles[i].az -= 1.0e-3 * r->particles[i].vz;   /* and a velocity term */
    }
}

static void case_call_stream(const char *label, const struct body *bs, size_t n,
                             double dt, double epsilon, size_t steps, int veldep){
    printf("%s (every call's r->t, r->dt, x, v and gravity, in order%s)\n",
           label, veldep ? "; force_is_velocity_dependent = 1" : "");
    rec_n[0] = rec_n[1] = 0;
    rec_overflow = 0;
    for (int side = 0; side < 2; side++){
        rec_side = side;
        struct reb_simulation *r = build(bs, n, dt, epsilon, side);
        r->additional_forces = f_record;
        r->force_is_velocity_dependent = (unsigned)veldep;
        reb_simulation_steps(r, steps);
        reb_simulation_free(r);
    }
    if (rec_overflow){
        printf("  FAIL %s: more than %d calls; raise REC_MAX\n", label, REC_MAX);
        failures++; return;
    }
    if (rec_n[0] != rec_n[1]){
        printf("  FAIL %s: REBOUND called the routine %zu times and this port %zu\n",
               label, rec_n[0], rec_n[1]);
        failures++; return;
    }
    if (rec_n[0] < 8 * steps){
        printf("  FAIL %s: only %zu calls in %zu steps; REBOUND issues at least\n"
               "       eight per step (one plus the seven nodes), so this case\n"
               "       is not exercising the node loop\n", label, rec_n[0], steps);
        failures++; return;
    }
    const char *names[5] = { "t", "dt", "x", "v", "a" };
    for (size_t k = 0; k < rec_n[0]; k++){
        const double A[5] = { rec[0][k].t, rec[0][k].dt, rec[0][k].x, rec[0][k].v, rec[0][k].a };
        const double B[5] = { rec[1][k].t, rec[1][k].dt, rec[1][k].x, rec[1][k].v, rec[1][k].a };
        for (int f = 0; f < 5; f++)
            if (bits_differ(A[f], B[f])){
                printf("  FAIL %s: call %zu of %zu, the routine was handed a different %s\n",
                       label, k, rec_n[0], names[f]);
                show(names[f], k, A[f], B[f]);
                failures++;
                return;
            }
    }
    printf("  ok   %s: %zu calls, 5 values each, bit for bit\n", label, rec_n[0]);
}

/* ------------------------------------------------------------------ */
/* The two controls                                                    */
/* ------------------------------------------------------------------ */
/* A routine that does nothing must leave the run bit-identical to no
 * routine at all. Both runs are ias15_cft, so this is not equivalence
 * with REBOUND but the statement that the mechanism is inert when it
 * should be: the whole round-trip - round x, v and a out, call, read
 * back, promote - runs at every node of the first run and not at all in
 * the second, and the two must land on the same bits. */
static void case_nop_force_is_free(const char *label, const struct body *bs, size_t n,
                                   double dt, double epsilon, size_t steps){
    printf("%s (a do-nothing routine against no routine at all, both ias15_cft)\n", label);
    struct snap with, without;
    if (run_cft(&with, label, bs, n, dt, epsilon, steps, f_nop, 0)) return;
    if (run_cft(&without, label, bs, n, dt, epsilon, steps, NULL, 0)) return;
    int bad = snap_diff_count(&with, &without);
    if (bad){
        printf("  FAIL %s: %d values differ, so registering a hook that does\n"
               "       nothing is not free and the no-force path is not unchanged\n",
               label, bad);
        failures++;
        return;
    }
    printf("  ok   %s: %zu particles x 9 values and the clock, bit for bit\n", label, with.N);
}

/* And the other direction, which is the one this repository has needed
 * twice: prove the force actually reached the arithmetic. Without this
 * every case above could be comparing two runs in which the routine was
 * never called, and they would all pass. */
static void case_force_changes_run(const char *label, const struct body *bs, size_t n,
                                   double dt, double epsilon, size_t steps){
    printf("%s (the constant force against no force, both ias15_cft)\n", label);
    struct snap forced, free_run;
    if (run_cft(&forced, label, bs, n, dt, epsilon, steps, f_constant, 0)) return;
    if (run_cft(&free_run, label, bs, n, dt, epsilon, steps, NULL, 0)) return;
    int bad = snap_diff_count(&forced, &free_run);
    if (!bad){
        printf("  FAIL %s: the two runs are bit-identical, so the force never\n"
               "       reached the arithmetic and every case above is a pass\n"
               "       against itself\n", label);
        failures++;
        return;
    }
    printf("  ok   %s: %d of the %zu values differ, so the force is in the run\n",
           label, bad, 9 * forced.N + 3);
}

/* ------------------------------------------------------------------ */
/* The timestep-modification hooks                                     */
/* ------------------------------------------------------------------ */
/* These need no mechanism and never did: REBOUND's DRIVER calls them,
 * not IAS15 - reb_simulation_step() at simulation.c:523 and :564,
 * either side of the integrator callback - and the shim already
 * re-promotes any coordinate that changed under it. What was missing
 * was a case saying so, which is this. The edit is deliberately a
 * velocity and not a position, because a position that changed would
 * also be caught by the collision and user-edit cases in
 * tools/cases_core.c and this must stand on its own. */
static void pre_mod(struct reb_simulation *r){
    r->particles[1].vy += 1.0e-9;
}
static void post_mod(struct reb_simulation *r){
    r->particles[1].vz -= 2.0e-10;
    r->particles[0].vx += 5.0e-11;
}

static void case_timestep_mods(const char *label, const struct body *bs, size_t n,
                               double dt, double epsilon, size_t steps,
                               int pre, int post){
    printf("%s (%s, %zu steps)\n", label,
           pre && post ? "both hooks" : (pre ? "pre_timestep_modifications"
                                             : "post_timestep_modifications"),
           steps);
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    for (int i = 0; i < 2; i++){
        struct reb_simulation *r = i ? rb : ra;
        if (pre)  r->pre_timestep_modifications  = pre_mod;
        if (post) r->post_timestep_modifications = post_mod;
    }
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* ------------------------------------------------------------------ */
/* What must no longer be refused                                      */
/* ------------------------------------------------------------------ */
/* accepted() names the row in cft_support_rows, and the coverage check
 * in check_dropin.c walks the table afterwards - so if a later parcel
 * deletes either row outright, these fail by name and say the case is
 * orphaned. The rows still exist: the SUBPROCESS path refuses both,
 * because it hands the run to another process which cannot call a
 * function pointer in this one. */
static void p_forces(struct reb_simulation *r){ r->additional_forces = f_nop; }
static void p_veldep(struct reb_simulation *r){
    r->additional_forces = f_drag;
    r->force_is_velocity_dependent = 1;
}
static void p_both(struct reb_simulation *r){
    r->pre_timestep_modifications = pre_mod;
    r->post_timestep_modifications = post_mod;
}

/* ------------------------------------------------------------------ */
/* The topic, in the order it is run                                   */
/* ------------------------------------------------------------------ */
void cases_forces(void){
    printf("\nr->additional_forces, called where REBOUND's IAS15 calls it:\n");

    /* First the control that makes the rest mean anything. */
    case_force_changes_run("a force changes the run", kepler, 2, 0.05, 1e-9, 200);
    case_nop_force_is_free("a do-nothing force is free", kepler, 2, 0.05, 1e-9, 200);

    /* Then equivalence with REBOUND running the identical routine. */
    case_force("kepler, constant force",        kepler, 2,      0.05, 1e-9, 400, f_constant, 0);
    case_force("kepler, constant, fixed step",  kepler, 2,      0.05, 0.0,  400, f_constant, 0);
    case_force("pythagorean, constant force",   pythagorean, 3, 0.01, 1e-9, 400, f_constant, 0);
    case_force("five bodies, constant force",   five, 5,        0.5,  1e-9, 300, f_constant, 0);

    printf("\na time-dependent force, which is what reads r->t at the nodes:\n");
    case_force("kepler, time-dependent",        kepler, 2,      0.05, 1e-9, 400, f_time, 0);
    case_force("pythagorean, time-dependent",   pythagorean, 3, 0.01, 1e-9, 400, f_time, 0);
    /* dt = 8 on an orbit of period 2 pi is rejected repeatedly at the
     * start, and a rejected attempt CHANGES r->dt - so the node times
     * of the next attempt move with it. A port that handed the routine
     * the r->dt the step began with would pass everything above and
     * fail here. */
    case_force("kepler, rejected steps, time",  kepler, 2,      8.0,  1e-9, 200, f_time, 0);

    printf("\na velocity-dependent force, which is what reads the predicted v:\n");
    case_force("kepler, drag",                  kepler, 2,      0.05, 1e-9, 400, f_drag, 1);
    case_force("pythagorean, drag",             pythagorean, 3, 0.01, 1e-9, 400, f_drag, 1);
    case_force("five bodies, drag",             five, 5,        0.5,  1e-9, 300, f_drag, 1);
    /* The same drag with the flag CLEAR. REBOUND then does not predict
     * the velocities and the routine sees the step-start ones, which is
     * a different problem from the case above - and a port that always
     * predicted them would fail this one. */
    case_force("kepler, drag, flag clear",      kepler, 2,      0.05, 1e-9, 400, f_drag, 0);

    printf("\nthe call stream itself: r->t, r->dt and what the routine was handed:\n");
    case_call_stream("kepler, call stream",       kepler, 2,      0.05, 1e-9, 40, 0);
    case_call_stream("kepler, call stream, veldep", kepler, 2,    0.05, 1e-9, 40, 1);
    case_call_stream("pythagorean, call stream",  pythagorean, 3, 0.01, 1e-9, 40, 0);
    case_call_stream("kepler, rejected, stream",  kepler, 2,      8.0,  1e-9, 40, 0);

    printf("\npre_ and post_timestep_modifications, which REBOUND's driver calls:\n");
    case_timestep_mods("kepler, post only",  kepler, 2,      0.05, 1e-9, 300, 0, 1);
    case_timestep_mods("kepler, pre only",   kepler, 2,      0.05, 1e-9, 300, 1, 0);
    case_timestep_mods("pythagorean, both",  pythagorean, 3, 0.01, 1e-9, 300, 1, 1);

    printf("\nand what must no longer be refused:\n");
    accepted("additional_forces", "r->additional_forces",          p_forces);
    accepted("veldep_forces",     "a velocity-dependent force",    p_veldep);
    accepted(NULL,                "the timestep-modification hooks", p_both);
    printf("\n");
}
