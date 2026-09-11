/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * THE GATE FOR THE DROP-IN, seen from REBOUND's side.
 *
 * The same program, twice: once with REBOUND's own "ias15" and once with
 * the registered "ias15_cft" at binary64. Every particle value must
 * agree BIT FOR BIT - not to a tolerance, not to a relative error, the
 * same 64 bits - and so must r->t, r->dt and r->dt_last_done.
 *
 * That is the existing equivalence gate (tools/check_equivalence.py)
 * seen from the other end. check_equivalence.py compares two record
 * files written by two standalone programs; this compares two live
 * struct reb_simulation, which is the only thing that can prove the shim
 * itself - the promotion, the step boundary, the clock, the view written
 * back - did not disturb the arithmetic.
 *
 * The last case is not an equivalence case and says so: it checks that
 * every unsupported feature is refused by name rather than computed.
 * The one before it, a particle added mid-run, IS an equivalence case -
 * it compares against REBOUND on the grown simulation - and is listed
 * separately only because what it exercises is the add hook.
 *
 *   build/check_dropin            all cases, at binary64
 *   build/check_dropin --wide     the same shim at binary128, which is
 *                                 the second half of the gate and is
 *                                 what `make check` runs after the first
 *   build/check_dropin -v         and print the first differing value
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

#include "rebound.h"
#include "cft_ias15.h"

static int verbose;
static int failures;

/* ------------------------------------------------------------------ */
/* Problems. Both simulations are built from the same doubles, so the   */
/* initial conditions cannot be the source of a difference.             */
/* ------------------------------------------------------------------ */
struct body { double m, x, y, z, vx, vy, vz; };

/* A Kepler pair at pericentre, a = 1, e = 1/2, barycentric and with
 * zero total momentum, so the pair stays where it is put and a round-off
 * difference is not buried under a linear drift. v_rel is
 * sqrt(G(M+m)(1+e)/(a(1-e))) = 1.7329166165744962, against an escape
 * velocity of 2.0009997501249219; the period is 6.2800460687587085. */
static const struct body kepler[] = {
    { 1.0,  -0.00049950049950049961, 0.0, 0.0,  0.0, -0.0017311854311433533, 0.0 },
    { 1e-3,  0.49950049950049957,    0.0, 0.0,  0.0,  1.7311854311433532,    0.0 },
};

/* Burrau's Pythagorean problem: three bodies at rest, masses 3, 4, 5 */
static const struct body pythagorean[] = {
    { 3.0,  1.0,  3.0, 0.0,  0.0, 0.0, 0.0 },
    { 4.0, -2.0, -1.0, 0.0,  0.0, 0.0, 0.0 },
    { 5.0,  1.0, -1.0, 0.0,  0.0, 0.0, 0.0 },
};

/* five bodies out of the plane, so that z is exercised too */
static const struct body five[] = {
    { 1.0,      0.0,  0.0,   0.0,     0.0,   0.0,      0.0     },
    { 9.55e-4,  5.2,  0.0,   0.0625,  0.0,   0.438,    0.00125 },
    { 2.86e-4, -9.5,  0.0,  -0.125,   0.0,  -0.324,    0.0025  },
    { 4.37e-5,  0.0, 19.2,   0.25,   -0.228, 0.0,     -0.00375 },
    { 5.18e-5,  0.0,-30.1,  -0.5,     0.182, 0.0,      0.005   },
};

static int wide_format = CFT_FP64;   /* --wide runs the cft side at binary128 */

static struct reb_simulation *build(const struct body *bs, size_t n, double dt, double epsilon,
                                    int use_cft){
    struct reb_simulation *r = reb_simulation_create();
    r->G = 1.0;
    r->dt = dt;
    r->exact_finish_time = 1;
    if (use_cft){
        struct cft_ias15_state *s = reb_simulation_set_integrator(r, "ias15_cft");
        if (!s){ fprintf(stderr, "check_dropin: ias15_cft not registered\n"); exit(2); }
        s->epsilon = epsilon;
        s->format = wide_format;
    }else{
        struct reb_integrator_ias15_state *s = reb_simulation_set_integrator(r, "ias15");
        if (!s){ fprintf(stderr, "check_dropin: ias15 missing\n"); exit(2); }
        s->epsilon = epsilon;
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
static int bits_differ(double a, double b){ return memcmp(&a, &b, sizeof a) != 0; }

static void show(const char *what, size_t i, double a, double b){
    uint64_t ua, ub;
    memcpy(&ua, &a, 8); memcpy(&ub, &b, 8);
    printf("      %-4s[%zu]  ias15 %.17g (%016llx)\n"
           "                 cft   %.17g (%016llx)\n",
           what, i, a, (unsigned long long)ua, b, (unsigned long long)ub);
}

static int compare(struct reb_simulation *ra, struct reb_simulation *rb, const char *label){
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

/* A particle added mid-run. REBOUND's ias15_alloc zeroes g, e, b, csb,
 * er, br, csx and csv whenever the array grows past its high-water mark
 * (realloc_dp7 zeroes the whole array); the shim's did_add_particle
 * invalidates and the next step does the same. So this IS an
 * equivalence case - on a grow. Removal is not, and is not claimed:
 * REBOUND keeps a polynomial that no longer describes the particle set
 * and the shim deliberately does not. */
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

/* ------------------------------------------------------------------ */
/* The refusals: named, not computed                                    */
/* ------------------------------------------------------------------ */
static void nop_forces(struct reb_simulation *r){ (void)r; }

/* The message itself goes to stderr, which is where a user would meet it;
 * what is asserted here is that the integration STOPPED and that the
 * clock did not move - nothing was computed. Positions are deliberately
 * not asserted: REBOUND's own boundary check and collision search run
 * after the integrator callback returns and may touch them, and that is
 * REBOUND's doing, not this integrator's. */
static int refused(const char *what, void (*poison)(struct reb_simulation *)){
    struct reb_simulation *r = build(kepler, 2, 0.05, 1e-9, 1);
    poison(r);
    double t_before = r->t, dtl_before = r->dt_last_done;
    /* the integrator's own callback, not reb_simulation_steps: what is
     * under test is this integrator's refusal, and REBOUND's post-step
     * boundary check and collision search have no business running with
     * a box that was never configured. */
    r->integrator.callbacks.step(r, r->integrator.state);
    int stopped = (r->status == REB_STATUS_GENERIC_ERROR);
    int still = (memcmp(&t_before, &r->t, sizeof(double)) == 0 &&
                 memcmp(&dtl_before, &r->dt_last_done, sizeof(double)) == 0);
    int status = (int)r->status;
    reb_simulation_free(r);
    if (stopped && still){ printf("  ok   %s refused, clock did not move\n", what); return 0; }
    printf("  FAIL %s: status %d (wanted %d), clock %s\n",
           what, status, (int)REB_STATUS_GENERIC_ERROR, still ? "still" : "MOVED");
    failures++;
    return 1;
}

static void p_softening(struct reb_simulation *r){ r->softening = 1e-6; }
static void p_forces(struct reb_simulation *r){ r->additional_forces = nop_forces; }
static void p_veldep(struct reb_simulation *r){ r->force_is_velocity_dependent = 1; }
static void p_ghost(struct reb_simulation *r){ r->N_ghost_x = 1; }
static void p_gravity(struct reb_simulation *r){ r->gravity = REB_GRAVITY_COMPENSATED; }
static void p_tree(struct reb_simulation *r){ r->gravity = REB_GRAVITY_TREE; }
static void p_collision(struct reb_simulation *r){ r->collision = REB_COLLISION_DIRECT; }
static void p_boundary(struct reb_simulation *r){ r->boundary = REB_BOUNDARY_PERIODIC; }
static void p_testp(struct reb_simulation *r){ r->N_active = 1; }
static void p_var(struct reb_simulation *r){ reb_simulation_add_variation_1st_order(r, -1); }
static void p_megno(struct reb_simulation *r){ r->calculate_megno = 1; }
/* The inverse of refused(): a feature this integrator must ACCEPT.
 * Worth its own helper because "the clock moved and nothing errored"
 * is the whole assertion, and writing it inline three times would
 * invite one of them to be written differently. */
static int accepted(const char *what, void (*prepare)(struct reb_simulation *)){
    struct reb_simulation *r = build(kepler, 2, 0.05, 1e-9, 1);
    prepare(r);
    double t_before = r->t;
    r->integrator.callbacks.step(r, r->integrator.state);
    int ok = (r->status != REB_STATUS_GENERIC_ERROR) &&
             (memcmp(&t_before, &r->t, sizeof(double)) != 0);
    int status = (int)r->status;
    reb_simulation_free(r);
    if (ok){ printf("  ok   %s accepted, the clock moved\n", what); return 0; }
    printf("  FAIL %s: status %d, and the clock did not move\n", what, status);
    failures++;
    return 1;
}

static void nop_gravity(struct reb_simulation *r){ (void)r; }
static void p_gravcustom(struct reb_simulation *r){ r->gravity_custom = nop_gravity; }
/* Through REBOUND's own constructor, so N_odes is set the way a user
 * would set it. reb_simulation_step() integrates it with its own
 * Bulirsch-Stoer state because this integrator is not named "bs",
 * exactly as it does under REBOUND's ias15 - so this must be
 * accepted, not refused. */
static void p_odes(struct reb_simulation *r){ reb_ode_create(r, 1); }
static void p_mindt(struct reb_simulation *r){ cft_ias15_get_state(r)->min_dt = 1e-8; }
static void p_mode(struct reb_simulation *r){ cft_ias15_get_state(r)->adaptive_mode = 0; }

static void case_refusals(void){
    printf("refusals: each unsupported feature named and the integration stopped\n");
    refused("non-zero softening",        p_softening);
    refused("additional_forces",         p_forces);
    refused("velocity-dependent forces", p_veldep);
    refused("ghost boxes",               p_ghost);
    refused("REB_GRAVITY_COMPENSATED",   p_gravity);
    refused("the tree code",             p_tree);
    refused("a custom gravity routine",  p_gravcustom);
    refused("collision detection",       p_collision);
    refused("periodic boundaries",       p_boundary);
    refused("test particles (N_active)", p_testp);
    refused("variational particles",     p_var);
    refused("MEGNO",                     p_megno);
    refused("min_dt",                    p_mindt);
    refused("adaptive_mode != PRS23",    p_mode);

    /* And one that must NOT be refused. It was, for one day, on a
     * premise that was not true. */
    accepted("an attached ODE set",      p_odes);
}

/* ------------------------------------------------------------------ */
int main(int argc, char **argv){
    int wide = 0;
    setvbuf(stdout, NULL, _IOLBF, 0);   /* so a pipe sees every case as it happens */
    for (int i = 1; i < argc; i++)
        if (!strcmp(argv[i], "-v")) verbose = 1;
        else if (!strcmp(argv[i], "--wide")){ wide = 1; wide_format = CFT_FP128; }
        else { fprintf(stderr, "usage: check_dropin [-v] [--wide]\n"); return 2; }

    cft_ias15_register("ias15_cft");
    /* Several of the step's scratch vectors are sized at their first use,
     * so the largest body count of the whole process is declared here. */
    cft_ias15_reserve(8);

    if (wide){
        printf("check_dropin --wide: the same shim at binary128, a smoke test\n\n");
        case_wide();
        printf("\n");
        if (failures){ printf("check_dropin --wide: %d FAILURES\n", failures); return 1; }
        printf("check_dropin --wide: passed\n");
        return 0;
    }

    printf("check_dropin: REBOUND's own ias15 against the registered "
           "ias15_cft at binary64, bit for bit\n\n");

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

    printf("\n");
    case_refusals();

    printf("\n");
    if (failures){ printf("check_dropin: %d FAILURES\n", failures); return 1; }
    printf("check_dropin: every case passed\n");
    return 0;
}
