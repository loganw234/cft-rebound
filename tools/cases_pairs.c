/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * WHICH PAIRS ARE COMPUTED: r->N_active, r->testparticle_type and
 * r->gravity_ignore_terms.
 *
 * reb_gravity_basic_calculate_acceleration() (gravity.c:167) is two
 * loops and three settings, not one loop over the triangle. All three
 * change the SET of addends a particle receives and the ORDER it
 * receives them in, so bit-identity is a claim about the sequence and
 * not about the sum - and the port reproduces the loops rather than the
 * total. The cases below are that claim.
 *
 * FOUR THINGS HERE WOULD EACH PASS A WEAKER GATE.
 *
 * 1. MASSLESS TEST PARTICLES MAKE N_active BIT-INERT, which is the trap
 *    this whole topic sets. Work it through: with m = 0 above N_active,
 *    every pair the skip removes contributes prefact*0 = +-0 to
 *    whichever particle would have received it, and the running sum
 *    starts at +0 and can never be -0 (a cancellation rounds to +0, and
 *    +0 + -0 is +0), so x + (+-0) returns x's bits. The real
 *    contributions also stay in the same relative order, because the
 *    massive particles are the low indices and every list here is
 *    ascending. So the run with N_active set is bit-identical to the run
 *    without it - and a gate that tested test particles the way most
 *    scripts write them would pass whether or not any of this was
 *    implemented. case_massless_is_inert asserts exactly that, so the
 *    fact is recorded rather than relied on, and every DIVERGENCE
 *    control below uses inactive particles that carry MASS.
 * 2. r->gravity_ignore_terms cannot be reached through IAS15 at all.
 *    reb_integrator_ias15_step() writes NONE over it at the top of every
 *    step (integrator_ias15.c:875), so no end-to-end case can tell the
 *    three values apart - REBOUND's own answer does not change either.
 *    case_clobbers_ignore_terms measures that on REBOUND rather than
 *    asserting it, and the pair selection itself is then tested one
 *    configuration at a time against REBOUND's gravity directly.
 * 3. Comparing integrations only would never see gravity alone, and
 *    would only ever see the handful of configurations someone thought
 *    to write down. The oracle cases put this port's gravity beside
 *    reb_simulation_update_acceleration()'s for one configuration and
 *    compare all 3N accelerations bit for bit, which localises a wrong
 *    pair set to the pair set; case_gravity_sweep then does that for
 *    EVERY configuration the engine will accept - N from 1, N_active
 *    from 0, both testparticle_type values, all three ignore_terms -
 *    and reports how many of them actually differ from the
 *    unrestricted answer, so a sweep that agreed by never restricting
 *    anything would say so.
 * 4. That any of it changed the answer. case_active_changes_run and
 *    case_tptype_changes_run require the restricted run to DIFFER from
 *    the unrestricted one and print by how much, because two runs that
 *    were never restricted agree for the wrong reason.
 *
 * WHAT IS NOT HERE. r->map is still refused, and tools/check_dropin.c
 * still poisons it: under IAS15 the map selects which particles are
 * INTEGRATED (integrator_ias15.c:287-289 and :314, the state is 3*N_map
 * long and indexed through map[]) while gravity still runs over all
 * r->N - so the engine would need its integrated set and its gravity
 * set to be different lengths, which every vector in step_attempt()
 * assumes they are not. The two in-tree users, MERCURIUS and TRACE, set
 * r->gravity = REB_GRAVITY_CUSTOM in the same breath and are refused by
 * that row anyway.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"
#include "ias15_engine.h"

/* REBOUND's enum by value; rebound.h:282-284. Spelled out so the tables
 * below can be indexed by it. */
#define IGN_N 3
static const char *ign_name[IGN_N] = { "NONE", "BETWEEN_0_AND_1", "INVOLVING_0" };

#define PAIRS_N 6     /* every problem here; check_dropin reserves 8 */

/* ------------------------------------------------------------------ */
/* Two six-body problems that differ only in the last three masses     */
/* ------------------------------------------------------------------ */
/* A star, two planets, and three test particles on bound near-circular
 * orbits. This is the configuration a planetesimal or debris script
 * actually has, and by the argument at the top of this file it is the
 * one in which N_active changes nothing - which is why it is here as a
 * CONTROL and not as the evidence. */
static const struct body dust[PAIRS_N] = {
    { 1.0,       0.0,   0.0,  0.0,     0.0,    0.0,    0.0     },
    { 9.55e-4,   5.2,   0.0,  0.0625,  0.0,    0.438,  0.00125 },
    { 2.86e-4,  -9.5,   0.0, -0.125,   0.0,   -0.324,  0.0025  },
    { 0.0,       2.0,   0.0,  0.0,     0.0,    0.70,   0.0     },
    { 0.0,       0.0,   3.1,  0.05,   -0.56,   0.0,    0.002   },
    { 0.0,      -7.3,   0.0, -0.2,     0.0,   -0.37,  -0.001   },
};

/* The same six bodies with mass on the last three: "semi-active"
 * particles, which feel the active ones and, at testparticle_type = 0,
 * do not pull on them or on each other. Both skips are live here - the
 * dropped test-test pairs and the dropped back-reaction - so this is the
 * problem every divergence control below runs. */
static const struct body semi[PAIRS_N] = {
    { 1.0,       0.0,   0.0,  0.0,     0.0,    0.0,    0.0     },
    { 9.55e-4,   5.2,   0.0,  0.0625,  0.0,    0.438,  0.00125 },
    { 2.86e-4,  -9.5,   0.0, -0.125,   0.0,   -0.324,  0.0025  },
    { 1.0e-5,    2.0,   0.0,  0.0,     0.0,    0.70,   0.0     },
    { 2.0e-5,    0.0,   3.1,  0.05,   -0.56,   0.0,    0.002   },
    { 3.0e-5,   -7.3,   0.0, -0.2,     0.0,   -0.37,  -0.001   },
};

/* ------------------------------------------------------------------ */
/* Snapshots: two ias15_cft runs cannot be alive at the same time      */
/* ------------------------------------------------------------------ */
struct psnap { size_t N; double p[PAIRS_N][9]; double t, dt, dtl; };

static int ptake(struct reb_simulation *r, struct psnap *s, const char *label){
    if (r->N > PAIRS_N){
        printf("  FAIL %s: %zu particles and PAIRS_N is %d\n", label, r->N, PAIRS_N);
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

/* How many of the 9N + 3 values differ by bit pattern, and - so that a
 * control can report a NUMBER rather than the word "differs" - the
 * largest absolute and relative gap and where it was.
 *
 * `body` stays -1 when nothing could be measured: only the clock
 * differed, or every gap was a NaN and so compared false against every
 * threshold. Both are worth saying out loud rather than indexing with,
 * which is what say_worst() is for. */
struct pdiff { int n; double worst_abs, worst_rel; int body, field; };

static const char *field_name[9] = { "x", "y", "z", "vx", "vy", "vz", "ax", "ay", "az" };

static void say_worst(const struct pdiff *d, double a, double b){
    if (d->body < 0){
        printf("       the differing values are the clock, or the gaps are not "
               "ordered (a NaN)\n");
        return;
    }
    printf("       worst %s[%d]: %.17g against %.17g\n",
           field_name[d->field], d->body, a, b);
    printf("       max |delta| %.6e, relative %.6e\n", d->worst_abs, d->worst_rel);
}

/* The value say_worst() names, or 0 when there is nothing to name. */
static double at_worst(const struct psnap *s, const struct pdiff *d){
    return d->body < 0 ? 0.0 : s->p[d->body][d->field];
}

static struct pdiff psnap_diff(const struct psnap *a, const struct psnap *b){
    struct pdiff d;
    memset(&d, 0, sizeof d);
    d.body = -1;
    if (a->N != b->N){ d.n = 1 << 20; return d; }
    for (size_t i = 0; i < a->N; i++)
        for (int k = 0; k < 9; k++){
            if (!bits_differ(a->p[i][k], b->p[i][k])) continue;
            d.n++;
            const double gap = fabs(a->p[i][k] - b->p[i][k]);
            const double scale = fabs(a->p[i][k]) > fabs(b->p[i][k])
                               ? fabs(a->p[i][k]) : fabs(b->p[i][k]);
            if (gap > d.worst_abs){
                d.worst_abs = gap;
                d.worst_rel = scale > 0 ? gap / scale : 0.0;
                d.body = (int)i; d.field = k;
            }
        }
    if (bits_differ(a->t, b->t)) d.n++;
    if (bits_differ(a->dt, b->dt)) d.n++;
    if (bits_differ(a->dtl, b->dtl)) d.n++;
    return d;
}

/* One ias15_cft run, snapshotted and freed so the next may have the
 * engine. */
static int run_pairs(struct psnap *out, const char *label,
                     const struct body *bs, size_t n, double dt, double epsilon,
                     size_t steps, size_t n_active, int tp_type){
    struct reb_simulation *r = build(bs, n, dt, epsilon, 1);
    r->N_active = n_active;
    r->testparticle_type = tp_type;
    r->testparticle_hidewarnings = 1;
    reb_simulation_steps(r, steps);
    int bad = ptake(r, out, label);
    reb_simulation_free(r);
    return bad;
}

/* ------------------------------------------------------------------ */
/* Gravity alone, one configuration at a time                          */
/* ------------------------------------------------------------------ */
/* The only case shape that can reach r->gravity_ignore_terms, because
 * IAS15 overwrites the field before an integration ever gets there. It
 * is also the sharpest form of the rest: a wrong pair set shows up as a
 * wrong acceleration on a named particle instead of as a trajectory
 * that drifted.
 *
 * The engine is one global instance opened by the shim on its first
 * step, so a run happens first and is thrown away; after that the
 * buffers are live and can be driven directly. Everything gravity reads
 * is written here - masses, G, softening, positions, and the three
 * settings - so nothing carries over from whatever ran before. */
static int engine_gravity(const struct body *bs, size_t n, double G,
                          size_t n_active, int tp_type, int ign, double *a_out){
    double m[PAIRS_N], x[3 * PAIRS_N], v[3 * PAIRS_N];
    for (size_t i = 0; i < n; i++){
        m[i] = bs[i].m;
        x[3*i] = bs[i].x; x[3*i+1] = bs[i].y; x[3*i+2] = bs[i].z;
        v[3*i] = bs[i].vx; v[3*i+1] = bs[i].vy; v[3*i+2] = bs[i].vz;
    }
    if (ias15_engine_set_bodies(n) != 0) return 1;
    ias15_engine_set_G_f64(G);
    ias15_engine_set_softening_f64(0.0);
    ias15_engine_set_masses_f64(m);
    ias15_engine_put_xv_f64(x, v);
    if (ias15_engine_set_active(n_active, tp_type, ign) != 0) return 1;
    ias15_engine_gravity();
    ias15_engine_get_a_f64(a_out);
    /* Leave the engine as the next case expects to find it. The shim
     * sets all three on every step anyway; this is so that a failure
     * here cannot become a failure somewhere else. */
    ias15_engine_set_active((size_t)-1, 0, 0);
    return 0;
}

/* REBOUND's own answer for the same configuration.
 * reb_simulation_update_acceleration() is what IAS15 calls, and with no
 * additional_forces attached it is reb_gravity_basic_calculate_
 * acceleration() and nothing else. */
static void rebound_gravity(const struct body *bs, size_t n, double G,
                            size_t n_active, int tp_type, int ign, double *a_out){
    struct reb_simulation *r = build(bs, n, 0.05, 1e-9, 0);
    r->G = G;
    r->N_active = n_active;
    r->testparticle_type = tp_type;
    r->gravity_ignore_terms = ign;
    r->testparticle_hidewarnings = 1;
    reb_simulation_update_acceleration(r);
    for (size_t i = 0; i < n; i++){
        a_out[3*i] = r->particles[i].ax;
        a_out[3*i+1] = r->particles[i].ay;
        a_out[3*i+2] = r->particles[i].az;
    }
    reb_simulation_free(r);
}

/* How many of the 3n accelerations differ, and the largest gap. */
static struct pdiff acc_diff(const double *a, const double *b, size_t n){
    struct pdiff d;
    memset(&d, 0, sizeof d);
    d.body = -1;
    for (size_t k = 0; k < 3 * n; k++){
        if (!bits_differ(a[k], b[k])) continue;
        d.n++;
        const double gap = fabs(a[k] - b[k]);
        const double scale = fabs(a[k]) > fabs(b[k]) ? fabs(a[k]) : fabs(b[k]);
        if (gap > d.worst_abs){
            d.worst_abs = gap;
            d.worst_rel = scale > 0 ? gap / scale : 0.0;
            d.body = (int)(k / 3); d.field = 6 + (int)(k % 3);
        }
    }
    return d;
}

/* One configuration: this port's gravity against REBOUND's, bit for
 * bit, AND the measured distance from the unrestricted answer, so that
 * a configuration whose restriction did nothing is reported as such
 * instead of passing quietly. */
static void case_gravity(const char *label, const struct body *bs, size_t n,
                         size_t n_active, int tp_type, int ign){
    double mine[3 * PAIRS_N], theirs[3 * PAIRS_N], plain[3 * PAIRS_N];
    printf("  %-46s", label);
    if (engine_gravity(bs, n, 1.0, n_active, tp_type, ign, mine)){
        printf("FAIL: the engine refused the configuration\n"); failures++; return;
    }
    rebound_gravity(bs, n, 1.0, n_active, tp_type, ign, theirs);
    struct pdiff d = acc_diff(mine, theirs, n);
    if (d.n){
        printf("FAIL\n");
        printf("    %d of the %zu accelerations differ from REBOUND's\n", d.n, 3 * n);
        if (d.body >= 0)
            show(field_name[d.field], (size_t)d.body,
                 theirs[3 * d.body + (d.field - 6)], mine[3 * d.body + (d.field - 6)]);
        else
            printf("    and no gap is ordered, so at least one side is a NaN\n");
        failures++;
        return;
    }
    /* and the distance from "no restriction at all", which is what says
     * this configuration is a different problem */
    rebound_gravity(bs, n, 1.0, (size_t)-1, 0, 0, plain);
    struct pdiff s = acc_diff(mine, plain, n);
    if (s.n)
        printf("ok   3N bit for bit; %d differ from unrestricted, max |d| %.3e (rel %.2e)\n",
               s.n, s.worst_abs, s.worst_rel);
    else
        printf("ok   3N bit for bit; identical to unrestricted (no term was dropped)\n");
}

/* Every configuration the engine will accept, in one line of output.
 *
 * The named cases above are readable evidence for the settings a user
 * would actually write; this is the coverage. N runs from 1, where
 * there are no pairs at all, and N_active from 0, where nothing is a
 * source - the two ends the loop bounds are most likely to get wrong,
 * and the place MAX(N_active, starti) stops being N_active. Every
 * combination is compared against REBOUND's own gravity bit for bit,
 * and the count of configurations that actually differ from the
 * unrestricted answer is reported beside it, so a sweep that agreed
 * with REBOUND by never restricting anything would be visible. */
static void case_gravity_sweep(const char *label, const struct body *bs, size_t nmax){
    int checked = 0, diverged = 0, bad = 0, shown = 0;
    for (size_t n = 1; n <= nmax; n++)
        for (int ign = 0; ign < IGN_N; ign++)
            for (int tp = 0; tp < 2; tp++)
                for (size_t na = 0; na <= n + 1; na++){
                    /* n + 1 stands for REBOUND's SIZE_MAX default; every
                     * value at or below n is a real N_active, and above
                     * n is the heap overread the support table refuses */
                    const size_t n_active = (na == n + 1) ? (size_t)-1 : na;
                    double mine[3 * PAIRS_N], theirs[3 * PAIRS_N], plain[3 * PAIRS_N];
                    if (engine_gravity(bs, n, 1.0, n_active, tp, ign, mine)){
                        printf("  FAIL %s: the engine refused N %zu, N_active %zu, "
                               "tp %d, ignore %s\n", label, n, na, tp, ign_name[ign]);
                        failures++; return;
                    }
                    rebound_gravity(bs, n, 1.0, n_active, tp, ign, theirs);
                    checked++;
                    struct pdiff d = acc_diff(mine, theirs, n);
                    if (d.n){
                        bad++;
                        if (shown++ < 3){
                            printf("  FAIL %s: N %zu, N_active %zu, tp %d, ignore %s - "
                                   "%d of %zu accelerations differ\n",
                                   label, n, na, tp, ign_name[ign], d.n, 3 * n);
                            if (d.body >= 0)
                                show(field_name[d.field], (size_t)d.body,
                                     theirs[3 * d.body + (d.field - 6)],
                                     mine[3 * d.body + (d.field - 6)]);
                        }
                        continue;
                    }
                    rebound_gravity(bs, n, 1.0, (size_t)-1, 0, 0, plain);
                    if (acc_diff(mine, plain, n).n) diverged++;
                }
    if (bad){
        printf("  FAIL %s: %d of %d configurations disagree with REBOUND\n",
               label, bad, checked);
        failures++;
        return;
    }
    printf("  ok   %s: %d configurations (N 1..%zu x N_active 0..N and unset x "
           "testparticle_type x\n       all three ignore_terms), every one bit for "
           "bit; %d of them differ from\n       the unrestricted answer\n",
           label, checked, nmax, diverged);
}

/* ------------------------------------------------------------------ */
/* Equivalence: whole integrations against REBOUND's own ias15         */
/* ------------------------------------------------------------------ */
static void case_pairs(const char *label, const struct body *bs, size_t n,
                       double dt, double epsilon, size_t steps,
                       size_t n_active, int tp_type){
    printf("%s (N = %zu, N_active = %zu, testparticle_type = %d, %zu steps)\n",
           label, n, n_active, tp_type, light_ran(steps));
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    for (int i = 0; i < 2; i++){
        struct reb_simulation *r = i ? rb : ra;
        r->N_active = n_active;
        r->testparticle_type = tp_type;
        r->testparticle_hidewarnings = 1;
    }
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);
    if (verbose) printf("    t = %.17g, dt = %.17g\n", ra->t, ra->dt);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* ------------------------------------------------------------------ */
/* The controls                                                        */
/* ------------------------------------------------------------------ */
/* Prove the skip reached the arithmetic: the same problem with and
 * without N_active must land on different bits. Run on `semi`, whose
 * inactive particles carry mass - on `dust` this cannot fail and cannot
 * pass, which is the next case. */
static void case_active_changes_run(const char *label, const struct body *bs, size_t n,
                                    double dt, double epsilon, size_t steps,
                                    size_t n_active){
    printf("%s (N_active = %zu against every particle active, both ias15_cft)\n",
           label, n_active);
    struct psnap restricted, all;
    if (run_pairs(&restricted, label, bs, n, dt, epsilon, steps, n_active, 0)) return;
    if (run_pairs(&all, label, bs, n, dt, epsilon, steps, (size_t)-1, 0)) return;
    struct pdiff d = psnap_diff(&restricted, &all);
    if (!d.n){
        printf("  FAIL %s: the two runs are bit-identical, so the skip never reached\n"
               "       the arithmetic and every case above is a pass against itself\n",
               label);
        failures++;
        return;
    }
    printf("  ok   %s: %d of the %zu values differ\n",
           label, d.n, 9 * restricted.N + 3);
    say_worst(&d, at_worst(&restricted, &d), at_worst(&all, &d));
}

/* The same for testparticle_type, which is the knob that was not in the
 * refusal table at all while N_active was refused, and which becomes
 * live the moment N_active is not. */
static void case_tptype_changes_run(const char *label, const struct body *bs, size_t n,
                                    double dt, double epsilon, size_t steps,
                                    size_t n_active){
    printf("%s (testparticle_type 1 against 0, N_active = %zu, both ias15_cft)\n",
           label, n_active);
    struct psnap felt, unfelt;
    if (run_pairs(&felt, label, bs, n, dt, epsilon, steps, n_active, 1)) return;
    if (run_pairs(&unfelt, label, bs, n, dt, epsilon, steps, n_active, 0)) return;
    struct pdiff d = psnap_diff(&felt, &unfelt);
    if (!d.n){
        printf("  FAIL %s: the two runs are bit-identical, so testparticle_type is\n"
               "       being ignored - which is the silent wrong answer this case\n"
               "       exists for\n", label);
        failures++;
        return;
    }
    printf("  ok   %s: %d of the %zu values differ\n",
           label, d.n, 9 * felt.N + 3);
    say_worst(&d, at_worst(&felt, &d), at_worst(&unfelt, &d));
}

/* And the finding the two controls above are shaped around: with
 * MASSLESS test particles the skip is exactly bit-inert, so a gate
 * written on the configuration most scripts actually use would pass
 * whether or not any of this existed. Asserted rather than assumed,
 * because if it ever stopped being true - a NaN leaking out of a pair
 * the list no longer names, say - this is where it would show. */
static void case_massless_is_inert(const char *label, const struct body *bs, size_t n,
                                   double dt, double epsilon, size_t steps,
                                   size_t n_active){
    printf("%s (massless test particles: N_active must change NOTHING)\n", label);
    struct psnap restricted, all;
    if (run_pairs(&restricted, label, bs, n, dt, epsilon, steps, n_active, 0)) return;
    if (run_pairs(&all, label, bs, n, dt, epsilon, steps, (size_t)-1, 0)) return;
    struct pdiff d = psnap_diff(&restricted, &all);
    if (d.n){
        printf("  FAIL %s: %d values differ. Every pair the skip removes contributes\n"
               "       prefact*0 = +-0, and adding +-0 to a sum that started at +0\n"
               "       cannot change its bits - so a difference here means a dropped\n"
               "       pair produced something that is not a zero\n", label, d.n);
        failures++;
        return;
    }
    printf("  ok   %s: bit for bit, as the arithmetic requires - which is why the\n"
           "       divergence controls run on massive inactive particles instead\n", label);
}

/* REBOUND's IAS15 writes NONE over r->gravity_ignore_terms at the top of
 * every step (integrator_ias15.c:875), so no value a user sets reaches
 * gravity under it. Measured on REBOUND ITSELF rather than asserted:
 * this case would fail if IAS15 ever honoured the field, and that is the
 * day the port would have to as well. */
static void case_clobbers_ignore_terms(const char *label, const struct body *bs, size_t n,
                                       double dt, double epsilon, size_t steps){
    printf("%s (REBOUND's own ias15, each value against NONE)\n", label);
    struct psnap base;
    {
        struct reb_simulation *r = build(bs, n, dt, epsilon, 0);
        reb_simulation_steps(r, steps);
        if (ptake(r, &base, label)){ reb_simulation_free(r); return; }
        reb_simulation_free(r);
    }
    for (int ign = 1; ign < IGN_N; ign++){
        struct psnap s;
        struct reb_simulation *r = build(bs, n, dt, epsilon, 0);
        r->gravity_ignore_terms = ign;
        reb_simulation_steps(r, steps);
        int bad = ptake(r, &s, label);
        const int left = (int)r->gravity_ignore_terms;
        reb_simulation_free(r);
        if (bad) return;
        struct pdiff d = psnap_diff(&s, &base);
        if (d.n){
            printf("  FAIL %s: %s changed REBOUND's own answer in %d values, so IAS15\n"
                   "       does honour the field and this port must too\n",
                   label, ign_name[ign], d.n);
            failures++;
            return;
        }
        printf("  ok   %s = %-15s no change, and the field reads %s after the step\n",
               label, ign_name[ign], ign_name[left]);
    }
}

/* ------------------------------------------------------------------ */
/* What must no longer be refused                                      */
/* ------------------------------------------------------------------ */
/* Names the row in cft_support_rows, so if a later parcel deletes it
 * outright the coverage check in check_dropin.c says this case is
 * orphaned. The row is still there and still refuses two things: any
 * inactive particle on the SUBPROCESS path, which has no way to carry
 * one, and N_active greater than N on either. */
static void p_active(struct reb_simulation *r){
    r->N_active = 1;
    r->testparticle_type = 1;
    r->testparticle_hidewarnings = 1;
}

/* ------------------------------------------------------------------ */
/* The topic, in the order it is run                                   */
/* ------------------------------------------------------------------ */
void cases_pairs(void){
    printf("\nwhich pairs are computed: r->N_active and r->testparticle_type\n");

    /* First the controls, so that what follows means something. */
    printf("\nthe controls - does the skip reach the arithmetic at all:\n");
    case_active_changes_run("massive inactive bodies, N_active",
                            semi, PAIRS_N, 0.3, 1e-9, 300, 3);
    case_tptype_changes_run("massive inactive bodies, tp_type",
                            semi, PAIRS_N, 0.3, 1e-9, 300, 3);
    case_massless_is_inert("massless test particles",
                           dust, PAIRS_N, 0.3, 1e-9, 300, 3);
    case_clobbers_ignore_terms("gravity_ignore_terms under ias15",
                               semi, PAIRS_N, 0.3, 1e-9, 60);

    /* Then gravity alone, one configuration at a time. This is where
     * ignore_terms is actually exercised: nothing downstream of IAS15
     * can see it. */
    printf("\ngravity alone, against reb_simulation_update_acceleration():\n");
    if (ias15_engine_format() != CFT_FP64){
        printf("  FAIL the engine is not at binary64; the oracle compares "
               "wide gravity against REBOUND's doubles\n");
        failures++;
    }else{
        /* the settings a user would write, named and with their
         * distance from the unrestricted answer beside them */
        for (int ign = 0; ign < IGN_N; ign++){
            char label[96];
            snprintf(label, sizeof label, "all active, ignore_terms = %s", ign_name[ign]);
            case_gravity(label, semi, PAIRS_N, (size_t)-1, 0, ign);
        }
        case_gravity("N_active 3, tp_type 0", semi, PAIRS_N, 3, 0, 0);
        case_gravity("N_active 3, tp_type 1", semi, PAIRS_N, 3, 1, 0);
        case_gravity("N_active 3, tp_type 1, INVOLVING_0", semi, PAIRS_N, 3, 1, 2);
        /* N_active = 0 is legal and means nothing is a source: both
         * loops are empty and every acceleration is +0. */
        case_gravity("N_active 0 (nothing is a source)", semi, PAIRS_N, 0, 0, 0);
        case_gravity("N_active 6 of 6 (the no-op)", semi, PAIRS_N, PAIRS_N, 0, 0);
        case_gravity("massless dust, N_active 3", dust, PAIRS_N, 3, 0, 0);
        case_gravity("massless dust, N_active 3, tp_type 1", dust, PAIRS_N, 3, 1, 0);
        /* and then all of them */
        case_gravity_sweep("every configuration", semi, PAIRS_N);
    }

    /* Then whole integrations against REBOUND's own ias15. */
    printf("\ntest particles, bit for bit against REBOUND's ias15:\n");
    case_pairs("massless dust, N_active 3",        dust, PAIRS_N, 0.3,  1e-9, 400, 3, 0);
    case_pairs("massless dust, N_active 1",        dust, PAIRS_N, 0.3,  1e-9, 400, 1, 0);
    case_pairs("massive semi-active, N_active 3",  semi, PAIRS_N, 0.3,  1e-9, 400, 3, 0);
    case_pairs("massive semi-active, N_active 1",  semi, PAIRS_N, 0.3,  1e-9, 400, 1, 0);
    case_pairs("semi-active, tp_type 1",           semi, PAIRS_N, 0.3,  1e-9, 400, 3, 1);
    case_pairs("semi-active, N_active 1, tp 1",    semi, PAIRS_N, 0.3,  1e-9, 400, 1, 1);
    /* fixed step: no adaptive controller to hide a difference in, and
     * REBOUND's predictor-corrector runs on exactly the accelerations
     * the pair list produced */
    case_pairs("semi-active, fixed step",          semi, PAIRS_N, 0.05, 0.0,  400, 3, 0);
    case_pairs("semi-active, fixed step, tp 1",    semi, PAIRS_N, 0.05, 0.0,  400, 3, 1);
    /* two bodies of which one is a test particle: the smallest problem
     * in which the active-active loop is empty and only the test-active
     * loop runs */
    case_pairs("one active, one test",             semi, 2,       0.05, 1e-9, 400, 1, 0);
    case_pairs("one active, one test, tp 1",       semi, 2,       0.05, 1e-9, 400, 1, 1);

    printf("\nand what must no longer be refused:\n");
    accepted("test_particles", "test particles (r->N_active)", p_active);
    printf("\n");
}
