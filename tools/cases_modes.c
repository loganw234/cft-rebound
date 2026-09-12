/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * IAS15's OTHER two step criteria: adaptive_mode 0 (INDIVIDUAL) and
 * 1 (GLOBAL).
 *
 * PRS23 (2) and AARSETH85 (3) share every sum - a0 squared and four
 * weighted sums of b - and differ in one expression, so AARSETH85 was a
 * branch. These two share none of those sums: they estimate the error
 * from the LAST term of the series against the acceleration at the end
 * of the sequence, and they change the CORRECTOR's convergence test as
 * well as the step control. That is why they were deferred rather than
 * added as one line, and it is also what this file has to prove was
 * done rather than ignored.
 *
 * Two assertions, and the second is the one that matters:
 *
 *   1. Equivalence. Whatever step sequence REBOUND's own ias15 chooses
 *      in each mode, the port chooses the same bits - the same assertion
 *      every other case in this gate makes.
 *
 *   2. The control. Each mode must choose a DIFFERENT step sequence from
 *      PRS23 on the same problem. Without it, a mode the port silently
 *      ignored would leave every case in (1) passing PRS23 against
 *      PRS23: a gate that cannot fail. AARSETH85 needed exactly this
 *      control (tools/cases_core.c says so at its own block).
 *
 * The control is for the STEP CONTROL half. It cannot speak for the
 * corrector half, because a port that took INDIVIDUAL's step estimate
 * and kept PRS23's corrector error would still choose a different
 * sequence and would still pass every one of these cases but one. That
 * one is named where it appears; it was found by measuring, and the
 * comment above it says what was measured.
 *
 * PARCELS.md, "P3 - The other two step criteria".
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"

static const char *mode_name(int m){
    switch (m){
        case 0:  return "INDIVIDUAL";
        case 1:  return "GLOBAL";
        case 2:  return "PRS23";
        case 3:  return "AARSETH85";
        default: return "?";
    }
}

/* ------------------------------------------------------------------ */
/* (1) equivalence, per mode                                            */
/* ------------------------------------------------------------------ */
static void case_mode_steps(const char *label, const struct body *bs, size_t n,
                            double dt, double epsilon, size_t steps, int mode){
    int save = adaptive_mode;
    adaptive_mode = mode;
    printf("%s (N = %zu, dt = %g, epsilon = %g, %zu steps, adaptive_mode %d %s)\n",
           label, n, dt, epsilon, light_ran(steps), mode, mode_name(mode));
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    reb_simulation_steps(ra, steps);
    reb_simulation_steps(rb, steps);
    if (verbose) printf("    t = %.17g, dt = %.17g\n", ra->t, ra->dt);
    compare(ra, rb, label);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
    adaptive_mode = save;
}

/* ------------------------------------------------------------------ */
/* (2) the control: does the mode actually change the step sequence?    */
/* ------------------------------------------------------------------ */
/* One step at a time, recording the step that was taken. REBOUND's
 * reb_simulation_steps(r, k) is a per-step loop plus a synchronize IAS15
 * does not implement, so the trajectory is the same as one call of k -
 * ref/ias15_ref.c relies on the same fact for its exact clock. */
#define SEQ_MAX 512
static void dt_sequence(const struct body *bs, size_t n, double dt, double epsilon,
                        int mode, int use_cft, size_t steps, double *out, double *t_end){
    int save = adaptive_mode;
    adaptive_mode = mode;
    struct reb_simulation *r = build(bs, n, dt, epsilon, use_cft);
    adaptive_mode = save;
    for (size_t i = 0; i < steps; i++){
        reb_simulation_steps(r, 1);
        out[i] = r->dt_last_done;
    }
    *t_end = r->t;
    reb_simulation_free(r);
}

/* Returns the number of steps whose size differs, and reports the first
 * one. Bit patterns, like every other comparison in this gate: two steps
 * that agree to every printed digit and differ in the last bit ARE a
 * different sequence, and saying so is the point. */
static size_t seq_differs(const double *a, const double *b, size_t steps,
                          size_t *first, double *fa, double *fb){
    size_t ndiff = 0;
    *first = steps; *fa = 0; *fb = 0;
    for (size_t i = 0; i < steps; i++)
        if (bits_differ(a[i], b[i])){
            if (!ndiff){ *first = i; *fa = a[i]; *fb = b[i]; }
            ndiff++;
        }
    return ndiff;
}

static void case_control(const char *label, const struct body *bs, size_t n,
                         double dt, double epsilon, size_t steps, int mode){
    static double prs_ref[SEQ_MAX], prs_cft[SEQ_MAX], m_ref[SEQ_MAX], m_cft[SEQ_MAX];
    double t_prs_ref, t_prs_cft, t_m_ref, t_m_cft;
    if (steps > SEQ_MAX){ printf("  FAIL %s: raise SEQ_MAX\n", label); failures++; return; }

    printf("%s: adaptive_mode %d (%s) against PRS23, %zu steps each\n",
           label, mode, mode_name(mode), steps);
    dt_sequence(bs, n, dt, epsilon, 2,    0, steps, prs_ref, &t_prs_ref);
    dt_sequence(bs, n, dt, epsilon, 2,    1, steps, prs_cft, &t_prs_cft);
    dt_sequence(bs, n, dt, epsilon, mode, 0, steps, m_ref,   &t_m_ref);
    dt_sequence(bs, n, dt, epsilon, mode, 1, steps, m_cft,   &t_m_cft);

    size_t first; double fa, fb;

    /* (a) REBOUND itself chooses differently, so the problem is one on
     *     which the two criteria can be told apart at all. A control on
     *     a problem where they happen to agree would prove nothing. */
    size_t nref = seq_differs(prs_ref, m_ref, steps, &first, &fa, &fb);
    if (!nref){
        printf("  FAIL %s: REBOUND's own ias15 takes the SAME %zu steps under PRS23 and\n"
               "       under %s, so this problem cannot tell the two criteria apart and\n"
               "       the equivalence cases above are PRS23 against PRS23\n",
               label, steps, mode_name(mode));
        failures++;
        return;
    }
    printf("  ref  REBOUND: %zu of %zu steps differ; first at step %zu, "
           "PRS23 %.17g vs %s %.17g (ratio %.6g)\n",
           nref, steps, first, fa, mode_name(mode), fb, fb / fa);

    /* (b) and so does the port. This is the assertion: a port that
     *     ignored the mode would take PRS23's steps here and ndiff
     *     would be 0. */
    size_t ncft = seq_differs(prs_cft, m_cft, steps, &first, &fa, &fb);
    if (!ncft){
        printf("  FAIL %s: the port takes the SAME %zu steps under PRS23 and under %s -\n"
               "       the mode was ignored, and every equivalence case above compared\n"
               "       PRS23 against PRS23\n", label, steps, mode_name(mode));
        failures++;
        return;
    }
    printf("  ok   port:    %zu of %zu steps differ; first at step %zu, "
           "PRS23 %.17g vs %s %.17g (ratio %.6g)\n",
           ncft, steps, first, fa, mode_name(mode), fb, fb / fa);
    printf("       after %zu steps: PRS23 reaches t = %.17g, %s reaches t = %.17g\n",
           steps, t_prs_cft, mode_name(mode), t_m_cft);

    /* (c) and the port's sequence is REBOUND's own, step for step. The
     *     equivalence cases compare the END state; this compares every
     *     step on the way, which is what "the same step sequence" means
     *     and what the two above have just been measured against. */
    size_t nseq = seq_differs(m_ref, m_cft, steps, &first, &fa, &fb);
    if (nseq){
        printf("  FAIL %s: under %s the port's step sequence differs from REBOUND's at\n"
               "       step %zu: ias15 %.17g, cft %.17g\n",
               label, mode_name(mode), first, fa, fb);
        failures++;
        return;
    }
    if (bits_differ(t_m_ref, t_m_cft)){
        printf("  FAIL %s: same steps, different clock: ias15 %.17g, cft %.17g\n",
               label, t_m_ref, t_m_cft);
        failures++;
        return;
    }
    printf("  ok   %s: all %zu steps bit for bit as REBOUND's ias15 chose them\n",
           label, steps);
}

/* ------------------------------------------------------------------ */
void cases_modes(void){
    printf("\nwith IAS15's other two step criteria, adaptive_mode 0 and 1:\n");

    /* The equivalence set, on the same three problems the PRS23 and
     * AARSETH85 blocks use. kepler and pythagorean are planar, so az is
     * +0 for all time and INDIVIDUAL divides by it on every z lane every
     * step - REBOUND's estimate does that too, isnormal() drops the
     * result, and this is the case that says the port drops it the same
     * way. five is out of the plane and has a body that starts at the
     * origin, which is where GLOBAL's own divide by |x|^2 is closest to
     * the edge. */
    for (int mode = 0; mode <= 1; mode++){
        case_mode_steps("kepler",       kepler, 2,      0.05, 1e-9, 400, mode);
        case_mode_steps("pythagorean",  pythagorean, 3, 0.01, 1e-9, 400, mode);
        case_mode_steps("five bodies",  five, 5,        0.5,  1e-9, 300, mode);
    }

    /* Rejections. dt = 8 on an orbit of period 2 pi is the rejection
     * configuration cases_core.c uses. */
    for (int mode = 0; mode <= 1; mode++)
        case_mode_steps("kepler, rejected steps", kepler, 2, 8.0, 1e-9, 200, mode);

    /* THE CORRECTOR'S OWN ERROR, which is the half of INDIVIDUAL that
     * is not in the step control at all (pc_error in src/ias15_cft.c).
     * It decides when the predictor-corrector loop stops, so it only
     * reaches the recorded state when a pass more or less changes the b
     * coefficients - and on a converged step it does not. Measured: with
     * INDIVIDUAL's corrector error disabled but its step control left
     * intact, every case above still passes, while the corrector's pass
     * count moves - 400 steps of data/problems/kepler.txt at dt = 0.01
     * run a mean of 2.595 passes against 2.570 and a peak of 5 against
     * 4. This case is the one that carries that through to the state,
     * and it was found by sweeping for it rather than assumed:
     *
     *   five bodies from dt = 500, which is half the outermost body's
     *   period, so the first attempt is rejected hard and the corrector
     *   is still far from converged when it stops. There, disabling the
     *   branch changes the RECORD.
     *
     * mode 1 runs the same configuration because GLOBAL shares PRS23's
     * corrector error and so must NOT move here - it is the control on
     * the control. */
    for (int mode = 0; mode <= 1; mode++)
        case_mode_steps("five bodies from dt = 500 (corrector far from converged)",
                        five, 5, 500.0, 1e-9, 300, mode);

    /* Two integrate() calls in a row, which is the shim path that resets
     * dt_last_done and so changes what a rejection predicts - run here
     * under a corrector whose error is a different quantity. */
    for (int mode = 0; mode <= 1; mode++){
        int save = adaptive_mode;
        adaptive_mode = mode;
        printf("kepler, two integrate calls (adaptive_mode %d %s)\n", mode, mode_name(mode));
        struct reb_simulation *ra = build(kepler, 2, 8.0, 1e-9, 0);
        struct reb_simulation *rb = build(kepler, 2, 8.0, 1e-9, 1);
        reb_simulation_integrate(ra, 20.0); reb_simulation_integrate(ra, 53.0);
        reb_simulation_integrate(rb, 20.0); reb_simulation_integrate(rb, 53.0);
        compare(ra, rb, "kepler, two integrate calls");
        reb_simulation_free(ra); reb_simulation_free(rb);
        adaptive_mode = save;
    }

    printf("\nthe control: each mode must choose a different step sequence from PRS23\n");
    case_control("kepler",      kepler, 2,      0.05, 1e-9, 200, 0);
    case_control("kepler",      kepler, 2,      0.05, 1e-9, 200, 1);
    case_control("pythagorean", pythagorean, 3, 0.01, 1e-9, 200, 0);
    case_control("pythagorean", pythagorean, 3, 0.01, 1e-9, 200, 1);
    case_control("five bodies", five, 5,        0.5,  1e-9, 200, 0);
    case_control("five bodies", five, 5,        0.5,  1e-9, 200, 1);

    /* And the two new modes must differ from EACH OTHER, or one of them
     * is the other under a different number. */
    {
        static double s0[SEQ_MAX], s1[SEQ_MAX];
        double t0, t1; size_t first; double fa, fb;
        const size_t steps = 200;
        dt_sequence(five, 5, 0.5, 1e-9, 0, 1, steps, s0, &t0);
        dt_sequence(five, 5, 0.5, 1e-9, 1, 1, steps, s1, &t1);
        size_t nd = seq_differs(s0, s1, steps, &first, &fa, &fb);
        if (!nd){
            printf("  FAIL five bodies: INDIVIDUAL and GLOBAL take the same %zu steps;\n"
                   "       one of the two modes is not reaching its own branch\n", steps);
            failures++;
        }else{
            printf("  ok   five bodies: INDIVIDUAL and GLOBAL differ too - %zu of %zu steps,\n"
                   "       first at step %zu, %.17g vs %.17g (ratio %.6g)\n",
                   nd, steps, first, fa, fb, fb / fa);
        }
    }
    printf("\n");
}
