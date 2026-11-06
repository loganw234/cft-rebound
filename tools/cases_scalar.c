/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * THE SCALAR ROUTE IN THE LIBRARY HALF (ROUND3.md, P5).
 *
 * A constant operand goes to the library as one element, through the
 * door's vrun_k() (src/ias15_cft.c, vmulk() and vfmak()). The family is
 * opt-in: unforced it is off, by=opt_in, even where the device publishes
 * the bit, because no card has measured it. So every case in check_dropin
 * holds the sites to REBOUND's own ias15 on the route this process opened
 * with: today's in the suite's own run, and the new one in
 * tools/check_scalar.py --library, which runs this whole suite with the
 * route forced on. That pair is the REBOUND gate covering the route. This
 * topic adds two things it cannot see.
 *
 * 1. THE ROUTE, READ. ias15_engine_route() is the library half's view of
 *    the program's `# route` line (V0, 2026-10-31). The route must be in
 *    the state this process's environment requires (expected(), below),
 *    and its counters must agree: on, calls and no fallbacks; off, the
 *    reverse. Without this, a library half that ignored the forced route
 *    would pass every case in the suite (P1's M4, lead.md 16:55 onward).
 * 2. THE SITES, COUNTED. A route line says calls went through the door,
 *    not that the right sites sent them. So each run here derives what
 *    the sites must have sent from the integrator's shape, seven nodes and
 *    seven levels, and from the steps it took:
 *      - a fixed-step run gives the pass count, and it must be whole;
 *      - a do-nothing velocity-dependent routine leaves every bit where it
 *        was and adds predict_velocities(), the library half's four sites
 *        and the program's none. So it must add exactly its calls for each
 *        of those same passes. That holds the only sites no standalone
 *        gate reaches;
 *      - the FMA form takes the other half of predict_positions() and
 *        correct(), and step_attempt()'s node and end-of-step sites. It
 *        has no REBOUND to equal, so its control is that it sends other
 *        calls than REBOUND's form of the same run.
 *
 * --light. Every case here scales: each takes a step count and reads the
 * steps actually taken back from REBOUND, r->steps_done, so a scaled run
 * asserts what it ran. No count here names a time or an event.
 *
 * Stated, not closed: the FMA form in the library half is held only by its
 * counts. Nothing compares it across the two routes, because the engine
 * opens once a process. tools/check_scalar.py holds the same sites across
 * the routes in the standalone program, at every format.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dropin_cases.h"
#include "ias15_engine.h"

/* The route's state, found by NAME among ias15_engine_route()'s rows. */
struct sroute { int found, on; const char *by; unsigned long long calls, fallbacks; };

static struct sroute scalar_route(void){
    struct sroute g;
    struct ias15_engine_route r;
    memset(&g, 0, sizeof g);
    for (int k = 0; ias15_engine_route(k, &r) == 0; k++)
        if (!strcmp(r.name, "scalar")){
            g.found = 1; g.on = r.on; g.by = r.by; g.calls = r.calls; g.fallbacks = r.fallbacks;
        }
    return g;
}

/* An item of a comma list, matched whole: the door's own reading of
 * CFT_REBOUND_HIDE_FEATURES (src/ias15_cft.c, engine_routes()). */
static int list_names(const char *list, const char *name){
    const size_t n = strlen(name);
    for (const char *p = list; p && *p; ){
        const size_t l = strcspn(p, ",");
        if (l == n && !strncmp(p, name, n)) return 1;
        p += l;
        if (*p == ',') p++;
    }
    return 0;
}

/* The state the door must have opened the route in, from the environment
 * alone, as tools/cases_gravity.c derives its two:
 *   - CFT_REBOUND_ROUTE_SCALAR on or off: that state, by=forced;
 *   - CFT_REBOUND_HIDE_FEATURES naming it: off by=hidden;
 *   - the software backend, with neither: off by=opt_in and nothing else.
 *     It publishes every bit, and the route is opt-in;
 *   - a card, with neither: off by=opt_in, or off by=absent on an image
 *     without the bit.
 * So unforced the route is never on, and `on by=device` is never right:
 * that would be auto taking an opt-in route. Forced on is the only state
 * in which this topic can see a library half silently on today's route.
 * Any other value of the variable is refused at open, and "auto" or empty
 * is the unforced case. */
static void expected(int *on, const char **by, const char **by_card_off){
    const char *want = getenv("CFT_REBOUND_ROUTE_SCALAR"), *hide = getenv("CFT_REBOUND_HIDE_FEATURES"),
               *art = getenv("CFT_REBOUND_ARTIFACT");
    *by_card_off = NULL;
    if (want && (!strcmp(want, "on") || !strcmp(want, "off"))){ *on = !strcmp(want, "on"); *by = "forced"; }
    else if (hide && list_names(hide, "scalar")){ *on = 0; *by = "hidden"; }
    else if (!art || !art[0]){ *on = 0; *by = "opt_in"; }
    else { *on = 0; *by = "opt_in"; *by_card_off = "absent"; }
}

/* The route held: in the state the environment requires, and its counters
 * agreeing over the stretch of work that moved them by dc and df. One line
 * either way, with "route scalar=STATE by=BY calls=N fallbacks=N" in it
 * word for word, which tools/check_scalar.py --library reads. */
static void route_held(const char *label, const struct sroute *g, unsigned long long dc, unsigned long long df){
    const char *scope = label[0] ? " in this case" : "";
    int want_on;
    const char *want_by, *card_off;
    if (!g->found){
        printf("  FAIL %sias15_engine_route() has no route named scalar\n", label);
        failures++;
        return;
    }
    expected(&want_on, &want_by, &card_off);
    const int state = (g->on == want_on && !strcmp(g->by, want_by)) ||
                      (card_off && !g->on && !strcmp(g->by, card_off));
    const int counts = g->on ? (dc > 0 && df == 0) : (df > 0 && dc == 0);
    if (state && counts)
        printf("  ok   %sroute scalar=%s by=%s calls=%llu fallbacks=%llu%s, as this environment requires\n",
               label, g->on ? "on" : "off", g->by, dc, df, scope);
    else{
        printf("  FAIL %sroute scalar=%s by=%s calls=%llu fallbacks=%llu%s, where this environment requires %s by=%s%s%s"
               ", with the counter of the route taken above zero and the other zero\n", label,
               g->on ? "on" : "off", g->by, dc, df, scope, want_on ? "on" : "off", want_by,
               card_off ? " or off by=" : "", card_off ? card_off : "");
        failures++;
    }
}

/* What one stretch of work sent through the door, the route held. */
static unsigned long long sent(const char *label, const struct sroute *before){
    char l[176];
    const struct sroute g = scalar_route();
    const unsigned long long dc = g.calls - before->calls, df = g.fallbacks - before->fallbacks;
    snprintf(l, sizeof l, "%s: ", label);
    route_held(l, &g, dc, df);
    return dc + df;
}

/* The integrator's shape: REBOUND's seven nodes, and the seven levels of
 * its b polynomial. Every count below is derived from these two. */
enum { NODES = 7, LEVELS = 7 };

/* Calls through the door a corrector pass, predict_positions() and
 * correct() over the nodes. Node n of the corrector has n - 1 of KC, and
 * under the FMA form n of KRINV. REBOUND's form of the predictor has KNUM
 * at each level, and H at each level and twice more, for dt*h/2 and dt*h. */
static unsigned long long door_a_pass(int fma){
    unsigned long long c = 0;
    for (int n = 1; n <= NODES; n++)
        c += fma ? (unsigned long long)(LEVELS + n + n - 1) : (unsigned long long)(LEVELS + LEVELS + 2 + n - 1);
    return c;
}

/* Calls a step accepted at its first attempt, outside the passes: g from
 * b (KD) and the next step's prediction (KBIN), one triangle of the levels
 * each, and under the FMA form both of a node's dt products, then the
 * end-of-step update's reciprocals (KPOSR with a0's, and KVELR). */
static unsigned long long door_a_step(int fma){
    const unsigned long long tri = LEVELS * (LEVELS - 1) / 2;
    return 2 * tri + (fma ? 2 * NODES + (LEVELS + 1) + LEVELS : 0);
}

/* predict_velocities() a pass: KVNUM at each level, and H at each level
 * and once more, for dt*h. */
static unsigned long long door_velocities_a_pass(void){ return (unsigned long long)NODES * (LEVELS + LEVELS + 1); }

/* The passes those calls imply: they must divide whole, and a step's
 * passes lie between one and max_iter. Zero if they do not. */
static unsigned long long passes_from(const char *label, unsigned long long calls, unsigned long long steps,
                                      int fma, int max_iter){
    const unsigned long long fixed = steps * door_a_step(fma), pass = door_a_pass(fma);
    const unsigned long long p = calls > fixed ? (calls - fixed) / pass : 0;
    if (!steps || calls <= fixed || (calls - fixed) % pass || p < steps || p > steps * (unsigned long long)max_iter){
        printf("  FAIL %s: %llu calls through the door in %llu steps, which is not %llu a step and a whole number of "
               "passes of %llu, between 1 and %d a step\n", label, calls, steps, door_a_step(fma), pass, max_iter);
        failures++;
        return 0;
    }
    printf("  ok   %s: %llu calls through the door are %llu steps of %llu and %llu passes of %llu\n", label, calls, steps,
           door_a_step(fma), p, pass);
    return p;
}

static void f_nop(struct reb_simulation *r){ (void)r; }

/* REBOUND's form on a fixed step, against REBOUND, and then the same run
 * with a do-nothing velocity-dependent routine. */
static void case_sites(const char *label, const struct body *bs, size_t n, double dt, size_t steps){
    char l2[160];
    struct sroute before;
    printf("%s (%zu bodies, dt = %g, fixed, %zu steps)\n", label, n, dt, light_ran(steps));
    struct reb_simulation *ra = build(bs, n, dt, 0.0, 0);
    struct reb_simulation *rb = build(bs, n, dt, 0.0, 1);
    reb_simulation_steps(ra, steps);
    before = scalar_route();
    reb_simulation_steps(rb, steps);
    compare(ra, rb, label);
    const unsigned long long plain = sent(label, &before);
    const unsigned long long taken = (unsigned long long)rb->steps_done;
    /* max_iter's 0 is resolved at the first step and written back */
    const unsigned long long passes = passes_from(label, plain, taken, 0, cft_ias15_get_state(rb)->max_iter);
    reb_simulation_free(ra);
    reb_simulation_free(rb);

    snprintf(l2, sizeof l2, "%s, a do-nothing velocity-dependent routine", label);
    ra = build(bs, n, dt, 0.0, 0);
    rb = build(bs, n, dt, 0.0, 1);
    ra->additional_forces = f_nop; ra->force_is_velocity_dependent = 1;
    rb->additional_forces = f_nop; rb->force_is_velocity_dependent = 1;
    reb_simulation_steps(ra, steps);
    before = scalar_route();
    reb_simulation_steps(rb, steps);
    compare(ra, rb, l2);
    const unsigned long long with = sent(l2, &before);
    const unsigned long long added = with - plain, want = passes * door_velocities_a_pass();
    if (!passes || with <= plain || added != want){
        printf("  FAIL %s: predict_velocities() added %llu calls, and %llu passes of %llu are %llu\n", l2,
               with > plain ? added : 0, passes, door_velocities_a_pass(), want);
        failures++;
    }else
        printf("  ok   %s: predict_velocities() added %llu calls, %llu passes of %llu\n", l2, added, passes,
               door_velocities_a_pass());
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* The FMA form. Its control is REBOUND's form of the same run. The two can
 * agree to the last bit at binary64 for many steps (measured on Kepler at
 * dt = 0.05: every value of a 40-step sample), so bits cannot say the form
 * took. The calls can. A form that never took sends what REBOUND's form
 * sends, exactly, because it is the same run. */
static void case_fma(const char *label, const struct body *bs, size_t n, double dt, size_t steps){
    struct sroute before;
    printf("%s (%zu bodies, dt = %g, fixed, %zu steps)\n", label, n, dt, light_ran(steps));
    struct reb_simulation *rb = build(bs, n, dt, 0.0, 1);
    before = scalar_route();
    reb_simulation_steps(rb, steps);
    const unsigned long long plain = sent(label, &before);
    reb_simulation_free(rb);
    rb = build(bs, n, dt, 0.0, 1);
    cft_ias15_get_state(rb)->arith_fma = 1;
    before = scalar_route();
    reb_simulation_steps(rb, steps);
    const unsigned long long fma = sent(label, &before);
    passes_from(label, fma, (unsigned long long)rb->steps_done, 1, cft_ias15_get_state(rb)->max_iter);
    if (fma == plain){
        printf("  FAIL %s: %llu calls, the same as REBOUND's form, so arith_fma never reached the arithmetic\n", label, fma);
        failures++;
    }else
        printf("  ok   %s: the form took, %llu calls against REBOUND's form's %llu\n", label, fma, plain);
    reb_simulation_free(rb);
}

void cases_scalar(void){
    printf("\nthe scalar route in the library half: a constant as one element\n");
    printf("the route every case so far took (ias15_engine_route()):\n");
    {   const struct sroute g = scalar_route();
        route_held("", &g, g.calls, g.fallbacks); }
    printf("\nthe sites, counted:\n");
    case_sites("two bodies", kepler, 2, 0.05, 100);
    case_sites("three bodies", pythagorean, 3, 0.01, 100);
    case_fma("two bodies, the FMA form", kepler, 2, 0.05, 100);
    printf("\n");
}
