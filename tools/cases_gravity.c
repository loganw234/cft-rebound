/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * GRAVITY'S ROUTES IN THE LIBRARY HALF (ROUND3.md, P1).
 *
 * Gravity runs as sequencer programs now, in both halves. Its two switches
 * are the door's: `indexed` (the accumulate fold) and `pair_gather` (the
 * pair program's six inputs), each decided once, at open, by the device's
 * bit or a CFT_REBOUND_ROUTE_* variable (src/ias15_cft.c, engine_routes()).
 * So every case in check_dropin runs gravity on the routes this process
 * opened with - on the software backend, unforced, both new routes - and
 * holds it to REBOUND's own ias15 bit for bit. That is the REBOUND gate
 * covering the routes. This topic adds three things it cannot see.
 *
 * 1. THE ROUTES, READ. ias15_engine_route() is the library half's view of
 *    the program's `# route` lines (V0, 2026-10-31). A route that is on
 *    must have carried calls and no fallbacks, and one that is off the
 *    reverse. On the software backend with no door variable set, both
 *    must be on `by=device`, because it publishes every bit. Without
 *    this, a library half that never took the new route would pass every
 *    case in the suite. Read here for the whole process so far, and
 *    around each case below.
 * 2. THE TABLES, REBUILT MID-RUN. Gravity's index tables are built for a
 *    pair list and a set of rows (build_scatter()), which change in two
 *    places: a body removed, through ias15_engine_set_bodies(), and
 *    N_active changed, through ias15_engine_set_active(). A table that
 *    outlived its N either reads past its source, which libcft and the
 *    host gather both refuse by name, or reads the wrong partner, which is
 *    silent. The N_active change is the one that can be silent: the pool
 *    keeps its size, so a stale table's indices all stay in range. Both
 *    happen mid-run here, against REBOUND doing the same, each with its
 *    physics counter asserted beside the comparison (CLAUDE.md, "Three
 *    ways this suite has lied"):
 *      - a merge that removed one body, "4 particles became 3";
 *      - an N_active change that measurably changed the run, against the
 *        same run unchanged.
 * 3. TEST PARTICLES THAT PULL BACK, on the routes, with the route counters
 *    asserted around them. Massive, so that the restriction is not
 *    bit-inert (tools/cases_pairs.c, point 1).
 *
 * --light. The merge case opts out by name, (reb_simulation_steps)(...):
 * its merge fires late, and a tenth of its steps is a run where nothing
 * happens. The N_active case scales, and its control is what says the
 * change still reached the arithmetic at the scaled length.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dropin_cases.h"
#include "ias15_engine.h"

/* A route's state, found by NAME among ias15_engine_route()'s rows. */
struct groute { int found, on; const char *by; unsigned long long calls, fallbacks; };

static struct groute route_named(const char *name){
    struct groute g;
    struct ias15_engine_route r;
    memset(&g, 0, sizeof g);
    for (int k = 0; ias15_engine_route(k, &r) == 0; k++)
        if (!strcmp(r.name, name)){
            g.found = 1; g.on = r.on; g.by = r.by; g.calls = r.calls; g.fallbacks = r.fallbacks;
        }
    return g;
}

static const char *const groutes[2] = { "indexed", "pair_gather" };

/* What one stretch of work moved: on, calls and no fallbacks; off, the
 * reverse. Prints the two route lines, and fails a route that carried
 * nothing or carried it the wrong way. */
static void routes_moved(const char *label, const struct groute before[2]){
    for (int k = 0; k < 2; k++){
        const struct groute g = route_named(groutes[k]);
        const unsigned long long dc = g.calls - before[k].calls, df = g.fallbacks - before[k].fallbacks;
        const int right = g.on ? (dc > 0 && df == 0) : (df > 0 && dc == 0);
        printf("  %s %s: route %s=%s by=%s carried %llu calls and %llu fallbacks\n", right ? "ok  " : "FAIL",
               label, groutes[k], g.on ? "on" : "off", g.by, dc, df);
        if (!right) failures++;
    }
}

static void routes_now(struct groute out[2]){
    for (int k = 0; k < 2; k++) out[k] = route_named(groutes[k]);
}

/* 1. The routes, for the whole process so far: every case before this one
 *    ran through them. */
static void case_routes_so_far(void){
    const char *art = getenv("CFT_REBOUND_ARTIFACT");
    const int software = !art || !art[0];
    const int unforced = !getenv("CFT_REBOUND_ROUTE_INDEXED") && !getenv("CFT_REBOUND_ROUTE_PAIR_GATHER") &&
                         !getenv("CFT_REBOUND_HIDE_FEATURES");
    printf("the routes every case so far took (ias15_engine_route()):\n");
    for (int k = 0; k < 2; k++){
        const struct groute g = route_named(groutes[k]);
        int right;
        if (!g.found){
            printf("  FAIL ias15_engine_route() has no route named %s\n", groutes[k]);
            failures++;
            continue;
        }
        right = g.on ? (g.calls > 0 && g.fallbacks == 0) : (g.fallbacks > 0 && g.calls == 0);
        if (software && unforced) right &= g.on && !strcmp(g.by, "device");
        printf("  %s route %s=%s by=%s calls=%llu fallbacks=%llu%s\n", right ? "ok  " : "FAIL", groutes[k],
               g.on ? "on" : "off", g.by, g.calls, g.fallbacks,
               software && unforced ? " (the software backend, unforced: on by=device is required)" : "");
        if (!right) failures++;
    }
}

/* 2a. A merge mid-run: ias15_engine_set_bodies() rebuilds the tables. A
 *     star and two small bodies that meet late - cases_core.c's
 *     case_collision, whose tangential velocities keep the remnant in
 *     orbit - and a fourth far out that keeps a system of three afterwards,
 *     so the rebuilt tables are tables of three pairs, not one. The fourth
 *     is light and far: at 1e-4 and r = 3 it moved the encounter off, and
 *     REBOUND kept all four (measured). */
static void case_merge(const char *label, double dt, double epsilon, size_t steps){
    static const struct body meet[4] = {
        { 1.0,     0.0,  0.0, 0.0,   0.0,  0.0,    0.0 },
        { 1.0e-3,  1.0,  0.0, 0.0,   0.0,  1.3,    0.0 },
        { 1.0e-3,  1.02, 0.0, 0.0,   0.0,  0.7,    0.0 },
        { 1.0e-10, -20.0, 0.0, 0.25, 0.0, -0.2236, 0.0 },
    };
    struct groute before[2];
    printf("%s (4 bodies, dt = %g, epsilon = %g, %zu steps unscaled, merge on contact)\n", label, dt, epsilon, steps);
    struct reb_simulation *ra = build(meet, 4, dt, epsilon, 0);
    struct reb_simulation *rb = build(meet, 4, dt, epsilon, 1);
    for (int i = 0; i < 2; i++){
        struct reb_simulation *r = i ? rb : ra;
        r->collision = REB_COLLISION_DIRECT;
        r->collision_resolve = reb_collision_resolve_merge;
        for (size_t k = 0; k < r->N; k++) r->particles[k].r = 0.02;
    }
    const size_t n_before = ra->N;
    /* NOT scaled by --light: the parentheses suppress the macro. The
     * merge fires late, and a tenth of the steps never reaches it. */
    (reb_simulation_steps)(ra, steps);
    routes_now(before);
    (reb_simulation_steps)(rb, steps);
    if (ra->N >= n_before || ra->N < 3){
        printf("  FAIL %s: REBOUND ended with %zu of %zu particles; the case needs one merge that leaves\n"
               "       three, or it tests no rebuild of the tables\n", label, ra->N, n_before);
        failures++;
        reb_simulation_free(ra); reb_simulation_free(rb);
        return;
    }
    printf("  (a merge removed one: %zu particles became %zu, on both sides: %zu and %zu)\n",
           n_before, ra->N, ra->N, rb->N);
    compare(ra, rb, label);
    routes_moved(label, before);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

/* Six bodies, every one with mass, so that an N_active restriction drops
 * real contributions (tools/cases_pairs.c's `semi`, restated here because
 * that file's is static). */
static const struct body grav6[6] = {
    { 1.0,       0.0,   0.0,  0.0,     0.0,    0.0,    0.0     },
    { 9.55e-4,   5.2,   0.0,  0.0625,  0.0,    0.438,  0.00125 },
    { 2.86e-4,  -9.5,   0.0, -0.125,   0.0,   -0.324,  0.0025  },
    { 1.0e-5,    2.0,   0.0,  0.0,     0.0,    0.70,   0.0     },
    { 2.0e-5,    0.0,   3.1,  0.05,   -0.56,   0.0,    0.002   },
    { 3.0e-5,   -7.3,   0.0, -0.2,     0.0,   -0.37,  -0.001   },
};

/* 2b. N_active changed mid-run, twice: every body active, then N_active 3,
 *     then every body again. Each change goes through
 *     ias15_engine_set_active() and rebuilds the rows and the fold's
 *     table, with the pool unchanged in size. Against REBOUND doing the
 *     same, and against the same cft run never changed, which must differ
 *     - or the change never reached the arithmetic and the comparison is a
 *     pass against itself. */
static void case_active_change(const char *label, double dt, double epsilon, size_t steps){
    struct groute before[2];
    printf("%s (6 massive bodies; N_active all, then 3, then all; %zu steps each)\n", label, light_ran(steps));
    struct reb_simulation *ra = build(grav6, 6, dt, epsilon, 0);
    struct reb_simulation *rc = NULL;
    reb_simulation_steps(ra, steps);
    ra->N_active = 3; ra->testparticle_type = 0; ra->testparticle_hidewarnings = 1;
    reb_simulation_steps(ra, steps);
    ra->N_active = SIZE_MAX;
    reb_simulation_steps(ra, steps);

    struct reb_simulation *rb = build(grav6, 6, dt, epsilon, 1);
    routes_now(before);
    reb_simulation_steps(rb, steps);
    rb->N_active = 3; rb->testparticle_type = 0; rb->testparticle_hidewarnings = 1;
    reb_simulation_steps(rb, steps);
    rb->N_active = SIZE_MAX;
    reb_simulation_steps(rb, steps);
    compare(ra, rb, label);
    routes_moved(label, before);
    /* the control: the same cft run with no change. One engine is one
     * simulation at a time, so rb is snapshotted and freed first. */
    double snap[6][6];
    const size_t n = rb->N;
    for (size_t i = 0; i < n && i < 6; i++){
        const struct reb_particle *q = &rb->particles[i];
        const double v[6] = { q->x, q->y, q->z, q->vx, q->vy, q->vz };
        memcpy(snap[i], v, sizeof v);
    }
    reb_simulation_free(rb);
    rc = build(grav6, 6, dt, epsilon, 1);
    reb_simulation_steps(rc, 3 * steps);
    int differ = 0;
    for (size_t i = 0; i < n && i < rc->N; i++){
        const struct reb_particle *q = &rc->particles[i];
        const double v[6] = { q->x, q->y, q->z, q->vx, q->vy, q->vz };
        for (int k = 0; k < 6; k++) differ += bits_differ(v[k], snap[i][k]);
    }
    if (!differ){
        printf("  FAIL %s: the run with N_active changed is bit-identical to the run without, so the\n"
               "       change never reached the arithmetic and the comparison above proves nothing\n", label);
        failures++;
    }else
        printf("  ok   %s: the change reached the arithmetic, %d of %zu values differ from the unchanged run\n",
               label, differ, 6 * n);
    reb_simulation_free(ra);
    reb_simulation_free(rc);
}

/* 3. Test particles that pull back, on the routes. */
static void case_pull_back(const char *label, double dt, double epsilon, size_t steps){
    struct groute before[2];
    printf("%s (6 massive bodies, N_active = 3, testparticle_type = 1, %zu steps)\n", label, light_ran(steps));
    struct reb_simulation *ra = build(grav6, 6, dt, epsilon, 0);
    struct reb_simulation *rb = build(grav6, 6, dt, epsilon, 1);
    for (int i = 0; i < 2; i++){
        struct reb_simulation *r = i ? rb : ra;
        r->N_active = 3; r->testparticle_type = 1; r->testparticle_hidewarnings = 1;
    }
    reb_simulation_steps(ra, steps);
    routes_now(before);
    reb_simulation_steps(rb, steps);
    compare(ra, rb, label);
    routes_moved(label, before);
    reb_simulation_free(ra);
    reb_simulation_free(rb);
}

void cases_gravity(void){
    printf("\ngravity's routes in the library half: indexed (the fold) and pair_gather (the pair inputs)\n");
    case_routes_so_far();
    printf("\nthe tables, rebuilt mid-run, against REBOUND's ias15:\n");
    case_merge("a merge mid-run", 0.01, 1e-9, 2000);
    case_active_change("N_active changed mid-run", 0.3, 1e-9, 200);
    printf("\ntest particles, on the routes:\n");
    case_pull_back("massive test particles that pull back", 0.3, 1e-9, 300);
    printf("\n");
}
