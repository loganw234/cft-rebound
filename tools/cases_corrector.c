/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The corrector's exit test through the drop-in (ROUND3.md, P2).
 *
 * The library half takes the exit test's route as the program does,
 * through the door: on a device that publishes CFT_FEAT_REDUCE_SEG - and
 * the software backend always does - pc_error()'s maxima are the
 * normalabs program and cft_reduce_seg, and every other case in this gate
 * already holds that route to REBOUND's own ias15, bit for bit. What no
 * other case can see is WHICH route ran, because the drop-in has no
 * trailer. So this reads ias15_engine_route() around integrations of its
 * own, in each of REBOUND's four modes (INDIVIDUAL takes the route's
 * other shape), and asserts that the route the engine says it took is
 * the one that carried the maxima: calls on the device route and no
 * fallback, or fallbacks and no call, and never neither. Each integration
 * is held to REBOUND's own as well, as every case here is.
 *
 * The engine opens once a process, so the route is whatever the
 * environment chose at the first bind. The case DERIVES which state that
 * must be from the environment, as the door's policy defines it
 * (src/ias15_cft.c, engine_routes()), and holds the engine to it, by= and
 * all:
 *   - CFT_REBOUND_ROUTE_REDUCE_SEG=on or =off: that state, by=forced;
 *   - otherwise CFT_REBOUND_HIDE_FEATURES naming reduce_seg: off,
 *     by=hidden;
 *   - otherwise the software backend: off, by=opt_in. It publishes every
 *     bit, and the route is opt-in (src/ias15_cft.c, engine_routes());
 *   - otherwise a card: off by=opt_in, or off by=absent.
 * `off by=device` is never right anywhere: a door that fell back without
 * saying so (V2's M4, 2026-11-02). Accepting whichever route opened would
 * pass a library half silently on today's route, and with the route
 * opt-in only a run that forces it on can tell: tools/check_corrector.py
 * --library runs this case so. It also runs both routes explicitly, and
 * derives each count, through the program's `# route` lines.
 *
 * Its steps are a count, scaled under --light: the property needs one
 * corrector pass and holds at any count.
 */
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "dropin_cases.h"
#include "ias15_engine.h"

static const char *const route = "reduce_seg";

static int reduce_seg_route(struct ias15_engine_route *out){
    for (int r = 0; ias15_engine_route(r, out) == 0; r++)
        if (!strcmp(out->name, route)) return 0;
    return -1;
}

/* Whether a comma-separated list names the route, item for item. */
static int list_names(const char *list, const char *name){
    size_t n = strlen(name);
    for (const char *p = list; p && *p; ){
        size_t l = strcspn(p, ",");
        if (l == n && !strncmp(p, name, n)) return 1;
        p += l;
        if (*p == ',') p++;
    }
    return 0;
}

/* The state the door must have opened in, from the environment: 1 on or 0
 * off, and the by= that alone can say so. A card may say either of two. */
static void expected(int *on, const char **by, const char **by_card_off){
    char var[64] = "CFT_REBOUND_ROUTE_";
    size_t l = strlen(var);
    for (const char *p = route; *p && l + 1 < sizeof var; p++) var[l++] = (char)toupper((unsigned char)*p);
    var[l] = 0;
    const char *want = getenv(var), *hide = getenv("CFT_REBOUND_HIDE_FEATURES"), *art = getenv("CFT_REBOUND_ARTIFACT");
    *by_card_off = NULL;
    if (want && (!strcmp(want, "on") || !strcmp(want, "off"))){ *on = !strcmp(want, "on"); *by = "forced"; }
    else if (hide && list_names(hide, route)){ *on = 0; *by = "hidden"; }
    else if (!art || !art[0]){ *on = 0; *by = "opt_in"; }
    else { *on = 0; *by = "opt_in"; *by_card_off = "absent"; }
}

static void case_route(const char *label, const struct body *bs, size_t n, double dt, double epsilon,
                       size_t steps, int mode){
    struct ias15_engine_route before, after;
    unsigned long long dc, df;
    int save = adaptive_mode;
    adaptive_mode = mode;
    printf("%s (N = %zu, dt = %g, epsilon = %g, %zu steps, adaptive_mode %d)\n",
           label, n, dt, epsilon, light_ran(steps), mode);
    struct reb_simulation *ra = build(bs, n, dt, epsilon, 0);
    struct reb_simulation *rb = build(bs, n, dt, epsilon, 1);
    if (reduce_seg_route(&before)){
        printf("  FAIL %s: the engine names no route reduce_seg\n", label);
        failures++;
    }else{
        /* An engine not yet opened reads zero, so the difference is this
         * case's alone either way. */
        reb_simulation_steps(ra, steps);
        reb_simulation_steps(rb, steps);
        compare(ra, rb, label);
        reduce_seg_route(&after);
        dc = after.calls - before.calls;
        df = after.fallbacks - before.fallbacks;
        int want_on;
        const char *want_by, *card_off;
        expected(&want_on, &want_by, &card_off);
        const int state = (after.on == want_on && !strcmp(after.by, want_by)) ||
                          (card_off && !after.on && !strcmp(after.by, card_off));
        const int counts = after.on ? (dc > 0 && df == 0) : (df > 0 && dc == 0);
        if (state && counts)
            printf("  ok   %s: %s=%s by=%s, as this environment requires, carried this case's %llu maxima %s\n", label,
                   route, after.on ? "on" : "off", after.by, after.on ? dc : df,
                   after.on ? "on the device route, and no fallback" : "on today's route, and no device call");
        else{
            printf("  FAIL %s: %s=%s by=%s, where this environment requires %s by=%s%s%s; the case moved its calls "
                   "by %llu and its fallbacks by %llu\n", label, route, after.on ? "on" : "off", after.by,
                   want_on ? "on" : "off", want_by, card_off ? " or off by=" : "", card_off ? card_off : "", dc, df);
            failures++;
        }
    }
    reb_simulation_free(ra);
    reb_simulation_free(rb);
    adaptive_mode = save;
}

void cases_corrector(void){
    printf("\ncorrector: the exit test's route through ias15_engine_route(), each integration against REBOUND's ias15\n");
    case_route("corrector kepler INDIVIDUAL", kepler, 2, 0.01, 1e-9, 100, 0);
    case_route("corrector kepler GLOBAL", kepler, 2, 0.01, 1e-9, 100, 1);
    case_route("corrector kepler PRS23", kepler, 2, 0.01, 1e-9, 100, 2);
    case_route("corrector kepler AARSETH85", kepler, 2, 0.01, 1e-9, 100, 3);
    case_route("corrector pythagorean PRS23", pythagorean, 3, 0.01, 1e-9, 100, 2);
}
