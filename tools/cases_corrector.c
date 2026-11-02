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
 * environment chose at the first bind: on in the suite, which runs on the
 * software backend; off under CFT_REBOUND_ROUTE_REDUCE_SEG=off, or on a
 * card without the bit. The case asserts whichever it finds, and says
 * which. tools/check_corrector.py runs both routes explicitly, and
 * derives each count, through the program's `# route` lines.
 *
 * Its steps are a count, scaled under --light: the property needs one
 * corrector pass and holds at any count.
 */
#include <stdio.h>
#include <string.h>

#include "dropin_cases.h"
#include "ias15_engine.h"

static int reduce_seg_route(struct ias15_engine_route *out){
    for (int r = 0; ias15_engine_route(r, out) == 0; r++)
        if (!strcmp(out->name, "reduce_seg")) return 0;
    return -1;
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
        if (after.on ? (dc > 0 && df == 0) : (df > 0 && dc == 0))
            printf("  ok   %s: reduce_seg=%s by=%s carried this case's %llu maxima %s\n", label,
                   after.on ? "on" : "off", after.by, after.on ? dc : df,
                   after.on ? "on the device route, and no fallback" : "on today's route, and no device call");
        else{
            printf("  FAIL %s: reduce_seg=%s by=%s, and this case moved its calls by %llu and its fallbacks by %llu\n",
                   label, after.on ? "on" : "off", after.by, dc, df);
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
