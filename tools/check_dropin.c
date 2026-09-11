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
 *
 * The file this comment heads is now main(), the refusals and the
 * coverage check. The cases live beside it, a file per topic:
 * tools/cases_core.c, tools/cases_wide.c, and whatever a parcel adds.
 * tools/dropin_cases.h says how to add one.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "dropin_cases.h"
#include "cft_supported.h"

static void coverage_report(void);
static void mark_covered(const char *row);

/* The message itself goes to stderr, which is where a user would meet it;
 * what is asserted here is that the integration STOPPED and that the
 * clock did not move - nothing was computed. Positions are deliberately
 * not asserted: REBOUND's own boundary check and collision search run
 * after the integrator callback returns and may touch them, and that is
 * REBOUND's doing, not this integrator's. */
int refused(const char *row, const char *what,
            void (*poison)(struct reb_simulation *)){
    mark_covered(row);
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

/* additional_forces and force_is_velocity_dependent were poisoned here
 * until they were implemented. They are now drop-in capabilities and
 * their rows are subprocess-only, so their cases moved to
 * tools/cases_forces.c - where they are accepted() and then run against
 * REBOUND's own ias15 with the same routine, which is the assertion
 * that retires a refusal rather than merely deleting it. */
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
int accepted(const char *row, const char *what,
             void (*prepare)(struct reb_simulation *)){
    mark_covered(row);
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
static void p_mode(struct reb_simulation *r){ cft_ias15_get_state(r)->adaptive_mode = 0; }


/* ------------------------------------------------------------------ *
 * Coverage over cft_support_rows
 *
 * refused() and accepted() each name the row they exercise, and this
 * walks the table afterwards. A drop-in row no case named is a FAILURE,
 * by name.
 *
 * That is the second half of what P0 bought. The first was that the
 * refusal list stopped being four copies; this is what stops the one
 * copy and its tests from drifting apart the way the four copies did.
 * It earned its place immediately: four rows - the particle map, the
 * ensemble count, the format and max_iter - had never been poisoned by
 * anything, and now are.
 */
#define MAX_COVERED 64
static const char *covered[MAX_COVERED];
static int n_covered;

static void mark_covered(const char *row){
    if (!row) return;
    for (int i = 0; i < n_covered; i++) if (!strcmp(covered[i], row)) return;
    if (n_covered < MAX_COVERED) covered[n_covered++] = row;
    else { printf("  FAIL coverage: more than %d rows exercised; raise MAX_COVERED\n",
                  MAX_COVERED); failures++; }
}

static int is_covered(const char *row){
    for (int i = 0; i < n_covered; i++) if (!strcmp(covered[i], row)) return 1;
    return 0;
}

static void coverage_report(void){
    int missing = 0, rows = 0;
    for (const struct cft_support_row *row = cft_support_rows; row->name; row++){
        if (!(row->paths & CFT_PATH_DROPIN)) continue;
        rows++;
        if (is_covered(row->name)) continue;
        printf("  FAIL coverage: cft_support_rows has \"%s\" and no case exercises it.\n"
               "       Add one to case_refusals(), or to your topic's file, naming that row.\n",
               row->name);
        failures++; missing++;
    }
    /* The other direction: a case naming a row that is gone. */
    for (int i = 0; i < n_covered; i++){
        int found = 0;
        for (const struct cft_support_row *row = cft_support_rows; row->name; row++)
            if (!strcmp(row->name, covered[i])) { found = 1; break; }
        if (!found){
            printf("  FAIL coverage: a case names row \"%s\", which is not in "
                   "cft_support_rows.\n       The refusal was removed; retire its case.\n",
                   covered[i]);
            failures++; missing++;
        }
    }
    if (!missing)
        printf("  ok   coverage: all %d drop-in rows of cft_support_rows exercised\n", rows);
}

/* Four poisons the coverage check asked for. Nothing exercised the
 * particle-map, ensemble, format or max_iter rows before the table
 * existed to be walked - which is the check earning its place on the
 * first run rather than in principle. */
static size_t one_index[1] = { 0 };
static void p_map(struct reb_simulation *r){ r->map = one_index; r->N_map = 1; }
static void p_ensemble(struct reb_simulation *r){ cft_ias15_get_state(r)->E = 2; }
static void p_format(struct reb_simulation *r){ cft_ias15_get_state(r)->format = 99; }
static void p_maxiter(struct reb_simulation *r){ cft_ias15_get_state(r)->max_iter = -1; }

static void case_refusals(void){
    printf("refusals: each unsupported feature named and the integration stopped\n");
    refused("ghost_boxes",    "ghost boxes",                  p_ghost);
    refused("gravity_module", "REB_GRAVITY_COMPENSATED",      p_gravity);
    refused("gravity_module", "the tree code",                p_tree);
    refused("gravity_custom", "a custom gravity routine",     p_gravcustom);
    refused("boundary",       "periodic boundaries",          p_boundary);
    refused("test_particles", "test particles (N_active)",    p_testp);
    refused("particle_map",   "r->map",                       p_map);
    refused("variational",    "variational particles",        p_var);
    refused("megno",          "MEGNO",                        p_megno);
    refused("adaptive_mode",  "adaptive_mode != PRS23",       p_mode);
    refused("ensemble_E",     "state->E != 1",                p_ensemble);
    refused("format",         "an out-of-range format",       p_format);
    refused("max_iter",       "a negative max_iter",          p_maxiter);

    /* And two that must NOT be refused. The first was, for one day, on
     * a premise that was not true. The second is the whole of the
     * collision row seen from binary64: the row exists for the drop-in
     * and its predicate must return 0 here - cases_wide.c proves it
     * returns 1 above binary64. */
    accepted(NULL,        "an attached ODE set",              p_odes);
    accepted("collision", "collision detection at binary64",  p_collision);

    coverage_report();
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
        cases_wide();
        if (failures){ printf("check_dropin --wide: %d FAILURES\n", failures); return 1; }
        printf("check_dropin --wide: passed\n");
        return 0;
    }

    printf("check_dropin: REBOUND's own ias15 against the registered "
           "ias15_cft at binary64, bit for bit\n\n");

    cases_core();
    /* A parcel's topic goes here, one line, beside its own file. */
    cases_forces();
    case_refusals();

    printf("\n");
    if (failures){ printf("check_dropin: %d FAILURES\n", failures); return 1; }
    printf("check_dropin: every case passed\n");
    return 0;
}
