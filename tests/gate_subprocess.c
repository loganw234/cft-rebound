/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * THE GATE FOR THE SUBPROCESS API, which did not have one.
 *
 * `libcft_rebound.a` and include/cft_rebound.h are a shipped public
 * entry point: cft_rebound_check() says whether the port can integrate
 * a simulation and cft_rebound_steps() integrates it, by writing the
 * particles out as an exact binary64 problem file, running the
 * standalone `ias15_cft` program on it, and rounding the wide result
 * back into r->particles.
 *
 * Until this file, nothing `make check` runs touched either of them.
 * The only things that did were examples/roundtrip.c and
 * examples/dropin.c, under `make example`. That is how a flag went
 * unforwarded for a month: cft_rebound_check() READ the simulation's
 * adaptive_mode to decide whether to refuse it, and the command builder
 * then dropped it on the floor, so an AARSETH85 simulation was accepted
 * and integrated at PRS23 with nothing saying so. A silent wrong
 * answer, which is the one outcome this repository is built to prevent,
 * on the one path it never tested.
 *
 * Two halves:
 *
 *   1. EQUIVALENCE, per step criterion. One call of N steps from a
 *      fresh state is exactly REBOUND's N steps from a fresh state -
 *      the header's warning about splitting a run is about the SECOND
 *      call, which starts with the b coefficients zeroed. So a single
 *      call must agree with REBOUND's own ias15 BIT FOR BIT, at every
 *      one of the four criteria.
 *
 *      And the control, without which that proves nothing: the four
 *      criteria must DISAGREE WITH EACH OTHER. Four modes all matching
 *      REBOUND is also what you get when REBOUND and the port both
 *      ignore the flag, so the assertion that matters is that the flag
 *      changes the answer at all.
 *
 *   2. REFUSALS. Every row of cft_support_rows that names this path
 *      must actually refuse, and the walk fails by name on any row no
 *      case exercised - the same coverage check tools/check_dropin.c
 *      applies to the drop-in rows, which until now covered the
 *      subprocess rows not at all.
 *
 * This program needs the ias15_cft executable. It sets
 * CFT_REBOUND_IAS15 from --program, or --build, before running.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rebound.h"
#include "cft_rebound.h"
#include "cft_supported.h"

static int failures;
static int verbose;

/* ------------------------------------------------------------------ */
/* A problem both sides are built from, so the initial conditions       */
/* cannot be the source of a difference.                                */
/* ------------------------------------------------------------------ */
struct body { double m, x, y, z, vx, vy, vz; };

/* The Kepler pair the rest of this repository's gates use: a = 1,
 * e = 1/2, barycentric, zero total momentum. */
static const struct body kepler[] = {
    { 1.0,  0.0,               0.0, 0.0,  0.0, -8.6645830828724809e-4, 0.0 },
    { 1e-3, 5.0e-1,            0.0, 0.0,  0.0,  8.6645830828724809e-1, 0.0 },
};

/* Three bodies, so a criterion that keys on the per-particle timescale
 * has more than one particle to choose between. */
static const struct body pythagorean[] = {
    { 3.0,  1.0,  3.0, 0.0,  0.0, 0.0, 0.0 },
    { 4.0, -2.0, -1.0, 0.0,  0.0, 0.0, 0.0 },
    { 5.0,  1.0, -1.0, 0.0,  0.0, 0.0, 0.0 },
};

static struct reb_simulation *build(const struct body *bs, size_t n,
                                    double dt, double epsilon, int mode){
    struct reb_simulation *r = reb_simulation_create();
    struct reb_integrator_ias15_state *s;
    r->G = 1.0;
    r->dt = dt;
    r->exact_finish_time = 1;
    s = reb_simulation_set_integrator(r, "ias15");
    if (!s){ fprintf(stderr, "gate_subprocess: ias15 missing\n"); exit(2); }
    s->epsilon = epsilon;
    s->adaptive_mode = mode;
    for (size_t i = 0; i < n; i++){
        struct reb_particle p;
        memset(&p, 0, sizeof p);
        p.m = bs[i].m;
        p.x = bs[i].x; p.y = bs[i].y; p.z = bs[i].z;
        p.vx = bs[i].vx; p.vy = bs[i].vy; p.vz = bs[i].vz;
        reb_simulation_add(r, p);
    }
    return r;
}

/* Bit patterns, never values - the same rule as every other gate. */
static int bits_differ(double a, double b){ return memcmp(&a, &b, sizeof a) != 0; }

/* How many of the 6N coordinates plus the clock differ. */
static int count_differences(struct reb_simulation *a, struct reb_simulation *b,
                             const char *label, int report){
    int bad = 0;
    if (a->N != b->N){
        if (report) printf("      N differs: %zu vs %zu\n", a->N, b->N);
        return 1000;
    }
    for (size_t i = 0; i < a->N; i++){
        const struct reb_particle *p = &a->particles[i], *q = &b->particles[i];
        const double pv[6] = { p->x, p->y, p->z, p->vx, p->vy, p->vz };
        const double qv[6] = { q->x, q->y, q->z, q->vx, q->vy, q->vz };
        static const char *nm[6] = { "x", "y", "z", "vx", "vy", "vz" };
        for (int c = 0; c < 6; c++)
            if (bits_differ(pv[c], qv[c])){
                bad++;
                if (report && (verbose || bad <= 2))
                    printf("      p%zu.%-2s  %a   vs %a\n", i, nm[c], pv[c], qv[c]);
            }
    }
    if (bits_differ(a->t, b->t)){
        bad++;
        if (report) printf("      t      %a   vs %a\n", a->t, b->t);
    }
    (void)label;
    return bad;
}

/* ------------------------------------------------------------------ */
/* 1. Equivalence, per criterion, with the control that says the flag   */
/*    travelled at all                                                  */
/* ------------------------------------------------------------------ */
static const char *mode_name(int m){
    switch (m){
        case 0: return "INDIVIDUAL";
        case 1: return "GLOBAL";
        case 2: return "PRS23";
        case 3: return "AARSETH85";
        default: return "?";
    }
}

static void case_modes(const char *label, const struct body *bs, size_t n,
                       double dt, double epsilon, long steps){
    struct reb_simulation *port[4] = { NULL, NULL, NULL, NULL };
    struct cft_rebound_options opt;
    int m;

    printf("%s (N = %zu, dt = %g, epsilon = %g, %ld steps in ONE call)\n",
           label, n, dt, epsilon, steps);

    memset(&opt, 0, sizeof opt);
    opt.format = CFT_REBOUND_FP64;
    opt.epsilon = epsilon;

    for (m = 0; m < 4; m++){
        struct reb_simulation *ref = build(bs, n, dt, epsilon, m);
        struct reb_simulation *sub = build(bs, n, dt, epsilon, m);
        int rc, bad;

        reb_simulation_steps(ref, steps);
        rc = cft_rebound_steps(sub, steps, &opt, NULL);
        if (rc != 0){
            printf("  FAIL %s at %s: cft_rebound_steps returned %d\n",
                   label, mode_name(m), rc);
            failures++;
            reb_simulation_free(ref); reb_simulation_free(sub);
            continue;
        }
        bad = count_differences(ref, sub, label, 1);
        if (bad){
            printf("  FAIL %s at %s (%d): %d values differ from REBOUND's own ias15\n",
                   label, mode_name(m), m, bad);
            failures++;
        }else{
            printf("  ok   %s at %-10s (%d): %zu particles x 6 values and the clock, "
                   "bit for bit\n", label, mode_name(m), m, sub->N);
        }
        reb_simulation_free(ref);
        port[m] = sub;                    /* kept for the control below */
    }

    /* THE CONTROL. Four criteria all matching REBOUND is also what a
     * dropped flag looks like, because REBOUND would be running its own
     * default on both sides. So the criteria have to disagree with each
     * other, and this is the assertion that a flag which never left this
     * process would fail. It is the check that was missing when
     * --adaptive-mode was not forwarded at all. */
    for (m = 0; m < 4; m++){
        int bad;
        if (m == 2 || !port[m] || !port[2]) continue;
        bad = count_differences(port[2], port[m], label, 0);
        if (bad == 0){
            printf("  FAIL %s: %s is bit-identical to PRS23, so the criterion did not\n"
                   "       reach the program - the run answered a different question\n"
                   "       from the one it was asked\n", label, mode_name(m));
            failures++;
        }else{
            printf("  ok   %s: %-10s differs from PRS23 in %d of %zu values\n",
                   label, mode_name(m), bad, 6 * port[2]->N + 1);
        }
    }
    for (m = 0; m < 4; m++) if (port[m]) reb_simulation_free(port[m]);
}

/* ------------------------------------------------------------------ */
/* 2. Refusals, with the coverage walk the drop-in gate already has     */
/* ------------------------------------------------------------------ */
#define MAX_COVERED 64
static const char *covered[MAX_COVERED];
static int n_covered;

static void mark_covered(const char *row){
    if (!row) return;
    for (int i = 0; i < n_covered; i++) if (!strcmp(covered[i], row)) return;
    if (n_covered < MAX_COVERED) covered[n_covered++] = row;
}

static void refused(const char *row, const char *what,
                    void (*poison)(struct reb_simulation *)){
    struct reb_simulation *r = build(kepler, 2, 0.05, 1e-9, 2);
    char why[512];
    int rc;
    mark_covered(row);
    poison(r);
    why[0] = 0;
    rc = cft_rebound_check(r, why, sizeof why);
    if (rc != 0 && why[0])
        printf("  ok   %s refused: %.90s%s\n", what, why, strlen(why) > 90 ? "..." : "");
    else{
        printf("  FAIL %s was NOT refused (rc = %d, why = \"%s\")\n", what, rc, why);
        failures++;
    }
    reb_simulation_free(r);
}

static void nop_forces(struct reb_simulation *r){ (void)r; }
static void nop_gravity(struct reb_simulation *r){ (void)r; }
static size_t one_index[1] = { 0 };

static void p_none(struct reb_simulation *r){ while (r->N) reb_simulation_remove_particle(r, 0); }
static void p_gravcustom(struct reb_simulation *r){ r->gravity_custom = nop_gravity; }
static void p_var(struct reb_simulation *r){ reb_simulation_add_variation_1st_order(r, -1); }
/* The FLAG alone, and no variational particles - that is the whole
 * point. reb_simulation_init_megno() adds them, so poisoning it the
 * ordinary way trips the variational row and says nothing about this
 * one. Setting the flag by hand used to be ACCEPTED here and then
 * integrated as plain gravity, reporting neither MEGNO nor a refusal. */
static void p_megno(struct reb_simulation *r){ r->calculate_megno = 1; }
static void p_forces(struct reb_simulation *r){ r->additional_forces = nop_forces; }
static void p_tsmods(struct reb_simulation *r){ r->post_timestep_modifications = nop_forces; }
static void p_veldep(struct reb_simulation *r){ r->force_is_velocity_dependent = 1; }
static void p_gravity(struct reb_simulation *r){ r->gravity = REB_GRAVITY_COMPENSATED; }
static void p_collision(struct reb_simulation *r){ r->collision = REB_COLLISION_DIRECT; }
static void p_ghost(struct reb_simulation *r){ r->N_ghost_x = 1; }
static void p_boundary(struct reb_simulation *r){ r->boundary = REB_BOUNDARY_PERIODIC; }
static void p_map(struct reb_simulation *r){ r->map = one_index; r->N_map = 1; }
static void p_testp(struct reb_simulation *r){ r->N_active = 1; }
static void p_odes(struct reb_simulation *r){ reb_ode_create(r, 1); }
static void p_integrator(struct reb_simulation *r){ reb_simulation_set_integrator(r, "leapfrog"); }
static void p_mode(struct reb_simulation *r){
    struct reb_integrator_ias15_state *s = r->integrator.state;
    if (s) s->adaptive_mode = 4;
}

/* body_count wants more particles than the cap, which is 1024 - too
 * slow to build one at a time here for no gain, so it is the one row
 * this gate names as deliberately untested rather than exercising. The
 * drop-in has no equivalent row; tools/check_bodycount.py covers the
 * program's own refusal at the cap directly.
 *
 * There is no ignore_terms case, and that is not an omission: the row
 * is gone. reb_integrator_ias15_step() writes NONE over
 * r->gravity_ignore_terms at the top of every step, so nothing a user
 * sets reaches gravity and refusing it refused a run that would have
 * been correct. This gate had such a case on its first run, written
 * against the table as it stood a merge earlier, and the coverage walk
 * below caught it - "a case names row ignore_terms, which is not in
 * cft_support_rows". That direction of the check exists for exactly
 * this and found its first victim immediately. */
static const char *const SKIPPED[] = { "body_count", NULL };

static void case_refusals(void){
    printf("\nrefusals on the subprocess path: every row that names it, by name\n");
    refused("no_particles",      "an empty simulation",            p_none);
    refused("gravity_custom",    "a custom gravity routine",       p_gravcustom);
    refused("variational",       "variational particles",          p_var);
    refused("megno",             "calculate_megno without them",   p_megno);
    refused("additional_forces", "r->additional_forces",           p_forces);
    refused("timestep_mods",     "timestep modifications",         p_tsmods);
    refused("veldep_forces",     "velocity-dependent forces",      p_veldep);
    refused("gravity_module",    "REB_GRAVITY_COMPENSATED",        p_gravity);
    refused("collision",         "collision detection",            p_collision);
    refused("ghost_boxes",       "ghost boxes",                    p_ghost);
    refused("boundary",          "periodic boundaries",            p_boundary);
    refused("particle_map",      "r->map",                         p_map);
    refused("test_particles",    "test particles (N_active)",      p_testp);
    refused("odes",              "an attached ODE set",            p_odes);
    refused("integrator_name",   "a non-IAS15 integrator",         p_integrator);
    refused("adaptive_mode",     "a mode naming no criterion",     p_mode);

    /* The coverage walk. A subprocess row no case above names is a
     * FAILURE, by name - the same rule tools/check_dropin.c applies to
     * the drop-in rows, and which nothing applied to these until this
     * file existed. */
    {
        int missing = 0, rows = 0;
        for (const struct cft_support_row *row = cft_support_rows; row->name; row++){
            int skipped = 0;
            if (!(row->paths & CFT_PATH_SUBPROCESS)) continue;
            rows++;
            for (const char *const *s = SKIPPED; *s; s++)
                if (!strcmp(*s, row->name)) skipped = 1;
            if (skipped){
                printf("  --   %s: deliberately untested here, see SKIPPED\n", row->name);
                continue;
            }
            for (int i = 0; i < n_covered; i++)
                if (!strcmp(covered[i], row->name)) goto found;
            printf("  FAIL coverage: cft_support_rows has \"%s\" on the subprocess path\n"
                   "       and no case here exercises it\n", row->name);
            failures++; missing++;
            found: ;
        }
        for (int i = 0; i < n_covered; i++){
            int exists = 0;
            for (const struct cft_support_row *row = cft_support_rows; row->name; row++)
                if (!strcmp(row->name, covered[i])){ exists = 1; break; }
            if (!exists){
                printf("  FAIL coverage: a case names row \"%s\", which is not in\n"
                       "       cft_support_rows. Retire the case.\n", covered[i]);
                failures++; missing++;
            }
        }
        if (!missing)
            printf("  ok   coverage: all %d subprocess rows of cft_support_rows accounted for\n",
                   rows);
    }
}

/* ------------------------------------------------------------------ */
int main(int argc, char **argv){
    const char *build_dir = NULL, *program = NULL;
    setvbuf(stdout, NULL, _IOLBF, 0);

    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "-v")) verbose = 1;
        else if (!strcmp(argv[i], "--build") && i + 1 < argc) build_dir = argv[++i];
        else if (!strcmp(argv[i], "--program") && i + 1 < argc) program = argv[++i];
        else { fprintf(stderr, "usage: gate_subprocess [-v] [--build DIR] "
                               "[--program PATH]\n"); return 2; }
    }

    /* The API finds the program through $CFT_REBOUND_IAS15 before it
     * tries the install path or PATH, so pointing it at the tree's own
     * build is one variable and needs no installed copy. */
    if (program || build_dir){
        static char buf[1024];
        if (program) snprintf(buf, sizeof buf, "CFT_REBOUND_IAS15=%s", program);
        else snprintf(buf, sizeof buf, "CFT_REBOUND_IAS15=%s/ias15_cft%s", build_dir,
#ifdef _WIN32
                      ".exe"
#else
                      ""
#endif
                      );
        putenv(buf);
    }

    printf("gate_subprocess: the subprocess API against REBOUND's own ias15,\n"
           "                 one call per run, at every step criterion\n\n");

    case_modes("kepler",      kepler, 2,      0.05, 1e-9, 120);
    case_modes("pythagorean", pythagorean, 3, 0.01, 1e-9, 120);
    case_refusals();

    printf("\n");
    if (failures){ printf("gate_subprocess: %d FAILURES\n", failures); return 1; }
    printf("gate_subprocess: PASS\n");
    return 0;
}
