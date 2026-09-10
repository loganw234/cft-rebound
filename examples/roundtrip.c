/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * cft-rebound, worked round trip: an ordinary REBOUND simulation, one
 * call at binary128, the answer back in REBOUND's own particles.
 *
 * Build it against an installed cft-rebound (see examples/Makefile):
 *
 *     make -C examples PREFIX=/where/you/installed
 *     ./examples/roundtrip
 *
 * or from this repository's root, `make examples`.
 *
 *
 * WHICH FORM THIS USES, AND WHY
 * -----------------------------
 * Everything between the two banner comments below is the port. Today
 * that is ONE CALL, cft_rebound_steps(), and behind it a subprocess:
 * the particles go out as an exact binary64 problem file, the
 * `ias15_cft` program integrates them at binary128, and the result is
 * correctly rounded back into r->particles. That is the only form that
 * exists as this is written.
 *
 * Note that the run below is ONE call for all 300 steps, deliberately.
 * The wide state does not persist between calls and IAS15 starts each
 * step from the previous step's b coefficients, so two calls of 150
 * steps do not give the same bits as one call of 300. Ask for all the
 * steps you want at once.
 *
 * ROADMAP parcel A is building the other one. REBOUND takes
 * user-provided integrators through reb_integrator_register() and
 * reb_simulation_set_integrator(), so the port becomes a registered
 * integrator, the wide state lives across steps instead of being
 * rebuilt per call, and the middle of cft_rebound_steps() is replaced
 * by three lines of REBOUND API.
 *
 * This example is written so that swap changes nothing here. It never
 * names the program, the problem file or the record; it sets options on
 * a struct and calls one function whose signature is the analogue of
 * reb_simulation_steps(). When the registration lands, this file still
 * compiles, still runs, and still prints the same numbers.
 *
 *
 * WHAT IT PRINTS
 * --------------
 * The same system, the same fixed steps, at binary64 and at binary128.
 * Fixed steps on purpose: two formats stepping identically can be
 * differenced, and the difference IS binary64's round-off rather than
 * an estimate of it. The energy drift of each run is the headline -
 * binary64's is round-off, binary128's is the method's - and the
 * position difference at the end is what the extra format bought.
 */

#include "cft_rebound.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define NSTEPS 300
#define DT     0.5

/* ------------------------------------------------------------------ */
/* An ordinary REBOUND simulation. Nothing here is cft-specific.       */
/* ------------------------------------------------------------------ */
/* Sun, Jupiter and Saturn, from orbital elements, in G = 1 units:
 * masses in solar masses, distances in AU, so a time unit is a year
 * over 2*pi and Jupiter's period is about 74.5. */
static struct reb_simulation *build_system(void){
    struct reb_simulation *r = reb_simulation_create();
    struct reb_particle sun = {0};
    r->G  = 1.0;
    r->dt = DT;
    sun.m = 1.0;
    reb_simulation_add(r, sun);
    /*                                   G      primary      m         a     e      inc  Omega omega  f */
    reb_simulation_add(r, reb_particle_from_orbit(r->G, r->particles[0], 9.55e-4, 5.20, 0.0489, 0.0, 0.0, 0.0, 0.0));
    reb_simulation_add(r, reb_particle_from_orbit(r->G, r->particles[0], 2.86e-4, 9.55, 0.0565, 0.0, 0.0, 0.0, 1.0));
    reb_simulation_move_to_com(r);
    return r;
}

/* ------------------------------------------------------------------ */

struct outcome {
    double e0, e1;        /* the energy before and after, at the run's own format */
    char   h0[96], h1[96];/* and the same two, unrounded, as exact hex text */
    double x[3 * 8];      /* the final positions, as REBOUND holds them */
    size_t n;
    double seconds;
    long   steps;
};

static int run_at(enum cft_rebound_format fmt, struct outcome *o){
    struct reb_simulation *r = build_system();
    struct cft_rebound_options opt;
    struct cft_rebound_result res;
    size_t i;
    int rc;

    /* ============ the port begins ============================== */
    memset(&opt, 0, sizeof opt);
    opt.format  = fmt;
    opt.epsilon = 0.0;      /* REBOUND's convention: 0 is a fixed step of
                             * r->dt. Set it to 1e-9 for IAS15's adaptive
                             * stepping, exactly as you would set
                             * ias15->epsilon on a stock simulation. */

    /* the energy this format sees before anything moves */
    rc = cft_rebound_steps(r, 0, &opt, &res);
    if (rc == 0){
        o->e0 = res.energy;
        strcpy(o->h0, res.energy_hex);
        /* ... and NSTEPS steps of it. The analogue of
         *     reb_simulation_steps(r, NSTEPS);
         * and it leaves r in the state that call would: r->particles,
         * r->t, r->dt and r->steps_done updated, in binary64. */
        rc = cft_rebound_steps(r, NSTEPS, &opt, &res);
    }
    /* ============ the port ends ================================ */

    if (rc == 0){
        o->e1      = res.energy;
        o->seconds = res.seconds;
        strcpy(o->h1, res.energy_hex);
        o->steps   = res.steps_done;
        o->n       = r->N;
        for (i = 0; i < r->N && i < 8; i++){
            o->x[3 * i]     = r->particles[i].x;
            o->x[3 * i + 1] = r->particles[i].y;
            o->x[3 * i + 2] = r->particles[i].z;
        }
        /* r is an ordinary REBOUND simulation again: every REBOUND tool
         * works on it. This one recomputes the energy at binary64 from
         * the rounded view, which is NOT res.energy - the gap between
         * them is the view's rounding, not the run's error. */
        printf("  (reb_simulation_energy at binary64 afterwards: %.17g)\n",
               reb_simulation_energy(r));
    }
    reb_simulation_free(r);
    return rc;
}

/* ------------------------------------------------------------------ */
/* And what a refusal looks like.                                      */
/* ------------------------------------------------------------------ */
/* The force model is basic pairwise gravity and nothing else, so a
 * simulation asking for more is refused before it runs rather than
 * quietly integrated as something it is not. This is the README's
 * "Scope" table demonstrated: each case sets one thing and prints the
 * message it gets back. */
static void a_force(struct reb_simulation *const r){ (void)r; }

static void show_refusals(void){
    static const char *const what[] = {
        "r->softening = 0.01",
        "r->collision = DIRECT",
        "r->gravity = TREE",
        "r->additional_forces set",
        "1025 particles"
    };
    size_t k;
    for (k = 0; k < sizeof what / sizeof what[0]; k++){
        struct reb_simulation *r = build_system();
        char why[256];
        switch (k){
            case 0: r->softening = 0.01; break;
            case 1: r->collision = REB_COLLISION_DIRECT; break;
            case 2: r->gravity = REB_GRAVITY_TREE; break;
            case 3: r->additional_forces = a_force; break;
            default: {
                struct reb_particle p = {0};
                p.m = 1e-12;
                while (r->N < 1025){ p.x = (double)r->N; reb_simulation_add(r, p); }
                break;
            }
        }
        if (cft_rebound_check(r, why, sizeof why))
            printf("  %-26s -> %s\n", what[k], why);
        else
            printf("  %-26s -> NOT REFUSED, which is a bug\n", what[k]);
        reb_simulation_free(r);
    }
}

int main(void){
    struct outcome a, b;
    double d2 = 0.0;
    size_t i;

    printf("cft-rebound worked round trip: Sun, Jupiter, Saturn;"
           " %d fixed steps of dt = %g\n", NSTEPS, DT);
    printf("(the program behind this call: %s)\n\n", cft_rebound_program_path(NULL));

    printf("binary64:\n");
    if (run_at(CFT_REBOUND_FP64, &a)) return 1;
    printf("binary128:\n");
    if (run_at(CFT_REBOUND_FP128, &b)) return 1;

    printf("\n  format     steps        |dE/E|        seconds\n");
    printf("  fp64      %6ld    %12.3e    %7.2f\n", a.steps, fabs((a.e1 - a.e0) / a.e0), a.seconds);
    printf("  fp128     %6ld    %12.3e    %7.2f\n", b.steps, fabs((b.e1 - b.e0) / b.e0), b.seconds);
    printf("\n  binary128's |dE/E| is 0 because res.energy is a DOUBLE - the\n");
    printf("  binary64 view of a binary128 number - and the drift is below its\n");
    printf("  last bit. res.energy_hex is the same number unrounded:\n");
    printf("    fp128 start %s\n", b.h0);
    printf("    fp128 end   %s\n", b.h1);
    printf("    fp64  start %s\n", a.h0);
    printf("    fp64  end   %s\n", a.h1);

    for (i = 0; i < 3 * a.n && i < 24; i++){
        double d = a.x[i] - b.x[i];
        d2 += d * d;
    }
    printf("\n  the two runs took the same steps and ended %.3e AU apart.\n", sqrt(d2));
    printf("  That distance is binary64's round-off over %d steps, measured\n", NSTEPS);
    printf("  rather than estimated, and the binary128 run is the one that is right.\n");
    printf("\n  Jupiter is now at x = %.17g AU (binary128 run, binary64 view).\n", b.x[3]);

    printf("\n  and the limits, so you meet them here rather than mid-run:\n");
    show_refusals();
    return 0;
}
