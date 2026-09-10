/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The drop-in, from outside the repository.
 *
 * examples/roundtrip.c uses the other way in: cft_rebound_steps(),
 * which writes a problem file and runs the ias15_cft program. This one
 * is the drop-in proper - the integrator registered inside this
 * process, selected by name, and driven by REBOUND's own
 * reb_simulation_steps(). Nothing here refers to a path inside the
 * source tree; it includes only installed headers and links only
 * installed libraries, which is the claim the README makes.
 *
 * It also takes a checkpoint, because Simulationarchive support is part
 * of the installed library rather than only of the gates, and an
 * installed library that cannot do what the gates do would not be the
 * same library.
 *
 *   make PREFIX=/usr/local
 *   ./dropin
 *
 * Two things this program does that a REBOUND program does not:
 * cft_ias15_register() once per process, and setting state->format.
 * Everything else is REBOUND's own API, which is the point.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rebound.h"
#include "cft_ias15.h"
#include "cft_archive.h"

#define ARCHIVE "dropin_checkpoint.bin"

static void add(struct reb_simulation *r, double m, double x, double vy){
    struct reb_particle p;
    memset(&p, 0, sizeof p);
    p.m = m; p.x = x; p.vy = vy;
    reb_simulation_add(r, p);
}

int main(void){
    /* Once per process, before any simulation selects it. The name is
     * what reb_simulation_set_integrator() takes, and it must be lower
     * case - REBOUND lower-cases the string from Python and compares
     * with strcmp, so a capital letter registers and can never be
     * selected. */
    cft_ias15_register("ias15_cft");

    struct reb_simulation *r = reb_simulation_create();
    r->G  = 1.0;
    r->dt = 0.01;
    add(r, 1.0,   0.0, 0.0);
    add(r, 1e-3,  1.0, 1.0);
    add(r, 5e-4,  2.5, 0.6324555320336759);

    struct cft_ias15_state *s = reb_simulation_set_integrator(r, "ias15_cft");
    if (!s){
        fprintf(stderr, "dropin: could not select ias15_cft. Was "
                        "cft_ias15_register() called?\n");
        return 1;
    }
    s->format  = CFT_FP256;   /* CFT_FP64 reproduces REBOUND bit for bit */
    s->epsilon = 0.0;         /* REBOUND's convention: 0 is a fixed step */

    /* An artifact, if one is wanted, goes here:
     *     cft_ias15_set_artifact("/path/to/cft_hw_quad.xclbin");
     * or in $CFT_REBOUND_ARTIFACT. Unset means the software backend,
     * and the two produce the same bits. */

    reb_simulation_steps(r, 40);
    printf("40 binary256 steps: t = %.17g\n", r->t);
    printf("  p1.x = %a\n", r->particles[1].x);

    /* The checkpoint. cft_archive_bind() points the simulation at the
     * descriptor list for the state's format, which is what puts the
     * wide state in the file; without it REBOUND writes a correct but
     * ordinary binary64 archive. */
    if (cft_archive_bind(r)){
        fprintf(stderr, "dropin: cft_archive_bind failed\n");
        return 1;
    }
    remove(ARCHIVE);
    reb_simulation_save_to_file(r, ARCHIVE);
    const double x_saved = r->particles[1].x;
    const double t_saved = r->t;
    reb_simulation_free(r);

    /* A separate program would start here, having registered the
     * integrator the same way. */
    enum cft_archive_status st;
    struct reb_simulation *b = cft_archive_load(ARCHIVE, 0, CFT_FP256, NULL, &st);
    if (!b){
        fprintf(stderr, "dropin: load failed: %s\n", cft_archive_status_str(st));
        return 1;
    }
    printf("reloaded (%s): t = %.17g, p1.x = %a  -> %s\n",
           cft_archive_status_str(st), b->t, b->particles[1].x,
           (b->t == t_saved && b->particles[1].x == x_saved)
               ? "the same bits it was saved with"
               : "DIFFERENT, which should not happen");

    reb_simulation_steps(b, 10);
    printf("ten more: t = %.17g\n", b->t);
    reb_simulation_free(b);
    remove(ARCHIVE);
    return 0;
}
