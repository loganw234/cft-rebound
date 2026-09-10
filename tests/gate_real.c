/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The gate nobody wrote: the REAL integrator through an archive.
 *
 * The four archive gates all drive tests/cft_shim_stub.c, whose own
 * header says "when parcel A lands, the gates should be re-pointed at
 * the real integrator". Parcel A landed. Nothing was re-pointed, and
 * nothing else in the suite archives anything, so the seam between the
 * two halves - a real IAS15 state written to a Simulationarchive and
 * read back - has never been executed.
 *
 * It is a seam worth a gate because the two halves hold their state
 * differently. The stub OWNS its 48 blobs. The real shim does not:
 * publish_state() in src/reb_integrator_cft.c points the state at the
 * engine's live buffers, which are file-scope statics in
 * src/ias15_cft.c. So a load that fills the state's pointers with
 * archived bytes is filling something the next step may simply
 * overwrite.
 *
 * The test is the one a checkpoint has to pass:
 *
 *     30 steps straight              ==     20 steps, save, exit,
 *                                           load, 10 steps
 *
 * bit for bit, in the particles and in r->t. Anything less than
 * identical is a checkpoint that silently changes the trajectory, which
 * for this repository is the same as a broken one.
 *
 * Two separate processes would be a stronger test still; one process is
 * enough to catch a state that is not carried, because the engine is
 * reset when a new simulation binds to it.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rebound.h"
#include "cft.h"
#include "cft_ias15.h"
#include "cft_archive.h"
#include "cft_ias15_fields.h"

#define STEPS_A   20     /* before the checkpoint */
#define STEPS_B   10     /* after it */
#define NBODY      4

static const char *ARCHIVE = "gate_real.bin";

/* A sun and three planets, spaced so the step control has something to
 * do without any close encounter. Exact binary64 values, so the two
 * runs start from identical bits. */
static void add_bodies(struct reb_simulation *r){
    struct reb_particle p;
    static const double m[NBODY]  = { 1.0,      1.0e-3,   5.0e-4,   2.0e-4 };
    static const double x[NBODY]  = { 0.0,      1.0,      2.5,      5.5    };
    static const double vy[NBODY] = { 0.0,      1.0,      0.6324555320336759, 0.4264014327112209 };
    for (int i = 0; i < NBODY; i++){
        memset(&p, 0, sizeof p);
        p.m = m[i]; p.x = x[i]; p.vy = vy[i];
        reb_simulation_add(r, p);
    }
}

static struct reb_simulation *build(int format, double dt){
    struct reb_simulation *r = reb_simulation_create();
    r->G  = 1.0;
    r->dt = dt;
    add_bodies(r);
    struct cft_ias15_state *s = reb_simulation_set_integrator(r, CFT_IAS15_INTEGRATOR_NAME);
    if (!s){ fprintf(stderr, "gate_real: set_integrator failed\n"); exit(2); }
    s->epsilon = 0.0;        /* fixed dt: the comparison is about state, not step control */
    s->format  = format;
    return r;
}

/* The particles as raw bytes, plus t and dt. This is what a user sees;
 * if these agree the checkpoint kept the trajectory. */
static unsigned char *snap(struct reb_simulation *r, size_t *len){
    size_t n = sizeof(double) * (size_t)(7 * r->N + 2);
    unsigned char *b = malloc(n);
    double *d = (double*)b;
    size_t k = 0;
    for (size_t i = 0; i < (size_t)r->N; i++){
        const struct reb_particle *p = &r->particles[i];
        d[k++] = p->m;  d[k++] = p->x;  d[k++] = p->y;  d[k++] = p->z;
        d[k++] = p->vx; d[k++] = p->vy; d[k++] = p->vz;
    }
    d[k++] = r->t; d[k++] = r->dt;
    *len = n;
    return b;
}

static int one_format(int format, const char *fname, double dt){
    int ok = 1;

    /* ---- the reference: STEPS_A + STEPS_B in one go ---- */
    struct reb_simulation *ref = build(format, dt);
    reb_simulation_steps(ref, STEPS_A + STEPS_B);
    size_t nref; unsigned char *sref = snap(ref, &nref);
    double tref = ref->t;
    reb_simulation_free(ref);          /* releases the engine for the next owner */

    /* ---- the checkpoint: STEPS_A, save, free ---- */
    struct reb_simulation *a = build(format, dt);
    reb_simulation_steps(a, STEPS_A);
    if (cft_archive_bind(a)){ printf("  %-6s FAIL  cft_archive_bind\n", fname); return 0; }
    remove(ARCHIVE);
    reb_simulation_save_to_file(a, ARCHIVE);
    double tsaved = a->t;
    reb_simulation_free(a);

    /* What actually reached the file. A gate that restores nothing and
     * compares nothing would pass; this is the guard against that. */
    struct cft_archive_info info;
    memset(&info, 0, sizeof info);
    if (cft_archive_probe(ARCHIVE, 0, &info)){
        printf("  %-6s FAIL  probe\n", fname); return 0;
    }
    printf("  %-6s archive: integrator=\"%s\" has_cft=%d cft_fields=%d blobs=%d/%d n_elem=%llu\n",
           fname, info.integrator, info.has_cft, info.n_cft_fields,
           info.n_blobs_seen, CFT_N_BLOBS, (unsigned long long)info.n_elem);
    if (!info.has_cft){
        printf("  %-6s FAIL  the real integrator wrote no cft_ fields; the wide state is not in the archive\n", fname);
        ok = 0;
    }

    /* ---- the restart ---- */
    enum cft_archive_status st;
    struct reb_simulation *b = cft_archive_load(ARCHIVE, 0, format, NULL, &st);
    if (!b){
        printf("  %-6s FAIL  load: %s\n", fname, cft_archive_status_str(st));
        return 0;
    }
    printf("  %-6s load: %s, t=%.17g (saved %.17g)\n",
           fname, cft_archive_status_str(st), b->t, tsaved);
    reb_simulation_steps(b, STEPS_B);
    size_t nb; unsigned char *sb = snap(b, &nb);

    if (nb != nref || memcmp(sb, sref, nb) != 0){
        printf("  %-6s FAIL  restart differs from the straight run\n", fname);
        printf("         t   straight %.17g   restarted %.17g\n", tref, b->t);
        const double *dr = (const double*)sref, *db = (const double*)sb;
        size_t nd = nref / sizeof(double), shown = 0;
        for (size_t i = 0; i < nd && shown < 4; i++){
            if (memcmp(&dr[i], &db[i], sizeof(double)) != 0){
                printf("         value %2lu  straight %.17g   restarted %.17g\n",
                       (unsigned long)i, dr[i], db[i]);
                shown++;
            }
        }
        ok = 0;
    } else {
        printf("  %-6s PASS  %lu values identical across the checkpoint\n",
               fname, (unsigned long)(nref / sizeof(double)));
    }

    free(sb); free(sref);
    reb_simulation_free(b);
    return ok;
}

int main(int argc, char **argv){
    int fmts[3] = { CFT_FP64, CFT_FP128, CFT_FP256 };
    const char *names[3] = { "fp64", "fp128", "fp256" };
    int only = -1;
    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "--fp64"))  only = 0;
        if (!strcmp(argv[i], "--fp128")) only = 1;
        if (!strcmp(argv[i], "--fp256")) only = 2;
    }

    cft_ias15_register(CFT_IAS15_INTEGRATOR_NAME);

    printf("gate_real: the registered integrator through a Simulationarchive\n");
    printf("           %d steps, checkpoint, %d more, against %d straight\n",
           STEPS_A, STEPS_B, STEPS_A + STEPS_B);

    /* One format per process. The engine opens once and sync_config()
     * refuses a format change, so a loop over all three would fail on
     * the second for a reason that is not about archives. The Makefile
     * runs this three times - which is also what a checkpoint is. */
    if (only < 0) only = 0;
    int ok = one_format(fmts[only], names[only], 0.01);
    printf("gate_real: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
