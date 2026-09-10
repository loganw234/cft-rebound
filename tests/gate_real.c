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

/* Long enough that the wide bits reach the binary64 output. At 20 + 10
 * steps only 2 of the 14 dumped lines differ between fp64 and fp256, so
 * a restart that had truncated the wide state to binary64 could pass on
 * the other 12; at 60 + 120 the difference is broad and the comparison
 * has margin. In software this is still milliseconds; on a card it is
 * about two minutes a format, which is the reason it is not larger. */
static long steps_a = 60;    /* before the checkpoint; --steps-a */
static long steps_b = 120;   /* after it;            --steps-b */
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

/* The same values snap() collects, as exact hex floats, one per line
 * and named. Two processes' dumps are compared with diff, so the
 * comparison is over bits: %a is exact and reversible, and a
 * difference in the last bit shows as a different line rather than
 * rounding to the same 17 digits. */
static void dump(struct reb_simulation *r){
    for (size_t i = 0; i < (size_t)r->N; i++){
        const struct reb_particle *q = &r->particles[i];
        printf("p%lu m  %a\n",  (unsigned long)i, q->m);
        printf("p%lu x  %a %a %a\n", (unsigned long)i, q->x,  q->y,  q->z);
        printf("p%lu v  %a %a %a\n", (unsigned long)i, q->vx, q->vy, q->vz);
    }
    printf("t  %a\n", r->t);
    printf("dt %a\n", r->dt);
}

/* One process of a multi-process checkpoint. Each returns 0 on
 * success so a shell can chain them; the comparison is the shell's,
 * because the point is that nothing is shared between the runs. */
static int phase_straight(int format, double dt){
    struct reb_simulation *r = build(format, dt);
    reb_simulation_steps(r, steps_a + steps_b);
    dump(r);
    reb_simulation_free(r);
    return 0;
}

static int phase_save(int format, double dt, const char *path){
    struct reb_simulation *r = build(format, dt);
    reb_simulation_steps(r, steps_a);
    if (cft_archive_bind(r)){ fprintf(stderr, "gate_real: bind failed\n"); return 2; }
    remove(path);
    reb_simulation_save_to_file(r, path);
    fprintf(stderr, "gate_real: wrote %s at t=%.17g after %ld steps\n",
            path, r->t, steps_a);
    reb_simulation_free(r);
    return 0;
}

static int phase_resume(int format, const char *path){
    enum cft_archive_status st;
    struct reb_simulation *r = cft_archive_load(path, 0, format, NULL, &st);
    if (!r){
        fprintf(stderr, "gate_real: load failed: %s\n", cft_archive_status_str(st));
        return 2;
    }
    fprintf(stderr, "gate_real: loaded %s (%s) at t=%.17g\n",
            path, cft_archive_status_str(st), r->t);
    reb_simulation_steps(r, steps_b);
    dump(r);
    reb_simulation_free(r);
    return 0;
}

static int one_format(int format, const char *fname, double dt){
    int ok = 1;

    /* ---- the reference: steps_a + steps_b in one go ---- */
    struct reb_simulation *ref = build(format, dt);
    reb_simulation_steps(ref, steps_a + steps_b);
    size_t nref; unsigned char *sref = snap(ref, &nref);
    double tref = ref->t;
    reb_simulation_free(ref);          /* releases the engine for the next owner */

    /* ---- the checkpoint: steps_a, save, free ---- */
    struct reb_simulation *a = build(format, dt);
    reb_simulation_steps(a, steps_a);
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
    reb_simulation_steps(b, steps_b);
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
    int phase = 0;                  /* 0 self-contained, 1 straight, 2 save, 3 resume */
    const char *path = NULL;
    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "--fp64"))  only = 0;
        if (!strcmp(argv[i], "--fp128")) only = 1;
        if (!strcmp(argv[i], "--fp256")) only = 2;
        if (!strcmp(argv[i], "--straight")) phase = 1;
        if (!strcmp(argv[i], "--save")   && i + 1 < argc){ phase = 2; path = argv[++i]; }
        if (!strcmp(argv[i], "--resume") && i + 1 < argc){ phase = 3; path = argv[++i]; }
        if (!strcmp(argv[i], "--steps-a") && i + 1 < argc) steps_a = atol(argv[++i]);
        if (!strcmp(argv[i], "--steps-b") && i + 1 < argc) steps_b = atol(argv[++i]);
    }
    if (only < 0) only = 0;

    cft_ias15_register(CFT_IAS15_INTEGRATOR_NAME);

    /* A phase writes only its dump to stdout, so a shell can diff two
     * processes' output directly. Everything else goes to stderr. */
    if (phase == 1) return phase_straight(fmts[only], 0.01);
    if (phase == 2) return phase_save(fmts[only], 0.01, path);
    if (phase == 3) return phase_resume(fmts[only], path);

    printf("gate_real: the registered integrator through a Simulationarchive\n");
    printf("           %ld steps, checkpoint, %ld more, against %ld straight\n",
           steps_a, steps_b, steps_a + steps_b);

    /* One format per process. The engine opens once and sync_config()
     * refuses a format change, so a loop over all three would fail on
     * the second for a reason that is not about archives. The Makefile
     * runs this three times - which is also what a checkpoint is. */
    int ok = one_format(fmts[only], names[only], 0.01);
    printf("gate_real: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
