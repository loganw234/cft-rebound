/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The gate nobody wrote: the REAL integrator through an archive.
 *
 * The four archive gates all drive tests/cft_shim_stub.c, whose own
 * header says "when parcel A lands, the gates should be re-pointed at
 * the real integrator". Parcel A landed. Nothing was re-pointed, and
 * nothing else in the suite archived anything, so the seam between the
 * two halves - a real IAS15 state written to a Simulationarchive and
 * read back - had never been executed. This file is that gate; `make
 * check` and `make check-quick` both run it at all three formats.
 *
 * It is a seam worth a gate because the two halves hold their state
 * differently. The stub OWNS its 50 blobs. The real shim does not:
 * publish_state() in src/reb_integrator_cft.c points the state at the
 * engine's live buffers, which are file-scope statics in
 * src/ias15_cft.c. So a load that fills the state's pointers with
 * archived bytes is filling something the next step may simply
 * overwrite.
 *
 * THREE SEQUENCES now, and they test three different things.
 *
 * 1. The plain checkpoint:
 *
 *     steps_a + steps_b straight     ==     steps_a steps, save, exit,
 *                                           load, steps_b steps
 *
 *    which by default is 180 against 60 + 120 (--steps-a, --steps-b).
 *    It was 20 + 10 when written; see the note beside the defaults
 *    below for why that was not enough margin.
 *
 * 2. A particle REMOVED, then a checkpoint, then one added back below
 *    the original high-water mark. REBOUND keeps its seven coefficient
 *    levels in one flat allocation that dpcast() re-slices at the
 *    current 3N, so a shrink strands a tail above the live region and
 *    a regrow under N_allocated reads it straight back. Within one
 *    process ias15_engine_alias_resize() reproduces that; across a
 *    save, until 2026-09-11, NEITHER the tail NOR the mark was carried,
 *    and a resumed run therefore zeroed its polynomial where REBOUND
 *    aliases it. ROADMAP.md, "What the collision work leaves open,
 *    across a checkpoint", is the record. --regrow selects it.
 *
 *    The engine is one process-wide instance, so a load in THIS process
 *    inherits whatever the previous simulation left in the shadow - and
 *    what the previous simulation left is very nearly the right answer.
 *    So the in-process form of this sequence scrubs the shadow with an
 *    unrelated shrink between the save and the load, and the sequence
 *    is ALSO run across three processes by tools/check_checkpoint.py,
 *    which inherits nothing at all. Its CFT_REBOUND_NO_ALIAS control
 *    drops the mark and the shadow on the restore and must diverge.
 *
 * 3. The corrector's pass cap, and REBOUND's warning off it. REBOUND
 *    counts every predictor-corrector loop that runs out of passes in
 *    ias15->iterations_max_exceeded, archives it as REB_UINT64 and
 *    warns once when it reaches ten (integrator_ias15.c:394). The
 *    engine counts it too and nothing carried it into the state, so the
 *    warning never fired and a checkpoint lost the count. --maxiter
 *    selects it; it runs with max_iter = 1 so that the cap is reached
 *    on every step and the counter is not a matter of luck.
 *
 * Sequences 1 and 2 compare bit for bit, in the particles and in r->t.
 * Anything less than identical is a checkpoint that silently changes
 * the trajectory, which for this repository is the same as a broken
 * one.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

#include "rebound.h"
#include "cft.h"
#include "cft_ias15.h"
#include "cft_archive.h"
#include "cft_ias15_fields.h"

/* How many strings REBOUND's r->messages holds. Taken from the symbol
 * rather than written down here: it is a plain const in rebound.c and
 * this file already links librebound. */
extern const size_t reb_messages_max_N;

/* Long enough that the wide bits reach the binary64 output. At 20 + 10
 * steps only 2 of the 14 dumped lines differ between fp64 and fp256, so
 * a restart that had truncated the wide state to binary64 could pass on
 * the other 12; at 60 + 120 the difference is broad and the comparison
 * has margin. In software this is still milliseconds; on a card it is
 * about two minutes a format, which is the reason it is not larger. */
static long steps_a = 60;    /* before the checkpoint; --steps-a */
static long steps_b = 120;   /* after it;            --steps-b */
/* Sequence 2 only: steps at the SHRUNKEN count, run once on each side
 * of the checkpoint so that both sides step a three-body system. */
static long steps_m = 40;    /* --steps-m */
#define NBODY      4
/* Four bodies so the mark is 4, a removal leaves 3, and the one added
 * back lands UNDER the mark - which is the re-reading path rather than
 * the zeroing one. tools/cases_core.c's case_remove_then_add makes the
 * same choice for the same reason. */
#define REMOVE_INDEX 2

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

/* The body that comes back after the removal. Exact binary64 values
 * again, and a different orbit from the one that left, so that a
 * survivor inheriting the departed body's polynomial - which is what
 * REBOUND does - is visible in the trajectory. */
static void add_extra(struct reb_simulation *r){
    struct reb_particle p;
    memset(&p, 0, sizeof p);
    p.m = 3.0e-4; p.x = -3.25; p.vy = -0.5547001962252291; p.vz = 0.00390625;
    reb_simulation_add(r, p);
}

static struct reb_simulation *build_cfg(int format, double dt, int max_iter){
    struct reb_simulation *r = reb_simulation_create();
    r->G  = 1.0;
    r->dt = dt;
    add_bodies(r);
    struct cft_ias15_state *s = reb_simulation_set_integrator(r, CFT_IAS15_INTEGRATOR_NAME);
    if (!s){ fprintf(stderr, "gate_real: set_integrator failed\n"); exit(2); }
    s->epsilon = 0.0;        /* fixed dt: the comparison is about state, not step control */
    s->format  = format;
    s->max_iter = max_iter;  /* 0 = the default for the format */
    return r;
}

static struct reb_simulation *build(int format, double dt){
    return build_cfg(format, dt, 0);
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

/* FNV-1a over every byte of the WIDE state: all CFT_N_BLOBS blobs at
 * n_elem elements, in the walker's order, plus the two counts.
 *
 * The particles are the binary64 VIEW, and above binary64 the view is
 * not sensitive to everything the state carries: a corrector that
 * starts from different b coefficients re-converges to within 1e-34 at
 * binary128, and the binary64 rounding of that is the same bits. Two
 * runs can therefore dump identical particles from states that differ
 * in thirty digits - which is exactly what the NO_ALIAS control did at
 * fp128 and fp256 before this line existed, leaving the comparison
 * passing and proving nothing. Hashing the wide bytes is what makes
 * the dump say what the state is rather than what it rounds to. */
static uint64_t wide_digest(struct reb_simulation *r){
    struct cft_ias15_state *s = (struct cft_ias15_state*)r->integrator.state;
    uint64_t h = 1469598103934665603ULL;
    size_t w;
    if (!s) return 0;
    w = cft_format_size((cft_format)s->format);
    for (int i = 0; i < CFT_N_BLOBS; i++){
        unsigned char **b = cft_archive_state_blob(s, i);
        if (!b || !*b) continue;
        for (size_t k = 0; k < s->n_elem * w; k++){ h ^= (*b)[k]; h *= 1099511628211ULL; }
    }
    for (size_t v = s->n_elem, k = 0; k < sizeof v; k++){ h ^= (v >> (8*k)) & 0xff; h *= 1099511628211ULL; }
    for (size_t v = s->hiwater_n_elem, k = 0; k < sizeof v; k++){ h ^= (v >> (8*k)) & 0xff; h *= 1099511628211ULL; }
    return h;
}

/* The same values snap() collects, as exact hex floats, one per line
 * and named, plus the wide state's digest. Two processes' dumps are
 * compared with diff, so the comparison is over bits: %a is exact and
 * reversible, and a difference in the last bit shows as a different
 * line rather than rounding to the same 17 digits. */
static void dump(struct reb_simulation *r){
    for (size_t i = 0; i < (size_t)r->N; i++){
        const struct reb_particle *q = &r->particles[i];
        printf("p%lu m  %a\n",  (unsigned long)i, q->m);
        printf("p%lu x  %a %a %a\n", (unsigned long)i, q->x,  q->y,  q->z);
        printf("p%lu v  %a %a %a\n", (unsigned long)i, q->vx, q->vy, q->vz);
    }
    printf("t  %a\n", r->t);
    printf("dt %a\n", r->dt);
    printf("wide %016llx\n", (unsigned long long)wide_digest(r));
}

/* ------------------------------------------------------------------ */
/* Sequence 2: the two halves of the remove-checkpoint-regrow run       */
/* ------------------------------------------------------------------ */
/* Up to the checkpoint: run, remove one, run at the smaller count. */
static int shrink_phase1(struct reb_simulation *r){
    reb_simulation_steps(r, steps_a);
    if (reb_simulation_remove_particle(r, REMOVE_INDEX) != 0){   /* 0 is success */
        fprintf(stderr, "gate_real: REBOUND refused the removal\n");
        return 2;
    }
    reb_simulation_steps(r, steps_m);
    return 0;
}

/* From the checkpoint on: more steps at the smaller count, then one
 * body back UNDER the mark, then the steps that make the divergence
 * visible in the particles. */
static int shrink_phase2(struct reb_simulation *r){
    reb_simulation_steps(r, steps_m);
    add_extra(r);
    if ((size_t)r->N != (size_t)NBODY){
        fprintf(stderr, "gate_real: N is %zu after the add, so it passed the mark and "
                        "this sequence exercises the zeroing path rather than the "
                        "re-reading one\n", r->N);
        return 2;
    }
    reb_simulation_steps(r, steps_b);
    return 0;
}

/* The engine is one process-wide instance and its coefficient shadow
 * outlives the simulation that filled it, so an in-process load would
 * otherwise inherit very nearly the tail it is supposed to be restoring
 * from the file. This writes a DIFFERENT tail over it: four bodies on
 * perturbed orbits, a removal, and enough steps that nothing of the
 * previous run is left below the mark either. */
static void scrub_shadow(int format, double dt){
    struct reb_simulation *z = build(format, dt);
    for (size_t i = 0; i < (size_t)z->N; i++){
        z->particles[i].vy *= 1.125;
        z->particles[i].z  += 0.015625 * (double)(i + 1);
    }
    reb_simulation_steps(z, 25);
    if (reb_simulation_remove_particle(z, 1) == 0)
        reb_simulation_steps(z, 25);
    reb_simulation_free(z);
}

/* ------------------------------------------------------------------ */
/* One process of a multi-process checkpoint                           */
/* ------------------------------------------------------------------ */
/* Each returns 0 on success so a shell can chain them; the comparison
 * is the shell's, because the point is that nothing is shared between
 * the runs. */
static int phase_straight(int format, double dt, int regrow){
    struct reb_simulation *r = build(format, dt);
    if (regrow){
        if (shrink_phase1(r) || shrink_phase2(r)){ reb_simulation_free(r); return 2; }
    }else{
        reb_simulation_steps(r, steps_a + steps_b);
    }
    dump(r);
    reb_simulation_free(r);
    return 0;
}

static int phase_save(int format, double dt, const char *path, int regrow){
    struct reb_simulation *r = build(format, dt);
    if (regrow){
        if (shrink_phase1(r)){ reb_simulation_free(r); return 2; }
    }else{
        reb_simulation_steps(r, steps_a);
    }
    if (cft_archive_bind(r)){ fprintf(stderr, "gate_real: bind failed\n"); return 2; }
    remove(path);
    reb_simulation_save_to_file(r, path);
    fprintf(stderr, "gate_real: wrote %s at t=%.17g with %zu particles after %ld steps\n",
            path, r->t, r->N, regrow ? steps_a + steps_m : steps_a);
    reb_simulation_free(r);
    return 0;
}

static int phase_resume(int format, const char *path, int regrow){
    enum cft_archive_status st;
    struct reb_simulation *r = cft_archive_load(path, 0, format, NULL, &st);
    if (!r){
        fprintf(stderr, "gate_real: load failed: %s\n", cft_archive_status_str(st));
        return 2;
    }
    struct cft_ias15_state *s = (struct cft_ias15_state*)r->integrator.state;
    fprintf(stderr, "gate_real: loaded %s (%s, provenance \"%s\") at t=%.17g\n",
            path, cft_archive_status_str(st),
            cft_ias15_provenance_str(s ? s->provenance : -1), r->t);
    if (regrow){
        if (shrink_phase2(r)){ reb_simulation_free(r); return 2; }
    }else{
        reb_simulation_steps(r, steps_b);
    }
    dump(r);
    reb_simulation_free(r);
    return 0;
}

/* ------------------------------------------------------------------ */
/* Sequence 1, in one process                                          */
/* ------------------------------------------------------------------ */
static int one_format(int format, const char *fname, double dt){
    int ok = 1;

    /* ---- the reference: steps_a + steps_b in one go ---- */
    struct reb_simulation *ref = build(format, dt);
    reb_simulation_steps(ref, steps_a + steps_b);
    size_t nref; unsigned char *sref = snap(ref, &nref);
    double tref = ref->t;
    uint64_t dref = wide_digest(ref);
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
    printf("  %-6s archive: integrator=\"%s\" has_cft=%d cft_fields=%d blobs=%d/%d "
           "alias=%d/%d n_elem=%llu hiwater_n_elem=%llu provenance=\"%s\"\n",
           fname, info.integrator, info.has_cft, info.n_cft_fields,
           info.n_blobs_seen, CFT_N_BLOBS,
           info.n_alias_blobs_seen, CFT_N_ALIAS_BLOBS,
           (unsigned long long)info.n_elem, (unsigned long long)info.hiwater_n_elem,
           cft_ias15_provenance_str(info.provenance));
    if (!info.has_cft){
        printf("  %-6s FAIL  the real integrator wrote no cft_ fields; the wide state is not in the archive\n", fname);
        ok = 0;
    }
    /* Nothing was ever removed here, so the mark is the live count and
     * the alias family must be absent - which is also what keeps this
     * archive byte for byte what the project wrote before that family
     * existed. */
    if (info.hiwater_n_elem || info.n_alias_blobs_seen){
        printf("  %-6s FAIL  a run that never shrank wrote %d high-water blobs and "
               "cft_hiwater_n_elem=%llu; it should carry neither\n",
               fname, info.n_alias_blobs_seen, (unsigned long long)info.hiwater_n_elem);
        ok = 0;
    }

    /* ---- the restart ---- */
    enum cft_archive_status st;
    struct reb_simulation *b = cft_archive_load(ARCHIVE, 0, format, NULL, &st);
    if (!b){
        printf("  %-6s FAIL  load: %s\n", fname, cft_archive_status_str(st));
        return 0;
    }
    struct cft_ias15_state *bs = (struct cft_ias15_state*)b->integrator.state;
    printf("  %-6s load: %s, t=%.17g (saved %.17g), state says \"%s\"\n",
           fname, cft_archive_status_str(st), b->t, tsaved,
           cft_ias15_provenance_str(bs ? bs->provenance : -1));
    /* The return value said "restored exactly". The STATE has to say so
     * too, because a caller may discard a return value and nothing else
     * in the state recorded it - ROADMAP.md line 404. */
    if (!bs || bs->provenance != CFT_PROV_EXACT){
        printf("  %-6s FAIL  the load returned %s and the state's provenance is \"%s\"\n",
               fname, cft_archive_status_str(st),
               cft_ias15_provenance_str(bs ? bs->provenance : -1));
        ok = 0;
    }
    reb_simulation_steps(b, steps_b);
    size_t nb; unsigned char *sb = snap(b, &nb);
    uint64_t db = wide_digest(b);

    /* The particles are the binary64 view; above binary64 they are not
     * sensitive to everything the state carries. See wide_digest(). */
    if (db != dref){
        printf("  %-6s FAIL  the WIDE state differs across the checkpoint: "
               "%016llx straight, %016llx restarted\n",
               fname, (unsigned long long)dref, (unsigned long long)db);
        ok = 0;
    }
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
    } else if (ok) {
        printf("  %-6s PASS  %lu values and the wide state (%016llx) identical across "
               "the checkpoint\n",
               fname, (unsigned long)(nref / sizeof(double)), (unsigned long long)dref);
    }

    free(sb); free(sref);
    reb_simulation_free(b);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Sequence 2, in one process                                          */
/* ------------------------------------------------------------------ */
static int one_format_regrow(int format, const char *fname, double dt){
    int ok = 1;

    /* ---- the reference: the whole sequence in one go ---- */
    struct reb_simulation *ref = build(format, dt);
    if (shrink_phase1(ref) || shrink_phase2(ref)){
        printf("  %-6s FAIL  the straight sequence would not run\n", fname);
        reb_simulation_free(ref); return 0;
    }
    size_t nref; unsigned char *sref = snap(ref, &nref);
    double tref = ref->t;
    uint64_t dref = wide_digest(ref);
    reb_simulation_free(ref);

    /* ---- run to the checkpoint and save ---- */
    struct reb_simulation *a = build(format, dt);
    if (shrink_phase1(a)){
        printf("  %-6s FAIL  the checkpoint run would not reach the save\n", fname);
        free(sref); reb_simulation_free(a); return 0;
    }
    if (cft_archive_bind(a)){ printf("  %-6s FAIL  cft_archive_bind\n", fname);
                              free(sref); reb_simulation_free(a); return 0; }
    remove(ARCHIVE);
    reb_simulation_save_to_file(a, ARCHIVE);
    double tsaved = a->t;
    size_t nsaved = a->N;
    reb_simulation_free(a);

    struct cft_archive_info info;
    memset(&info, 0, sizeof info);
    if (cft_archive_probe(ARCHIVE, 0, &info)){
        printf("  %-6s FAIL  probe\n", fname); free(sref); return 0;
    }
    unsigned long long mark = cft_archive_hiwater_N(&info);
    printf("  %-6s archive: %zu live particles (n_elem %llu), mark %llu, alias=%d/%d "
           "blobs of %llu elements (cft_hiwater_n_elem=%llu)\n",
           fname, nsaved, (unsigned long long)info.n_elem, mark,
           info.n_alias_blobs_seen, CFT_N_ALIAS_BLOBS,
           (unsigned long long)(info.alias_blob_bytes /
               (info.format > 0 ? cft_format_size((cft_format)info.format) : 1)),
           (unsigned long long)info.hiwater_n_elem);

    /* This sequence is only a test of anything if the file actually
     * records a mark ABOVE the live count. If it does not, the resumed
     * run would take r->N for its mark and be right by accident. */
    if (mark != (unsigned long long)NBODY || info.n_alias_blobs_seen != CFT_N_ALIAS_BLOBS){
        printf("  %-6s FAIL  the archive records a high-water mark of %llu and %d of %d "
               "coefficient blobs; without both this sequence proves nothing\n",
               fname, mark, info.n_alias_blobs_seen, CFT_N_ALIAS_BLOBS);
        ok = 0;
    }

    /* Make sure the load has to get the tail from the FILE. */
    scrub_shadow(format, dt);

    /* ---- resume, step, and grow back under the mark ---- */
    enum cft_archive_status st;
    struct reb_simulation *b = cft_archive_load(ARCHIVE, 0, format, NULL, &st);
    if (!b){
        printf("  %-6s FAIL  load: %s\n", fname, cft_archive_status_str(st));
        free(sref); return 0;
    }
    printf("  %-6s load: %s, t=%.17g (saved %.17g), %zu particles\n",
           fname, cft_archive_status_str(st), b->t, tsaved, b->N);
    if (shrink_phase2(b)){
        printf("  %-6s FAIL  the resumed sequence would not run\n", fname);
        free(sref); reb_simulation_free(b); return 0;
    }
    size_t nb; unsigned char *sb = snap(b, &nb);
    uint64_t db = wide_digest(b);

    if (db != dref){
        printf("  %-6s FAIL  the WIDE state differs across the regrown checkpoint: "
               "%016llx straight, %016llx restarted\n",
               fname, (unsigned long long)dref, (unsigned long long)db);
        ok = 0;
    }
    if (nb != nref || memcmp(sb, sref, nb) != 0){
        printf("  %-6s FAIL  the regrown restart differs from the straight run\n", fname);
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
    } else if (ok) {
        printf("  %-6s PASS  %lu values and the wide state (%016llx) identical across a "
               "removal, a checkpoint and a regrow under the mark\n",
               fname, (unsigned long)(nref / sizeof(double)), (unsigned long long)dref);
    }

    /* And the corner an APPENDED snapshot opens, which is what a real
     * Simulationarchive is: snapshot 0 here holds three particles with
     * a mark of four, and the snapshot appended now holds four of each.
     * Two things could go wrong and both were live possibilities.
     *
     * The high-water family must still be there. REBOUND stores an
     * append as a DIFF, and binarydata.c writes a header and a name but
     * no DATA for a field that the new snapshot no longer has, so a
     * family that appeared and then vanished would corrupt the append -
     * which is why publish_state keeps publishing it once the mark has
     * been above the live count, even after a regrow levels the two.
     *
     * And the snapshot must LOAD. Its fifty live blobs are twelve
     * elements against snapshot 0's nine, and reading that as a file
     * disagreeing with itself refused every archive appended after a
     * particle was removed. */
    {
        struct cft_archive_info app;
        reb_simulation_save_to_file(b, ARCHIVE);
        memset(&app, 0, sizeof app);
        if (cft_archive_probe(ARCHIVE, -1, &app)){
            printf("  %-6s FAIL  probe of the appended snapshot\n", fname);
            ok = 0;
        }else{
            printf("  %-6s appended snapshot: %lld snapshots, n_elem %llu (base 9), "
                   "has_hiwater=%d cft_hiwater_n_elem=%llu, alias=%d/%d\n",
                   fname, (long long)app.n_snapshots, (unsigned long long)app.n_elem,
                   app.has_hiwater, (unsigned long long)app.hiwater_n_elem,
                   app.n_alias_blobs_seen, CFT_N_ALIAS_BLOBS);
            if (!app.has_hiwater || app.n_alias_blobs_seen != CFT_N_ALIAS_BLOBS){
                printf("  %-6s FAIL  the high-water family left the archive when the "
                       "regrow levelled the mark (%d of %d blobs, cft_hiwater_n_elem "
                       "present=%d); a diff cannot express a field that has gone\n",
                       fname, app.n_alias_blobs_seen, CFT_N_ALIAS_BLOBS, app.has_hiwater);
                ok = 0;
            }
            enum cft_archive_status sa;
            struct reb_simulation *d = cft_archive_load(ARCHIVE, -1, format, NULL, &sa);
            if (!d || sa != CFT_ARCHIVE_EXACT){
                printf("  %-6s FAIL  the appended snapshot would not load: %s\n",
                       fname, cft_archive_status_str(sa));
                ok = 0;
            }else{
                struct cft_ias15_state *ds = (struct cft_ias15_state*)d->integrator.state;
                unsigned long long dm = cft_archive_hiwater_N(&app);
                printf("  %-6s the appended snapshot loads: %s, %zu particles, mark %llu\n",
                       fname, cft_archive_status_str(sa), d->N, dm);
                if (!ds || ds->n_elem != 3 * d->N || dm != (unsigned long long)NBODY){
                    printf("  %-6s FAIL  it came back with n_elem %zu for %zu particles "
                           "and a mark of %llu\n", fname, ds ? ds->n_elem : 0, d->N, dm);
                    ok = 0;
                }
            }
            if (d) reb_simulation_free(d);
        }
    }

    free(sb); free(sref);
    reb_simulation_free(b);
    return ok;
}

/* ------------------------------------------------------------------ */
/* Sequence 3: iterations_max_exceeded and REBOUND's warning           */
/* ------------------------------------------------------------------ */
/* How many stored messages carry REBOUND's non-convergence sentence.
 * r->messages[i][0] is a type byte, so the text starts at + 1. */
static int count_warnings(struct reb_simulation *r){
    int n = 0;
    if (!r->messages) return 0;
    for (size_t i = 0; i < reb_messages_max_N; i++){
        const char *m = r->messages[i];
        if (!m) break;
        if (strstr(m + 1, "predictor corrector loops in IAS15 did not converge")) n++;
    }
    return n;
}

static int one_format_maxiter(int format, const char *fname, double dt){
    int ok = 1;

    /* max_iter = 1 reaches the cap on every step, so the counter is
     * arithmetic rather than luck and the warning's threshold of ten is
     * crossed inside steps_a. The trajectory is a truncated corrector's
     * and this sequence does not compare it to anything. */
    struct reb_simulation *a = build_cfg(format, dt, 1);
    a->save_messages = 1;              /* collect rather than print */
    reb_simulation_steps(a, steps_a);
    struct cft_ias15_state *as = (struct cft_ias15_state*)a->integrator.state;
    unsigned long long counted = as ? as->iterations_max_exceeded : 0;
    int warned = count_warnings(a);
    printf("  %-6s max_iter=1: %ld steps counted %llu non-convergences, %d warning(s)\n",
           fname, steps_a, counted, warned);
    if (counted < 10){
        printf("  %-6s FAIL  the cap was reached %llu times in %ld steps; this sequence "
               "needs it above ten or the warning's threshold is never crossed\n",
               fname, counted, steps_a);
        ok = 0;
    }
    if (warned != 1){
        printf("  %-6s FAIL  REBOUND warns exactly once, on the increment that reaches "
               "ten; this run warned %d times\n", fname, warned);
        ok = 0;
    }

    if (cft_archive_bind(a)){ printf("  %-6s FAIL  cft_archive_bind\n", fname);
                              reb_simulation_free(a); return 0; }
    remove(ARCHIVE);
    reb_simulation_save_to_file(a, ARCHIVE);
    reb_simulation_free(a);

    struct cft_archive_info info;
    memset(&info, 0, sizeof info);
    if (cft_archive_probe(ARCHIVE, 0, &info)){
        printf("  %-6s FAIL  probe\n", fname); return 0;
    }
    if (info.iterations_max_exceeded != counted){
        printf("  %-6s FAIL  the state counted %llu and the archive says %llu\n",
               fname, counted, (unsigned long long)info.iterations_max_exceeded);
        ok = 0;
    }

    enum cft_archive_status st;
    struct reb_simulation *b = cft_archive_load(ARCHIVE, 0, format, NULL, &st);
    if (!b){
        printf("  %-6s FAIL  load: %s\n", fname, cft_archive_status_str(st));
        return 0;
    }
    b->save_messages = 1;
    struct cft_ias15_state *bs = (struct cft_ias15_state*)b->integrator.state;
    if (!bs || bs->iterations_max_exceeded != counted){
        printf("  %-6s FAIL  the count did not survive the checkpoint: %llu, wanted %llu\n",
               fname, bs ? bs->iterations_max_exceeded : 0, counted);
        ok = 0;
    }
    reb_simulation_steps(b, 10);
    unsigned long long after = bs ? bs->iterations_max_exceeded : 0;
    int warned_again = count_warnings(b);
    printf("  %-6s after a checkpoint and 10 more steps: %llu (was %llu), %d new warning(s)\n",
           fname, after, counted, warned_again);
    if (after <= counted){
        printf("  %-6s FAIL  the counter did not advance after the restart; it is the "
               "SIMULATION's and monotonic, not the process's\n", fname);
        ok = 0;
    }
    /* REBOUND fires on the increment that EQUALS ten, so a run resumed
     * above ten never fires again. */
    if (warned_again != 0){
        printf("  %-6s FAIL  the warning fired %d more time(s) after a restart already "
               "above the threshold\n", fname, warned_again);
        ok = 0;
    }
    if (ok) printf("  %-6s PASS  the count is the simulation's, survives a checkpoint, "
                   "and the warning fires once\n", fname);
    reb_simulation_free(b);
    return ok;
}

int main(int argc, char **argv){
    int fmts[3] = { CFT_FP64, CFT_FP128, CFT_FP256 };
    const char *names[3] = { "fp64", "fp128", "fp256" };
    int only = -1;
    int phase = 0;                  /* 0 self-contained, 1 straight, 2 save, 3 resume */
    int regrow = 0;                 /* --regrow: sequence 2 rather than sequence 1 */
    int maxiter = 0;                /* --maxiter: sequence 3 only */
    const char *path = NULL;
    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "--fp64"))  only = 0;
        if (!strcmp(argv[i], "--fp128")) only = 1;
        if (!strcmp(argv[i], "--fp256")) only = 2;
        if (!strcmp(argv[i], "--straight")) phase = 1;
        if (!strcmp(argv[i], "--regrow"))   regrow = 1;
        if (!strcmp(argv[i], "--maxiter"))  maxiter = 1;
        if (!strcmp(argv[i], "--save")   && i + 1 < argc){ phase = 2; path = argv[++i]; }
        if (!strcmp(argv[i], "--resume") && i + 1 < argc){ phase = 3; path = argv[++i]; }
        if (!strcmp(argv[i], "--steps-a") && i + 1 < argc) steps_a = atol(argv[++i]);
        if (!strcmp(argv[i], "--steps-b") && i + 1 < argc) steps_b = atol(argv[++i]);
        if (!strcmp(argv[i], "--steps-m") && i + 1 < argc) steps_m = atol(argv[++i]);
    }
    if (only < 0) only = 0;

    cft_ias15_register(CFT_IAS15_INTEGRATOR_NAME);
    /* Sequence 2 resumes a THREE-body archive and then grows back to
     * four, and several of the engine's scratch vectors - the
     * coefficient shadow among them - are sized once at their first
     * use. A resumed process that has not reserved the mark cannot hold
     * it, and is refused by name rather than silently mis-sized. */
    cft_ias15_reserve(NBODY);

    /* A phase writes only its dump to stdout, so a shell can diff two
     * processes' output directly. Everything else goes to stderr. */
    if (phase == 1) return phase_straight(fmts[only], 0.01, regrow);
    if (phase == 2) return phase_save(fmts[only], 0.01, path, regrow);
    if (phase == 3) return phase_resume(fmts[only], path, regrow);

    printf("gate_real: the registered integrator through a Simulationarchive\n");
    printf("           1. %ld steps, checkpoint, %ld more, against %ld straight\n",
           steps_a, steps_b, steps_a + steps_b);
    printf("           2. %ld steps, a particle removed, %ld, checkpoint, %ld, one added "
           "back under the mark, %ld\n", steps_a, steps_m, steps_m, steps_b);
    printf("           3. max_iter = 1: the corrector's count and REBOUND's warning, "
           "across the same checkpoint\n");

    /* One format per process. The engine opens once and sync_config()
     * refuses a format change, so a loop over all three would fail on
     * the second for a reason that is not about archives. The Makefile
     * runs this three times - which is also what a checkpoint is. */
    int ok = 1;
    if (!maxiter){
        ok = one_format(fmts[only], names[only], 0.01) && ok;
        ok = one_format_regrow(fmts[only], names[only], 0.01) && ok;
    }
    ok = one_format_maxiter(fmts[only], names[only], 0.01) && ok;
    printf("gate_real: %s\n", ok ? "PASS" : "FAIL");
    return ok ? 0 : 1;
}
