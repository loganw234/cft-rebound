/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Gate 3: open a STOCK REBOUND archive here and confirm it is promoted
 * and says so; and the loud refusal - a binary256 archive loaded into a
 * binary64 run.
 *
 * Usage: gate_promote [dir]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cft_shim_stub.h"
#include "integrator_ias15.h"

static const char *dir = ".";
static char *path(const char *stem){
    static char buf[1024];
    snprintf(buf, sizeof buf, "%s/%s", dir, stem);
    return buf;
}

static void add_bodies(struct reb_simulation *r){
    struct reb_particle p = {0};
    p.m = 1.0;    p.x = 0;    p.y = 0;   p.z = 0;    p.vx = 0;    p.vy = 0;    p.vz = 0;
    reb_simulation_add(r, p);
    p.m = 1e-3;   p.x = 1.0;  p.y = 0;   p.z = 0;    p.vx = 0;    p.vy = 1.0;  p.vz = 0;
    reb_simulation_add(r, p);
    p.m = 1e-4;   p.x = 0;    p.y = 2.0; p.z = 0.05; p.vx = -0.7; p.vy = 0;    p.vz = 0.01;
    reb_simulation_add(r, p);
}

int main(int argc, char **argv){
    if (argc > 1) dir = argv[1];
    int fail = 0;
    cft_shim_register();

    /* ---- a stock archive: REBOUND's own IAS15, its own arithmetic --- */
    remove(path("gate3_stock.bin"));
    struct reb_simulation *s = reb_simulation_create();     /* integrator "ias15" */
    s->G = 1.0; s->dt = 0.01; s->exact_finish_time = 0;
    add_bodies(s);
    reb_simulation_steps(s, 12);
    reb_simulation_save_to_file(s, path("gate3_stock.bin"));
    reb_simulation_free(s);

    struct cft_archive_info info;
    if (cft_archive_probe(path("gate3_stock.bin"), 0, &info)){
        printf("gate 3: FAIL  could not probe the stock archive\n");
        return 1;
    }
    printf("gate 3: stock archive integrator \"%s\", has_cft = %d, cft_ fields = %d\n",
           info.integrator, info.has_cft, info.n_cft_fields);
    if (info.has_cft){ printf("gate 3: FAIL  a stock archive should carry no cft_ fields\n"); fail = 1; }

    /* the same archive read plainly, to have REBOUND's own numbers */
    struct reb_simulation *ref = reb_simulation_create_from_file(path("gate3_stock.bin"), 0);
    if (!ref || !ref->integrator.name || strcmp(ref->integrator.name, "ias15")){
        printf("gate 3: FAIL  the reference load did not give REBOUND's ias15\n");
        return 1;
    }
    struct reb_integrator_ias15_state *ia = (struct reb_integrator_ias15_state*)ref->integrator.state;

    /* and read through the archive module, which must promote */
    enum cft_archive_status st;
    struct reb_simulation *r = cft_archive_load(path("gate3_stock.bin"), 0, CFT_FP128, NULL, &st);
    if (!r){
        printf("gate 3: FAIL  load: %s\n", cft_archive_status_str(st));
        return 1;
    }
    printf("  load said: %s (%d)\n", cft_archive_status_str(st), (int)st);
    if (st != CFT_ARCHIVE_PROMOTED){ printf("gate 3: FAIL  expected CFT_ARCHIVE_PROMOTED\n"); fail = 1; }
    if (!r->integrator.name || strcmp(r->integrator.name, CFT_IAS15_INTEGRATOR_NAME)){
        printf("gate 3: FAIL  integrator after promotion is \"%s\"\n",
               r->integrator.name ? r->integrator.name : "(none)");
        fail = 1;
    }

    struct cft_ias15_state *cs = (struct cft_ias15_state*)r->integrator.state;
    size_t N3 = 3 * r->N;
    if (cs->n_elem != N3 || cs->format != CFT_FP128){
        printf("gate 3: FAIL  promoted state is n_elem %zu format %d, wanted %zu / %d\n",
               cs->n_elem, cs->format, N3, CFT_FP128);
        fail = 1;
    }
    /* The return value said "promoted from binary64" and the STATE has
     * to say so too: a caller may discard a return value, and a state
     * that was narrowed to binary64 and one restored bit for bit are
     * otherwise indistinguishable from the inside. ROADMAP.md line 404. */
    if (cs->provenance != CFT_PROV_PROMOTED){
        printf("gate 3: FAIL  the load returned %s and the state's provenance is \"%s\"\n",
               cft_archive_status_str(st), cft_ias15_provenance_str(cs->provenance));
        fail = 1;
    }
    /* A promoted state is at its own mark by definition. */
    if (cs->hiwater_n_elem){
        printf("gate 3: FAIL  a promoted state carries cft_hiwater_n_elem %zu\n",
               cs->hiwater_n_elem);
        fail = 1;
    }

    /* Promotion must be EXACT: binary64 is a subset of binary128, so
     * rounding the wide bytes back must give the very same doubles. */
    double *back = malloc(7 * N3 * sizeof(double));
    int bad = 0;
    struct { const char *name; unsigned char *wide; const double *ref; size_t n; } chk[] = {
        { "x0",   cs->x0,   NULL, N3 },
        { "v0",   cs->v0,   NULL, N3 },
        { "a0",   cs->a0,   ia->a0,   N3 },
        { "csx",  cs->csx,  ia->csx,  N3 },
        { "csv",  cs->csv,  ia->csv,  N3 },
        { "csa0", cs->csa0, ia->csa0, N3 },
    };
    double *pos = malloc(N3 * sizeof(double)), *vel = malloc(N3 * sizeof(double));
    for (size_t i = 0; i < r->N; i++){
        pos[3*i] = ref->particles[i].x; pos[3*i+1] = ref->particles[i].y; pos[3*i+2] = ref->particles[i].z;
        vel[3*i] = ref->particles[i].vx; vel[3*i+1] = ref->particles[i].vy; vel[3*i+2] = ref->particles[i].vz;
    }
    chk[0].ref = pos; chk[1].ref = vel;
    for (size_t k = 0; k < sizeof chk / sizeof chk[0]; k++){
        cft_archive_demote_doubles(NULL, cs->format, chk[k].wide, back, chk[k].n);
        for (size_t i = 0; i < chk[k].n; i++)
            if (memcmp(&back[i], &chk[k].ref[i], sizeof(double))){
                printf("  %s[%zu]: %a promoted, %a back\n", chk[k].name, i, chk[k].ref[i], back[i]);
                bad = 1;
            }
    }
    struct { const char *name; unsigned char **wide; const double *ref; } lv[] = {
        { "g",   cs->g,   ia->g   }, { "b",  cs->b,  ia->b  }, { "csb", cs->csb, ia->csb },
        { "e",   cs->e,   ia->e   }, { "br", cs->br, ia->br }, { "er",  cs->er,  ia->er  },
    };
    for (size_t k = 0; k < sizeof lv / sizeof lv[0]; k++)
        for (int j = 0; j < 7; j++){
            cft_archive_demote_doubles(NULL, cs->format, lv[k].wide[j], back, N3);
            for (size_t i = 0; i < N3; i++)
                if (memcmp(&back[i], &lv[k].ref[(size_t)j * N3 + i], sizeof(double))){
                    printf("  %s%d[%zu] differs\n", lv[k].name, j, i);
                    bad = 1;
                }
        }
    if (bad) fail = 1;
    printf("  promotion of REBOUND's binary64 IAS15 state into %s: %s "
           "(epsilon %g, adaptive_mode %d, abi \"%s\", provenance \"%s\")\n",
           cft_format_name((cft_format)cs->format),
           bad ? "NOT EXACT" : "exact, every value round-trips",
           cs->epsilon, cs->adaptive_mode, cs->cft_abi,
           cft_ias15_provenance_str(cs->provenance));
    free(back); free(pos); free(vel);
    reb_simulation_free(ref);
    reb_simulation_free(r);

    /* ---- the loud refusal ------------------------------------------ */
    remove(path("gate3_fp256.bin"));
    struct reb_simulation *w = reb_simulation_create();
    w->G = 1.0; w->dt = 0.01; w->exact_finish_time = 0;
    add_bodies(w);
    if (cft_shim_setup(w, CFT_FP256)) return 2;
    reb_simulation_steps(w, 3);
    reb_simulation_save_to_file(w, path("gate3_fp256.bin"));
    reb_simulation_free(w);

    printf("  a binary256 archive loaded into a binary64 run says:\n");
    fflush(stdout);
    enum cft_archive_status st2;
    struct reb_simulation *bad_r = cft_archive_load(path("gate3_fp256.bin"), 0, CFT_FP64, NULL, &st2);
    fflush(stderr);
    if (bad_r || st2 != CFT_ARCHIVE_FORMAT_MISMATCH){
        printf("gate 3: FAIL  a binary256 archive was accepted into a binary64 run (%s)\n",
               cft_archive_status_str(st2));
        if (bad_r) reb_simulation_free(bad_r);
        fail = 1;
    }else{
        printf("  refused with %s and a NULL simulation, as it should be\n",
               cft_archive_status_str(st2));
    }
    /* and the same file, loaded at the format it was written with */
    enum cft_archive_status st3;
    struct reb_simulation *good_r = cft_archive_load(path("gate3_fp256.bin"), 0, CFT_FP256, NULL, &st3);
    if (!good_r || st3 != CFT_ARCHIVE_EXACT){
        printf("gate 3: FAIL  the same archive at its own format: %s\n", cft_archive_status_str(st3));
        fail = 1;
    }else{
        printf("  the same file at binary256: %s\n", cft_archive_status_str(st3));
        reb_simulation_free(good_r);
    }

    printf("gate 3: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}
