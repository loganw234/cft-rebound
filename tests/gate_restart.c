/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Gate 1: write an archive mid-run, restart from it, and continue
 * bit-identically to an uninterrupted run.
 *
 * Two shapes, at every format:
 *   a) a fresh single-snapshot archive written after 10 of 20 steps
 *      (the checkpoint case);
 *   b) a two-snapshot archive, the second appended - which REBOUND
 *      stores as a DIFF against the first, so the wide blobs travel
 *      through reb_binarydata_diff as well as through the writer.
 *
 * Usage: gate_restart [dir]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cft_shim_stub.h"

static const char *dir = ".";

static char *path(const char *stem){
    static char buf[1024];
    snprintf(buf, sizeof buf, "%s/%s", dir, stem);
    return buf;
}

static struct reb_simulation *mk(int format){
    struct reb_simulation *r = reb_simulation_create();
    r->G = 1.0;
    r->dt = 0.01;
    r->exact_finish_time = 0;
    struct reb_particle p = {0};
    p.m = 1.0;   p.x =  1.0; p.y = 0.25; p.z = -0.5;  p.vx = 0.0;  p.vy = 0.75; p.vz = 0.125;
    reb_simulation_add(r, p);
    p.m = 1e-3;  p.x = -0.5; p.y = 1.5;  p.z =  0.25; p.vx = -0.5; p.vy = 0.25; p.vz = -0.75;
    reb_simulation_add(r, p);
    p.m = 1e-4;  p.x =  2.0; p.y = -1.0; p.z =  0.75; p.vx = 0.25; p.vy = -0.5; p.vz = 0.5;
    reb_simulation_add(r, p);
    if (cft_shim_setup(r, format)){ fprintf(stderr, "setup failed\n"); exit(2); }
    return r;
}

static int fail = 0;

static void run_format(int format){
    const char *fn = cft_format_name((cft_format)format);

    /* the uninterrupted run */
    struct reb_simulation *a = mk(format);
    reb_simulation_steps(a, 20);
    size_t la; unsigned char *sa = cft_shim_snapshot(a, &la);
    int nz = cft_shim_blobs_nonzero(a);

    /* (a) one snapshot, written after 10 steps */
    remove(path("gate1_one.bin"));
    struct reb_simulation *b = mk(format);
    reb_simulation_steps(b, 10);
    reb_simulation_save_to_file(b, path("gate1_one.bin"));
    reb_simulation_free(b);

    enum cft_archive_status st;
    b = cft_archive_load(path("gate1_one.bin"), 0, format, NULL, &st);
    if (!b){
        printf("  %-9s FAIL  load: %s\n", fn, cft_archive_status_str(st));
        fail = 1; free(sa); reb_simulation_free(a); return;
    }
    reb_simulation_steps(b, 10);
    size_t lb; unsigned char *sb = cft_shim_snapshot(b, &lb);

    /* (b) two snapshots, the second an appended diff */
    remove(path("gate1_two.bin"));
    struct reb_simulation *c = mk(format);
    reb_simulation_steps(c, 5);
    reb_simulation_save_to_file(c, path("gate1_two.bin"));
    reb_simulation_steps(c, 5);
    reb_simulation_save_to_file(c, path("gate1_two.bin"));
    reb_simulation_free(c);

    enum cft_archive_status st2;
    c = cft_archive_load(path("gate1_two.bin"), -1, format, NULL, &st2);
    if (!c){
        printf("  %-9s FAIL  load(-1): %s\n", fn, cft_archive_status_str(st2));
        fail = 1; free(sa); free(sb); reb_simulation_free(a); reb_simulation_free(b); return;
    }
    reb_simulation_steps(c, 10);
    size_t lc; unsigned char *sc = cft_shim_snapshot(c, &lc);

    struct cft_archive_info info;
    cft_archive_probe(path("gate1_two.bin"), -1, &info);

    int ok1 = (la == lb) && memcmp(sa, sb, la) == 0;
    int ok2 = (la == lc) && memcmp(sa, sc, la) == 0;
    printf("  %-9s %s  %zu state bytes, %d/48 blobs non-zero, %d cft_ fields, "
           "n_elem %llu, %llu B a blob, %lld snapshots; one-snapshot %s, appended-diff %s (%s, %s)\n",
           fn, (ok1 && ok2) ? "PASS" : "FAIL", la, nz, info.n_cft_fields,
           (unsigned long long)info.n_elem, (unsigned long long)info.blob_bytes,
           (long long)info.n_snapshots,
           ok1 ? "identical" : "DIFFERS", ok2 ? "identical" : "DIFFERS",
           cft_archive_status_str(st), cft_archive_status_str(st2));
    if (!ok1 || !ok2 || nz != 48) fail = 1;
    if (nz != 48) printf("  %-9s FAIL  only %d of 48 blobs carry a non-zero byte\n", fn, nz);

    free(sa); free(sb); free(sc);
    reb_simulation_free(a); reb_simulation_free(b); reb_simulation_free(c);
}

int main(int argc, char **argv){
    if (argc > 1) dir = argv[1];
    cft_shim_register();
    printf("gate 1: mid-run archive, restart, bit-identical continuation\n");
    int shape = cft_archive_selftest();
    printf("  descriptor lists: %s\n", shape ?
           "MALFORMED" :
           "48 blobs + 11 scalars at each of fp64/fp128/fp256, every name cft_-prefixed and unique, cft_n_elem last");
    if (shape) fail = 1;
    run_format(CFT_FP64);
    run_format(CFT_FP128);
    run_format(CFT_FP256);
    printf("gate 1: %s\n", fail ? "FAIL" : "PASS");
    return fail;
}
