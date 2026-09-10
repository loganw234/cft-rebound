/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Gate 2, part 1: write a cft archive at binary128 and record, exactly,
 * the binary64 view a stock REBOUND reader must come back with.
 *
 * Usage: gate_write [dir]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cft_shim_stub.h"

static void expect(const char *path, struct reb_simulation *r){
    FILE *f = fopen(path, "wb");
    if (!f) exit(2);
    size_t N = r->N;
    fwrite(&N, sizeof N, 1, f);
    for (size_t i = 0; i < N; i++){
        double v[7] = { r->particles[i].x, r->particles[i].y, r->particles[i].z,
                        r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
                        r->particles[i].m };
        fwrite(v, sizeof v, 1, f);
    }
    fwrite(&r->t, sizeof(double), 1, f);
    fwrite(&r->dt, sizeof(double), 1, f);
    fclose(f);
}

int main(int argc, char **argv){
    const char *dir = (argc > 1) ? argv[1] : ".";
    char bin[1024], exp[1024], bin2[1024], exp2[1024];
    snprintf(bin,  sizeof bin,  "%s/gate2_cft.bin", dir);
    snprintf(exp,  sizeof exp,  "%s/gate2_expect.bin", dir);
    snprintf(bin2, sizeof bin2, "%s/gate2_cft_multi.bin", dir);
    snprintf(exp2, sizeof exp2, "%s/gate2_expect_multi.bin", dir);

    cft_shim_register();
    struct reb_simulation *r = reb_simulation_create();
    r->G = 1.0;
    r->dt = 0.01;
    r->exact_finish_time = 0;
    struct reb_particle p = {0};
    p.m = 1.0;  p.x =  1.0; p.y = 0.25; p.z = -0.5;  p.vx = 0.0;  p.vy = 0.75; p.vz = 0.125;
    reb_simulation_add(r, p);
    p.m = 1e-3; p.x = -0.5; p.y = 1.5;  p.z =  0.25; p.vx = -0.5; p.vy = 0.25; p.vz = -0.75;
    reb_simulation_add(r, p);
    p.m = 1e-4; p.x =  2.0; p.y = -1.0; p.z =  0.75; p.vx = 0.25; p.vy = -0.5; p.vz = 0.5;
    reb_simulation_add(r, p);
    if (cft_shim_setup(r, CFT_FP128)) return 2;

    reb_simulation_steps(r, 7);

    remove(bin);
    reb_simulation_save_to_file(r, bin);
    expect(exp, r);

    /* and a three-snapshot archive, so the stock reader has to index a
     * file whose blobs carry unknown fields and read the last one, which
     * REBOUND stores as a diff. */
    remove(bin2);
    reb_simulation_save_to_file(r, bin2);
    reb_simulation_steps(r, 4);
    reb_simulation_save_to_file(r, bin2);
    reb_simulation_steps(r, 4);
    reb_simulation_save_to_file(r, bin2);
    expect(exp2, r);

    struct cft_archive_info info, info2;
    if (cft_archive_probe(bin, 0, &info)) return 2;
    if (cft_archive_probe(bin2, -1, &info2)) return 2;
    printf("gate 2 (write): %s at %s, integrator \"%s\", %d cft_ fields, "
           "%d/48 blobs, %llu B a blob, n_elem %llu, abi \"%s\"\n",
           bin, cft_format_name(CFT_FP128), info.integrator, info.n_cft_fields,
           info.n_blobs_seen, (unsigned long long)info.blob_bytes,
           (unsigned long long)info.n_elem, info.abi);
    printf("                %s, %lld snapshots, last one probes to n_elem %llu at %s\n",
           bin2, (long long)info2.n_snapshots, (unsigned long long)info2.n_elem,
           cft_format_name((cft_format)info2.format));
    reb_simulation_free(r);
    return 0;
}
