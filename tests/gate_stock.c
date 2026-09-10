/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Gate 2, part 2: STOCK REBOUND opens the cft archive.
 *
 * This program links the pinned upstream librebound and nothing of
 * cft-rebound's: no integrator is registered, no cft_ field descriptor
 * exists in this process. That is exactly the position a collaborator
 * without cft is in. It reports every warning bit REBOUND raised and
 * compares the recovered binary64 state, bit for bit, with what the
 * writer recorded - for a single-snapshot archive and for the last
 * snapshot of a three-snapshot one.
 *
 * Usage: gate_stock [dir]
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "rebound.h"
#include "binarydata.h"
#include "simulationarchive.h"

static const struct { int bit; const char *name; } warnbits[] = {
    { REB_BINARYDATA_ERROR_NOFILE,             "ERROR_NOFILE" },
    { REB_BINARYDATA_WARNING_VERSION,          "WARNING_VERSION" },
    { REB_BINARYDATA_WARNING_POINTERS,         "WARNING_POINTERS" },
    { REB_BINARYDATA_WARNING_PARTICLES,        "WARNING_PARTICLES" },
    { REB_BINARYDATA_ERROR_FILENOTOPEN,        "ERROR_FILENOTOPEN" },
    { REB_BINARYDATA_ERROR_OUTOFRANGE,         "ERROR_OUTOFRANGE" },
    { REB_BINARYDATA_ERROR_SEEK,               "ERROR_SEEK" },
    { REB_BINARYDATA_WARNING_FIELD_UNKNOWN,    "WARNING_FIELD_UNKNOWN" },
    { REB_BINARYDATA_ERROR_INTEGRATOR,         "ERROR_INTEGRATOR" },
    { REB_BINARYDATA_WARNING_CORRUPTFILE,      "WARNING_CORRUPTFILE" },
    { REB_BINARYDATA_ERROR_OLD,                "ERROR_OLD" },
    { REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR,"WARNING_CUSTOM_INTEGRATOR" },
    { 0, NULL }
};

static int one(const char *bin, const char *exp, int64_t snapshot, const char *what){
    enum REB_BINARYDATA_ERROR_CODE warnings = REB_BINARYDATA_WARNING_NONE;
    struct reb_simulationarchive *sa = reb_simulationarchive_create_from_file_with_messages(bin, &warnings);
    if (!sa || !sa->inf){
        printf("  %s: FAIL  stock REBOUND could not open %s\n", what, bin);
        return 1;
    }
    int64_t nblobs = sa->nblobs;
    struct reb_simulation *r = reb_simulation_create();
    reb_simulation_init_from_simulationarchive_with_messages(r, sa, snapshot, &warnings);
    reb_simulationarchive_free(sa);

    printf("  %s: %lld snapshots indexed; warnings = 0x%x:", what,
           (long long)nblobs, (unsigned)warnings);
    if (!warnings) printf(" (none)");
    for (int i = 0; warnbits[i].name; i++)
        if (warnings & warnbits[i].bit) printf(" %s", warnbits[i].name);
    printf("\n    integrator.name in this reader after the load: \"%s\"\n",
           r->integrator.name ? r->integrator.name : "(none)");

    FILE *f = fopen(exp, "rb");
    if (!f){ printf("  %s: FAIL  missing %s\n", what, exp); reb_simulation_free(r); return 1; }
    size_t N = 0;
    int bad = 0;
    if (fread(&N, sizeof N, 1, f) != 1){ fclose(f); return 1; }
    if (N != r->N){ printf("    N: expected %zu, got %zu\n", N, r->N); bad = 1; }
    for (size_t i = 0; i < N && i < r->N; i++){
        double v[7];
        if (fread(v, sizeof v, 1, f) != 1){ bad = 1; break; }
        double g[7] = { r->particles[i].x, r->particles[i].y, r->particles[i].z,
                        r->particles[i].vx, r->particles[i].vy, r->particles[i].vz,
                        r->particles[i].m };
        for (int k = 0; k < 7; k++)
            if (memcmp(&v[k], &g[k], sizeof(double)) != 0){
                printf("    particle %zu field %d: expected %a, got %a\n", i, k, v[k], g[k]);
                bad = 1;
            }
    }
    double t, dt;
    if (fread(&t, sizeof t, 1, f) == 1 && fread(&dt, sizeof dt, 1, f) == 1){
        if (memcmp(&t, &r->t, sizeof t)) { printf("    t: expected %a, got %a\n", t, r->t); bad = 1; }
        if (memcmp(&dt, &r->dt, sizeof dt)) { printf("    dt: expected %a, got %a\n", dt, r->dt); bad = 1; }
    }else bad = 1;
    fclose(f);

    printf("    binary64 state recovered: %s (%zu particles, t = %a)\n",
           bad ? "WRONG" : "bit for bit as written", r->N, r->t);
    reb_simulation_free(r);
    return bad;
}

int main(int argc, char **argv){
    const char *dir = (argc > 1) ? argv[1] : ".";
    char bin[1024], exp[1024], bin2[1024], exp2[1024];
    snprintf(bin,  sizeof bin,  "%s/gate2_cft.bin", dir);
    snprintf(exp,  sizeof exp,  "%s/gate2_expect.bin", dir);
    snprintf(bin2, sizeof bin2, "%s/gate2_cft_multi.bin", dir);
    snprintf(exp2, sizeof exp2, "%s/gate2_expect_multi.bin", dir);

    printf("gate 2: a cft archive opened by a stock REBOUND reader\n");
    int bad = one(bin, exp, 0, "one snapshot");
    bad |= one(bin2, exp2, -1, "last of three");
    printf("gate 2: %s\n", bad ? "FAIL" : "PASS");
    return bad;
}
