/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * cft-rebound as a REBOUND drop-in.
 *
 * REBOUND accepts user-provided integrators through a public API, so
 * this is a registration and not a patch. A user writes their simulation
 * exactly as they do now and changes one line:
 *
 *     #include "cft_ias15.h"
 *     ...
 *     cft_ias15_register("ias15_cft");
 *     struct cft_ias15_state* st = reb_simulation_set_integrator(r, "ias15_cft");
 *     st->format = CFT_FP128;          // optional; the default is CFT_FP64
 *
 * Everything else in REBOUND - particles, outputs, callbacks, the
 * visualisation - keeps working, because none of it is replaced.
 *
 * THE DIVISION OF STATE. struct reb_particle holds double x, y, z, vx,
 * vy, vz and always will, so a wide integration cannot live in it:
 *
 *   - REBOUND's particles are the binary64 VIEW, correctly rounded from
 *     the wide state after every step;
 *   - the wide state is the integrator's own, in struct cft_ias15_state.
 *
 * At CFT_FP64 the two are the same bits and the view costs nothing,
 * which is what the equivalence gate (tools/check_dropin.c) checks: the
 * same program with "ias15" and with "ias15_cft" must agree on every
 * particle value bit for bit.
 *
 * WHAT IS REFUSED. Additional forces, velocity-dependent forces,
 * collisions, ghost boxes, non-trivial boundaries, any gravity module
 * but REB_GRAVITY_BASIC, non-zero softening, test particles (N_active),
 * r->map subsets, variational particles and MEGNO are out of scope. The
 * step detects each of them, names it in a REBOUND error message and
 * stops the integration rather than compute something wrong.
 */
#ifndef CFT_IAS15_H
#define CFT_IAS15_H

#include <stddef.h>
#include <stdint.h>
#include "cft.h"

struct reb_simulation;
struct reb_integrator;

#ifdef __cplusplus
extern "C" {
#endif

/* --------------------------------------------------------------------
 * The state, as ROADMAP.md defines it.
 *
 * This struct is the contract between the integrator (parcel A) and the
 * Simulationarchive support (parcel B): parcel B writes
 * cft_ias15_field_descriptor_list against these offsets. The fields
 * below, their names, their types and their order are ROADMAP.md's, and
 * nothing may be inserted among them.
 *
 * W is the wide element width in bytes - cft_format_size(format), 8, 16
 * or 32. The byte blobs hold n_elem elements each; the seven-fold arrays
 * hold n_elem elements per level.
 * -------------------------------------------------------------------- */
struct cft_ias15_state {
    /* configuration, plain doubles, archived as REB_DOUBLE */
    double   epsilon;          /* as REBOUND's */
    double   min_dt;
    int      adaptive_mode;
    int      format;           /* CFT_FP64 | CFT_FP128 | CFT_FP256 */
    int      max_iter;         /* 12 in REBOUND; binary256 needs ~22 */
    int      arith_fma;        /* 0 = REBOUND's roundings, 1 = the FMA form */

    /* the wide state: byte blobs, archived as REB_POINTER with
     * element_size = W. Lengths are 3N or 3N*E elements. */
    unsigned char *x, *v;                  /* the live coordinates: what a
                                            * restart continues FROM. x0/v0/a0
                                            * are the step's starting copy and
                                            * are refreshed from these at the
                                            * top of every step. */
    unsigned char *x0, *v0, *a0;      /* position, velocity, acceleration */
    unsigned char *csx, *csv, *csa0;  /* compensated-summation carries */
    unsigned char *g[7], *b[7], *e[7], *br[7], *er[7], *csb[7];
    size_t   n_elem;           /* 3N, or 3NE for an ensemble */
    size_t   E;                /* ensemble members, 1 for a plain sim */

    /* provenance, so a reader knows what produced the wide bytes */
    char     cft_abi[16];      /* e.g. "0.11" */
    uint64_t constants_digest; /* the Gauss-Radau set the run used */
};

/* --------------------------------------------------------------------
 * The integrator
 * -------------------------------------------------------------------- */

/* The callbacks. Pass this to reb_integrator_register() yourself, or
 * call cft_ias15_register() below, which does it once. */
extern const struct reb_integrator cft_ias15_integrator;

/* Register under `name` (NULL means "ias15_cft"). Registering the same
 * name twice is a REBOUND error, so this is a no-op after the first
 * call with the same name.
 *
 * ONE CUSTOM INTEGRATOR PER PROCESS, and that is upstream's limit, not
 * this one. rebound.c's reb_integrator_register scans the existing list
 * with
 *
 *     while (list[N].name){ N++; if (strcmp(list[N].name, name)==0) ... }
 *
 * which reads list[N].name AFTER the increment, so a SECOND registration
 * reaches the {0} terminator and passes its NULL name to strcmp. On this
 * host (REBOUND bdfda4bd, mingw64) the call never returns. So do not
 * register another custom integrator alongside this one until upstream
 * moves the increment. */
void cft_ias15_register(const char *name);

/* The name this integrator is registered under by everything in this
 * repository, and the one the Simulationarchive's field prefix is
 * built from. Lower case is not a style choice: REBOUND's Python layer
 * lower-cases what you assign to sim.integrator before it reaches the
 * C strcmp, so a name with a capital in it registers fine and can
 * never be selected. */
#define CFT_IAS15_INTEGRATOR_NAME "ias15_cft"

/* The state of a simulation using this integrator, or NULL if it is
 * using another one. Same pointer reb_simulation_set_integrator returns. */
struct cft_ias15_state *cft_ias15_get_state(struct reb_simulation *r);

/* --------------------------------------------------------------------
 * Two things the state struct cannot carry
 *
 * The engine in src/ias15_cft.c is one global instance - every buffer in
 * it is a file-scope static - so the settings below are the process's,
 * not a simulation's, and they must be made before the first step.
 * -------------------------------------------------------------------- */

/* Reserve room for up to n bodies. Several of the step's scratch vectors
 * are allocated at their first use and never resized, so a simulation
 * that will grow past its starting N must say so before it starts. The
 * default is the N of the first step; a later reb_simulation_add() past
 * the reserved count is refused, not silently mis-sized. */
void cft_ias15_reserve(size_t n);

/* The libcft artifact to open, or NULL (the default) for the software
 * backend. Not exercised in this parcel. */
void cft_ias15_set_artifact(const char *path);

/* The predictor-corrector tolerance exponent: the loop stops under
 * 1e-16 * 2^-shift. -1 (the default) means REBOUND's 1e-16 at binary64
 * and the same number of ulps in the wider formats (53-p). */
void cft_ias15_set_pc_tol_shift(int shift);

/* Diagnostics from the engine, for a caller that wants them. */
uint32_t           cft_ias15_flags_seen(void);
unsigned long long cft_ias15_library_calls(void);

/* --------------------------------------------------------------------
 * Parcel B's symbol
 *
 * The archive field descriptors. Parcel A ships a placeholder holding
 * only the terminator (src/cft_ias15_fields.c), which REBOUND reads as
 * "this integrator adds no fields"; parcel B replaces that file. The
 * names it defines are prefixed cft_ and reach the file as
 * "integrator.<name>.cft_<field>", REBOUND prefixing them itself.
 * -------------------------------------------------------------------- */
struct reb_binarydata_field_descriptor;
extern const struct reb_binarydata_field_descriptor cft_ias15_field_descriptor_list[];

#ifdef __cplusplus
}
#endif
#endif /* CFT_IAS15_H */
