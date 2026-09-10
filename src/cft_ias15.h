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
 * stops the integration rather than compute something wrong. So are
 * four of the state's own settings: an adaptive_mode other than PRS23
 * (2), a non-zero min_dt, a format that is not one of the three, and a
 * max_iter below 0 (0 itself means the default for the format); and an
 * E other than 1, because an ensemble is E
 * independent systems and a reb_simulation is one (docs/ENSEMBLE.md).
 *
 * Three things the subprocess API refuses and this does NOT check:
 * pre_/post_timestep_modifications (deliberate - the step re-promotes
 * a coordinate a callback edited), r->gravity_custom, and r->N_odes.
 * The README's scope table says so beside the table.
 */
#ifndef CFT_IAS15_H
#define CFT_IAS15_H

#include <stddef.h>
#include <stdint.h>
#include "rebound.h"       /* and binarydata.h below needs its REB_API
                            * macro, so the order here is not cosmetic */
#include "binarydata.h"    /* struct reb_binarydata_field_descriptor, which
                            * the field list below is an array of - so this
                            * header compiles on its own rather than only
                            * after whatever included it first */
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
 * Simulationarchive support (parcel B): src/cft_ias15_fields.c writes
 * the descriptor lists against these offsets with offsetof, so a field
 * may be added but the meaning of an existing one may not change under
 * an archive already written.
 *
 * ROADMAP.md's first version of this struct was wrong in one place and
 * both parcels inherited it: it had no x and v, so a checkpoint carried
 * the previous step's starting copy instead of the live coordinates and
 * silently truncated position and velocity to binary64 above binary64.
 * They are here now and the blob count is 50 (CFT_N_BLOBS in
 * src/cft_ias15_fields.h). docs/VALIDATION.md entry 25.
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
    int      max_iter;         /* the corrector's pass cap. 0 means the
                                * default for this format, resolved on the
                                * first step and written back here; see
                                * cft_ias15_default_max_iter() */
    int      arith_fma;        /* 0 = REBOUND's roundings, 1 = the FMA form */

    /* the wide state: byte blobs, archived as REB_POINTER with
     * element_size = W. Lengths are 3N or 3N*E elements. */
    unsigned char *x, *v;                  /* the live coordinates: what a
                                            * restart continues FROM. x0/v0/a0
                                            * are the step's starting copy and
                                            * are refreshed from these at the
                                            * top of every step. */
    unsigned char *x0, *v0, *a0;      /* NOT the live position, velocity and
                                       * acceleration, whatever their names
                                       * suggest: step_attempt() sets these
                                       * FROM x/v at the top of a step, so
                                       * after a completed step x0 is the
                                       * PREVIOUS step's start. Calling them
                                       * "position, velocity, acceleration"
                                       * here is what made two parcels archive
                                       * them and not x/v - VALIDATION 25. */
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

/* The same settings by value, for a caller that holds a simulation and
 * not the state - which is every ctypes caller. REBOUND resolves
 * sim.integrator.<field> against the registered field_descriptor_list,
 * and every name in ours is cft_-prefixed, so the spellings a REBOUND
 * user knows (sim.integrator.epsilon) do not resolve here. Without
 * these, Python can select the integrator and never reach the format.
 *
 *   format    CFT_FP64 | CFT_FP128 | CFT_FP256
 *   epsilon   as REBOUND's: 0 is a fixed step
 *   max_iter  the corrector's pass cap; 0 means the default for the
 *             format, which is REBOUND's 12 at binary64 and more above
 *             it because the corrector needs about 6 passes at
 *             binary128 and 18-22 at binary256 - see docs/VALIDATION.md
 *
 * Returns 0, or non-zero if r is not using this integrator or a value
 * is not one this integrator accepts. Settings that are refused leave
 * the state unchanged; nothing is applied partially. */
int cft_ias15_configure(struct reb_simulation *r, int format,
                        double epsilon, int max_iter);

/* "fp64", "fp128", "fp256" to the matching CFT_FP* value, or -1.
 * So a ctypes caller need not hard-code an enumerator that belongs to
 * libcft and could be renumbered there. */
int cft_ias15_format_code(const char *name);

/* The corrector's pass cap for a format when the caller has not chosen
 * one: 12 at binary64, which is REBOUND's own and part of what makes
 * the binary64 run bit-identical to it, and more above it because the
 * corrector needs more passes there. state->max_iter of 0 means this,
 * resolved on the first step and written back. */
int cft_ias15_default_max_iter(int format);

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
 * backend. An explicit setting here wins; otherwise $CFT_REBOUND_ARTIFACT
 * names one, which is what lets the gate suite reach a card without a
 * flag of its own. Empty is treated as unset. The bits are the same
 * either way - docs/VALIDATION.md entry 27 ran check_dropin and
 * gate_real on a U50C quad tile at every format and diffed the output
 * against the software backend's: identical. */
void cft_ias15_set_artifact(const char *path);

/* The predictor-corrector tolerance exponent: the loop stops under
 * 1e-16 * 2^-shift. -1 (the default) means REBOUND's 1e-16 at binary64
 * and the same number of ulps in the wider formats (53-p). */
void cft_ias15_set_pc_tol_shift(int shift);

/* Diagnostics from the engine, for a caller that wants them. */
uint32_t           cft_ias15_flags_seen(void);
unsigned long long cft_ias15_library_calls(void);

/* --------------------------------------------------------------------
 * The archive field descriptors
 *
 * Defined once, in src/cft_ias15_fields.c, from the macros in
 * src/cft_ias15_fields.h. The names are prefixed cft_ and reach the file
 * as "integrator.<name>.cft_<field>", REBOUND prefixing them itself.
 *
 * This file held only a terminator until 2026-09-10 - "parcel A's
 * placeholder, for parcel B to replace" - while parcel B kept the real
 * lists as file-scope statics in src/cft_archive.c, which is not linked
 * into every drop-in target. REBOUND therefore resolved no cft_ field on
 * read and a checkpoint came back from create() defaults. See
 * docs/VALIDATION.md entry 25 and tests/gate_real.c, the gate that
 * caught it.
 * -------------------------------------------------------------------- */
struct reb_binarydata_field_descriptor;
extern const struct reb_binarydata_field_descriptor cft_ias15_field_descriptor_list[];

#ifdef __cplusplus
}
#endif
#endif /* CFT_IAS15_H */
