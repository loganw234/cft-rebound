/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The refusal list, as data. src/cft_supported.h says why it is data.
 *
 * The force model is basic pairwise Newtonian gravity, full stop.
 * Everything here computes something DIFFERENT from what the port
 * would compute, so each is refused by name rather than silently
 * ignored or approximated.
 *
 * A row's message carries no prefix and no trailing newline. The
 * drop-in adds "ias15_cft: "; the subprocess API's caller prints it
 * however it likes.
 */
#include <stdio.h>
#include <string.h>

#include "cft_supported.h"
#include "cft.h"

/* Every row is a predicate and a sentence. They are written as pairs
 * rather than as one function returning a message because the walk
 * asks the question far more often than it answers it: a supported
 * simulation evaluates every predicate and formats nothing. */

/* ---- the simulation ---------------------------------------------- */

static int hit_gravity_custom(const struct cft_support_ctx *c){
    return c->r->gravity_custom != NULL;
}
static void say_gravity_custom(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "a custom gravity routine is not supported. The engine issues "
                   "REBOUND's basic pairwise gravity itself, so r->gravity_custom "
                   "would never be called and the run would answer a different "
                   "problem.");
}

static int hit_variational(const struct cft_support_ctx *c){
    return c->r->N_var != 0 || c->r->particles_var != NULL;
}
static void say_variational(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "variational particles are not supported (r->N_var = %zu). "
                   "The wide state carries the real particles only.", c->r->N_var);
}

static int hit_megno(const struct cft_support_ctx *c){
    return c->r->calculate_megno != 0;
}
static void say_megno(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "MEGNO is not supported; it needs variational particles.");
}

static int hit_additional_forces(const struct cft_support_ctx *c){
    return c->r->additional_forces != NULL;
}
static void say_additional_forces(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "r->additional_forces is not supported. The engine issues "
                   "REBOUND's basic pairwise gravity and nothing else.");
}

/* Refused on the subprocess path and NOT on the drop-in, which is not
 * an oversight: the drop-in's step re-promotes any coordinate that
 * changed under it, which is how a callback that edits a particle is
 * meant to work. A REBOUNDx FORCE still reaches the drop-in unrefused
 * and uncomputed, which is the open item README names. */
static int hit_timestep_mods(const struct cft_support_ctx *c){
    return c->r->pre_timestep_modifications != NULL ||
           c->r->post_timestep_modifications != NULL;
}
static void say_timestep_mods(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "timestep modifications are not supported (REBOUNDx and the like).");
}

static int hit_veldep(const struct cft_support_ctx *c){
    return c->r->force_is_velocity_dependent != 0;
}
static void say_veldep(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "velocity-dependent forces are not supported "
                   "(r->force_is_velocity_dependent = 1).");
}

static int hit_gravity_module(const struct cft_support_ctx *c){
    return c->r->gravity != REB_GRAVITY_BASIC;
}
static void say_gravity_module(const struct cft_support_ctx *c, char *b, size_t n){
    const char *which =
        c->r->gravity == REB_GRAVITY_TREE
            ? "the tree code; this port is direct summation"
            : (c->r->gravity == REB_GRAVITY_COMPENSATED
                ? "COMPENSATED is a different summation from the one ported"
                : "the Jacobi module and a custom routine are not ported either");
    snprintf(b, n, "only REB_GRAVITY_BASIC is supported, r->gravity is %d (%s).",
             (int)c->r->gravity, which);
}

static int hit_ignore_terms(const struct cft_support_ctx *c){
    return c->r->gravity_ignore_terms != REB_GRAVITY_IGNORE_TERMS_NONE;
}
static void say_ignore_terms(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "r->gravity_ignore_terms must be NONE; this port computes every pair.");
}

/* Two behaviours, one capability, and the row keeps them together
 * rather than splitting into two rows the gate would have to know
 * about separately.
 *
 * The drop-in supports collisions: REBOUND's driver runs the search and
 * the resolver between steps, and the integrator's whole obligation is
 * to survive the removal the way REBOUND does - which
 * ias15_engine_alias_resize() reproduces, coefficient aliasing and all.
 * Above binary64 the removal shifts r->particles, which are binary64,
 * so accurate = 0 would re-promote the survivors from them and lose
 * every wide tail; that is refused rather than approximated.
 *
 * The subprocess API refuses any collision, because REBOUND's driver
 * does not run between ITS steps. */
static int hit_collision(const struct cft_support_ctx *c){
    if (c->r->collision == REB_COLLISION_NONE) return 0;
    if (!c->have_state) return 1;
    return c->format != CFT_FP64 && !c->accurate;
}
static void say_collision(const struct cft_support_ctx *c, char *b, size_t n){
    if (!c->have_state){
        snprintf(b, n, "collision detection is not supported by this API "
                       "(r->collision != NONE); the drop-in integrator supports it, "
                       "because REBOUND's own driver runs the search between steps "
                       "there and this one does not.");
        return;
    }
    snprintf(b, n, "collision detection at %s needs state->accurate = 1. A removal "
                   "shifts r->particles, which are binary64, so at accurate = 0 - "
                   "where this port reproduces REBOUND exactly, aliased coefficient "
                   "levels and all - the wide state is re-promoted from them and "
                   "every coordinate loses its tail. accurate = 1 shifts the wide "
                   "state instead, which is more accurate than REBOUND and therefore "
                   "not bit-identical to it.",
             cft_format_name((cft_format)c->format));
}

static int hit_ghost_boxes(const struct cft_support_ctx *c){
    return c->r->N_ghost_x || c->r->N_ghost_y || c->r->N_ghost_z;
}
static void say_ghost_boxes(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "ghost boxes are not supported (N_ghost = %d, %d, %d). The "
                   "engine's pair term uses a ghost-box offset of exactly +0.",
             c->r->N_ghost_x, c->r->N_ghost_y, c->r->N_ghost_z);
}

/* No check on r->OMEGA. It reaches the dynamics only through the
 * shearing sheet (boundary SHEAR) or the SEI integrator, both refused
 * here, and REBOUND initialises r->OMEGAZ to -1 as a sentinel meaning
 * "use OMEGA" - so a test on either refuses every default simulation,
 * which this list did until the worked example ran. */
static int hit_boundary(const struct cft_support_ctx *c){
    return c->r->boundary != REB_BOUNDARY_NONE;
}
static void say_boundary(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "only REB_BOUNDARY_NONE is supported, r->boundary is %d.",
             (int)c->r->boundary);
}

static int hit_particle_map(const struct cft_support_ctx *c){
    return c->r->map != NULL || c->r->N_map != 0;
}
static void say_particle_map(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "r->map (integrating a subset of the particles) is not supported.");
}

static int hit_test_particles(const struct cft_support_ctx *c){
    return c->r->N_active != (size_t)-1 && c->r->N_active != c->r->N;
}
static void say_test_particles(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "test particles are not supported (r->N_active = %zu of %zu). "
                   "Every particle in the engine's pair list is active.",
             c->r->N_active, c->r->N);
}

/* Refused on the subprocess path and NOT on the drop-in, and the
 * drop-in is right: reb_simulation_step() integrates attached ODE sets
 * itself, with a private Bulirsch-Stoer state, for every integrator
 * whose name is not "bs" (simulation.c, "Integrate other ODEs"). So
 * under the drop-in they are integrated exactly as they are under
 * REBOUND's own ias15. The subprocess API hands the problem to another
 * process, which never sees them. */
static int hit_odes(const struct cft_support_ctx *c){
    return c->r->N_odes != 0;
}
static void say_odes(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "attached ODE sets are not supported by this API.");
}

static int hit_no_particles(const struct cft_support_ctx *c){
    return c->r->N < 1;
}
static void say_no_particles(const struct cft_support_ctx *c, char *b, size_t n){
    (void)c;
    snprintf(b, n, "no particles.");
}

static int hit_body_count(const struct cft_support_ctx *c){
    return c->r->N > CFT_REBOUND_MAX_BODIES;
}
static void say_body_count(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "%lu particles, and the port's limit is %d per system "
                   "(memory and time grow as N^2).",
             (unsigned long)c->r->N, CFT_REBOUND_MAX_BODIES);
}

static int hit_integrator_name(const struct cft_support_ctx *c){
    return c->r->integrator.name != NULL &&
           strcmp(c->r->integrator.name, "ias15") != 0;
}
static void say_integrator_name(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "the simulation's integrator is \"%.80s\"; this port is IAS15 only.",
             c->r->integrator.name);
}

/* ---- the integrator's own settings -------------------------------- */

static int hit_ensemble(const struct cft_support_ctx *c){
    return c->have_state && c->E != 1;
}
static void say_ensemble(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "state->E is %zu. An ensemble is E independent systems in one run "
                   "and a reb_simulation is one system; use the standalone ias15_cft "
                   "program for ensembles (docs/ENSEMBLE.md).", c->E);
}

/* All four of REBOUND's criteria are implemented for the drop-in, so
 * what is left to refuse there is a number that is not one of them.
 *
 * The subprocess path is NOT widened with it, and that is the row
 * carrying a fact rather than an oversight: cft_rebound_steps() hands
 * the run to the standalone ias15_cft program and does not pass
 * --adaptive-mode, so whatever the simulation asks for, the program
 * runs its own default. Accepting 0 or 1 there would mean accepting
 * them and then silently running PRS23. (The same is already true of
 * AARSETH85; this row leaves that exactly as it found it, because
 * changing it is a change to cft_rebound_run.c's command line and not
 * to this table. PARCELS.md P3's report says so.) */
static int hit_adaptive_mode(const struct cft_support_ctx *c){
    if (c->have_state) return c->adaptive_mode < 0 || c->adaptive_mode > 3;
    return c->adaptive_mode != 2 && c->adaptive_mode != 3;
}
static void say_adaptive_mode(const struct cft_support_ctx *c, char *b, size_t n){
    if (!c->have_state){
        snprintf(b, n, "adaptive_mode %d is not supported by this API. This path runs "
                       "the standalone ias15_cft program, whose command line does not "
                       "carry the step criterion, so only the criteria that program "
                       "defaults to can be honoured; the drop-in integrator implements "
                       "all four of REBOUND's.", c->adaptive_mode);
        return;
    }
    snprintf(b, n, "adaptive_mode %d is not one of REBOUND's: INDIVIDUAL (0), "
                   "GLOBAL (1), PRS23 (2, its default since January 2024) or "
                   "AARSETH85 (3).", c->adaptive_mode);
}

static int hit_format(const struct cft_support_ctx *c){
    return c->format != CFT_FP64 && c->format != CFT_FP128 && c->format != CFT_FP256;
}
static void say_format(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "format %d is not one of CFT_FP64 (%d), CFT_FP128 (%d) or "
                   "CFT_FP256 (%d).", c->format, CFT_FP64, CFT_FP128, CFT_FP256);
}

static int hit_max_iter(const struct cft_support_ctx *c){
    return c->max_iter < 1;
}
static void say_max_iter(const struct cft_support_ctx *c, char *b, size_t n){
    snprintf(b, n, "max_iter = %d. Use 0 for the default at this format, or a "
                   "positive cap; REBOUND uses 12 and binary256 needs about 22.",
             c->max_iter);
}

/* ------------------------------------------------------------------ *
 * The table.
 *
 * Order is the order a doubly-invalid simulation is told about, and
 * nothing else depends on it: every row here stops the run.
 *
 * The `paths` column is the whole of the difference between the two
 * entry points. Five rows are one-sided, and each says why above its
 * predicate: the drop-in supports collisions and attached ODE sets,
 * the subprocess API refuses those and also refuses timestep
 * modifications, ignore-terms, a foreign integrator and a body count
 * past the cap; MEGNO, the ensemble count, the format and max_iter are
 * settings only the drop-in has.
 */
const struct cft_support_row cft_support_rows[] = {
    /* the simulation */
    { "no_particles",      CFT_PATH_SUBPROCESS, hit_no_particles,      say_no_particles      },
    { "body_count",        CFT_PATH_SUBPROCESS, hit_body_count,        say_body_count        },
    { "gravity_custom",    CFT_PATH_BOTH,       hit_gravity_custom,    say_gravity_custom    },
    { "variational",       CFT_PATH_BOTH,       hit_variational,       say_variational       },
    { "megno",             CFT_PATH_DROPIN,     hit_megno,             say_megno             },
    { "additional_forces", CFT_PATH_BOTH,       hit_additional_forces, say_additional_forces },
    { "timestep_mods",     CFT_PATH_SUBPROCESS, hit_timestep_mods,     say_timestep_mods     },
    { "veldep_forces",     CFT_PATH_BOTH,       hit_veldep,            say_veldep            },
    { "gravity_module",    CFT_PATH_BOTH,       hit_gravity_module,    say_gravity_module    },
    { "ignore_terms",      CFT_PATH_SUBPROCESS, hit_ignore_terms,      say_ignore_terms      },
    { "collision",         CFT_PATH_BOTH,       hit_collision,         say_collision         },
    { "ghost_boxes",       CFT_PATH_BOTH,       hit_ghost_boxes,       say_ghost_boxes       },
    { "boundary",          CFT_PATH_BOTH,       hit_boundary,          say_boundary          },
    { "particle_map",      CFT_PATH_BOTH,       hit_particle_map,      say_particle_map      },
    { "test_particles",    CFT_PATH_BOTH,       hit_test_particles,    say_test_particles    },
    { "odes",              CFT_PATH_SUBPROCESS, hit_odes,              say_odes              },
    { "integrator_name",   CFT_PATH_SUBPROCESS, hit_integrator_name,   say_integrator_name   },
    /* the integrator's own settings */
    { "ensemble_E",        CFT_PATH_DROPIN,     hit_ensemble,          say_ensemble          },
    { "adaptive_mode",     CFT_PATH_BOTH,       hit_adaptive_mode,     say_adaptive_mode     },
    { "format",            CFT_PATH_DROPIN,     hit_format,            say_format            },
    { "max_iter",          CFT_PATH_DROPIN,     hit_max_iter,          say_max_iter          },
    { NULL, 0, NULL, NULL }
};

const struct cft_support_row *
cft_support_first_refusal(const struct cft_support_ctx *c, unsigned path){
    if (!c || !c->r) return NULL;
    for (const struct cft_support_row *row = cft_support_rows; row->name; row++){
        if (!(row->paths & path)) continue;
        if (row->hit(c)) return row;
    }
    return NULL;
}
