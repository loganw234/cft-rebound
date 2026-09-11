/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The IAS15 engine of src/ias15_cft.c, exposed to a caller that holds
 * the particles itself. Compile ias15_cft.c with -DIAS15_CFT_LIBRARY to
 * get these entry points and no main().
 *
 * The engine is ONE GLOBAL INSTANCE: every buffer in ias15_cft.c is a
 * file-scope static, deliberately, because the standalone program is one
 * run of one problem. A caller must therefore use it from one place at a
 * time; src/reb_integrator_cft.c enforces that and refuses rather than
 * let two simulations share the vectors.
 *
 * Nothing here rounds except eng_promote/eng_round, which are 754-2019
 * 5.4.2 convertFormat: binary64 into the wide format is exact (binary64
 * is a subset of every wider one) and the wide format back to binary64
 * is the correctly rounded view.
 *
 * The order of calls:
 *
 *     ias15_engine_open(fmt, artifact)
 *     ias15_engine_alloc(n_cap, max_iter, arith_fma, cs_aug, tol_shift, quiet)
 *     ias15_engine_set_epsilon_f64(eps)          once, before the first step
 *     ias15_engine_set_bodies(n)                 n <= n_cap
 *     per step:  set_G, set_masses, put_dt, put_dt_last, [put_xv], step, get_xv
 */
#ifndef IAS15_ENGINE_H
#define IAS15_ENGINE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

/* The engine's own vectors, borrowed. Each is a byte blob of `n_elem`
 * elements of `elem_size` bytes, in the run's format. The names are
 * REBOUND's IAS15 names one for one (integrator_ias15.c), which is what
 * lets the archive descriptors be checked against it line by line. */
struct ias15_engine_view {
    unsigned char *x, *v, *a, *at;
    unsigned char *x0, *v0, *a0;
    unsigned char *csx, *csv, *csa0;
    unsigned char *g[7], *b[7], *e[7], *csb[7], *er[7], *br[7];
    size_t n_elem;      /* 3 * n_bodies */
    size_t elem_size;   /* 8, 16 or 32 */
    size_t n_bodies;
};

/* fmt is a cft_format: CFT_FP64 = 1, CFT_FP128 = 2, CFT_FP256 = 3.
 * artifact NULL selects the software backend.
 * 0 on success; -1 already open, -2 bad format, -3 cft_open failed. */
int    ias15_engine_open(int fmt, const char *artifact);
int    ias15_engine_is_open(void);
int    ias15_engine_format(void);

/* Allocate for up to n_cap bodies. Every later body count must be at or
 * under n_cap: several of the step's scratch vectors are allocated at
 * their first use and never resized, so the capacity is what makes a
 * REBOUND particle addition safe. 0 on success. */
int    ias15_engine_alloc(size_t n_cap, int max_iter, int arith_fma,
                          int cs_augmented, int tol_shift, int quiet);
size_t ias15_engine_capacity(void);

/* The live body count. Rebuilds the pair list. -2 if n > capacity. */
int    ias15_engine_set_bodies(size_t n);

/* epsilon, and the two knobs that are plain switches. Setting epsilon
 * marks sqrt7(epsilon*5040) stale so that the next step derives it
 * again; within one run it is derived once, as the standalone program
 * derives it once. */
int    ias15_engine_set_epsilon_f64(double eps);
int    ias15_engine_adaptive(void);       /* 1 if epsilon > 0 */
void   ias15_engine_set_max_iter(int n);
void   ias15_engine_set_arith_fma(int on);

void   ias15_engine_set_G_f64(double G);

/* ---- which pairs gravity computes ----------------------------------
 *
 * REBOUND's reb_gravity_basic_calculate_acceleration() (gravity.c:167)
 * is two loops and three settings, not one loop over the triangle, and
 * all three change WHICH addends a particle receives and in what order
 * - so they cannot be had by zeroing a mass, and this is where they
 * live. Set them together; the pair-to-particle table is rebuilt once.
 *
 *   n_active           r->N_active. (size_t)-1 - REBOUND's SIZE_MAX
 *                      default - means every particle is active, which
 *                      is what the engine starts at. Clamped to N.
 *   testparticle_type  r->testparticle_type. 0: the particles at and
 *                      above n_active feel the active ones and do not
 *                      pull back. 1: they pull back.
 *   ignore_terms       r->gravity_ignore_terms, by value:
 *                      0 NONE, 1 BETWEEN_0_AND_1, 2 INVOLVING_0.
 *
 * A caller reproducing IAS15 will only ever pass 0 for ignore_terms:
 * reb_integrator_ias15_step() writes NONE over the field at the top of
 * every step (integrator_ias15.c:875), so no value a user sets survives
 * to reach gravity under that integrator. The argument exists so that
 * the value is read from the simulation rather than assumed.
 *
 * Cheap on every step - it returns having touched nothing when the
 * three are unchanged. 0 on success, -1 if not allocated, -2 if
 * ignore_terms is not one of the three. */
int    ias15_engine_set_active(size_t n_active, int testparticle_type,
                               int ignore_terms);

/* One gravity evaluation over the x, masses and restrictions the engine
 * holds now; the answer is where ias15_engine_get_a_f64() reads it.
 * Not part of a step: it is how a gate puts this gravity beside
 * REBOUND's own for one configuration, which is the only way to reach a
 * setting REBOUND's IAS15 overwrites before gravity ever sees it. */
void   ias15_engine_gravity(void);

/* Remove one body from the wide state, shifting every vector down so
 * that each survivor keeps its own polynomial and its own wide tail.
 * This is NOT what REBOUND does - REBOUND shifts its particles and
 * leaves its polynomial arrays alone - so it is only for a caller
 * that has asked for accuracy over equivalence. 0 on success. */
int    ias15_engine_remove_body(size_t index);

/* The opposite choice, and REBOUND's: do not move anything, and read
 * the seven coefficient levels at the new stride. REBOUND keeps them
 * in one flat buffer that dpcast() slices at the current 3N, so a
 * count change below its high-water mark silently re-reads every
 * level at a different offset; each level is its own allocation here,
 * so that has to be performed. Call it BEFORE set_bodies, while the
 * old count is still the engine's. A no-op for an ensemble. */
void   ias15_engine_alias_resize(size_t n_new);

/* ---- the force hook: r->additional_forces at every substage --------
 *
 * REBOUND's IAS15 does not call gravity; it calls
 * reb_simulation_update_acceleration(), which is gravity FOLLOWED BY
 * r->additional_forces(r), at the top of every step attempt and again
 * at each of the seven Gauss-Radau nodes (integrator_ias15.c:461). The
 * engine therefore has to offer the same seam.
 *
 * It is a plain callback over binary64 buffers and not a
 * struct reb_simulation, deliberately: src/ias15_cft.c builds both as
 * the standalone program and, under -DIAS15_CFT_LIBRARY, as a library
 * half that knows nothing about REBOUND, and that separation is load
 * bearing. cft_force_hook() in src/reb_integrator_cft.c is what turns
 * these buffers into particles and back.
 *
 * The engine calls fn with:
 *   n    0 for the call at the top of the attempt, 1..7 for the node
 *   h_n  REBOUND's h[n] as a binary64, +0 at n = 0. The caller wants
 *        this rather than the engine's clock because REBOUND computes
 *        r->t = t_beginning + r->dt*h[n] in binary64 and a
 *        time-dependent force that saw a different t would compute a
 *        different acceleration. The port's derived KH[] rounds to
 *        REBOUND's h literals bit for bit at binary64, so handing the
 *        node's h over is enough for the caller to reproduce the
 *        expression exactly.
 *   dt   the CURRENT attempt's step, rounded to binary64. Not r->dt as
 *        the caller last saw it: a rejected attempt changes it, and the
 *        node times move with it.
 *   x,v  the node's positions and velocities, 3n values each, the
 *        correctly rounded view of the wide state
 *   a    in: the wide gravity, correctly rounded. out: whatever the
 *        routine leaves there.
 *
 * A COMPONENT THE ROUTINE CHANGES BECOMES BINARY64. The engine compares
 * `a` bit for bit against what it handed over and promotes back only
 * where it differs, so an untouched component - and a do-nothing
 * routine - costs nothing and leaves a wide run wide. Where the routine
 * did write, the value it wrote is all there is: r->particles are
 * binary64 and that is the only interface REBOUND offers a callback.
 * Gravity stays wide; the user's force does not. At CFT_FP64, which is
 * where the equivalence claim lives, the two are the same bits and this
 * is REBOUND's own `at[k] = particles[mk].ax` exactly.
 *
 * velocity_dependent is REBOUND's r->force_is_velocity_dependent, and
 * it selects one thing: whether the node's velocities are predicted
 * from the b polynomial before the call. REBOUND predicts them only
 * under `r->calculate_megno || (r->additional_forces &&
 * r->force_is_velocity_dependent)` (integrator_ias15.c:434), and a port
 * that predicted them always would hand a velocity-independent force a
 * different v from the one REBOUND hands it - not a difference the
 * force uses, but not one to leave to luck either.
 *
 * fn NULL removes the hook, which is the no-force path and pays for
 * nothing. */
typedef void (*ias15_engine_force_fn)(void *user, int n, double h_n, double dt,
                                      const double *x, const double *v, double *a);
void   ias15_engine_set_force_hook(ias15_engine_force_fn fn, void *user,
                                   int velocity_dependent);

/* The shadow ias15_engine_alias_resize() performs that re-reading
 * through: six flat arrays, one per seven-level coefficient family,
 * holding exactly what REBOUND's own six allocations hold - including
 * the tail above the live region that a shrink strands and that a
 * regrow under the high-water mark reads straight back.
 *
 * Borrowed, like ias15_engine_view()'s vectors, and for the same
 * caller: src/cft_archive.c has to put this in a Simulationarchive and
 * take it back out, because within one process the shadow is what makes
 * a remove-then-regrow bit-identical to REBOUND and across a checkpoint
 * nothing carried it (ROADMAP.md, "What the collision work leaves open,
 * across a checkpoint").
 *
 * The shadow does not exist until a body count change needs it, so
 * `create` says what to do about that: 0 reports what is there, and
 * `flat[0]` is NULL when nothing is - which is how a caller asks "has
 * this run ever resized". 1 allocates it, zeroed, as alias_resize()
 * would. n_elem is its full length, 7*3*capacity; a caller writing
 * fewer elements than that must zero the rest itself. */
struct ias15_engine_alias_view {
    unsigned char *flat[6];   /* g, b, e, csb, er, br - alias_resize()'s order */
    size_t n_elem;            /* elements in EACH: 7 * 3 * ias15_engine_capacity() */
    size_t elem_size;         /* 8, 16 or 32 */
};
void   ias15_engine_alias_view(struct ias15_engine_alias_view *out, int create);

/* IAS15's adaptive_mode, all four of REBOUND's: 0 INDIVIDUAL,
 * 1 GLOBAL, 2 PRS23 (REBOUND's default since January 2024) and
 * 3 AARSETH85. 2 and 3 share every sum and differ in one
 * expression; 0 and 1 take REBOUND's other error estimate, and 0
 * also changes the CORRECTOR's convergence test - see pc_error()
 * in src/ias15_cft.c, which is the half the step control does not
 * reach. */
void   ias15_engine_set_adaptive_mode(int mode);
/* IAS15's min_dt, a floor on |dt|. 0 disables it, which is REBOUND's
 * default: at 0 the comparison is false and nothing is selected. */
void   ias15_engine_set_min_dt_f64(double min_dt);
/* r->softening. Squared internally in the run's format; 0 is the
 * default and issues the same +0 addend the port always issued. */
void   ias15_engine_set_softening_f64(double softening);
void   ias15_engine_set_masses_f64(const double *m);      /* n_bodies values */

void   ias15_engine_put_xv_f64(const double *x, const double *v);  /* 3n each */
void   ias15_engine_get_xv_f64(double *x, double *v);
void   ias15_engine_get_a_f64(double *a);

void   ias15_engine_put_dt_f64(double dt);
void   ias15_engine_put_dt_last_f64(double dt_last);
double ias15_engine_get_dt_f64(void);
double ias15_engine_get_dt_last_f64(void);
void   ias15_engine_put_t_f64(double t);
double ias15_engine_get_t_f64(void);

/* One accepted IAS15 step: step_attempt() repeated until the step is
 * accepted, which is what reb_integrator_ias15_step does with
 * reb_integrator_ias15_step_try. Returns the number of attempts. */
long   ias15_engine_step(double *dt_done_out);

/* Zero b, e, g, er, br, csb, csx, csv and the step counters: the
 * polynomial describes a particle set that no longer exists. */
void   ias15_engine_reset_state(void);

void   ias15_engine_view(struct ias15_engine_view *out);

uint32_t           ias15_engine_flags(void);       /* the OR of every operation's flags */
unsigned long long ias15_engine_calls(void);
unsigned long long ias15_engine_max_exceeded(void);
uint64_t           ias15_engine_constants_digest(void);

#ifdef __cplusplus
}
#endif
#endif /* IAS15_ENGINE_H */
