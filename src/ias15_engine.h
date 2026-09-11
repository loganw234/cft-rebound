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

/* IAS15's adaptive_mode: 2 (PRS23) or 3 (AARSETH85). They share every
 * sum and differ in one expression; 0 and 1 are REBOUND's other
 * branch and are not implemented. */
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
