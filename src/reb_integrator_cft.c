/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The REBOUND integrator shim: struct reb_integrator over the IAS15
 * engine of src/ias15_cft.c.
 *
 * One step is
 *
 *   1. promote r->particles into the wide format - exact, because
 *      binary64 is a subset of every wider one - but only where they
 *      changed under us, so that a wide run is not truncated back to
 *      binary64 once a step;
 *   2. one IAS15 step through the engine, which is REBOUND's
 *      reb_integrator_ias15_step_try repeated until a step is accepted;
 *   3. round the wide result back into r->particles and advance r->t,
 *      r->dt and r->dt_last_done the way REBOUND does.
 *
 * Nothing here performs an arithmetic operation on a coordinate. The
 * only floating-point work is 754-2019 5.4.2 convertFormat, inside the
 * engine; the comparisons below are memcmp on bit patterns, which is
 * why a signalling value or a -0 cannot be misread as "unchanged".
 *
 * THE ENGINE IS ONE GLOBAL INSTANCE. Every buffer in ias15_cft.c is a
 * file-scope static. Two simulations therefore cannot use this
 * integrator at the same time, and a second one is refused with a
 * message that says so. A simulation that has been freed releases the
 * engine and the next one adopts it, at the same format - the Gauss-
 * Radau constants are derived at that format and every buffer is sized
 * for its element width - and within the reserved capacity. epsilon,
 * max_iter and the arithmetic form are reconciled on every step, so a
 * program may run one simulation after another at different tolerances,
 * and a user who edits the state mid-run gets what they asked for.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "rebound.h"
#include "cft.h"
#include "cft_ias15.h"
#include "cft_supported.h"
#include "ias15_engine.h"

/* ------------------------------------------------------------------ */
/* Refusals                                                            */
/* ------------------------------------------------------------------ */
/* A refusal is a REBOUND error message plus REB_STATUS_GENERIC_ERROR,
 * which stops reb_simulation_steps() and reb_simulation_integrate()
 * at the top of their loops. The user's process is not killed and the
 * simulation is left exactly as it was: nothing was computed. */
static void refuse(struct reb_simulation *r, const char *fmt, ...){
    /* Wide enough for the longest row in src/cft_supported.c plus the
     * "ias15_cft: " prefix. The collision row alone is over 400. */
    char buf[1024];
    va_list ap; va_start(ap, fmt);
    vsnprintf(buf, sizeof buf, fmt, ap);
    va_end(ap);
    reb_simulation_error(r, buf);
    if (r) r->status = REB_STATUS_GENERIC_ERROR;
}

/* ------------------------------------------------------------------ */
/* Process-wide settings and the binary64 view                         */
/* ------------------------------------------------------------------ */
static size_t      reserve_N;          /* cft_ias15_reserve(); 0 = the N of the first step */
static const char *artifact_path;

/* The artifact this process will open, or NULL for the software
 * backend. An explicit cft_ias15_set_artifact() wins; otherwise
 * CFT_REBOUND_ARTIFACT names one, which is what lets the gate suite
 * run against a card without a single gate being modified. Empty is
 * treated as unset, so CFT_REBOUND_ARTIFACT= forces software. */
static const char *resolve_artifact(void){
    const char *p;
    if (artifact_path && artifact_path[0]) return artifact_path;
    p = getenv("CFT_REBOUND_ARTIFACT");
    if (p && p[0]) return p;
    return NULL;
}
static int         pc_tol_shift = -1;

static struct cft_ias15_state *engine_owner;   /* the state the engine is bound to */
static int    engine_ready;                    /* the engine has been allocated */
/* REBOUND's ias15->N_allocated, in particles rather than in 3N.
 * reb_integrator_ias15_alloc() zeroes the polynomial only when 3N
 * exceeds it, so a removal and a regrow under the mark leave the
 * polynomial alone. Reproduced here because the equivalence claim is
 * to REBOUND's behaviour, not to the behaviour REBOUND ought to
 * have. */
static size_t engine_hiwater_N;
static int    engine_fmt;
static double engine_eps;
static int    engine_max_iter;
static int    engine_arith_fma;

/* the binary64 view we last wrote into r->particles, r->t and r->dt.
 * If REBOUND or the user has changed one of them since, that value is
 * theirs and is promoted back in; if not, the wide state is the truth
 * and is left alone. At CFT_FP64 both branches are the same bits. */
static double *view_x, *view_v;
static double *tmp_x, *tmp_v, *tmp_a, *tmp_m;   /* the step's scratch, sized once */
static size_t  view_N_allocated;
static size_t  view_N;
static double  view_t, view_dt, view_dt_last;
static int     view_valid;

static int same_bits(const double *a, const double *b, size_t n){
    return memcmp(a, b, n * sizeof(double)) == 0;
}

static void view_alloc(size_t n){
    if (n <= view_N_allocated && view_x) return;
    free(view_x); free(view_v); free(tmp_x); free(tmp_v); free(tmp_a); free(tmp_m);
    view_x = malloc(3 * n * sizeof(double));
    view_v = malloc(3 * n * sizeof(double));
    tmp_x  = malloc(3 * n * sizeof(double));
    tmp_v  = malloc(3 * n * sizeof(double));
    tmp_a  = malloc(3 * n * sizeof(double));
    tmp_m  = malloc(n * sizeof(double));
    if (!view_x || !view_v || !tmp_x || !tmp_v || !tmp_a || !tmp_m){
        fprintf(stderr, "ias15_cft: out of memory\n"); exit(1); }
    view_N_allocated = n;
}

/* ------------------------------------------------------------------ */
/* What this integrator does not do                                    */
/* ------------------------------------------------------------------ */
/* The list is src/cft_supported.h - one table both entry points walk,
 * so that removing a refusal is one row rather than the same edit in
 * four files. ROADMAP.md's last section is the record of what the other
 * arrangement costs.
 *
 * What stays here is this path's half: the state, the one fixup that is
 * a write rather than a question, and the "ias15_cft: " prefix, which
 * belongs to the integrator and not to the fact. */
static int supported(struct reb_simulation *r, struct cft_ias15_state *st){
    if (!st){ refuse(r, "ias15_cft: no integrator state (was "
                        "reb_simulation_set_integrator called?)"); return 0; }

    /* Not a question, so not a row: 0 means "the default for this
     * format", which create() cannot resolve because the caller chooses
     * the format afterwards. Written back, so the state - and any
     * archive of it - carries the number actually used rather than the
     * sentinel. Guarded on a valid format because it used to sit after
     * the format check and an invalid format must still be refused by
     * name rather than by max_iter. */
    if (st->max_iter == 0 &&
        (st->format == CFT_FP64 || st->format == CFT_FP128 || st->format == CFT_FP256))
        st->max_iter = cft_ias15_default_max_iter(st->format);

    struct cft_support_ctx c;
    memset(&c, 0, sizeof c);
    c.r             = r;
    c.have_state    = 1;
    c.format        = st->format;
    c.accurate      = st->accurate;
    c.max_iter      = st->max_iter;
    c.E             = st->E;
    c.adaptive_mode = st->adaptive_mode;

    const struct cft_support_row *row =
        cft_support_first_refusal(&c, CFT_PATH_DROPIN);
    if (row){
        char why[768];
        row->say(&c, why, sizeof why);
        refuse(r, "ias15_cft: %s", why);
        return 0;
    }
    return 1;
}

/* ------------------------------------------------------------------ */
/* Binding the one engine to this simulation                           */
/* ------------------------------------------------------------------ */
static void publish_state(struct cft_ias15_state *st){
    struct ias15_engine_view w;
    ias15_engine_view(&w);
    st->x = w.x; st->v = w.v;
    st->x0 = w.x0; st->v0 = w.v0; st->a0 = w.a0;
    st->csx = w.csx; st->csv = w.csv; st->csa0 = w.csa0;
    for (int m = 0; m < 7; m++){
        st->g[m] = w.g[m]; st->b[m] = w.b[m]; st->e[m] = w.e[m];
        st->br[m] = w.br[m]; st->er[m] = w.er[m]; st->csb[m] = w.csb[m];
    }
    st->n_elem = w.n_elem;
    st->E = 1;
    st->constants_digest = ias15_engine_constants_digest();
}

/* A state that came out of an archive owns its blobs: REBOUND's loader
 * realloc'd each one and read the file into it. publish_state would
 * overwrite those pointers with the engine's, dropping the whole
 * restored state on the first step and continuing from an engine
 * reset_state() had just zeroed - which is what tests/gate_real.c
 * measured before this existed. Copy them in, release them, and let
 * publish_state repoint the state at the engine like any other.
 *
 * "Came out of an archive" needs no flag: publish_state only ever sets
 * these to the engine's own buffers, and bind_engine refuses a second
 * simulation while the first owns the engine, so a non-NULL pointer
 * that is not the engine's can only have been loaded.
 *
 * Returns 0 having refused; the caller must not go on to step. */
static int adopt_loaded_state(struct reb_simulation *r, struct cft_ias15_state *st){
    struct ias15_engine_view w;
    ias15_engine_view(&w);
    if (!st->x || st->x == w.x) return 1;          /* nothing was loaded */
    /* The negative control for tools/check_checkpoint.py. With this set,
     * the loaded state is deliberately NOT carried into the engine - the
     * exact defect gate_real found - so the checker can require its own
     * comparison to FAIL and thereby prove it would notice. A gate that
     * cannot fail on demand is not evidence. */
    if (getenv("CFT_REBOUND_NO_ADOPT")){
        /* Discard rather than leak: publish_state is about to overwrite
         * these pointers either way, and nothing else frees them. */
        free(st->x);   free(st->v);
        free(st->x0);  free(st->v0);  free(st->a0);
        free(st->csx); free(st->csv); free(st->csa0);
        for (int m = 0; m < 7; m++){
            free(st->g[m]);  free(st->b[m]);  free(st->e[m]);
            free(st->br[m]); free(st->er[m]); free(st->csb[m]);
        }
        return 1;
    }

    size_t width = cft_format_size((cft_format)st->format);
    if (!width) return 1;
    if (st->n_elem != w.n_elem){
        refuse(r, "ias15_cft: the archive holds %zu wide elements and the engine "
                  "is allocated for %zu (%zu bodies). Load into a simulation with "
                  "the same particle count, or call cft_ias15_reserve() first.",
               st->n_elem, w.n_elem, w.n_bodies);
        return 0;
    }
    const size_t bytes = st->n_elem * width;

#define CFT_ADOPT(dst, src) do{ \
        if ((src) && (src) != (dst)){ memcpy((dst), (src), bytes); free(src); (src) = NULL; } \
    }while(0)
    CFT_ADOPT(w.x,    st->x);    CFT_ADOPT(w.v,    st->v);
    CFT_ADOPT(w.x0,   st->x0);   CFT_ADOPT(w.v0,   st->v0);   CFT_ADOPT(w.a0,   st->a0);
    CFT_ADOPT(w.csx,  st->csx);  CFT_ADOPT(w.csv,  st->csv);  CFT_ADOPT(w.csa0, st->csa0);
    for (int m = 0; m < 7; m++){
        CFT_ADOPT(w.g[m],  st->g[m]);  CFT_ADOPT(w.b[m],  st->b[m]);
        CFT_ADOPT(w.e[m],  st->e[m]);  CFT_ADOPT(w.br[m], st->br[m]);
        CFT_ADOPT(w.er[m], st->er[m]); CFT_ADOPT(w.csb[m], st->csb[m]);
    }
#undef CFT_ADOPT

    /* The clock is REBOUND's own, restored from its own fields, and the
     * engine's copy was not reset with the rest - so write it here
     * rather than leave the step's "has it changed" guard to notice. */
    ias15_engine_put_t_f64(r->t);
    ias15_engine_put_dt_f64(r->dt);
    ias15_engine_put_dt_last_f64(r->dt_last_done);
    view_t = r->t; view_dt = r->dt; view_dt_last = r->dt_last_done;

    /* The wide x and v ARE the truth now, so the step must not
     * re-promote r->particles over them: that would round a binary256
     * restart back to binary64 in its first step. Publishing the view
     * the engine now holds is what makes the guard agree. */
    ias15_engine_get_xv_f64(view_x, view_v);
    view_valid = 1;
    return 1;
}

/* epsilon, max_iter and the arithmetic form may change between steps -
 * the format may not. Applied on every step, so a user who edits the
 * state gets what they asked for rather than what the first step saw. */
static int sync_config(struct reb_simulation *r, struct cft_ias15_state *st){
    if (st->format != engine_fmt){
        refuse(r, "ias15_cft: this process opened the engine at %s and the format cannot "
                  "be changed afterwards (the Gauss-Radau constants are derived at that "
                  "format and every buffer is sized for its element width). Requested %s.",
               cft_format_name((cft_format)engine_fmt), cft_format_name((cft_format)st->format));
        return 0;
    }
    if (memcmp(&st->epsilon, &engine_eps, sizeof(double)) != 0){
        if (ias15_engine_set_epsilon_f64(st->epsilon) != 0){
            refuse(r, "ias15_cft: epsilon could not be set."); return 0; }
        engine_eps = st->epsilon;
    }
    if (st->max_iter != engine_max_iter){ ias15_engine_set_max_iter(st->max_iter); engine_max_iter = st->max_iter; }
    if (st->arith_fma != engine_arith_fma){ ias15_engine_set_arith_fma(st->arith_fma); engine_arith_fma = st->arith_fma; }
    return 1;
}

static int bind_engine(struct reb_simulation *r, struct cft_ias15_state *st){
    if (engine_owner == st) return sync_config(r, st);
    if (engine_owner){
        refuse(r, "ias15_cft: another simulation is already using this integrator. The "
                  "engine in src/ias15_cft.c is one global instance (every buffer in it "
                  "is a file-scope static), so one simulation may use ias15_cft at a "
                  "time; free the other simulation first.");
        return 0;
    }
    if (!engine_ready){
        size_t cap = reserve_N > r->N ? reserve_N : r->N;
        const char *art = resolve_artifact();
        if (ias15_engine_open(st->format, art) != 0){
            refuse(r, "ias15_cft: cft_open(%s) failed: %s",
                   art ? art : "software backend", cft_last_error());
            return 0;
        }
        int shift = pc_tol_shift;
        if (shift < 0){
            shift = (st->format == CFT_FP128) ? 113 - 53
                  : (st->format == CFT_FP256) ? 237 - 53 : 0;
        }
        if (ias15_engine_alloc(cap, st->max_iter, st->arith_fma, 0, shift, 1) != 0){
            refuse(r, "ias15_cft: the engine could not be allocated for %zu bodies.", cap);
            return 0;
        }
        if (ias15_engine_set_epsilon_f64(st->epsilon) != 0){
            refuse(r, "ias15_cft: epsilon could not be set."); return 0;
        }
        engine_ready = 1;
        engine_fmt = st->format; engine_eps = st->epsilon;
        engine_max_iter = st->max_iter; engine_arith_fma = st->arith_fma;
    }else{
        /* The engine outlives the simulation that opened it; a later one
         * adopts it. sync_config decides what may still change. */
        if (!sync_config(r, st)) return 0;
    }
    if (r->N > ias15_engine_capacity()){
        refuse(r, "ias15_cft: %zu particles, but the engine was allocated for %zu. "
                  "Several of the step's scratch vectors are sized at their first use, "
                  "so call cft_ias15_reserve(N) before the first step.",
               r->N, ias15_engine_capacity());
        return 0;
    }
    if (ias15_engine_set_bodies(r->N) != 0){
        refuse(r, "ias15_cft: the engine refused %zu bodies.", r->N); return 0; }
    /* The first bind is REBOUND's first alloc: N_allocated is 0, so
     * 3N always exceeds it and the polynomial starts from zero. */
    ias15_engine_reset_state();
    engine_hiwater_N = r->N;
    view_alloc(ias15_engine_capacity());
    view_N = r->N;
    view_valid = 0;
    engine_owner = st;
    if (!adopt_loaded_state(r, st)){ engine_owner = NULL; return 0; }
    publish_state(st);
    return 1;
}

/* ------------------------------------------------------------------ */
/* r->additional_forces, at every substage                             */
/* ------------------------------------------------------------------ */
/* REBOUND's IAS15 calls reb_simulation_update_acceleration(), not
 * gravity: simulation.c:643 is gravity followed by
 * r->additional_forces(r), and integrator_ias15.c calls it at :290 for
 * the step-start acceleration and again at :461 for each of the seven
 * Gauss-Radau nodes, reading particles[mk].ax straight back into at[]
 * at :467. The engine offers that seam as a callback over binary64
 * buffers (src/ias15_engine.h) because it may not include rebound.h;
 * this is the half that knows what a particle is.
 *
 * It lives here rather than in a src/cft_forces.c of its own for one
 * reason: a new object in src/ would have to be added to DROPIN_OBJ,
 * to PIC_OBJ and to the install rules, which is a larger edit to the
 * Makefile than the whole of this is to this file.
 *
 * THE USER'S FORCE IS EVALUATED AT BINARY64, even in a binary256 run.
 * r->particles are binary64 and they are the only interface REBOUND
 * offers a callback, so there is nowhere wider to hand the routine its
 * inputs or to take its answer back. Gravity stays wide; a component
 * the routine writes is binary64 from that node on. That is a real
 * ceiling of this feature and not a detail: docs and README belong to
 * the integrator, so it is stated here and in src/ias15_engine.h. */
static double force_t_beginning;   /* REBOUND's t_beginning for this step */

static void cft_force_hook(void *user, int n, double h_n, double dt,
                           const double *x, const double *v, double *a){
    struct reb_simulation *r = user;
    const size_t N = r->N;

    /* r->t at the node, as REBOUND's own binary64 expression:
     *     r->t = t_beginning + r->dt * h[n];      (integrator_ias15.c:408)
     * NOT the engine's wide clock rounded - the two differ in the last
     * bits, and a time-dependent force that saw a different time would
     * compute a different acceleration. h_n is the port's derived KH[n]
     * rounded to binary64, which tools/gen_constants.py checks is
     * REBOUND's h[] literal bit for bit. `dt` is the attempt's step and
     * not r->dt as it stood when this step began, because a rejected
     * attempt changes it and the node times move with it.
     *
     * n = 0 is the call at the top of the attempt, where REBOUND has
     * not touched r->t and it still holds t_beginning - which is also
     * what integrator_ias15.c:608 restores it to after the loop, so
     * writing it here is that restore for the next attempt. */
    const double dth = dt * h_n;
    r->dt = dt;
    r->t = n ? force_t_beginning + dth : force_t_beginning;

    for (size_t i = 0; i < N; i++){
        r->particles[i].x  = x[3*i]; r->particles[i].y  = x[3*i+1]; r->particles[i].z  = x[3*i+2];
        r->particles[i].vx = v[3*i]; r->particles[i].vy = v[3*i+1]; r->particles[i].vz = v[3*i+2];
        r->particles[i].ax = a[3*i]; r->particles[i].ay = a[3*i+1]; r->particles[i].az = a[3*i+2];
    }
    r->additional_forces(r);
    for (size_t i = 0; i < N; i++){
        a[3*i] = r->particles[i].ax; a[3*i+1] = r->particles[i].ay; a[3*i+2] = r->particles[i].az;
    }
    /* The particles are deliberately left at the node, positions and
     * all: REBOUND leaves them there too until the end of step_try, and
     * the step below writes the accepted view over them afterwards. */
}

/* ------------------------------------------------------------------ */
/* The step                                                            */
/* ------------------------------------------------------------------ */
static void cft_ias15_step(struct reb_simulation *r, void *p){
    struct cft_ias15_state *st = p;

    r->gravity_ignore_terms = REB_GRAVITY_IGNORE_TERMS_NONE;   /* as REBOUND's ias15 step does */

    if (!supported(r, st)) return;

    if (r->N == 0){                 /* REBOUND's own empty branch, verbatim */
        r->t += r->dt;
        r->dt_last_done = r->dt;
        return;
    }
    if (!bind_engine(r, st)) return;

    /* the body count changed under us: did_add_particle /
     * will_remove_particle set view_valid to 0 and the state was
     * invalidated there. Re-size and start the polynomial again. */
    if (view_N != r->N){
        if (r->N > ias15_engine_capacity()){
            refuse(r, "ias15_cft: the simulation grew to %zu particles and the engine was "
                      "allocated for %zu. Call cft_ias15_reserve(N) before the first step.",
                   r->N, ias15_engine_capacity());
            return;
        }
        /* REBOUND does not move its coefficient levels when the count
         * changes; it re-reads them at the new stride, because all
         * seven share one buffer that dpcast() slices at the current
         * N3. Here they are seven allocations, so that reading is
         * performed - before set_bodies, while the old count is still
         * the engine's. At accurate = 1 the state was already shifted
         * in will_remove_particle and must not be re-read. */
        if (!st->accurate) ias15_engine_alias_resize(r->N);
        if (ias15_engine_set_bodies(r->N) != 0){
            refuse(r, "ias15_cft: the engine refused %zu bodies.", r->N); return; }
        /* And zero only past the high-water mark, which is what
         * REBOUND's `if (N3 > N_allocated)` amounts to: realloc_dp7
         * zeroes the whole array, and nothing below the mark does. */
        if (r->N > engine_hiwater_N){
            ias15_engine_reset_state();
            engine_hiwater_N = r->N;
        }
        view_N = r->N;
        view_valid = 0;
        publish_state(st);
    }

    const size_t N = r->N, N3 = 3 * N;

    /* --- masses and G, every step: cheap, and a user may change them - */
    for (size_t i = 0; i < N; i++) tmp_m[i] = r->particles[i].m;
    ias15_engine_set_masses_f64(tmp_m);
    ias15_engine_set_G_f64(r->G);
    /* Set every step for the same reason G and the masses are: a user
     * may change it between steps and should get what they asked for
     * rather than what the first step saw. */
    ias15_engine_set_softening_f64(r->softening);
    ias15_engine_set_min_dt_f64(st->min_dt);
    ias15_engine_set_adaptive_mode(st->adaptive_mode);
    /* And the three settings of the SUM, read where REBOUND's gravity
     * reads them - after the line at the top of this function that
     * writes NONE over gravity_ignore_terms, because that line is
     * reb_integrator_ias15_step() (integrator_ias15.c:875) and gravity
     * runs downstream of it. So r->gravity_ignore_terms is always NONE
     * by the time it is passed here, which is exactly what REBOUND's
     * IAS15 computes; it is read rather than hard-coded so that the
     * port follows the field rather than a claim about it. */
    ias15_engine_set_active(r->N_active, r->testparticle_type,
                            (int)r->gravity_ignore_terms);

    /* --- the particles ------------------------------------------------
     * The wide state is the truth. r->particles are re-promoted only if
     * they differ, bit for bit, from the view this integrator last wrote
     * into them - which happens on the first step, after a particle was
     * added or removed, and whenever the user or a
     * pre_timestep_modifications callback edited a coordinate. At
     * CFT_FP64 the two paths are the same bits, so the equivalence gate
     * cannot tell them apart; at CFT_FP128 and above, re-reading a view
     * that nobody touched would truncate the run to binary64 once a
     * step, which is exactly the mistake this guard exists to avoid. */
    for (size_t i = 0; i < N; i++){
        tmp_x[3*i] = r->particles[i].x; tmp_x[3*i+1] = r->particles[i].y; tmp_x[3*i+2] = r->particles[i].z;
        tmp_v[3*i] = r->particles[i].vx; tmp_v[3*i+1] = r->particles[i].vy; tmp_v[3*i+2] = r->particles[i].vz;
    }
    if (!view_valid || !same_bits(tmp_x, view_x, N3) || !same_bits(tmp_v, view_v, N3))
        ias15_engine_put_xv_f64(tmp_x, tmp_v);

    /* --- the clock and the step ---------------------------------------
     * r->dt is REBOUND's, and the driver rewrites it to land exactly on
     * tmax when exact_finish_time is 1; r->dt_last_done is reset to 0 at
     * the top of every reb_simulation_integrate(), which is what tells
     * step_try not to predict e and b after a first-attempt rejection.
     * Both are adopted whenever they differ from what we last wrote. */
    if (!view_valid || memcmp(&r->t, &view_t, sizeof(double)) != 0)
        ias15_engine_put_t_f64(r->t);
    if (!view_valid || memcmp(&r->dt, &view_dt, sizeof(double)) != 0)
        ias15_engine_put_dt_f64(r->dt);
    if (!view_valid || memcmp(&r->dt_last_done, &view_dt_last, sizeof(double)) != 0)
        ias15_engine_put_dt_last_f64(r->dt_last_done);

    /* --- the force hook ------------------------------------------------
     * Registered every step, for the same reason G and the masses are
     * set every step: a user may attach or detach a routine between
     * steps, or flip force_is_velocity_dependent, and should get what
     * they asked for rather than what the first step saw. A NULL
     * routine removes the hook, so the no-force path costs one test of
     * a null pointer per Gauss-Radau node and nothing else. */
    force_t_beginning = r->t;
    ias15_engine_set_force_hook(r->additional_forces ? cft_force_hook : NULL, r,
                                r->force_is_velocity_dependent != 0);

    /* --- one accepted step -------------------------------------------- */
    double dt_done = 0;
    ias15_engine_step(&dt_done);

    /* --- the binary64 view -------------------------------------------- */
    ias15_engine_get_xv_f64(view_x, view_v);
    ias15_engine_get_a_f64(tmp_a);
    for (size_t i = 0; i < N; i++){
        r->particles[i].x = view_x[3*i]; r->particles[i].y = view_x[3*i+1]; r->particles[i].z = view_x[3*i+2];
        r->particles[i].vx = view_v[3*i]; r->particles[i].vy = view_v[3*i+1]; r->particles[i].vz = view_v[3*i+2];
        r->particles[i].ax = tmp_a[3*i]; r->particles[i].ay = tmp_a[3*i+1]; r->particles[i].az = tmp_a[3*i+2];
    }
    r->t = ias15_engine_get_t_f64();
    r->dt = ias15_engine_get_dt_f64();
    r->dt_last_done = dt_done;
    view_t = r->t; view_dt = r->dt; view_dt_last = r->dt_last_done;
    view_valid = 1;
    publish_state(st);
}

/* ------------------------------------------------------------------ */
/* create / free / the particle hooks                                  */
/* ------------------------------------------------------------------ */
/* REBOUND's create() takes no arguments - it cannot see r, and so it
 * cannot see r->N. Everything that depends on the body count is
 * therefore allocated at the first step, which is the first moment the
 * particle count is knowable. */
static void *cft_ias15_create(void){
    struct cft_ias15_state *st = calloc(1, sizeof *st);
    if (!st) return NULL;
    st->epsilon = 1e-9;              /* REBOUND's default */
    st->min_dt = 0.0;
    st->adaptive_mode = 2;           /* PRS23 */
    st->format = CFT_FP64;
    st->max_iter = 0;                /* 0 = the default for the format the
                                      * caller is about to choose; resolved
                                      * on the first step */
    st->arith_fma = 0;               /* REBOUND's sequence of roundings */
    st->E = 1;
    snprintf(st->cft_abi, sizeof st->cft_abi, "%u.%u",
             (unsigned)(cft_abi_version() >> 16), (unsigned)(cft_abi_version() & 0xffff));
    return st;
}

static void cft_ias15_free(void *p){
    struct cft_ias15_state *st = p;
    if (!st) return;
    if (engine_owner == st){
        /* The engine's buffers are the process's and stay allocated;
         * releasing ownership is what lets the next simulation adopt
         * them. The state's blobs point into them and must not be
         * freed here. */
        engine_owner = NULL;
        view_valid = 0;
        view_N = 0;
        /* The hook holds the simulation that is being freed. Nothing
         * would call it - every step re-registers before it steps - but
         * a dangling pointer in a process-global is not something to
         * leave lying around for the next simulation to adopt. */
        if (ias15_engine_is_open()) ias15_engine_set_force_hook(NULL, NULL, 0);
    }else{
        /* Never bound, so any blob it holds is one REBOUND's loader
         * allocated - see adopt_loaded_state for why that is the only
         * way a non-owner has one. Freeing them here is what keeps a
         * probe-and-discard load from leaking the whole wide state. */
        free(st->x);   free(st->v);
        free(st->x0);  free(st->v0);  free(st->a0);
        free(st->csx); free(st->csv); free(st->csa0);
        for (int m = 0; m < 7; m++){
            free(st->g[m]);  free(st->b[m]);  free(st->e[m]);
            free(st->br[m]); free(st->er[m]); free(st->csb[m]);
        }
    }
    free(st);
}

/* Adding or removing a particle invalidates the view: the next step
 * resizes the engine and re-promotes every coordinate from
 * r->particles.
 *
 * What happens to the polynomial is REBOUND's, in both directions.
 * Past the high-water mark reb_integrator_ias15_alloc() reallocates
 * and realloc_dp7 zeroes the whole array, which the step reproduces
 * with reset_state. Below the mark REBOUND reallocates nothing and
 * re-reads its seven levels at the new stride, which the step
 * reproduces with ias15_engine_alias_resize. So the binary64 gate's
 * add-a-particle case and its merge case are both bit-identical, and
 * for the same reason: the behaviour is performed, not guessed at.
 *
 * At a wide format accurate = 0 still costs the survivors their wide
 * coordinates, because after the shift the only thing describing them
 * is the binary64 view. Nothing is lost at CFT_FP64, where the view is
 * the state. Above it, set state->accurate = 1, which shifts the wide
 * state itself - a removal only; an add there still re-promotes. */
static void cft_ias15_did_add_particle(struct reb_simulation *r){
    (void)r;
    view_valid = 0;
}

static void cft_ias15_will_remove_particle(struct reb_simulation *r, size_t index){
    /* REBOUND calls this BEFORE its own range check, so an index past
     * the end arrives here and is then refused without a particle ever
     * being removed (particle.c: the hook, then `if (index >= r->N)`).
     * Both built-in integrators that implement this hook guard against
     * it and so does this one: acting on it would shift the state for a
     * removal that never happened. */
    if (!r || index >= r->N) return;

    struct cft_ias15_state *st = cft_ias15_get_state(r);
    if (st && st->accurate && engine_owner == st && engine_ready){
        /* Shift the wide state with the particles, so each survivor
         * keeps its own polynomial. Not what REBOUND does. */
        if (ias15_engine_remove_body(index) == 0){
            /* The view has to move with it, or the step compares
             * REBOUND's shifted particles against an unshifted view,
             * decides they changed, and re-promotes the whole vector -
             * which is the truncation this mode exists to avoid. */
            ias15_engine_get_xv_f64(view_x, view_v);
            return;
        }
    }
    view_valid = 0;
}

const struct reb_integrator cft_ias15_integrator = {
    .documentation =
        "IAS15 with every floating-point operation issued through libcft, the "
        "IEEE 754-2019 binary32/64/128/256 library of the cft-fp256 project, so "
        "that the same integrator runs at binary64, binary128 and binary256 with "
        "the same bits everywhere.\n\n"
        "At binary64 it is REBOUND's own IAS15 bit for bit: the sequence of "
        "roundings is integrator_ias15.c's, operation by operation. The wide "
        "state lives in struct cft_ias15_state; r->particles are the binary64 "
        "view of it, correctly rounded after every step, so every REBOUND "
        "output, callback and visualisation works unchanged.\n\n"
        "Set state->format to CFT_FP128 or CFT_FP256 for a wide run. Collisions, "
        "ghost boxes, boundaries, gravity modules other than REB_GRAVITY_BASIC, "
        "non-zero softening, test particles, r->map and variational particles are "
        "refused rather than approximated.\n\n"
        "r->additional_forces is called where REBOUND's IAS15 calls it: after "
        "gravity at the top of every step attempt and at each of the seven "
        "Gauss-Radau nodes, with r->t set to t_beginning + r->dt*h[n] and, when "
        "r->force_is_velocity_dependent is set, with the node's predicted "
        "velocities. The routine is evaluated at binary64 even in a wide run, "
        "because r->particles are binary64 and are the only interface REBOUND "
        "offers a callback: gravity stays wide, the user's force does not.",
    .step = cft_ias15_step,
    .synchronize = NULL,       /* IAS15 is not a DKD scheme; REBOUND's own is NULL too */
    .create = cft_ias15_create,
    .free = cft_ias15_free,
    .did_add_particle = cft_ias15_did_add_particle,
    .will_remove_particle = cft_ias15_will_remove_particle,
    .field_descriptor_list = cft_ias15_field_descriptor_list,
};

/* ------------------------------------------------------------------ */
/* The registration and the process-wide settings                      */
/* ------------------------------------------------------------------ */
static char registered_name[64];

void cft_ias15_register(const char *name){
    if (!name) name = "ias15_cft";
    if (registered_name[0] && strcmp(registered_name, name) == 0) return;   /* already done */
    reb_integrator_register(cft_ias15_integrator, name);
    if (!registered_name[0]) snprintf(registered_name, sizeof registered_name, "%s", name);
}

struct cft_ias15_state *cft_ias15_get_state(struct reb_simulation *r){
    if (!r || !r->integrator.name) return NULL;
    if (r->integrator.callbacks.step != cft_ias15_step) return NULL;
    return r->integrator.state;
}

int cft_ias15_default_max_iter(int format){
    /* REBOUND's 12 is right at binary64 and is part of what makes the
     * binary64 run REBOUND's own, bit for bit. Above it the measured
     * pass counts are 4-9 at binary128 and 10-20 at binary256
     * (docs/VALIDATION.md), so these sit clear of them. A converged
     * step exits early, so a cap costs nothing until it is reached. */
    switch (format){
        case CFT_FP256: return 60;
        case CFT_FP128: return 24;
        default:        return 12;
    }
}

int cft_ias15_format_code(const char *name){
    if (!name) return -1;
    if (!strcmp(name, "fp64"))  return CFT_FP64;
    if (!strcmp(name, "fp128")) return CFT_FP128;
    if (!strcmp(name, "fp256")) return CFT_FP256;
    return -1;
}

int cft_ias15_configure(struct reb_simulation *r, int format,
                        double epsilon, int max_iter){
    struct cft_ias15_state *st = cft_ias15_get_state(r);
    if (!st) return 1;
    if (format != CFT_FP64 && format != CFT_FP128 && format != CFT_FP256) return 2;
    if (epsilon < 0.0) return 3;
    if (max_iter < 0) return 4;
    if (max_iter == 0) max_iter = cft_ias15_default_max_iter(format);
    st->format   = format;
    st->epsilon  = epsilon;
    st->max_iter = max_iter;
    return 0;
}

void cft_ias15_reserve(size_t n){ reserve_N = n; }
void cft_ias15_set_artifact(const char *path){ artifact_path = path; }
void cft_ias15_set_pc_tol_shift(int shift){ pc_tol_shift = shift; }
uint32_t cft_ias15_flags_seen(void){ return ias15_engine_is_open() ? ias15_engine_flags() : 0; }
unsigned long long cft_ias15_library_calls(void){ return ias15_engine_is_open() ? ias15_engine_calls() : 0; }
