/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * The packaging layer: cft_rebound_check() and cft_rebound_steps(),
 * declared in include/cft_rebound.h.
 *
 * THIS IS THE SUBPROCESS FORM. It writes the simulation's particles out
 * as an exact binary64 problem file, runs the `ias15_cft` program on it
 * at the requested format, and rounds the last recorded state back into
 * r->particles. Nothing is approximated on the way in - binary64 is a
 * subset of every wider format, and the problem file carries exact hex
 * floats - and the way out is one correctly rounded conversion per
 * value, done by libcft, which is the definition of the bits in this
 * project.
 *
 * IT IS STILL THE SUBPROCESS FORM, ON PURPOSE. ROADMAP parcel A landed
 * and the registration form exists - src/reb_integrator_cft.c,
 * src/cft_ias15.h, -lcft_ias15 - but this file was NOT rewritten as a
 * wrapper around it. The two are different products: there the engine
 * runs in the caller's address space and the wide state lives across
 * steps; here the caller keeps the engine out of its process and pays
 * one exec per call. A caller who wants the wide state to persist
 * should use the drop-in, which is what examples/dropin.c shows and
 * what the README's "The drop-in" section documents.
 *
 * What that means for a caller of THIS file is unchanged and is the
 * warning in include/cft_rebound.h: the wide state does not survive the
 * call, so ask for a whole run in one call.
 */

#include "cft_rebound.h"
#include "hexfloat.h"
#include "cft.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
#include <process.h>
#define CFT_GETPID _getpid
#else
#include <unistd.h>
#define CFT_GETPID getpid
#endif

/* Where `make install` puts the program: $(PREFIX)/bin/ias15_cft. The
 * root Makefile defines this; if the file is not actually there - the
 * library was built for one prefix and moved, or is being used straight
 * out of the build tree - the search falls through to PATH. */
#ifndef CFT_REBOUND_IAS15_DEFAULT
#define CFT_REBOUND_IAS15_DEFAULT "ias15_cft"
#endif

/* Kept in step with CFT_MAX_BODIES in src/ias15_cft.c. The program is
 * the authority and refuses on its own; this is here so that the
 * refusal arrives before a temporary file is written. */
#define CFT_REBOUND_MAX_BODIES 1024

const char *cft_rebound_format_name(enum cft_rebound_format f){
    switch (f){
        case CFT_REBOUND_FP64:  return "fp64";
        case CFT_REBOUND_FP128: return "fp128";
        case CFT_REBOUND_FP256: return "fp256";
    }
    return "?";
}

static int exists(const char *path){
    FILE *f = fopen(path, "rb");
    if (!f) return 0;
    fclose(f);
    return 1;
}

const char *cft_rebound_program_path(const struct cft_rebound_options *opt){
    const char *p;
    if (opt && opt->program && opt->program[0]) return opt->program;
    p = getenv("CFT_REBOUND_IAS15");
    if (p && p[0]) return p;
    if (exists(CFT_REBOUND_IAS15_DEFAULT)) return CFT_REBOUND_IAS15_DEFAULT;
    return "ias15_cft";     /* let PATH decide */
}

/* opt->artifact, then $CFT_REBOUND_ARTIFACT, then NULL for software.
 * src/reb_integrator_cft.c resolve_artifact() is the same rule for the
 * registered integrator; keep them together if either changes. */
static const char *artifact_of(const struct cft_rebound_options *opt){
    const char *p;
    if (opt && opt->artifact && opt->artifact[0]) return opt->artifact;
    p = getenv("CFT_REBOUND_ARTIFACT");
    if (p && p[0]) return p;
    return NULL;
}

static const char *workdir_of(const struct cft_rebound_options *opt){
    const char *p;
    if (opt && opt->workdir && opt->workdir[0]) return opt->workdir;
    if ((p = getenv("TMPDIR")) && p[0]) return p;
    if ((p = getenv("TEMP"))   && p[0]) return p;
    if ((p = getenv("TMP"))    && p[0]) return p;
    return ".";
}

static void say(char *why, size_t n, const char *msg){
    if (why && n){
        size_t k = strlen(msg);
        if (k >= n) k = n - 1;
        memcpy(why, msg, k);
        why[k] = 0;
    }
}

/* ------------------------------------------------------------------ */
/* The refusal list                                                    */
/* ------------------------------------------------------------------ */
/* The force model here is basic pairwise Newtonian gravity, full stop.
 * Everything below computes something DIFFERENT from what the port
 * would compute, so each is refused by name rather than silently
 * ignored. README "Scope" is this list in prose. */
int cft_rebound_check(const struct reb_simulation *r, char *why, size_t n){
    if (!r){ say(why, n, "no simulation"); return 1; }

    if (r->N_var != 0){
        say(why, n, "variational particles are not supported (r->N_var != 0)"); return 1; }
    {   size_t N = r->N;
        if (N < 1){ say(why, n, "no particles"); return 1; }
        if (N > CFT_REBOUND_MAX_BODIES){
            char b[160];
            sprintf(b, "%lu particles, and the port's limit is %d per system"
                       " (memory and time grow as N^2)",
                    (unsigned long)N, CFT_REBOUND_MAX_BODIES);
            say(why, n, b); return 1;
        }
    }
    if (r->N_active != (size_t)-1 && r->N_active != r->N){
        say(why, n, "test particles are not supported (r->N_active is set):"
                    " they change which pairs are computed"); return 1; }
    if (r->map != NULL || r->N_map != 0){
        say(why, n, "a particle map is not supported (r->map != NULL)"); return 1; }
    if (r->additional_forces != NULL){
        say(why, n, "additional forces are not supported (r->additional_forces is set)"); return 1; }
    if (r->pre_timestep_modifications != NULL || r->post_timestep_modifications != NULL){
        say(why, n, "timestep modifications are not supported (REBOUNDx and the like)"); return 1; }
    if (r->force_is_velocity_dependent){
        say(why, n, "velocity-dependent forces are not supported"); return 1; }
    if (r->gravity != REB_GRAVITY_BASIC){
        say(why, n, r->gravity == REB_GRAVITY_TREE
                ? "the tree code is not supported (r->gravity = TREE); this port is direct summation"
                : (r->gravity == REB_GRAVITY_COMPENSATED
                   ? "r->gravity = COMPENSATED is a different summation from the one ported; use BASIC"
                   : "only r->gravity = BASIC is supported"));
        return 1; }
    if (r->gravity_custom != NULL){
        say(why, n, "a custom gravity routine is not supported"); return 1; }
    if (r->gravity_ignore_terms != REB_GRAVITY_IGNORE_TERMS_NONE){
        say(why, n, "r->gravity_ignore_terms must be NONE; this port computes every pair"); return 1; }
    if (r->collision != REB_COLLISION_NONE){
        say(why, n, "collision detection is not supported (r->collision != NONE)"); return 1; }
    if (r->boundary != REB_BOUNDARY_NONE){
        say(why, n, "boundary conditions are not supported (r->boundary != NONE)"); return 1; }
    if (r->N_ghost_x || r->N_ghost_y || r->N_ghost_z){
        say(why, n, "ghost boxes are not supported (r->N_ghost_* != 0)"); return 1; }
    /* No check on r->OMEGA. It reaches the dynamics only through the
     * shearing sheet (boundary SHEAR) or the SEI integrator, both
     * refused above, and REBOUND initialises r->OMEGAZ to -1 as a
     * sentinel meaning "use OMEGA" - so a test on either of them
     * refuses every default simulation, which this one did until the
     * worked example ran. */
    if (r->N_odes != 0){
        say(why, n, "attached ODE sets are not supported"); return 1; }
    if (r->integrator.name && strcmp(r->integrator.name, "ias15") != 0){
        char b[160];
        sprintf(b, "the simulation's integrator is \"%.80s\"; this port is IAS15 only",
                r->integrator.name);
        say(why, n, b); return 1; }
    if (r->integrator.state && r->integrator.name && strcmp(r->integrator.name, "ias15") == 0){
        const struct reb_integrator_ias15_state *s = r->integrator.state;
        if (s->adaptive_mode != 2){
            say(why, n, "only IAS15's PRS23 step criterion (adaptive_mode 2) is ported"); return 1; }
        if (s->min_dt != 0.0){
            say(why, n, "IAS15's min_dt is not ported; leave it at 0"); return 1; }
    }
    say(why, n, "");
    return 0;
}

/* ------------------------------------------------------------------ */
/* Reading a record back: wide hex text -> correctly rounded binary64  */
/* ------------------------------------------------------------------ */
/* src/hexfloat.h reads binary64 hex text only, and refuses anything
 * that would need rounding - deliberately, because a binary64 record
 * must never hold a value that does. A binary128 or binary256 record
 * does, on almost every line, so the conversion here goes through
 * libcft's own cft_from_hex_char at CFT_FP64 with roundTiesToEven. That
 * is the "binary64 VIEW" of ROADMAP.md, produced by the library that
 * defines what the bits mean rather than by a parser written twice. */
static int widehex_to_double(cft_device *dev, const char *const *in, double *out,
                             size_t n, const char **bad){
    size_t bad_index = 0;
    uint32_t flags = 0;
    cft_status st = cft_from_hex_char(dev, CFT_FP64, CFT_RNE, in, out, n, &bad_index, &flags);
    if (st != CFT_OK){ if (bad) *bad = in[bad_index]; return 1; }
    return 0;
}

/* the tokens of the last "sample" line, and the trailing summary */
struct record {
    char  **tok;      /* tokens of the last sample line */
    size_t  ntok;
    char    tail[512];
};

static void record_free(struct record *rec){
    size_t i;
    if (!rec->tok) return;
    for (i = 0; i < rec->ntok; i++) free(rec->tok[i]);
    free(rec->tok);
    rec->tok = NULL; rec->ntok = 0;
}

static int record_read(const char *path, struct record *rec, size_t nfields){
    FILE *f = fopen(path, "r");
    char *line;
    size_t cap = 1 << 16;
    memset(rec, 0, sizeof *rec);
    if (!f) return 1;
    line = malloc(cap);
    if (!line){ fclose(f); return 1; }
    while (fgets(line, (int)cap, f)){
        if (!strncmp(line, "sample ", 7)){
            char *save = line, *t;
            size_t k = 0;
            record_free(rec);
            rec->tok = malloc(nfields * sizeof *rec->tok);
            if (!rec->tok) break;
            for (t = strtok(save, " \t\r\n"); t; t = strtok(NULL, " \t\r\n")){
                if (!strcmp(t, "|")) break;          /* the exact-time pair follows */
                if (k >= nfields) continue;
                rec->tok[k] = malloc(strlen(t) + 1);
                if (!rec->tok[k]) break;
                strcpy(rec->tok[k], t);
                k++;
            }
            /* the two exact-time words after the bar, appended */
            for (t = strtok(NULL, " \t\r\n"); t && k < nfields; t = strtok(NULL, " \t\r\n")){
                rec->tok[k] = malloc(strlen(t) + 1);
                if (!rec->tok[k]) break;
                strcpy(rec->tok[k], t);
                k++;
            }
            rec->ntok = k;
        }else if (!strncmp(line, "# steps_done=", 13)){
            strncpy(rec->tail, line + 2, sizeof rec->tail - 1);
            rec->tail[sizeof rec->tail - 1] = 0;
        }
    }
    free(line);
    fclose(f);
    return rec->ntok ? 0 : 1;
}

static long tail_long(const char *tail, const char *key, long dflt){
    const char *p = strstr(tail, key);
    return p ? atol(p + strlen(key)) : dflt;
}
static double tail_double(const char *tail, const char *key, double dflt){
    const char *p = strstr(tail, key);
    return p ? atof(p + strlen(key)) : dflt;
}

/* ------------------------------------------------------------------ */
/* The call                                                            */
/* ------------------------------------------------------------------ */
int cft_rebound_steps(struct reb_simulation *r, long nsteps,
                      const struct cft_rebound_options *opt,
                      struct cft_rebound_result *res){
    struct cft_rebound_options defaults;
    char why[256];
    char probpath[1024], recpath[1024], cmd[4096];
    const char *dir, *prog, *art;
    FILE *f;
    size_t i, N, nfields, nvals;
    cft_device *dev = NULL;
    cft_status st;
    struct record rec;
    double *vals = NULL;
    const char **in = NULL;
    const char *bad = NULL;
    int rc = 1;

    if (!opt){ memset(&defaults, 0, sizeof defaults); opt = &defaults; }
    if (nsteps < 0){ fprintf(stderr, "cft_rebound: nsteps must not be negative\n"); return 1; }
    if (cft_rebound_check(r, why, sizeof why)){
        fprintf(stderr, "cft_rebound: %s\n", why);
        return 1;
    }
    N = r->N;
    dir  = workdir_of(opt);
    prog = cft_rebound_program_path(opt);
    art  = artifact_of(opt);
    sprintf(probpath, "%.900s/cft_rebound_%d.problem", dir, (int)CFT_GETPID());
    sprintf(recpath,  "%.900s/cft_rebound_%d.record",  dir, (int)CFT_GETPID());
    /* prog comes from the caller or the environment and is otherwise
     * unbounded; 200 covers the flags and their arguments. Fail rather
     * than truncate - a truncated path is a wrong path. */
    if (strlen(prog) + strlen(probpath) + strlen(recpath)
        + (art ? strlen(art) : 0) + 200 > sizeof cmd){
        fprintf(stderr, "cft_rebound: the command line would not fit in %lu bytes;"
                        " use a shorter program or workdir path\n", (unsigned long)sizeof cmd);
        return 1;
    }

    /* the problem: exact binary64 hex, so nothing is rounded going in */
    f = fopen(probpath, "w");
    if (!f){ fprintf(stderr, "cft_rebound: cannot write %s\n", probpath); return 1; }
    fprintf(f, "# written by cft_rebound_steps()\n");
    { char b[40]; hexfloat_print(r->G, b); fprintf(f, "name rebound\nG %s\nN %lu\n", b, (unsigned long)N); }
    for (i = 0; i < N; i++){
        const struct reb_particle *p = &r->particles[i];
        char m[40], x[40], y[40], z[40], vx[40], vy[40], vz[40];
        hexfloat_print(p->m, m);   hexfloat_print(p->x, x);   hexfloat_print(p->y, y);
        hexfloat_print(p->z, z);   hexfloat_print(p->vx, vx); hexfloat_print(p->vy, vy);
        hexfloat_print(p->vz, vz);
        fprintf(f, "body p%lu %s %s %s %s %s %s %s\n", (unsigned long)i, m, x, y, z, vx, vy, vz);
    }
    fclose(f);

    /* the run. dt and epsilon go across as exact hex floats, so the
     * program starts from the same bits this process holds. */
    {
        char dtb[40], epsb[40];
        int k;
        hexfloat_print(r->dt, dtb);
        hexfloat_print(opt->epsilon, epsb);
        k = sprintf(cmd, "\"%s\" --format %s --problem \"%s\" --dt %s --epsilon %s"
                         " --steps %ld --sample %ld --quiet",
                    prog, cft_rebound_format_name(opt->format), probpath, dtb, epsb,
                    nsteps, nsteps > 0 ? nsteps : 1);
        /* The corrector's pass cap. ias15_cft defaults to REBOUND's 12 at
         * every format, and that is right at binary64 - it is part of what
         * makes the binary64 run REBOUND's own, bit for bit - but wrong
         * above it: docs/VALIDATION.md measures 2-3 passes a step at
         * binary64, 4-9 at binary128 and 10-20 at binary256, so at
         * binary256 a cap of 12 is hit on essentially every step and the
         * corrector is truncated rather than converged. The caps below sit
         * clear of those measured counts (60 is what this repository's own
         * binary256 runs use); a converged step exits early, so a cap only
         * costs anything when it is actually reached, and the result says
         * when it was. */
        {
            int mi = opt->max_iter;
            if (mi <= 0) mi = (opt->format == CFT_REBOUND_FP256) ? 60
                            : (opt->format == CFT_REBOUND_FP128) ? 24
                            : 0;    /* binary64: leave REBOUND's 12 alone */
            if (mi > 0) k += sprintf(cmd + k, " --max-iter %d", mi);
        if (r->softening != 0.0){
            char sb[40];
            hexfloat_print(r->softening, sb);
            k += sprintf(cmd + k, " --softening %s", sb);
        }
        }
        if (art) k += sprintf(cmd + k, " --artifact \"%s\"", art);
        k += sprintf(cmd + k, " > \"%s\"", recpath);
        (void)k;
    }
    if (opt->verbose) fprintf(stderr, "cft_rebound: %s\n", cmd);
    {
        /* cmd.exe eats the outer pair of quotes when the string starts
         * with one, so a quoted program path needs a second pair. */
        int status;
#ifdef _WIN32
        char wrapped[4200];
        sprintf(wrapped, "\"%s\"", cmd);
        status = system(wrapped);
#else
        status = system(cmd);
#endif
        if (status != 0){
            fprintf(stderr, "cft_rebound: %s failed (see the record file %s)\n", prog, recpath);
            fprintf(stderr, "cft_rebound: is ias15_cft on PATH, or CFT_REBOUND_IAS15 set?\n");
            goto out;
        }
    }

    /* "sample" k step | t dt_next dt_last E | 6 per body | t_hi t_lo:
     * three words of label, then every number of the line. */
    nfields = 3 + 4 + 6 * N + 2;
    nvals   = nfields - 3;
    if (record_read(recpath, &rec, nfields)){
        fprintf(stderr, "cft_rebound: no record in %s\n", recpath);
        goto out;
    }
    if (rec.ntok != nfields){
        fprintf(stderr, "cft_rebound: record has %lu fields, expected %lu\n",
                (unsigned long)rec.ntok, (unsigned long)nfields);
        goto out;
    }

    st = cft_open(NULL, 0, &dev);   /* software backend: no card is touched */
    if (st != CFT_OK){
        fprintf(stderr, "cft_rebound: cft_open: %s\n", cft_strerror(st));
        goto out;
    }
    /* every hex word of the line at once, correctly rounded to binary64 */
    vals = malloc(nvals * sizeof *vals);
    in   = malloc(nvals * sizeof *in);
    if (!vals || !in){ fprintf(stderr, "cft_rebound: out of memory\n"); goto out; }
    for (i = 3; i < nfields; i++) in[i - 3] = rec.tok[i];
    if (widehex_to_double(dev, in, vals, nvals, &bad)){
        fprintf(stderr, "cft_rebound: cannot read \"%s\" from the record\n", bad ? bad : "?");
        goto out;
    }

    /* the binary64 view, back into REBOUND */
    for (i = 0; i < N; i++){
        struct reb_particle *p = &r->particles[i];
        const double *b = vals + 4 + 6 * i;
        p->x  = b[0]; p->y  = b[1]; p->z  = b[2];
        p->vx = b[3]; p->vy = b[4]; p->vz = b[5];
    }
    r->t            = vals[0];
    r->dt           = vals[1];
    r->dt_last_done = vals[2];
    r->steps_done  += (uint64_t)atol(rec.tok[2]);

    if (res){
        memset(res, 0, sizeof *res);
        res->t        = vals[0];
        res->dt_next  = vals[1];
        res->dt_last  = vals[2];
        res->energy   = vals[3];
        strncpy(res->energy_hex, rec.tok[6], sizeof res->energy_hex - 1);
        res->energy_hex[sizeof res->energy_hex - 1] = 0;
        res->t_hi     = vals[nvals - 2];
        res->t_lo     = vals[nvals - 1];
        res->steps_done          = tail_long(rec.tail, "steps_done=", 0);
        res->steps_rejected      = tail_long(rec.tail, "steps_rejected=", 0);
        res->iterations_max_exceeded = tail_long(rec.tail, "iterations_max_exceeded=", 0);
        res->mean_pc_iterations  = tail_double(rec.tail, "mean_pc_iterations=", 0.0);
        res->max_pc_iterations   = (int)tail_long(rec.tail, "max_pc_iterations=", 0);
        res->seconds             = tail_double(rec.tail, "seconds=", 0.0);
    }
    rc = 0;

out:
    if (dev) cft_close(dev);
    free(vals);
    free((void *)in);
    record_free(&rec);
    remove(probpath);
    if (rc == 0) remove(recpath);
    return rc;
}
