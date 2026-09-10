/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Phase 0: REBOUND's own IAS15, in plain binary64, on the reference
 * problems - the ground truth every later claim is measured against.
 *
 * Reads a problem file (data/problems/, exact binary64 hex floats),
 * runs REBOUND's IAS15 for a number of steps, and writes one
 * record line per sample with the full state as exact hex floats, so
 * that the libcft port can be compared bit for bit and an mpmath
 * oracle can score the energy from the exact bits.
 *
 *   ias15_ref --problem FILE [--dt DT] [--epsilon EPS] [--steps N]
 *             [--sample K] [--adaptive-mode M]
 *
 * epsilon = 0 turns adaptive stepping off (REBOUND's own convention),
 * so --dt is then the fixed step. Records go to stdout. */
#include "rebound.h"
#include "../src/hexfloat.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

struct body { char name[32]; double m, x, y, z, vx, vy, vz; };

static double need_hex(const char *tok, const char *what){
    double v;
    if (!tok || !hexfloat_parse(tok, &v)){
        fprintf(stderr, "ias15_ref: bad hex float for %s: %s\n", what, tok ? tok : "(null)");
        exit(2);
    }
    return v;
}

static int read_problem(const char *path, char *name, double *G, struct body *bodies, int max){
    FILE *f = fopen(path, "r");
    if (!f){ perror(path); exit(2); }
    char line[1024];
    int n = 0, N = -1;
    *G = 1.0;
    name[0] = 0;
    while (fgets(line, sizeof line, f)){
        if (line[0] == '#' || line[0] == '\n') continue;
        char *tok = strtok(line, " \t\r\n");
        if (!tok) continue;
        if (strcmp(tok, "name") == 0){ strncpy(name, strtok(NULL, " \t\r\n"), 31); name[31] = 0; }
        else if (strcmp(tok, "G") == 0){ *G = need_hex(strtok(NULL, " \t\r\n"), "G"); }
        else if (strcmp(tok, "N") == 0){ N = atoi(strtok(NULL, " \t\r\n")); }
        else if (strcmp(tok, "body") == 0){
            if (n >= max){ fprintf(stderr, "too many bodies\n"); exit(2); }
            struct body *b = &bodies[n++];
            strncpy(b->name, strtok(NULL, " \t\r\n"), 31); b->name[31] = 0;
            b->m  = need_hex(strtok(NULL, " \t\r\n"), "m");
            b->x  = need_hex(strtok(NULL, " \t\r\n"), "x");
            b->y  = need_hex(strtok(NULL, " \t\r\n"), "y");
            b->z  = need_hex(strtok(NULL, " \t\r\n"), "z");
            b->vx = need_hex(strtok(NULL, " \t\r\n"), "vx");
            b->vy = need_hex(strtok(NULL, " \t\r\n"), "vy");
            b->vz = need_hex(strtok(NULL, " \t\r\n"), "vz");
        }else{
            fprintf(stderr, "ias15_ref: unknown line '%s' in %s\n", tok, path);
            exit(2);
        }
    }
    fclose(f);
    if (N != n){ fprintf(stderr, "ias15_ref: N says %d, %d bodies read\n", N, n); exit(2); }
    return n;
}

static void put(const char *label, double x){
    char buf[40]; hexfloat_print(x, buf); printf(" %s%s", label, buf);
}

int main(int argc, char **argv){
    const char *problem = NULL;
    double dt = 0.01, epsilon = 1e-9;
    long steps = 1000, sample = 100;
    int adaptive_mode = 2;
    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "--problem") && i + 1 < argc) problem = argv[++i];
        else if (!strcmp(argv[i], "--dt") && i + 1 < argc) dt = atof(argv[++i]);
        else if (!strcmp(argv[i], "--epsilon") && i + 1 < argc) epsilon = atof(argv[++i]);
        else if (!strcmp(argv[i], "--steps") && i + 1 < argc) steps = atol(argv[++i]);
        else if (!strcmp(argv[i], "--sample") && i + 1 < argc) sample = atol(argv[++i]);
        else if (!strcmp(argv[i], "--adaptive-mode") && i + 1 < argc) adaptive_mode = atoi(argv[++i]);
        else { fprintf(stderr, "ias15_ref: unknown argument %s\n", argv[i]); return 2; }
    }
    if (!problem){ fprintf(stderr, "usage: ias15_ref --problem FILE [--dt DT] [--epsilon EPS] [--steps N] [--sample K]\n"); return 2; }

    char name[32]; double G; struct body bodies[64];
    int N = read_problem(problem, name, &G, bodies, 64);

    struct reb_simulation *r = reb_simulation_create();
    r->G = G;
    r->dt = dt;
    reb_simulation_set_integrator(r, "ias15");
    struct reb_integrator_ias15_state *ias15 = r->integrator.state;
    ias15->epsilon = epsilon;
    ias15->adaptive_mode = adaptive_mode;
    for (int i = 0; i < N; i++){
        struct reb_particle p = {0};
        p.m = bodies[i].m; p.x = bodies[i].x; p.y = bodies[i].y; p.z = bodies[i].z;
        p.vx = bodies[i].vx; p.vy = bodies[i].vy; p.vz = bodies[i].vz;
        reb_simulation_add(r, p);
    }

    printf("# cft-rebound ias15 record v1\n");
    printf("# program=ref impl=rebound format=fp64 problem=%s N=%d steps=%ld sample=%ld adaptive_mode=%d\n",
           name, N, steps, sample, adaptive_mode);
    { char a[40], b[40], c[40]; hexfloat_print(dt, a); hexfloat_print(epsilon, b); hexfloat_print(G, c);
      printf("# dt0=%s epsilon=%s G=%s\n", a, b, c); }
    for (int i = 0; i < N; i++){ char m[40]; hexfloat_print(bodies[i].m, m); printf("# body %d %s m=%s\n", i, bodies[i].name, m); }
    printf("# columns: sample step t dt_next dt_last E_rebound then per body x y z vx vy vz ; then exact-time pair t_hi t_lo\n");

    /* The exact elapsed time, kept beside REBOUND's own r->t. r->t is a
     * plain running sum of the steps taken and its rounding error grows
     * with the step count; it never feeds back into the trajectory (the
     * force and the step control are time-independent), so it is only a
     * label - but a label the oracle would mistake for a position error
     * once a run reaches 1e8 steps. The pair (t_hi, t_lo) accumulates
     * every dt_last_done with an exact TwoSum (Knuth), the same idea the
     * libcft port implements with 754-2019's augmentedAddition; -std=c99
     * pins -ffp-contract=off, so no FMA contraction disturbs it. Stepping
     * one step at a time is what reb_simulation_steps(r, todo) does
     * internally (a per-step loop, plus a synchronize IAS15 does not
     * implement), so the trajectory is unchanged. */
    double t_hi = 0.0, t_lo = 0.0;
    long done = 0, k = 0;
    while (1){
        printf("sample %ld %ld", k, done);
        put("", r->t); put("", r->dt); put("", r->dt_last_done); put("", reb_simulation_energy(r));
        for (int i = 0; i < N; i++){
            struct reb_particle *p = &r->particles[i];
            put("", p->x); put("", p->y); put("", p->z); put("", p->vx); put("", p->vy); put("", p->vz);
        }
        printf(" |"); put("", t_hi); put("", t_lo);
        printf("\n");
        fflush(stdout);
        if (done >= steps) break;
        long todo = sample; if (done + todo > steps) todo = steps - done;
        for (long s = 0; s < todo; s++){
            reb_simulation_steps(r, 1);
            double d = r->dt_last_done;
            double sum = t_hi + d;
            double bv = sum - t_hi;
            double err = (t_hi - (sum - bv)) + (d - bv);
            t_hi = sum;
            t_lo += err;
        }
        done += todo;
        k++;
    }
    printf("# steps_done=%llu iterations_max_exceeded=%llu\n",
           (unsigned long long)r->steps_done, (unsigned long long)ias15->iterations_max_exceeded);
    reb_simulation_free(r);
    return 0;
}
