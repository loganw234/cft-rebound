/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * IAS15 with every floating-point operation issued through libcft.
 *
 * This is a re-implementation of the integrator published in REBOUND's
 * src/integrator_ias15.c (Rein & Spiegel 2015; Everhart 1985; the PRS23
 * step criterion of Pham, Rein & Spiegel 2024) in which no C arithmetic
 * operator ever touches a floating-point value. Every add, multiply,
 * divide, square root, comparison and classification is a cft.h call at
 * one of binary64, binary128 or binary256, on libcft's software backend
 * or, given an artifact, on the tile - the same bits either way.
 *
 * THE SEQUENCE OF ROUNDINGS IS REBOUND'S. The order of every operation
 * follows integrator_ias15.c and gravity.c exactly: C evaluates
 * `b6*7.*h/9.` as ((b6*7)*h)/9 and the code below issues those three
 * roundings in that order, no FMA is ever used where REBOUND has a
 * product and a sum, and the pairwise-gravity accumulation adds each
 * particle's partners in the order REBOUND's half-N^2 loop visits them.
 * The reason is the gate: at binary64 this program must agree with
 * REBOUND's own IAS15 bit for bit (tools/check_equivalence.py), and only
 * then is raising the format a controlled experiment. Where REBOUND's
 * double literals are integers or dyadic rationals they are recreated
 * here from their definitions (binomial coefficients, (j+1)/(j+3), and
 * so on), not typed; where they are the Gauss-Radau constants they come
 * from src/ias15_constants.h, derived at 130 digits.
 *
 * ENSEMBLES. A problem file may hold E independent systems of N bodies
 * (`system NAME` lines, docs/ENSEMBLE.md). They are integrated in ONE
 * run with their coordinates side by side in every vector - 3NE lanes
 * for the predictor and corrector, E N(N-1)/2 lanes for the gravity -
 * so that a two-body problem, one pair on its own, becomes E pairs.
 * Every system keeps its own adaptive step, its own corrector exit and
 * its own accept/reject decision: the vectors are shared, the scalar
 * control is per system, and a system that has left the corrector (or
 * rejected its step, or finished its sample block) is masked by a
 * snapshot taken at that moment and restored after the shared vector
 * work, byte for byte. Elementwise arithmetic is per element, so the
 * gate is exact: an ensemble of E systems must reproduce, bit for bit,
 * the records of those E systems run one at a time (`--member k`,
 * tools/check_ensemble.py). What the shared vector costs is idle lanes
 * whenever the systems disagree about how many corrector passes they
 * need, reported as pc_lane_efficiency in the trailer.
 *
 * WHAT IS NOT MIRRORED: variational particles, MEGNO, velocity-dependent
 * forces, ghost boxes, the tree code, softening other than zero, and the
 * pre-2024 GLOBAL/INDIVIDUAL step criteria. None of those touch the
 * arithmetic of the scheme itself.
 *
 * TWO KNOBS ARE FORMAT-DEPENDENT BY NECESSITY. The predictor-corrector
 * loop stops when its error estimate falls under 1e-16, REBOUND's
 * "machine precision"; here that threshold is 1e-16 * 2^(53-p), i.e. the
 * same number of ulps in every format (--pc-tol-shift overrides the
 * exponent). And --cs augmented replaces Kahan's compensated add_cs()
 * with the exact error term of IEEE 754-2019 9.5's augmentedAddition.
 *
 * Usage:
 *   ias15_cft --format fp64|fp128|fp256 --problem FILE [--dt DT]
 *             [--epsilon EPS] [--steps N] [--sample K] [--cs kahan|augmented]
 *             [--pc-tol-shift S] [--max-iter M] [--no-flag-abort]
 *             [--artifact PATH] [--dump-constants] [--quiet]
 *             [--member K] [--dt-file FILE] [--dt-out FILE]
 */
#include "cft.h"
#include "ias15_constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdint.h>
#include <time.h>

/* ------------------------------------------------------------------ */
/* The device, the format, and the accounting                          */
/* ------------------------------------------------------------------ */
static cft_device *dev;
static cft_format  F;
static int         FI;          /* 0 fp64, 1 fp128, 2 fp256 */
static size_t      ESZ;
static uint32_t    flags_union;
static int         flag_abort = 1;
static unsigned long long ncalls, ncalls_divsqrt;

static void die(const char *fmt, ...){
    va_list ap; va_start(ap, fmt);
    fprintf(stderr, "ias15_cft: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    exit(3);
}

#define CERT_FLAGS (CFT_FLAG_INVALID | CFT_FLAG_DIVBYZERO | CFT_FLAG_OVERFLOW | CFT_FLAG_UNDERFLOW)

static void note(cft_status st, uint32_t fl, const char *what){
    ncalls++;
    if (st != CFT_OK) die("%s: %s (%s)", what, cft_strerror(st), cft_last_error());
    flags_union |= fl;
    if (flag_abort && (fl & CERT_FLAGS))
        die("%s raised 0x%02x - this workload can only raise inexact; invalid, "
            "divide-by-zero, overflow or underflow means the integration or the "
            "port is wrong (pass --no-flag-abort to continue anyway)", what, fl);
}

typedef unsigned char *V;
#define E(v, i) ((v) + (size_t)(i) * ESZ)

static V valloc(size_t n){
    V v = calloc(n ? n : 1, ESZ);
    if (!v) die("out of memory");
    return v;
}

/* ------------------------------------------------------------------ */
/* Elementwise operations. d may alias any input (cft_run's rule).      */
/* ------------------------------------------------------------------ */
static void vfma(V d, const V a, const V b, const V c, size_t n){
    uint32_t fl = 0; cft_status st = cft_run(dev, CFT_FMA, F, CFT_RNE, a, b, c, d, n, &fl, NULL); note(st, fl, "fma"); }
static void vadd(V d, const V a, const V c, size_t n){
    uint32_t fl = 0; cft_status st = cft_run(dev, CFT_ADD, F, CFT_RNE, a, NULL, c, d, n, &fl, NULL); note(st, fl, "add"); }
static void vsub(V d, const V a, const V c, size_t n){
    uint32_t fl = 0; cft_status st = cft_run(dev, CFT_SUB, F, CFT_RNE, a, NULL, c, d, n, &fl, NULL); note(st, fl, "sub"); }
static void vmul(V d, const V a, const V b, size_t n){
    uint32_t fl = 0; cft_status st = cft_run(dev, CFT_MUL, F, CFT_RNE, a, b, NULL, d, n, &fl, NULL); note(st, fl, "mul"); }
static void vneg(V d, const V a, size_t n){
    uint32_t fl = 0; cft_status st = cft_run(dev, CFT_NEG, F, CFT_RNE, a, NULL, NULL, d, n, &fl, NULL); note(st, fl, "neg"); }
static void vabs(V d, const V a, size_t n){
    uint32_t fl = 0; cft_status st = cft_run(dev, CFT_ABS, F, CFT_RNE, a, NULL, NULL, d, n, &fl, NULL); note(st, fl, "abs"); }
static void vdiv(V d, const V a, const V b, size_t n){
    uint32_t fl = 0; ncalls_divsqrt++; cft_status st = cft_div(dev, F, CFT_RNE, a, b, d, n, &fl, NULL); note(st, fl, "div"); }
static void vsqrt(V d, const V a, size_t n){
    uint32_t fl = 0; ncalls_divsqrt++; cft_status st = cft_sqrt(dev, F, CFT_RNE, a, d, n, &fl, NULL); note(st, fl, "sqrt"); }
static void vclass(uint8_t *cls, const V a, size_t n){
    ncalls++; cft_status st = cft_class(dev, F, a, cls, n);
    if (st != CFT_OK) die("class: %s", cft_strerror(st)); }
static void vaugadd(V r, V e, const V a, const V b, size_t n){
    uint32_t fl = 0; cft_status st = cft_augmented_add(dev, F, a, b, r, e, n, &fl); note(st, fl, "augmented_add"); }

/* Scalars are vectors of one element. A predicate result is 1.0 or
 * +0.0; "nonzero bits" reads it without arithmetic. */
static int is_nonzero_bits(const V a){
    for (size_t i = 0; i < ESZ; i++) if (a[i]) return 1;
    return 0;
}
static int s_lt(const V a, const V b){ /* a < b */
    static V t; if (!t) t = valloc(1);
    uint32_t fl = 0; note(cft_run(dev, CFT_CMPLT, F, CFT_RNE, a, b, NULL, t, 1, &fl, NULL), fl, "cmplt");
    return is_nonzero_bits(t);
}
static int s_le(const V a, const V b){
    static V t; if (!t) t = valloc(1);
    uint32_t fl = 0; note(cft_run(dev, CFT_CMPLE, F, CFT_RNE, a, b, NULL, t, 1, &fl, NULL), fl, "cmple");
    return is_nonzero_bits(t);
}
static int s_eq(const V a, const V b){
    static V t; if (!t) t = valloc(1);
    uint32_t fl = 0; note(cft_run(dev, CFT_CMPEQ, F, CFT_RNE, a, b, NULL, t, 1, &fl, NULL), fl, "cmpeq");
    return is_nonzero_bits(t);
}
static int is_normal_class(uint8_t c){
    return c == CFT_CLASS_NEG_NORM || c == CFT_CLASS_POS_NORM;
}
static int s_isnormal(const V a){
    uint8_t c; vclass(&c, a, 1); return is_normal_class(c);
}

/* Host-side data movement: byte copies, never arithmetic. */
static void vcopy(V d, const V s, size_t n){ memcpy(d, s, n * ESZ); }
static void vbcast(V d, const V s, size_t n){ for (size_t i = 0; i < n; i++) memcpy(E(d, i), s, ESZ); }
static void vzero(V d, size_t n){ memset(d, 0, n * ESZ); }  /* +0 in every format */

/* ------------------------------------------------------------------ */
/* Constants: parsed, never typed, and cross-checked                   */
/* ------------------------------------------------------------------ */
static size_t NMAX;  /* every broadcast constant has this many elements */

static V from_text(const char *s, int hex){
    V v = valloc(1);
    const char *in[1] = { s };
    size_t bad = 0; uint32_t fl = 0;
    cft_status st = hex ? cft_from_hex_char(dev, F, CFT_RNE, in, v, 1, &bad, &fl)
                        : cft_from_decimal_char(dev, F, CFT_RNE, in, v, 1, &bad, &fl);
    if (st != CFT_OK) die("cannot parse '%s': %s", s, cft_strerror(st));
    return v;
}
static int is_hex_text(const char *s){
    if (s && (*s == '-' || *s == '+')) s++;
    return s && s[0] == '0' && (s[1] == 'x' || s[1] == 'X');
}
static V kdec(const char *dec){  /* a broadcast constant from its decimal text */
    V s = from_text(dec, 0);
    V v = valloc(NMAX); vbcast(v, s, NMAX); free(s);
    return v;
}
static V kint(long x){ char buf[32]; sprintf(buf, "%ld", x); return kdec(buf); }

/* The derived Gauss-Radau constants: hex is the exact form, decimal the
 * 100-digit form; both must round to the same bits or the header lies. */
static V kderived(const char *hex, const char *dec, const char *name, int i, int *nchecked){
    V a = from_text(hex, 1);
    V b = from_text(dec, 0);
    if (memcmp(a, b, ESZ) != 0)
        die("constant %s[%d]: hex form %s and decimal form %s round to different %s bits",
            name, i, hex, dec, cft_format_name(F));
    (*nchecked)++;
    V v = valloc(NMAX); vbcast(v, a, NMAX); free(a); free(b);
    return v;
}

static V KH[8], KRR[28], KC[21], KD[21];
static V K0, K1, K2, K3, K4, K5, K6, K7, K10, K20, K0_1, K0_25, K1E2, K1E7, K1EM7, K5040, K1E300, K1E16, KINF;
static V KNUM[8], KDEN[8];      /* the predictor's (j+1)/(j+3) per level, j = 6..0 and the a0 level */
static V KPOS[8], KVEL[8];      /* the end-of-step divisors */
static V KY3[7], KY4[7], KY5[7];/* PRS23 coefficients (m+1), m(m+1), (m-1)m(m+1) */
static V KBIN[7][7];            /* binomial(j+1, m+1) for predict_next_step */
static V TOL;                   /* 1e-16 * 2^(53-p) */
/* the FMA form (--arith fma): reciprocals and per-substep factors, each
 * an exact rational rounded once */
static V KRINV[28];             /* 1/rr */
static V KHF[7][7];             /* h_n (j+1)/(j+3), [n-1][lvl], j = 6 - lvl */
static V KPOSR[8], KVELR[8];    /* 1/((m+2)(m+3)), 1/(m+2); [7] is a0's 1/2 and 1 */
static V KHALF;
static int arith_fma = 0;       /* 0: REBOUND's divisions; 1: the FMA form */
static int engine_program = 0;  /* 1: predictor and corrector as sequencer programs */
static const char *programs_dir = "programs/out";

static long binom(long n, long k){ long r = 1; for (long i = 1; i <= k; i++) r = r * (n - k + i) / i; return r; }

static void make_constants(int tol_shift, int quiet){
    int nchecked = 0;
    for (int i = 0; i < 8;  i++) KH[i]  = kderived(ias15_h_hex[FI][i],  ias15_h_dec[i],  "h",  i, &nchecked);
    for (int i = 0; i < 28; i++) KRR[i] = kderived(ias15_rr_hex[FI][i], ias15_rr_dec[i], "rr", i, &nchecked);
    for (int i = 0; i < 21; i++) KC[i]  = kderived(ias15_c_hex[FI][i],  ias15_c_dec[i],  "c",  i, &nchecked);
    for (int i = 0; i < 21; i++) KD[i]  = kderived(ias15_d_hex[FI][i],  ias15_d_dec[i],  "d",  i, &nchecked);
    for (int i = 0; i < 28; i++) KRINV[i] = kderived(ias15_rinv_hex[FI][i], ias15_rinv_dec[i], "rinv", i, &nchecked);
    for (int n = 0; n < 7; n++) for (int l = 0; l < 7; l++)
        KHF[n][l] = kderived(ias15_hf_hex[FI][7 * n + l], ias15_hf_dec[7 * n + l], "hf", 7 * n + l, &nchecked);
    K0 = kint(0); K1 = kint(1); K2 = kint(2); K3 = kint(3); K4 = kint(4); K5 = kint(5); K6 = kint(6); K7 = kint(7);
    K10 = kint(10); K20 = kint(20); K5040 = kint(5040);
    K0_1 = kdec("0.1"); K0_25 = kdec("0.25"); K1E2 = kdec("1e2"); K1E7 = kdec("1e7"); K1EM7 = kdec("1e-7");
    K1E300 = kdec("1e300"); K1E16 = kdec("1e-16"); KINF = kdec("inf");
    /* predictor levels: index 0 is the b6 level ... 6 is the b0 level, 7 the a0 level.
     * The Taylor factor at the b_j level is (j+1)/(j+3); REBOUND writes it
     * reduced (3/4 for 6/8, 2/3 for 4/6, 1/2 for 2/4), so reduce by the gcd
     * here too and the sequence of roundings is literally the same. */
    for (int lvl = 0; lvl < 8; lvl++){
        long j = 6 - lvl;              /* b_j at this level; j = -1 is the a0 level */
        long num = j + 1, den = j + 3, g2 = num, h2 = den;
        while (h2){ long r = g2 % h2; g2 = h2; h2 = r; }
        KNUM[lvl] = kint(num / g2);
        KDEN[lvl] = kint(den / g2);
    }
    for (int m = 0; m < 7; m++){
        KPOS[m] = kint((m + 2) * (m + 3));   /* b_m / ((m+2)(m+3)) dt^2 */
        KVEL[m] = kint(m + 2);               /* b_m / (m+2) dt */
        KY3[m] = kint(m + 1);
        KY4[m] = kint((long)m * (m + 1));
        KY5[m] = kint((long)(m - 1) * m * (m + 1));
    }
    KPOS[7] = kint(2);   /* a0 / 2 dt^2 */
    KVEL[7] = kint(1);   /* a0 dt */
    /* the reciprocals of those small integers, correctly rounded by the
     * library's own divide, and one half */
    for (int m = 0; m < 8; m++){
        V s = valloc(1);
        vdiv(s, K1, KPOS[m], 1); KPOSR[m] = valloc(NMAX); vbcast(KPOSR[m], s, NMAX);
        vdiv(s, K1, KVEL[m], 1); KVELR[m] = valloc(NMAX); vbcast(KVELR[m], s, NMAX);
        free(s);
    }
    { V s = valloc(1); vdiv(s, K1, K2, 1); KHALF = valloc(NMAX); vbcast(KHALF, s, NMAX); free(s); }
    for (int m = 0; m < 7; m++)
        for (int j = 0; j < 7; j++)
            KBIN[m][j] = (j > m) ? kint(binom(j + 1, m + 1)) : NULL;
    /* the convergence tolerance: 1e-16 scaled by an exact power of two */
    TOL = valloc(NMAX);
    { uint32_t fl = 0; V s = valloc(1);
      cft_status st = cft_scaleb(dev, F, CFT_RNE, K1E16, (int64_t)(-tol_shift), s, 1, &fl, NULL);
      if (st != CFT_OK) die("scaleb: %s", cft_strerror(st));
      if (fl) die("the tolerance scaling was inexact (0x%x) - that must never happen", fl);
      vbcast(TOL, s, NMAX); free(s); }
    if (!quiet) fprintf(stderr, "constants: %d derived values, hex and decimal forms agree at %s\n",
                        nchecked, cft_format_name(F));
}

/* ------------------------------------------------------------------ */
/* The problem: E independent systems of N bodies each                 */
/* ------------------------------------------------------------------ */
static size_t E = 1;       /* systems in this run */
static size_t N, NB, N3;   /* bodies per system, bodies in all, coordinates in all */
static size_t L;           /* lanes per system, 3N */
static size_t PS, P;       /* pairs per system, pairs in all */
static int member = -1;    /* --member k: run system k of an ensemble file alone */
static V mass;             /* NB */
static V G;                /* broadcast */
static V X0, V0;           /* the initial condition as read (N3) */
static char problem_name[64];
static char (*body_names)[32];   /* NB */
static char (*sys_names)[32];    /* E */

static void read_problem(const char *path){
    FILE *f = fopen(path, "r");
    if (!f){ perror(path); exit(2); }
    char line[2048];
    /* first pass: count bodies and systems, find N */
    long nb = 0, ns = 0, nline = -1;
    while (fgets(line, sizeof line, f)){
        if (!strncmp(line, "body ", 5)) nb++;
        else if (!strncmp(line, "system ", 7) || !strcmp(line, "system\n")) ns++;
        else if (!strncmp(line, "N ", 2)) nline = atol(line + 2);
    }
    rewind(f);
    size_t Efile = ns ? (size_t)ns : 1;
    if (nline < 0){ if (Efile > 1) die("an ensemble file needs an N line (bodies per system)"); nline = nb; }
    N = (size_t)nline;
    if (N < 1 || N > 64) die("N = %zu bodies per system is outside 1..64", N);
    if ((size_t)nb != Efile * N) die("%ld body lines for %zu system(s) of %zu bodies in %s", nb, Efile, N, path);
    if (member >= 0){
        if ((size_t)member >= Efile) die("--member %d, but %s has %zu system(s)", member, path, Efile);
        E = 1;
    }else E = Efile;
    NB = E * N; N3 = 3 * NB; L = 3 * N; PS = N * (N - 1) / 2; P = E * PS;
    /* everything is allocated at the largest count any call uses */
    NMAX = N3 > P ? N3 : P; if (NMAX < 8) NMAX = 8;
    mass = valloc(NB);
    G = valloc(NMAX);
    X0 = valloc(N3); V0 = valloc(N3);
    body_names = calloc(NB, sizeof *body_names);
    sys_names = calloc(E, sizeof *sys_names);
    if (!body_names || !sys_names) die("out of memory");
    /* second pass: parse. Hex floats are exact in every format. */
    long cur = -1;          /* the system the bodies being read belong to */
    long in_cur = 0;        /* bodies read in that system */
    size_t i = 0;           /* bodies kept */
    int have_G = 0;
    while (fgets(line, sizeof line, f)){
        if (line[0] == '#' || line[0] == '\n') continue;
        char *tok = strtok(line, " \t\r\n");
        if (!tok) continue;
        if (!strcmp(tok, "name")){ strncpy(problem_name, strtok(NULL, " \t\r\n"), 63); }
        else if (!strcmp(tok, "G")){ V s = from_text(strtok(NULL, " \t\r\n"), 1); vbcast(G, s, NMAX); free(s); have_G = 1; }
        else if (!strcmp(tok, "N")){ /* counted in the first pass */ }
        else if (!strcmp(tok, "E")){ if ((size_t)atol(strtok(NULL, " \t\r\n")) != Efile) die("E line disagrees with the system count"); }
        else if (!strcmp(tok, "system")){
            if (cur >= 0 && in_cur != (long)N) die("system %ld has %ld bodies, not %zu", cur, in_cur, N);
            cur++; in_cur = 0;
            const char *nm = strtok(NULL, " \t\r\n");
            if (member < 0) strncpy(sys_names[cur], nm ? nm : "", 31);
            else if (cur == member) strncpy(sys_names[0], nm ? nm : "", 31);
        }
        else if (!strcmp(tok, "body")){
            if (ns && cur < 0) die("body line before the first system line in %s", path);
            in_cur++;
            if (member >= 0 && cur != member){ continue; }
            if (i >= NB) die("too many body lines");
            strncpy(body_names[i], strtok(NULL, " \t\r\n"), 31);
            const char *fields[7]; for (int k = 0; k < 7; k++){ fields[k] = strtok(NULL, " \t\r\n"); if (!fields[k]) die("short body line"); }
            V s;
            s = from_text(fields[0], 1); memcpy(E(mass, i), s, ESZ); free(s);
            for (int k = 0; k < 3; k++){ s = from_text(fields[1 + k], 1); memcpy(E(X0, 3 * i + k), s, ESZ); free(s); }
            for (int k = 0; k < 3; k++){ s = from_text(fields[4 + k], 1); memcpy(E(V0, 3 * i + k), s, ESZ); free(s); }
            i++;
        }else die("unknown line '%s' in %s", tok, path);
    }
    fclose(f);
    if (ns && in_cur != (long)N) die("the last system has %ld bodies, not %zu", in_cur, N);
    if (i != NB) die("read %zu bodies, expected %zu", i, NB);
    if (!have_G) die("no G line in %s", path);
    if (!ns) snprintf(sys_names[0], sizeof sys_names[0], "%.31s", problem_name);
}

/* ------------------------------------------------------------------ */
/* State                                                                */
/* ------------------------------------------------------------------ */
static V x, v, a;                      /* the particles: position, velocity, acceleration (N3) */
static V x0, v0, a0, csx, csv, csa0, at; /* IAS15's step-start copies and compensations (N3) */
static V g[7], b[7], e[7], csb[7], er[7], br[7];
static V gcs;                            /* gravity_cs: always +0 for basic gravity */
static V T1, T2, T3, T4;                 /* scratch, NMAX */
/* per system, one scalar each (E elements) */
static V SDT, SDTLAST, STPLAIN, STHI, STLO, SPCE, SPCELAST, SDTNEW, SRATIO, SDTDONE;
/* per lane (N3): the owning system's dt, dt_done and step ratio */
static V DTB, DTDB, RB;
static V PE[7], PB[7];                   /* predict_next_step's inputs, gathered per system */
static int *active, *pc_active, *accepted, *pns_skip, *iters, *pc_max;
static long *done, *rejected, *max_exceeded;
static unsigned long long *pc_total;
static int cs_augmented = 0;
static int max_iter = 12;
static int trace_pc = 0;          /* --trace-pc N: print the corrector's error per iteration for N steps (system 0) */
static long steps_traced = 0;
static unsigned long long pc_lane_sum, pc_lane_max;   /* lane efficiency of the shared corrector loop */
static int adaptive;              /* epsilon > 0 */
/* --dt-file: a prescribed step sequence, shared by every system */
static V *dtlist; static long ndtlist;
static FILE *dt_out;

/* the per-lane state, as slices for masking */
#define NSLICE 52
static V slice_ptr[NSLICE];
static V *snap;                   /* E snapshot blocks of NSLICE * L elements */

/* pairs, in REBOUND's visiting order within each system: i = 1..N-1, j = 0..i-1 */
static size_t *pair_i, *pair_j;
static V pxi, pxj, pmi, pmj, pdx, pdy, pdz, pt1, pt2, ps, pr, pr3, ppf, ppfi, ppfj;
static V pcxi, pcyi, pczi, pcxj, pcyj, pczj;
static V addend;

static void alloc_state(void){
    x = valloc(N3); v = valloc(N3); a = valloc(N3);
    x0 = valloc(N3); v0 = valloc(N3); a0 = valloc(N3); csx = valloc(N3); csv = valloc(N3); csa0 = valloc(N3); at = valloc(N3);
    for (int m = 0; m < 7; m++){ g[m] = valloc(N3); b[m] = valloc(N3); e[m] = valloc(N3); csb[m] = valloc(N3); er[m] = valloc(N3); br[m] = valloc(N3); }
    gcs = valloc(N3);
    T1 = valloc(NMAX); T2 = valloc(NMAX); T3 = valloc(NMAX); T4 = valloc(NMAX);
    SDT = valloc(E); SDTLAST = valloc(E); STPLAIN = valloc(E); STHI = valloc(E); STLO = valloc(E);
    SPCE = valloc(E); SPCELAST = valloc(E); SDTNEW = valloc(E); SRATIO = valloc(E); SDTDONE = valloc(E);
    DTB = valloc(NMAX); DTDB = valloc(NMAX); RB = valloc(NMAX);
    for (int m = 0; m < 7; m++){ PE[m] = valloc(N3); PB[m] = valloc(N3); }
    active = calloc(E, sizeof *active); pc_active = calloc(E, sizeof *pc_active); accepted = calloc(E, sizeof *accepted);
    pns_skip = calloc(E, sizeof *pns_skip); iters = calloc(E, sizeof *iters); pc_max = calloc(E, sizeof *pc_max);
    done = calloc(E, sizeof *done); rejected = calloc(E, sizeof *rejected); max_exceeded = calloc(E, sizeof *max_exceeded);
    pc_total = calloc(E, sizeof *pc_total);
    if (!active || !pc_active || !accepted || !pns_skip || !iters || !pc_max || !done || !rejected || !max_exceeded || !pc_total) die("out of memory");
    { int k = 0;
      slice_ptr[k++] = x; slice_ptr[k++] = v; slice_ptr[k++] = a; slice_ptr[k++] = x0; slice_ptr[k++] = v0; slice_ptr[k++] = a0;
      slice_ptr[k++] = csx; slice_ptr[k++] = csv; slice_ptr[k++] = csa0; slice_ptr[k++] = at;
      for (int m = 0; m < 7; m++){ slice_ptr[k++] = g[m]; slice_ptr[k++] = b[m]; slice_ptr[k++] = e[m]; slice_ptr[k++] = csb[m]; slice_ptr[k++] = er[m]; slice_ptr[k++] = br[m]; }
      if (k != NSLICE) die("slice table"); }
    snap = calloc(E, sizeof *snap); if (!snap) die("out of memory");
    for (size_t s = 0; s < E; s++) snap[s] = valloc(NSLICE * L);
    pair_i = calloc(P ? P : 1, sizeof *pair_i); pair_j = calloc(P ? P : 1, sizeof *pair_j);
    size_t l = 0;
    for (size_t s = 0; s < E; s++)
        for (size_t i = 1; i < N; i++) for (size_t j = 0; j < i; j++){ pair_i[l] = s * N + i; pair_j[l] = s * N + j; l++; }
    pxi = valloc(P); pxj = valloc(P); pmi = valloc(P); pmj = valloc(P);
    pdx = valloc(P); pdy = valloc(P); pdz = valloc(P); pt1 = valloc(P); pt2 = valloc(P); ps = valloc(P);
    pr = valloc(P); pr3 = valloc(P); ppf = valloc(P); ppfi = valloc(P); ppfj = valloc(P);
    pcxi = valloc(P); pcyi = valloc(P); pczi = valloc(P); pcxj = valloc(P); pcyj = valloc(P); pczj = valloc(P);
    addend = valloc(N3);
    for (l = 0; l < P; l++){ memcpy(E(pmi, l), E(mass, pair_i[l]), ESZ); memcpy(E(pmj, l), E(mass, pair_j[l]), ESZ); }
}

/* Masking. A system that must not change while the shared vectors do
 * is snapshotted before and restored after: byte copies, so what comes
 * back is exactly what went in. */
static void save_system(size_t s){
    for (int k = 0; k < NSLICE; k++) memcpy(E(snap[s], (size_t)k * L), E(slice_ptr[k], L * s), L * ESZ);
}
static void restore_system(size_t s){
    for (int k = 0; k < NSLICE; k++) memcpy(E(slice_ptr[k], L * s), E(snap[s], (size_t)k * L), L * ESZ);
}
/* a per-system scalar into every lane of that system */
static void bcast_sys(V dst, const V src){
    for (size_t s = 0; s < E; s++) for (size_t k = 0; k < L; k++) memcpy(E(dst, L * s + k), E(src, s), ESZ);
}
/* copy one system's slice of a per-lane vector */
static void copy_sys(V dst, const V src, size_t s){ memcpy(E(dst, L * s), E(src, L * s), L * ESZ); }
static void zero_sys(V dst, size_t s){ memset(E(dst, L * s), 0, L * ESZ); }

/* ------------------------------------------------------------------ */
/* Gravity, as REBOUND's reb_gravity_basic_calculate_acceleration       */
/* ------------------------------------------------------------------ */
static void gravity(void){
    if (P){
        /* per pair, in the order of REBOUND's loop, vectorised over pairs */
        for (size_t l = 0; l < P; l++){ memcpy(E(pxi, l), E(x, 3 * pair_i[l]), ESZ); memcpy(E(pxj, l), E(x, 3 * pair_j[l]), ESZ); }
        vsub(pdx, pxi, pxj, P);                                   /* dx = x_i - x_j   (ghost box offset is +0: gb.x + x_i == x_i) */
        for (size_t l = 0; l < P; l++){ memcpy(E(pxi, l), E(x, 3 * pair_i[l] + 1), ESZ); memcpy(E(pxj, l), E(x, 3 * pair_j[l] + 1), ESZ); }
        vsub(pdy, pxi, pxj, P);
        for (size_t l = 0; l < P; l++){ memcpy(E(pxi, l), E(x, 3 * pair_i[l] + 2), ESZ); memcpy(E(pxj, l), E(x, 3 * pair_j[l] + 2), ESZ); }
        vsub(pdz, pxi, pxj, P);
        vmul(pt1, pdx, pdx, P);                                   /* dx*dx + dy*dy + dz*dz + softening2, left to right */
        vmul(pt2, pdy, pdy, P);
        vadd(ps, pt1, pt2, P);
        vmul(pt1, pdz, pdz, P);
        vadd(ps, ps, pt1, P);
        vadd(ps, ps, K0, P);                                      /* + softening2, which is 0*0 = +0 */
        vsqrt(pr, ps, P);                                         /* _r */
        vmul(pr3, pr, pr, P);                                     /* _r*_r*_r */
        vmul(pr3, pr3, pr, P);
        vdiv(ppf, G, pr3, P);                                     /* prefact = G/(_r*_r*_r) */
        vneg(pt1, ppf, P);                                        /* prefactj = -prefact*m_j */
        vmul(ppfj, pt1, pmj, P);
        vmul(ppfi, ppf, pmi, P);                                  /* prefacti = prefact*m_i */
        vmul(pcxi, ppfj, pdx, P); vmul(pcyi, ppfj, pdy, P); vmul(pczi, ppfj, pdz, P);   /* what particle i receives */
        vmul(pcxj, ppfi, pdx, P); vmul(pcyj, ppfi, pdy, P); vmul(pczj, ppfi, pdz, P);   /* what particle j receives */
    }
    /* ax = 0, then each particle's partners in ascending order, which is
     * the order REBOUND's (i, j<i) loop delivers them; every system's
     * particles at once, each with its own partners */
    vzero(a, N3);
    for (size_t t = 0; t + 1 < N; t++){
        for (size_t pg = 0; pg < NB; pg++){
            size_t s = pg / N, p = pg % N;
            size_t q = (t < p) ? t : t + 1;   /* the t-th partner of p, within the system */
            size_t l; V cx, cy, cz;
            if (q < p){ l = s * PS + p * (p - 1) / 2 + q; cx = pcxi; cy = pcyi; cz = pczi; }
            else      { l = s * PS + q * (q - 1) / 2 + p; cx = pcxj; cy = pcyj; cz = pczj; }
            memcpy(E(addend, 3 * pg),     E(cx, l), ESZ);
            memcpy(E(addend, 3 * pg + 1), E(cy, l), ESZ);
            memcpy(E(addend, 3 * pg + 2), E(cz, l), ESZ);
        }
        vadd(a, a, addend, N3);
    }
}

/* ------------------------------------------------------------------ */
/* Compensated summation: REBOUND's add_cs, or the exact 9.5 primitive  */
/* ------------------------------------------------------------------ */
static void add_cs(V p, V cs, const V inp, size_t n){
    /* y = inp - cs; t = p + y; cs = (t - p) - y; p = t */
    vsub(T3, inp, cs, n);
    if (!cs_augmented){
        vadd(T4, p, T3, n);
        vsub(T2, T4, p, n);
        vsub(cs, T2, T3, n);
        vcopy(p, T4, n);
    }else{
        /* (t, err) = augmentedAddition(p, y): t + err == p + y exactly,
         * so the compensation is exactly -err instead of Kahan's
         * estimate (t - p) - y, which is exact only when |p| >= |y|. */
        vaugadd(T4, T2, p, T3, n);
        vneg(cs, T2, n);
        vcopy(p, T4, n);
    }
}

/* ------------------------------------------------------------------ */
/* The predictor: positions at substep n from the b polynomial          */
/* ------------------------------------------------------------------ */
static V D1B[8], D2B[8];   /* fl(dt h_n) and its half, per lane, for the FMA form */

static void predict_positions(int n){
    V H = KH[n];
    if (arith_fma){
        /* t = b6; t = t*K_lvl + b_j down the levels; then + a0, then
         * fl(dt h/2) and v0, then * fl(dt h) - one rounding per level
         * instead of REBOUND's three */
        vcopy(T1, b[6], N3);
        for (int lvl = 0; lvl < 6; lvl++) vfma(T1, T1, KHF[n - 1][lvl], b[5 - lvl], N3);
        vfma(T1, T1, KHF[n - 1][6], a0, N3);
        vfma(T1, T1, D2B[n], v0, N3);
        vmul(T1, T1, D1B[n], N3);
        vsub(T1, T1, csx, N3);
        vadd(x, T1, x0, N3);
        return;
    }
    /* ((((((((b6*7*h/9 + b5)*3*h/4 + b4)*5*h/7 + b3)*2*h/3 + b2)*3*h/5 + b1)*h/2 + b0)*h/3 + a0)*dt*h/2 + v0)*dt*h,
     * each level as ((T*num)*h)/den with num/den = (j+1)/(j+3); the
     * levels REBOUND writes reduced (3/4, 2/3, 1/2) differ only by a
     * power of two, which rounds identically. The binary64 gate is what
     * proves that sentence. */
    /* the b6 level: (b6*7)*h/9; the numerator 1 of the b0 level is a
     * multiply by one, exact, kept so every level is the same code */
    vmul(T1, b[6], KNUM[0], N3);
    for (int lvl = 0; lvl < 7; lvl++){
        if (lvl > 0) vmul(T1, T1, KNUM[lvl], N3);
        vmul(T1, T1, H, N3); vdiv(T1, T1, KDEN[lvl], N3);
        vadd(T1, T1, lvl < 6 ? b[5 - lvl] : a0, N3);
    }
    vmul(T1, T1, DTB, N3); vmul(T1, T1, H, N3); vdiv(T1, T1, K2, N3); vadd(T1, T1, v0, N3);   /* *dt*h/2 + v0 */
    vmul(T1, T1, DTB, N3); vmul(T1, T1, H, N3);                                              /* *dt*h */
    vsub(T1, T1, csx, N3);          /* xk = -csx + (...) */
    vadd(x, T1, x0, N3);            /* particle.x = xk + x0 */
}

/* ------------------------------------------------------------------ */
/* The corrector at substep n: improve g and b                          */
/* ------------------------------------------------------------------ */
static void pc_error(const V tmp);

static void correct(int n){
    /* gk = at; gk_cs = gravity_cs; add_cs(gk, gk_cs, -a0); add_cs(gk, gk_cs, csa0) */
    static V gk, gk_cs, tmp, neg_a0, gnew, told, y, t, u, prod;
    if (!gk){ gk = valloc(N3); gk_cs = valloc(N3); tmp = valloc(N3); neg_a0 = valloc(N3); gnew = valloc(N3); told = valloc(N3);
              y = valloc(N3); t = valloc(N3); u = valloc(N3); prod = valloc(N3); }
    vcopy(gk, at, N3);
    vcopy(gk_cs, gcs, N3);
    vneg(neg_a0, a0, N3);
    /* REBOUND applies add_cs twice here, Kahan's form regardless of --cs:
     * the compensation of the FORCE (gravity_cs) is a separate mechanism
     * from the compensation of the b coefficients, and with basic gravity
     * both inputs are +0. Written out so the roundings are its. */
    vsub(y, neg_a0, gk_cs, N3); vadd(t, gk, y, N3); vsub(u, t, gk, N3); vsub(gk_cs, u, y, N3); vcopy(gk, t, N3);
    vsub(y, csa0, gk_cs, N3);   vadd(t, gk, y, N3); vsub(u, t, gk, N3); vsub(gk_cs, u, y, N3); vcopy(gk, t, N3);
    /* g[n-1] = (((gk/rr[base] - g0)/rr[base+1] - g1)/... - g[n-2])/rr[base+n-1] */
    size_t base = (size_t)n * (n - 1) / 2;
    vcopy(told, g[n - 1], N3);
    if (arith_fma) vmul(gnew, gk, KRINV[base], N3); else vdiv(gnew, gk, KRR[base], N3);
    for (int m = 1; m < n; m++){
        vsub(gnew, gnew, g[m - 1], N3);
        if (arith_fma) vmul(gnew, gnew, KRINV[base + m], N3); else vdiv(gnew, gnew, KRR[base + m], N3);
    }
    vcopy(g[n - 1], gnew, N3);
    vsub(tmp, gnew, told, N3);                 /* tmp = g_new - g_old */
    size_t cbase = (size_t)(n - 1) * (n - 2) / 2;
    for (int m = 0; m < n - 1; m++){ vmul(prod, tmp, KC[cbase + m], N3); add_cs(b[m], csb[m], prod, N3); }
    add_cs(b[n - 1], csb[n - 1], tmp, N3);
    if (n == 7) pc_error(tmp);
}

/* predictor_corrector_error = max|tmp| / max|at| over normal values
 * (GLOBAL/PRS23 modes), per system over its own coordinates */
static void pc_error(const V tmp){
    static V maxak, maxb6, aabs, tabs; static uint8_t *ca, *ct;
    if (!maxak){ maxak = valloc(1); maxb6 = valloc(1); aabs = valloc(N3); tabs = valloc(N3); ca = malloc(N3); ct = malloc(N3); }
    vabs(aabs, at, N3); vabs(tabs, tmp, N3);
    vclass(ca, aabs, N3); vclass(ct, tabs, N3);
    for (size_t s = 0; s < E; s++){
        if (!pc_active[s]) continue;
        vzero(maxak, 1); vzero(maxb6, 1);
        for (size_t k = L * s; k < L * (s + 1); k++){
            if (is_normal_class(ca[k]) && s_lt(maxak, E(aabs, k))) vcopy(maxak, E(aabs, k), 1);
            if (is_normal_class(ct[k]) && s_lt(maxb6, E(tabs, k))) vcopy(maxb6, E(tabs, k), 1);
        }
        vdiv(E(SPCE, s), maxb6, maxak, 1);
    }
}

/* ------------------------------------------------------------------ */
/* The same two pieces as sequencer programs (--engine program)         */
/* ------------------------------------------------------------------ */
static cft_program *prog_predict, *prog_predict_ens, *prog_correct[8];
static V prog_bank, prog_sin, prog_sout, prog_dep;
static int use_ens_predict;   /* E > 1 with per-system steps: dt rides in the scratch block, not the bank */

static void *read_file(const char *path, size_t *len){
    FILE *f = fopen(path, "rb");
    if (!f) die("cannot open program image %s", path);
    fseek(f, 0, SEEK_END); long n = ftell(f); fseek(f, 0, SEEK_SET);
    void *buf = malloc(n > 0 ? (size_t)n : 1);
    if (fread(buf, 1, (size_t)n, f) != (size_t)n) die("short read on %s", path);
    fclose(f); *len = (size_t)n; return buf;
}

static cft_program *load_one(const char *stem){
    char path[512]; size_t len; void *img; cft_program *p;
    snprintf(path, sizeof path, "%s/%s-%s.cftp", programs_dir, stem, cft_format_name(F));
    img = read_file(path, &len);
    cft_status st = cft_program_load(dev, img, len, &p);
    if (st != CFT_OK) die("cft_program_load(%s): %s (%s)", path, cft_strerror(st), cft_last_error());
    free(img);
    return p;
}

static void load_programs(void){
    prog_predict = load_one("predict");
    if (use_ens_predict) prog_predict_ens = load_one("predict-ens");
    for (int n = 1; n < 8; n++){ char stem[16]; snprintf(stem, sizeof stem, "correct%d", n); prog_correct[n] = load_one(stem); }
    prog_bank = valloc(16); prog_sin = valloc(N3 * 22); prog_sout = valloc(N3 * 22); prog_dep = valloc(N3);
}

static void run_program(cft_program *prog, const V a_, const V bb, const V c, size_t nbank,
                        size_t nsin, size_t nsout, const char *what){
    cft_run_args ra; memset(&ra, 0, sizeof ra);
    ra.struct_size = sizeof ra;
    ra.a = a_; ra.b = bb; ra.c = c; ra.n = N3;
    ra.bank = prog_bank; ra.bank_bytes = nbank * ESZ;
    ra.scratch_in = nsin ? prog_sin : NULL; ra.scratch_in_bytes = nsin * N3 * ESZ;
    ra.scratch_out = nsout ? prog_sout : NULL; ra.scratch_out_bytes = nsout * N3 * ESZ;
    ra.deposits = prog_dep; ra.counts = NULL;
    uint32_t fl = 0, bus = 0; ra.flags_out = &fl; ra.bus_out = &bus;
    cft_status st = cft_program_run_ex(prog, &ra);
    if (st == CFT_OK && (bus & CFT_STATUS_DEPOSIT_OVERFLOW)) die("%s: deposit overflow", what);
    note(st, fl, what);
}

static void predict_positions_program(int n){
    for (int l = 0; l < 7; l++) memcpy(E(prog_bank, l), KHF[n - 1][l], ESZ);
    if (!use_ens_predict){
        /* bank: K0..K6, D2, D1 - one dt for the whole run, lane 0's */
        memcpy(E(prog_bank, 7), E(D2B[n], 0), ESZ);
        memcpy(E(prog_bank, 8), E(D1B[n], 0), ESZ);
        /* scratch in, lane-major: b0..b6, csx */
        for (size_t k = 0; k < N3; k++){
            for (int m = 0; m < 7; m++) memcpy(E(prog_sin, 8 * k + m), E(b[m], k), ESZ);
            memcpy(E(prog_sin, 8 * k + 7), E(csx, k), ESZ);
        }
        run_program(prog_predict, x0, v0, a0, 9, 8, 0, "predict program");
    }else{
        /* bank: K0..K6; scratch in: b0..b6, csx, then the lane's own D2, D1 */
        for (size_t k = 0; k < N3; k++){
            for (int m = 0; m < 7; m++) memcpy(E(prog_sin, 10 * k + m), E(b[m], k), ESZ);
            memcpy(E(prog_sin, 10 * k + 7), E(csx, k), ESZ);
            memcpy(E(prog_sin, 10 * k + 8), E(D2B[n], k), ESZ);
            memcpy(E(prog_sin, 10 * k + 9), E(D1B[n], k), ESZ);
        }
        run_program(prog_predict_ens, x0, v0, a0, 7, 10, 0, "predict-ens program");
    }
    vcopy(x, prog_dep, N3);
}

static void correct_program(int n){
    size_t base = (size_t)n * (n - 1) / 2, cbase = (size_t)(n - 1) * (n - 2) / 2;
    for (int m = 0; m < n; m++) memcpy(E(prog_bank, m), KRINV[base + m], ESZ);
    for (int m = 0; m < n - 1; m++) memcpy(E(prog_bank, n + m), KC[cbase + m], ESZ);
    for (size_t k = 0; k < N3; k++){
        memcpy(E(prog_sin, 22 * k), E(a0, k), ESZ);
        for (int m = 0; m < 7; m++){
            memcpy(E(prog_sin, 22 * k + 1 + m), E(g[m], k), ESZ);
            memcpy(E(prog_sin, 22 * k + 8 + m), E(b[m], k), ESZ);
            memcpy(E(prog_sin, 22 * k + 15 + m), E(csb[m], k), ESZ);
        }
    }
    run_program(prog_correct[n], at, NULL, NULL, (size_t)(2 * n - 1), 22, 22, "correct program");
    for (size_t k = 0; k < N3; k++){
        for (int m = 0; m < 7; m++){
            memcpy(E(g[m], k),   E(prog_sout, 22 * k + 1 + m), ESZ);
            memcpy(E(b[m], k),   E(prog_sout, 22 * k + 8 + m), ESZ);
            memcpy(E(csb[m], k), E(prog_sout, 22 * k + 15 + m), ESZ);
        }
    }
    if (n == 7) pc_error(prog_dep);
}

/* ------------------------------------------------------------------ */
/* predict_next_step, verbatim in structure, over every lane at once   */
/* ------------------------------------------------------------------ */
/* rb holds each lane's own step ratio; the powers q2..q7 are the same
 * six multiplications REBOUND does on its scalar, issued per lane. The
 * ratio > 20 branch (zero the prediction) is applied per system by the
 * caller after the fact - it is unreachable in REBOUND's own flow. */
static void predict_next_step_vec(const V rb, V *_e, V *_b, V *eo, V *bo){
    static V q[8], be[7], s, u;
    if (!s){ for (int i = 1; i < 8; i++) q[i] = valloc(N3); for (int m = 0; m < 7; m++) be[m] = valloc(N3); s = valloc(N3); u = valloc(N3); }
    vcopy(q[1], rb, N3);
    vmul(q[2], q[1], q[1], N3); vmul(q[3], q[1], q[2], N3); vmul(q[4], q[2], q[2], N3);
    vmul(q[5], q[2], q[3], N3); vmul(q[6], q[3], q[3], N3); vmul(q[7], q[3], q[4], N3);
    for (int m = 0; m < 7; m++) vsub(be[m], _b[m], _e[m], N3);
    for (int m = 0; m < 7; m++){
        /* e_m = q_{m+1} * (b6*C(7,m+1) + b5*C(6,m+1) + ... + b_{m+1}*C(m+2,m+1) + b_m), left to right */
        int first = 1;
        for (int j = 6; j > m; j--){
            vmul(u, _b[j], KBIN[m][j], N3);
            if (first){ vcopy(s, u, N3); first = 0; } else vadd(s, s, u, N3);
        }
        if (first) vcopy(s, _b[m], N3); else vadd(s, s, _b[m], N3);
        vmul(eo[m], q[m + 1], s, N3);
    }
    for (int m = 0; m < 7; m++) vadd(bo[m], eo[m], be[m], N3);
}

/* ------------------------------------------------------------------ */
/* sqrt7 and the PRS23 step criterion                                   */
/* ------------------------------------------------------------------ */
static void sqrt7(V out, const V ain){
    V aa = valloc(1), scale = valloc(1), xx = valloc(1), x6 = valloc(1), q = valloc(1), d = valloc(1);
    vcopy(aa, ain, 1); vcopy(scale, K1, 1);
    while (s_lt(aa, K1EM7) && s_isnormal(aa)){ vmul(scale, scale, K0_1, 1); vmul(aa, aa, K1E7, 1); }
    while (s_lt(K1E2, aa) && s_isnormal(aa)){ vmul(scale, scale, K10, 1); vmul(aa, aa, K1EM7, 1); }
    vcopy(xx, K1, 1);
    for (int k = 0; k < 20; k++){
        vmul(x6, xx, xx, 1); vmul(x6, x6, xx, 1); vmul(x6, x6, xx, 1); vmul(x6, x6, xx, 1); vmul(x6, x6, xx, 1);
        vdiv(q, aa, x6, 1); vsub(d, q, xx, 1); vdiv(d, d, K7, 1); vadd(xx, xx, d, 1);
    }
    vmul(out, xx, scale, 1);
    free(aa); free(scale); free(xx); free(x6); free(q); free(d);
}

static V EPS, EPS5040;

/* For every active system: accepted[s] and SDTNEW[s], from SDTDONE[s].
 * The per-particle timescale is REBOUND's, particle by particle; the
 * minimum is taken within each system. */
static void choose_timestep(void){
    static V sq, tmp, y[6], a0i, ts2, mints2, num, den, s7; static uint8_t *cls;
    if (!sq){ sq = valloc(N3); tmp = valloc(N3); for (int i = 0; i < 6; i++) y[i] = valloc(NB); a0i = valloc(NB); ts2 = valloc(NB); mints2 = valloc(E); num = valloc(NB); den = valloc(NB); cls = malloc(NB);
              /* sqrt7(epsilon*5040) depends on nothing that changes: once, the same bits every step */
              s7 = valloc(1); sqrt7(s7, EPS5040); }
    /* per component: a0^2; (a0+b0+...+b6)^2; (sum (m+1) b_m)^2; (sum m(m+1) b_m)^2; (sum (m-1)m(m+1) b_m)^2 */
    V sums[5];
    static V S[5]; if (!S[0]) for (int i = 0; i < 5; i++) S[i] = valloc(N3);
    vmul(S[0], a0, a0, N3);
    vadd(tmp, a0, b[0], N3); for (int m = 1; m < 7; m++) vadd(tmp, tmp, b[m], N3);
    vmul(S[1], tmp, tmp, N3);
    vcopy(tmp, b[0], N3); for (int m = 1; m < 7; m++){ vmul(sq, b[m], KY3[m], N3); vadd(tmp, tmp, sq, N3); }
    vmul(S[2], tmp, tmp, N3);
    vmul(tmp, b[1], KY4[1], N3); for (int m = 2; m < 7; m++){ vmul(sq, b[m], KY4[m], N3); vadd(tmp, tmp, sq, N3); }
    vmul(S[3], tmp, tmp, N3);
    vmul(tmp, b[2], KY5[2], N3); for (int m = 3; m < 7; m++){ vmul(sq, b[m], KY5[m], N3); vadd(tmp, tmp, sq, N3); }
    vmul(S[4], tmp, tmp, N3);
    /* per particle: sum of the three components, starting from 0, in order */
    for (int i = 0; i < 5; i++){
        sums[i] = (i == 0) ? a0i : y[i];
        vzero(sums[i], NB);
        for (int c = 0; c < 3; c++){
            for (size_t p = 0; p < NB; p++) memcpy(E(tmp, p), E(S[i], 3 * p + c), ESZ);
            vadd(sums[i], sums[i], tmp, NB);
        }
    }
    /* timescale2 = 2*y2/(y3 + sqrt(y4*y2)), only for particles whose a0^2
     * is normal - REBOUND skips the others before dividing, and a 0/0
     * here would raise a certificate flag for a particle it never looks at */
    vclass(cls, a0i, NB);
    for (size_t s = 0; s < E; s++) vcopy(E(mints2, s), KINF, 1);
    for (size_t p = 0; p < NB; p++){
        size_t s = p / N;
        if (!active[s]) continue;
        if (!is_normal_class(cls[p])) continue;    /* no acceleration, or not finite: skip */
        vmul(E(num, p), K2, E(y[1], p), 1);
        vmul(E(tmp, p), E(y[3], p), E(y[1], p), 1); vsqrt(E(tmp, p), E(tmp, p), 1); vadd(E(den, p), E(y[2], p), E(tmp, p), 1);
        vdiv(E(ts2, p), E(num, p), E(den, p), 1);
        if (s_isnormal(E(ts2, p)) && s_lt(E(ts2, p), E(mints2, s))) vcopy(E(mints2, s), E(ts2, p), 1);
    }
    for (size_t s = 0; s < E; s++){
        if (!active[s]) continue;
        V dt_done = E(SDTDONE, s), dt_new = E(SDTNEW, s), mi = E(mints2, s);
        if (s_isnormal(mi)){
            V r = valloc(1);
            vsqrt(r, mi, 1); vmul(r, r, dt_done, 1);
            vmul(dt_new, r, s7, 1);
            free(r);
        }else{
            vdiv(dt_new, dt_done, K0_25, 1);
        }
        /* min_dt is 0: fabs(dt_new) < 0 is never true, but issue it */
        { V ad = valloc(1); vabs(ad, dt_new, 1); (void)s_lt(ad, K0); free(ad); }
        /* fabs(dt_new/dt_done) < safety_factor -> reject */
        { V rr = valloc(1), ar = valloc(1); int reject, larger;
          vdiv(rr, dt_new, dt_done, 1); vabs(ar, rr, 1);
          reject = s_lt(ar, K0_25);
          larger = s_lt(K1, ar);
          if (!reject && larger){
              V inv = valloc(1); vdiv(inv, K1, K0_25, 1);    /* 1./safety_factor */
              if (s_lt(inv, rr)) vdiv(dt_new, dt_done, K0_25, 1);
              free(inv);
          }
          free(rr); free(ar);
          accepted[s] = !reject; }
    }
}

/* ------------------------------------------------------------------ */
/* One attempt at a step, for every active system at once:             */
/* REBOUND's reb_integrator_ias15_step_try                              */
/* ------------------------------------------------------------------ */
static void step_attempt(void){
    for (size_t s = 0; s < E; s++) if (!active[s]) save_system(s);
    gravity();                                   /* reb_simulation_update_acceleration */
    vcopy(x0, x, N3); vcopy(v0, v, N3); vcopy(a0, a, N3);
    vzero(csa0, N3);                             /* basic gravity: gravity_cs is csa0, always 0 */
    for (int m = 0; m < 7; m++) vzero(csb[m], N3);
    /* g from b */
    for (int m = 0; m < 6; m++){
        int first = 1;
        for (int j = 6; j > m; j--){
            vmul(T2, b[j], KD[(size_t)j * (j - 1) / 2 + m], N3);
            if (first){ vcopy(T1, T2, N3); first = 0; } else vadd(T1, T1, T2, N3);
        }
        vadd(g[m], T1, b[m], N3);
    }
    vcopy(g[6], b[6], N3);
    if (arith_fma){
        /* fl(dt h_n) and its exact half, each lane with its system's dt */
        for (int n = 1; n < 8; n++){
            if (!D1B[n]){ D1B[n] = valloc(NMAX); D2B[n] = valloc(NMAX); }
            vmul(D1B[n], DTB, KH[n], N3); vmul(D2B[n], D1B[n], KHALF, N3);
        }
    }

    for (size_t s = 0; s < E; s++){
        vcopy(E(SPCE, s), K1E300, 1);
        vcopy(E(SPCELAST, s), K2, 1);
        iters[s] = 0;
        pc_active[s] = active[s];
    }
    int npass = 0;
    while (1){
        int any = 0;
        for (size_t s = 0; s < E; s++){
            if (!pc_active[s]) continue;
            V pce = E(SPCE, s), pce_last = E(SPCELAST, s);
            if (trace_pc && s == 0 && iters[s] > 0){
                char buf[160]; size_t len = 0;
                cft_to_decimal_char(dev, F, CFT_RNE, pce, 6, buf, sizeof buf, &len, NULL);
                fprintf(stderr, "  step %ld iteration %d: predictor_corrector_error %s\n", steps_traced, iters[s], buf);
            }
            int stop = 0;
            if (s_lt(pce, TOL)) stop = 1;
            else if (iters[s] > 2 && s_le(pce_last, pce)) stop = 1;
            else if (iters[s] >= max_iter){ max_exceeded[s]++; stop = 1; }
            if (stop){ pc_active[s] = 0; save_system(s); continue; }
            vcopy(pce_last, pce, 1);
            vzero(pce, 1);
            iters[s]++;
            any = 1;
        }
        if (!any) break;
        npass++;
        for (int n = 1; n < 8; n++){
            if (engine_program) predict_positions_program(n); else predict_positions(n);
            gravity();
            vcopy(at, a, N3);
            if (engine_program) correct_program(n); else correct(n);
        }
        /* a system that left the corrector keeps the state it left with */
        for (size_t s = 0; s < E; s++) if (active[s] && !pc_active[s]) restore_system(s);
    }
    if (trace_pc){ steps_traced++; if (steps_traced >= trace_pc) trace_pc = 0; }
    for (size_t s = 0; s < E; s++){
        if (!active[s]) continue;
        pc_total[s] += (unsigned long long)iters[s];
        if (iters[s] > pc_max[s]) pc_max[s] = iters[s];
        pc_lane_sum += (unsigned long long)iters[s];
        pc_lane_max += (unsigned long long)npass;
    }

    for (size_t s = 0; s < E; s++) vcopy(E(SDTDONE, s), E(SDT, s), 1);
    if (adaptive){
        choose_timestep();
        for (size_t s = 0; s < E; s++){
            if (!active[s]) continue;
            if (!accepted[s]){
                /* reset particles to the step start */
                copy_sys(x, x0, s); copy_sys(v, v0, s); copy_sys(a, a0, s);
                rejected[s]++;
            }
            vcopy(E(SDT, s), E(SDTNEW, s), 1);
        }
    }else{
        for (size_t s = 0; s < E; s++) accepted[s] = active[s];
    }
    /* positions and velocities at the end of the step, every lane; a
     * rejected system is restored afterwards, which is REBOUND's "skip" */
    for (size_t s = 0; s < E; s++) if (active[s] && !accepted[s]) save_system(s);
    bcast_sys(DTDB, SDTDONE);
    for (int m = 6; m >= 0; m--){
        if (arith_fma) vmul(T1, b[m], KPOSR[m], N3); else vdiv(T1, b[m], KPOS[m], N3);
        vmul(T1, T1, DTDB, N3); vmul(T1, T1, DTDB, N3); add_cs(x0, csx, T1, N3);
    }
    if (arith_fma) vmul(T1, a0, KPOSR[7], N3); else vdiv(T1, a0, KPOS[7], N3);
    vmul(T1, T1, DTDB, N3); vmul(T1, T1, DTDB, N3); add_cs(x0, csx, T1, N3);
    vmul(T1, v0, DTDB, N3); add_cs(x0, csx, T1, N3);
    for (int m = 6; m >= 0; m--){
        if (arith_fma) vmul(T1, b[m], KVELR[m], N3); else vdiv(T1, b[m], KVEL[m], N3);
        vmul(T1, T1, DTDB, N3); add_cs(v0, csv, T1, N3);
    }
    vmul(T1, a0, DTDB, N3); add_cs(v0, csv, T1, N3);
    for (size_t s = 0; s < E; s++) if (active[s] && !accepted[s]) restore_system(s);

    for (size_t s = 0; s < E; s++){
        if (!(active[s] && accepted[s])) continue;
        V dt_done = E(SDTDONE, s);
        /* t += dt_done: the plain sum REBOUND keeps, and an exact one beside it */
        vadd(E(STPLAIN, s), E(STPLAIN, s), dt_done, 1);
        { V r = valloc(1), err = valloc(1); vaugadd(r, err, E(STHI, s), dt_done, 1); vadd(E(STLO, s), E(STLO, s), err, 1); vcopy(E(STHI, s), r, 1); free(r); free(err); }
        vcopy(E(SDTLAST, s), dt_done, 1);
        if (dt_out) { char buf[160]; size_t len = 0; cft_to_hex_char(dev, F, dt_done, buf, sizeof buf, &len); fprintf(dt_out, "%s\n", buf); }
        done[s]++;
        if (dtlist){
            if (done[s] < ndtlist) vcopy(E(SDT, s), dtlist[done[s]], 1);
        }
    }
    vcopy(x, x0, N3); vcopy(v, v0, N3);
    for (size_t s = 0; s < E; s++){
        if (!(active[s] && accepted[s])) continue;
        for (int m = 0; m < 7; m++){ copy_sys(er[m], e[m], s); copy_sys(br[m], b[m], s); }
    }
    /* the prediction of the next step's e and b: from e, b with ratio
     * dt/dt_done after an accepted step; from er, br with ratio
     * dt/dt_last_done after a rejected one (not at all if it was the
     * first step); every active system's lanes in one pass */
    int any_pns = 0;
    for (size_t s = 0; s < E; s++){
        pns_skip[s] = 1;
        if (!active[s]) continue;
        if (accepted[s]){
            vdiv(E(SRATIO, s), E(SDT, s), E(SDTDONE, s), 1);
            for (int m = 0; m < 7; m++){ copy_sys(PE[m], e[m], s); copy_sys(PB[m], b[m], s); }
            pns_skip[s] = 0;
        }else if (!s_eq(E(SDTLAST, s), K0)){
            vdiv(E(SRATIO, s), E(SDT, s), E(SDTLAST, s), 1);
            for (int m = 0; m < 7; m++){ copy_sys(PE[m], er[m], s); copy_sys(PB[m], br[m], s); }
            pns_skip[s] = 0;
        }
        if (!pns_skip[s]) any_pns = 1;
    }
    if (any_pns){
        bcast_sys(RB, SRATIO);
        predict_next_step_vec(RB, PE, PB, e, b);
        for (size_t s = 0; s < E; s++){
            if (pns_skip[s]) continue;
            if (s_lt(K20, E(SRATIO, s)))           /* ratio > 20: do not predict */
                for (int m = 0; m < 7; m++){ zero_sys(e[m], s); zero_sys(b[m], s); }
        }
        for (size_t s = 0; s < E; s++) if (active[s] && pns_skip[s]) restore_system(s);
    }
    for (size_t s = 0; s < E; s++) if (!active[s]) restore_system(s);
    bcast_sys(DTB, SDT);
}

/* ------------------------------------------------------------------ */
/* Energy, as REBOUND's reb_simulation_energy, every system at once     */
/* ------------------------------------------------------------------ */
static void energy_all(V out){
    static V bvx, bvy, bvz, bt, bu, ekin, epot, half, halfb, tt, uu;
    if (!bvx){ bvx = valloc(NB); bvy = valloc(NB); bvz = valloc(NB); bt = valloc(NB); bu = valloc(NB);
               ekin = valloc(E); epot = valloc(E); half = valloc(1); halfb = valloc(NB); tt = valloc(P); uu = valloc(P); }
    vdiv(half, K1, K2, 1); vbcast(halfb, half, NB);
    /* 0.5 * m * (vx*vx + vy*vy + vz*vz), every body at once */
    for (size_t i = 0; i < NB; i++){ memcpy(E(bvx, i), E(v, 3 * i), ESZ); memcpy(E(bvy, i), E(v, 3 * i + 1), ESZ); memcpy(E(bvz, i), E(v, 3 * i + 2), ESZ); }
    vmul(bt, bvx, bvx, NB); vmul(bu, bvy, bvy, NB); vadd(bt, bt, bu, NB);
    vmul(bu, bvz, bvz, NB); vadd(bt, bt, bu, NB);
    vmul(bu, halfb, mass, NB); vmul(bu, bu, bt, NB);
    /* - G m_i m_j / r for every pair at once (the pair list is (i > j);
     * the energy loop's (i < j) gives the same dx up to sign, and the
     * square is the same bits) */
    if (P){
        for (size_t l = 0; l < P; l++){ memcpy(E(pxi, l), E(x, 3 * pair_i[l]), ESZ); memcpy(E(pxj, l), E(x, 3 * pair_j[l]), ESZ); }
        vsub(pdx, pxi, pxj, P);
        for (size_t l = 0; l < P; l++){ memcpy(E(pxi, l), E(x, 3 * pair_i[l] + 1), ESZ); memcpy(E(pxj, l), E(x, 3 * pair_j[l] + 1), ESZ); }
        vsub(pdy, pxi, pxj, P);
        for (size_t l = 0; l < P; l++){ memcpy(E(pxi, l), E(x, 3 * pair_i[l] + 2), ESZ); memcpy(E(pxj, l), E(x, 3 * pair_j[l] + 2), ESZ); }
        vsub(pdz, pxi, pxj, P);
        vmul(tt, pdx, pdx, P); vmul(uu, pdy, pdy, P); vadd(tt, tt, uu, P); vmul(uu, pdz, pdz, P); vadd(tt, tt, uu, P); vsqrt(tt, tt, P);
        /* G*m_j*m_i with j the larger index, as reb_simulation_energy
         * multiplies it; the list's pair_i IS the larger index */
        vmul(uu, G, pmi, P); vmul(uu, uu, pmj, P); vdiv(uu, uu, tt, P);
    }
    /* the sums, in REBOUND's order within each system */
    for (size_t s = 0; s < E; s++){
        vzero(E(ekin, s), 1); vzero(E(epot, s), 1);
        for (size_t i = 0; i < N; i++) vadd(E(ekin, s), E(ekin, s), E(bu, s * N + i), 1);
        for (size_t i = 0; i < N; i++) for (size_t j = i + 1; j < N; j++)
            vsub(E(epot, s), E(epot, s), E(uu, s * PS + j * (j - 1) / 2 + i), 1);   /* pair (j > i) of the list */
        vadd(E(out, s), E(ekin, s), E(epot, s), 1);
        vadd(E(out, s), E(out, s), K0, 1);   /* + energy_offset, which is 0 */
    }
}

/* ------------------------------------------------------------------ */
/* Records                                                              */
/* ------------------------------------------------------------------ */
static void puthex(const V a_){
    char buf[160]; size_t len = 0;
    cft_status st = cft_to_hex_char(dev, F, a_, buf, sizeof buf, &len);
    if (st != CFT_OK) die("to_hex: %s", cft_strerror(st));
    printf(" %s", buf);
}

static void dump_constants(void){
    printf("# constants at %s (hex = exact bits)\n", cft_format_name(F));
    for (int i = 0; i < 8; i++){ printf("h[%d]", i); puthex(KH[i]); printf("\n"); }
    for (int i = 0; i < 28; i++){ printf("rr[%d]", i); puthex(KRR[i]); printf("\n"); }
    for (int i = 0; i < 21; i++){ printf("c[%d]", i); puthex(KC[i]); printf("\n"); }
    for (int i = 0; i < 21; i++){ printf("d[%d]", i); puthex(KD[i]); printf("\n"); }
    printf("pc_tol"); puthex(TOL); printf("\n");
}

static void read_dt_file(const char *path, long steps){
    FILE *f = fopen(path, "r");
    if (!f){ perror(path); exit(2); }
    char line[256]; long cap = 1024, n = 0;
    dtlist = malloc(cap * sizeof *dtlist);
    while (fgets(line, sizeof line, f)){
        if (line[0] == '#' || line[0] == '\n') continue;
        char *tok = strtok(line, " \t\r\n"); if (!tok) continue;
        if (n == cap){ cap *= 2; dtlist = realloc(dtlist, cap * sizeof *dtlist); if (!dtlist) die("out of memory"); }
        dtlist[n++] = from_text(tok, 1);
    }
    fclose(f);
    if (n < steps) die("--dt-file %s has %ld steps, the run needs %ld", path, n, steps);
    ndtlist = n;
}

int main(int argc, char **argv){
    const char *problem = NULL, *artifact = NULL, *fmtname = "fp64", *dt_file = NULL, *dt_out_path = NULL;
    const char *dt_txt = "0.01", *eps_txt = "1e-9";
    long steps = 1000, sample = 100;
    int tol_shift = -1, quiet = 0, do_dump = 0;
    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "--format") && i + 1 < argc) fmtname = argv[++i];
        else if (!strcmp(argv[i], "--problem") && i + 1 < argc) problem = argv[++i];
        else if (!strcmp(argv[i], "--dt") && i + 1 < argc) dt_txt = argv[++i];
        else if (!strcmp(argv[i], "--epsilon") && i + 1 < argc) eps_txt = argv[++i];
        else if (!strcmp(argv[i], "--steps") && i + 1 < argc) steps = atol(argv[++i]);
        else if (!strcmp(argv[i], "--sample") && i + 1 < argc) sample = atol(argv[++i]);
        else if (!strcmp(argv[i], "--cs") && i + 1 < argc){ const char *m = argv[++i]; if (!strcmp(m, "kahan")) cs_augmented = 0; else if (!strcmp(m, "augmented")) cs_augmented = 1; else die("--cs kahan|augmented"); }
        else if (!strcmp(argv[i], "--pc-tol-shift") && i + 1 < argc) tol_shift = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--max-iter") && i + 1 < argc) max_iter = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--trace-pc") && i + 1 < argc) trace_pc = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--arith") && i + 1 < argc){ const char *m = argv[++i]; if (!strcmp(m, "rebound")) arith_fma = 0; else if (!strcmp(m, "fma")) arith_fma = 1; else die("--arith rebound|fma"); }
        else if (!strcmp(argv[i], "--engine") && i + 1 < argc){ const char *m = argv[++i]; if (!strcmp(m, "loop")) engine_program = 0; else if (!strcmp(m, "program")) engine_program = 1; else die("--engine loop|program"); }
        else if (!strcmp(argv[i], "--programs") && i + 1 < argc) programs_dir = argv[++i];
        else if (!strcmp(argv[i], "--artifact") && i + 1 < argc) artifact = argv[++i];
        else if (!strcmp(argv[i], "--member") && i + 1 < argc) member = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--dt-file") && i + 1 < argc) dt_file = argv[++i];
        else if (!strcmp(argv[i], "--dt-out") && i + 1 < argc) dt_out_path = argv[++i];
        else if (!strcmp(argv[i], "--no-flag-abort")) flag_abort = 0;
        else if (!strcmp(argv[i], "--dump-constants")) do_dump = 1;
        else if (!strcmp(argv[i], "--quiet")) quiet = 1;
        else die("unknown argument %s", argv[i]);
    }
    if (!strcmp(fmtname, "fp64")){ F = CFT_FP64; FI = 0; if (tol_shift < 0) tol_shift = 0; }
    else if (!strcmp(fmtname, "fp128")){ F = CFT_FP128; FI = 1; if (tol_shift < 0) tol_shift = 113 - 53; }
    else if (!strcmp(fmtname, "fp256")){ F = CFT_FP256; FI = 2; if (tol_shift < 0) tol_shift = 237 - 53; }
    else die("--format fp64|fp128|fp256");
    ESZ = cft_format_size(F);
    cft_status st = cft_open(artifact, 0, &dev);
    if (st != CFT_OK) die("cft_open(%s): %s (%s)", artifact ? artifact : "software", cft_strerror(st), cft_last_error());
    if (!problem){
        if (do_dump){ NMAX = 8; make_constants(tol_shift, quiet); dump_constants(); return 0; }
        die("usage: ias15_cft --format F --problem FILE [--dt DT] [--epsilon EPS] [--steps N] [--sample K] [--cs kahan|augmented] [--member K] [--dt-file FILE]");
    }
    if (engine_program && !arith_fma) die("--engine program needs --arith fma: a correctly rounded divide is not a sequencer program (cft-fp256 docs/ORBITS.md)");
    if (engine_program && cs_augmented) die("--engine program carries Kahan's add_cs only; --cs augmented is a host operation");
    read_problem(problem);
    make_constants(tol_shift, quiet);
    alloc_state();

    /* dt and epsilon from their decimal text, correctly rounded in the
     * format - at binary64 the same bits a C literal would be - or from
     * a hex float, which is exact (a record's own dt0 fed back in) */
    EPS = from_text(eps_txt, is_hex_text(eps_txt));
    EPS5040 = valloc(1); vmul(EPS5040, EPS, K5040, 1);
    adaptive = s_lt(K0, EPS);
    if (dt_file){
        if (adaptive) die("--dt-file prescribes every step; it needs --epsilon 0");
        read_dt_file(dt_file, steps);
        for (size_t s = 0; s < E; s++) vcopy(E(SDT, s), dtlist[0], 1);
    }else{
        V s0 = from_text(dt_txt, is_hex_text(dt_txt));
        for (size_t s = 0; s < E; s++) vcopy(E(SDT, s), s0, 1);
        free(s0);
    }
    bcast_sys(DTB, SDT);
    if (dt_out_path){ if (E != 1) die("--dt-out records one system's steps; use --member"); dt_out = fopen(dt_out_path, "w"); if (!dt_out){ perror(dt_out_path); exit(2); } }
    use_ens_predict = engine_program && E > 1 && adaptive;
    if (engine_program) load_programs();
    if (do_dump) dump_constants();
    vcopy(x, X0, N3); vcopy(v, V0, N3);
    vzero(csx, N3); vzero(csv, N3);
    vzero(SDTLAST, E); vzero(STPLAIN, E); vzero(STHI, E); vzero(STLO, E);
    vzero(gcs, N3);
    for (int m = 0; m < 7; m++){ vzero(g[m], N3); vzero(b[m], N3); vzero(e[m], N3); vzero(csb[m], N3); vzero(er[m], N3); vzero(br[m], N3); }

    if (E == 1) printf("# cft-rebound ias15 record v1\n");
    else        printf("# cft-rebound ias15 ensemble record v1\n");
    printf("# program=cft impl=libcft-%u.%u format=%s backend=%s problem=%s",
           (unsigned)(cft_abi_version() >> 16), (unsigned)(cft_abi_version() & 0xffff), cft_format_name(F),
           artifact ? artifact : "software", problem_name);
    if (E > 1) printf(" E=%zu", E);
    printf(" N=%zu steps=%ld sample=%ld cs=%s arith=%s engine=%s pc_tol_shift=%d max_iter=%d",
           N, steps, sample, cs_augmented ? "augmented" : "kahan",
           arith_fma ? "fma" : "rebound", engine_program ? "program" : "loop", tol_shift, max_iter);
    if (member >= 0) printf(" member=%d system=%s", member, sys_names[0]);
    if (dt_file) printf(" dt_file=%s", dt_file);
    printf("\n");
    printf("# dt0="); puthex(E(SDT, 0)); printf(" epsilon="); puthex(EPS); printf(" G="); puthex(G); printf(" pc_tol="); puthex(TOL); printf("\n");
    if (E == 1){
        for (size_t i = 0; i < N; i++){ printf("# body %zu %s m=", i, body_names[i]); puthex(E(mass, i)); printf("\n"); }
        printf("# columns: sample step t dt_next dt_last E then per body x y z vx vy vz ; then exact-time pair t_hi t_lo\n");
    }else{
        for (size_t s = 0; s < E; s++){
            printf("# system %zu %s\n", s, sys_names[s]);
            for (size_t i = 0; i < N; i++){ printf("# body %zu %s m=", i, body_names[s * N + i]); puthex(E(mass, s * N + i)); printf("\n"); }
        }
        printf("# columns: system s sample step t dt_next dt_last E then per body x y z vx vy vz ; then exact-time pair t_hi t_lo\n");
    }

    V Eo = valloc(E);
    clock_t c0 = clock();
    long k = 0, blocks_done = 0;
    while (1){
        energy_all(Eo);
        for (size_t s = 0; s < E; s++){
            if (E > 1) printf("system %zu ", s);
            printf("sample %ld %ld", k, done[s]);
            puthex(E(STPLAIN, s)); puthex(E(SDT, s)); puthex(E(SDTLAST, s)); puthex(E(Eo, s));
            for (size_t i = L * s; i < L * (s + 1); i += 3){ puthex(E(x, i)); puthex(E(x, i + 1)); puthex(E(x, i + 2)); puthex(E(v, i)); puthex(E(v, i + 1)); puthex(E(v, i + 2)); }
            printf(" |"); puthex(E(STHI, s)); puthex(E(STLO, s));
            printf("\n");
        }
        fflush(stdout);
        if (blocks_done >= steps) break;
        long todo = sample; if (blocks_done + todo > steps) todo = steps - blocks_done;
        long target = blocks_done + todo;
        while (1){
            int need = 0;
            for (size_t s = 0; s < E; s++){ active[s] = done[s] < target; need |= active[s]; }
            if (!need) break;
            step_attempt();
        }
        blocks_done = target;
        k++;
    }
    double secs = (double)(clock() - c0) / CLOCKS_PER_SEC;
    long rej_all = 0, mx_all = 0, att_all = 0; unsigned long long pc_all = 0; int pcmax_all = 0;
    for (size_t s = 0; s < E; s++){
        rej_all += rejected[s]; mx_all += max_exceeded[s]; att_all += done[s] + rejected[s]; pc_all += pc_total[s];
        if (pc_max[s] > pcmax_all) pcmax_all = pc_max[s];
        if (E > 1)
            printf("# system %zu steps_done=%ld iterations_max_exceeded=%ld steps_rejected=%ld mean_pc_iterations=%.3f max_pc_iterations=%d\n",
                   s, done[s], max_exceeded[s], rejected[s],
                   done[s] + rejected[s] ? (double)pc_total[s] / (double)(done[s] + rejected[s]) : 0.0, pc_max[s]);
    }
    printf("# steps_done=%ld iterations_max_exceeded=%ld steps_rejected=%ld mean_pc_iterations=%.3f max_pc_iterations=%d flags_seen=0x%02x calls=%llu divsqrt_calls=%llu seconds=%.3f steps_per_s=%.2f",
           done[0], mx_all, rej_all,
           att_all ? (double)pc_all / (double)att_all : 0.0, pcmax_all,
           flags_union, ncalls, ncalls_divsqrt, secs, secs > 0 ? done[0] / secs : 0.0);
    if (E > 1) printf(" E=%zu system_steps_per_s=%.2f pc_lane_efficiency=%.3f", E, secs > 0 ? (double)done[0] * (double)E / secs : 0.0,
                      pc_lane_max ? (double)pc_lane_sum / (double)pc_lane_max : 0.0);
    printf("\n");
    if (!quiet) fprintf(stderr, "done: %ld steps%s, %llu library calls, %.1f s, flags 0x%02x\n", done[0], E > 1 ? " per system" : "", ncalls, secs, flags_union);
    if (dt_out){
        /* the step the run would take next, so that a replay of this
         * file forms the same final prediction ratio and reports the
         * same dt_next: steps + 1 lines in all */
        char buf[160]; size_t len = 0; cft_to_hex_char(dev, F, E(SDT, 0), buf, sizeof buf, &len); fprintf(dt_out, "%s\n", buf);
        fclose(dt_out);
    }
    cft_close(dev);
    return 0;
}
