/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The wave program (ROUND5.md, P2): src/wave/'s step, one system at a
 * time on the host (the twin, `make wave-host`) or one system to a
 * workgroup on a GPU (`make HIP=1 wave`, the same file built by hipcc,
 * with --device gpu), writing this repository's record format so that
 * tools/compare_run.py holds it to `ias15_cft --arith fma --format fp64`
 * bit for bit, every column, member by member.
 *
 *   wave --problem FILE [--dt DT] [--epsilon EPS] [--steps N] [--sample K]
 *        [--max-iter M] [--pc-tol-shift S] [--member K] [--no-flag-abort]
 *        [--device host|gpu] [--quiet] [--attempts-dir DIR]
 *        [--format fp64] [--arith fma] [--cs kahan]
 *   wave --dump-constants
 *   wave --zonal-selftest-from FILE      (ias15_cft --zonal-selftest's output)
 *
 * What it takes is the port's: the problem file's lines (name, G, N, E,
 * system, body, active, zonal), --dt and --epsilon as a decimal (strtod,
 * correctly rounded) or an exact hex float, and the run's shape. What it
 * refuses by name is everything the port has and this does not: any
 * format but fp64, REBOUND's form, --cs augmented, adaptive modes other
 * than PRS23, softening, a step floor, dt files, ensembles of more than
 * WAVE_MAXN bodies a system, and a non-finite input value (the port would
 * take one and abort on its first flag).
 *
 * The record names its producer: program=wave, engine=wave, backend=host
 * or the GPU's, and a `# device` line with the GPU's own description.
 * tools/compare_run.py reruns it with --rerun-engine loop.
 */
#include "wave_step.h"
#include "wave_io.h"
#include "ias15_constants.h"
#include "ias15_limits.h"
#include <stdarg.h>
#include <time.h>
#include <string>
#include <vector>

using namespace wave;

static_assert(IAS15_ZONAL_J2_NUM == ZJ2_NUM && IAS15_ZONAL_J2_DEN == ZJ2_DEN && IAS15_ZONAL_J2_S1 == ZJ2_S1 &&
              IAS15_ZONAL_J2_S0 == ZJ2_S0 && IAS15_ZONAL_J2_Z0 == ZJ2_Z0 && IAS15_ZONAL_J4_NUM == ZJ4_NUM &&
              IAS15_ZONAL_J4_DEN == ZJ4_DEN && IAS15_ZONAL_J4_S2 == ZJ4_S2 && IAS15_ZONAL_J4_S1 == ZJ4_S1 &&
              IAS15_ZONAL_J4_S0 == ZJ4_S0 && IAS15_ZONAL_J4_Z0 == ZJ4_Z0 && IAS15_ZONAL_J4_Z1 == ZJ4_Z1 &&
              IAS15_ZONAL_P2_S1 == ZP2_S1 && IAS15_ZONAL_P2_S0 == ZP2_S0 && IAS15_ZONAL_P2_DEN == ZP2_DEN &&
              IAS15_ZONAL_P4_S2 == ZP4_S2 && IAS15_ZONAL_P4_S1 == ZP4_S1 && IAS15_ZONAL_P4_S0 == ZP4_S0 &&
              IAS15_ZONAL_P4_DEN == ZP4_DEN && IAS15_N_ZONAL == ZN,
              "src/wave/wave_step.h's zonal indices are not src/ias15_constants.h's");
static_assert(IAS15_N_H == 8 && IAS15_N_RR == 28 && IAS15_N_C == 21 && IAS15_N_D == 21 && IAS15_N_RINV == 28 &&
              IAS15_N_HF == 49, "src/ias15_constants.h's array sizes are not the step's");

static int quiet;

static void die(int rc, const char *fmt, ...){
    va_list ap; va_start(ap, fmt);
    fprintf(stderr, "wave: ");
    vfprintf(stderr, fmt, ap);
    fprintf(stderr, "\n");
    va_end(ap);
    fflush(stdout);
    exit(rc);
}

static double now_s(void){
    struct timespec ts; clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

static const char *const phase_name[PH_N] = {
    "nothing", "gravity", "J2 and J4", "the step's start (g from b)", "the predictor", "the corrector",
    "the corrector's error", "the step control (PRS23)", "the end-of-step update", "the exact-time clock",
    "the next step's prediction", "the energy", "sqrt7", "the constants" };
static const char *const kind_name[8] = { "?", "add", "sub", "mul", "div", "sqrt", "fma", "augmented_add" };
static std::string flag_names(unsigned f){
    std::string s;
    const char *n[4] = { "INVALID", "DIVBYZERO", "OVERFLOW", "UNDERFLOW" };
    for (int i = 0; i < 4; i++) if (f >> i & 1u){ if (!s.empty()) s += "|"; s += n[i]; }
    return s.empty() ? "none" : s;
}

/* ---- the problem ------------------------------------------------------- */
struct Body { std::string name; double m, x[3], v[3]; };
struct Problem {
    std::string name;
    double G = 0.0;
    int N = 0, E = 1, Efile = 1;
    std::vector<std::string> sysnames;
    std::vector<Body> bodies;               /* E * N, the systems' bodies in order */
    long active = -1;
    int zonal_on = 0, zonal_k = 0;
    double zr = 0, zj2 = 0, zj4 = 0;
};

static double hexval_or_die(const char *s, const char *what, const char *path){
    double v;
    if (!s) die(2, "a short %s line in %s", what, path);
    if (!wave_from_hex(s, &v)) die(3, "cannot parse '%s' (%s in %s)", s, what, path);
    if (!is_finite(v)) die(3, "'%s' (%s in %s) is not finite: the port would take it and abort on its first flag; "
                           "this program refuses it before anything runs", s, what, path);
    return v;
}

/* read_problem() (src/ias15_cft.c:1027), line for line in what it takes and refuses */
static Problem read_problem(const char *path, int member){
    FILE *f = fopen(path, "r");
    if (!f){ perror(path); exit(2); }
    char line[2048];
    long nb = 0, ns = 0, nline = -1;
    while (fgets(line, sizeof line, f)){
        if (!strncmp(line, "body ", 5)) nb++;
        else if (!strncmp(line, "system ", 7) || !strcmp(line, "system\n")) ns++;
        else if (!strncmp(line, "N ", 2)) nline = atol(line + 2);
    }
    rewind(f);
    Problem P;
    P.Efile = ns ? (int)ns : 1;
    if (nline < 0){ if (P.Efile > 1) die(3, "an ensemble file needs an N line (bodies per system)"); nline = nb; }
    P.N = (int)nline;
    if (P.N < 1 || P.N > CFT_MAX_BODIES) die(3, "N = %d bodies per system is outside 1..%d", P.N, CFT_MAX_BODIES);
    if (P.N > WAVE_MAXN) die(3, "N = %d bodies per system: this program's systems hold at most %d (WAVE_MAXN, one workgroup)",
                             P.N, WAVE_MAXN);
    if (nb != (long)P.Efile * P.N) die(3, "%ld body lines for %d system(s) of %d bodies in %s", nb, P.Efile, P.N, path);
    if (member >= 0){
        if (member >= P.Efile) die(3, "--member %d, but %s has %d system(s)", member, path, P.Efile);
        P.E = 1;
    }else P.E = P.Efile;
    P.sysnames.assign(P.E, "");
    long cur = -1, in_cur = 0;
    int have_G = 0;
    while (fgets(line, sizeof line, f)){
        if (line[0] == '#' || line[0] == '\n') continue;
        char *tok = strtok(line, " \t\r\n");
        if (!tok) continue;
        if (!strcmp(tok, "name")){ const char *t = strtok(NULL, " \t\r\n"); P.name = std::string(t ? t : "").substr(0, 63); }
        else if (!strcmp(tok, "G")){ P.G = hexval_or_die(strtok(NULL, " \t\r\n"), "G", path); have_G = 1; }
        else if (!strcmp(tok, "N")){ }
        else if (!strcmp(tok, "E")){
            const char *t = strtok(NULL, " \t\r\n");
            if (!t || atol(t) != P.Efile) die(3, "E line disagrees with the system count");
        }
        else if (!strcmp(tok, "system")){
            if (cur >= 0 && in_cur != P.N) die(3, "system %ld has %ld bodies, not %d", cur, in_cur, P.N);
            cur++; in_cur = 0;
            const char *nm = strtok(NULL, " \t\r\n");
            const std::string s = std::string(nm ? nm : "").substr(0, 31);
            if (member < 0) P.sysnames[cur] = s;
            else if (cur == member) P.sysnames[0] = s;
        }
        else if (!strcmp(tok, "body")){
            if (ns && cur < 0) die(3, "body line before the first system line in %s", path);
            in_cur++;
            if (member >= 0 && cur != member) continue;
            if ((long)P.bodies.size() >= (long)P.E * P.N) die(3, "too many body lines");
            Body b;
            const char *nm = strtok(NULL, " \t\r\n");
            b.name = std::string(nm ? nm : "").substr(0, 31);
            const char *fld[7];
            for (int k = 0; k < 7; k++){ fld[k] = strtok(NULL, " \t\r\n"); if (!fld[k]) die(3, "short body line"); }
            b.m = hexval_or_die(fld[0], "a body's mass", path);
            for (int k = 0; k < 3; k++) b.x[k] = hexval_or_die(fld[1 + k], "a body's position", path);
            for (int k = 0; k < 3; k++) b.v[k] = hexval_or_die(fld[4 + k], "a body's velocity", path);
            P.bodies.push_back(b);
        }
        else if (!strcmp(tok, "active")){
            const char *t = strtok(NULL, " \t\r\n");
            const char *more = t ? strtok(NULL, " \t\r\n") : NULL;
            char *end = NULL;
            const long n = t ? strtol(t, &end, 10) : -1;
            if (P.active >= 0) die(3, "a second active line in %s: a problem has one N_active", path);
            if (more) die(3, "'active %s %s' in %s: the line is 'active n', one field", t, more, path);
            if (!t || *end || n < 1 || n > P.N)
                die(3, "'active %s' in %s: it takes the count of massive bodies a system, 1 to N = %d", t ? t : "", path, P.N);
            P.active = n;
        }
        else if (!strcmp(tok, "zonal")){
            if (P.zonal_on) die(3, "a second zonal line in %s: a problem has one oblate body, and the first zonal line named body %d",
                                path, P.zonal_k);
            const char *fld[4];
            for (int k = 0; k < 4; k++){ fld[k] = strtok(NULL, " \t\r\n"); if (!fld[k]) die(3, "a short zonal line in %s: it is 'zonal k R J2 J4'", path); }
            const char *more = strtok(NULL, " \t\r\n");
            if (more) die(3, "a zonal line in %s with a fifth field '%s': it is 'zonal k R J2 J4'", path, more);
            char *end = NULL;
            const long k = strtol(fld[0], &end, 10);
            if (*end || k < 0 || k >= P.N) die(3, "'zonal %s' in %s: k is a body of each system, 0 to %d", fld[0], path, P.N - 1);
            P.zonal_k = (int)k;
            P.zr = hexval_or_die(fld[1], "the zonal R", path);
            P.zj2 = hexval_or_die(fld[2], "the zonal J2", path);
            P.zj4 = hexval_or_die(fld[3], "the zonal J4", path);
            P.zonal_on = 1;
        }
        else die(3, "unknown line '%s' in %s", tok, path);
    }
    fclose(f);
    if (ns && in_cur != P.N) die(3, "the last system has %ld bodies, not %d", in_cur, P.N);
    if ((long)P.bodies.size() != (long)P.E * P.N) die(3, "read %zu bodies, expected %d", P.bodies.size(), P.E * P.N);
    if (!have_G) die(3, "no G line in %s", path);
    if (!ns) P.sysnames[0] = P.name.substr(0, 31);
    return P;
}

/* ---- the constants ----------------------------------------------------- */
static double kconst(const char *hex, const char *dec, const char *name, int i){
    double h, d;
    if (!wave_from_hex(hex, &h)) die(3, "constant %s[%d]: cannot parse %s", name, i, hex);
    d = strtod(dec, NULL);
    if (d2u(h) != d2u(d))
        die(3, "constant %s[%d]: hex form %s and decimal form %s round to different binary64 bits", name, i, hex, dec);
    return h;
}

static void make_constants(Consts &K, int tol_shift){
    memset(&K, 0, sizeof K);
    for (int i = 0; i < 8; i++) K.h[i] = kconst(ias15_h_hex[0][i], ias15_h_dec[i], "h", i);
    for (int i = 0; i < 28; i++) K.rr[i] = kconst(ias15_rr_hex[0][i], ias15_rr_dec[i], "rr", i);
    for (int i = 0; i < 21; i++) K.c[i] = kconst(ias15_c_hex[0][i], ias15_c_dec[i], "c", i);
    for (int i = 0; i < 21; i++) K.d[i] = kconst(ias15_d_hex[0][i], ias15_d_dec[i], "d", i);
    for (int i = 0; i < 28; i++) K.rinv[i] = kconst(ias15_rinv_hex[0][i], ias15_rinv_dec[i], "rinv", i);
    for (int i = 0; i < 49; i++) K.hf[i] = kconst(ias15_hf_hex[0][i], ias15_hf_dec[i], "hf", i);
    for (int i = 0; i < ZN; i++) K.zk[i] = kconst(ias15_zonal_hex[0][i], ias15_zonal_dec[i], "zonal", i);
    /* the reciprocals by the library's divide (make_constants(), :917): here the twin's own */
    Fl fl; fl.clear(); fl.phase = PH_SETUP;
    for (int m = 0; m < 7; m++){
        K.kposr[m] = w_div(1.0, (double)((m + 2) * (m + 3)), fl);
        K.kvelr[m] = w_div(1.0, (double)(m + 2), fl);
    }
    K.kposr[7] = w_div(1.0, 2.0, fl);
    K.kvelr[7] = w_div(1.0, 1.0, fl);
    K.khalf = w_div(1.0, 2.0, fl);
    K.invb = w_div(1.0, 0.25, fl);                               /* 1/safety_factor (choose_timestep()) */
    if (fl.bits) die(3, "the constants raised %s", flag_names(fl.bits).c_str());
    /* the tolerance: 1e-16 scaled by an exact power of two (cft_scaleb) */
    K.tol = ldexp(1e-16, -tol_shift);
    if (ldexp(K.tol, tol_shift) != 1e-16) die(3, "the tolerance scaling was inexact - that must never happen");
}

static void dump_constants(const Consts &K){
    char b[40];
    printf("# constants at fp64 (hex = exact bits)\n");
    for (int i = 0; i < 8; i++){ wave_to_hex(K.h[i], b); printf("h[%d] %s\n", i, b); }
    for (int i = 0; i < 28; i++){ wave_to_hex(K.rr[i], b); printf("rr[%d] %s\n", i, b); }
    for (int i = 0; i < 21; i++){ wave_to_hex(K.c[i], b); printf("c[%d] %s\n", i, b); }
    for (int i = 0; i < 21; i++){ wave_to_hex(K.d[i], b); printf("d[%d] %s\n", i, b); }
    wave_to_hex(K.tol, b); printf("pc_tol %s\n", b);
}

static double text_value(const char *s, const char *what){
    double v;
    const char *t = s + (*s == '-' || *s == '+');
    if (t[0] == '0' && (t[1] == 'x' || t[1] == 'X')){
        if (!wave_from_hex(s, &v)) die(3, "cannot parse %s '%s'", what, s);
    }else{
        char *end = NULL;
        v = strtod(s, &end);
        if (!end || *end || end == s) die(3, "cannot parse %s '%s'", what, s);
    }
    return v;
}

/* ---- the run key (src/ias15_cft.c:4420, sched_key()) --------------------
 * What a member's attempts record and schedule are filed under, for the
 * driver's particle replay. A copy, held by tools/check_wave.py to
 * `ias15_cft --run-key --member k` for every member it writes. */
static std::string run_key(const Problem &P, const Consts &K, int member, long steps, long sample, int tol_shift,
                           const char *engine, double dt0, const Body *b){
    char h[5][40], z[3][40];
    const double val[5] = { dt0, K.eps, K.soft2, K.mindt, K.G };
    for (int i = 0; i < 5; i++) wave_to_hex(val[i], h[i]);
    uint64_t hs = CFT_DIGEST_SEED;
    for (int i = 0; i < P.N; i++){ uint64_t u = d2u(b[i].m); for (int k = 0; k < 8; k++){ hs ^= (u >> (8 * k)) & 0xff; hs *= 1099511628211ULL; } }
    for (int i = 0; i < P.N; i++)
        for (int c = 0; c < 3; c++){
            const uint64_t ux = d2u(b[i].x[c]), uv = d2u(b[i].v[c]);
            for (int k = 0; k < 8; k++){
                hs ^= (ux >> (8 * k)) & 0xff; hs *= 1099511628211ULL;
                hs ^= (uv >> (8 * k)) & 0xff; hs *= 1099511628211ULL;
            }
        }
    char out[1536];
    int w = snprintf(out, sizeof out, "format=fp64 problem=%s E=1 N=%d member=%d steps=%ld sample=%ld adaptive_mode=2 "
                     "max_iter=%d pc_tol_shift=%d cs=kahan arith=fma engine=%s dt0=%s epsilon=%s softening2=%s "
                     "min_dt=%s G=%s dt_file=no state=%016llx", P.name.c_str(), P.N, member, steps, sample, K.max_iter,
                     tol_shift, engine, h[0], h[1], h[2], h[3], h[4], (unsigned long long)hs);
    if (P.active >= 0) w += snprintf(out + w, sizeof out - w, " active=%ld", P.active);
    if (P.zonal_on){
        wave_to_hex(P.zr, z[0]); wave_to_hex(P.zj2, z[1]); wave_to_hex(P.zj4, z[2]);
        w += snprintf(out + w, sizeof out - w, " zonal=%d,%s,%s,%s", P.zonal_k, z[0], z[1], z[2]);
    }
    return std::string(out);
}

/* ---- what a block leaves of each system --------------------------------- */
struct SampleOut {
    double t, dt, dtlast, energy, thi, tlo;
    double x[3][WAVE_MAXN], v[3][WAVE_MAXN];
    long long done, rejected, max_exceeded, pc_total, attempts;
    int pc_max, stopped;
    unsigned flag_bits, flag_first;
};

__host__ __device__ inline void fill_sample(const Sys &S, const Shared &sh, int N, SampleOut &o){
    o.t = S.tplain; o.dt = S.dt; o.dtlast = S.dtlast; o.energy = sh.energy; o.thi = S.thi; o.tlo = S.tlo;
    for (int c = 0; c < 3; c++) for (int p = 0; p < N; p++){ o.x[c][p] = S.x[c][p]; o.v[c][p] = S.v[c][p]; }
    o.done = S.done; o.rejected = S.rejected; o.max_exceeded = S.max_exceeded; o.pc_total = S.pc_total;
    o.attempts = S.attempts; o.pc_max = S.pc_max; o.stopped = S.stopped; o.flag_bits = S.flag_bits; o.flag_first = S.flag_first;
}

/* One attempt's line for a member's attempts record and schedule */
struct Attempt { double dt; int accepted, iters, capped; };

#ifdef WAVE_HIP
#include "wave_gpu.inc"
#endif

/* ---- the host twin -------------------------------------------------------- */
struct HostRunner {
    const Consts &K;
    std::vector<Sys> &S;
    std::vector<std::vector<Attempt>> *log;
    Shared sh;
    HostRunner(const Consts &K_, std::vector<Sys> &S_, std::vector<std::vector<Attempt>> *log_) : K(K_), S(S_), log(log_){}
    void block(long long target, std::vector<SampleOut> &out){
        for (size_t s = 0; s < S.size(); s++){
            HostExec ex{K.N};
            Fl fl; fl.clear();
            Sys &y = S[s];
            while (y.done < target && !y.stopped){
                step_attempt(ex, K, y, sh, fl);
                gather_flags(ex, K, y, sh, fl);
                if (log) (*log)[s].push_back(Attempt{ y.last_dt, y.last_accepted, y.last_iters, y.last_capped });
            }
            if (!y.stopped){ energy(ex, K, y, sh, fl); gather_flags(ex, K, y, sh, fl); }
            fill_sample(y, sh, K.N, out[s]);
        }
    }
};

/* ---- --zonal-selftest-from: the port's own points, this program's term ----
 * On the host, or with --device gpu in a kernel, a workgroup of eight a
 * system. */
static int zonal_selftest_from(const char *path, Consts &K, int gpu){
    FILE *f = fopen(path, "r");
    if (!f){ perror(path); return 2; }
    char line[4096];
    int E = 0, N = 0, k = -1, act = -1;
    double G = 0, R = 0, J2 = 0, J4 = 0;
    static double m[2][8], x[2][3][8];
    int got = 0;
    auto bits = [](const char *s){ return u2d(strtoull(s, NULL, 16)); };
    while (fgets(line, sizeof line, f)){
        int s, j;
        char a[40], b[40], c[40], d[40];
        if (sscanf(line, "zonal selftest format=fp64 E=%d N=%d k=%d active=%d", &E, &N, &k, &act) == 4) continue;
        if (sscanf(line, "zst G %39s", a) == 1){ G = bits(a); continue; }
        if (sscanf(line, "zst R %39s", a) == 1){ R = bits(a); continue; }
        if (sscanf(line, "zst J2 %39s", a) == 1){ J2 = bits(a); continue; }
        if (sscanf(line, "zst J4 %39s", a) == 1){ J4 = bits(a); continue; }
        if (sscanf(line, "zst body %d %d %39s %39s %39s %39s", &s, &j, a, b, c, d) == 6){
            if (s < 0 || s > 1 || j < 0 || j > 7) die(3, "--zonal-selftest-from: body %d %d is outside 2 systems of 8", s, j);
            m[s][j] = bits(a); x[s][0][j] = bits(b); x[s][1][j] = bits(c); x[s][2][j] = bits(d); got++;
        }
    }
    fclose(f);
    if (E != 2 || N != 8 || k < 0 || act < 1 || got != 16)
        die(3, "--zonal-selftest-from %s: not ias15_cft --zonal-selftest's fp64 output (E=%d N=%d k=%d active=%d, %d bodies)",
            path, E, N, k, act, got);
    K.N = N; K.na = act; K.zonal_on = 1; K.zonal_k = k; K.G = G; K.zr = R; K.zj2 = J2; K.zj4 = J4;
    std::vector<Sys> S(E);
    for (int s = 0; s < E; s++){
        memset(&S[s], 0, sizeof S[s]);
        for (int j = 0; j < N; j++){ S[s].m[j] = m[s][j]; for (int c = 0; c < 3; c++) S[s].x[c][j] = x[s][c][j]; }
    }
    double pot[2] = { 0.0, 0.0 };
    unsigned flags = 0;
    std::string dline;
#ifdef WAVE_HIP
    if (gpu){ std::string b; dline = gpu_open(b); flags = gpu_zonal_selftest(K, S, pot); }
    else
#endif
    {
        (void)gpu;
        static Shared sh;
        for (int s = 0; s < E; s++){
            Fl fl; fl.clear();
            HostExec ex{N};
            zonal_selftest_one(ex, K, S[s], sh, fl, &pot[s]);
            flags |= fl.bits;
        }
    }
    printf("wave zonal selftest format=fp64 E=%d N=%d k=%d active=%d on the %s\n", E, N, k, act, gpu ? "GPU" : "host");
    if (!dline.empty()) printf("# device %s\n", dline.c_str());
    for (int s = 0; s < E; s++)
        for (int j = 0; j < N; j++){
            printf("zst acc %d %d", s, j);
            for (int c = 0; c < 3; c++) printf(" %016llx", (unsigned long long)d2u(S[s].a[c][j]));
            printf("\n");
        }
    for (int s = 0; s < E; s++) printf("zst pot %d %016llx\n", s, (unsigned long long)d2u(pot[s]));
    printf("wave zonal selftest: flags raised 0x%x\n", flags);
    return 0;
}

/* ---- main -------------------------------------------------------------- */
int main(int argc, char **argv){
    const char *problem = NULL, *dt_txt = "0.01", *eps_txt = "1e-9", *device = "host", *att_dir = NULL;
    const char *zst = NULL, *key_engine = "loop";
    long steps = 1000, sample = 100;
    int max_iter = 12, tol_shift = 0, member = -1, flag_abort = 1, do_dump = 0;
    for (int i = 1; i < argc; i++){
        const char *a = argv[i];
        const int more = i + 1 < argc;
        if (!strcmp(a, "--problem") && more) problem = argv[++i];
        else if (!strcmp(a, "--dt") && more) dt_txt = argv[++i];
        else if (!strcmp(a, "--epsilon") && more) eps_txt = argv[++i];
        else if (!strcmp(a, "--steps") && more) steps = atol(argv[++i]);
        else if (!strcmp(a, "--sample") && more) sample = atol(argv[++i]);
        else if (!strcmp(a, "--max-iter") && more) max_iter = atoi(argv[++i]);
        else if (!strcmp(a, "--pc-tol-shift") && more) tol_shift = atoi(argv[++i]);
        else if (!strcmp(a, "--member") && more) member = atoi(argv[++i]);
        else if (!strcmp(a, "--device") && more) device = argv[++i];
        else if (!strcmp(a, "--attempts-dir") && more) att_dir = argv[++i];
        else if (!strcmp(a, "--key-engine") && more) key_engine = argv[++i];
        else if (!strcmp(a, "--zonal-selftest-from") && more) zst = argv[++i];
        else if (!strcmp(a, "--no-flag-abort")) flag_abort = 0;
        else if (!strcmp(a, "--dump-constants")) do_dump = 1;
        else if (!strcmp(a, "--quiet")) quiet = 1;
        else if (!strcmp(a, "--format") && more){ if (strcmp(argv[++i], "fp64")) die(3, "--format %s: the wave kernel is the hardware's binary64 only (fp64); binary128 and binary256 wait for a device divide (ROUND5.md)", argv[i]); }
        else if (!strcmp(a, "--arith") && more){ if (strcmp(argv[++i], "fma")) die(3, "--arith %s: the wave kernel is the port's FMA form only (--arith fma)", argv[i]); }
        else if (!strcmp(a, "--cs") && more){ if (strcmp(argv[++i], "kahan")) die(3, "--cs %s: the wave kernel carries Kahan's add_cs only (--cs kahan)", argv[i]); }
        else if (!strcmp(a, "--adaptive-mode") && more){ if (strcmp(argv[++i], "2")) die(3, "--adaptive-mode %s: the wave kernel is PRS23 only (2)", argv[i]); }
        else if (!strcmp(a, "--softening") || !strcmp(a, "--min-dt") || !strcmp(a, "--dt-file") || !strcmp(a, "--engine"))
            die(3, "%s: not built in the wave kernel", a);
        else die(3, "unknown argument %s", a);
    }
    Consts K;
    if (tol_shift < 0) die(3, "--pc-tol-shift %d: negative", tol_shift);
    make_constants(K, tol_shift);
    if (do_dump && !problem){ dump_constants(K); return 0; }
    const int gpu = !strcmp(device, "gpu");
    if (!gpu && strcmp(device, "host")) die(3, "--device %s: host or gpu", device);
#ifndef WAVE_HIP
    if (gpu) die(3, "--device gpu: this build has no GPU (make HIP=1 wave)");
#endif
    if (zst) return zonal_selftest_from(zst, K, gpu);
    if (!problem) die(2, "usage: wave --problem FILE [--dt DT] [--epsilon EPS] [--steps N] [--sample K] [--member K] [--device host|gpu]");
    if (steps < 0 || sample < 1) die(3, "--steps %ld --sample %ld", steps, sample);
    if (max_iter < 1) die(3, "--max-iter %d", max_iter);

    Problem P = read_problem(problem, member);
    K.N = P.N;
    K.na = P.active >= 0 ? (int)P.active : P.N;
    K.G = P.G; K.soft2 = 0.0; K.mindt = 0.0;
    K.zonal_on = P.zonal_on; K.zonal_k = P.zonal_k; K.zr = P.zr; K.zj2 = P.zj2; K.zj4 = P.zj4;
    if (K.zonal_on && K.zonal_k >= K.na)
        die(3, "the zonal line's body %d is a test particle (active %d): J2 and J4 belong to one of the first "
            "N_active bodies, which act on the others", K.zonal_k, K.na);
    K.max_iter = max_iter;
    K.flag_abort = flag_abort;
    K.abort_mask = F_INVALID | F_DIVBYZERO | F_OVERFLOW | F_UNDERFLOW;
    K.eps = text_value(eps_txt, "--epsilon");
    {
        Fl fl; fl.clear(); fl.phase = PH_SETUP;
        K.eps5040 = w_mul(K.eps, 5040.0, fl);
        if (fl.bits && flag_abort) die(3, "epsilon*5040 raised %s", flag_names(fl.bits).c_str());
    }
    K.adaptive = 0.0 < K.eps;
    const double dt0 = text_value(dt_txt, "--dt");
    if (do_dump) dump_constants(K);

    const int E = P.E, N = P.N;
    std::vector<Sys> S(E);
    for (int s = 0; s < E; s++){
        Sys &y = S[s];
        memset(&y, 0, sizeof y);
        for (int p = 0; p < N; p++){
            const Body &b = P.bodies[(size_t)s * N + p];
            y.m[p] = b.m;
            for (int c = 0; c < 3; c++){ y.x[c][p] = b.x[c]; y.v[c][p] = b.v[c]; }
        }
        y.dt = dt0;
        y.sampled_at = -1;
    }

    /* ---- the header: main()'s, line for line (src/ias15_cft.c:4696) ---- */
    char hb[40];
    auto puthex = [&](double v){ wave_to_hex(v, hb); printf(" %s", hb); };
    std::string backend = gpu ? "gpu" : "host";
    std::string devline;
#ifdef WAVE_HIP
    if (gpu) devline = gpu_open(backend);
#endif
    printf(E == 1 ? "# cft-rebound ias15 record v1\n" : "# cft-rebound ias15 ensemble record v1\n");
    printf("# program=wave impl=cft-rebound-wave-1 format=fp64 backend=%s problem=%s", backend.c_str(), P.name.c_str());
    if (E > 1) printf(" E=%d", E);
    printf(" N=%d steps=%ld sample=%ld cs=kahan arith=fma engine=wave pc_tol_shift=%d max_iter=%d", N, steps, sample,
           tol_shift, max_iter);
    if (member >= 0) printf(" member=%d system=%s", member, P.sysnames[0].c_str());
    printf("\n");
    {
        std::string cv = __VERSION__;
        for (auto &ch : cv) if (ch == ' ') ch = ',';
        printf("# identity wave_source=src/wave compiler=%s contract=off fma=hardware flags=checked\n", cv.c_str());
    }
    if (!devline.empty()) printf("# device %s\n", devline.c_str());
    printf("# dt0="); puthex(dt0); printf(" epsilon="); puthex(K.eps); printf(" G="); puthex(K.G); printf(" pc_tol="); puthex(K.tol);
    printf("\n");
    if (P.active >= 0){
        /* build_scatter()'s count for system 0: active-active, then test-active */
        const long na = P.active;
        long np = 0;
        for (long i = 1; i < na; i++) np += i;
        for (long i = na > 1 ? na : 1; i < N; i++) np += na;
        printf("# active %ld of %d bodies a system: %ld pairs a force evaluation, a system\n", P.active, N, np);
    }
    if (P.zonal_on){
        int zreact = 0;
        for (int j = 0; j < N; j++) zreact += zonal_gives(K, j);
        printf("# zonal body %d of every system, about z: R=", P.zonal_k); puthex(P.zr);
        printf(" J2="); puthex(P.zj2); printf(" J4="); puthex(P.zj4);
        printf("; it acts on %d bodies a system, and body %d takes the reaction of %d\n", N - 1, P.zonal_k, zreact);
    }
    if (E == 1){
        for (int i = 0; i < N; i++){ printf("# body %d %s m=", i, P.bodies[i].name.c_str()); puthex(P.bodies[i].m); printf("\n"); }
        printf("# columns: sample step t dt_next dt_last E then per body x y z vx vy vz ; then exact-time pair t_hi t_lo\n");
    }else{
        for (int s = 0; s < E; s++){
            printf("# system %d %s\n", s, P.sysnames[s].c_str());
            for (int i = 0; i < N; i++){
                printf("# body %d %s m=", i, P.bodies[(size_t)s * N + i].name.c_str()); puthex(P.bodies[(size_t)s * N + i].m);
                printf("\n");
            }
        }
        printf("# columns: system s sample step t dt_next dt_last E then per body x y z vx vy vz ; then exact-time pair t_hi t_lo\n");
    }
    fflush(stdout);

    std::vector<std::vector<Attempt>> log(E);
    std::vector<SampleOut> out(E);
    HostRunner host(K, S, att_dir ? &log : NULL);
#ifdef WAVE_HIP
    GpuRunner *dev = gpu ? new GpuRunner(K, S, att_dir ? &log : NULL) : NULL;
#endif
    const double w0 = now_s();
    double t_kernel = 0.0, t_first = 0.0;   /* every launch's time, and the sample-0 launch's (the device's first) */
    long long launches = 0;
    long k = 0, blocks_done = 0;
    for (;;){
        const long long target = blocks_done;
#ifdef WAVE_HIP
        if (dev){
            const double t0 = now_s();
            dev->block(target, out);
            const double d = now_s() - t0;
            t_kernel += d;
            if (target == 0) t_first += d;
        }
        else
#endif
            host.block(target, out);
        for (int s = 0; s < E; s++){
            const SampleOut &o = out[s];
            if (o.stopped){
                const unsigned ff = o.flag_first;
                die(3, "system %d (%s) stopped: %s raised %s (0x%02x) in %s; every flag this system raised: %s. This workload "
                    "can only raise inexact; invalid, divide-by-zero, overflow or underflow means the integration or the "
                    "port is wrong (pass --no-flag-abort to continue anyway)", s, P.sysnames[s].c_str(),
                    kind_name[(ff >> 4) & 0xff & 7], flag_names(ff & 0xf).c_str(), ff & 0xf,
                    phase_name[(ff >> 12) < PH_N ? (ff >> 12) : 0], flag_names(o.flag_bits).c_str());
            }
        }
        for (int s = 0; s < E; s++){
            const SampleOut &o = out[s];
            if (E > 1) printf("system %d ", s);
            printf("sample %ld %lld", k, o.done);
            puthex(o.t); puthex(o.dt); puthex(o.dtlast); puthex(o.energy);
            for (int p = 0; p < N; p++){ puthex(o.x[0][p]); puthex(o.x[1][p]); puthex(o.x[2][p]); puthex(o.v[0][p]); puthex(o.v[1][p]); puthex(o.v[2][p]); }
            printf(" |"); puthex(o.thi); puthex(o.tlo);
            printf("\n");
        }
        fflush(stdout);
        if (blocks_done >= steps) break;
        long todo = sample; if (blocks_done + todo > steps) todo = steps - blocks_done;
        blocks_done += todo;
        k++;
    }
    const double wall = now_s() - w0;
#ifdef WAVE_HIP
    if (dev) launches = dev->launches;
#endif
    long long rej_all = 0, mx_all = 0, att_all = 0, pc_all = 0; int pcmax_all = 0; unsigned fl_all = 0;
    for (int s = 0; s < E; s++){
        const SampleOut &o = out[s];
        rej_all += o.rejected; mx_all += o.max_exceeded; att_all += o.done + o.rejected; pc_all += o.pc_total;
        if (o.pc_max > pcmax_all) pcmax_all = o.pc_max;
        fl_all |= o.flag_bits;
        if (E > 1)
            printf("# system %d steps_done=%lld iterations_max_exceeded=%lld steps_rejected=%lld mean_pc_iterations=%.3f max_pc_iterations=%d\n",
                   s, o.done, o.max_exceeded, o.rejected,
                   o.done + o.rejected ? (double)o.pc_total / (double)(o.done + o.rejected) : 0.0, o.pc_max);
    }
    printf("# steps_done=%lld iterations_max_exceeded=%lld steps_rejected=%lld mean_pc_iterations=%.3f max_pc_iterations=%d "
           "flags_checked=0x0f flags_raised=0x%02x wall_seconds=%.3f kernel_seconds=%.3f first_launch_seconds=%.3f launches=%lld E=%d "
           "system_steps=%lld system_steps_per_s=%.2f\n",
           out[0].done, mx_all, rej_all, att_all ? (double)pc_all / (double)att_all : 0.0, pcmax_all, fl_all, wall,
           t_kernel, t_first, launches, E, (long long)steps * E, wall > 0 ? (double)steps * E / wall : 0.0);
    fflush(stdout);
    if (att_dir){
        for (int s = 0; s < E; s++){
            /* the file is member<k>: system k of an ensemble, the --member k
             * run, 0 for one system (P3's names); the key is the one the
             * port's rerun of that member prints, which for a problem run
             * whole as one system has no --member, member=-1 */
            const int mk = member >= 0 ? member : s;
            const int keym = member >= 0 ? member : (E > 1 ? s : -1);
            const std::string key = run_key(P, K, keym, steps, sample, tol_shift, key_engine, dt0, &P.bodies[(size_t)s * N]);
            char path[1024], hx[40];
            snprintf(path, sizeof path, "%s/member%d.att", att_dir, mk);
            FILE *fa = fopen(path, "w");
            if (!fa){ perror(path); return 2; }
            fprintf(fa, "# cft-rebound attempts v1\n# run %s\n# a line a step attempt, in the order they were tried: "
                    "the dt it tried, then a where it was accepted or r where it was rejected; the last line is the dt "
                    "the run takes next, then next\n", key.c_str());
            for (const Attempt &a : log[s]){ wave_to_hex(a.dt, hx); fprintf(fa, "%s %c\n", hx, a.accepted ? 'a' : 'r'); }
            wave_to_hex(out[s].dt, hx); fprintf(fa, "%s next\n", hx);
            fclose(fa);
            snprintf(path, sizeof path, "%s/member%d.sched", att_dir, mk);
            FILE *fs = fopen(path, "w");
            if (!fs){ perror(path); return 2; }
            fprintf(fs, "# cft-rebound corrector schedule v1\n# run %s\n# a line a step attempt, a field a system: the "
                    "corrector passes it took, c where the cap stopped it, - where it was not in the attempt\n", key.c_str());
            for (const Attempt &a : log[s]) fprintf(fs, "%d%s\n", a.iters, a.capped ? "c" : "");
            fclose(fs);
            snprintf(path, sizeof path, "%s/member%d.key", att_dir, mk);
            FILE *fk = fopen(path, "w");
            if (fk){ fprintf(fk, "%s\n", key.c_str()); fclose(fk); }
        }
    }
    if (!quiet) fprintf(stderr, "wave: %ld steps a system, %d system(s), %.3f s on the %s, flags raised 0x%02x\n",
                        steps, E, wall, gpu ? "GPU" : "host", fl_all);
    return 0;
}
