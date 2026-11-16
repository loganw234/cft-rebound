/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * GATE 4 OF ROUND 5's P2: THE HARDWARE'S ARITHMETIC, MEASURED FIRST.
 *
 * Before anything rests on the hardware's add, subtract, multiply, fma,
 * divide and square root being libcft's correctly rounded binary64, this
 * runs a corpus through src/wave/wave_arith.h's operations - on the host
 * CPU (the twin's), and with --device gpu in a kernel on the GPU - and
 * holds every result to libcft's SOFTWARE binary64 (cft_open(NULL)), one
 * element a call so that each element's flags are its own:
 *   - the result's bits (a NaN is held as a NaN; its sign and payload,
 *     which no libcft contract fixes for arithmetic, are counted apart);
 *   - the four flags the integrator aborts on, INVALID, DIVBYZERO,
 *     OVERFLOW and UNDERFLOW, as wave_arith.h derives them from the
 *     operands and the result, against libcft's own;
 *   - augmentedAddition's two results and flags against
 *     cft_augmented_add, the exact-time clock's call.
 * On the host only, also the text: wave_to_hex() against cft_to_hex_char
 * and wave_from_hex() against cft_from_hex_char, and every decimal literal
 * the step uses against cft_from_decimal_char.
 *
 * The corpus is fixed by its seed: random bit patterns of every class;
 * random operands of like magnitude (cancellation) and of moderate
 * magnitude; ties - sums exactly halfway between two binary64 values,
 * where round-to-even and the clock's ties-toward-zero part, and products
 * exactly halfway; subnormal operands and results, exact and inexact
 * either side of 2^-1022 (tininess after rounding); results either side
 * of overflow, the clock's overflow tie included; and the specials
 * crossed with each other.
 *
 * Every mismatch is printed by name, and any one fails the gate:
 *
 *   wave_corpus [--device host|gpu] [--seed S] [--n N] [--quiet]
 */
#include "wave_arith.h"
#include "wave_io.h"
#include "cft.h"
#include <stdarg.h>
#include <string>
#include <vector>

using namespace wave;

static_assert((unsigned)F_INVALID == (unsigned)CFT_FLAG_INVALID && (unsigned)F_DIVBYZERO == (unsigned)CFT_FLAG_DIVBYZERO &&
              (unsigned)F_OVERFLOW == (unsigned)CFT_FLAG_OVERFLOW && (unsigned)F_UNDERFLOW == (unsigned)CFT_FLAG_UNDERFLOW,
              "wave_arith.h's flag bits are not cft.h's");

enum { OP_ADD, OP_SUB, OP_MUL, OP_DIV, OP_SQRT, OP_FMA, OP_AUG, NOPS };
static const char *const op_name[NOPS] = { "add", "sub", "mul", "div", "sqrt", "fma", "augmented_add" };
static const unsigned CHECKED = F_INVALID | F_DIVBYZERO | F_OVERFLOW | F_UNDERFLOW;

struct Case { double a, b, c; };
struct Res { double r, e; unsigned fl; };

/* The one evaluator both sides run: the operation through wave_arith.h */
__host__ __device__ inline Res eval_one(int op, const Case &k){
    Fl fl; fl.clear();
    Res o; o.e = 0.0;
    switch (op){
    case OP_ADD:  o.r = w_add(k.a, k.b, fl); break;
    case OP_SUB:  o.r = w_sub(k.a, k.b, fl); break;
    case OP_MUL:  o.r = w_mul(k.a, k.b, fl); break;
    case OP_DIV:  o.r = w_div(k.a, k.b, fl); break;
    case OP_SQRT: o.r = w_sqrt(k.a, fl); break;
    case OP_FMA:  o.r = w_fma(k.a, k.b, k.c, fl); break;
    default:      w_augadd(k.a, k.b, o.r, o.e, fl); break;
    }
    o.fl = fl.bits;
    return o;
}

#ifdef WAVE_HIP
__global__ void corpus_kernel(int op, const Case *cs, Res *out, size_t n){
    const size_t i = (size_t)blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = eval_one(op, cs[i]);
}
#define HIPCHECK(x) do { hipError_t e_ = (x); if (e_ != hipSuccess){ fprintf(stderr, "wave_corpus: %s: %s\n", #x, hipGetErrorString(e_)); exit(2); } } while (0)
static std::string gpu_name;
static void eval_gpu(int op, const std::vector<Case> &cs, std::vector<Res> &out){
    Case *dc; Res *dr;
    const size_t n = cs.size();
    HIPCHECK(hipMalloc((void **)&dc, n * sizeof(Case)));
    HIPCHECK(hipMalloc((void **)&dr, n * sizeof(Res)));
    HIPCHECK(hipMemcpy(dc, cs.data(), n * sizeof(Case), hipMemcpyHostToDevice));
    const unsigned threads = 128;
    hipLaunchKernelGGL(corpus_kernel, dim3((unsigned)((n + threads - 1) / threads)), dim3(threads), 0, 0, op, dc, dr, n);
    HIPCHECK(hipGetLastError());
    HIPCHECK(hipDeviceSynchronize());
    out.resize(n);
    HIPCHECK(hipMemcpy(out.data(), dr, n * sizeof(Res), hipMemcpyDeviceToHost));
    HIPCHECK(hipFree(dc)); HIPCHECK(hipFree(dr));
}
#endif

static void eval_host(int op, const std::vector<Case> &cs, std::vector<Res> &out){
    out.resize(cs.size());
    for (size_t i = 0; i < cs.size(); i++) out[i] = eval_one(op, cs[i]);
}

/* ---- the reference: libcft's software binary64, one element a call ---- */
static cft_device *sw;
static Res eval_ref(int op, const Case &k){
    Res o; o.e = 0.0;
    uint32_t fl = 0;
    cft_status st = CFT_OK;
    switch (op){
    case OP_ADD:  st = cft_run(sw, CFT_ADD, CFT_FP64, CFT_RNE, &k.a, NULL, &k.b, &o.r, 1, &fl, NULL); break;
    case OP_SUB:  st = cft_run(sw, CFT_SUB, CFT_FP64, CFT_RNE, &k.a, NULL, &k.b, &o.r, 1, &fl, NULL); break;
    case OP_MUL:  st = cft_run(sw, CFT_MUL, CFT_FP64, CFT_RNE, &k.a, &k.b, NULL, &o.r, 1, &fl, NULL); break;
    case OP_DIV:  st = cft_div(sw, CFT_FP64, CFT_RNE, &k.a, &k.b, &o.r, 1, &fl, NULL); break;
    case OP_SQRT: st = cft_sqrt(sw, CFT_FP64, CFT_RNE, &k.a, &o.r, 1, &fl, NULL); break;
    case OP_FMA:  st = cft_run(sw, CFT_FMA, CFT_FP64, CFT_RNE, &k.a, &k.b, &k.c, &o.r, 1, &fl, NULL); break;
    default:      st = cft_augmented_add(sw, CFT_FP64, &k.a, &k.b, &o.r, &o.e, 1, &fl); break;
    }
    if (st != CFT_OK){ fprintf(stderr, "wave_corpus: libcft %s refused: %s (%s)\n", op_name[op], cft_strerror(st), cft_last_error()); exit(2); }
    o.fl = fl & CHECKED;
    return o;
}

/* ---- the corpus ------------------------------------------------------------ */
struct Rng {
    uint64_t s;
    uint64_t next(){ s += 0x9e3779b97f4a7c15ull; uint64_t z = s; z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ull;
                     z = (z ^ (z >> 27)) * 0x94d049bb133111ebull; return z ^ (z >> 31); }    /* splitmix64 */
    double unit(){ return (double)(next() >> 11) * 0x1p-53; }                               /* [0, 1), exact */
    int range(int lo, int hi){ return lo + (int)(next() % (uint64_t)(hi - lo + 1)); }
    /* a normal of random sign, significand and exponent in [lo, hi] */
    double normal(int lo, int hi){
        const uint64_t frac = next() & 0x000fffffffffffffull;
        const int e = range(lo, hi);
        return u2d((next() & 1 ? 0x8000000000000000ull : 0) | ((uint64_t)(e + 1023) << 52) | frac);
    }
    double any(){ return u2d(next()); }
    double sub(){ return u2d((next() & 0x800fffffffffffffull) | 1); }                       /* a subnormal, never 0 */
};

static std::vector<double> specials(){
    const double v[] = { 0.0, 0x1p-1074, 0x1.8p-1073, 0x0.fffffffffffffp-1022, 0x1p-1022, 0x1.0000000000001p-1022,
                         0x1.fffffffffffffp-1022, 0x1p-537, 0x1p-511, 0x1p-500, 0.5, 1.0, 0x1.0000000000001p+0,
                         0x1.fffffffffffffp+0, 3.0, 0x1p511, 0x1p970, 0x1p1023, 0x1.fffffffffffffp+1023, 1e-7,
                         0.1, 7.0, 1e300, 1e-300 };
    std::vector<double> s;
    for (double x : v){ s.push_back(x); s.push_back(-x); }
    s.push_back(HUGE_VAL); s.push_back(-HUGE_VAL);
    s.push_back(u2d(0x7ff8000000000000ull));
    return s;
}

static void add_case(std::vector<Case> &c, double a, double b, double d = 0.0){ c.push_back(Case{ a, b, d }); }

/* a random binary64 x with an odd last bit and |x| in binade e, and a b
 * that puts x + b exactly halfway between x and a neighbour */
static void ties_add(Rng &R, std::vector<Case> &c, int n){
    for (int i = 0; i < n; i++){
        const int e = R.range(-1000, 1000);
        uint64_t M = (1ull << 52) | (R.next() & 0x000fffffffffffffull) | 1ull;     /* odd */
        if (R.next() & 1) M &= ~1ull;                                                /* sometimes even */
        const double x = ldexp((double)M, e - 52) * (R.next() & 1 ? -1.0 : 1.0);
        const uint64_t k = R.next() % (1ull << 20);
        double b = ldexp((double)(2 * k + 1), e - 53);                               /* an odd multiple of half an ulp */
        if (R.next() & 1) b = -b;
        if (R.next() & 3) b = w_copysign(ldexp(1.0, e - 53), R.next() & 1 ? x : -x); /* exactly half an ulp */
        add_case(c, x, b);
        add_case(c, b, x);
    }
}

/* products exactly halfway: an odd 53-bit M with 3M of 54 bits, times 3 */
static void ties_mul(Rng &R, std::vector<Case> &c, int n){
    for (int i = 0; i < n; i++){
        uint64_t M;
        do { M = (1ull << 52) | (R.next() & 0x000fffffffffffffull) | 1ull; } while (3 * M >= (1ull << 54));
        const int e = R.range(-600, 600), f = R.range(-400, 400);
        const double a = ldexp((double)M, e - 52), b = ldexp(3.0, f);
        add_case(c, R.next() & 1 ? -a : a, b);
    }
}

static std::vector<Case> corpus(int op, Rng &R, int n){
    std::vector<Case> c;
    const std::vector<double> sp = specials();
    /* the specials crossed */
    if (op == OP_SQRT){ for (double a : sp) add_case(c, a, 0.0); }
    else if (op == OP_FMA){
        for (double a : sp) for (double b : sp) add_case(c, a, b, sp[(c.size() * 7) % sp.size()]);
        for (double a : sp) for (double b : sp) add_case(c, a, 1.0, b);
    }else for (double a : sp) for (double b : sp) add_case(c, a, b);
    for (int i = 0; i < n; i++){
        switch (op){
        case OP_ADD: case OP_SUB: case OP_AUG: {
            add_case(c, R.any(), R.any());
            add_case(c, R.normal(-30, 30), R.normal(-30, 30));
            const double a = R.normal(-1022, 1023);
            add_case(c, a, w_neg(a) * (1.0 + R.unit() * 0x1p-20));                  /* cancellation */
            add_case(c, R.sub(), R.sub());
            add_case(c, R.sub(), R.normal(-1022, -1000));
            add_case(c, R.normal(1015, 1023), R.normal(1015, 1023));                /* near overflow */
            add_case(c, R.normal(-1022, 1023), R.normal(-1022, 1023));
            break; }
        case OP_MUL: {
            add_case(c, R.any(), R.any());
            add_case(c, R.normal(-30, 30), R.normal(-30, 30));
            add_case(c, R.normal(-600, -400), R.normal(-700, -500));                /* subnormal products */
            add_case(c, R.normal(-560, -500), R.normal(-560, -500));
            const double t = 1.0 + R.unit();                                        /* products about 2^-1022 */
            add_case(c, ldexp(t, -511), ldexp(1.0 / t, -511));
            add_case(c, ldexp(1.0, R.range(-1074, -1000)), (double)(R.next() % (1u << 20)) * 0x1p-40);   /* exact tiny */
            add_case(c, R.sub(), R.normal(-60, 60));
            add_case(c, R.normal(500, 530), R.normal(490, 530));                    /* near overflow */
            break; }
        case OP_DIV: {
            add_case(c, R.any(), R.any());
            add_case(c, R.normal(-30, 30), R.normal(-30, 30));
            add_case(c, R.normal(-1000, -900), R.normal(100, 200));                 /* subnormal quotients */
            const double t = 1.0 + R.unit();
            add_case(c, ldexp(t, -1022), t * (1.0 + R.unit() * 0x1p-40));           /* about 2^-1022 */
            add_case(c, ldexp((double)(R.next() % (1u << 20)), -1074), ldexp(1.0, R.range(0, 30)));   /* exact tiny */
            add_case(c, R.sub(), R.normal(-20, 20));
            add_case(c, R.normal(1000, 1023), R.normal(-100, -1));                  /* overflow */
            add_case(c, R.normal(-30, 30), R.next() & 1 ? 0.0 : -0.0);              /* by zero */
            break; }
        case OP_SQRT: {
            add_case(c, R.any(), 0.0);
            add_case(c, w_abs(R.normal(-1022, 1023)), 0.0);
            add_case(c, w_abs(R.sub()), 0.0);
            const double q = (double)(R.next() % (1u << 26));
            add_case(c, q * q, 0.0);                                                /* exact squares */
            add_case(c, w_neg(w_abs(R.normal(-30, 30))), 0.0);
            break; }
        case OP_FMA: {
            add_case(c, R.any(), R.any(), R.any());
            add_case(c, R.normal(-30, 30), R.normal(-30, 30), R.normal(-60, 60));
            const double a = R.normal(-500, -470), b = R.normal(-500, -470);
            const double p = a * b;                                                 /* about 2^-1000..2^-940 */
            add_case(c, a, b, w_neg(p));                                            /* the product's error: tiny */
            add_case(c, a, b, w_neg(p) + R.sub());
            add_case(c, R.normal(-560, -500), R.normal(-560, -500), R.sub());
            add_case(c, R.normal(-540, -530), R.normal(-540, -530), R.normal(-1022, -1000));
            add_case(c, R.normal(500, 530), R.normal(490, 530), R.normal(1000, 1023));
            const double x = R.normal(-60, 60), y = R.normal(-60, 60);
            add_case(c, x, y, w_neg(x * y));                                        /* exact cancellations */
            break; }
        }
    }
    if (op == OP_ADD || op == OP_SUB || op == OP_AUG || op == OP_FMA){
        std::vector<Case> t;
        ties_add(R, t, n / 4 + 64);
        for (const Case &k : t){
            if (op == OP_FMA) add_case(c, k.a, 1.0, k.b);
            else if (op == OP_SUB) add_case(c, k.a, w_neg(k.b));
            else add_case(c, k.a, k.b);
        }
    }
    if (op == OP_FMA){
        /* fma's underflow at its edge: a*b + c within a few units of
         * 2^-1022 - 2^-1076, the largest tiny value, with the product's
         * bits below 2^-1078 (a sticky remainder) or not, either sign:
         * a = j (1 + u 2^-52) 2^-538, b = 2^-538, c = k 2^-1074 */
        for (int64_t k = (1ll << 52) - 3; k <= (1ll << 52) + 3; k++)
            for (int j = 1; j <= 8; j++)
                for (int u = 0; u < 2; u++)
                    for (int sg = 0; sg < 4; sg++){
                        const double a = (sg & 1 ? -1.0 : 1.0) * ldexp((double)j * (1.0 + u * 0x1p-52), -538);
                        const double cc = (sg & 2 ? -1.0 : 1.0) * ldexp((double)k, -1074);
                        add_case(c, a, 0x1p-538, cc);
                    }
    }
    if (op == OP_MUL || op == OP_FMA){
        std::vector<Case> t;
        ties_mul(R, t, n / 8 + 64);
        for (const Case &k : t){ if (op == OP_FMA) add_case(c, k.a, k.b, 0.0); else add_case(c, k.a, k.b); }
    }
    if (op == OP_AUG){
        const double M = u2d(0x7fefffffffffffffull);
        add_case(c, M, 0x1p970); add_case(c, -M, -0x1p970);                         /* the RTTZ overflow tie */
        add_case(c, M, 0x1p970 + 0x1p918); add_case(c, M, 0x1.fffffffffffffp969);
        add_case(c, 0x1p1023, 0x1.fffffffffffffp1022); add_case(c, M, M);
    }
    return c;
}

/* ---- comparison ------------------------------------------------------------ */
static int same_bits(double x, double y){
    if (is_nan(x) && is_nan(y)) return 1;
    return d2u(x) == d2u(y);
}
static int nan_pattern_differs(double x, double y){ return is_nan(x) && is_nan(y) && d2u(x) != d2u(y); }

static int quiet_mode;
static unsigned long long printed;
static void report(const char *side, int op, const Case &k, const Res &ref, const Res &got, const char *what){
    char a[40], b[40], c[40], r1[40], r2[40], e1[40], e2[40];
    if (printed++ > 2000) return;
    wave_to_hex(k.a, a); wave_to_hex(k.b, b); wave_to_hex(k.c, c);
    wave_to_hex(ref.r, r1); wave_to_hex(got.r, r2); wave_to_hex(ref.e, e1); wave_to_hex(got.e, e2);
    printf("MISMATCH %s %s (%s): a=%s b=%s c=%s; libcft r=%s e=%s flags=0x%02x; %s r=%s e=%s flags=0x%02x\n",
           side, op_name[op], what, a, b, c, r1, e1, ref.fl, side, r2, e2, got.fl);
}

/* ---- the text, host only ---------------------------------------------------- */
static int text_checks(Rng &R, int n){
    int bad = 0;
    unsigned long long nfmt = 0, nparse = 0, nlit = 0;
    std::vector<double> vals = specials();
    for (int i = 0; i < n; i++){ vals.push_back(R.any()); vals.push_back(R.sub()); vals.push_back(R.normal(-1022, 1023)); }
    for (double x : vals){
        char mine[40], ref[160]; size_t len = 0;
        wave_to_hex(x, mine);
        if (cft_to_hex_char(sw, CFT_FP64, &x, ref, sizeof ref, &len) != CFT_OK){ printf("MISMATCH hex: libcft refused\n"); bad++; continue; }
        nfmt++;
        if (strcmp(mine, ref)){ printf("MISMATCH to_hex %016llx: libcft %s, wave %s\n", (unsigned long long)d2u(x), ref, mine); bad++; }
        if (!is_nan(x)){
            double back;
            nparse++;
            if (!wave_from_hex(ref, &back) || d2u(back) != d2u(x)){ printf("MISMATCH from_hex %s: wave gives %016llx\n", ref, (unsigned long long)d2u(back)); bad++; }
        }
    }
    /* sequences with more bits than binary64 holds, rounded */
    std::vector<std::string> seq = { "0x1.00000000000008p+0", "0x1.00000000000018p+0", "0x1.000000000000080000001p+0",
        "0x1.fffffffffffff8p+0", "0x1.fffffffffffff7ffffp+1023", "0x1.fffffffffffff8p+1023", "0x1p-1075", "0x1.0000001p-1075",
        "0x1.8p-1074", "0x1.4p-1074", "0x0.0000000000001p-1022", "0x0.00000000000008p-1022", "0x10p-4", "0x.8p1",
        "0x000000000000000000001p+0", "0x1.23456789abcdef0123456789p+100", "-0x0p+0", "0x0.0p+0", "-0x1.8p-1074" };
    for (int i = 0; i < n; i++){
        char b[96];
        snprintf(b, sizeof b, "%s0x%llx.%llxp%d", R.next() & 1 ? "-" : "", (unsigned long long)(R.next() >> R.range(0, 63)),
                 (unsigned long long)R.next(), R.range(-1200, 1100));
        seq.push_back(b);
    }
    for (const std::string &s : seq){
        double mine = 0, ref = 0; size_t bi = 0; uint32_t fl = 0;
        const char *in[1] = { s.c_str() };
        const cft_status st = cft_from_hex_char(sw, CFT_FP64, CFT_RNE, in, &ref, 1, &bi, &fl);
        const int ok = wave_from_hex(s.c_str(), &mine);
        nparse++;
        if ((st == CFT_OK) != (ok != 0) || (ok && d2u(mine) != d2u(ref))){
            printf("MISMATCH from_hex \"%s\": libcft %s %016llx, wave %s %016llx\n", s.c_str(), st == CFT_OK ? "ok" : "refused",
                   (unsigned long long)d2u(ref), ok ? "ok" : "refused", (unsigned long long)d2u(mine));
            bad++;
        }
    }
    /* the step's decimal literals (src/ias15_cft.c make_constants() and
     * main()'s defaults): C++'s, as the compiler rounded them */
    struct { const char *text; double v; } lit[] = { { "0.1", 0.1 }, { "0.25", 0.25 }, { "1e2", 1e2 }, { "1e7", 1e7 },
        { "1e-7", 1e-7 }, { "1e300", 1e300 }, { "1e-16", 1e-16 }, { "5040", 5040.0 }, { "20", 20.0 }, { "10", 10.0 },
        { "7", 7.0 }, { "2", 2.0 }, { "0.01", 0.01 }, { "1e-9", 1e-9 } };
    for (auto &l : lit){
        double ref = 0; size_t bi = 0; uint32_t fl = 0;
        const char *in[1] = { l.text };
        nlit++;
        if (cft_from_decimal_char(sw, CFT_FP64, CFT_RNE, in, &ref, 1, &bi, &fl) != CFT_OK || d2u(ref) != d2u(l.v)){
            printf("MISMATCH literal %s: libcft %016llx, C++ %016llx\n", l.text, (unsigned long long)d2u(ref), (unsigned long long)d2u(l.v));
            bad++;
        }
        const double s = strtod(l.text, NULL);
        if (d2u(s) != d2u(ref)){ printf("MISMATCH strtod %s: libcft %016llx, strtod %016llx\n", l.text, (unsigned long long)d2u(ref), (unsigned long long)d2u(s)); bad++; }
    }
    printf("text: %llu values to hex against cft_to_hex_char, %llu sequences from hex against cft_from_hex_char, "
           "%llu decimal literals (and strtod) against cft_from_decimal_char: %d mismatches\n", nfmt, nparse, nlit, bad);
    return bad;
}

int main(int argc, char **argv){
    const char *device = "host";
    uint64_t seed = 20261116;
    int n = 20000;
    for (int i = 1; i < argc; i++){
        if (!strcmp(argv[i], "--device") && i + 1 < argc) device = argv[++i];
        else if (!strcmp(argv[i], "--seed") && i + 1 < argc) seed = strtoull(argv[++i], NULL, 10);
        else if (!strcmp(argv[i], "--n") && i + 1 < argc) n = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--quiet")) quiet_mode = 1;
        else { fprintf(stderr, "wave_corpus: unknown argument %s\n", argv[i]); return 2; }
    }
    const int gpu = !strcmp(device, "gpu");
#ifndef WAVE_HIP
    if (gpu){ fprintf(stderr, "wave_corpus: --device gpu: this build has no GPU (make HIP=1 wave)\n"); return 2; }
#else
    if (gpu){
        hipDeviceProp_t pr; int dv = 0;
        HIPCHECK(hipGetDevice(&dv));
        HIPCHECK(hipGetDeviceProperties(&pr, dv));
        char line[512];
#if defined(__HIP_PLATFORM_AMD__)
        const char *arch = pr.gcnArchName;
#else
        char archb[32]; snprintf(archb, sizeof archb, "sm_%d%d", pr.major, pr.minor);
        const char *arch = archb;
#endif
        /* the device as the kernel's record names it (wave_gpu.inc, gpu_open()) */
        snprintf(line, sizeof line, "hip:%d = %s (%s, %d multiprocessors, bus %04x:%02x:%02x.0, %zu MiB, warp %d)", dv,
                 pr.name, arch, pr.multiProcessorCount, pr.pciDomainID, pr.pciBusID, pr.pciDeviceID,
                 (size_t)(pr.totalGlobalMem >> 20), pr.warpSize);
        gpu_name = line;
        printf("device: %s\n", line);
    }
#endif
    if (cft_open(NULL, 0, &sw) != CFT_OK){ fprintf(stderr, "wave_corpus: cft_open(software): %s\n", cft_last_error()); return 2; }
    printf("wave_corpus: src/wave/wave_arith.h on the %s against libcft's software binary64 (cft-fp256 at the pin), "
           "seed %llu, n %d\n", gpu ? "GPU" : "host CPU", (unsigned long long)seed, n);
    Rng R; R.s = seed;
    unsigned long long total = 0, bad = 0;
    for (int op = 0; op < NOPS; op++){
        Rng Ro; Ro.s = seed * 1000003ull + (uint64_t)op;
        const std::vector<Case> cs = corpus(op, Ro, n);
        std::vector<Res> got;
#ifdef WAVE_HIP
        if (gpu) eval_gpu(op, cs, got); else
#endif
        eval_host(op, cs, got);
        unsigned long long nb = 0, nf = 0, nn = 0, flagged[4] = { 0, 0, 0, 0 }, sub = 0, ties = 0;
        for (size_t i = 0; i < cs.size(); i++){
            const Res ref = eval_ref(op, cs[i]);
            const Res &g = got[i];
            for (int f = 0; f < 4; f++) flagged[f] += (ref.fl >> f) & 1u;
            if (!is_zero(ref.r) && w_abs(ref.r) < DBL_MIN) sub++;
            if (op == OP_AUG && !is_zero(ref.e) && is_finite(ref.e)) ties += (d2u(ref.r) != d2u(cs[i].a + cs[i].b));
            const int rb = same_bits(ref.r, g.r) && (op != OP_AUG || same_bits(ref.e, g.e));
            if (!rb){ nb++; report(gpu ? "gpu" : "host", op, cs[i], ref, g, "result bits"); }
            if (ref.fl != (g.fl & CHECKED)){ nf++; report(gpu ? "gpu" : "host", op, cs[i], ref, g, "flags"); }
            if (nan_pattern_differs(ref.r, g.r)) nn++;
        }
        total += cs.size(); bad += nb + nf;
        printf("%-14s %8zu cases: %llu result mismatches, %llu flag mismatches; libcft raised INVALID %llu, DIVBYZERO %llu, "
               "OVERFLOW %llu, UNDERFLOW %llu; %llu subnormal results%s; %llu NaN results of another sign or payload "
               "(counted, not compared)\n",
               op_name[op], cs.size(), nb, nf, flagged[0], flagged[1], flagged[2], flagged[3], sub,
               op == OP_AUG ? (std::string("; ") + std::to_string(ties) + " ties where round-to-even and the clock's ties-toward-zero differ").c_str() : "",
               nn);
    }
    int tbad = 0;
    if (!gpu) tbad = text_checks(R, n / 4);
    printf("RESULT: %s: %llu cases, %llu mismatches%s\n", bad || tbad ? "FAIL" : "PASS", total, bad,
           gpu ? "" : (std::string(", text ") + std::to_string(tbad)).c_str());
    cft_close(sw);
    return bad || tbad ? 1 : 0;
}
