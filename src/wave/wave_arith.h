/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The step's arithmetic: binary64 add, subtract, multiply, divide, square
 * root and fused multiply-add, each the hardware's own correctly rounded
 * operation, round to nearest even - and IEEE 754-2019's augmentedAddition
 * (9.5) for the exact-time clock. Gate 4 (wave_corpus.cpp) measures every
 * one against libcft's software binary64 before anything rests on it.
 *
 * FLAGS. A GPU has no IEEE status flags, and the port aborts when a call
 * raises INVALID, DIVBYZERO, OVERFLOW or UNDERFLOW (src/ias15_cft.c,
 * note(); --no-flag-abort turns that off). So each operation here decides
 * from its operands and its result which of those four libcft would have
 * raised, by libcft's own definitions (cft-fp256 host/src/softfloat.c and
 * augmented.c, read at the pin):
 *   - INVALID: a NaN from operands that hold none (0/0, inf-inf, 0*inf,
 *     the square root of a negative);
 *   - DIVBYZERO: a finite non-zero numerator over a zero;
 *   - OVERFLOW: an infinity from finite operands;
 *   - UNDERFLOW: TINY AFTER ROUNDING AND INEXACT. Tiny means the result
 *     rounded to 53 bits with the exponent unbounded lies below 2^-1022;
 *     inexact means the exact result is not the delivered one. A result
 *     rounded to zero from a non-zero exact value is both.
 * The common case costs one comparison: a result that is finite and
 * larger in magnitude than 2^-1022 raised none of the four, whatever its
 * operands. Only the rest goes to the exact analysis below, which is done
 * in integers and in correctly rounded binary64 operations that cannot
 * themselves underflow, never by asking the hardware.
 *
 * NaN operands: a quiet NaN raises nothing and a signaling one INVALID,
 * as in libcft. (The hardware's NaN results may differ from libcft's in
 * sign and payload; no contract fixes those, and a NaN reaches a record
 * only past an INVALID, which stops the run by default.)
 *
 * Each operation takes the thread's Fl: the union of what was raised, and
 * the first raise by phase and operation, so that a system that stops
 * says where.
 */
#ifndef WAVE_ARITH_H
#define WAVE_ARITH_H

#include "wave_hd.h"

namespace wave {

/* libcft's bit values (cft.h's CFT_FLAG_*): wave_corpus.cpp, which links
 * libcft, asserts each against cft.h's own. */
enum : unsigned { F_INVALID = 1u, F_DIVBYZERO = 2u, F_OVERFLOW = 4u, F_UNDERFLOW = 8u };

enum Kind : unsigned { K_ADD = 1, K_SUB, K_MUL, K_DIV, K_SQRT, K_FMA, K_AUG };

struct Fl {
    unsigned bits;      /* every flag raised */
    unsigned first;     /* phase << 12 | kind << 4 | flags, of the first raise; 0 for none */
    unsigned phase;     /* where the step is, set by the step as it goes */
    __host__ __device__ WAVE_INLINE void clear(){ bits = 0; first = 0; phase = 0; }
    __host__ __device__ WAVE_INLINE void raise(unsigned f, unsigned kind){
        if (!f) return;
        if (!first) first = (phase << 12) | (kind << 4) | f;
        bits |= f;
    }
};

/* ---- the two operations a compiler might hand to a library ------------
 * On an x86-64 host, fma() and sqrt() are what the C runtime defines them
 * as unless the compiler inlines them, and it does so only when it
 * optimises: MEASURED with MinGW's g++ 16.1 at -O0, three calls to
 * libmingwex's fma, a software emulation, and a ring record that differed
 * from the port's. So the host takes the instructions themselves through
 * their intrinsics, vfmadd and sqrtsd, at every optimisation level; the
 * device takes its own fma and sqrt, which gate 4 holds to libcft. */
#if !WAVE_DEVICE_PASS && (defined(__x86_64__) || defined(_M_X64))
}   /* namespace wave */
#  include <immintrin.h>
namespace wave {
__host__ __device__ WAVE_INLINE double hw_fma(double a, double b, double c){
    return _mm_cvtsd_f64(_mm_fmadd_sd(_mm_set_sd(a), _mm_set_sd(b), _mm_set_sd(c)));
}
__host__ __device__ WAVE_INLINE double hw_sqrt(double a){
    const __m128d v = _mm_set_sd(a);
    return _mm_cvtsd_f64(_mm_sqrt_sd(v, v));
}
#else
__host__ __device__ WAVE_INLINE double hw_fma(double a, double b, double c){ return fma(a, b, c); }
__host__ __device__ WAVE_INLINE double hw_sqrt(double a){ return sqrt(a); }
#endif

/* ---- bits ------------------------------------------------------------ */
__host__ __device__ WAVE_INLINE uint64_t d2u(double x){ uint64_t u; memcpy(&u, &x, 8); return u; }
__host__ __device__ WAVE_INLINE double u2d(uint64_t u){ double x; memcpy(&x, &u, 8); return x; }
__host__ __device__ WAVE_INLINE int is_nan(double x){ return (d2u(x) & 0x7fffffffffffffffull) > 0x7ff0000000000000ull; }
/* a signaling NaN: a NaN whose quiet bit (the fraction's top bit) is clear */
__host__ __device__ WAVE_INLINE int is_snan(double x){
    const uint64_t u = d2u(x) & 0x7fffffffffffffffull;
    return u > 0x7ff0000000000000ull && !(u & 0x0008000000000000ull);
}
__host__ __device__ WAVE_INLINE int is_inf(double x){ return (d2u(x) & 0x7fffffffffffffffull) == 0x7ff0000000000000ull; }
__host__ __device__ WAVE_INLINE int is_zero(double x){ return (d2u(x) & 0x7fffffffffffffffull) == 0; }
__host__ __device__ WAVE_INLINE int is_finite(double x){ return (d2u(x) & 0x7ff0000000000000ull) != 0x7ff0000000000000ull; }
__host__ __device__ WAVE_INLINE int sign_of(double x){ return (int)(d2u(x) >> 63); }
/* libcft's CFT_CLASS_*_NORM: neither zero, subnormal, infinite nor NaN */
__host__ __device__ WAVE_INLINE int is_normal(double x){
    uint64_t e = (d2u(x) >> 52) & 0x7ff;
    return e != 0 && e != 0x7ff;
}
/* |x|, a sign bit cleared: libcft's CFT_ABS, which raises nothing */
__host__ __device__ WAVE_INLINE double w_abs(double x){ return u2d(d2u(x) & 0x7fffffffffffffffull); }
/* -x, the sign bit flipped: CFT_NEG, exact and raising nothing */
__host__ __device__ WAVE_INLINE double w_neg(double x){ return u2d(d2u(x) ^ 0x8000000000000000ull); }
/* CFT_COPYSIGN: a's magnitude with b's sign */
__host__ __device__ WAVE_INLINE double w_copysign(double a, double b){
    return u2d((d2u(a) & 0x7fffffffffffffffull) | (d2u(b) & 0x8000000000000000ull));
}

/* A finite non-zero x as M * 2^E with M an integer of at most 53 bits:
 * the encoding's own fields (a subnormal's M is its fraction, E -1074). */
__host__ __device__ WAVE_INLINE void decomp(double x, uint64_t &M, int &E){
    uint64_t u = d2u(x);
    int be = (int)((u >> 52) & 0x7ff);
    M = u & 0x000fffffffffffffull;
    if (be){ M |= 0x0010000000000000ull; E = be - 1075; }
    else     E = -1074;
}
__host__ __device__ WAVE_INLINE int bitlen64(uint64_t v){ int n = 0; while (v){ n++; v >>= 1; } return n; }

/* ---- 128-bit integers, for the one check that needs them ------------- */
struct U128 { uint64_t hi, lo; };
__host__ __device__ WAVE_INLINE U128 u128_mul(uint64_t a, uint64_t b){
    const uint64_t a0 = a & 0xffffffffu, a1 = a >> 32, b0 = b & 0xffffffffu, b1 = b >> 32;
    const uint64_t p00 = a0 * b0, p01 = a0 * b1, p10 = a1 * b0, p11 = a1 * b1;
    const uint64_t mid = (p00 >> 32) + (p01 & 0xffffffffu) + (p10 & 0xffffffffu);
    U128 r;
    r.lo = (mid << 32) | (p00 & 0xffffffffu);
    r.hi = p11 + (p01 >> 32) + (p10 >> 32) + (mid >> 32);
    return r;
}
__host__ __device__ WAVE_INLINE int u128_cmp(U128 a, U128 b){
    if (a.hi != b.hi) return a.hi < b.hi ? -1 : 1;
    if (a.lo != b.lo) return a.lo < b.lo ? -1 : 1;
    return 0;
}
__host__ __device__ WAVE_INLINE U128 u128_add(U128 a, U128 b){
    U128 r; r.lo = a.lo + b.lo; r.hi = a.hi + b.hi + (r.lo < a.lo); return r;
}
__host__ __device__ WAVE_INLINE U128 u128_sub(U128 a, U128 b){   /* a >= b */
    U128 r; r.lo = a.lo - b.lo; r.hi = a.hi - b.hi - (a.lo < b.lo); return r;
}
/* v << s, for 0 <= s and a result that fits */
__host__ __device__ WAVE_INLINE U128 u128_shl(U128 v, int s){
    U128 r;
    if (s <= 0) return v;
    if (s >= 128){ r.hi = r.lo = 0; return r; }
    if (s >= 64){ r.hi = v.lo << (s - 64); r.lo = 0; return r; }
    r.hi = (v.hi << s) | (v.lo >> (64 - s)); r.lo = v.lo << s;
    return r;
}
/* v >> s, s >= 0, with *sticky set when a one bit is shifted out */
__host__ __device__ WAVE_INLINE U128 u128_shr(U128 v, int s, int *sticky){
    U128 r;
    if (s <= 0) return v;
    if (s >= 128){ *sticky |= (v.hi | v.lo) != 0; r.hi = r.lo = 0; return r; }
    if (s >= 64){
        const int t = s - 64;
        *sticky |= v.lo != 0 || (t && (v.hi << (64 - t)) != 0);
        r.lo = t ? v.hi >> t : v.hi; r.hi = 0;
        return r;
    }
    *sticky |= (v.lo << (64 - s)) != 0;
    r.lo = (v.lo >> s) | (v.hi << (64 - s)); r.hi = v.hi >> s;
    return r;
}

/* ---- the underflow decision ------------------------------------------ */
/* For a result whose exact value is M * 2^S with 2^52 <= M < 2^53 its
 * correct 53-bit rounding (exponent unbounded): tiny iff that lies below
 * 2^-1022, i.e. S <= -1075. Inexact on the binary64 grid unless the exact
 * value was exactly this one (exact53) and no bit of M lies below 2^-1074. */
__host__ __device__ WAVE_INLINE int uf_from53(uint64_t M, int S, int exact53){
    if (S > -1075) return 0;                          /* not tiny */
    if (!exact53) return 1;                           /* tiny and inexact */
    const int k = -1074 - S;                          /* bits of M below the subnormal quantum */
    if (k <= 0) return 0;
    if (k >= 53) return 1;
    return (M & ((1ull << k) - 1)) != 0;
}

/* a*b, both finite and non-zero, its rounded result not normal: whether
 * libcft raises UNDERFLOW. fa and fb are the integer significands as
 * doubles (exact); their product is a normal number and its rounding error
 * is exact, so nothing here underflows. */
__host__ __device__ WAVE_INLINE int uf_mul(double a, double b){
    uint64_t Ma, Mb, Mh; int Ea, Eb, Eh;
    decomp(a, Ma, Ea); decomp(b, Mb, Eb);
    const double fa = (double)Ma, fb = (double)Mb;
    const double hi = fa * fb;
    const double lo = hw_fma(fa, fb, -hi);
    decomp(hi, Mh, Eh);                               /* hi >= 1 is normal: 2^52 <= Mh < 2^53 */
    return uf_from53(Mh, Eh + Ea + Eb, lo == 0.0);
}

/* a/b, both finite and non-zero: the same question. q = fa/fb lies in
 * [2^-53, 2^53], normal, correctly rounded; its remainder fa - q fb is
 * exact for a correctly rounded quotient, and zero iff the quotient was. */
__host__ __device__ WAVE_INLINE int uf_div(double a, double b){
    uint64_t Ma, Mb, Mq; int Ea, Eb, Eq;
    decomp(a, Ma, Ea); decomp(b, Mb, Eb);
    const double fa = (double)Ma, fb = (double)Mb;
    const double q = fa / fb;
    const double rem = hw_fma(-q, fb, fa);
    decomp(q, Mq, Eq);
    return uf_from53(Mq, Eq + Ea - Eb, rem == 0.0);
}

/* fma(a, b, c), all finite, a and b non-zero, c non-zero, its rounded
 * result not normal (|r| <= 2^-1022). When either term reaches 2^-966 the
 * small result needs a cancellation, and then both terms' quanta are at or
 * above 2^-1074, so the sum is exact and nothing is raised. Below that
 * both terms fit a 128-bit window whose unit is 2^-1078, where the sum is
 * formed exactly (the product's bits below the unit kept as a sticky
 * remainder with the product's sign). Tiny after rounding is then
 * |x| < 2^-1022 - 2^-1076, which is (2^56 - 4) units; the subnormal
 * quantum is 16 units. */
__host__ __device__ WAVE_INLINE int uf_fma(double a, double b, double c){
    uint64_t Ma, Mb, Mc; int Ea, Eb, Ec;
    decomp(a, Ma, Ea); decomp(b, Mb, Eb); decomp(c, Mc, Ec);
    const int top_p = Ea + Eb + bitlen64(Ma) + bitlen64(Mb);   /* |a b| < 2^top_p */
    const int top_c = Ec + bitlen64(Mc);                         /* |c| < 2^top_c */
    if (top_p > -966 || top_c > -966) return 0;
    const int sp = sign_of(a) ^ sign_of(b), sc = sign_of(c);
    int sticky = 0;
    U128 P = u128_mul(Ma, Mb);
    const int sh = Ea + Eb + 1078;
    if (sh >= 0) P = u128_shl(P, sh); else P = u128_shr(P, -sh, &sticky);
    U128 C; C.hi = 0; C.lo = Mc;
    C = u128_shl(C, Ec + 1078);                                   /* Ec >= -1074: exact */
    U128 S; int ss;                                               /* |S| and its sign */
    if (sp == sc){ S = u128_add(P, C); ss = sp; }
    else if (u128_cmp(P, C) >= 0){ S = u128_sub(P, C); ss = sp; }
    else { S = u128_sub(C, P); ss = sc; }
    const int s_zero = (S.hi | S.lo) == 0;
    if (s_zero && !sticky) return 0;                              /* exactly zero */
    U128 T0; T0.hi = 0; T0.lo = (1ull << 56) - 4;
    const int cmp = u128_cmp(S, T0);
    int tiny;
    if (cmp < 0) tiny = 1;
    else if (cmp > 0) tiny = 0;
    else tiny = sticky && (sp != ss);                             /* the remainder pulls |x| below */
    const int inexact = sticky || (S.lo & 15u) != 0;
    return tiny && inexact;
}

/* ---- the operations --------------------------------------------------- */
/* A result is quiet - none of the four raised - when it is finite and its
 * magnitude exceeds 2^-1022, whatever the operation (an add's subnormal
 * result is exact, so an add can never underflow). Anything else goes to
 * the slow path, which decides exactly. */
__host__ __device__ WAVE_INLINE int no_flag(double r){
    const double ar = w_abs(r);
    return ar > DBL_MIN && ar <= DBL_MAX;
}

/* The slow path of an add or a subtract reads its operands' classes only -
 * NaN, signaling, finite - and none of them depends on a sign, so w_sub
 * hands it b as it is. It once handed it w_neg(b): a sign flip in bits,
 * which nvcc compiles as a negation that quiets a signaling NaN, so a - sNaN
 * raised no INVALID on the RTX 5060 Ti (gate 4, MEASURED: 3 cases of
 * 1,064,352; the host twin's g++ keeps the bits, and found none). */
__host__ __device__ inline void slow_add(double a, double b, double r, Fl &fl, unsigned kind){
    if (is_nan(r)){
        if ((!is_nan(a) && !is_nan(b)) || is_snan(a) || is_snan(b)) fl.raise(F_INVALID, kind);   /* inf - inf */
        return;
    }
    if (is_inf(r) && is_finite(a) && is_finite(b)) fl.raise(F_OVERFLOW, kind);
}
/* WAVE_UNCHECKED compiles every check out, for one build only: the
 * gate's contraction control (tools/check_wave.py). A result's check reads
 * the product, so with the checks in a compiler may not fuse a*b + c even
 * when asked to (MEASURED with MinGW's g++ 16.1: -ffp-contract=fast on the
 * checked source changed no bit). Without them it can, and that build
 * must fail the gate, while the unchecked build with contraction off must
 * pass it: the checks move no bit, contraction does. */
#ifdef WAVE_UNCHECKED
#  define WAVE_CHECK(cond, call) ((void)fl)
#else
#  define WAVE_CHECK(cond, call) do { if (cond) call; } while (0)
#endif

__host__ __device__ WAVE_INLINE double w_add(double a, double b, Fl &fl){
    const double r = a + b;
    WAVE_CHECK(!(w_abs(r) <= DBL_MAX), slow_add(a, b, r, fl, K_ADD));
    return r;
}
__host__ __device__ WAVE_INLINE double w_sub(double a, double b, Fl &fl){
    const double r = a - b;
    WAVE_CHECK(!(w_abs(r) <= DBL_MAX), slow_add(a, b, r, fl, K_SUB));
    return r;
}

__host__ __device__ inline void slow_mul(double a, double b, double r, Fl &fl){
    if (is_nan(a) || is_nan(b)){ if (is_snan(a) || is_snan(b)) fl.raise(F_INVALID, K_MUL); return; }
    if (is_nan(r)){ fl.raise(F_INVALID, K_MUL); return; }                 /* 0 * inf */
    if (!is_finite(a) || !is_finite(b)) return;                         /* an exact infinity */
    if (is_inf(r)){ fl.raise(F_OVERFLOW, K_MUL); return; }
    if (is_zero(a) || is_zero(b)) return;                               /* an exact zero */
    if (uf_mul(a, b)) fl.raise(F_UNDERFLOW, K_MUL);
}
__host__ __device__ WAVE_INLINE double w_mul(double a, double b, Fl &fl){
    const double r = a * b;
    WAVE_CHECK(!no_flag(r), slow_mul(a, b, r, fl));
    return r;
}

__host__ __device__ inline void slow_div(double a, double b, double r, Fl &fl){
    if (is_nan(a) || is_nan(b)){ if (is_snan(a) || is_snan(b)) fl.raise(F_INVALID, K_DIV); return; }
    if (is_zero(b)){
        if (is_zero(a)) fl.raise(F_INVALID, K_DIV);                     /* 0/0 */
        else if (is_finite(a)) fl.raise(F_DIVBYZERO, K_DIV);
        return;                                                         /* inf/0 is an exact infinity */
    }
    if (is_inf(a)){
        if (is_inf(b)) fl.raise(F_INVALID, K_DIV);                      /* inf/inf */
        return;
    }
    if (is_inf(b) || is_zero(a)) return;                                /* an exact zero */
    if (is_inf(r)){ fl.raise(F_OVERFLOW, K_DIV); return; }
    (void)r;
    if (uf_div(a, b)) fl.raise(F_UNDERFLOW, K_DIV);
}
__host__ __device__ WAVE_INLINE double w_div(double a, double b, Fl &fl){
    const double r = a / b;
    WAVE_CHECK(!no_flag(r), slow_div(a, b, r, fl));
    return r;
}

/* A square root never overflows, and the root of a positive binary64 is at
 * least 2^-537, so it never underflows: only the root of a negative. */
__host__ __device__ WAVE_INLINE double w_sqrt(double a, Fl &fl){
    const double r = hw_sqrt(a);
    WAVE_CHECK((sign_of(a) && !is_zero(a) && !is_nan(a)) || is_snan(a), fl.raise(F_INVALID, K_SQRT));
    return r;
}

__host__ __device__ inline void slow_fma(double a, double b, double c, double r, Fl &fl){
    if (is_nan(a) || is_nan(b) || is_nan(c)){ if (is_snan(a) || is_snan(b) || is_snan(c)) fl.raise(F_INVALID, K_FMA); return; }
    if ((is_inf(a) && is_zero(b)) || (is_zero(a) && is_inf(b))){ fl.raise(F_INVALID, K_FMA); return; }
    if (is_inf(a) || is_inf(b)){
        if (is_inf(c) && (sign_of(a) ^ sign_of(b)) != sign_of(c)) fl.raise(F_INVALID, K_FMA);
        return;
    }
    if (is_inf(c)) return;                                              /* an exact infinity */
    if (is_inf(r)){ fl.raise(F_OVERFLOW, K_FMA); return; }
    if (is_zero(a) || is_zero(b)) return;                               /* r is c, exactly */
    if (is_zero(c)){ if (uf_mul(a, b)) fl.raise(F_UNDERFLOW, K_FMA); return; }
    if (uf_fma(a, b, c)) fl.raise(F_UNDERFLOW, K_FMA);
}
/* The port's FMA form issues fma at exactly the sites the step calls this. */
__host__ __device__ WAVE_INLINE double w_fma(double a, double b, double c, Fl &fl){
    const double r = hw_fma(a, b, c);
    WAVE_CHECK(!no_flag(r), slow_fma(a, b, c, r, fl));
    return r;
}

/* ---- augmentedAddition (IEEE 754-2019 9.5), binary64 ------------------
 * r = roundTiesTowardZero(a + b), e = (a + b) - r exactly, e signed as r
 * when it is zero; libcft's cft_augmented_add (host/src/augmented.c, read
 * at the pin), which the port's exact-time clock calls every step.
 *
 * TwoSum gives s = RNE(a + b) and its exact error t. The two roundings
 * differ only on a tie - a + b exactly halfway between s and its
 * neighbour n on t's side - where RNE took the even one and RTTZ takes
 * the one nearer zero. So: n is s's neighbour toward t; a tie is
 * 2t == n - s (both exact); and if n is nearer zero than s, r = n and
 * e = -t. Ties cannot occur where the two neighbours are subnormal (the
 * operands' sum is on their grid), so the halving is exact.
 *
 * The overflow edge: RNE overflows on the tie DBL_MAX + ulp/2 and RTTZ
 * does not; there r = DBL_MAX and e = ulp/2 = 2^970. Past it both are
 * infinite, r = e = the infinity, with OVERFLOW. Flags as augmented.c
 * raises them: INVALID for inf + -inf, OVERFLOW above, and UNDERFLOW when
 * e is non-zero and below 2^-1022 in magnitude.
 */
__host__ __device__ inline void w_augadd(double a, double b, double &r_out, double &e_out, Fl &fl){
    if (is_nan(a) || is_nan(b)){
        if (is_snan(a) || is_snan(b)) fl.raise(F_INVALID, K_AUG);
        r_out = e_out = a + b; return;
    }
    if (is_inf(a) || is_inf(b)){
        if (is_inf(a) && is_inf(b) && sign_of(a) != sign_of(b)){
            fl.raise(F_INVALID, K_AUG); r_out = e_out = a + b; return;
        }
        r_out = e_out = is_inf(a) ? a : b;
        return;
    }
    if (is_zero(a) && is_zero(b)){ r_out = e_out = a + b; return; }     /* 6.3's zero, e the same */
    if (is_zero(a) || is_zero(b)){
        const double keep = is_zero(a) ? b : a;
        r_out = keep; e_out = w_copysign(0.0, keep);
        return;
    }
    const double s = a + b;
    if (is_inf(s)){
        /* The only finite RTTZ result past RNE's overflow is the tie at
         * DBL_MAX + 2^970. Both operands carry s's sign here (with
         * opposite signs |a + b| would be under the larger, which is
         * finite), so the larger, big, lies within a factor two of DBL_MAX
         * and big - DBL_MAX is exact (Sterbenz); TwoSum of that and the
         * smaller then says exactly whether a + b - DBL_MAX == 2^970. */
        const double dm = w_copysign(u2d(0x7fefffffffffffffull), s);
        const int a_big = w_abs(a) >= w_abs(b);
        const double big = a_big ? a : b, small = a_big ? b : a;
        const double d = big - dm;
        const double u = d + small, ub = u - d;
        const double err = (d - (u - ub)) + (small - ub);
        if (u == w_copysign(0x1p970, s) && is_zero(err)){
            r_out = dm; e_out = u;
            return;
        }
        fl.raise(F_OVERFLOW, K_AUG);
        r_out = e_out = s;
        return;
    }
    const double bb = s - a;
    const double t = (a - (s - bb)) + (b - bb);                /* TwoSum: a + b == s + t exactly */
    double r = s, e = t;
    if (!is_zero(t)){
        uint64_t us = d2u(s);
        /* n: s's neighbour toward t - one step down in magnitude when t has
         * the other sign, up when the same; s is non-zero here (s == 0
         * would make a + b == 0 exactly, and then t == 0) */
        const int toward_zero = sign_of(t) != sign_of(s);
        const double n = u2d(toward_zero ? us - 1 : us + 1);
        if (2.0 * t == n - s && toward_zero){ r = n; e = w_neg(t); }
    }
    if (is_zero(e)) e = w_copysign(0.0, r);
    if (!is_zero(e) && w_abs(e) < DBL_MIN) fl.raise(F_UNDERFLOW, K_AUG);
    r_out = r; e_out = e;
}

} /* namespace wave */

#endif /* WAVE_ARITH_H */
