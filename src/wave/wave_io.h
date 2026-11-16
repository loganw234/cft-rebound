/* SPDX-License-Identifier: GPL-3.0-or-later
 * Copyright 2026 the cft-rebound contributors.
 *
 * The host's text: binary64 to and from libcft's hex spelling, and a
 * problem file read as src/ias15_cft.c's read_problem() reads it. Host
 * only; the kernel never sees text.
 *
 *   wave_to_hex()    cft_to_hex_char at binary64 (cft-fp256 host/src/
 *                    chars.c, to_hex_one(), read at the pin): a leading 1
 *                    whatever the binade, so the smallest subnormal is
 *                    0x1p-1074; trailing zero nibbles dropped; zero
 *                    0x0p+0; inf, nan, and a NaN's payload in brackets.
 *   wave_from_hex()  cft_from_hex_char at binary64, round to nearest even
 *                    (from_hex_one()): [+-]0x<hex>[.<hex>]p[+-]<dec>, the
 *                    exponent required, or inf, infinity, nan.
 * wave_corpus.cpp holds both to libcft's own, on the host, value by value.
 */
#ifndef WAVE_IO_H
#define WAVE_IO_H

#include "wave_arith.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

namespace wave {

/* buf needs 32 bytes */
inline void wave_to_hex(double x, char *buf){
    const uint64_t u = d2u(x);
    char *o = buf;
    if (u >> 63) *o++ = '-';
    const uint64_t be = (u >> 52) & 0x7ff, frac = u & 0x000fffffffffffffull;
    if (be == 0x7ff){
        if (!frac){ strcpy(o, "inf"); return; }
        const int signaling = !((frac >> 51) & 1);
        const uint64_t pay = frac & 0x0007ffffffffffffull;
        o += sprintf(o, "%s", signaling ? "snan" : "nan");
        if (pay) sprintf(o, "(0x%llx)", (unsigned long long)pay);
        return;
    }
    if (be == 0 && frac == 0){ strcpy(o, "0x0p+0"); return; }
    uint64_t m; int e;
    decomp(x, m, e);
    const int len = bitlen64(m);
    const int exp = e + len - 1;
    uint64_t f = m & ~(1ull << (len - 1));
    int nib = (len + 2) / 4;
    if (nib) f <<= (4 * nib - (len - 1));
    while (nib > 0 && (f & 0xf) == 0){ f >>= 4; nib--; }
    o += sprintf(o, "0x1");
    if (nib){
        *o++ = '.';
        for (int i = nib; i-- > 0; ) *o++ = "0123456789abcdef"[(f >> (4 * i)) & 0xf];
    }
    sprintf(o, "p%c%d", exp < 0 ? '-' : '+', exp < 0 ? -exp : exp);
}

inline int hexval(int c){
    if (c >= '0' && c <= '9') return c - '0';
    c = tolower(c);
    if (c >= 'a' && c <= 'f') return c - 'a' + 10;
    return -1;
}

inline int ci_eq(const char *s, const char *lit){
    while (*s && *lit){ if (tolower((unsigned char)*s) != *lit) return 0; s++; lit++; }
    return !*s && !*lit;
}

/* 1 and *out on success, 0 on a sequence libcft would refuse */
inline int wave_from_hex(const char *s, double *out){
    int neg = 0;
    if (*s == '+' || *s == '-'){ neg = *s == '-'; s++; }
    if (!*s) return 0;
    if (ci_eq(s, "inf") || ci_eq(s, "infinity")){ *out = neg ? -HUGE_VAL : HUGE_VAL; return 1; }
    if (ci_eq(s, "nan")){ *out = u2d((neg ? 0xfff8000000000000ull : 0x7ff8000000000000ull)); return 1; }
    if (s[0] != '0' || (s[1] != 'x' && s[1] != 'X')) return 0;
    s += 2;
    uint64_t acc = 0;         /* the leading significant digits, at most 60 bits */
    int sticky = 0, digits = 0, frac_digits = 0, started = 0, in_frac = 0;
    long dropped = 0;         /* significant digits past acc's capacity */
    for (;; s++){
        if (*s == '.'){
            if (in_frac) return 0;
            in_frac = 1;
            continue;
        }
        const int h = hexval((unsigned char)*s);
        if (h < 0) break;
        digits++;
        if (in_frac) frac_digits++;
        if (!started && h == 0) continue;
        started = 1;
        if (acc >> 56){ dropped++; if (h) sticky = 1; }
        else acc = acc * 16 + (uint64_t)h;
    }
    if (!digits) return 0;
    if (*s != 'p' && *s != 'P') return 0;
    s++;
    int eneg = 0;
    if (*s == '+' || *s == '-'){ eneg = *s == '-'; s++; }
    if (!isdigit((unsigned char)*s)) return 0;
    long long ev = 0;
    for (; isdigit((unsigned char)*s); s++){ if (ev < 100000000) ev = ev * 10 + (*s - '0'); }
    if (*s) return 0;
    if (eneg) ev = -ev;
    if (!acc && !sticky){ *out = neg ? -0.0 : 0.0; return 1; }
    /* value = (acc + sticky) * 2^q0 */
    long long q0 = ev - 4LL * frac_digits + 4LL * dropped;
    const int L = bitlen64(acc);
    const long long top = q0 + L - 1;                 /* the leading bit's exponent */
    if (top > 1023){ *out = neg ? -HUGE_VAL : HUGE_VAL; return 1; }
    long long q = top - 52;                           /* the result's quantum */
    if (q < -1074) q = -1074;
    uint64_t M;
    if (q <= q0){
        M = acc << (q0 - q);                          /* exact: no digit was dropped (acc < 2^57) */
    }else{
        /* (acc + d) * 2^q0, 0 <= d < 1 and d > 0 iff sticky, to a multiple
         * of 2^q, ties to even. acc < 2^60, so past a shift of 63 it is
         * under half a quantum and rounds to zero. */
        const long long sh = q - q0;
        if (sh >= 64) M = 0;
        else{
            M = acc >> sh;
            const uint64_t rem = acc & ((1ull << sh) - 1);
            const uint64_t half = 1ull << (sh - 1);
            if (rem > half || (rem == half && (sticky || (M & 1)))) M++;
        }
        if (M >> 53){ M >>= 1; q++; }                 /* the rounding carried out */
        if (q + 52 > 1023 && (M >> 52)){ *out = neg ? -HUGE_VAL : HUGE_VAL; return 1; }
    }
    double v = ldexp((double)M, (int)q);
    *out = neg ? -v : v;
    return 1;
}

} /* namespace wave */

#endif /* WAVE_IO_H */
