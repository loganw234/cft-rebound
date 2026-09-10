/* SPDX-License-Identifier: GPL-3.0-or-later
 *
 * Exact hexadecimal text for binary64, independent of the C library.
 *
 * MSVCRT's printf has no %a and its strtod reads no hex floats, and
 * this repository's records must be the same bytes on every host. So
 * both directions are done here on the bits, and the text form is the
 * canonical one libcft's cft_to_hex_char writes: a single leading 1,
 * no trailing zeros in the fraction, an explicitly signed binary
 * exponent, "0x0p+0" for zero, "inf"/"nan" for the specials. What
 * one side writes the other side reads back bit for bit. */
#ifndef CFT_REBOUND_HEXFLOAT_H
#define CFT_REBOUND_HEXFLOAT_H

#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <ctype.h>

static inline uint64_t hexfloat_bits(double x){
    uint64_t u; memcpy(&u, &x, sizeof u); return u;
}
static inline double hexfloat_from_bits(uint64_t u){
    double x; memcpy(&x, &u, sizeof x); return x;
}

/* Writes at most 32 bytes including the NUL. */
static inline void hexfloat_print(double x, char *out){
    uint64_t u = hexfloat_bits(x);
    int sign = (int)(u >> 63);
    int bexp = (int)((u >> 52) & 0x7ff);
    uint64_t frac = u & ((1ull << 52) - 1);
    char *p = out;
    if (sign) *p++ = '-';
    if (bexp == 0x7ff){
        strcpy(p, frac ? "nan" : "inf");
        return;
    }
    if (bexp == 0 && frac == 0){
        strcpy(p, "0x0p+0");
        return;
    }
    int e;
    if (bexp == 0){
        /* subnormal: normalise so the leading 1 is explicit */
        e = -1022;
        while (!(frac & (1ull << 52))){ frac <<= 1; e--; }
        frac &= (1ull << 52) - 1;
    }else{
        e = bexp - 1023;
    }
    char digs[14];
    int n = 0;
    for (int i = 12; i >= 0; i--){
        int d = (int)((frac >> (4 * i)) & 0xf);
        digs[n++] = (char)(d < 10 ? '0' + d : 'a' + d - 10);
    }
    while (n > 0 && digs[n - 1] == '0') n--;
    digs[n] = 0;
    if (n) sprintf(p, "0x1.%sp%+d", digs, e);
    else   sprintf(p, "0x1p%+d", e);
}

/* Reads the form above (also "0x1.8000000000000p-1" with trailing
 * zeros and "0x0.0p+0", which Python's float.hex() writes). Returns 1
 * on success. Only values exactly representable in binary64 are
 * accepted - more than 13 significant hex digits after the point is
 * refused rather than rounded, because a record file is never
 * supposed to hold a value that needs rounding. */
static inline int hexfloat_parse(const char *s, double *out){
    int sign = 0;
    if (*s == '-'){ sign = 1; s++; } else if (*s == '+') s++;
    if (strcmp(s, "inf") == 0){ *out = hexfloat_from_bits((uint64_t)sign << 63 | 0x7ffull << 52); return 1; }
    if (strcmp(s, "nan") == 0){ *out = hexfloat_from_bits((uint64_t)sign << 63 | 0x7ffull << 52 | 1ull << 51); return 1; }
    if (!(s[0] == '0' && (s[1] == 'x' || s[1] == 'X'))) return 0;
    s += 2;
    uint64_t mant = 0;
    int nd = 0, seen_point = 0, fracdigits = 0;
    for (; *s && *s != 'p' && *s != 'P'; s++){
        if (*s == '.'){ if (seen_point) return 0; seen_point = 1; continue; }
        int c = tolower((unsigned char)*s);
        int d;
        if (c >= '0' && c <= '9') d = c - '0';
        else if (c >= 'a' && c <= 'f') d = c - 'a' + 10;
        else return 0;
        if (mant == 0 && d == 0){ if (seen_point) fracdigits++; continue; }
        if (nd >= 16) return 0;
        mant = mant << 4 | (uint64_t)d;
        nd++;
        if (seen_point) fracdigits++;
    }
    if (*s != 'p' && *s != 'P') return 0;
    s++;
    int esign = 1;
    if (*s == '-'){ esign = -1; s++; } else if (*s == '+') s++;
    if (!isdigit((unsigned char)*s)) return 0;
    long e = 0;
    for (; isdigit((unsigned char)*s); s++) e = e * 10 + (*s - '0');
    if (*s) return 0;
    e = esign * e - 4L * fracdigits;
    if (mant == 0){ *out = sign ? -0.0 : 0.0; return 1; }
    /* value = mant * 2^e, exact in binary64 if mant fits 53 bits after
     * stripping trailing zero bits, and the exponent is in range */
    while (!(mant & 1)){ mant >>= 1; e++; }
    int bits = 0; for (uint64_t t = mant; t; t >>= 1) bits++;
    if (bits > 53) return 0;
    /* shift so that mant has exactly 53 bits: value = m53 * 2^(e - (53-bits)) */
    uint64_t m53 = mant << (53 - bits);
    long e2 = e - (53 - bits);
    long bexp = e2 + 52 + 1023;   /* exponent field for 1.f * 2^(e2+52) */
    if (bexp >= 0x7ff) return 0;
    if (bexp <= 0){
        /* subnormal: only exact if the dropped bits are zero */
        long shift = 1 - bexp;
        if (shift > 53) return 0;
        if (m53 & ((1ull << shift) - 1)) return 0;
        uint64_t u = (uint64_t)sign << 63 | (m53 >> shift);
        *out = hexfloat_from_bits(u);
        return 1;
    }
    uint64_t u = (uint64_t)sign << 63 | (uint64_t)bexp << 52 | (m53 & ((1ull << 52) - 1));
    *out = hexfloat_from_bits(u);
    return 1;
}

#endif
