#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Derive the IAS15 constants, never transcribe them.
#
# REBOUND's src/integrator_ias15.c hard-codes the Gauss-Radau spacings
# h[8] and the derived arrays rr[28], c[21], d[21] as double literals
# with 25 to 31 significant decimal digits. binary128 needs 36 and
# binary256 needs 73 (754-2019 5.12.2's Pmin), so a port that copied
# those literals would round every constant to binary64 and cap the
# whole integrator at binary64 accuracy while claiming binary256.
#
# This script derives them from their definitions at 130 decimal digits
# with mpmath and emits them per target format, correctly rounded, in
# hexadecimal (exact) and decimal (100 significant digits) character
# form, then checks the derivation two ways:
#
#   1. every non-zero node h[n] satisfies its defining polynomial to
#      the working precision, under two independent definitions of the
#      Radau nodes (the (P7+P8)/(1+x) quotient and a Jacobi polynomial);
#   2. the derived values, correctly rounded to binary64, reproduce
#      REBOUND's published double literals bit for bit - all 78 of
#      them - and the binary64 rounding is done twice by different
#      code (this file's integer rounding and CPython's Fraction ->
#      float), which must agree.
#
# A mismatch in either check means the wrong thing was derived, and
# the script exits non-zero without writing anything.
#
# Usage:  python tools/gen_constants.py [--header PATH] [--json PATH]
#         [--rebound-src PATH/integrator_ias15.c]
#
# The rebound source is only READ, to extract the published literals
# for check 2. Nothing in the emitted constants comes from it.

import argparse
import json
import os
import re
import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

DPS = 130
mp.dps = DPS

FORMATS = [
    # name,   precision p (significand bits incl. hidden), emin
    ("fp64",  53,  -1022),
    ("fp128", 113, -16382),
    ("fp256", 237, -262142),
]
DEC_DIGITS = 100  # > Pmin(fp256) = 73, so the decimal form round-trips


# --------------------------------------------------------------------
# Exact Legendre polynomials as Fraction coefficient lists (low->high)
# --------------------------------------------------------------------
def poly_add(a, b):
    n = max(len(a), len(b))
    return [(a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0) for i in range(n)]

def poly_scale(a, s):
    return [s * x for x in a]

def poly_mul_x(a):
    return [Fraction(0)] + list(a)

def legendre_coeffs(n):
    """Exact coefficients of P_n via the Bonnet recurrence."""
    P0 = [Fraction(1)]
    if n == 0:
        return P0
    P1 = [Fraction(0), Fraction(1)]
    for k in range(1, n):
        # (k+1) P_{k+1} = (2k+1) x P_k - k P_{k-1}
        t = poly_add(poly_scale(poly_mul_x(P1), Fraction(2 * k + 1)),
                     poly_scale(P0, Fraction(-k)))
        P0, P1 = P1, poly_scale(t, Fraction(1, k + 1))
    return P1

def poly_div_linear(a, r):
    """Divide polynomial a (low->high) by (x - r) exactly; return (q, rem)."""
    n = len(a) - 1
    q = [Fraction(0)] * n
    rem = a[n]
    for i in range(n - 1, -1, -1):
        q[i] = rem
        rem = a[i] + rem * r
    return q, rem

def poly_eval_mp(a, x):
    acc = mpf(0)
    for coef in reversed(a):
        acc = acc * x + mpf(coef.numerator) / mpf(coef.denominator)
    return acc

def poly_deriv(a):
    return [a[i] * i for i in range(1, len(a))]


# --------------------------------------------------------------------
# The nodes
# --------------------------------------------------------------------
def radau_nodes():
    """The seven free nodes of the 8-point Gauss-Radau rule on [-1,1]
    with the fixed node at x = -1, mapped to [0,1] by h = (1+x)/2.

    Definition A: the zeros of (P_7(x) + P_8(x)) / (1 + x).
    """
    p = poly_add(legendre_coeffs(7), legendre_coeffs(8))
    q, rem = poly_div_linear(p, Fraction(-1))
    if rem != 0:
        raise SystemExit("P7+P8 is not divisible by (x+1): derivation wrong")
    # polyroots wants high->low coefficients
    coeffs = [mpf(c.numerator) / mpf(c.denominator) for c in reversed(q)]
    roots = mpmath.polyroots(coeffs, maxsteps=200, extraprec=400)
    roots = sorted(mpmath.re(r) for r in roots)
    if len(roots) != 7:
        raise SystemExit("expected 7 roots, got %d" % len(roots))
    # Polish each root by Newton on the exact polynomial.
    dq = poly_deriv(q)
    polished = []
    for r in roots:
        x = r
        for _ in range(8):
            x = x - poly_eval_mp(q, x) / poly_eval_mp(dq, x)
        polished.append(x)
    return q, polished


def check_nodes(q, xs):
    """Check 1: each node satisfies its defining polynomial(s)."""
    ok = True
    tol = mpf(10) ** (-(DPS - 10))
    # Definition A: the quotient polynomial itself.
    for i, x in enumerate(xs):
        v = abs(poly_eval_mp(q, x))
        if v > tol:
            print("FAIL node %d: |(P7+P8)/(1+x)| = %s" % (i + 1, mpmath.nstr(v, 5)))
            ok = False
    # Definition B: an independent one. The free Radau nodes with the
    # endpoint at -1 are also the zeros of the Jacobi polynomial
    # P_7^(0,1)(x) (Hildebrand, Introduction to Numerical Analysis,
    # 8.9). We test BOTH orderings of the Jacobi parameters and report
    # which one vanishes, so this is a measurement rather than an
    # assertion.
    which = None
    for (alpha, beta) in ((0, 1), (1, 0)):
        worst = mpf(0)
        for x in xs:
            worst = max(worst, abs(mpmath.jacobi(7, alpha, beta, x)))
        print("  Jacobi P_7^(%d,%d) at the nodes: max |value| = %s" %
              (alpha, beta, mpmath.nstr(worst, 5)))
        if worst < tol:
            which = (alpha, beta)
    if which is None:
        print("FAIL: no Jacobi polynomial vanishes at the nodes")
        ok = False
    else:
        print("  definition B satisfied: nodes are the zeros of P_7^(%d,%d)" % which)
    # And the nodes must be distinct, in (-1, 1), and there must be 7.
    for i in range(len(xs) - 1):
        if not (xs[i] < xs[i + 1]):
            ok = False
            print("FAIL: nodes not increasing")
    if not (xs[0] > -1 and xs[-1] < 1):
        ok = False
        print("FAIL: node outside (-1,1)")
    return ok


# --------------------------------------------------------------------
# The derived arrays, by REBOUND's own recurrence (GENERATE_CONSTANTS)
# --------------------------------------------------------------------
def derive_arrays(h):
    rr = []
    for j in range(1, 8):
        for k in range(j):
            rr.append(h[j] - h[k])
    assert len(rr) == 28
    c = [mpf(0)] * 21
    d = [mpf(0)] * 21
    c[0] = -h[1]
    d[0] = h[1]
    l = 0
    for j in range(2, 7):
        l += 1
        c[l] = -h[j] * c[l - j + 1]
        d[l] = h[1] * d[l - j + 1]
        for k in range(2, j):
            l += 1
            c[l] = c[l - j] - h[j] * c[l - j + 1]
            d[l] = d[l - j] + h[k] * d[l - j + 1]
        l += 1
        c[l] = c[l - j] - h[j]
        d[l] = d[l - j] + h[j]
    assert l == 20
    return rr, c, d


# --------------------------------------------------------------------
# Correct rounding to a binary format, done on exact integers
# --------------------------------------------------------------------
def mpf_to_fraction(x):
    # mpmath may be running on gmpy2, in which case man is an mpz and
    # Fraction arithmetic refuses it - normalise to Python ints.
    sign, man, exp, bc = x._mpf_
    man = int(man)
    exp = int(exp)
    if man == 0:
        return Fraction(0)
    f = Fraction(man) * (Fraction(2) ** exp)
    return -f if sign else f

def round_to_binary(F, p, emin):
    """Round the exact rational F to the nearest binary float with a
    p-bit significand (ties to even). Returns (sign, m, e) with value
    (-1)^sign * m * 2^e and 2^(p-1) <= m < 2^p, or m == 0. Values are
    assumed normal (all IAS15 constants are); the assertion checks."""
    if F == 0:
        return 0, 0, 0
    sign = 1 if F < 0 else 0
    F = abs(F)
    # find e with 2^(p-1) <= F / 2^e < 2^p
    e = F.numerator.bit_length() - F.denominator.bit_length() - p
    # adjust so the integer part is in range
    while F / (Fraction(2) ** e) >= 2 ** p:
        e += 1
    while F / (Fraction(2) ** e) < 2 ** (p - 1):
        e -= 1
    scaled = F / (Fraction(2) ** e)
    m = scaled.numerator // scaled.denominator
    rem = scaled - m
    half = Fraction(1, 2)
    if rem > half or (rem == half and (m & 1)):
        m += 1
        if m == 2 ** p:
            m = 2 ** (p - 1)
            e += 1
    assert e + p - 1 >= emin, "subnormal constant, not expected"
    return sign, m, e

def hex_float(sign, m, e, p):
    """Canonical hex-float text: one leading 1, no trailing zeros, an
    explicitly signed binary exponent - the form cft_to_hex_char
    writes and cft_from_hex_char reads exactly."""
    if m == 0:
        return "0x0p+0"
    # value = m * 2^e with m in [2^(p-1), 2^p): 1.frac * 2^(e+p-1)
    frac = m - (1 << (p - 1))
    nbits = p - 1
    ndig = (nbits + 3) // 4
    frac <<= (ndig * 4 - nbits)
    s = ("%0" + str(ndig) + "x") % frac
    s = s.rstrip("0")
    bexp = e + p - 1
    body = "0x1" + ("." + s if s else "") + ("p%+d" % bexp)
    return ("-" if sign else "") + body

def fraction_of_hex(sign, m, e):
    f = Fraction(m) * (Fraction(2) ** e)
    return -f if sign else f

def decimal_text(F, digits):
    """digits significant decimal digits of the rational F, correctly
    rounded (ties to even), in the "d.ddd...e+X" form."""
    if F == 0:
        return "0"
    sign = "-" if F < 0 else ""
    F = abs(F)
    # find k with 10^(digits-1) <= F * 10^k < 10^digits
    k = 0
    ten = Fraction(10)
    lo = ten ** (digits - 1)
    hi = ten ** digits
    while F * ten ** k >= hi:
        k -= 1
    while F * ten ** k < lo:
        k += 1
    scaled = F * ten ** k
    n = scaled.numerator // scaled.denominator
    rem = scaled - n
    if rem > Fraction(1, 2) or (rem == Fraction(1, 2) and (n & 1)):
        n += 1
        if n == hi:
            n //= 10
            k -= 1
    s = str(n)
    assert len(s) == digits
    exp10 = -k + digits - 1
    return "%s%s.%se%+d" % (sign, s[0], s[1:], exp10)


# --------------------------------------------------------------------
# REBOUND's published literals, for check 2
# --------------------------------------------------------------------
def read_rebound_literals(path):
    src = open(path, encoding="utf-8").read()
    def arr(name, n):
        m = re.search(r"static const double %s\[%d\]\s*=\s*\{([^}]*)\}" % (name, n), src)
        if not m:
            raise SystemExit("could not find %s[%d] in %s" % (name, n, path))
        vals = [v.strip() for v in m.group(1).split(",")]
        if len(vals) != n:
            raise SystemExit("%s has %d entries, expected %d" % (name, len(vals), n))
        return vals
    return {"h": arr("h", 8), "rr": arr("rr", 28), "c": arr("c", 21), "d": arr("d", 21)}


def main():
    ap = argparse.ArgumentParser()
    here = os.path.dirname(os.path.abspath(__file__))
    root = os.path.dirname(here)
    ap.add_argument("--header", default=os.path.join(root, "src", "ias15_constants.h"))
    ap.add_argument("--json", default=os.path.join(root, "data", "ias15_constants.json"))
    ap.add_argument("--rebound-src", default=os.path.join(
        root, "third_party", "rebound", "src", "integrator_ias15.c"))
    ap.add_argument("--no-write", action="store_true")
    args = ap.parse_args()

    print("cft-rebound gen_constants: mpmath %s at %d digits" % (mpmath.__version__, DPS))
    q, xs = radau_nodes()
    print("check 1: the nodes against their defining polynomials")
    ok1 = check_nodes(q, xs)
    h = [mpf(0)] + [(x + 1) / 2 for x in xs]
    rr, c, d = derive_arrays(h)

    arrays = {"h": h, "rr": rr, "c": c, "d": d}
    lits = read_rebound_literals(args.rebound_src)

    print("check 2: correctly rounded to binary64 against REBOUND's %d literals" %
          sum(len(v) for v in lits.values()))
    ok2 = True
    n_checked = 0
    for name in ("h", "rr", "c", "d"):
        for i, v in enumerate(arrays[name]):
            F = mpf_to_fraction(v)
            s, m, e = round_to_binary(F, 53, -1022)
            mine = float(fraction_of_hex(s, m, e))     # exact: fits a double
            cpython = float(F)                         # CPython rounds Fraction->float correctly
            published = float(lits[name][i])           # C literal, parsed by CPython
            n_checked += 1
            if mine != cpython:
                ok2 = False
                print("FAIL %s[%d]: my rounding %a differs from CPython's %a" % (name, i, mine, cpython))
            if mine != published:
                ok2 = False
                print("FAIL %s[%d]: derived %a, REBOUND literal %a (%s)" %
                      (name, i, mine, published, lits[name][i]))
    print("  %d constants compared, %s" % (n_checked, "all match" if ok2 else "MISMATCHES"))

    # A third, weaker sanity line: the published 31-digit h literals
    # against the derived nodes at their own digit count.
    worst = mpf(0)
    for i in range(1, 8):
        worst = max(worst, abs(h[i] - mpf(lits["h"][i])))
    print("  max |h_derived - h_literal(31 digits)| = %s" % mpmath.nstr(worst, 3))

    if not (ok1 and ok2):
        print("DERIVATION CHECK FAILED - nothing written")
        return 1

    # ---------------- emission ----------------
    out = {"dps": DPS, "dec_digits": DEC_DIGITS, "formats": {}, "decimal": {}}
    for name, vals in arrays.items():
        out["decimal"][name] = [decimal_text(mpf_to_fraction(v), DEC_DIGITS) for v in vals]
    for fname, p, emin in FORMATS:
        fmt = {}
        for name, vals in arrays.items():
            fmt[name] = [hex_float(*round_to_binary(mpf_to_fraction(v), p, emin), p=p) for v in vals]
        out["formats"][fname] = fmt

    if args.no_write:
        print("(--no-write: not emitting)")
        return 0

    os.makedirs(os.path.dirname(args.json), exist_ok=True)
    with open(args.json, "w", newline="\n") as f:
        json.dump(out, f, indent=1)
        f.write("\n")

    os.makedirs(os.path.dirname(args.header), exist_ok=True)
    with open(args.header, "w", newline="\n") as f:
        w = f.write
        w("/* Generated by tools/gen_constants.py - DO NOT EDIT.\n")
        w(" * SPDX-License-Identifier: GPL-3.0-or-later\n")
        w(" *\n")
        w(" * IAS15's Gauss-Radau spacings h[8] and the derived arrays rr[28],\n")
        w(" * c[21], d[21], derived at %d decimal digits and correctly rounded\n" % DPS)
        w(" * to each format (ties to even). The hexadecimal forms are exact\n")
        w(" * and are what the integrator loads through cft_from_hex_char; the\n")
        w(" * decimal forms carry %d significant digits and exist so the\n" % DEC_DIGITS)
        w(" * program can check, at start-up, that cft_from_decimal_char rounds\n")
        w(" * them to the same bits - two conversions, one answer.\n")
        w(" *\n")
        w(" * Format index: 0 = binary64, 1 = binary128, 2 = binary256\n")
        w(" * (CFT_FP64 - 1, CFT_FP128 - 1, CFT_FP256 - 1). */\n")
        w("#ifndef IAS15_CONSTANTS_H\n#define IAS15_CONSTANTS_H\n\n")
        w("#define IAS15_NFMT 3\n")
        for name, vals in arrays.items():
            n = len(vals)
            w("#define IAS15_N_%s %d\n" % (name.upper(), n))
        w("\n")
        for name, vals in arrays.items():
            n = len(vals)
            w("static const char *const ias15_%s_hex[IAS15_NFMT][%d] = {\n" % (name, n))
            for fname, p, emin in FORMATS:
                w("    { /* %s */\n" % fname)
                for hx in out["formats"][fname][name]:
                    w("        \"%s\",\n" % hx)
                w("    },\n")
            w("};\n")
            w("static const char *const ias15_%s_dec[%d] = {\n" % (name, n))
            for dec in out["decimal"][name]:
                w("    \"%s\",\n" % dec)
            w("};\n\n")
        w("#endif /* IAS15_CONSTANTS_H */\n")
    print("wrote %s and %s" % (args.header, args.json))
    return 0


if __name__ == "__main__":
    sys.exit(main())
