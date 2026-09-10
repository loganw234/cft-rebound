#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Divergence within an ensemble record: at every sample, the distance
# of each member's state from member 0's, from the exact bits - the
# measurement an ensemble of perturbed initial conditions exists for
# (sensitivity to initial conditions, the finite-time Lyapunov rate).
# Members must share the sampling times, which a --dt-file run
# guarantees; with per-system adaptive steps the samples of different
# members are at different times and the tool refuses.
#
#   python tools/divergence.py ENSEMBLE.rec [--csv OUT] [--reference K]
#
# Prints per sample: t, then for each member k != reference the
# largest relative coordinate difference |x_k - x_ref| / L (L the
# largest coordinate of the reference at that sample) and the relative
# speed difference; and a summary of the exponential growth rate
# between consecutive samples for each member.

import argparse
import math
import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from oracle import hex_to_fraction, F2m  # noqa: E402

mp.dps = 40


def parse(path):
    E = None
    samples = {}   # (system, k) -> (step, t, state)
    for line in open(path):
        if line.startswith("# program="):
            for kv in line[2:].split():
                k, _, v = kv.partition("=")
                if k == "E":
                    E = int(v)
        elif line.startswith("system ") and " sample " in line:
            toks = line.split()
            s = int(toks[1]); k = int(toks[3]); step = int(toks[4])
            body = toks[5:]
            t = None
            if "|" in body:
                i = body.index("|")
                t = hex_to_fraction(body[i + 1]) + hex_to_fraction(body[i + 2])
                body = body[:i]
            vals = [hex_to_fraction(x) for x in body]
            samples[(s, k)] = (step, t if t is not None else vals[0], vals[4:])
    if E is None:
        raise SystemExit("not an ensemble record")
    return E, samples


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("record")
    ap.add_argument("--csv")
    ap.add_argument("--reference", type=int, default=0)
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    E, S = parse(args.record)
    ref = args.reference
    ks = sorted(k for (s, k) in S if s == ref)
    out = open(args.csv, "w", newline="\n") if args.csv else None
    if out:
        out.write("sample,step,t," + ",".join("dx_m%d,dv_m%d" % (k, k) for k in range(E) if k != ref) + "\n")
    prev = None
    rates = {k: [] for k in range(E)}
    for k in ks:
        step, t, st = S[(ref, k)]
        N = len(st) // 6
        L = max(abs(q) for i in range(N) for q in st[6 * i:6 * i + 3]) or Fraction(1)
        V = max(abs(q) for i in range(N) for q in st[6 * i + 3:6 * i + 6]) or Fraction(1)
        row = []
        for m in range(E):
            if m == ref:
                continue
            stm = S[(m, k)][2]
            if S[(m, k)][1] != t:
                raise SystemExit("member %d is at t=%s when member %d is at t=%s (sample %d): the members do not share their times" %
                                 (m, float(S[(m, k)][1]), ref, float(t), k))
            dx = max(abs(stm[6 * i + c] - st[6 * i + c]) for i in range(N) for c in range(3)) / L
            dv = max(abs(stm[6 * i + 3 + c] - st[6 * i + 3 + c]) for i in range(N) for c in range(3)) / V
            row.append((float(dx), float(dv)))
            if prev is not None and prev[1][m - (1 if m > ref else 0)][0] > 0 and dx > 0 and t > prev[0]:
                rates[m].append(math.log(float(dx) / prev[1][m - (1 if m > ref else 0)][0]) / float(t - prev[0]))
        if out:
            out.write("%d,%d,%s," % (k, step, mpmath.nstr(F2m(t), 18)) + ",".join("%.6e,%.6e" % r for r in row) + "\n")
        if not args.quiet:
            print("sample %4d t %14.8g  " % (k, float(t)) + "  ".join("m%d %.3e" % (m, row[j][0]) for j, m in enumerate(mm for mm in range(E) if mm != ref)))
        prev = (t, row)
    if out:
        out.close()
    print("DIVERGENCE %s: %d members, %d samples, reference member %d; largest final |dx|/L %s" %
          (args.record, E, len(ks), ref, "%.3e" % max(r[0] for r in row) if row else "-"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
