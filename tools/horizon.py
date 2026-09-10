#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The horizon: from the oracle CSVs of an ensemble of runs (one CSV per
# member, the same sampling grid), the distribution of the error at
# every sample - median, 90th percentile, maximum - and the time at
# which a chosen percentile of the phase error crosses each threshold,
# measured where the data reach and extrapolated by a power-law fit of
# the percentile curve where they do not (marked as such).
#
# A single run of a round-off random walk does not fit a power law
# (docs/VALIDATION.md, the horizon entry); the percentile of an
# ensemble does. The phase error is the along-track timing error in
# units of the period (tools/oracle.py), i.e. "how late is the planet",
# which is the error a user of a long integration cares about.
#
#   python tools/horizon.py CSV [CSV ...] [--column rel_phase_err] [--percentile 90]
#          [--thresholds 1e-9,1e-6,1e-3,1e-1] [--fit-from ORBITS] [--table OUT.md]

import argparse
import csv
import glob
import math
import os
import sys


def load(paths):
    members = []
    for p in paths:
        rows = list(csv.DictReader(open(p)))
        members.append([(float(r["orbits"]), {k: float(v) for k, v in r.items() if k not in ("sample", "step")}) for r in rows])
    n = min(len(m) for m in members)
    if any(len(m) != n for m in members):
        print("note: members have different lengths; using the first %d samples of each" % n)
    return [m[:n] for m in members]


def percentile(xs, p):
    xs = sorted(xs)
    if len(xs) == 1:
        return xs[0]
    k = (len(xs) - 1) * p / 100.0
    lo = int(math.floor(k)); hi = min(lo + 1, len(xs) - 1)
    return xs[lo] + (xs[hi] - xs[lo]) * (k - lo)


def fit(pairs):
    xs = [math.log10(t) for t, y in pairs if t > 0 and y > 0]
    ys = [math.log10(y) for t, y in pairs if t > 0 and y > 0]
    n = len(xs)
    if n < 3:
        return None
    mx, my = sum(xs) / n, sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return None
    a = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    b = my - a * mx
    r = math.sqrt(sum((y - (a * x + b)) ** 2 for x, y in zip(xs, ys)) / n)
    return a, b, n, r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csvs", nargs="+")
    ap.add_argument("--column", default="rel_phase_err")
    ap.add_argument("--percentile", type=float, default=90.0)
    ap.add_argument("--thresholds", default="1e-9,1e-6,1e-3,1e-1")
    ap.add_argument("--fit-from", type=float, default=0.0, help="fit the percentile curve from this many orbits on")
    ap.add_argument("--rows", type=int, default=12, help="how many rows of the distribution to print")
    ap.add_argument("--table", default=None)
    args = ap.parse_args()
    paths = []
    for c in args.csvs:
        paths.extend(sorted(glob.glob(c)) or [c])
    members = load(paths)
    E = len(members)
    n = len(members[0])
    col = args.column
    thresholds = [float(t) for t in args.thresholds.split(",")]
    print("horizon: %d members, %d samples, column %s, percentile %g" % (E, n, col, args.percentile))
    dist = []
    for i in range(n):
        t = members[0][i][0]
        vals = [m[i][1][col] for m in members]
        dist.append((t, percentile(vals, 50), percentile(vals, args.percentile), max(vals), min(vals)))
    step = max(1, (n - 1) // args.rows)
    lines = []
    lines.append("| orbits | median | p%g | max | min |" % args.percentile)
    lines.append("|---|---|---|---|---|")
    for t, med, pp, mx, mn in dist[1::step] + ([dist[-1]] if (n - 1) % step else []):
        lines.append("| %.4g | %.3e | %.3e | %.3e | %.3e |" % (t, med, pp, mx, mn))
    print("\n".join(lines))
    # measured crossings of the percentile curve
    print("crossings of the p%g curve (first sample at or above the threshold):" % args.percentile)
    measured = {}
    for thr in thresholds:
        hit = next((t for t, med, pp, mx, mn in dist[1:] if pp >= thr), None)
        measured[thr] = hit
        print("   %g: %s" % (thr, ("%.4g orbits" % hit) if hit else "not reached in %.4g orbits" % dist[-1][0]))
    f = fit([(t, pp) for t, med, pp, mx, mn in dist[1:] if t >= args.fit_from])
    if f:
        a, b, nn, r = f
        print("fit of the p%g curve over %d samples (orbits >= %g): %s = %.3e * orbits^%.3f, rms residual %.2f dex" % (args.percentile, nn, args.fit_from, col, 10 ** b, a, r))
        print("extrapolated crossings (power law):")
        for thr in thresholds:
            if a > 0:
                print("   %g: %.3g orbits%s" % (thr, 10 ** ((math.log10(thr) - b) / a), "" if measured[thr] else " (extrapolated)"))
    fm = fit([(t, med) for t, med, pp, mx, mn in dist[1:] if t >= args.fit_from])
    if fm:
        print("fit of the median: %s = %.3e * orbits^%.3f, rms residual %.2f dex" % (col, 10 ** fm[1], fm[0], fm[3]))
    if args.table:
        with open(args.table, "w", newline="\n") as out:
            out.write("\n".join(lines) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
