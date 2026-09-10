#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The round-off floor of a format, measured directly: two records of
# the SAME fixed-step run at two formats share the same truncation
# error (same scheme, same dt, same constants to within their rounding),
# so their difference at a common step is the narrower format's
# round-off - not the total error the oracle reports, which the
# truncation dominates whenever the step is one IAS15 would normally
# take. The wider record is the reference; at binary256 its own
# round-off is 2^-184 of the narrower's and drops out.
#
#   python tools/compare_formats.py NARROW.rec WIDE.rec [--csv OUT]
#
# Prints per common sample: the relative energy difference and the
# largest relative coordinate difference, and a summary.

import argparse
import os
import sys

import mpmath
from mpmath import mp, mpf

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from oracle import parse, hex_to_fraction, F2m, invariants  # noqa: E402

mp.dps = 60


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("narrow")
    ap.add_argument("wide")
    ap.add_argument("--csv")
    ap.add_argument("--quiet", action="store_true")
    ap.add_argument("--allow-dt-mismatch", action="store_true",
                    help="compare even when the two formats rounded the decimal dt differently (the result then includes the time offset)")
    args = ap.parse_args()
    hn, mn, sn, _ = parse(args.narrow)
    hw, mw, sw, _ = parse(args.wide)
    if hn.get("problem") != hw.get("problem"):
        print("records are not the same problem: %s/%s" % (hn.get("problem"), hw.get("problem")))
        return 2
    if hn.get("dt0") != hw.get("dt0"):
        print("WARNING: the two formats rounded dt differently (%s vs %s): a non-dyadic decimal step; the difference below includes the resulting time offset" % (hn.get("dt0"), hw.get("dt0")))
        if not args.allow_dt_mismatch:
            print("refusing; pass --allow-dt-mismatch to compare anyway, or use a dyadic dt")
            return 2
    G = F2m(hex_to_fraction(hw["G"]))
    wide_by_step = {s["step"]: s for s in sw}
    out = open(args.csv, "w", newline="\n") if args.csv else None
    if out:
        out.write("step,rel_energy_diff,max_rel_coord_diff\n")
    scale = None
    worstE = mpf(0); worstX = mpf(0); last = None; ncommon = 0
    for s in sn:
        w = wide_by_step.get(s["step"])
        if w is None:
            continue
        ncommon += 1
        En, _, pn, vn = invariants(G, mn, s["state"])
        Ew, _, pw, vw = invariants(G, mw, w["state"])
        dE = abs((En - Ew) / Ew) if Ew != 0 else mpf(0)
        if scale is None:
            scale = max(abs(q) for x in pw for q in x)   # a length scale from the wide state
        dX = mpf(0)
        for i in range(len(pn)):
            for c in range(3):
                dX = max(dX, abs(pn[i][c] - pw[i][c]) / scale)
        worstE = max(worstE, dE); worstX = max(worstX, dX); last = (s["step"], dE, dX)
        if out:
            out.write("%d,%s,%s\n" % (s["step"], mpmath.nstr(dE, 6), mpmath.nstr(dX, 6)))
        if not args.quiet:
            print("step %8d  |dE/E| %s  max|dx|/L %s" % (s["step"], mpmath.nstr(dE, 4), mpmath.nstr(dX, 4)))
    if out:
        out.close()
    if last is None:
        print("no common samples")
        return 1
    print("FLOOR %s vs %s: %d common samples, last step %d, max |dE/E| %s (last %s), max |dx|/L %s (last %s)" %
          (hn.get("format"), hw.get("format"), ncommon, last[0], mpmath.nstr(worstE, 4), mpmath.nstr(last[1], 4),
           mpmath.nstr(worstX, 4), mpmath.nstr(last[2], 4)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
