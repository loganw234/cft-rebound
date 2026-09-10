#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The regression gate: the first samples of committed records under
# results/raw/ are recomputed by the current build and must match bit
# for bit. A record's sampling does not perturb its trajectory, so a
# short run on the same sampling grid reproduces the long run's early
# samples exactly - at every format, fixed and adaptive, both problems.
# This is what protects the single-system path when the integrator
# changes shape (the ensemble rewrite, for one).
#
#   python tools/check_records.py [--build build] [--quick]

import argparse
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import parse_record, run  # noqa: E402

root = os.path.dirname(here)


def header_args(text):
    """The run's own arguments, from its header lines."""
    h = {}
    for line in text.splitlines():
        if line.startswith("# program="):
            for kv in line[2:].split():
                k, _, v = kv.partition("=")
                h[k] = v
        elif line.startswith("# dt0="):
            toks = line[2:].split()
            keys = [t for t in toks if t.endswith("=")]
            vals = [t for t in toks if not t.endswith("=") and "=" not in t]
            for k, v in zip(keys, vals):
                h[k[:-1]] = v
    return h


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--exe", default="ias15_cft")
    ap.add_argument("--quick", action="store_true", help="fp64 records only")
    ap.add_argument("--samples", type=int, default=2, help="how many samples after the first to recompute")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, args.exe + exe)
    raw = os.path.join(root, "results", "raw")
    problems = {"kepler": os.path.join(root, "data", "problems", "kepler.txt"),
                "outer": os.path.join(root, "data", "problems", "outer.txt")}
    names = ["kepler_fixed_dt0.05_fp64", "kepler_adaptive_eps1e-9_fp64", "outer_fixed_dt40_fp64", "kepler_fixed_dt0.05_fp64_aug",
             "kepler_fixed_dt0.05_fp128", "kepler_adaptive_eps1e-9_fp128", "outer_fixed_dt40_fp128",
             "kepler_fixed_dt0.05_fp256", "kepler_adaptive_eps1e-9_fp256"]
    if args.quick:
        names = [n for n in names if n.endswith("fp64") or n.endswith("fp64_aug")]
    allok = True
    print("check_records: the first %d samples of the committed records, recomputed" % args.samples)
    for name in names:
        path = os.path.join(raw, name + ".rec")
        text = open(path).read()
        h = header_args(text)
        samples, _ = parse_record(text)
        sample = int(h["sample"])
        steps = sample * args.samples
        # dt and epsilon: the header carries them as exact hex, which the
        # port's parser takes as-is at any format
        cmd = [cft, "--format", h["format"], "--problem", problems[h["problem"]], "--dt", h["dt0"], "--epsilon", h["epsilon"],
               "--steps", str(steps), "--sample", str(sample), "--max-iter", h["max_iter"], "--cs", h.get("cs", "kahan"),
               "--arith", h.get("arith", "rebound"), "--quiet"]
        # the record's dt0 is the decimal's rounding at that format, and
        # from_text on the hex form reproduces those bits exactly
        out = run(cmd)
        if out is None:
            allok = False
            continue
        new, _ = parse_record(out)
        old = [s for s in samples if s[1] <= steps]
        if len(new) != len(old):
            print("  FAIL %s: %d samples recomputed, %d in the record up to step %d" % (name, len(new), len(old), steps))
            allok = False
            continue
        nvals = 0
        bad = None
        for (ok_, ostep, ov), (nk, nstep, nv) in zip(old, new):
            if ostep != nstep or len(ov) != len(nv):
                bad = "sample %d misaligned (%d/%d, %d/%d values)" % (ok_, ostep, nstep, len(ov), len(nv))
                break
            for j, (a, b) in enumerate(zip(ov, nv)):
                nvals += 1
                if a != b:
                    bad = "sample %d (step %d) column %d: record %s now %s" % (ok_, ostep, j, a[0], b[0])
                    break
            if bad:
                break
        if bad:
            print("  FAIL %s: %s" % (name, bad))
            allok = False
        else:
            print("  ok   %s: %d samples, %d values identical" % (name, len(new), nvals))
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
