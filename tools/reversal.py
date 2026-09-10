#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The reversal test: integrate N fixed steps forward, restart from the
# exact recorded state, integrate N steps back with -dt, and compare
# the returned state with the initial one - bit for bit, and as a
# relative deviation.
#
# IAS15 is not time-symmetric: its Gauss-Radau nodes include the start
# of the step and not the end, so the backward polynomial is a
# different polynomial and even exact arithmetic would not return
# exactly. What the test measures is therefore the method's own
# reversal error at a format where the arithmetic is far below it,
# and the arithmetic's contribution where it is not. A bit-exact
# return would need a symmetric (or a lattice, JANUS-style) scheme;
# what IS exact here is that the run repeats bit for bit.
#
#   python tools/reversal.py --problem FILE --format F --dt DT --steps N [--build build] [--keep DIR]

import argparse
import os
import subprocess
import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import hex_to_fraction  # noqa: E402

root = os.path.dirname(here)
mp.dps = 60


def read_problem(path):
    G, bodies = None, []
    for line in open(path):
        toks = line.split()
        if not toks or toks[0] == "#":
            continue
        if toks[0] == "G":
            G = hex_to_fraction(toks[1])[0]
        elif toks[0] == "body":
            bodies.append((toks[1], [hex_to_fraction(t)[0] for t in toks[2:9]]))
    return G, bodies


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        print("  command failed (%d): %s\n%s" % (p.returncode, " ".join(cmd), p.stderr[-1500:]))
        sys.exit(1)
    return p.stdout


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--problem", required=True)
    ap.add_argument("--format", default="fp64")
    ap.add_argument("--dt", required=True)
    ap.add_argument("--steps", type=int, required=True)
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--keep", default=None)
    ap.add_argument("--extra", default="", help="extra ias15_cft arguments, e.g. '--arith fma'")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, "ias15_cft" + exe)
    keep = args.keep or os.path.join(args.build, "reversal")
    os.makedirs(keep, exist_ok=True)
    tag = "%s_%s_dt%s_n%d" % (os.path.splitext(os.path.basename(args.problem))[0], args.format, args.dt.replace("/", "_"), args.steps)
    extra = args.extra.split()
    common = [cft, "--format", args.format, "--epsilon", "0", "--steps", str(args.steps), "--sample", str(args.steps), "--max-iter", "60", "--quiet"] + extra
    fwd = run(common + ["--problem", args.problem, "--dt", args.dt])
    fwd_rec = os.path.join(keep, tag + "_fwd.rec")
    open(fwd_rec, "w", newline="\n").write(fwd)
    turn = os.path.join(keep, tag + "_turn.txt")
    run([sys.executable, os.path.join(here, "state_to_problem.py"), fwd_rec, "--out", turn, "--last"])
    back = run(common + ["--problem", turn, "--dt", "-" + args.dt.lstrip("+")])
    back_rec = os.path.join(keep, tag + "_back.rec")
    open(back_rec, "w", newline="\n").write(back)
    ret = os.path.join(keep, tag + "_return.txt")
    run([sys.executable, os.path.join(here, "state_to_problem.py"), back_rec, "--out", ret, "--last"])

    G, b0 = read_problem(args.problem)
    _, b1 = read_problem(ret)
    identical = all(x == y for (_, u), (_, v) in zip(b0, b1) for x, y in zip(u, v))
    L = max(abs(q) for _, u in b0 for q in u[1:4])
    V = max(abs(q) for _, u in b0 for q in u[4:7]) or Fraction(1)
    dx = max(abs(u[c] - v[c]) for (_, u), (_, v) in zip(b0, b1) for c in (1, 2, 3))
    dv = max(abs(u[c] - v[c]) for (_, u), (_, v) in zip(b0, b1) for c in (4, 5, 6))
    # the forward run's own error against the closed form is not needed
    # here; the energy at the return point against the initial energy is
    def energy(bodies):
        ek = mpf(0); ep = mpf(0)
        for i, (_, u) in enumerate(bodies):
            m = mpf(u[0].numerator) / u[0].denominator
            ek += m * sum((mpf(q.numerator) / q.denominator) ** 2 for q in u[4:7]) / 2
            for j in range(i + 1, len(bodies)):
                w = bodies[j][1]
                d = mpmath.sqrt(sum((mpf((u[c] - w[c]).numerator) / (u[c] - w[c]).denominator) ** 2 for c in (1, 2, 3)))
                ep -= (mpf(G.numerator) / G.denominator) * m * (mpf(w[0].numerator) / w[0].denominator) / d
        return ek + ep
    E0 = energy(b0); E1 = energy(b1)
    print("REVERSAL problem=%s format=%s dt=%s steps=%d bits_identical=%s max|dx|/L=%s max|dv|/V=%s |dE/E|=%s" %
          (os.path.basename(args.problem), args.format, args.dt, args.steps, "yes" if identical else "no",
           mpmath.nstr(mpf(dx.numerator) / dx.denominator / (mpf(L.numerator) / L.denominator), 4),
           mpmath.nstr(mpf(dv.numerator) / dv.denominator / (mpf(V.numerator) / V.denominator), 4),
           mpmath.nstr(abs((E1 - E0) / E0), 4)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
