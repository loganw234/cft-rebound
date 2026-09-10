#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of the sequencer-program prototype: with --arith fma, the
# predictor and corrector run as orbit-sequencer programs
# (--engine program, tools/gen_programs.py) must produce the SAME BITS
# as the host loop issuing the same operations one cft_run at a time
# (--engine loop), at every format, on both problems, with fixed and
# adaptive steps. Also records how far the FMA form drifts from
# REBOUND's own rounding sequence, which it must - it is a different
# sequence - but only at the round-off level.
#
#   python tools/check_program_engine.py [--build build] [--exe ias15_cft]
#                                        [--formats fp64,fp128,fp256]
#   (`make check-quick` passes --formats fp64)

import argparse
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import parse_record, compare, run  # noqa: E402

root = os.path.dirname(here)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--exe", default="ias15_cft")
    ap.add_argument("--formats", default="fp64,fp128,fp256")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, args.exe + exe)
    kepler = os.path.join(root, "data", "problems", "kepler.txt")
    outer = os.path.join(root, "data", "problems", "outer.txt")
    progs = os.path.join(root, "programs", "out")
    cases = [
        ("kepler fixed dt=0.05",     kepler, 2, ["--dt", "0.05", "--epsilon", "0",    "--steps", "200", "--sample", "50"]),
        ("kepler adaptive eps=1e-9", kepler, 2, ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "200", "--sample", "50"]),
        ("outer fixed dt=40",        outer,  6, ["--dt", "40",   "--epsilon", "0",    "--steps", "60",  "--sample", "20"]),
    ]
    allok = True
    for fmt in args.formats.split(","):
        print("check_program_engine at %s: --engine program against --engine loop, both --arith fma" % fmt)
        for name, prob, nb, extra in cases:
            common = [cft, "--format", fmt, "--problem", prob, "--arith", "fma", "--max-iter", "60", "--quiet"] + extra
            loop = run(common + ["--engine", "loop"])
            prog = run(common + ["--engine", "program", "--programs", progs])
            if loop is None or prog is None:
                allok = False
                continue
            allok &= compare(name, loop, prog, 6, nb)
            # the FMA form against REBOUND's sequence: must differ (it is a different rounding
            # sequence) and only at round-off level - report the first sample's relative gap
            reb = run([cft, "--format", fmt, "--problem", prob, "--arith", "rebound", "--max-iter", "60", "--quiet"] + extra)
            if reb is not None:
                rs, _ = parse_record(reb)
                ls, _ = parse_record(loop)
                worst = 0.0
                for (_, _, rv), (_, _, lv) in zip(rs, ls):
                    for a, b in zip(rv[4:], lv[4:]):
                        if a[0] != 0:
                            worst = max(worst, abs(float((b[0] - a[0]) / a[0])))
                print("       fma form vs REBOUND's sequence: largest relative difference over the run %.3e" % worst)
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
