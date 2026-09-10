#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The body-count gate. `ias15_cft` refused more than 64 bodies for no
# reason in the code; the cap is now CFT_MAX_BODIES (1024) and this
# checks that the number means something:
#
#   1. at a few hundred bodies the binary64 port is still REBOUND's own
#      IAS15 bit for bit, fixed step and adaptive. This is the same gate
#      as tools/check_equivalence.py, run where it could actually break:
#      gravity is summed per particle over its partners in ascending
#      order, and if that order drifted from REBOUND's (i, j<i) loop it
#      would show as a last-bit difference at large N and at no other N.
#   2. the cap itself is reachable - a problem at exactly CFT_MAX_BODIES
#      is read, allocated and evaluated, not just accepted;
#   3. one body past it is refused with a message that says the limit.
#
# The problems come from tools/make_nbody.py, which is deterministic, so
# nothing large is committed to the repository.
#
# Both the memory and the time of this integrator grow as E*N^2, so the
# default here is deliberately small - N = 256, two steps, about a
# minute - and --quick is one step. `--n 512 --steps 1` is the same gate
# where the ROADMAP said the port had been exercised; it costs about
# 50 s a step at binary64 on the host this was written on, which is why
# it is not the default.
#
#   python tools/check_bodycount.py [--build build] [--quick]
#                                   [--n N] [--steps K]

import argparse
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
sys.path.insert(0, here)

from check_equivalence import compare, run   # noqa: E402


def cap_from_source():
    """CFT_MAX_BODIES as src/ias15_cft.c defines it - read, not
    transcribed, so this test cannot drift from the program."""
    path = os.path.join(root, "src", "ias15_cft.c")
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#define CFT_MAX_BODIES"):
                return int(line.split()[2])
    raise SystemExit("check_bodycount: no #define CFT_MAX_BODIES in src/ias15_cft.c")


def make_problem(n, outdir):
    path = os.path.join(outdir, "nbody%d.txt" % n)
    if not os.path.exists(path):
        p = subprocess.run([sys.executable, os.path.join(here, "make_nbody.py"), str(n), path],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if p.returncode != 0:
            print(p.stderr[-2000:])
            raise SystemExit("check_bodycount: make_nbody.py failed")
    return path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--n", type=int, default=None)
    ap.add_argument("--steps", type=int, default=None)
    args = ap.parse_args()

    exe = ".exe" if os.name == "nt" else ""
    ref = os.path.join(args.build, "ias15_ref" + exe)
    cft = os.path.join(args.build, "ias15_cft" + exe)
    outdir = os.path.join(args.build, "problems")
    os.makedirs(outdir, exist_ok=True)

    cap = cap_from_source()
    n = args.n if args.n is not None else 256
    steps = args.steps if args.steps is not None else (1 if args.quick else 2)
    print("check_bodycount: cap = %d, equivalence at N = %d, %d step(s)" % (cap, n, steps))

    prob = make_problem(n, outdir)
    allok = True
    cases = [
        ("nbody%d fixed dt=0.005" % n, ["--dt", "0.005", "--epsilon", "0"]),
        ("nbody%d adaptive eps=1e-9" % n, ["--dt", "0.005", "--epsilon", "1e-9"]),
    ]
    for name, extra in cases:
        tail = ["--problem", prob, "--steps", str(steps), "--sample", "1"] + extra
        r = run([ref] + tail)
        c = run([cft, "--format", "fp64", "--quiet"] + tail)
        if r is None or c is None:
            allok = False
            continue
        allok &= compare(name, r, c, 6, n)

    # the cap is reachable: read, allocate and evaluate at exactly N = cap.
    # --steps 0 writes the initial sample and stops, which is the whole
    # allocation and one energy evaluation over all N(N-1)/2 pairs.
    big = make_problem(cap, outdir)
    p = subprocess.run([cft, "--format", "fp64", "--quiet", "--problem", big,
                        "--dt", "0.005", "--epsilon", "0", "--steps", "0", "--sample", "1"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0 or "sample 0 0" not in p.stdout:
        print("  FAIL N = %d (the cap) did not run: exit %d\n%s" % (cap, p.returncode, p.stderr[-2000:]))
        allok = False
    else:
        print("  ok   N = %d (the cap): allocated and evaluated" % cap)

    # one past it is refused, and the message says the limit
    over = make_problem(cap + 1, outdir)
    p = subprocess.run([cft, "--format", "fp64", "--quiet", "--problem", over,
                        "--dt", "0.005", "--epsilon", "0", "--steps", "0", "--sample", "1"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    msg = (p.stderr or "") + (p.stdout or "")
    if p.returncode == 0:
        print("  FAIL N = %d was accepted; the cap is %d" % (cap + 1, cap))
        allok = False
    elif ("1..%d" % cap) not in msg:
        print("  FAIL N = %d was refused, but the message does not name 1..%d:\n    %s"
              % (cap + 1, cap, msg.strip()[:300]))
        allok = False
    else:
        print("  ok   N = %d refused: %s" % (cap + 1, msg.strip().splitlines()[-1][:160]))

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
