#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The correctness gate of the port: at binary64, ias15_cft must agree
# with REBOUND's own IAS15 (ias15_ref) BIT FOR BIT on every recorded
# value of every sample - time, next and last step, and the whole state
# of every body - on both problems, with fixed and with adaptive steps.
#
# Values are compared as exact rationals parsed from the hex-float
# records, so "agree" means the same bits (a -0 against a +0 would be
# reported by the sign check below).
#
#   python tools/check_equivalence.py [--build build] [--quick]

import argparse
import os
import subprocess
import sys
from fractions import Fraction

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)


def hex_to_fraction(s):
    """Exact value of a canonical hex float of any width, and its sign
    bit (so that -0 and +0 are told apart)."""
    t = s.strip()
    neg = t.startswith("-")
    if neg or t.startswith("+"):
        t = t[1:]
    if t in ("inf", "nan"):
        return (t, neg)
    if not t[:2].lower() == "0x":
        raise ValueError("not a hex float: %r" % s)
    body = t[2:]
    mant, _, exp = body.partition("p")
    if not exp:
        mant, _, exp = body.partition("P")
    ip, _, fp = mant.partition(".")
    digits = (ip + fp) or "0"
    val = Fraction(int(digits, 16), 16 ** len(fp)) * (Fraction(2) ** int(exp))
    return (-val if neg else val, neg)


def parse_record(text):
    samples = []
    tail = {}
    for line in text.splitlines():
        if line.startswith("sample "):
            toks = line.split()
            body = toks[3:]
            if "|" in body:
                body = body[:body.index("|")]
            samples.append((int(toks[1]), int(toks[2]), [hex_to_fraction(x) for x in body]))
        elif line.startswith("# steps_done="):
            for kv in line[2:].split():
                k, _, v = kv.partition("=")
                tail[k] = v
    return samples, tail


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0:
        print("  command failed (%d): %s" % (p.returncode, " ".join(cmd)))
        print(p.stderr[-2000:])
        return None
    return p.stdout


def compare(name, ref_txt, cft_txt, ncols_per_body, nbodies):
    rs, rt = parse_record(ref_txt)
    cs, ct = parse_record(cft_txt)
    if len(rs) != len(cs):
        print("  FAIL %s: %d samples from ref, %d from cft" % (name, len(rs), len(cs)))
        return False
    labels = ["t", "dt_next", "dt_last", "E"]
    for i in range(nbodies):
        for c in ("x", "y", "z", "vx", "vy", "vz"):
            labels.append("%s[%d]" % (c, i))
    worst = None
    nvals = 0
    for (rk, rstep, rv), (ck, cstep, cv) in zip(rs, cs):
        if rk != ck or rstep != cstep:
            print("  FAIL %s: sample %d/%d step %d/%d misaligned" % (name, rk, ck, rstep, cstep))
            return False
        if len(rv) != len(cv):
            print("  FAIL %s: sample %d has %d values from ref, %d from cft" % (name, rk, len(rv), len(cv)))
            return False
        for j, (a, b) in enumerate(zip(rv, cv)):
            nvals += 1
            if a != b:
                if worst is None:
                    worst = (rk, rstep, labels[j] if j < len(labels) else str(j), a, b)
    if worst:
        k, step, lab, a, b = worst
        av = a[0]; bv = b[0]
        rel = ""
        if isinstance(av, Fraction) and isinstance(bv, Fraction) and av != 0:
            rel = " rel %.3e" % float((bv - av) / av)
        print("  FAIL %s: first difference at sample %d (step %d) %s: ref %s cft %s%s" %
              (name, k, step, lab, float(av) if isinstance(av, Fraction) else av,
               float(bv) if isinstance(bv, Fraction) else bv, rel))
        return False
    print("  ok   %s: %d samples, %d values identical (steps_rejected=%s, max_exceeded=%s)" %
          (name, len(rs), nvals, ct.get("steps_rejected", "?"), ct.get("iterations_max_exceeded", "?")))
    return True


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--cs", default="kahan")
    ap.add_argument("--exe", default="ias15_cft")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    ref = os.path.join(args.build, "ias15_ref" + exe)
    cft = os.path.join(args.build, args.exe + exe)
    kepler = os.path.join(root, "data", "problems", "kepler.txt")
    outer = os.path.join(root, "data", "problems", "outer.txt")
    scale = 1 if not args.quick else 4
    cases = [
        ("kepler fixed dt=0.05",     kepler, 2, ["--dt", "0.05", "--epsilon", "0",    "--steps", str(800 // scale),  "--sample", "50"]),
        ("kepler adaptive eps=1e-9", kepler, 2, ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(800 // scale),  "--sample", "50"]),
        ("outer fixed dt=40",        outer,  6, ["--dt", "40",   "--epsilon", "0",    "--steps", str(200 // scale),  "--sample", "25"]),
        ("outer adaptive eps=1e-9",  outer,  6, ["--dt", "40",   "--epsilon", "1e-9", "--steps", str(200 // scale),  "--sample", "25"]),
        # a far too large first step forces the rejected-step path and
        # REBOUND's cap of 12 corrector passes, once each
        ("kepler adaptive from dt0=1.0 (rejection)",  kepler, 2, ["--dt", "1.0",  "--epsilon", "1e-9", "--steps", str(300 // scale),  "--sample", "50"]),
        ("outer adaptive from dt0=2000 (rejection)",  outer,  6, ["--dt", "2000", "--epsilon", "1e-9", "--steps", str(100 // scale),  "--sample", "25"]),
    ]
    allok = True
    print("check_equivalence: ias15_cft --format fp64 --cs %s against REBOUND's IAS15" % args.cs)
    for name, prob, nb, extra in cases:
        r = run([ref, "--problem", prob] + extra)
        c = run([cft, "--format", "fp64", "--problem", prob, "--cs", args.cs, "--quiet"] + extra)
        if r is None or c is None:
            allok = False
            continue
        allok &= compare(name, r, c, 6, nb)
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
