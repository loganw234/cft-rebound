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
# predict_next_step is a program as well (nextstep-<fmt>, once a step),
# and two more things are held for it beside the records:
#   - that it RAN. The trailer's nextstep_calls is zero on the loop
#     engine and, on the program engine, one for every step attempt that
#     predicts: at least the steps some system accepted, at most every
#     attempt. Records compared between two engines that both skipped it
#     would agree just as well.
#   - that each case gives it the step ratio the case is there for. A
#     fixed step's ratio is exactly 1, so on a fixed case the program's
#     powers of q are all 1 and a fault in them is inert by construction
#     (METHOD.md 5); an adaptive case must show a ratio that is not 1.
#     And the ratio is per SYSTEM, since each keeps its own step: a
#     program that took it from the constant bank would be right for one
#     system, and for an ensemble whose systems happen to share a ratio,
#     and wrong otherwise. Only the ensemble case can see that, so it must
#     show its systems' ratios differing.
# The trailer's step counters and flags_seen must also be the loop's;
# nothing else compares them.
#
#   python tools/check_program_engine.py [--build build] [--exe ias15_cft]
#                                        [--formats fp64,fp128,fp256]
#   (`make check-quick` passes --formats fp64)

import argparse
import os
import sys
from fractions import Fraction

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import parse_record, compare, run, hex_to_fraction  # noqa: E402

root = os.path.dirname(here)

# the trailer's step counters: the two engines must agree on every one
COUNTERS = ("steps_done", "steps_rejected", "iterations_max_exceeded",
            "mean_pc_iterations", "max_pc_iterations", "flags_seen")


def trailers(text):
    """The run's trailer, and each system's (an ensemble prints a line
    per system too), as dicts of the printed text."""
    whole, systems = {}, {}
    for line in text.splitlines():
        if line.startswith("# steps_done="):
            whole = dict(kv.partition("=")[::2] for kv in line[2:].split())
        elif line.startswith("# system ") and " steps_done=" in line:
            toks = line.split()
            systems[int(toks[2])] = dict(kv.partition("=")[::2] for kv in toks[3:])
    return whole, systems


def sampled_ratios(text):
    """For every sample, each system's dt_next / dt_last, exactly. That is
    the ratio predict_next_step was given on the system's last attempt,
    which was accepted: a system is sampled when its block of steps is
    done. None where dt_last is zero, before the first step."""
    out = {}
    for line in text.splitlines():
        toks = line.split()
        if line.startswith("sample "):
            s, k, dn, dl = 0, int(toks[1]), toks[4], toks[5]
        elif line.startswith("system ") and len(toks) > 7 and toks[2] == "sample":
            s, k, dn, dl = int(toks[1]), int(toks[3]), toks[6], toks[7]
        else:
            continue
        num, den = hex_to_fraction(dn)[0], hex_to_fraction(dl)[0]
        out.setdefault(k, {})[s] = None if den == 0 else num / den
    return [out[k] for k in sorted(out)]


def compare_lines(name, loop_txt, prog_txt, nsys):
    """An ensemble's records, every line, as text: canonical hex, so the
    same text is the same bits. check_equivalence's compare() reads only
    single-system lines and would find nothing here to compare."""
    a = [l for l in loop_txt.splitlines() if l and not l.startswith("#")]
    b = [l for l in prog_txt.splitlines() if l and not l.startswith("#")]
    if not a or len(a) % nsys or len(a) != len(b):
        print("  FAIL %s: %d record lines on the loop, %d on the program engine, for %d systems" % (name, len(a), len(b), nsys))
        return False
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            xs, ys = x.split(), y.split()
            j = next((j for j, (p, q) in enumerate(zip(xs, ys)) if p != q), min(len(xs), len(ys)))
            print("  FAIL %s: first difference in record line %d (%s), token %d: loop %s, program %s" %
                  (name, i, " ".join(xs[:4]), j, xs[j] if j < len(xs) else "(none)", ys[j] if j < len(ys) else "(none)"))
            return False
    print("  ok   %s: %d record lines identical, %d systems x %d samples" % (name, len(a), nsys, len(a) // nsys))
    return True


def held(name, loop_txt, prog_txt, nsys, kind):
    """Everything beside the records; prints what it held, or a FAIL."""
    lw, ls = trailers(loop_txt)
    pw, ps = trailers(prog_txt)
    for k in COUNTERS:
        if lw.get(k) is None or lw.get(k) != pw.get(k):
            print("  FAIL %s: trailer %s is %s on the loop engine, %s on the program engine" % (name, k, lw.get(k), pw.get(k)))
            return False
    if nsys > 1:
        if sorted(ls) != list(range(nsys)) or ls != ps:
            print("  FAIL %s: the %d systems' trailer lines differ between the engines" % (name, nsys))
            return False
        done = [int(ls[s]["steps_done"]) for s in range(nsys)]
        rejected = [int(ls[s]["steps_rejected"]) for s in range(nsys)]
    else:
        done, rejected = [int(lw["steps_done"])], [int(lw["steps_rejected"])]
    if "nextstep_calls" not in lw or "nextstep_calls" not in pw:
        print("  FAIL %s: a trailer has no nextstep_calls" % name)
        return False
    nl, npr = int(lw["nextstep_calls"]), int(pw["nextstep_calls"])
    lo, hi = max(done), sum(done) + sum(rejected)
    if nl != 0 or not lo <= npr <= hi:
        print("  FAIL %s: nextstep_calls is %d on the loop engine (must be 0) and %d on the program engine (must be %d to %d)" %
              (name, nl, npr, lo, hi))
        return False
    said = ["counters and flags the loop's", "nextstep ran %d times" % npr]
    rs = [r for r in sampled_ratios(loop_txt) if r and all(v is not None for v in r.values())]
    if not rs:
        print("  FAIL %s: no sampled step carries a ratio" % name)
        return False
    if kind == "fixed":
        if any(v != 1 for r in rs for v in r.values()):
            print("  FAIL %s: a fixed step's ratio is not 1" % name)
            return False
        said.append("ratio 1 at all %d sampled steps, so q is inert here" % len(rs))
    elif kind == "adaptive":
        n = sum(1 for r in rs if any(v != 1 for v in r.values()))
        if n == 0:
            print("  FAIL %s: every sampled ratio is 1, so nothing here exercises q" % name)
            return False
        said.append("ratio not 1 at %d of %d sampled steps" % (n, len(rs)))
    else:
        # the systems' ratios must differ by more than any rounding could
        # hide: then the formats' roundings of them differ as well
        n, last = 0, None
        for r in rs:
            v = sorted(r.values())
            if v[-1] - v[0] > v[0] * Fraction(1, 2 ** 40):
                n, last = n + 1, (v[0], v[-1])
        if n == 0:
            print("  FAIL %s: at no sampled step do the %d systems' ratios differ" % (name, nsys))
            return False
        said.append("the systems' ratios differ at %d of %d sampled steps (%.6f to %.6f at the last)" %
                    (n, len(rs), float(last[0]), float(last[1])))
        if not any(rejected) or all(rejected):
            print("  FAIL %s: expected a mixed attempt, some systems rejecting and some accepting, but steps_rejected=%s" % (name, rejected))
            return False
        said.append("a mixed attempt: steps_rejected=%s" % rejected)
    print("       %s" % "; ".join(said))
    return True


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
    # The per-system ratio's case (see the top): tools/check_ensemble.py's
    # "spread" Kepler family, its planet starting at 0.5, 0.9, ... 2.1, so
    # an 8.6x spread in natural step. From a first step of 0.6 the four
    # inner systems reject and the outermost accepts, so its first attempt
    # also runs the program over systems that predict nothing (their first
    # step, rejected) beside one that predicts, and keeps only the latter.
    ensdir = os.path.join(args.build, "program_engine")
    os.makedirs(ensdir, exist_ok=True)
    spread = os.path.join(ensdir, "kepler_spread_E5.txt")
    if run([sys.executable, os.path.join(here, "make_ensemble.py"), kepler, "--members", "5", "--out", spread,
            "--offset", "planet", "x", "0.4"]) is None:
        print("RESULT: FAIL (could not write the ensemble file)")
        return 1
    # (name, problem, bodies, systems, what the ratio must show, arguments).
    # None is scaled in any tier: this script takes no --light, and the
    # ensemble's step count is part of what it asserts.
    cases = [
        ("kepler fixed dt=0.05",     kepler, 2, 1, "fixed",    ["--dt", "0.05", "--epsilon", "0",    "--steps", "200", "--sample", "50"]),
        ("kepler adaptive eps=1e-9", kepler, 2, 1, "adaptive", ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "200", "--sample", "50"]),
        ("outer fixed dt=40",        outer,  6, 1, "fixed",    ["--dt", "40",   "--epsilon", "0",    "--steps", "60",  "--sample", "20"]),
        ("kepler spread E=5 adaptive from dt0=0.6 (mixed rejection)", spread, 2, 5, "systems",
         ["--dt", "0.6", "--epsilon", "1e-9", "--steps", "120", "--sample", "40"]),
    ]
    allok = True
    for fmt in args.formats.split(","):
        print("check_program_engine at %s: --engine program against --engine loop, both --arith fma" % fmt)
        for name, prob, nb, nsys, kind, extra in cases:
            common = [cft, "--format", fmt, "--problem", prob, "--arith", "fma", "--max-iter", "60", "--quiet"] + extra
            loop = run(common + ["--engine", "loop"])
            prog = run(common + ["--engine", "program", "--programs", progs])
            if loop is None or prog is None:
                allok = False
                continue
            ok = compare(name, loop, prog, 6, nb) if nsys == 1 else compare_lines(name, loop, prog, nsys)
            ok = ok and held(name, loop, prog, nsys, kind)
            allok &= ok
            if nsys > 1:
                continue
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
