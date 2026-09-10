#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The correctness gate of the ensemble mode: an ensemble of E systems
# integrated in one run must reproduce, BIT FOR BIT, the records of
# those E systems run one at a time (--member k) - every recorded value
# of every sample of every member, and each member's step statistics
# (steps rejected, corrector passes, cap hits). With fixed and with
# per-system adaptive steps, on the loop and the program engines, at
# every format, with a mixed attempt in which some members reject
# their step while others accept it.
#
#   python tools/check_ensemble.py [--build build] [--formats fp64,fp128,fp256]
#                                  [--quick] [--exe ias15_cft] [--keep]

import argparse
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import hex_to_fraction, run  # noqa: E402

root = os.path.dirname(here)


def parse(text, system=None):
    """Sample lines of one system (or of a single-system record) as
    (sample, step, values), and its trailer as a dict."""
    samples = []
    tail = {}
    want = None if system is None else "system %d sample " % system
    for line in text.splitlines():
        if want is None and line.startswith("sample "):
            toks = line.split()[1:]
        elif want is not None and line.startswith(want):
            toks = line.split()[3:]
        elif line.startswith("# steps_done=") and system is None:
            for kv in line[2:].split():
                k, _, v = kv.partition("=")
                tail[k] = v
            continue
        elif system is not None and line.startswith("# system %d steps_done=" % system):
            for kv in line.split()[3:]:
                k, _, v = kv.partition("=")
                tail[k] = v
            continue
        else:
            continue
        body = toks[2:]
        if "|" in body:
            i = body.index("|")
            body = body[:i] + body[i + 1:]     # the exact-time pair is compared too
        samples.append((int(toks[0]), int(toks[1]), [hex_to_fraction(x) for x in body]))
    return samples, tail


def compare(name, ens_txt, member_txt, k):
    es, et = parse(ens_txt, k)
    ms, mt = parse(member_txt)
    if not es or not ms:
        print("  FAIL %s member %d: %d ensemble samples, %d member samples" % (name, k, len(es), len(ms)))
        return False
    if len(es) != len(ms):
        print("  FAIL %s member %d: %d samples in the ensemble, %d alone" % (name, k, len(es), len(ms)))
        return False
    nvals = 0
    for (ek, estep, ev), (mk, mstep, mv) in zip(es, ms):
        if ek != mk or estep != mstep:
            print("  FAIL %s member %d: sample %d/%d step %d/%d misaligned" % (name, k, ek, mk, estep, mstep))
            return False
        if len(ev) != len(mv):
            print("  FAIL %s member %d: sample %d has %d values in the ensemble, %d alone" % (name, k, ek, len(ev), len(mv)))
            return False
        for j, (a, b) in enumerate(zip(ev, mv)):
            nvals += 1
            if a != b:
                print("  FAIL %s member %d: first difference at sample %d (step %d) column %d: ensemble %s alone %s" %
                      (name, k, ek, estep, j, a[0], b[0]))
                return False
    for key in ("steps_done", "steps_rejected", "iterations_max_exceeded", "mean_pc_iterations", "max_pc_iterations"):
        if et.get(key) != mt.get(key):
            print("  FAIL %s member %d: trailer %s is %s in the ensemble, %s alone" % (name, k, key, et.get(key), mt.get(key)))
            return False
    return nvals, et


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--exe", default="ias15_cft")
    ap.add_argument("--formats", default="fp64,fp128,fp256")
    ap.add_argument("--quick", action="store_true", help="fp64 only, fewer steps")
    ap.add_argument("--keep", default=None, help="directory for the ensemble files (default build/ensembles)")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, args.exe + exe)
    py = sys.executable
    progs = os.path.join(root, "programs", "out")
    kepler = os.path.join(root, "data", "problems", "kepler.txt")
    outer = os.path.join(root, "data", "problems", "outer.txt")
    ensdir = args.keep or os.path.join(args.build, "ensembles")
    os.makedirs(ensdir, exist_ok=True)
    formats = ["fp64"] if args.quick else args.formats.split(",")
    scale = 4 if args.quick else 1

    # the ensembles: a Kepler family whose members have different
    # velocities (so different eccentricities, periods and step
    # sequences), a Kepler family a few ulps apart, and an outer solar
    # system family with Jupiter displaced by a small dyadic offset
    kep_v = os.path.join(ensdir, "kepler_vscale_E5.txt")
    kep_u = os.path.join(ensdir, "kepler_ulps_E4.txt")
    kep_s = os.path.join(ensdir, "kepler_spread_E5.txt")
    out_j = os.path.join(ensdir, "outer_jup_E3.txt")
    gen = [
        [py, os.path.join(here, "make_ensemble.py"), kepler, "--members", "5", "--out", kep_v, "--vscale", "planet", "-0.08"],
        [py, os.path.join(here, "make_ensemble.py"), kepler, "--members", "4", "--out", kep_u, "--offset", "planet", "x", "0x1p-40"],
        # the planet started at 0.5, 0.9, ... 2.1 with the pericentre speed
        # of the base orbit: an 8.6x spread in natural step, so that a first
        # step of 0.6 is rejected by the inner members and accepted by the
        # outermost in the same attempt
        [py, os.path.join(here, "make_ensemble.py"), kepler, "--members", "5", "--out", kep_s, "--offset", "planet", "x", "0.4"],
        [py, os.path.join(here, "make_ensemble.py"), outer, "--members", "3", "--out", out_j, "--offset", "jupiter", "x", "0x1p-20"],
    ]
    for g in gen:
        if run(g) is None:
            print("RESULT: FAIL (could not write an ensemble file)")
            return 1

    # (name, file, E, extra args, engine, wants a mixed accept/reject attempt)
    cases = [
        ("kepler vscale E=5 fixed dt=0.05",             kep_v, 5, ["--dt", "0.05", "--epsilon", "0",    "--steps", str(200 // scale), "--sample", "50"],  "loop", False),
        ("kepler vscale E=5 adaptive eps=1e-9",         kep_v, 5, ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(200 // scale), "--sample", "50"],  "loop", False),
        ("kepler spread E=5 adaptive from dt0=0.6 (mixed rejection)", kep_s, 5, ["--dt", "0.6", "--epsilon", "1e-9", "--steps", str(120 // scale), "--sample", "40"], "loop", True),
        ("kepler spread E=5 adaptive from dt0=0.6, program engine (mixed rejection)", kep_s, 5, ["--dt", "0.6", "--epsilon", "1e-9", "--steps", str(120 // scale), "--sample", "40"], "program", True),
        ("kepler ulps E=4 adaptive eps=1e-9",           kep_u, 4, ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(200 // scale), "--sample", "50"],  "loop", False),
        ("kepler vscale E=5 adaptive eps=1e-9, program engine", kep_v, 5, ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(120 // scale), "--sample", "40"], "program", False),
        ("kepler vscale E=5 fixed dt=0.05, program engine",     kep_v, 5, ["--dt", "0.05", "--epsilon", "0",    "--steps", str(120 // scale), "--sample", "40"], "program", False),
        ("outer jupiter E=3 fixed dt=40",               out_j, 3, ["--dt", "40",   "--epsilon", "0",    "--steps", str(60 // scale),  "--sample", "20"],  "loop", False),
        ("outer jupiter E=3 adaptive eps=1e-9",         out_j, 3, ["--dt", "40",   "--epsilon", "1e-9", "--steps", str(60 // scale),  "--sample", "20"],  "loop", False),
    ]
    allok = True
    for fmt in formats:
        print("check_ensemble at %s: the ensemble against its members run alone" % fmt)
        for name, prob, E, extra, engine, mixed in cases:
            common = [cft, "--format", fmt, "--problem", prob, "--max-iter", "60", "--quiet"] + extra
            if engine == "program":
                common += ["--arith", "fma", "--engine", "program", "--programs", progs]
            ens = run(common)
            if ens is None:
                allok = False
                continue
            ok = True
            total = 0
            rejected = []
            for k in range(E):
                mem = run(common + ["--member", str(k)])
                if mem is None:
                    ok = False
                    break
                r = compare(name, ens, mem, k)
                if r is False:
                    ok = False
                    break
                total += r[0]
                rejected.append(int(r[1].get("steps_rejected", "0")))
            if ok and mixed and not (any(rejected) and not all(rejected)):
                print("  FAIL %s: expected a mixed attempt (some members rejecting, some accepting) but steps_rejected=%s" % (name, rejected))
                ok = False
            if ok:
                eff = ""
                for line in ens.splitlines():
                    if line.startswith("# steps_done="):
                        for kv in line[2:].split():
                            if kv.startswith("pc_lane_efficiency=") or kv.startswith("mean_pc_iterations="):
                                eff += " " + kv
                print("  ok   %s: %d members, %d values identical, steps_rejected=%s,%s" % (name, E, total, rejected, eff))
            allok &= ok
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
