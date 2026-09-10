#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The precision-versus-error measurement: run the port at binary64,
# binary128 and binary256 on the same problems with the same steps,
# score every record with tools/oracle.py, and tabulate.
#
#   python tools/sweep.py [--jobs J] [--only PATTERN] [--out results]
#
# Each job's record goes to results/raw/<name>.rec, its oracle table
# to results/<name>.csv, and one summary line per job to
# results/summary.txt. Jobs run in a pool; the longest first.

import argparse
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
exe = ".exe" if os.name == "nt" else ""
CFT = os.path.join(root, "build", "ias15_cft" + exe)
KEPLER = os.path.join(root, "data", "problems", "kepler.txt")
OUTER = os.path.join(root, "data", "problems", "outer.txt")
T_KEPLER = 6.2800460687587205   # from the oracle, for step counts only

FORMATS = ["fp64", "fp128", "fp256"]


def followup_jobs():
    """Where the first sweep left the question: the wider formats were
    truncation-limited at every step size it used. Smaller steps over
    shorter horizons, to find the binary128 floor and see whether
    binary256 goes under it."""
    J = []
    for dt, orbits, fmts in ((0.0125, 20, FORMATS), (0.00625, 10, ("fp128", "fp256"))):
        steps = int(round(orbits * T_KEPLER / dt))
        for f in fmts:
            J.append(("kepler_fixed_dt%g_%s" % (dt, f), KEPLER, f,
                      ["--dt", str(dt), "--epsilon", "0", "--steps", str(steps), "--sample", str(steps // 20)]))
    for dt, fmts in ((20, ("fp256",)), (10, ("fp128", "fp256"))):
        steps = int(round(300 * 365.25 / dt))
        for f in fmts:
            J.append(("outer_fixed_dt%d_%s_300y" % (dt, f), OUTER, f,
                      ["--dt", str(dt), "--epsilon", "0", "--steps", str(steps), "--sample", str(steps // 15)]))
    return J


def followup2_jobs():
    """The same question with DYADIC step sizes, so that every format
    integrates literally the same step (0.05 rounds differently at
    each format, which pollutes a record-to-record comparison at the
    2^-p level; 1/16, 1/64 and 1/128 do not)."""
    J = []
    for dt, orbits, fmts in ((1.0 / 16, 100, FORMATS), (1.0 / 64, 20, FORMATS), (1.0 / 128, 10, ("fp128", "fp256"))):
        steps = int(round(orbits * T_KEPLER / dt))
        for f in fmts:
            J.append(("kepler_fixed_dt%s_%s" % (repr(dt), f), KEPLER, f,
                      ["--dt", repr(dt), "--epsilon", "0", "--steps", str(steps), "--sample", str(steps // 20)]))
    return J


def jobs():
    J = []
    # Kepler, fixed step, 200 orbits at three step sizes
    for dt in (0.1, 0.05, 0.025):
        steps = int(round(200 * T_KEPLER / dt))
        for f in FORMATS:
            J.append(("kepler_fixed_dt%g_%s" % (dt, f), KEPLER, f,
                      ["--dt", str(dt), "--epsilon", "0", "--steps", str(steps), "--sample", str(steps // 40)]))
    # the augmented-summation variant, at the middle step size
    for f in ("fp64", "fp256"):
        steps = int(round(200 * T_KEPLER / 0.05))
        J.append(("kepler_fixed_dt0.05_%s_aug" % f, KEPLER, f,
                  ["--dt", "0.05", "--epsilon", "0", "--steps", str(steps), "--sample", str(steps // 40), "--cs", "augmented"]))
    # Kepler, adaptive, three tolerances (steps scale like eps^(-1/7))
    for eps, steps in (("1e-9", 10000), ("1e-12", 27000), ("1e-16", 100000)):
        for f in FORMATS:
            if eps == "1e-16" and f == "fp256":
                continue   # ~100k steps at fp256 is several hours; run separately if the others warrant it
            J.append(("kepler_adaptive_eps%s_%s" % (eps, f), KEPLER, f,
                      ["--dt", "0.01", "--epsilon", eps, "--steps", str(steps), "--sample", str(steps // 40)]))
    # outer solar system, fixed step: 1200 years at 40 and 20 days
    for dt, fmts in ((40, FORMATS), (20, ("fp64", "fp128"))):
        steps = int(round(1200 * 365.25 / dt))
        for f in fmts:
            J.append(("outer_fixed_dt%d_%s" % (dt, f), OUTER, f,
                      ["--dt", str(dt), "--epsilon", "0", "--steps", str(steps), "--sample", str(steps // 30)]))
    return J


def cost(job):
    """A rough ordering key: fp256 first, then by step count."""
    name, prob, f, args = job
    steps = int(args[args.index("--steps") + 1])
    w = {"fp64": 1, "fp128": 3, "fp256": 9}[f]
    if "outer" in name:
        w *= 3
    return -steps * w


def run_job(job, outdir):
    name, prob, f, args = job
    rec = os.path.join(outdir, "raw", name + ".rec")
    log = os.path.join(outdir, "raw", name + ".log")
    csv = os.path.join(outdir, name + ".csv")
    t0 = time.time()
    with open(rec, "w") as fo, open(log, "w") as fe:
        p = subprocess.run([CFT, "--format", f, "--problem", prob, "--max-iter", "60", "--quiet"] + args,
                           stdout=fo, stderr=fe, text=True)
    wall = time.time() - t0
    if p.returncode != 0:
        return "%s: FAILED rc=%d after %.0f s (see %s)" % (name, p.returncode, wall, log)
    q = subprocess.run([sys.executable, os.path.join(here, "oracle.py"), rec, "--csv", csv, "--quiet"],
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if q.returncode != 0:
        return "%s: oracle FAILED: %s" % (name, q.stderr[-500:])
    summary = [l for l in q.stdout.splitlines() if l.startswith("SUMMARY")]
    return "%s: wall=%.0fs %s" % (name, wall, summary[0] if summary else q.stdout.strip())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 1))
    ap.add_argument("--only", default=None)
    ap.add_argument("--out", default=os.path.join(root, "results"))
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--set", default="main", choices=["main", "followup", "followup2"])
    ap.add_argument("--exe", default=None, help="binary to run (default build/ias15_cft)")
    args = ap.parse_args()
    global CFT
    if args.exe:
        CFT = args.exe
    J = {"main": jobs, "followup": followup_jobs, "followup2": followup2_jobs}[args.set]()
    if args.only:
        J = [j for j in J if args.only in j[0]]
    J.sort(key=cost)
    if args.list:
        for j in J:
            print(j[0], " ".join(j[3]))
        return 0
    os.makedirs(os.path.join(args.out, "raw"), exist_ok=True)
    summary_path = os.path.join(args.out, "summary.txt")
    print("sweep: %d jobs, %d at a time, records under %s" % (len(J), args.jobs, args.out), flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=args.jobs) as pool, open(summary_path, "a") as sf:
        futs = {pool.submit(run_job, j, args.out): j for j in J}
        for fut in futs:
            pass
        done = 0
        from concurrent.futures import as_completed
        for fut in as_completed(futs):
            line = fut.result()
            done += 1
            print("[%d/%d %.0fs] %s" % (done, len(J), time.time() - t0, line), flush=True)
            sf.write(line + "\n")
            sf.flush()
    print("sweep done in %.0f s" % (time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
