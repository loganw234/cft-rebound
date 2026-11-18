#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# What this port's default path costs on a GPU (round 5, P1): a call and a
# step, at each format, on both engines, on hip:N and on the software
# backend of the same host, from the trailer's wall_seconds.
#
#   python3 hw/gpu/costs.py [--build build] [--artifact hip:0]
#                           [--formats fp64,fp128,fp256] [--problems kepler,ring]
#
# MEASURED, one run a row, and unranked within CLAUDE.md's noise (two runs
# of one configuration differ by up to 6% on the small problems): a row says
# what one run took, not which route is faster by less than that.
#
# The clocks. wall_seconds is the program's own, from its start to its
# trailer, and is the one used here. The trailer's seconds and steps_per_s
# are CPU time, clock(), which the HIP runtime's own threads make larger
# than the wall time on a GPU (V0 measured 1.5 to 2.0 over 24 runs); the
# CPU-to-wall ratio is printed beside each GPU row so that no one reads them
# as time. A call is wall_seconds over the trailer's calls= (every library
# call the port makes, the host's augmented sums among them); a step is
# wall_seconds over steps_done. Both include the energy at each sample and
# the start-up the trailer's clock covers.
#
# Every GPU row's record is held to its software row's (the header's
# backend=, the trailer's timings and the by= of a route the device does not
# publish set aside): a cost of a run that did not compute the same thing
# would be a cost of something else.
#
# The problems:
#   kepler: data/problems/kepler.txt, --dt 0.01 --epsilon 1e-9 --steps 50
#           --sample 10, P0's smoke (VALIDATION entry 54), so its rows meet
#           P0's second time;
#   ring:   one wave of tools/make_ring.py's ring, Saturn with its zonal
#           line, six moons and 57 test particles (64 bodies, GPU-E2's
#           shape a2, what P2's kernel replaces), --dt 100 --epsilon 1e-9,
#           8 steps at binary64, 6 at binary128 and 4 at binary256.
# Each on the loop engine in REBOUND's order and on the program engine in
# the FMA form (`--engine program --arith fma`).

import argparse
import os
import re
import subprocess
import sys
import tempfile
import time

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(os.path.dirname(here))
EXE = ".exe" if os.name == "nt" else ""


def run(cmd, env):
    t0 = time.time()
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env, cwd=root,
                       stdin=subprocess.DEVNULL)
    return p.returncode, p.stdout, p.stderr, time.time() - t0


def trailer(text):
    m = re.search(r"^# steps_done=.*$", text, re.M)
    return dict(t.partition("=")[::2] for t in m.group(0)[2:].split()) if m else {}


def norm(text):
    out = []
    for l in text.splitlines():
        if l.startswith("# program="):
            l = re.sub(r" backend=\S+", " backend=*", l)
        elif l.startswith("# steps_done="):
            l = re.sub(r" (?:seconds|steps_per_s|wall_seconds|system_steps_per_s|t_[a-z_]+)=\S+", "", l)
        l = re.sub(r"\b(indexed|pair_gather)=off by=(absent|opt_in)\b", r"\1=off by=~", l)
        out.append(l)
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--artifact", default="hip:0")
    ap.add_argument("--formats", default="fp64,fp128,fp256")
    ap.add_argument("--problems", default="kepler,ring")
    ap.add_argument("--engines", default="loop,program")
    ap.add_argument("--work", default=None)
    args = ap.parse_args()
    real = os.path.join(os.path.abspath(args.build), "ias15_cft" + EXE)
    art = args.artifact
    work = args.work or tempfile.mkdtemp(prefix="costs-")
    os.makedirs(work, exist_ok=True)
    base = {k: v for k, v in os.environ.items() if not k.upper().startswith("CFT_")}
    env_dev = dict(base, CFT_REBOUND_ARTIFACT=art, CFT_HIP_INFO="1")
    env_sw = dict(base)
    ring = os.path.join(work, "ring1w_z.txt")
    if "ring" in args.problems:
        rc, o, e, _ = run([sys.executable, "-B", os.path.join(root, "tools", "make_ring.py"), ring, "--one-system",
                           "--particles", "57", "--seed", "1"], env_sw)
        if rc:
            print("costs: make_ring.py failed: %s" % (e or o)[-300:])
            return 2
    problems = {
        "kepler": (os.path.join(root, "data", "problems", "kepler.txt"),
                   {f: ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "50", "--sample", "10"]
                    for f in ("fp64", "fp128", "fp256")}),
        "ring": (ring, {"fp64": ["--dt", "100", "--epsilon", "1e-9", "--steps", "8", "--sample", "4"],
                       "fp128": ["--dt", "100", "--epsilon", "1e-9", "--steps", "6", "--sample", "3"],
                       "fp256": ["--dt", "100", "--epsilon", "1e-9", "--steps", "4", "--sample", "2"]}),
    }
    engines = {"loop": [], "program": ["--engine", "program", "--arith", "fma"]}
    rc, o, e, _ = run([real, "--probe"], env_dev)
    dev = [l for l in e.splitlines() if l.startswith("cft hip: ")]
    print("costs: this port's default path on %s, MEASURED, one run a row, unranked within the noise" % art)
    print("  %s, built %s" % (real, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(real)))))
    print("  device: %s" % (dev[0] if dev else "NONE: %s" % e.strip()[-300:]))
    if rc or not dev:
        return 2
    print("  started %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    head = ("problem", "format", "engine", "steps", "rej", "calls", "divsqrt", "programs", "gpu wall s", "us/call",
            "ms/step", "gpu cpu/wall", "sw wall s", "gpu/sw", "record")
    rows = []
    for pname in args.problems.split(","):
        path, per_fmt = problems[pname]
        for fmt in args.formats.split(","):
            for eng in args.engines.split(","):
                cmd = [real, "--format", fmt, "--problem", path] + per_fmt[fmt] + engines[eng] + ["--quiet"]
                rd, od, ed, _ = run(cmd, env_dev)
                rs, os_, es, _ = run(cmd, env_sw)
                td, ts = trailer(od), trailer(os_)
                lines = [l for l in ed.splitlines() if l.startswith("cft hip: ")]
                same = rd == 0 and rs == 0 and norm(od) == norm(os_) and lines == dev
                try:
                    wall, calls = float(td["wall_seconds"]), int(td["calls"])
                    steps, cpu = int(td["steps_done"]), float(td["seconds"])
                    sw = float(ts["wall_seconds"])
                    row = (pname, fmt, eng, str(steps), td["steps_rejected"], str(calls), td["divsqrt_calls"],
                           td.get("program_calls", "-"), "%.2f" % wall, "%.0f" % (1e6 * wall / calls),
                           "%.0f" % (1e3 * wall / steps), "%.2f" % (cpu / wall), "%.2f" % sw, "%.1f" % (wall / sw),
                           "identical" if same else "DIFFERS")
                except (KeyError, ValueError, ZeroDivisionError):
                    row = (pname, fmt, eng, "-", "-", "-", "-", "-", "-", "-", "-", "-", "-", "-",
                           "rc %d/%d: %s" % (rd, rs, (ed or es).strip()[-120:]))
                rows.append(row)
                print("  " + " | ".join(row))
                sys.stdout.flush()
    print("  ended %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
    print()
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for r in rows:
        print("| " + " | ".join(r) + " |")
    bad = [r for r in rows if r[-1] != "identical"]
    print("RESULT: %s" % ("every GPU row's record identical to its software row's" if not bad else
                          "%d row(s) NOT identical or not run" % len(bad)))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
