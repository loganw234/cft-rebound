#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Round 3's seam tests. Each holds two parcels together, where neither
# parcel's own gate looks (ROUND3.md, "What the integrator keeps").
#
# 1. A ROUTE ACROSS A CHECKPOINT (P2 x P3, P1 x P3). A route never changes
#    a record line (the door's claim, held by each routed parcel's gate).
#    A checkpoint resumes bit for bit (P3's claim, held by
#    tools/check_checkpoint.py). Each gate holds its claim with the
#    other's variable at its default. This part crosses them, with three
#    processes as check_checkpoint's sequences use:
#
#      gate_real --fpNNN --straight              -> the reference dump
#      gate_real --fpNNN --save   ckpt.bin       (the route forced one way)
#      gate_real --fpNNN --resume ckpt.bin       (the route forced the other)
#
#    The resumed dump must still be the straight run's, in both
#    directions. Two defects show here and nowhere else: a route's state
#    that an archive does not carry, and an archive whose restore depends
#    on the route that wrote it.
#
#      P2 x P3   reduce_seg, the corrector's device exit test.
#      P1 x P3   indexed and pair_gather, gravity's routes. The regrow
#                sequence removes a body before the save and adds one
#                back after the load, so gravity's tables are rebuilt on
#                both sides of the checkpoint.
#      P5 x P3   scalar, the one-element constants.
#
#    ROUND3.md named the P2 x P3 seam "a replay across a checkpoint". The
#    two do not meet. The replay is the standalone program's (main()'s
#    parser), and the checkpoint is the drop-in's Simulationarchive. The
#    route is what they share, so the route is the seam.
#
#    The control: the same resume with CFT_REBOUND_NO_ADOPT, the shim
#    discarding what it loaded, must differ from the straight run. If it
#    does not, the comparison proves nothing. The sequences, the process
#    runner and the comparison are check_checkpoint's own, imported, never
#    copied: a new sequence there reaches this gate unedited.
#
# 2. A REPLAYED SCHEDULE OVER THE OPT-IN ROUTES (P1 x P2, P5 x P2). The
#    program records a schedule (--record-iters) with every opt-in route
#    forced one way, and replays it (--replay-iters) with all forced the
#    other.
#    The replay must give the recording's record lines and every counter
#    a route may not move, system by system, and it must have computed no
#    exit test - reduce_seg carried nothing - so the schedule drove it.
#
# 3. THE STEP-CONTROL PROGRAM OVER THE OPT-IN ROUTES (P1 x P4, P5 x P4).
#    The program engine, P4's step control on the device, with every
#    opt-in route forced on and then off, against the loop engine with
#    them off: the loop engine's record and counters, and the
#    step-control program run in both, nextstep_calls above zero and the
#    same in each.
#
# THE PHYSICS COUNTER, in every part. Each forced route must be the one
# taken, by=forced, and must have carried its run: calls and no fallback
# when on, fallbacks and no call when off. A crossing whose route carried
# nothing compares two runs of today's route, and passes for the wrong
# reason. gate_real says what its routes carried when given --routes; the
# program says it in its `# route` lines. The cases that are there for a
# rejection assert that they had one.
#
# No door variable is inherited, so the environment cannot choose a route
# behind the gate's back. CFT_REBOUND_ARTIFACT stays, and on a card the
# runs go one at a time.
#
#   python tools/check_seams.py [--build build] [--quick]
#   (--quick: binary64 only, the light tiers)

import argparse
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
sys.path.insert(0, here)
import check_checkpoint as CC  # noqa: E402

# The routes a merged parcel carries. Add a name here when its parcel
# merges.
ROUTES = ["reduce_seg", "indexed", "pair_gather", "scalar"]
# What parts 2 and 3 force: the routes auto leaves off.
OPT_IN = ("indexed", "pair_gather", "scalar")
SEQUENCES = CC.SEQUENCES
# What a route may not move, and a replay must keep, run and system alike.
SYSTEM = ("steps_done", "iterations_max_exceeded", "steps_rejected", "mean_pc_iterations", "max_pc_iterations")
JOBS = 1 if os.environ.get("CFT_REBOUND_ARTIFACT") else 4

# What gate_real --routes prints (tests/gate_real.c, print_routes()), and
# the program's trailer lines.
GATE_ROUTE = re.compile(r"^gate_real: route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)\s*$", re.M)
PROGRAM_ROUTE = re.compile(r"^# route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)\s*$", re.M)


def door_var(route):
    return "CFT_REBOUND_ROUTE_" + route.upper()


def routes_in(text, pattern):
    return {m.group(1): (m.group(2), m.group(3), int(m.group(4)), int(m.group(5))) for m in pattern.finditer(text)}


def carried(routes, route, forced):
    """Whether a run took the route it was forced to and that route
    carried the run, with a phrase that says what it carried."""
    if route not in routes:
        return False, "no route line names %s" % route
    state, by, calls, falls = routes[route]
    if state != forced or by != "forced":
        return False, "%s=%s by=%s, though it was forced %s" % (route, state, by, forced)
    if forced == "on" and (calls == 0 or falls != 0):
        return False, "%s forced on carried %d calls and %d fallbacks" % (route, calls, falls)
    if forced == "off" and (falls == 0 or calls != 0):
        return False, "%s forced off carried %d calls and %d fallbacks" % (route, calls, falls)
    return True, ("%d calls on the new route" % calls) if forced == "on" else ("%d fallbacks to today's" % falls)


def opt_in_carried(routes, forced):
    oks, said = zip(*(carried(routes, g, forced) for g in OPT_IN))
    return all(oks), ", ".join("%s %s" % (g, s) if ok else s for g, ok, s in zip(OPT_IN, oks, said))


def kv(text):
    return dict(x.partition("=")[::2] for x in text.split())


def program(cmd, forced=None):
    """One run of the standalone program, the opt-in routes forced as
    asked: its rc, stderr, record lines, trailer, each system's trailer
    and its route lines."""
    env = dict(os.environ)
    env.update({door_var(g): forced for g in OPT_IN} if forced else {})
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env, cwd=root)
    tail, systems = {}, {}
    for line in p.stdout.splitlines():
        if line.startswith("# steps_done="):
            tail = kv(line[2:])
        m = re.match(r"^# system (\d+) (steps_done=.*)$", line)
        if m:
            systems[int(m.group(1))] = kv(m.group(2))
    return {"rc": p.returncode, "err": p.stderr, "recs": [l for l in p.stdout.splitlines() if l and not l.startswith("#")],
            "tail": tail, "systems": systems, "routes": routes_in(p.stdout, PROGRAM_ROUTE)}


def same_counters(a, b):
    """The run's counters a route may not move, and each system's."""
    return (all(a["tail"].get(k) is not None and a["tail"].get(k) == b["tail"].get(k) for k in SYSTEM)
            and a["tail"].get("pc_lane_efficiency") == b["tail"].get("pc_lane_efficiency")
            and set(a["systems"]) == set(b["systems"])
            and all(a["systems"][s].get(k) == b["systems"][s].get(k) for s in a["systems"] for k in SYSTEM))


def across_checkpoint(exe, formats, wd, pool, verdict):
    """Part 1. Returns how many crossings had their route carry both
    halves, and how many there were to cross."""
    keys = [(fmt, seq, extra) for fmt in formats for seq, extra, _ in SEQUENCES]
    straight = dict(zip([k[:2] for k in keys],
                        pool.map(lambda k: CC.run(exe, ["--" + k[0]] + k[2] + ["--straight"], True), keys)))
    jobs = [(fmt, seq, extra, route, wrote, read) for fmt, seq, extra in keys for route in ROUTES
            for wrote, read in (("on", "off"), ("off", "on"))]

    def crossing(job):
        fmt, seq, extra, route, wrote, read = job
        flag, var = "--" + fmt, door_var(route)
        head = "%s %s: %s %s at the save" % (fmt, seq, route, wrote)
        ok, ref, err = straight[(fmt, seq)]
        if not ok:
            return False, False, "%s: the straight run failed: %s" % (head, err.strip()[-200:])
        ckpt = os.path.join(wd, "seam_%s_%s_%s_%s.bin" % (fmt, seq, route, wrote))
        ok, _, err = CC.run(exe, [flag] + extra + ["--routes", "--save", ckpt], False, {var: wrote})
        if not ok or not os.path.exists(ckpt):
            return False, False, "%s: the save failed: %s" % (head, err.strip()[-200:])
        okw, at_save = carried(routes_in(err, GATE_ROUTE), route, wrote)
        ok, resumed, err = CC.run(exe, [flag] + extra + ["--routes", "--resume", ckpt], True, {var: read})
        if not ok:
            return False, False, "%s: the resume with %s=%s failed: %s" % (head, var, read, err.strip()[-200:])
        okr, at_resume = carried(routes_in(err, GATE_ROUTE), route, read)
        okc, control, _ = CC.run(exe, [flag] + extra + ["--resume", ckpt], True,
                                 {var: read, "CFT_REBOUND_NO_ADOPT": "1"})
        sensitive = okc and control != ref
        same = resumed == ref
        return (same and sensitive and okw and okr, okw and okr,
                "%s (%s), %s at the resume (%s): %s; %s"
                % (head, at_save, read, at_resume,
                   "%d values identical to the straight run" % len(ref.splitlines())
                   if same else "the resume DIFFERS - %s" % CC.first_difference(ref, resumed),
                   "CFT_REBOUND_NO_ADOPT differs, so the comparison bites" if sensitive
                   else "NOT SENSITIVE: a resume that discards the archive still matches"))

    print("check_seams 1: a route across a checkpoint, routes %s, formats %s, sequences %s"
          % (", ".join(ROUTES), ", ".join(formats), ", ".join(s[0] for s in SEQUENCES)))
    n = 0
    for ok, both, text in pool.map(crossing, jobs):
        verdict(ok, text)
        n += both
    return n, len(jobs)


def replay_over_routes(cft, formats, problems, wd, pool, verdict):
    """Part 2: a schedule recorded on one side of the opt-in routes,
    replayed on the other."""
    jobs = [(fmt, i, rec, rep) for fmt in formats for i, p in enumerate(problems) if not p[3] or fmt in p[3]
            for rec, rep in (("on", "off"), ("off", "on"))]

    def one(job):
        fmt, i, rec, rep = job
        name, prob, args, _, rejects = problems[i]
        sched = os.path.join(wd, "replay_%s_%d_%s.sched" % (fmt, i, rec))
        base = [cft, "--format", fmt, "--problem", prob, "--quiet"] + args
        a = program(base + ["--record-iters", sched], rec)
        if a["rc"]:
            return False, 0, "%s %s: the recording did not run (rc %d): %s" % (fmt, name, a["rc"], a["err"].strip()[-200:])
        b = program(base + ["--replay-iters", sched], rep)
        if b["rc"]:
            return False, 0, "%s %s: the replay did not run (rc %d): %s" % (fmt, name, b["rc"], b["err"].strip()[-200:])
        oka, at_rec = opt_in_carried(a["routes"], rec)
        okb, at_rep = opt_in_carried(b["routes"], rep)
        test = b["routes"].get("reduce_seg")
        no_test = test is not None and test[2:] == (0, 0)
        same = b["recs"] == a["recs"] and len(a["recs"]) > 0 and same_counters(a, b)
        rej = a["tail"].get("steps_rejected", "?")
        physics = not rejects or rej not in ("0", "?")
        return (same and no_test and oka and okb and physics, oka and okb,
                "%s %s: recorded with the opt-in routes %s (%s), replayed %s (%s): %s; %s; steps_rejected=%s%s"
                % (fmt, name, rec, at_rec, rep, at_rep,
                   "%d record lines and every counter a route may not move identical, system by system"
                   % len(a["recs"]) if same else "the replay DIFFERS from its recording",
                   "it computed no exit test, so the schedule drove it" if no_test
                   else "IT COMPUTED AN EXIT TEST (reduce_seg %s), so the schedule did not drive it" % (test,),
                   rej, "" if not rejects else (", the rejection the case needs" if physics else ", AND THE CASE NEEDS ONE")))

    print("check_seams 2: a replayed schedule over the opt-in routes, formats %s" % ", ".join(formats))
    n = 0
    for ok, both, text in pool.map(one, jobs):
        verdict(ok, text)
        n += both
    return n, len(jobs)


def program_over_routes(cft, progs, formats, problems, pool, verdict):
    """Part 3: the step-control program over the opt-in routes, on and
    off, against the loop engine with them off."""
    jobs = [(fmt, i) for fmt in formats for i, p in enumerate(problems) if not p[3] or fmt in p[3]]

    def one(job):
        fmt, i = job
        name, prob, args, _, rejects = problems[i]
        common = [cft, "--format", fmt, "--problem", prob, "--quiet", "--arith", "fma", "--max-iter", "60"] + args
        loop = program(common + ["--engine", "loop"], "off")
        runs = {side: program(common + ["--engine", "program", "--programs", progs], side) for side in ("on", "off")}
        bad = ["%s rc %d: %s" % (k, r["rc"], r["err"].strip()[-160:]) for k, r in (("loop", loop),) + tuple(runs.items())
               if r["rc"]]
        if bad:
            return False, 0, "%s %s: %s" % (fmt, name, "; ".join(bad))
        okl, at_loop = opt_in_carried(loop["routes"], "off")
        said, oks, steps = [], [], []
        for side, r in runs.items():
            ok, at = opt_in_carried(r["routes"], side)
            oks.append(ok)
            said.append("%s (%s)" % (side, at))
            steps.append(r["tail"].get("nextstep_calls", "0"))
        same = all(r["recs"] == loop["recs"] and len(loop["recs"]) > 0 and same_counters(r, loop) for r in runs.values())
        stepped = steps[0] == steps[1] and int(steps[0]) > 0
        rej = loop["tail"].get("steps_rejected", "?")
        physics = not rejects or rej not in ("0", "?")
        return (same and stepped and okl and all(oks) and physics, okl and all(oks),
                "%s %s: the program engine with the opt-in routes %s, against the loop engine with them off (%s): %s; "
                "nextstep_calls %s; steps_rejected=%s%s"
                % (fmt, name, " and ".join(said), at_loop,
                   "%d record lines and every counter a route may not move identical" % len(loop["recs"])
                   if same else "the program engine DIFFERS from the loop engine",
                   " and ".join(steps) + ("" if stepped else ", SO THE STEP CONTROL PROGRAM DID NOT RUN IN BOTH"),
                   rej, "" if not rejects else (", the rejection the case needs" if physics else ", AND THE CASE NEEDS ONE")))

    print("check_seams 3: the step-control program over the opt-in routes, formats %s" % ", ".join(formats))
    n = 0
    for ok, both, text in pool.map(one, jobs):
        verdict(ok, text)
        n += both
    return n, len(jobs)


def main():
    ap = argparse.ArgumentParser(description="round 3's seam tests")
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true", help="binary64 only")
    a = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    gate = os.path.join(a.build, "gate_real")
    cft = os.path.join(a.build, "ias15_cft" + exe)
    progs = os.path.join(root, "programs", "out")
    formats = ["fp64"] if a.quick else CC.FORMATS
    for k in sorted(os.environ):
        if k.upper().startswith("CFT_REBOUND_ROUTE_") or k.upper().startswith("CFT_REBOUND_HIDE"):
            print("check_seams: %s=%s is not inherited" % (k, os.environ.pop(k)))
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= ok
        print("  %s %s" % ("ok  " if ok else "FAIL", text))

    with tempfile.TemporaryDirectory(prefix="cft_seams_") as wd, ThreadPoolExecutor(max_workers=JOBS) as pool:
        ens = os.path.join(wd, "outer4.txt")
        rc = subprocess.run([sys.executable, os.path.join(here, "make_ensemble.py"),
                             os.path.join(root, "data", "problems", "outer.txt"), "--members", "4", "--out", ens,
                             "--geometric", "2", "x", "0.25"], stdout=subprocess.DEVNULL, cwd=root).returncode
        pdir = os.path.join(root, "data", "problems")
        # name, problem, arguments, formats (None: every one), whether the
        # case is there for its rejection
        problems = [
            ("pythagorean, three bodies", os.path.join(pdir, "pythagorean.txt"),
             ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "40", "--sample", "20"], None, False),
            ("kepler e=0.99, a rejected step", os.path.join(pdir, "kepler_e099.txt"),
             ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "30", "--sample", "10"], None, True),
            ("outer x4, an ensemble", ens,
             ["--dt", "40", "--epsilon", "1e-9", "--steps", "10", "--sample", "10"], None, False),
        ]
        if rc:
            verdict(False, "make_ensemble.py failed (rc %d), so the ensemble cases cannot run" % rc)
        counts = [across_checkpoint(gate, formats, wd, pool, verdict),
                  replay_over_routes(cft, formats, problems, wd, pool, verdict),
                  program_over_routes(cft, progs, formats, problems, pool, verdict)]
    # A gate that crossed nothing passes nothing: every crossing's routes
    # must have carried both of its halves.
    for (n, want), what in zip(counts, ("route crossings of a checkpoint", "replays across the opt-in routes",
                                        "program-engine pairs over the opt-in routes")):
        verdict(n == want and want > 0, "%d of %d %s carried by the forced routes on both sides" % (n, want, what))
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
