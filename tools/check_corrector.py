#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of the corrector's exit test (ROUND3.md, P2): the test on the
# device, the recorded schedule, and the self-test that holds pc_error().
#
#   python tools/check_corrector.py [--build build] [--exe ias15_cft]
#                                   [--formats fp64,fp128,fp256] [--quick]
#                                   [--keep DIR]
#   (`make check-quick` and `make check-light` pass --quick: binary64 only,
#    and shorter runs)
#
# What it asserts, and what lets each assertion fail:
#   1. --corrector-selftest at every format. pc_error() on both routes
#      equals a reference computed outside it, in every mode, on systems
#      crafted for each property of the device route (the +0 identity, the
#      normal-only mask, the 1/1 a system out of the corrector is given)
#      and on random ones; the data holds every class of value; calls= and
#      divsqrt_calls= move the same on both routes; both controls see their
#      difference; and the image the engine loaded is the pinned tree's, by
#      the SHA-256 that tree's own programs/MANIFEST gives it. On the
#      software backend the route must be published: nothing passes
#      uncompared.
#   2. Both routes, run explicitly, on every problem in data/problems/ in
#      every mode at binary64, on a set at the wider formats, and on an
#      ensemble and the program engine. The route forced on and forced off
#      must give the same record lines, the same counters (the eight
#      check_door.py holds, program_calls, and each system's), the same passes in
#      every step attempt: their schedules, byte for byte. Each run's
#      `# route` line is asserted: on by=forced with no fallback, off
#      by=forced with no call, and the maxima one carried as calls the
#      other carried as fallbacks - a count DERIVED from the schedule, its
#      passes times the maxima a pass the self-test measured.
#   3. The replay. --replay-iters of each of those schedules gives its run's
#      record lines and counters, and computes no test: its route line
#      carries nothing, and divsqrt_calls= and calls= fall by exactly what
#      the test costs, from the schedule and the self-test's counts. The
#      physics that makes this a test is asserted, not assumed: at every
#      format a step rejected, a cap hit, and an ensemble whose systems
#      leave the corrector at different passes in one attempt.
#   4. The replay drives the run, at every format. A count moved by one
#      changes no record line on most attempts (P2's ledger entries), so
#      the controls are: the first attempt cut to one pass must change the
#      record; a moved count's passes must show in the replay's own
#      counters; a capped count moved by one is refused by name; and in a
#      mixed attempt of the ensemble, ONE system that is not the longest cut
#      to one pass must change its own lines and no other system's, and
#      those lines must be the same cut of that system run alone (--member),
#      whose schedule must be the ensemble's column for it attempt by
#      attempt. That holds each system's own exit exactly: a replay that
#      kept every system in the corrector to the attempt's last count (V2's
#      M2), or froze one a pass late (M2b), passed everything else. Each of
#      the three cut replays must RUN. On these problems every one does, at
#      every format, so a refusal by name says the cut desynced the replay,
#      and fails: when it passed, a replay that refused every count of 1 by
#      name passed all three (V2's Rboth).
#   5. Refusals by name: a schedule from another run, one cut short, one
#      too long, a malformed field, a file that is not a schedule, and
#      --trace-pc or --record-iters beside --replay-iters.
#   6. The run key holds every field. A schedule replayed under a run that
#      differs in one thing - epsilon, each of the four things the state's
#      hash covers (a position, a velocity, a mass, a dt-file step), and
#      every other field the key carries, read from the key - is refused by
#      name, and a field that no case moves fails by name (V2's M3, M3s,
#      M3v, M3m and M3dt).
#
# Each run works in <build>/corrector/<pid>, so two runs sharing a --build
# directory do not clobber each other; a passing run removes it, a failing
# one keeps it.
#
# The route's name and variable are read from the program, never typed.

import argparse
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)

# check_door.py's eight, and program_calls: the card-day normalisers keep
# it, so a route must not move it either (the normalabs run is not one)
COUNTERS = ("steps_done", "iterations_max_exceeded", "steps_rejected", "mean_pc_iterations",
            "max_pc_iterations", "flags_seen", "calls", "divsqrt_calls", "program_calls")
SYSTEM = ("steps_done", "iterations_max_exceeded", "steps_rejected", "mean_pc_iterations", "max_pc_iterations")
ok_all = True


def verdict(ok, text):
    global ok_all
    ok_all &= bool(ok)
    print("  %s %s" % ("ok  " if ok else "FAIL", text))
    return ok


def run(cmd, extra_env=None):
    """No door variable is inherited, so the environment cannot choose a
    route behind the gate's back. CFT_REBOUND_ARTIFACT stays."""
    env = {k: v for k, v in os.environ.items()
           if not k.upper().startswith("CFT_REBOUND_ROUTE_") and not k.upper().startswith("CFT_REBOUND_HIDE")}
    env.update(extra_env or {})
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    return p.returncode, p.stdout, p.stderr


def kv(text):
    return dict(x.partition("=")[::2] for x in text.split())


def parse(out, route):
    """Record lines, the trailer, each system's trailer, and the route's line."""
    rec = [l for l in out.splitlines() if not l.startswith("#")]
    tail, systems, rline = {}, {}, None
    for l in out.splitlines():
        if l.startswith("# steps_done="):
            tail = kv(l[2:])
        m = re.match(r"^# system (\d+) (steps_done=.*)$", l)
        if m:
            systems[int(m.group(1))] = kv(m.group(2))
        m = re.match(r"^# route %s=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)$" % re.escape(route), l)
        if m:
            rline = (m.group(1), m.group(2), int(m.group(3)), int(m.group(4)))
    return rec, tail, systems, rline


def schedule(path):
    with open(path, encoding="utf-8") as f:
        lines = f.read().splitlines()
    return [l for l in lines if l.startswith("#")], [l.split() for l in lines if l and not l.startswith("#")]


def tally(body):
    """Attempts, the shared passes (one pc_error() a pass), the capped
    fields, and the attempts whose systems left at different passes."""
    passes = capped = mixed = 0
    for fields in body:
        counts = [int(f.rstrip("c")) for f in fields if f != "-"]
        passes += max(counts) if counts else 0
        capped += sum(f.endswith("c") for f in fields)
        mixed += len(set(counts)) > 1
    return len(body), passes, capped, mixed


def write(path, head, body):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(head + [" ".join(b) for b in body]) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--exe", default="ias15_cft")
    ap.add_argument("--formats", default="fp64,fp128,fp256")
    ap.add_argument("--quick", action="store_true", help="binary64 only, shorter runs")
    ap.add_argument("--library", action="store_true", help="also run check_dropin --light with the route forced on")
    ap.add_argument("--keep", default=None, help="directory for schedules and ensembles, kept "
                    "(default <build>/corrector/<pid>, removed after a pass)")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, args.exe + exe)
    formats = ["fp64"] if args.quick else args.formats.split(",")
    # A directory of this process's own, so two runs sharing one --build
    # cannot clobber each other's schedules (V2, 2026-11-02: one read the
    # other's fp256 schedules as fp64 ones). Removed after a pass, kept after
    # a failure for what it holds.
    work = args.keep or os.path.join(args.build, "corrector", str(os.getpid()))
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)
    software = not os.environ.get("CFT_REBOUND_ARTIFACT")
    progs = os.path.join(root, "programs", "out")
    manifest = os.path.join(root, "third_party", "cft-fp256", "programs", "MANIFEST")
    pinned = {}
    with open(manifest, encoding="utf-8") as f:
        for line in f:
            p = line.split()
            if len(p) == 2 and not line.startswith("#"):
                pinned[p[1]] = p[0]

    # 1. The self-test, which also names the route and measures what a
    #    pass of pc_error() costs in each mode.
    route = env = None
    cost = {}                        # format -> mode -> (calls, divides) a pass
    published = True
    for fmt in formats:
        print("check_corrector at %s: pc_error() on both routes against a reference computed outside it" % fmt)
        rc, out, err = run([cft, "--corrector-selftest", "--format", fmt])
        m = re.search(r"^corrector route (\w+)=(on|off) by=(\w+) env=(\S+)$", out, re.M)
        s = re.search(r"^corrector: (\d+) identical, (\d+) differ, (\d+) not compared, (\d+) flag differences, "
                      r"(\d+) accounting errors; controls (\d+) of (\d+)$", out, re.M)
        if not m or not s:
            verdict(False, "--corrector-selftest did not run (rc %d): %s" % (rc, (err or out)[-400:].strip()))
            continue
        route, env, by = m.group(1), m.group(4), m.group(3)
        same, differ, skipped, fbad, cbad, cok, cn = (int(g) for g in s.groups())
        for line in out.splitlines():
            if line.startswith(("DIFF ", "COUNT ", "CONTROL FAILED", "skip ")):
                print("       " + line)
        verdict(rc == 0 and differ == 0 and fbad == 0 and same > 0,
                "%d of %d mode-and-route results equal the reference in every system, bits and flags (rc %d)"
                % (same, same + differ + skipped, rc))
        verdict(cbad == 0, "each route's accounting: calls= and divsqrt_calls= alike, one reduce_seg operation a "
                "maximum on the counter meant (%d errors)" % cbad)
        verdict(cn > 0 and cok == cn, "controls: %d of %d saw their difference" % (cok, cn))
        d = re.search(r"^corrector data: normal (\d+), \+0 (\d+), -0 (\d+), subnormal (\d+), infinite (\d+), NaN (\d+); "
                      r"(\d+) of (\d+) systems out of the corrector$", out, re.M)
        counts = [int(g) for g in d.groups()] if d else [0] * 8
        verdict(d and all(c > 0 for c in counts[:7]),
                "the data holds every class, as libcft classifies it: normal %d, +0 %d, -0 %d, subnormal %d, "
                "infinite %d, NaN %d, and %d of %d systems out of the corrector" % tuple(counts))
        if by == "absent" and not software:
            published = False
            print("  note %s: the device does not publish it, so only today's route was compared" % route)
        else:
            verdict(skipped == 0 and (by == "opt_in" or not software),
                    "%s by=%s, and nothing left uncompared%s"
                    % (route, by, " (the software backend publishes every bit, and auto leaves this route off)"
                       if software else ""))
            g = re.search(r"^corrector image (\S+) sha256=([0-9a-f]{64})$", out, re.M)
            name = (g.group(1) + ".cftp") if g else None
            verdict(g and name in pinned and pinned[name] == g.group(2),
                    "the image the engine loaded, %s, is the pinned tree's: sha256 %s, and its MANIFEST says %s"
                    % (name, g.group(2)[:16] if g else "?", pinned.get(name, "nothing")[:16]))
        cost[fmt] = {}
        for c in re.finditer(r"^ok   mode (\d): calls= moved (\d+) and divsqrt_calls= (\d+) on each route$", out, re.M):
            cost[fmt][int(c.group(1))] = (int(c.group(2)), int(c.group(3)))
        if published:
            verdict(sorted(cost[fmt]) == [0, 1, 2, 3], "a pass of pc_error() costs, by mode: %s"
                    % ", ".join("%d: %d calls and %d divide" % (k, v[0], v[1]) for k, v in sorted(cost[fmt].items())))
    if route is None:
        print("RESULT: FAIL")
        return 1

    # The problems: every file in data/problems/, at every mode; and an
    # ensemble whose members reject and accept in the same attempt.
    ens = os.path.join(work, "kepler_spread_E5.txt")
    rc, out, err = run([sys.executable, os.path.join(here, "make_ensemble.py"), os.path.join(root, "data", "problems", "kepler.txt"),
                        "--members", "5", "--out", ens, "--offset", "planet", "x", "0.4"])
    if rc:
        verdict(False, "make_ensemble.py: %s" % (err or out)[-300:])
    pdir = os.path.join(root, "data", "problems")
    probs = sorted(f for f in os.listdir(pdir) if f.endswith(".txt"))
    long_dt = {"outer.txt": "40"}                   # the outer planets' natural step, in days
    steps = "30" if args.quick else "60"

    def config(fmt, label, prob, extra, mode, engine="loop"):
        a = ["--format", fmt, "--problem", prob, "--quiet", "--adaptive-mode", str(mode)] + extra
        if mode in (0, 1):
            a.append("--no-flag-abort")             # REBOUND's own estimate divides by zero there
        if engine == "program":
            a += ["--arith", "fma", "--engine", "program", "--programs", progs]
        return {"fmt": fmt, "label": "%s, mode %d%s" % (label, mode, ", program engine" if engine == "program" else ""),
                "args": a, "mode": mode}

    configs = []
    for fmt in formats:
        wide = fmt != "fp64"
        n = ("12" if fmt == "fp256" else "20") if wide else steps
        cap = ["--max-iter", "60"] if fmt == "fp256" else []    # so its counts vary; one run keeps the cap at 12
        run_args = lambda dt, eps: ["--dt", dt, "--epsilon", eps, "--steps", n, "--sample", str(int(n) // 3)]
        for f in probs:
            p = os.path.join(pdir, f)
            dt = long_dt.get(f)
            if wide and f not in ("kepler.txt", "kepler_e099.txt", "pythagorean.txt"):
                continue
            modes_a = (0, 1, 2, 3) if not wide else ((2, 0) if f == "kepler_e099.txt" else (2,))
            for mode in modes_a:
                configs.append(config(fmt, "%s adaptive" % f[:-4], p, run_args(dt or "0.01", "1e-9") + cap, mode))
            for mode in ((0, 2) if not wide else ()):
                configs.append(config(fmt, "%s fixed" % f[:-4], p, run_args(dt or "0.05", "0"), mode))
        kepler = os.path.join(pdir, "kepler.txt")
        configs.append(config(fmt, "kepler from dt0=1.0 (a rejection)", kepler, run_args("1.0", "1e-9") + cap, 2))
        if fmt == "fp256":
            configs.append(config(fmt, "kepler adaptive at --max-iter 12 (every attempt capped)", kepler,
                                  run_args("0.01", "1e-9"), 2))
        configs.append(config(fmt, "kepler spread E=5 from dt0=0.6 (mixed rejection)", ens,
                              ["--dt", "0.6", "--epsilon", "1e-9", "--steps", n if wide else str(2 * int(steps)),
                               "--sample", "10"] + (cap or ["--max-iter", "60"]), 2))
        if not wide:
            configs.append(config(fmt, "kepler spread E=5 from dt0=0.6 (mixed rejection)", ens,
                                  ["--dt", "0.6", "--epsilon", "1e-9", "--steps", str(2 * int(steps)), "--sample", "10",
                                   "--max-iter", "60"], 0))
        configs.append(config(fmt, "kepler adaptive", kepler, run_args("0.01", "1e-9") + cap, 2, "program"))
        if not wide:
            configs.append(config(fmt, "outer fixed", os.path.join(pdir, "outer.txt"), run_args("40", "0"), 2, "program"))

    # 2 and 3. Each configuration: the route on, the route off, each with its
    # schedule recorded, then a replay of the first.
    def job(i):
        c = configs[i]
        s_on, s_off = os.path.join(work, "c%d-on.sched" % i), os.path.join(work, "c%d-off.sched" % i)
        cmd = [cft] + c["args"]
        on = run(cmd + ["--record-iters", s_on], {env: "on"}) if published else None
        off = run(cmd + ["--record-iters", s_off], {env: "off"})
        rep = run(cmd + ["--replay-iters", s_on if published else s_off])
        return on, off, rep, s_on, s_off

    # One at a time on a card: four would share the one tile.
    with ThreadPoolExecutor(max_workers=1 if os.environ.get("CFT_REBOUND_ARTIFACT") else 4) as pool:
        results = list(pool.map(job, range(len(configs))))

    seen = {fmt: {"rejected": 0, "capped": 0, "mixed": 0, "mixed_reject": 0} for fmt in formats}
    for c, (on, off, rep, s_on, s_off) in zip(configs, results):
        fmt, mode = c["fmt"], c["mode"]
        print("check_corrector at %s: %s" % (fmt, c["label"]))
        if (on and on[0]) or off[0]:
            verdict(False, "did not run (rc %s/%d): %s" % (on[0] if on else "-", off[0], ((on[2] if on else "") + off[2])[-300:].strip()))
            continue
        r1, t1, y1, l1 = parse(off[1], route)
        head, body = schedule(s_off)
        attempts, passes, capped, mixed = tally(body)
        per = 1 if mode == 0 else 2                             # maxima a pass: INDIVIDUAL's one, or two
        if published:
            r0, t0, y0, l0 = parse(on[1], route)
            same_sched = open(s_on, encoding="utf-8").read() == open(s_off, encoding="utf-8").read()
            missing = [k for k in COUNTERS if k not in t0 or k not in t1]
            verdict(r0 == r1 and len(r0) > 0 and not missing and all(t0[k] == t1[k] for k in COUNTERS)
                    and all(y0.get(s, {}).get(k) == y1.get(s, {}).get(k) for s in set(y0) | set(y1) for k in SYSTEM)
                    and same_sched,
                    "device and today's route: %d record lines, the counters%s and every attempt's passes identical "
                    "(%d attempts, %d passes)" % (len(r0), " (MISSING %s)" % missing if missing else "", attempts, passes))
            verdict(l0 == ("on", "forced", per * passes, 0) and l1 == ("off", "forced", 0, per * passes),
                    "route lines: %s on %s, off %s; wanted %d maxima (%d passes x %d) as calls, then as fallbacks"
                    % (route, l0, l1, per * passes, passes, per))
        else:
            verdict(l1 == ("off", "forced", 0, per * passes), "route line off: %s; wanted %d fallbacks" % (l1, per * passes))
        # the replay
        r2, t2, y2, l2 = parse(rep[1], route)
        calls_pass, div_pass = cost.get(fmt, {}).get(mode, (None, None))
        if rep[0]:
            verdict(False, "the replay did not run (rc %d): %s" % (rep[0], rep[2][-300:].strip()))
            continue
        if published:                       # held to the run whose schedule it replays
            r1, t1, y1 = r0, t0, y0
            attempts, passes, capped, mixed = tally(schedule(s_on)[1])
        same = (r2 == r1 and all(t2.get(k) == t1.get(k) for k in SYSTEM)
                and all(y2.get(s, {}).get(k) == y1.get(s, {}).get(k) for s in set(y1) | set(y2) for k in SYSTEM)
                and t2.get("pc_lane_efficiency") == t1.get("pc_lane_efficiency"))
        verdict(same, "the replay: %d record lines and every counter of the run, iterations_max_exceeded=%s "
                "steps_rejected=%s mean_pc_iterations=%s" % (len(r2), t2.get("iterations_max_exceeded"),
                                                             t2.get("steps_rejected"), t2.get("mean_pc_iterations")))
        dcalls = int(t1["calls"]) - int(t2["calls"])
        ddiv = int(t1["divsqrt_calls"]) - int(t2["divsqrt_calls"])
        want_calls = None if calls_pass is None else calls_pass * passes + 2 * (passes + attempts)
        want_div = None if div_pass is None else div_pass * passes
        flags_sub = (int(t1["flags_seen"], 16) | int(t2["flags_seen"], 16)) == int(t1["flags_seen"], 16)
        verdict(l2 is not None and l2[2:] == (0, 0) and ddiv == want_div and dcalls == want_calls and flags_sub,
                "it computed no test: %s carried %s; divsqrt_calls %d fewer (%d passes x %s); calls %d fewer "
                "(wanted %s: %s a pass, and two compares at each of %d tests); flags 0x%s within 0x%s"
                % (route, "nothing" if l2 and l2[2:] == (0, 0) else str(l2), ddiv, passes, div_pass, dcalls, want_calls,
                   calls_pass, passes + attempts, t2["flags_seen"][2:], t1["flags_seen"][2:]))
        seen[fmt]["rejected"] += int(t1["steps_rejected"]) > 0
        seen[fmt]["capped"] += capped > 0 and int(t1["iterations_max_exceeded"]) == capped
        if y1:
            rej = [int(y1[s]["steps_rejected"]) for s in sorted(y1)]
            seen[fmt]["mixed"] += mixed > 0
            seen[fmt]["mixed_reject"] += any(rej) and not all(rej)
            print("       ensemble: %d of %d attempts had systems leave the corrector at different passes; "
                  "steps_rejected by system %s" % (mixed, attempts, rej))
    for fmt in formats:
        s = seen[fmt]
        verdict(s["rejected"] > 0 and s["capped"] > 0 and s["mixed"] > 0 and s["mixed_reject"] > 0,
                "%s: the replays crossed %d runs with a rejected step, %d with a capped attempt (each counted as "
                "iterations_max_exceeded), %d ensembles whose systems left at different passes, %d with some members "
                "rejecting and some not" % (fmt, s["rejected"], s["capped"], s["mixed"], s["mixed_reject"]))

    # 4. The replay drives the run, at every format, on its rejection run.
    def schedule_of(i):
        return results[i][3] if published else results[i][4]

    def desync(rc, err, who="the replay"):
        # A cut replay must RUN: on these problems every one does, at every
        # format. A refusal by name says the cut desynced the replay from its
        # schedule, and it fails - accepted, it let a replay that refused
        # every count of 1 by name pass all three cuts (V2's Rboth).
        tail = err.strip()[-140:]
        if rc == 3 and "--replay-iters" in err:
            return "the cut desynced %s, which refused it by name (rc 3): %s" % (who, tail)
        return "%s failed (rc %d): %s" % (who, rc, tail)

    for fmt in formats:
        k = next(i for i, c in enumerate(configs) if c["fmt"] == fmt and "dt0=1.0" in c["label"])
        c = configs[k]
        base = [cft] + c["args"]
        if results[k][1][0] or (published and results[k][0][0]):
            verdict(False, "%s did not run at %s, so the replay's controls have no schedule to move" % (c["label"], fmt))
            continue
        head, body = schedule(schedule_of(k))
        rec0, t0, _, _ = parse(results[k][1][1], route)
        print("check_corrector at %s: the replay drives the run (%s)" % (fmt, c["label"]))
        first_free = next(i for i, f in enumerate(body) if not f[0].endswith("c"))
        # (a) one count moved by one runs, and the replay's own counters show
        #     exactly the passes it was given, where no record bit moves
        for d in (1, -1):
            mut = [list(f) for f in body]
            mut[first_free][0] = str(int(body[first_free][0]) + d)
            path = os.path.join(work, "%s-moved%+d.sched" % (fmt, d))
            write(path, head, mut)
            rc, out, err = run(base + ["--replay-iters", path])
            rec, t, _, _ = parse(out, route)
            want = "%.3f" % (sum(int(f[0].rstrip("c")) for f in mut) / len(mut))
            verdict(rc == 0 and t.get("mean_pc_iterations") == want and want != t0.get("mean_pc_iterations"),
                    "attempt %d's %s passes moved %+d: the replay ran them, mean_pc_iterations %s (wanted %s; the "
                    "recording's %s), and %d of %d record lines differ"
                    % (first_free, body[first_free][0], d, t.get("mean_pc_iterations"), want, t0.get("mean_pc_iterations"),
                       sum(a != b for a, b in zip(rec, rec0)) + abs(len(rec) - len(rec0)), len(rec0)))
        # (b) the first attempt cut to one pass: the replay runs, and the
        #     record changes
        mut = [list(f) for f in body]
        mut[0][0] = "1"
        path = os.path.join(work, "%s-first-one.sched" % fmt)
        write(path, head, mut)
        rc, out, err = run(base + ["--replay-iters", path])
        rec, _, _, _ = parse(out, route)
        nd = sum(a != b for a, b in zip(rec, rec0)) + abs(len(rec) - len(rec0))
        verdict(rc == 0 and nd > 0,
                "the first attempt's %s passes cut to 1: %s" % (body[0][0], ("%d of %d record lines differ" % (nd, len(rec0)))
                                                                  if rc == 0 else desync(rc, err)))
        # (c) a count the cap stopped, moved by one, is refused by name: on
        #     the first run of this format that has one
        hit = None
        for i, cc in enumerate(configs):
            if cc["fmt"] != fmt or results[i][1][0] or (published and results[i][0][0]):
                continue
            h2, b2 = schedule(schedule_of(i))
            at = next((j for j, f in enumerate(b2) if len(f) == 1 and f[0].endswith("c")), None)
            if at is not None:
                hit = (cc, h2, b2, at)
                break
        if hit is None:
            verdict(False, "no single-system run at %s had a capped attempt, so the cap's refusal is not exercised" % fmt)
            continue
        cc, h2, b2, at = hit
        mut = [list(f) for f in b2]
        mut[at][0] = "%dc" % (int(b2[at][0].rstrip("c")) + 1)
        path = os.path.join(work, "%s-capped.sched" % fmt)
        write(path, h2, mut)
        rc, out, err = run([cft] + cc["args"] + ["--replay-iters", path])
        verdict(rc == 3 and "the cap stopped system 0" in err,
                "%s: attempt %d, which the cap stopped at %s, moved by one, is refused by name (rc %d)"
                % (cc["label"], at, b2[at][0], rc))
        # (d) Each system's own exit. In the mixed ensemble's first attempt
        #     whose systems leave at different passes, the system with the
        #     most passes short of the attempt's longest is cut to ONE. The
        #     replay must run, and change that system's record lines and no
        #     other system's. A replay that kept every system in the
        #     corrector until the attempt's last count, while counting each
        #     one's passes exactly, reproduces every recording - the extra
        #     passes land on a fixed point - and leaves this cut converging
        #     with the rest (V2's M2, 2026-11-02).
        ke = next((i for i, c2 in enumerate(configs) if c2["fmt"] == fmt and "mixed rejection" in c2["label"]
                   and c2["mode"] == 2), None)
        if ke is None or results[ke][1][0] or (published and results[ke][0][0]):
            verdict(False, "the mixed ensemble did not run at %s, so no system's own exit is held" % fmt)
            continue
        he, be = schedule(schedule_of(ke))
        pick = None
        for a, fields in enumerate(be):
            counts = [(s, int(f.rstrip("c"))) for s, f in enumerate(fields) if f != "-"]
            longest = max((n for _, n in counts), default=0)
            short = [(n, -s) for s, n in counts if 2 <= n < longest]
            if short:
                n, s = max(short)
                pick = (a, -s, n, longest)
                break
        if pick is None:
            verdict(False, "no attempt of the mixed ensemble at %s had a system leave two or more passes before "
                    "the longest, so no system's own exit is held" % fmt)
            continue
        a, s, n, longest = pick
        mut = [list(f) for f in be]
        mut[a][s] = "1"
        path = os.path.join(work, "%s-one-system.sched" % fmt)
        write(path, he, mut)
        rc, out, err = run([cft] + configs[ke]["args"] + ["--replay-iters", path])
        def by_system(text):
            lines = {}
            for l in text.splitlines():
                if l.startswith("system "):
                    lines.setdefault(int(l.split()[1]), []).append(l)
            return lines
        ref = by_system(results[ke][0][1] if published else results[ke][1][1])
        new = by_system(out)
        moved = sorted(k for k in set(ref) | set(new) if ref.get(k) != new.get(k))
        verdict(rc == 0 and moved == [s],
                "%s: attempt %d, system %d's %d passes cut to 1 (the attempt's longest ran %d): %s"
                % (configs[ke]["label"], a, s, n, longest,
                   ("the systems whose record lines moved are %s, and only %d's may" % (moved, s)) if rc == 0 else
                   desync(rc, err)))
        #     "Only that system moved" cannot see a system frozen a pass late,
        #     because it moves too (V2's M2b). The exact reference can: the
        #     same system run ALONE (--member), whose schedule column must be
        #     the ensemble's attempt by attempt, with the same cut at its own
        #     attempt. An ensemble's system is its member alone, bit for bit
        #     (tools/check_ensemble.py), so the cut lines must be the same.
        margs = configs[ke]["args"] + ["--member", str(s)]
        msched = os.path.join(work, "%s-member.sched" % fmt)
        mrc, mout, merr = run([cft] + margs + ["--record-iters", msched])
        mh, mb = schedule(msched) if mrc == 0 else ([], [])
        own = [f[s] for f in be if f[s] != "-"]
        verdict(mrc == 0 and [f[0] for f in mb] == own,
                "system %d's passes in the ensemble's schedule are its own run alone (--member %d), attempt by attempt: "
                "%d of %d attempts" % (s, s, sum(1 for x, y in zip(mb, own) if x[0] == y), len(own)))
        am = sum(1 for f in be[:a] if f[s] != "-")      # the cut attempt, counted in the member's own attempts
        exact, detail = False, "the member's schedule has no attempt %d" % am
        if am < len(mb):
            mcut = [list(f) for f in mb]
            mcut[am][0] = "1"
            mpath = os.path.join(work, "%s-member-cut.sched" % fmt)
            write(mpath, mh, mcut)
            rrc, rout, rerr = run([cft] + margs + ["--replay-iters", mpath])
            alone = [l for l in rout.splitlines() if l.startswith("sample ")]
            mine = [l.split(" ", 2)[2] for l in new.get(s, [])]      # its lines, "system <s> " taken off
            if rc == 0 and rrc == 0:
                diff = sum(1 for x, y in zip(mine, alone) if x != y) + abs(len(mine) - len(alone))
                exact = diff == 0 and len(alone) > 0
                detail = "%d of %d lines differ" % (diff, len(alone))
            else:
                detail = "; ".join(desync(c, e, who) for who, c, e in (("the ensemble's replay", rc, err),
                                                                      ("the member's replay", rrc, rerr)) if c != 0)
        verdict(exact, "and its cut lines are the same cut of --member %d run alone (its attempt %d): %s" % (s, am, detail))

    # 5. Refusals by name, at the first format.
    fmt0 = formats[0]
    k = next(i for i, c in enumerate(configs) if c["fmt"] == fmt0 and "dt0=1.0" in c["label"])
    base = [cft] + configs[k]["args"]
    src = schedule_of(k)
    head, body = schedule(src)
    print("check_corrector at %s: refusals" % fmt0)
    other = [a if a != fmt0 else ("fp128" if fmt0 != "fp128" else "fp64") for a in base]
    for name, cmd, path, phrase in [
        ("a schedule from another run (another format)", other, src, "was recorded from a different run"),
        ("a schedule cut short", base, ("cut", body[:-1]), "the schedule ends after"),
        ("a schedule too long", base, ("long", body + [body[-1]]), "the schedule holds more"),
        ("a malformed field", base, ("bad", [["x"]] + body[1:]), "is not a pass count"),
        ("a system marked out of an attempt it is in", base, ("dash", [["-"]] + body[1:]), "was not in step attempt"),
    ]:
        if isinstance(path, tuple):
            p = os.path.join(work, "refuse-%s.sched" % path[0])
            write(p, head, path[1])
            path = p
        rc, out, err = run(cmd + ["--replay-iters", path])
        verdict(rc == 3 and phrase in err and not [l for l in out.splitlines() if l.startswith("# steps_done=")],
                "%s is refused by name, with no trailer (rc %d)" % (name, rc))
    rc, out, err = run(base + ["--replay-iters", src, "--trace-pc", "2"])
    verdict(rc == 3 and "--trace-pc and --replay-iters" in err, "--trace-pc beside --replay-iters is refused by name (rc %d)" % rc)
    rc, out, err = run(base + ["--replay-iters", src, "--record-iters", os.path.join(work, "x.sched")])
    verdict(rc == 3 and "--record-iters and --replay-iters" in err, "--record-iters beside --replay-iters is refused by name (rc %d)" % rc)
    rc, out, err = run(base + ["--replay-iters", os.path.join(pdir, "kepler.txt")])
    verdict(rc == 3 and "is not a corrector schedule" in err, "a file that is not a schedule is refused by name (rc %d)" % rc)

    # 6. The run key holds every field. A schedule replayed under a run that
    #    differs in one thing is refused by name, and the refusal's two key
    #    lines say which fields differ. Every field the key carries must be
    #    moved by some case here, read from the key itself, so a field added
    #    to it untested fails by name. Epsilon and the initial state each
    #    have a case that moves nothing else (V2's M3 and M3s: a key without
    #    either passed every other check, and replayed a schedule recorded
    #    at epsilon 1e-9 under 1e-8 with rc 0).
    print("check_corrector at %s: the run key holds every field" % fmt0)
    kep = os.path.join(pdir, "kepler.txt")
    kbase = ["--format", fmt0, "--quiet", "--dt", "0.01", "--epsilon", "1e-9", "--steps", "4", "--sample", "2"]
    ksched = os.path.join(work, "key.sched")
    rc, out, err = run([cft, "--problem", kep] + kbase + ["--record-iters", ksched])
    khead, _ = schedule(ksched)
    keyline = next((l[6:] for l in khead if l.startswith("# run ")), "")
    fields = set(kv(keyline))
    with open(kep, encoding="utf-8") as f:
        ktext = f.read()

    def variant(tag, edit):
        p = os.path.join(work, "key-%s.txt" % tag)
        with open(p, "w", encoding="utf-8", newline="\n") as f:
            f.write(edit(ktext))
        return p

    def nudge(tok):
        # the last body's field tok (2 its mass, 3 its x, 7 its vy), its final
        # hex digit moved: one value, about one ulp's worth of change
        def edit(text):
            lines = text.splitlines()
            i = max(k for k, l in enumerate(lines) if l.startswith("body "))
            t = lines[i].split()
            m, e = t[tok].split("p")
            t[tok] = m[:-1] + ("0" if m[-1] != "0" else "1") + "p" + e
            lines[i] = " ".join(t)
            return "\n".join(lines) + "\n"
        return edit

    # The state field is ONE hash, of the masses, the positions, the
    # velocities and a dt-file's steps (src/ias15_cft.c, sched_key()), so
    # each of the four is moved by a case of its own (V2's M3v, M3m, M3dt: a
    # hash without any one of them passed a gate that moved a position only).
    dtf, dtf2 = os.path.join(work, "key-dt.txt"), os.path.join(work, "key-dt2.txt")
    with open(dtf, "w", encoding="utf-8", newline="\n") as f:
        f.write("0x1.47ae147ae147bp-7\n" * 5)
    with open(dtf2, "w", encoding="utf-8", newline="\n") as f:
        f.write("0x1.47ae147ae147bp-7\n" * 2 + "0x1.47ae147ae147ap-7\n" + "0x1.47ae147ae147bp-7\n" * 2)
    replace = lambda a, b, over=kbase: [b if x == a else x for x in over]
    kbase_dt = replace("1e-9", "0") + ["--dt-file", dtf]
    ksched_dt = os.path.join(work, "key-dt.sched")
    rc, out, err = run([cft, "--problem", kep] + kbase_dt + ["--record-iters", ksched_dt])
    # (what the case is, the field it must move, the problem, the run, the
    #  only fields it may move or None, the schedule it replays)
    cases = [
        ("epsilon", "epsilon", kep, replace("1e-9", "1e-8"), {"epsilon"}, ksched),
        ("dt0", "dt0", kep, replace("0.01", "0.02"), {"dt0"}, ksched),
        ("steps", "steps", kep, replace("4", "5"), {"steps"}, ksched),
        ("sample", "sample", kep, replace("2", "1"), {"sample"}, ksched),
        ("adaptive_mode", "adaptive_mode", kep, kbase + ["--adaptive-mode", "3"], {"adaptive_mode"}, ksched),
        ("max_iter", "max_iter", kep, kbase + ["--max-iter", "11"], {"max_iter"}, ksched),
        ("pc_tol_shift", "pc_tol_shift", kep, kbase + ["--pc-tol-shift", "1"], {"pc_tol_shift"}, ksched),
        ("cs", "cs", kep, kbase + ["--cs", "augmented"], {"cs"}, ksched),
        ("arith", "arith", kep, kbase + ["--arith", "fma"], {"arith"}, ksched),
        ("engine", "engine", kep, kbase + ["--arith", "fma", "--engine", "program", "--programs", progs], None, ksched),
        ("softening2", "softening2", kep, kbase + ["--softening", "0.001"], {"softening2"}, ksched),
        ("min_dt", "min_dt", kep, kbase + ["--min-dt", "0.001"], {"min_dt"}, ksched),
        ("format", "format", kep, replace(fmt0, "fp128" if fmt0 != "fp128" else "fp64"), None, ksched),
        ("state: a position", "state", variant("x", nudge(3)), kbase, {"state"}, ksched),
        ("state: a velocity", "state", variant("v", nudge(7)), kbase, {"state"}, ksched),
        ("state: a mass", "state", variant("m", nudge(2)), kbase, {"state"}, ksched),
        ("state: a dt-file step", "state", kep, replace(dtf, dtf2, kbase_dt), {"state"}, ksched_dt),
        ("G", "G", variant("G", lambda t: re.sub(r"(?m)^G \S+$", "G 0x1.0000000000001p+0", t)), kbase, {"G"}, ksched),
        ("problem", "problem", variant("name", lambda t: re.sub(r"(?m)^name \S+$", "name kepler_renamed", t)), kbase,
         {"problem"}, ksched),
        ("N", "N", os.path.join(pdir, "outer.txt"), kbase, None, ksched),
        ("E", "E", ens, kbase, None, ksched),
        ("member", "member", ens, kbase + ["--member", "0"], None, ksched),
        ("dt_file", "dt_file", kep, kbase_dt, None, ksched),
    ]
    moved_all = set()
    for label, field, prob, a, only, sched in cases:
        rc, out, err = run([cft, "--problem", prob] + a + ["--replay-iters", sched])
        m = re.search(r"it was recorded from: (.*)\n\s*this run is:\s+(.*)", err)
        moved = sorted(k for k in set(kv(m.group(1))) | set(kv(m.group(2)))
                       if kv(m.group(1)).get(k) != kv(m.group(2)).get(k)) if m else []
        moved_all |= set(moved)
        verdict(rc == 3 and "was recorded from a different run" in err and field in moved
                and (only is None or set(moved) == only),
                "a schedule replayed under another %s is refused by name (rc %d); the keys differ in %s%s"
                % (label, rc, ", ".join(moved) or "nothing", "" if only is None else ", and in nothing else"))
    missing = sorted(fields - moved_all)
    verdict(fields and not missing, "every field of the run key is held by a case: %d of %d%s"
            % (len(fields) - len(missing), len(fields), " (UNHELD: %s)" % ", ".join(missing) if missing else ""))

    # 7. The library half on the device route, forced on. The route is
    # opt-in, so the suite's own check_dropin leg runs today's, and this is
    # the one run that takes the device's test there. It is also the only
    # state in which a library half silently on today's route differs from
    # a correct one. tools/cases_corrector.c asserts the route and its counts.
    if args.library and published:
        print("check_corrector: the drop-in's suite with the exit test forced onto the device")
        rc, out, err = run([os.path.join(args.build, "check_dropin" + exe), "--light"], {env: "on"})
        maxima = [int(c) for c in re.findall(r"^\s+ok\s+.*%s=on by=forced.*carried this case's (\d+) maxima" % re.escape(route),
                                             out, re.M)]
        passed = "check_dropin: every case passed" in out
        verdict(rc == 0 and passed and maxima and all(maxima),
                "check_dropin --light, the exit test forced onto the device: rc %d, %s; %d cases carried %d maxima on it"
                % (rc, "every case passed" if passed else "NOT every case passed", len(maxima), sum(maxima)))

    if ok_all and not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    elif not ok_all:
        print("the schedules and files of this run are kept in %s" % work)
    print("RESULT: %s" % ("PASS" if ok_all else "FAIL"))
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
