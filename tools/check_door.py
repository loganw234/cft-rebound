#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of the door (ROUND3.md, P0). ABI 0.12 to 0.14 opened three
# routes: a scalar operand, an indexed gather and a segmented reduction.
# Each is held against the route this port already takes, which the
# contract makes the same bits. Also held is the policy that chooses
# between them, which every round-3 parcel shares.
#
#   python tools/check_door.py [--build build] [--exe ias15_cft]
#                              [--formats fp64,fp128,fp256]
#   (`make check-quick` and `make check-light` pass --formats fp64)
#
# What it asserts, and what lets each assertion fail:
#   1. --door-selftest at each format:
#        - every case identical, and every control seeing its difference;
#        - BOTH routes of every published family counted as having run.
#          A self-test that took one route twice would be comparing that
#          route against itself.
#        - a family the device does not publish has every case reported
#          as not compared, and never as a pass.
#   2. With every bit hidden, no case is identical and every case is
#      not compared. A hidden route must never count as a pass.
#   3. A route forced on against a hidden bit is refused, and the
#      refusal names the variable and the bit. So are a bad value and a
#      bad route name.
#   4. A run with every route forced off gives the same record as the
#      default run. That is the fallback's own gate. The route lines say
#      how many calls each family carried, so a comparison that nothing
#      exercised shows up as one. In P0 nothing calls the door yet.
#   5. The record's `# identity` line names the library build.
#
# The route names, variables and bits are read from the program's own
# output, never typed here.

import argparse
import os
import re
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)

ROUTE_LINE = re.compile(r"^door route (\w+)=(on|off) by=(\w+) env=(\S+) bit=(\S+)$")
COUNT_LINE = re.compile(r"^door route (\w+) calls=(\d+) fallbacks=(\d+)$")
SUMMARY = re.compile(r"^door: (\d+) identical, (\d+) differ, (\d+) not compared; controls (\d+) of (\d+)$")
TRAILER_ROUTE = re.compile(r"^# route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)$")
COUNTERS = ("steps_done", "iterations_max_exceeded", "steps_rejected", "mean_pc_iterations",
            "max_pc_iterations", "flags_seen", "calls", "divsqrt_calls")


def run(cmd, extra_env=None):
    """Run with no door variable inherited, so the environment cannot
    choose a route behind the gate's back. CFT_REBOUND_ARTIFACT stays: a
    card run is the same gate."""
    env = {k: v for k, v in os.environ.items()
           if not k.startswith("CFT_REBOUND_ROUTE_") and k != "CFT_REBOUND_HIDE_FEATURES"}
    env.update(extra_env or {})
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    return p.returncode, p.stdout, p.stderr


def selftest(cft, fmt, extra_env=None):
    rc, out, err = run([cft, "--door-selftest", "--format", fmt], extra_env)
    routes, counts, summary, skips = {}, {}, None, {}
    for line in out.splitlines():
        m = ROUTE_LINE.match(line)
        if m:
            routes[m.group(1)] = {"on": m.group(2), "by": m.group(3), "env": m.group(4), "bit": m.group(5)}
            continue
        m = COUNT_LINE.match(line)
        if m:
            counts[m.group(1)] = (int(m.group(2)), int(m.group(3)))
            continue
        m = SUMMARY.match(line)
        if m:
            summary = tuple(int(g) for g in m.groups())
            continue
        if line.startswith("skip "):
            fam = line.split()[1]
            skips[fam] = skips.get(fam, 0) + 1
        elif line.startswith("DIFF ") or line.startswith("CONTROL FAILED"):
            print("       " + line)
    return rc, out, err, routes, counts, summary, skips


def record(cft, fmt, problem, extra, extra_env=None):
    rc, out, err = run([cft, "--format", fmt, "--problem", problem, "--quiet"] + extra, extra_env)
    samples = [line for line in out.splitlines() if line.startswith("sample ")]
    tail, routes, identity = {}, {}, None
    for line in out.splitlines():
        if line.startswith("# steps_done="):
            for kv in line[2:].split():
                k, _, v = kv.partition("=")
                tail[k] = v
        m = TRAILER_ROUTE.match(line)
        if m:
            routes[m.group(1)] = (m.group(2), m.group(3), int(m.group(4)), int(m.group(5)))
        if line.startswith("# identity "):
            identity = line
    return rc, err, samples, tail, routes, identity


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--exe", default="ias15_cft")
    ap.add_argument("--formats", default="fp64,fp128,fp256")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, args.exe + exe)
    formats = args.formats.split(",")
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= ok
        print("  %s %s" % ("ok  " if ok else "FAIL", text))

    # 1. The self-test at every format.
    routes = None
    cases = None
    for fmt in formats:
        print("check_door at %s: each route against today's route" % fmt)
        rc, out, err, routes, counts, summary, skips = selftest(cft, fmt)
        if rc not in (0, 1) or summary is None or not routes:
            verdict(False, "--door-selftest did not run (rc %d): %s" % (rc, (err or out)[-400:].strip()))
            continue
        same, differ, skipped, cok, cn = summary
        cases = same + differ + skipped
        present = [r for r, d in routes.items() if d["by"] != "absent"]
        verdict(rc == 0 and differ == 0, "%d of %d cases identical, %d differ (rc %d)" % (same, cases, differ, rc))
        verdict(cok == cn and cn >= len(present),
                "controls: %d of %d saw their difference, one or more for each of %d published families" % (cok, cn, len(present)))
        for r, d in routes.items():
            calls, falls = counts.get(r, (0, 0))
            if d["by"] == "absent":
                print("  note %s: %s is not published on this device, so its %d cases ran today's route only"
                      " and were not compared" % (r, d["bit"], skips.get(r, 0)))
                continue
            verdict(skips.get(r, 0) == 0 and calls > 0 and falls > 0,
                    "%s: both routes ran (%d calls through %s, %d through today's route), nothing skipped"
                    % (r, calls, d["bit"], falls))

    if routes is None or cases is None:
        print("RESULT: FAIL")
        return 1
    names = sorted(routes)
    f0 = formats[0]

    # 2. Every bit hidden: nothing may count as a pass.
    print("check_door at %s: every bit hidden" % f0)
    rc, out, err, hroutes, _, hsum, _ = selftest(cft, f0, {"CFT_REBOUND_HIDE_FEATURES": ",".join(names)})
    ok = (rc == 0 and hsum is not None and hsum[0] == 0 and hsum[2] == cases
          and all(d["by"] == "hidden" and d["on"] == "off" for d in hroutes.values()))
    verdict(ok, "hidden: %s identical, %s not compared of %d, every route off by=hidden"
            % (hsum[0] if hsum else "?", hsum[2] if hsum else "?", cases))

    # 3. Refusals, each by name.
    print("check_door at %s: refusals" % f0)
    for r in names:
        env, bit = routes[r]["env"], routes[r]["bit"]
        rc, out, err = run([cft, "--door-selftest", "--format", f0],
                           {"CFT_REBOUND_HIDE_FEATURES": r, env: "on"})
        verdict(rc != 0 and (env + "=on") in err and bit in err and "does not publish" in err,
                "%s=on against a hidden %s is refused by name (rc %d)" % (env, bit, rc))
    env0 = routes[names[0]]["env"]
    rc, out, err = run([cft, "--door-selftest", "--format", f0], {env0: "yes"})
    verdict(rc != 0 and "expected on, off or auto" in err, "%s=yes is refused by name (rc %d)" % (env0, rc))
    rc, out, err = run([cft, "--door-selftest", "--format", f0], {"CFT_REBOUND_HIDE_FEATURES": "no_such_route"})
    verdict(rc != 0 and "is not a route" in err, "CFT_REBOUND_HIDE_FEATURES=no_such_route is refused by name (rc %d)" % rc)

    # 4. Every route forced off against the default: the same record.
    kepler = os.path.join(root, "data", "problems", "kepler.txt")
    outer = os.path.join(root, "data", "problems", "outer.txt")
    runs = [("kepler adaptive", kepler, ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "50", "--sample", "25"]),
            ("outer fixed", outer, ["--dt", "40", "--epsilon", "0", "--steps", "20", "--sample", "10"])]
    off = {routes[r]["env"]: "off" for r in names}
    for fmt in formats:
        print("check_door at %s: every route forced off against the default" % fmt)
        for name, prob, extra in runs:
            rc0, err0, s0, t0, r0, ident = record(cft, fmt, prob, extra)
            rc1, err1, s1, t1, r1, _ = record(cft, fmt, prob, extra, off)
            if rc0 or rc1:
                verdict(False, "%s did not run (rc %d/%d): %s" % (name, rc0, rc1, (err0 or err1)[-300:].strip()))
                continue
            carried = sum(v[2] + v[3] for v in list(r0.values()) + list(r1.values()))
            same = s0 == s1 and len(s0) > 0 and all(t0.get(k) == t1.get(k) for k in COUNTERS)
            forced = sorted(r1) == names and all(v[0] == "off" and v[1] == "forced" for v in r1.values())
            verdict(same and forced and sorted(r0) == names,
                    "%s: %d sample lines and the counters identical; forced-off run reports every route off by=forced%s"
                    % (name, len(s0), "" if carried else
                       " (vacuous for now: no route carried a call in either run, which is P0's state)"))
            if fmt == f0 and name == runs[0][0]:
                verdict(ident is not None and "libcft_build=" in ident and (" image=none" in ident or " image_sha256=" in ident),
                        "the record names its library build: %s" % (ident or "no # identity line"))

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
