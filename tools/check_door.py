#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of the door (ROUND3.md, P0). ABI 0.12 to 0.14 opened three
# routes: a scalar operand, an indexed gather and a segmented reduction.
# Each route, and the route this port already takes, is held against a
# reference computed outside the wrappers. Also held is the policy that
# chooses between them, which every round-3 parcel shares.
#
#   python tools/check_door.py [--build build] [--exe ias15_cft]
#                              [--formats fp64,fp128,fp256]
#   (`make check-quick` and `make check-light` pass --formats fp64)
#
# What it asserts, and what lets each assertion fail:
#   1. --door-selftest at each format:
#        - every case identical, and every control seeing its difference;
#        - every wrapper call accounted exactly, both in calls= and in
#          the route's own counters;
#        - something compared at all. A device that publishes nothing
#          passes nothing.
#        - on the software backend, every bit published, because it
#          publishes them all (cft.h). A probe that reads none is broken,
#          and it passed every check here before 2026-10-31 (V0's N4).
#          So every route is by=device there, except the opt-in routes
#          (OPT_IN), which auto leaves off and which say by=opt_in.
#        - a family that a card does not publish has every case reported
#          as not compared, and never as a pass.
#   2. With every bit hidden, no case is identical and every case is not
#      compared.
#   3. Refusals, each by name:
#        - a route forced on against a hidden bit;
#        - a bad value, a bad route name, and a misspelled variable;
#        - every --door-refuse call, on BOTH routes, so that today's route
#          never computes what the new one refuses. Where the two refuse
#          in different words, each route is held to its own, which says
#          which route ran.
#   4. A run with every route forced off gives the same record as the
#      default run, with every counter present on both sides. That is
#      the fallback's own gate. The route lines say how many calls each
#      family carried, so a comparison that nothing exercised shows up as
#      one. In P0 nothing calls the door yet.
#      The opt-in routes are off by default, so that comparison no longer
#      reaches them, and the opposite one does: a run with every opt-in
#      route forced on, against the default, compared the same way, with
#      each forced route on by=forced and carrying calls, and each off in
#      the default by=opt_in and carrying fallbacks. Where the device does
#      not publish an opt-in route's bit, forcing it on is refused, and the
#      comparison is skipped and says so.
#   5. The record's `# identity` line names the library build.
#
# The route names, variables, bits and refusal cases are read from the
# program's own output and source, never typed here, except OPT_IN: which
# routes are opt-in is the owner's decision (2026-11-02), and a program
# that changed it must fail here rather than be read.

import argparse
import os
import re
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)

ROUTE_LINE = re.compile(r"^door route (\w+)=(on|off) by=(\w+) env=(\S+) bit=(\S+)$")
COUNT_LINE = re.compile(r"^door route (\w+) calls=(\d+) fallbacks=(\d+)$")
SUMMARY = re.compile(r"^door: (\d+) identical, (\d+) differ, (\d+) not compared, (\d+) accounting errors; "
                     r"controls (\d+) of (\d+)$")
TRAILER_ROUTE = re.compile(r"^# route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)$")
COUNTERS = ("steps_done", "iterations_max_exceeded", "steps_rejected", "mean_pc_iterations",
            "max_pc_iterations", "flags_seen", "calls", "divsqrt_calls")
# The opt-in routes: auto leaves each off even where the device publishes
# its bit, by=opt_in, and CFT_REBOUND_ROUTE_<NAME>=on takes it. They are
# gravity's two switches (ROUND3.md, P1). On the rev8a tile each made a
# step 22-40% slower at 64 to 512 bodies with the same bits (lead.md,
# 2026-11-02 16:55), so the owner made them opt-in.
OPT_IN = ("indexed", "pair_gather")
# What each --door-refuse case must say, and which route it belongs to.
# The cases themselves are read from the program's own refusal of an
# unknown one.
REFUSALS = {"ix-alias-dense": ("indexed", "once any operand has an index table"),
            "ix-alias-source": ("indexed", "indexed source"),
            "ix-table-on-null": ("indexed", "has an index table and no source"),
            "ix-zero-source": ("indexed", "beside a source of zero elements"),
            "reduce-partial": ("reduce_seg", "are not whole segments"),
            "px-alias-pool": ("indexed", "may write over none of its inputs"),
            "px-stream-a-null": ("indexed", "stream a is NULL, and a program run reads it into r0"),
            "px-table-on-null": ("indexed", "has an index table and no source, which cft_program_run_ex refuses"),
            "px-zero-source": ("indexed", "beside a source of zero elements, which cft_program_run_ex refuses"),
            "px-table-no-block": ("indexed", "declares no scratch input, which cft_program_run_ex refuses"),
            # An index past its source: the new route refuses it in libcft's
            # words and today's in the engine's, each impossible on the other
            # route. So a route that says it is on while today's runs fails
            # here (V1's m4_lie, 2026-11-02), and so does today's route with
            # its bound check gone.
            "px-index-past-source": ("indexed", {"on": "idx_scratch_src says the pool holds",
                                                 "off": "is past its source of"}),
            "px-pair-past-source": ("pair_gather", {"on": "idx_scratch_src says the pool holds",
                                                    "off": "is past its source of"})}
DATA = re.compile(r"^door data: (\d+) of (\d+) nonzero values carry bits below binary64's reach$", re.M)
POSITIONS = re.compile(r"^door data: (\d+) of (\d+) bit positions below binary64's reach are set in some value$", re.M)
ZEROS = re.compile(r"^door data: zeros (\d+) positive and (\d+) negative$", re.M)


def run(cmd, extra_env=None):
    """Run with no door variable inherited, so the environment cannot
    choose a route behind the gate's back. CFT_REBOUND_ARTIFACT stays: a
    card run is the same gate."""
    env = {k: v for k, v in os.environ.items()
           if not k.upper().startswith("CFT_REBOUND_ROUTE_") and not k.upper().startswith("CFT_REBOUND_HIDE")}
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
        elif line.startswith(("DIFF ", "COUNT ", "CONTROL FAILED")):
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
    software = not os.environ.get("CFT_REBOUND_ARTIFACT")
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= ok
        print("  %s %s" % ("ok  " if ok else "FAIL", text))

    # 1. The self-test at every format.
    routes = None
    cases = None
    for fmt in formats:
        print("check_door at %s: each route and today's against a reference computed outside both" % fmt)
        rc, out, err, routes, counts, summary, skips = selftest(cft, fmt)
        if rc not in (0, 1) or summary is None or not routes:
            verdict(False, "--door-selftest did not run (rc %d): %s" % (rc, (err or out)[-400:].strip()))
            routes = None
            continue
        same, differ, skipped, bad, cok, cn = summary
        cases = same + differ + skipped
        present = [r for r, d in routes.items() if d["by"] != "absent"]
        verdict(rc == 0 and differ == 0 and same > 0,
                "%d of %d cases identical to the reference on both routes, %d differ (rc %d)" % (same, cases, differ, rc))
        verdict(bad == 0, "every wrapper call accounted: one operation in calls=, and the one route counter meant (%d errors)" % bad)
        # The data's own check. Below binary64 there is nothing to carry;
        # above it, a value without low bits cannot catch a route that
        # drops them (V0's F4).
        dm = DATA.search(out)
        beyond, values = (int(dm.group(1)), int(dm.group(2))) if dm else (-1, 0)
        wide = fmt != "fp64"
        verdict(values > 0 and beyond == (values if wide else 0),
                "the data: %d of %d nonzero values carry bits below binary64's reach (wanted %s)"
                % (beyond, values, "all of them" if wide else "none, at binary64"))
        # Every such bit position set somewhere, the positions measured by
        # the program (V0's F15 and F16), and zeros of both signs present
        # (V0's F18).
        pm, zm = POSITIONS.search(out), ZEROS.search(out)
        seen, posns = (int(pm.group(1)), int(pm.group(2))) if pm else (-1, -2)
        zpos, zneg = (int(zm.group(1)), int(zm.group(2))) if zm else (0, 0)
        verdict(seen == posns and (posns > 0) == wide,
                "the data: %d of %d bit positions below binary64's reach are set in some value" % (seen, posns))
        verdict(zpos > 0 and zneg > 0, "the data: zeros of both signs, %d positive and %d negative" % (zpos, zneg))
        verdict(cok == cn and cn >= len(present),
                "controls: %d of %d saw their difference, one or more for each of %d published families" % (cok, cn, len(present)))
        if software:
            # every bit published: each route on by=device, but an opt-in
            # route off by=opt_in, which only a published bit can say
            want = {r: ("off", "opt_in") if r in OPT_IN else ("on", "device") for r in routes}
            verdict(all((d["on"], d["by"]) == want[r] for r, d in routes.items()) and all(g in routes for g in OPT_IN),
                    "the software backend publishes every bit, and auto leaves only the opt-in routes off: %s"
                    % ", ".join("%s=%s by=%s" % (r, d["on"], d["by"]) for r, d in sorted(routes.items())))
        for r, d in sorted(routes.items()):
            calls, falls = counts.get(r, (0, 0))
            if d["by"] == "absent":
                print("  note %s: %s is not published on this device, so its %d cases ran today's route only"
                      " and were not compared" % (r, d["bit"], skips.get(r, 0)))
                continue
            verdict(skips.get(r, 0) == 0 and calls > 0 and falls > 0,
                    "%s: both routes ran (%d calls through %s, %d through today's route), nothing left uncompared"
                    % (r, calls, d["bit"], falls))

    if routes is None or cases is None:
        print("RESULT: FAIL")
        return 1
    names = sorted(routes)
    f0 = formats[0]

    # 2. Every bit hidden: nothing may count as a pass.
    print("check_door at %s: every bit hidden" % f0)
    rc, out, err, hroutes, _, hsum, _ = selftest(cft, f0, {"CFT_REBOUND_HIDE_FEATURES": ",".join(names)})
    ok = (rc == 0 and hsum is not None and hsum[0] == 0 and hsum[2] == cases and hsum[3] == 0
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
    typo = env0[:-1]          # a near miss: the route's own variable, one letter short
    rc, out, err = run([cft, "--door-selftest", "--format", f0], {typo: "on"})
    verdict(rc != 0 and typo in err and "is not one of the door's variables" in err,
            "a misspelled %s=on is refused by name, not ignored (rc %d)" % (typo, rc))
    hide_typo = "CFT_REBOUND_HIDE_FEATURES"[:-1]      # one letter short, as above
    rc, out, err = run([cft, "--door-selftest", "--format", f0], {hide_typo: names[0]})
    verdict(rc != 0 and hide_typo in err and "is not one of the door's variables" in err,
            "a misspelled %s is refused by name, not ignored (rc %d)" % (hide_typo, rc))
    rc, out, err = run([cft, "--door-refuse", "no_such_case", "--format", f0])
    listed = set(re.findall(r"[a-z]+(?:-[a-z]+)+", err.split("the cases are", 1)[-1])) if "the cases are" in err else set()
    verdict(rc != 0 and listed == set(REFUSALS),
            "--door-refuse names its cases, and they are the ones this gate holds: %s" % ", ".join(sorted(listed)))
    for case, (route, phrase) in sorted(REFUSALS.items()):
        if routes.get(route, {}).get("by") == "absent":
            print("  note %s: %s is absent here, so only today's route was asked to refuse it" % (case, route))
        for state in ("off", "on"):
            if state == "on" and routes.get(route, {}).get("by") == "absent":
                continue
            rc, out, err = run([cft, "--door-refuse", case, "--format", f0], {routes[route]["env"]: state})
            words = phrase[state] if isinstance(phrase, dict) else phrase
            verdict(rc == 3 and words in err and "NOT refused" not in out,
                    "%s is refused by name with %s=%s%s (rc %d)%s"
                    % (case, routes[route]["env"], state,
                       ", in the words of the route that must have run" if isinstance(phrase, dict) else "", rc,
                       "" if rc == 3 and words in err else ": wanted \"%s\", got: %s" % (words, err.strip()[-200:])))

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
            missing = [k for k in COUNTERS if k not in t0 or k not in t1]
            carried = sum(v[2] + v[3] for v in list(r0.values()) + list(r1.values()))
            same = s0 == s1 and len(s0) > 0 and not missing and all(t0[k] == t1[k] for k in COUNTERS)
            forced = sorted(r1) == names and all(v[0] == "off" and v[1] == "forced" for v in r1.values())
            verdict(same and forced and sorted(r0) == names,
                    "%s: %d sample lines and all %d counters identical%s; forced-off run reports every route off by=forced%s"
                    % (name, len(s0), len(COUNTERS), " (MISSING: %s)" % ", ".join(missing) if missing else "",
                       "" if carried else
                       " (vacuous for now: no route carried a call in either run, which is P0's state)"))
            if fmt == f0 and name == runs[0][0]:
                verdict(ident is not None and "libcft_build=" in ident and (" image=none" in ident or " image_sha256=" in ident),
                        "the record names its library build: %s" % (ident or "no # identity line"))

    # 4b. The opt-in routes forced on against the default. The default
    # leaves them off, so part 4 compares today's route with today's route
    # for them; this is the comparison that reaches the new one. The other
    # routes must say the same in both runs.
    opt = [r for r in names if r in OPT_IN]
    on_env = {routes[r]["env"]: "on" for r in opt}
    unpublished = [r for r in opt if routes[r]["by"] == "absent"]
    for fmt in formats if opt else []:
        print("check_door at %s: every opt-in route forced on against the default" % fmt)
        if unpublished:
            print("  note %s: not published on this device, so forcing on is refused, and today's route is all "
                  "there is to compare" % ", ".join(unpublished))
            continue
        for name, prob, extra in runs:
            rc0, err0, s0, t0, r0, _ = record(cft, fmt, prob, extra)
            rc1, err1, s1, t1, r1, _ = record(cft, fmt, prob, extra, on_env)
            if rc0 or rc1:
                verdict(False, "%s did not run (rc %d/%d): %s" % (name, rc0, rc1, (err0 or err1)[-300:].strip()))
                continue
            missing = [k for k in COUNTERS if k not in t0 or k not in t1]
            differ = [k for k in COUNTERS if k not in missing and t0[k] != t1[k]]
            same = s0 == s1 and len(s0) > 0 and not missing and not differ
            dflt = all(r in r0 and r0[r][:2] == ("off", "opt_in") and r0[r][2] == 0 and r0[r][3] > 0 for r in opt)
            forced = all(r in r1 and r1[r][:2] == ("on", "forced") and r1[r][2] > 0 and r1[r][3] == 0 for r in opt)
            others = sorted(r0) == names and sorted(r1) == names and all(r0[r] == r1[r] for r in names if r not in opt)
            verdict(same and dflt and forced and others,
                    "%s: %d sample lines and all %d counters identical%s%s%s; the default leaves %s off by=opt_in "
                    "(fallbacks %s), forced on they carry calls %s and no fallback%s"
                    % (name, len(s0), len(COUNTERS), " (MISSING: %s)" % ", ".join(missing) if missing else "",
                       "" if s0 == s1 else " (THE SAMPLE LINES DIFFER)",
                       "" if not differ else " (DIFFER: %s)" % ", ".join("%s %s against %s" % (k, t0[k], t1[k]) for k in differ),
                       " and ".join(opt), "/".join(str(r0[r][3]) if r in r0 else "?" for r in opt),
                       "/".join(str(r1[r][2]) if r in r1 else "?" for r in opt),
                       "" if dflt and forced and others else
                       " (ROUTES: default %s; forced %s)" % (", ".join("%s=%s by=%s calls=%d fallbacks=%d" % ((r,) + r0[r]) for r in sorted(r0)),
                                                             ", ".join("%s=%s by=%s calls=%d fallbacks=%d" % ((r,) + r1[r]) for r in sorted(r1)))))

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
