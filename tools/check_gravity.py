#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of gravity's routes (ROUND3.md, P1).
#
# Gravity is a handful of sequencer program runs a force evaluation now,
# with two switches on the door: `indexed`, whether the accumulate fold's
# addends are gathered by the device through a table or by the host, and
# `pair_gather`, the same choice for the pair program's six coordinates.
# On either route the same programs run on the same bits. Both switches
# are opt-in (src/ias15_cft.c, engine_routes()): unforced, each is off,
# by=opt_in, even where the device publishes the bit, because on the
# rev8a tile each made a step slower (lead.md, 2026-11-02 16:55). The
# device's gather runs only where a variable asks for it.
#
#   python tools/check_gravity.py [--build build] [--quick] [--library]
#   (`make check` passes --library; check-quick and check-light pass --quick
#   --library. The routes are opt-in, so --library is the one leg in any tier
#   that can see a library half silently on today's route.)
#
# What it asserts, and what lets each assertion fail:
#   1. The programs this binary carries are the ones programs/out/ holds,
#      byte for byte, by SHA-256 and in both directions, at every format:
#      the binary reports them itself (--door-selftest's `door image`
#      lines, through cft_program_digest), and a binary built before
#      `make programs` regenerated them fails here by name.
#   2. Every combination of the two switches gives the same record at
#      every format, on problems whose rows differ - two bodies, three,
#      six with a massless one, an ensemble, 33 bodies, a rejected step,
#      and (at binary64) rows one past the widest fold, which take two
#      runs a force evaluation. Every run's two `# route` lines are
#      asserted: on or off, by= (forced; for the default run, off by=opt_in,
#      or on a card off by=absent), calls= above zero where the route ran
#      and zero where it did not, and fallbacks= the reverse. The physics
#      below reads its counts from the default run, today's routes, whose
#      fallbacks count the same runs. The trailer's counters must agree
#      across the combinations too, calls= and program_calls= included:
#      both routes run one program per call.
#      Beside each, the physics it had to have: steps taken, and the fold
#      runs a force evaluation, which must be the rows' length over the
#      widest fold, rounded up - so the chained case is chained - and, for
#      the rejected step's case, steps_rejected above zero. One problem
#      has test particles (an `active` line): its rows are N_active long,
#      read from the record's own `# active` line.
#   2b. Test particles (round 4's P0, ROUND4.md): outer with two more
#      massless bodies, run with `active 5` and with every body active,
#      gives the same record at every format, on the loop engine in both
#      forms and on the program engine. The record says `active 5 of 8`
#      with 25 pairs a force evaluation (REBOUND's two loops: 10 active
#      pairs, then 3 x 5), the test particles move, and the comparison can
#      fail: with one of them given a mass of 2^-40 the two runs differ.
#      testparticle_type 0 is held too: under `active 5` that 2^-40 mass
#      must change NOTHING, the energy column included - a test particle's
#      mass moves no active body (type 1 would), and REBOUND's energy
#      counts only the active bodies (reb_simulation_energy, N_interact).
#      V0 found both unheld in P0's first version, whose test particles
#      were all massless. And an ensemble of two systems with `active 5`
#      against every body active, and its member against the system alone.
#   2c. The problem lines' refusals: an active count of 0 or past N, a
#      second active line, a zonal line (until round 4's P1 builds the
#      term), and a corrector schedule recorded with `active 5` replayed
#      without it - the run key names the line, so a schedule replays only
#      its own N_active.
#   3. The body cap's refusal still holds one past it. The cap is read
#      from src/ias15_limits.h, never typed here. (Its one force
#      evaluation, at the cap, on the routes against the elementwise
#      reference, is --door-selftest's, held by tools/check_door.py.)
#   4. --library: the drop-in's whole suite (check_dropin --light) with
#      both gravity routes forced ON passes against REBOUND, and its own
#      gravity topic reads both routes on by=forced through
#      ias15_engine_route(), with calls and no fallback. The suite's own
#      check_dropin leg is the today's-routes half, since the routes are
#      opt-in. A library half silently on today's route (M4) passes that
#      leg, which is right for it, and fails here: forced on is the one
#      state that tells them apart. A device without the bit refuses the
#      routes forced on, here as in part 2, so this gate needs one that
#      publishes it.
#
# The route names, their variables and the hide list are read from the
# program's own output, never typed here, except the two gravity names
# this gate exists for, which it checks the program has.

import argparse
import hashlib
import os
import re
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)

ROUTE_LINE = re.compile(r"^door route (\w+)=(on|off) by=(\w+) env=(\S+) bit=(\S+)$")
IMAGE_LINE = re.compile(r"^door image (\S+)-(fp\d+) bytes=(\d+) sha256=([0-9a-f]{64})$")
TRAILER_ROUTE = re.compile(r"^# route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)$")
COUNTERS = ("steps_done", "iterations_max_exceeded", "steps_rejected", "mean_pc_iterations",
            "max_pc_iterations", "flags_seen", "calls", "divsqrt_calls", "program_calls")
GRAVITY = ("indexed", "pair_gather")
COMBOS = (("on", "on"), ("off", "off"), ("on", "off"), ("off", "on"))


def clean_env(extra=None):
    """No door variable inherited: the environment may not choose a route
    behind the gate's back. CFT_REBOUND_ARTIFACT stays."""
    env = {k: v for k, v in os.environ.items()
           if not k.upper().startswith("CFT_REBOUND_ROUTE_") and not k.upper().startswith("CFT_REBOUND_HIDE")}
    env.update(extra or {})
    return env


def run(cmd, extra=None):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=clean_env(extra), cwd=root)
    return p.returncode, p.stdout, p.stderr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true", help="binary64 only, and fewer steps")
    ap.add_argument("--library", action="store_true", help="also run check_dropin --light with the gravity routes off")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, "ias15_cft" + exe)
    formats = ["fp64"] if args.quick else ["fp64", "fp128", "fp256"]
    software = not os.environ.get("CFT_REBOUND_ARTIFACT")
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= ok
        print("  %s %s" % ("ok  " if ok else "FAIL", text))

    # The routes, from the program: --probe names them and the device's
    # depth, and a self-test with every bit hidden runs no case at all but
    # still prints each route's variable and every image it carries.
    rc, out, err = run([cft, "--probe", "--format", "fp64"])
    names = re.findall(r"^route (\w+): ", out, re.M)
    depth = re.search(r"^max_scratch: (\d+)$", out, re.M)
    depth = int(depth.group(1)) if depth else 0
    hide_all = {"CFT_REBOUND_HIDE_FEATURES": ",".join(sorted(names))}
    rc, out, err = run([cft, "--door-selftest", "--format", "fp64"], hide_all)
    routes = {}
    for line in out.splitlines():
        m = ROUTE_LINE.match(line)
        if m:
            routes[m.group(1)] = {"env": m.group(4), "bit": m.group(5)}
    if not names or sorted(routes) != sorted(names) or any(g not in routes for g in GRAVITY):
        print("check_gravity: the program does not name the routes %s (--probe names %s, the self-test %s; rc %d): %s"
              % (", ".join(GRAVITY), ", ".join(names) or "none", ", ".join(sorted(routes)) or "none", rc,
                 (err or out)[-300:].strip()))
        print("RESULT: FAIL")
        return 1
    widths = [int(m.group(1)) for m in re.finditer(r"^door image gravity-fold(\d+)-fp64 ", out, re.M)]
    widths = [w for w in widths if not depth or w <= depth]
    wmax = max(widths) if widths else 0
    print("check_gravity: routes %s; the fold's widest admitted width is %d (%d carried, a scratch depth of %d)"
          % (", ".join("%s (%s)" % (g, routes[g]["env"]) for g in GRAVITY), wmax, len(widths), depth))

    # 1. The images.
    print("check_gravity: the gravity programs compiled in, against programs/out/")
    outdir = os.path.join(root, "programs", "out")
    for fmt in formats:
        rc, out, err = run([cft, "--door-selftest", "--format", fmt], hide_all)
        carried = {}
        for line in out.splitlines():
            m = IMAGE_LINE.match(line)
            if m and m.group(2) == fmt:
                carried[m.group(1)] = (int(m.group(3)), m.group(4))
        onfile = {}
        for name in sorted(os.listdir(outdir)):
            if name.startswith("gravity-") and name.endswith("-%s.cftp" % fmt):
                with open(os.path.join(outdir, name), "rb") as f:
                    data = f.read()
                onfile[name[:-len("-%s.cftp" % fmt)]] = (len(data), hashlib.sha256(data).hexdigest())
        missing = sorted(set(onfile) - set(carried))
        extra = sorted(set(carried) - set(onfile))
        stale = sorted(s for s in set(onfile) & set(carried) if onfile[s] != carried[s])
        verdict(rc == 0 and onfile and not missing and not extra and not stale,
                "%s: %d images carried, %d in programs/out, all %s by SHA-256%s"
                % (fmt, len(carried), len(onfile), "equal" if not stale else "NOT equal",
                   "".join(" (%s: %s)" % (k, ", ".join(v)) for k, v in
                           (("not carried", missing), ("not on file", extra), ("stale, rebuild", stale)) if v)))

    # 2. The route matrix.
    probs = os.path.join(root, "data", "problems")
    scratch = os.path.join(args.build, "check_gravity")
    os.makedirs(scratch, exist_ok=True)
    nb33 = os.path.join(scratch, "nbody33.txt")
    nb258 = os.path.join(scratch, "nbody258.txt")
    ens = os.path.join(scratch, "outer4.txt")
    py = sys.executable
    for n, path in ((33, nb33), (258, nb258)):
        if not os.path.exists(path):
            subprocess.run([py, os.path.join(here, "make_nbody.py"), str(n), path], check=True, cwd=root,
                           stdout=subprocess.DEVNULL)
    if not os.path.exists(ens):
        subprocess.run([py, os.path.join(here, "make_ensemble.py"), os.path.join(probs, "outer.txt"),
                        "--members", "4", "--out", ens, "--geometric", "2", "x", "0.25"], check=True, cwd=root,
                       stdout=subprocess.DEVNULL)
    # Round 4's P0 problems: outer with two more massless bodies, Pluto
    # turned a quarter and a half turn about z - a swap and negations, so
    # every value is exactly outer's own - with and without `active 5`,
    # with one of them given a mass of 2^-40 (the control), and two that
    # must be refused.
    with open(os.path.join(probs, "outer.txt")) as f:
        outer_lines = f.read().splitlines()
    pl = [l for l in outer_lines if l.startswith("body pluto ")][0].split()
    px, py_, pz, pvx, pvy, pvz = pl[3:9]
    neg = lambda h: h[1:] if h.startswith("-") else "-" + h

    def outer8_lines(tp2_mass="0x0.0p+0"):
        lines = ["N 8" if l.startswith("N ") else l for l in outer_lines]
        lines.append("body tp1 0x0.0p+0 %s %s %s %s %s %s" % (neg(py_), px, pz, neg(pvy), pvx, pvz))
        lines.append("body tp2 %s %s %s %s %s %s %s" % (tp2_mass, neg(px), neg(py_), pz, neg(pvx), neg(pvy), pvz))
        return lines

    def write(name, lines):
        path = os.path.join(scratch, name)
        with open(path, "w") as f:
            f.write("\n".join(lines) + "\n")
        return path

    def outer8(name, active=None, tp2_mass="0x0.0p+0", extra=()):
        return write(name, outer8_lines(tp2_mass) + (["active %d" % active] if active is not None else []) + list(extra))

    def outer8x2(name, active=None):
        # two systems of the same eight bodies: the header lines, then each
        # system's body lines under its own system line
        lines = outer8_lines()
        head = [l for l in lines if not l.startswith("body ")]
        bodies = [l for l in lines if l.startswith("body ")]
        return write(name, head + ["E 2", "system one"] + bodies + ["system two"] + bodies
                     + (["active %d" % active] if active is not None else []))
    p_a5 = outer8("outer8_active5.txt", 5)
    p_all = outer8("outer8_all.txt")
    p_a5m = outer8("outer8_active5_mass.txt", 5, "0x1.0000000000000p-40")
    p_allm = outer8("outer8_all_mass.txt", None, "0x1.0000000000000p-40")
    p_ens5 = outer8x2("outer8x2_active5.txt", 5)
    p_ensall = outer8x2("outer8x2_all.txt")
    p_a9 = outer8("outer8_active9.txt", 9)
    p_a0 = outer8("outer8_active0.txt", 0)
    p_twice = outer8("outer8_active_twice.txt", 5, extra=("active 5",))
    p_zonal = outer8("outer8_zonal.txt", 5, extra=("zonal 0 0x1.0000000000000p+0 0x1.0000000000000p-10 0x0.0p+0",))

    s = 1 if args.quick else 2
    rejecting = "kepler e=0.99, a rejected step"     # here for its rejection, which is asserted below
    problems = [
        ("kepler, two bodies, adaptive", os.path.join(probs, "kepler.txt"), ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(40 * s), "--sample", "20"], None),
        (rejecting, os.path.join(probs, "kepler_e099.txt"), ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(30 * s), "--sample", "10"], None),
        ("pythagorean, three bodies", os.path.join(probs, "pythagorean.txt"), ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(40 * s), "--sample", "20"], None),
        ("outer, six bodies, pluto massless", os.path.join(probs, "outer.txt"), ["--dt", "40", "--epsilon", "0", "--steps", str(20 * s), "--sample", "10"], None),
        ("outer and two more, three test particles (active 5)", p_a5, ["--dt", "40", "--epsilon", "0", "--steps", str(20 * s), "--sample", "10"], None),
        ("outer x4, an ensemble", ens, ["--dt", "40", "--epsilon", "1e-9", "--steps", str(10 * s), "--sample", "10"], None),
        ("33 bodies, rows of 32", nb33, ["--dt", "0.01", "--epsilon", "1e-9", "--steps", str(2 * s), "--sample", "2"], None),
        ("258 bodies, rows one past the widest fold", nb258, ["--dt", "0.01", "--epsilon", "0", "--steps", "1", "--sample", "1", "--max-iter", "2"], ("fp64",)),
    ]
    for fmt in formats:
        print("check_gravity at %s: every combination of indexed and pair_gather, the same record" % fmt)
        for name, prob, extra, only in problems:
            if only and fmt not in only:
                continue
            base = [cft, "--format", fmt, "--problem", prob, "--quiet"] + extra
            results = []
            for combo in (None,) + COMBOS:
                env = {} if combo is None else {routes["indexed"]["env"]: combo[0], routes["pair_gather"]["env"]: combo[1]}
                rc, out, err = run(base, env)
                recs = [l for l in out.splitlines() if l and not l.startswith("#")]
                tail, rl = {}, {}
                for line in out.splitlines():
                    if line.startswith("# steps_done="):
                        for kv in line[2:].split():
                            k, _, v = kv.partition("=")
                            tail[k] = v
                    m = TRAILER_ROUTE.match(line)
                    if m:
                        rl[m.group(1)] = (m.group(2), m.group(3), int(m.group(4)), int(m.group(5)))
                results.append((combo, rc, err, recs, tail, rl, out))
            bad = []
            ref = results[0]
            for combo, rc, err, recs, tail, rl, _ in results:
                label = "default" if combo is None else "indexed=%s pair_gather=%s" % combo
                if rc != 0:
                    bad.append("%s: rc %d: %s" % (label, rc, err.strip()[-200:]))
                    continue
                if recs != ref[3] or not recs:
                    bad.append("%s: the record differs from the default run's" % label)
                for k in COUNTERS:
                    if tail.get(k) is None or tail.get(k) != ref[4].get(k):
                        bad.append("%s: %s=%s against the default run's %s" % (label, k, tail.get(k), ref[4].get(k)))
                for i, g in enumerate(GRAVITY):
                    if g not in rl:
                        bad.append("%s: no # route %s line" % (label, g))
                        continue
                    on, by, calls, falls = rl[g]
                    if combo is not None:
                        want_on, want_by = combo[i] == "on", ("forced",)
                    elif software:
                        # every bit published (cft.h), and both routes opt-in
                        want_on, want_by = False, ("opt_in",)
                    else:
                        want_on, want_by = False, ("opt_in", "absent")
                    ok = (on == ("on" if want_on else "off") and by in want_by and
                          ((calls > 0 and falls == 0) if want_on else (falls > 0 and calls == 0)))
                    if not ok:
                        bad.append("%s: # route %s=%s by=%s calls=%d fallbacks=%d (wanted %s by=%s, and the counter "
                                   "of the route that ran above zero, the other zero)"
                                   % (label, g, on, by, calls, falls, "on" if want_on else "off", "/".join(want_by)))
            # the physics: steps, and fold runs a force evaluation, which
            # must be the rows' length over the widest fold, rounded up.
            # Every body active, a row is N - 1 long, N read from the
            # record's own header. With an `active` line the longest row is
            # a test particle's, N_active long (testparticle_type 0), read
            # from the record's `# active` line.
            r0 = ref[5]
            evals = sum(r0.get("pair_gather", (0, 0, 0, 0))[2:])
            folds = sum(r0.get("indexed", (0, 0, 0, 0))[2:])
            runs = folds / evals if evals else 0
            m = re.search(r"^# program=.* N=(\d+) ", ref[6], re.M)
            rows = int(m.group(1)) - 1 if m else -1
            ma = re.search(r"^# active (\d+) of (\d+) bodies a system", ref[6], re.M)
            if ma and int(ma.group(1)) < int(ma.group(2)):
                rows = int(ma.group(1))
            want = (rows + wmax - 1) // wmax if rows > 0 and wmax else -1
            if evals == 0 or runs != want:
                bad.append("the physics: %d force evaluations and %.3f fold runs each, for rows of %d and a widest "
                           "fold of %d (wanted %d)" % (evals, runs, rows, wmax, want))
            rejected = ref[4].get("steps_rejected", "")
            if name == rejecting and not (rejected.isdigit() and int(rejected) > 0):
                bad.append("the physics: steps_rejected=%s, and this case is here for a rejected step (V1, 2026-11-02: "
                           "it was printed and not asserted)" % (rejected or "missing"))
            physics = ("%s steps, %s rejected, %d force evaluations, %.3g fold runs each"
                       % (ref[4].get("steps_done"), ref[4].get("steps_rejected"), evals, runs))
            if not bad:
                verdict(True, "%s: %d sample lines, all 5 runs identical, counters and route lines as they must be; %s"
                        % (name, len(ref[3]), physics))
            else:
                verdict(False, "%s: %d of the checks failed (%s):%s" % (name, len(bad), physics,
                                                                       "".join("\n         " + b for b in bad[:8])))

    # 2b. Test particles: the active line against the same bodies all
    # active (round 4's P0). Massless bodies add exactly +0 to every sum
    # (tools/cases_pairs.c, case_massless_is_inert, for the drop-in), so
    # the two runs must agree bit for bit - and the control must not.
    print("check_gravity: test particles (`active 5`) against the same bodies all active, the same record")
    want_pairs = 5 * 4 // 2 + (8 - 5) * 5     # REBOUND's two loops: the active pairs, then test-active
    active_line = "# active 5 of 8 bodies a system: %d pairs a force evaluation, a system" % want_pairs
    for fmt in formats:
        for engine, arith in (("loop", "rebound"), ("loop", "fma"), ("program", "fma")):
            base = [cft, "--format", fmt, "--quiet", "--engine", engine, "--arith", arith,
                    "--dt", "40", "--epsilon", "0", "--steps", str(20 * s), "--sample", "10"]
            res = {}
            for key, path in (("a5", p_a5), ("all", p_all), ("a5m", p_a5m), ("allm", p_allm)):
                rc, out, err = run(base + ["--problem", path])
                res[key] = (rc, [l for l in out.splitlines() if l and not l.startswith("#")], out, err)
            for key, path in (("e5", p_ens5), ("eall", p_ensall)):
                rc, out, err = run(base + ["--problem", path])
                res[key] = (rc, [l for l in out.splitlines() if l and not l.startswith("#")], out, err)
            rc, out, err = run(base + ["--problem", p_ens5, "--member", "1"])
            res["m1"] = (rc, [l for l in out.splitlines() if l and not l.startswith("#")], out, err)
            bad = ["%s: rc %d: %s" % (k, r[0], r[3].strip()[-200:]) for k, r in res.items() if r[0] != 0]
            if not bad:
                if not res["a5"][1] or res["a5"][1] != res["all"][1]:
                    bad.append("the records differ, or are empty (%d and %d sample lines)" % (len(res["a5"][1]), len(res["all"][1])))
                if active_line not in res["a5"][2].splitlines():
                    bad.append("no '%s' line" % active_line)
                if any(l.startswith("# active") for l in res["all"][2].splitlines()):
                    bad.append("an `# active` line in the run with no active line")
                if res["a5m"][1] == res["allm"][1]:
                    bad.append("the control did not bite: with tp2 given a mass of 2^-40 the two runs still agree")
                if res["a5m"][1] != res["a5"][1]:
                    bad.append("testparticle_type 0 does not hold: under active 5 a test particle's 2^-40 mass changed "
                               "the record (an active body moved, or the energy column counted it)")
                if not res["e5"][1] or res["e5"][1] != res["eall"][1]:
                    bad.append("the two-system ensemble with active 5 differs from it all active, or is empty")
                # the member's lines are the system's own, without the ensemble's leading system field
                mine = [l.split(None, 2)[2] for l in res["e5"][1] if l.split(None, 2)[:2] == ["system", "1"]]
                if not mine or mine != res["m1"][1]:
                    bad.append("--member 1 of the active 5 ensemble differs from that system's lines in the ensemble "
                               "(%d against %d)" % (len(res["m1"][1]), len(mine)))
                cols = re.search(r"^# columns: (.*?) then per body", res["a5"][2], re.M)
                lead = len(cols.group(1).split()) if cols else -1
                first, last = res["a5"][1][0].split(), res["a5"][1][-1].split()
                tp1x = lead + 6 * 6       # the record's own leading columns, then six a body; tp1 is body 6
                if lead < 0 or len(first) <= tp1x or first[tp1x] == last[tp1x]:
                    bad.append("the physics: tp1's x did not move between the first and last sample (or no # columns line)")
            label = "%s, %s engine, %s form" % (fmt, engine, arith)
            verdict(not bad, "%s: %s" % (label, "%d sample lines identical, %d pairs a force evaluation against %d all "
                                         "active, the test particles move, a 2^-40 mass changes nothing under active 5 "
                                         "(energy included) and everything without it, and the two-system ensemble and "
                                         "its member agree"
                                         % (len(res["a5"][1]), want_pairs, 8 * 7 // 2) if not bad else "; ".join(bad[:4])))

    # 2c. The problem lines' refusals (round 4's P0).
    print("check_gravity: the problem lines' refusals")
    for path, what, needle in ((p_a9, "active 9 in a system of 8", "'active 9'"),
                               (p_a0, "active 0, no massive body", "'active 0'"),
                               (p_twice, "a second active line", "a second active line")):
        rc, out, err = run([cft, "--format", "fp64", "--problem", path, "--steps", "1", "--quiet"])
        verdict(rc != 0 and needle in err and ("1 to N = 8" in err or needle == "a second active line"),
                "%s is refused by name (rc %d): %s" % (what, rc, err.strip()[-160:]))
    rc, out, err = run([cft, "--format", "fp64", "--problem", p_zonal, "--steps", "1", "--quiet"])
    verdict(rc != 0 and "round 4's P1" in err,
            "a zonal line is refused by name until round 4's P1 builds the term (rc %d): %s" % (rc, err.strip()[-160:]))
    sched = os.path.join(scratch, "outer8_active5.sched")
    common = ["--format", "fp64", "--dt", "40", "--epsilon", "1e-9", "--steps", "4", "--sample", "2", "--quiet"]
    rc1, _, err1 = run([cft, "--problem", p_a5, "--record-iters", sched] + common)
    rc2, _, err2 = run([cft, "--problem", p_all, "--replay-iters", sched] + common)
    rc3, _, err3 = run([cft, "--problem", p_a5, "--replay-iters", sched] + common)
    verdict(rc1 == 0 and rc3 == 0 and rc2 != 0 and "different run" in err2 and " active=5" in err2,
            "a schedule recorded with active 5 replays with it (rc %d, %d) and is refused without it (rc %d)%s"
            % (rc1, rc3, rc2, "" if rc2 != 0 and " active=5" in err2 else
               ": %s" % (err2 or err1 or err3).strip()[-200:]))

    # 3. The body cap, and one past it.
    print("check_gravity: the body cap's refusal")
    cap = None
    with open(os.path.join(root, "src", "ias15_limits.h")) as f:
        m = re.search(r"^#define\s+CFT_MAX_BODIES\s+(\d+)", f.read(), re.M)
        cap = int(m.group(1)) if m else None
    if cap is None:
        verdict(False, "no #define CFT_MAX_BODIES in src/ias15_limits.h")
    else:
        past = os.path.join(scratch, "nbody%d.txt" % (cap + 1))
        if not os.path.exists(past):
            subprocess.run([py, os.path.join(here, "make_nbody.py"), str(cap + 1), past], check=True, cwd=root,
                           stdout=subprocess.DEVNULL)
        rc, out, err = run([cft, "--format", "fp64", "--problem", past, "--steps", "1", "--quiet"])
        verdict(rc != 0 and ("outside 1..%d" % cap) in err,
                "%d bodies, one past CFT_MAX_BODIES = %d, is refused by name (rc %d): %s"
                % (cap + 1, cap, rc, err.strip()[-160:]))

    # 4. The library half on the new routes, forced on: the routes are
    # opt-in, so the suite's own check_dropin leg runs today's, and this
    # is the one run of the suite that takes the new ones. It is also the
    # only state in which a library half silently on today's route (M4)
    # differs from a correct one.
    # A device without the bit refuses the routes forced on, here as in the
    # matrix above, so this gate needs one that publishes it.
    if args.library:
        print("check_gravity: the drop-in's suite with both gravity routes forced on")
        dropin = os.path.join(args.build, "check_dropin" + exe)
        env = {routes[g]["env"]: "on" for g in GRAVITY}
        rc, out, err = run([dropin, "--light"], env)
        lines = out.splitlines()
        seen = {}
        for line in lines:
            m = re.match(r"^\s+(ok  |FAIL) route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)", line)
            if m:
                seen[m.group(2)] = (m.group(1).strip(), m.group(3), m.group(4), int(m.group(5)), int(m.group(6)))
        right = all(g in seen and seen[g][:3] == ("ok", "on", "forced") and seen[g][3] > 0 and seen[g][4] == 0
                    for g in GRAVITY)
        verdict(rc == 0 and "check_dropin: every case passed" in out and right,
                "check_dropin --light, gravity forced onto the new routes: rc %d, %s; %s"
                % (rc, "every case passed" if "check_dropin: every case passed" in out else "NOT every case passed",
                   ", ".join("%s=%s by=%s calls=%d fallbacks=%d" % ((g,) + seen[g][1:]) if g in seen else g + " unread"
                             for g in GRAVITY)))
        if rc != 0:
            for line in lines:
                if "FAIL" in line:
                    print("       " + line.strip())

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
