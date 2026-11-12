#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Round 4's seam tests (ROUND4.md, "What the integrator keeps"). Each holds
# two parcels together where neither parcel's own gate looks. The lead owns
# this file; a part is added as each parcel lands.
#
# 1. THE RING IN ITS OWN FIELD (P1 x P2). P2's generator writes every file
#    circular in its own field: with its zonal line, each particle at the
#    circular speed of Saturn's J2 and J4 in the equatorial plane,
#        v^2 = (G M / r) (1 + 3/2 J2 (R/r)^2 - 15/8 J4 (R/r)^4).
#    P1's term is that field. Neither gate runs the two together:
#    check_ring runs --no-zonal rings, and check_zonal its own problem.
#    So this part runs a ring WITH its zonal line, and the same file with
#    the line deleted (P2: delete it, do not regenerate with --no-zonal,
#    or the speeds change too):
#    (a) THE PHYSICS. Without the line each particle starts too fast by
#        delta = 3/2 J2 (R/r)^2 - 15/8 J4 (R/r)^4, so it is at the
#        pericentre of an orbit of eccentricity delta, and its distance
#        from Saturn grows as delta r (1 - cos M), M the mean anomaly the
#        run reaches. That is about 1,000 km in two thirds of an orbit.
#        Each particle's measured growth must be outward and within 10% of
#        that theory, which holds P2's speeds to the J2 they name. (J4's
#        share of delta is about 1.4%, under this check's 10%: J4 is held
#        by the circular check below and its J4-zeroed control.)
#        With the line, P1's field must keep every particle's distance from
#        Saturn within 1e-3 of the same theory's smallest excursion: a J2
#        or J4 radial force wrong by more than about 0.1% of J2's, in P1's
#        term or in P2's speeds, fails it (measured: 3.5e-5 over the 12
#        steps of --quick, 1.5e-4 over the 24 of the full run).
#        THE CONTROL, permanent: the same file with J4 zeroed in its zonal
#        line. The particles' speeds still carry J4's share of delta,
#        -15/8 J4 (R/r)^4, so each must move outward by that share's
#        delta r (1 - cos M), within 10%, and so leave the bound above by
#        more than ten times: the bound resolves J4, not only J2.
#        The energy, the zonal potential included, must hold to binary64's
#        level over the run.
#        THE NODES (ROUND4.md's precession counter). The run with the line
#        goes on for about ten orbits. Each particle's node, from h = r x v
#        relative to Saturn and fitted over the run, must regress at
#        n - nu, first-order theory in J2 and J4:
#          n^2  = (G M / r^3) (1 + 3/2 J2 (R/r)^2 - 15/8 J4 (R/r)^4)
#          nu^2 = (G M / r^3) (1 + 9/2 J2 (R/r)^2 - 75/8 J4 (R/r)^4),
#        to 0.5% (measured: under 0.1%, the size of J2^2), and the same
#        theory without J4 must miss by at least five times as much: the
#        run resolves J4's vertical force, which the circular check above,
#        in the plane, does not reach. In --quick, three orbits, the bound
#        is 2% and J4 is not claimed.
#    (b) THE BITS. With the zonal line: the program engine against the loop
#        engine in the FMA form; `active 7` against every body active (the
#        test particles' reaction is not taken under `active`, and is an
#        exact zero for a massless body when every body is active, so the
#        samples are identical while the record's counts differ: 420 pairs
#        against none, a reaction of 6 against 63, each derived from the
#        file); a two-wave ensemble against each wave alone (--member).
#    THE RECORD. The `# zonal` line must name R, J2 and J4 as the file
#    gives them, and its counts must be the file's: N - 1 bodies feel the
#    term, and active - 1 give body 0 their reaction.
#
# 2. A RING PARTICLE REPLAYED OUT OF ITS WAVE, ITS ZONAL LINE RIDING ALONG
#    (P2 x P3, P1 x P3). P3's gate replays particles out of waves it makes
#    by hand, with no zonal line, since P1 had not landed. Here the waves
#    are P2's generator's, 64 bodies with the moons at their own masses,
#    and they carry P1's zonal line. The sub-system keeps that line
#    (extract_particle copies it verbatim, and the run key names it), so
#    P1's term is evaluated over 8 bodies in the replay and 64 in the wave.
#    The term's reaction sums only the active bodies, ascending, and a
#    massless body's acceleration reads only the active ones, so the
#    replay must give the wave's bits: the particle, all seven massive
#    bodies, every sample's head and exact-time pair, and the counters.
#    Three particles a wave (its innermost, middle and outermost), from a
#    sane first step and from one far too large, whose first attempt must
#    be rejected: the one kind of rejection that leaves state behind (P3,
#    measured). P3's extract() and its gate's compare() are imported, not
#    copied.
#    THE CONTROL, permanent: one sub-system with its zonal line deleted,
#    its schedule re-keyed by hand through the program's own --run-key (the
#    tool would refuse it), must replay to other bits in the massive
#    bodies. So the line is load-bearing in the replay, and a tool that
#    dropped it would fail this part.
#
# 3. THE DRIVER ON THE RING'S OWN SHAPE (P2 x P4, P1 x P4; VL, 2026-11-12:
#    no gate ran it). tools/compare_run.py's own gate drives hand-made
#    waves with no zonal line, and a GPU run's record will be P2's waves
#    with P1's line. So a two-wave ensemble of P2's ring, its zonal line
#    kept, is recorded from a first step far too large (each wave rejects
#    its first attempt, so both are sampled), and the driver holds it: both
#    members IDENTICAL, each one's two particles replayed alone with the
#    zonal line and IDENTICAL to the record. THE CONTROL, permanent: one
#    hex digit of a moon's x on wave 0's last sample line. Wave 0 FAILs, and
#    both its particles' replays FAIL naming that line and column, while
#    wave 1 and its particles stay IDENTICAL.
#
#   python tools/check_ring_seams.py [--build build] [--quick] [--keep]

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_ensemble import parse  # noqa: E402

root = os.path.dirname(here)
PY = sys.executable

FAILS = []


def verdict(ok, text):
    print("  %s %s" % ("ok  " if ok else "FAIL", text))
    if not ok:
        FAILS.append(text)
    return ok


def run(cmd):
    p = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def columns(text):
    """The record's own `# columns:` line: how many head values precede the
    bodies, and the six names of each body's block, asserted."""
    m = re.search(r"^# columns: sample step (.*?) then per body (.*?) ; then exact-time pair", text, re.M)
    if not m:
        raise SystemExit("check_ring_seams: no # columns: line of the expected form")
    head = m.group(1).split()
    body = m.group(2).split()
    if body != ["x", "y", "z", "vx", "vy", "vz"] or "E" not in head:
        raise SystemExit("check_ring_seams: the columns line is %r" % m.group(0))
    return len(head), head.index("E"), head.index("t")


def problem_lines(path):
    with open(path) as f:
        return f.read().splitlines()


def field(lines, key):
    for l in lines:
        if l.split()[:1] == [key]:
            return l.split()[1:]
    return None


def hexf(s):
    return float.fromhex(s)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="build")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()
    exe = os.path.join(root, args.build, "ias15_cft.exe" if os.name == "nt" else "ias15_cft")
    if not os.path.exists(exe):
        exe = os.path.join(root, args.build, "ias15_cft")
    wd = os.path.join(root, args.build, "check_ring_seams")
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    print("check_ring_seams: %s, built %s" % (os.path.relpath(exe, root),
          time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(exe)))))

    # ---------------------------------------------------------------- 1
    print("check_ring_seams 1: the ring in its own field (P1 x P2)")
    small = os.path.join(wd, "ring15_z.txt")
    waves = os.path.join(wd, "ring2w_z.txt")
    for out, extra in ((small, ["--one-system", "--particles", "8"]), (waves, ["--particles", "114"])):
        rc, o, e = run([PY, "-B", os.path.join(here, "make_ring.py"), out, "--seed", "1"] + extra)
        if rc != 0:
            raise SystemExit("check_ring_seams: make_ring.py failed: %s" % (e or o).strip()[-300:])
    def without_zonal(src, dst):
        lines = problem_lines(src)
        kept = [l for l in lines if l.split()[:1] != ["zonal"]]
        if len(kept) != len(lines) - 1:
            raise SystemExit("check_ring_seams: %s has %d zonal lines, not 1" % (src, len(lines) - len(kept)))
        with open(dst, "w", newline="\n") as f:
            f.write("\n".join(kept) + "\n")
        return dst
    small_nz = without_zonal(small, os.path.join(wd, "ring15_nz.txt"))
    def j4_zeroed(src, dst):
        lines = problem_lines(src)
        out = []
        for l in lines:
            f = l.split()
            if f[:1] == ["zonal"]:
                l = " ".join(f[:4] + ["0x0.0p+0"])
            out.append(l)
        with open(dst, "w", newline="\n") as fh:
            fh.write("\n".join(out) + "\n")
        return dst
    small_j4 = j4_zeroed(small, os.path.join(wd, "ring15_noj4.txt"))
    def all_active(src, dst):
        lines = problem_lines(src)
        kept = [l for l in lines if l.split()[:1] != ["active"]]
        with open(dst, "w", newline="\n") as f:
            f.write("\n".join(kept) + "\n")
        return dst

    pl = problem_lines(small)
    N = int(field(pl, "N")[0]) if field(pl, "N") else sum(1 for l in pl if l.startswith("body "))
    n_active = int(field(pl, "active")[0])
    zk, zR, zJ2, zJ4 = field(pl, "zonal")
    zk = int(zk)
    GM = hexf(field(pl, "G")[0]) * hexf([l for l in pl if l.startswith("body ")][zk].split()[2])
    R, J2, J4 = hexf(zR), hexf(zJ2), hexf(zJ4)

    common = ["--dt", "100", "--epsilon", "1e-9", "--max-iter", "60"]
    steps = 12 if args.quick else 24
    # the run with the zonal line goes on for the nodes: about ten orbits,
    # three in quick mode; its circularity is read over the first `steps`
    psteps = 120 if args.quick else 360
    phys = {}
    for tag, path, ns in (("with", small, psteps), ("without", small_nz, steps), ("noj4", small_j4, steps)):
        rc, out, err = run([exe, "--format", "fp64", "--problem", path, "--steps", str(ns), "--sample", "1"] + common)
        if rc != 0:
            raise SystemExit("check_ring_seams: the %s run failed (rc %d): %s" % (tag, rc, err.strip()[-300:]))
        phys[tag] = out
    nhead, iE, it = columns(phys["with"])

    # the record's zonal line, against the file
    zl = re.search(r"^# zonal body (\d+) of every system, about z: R= (\S+) J2= (\S+) J4= (\S+); it acts on (\d+) "
                   r"bodies a system, and body (\d+) takes the reaction of (\d+)$", phys["with"], re.M)
    verdict(zl is not None and int(zl.group(1)) == zk and hexf(zl.group(2)) == R and hexf(zl.group(3)) == J2
            and hexf(zl.group(4)) == J4 and int(zl.group(5)) == N - 1 and int(zl.group(7)) == n_active - 1,
            "the record's zonal line is the file's: body %d, R %s, J2 %s, J4 %s, felt by N - 1 = %d bodies and "
            "body %d's reaction from active - 1 = %d%s"
            % (zk, zR, zJ2, zJ4, N - 1, zk, n_active - 1, "" if zl else " (no # zonal line)"))
    verdict(not re.search(r"^# zonal ", phys["without"], re.M),
            "the file with its zonal line deleted runs without the term (no # zonal line)")

    # (a) the physics
    def radii(text):
        out = []
        for _, step, v in parse(text)[0]:
            t = float(v[it][0])
            xs = [float(v[nhead + 6 * b + c][0]) for b in range(N) for c in range(3)]
            sx = xs[3 * zk:3 * zk + 3]
            r = [math.dist(xs[3 * b:3 * b + 3], sx) for b in range(N)]
            out.append((t, r, float(v[iE][0])))
        return out
    rw_all, rn, r4 = radii(phys["with"]), radii(phys["without"]), radii(phys["noj4"])
    if len(rw_all) != psteps + 1 or len(rn) != steps + 1 or len(r4) != steps + 1:
        raise SystemExit("check_ring_seams: %d, %d and %d samples, not %d, %d and %d"
                         % (len(rw_all), len(rn), len(r4), psteps + 1, steps + 1, steps + 1))
    rw = rw_all[:steps + 1]
    BOUND = 1e-3
    def excursion(rs, share):
        """Each particle's growth against share(q) r (1 - cos M): whether all
        moved outward, the worst relative deviation, and the theory's range."""
        tend = rs[-1][0]
        out, worst, lo, hi = True, 0.0, float("inf"), 0.0
        for b in range(n_active, N):
            r0 = rs[0][1][b]
            q = (R / r0) ** 2
            M = 2 * math.pi * tend / (2 * math.pi * math.sqrt(r0 ** 3 / GM))
            theory = share(q) * r0 * (1 - math.cos(M))
            grew = rs[-1][1][b] - r0
            out &= grew > 0
            worst = max(worst, abs(grew - theory) / theory)
            lo, hi = min(lo, theory), max(hi, theory)
        return out, worst, lo, hi, tend
    delta = lambda q: 1.5 * J2 * q - 15.0 / 8.0 * J4 * q * q
    j4share = lambda q: -15.0 / 8.0 * J4 * q * q
    out_n, dev_n, lo_n, hi_n, tend = excursion(rn, delta)
    verdict(out_n and dev_n < 0.10,
            "without the zonal line each of the %d particles moves outward as an orbit of eccentricity delta "
            "from its pericentre: %.1f to %.1f km by t = %.0f s, within %.1f%% of delta r (1 - cos M) (the bound is 10%%)"
            % (N - n_active, lo_n, hi_n, tend, 100 * dev_n))
    worst_with = max(max(abs(s[1][b] - rw[0][1][b]) for s in rw) for b in range(n_active, N))
    verdict(worst_with < BOUND * lo_n,
            "with it, P1's field keeps every particle circular: its distance from Saturn moves %.4f km at most, "
            "%.2e of the smallest excursion the field's absence gives (the bound is %.0e)"
            % (worst_with, worst_with / lo_n, BOUND))
    out_4, dev_4, lo_4, hi_4, _ = excursion(r4, j4share)
    worst_4 = max(max(abs(s[1][b] - r4[0][1][b]) for s in r4) for b in range(n_active, N))
    verdict(out_4 and dev_4 < 0.10 and worst_4 > 10 * BOUND * lo_n,
            "the control, J4 zeroed in the zonal line: each particle moves outward by J4's share of delta, %.2f to "
            "%.2f km, within %.1f%% of -15/8 J4 (R/r)^4 r (1 - cos M), and %.1f times the bound above, so the bound "
            "resolves J4" % (lo_4, hi_4, 100 * dev_4, worst_4 / (BOUND * lo_n)))
    e0, e1 = rw_all[0][2], rw_all[-1][2]
    verdict(abs(e1 - e0) <= 1e-13 * abs(e0),
            "the energy with the zonal potential moves %.2e of itself over %d steps at binary64 (the bound is 1e-13)"
            % (abs(e1 - e0) / abs(e0), psteps))

    # the nodes: each particle's node, from h = r x v relative to Saturn,
    # fitted by least squares over the run, against n - nu
    def node_rates(text):
        series = {b: ([], []) for b in range(n_active, N)}
        prev = {}
        for _, step, v in parse(text)[0]:
            t = float(v[it][0])
            s = [float(v[nhead + 6 * zk + c][0]) for c in range(6)]
            for b in range(n_active, N):
                p = [float(v[nhead + 6 * b + c][0]) for c in range(6)]
                r = [p[i] - s[i] for i in range(3)]
                u = [p[3 + i] - s[3 + i] for i in range(3)]
                h = (r[1] * u[2] - r[2] * u[1], r[2] * u[0] - r[0] * u[2], r[0] * u[1] - r[1] * u[0])
                om = math.atan2(h[0], -h[1])
                if b in prev:
                    om += 2 * math.pi * round((prev[b] - om) / (2 * math.pi))
                prev[b] = om
                series[b][0].append(t)
                series[b][1].append(om)
        out = {}
        for b, (ts, oms) in series.items():
            mt, mo = sum(ts) / len(ts), sum(oms) / len(oms)
            out[b] = sum((a - mt) * (c - mo) for a, c in zip(ts, oms)) / sum((a - mt) ** 2 for a in ts)
        return out, ts[-1] - ts[0]
    rates, span = node_rates(phys["with"])
    worst, least_nj4 = 0.0, float("inf")
    for b, rate in rates.items():
        r0 = sum(s[1][b] for s in rw_all) / len(rw_all)
        q, n0 = (R / r0) ** 2, math.sqrt(GM / r0 ** 3)
        theory = n0 * (math.sqrt(1 + 1.5 * J2 * q - 15 / 8 * J4 * q * q) - math.sqrt(1 + 4.5 * J2 * q - 75 / 8 * J4 * q * q))
        nj4 = n0 * (math.sqrt(1 + 1.5 * J2 * q) - math.sqrt(1 + 4.5 * J2 * q))
        worst = max(worst, abs(rate / theory - 1))
        least_nj4 = min(least_nj4, abs(rate / nj4 - 1))
    r_in = min(rw_all[0][1][b] for b in range(n_active, N))
    orbits = span * math.sqrt(GM / r_in ** 3) / (2 * math.pi)
    NODE_BOUND = 0.02 if args.quick else 0.005
    resolved = least_nj4 >= 5 * worst
    verdict(rates and all(r < 0 for r in rates.values()) and worst < NODE_BOUND and (args.quick or resolved),
            "the nodes: every particle's node regresses at n - nu, first-order theory in J2 and J4, to %.2e at worst "
            "over %d steps (%.1f orbits of the innermost; the bound is %g); the same theory without J4 misses by %.2e "
            "at least, %.0f times as far%s"
            % (worst, psteps, orbits, NODE_BOUND, least_nj4, least_nj4 / worst if worst else float("inf"),
               "" if not args.quick else " (in --quick, three orbits, J4 is not claimed resolved)"))

    # (b) the bits
    fmts = ["fp64"] if args.quick else ["fp64", "fp128", "fp256"]
    for fmt in fmts:
        ns = {"fp64": 6, "fp128": 4, "fp256": 2}[fmt]
        base = [exe, "--format", fmt, "--problem", small, "--steps", str(ns), "--sample", "1", "--arith", "fma"] + common
        r1, o1, e1_ = run(base + ["--engine", "loop"])
        r2, o2, e2_ = run(base + ["--engine", "program"])
        s1, t1 = parse(o1)
        s2, t2 = parse(o2)
        same = r1 == 0 and r2 == 0 and s1 and s1 == s2 and all(t1.get(k) == t2.get(k) for k in
                                                                  ("steps_done", "steps_rejected", "mean_pc_iterations"))
        verdict(same, "%s, with the zonal line: the program engine against the loop engine in the FMA form, %d sample "
                "lines identical, the counters too (%s steps, %s passes a step)"
                % (fmt, len(s1), t1.get("steps_done"), t1.get("mean_pc_iterations")))
    wave1 = os.path.join(wd, "wave1_z.txt")
    wl = problem_lines(waves)
    head = [l for l in wl if not l.startswith(("system ", "body ", "E "))]
    first = []
    insys = -1
    for l in wl:
        if l.startswith("system "):
            insys += 1
        elif l.startswith("body ") and insys == 0:
            first.append(l)
    with open(wave1, "w", newline="\n") as f:
        f.write("\n".join(head + first) + "\n")
    wave1_all = all_active(wave1, os.path.join(wd, "wave1_z_all.txt"))
    wn = int(field(problem_lines(wave1), "N")[0])
    wa = int(field(problem_lines(wave1), "active")[0])
    ra, oa, ea = run([exe, "--format", "fp64", "--problem", wave1, "--steps", "4", "--sample", "1"] + common)
    rb, ob, eb = run([exe, "--format", "fp64", "--problem", wave1_all, "--steps", "4", "--sample", "1"] + common)
    sa, ta = parse(oa)
    sb, tb = parse(ob)
    pairs_a = re.search(r"^# active (\d+) of (\d+) bodies a system: (\d+) pairs", oa, re.M)
    react_a = re.search(r"takes the reaction of (\d+)$", oa, re.M)
    react_b = re.search(r"takes the reaction of (\d+)$", ob, re.M)
    want_pairs = wa * (wa - 1) // 2 + (wn - wa) * wa
    ok = (ra == 0 and rb == 0 and sa and sa == sb and pairs_a and int(pairs_a.group(3)) == want_pairs
          and not re.search(r"^# active ", ob, re.M) and react_a and react_b
          and int(react_a.group(1)) == wa - 1 and int(react_b.group(1)) == wn - 1)
    verdict(ok, "fp64, a wave with its zonal line: `active %d` against all %d bodies active, %d sample lines "
            "identical; %s pairs a force evaluation against none named (n(n-1)/2 + (N-n)n = %d from the file), "
            "and a reaction of %s against %s (active - 1 = %d, N - 1 = %d)"
            % (wa, wn, len(sa), pairs_a.group(3) if pairs_a else "no", want_pairs,
               react_a.group(1) if react_a else "?", react_b.group(1) if react_b else "?", wa - 1, wn - 1))
    re_, oe, ee = run([exe, "--format", "fp64", "--problem", waves, "--steps", "4", "--sample", "1"] + common)
    for k in (0, 1):
        rm, om, em = run([exe, "--format", "fp64", "--problem", waves, "--member", str(k), "--steps", "4",
                          "--sample", "1"] + common)
        es, et = parse(oe, k)
        ms, mt = parse(om)
        verdict(re_ == 0 and rm == 0 and es and es == ms and all(et.get(x) == mt.get(x) for x in
                                                                  ("steps_done", "steps_rejected", "mean_pc_iterations")),
                "fp64, the two-wave ensemble with its zonal line: wave %d's %d sample lines and counters are its run "
                "alone (--member %d)" % (k, len(es), k))

    part2(exe, wd, waves, args.quick)
    part3(exe, wd, waves)

    if not args.keep:
        shutil.rmtree(wd, ignore_errors=True)
    print("RESULT: %s" % ("FAIL" if FAILS else "PASS"))
    return 1 if FAILS else 0


def part3(exe, wd, waves):
    print("check_ring_seams 3: the driver on the ring's own shape, P2's waves with P1's zonal line (P2 x P4, P1 x P4)")
    pl = problem_lines(waves)
    n_active = int(field(pl, "active")[0])
    E = int(field(pl, "E")[0])
    rc, rec_text, err = run([exe, "--format", "fp64", "--problem", waves, "--dt", "1e5", "--epsilon", "1e-9",
                             "--max-iter", "60", "--steps", "8", "--sample", "4", "--quiet"])
    if rc != 0:
        verdict(False, "the ring's ensemble did not run (rc %d): %s" % (rc, err.strip()[-200:]))
        return
    rec = os.path.join(wd, "ring2w_z.rec")
    with open(rec, "w", newline="\n") as f:
        f.write(rec_text)
    rejected = [int(m.group(1)) for m in re.finditer(r"^# system \d+ steps_done=\S+ .*?steps_rejected=(\d+)",
                                                     rec_text, re.M)]
    drv = [PY, "-B", os.path.join(here, "compare_run.py"), None, "--problem", waves, "--format", "fp64",
           "--arith", "rebound", "--exe", exe, "--strata", str(E), "--seed", "1"]

    def drive(path):
        cmd = list(drv)
        cmd[3] = path
        return run(cmd)

    def members(out):
        return {int(m.group(1)): m.group(2) for m in re.finditer(r"^member (\d+) \(\S+\): (IDENTICAL|FAIL): ", out, re.M)}

    def particles(out):
        return [(int(m.group(1)), m.group(2), m.group(3), m.group(4)) for m in re.finditer(
            r"^  member (\d+): particle \S+ \(body (\d+), [^)]*\): replayed alone with the member's \d+ massive bodies"
            r"( and its zonal line)? .*?held to the record: (IDENTICAL|FAIL)", out, re.M)]

    rc, out, err = drive(rec)
    mem, parts = members(out), particles(out)
    ok = (rc == 0 and len(rejected) == E and all(r > 0 for r in rejected) and sorted(mem) == list(range(E))
          and all(v == "IDENTICAL" for v in mem.values()) and len(parts) == 2 * E
          and all(p[2] and p[3] == "IDENTICAL" for p in parts) and "\nRESULT: PASS: " in out)
    verdict(ok, "fp64, the driver on a %d-wave ensemble of P2's ring with P1's zonal line, from a first step of 1e5 s "
            "(rejections %s, so every wave is sampled): %d members IDENTICAL, and %d particles replayed alone with the "
            "zonal line, IDENTICAL to the record (rc %d)"
            % (E, rejected, sum(v == "IDENTICAL" for v in mem.values()), sum(p[3] == "IDENTICAL" for p in parts), rc))

    # the control: one hex digit of a moon's x on wave 0's last sample line,
    # its token found by P4's own reading of the record's `# columns:` line
    from compare_run import Columns
    rows = rec_text.splitlines()
    N = int(field(pl, "N")[0])
    cols = Columns(next(l for l in rows if l.startswith("# columns:")), N)
    i = max(i for i, l in enumerate(rows) if l.startswith("system 0 sample "))
    toks = rows[i].split()
    moon = 1
    t = cols.body_tok(moon, "x")
    old = toks[t]
    mant = old.split("p")[0]
    toks[t] = mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + old[len(mant):]
    bad = os.path.join(wd, "ring2w_z_moon.rec")
    with open(bad, "w", newline="\n") as f:
        f.write("\n".join(rows[:i] + [" ".join(toks)] + rows[i + 1:]) + "\n")
    rc, out, err = drive(bad)
    mem, parts = members(out), particles(out)
    hit = re.findall(r"^  member 0: particle .*?held to the record: FAIL: record line %d \(system 0 sample \d+, step \d+\): "
                     r"x of body %d " % (i + 1, moon), out, re.M)
    ok = (rc == 1 and mem.get(0) == "FAIL" and all(mem.get(k) == "IDENTICAL" for k in range(1, E)) and len(hit) == 2
          and all(p[3] == "IDENTICAL" for p in parts if p[0] != 0)
          and re.search(r"^member 0 \(\S+\): FAIL: first differing line: record line %d .*?x of body %d " % (i + 1, moon),
                        out, re.M) is not None)
    verdict(ok, "the control: one hex digit of a moon's x (body %d) on wave 0's line %d: wave 0 FAILs naming it, and so "
            "do %d of its 2 particles' replays, while the other wave%s and %s particles stay IDENTICAL (rc %d)"
            % (moon, i + 1, len(hit), "s" if E > 2 else "", "their" if E > 2 else "its", rc))


def part2(exe, wd, waves, quick):
    from extract_particle import extract, Refused, program_key
    from check_replay import compare, record, outcomes, body_mismatches
    control_done = False
    print("check_ring_seams 2: a ring particle replayed out of its wave, its zonal line riding along (P2 x P3, P1 x P3)")
    head = problem_lines(waves)
    zline = [l for l in head if l.split()[:1] == ["zonal"]]
    n_active = int(field(head, "active")[0])
    N = int(field(head, "N")[0])
    picks = [n_active, (n_active + N - 1) // 2, N - 1]
    cases = [("fp64", "100", 8, 0), ("fp64", "1e5", 8, 1)]
    if not quick:
        cases += [("fp128", "100", 6, 1), ("fp128", "1e5", 6, 0), ("fp256", "1e5", 3, 0)]
    for fmt, dt0, steps, member in cases:
        tag = "%s-dt%s-m%d" % (fmt, dt0, member)
        att = os.path.join(wd, tag + ".att")
        sch = os.path.join(wd, tag + ".sched")
        run_args = ["--format", fmt, "--dt", dt0, "--epsilon", "1e-9", "--max-iter", "60", "--steps", str(steps),
                    "--sample", "1"]
        rc, wout, werr = run([exe, "--problem", waves, "--member", str(member)] + run_args
                             + ["--attempts-out", att, "--record-iters", sch])
        if rc != 0:
            verdict(False, "%s: wave %d did not run (rc %d): %s" % (tag, member, rc, werr.strip()[-200:]))
            continue
        wrec, wtail = record(wout)
        outc = outcomes(att)
        rejected = int(wtail.get("steps_rejected", "0"))
        first_rejected = bool(outc) and outc[0] == "r"
        if dt0 == "1e5":
            verdict(first_rejected and rejected > 0,
                    "%s: wave %d from a first step of %s s rejects its first attempt (%d rejected of %d attempts), the "
                    "rejection that leaves state behind" % (fmt, member, dt0, rejected, len(outc)))
        for j in picks:
            try:
                res = extract(exe, waves, att, sch, j, os.path.join(wd, tag), run_args, member=member)
            except Refused as e:
                verdict(False, "%s: body %d of wave %d: extract_particle refused: %s" % (fmt, j, member, e))
                continue
            with open(res["problem"]) as f:
                sub_zonal = [l for l in f.read().splitlines() if l.split()[:1] == ["zonal"]]
            rc, sout, serr = run([exe] + res["args"])
            if rc != 0:
                verdict(False, "%s: body %d of wave %d: the replay did not run (rc %d): %s"
                        % (fmt, j, member, rc, serr.strip()[-200:]))
                continue
            srec, stail = record(sout)
            nsamp, nvals, bad = compare(wrec, srec, res["bodies"])
            keys = ("steps_done", "steps_rejected", "iterations_max_exceeded", "mean_pc_iterations", "max_pc_iterations")
            tails = [k for k in keys if wtail.get(k) != stail.get(k)]
            zrec = re.search(r"it acts on (\d+) bodies a system, and body \d+ takes the reaction of (\d+)$", sout, re.M)
            ok = (not bad and not tails and sub_zonal == zline and zrec is not None
                  and int(zrec.group(1)) == n_active and int(zrec.group(2)) == n_active - 1)
            verdict(ok, "%s, first step %s s: body %d of wave %d replayed alone with the wave's %d massive bodies and "
                    "its zonal line (felt by %s bodies, %s reacting): %d samples, %d values, heads and exact-time pairs "
                    "the wave's, and its counters (steps_rejected=%s, mean_pc_iterations=%s)%s"
                    % (fmt, dt0, j, member, n_active, zrec.group(1) if zrec else "?", zrec.group(2) if zrec else "?",
                       nsamp, nvals, stail.get("steps_rejected"), stail.get("mean_pc_iterations"),
                       "" if ok else ": " + "; ".join((bad + ["trailer " + k for k in tails]
                                                        + ([] if sub_zonal == zline else ["the zonal line differs"]))[:3])))
            if ok and not control_done:
                control_done = True
                # the control: the same sub-system without its zonal line,
                # re-keyed by hand, must replay to other bits
                sub_nz = res["problem"][:-4] + "-nozonal.txt"
                with open(res["problem"]) as f:
                    kept = [l for l in f.read().splitlines() if l.split()[:1] != ["zonal"]]
                with open(sub_nz, "w", newline="\n") as f:
                    f.write("\n".join(kept) + "\n")
                a = list(res["args"])
                ri = a.index("--replay-iters")
                sched_nz = res["schedule"][:-6] + "-nozonal.sched"
                a = [sub_nz if x == res["problem"] else x for x in a[:ri]]
                key = program_key(exe, a)
                with open(res["schedule"]) as f:
                    sl = f.read().splitlines()
                with open(sched_nz, "w", newline="\n") as f:
                    f.write("\n".join([sl[0], "# run " + key] + sl[2:]) + "\n")
                rc, cout, cerr = run([exe] + a + ["--replay-iters", sched_nz])
                crec, _ = record(cout)
                nb = body_mismatches(wrec, crec, res["bodies"]) if rc == 0 else -1
                verdict(rc == 0 and nb > 0 and key != res["sub_key"],
                        "the control, %s: body %d's sub-system with its zonal line deleted, re-keyed by hand (its key "
                        "differs from the tool's, which names the line), replays (rc %d) to other bits in %d of its %d "
                        "sample-and-body pairs: the line is load-bearing"
                        % (fmt, j, rc, nb, len(wrec) * len(res["bodies"])))


if __name__ == "__main__":
    sys.exit(main())
