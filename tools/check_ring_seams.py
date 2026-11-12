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
#        that theory, which holds P2's speeds to the J2 and J4 they name.
#        With the line, P1's field must keep every particle's distance from
#        Saturn within 1e-3 of the same theory's smallest excursion: a J2
#        or J4 radial force wrong by more than about 0.1% of J2's, in P1's
#        term or in P2's speeds, fails it (measured: about 3e-5).
#        THE CONTROL, permanent: the same file with J4 zeroed in its zonal
#        line. The particles' speeds still carry J4's share of delta,
#        -15/8 J4 (R/r)^4, so each must move outward by that share's
#        delta r (1 - cos M), within 10%, and so leave the bound above by
#        more than ten times: the bound resolves J4, not only J2.
#        The energy, the zonal potential included, must hold to binary64's
#        level over the run.
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
    phys = {}
    for tag, path in (("with", small), ("without", small_nz), ("noj4", small_j4)):
        rc, out, err = run([exe, "--format", "fp64", "--problem", path, "--steps", str(steps), "--sample", "1"] + common)
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
    rw, rn, r4 = radii(phys["with"]), radii(phys["without"]), radii(phys["noj4"])
    if len(rw) != steps + 1 or len(rn) != steps + 1 or len(r4) != steps + 1:
        raise SystemExit("check_ring_seams: %d, %d and %d samples, not %d" % (len(rw), len(rn), len(r4), steps + 1))
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
    e0, e1 = rw[0][2], rw[-1][2]
    verdict(abs(e1 - e0) <= 1e-13 * abs(e0),
            "the energy with the zonal potential moves %.2e of itself over %d steps at binary64 (the bound is 1e-13)"
            % (abs(e1 - e0) / abs(e0), steps))

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

    if not args.keep:
        shutil.rmtree(wd, ignore_errors=True)
    print("RESULT: %s" % ("FAIL" if FAILS else "PASS"))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
