#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The wave kernel's gate (ROUND5.md, P2): src/wave/'s step, on the host
# (the twin, build/wave_host) or on a GPU (build/wave --device gpu), held
# to this tree's port, `ias15_cft --arith fma --format fp64`, BIT FOR BIT,
# every column, member by member, through tools/compare_run.py - the
# driver a GPU run's record is held by - with particles sampled through
# round 4's replay. Beside every verdict, the physics.
#
#   python tools/check_wave.py [--build build] [--quick] [--keep DIR]
#          [--device host|gpu] [--wave PATH] [--twin PATH] [--corpus PATH]
#   (`make check-quick` and `make check-light` pass --quick: fewer and
#    shorter runs, every part still run.)
#
# THE PARTS, each asserted, each named in its verdict:
#   1. Gate 4: the arithmetic, measured. build/wave_corpus (or, on a GPU,
#      build/wave_corpus_gpu --device gpu) runs wave_arith.h's add, sub,
#      mul, div, sqrt, fma and augmentedAddition over a corpus - random
#      operands, ties, subnormals, extremes - against libcft's software
#      binary64, every result's bits and the four flags the port aborts
#      on; on the host also the hex text both ways and the step's decimal
#      literals. Any mismatch fails, and is printed by name.
#   2. The constants: the program's --dump-constants against the port's.
#   3. J2 and J4: the port's own --zonal-selftest points (two systems of
#      eight, an oblate body 2, `active 6`, a test particle with a mass),
#      read from its output and evaluated here; every acceleration and both
#      potentials must be the port's bits. tools/check_zonal.py holds the
#      port's to cft-fp256's golden model and to mpmath.
#   4. Records, through tools/compare_run.py --rerun-engine loop, every
#      member sampled (--strata E), and every member's particle replays:
#      A. a ring from tools/make_ring.py with its zonal line, 3 waves of 64
#         (2 in --quick), from dt0 600 s;
#      B. the same ring from dt0 1e5 s, far too large: every wave must
#         reject its first attempt, before any step is accepted (the one
#         rejection that leaves state; CLAUDE.md);
#      C. the ring with one test particle of wave 1 aimed past Janus: that
#         wave must reject attempts AFTER accepted steps;
#      D. this tree's problems: kepler_e099 (rejections and a cap hit),
#         zonal_saturn, outer and pythagorean (full only), and kepler at
#         fixed steps;
#      E. one member of A run alone (--member 1), a one-system record;
#      F. A's ring with every test particle given a mass of 2^-20: under
#         `active 7` it moves nothing, so every sample line must be A's. A
#         massless particle hides a reaction or a kinetic energy taken from
#         it (MEASURED: two such mutants of the step pass A and fail F);
#         the reactions' ORDER, below half an ulp of Saturn's acceleration
#         in a ring, is part 3's (its selftest catches it).
#      Each run also writes every member's attempts record and corrector
#      schedule (--attempts-dir), which must be the port's own for that
#      member, line for line, under the port's --run-key.
#      Beside each: steps, rejections, passes and cap hits from the
#      record's trailers, and for a ring the physics the bits cannot show,
#      asserted (ring_physics(), below): every particle started on its
#      circular orbit in Saturn's field with J2 and J4, so its distance from
#      Saturn moves by under 1e-2 of what the field's absence would move it
#      (check_ring_seams' delta r (1 - cos M)), it sweeps its orbit's n t in
#      that field to 1e-4 (Kepler's n misses by 2.2e-3 or more), and the
#      energy holds to 1e-13; the worst particle is named, with the
#      closest a moon came to it at a sample.
#   5. Controls, each asserted to fail by name:
#      - contraction: the same source with its flag checks compiled out
#        passes, and with a*b + c fused (-ffp-contract=fast; --fmad=true
#        under nvcc) fails compare_run;
#      - one hex digit planted in a sampled member's record, a column on
#        two lines, fails compare_run naming the member and the column;
#      - a trailer counter changed fails it;
#      - the physics: wave 0 of A's ring with its zonal line deleted runs,
#        its sample lines are the port's bit for bit, and the physics
#        above must fail it (bit-identity is not validity).
#   6. Flags: a problem made to raise each of INVALID, DIVBYZERO, OVERFLOW
#      and UNDERFLOW stops the port (rc 3) and stops this program (rc 3),
#      naming the flag; and with --no-flag-abort both run on and report
#      the same four flags.
# On a GPU (--device gpu): the same parts on build/wave --device gpu, and
# every record also against the host twin's, line for line.
#
# Runs work in <build>/wave/<pid> unless --keep; a passing run removes it.

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
sys.path.insert(0, here)
from check_ensemble import parse as parse_samples   # noqa: E402
from check_replay import FLYBY                       # noqa: E402

EXE = ".exe" if os.name == "nt" else ""
FAILS = []
COUNTERS = ("steps_done", "steps_rejected", "iterations_max_exceeded", "mean_pc_iterations", "max_pc_iterations")


def verdict(ok, text):
    print("  %s %s" % ("ok  " if ok else "FAIL", text))
    if not ok:
        FAILS.append(text)
    return ok


def run(cmd, env=None):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env)
    return p.returncode, p.stdout, p.stderr


def built(path):
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(path)))


def write(path, text):
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write(text)


def trailers(text):
    """{system: counters} from a record's trailers (an ensemble's per-system
    lines, or a one-system record's totals as system 0), and the totals."""
    per, tot = {}, {}
    for l in text.splitlines():
        m = re.match(r"^# system (\d+) (steps_done=.*)$", l)
        if m:
            per[int(m.group(1))] = dict(t.split("=", 1) for t in m.group(2).split() if "=" in t)
        elif l.startswith("# steps_done="):
            tot = dict(t.split("=", 1) for t in l[2:].split() if "=" in t)
    if not per and tot:
        per[0] = tot
    return per, tot


def sample_lines(text):
    return [l for l in text.splitlines() if l.startswith("sample ") or l.startswith("system ")]


# ---- ring files ---------------------------------------------------------------
def make_ring(path, particles, seed):
    rc, out, err = run([sys.executable, os.path.join(here, "make_ring.py"), path, "--particles", str(particles),
                        "--seed", str(seed)])
    if rc != 0:
        raise SystemExit("check_wave: make_ring failed: %s" % (err or out))


def massive_test_particles(src, dst, n_active, mass):
    """src with every test particle of every system given `mass`: under the
    `active` line it still moves nothing (testparticle_type 0), so the record's
    sample lines must not move, while a step that took its reaction or its
    energy would. A massless particle hides both."""
    out, k = [], 0
    for l in open(src).read().splitlines():
        t = l.split()
        if t[:1] == ["system"]:
            k = 0
        elif t[:1] == ["body"]:
            if k >= n_active:
                t[2] = float(mass).hex()
                l = " ".join(t)
            k += 1
        out.append(l)
    write(dst, "\n".join(out) + "\n")


def flyby_ring(src, dst, wave):
    """src with the last test particle of system `wave` replaced by one aimed
    past Janus: tools/check_replay.py's FLYBY (km out from it along its
    radius, km/s toward it, km off its plane), in a make_ring wave. Plain
    binary64 arithmetic, each value written as exact hex."""
    lines = open(src).read().splitlines()
    sysidx, bodies = -1, {}
    for i, l in enumerate(lines):
        t = l.split()
        if t[:1] == ["system"]:
            sysidx += 1
        elif t[:1] == ["body"] and sysidx == wave:
            bodies.setdefault("list", []).append(i)
    rows = bodies["list"]
    jan = [l for l in (lines[i] for i in rows) if l.split()[1] == "janus"]
    sat = [l for l in (lines[i] for i in rows) if l.split()[1] == "saturn"]
    if len(jan) != 1 or len(sat) != 1:
        raise SystemExit("check_wave: wave %d of %s has no single janus and saturn" % (wave, src))
    jt, st = jan[0].split(), sat[0].split()
    jx = [float.fromhex(v) for v in jt[3:6]]
    jv = [float.fromhex(v) for v in jt[6:9]]
    sx = [float.fromhex(v) for v in st[3:6]]
    rel = [jx[k] - sx[k] for k in range(3)]
    r = math.sqrt(sum(c * c for c in rel))
    u = [c / r for c in rel]
    d, w, b = FLYBY
    x = [jx[0] + d * u[0], jx[1] + d * u[1], jx[2] + d * u[2] + b]
    v = [jv[k] - w * u[k] for k in range(3)]
    last = rows[-1]
    old = lines[last].split()
    lines[last] = "body %s %s %s" % (old[1], old[2], " ".join(float(q).hex() for q in x + v))
    write(dst, "\n".join(lines) + "\n")
    return old[1]


# ---- the physics beside the bits ---------------------------------------------
def ring_field(prob):
    """From a ring's problem file: the oblate body's index and GM (G times
    its mass), its zonal line's R, J2 and J4, and each system's body names."""
    def num(x):
        return float.fromhex(x) if x.lstrip("+-").startswith("0x") else float(x)
    G, zl, systems = 1.0, None, []
    for l in open(prob):
        t = l.split()
        if not t or t[0].startswith("#"):
            continue
        if t[0] == "G":
            G = num(t[1])
        elif t[0] == "zonal":
            zl = (int(t[1]), num(t[2]), num(t[3]), num(t[4]))
        elif t[0] == "system":
            systems.append([])
        elif t[0] == "body":
            if not systems:
                systems.append([])
            systems[-1].append((t[1], num(t[2])))
    if zl is None:
        raise SystemExit("check_wave: %s has no zonal line, so no field to hold a ring's physics to" % prob)
    k, R, J2, J4 = zl
    return k, G * systems[0][k][1], R, J2, J4, [[n for n, _ in s] for s in systems]


# Every particle of a make_ring ring starts on the circular orbit of the
# oblate body's field, v^2 = (G M / r)(1 + delta), delta = 3/2 J2 (R/r)^2
# - 15/8 J4 (R/r)^4 (tools/make_ring.py). So, the bits aside:
#   circ - its distance from the body moves only by what the moons do: the
#          largest change over the samples, over the largest excursion the
#          field's ABSENCE would give at a sample, delta r (1 - cos M)
#          (check_ring_seams' theory: without the field it starts at the
#          pericentre of an orbit of eccentricity delta);
#   mm   - it sweeps its circular orbit's n t about the body, n = sqrt(G M
#          / r^3 (1 + delta)): the angle swept over n t, less 1. Kepler's n,
#          without J2 and J4, misses by 2.2e-3 to 2.7e-3 on the gate's ring
#          (MEASURED); so mm resolves J2, and not J4's share (about 4e-5);
#   and the energy, the zonal potential included, over the run.
# What the bounds stand on, MEASURED on the full tier's 3-wave ring (round
# 5, P2, at 7dec985, where a fixed 1 km bound on the change failed): circ
# at most 1.43e-3 in case A, from p0000159, 679 km from Pandora at a
# sample, which moved 1.76 km (and p0000131, 1,911 km from Prometheus,
# 1.46 km; every other particle 0.26 km at most); 3.4e-4 in B, and 2.34e-3
# in C's waves but the flyby's; |mm| at most 9e-6; the energy 6.5e-15.
# The field's absence gives circ about 1 and mm about -2 delta, -9e-3:
# part 5's control runs it, and both must catch it.
# circ reads the samples, so it holds a run whose samples are at most
# GAP_MAX radians of the orbit apart: their largest theory excursion is then
# at least (1 + cos 1) / 2, 0.77, of the run's. mm tracks each particle's phase
# from sample to sample, and so holds any run in which a particle stays
# within half a turn of its circular orbit between samples.
CIRC_BOUND, MM_BOUND, E_BOUND, GAP_MAX = 1e-2, 1e-4, 1e-13, 2.0


def ring_physics(text, prob, n_active, member=None):
    """For each system of a ring's record: the largest change of any test
    particle's distance from the oblate body (km), circ and mm at their
    worst and whose, the closest a massive body came to circ's worst at a
    sample, Kepler's n's range, the least angle swept, and |dE/E|."""
    zk, GM, R, J2, J4, names = ring_field(prob)
    out = {}
    systems = sorted({int(l.split()[1]) for l in text.splitlines() if l.startswith("system ")}) or [None]
    for s in systems:
        samples, _ = parse_samples(text, s)
        if not samples:
            continue
        nb = (len(samples[0][2]) - 4 - 2) // 6
        rows = []
        for _, _, v in samples:
            xs = [float(v[4 + 6 * b + c][0]) for b in range(nb) for c in range(3)]
            rows.append((float(v[0][0]), xs, float(v[3][0])))
        nm = names[member if member is not None else (s or 0)]
        t0, tend = rows[0][0], rows[-1][0] - rows[0][0]
        res = {"dr": 0.0, "circ": 0.0, "circ_who": "", "mm": 0.0, "mm_who": "", "kepler": [math.inf, -math.inf],
               "swept": math.inf, "near": ("", math.inf), "t": tend, "gap": 0.0}
        for b in range(n_active, nb):
            r0 = None
            dr, th, ang, near = 0.0, 0.0, 0.0, ("", math.inf)
            for t, xs, _ in rows:
                p = [xs[3 * b + c] - xs[3 * zk + c] for c in range(3)]
                r = math.sqrt(sum(c * c for c in p))
                phi = math.atan2(p[1], p[0])          # about the oblate body's pole, the z axis
                if r0 is None:
                    r0 = r
                    q = (R / r0) ** 2
                    delta = 1.5 * J2 * q - 15.0 / 8.0 * J4 * q * q
                    nk = math.sqrt(GM / r0 ** 3)
                    n = nk * math.sqrt(1.0 + delta)
                else:
                    # the angle swept since the last sample: n dt, and what
                    # the sample says beyond it, within half a turn - so the
                    # samples may be any distance apart that leaves the
                    # particle within half a turn of its circular orbit
                    adv = n * (t - tp)
                    ang += adv + ((phi - php - adv + math.pi) % (2.0 * math.pi) - math.pi)
                    res["gap"] = max(res["gap"], adv)
                tp, php = t, phi
                dr = max(dr, abs(r - r0))
                th = max(th, delta * r0 * (1.0 - math.cos(nk * (t - t0))))
                for m in range(n_active):
                    if m != zk:
                        d = math.dist(xs[3 * b:3 * b + 3], xs[3 * m:3 * m + 3])
                        if d < near[1]:
                            near = (nm[m], d)
            who = nm[b] if b < len(nm) else str(b)
            circ = dr / th if th > 0 else math.inf
            mm = ang / (n * tend) - 1.0 if tend > 0 else math.inf
            kep = ang / (nk * tend) - 1.0 if tend > 0 else math.inf
            res["dr"] = max(res["dr"], dr)
            if circ >= res["circ"]:
                res["circ"], res["circ_who"], res["near"] = circ, who, near
            if abs(mm) >= abs(res["mm"]):
                res["mm"], res["mm_who"] = mm, who
            res["kepler"] = [min(res["kepler"][0], kep), max(res["kepler"][1], kep)]
            res["swept"] = min(res["swept"], ang)
        e0, e1 = rows[0][2], rows[-1][2]
        res["dE"] = abs(e1 - e0) / abs(e0) if e0 else 0.0
        out[s if s is not None else 0] = res
    return out


def physics_verdict(ph, skip=()):
    """circ, mm and the energy at their worst over a record's systems (less
    those in skip), each against its bound: (ok, text)."""
    sel = [v for s, v in sorted(ph.items()) if s not in skip]
    w = max(sel, key=lambda v: v["circ"])
    m = max(sel, key=lambda v: abs(v["mm"]))
    de = max(v["dE"] for v in sel)
    klo, khi = min(v["kepler"][0] for v in sel), max(v["kepler"][1] for v in sel)
    swept = min(v["swept"] for v in sel)
    gap = max(v["gap"] for v in sel)
    ok = gap <= GAP_MAX and w["circ"] <= CIRC_BOUND and abs(m["mm"]) <= MM_BOUND and de <= E_BOUND and swept > 0.0
    if gap <= GAP_MAX:
        circ = ("every test particle's distance from Saturn moved %.4f km at most, %.2e of the largest excursion J2 "
                "and J4's absence would give, delta r (1 - cos M) (the bound is %.0e; the worst, %s, %.0f km from %s "
                "at a sample)" % (max(v["dr"] for v in sel), w["circ"], CIRC_BOUND, w["circ_who"], w["near"][1],
                                  w["near"][0]))
    else:
        circ = ("the distance from Saturn NOT HELD: the samples are %.2f rad of an orbit apart, more than the %.0f "
                "it reads (it moved %.4f km at most at the samples)" % (gap, GAP_MAX, max(v["dr"] for v in sel)))
    text = ("%s; each swept its circular orbit's n t in that field to %.1e (%s; the bound is %.0e; Kepler's n "
            "misses by %.1e to %.1e), at least %.3f rad; the energy moved %.2e of itself (the bound is %.0e)"
            % (circ, m["mm"], m["mm_who"], MM_BOUND, klo, khi, swept, de, E_BOUND))
    return ok, text


# ---- compare_run ----------------------------------------------------------------
def compare(record, problem, cft, E, keep=None):
    cmd = [sys.executable, os.path.join(here, "compare_run.py"), record, "--problem", problem, "--format", "fp64",
           "--arith", "fma", "--exe", cft, "--rerun-engine", "loop", "--strata", str(max(1, E))]
    if keep:
        cmd += ["--keep", keep]
    rc, out, err = run(cmd)
    return rc, out + err


def compare_summary(out):
    members = re.findall(r"^member (\d+) \(([^)]*)\): (IDENTICAL|FAIL)", out, re.M)
    parts = re.findall(r"^  member \d+: particle (\S+) .*?: (IDENTICAL|FAIL)", out, re.M)
    res = re.search(r"^RESULT: (\w+)", out, re.M)
    return members, parts, res.group(1) if res else "none"


# ---- the gate ---------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--keep", default=None)
    ap.add_argument("--device", default="host", choices=("host", "gpu"))
    ap.add_argument("--wave", default=None, help="the program under test (default build/wave_host, or build/wave on a GPU)")
    ap.add_argument("--twin", default=None, help="on a GPU: the host twin each GPU record is also held to "
                                                 "(default: the GPU build's own --device host)")
    ap.add_argument("--corpus", default=None)
    ap.add_argument("--expect", default=None, help="on a GPU: text its device line must contain, or nothing runs")
    args = ap.parse_args()
    B = os.path.abspath(args.build)
    gpu = args.device == "gpu"
    wave = os.path.abspath(args.wave or os.path.join(B, ("wave" if gpu else "wave_host") + EXE))
    unchecked = os.path.join(B, ("wave_unchecked" if gpu else "wave_host_unchecked") + EXE)
    contract = os.path.join(B, ("wave_contract" if gpu else "wave_host_contract") + EXE)
    corpus = os.path.abspath(args.corpus or os.path.join(B, ("wave_corpus_gpu" if gpu else "wave_corpus") + EXE))
    cft = os.path.join(B, "ias15_cft" + EXE)
    dev = ["--device", "gpu"] if gpu else []
    work = os.path.abspath(args.keep or os.path.join(B, "wave", str(os.getpid())))
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)
    print("check_wave: the wave kernel on the %s, held to this tree's port" % ("GPU" if gpu else "host (the twin)"))
    for p in (wave, unchecked, contract, corpus, cft):
        if not os.path.isfile(p):
            raise SystemExit("check_wave: %s is not built (make wave-host%s)" % (p, "; make HIP=1 wave" if gpu else ""))
        print("  %s, built %s" % (p, built(p)))
    twin = None
    if gpu:
        twin = [os.path.abspath(args.twin)] if args.twin else [wave, "--device", "host"]
        # The proof the GPU ran is its own description, not a backend= the
        # program writes (V0, VALIDATION entry 54): every GPU record carries
        # a `# device` line from the runtime, and --expect names the GPU it
        # must be, so a shell that would open another device runs nothing.
        rc, out, err = run([wave] + dev + ["--problem", os.path.join(root, "data", "problems", "kepler.txt"), "--dt",
                                           "0.01", "--steps", "2", "--sample", "2", "--quiet"])
        dline = next((l for l in out.splitlines() if l.startswith("# device ")), "")
        print("  " + (dline or "no device line: %s" % (err or out).strip()[-200:]))
        if rc != 0 or not dline or (args.expect and args.expect not in dline):
            raise SystemExit("check_wave: REFUSED: the GPU is not the one expected (--expect %r): rc %d, %r"
                             % (args.expect, rc, dline))

    # 1. gate 4 ------------------------------------------------------------
    print("1. gate 4: the arithmetic against libcft's software binary64")
    rc, out, err = run([corpus] + dev + ["--n", "4000" if args.quick else "20000"])
    for l in out.splitlines():
        if l.startswith(("device:", "wave_corpus:", "text:", "RESULT")) or re.match(r"^\S+ +\d+ cases", l):
            print("    " + l)
    mism = [l for l in out.splitlines() if l.startswith("MISMATCH")]
    for l in mism[:50]:
        print("    " + l)
    verdict(rc == 0 and "RESULT: PASS" in out and not mism,
            "every result's bits and every flag of %d operations are libcft's (rc %d, %d mismatches)"
            % (sum(int(n) for n in re.findall(r"^\S+ +(\d+) cases", out, re.M)), rc, len(mism)))

    # 2. the constants ------------------------------------------------------
    print("2. the constants")
    rc1, o1, _ = run([cft, "--format", "fp64", "--dump-constants"])
    rc2, o2, _ = run([wave, "--dump-constants"])
    c1 = [l for l in o1.splitlines() if l and not l.startswith("constants:")]
    c2 = [l for l in o2.splitlines() if l]
    verdict(rc1 == rc2 == 0 and c1 == c2 and len(c1) > 70,
            "--dump-constants gives the port's %d lines (h, rr, c, d and pc_tol at binary64)" % len(c1))

    # 3. J2 and J4 ----------------------------------------------------------
    print("3. J2 and J4 at the port's --zonal-selftest points")
    rc1, o1, _ = run([cft, "--format", "fp64", "--zonal-selftest"])
    zst = os.path.join(work, "zonal_selftest.txt")
    write(zst, o1)
    rc2, o2, e2 = run([wave] + dev + ["--zonal-selftest-from", zst])
    z1 = [l for l in o1.splitlines() if l.startswith(("zst acc", "zst pot"))]
    z2 = [l for l in o2.splitlines() if l.startswith(("zst acc", "zst pot"))]
    verdict(rc1 == rc2 == 0 and z1 == z2 and len(z1) == 18,
            "%d accelerations and potentials, every one the port's bits (16 bodies of 2 systems, k = 2, active 6, "
            "a test particle with a mass)%s%s" % (len(z1), " in a kernel, a workgroup a system" if gpu else "",
                                                 "" if rc2 == 0 else ": " + e2.strip()[-200:]))

    # 4. records ------------------------------------------------------------
    print("4. records, through tools/compare_run.py --rerun-engine loop, every member sampled")
    nwaves = 2 if args.quick else 3
    steps, sample = (12, 4) if args.quick else (24, 6)
    ring = os.path.join(work, "ring.txt")
    make_ring(ring, 57 * nwaves, 1)
    fly = os.path.join(work, "ring_flyby.txt")
    flyname = flyby_ring(ring, fly, 1)
    heavy = os.path.join(work, "ring_massive_tp.txt")
    massive_test_particles(ring, heavy, 7, 2.0 ** -20)
    probs = os.path.join(root, "data", "problems")
    cases = [
        ("A", "the ring, %d waves with the zonal line, from dt0 600 s" % nwaves, ring,
         ["--dt", "600", "--steps", str(steps), "--sample", str(sample)], nwaves, 7),
        ("B", "the same ring from dt0 1e5 s, far too large", ring,
         ["--dt", "1e5", "--steps", str(steps // 2), "--sample", str(sample // 2)], nwaves, 7),
        ("C", "the ring with %s of wave 1 aimed past Janus, from dt0 100 s" % flyname, fly,
         ["--dt", "100", "--steps", str(steps if args.quick else 30), "--sample", "1"], nwaves, 7),
        ("D1", "kepler_e099 from dt0 0.01", os.path.join(probs, "kepler_e099.txt"),
         ["--dt", "0.01", "--steps", "80", "--sample", "20"], 1, None),
        ("D2", "zonal_saturn from dt0 100 s", os.path.join(probs, "zonal_saturn.txt"),
         ["--dt", "100", "--steps", "60", "--sample", "20"], 1, 2),
        ("D3", "kepler at a fixed step (epsilon 0)", os.path.join(probs, "kepler.txt"),
         ["--dt", "0.01", "--epsilon", "0", "--steps", "30", "--sample", "10"], 1, None),
        ("E", "member 1 of the ring alone (--member 1)", ring,
         ["--member", "1", "--dt", "600", "--steps", str(steps), "--sample", str(sample)], 1, 7),
        ("F", "the ring with every test particle given a mass of 2^-20, from dt0 600 s", heavy,
         ["--dt", "600", "--steps", str(steps), "--sample", str(sample)], nwaves, 7),
    ]
    if not args.quick:
        cases += [("D4", "outer from dt0 10", os.path.join(probs, "outer.txt"), ["--dt", "10", "--steps", "40", "--sample", "10"], 1, None),
                  ("D5", "pythagorean from dt0 0.01", os.path.join(probs, "pythagorean.txt"),
                   ["--dt", "0.01", "--steps", "200", "--sample", "50"], 1, None)]
    records = {}
    for name, what, prob, a, E, n_active in cases:
        rec = os.path.join(work, "%s.rec" % name)
        attd = os.path.join(work, "att_%s" % name)
        os.makedirs(attd, exist_ok=True)
        rc, out, err = run([wave] + dev + ["--problem", prob, "--quiet", "--attempts-dir", attd] + a)
        write(rec, out)
        records[name] = (rec, prob, a, E, out)
        if not verdict(rc == 0, "%s %s: the program ran (rc %d)%s" % (name, what, rc, "" if rc == 0 else ": " + err.strip()[-300:])):
            continue
        if gpu:
            d = [l for l in out.splitlines() if l.startswith("# device ")]
            verdict(len(d) == 1 and d[0] == dline, "%s: the record's one device line is the GPU's" % name)
        rc, cout = compare(rec, prob, cft, E)
        members, parts, res = compare_summary(cout)
        ident = sum(1 for m in members if m[2] == "IDENTICAL")
        pid = sum(1 for p in parts if p[1] == "IDENTICAL")
        per, tot = trailers(out)
        cnt = "; ".join("member %d: steps %s, rejected %s, mean passes %s, max %s, cap hits %s" % (
            s, c.get("steps_done"), c.get("steps_rejected"), c.get("mean_pc_iterations"), c.get("max_pc_iterations"),
            c.get("iterations_max_exceeded")) for s, c in sorted(per.items()))
        ok = rc == 0 and res == "PASS" and ident == E and pid == len(parts)
        verdict(ok, "%s %s: compare_run %s: %d of %d members IDENTICAL%s; %s" % (
            name, what, res, ident, E, ", %d of %d particles replayed IDENTICAL" % (pid, len(parts)) if parts else "", cnt))
        if not ok:
            for l in cout.splitlines():
                if "FAIL" in l or "REFUSED" in l:
                    print("    " + l)
        if gpu:
            r2, o2, e2 = run(twin + ["--problem", prob, "--quiet"] + a)
            verdict(r2 == 0 and sample_lines(o2) == sample_lines(out) and
                    {s: {k: c.get(k) for k in COUNTERS} for s, c in trailers(o2)[0].items()} ==
                    {s: {k: c.get(k) for k in COUNTERS} for s, c in per.items()},
                    "%s: the GPU's %d sample lines and counters are the host twin's, line for line"
                    % (name, len(sample_lines(out))))
        # the attempts records and schedules: the port's own, under its key
        same, keys, total = 0, 0, 0
        first_rej, late_rej = 0, 0
        for s in range(E):
            k = s if "--member" not in a else int(a[a.index("--member") + 1])
            att, sch = os.path.join(attd, "member%d.att" % k), os.path.join(attd, "member%d.sched" % k)
            pa, ps = os.path.join(work, "port_%s_%d.att" % (name, k)), os.path.join(work, "port_%s_%d.sched" % (name, k))
            base = list(a)
            if "--member" in a:
                i = a.index("--member")
                base = a[:i] + a[i + 2:]
            # the member run alone, as compare_run reruns it and as P3's
            # driver keys a producer's files: --member k; a problem run
            # whole as one system is rerun whole, with no --member (the port
            # refuses --member on a file with no system lines)
            pargs = [cft, "--format", "fp64", "--arith", "fma", "--engine", "loop", "--problem", prob, "--quiet"] +                     (["--member", str(k)] if E > 1 or "--member" in a else []) + base
            r1, _, e1 = run(pargs + ["--attempts-out", pa, "--record-iters", ps])
            r2, kout, _ = run(pargs + ["--run-key"])
            total += 1
            try:
                mine_a = open(att).read().splitlines()
                mine_s = open(sch).read().splitlines()
                port_a = [l.rstrip("\r") for l in open(pa).read().splitlines()]
                port_s = [l.rstrip("\r") for l in open(ps).read().splitlines()]
            except OSError as e:
                print("    %s member %d: %s" % (name, k, e))
                continue
            same += mine_a == port_a and mine_s == port_s
            keys += r2 == 0 and mine_a[1:2] == ["# run " + kout.strip()]
            outcomes = [l.split()[1] for l in mine_a if l and not l.startswith("#")][:-1]
            first_rej += bool(outcomes) and outcomes[0] == "r"
            late_rej += any(o == "r" for o in outcomes[outcomes.index("a"):]) if "a" in outcomes else 0
        verdict(same == total and keys == total,
                "%s: every member's attempts record and corrector schedule are the port's own, line for line (%d of %d), "
                "under the port's --run-key (%d of %d)" % (name, same, total, keys, total))
        if name == "B":
            verdict(first_rej == E, "B: every wave rejected its first attempt, before any step was accepted (%d of %d)"
                    % (first_rej, E))
        if name == "C":
            verdict(late_rej >= 1, "C: %d wave(s) rejected an attempt AFTER an accepted step" % late_rej)
        if name == "F":
            verdict(sample_lines(out) == sample_lines(records["A"][4]),
                    "F: the test particles' masses move no bit: every sample line is A's (%d lines), as the port's "
                    "testparticle_type 0 makes it; a step that took their reaction on Saturn or their kinetic energy "
                    "would move them" % len(sample_lines(out)))
        if name == "D1":
            c = per.get(0, {})
            verdict(int(c.get("steps_rejected", 0)) >= 1 and int(c.get("iterations_max_exceeded", 0)) >= 1,
                    "D1: kepler_e099 rejects (%s) and hits the corrector's cap (%s)"
                    % (c.get("steps_rejected"), c.get("iterations_max_exceeded")))
        if n_active is not None and name in ("A", "B", "C", "E"):
            member = int(a[a.index("--member") + 1]) if "--member" in a else None
            phys = ring_physics(out, prob, n_active, member)
            ok, text = physics_verdict(phys, skip=(1,) if name == "C" else ())
            verdict(ok, "%s: the physics%s: %s" % (name, ", the flyby's wave aside" if name == "C" else "", text))
            if name == "C":
                f = phys[1]
                print("    C: the flyby's wave: its particles' distance from Saturn moved %.1f km at most (%s, %.0f km "
                      "from %s at a sample), each swept at least %.4f rad"
                      % (f["dr"], f["circ_who"], f["near"][1], f["near"][0], f["swept"]))

    # 5. controls -------------------------------------------------------------
    print("5. controls, each of which must fail")
    recA, probA, aA, EA, outA = records["A"]
    for label, exe, want in (("the flag checks compiled out, contraction off", unchecked, "PASS"),
                             ("the flag checks compiled out, a*b + c fused (contraction on)", contract, "FAIL")):
        rec = os.path.join(work, "A_%s.rec" % os.path.basename(exe).split(".")[0])
        rc, out, err = run([exe] + dev + ["--problem", probA, "--quiet"] + aA)
        write(rec, out)
        rc, cout = compare(rec, probA, cft, EA)
        members, parts, res = compare_summary(cout)
        fails = [m for m in members if m[2] == "FAIL"]
        first = next((l for l in cout.splitlines() if re.match(r"^member \d+ .*: FAIL:", l)), "")
        verdict(res == want, "contraction control, %s: compare_run %s%s" % (
            label, res, (" (%s)" % first.strip()[:240]) if first else ", %d of %d members IDENTICAL" % (
                sum(1 for m in members if m[2] == "IDENTICAL"), EA)))
    # one hex digit planted in a sampled member's record: a column on two lines
    recE, probE, aE, _, outE = records["E"]
    lines = outE.splitlines()
    sl = [i for i, l in enumerate(lines) if l.startswith("sample ")]
    cols = [("t", 3), ("E", 6), ("vz of the last body", -4), ("t_lo", -1)]
    if args.quick:
        cols = cols[1:2] + cols[3:]
    for cname, tok in cols:
        mut = list(lines)
        planted = 0
        for i in sl[1:3]:
            t = mut[i].split()
            v = t[tok]
            j = len(v) - 1 if "p" not in v else v.index("p") - 1
            if not v[j].isalnum() or v[j] in "x.":
                j = 2
            c = v[j]
            nc = "1" if c != "1" else "2"
            t[tok] = v[:j] + nc + v[j + 1:]
            mut[i] = " ".join(t)
            planted += 1
        rec = os.path.join(work, "E_planted.rec")
        write(rec, "\n".join(mut) + "\n")
        rc, cout = compare(rec, probE, cft, 1)
        line = next((l for l in cout.splitlines() if re.match(r"^member \d+ .*: FAIL:", l)), "")
        verdict(rc == 1 and "FAIL" in line and planted == 2,
                "a digit planted in %s on %d sample lines of the member's record: compare_run rc %d: %s"
                % (cname, planted, rc, line.strip()[:260]))
    # a trailer counter changed
    recB, probB, aB, EB, outB = records["B"]
    mut = re.sub(r"^(# system 0 steps_done=\S+ iterations_max_exceeded=\S+ steps_rejected=)(\d+)",
                 lambda m: m.group(1) + str(int(m.group(2)) + 1), outB, flags=re.M)
    mut = re.sub(r"^(# steps_done=\S+ iterations_max_exceeded=\S+ steps_rejected=)(\d+)",
                 lambda m: m.group(1) + str(int(m.group(2)) + 1), mut, flags=re.M)
    rec = os.path.join(work, "B_counter.rec")
    write(rec, mut)
    rc, cout = compare(rec, probB, cft, EB)
    line = next((l for l in cout.splitlines() if "FAIL" in l and "trailer" in l), "")
    verdict(rc == 1 and mut != outB and "steps_rejected" in line,
            "member 0's steps_rejected one more, the total with it: compare_run rc %d: %s" % (rc, line.strip()[:240]))
    # bit-identity is not validity (CLAUDE.md): wave 0 of A's ring with its
    # zonal line deleted (not regenerated, so the speeds are still its
    # field's) runs, and this program's record is the port's bit for bit;
    # the physics, held to the field the speeds were made for, must fail it
    noz = os.path.join(work, "ring_no_zonal.txt")
    write(noz, "".join(l for l in open(probA) if not l.startswith("zonal ")))
    za = ["--member", "0"] + aA
    rc1, out1, err1 = run([wave] + dev + ["--problem", noz, "--quiet"] + za)
    rc2, out2, err2 = run([cft, "--format", "fp64", "--arith", "fma", "--engine", "loop", "--problem", noz, "--quiet"] + za)
    write(os.path.join(work, "A_no_zonal.rec"), out1)
    same = rc1 == 0 and rc2 == 0 and sample_lines(out1) == sample_lines(out2) and len(sample_lines(out1)) > 1
    caught, text = False, err1.strip()[-200:]
    if rc1 == 0:
        php = ring_physics(out1, probA, 7, 0)
        _, text = physics_verdict(php)
        p0 = php[0]
        caught = p0["gap"] <= GAP_MAX and p0["circ"] > CIRC_BOUND and abs(p0["mm"]) > MM_BOUND   # each, alone
    verdict(same and caught,
            "the physics control: wave 0 of A's ring with its zonal line deleted: this program's %d sample lines are "
            "the port's (%s), and the physics FAILS it, the distance and the angle each: %s"
            % (len(sample_lines(out1)), "IDENTICAL" if same else "NOT identical, rc %d/%d" % (rc1, rc2), text))

    # 6. flags ------------------------------------------------------------------
    print("6. flags: the port aborts, and this program stops by name")
    flagcases = [
        ("DIVBYZERO", "two bodies at one position", "G 0x1p+0\nbody a 0x1p+0 0x1p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0\n"
         "body b 0x1p+0 0x1p+0 0x0p+0 0x0p+0 0x0p+0 0x1p+0 0x0p+0\n", 2),
        ("INVALID", "G = 0 and two bodies at one position", "G 0x0p+0\nbody a 0x1p+0 0x1p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0\n"
         "body b 0x1p+0 0x1p+0 0x0p+0 0x0p+0 0x0p+0 0x1p+0 0x0p+0\n", 1),
        ("OVERFLOW", "two bodies 2^600 apart", "G 0x1p+0\nbody a 0x1p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0\n"
         "body b 0x1p+0 0x1p+600 0x0p+0 0x0p+0 0x0p+0 0x1p+0 0x0p+0\n", 4),
        ("UNDERFLOW", "two bodies 2^-540 apart", "G 0x1p+0\nbody a 0x1p-40 0x0p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0 0x0p+0\n"
         "body b 0x1p-40 0x1p-540 0x0p+0 0x0p+0 0x0p+0 0x1p+0 0x0p+0\n", 8),
    ]
    for fname, what, body, bit in flagcases:
        prob = os.path.join(work, "flag_%s.txt" % fname)
        write(prob, "name flag_%s\n%s" % (fname.lower(), body))
        a = ["--problem", prob, "--dt", "0.01", "--steps", "4", "--sample", "2", "--quiet"]
        r1, o1, e1 = run([cft, "--format", "fp64", "--arith", "fma"] + a)
        r2, o2, e2 = run([wave] + dev + a)
        m1 = re.search(r"raised 0x([0-9a-f]+)", e1)
        verdict(r1 == 3 and r2 == 3 and m1 and int(m1.group(1), 16) & bit and fname in e2,
                "%s (%s): the port stops, rc %d, \"%s\"; this program stops, rc %d, \"%s\"" % (
                    fname, what, r1, (m1.group(0) if m1 else e1.strip()[-80:]), r2,
                    re.sub(r"\. This workload.*", "", e2.strip())[:200]))
        r1, o1, e1 = run([cft, "--format", "fp64", "--arith", "fma", "--no-flag-abort"] + a)
        r2, o2, e2 = run([wave] + dev + ["--no-flag-abort"] + a)
        f1 = re.search(r"flags_seen=0x([0-9a-f]+)", o1)
        f2 = re.search(r"flags_raised=0x([0-9a-f]+)", o2)
        verdict(r1 == 0 and r2 == 0 and f1 and f2 and int(f1.group(1), 16) & 0xf == int(f2.group(1), 16),
                "%s with --no-flag-abort: both run on; the port's flags_seen=0x%s (less inexact) is this program's "
                "flags_raised=0x%s" % (fname, f1.group(1) if f1 else "?", f2.group(1) if f2 else "?"))

    if FAILS:
        print("RESULT: FAIL: %d check(s) failed; the work is kept in %s" % (len(FAILS), work))
        return 1
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    print("RESULT: PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
