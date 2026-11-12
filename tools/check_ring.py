#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of the ring (ROUND4.md, P2): tools/make_ring.py's Saturn ring,
# run at every format and held to bits AND to physics.
#
#   python tools/check_ring.py [--build build] [--quick]
#   (`make check` runs it whole; check-quick and check-light pass --quick:
#   binary64 only, and fewer steps.)
#
# The ring is generated here, into the build tree (build/check_ring/),
# never committed: 228 particles, seed 1, --no-zonal, so 4 waves of 64
# bodies - Saturn and six moons, then 57 test particles each. The zonal
# line is refused until round 4's P1 lands; the lead's P1 x P2 seam test
# runs the ring with it.
#
# What it asserts, and what lets each assertion fail:
#   1. The generator. Twice the same bytes, and different bytes from
#      another seed; the same bytes again on mpmath's pure-Python backend
#      (MPMATH_NOGMPY), the closest this host comes to a second host.
#      Then the file, read back as the program reads it: N, E and the
#      active line, N the generator's WAVE (64) and active its MASSIVE
#      (7); Saturn and the moons first in every system, bit for bit the
#      same in every wave, and the only massive bodies; the particles
#      sorted by radius across the waves and inside the annulus; the
#      barycentre of the massive bodies at the origin; every particle
#      circular in the file's own field (|v|^2 r / GM is the field's
#      factor, 1 with no zonal line, and r.v is 0); every moon inclined:
#      sin i, from r x v relative to Saturn, above zero (not |z|: Pan
#      starts at its node, z = 0; a flat moon would leave its z terms and
#      the z of its reaction on Saturn exactly zero, V2); no start within
#      the stated Hill radii of a moon, the Hill radii derived from the
#      file's own masses and states. The zonal line, written without
#      --no-zonal, holds the generator's constants, and its particles and
#      its circular moon (Pan) the oblate field's circular speed. The
#      exclusion bites (a wide one refuses draws) and is bounded (one that
#      leaves none of the annulus is refused by name), a count that does
#      not fill waves of 57 is refused by name, and --one-system writes
#      the same body lines as one wave of the same draws.
#   2. The engines: one wave at every format. The loop and the program
#      engine in the FMA form agree bit for bit, every record line and
#      every trailer counter, and the program engine's predictor ran. The
#      loop engine in REBOUND's form runs it too: it must differ from the
#      FMA form (a different rounding sequence; an --arith that did
#      nothing would pass everything else) and agree with it to within
#      2^(-p/2) relative, p the format's precision. MEASURED at binary64
#      on w00001: the two forms' dt_next differ by 1.7e-11 relative after
#      three steps, far above the positions' rounding, and 2^-26.5 clears
#      it 650 times. BELIEVED the reason: PRS23 builds the step from the
#      high b coefficients, where the corrector's divided differences
#      amplify round-off.
#   3. Test particles are inert on the real problem: the wave with its
#      `active 7` line deleted - every one of its 64 bodies active, the 57
#      particles massless - gives the same record at every format, on the
#      loop engine in REBOUND's form, and at binary64 in the FMA form on
#      both engines too. And REBOUND's own IAS15 (build/ias15_ref, its
#      build time printed, its record's header held to impl=rebound, so a
#      stand-in at that path fails) runs the wave as one system with every
#      body active (it reads no active line): bit for bit ias15_cft's at
#      binary64 in REBOUND's form and at the program's defaults (REBOUND's
#      12 passes), which are the same samples as the gate's all-active
#      run, so REBOUND has the bits `active 7` has.
#   4. The ensemble: each of the 4 waves, run alone (--member), gets the
#      bits it gets in the ensemble, every sample and its trailer, on the
#      program engine (the FMA form; its ensemble path predicts each
#      system with its own step ratio, and here the waves' steps differ).
#      And the comparison can fail: in a copy with one particle of one
#      wave moved by one ulp in x, that wave differs from its own run
#      alone and every other wave still matches. Then the ring as one
#      system of 37 bodies runs, at binary64.
#   5. The physics, printed beside the bits and asserted:
#      - 420 pairs a force evaluation, read from every active run's
#        `# active` line, against REBOUND's two loops' count from that
#        run's own problem file (n(n-1)/2 + (N-n)n), and no such line in
#        a run with every body active;
#      - the moons act: the wave with the six moons' masses set to zero
#        (the same starting bits) differs from the wave with them, by
#        more than round-off could;
#      - no particle leaves the annulus or reaches Saturn, at any sample
#        of any run;
#      - every wave's largest |z| from Saturn is above zero at every
#        sample (a flat ring leaves every z term of J2 and J4 exactly
#        zero, ROUND4.md "the trap"); this counter, and the moons'
#        inclination counter of part 1, read zero on a flat copy of the
#        file, so each can fail;
#      - the far wave (the innermost, wholly inside every moon's
#        pericentre, which the file is held to) takes within 10% of 36.2
#        steps an orbit at epsilon 1e-9: cft-fp256's GPU-E section 3,
#        MEASURED there on a circular orbit about a point mass, REBOUND's
#        form, --max-iter 60. A miss is a finding, not a bound to loosen.
#        The orbit is its shortest, from the file's own state; the steps
#        and the time are the record's.
#      Beside them, MEASURED and printed, not asserted: each wave's step
#      in the ensemble against its shortest orbit and its innermost
#      particle's. The shortest orbit in a system sets its step (PRS23's
#      minimum runs over every particle), so a wave outside Pan's orbit
#      steps at Pan's rate.
#
# Every case is a step-count case (--steps); none integrates to a time,
# none is a step sequence (--dt-file) and none merges, so none has a
# CLAUDE.md --light rule to keep; --quick scales them itself.

import argparse
import math
import os
import re
import subprocess
import sys
import time
from fractions import Fraction

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import compare as equivalent, hex_to_fraction  # noqa: E402
import make_ring  # noqa: E402 - the wave's size, the massive count and the zonal constants, from their one source

root = os.path.dirname(here)

RING = ["--particles", "228", "--seed", "1", "--no-zonal"]
RING_SEED2 = ["--particles", "228", "--seed", "2", "--no-zonal"]
ENGINE_WAVE = 1          # straddles Pan's orbit: the moons act most here
FAR_WAVE = 0             # the innermost wave
SPO_REF, SPO_TOL = 36.2, 0.10   # GPU-E section 3, MEASURED there; the brief's bound
EPS, DT0, MAX_ITER = "1e-9", "100", "60"
# IEEE 754's precisions, by definition; checked against each record's own
# pc_tol_shift, which the program sets to the format's p - 53
PREC = {"fp64": 53, "fp128": 113, "fp256": 237}
COUNTERS = ("steps_done", "steps_rejected", "iterations_max_exceeded",
            "mean_pc_iterations", "max_pc_iterations", "flags_seen")
SYS_COUNTERS = ("steps_done", "steps_rejected", "iterations_max_exceeded",
                "mean_pc_iterations", "max_pc_iterations")
ACTIVE_LINE = re.compile(r"^# active (\d+) of (\d+) bodies a system: (\d+) pairs a force evaluation, a system$", re.M)


# ---------------------------------------------------------------- files

def read_problem(path):
    """The problem as ias15_cft reads it: G, N, E, active, zonal, the
    systems (name and body lines, as tokens), and the `# ring` lines."""
    p = {"G": None, "N": None, "E": None, "active": None, "zonal": None, "systems": [], "ring": {}}
    cur = None
    for line in open(path):
        t = line.split()
        if not t:
            continue
        if t[0] == "#":
            if len(t) > 2 and t[1] == "ring":
                p["ring"][t[2]] = t[3:]
            continue
        if t[0].startswith("#"):
            continue
        if t[0] in ("G",):
            p["G"] = float.fromhex(t[1])
        elif t[0] in ("N", "E", "active"):
            p[t[0]] = int(t[1])
        elif t[0] == "zonal":
            p["zonal"] = (int(t[1]),) + tuple(float.fromhex(x) for x in t[2:5])
        elif t[0] == "system":
            cur = {"name": t[1], "bodies": []}
            p["systems"].append(cur)
        elif t[0] == "body":
            if cur is None:
                cur = {"name": None, "bodies": []}
                p["systems"].append(cur)
            cur["bodies"].append(t[1:9])
    return p


def state(tokens):
    """A body line's (m, x, y, z, vx, vy, vz) as binary64 values."""
    return [float.fromhex(x) for x in tokens[1:8]]


def rel(b, s):
    """Body b's state relative to body s (both state() lists), position then velocity."""
    return [b[1 + c] - s[1 + c] for c in range(6)]


def norm(v):
    return math.sqrt(sum(x * x for x in v))


def elements(b, s, G):
    """b's osculating semi-major axis, eccentricity and period about s."""
    d = rel(b, s)
    rv, vv = d[:3], d[3:]
    r, v2, mu = norm(rv), sum(x * x for x in vv), G * (s[0] + b[0])
    a = 1.0 / (2.0 / r - v2 / mu)
    rdotv = sum(p * q for p, q in zip(rv, vv))
    ev = [((v2 - mu / r) * rv[c] - rdotv * vv[c]) / mu for c in range(3)]
    return a, norm(ev), 2.0 * math.pi * math.sqrt(a ** 3 / mu)


def largest_z(sat, particles):
    """The flat-ring counter: the largest |z| of any particle from Saturn
    (positions as (x, y, z, ...) lists)."""
    return max(abs(p[2] - sat[2]) for p in particles)


def inclination(b, s):
    """b's orbit about s (state() lists), from h = r x v: (sin i, i in
    radians). sin i is |(h_x, h_y)| / |h|: exactly zero for an orbit in the
    plane, prograde or retrograde, and unlike |z| not zero for a body that
    starts on its node."""
    d = rel(b, s)
    h = [d[1] * d[5] - d[2] * d[4], d[2] * d[3] - d[0] * d[5], d[0] * d[4] - d[1] * d[3]]
    tilt = math.hypot(h[0], h[1])
    return tilt / norm(h), math.atan2(tilt, h[2])


def edit_copy(src, dst, fn):
    """Write a copy of a problem file with each line passed through fn
    (which returns the line, a replacement, or None to drop it)."""
    with open(src) as f, open(dst, "w", newline="\n") as g:
        for line in f:
            out = fn(line)
            if out is not None:
                g.write(out)
    return dst


# ---------------------------------------------------------------- records

def parse(out):
    """A record: {system: [sample tokens]} (each list from the sample
    number on, `|` included), the whole trailer, and the per-system
    trailers."""
    recs, whole, systems = {}, {}, {}
    for line in out.splitlines():
        if line.startswith("sample "):
            recs.setdefault(0, []).append(line.split()[1:])
        elif line.startswith("system ") and " sample " in line:
            t = line.split()
            recs.setdefault(int(t[1]), []).append(t[3:])
        elif line.startswith("# steps_done="):
            whole = dict(kv.partition("=")[::2] for kv in line[2:].split())
        elif line.startswith("# system ") and " steps_done=" in line:
            t = line.split()
            systems[int(t[2])] = dict(kv.partition("=")[::2] for kv in t[3:])
    return recs, whole, systems


def bodies_at(sample, n):
    """A sample's body states as floats: [(x, y, z, vx, vy, vz)] * n."""
    vals = [float.fromhex(x) for x in sample[6:6 + 6 * n]]
    return [vals[6 * i:6 * i + 6] for i in range(n)]


def first_difference(a, b):
    """Where two sample-token lists first differ: (sample index, token
    index, a's token, b's token), or None."""
    if len(a) != len(b):
        return (min(len(a), len(b)), -1, "%d samples" % len(a), "%d samples" % len(b))
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            j = next((j for j, (p, q) in enumerate(zip(x, y)) if p != q), min(len(x), len(y)))
            return (i, j, x[j] if j < len(x) else "(none)", y[j] if j < len(y) else "(none)")
    return None


def column(j, n):
    """A sample token's name, for a FAIL line."""
    lead = ["sample", "step", "t", "dt_next", "dt_last", "E"]
    if j < 0:
        return "the sample count"
    if j < len(lead):
        return lead[j]
    k = j - len(lead)
    if k < 6 * n:
        return "body %d %s" % (k // 6, ("x", "y", "z", "vx", "vy", "vz")[k % 6])
    return "the exact-time pair"


def largest_relative(a, b, n):
    """The largest relative difference between two runs' samples: per body,
    |dposition| / |position| and |dvelocity| / |velocity|, and t, the
    steps and E against their own size, all exact (fractions of the hex)."""
    worst = Fraction(0)
    for x, y in zip(a, b):
        for j in range(2, 6):
            p, q = hex_to_fraction(x[j])[0], hex_to_fraction(y[j])[0]
            if p != q and p != 0:
                worst = max(worst, abs((q - p) / p))
        for i in range(n):
            for half in (0, 3):
                u = [hex_to_fraction(x[6 + 6 * i + half + c])[0] for c in range(3)]
                w = [hex_to_fraction(y[6 + 6 * i + half + c])[0] for c in range(3)]
                d2 = sum((p - q) ** 2 for p, q in zip(u, w))
                n2 = sum(p * p for p in u)
                if d2 and n2:
                    worst = max(worst, Fraction(math.sqrt(d2 / n2)))
    return float(worst)


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true", help="binary64 only, and fewer steps")
    args = ap.parse_args()
    t_start = time.time()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, "ias15_cft" + exe)
    progs = os.path.join(root, "programs", "out")
    gen = os.path.join(here, "make_ring.py")
    py = sys.executable
    scratch = os.path.join(args.build, "check_ring")
    os.makedirs(scratch, exist_ok=True)
    formats = ["fp64"] if args.quick else ["fp64", "fp128", "fp256"]
    # a step-count case at every format; the wide formats' steps cost more
    steps = {"fp64": 4, "fp128": 3, "fp256": 2} if args.quick else {"fp64": 6, "fp128": 3, "fp256": 2}
    spo_steps = 16 if args.quick else 24
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= bool(ok)
        print("  %s %s" % ("ok  " if ok else "FAIL", text))
        return ok

    ref_exe = os.path.join(args.build, "ias15_ref" + exe)
    for binary in (cft, ref_exe):
        if not os.path.exists(binary):
            print("check_ring: no %s" % binary)
            print("RESULT: FAIL")
            return 1
    images = [os.path.join(progs, f) for f in os.listdir(progs) if f.endswith(".cftp")] if os.path.isdir(progs) else []
    stamp = lambda t: time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(t))
    print("check_ring: %s, built %s; %s, built %s; %d program images in %s, the newest %s" % (
        cft, stamp(os.path.getmtime(cft)), ref_exe, stamp(os.path.getmtime(ref_exe)), len(images), progs,
        stamp(max(os.path.getmtime(f) for f in images)) if images else "none"))

    def generate(out, extra, env=None):
        e = {k: v for k, v in os.environ.items() if k != "MPMATH_NOGMPY"}   # mpmath's own choice, unless asked
        e.update(env or {})
        p = subprocess.run([py, "-B", gen, out] + extra, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           text=True, cwd=root, env=e)
        return p.returncode, p.stderr.strip()

    def readbytes(path):
        with open(path, "rb") as f:
            return f.read()

    # ------------------------------------------------------------ 1. the generator
    print("check_ring: the generator, tools/make_ring.py %s" % " ".join(RING))
    ring = os.path.join(scratch, "ring.txt")
    again = os.path.join(scratch, "ring_again.txt")
    seed2 = os.path.join(scratch, "ring_seed2.txt")
    pure = os.path.join(scratch, "ring_purepython.txt")
    rc1, say1 = generate(ring, RING)
    rc2, _ = generate(again, RING)
    rc3, _ = generate(seed2, RING_SEED2)
    rc4, say4 = generate(pure, RING, {"MPMATH_NOGMPY": "1"})
    if rc1 or rc2 or rc3 or rc4:
        print("  FAIL the generator failed (rc %d %d %d %d): %s" % (rc1, rc2, rc3, rc4, say1[-300:]))
        print("RESULT: FAIL")
        return 1
    print("       %s" % say1.replace("make_ring: wrote ", "")[-260:])
    b1, b2, b3, b4 = readbytes(ring), readbytes(again), readbytes(seed2), readbytes(pure)
    backends = (re.search(r"(\w+) backend\)", say1), re.search(r"(\w+) backend\)", say4))
    backends = tuple(m.group(1) if m else "?" for m in backends)
    verdict(b1 == b2 and b1 != b3, "the same arguments twice give the same %d bytes, and seed 2 gives other bytes (%s)"
            % (len(b1), "differing" if b1 != b3 else "IDENTICAL"))
    verdict(b1 == b4, "mpmath's %s backend and its %s backend give the same bytes%s"
            % (backends[0], backends[1], "" if backends[0] != backends[1] else " (one backend only on this host)"))

    prob = read_problem(ring)
    G, N, E, na = prob["G"], prob["N"], prob["E"], prob["active"]
    systems = prob["systems"]
    annulus = [float(x) for x in prob["ring"].get("annulus_km", [])]
    r_saturn = float(prob["ring"].get("saturn_radius_km", ["nan"])[0])
    excl = dict(kv.partition("=")[::2] for kv in prob["ring"].get("exclusion", []))
    k_hill = float(excl.get("hill", "nan"))
    ok = (G and N and E and na and len(systems) == E and len(annulus) == 2 and r_saturn > 0 and k_hill > 0
          and all(len(s["bodies"]) == N for s in systems) and prob["zonal"] is None
          and N == make_ring.WAVE and na == make_ring.MASSIVE)
    if not verdict(ok, "the file: N=%s E=%s active=%s (make_ring's WAVE %d and MASSIVE %d), %d systems of %s bodies, "
                       "no zonal line under --no-zonal, its annulus %s km, Saturn's radius %s km and the exclusion's "
                       "%s Hill radii read from its head"
                       % (N, E, na, make_ring.WAVE, make_ring.MASSIVE, len(systems),
                          sorted(set(len(s["bodies"]) for s in systems)), annulus, r_saturn, k_hill)):
        print("RESULT: FAIL")
        return 1
    massive = [systems[0]["bodies"][i] for i in range(na)]
    same = all(s["bodies"][:na] == massive for s in systems)
    masses = [state(b)[0] for b in massive]
    pmass = set(state(b)[0] for s in systems for b in s["bodies"][na:])
    verdict(same and massive[0][0] == "saturn" and all(m > 0 for m in masses) and pmass == {0.0},
            "Saturn first, then the moons (%s): the %d massive bodies, bit for bit the same in all %d waves; "
            "every particle massless" % (", ".join(b[0] for b in massive[1:]), na, E))
    sat = state(massive[0])
    moons = [(b[0], state(b)) for b in massive[1:]]
    # the barycentre of the massive bodies, exactly
    mt = sum(Fraction(m) for m in masses)
    cen = [sum(Fraction(state(b)[0]) * Fraction(state(b)[1 + c]) for b in massive) / mt for c in range(6)]
    cpos, cvel = float(max(abs(c) for c in cen[:3])), float(max(abs(c) for c in cen[3:]))
    verdict(cpos < 1e-9 and cvel < 1e-12, "the frame: the massive bodies' barycentre is at %.2e km and moves at "
            "%.2e km/s" % (cpos, cvel))
    # sorted by radius, inside the annulus, circular, inclined
    radii, worst_circ, worst_radial, zmin_wave, order_ok = [], 0.0, 0.0, None, True
    for s in systems:
        rs = []
        for b in s["bodies"][na:]:
            d = rel(state(b), sat)
            r = norm(d[:3])
            rs.append(r)
            worst_circ = max(worst_circ, abs(sum(x * x for x in d[3:]) * r / (G * sat[0]) - 1.0))
            worst_radial = max(worst_radial, abs(sum(d[c] * d[3 + c] for c in range(3))) / (r * norm(d[3:])))
        if rs != sorted(rs) or (radii and rs[0] < radii[-1]):
            order_ok = False
        radii += rs
        z = largest_z(state(s["bodies"][0])[1:], [state(b)[1:] for b in s["bodies"][na:]])
        zmin_wave = z if zmin_wave is None else min(zmin_wave, z)
    verdict(order_ok and annulus[0] <= min(radii) and max(radii) < annulus[1],
            "the particles: sorted by radius across the %d waves, from %.3f to %.3f km, inside the annulus "
            "[%g, %g)" % (E, min(radii), max(radii), annulus[0], annulus[1]))
    verdict(worst_circ < 2.0 ** -45 and worst_radial < 2.0 ** -40,
            "circular in the file's own field (a point mass's): |v|^2 r / GM - 1 at most %.1e, r.v / (|r||v|) at "
            "most %.1e" % (worst_circ, worst_radial))
    verdict(zmin_wave > 0, "not flat: every wave's largest |z| from Saturn is at least %.3f km" % zmin_wave)
    # The moons out of the plane, by inclination. A flat moon leaves its J2
    # and J4 z terms, and the z of its reaction on Saturn, exactly zero, and
    # the moons carry all of the ring's reactions, the particles being
    # massless (V2, 2026-11-12). Not by |z|: Pan starts at its node, z = 0.
    incl = [(name, inclination(m, sat), m[3] - sat[3]) for name, m in moons]
    verdict(all(sin_i > 0 for _, (sin_i, _), _ in incl),
            "the moons are not flat: each one's inclination from r x v relative to Saturn is above zero (sin i > 0): "
            "%s degrees%s" % (", ".join("%s %.4f" % (n, math.degrees(i)) for n, (_, i), _ in incl),
                              "".join("; %s's z from Saturn is %r km, at its node, which |z| would call flat" % (n, z)
                                      for n, (sin_i, _), z in incl if z == 0 and sin_i > 0)))
    # the Hill radii from the file's own state, and no start within the
    # exclusion's reach of a moon
    hill = {}
    for name, m in moons:
        a = elements(m, sat, G)[0]
        hill[name] = a * (m[0] / (3.0 * sat[0])) ** (1.0 / 3.0)
    said = dict(kv.partition("=")[::2] for kv in prob["ring"].get("hill_km", []))
    agree = all(name in said and abs(float(said[name]) / hill[name] - 1) < 1e-9 for name in hill)
    closest = min((norm(rel(state(b), m)[:3]) / hill[name], name)
                  for s in systems for b in s["bodies"][na:] for name, m in moons)
    verdict(agree and closest[0] >= k_hill,
            "Hill radii from the file's states (%s km)%s; the closest start to a moon is %.1f Hill radii from %s, "
            "past the %g refused" % (", ".join("%s %.1f" % kv for kv in hill.items()),
                                     ", as the head says" if agree else ", NOT as the head says",
                                     closest[0], closest[1], k_hill))
    # the zonal line, the exclusion biting, and a count that is refused
    zring = os.path.join(scratch, "ring_zonal.txt")
    rc, say = generate(zring, ["--particles", "57", "--seed", "3"])
    zp = read_problem(zring) if rc == 0 else None
    want = (0, float(make_ring.ZONAL_R_KM), float(make_ring.ZONAL_J2), float(make_ring.ZONAL_J4))
    zworst, circular = None, [m[0] for m in make_ring.MOONS if float(m[3]) == 0]
    if zp and zp["zonal"] == want:
        k, R, J2, J4 = zp["zonal"]
        zs = state(zp["systems"][0]["bodies"][k])
        zworst = 0.0
        # every particle, and every moon on a circular orbit (Pan), which
        # the oblate factor must keep circular
        for b in zp["systems"][0]["bodies"]:
            if b[0] not in circular and state(b)[0] != 0:
                continue
            d = rel(state(b), zs)
            r = norm(d[:3])
            q = (R / r) ** 2
            factor = 1.0 + 1.5 * J2 * q - 1.875 * J4 * q * q
            zworst = max(zworst, abs(sum(x * x for x in d[3:]) * r / (G * (zs[0] + state(b)[0])) / factor - 1.0))
    verdict(zworst is not None and zworst < 2.0 ** -45,
            "without --no-zonal: `zonal %s` holds R, J2 and J4, and every particle, and %s, has the oblate field's "
            "circular speed, 1 + 3/2 J2 (R/r)^2 - 15/8 J4 (R/r)^4 over the point mass's, to %s"
            % (" ".join("%r" % x for x in want), " and ".join(circular) or "no moon",
               "%.1e" % zworst if zworst is not None else "nothing (rc %d, zonal line %s)" % (rc, zp and zp["zonal"])))
    wide = os.path.join(scratch, "ring_exclusion.txt")
    rc, say = generate(wide, ["--particles", "57", "--seed", "1", "--no-zonal", "--exclude-hill", "300"])
    wp = read_problem(wide) if rc == 0 else None
    refused, wclosest = -1, -1.0
    if wp:
        refused = int(dict(kv.partition("=")[::2] for kv in wp["ring"]["exclusion"])["refused"])
        wclosest = min(norm(rel(state(b), state(m))[:3]) / hill[m[0]]
                       for b in wp["systems"][0]["bodies"][wp["active"]:] for m in wp["systems"][0]["bodies"][1:wp["active"]])
    verdict(refused > 0 and wclosest >= 300,
            "the exclusion bites: at 300 Hill radii it refused %d draws, and the closest start kept is %.1f Hill "
            "radii from a moon" % (refused, wclosest))
    # and it is bounded: one that leaves none of the annulus is refused by
    # name, where it once drew forever (V2)
    unbounded = os.path.join(scratch, "ring_unbounded.txt")
    if os.path.exists(unbounded):
        os.remove(unbounded)
    rc, say = generate(unbounded, ["--particles", "57", "--seed", "1", "--no-zonal", "--exclude-hill", "1000"])
    verdict(rc != 0 and "Lower --exclude-hill" in say and not os.path.exists(unbounded),
            "an exclusion that leaves none of the annulus, 1000 Hill radii, is refused by name and writes nothing "
            "(rc %d): %s" % (rc, say[-170:]))
    rc, say = generate(os.path.join(scratch, "ring_56.txt"), ["--particles", "56"])
    verdict(rc != 0 and "multiple of 57" in say, "56 particles, which do not fill a wave, are refused by name (rc %d): %s"
            % (rc, say[-110:]))
    # one system: the same draws as one wave are the same body lines
    wave57, one57, one30 = (os.path.join(scratch, f) for f in ("ring_wave57.txt", "ring_one57.txt", "ring_one30.txt"))
    rcs = [generate(wave57, ["--particles", "57", "--seed", "5", "--no-zonal"])[0],
           generate(one57, ["--particles", "57", "--seed", "5", "--no-zonal", "--one-system"])[0],
           generate(one30, ["--particles", "30", "--seed", "5", "--no-zonal", "--one-system"])[0]]
    pa, pb, pc = (read_problem(f) if r == 0 else None for f, r in zip((wave57, one57, one30), rcs))
    one_ok = bool(pa and pb and pc) and pa["E"] == 1 and pb["E"] is None and pc["E"] is None and \
        pa["systems"][0]["bodies"] == pb["systems"][0]["bodies"] and pb["N"] == len(pb["systems"][0]["bodies"]) and \
        len(pc["systems"]) == 1 and pc["N"] == len(pc["systems"][0]["bodies"]) and pc["active"] == na
    if one_ok:
        cs = state(pc["systems"][0]["bodies"][0])
        rr = [norm(rel(state(b), cs)[:3]) for b in pc["systems"][0]["bodies"][na:]]
        one_ok = rr == sorted(rr)
    verdict(one_ok, "--one-system: 57 particles as one system are the same %s body lines as the same draws as one "
            "wave, with no E or system line; 30 particles are one system of N = %s, `active %s`, sorted by radius"
            % (pb and pb["N"], pc and pc["N"], pc and pc["active"]))

    # ------------------------------------------------------------ the runs
    cache = {}

    def run(key, problem, fmt, engine, arith, nsteps, sample, member=None, dt=DT0, max_iter=MAX_ITER):
        """ias15_cft on a problem; each distinct run once. max_iter None
        leaves the program's default, REBOUND's 12 passes."""
        cmd = [cft, "--format", fmt, "--problem", problem, "--dt", dt, "--epsilon", EPS,
               "--steps", str(nsteps), "--sample", str(sample)] + \
              (["--max-iter", max_iter] if max_iter is not None else []) + \
              ["--quiet", "--engine", engine, "--arith", arith]
        if engine == "program":
            cmd += ["--programs", progs]
        if member is not None:
            cmd += ["--member", str(member)]
        ck = tuple(cmd)
        if ck not in cache:
            t0 = time.time()
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=root)
            cache[ck] = (p.returncode, p.stdout, p.stderr, time.time() - t0, problem, member, key)
        rc, out, err = cache[ck][:3]
        if rc:
            print("  FAIL %s: rc %d: %s" % (key, rc, err.strip()[-300:]))
        return rc, out

    def pairs_of(p):
        """REBOUND's two loops' pair count, from a problem's own N and
        active: n(n-1)/2 active pairs, then (N-n)n test-active ones."""
        return p["active"] * (p["active"] - 1) // 2 + (p["N"] - p["active"]) * p["active"]

    want_pairs = pairs_of(prob)
    phys = {"rmin": float("inf"), "rmax": 0.0, "excursion": 0.0, "zmin": float("inf"), "samples": 0}
    problems = {}

    def physics(label, out, problem, member):
        """Every sample of a record against the annulus, Saturn and the
        plane, each particle's change of radius against its start in the
        run's own problem file. A miss fails the gate here, by name."""
        nonlocal allok
        if problem not in problems:
            problems[problem] = read_problem(problem)
        p = problems[problem]
        mine = p["systems"] if member is None else [p["systems"][member]]
        # with every body active, the particles are still the ring's, after its na massive bodies
        n, act = p["N"], p["active"] if p["active"] is not None else na
        recs, _, _ = parse(out)
        bad = []
        for s, samples in recs.items():
            start = [norm(rel(state(b), state(mine[s]["bodies"][0]))[:3]) for b in mine[s]["bodies"]]
            for smp in samples:
                bs = bodies_at(smp, n)
                for i in range(act, n):
                    r = norm([bs[i][c] - bs[0][c] for c in range(3)])
                    phys["rmin"], phys["rmax"] = min(phys["rmin"], r), max(phys["rmax"], r)
                    phys["excursion"] = max(phys["excursion"], abs(r - start[i]))
                    if not (annulus[0] <= r < annulus[1]) or r <= r_saturn:
                        bad.append("%s: %s at %.3f km from Saturn at sample %s" % (label, mine[s]["bodies"][i][0], r, smp[0]))
                zmax = largest_z(bs[0], bs[act:])
                phys["zmin"] = min(phys["zmin"], zmax)
                if zmax <= 0:
                    bad.append("%s: system %d is flat at sample %s" % (label, s, smp[0]))
                phys["samples"] += 1
        for b in bad[:3]:
            print("  FAIL the physics, %s" % b)
        allok &= not bad
        return not bad

    names = [[b[0] for b in s["bodies"]] for s in systems]
    r0 = {b[0]: norm(rel(state(b), sat)[:3]) for s in systems for b in s["bodies"][na:]}

    # ------------------------------------------------------------ 2. the engines
    # 3. and test particles inert, beside them
    allact = edit_copy(ring, os.path.join(scratch, "ring_all_active.txt"),
                       lambda l: None if l.startswith("active ") else l)
    print("check_ring: one wave (%s, radii %.0f to %.0f km) on both engines, and every body active against `active %d`"
          % (systems[ENGINE_WAVE]["name"], r0[names[ENGINE_WAVE][na]], r0[names[ENGINE_WAVE][-1]], na))
    for fmt in formats:
        n, smp = steps[fmt], 1
        lab = lambda e, a, w="": "%s %s engine %s form%s" % (fmt, e, a, w)
        rc1, lf = run(lab("loop", "fma"), ring, fmt, "loop", "fma", n, smp, ENGINE_WAVE)
        rc2, pf = run(lab("program", "fma"), ring, fmt, "program", "fma", n, smp, ENGINE_WAVE)
        rc3, lr = run(lab("loop", "rebound"), ring, fmt, "loop", "rebound", n, smp, ENGINE_WAVE)
        if rc1 or rc2 or rc3:
            allok = False
            continue
        (lrec, lw, _), (prec_, pw, _), (rrec, rw, _) = parse(lf), parse(pf), parse(lr)
        d = first_difference(lrec.get(0, []), prec_.get(0, []))
        bad = [] if d is None else ["first difference at sample %d, %s: loop %s, program %s" % (d[0], column(d[1], N), d[2], d[3])]
        bad += ["trailer %s is %s on the loop engine, %s on the program engine" % (k, lw.get(k), pw.get(k))
                for k in COUNTERS if lw.get(k) is None or lw.get(k) != pw.get(k)]
        done, nl, npg = int(lw.get("steps_done", -1)), int(lw.get("nextstep_calls", -1)), int(pw.get("nextstep_calls", -1))
        rej = int(lw.get("steps_rejected", 0))
        if nl != 0 or not done <= npg <= done + rej:
            bad.append("nextstep_calls %d on the loop engine (must be 0) and %d on the program engine (must be %d to %d)"
                       % (nl, npg, done, done + rej))
        if not lrec.get(0):
            bad.append("no samples")
        verdict(not bad, "%s: %s" % (fmt, "; ".join(bad) if bad else
                "the loop and program engines in the FMA form: %d sample lines identical over %d steps (%s passes a "
                "step), counters and flags the loop's, the program's predictor ran %d times"
                % (len(lrec[0]), done, lw.get("mean_pc_iterations"), npg)))
        # REBOUND's form: runs, differs, and only at round-off
        dr = first_difference(lrec[0], rrec.get(0, []))
        rel_ = largest_relative(lrec[0], rrec.get(0, []), N) if rrec.get(0) and len(rrec[0]) == len(lrec[0]) else float("inf")
        bound = 2.0 ** (-PREC[fmt] / 2)
        shift = re.search(r" pc_tol_shift=(\d+) ", lr)
        ok = dr is not None and rel_ < bound and rw.get("steps_done") == lw.get("steps_done") and \
            shift is not None and int(shift.group(1)) == PREC[fmt] - 53
        verdict(ok, "%s: the loop engine in REBOUND's form runs it, %s the FMA form (%s), its largest relative "
                "difference %.3e against a bound of 2^-%g = %.1e (p = %d, its record's pc_tol_shift %s + 53), %s steps "
                "as the FMA form's"
                % (fmt, "differing from" if dr is not None else "IDENTICAL to",
                   "first at sample %d, %s" % (dr[0], column(dr[1], N)) if dr is not None else "an --arith that did nothing",
                   rel_, PREC[fmt] / 2, bound, PREC[fmt], shift.group(1) if shift else "missing", rw.get("steps_done")))
        # 3. inert: every body active, on the loop engine in REBOUND's form,
        # and at binary64 on both FMA rows as well
        rows = [("loop", "rebound", lr)] + ([("loop", "fma", lf), ("program", "fma", pf)] if fmt == "fp64" else [])
        for engine, arith, act_out in rows:
            rc, all_out = run(lab(engine, arith, ", all active"), allact, fmt, engine, arith, n, smp, ENGINE_WAVE)
            if rc:
                allok = False
                continue
            arec, aw, _ = parse(all_out)
            crec, cw, _ = parse(act_out)
            d = first_difference(crec.get(0, []), arec.get(0, []))
            bad = [] if d is None else ["first difference at sample %d, %s: active %s, all active %s" % (d[0], column(d[1], N), d[2], d[3])]
            bad += ["trailer %s is %s with `active`, %s all active" % (k, cw.get(k), aw.get(k))
                    for k in COUNTERS if cw.get(k) != aw.get(k)]
            ma = ACTIVE_LINE.search(act_out)
            if not ma or int(ma.group(3)) != want_pairs:
                bad.append("the active run's `# active` line says %s pairs, and the file's N = %d with active %d give %d"
                           % (ma.group(3) if ma else "nothing: no line", N, na, want_pairs))
            if ACTIVE_LINE.search(all_out):
                bad.append("an `# active` line in the run with every body active")
            if not crec.get(0):
                bad.append("no samples")
            verdict(not bad, "%s %s engine %s form: `active %d` against all %d active, %s" % (
                fmt, engine, arith, na, N, "; ".join(bad) if bad else
                "%d sample lines identical, the trailer too: %s pairs a force evaluation against %d"
                % (len(crec[0]), ma.group(3), N * (N - 1) // 2)))

    # REBOUND's own IAS15 on the wave, at binary64. ias15_ref reads no
    # active, E or system line, so the wave goes to it as one system with
    # every body active: by part 3, the same record as `active 7`.
    single = os.path.join(scratch, "ring_%s_one_system.txt" % systems[ENGINE_WAVE]["name"])
    with open(single, "w", newline="\n") as f:
        f.write("# %s of %s as one system, every body active, for ias15_ref\n" % (systems[ENGINE_WAVE]["name"], ring))
        f.write("name ring_%s\nG %s\nN %d\n" % (systems[ENGINE_WAVE]["name"], G.hex(), N))
        f.writelines("body %s\n" % " ".join(b) for b in systems[ENGINE_WAVE]["bodies"])
    n64 = steps["fp64"]
    p = subprocess.run([ref_exe, "--problem", single, "--dt", DT0, "--epsilon", EPS, "--steps", str(n64),
                        "--sample", "1"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=root)
    # the port at its defaults, REBOUND's: the loop engine, REBOUND's form, 12 passes at most
    rc, cref = run("fp64 the wave as one system, REBOUND's defaults", single, "fp64", "loop", "rebound", n64, 1,
                   max_iter=None)
    if p.returncode or rc:
        verdict(False, "REBOUND's own IAS15 on the wave: ias15_ref rc %d (%s), ias15_cft rc %d"
                % (p.returncode, p.stderr.strip()[-200:], rc))
    else:
        # METHOD.md 5: the binary at that path must be REBOUND's, not a stand-in
        who = re.search(r"^# program=(\S+) impl=(\S+)", p.stdout, re.M)
        verdict(who is not None and who.group(2) == "rebound",
                "%s (built %s) says `%s` in its record's header: REBOUND's own IAS15, not a stand-in at that path"
                % (os.path.basename(ref_exe), stamp(os.path.getmtime(ref_exe)),
                   who.group(0)[2:] if who else "no # program= line"))
        allok &= equivalent("fp64: %s as one system of %d, every body active: REBOUND's own IAS15 (ias15_ref) against "
                            "ias15_cft in REBOUND's form at its defaults" % (systems[ENGINE_WAVE]["name"], N),
                            p.stdout, cref, 6, N)
        # and that record is the gate's own: the all-active member run's, at --max-iter 60
        rc, mem = run("fp64 loop engine rebound form, all active", allact, "fp64", "loop", "rebound", n64, 1, ENGINE_WAVE)
        d = first_difference(parse(cref)[0].get(0, []), parse(mem)[0].get(0, [])) if rc == 0 else (0, -1, "rc", rc)
        verdict(d is None, "the same %d sample lines as the gate's all-active run of %s (--member %d, --max-iter %s), "
                "so REBOUND's IAS15 has the bits `active %d` has%s"
                % (n64 + 1, systems[ENGINE_WAVE]["name"], ENGINE_WAVE, MAX_ITER, na,
                   "" if d is None else ": first difference at sample %s, %s: %s against %s" % (d[0], column(d[1], N), d[2], d[3])))

    # ------------------------------------------------------------ 4. the ensemble
    print("check_ring: the ensemble of %d waves against each wave run alone (--member), fp64, the program engine" % E)
    # the same run as the engines' binary64 program row, so wave ENGINE_WAVE's is shared
    ens_args = ("fp64", "program", "fma", steps["fp64"], 1)
    rc, ens = run("the ensemble", ring, *ens_args)
    erec, ew, esys = parse(ens) if rc == 0 else ({}, {}, {})
    alone = {}
    for w in range(E):
        rc, out = run("wave %d alone" % w, ring, *ens_args, member=w)
        alone[w] = parse(out) if rc == 0 else ({}, {}, {})

    def against_alone(rec, sysw, w):
        """System w of an ensemble record against wave w run alone: None, or
        where they first differ."""
        mrec, mw, _ = alone[w]
        d = first_difference(rec.get(w, []), mrec.get(0, []))
        if d is not None:
            return "first difference at sample %d, %s: ensemble %s, alone %s" % (d[0], column(d[1], N), d[2], d[3])
        for k in SYS_COUNTERS:
            if sysw.get(w, {}).get(k) is None or sysw[w].get(k) != mw.get(k):
                return "trailer %s is %s in the ensemble, %s alone" % (k, sysw.get(w, {}).get(k), mw.get(k))
        if not rec.get(w):
            return "no samples"
        return None

    for w in range(E):
        why = against_alone(erec, esys, w)
        verdict(why is None, "%s: %s" % (systems[w]["name"], why or "%d sample lines identical, its trailer too (%s steps, "
                                         "%s passes a step)" % (len(erec[w]), esys[w].get("steps_done"),
                                                                 esys[w].get("mean_pc_iterations"))))
    # can it fail: one ulp in one particle's x, in a copy
    jw, jb = E - 2, na + 29
    target = names[jw][jb]
    moved = {}

    def nudge(line):
        t = line.split()
        if line.startswith("body ") and t[1] == target:
            x = float.fromhex(t[3])
            t[3] = math.nextafter(x, math.inf).hex()
            moved["x"] = (x, float.fromhex(t[3]))
            return " ".join(t) + "\n"
        return line
    nud = edit_copy(ring, os.path.join(scratch, "ring_one_ulp.txt"), nudge)
    rc, nout = run("the copy with one ulp moved", nud, *ens_args)
    if rc == 0:
        nrec, _, nsys = parse(nout)
        whys = {w: against_alone(nrec, nsys, w) for w in range(E)}
        others = [w for w in range(E) if w != jw and whys[w] is not None]
        ndiff = sum(1 for x, y in zip(nrec.get(jw, []), alone[jw][0].get(0, [])) if x != y)
        verdict(whys[jw] is not None and not others,
                "the comparison can fail: %s's x moved by one ulp (%s to %s) in a copy, %s differs from its run "
                "alone (%s; %d of its %d samples differ) and %s"
                % (target, moved["x"][0].hex(), moved["x"][1].hex(), systems[jw]["name"], whys[jw] or "it did NOT",
                   ndiff, len(nrec.get(jw, [])),
                   "every other wave still matches" if not others else "other waves differ too: %s" % others))
    # MEASURED beside it: each wave's last step against its shortest orbit
    rates = []
    for w in range(E):
        smp = erec.get(w, [])
        if len(smp) < 2:
            continue
        dt = (float.fromhex(smp[-1][2]) - float.fromhex(smp[-2][2])) / (int(smp[-1][1]) - int(smp[-2][1]))
        sw = state(systems[w]["bodies"][0])
        shortest = min((elements(state(b), sw, G)[2], b[0]) for b in systems[w]["bodies"][1:])
        inner = elements(state(systems[w]["bodies"][na]), sw, G)[2]
        rates.append("%s dt %.2f s, %.2f an orbit of %s, its shortest (%.2f of its innermost particle's)"
                     % (systems[w]["name"], dt, shortest[0] / dt, shortest[1], inner / dt))
    print("       MEASURED, not asserted, each wave's last step: %s" % "; ".join(rates))
    # the ring as one system (--one-system), which the program runs as it is
    if pc:
        rc, out = run("one system of %d bodies" % pc["N"], one30, "fp64", "loop", "rebound", 2, 1)
        orec, ow, _ = parse(out) if rc == 0 else ({}, {}, {})
        ma = ACTIVE_LINE.search(out) if rc == 0 else None
        verdict(rc == 0 and len(orec.get(0, [])) == 3 and ow.get("steps_done") == "2" and ma
                and int(ma.group(3)) == pairs_of(pc),
                "one system of %d bodies, fp64, REBOUND's form: %s steps, %s pairs a force evaluation against the %d "
                "its file's N and active give" % (pc["N"], ow.get("steps_done"), ma.group(3) if ma else "no line",
                                                   pairs_of(pc)))

    # ------------------------------------------------------------ 5. the physics
    print("check_ring: the physics")
    # the moons act: their masses set to zero, the same starting bits
    count = {"n": 0}

    def no_moons(line):
        t = line.split()
        if line.startswith("system "):
            count["n"] = 0
        if line.startswith("body "):
            count["n"] += 1
            if 2 <= count["n"] <= na:
                t[2] = (0.0).hex()
                return " ".join(t) + "\n"
        return line
    nomoons = edit_copy(ring, os.path.join(scratch, "ring_no_moons.txt"), no_moons)
    n64 = steps["fp64"]
    rc1, with_m = run("fp64 loop engine rebound form", ring, "fp64", "loop", "rebound", n64, 1, ENGINE_WAVE)
    rc2, without = run("fp64 loop engine rebound form, the moons massless", nomoons, "fp64", "loop", "rebound", n64, 1, ENGINE_WAVE)
    if rc1 == 0 and rc2 == 0:
        a, b = parse(with_m)[0].get(0, []), parse(without)[0].get(0, [])
        last_a, last_b = bodies_at(a[-1], N), bodies_at(b[-1], N)
        dmax, who = max((norm([last_a[i][c] - last_b[i][c] for c in range(3)]), names[ENGINE_WAVE][i])
                        for i in range(na, N))
        reld = dmax / r0[who]
        verdict(first_difference(a, b) is not None and reld > 2.0 ** -40,
                "the moons act: with their masses zeroed the wave differs; after %s steps a particle (%s) is %.3e km "
                "from where the moons put it, %.1e of its radius, against a binary64 rounding of %.1e"
                % (a[-1][1], who, dmax, reld, 2.0 ** -53))
    # the flat-ring counters on a flat copy of the file: each must read zero
    zflat, iflat = 0.0, 0.0
    for s in systems:
        flat = [state(b) for b in s["bodies"]]          # m, x, y, z, vx, vy, vz
        for p in flat:
            p[3], p[6] = 0.0, 0.0
        zflat = max(zflat, largest_z(flat[0][1:], [p[1:] for p in flat[na:]]))
        iflat = max(iflat, max(inclination(p, flat[0])[0] for p in flat[1:na]))
    verdict(zflat == 0 and iflat == 0,
            "the |z| and the moons' inclination counters can fail: on a flat copy of the file (every z and vz zeroed) "
            "the largest |z| reads %r and the largest moon's sin i %r" % (zflat, iflat))
    # the far wave: wholly inside every moon's pericentre by the exclusion's
    # Hill radii, so its own particles set its step; and its steps an orbit
    fw = systems[FAR_WAVE]
    fs = state(fw["bodies"][0])
    peri = min((elements(m, sat, G)[0] * (1 - elements(m, sat, G)[1]) - k_hill * hill[name], name) for name, m in moons)
    fr = max(r0[b[0]] for b in fw["bodies"][na:])
    P, who = min((elements(state(b), fs, G)[2], b[0]) for b in fw["bodies"][1:])
    dt0 = "%.17g" % (0.01 * P / (2 * math.pi))     # GPU-E's first step: 0.01 in units of 1/n
    rc, out = run("the far wave", ring, "fp64", "loop", "rebound", spo_steps, spo_steps // 2, FAR_WAVE, dt0)
    if rc == 0:
        smp = parse(out)[0].get(0, [])
        nst = int(smp[2][1]) - int(smp[1][1])
        span = float.fromhex(smp[2][2]) - float.fromhex(smp[1][2])
        spo = nst * P / span
        verdict(fr < peri[0] and abs(spo / SPO_REF - 1) <= SPO_TOL,
                "the far wave %s (radii to %.0f km, wholly inside %s's pericentre less %g Hill radii, %.0f km): %.2f "
                "steps an orbit of %s, its shortest (%.0f s), over steps %s to %s at epsilon %s, against GPU-E's "
                "MEASURED %.1f: %+.1f%% (the bound is %d%%)"
                % (fw["name"], fr, peri[1], k_hill, peri[0], spo, who, P, smp[1][1], smp[2][1], EPS, SPO_REF,
                   100 * (spo / SPO_REF - 1), round(100 * SPO_TOL)))
    # every run, last: its `# active` line against its own problem's N and
    # active, and every sample against the annulus, Saturn and the plane
    seen, without_line, bad = {}, 0, []
    for rc, out, err, secs, problem, member, key in cache.values():
        if rc:
            continue
        physics(key, out, problem, member)
        ma = ACTIVE_LINE.search(out)
        p = problems[problem]
        if p["active"] is None:
            without_line += 1
            if ma:
                bad.append("%s: an `# active` line, and its problem has no active line" % key)
        else:
            got = (int(ma.group(1)), int(ma.group(2)), int(ma.group(3))) if ma else None
            if got != (p["active"], p["N"], pairs_of(p)):
                bad.append("%s: %s, and its problem's N = %d with active %d give %d" % (
                    key, ma.group(0) if ma else "no `# active` line", p["N"], p["active"], pairs_of(p)))
            else:
                seen[got] = seen.get(got, 0) + 1
    verdict(seen and not bad,
            "the pairs, from every run's `# active` line, each against REBOUND's two loops' n(n-1)/2 + (N-n)n from its "
            "own problem file: %s; and no line in the %d runs with every body active"
            % ("; ".join("%d pairs a force evaluation in %d run%s of `active %d of %d` (%d with every body active)"
                         % (k[2], v, "" if v == 1 else "s", k[0], k[1], k[1] * (k[1] - 1) // 2)
                         for k, v in sorted(seen.items(), key=lambda kv: kv[0][1])),
               without_line) if not bad else "the pairs: %s" % "; ".join(bad[:3]))
    verdict(phys["samples"] > 0 and annulus[0] <= phys["rmin"] and phys["rmax"] < annulus[1] and phys["rmin"] > r_saturn
            and phys["zmin"] > 0,
            "over all %d runs, %d sampled systems: every particle stayed in the annulus [%g, %g) and off Saturn (%g km), "
            "at %.3f to %.3f km, the largest change of a radius %.3f km; every wave's largest |z| at least %.3f km"
            % (len(cache), phys["samples"], annulus[0], annulus[1], r_saturn, phys["rmin"], phys["rmax"],
               phys["excursion"], phys["zmin"]))
    secs = sum(v[3] for v in cache.values())
    print("check_ring: %d runs of ias15_cft, %.0f s of them; %.0f s in all" % (len(cache), secs, time.time() - t_start))
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
