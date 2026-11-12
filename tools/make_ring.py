#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# A Saturn ring as a problem file (ROUND4.md, P2): Saturn with its zonal
# line, six moons and N test particles, as one system, or as an ensemble
# of 64-body waves - the 7 massive bodies, then 57 particles - which is
# cft-fp256's study GPU-E's shape a2, one 64-lane wave a system.
#
#   python tools/make_ring.py OUT.txt --particles N [--seed S] [--one-system]
#          [--no-zonal] [--rmin KM] [--rmax KM] [--exclude-hill K] [--name NAME]
#
# The same arguments always give the same bytes; the output's path is not
# one of them, and nothing else (a time, a host, a version) is written.
#
# UNITS. km and s, with G = 1, so every mass is a GM in km^3 s^-2.
#
# THE MEASURED CONSTANTS, every one READ in cft-fp256's
# docs/studies/GPU-E-mi300x-ring.md section 4, which names its sources;
# no web page was read for this file (ROUND4.md, decision 5):
#   Saturn  GM 37.931e6, equatorial radius (1 bar) 60,268 km
#           - NASA's Saturn fact sheet;
#           J2 16,290.573e-6 and J4 -935.314e-6 - Cassini (Iess et al.
#           2019), in Fortney et al., arXiv:2304.09215, Table 4.1.
#           The study's J4 table states no reference radius and BELIEVES
#           60,330 km. That pair - those J2 and J4 with R = 60,330 - is
#           the zonal line, and the file's head says so. The 60,268 is
#           used only as "Saturn's surface", for tools/check_ring.py.
#   moons   GM, a and e - JPL SSD's physical parameters and mean elements
#           (Laplace-plane elements, "not intended for ephemeris
#           computation"):
#             Mimas       2.50349   186,000 km  0.020
#             Pan         0.00028   133,600 km  0
#             Prometheus  0.01071   139,400 km  0.002
#             Pandora     0.00926   141,700 km  0.004
#             Epimetheus  0.03514   151,400 km  0.020
#             Janus       0.12662   151,500 km  0.007
#   annulus 130,000 to 142,000 km - the study's recommended ring (its
#           answer 1): the Encke gap with Pan, the Keeler gap, the A
#           ring's edge, the F ring between Prometheus and Pandora.
# Each decimal is rounded once to binary64 (Python's float() is correctly
# rounded), and those binary64 values are the ones written and the ones
# everything below is computed from: the speeds are the field's that the
# file itself describes.
#
# THE SYNTHETIC CONSTANTS (ROUND4.md, decision 5: an epoch's would need
# JPL's tables read from the web). The moons' inclinations, nodes,
# pericentres and true anomalies are MOON_ANGLES below: fixed, small, not
# zero (a flat ring would leave every z term of J2 and J4 exactly zero,
# ROUND4.md "the trap"), and NOT any epoch's or any table's. Epimetheus and
# Janus share an orbit 100 km apart, so their true longitudes are put
# 180 degrees apart: phases near each other would start a close encounter.
#
# THE PARTICLES. Each draw takes four numbers from Python's
# random.Random(seed) (Mersenne Twister; random() is guaranteed to give
# the same sequence for the same seed across versions and hosts), in
# this order: the radius, uniform in [rmin, rmax); the argument of
# latitude and the node, uniform in [0, 2 pi); the inclination, uniform
# in [INC_LO, INC_HI) radians - small, and never zero.
#   - A draw whose start lies within EXCLUDE_HILL Hill radii of a moon's
#     start is refused and drawn again, because a particle started there
#     would make its wave's steps collapse. The Hill radius is
#     a (GM_moon / (3 GM_Saturn))^(1/3), computed from the constants above
#     (the study's column: Pan 18 km, Prometheus 63, Pandora 61). How many
#     draws were refused, and near which moon, is in the file's head.
#   - The orbit is circular in the field the file describes. With the
#     zonal line, its speed is the circular speed in Saturn's equatorial
#     plane with J2 and J4 at the particle's radius,
#         v^2 = (GM / r) (1 + 3/2 J2 (R/r)^2 - 15/8 J4 (R/r)^4),
#     from Phi = -(GM/r) (1 - J2 (R/r)^2 P2(u) - J4 (R/r)^4 P4(u)) at u = 0
#     (ROUND4.md, P1). With --no-zonal there is no oblate field in the
#     problem, so the speed is the point mass's, sqrt(GM / r): a particle
#     given the oblate speed in a point mass's field is on an e ~ 0.005
#     orbit and swings out by about 1,400 km.
#   - The particles are sorted by radius, then cut into waves of 57.
# THE MOONS are Saturn-centric Kepler orbits with mu = G (M + m), in
# REBOUND's own form (third_party/rebound/src/tools.c, reb_particle_from_
# orbit: Murray & Dermott eq. 2.122 and 2.36). With the zonal line each
# velocity is scaled by the same oblate factor at the moon's start radius,
# so a circular one (Pan) stays circular in the oblate field; without it
# the factor is 1.
# THE FRAME is the barycentre of the 7 massive bodies, position and
# velocity, subtracted from every body. In an ensemble every wave holds
# the same Saturn and moons, bit for bit, so `active 7` is right for every
# system and every wave has the same frame.
#
# EXACTNESS. Everything is computed in mpmath at PREC bits from the
# binary64 constants and the binary64 draws, and every written value is
# rounded once to binary64, to nearest, by b64() below, then written as
# an exact hex float. Python's math.sin and math.cos are not guaranteed to
# be correctly rounded and are not used. That the file is then the same on
# every host is BELIEVED (mpmath is integer arithmetic in pure Python or
# gmpy); tools/check_ring.py holds the two mpmath backends of one host to
# the same bytes, and the check across hosts waits for a WSL slot.

import argparse
import math
import os
import random
import sys

import mpmath
from mpmath.libmp import normalize, round_nearest

PREC = 160                      # mpmath's working precision, bits

# Read in GPU-E section 4 (see the head). Decimal strings, rounded once.
SATURN_NAME = "saturn"
SATURN_GM = "37.931e6"          # km^3 s^-2, NASA's Saturn fact sheet
SATURN_REQ_KM = "60268"         # equatorial radius, 1 bar, the same sheet
ZONAL_R_KM = "60330"            # BELIEVED by the study (its J4 table states none)
ZONAL_J2 = "16290.573e-6"       # Cassini, Iess et al. 2019
ZONAL_J4 = "-935.314e-6"        # Cassini, Iess et al. 2019
RMIN_KM, RMAX_KM = "130000", "142000"   # the study's annulus
#           name          GM, km^3 s^-2   a, km      e
MOONS = (("mimas",       "2.50349",     "186000",  "0.020"),
         ("pan",         "0.00028",     "133600",  "0"),
         ("prometheus",  "0.01071",     "139400",  "0.002"),
         ("pandora",     "0.00926",     "141700",  "0.004"),
         ("epimetheus",  "0.03514",     "151400",  "0.020"),
         ("janus",       "0.12662",     "151500",  "0.007"))
# SYNTHETIC (decision 5), degrees: inclination, node, argument of
# pericentre, true anomaly. True longitudes (node + pericentre + anomaly):
# Mimas 90, Pan 150, Prometheus 270, Pandora 330, Epimetheus 45, Janus 225.
MOON_ANGLES = {"mimas":      ("1.0",  "30",  "60", "0"),
               "pan":        ("0.05", "150", "0",  "0"),
               "prometheus": ("0.05", "210", "0",  "60"),
               "pandora":    ("0.05", "300", "30", "0"),
               "epimetheus": ("0.2",  "0",   "0",  "45"),
               "janus":      ("0.2",  "180", "0",  "45")}
INC_LO, INC_HI = "1e-4", "1e-3"  # the particles' inclinations, radians: synthetic
EXCLUDE_HILL = "5"               # Hill radii: "a few", ROUND4.md P2
PER_WAVE = 57                    # particles a wave: 64 lanes less Saturn and the moons
MASSIVE = 1 + len(MOONS)         # bodies 0..6 of every system


def b64(x):
    """x, an mpf, rounded once to the nearest binary64 (ties to even)."""
    sign, man, exp, bc = x._mpf_
    if not man:
        if x == 0:
            return -0.0 if sign else 0.0
        raise SystemExit("make_ring: a value is not finite: %r" % x)
    sign, man, exp, bc = normalize(sign, man, exp, bc, 53, round_nearest)
    if not -1022 <= exp + bc - 1 <= 1023:
        raise SystemExit("make_ring: %r is outside binary64's normal range" % x)
    v = math.ldexp(man, exp)
    v = -v if sign else v
    if mp.mpf(v) != mp.mpf((sign, man, exp, bc)):
        raise SystemExit("make_ring: %r did not round exactly" % x)
    return v


mp = mpmath.MPContext()
mp.prec = PREC


def f64(text):
    """A decimal constant as the binary64 the file holds, and as an mpf."""
    v = float(text)
    return v, mp.mpf(v)


def rot(r, v0, e, inc, node, peri, f):
    """REBOUND's reb_particle_from_orbit (tools.c), relative to the
    primary: position and velocity in mpmath."""
    cO, sO = mp.cos(node), mp.sin(node)
    co, so = mp.cos(peri), mp.sin(peri)
    cf, sf = mp.cos(f), mp.sin(f)
    ci, si = mp.cos(inc), mp.sin(inc)
    x = r * (cO * (co * cf - so * sf) - sO * (so * cf + co * sf) * ci)
    y = r * (sO * (co * cf - so * sf) + cO * (so * cf + co * sf) * ci)
    z = r * (so * cf + co * sf) * si
    vx = v0 * ((e + cf) * (-ci * co * sO - cO * so) - sf * (co * cO - ci * so * sO))
    vy = v0 * ((e + cf) * (ci * co * cO - sO * so) - sf * (co * sO + ci * so * cO))
    vz = v0 * ((e + cf) * co * si - sf * si * so)
    return [x, y, z, vx, vy, vz]


class Ring:
    def __init__(self, zonal):
        self.G = mp.mpf(1)
        self.gm_s_f, self.gm_s = f64(SATURN_GM)
        self.zr_f, self.zr = f64(ZONAL_R_KM)
        self.j2_f, self.j2 = f64(ZONAL_J2)
        self.j4_f, self.j4 = f64(ZONAL_J4)
        self.zonal = zonal
        self.moons = []
        deg = mp.pi / 180
        for name, gm, a, e in MOONS:
            gm_f, gm_m = f64(gm)
            _, a_m = f64(a)
            _, e_m = f64(e)
            inc, node, peri, f = (mp.mpf(float(t)) * deg for t in MOON_ANGLES[name])
            r = a_m * (1 - e_m * e_m) / (1 + e_m * mp.cos(f))
            v0 = mp.sqrt(self.G * (gm_m + self.gm_s) / a_m / (1 - e_m * e_m)) * self.oblate(r)
            hill = a_m * mp.cbrt(gm_m / (3 * self.gm_s))
            self.moons.append({"name": name, "gm_f": gm_f, "gm": gm_m, "a": a_m, "e": e_m, "hill": hill,
                               "state": rot(r, v0, e_m, inc, node, peri, f)})

    def oblate(self, r):
        """The circular speed's factor over the point mass's at r in the
        equatorial plane: 1 without the zonal line."""
        if not self.zonal:
            return mp.mpf(1)
        q2 = (self.zr / r) ** 2
        return mp.sqrt(1 + mp.mpf(3) / 2 * self.j2 * q2 - mp.mpf(15) / 8 * self.j4 * q2 * q2)

    def particle(self, r, lat, node, inc):
        v0 = mp.sqrt(self.G * self.gm_s / r) * self.oblate(r)
        return rot(r, v0, mp.mpf(0), inc, node, mp.mpf(0), lat)


def draw(ring, n, seed, rmin, rmax, k_hill, place):
    """n particles, each as (u0, place(state)) - u0 the draw the radius came
    from, which orders them as the radius does - sorted by radius; the
    refused draws, by moon; the draws in all; and the largest |z| from
    Saturn."""
    rng = random.Random(seed)
    two_pi = 2 * mp.pi
    inc_lo, inc_hi = mp.mpf(float(INC_LO)), mp.mpf(float(INC_HI))
    lo, hi = mp.mpf(rmin), mp.mpf(rmax)
    for m in ring.moons:
        m["r0"] = mp.sqrt(sum(m["state"][c] ** 2 for c in range(3)))
        m["reach"] = k_hill * m["hill"]
    kept, refused, draws, zmax = [], {m["name"]: 0 for m in ring.moons}, 0, mp.mpf(0)
    while len(kept) < n:
        u = [rng.random() for _ in range(4)]
        draws += 1
        r = lo + (hi - lo) * mp.mpf(u[0])
        inc = inc_lo + (inc_hi - inc_lo) * mp.mpf(u[3])
        s = ring.particle(r, two_pi * mp.mpf(u[1]), two_pi * mp.mpf(u[2]), inc)
        near = None
        for m in ring.moons:
            # |p - m| >= | |p| - |m| |, so a moon whose orbit's radius is
            # out of reach cannot refuse the draw: skip the 3-D distance
            if abs(r - m["r0"]) >= m["reach"]:
                continue
            if sum((s[c] - m["state"][c]) ** 2 for c in range(3)) < m["reach"] ** 2:
                near = m["name"]
                break
        if near:
            refused[near] += 1
            continue
        zmax = max(zmax, abs(s[2]))
        kept.append((u[0], place(s)))
    kept.sort(key=lambda t: t[0])
    return kept, refused, draws, zmax


def main():
    ap = argparse.ArgumentParser(description="Write a Saturn ring problem (ROUND4.md, P2).")
    ap.add_argument("out")
    ap.add_argument("--particles", type=int, required=True, help="test particles in all")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--one-system", action="store_true", help="one system of 7 + N bodies, not waves of 64")
    ap.add_argument("--no-zonal", action="store_true",
                    help="no zonal line, and the point mass's speeds (the line is refused until round 4's P1)")
    ap.add_argument("--rmin", default=RMIN_KM, help="the annulus's inner radius, km (default %s)" % RMIN_KM)
    ap.add_argument("--rmax", default=RMAX_KM, help="the annulus's outer radius, km (default %s)" % RMAX_KM)
    ap.add_argument("--exclude-hill", default=EXCLUDE_HILL,
                    help="refuse a start within this many Hill radii of a moon (default %s)" % EXCLUDE_HILL)
    ap.add_argument("--name", default=None)
    args = ap.parse_args()
    n = args.particles
    if n < 1:
        raise SystemExit("make_ring: --particles must be at least 1")
    if not args.one_system and n % PER_WAVE:
        raise SystemExit("make_ring: %d particles do not fill waves of %d: every system of an ensemble has "
                         "N bodies, so the count must be a multiple of %d" % (n, PER_WAVE, PER_WAVE))
    rmin, rmax, k_hill = float(args.rmin), float(args.rmax), float(args.exclude_hill)
    if not 0 < rmin < rmax:
        raise SystemExit("make_ring: the annulus %r to %r km is empty" % (rmin, rmax))
    if not k_hill >= 0:
        raise SystemExit("make_ring: --exclude-hill must be at least 0")
    zonal = not args.no_zonal
    ring = Ring(zonal)

    # the massive bodies, Saturn at the Saturn-centric origin, then the
    # frame: their barycentre, which no particle moves (they are massless)
    massive = [(SATURN_NAME, ring.gm_s_f, ring.gm_s, [mp.mpf(0)] * 6)]
    massive += [(m["name"], m["gm_f"], m["gm"], m["state"]) for m in ring.moons]
    mtot = sum(b[2] for b in massive)
    centre = [sum(b[2] * b[3][c] for b in massive) / mtot for c in range(6)]

    def place(mass_f, state):
        """A body's mass and its state in the frame, rounded once, as hex."""
        return " ".join(float(v).hex() for v in [mass_f] + [b64(state[c] - centre[c]) for c in range(6)])

    head_massive = ["body %s %s\n" % (nm, place(mf, st)) for nm, mf, _, st in massive]
    kept, refused, draws, zmax = draw(ring, n, args.seed, rmin, rmax, mp.mpf(k_hill), lambda s: place(0.0, s))
    plines = ["body p%07d %s\n" % (j, vals) for j, (_, vals) in enumerate(kept)]
    r_lo = mp.mpf(rmin) + (mp.mpf(rmax) - mp.mpf(rmin)) * mp.mpf(kept[0][0])
    r_hi = mp.mpf(rmin) + (mp.mpf(rmax) - mp.mpf(rmin)) * mp.mpf(kept[-1][0])
    nwave = 1 if args.one_system else n // PER_WAVE
    per_sys = MASSIVE + (n if args.one_system else PER_WAVE)
    name = args.name or "ring_n%d_s%d%s%s" % (n, args.seed, "_one" if args.one_system else "", "" if zonal else "_nozonal")
    argv = ["--particles", str(n), "--seed", str(args.seed)] + (["--one-system"] if args.one_system else []) + \
           (["--no-zonal"] if args.no_zonal else []) + ["--rmin", args.rmin, "--rmax", args.rmax,
                                                       "--exclude-hill", args.exclude_hill]
    with open(args.out, "w", newline="\n") as f:
        f.write("# cft-rebound problem file, binary64 hex floats\n")
        f.write("# GENERATED by tools/make_ring.py %s\n" % " ".join(argv))
        f.write("# A Saturn ring (ROUND4.md, P2): Saturn, six moons and %d test particles, %s.\n"
                % (n, "one system" if args.one_system else "%d waves of %d bodies, %d massive then %d particles"
                   % (nwave, per_sys, MASSIVE, PER_WAVE)))
        f.write("# Units km and s, G = 1, every mass a GM in km^3 s^-2. Values READ in cft-fp256's\n"
                "# docs/studies/GPU-E-mi300x-ring.md section 4 (no web page read here): Saturn's GM %s and\n"
                "# equatorial radius %s km (NASA's fact sheet); J2 %s and J4 %s (Cassini, Iess et al. 2019,\n"
                "# in Fortney et al. arXiv:2304.09215 Table 4.1); the moons' GM, a and e (JPL SSD).\n"
                % (SATURN_GM, SATURN_REQ_KM, ZONAL_J2, ZONAL_J4))
        f.write("# The zonal line's reference radius, %s km, is the study's BELIEF: its J4 table states none.\n"
                % ZONAL_R_KM)
        if zonal:
            f.write("# zonal: body 0 is oblate with R = %s km, J2 = %s, J4 = %s; every orbit's speed is the\n"
                    "# circular speed in the oblate field at its start radius, in the equatorial plane.\n"
                    % (ZONAL_R_KM, ZONAL_J2, ZONAL_J4))
        else:
            f.write("# --no-zonal: no zonal line, so the field is the point mass's and so are the speeds.\n")
        f.write("# SYNTHETIC (ROUND4.md decision 5), not an epoch's: the moons' inclination, node, pericentre\n"
                "# and true anomaly, degrees: %s;\n" % "; ".join("%s %s" % (k, " ".join(MOON_ANGLES[k])) for k, *_ in MOONS))
        f.write("# the particles' radius uniform in [%s, %s) km, argument of latitude and node uniform, inclination\n"
                "# uniform in [%s, %s) rad, from random.Random(%d); circular orbits about Saturn.\n"
                % (args.rmin, args.rmax, INC_LO, INC_HI, args.seed))
        f.write("# Hill radii, km: %s.\n" % ", ".join("%s %.1f" % (m["name"], float(m["hill"])) for m in ring.moons))
        f.write("# A start within %s Hill radii of a moon's start was refused and drawn again: %d of %d draws (%s).\n"
                % (args.exclude_hill, draws - n, draws, ", ".join("%s %d" % kv for kv in refused.items())))
        f.write("# Frame: the barycentre of the %d massive bodies, applied to every body. Particles sorted by radius%s.\n"
                % (MASSIVE, "" if args.one_system else ", then cut into waves"))
        # key=value lines for tools/check_ring.py; ias15_cft skips every # line
        f.write("# ring annulus_km %r %r\n" % (rmin, rmax))
        f.write("# ring saturn_radius_km %r\n" % float(SATURN_REQ_KM))
        f.write("# ring exclusion hill=%r refused=%d draws=%d\n" % (k_hill, draws - n, draws))
        f.write("# ring hill_km %s\n" % " ".join("%s=%r" % (m["name"], float(m["hill"])) for m in ring.moons))
        f.write("name %s\n" % name)
        f.write("G %s\n" % (1.0).hex())
        f.write("N %d\n" % per_sys)
        if not args.one_system:
            f.write("E %d\n" % nwave)
        f.write("active %d\n" % MASSIVE)
        if zonal:
            f.write("zonal 0 %s %s %s\n" % (ring.zr_f.hex(), ring.j2_f.hex(), ring.j4_f.hex()))
        if args.one_system:
            f.writelines(head_massive + plines)
        else:
            for w in range(nwave):
                f.write("system w%05d\n" % w)
                f.writelines(head_massive + plines[w * PER_WAVE:(w + 1) * PER_WAVE])
    print("make_ring: wrote %s: %s, %d particles from %.3f to %.3f km, largest |z| %.3f km; refused %d of %d draws "
          "within %s Hill radii of a moon%s (mpmath %s, %s backend)"
          % (args.out, "one system of %d bodies" % per_sys if args.one_system else "%d waves of %d bodies" % (nwave, per_sys),
             n, float(r_lo), float(r_hi), float(zmax), draws - n, draws, args.exclude_hill,
             "" if zonal else ", no zonal line", mpmath.__version__, mpmath.libmp.BACKEND), file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
