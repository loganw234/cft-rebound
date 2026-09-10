#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Write the two reference problems as exact binary64 bit patterns, so
# that REBOUND, the plain-double reference and the libcft port at every
# format all integrate the SAME initial state - widened exactly, never
# re-rounded.
#
#   kepler   a two-body orbit, G = 1, m0 = 1, m1 = 1/1000, a = 1,
#            e = 1/2, started at pericentre in the barycentric frame.
#            e = 1/2 is dyadic, so the only rounded initial values are
#            the barycentric weights m/(m0+m1) and the speed sqrt(3.003).
#   outer    the outer solar system of REBOUND's own example
#            examples/outer_solar_system/problem.c (Applegate et al.
#            1986 values; Sun + inner planets, Jupiter, Saturn, Uranus,
#            Neptune, Pluto as a test particle), in AU, days and solar
#            masses with G = k^2. Those numbers are measurements and
#            cannot be derived, so they are COPIED BY PROGRAM out of the
#            pinned REBOUND source rather than retyped.
#
# Python floats are IEEE binary64 with the same round-to-nearest-even
# arithmetic as C doubles on x86-64, so the divisions and products
# below produce the bits REBOUND's example would produce itself.

import math
import os
import re
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
outdir = os.path.join(root, "data", "problems")


def hx(x):
    return float(x).hex()


def write_problem(path, name, G, bodies, note):
    with open(path, "w", newline="\n") as f:
        f.write("# cft-rebound problem file, binary64 hex floats\n")
        f.write("# %s\n" % note)
        f.write("name %s\n" % name)
        f.write("G %s\n" % hx(G))
        f.write("N %d\n" % len(bodies))
        for b in bodies:
            f.write("body %s %s %s %s %s %s %s %s\n" % (b[0], hx(b[1]), hx(b[2]), hx(b[3]),
                                                       hx(b[4]), hx(b[5]), hx(b[6]), hx(b[7])))
    print("wrote", path)


def kepler(e=0.5, fname="kepler.txt", pname="kepler"):
    G = 1.0
    m0 = 1.0
    m1 = 1.0 / 1000.0
    a = 1.0
    mu = G * (m0 + m1)
    q = a * (1.0 - e)                      # pericentre distance (exact for e = 1/2)
    v = math.sqrt(mu * (1.0 + e) / q)      # pericentre speed
    mt = m0 + m1
    w0 = m1 / mt
    w1 = m0 / mt
    bodies = [
        ("star",   m0, -q * w0, 0.0, 0.0, 0.0, -v * w0, 0.0),
        ("planet", m1,  q * w1, 0.0, 0.0, 0.0,  v * w1, 0.0),
    ]
    write_problem(os.path.join(outdir, fname), pname, G, bodies,
                  "two-body, G=1, m0=1, m1=1e-3, a=1, e=%s, pericentre, barycentric frame" % ("1/2" if e == 0.5 else repr(e)))


def pythagorean():
    """Burrau's problem (Burrau 1913; Szebehely & Peters 1967): masses 3,
    4, 5 at rest at the vertices of a 3-4-5 right triangle, G = 1. Every
    initial value is an integer, so it is exact in every format, and the
    centre of mass is at the origin already: (3*1 - 4*2 + 5*1, 3*3 - 4 - 5)
    = (0, 0). The system is chaotic with repeated close encounters and
    ends, near t = 70, with the 4-5 binary ejecting body 3."""
    bodies = [
        ("m3", 3.0,  1.0,  3.0, 0.0, 0.0, 0.0, 0.0),
        ("m4", 4.0, -2.0, -1.0, 0.0, 0.0, 0.0, 0.0),
        ("m5", 5.0,  1.0, -1.0, 0.0, 0.0, 0.0, 0.0),
    ]
    write_problem(os.path.join(outdir, "pythagorean.txt"), "pythagorean", 1.0, bodies,
                  "Burrau's Pythagorean three-body problem: m = 3, 4, 5 at (1,3), (-2,-1), (1,-1), at rest, G = 1; all values exact")


def outer():
    src = os.path.join(root, "third_party", "rebound", "examples", "outer_solar_system", "problem.c")
    text = open(src, encoding="utf-8").read()

    def block(name):
        m = re.search(r"double %s\[6\](?:\[3\])?\s*=\s*\{(.*?)\};" % name, text, re.S)
        if not m:
            raise SystemExit("could not find %s in %s" % (name, src))
        return m.group(1)

    def rows(name):
        body = block(name)
        out = []
        for line in body.splitlines():
            line = line.split("//")[0].strip().strip(",").strip()
            if not line:
                continue
            line = line.strip("{}")
            vals = [float(eval(v.strip())) for v in line.split(",") if v.strip()]
            out.append(vals)
        return out

    pos = rows("ss_pos")
    vel = rows("ss_vel")
    mass_body = block("ss_mass")
    masses = []
    for line in mass_body.splitlines():
        line = line.split("//")[0].strip().strip(",").strip()
        if line:
            masses.append(float(eval(line)))
    if not (len(pos) == 6 and len(vel) == 6 and len(masses) == 6):
        raise SystemExit("parsed %d/%d/%d rows, expected 6" % (len(pos), len(vel), len(masses)))
    mk = re.search(r"const double k = ([0-9.eE+-]+);", text)
    if not mk:
        raise SystemExit("could not find the Gaussian constant k")
    k = float(mk.group(1))
    G = k * k
    names = ["sun", "jupiter", "saturn", "uranus", "neptune", "pluto"]
    bodies = []
    for i in range(6):
        bodies.append((names[i], masses[i], pos[i][0], pos[i][1], pos[i][2],
                       vel[i][0], vel[i][1], vel[i][2]))
    write_problem(os.path.join(outdir, "outer.txt"), "outer", G, bodies,
                  "REBOUND examples/outer_solar_system/problem.c values, copied by program; AU, days, Msun; G = k*k, k = %r" % k)


def main():
    os.makedirs(outdir, exist_ok=True)
    kepler()
    kepler(0.9, "kepler_e09.txt", "kepler_e09")
    kepler(0.99, "kepler_e099.txt", "kepler_e099")
    pythagorean()
    outer()
    return 0


if __name__ == "__main__":
    sys.exit(main())
