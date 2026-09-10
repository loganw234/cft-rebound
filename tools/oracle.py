#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Score a record file against the truth, from the exact bits.
#
# Every value in a record is an exact hex float, so the analysis never
# adds rounding of its own: positions, velocities and the exact elapsed
# time (t_hi + t_lo) are converted to rationals and the invariants are
# evaluated with mpmath at 60 digits.
#
#   kepler  the two-body problem has a closed form. From the exact
#           initial state the script derives mu, a, e, the period T and
#           the initial eccentric anomaly, then for every sample solves
#           Kepler's equation at the exact elapsed time and reports the
#           position error of the planet against the closed form, the
#           energy error and the angular-momentum error, all relative.
#   outer   no closed form: energy and angular momentum only.
#
#   python tools/oracle.py RECORD [--csv OUT]
#
# Prints one line per sample: orbits (t/T for kepler, years for outer),
# then relative errors; and a summary line with the maxima.

import argparse
import sys
from fractions import Fraction

import mpmath
from mpmath import mp, mpf

mp.dps = 60


def hex_to_fraction(s):
    t = s.strip()
    neg = t.startswith("-")
    if neg or t.startswith("+"):
        t = t[1:]
    if t in ("inf", "nan"):
        raise ValueError("non-finite value in record: %s" % s)
    body = t[2:]
    mant, _, exp = body.partition("p")
    ip, _, fp = mant.partition(".")
    digits = (ip + fp) or "0"
    val = Fraction(int(digits, 16), 16 ** len(fp)) * (Fraction(2) ** int(exp))
    return -val if neg else val


def F2m(x):
    return mpf(x.numerator) / mpf(x.denominator)


def parse(path):
    header = {}
    masses = []
    samples = []
    tail = {}
    with open(path) as f:
        for line in f:
            if line.startswith("# program="):
                for kv in line[2:].split():
                    k, _, v = kv.partition("=")
                    header[k] = v
            elif line.startswith("# dt0="):
                for kv in line[2:].split():
                    k, _, v = kv.partition("=")
                    if v:
                        header[k] = v
                # the cft program prints "dt0= 0x..." with a space; recover
                toks = line[2:].split()
                keys = [t for t in toks if t.endswith("=")]
                vals = [t for t in toks if not t.endswith("=") and "=" not in t]
                if keys and len(keys) == len(vals):
                    for k, v in zip(keys, vals):
                        header[k[:-1]] = v
            elif line.startswith("# body"):
                toks = line.split()
                mtok = [t for t in toks if t.startswith("m=") or t.startswith("0x") or t.startswith("-0x")]
                m = mtok[-1]
                if m.startswith("m="):
                    m = m[2:]
                masses.append(hex_to_fraction(m))
            elif line.startswith("sample "):
                toks = line.split()
                k, step = int(toks[1]), int(toks[2])
                body = toks[3:]
                exact_t = None
                if "|" in body:
                    i = body.index("|")
                    thi, tlo = body[i + 1], body[i + 2]
                    exact_t = hex_to_fraction(thi) + hex_to_fraction(tlo)
                    body = body[:i]
                vals = [hex_to_fraction(x) for x in body]
                samples.append({"k": k, "step": step, "t": vals[0], "dt_next": vals[1], "dt_last": vals[2],
                                "E_prog": vals[3], "state": vals[4:], "t_exact": exact_t})
            elif line.startswith("# steps_done="):
                for kv in line[2:].split():
                    k, _, v = kv.partition("=")
                    tail[k] = v
    return header, masses, samples, tail


def invariants(G, masses, state):
    N = len(masses)
    ekin = mpf(0)
    epot = mpf(0)
    L = [mpf(0)] * 3
    pos = []
    vel = []
    for i in range(N):
        x = [F2m(state[6 * i + c]) for c in range(3)]
        v = [F2m(state[6 * i + 3 + c]) for c in range(3)]
        pos.append(x)
        vel.append(v)
        m = F2m(masses[i])
        ekin += m * (v[0] ** 2 + v[1] ** 2 + v[2] ** 2) / 2
        L[0] += m * (x[1] * v[2] - x[2] * v[1])
        L[1] += m * (x[2] * v[0] - x[0] * v[2])
        L[2] += m * (x[0] * v[1] - x[1] * v[0])
    for i in range(N):
        for j in range(i + 1, N):
            d = mpmath.sqrt(sum((pos[i][c] - pos[j][c]) ** 2 for c in range(3)))
            epot -= G * F2m(masses[i]) * F2m(masses[j]) / d
    return ekin + epot, L, pos, vel


class Kepler:
    """Closed-form two-body motion from the exact initial state."""

    def __init__(self, G, masses, state):
        m0, m1 = F2m(masses[0]), F2m(masses[1])
        self.M = m0 + m1
        self.w0, self.w1 = m0 / self.M, m1 / self.M
        x0 = [F2m(state[c]) for c in range(3)]
        v0 = [F2m(state[3 + c]) for c in range(3)]
        x1 = [F2m(state[6 + c]) for c in range(3)]
        v1 = [F2m(state[9 + c]) for c in range(3)]
        self.Xcm = [self.w0 * x0[c] + self.w1 * x1[c] for c in range(3)]
        self.Vcm = [self.w0 * v0[c] + self.w1 * v1[c] for c in range(3)]
        self.r0 = [x1[c] - x0[c] for c in range(3)]
        self.v0 = [v1[c] - v0[c] for c in range(3)]
        self.mu = G * self.M
        r = mpmath.sqrt(sum(q * q for q in self.r0))
        v2 = sum(q * q for q in self.v0)
        eps = v2 / 2 - self.mu / r
        self.a = -self.mu / (2 * eps)
        rv = sum(self.r0[c] * self.v0[c] for c in range(3))
        ecosE = 1 - r / self.a
        esinE = rv / mpmath.sqrt(self.mu * self.a)
        self.e = mpmath.sqrt(ecosE ** 2 + esinE ** 2)
        self.E0 = mpmath.atan2(esinE, ecosE)
        self.n = mpmath.sqrt(self.mu / self.a ** 3)
        self.T = 2 * mpmath.pi / self.n
        self.rmag0 = r

    def at(self, t):
        """Relative position/velocity at time t via the f and g functions."""
        M = self.n * t   # mean anomaly advanced from E0: solve dE - e(sin(E0+dE) - sin E0) = n t
        e, E0 = self.e, self.E0
        dE = M
        for _ in range(100):
            f = dE - e * (mpmath.sin(E0 + dE) - mpmath.sin(E0)) - M
            fp = 1 - e * mpmath.cos(E0 + dE)
            step = f / fp
            dE -= step
            if abs(step) < mpf(10) ** (-(mp.dps - 5)):
                break
        a, r0 = self.a, self.rmag0
        f = 1 - (a / r0) * (1 - mpmath.cos(dE))
        g = t - (dE - mpmath.sin(dE)) / self.n
        r = [f * self.r0[c] + g * self.v0[c] for c in range(3)]
        rmag = mpmath.sqrt(sum(q * q for q in r))
        fdot = -mpmath.sqrt(self.mu * a) / (r0 * rmag) * mpmath.sin(dE)
        gdot = 1 - (a / rmag) * (1 - mpmath.cos(dE))
        v = [fdot * self.r0[c] + gdot * self.v0[c] for c in range(3)]
        return r, v

    def bodies_at(self, t):
        r, v = self.at(t)
        X = [self.Xcm[c] + self.Vcm[c] * t for c in range(3)]
        x0 = [X[c] - self.w1 * r[c] for c in range(3)]
        x1 = [X[c] + self.w0 * r[c] for c in range(3)]
        v0 = [self.Vcm[c] - self.w1 * v[c] for c in range(3)]
        v1 = [self.Vcm[c] + self.w0 * v[c] for c in range(3)]
        return x0, v0, x1, v1


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("record")
    ap.add_argument("--csv")
    ap.add_argument("--quiet", action="store_true")
    args = ap.parse_args()
    header, masses, samples, tail = parse(args.record)
    G = F2m(hex_to_fraction(header["G"]))
    problem = header.get("problem", "?")
    fmt = header.get("format", "?")
    E0, L0, _, _ = invariants(G, masses, samples[0]["state"])
    L0mag = mpmath.sqrt(sum(q * q for q in L0))
    kep = Kepler(G, masses, samples[0]["state"]) if problem == "kepler" else None
    fixed = header.get("epsilon", "") in ("0x0p+0", "0x0.0p+0", "0")
    dt0 = hex_to_fraction(header["dt0"]) if "dt0" in header else None
    rows = []
    out = open(args.csv, "w", newline="\n") if args.csv else None
    if out:
        out.write("sample,step,orbits,t,rel_energy_err,rel_angmom_err,rel_pos_err,rel_vel_err\n")
    maxE = mpf(0); maxL = mpf(0); maxP = mpf(0); lastP = mpf(0); lastE = mpf(0)
    for s in samples:
        if s["t_exact"] is not None:
            t = F2m(s["t_exact"])
        elif fixed and dt0 is not None:
            t = F2m(dt0 * s["step"])
        else:
            t = F2m(s["t"])
        E, L, pos, vel = invariants(G, masses, s["state"])
        dE = abs((E - E0) / E0)
        Lmag = mpmath.sqrt(sum(q * q for q in L))
        dL = abs((Lmag - L0mag) / L0mag)
        dP = mpf(0); dV = mpf(0)
        if kep:
            x0, v0, x1, v1 = kep.bodies_at(t)
            dP = mpmath.sqrt(sum((pos[1][c] - x1[c]) ** 2 for c in range(3))) / kep.a
            vscale = mpmath.sqrt(kep.mu / kep.a)
            dV = mpmath.sqrt(sum((vel[1][c] - v1[c]) ** 2 for c in range(3))) / vscale
            orbits = t / kep.T
        else:
            orbits = t / mpf(365.25)   # years, for the outer solar system in days
        maxE = max(maxE, dE); maxL = max(maxL, dL); maxP = max(maxP, dP); lastP = dP; lastE = dE
        rows.append((s["k"], s["step"], orbits, t, dE, dL, dP, dV))
        if out:
            out.write("%d,%d,%s,%s,%s,%s,%s,%s\n" % (s["k"], s["step"], mpmath.nstr(orbits, 12), mpmath.nstr(t, 20),
                                                    mpmath.nstr(dE, 6), mpmath.nstr(dL, 6), mpmath.nstr(dP, 6), mpmath.nstr(dV, 6)))
        if not args.quiet:
            print("sample %4d step %8d %s %12s  dE/E %s  dL/L %s  dx/a %s" %
                  (s["k"], s["step"], "orbits" if kep else "years ", mpmath.nstr(orbits, 8),
                   mpmath.nstr(dE, 4), mpmath.nstr(dL, 4), mpmath.nstr(dP, 4) if kep else "-"))
    if out:
        out.close()
    last = rows[-1]
    print("SUMMARY problem=%s format=%s cs=%s steps=%d %s=%s max_dE=%s last_dE=%s max_dL=%s max_dx=%s last_dx=%s rejected=%s max_exceeded=%s mean_pc=%s seconds=%s" %
          (problem, fmt, header.get("cs", "-"), last[1], "orbits" if kep else "years", mpmath.nstr(last[2], 8),
           mpmath.nstr(maxE, 4), mpmath.nstr(lastE, 4), mpmath.nstr(maxL, 4),
           mpmath.nstr(maxP, 4) if kep else "-", mpmath.nstr(lastP, 4) if kep else "-",
           tail.get("steps_rejected", "-"), tail.get("iterations_max_exceeded", "-"),
           tail.get("mean_pc_iterations", "-"), tail.get("seconds", "-")))
    if kep and not args.quiet:
        print("kepler: a=%s e=%s T=%s mu=%s" % (mpmath.nstr(kep.a, 20), mpmath.nstr(kep.e, 20), mpmath.nstr(kep.T, 20), mpmath.nstr(kep.mu, 20)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
