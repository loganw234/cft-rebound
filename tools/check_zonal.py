#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of the zonal term: J2 and J4 of one body of every system, in
# libcft's format (ROUND4.md, P1; src/ias15_zonal.inc).
#
#   python tools/check_zonal.py [--build build] [--quick]
#   (`make check` runs it whole; check-quick and check-light pass --quick,
#   which integrates fewer orbits at binary128 and binary256. Every case
#   integrates to a count of steps; none is a check_dropin case, so
#   --light's scaling never reaches one.)
#
# What it asserts, and what lets each assertion fail:
#   1. ARITHMETIC. `ias15_cft --zonal-selftest --format F`, at all three
#      formats, evaluates the term and the potential at its fixed points and
#      traces every library operation it makes (operands and result as bits,
#      and each element's flags) and every element it moves between them.
#      The golden model (third_party/cft-fp256/python/cft_golden) recomputes
#      each traced element: bits and flags must be its. Every operand must be
#      a declared input or what the trace last wrote there, and the results
#      the selftest prints must be the trace's: so no value reaches the
#      result by a path the trace does not show. A value corrupted by one ulp
#      is named by its operation. Beside it, the physics: at every point off
#      the plane each component of the term is nonzero, and on it z's is +0.
#   2. FORMULA. The same points against minus the gradient of the potential,
#      Phi = G M [J2 R^2 r^-3 P2(u) + J4 R^4 r^-5 P4(u)], u = dz/r, taken by
#      mpmath.diff at 200 digits from mpmath's own Legendre polynomials,
#      never from the force: each body's acceleration, body k's as the sum
#      of every interacting body's reaction (m_j/M times its term, under
#      testparticle_type 0), and each system's potential.
#      The bound is derived per component from the trace's own dataflow
#      (Higham, Accuracy and Stability, 3.1): a rounded result is its exact
#      operation times (1 + d), |d| <= u = 2^-p; so a value reached through k
#      roundings, counted on the path with multiplicity (one more for an add
#      over its deeper operand, one more for a product or quotient over both
#      operands' counts), is within gamma_k = k u / (1 - k u) of its exact
#      value times the same dataflow evaluated in absolute values (each sum
#      as the sum of its operands' absolute evaluations, or, of two exact
#      operands, as its own exact value, one monomial; each divisor, here
#      a sum of squares, its root, or an exact constant, by its own exact
#      value, which the gate checks it is). So
#          |computed - reference| <= gamma_k * |abs evaluation| + 1e-150 * scale,
#      the last mpmath's own allowance at 200 digits, and the dataflow
#      evaluated exactly must equal the reference within that allowance:
#      that is the formula, held apart from its rounding, at every point,
#      for body k's reaction and for the potential alike. At binary256 the
#      bound must itself be 1e-68 of the term or tighter at every point.
#      (The first version added |exact - reference| to the bound and held
#      only the points to it, so a reaction left out passed its own line:
#      the dropped-reaction control found it.) Beside it, the physics: J4's z
#      part, the potential's gradient with J4 and without, is many bounds
#      above the rounding at every point off the plane, so a dropped or
#      wrong J4 z term cannot hide.
#   3. PHYSICS, on data/problems/zonal_saturn.txt: Saturn with Cassini's J2
#      and J4 at R = 60330 km, a Mimas-like moon inclined 3 degrees and a
#      massless ring particle at the Encke gap inclined 2 degrees, both
#      slightly eccentric, `active 2`. The file is derived here
#      (--write-problem, from its sources' values and the orbits' elements in
#      mpmath) and held to that byte for byte.
#      (a) The record's energy, its zonal potential included, moves no more
#          than ten times what the same orbits move without the zonal line, at
#          every format: the integrator's own level, measured beside it.
#      (c) So does the total momentum, which only the reaction on Saturn can
#          keep.
#      (b) The particle's periapse and node longitudes, from the cubics
#          through r and dr/dt and through z and dz/dt at every step, give the
#          apsidal and nodal rates by least squares. First-order theory, from
#          the potential alone (n^2 = r^-1 dPhi/dr and nu^2 = d2Phi/dz2 in the
#          plane, kappa^2 = r^-3 d(r^4 n^2)/dr, each made linear in J2 and J4
#          by mpmath's derivative in them), at the radius whose mean motion is
#          the one measured, gives n - kappa and n - nu. They must agree within
#          the second order's size: the J terms squared, computed (the exact
#          epicyclic rate less the linear one), plus 2 (e^2 + sin^2 i) of the
#          rate, the size of the classical J2 rates' e^2 and i^2 terms
#          (BELIEVED, from (1 - e^2)^-2 and the inclination factors of the
#          secular J2 rates), e and i measured from the run. And the same
#          theory without J4 must miss: the run resolves J4.
#   4. FORMAT DIFFERENCE. binary128 at epsilon 1e-9 writes its steps
#      (--dt-out); binary128 and binary256 take them (--dt-file, --epsilon 0).
#      Their positions must differ, and by no more than u128 * steps^2 of the
#      distance from Saturn: binary128's round-off, eleven decades under
#      binary64's. Beside it, how far the term itself moved them.
#      ITS SCOPE (the lead's ruling, 2026-11-12): it can fail on a defect that
#      differs between the two formats, a constant or an operation narrowed
#      at one format only, and it cannot fail on a narrowing both formats
#      share. A term rounded to binary64 in both runs is the same binary64
#      value in both, and cancels in their difference (measured: 3.4e-33 on
#      120 steps, against the correct build's 3.4e-33). Parts 1 and 2 hold that
#      one, against references without it: the golden model and mpmath. Its
#      negative control is the term rounded to binary64 at binary128 only.
#   Then: the term reaches both engines and both forms (the program engine's
#      record is the loop engine's in the FMA form, and every engine and form
#      differs from its run without the zonal line); an ensemble of two
#      systems gives each its own run alone (--member); and a test particle's
#      2^-40 mass changes nothing under `active 2`, the energy included.
#   5. THE LINE. Accepted, and named by the record's `# zonal` line with its
#      values and its counts; a corrector schedule recorded under it carries
#      `zonal=` and replays, and is refused under another J2, J4 or R, or
#      without the line, with the keys differing in zonal= and nothing else
#      (the run key's zonal= field is held here; tools/check_corrector.py's
#      walk sees only the fields a run prints). And refused by name, before
#      any record: a zonal body that is a test particle, one past N, a short
#      line.
#
# The selftest's order of operations is the C's (src/ias15_zonal.inc says
# which and why). Nothing here restates it: the dataflow is read from the
# trace, the reference from the potential.

import argparse
import os
import re
import subprocess
import sys

sys.dont_write_bytecode = True        # nothing lands in the pinned clone
here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
sys.path.insert(0, os.path.join(root, "third_party", "cft-fp256", "python"))

import mpmath                          # noqa: E402
from mpmath import mp, mpf             # noqa: E402
from cft_golden import softfloat as sf  # noqa: E402
from cft_golden.formats import FORMATS  # noqa: E402

mp.dps = 200
FMTS = ("fp64", "fp128", "fp256")
CERT = sf.FLAG_INVALID | sf.FLAG_DIVZERO | sf.FLAG_OVERFLOW | sf.FLAG_UNDERFLOW
GOLDEN = {"add": sf.add, "sub": sf.sub, "mul": sf.mul, "div": sf.div, "sqrt": sf.sqrt}


def value(fmt, bits):
    """A finite element's bits as its exact value."""
    u = sf.unpack(fmt, bits)
    if u.kind == sf.ZERO:
        return mpf(0)
    if u.kind not in (sf.NORM, sf.SUB):
        raise ValueError("not finite: %x" % bits)
    v = mpf(u.m) * mpf(2) ** u.e
    return -v if u.sign else v


class Node:
    """One value of the dataflow: its bits, its exact value, its absolute
    evaluation and its rounding count."""
    __slots__ = ("bits", "exact", "absv", "k")

    def __init__(self, bits, exact, absv, k):
        self.bits, self.exact, self.absv, self.k = bits, exact, absv, k


def read_trace(fmtname, text):
    """The selftest's output, checked op by op: (failures, the problem, the
    results as nodes, counts)."""
    fmt = FORMATS[fmtname]
    bad, cur, ops, els = [], {}, 0, 0
    prob = {"bodies": {}}
    acc, pot = {}, {}
    flags_all = 0
    lines = text.splitlines()
    i = 0
    while i < len(lines):
        t = lines[i].split()
        i += 1
        if not t:
            continue
        if t[0] == "zonal" and t[1] == "selftest" and "=" in t[2]:
            kv = dict(x.split("=") for x in t[2:])
            prob.update(E=int(kv["E"]), N=int(kv["N"]), k=int(kv["k"]), active=int(kv["active"]), tptype=int(kv["tptype"]))
        elif t[0] == "zst" and t[1] in ("G", "R", "J2", "J4"):
            prob[t[1]] = value(fmt, int(t[2], 16))
        elif t[0] == "zst" and t[1] == "body":
            s, j = int(t[2]), int(t[3])
            prob["bodies"][(s, j)] = tuple(value(fmt, int(b, 16)) for b in t[4:8])
        elif t[0] == "zin":
            v = value(fmt, int(t[3], 16))
            cur[(t[1], int(t[2]))] = Node(int(t[3], 16), v, abs(v), 0)
        elif t[0] == "zcp":
            src = (t[3], int(t[4]))
            if src not in cur:
                bad.append("a move from %s[%s], which no input declares and no traced operation wrote" % src)
                continue
            cur[(t[1], int(t[2]))] = cur[src]
        elif t[0] == "zop":
            seq, op, dn, an, bn, n = t[1], t[2], t[3], t[4], t[5], int(t[6])
            ops += 1
            for e in range(n):
                z = lines[i].split()
                i += 1
                if z[0] != "zel" or z[1] != seq or int(z[2]) != e:
                    bad.append("operation %s: element line %d is not its element %d" % (seq, i, e))
                    break
                unary = bn == "-"
                ab = int(z[3], 16)
                bb = None if unary else int(z[4], 16)
                rb, fl = int(z[-2], 16), int(z[-1], 16)
                els += 1
                flags_all |= fl
                srcs = [(an, e)] + ([] if unary else [(bn, e)])
                given = [ab] + ([] if unary else [bb])
                nodes = []
                for src, g in zip(srcs, given):
                    nd = cur.get(src)
                    if nd is None:
                        bad.append("operation %s (%s) element %d: its operand %s[%d] is no declared input and no "
                                   "traced result" % (seq, op, e, src[0], src[1]))
                    elif nd.bits != g:
                        bad.append("operation %s (%s) element %d: its operand %s[%d] is %x, and the trace last left "
                                   "%x there" % (seq, op, e, src[0], src[1], g, nd.bits))
                    nodes.append(nd)
                gb, gf = GOLDEN[op](fmt, ab) if unary else GOLDEN[op](fmt, ab, bb)
                if (gb, gf) != (rb, fl):
                    bad.append("operation %s (%s %s = %s %s %s) element %d: the library gave %x flags %x, the golden "
                               "model %x flags %x" % (seq, op, dn, an, op, bn, e, rb, fl, gb, gf))
                if len(nodes) != len(srcs) or any(nd is None for nd in nodes):
                    continue
                cur[(dn, e)] = combine(op, nodes, rb, fmt)
        elif t[0] == "zst" and t[1] == "acc":
            s, j = int(t[2]), int(t[3])
            acc[(s, j)] = [int(b, 16) for b in t[4:7]]
        elif t[0] == "zst" and t[1] == "pot":
            pot[int(t[2])] = int(t[3], 16)
    # the results printed must be the trace's
    res_acc, res_pot = {}, {}
    for (s, j), bits in acc.items():
        nds = []
        for c in range(3):
            nd = cur.get(("a", 3 * (s * prob["N"] + j) + c))
            if nd is None or nd.bits != bits[c]:
                bad.append("the printed acceleration of system %d body %d, component %d, is not the trace's" % (s, j, c))
            nds.append(nd)
        res_acc[(s, j)] = nds
    for s, bits in pot.items():
        nd = cur.get(("epot", s))
        if nd is None or nd.bits != bits:
            bad.append("the printed potential of system %d is not the trace's" % s)
        res_pot[s] = nd
    if flags_all & CERT:
        bad.append("the term raised flags 0x%x: only inexact is a term's to raise" % flags_all)
    return bad, prob, res_acc, res_pot, ops, els, flags_all


def combine(op, nodes, rbits, fmt):
    """The node an operation makes: its exact value and absolute evaluation
    from its operands', and its rounding count (Higham's rules)."""
    a = nodes[0]
    if op == "sqrt":
        if a.absv != a.exact:
            raise SystemExit("check_zonal: a square root of a value with cancellation in it; the bound below does not apply")
        return Node(rbits, mpmath.sqrt(a.exact), mpmath.sqrt(a.absv), a.k + 1)
    b = nodes[1]
    if op in ("add", "sub"):
        ex = a.exact + b.exact if op == "add" else a.exact - b.exact
        if a.k == 0 and b.k == 0:          # two exact operands: (a + b)(1 + d), one monomial
            return Node(rbits, ex, abs(ex), 1)
        return Node(rbits, ex, a.absv + b.absv, max(a.k, b.k) + 1)
    if op == "mul":
        return Node(rbits, a.exact * b.exact, a.absv * b.absv, a.k + b.k + 1)
    # div: the divisor's own exact value, which must carry no cancellation
    if abs(b.absv - abs(b.exact)) > abs(b.exact) * mpf(10) ** -150:
        raise SystemExit("check_zonal: a divisor with cancellation in it; the bound below does not apply")
    return Node(rbits, a.exact / b.exact, a.absv / abs(b.exact), a.k + b.k + 1)


def gamma(k, p):
    u = mpf(2) ** -p
    return k * u / (1 - k * u)


def reference(prob, s, j, j4=True):
    """Minus the gradient of the zonal potential at body j of system s, from
    the potential alone: mpmath.diff of G M [J2 R^2 r^-3 P2(u) + J4 R^4 r^-5 P4(u)]."""
    G, R, J2, J4 = prob["G"], prob["R"], prob["J2"], prob["J4"] if j4 else mpf(0)
    M = prob["bodies"][(s, prob["k"])][0]
    xk = prob["bodies"][(s, prob["k"])][1:]
    xj = prob["bodies"][(s, j)][1:]

    def phi(x, y, z):
        r = mpmath.sqrt(x * x + y * y + z * z)
        u = z / r
        return G * M * (J2 * R ** 2 * r ** -3 * mpmath.legendre(2, u) + J4 * R ** 4 * r ** -5 * mpmath.legendre(4, u))

    d = [xj[c] - xk[c] for c in range(3)]
    grad = [-mpmath.diff(phi, d, o) for o in ((1, 0, 0), (0, 1, 0), (0, 0, 1))]
    return grad, phi(*d), d


def held(fmtname, p, nd, ref, scale):
    """One computed component against its reference: (|error|, its bound,
    whether the dataflow evaluated exactly is the reference). The bound is
    gamma_k times the absolute evaluation, plus 1e-150 of the scale for
    mpmath's own error at 200 digits; the exact dataflow must be within
    that allowance of the reference, or the formula itself is wrong."""
    allow = mpf(10) ** -150 * (nd.absv + scale)
    err = abs(value(FORMATS[fmtname], nd.bits) - ref)
    return err, gamma(nd.k, p) * nd.absv + allow, abs(nd.exact - ref) <= allow


def gate_formula(fmtname, prob, res_acc, res_pot, verdict, quiet=False):
    """Gate 2 at one format: a line a point, and the system's sums."""
    p = FORMATS[fmtname].prec
    N, k, na, tp = prob["N"], prob["k"], prob["active"], prob["tptype"]
    ni = N if tp else na
    worst_rel_bound, least_j4 = mpf(0), None
    ok_all = True
    lines = []
    for s in range(prob["E"]):
        M = prob["bodies"][(s, k)][0]
        react = [mpf(0)] * 3
        upot = mpf(0)
        for j in range(N):
            if j == k:
                continue
            ref, ph, d = reference(prob, s, j)
            ref0, _, _ = reference(prob, s, j, j4=False)
            m = prob["bodies"][(s, j)][0]
            if j < ni:
                for c in range(3):
                    react[c] -= m / M * ref[c]
                upot += m * ph
            nds = res_acc[(s, j)]
            norm = mpmath.sqrt(sum(v * v for v in ref))
            plane = d[2] == 0
            errs, bounds, exacts = [], [], []
            for c in range(3):
                e_, b_, x_ = held(fmtname, p, nds[c], ref[c], norm)
                errs.append(e_), bounds.append(b_), exacts.append(x_)
            ok = all(e <= b for e, b in zip(errs, bounds)) and all(exacts)
            relb = max(bounds) / norm
            worst_rel_bound = max(worst_rel_bound, relb)
            j4z = abs(ref[2] - ref0[2])
            if not plane:
                r = j4z / bounds[2]
                least_j4 = r if least_j4 is None or r < least_j4 else least_j4
            ok_all &= ok
            lines.append((ok, "system %d body %d, %s, r = %s R: |error| %s of the term, the bound %s (k up to %d)%s%s"
                          % (s, j, "on the plane" if plane else "off the plane, u = %s" % mpmath.nstr(d[2] / mpmath.sqrt(sum(v * v for v in d)), 3),
                             mpmath.nstr(mpmath.sqrt(sum(v * v for v in d)) / prob["R"], 3),
                             mpmath.nstr(max(errs) / norm, 3), mpmath.nstr(relb, 3), max(nd.k for nd in nds),
                             "" if plane else "; J4's z part %s bounds above it" % mpmath.nstr(j4z / bounds[2], 3),
                             "" if all(exacts) else "; THE FORMULA ITSELF differs from the gradient")))
        # body k: the interacting bodies' reactions; the potential
        nds = res_acc[(s, k)]
        norm = mpmath.sqrt(sum(v * v for v in react))
        hk = [held(fmtname, p, nds[c], react[c], norm) for c in range(3)]
        ok = all(e <= b and x for e, b, x in hk) and norm > 0
        ok_all &= ok
        lines.append((ok, "system %d body %d (the oblate one): the reaction of the %d bodies that act on it, |error| %s of "
                      "it, the bound %s%s" % (s, k, ni - 1, mpmath.nstr(max(h[0] for h in hk) / norm, 3),
                                              mpmath.nstr(max(h[1] for h in hk) / norm, 3),
                                              "" if all(h[2] for h in hk) else "; THE FORMULA ITSELF differs from the gradient")))
        perr, pb, px = held(fmtname, p, res_pot[s], upot, abs(upot))
        ok = perr <= pb and px
        ok_all &= ok
        lines.append((ok, "system %d potential, the J terms of the %d interacting bodies: |error| %s of it, the bound %s%s"
                      % (s, ni - 1, mpmath.nstr(perr / abs(upot), 3), mpmath.nstr(pb / abs(upot), 3),
                         "" if px else "; THE FORMULA ITSELF differs from the potential")))
    for ok, text in lines:
        if not quiet or not ok:
            verdict(ok, "%s: %s" % (fmtname, text))
    return ok_all, worst_rel_bound, least_j4


# ---------------------------------------------------------------------
# The physics problem: derived here, and the committed file held to it
# ---------------------------------------------------------------------
PROBLEM = os.path.join(root, "data", "problems", "zonal_saturn.txt")
# Each value as its source gives it; the file holds the nearest binary64.
# Units km and s, G = 1, a mass is a GM in km^3 s^-2.
SATURN = {"GM": "37.931e6", "R": "60330", "J2": "16290.573e-6", "J4": "-935.314e-6"}
MOON = {"GM": "2.50349", "a": "186000", "e": "0.02", "i": "3", "node": "40", "peri": "70", "M": "10"}
RING = {"a": "133410", "e": "0.02", "i": "2", "node": "200", "peri": "30", "M": "250"}
HEAD = """\
# Saturn with J2 and J4, a Mimas-like moon and one ring particle: the
# physics problem of the zonal term's gate (round 4's P1). Written by
# `python tools/check_zonal.py --write-problem`, which derives every value
# below and which the gate holds this file to, byte for byte.
#
# Units km and s, G = 1, each mass a GM in km^3 s^-2; each number is the
# nearest binary64 to the source's decimal (or to the state its elements
# give, computed in mpmath at 50 digits), written exactly as hex.
#
#   Saturn's GM, 37.931e6: NASA's Saturn fact sheet
#     (nssdc.gsfc.nasa.gov/planetary/factsheet/saturnfact.html).
#   J2 = 16290.573e-6 and J4 = -935.314e-6: Iess et al. 2019, Science 364,
#     eaat2965 (Cassini's Grand Finale), as tabulated by Fortney et al.,
#     arXiv:2304.09215, Table 4.1. Both read there by cft-fp256's study
#     GPU-E, section 4, which names no reference radius for them.
#   R = 60330 km, the radius those J2 and J4 are referred to: BELIEVED (the
#     study's belief, and this file's), the reference radius of Cassini's
#     gravity solutions. It is not Saturn's 1-bar equatorial radius,
#     60,268 km: a J_n means something only with its own radius.
#   The moon: Mimas's GM, 2.50349, and its a = 186,000 km and e = 0.02
#     (JPL SSD, via the same study). Its inclination of 3 degrees, node,
#     periapse and phase are synthetic.
#   The ring particle: massless, a = 133,410 km (the Encke gap, NASA's ring
#     fact sheet, via the study), e = 0.02 and i = 2 degrees, synthetic.
#   Both orbits are Keplerian about Saturn alone at t = 0; in the oblate
#   field they precess. `active 2`: the particle is a test particle.
"""


def to_hex(v):
    """v, an mpf, rounded once to the nearest binary64, as exact hex."""
    from fractions import Fraction
    s, m, e, _ = v._mpf_
    f = Fraction(int(m)) * Fraction(2) ** int(e) if m else Fraction(0)
    return float.hex(float(-f if s else f))


def kepler_state(mu, el):
    """Position and velocity from Keplerian elements, in mpmath."""
    with mp.workdps(50):
        a, e = mpf(el["a"]), mpf(el["e"])
        i, node, peri, M = (mpmath.radians(mpf(el[k])) for k in ("i", "node", "peri", "M"))
        E = mpmath.findroot(lambda E: E - e * mpmath.sin(E) - M, M)
        b = a * mpmath.sqrt(1 - e * e)
        n = mpmath.sqrt(mu / a ** 3)
        Ed = n / (1 - e * mpmath.cos(E))
        p = (a * (mpmath.cos(E) - e), b * mpmath.sin(E))
        q = (-a * mpmath.sin(E) * Ed, b * mpmath.cos(E) * Ed)
        co, so, ci, si, cw, sw = (mpmath.cos(node), mpmath.sin(node), mpmath.cos(i), mpmath.sin(i),
                                  mpmath.cos(peri), mpmath.sin(peri))
        rot = ((co * cw - so * sw * ci, -co * sw - so * cw * ci),
               (so * cw + co * sw * ci, -so * sw + co * cw * ci),
               (sw * si, cw * si))
        return ([rot[c][0] * p[0] + rot[c][1] * p[1] for c in range(3)],
                [rot[c][0] * q[0] + rot[c][1] * q[1] for c in range(3)])


def problem_text(zonal=True, extra_active=True):
    """The problem, as the file holds it; zonal=False drops the zonal line,
    which is the same orbits without the term."""
    with mp.workdps(50):
        gm = mpf(SATURN["GM"])
        mx, mv = kepler_state(gm + mpf(MOON["GM"]), MOON)
        rx, rv = kepler_state(gm, RING)
        lines = [HEAD.rstrip("\n"), "name zonal_saturn", "G %s" % to_hex(mpf(1))]
        if extra_active:
            lines.append("active 2")
        if zonal:
            lines.append("zonal 0 %s %s %s" % tuple(to_hex(mpf(SATURN[k])) for k in ("R", "J2", "J4")))
        z = to_hex(mpf(0))
        lines.append("body saturn %s %s" % (to_hex(gm), " ".join([z] * 6)))
        lines.append("body mimas %s %s" % (to_hex(mpf(MOON["GM"])), " ".join(to_hex(v) for v in mx + mv)))
        lines.append("body ring %s %s" % (z, " ".join(to_hex(v) for v in rx + rv)))
    return "\n".join(lines) + "\n"


# ---------------------------------------------------------------------
# Records, and the physics read from them
# ---------------------------------------------------------------------
def hexval(s):
    """A record's hex float, exactly."""
    neg = s.startswith("-")
    s = s.lstrip("+-")
    if not s.startswith("0x"):
        raise ValueError("not a hex float: %s" % s)
    m, e = s[2:].split("p")
    ip, _, fp = m.partition(".")
    v = mpf(int(ip + fp, 16)) * mpf(2) ** (int(e) - 4 * len(fp))
    return -v if neg else v


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=root)
    return p.returncode, p.stdout, p.stderr


def record(out):
    """A one-system record: its masses, its samples [(t, E, [[x y z vx vy vz]
    a body])], and its trailer as a dict."""
    masses, samples, tail = [], [], {}
    for line in out.splitlines():
        m = re.match(r"^# body (\d+) \S+ m= (\S+)$", line)
        if m:
            masses.append(hexval(m.group(2)))
        elif line.startswith("# steps_done="):
            tail = dict(kv.split("=", 1) for kv in line[2:].split())
        elif line.startswith("sample "):
            t = line.split()
            bar = t.index("|")
            v = t[3:bar]
            body = [hexval(x) for x in v[4:]]
            samples.append((hexval(t[bar + 1]) + hexval(t[bar + 2]), hexval(v[3]),
                            [body[6 * b:6 * b + 6] for b in range(len(body) // 6)]))
    return masses, samples, tail


def energy_drift(samples):
    E0 = samples[0][1]
    return max(abs(s[1] - E0) for s in samples) / abs(E0)


def momentum_drift(samples, masses):
    """The total momentum's largest excursion from its start, over the sum of
    the bodies' |m v| at the start (masses are GM: the same ratio)."""
    def P(st):
        return [sum(m * b[3 + c] for m, b in zip(masses, st)) for c in range(3)]
    p0 = P(samples[0][2])
    scale = sum(abs(m) * mpmath.sqrt(sum(b[3 + c] ** 2 for c in range(3))) for m, b in zip(masses, samples[0][2]))
    return max(mpmath.sqrt(sum((a - b) ** 2 for a, b in zip(P(s[2]), p0))) for s in samples) / scale


def hermite(t0, t1, f0, d0, f1, d1, t):
    """The cubic through f and its derivative d at both ends."""
    h = t1 - t0
    s = (t - t0) / h
    return ((2 * s ** 3 - 3 * s ** 2 + 1) * f0 + (s ** 3 - 2 * s ** 2 + s) * h * d0 +
            (-2 * s ** 3 + 3 * s ** 2) * f1 + (s ** 3 - s ** 2) * h * d1)


def root_in(fn, a, b):
    """fn's sign change in [a, b], by bisection to the working precision."""
    fa = fn(a)
    for _ in range(200):
        m = (a + b) / 2
        fm = fn(m)
        if (fm >= 0) == (fa >= 0):
            a, fa = m, fm
        else:
            b = m
    return (a + b) / 2


def apsides_and_nodes(samples, ib, ik=0):
    """Body ib about body ik: each periapse (a minimum of r, from the cubic
    through r and dr/dt) and each ascending node (z upward through 0, from
    the cubic through z and dz/dt), with the body's longitude there. Returns
    the rates of the periapse's and the node's longitudes by least squares
    over the events, the mean motion between the first and last periapse,
    and the orbit's e and sin i."""
    rel = []
    for t, _, st in samples:
        a, b = st[ib], st[ik]
        rel.append((t, [a[c] - b[c] for c in range(6)]))
    th_prev, turns, unwrapped = None, 0, []
    rmin, rmax, zmax = None, None, mpf(0)
    for t, q in rel:
        th = mpmath.atan2(q[1], q[0])
        if th_prev is not None and th < th_prev - mpmath.pi:
            turns += 1
        th_prev = th
        unwrapped.append(th + 2 * mpmath.pi * turns)
        r = mpmath.sqrt(q[0] ** 2 + q[1] ** 2 + q[2] ** 2)
        rmin = r if rmin is None or r < rmin else rmin
        rmax = r if rmax is None or r > rmax else rmax
        zmax = max(zmax, abs(q[2]))

    def lon_at(i, t):
        (t0, a), (t1, b) = rel[i], rel[i + 1]
        x = hermite(t0, t1, a[0], a[3], b[0], b[3], t)
        y = hermite(t0, t1, a[1], a[4], b[1], b[4], t)
        th = mpmath.atan2(y, x)
        base = unwrapped[i]
        return th + 2 * mpmath.pi * mpmath.nint((base - th) / (2 * mpmath.pi))

    peri, node = [], []
    for i in range(len(rel) - 1):
        (t0, a), (t1, b) = rel[i], rel[i + 1]
        r0, r1 = mpmath.sqrt(sum(a[c] ** 2 for c in range(3))), mpmath.sqrt(sum(b[c] ** 2 for c in range(3)))
        v0, v1 = sum(a[c] * a[3 + c] for c in range(3)) / r0, sum(b[c] * b[3 + c] for c in range(3)) / r1
        if v0 < 0 <= v1:
            h = t1 - t0

            def drdt(t, t0=t0, r0=r0, v0=v0, r1=r1, v1=v1, h=h):
                s = (t - t0) / h
                return ((6 * s ** 2 - 6 * s) * r0 + (3 * s ** 2 - 4 * s + 1) * h * v0 +
                        (-6 * s ** 2 + 6 * s) * r1 + (3 * s ** 2 - 2 * s) * h * v1) / h
            tp = root_in(drdt, t0, t1)
            peri.append((tp, lon_at(i, tp)))
        if a[2] < 0 <= b[2]:
            tn = root_in(lambda t, t0=t0, t1=t1, a=a, b=b: hermite(t0, t1, a[2], a[5], b[2], b[5], t), t0, t1)
            node.append((tn, lon_at(i, tn)))

    def rate(ev):
        # the longitude less a turn per event, against time, by least squares
        ys = [lon - 2 * mpmath.pi * k for k, (_, lon) in enumerate(ev)]
        ts = [t for t, _ in ev]
        n = len(ev)
        tm, ym = sum(ts) / n, sum(ys) / n
        return sum((t - tm) * (y - ym) for t, y in zip(ts, ys)) / sum((t - tm) ** 2 for t in ts)
    if len(peri) < 3 or len(node) < 3:
        return None
    nmean = (peri[-1][1] - peri[0][1]) / (peri[-1][0] - peri[0][0])
    rg = (rmin + rmax) / 2
    return {"apse": rate(peri), "node": rate(node), "n": nmean, "e": (rmax - rmin) / (rmax + rmin),
            "sini": zmax / rg, "periapses": len(peri), "nodes": len(node), "span": peri[-1][0] - peri[0][0]}


def theory(GM, R, J2, J4, n_meas):
    """First-order theory in J2 and J4, from the potential alone: the
    frequencies of a near-circular, near-equatorial orbit at the radius
    whose mean motion is the one measured - n^2 = (1/r) dPhi/dr and
    nu^2 = d2Phi/dz2 in the plane, kappa^2 = r^-3 d(r^4 n^2)/dr - each made
    linear in J2 and J4 by mpmath's derivative in them at zero."""
    with mp.workdps(60):
        def phi(x, y, z, j2, j4):
            r = mpmath.sqrt(x * x + y * y + z * z)
            u = z / r
            return -GM / r + GM * (j2 * R ** 2 * r ** -3 * mpmath.legendre(2, u) + j4 * R ** 4 * r ** -5 * mpmath.legendre(4, u))

        def n2(r, j2, j4):
            return mpmath.diff(lambda q: phi(q, 0, 0, j2, j4), r) / r

        def freqs(r, j2, j4):
            n = mpmath.sqrt(n2(r, j2, j4))
            k = mpmath.sqrt(mpmath.diff(lambda q: q ** 4 * n2(q, j2, j4), r) / r ** 3)
            nu = mpmath.sqrt(mpmath.diff(lambda z: phi(r, 0, z, j2, j4), 0, 2))
            return n, k, nu

        def linear(f, r):
            # f(r, J2, J4) to first order in J2 and J4
            return (f(r, 0, 0) + J2 * mpmath.diff(lambda j: f(r, j, 0), 0) +
                    J4 * mpmath.diff(lambda j: f(r, 0, j), 0))
        n1 = lambda r: linear(lambda r, a, b: freqs(r, a, b)[0], r)
        rg = mpmath.findroot(lambda r: n1(r) - n_meas, (GM / n_meas ** 2) ** (mpf(1) / 3))
        apse = linear(lambda r, a, b: freqs(r, a, b)[0] - freqs(r, a, b)[1], rg)
        node = linear(lambda r, a, b: freqs(r, a, b)[0] - freqs(r, a, b)[2], rg)
        n, k, nu = freqs(rg, J2, J4)
        return {"rg": rg, "apse": apse, "node": node, "apse_exact": n - k, "node_exact": n - nu,
                "eps2": J2 * (R / rg) ** 2}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true", help="fewer orbits at binary128 and binary256")
    ap.add_argument("--write-problem", action="store_true", help="write data/problems/zonal_saturn.txt and stop")
    args = ap.parse_args()
    if args.write_problem:
        with open(PROBLEM, "w", newline="\n") as f:
            f.write(problem_text())
        print("wrote %s" % os.path.relpath(PROBLEM, root))
        return 0
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, "ias15_cft" + exe)
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= bool(ok)
        print("  %s %s" % ("ok  " if ok else "FAIL", text))

    # 1 and 2. The arithmetic and the formula, at every format.
    for fmtname in FMTS:
        p = subprocess.run([cft, "--zonal-selftest", "--format", fmtname], stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE, text=True, cwd=root)
        print("check_zonal at %s: the term at the selftest's points, every operation against the golden model"
              % fmtname)
        if p.returncode:
            verdict(False, "--zonal-selftest rc %d: %s" % (p.returncode, p.stderr.strip()[-300:]))
            continue
        bad, prob, res_acc, res_pot, ops, els, fl = read_trace(fmtname, p.stdout)
        # the physics beside it: z's term is +0 on the plane and nonzero off it
        live, dead, wrongs = 0, 0, []
        for (s, j), nds in sorted(res_acc.items()):
            d = [prob["bodies"][(s, j)][1 + c] - prob["bodies"][(s, prob["k"])][1 + c] for c in range(3)]
            zb = nds[2].bits if nds[2] is not None else None
            if j == prob["k"]:
                continue
            if d[2] == 0:
                dead += 1
                if zb is None or value(FORMATS[fmtname], zb) != 0:
                    wrongs.append("system %d body %d is on the plane and its z term is not 0" % (s, j))
            else:
                live += 1
                if any(nd is None or value(FORMATS[fmtname], nd.bits) == 0 for nd in nds):
                    wrongs.append("system %d body %d is off the plane and a component of its term is 0" % (s, j))
        verdict(not bad and not wrongs and ops > 0,
                "%d operations, %d elements, each against the golden model's bits and flags (flags seen 0x%x), each "
                "operand against a declared input or the trace's own last result, each printed result against the "
                "trace: %d mismatches; %d points off the plane (every component wanted nonzero) and %d on it (z's "
                "wanted exactly 0): %d wrong%s"
                % (ops, els, fl, len(bad), live, dead, len(wrongs), "".join("\n         " + b for b in (bad + wrongs)[:12])))
        print("check_zonal at %s: the term against the potential's gradient, mpmath at %d digits" % (fmtname, mp.dps))
        ok, worst, least = gate_formula(fmtname, prob, res_acc, res_pot, verdict)
        if fmtname == "fp256":
            verdict(worst <= mpf("1e-68"), "fp256: the bound is %s of the term at its widest, against the 1e-68 asked"
                    % mpmath.nstr(worst, 3))
        verdict(least is not None and least > 1000, "%s: J4's z part over the bound, at its least off the plane: %s "
                "(wanted above 1000, so that it cannot hide in the rounding)"
                % (fmtname, mpmath.nstr(least, 3) if least else "no point"))

    work = os.path.join(args.build, "check_zonal")
    os.makedirs(work, exist_ok=True)

    def write(name, text):
        path = os.path.join(work, name)
        with open(path, "w", newline="\n") as f:
            f.write(text)
        return path

    # The problem: the committed file is the derived one, byte for byte.
    print("check_zonal: the physics problem, Saturn with J2 and J4, a moon and a ring particle")
    try:
        with open(PROBLEM, encoding="utf-8", newline="") as f:
            committed = f.read()
    except OSError:
        committed = None
    verdict(committed == problem_text(), "%s against what --write-problem derives from its sources' values, byte for "
            "byte: %s" % (os.path.relpath(PROBLEM, root).replace(os.sep, "/"),
                          "the same" if committed == problem_text() else "DIFFERENT (regenerate it, and review the diff)"))
    nozonal = write("zonal_saturn_nozonal.txt", problem_text(zonal=False))
    head = problem_text()
    R_hex, J2_hex, J4_hex = re.search(r"^zonal 0 (\S+) (\S+) (\S+)$", head, re.M).groups()
    s64 = 600 if args.quick else 1100          # binary64: about 16 or 30 of the particle's orbits
    sw = 40 if args.quick else 120             # binary128 and binary256

    def integrate(fmt, prob, steps, extra=()):
        a = [cft, "--format", fmt, "--problem", prob, "--steps", str(steps), "--sample", "1", "--quiet"] + list(extra)
        rc, out, err = run(a)
        if rc:
            verdict(False, "%s: rc %d: %s" % (" ".join(os.path.basename(x) for x in a[1:6]), rc, err.strip()[-300:]))
            return None
        return record(out) + (out,)

    adapt = ["--dt", "100", "--epsilon", "1e-9"]
    dtf = os.path.join(work, "zonal_saturn_fp128.dt")
    runs = {
        ("fp64", "zonal"): integrate("fp64", PROBLEM, s64, adapt),
        ("fp64", "none"): integrate("fp64", nozonal, s64, adapt),
        ("fp128", "zonal"): integrate("fp128", PROBLEM, sw, adapt + ["--dt-out", dtf]),
        ("fp128", "none"): integrate("fp128", nozonal, sw, adapt),
    }
    replay = ["--dt-file", dtf, "--epsilon", "0"]
    runs[("fp128", "replay")] = integrate("fp128", PROBLEM, sw, replay)
    runs[("fp256", "zonal")] = integrate("fp256", PROBLEM, sw, replay)
    runs[("fp256", "none")] = integrate("fp256", nozonal, sw, replay)
    if any(v is None for v in runs.values()):
        print("RESULT: FAIL")
        return 1

    # 3(a) and 3(c): energy and momentum, against the same orbits without the term
    print("check_zonal: 3(a) the energy with the zonal potential, and 3(c) the momentum with the reaction, at every format")
    for fmt in FMTS:
        z, n = runs[(fmt, "zonal")], runs[(fmt, "none")]
        u = mpf(2) ** -FORMATS[fmt].prec
        ez, en = energy_drift(z[1]), energy_drift(n[1])
        pz, pn = momentum_drift(z[1], z[0]), momentum_drift(n[1], n[0])
        how = "adaptive" if fmt != "fp256" else "on binary128's steps"
        verdict(ez <= 10 * max(en, u), "%s, %d steps %s: the energy, Phi's J terms included, moves %s of itself; the same "
                "orbits without the term, %s (the integrator's own level; the bound is 10 times it)"
                % (fmt, len(z[1]) - 1, how, mpmath.nstr(ez, 3), mpmath.nstr(en, 3)))
        verdict(pz <= 10 * max(pn, u), "%s: the total momentum moves %s of the bodies' |m v|, and %s without the term "
                "(the bound is 10 times it); body 0 takes the moon's reaction" % (fmt, mpmath.nstr(pz, 3), mpmath.nstr(pn, 3)))

    # 3(b): the apsides and the node, against first-order theory
    print("check_zonal: 3(b) the ring particle's apsides and node, against first-order theory in J2 and J4")
    masses, samples = runs[("fp64", "zonal")][0], runs[("fp64", "zonal")][1]
    meas = apsides_and_nodes(samples, 2)
    if meas is None:
        verdict(False, "fewer than three periapses or nodes in %d samples" % len(samples))
    else:
        GM, R, J2, J4 = masses[0], hexval(R_hex), hexval(J2_hex), hexval(J4_hex)
        th = theory(GM, R, J2, J4, meas["n"])
        th0 = theory(GM, R, J2, mpf(0), meas["n"])
        second = 2 * (meas["e"] ** 2 + meas["sini"] ** 2)
        orbits = meas["span"] * meas["n"] / (2 * mpmath.pi)
        for key, name, sense in (("apse", "the apsidal rate (wanted positive: the apsides advance)", 1),
                                 ("node", "the nodal rate (wanted negative: the node regresses)", -1)):
            first, exact, got, nj4 = th[key], th[key + "_exact"], meas[key], th0[key]
            tol = abs(exact - first) + abs(first) * second
            ok = abs(got - first) <= tol and got * sense > 0
            verdict(ok, "%s: %s rad/s over %s orbits; first-order theory %s (n - %s); off by %s, against the second "
                    "order's %s (the J terms squared %s, and 2(e^2 + sin^2 i) = %s with e = %s, i = %s deg)"
                    % (name, mpmath.nstr(got, 6), mpmath.nstr(orbits, 3), mpmath.nstr(first, 6),
                       "kappa" if key == "apse" else "nu", mpmath.nstr(abs(got / first - 1), 3), mpmath.nstr(tol / abs(first), 3),
                       mpmath.nstr(abs(exact / first - 1), 3), mpmath.nstr(second, 3), mpmath.nstr(meas["e"], 3),
                       mpmath.nstr(mpmath.degrees(mpmath.asin(meas["sini"])), 3)))
            verdict(abs(got - nj4) > tol, "%s, J4 resolved: the same theory without J4 gives %s, %s off the measured "
                    "rate (wanted past the second order's %s)"
                    % (key == "apse" and "the apsidal rate" or "the nodal rate", mpmath.nstr(nj4, 6),
                       mpmath.nstr(abs(got / nj4 - 1), 3), mpmath.nstr(tol / abs(first), 3)))

    # 4: binary128 against binary256 on the same steps
    print("check_zonal: 4. binary128 against binary256 on binary128's own steps (--dt-out, then --dt-file at epsilon 0)")
    a, b, nz = runs[("fp128", "replay")][1], runs[("fp256", "zonal")][1], runs[("fp256", "none")][1]
    worst, moved = mpf(0), mpf(0)
    for (_, _, sa), (_, _, sb), (_, _, sn) in zip(a, b, nz):
        for i in range(1, len(sb)):
            r = mpmath.sqrt(sum((sb[i][c] - sb[0][c]) ** 2 for c in range(3)))
            worst = max(worst, mpmath.sqrt(sum((sa[i][c] - sb[i][c]) ** 2 for c in range(3))) / r)
            moved = max(moved, mpmath.sqrt(sum((sn[i][c] - sb[i][c]) ** 2 for c in range(3))) / r)
    bound = mpf(2) ** -113 * sw ** 2
    verdict(len(a) == len(b) == sw + 1 and 0 < worst <= bound,
            "%d samples on the same %d steps: the positions differ by %s of the distance from Saturn at most, against "
            "binary128's round-off, u128 * steps^2 = %s (binary64's: %s), and must differ at all; the term itself "
            "moved them %s. Its scope: a precision that differs between the formats, never a narrowing both share "
            "(parts 1 and 2 hold that)" % (len(a), sw, mpmath.nstr(worst, 3), mpmath.nstr(bound, 3),
                                           mpmath.nstr(mpf(2) ** -53 * sw ** 2, 3), mpmath.nstr(moved, 3)))

    # Both engines and both forms see the term
    print("check_zonal: both engines and both forms, and an ensemble against its members")
    eng = (("loop", "rebound"), ("loop", "fma"), ("program", "fma"))
    for fmt in (("fp64",) if args.quick else FMTS):
        recs = {}
        for engine, arith in eng:
            for prob, tag in ((PROBLEM, "zonal"), (nozonal, "none")):
                extra = ["--arith", arith, "--engine", engine] + (["--programs", "programs/out"] if engine == "program" else [])
                rc, out, err = run([cft, "--format", fmt, "--problem", prob, "--steps", "12", "--sample", "4", "--quiet"]
                                   + adapt + extra)
                recs[(engine, arith, tag)] = [l for l in out.splitlines() if l.startswith("sample ")] if rc == 0 else None
        same = recs[("loop", "fma", "zonal")] is not None and recs[("loop", "fma", "zonal")] == recs[("program", "fma", "zonal")]
        seen = all(recs[(e, a_, "zonal")] and recs[(e, a_, "zonal")] != recs[(e, a_, "none")] for e, a_ in eng)
        verdict(same and seen, "%s: the program engine against the loop engine's FMA-form record, %d samples: %s; each "
                "engine and form against its run without the zonal line: %s"
                % (fmt, len(recs[("loop", "fma", "zonal")] or []), "identical" if same else "NOT identical",
                   "every one differs" if seen else "NOT every one differs, or one did not run"))
    lines = problem_text().splitlines()
    heads = [l for l in lines if l and not l.startswith("body ")]
    bodies = [l for l in lines if l.startswith("body ")]
    moon2 = bodies[1].split()
    moon2[2] = to_hex(hexval(moon2[2]) * 3)                 # system two's moon three times as massive
    ens = write("zonal_saturn_x2.txt", "\n".join(heads + ["N 3", "E 2", "system one"] + bodies + ["system two"]
                                                 + bodies[:1] + [" ".join(moon2)] + bodies[2:]) + "\n")
    common = ["--format", "fp64", "--problem", ens, "--steps", "12", "--sample", "4", "--quiet"] + adapt
    rc, out, err = run([cft] + common)
    whole = [l.split(None, 2) for l in out.splitlines() if l.startswith("system ")] if rc == 0 else []
    members = []
    for k in (0, 1):
        rc_k, out_k, _ = run([cft] + common + ["--member", str(k)])
        mine = [w[2] for w in whole if w[1] == str(k)]
        members.append(rc_k == 0 and mine and mine == [l for l in out_k.splitlines() if l.startswith("sample ")])
    differ = [w[2] for w in whole if w[1] == "0"] != [w[2] for w in whole if w[1] == "1"]
    verdict(rc == 0 and all(members) and differ, "an ensemble of two (the second moon three times as massive): each "
            "system against its member run alone: %s; the two systems against each other: %s"
            % ("both identical" if all(members) else "NOT both identical (rc %d)" % rc,
               "they differ" if differ else "THE SAME, so the comparison cannot fail"))

    # A test particle's mass moves no active body here either
    t40 = problem_text().replace("body ring 0x0.0p+0 ", "body ring 0x1.0000000000000p-40 ")
    m40 = write("zonal_saturn_ring_m40.txt", t40)
    rc1, o1, _ = run([cft, "--format", "fp64", "--problem", PROBLEM, "--steps", "12", "--sample", "4", "--quiet"] + adapt)
    rc2, o2, _ = run([cft, "--format", "fp64", "--problem", m40, "--steps", "12", "--sample", "4", "--quiet"] + adapt)
    s1 = [l for l in o1.splitlines() if l.startswith("sample ")]
    same40 = bool(s1) and s1 == [l for l in o2.splitlines() if l.startswith("sample ")]
    verdict(rc1 == 0 and rc2 == 0 and same40 and t40 != problem_text(),
            "under active 2, the ring particle with a mass of 2^-40 against massless, %d samples, the energy column "
            "included: %s (a test particle takes the term and gives no reaction; rc %d, %d)"
            % (len(s1), "identical" if same40 else "THEY DIFFER", rc1, rc2))

    # 5. The zonal line: accepted, its record line, its run key, its refusals
    print("check_zonal: 5. the zonal line, its record line, the run key's zonal= field, and the refusals")
    rc, out, err = run([cft, "--format", "fp64", "--problem", PROBLEM, "--steps", "2", "--sample", "1", "--quiet"] + adapt)
    m = re.search(r"^# zonal body (\d+) of every system, about z: R= (\S+) J2= (\S+) J4= (\S+); it acts on (\d+) bodies a "
                  r"system, and body (\d+) takes the reaction of (\d+)$", out, re.M)
    want = (0, hexval(R_hex), hexval(J2_hex), hexval(J4_hex), 2, 0, 1)
    got = (int(m.group(1)), hexval(m.group(2)), hexval(m.group(3)), hexval(m.group(4)), int(m.group(5)), int(m.group(6)),
           int(m.group(7))) if m else None
    verdict(rc == 0 and got == want, "the zonal line: rc %d, and the record's # zonal line against the file's body 0, "
            "R, J2 and J4, acting on 2 bodies (N - 1) and taking the reaction of 1 (active 2, less body 0): %s"
            % (rc, "all of them" if got == want else (m.group(0) if m else "no # zonal line")))

    def kv(line):
        return dict(t.split("=", 1) for t in line.split() if "=" in t)
    sched = os.path.join(work, "zonal_saturn.sched")
    kcommon = ["--format", "fp64", "--steps", "4", "--sample", "2", "--quiet"] + adapt
    rc0, out0, err0 = run([cft, "--problem", PROBLEM, "--record-iters", sched] + kcommon)
    key = ""
    if rc0 == 0:
        with open(sched) as f:
            key = next((l[6:].strip() for l in f if l.startswith("# run ")), "")
    zfield = kv(key).get("zonal", "")
    rcs, outs, _ = run([cft, "--problem", PROBLEM, "--replay-iters", sched] + kcommon)
    attempts = sum(1 for l in open(sched) if not l.startswith("#")) if rc0 == 0 else 0
    lines0 = [l for l in out0.splitlines() if l.startswith("sample ")]
    lines_r = [l for l in outs.splitlines() if l.startswith("sample ")]
    verdict(rc0 == 0 and zfield.startswith("0,") and len(zfield.split(",")) == 4 and rcs == 0 and attempts > 0 and lines0
            and lines0 == lines_r,
            "a schedule recorded with the zonal line, %d attempts, its key's zonal=%s; replayed under it: rc %d, %d of %d "
            "sample lines the recording's" % (attempts, zfield or "(none)", rcs,
                                              sum(1 for x, y in zip(lines0, lines_r) if x == y), len(lines0)))

    def nudge(text, idx):
        # the zonal line's field idx (1 R, 2 J2, 3 J4), its last hex digit moved
        def f(mm):
            t = mm.group(0).split()
            mant, ex = t[idx + 1].split("p")
            t[idx + 1] = mant[:-1] + ("0" if mant[-1] != "0" else "1") + "p" + ex
            return " ".join(t)
        return re.sub(r"(?m)^zonal .*$", f, text)
    for label, path in (("J2", write("zonal_j2.txt", nudge(head, 2))), ("J4", write("zonal_j4.txt", nudge(head, 3))),
                        ("R", write("zonal_r.txt", nudge(head, 1))), ("no zonal line", nozonal)):
        rc, out, err = run([cft, "--problem", path, "--replay-iters", sched] + kcommon)
        mm = re.search(r"it was recorded from: (.*)\n\s*this run is:\s+(.*)", err)
        moved = sorted(k for k in set(kv(mm.group(1))) | set(kv(mm.group(2)))
                       if kv(mm.group(1)).get(k) != kv(mm.group(2)).get(k)) if mm else []
        verdict(rc == 3 and "was recorded from a different run" in err and moved == ["zonal"],
                "the schedule replayed under another %s: rc %d, %s; the keys differ in %s (wanted refused by name, the "
                "keys differing in zonal= alone)" % (label, rc, "refused by name" if "was recorded from a different run"
                                                     in err else "not refused by name", ", ".join(moved) or "nothing"))
    refusals = (("a zonal body that is a test particle", head.replace("zonal 0 ", "zonal 2 "), "is a test particle"),
                ("a zonal body past N", head.replace("zonal 0 ", "zonal 3 "), "k is a body of each system"),
                ("a short zonal line", re.sub(r"(?m)^(zonal \S+ \S+ \S+) \S+$", r"\1", head), "a short zonal line"))
    for label, text, needle in refusals:
        rc, out, err = run([cft, "--format", "fp64", "--problem", write("zonal_refused.txt", text), "--steps", "1", "--quiet"])
        verdict(rc != 0 and needle in err and not out.strip(), "%s: rc %d, %s of stdout (wanted a refusal by name, "
                "before any record): %s" % (label, rc, "nothing" if not out.strip() else "%d lines" % len(out.splitlines()),
                                             err.strip()[-150:]))

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
