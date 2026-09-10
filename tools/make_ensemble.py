#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Write an ensemble file: E copies of a base problem, member k
# perturbed by a rule, every value an exact binary64 bit pattern so
# that the ensemble run, the members run alone (--member k) and every
# format integrate literally the same initial state.
#
#   python tools/make_ensemble.py BASE.txt --members E --out ENS.txt
#          [--offset BODY COORD DELTA]   member k adds k*DELTA to that coordinate
#          [--vscale BODY FACTOR]        member k multiplies that body's velocity by (1 + k*FACTOR)
#          [--geometric BODY COORD DELTA] member 0 untouched, member k adds DELTA*2^(k-1)
#          [--name NAME]
#
# COORD is one of x y z vx vy vz. Python floats are binary64 with
# round-to-nearest-even, so k*DELTA and the sum round exactly as a C
# double would; with a dyadic DELTA the offset itself is exact and only
# the sum rounds. The file format is the problem format plus `E` and
# `system NAME` lines; a file with no `system` line is one system.

import argparse
import os
import sys

COORDS = {"x": 2, "y": 3, "z": 4, "vx": 5, "vy": 6, "vz": 7}


def read_problem(path):
    name, G, N, bodies = None, None, None, []
    for line in open(path):
        if line.startswith("#") or not line.strip():
            continue
        toks = line.split()
        if toks[0] == "name":
            name = toks[1]
        elif toks[0] == "G":
            G = toks[1]
        elif toks[0] == "N":
            N = int(toks[1])
        elif toks[0] == "body":
            bodies.append([toks[1], float.fromhex(toks[2])] + [float.fromhex(t) for t in toks[3:9]])
        elif toks[0] in ("system", "E"):
            raise SystemExit("%s is already an ensemble file" % path)
        else:
            raise SystemExit("unknown line %r in %s" % (line, path))
    if N is None:
        N = len(bodies)
    if N != len(bodies):
        raise SystemExit("N=%d but %d bodies in %s" % (N, len(bodies), path))
    return name, G, bodies


def num(s):
    """A decimal or a hex float (0x1p-40) as a binary64."""
    return float.fromhex(s) if s.lower().startswith(("0x", "-0x")) else float(s)


def body_index(bodies, key):
    if key.isdigit():
        return int(key)
    for i, b in enumerate(bodies):
        if b[0] == key:
            return i
    raise SystemExit("no body %r" % key)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("base")
    ap.add_argument("--members", type=int, required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--name", default=None)
    ap.add_argument("--offset", nargs=3, action="append", default=[], metavar=("BODY", "COORD", "DELTA"))
    ap.add_argument("--vscale", nargs=2, action="append", default=[], metavar=("BODY", "FACTOR"))
    ap.add_argument("--geometric", nargs=3, action="append", default=[], metavar=("BODY", "COORD", "DELTA"))
    ap.add_argument("--members-dir", default=None, help="also write each member as a single-system problem file m<k>.txt here (for ias15_ref, which reads no ensemble files)")
    args = ap.parse_args()
    name, G, bodies = read_problem(args.base)
    E = args.members
    if E < 1:
        raise SystemExit("--members must be at least 1")
    ename = args.name or (name + "_ens")
    rules = []
    with open(args.out, "w", newline="\n") as f:
        f.write("# cft-rebound ensemble file, binary64 hex floats\n")
        f.write("# base %s, %d members\n" % (os.path.basename(args.base), E))
        for body, coord, delta in args.offset:
            rules.append("member k: %s.%s += k * %r" % (body, coord, num(delta)))
        for body, factor in args.vscale:
            rules.append("member k: %s.v *= 1 + k * %r" % (body, num(factor)))
        for body, coord, delta in args.geometric:
            rules.append("member k > 0: %s.%s += %r * 2^(k-1)" % (body, coord, num(delta)))
        for r in rules:
            f.write("# %s\n" % r)
        f.write("name %s\n" % ename)
        f.write("G %s\n" % G)
        f.write("N %d\n" % len(bodies))
        f.write("E %d\n" % E)
        if args.members_dir:
            os.makedirs(args.members_dir, exist_ok=True)
        for k in range(E):
            f.write("system m%d\n" % k)
            mf = open(os.path.join(args.members_dir, "m%d.txt" % k), "w", newline="\n") if args.members_dir else None
            if mf:
                mf.write("# member %d of %s\nname %s_m%d\nG %s\nN %d\n" % (k, os.path.basename(args.out), ename, k, G, len(bodies)))
            for bi, b in enumerate(bodies):
                b = list(b)
                for body, coord, delta in args.offset:
                    if body_index(bodies, body) == bi:
                        b[COORDS[coord]] = b[COORDS[coord]] + k * num(delta)
                for body, factor in args.vscale:
                    if body_index(bodies, body) == bi:
                        s = 1.0 + k * num(factor)
                        for c in (5, 6, 7):
                            b[c] = b[c] * s
                for body, coord, delta in args.geometric:
                    if k > 0 and body_index(bodies, body) == bi:
                        b[COORDS[coord]] = b[COORDS[coord]] + num(delta) * (2.0 ** (k - 1))
                f.write("body %s %s\n" % (b[0], " ".join(float(v).hex() for v in b[1:])))
                if mf:
                    mf.write("body %s %s\n" % (b[0], " ".join(float(v).hex() for v in b[1:])))
            if mf:
                mf.close()
    print("wrote %s: %d systems of %d bodies%s" % (args.out, E, len(bodies), ("; " + "; ".join(rules)) if rules else ""))
    return 0


if __name__ == "__main__":
    sys.exit(main())
