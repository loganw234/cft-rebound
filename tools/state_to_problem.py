#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Turn one sample of a record into a problem file: the state at that
# sample, bit for bit, as the initial condition of a new run. The record
# holds every coordinate as an exact hex float, and the problem parser
# reads hex floats exactly at any format, so a run started from this
# file continues from precisely the recorded state - which is what a
# reversal test (integrate forward, then back with -dt) needs.
#
#   python tools/state_to_problem.py RECORD --out PROBLEM.txt [--sample K | --last] [--system S] [--name NAME]

import argparse
import sys


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("record")
    ap.add_argument("--out", required=True)
    ap.add_argument("--sample", type=int, default=None)
    ap.add_argument("--last", action="store_true")
    ap.add_argument("--system", type=int, default=None)
    ap.add_argument("--name", default=None)
    args = ap.parse_args()
    name, G, bodies, masses, lines = None, None, [], [], []
    want = None if args.system is None else "system %d sample " % args.system
    cur_sys = None
    for line in open(args.record):
        if line.startswith("# program="):
            for kv in line[2:].split():
                k, _, v = kv.partition("=")
                if k == "problem":
                    name = v
        elif line.startswith("# dt0="):
            toks = line[2:].split()
            for i, t in enumerate(toks):
                if t == "G=":
                    G = toks[i + 1]
        elif line.startswith("# system ") and "steps_done" not in line:
            cur_sys = int(line.split()[2])
        elif line.startswith("# body"):
            if args.system is None or cur_sys == args.system:
                toks = line.split()
                bodies.append(toks[3])
                masses.append(toks[-1])
        elif (want is None and line.startswith("sample ")) or (want is not None and line.startswith(want)):
            lines.append(line.split())
    if not lines:
        raise SystemExit("no sample lines")
    if args.last or args.sample is None:
        toks = lines[-1]
    else:
        toks = next(l for l in lines if int(l[1 if want is None else 3]) == args.sample)
    if want is not None:
        toks = toks[2:]
    step = toks[2]
    body = toks[3:]
    if "|" in body:
        body = body[:body.index("|")]
    state = body[4:]
    N = len(bodies)
    if len(state) != 6 * N:
        raise SystemExit("%d state values for %d bodies" % (len(state), N))
    with open(args.out, "w", newline="\n") as f:
        f.write("# cft-rebound problem file, exact hex floats\n")
        f.write("# the state of %s at step %s of %s, t = %s\n" % (name, step, args.record, toks[3]))
        f.write("name %s\n" % (args.name or (name + "_at" + step)))
        f.write("G %s\n" % G)
        f.write("N %d\n" % N)
        for i in range(N):
            f.write("body %s %s %s\n" % (bodies[i], masses[i], " ".join(state[6 * i:6 * i + 6])))
    print("wrote %s: %d bodies at step %s" % (args.out, N, step))
    return 0


if __name__ == "__main__":
    sys.exit(main())
