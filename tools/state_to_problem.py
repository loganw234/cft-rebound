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
# The problem is the record's, every line of it: name, G, N and the
# bodies, and the `active` and `zonal` lines (ROUND4.md, P0) from the
# record's `# active` and `# zonal` lines, which the program writes only
# for a problem that has them, so a record without them gives the file it
# always gave. Until VALIDATION entry 52 this tool dropped those two lines,
# and tools/reversal.py ran a zonal or `active` problem's backward leg on
# another problem. A problem line added later must be carried here too:
# reversal.py refuses a turn problem whose lines are not the original's,
# and tools/check_reversal.py holds both tools to it.
#
#   python tools/state_to_problem.py RECORD --out PROBLEM.txt [--sample K | --last] [--system S] [--name NAME]

import argparse
import re
import sys

# The record's two problem lines, as src/ias15_cft.c writes them (the
# header, and zonal_record() in src/ias15_zonal.inc): each value is the
# run's own exact hex, after a "KEY= ".
ACTIVE = re.compile(r"# active (\d+) of (\d+) bodies a system:")
ZONAL = re.compile(r"# zonal body (\d+) of every system, about z: R= (\S+) J2= (\S+) J4= ([^\s;]+);")


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
    active, zonal = None, None
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
        elif line.startswith("# active "):
            m = ACTIVE.match(line)
            if not m or active is not None:
                raise SystemExit("state_to_problem: %s: %s `# active` line: %s"
                                 % (args.record, "a second" if m else "an unreadable", line.strip()))
            active = m.groups()
        elif line.startswith("# zonal "):
            m = ZONAL.match(line)
            if not m or zonal is not None:
                raise SystemExit("state_to_problem: %s: %s `# zonal` line: %s"
                                 % (args.record, "a second" if m else "an unreadable", line.strip()))
            zonal = m.groups()
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
    # the two lines name bodies of each system: hold them to this one's N
    if active is not None and (int(active[1]) != N or not 1 <= int(active[0]) <= N):
        raise SystemExit("state_to_problem: %s: `# active %s of %s` for %d bodies a system" % (args.record, active[0], active[1], N))
    if zonal is not None and not int(zonal[0]) < N:
        raise SystemExit("state_to_problem: %s: `# zonal body %s` for %d bodies a system" % (args.record, zonal[0], N))
    with open(args.out, "w", newline="\n") as f:
        f.write("# cft-rebound problem file, exact hex floats\n")
        f.write("# the state of %s at step %s of %s, t = %s\n" % (name, step, args.record, toks[3]))
        f.write("name %s\n" % (args.name or (name + "_at" + step)))
        f.write("G %s\n" % G)
        f.write("N %d\n" % N)
        if active is not None:
            f.write("active %s\n" % active[0])
        if zonal is not None:
            f.write("zonal %s %s %s %s\n" % zonal)
        for i in range(N):
            f.write("body %s %s %s\n" % (bodies[i], masses[i], " ".join(state[6 * i:6 * i + 6])))
    carried = "".join(("" if active is None else ", active %s" % active[0],
                       "" if zonal is None else ", zonal body %s" % zonal[0]))
    print("wrote %s: %d bodies at step %s%s" % (args.out, N, step, carried))
    return 0


if __name__ == "__main__":
    sys.exit(main())
