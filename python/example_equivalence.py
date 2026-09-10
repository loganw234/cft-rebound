# SPDX-License-Identifier: GPL-3.0-or-later
"""The binary64 equivalence gate, from Python.

The same problem, the same fixed step, twice: once through REBOUND's own
``ias15`` and once through the integrator this repository registers. At
binary64 the two must agree to the last bit - that is the gate the C
side runs (tools/check_equivalence.py) and this is the same claim seen
from a ``rebound.Simulation``.

    python python/example_equivalence.py --library build/libias15_cft.so

The initial conditions are read from data/problems/kepler.txt as exact
hex floats, so this run and ``build/ias15_ref`` start from identical
bits.
"""

import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, HERE)

import cft_rebound  # noqa: E402


def read_problem(path):
    """A cft-rebound problem file: G, then one 'body' line per particle."""
    G = 1.0
    bodies = []
    with open(path) as f:
        for line in f:
            tok = line.split()
            if not tok or tok[0].startswith("#"):
                continue
            if tok[0] == "G":
                G = float.fromhex(tok[1])
            elif tok[0] == "body":
                # body NAME m x y z vx vy vz
                m, x, y, z, vx, vy, vz = (float.fromhex(t) for t in tok[2:9])
                bodies.append((m, x, y, z, vx, vy, vz))
    return G, bodies


def run(integrator, G, bodies, dt, steps):
    import rebound

    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = integrator
    sim.dt = dt
    try:
        sim.integrator.epsilon = 0.0  # REBOUND's convention for a fixed step
    except AttributeError:
        print("  note: %r has no 'epsilon' field; step control is its own"
              % str(sim.integrator))
    for m, x, y, z, vx, vy, vz in bodies:
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    sim.steps(steps)
    out = [sim.t]
    for p in sim.particles:
        out += [p.x, p.y, p.z, p.vx, p.vy, p.vz]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library", default=os.environ.get("CFT_REBOUND_LIB"),
                    help="the shared library that registers the integrator")
    ap.add_argument("--integrator", default="ias15_cft")
    ap.add_argument("--problem",
                    default=os.path.join(ROOT, "data", "problems", "kepler.txt"))
    ap.add_argument("--dt", default="0.05")
    ap.add_argument("--steps", type=int, default=252)
    args = ap.parse_args()

    if not args.library:
        ap.error("give --library (or set CFT_REBOUND_LIB): the library built "
                 "from this repository that registers %r" % args.integrator)

    import rebound
    print("rebound %s (%s)" % (rebound.__version__, rebound.__githash__))
    print("librebound %s" % cft_rebound.library_path())

    cft_rebound.load(args.library)
    names = cft_rebound.registered()
    print("loaded %s; registered integrators: %s" % (args.library, " ".join(names)))
    if args.integrator not in names:
        print("\n%r is not registered. This example needs the library that "
              "registers it." % args.integrator)
        return 2
    if args.integrator != args.integrator.lower():
        print("\n%r has an upper-case letter. REBOUND's Python layer "
              "lowercases the name before it reaches C, so it can never be "
              "selected from Python." % args.integrator)
        return 2

    dt = float.fromhex(args.dt) if args.dt.startswith(("0x", "-0x")) else float(args.dt)
    G, bodies = read_problem(args.problem)
    print("problem %s: N = %d, G = %s, dt = %s, %d steps\n"
          % (os.path.basename(args.problem), len(bodies), G.hex(), dt.hex(), args.steps))

    print("REBOUND's own ias15:")
    a = run("ias15", G, bodies, dt, args.steps)
    print("  %s:" % args.integrator)
    b = run(args.integrator, G, bodies, dt, args.steps)

    labels = ["t"] + ["p%d.%s" % (i, c) for i in range(len(bodies))
                      for c in ("x", "y", "z", "vx", "vy", "vz")]
    worst = 0.0
    differing = 0
    print("\n%-8s %-24s %-24s" % ("field", "ias15", args.integrator))
    for name, u, v in zip(labels, a, b):
        same = (u == v)
        if not same:
            differing += 1
            worst = max(worst, abs(u - v))
        print("%-8s %-24s %-24s %s" % (name, u.hex(), v.hex(), "" if same else "<-- differs"))

    print()
    if differing == 0:
        print("IDENTICAL: %d values, bit for bit." % len(a))
        return 0
    print("DIFFERS: %d of %d values, largest absolute difference %.3e"
          % (differing, len(a), worst))
    print("At binary64 this should be zero. A non-zero result is the gate "
          "failing, not a tolerance to widen.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
