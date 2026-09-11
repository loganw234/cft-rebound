# SPDX-License-Identifier: GPL-3.0-or-later
"""The binary64 equivalence gate, from Python.

The same problem, the same fixed step, twice: once through REBOUND's own
``ias15`` and once through the integrator this repository registers. At
binary64 the two must agree to the last bit - that is the gate the C
side runs (tools/check_equivalence.py) and this is the same claim seen
from a ``rebound.Simulation``.

    python python/example_equivalence.py --library build/libcft_ias15.so

(``.dylib`` on macOS, ``.dll`` on Windows - what ``make python-lib``
wrote. ``make check-python`` runs exactly this.)

The initial conditions are read from data/problems/kepler.txt as exact
hex floats, so this run and ``build/ias15_ref`` start from identical
bits.

``--expect-differ`` inverts the verdict, which is what makes this a
gate rather than a demonstration: run it at ``--format fp128`` and the
13 values must MOVE. See ``make check-python-control``.
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


def run(integrator, G, bodies, dt, steps, lib=None, format="fp64"):
    import rebound

    sim = rebound.Simulation()
    sim.G = G
    sim.integrator = integrator
    sim.dt = dt

    # Both runs must have the SAME step control or the comparison is
    # between two different integrations and says nothing about
    # arithmetic. REBOUND's own integrators expose it as an attribute
    # under the name a user knows; this one is configured through the
    # library, because REBOUND resolves sim.integrator.<field> against
    # the registered field descriptor list and every name in ours is
    # cft_-prefixed - so sim.integrator.epsilon does not resolve here.
    if lib is not None:
        cft_rebound.configure(lib, sim, format=format, epsilon=0.0)
        print("  configured through the library: %s, fixed step" % format)
    else:
        # REBOUND's own: sim.integrator is the configuration object for
        # whichever integrator is selected, and ias15's descriptor list
        # spells this field "epsilon". 0.0 is REBOUND's convention for
        # a fixed step.
        sim.integrator.epsilon = 0.0
    for m, x, y, z, vx, vy, vz in bodies:
        sim.add(m=m, x=x, y=y, z=z, vx=vx, vy=vy, vz=vz)
    sim.steps(steps)
    out = [sim.t]
    for p in sim.particles:
        out += [p.x, p.y, p.z, p.vx, p.vy, p.vz]
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--library",
                    default=(os.environ.get("CFT_REBOUND_LIB")
                             or cft_rebound.default_library()),
                    help="the shared library that registers the integrator. "
                         "Defaults to $CFT_REBOUND_LIB, then to whatever "
                         "`make python-lib` left in build/.")
    ap.add_argument("--integrator", default="ias15_cft")
    ap.add_argument("--format", default="fp64",
                    choices=("fp64", "fp128", "fp256"),
                    help="the wide format for the registered integrator. "
                         "Only fp64 can be identical to REBOUND; the wider "
                         "ones are expected to differ, and by how much is "
                         "the point of running them.")
    ap.add_argument("--problem",
                    default=os.path.join(ROOT, "data", "problems", "kepler.txt"))
    ap.add_argument("--dt", default="0.05")
    ap.add_argument("--steps", type=int, default=252)
    ap.add_argument("--expect-differ", action="store_true",
                    help="invert the verdict: succeed only if the two runs "
                         "DIFFER. This is the control for the gate. At "
                         "--format fp128 the arithmetic has to reach the "
                         "integrator, so 13 identical values would mean the "
                         "format was accepted and ignored - and that the "
                         "default run is binary64 against binary64, passing "
                         "against itself.")
    args = ap.parse_args()

    if not args.library:
        ap.error("give --library (or set CFT_REBOUND_LIB): the library built "
                 "from this repository that registers %r" % args.integrator)

    import rebound
    print("rebound %s (%s)" % (rebound.__version__, rebound.__githash__))
    print("librebound %s" % cft_rebound.library_path())

    lib = cft_rebound.load(args.library)
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
    b = run(args.integrator, G, bodies, dt, args.steps, lib=lib, format=args.format)

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
        if args.expect_differ:
            print("\nBut --expect-differ was given, so this is a FAILURE: at "
                  "%s the arithmetic must reach the integrator and move the "
                  "answer. Identical values mean the format was accepted and "
                  "ignored, which would also mean the plain run compares "
                  "binary64 against binary64 and passes against itself."
                  % args.format)
            return 1
        return 0
    print("DIFFERS: %d of %d values, largest absolute difference %.3e"
          % (differing, len(a), worst))
    if args.expect_differ:
        print("\nExpected, and this is the control passing: %s changed the "
              "answer, so the format reaches the arithmetic and the binary64 "
              "comparison is a real one." % args.format)
        return 0
    print("At binary64 this should be zero. A non-zero result is the gate "
          "failing, not a tolerance to widen.")
    return 1


if __name__ == "__main__":
    sys.exit(main())
