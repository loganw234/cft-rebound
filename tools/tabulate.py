#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Turn the sweep's per-sample CSVs into the tables and the figure the
# documentation quotes: energy error against orbit count per format
# and step size, for the fixed-step Kepler runs, the adaptive runs and
# the outer solar system.
#
#   python tools/tabulate.py [--results results] [--png results/energy_vs_orbits.png]

import argparse
import csv
import glob
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)


def read_csv(path):
    rows = []
    with open(path) as f:
        for r in csv.DictReader(f):
            rows.append({k: (float(v) if k not in ("sample", "step") else int(v)) for k, v in r.items()})
    return rows


def read_summary(path):
    out = {}
    if not os.path.exists(path):
        return out
    for line in open(path):
        name, _, rest = line.partition(":")
        name = name.strip()
        kv = {}
        for tok in rest.split():
            k, _, v = tok.partition("=")
            if v:
                kv[k] = v
        out[name] = kv
    return out


def fmt_e(x):
    return "%.2e" % x if x else "0"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(root, "results"))
    ap.add_argument("--png", default=None)
    args = ap.parse_args()
    R = args.results
    summ = read_summary(os.path.join(R, "summary.txt"))
    runs = {}
    for p in sorted(glob.glob(os.path.join(R, "*.csv"))):
        name = os.path.basename(p)[:-4]
        runs[name] = read_csv(p)

    def row(name):
        rr = runs.get(name)
        if not rr:
            return None
        last = rr[-1]
        maxE = max(r["rel_energy_err"] for r in rr)
        maxL = max(r["rel_angmom_err"] for r in rr)
        s = summ.get(name, {})
        return (last["orbits"], last["step"], maxE, last["rel_energy_err"], maxL, last["rel_pos_err"],
                s.get("mean_pc", "-"), s.get("rejected", "-"), s.get("seconds", "-"), s.get("wall", "-"))

    print("## Kepler, fixed step, 200 orbits\n")
    print("| dt | steps | format | max |dE/E| | final |dE/E| | max |dL/L| | final |dx|/a | corrector passes | seconds |")
    print("|---|---|---|---|---|---|---|---|---|")
    for dt in ("0.1", "0.05", "0.025"):
        for f in ("fp64", "fp128", "fp256"):
            r = row("kepler_fixed_dt%s_%s" % (dt, f))
            if r:
                print("| %s | %d | %s | %s | %s | %s | %s | %s | %s |" % (dt, r[1], f, fmt_e(r[2]), fmt_e(r[3]), fmt_e(r[4]), fmt_e(r[5]), r[6], r[8]))
    print("\n## Kepler, fixed step dt=0.05, augmented summation\n")
    print("| format | cs | max |dE/E| | final |dE/E| | max |dL/L| | final |dx|/a |")
    print("|---|---|---|---|---|---|")
    for f in ("fp64", "fp256"):
        for suffix, cs in (("", "kahan"), ("_aug", "augmented")):
            r = row("kepler_fixed_dt0.05_%s%s" % (f, suffix))
            if r:
                print("| %s | %s | %s | %s | %s | %s |" % (f, cs, fmt_e(r[2]), fmt_e(r[3]), fmt_e(r[4]), fmt_e(r[5])))
    print("\n## Kepler, adaptive step\n")
    print("| epsilon | format | steps | orbits | max |dE/E| | final |dE/E| | final |dx|/a | corrector passes | rejected | seconds |")
    print("|---|---|---|---|---|---|---|---|---|---|")
    for eps in ("1e-9", "1e-12", "1e-16"):
        for f in ("fp64", "fp128", "fp256"):
            r = row("kepler_adaptive_eps%s_%s" % (eps, f))
            if r:
                print("| %s | %s | %d | %.1f | %s | %s | %s | %s | %s | %s |" % (eps, f, r[1], r[0], fmt_e(r[2]), fmt_e(r[3]), fmt_e(r[5]), r[6], r[7], r[8]))
    print("\n## Outer solar system, fixed step, 1200 years\n")
    print("| dt (days) | steps | format | max |dE/E| | final |dE/E| | max |dL/L| | corrector passes | seconds |")
    print("|---|---|---|---|---|---|---|---|")
    for dt in ("40", "20"):
        for f in ("fp64", "fp128", "fp256"):
            r = row("outer_fixed_dt%s_%s" % (dt, f))
            if r:
                print("| %s | %d | %s | %s | %s | %s | %s | %s |" % (dt, r[1], f, fmt_e(r[2]), fmt_e(r[3]), fmt_e(r[4]), r[6], r[8]))

    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
        colors = {"fp64": "#1f77b4", "fp128": "#2ca02c", "fp256": "#d62728"}
        styles = {"0.1": ":", "0.05": "-", "0.025": "--"}
        ax = axes[0]
        for dt in ("0.1", "0.05", "0.025"):
            for f in ("fp64", "fp128", "fp256"):
                rr = runs.get("kepler_fixed_dt%s_%s" % (dt, f))
                if rr:
                    xs = [r["orbits"] for r in rr[1:]]
                    ys = [max(r["rel_energy_err"], 1e-80) for r in rr[1:]]
                    ax.plot(xs, ys, styles[dt], color=colors[f], label="%s dt=%s" % (f, dt))
        ax.set_yscale("log"); ax.set_xlabel("orbits"); ax.set_ylabel("|dE/E|")
        ax.set_title("Kepler e=1/2, fixed step, IAS15 through libcft"); ax.legend(fontsize=7, ncol=3)
        ax = axes[1]
        for dt in ("40", "20"):
            for f in ("fp64", "fp128", "fp256"):
                rr = runs.get("outer_fixed_dt%s_%s" % (dt, f))
                if rr:
                    xs = [r["orbits"] for r in rr[1:]]
                    ys = [max(r["rel_energy_err"], 1e-80) for r in rr[1:]]
                    ax.plot(xs, ys, "-" if dt == "40" else "--", color=colors[f], label="%s dt=%s d" % (f, dt))
        ax.set_yscale("log"); ax.set_xlabel("years"); ax.set_ylabel("|dE/E|")
        ax.set_title("outer solar system, fixed step"); ax.legend(fontsize=7, ncol=2)
        fig.tight_layout()
        fig.savefig(args.png, dpi=110)
        print("\nwrote %s" % args.png)
    return 0


if __name__ == "__main__":
    sys.exit(main())
