#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Turn the sweeps' per-sample CSVs and summary lines into the tables
# and the figure the documentation quotes.
#
#   python tools/tabulate.py [--results results] [--png results/energy_vs_orbits.png]
#
# Reads results/summary.txt (one line per finished job; a job that ran
# twice keeps its last line), results/<job>.csv (the oracle's
# per-sample table) and results/floor_*.csv (compare_formats output).

import argparse
import csv
import glob
import os
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
FORMATS = ("fp64", "fp128", "fp256")


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
        if "max_dE" in kv:
            out[name] = kv
    return out


def e(x):
    if x is None:
        return "-"
    return "%.2e" % x if x else "0"


class Runs:
    def __init__(self, R):
        self.summ = read_summary(os.path.join(R, "summary.txt"))
        self.csv = {}
        for p in sorted(glob.glob(os.path.join(R, "*.csv"))):
            name = os.path.basename(p)[:-4]
            if not name.startswith("floor_"):
                self.csv[name] = read_csv(p)

    def stat(self, name):
        rr = self.csv.get(name)
        s = self.summ.get(name)
        if not rr or not s:
            return None
        last = rr[-1]
        return {"orbits": last["orbits"], "steps": last["step"],
                "maxE": max(r["rel_energy_err"] for r in rr), "lastE": last["rel_energy_err"],
                "maxL": max(r["rel_angmom_err"] for r in rr), "lastX": last["rel_pos_err"],
                "pc": s.get("mean_pc", "-"), "sec": s.get("seconds", "-"), "rej": s.get("rejected", "-")}


def fixed_table(runs, title, rows, unit):
    print("### %s\n" % title)
    print("| dt | steps | %s | max dE/E fp64 | fp128 | fp256 | final dx/a fp64 | fp128 | fp256 | corrector passes fp64/fp128/fp256 | seconds fp64/fp128/fp256 |" % unit)
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for label, names in rows:
        st = [runs.stat(n) for n in names]
        if not any(st):
            continue
        ref = next(s for s in st if s)
        print("| %s | %d | %.1f | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            label, ref["steps"], ref["orbits"],
            e(st[0]["maxE"]) if st[0] else "-", e(st[1]["maxE"]) if st[1] else "-", e(st[2]["maxE"]) if st[2] else "-",
            e(st[0]["lastX"]) if st[0] else "-", e(st[1]["lastX"]) if st[1] else "-", e(st[2]["lastX"]) if st[2] else "-",
            "/".join(s["pc"] if s else "-" for s in st),
            "/".join(("%.0f" % float(s["sec"])) if s and s["sec"] != "-" else "-" for s in st)))
    print()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default=os.path.join(root, "results"))
    ap.add_argument("--png", default=None)
    args = ap.parse_args()
    R = args.results
    runs = Runs(R)

    fixed_table(runs, "Kepler (e = 1/2), fixed step", [
        ("0.1", ["kepler_fixed_dt0.1_%s" % f for f in FORMATS]),
        ("1/16", ["kepler_fixed_dt0.0625_%s" % f for f in FORMATS]),
        ("0.05", ["kepler_fixed_dt0.05_%s" % f for f in FORMATS]),
        ("0.025", ["kepler_fixed_dt0.025_%s" % f for f in FORMATS]),
        ("1/64", ["kepler_fixed_dt0.015625_%s" % f for f in FORMATS]),
        ("1/128", ["kepler_fixed_dt0.0078125_%s" % f for f in FORMATS]),
    ], "orbits")

    print("### Kepler, adaptive step (PRS23)\n")
    print("| epsilon | steps | orbits | max dE/E fp64 | fp128 | fp256 | final dx/a fp64 | fp128 | fp256 | corrector passes | seconds |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for eps in ("1e-9", "1e-12", "1e-16"):
        st = [runs.stat("kepler_adaptive_eps%s_%s" % (eps, f)) for f in FORMATS]
        if not any(st):
            continue
        ref = next(s for s in st if s)
        print("| %s | %d | %.1f | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            eps, ref["steps"], ref["orbits"],
            *[e(s["maxE"]) if s else "-" for s in st], *[e(s["lastX"]) if s else "-" for s in st],
            "/".join(s["pc"] if s else "-" for s in st),
            "/".join(("%.0f" % float(s["sec"])) if s and s["sec"] != "-" else "-" for s in st)))
    print()

    print("### Outer solar system, fixed step\n")
    print("| dt (days) | steps | years | max dE/E fp64 | fp128 | fp256 | max dL/L fp64 | fp128 | fp256 | corrector passes | seconds |")
    print("|---|---|---|---|---|---|---|---|---|---|---|")
    for label, names in (("40", ["outer_fixed_dt40_%s" % f for f in FORMATS]),
                         ("20", ["outer_fixed_dt20_%s" % f for f in FORMATS]),
                         ("20, 300 y", ["outer_fixed_dt20_%s_300y" % f for f in FORMATS]),
                         ("10, 300 y", ["outer_fixed_dt10_%s_300y" % f for f in FORMATS])):
        st = [runs.stat(n) for n in names]
        if not any(st):
            continue
        ref = next(s for s in st if s)
        print("| %s | %d | %.0f | %s | %s | %s | %s | %s | %s | %s | %s |" % (
            label, ref["steps"], ref["orbits"],
            *[e(s["maxE"]) if s else "-" for s in st], *[e(s["maxL"]) if s else "-" for s in st],
            "/".join(s["pc"] if s else "-" for s in st),
            "/".join(("%.0f" % float(s["sec"])) if s and s["sec"] != "-" else "-" for s in st)))
    print()

    print("### Compensated summation: Kahan against 9.5 augmentedAddition, Kepler dt = 0.05, 200 orbits\n")
    print("| format | cs | max dE/E | final dE/E | max dL/L | final dx/a |")
    print("|---|---|---|---|---|---|")
    for f in ("fp64", "fp256"):
        for suffix, cs in (("", "kahan"), ("_aug", "augmented")):
            s = runs.stat("kepler_fixed_dt0.05_%s%s" % (f, suffix))
            if s:
                print("| %s | %s | %s | %s | %s | %s |" % (f, cs, e(s["maxE"]), e(s["lastE"]), e(s["maxL"]), e(s["lastX"])))
    print()

    print("### Round-off floors measured directly (narrow record minus wide record, same dyadic step)\n")
    print("| pair | steps | max dE/E | final dE/E | max dx/L | final dx/L |")
    print("|---|---|---|---|---|---|")
    for p in sorted(glob.glob(os.path.join(R, "floor_*.csv"))):
        rows = read_csv(p)
        if not rows:
            continue
        name = os.path.basename(p)[6:-4]
        print("| %s | %d | %s | %s | %s | %s |" % (name, rows[-1]["step"],
              e(max(r["rel_energy_diff"] for r in rows)), e(rows[-1]["rel_energy_diff"]),
              e(max(r["max_rel_coord_diff"] for r in rows)), e(rows[-1]["max_rel_coord_diff"])))
    print()

    if args.png:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, axes = plt.subplots(1, 2, figsize=(13, 5))
        colors = {"fp64": "#1f77b4", "fp128": "#2ca02c", "fp256": "#d62728"}
        ax = axes[0]
        for dt, style in (("0.1", ":"), ("0.05", "-"), ("0.025", "--"), ("0.0625", "-."), ("0.015625", (0, (1, 1))), ("0.0078125", (0, (5, 1, 1, 1)))):
            for f in FORMATS:
                rr = runs.csv.get("kepler_fixed_dt%s_%s" % (dt, f))
                if rr:
                    ax.plot([r["orbits"] for r in rr[1:]], [max(r["rel_energy_err"], 1e-80) for r in rr[1:]],
                            linestyle=style, color=colors[f], label="%s dt=%s" % (f, dt), linewidth=1.2)
        ax.set_yscale("log"); ax.set_xscale("log"); ax.set_xlabel("orbits"); ax.set_ylabel("|dE/E| against the exact two-body solution")
        ax.set_title("Kepler e=1/2, fixed step: IAS15 through libcft"); ax.legend(fontsize=6, ncol=3)
        ax = axes[1]
        for dt, style in (("40", "-"), ("20", "--")):
            for f in FORMATS:
                for suffix in ("", "_300y"):
                    rr = runs.csv.get("outer_fixed_dt%s_%s%s" % (dt, f, suffix))
                    if rr:
                        ax.plot([r["orbits"] for r in rr[1:]], [max(r["rel_energy_err"], 1e-80) for r in rr[1:]],
                                linestyle=style, color=colors[f], label="%s dt=%s d%s" % (f, dt, suffix.replace("_", " ")), linewidth=1.2)
        rr = runs.csv.get("outer_fixed_dt10_fp128_300y"); rr2 = runs.csv.get("outer_fixed_dt10_fp256_300y")
        for f, r in (("fp128", rr), ("fp256", rr2)):
            if r:
                ax.plot([x["orbits"] for x in r[1:]], [max(x["rel_energy_err"], 1e-80) for x in r[1:]], linestyle=":", color=colors[f], label="%s dt=10 d 300y" % f, linewidth=1.2)
        ax.set_yscale("log"); ax.set_xscale("log"); ax.set_xlabel("years"); ax.set_ylabel("|dE/E|")
        ax.set_title("outer solar system (Sun, 4 giants, Pluto), fixed step"); ax.legend(fontsize=6, ncol=2)
        fig.tight_layout()
        fig.savefig(args.png, dpi=110)
        print("wrote %s" % args.png)
    return 0


if __name__ == "__main__":
    sys.exit(main())
