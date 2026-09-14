#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
"""A real REBOUND workload across every operating mode, priced and scored.

Not a microbenchmark. Every row here is an actual IAS15 integration of an
actual problem, run by a real program, and the baseline row is REBOUND's
own IAS15 in plain double - the thing a person runs today.

    python3 hw/bench-workload.py --image "1t-f128=$HOME/cardday-f128s/cft_hw_f128_1x.xclbin" \
                                 --image "4t-full=$HOME/cardday-rev4/cft_hw_quad.xclbin" \
                                 --image "6t-f128=$HOME/cardday-f128x6b/cft_hw_f128_6x.xclbin"

THE MODES

    ref/fp64        build/ias15_ref - REBOUND's IAS15, hardware double.
                    The traditional option, and the speed nothing here beats.
    sw/fp64         the port on libcft's softfloat. Bit-for-bit REBOUND
                    (gate 1), so this row prices the CONTRACT, not the physics.
    sw/fp128        the precision people actually want, in software. This is
                    the row the card has to beat to be worth owning.
    <img>/fp64      the port on a tile.
    <img>/fp128     the same, at binary128.

WHY BOTH A CLOCK AND AN ERROR COLUMN. Speed alone cannot answer "is this
worth it", because the rows do not compute the same thing: binary64
round-off floors the energy error near 1e-15 no matter how small the step,
and binary128 does not. A row that is ten times slower and nineteen orders
more accurate is not a loss. So every row carries the relative energy
drift it actually achieved, and the comparison a reader should make is
down a column of EQUAL accuracy - which is why sw/fp128 and card/fp128 sit
on the same chart.

THE CORRECTNESS CHECK IS FREE AND IS TAKEN. ref/fp64 and sw/fp64 must
produce byte-identical records - that is this repository's gate 1 - and so
must every card row at the same format. Any row whose record differs from
its software counterpart is marked and the run is a failure, because a
faster wrong answer is not a result.

The corrector's mean pass count is reported per row for the reason
docs/VALIDATION.md keeps repeating: identical bits prove the backends
agree, never that the problem posed was the one meant.
"""
import argparse
import os
import pathlib
from fractions import Fraction
import re
import shutil
import subprocess
import sys
import time

REPO = pathlib.Path(__file__).resolve().parent.parent


def run_timed(cmd, env=None):
    """-> (seconds, returncode, stdout)"""
    t0 = time.monotonic()
    p = subprocess.run(cmd, capture_output=True, text=True, env=env)
    return time.monotonic() - t0, p.returncode, p.stdout + p.stderr


def data_only(record: str) -> str:
    """The values the equivalence claim is actually about, for comparing
    across PROGRAMS - defined the way tools/check_equivalence.py defines
    it, not a fresh way.

    Two things are excluded and both would otherwise report a divergence
    that is not one. The header lines differ because the programs differ
    (one carries REBOUND's settings, the other cs, arith, engine and the
    tolerance shift). And everything after the `|` is the exact-time
    pair, a double-double decomposition of the instant the sample was
    taken: ias15_ref and ias15_cft split the same time differently
    between the high and low word, so those two fields differ while `t`
    itself and every coordinate are bit-identical. The gate's own parser
    truncates at the `|` for exactly this reason."""
    out = []
    for line in record.splitlines():
        if not line.startswith("sample "):
            continue
        toks = line.split()
        if "|" in toks:
            toks = toks[:toks.index("|")]
        out.append(" ".join(toks))
    return "\n".join(out)


def normalise(record: str) -> str:
    """A record with the two fields that cannot match across runs removed:
    which program and backend produced it, and how long it took. For
    comparing two runs of the SAME program, where the counters are
    comparable and worth keeping in the comparison."""
    out = []
    for line in record.splitlines():
        line = re.sub(r"backend=\S+", "backend=X", line)
        line = re.sub(r"program=\S+", "program=X", line)
        line = re.sub(r"impl=\S+", "impl=X", line)
        line = re.sub(r" seconds=[0-9.]+", "", line)
        line = re.sub(r" steps_per_s=[0-9.]+", "", line)
        # The instrumented trailer's own timings (2026-09-14) and the
        # ensemble's rate: all wall-clock, none of them a value.
        line = re.sub(r" (wall_seconds|t_program|t_elem|t_divsqrt|t_gravity|system_steps_per_s)=[0-9.]+", "", line)
        out.append(line)
    return "\n".join(out)


def energy_drift(record: str):
    """Relative |E - E0| / |E0| from the record's own energy column.

    The record's columns are `sample step t dt_next dt_last E ...`, so the
    energy is field 6 (zero-based) of every data line, written as a hex
    float. Taken from the FIRST and LAST samples of the same run, so it is
    this integration's own drift and not a comparison against a constant
    somebody typed."""
    vals = []
    for line in record.splitlines():
        if not line.startswith("sample "):
            continue
        f = line.split()
        if len(f) > 6:
            try:
                vals.append(hex_fraction(f[6]))
            except ValueError:
                pass
    if len(vals) < 2 or vals[0] == 0:
        return None
    # EXACT rational arithmetic, converted to float only for printing.
    # float.fromhex rounds a binary128 energy to a double, and two values
    # that differ at the 1e-30 level then round to the SAME double, so
    # the drift reads exactly zero - which the first run of this table
    # printed for every binary128 row. The accuracy column is half the
    # point of the table; it cannot be measured with the instrument
    # whose limitation it exists to show.
    return float(abs(vals[-1] - vals[0]) / abs(vals[0]))


def hex_fraction(tok: str):
    """A hex float of any width as an exact Fraction, the way
    tools/check_equivalence.py reads them."""
    s = tok.strip()
    neg = s.startswith("-")
    if neg:
        s = s[1:]
    if not s.lower().startswith("0x"):
        raise ValueError(tok)
    body = s[2:]
    mant, _, exp = body.lower().partition("p")
    ip, _, fp = mant.partition(".")
    digits = (ip + fp) or "0"
    val = Fraction(int(digits, 16), 16 ** len(fp)) * (Fraction(2) ** int(exp or "0"))
    return -val if neg else val


def passes(record: str):
    m = re.search(r"mean_pc_iterations=([0-9.]+)", record)
    return m.group(1) if m else "?"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--image", action="append", default=[],
                    metavar="LABEL=PATH", help="repeatable")
    ap.add_argument("--formats", default="fp64,fp128")
    ap.add_argument("--csv", default=str(REPO / "bench-modes" / "workload.csv"))
    ap.add_argument("--build", default=str(REPO / "build"))
    ap.add_argument("--quick", action="store_true",
                    help="the small problems only")
    a = ap.parse_args()

    build = pathlib.Path(a.build)
    ref, cft = build / "ias15_ref", build / "ias15_cft"
    for p in (ref, cft):
        if not p.exists():
            print(f"missing {p}; run make first", file=sys.stderr)
            return 2

    images = []
    for spec in a.image:
        if "=" not in spec:
            print(f"--image wants LABEL=PATH, got {spec}", file=sys.stderr)
            return 2
        label, path = spec.split("=", 1)
        if not os.path.exists(path):
            print(f"no such image: {path}", file=sys.stderr)
            return 2
        info = subprocess.run([str(cft), "--probe", "--format", "fp64",
                               "--artifact", path], capture_output=True, text=True).stdout
        tiles = re.search(r"^tiles: (\d+)", info, re.M)
        fmts = re.search(r"^formats: (.+)$", info, re.M)
        images.append((label, path, tiles.group(1) if tiles else "?",
                       (fmts.group(1) if fmts else "").split()))
        print(f"image {label}: {images[-1][2]} tile(s), formats {' '.join(images[-1][3])}")

    # Real problems. The two this project has always timed, the chaotic
    # one, and body counts spanning the measured crossover. Step counts
    # fall as N rises because gravity is N^2/2 pairs at eight nodes of
    # every corrector pass; the count is reported per row so no row can
    # quietly mean a different amount of work.
    #
    # THE COUNTS ARE SIZED FROM THE SLOWEST ROW, not the fastest. Every
    # row of a problem runs the same number of steps, and the slowest is
    # the card at binary128 - so a count chosen against REBOUND-in-double,
    # which is three orders faster, produces a run nobody waits for. These
    # are calibrated from the per-step costs measured on 2026-09-14
    # (docs/VALIDATION.md entry 36) to put the slowest row of each problem
    # in the tens of seconds. The first draft of this file asked for 2000
    # Kepler steps, which was five and a half hours of card time for one
    # table, and was killed seventeen minutes in.
    probs = [("kepler", REPO / "data/problems/kepler.txt", 200),
             ("outer", REPO / "data/problems/outer.txt", 100),
             ("pythagorean", REPO / "data/problems/pythagorean.txt", 100)]
    if not a.quick:
        gen = pathlib.Path(a.csv).parent
        gen.mkdir(parents=True, exist_ok=True)
        for n, steps in ((64, 20), (256, 5), (512, 3)):
            p = gen / f"nbody{n}.txt"
            if not p.exists():
                subprocess.run([sys.executable, str(REPO / "tools/make_nbody.py"),
                                str(n), str(p)], capture_output=True)
            if p.exists():
                probs.append((f"nbody{n}", p, steps))

    formats = a.formats.split(",")
    rows, failures = [], 0
    outdir = pathlib.Path(a.csv).parent
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"\n{'problem':<12} {'steps':>6} {'mode':<14} {'format':<6} "
          f"{'seconds':>9} {'steps/s':>9} {'vs ref':>8} {'energy drift':>13} "
          f"{'record':<10} passes")
    for pname, ppath, steps in probs:
        if not pathlib.Path(ppath).exists():
            print(f"{pname}: no such problem file, skipped by name")
            continue
        base = ["--problem", str(ppath), "--steps", str(steps), "--quiet"]

        # The traditional option: REBOUND, hardware double.
        secs, rc, out = run_timed([str(ref), "--problem", str(ppath),
                                   "--steps", str(steps)])
        if rc != 0:
            print(f"{pname}: ias15_ref rc={rc}: {out[-160:]}")
            failures += 1
            continue
        ref_secs, ref_rec = secs, out
        ref_drift = energy_drift(out)
        rows.append((pname, steps, "ref", "fp64", secs, 1.0, ref_drift, "reference", "-"))
        print(f"{pname:<12} {steps:>6} {'ref (REBOUND)':<14} {'fp64':<6} {secs:>9.2f} "
              f"{steps/secs:>9.1f} {1.0:>8.2f} {fmt_drift(ref_drift):>13} "
              f"{'reference':<10} -")

        # The port at binary64 in REBOUND'S OWN rounding sequence, which
        # is where this repository's equivalence claim lives: gate 1 says
        # these values are REBOUND's bit for bit. It is a separate row
        # from the one below on purpose - `--arith fma --engine program`
        # is a DIFFERENT sequence of roundings, documented to differ from
        # REBOUND at round-off (3e-15 at binary64), and holding it to
        # bit-identity with REBOUND reports a failure that is not one.
        # The first run of this table did exactly that.
        secs, rc, out = run_timed([str(cft), "--format", "fp64"] + base)
        if rc == 0:
            same = data_only(out) == data_only(ref_rec)
            rec = "identical" if same else "DIFFERS from REBOUND"
            if not same:
                failures += 1
            d = energy_drift(out)
            rows.append((pname, steps, "sw-rebound-form", "fp64", secs,
                         ref_secs / secs, d, rec, passes(out)))
            print(f"{pname:<12} {steps:>6} {'sw REBOUND fm':<14} {'fp64':<6} {secs:>9.2f} "
                  f"{steps/secs:>9.1f} {ref_secs/secs:>8.2f} {fmt_drift(d):>13} "
                  f"{rec:<10} {passes(out)}")
        else:
            print(f"{pname} sw-rebound-form: rc={rc}: {out[-160:]}")
            failures += 1

        sw_rec = {}
        for f in formats:
            # The port in the FMA form, which is the form the card runs:
            # a correctly rounded divide is not a sequencer program, so
            # --engine program requires --arith fma. This row is the
            # reference the card rows are held to, and it is not held to
            # REBOUND - the row above is.
            secs, rc, out = run_timed([str(cft), "--format", f, "--engine", "program",
                                       "--arith", "fma", "--max-iter", "60"] + base)
            if rc != 0:
                print(f"{pname} sw/{f}: rc={rc}: {out[-160:]}")
                failures += 1
                continue
            sw_rec[f] = out
            d = energy_drift(out)
            rows.append((pname, steps, "sw", f, secs, ref_secs / secs, d,
                         "fma reference", passes(out)))
            print(f"{pname:<12} {steps:>6} {'sw FMA form':<14} {f:<6} {secs:>9.2f} "
                  f"{steps/secs:>9.1f} {ref_secs/secs:>8.2f} {fmt_drift(d):>13} "
                  f"{'fma ref':<10} {passes(out)}")

            for label, path, tiles, have in images:
                if f not in have:
                    print(f"{pname:<12} {steps:>6} {label:<14} {f:<6} "
                          f"{'-':>9} {'-':>9} {'-':>8} {'-':>13} "
                          f"{'no rung':<10} the image does not carry {f}")
                    continue
                secs, rc, out = run_timed([str(cft), "--format", f, "--engine", "program",
                                           "--arith", "fma", "--max-iter", "60",
                                           "--artifact", path] + base)
                if rc != 0:
                    print(f"{pname} {label}/{f}: rc={rc}: {out[-160:]}")
                    failures += 1
                    continue
                d = energy_drift(out)
                same = normalise(out) == normalise(sw_rec.get(f, ""))
                rec = "identical" if same else "DIFFERS from software"
                if not same:
                    failures += 1
                rows.append((pname, steps, label, f, secs, ref_secs / secs, d, rec, passes(out)))
                print(f"{pname:<12} {steps:>6} {label:<14} {f:<6} {secs:>9.2f} "
                      f"{steps/secs:>9.1f} {ref_secs/secs:>8.2f} {fmt_drift(d):>13} "
                      f"{rec:<10} {passes(out)}")

    with open(a.csv, "w", encoding="utf-8", newline="\n") as fh:
        fh.write("problem,steps,mode,format,seconds,steps_per_s,vs_ref,energy_drift,record,passes\n")
        for r in rows:
            fh.write(f"{r[0]},{r[1]},{r[2]},{r[3]},{r[4]:.4f},{r[1]/r[4]:.4f},"
                     f"{r[5]:.4f},{'' if r[6] is None else f'{r[6]:.6e}'},{r[7]},{r[8]}\n")
    print(f"\n{len(rows)} rows -> {a.csv}; {failures} failure(s)")
    print("vs ref = REBOUND-in-double seconds / this row's seconds; above 1.00 beats REBOUND.")
    print("Compare DOWN a column of equal energy drift: a slower row that is many orders")
    print("more accurate is not a loss, and binary64 cannot buy that accuracy at any speed.")
    return 1 if failures else 0


def fmt_drift(d):
    return "-" if d is None else f"{d:.3e}"


if __name__ == "__main__":
    sys.exit(main())
