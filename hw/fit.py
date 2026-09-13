#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
"""How many tiles of a measured kind fit the U50, from a linked single.

Derived from measurements, never typed: give it the placed CLB LUT total
of a ONE-tile image (its manifest's util_clb_luts, or full_util_placed.rpt)
and it prints, for every tile count the HBM wall allows, the modelled
LUT total and whether cft-fp256's layout rules call it a fit. Give it a
linked four-tile image of the SAME kind too and the per-compute-unit
crossbar cost is derived from the pair instead of taken from
cft-fp256's hw/gen_layouts.py.

    python3 hw/fit.py --single-luts 252733                   # the rev-4 full tile
    python3 hw/fit.py --single-luts N --quad-luts M          # calibrated on a pair

Constants and their provenance (cft-fp256 hw/gen_layouts.py, 2026-09-02):
device 870,720 LUT is the "Available" column of a routed report; the
one-CU shell 123,897 was differenced from routed builds; each further CU
adds (161,775 - 123,897) / 3 of crossbar; four masters a tile on 32 HBM
pseudo-channels walls at eight; 80% is "fits" (the quad closed at 80.6%),
85% the practical routing limit.
"""
import argparse

LUT_DEVICE = 870_720
SHELL_1CU = 123_897
SHELL_PER_EXTRA_CU = (161_775 - 123_897) // 3
MAX_TILES = 32 // 4
FITS_PCT, LIMIT_PCT = 80.0, 85.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--single-luts", type=int, required=True,
                    help="placed CLB LUTs of the one-tile image, shell included")
    ap.add_argument("--quad-luts", type=int,
                    help="placed CLB LUTs of a four-tile image of the same tile; "
                         "derives the per-CU crossbar cost from the pair")
    ap.add_argument("--shell", type=int, default=SHELL_1CU)
    ap.add_argument("--device", type=int, default=LUT_DEVICE)
    a = ap.parse_args()

    tile = a.single_luts - a.shell
    per_cu = SHELL_PER_EXTRA_CU
    src = "cft-fp256 hw/gen_layouts.py"
    if a.quad_luts:
        per_cu = (a.quad_luts - a.shell - 4 * tile) // 3
        src = "derived from the single/quad pair"
    print(f"tile: {tile:,} LUT (single {a.single_luts:,} less shell {a.shell:,}); "
          f"each further CU {per_cu:,} ({src}); device {a.device:,}")
    print(f"{'tiles':>5} {'HBM PCs':>7} {'model LUT':>10} {'of device':>9}  verdict")
    for n in range(1, MAX_TILES + 1):
        total = a.shell + (n - 1) * per_cu + n * tile
        pct = 100.0 * total / a.device
        verdict = "fits" if pct <= FITS_PCT else "tight" if pct <= LIMIT_PCT else "no"
        print(f"{n:>5} {4 * n:>7} {total:>10,} {pct:>8.1f}%  {verdict}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
