#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The ensemble across every operating mode: the one workload where a
# multi-tile image gets long vectors, and therefore the fair test the
# six-tile image never had.
#
#   hw/bench-ensemble.sh --image "1t-f128=$HOME/cardday-f128s/cft_hw_f128_1x.xclbin" \
#                        --image "4t-full=$HOME/cardday-rev4/cft_hw_quad.xclbin" \
#                        --image "6t-f128=$HOME/cardday-f128x6b/cft_hw_f128_6x.xclbin"
#
# E independent Kepler systems in one run, members one ulp apart on the
# planet's x (a dyadic --offset, so every member's state is an exact
# binary64 and the perturbation cannot run away - docs/VALIDATION.md's
# ensemble entries record what --geometric did to a run that used it).
# The same integration the ledger's software throughput scan used:
# fixed dt = 0.05, 20 steps, the FMA form on the program engine. At
# E members the predictor and corrector calls carry 6E elements and the
# gravity calls E, so this is where the crossover in elements per call
# is crossed by width rather than by body count.
#
# Reported per (E, format, mode): seconds, system-steps a second (the
# ledger's unit), the ratio against the software backend, whether the
# record is byte-identical to software, and pc_lane_efficiency - the
# fraction of lane-passes that did useful work, which is the physics-
# side counter here: an ensemble whose members leave the corrector at
# different passes wastes lanes, and a ratio quoted without it could be
# a ratio on idle silicon.
set -uo pipefail
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO" || exit 1
PY=${PYTHON:-python3}
OUT="$REPO/bench-modes"; STEPS=20; DT=0.05
FORMATS="fp64 fp128"
MEMBERS="1 8 64 512 4096"
LABELS=(); PATHS=()
die () { echo "FATAL: $*" >&2; exit 2; }
while [ $# -gt 0 ]; do
  case "$1" in
    --image)   LABELS+=("${2%%=*}"); PATHS+=("${2#*=}"); shift 2;;
    --formats) FORMATS=${2:?}; shift 2;;
    --members) MEMBERS=${2:?}; shift 2;;
    --steps)   STEPS=${2:?}; shift 2;;
    --out)     OUT=${2:?}; shift 2;;
    -h|--help) sed -n '2,30p' "$0"; exit 0;;
    *) die "unknown option $1";;
  esac
done
[ -x build/ias15_cft ] || die "build/ias15_cft missing; make first"
for p in ${PATHS+"${PATHS[@]}"}; do [ -f "$p" ] || die "no such image: $p"; done
# shellcheck disable=SC1091
source /opt/xilinx/xrt/setup.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=${XILINX_XRT:-/opt/xilinx/xrt}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
mkdir -p "$OUT"
CSV="$OUT/ensemble.csv"; T=$(mktemp -d)
echo "members,format,mode,tiles,seconds,system_steps_per_s,vs_sw,record,lane_efficiency,passes" > "$CSV"

TILES=(); FMTS=()
for i in ${PATHS+"${!PATHS[@]}"}; do
  info=$(./build/ias15_cft --probe --format fp64 --artifact "${PATHS[$i]}" 2>/dev/null)
  TILES+=("$(printf '%s' "$info" | sed -n 's/^tiles: //p')")
  FMTS+=("$(printf '%s' "$info" | sed -n 's/^formats: //p')")
  echo "image ${LABELS[$i]}: ${TILES[$i]} tile(s), formats ${FMTS[$i]}"
done
norm () { sed -E 's#backend=[^ ]+#backend=X#; s/ seconds=[0-9.]+//; s/ steps_per_s=[0-9.]+//; s/ (wall_seconds|t_program|t_elem|t_divsqrt|t_gravity|system_steps_per_s)=[0-9.]+//g' "$1"; }
eff () { grep -oE 'pc_lane_efficiency=[0-9.]+' "$1" | tail -1 | cut -d= -f2; }
pss () { grep -oE 'mean_pc_iterations=[0-9.]+' "$1" | tail -1 | cut -d= -f2; }
run_one () {  # <out> <fmt> <problem> [--artifact X]
  local out=$1 fmt=$2 prob=$3; shift 3
  local t0 t1
  t0=$(date +%s.%N)
  ./build/ias15_cft --format "$fmt" --problem "$prob" --dt "$DT" --epsilon 0 \
      --steps "$STEPS" --engine program --arith fma --max-iter 60 --quiet "$@" > "$out" 2>&1
  local rc=$?
  t1=$(date +%s.%N)
  awk -v a="$t0" -v b="$t1" -v r="$rc" 'BEGIN { printf "%.2f %d", b - a, r }'
}

printf '\n%-7s %-6s %-10s %5s %9s %11s %7s %-10s %6s %s\n' E format mode tiles seconds sys-steps/s vs-sw record lane-eff passes
for E in $MEMBERS; do
  prob="$T/ens$E.txt"
  $PY tools/make_ensemble.py data/problems/kepler.txt --members "$E" --out "$prob" \
      --offset planet x 0x1p-52 > /dev/null 2>&1 || { echo "E=$E: make_ensemble failed"; continue; }
  for fmt in $FORMATS; do
    read -r sw_s sw_rc <<< "$(run_one "$T/sw.txt" "$fmt" "$prob")"
    if [ "$sw_rc" -ne 0 ]; then
      printf '%-7s %-6s %-10s %5s %9s %11s %7s %-10s %6s %s\n' "$E" "$fmt" software 0 "$sw_s" - - "rc=$sw_rc" - "$(tail -c 120 "$T/sw.txt" | tr '\n' ' ')"
      continue
    fi
    sps=$(awk -v e="$E" -v s="$STEPS" -v t="$sw_s" 'BEGIN { printf "%.1f", (t > 0) ? e * s / t : 0 }')
    printf '%-7s %-6s %-10s %5s %9s %11s %7s %-10s %6s %s\n' "$E" "$fmt" software 0 "$sw_s" "$sps" 1.00 reference "$(eff "$T/sw.txt")" "$(pss "$T/sw.txt")"
    echo "$E,$fmt,software,0,$sw_s,$sps,1.00,reference,$(eff "$T/sw.txt"),$(pss "$T/sw.txt")" >> "$CSV"
    for i in ${PATHS+"${!PATHS[@]}"}; do
      case " ${FMTS[$i]} " in *" $fmt "*) ;; *)
        printf '%-7s %-6s %-10s %5s %9s %11s %7s %-10s %6s %s\n' "$E" "$fmt" "${LABELS[$i]}" "${TILES[$i]}" - - - "no rung" - -
        continue;; esac
      read -r c_s c_rc <<< "$(run_one "$T/card.txt" "$fmt" "$prob" --artifact "${PATHS[$i]}")"
      if [ "$c_rc" -ne 0 ]; then
        printf '%-7s %-6s %-10s %5s %9s %11s %7s %-10s %6s %s\n' "$E" "$fmt" "${LABELS[$i]}" "${TILES[$i]}" "$c_s" - - "rc=$c_rc" - "$(tail -c 120 "$T/card.txt" | tr '\n' ' ')"
        continue
      fi
      n=$(diff <(norm "$T/sw.txt") <(norm "$T/card.txt") | grep -c '^<')
      rec=$([ "$n" -eq 0 ] && echo identical || echo "DIFFER:$n")
      csps=$(awk -v e="$E" -v s="$STEPS" -v t="$c_s" 'BEGIN { printf "%.1f", (t > 0) ? e * s / t : 0 }')
      ratio=$(awk -v a="$sw_s" -v b="$c_s" 'BEGIN { printf "%.2f", (b > 0) ? a / b : 0 }')
      printf '%-7s %-6s %-10s %5s %9s %11s %7s %-10s %6s %s\n' "$E" "$fmt" "${LABELS[$i]}" "${TILES[$i]}" "$c_s" "$csps" "$ratio" "$rec" "$(eff "$T/card.txt")" "$(pss "$T/card.txt")"
      echo "$E,$fmt,${LABELS[$i]},${TILES[$i]},$c_s,$csps,$ratio,$rec,$(eff "$T/card.txt"),$(pss "$T/card.txt")" >> "$CSV"
    done
  done
done
echo
echo "-> $CSV. vs-sw above 1.00 beats the software backend; lane-eff is the share of"
echo "   lane-passes that did useful work, and a ratio without it may be a ratio on idle lanes."
rm -rf "$T"
