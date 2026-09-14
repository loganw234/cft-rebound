#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Where does this hardware start beating what a person would otherwise
# run? One dataset, every operating mode, one element-count ladder.
#
#   hw/bench-modes.sh --image "1t-f128=$HOME/cardday-f128s/cft_hw_f128_1x.xclbin" \
#                     --image "4t-full=$HOME/cardday-rev4/cft_hw_quad.xclbin" \
#                     --image "6t-f128=$HOME/cardday-f128x6b/cft_hw_f128_6x.xclbin"
#
# THE MODES, which is the axis that matters rather than the format ladder:
#
#   cpu-hw     the CPU's own double, with a hardware fma. The traditional
#              option at binary64 and the ceiling nothing in software beats.
#   quadmath   gcc __float128 + fmaq. The traditional option at binary128
#              on x86, and the thing a tile has to beat to be worth owning.
#   mpfr       GNU MPFR at the format's precision, and mpfr+754 emulating
#              the binary format exactly. What people reach for when they
#              want arbitrary precision at all.
#   software   libcft's own softfloat, same bits as the card. Not a rival -
#              it is the reference, and it is on the chart so the cost of
#              the CONTRACT can be told apart from the cost of the format.
#   <img>/host      a tile, operands staged across PCIe on every call.
#                   What a first port gets, and what IAS15 does today.
#   <img>/resident  the same tile with cft_alloc'd buffers filled once.
#                   What the engine can do when the bus is not in the way.
#
# The four peer rows come from cft-fp256's own host/cft-bench-peers,
# which is the validated tool for exactly this; this script drives it
# rather than reimplementing it.
#
# THIS IS THE MECHANISM, NOT THE RESULT. It prices one elementwise
# operation, which is a microbenchmark. What an actual REBOUND
# integration does with these modes is hw/bench-workload.py, and that is
# the one to read first; this exists to explain its shape.
#
# Each image is labelled and reports its own compute-unit count, so a
# one-tile and a six-tile image are different rows rather than different
# runs, and the partitioning cost is visible as a line on the same chart.
#
# cft-fp256's own hw/bench-sweep.sh does the format-ladder version of
# this and is the better tool for that question; it takes exactly one
# --single and one --quad, which is why a six-tile image needs this one.
#
# Formats default to fp64 and fp128: binary64 is what people run, and
# binary128 is the rung this project's precision study says is worth
# having (docs/HORIZON.md). Pass --formats to widen it; an image that
# does not carry a format is skipped BY NAME, never silently.
set -uo pipefail

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO" || exit 1
CFT="$REPO/third_party/cft-fp256"
BENCH="$CFT/host/cft-bench"
PEERS="$CFT/host/cft-bench-peers"
OUT="$REPO/bench-modes"; TSEC=0.15; QUICK=0
FORMATS="fp64 fp128"
LABELS=(); PATHS=()

die () { echo "FATAL: $*" >&2; exit 2; }

while [ $# -gt 0 ]; do
  case "$1" in
    --image)   LABELS+=("${2%%=*}"); PATHS+=("${2#*=}"); shift 2;;
    --formats) FORMATS=${2:?}; shift 2;;
    --out)     OUT=${2:?}; shift 2;;
    --time)    TSEC=${2:?}; shift 2;;
    --quick)   QUICK=1; shift;;
    -h|--help) sed -n '2,40p' "$0"; exit 0;;
    *) die "unknown option $1";;
  esac
done

[ -x "$BENCH" ] || die "no cft-bench; make -C $CFT/host XRT=1 cft-bench"
[ -x "$PEERS" ] || die "no cft-bench-peers; make -C $CFT/host XRT=1 cft-bench-peers (needs libmpfr-dev + libgmp-dev)"
for p in ${PATHS+"${PATHS[@]}"}; do [ -f "$p" ] || die "no such image: $p"; done
# shellcheck disable=SC1091
source /opt/xilinx/xrt/setup.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=${XILINX_XRT:-/opt/xilinx/xrt}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}

mkdir -p "$OUT"
CSV="$OUT/modes.csv"; LOG="$OUT/modes.log"; : > "$LOG"
echo "mode,label,tiles,path,format,op,n,ns_per_elem,elems_per_s,reps,seconds" > "$CSV"

# What each image says it is, asked once rather than assumed.
TILES=(); FMTS=()
for i in ${PATHS+"${!PATHS[@]}"}; do
  info=$("$REPO/build/ias15_cft" --probe --format fp64 --artifact "${PATHS[$i]}" 2>/dev/null)
  TILES+=("$(printf '%s' "$info" | sed -n 's/^tiles: //p')")
  FMTS+=("$(printf '%s' "$info" | sed -n 's/^formats: //p')")
  echo "image ${LABELS[$i]}: ${TILES[$i]} tile(s), formats ${FMTS[$i]}" | tee -a "$LOG"
done

ladder_for () {   # <bytes per element>
  local esz=$1 cap n out=""
  cap=$(( 134217728 / esz ))
  for n in 1 4 16 64 256 1024 4096 16384 65536 262144 1048576 4194304 16777216; do
    [ "$n" -le "$cap" ] || break
    if [ "$QUICK" = 1 ]; then
      case "$n" in 1|64|4096|262144|4194304) ;; *) continue;; esac
    fi
    out="$out $n"
  done
  echo "$out"
}
esz_of () { case "$1" in fp32) echo 4;; fp64) echo 8;; fp128) echo 16;; fp256) echo 32;; *) echo 0;; esac; }

rows=0; failed=0
# cft-bench row -> our CSV. Its columns are
# format,bytes_per_elem,op,ns_per_elem,elems_per_s,mb_per_s,reps,seconds
# and only the fma op is kept: one op per chart, and fma is the one the
# integrator actually issues by the thousand.
bench_point () {  # <mode> <label> <tiles> <path> <format> <n> [artifact] [--resident]
  local mode=$1 label=$2 tiles=$3 pth=$4 fmt=$5 n=$6; shift 6
  local out rc
  out=$("$BENCH" "$@" -f "$fmt" -n "$n" -t "$TSEC" --csv 2>>"$LOG"); rc=$?
  if [ "$rc" -ne 0 ] || [ -z "$out" ]; then
    failed=$((failed+1)); echo "  FAILED rc=$rc: $mode/$label $fmt n=$n" >> "$LOG"; return 1
  fi
  printf '%s\n' "$out" | tail -n +2 | awk -F, -v m="$mode" -v l="$label" -v t="$tiles" \
      -v p="$pth" -v n="$n" '$3=="fma" && NF>=8 {
        printf "%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s\n", m,l,t,p,$1,$3,n,$4,$5,$7,$8 }' >> "$CSV"
  rows=$((rows+1)); return 0
}

peer_point () {   # <format> <n>
  local fmt=$1 n=$2 out rc
  out=$("$PEERS" -f "$fmt" -n "$n" -t "$TSEC" --csv 2>>"$LOG"); rc=$?
  [ "$rc" -eq 0 ] && [ -n "$out" ] || { failed=$((failed+1)); return 1; }
  # cft-bench-peers --csv is format,impl,op,ns_per_elem,elems_per_s,rel_time,...
  # Our own is mode,label,tiles,path,format,op,n,ns_per_elem,elems_per_s,reps,seconds.
  # fma only, one op per chart; and libcft's own row is dropped because
  # this script measures that itself through cft-bench, and one number
  # reported twice under two names is how a chart grows a rival that
  # does not exist.
  printf '%s\n' "$out" | tail -n +2 | awk -F, -v n="$n" \
      '$3=="fma" && $2!="libcft" { printf "peer,%s,0,-,%s,%s,%s,%s,%s,,\n", $2,$1,$3,n,$4,$5 }' >> "$CSV"
  rows=$((rows+1)); return 0
}

echo "== modes sweep $(date -Is); host $(uname -n)" | tee -a "$LOG"
for fmt in $FORMATS; do
  esz=$(esz_of "$fmt"); [ "$esz" -gt 0 ] || die "unknown format $fmt"
  for n in $(ladder_for "$esz"); do
    printf '%-6s n=%-9s' "$fmt" "$n" | tee -a "$LOG"
    peer_point "$fmt" "$n"                                       && printf ' peer' | tee -a "$LOG"
    bench_point software software 0 - "$fmt" "$n"                && printf ' sw'   | tee -a "$LOG"
    for i in ${PATHS+"${!PATHS[@]}"}; do
      case " ${FMTS[$i]} " in
        *" $fmt "*) ;;
        *) printf ' %s:no-rung' "${LABELS[$i]}" | tee -a "$LOG"; continue;;
      esac
      bench_point staged   "${LABELS[$i]}" "${TILES[$i]}" "${PATHS[$i]}" "$fmt" "$n" "${PATHS[$i]}" \
        && printf ' %s:st' "${LABELS[$i]}" | tee -a "$LOG"
      bench_point resident "${LABELS[$i]}" "${TILES[$i]}" "${PATHS[$i]}" "$fmt" "$n" "${PATHS[$i]}" --resident \
        && printf ' %s:res' "${LABELS[$i]}" | tee -a "$LOG"
    done
    echo | tee -a "$LOG"
  done
done
echo "== done $(date -Is): $rows rows, $failed failed -> $CSV" | tee -a "$LOG"
[ "$rows" -gt 0 ] || { echo "NO ROWS - the dataset is empty" | tee -a "$LOG"; exit 1; }
