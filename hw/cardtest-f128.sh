#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The card series for an f128 image, on the box that holds the card
# (amd-arc-box). Detached: `setsid nohup hw/cardtest-f128.sh <xclbin>
# <name> > /dev/null 2>&1 &`; the log is ~/f128-logs-<name>/cardtest.log
# and ends with CARDTEST-<name>-COMPLETE.
#
#   hw/cardtest-f128.sh <xclbin> <name>
#
# Needs: this repo built here with XRT sourced (make all programs
# build/gate_real - the Makefile links the card backend when XILINX_XRT
# is set), the pinned cft-fp256's host tools built XRT=1, and cft-fp256's
# published vectors at $VECTORS (default ~/cft-fp256/vectors/out).
#
# What it proves, in order. (0) What the image says about itself.
# (1) The refusal, which is the point of the image: binary256 refused BY
# NAME at every cft-rebound entry point (the standalone program, the
# drop-in through gate_real) with nothing computed, and by the TILE
# itself when the CAPS check is bypassed (cft-resident issues MODE
# precision 3 straight at the kernel: STATUS[3], nothing written).
# (2) cft-fp256's own matrix and published sets on the rungs the image
# carries, skipping binary256 by name. (3) The engine's rate. (4) This
# project's gates at binary64 and binary128 on the card, and the
# program engine's records card-vs-software. (5) Wall-clock, card
# against the software backend, on the problems docs/HARDWARE.md timed.
set -u
X=${1:?xclbin}; NAME=${2:?name}
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO" || exit 1
H=$REPO/third_party/cft-fp256/host
HW=$REPO/third_party/cft-fp256/hw
VECTORS=${VECTORS:-$HOME/cft-fp256/vectors/out}
PY=${PYTHON:-python3}
# shellcheck disable=SC1091
source /opt/xilinx/xrt/setup.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=/opt/xilinx/xrt/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
L=$HOME/f128-logs-$NAME; T=/tmp/f128-$NAME
mkdir -p "$L" "$T"
LOG=$L/cardtest.log
log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
# A record with the two fields that cannot match removed: which backend
# produced it, and how long it took. Nothing else is touched.
norm_record() { sed -E 's#backend=[^ ]+#backend=X#; s/ seconds=[0-9.]+//; s/ steps_per_s=[0-9.]+//; s/ (wall_seconds|t_program|t_elem|t_divsqrt|t_gravity|system_steps_per_s)=[0-9.]+//g' "$1"; }
therm() { xbutil examine -d 0000:02:00.1 -r thermal -r electrical 2>/dev/null | grep -E "FPGA|Int Vcc|^\s*Power  " | sed 's/^\s*//; s/\s\+/ /g' | tr '\n' ';'; }
last() { tail -n "${2:-1}" "$1" 2>/dev/null | tr '\n' ' ' | cut -c1-"${3:-300}"; }
{
  log "=== f128 card series on $NAME ($X); sha $(sha256sum "$X" | cut -c1-16)...; cft-rebound $(git rev-parse --short HEAD), cft-fp256 $(git -C third_party/cft-fp256 rev-parse --short HEAD); before: $(therm)"
  for b in build/ias15_cft build/gate_real build/check_dropin "$H/device-test" "$H/cft-selftest" "$H/cft-resident" "$H/cft-bench"; do
    [ -x "$b" ] || log "MISSING: $b"
  done
  [ -d "$VECTORS" ] || log "MISSING: vectors at $VECTORS"

  # ---- 0. the image's own word ----------------------------------------------
  build/ias15_cft --probe --format fp256 --artifact "$X" > "$T/probe.txt" 2>&1
  log "probe (rc=$?): $(tr '\n' ';' < "$T/probe.txt" | cut -c1-400)"

  # ---- 1. the refusal, by name, at every entry point --------------------------
  for fmt in fp64 fp128 fp256; do
    build/ias15_cft --format "$fmt" --artifact "$X" --problem data/problems/kepler.txt --steps 2 --quiet > "$T/open-$fmt.txt" 2>&1; rc=$?
    log "ias15_cft --format $fmt --steps 2: rc=$rc; $(grep -m1 -E 'carries|ias15_cft:' "$T/open-$fmt.txt" | cut -c1-360)"
  done
  CFT_REBOUND_ARTIFACT=$X build/gate_real --fp256 > "$T/gate_real-fp256.txt" 2>&1; rc=$?
  log "gate_real --fp256 (drop-in path): rc=$rc; $(grep -m1 -E 'carries|ERROR|ias15_cft' "$T/gate_real-fp256.txt" | cut -c1-360)"
  # libcft's own word, with no help from this repo: the status, the
  # sentence (or its absence) and whether the output was touched.
  cc -std=c99 -O2 -Wall -I"$H/include" -o "$T/refusal_probe" hw/refusal_probe.c "$H/libcft.a" \
     -lm -L"$XILINX_XRT/lib" -lxrt_coreutil -lstdc++ -lpthread -luuid > "$T/refusal_probe.build" 2>&1 \
     || log "refusal_probe did not build: $(last "$T/refusal_probe.build" 3 300)"
  for fmt in fp128 fp256; do
    "$T/refusal_probe" "$X" "$fmt" > "$T/refusal-$fmt.txt" 2>&1; rc=$?
    log "refusal_probe $fmt (libcft directly): rc=$rc; $(tr '\n' ';' < "$T/refusal-$fmt.txt" | cut -c1-420)"
  done
  # The tile's own refusal: cft-resident reads no CAPS and issues MODE
  # precision 3 to the kernel; the run must end with STATUS[3] and
  # nothing written. Its non-zero rc IS the pass here.
  "$H/cft-resident" "$X" -f fp256 --op fma -n 4096 -r 2 > "$T/resident-fp256.txt" 2>&1; rc=$?
  log "cft-resident -f fp256 (MODE prec 3 straight at the tile, no CAPS check): rc=$rc; $(grep -m2 -E 'status|STATUS|refus' "$T/resident-fp256.txt" | tr '\n' ';' | cut -c1-300)"

  # ---- 2. cft-fp256's matrix and published sets on the rungs carried ---------
  for args in "-q -n 8" "-n 4096" "-r"; do
    tag=$(echo "$args" | tr -d ' -'); t0=$(date +%s)
    bash "$HW/run-device-test.sh" "$X" $args > "$T/dt-$tag.log" 2>&1; rc=$?
    log "device-test $args: rc=$rc in $(( $(date +%s)-t0 )) s - $(grep -aE 'checks, [0-9]+ failed$|^device:|not on this device' "$T/dt-$tag.log" | tail -n 3 | tr '\n' ';' | cut -c1-320)"
  done
  t0=$(date +%s)
  "$H/cft-selftest" "$VECTORS" "$X" > "$T/selftest.log" 2>&1; rc=$?
  log "cft-selftest: rc=$rc in $(( $(date +%s)-t0 )) s - $(grep -aE 'cases checked|all matching|MISMATCH|differ' "$T/selftest.log" | tail -n 2 | tr '\n' ';' | cut -c1-200); binary256 sets skipped by name: $(grep -ac 'fp256 not on this device' "$T/selftest.log")"

  # ---- 3. the engine's rate ---------------------------------------------------
  for fmt in fp32 fp64 fp128; do
    "$H/cft-resident" "$X" -f "$fmt" --op fma -n 1048576 -r 20 > "$T/resident-$fmt.txt" 2>&1; rc=$?
    log "cft-resident -f $fmt fma n=1M: rc=$rc; $(grep -aE 'M/s|Mbeat|elem/s|GB/s' "$T/resident-$fmt.txt" | head -n 2 | tr '\n' ';' | cut -c1-260)"
  done
  "$H/cft-bench" "$X" -n 1048576 -t 1 --csv > "$T/bench.csv" 2> "$T/bench.err"; rc=$?
  log "cft-bench n=1M (staged): rc=$rc; fma rows (fmt,ns/elem,elems/s,MB/s): $(awk -F, '$3=="fma"{print $1","$4","$5","$6}' "$T/bench.csv" | tr '\n' ';' | cut -c1-240)"
  log "  thermals: $(therm)"

  # ---- 4. this project's gates on the card --------------------------------------
  for fmt in fp64 fp128; do
    t0=$(date +%s)
    CFT_REBOUND_ARTIFACT=$X build/gate_real --$fmt > "$T/gate_real-$fmt.txt" 2>&1; rc=$?
    log "gate_real --$fmt on the card: rc=$rc in $(( $(date +%s)-t0 )) s; $(last "$T/gate_real-$fmt.txt" 2 240)"
  done
  t0=$(date +%s)
  CFT_REBOUND_ARTIFACT=$X build/check_dropin --light > "$T/check_dropin-light.txt" 2>&1; rc=$?
  log "check_dropin --light on the card: rc=$rc in $(( $(date +%s)-t0 )) s; $(last "$T/check_dropin-light.txt" 3 300)"
  t0=$(date +%s)
  CFT_REBOUND_ARTIFACT=$X $PY tools/check_program_engine.py --build build --formats fp128 > "$T/program_engine-fp128.txt" 2>&1; rc=$?
  log "check_program_engine fp128 (program engine on the card vs the host loop): rc=$rc in $(( $(date +%s)-t0 )) s; $(last "$T/program_engine-fp128.txt" 2 300)"

  # ---- 5. wall clock, card against the software backend ------------------------
  [ -f "$T/n64.txt" ] || $PY tools/make_nbody.py 64 "$T/n64.txt" > /dev/null 2>&1
  time_one() {  # <label> <ias15_cft args...>
    local label=$1; shift
    local t0 t1 sw card
    t0=$(date +%s.%N); build/ias15_cft "$@" --quiet > "$T/time-sw.txt" 2>&1; t1=$(date +%s.%N)
    sw=$(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.2f", b - a }')
    t0=$(date +%s.%N); build/ias15_cft "$@" --artifact "$X" --quiet > "$T/time-card.txt" 2>&1; t1=$(date +%s.%N)
    card=$(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.2f", b - a }')
    # Compare the RECORD, not the file. Two comment lines cannot match
    # and must not be allowed to say the arithmetic did not: the header
    # names the backend, and the trailer carries the elapsed time. A
    # plain cmp reported "records differ" on every row of this table on
    # 2026-09-13 while every data line was byte-identical. Everything
    # else is kept, the physics counters included, so a real divergence
    # in steps_done, flags_seen, calls or the values themselves still
    # shows.
    local same="records DIFFER"
    if diff -q <(norm_record "$T/time-sw.txt") <(norm_record "$T/time-card.txt") > /dev/null; then
      same="records identical"
    else
      same="records DIFFER in $(diff <(norm_record "$T/time-sw.txt") <(norm_record "$T/time-card.txt") | grep -c '^<') lines"
    fi
    # The physics counter beside the timing: identical records prove
    # the backends agree, and the corrector's pass count proves the
    # problem posed was the one meant (bit-identity is not validity).
    log "  $label: software ${sw}s, card ${card}s, ratio $(awk -v a="$sw" -v b="$card" 'BEGIN { printf "%.2f", (b > 0) ? a / b : 0 }'); $same; $(grep -oE 'mean_pc_iterations=[0-9.]+' "$T/time-card.txt" | head -n 1)"
  }
  log "wall clock, --engine program --arith fma (docs/HARDWARE.md's form), card vs software:"
  time_one "kepler N=2 fp64 200 steps"   --format fp64  --problem data/problems/kepler.txt --steps 200 --engine program --arith fma --max-iter 60
  time_one "kepler N=2 fp128 200 steps"  --format fp128 --problem data/problems/kepler.txt --steps 200 --engine program --arith fma --max-iter 60
  time_one "outer N=6 fp128 100 steps"   --format fp128 --problem data/problems/outer.txt  --steps 100 --engine program --arith fma --max-iter 60
  time_one "nbody N=64 fp64 20 steps"    --format fp64  --problem "$T/n64.txt" --steps 20 --engine program --arith fma --max-iter 60
  time_one "nbody N=64 fp128 20 steps"   --format fp128 --problem "$T/n64.txt" --steps 20 --engine program --arith fma --max-iter 60
  log "=== CARDTEST-$NAME-COMPLETE; after: $(therm)"
} >> "$LOG" 2>&1
