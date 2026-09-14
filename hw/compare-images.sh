#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The question this project builds a bitstream to answer: does the
# specialised image beat the general one on THIS workload, and does it
# still produce the same bits?
#
#   hw/compare-images.sh <label>=<xclbin> [<label>=<xclbin> ...]
#
# Runs the same IAS15 problems on the software backend and on every
# named image, at every format all of them carry, and reports per image:
# wall clock, the ratio against software, whether the record is
# byte-identical to the software backend's, and the corrector's mean
# pass count.
#
# THE RECORD COMPARISON IS THE POINT, not a courtesy. A faster image
# that computes different bits is not a faster image, it is a different
# integrator; and the whole claim of this port is that the answer does
# not depend on which backend produced it. Every row prints both.
#
# THE PASS COUNT IS BESIDE EVERY RATIO for the reason docs/VALIDATION.md
# records: bit-identity proves the two backends agree, never that the
# problem posed was the one meant. An ensemble perturbed into nonsense
# once produced perfectly bit-identical runs at a tenth of the corrector
# work and a speedup that reversed when the setup was fixed. If the pass
# count is not what the format needs - about 2.4 at binary64, 6 at
# binary128 - the timing below is about some other problem.
#
# Formats are chosen per image from what it publishes, so an image
# without binary256 contributes rows at the rungs it has and is not
# silently credited with the one it lacks.
set -u
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO" || exit 1
PY=${PYTHON:-python3}
# shellcheck disable=SC1091
source /opt/xilinx/xrt/setup.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=$XILINX_XRT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
T=${TMPDIR:-/tmp}/cmp-$$
mkdir -p "$T"
STEPS=${STEPS:-100}

LABELS=(); IMAGES=()
for spec in "$@"; do
  [ "${spec%%=*}" != "$spec" ] || { echo "argument '$spec' is not label=path" >&2; exit 2; }
  LABELS+=("${spec%%=*}"); IMAGES+=("${spec#*=}")
done
[ ${#IMAGES[@]} -gt 0 ] || { echo "usage: compare-images.sh <label>=<xclbin> ..." >&2; exit 2; }

echo "== images"
printf 'software backend\n'
for i in "${!IMAGES[@]}"; do
  [ -f "${IMAGES[$i]}" ] || { echo "no such image: ${IMAGES[$i]}" >&2; exit 2; }
  fmts=$(./build/ias15_cft --probe --format fp64 --artifact "${IMAGES[$i]}" 2>/dev/null \
         | sed -n 's/^formats: //p')
  tiles=$(./build/ias15_cft --probe --format fp64 --artifact "${IMAGES[$i]}" 2>/dev/null \
          | sed -n 's/^tiles: //p')
  echo "${LABELS[$i]}: ${tiles:-?} tile(s), formats ${fmts:-unreadable}, $(basename "${IMAGES[$i]}")"
  eval "FMTS_$i=\"\$fmts\""
done

# Problems: the two docs/HARDWARE.md timed, plus body counts spanning the
# measured crossover (the card beats a core between 8 and 32 bodies and
# saturates by 64-512), because tile COUNT can only pay where there are
# coordinates to spread over them.
#
# The step count is per problem and falls as the body count rises,
# because gravity is N^2/2 pairs at each of eight nodes of every
# corrector pass: the software backend does about 120 Kepler steps a
# second at binary64 and a great deal fewer at N=512, and a comparison
# nobody waits for is a comparison nobody makes. Ratios are unaffected -
# both sides run the same count - as long as the count is above the
# per-call overhead, which the seconds column shows.
PROBS=(data/problems/kepler.txt data/problems/outer.txt)
STEPS_OF=("$STEPS" "$STEPS")
for n in 64 256 512; do
  p="$T/n$n.txt"
  if $PY tools/make_nbody.py "$n" "$p" > /dev/null 2>&1; then
    PROBS+=("$p")
    case $n in 64) STEPS_OF+=(20) ;; 256) STEPS_OF+=(5) ;; *) STEPS_OF+=(3) ;; esac
  fi
done

# A record with the two fields that cannot match removed: which backend
# produced it, and how long it took. Everything else is kept - the data
# lines and the physics counters - so a real divergence still shows.
# Comparing the raw files instead reports "differ" on every row while
# every value is byte-identical, which is a comparison that can only
# ever fail (measured, 2026-09-13).
norm_record() { sed -E 's#backend=[^ ]+#backend=X#; s/ seconds=[0-9.]+//; s/ steps_per_s=[0-9.]+//; s/ (wall_seconds|t_program|t_elem|t_divsqrt|t_gravity|system_steps_per_s)=[0-9.]+//g' "$1"; }

run_one() {  # <out> <artifact|""> <fmt> <problem> <steps>
  local out=$1 art=$2 fmt=$3 prob=$4 nsteps=$5 t0 t1
  local -a extra=()
  [ -n "$art" ] && extra=(--artifact "$art")
  t0=$(date +%s.%N)
  ./build/ias15_cft --format "$fmt" --problem "$prob" --steps "$nsteps" \
      --engine program --arith fma --max-iter 60 --quiet \
      "${extra[@]}" > "$out" 2>&1
  local rc=$?
  t1=$(date +%s.%N)
  awk -v a="$t0" -v b="$t1" -v r="$rc" 'BEGIN { printf "%.2f %d", b - a, r }'
}

printf '\n%-22s %-6s %-10s %9s %7s %-10s %s\n' \
       problem/steps format image seconds vs-sw record passes
for pi in "${!PROBS[@]}"; do
  prob=${PROBS[$pi]}; nsteps=${STEPS_OF[$pi]}
  pname="$(basename "$prob" .txt)/$nsteps"
  for fmt in fp64 fp128 fp256; do
    read -r sw_s sw_rc <<< "$(run_one "$T/sw.txt" "" "$fmt" "$prob" "$nsteps")"
    [ "$sw_rc" -eq 0 ] || { printf '%-22s %-6s %-10s %9s %7s %-10s %s\n' \
        "$pname" "$fmt" software "$sw_s" - "rc=$sw_rc" "$(tail -c 120 "$T/sw.txt" | tr '\n' ' ')"; continue; }
    passes=$(grep -oE 'mean_pc_iterations=[0-9.]+' "$T/sw.txt" | head -1 | cut -d= -f2)
    printf '%-22s %-6s %-10s %9s %7s %-10s %s\n' "$pname" "$fmt" software "$sw_s" 1.00 "reference" "${passes:-?}"
    for i in "${!IMAGES[@]}"; do
      eval "have=\$FMTS_$i"
      case " $have " in *" $fmt "*) ;; *)
        printf '%-22s %-6s %-10s %9s %7s %-10s %s\n' "$pname" "$fmt" "${LABELS[$i]}" - - "no rung" "the image does not carry $fmt"
        continue ;;
      esac
      read -r c_s c_rc <<< "$(run_one "$T/card-$i.txt" "${IMAGES[$i]}" "$fmt" "$prob" "$nsteps")"
      if [ "$c_rc" -ne 0 ]; then
        printf '%-22s %-6s %-10s %9s %7s %-10s %s\n' "$pname" "$fmt" "${LABELS[$i]}" "$c_s" - "rc=$c_rc" "$(tail -c 140 "$T/card-$i.txt" | tr '\n' ' ')"
        continue
      fi
      ndiff=$(diff <(norm_record "$T/sw.txt") <(norm_record "$T/card-$i.txt") | grep -c '^<')
      if [ "$ndiff" -eq 0 ]; then local_rec="identical"; else local_rec="DIFFER:$ndiff"; fi
      cpasses=$(grep -oE 'mean_pc_iterations=[0-9.]+' "$T/card-$i.txt" | head -1 | cut -d= -f2)
      printf '%-22s %-6s %-10s %9s %7s %-10s %s\n' "$pname" "$fmt" "${LABELS[$i]}" "$c_s" \
             "$(awk -v a="$sw_s" -v b="$c_s" 'BEGIN { printf "%.2f", (b > 0) ? a / b : 0 }')" \
             "$local_rec" "${cpasses:-?}"
    done
  done
done
echo
echo "record 'identical' = byte-identical to the software backend's record for that"
echo "problem and format; the passes column is the corrector's mean per step, which"
echo "must match what the format needs or the timings are about another problem."
rm -rf "$T"
