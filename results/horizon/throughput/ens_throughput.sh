#!/bin/bash
# Ensemble throughput on the software backend: the same Kepler
# integration the tile ran (fp256, --arith fma --engine program, fixed
# dt = 0.05, 20 steps), at E = 1 .. 1024 members a few ulps apart.
# Sequential, one job at a time, on an otherwise idle box.
R=/c/Users/logan/source/repos/cft-rebound
O=/c/Users/logan/AppData/Local/Temp/claude/C--Users-logan-source-repos/ee78a62e-821f-4dd3-a7cd-6d09b5670eed/scratchpad/rebound2/throughput
PY=${PYTHON:-python3}
mkdir -p "$O"
cd "$R"
for E in 1 4 16 64 256 1024; do
  $PY tools/make_ensemble.py data/problems/kepler.txt --members $E --out "$O/kepler_E$E.txt" --offset planet x 0x1p-52 > /dev/null
  for eng in program loop; do
    if [ "$eng" = loop ] && [ $E -gt 64 ]; then continue; fi
    ./build/ias15_cft.exe --format fp256 --problem "$O/kepler_E$E.txt" --dt 0.05 --epsilon 0 --steps 20 --sample 20 --max-iter 60 --arith fma --engine $eng --programs programs/out --quiet > "$O/E${E}_$eng.rec" 2> "$O/E${E}_$eng.log"
    echo "E=$E $eng: $(grep '^# steps_done' "$O/E${E}_$eng.rec")"
  done
done
# the same at binary64, program engine, for the per-format scaling
for E in 1 64 1024; do
  ./build/ias15_cft.exe --format fp64 --problem "$O/kepler_E$E.txt" --dt 0.05 --epsilon 0 --steps 20 --sample 20 --max-iter 60 --arith fma --engine program --programs programs/out --quiet > "$O/E${E}_program_fp64.rec" 2> /dev/null
  echo "fp64 E=$E program: $(grep '^# steps_done' "$O/E${E}_program_fp64.rec")"
done
echo "throughput done $(date +%T)"
