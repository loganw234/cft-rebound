#!/bin/bash
# 64-member binary64 ensembles through REBOUND itself: the members are
# the Kepler e = 1/2 problem with the planet's x shifted by k * 2^-52
# (2k ulps), so their spread at any time is the distribution of the
# round-off random walk. Adaptive at REBOUND's default epsilon, and
# fixed dt = 1/16; 1e7 steps each (2e5 and 1e5 orbits).
R=/c/Users/logan/source/repos/cft-rebound
O=/c/Users/logan/AppData/Local/Temp/claude/C--Users-logan-source-repos/ee78a62e-821f-4dd3-a7cd-6d09b5670eed/scratchpad/rebound2/ens64
PY=C:/Users/logan/AppData/Local/Programs/Miniconda3/python.exe
mkdir -p "$O/members" "$O/adapt" "$O/fixed"
cd "$R"
$PY tools/make_ensemble.py data/problems/kepler.txt --members 64 --out "$O/kepler_ulps_E64.txt" --offset planet x 0x1p-52 --members-dir "$O/members"
run_adapt(){ k=$1; "$R/build/ias15_ref.exe" --problem "$O/members/m$k.txt" --dt 0.01 --epsilon 1e-9 --steps 10000000 --sample 100000 > "$O/adapt/m$k.rec" 2> /dev/null; $PY "$R/tools/oracle.py" "$O/adapt/m$k.rec" --csv "$O/adapt/m$k.csv" --quiet > "$O/adapt/m$k.summary"; echo "adapt m$k done $(date +%T)"; }
run_fixed(){ k=$1; "$R/build/ias15_ref.exe" --problem "$O/members/m$k.txt" --dt 0.0625 --epsilon 0 --steps 10000000 --sample 100000 > "$O/fixed/m$k.rec" 2> /dev/null; $PY "$R/tools/oracle.py" "$O/fixed/m$k.rec" --csv "$O/fixed/m$k.csv" --quiet > "$O/fixed/m$k.summary"; echo "fixed m$k done $(date +%T)"; }
export -f run_adapt run_fixed
export R O PY
seq 0 63 | xargs -P 4 -I{} bash -c 'run_adapt {}'
seq 0 63 | xargs -P 4 -I{} bash -c 'run_fixed {}'
echo "ens64 done $(date +%T)"
