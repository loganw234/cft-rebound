#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Round 3's first card session (ROUND3.md, "The card day"): what P0 alone
# makes measurable, at P0's commit. The parcels' routes are a second
# session, after they merge. Run on the box that holds the card
# (amd-arc-box), in a clone of its own:
#
#   bash hw/cardday-round3.sh --build          # no card: the pins, then make
#   bash hw/cardday-round3.sh --list           # the legs, in order
#   setsid nohup bash hw/cardday-round3.sh <tag> [LEG ...] > /dev/null 2>&1 &
#
# The log is ~/r3-logs-<tag>/cardday.log, and each leg's own output sits
# beside it. The log ends with CARDDAY-R3-<tag>-COMPLETE when every leg
# ran, and with CARDDAY-R3-<tag>-ENDED when the session stopped early.
#
# The legs run in the order below, which is the order of what they
# decide, so a session cut short keeps what mattered most. Naming legs
# runs only those. `touch ~/r3-logs-<tag>/STOP` ends the session before
# its next leg; nothing is killed.
#
# THE GUARD. The card is shared: cft-fp256's scheduler and nightly runs
# use it too. Before every leg this script looks for another card user.
# - It uses the pinned cft-fp256's own definition (hw/card_stress.py,
#   card_user), widened by two of this port's: any process with an
#   .xclbin argument, and any with CFT_REBOUND_ARTIFACT in its
#   environment.
# - Its own process tree is not "another".
# - It waits up to GUARD_WAIT seconds (default 600) for the card, then
#   stops, naming who holds it. GUARD_IGNORE is a regex of command lines
#   it may pass over.
# - It never signals anything.
# - If it cannot look, that is a stop, not a free card.
# The reverse does not hold. cft-fp256's guard does not count this
# port's programs as card users: its XCLBIN_TOOLS names its own tools.
# So a session needs its scheduler to hold the card for the length of
# the window (`cft-sched hold u50-0`).
#
# THE IMAGES are one tile each, because on a card tile count is a cost
# (CLAUDE.md):
#   F128S  entry 37's image: one f128 tile at 150 MHz, no binary256. It
#          is the "before", comparable with entries 34 to 37. Believed to
#          publish none of the door's three bits, since it predates them.
#   REV8A  revision 8's single tile, every rung. cft-fp256 read its
#          seq_features on this card as 0x7ff1f, which holds all three
#          door bits: scalar 0x800, reduce_seg 0x1000, indexed 0x2000.
#   REV8Q  revision 8's quad, opened as one tile (CFT_XRT_TILES=1), for
#          the tile-count row only.
# The probe leg prints what each image publishes. Read it before you
# believe the "believed" and the "read" above.
#
# THE LEGS
#   probe        Each image's own word: VERSION, CAPS and the door's
#                bits.
#   pin-f128s    Does libcft 0.18 open entry 37's image (ROUND3.md, "Read
#                the provider's code", the third claim)? Records card
#                against software at every rung the image carries,
#                gate_real at binary64 and binary128, and binary256
#                refused by name.
#   split-f128s  The gravity split's "before" (ROADMAP item 4), by
#                hw/bench-workload.py, whose decompose.csv carries every
#                trailer key. It is what orders P1's routes.
#   door-rev8a   tools/check_door.py on a tile that publishes the bits.
#                On software each new route is the same primitive as
#                today's; on this tile it is not, so this is the first
#                test of the tile's routes.
#   door-f128s   The same on a tile without the bits: each route
#                reported absent, and today's route only.
#   split-rev8a  The gravity split on the image that the parcels' routes
#                will run on.
#   tiles        CFT_XRT_TILES: REV8Q opened as one tile, against REV8A,
#                each through hw/compare-images.sh.
#   gates-rev8a  The gates through the new library at every format:
#                gate_real, check_dropin and check_dropin --wide, and the
#                ensemble gate (docs/HARDWARE.md, "What to do first",
#                item 1).
#   mixed        Idle lanes priced (the same list's item 4), by
#                hw/bench-ensemble.sh on REV8A: first with members that
#                leave the corrector together, then with members that do
#                not.
# shellcheck disable=SC1091
set -u
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO" || exit 1
PY=${PYTHON:-python3}
F128S=${F128S:-$HOME/cardday-f128s/cft_hw_f128_1x.xclbin}
REV8A=${REV8A:-$HOME/cardday-rev8a/cft_hw_single.xclbin}
REV8Q=${REV8Q:-$HOME/cardday-rev8q/cft_hw_quad.xclbin}
GUARD_WAIT=${GUARD_WAIT:-600}
LEGS=(probe pin-f128s split-f128s door-rev8a door-f128s split-rev8a tiles gates-rev8a mixed)
log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*"; }
die() { echo "FATAL: $*" >&2; exit 2; }
last() { tail -n "${2:-1}" "$1" 2>/dev/null | tr '\n' ' ' | cut -c1-"${3:-300}"; }
therm() { xbutil examine -d 0000:02:00.1 -r thermal -r electrical 2>/dev/null | grep -E "FPGA|Int Vcc|^\s*Power  " | sed 's/^\s*//; s/\s\+/ /g' | tr '\n' ';'; }
# What a reader decides on, counted from the leg's own output, never
# from an exit code alone: the gates' verdict lines and the benches'
# record column.
tally() {
    printf 'ok %s, FAIL %s, identical %s, DIFFER %s' \
        "$(grep -c -E '^\s*ok\b' "$1")" "$(grep -c -E '\bFAIL\b' "$1")" \
        "$(grep -c -E ' identical( |$)' "$1")" "$(grep -c -E '\bDIFFER' "$1")"
}
# The card runtime build/ias15_cft links, printed; nothing, and status 1,
# when it links none.
xrt_linked() { ldd build/ias15_cft 2>/dev/null | grep -o -m1 '/[^ ]*libxrt_coreutil[^ ]*'; }

case "${1:-}" in
  ""|-h|--help) sed -n '2,80p' "$0"; exit 0;;
  --list) printf '%s\n' "${LEGS[@]}"; exit 0;;
  --build)
    bash tools/check_pins.sh || die "the clones are not at the pins (tools/check_pins.sh)"
    [ -z "$(git status --porcelain)" ] || die "the tree is not clean: a session runs a commit, not a working copy"
    source /opt/xilinx/xrt/setup.sh >/dev/null 2>&1
    [ -n "${XILINX_XRT:-}" ] || die "no XRT: /opt/xilinx/xrt/setup.sh set nothing"
    nice -n 10 make all programs build/gate_real || die "make failed"
    # make's "up to date" is not "built with XRT": a tree first built
    # without it stays that way, which the first dry run of this script
    # showed. So the link is checked, not assumed.
    xrt_linked || die "build/ias15_cft does not link XRT, so this tree was first built without it. make clean keeps the pinned libcft.a; build a fresh clone with XRT sourced"
    log "built $(git rev-parse --short HEAD), linked against $(xrt_linked)"
    exit 0;;
  -*) die "unknown option $1";;
esac
TAG=$1; shift
for l in "$@"; do printf '%s\n' "${LEGS[@]}" | grep -qx -- "$l" || die "unknown leg $l (--list names them)"; done
if [ $# -gt 0 ]; then RUN=("$@"); else RUN=("${LEGS[@]}"); fi
for b in build/ias15_cft build/ias15_ref build/gate_real build/check_dropin; do
    [ -x "$b" ] || die "$b missing: run --build first"
done
xrt_linked > /dev/null || die "build/ias15_cft does not link XRT: build a fresh clone with --build, XRT sourced"
for x in "$F128S" "$REV8A" "$REV8Q"; do [ -f "$x" ] || die "no such image: $x"; done
source /opt/xilinx/xrt/setup.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=${XILINX_XRT:-/opt/xilinx/xrt}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
L=$HOME/r3-logs-$TAG
mkdir -p "$L"
LOG=$L/cardday.log

# Every other card user, one line each, then "#ok". Without the "#ok"
# the guard could not look, and the caller must not read that as a free
# card.
foreign() {
    "$PY" - "$REPO/third_party/cft-fp256/hw" <<'EOF'
import os, sys
sys.path.insert(0, sys.argv[1])
import card_stress as CS
procs = CS.list_processes()
skip = CS.descendants_of(os.getpid(), procs)
by_pid = {p["pid"]: p for p in procs}
cur = os.getpid()
while cur in by_pid:                       # our ancestors: this script, its shell, the launcher
    skip.add(cur)
    cur = by_pid[cur]["ppid"]
    if cur in skip:
        break
for p in procs:
    if p["pid"] in skip:
        continue
    why = CS.card_user(p)
    if not why and any(a.lower().endswith(".xclbin") for a in p["args"]):
        why = "an .xclbin argument"
    if not why and "CFT_REBOUND_ARTIFACT" in CS.proc_environ(p["pid"]):
        why = "CFT_REBOUND_ARTIFACT in its environment"
    if why:
        print(p["pid"], why + ":", " ".join(p["args"])[:200])
print("#ok")
EOF
}
guard() {
    local waited=0 out users
    while :; do
        out=$(foreign 2>&1)
        if [ "$(printf '%s\n' "$out" | tail -n 1)" != "#ok" ]; then
            log "STOPPED: the guard could not look for other card users: $(printf '%s' "$out" | tr '\n' ';' | cut -c1-300)"
            return 1
        fi
        users=$(printf '%s\n' "$out" | sed '$d' | { if [ -n "${GUARD_IGNORE:-}" ]; then grep -v -E -- "$GUARD_IGNORE"; else cat; fi; } | paste -sd ';' -)
        [ -z "$users" ] && return 0
        if [ "$waited" -ge "$GUARD_WAIT" ]; then
            log "STOPPED: the card is still in use outside this session after ${waited} s: $users"
            return 1
        fi
        [ "$waited" -eq 0 ] && log "waiting for the card, in use outside this session: $users"
        sleep 30; waited=$((waited + 30))
    done
}

leg_probe() {
    local name x
    for name in F128S REV8A REV8Q; do
        x=${!name}
        build/ias15_cft --probe --artifact "$x" > "$1/probe-$name.txt" 2>&1
        log "probe $name (rc=$?; $(basename "$x"), sha256 $(sha256sum "$x" | cut -c1-16)...): $(tr '\n' ';' < "$1/probe-$name.txt" | cut -c1-900)"
    done
    CFT_XRT_TILES=1 build/ias15_cft --probe --artifact "$REV8Q" > "$1/probe-REV8Q-tile1.txt" 2>&1
    log "probe REV8Q with CFT_XRT_TILES=1 (rc=$?): $(tr '\n' ';' < "$1/probe-REV8Q-tile1.txt" | cut -c1-900)"
}

leg_pin_f128s() {
    local fmt rc
    hw/compare-images.sh "1t-f128=$F128S" > "$1/compare-images.txt" 2>&1
    log "compare-images on F128S (rc=$?): $(tally "$1/compare-images.txt")"
    for fmt in fp64 fp128; do
        CFT_REBOUND_ARTIFACT=$F128S build/gate_real --$fmt > "$1/gate_real-$fmt.txt" 2>&1
        log "gate_real --$fmt on F128S (rc=$?): $(tally "$1/gate_real-$fmt.txt"); $(last "$1/gate_real-$fmt.txt" 2)"
    done
    build/ias15_cft --format fp256 --artifact "$F128S" --problem data/problems/kepler.txt --steps 2 --quiet \
        > "$1/fp256-refused.txt" 2>&1; rc=$?
    if [ "$rc" -ge 1 ] && [ "$rc" -le 2 ] && grep -q -i -E 'binary256|fp256' "$1/fp256-refused.txt"; then
        log "binary256 on F128S: refused by name (rc=$rc): $(grep -m1 -i -E 'binary256|fp256' "$1/fp256-refused.txt" | cut -c1-300)"
    else
        log "FAIL binary256 on F128S was not refused by name (rc=$rc): $(last "$1/fp256-refused.txt" 2)"
    fi
}

leg_split() {  # <dir> <label> <xclbin>
    "$PY" hw/bench-workload.py --image "$2=$3" --formats fp64,fp128 --csv "$1/workload.csv" \
        > "$1/bench-workload.txt" 2>&1
    log "bench-workload $2 (rc=$?): $(tally "$1/bench-workload.txt"); decompose.csv $(wc -l < "$1/decompose.csv" 2>/dev/null || echo no) lines; $(last "$1/bench-workload.txt" 2)"
}
leg_split_f128s() { leg_split "$1" 1t-f128 "$F128S"; }
leg_split_rev8a() { leg_split "$1" 1t-rev8a "$REV8A"; }

leg_door() {  # <dir> <xclbin> <formats>
    CFT_REBOUND_ARTIFACT=$2 "$PY" tools/check_door.py --build build --formats "$3" > "$1/check_door.txt" 2>&1
    log "check_door (rc=$?): $(tally "$1/check_door.txt"); $(grep -E 'identical to the reference|is not published' "$1/check_door.txt" | tr '\n' ';' | cut -c1-700)"
}
leg_door_rev8a() { leg_door "$1" "$REV8A" fp64,fp128,fp256; }
leg_door_f128s() { leg_door "$1" "$F128S" fp64,fp128; }

leg_tiles() {
    hw/compare-images.sh "1t-rev8a=$REV8A" > "$1/rev8a.txt" 2>&1
    log "compare-images REV8A (rc=$?): $(tally "$1/rev8a.txt")"
    CFT_XRT_TILES=1 hw/compare-images.sh "rev8q-tile1=$REV8Q" > "$1/rev8q-tile1.txt" 2>&1
    log "compare-images REV8Q with CFT_XRT_TILES=1 (rc=$?): $(tally "$1/rev8q-tile1.txt")"
}

leg_gates_rev8a() {
    local fmt
    for fmt in fp64 fp128 fp256; do
        CFT_REBOUND_ARTIFACT=$REV8A build/gate_real --$fmt > "$1/gate_real-$fmt.txt" 2>&1
        log "gate_real --$fmt (rc=$?): $(tally "$1/gate_real-$fmt.txt"); $(last "$1/gate_real-$fmt.txt" 2)"
    done
    CFT_REBOUND_ARTIFACT=$REV8A build/check_dropin > "$1/check_dropin.txt" 2>&1
    log "check_dropin (rc=$?): $(tally "$1/check_dropin.txt"); $(last "$1/check_dropin.txt" 2)"
    CFT_REBOUND_ARTIFACT=$REV8A build/check_dropin --wide > "$1/check_dropin-wide.txt" 2>&1
    log "check_dropin --wide (rc=$?): $(tally "$1/check_dropin-wide.txt"); $(last "$1/check_dropin-wide.txt" 2)"
    CFT_REBOUND_ARTIFACT=$REV8A "$PY" tools/check_ensemble.py --build build > "$1/check_ensemble.txt" 2>&1
    log "check_ensemble (rc=$?): $(tally "$1/check_ensemble.txt"); $(last "$1/check_ensemble.txt" 2)"
}

leg_mixed() {
    hw/bench-ensemble.sh --image "1t-rev8a=$REV8A" --formats "fp128 fp256" --members "64 512" \
        --out "$1/together" > "$1/together.txt" 2>&1
    log "bench-ensemble, members together (rc=$?): $(tally "$1/together.txt"); $(last "$1/together.txt" 2)"
    hw/bench-ensemble.sh --image "1t-rev8a=$REV8A" --formats "fp128 fp256" --members "64 512" \
        --perturb "--vscale planet 0.0005" --out "$1/apart" > "$1/apart.txt" 2>&1
    log "bench-ensemble, members apart (rc=$?): $(tally "$1/apart.txt"); $(last "$1/apart.txt" 2)"
}

{
    log "=== round 3 card session $TAG: cft-rebound $(git rev-parse --short HEAD), cft-fp256 $(git -C third_party/cft-fp256 rev-parse --short HEAD); legs: ${RUN[*]}; before: $(therm)"
    ran=0
    for leg in "${RUN[@]}"; do
        if [ -e "$L/STOP" ]; then log "STOPPED by $L/STOP, before $leg"; break; fi
        guard || break
        d=$L/$leg; mkdir -p "$d"; t0=$(date +%s)
        log "--- $leg"
        "leg_${leg//-/_}" "$d"
        log "--- $leg done in $(( $(date +%s) - t0 )) s; $(therm)"
        ran=$((ran + 1))
    done
    if [ "$ran" -eq "${#RUN[@]}" ]; then log "CARDDAY-R3-$TAG-COMPLETE"; else log "CARDDAY-R3-$TAG-ENDED: $ran of ${#RUN[@]} legs ran"; fi
} >> "$LOG" 2>&1
