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
# THE WINDOW. UNTIL=HH:MM is when the card must be free again. A leg is
# skipped, by name, when its estimate in EST below would end past that
# time. A later, cheaper leg can still run. The guard and STOP act only
# between legs. A card process killed mid-run can leave the tile
# unusable until its image is reloaded (cft-fp256's scheduler notes,
# decision N-5). So a leg is never started that would need killing. The
# estimates are believed, from entries 27 and 37 and the dry runs, and
# each leg logs its actual time beside its estimate, which is what
# corrects them.
#
# THE GUARD. The card is shared: cft-fp256's scheduler and nightly runs
# use it too. Before every leg this script looks for another card user.
# - It uses the pinned cft-fp256's own definition (hw/card_stress.py,
#   card_user). It widens it to this port's own programs and scripts
#   (ias15_cft, the gates, the check and bench scripts), counted when
#   they have an .xclbin argument or CFT_REBOUND_ARTIFACT in their
#   environment. As card_stress has it, a cp or an editor naming an
#   image does not hold the card, and neither does an idle shell with
#   the variable exported.
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
#          is the "before", comparable with entries 34 to 37. Its probe
#          read seq_features 0xf1f, contract 0x800 (docs/BITSTREAM.md,
#          "The tile count, derived"). That is the scalar bit, but not
#          indexed or reduce_seg, which came after it. So P1's and P2's
#          routes cannot run on it.
#   REV8A  revision 8's single tile, every rung. cft-fp256 read its
#          seq_features on this card as 0x7ff1f, which holds all three
#          door bits: scalar 0x800, reduce_seg 0x1000, indexed 0x2000.
#   REV8Q  revision 8's quad, opened as one tile (CFT_XRT_TILES=1), for
#          the tile-count row only.
#   DEEP   revision 8's deep single: one tile, the single's seq_features,
#          and 8,192 scratch slots a lane against the single's 4,096
#          (cft-fp256 docs/LANGUAGE.md, scratch-capacity). The port's
#          programs use 10 to 22 slots, so the depth itself should buy
#          nothing here. The deep leg asks whether the image costs or
#          saves per call anyway. The owner pointed it out on 2026-11-02.
# The probe leg prints what each image publishes on this card, at this
# pin. Read it before you trust the two readings above, which were taken
# at older pins.
#
# THE LEGS
#   probe        Each image's own word: VERSION, CAPS and the door's
#                bits.
#   pin-f128s    Does libcft 0.18 open entry 37's image (ROUND3.md, "Read
#                the provider's code", the third claim)? gate_real at
#                binary64 and binary128, and binary256 refused by name.
#                The records, card against software, are split-f128s's.
#   split-f128s  The gravity split's "before" (ROADMAP item 4), by
#                hw/bench-workload.py, whose decompose.csv carries every
#                trailer key. It is what orders P1's routes. Every card
#                row's record is held byte-identical to its software
#                row's, which is the pin's record check too.
#   door-rev8a   tools/check_door.py on a tile that publishes the bits.
#                On software each new route is the same primitive as
#                today's; on this tile it is not, so this is the first
#                test of the tile's routes.
#   door-f128s   The same on a tile with only the scalar bit: the scalar
#                route on the tile, and indexed and reduce_seg reported
#                absent, on today's route only.
#   split-rev8a  The gravity split on the image that the parcels' routes
#                will run on.
#   tiles        CFT_XRT_TILES: REV8Q opened as one tile, through
#                hw/bench-workload.py --quick. That runs the small
#                problems, where tile count cost the most (CLAUDE.md).
#                Compare its rows with split-rev8a's, which ran the same
#                problems and step counts.
#
# The software rows run inside each bench, so they cost card time:
# minutes for the 512-body rows at binary128, measured in the dry run.
# So no leg runs a software reference that another leg already ran. That
# is why the session uses no hw/compare-images.sh, which recomputes its
# own software reference for every problem at every format.
#   gates-rev8a  The gates through the new library at every format:
#                gate_real, check_dropin and check_dropin --wide, and the
#                ensemble gate (docs/HARDWARE.md, "What to do first",
#                item 1).
#   mixed        Idle lanes priced (the same list's item 4), by
#                hw/bench-ensemble.sh on REV8A: first with members that
#                leave the corrector together, then with members that do
#                not.
#   deep         DEEP's probe, then hw/bench-workload.py --quick on it.
#                Its rows compare with split-rev8a's and tiles', which
#                ran the same problems and step counts.
# shellcheck disable=SC1091
set -u
REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO" || exit 1
PY=${PYTHON:-python3}
F128S=${F128S:-$HOME/cardday-f128s/cft_hw_f128_1x.xclbin}
REV8A=${REV8A:-$HOME/cardday-rev8a/cft_hw_single.xclbin}
REV8Q=${REV8Q:-$HOME/cardday-rev8q/cft_hw_quad.xclbin}
DEEP=${DEEP:-$HOME/cardday-rev8d/cft_hw_deep.xclbin}
GUARD_WAIT=${GUARD_WAIT:-600}
UNTIL=${UNTIL:-}
LEGS=(probe pin-f128s split-f128s door-rev8a door-f128s split-rev8a tiles gates-rev8a mixed deep)
# Minutes each leg is believed to take on one tile, on the generous side.
declare -A EST=([probe]=2 [pin-f128s]=5 [split-f128s]=35 [door-rev8a]=25 [door-f128s]=15
                [split-rev8a]=35 [tiles]=10 [gates-rev8a]=110 [mixed]=50 [deep]=8)
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
if [ -n "$UNTIL" ]; then
    until_s=$(date -d "$UNTIL" +%s 2>/dev/null) || die "UNTIL=$UNTIL is not a time of day (HH:MM)"
fi
for b in build/ias15_cft build/ias15_ref build/gate_real build/check_dropin; do
    [ -x "$b" ] || die "$b missing: run --build first"
done
for x in "$F128S" "$REV8A" "$REV8Q"; do [ -f "$x" ] || die "no such image: $x"; done
# XRT first, then the link check. ldd resolves libxrt_coreutil only
# through XRT's own lib directory, so checking before sourcing refused a
# good build. The second dry run measured exactly that.
source /opt/xilinx/xrt/setup.sh >/dev/null 2>&1
export LD_LIBRARY_PATH=${XILINX_XRT:-/opt/xilinx/xrt}/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
xrt_linked > /dev/null || die "build/ias15_cft does not link XRT: build a fresh clone with --build, XRT sourced"
L=$HOME/r3-logs-$TAG
mkdir -p "$L"
LOG=$L/cardday.log

# Every other card user, one line each, then "#ok". Without the "#ok"
# the guard could not look, and the caller must not read that as a free
# card.
foreign() {
    "$PY" - "$REPO/third_party/cft-fp256/hw" <<'EOF'
import os, re, sys
sys.path.insert(0, sys.argv[1])
import card_stress as CS
# This port's programs and scripts. As card_stress has it for its own
# tools, an image named by a cp, a sha256sum or an editor does not hold
# the card, and neither does an idle shell with the variable exported.
PORT = re.compile(r"^(ias15_cft|ias15_ref|gate_[a-z_]+|check_[a-z_]+(\.py)?|bench-[a-z]+\.(py|sh)"
                  r"|cardtest-[a-z0-9]+\.sh|compare-images\.sh)$")
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
    names = [os.path.basename(a) for a in p["args"]]
    ours = any(PORT.match(n[:-4] if n.endswith(".exe") else n) for n in names)
    if not why and ours and any(a.lower().endswith(".xclbin") for a in p["args"]):
        why = "this port's program, with an .xclbin argument"
    if not why and ours and "CFT_REBOUND_ARTIFACT" in CS.proc_environ(p["pid"]):
        why = "this port's program, with CFT_REBOUND_ARTIFACT in its environment"
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
    for fmt in fp64 fp128; do
        CFT_REBOUND_ARTIFACT=$F128S build/gate_real --$fmt > "$1/gate_real-$fmt.txt" 2>&1
        log "gate_real --$fmt on F128S (rc=$?): $(tally "$1/gate_real-$fmt.txt"); $(last "$1/gate_real-$fmt.txt" 2)"
    done
    build/ias15_cft --format fp256 --artifact "$F128S" --problem data/problems/kepler.txt --steps 2 --quiet \
        > "$1/fp256-refused.txt" 2>&1; rc=$?
    # die() exits 3, and so does an image that never opened. So the exit
    # status cannot tell the refusal apart; only engine_check_format()'s
    # own sentence can.
    if [ "$rc" -ne 0 ] && [ "$rc" -lt 128 ] && grep -q 'and this run asked for' "$1/fp256-refused.txt" \
       && grep -q 'nothing was computed' "$1/fp256-refused.txt"; then
        log "binary256 on F128S: refused by name (rc=$rc): $(grep -m1 -o 'carries [^(]*' "$1/fp256-refused.txt" | cut -c1-200)"
    else
        log "FAIL binary256 on F128S was not refused by name (rc=$rc): $(last "$1/fp256-refused.txt" 2)"
    fi
}

leg_split() {  # <dir> <label> <xclbin> [bench-workload option ...]
    local d=$1 label=$2 x=$3
    shift 3
    local rc
    "$PY" hw/bench-workload.py --image "$label=$x" --formats fp64,fp128 --csv "$d/workload.csv" "$@" \
        > "$d/bench-workload.txt" 2>&1; rc=$?
    log "bench-workload $label (rc=$rc): $(card_rows "$d/workload.csv" "$label" 9); $(grep -c 'no rung' "$d/bench-workload.txt") no-rung rows; decompose.csv $(wc -l < "$d/decompose.csv" 2>/dev/null || echo no) lines"
}
# A bench's card rows, counted from its CSV: <csv> <label> <record column>.
# An image that does not open reads to the benches as one that carries
# no format. They print "no rung" rows, write no card row, and exit 0,
# as the third dry run measured. So a card row that never ran is a FAIL
# here, never a silence.
card_rows() {
    local ran ident
    ran=$(awk -F, -v m="$2" '$3 == m' "$1" 2>/dev/null | wc -l)
    ident=$(awk -F, -v m="$2" -v c="$3" '$3 == m && $c == "identical"' "$1" 2>/dev/null | wc -l)
    if [ "$ran" -eq 0 ]; then echo "FAIL no card row ran: the image did not open, or carries no format asked for"
    elif [ "$ident" -ne "$ran" ]; then echo "FAIL $((ran - ident)) of $ran card rows are not identical to software"
    else echo "$ran card rows, every one identical to software"; fi
}
leg_split_f128s() { leg_split "$1" 1t-f128 "$F128S"; }
leg_split_rev8a() { leg_split "$1" 1t-rev8a "$REV8A"; }

leg_door() {  # <dir> <xclbin> <formats>
    CFT_REBOUND_ARTIFACT=$2 "$PY" tools/check_door.py --build build --formats "$3" > "$1/check_door.txt" 2>&1
    log "check_door (rc=$?): $(tally "$1/check_door.txt"); $(grep -E 'identical to the reference|is not published' "$1/check_door.txt" | tr '\n' ';' | cut -c1-700)"
}
leg_door_rev8a() { leg_door "$1" "$REV8A" fp64,fp128,fp256; }
leg_door_f128s() { leg_door "$1" "$F128S" fp64,fp128; }

leg_tiles() { ( export CFT_XRT_TILES=1; leg_split "$1" rev8q-tile1 "$REV8Q" --quick ); }

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

leg_deep() {
    [ -f "$DEEP" ] || { log "FAIL no such image: $DEEP"; return; }
    build/ias15_cft --probe --artifact "$DEEP" > "$1/probe-DEEP.txt" 2>&1
    log "probe DEEP (rc=$?; $(basename "$DEEP"), sha256 $(sha256sum "$DEEP" | cut -c1-16)...): $(tr '\n' ';' < "$1/probe-DEEP.txt" | cut -c1-900)"
    leg_split "$1" 1t-rev8d "$DEEP" --quick
}

leg_mixed() {
    hw/bench-ensemble.sh --image "1t-rev8a=$REV8A" --formats "fp128 fp256" --members "64 512" \
        --out "$1/together" > "$1/together.txt" 2>&1
    log "bench-ensemble, members together (rc=$?): $(card_rows "$1/together/ensemble.csv" 1t-rev8a 8)"
    hw/bench-ensemble.sh --image "1t-rev8a=$REV8A" --formats "fp128 fp256" --members "64 512" \
        --perturb "--vscale planet 0.0005" --out "$1/apart" > "$1/apart.txt" 2>&1
    log "bench-ensemble, members apart (rc=$?): $(card_rows "$1/apart/ensemble.csv" 1t-rev8a 8)"
}

{
    log "=== round 3 card session $TAG: cft-rebound $(git rev-parse --short HEAD), cft-fp256 $(git -C third_party/cft-fp256 rev-parse --short HEAD); legs: ${RUN[*]}; ${UNTIL:+the card free again at $UNTIL; }before: $(therm)"
    ran=0
    for leg in "${RUN[@]}"; do
        if [ -e "$L/STOP" ]; then log "STOPPED by $L/STOP, before $leg"; break; fi
        guard || break
        if [ -n "$UNTIL" ] && [ $(( $(date +%s) + ${EST[$leg]} * 60 )) -gt "$until_s" ]; then
            log "SKIPPED $leg: its estimate, ${EST[$leg]} min, would end past $UNTIL"
            continue
        fi
        d=$L/$leg; mkdir -p "$d"; t0=$(date +%s)
        log "--- $leg (estimate ${EST[$leg]} min)"
        "leg_${leg//-/_}" "$d"
        log "--- $leg done in $(( $(date +%s) - t0 )) s (estimate ${EST[$leg]} min); $(therm)"
        ran=$((ran + 1))
    done
    if [ "$ran" -eq "${#RUN[@]}" ]; then log "CARDDAY-R3-$TAG-COMPLETE"; else log "CARDDAY-R3-$TAG-ENDED: $ran of ${#RUN[@]} legs ran"; fi
} >> "$LOG" 2>&1
