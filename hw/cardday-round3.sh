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
#
# SESSION B: the parcels' routes, on REV8A, the image that publishes
# their bits. Each leg runs in a clone of the parcel's verified commit,
# merged with this script's branch, and first checks that the build is
# that parcel's. Name its legs to run them; with no legs named the
# script runs session A, above.
#   p4-program   P4's tools/check_program_engine.py on the card: the
#                step-control program against the loop engine.
#   p3-identity  P3's gate_real at every format on the card: a
#                card-written archive, whose identity names the image.
#   p2-corrector P2's tools/check_corrector.py on the card, its runs one
#                at a time.
#   p2-bench     A bench with the device test forced off, then on. The
#                second is also the merged tree's default route on a
#                card, so it is the round's "after" beside split-rev8a's
#                "before".
#   p2-replay    A schedule recorded on the card and replayed there, and
#                one recorded on the software backend and replayed on
#                the card, each held to the card's recording.
#   p1-gate      P1's tools/check_gravity.py on the card.
#   p1-gravity   p1-gate, then a bench for each of the four combinations
#                of indexed and pair_gather.
#   p5-scalar    A bench with the scalar route forced off, then on.
# A route bench reads its build's software rows from a cache
# (bench-workload.py --software-cache), so only the first bench of a
# build pays for them. BENCH_ARGS adds arguments to every bench, for
# example --quick when the window is short.
#
# SESSION C: the merged tree, 2026-11-06. What session B could not fit,
# and what no parcel's own clone could show.
#   ensemble     tools/check_ensemble.py on the card, gates-rev8a's last
#                gate on its own.
#   dropin-light check_dropin --light: the merged engine's drop-in on a
#                tile, at a quarter of check_dropin's card time.
#   dropin-wide  check_dropin --wide, gates-rev8a's binary128 gate on its
#                own.
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
SESSION_A=(probe pin-f128s split-f128s door-rev8a door-f128s split-rev8a tiles gates-rev8a mixed deep)
LEGS=("${SESSION_A[@]}" p4-program p3-identity p2-corrector p2-bench p2-replay p1-gate p1-gravity p5-scalar
      ensemble dropin-light dropin-wide)
# Minutes each leg takes on one tile, on the generous side. Where a leg
# has run, this is its time on 2026-11-02 and a margin; those days' legs
# ran up to four times the first estimates.
declare -A EST=([probe]=2 [pin-f128s]=10 [split-f128s]=45 [door-rev8a]=15 [door-f128s]=15
                [split-rev8a]=45 [tiles]=10 [gates-rev8a]=300 [mixed]=60 [deep]=8
                [p4-program]=65 [p3-identity]=35 [p2-corrector]=150 [p2-bench]=55 [p2-replay]=20
                [p1-gate]=70 [p1-gravity]=130 [p5-scalar]=55
                [ensemble]=90 [dropin-light]=80 [dropin-wide]=150)
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
if [ $# -gt 0 ]; then RUN=("$@"); else RUN=("${SESSION_A[@]}"); fi
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
    "$PY" hw/bench-workload.py --image "$label=$x" --formats fp64,fp128 --csv "$d/workload.csv" \
        --software-cache "$HOME/r3-swcache" ${BENCH_ARGS:-} "$@" \
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

# A card run's route lines, joined. The benches drop them, because they
# differ by design between a card and software, so the route benches
# show them here instead.
routes() { grep '^# route ' "$1" | sed -E 's/^# route ([a-z_]+)=([a-z]+) by=([a-z_]+) calls=([0-9]+) fallbacks=([0-9]+)$/\1=\2 by \3 (\4 calls, \5 fallbacks)/' | paste -sd ';' -; }
# Data lines, the comparison tools/check_records.py makes: identical, or
# how many differ. Only for two runs that both exited 0. Two refusals
# print the same message, and session B's dry run measured this function
# calling them identical.
same_data() {
    local n
    if [ "$(grep -c -v '^#' "$1")" -eq 0 ] || [ "$(grep -c -v '^#' "$2")" -eq 0 ]; then echo EMPTY; return; fi
    n=$(diff <(grep -v '^#' "$1") <(grep -v '^#' "$2") | grep -c '^<')
    [ "$n" -eq 0 ] && echo identical || echo "DIFFER:$n"
}
# The trailer's step counters, which a replay must reproduce as well as
# the bits (ROUND3.md, P2's trap).
counters() { grep -o -E '(steps_done|steps_rejected|iterations_max_exceeded|mean_pc_iterations)=[0-9.]+' "$1" | paste -sd ' ' -; }

# One bench per combination of forced routes, each in its own directory.
# A combination is a comma-separated list of NAME=on|off, for
# CFT_REBOUND_ROUTE_<NAME>. A short run first shows which routes ran,
# because a bench's own output drops the route lines.
route_benches() {
    local d=$1 combo sub kv
    shift
    for combo in "$@"; do
        sub=$d/${combo//[,=]/-}
        mkdir -p "$sub"
        (
            for kv in ${combo//,/ }; do export "CFT_REBOUND_ROUTE_${kv%%=*}=${kv#*=}"; done
            build/ias15_cft --format fp64 --problem data/problems/kepler.txt --steps 2 --artifact "$REV8A" \
                --engine program --arith fma --max-iter 60 --quiet > "$sub/routes.txt" 2>&1
            local rc=$? r
            r=$(routes "$sub/routes.txt")
            if [ "$rc" -ne 0 ] || [ -z "$r" ]; then
                log "FAIL $combo: the route check did not run (rc=$rc): $(last "$sub/routes.txt" 1)"
            else
                log "$combo: routes $r"
            fi
            leg_split "$sub" 1t-rev8a "$REV8A"
        )
    done
}

leg_p4_program() {
    [ -f programs/nextstep-fp64.cfta ] || { log "FAIL this build is not P4's: no programs/nextstep-fp64.cfta"; return; }
    CFT_REBOUND_ARTIFACT=$REV8A "$PY" tools/check_program_engine.py --build build > "$1/check_program_engine.txt" 2>&1
    log "check_program_engine on REV8A (rc=$?): $(tally "$1/check_program_engine.txt"); $(grep -h '^RESULT' "$1/check_program_engine.txt" | tail -n 1)"
}

leg_p3_identity() {
    local fmt
    grep -q CFT_FD_IDENTITY src/cft_ias15_fields.h 2>/dev/null || { log "FAIL this build is not P3's: no identity field"; return; }
    for fmt in fp64 fp128 fp256; do
        CFT_REBOUND_ARTIFACT=$REV8A build/gate_real --$fmt > "$1/gate_real-$fmt.txt" 2>&1
        log "gate_real --$fmt on REV8A (rc=$?): $(tally "$1/gate_real-$fmt.txt"); $(last "$1/gate_real-$fmt.txt" 1)"
    done
    # gate_real holds the words to a 64-digit digest, and prints them only
    # when a run resumes. So one save and resume puts them in the log,
    # beside the image file's own digest.
    CFT_REBOUND_ARTIFACT=$REV8A build/gate_real --fp128 --save "$1/ckpt.bin" > "$1/save.txt" 2>&1
    CFT_REBOUND_ARTIFACT=$REV8A build/gate_real --fp128 --resume "$1/ckpt.bin" > /dev/null 2> "$1/resume.txt"
    log "a card-written archive (rc=$?): $(grep -m1 -o 'written on "[^"]*"' "$1/resume.txt" | cut -c1-220); the image file's sha256 $(sha256sum "$REV8A" | cut -c1-64)"
}

leg_p2_corrector() {
    [ -f tools/check_corrector.py ] || { log "FAIL this build is not P2's: no tools/check_corrector.py"; return; }
    CFT_REBOUND_ARTIFACT=$REV8A "$PY" tools/check_corrector.py --build build > "$1/check_corrector.txt" 2>&1
    log "check_corrector on REV8A (rc=$?): $(tally "$1/check_corrector.txt"); $(grep -h '^RESULT' "$1/check_corrector.txt" | tail -n 1)"
}

leg_p2_bench() {
    [ -f tools/check_corrector.py ] || { log "FAIL this build is not P2's: no tools/check_corrector.py"; return; }
    route_benches "$1" "REDUCE_SEG=off" "REDUCE_SEG=on"
}

leg_p2_replay() {
    local d=$1 spec pname prob fmt steps t0 t1 t2 rc1 rc2 rc3 rc4 c1 c2 c3 v1 v2 verdict
    local -a common
    [ -f tools/check_corrector.py ] || { log "FAIL this build is not P2's: no tools/check_corrector.py"; return; }
    "$PY" tools/make_nbody.py 64 "$d/nbody64.txt" > /dev/null 2>&1
    "$PY" tools/make_nbody.py 256 "$d/nbody256.txt" > /dev/null 2>&1
    for spec in "kepler data/problems/kepler.txt fp128 200" "outer data/problems/outer.txt fp128 100" \
                "nbody64 $d/nbody64.txt fp128 20" "nbody256 $d/nbody256.txt fp128 5"; do
        read -r pname prob fmt steps <<< "$spec"
        common=(--format "$fmt" --problem "$prob" --steps "$steps" --engine program --arith fma --max-iter 60 --quiet)
        t0=$(date +%s.%N)
        build/ias15_cft "${common[@]}" --artifact "$REV8A" --record-iters "$d/$pname.card.iters" > "$d/$pname.card.txt" 2>&1; rc1=$?
        t1=$(date +%s.%N)
        build/ias15_cft "${common[@]}" --artifact "$REV8A" --replay-iters "$d/$pname.card.iters" > "$d/$pname.replay.txt" 2>&1; rc2=$?
        t2=$(date +%s.%N)
        build/ias15_cft "${common[@]}" --record-iters "$d/$pname.sw.iters" > "$d/$pname.sw.txt" 2>&1; rc3=$?
        build/ias15_cft "${common[@]}" --artifact "$REV8A" --replay-iters "$d/$pname.sw.iters" > "$d/$pname.swreplay.txt" 2>&1; rc4=$?
        if [ "$rc1$rc2$rc3$rc4" != 0000 ]; then
            log "FAIL replay $pname $fmt: a run did not complete (the card's recording rc=$rc1, its replay $rc2, the software recording $rc3, its replay on the card $rc4): $(last "$d/$pname.card.txt" 1)"
            continue
        fi
        c1=$(counters "$d/$pname.card.txt"); c2=$(counters "$d/$pname.replay.txt"); c3=$(counters "$d/$pname.swreplay.txt")
        v1=$(same_data "$d/$pname.card.txt" "$d/$pname.replay.txt"); v2=$(same_data "$d/$pname.card.txt" "$d/$pname.swreplay.txt")
        if [ -z "$c1" ] || [ "$c1" != "$c2" ] || [ "$c1" != "$c3" ] || [ "$v1" != identical ] || [ "$v2" != identical ]; then
            verdict=FAIL
        else
            verdict=ok
        fi
        log "$verdict replay $pname $fmt $steps steps: the card's recording $(awk -v a="$t0" -v b="$t1" 'BEGIN { printf "%.2f", b - a }') s, its replay $(awk -v a="$t1" -v b="$t2" 'BEGIN { printf "%.2f", b - a }') s; the replay against the recording: data $v1, counters [$c1] against [$c2]; the software recording replayed on the card: data $v2, counters [$c3]"
    done
}

leg_p1_gate() {
    [ -f tools/check_gravity.py ] || { log "FAIL this build is not P1's: no tools/check_gravity.py"; return; }
    CFT_REBOUND_ARTIFACT=$REV8A "$PY" tools/check_gravity.py --build build > "$1/check_gravity.txt" 2>&1
    log "check_gravity on REV8A (rc=$?): $(tally "$1/check_gravity.txt"); $(grep -h '^check_gravity: routes' "$1/check_gravity.txt" | cut -c1-240)"
}

leg_p1_gravity() {
    leg_p1_gate "$1"
    route_benches "$1" "INDEXED=off,PAIR_GATHER=off" "INDEXED=on,PAIR_GATHER=off" \
                       "INDEXED=off,PAIR_GATHER=on" "INDEXED=on,PAIR_GATHER=on"
}

leg_p5_scalar() {
    [ -f tools/check_scalar.py ] || { log "FAIL this build is not P5's: no tools/check_scalar.py"; return; }
    route_benches "$1" "SCALAR=off" "SCALAR=on"
}

leg_ensemble() {
    CFT_REBOUND_ARTIFACT=$REV8A "$PY" tools/check_ensemble.py --build build > "$1/check_ensemble.txt" 2>&1
    log "check_ensemble (rc=$?): $(tally "$1/check_ensemble.txt"); $(last "$1/check_ensemble.txt" 2)"
}

leg_dropin_light() {
    CFT_REBOUND_ARTIFACT=$REV8A build/check_dropin --light > "$1/check_dropin-light.txt" 2>&1
    log "check_dropin --light (rc=$?): $(tally "$1/check_dropin-light.txt"); $(last "$1/check_dropin-light.txt" 2); $(grep -h -E 'route (indexed|pair_gather|reduce_seg)=' "$1/check_dropin-light.txt" | tail -n 3 | cut -c1-120 | tr '\n' ';')"
}

leg_dropin_wide() {
    CFT_REBOUND_ARTIFACT=$REV8A build/check_dropin --wide > "$1/check_dropin-wide.txt" 2>&1
    log "check_dropin --wide (rc=$?): $(tally "$1/check_dropin-wide.txt"); $(last "$1/check_dropin-wide.txt" 2)"
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
