#!/bin/bash
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Build the specialised tile for this project: cft-fp256's cft_krnl with
# the binary256 rung LEFT OUT (EN_FP256=0 - the "f128" variant, fp32/
# fp64/fp128), N compute units, at a named kernel clock, linked for the
# Alveo U50 with Vitis 2022.2. docs/BITSTREAM.md says why this image
# exists and what it measured; this file is the recipe as an executable,
# because cft-fp256's docs/BITSTREAM-BUILDS.md records what happens to a
# recipe kept in memory.
#
#   hw/build-variant.sh --tiles N --freq HZ --tag NAME [options]
#
#   --tiles N          compute units, 1..8 (hw/gen_layout.py writes the
#                      link config; 8 is the HBM wall)
#   --freq HZ          the kernel clock. REQUIRED, no default: cft-fp256's
#                      rebuild-2022.sh defaults to 10 MHz and an image built
#                      there meets timing trivially and is useless. Below
#                      --min-freq is refused.
#   --tag NAME         names the build dir (build-hw-NAME), the staging
#                      dir (~/cardday-NAME) and the logs (~/f128-logs-NAME)
#   --variant f128     EN_FP256=0 (default). f64: EN_FP128=0 EN_FP256=0.
#   --cft DIR          the cft-fp256 checkout, default third_party/cft-fp256;
#                      its HEAD must be the commit third_party/MANIFEST pins
#                      and its rtl/ and hw/kernel.xml must be clean
#   --build DIR        default <repo>/build-hw-NAME
#   --stage DIR        default ~/cardday-NAME
#   --min-freq HZ      default 100000000
#   --min-free-gb N    refuse to link under this much free memory (20)
#   --no-retiming      kernel synthesis without register retiming
#   --no-phys-opt      no physical optimisation passes
#   --dry-run          every assertion, packaging and the wrapper check,
#                      then stop before v++ (minutes, not hours)
#
# The recipe is cft-fp256's card-day one (135 MHz images: RETIMING=1
# PHYS_OPT=1, default directives) with the clock and the generic changed.
# The .xo is packaged by hw/package_variant.tcl, which sets the generic
# on the packaged IP - cft-fp256's own packager strips every parameter,
# so no narrow tile could be packaged by it. hw/verify_variant_xo.tcl
# then instantiates the packaged IP and reads back what its synthesis
# wrapper passes to cft_krnl, BEFORE the hours are spent: a generic that
# did not survive packaging is the full tile under a variant's name, and
# nothing in a successful link would say so.
#
# What is asserted, and refused rather than warned: the cft-fp256 commit
# (the MANIFEST pin, and rtl/ clean), the clock, free memory, the generic
# in the packaged IP and in the wrapper, and - after the link - that the
# image matches its manifest (cft-fp256's hw/verify-image.sh) before it
# is staged. A negative kernel WNS is printed in capitals and the image
# is still staged, marked: determinism does not depend on the clock, and
# a slower image that is right is a result; it is just not this one.
#
# Run it detached: `setsid nohup hw/build-variant.sh ... > log 2>&1 &`.
set -uo pipefail

TILES=""; FREQ=""; TAG=""; VARIANT="f128"; CFT=""; BUILD=""; STAGE=""
MIN_FREQ=100000000; MIN_FREE_GB=20; RETIMING=1; PHYS_OPT=1; DRY_RUN=0
PLATFORM=${PLATFORM:-xilinx_u50_gen3x16_xdma_5_202210_1}
PART=${PART:-xcu50-fsvh2104-2-e}

die () { echo "FATAL: $*" >&2; exit 1; }

while [ $# -gt 0 ]; do
  case "$1" in
    --tiles)       TILES="${2:-}"; shift 2 ;;
    --freq)        FREQ="${2:-}"; shift 2 ;;
    --tag)         TAG="${2:-}"; shift 2 ;;
    --variant)     VARIANT="${2:-}"; shift 2 ;;
    --cft)         CFT="${2:-}"; shift 2 ;;
    --build)       BUILD="${2:-}"; shift 2 ;;
    --stage)       STAGE="${2:-}"; shift 2 ;;
    --min-freq)    MIN_FREQ="${2:-}"; shift 2 ;;
    --min-free-gb) MIN_FREE_GB="${2:-}"; shift 2 ;;
    --no-retiming) RETIMING=0; shift ;;
    --no-phys-opt) PHYS_OPT=0; shift ;;
    --dry-run)     DRY_RUN=1; shift ;;
    -h|--help)     sed -n '2,60p' "$0"; exit 0 ;;
    *)             die "unknown option $1" ;;
  esac
done

[ -n "$TILES" ] || die "--tiles is required"
[ -n "$FREQ" ]  || die "--freq is required (there is deliberately no default)"
[ -n "$TAG" ]   || die "--tag is required"
case "$TILES" in ''|*[!0-9]*) die "--tiles $TILES is not a number" ;; esac
case "$FREQ"  in ''|*[!0-9]*) die "--freq $FREQ is not a number" ;; esac
[ "$FREQ" -ge "$MIN_FREQ" ] || die "--freq $FREQ is below --min-freq $MIN_FREQ; a slow image meets timing trivially and answers nothing"

case "$VARIANT" in
  f128) GENERICS="EN_FP256=0";            RUNGS="fp32 fp64 fp128" ;;
  f64)  GENERICS="EN_FP128=0 EN_FP256=0"; RUNGS="fp32 fp64" ;;
  *)    die "--variant $VARIANT: f128 or f64" ;;
esac

REPO=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
cd "$REPO" || die "cannot enter $REPO"
[ -f hw/package_variant.tcl ] && [ -f hw/verify_variant_xo.tcl ] && [ -f hw/gen_layout.py ] \
  || die "$REPO does not look like cft-rebound (hw/ incomplete)"
[ -n "$CFT" ]   || CFT="$REPO/third_party/cft-fp256"
[ -n "$BUILD" ] || BUILD="$REPO/build-hw-$TAG"
[ -n "$STAGE" ] || STAGE="$HOME/cardday-$TAG"
LOGDIR="$HOME/f128-logs-$TAG"
CFT=$(cd "$CFT" 2>/dev/null && pwd) || die "no cft-fp256 checkout at $CFT"

# ---- the RTL: which commit, and is it really that commit -------------
PIN=$(awk '$1 == "cft-fp256" {print $3}' third_party/MANIFEST)
[ -n "$PIN" ] || die "third_party/MANIFEST names no cft-fp256 pin"
CFT_HEAD=$(git -C "$CFT" rev-parse HEAD 2>/dev/null) || die "$CFT is not a git checkout"
[ "$CFT_HEAD" = "$PIN" ] || die "cft-fp256 checkout is at $CFT_HEAD, MANIFEST pins $PIN (tools/fetch_third_party.sh)"
# The pointer can be right while the tree is not: assert the content.
CFT_DIRT=$(git -C "$CFT" status --porcelain --untracked-files=all -- rtl/ hw/kernel.xml 2>/dev/null)
[ -z "$CFT_DIRT" ] || die "cft-fp256 rtl/ or hw/kernel.xml differs from $PIN:"$'\n'"$CFT_DIRT"
CFT_RTL_TREE=$(git -C "$CFT" rev-parse "HEAD:rtl")
CFT_XML_BLOB=$(git -C "$CFT" rev-parse "HEAD:hw/kernel.xml")
for f in rtl/cft_krnl.sv hw/kernel.xml hw/verify-image.sh; do
  [ -f "$CFT/$f" ] || die "$CFT lacks $f"
done
# The generic must exist in the RTL this checkout carries; a parameter
# name that does not is packaged as nothing and warns about it.
for g in $GENERICS; do
  grep -q "parameter bit ${g%%=*}" "$CFT/rtl/cft_krnl.sv" || die "rtl/cft_krnl.sv declares no 'parameter bit ${g%%=*}'"
done
REB_HEAD=$(git rev-parse HEAD 2>/dev/null || echo unknown)
REB_DIRT=$(git status --porcelain -- hw/ third_party/MANIFEST 2>/dev/null)

# ---- tools ------------------------------------------------------------
for root in /data/Xilinx /opt/Xilinx /tools/Xilinx; do
  if [ -f "$root/Vitis/2022.2/settings64.sh" ]; then
    # shellcheck disable=SC1091
    set +u; source "$root/Vitis/2022.2/settings64.sh"; set -u
    break
  fi
done
command -v v++ >/dev/null || die "Vitis 2022.2 not found (this platform's shell IP is 2022-era encrypted; docs/BITSTREAM-BUILDS.md)"
command -v vivado >/dev/null || die "vivado not on PATH after sourcing Vitis"
export PLATFORM_REPO_PATHS=${PLATFORM_REPO_PATHS:-$HOME/platforms-root/opt/xilinx/platforms}
[ -d "$PLATFORM_REPO_PATHS/$PLATFORM" ] || die "platform $PLATFORM not under $PLATFORM_REPO_PATHS"
command -v python3 >/dev/null || die "python3 is needed for hw/gen_layout.py"

echo "cft-rebound $VARIANT tile: $TILES x cft_krnl ($RUNGS), $FREQ Hz, retiming=$RETIMING phys_opt=$PHYS_OPT"
echo "cft-fp256 $PIN (rtl tree $CFT_RTL_TREE), cft-rebound $REB_HEAD${REB_DIRT:+ (hw/ or MANIFEST modified)}"
echo "host $(hostname), $(nproc) threads, $(free -g | awk '/^Mem:/ {print $7}') GB free; build $BUILD, stage $STAGE, logs $LOGDIR"

mkdir -p "$BUILD" "$LOGDIR"
S="$LOGDIR/00-summary.txt"
{ echo "cft-rebound $VARIANT: $TILES tiles at $FREQ Hz from cft-fp256 $PIN"; date -Is; } | tee -a "$S"

# ---- link config, and the clock argument derived from it -----------------
LINK_CFG="$BUILD/link.cfg"
python3 hw/gen_layout.py --tiles "$TILES" --out "$LINK_CFG" | tee -a "$S" || die "gen_layout failed"
if [ "$TILES" -eq 1 ] && [ -f "$CFT/hw/link.cfg" ]; then
  python3 hw/gen_layout.py --tiles 1 --check "$CFT/hw/link.cfg" | tee -a "$S" || die "one-tile layout differs from cft-fp256's hw/link.cfg"
fi
if [ "$TILES" -eq 4 ] && [ -f "$CFT/hw/link_quad.cfg" ]; then
  python3 hw/gen_layout.py --tiles 4 --check "$CFT/hw/link_quad.cfg" | tee -a "$S" || die "four-tile layout differs from cft-fp256's hw/link_quad.cfg"
fi
# Exactly rebuild-2022.sh's derivation, for exactly its reason: every CU
# in nk= must be named by --clock.freqHz or it runs at the platform
# default, roughly double this design's ceiling, and nothing says so.
CLOCK_CUS=$(sed -n 's/^[[:space:]]*nk=[^:]*:[0-9]*:\(.*\)$/\1/p' "$LINK_CFG" | tr -d ' ' | head -1)
case "$CLOCK_CUS" in
  "" | *[!A-Za-z0-9_.]* ) die "no usable CU names from the nk= line of $LINK_CFG: $(printf '%q' "$CLOCK_CUS")" ;;
esac
NCU=$(echo "$CLOCK_CUS" | tr '.' '\n' | wc -l)
[ "$NCU" -eq "$TILES" ] || die "nk= names $NCU CUs, --tiles says $TILES"
CLOCK_ARG=$(printf '%s\n' "$CLOCK_CUS" | tr '.' '\n' | sed 's/$/.ap_clk/' | paste -sd,)
echo "Clock constraint targets: $CLOCK_CUS" | tee -a "$S"
echo "Clock constraint argument: ${FREQ}:${CLOCK_ARG}" | tee -a "$S"

# ---- package the variant .xo ------------------------------------------------
rm -f "$BUILD/cft_krnl.xo"
rm -rf "$BUILD/packaged_kernel" "$BUILD/tmp_kernel_pack" "$BUILD/xo_verify"
echo "== package_xo ($(vivado -version | head -1)) generics: $GENERICS" | tee -a "$S"
CFT_GENERICS="$GENERICS" vivado -mode batch -nolog -nojournal \
    -source hw/package_variant.tcl \
    -tclargs "$PART" "$BUILD" "$CFT/rtl" "$CFT/hw/kernel.xml" > "$LOGDIR/package.log" 2>&1
rc=$?
grep -E "^(VARIANT|GENERIC|HDLPARAM):" "$LOGDIR/package.log" | sed 's/^/    /' | tee -a "$S"
[ "$rc" -eq 0 ] && [ -f "$BUILD/cft_krnl.xo" ] || die "packaging failed (rc=$rc): $LOGDIR/package.log"
HDL_PARAMS=$(grep -E "^HDLPARAM:" "$LOGDIR/package.log" | sed 's/^HDLPARAM: //' | tr -d '"' | paste -sd' ')
for g in $GENERICS; do
  n=${g%%=*}; v=${g#*=}
  echo "$HDL_PARAMS" | grep -qE "(^| )$n = $v( |$)" || die "the packaged IP does not carry $n = $v: $HDL_PARAMS"
done

# ---- and what the packaged IP would instantiate ------------------------------
echo "== wrapper check" | tee -a "$S"
vivado -mode batch -nolog -nojournal -source hw/verify_variant_xo.tcl \
    -tclargs "$PART" "$BUILD/packaged_kernel" "$BUILD/xo_verify" > "$LOGDIR/verify_xo.log" 2>&1
rc=$?
grep -E "^WRAPPER" "$LOGDIR/verify_xo.log" | sed 's/^/    /' | tee -a "$S"
[ "$rc" -eq 0 ] || die "wrapper check failed (rc=$rc): $LOGDIR/verify_xo.log"
for g in $GENERICS; do
  n=${g%%=*}; v=${g#*=}
  # .EN_FP256("0") / .EN_FP256(1'b0) / .EN_FP256(1'B0) / .EN_FP256(0):
  # strip the spellings. Vivado 2022.2 writes the sized form with an
  # UPPERCASE B, which the first version of this line did not know, and
  # it refused a correctly packaged .xo - the right outcome for a check
  # whose false negatives cost minutes and whose false positives cost
  # a two-hour link of the wrong tile.
  got=$(grep -E "^WRAPPER_PARAM: *\.$n *\(" "$LOGDIR/verify_xo.log" | head -1 \
        | sed -E "s/^WRAPPER_PARAM: *\.$n *\(([^)]*)\).*/\1/" | tr -d '" ' | sed -E "s/^[0-9]+'[bB]//")
  [ -n "$got" ] || die "the synthesis wrapper names no $n override at all - the generic did not survive packaging"
  [ "$got" = "$v" ] || die "the synthesis wrapper passes .$n($got), wanted $v"
  echo "    wrapper passes .$n($got) - the generic survived packaging" | tee -a "$S"
done

if [ "$DRY_RUN" -eq 1 ]; then
  echo "--dry-run: every assertion passed, the .xo is packaged and its wrapper carries the generic; not linking" | tee -a "$S"
  exit 0
fi

# ---- the link ----------------------------------------------------------------
free_gb=$(free -g | awk '/^Mem:/ {print $7}')
[ "$free_gb" -ge "$MIN_FREE_GB" ] || die "$free_gb GB free, under $MIN_FREE_GB: an OOM-killed place_design reads like a design failure"
vprop=()
if [ "$PHYS_OPT" = 1 ]; then
  vprop+=(--vivado.prop "run.impl_1.{STEPS.PHYS_OPT_DESIGN.IS_ENABLED}=true")
  vprop+=(--vivado.prop "run.impl_1.{STEPS.POST_ROUTE_PHYS_OPT_DESIGN.IS_ENABLED}=true")
fi
if [ "$RETIMING" = 1 ]; then
  # The kernel IP's own synthesis run inside the vpl project; the name
  # follows the kernel name, which is why the kernel is still cft_krnl.
  vprop+=(--vivado.prop "run.ulp_cft_krnl_1_0_synth_1.{STEPS.SYNTH_DESIGN.ARGS.RETIMING}=true")
fi
STARTED=$(date -Is); t0=$(date +%s)
echo "== v++ link -t hw ($TILES CUs, ${FREQ} Hz) started $STARTED" | tee -a "$S"
v++ -l -t hw --platform "$PLATFORM" --config "$LINK_CFG" \
    --clock.freqHz "${FREQ}:${CLOCK_ARG}" "${vprop[@]}" \
    --save-temps --temp_dir "$BUILD/_x_hw" \
    -o "$BUILD/cft_hw.xclbin" "$BUILD/cft_krnl.xo" > "$LOGDIR/link.log" 2>&1
LINK_RC=$?
t1=$(date +%s)
echo "    v++ rc=$LINK_RC in $(( (t1 - t0) / 60 )) min -> $LOGDIR/link.log" | tee -a "$S"

# ---- manifest: rebuild-2022.sh's keys, so hw/verify-image.sh reads it, plus ours
XB="$BUILD/cft_hw.xclbin"
MAN="$BUILD/cft_hw.manifest.txt"
imp="$BUILD/_x_hw/link/vivado/vpl/prj/prj.runs/impl_1"
{
  echo "artifact:      cft_hw.xclbin"
  if [ -f "$XB" ]; then
    echo "sha256:        $(sha256sum "$XB" | cut -d' ' -f1)"
    echo "bytes:         $(stat -c%s "$XB")"
  else
    echo "sha256:        none - the link produced no image"
    echo "bytes:         0"
  fi
  echo "built:         $(date -Is)"
  echo "started:       $STARTED"
  echo "host:          $(hostname)"
  echo "commit:        $PIN"
  echo "describe:      cft-fp256 $(git -C "$CFT" describe --tags --always 2>/dev/null || echo unknown), packaged by cft-rebound $REB_HEAD"
  echo "tree:          rtl/ and hw/kernel.xml clean at $PIN"
  echo "bitstream_sources: rtl/ tree $CFT_RTL_TREE, hw/kernel.xml blob $CFT_XML_BLOB, plus cft-rebound hw/package_variant.tcl"
  echo "variant:       $VARIANT ($RUNGS)"
  echo "generics:      $GENERICS"
  echo "hdl_params:    $HDL_PARAMS"
  echo "tiles:         $TILES"
  echo "cft_fp256_commit:  $PIN"
  echo "cft_rebound_commit: $REB_HEAD${REB_DIRT:+ (hw/ or MANIFEST modified)}"
  echo "link_rc:       $LINK_RC"
  echo "target:        hw"
  echo "platform:      $PLATFORM"
  echo "part:          $PART"
  echo "link_cfg:      $LINK_CFG"
  echo "kernel_freq:   $FREQ"
  echo "retiming:      $RETIMING"
  echo "place_directive: default"
  echo "route_directive: default"
  echo "phys_opt:      $PHYS_OPT"
  echo "vpp_props:     none"
  echo "clock_cus:     $CLOCK_CUS"
  echo "clock_arg:     ${FREQ}:${CLOCK_ARG}"
  echo "vivado:        $(vivado -version 2>/dev/null | head -1)"
  rpt=""
  for cand in "$imp/hw_bb_locked_timing_summary_routed.rpt" "$imp/dr_timing_summary.rpt"; do
    [ -f "$cand" ] && { rpt="$cand"; break; }
  done
  if [ -n "$rpt" ]; then
    echo "timing_report:  $(basename "$rpt")"
    wns=$(awk '/Design Timing Summary/{f=1} f && $1 ~ /^-?[0-9]+\.[0-9]+$/ {print $1; exit}' "$rpt")
    echo "routed_wns_ns: ${wns:-unknown}   # whole design, shell included"
    kclk=${KERNEL_CLK:-clk_out1_ulp_clk_wiz_0}
    kline=$(awk -v c="$kclk" '/Intra Clock Table/{f=1} f && /Inter Clock Table/{exit} f && $1 == c {print $2, $4, $5; exit}' "$rpt")
    if [ -n "$kline" ]; then
      set -- $kline
      echo "kernel_clock:  $kclk"
      echo "kernel_wns_ns: $1   # THIS is the design's margin"
      echo "kernel_failing_endpoints: $2 of $3"
    else
      echo "kernel_clock:  $kclk (not found in Intra Clock Table)"
      echo "kernel_wns_ns: unknown"
    fi
    worst=$(awk '/Intra Clock Table/{f=1} f && $2 ~ /^-[0-9]+\.[0-9]+$/ { printf "%s%s %s %s/%s", sep, $1, $2, $4, $5; sep="; " }' "$rpt")
    [ -n "$worst" ] && printf 'violating_clocks: %s\n' "$worst"
    [ -f "$imp/_new_clk_freq" ] && printf 'final_clocks:  %s\n' "$(tr '\n' ' ' < "$imp/_new_clk_freq")"
  else
    echo "timing_report:  none found under $imp"
  fi
  # Whole-design utilisation after placement: the tile's cost is this
  # minus the shell's, and the shell is the same in every image.
  util=""
  for cand in "$imp/full_util_routed.rpt" "$imp/full_util_placed.rpt"; do
    [ -f "$cand" ] && { util="$cand"; break; }
  done
  if [ -n "$util" ]; then
    echo "util_report:   $(basename "$util")"
    for k in "CLB LUTs" "CLB Registers" "Block RAM Tile" "URAM" "DSPs"; do
      v=$(grep -m1 -E "^\| $k +\|" "$util" | awk -F'|' '{gsub(/ /,"",$3); print $3}')
      printf 'util_%s: %s\n' "$(echo "$k" | tr 'A-Z ' 'a-z_')" "${v:-unknown}"
    done
  fi
} > "$MAN"
echo "== manifest: $MAN" | tee -a "$S"
grep -E "^(sha256|kernel_wns_ns|routed_wns_ns|kernel_failing|violating|final_clocks|util_clb_luts|util_dsps|util_block|util_uram):" "$MAN" | sed 's/^/    /' | tee -a "$S"
WNS_K=$(grep -m1 "kernel_wns_ns:" "$MAN" | awk '{print $2}')
case "$WNS_K" in
  -*) echo "    *** NEGATIVE kernel WNS $WNS_K - this image does NOT close at ${FREQ} Hz ***" | tee -a "$S" ;;
esac

[ "$LINK_RC" -eq 0 ] && [ -f "$XB" ] || { echo "== link failed; nothing staged (manifest kept for the timing it reports)" | tee -a "$S"; exit 1; }

# ---- verify, then stage ------------------------------------------------------
bash "$CFT/hw/verify-image.sh" "$XB" "$MAN" > "$LOGDIR/verify-image.log" 2>&1
vrc=$?
grep -E "^(PASS|FAIL|SKIP|== )" "$LOGDIR/verify-image.log" | sed 's/^/    /' | tee -a "$S"
[ "$vrc" -eq 0 ] || { echo "== verify-image FAILED - not staging" | tee -a "$S"; exit 1; }

NAME="cft_hw_${VARIANT}_${TILES}x"
mkdir -p "$STAGE"
cp "$XB" "$STAGE/$NAME.xclbin"
cp "$MAN" "$STAGE/$NAME.manifest.txt"
cp "$LINK_CFG" "$STAGE/$NAME.link.cfg"
before=$(grep -m1 "^sha256:" "$MAN" | awk '{print $2}')
after=$(sha256sum "$STAGE/$NAME.xclbin" | awk '{print $1}')
[ "$before" = "$after" ] || { echo "*** COPY MISMATCH: manifest $before, copy $after ***" | tee -a "$S"; exit 1; }
( cd "$STAGE" && sha256sum "$NAME.xclbin" "$NAME.manifest.txt" >> SHA256SUMS && sort -u -k2 SHA256SUMS -o SHA256SUMS )
{
  echo "$NAME.xclbin: cft-rebound $VARIANT tile ($RUNGS; binary256 refused by CAPS[3:0] and STATUS[3]), $TILES CU(s), ${FREQ} Hz,"
  echo "  retiming=$RETIMING phys_opt=$PHYS_OPT, kernel WNS ${WNS_K:-?}; cft-fp256 $PIN (rtl tree $CFT_RTL_TREE), cft-rebound $REB_HEAD"
  echo "  sha256 $after; sha256sum -c SHA256SUMS before loading"
} >> "$STAGE/README"
echo "== staged $STAGE/$NAME.xclbin ($after)" | tee -a "$S"
{ echo "=== done ==="; date -Is; } | tee -a "$S"
