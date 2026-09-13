# The f128 image: the tile this project uses, with binary256 left out

A bitstream for the Alveo U50 built from cft-fp256's tile with the
binary256 rung removed (`EN_FP256=0`): fp32, fp64 and fp128, more
compute units than the full tile fits, at a higher clock than the full
tile closes. It is a **performance** product for a workload whose
correctness was settled on the full image: every gate in this repo
passed on the revision-4 pair first (docs/VALIDATION.md), and the
precision study said what a wide format is worth here - binary128
recovers the four orders binary64 loses to round-off and binary256
matches binary128 to every printed digit at every step a user would
take (docs/HORIZON.md). So binary256 is the rung this integrator can do
without, and it is the most expensive one on the tile: 35k of a full
tile's 131k LUTs by cft-fp256's own accounting (its docs/LAYOUTS.md),
and for most of the design's life the path that set its clock.

The image refuses what it lacks, loudly, at three layers, and this
document says how that was proved. It also records what the image
measured, how it was built, how to build it again, and what it asked
of cft-fp256 on the way.

## What it is

One RTL, one generic. cft-fp256's `rtl/cft_krnl.sv` has carried
`EN_FP64`, `EN_FP128` and `EN_FP256` since 2026-08-29; its `PREC_CAPS`
localparam is `{EN_FP256, EN_FP128, EN_FP64, 1'b1}`, published as
`CAPS[3:0]`, and `prec_ok` refuses a run whose `MODE` names a rung the
build lacks - `STATUS[3]`, nothing started, no memory touched. fp32 is
the baseline and has no generic; an fp64/fp128-only tile would need
one (see the asks below). So "fp64/128" here means **fp32 + fp64 +
fp128**, `CAPS[3:0] = 0x7`, contract `0x800`, every sequencer feature
of revision 4 intact: the sequencer, engine, FIFOs, reduction
accumulator and CSR are `BEAT_BITS`-wide and do not shrink with the
rungs.

The kernel keeps the name `cft_krnl`. libcft's XRT backend opens
compute units by the literal name `cft_krnl:{cft_krnl_N}`
(cft-fp256 `host/src/backend_xrt.cpp`), so a kernel renamed after its
variant would open zero tiles. What the image lacks is not in its
name; it is in the tile's own word, `CAPS`, which is what every layer
below reads.

## What it refuses, and where

Three layers, each proven on the card (the "Measured" section carries
the transcripts):

1. **The tile.** `MODE` precision code 3 on this image ends the run at
   start with `STATUS[3]` set and `D` untouched. cft-fp256's
   `cft-resident` issues `MODE` straight at the kernel with no `CAPS`
   check, so `cft-resident <image> -f fp256` is the way to see the
   hardware refuse for itself.
2. **libcft.** `cft_run`, `cft_reduce` and the program path consult
   `format_mask` (`CAPS[3:0]`) and answer `CFT_ERR_UNSUPPORTED` before
   a byte is issued; cft-fp256's own tools skip the format by name
   (`device-test`: "fp256 not on this device, skipped"; `cft-selftest`:
   one line per skipped set).
3. **This repo.** Both open sites - the standalone `ias15_cft` and
   `ias15_engine_open()` behind the drop-in integrator and the
   subprocess API - ask `cft_get_caps` and `cft_supports` at open and
   refuse with a sentence naming the artifact, the formats it carries,
   its `CAPS` mask, tile count and contract, and the format asked for.
   Before 2026-09-13 the refusal arrived at the first operation as a
   bare status: `fma: operation or format not available on this device
   ()`. `ias15_cft --probe [--format F] [--artifact X]` prints the
   device's own word and never refuses, so a reader can ask an image
   what it would refuse before running anything.

`make check` on this image fails its three binary256 legs, by that
refusal, and that is correct: `make check` never skips (CLAUDE.md), and
an image that cannot do binary256 does not pass a binary256 gate. The
card series in `hw/cardtest-f128.sh` runs the binary64 and binary128
gates and demonstrates the binary256 refusal at every entry point
instead.

## How it is built

cft-fp256's packager (`hw/package_kernel.tcl`) removes every user
parameter from the packaged IP, so a bitstream from it can only carry
the RTL defaults - the full tile; its `docs/LAYOUTS.md` names this as
the first of two changes that would turn its narrow layouts from
placeholders into builds. This repo carries that half of the flow
rather than patching cft-fp256:

| file | what |
|---|---|
| `hw/package_variant.tcl` | packages `cft_krnl` from a cft-fp256 checkout with generics from `CFT_GENERICS`, set on the packaged IP's HDL parameters before the user parameters are removed; prints every parameter the `.xo` carries |
| `hw/verify_variant_xo.tcl` | instantiates the packaged IP in a throwaway project and prints the parameter overrides its synthesis wrapper hands to `cft_krnl` - a generic that did not survive packaging shows here as its default, in minutes |
| `hw/gen_layout.py` | the N-tile link config, derived from the `HBM[4(t-1)+m]` rule; `--check` holds it to the one- and four-tile configs cft-fp256 has linked |
| `hw/build-variant.sh` | the recipe: pin asserted by content, clock required, packaging and the wrapper check, `v++` with retiming and phys_opt (cft-fp256's card-day recipe), a manifest in `rebuild-2022.sh`'s shape so `hw/verify-image.sh` verifies it, then staging with `SHA256SUMS` |
| `hw/fit.py` | how many tiles fit, from a linked single's placed LUTs |
| `hw/cardtest-f128.sh` | the card series (below) |

On the box that holds the card, from a checkout of this repo whose
`third_party/cft-fp256` is at the MANIFEST pin (`tools/fetch_third_party.sh`):

```bash
setsid nohup hw/build-variant.sh --tiles 1 --freq 150000000 --tag f128s > ~/f128s.out 2>&1 &
```

The script refuses rather than warns on: a cft-fp256 checkout not at
the pin or with `rtl/` or `hw/kernel.xml` modified; a clock under
100 MHz (cft-fp256's `rebuild-2022.sh` defaults to 10 MHz, and an image
built there meets timing trivially and answers nothing); under 20 GB
free; a packaged IP or a wrapper not carrying the generic; and, after
the link, an image that does not match its manifest. The manifest
records the cft-fp256 commit and its `rtl/` tree hash, the cft-rebound
commit, the generics and every HDL parameter, the clock, both WNS
figures (kernel and whole-design), and the placed utilisation. Build
trees land in `build-hw-<tag>/` (ignored), images in `~/cardday-<tag>/`.

The first run of the flow refused its own correctly packaged `.xo`:
Vivado 2022.2's wrapper spells the bit `1'B0` and the check knew only
`1'b0`. That is the failure mode wanted from a check whose false
negatives cost minutes and whose false positives cost a link of the
wrong tile.

## Measured

TBD - filled from the single at 150 MHz and the multi-tile image.

## The tile count, derived

TBD - `hw/fit.py` on the linked single.

## What it asked of cft-fp256

Reported to cft-fp256 rather than changed there (its fixes were in
flight while this was built):

1. **Packaging strips generics.** `hw/package_kernel.tcl` removes every
   user parameter, so no narrow variant can be packaged by it; the
   RTL's trim generics have never reached a bitstream. Setting the HDL
   parameter's value before the removal is enough (`hw/package_variant.tcl`
   here, a candidate for a `CFT_GENERICS` hook in `package_kernel.tcl`),
   and a wrapper read-back (`hw/verify_variant_xo.tcl`) is the cheap
   proof that it worked.
2. **The XRT backend opens tiles by name.** `backend_xrt.cpp` enumerates
   `cft_krnl:{cft_krnl_N}` literally, so the variant kernel names in
   `hw/layouts/*.cfg` (`cft_krnl_f128_N`) would open no tiles;
   `docs/LAYOUTS.md`'s "homogeneous: works with the current host" is
   true only if the kernel keeps its name. A name-independent lookup
   (by VLNV, or by reading `MAGIC`) would make the catalogue's names
   usable.
3. **A format refusal carries no sentence.** `cft_run`, `cft_reduce`
   and the program path return `CFT_ERR_UNSUPPORTED` for a format the
   device lacks without `cft_set_error`, so `cft_last_error()` is empty
   or stale at the moment a caller most wants it to say "this device
   carries fp32 fp64 fp128". This repo works around it at open.
4. **No bench simulates a trimmed full-beat tile.** `tb/test_krnl.py`
   issues fp256 operations unconditionally and the one trimmed bench,
   `test_krnl_quarter`, is `BEAT_BITS=64`. An `EN_FP256=0` run at
   `BEAT_BITS=256` - `CAPS = 0x7`, the three rungs bit-exact, precision
   3 refused with `STATUS[3]` - would have proved this image's RTL
   before a link did.
5. **fp32 has no generic.** It is the baseline; an fp64/fp128-only tile
   (this project's actual demand) would save its bank - 26.6k LUT by
   `docs/LAYOUTS.md` - per tile.

Not a defect, but noted: cft-fp256's own gate for `make sim` and the
K325T open-core numbers are unaffected; nothing in cft-fp256 was
modified for this image, and the pinned commit's `rtl/` tree hash is
in every manifest.
