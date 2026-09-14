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

## Measured: the single tile

`~/cardday-f128s/cft_hw_f128_1x.xclbin` on amd-arc-box, one compute
unit at **150 MHz**, built 2026-09-13 from cft-fp256 `ca19fe3` (`rtl/`
tree `0aef8eea`, byte-identical to the revision-4 pair's) by
cft-rebound `247f5c9`, sha256 `90457667...`, 36,646,752 bytes,
`hw/verify-image.sh` 8 of 8 PASS.

**It closes at a clock the full tile has never been asked for.** The
full tile's shipped images are 135 MHz; this one closed 150 with every
endpoint met.

| | full tile (rev-4 single, 135 MHz) | f128 (150 MHz) |
|---|---|---|
| kernel WNS | +0.210 ns | +0.027 ns |
| failing endpoints | 0 of 129,804 | 0 of 105,669 |
| implied path delay | 7.197 ns | 6.640 ns |
| CLB LUTs (placed, shell included) | 252,733 | 220,335 |
| CLB registers | 230,610 | 216,079 |
| DSPs | 311 | 168 |
| BRAM tiles | 258 | 258 |
| URAM | 12 | 12 |

Net of the 123,897-LUT shell the tile is **96,438 LUTs against
128,836**, so the binary256 rung was 32,398 LUTs and 143 DSPs - a
quarter of the tile, close to cft-fp256's model of 35,333 (its
`docs/LAYOUTS.md`). Block RAM and UltraRAM do not move: they belong to
the engine and the sequencer, which are `BEAT_BITS`-wide and do not
shrink with the rungs. `docs/LAYOUTS.md` guessed 150 MHz for this
variant and marked it "target, unmeasured"; it is measured now, and
the guess was right to the megahertz.

**+0.027 ns is not headroom, it is arrival.** Vivado works a path
exactly as hard as the constraint asks and then stops, so a closing
build always lands just above zero and the slack says nothing about
what is left. The number that transfers is the path delay: 6.640 ns
here against the full tile's 7.197 at its own ask. Read 150 MHz as at
or very near this tile's single-CU ceiling, not as a floor.

### What it computes, on the card

Everything this repo and cft-fp256 can ask of the rungs it carries
(`hw/cardtest-f128.sh`, 2026-09-13, log `~/f128-logs-f128s/cardtest.log`):

| gate | result |
|---|---|
| `cft-selftest` over cft-fp256's published vectors | 111 sets, **892,548 cases, all matching**, in 168 s; the binary256 sets skipped by name |
| `device-test`, three modes | 813, 2004 and 687 checks, **0 failed** |
| `gate_real --fp64` | **PASS**, 170 s |
| `gate_real --fp128` | **PASS**, 253 s |
| `check_dropin --light` | **every case passed**, all 14 drop-in refusal rows exercised, 3,062 s |
| `check_program_engine fp128` | **PASS**; the FMA form against REBOUND's own rounding sequence differs by 3.245e-34 at most |

892,548 rather than the full tile's 1,071,635 because one format of
four is gone; nothing else in the census changed.

### The rate

The engine's own, with the operands already resident on the card
(`cft-resident`, n = 1,048,576, 20 reps), against the staged path the
library uses by default (`cft-bench`, same n):

| format | resident M elem/s | Mbeat/s per tile | staged M elem/s | staged GB/s |
|---|---|---|---|---|
| fp32 | 883.9 | 110.5 | 169.9 | 2.72 |
| fp64 | 458.3 | 114.6 | 95.6 | 3.06 |
| fp128 | 234.2 | 117.1 | 46.2 | 2.96 |

Every resident row came back `status 0x0`, with the hardware's bytes
identical to the software backend's, identical across repeats, and
identical across compute units. Power at the wall of the tile was
15.0 W with the card at 37 C.

**The clock reaches the arithmetic.** Measured against the full tile
back to back, same card, same session, same tool, rather than against a
figure quoted from cft-fp256's record:

| format | full tile @135 MHz | f128 @150 MHz | ratio |
|---|---|---|---|
| fp32 | 808.92 M/s (101.11 Mbeat/s) | 897.27 M/s (112.16) | 1.109 |
| fp64 | 417.23 M/s (104.31 Mbeat/s) | 462.03 M/s (115.51) | 1.107 |
| fp128 | 211.96 M/s (105.98 Mbeat/s) | 234.95 M/s (117.47) | 1.108 |

The clock ratio is 150/135 = 1.1111 and the measured gain is 99.7 to
99.8% of it at all three rungs, so the extra megahertz turn into
arithmetic and removing the rung introduced no new wall. The per-format
digests are **identical between the two images**
(`0ab39ba8b1e1ef3a`, `b4dd984b36bc7e7c`, `ca3e223f7f0acad8`): the
determinism claim holds across a change of silicon, not only across
backends. Run-to-run spread is about 1.5%, so the back-to-back pair is
the one the ratio is taken from.

**The staged path does not move and cannot.** All three formats sit
between 2.7 and 3.1 GB/s, the same band every image of this project has
measured, because that path stages operands across PCIe on every call.
A faster tile cannot show there. This matters for reading the wall-clock
table below: `cft_run` is what IAS15 issues.

### What it does for IAS15

Wall clock against the software backend on one core, the program engine
in the FMA form, `--max-iter 60`. The pass count is printed because
bit-identity alone would not prove the problem posed was the one meant.

| problem | steps | format | software | card | ratio | corrector passes |
|---|---|---|---|---|---|---|
| Kepler, N=2 | 200 | fp64 | 2.61 s | 22.02 s | **0.12x** | 3.43 |
| Kepler, N=2 | 200 | fp128 | 6.67 s | 40.25 s | **0.17x** | 8.56 |
| outer, N=6 | 100 | fp128 | 8.67 s | 24.20 s | **0.36x** | 8.42 |
| nbody, N=64 | 20 | fp64 | 24.19 s | 8.77 s | **2.76x** | 3.20 |
| nbody, N=64 | 20 | fp128 | 69.63 s | 22.49 s | **3.10x** | 7.95 |

Every one of those records is **byte-identical to the software
backend's** in every value and every physics counter. The crossover and
the plateau are exactly where docs/HARDWARE.md put them from the full
tile: a card loses badly on a handful of coordinates, where per-call
cost is the whole story, and wins by about 3x from a few dozen bodies
up. **Removing a rung did not move the crossover**, which is the
expected result and worth stating: the crossover is set by per-call
overhead and the scatter in gravity, not by how many rungs the silicon
carries.

So the honest summary of the single tile is that it is 1.10x the full
tile where the tile's own rate is what is being measured, and
indistinguishable from it where IAS15's wall clock is.

**A trap this table set for me.** The first run of it reported "records
differ" on all five rows while every data line was byte-identical: the
comparison was `cmp` on the whole file, and a record carries a header
naming the backend and a trailer carrying elapsed seconds, neither of
which can ever match. Both scripts now normalise exactly those two
fields and compare everything else, counters included.

### What it refuses, on the card

Every layer, measured on this image 2026-09-13 (`hw/cardtest-f128.sh`,
log `~/f128-logs-f128s/cardtest.log`):

- **The tile.** `cft-resident -f fp256` issues `MODE` precision 3 with
  no capability check and the run comes back `status 0x8` - `STATUS[3]`,
  the refusal bit - with the hardware's flags clean and its output
  differing from software, because the tile computed nothing.
- **libcft.** `cft_run(FMA, fp256)` and `cft_reduce(SUM, fp256)` both
  answer status 2, "operation or format not available on this device",
  and the output buffer is **untouched** - not zeroed, which would be
  the worst shape of a wrong answer. `cft_last_error()` is empty, which
  is ask 3 below.
- **cft-fp256's own tools.** `device-test` prints "fp256 not on this
  device, skipped" and passes 813, 2004 and 687 checks with 0 failed in
  its three modes; the counts are lower than a full tile's because one
  format of four is gone.
- **This repo.** The standalone program exits 3 and the drop-in fails
  the simulation, both with: *"carries fp32 fp64 fp128 (CAPS[3:0] =
  0x7, 1 tile, contract 0x00000800) and this run asked for fp256"*.

`ias15_cft --probe` on it reports `formats: fp32 fp64 fp128`,
`format_mask: 0x7`, contract `0x00000800`, and every revision-4
sequencer capacity intact: 64 deposits, 16,384 instructions, 512
constants, 256 scratch slots, `seq_features 0xf1f`, resident buffers.

## The tile count, derived

From the linked single, not from a model of one. `hw/fit.py` takes the
placed LUT total of a one-tile image and the per-compute-unit cost, and
the per-CU cost is itself derived from a linked pair. cft-fp256's
revision-4 pair calibrates it: single 252,733 and quad 640,500 placed,
so three extra compute units cost 387,767, i.e. **419 LUTs of crossbar
per CU beyond the tile itself**. (cft-fp256's `hw/gen_layouts.py` uses
12,626, differenced from 2026-09-02 builds; the revision-4 pair is the
same RTL generation and the same tools as this image, so it is the
better calibration. Both are in `hw/fit.py`.)

    # the per-CU cost, from cft-fp256's OWN linked pair, once:
    python3 hw/fit.py --single-luts 252733 --quad-luts 640500   # -> 419
    # then this tile's table, from its linked single and that cost:
    python3 hw/fit.py --single-luts 220335 --per-cu 419

Passing a `--quad-luts` this script itself predicted would derive the
per-CU cost from the per-CU cost, and the table would then say only
that the arithmetic is self-consistent. `--per-cu` exists so that
cannot happen by accident.

| tiles | HBM PCs | model LUT | of device | verdict |
|---|---|---|---|---|
| 4 | 16 | 510,909 | 58.7% | fits |
| 5 | 20 | 607,767 | 69.8% | fits |
| **6** | **24** | **704,625** | **80.9%** | **tight** |
| 7 | 28 | 801,483 | 92.0% | no |

**Six is the attempt.** It sits at 80.9% of the device, the same
neighbourhood as the full-tile quad that closed at 80.6% on 2026-09-02,
and inside the HBM wall at 24 of 32 pseudo-channels. Seven is 92%,
past cft-fp256's 85% practical routing limit. Against the shipped
four-tile full image that is 50% more compute units.

The clock for it is **140 MHz**, not 150. The full tile lost 0.188 ns
of path going from one compute unit to four at the same ask (+0.210 to
+0.022 at 135 MHz); six units at a higher occupancy should cost at
least as much, which puts this tile's six-CU ceiling near 145 MHz.
140 asks for 7.143 ns against a single-CU path of 6.640 and leaves
about 0.2 ns for the crowding. Six tiles at 140 MHz against four at
135 is 1.56x the aggregate tile-clock product.

### Measured: the multi-tile image

**Six tiles do not close at 140 MHz.** 389 minutes, kernel WNS
**-0.558 ns**, 3,015 of 632,233 endpoints failing, no image. Details
and the lesson in docs/VALIDATION.md entry 35; the short version is
that the **area** model above was right to within 0.9% (698,188 LUTs
placed against 704,620 predicted, with block RAM, UltraRAM and DSPs
exact) and the **timing** estimate was not.

The estimate scaled the full tile's one-unit-to-four penalty of
0.188 ns, measured at 73.6% occupancy, to six units at 80.19%. The real
penalty was about 1.06 ns. Area is a sum; timing is not. A per-unit
timing penalty measured at one occupancy does not predict a higher one,
because congestion cost is flat while the router has room and steep
once it does not. The measured six-unit ceiling is a 7.701 ns path,
about **130 MHz**.

Relaunched at **125 MHz** (8.000 ns, so 0.3 ns of real margin over a
path the tools have demonstrated) rather than 130 (7.692 ns, nine
picoseconds under it).

**It closed:** kernel WNS **+0.061 ns**, 322 minutes, `verify-image`
8 of 8 including *24 masters, all channels HBM, no channel shared*,
sha256 `51f00fb6...`, staged at
`~/cardday-f128x6b/cft_hw_f128_6x.xclbin`. Placed at 696,806 LUTs
(80.03%).

**The arithmetic scaled exactly as asked:** 4,476 / 2,309 / 1,178 M
elem/s at fp32/fp64/fp128 with all six units engaged, 71.6 to 75.4
GB/s, against the shipped quad's 3,224 / 1,664 / 845. That is 1.387 to
1.393 against 6x125 / 4x135 = 1.389 predicted. Per-unit rate tracks the
clock ratio exactly, and the digests match the single-tile images and
software at every format.

## What it is worth to IAS15, which is 2%

The part that matters, and it is not what this image was built for.
Full tables in docs/VALIDATION.md entry 36.

**More tiles is slower, monotonically, in all ten card rows measured** -
one tile beats four beats six, at every body count and both formats,
with every record byte-identical to software. The six-tile image carries
1.39x the aggregate arithmetic of the shipped quad and integrates about
10% slower than it. The library partitions every `cft_run` across all
tiles a device presents, so one call becomes a launch and a staging
round per tile, and IAS15 issues thousands of small calls per step.
Tile count is a divisor on an already short vector and a multiplier on
a fixed cost.

**And the specialised image itself is worth 1 to 3%.** Both
one-compute-unit images head to head on the same card, which is the only
comparison that separates the clock from the partitioning:

| problem | steps | format | f128 @150 | full tile @135 | ratio |
|---|---|---|---|---|---|
| nbody N=64 | 20 | fp128 | 29.21 s | 29.72 s | 1.02 |
| nbody N=256 | 5 | fp128 | 63.61 s | 64.90 s | 1.02 |
| nbody N=512 | 3 | fp128 | 148.09 s | 151.51 s | 1.02 |

So the 1.24x to 2.67x by which this image's single tile beats the
shipped quad is **tile count, not the specialised bitstream**. Loading
cft-fp256's existing one-tile full image recovers 97 to 99% of it and
keeps binary256.

An 11% faster engine moving the integration by 2% puts the arithmetic at
about a fifth of the wall clock. The rest is per-call staging across
PCIe and the scatter in gravity's accumulate half, and no bitstream
addresses either. **The lever is the library** - state resident on the
card across a step, and a device-side scatter - which docs/HARDWARE.md
already ranks.

### What to load

- **IAS15 on this card: a one-tile image**, either kind. Prefer the
  full tile if binary256 is ever wanted, this one for 2% and a loud
  refusal instead.
- **A multi-tile image for long vectors or independent per-tile work**,
  where the 1.39x scaling above is real and perfect. IAS15 is not that
  workload at any body count reachable here.

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
