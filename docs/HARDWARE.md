# IAS15 as a resident workload on the tile: the design, the numbers, and what is unverified

**Read this before any number below.** This document was written as a
design plus a software-backend prototype, at a time when no card was
opened, no XRT call was made and no `cft://` server was contacted for
it, because both were in use. Every card number in the design sections
is therefore either quoted from cft-fp256's own measurements (its
docs/INTEGRATION.md and docs/BENCHMARKS.md, in the pinned clone under
third_party/, dated 2026-09-09) or is a projection built from those
and from operation counts measured here on the software backend.
Projections are marked as such, each time.

Two sections are **not** projections and were added later, from runs on
an Alveo U50C: "What the card said, and the crossover", and item 0 of
"What to do first". docs/VALIDATION.md entries 27 (the gate suite on
the quad tile) and the "first run on a tile" entry are their record.
Where a projection and a measurement disagree, the measurement is the
one to believe; the projections are kept because the model they came
from is what the measurements are a check on.

## Why IAS15 is unusually favourable, and where it is not

An IAS15 step is seven substeps, each a prediction of every
coordinate, one force evaluation, and a correction of the `g` and `b`
coefficients of every coordinate, repeated until the corrector
converges - 2 to 3 passes at binary64, 6 at binary128 and about 18 at
binary256 for the Kepler problem at a step of 0.05 (docs/VALIDATION.md).
Between those passes nothing leaves the step: the same positions,
velocities, accelerations and the seven `b`, `g`, `e` arrays are read
and rewritten, over and over. Per coordinate and per substep-pass the
loop engine issues about 110 elementwise operations; per binary256
step that is about 14,000 operations per coordinate on data that is
47 format-width values per coordinate. That ratio is the case for
keeping the state on the card.

What is NOT elementwise across coordinates is the gravity, which
couples every pair of particles, and the convergence test, which is a
maximum over all coordinates. Those two are where the host, or a
missing device mechanism, re-enters the loop.

## Operation counts, measured on the software backend

`ias15_cft` counts every library call. Per step, both problems, the
default settings of the sweep (`--max-iter 60`, corrector tolerance
1e-16 scaled by 2^(53-p)):

| problem | format | form | engine | calls per step | of which div/sqrt | corrector passes per step |
|---|---|---|---|---|---|---|
| Kepler, N=2 (6 coordinates, 1 pair) | fp64 | REBOUND's | loop | 1,868 | 255 | 2.39 |
| Kepler | fp128 | REBOUND's | loop | 4,321 | 612 | 6.0 |
| Kepler | fp256 | REBOUND's | loop | ~9,400 | ~1,200 | ~18 |
| Kepler | fp256 | FMA | loop | 9,400 | 276 | 18.1 |
| Kepler | fp256 | FMA | program | 3,698 | 276 | 18.1 |
| outer, N=6 (18 coordinates, 15 pairs) | fp64, adaptive | REBOUND's | loop | 2,804 | 380 | 3.1 |
| outer | fp256 | FMA | loop | 9,263 | 251 | 16.0 |
| outer | fp256 | FMA | program | 4,223 | 251 | 16.0 |

Two things in that table matter for the card. The corrector passes
grow with the format because the tolerance does: each pass gains about
2.7 decimal digits at this step size (the trace in docs/VALIDATION.md),
so binary256 needs seven times the passes of binary64 to converge to
its own precision. And the program engine removes 60 percent of the
calls: everything that was per-coordinate becomes two `cft_program_run_ex`
calls per substep-pass, and what remains is the gravity (about 20
vectorised calls per pass including the one square root and the one
divide), the end-of-step update, the step control and the scalar
bookkeeping.

The FMA form is what makes the program possible: REBOUND divides by
`rr[]` and by small integers at every level, and a correctly rounded
divide is not a sequencer program (cft-fp256 docs/ORBITS.md, obstacle
2). `--arith fma` multiplies by reciprocals and fractions that are
exact rationals rounded once (src/ias15_constants.h). It is a different
rounding sequence from REBOUND's - it differs at 3e-15 relative at
binary64 and 1e-70 at binary256 on the equivalence cases, at round-off
- and the records say which form ran.

## The resident design

One lane per coordinate for the predictor and corrector; one lane per
pair for the gravity; the ensemble dimension - E independent systems
- multiplies the lane count, which is what fills a 16-deep pipeline
and a four-tile card.

**State, per coordinate, held on the card across substeps and steps:**
`x0, v0, a0, csx, csv` and the seven-fold `g, b, e, csb, er, br` - 47
format-width values, 1,504 bytes at binary256; plus the working
`x, at` and the pair intermediates. For E = 1,000 six-body systems
that is 18,000 coordinates and 27 MB, a fraction of the U50's 8 GB.

*Do not read 47 as the archive's blob count, which is 50
(`CFT_N_BLOBS`). They count different things and neither is the other
plus a constant: the archive carries `csa0` and the live `x` and `v`,
which this list either omits or counts as working state, and the
sizing above has not been re-derived since `x` and `v` were added on
2026-09-10. Whether a resident design needs 47, 48 or 50 slots a lane
is a question for whoever builds it, and no card has run the resident
path.*

**Per substep-pass:**

1. `predict` - one program run over all coordinates. Streams `a, b, c`
   = `x0, v0, a0`; the lane's `b0..b6, csx` from its scratch block;
   `K0..K6` (the substep's `h_n (j+1)/(j+3)`), `fl(dt h_n)/2` and
   `fl(dt h_n)` from the bank; deposits `x`. 21 instructions
   (programs/predict-*.cfta).
2. gravity over pairs - one program run over all pairs, given each
   pair's two positions and two masses: `dx, dy, dz`, `r^2`, `1/r^3`,
   the two prefactors and the six contributions, deposited. **Its
   inputs are a gather** - pair l needs `x[i_l]` and `x[j_l]` from the
   predictor's deposits - which nothing on the device does today; the
   prototype's host does it with byte copies, and on the card it is a
   staged buffer per substep-pass or a device-side gather that does
   not exist (asks, below). `1/r^3` is the second fork: the correctly
   rounded route is `cft_sqrt` then `cft_div`, each composed of a
   program core with host prep and finish on the operand bits, so it
   is three round trips and a staged pass; the seed-and-Newton route
   (`CFT_RSQRT_SEED` plus a fixed number of FMA passes, cft-fp256's
   `--rsqrt newton`) is one program and different bits. The prototype
   keeps the correctly rounded route so that the FMA form's only
   departure from REBOUND is the reciprocal constants.
3. accumulate - each particle's acceleration is the sum of its
   partners' contributions in REBOUND's partner order (N-1 adds per
   coordinate). One program run over coordinates whose scratch block
   holds the gathered contributions; the same gather problem as step 2,
   in the other direction. The prototype issues N-1 vectorised
   `cft_run` adds with host gathers between.
4. `correct_n` - one program run over all coordinates. Stream `a` =
   `at`; the lane's `a0, g0..g6, b0..b6, csb0..csb6` from its 22-slot
   scratch block, written back; `1/rr` and `c` for the substep from
   the bank; deposits `tmp`, the change in `g_{n-1}`. 22 to 106
   instructions (programs/correct{1..7}-*.cfta).
5. at n = 7, the convergence test: `max |tmp| / max |at|` over all
   coordinates. A maximum is not one of `cft_reduce`'s trees (sum, dot,
   sumsq, sumabs), so the host reads the deposits - 3N x E values - and
   decides. On the card that is one `cft_buffer_from_device` of the
   deposit buffer per pass.

**Per step:** the end-of-step update (16 compensated additions per
coordinate, expressible as one more program over the scratch block;
the prototype issues it as elementwise calls), the PRS23 step control
(per-particle sums and a 20-iteration scalar Newton for the seventh
root; a few dozen scalar operations the host does in software, which
is the same bits by contract), and `predict_next_step`, elementwise
again.

**What the prototype proves.** `ias15_cft --arith fma --engine
program` runs steps 1 and 4 as the assembled programs through
`cft_program_run_ex` with the state riding in the per-lane scratch
block and the constants in a per-run bank, on the software backend,
and tools/check_program_engine.py shows the records bit-identical to
the host loop issuing the same operations one `cft_run` at a time, at
all three formats, on both problems, with fixed and adaptive steps.
The shape - state per lane, constants per run, the substep as two
program runs around a force evaluation - is therefore known to
compute the right bits. Its speed on the software backend is not the
point and is worse (the software program executor interprets each
lane's scratch traffic): 3.3 steps/s against 13.3 for the loop at
binary256 on the Kepler problem.

## Arithmetic intensity and the break-even

cft-fp256's measured constants (its docs/INTEGRATION.md, 2026-09-09,
one U50 tile at 135 MHz): a resident `cft_run` costs about 35 us fixed
and the engine moves about 107 M beats/s a tile, format-blind, so
107 M binary256 or 428 M binary64 element-operations a second; the
staged path moves 2.3 to 3.3 GB/s; four tiles are four times the
resident rate. The software backend on this host does 3.2 M binary64
FMA elements/s in bulk, and at the small n of one system about one
elementwise call a microsecond.

**Intensity.** With the state resident, what crosses the link per step
is the convergence deposit each pass (3N x E values), and the sample
when one is asked for. Per binary256 Kepler step at 18 passes that is
about 84,000 element-operations per system against 18 x 6 x 32 = 3.5 KB
of deposits read back: about 25 operations per byte, and 14,000 per
coordinate-step against INTEGRATION.md's break-even of about 4
instructions per element on PCIe. The staged elementwise route - every
`cft_run` on host pointers - is the opposite: 128 bytes moved per
element-operation, 1/128 of an operation per byte, bus-bound by two
orders of magnitude, which is why the loop engine as it stands must
never be pointed at a card.

**Time per step, projected.** Round trips do not scale with E; element
work does. Counting the design's round trips per substep-pass at four
(predict, gravity, accumulate, correct) plus the correctly rounded
`1/r^3` at three, and about ten per step outside the passes:

| format | passes per step | round trips per step | fixed cost at 35 us | element work per system-step | one tile, E = 1 | one tile, E = 1,000 | software backend, per system |
|---|---|---|---|---|---|---|---|
| fp64 | 2.4 | ~130 | 4.5 ms | 11,000 ops = 0.026 ms | ~4.5 ms | ~30 ms | 8.3 ms |
| fp256 | 18 | ~890 | 31 ms | 84,000 ops = 0.79 ms | ~32 ms | ~0.8 s | 75 to 106 ms |

Read across: a single system on the tile is 2 to 3 times faster than
the software backend at either format, and only because the fixed cost
is small; the gain the card exists for arrives with the ensemble, where
the projection is 250x (fp64) to 100x (fp256) over one core of the
software backend at E = 1,000, and the four-tile card four times that.
The binary256 column is engine-bound already at E of a few hundred -
84,000 operations per system-step at 107 M/s is the floor, 1,300
system-steps a second a tile - and that floor is the honest capacity of
this hardware for this integrator: a few thousand binary256 IAS15 steps
a second per tile, whatever the ensemble.

**Against a CPU running REBOUND itself.** REBOUND's own binary64 IAS15
does the Kepler step in 9.2 us on this host (100,000 steps in 0.92 s,
`build/ias15_ref`). The tile's fp64 floor is 26 us per system-step.
So at binary64 the tile is about three times SLOWER than one CPU core,
even fully loaded, and no design changes that: the tile is a binary256
engine, and IAS15 at binary64 is not what it is for. At binary256 the
comparison is against an MPFR-class CPU implementation, which this
repository did not build or time; the software backend's 75 to 106 ms
per step is a lower bound on how slow a naive one is, and the tile's
0.79 ms per system-step is 100x that.

## What could not be verified without hardware

- **The fixed cost of `cft_program_run_ex`.** The 35 us is measured
  for a resident `cft_run`; a program run stages its bank and its
  scratch blocks every time (cft.h, ABI 0.11: "its image, constant
  bank, counts and scratch blocks are staged always"), so each
  `correct` run moves 22 x 3N x E values in and out over PCIe. For
  E = 1,000 at binary256 that is 8.4 MB per run, 3.7 ms at 2.3 GB/s,
  per substep-pass: comparable to the element work itself. **The
  scratch block is where IAS15's state lives, and today it is not
  resident.** The single most valuable library change for this
  workload is to let `scratch_in` and `scratch_out` bind `cft_alloc`
  buffers the way `a, b, c` and `deposits` do; they are operand-shaped
  (n x count elements, lane-major, dense) so the mechanism exists.
- **The sequencer's rate on real programs.** 107 M beats/s is the
  elementwise engine's read path. A program instruction is one
  beat-operation per lane, but `LDL`/`STL` traffic, deposits and the
  16-lane block floor (a dependent chain runs at pipeline speed below
  16 binary256 lanes in flight, cft-fp256's docs/SEQUENCER.md) are not
  in that number. Not measured.
- **Gathers.** Steps 2 and 3 need a lane to read another lane's
  result. No device-side gather exists; the prototype's host does it.
  The asks are a device-side index-table copy, or a program-model
  change that lets a lane read a neighbour's deposit.
- **The maximum.** A `CFT_MAX` reduction would keep the convergence
  test on the card; today the host reads 3N x E deposits per pass.
- **Correct rounding of `1/r^3` on the card.** The composed route
  needs the host between its halves; the on-chip route is the
  seed-and-Newton composition and gives different bits. Which one a
  card run uses is a decision the bit-for-bit record must name.
- **Multi-tile.** The library partitions a `cft_run` over four
  compute units and the same holds for programs; the ensemble splits
  by lane index and each tile's bits are the software backend's, but
  the four-way projection above assumes the partition costs nothing
  and was not run.
- **Everything in the tables marked projected.** The operation counts
  and the software-backend times are measured; the card times are
  arithmetic on cft-fp256's published constants.

## What the card said, and the crossover

The record-diff task - `ias15_cft --engine program --artifact <xclbin>`
against the software backend, byte for byte - was run by the integrator
on the read-ahead quad (docs/VALIDATION.md, "the first run on a tile"):
the records are identical bit for bit, and the tile was 4.5x slower than
one core of the software backend on the two-body problem. (It was item 2
of this document's first to-do list; that list has since been rewritten
and renumbered, so the items below are not the ones those numbers named.) The
crossover was then measured on the same integration at three widths,
fp256, `--arith fma --engine program`:

| N | pairs per force call | software backend | card | card / software |
|---|---|---|---|---|
| 2 | 1 | 5.87 s | 26.93 s | 4.59x slower |
| 8 | 28 | 10.80 s | 41.80 s | 3.87x slower |
| 32 | 496 | 29.01 s | 22.65 s | 0.78x - the card wins |

So the card's deficit is per-call overhead and width amortises it,
with the crossover near 500 elements per call - well under the ~4,000
the naive overhead estimate above suggested. And difficulty alone does
not move the ratio: Burrau's Pythagorean problem drives the adaptive
step over three orders of magnitude and needs 18.2 corrector passes
against 15.07, but the software backend did 51,900 library calls a
second on it against 51,200 on the easy problem. Difficulty buys more
calls at the same rate; width is the lever.

## The ensemble on the tile

The ensemble - the third item of the first to-do list, and the step that
turns the design's lane count from 6 into 6,000 - is implemented
(docs/ENSEMBLE.md): E systems in one run,
3NE lanes for the predictor and corrector, E N(N-1)/2 for the gravity,
every member bit for bit its solo run (tools/check_ensemble.py, gated
at all three formats on the software backend). What it maps to on the
card, and what needs the card to confirm:

- **Width.** A Kepler ensemble of E members puts 6E elements into
  every predictor and corrector call and E into every gravity call.
  Against the measured crossover of ~500 elements per call, the
  predictor and corrector cross at E ~ 80 and the gravity's pair calls
  at E ~ 500; at E = 1,000 every call in the step is past it. The
  projection is the one above - 100x the software backend at E = 1,000
  for one tile - and it is a projection until the ensemble gate is run
  with `--artifact`. The software backend's own throughput against E
  is measured in the ledger (the ensemble throughput entry): it too
  gains from width, since its per-call cost has a fixed part, and that
  curve is the baseline the card's must beat.
- **Masking.** The prototype masks a member that has left the
  corrector, rejected its step or finished its block by a byte
  snapshot and restore on the host. On the tile that is a `SETACT`
  lane mask, which the sequencer has; but the host's snapshot is not
  how a card should do it, because the masked lanes still compute and
  their state still crosses the bus in the scratch block. The wasted
  work is measured by `pc_lane_efficiency` in every ensemble trailer:
  0.77 to 0.92 at binary64 on the gate's families, and at binary256,
  where the pass count runs from 10 to 23, whatever the ledger says.
  The ask is a per-run lane mask in `cft_run_args` (a host-supplied
  bitmap the engine honours), so that an idle lane costs neither a
  beat nor a byte.
- **Per-lane dt.** `predict-ens-<fmt>` carries `fl(dt h_n)` and its
  half in scratch slots 8 and 9, so an adaptive ensemble's predictor
  stages 10 slots a lane per run instead of 8 - 25 percent more
  scratch traffic on the predictor, none on the corrector, until the
  scratch block can be bound resident (the first ask below). A
  fixed-step ensemble runs the original, card-verified `predict-<fmt>`
  with more lanes and needs nothing new.
- **The scalar control.** The first ensemble build kept REBOUND's
  scalar control scalar: two compare round trips per lane per
  corrector pass for the `max` of the convergence test, one divide
  per system per pass, and about twelve scalar operations per system
  per step for the step control - some 220 library calls per
  system-step at binary256, which the software backend's throughput
  scan showed as a call count growing from 74,000 to 4.2 million
  between E = 1 and E = 1,024 (the ledger). On the software backend
  that was 3 percent of the time; on the card it would have been
  220,000 round trips of 35 us per step at E = 1,000, eight seconds
  against the 0.8 s of element work. It is gone: the maximum over a
  system's lanes and the minimum over its particles are selections
  by bit pattern on the host (the unsigned order of non-negative IEEE
  encodings is their numeric order - the same reading of bits the
  predicates already use, not arithmetic), the convergence quotient
  and the two exit tests are one E-wide divide and two E-wide
  compares per pass, and the step control is a dozen E-wide calls
  per step with masked inputs for the particles REBOUND skips and
  both branches computed for every system, each keeping its own -
  every element still the operation the system would issue alone,
  which the gates re-proved at every format. The step's round-trip
  count is now independent of E: 71,213 calls for 20 binary256 steps
  at E = 1 and 74,671 at E = 1,024 on the final build (the ledger's
  throughput entry), against 4,171,889 on the first. What remains for
  the card is that the convergence decision still reads 3NE deposits
  per pass: a `CFT_MAX` reduction over the deposits would keep it on
  the tile, and that ask stands. And the baseline is now measured:
  the software backend does 62 binary256 system-steps a second at
  E = 1,024 on one core (it gains 3.5x from width itself), so the
  engine-bound projection of 1,300 per tile is 21x that, not the 100x
  quoted against the E = 1 rate above.
- **Memory.** About 350 vectors of `max(3NE, E N(N-1)/2)` elements
  (52 state, ~300 broadcast constants) plus E snapshots: 60 MB at
  E = 1,000 and binary256, on the host; on the card only the operands
  of each call cross, so the constant width is bus traffic, not
  device memory, until `cft_run` gains a scalar-broadcast operand -
  a fourth ask, small, and worth 300 staged vectors a step.

## What to do first, if a card is available

0. Point the whole suite at it, which is now one variable:

       export CFT_REBOUND_ARTIFACT=/path/to/cft_hw_quad.xclbin
       make check

   No gate takes a flag for the card; every entry point falls back to
   this variable. Until 2026-09-10 nothing in the suite could open an
   artifact at all - the drop-in had a setter no gate called, and
   `struct cft_rebound_options` had no artifact field - so every
   hardware number in this document came from the standalone program.

1. Run the ensemble gate on the card, all three formats. The contract
   says the records are the software backend's; the record of a
   1,000-member ensemble on a tile is the proof, and it does not exist
   yet. Note that the python checkers spawn a fresh `ias15_cft` per
   case and pay a device open every time, so they are the expensive way
   to use a card; `build/check_dropin` and `build/gate_real` open the
   engine once per process.
2. Time the same Kepler ensemble at E = 1, 8, 64, 512, 4,096 members,
   fp256, `--engine program`, 20 steps, card against software backend,
   which is the ledger's throughput entry with the card column filled
   in. That curve is the whole economic case, measured.
3. Time `cft_program_run_ex` on `correct7-fp256.cftp` over 6,000 lanes
   with the 22-slot scratch block, staged, against the 35 us + 3.7 ms
   projection above; and `predict-ens-fp256.cftp` with its 10. Those
   numbers decide whether scratch residency is the first library job.
4. A mixed ensemble (members needing 10 and 23 passes in the same
   step) to price idle lanes on the tile in wall time, which is what
   decides whether the lane-mask ask is worth its library change.

## The specialised bitstream: tried, measured, and worth 2%

**2026-09-13/14.** Everything above treats the tile as given. It is not:
cft-fp256's kernel carries trim generics, so an image can be built with
the binary256 rung left out - smaller, faster, and refusing binary256 by
name. That was built, on the reasoning that this integrator's own
precision study (docs/HORIZON.md) says binary128 is the rung worth
having. Three images now exist and have been measured against each other
on the card. **docs/BITSTREAM.md is the full account; docs/VALIDATION.md
entries 34 to 36 are the measurements.** The conclusion belongs here,
because it changes what this document's to-do list is for.

What the hardware did, all of it as designed:

- one tile, no binary256, closes **150 MHz** where the full tile ships at
  135, with a quarter of the tile gone (32,398 LUTs and 143 DSPs);
- its engine runs **1.108x** the full tile's rate, measured back to back,
  which is 99.7% of the clock ratio;
- six such tiles fit and close at 125 MHz and deliver **1.39x** the
  shipped quad's aggregate arithmetic, against 1.389 predicted;
- every record byte-identical to software across one, four and six tiles
  and two different bitstreams.

What it did for IAS15:

- **more tiles is slower, monotonically** - one tile beat four beat six
  in all ten card rows, at every body count and both formats. The
  library partitions every `cft_run` across every tile a device
  presents, so one call becomes a kernel launch and a staging round per
  tile, and this integrator issues thousands of small calls per step.
  Tile count is a divisor on an already short vector and a multiplier on
  a fixed cost;
- **the specialised image itself is worth 1 to 3%**, isolated by running
  both one-compute-unit images head to head.

**So the arithmetic is about a fifth of the wall clock.** An 11% faster
engine moving the integration by 2% puts it there directly, and that is
the single most useful number this document now contains: **four fifths
of a card-side IAS15 step is not arithmetic**, and no bitstream
addresses any of it. The items in "What to do first" below are therefore
not a list of nice-to-haves beside a hardware lever. They are the whole
lever.

Read the ranked asks in that light:

1. **State resident on the card across a step.** The engine's vectors are
   plain `calloc` today (`valloc` in src/ias15_cft.c), so every operand
   of every call is staged across PCIe and staged back. libcft has
   carried device-resident `cft_alloc` buffers since ABI 0.11, and
   cft-fp256 measured 3.0 to 3.3x on the raw path from exactly this
   change. Nothing in this repository uses it yet. **This is the
   cheapest large item and it needs no hardware.**
2. **Fewer, larger calls.** The program engine already removed 60% of
   them; the remainder is gravity, the end-of-step update and the step
   control. Every call removed is a launch and a staging round removed
   per tile.
3. **A device-side scatter for gravity's accumulate half**, which is
   36-42% of card wall clock as N-1 narrow vector adds, and projects
   1.99x to 3.1x. This is the one that needs a sequencer feature rather
   than a library change.
4. **The convergence test, removed exactly** via a recorded corrector
   schedule - 13-19% of calls, worth 4-9% on the card and nothing in
   software.

And one negative result worth carrying: **do not reach for a wider
image to make this faster.** It was tried. Measure the partitioning cost
with the images that already exist before building another one; two
existing images predict the whole outcome in an hour, and two multi-tile
links cost thirteen.

## Moving the state onto the card: what it would actually take

The section above designs a resident integrator. This one prices the
gap between that design and `--engine program` as it stands, because
the 2% result makes this the work that matters. Written 2026-09-14
after reading both sides of the seam.

### What the library already does, verified rather than assumed

The pinned libcft (`ca19fe3`, ABI 0.12) binds **every** operand of
`cft_program_run_ex` through the device-resident buffer registry - the
three streams, the deposit window, **and both scratch blocks**
(`host/src/device.c`, `bind_role` for `CFT_ROLE_SI` and `CFT_ROLE_SO`).
Its comment names this document's first ask as the reason. So the
library-side item this file has carried since 2026-09-10 is **done**,
and nothing below is blocked on cft-fp256.

What is missing is on this side: `valloc` in src/ias15_cft.c is
`calloc`, so not one buffer in this port is a `cft_alloc` buffer, and
the registry has nothing to recognise. Every call stages.

### What a substep-pass moves today, in elements of `N3`

`run_program` rebuilds its scratch block on the host from `b[]`, `csx`,
`a0`, `g[]` and `csb[]` before every run, and reads it back after:

| per run | host memcpy | staged over PCIe |
|---|---|---|
| predictor | 8 x N3 in, N3 out | 12 x N3 |
| corrector (each of 7) | 22 x N3 in, 22 x N3 out | 46 x N3 |

A corrector pass is seven of each, so about **406 x N3 elements cross
the bus per pass** and roughly the same is memcpy'd on the host first.
At 512 bodies and binary128 that is 10 MB a pass and about 80 MB a step,
each way - and the memcpy is pure waste, because the values it moves
came from the device and are going straight back to it.

### The three changes, cheapest first

1. **Allocate through `cft_alloc` (hours).** `x0`, `v0`, `a0` are
   written once a step and read fifty-odd times; `at` once a substep.
   Making those and the four program buffers resident stops the streams
   being re-staged. It is perhaps 7% of the traffic on its own - the
   scratch block is the bulk - but it is the prerequisite for everything
   else and it makes residency *measurable*, through
   `cft_buffer_get_info`'s resident and staged bind counts.

2. **Keep the scratch block on the device between passes (the real
   one).** The corrector reads `g/b/csb` out of `scratch_out` and writes
   the same values back into `scratch_in` for the next pass, through the
   host, 44 x N3 elements at a time. If the program's output layout is
   its own input layout, `scratch_out` **is** next pass's `scratch_in`
   and the host never touches it: cft-fp256's `resume-fp64` program
   demonstrates exactly that shape, two runs through one block with
   identical hashes. The host would then read the block once at the end
   of a step rather than fourteen times a pass, for the end-of-step
   update, the step control, `dpcast` and the checkpoint. This removes
   most of the 406 and all of the matching memcpy, and it is the change
   docs/HARDWARE.md's resident design has always described.

3. **The convergence test, which survives both of the above.**
   `pc_error` reads `tmp` and `at` in full, element by element on the
   host, every pass, so even a fully resident corrector pays one
   `from_device` of 2 x N3 per pass. Two ways out, and they are
   different in kind: the **recorded corrector schedule** replays a
   known pass count and removes the test exactly, bit for bit (13-19%
   of calls, and it has been costed here before); or a device-side
   maximum, for which **ABI 0.12 added `CFT_MAXALL`** - though
   `pc_error` maximises over *normal* values only and per system, which
   that reduction does not express, so this route needs thought rather
   than a substitution.

Gravity's gather and scatter (steps 2 and 3 of the resident design) are
untouched by any of this and still need a device-side mechanism that
does not exist. They are the next wall, not this one.

### The measurement that should come first

None of the above says which of per-call latency, staged bytes and the
host's own gather actually dominates - the byte counts here are
analytic, and the wall clock is measured, and nothing in this repository
connects them. `ias15_cft` already counts library calls; adding staged
bytes beside that count, and timing a step with the gravity stubbed out,
would turn this ranking from an argument into a number. Do that before
writing any of the three, because a day spent on residency when the wall
is the gather would be the same mistake the six-tile image already made
once.

### Measured, later the same day: the ranking above was wrong

The section above asked for one measurement before anything was
built, and the measurement came back the same afternoon
(`hw/bench-workload.py`'s trailer instrumentation, docs/VALIDATION.md
entry 37): every library call timed on a wall clock and bucketed,
gravity timed inclusively, bytes presented per call summed. One f128
tile, the FMA form on the program engine, 2026-09-14.

Share of the card's wall clock, per bucket:

| problem | format | wall | program runs | elementwise | div and sqrt | gravity, inclusive |
|---|---|---|---|---|---|---|
| Kepler, 2 bodies | binary64 | 24.1 s | 8.1% | 77.3% | 13.0% | 59.8% |
| Kepler, 2 bodies | binary128 | 46.8 s | 10.5% | 75.9% | 12.2% | 71.7% |
| n-body, 64 | binary64 | 9.4 s | 3.7% | 63.6% | 30.6% | 82.7% |
| n-body, 64 | binary128 | 25.3 s | 4.3% | 54.8% | 39.0% | 89.9% |
| n-body, 256 | binary64 | 27.3 s | 0.9% | 48.5% | 48.0% | 73.5% |
| n-body, 256 | binary128 | 59.3 s | 1.1% | 37.4% | 59.2% | 85.8% |

Gravity overlaps the other three - it is timed around its own calls -
which is why it is a column and not a segment.

**The program runs are at most a tenth of the step, and one percent at
256 bodies.** Every byte the residency items above would keep on the
card moves inside those runs. Items 1 and 2 of the ranking above are
therefore worth 1 to 10%, and the argument that put them first - that
406 x N3 elements cross the bus per pass - was true and beside the
point: at 256 bodies and binary128 the run presents 2.4 GB a step, and
at the bus rate this project has measured that is about 0.8 s of an
11.9 s step. Bytes are not the wall. Two other things are.

1. **Per-call cost, at every size.** A Kepler step issues 172,000
   elementwise calls of six elements each over 200 steps, and on the
   card they cost **108 microseconds apiece** against 0.73 in software -
   a launch, a staging round and a read-back, for six numbers. That is
   77% of the step at two bodies. It is not a residency problem; a
   resident six-element buffer still costs the launch. It is a
   call-count problem, and the fix is the one the design section has
   always named: more of the substep as a program, so that gravity and
   the step control stop being hundreds of separate calls.
2. **The correctly rounded divide and square root, once the vectors are
   long.** At 256 bodies and binary128 the 674 `cft_div` and `cft_sqrt`
   calls of a five-step run - `1/r^3` for every pair, once per force
   evaluation - take **35.1 of 59.3 seconds**, 52 ms a call over 32,640
   pairs, **1.6 microseconds an element**, against 4.3 nanoseconds for
   an FMA on the same tile. In software the same route is 7.6
   microseconds an element and 71% of the step. cft-fp256's own peers
   tool prices the software square root at binary128 at about 10,900 ns
   against MPFR's 49, so this is a property of the composed
   correctly-rounded route - a program core with host prep and finish
   on every operand - and not of the tile. The seed-and-Newton route
   (`CFT_RSQRT_SEED` plus fixed FMA passes) is one program and different
   bits; a device-side correctly rounded divide is the ask that keeps
   the bits.

So the order for the card, measured rather than argued, is: **cut the
call count, then fix the divide, then residency.** The residency work
was the right thing to cost and the wrong thing to start with, and the
paragraph above that said "do that before writing any of the three"
is the one part of the section that was right.

Gravity's own gather and scatter, which the design section treats as
the next wall, are inside the inclusive column; at 256 bodies the host
byte-copies and the accumulate adds are the part of that 86% which is
neither the divide nor the launches, and this instrumentation does not
split them out. That is the next measurement.
