# IAS15 as a resident workload on the tile: the design, the numbers, and what is unverified

Everything below is a design plus a software-backend prototype. **No
card was opened, no XRT call was made, and no `cft://` server was
contacted for this document**; both were in use. Every card number
here is either quoted from cft-fp256's own measurements (docs/
INTEGRATION.md and docs/BENCHMARKS.md of the pinned clone, dated
2026-09-09) or is a projection built from those and from operation
counts measured here on the software backend. Projections are marked
as such, each time.

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
  16 binary256 lanes in flight, docs/SEQUENCER.md) are not in that
  number. Not measured.
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

Item 2 of the original list below was run by the integrator on the
read-ahead quad (docs/VALIDATION.md, "the first run on a tile"): the
records are identical bit for bit, and the tile was 4.5x slower than
one core of the software backend on the two-body problem. The
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

Item 3 is now implemented (docs/ENSEMBLE.md): E systems in one run,
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

1. Run tools/check_ensemble.py with `--artifact`: the ensemble gate on
   the card, all three formats. The contract says the records are the
   software backend's; the record of a 1,000-member ensemble on a tile
   is the proof, and it does not exist yet.
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
