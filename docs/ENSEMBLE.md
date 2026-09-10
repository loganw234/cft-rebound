# Ensembles: E systems in one run

The design of `ias15_cft`'s ensemble mode, the decision about step
sharing and what it costs, and the gate that says it is right. The
measurements are in docs/VALIDATION.md (the ensemble entry) and the
hardware projection in docs/HARDWARE.md.

## Why

The first run on a tile (docs/VALIDATION.md, 2026-09-09) settled
correctness and said the expected thing about speed: a two-body
problem is six coordinates and one pair, the engine's width is its
whole advantage, and six lanes leave it idle while every program run
pays its fixed cost. The integrator then measured the crossover on the
card: the tile loses to the software backend by 4.6x at one pair, by
3.9x at 28 pairs, and wins at 496 pairs (N = 32). Width is the lever.
Difficulty is not - a harder problem buys more calls at the same rate
on both sides.

An ensemble is the width a two-body problem does not have. E
independent systems integrated in one run put 3NE coordinates in the
predictor's and corrector's vectors and E N(N-1)/2 pairs in the
gravity's, so the one-pair problem becomes an E-pair call. And it is
not a benchmark trick: an ensemble of nearby initial conditions is how
chaotic systems are studied, how sensitivity to initial conditions is
measured, and how a Lyapunov rate is estimated; tools/divergence.py is
that measurement and docs/VALIDATION.md runs it on the Pythagorean
problem.

## The layout

A problem file may carry `E` and `system NAME` lines
(tools/make_ensemble.py writes them; a file without `system` lines is
one system). System s occupies coordinates `3N s .. 3N (s+1) - 1` of
every per-lane vector (`x, v, a, x0, v0, a0, csx, csv, csa0, at` and
the seven-fold `g, b, e, csb, er, br`) and pairs `PS s .. PS (s+1) - 1`
of the pair list, where PS = N(N-1)/2 and the pairs within a system
are visited in REBOUND's order `i = 1..N-1, j < i`. The gravity's
accumulation walks every particle of every system at once, each
particle adding its own system's partners in REBOUND's partner order.

Everything scalar in REBOUND's step - the step `dt`, `dt_last_done`,
the time, the corrector's error estimate and its previous value, the
iteration count, the accept/reject decision and the ratio for
`predict_next_step` - is an array of E scalars. Where the scalar used
to be broadcast into every lane (`DTB`), it is now broadcast per
system: each lane carries its own system's `dt`, `dt_done` and step
ratio. Because `cft_run` is elementwise and each element depends on
its own inputs only, an operation over the shared vector is, lane by
lane, the operation the system would have done alone with the same
scalar - the same bits.

The FMA form's per-substep factors `fl(dt h_n)` and `fl(dt h_n)/2`
become per-lane vectors too (one multiplication over the lanes instead
of one scalar multiplication and a broadcast: the same bits per lane).
The sequencer's predictor kept them in the constant bank, which is
per run, so an ensemble with per-system steps uses a variant program,
`predict-ens-<fmt>`, that reads them from two more scratch slots (10
instead of 8) and keeps only the seven format-wide constants in the
bank; the corrector programs do not touch `dt` and are unchanged. A
fixed-step ensemble shares one `dt` and runs the original, tile-
verified `predict-<fmt>` with more lanes.

## The step-sharing decision

The point of IAS15 is its adaptive step, and the members of a real
ensemble want different steps: an ensemble spanning eccentricities, or
one that is watched through a close encounter that some members meet
and others do not. Three designs were on the table.

1. **One shared step for the ensemble**, chosen from all members
   (the minimum, say), with one shared corrector exit. Simplest, and
   the one a fixed-step run reduces to. It changes the science: each
   member is integrated with a step that is not its own, so its record
   is not the record of that system run alone, and the whole ensemble
   crawls at the pace of its hardest member. Its corrector exit is
   wrong for every member but the slowest to converge.

2. **Per-system steps by running the members through separate
   vectors** - no shared vector, no gain. That is what exists already.

3. **Per-system steps in the shared vector, with masking.** The vectors
   are shared, the scalar control is per system, and a system that
   has left the corrector, rejected its step, or finished its sample
   block is snapshotted at that moment (byte copies of its slices of
   every per-lane vector, 52 of them) and restored after the shared
   vector work. The arithmetic on a masked system's lanes still
   happens - the vector is the vector - and is thrown away.

The third is what is built, because it keeps the science exactly:
every member's record is bit for bit the record of that member run
alone, adaptive step and all, and that is a gate rather than a hope
(below). What it costs is idle lanes. The shared corrector loop runs
until every active system has converged, so a pass is issued as long
as any system still needs one; a system that converged in 13 passes
rides along, masked, while a neighbour takes 23. The trailer reports
this as `pc_lane_efficiency`, the fraction of issued lane-passes that
were live. On the gate's Kepler families at binary64 it is 0.77 to
0.92; at binary256, where the pass count is 10 to 23, it is measured
in the ledger. On the tile the same masking is a `SETACT` lane mask
rather than a snapshot, and the wasted passes are wasted engine beats,
not bus traffic. A rejected step in one member is the same mechanism:
that member is restored to its step start and re-predicts from its
`er, br` with `dt_new / dt_last_done` in the same vector pass in
which the accepted members predict from `e, b` with `dt / dt_done`,
its lanes carrying its own ratio.

Members are synchronised by accepted step count, not by time: every
member does `steps` accepted steps and reaches its own time. That is
what makes the per-member record identical to the solo record (which
is also `steps` accepted steps), and it is the honest unit for an
adaptive integrator; a member that takes small steps through an
encounter covers less simulated time than its neighbours in the same
run. A user who wants every member at the same final time reads the
per-system times in the record and runs the block again, or uses a
prescribed shared sequence (`--dt-file`), which is what the Lyapunov
measurement wants anyway: members that share their sampling times.

The first version was going to be a fixed shared step. It turned out
that per-system steps cost about as much code as the masking the
shared corrector exit needed anyway, and the fixed-step ensemble is
now simply the case where all E steps are equal.

## What the ensemble does not do

- It does not share anything across members: no common centre of
  mass, no coupling. E systems of N bodies is not one system of EN
  bodies (that would be `N` = EN in one system, which the code also
  runs, one pair list of EN(EN-1)/2).
- `G` is per file, shared by every member; masses are per body and
  may differ between members.
- The step control's per-particle timescale is still issued particle
  by particle (four scalar calls and two compares per particle per
  step, as REBOUND's loop does it), so an ensemble of E two-body
  systems issues about 12E scalar library calls per step for the step
  control on top of the ~3,700 to 9,400 vector calls of the step
  itself. That is under 1 percent of the software backend's time at
  E = 1,000 and it is the first thing to vectorise if it ever shows;
  on the tile the honest fix is a `CFT_MIN` reduction (docs/HARDWARE.md).
- Memory is 52 per-lane vectors plus about 300 broadcast constants,
  each `max(3NE, E N(N-1)/2)` elements wide, plus the E snapshots:
  about 60 MB at E = 1,000, N = 2, binary256; 600 MB at E = 10,000.
  `cft_run` has no scalar-broadcast operand, so the constants are
  really that wide.

## The gate

tools/check_ensemble.py: for each case, the ensemble is run once and
each member is run alone with `--member k` on the same file with the
same arguments, and every recorded value of every sample - time, next
and last step, energy, every coordinate, the exact-time pair - is
compared as an exact rational, together with each member's step
statistics (steps rejected, corrector passes mean and maximum, cap
hits). Cases: a five-member Kepler family whose members have
different eccentricities (planet velocity scaled by 1 - 0.08k), fixed
and adaptive, on the loop engine and on the program engine (the
`predict-ens` program); a four-member family a few ulps apart; a
five-member family spread in starting radius from 0.5 to 2.1 started
with a first step of 0.6, which the four inner members reject in the
same attempt in which the outermost accepts it, on both engines; and
a three-member outer solar system with Jupiter displaced by 2^-20 AU,
fixed and adaptive. The result at every format is in docs/VALIDATION.md.

tools/check_records.py guards the other side: the single-system path
through the same code must still reproduce the committed records of
the sweep (results/raw/) bit for bit, which it does at all three
formats.

## Using it

    python tools/make_ensemble.py data/problems/kepler.txt --members 64 \
        --out ens.txt --offset planet x 0x1p-52
    build/ias15_cft --format fp256 --problem ens.txt --epsilon 1e-9 --dt 0.01 \
        --steps 10000 --sample 250 --max-iter 60 > ens.rec
    build/ias15_cft --format fp256 --problem ens.txt --member 17 ...    # member 17 alone
    python tools/oracle.py ens.rec        # scores one system: see --system in tools/state_to_problem.py
    python tools/divergence.py ens.rec    # members against member 0 (shared times: use --dt-file)

The record of an ensemble carries one `system s sample ...` line per
member per sample and one `# system s steps_done=...` trailer per
member, then the usual trailer with the totals and
`pc_lane_efficiency`; a `--member` run writes the ordinary single-
system record, which is what the gate diffs against.
