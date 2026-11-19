<!-- SPDX-License-Identifier: GPL-3.0-or-later
     Copyright 2026 the cft-rebound contributors. -->

# The ring run's CPU half: what a GPU wave will be held to

**As of 2026-11-12, at `d9c8c8c` and the docs that land with it; round
5's sections as of 2026-11-18, at `fda3b40`.**
cft-fp256's study GPU-E (its docs/studies/GPU-E-mi300x-ring.md, merged at
`56006987`) sized an IAS15 run of a Saturn ring on an MI300X: 10^6 test
particles in 64-body waves, each wave a system of Saturn with J2 and J4,
six moons and 57 particles. Round 4 (ROUND4.md) built on the CPU what
that run computes and what it will be checked against. Round 5
(ROUND5.md) taught the driver to take a producer's own files, and held
this tree on a GPU as it is. This file says what exists, what holds it,
and what a GPU run must carry to be held. docs/VALIDATION.md entries 46
to 51 are round 4's measurements, and entries 54 on are round 5's.

Every figure here is MEASURED on this desktop's CPU, on the software
backend, unless it is labelled: READ (in a named source), COMPUTED
(arithmetic on the others) or BELIEVED (taken on trust, not checked).
Nothing in round 4 ran on a GPU, a card or WSL. A figure from round 5
names the GPU it was measured on.

## The ring

### What it holds

| | value | source |
|---|---|---|
| Saturn | GM 37.931e6 km^3 s^-2; equatorial radius 60,268 km | READ in GPU-E section 4: NASA's Saturn fact sheet |
| J2, J4 | 16,290.573e-6 and -935.314e-6 | READ there: Cassini, Iess et al. 2019, in Fortney et al., arXiv:2304.09215, Table 4.1 |
| R, the radius J2 and J4 refer to | 60,330 km | BELIEVED: GPU-E's belief, since its table names none |
| six moons | Mimas, Pan, Prometheus, Pandora, Epimetheus, Janus: GM, a and e | READ in GPU-E section 4: JPL SSD |
| the moons' inclinations, nodes, pericentres and phases | small, not zero, fixed in `tools/make_ring.py` | synthetic, ROUND4.md decision 5 |
| the annulus | 130,000 to 142,000 km | GPU-E's answer 1 |
| particles | massless, inclined 1e-4 to 1e-3 rad, from `random.Random(seed)` | synthetic |

Units are km and s with G = 1, so every mass is a GM.

### How it is generated

    python tools/make_ring.py OUT.txt --particles N [--seed S] [--one-system]
           [--no-zonal] [--rmin KM] [--rmax KM] [--exclude-hill K] [--name NAME]

- **Saturn, then the six moons, then the particles,** so `active 7` is
  right for every system. The file carries `active 7` and
  `zonal 0 R J2 J4` unless `--no-zonal` is given.
- **Waves.** The particles are sorted by radius and cut into waves of 64
  bodies, 7 massive and 57 particles, GPU-E's shape a2; `--one-system`
  writes them as one system instead. For waves, a count that is not a
  multiple of 57 is refused by name; `--one-system` takes any count.
- **Exact.** Everything is computed in mpmath at 160 bits from the
  binary64 constants, and each written value is rounded to binary64 once
  and written as exact hex. Python's trig is not used. The moons'
  states are REBOUND's own `reb_particle_from_orbit`, and the frame is
  the barycentre of the 7 massive bodies.
- **Circular in its own field.** With the zonal line, each particle has
  the circular speed in Saturn's equatorial plane,
  v^2 = (GM/r)(1 + 3/2 J2 (R/r)^2 - 15/8 J4 (R/r)^4), and each moon's
  velocity is scaled by the same factor at its radius, so Pan stays
  circular. With `--no-zonal`, a point mass's speed. So a file's zonal
  line is isolated by deleting it from that file, never by comparing a
  `--no-zonal` file: COMPUTED by P2, a with-zonal file read in a point
  mass's field puts its particles on orbits of e 0.0045 to 0.0053, out
  to 143,253 km, past the annulus.
- **The exclusion.** A start within 5 Hill radii of a moon is drawn
  again, and the file's head counts the refusals by moon. More than 99
  refusals a particle, and 1,000 beside, is refused by name.
- **The same bytes from the same arguments.** V2 ran it ten times, under
  another PYTHONHASHSEED, without gmpy and to another path: one sha256.
  An independent 400-bit recomputation reproduced every value, 399,049
  of them in a 57,000-particle ring with its zonal line. BELIEVED: the
  same bytes on every host; the check across hosts waits for a WSL slot.

### At the study's size

`--particles 1000000` is refused by name: 10^6 is not a multiple of 57.
`--particles 1000008`, 17,544 waves of 57, wrote 174,691,521 bytes. Its
particles run from 130,000.025 to 141,999.987 km, the largest |z| from
Saturn is 141.850 km, and 49 of 1,000,057 draws were refused: Pan 2,
Prometheus 19, Pandora 28 (P2, 2026-11-11). An ensemble's systems all
have N bodies, so GPU-E's last wave of 49 particles cannot be written
here. The run is 17,544 full waves and 1,000,008 particles, 17,543 waves
and 999,951, or a padded last wave. Nothing integrated the ring at that
size.

## What holds each part

Each gate is one leg in all three suite tiers.

| part | gate | what it holds | entry |
|---|---|---|---|
| a problem's `active n` line | `tools/check_gravity.py`, parts 2b and 2c | `active 5` gives every body active's record at every format, on both engines and in both forms, at 25 pairs a force evaluation against 28; a test particle's 2^-40 mass moves nothing, the energy included, and moves everything without the line; the refusals | 46 |
| J2 and J4, a problem's `zonal k R J2 J4` line (`src/ias15_zonal.inc`) | `tools/check_zonal.py` | every traced operation against cft-fp256's golden model, bits and flags; the term against mpmath's gradient of the potential at 200 digits; energy, momentum, and the apsidal and nodal rates of first-order theory, J4 resolved; binary128 against binary256, with its scope; both engines, an ensemble, the run key | 47 |
| the ring generator | `tools/check_ring.py` | the same bytes; circular in its own field; every wave and every moon out of the plane; one wave at every format on both engines; `active 7` against all 64 active, 420 pairs against 2,016; REBOUND's own IAS15 on a wave, bit for bit at binary64; the ensemble against its members; the far wave's steps an orbit | 48 |
| one particle out of its wave | `tools/check_replay.py` | three waves at every format, rejecting nothing, first, and after accepted steps; each sampled particle and every massive body at the wave's bits, with its counters; the records load bearing; refusals; a dt file through a pipe | 49 |
| the driver | `tools/check_compare_run.py` | its PASS on a nine-wave ensemble, and a change planted in every column it compares, for members and for particles, each failing by name: bits, not values; a record of another format, arithmetic or problem refused, and one lacking a field the rerun needs (as REBOUND's own ias15_ref's does); another producer's record that carries every field is held like this tree's; and, with `--producer-files`, particles replayed from a producer's own files, the same as from a rerun's, a producer whose exit test is not this tree's, each file refused by its run key when it is another wave's, and the limit's line word for word | 50, 51, 57 |
| the parts together | `tools/check_ring_seams.py` | part 1, the ring in its own field, with the node counter (P1 x P2); part 2, a ring particle replayed with its zonal line (P2 x P3, P1 x P3); part 3, the driver on the ring's own shape (P2 x P4, P1 x P4) | 51 |

Where a gate compares one run of this tree with another, a defect both
runs share cancels; of the runs above, only REBOUND's own IAS15, at
binary64, is another build. So the zonal term is held to references of
its own, the golden model and mpmath, the physics counters compare with
theory, and the format difference below says what it cannot see.

## What a GPU run's record must carry to be held

    python tools/compare_run.py RECORD --problem FILE --format F --arith A
           [--exe build/ias15_cft] [--programs programs/out]
           [--seed 1] [--strata 8] [--format-diff F2] [--plan]
           [--rerun-engine loop|program] [--keep DIR]
           [--producer-files DIR [--particle-strata P]]

The driver reads only the record format, the lines `ias15_cft`'s
`main()` writes, from its first line to the trailer. An engine that
writes the same lines is held with no change to the driver: V4 rewrote a
record as `program=gpu impl=hip-6.2 backend=hip:0`, with one extra header
line and one extra trailer line, and it passed.

**What it must carry.**
- **The header the rerun is built from.** On `# program=`: format,
  arith, cs, engine, pc_tol_shift, max_iter, steps, sample, N and
  problem. On `# dt0=`: dt0, epsilon, G and pc_tol. The `# columns:`
  line, and a problem's `# active` and `# zonal` lines. A header missing
  a field is refused by name.
- **The format and the arithmetic,** which the caller names with
  `--format` and `--arith`. A record of another one is refused: a genuine
  binary128 record reruns at binary128 and agrees, so the caller says
  which run it means. Records in the two forms are not comparable: on a
  ring wave at binary64 the forms' dt_next differ by 1.7e-11 after 3
  steps (P2).
- **The problem file it ran,** held to the record before anything runs:
  its name, G, systems, every body's name and mass, its `active` line,
  and every system's sample 0. A test particle's tiny mass, or the active
  line dropped from massless ones, moves no bit, so only this check sees
  it.
- **Every member whole:** its system line, bodies, a sample line at
  every step steps and sample make, and a trailer with steps_done,
  steps_rejected, iterations_max_exceeded, mean_pc_iterations and
  max_pc_iterations; the members' rejections must sum to the record's
  total.
- **Bits.** A sampled member's lines are compared token for token: the
  sign of a zero and a digit past binary64 are differences. The
  exact-time pair, t_hi and t_lo, is compared too. It comes from
  `vaugadd`, the API call `cft_augmented_add`, which at the pin runs on
  the host on every device, at every step attempt: no `hip:0` run meets
  AUGADD (VALIDATION entries 54 and 55; this bullet said BELIEVED, refused
  as AUGADD on hip:N, until round 5).

**The run settings a record does not name.** `--adaptive-mode`,
`--softening`, `--min-dt` and `--no-flag-abort` are not in the record,
though the corrector schedule's run key names the first three (P4, read
in `main()`). The rerun takes their defaults, so a run made otherwise
passes only where the setting moved no bit. A GPU run's record must name
any non-default itself, and the driver does not read one yet: NOT BUILT.
REBOUND's `ias15_ref` prints adaptive_mode=2 in its record; `ias15_cft`
does not.

**The producer.** A record's first line does not name it: `ias15_ref`
writes the same `# cft-rebound ias15 record v1`. The driver refuses a
record that lacks a field it reruns from, as `ias15_ref`'s does, and
`--rerun-engine` reruns a record whose engine this tree does not have.

**What the sample trusts.** Members are sampled by seed, one from each
stratum of member index, which is radius, plus every member whose
trailer shows a rejection. MEASURED by V4: a member whose trailer
under-reports its rejections, the record's total adjusted to match,
leaves the sample, and an unsampled member's lines are never rerun. The
driver prints that beside the sample.

**Particles from a rerun.** Without `--producer-files` the driver
samples two test particles of each member it reruns, the least closest
approach to an active body at the samples and one by the seed, and
replays each alone through `tools/extract_particle.py` from the attempts
record and corrector schedule this tree's rerun wrote. So a particle is
sampled only in a member that was rerun.

**Particles without a rerun: a producer's own files** (round 5, P3;
VALIDATION entry 57). With `--producer-files DIR`, particles are sampled
in members the driver does not rerun. For every wave, the producer
writes in DIR:
- `member<k>.att`, the attempts record in `--attempts-out`'s form (each
  attempt's dt and outcome, then the next dt), and `member<k>.sched`,
  the corrector schedule in `--record-iters`'s (each attempt's passes),
  k being the member's index in the problem file, and 0 for a single
  problem;
- each as that member run alone writes it, keyed on its `# run` line
  with what `ias15_cft --run-key` prints for the arguments the driver
  reruns the member with: its engine the record's own if this tree has
  it, and `--rerun-engine`'s if not. A single problem with no system
  lines keys `member=-1`.

The particles sampled are the member sample's, and one member's from
each of `--particle-strata` ranges of member index (default:
`--strata`'s).

**What the files let the driver hold.**
- **That each file is its member's run's,** by its key, before anything
  runs. A file of another wave, format, `active` line or `zonal` line is
  refused, naming the member, the file and every key field that differs;
  a file not of its form is refused by its first line.
- **The arithmetic that reaches each sampled particle's samples and
  trailer,** with no rerun of its member: every sample's head fields and
  exact-time pair, the particle, every massive body and the trailer's
  counters, held to the record's own lines. Only sample lines and the
  trailer are compared, and a rejected attempt after an accepted step
  moves no state bit (VALIDATION entry 49): its dt shows only as the
  dt_next of a sample line printed just before it. MEASURED with
  ias15_cft as the producer: the same 18 particles from its files as
  from a rerun's, at every format, 1,512 tokens each time.

**What they do not.**
- **The producer's decisions.** A replay takes each attempt's dt and
  outcome, and its passes, as given. Where a producer's exit test stopped
  two passes early on a first attempt, the driver holds that member, NOT
  RERUN, both its particles IDENTICAL (P3's GATE 1b, MEASURED; that run
  ends rc 1 because another member, rerun, FAILs). Where nothing else
  differs the run says RESULT: PASS (V3, MEASURED on a three-wave ring).
  The same member FAILs once it is rerun. Decisions are held only in a
  member that is rerun: one a stratum, and every member with a rejected
  step.
- **A NOT RERUN member's unsampled particles,** which are held by
  nothing.
- **The producer's settings.** The key holds the files to the arguments
  the record names, not to those the producer ran with. Its epsilon,
  pc_tol_shift and max_iter, and the settings no record carries, are held
  only in a member that is rerun, and there only through what they move
  in its lines. MEASURED by V3: produced at pc_tol_shift 2 or 8, or at
  max_iter 12, a member whose header names the defaults is IDENTICAL on
  its rerun.
- **Passes, in the bits.** One pass more or fewer moved no bit at the
  samples in 7 of 7 tries at binary64; the trailer's mean_pc_iterations,
  to 3 decimals, sees it, and on a first attempt max_pc_iterations too.
  COMPUTED: past about 1,000 attempts one pass is not certain to show in
  the mean, which a ring wave reaches past about 28 orbits (V3).

The driver prints that limit once, beside its plan, whenever it reads a
producer's files, and its gate holds the line word for word. Every
member with a rejected step stays in the member sample, so a ring
started from too large a first step is rerun whole, files or not.

**Also not built:** sampling by a close approach to a moon, which ROUND4.md
named and which has no threshold; steps and passes an orbit, since a
record carries no period and a wave's step is its shortest orbit's; both
format pairs in one call.

**The GPU seam test** (ROUND4.md, "What the integrator keeps"): built in
round 5 as part 3 of check-gpu, below. The driver holds a `hip:0` record
of the two-wave ring to this tree's CPU reruns, and one digit is planted
a run in each wave, in both forms, in every class a particle's replay
carries, each caught by line and column.

## What check-gpu holds of a GPU run of this tree

Round 5's P1 (VALIDATION entry 55): this tree as it is, its arithmetic
through libcft's `hip:N` device at the pin. READ in cft-fp256's
docs/HOSTAPI.md, "The GPU backend", at `eae74e3e`: an element is
`cft_sf_compute` and a program lane is `seq_run_block`, one thread a
lane; the device refuses an index table, augmentedAddition's two
control codes, and any non-zero scratch depth at open, 256 included,
since it fixes its own at 256 slots.

    make clean
    make HIP=1 [HIP_PLATFORM=nvidia HIP_ROOT=... HIP_ARCH=...] all programs
    make HIP=1 [the same variables as the build] [GPU_EXPECT=TEXT] \
         [GPU_LEGS=a,b | GPU_SKIP=a,b] check-gpu

It is in none of the three tiers, does not use the gate cache, and
refuses by name where there is no `hip:0`. It runs eleven legs, ten in
their quick form and check_dropin in its `--light` form: check_records,
check_equivalence, check_program_engine, check_ensemble, check_ring,
check_replay, check_compare_run, check_zonal, check_reversal,
check_corrector and check_dropin_light.

**What it holds,** for every process a leg starts, whatever starts it:
- **That it ran on the device.** Every process prints the probe's
  `CFT_HIP_INFO` device line, or it is a stand-in and its leg FAILs. A
  process refused by its own arguments before it opens a device, as the
  software backend refuses the same arguments, is the exception: it is
  counted apart (tools/check_gpu.py:565-570), and still as a device
  process (:562). The record's `backend=hip:0` only echoes the artifact
  string (entry 54).
- **That its record is the software backend's.** The same arguments run
  again on the software backend, and stdout and every output file are
  compared byte for byte, with the header's `backend=` value, the
  self-tests' banners, the trailer's 13 timings by name and an
  unpublished route's `by=` set aside. MEASURED on both GPUs: every
  process of every leg, 774 of them in check_compare_run alone.
- **That nothing ran past the witness,** where strace is installed:
  every exec accounted for, on nvidia-box's RX 7600.
- **The routes:** forced on, `scalar` and `reduce_seg` give the software
  backend's bits, and `indexed` and `pair_gather` are refused by name.
- **The GPU seam:** compare_run.py holds a `hip:0` record of the
  two-wave ring to this tree's CPU reruns. A plant in each wave of every
  class a particle's replay carries, the schedule, the trailer's
  counters, a massive body and a sampled particle's own body, and in
  REBOUND's form of a body no replay carries, fails as its class
  requires. P1's two member mutants of compare_run.py, and V1's four (a
  particle never compared, its own body skipped, a member's or a
  particle's counters never compared), each FAIL it.

**What runs where, measured at the pin** (entry 55):
- the clock's augmented sum on the host, at every step attempt;
- the divides and square roots by the program route, their rounding on
  the host;
- no reduction on the device by default, the corrector's maxima being
  the host's;
- every program within the device's 256 scratch slots and 512
  constants, the widest with no margin.

**What it does not hold.**
- **That every operation ran on the GPU.** The device line proves the
  process opened it (V1).
- **A bypass on the RTX 5060 Ti.** WSL cft-gpu has no strace.
- **stderr,** where a difference is printed and never failed on.
- **A refused process as other than a device process.** A process
  refused by its own arguments counts as one, and its record among the
  identical ones, so a leg whose processes were all refused passes the
  driver's check that a process ran on the device; only the leg's own
  verdict would fail it (V1; the second of the merge message's three
  open notes).
- **The legs it does not run.** gate_real's self-contained mode refuses
  `hip:0` on its identity check, and check_seams fails there on the
  forced routes the device does not publish. check_checkpoint passes on
  `hip:0`, outside check-gpu.

LEAD-FILL: the wave kernel's section (round 5, P2, entry 56): what check_wave holds of a kernel's records and files

## The format difference's scope

`--format-diff F2` reruns a member with `--dt-out`, then at F2 with
`--dt-file` and `--epsilon 0`, so on the same accepted steps, and
measures the two with `tools/compare_formats.py`. **It shows the narrower
format's round-off on those steps, and a narrowing both runs share passes
it.** That is the lead's ruling of 2026-11-12, and the driver prints it
beside every number.
- **MEASURED by P1.** The zonal term rounded to binary64 in both runs
  leaves binary128 against binary256 at 3.371e-33 of the distance from
  Saturn over 120 steps; the correct build gives 3.402e-33. Rounded at
  binary128 only: 4.029e-18.
- **What sees a shared narrowing** is a reference without it:
  cft-fp256's golden model and mpmath, which fail that build in
  check_zonal's gates 1 and 2, or another build, such as REBOUND's own
  IAS15, at the one width it has.
- **So on the GPU,** a wave's binary64 on its binary128 steps measures
  binary64's round-off, and binary256 on a tenth of the waves measures
  binary128's. Neither can show that the GPU's arithmetic is the
  format's. The bit comparison against this tree can.
- **GPU-E's section 6, item 3(b)** proposes the format difference as "the
  check that fails on exactly the defect the README names", and its
  question 4 leans on it. BELIEVED, by message from cft-fp256's lead on
  2026-11-12: both are to be corrected in a new commit on its main. At
  cft-fp256's `be40005b` the study is unchanged since `56006987` (MEASURED by
  the lead, read only: `git diff 56006987 be40005b` of the study is empty).
- **A rejection,** beside the number. The replay at F2 takes no rejected
  attempt. A rejection before a member's first accepted step leaves
  state, the rejected attempt's corrected b, so that member's difference
  includes it; one after an accepted step leaves no state bit (entry 49).
  The note beside each number says which.

## The numbers GPU-E's cost model needs, from this tree

| quantity | this tree | GPU-E | source |
|---|---|---|---|
| waves for 10^6 particles | 17,544 full waves of 57: 1,000,008 particles | 17,544, its last of 49 (BELIEVED, by message) | P2; entry 48 |
| pairs a force evaluation, a wave | 420 under `active 7`; 2,016 with every body active | 420 | P2, read from the record's `# active` line |
| steps an orbit, a wave inside Pan's orbit | 35.89, for the wave with radii to 132,941 km, at epsilon 1e-9 | 36.2, MEASURED on a circular orbit about a point mass | P2; -0.9% |
| a wave outside Pan's orbit | steps at Pan's rate, 1,388.09 s: Pan's period, 49,819 s, over 35.89; 38.4 steps per orbit of the outermost wave's innermost particle | not in the model | P2, a 228-particle ring at binary64 |
| more steps over the whole ring | about 3% | not in the model | COMPUTED by P2 from a uniform radius, not run |
| J2 and J4, a body a force evaluation | 3 sub, 34 mul, 14 add, 1 sqrt, 10 div; no FMA; two of the divides by 2 and 8, exact | ESTIMATE: 1 DIV, 21 MUL, 6 FMA, 6 ADD | V1, counted by hand, equal to the trace; entry 47 |

- **The 34 multiplies are every body's,** test particles included: the
  13 that form the reaction on Saturn are as wide as the bodies, and a
  test particle's are thrown away (V1).
- **Beside them, a call and not a body,** 8 multiplies forming 3 J2 R^2
  and 5 J4 R^4, and one subtract a reacting body: P1's count from the
  trace, which V1's per-body count does not include.
- **The order is REBOUNDx's as recalled,** BELIEVED. A one-divide order
  is possible; it would give up decision 4's operation-for-operation
  comparison with REBOUNDx (P1).
- **Under `active 7`** the term acts on the 63 other bodies of a wave, and
  Saturn takes the reaction of the 6 moons.
- **A force evaluation is gravity, then the term;** a step attempt of p
  passes makes 7p + 1 of them (GPU-E section 2, READ; unchanged here).

## What is BELIEVED

- **R = 60,330 km** as the radius Cassini's J2 and J4 refer to. GPU-E's
  table names none.
- **REBOUNDx's order** in `src/ias15_zonal.inc`, recalled and not read.
  REBOUNDx is not pinned (decision 4). The formula and the physics gates
  hold the term either way.
  READ by study DEBRIS-A (2026-11-12, from REBOUNDx's repository on the
  web, nothing run against REBOUNDx): the recalled order is REBOUNDx's as
  of 5406ae5 (2020-04-13). REBOUNDx rewrote gravitational_harmonics.c on
  2024-10-20 (8402022, a general spin axis), with another order. So
  decision 4's comparison would match the 2020 order only, and a current
  REBOUNDx is BELIEVED to differ from the port in bits. The study did
  MEASURE the port at binary64 bit-identical to REBOUND 5.1.1's IAS15
  with J2 written in Python in the 2020 order (low Earth orbit, J4 zero).
- **The same ring file on every host.** Measured on one host, with
  mpmath's two backends.
- ~~**The exact-time pair from `vaugadd`,** and G1B refusing AUGADD on
  hip:N, so a run on `hip:0` as it is would very likely meet it.~~
  MEASURED false in round 5 (VALIDATION entry 54). `vaugadd` is the API
  call `cft_augmented_add`, which runs on the host on every device, and
  every `hip:0` run carries its exact-time pair. AUGADD is refused only
  inside a program, and none of this port's programs uses it.
- **About 3% more steps from Pan's cap,** computed from a uniform radius.
- **Why the two forms' steps differ** at 1.7e-11: PRS23 builds the step
  from the high b coefficients, which amplify round-off (P2).
- **The POSIX pipe branch** of `check_replay`, `/dev/stdin`, which has not
  run here.
- **What cft-fp256's lead said by message:** G1B's merge, GPU-E's coming
  correction, and its cost model's last wave of 49.

## The decisions that are Logan's

ROUND4.md's list, as the round left it, and Logan's answers of
2026-11-13, relayed by cft-fp256's lead in his words:

| | decision | as the round left it | Logan, 2026-11-13 |
|---|---|---|---|
| 1 | merges and pushes for this round | delegated to the lead on 2026-11-11, on round 3's bar, for round 4's commits | unchanged |
| 2 | the study's questions 1, 3 and 5: the workload, the form and the formats | open; the CPU half serves every answer | open |
| 3 | its question 2, a device divide before cft-fp256's exp family | open; it decides when the GPU half's step 3 starts. The term adds 10 divides and a square root a body a force evaluation to gravity's one each a pair | "Yes with GPU round 2": cft-fp256's, expected the week of 2026-11-16 |
| 4 | REBOUNDx as a pinned upstream | open; one operation-for-operation comparison at binary64 would confirm or correct the order. Study DEBRIS-A read that the recalled order is REBOUNDx's of 2020 (5406ae5), rewritten in 2024 (8402022): the comparison needs that older pin | "Keep default": no pin |
| 5 | the moons' phases | synthetic; an epoch's needs JPL's tables from the web | open |
| 6 | P5, the wave kernel at the hardware's binary64 | open, and then a GPU slot from cft-fp256 | "Yes"; its GPU slot is to be booked in the TimeSheet, early (none yet) |
| 7 | the pin's move to G1B | open, once G1B merges | "Yes, to G1C": `eae74e3`, VALIDATION entry 53 |

GPU-E's question 4 asked for J2 and J4 "verified by the golden model and
the format-difference check". The round answers with the golden model,
mpmath and the physics; the format difference cannot see the defect it
was meant for. Logan answered "Yes" on 2026-11-13 (relayed): on the GPU,
J2 and J4 are held as tools/check_zonal.py holds them here, with the
format difference kept for a narrowing in one format.

## The record

docs/VALIDATION.md:
- **46:** P0, the `active` and `zonal` lines and the seam;
- **47:** P1, J2 and J4 in libcft's format, and the format difference's
  blind spot;
- **48:** P2, the ring generator, Pan's cap and 1,000,008 particles;
- **49:** P3, a particle replayed out of its wave, and which rejections
  leave state;
- **50:** P4, the driver;
- **51:** the seams, the P3 x P4 wiring, VL, the merges and the suites;
- **54:** round 5's P0, this tree built with libcft's GPU backend;
- **55:** round 5's P1, this tree on `hip:0` as it is, and check-gpu;
- **57:** round 5's P3, the driver taking a producer's own files.
