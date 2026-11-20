<!-- SPDX-License-Identifier: GPL-3.0-or-later
     Copyright 2026 the cft-rebound contributors. -->

# Round 5: the GPU half begins - this tree on a GPU as it is, then the wave kernel at the hardware's binary64

**Status, 2026-11-20: all four parcels are on main, and the round closes
when these records land.** P3 landed at `4a8fe99` and P1 at `fda3b40` on
2026-11-18, PL, the lead's own gate fixes, at `620c98b` on 2026-11-19,
and P2 at `61e45bc` on 2026-11-20, each on its verifier's READY and a
full make check (docs/VALIDATION.md entries 55 to 58; P0 is entry 54).
The pin stays at cft-fp256 `eae74e3`; its move waits on Logan's word.

Logan approved the round on 2026-11-16 ("the plan is also approved"),
after its outline: the pin to G1C, this tree on `hip:0`, then P5 at
binary64. His answers of 2026-11-13, relayed by cft-fp256's lead
(ROADMAP.md, "The decisions that are Logan's"):
- decision 6, P5: "Yes";
- decision 7, the pin: "Yes, to G1C", done in docs/VALIDATION.md entry
  53;
- GPU-E's question 4: "Yes", J2 and J4 held on the GPU as
  tools/check_zonal.py holds them on the CPU;
- decision 3, a device divide: "Yes with GPU round 2", which is
  cft-fp256's work, not this round's.

Round 4 built on the CPU what a GPU wave will be held to: the ring
generator, tools/compare_run.py, the particle replay and J2 and J4 in
libcft's format (docs/RING.md). This round takes the first two GPU steps
of ROADMAP.md's "The GPU half": this tree as it stands, run on a GPU
through libcft's new `hip:N` device, and the wave kernel at binary64,
the hardware's own.

## Where it starts

**The pin.** cft-fp256 `eae74e3`, its main with G1B, the HIP backend,
and G1C, its GPU gates (entry 53). The same calls return the same bits
through it on the software backend.

**What `hip:N` does with this port's calls.** READ, from cft-fp256's
lead on 2026-11-16 and its docs/HOSTAPI.md, "The GPU backend":
- **On the GPU:**
  - elementwise `cft_run`, one `cft_sf_compute` an element;
  - `cft_reduce`, softfloat.c's tree, one level a launch. This port's
    default path never calls it: the corrector's maxima are `cft_class`
    and a host loop unless `reduce_seg` is forced on (src/ias15_cft.c:1835).
    Under gdb on the RX 7600, P1 found no backend reduction call, and V1
    no entry to `cft_reduce` or `cft_reduce_seg`;
  - `cft_program_run_ex`, one thread a lane. Scratch depth is fixed at
    256 slots and banks at 512. This port's widest program,
    gravity-fold256-fp128, uses all 256 slots: it fits, with no margin
    (P1, measured on both GPUs). A program using augmentedAddition's codes
    AUGADD and AUGERR is refused by its loader, by name. An index table is
    refused at the call (device.c:861 and :2281; V0). `seq_features` is
    0x75f1f: the software backend's less INDEXED and AUGADD.
- **Partly on the host:** `cft_div` and `cft_sqrt` take the program
  route. One sequencer program a chunk runs the floating core on the
  GPU; classification, centring and rounding run on the host
  (divsqrt.c's `divsqrt_via_program`).
- **On the host, on every device:** `cft_augmented_add`, which this
  port's exact-time clock calls at every step attempt, rejected ones
  included (src/ias15_cft.c:3049; P1 and V1 under gdb, 5 calls in 5
  attempts), and its compensated sums under `--cs augmented` (:1680).
- **Refused:** `cft_get_image_id`, since a GPU loads no image; and a
  non-zero `cft_open_ex` depth, which this port never asks for (it calls
  `cft_open`).
- **G1C** held `hip:N` to the software backend's bits on both GPUs at
  every format: 14 corpora, 125,032 vector cases, and a replay of
  1,228,415 cases.
- **A call's cost.** G1B's first part measured about 87 µs a call on the
  RTX 5060 Ti and 234 µs on the RX 7600, staged. G1C measured 107 µs on
  the 5060 Ti. V0 measured this port's own on the RX 7600, `wall_seconds`
  over the trailer's calls: 286 µs on kepler's loop engine at binary64,
  and 218 to 769 µs over its 24 runs (binary64 218 to 364, binary128 336
  to 472, binary256 364 to 769).

**P0's measurements, 2026-11-16.**
- `make HIP=1 all programs` built this tree on nvidia-box in 31 s, and
  ias15_cft links ROCm's libamdhip64.
- One binary, the software backend against `hip:0`, on kepler at
  binary64, binary128 and binary256, on the loop engine and on the
  program engine in the FMA form. Every sample line and every counter is
  the same in all 6 runs: 50 steps each, 3.44, 8.52 or 8.58, and 12.00
  passes a step. The only differing lines are the routes: `indexed` and
  `pair_gather` read `by=absent` on `hip:0` against `by=opt_in` in
  software, since the device publishes no INDEXED.
- Wall seconds, one run each (MEASURED, not ranked):

  | format | engine | software | hip:0 |
  |---|---|---|---|
  | binary64 | loop | 1.10 | 31.6 |
  | binary64 | program | 0.85 | 7.98 |
  | binary128 | loop | 3.29 | 85.9 |
  | binary128 | program | 2.40 | 18.8 |
  | binary256 | loop | 7.13 | 154.0 |
  | binary256 | program | 3.78 | 35.1 |

  So this port's default path is bound by its calls on a GPU, as on the
  card (entry 37). V0 measured 218 to 364 µs a call at binary64 from 2
  to 64 bodies, whatever the size; P2's waves are larger, and nothing
  above 64 bodies is measured. How that compares with software depends on the
  problem. V0 found a 64-body ring wave faster on `hip:0` than in
  software: 6.31 against 7.19 s on the loop engine, and 4.85 against
  7.26 s on the program engine, one run each. Its GPU time is still its
  calls (244 and 271 µs a call). What changes at 64 bodies is that
  software is slower.
- `CFT_HIP_INFO=1` names the device: "AMD Radeon RX 7600 (gfx1102, 16
  multiprocessors, bus 0000:03:00.0, 8176 MiB), HIP platform amd,
  runtime 71626385".

**Study GPU-E2** (cft-fp256's docs/studies/GPU-E-mi300x-ring.md, READ):
- **The wave.** Its shape a2 is one 64-lane wave a system: Saturn with J2
  and J4 and six moons in seven lanes, and 57 test particles. Each wave
  takes one step and one pass count, and holds its maxima across the
  wave. The study runs one kernel a substep, each thread a body.
- **What it leaves open.** It states no mapping onto NVIDIA's 32-wide
  warps or RDNA3's wave32.
- **Hardware binary64.** Every binary64 figure in it is libcft's own,
  and it says nothing about the hardware's binary64. cft-fp256's
  ROADMAP believes the native FMA correctly rounded, but notes that
  libcft's contract also needs flags and subnormals.
- **Its estimates for the ring at binary64,** 10 orbits:
  - the RX 7600, 3 x 10^4 particles: 13 minutes;
  - the RTX 5060 Ti, 10^5 particles: 14 minutes.

**cft-fp256's GPU round 2,** READ from its lead: augmentedAddition on
the device, index tables and deeper scratch, and the in-program divide
and square root. It is drafted, dispatched from Tuesday 2026-11-17 at
the earliest, and shares both GPUs.

**The GPUs:**
- **The RX 7600:** nvidia-box (192.168.0.2, by SSH), ROCm with HIP 7.16
  and hipcc for gfx1102. This round works in ~/cft-rebound-r5* there,
  and never in another session's directories.
- **The RTX 5060 Ti:** WSL `cft-gpu`, through HIP 7.2.4's NVIDIA layer
  (`HIP_PLATFORM=nvidia HIP_ROOT=/opt/rocm-7.2.4 HIP_ARCH=120`).
  - It is CUDA device 1 in PCI order. Set `CUDA_DEVICE_ORDER=PCI_BUS_ID
    CUDA_VISIBLE_DEVICES=1`, so that `hip:0` is the 5060 Ti and never
    the 4070.
  - The RTX 4070 drives the displays and stays unused.
- **Booking:** every run of 10 minutes or more is booked in the
  TimeSheet (`C:\Users\logan\source\repos\TimeSheet`), with its mode
  stated, and released when it ends.

## P0 - a HIP build of this tree, and a GPU host it runs on *(the lead's)*

- **Done:** the Makefile's `HIP=1` (`d95ac0d`). It passes `HIP_ROOT`,
  `HIP_PLATFORM`, `HIP_CUDA` and `HIP_ARCH` to libcft's Makefile and
  links the GPU runtime: libamdhip64 on ROCm, the CUDA runtime under
  HIP's NVIDIA layer. `HIP=1` with `XRT=1` is refused, as libcft refuses
  it.
- **Done:** ~/cft-rebound-r5 on nvidia-box, cloned from git bundles of
  this repository and its two pins, carried over SSH. Nothing was
  fetched from the web, and check-pins passes there.
- **Done:** the smoke above.
- **Done:** VALIDATION entry 54, the round's ledger, and P0's full
  `make check` on the desktop at `d7b5047`, where `HIP=0` changes
  nothing: rc 0, 25 legs, 1,463 lines beginning `  ok `, 0 FAIL, 06:06
  to 09:38 on 2026-11-16, on a shared host (r5-check-d7b5047.log). main
  fast-forwarded to `5b1cac6` and was pushed that morning. P1 built the
  tree in WSL `cft-gpu` (VALIDATION entry 55).
- **Found, not fixed:**
  - GCC 13 on Ubuntu warns at src/ias15_cft.c:2155 that `vcopy`'s
    source is null. It does so only in the library object
    (`-DIAS15_CFT_LIBRARY`), where `load_programs()`, which holds
    `prog_dep`'s only allocation, is compiled out. That makes `prog_dep`
    null there outright. It is harmless only because the library sets
    `engine_program = 0` and nothing sets it to 1, so the call is never
    reached; a library path that turned the program engine on would copy
    from null. The program's own object does not warn (V0).
  - `HIP_PLATFORM=cuda`, libcft's plain-nvcc platform, falls into the
    Makefile's amd branch. `make -n` shows `-lamdhip64` on its link, so
    it is expected to fail there; not run, since nvidia-box has no CUDA.
    `make HIP=1 example` fails (run): examples/Makefile has no HIP block.
    Neither is this round's path (V0).

## The parcels

### P1 - this tree on `hip:0` as it is *(ROADMAP's GPU half, step 2)*

Every record the suite makes, made on `hip:0`, must be the software
backend's, or its leg refuses by name. Beside the bits go what the
default path costs on a GPU, so that P2 knows what it replaces.

Measured by P1: a whole-output comparison sets aside the header's
`backend=` value, the trailer's timings and the `by=` of a route the
device does not publish, and also the first lines of the door and
corrector self-tests, which name the backend they ran on.

- **Where:** the RX 7600 first, then the RTX 5060 Ti.
- **What:** the suite's record-making legs in `--quick`, with
  `CFT_REBOUND_ARTIFACT=hip:0`:
  - check_records, check_equivalence, check_program_engine at binary64;
  - check_ensemble, check_ring, check_replay and check_compare_run;
  - check_zonal and check_reversal;
  - check_corrector and check_dropin. check_dropin goes through the
    registered integrator; check_corrector does only in its `--library`
    step, and otherwise runs the standalone program (P1's report, item 3;
    this line first said both went through it).

  A full `make check` is out of reach. On kepler the calls make a run 7
  to 30 times slower, and V0 measured `check_dropin --light` at 449 s in
  software on nvidia-box, so roughly 1 to 3.5 hours on `hip:0`. Book
  such windows long, or hold a named subset. Measured: inside
  check_corrector's `--library` step, `check_dropin --light` took 62 min
  alone on the RTX 5060 Ti (P1, at `20dad58`, 2026-11-16) and 66 min
  there (V1, at `ab75fbe`, 2026-11-18), each time with its host, the
  desktop, shared.
- **Prove the GPU ran.** A leg that quietly ran the software backend is
  a stand-in, and is a finding. Each run's `CFT_HIP_INFO` device line is
  required beside its verdict. The record's `backend=hip:0` only echoes
  the artifact string the run was given (src/ias15_cft.c:4697-4699). It
  does not prove the GPU ran; the device line does (V0).
- **The routes on `hip:0`** (V0, at binary64 and binary128): forced on,
  `scalar` and `reduce_seg` run on the GPU in the software backend's
  bits, and `indexed` and `pair_gather` are refused by name.
  `--door-selftest` on `hip:0` at binary64 holds 116 cases identical,
  with 48 not compared because their routes are absent. P1 holds them
  at every format.
- **The GPU seam test** (ROUND4.md, "What the integrator keeps"):
  tools/compare_run.py pointed at P1's first `hip:0` ring record, with a
  difference planted in that record and shown caught by name.
- **What runs where, confirmed by measurement:** the clock's augmented
  sum on the host at every step; the divides and square roots by the
  program route; the index-table routes absent; programs within 256
  scratch slots and 512 bank slots. P1 measured the first at every step
  attempt, rejected ones included, and the last with no margin: one
  program uses all 256 slots (VALIDATION entry 55).
- **Costs:** per call and per step, at each format, on each GPU,
  MEASURED and labelled, on the trailer's `wall_seconds`. Its `seconds`
  and `steps_per_s` are CPU time, 1.5 to 2.0 times the wall time on
  `hip:0` (V0, over 24 runs; BELIEVED, not measured: HIP's own threads
  are the cause). That ratio is ROCm's: P1 measured 1.54 to 1.92 on the
  RX 7600, and 1.02 to 1.03 on the RTX 5060 Ti through HIP's NVIDIA
  layer, in all 12 of its cost rows there. One run a row ranks nothing
  within the noise (CLAUDE.md).
- **Physics counters beside the bits,** as every gate prints them.
- **Not P1's:** any change to the port's source. A refusal that would
  need one is reported first.

### P2 - the wave kernel at the hardware's binary64 *(decision 6, P5)*

One system a 64-thread workgroup, in the hardware's own binary64. It
writes this repository's record format, so that tools/compare_run.py
holds it to the CPU port bit for bit.

- **The shape,** GPU-E2's a2: Saturn with J2 and J4 and six moons in
  seven lanes, and 57 test particles. One step and one pass count a
  system. The system's maxima and minimum are reduced in shared memory,
  so that the kernel is the same on NVIDIA's 32-wide warps and on
  RDNA3's wave32 and wave64; a wave-level intrinsic is an optimisation
  for later.
- **The arithmetic is the hardware's binary64:** add, multiply, fma,
  divide and square root. They are BELIEVED correctly rounded, with
  subnormals, on both toolchains by default (no fast-math). P2's first
  gate measures that on a corpus against libcft's software binary64,
  `cft_div` and `cft_sqrt` included, before anything rests on it.
  - **Contraction, measured by P2 (2026-11-16).** The control P2's brief
    named, the checked source built with contraction on, cannot fail as
    worded: with the flag checks in, every product is read by its check,
    and that build gave the shipped build's records on g++ and on nvcc.
    The control that bites compiles the checks out: with contraction off
    that source passes compare_run, and with it on it fails.
- **Flags.** A GPU has no IEEE flags. Where the port aborts on OVERFLOW,
  DIVBYZERO or UNDERFLOW, the kernel checks the results itself and
  stops by name. Anything it cannot check is a stated limit.
- **The form is the port's FMA form at binary64** (`--arith fma`),
  operation for operation:
  - the predictor and the corrector's passes;
  - gravity in ascending pair order, with the `active` split;
  - J2 and J4 in src/ias15_zonal.inc's order;
  - the step control;
  - the exact-time clock, an augmented sum as libcft's
    `cft_augmented_add` computes it. Not TwoSum in binary64, as this line
    first said: libcft's augmentedAddition rounds an exact tie toward
    zero, where TwoSum rounds it to even. P2 measured 7,259 such ties in
    152,735 cases, and a TwoSum clock's records differ from the port's in
    two of the three cases it tried (2026-11-16).

  So its records are `ias15_cft --arith fma --format fp64`'s, bit for
  bit, every column, member by member.
- **How it is built:**
  - The step's arithmetic is `__host__ __device__` functions over
    `double`. They run first on the host, with no GPU slot, held to the
    CPU port's records.
  - Then the same source runs as a HIP kernel: hipcc for gfx1102 on
    nvidia-box, and HIP's NVIDIA layer for sm_120 in cft-gpu.
  - A `make HIP=1` target builds it.
- **What it is held to:**
  - tools/compare_run.py against CPU reruns, on the ring generator's
    waves with their zonal line: each member, and particles through
    round 4's replay;
  - the GPU seam test on its first record;
  - J2 and J4 at binary64 to the CPU port's bits, which check_zonal
    holds to the golden model and mpmath (Logan's "Yes" to GPU-E's
    question 4).
- **What it measures:** system-steps a second on each GPU at a few wave
  counts, MEASURED and labelled, set beside GPU-E2's estimates. No other
  speed claim.
- **Not P2's:** binary128 and binary256, which need the device divide
  (cft-fp256's round 2); more than one GPU; the ladder's full sizes; the
  MI300X.

### P3 - the driver takes a producer's own attempts records and schedules *(CPU only)*

docs/RING.md: "To sample its particles without a rerun, a GPU run must
write each wave's attempts record and corrector schedule, and the driver
must take them: not built."

- **What it builds:** tools/compare_run.py learns to take, beside a
  producer's record, each member's attempts record and corrector
  schedule, in the formats ias15_cft's `--attempts-out` and
  `--record-iters` write. It samples particles by replaying them against
  the producer's own schedule, with no rerun of the member.
  - **Measured by P3: a producer's files hold its arithmetic, not its
    decisions.** A replay takes each attempt's dt and outcome, and its
    passes, as given. A producer whose first attempt took 5 passes where
    this tree's takes 7, so that 2 of its 3 sample lines are not this
    tree's, has that member held, NOT RERUN, both its particles
    IDENTICAL (P3's GATE 1b, whose run ends rc 1 because another member,
    rerun, FAILs; VALIDATION entry 57). Where nothing else differs, the
    run says RESULT: PASS (V3, on a three-wave ring, wave 0 at 2 passes
    where this tree takes 4). The same member FAILs once it is rerun. The
    decisions are held only in a member that is rerun.
- **What it is held to:**
  - with ias15_cft as the producer, particles sampled from its own files
    are identical to those sampled through a rerun. This test cannot fail
    by itself: ias15_cft's files are a rerun's to the byte, and a driver
    that replays rerun members from the rerun's schedule instead, while
    naming the producer's, passed it (P3's mutant M1). GATE 1b's producer,
    whose exit test is not this tree's, fails such a driver;
  - one change planted in each file, the attempts record and the
    schedule, is caught by name;
  - a schedule from another wave is refused by its run key.
- **No GPU.** P2's kernel writes these files when it can, and the seam
  test below joins them.

## Rules for every parcel

Round 4's rules (ROUND4.md, "Rules for every parcel") hold, except
where this round's GPUs change them:
- **The GPU hosts.** A parcel may use a GPU only in a window it has
  booked in the TimeSheet, as `--who cft-REBOUND/<parcel>`, with its mode
  stated (`exclusive` for gates and stage runs), and must release it.
  - The TimeSheet is on the desktop, so a run is booked there before
    its SSH session, with that session's Windows pid attached (`book ...
    --pid`).
  - It runs only in its own directory there, ~/cft-rebound-r5-<parcel>
    on nvidia-box and /root/cft-rebound-r5-<parcel> in cft-gpu.
  - It never touches another session's directories.
  - One agent at a time on a GPU host.
- **Pins come by bundle.** Every clone, here and on a GPU host, is made
  from this repository's local clones, and nothing is fetched from the
  web. P0's bundles are in the round's ledger.
- **The device trap.** In cft-gpu, `CUDA_DEVICE_ORDER=PCI_BUS_ID
  CUDA_VISIBLE_DEVICES=1` comes before any `hip:0`. Keep `CFT_HIP_INFO`'s
  line in every record set.
- **The stale-library trap** (V0). In a tree first built without HIP,
  `make HIP=1 all` returns 0 but leaves libcft.a software-only, and
  `hip:0` is then refused by name. `make clean` first.
- **No speed claims** beyond MEASURED rows, labelled, one run each,
  unranked within the noise.
- **No card, and never the RTX 4070.**

## What the integrator keeps

- **P0,** before anything is dispatched.
- **The seam tests:**
  - **The GPU seam,** in P1's brief, which is the GPU half's first: a
    difference planted in the first GPU record, shown caught.
  - **P1 x P2:** one ring wave at binary64 in the FMA form three ways,
    on the CPU, on `hip:0` through libcft, and in P2's kernel. All three
    are the same bits, beside the physics counters.
  - **P2 x P3:** the driver samples P2's particles from P2's own
    attempts records and schedules, if P2 writes them.
- **The integrator-only files,** as in round 4, and this file.
- **Every merge serially,** each with a full `make check` on the desktop
  after it, never while a suite runs. The GPU legs are not in the
  desktop's tiers. Each parcel's GPU results are read from its logs at
  the merge.
- **A verifier per parcel,** V1 to V3, re-running from a clean build
  and fixing nothing. A GPU verifier books its own window.
- **VALIDATION entries** 54 (P0) onward.

## Order of the week

| when | who | what |
|---|---|---|
| Monday | lead | P0: done on the RX 7600 except the launcher, the entry and its suite. P3 dispatched (CPU). P2's host half dispatched (CPU). |
| Monday afternoon | P1 | On the RX 7600, in the lead's window (09:00 to 17:00), or its own. |
| from Tuesday | P1, P2 | The 5060 Ti and the RX 7600, in windows clear of cft-fp256's round 2. |
| as each returns | V1 to V3 | Verify. |
| as each lands | lead | Merge, seam test, full `make check`. |

## Decisions that are Logan's

1. **Merges and pushes for this round.** Delegated as in round 4, or
   each at his word.
2. **mpmath on the GPU hosts.** nvidia-box's Python has none, and the
   gates that check exactly need it.
   - (a) Copy the desktop's mpmath 1.3.0, which is pure Python, into a
     directory beside the clone, on PYTHONPATH.
   - (b) pip-install it, which is a download.

   The plan takes (a) unless he says otherwise.
3. **P2's flags at the hardware's binary64.** By default the kernel
   checks its own results for what the port's flags catch, and states
   whatever it cannot.

## Not in this round, and why

- **Binary128 and binary256 on the GPU,** in a kernel. They need a
  divide and square root inside it: cft-fp256's GPU round 2.
- **The ladder at full size** (3 x 10^4 and 10^5 particles) **and the
  MI300X:** after P2's kernel is held at small sizes. The MI300X waits on
  AMD's time.
- **The card:** nothing in this round needs it.

## What the round found against this plan

Each item was measured unless it says otherwise. Its record is in
docs/VALIDATION.md, by entry.

**Where it starts, and P1** (entries 54 and 55):
- `cft_reduce` runs on the GPU, but this port's default path never calls
  it.
- The 256 scratch slots fit with no margin: gravity-fold256-fp128 uses
  all of them.
- The clock's augmented sum runs on the host at every step attempt,
  rejected ones included, not at every step.
- The trailer's CPU time is 1.5 to 2.0 times its wall time under ROCm,
  and 1.02 to 1.03 under HIP's NVIDIA layer.
- A whole-output comparison must also set aside the door and corrector
  self-tests' first lines, which name the backend.
- check_corrector goes through the registered integrator only in its
  `--library` step.
- `check_dropin --light` took 62 and 66 min alone on the RTX 5060 Ti,
  inside the plan's estimate of 1 to 3.5 hours.
- Of the legs outside P1's list, gate_real's self-contained mode refuses
  `hip:0` on its identity check, check_checkpoint passes there as it
  is, and check_seams fails only on the forced routes the device does
  not publish. P1 first read gate_real's refusal as binding the other
  two.
- The GPU seam as worded, one difference planted and shown caught, was
  not enough. P1's seam went back three times, each time for a copy of
  compare_run that compared less and passed.

**P2** (entry 56):
- The clock is not TwoSum in binary64: libcft rounds an exact tie toward
  zero.
- The contraction control as worded cannot fail on the checked source.
  The control that bites compiles the flag checks out.
- The flags the port aborts on include INVALID, which the plan's list
  leaves out.
- Under nvcc, a bit-level negation or abs is an add, which quiets a
  signaling NaN.
- Gate 4's corpus depended on the compiler's order of argument draws.
- MinGW's library fma and sqrt are not correctly rounded, and g++ calls
  them at -O0.
- GPU-E2's estimates priced libcft's software binary64, not the
  hardware's, so they are not the kernel's.
- The port computes a ring wave's test-test pairs and the kernel does
  not, so two coincident test particles stop only the port.
- check_wave went back six times for copies of compare_run, or of
  check_wave, that did less and passed, and once for its work directory,
  which was the kernel's own path (and entry 58).
- The plan left the ladder's full sizes out of the round. P2's rate rows
  ran GPU-E2's rungs over 10 orbits anyway, 527 waves on the RX 7600 and
  1,755 on the 5060 Ti, one run each, with two members of each held
  IDENTICAL by compare_run. The rungs were not compared with each other.

**P3** (entry 57):
- A producer's files hold its arithmetic, not its decisions; those are
  held only in a member that is rerun.
- The plan's first P3 test cannot fail by itself. GATE 1b, a producer
  whose exit test is not this tree's, can.
- A planted pass count is caught by the trailer's counters, not by the
  bits.
- "With no rerun of the member" never applies to a rejecting wave, which
  stays in the member sample.
- compare_run's output at first said a NOT RERUN member was held, and
  passed a producer whose exit test stopped two passes early. The limit
  is now printed, and held word for word.

**The integration** (entry 58):
- The gate cache's keys missed files that legs execute or read: the
  program images in a fresh clone, check_ensemble's images, 6 inputs in 3
  legs, and the constants leg's REBOUND file. PL closed them.
- A GPU run with other arguments than the gate's is not the gate's run.
  check_wave's own line in check-gpu first ran at the merge, and crashed
  on its work directory; the lead's seam script first met a GPU record's
  header there too.
- cft-fp256's GPU round 2 landed during the round. Round 5 stayed on
  `eae74e3` at Logan's word, and the pin moves only on that word, after
  the close.

**Not in this round, as planned:** binary128 and binary256 on a GPU, the
ladder compared rung against rung, and the MI300X.
