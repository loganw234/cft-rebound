<!-- SPDX-License-Identifier: GPL-3.0-or-later
     Copyright 2026 the cft-rebound contributors. -->

# Round 4: the ring's CPU half - what a GPU wave will be held to, built here first

**Status, 2026-11-12: built on the CPU, and every commit READY by its
verifier. P0 is on main at `9a42b1d`, P1 at `291c99c`, P2 and P4 at
`7d50392`, and P3 with the seams and the wiring at `d9c8c8c`, each on the
suite docs/VALIDATION.md entry 51 records.
Entries 46 to 51 are the round's record, docs/RING.md is the ring's CPU
half as it stands, and "What the round found against this plan" is at
the end. The plan below is as it was written.**

> "When time comes to run the MI300X trials, I'd also like to do an IAS15
> run, but utilizing as many resources from the MI300X as possible, perhaps
> something 'simple' but highly parallel like hundreds of thousands of
> bodies in a ring around a world? Similar to Saturn." (Logan, 2026-11-11)

cft-fp256's study GPU-E answered that on paper
(its docs/studies/GPU-E-mi300x-ring.md, merged at `56006987`). The same
day its lead handed the work here: the port derives from REBOUND's IAS15,
so it lives under REBOUND's licence, in this repository (README.md:163-170,
and the study's question 6, answered).

**The study's recommendation, which is Logan's to decide:**
- a test-particle ring, one 64-lane wave a system: 57 massless particles,
  Saturn with J2 and J4, and six moons, each wave on its own adaptive step;
- 10^6 particles from 130,000 to 142,000 km, in 17,544 waves sorted by
  radius;
- binary128 at epsilon 1e-9 for 10 orbits, binary64 on the same steps,
  binary256 on a tenth of the waves as the difference reference;
- the FMA form first; about 5 to 13 device-hours on the generic tier (its
  ESTIMATE).

**Why this round runs on the CPU only.** A GPU wave can be held to this
port only if everything it computes exists here first, op for op. Three
things that run needs are missing here, and none of them needs a GPU:
- **J2 and J4 in libcft's format.** `additional_forces` is binary64 even
  in a binary256 run (README.md:122-128). The standalone program has no
  force but gravity: its `apply_force_hook` compiles to nothing
  (src/ias15_cft.c:2664).
- **Test particles in the standalone program.** The engine has them, and
  the drop-in sets them every step (src/reb_integrator_cft.c:943). But
  no line of a problem file sets them for `ias15_cft`.
- **A ring, and the records a GPU run will be compared against.**

Three more come from cft-fp256, and none has shipped:
- **G1B**, a GPU as a libcft device. It is on that repository's branches
  `g1b` and `g1b2`, due to merge on 2026-11-12 (per its lead).
- **A divide and square root callable inside a kernel.** No round plans
  one; that is the study's question 2.
- **Its gates through `hip:N`**, G1C.

So this round builds the CPU half. The GPU half (below, "After this round")
starts as those land.

**The target.** On the CPU, at binary64, binary128 and binary256, on
both engines and in both forms where they exist:
- **A Saturn ring integrates:** Saturn with J2 and J4 in the run's own
  format, six moons, and 57 test particles a wave computed as test
  particles, as one system or as an ensemble of waves.
- **J2 and J4 are held three ways:**
  - their arithmetic, operation by operation, against the golden model;
  - their formula, against the potential's own gradient in mpmath;
  - their physics: energy conserved, and apsides and nodes moving at
    theory's rates.
- **One particle replays out of its wave bit for bit,** rejected attempts
  included. That is the cheap check a GPU run's particles will get.
- **A driver compares a run's records with this tree.** It samples waves
  and particles, reruns them here, and says which match, with the physics
  counters beside the bits.
- **Every existing bit stays where it is:** no record line moves.

`METHOD.md` below is `../ParcelRound`'s, the canonical copy.

---

## What reading the code already turned up

Two read-only passes went over this tree and cft-fp256 on 2026-11-11
before a brief was written. Verify these before relying on them.

**In this tree:**
- **Test particles exist in the engine, and are gated.**
  `build_scatter()` (src/ias15_cft.c:1158-1184) builds REBOUND's two loops
  from `eng_active` and `eng_tptype`. A test particle's row holds its
  N_active partners in ascending order. An active body gets the reverse
  term only under `testparticle_type`. `ias15_cft` never sets
  `eng_active`: only `--door-selftest` does.
- **The study's massless claim holds and is gated, but not where it
  cites.** src/ias15_cft.c:1542-1546 is about `gravity_reference`'s +0
  padding. The claim itself is in tools/cases_pairs.c: a massless active
  body gets a test particle's bits, asserted by `case_massless_is_inert`.
  - It holds under round-to-nearest.
  - It fails where two massless bodies coincide: r = 0 gives inf times 0,
    a NaN.
  - It costs the extra pairs: 2,016 a wave, against 420.
- **Gravity is one routine for both engines,** `gravity_body`, and has
  no FMA in it. Its arithmetic, in order:
  - dx, then ((dx^2 + dy^2) + dz^2) + soft2;
  - a correctly rounded square root, then r^3 = (r r) r;
  - a correctly rounded divide;
  - each body's left fold from +0, in ascending partner index.

  A term added after the fold reaches every call site in both engines,
  because predict and correct read only `a0` and `at`. INFERRED from the
  call sites.
- **The program engine needs `--arith fma`; the loop engine takes both
  forms.** The forms differ only in three places: the predictor, the
  corrector's reciprocals and the end-of-step update. Gravity, `pc_error`,
  the step control and `predict_next_step` are shared.
- **A problem is a file of hex lines.** `read_problem`
  (src/ias15_cft.c:1012-1046) knows `name`, `G`, `N`, `E`, `system` and
  `body`, and anything else dies by name. Nothing can say "test
  particles" or "J2" today.
- **The replay cannot carry a rejected step.**
  - `--dt-file` needs `--epsilon 0` (src/ias15_cft.c:4412), so a replay
    never rejects.
  - `--dt-out` writes accepted steps only, and only for one system
    (:4421).
  - A rejection changes the next step's prediction, from `er` and `br`.
  - `--replay-iters` does carry rejected attempts.
  - The two options together are allowed (the run key hashes the dt
    list). No gate replays them together: tools/check_corrector.py runs
    the pair only to show it refused as a different run.
  - The run key hashed the masses, the initial state and the dt list, and
    nothing that would have said "test particles" or "J2". P0 keys both
    lines.
- **What a GPU wave must reproduce,** in order (`step_attempt`):
  1. gravity, then `a0 = a`, then `g` from `b`;
  2. passes of seven nodes (predict, gravity, `at = a`, correct), until
     `pc_error` exits;
  3. `choose_timestep` (PRS23, `dtnew_prs23`);
  4. the end-of-step update, with Kahan `add_cs`;
  5. `predict_next_step`.

  Its reductions are these:
  - `pc_error`'s two maxima over the system's coordinates, taken over
    normal values only and by bitwise selection, so order-free, then one
    divide;
  - PRS23's sums over each body's three coordinates, x then y then z,
    from +0;
  - PRS23's minimum over the bodies, which is order-free;
  - the clock.

  Everything else is per coordinate.
- **`CFT_REBOUND_ARTIFACT` goes straight to `cft_open()`**
  (src/ias15_cft.c:4347-4351). Once the pin carries G1B,
  `CFT_REBOUND_ARTIFACT=hip:0` opens a GPU with no change here.
- **`divfull` and `sqrtfull` are in the pinned library already**
  (third_party/cft-fp256/programs/divfull-*.cfta, at `922c2ad`).
  ROADMAP.md says "Neither is built" beside their prices. What is not
  built is their use in gravity, and the docs sweep says so.

**In cft-fp256** (read only):
- **G1A is on its main, at `fee1b9e2`.** libcft's softfloat, its bigint
  and the sequencer's executor compile as device code (`CFT_HD`).
  - `cft_sf_compute` does FMA, ADD, SUB and MUL; the compares, selects
    and integer operations; and the reciprocal and root seeds, in every
    format.
  - **It has no DIV and no SQRT.** divsqrt.c carries no `CFT_HD`.
- **G1B is not merged.** On its branches, `cft_open("hip:N")` runs
  elementwise calls, reductions and program images, one thread a lane,
  with resident `cft_alloc` buffers.
  - It refuses by name: index tables, AUGADD/AUGERR, any scratch depth,
    and `cft_get_image_id`. All of those are that backend's round 2.
  - Its `cft_div` and `cft_sqrt` are round trips orchestrated by the
    host.
  - It was checked bit for bit on the RTX 5060 Ti and the RX 7600.
- **No GPU has run `divfull` and `sqrtfull` yet.** On the host their
  images are byte-identical to the golden model's and held to its
  `softfloat.div` and `softfloat.sqrt`. No round plans a per-lane entry
  point for them.
- **The golden model, `python/cft_golden`, is in this tree's pinned
  clone.** It evaluates add, mul, fma, div and sqrt on bit patterns, in
  every format, and returns bits and flags. It runs here, without
  touching ../cft-fp256.

**The trap that would be a silent wrong answer: a flat ring.** Suppose
the particles and the moons all lie in the equatorial plane. Then dz = 0,
every z term of J2 and J4 is exactly zero, and a wrong z term passes
every gate at every format, bit for bit (METHOD.md §5, "inert by
construction"). So the ring's particles and moons are inclined (small,
but not zero), and P1's points leave the plane.

---

## P0 - two problem lines, and a seam for the zonal term *(the lead's)*

Every parcel reads or writes problem files, and three of them write
src/ias15_cft.c. So the lead lands the shared part first.

- **`active n`,** with n from 1 to N. The first n bodies of every system
  are massive and the rest are test particles, with `testparticle_type`
  0, through the engine's own `eng_active`. The record's energy follows
  it as REBOUND's does.
  - It is built in P0: it is one assignment, on a path the drop-in
    already gates.
  - The record says so in one line, with the pair count. That count is
    derived from the rows the engine built, not computed beside them.
- **`zonal k R J2 J4`.** Body k of each system is oblate about the z
  axis, with reference radius R.
  - P0 parses and stores it, and **refuses it by name until P1**.
  - The refusal's case belongs to P0, and P1 replaces it.
- **The seam.** `src/ias15_zonal.inc` is included by src/ias15_cft.c.
  - It holds three functions: `zonal_accel()`, called at the end of
    `gravity_body` when a zonal line was read; `zonal_potential()`, for
    `energy_all`; and `zonal_record()`, for the record.
  - P0's bodies are empty and cannot be reached, because the line is
    refused. P1 owns the file.
- **The run key.** The corrector schedule's key names both lines,
  appended only when a problem has them, so no older key moves. A
  schedule recorded with one N_active or one J2 is then refused under
  another.
- **The gate,** in tools/check_gravity.py. Its route matrix already runs
  problems whose rows differ, and `active` makes them differ again. (The
  drop-in's own N_active rows are tools/cases_gravity.c's and
  tools/cases_pairs.c's.)
  - **The problem:** outer with two more massless bodies, made from Pluto
    by a quarter and a half turn about z, so every value is exact.
  - **The route matrix** runs it with `active 5`, its rows read from the
    record's `# active` line.
  - **Inertness.** `active 5` and every body active must give the same
    record at every format, on the loop engine in both forms and on the
    program engine. The record says 25 pairs a force evaluation, and the
    test particles move. That is the inert claim, in the standalone
    program.
  - **The control:** one test particle given a mass of 2^-40, after which
    the two runs must differ. It stays in the gate. (A mass of 1e-30 would
    be inert at binary64, below its last bit.)
  - **testparticle_type 0, held.** Under `active 5`, a test particle's
    2^-40 mass must change nothing, the energy column included: it moves
    no active body, and REBOUND's energy (`reb_simulation_energy`)
    counts only the active bodies. V0 found both unheld in P0's first
    version, whose test particles were all massless. A copy at type 1
    fails this check, and so does one whose energy counts every body.
  - **An ensemble.** Two systems with `active 5` match every body active,
    and `--member 1` matches the system run alone.
  - **The refusals:** an active count of 0 or past N; a second active
    line; the zonal line; and a schedule recorded with `active 5`,
    replayed without it.
- **Behaviour is preserved.** No existing problem has either line, so
  `make check` must give merge 6's result: 841 ok and six digests
  unchanged (docs/VALIDATION.md entry 45). That is read from the log.
- **The round's ledger** goes in .ledger/, from ../ParcelRound's
  template, together with this file.

P0 lands and is pushed before anything is dispatched, so every worktree
branches from it.

---

## The parcels

### P1 - J2 and J4 in libcft's format *(the long pole; dispatch first)*

**Job.** Fill `src/ias15_zonal.inc`, then lift P0's refusal.
- **The acceleration:** each body's acceleration from body k's J2 and
  J4, computed in the run's format with the library's correctly rounded
  operations. It is added after `gravity_body`'s fold, in ascending body
  order.
- **The reaction on body k:** each body's term times its mass over k's.
  It is zero for a test particle.
- **The potential,** in `energy_all`, and **a record line.**

There is no FMA in the term, as there is none in gravity. So one order
serves both forms, and a later binary64 comparison with REBOUNDx's
harmonics can be op for op (decision 4).

**The order is one source.** The C is the definition, and its head says
which order it chose and why. The gates hold the C to the arithmetic and
to the physics; they never restate it.

The potential of body k (mass M, reference radius R) at r = (dx, dy, dz)
from it, with u = dz/r:

    Phi = -(G M / r) [1 - J2 (R/r)^2 P2(u) - J4 (R/r)^4 P4(u)]
    P2 = (3u^2 - 1)/2,   P4 = (35u^4 - 30u^2 + 3)/8

**Gates** (tools/check_zonal.py, one leg):
1. **Arithmetic.** `ias15_cft --zonal-selftest` evaluates the term at a
   fixed set of points: on and off the plane, near and far, both
   hemispheres. It traces every operation it performs, with its operands
   and result in hex.
   - The golden model checks each traced operation, bits and flags, at
     all three formats.
   - The trace comes from the C itself, so nothing is transcribed.
2. **Formula.** The same points against -grad Phi, evaluated in mpmath at
   200 digits from the potential above, not from the force's code.
   - At binary256 they must agree to the format's own rounding, a
     relative 1e-68 or better.
   - A wrong coefficient misses by tens of decades, since the coefficients
     are small exact rationals.
3. **Physics.** Saturn and one inclined, slightly eccentric particle:
   - **(a)** the energy including Phi is conserved to the integrator's own
     level, at every format;
   - **(b)** the apsides advance at n - kappa and the node regresses at
     n - nu, first-order theory in J2 and J4, within the size of the
     second-order terms (about a percent). The rates are printed beside
     the comparison and asserted.
4. **Format difference.** One run at binary128 and at binary256 on the
   same steps (`--dt-out`, then `--dt-file`). The positions must differ
   at binary128's round-off, not binary64's.

**Negative controls.** Build each from a copy, run it, show it fails,
delete it, and report it with its numbers:
- the term rounded to binary64 before it is added: gate 4 fails, at
  about 1e-16;
- J4's z part dropped: gate 2 fails off the plane and passes on it. That
  is why the points leave the plane;
- the force's z sign flipped and the potential left alone: gate 3(a)
  fails;
- one traced operation's result corrupted by one ulp: gate 1 names it.

**Files you own:**
- `src/ias15_zonal.inc`;
- tools/check_zonal.py and data/problems/zonal_*.txt;
- the constants it needs, in data/ias15_constants.json and the header
  tools/gen_constants.py makes;
- tools/check_scalar.py's rows, if its recount of call sites needs them.
  It recounts src/ias15_cft.c alone, so calls in the include are
  invisible to it: P1 extends the recount, or says why the term has no
  site it would count;
- in src/ias15_cft.c, only P0's refusal line and the `--zonal-selftest`
  option;
- one leg in each of the Makefile's three tiers.

**The trap.** A test with dz = 0 is inert for every z term, so gates 1
and 2 sample off the plane and gate 3 uses an inclined orbit. And J2 and
J4 mean something only with their reference radius. The study's J4 has
none stated (it believes 60,330 km), so P1 states the pair it used, with
its source.

### P2 - The ring: a generator, waves and ensembles, run at every format

**Job.** `tools/make_ring.py` writes ring problems: Saturn with its zonal
line, six moons and N test particles. It writes them as one system, or as
an ensemble of 64-body waves (7 massive bodies and 57 particles each),
with the particles sorted by radius before they are cut into waves. A leg
then runs them.

**What it writes:**
- **Saturn, the moons and the annulus** come from the study's READ table
  (GPU-E §4: NASA's fact sheets, Cassini's J2 and J4, JPL's moons), each
  value with its source in the file's head. Units are km and s, with
  G = 1 and masses as GM.
- **The moons' phases and inclinations are synthetic:** fixed, documented
  values, unless Logan asks for an epoch's (decision 5). So this round
  needs no web read. The inclinations are small and not zero. Epimetheus
  and Janus are co-orbital, 100 km apart, so their phases stay far
  apart. Particles that would start within a few Hill radii of Pan,
  Prometheus or Pandora (all three are inside the annulus) are excluded
  and counted.
- **Particles** go on near-circular orbits in the oblate field, at the
  circular speed with J2 and J4. Their inclinations are small and not
  zero, and their radii and phases come from a seeded generator. Every
  value is written exactly, as hex.
  - mpmath does anything binary64 does not round correctly (the trig), so
    the file should be the same on every host. That is BELIEVED; the
    check across hosts waits for a WSL slot.
- **The massive bodies come first,** so `active 7` is right for every
  system.

**Gates** (tools/check_ring.py, one leg):
- the generator twice gives the same bytes, and a different seed gives
  different ones;
- one wave runs at binary64, binary128 and binary256: on both engines in
  the FMA form, and on the loop engine in REBOUND's. The engines agree
  bit for bit, as check_program_engine holds them;
- in an ensemble of k waves, each member gets the bits it gets alone
  (`--member`, as check_ensemble does);
- `active 7` against all 64 bodies active gives the same bits, at every
  format. That is the inert claim on the real problem;
- **physics counters, each asserted:**
  - 420 pairs a force evaluation;
  - the moons do something: a run with them differs from one without;
  - no particle leaves the annulus or reaches Saturn over the run;
  - the steps an orbit for a wave far from the moons are within 10% of
    the study's MEASURED 36.2 at epsilon 1e-9. That was measured on a
    circular orbit about a point mass, with no J2 and no moons, so a
    miss is a finding to report, not a reason to loosen the bound.

The rings run without the zonal line until P1 lands: the generator has a
switch. The zonal line's refusal is P0's case, which P1 replaces, and the
lead's P1 x P2 seam test runs the ring with the line.

**Negative control.** In a copy of an ensemble, move one particle's
starting x by one ulp. The comparison against `--member` must fail for
that member, and only for it.

**Files you own:** tools/make_ring.py, tools/check_ring.py, and one leg.

### P3 - A particle replayed out of its wave, rejected attempts included

**Job.** Make one particle's trajectory reproducible, bit for bit, from
its wave's record, in a system of Saturn, the moons and that particle
alone, whatever the wave did. That is how a GPU run's particles get
checked at scale: a wave makes 420 pairs a force evaluation, and the
particle's own system makes 28.

**What it needs:**
- **The attempts.** It needs a record of every attempt, not only of every
  accepted step: the dt tried, whether it was accepted, and the passes.
  Add an attempts form to `--dt-out`, and a replay that follows it,
  rejections included. On a recorded rejection that means the same
  restore, and the same prediction from `er` and `br`.
- **A run key that admits the sub-system.** `--replay-iters` refuses a
  run with a different key, and the particle's system is a different
  problem by construction. P3 picks one way through and says why:
  - the key gains a part that names the massive bodies and the schedule,
    and leaves the test particles out; or
  - a tool re-keys the schedule for the sub-system, after asserting that
    the sub-system's massive bodies are the wave's, line for line.
- **The extraction.** A tool writes the sub-system's problem (the wave's
  massive bodies and zonal line, and particle j) and its schedule.

**Gates** (tools/check_replay.py, one leg):
- On a small ring wave at every format, each of k sampled particles,
  replayed alone, gives the wave's bits for that particle and for every
  massive body, at every sample.
- A wave forced to reject, by a first dt far too large, passes the same
  test.
- **Physics counter:** the replayed schedule's rejected attempts,
  asserted not zero in the rejection case.

**Negative controls:**
- one attempt's dt moved by one ulp: the particle's bits differ;
- the schedule with its rejected attempt deleted: the replay is refused
  by name, or differs.

Round 3 found a pass count moved by one inert on most attempts (docs/
VALIDATION.md entry 41). So these controls move a dt, not a count.

**Files you own:**
- in src/ias15_cft.c, only the dt record's writer and reader, the step
  loop's replay branch and the run key;
- tools/extract_particle.py and tools/check_replay.py;
- one leg.

P3 starts before P2 lands. Its waves are hand-made (`active` is P0's), and
the lead's seam test uses P2's.

**The trap.** The sub-system's massive bodies match the wave's only
because test particles add exactly nothing to them. A "test" particle
with a tiny mass would make the massive bodies differ, while its own bits
could still match. So the gate compares the massive bodies too.

### P4 - A run's records against this tree

**Job.** `tools/compare_run.py` takes a ring run's records: today this
tree's own ensemble run, and later a GPU's.
1. **It samples waves:** seeded and stratified by radius, plus every wave
   whose record shows a rejection or a close approach to a moon.
2. **It samples particles** from those waves.
3. **It reruns them here:** each wave as `--member`, and each particle
   through P3's replay once that lands.
4. **It reports each one** as identical or not, by name, with the physics
   counters beside the bits:
   - steps and passes an orbit, and rejections;
   - each particle's closest approach to each moon;
   - the format differences: binary64 against binary128 on the same
     steps, and binary128 against binary256.

**The seam this parcel must not leave to nobody.** Until a GPU engine
exists, the run being compared is this tree's own, so every case compares
the CPU with the CPU. That is METHOD.md §1's canonical failure: a gate
driving a stand-in.
- The driver reads only the record format, never ias15_cft's internals.
  A GPU engine that writes the same records is then held without
  changing the driver.
- The lead owns the test that points the driver at the first GPU run
  ("What the integrator keeps").

**Gates** (tools/check_compare_run.py, one leg): the driver on a small
ensemble run reports every sampled member identical. Negative controls:
- one byte changed in one member's record: FAIL, naming that member and
  line;
- a member missing: FAIL, naming it;
- a record from another format: refused by name.

**Files you own:** tools/compare_run.py, tools/check_compare_run.py, and
one leg.

### P5 - The wave kernel at the hardware's binary64 *(proposed; needs Logan's yes and a GPU slot from cft-fp256)*

P5 is not dispatched with the rest. The kernel the study prices (its §6,
items 2 and 7) needs a divide inside a kernel, and none exists.
- **Why binary64 now.** At binary64 the hardware's own divide and square
  root are IEEE's correctly rounded operations, on both vendors. That is
  BELIEVED; P5's first gate measures it on a corpus before anything relies
  on it.
- **What it would build.** The wave kernel's structure could be built now
  and held to this round's CPU references at binary64: the state laid
  out by lane, one thread a body, the wave's maxima and minimum, the
  masks, J2 and J4, and the test particles. It is written over an
  arithmetic type that libcft's device arithmetic replaces later.
- **What it would measure:** what the study could not, a real IAS15
  kernel's registers and rate (on the hardware's doubles, not
  softfloat's).
- **Where it runs.** Its machines are cft-fp256's: the 5060 Ti in WSL
  `cft-gpu`, and the RX 7600 on nvidia-box. So it runs only in a slot
  cft-fp256 gives, one agent at a time.

---

## Rules for every parcel

- **Work in a worktree,** in the background. **Never push, merge or
  rebase**; the integrator does all three. Commit in your worktree only.
- **Check your base** with `git log --oneline -1`, against the commit in
  your dispatch message.
- **Build in your own tree.** On Windows the PATH line is not optional:

  ```bash
  export PATH="/c/msys64/mingw64/bin:$PATH"
  make OS=Windows_NT TMP="$TEMP" TEMP="$TEMP" \
       PYTHON="/c/Users/logan/AppData/Local/Programs/Miniconda3/python.exe" <target>
  ```

  third_party/ is not in a worktree: clone it from the main checkout's
  third_party/ at the pins (a local clone, not a download).
- **Run `make check-light` while you work** (31 minutes cold). Before you
  report, run a full `make check` (2 h 36 min on 2026-11-06, longer
  today: this host is shared and was at 100% load at 17:22). Start it from
  PowerShell (CLAUDE.md, "Launching a long suite from an agent"). **Read
  the log, not the exit code.**
- **A new leg names every artifact it executes as an `--input`**
  (tools/gate_cache.py).
- **Follow the `--light` scaling rules,** and say which kind each new case
  is (CLAUDE.md).
- **Bit-identity is not validity.** Print the physics counter beside the
  comparison, and assert it.
- **Nothing compares the trailer for you.** Derive; never transcribe.
- **A gate that cannot fail is not a gate.** Run your named control,
  confirm it fails, delete its artifact, and report both results.
- **READY means "a gate, or a stated limit".**
- **The ledger.** Read every file in it at three points: before you
  start, before you design anything this plan calls shared, and before
  you report. Arm the `urgent/` watch first. Append to your own file
  only.
- **Label every claim** as measured or believed.
- **Never run bare `git stash`.**
- **Read ../cft-fp256; never write it.** Run the golden model from
  third_party/cft-fp256/python with `python -B`. Touch nothing in
  ../ParcelRound.
- **No GPU, no WSL, no card.** The GPUs and WSL are cft-fp256's to
  schedule, and WSL work there is one agent at a time. amd-arc-box's U50
  has cft-fp256's nightly card run, from 18:00 to about 08:00.
- **No speed claims.**
- **Report:**
  - what changed, file by file;
  - the gate output, as actual lines;
  - the controls, with their numbers;
  - what you did not do;
  - anything you found that this plan got wrong.

## What the integrator keeps

- **P0,** before anything is dispatched.
- **The seam tests.** Each merge gets one that exercises two parcels
  together:
  - **P1 x P2:** a ring with J2 and J4, at every format, on both engines,
    with the precession counter on a ring particle.
  - **P2 x P3:** a particle of a generated wave, replayed out of it.
  - **P1 x P3:** the same, with the zonal line. The sub-system keeps
    Saturn's J2, and a particle's reaction on Saturn is exactly zero.
  - **P3 x P4:** the driver's particle sampling, through the replay.
  - **The GPU seam:** when a GPU engine first writes a record, the driver
    is pointed at it. Its first case is a difference planted in the GPU's
    record, shown caught. This goes into the GPU half's first brief, and
    is not left to that parcel.
- **The integrator-only files:**
  - README.md, ROADMAP.md and CLAUDE.md;
  - docs/VALIDATION.md and docs/HARDWARE.md;
  - PARCELS.md, ROUND3.md and this file;
  - third_party/MANIFEST;
  - a new docs/RING.md, at the round's end.

  A parcel that needs one changed says so.
- **Every merge, serially, with a full `make check` after it.** Never
  merge while a suite runs. Every parcel adds a leg to the Makefile's
  three tiers, and the lead resolves those conflicts at each merge. P2
  and P4 share no file beyond a leg each, so the lead may merge both and
  run one suite, and will say so.
- **A verifier per parcel, V1-V4.** Each re-runs the gate from a clean
  build, re-runs the control, diffs the change against the ownership
  lists, and fixes nothing.
- **The VALIDATION entries.** Entry 46 is P0's; the parcels' start at 47.
- **The summary to the owner,** with the cost per agent and the
  verifiers' share.

## Order of the day

| when | who | what |
|---|---|---|
| tonight | lead | P0 built, its `make check` run and read; the ledger and its `urgent/` watch armed. |
| before dispatch | Logan | The merge and push of P0. Round 3's delegation ended (decision 1). |
| dispatch | lead | **P1 first** (the long pole), then P2, P3 and P4. |
| as each returns | V1-V4 | Verify. A send-back resumes the parcel with its context. |
| as each lands | lead | Merge, run the seam test and the full `make check`, and read the log. |
| last | lead | The entries, docs/RING.md, ROADMAP, README and CLAUDE.md. |

If the round runs short, **P3's rejected attempts drop first, to a stated
limit.** P4 then reruns the whole wave wherever a rejection happened.

---

## After this round: the GPU half

1. **The pin moves to G1B** once it merges, at Logan's word. `make check`
   must show the software bits unchanged.
2. **This tree on `hip:0` as it is.** This is the study's route with no
   port. It is slow by design: the per-call structure HARDWARE.md measured
   on the card.
   - It finds what G1B refuses on this engine's default path. Index
     tables are opt-in here, so they should not reach it. The clock's
     augmented sum calls `vaugadd` at every step and G1B refuses AUGADD,
     so that very likely does. Whether the scratch or the image id do is
     not known.
   - It holds a GPU's bits to software's on the existing gates, at small
     sizes.
   - It would be the first IAS15 run on a GPU, and it needs a slot from
     cft-fp256.
3. **A divide and square root inside a kernel.** This is cft-fp256's work
   (the study's question 2). The images and the executor exist; what is
   missing is a per-lane entry point and its gates on a GPU. This tree
   asks for:
   - binary64, binary128 and binary256;
   - round to nearest;
   - the five flags, since the engine aborts on four of them;
   - host `cft_div`'s and `cft_sqrt`'s bits.
4. **The wave kernel** (the study's a2; its §6, items 2 and 7): one wave a
   system, held to this round's CPU references through P4's driver. Its
   prototype at binary128 measures registers, scratch and rate.
5. **The ladder:** the RX 7600 at 3 x 10^4 particles, then the RTX 5060
   Ti at 10^5. The study's ESTIMATEs, generic tier, by format: 10 min to
   4.3 h on the 7600, and 11 min to 7.7 h on the 5060 Ti. Each rung's
   records are compared with the others' and with the CPU's.
6. **The MI300X,** when AMD's time is approved: the runbook, its
   checkpoints and records, and the sampled comparison.

## Decisions that are Logan's

1. **Merges and pushes for this round,** since round 3's delegation
   ended. P0 must land before anything is dispatched. Meanwhile P0 is
   built and its suite run, and nothing merges.
2. **The study's questions 1, 3 and 5:** the workload, the form and the
   formats. This round's CPU half serves every answer; the GPU half needs
   them.
3. **Its question 2, the device divide before cft-fp256's exp family.**
   It decides when the GPU half's step 3 starts.
4. **REBOUNDx as a pinned upstream** (a download), for an op-for-op check
   of J2 and J4 at binary64 against what REBOUND's users run. The default
   is not to pin it: P1's formula gate holds the formula without it.
5. **The moons' phases.** Synthetic by default; an epoch's needs JPL's
   tables read from the web.
6. **P5:** yes or no, and then a GPU slot from cft-fp256.
7. **The pin's move to G1B,** once it merges.

## Not in this round, and why

- **Anything on a GPU or the card:** that is the GPU half, above.
- **J2 and J4 for the drop-in.** REBOUND's users have REBOUNDx's at
  binary64. A wide one would need an API of this port's own, and nothing
  has asked for it.
- **J6 and beyond.** Saturn's J6 term is about 1e-6 of the central force
  at the Encke gap. The line can grow when a run needs it.
- **A self-gravitating ring, or collisions.** The study prices them at
  10^3 to 5 x 10^4 times the test-particle ring. This port refuses the
  tree code and shear boundaries by name (README.md:47-58).
- **Round 3's list in ROADMAP.md.** It waits. Its item 3, the divide and
  square root in gravity, is the same mechanism as the GPU half's step 3,
  `divfull` and `sqrtfull` as programs, so whichever lands first serves
  the other.
- **The software backend's cost** (ROADMAP.md item 2) has a proposal of
  its own.

---

## What the round found against this plan

Each item is measured unless it says otherwise. The urgent files and
ledger entries named are in the round's ledger, archived as
`archive/round4-ledger.zip`; docs/VALIDATION.md carries the numbers.

1. **P1's gate-4 control could not fail as written.** The term rounded
   to binary64 before it is added, in both runs, leaves binary128 against
   binary256 on the same 120 steps at 3.371e-33; the correct build gives
   3.402e-33. Rounded at binary128 only: 4.029e-18 (P1,
   `urgent/P1-gate4-blind.md`). The lead's ruling of 2026-11-12 01:33
   (`urgent/lead-gate4-ruling.md`): gate 4 stays with its scope stated,
   its control is the binary128-only rounding, the control above is
   shown to fail gates 1 and 2, and P4's format difference prints its
   scope. Both controls were run from copies, by P1 and again by V1, and
   deleted: neither is a permanent control in the gate. Entry 47.
2. **"The replay cannot carry a rejected step" holds only for a
   rejection before the first accepted step.** After an accepted step a
   rejection leaves no state bit: deleted from the accepted-steps replay,
   0 of 7 record lines differ at binary64 and binary128, where one
   first-attempt rejection moves 5 of 5, 7 of 7 and 7 of 7 on three
   problems (P3's ledger, 2026-11-11 21:58; P4's member 8; V3). P3 built
   the attempts form all the same, for those waves and for the counters,
   so the fallback above never ran, and it would have needed only waves
   that reject first. Entry 49.
3. **10^6 particles is 17,544 waves of 57, 1,000,008 particles.**
   `--particles 1000000` is refused by name, and an ensemble's systems all
   have N bodies, so a 49-particle last wave cannot be written (P2's
   ledger, 2026-11-11 22:35; the lead's, 2026-11-12 05:59). And "each wave
   on its own adaptive step" hides a cap: a wave outside Pan's orbit steps
   at Pan's rate, about 3% more steps over the ring, computed by P2 from a
   uniform radius. Entry 48.
4. **The parcels' clone loop was refused by the worktree guard.** The
   briefs' `while read ... git -C` loop over third_party/ was refused by
   the worktree-isolation check; plain, separate `git clone` and
   `git checkout` commands worked (P2's and P4's reports).
5. **Decision 1 was delegated.** Logan delegated round 4's merges and
   pushes to the lead on 2026-11-11, on round 3's bar (the lead's ledger,
   21:15). Entry 51 has the merges.
6. **The lead's own seams and wiring got a verifier, VL.** No parcel's
   verifier saw the fix after V1, the seams, the a5 hardening, the P3 x P4
   wiring or the merges, so VL was dispatched on them (the lead's ledger,
   2026-11-12 05:52). It sent them back twice, each time for a control
   that held less than its comparison, and was READY on `d9c8c8c`.
   Entry 51.
