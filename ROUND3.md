<!-- SPDX-License-Identifier: GPL-3.0-or-later
     Copyright 2026 the cft-rebound contributors. -->

# Round 3: up to date with the library that was built for this port

**Status, 2026-10-31: P0 is built and verified twice, and nothing is
dispatched.**
- The re-pin, P0(a), is merged and pushed (`ebaeaa1`).
- P0's first build is on branch `round3-p0`. V0 sent it back, and the
  fixes are on `round3-p0-v1`.
- V0's scoped re-check confirmed those fixes, and found four gaps and
  one more route divergence. They are fixed on the same branch.
- P0 lands when its full suite is green and the owner says so.
- The rest of this plan is the owner's to approve, change or cut.

Round 2 ([PARCELS.md](PARCELS.md), which stays as that round's record
because the tree cites it) took the port from correct to usable. Since
then this repository has stood still and its library has not. The last
commit here is from 2026-09-14: entry 37 and the work it scoped. In the
six weeks after that, cft-fp256 went from ABI 0.12 to 0.18:

- it built all eight asks this port filed (its docs/ROADMAP.md, "What
  a real workload asked for"), including the gather and scatter that
  docs/VALIDATION.md entry 37 ranked first;
- it built certificates, which nobody here asked for (its
  docs/CERTIFICATES.md, from 2026-09-28);
- it cut program runs across tiles (2026-09-25) and let a caller open a
  subset of an image's tiles (`CFT_XRT_TILES`);
- it took the divide ask, built it, and measured the result slower
  (ask 8).

Entry 38 moved the pin to `922c2ad` (ABI 0.18) on 2026-10-30. It showed
that the calls the port already makes return the same bits, and nothing
in the port calls anything new. This round is the adoption, in entry
37's order. It also covers the two pieces of state the archive still
does not carry, and the verification entry 38 left open.

**The target.** On the card:
- a step whose gravity is a handful of program runs, not hundreds of
  calls;
- a corrector exit test that does not loop on the host;
- a checkpoint that is complete at every format;
- a record that says which library build and device image made it;
- every bit where it is now.

`METHOD.md` and `CASE-STUDY-4.md` below are `../ParcelRound`'s. That
copy of the method is canonical, and `docs/PARCEL-ROUNDS.md` is an
older copy that loses wherever the two differ.

---

## What this round is built around

### The shared thing: new entry points, mostly behind capability bits

Four of this round's five parcels rely on a library entry point the
port has never called:

| parcel | entry point | ABI | capability bit |
|---|---|---|---|
| P1 | an indexed program run (`cft_run_args.idx_*`) | 0.14 | `CFT_SEQ_FEAT_INDEXED`, CAPS2[9] |
| P2 | a segmented reduction (`cft_reduce_seg`) | 0.13 | `CFT_FEAT_REDUCE_SEG`, CAPS2[8] |
| P5 | a scalar-broadcast elementwise call (`cft_run_ex` with `scalar_mask`) | 0.12 | `CFT_SEQ_FEAT_SCALAR`, CAPS2[7] |
| P3 | the identity calls `cft_build_id()` and `cft_get_image_id()` | 0.15, 0.18 | none; refused by name on the software backend |

Each of the first three needs the same four things:
- a probe of the bit;
- a fallback to the route the port takes today, which gives the same
  bits;
- a wrapper that times and counts the call in the trailer's buckets;
- a trailer line saying which route ran.

The last is the one that matters most. The trailer is a `#` line, and
records are compared on non-comment lines only (entry 38), so a card run
that silently took the fallback would pass every gate. Built four
times, these drift. So **P0 builds them once, and lands before any
parcel starts.**

The trailer held a second shared fact, and it was already stale.
[hw/report.py:298](hw/report.py:298) names the four keys its figure
reads (`t_elem`, `t_program`, `t_divsqrt`, `t_gravity`). The gravity
split that ROADMAP item 4 says decides P1's shape (`t_g_copy`,
`t_g_accum`, `t_g_lib`, `t_g_divsqrt`) was in no list any script read.
No script in the tree wrote `docs/bench/decompose.csv` at all. That is
METHOD.md §1's failure in miniature: a transcribed list that a new key
walked past. P0(d) below is the fix.

### Measure what decides, and build what the measurement will choose between

Entry 37 found the residency plan "costed correctly and ordered
wrongly". Its instruction to measure first was the part that was right.
The measurement that orders P1 is ROADMAP item 4, the gravity split,
and it needs the card. The card is cft-fp256's until Mon 2026-11-02
07:45 (its weekend window), with its nightly runs after that.

So this round neither waits for the card nor guesses. P1 builds both
pair-gather routes behind P0's route switch, and the card day picks
between them without a rebuild. **No parcel claims a speed-up.** A
software timing says nothing about a tile that charges 108 µs a call,
so every speed figure this round produces is the lead's, from the card
day. That is the lead's choice and the owner's to overrule. Holding P1
until the split is measured would cost a day and save building one
route that may lose.

### Read the provider's code, not its record

METHOD.md §1 says a round that answers a request reads the requester's
code. This round is the requester, so the rule turns around. Every
claim below that comes from cft-fp256's docs is **believed, not
measured here**, and any brief that relies on one says "verify it, do
not trust it". Three of them are load-bearing:
- that `normalabs-<fmt>` with `cft_reduce_seg` is `pc_error` exactly,
  in both modes (P2);
- that the indexed run folds in table order from +0 (P1);
- that the new library opens this repository's f128 image (the card
  day). Contract 0x800 is in `KNOWN_VERSIONS` at the pinned
  `third_party/cft-fp256/host/src/backend_xrt.cpp:267`, but nothing
  here has opened the image through 0.18.

---

## Before dispatch

| # | who | what | state |
|---|---|---|---|
| 1 | owner | Merge and push `repin-cft-fp256` (`ebaeaa1`). This is P0(a). | **done** 2026-10-31: the lead fast-forwarded main and pushed it, on the owner's word |
| 2 | lead | `make check` on Linux, in the cft2204 WSL distro. It closes entry 38's Linux gap and gives entry 29's open item its first Linux result. amd-arc-box stays untouched while cft-fp256 has it. | **done**: see below |
| 3 | owner | `make check-python` needs a `rebound` wheel in the Miniconda env. | **done**: see below |
| 4 | lead | Build P0(b)-(e) below. The suite must be green before and after, and entry 38's 27-integration A/B must show the records unchanged. | **done**, twice: see below |
| 5 | verifier V0 | Verify P0 before any parcel that reads it is dispatched. This is CASE-STUDY-4's proposal, which the owner made the precondition of that study's first round. | **done**: sent P0 back, then re-checked the fixes; see below |
| 6 | lead | Archive round 2's `.ledger/` and stand up round 3's, `urgent/` included. | **done**: see below |
| 7 | lead | Copy `~/patch_fixed_schedule.py` from amd-arc-box into the ledger for P2. Read only; no card. | before P2 is dispatched |
| 8 | owner | Prune the ten stale agent worktrees under `.claude/worktrees/`. | **nothing to prune**: see below |

Rows 2, 4, 5, 3, 6 and 8 in more detail:
- **Row 2, Linux.** This was `make check` at `ebaeaa1` in cft2204 (gcc
  11.4, glibc), from local clones only.
  - It passed: rc 0, 14 legs and 222 `ok` lines, the same count as
    Windows at that pin. `make example` passed too.
  - All six wide-state digests are identical to Windows'.
  - It found two warnings that mingw hides, both in this repository's
    tests, and both are fixed on `round3-p0-v1`. `tests/gate_real.c`
    printed a `uint64_t` with `%llu`. `tests/gate_subprocess.c` called
    `putenv()` undeclared under `-std=c99`, which gcc 14 and later
    reject.
- **Row 4, P0.** The first build, `971f93d`, passed the full `make
  check` on Windows: 15 legs, 250 `ok` lines and 0 FAIL, with the
  digests unchanged.
  - An A/B of 27 integrations against the pre-P0 binaries showed 27 of
    27 identical, with the control seeing its difference and the 27
    program images byte-identical.
  - `round3-p0-v1` repeats both.
- **Row 5, V0.** V0 rebuilt both trees from clean and confirmed that no
  record line moved, across 51 integrations of its own choosing. It
  sent P0 back on what P0 hands the parcels:
  - five ways to defeat the gate that it did not catch, one passing
    with nothing compared;
  - two wrappers that answer differently by route on aliased calls;
  - three false claims;
  - smaller items.

  The fixes are listed in P0(c) and (d) below. V0's entry is
  `.ledger/verifier-P0.md`.

  Its scoped re-check of `90e6f3f` re-ran the records (51 of 51), the
  gate, `check-light` (15 legs, 0 FAIL) and the Python legs. It
  confirmed each fix against its finding, and found reading a scalar
  once correct across five aliasing shapes. Four of its ten new
  sabotages passed the gate:
  - F4: the data's low bytes never checked;
  - F6: the hide-variable typo never tested;
  - F11: no two-operand reduction;
  - F4 combined with N8.

  It also found `vrun_ix` still splitting by route on a table on a
  NULL operand, and two limits stated below. All are dealt with at this
  branch's tip. Its entry is `.ledger/verifier-P0-recheck.md`.
- **Row 3, Python.** The owner allowed the download, and the lead
  installed rebound 5.1.1, the `cp312-win_amd64` wheel at githash
  `33549d1d`.
  - `check-python` gives "IDENTICAL: 13 values, bit for bit", and its
    binary128 control differs in 9 of 13 values.
  - That was true only after a fix. From `24900ee` (2026-09-14) the
    Windows DLL imported `libwinpthread-1.dll`, which no Python process
    can find, through the trailer's `clock_gettime()`. `47f0bc3` links
    it static and refuses any `lib*.dll` import; docs/PYTHON.md, "The
    Windows rule".
- **Row 6, the ledger.** The owner said to archive it in this
  repository.
  - Round 2's is `archive/round2-ledger.zip` (`21579ee`): six files,
    each with its original timestamp, verified byte for byte.
  - Its working copy went to the Recycle Bin, not deleted.
  - Round 3's is `C:\Users\logan\source\repos\cft-rebound\.ledger\`:
    the template's README, `urgent/`, and three seed entries from the
    lead.
- **Row 8, the worktrees.** All ten were already gone, with their
  registry entries, by 2026-10-31. Every local branch but `round3-p0`
  is merged into main.

Parcels branch from main **at or after** P0's merge. The exact tip goes
in each dispatch message, not here (METHOD.md §2).

---

## P0 - one door to the new library *(the lead's)*

**(a) The pin. Built:** `ebaeaa1`, docs/VALIDATION.md entry 38.
- `make check` passes at both pins, 222 `ok` lines each, with identical
  wide-state digests.
- 27 of 27 integrations are bit-identical across the pin.
- Checkpoints cross the pin in both directions, 12 of 12.

It has not been run on the card, on Linux or under Python; see rows 2
and 3 above, and the card day.

This section describes P0 **as built**. Where building it corrected the
plan, the correction is marked. Every line number below is at
`round3-p0-v1`'s tip, which parcels branch from once P0 lands.

**(b) The routes, probed and switched in one place: `engine_routes()`
([src/ias15_cft.c:214](src/ias15_cft.c:214)).**

*Correction.* The plan put the probe in `probe_device()`. That is
`--probe`'s reporter, and it runs only in the program. The function both
paths already call at open is `engine_check_format()` (:114). So
`engine_routes()` runs right after it, in `main()` (:2689) and in
`ias15_engine_open()` (:2978) alike.

It reads `cft_caps.seq_features` once for three bits (`cft.h` here means
the pinned `third_party/cft-fp256/host/include/cft.h`):
- `CFT_SEQ_FEAT_SCALAR`, at `cft.h:891`;
- `CFT_FEAT_REDUCE_SEG`, at :931;
- `CFT_SEQ_FEAT_INDEXED`, at :974.

The policy is stated once, in `route_def[]`, for the routes `scalar`,
`indexed` and `reduce_seg`:
- `CFT_REBOUND_ROUTE_<NAME>` takes `auto` (the default: the device
  decides), `off` (today's route, whatever the device publishes) or `on`.
- `on` against a device without the bit is **refused by name**. The
  program refuses with `die()`. The library returns `-4` with the
  sentence in `ias15_engine_error()`, as it does for an absent format,
  so the shim needed no change.
- `CFT_REBOUND_HIDE_FEATURES=<name>[,<name>]` is the test-only hook
  that makes a published bit look absent.
- A bad value or a bad name is refused by name too. So is a
  **misspelled variable**: any `CFT_REBOUND_ROUTE_*` or
  `CFT_REBOUND_HIDE*` the door does not know. The prefix is matched
  without case, as Windows' `getenv()` does. *V0 found that a typo used
  to leave its route on `auto` with no sign.*
- **The library half can see its routes** through `ias15_engine_route()`
  (`src/ias15_engine.h:318`). A library-path gate for a routed call
  reads it, as the program's gates read the `# route` lines.

Each family's calls and fallbacks are counted, and the trailer prints
one line per family: `# route <name>=on|off by=device|absent|hidden|forced
calls=N fallbacks=N`. `--probe` names the three bits as the device
publishes them.

**(c) The wrappers that the contract alone defines.**

*Correction.* The plan had P0 build a program-run wrapper with index
tables. Its shape is P1's design: which streams and blocks are indexed,
the lane count, the deposits. Building it ahead would have fixed an API
before its only consumer existed. **It moves to P1.** P0 built the three
wrappers whose behaviour the contract alone fixes:
- `vrun_k()` (:506): a scalar-broadcast elementwise call over
  `cft_run_ex`, or a broadcast vector on today's route;
- `vrun_ix()` (:572): an indexed elementwise call through `idx_*`, or
  the host's gather on today's route, with `CFT_IDX_NONE` reading +0;
- `vreduce_seg()` (:641): `cft_reduce_seg`, or one `cft_reduce` per
  segment, which is its definition.

Each is timed into its neighbours' buckets: `t_elem`, the new `t_reduce`,
and `t_g_lib` inside gravity. Each counts **one** operation in `calls=`,
whichever route carries it. So `calls=` compares equal between a card
without a bit and a backend with it, and the route line says which route
ran.

They are marked `DOOR_UNUSED` until a parcel's call lands: P5 uses
`vrun_k`, P1 may use `vrun_ix`, and P2 uses `vreduce_seg`.

**Two rules that V0's send-back added, which a parcel inherits:**
- **A scalar operand is read once, before any element is written**, on
  both routes. libcft's software backend reads element 0 again after
  writing it (cft-fp256 `host/src/device.c`, 2068-2086). With the
  result aliasing the scalar, V0 measured 3,838 to 4,095 of 4,096 lanes
  differing between routes. Whether re-reading is libcft's intended
  contract is a question for cft-fp256's owner.
- **Once any operand has an index table, the result may overlap no
  operand**, dense or indexed, on both routes. That is `cft.h`'s rule.
  Before the fix, today's route computed what the new route refused.
  - The same goes for a table on a NULL operand and a table beside a
    zero-element source, which V0's re-check found split by route.
  - One of `cft.h`'s rules is a **stated limit**: a table on an operand
    the opcode does not read. libcft exports no arity, and one typed
    here would drift from it. The new route refuses such a call by
    name. Today's route ignores the operand, so its result is the one
    the opcode computes anyway.

**They are held by `--door-selftest` (:2478) and its gate,
`tools/check_door.py`, in all three suite tiers.** The self-test runs
131 cases at each format: 98 scalar (two of them with the result
aliasing the scalar), 15 indexed, and 18 segmented, six of those a
two-operand dot product.
- **Both routes are compared against a reference computed outside the
  wrappers**, straight from libcft. A defect both routes share is then
  caught too. *V0's N1, the mask ignored on both routes, passed a
  route-against-route comparison.*
- **Every wrapper call's accounting is checked:** one operation in
  `calls=`, and exactly the counter meant. *V0's N2 and N3, a swapped
  counter and a fallback counted per segment, passed before.*
- **The data uses full-width mantissas at binary128 and binary256.**
  Widened binary64 left every operand's low bytes zero there. That
  layout, little-endian with the sign in the last byte's top bit, was
  measured on 2026-10-31. Every scalar case also runs with a scalar of
  +0 and of -0. *V0's N7 and N8, dropped low words, passed before.*
  - **The data checks itself.** Each nonzero value goes through
    binary64 and back in libcft, and the gate requires all of them to
    come back different at binary128 and binary256, and none at
    binary64. Measured: 9,422,984 of 9,422,984 and 9,425,290 of
    9,425,290. *V0's F4 had turned the randomising off unseen.*
- *Correction:* **there is no "largest shape the integrator presents".**
  `CFT_MAX_BODIES` bounds the bodies of one system, not the number of
  systems. The shapes go past every width run so far: entry 37's 24,576
  lanes, and V0's 49,152 lanes, 524,272 pairs and 101,472-element
  gather source.
- **The gate fails when nothing was compared.** On the software
  backend, every bit must be published, since it publishes them all.
  *V0's N4, a probe that reads no bit, passed every check before.*
- Each family has a control that must see a difference. Hidden bits
  count as "not compared", never as passes.
- `--door-refuse` (:2667) makes five calls that both routes must
  refuse, and the gate runs each on both routes.

The checks also cover every refusal by name, including a misspelled
route variable and a misspelled hide variable. A forced-off run gives
the same record as the default, with every counter required on both
sides. That comparison is vacuous until a parcel routes a call, and the
gate says so.

Two limits on the variables are stated, not closed:
- **A typo inside the prefix** makes an unrelated variable, which
  nothing reads, for example `CFT_REBOUND_ROTUE_INDEXED`.
- **On Windows, the library half inside Python** sees only the
  environment as the process started. A variable set through
  `os.environ` is not seen, and `CFT_REBOUND_ARTIFACT` was already like
  that. docs/PYTHON.md, "The Windows rule", has V0's measurement.

Measured on 2026-10-31 at `round3-p0-v1`: 131 of 131 cases identical to
the reference on both routes, at binary64, binary128 and binary256, with
0 accounting errors and 3 of 3 controls. Seventeen sabotaged copies of
the engine all made the gate fail:
- the five the first gate caught: broadcasting element 1, `CFT_IDX_NONE`
  reading element 0, a segment dropping an element, the hide hook
  ignored, and a forced route not refused;
- the five V0 found it passing, with N7 run at both binary128 and
  binary256: N1, N2, N3, N4 and N7/N8;
- the four V0's re-check found it passing: F4, F4 with N8, F6 and F11;
- the two rules the re-check added, each removed: F12 and F13.

On the software backend these test the wrappers. There, each new route
is the same primitive today's route uses. On a card publishing a bit,
they test the tile's route as well.

**Nothing in the integration calls the door yet.** That is what makes
P0 a refactor. The suite and the A/B hold it to "no record line moved";
new comment lines are not changed record lines.

**(d) The trailer's keys, and the decomposition pass, derived.**
- The instrumentation keys print from one table, `trailer_keys[]`
  (:2291), in their old order and formats, plus `t_reduce`.
- *Correction.* Identity is not in the trailer. It is a `# identity`
  header line (:2852) holding `cft_build_id()` and, where the device
  answers, the image's SHA-256. It is a line of its own because
  `tools/check_records.py` reads `# program=` back as the run's
  arguments.
- The four card-day normalisers (`hw/bench-ensemble.sh`,
  `hw/cardtest-f128.sh`, `hw/compare-images.sh` and
  `hw/bench-workload.py`) drop the `# route` and `# identity` lines.
  Both say where a run ran, not what it computed, and both differ
  between a card and the software backend by design.
  - The normalisers also drop `staged_bytes`. *V0 found that it depends
    on the route:* a scalar stages one element on the new route and n
    on today's.
  - The three shell scripts **print each card run's routes** beside
    their verdict. They used to delete those lines unseen.
- **The gate cache keys on every `CFT_REBOUND_` variable, by prefix**
  (`tools/gate_cache.py --env-prefix`). *V0 found a pass with every
  route hidden satisfying the next default run.*
- `hw/bench-workload.py` now writes `decompose.csv` itself: one row per
  port run, with every key its trailer carried as a column, read by
  `trailer()` (:169) and never listed.
- `hw/report.py` (:321) adds a table of every column its own table
  doesn't show. A parcel's new key reaches the bench sheet with no edit
  to either file.
- `bench-workload.py` also runs on Windows now, taking `.exe` as
  `tools/check_program_engine.py` does. That is how this was tested here:
  `--quick` at binary64, 9 rows, 0 failures, 36 columns.

**For P3:** `struct ias15_engine_view` gains `t`, `t_hi` and `t_lo`, the
wide clock read and written in place (`src/ias15_engine.h:55`). The
engine also gains `ias15_engine_build_id()` and
`ias15_engine_image_digest()` (:305).

**(e) The gates, globbed: not done, and why.**
- Only P1 and P2 add drop-in topics, which is below METHOD.md §2's
  three-column bar.
- `main()`'s order carries meaning. `cases_wide()` is a separate mode,
  and `cases_seam()` runs last.
- A glob would have to special-case both.

So round 2's one-line edits stand, and each brief lists them as
expected:
- `CASES_SRC` (Makefile:470);
- one declaration in `tools/dropin_cases.h`;
- one call in `tools/check_dropin.c`'s `main()`;
- for a new Python leg, one line in each of `check` (:587),
  `check-quick` (:612) and `check-light` (:640).

---

## The parcels

Each parcel lists:
- the files it owns, and the files it must not touch;
- what reading has already turned up, and the trap;
- the gate, and the negative control, named by the property that makes
  it bite (CASE-STUDY-4).

A parcel that cannot pass its own gate is not done. Each brief is
written at dispatch from `../ParcelRound/templates/brief.md` and the
parcel's section here. It adds the ledger's path, the exact base
commit, and the scratch directory.

### P1 - Gravity as indexed program runs *(the long pole; dispatch first)*

**Owns:**
- `gravity_body()` ([src/ias15_cft.c:1150](src/ias15_cft.c:1150)) and
  `build_scatter()` (:1053);
- the pair-table part of `alloc_state()` (:1083), plus the rebuild in
  `ias15_engine_set_active()` (:3215) if the tables change shape;
- **the indexed program-run wrapper**, a new function beside
  `run_program()` (:1393), which moved here from P0's plan. It takes its
  route from P0's `route_on[ROUTE_INDEXED]`, and its fallback is the
  host's gather and a dense run. It counts into `route_calls[]` and
  `route_fallbacks[]` like `vrun_ix()` (:572), and is held by its own
  cases in `--door-selftest` (:2478) beside P0's;
- a new gravity section in `tools/gen_programs.py`, with one entry in
  its `main()`;
- one `load_one()` line in `load_programs()` (:1385);
- a new `tools/cases_gravity.c`.

**Forbidden:**
- `pc_error()`, `correct*()` and `step_attempt()` (P2's);
- `predict_next_step_vec()`, `dtnew_*()` and `choose_timestep()` (P4's,
  or nobody's);
- the archive files (P3's);
- P0's probe and wrappers: use them, and report it if they are wrong;
- the integrator-only files.

**What reading turned up. Verify it; do not trust it.**

- **Entry 37's shares, on the card's one f128 tile.**
  - 172,000 six-element calls at 108 µs each are 77% of a Kepler step.
  - Gravity, inclusive, is 74-90% of a step at 64-256 bodies.
  - Divide and square root, which this parcel does not touch, are 39%
    of a step at 64 bodies and 59% at 256 (binary128).
  - So the ceiling, before the gather is paid for, is about 2x at 64
    bodies and 1.36x at 256. *Believed*: this is arithmetic on entry
    37's shares.
- **cft-fp256's measurement.** cft-fp256 read `build_scatter()` as a
  gather by a static table followed by a fixed-order fold from +0, and
  built ABI 0.14's index tables for that shape. On its card at 128
  bodies, one indexed run beat the 127 dense calls it replaces by 1.51x
  to 1.71x, by format, with the 127 host gathers gone as well (its
  docs/ROADMAP.md, ask 1). *That is their measurement, of that piece
  only.*
- **The device gather is not free.** It costs 310-340 ns a gathered
  element on their card at every format. For the pair inputs, the host
  loops at :1153-1157 (`t_g_copy`) may well be cheaper.
- **The pair gather exists twice.** `energy_all()` (:2208) has its own
  copy at :2222-2226, for the record's energy column. It is not gravity
  and is not on the hot path, so it stays as it is, and the energy
  column is part of every gate's record. The card day
  decides.
- **A finding to confirm before relying on it.** The port's fold is
  believed to start from +0 and run in one static order for every body.
  Confirm that, including under `N_active` and `gravity_ignore_terms`,
  where round 2's P2 reproduced REBOUND's `(i, j<i)` loop under the
  skip. If it does not, the indexed fold is a different sum. Design
  from what is there, and say so.

**The trap that would be a silent wrong answer: the fold's order and
its zero.**
- A fold that runs as a tree, or starts from its first addend instead of
  +0, or lets `CFT_IDX_NONE` padding into a sum REBOUND never forms,
  changes bits only where the reordering shows. Plenty of problems hide
  it.
- +0 + -0 is +0. So a fold that starts from the first addend differs
  from REBOUND's `ax = 0; ax += ...` exactly when that addend is -0.

**The routes.**
1. The accumulate half as one indexed run per force evaluation. This is
   the shape with a card measurement behind it.
2. The pair arithmetic as one program over P lanes. Its six inputs are
   gathered **either** by today's host loops **or** by index tables,
   chosen through P0's switch.

The dense route stays, and is the fallback.

**Also priced, not built: a whole-program divide and square root.**
They would give the contract's bits *inside* a program, which could put
gravity's divide and square root inside the run instead of costing a
round trip each.
- The images are `third_party/cft-fp256/programs/divfull-<fmt>.cfta`
  and `sqrtfull-<fmt>.cfta`.
- The route is opt-in through `CFT_DIVSQRT_FULL`
  (`host/src/divsqrt.c:1459` in the pinned tree).
- cft-fp256 measured it slower per element: 2.33 against 1.70 µs, and
  1.185 against 1.01 after sequencer revision 5.

Report the instruction count and what adopting it would take. Do not
adopt it.

**Gate.**
- Every existing gate passes unchanged: REBOUND bit for bit at binary64
  on every problem, the wide formats, and `cases_pairs`.
- A route-equivalence case shows every route gives identical records at
  every format, with test particles. It includes one force evaluation
  at the largest body count the engine claims (`CFT_MAX_BODIES`, 1024,
  [src/ias15_limits.h:29](src/ias15_limits.h:29)), and the refusal one
  past it must still hold. It also asserts the route line.

**Negative controls.**
- *Order.* A table with two partners swapped must change the binary64
  record. Use a problem first shown to be order-sensitive, and say
  which. A control on a problem where the swap is inert proves nothing
  (METHOD.md §5, "inert by construction").
- *Route.* One index pointed at the wrong partner must change the
  answer. That proves the indexed route ran, not the fallback.

**Stop line.** If the library's indexed run cannot express the fold in
the port's order, stop and report. Do not invent a different
summation.

### P2 - The corrector's exit test

**Owns:**
- `pc_error()` ([src/ias15_cft.c:1311](src/ias15_cft.c:1311)) and its
  two call sites, `if (n == 7) pc_error(...)` at :1297 and :1457;
- the corrector's exit test inside `step_attempt()` (:2007), and only
  that;
- the two new flags in `main()`'s parser (:2689);
- one Makefile rule that assembles `normalabs-<fmt>` **from the pinned
  tree** (`third_party/cft-fp256/programs/`), so the pin governs it and
  no copy drifts;
- a new `tools/cases_corrector.c` or check script.

**Forbidden:** gravity (P1), the step control (P4), the archive (P3),
and P0's probe and wrappers.

**Two mechanisms aim at one 4-9%.** That is the share of a card step
cft-fp256 measured for removing the test entirely. On software it
measured nothing at all (its ask 3; ROADMAP item 3).

- **(a) The recorded schedule.** `--record-iters FILE` writes one
  corrector pass count per step attempt. `--replay-iters FILE` runs
  exactly those counts, and skips `pc_error()` at both call sites and
  the exit test. The prototype is `~/patch_fixed_schedule.py` on
  amd-arc-box, dated 2026-09-10 and 123 commits behind; the lead copies
  it into the ledger. ROADMAP says its four anchors still exist.
  *Verify.* The replay's byte-identical record **is** the proof that the
  optimisation is exact, a claim only a bit-reproducible integrator can
  make.
- **(b) The test on the device.** `cft_reduce_seg` (`cft.h:1321`), through
  P0's `vreduce_seg()` (:641), with `seg` = 3N, over the output of
  `normalabs-<fmt>`. cft-fp256 states
  that this pair is `pc_error` exactly in both modes, the +0 identity
  included (its ROADMAP, ask 7). *A finding to confirm* by reading
  `pc_error()`:
  - INDIVIDUAL mode is the maximum of `|tmp/at|`, a divide per
    coordinate before the maximum;
  - the other modes take `max|tmp| / max|at|`;
  - a system that has left the corrector is given 1/1.

  If the pair cannot express INDIVIDUAL's divide per coordinate, that
  mode keeps the host route, and the report says so.

**The trap: a replay that reproduces the state but not the counters.**
- `iterations_max_exceeded` (round 2's P4), `steps_rejected` and the
  pass counts are printed in the trailer, and they must replay too.
- So must an ensemble whose systems leave the corrector at different
  passes.

**Gate.**
- The replay is byte-identical to the recording run at every format,
  and in a mixed-pass-count ensemble. Its own case asserts the
  trailer's counters.
- The device test gives the host route's pass counts and records on
  every problem, in every mode it claims.

**Negative controls.**
- *Replay.* A schedule with one count changed by one must change the
  record, or be refused by name. That shows the replay drives the run.
- *The device route.* There is one case per property, and each must
  fail with its property removed from the route:
  - a system with no normal coordinate (the +0 identity);
  - a subnormal coordinate (normal-only);
  - a system that has left the corrector (1/1).

**Stop line.** If the device route needs per-coordinate host work of
the kind the reduction was meant to remove, stop and report.

### P3 - The archive's last two members, and who made it

**Owns:**
- `src/cft_ias15.h` (append fields only), `src/cft_ias15_fields.h` and
  `src/cft_ias15_fields.c`;
- `src/cft_archive.c`;
- `src/reb_integrator_cft.c`, only where the state is filled for a save
  and adopted on load (`adopt_loaded_state()`,
  [:314](src/reb_integrator_cft.c:314));
- `tests/gate_*.c`, `tools/check_checkpoint.py` and
  `tools/check_archive.py`.

**Forbidden:** `src/ias15_cft.c`, entirely. P0 provides what this parcel
reads:
- `struct ias15_engine_view`'s `t`, `t_hi` and `t_lo`, the wide clock,
  read and written in place (`src/ias15_engine.h:55`);
- `ias15_engine_build_id()` and `ias15_engine_image_digest()` (:305).

If another accessor is needed, ask in the ledger and the lead adds it.

**Three items, one on-disk change.**

1. **The wide clock** (ROADMAP, "What round 2 leaves open: the wide
   clock"). A plain checkpoint differs from the straight run by one
   binary64 ulp in `t` on 2 of 40 step splits, at binary128 and at
   binary256, because the engine's wide `t` is not archived. It
   reproduces on builds from before round 2's P4.
2. **`provenance`** (ROADMAP, "Corrections to the state struct above";
   PARCELS.md, P4 item 3). Nothing in the state records that it was
   promoted from binary64 rather than restored exactly. Only the return
   value of `cft_archive_finish_load()` carries that, and a caller may
   discard it.

   The name now collides with cft-fp256's certificate vocabulary, where
   "Provenance" means the device lines. Propose a name; the lead
   settles it.
3. **Identity.** The library build and the device image, as P0's door
   reads them, go into the archive. That is the "producer says it ran
   on" line, which a certificate's audit reports as stated rather than
   checked (cft-fp256 docs/CERTIFICATES.md, "What a certificate says,
   and what an audit proves").

`CFT_N_SCALARS` and `CFT_N_BLOBS` stay derived, not transcribed (round
2's rule).

**Old archives still load.** The Simulationarchive tolerates strangers
(ROADMAP, "The Simulationarchive is extensible, by name, and tolerates
strangers"). A resumed run whose file has no wide clock behaves as it
does today. Whether it says so is the parcel's call, recorded with its
reason.

**The trap.** A checkpoint case that compares `t`'s binary64 view is
vacuous through the observable (METHOD.md §5). Compare the wide bytes.

**Gate.**
- The 2-of-40 case at binary128 and binary256 is **watched failing** on
  the tree without the fix, then passes with it.
- `provenance` is asserted on both load paths.
- A software-backend archive carries the build id and records the image
  as refused. A card-written archive is the card day's.

**Negative control.** With the new field's load removed, the 2-of-40
case fails again.

### P4 - The step control as a program *(the first to drop if the day runs short)*

This is ROADMAP item 1, second bullet. The design names
`predict_next_step`'s sixteen compensated additions per coordinate as
"one more program over the scratch block", but they are still
elementwise calls. The work is mechanical, in the shape of `correct_n`.

**Owns:**
- `predict_next_step_vec()`
  ([src/ias15_cft.c:1467](src/ias15_cft.c:1467));
- a new section in `tools/gen_programs.py`, with one entry in its
  `main()`;
- one `load_one()` line;
- its case in `tools/check_program_engine.py`.

**Shares with P1** the two one-line entries, in `gen_programs.py`'s
`main()` and in `load_programs()`. The functions are disjoint, so the
conflict is one line, and the lead resolves it.

**Forbidden:** `dtnew_*()` and `choose_timestep()`, which are nobody's
this round, so leave them; and everything P1, P2 and P3 own.

**Gate.** The program route's records equal the elementwise route's at
every format, through `check_program_engine.py`'s existing comparison.

**Negative control.** A program with one compensated addition dropped
must change the record. Find a problem where it does first.

### P5 - `scalar_mask` at the broadcast sites *(second wave)*

This is ROADMAP item 2. As of 2026-09-14 there were 24 call sites that
pass a broadcast constant vector (`KRINV`, `KHF`, `KC`, `KY4`, `KY5`,
`KNUM`, `KVNUM`, `KPOSR`, `K5040`, ...), each staging N3 copies of one
value. Recount them, and use P0's `vrun_k()`
([src/ias15_cft.c:506](src/ias15_cft.c:506)), taking its `DOOR_UNUSED`
off when the first call lands.

The sites run through functions P1, P2 and P4 own, so **P5 is
dispatched after those three merge.**

**Value.** A few percent of a card step (entry 37). Do it, as ROADMAP
says, "for the memory and the honesty of having asked".

**Gate.** Records are unchanged. The contract says a scalar operand
applies element 0 to the whole run, so the bits must not move.

**Negative control.** `staged_bytes` in the trailer drops by an amount
**derived** from the recounted sites and their widths, not transcribed,
while the bits do not move.

---

## Rules for every parcel

- Work in a **worktree**, in the background. **Never push, merge or
  rebase**; the integrator does all three. Commit in your worktree
  only.
- **Check your base** with `git log --oneline -1`, against the commit
  in your dispatch message.
- Build in your own tree; `build[0-9]*/` is gitignored for exactly
  this. On Windows the PATH line is not optional:

  ```bash
  export PATH="/c/msys64/mingw64/bin:$PATH"
  make OS=Windows_NT TMP="$TEMP" TEMP="$TEMP" \
       PYTHON="/c/Users/logan/AppData/Local/Programs/Miniconda3/python.exe" <target>
  ```

  Without it, bare `cc` is the 32-bit toolchain, and the link fails on a
  good object with `file format not recognized`.
- Run `make check-light` while you work. Run a full `make check` (about
  70 minutes here) in the background before you report. **Read the log,
  not the exit code.**
- **A new leg names every artifact it executes as an `--input`.** The
  gate cache keys on executed artifacts, never on source
  (`tools/gate_cache.py`). A missing input lets `check-quick` hold a
  stale pass across your change.
- **Follow the `--light` scaling rules.**
  - Cases that integrate to a time stay unscaled.
  - So do step-sequence cases.
  - A case that must not be scaled opts out by name, through
    `(reb_simulation_steps)(...)`.
  - A new case says which kind it is (CLAUDE.md).
- **Bit-identity is not validity.** Print the physics counter beside
  the comparison and assert it (CLAUDE.md, "Three ways this suite has
  lied").
- **Nothing compares the trailer for you.** Whatever you rely on there,
  a route line, a counter or a byte count, your own case asserts.
- **Derive; never transcribe.**
- **A gate that cannot fail is not a gate.** Run your named control,
  confirm it fails, delete its artifact, and report both results.
- **READY is "a gate, or a stated limit"**, and a limit is stated by the
  behaviour it concedes (CASE-STUDY-4).
- **The ledger.** Read every file in it at three points: before you
  start, before you design anything this plan calls shared or
  forbidden, and before you report. Arm the `urgent/` watch first.
  Append to your own file only.
- **Label every claim** as measured or believed.
- **Never run bare `git stash`.** The stash stack is shared across
  worktrees.
- **Read `../cft-fp256`; never write it.** Touch nothing in
  `../ParcelRound`.
- **No card.** amd-arc-box and the U50 belong to the lead and the
  owner, on the card day.
- **No speed claims.** No parcel makes one; the card day measures.
- **Report:**
  - what changed, file by file;
  - the gate output, as actual lines;
  - the controls, with their numbers;
  - what you did not do;
  - **anything you found that this plan got wrong.**

## What the integrator keeps

- **P0**, before anything is dispatched.
- **The seam tests.** They belong to no parcel. Each merge gets one
  that exercises two parcels together:
  - a replayed schedule over the indexed gravity route (P1 × P2);
  - a checkpoint taken on the indexed route and resumed on the dense
    one, bit for bit (P1 × P3);
  - a replay across a checkpoint (P2 × P3);
  - the step-control program over indexed gravity (P1 × P4).
- **The integrator-only files:**
  - `README.md`, `ROADMAP.md` and `CLAUDE.md`;
  - `docs/VALIDATION.md` and `docs/HARDWARE.md`;
  - `PARCELS.md` and this file;
  - `third_party/MANIFEST`.

  A parcel that needs one changed says so in its report.
- **Every merge, serially, with the full `make check` after each one.**
  At about 70 minutes a run, four merges are close to five hours of
  suite. P3 and P4 share no file, so the lead may merge both and run one
  suite, and will say so (METHOD.md §8). Never merge while a suite is
  running.
- **A verifier per parcel:** V1-V4, and V5 for the second wave. Each
  verifier:
  - re-runs the gate from a clean build;
  - re-runs the control;
  - diffs the change against the ownership lists above;
  - fixes nothing.

  "Found nothing" is an acceptable answer.
- **The VALIDATION entries**, from 39. Parcels report; the lead writes.
- **The docs sweep, in idle time.**
  - HARDWARE.md's sentences that call a mechanism missing now mean
    "missing at the old pin" (its 2026-10-30 note).
  - CLAUDE.md's tile advice waits for the card day's `CFT_XRT_TILES`
    row.
- **The summary to the owner**, with the cost per agent and the
  verifiers' share.

## Order of the day

| when | who | what |
|---|---|---|
| before dispatch | owner, lead, V0 | Everything in the table above: the pin merged, the Linux leg run, P0 built, verified and pushed. |
| dispatch | lead | **P1 first** (the long pole), then P2, P3 and P4. The lead arms its own `urgent/` watch with the first dispatch. |
| as each returns | V1-V4 | Verify. A send-back resumes the parcel with its context, so it needs no re-brief. |
| as each lands | lead | Merge, run the seam test and the full `make check`, and read the log. |
| second wave | P5, V5 | After P1, P2 and P4 are merged. |
| the card day | lead, owner | See below. |
| last | lead | The VALIDATION entries, then ROADMAP, README and CLAUDE.md. |

If the day runs short, **P4 is dropped first**: a few dozen calls a
step, against gravity's hundreds. P5 goes next: it is worth a few
percent, and its sites move under every other parcel.

## The card day

This comes after the build phase, on a day the owner books outside
cft-fp256's runs. It runs against the **one-tile** f128 image, since on
a card tile count is a cost (CLAUDE.md). Every card row is held
byte-identical to its software counterpart.

1. **The pin on the card.** Run the suite against the image through the
   new library. This is entry 38's missing row, and it settles whether
   0.18 opens the image.
2. **The gravity split (ROADMAP item 4), as the "before".**
   - Build P0's commit on amd-arc-box beside the round's tip.
   - Run Kepler, the outer solar system, the Pythagorean problem, and
     n-body at 64, 256 and 512 bodies, at binary64 and binary128.
3. **P1's routes.**
   - The accumulate fold, indexed against dense.
   - The pair gather, host against device.

   Run both at 64, 128, 256 and 512 bodies. The split says which
   gather wins, and by how much.
4. **P2:** the device test against the host loop, and a replay against
   its recording.
5. **P4's program and P5's staged bytes**, if merged.
6. **`CFT_XRT_TILES`:** a quad opened as one tile, against the one-tile
   image. If they are equal, CLAUDE.md's advice gains "or the quad
   with one tile opened".
7. **Still owed from before this round:**
   - the ensemble gate with an artifact;
   - the mixed-pass-count ensemble that prices idle lanes.

   Both are docs/HARDWARE.md, "What to do first", items 1 and 4.
8. **P3's identity on a card-written archive:** the image digest is
   present.

## Not in this round, and why

- **Residency (ROADMAP item 5).** Entry 37 prices it at 1-10% of a step
  and puts it after items 1 and 3, which are this round.
- **The divide.** cft-fp256 built the ask and measured it slower. The
  number moves with its sequencer RTL, not with anything here.
- **Certifying an IAS15 run.**
  - cft-fp256's audit re-runs segments of one program image. By its own
    stated limits, a route whose divide and square root are host calls
    between program runs cannot be certified (docs/CERTIFICATES.md,
    "What version 1 does not do"). IAS15 here has that shape, and more.
  - The natural segment for this port is the interval between two
    archive snapshots, re-run on another conforming build. Entry 38's
    cross-build checkpoint check is already that, 12 of 12. It needs
    P3's complete archive first, and then a certificate kind that
    cft-fp256 would have to define. That would be a ninth ask, the
    owner's to file.
  - At binary64, REBOUND itself is an independent auditor. Above
    binary64, the only independent auditor would be an IAS15 in
    cft-fp256's golden model.
- **Bitstreams and image builds.** Arithmetic is about a fifth of card
  wall clock, and an 11% faster engine moved the integration 2%
  (CLAUDE.md).
- **cft-fp256's `exp2m1` fix** (`48e2aa22`, after the pin). This port
  never calls `exp2m1`.
- **Two defects on the drop-in's refusal path.** V0 found both on
  2026-10-31. They were there before P0: pre-P0 did the same for a
  missing artifact.
  - The refusal sentence prints twice, and `steps_done` becomes 1.
  - From Python, REBOUND's message store cuts the sentence at 254
    characters. The door's sentences put what matters first for that
    reason.

  The path is `src/reb_integrator_cft.c`'s, which no parcel owns this
  round, so both are filed for the next round.
- **Variational particles, MEGNO, WHFast, ghost boxes and the tree
  code.** Unchanged from round 2's list (PARCELS.md).
