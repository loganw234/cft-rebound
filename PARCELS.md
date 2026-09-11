<!-- SPDX-License-Identifier: GPL-3.0-or-later
     Copyright 2026 the cft-rebound contributors. -->

# Round 2: from a correct port to one a stranger can use

Round 1 (ROADMAP.md, "The work, in four parcels") made the drop-in
exist and made it bit-identical. This round is about the gap between
*correct* and *usable*: the refusals a real simulation actually hits.

**The target.** Someone with an existing REBOUND script points it at
`ias15_cft`, and either it runs or it refuses with a reason they can
act on. Today a large fraction of real scripts refuse, and the reasons
are not exotic — they are REBOUNDx and massless test particles.

Ranked by how often a real script hits them:

| refusal | who hits it |
|---|---|
| `additional_forces`, `pre_`/`post_timestep_modifications` | **REBOUNDx**: GR corrections, tides, migration, radiation forces. The single largest blocker |
| `N_active` (test particles) | Anything with massless particles — most planetesimal and debris work |
| `REB_GRAVITY_COMPENSATED` | People who already care about round-off, which is precisely this project's audience |
| `adaptive_mode` 0 and 1 | Narrower since PRS23 became the default in January 2024, but the port claims IAS15 and covers two of four |
| variational particles, MEGNO | A real community, but a second particle set is a round of its own |

Plus two things that are not refusals and still stop people: **Python
on Windows does not work**, and **nothing has ever run on macOS**.

---

## The rule this round is built around

ROADMAP.md's last section is the post-mortem of round 1, and its
finding was not "a parcel was wrong". Every parcel was right. The
finding was:

> Where a shared fact appeared twice, the copies drifted.

The blob count was hand-written in four places and the state struct in
two, and growing the list broke all six in different ways.

**This round has exactly the same shape of shared fact, and more
copies of it.** "What is supported" is currently stated in four
places:

1. `supported()` in [src/reb_integrator_cft.c:126](src/reb_integrator_cft.c:126) — 16 refusals
2. `cft_rebound_check()` in [src/cft_rebound_run.c](src/cft_rebound_run.c) — a deliberately different list
3. `case_refusals()` in [tools/check_dropin.c](tools/check_dropin.c)
4. The support table in README.md

Every parcel below removes a refusal, so without a fix all five would
edit all four files, and the drift would be guaranteed rather than
merely likely.

So **P0 comes first, is mine, and lands before any parcel starts.**

---

## P0 — one support table — **DONE**

Landed before dispatch, with the full suite green on both sides of it.
What a parcel now finds:

**(a) The refusal list is one table.** [src/cft_supported.h](src/cft_supported.h)
declares it and [src/cft_supported.c](src/cft_supported.c) holds the 21
rows: a name, a mask of which paths refuse it, a predicate, a message.
`supported()` and `cft_rebound_check()` are each now a context fill plus
one call to `cft_support_first_refusal()` — the loop exists once, and so
does every message.

**Unlocking a capability is deleting one row.** Nothing else in `src/`
states the list.

`case_refusals()` keeps its own poison functions, but every
`refused()` / `accepted()` call names the row it exercises, and
`coverage_report()` walks the table afterwards and **fails by name** on
any drop-in row no case named — and on any case naming a row that is
gone, which is what catches a test left behind after a refusal is
removed. Verified with a negative control: a dummy row makes the gate
fail with `coverage: cft_support_rows has "control_row" and no case
exercises it`.

It earned its keep on the first run. Four rows — the particle map, the
ensemble count, the format and `max_iter` — had never been poisoned by
anything, and are now.

**(b) The gate is a file per topic.**

| file | what |
|---|---|
| [tools/dropin_cases.h](tools/dropin_cases.h) | the scaffolding, declared once, and the instructions for adding a topic |
| [tools/dropin_common.c](tools/dropin_common.c) | the problems, `build()`, `compare()` |
| [tools/cases_core.c](tools/cases_core.c) | the binary64 equivalence set |
| [tools/cases_wide.c](tools/cases_wide.c) | the binary128 smoke tests |
| [tools/check_dropin.c](tools/check_dropin.c) | `main()`, the refusals, the coverage check |

**A parcel writes `tools/cases_<topic>.c`, declares one entry point in
`dropin_cases.h`, adds it to `CASES_SRC` in the Makefile, and adds one
call in `main()`.** Four small edits outside its own file, none of them
in another parcel's territory.

One more shared fact collapsed on the way past: the body-count cap was
1024 in two files with a comment saying they were "kept in step". It is
[src/ias15_limits.h](src/ias15_limits.h) now.

---

## The parcels

Each one is: owned files, forbidden files, what it unlocks, and the
gate it must pass. A parcel that cannot pass its own gate is not done —
ROADMAP.md line 286, and it still holds.

### P1 — Forces at every substage *(largest; the REBOUNDx unlock)*

**Owns:** `step_attempt()`'s substage force path in
[src/ias15_cft.c:975](src/ias15_cft.c:975) (the two `gravity()` call
sites, lines 977 and 1035); a new `ias15_engine_set_force_hook()` in
`src/ias15_engine.h`; a new `src/cft_forces.c` for the shim-side
plumbing.

**Forbidden:** `gravity()` itself (P2 owns it), `choose_timestep()`
(P3), the archive files (P4), anything under `python/` (P5).

**Unlocks:** `r->additional_forces`, `r->force_is_velocity_dependent`,
`r->pre_timestep_modifications`, `r->post_timestep_modifications`.

**Scope is smaller than it looks — two of the four already work.**
`pre_` and `post_timestep_modifications` are called by REBOUND's
*driver*, not by IAS15: `simulation.c:523` and `:564`, around the step
callback. They edit `r->particles` between steps, and the shim already
re-promotes any coordinate that changed under it. So they need a **gate
case, not a mechanism** — and the row that refuses them is on the
subprocess path only, so there is nothing to delete for them either.

The mechanism is for `additional_forces`, and velocity-dependent forces
come with it for free: IAS15's predictor supplies `v` at every node, so
that refusal is inherited rather than required.

**The shape.** `reb_simulation_update_acceleration()` (`simulation.c:643`)
is gravity followed by `r->additional_forces(r)`, and IAS15 calls it at
every Gauss-Radau node — `integrator_ias15.c:461` — then reads
`particles[mk].ax/ay/az` straight back into `at[]` at `:467`. So the
port needs, after each wide `gravity()`: round the node's `x`, `v` and
`a` out to `r->particles`, call the routine, promote `ax/ay/az` back.
Promotion from binary64 is exact, so the total is wide gravity plus an
exact binary64 addend, correctly rounded.

**`r->t` moves at every node, and this is the part to get right.**
`integrator_ias15.c:408` sets `r->t = t_beginning + r->dt * h[n]` before
each call and restores `t_beginning` at `:608`. A time-dependent force —
most of REBOUNDx — reads it. For bit-identity the shim must compute that
same binary64 expression from REBOUND's own `h[n]`, **not** round the
wide clock: the two differ in the last bits, and a force that sees a
different time computes a different acceleration. The port's derived
`KH[8]` rounds to REBOUND's `h` literals bit for bit at binary64
(`tools/gen_constants.py` checks all 78), so exposing `h[n]` as a double
from the engine is enough.

**The engine may not include `rebound.h`.** `src/ias15_cft.c` builds
both as the standalone program and, under `-DIAS15_CFT_LIBRARY`, as a
library half that knows nothing about REBOUND — that separation is load
bearing and must survive. So the hook is a plain callback with the
rounded buffers, registered by the shim; the engine never sees a
`struct reb_simulation`. `ias15_engine_remove_body()` and
`ias15_engine_alias_resize()` are the shape to follow.

**Gate** — a new `tools/cases_forces.c`, wired as P0's table above describes:
- a constant additional acceleration, against REBOUND's own ias15 with
  the same routine — bit for bit at binary64;
- a **time-dependent** force, which is the case that catches a wrong
  `r->t` at the nodes and which nothing else here would;
- a **velocity-dependent** drag force, same assertion, which is what
  retires that refusal rather than merely deleting it;
- a `post_timestep_modifications` routine that edits a coordinate;
- **the no-force path is unchanged**: every existing case in the gate
  still passes, which is the regression this parcel is most likely to
  cause — 8 round-trips per attempt land on the hot path, and the
  no-force branch must not pay for them.

**Negative control.** A force routine that does nothing must leave the
run bit-identical to no force routine at all; and the constant-force
case must differ from the no-force run, or it is a pass against itself.

**State in the brief:** the user's force is evaluated at **binary64**
even in a binary256 run, because `r->particles` are binary64 and that is
the only interface REBOUND offers a callback. Gravity stays wide. This
is a real ceiling and the parcel must document it, not discover it.

### P2 — Which pairs are computed

**Owns:** `gravity()` at [src/ias15_cft.c:549](src/ias15_cft.c:549);
the pair-list construction in `alloc_state()` and
`ias15_engine_set_bodies()`; a new `ias15_engine_set_active()`.

**Forbidden:** `step_attempt()` (P1), `choose_timestep()` (P3), the
archive files (P4).

**Unlocks:** `r->N_active` (test particles), `r->gravity_ignore_terms`,
`r->map`.

**The trap.** `N_active` does not merely zero a mass — REBOUND *skips*
those pair terms, so both the set of addends and their summation order
change. Bit-identity requires reproducing its `(i, j<i)` loop under the
skip, not reproducing the sum. The existing ascending-partner scatter
at line 576 is where that lands.

**Gate** — a new `tools/cases_pairs.c`: test particles bit-identical at
binary64 on at least two problems; each `gravity_ignore_terms` value;
a particle map. Plus the control this repo has needed twice before —
**prove the skip actually changed something**, or the case is a
pass against itself.

### P3 — The other two step criteria

**Owns:** `choose_timestep()` at
[src/ias15_cft.c:869](src/ias15_cft.c:869), and nothing else.

**Forbidden:** everything else. This is the cleanest seam in the round.

**Unlocks:** `adaptive_mode` 0 (INDIVIDUAL) and 1 (GLOBAL), which take
REBOUND's other error estimate entirely. PRS23 (2) and AARSETH85 (3)
are already there and must not move.

**Gate** — a new `tools/cases_modes.c`: both modes bit-identical against
REBOUND on the existing problem set — **and a control that each mode
chooses a different step sequence from PRS23 on the same problem.**
AARSETH85 needed exactly this control: a mode the port silently ignored
would leave every case passing PRS23 against PRS23.

### P4 — State completeness across a checkpoint

**Owns:** `src/cft_ias15.h` (append fields only), `src/cft_ias15_fields.h`,
`src/cft_ias15_fields.c`, `src/cft_archive.c`, `tests/gate_*.c`,
`tools/check_checkpoint.py`.

**Forbidden:** `src/ias15_cft.c` except for read-only accessors it may
*request* from the integrator; the shim's `supported()`.

**Three known gaps, all recorded:**

1. **`iterations_max_exceeded` has no home in the state.** REBOUND
   archives it as `REB_UINT64` and fires its "at least 10 predictor
   corrector loops did not converge" warning off it. The engine counts
   it (`ias15_engine_max_exceeded()`) and the shim drops it, so the
   warning never fires. Accumulate as a delta into the state — the
   engine's counter is process-global and `reset_state()` zeroes it,
   REBOUND's is per-simulation and monotonic.
2. **The checkpoint corner opened on 2026-09-11** (ROADMAP.md, "What
   the collision work leaves open"). REBOUND archives its whole flat
   coefficient buffer at `N_allocated`, so the stranded tail and the
   mark both survive a save; this port archives each level at `3N` and
   has no `N_allocated`. Remove → save → resume → regrow diverges.
3. **`provenance`** — ROADMAP.md line 404 says it is still owed:
   nothing in the state records that it was promoted from binary64
   rather than restored exactly, and only `cft_archive_finish_load()`'s
   return value carries that, which a caller may discard.

**Gate:** extend `tests/gate_real.c` with the remove → save → resume →
regrow sequence, which **must fail before the fix**. `CFT_N_SCALARS`
and `CFT_N_BLOBS` are derived, not transcribed — keep it that way.

### P5 — Python on Windows, and packaging

**Owns:** `python/`, `docs/PYTHON.md`, the Makefile's `python-lib` /
`check-python` / `install` targets, and a `pyproject.toml` if it
earns one.

**Forbidden:** every `src/` file. This parcel is pure plumbing and must
stay that way.

**The known blocker,** from docs/PYTHON.md line 184: a Windows DLL may
not carry undefined symbols, so it cannot leave REBOUND's symbols to
the loader the way the POSIX `.so` does — it has to link
`rebound.__libpath__`, which ties the build to a Python minor version.
That is the design question to answer, not a detail to code around.

**Gate:** `make check-python` on Windows, getting the same 13 values
bit for bit that the Linux leg gets. Plus `check-rebound-match` still
refusing a mismatched wheel.

---

## What the integrator keeps

- **P0**, before anything is dispatched.
- **The seam tests.** ROADMAP.md's rule, learned the expensive way:
  they belong to no parcel, so they belong to the integrator. Each
  merge gets a test that exercises *two* parcels together — forces plus
  test particles, a mode change across a checkpoint.
- **The integrator-only files**, which no parcel may touch: `README.md`,
  `docs/COMPATIBILITY.md`, `docs/COMPLIANCE.md`, `docs/BENCHMARKS.md`,
  `docs/SCALING.md`, `.github/workflows/gates.yml`.
- **Every merge, and the full `make check` after each one** — including
  the fp128 and fp256 legs, which parcels are not expected to run.
- **The VALIDATION entries.**

## Order of the day

| when | who | what |
|---|---|---|
| ~~first~~ | integrator | ~~P0~~ **done**, suite green before and after |
| then, in parallel | P1 · P2 · P3 · P4 · P5 | five worktrees, no shared owned file |
| as each lands | integrator | merge, seam test, full suite |
| last | integrator | README support table, COMPATIBILITY, VALIDATION |

P1 is the long pole and should be dispatched first of the five. P3 and
P4 are the most likely to finish early. If the day runs short, **P3 and
P5 are the ones to drop** — P3 because PRS23 is REBOUND's default and
covers most users, P5 because the POSIX path works and Windows Python
is a packaging problem rather than a correctness one.

## Rules for every parcel

- Work in a **worktree or clone**, in the background. **Never push,
  never merge** — the integrator does both.
- Build in your own tree. `build[0-9]*/` is gitignored for exactly
  this; do not build into another tree's `build/` while a suite runs
  there.
- **Derive, never transcribe.** Counts and name lists come from the
  definition that owns them.
- **A gate that cannot fail is not a gate.** Every parcel above names
  its negative control; run it, then delete the control binary.
- Touch nothing under `atlas-engine` or `cft-fp256`. Different repos,
  active work.

## Not in this round, and why

- **Variational particles and MEGNO.** A second particle set through
  the same integrator is a round of its own, and the audience is
  narrower than REBOUNDx's.
- **WHFast.** A different integrator, not a gap in this one. It is
  ranked first to follow in docs/INTEGRATORS.md and stays there.
- **Ghost boxes, periodic and shear boundaries, the tree code.**
  Genuinely out of scope for a direct-summation port, and the refusals
  are correct rather than incomplete.
- **The cft-fp256 threads** — sequencer revision 4's remaining layers
  (RTL, host, `CAPS2[6]`, fuzz corpus, docs), the scratch capacity
  increase, and the ordered segmented accumulate. Different repo, zero
  file overlap with anything above; it can run as a sixth lane, but it
  does not make cft-rebound more usable and so is not part of this day.
