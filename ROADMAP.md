# Roadmap: from a correct experiment to a usable integrator

Written 2026-09-10, after the port ran on the FP256 tile.

> **STATUS: all four parcels landed the same day this was written.**
> Read the body below as the plan and the reasoning, not as a statement
> of what is missing - every "will", "should" and "must" in parcels A
> through D describes work that has since been done. What is current:
>
> | the plan said | where it is now |
> |---|---|
> | A. the integrator shim | `src/reb_integrator_cft.c`, `src/cft_ias15.h`, `-lcft_ias15`; gated by `tools/check_dropin.c` |
> | B. Simulationarchive support | `src/cft_archive.c`, `src/cft_ias15_fields.c`; gated by `tests/gate_{restart,write,stock,promote,real}.c` |
> | C. scope, limits, packaging | the README's scope table, the 1024-body cap, `examples/roundtrip.c`, `make install` |
> | D. Python | `make python-lib`, `make check-python`, `python/cft_rebound.py`, `cft_ias15_configure()` - closed, except Windows |
>
> One thing in the body is **still open** and is marked where it
> appears: the `provenance` member (see "Corrections to the state
> struct"). `max_iter` by format in the shim, listed here as open until
> 2026-09-14, is not: `reb_integrator_cft.c` sets it from
> `cft_ias15_default_max_iter(format)` when the caller leaves it at 0
> (see "A footgun found by building on it" for why it once was not).
> The last section, "What integration found",
> is the newest text in this file and is the one to read first if you
> only read one. docs/VALIDATION.md entries 25 through 28 are the
> measurements.

The hard part is done and it is the part that is normally hard:
**correctness**. At binary64 this is REBOUND's own IAS15 bit for bit;
on the card every record is bit-identical to the software backend,
including an ensemble reproducing its members run one at a time. What
is missing is not accuracy. It is that nobody can use it.

Today `ias15_cft` is a standalone program that reads a bespoke problem
file and writes a bespoke record. A REBOUND user has a
`struct reb_simulation` in C, or a `rebound.Simulation` in Python, and
there is no path from one to the other. That is the gap.

*(Those two paragraphs are the problem statement this document was
written to answer. There are three paths now: the drop-in, the
subprocess API, and Python. The README's "Using it from your own
REBOUND program" is the current account.)*

---

## What the investigation found, and it is good news twice

Two questions decided the shape of this roadmap, and REBOUND answers
both itself. **Neither needs a rewrite, and neither needs a fork.**

### 1. REBOUND accepts user-provided integrators, officially

```c
REB_API void reb_integrator_register(const struct reb_integrator integrator, const char* name);
REB_API void* reb_simulation_set_integrator(struct reb_simulation* r, const char* name);
```

`struct reb_integrator` (rebound.h) takes `step`, `synchronize`,
`create`, `free`, the particle add/remove hooks, **and** a
`field_descriptor_list` described in its own comment as "Information on
how to safe/load the integrator state from restart files for this
integrator."

So the drop-in is not a patch to REBOUND. It is a registration. A user
writes their simulation exactly as they do now and changes one line:

```c
reb_integrator_register(cft_ias15_integrator, "ias15_cft");
reb_simulation_set_integrator(r, "ias15_cft");
```

Everything else in REBOUND - particles, outputs, callbacks, the
archive - keeps working, because we are not replacing any of it.

### 2. The Simulationarchive is extensible, by name, and tolerates strangers

The on-disk record is:

```c
struct reb_binarydata_field {
    uint64_t size_name;     // including the \0
    uint64_t size_data;
};
```

Fields are identified **by name**, not by a numeric id, and a reader
that meets a name it does not know does this (binarydata.c):

```c
case REB_FIELD_NOT_FOUND:
    *warnings |= REB_BINARYDATA_WARNING_FIELD_UNKNOWN;
    int err = fseek(inf, field.size_data, SEEK_CUR);
```

It warns, seeks past the field - **and then `goto finish_fields`,
abandoning the rest of the snapshot.** That last clause was missed when
this roadmap was first written, and two parcels found it independently
by running it. What follows is the corrected account.

**What is true.** A stock REBOUND reader opens a cft archive and
recovers every particle coordinate, mass, `t` and `dt` **bit for bit**
(parcel B verified it, on one snapshot and on the last of three). Names
are unique strings rather than indices, so a `cft_` field cannot
collide with any present or future REBOUND field.

**What is not.** Three corrections, all of which change what we should
promise:

- The binary64 state survives **incidentally**, because REBOUND writes
  every simulation field before any integrator field. It is not a
  designed skip and it is not robust to upstream reordering.
- A stock reader never reaches the unknown-field path at all.
  Integrator fields are resolved through the integrator's NAME first
  (binarydata.c:152), so an unregistered `ias15_cft` produces
  `Error! Integrator not found.` and `WARNING_CORRUPTFILE` - not
  `WARNING_FIELD_UNKNOWN`, which for this case is unreachable, as is
  `WARNING_CUSTOM_INTEGRATOR`, which nothing in REBOUND ever sets.
  This is structural: any custom integrator name does it, and there is
  no writer-side fix short of forking REBOUND.
- A stock reader therefore restarts from a **cold** IAS15 series, not a
  bit-identical continuation.

Direct compatibility and native wide precision are still the same file,
and that remains the right design. But the barrier to entry is lowered
less than first claimed: a collaborator without cft gets a usable
simulation AND is told a perfectly good file seems to be corrupted.
**That sentence belongs in the README beside the compatibility claim**,
not in a footnote.

### The one real constraint, stated plainly

`struct reb_particle` holds `double x, y, z, vx, vy, vz`. REBOUND's
particle array is binary64 and will remain so. So a wide-format
integration cannot live in it. The division is:

- **REBOUND's particles are the binary64 VIEW.** Correctly rounded from
  the wide state after each step, so every REBOUND output, callback and
  visualisation works unchanged.
- **The wide state is the integrator's own**, carried in its state
  struct and archived through the `cft_` fields.

A run therefore restarts from an archive *exactly* if the reader
understands the `cft_` fields, and *approximately* - at binary64, which
is what REBOUND itself would have given - if it does not. Both are
correct behaviours and the file says which you got.

---

## The state struct both parcels build against

Defined here so the integrator and the archive can be written in
parallel. `W` is the wide element width in bytes: 8, 16 or 32.

```c
struct cft_ias15_state {
    /* configuration, plain doubles, archived as REB_DOUBLE */
    double   epsilon;          /* as REBOUND's */
    double   min_dt;
    int      adaptive_mode;
    int      format;           /* CFT_FP64 | CFT_FP128 | CFT_FP256 */
    int      max_iter;         /* 12 in REBOUND; binary256 needs ~22 */
    int      arith_fma;        /* 0 = REBOUND's roundings, 1 = the FMA form */

    /* the wide state: byte blobs, archived as REB_POINTER with
     * element_size = W. Lengths are 3N or 3N*E elements. */
    unsigned char *x, *v;             /* the LIVE coordinates: what a restart
                                       * continues from. Missing from the first
                                       * version of this struct - see the note
                                       * at the end of this file. */
    unsigned char *x0, *v0, *a0;      /* the step's starting copy, refreshed
                                       * from x/v at the top of every step */
    unsigned char *csx, *csv, *csa0;  /* compensated-summation carries */
    unsigned char *g[7], *b[7], *e[7], *br[7], *er[7], *csb[7];
    size_t   n_elem;           /* 3N, or 3NE for an ensemble */
    size_t   E;                /* ensemble members, 1 for a plain sim */

    /* provenance, so a reader knows what produced the wide bytes */
    char     cft_abi[16];      /* e.g. "0.11" */
    uint64_t constants_digest; /* the Gauss-Radau set the run used */
};
```

The field names mirror REBOUND's own IAS15 descriptor list - `at`,
`x0`, `v0`, `a0`, `csx`, `csv`, `csa0`, `g`, `b`, `csb`, `e`, `br`,
`er` - so the mapping is one to one and a reviewer can check it against
`integrator_ias15.c` line by line.

---

## The work, in four parcels

**All four landed on 2026-09-10.** Kept as written, because the last
section of this file is about what this way of splitting the work cost,
and that argument needs the briefs it is arguing about. Each is
independently reviewable and they share only the struct above - which
is exactly the problem, as it turned out.

### A. The REBOUND integrator shim (the drop-in)

Implement `struct reb_integrator` and register it.

- `create` / `free`: allocate the state for `r->N` particles.
- `step`: promote `r->particles` to the wide format (exact - binary64
  is a subset of every wider one), take one IAS15 step through the
  existing engine, round the result back into `r->particles`, and
  advance `r->t` and `r->dt` as REBOUND expects.
- `did_add_particle` / `will_remove_particle`: resize and invalidate -
  and, below the high-water mark, re-read the coefficient levels at
  the new stride, which is what REBOUND does rather than what it was
  first thought to do (see the upstream-defects section).
- Honour `r->exact_finish_time`, and `synchronize` where it applies.

**Gate:** a C program that builds a simulation, registers this
integrator at binary64, integrates, and gets bit-identical particle
values to the same program using REBOUND's own `"ias15"`. That is the
existing equivalence gate seen from REBOUND's side, and it is the one
result that proves the shim did not change the arithmetic.

Deliberately out of scope: additional forces, ghost boxes, the tree
code, variational particles. Detect them and refuse with a clear
message rather than computing something wrong.

Since settled the other way, and each with a gate: non-zero
`softening`, `min_dt`, the AARSETH85 step criterion, and **collision
detection** - the driver runs the search and the resolver between
steps, so the integrator's whole obligation is to survive the
removal, which is reproducible once you know what REBOUND actually
does across one. Above binary64 it needs `state->accurate = 1` and
is then deliberately not bit-identical; the README's note under the
support table has the rule.

### B. Simulationarchive support

- Fill `field_descriptor_list` for the configuration and, as
  `REB_POINTER` with `element_size = W`, the wide arrays.
- Prefix every added name `cft_`.
- On load: if the `cft_` fields are present and the format and element
  count agree, restore exactly; otherwise promote REBOUND's binary64
  state and say so.
- Refuse a mismatch loudly - an archive whose `cft_format` is binary256
  loaded into a binary64 run is a user error, not something to paper
  over.

**Gates, all short:** write an archive mid-run, restart, and continue
bit-identically to an uninterrupted run; open a cft archive in stock
REBOUND and confirm it reads with the unknown-field warning and the
right binary64 state; open a stock archive here and confirm promotion.

### C. Scope, limits and packaging

- Raise the 64-body cap. **CORRECTED after parcel C built it:** the
  reason given here was wrong. `body_names` is already
  `static char (*body_names)[32]`, `calloc`'d from the body count, and
  so is `sys_names`. Nothing needed resizing - the cap was a bare
  literal with nothing behind it. What actually bounds N is that
  gravity is an explicit pair list, so memory and time both grow as
  `E*N^2`: at binary256 one system needs 1.1 GB at N=512 and 4.4 GB at
  N=1024, and one fixed binary64 step takes 51 s at N=512 and 196 s at
  N=1024. "Exercised to 512 with no trouble" was true and incomplete.
  Set at 1024. `ref/ias15_ref.c` had its own `struct body bodies[64]`
  and had to grow too, or the equivalence gate cannot reach the sizes
  the port now allows.
- A README section stating what is supported and what is refused,
  above the fold, so nobody discovers the limits by hitting them.
- A worked round trip: a REBOUND C program, run at binary128, result
  back in REBOUND, in the repository as a runnable example.
- `make install` or an equivalent that puts a header and a library
  where a REBOUND build can find them.

### D. Python, if it is cheap

Investigate first, implement only if the answer is short. REBOUND's
Python layer wraps the C library with ctypes, so a registered C
integrator may already be reachable by name from Python with little
more than a loader. If it needs its own binding layer, write down what
that would cost and stop.

---

## Testing policy for this round (as it was set)

**Quick tests only.** Every parcel runs the fast gates - the binary64
equivalence check, a short archive round trip, `make check-quick` - and
does not wait on the long suites. The full-format runs, the census
replays and the hardware sweeps belong to integration, once the pieces
are together, because running them per-parcel serialises four agents
behind a queue that proves nothing about their interfaces.

A parcel that cannot pass its own quick gate is not done, and saying so
is better than a long run that hides it.

---

## Not in this round, and why

- **Performance.** The measurements are in docs/HARDWARE.md and the
  integrator's log. In short: the tile beats one core from about 100
  ensemble members, reaching 13.3x at 16,384; four tiles are worth
  nothing over one because the workload runs the staged path and is
  bus-bound, where resident data would be worth 4x. None of that
  changes what the interfaces above should look like, and chasing it
  first would mean designing them around today's bottleneck.
- **Other integrators.** WHFast is ranked first in
  docs/INTEGRATORS.md and should follow the same shim, not precede it.
- **The device-side gather, scatter and scalar broadcast.** They are
  asks on cft-fp256, recorded in docs/HARDWARE.md and in that project's
  ROADMAP, and they are somebody else's parcel.

---

## Upstream defects found while building this

Recorded because they constrain the design, and two of them were found
twice, independently, which is why they are stated as facts rather than
suspicions.

- **`reb_integrator_register` faults on a SECOND custom integrator in
  one process.** The scan increments before testing and `strcmp`s the
  `{0}` terminator's NULL name. At `-O0` that segfaults; at `-O2` and
  `-O3`, which is how every published wheel is built, the compiler
  deletes the loop's exit via `strcmp`'s `nonnull` attribute and the
  call never returns. **One custom integrator per process**, and a
  stale second copy of our library is a hang rather than an error.
  Colliding with a built-in name is clean.
- **`element_size = W` does not survive a load.** The reader computes
  `n_elem` as `size_data / element_size` using the descriptor list the
  integrator was REGISTERED with, so a binary128 or binary256 archive
  read back through a binary64 registration gets the wrong count.
  Putting the count field last rescues a whole snapshot but not an
  appended one, because those are stored as a diff that omits unchanged
  fields. The count must be probed from the file.
- **`offset_N = SIZE_MAX` is broken by the address fixup** - the
  sentinel becomes `base - 1` and the read path writes eight bytes
  below the object. REBOUND's own `display_settings` uses it; nothing
  here does.
- **Integrator names must be lowercase.** The Python layer lowercases
  on assignment. `ias15_cft` is safe.
- **`reb_binarydata_diff` writes a field header with no data for a field
  that has disappeared.** The not-found branch (`binarydata.c:325-333`)
  emits `write_to_stream(&field1)` and `write_to_stream(name1)` and
  stops: no third call for the data. So a field present in the old
  snapshot and absent from the new one lands in the diff claiming a
  `size_data` it did not write, and a reader trusting that length walks
  into the middle of the next record's header. Reproduced minimally
  against the function, and end to end as an **0xC0000005** in
  `gate_real --fp64` when this port let its alias family vanish once a
  regrow levelled the mark. The fix here is that a family, once
  published, keeps being published for the life of the binding.

  **A consequence for this repository's own archives.** The same defect
  makes mixing builds across the round-2 commit unsafe in one
  direction: an archive written by a build that has the three new
  scalars and appended by one that does not loses them from the diff,
  emits `cft_iterations_max_exceeded` with `size_data=8` and no data,
  and the snapshot will not load — REBOUND itself reports "The binary
  file seems to be corrupted." Old-written and new-appended is fine.
  Do not mix builds across that commit within one archive.
- **A body-count change re-reads IAS15's coefficient levels at a
  different offset.** All seven share one flat allocation, and
  `dpcast(ias15->b, N3)` slices it at the *current* `3N` on every step
  (`integrator_ias15.c:217`): `p[m] = dp + m*N3`. Below
  `N_allocated` - a removal, or a regrow under the high-water mark -
  `reb_integrator_ias15_alloc()` reallocates nothing, so the buffer
  keeps its contents and every level above 0 is read short by
  `old_3N - new_3N`; level 1 begins in the tail of level 0, and so on
  up. Measured on a merge from three bodies to two: level 0 unchanged,
  and REBOUND's `b1[3]` is bit for bit what its `b1[0]` becomes.
  This project believed, and wrote down, that a removal simply left a
  stale polynomial behind - each survivor inheriting its neighbour's.
  It does not, and the difference matters here, because every level is
  its own allocation in this port and the aliasing had to be
  *performed* to reproduce it (`ias15_engine_alias_resize()`). A
  regrow reads back the tail the shrink stranded above the live region,
  which is why that is done through a persistent flat shadow rather
  than in place; `case_remove_then_add` in tools/check_dropin.c fails
  without it and passes with it.

## The cost, which was nowhere in this repository

Parcel C measured it and it belongs at the top of any expectation.
**The software library is about 2,300x slower than REBOUND's own IAS15
at the same binary64** - 40.2 steps a second against 92,400 on the
outer solar system. The tile's best measured ensemble speedup is 13.3x,
so on hardware this is still roughly two orders of magnitude slower
than stock REBOUND at binary64.

That is not a defect and it does not want fixing. It is the same
statement cft-fp256's own README makes: this is not for binary64
throughput. What it is for is the precision commodity hardware does not
offer and the guarantee that the answer does not move. Anyone choosing
this over REBOUND for a binary64 run has chosen wrong, and the README
now says so in those terms.

## A footgun found by building on it

`ias15_cft --max-iter` defaults to REBOUND's 12 **at every format**, so
a binary256 run silently truncates its corrector: 20 steps measured
`mean_pc = 12.000, max_pc = 12` - the cap on every single step. Given a
format-appropriate limit the same run converges at `mean_pc = 14.450,
max_pc = 18`, and the energy differs in the last twelve hex digits.

Nothing is wrong with the bits it produced; they are the correct answer
to a truncated iteration. But it is silent, and anything wrapping this
program must set the limit by format. Parcel C's packaging layer picks
12/24/60. **The shim in parcel A must do the same.**

## Corrections to the state struct above

From building against it: there is no `at` member though the mirror
list names one (safe - it is per-substep scratch); `g[7]` as seven
pointers cannot be a single field with `element_size = 7*8`, so it
archives as seven; `char cft_abi[16]` cannot be archived by any dtype
and goes out as two `uint64`; and nothing records that a state was
PROMOTED rather than restored exactly, which wants a one-word
`provenance` member.

Two members are **aspirational, not ports of existing behaviour**:
`adaptive_mode` (the port has no mode switch - PRS23 only) and
`min_dt` (no such option exists). Implement them or drop them
deliberately; do not assume there is behaviour to wrap.

**Settled:** neither was implemented, and both are refused rather than
ignored. `src/reb_integrator_cft.c` stops the integration with a named
REBOUND error for any `adaptive_mode` but PRS23 (2) and for a non-zero
`min_dt`. The members stay in the struct because they are archived and
removing them would change the on-disk field list. `provenance` was
**not** added and is still owed: nothing in the state records that it
was promoted from binary64 rather than restored exactly, and only the
return value of `cft_archive_finish_load()` carries that, which a
caller may discard.

## Two REBOUND defaults that will bite the refusal checks

Parcel C's first refusal list rejected every default simulation.
`r->OMEGAZ` is initialised to **-1.0** as a sentinel, not 0, and
`r->N_active` defaults to **SIZE_MAX**, not -1. Any check for "the user
set this" must test against those, and parcel A is writing the same
checks.

## Splitting a run changes the answer, measured

IAS15 starts each step's corrector from the previous step's `b`
coefficients. One call of 20 Kepler steps against two calls of 10, with
the state fed back exactly, differ by **1 to 6 ulps in four of twelve
coordinates**.

That is the empirical case for everything above: a restart that does
not carry `b` and `e` is not the same run, it is a nearby one. It is
why the state struct lists them, why the archive has to persist them,
and why parcel B's restart gate compares all the blobs rather than just
the particles. (It said "all 48 blobs" when written. There are 50; the
gate reads the count from `CFT_N_BLOBS` and no longer has a number of
its own - see the last section of this file.)

---

## What integration found, and what it costs to state a struct twice

The four parcels landed and every gate passed, and a checkpoint through
the registered integrator still did not work. Recorded here because the
cause is a property of how this document split the work, not of any one
parcel.

**The struct above was the shared contract, and it was wrong in one
place.** It lists `x0`, `v0`, `a0` as "position, velocity,
acceleration". They are not: `step_attempt()` copies `x` into `x0` at
the top of a step and leaves the advanced position in `x`, so after a
completed step `x0` is the *previous* step's start. The live coordinates
are `x` and `v`, and this document never mentioned them, so neither
parcel archived them. At binary64 that is invisible - `r->particles`
carry the same bits and the step re-promotes them - and at binary128 and
binary256 it silently truncated a checkpoint to binary64. The struct
now carries `x` and `v`, and the blob count is 50.

**Two parcels each implemented half of one thing, and the halves never
met.** Parcel A's brief said `src/cft_ias15_fields.c` is "PARCEL B's
file"; parcel B put its descriptor lists in `src/cft_archive.c` as
file-scope statics and never touched A's file. Both were faithful to
this document. The result was an integrator registered with a list
holding nothing but its terminator, so REBOUND could not resolve a
single `cft_` field on read - and because `cft_archive_finish_load()`
validates the FILE rather than what was restored, it returned "restored
exactly" for a restore that had not happened. A silent wrong answer,
which is the one outcome this repository is built to prevent.

**Where a shared fact appeared twice, the copies drifted.** The state
struct was defined in both A's and B's headers (caught at merge, found
byte-identical). The descriptor lists were in two files. The
index-to-member walker over the blobs was in two files. The blob count
48 was written out by hand in four places, and the set of blob names a
fifth. Growing the list from 48 to 50 broke every one of those copies
in a different way: a refused load, three crashed gates, and a passing
gate that printed "only 50 of 48".

**The rule this suggests for a next round.** A parcel brief may name a
shared fact, but exactly one parcel must own the file that states it,
and the others must be told to include it rather than restate it. Where
a count or a name list can be derived from a definition, derive it: blob
identity now comes from the descriptor list, and the count from
`CFT_N_BLOBS` in the header that defines the macros.

**The gate that would have caught all of it** is short and takes under
a second: N steps straight against a save and a resume that add to N,
compared bit for bit. It is `tests/gate_real.c`. (It was 30 = 20 + 10
when written; the defaults are 180 = 60 + 120 now, because at 20 + 10
only 2 of the 14 dumped lines differ between binary64 and binary256 and
a truncated restart could pass on the other 12. `--steps-a`/`--steps-b`
set them, and `tools/check_checkpoint.py` drives the same binary in
three separate processes.) Neither
parcel was asked for it, because it belongs to no parcel - it tests the
seam. A round split into parcels should name the seam tests too, and
give them to the integrator.

See docs/VALIDATION.md entry 25 for the measurements.

## What the collision work leaves open, across a checkpoint

**Closed by round 2's parcel P4 — kept for the reasoning, which still
explains why the archive is shaped as it is.** The high-water mark, the
stranded tail and the counters all travel now; `gate_real`'s second
sequence is remove → save → resume → regrow, and it fails without the
fix. Two things the section below did not know, both measured while
closing it: the stranded tail is **not** confined to the seven
coefficient levels — `csx` and `csv` strand one too, and the predictor
reads it before writing it — and a remove → checkpoint → regrow at
`state->accurate = 1` is still not exact, because that mode shifts the
wide state rather than aliasing it, so the tail sits in the live arrays
rather than in the shadow. What follows is the original statement of the
problem.


Reproducing REBOUND's coefficient re-slice (see the upstream-defects
section) needs the part of its buffer that a shrink strands above the
live region: no stride reaches it while the count is small, and a
regrow under the high-water mark reads it straight back. Within one
process `ias15_engine_alias_resize()` keeps it in a flat shadow, and
`case_remove_then_add` in tools/check_dropin.c holds that down - the
case fails if the shadow is dropped between resizes.

**Across a checkpoint it is not kept, and neither is the mark.**
REBOUND stores its coefficient arrays as one `REB_POINTER` field of
`N_allocated` elements of `7*sizeof(double)`, so the stranded tail and
`N_allocated` itself both survive a save. This port stores each level
as its own blob of `3N` elements and has no `N_allocated` at all: the
shim sets its high-water mark from `r->N` when it binds. So a run that
removes a particle, checkpoints, resumes, and then adds one back below
the *original* mark diverges twice over - the port zeroes where
REBOUND aliases, because its mark is now the smaller count, and it
would read zeros where REBOUND reads the tail even if the mark were
right.

Closing it means archiving both, which adds blobs and a scalar and so
changes the on-disk field list. Left open deliberately: the sequence
needs a removal, a save, a resume and a regrow that stays under a mark
the resumed process cannot see, and every other removal path - a merge,
a bare `reb_simulation_remove_particle`, a regrow within one process -
is bit-identical and gated.

## What round 2 leaves open: the wide clock

`t` is carried across a checkpoint as the binary64 view and as the
unevaluated pair the record prints, but **the engine's wide clock is
not archived**. Measured while verifying parcel P4: a plain checkpoint
at binary128 differs from the straight run by one binary64 ulp in `t`
for 2 of 40 step splits tried, and the same at binary256. It reproduces
on a build from before that parcel, so it is not something the parcel
introduced — but for work titled "state completeness across a
checkpoint" it is a piece of wide state that still does not survive one,
and it was recorded nowhere until now.

It is narrow: the particles, the coefficients and the compensations all
survive exactly, so a resumed run continues on the right trajectory and
only its clock differs, in the last bit, on some splits. Closing it
means archiving the wide `t` as a blob, which is one more field and one
more line in the descriptor list.


## What entry 37 leaves to build, scoped

**2026-09-14.** The measurement in docs/VALIDATION.md entry 37 ordered
the card-side work: cut the call count, then fix the divide, then
residency. Each item, sized from the tree as it stands rather than
from the design document, so the next round can pick one up without
re-deriving it.

### 1. Cut the call count (the biggest lever, rebound-side)

172,000 elementwise calls of six elements in a 200-step Kepler run, at
108 us each on the card: 77% of the step. The predictor and the seven
correctors are already programs (two `cft_program_run_ex` per substep
pass); what remains as hundreds of separate calls is **gravity** (about
twenty vectorised calls per force evaluation, `gravity_body()` in
src/ias15_cft.c) and the **end-of-step update and step control** (a
few dozen scalar-width calls a step, `pc_error()` and `dtnew_*()`).

- Gravity as one program needs its inputs gathered per pair - `x[i_l]`
  and `x[j_l]` from the predictor's deposits - and its outputs scattered
  per particle. **Nothing on the device does that today** (docs/
  HARDWARE.md, "The resident design", steps 2 and 3): the host does it
  with the byte-copy loops entry 37's split will price. So this item
  splits in two: the pair arithmetic as a program (possible now - a
  program over P lanes whose inputs are the six gathered coordinate
  vectors; removes ~15 of gravity's ~20 calls per evaluation, and the
  gather stays on the host as it is), and the gather/scatter on the
  device (a cft-fp256 sequencer ask: an indexed operand read, or an
  in-program cross-lane reduction, both already on cft-fp256's asks
  list from the workload round).
- The step control as a program: `predict_next_step`'s sixteen
  compensated additions per coordinate are named in the design as "one
  more program over the scratch block" and are still elementwise calls.
  Mechanical; the same shape as `correct_n`.

### 2. Adopt `cft_run_ex` with `scalar_mask` (delivered, unadopted)

cft-fp256 shipped the scalar-broadcast operand in ABI 0.12 because this
project asked for it (docs/HARDWARE.md's fourth ask), behind
`CAPS2[7]`; both f128 images and the revision-4 pair publish it, and the
software backend always carries it. `grep` finds **24 call sites** in
src/ias15_cft.c that pass a broadcast constant vector as an operand
(`vmul`/`vfma`/`vadd`/`vsub` with a `K*` vector: `KRINV`, `KHF`, `KC`,
`KY4`, `KY5`, `KNUM`, `KVNUM`, `KPOSR`, `K5040`, ...), each staging N3
copies of one value. Adopting it: a `vmulk(d, a, k)` family that calls
`cft_run_ex` with the constant as a one-element operand and
`scalar_mask` set, used at those 24 sites; the bits are the contract's
same bits (a scalar operand applies element 0 to the whole run), so the
gates hold it. Worth the ~300 staged vectors a step the design counted,
which entry 37 prices at a few percent - bytes are not the wall. Do it
for the memory and the honesty of having asked, not for speed.

### 3. The recorded corrector schedule (prototype exists)

`~/patch_fixed_schedule.py` on amd-arc-box (2026-09-10, against a tree
123 commits behind) adds `--record-iters FILE` and `--replay-iters
FILE`: the first writes one corrector pass count per step attempt, the
second runs exactly those counts and skips `pc_error()` at both of its
call sites (loop and program engines) and the corrector's own exit
test. The replay must produce a byte-identical record, and that
identity is the proof the optimisation is exact - it is a claim only a
bit-reproducible integrator can make. Entry 37's measurement says what
it is worth: `pc_error` is inside the elementwise bucket at small N and
is 13-19% of calls; on the card that is 4-9% of a step. Port the
prototype's four hunks to the current tree (its anchors are the
`if (n == 7) pc_error(...)` lines, the corrector loop and the argument
parser; all four still exist), gate it on the identity, and record it.

### 4. The measurement before the next round (running as this is written)

Entry 37's split cannot say how much of gravity's 86% is the host's
byte-copies against the accumulate adds against the divide against the
launches. `src/ias15_cft.c` now times the three pair gathers and the
scatter loop as `t_g_copy`, the accumulate half as `t_g_accum`, and
every library call inside gravity as `t_g_lib` with `t_g_divsqrt`
separately; the decomposition pass re-runs with them and entry 38 will
carry the table. Whether item 1 starts with the pair program or with
the device gather depends on it.

### 5. Residency, demoted to last, still real

`cft_alloc` for `x0`, `v0`, `a0`, `at` and the four program buffers,
then the scratch block resident between corrector passes. docs/
HARDWARE.md's "Moving the state onto the card" has the byte accounting
and the library side is verified done; entry 37 prices the whole item
at 1-10% of a step. Worth doing after 1 and 3, not before.

### Still owed from before this round, unchanged

- The `provenance` member (this file, "Corrections to the state struct
  above"; PARCELS.md's parcel-C note).
- The wide clock `t` across a checkpoint: 1 ulp on 2 of 40 splits
  ("What round 2 leaves open: the wide clock").
- The full software suite on Linux, entry 29's open item: `make check`
  on amd-arc-box, 2026-09-14, result in entry 38.
- docs/HARDWARE.md "What to do first" items 1 and 4: the ensemble gate
  with an artifact, and the mixed-pass-count ensemble that prices idle
  lanes - both in the chain running as this is written, results in
  entry 38.
