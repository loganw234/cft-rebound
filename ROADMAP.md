# Roadmap: from a correct experiment to a usable integrator

Written 2026-09-10, after the port ran on the FP256 tile.

The hard part is done and it is the part that is normally hard:
**correctness**. At binary64 this is REBOUND's own IAS15 bit for bit;
on the card every record is bit-identical to the software backend,
including an ensemble reproducing its members run one at a time. What
is missing is not accuracy. It is that nobody can use it.

Today `ias15_cft` is a standalone program that reads a bespoke problem
file and writes a bespoke record. A REBOUND user has a
`struct reb_simulation` in C, or a `rebound.Simulation` in Python, and
there is no path from one to the other. That is the gap.

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

Each is independently reviewable and they share only the struct above.

### A. The REBOUND integrator shim (the drop-in)

Implement `struct reb_integrator` and register it.

- `create` / `free`: allocate the state for `r->N` particles.
- `step`: promote `r->particles` to the wide format (exact - binary64
  is a subset of every wider one), take one IAS15 step through the
  existing engine, round the result back into `r->particles`, and
  advance `r->t` and `r->dt` as REBOUND expects.
- `did_add_particle` / `will_remove_particle`: resize and invalidate.
- Honour `r->exact_finish_time`, and `synchronize` where it applies.

**Gate:** a C program that builds a simulation, registers this
integrator at binary64, integrates, and gets bit-identical particle
values to the same program using REBOUND's own `"ias15"`. That is the
existing equivalence gate seen from REBOUND's side, and it is the one
result that proves the shim did not change the arithmetic.

Deliberately out of scope: additional forces, collisions, ghost boxes,
the tree code, non-zero softening, variational particles. Detect them
and refuse with a clear message rather than computing something wrong.

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

## Testing policy for this round

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
and why parcel B's restart gate compares all 48 blobs rather than just
the particles.

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

**The gate that would have caught all of it** is 190 lines and takes
under a second: 30 steps straight against 20 steps, save, free, load,
10 more, compared bit for bit. It is `tests/gate_real.c`. Neither
parcel was asked for it, because it belongs to no parcel - it tests the
seam. A round split into parcels should name the seam tests too, and
give them to the integrator.

See docs/VALIDATION.md entry 25 for the measurements.
