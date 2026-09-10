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

It **warns and seeks past it**. Which means:

- a `cft_`-prefixed field cannot collide with any present or future
  REBOUND field, because names are unique strings rather than indices;
- an archive written by cft-rebound **opens in stock REBOUND**, which
  recovers the full binary64 state and warns about the extra fields;
- an archive written by stock REBOUND opens here and is promoted.

Direct compatibility and native wide-precision support are therefore
the same file, not a choice between them. That removes the barrier to
entry the long way round: a collaborator without cft can still read
your archives.

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
    unsigned char *x0, *v0, *a0;      /* position, velocity, acceleration */
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

- Raise the 64-body cap. It is arbitrary: `valloc` is `calloc`, every
  array is sized from N at run time, and the only fixed object is
  `body_names[64][32]`, which is display text. It has been exercised to
  512 bodies on a scratch clone with no trouble.
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
