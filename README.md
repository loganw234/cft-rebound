# cft-rebound

REBOUND's IAS15 integrator with its arithmetic routed through libcft,
the IEEE 754-2019 binary32/64/128/256 library of the cft-fp256
project, so that the same integrator runs at binary64, binary128 and
binary256 - and, given an artifact, on the FP256 tile - with the same
bits everywhere.

It exists to answer one question with measurements: IAS15 was
designed by Rein & Spiegel to sit at the binary64 round-off floor.
Does a wider format push that floor down and extend the usable
integration horizon, and at what cost? docs/VALIDATION.md has the
numbers; the short answer is in "What was found" below.

## Scope: what it integrates, and what it refuses

Read this first; it decides in about a minute whether this is any use
to you.

**What it does.** REBOUND's IAS15, unchanged as an algorithm, with
every floating-point operation carried out at binary64, binary128 or
binary256, up to **1024 bodies** per system. At binary64 it is
REBOUND's own IAS15 bit for bit, which is a gate the port passes on
every build. An ensemble of independent systems integrates in one run,
each with its own step (docs/ENSEMBLE.md).

Gravity is REBOUND's `REB_GRAVITY_BASIC` - every pair, nothing
approximated - and on the **drop-in integrator** it now comes with most
of what a real script puts around it:

| | |
|---|---|
| `r->additional_forces` | called after gravity at **every Gauss-Radau node**, exactly where REBOUND's IAS15 calls it, with `r->t` set to REBOUND's own `t_beginning + dt*h[n]`. This is the REBOUNDx path |
| velocity-dependent forces | the node's velocities are predicted first, under REBOUND's own condition for doing so |
| `pre_`/`post_timestep_modifications` | REBOUND's driver calls these between steps and the port carries them |
| collision detection and resolution | including REBOUND's own resolvers. Bit-identical across a removal |
| test particles (`N_active`), `testparticle_type` | the pair set follows REBOUND's two loops, in REBOUND's order |
| all four step criteria | INDIVIDUAL (0), GLOBAL (1), PRS23 (2, REBOUND's default), AARSETH85 (3) |
| non-zero `softening`, IAS15's `min_dt` | as REBOUND computes them |
| `r->gravity_ignore_terms` | accepted because IAS15 overwrites it with `NONE` at the top of every step, so nothing a caller sets ever reaches gravity - refusing it refused a run that would have been right |

Everything below is refused, **by name, with a message** - never
silently ignored, never approximated. The list is one table in code
(`cft_support_rows`, src/cft_supported.c) that both entry points walk,
and the gate walks it too and fails on any row no case exercises.

**Refused on both paths:**

| | |
|---|---|
| variational particles | the wide state carries the real particles only. MEGNO needs them, so the ordinary route into it is refused here too |
| ghost boxes, periodic or shear boundaries | not ported |
| the tree code (`REB_GRAVITY_TREE`) | this is direct summation. `REB_GRAVITY_COMPENSATED` is refused too: a different summation from the one ported |
| a custom gravity routine (`r->gravity_custom`) | the engine issues gravity itself, so it would never be called |
| particle maps (`r->map`) | a setting of the *integrator*, not of gravity: IAS15 runs its state over `N_map` elements while gravity still runs over all `N` |
| `N_active > N` | a heap overread in REBOUND's own gravity loop. There is no bit-identity to claim against undefined behaviour |
| an `adaptive_mode` naming none of the four | |
| every integrator except IAS15 | WHFast is ranked first to follow (docs/INTEGRATORS.md) and has not been done |

**Refused on the subprocess API only**, because it hands the run to
another process that cannot call a function pointer in this one, and
because REBOUND's driver does not run between its steps:
`r->additional_forces`, velocity-dependent forces, the timestep
modification hooks, and attached ODE sets. All four work on the
drop-in.

**Refused differently by each path.** Two rows apply to both but ask a
different question on each, which is why they are neither in the table
above nor in the list before this one:

| | drop-in | subprocess |
|---|---|---|
| collision detection | refused only above binary64 without `state->accurate = 1` | refused outright — REBOUND's driver runs the search between steps there and not here |
| test particles (`N_active`) | refused only for `N_active > N`, which is undefined behaviour upstream | refused for any inactive particle |

**Refused on the drop-in only**, because they are its own settings:
MEGNO set directly (`r->calculate_megno` without variational
particles), a `format` or `max_iter` out of range, a state whose `E` is not
1 (an ensemble is E independent systems and a `reb_simulation` is one -
use the standalone program, docs/ENSEMBLE.md), and `arith_fma` combined
with a velocity-dependent force, which would need an FMA form of the
velocity predictor that has no oracle to check it against.

**And one refusal that is not about a setting:** a force routine may
write `ax`, `ay` and `az` and nothing else. REBOUND carries more out of
that call - a position or velocity written at the top of a step attempt
becomes the step's initial condition, a velocity written at a node
survives into the next one while `force_is_velocity_dependent` is
clear, and a mass is picked up by gravity at the next node - so a
routine that writes them would answer a different problem here. The
port refuses and **discards the step in progress**, leaving the
simulation exactly as it was. Edit coordinates from
`r->pre_timestep_modifications` instead, which REBOUND's own driver
calls between steps and which this port does carry.

Both entry points walk the same table. `cft_rebound_check()`
(src/cft_rebound_run.c) is the subprocess API's half, and
`cft_rebound_steps()` calls it before doing anything, so a refusal
arrives at the start of your run rather than in the middle of it;
`supported()` (src/reb_integrator_cft.c) is the drop-in's, checked at
the top of every step. The rows above say which path each applies to,
and that column is the whole of the difference between them - it used
to be the emergent result of two functions written months apart.

Two details worth knowing before you pick a path.

**Collisions above binary64 need `state->accurate = 1`.** What the
integrator owes across a removal is to survive it the way REBOUND's own
IAS15 does - which is not by keeping a stale polynomial, as this project
believed until it measured one. REBOUND's seven coefficient levels share
a single buffer that `dpcast()` re-slices at the *current* `3N` on every
step, so a removal below its high-water mark leaves the buffer alone and
re-reads every level at an offset short by `old_3N - new_3N`. The port
performs that reading rather than approximating it, so a merge at
binary64 is bit-identical across the removal (`case_collision` in
tools/check_dropin.c). Above binary64 the cost is elsewhere: a removal
shifts `r->particles`, which are binary64, and the default re-promotes
the survivors from them and loses every wide tail. `accurate = 1` shifts
the wide state instead, which is more accurate than REBOUND and so
deliberately not bit-identical to it.

**A force routine is evaluated at binary64 even in a binary256 run.**
`r->particles` are binary64 and they are the only interface REBOUND
offers a callback, so there is nowhere wider to hand a routine its
inputs or take its answer back. Gravity stays wide; a component the
routine writes is binary64 from that node on, and one it leaves alone
keeps its wide value exactly. That is a real ceiling of the feature, not
a detail.

**What it costs, and this is the part to weigh.** The arithmetic is a
software library, not the hardware's doubles. On the outer solar
system (N = 6, fixed step) REBOUND's own IAS15 does 92,000 steps a
second on this host and the port does **40** at binary64, 11 at
binary128 and 3.6 at binary256. The first factor of 2,300 buys you
nothing but the ability to change format; only the rest buys accuracy.

Gravity is an explicit pair list, so time and memory then grow as
N^2 - per system, and an ensemble of E systems multiplies both by E:

| one system | N = 64 | N = 256 | N = 512 | N = 1024 |
|---|---|---|---|---|
| memory, binary64 | 8 MB | 79 MB | 283 MB | 1.1 GB |
| memory, binary256 | 23 MB | 283 MB | 1.1 GB | 4.4 GB |
| one fixed binary64 step | 1.1 s | 14 s | 51 s | 196 s |

The 1024-body cap is set where the worst case, binary256, still fits
an ordinary workstation: 2048 would want 17.7 GB. Time stops you well
before memory does, so a few hundred bodies is a demonstration rather
than a survey. What this is for is a *small* system integrated
*accurately*; docs/HORIZON.md measures where that is worth doing, and
is equally clear about the many cases where binary64 is entirely
sufficient.

**What you get back.** REBOUND's `struct reb_particle` holds `double`
and will keep doing so, so `r->particles` is the binary64 **view** of
the run, correctly rounded from the wide state, and every REBOUND
output, callback and tool works on it unchanged. The wide state is the
integrator's own. `struct cft_rebound_result` reports the run's energy
both as a rounded `double` and as exact hexadecimal text, because at
binary128 and above the digits that moved are past binary64's last
bit.

## Licence and provenance

This repository is **GPL-3.0-or-later** (LICENSE, NOTICE). REBOUND is
GPL-3.0-or-later; libcft is Apache-2.0 and the cft-fp256 repository
states a permissive-only dependency rule. Apache-2.0 combines into a
GPL-3.0 work and not the reverse, so the combined work lives here,
under REBOUND's licence, and no REBOUND source is ever copied into
cft-fp256.

Both upstreams are **pinned clones** under `third_party/`, recorded by
commit in `third_party/MANIFEST` and reproduced by
`tools/fetch_third_party.sh`, which fails hard on a wrong commit.
Nothing is vendored by hand.

| upstream | commit | licence |
|---|---|---|
| REBOUND (hannorein/rebound) | bdfda4bd (5.1.1 + 6) | GPL-3.0-or-later |
| cft-fp256 (loganw234/cft-fp256) | 40ce35e3 (libcft ABI 0.11) | Apache-2.0 |

## Layout

    include/cft_rebound.h    the subprocess API: what a REBOUND program includes
    src/cft_rebound_run.c    and what it links: the refusal list and the run call
    src/cft_ias15.h          the DROP-IN's public header: the state struct, the
                             registration, cft_ias15_configure()
    src/reb_integrator_cft.c the struct reb_integrator shim over the engine, and
                             the DROP-IN's refusal list (supported())
    examples/roundtrip.c     the worked round trip, runnable (`make example`)
    examples/dropin.c        the same from the drop-in side: register,
                             integrate at binary256, checkpoint, reload
    examples/Makefile        a REBOUND program's makefile with the added lines,
                             marked, for each of the two ways in
    src/ias15_cft.c          the port: every floating-point operation is a cft.h call
    src/ias15_engine.h       the same file as a library: the entry points the shim
                             drives (-DIAS15_CFT_LIBRARY, no main)
    src/ias15_constants.h    GENERATED: the Gauss-Radau constants at every format
    src/hexfloat.h           exact hex-float text for binary64, libc-independent
    src/cft_ias15_state.h    the width helper, and the include that reaches the
                             state struct in cft_ias15.h
    src/cft_archive.h/.c     Simulationarchive: the cft_ fields, the probe, the load
    src/cft_ias15_fields.h/.c  the ONE definition of the field descriptors, and
                             CFT_N_BLOBS - a literal beside the list it
                             counts, checked against it by the selftest.
                             CFT_N_ALIAS_BLOBS and CFT_N_SCALARS are two
                             more on the same terms
    tests/cft_shim_stub.c    a stand-in integrator, for four of the archive gates
    tests/gate_restart.c     gate 1: checkpoint and restart, bit for bit
    tests/gate_write.c       gate 2a: write a binary128 archive, and record the
                             exact binary64 view a stock reader must recover
    tests/gate_stock.c       gate 2b: that archive opened by a program linked
                             against upstream REBOUND and nothing of ours
    tests/gate_promote.c     gate 3: a stock archive promoted, and the refusal
    tests/gate_real.c        the seam: the REAL integrator through an archive
    ref/ias15_ref.c          REBOUND's own IAS15 on the same problems, plain double
    ref/whfast512_stub.c     why REBOUND's WHFast512 is not built here (MinGW)
    tools/gen_constants.py   derives and CHECKS the constants (mpmath, 130 digits)
    tools/gen_programs.py    the predictor and corrector as sequencer programs
    tools/make_problems.py   the problems as exact binary64 bit patterns
    tools/make_ensemble.py   E perturbed copies of a problem as one ensemble file
    tools/make_nbody.py      a deterministic N-body problem of any size
    tools/check_equivalence.py     the gate: binary64 port == REBOUND, bit for bit
    tools/check_program_engine.py  the gate: programs == host loop, bit for bit
    tools/check_ensemble.py        the gate: ensemble == its members run alone, bit for bit
    tools/check_records.py         the gate: the committed records, recomputed, bit for bit
    tools/check_archive.py         the gates: an archive round trip, and both readers
    tools/check_bodycount.py       the gate: the port == REBOUND at a few hundred bodies
    tools/check_checkpoint.py      the gate: a checkpoint ACROSS processes, with a
                                   negative control that must differ
    tools/check_dropin.c           the gate: the registered integrator == REBOUND's
                                   own ias15, and every refusal refused
    tools/oracle.py          scores a record from its exact bits (mpmath)
    tools/compare_formats.py the round-off floor: one run at two formats, differenced
    tools/horizon.py         percentiles of an ensemble's error against time, and crossings
    tools/divergence.py      members of an ensemble against member 0: the Lyapunov measurement
    tools/reversal.py        forward N steps, back N steps, from the exact recorded state
    tools/state_to_problem.py      a record's sample as a problem file, bit for bit
    tools/sweep.py           the precision-versus-error measurement
    tools/tabulate.py        its tables and figure
    programs/*.cfta          GENERATED: 27 orbit-sequencer programs
    data/problems/           kepler.txt (e = 1/2), kepler_e09.txt, kepler_e099.txt,
                             pythagorean.txt (Burrau's problem), outer.txt
    data/ias15_constants.json
    results/                 the sweep's per-sample tables and summary
    docs/VALIDATION.md       what ran and what it said, failures included
    docs/ENSEMBLE.md         E systems in one run: the layout, the step decision, the gate
    docs/HORIZON.md          when binary64 stops being enough, measured
    docs/HARDWARE.md         the resident design, its break-even, what is unverified
    docs/INTEGRATORS.md      which other REBOUND integrators port well, ranked
    docs/PYTHON.md           reaching the integrator from rebound.Simulation, and
                             the four sharp edges that come with it
    python/cft_rebound.py    the loader: puts a registered integrator within reach
                             of REBOUND's Python package
    python/example_equivalence.py  the binary64 equivalence gate, from Python
    src/cft_ias15_shared.c   the constructor that registers the integrator
                             when the shared library is loaded
    tools/fetch_third_party.sh     clones the pinned upstreams, hard-fails on a
                                   wrong commit
    ROADMAP.md               the plan that got this here, all four parcels landed,
                             plus what integration found and the upstream defects

## Building

Plain C99 and GNU make. On Linux or macOS, with git, gcc and a Python
that has mpmath:

    make third-party        # clone and verify the pinned upstreams
    make libcft             # libcft from the pinned clone
    make                    # librebound (static), ias15_ref, ias15_cft,
                            # check_dropin, libcft_rebound.a, libcft_ias15.a
    make programs           # assemble the sequencer programs (needs cft-asm)
    make check              # the gates (make check-quick for binary64 only)

On Windows the toolchain is MSYS2's mingw64 gcc from Git Bash, and the
same three traps cft-fp256's host/Makefile documents apply; pass them
as make variables:

    PATH="/c/msys64/mingw64/bin:$PATH" make CC=gcc OS=Windows_NT \
        TMP=C:/Users/you/AppData/Local/Temp TEMP=C:/Users/you/AppData/Local/Temp \
        PYTHON=C:/path/to/python.exe

If libcft was built with the XRT backend - the one that talks to a
card - its objects are C++ and need the XRT runtime at link time.
Sourcing XRT's own setup script is enough:

    . /opt/xilinx/xrt/setup.sh    # sets XILINX_XRT
    make                          # the link flags follow from it

`XRT=0` turns that off for a software-only libcft in a shell that
happens to have XRT sourced, and `LIBS=` overrides the lot. Without
it a plain `make` against such a library stops on a page of
undefined references to `xrt::bo`, which names neither the cause nor
the fix.

REBOUND's own build assumes MSVC on Windows; here its library is
compiled with gcc from its sources minus `integrator_whfast512.c`,
whose private `__m512d` typedef collides with the one mingw's
`windows.h` brings in, plus a stub (ref/whfast512_stub.c). IAS15 does
not touch WHFast512. `-std=c99` pins `-ffp-contract=off`, so no FMA
contraction changes REBOUND's doubles.

## Using it from your own REBOUND program

    make install PREFIX=/usr/local        # DESTDIR=... also honoured

installs both ways in. Add the include path, then pick a library. If
your libcft was built for a card, a third line is needed - see the
block marked in `examples/Makefile`, or source XRT's `setup.sh` and let
that Makefile derive it.

### The drop-in

    CFLAGS  += -I/usr/local/include
    LDLIBS  += -L/usr/local/lib -lcft_ias15 -lcft

The integrator runs inside your process and REBOUND's own driver
calls it, so your program keeps the API it has. Two lines are new:
register once, and say which format.

    #include "rebound.h"
    #include "cft_ias15.h"

    cft_ias15_register("ias15_cft");            /* once per process */

    struct reb_simulation *r = reb_simulation_create();
    /* ... your particles, G and dt, exactly as now ... */
    struct cft_ias15_state *s =
        reb_simulation_set_integrator(r, "ias15_cft");
    s->format  = CFT_FP256;      /* CFT_FP64 is REBOUND, bit for bit */
    s->epsilon = 1e-9;           /* 0 is REBOUND's fixed step */
    s->max_iter = 60;            /* see below - the struct default is 12 */

    reb_simulation_steps(r, 1000);               /* REBOUND's own call */

**Set `max_iter` above binary64.** The state a fresh
`reb_simulation_set_integrator()` hands you has `max_iter = 12`,
REBOUND's own number, at *every* format, and binary256 needs about 18
to 22 passes to converge its corrector - so a binary256 run that leaves
it alone hits the cap on every step and truncates the iteration
silently. The bits it produces are the correct answer to a truncated
iteration, not a wrong one, but they are not the answer you asked for.
`cft_ias15_configure(r, format, epsilon, 0)` picks 12/24/60 by format;
setting the member yourself does not.

Add `cft_archive.h` and a `cft_archive_bind(r)` before
`reb_simulation_save_to_file()` and the Simulationarchive carries the
wide state, so a run can stop and be picked up later at full width.
`examples/dropin.c` is all of that in one runnable file, built against
the install rather than against the source tree.

### From Python

    make python-lib          # build/libcft_ias15.so (.dylib on macOS)
    make check-python        # REBOUND's own ias15 against it, from Python

Then, with the `rebound` package installed:

    import sys, rebound
    sys.path.insert(0, "python")
    import cft_rebound

    lib = cft_rebound.load("build/libcft_ias15.so")

    sim = rebound.Simulation()
    sim.dt = 0.05
    sim.add(m=1.0)
    sim.add(m=1e-3, x=1.0, vy=1.0)
    sim.integrator = "ias15_cft"
    cft_rebound.configure(lib, sim, format="fp256", epsilon=0.0)
    sim.steps(1000)

The library registers itself when it loads, so the name is selectable
immediately. The `configure` call is not optional if you want anything
but the default: REBOUND builds `sim.integrator.epsilon` and its
siblings from structs it knows and does not know this integrator's, so
that is where the format is chosen. At binary64 the run is REBOUND's
own, bit for bit, which is what `make check-python` checks.

The pinned REBOUND and the installed wheel must be the same source, or
the two sides of every call disagree about how `struct reb_simulation`
is laid out. `make python-lib` compares them and refuses if they
differ.

**Windows works**, by the opposite loading rule. A DLL may not carry
undefined symbols, so where the POSIX shared object leaves REBOUND's
symbols to the loader, the Windows build links the wheel's own
`librebound` - asked of the interpreter that will load it. That ties
the library to a Python minor version, and a mismatch is worse than it
sounds: with a *different* wheel's librebound findable beside it the
load **succeeds**, registers into a second REBOUND, and surfaces only
as `RuntimeError: Integrator not found`. So `load()` reads the import
table and refuses by name first. docs/PYTHON.md has the detail, and
why there is no `pyproject.toml`.

### The subprocess API

    CFLAGS  += -I/usr/local/include
    LDLIBS  += -L/usr/local/lib -lcft_rebound -lcft

For a caller that does not want the engine in its own address space.
`cft_rebound_steps()` writes a problem file and runs the `ias15_cft`
program, which the install also puts on the system. Include
`cft_rebound.h` instead of `rebound.h`, build your simulation exactly
as you do now, and call it where you would have called
`reb_simulation_steps()`:

    struct cft_rebound_options opt = {0};
    opt.format  = CFT_REBOUND_FP128;
    opt.epsilon = 1e-9;                   /* 0 is REBOUND's fixed step */
    struct cft_rebound_result res;
    cft_rebound_steps(r, 1000, &opt, &res);
    /* r->particles are now the binary64 view of a binary128 run */

`examples/roundtrip.c` is that in full - the same system at binary64
and binary128 on identical steps, differenced - and
`examples/Makefile` is an ordinary REBOUND program's makefile with the
two lines above marked. `make example` stage-installs into
`build/stage`, builds it against that, and runs it, so the example
also tests the install.

`cft_rebound_steps()` runs the `ias15_cft` program in a subprocess and
the wide state does not survive the call, so **ask for a whole run in
one call**: IAS15 starts each step's corrector from the previous step's
b coefficients, and a call that begins with them zeroed lands on
different bits - one to six ulps over twenty Kepler steps, measured in
docs/VALIDATION.md. That restriction is the reason to prefer the
drop-in above, where the state persists across steps and across a
checkpoint. This form was deliberately kept rather than rewritten as a
wrapper around it, for a caller who does not want the engine in their
own address space; `make example` builds and runs both, and each case
declares the answer it expects, so a change that breaks either is a
build failure.

That sentence was false until 2026-09-11. `roundtrip.c` used to *print*
what the library answered rather than assert what it should answer, and
`main()` returned 0 regardless — so when `softening` became supported,
the example went on shipping the line `-> NOT REFUSED, which is a bug`
about correct behaviour, and the build passed through it. A case that
prints rather than asserts will print whatever it is given, forever, and
look fine.

## Running

    build/ias15_ref --problem data/problems/kepler.txt --dt 0.05 --epsilon 0 --steps 2520 --sample 252
    build/ias15_cft --format fp256 --problem data/problems/kepler.txt --dt 0.05 --epsilon 0 \
                    --steps 2520 --sample 252 --max-iter 60
    python tools/oracle.py record.txt

Both programs write the same record: one line per sample with the
time, the next and last step, the energy and every coordinate as an
exact hex float, so that records are diffable bit for bit and the
oracle scores them from their exact bits. `epsilon 0` is REBOUND's
convention for a fixed step.

`ias15_cft` options that matter:

- `--format fp64|fp128|fp256`; `--artifact PATH` opens a tile instead
  of the software backend, and `$CFT_REBOUND_ARTIFACT` is the same
  thing for a program that takes no flag of its own.
- `--arith rebound|fma`: REBOUND's exact sequence of roundings
  (divisions and all), or the FMA form with reciprocal constants that a
  sequencer program can run. `--engine program` runs the predictor and
  corrector as the assembled programs (`--programs DIR`).
- `--cs kahan|augmented`: REBOUND's compensated summation, or the exact
  error term of 754-2019 9.5's augmentedAddition.
- `--max-iter N` (12 in REBOUND; the corrector needs about 18 passes at
  binary256), `--pc-tol-shift S` (the corrector tolerance 1e-16 is
  scaled by 2^-S, default 53-p), `--trace-pc N`.
- An ensemble file (`tools/make_ensemble.py`; `E` and `system` lines)
  integrates all its systems in one run, each with its own adaptive
  step (docs/ENSEMBLE.md); `--member K` runs system K of it alone,
  which is what the ensemble gate diffs against.
- `--dt-file FILE` replays a prescribed step sequence (one hex float
  per step, shared by every system; needs `--epsilon 0`) and
  `--dt-out FILE` records one; `--dt` and `--epsilon` take a decimal
  or an exact hex float.

### The two suite targets, and why one of them skips

`make check` runs every gate at every format. On this host that is about
an hour, and **it never skips anything** - it is the target whose
meaning must not change, so it always executes.

`make check-quick` runs the same gates with the cheaper arguments and
consults a cache: a leg whose inputs are byte for byte what they were
when that leg last passed is skipped. Measured on this host: **eighteen
minutes with an empty cache, five seconds when nothing has changed.**

The key is over the **artifacts a leg executes** - the binaries, its own
script, the data files it reads, and `$CFT_REBOUND_ARTIFACT` - and never
over source files. A source list has to enumerate every header and every
compiler flag, and the cost of missing one is skipping a leg that would
have failed, which is the single outcome this repository exists to
prevent. A binary is the closed-over result of all of it: if
`check_dropin` is byte-identical, nothing that could change its answer
changed. It works here because every verdict in this suite is
bit-identity - no tolerance, no timing, no seed - so a rerun of the same
bytes cannot answer differently.

Four things keep it from becoming a gate that cannot fail:

- **`make check` never reads it.** "The full suite passed" means what it
  has always meant.
- **Every skip prints itself**, by name and key, with a tally at the
  end. `make gate-cache-report` shows what is held.
- **A failure retires the stamp**, so a leg that failed cannot be
  skipped tomorrow on yesterday's pass.
- **The cache lives under the build tree** and is written by atomic
  rename, so it cannot travel between `build/` and `build2/` or between
  machines, and a half-written stamp cannot be read as valid.

`$CFT_REBOUND_ARTIFACT` is part of every key, so a pass on the software
backend can never satisfy a run against a card.

`make check-light` is the third tier, for when you have changed
something and want a fast signal from *everything* rather than a
thorough one from part of it. It runs every case with the step counts
scaled down. Measured on `check_dropin`: **159 seconds against 585, the
same 130 cases, zero failures.**

One case opts out by name, and the reason is the interesting part. The
merge case fires its collision at step 1291 of 2000; a tenth of the
steps is a run in which nothing happens. It does not pass vacuously - it
reports "no collision occurred ... so this case proves nothing" and
fails - because that guard has been there since the case was written.
Scaling it would have turned a real case into a no-op, silently, if the
guard had not been. It steps through `(reb_simulation_steps)(...)`,
where the parentheses suppress the scaling macro.

`check-light` is never a substitute for `make check`. It says whether a
change is obviously wrong, not whether it is right.

One honest consequence: a documentation-only change invalidates nothing,
so `check-quick` will skip every leg and pass having executed none. That
is correct, and it is why the skip lines are loud.

### Running the gates on a card

    export CFT_REBOUND_ARTIFACT=/path/to/cft_hw_quad.xclbin
    make check

One variable, and every gate that opens the engine opens the tile
instead - no gate takes a flag for it. The same variable is the
fallback for all three ways into this library, so they cannot
disagree about which device a run used:

- the registered integrator, after `cft_ias15_set_artifact()`;
- `cft_rebound_steps()`, after `opt.artifact`;
- the `ias15_cft` program, after `--artifact`.

An explicit setting always wins; unset or empty means the software
backend. Whichever is used, the bits are the same - that is what
the gates check - so the card is a speed choice, not a numerical
one.

Note that a gate is not a benchmark. The python checkers spawn a
fresh `ias15_cft` per case and would pay a device open every time;
`build/check_dropin` and `build/gate_real` open the engine once for
the whole program and are the ones worth pointing at a card.

Budget for it. On a U50C quad tile, `check_dropin` takes 2,030 s at
binary64 and 2,195 s for its binary128 pass, against a few seconds
each in software. That is the per-call cost at two and five bodies,
not a surprise - docs/HARDWARE.md predicts it and the ensemble numbers
are where the card earns its keep. What the two runs establish is that
the answer does not depend on the backend.

`make check` runs every gate at every format (about half an hour);

`make check-quick` runs them at binary64 in a few minutes. Both
include `tests/gate_real.c`, which takes a checkpoint through the
registered integrator and requires the restarted run to be bit
identical to an uninterrupted one, and `tools/check_bodycount.py`,
which re-runs the equivalence gate
at a few hundred bodies - where the per-particle summation order could
drift from REBOUND's and nowhere else - and checks that the 1024-body
cap is reachable and that 1025 is refused.

## What was found

See docs/VALIDATION.md for every number. In brief:

- **The constants.** REBOUND's `h[]`, `rr[]`, `c[]`, `d[]` are 25-31
  digit doubles; binary256 needs 73. tools/gen_constants.py derives the
  Gauss-Radau nodes at 130 digits (zeros of (P7+P8)/(1+x), confirmed as
  the zeros of the Jacobi polynomial P_7^(0,1) to 6e-130) and the
  arrays by REBOUND's own recurrence, and all 78 values rounded to
  binary64 reproduce REBOUND's literals bit for bit.
- **The port.** At binary64, `ias15_cft` is bit-for-bit REBOUND's IAS15
  on both problems, with fixed and adaptive steps: tools/
  check_equivalence.py, 1,264 recorded values identical.
- **Precision.** Measured on 41 runs (results/, tools/sweep.py):
  **binary128 is worth having and binary256 mostly is not.** At
  REBOUND's default epsilon the binary64 Kepler run is round-off
  limited by four orders of magnitude (5.3e-15 against the 4.6e-19 the
  method delivers), so binary128 buys those four orders at once, ten at
  epsilon 1e-12, and eighteen at 1e-16 (8.6e-34 against 4.7e-16 over
  195 orbits, its own floor). binary256 is identical to binary128 at
  every step size a user would choose, because a fifteenth-order
  method's truncation error - which scales as dt^15, measured - sits
  far above the binary128 floor until 800 steps an orbit (Kepler) or
  10-20 day steps (outer solar system, where it is then worth one to
  six orders more). The cost is the corrector: 2-3 passes a step at
  binary64, 4-9 at binary128, 10-20 at binary256. The wide-format
  error grows linearly in time - a ruler, not a random walk.
  Replacing Kahan's `add_cs` with 9.5's exact augmentedAddition lowers
  the binary64 floor 3-4x on the one run tried and changes nothing at
  binary256.
- **The programs.** The predictor and corrector run as orbit-sequencer
  programs with the state in the per-lane scratch block, bit-identical
  to the host loop at every format: one predictor and seven correctors
  at each of three formats, plus a per-lane-`dt` predictor for an
  adaptive ensemble, which is the 27 files in `programs/`.
- **The hardware.** docs/HARDWARE.md: the design keeps 47 values per
  coordinate resident; the projection is that a single system on the
  tile is 2-3x the software backend and an ensemble of a thousand is
  100-250x, that at binary64 the tile is three times slower than one
  CPU core running REBOUND itself, and that the scratch block - where
  the state lives - is staged every run under today's contract, which
  is the first thing to change. The first run on a tile then gave
  identical records to the software backend and was 4.5x slower on
  the two-body problem; the crossover, measured by the integrator, is
  width - the card wins at 496 pairs a call and loses at 28 - and
  difficulty alone does not move the ratio. The whole drop-in gate
  has since run on the quad tile at binary64: every case bit for bit
  against REBOUND's own IAS15, and every refusal refused, in 2,030 s
  against a few seconds in software. That ratio is the per-call cost
  at N = 2 and is exactly what the paragraph above predicts; what it
  establishes is that the backend does not change the answer, which
  is the claim the card exists to keep.
- **Ensembles.** E systems integrate in one run with their
  coordinates side by side in every vector, each with its own
  adaptive step, corrector exit and accept/reject decision, masked by
  snapshot where they part; the one-pair problem becomes E pairs.
  The gate: an ensemble reproduces, bit for bit, the records of its
  members run alone - 27 cases at three formats, fixed and adaptive,
  both engines, including an attempt in which four members reject a
  step while the fifth accepts it (tools/check_ensemble.py,
  docs/ENSEMBLE.md). The cost is idle lanes when members disagree
  about their corrector pass count: 1 to 31 percent, measured. The
  step's library-call count is independent of E (71,213 calls for 20
  binary256 steps at E = 1, 74,671 at E = 1,024), which is what the
  tile's per-call overhead needs; the software backend, element-bound,
  does 18 binary256 system-steps a second at E = 1 and 62 at 1,024 on
  one core, the baseline a card run must beat.
- **The horizon** (docs/HORIZON.md). Binary64's error is a random walk
  and a single run is one draw of it; a 64-member binary64 ensemble
  through REBOUND itself gives Brouwer's law to 0.02 dex, the
  90th-percentile phase error 2.1e-16 orbits^1.5, which reaches 1e-6
  of an orbit near 2e6 orbits at any setting (a billion-step single
  run landed at 8e-7 after 1e7 orbits). The wide formats' error is
  the method's, a t^2 ruler to 0.00 dex, 1.7e-21 orbits^2 at REBOUND's
  default tolerance and 3.9e-28 at 1e-12: binary128 reaches 1e-6 at
  2.4e7 orbits at the default and 5e10 at 1e-12, for 2.5 to 12 times
  the force evaluations. On Burrau's chaotic Pythagorean problem the
  binary64 solution is lost by t = 66 and the binary128 solution runs
  eighteen decades under it on the same steps, the two divergence
  curves being one curve 2^60 apart; the method's own error, measured
  separately, is nine decades under binary64's round-off at t = 60.
  An eight-member ensemble one to sixty-four ulps apart keeps its
  members in the ratios 2, 4, ... 64 through the whole evolution at
  binary256 and in no order at binary64: at binary64 an ulp-scale
  ensemble measures the arithmetic, at binary256 the dynamics.
  Eccentricity multiplies every floor by about fifty at e = 0.99.
  Binary64 is entirely sufficient for regular systems under a million
  orbits at 1e-6 of an orbit, for chaotic systems past their horizon
  where the science is statistical, and for anything limited by its
  physics; the niche is phases of a regular system past a few million
  orbits, close encounters passed with the energy still the method's,
  and results demanded exact rather than close.
- **What exactness permits.** Two formats differenced on identical
  steps measure a round-off floor rather than estimate it; an
  ensemble is gated bit for bit against its members; a recorded step
  sequence replayed reproduces its run to the bit. What it does not
  give IAS15 is an exact return: forward N steps and back N, no
  format comes home bit for bit (the Radau nodes are asymmetric), but
  binary128 and binary256 return to the same 7.8e-19, the method's
  own reversal error, where binary64 returns to 1.2e-13.
