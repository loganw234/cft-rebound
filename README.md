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
binary256. The force model is **basic pairwise Newtonian gravity** -
REBOUND's `REB_GRAVITY_BASIC`, every pair, nothing approximated -
with IAS15's PRS23 adaptive step or a fixed one, up to **1024 bodies**
per system. At binary64 it is REBOUND's own IAS15 bit for bit, which
is a gate the port passes on every build. An ensemble of independent
systems integrates in one run, each with its own step
(docs/ENSEMBLE.md).

That is all of it. Everything below is refused, by name, with a
message - never silently ignored, never approximated:

| refused | |
|---|---|
| additional forces, `pre_`/`post_timestep_modifications` (REBOUNDx) | only gravity is ported |
| collision detection and resolution | not ported |
| ghost boxes, periodic or shear boundaries | not ported |
| the tree code (`REB_GRAVITY_TREE`) | this is direct summation. `REB_GRAVITY_COMPENSATED` is refused too: a different summation from the one ported |
| non-zero `softening` | not ported |
| test particles (`N_active`), `gravity_ignore_terms`, particle maps | they change which pairs are computed |
| variational particles | not ported |
| velocity-dependent forces | not ported |
| every integrator except IAS15 | WHFast is ranked first to follow (docs/INTEGRATORS.md) and has not been done |
| IAS15's `min_dt`, and adaptive modes other than PRS23 | not ported |

`cft_rebound_check()` is that table in code, and `cft_rebound_steps()`
calls it before doing anything, so a refusal arrives at the start of
your run rather than in the middle of it.

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

    include/cft_rebound.h    the public surface: what a REBOUND program includes
    src/cft_rebound_run.c    and what it links: the refusal list and the run call
    examples/roundtrip.c     the worked round trip, runnable (`make example`)
    examples/Makefile        a REBOUND program's makefile with the two added lines
    src/ias15_cft.c          the port: every floating-point operation is a cft.h call
    src/ias15_constants.h    GENERATED: the Gauss-Radau constants at every format
    src/hexfloat.h           exact hex-float text for binary64, libc-independent
    src/cft_ias15_state.h    the state the integrator shim and the archive share
    src/cft_archive.c        Simulationarchive: the cft_ fields, the probe, the load
    tests/cft_shim_stub.c    a stand-in integrator, for the archive gates only
    tests/gate_*.c           those gates: restart, stock reader, promotion, refusal
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
    ROADMAP.md               what is left before this is usable, and the two
                             REBOUND extension points that make it cheap

## Building

Plain C99 and GNU make. On Linux or macOS, with git, gcc and a Python
that has mpmath:

    make third-party        # clone and verify the pinned upstreams
    make libcft             # libcft from the pinned clone
    make                    # librebound (static), ias15_ref, ias15_cft
    make programs           # assemble the sequencer programs (needs cft-asm)
    make check              # the gates

On Windows the toolchain is MSYS2's mingw64 gcc from Git Bash, and the
same three traps cft-fp256's host/Makefile documents apply; pass them
as make variables:

    PATH="/c/msys64/mingw64/bin:$PATH" make CC=gcc OS=Windows_NT \
        TMP=C:/Users/you/AppData/Local/Temp TEMP=C:/Users/you/AppData/Local/Temp \
        PYTHON=C:/path/to/python.exe

REBOUND's own build assumes MSVC on Windows; here its library is
compiled with gcc from its sources minus `integrator_whfast512.c`,
whose private `__m512d` typedef collides with the one mingw's
`windows.h` brings in, plus a stub (ref/whfast512_stub.c). IAS15 does
not touch WHFast512. `-std=c99` pins `-ffp-contract=off`, so no FMA
contraction changes REBOUND's doubles.

## Using it from your own REBOUND program

    make install PREFIX=/usr/local        # DESTDIR=... also honoured

installs `cft_rebound.h` and `libcft_rebound.a` (plus libcft and the
`ias15_cft` program the library runs) under that prefix. Two lines go
into your own build:

    CFLAGS  += -I/usr/local/include
    LDLIBS  += -L/usr/local/lib -lcft_rebound -lcft

Then include `cft_rebound.h` instead of `rebound.h`, build your
simulation exactly as you do now, and call `cft_rebound_steps()` where
you would have called `reb_simulation_steps()`:

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

Today `cft_rebound_steps()` runs the `ias15_cft` program in a
subprocess and the wide state does not survive the call, so **ask for a
whole run in one call**: IAS15 starts each step's corrector from the
previous step's b coefficients, and a call that begins with them zeroed
lands on different bits - one to six ulps over twenty Kepler steps,
measured in docs/VALIDATION.md. ROADMAP parcel A replaces that
middle with `reb_integrator_register()`; the call above does not
change when it does, which is why the example is written against it.
Both `include/cft_rebound.h` and `examples/roundtrip.c` say which form
they are and what replaces it.

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
  of the software backend (untested here; the card was in use).
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

`make check` runs every gate at every format (about half an hour);
`make check-quick` runs them at binary64 in a few minutes. Both
include `tools/check_bodycount.py`, which re-runs the equivalence gate
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
- **The programs.** The predictor and corrector run as 24
  orbit-sequencer programs with the state in the per-lane scratch
  block, bit-identical to the host loop at every format.
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
  difficulty alone does not move the ratio.
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
