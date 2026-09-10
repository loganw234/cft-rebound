# Validation ledger

What ran, against which commits, and exactly what it said - passes and
failures alike, numbers verbatim from the runs. One entry per
campaign, newest last; a correction is a new entry naming the old one.
The style is cft-fp256's docs/VALIDATION.md, on purpose.

Host for every entry: DESKTOP-class Windows 11 box, MSYS2 mingw64 gcc
16.1.0, `-std=c99 -O2`, Miniconda Python 3.12.9 with mpmath 1.3.0
(on gmpy2), libcft from cft-fp256 at 40ce35e3 (ABI 0.11, software
backend only), REBOUND at bdfda4bd (5.1.1 + 6). No card, no XRT, no
`cft://` server: both were in use, and everything here is the software
backend, which is the contract's definition of the bits.

---

## 2026-09-09 - Phase 0: REBOUND's IAS15 in plain double, on Windows with gcc

REBOUND's `src/Makefile.defs` has three branches: Linux, Darwin, and
a Windows branch that assumes MSVC (`cl`, `cmd.exe`, `/Ox /fp:precise`).
Under Git Bash `uname` says `MINGW64_NT-...`, which lands in the
generic branch (`-std=c99 -Wpointer-arith -D_GNU_SOURCE -fPIC`, no
optimisation, `librebound.so`). Two things stop that build:

- `Cannot create temporary file in C:\WINDOWS\: Permission denied`
  unless `TMP`/`TEMP` are passed as make variables (the same trap
  cft-fp256's host build documents; a shell export does not reach
  gcc's subprocesses here).
- `integrator_whfast512.c:60:28: error: conflicting types for
  '__m512d'`: REBOUND defines its own `__m512d` struct so the file need
  not be compiled with AVX-512 flags, and mingw's `windows.h` (pulled
  in by `rebound_internal.h` for `winsock2.h`) includes `x86intrin.h`,
  which already defines it. MSVC's `windows.h` does not, so REBOUND's
  own Windows build never meets this.

So this repository builds REBOUND's library itself (root Makefile):
every `src/*.c` except `integrator_whfast512.c`, plus
`ref/whfast512_stub.c` so the built-in integrator table still resolves,
compiled `-std=c99 -O2 -D_GNU_SOURCE -DBUILDINGLIBREBOUND` (the last
because `rebound.h` marks its API `dllimport` on Windows otherwise, and
the library is linked statically) into `build/librebound.a`. Everything
IAS15 touches is in that archive unchanged; `-std=c99` pins
`-ffp-contract=off` so no FMA contraction alters the doubles.

`build/ias15_ref` reads a problem file of exact binary64 hex floats
(`data/problems/*.txt`, written by tools/make_problems.py - the outer
solar system's table copied by program out of REBOUND's own
`examples/outer_solar_system/problem.c` rather than retyped) and runs
`reb_simulation_steps` with IAS15 at `adaptive_mode = PRS23`, the
2024 default. Smoke run, Kepler (m0 = 1, m1 = 1/1000, a = 1, e = 1/2,
barycentric), 2,000 adaptive steps at epsilon 1e-9 from dt0 = 0.01:
t reaches 245.0 (39 orbits), the energy moves from
-0x1.0624dd2f1a9f4p-11 to -0x1.0624dd2f1aap-11, i.e. a few parts in
1e15, `iterations_max_exceeded = 0`. Native speed on this host,
Kepler at a fixed dt = 0.05: 100,000 steps in 0.92 s, 9.2 us a step.

## 2026-09-09 - Phase 1: the constants derived, and checked two ways

`src/integrator_ias15.c` hard-codes `h[8]` at 31 significant digits
and `rr[28]`, `c[21]`, `d[21]` at 25. binary128 needs 36 and binary256
73 (5.12.2's Pmin, which `cft_format_decimal_digits` reports as 17,
36, 73). tools/gen_constants.py derives them instead, at
`mp.dps = 130`:

- the seven non-zero nodes as the zeros of the exact polynomial
  `(P_7(x) + P_8(x)) / (1 + x)` (Legendre coefficients as Fractions,
  synthetic division with remainder checked to be 0, `polyroots` at
  extra precision, then eight Newton polishes on the exact
  polynomial), mapped to `h = (1 + x)/2`;
- `rr`, `c`, `d` by exactly the recurrence of REBOUND's own
  `GENERATE_CONSTANTS` block, in mpmath instead of GMP.

**Check 1, the defining polynomial.** All seven nodes satisfy the
quotient polynomial to below 1e-120. An independent definition was
tested rather than assumed: the Jacobi polynomial `P_7^(0,1)` vanishes
at the nodes with `max |value| = 5.8015e-130`, while `P_7^(1,0)` does
not (`3.4843`), which is the right way round for a Radau rule with its
fixed node at the left end.

**Check 2, REBOUND's published doubles.** Every one of the 78 values,
rounded to binary64 by this file's own integer rounding, equals the
literal REBOUND's source carries, bit for bit, and equals CPython's
independent `Fraction -> float` rounding of the same rational. The
31-digit `h` literals differ from the derived nodes by at most
1.18e-26, i.e. they were correctly rounded to their own length.

**A bug the check caught on its first run.** mpmath here runs on
gmpy2, whose mantissas are `mpz`; `Fraction` arithmetic refuses them
(`SystemError: Object does not appear to be Fraction`). Converting to
`int` fixed it; nothing was written until both checks passed.

Emitted: `src/ias15_constants.h` and `data/ias15_constants.json`, each
value as a canonical hex float per format (exact) and a 100-digit
decimal (so the C program can check that `cft_from_hex_char` and
`cft_from_decimal_char` round to the same bits at start-up - it does,
"constants: 78 derived values, hex and decimal forms agree" at every
format, later 155 with the FMA-form tables). Sample, `h[1]`:
binary64 `0x1.cce7242fd9812p-5`; binary256
`0x1.cce7242fd98120b32ee8a30ea77fda2d4227ed3a40ed4a94ce5a6a68bbb6f8p-5`.

## 2026-09-09 - the software backend's per-call cost, which sizes everything after it

A 20,000-call microbenchmark against `libcft.a` (scratch
`callbench.c`), one call at a time, n elements each, idle host:

| n | fp64 fma | fp64 add | fp64 div | fp64 sqrt | fp256 fma | fp256 div | fp256 sqrt |
|---|---|---|---|---|---|---|---|
| 3 | 0.40 us | 0.65 us | 16.2 us | 21.0 us | 0.90 us | 30.0 us | 40.4 us |
| 18 | 2.75 us | 3.55 us | 60.4 us | 87.6 us | 5.40 us | 211 us | 217 us |
| 1000 | 146 us | 200 us | 3,296 us | 4,352 us | 304 us | 7,248 us | 9,824 us |

So an elementwise call costs about 0.15 us an element at every n
above three (the fixed overhead is under half a microsecond), and a
composed divide or square root costs 25 to 30 elementwise passes, as
cft.h says it does. An IAS15 step that follows REBOUND's operation
order divides at every level of the predictor and the corrector; on
this backend the divides are half its time.

## 2026-09-09 - Phase 3 gate: the port at binary64 is REBOUND's IAS15, bit for bit

`src/ias15_cft.c` issues every floating-point operation through cft.h,
in REBOUND's order: `((b6*7)*h)/9` as three roundings, no FMA where
REBOUND has a product and a sum, Kahan's `add_cs` as its four
operations, the half-N^2 gravity summed partner by partner in the
order REBOUND's loop delivers them (vectorised over pairs, then one
add per accumulation round). The integer constants of the scheme are
derived in code - Taylor factors `(j+1)/(j+3)` reduced by their gcd
so the literal sequence matches REBOUND's `3/4`, `2/3`, `1/2`; the
binomial table of `predict_next_step`; the PRS23 weights `m+1`,
`m(m+1)`, `(m-1)m(m+1)`; the end-of-step divisors `(m+2)(m+3)` and
`m+2` - and every decimal (`1e-16`, `1e300`, `0.25`, `1e-7`, ...)
arrives through `cft_from_decimal_char`, which at binary64 gives the
bits a C literal gives.

tools/check_equivalence.py runs `ias15_ref` and `ias15_cft --format
fp64` on four cases and compares every recorded value as an exact
rational:

    kepler fixed dt=0.05:      17 samples, 272 values identical  (800 steps)
    kepler adaptive eps=1e-9:  17 samples, 272 values identical  (800 steps, 0 rejected)
    outer fixed dt=40:          9 samples, 360 values identical  (200 steps)
    outer adaptive eps=1e-9:    9 samples, 360 values identical  (200 steps, 0 rejected)
    RESULT: PASS

1,264 values - time, next and last step, energy, and 6 coordinates of
every body - identical on the first run of the full gate; the only
earlier failure was a compile error (a missing `n` argument in the
predictor, caught by the compiler). Those four cases do not exercise
the rejected-step path (`dt_new/dt_done < 0.25`) or REBOUND's cap of
12 corrector passes, because neither occurs in 1,000 steps of these
problems at a sane first step. Two more cases were added to force
them - a first step of 1.0 on the Kepler problem and of 2,000 days on
the outer solar system, each of which REBOUND rejects once and, on the
retry from the far too large step, runs the corrector to its cap
once (`steps_rejected=1, iterations_max_exceeded=1` in both
programs):

    kepler adaptive from dt0=1.0:   7 samples, 112 values identical
    outer adaptive from dt0=2000:   5 samples, 200 values identical

so the rejection path, the restore of the particles, the predictor
re-run with `ratio = dt/dt_last_done`, and the unconverged-step exit
are REBOUND's bit for bit as well. What remains unexercised is the
`ratio > 20` branch of `predict_next_step`, which the 4x growth clamp
makes unreachable in REBOUND's own flow; it is written to mirror
REBOUND and is untested.

Throughput at binary64 on the Kepler problem: 1,868 library calls a
step (255 of them div or sqrt), 120 steps a second on an idle host,
8.3 ms a step - 900 times REBOUND's own 9.2 us.

## 2026-09-09 - the corrector at wider formats: convergence per pass, and REBOUND's cap of 12

The predictor-corrector loop stops when `max|db6| / max|a| < 1e-16`.
Here that threshold is `1e-16 * 2^(53-p)` - the same number of ulps
at every format - and `--trace-pc` prints the error per pass. Kepler,
fixed dt = 0.05, the first step from a cold start (`b = 0`):

    pass   fp64        fp128       fp256
      1    1.92e-4     1.92e-4     1.92e-4
      2    1.93e-4     1.93e-4     1.93e-4
      3    6.40e-7     6.40e-7     6.40e-7
      4    2.80e-10    2.80e-10    2.80e-10
      5    8.77e-14    1.04e-13    1.04e-13
      6    0           5.88e-17    5.88e-17
      7                2.81e-20    2.81e-20
     ...                ...         ...
     10                4.45e-30    4.38e-30
     11                0           2.31e-33
     ...                            ...
     21                            4.06e-66
     22                            0

The three formats agree to the digits printed until each one's
round-off floor, where the error drops to exactly 0 (the iteration
reaches a fixed point of the finite-precision map). Each pass gains a
factor of about 500, 2.7 decimal digits, at this step size. Once the
predictor is warm (step 3 onward) a pass starts near 7e-7 and needs
about 2.4 passes at binary64, 6 at binary128 and 16 to 18 at
binary256. **REBOUND's cap of 12 passes is therefore hit on every
binary256 step** (`iterations_max_exceeded=100` of 100 in the first
probe) and the step proceeds unconverged; the sweep runs with
`--max-iter 60`, which no binary64 step reaches (`max_pc_iterations`
6 to 7 on the gate cases, under REBOUND's 12, so the gate results are
unaffected by the raised cap). This is the first cost of a wider
format in IAS15, and it is structural: the corrector converges
geometrically, so digits cost passes.

## 2026-09-09 - Phase 4 prototype: the predictor and corrector as sequencer programs, bit-identical to the loop

Two changes, each gated.

**The FMA form** (`--arith fma`). A correctly rounded divide is not a
sequencer program (cft-fp256 docs/ORBITS.md, obstacle 2), and REBOUND
divides at every level. The FMA form replaces each constant division
by one multiplication by an exact rational rounded once: `1/rr[i]` and
`h_n (j+1)/(j+3)` join the derived tables (155 hex/decimal pairs now
checked at start-up), `1/((m+2)(m+3))` and `1/(m+2)` come from
`cft_div` at run time, and the predictor collapses to one FMA per
level. It is a different rounding sequence and it differs from
REBOUND's at round-off: largest relative difference over 200 steps,
Kepler fixed dt = 0.05, `2.6e-15` at binary64, `9.0e-33` at binary128,
`2.1e-70` at binary256; on the outer solar system at binary64 the
60-step gate records were identical to the last bit - the differing
terms sit below the position's ulp for that long - and a 600-step run
shows the two forms parting from the first 100-step sample on; at the
wider formats the 60-step gap is `3.2e-34` and `9.1e-72`. The record
header names the form.

**Determinism, run to run and build to build.** The same binary256
adaptive Kepler run (30 steps, epsilon 1e-9, cap 60) twice on the
current build and once on the previous build (before the FMA form and
the flag fix existed) gives one SHA-256 over the sample lines,
`8f6d73f2...0ad697`, three times.

**The programs** (`--engine program`). tools/gen_programs.py emits 24
`.cfta` files - `predict-<fmt>` (21 instructions, a 9-entry bank:
seven `h_n (j+1)/(j+3)`, `fl(dt h_n)/2`, `fl(dt h_n)`; `x0, v0, a0` as
the three streams, `b0..b6, csx` from an 8-slot scratch block) and
`correct<n>-<fmt>` for n = 1..7 (22 to 106 instructions, a bank of
`2n-1` entries, a 22-slot scratch block in and out holding `a0, g, b,
csb`, one deposit) - and `cft-asm` from the pinned clone assembles all
24 (`correct7-fp256`: 106 instructions, 13 external constants, flags
`BANK_EXT SCRATCH_IO`, features `REGS32 BANK_PTR SCRATCH SCRATCH_IO`,
sha256 d96c1300...). Two assembler refusals on the way, both the
namespace rule docs/PROGRAMS.md warns about: a constant named `R0`
"would shadow a register", and a slot named `A0` collided with the
register alias `a0`.

tools/check_program_engine.py, three cases at each format, the program
engine against the host loop in the same FMA form:

    fp64:   kepler fixed 80/80, kepler adaptive 80/80, outer fixed 160/160 values identical
    fp128:  the same three, identical
    fp256:  the same three, identical
    RESULT: PASS

Library calls per binary256 Kepler step fall from 9,400 to 3,698 (the
remaining ones are the gravity, the end-of-step update, the step
control and scalar bookkeeping); on the software backend the program
engine is slower in wall time (3.3 against 13.3 steps a second), which
says nothing about the tile and everything about interpreting scratch
traffic in software.

## 2026-09-09 - a bug in the port's flag capture, found by an inconsistency between the two engines

Every record had printed `flags_seen=0x00`; the program engine's
records printed `0x10`. A standalone test showed `cft_run` reporting
inexact correctly all along. The wrappers were written
`note(cft_run(..., &fl, NULL), fl, "add")`, and C does not order the
evaluation of the two arguments - gcc read `fl` before the call wrote
it. Consequence: the four certificate flags (invalid, divide-by-zero,
overflow, underflow) could never have aborted a run of the loop
engine. Fixed by calling, then reading. The computed bits of every
record are unaffected (flags are outputs, not inputs); the
equivalence gate was re-run on the fixed build (`--quick`, 400 values)
and passes; every run since prints `0x10`, inexact and nothing else,
which is what a healthy IAS15 step must raise. The sweep below was
already running on the earlier build: its bits are the same and its
`flags_seen` field is wrong, and the re-run of its configurations on
the fixed build is recorded in the sweep entry.

## 2026-09-09 - Phase 2: does a wider format buy anything? Measured, both problems, three formats

The sweep (tools/sweep.py, sets `main`, `followup2`, and the outer
rows of `followup`): the port at binary64, binary128 and binary256 on
the two problems, every run scored by tools/oracle.py from the exact
bits of its record - the Kepler runs against the closed-form two-body
solution at the exact elapsed time, both problems' energy and angular
momentum evaluated at 60 digits. All 41 runs raised inexact and
nothing else; `--max-iter 60` throughout (REBOUND's 12 is hit on every
binary256 step, above). The box was running 11 to 20 of these at once,
so the `seconds` columns are relative to each other and about 2x an
idle host. tools/tabulate.py regenerates every table below from
results/; the per-sample CSVs are committed beside them.

### Kepler, e = 1/2, fixed step: max |dE/E| over the run, final |dx|/a, corrector passes per step

| dt | steps/orbit | orbits | fp64 | fp128 | fp256 | dx fp64 | dx fp128 | dx fp256 | passes fp64/128/256 | seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.1 | 63 | 200 | 3.04e-15 | 5.92e-17 | 5.92e-17 | 3.1e-12 | 9.6e-14 | 9.6e-14 | 3.3 / 7.9 / 18.3 | 265 / 1055 / 4339 |
| 1/16 | 100 | 100 | 1.28e-15 | 2.65e-20 | 2.65e-20 | 3.1e-13 | 2.2e-17 | 2.2e-17 | 2.7 / 6.8 / 16.1 | 186 / 856 / 3123 |
| 0.05 | 126 | 200 | 5.51e-15 | 1.88e-21 | 1.88e-21 | 8.4e-12 | 3.1e-18 | 3.1e-18 | 2.5 / 6.3 / 15.2 | 417 / 1910 / 6349 |
| 0.025 | 251 | 200 | 2.52e-15 | 5.79e-26 | 5.79e-26 | 3.5e-12 | 9.4e-23 | 9.4e-23 | 2.2 / 5.1 / 12.9 | 909 / 3232 / 8089 |
| 1/64 | 402 | 20 | 7.62e-16 | 5.06e-30 | 5.06e-30 | 8.2e-14 | 8.2e-28 | 8.2e-28 | 2.0 / 4.6 / 11.6 | 147 / 490 / 2009 |
| 1/128 | 804 | 10 | - | 4.55e-34 | 7.72e-35 | - | 1.6e-32 | 6.3e-33 | - / 3.8 / 10.1 | - / 353 / 1787 |

Read the binary128 and binary256 columns together: **they are equal
to every digit printed at every step size down to 1/64**, and only at
1/128 (804 steps per orbit) does binary256 fall below binary128, by a
factor of six. binary128 against binary64 is another matter: 5 orders
at dt = 0.1, 10 at 0.05, 11 at 0.025, 14 at 1/64.

The reason is in the step-size scaling. Halving the step from 0.05 to
0.025 divides the binary128 error by 3.2e4 (2^15 = 32,768: the method
is fifteenth order and this is it, measured), from 1/64 to 1/128 by
another 3.3e4 (5.06e-30 x 10/20 orbits / 7.72e-35). So at every step
size REBOUND's users would choose, the wide-format run is limited by
the method's own truncation error, which the format cannot touch, and
binary128's round-off floor is reached only at 800 steps an orbit.
binary256's floor is 2^-184 below that again, i.e. 55 more decimal
orders, which at dt^15 is a step 10^(55/15) = 4,600 times smaller
still. Nobody will take it.

**The round-off floors, measured directly.** Two records of the same
fixed-step run at two formats share their truncation error, so their
difference is the narrower format's round-off (tools/compare_formats.py;
these are the dyadic steps, so both formats integrate literally the
same dt - 0.05 rounds differently at each format and would put a
1e-16 time offset into the comparison):

| pair, step, horizon | max dE/E | final dE/E | max dx/L | final dx/L |
|---|---|---|---|---|
| fp64 - fp256, dt 1/16, 100 orbits | 1.28e-15 | 1.86e-16 | 6.2e-13 | 6.2e-13 |
| fp64 - fp256, dt 1/64, 20 orbits | 7.62e-16 | 1.98e-16 | 1.6e-13 | 1.6e-13 |
| fp128 - fp256, dt 1/16, 100 orbits | 1.34e-33 | 1.30e-33 | 1.3e-30 | 1.3e-30 |
| fp128 - fp256, dt 1/64, 20 orbits | 7.15e-34 | 2.67e-34 | 7.7e-32 | 7.7e-32 |
| fp128 - fp256, dt 1/128, 10 orbits | 4.81e-34 | 5.64e-35 | 4.4e-32 | 4.4e-32 |

The energy floor sits at a few ulps of the format (binary64's ulp of 1
is 2.2e-16, binary128's 1.9e-34) and does not grow much over these
horizons; the position floor grows roughly linearly with time, as a
phase error does. The two floors are 18 decimal orders apart, which is
2^(113-53) = 1.15e18.

**How the error grows** (results/*.csv). Where the wide format is
truncation-limited the energy error grows LINEARLY in time - binary128
at dt = 0.025: 1.16e-26, 2.32e-26, 3.47e-26, 4.63e-26, 5.79e-26 at 40,
80, 120, 160, 200 orbits, a straight line - and the position error
quadratically. That is a systematic error, the signature of a
non-symplectic method, and it does not average out: a wide format
turns IAS15's error from a random walk into a ruler. Where the run is
round-off-limited (binary64 at any of these steps) the energy error
wanders between 2e-17 and 5e-15 with no trend over 200 orbits and the
position error grows about linearly.

### Kepler, adaptive step (PRS23), about 195 orbits

| epsilon | steps | fp64 | fp128 | fp256 | dx fp64 | dx fp128 | dx fp256 | passes | seconds |
|---|---|---|---|---|---|---|---|---|---|
| 1e-9 (REBOUND's default) | 10,000 | 5.33e-15 | 4.59e-19 | 4.59e-19 | 6.7e-12 | 4.6e-16 | 4.6e-16 | 3.4 / 8.6 / 19.7 | 223 / 699 / 3822 |
| 1e-12 | 27,000 | 8.04e-16 | 1.02e-25 | 1.02e-25 | 5.6e-14 | 7.4e-23 | 7.4e-23 | 2.3 / 6.1 / 15.2 | 467 / 2389 / 6596 |
| 1e-16 | 100,000 | 4.67e-16 | 8.55e-34 | not run | 1.3e-13 | 8.0e-31 | | 2.0 / 4.2 / - | 2486 / 5336 / - |

This is the answer to the question as REBOUND's users meet it. **At
the default epsilon the binary64 run is round-off-limited by four
orders of magnitude** - 5.3e-15 against the 4.6e-19 the method
actually delivers at that step sequence - so binary128 buys those four
orders at once, for 2.5x the corrector passes and no change of step.
Tightening epsilon to 1e-12 (2.7x the steps) makes it ten orders,
8.0e-16 against 1.0e-25; at 1e-16 (10x the steps) binary128 reaches
its own floor, 8.6e-34, eighteen orders under binary64's 4.7e-16 at
the same step count, with the position error 8.0e-31 against 1.3e-13.
binary256 is again identical to binary128 at every epsilon it was run
at, and was not run at 1e-16 (100,000 binary256 steps is about six
hours here and the fixed-step rows already say what it would show).

### The outer solar system (Sun + inner planets, Jupiter, Saturn, Uranus, Neptune, Pluto as a test particle)

| dt (days) | steps | years | fp64 | fp128 | fp256 | dL/L fp64 | fp128 | fp256 | passes | seconds |
|---|---|---|---|---|---|---|---|---|---|---|
| 40 | 10,958 | 1200 | 1.21e-15 | 3.97e-31 | 3.98e-31 | 4.1e-16 | 4.0e-32 | 4.0e-32 | 2.2 / 6.0 / 15.2 | 493 / 2287 / 7274 |
| 20 | 21,915 | 1200 | 9.17e-16 | 7.77e-34 | - | 3.3e-16 | 3.1e-34 | - | 2.0 / 5.0 / - | 1103 / 4500 / - |
| 20 | 5,479 | 300 | - | - | 4.91e-36 | - | - | 2.8e-37 | - / - / 12.9 | - / - / 4323 |
| 10 | 10,958 | 300 | - | 2.16e-34 | 1.49e-40 | - | 7.3e-35 | 8.6e-42 | - / 4.0 / 11.0 | - / 2190 / 5458 |

and the floors measured directly: fp64 - fp128 at dt 40 (1200 y)
1.21e-15 in energy, 8.5e-14 in position; at dt 20, 9.2e-16 and
6.4e-14; fp128 - fp256 at dt 20 (300 y) 2.36e-34 and 2.4e-33; at dt 10
(300 y) 2.16e-34 and 8.0e-34.

**Here binary256 does pay.** This system's truncation error at a
given step count per orbit is far smaller than the eccentric Kepler
orbit's (nearly circular orbits, small accelerations), so binary128
hits its floor already at 20 days (7.8e-34 over 1200 years, 2.4e-34
measured against binary256 over 300) and binary256 goes under it: 48x
at 20 days, and at 10 days 1.49e-40 against binary128's 2.16e-34, six
orders of magnitude, with the fifteenth-order scaling intact
(4.91e-36 / 2^15 = 1.50e-40). At 40 days binary128 is still
truncation-limited at 4.0e-31 and binary256 lands on the same number, 3.98e-31 against 3.97e-31, the 6e-34 between them being binary128 own floor.

### Compensated summation: REBOUND's add_cs against 9.5 augmentedAddition

`--cs augmented` replaces Kahan's `cs = (t - p) - y` with the exact
error term of `cft_augmented_add`. Kepler, dt = 0.05, 200 orbits:

| format | cs | max dE/E | final dE/E | max dL/L | final dx/a |
|---|---|---|---|---|---|
| fp64 | kahan | 5.51e-15 | 2.07e-15 | 8.96e-16 | 8.4e-12 |
| fp64 | augmented | 1.69e-15 | 6.27e-16 | 4.70e-16 | 1.9e-12 |
| fp256 | kahan | 1.88e-21 | 1.88e-21 | 6.52e-23 | 3.1e-18 |
| fp256 | augmented | 1.88e-21 | 1.88e-21 | 6.52e-23 | 3.1e-18 |

At binary64, where the run sits on the round-off floor, the exact
primitive lowers it: 3.3x in the energy maximum, 4.4x in the final
position error, on this one run. Kahan's form is Fast2Sum, exact only
when the accumulator is at least as large as the increment, and IAS15
feeds it increments larger than the accumulator whenever a coordinate
crosses zero or a `b` coefficient starts from nothing; the augmented
form has no such condition. At binary256 the two runs' recorded
states are IDENTICAL at all 41 samples (tools/compare_formats.py: 0.0
throughout 25,120 steps), while the mean corrector pass count differs
in the fourth digit (15.171 against 15.169) - so a few compensation
terms did differ, and none of them reached a recorded coordinate,
because at this step the run is truncation-limited by twelve orders
of magnitude. One run per format; the binary64 factor is a single
sample, not a distribution.

### What it costs

The corrector's pass count is the price: 2.0 to 3.4 passes a step at
binary64, 4 to 9 at binary128, 10 to 20 at binary256 (each pass is a
force evaluation and a full b/g update), on top of the software
backend's per-operation cost, which is 1.3x from binary64 to binary128
and 2x to binary256 on this workload (docs/HARDWARE.md). The
`seconds` columns carry both: binary128 costs 3 to 5 times binary64
per step here, binary256 10 to 17 times. The step count needed to
reach a floor is the other price, and it is the one that decides:
binary128's floor needs 10x the steps of REBOUND's default on the
Kepler problem; binary256's would need 4,600x more again.

### The answer

For IAS15 on these problems: **binary128 is worth having and binary256
mostly is not.** binary128 removes four orders of error at REBOUND's
default settings and up to eighteen at ten times the steps, and turns
the error from a random walk into a linear drift whose slope is the
method's own; binary256 adds nothing over binary128 until the step is
small enough that a fifteenth-order method's truncation error falls
under 1e-34, which on the eccentric Kepler problem is 800 steps an
orbit and on the outer solar system 10 to 20 days (where it is worth
one to six orders more). The wide format does not extend the usable
horizon in the sense a symplectic method's does - the error still
grows, and it grows linearly - but it moves the floor that horizon is
measured against by 2^(p-53), which for binary128 over 200 orbits is
the difference between 1e-15 and 1e-33.

**Bit-for-bit, again, at length.** Three of the sweep's binary64
records were compared with REBOUND's own IAS15 run for the same steps:
the adaptive eps = 1e-9 Kepler run over 10,000 steps (41 samples, 656
values), the outer solar system at 40 days over 10,958 steps (32
samples, 1,280 values) and the fixed dt = 0.025 Kepler run over 50,240
steps (41 samples, 656 values): identical. Two binary64 jobs ran twice
by accident of scheduling, once on the sweep's original build and once
on the flag-fixed build: identical summaries to every digit. The
duplicates are why results/summary.txt has more lines than runs.

## 2026-09-09 - the shipped build, gated once more

After the sweeps, `build/ias15_cft` was rebuilt from the final source
and `make check` run on it: the constants derivation re-checked
(78 of 78 against REBOUND, nodes to 6e-130), the six binary64
equivalence cases (1,576 values identical, rejection cases included),
and the program engine at all three formats (nine cases, 960 values
identical). The results directory holds every per-sample CSV, the
direct floor measurements, results/tables.md as tools/tabulate.py
writes it, and results/energy_vs_orbits.png. results/raw/ (the hex
records, about 1 MB in all) is committed too; every number above is
recomputable from the commands in tools/sweep.py.


## 2026-09-09 - the first run on a tile

Everything above this entry was measured on the software backend,
because the card was in use while this port was written and the brief
forbade touching it. This is the port meeting hardware for the first
time, run by the integrator on the box that owns the card.

### What ran

The pinned tree, cloned from a bundle onto the Linux box that has XRT
and the Alveo U50. libcft rebuilt from the same pinned cft-fp256 clone
with `XRT=1`, and `ias15_cft` linked against that archive. The twenty-
four sequencer programs were assembled from the **committed** `.cfta`
files rather than regenerated, so what ran on the tile is exactly what
the software-backend gates checked.

The same integration twice, changing one argument:

      ias15_cft --format fp256 --arith fma --engine program
                --programs programs/out --problem kepler.txt
                --dt 0.05 --epsilon 0 --steps 200 --sample 20 --max-iter 60
      ... and again with --artifact ~/cardday-ra/cft_hw_quad.xclbin

The artifact is the read-ahead quad, four `cft_krnl` tiles at 135 MHz.

### The result

**The records are identical.** Every data line matches bit for bit -
eleven samples, each carrying the time, the next and last step, the
energy and every coordinate as an exact hex float. So does every
derived quantity in the trailer:

      steps_done                 200      200
      iterations_max_exceeded      0        0
      steps_rejected               0        0
      mean_pc_iterations      15.070   15.070
      max_pc_iterations           23       23
      flags_seen                0x10     0x10
      library calls          624,006  624,006
      divsqrt calls           45,860   45,860

Two lines differ in the whole file and neither is arithmetic: the
provenance header, which names the backend, and the wall clock.

That the predictor-corrector took the same number of passes on both -
15.070 on average, 23 at worst, over 200 steps of an adaptive scheme
whose iteration count is decided by a tolerance test on computed values
- is the part worth noticing. A single differing bit anywhere would
have moved it.

### The wall clock, which is the honest half

      software backend, one core     12.19 s    16.41 steps a second
      the tile                       55.37 s     3.61 steps a second

**The tile is 4.5x slower here, and that is the expected result rather
than a disappointment.** This problem is two bodies: six coordinates.
The engine's width is its whole advantage and six lanes leave nearly
all of it idle, while every program run still pays its fixed cost -
two runs per substep-pass, 624,006 library calls for 200 steps. The
projection in docs/HARDWARE.md said as much before this ran: one system
is 2-3x the software backend at best and the win needs an ensemble,
where a thousand systems put roughly a thousand coordinates across the
lanes.

So this run settles correctness on hardware and says nothing good about
throughput, which is what a two-body problem was always going to say.
The ensemble mode named as next step 3 is the one that would move it.

### What this does not cover

Single tile's worth of work on a four-tile artifact; no ensemble; one
problem; 200 steps; binary256 only; `--arith fma` and `--engine
program` only, since the program engine requires them. The three asks
of cft-fp256 in docs/HARDWARE.md are untouched and still asks: binding
`cft_alloc` buffers to `scratch_in`/`scratch_out`, a device-side
gather, and a `CFT_MAX` reduction.


## 2026-09-10 - the rewrite for ensembles, and the single-system path gated three ways

`src/ias15_cft.c` was rewritten so that E systems ride in one vector
(docs/ENSEMBLE.md). The single-system path goes through the same code,
so before anything else it had to be shown unchanged, and it was shown
three ways on the rebuilt binary:

- tools/check_equivalence.py, full: the six binary64 cases against
  REBOUND's own IAS15, 1,576 values identical, the two rejection cases
  included (`steps_rejected=1, max_exceeded=1` in both programs).
- tools/check_program_engine.py at all three formats: nine cases, 960
  values identical, and the FMA form's gap from REBOUND's rounding
  sequence at `2.625e-15`, `9.004e-33` and `2.092e-70` on the fixed
  Kepler case - the same digits the 2026-09-09 entry printed, which is
  its own regression check.
- tools/check_records.py, new: the first two samples after the start
  of nine committed sweep records (`results/raw/`: Kepler fixed and
  adaptive and the outer solar system at binary64 and binary128, Kepler
  fixed and adaptive at binary256, and the augmented-summation run)
  recomputed by the new binary on the same sampling grid and compared
  as exact rationals: 27 samples, 576 values identical. A record's
  sampling does not perturb its trajectory, so a short run reproduces
  a long run's early samples exactly, and this is what protects the
  E = 1 path when the integrator changes shape.

Two things changed on purpose in the E = 1 path and are visible only
in the call count, not in any recorded value: the energy is now
vectorised over bodies and pairs (the same operations in REBOUND's
order, issued once for all systems), and `sqrt7(epsilon * 5040)`,
whose input never changes, is computed once instead of once per step
(20 Newton iterations, about 200 scalar calls a step). The
`calls=` trailer of a Kepler run therefore no longer reads what the
2026-09-09 entries report; the data lines do.

One bug was caught by the new binary before any gate ran: the
vectorised energy first multiplied `G * m_smaller * m_larger` where
`reb_simulation_energy` multiplies `G * m_j * m_i` with `j` the larger
index. Identical for G = 1 and irrelevant for the Kepler gates; not
identical for the outer solar system's `G = k^2`, where the quick
equivalence gate would have failed on the energy column. Found by
reading, fixed before building.

## 2026-09-10 - the ensemble mode and its gate

tools/check_ensemble.py (docs/ENSEMBLE.md): each case runs an ensemble
once and every member alone with `--member k` on the same file with
the same arguments, and compares every recorded value of every sample
- time, next and last step, energy, every coordinate, the exact-time
pair - as an exact rational, plus each member's own trailer (steps
done and rejected, corrector passes mean and maximum, cap hits). The
families: five Kepler orbits with the planet's velocity scaled by
1 - 0.08k (eccentricities 0.5 down to 0.03, different periods and
step sequences); four Kepler orbits 2^-40 apart in x; five Kepler
orbits with the planet started at 0.5, 0.9, 1.3, 1.7 and 2.1 with the
base orbit's pericentre speed (an 8.6x spread in natural step, the
outer two unbound) from a first step of 0.6, which the inner four
reject in the same attempt in which the outermost accepts; and three
outer solar systems with Jupiter displaced by 2^-20 AU. Fixed and
adaptive steps, the loop engine and the program engine (the latter
with the new `predict-ens` program when the steps are per system).

At binary64, `--quick` (50 to 30 steps a case, 37 s in all), nine
cases:

    kepler vscale E=5 fixed dt=0.05:                 5 members, 180 values identical, rejected [0,0,0,0,0], lane efficiency 0.766
    kepler vscale E=5 adaptive eps=1e-9:             5 members, 180 values identical, rejected [0,0,0,0,0], 0.868
    kepler spread E=5 adaptive from dt0=0.6:         5 members, 180 values identical, rejected [1,1,1,1,0], 0.862
    kepler spread E=5, the same, program engine:     5 members, 180 values identical, rejected [1,1,1,1,0], 0.836
    kepler ulps E=4 adaptive eps=1e-9:               4 members, 144 values identical, rejected [0,0,0,0], 0.919
    kepler vscale E=5 adaptive eps=1e-9, program:    5 members, 180 values identical, 0.871
    kepler vscale E=5 fixed dt=0.05, program:        5 members, 180 values identical, 0.806
    outer jupiter E=3 fixed dt=40:                   3 members, 252 values identical, 0.927
    outer jupiter E=3 adaptive eps=1e-9:             3 members, 252 values identical, 0.974
    RESULT: PASS

The first run of the gate failed one case by design of the case, not
of the code: the mixed-rejection family was first the velocity-scaled
one from dt0 = 0.6, and all five members rejected the step
(`steps_rejected=[1,1,1,1,1]`), so the attempt was not mixed and the
gate said so. The spread family replaced it; the trace above is the
second run. `pc_lane_efficiency` is the fraction of issued lane-passes
that belonged to a system still iterating: with 2 to 6 passes a step
at binary64 the members disagree often enough to idle a quarter of
the lanes on the fixed-step Kepler family and 3 percent on the outer
solar system, whose members are nearly identical.

The full gate, 200 to 60 steps a case, at all three formats, 27 cases
in all (59 minutes on a box running ten other jobs):

    fp64:   nine cases, 3,708 values identical; lane efficiency 0.750 to 0.972
    fp128:  nine cases, 3,708 values identical; 0.721 to 0.984
    fp256:  nine cases, 3,708 values identical; 0.687 to 0.996
    RESULT: PASS

with the mixed-rejection family rejecting `[1, 1, 1, 1, 0]` at every
format on both engines, and the corrector passes per step rising from
2.3-3.5 at binary64 through 6.0-8.9 at binary128 to 13.7-20.8 at
binary256. The lane efficiency is what the shared corrector loop
costs: at binary256 the members of the velocity-scaled family
disagree enough about their pass count to idle 18 percent of the
lanes on the fixed step and 6 percent on the adaptive one; the
spread family, whose members range from a bound e = 1/2 orbit to an
unbound one, idles 31 percent; the ulp family and the outer solar
system, whose members are nearly the same system, idle under 1
percent. So an ensemble of a chaotic system's nearby copies - the
scientific case - wastes almost nothing, and an ensemble spanning
regimes wastes a third.

## 2026-09-10 - prescribed steps: a replay is the run, bit for bit

`--dt-file FILE` feeds a run its step sequence (one exact hex float
per step) and `--dt-out FILE` records one; the replay sets the next
step before the prediction ratio `dt / dt_done` is formed, exactly as
the adaptive path does after an accepted step. The consequence is a
checkable property: an adaptive run that rejected no step, replayed
from its own recorded sequence, must reproduce its record to the bit.
Kepler and the Pythagorean problem, 2,000 adaptive steps at
binary64 (`steps_rejected=0` in both), replayed:

    kepler:       2001 dt lines; replay IDENTICAL on all 21 sample lines
    pythagorean:  2001 dt lines; replay IDENTICAL on all 21 sample lines

The first version of `--dt-out` wrote `steps` lines - one per accepted
step - and the replays differed in exactly one value: the last
sample's `dt_next`, because the adaptive run had already chosen its
2,001st step and the replay's sequence had ended. Every position,
velocity, time and energy was identical. `--dt-out` now writes the
step the run would take next as its last line (`steps + 1` in all),
and the difference is gone. The point of the mechanism is the horizon
campaign below: a binary256 run's sequence, rounded to binary64, lets
every format take literally the same steps, so that the difference
between two formats' records is the narrower one's round-off and
nothing else.

## 2026-09-10 - the reversal test: no format returns exactly, and the wide ones return to the same error

tools/reversal.py: N fixed steps forward, the exact recorded state
written as a problem file (tools/state_to_problem.py), N steps back
with `-dt`, and the returned state against the initial one as exact
rationals. Kepler at dt = 1/16, 1,000 steps out and back (10 orbits
each way); the outer solar system at 40 days, 500 steps each way:

    problem  format  bits identical  max|dx|/L   max|dv|/V   |dE/E|
    kepler   fp64    no              1.243e-13   8.254e-14   1.431e-16
    kepler   fp128   no              7.813e-19   5.207e-19   5.040e-21
    kepler   fp256   no              7.813e-19   5.207e-19   5.040e-21
    outer    fp64    no              7.893e-16   6.972e-15   2.489e-16
    outer    fp128   no              1.162e-31   9.595e-31   2.106e-31
    outer    fp256   no              1.154e-31   9.531e-31   2.104e-31

IAS15 is not time-symmetric - its Gauss-Radau nodes include the start
of the step and not the end, so the backward polynomial is a different
polynomial - and no exact arithmetic would bring it back bit for bit;
the test was run to see which error is which, and it does. On Kepler
the binary128 and binary256 returns are identical to every digit
printed: 7.813e-19 of the orbit is the method's own asymmetry at this
step, with the arithmetic invisible under it, and binary64's
1.243e-13 is a factor of 1.6e5 above it and is all arithmetic. On the
outer solar system, whose 40-day step is so smooth that the method's
asymmetry is 1.15e-31, binary128 sits 8e-34 above binary256 - its own
floor showing - and binary64 is 7e15 times higher. A bit-exact return
is a property of a symmetric or a lattice scheme (JANUS, ranked for it
in docs/INTEGRATORS.md), not of this one; what is exact here is the
run itself, repeated.

## 2026-09-10 - the horizon of a regular orbit at binary64, measured with REBOUND itself

The question (docs/HORIZON.md): after how much simulated time does a
binary64 answer stop being worth having. Every binary64 run in this
entry is REBOUND's own IAS15 through `build/ias15_ref`, which the
equivalence gate shows is the port at binary64 bit for bit and which
does a Kepler step in 5 to 9 us here; `ias15_ref` now writes an exact
TwoSum time pair beside REBOUND's `r->t` (a label, not an input: the
force and the step control are time-independent, and the quick
equivalence gate passed unchanged), because `r->t` alone mislabels a
1e8-step run at the 1e-5 level, which is the size of the effect being
measured. The error is the along-track phase error of tools/oracle.py
in units of the period; the oracle's Kepler solver was made robust on
the way (below).

**Single runs.** Kepler e = 1/2 at the fixed dyadic step 1/16
(100.48 steps an orbit, the time exact by construction):

    orbits   steps         max |dE/E|   phase error (last sample)
    1e4      1,004,807     1.109e-14    2.09e-11
    1e5      10,048,074    5.027e-14    1.879e-9
    1e6      100,480,737   1.271e-13    3.486e-8

and at REBOUND's default adaptive settings (epsilon 1e-9, 1e8 steps
each, 12 to 17 minutes a run):

    e     steps/orbit  orbits      max |dE/E|   phase (max / last)
    0.5   51.2         1,953,174   6.546e-13    9.857e-7 / 9.857e-7
    0.9   97.2         1,029,258   5.089e-13    8.944e-8 / 6.554e-9
    0.99  159          627,838     6.824e-12    1.865e-6 / 1.792e-6

The e = 0.99 run first reported a maximum position error of 8.3
semi-major axes, which is geometrically impossible for a bound orbit
(two points on it are at most 4a apart): the oracle started Newton on
Kepler's equation at dE = M with M ~ 4e6 radians and e = 0.99 and did
not converge. The solver now reduces the mean anomaly modulo 2 pi at
the working precision and runs Newton inside a bisection bracket on
the monotonic function; re-scoring the committed
`kepler_fixed_dt0.05_fp64` record reproduces its committed CSV in
every column of every sample but one, the first sample's position
error, which reads 1.68e-72 instead of 0.0 (the 60-digit reduction's
own round-off). With the robust solver the e = 0.99 run's maximum
`|dx|/a` is 1.643e-4 and its phase error 1.865e-6 of an orbit.

**Single runs do not measure the law.** Fitting `phase = C t^alpha`
to those runs gives alpha = 0.59, 1.38 and 1.30 on the three
fixed-step spans, 1.92 for the adaptive e = 1/2 run (0.06 dex
residual, on top of an energy error growing as t^0.85 - a drift, not
a walk, in that realisation), 0.75 for e = 0.9 and 1.90 for e = 0.99,
with 0.2 to 0.5 dex of scatter in most; the e = 0.9 run's phase error
at a million orbits (6.6e-9) is five times SMALLER than the e = 1/2
run's at the same count. A round-off random walk is a distribution
and one run is one draw from it.

**Ensembles measure it.** tools/make_ensemble.py wrote 64 copies of
the e = 1/2 problem with the planet's x shifted by 2k ulps
(`--offset planet x 0x1p-52`) and their member files; each was run
through `ias15_ref` for 1e7 steps at epsilon 1e-9 (2e5 orbits) and
again at the fixed step 1/16 (1e5 orbits), 128 runs, four at a time,
and tools/horizon.py took the distribution at every sample:

    adaptive, epsilon 1e-9, 64 members, phase error:
      orbits    median      p90         max         min
      1953      1.042e-11   2.350e-11   3.723e-11   1.259e-13
      2.15e4    3.458e-10   7.513e-10   1.268e-9    8.371e-12
      6.06e4    1.793e-9    4.151e-9    7.340e-9    1.725e-10
      1.19e5    4.410e-9    1.254e-8    1.582e-8    8.571e-12
      1.78e5    8.334e-9    1.995e-8    3.169e-8    3.233e-10
      p90    = 2.110e-16 * orbits^1.525   (90 samples from 2e4 orbits, rms residual 0.02 dex)
      median = 2.073e-16 * orbits^1.446   (0.02 dex)
      energy p90 1.489e-14 at 1953 orbits, 1.314e-13 at 1.78e5: ~ orbits^0.5
      the p90 curve crosses 1e-9 of an orbit at 2.73e4 orbits (the fit says 2.39e4)

    fixed dt = 1/16, 64 members, phase error:
      p90    = 3.212e-16 * orbits^1.471   (0.01 dex)
      median = 9.761e-17 * orbits^1.494   (0.02 dex)
      energy p90 = 1.729e-16 * orbits^0.540 (0.03 dex); median ~ orbits^0.477
      crosses 1e-9 at 2.59e4 orbits (fit 2.6e4)

That is Brouwer's law - the energy error a t^1/2 walk, the phase error
its t^3/2 integral - measured to 0.01-0.02 dex, and the fixed-step
coefficient is 1.5x the adaptive one where 100 steps an orbit against
51 is a walk 1.4x longer. The single long run (2e6 orbits, 9.86e-7)
lands where the adaptive p90 law says 1.1e-6. One more single run
went a decade further: the fixed step 1/16 for 1e7 orbits,
1,004,807,371 steps, 2.6 hours of REBOUND, phase error 3.49e-8 at
1e6 orbits, 3.24e-7 at 2e6, 6.07e-7 at 3e6 and 8.19e-7 at 1e7 (the
walk stalled between 3e6 and 9e6 orbits, as walks do), against the
fixed-step median law's 2.6e-6 and p90's 6.6e-6 at 1e7: a draw on the
low side, within the ensemble's spread (its members differ by a
factor of 100 at every sample). The energy error of that run never
left the 2e-14 to 3e-13 band in a billion steps. Extrapolated on the
law beyond the data, the p90 phase error reaches:

    threshold   adaptive eps 1e-9   fixed dt 1/16
    1e-6        2.21e6 orbits       2.86e6
    1e-3        2.05e8              3.13e8
    0.1         4.2e9               7.2e9

The extrapolation rests on the law, not on the data, past 2e5 orbits
(2e6 for the one long run); the law is Brouwer's and the fit is clean,
but it is an extrapolation and is marked as one in docs/HORIZON.md.

**The other law.** Re-scoring the committed sweep records with the
phase column: every truncation-limited wide-format record has energy
error proportional to t^1.00 and phase error to t^2.00, both to 0.00
dex over 200 orbits - the systematic per-orbit error of the method,
integrated twice - with coefficients set by the step control:

    binary128 (= binary256 to every digit), phase = C * orbits^2:
      adaptive eps 1e-9:   C = 1.71e-21     (energy 4.589e-19 at 195 orbits, ~ t^1.00)
      adaptive eps 1e-12:  C = 3.86e-28     (energy 1.025e-25, ~ t^1.00)
      adaptive eps 1e-16:  C = 1.92e-36, alpha 2.07 (energy 1.57e-34, flat: the binary128 floor)
      fixed dt 0.05:       C = 7.04e-24;  dt 1/16: 1.99e-22;  dt 0.025: 2.17e-28;  dt 1/64: 1.90e-31

Against binary64's 2.11e-16 t^1.5 walk, the two laws cross at 1.5e10
orbits at epsilon 1e-9: binary128 at REBOUND's default tolerance
reaches 1e-6 of an orbit at 2.4e7 orbits (eleven times binary64's
2.2e6) and 0.1 at 7.6e9 (1.8 times binary64's 4.2e9), because the
method's t^2 overtakes the arithmetic's t^1.5; at epsilon 1e-12 it
reaches 1e-6 at 5.1e10 orbits, 23,000 times later than binary64 can
at any setting. The table with costs is in docs/HORIZON.md, and the
sentence it supports is: binary64's horizon is fixed by the
arithmetic and no setting moves it; binary128's is set by the method
and is for sale at 2.5 to 12.5 times the force evaluations an orbit.

**The outer solar system**, fixed 40-day step, binary64 through
`ias15_ref`: 1e6 years (9,131,250 steps) max |dE/E| 1.878e-14, max
|dL/L| 6.072e-15; 1e7 years (91,312,500 steps, 16 minutes) 5.091e-14
and 1.799e-14, the energy error between 2e-14 and 5e-14 throughout
and fitted as t^0.63 with 0.5 dex of scatter - one draw of a walk
again. No closed form, so no phase; the sweep's binary128 record at
the same step gives 3.97e-31 over 1,200 years, a t^1.00 ruler, so at
40 days the two formats are 16 orders apart in energy at 1,200 years
and, on their laws, 15 orders at 1e7 years.


## 2026-09-10 - the chaotic horizon: Burrau's Pythagorean problem at three formats on one step sequence

Burrau's problem (masses 3, 4, 5 at rest at the vertices of a 3-4-5
triangle, G = 1; `data/problems/pythagorean.txt`, every value an
integer and exact in every format) is the classic sensitive three-body
system: repeated close encounters, an adaptive step spanning four
orders of magnitude, and near t = 60-70 the formation of a 4-5 binary
that ejects body 3. The port takes it in its stride - a binary256
adaptive run at REBOUND's default epsilon 1e-9 from dt0 = 0.01 did
6,000 steps to t = 73.4 in 29.6 minutes on a loaded box, 16.94
corrector passes a step (22 at worst), no step rejected, inexact the
only flag, its step falling to 5e-6 at the closest passage - but the
question is what its arithmetic is worth against the dynamics, and
that needs three formats on *literally the same steps*: that run's
sequence, rounded to binary64 (`--dt-out`, tools, then `--dt-file`),
replayed at binary64, binary128 and binary256, sampled every 20
steps. Every replay then has the same 301 sample times to the bit and
the difference between two of them is the narrower format's round-off,
amplified by the dynamics, and nothing else.

**The two divergence curves are one curve, 18 decades apart.** The
largest relative coordinate difference (`|dx|/L`, tools/
compare_formats.py) first exceeds

    threshold  binary64 - binary256   |  threshold  binary128 - binary256
    1e-12      t = 6.89  (step 360)   |  1e-30      t = 6.89  (step 360)
    1e-9       t = 22.96 (step 1340)  |  1e-27      t = 22.96 (step 1360)
    1e-6       t = 46.00 (step 2480)  |  1e-24      t = 46.01 (step 2500)
    1e-3       t = 59.78 (step 3560)  |  1e-21      t = 60.63 (step 3700)
    1e-2       t = 63.21 (step 4220)  |
    1          t = 66.03 (step 4740)  |  (3.1e-21 at t = 73.4, the end; 3.5e-20 at most)

The binary128 curve crosses each threshold 1e-18 = 2^-59.8 lower at
the same time as the binary64 curve, to within one sample: the seed
is the format's floor, 2^-53 against 2^-113, and the amplification is
the problem's and the same for both. It is a staircase - flat between
encounters, a jump at each: 1e-12 at t = 7, 2e-11 at t = 16 after the
first close passage (where the step is 7e-6), 7e-10 at t = 20, 4e-9
at t = 28-34, 7e-8 at t = 42, 6e-5 at t = 50, 1e-3 at t = 60 - which
averages a decade of divergence every four time units, an e-folding
time of about 1.7. **The binary64 solution of this problem is wrong at
the 1e-6 level by t = 46, at the 1e-3 level by t = 60, and is unrelated
to the true trajectory from t = 66, before the ejection that is the
problem's whole point; the binary128 solution reaches the end 3e-21
from binary256's**, with 18 decades of the same staircase still to
climb - about 70 more time units at the measured rate, if the
dynamics stayed chaotic that long, which they do not once the binary
has formed and body 3 has gone.

**The replay artefact, and why it is not the integrator.** From
t = 66 the binary64 replay's energy error is 698 (yes, 698: the
relative error, not its logarithm) and its distance from binary256 is
124 length units. That is the replay, not binary64 IAS15: the
prescribed sequence was chosen for the binary256 trajectory, and once
the binary64 trajectory has left it by 1e-2 (t = 63) an encounter
arrives where the sequence has a large step, and a fifteenth-order
method through a close passage at a step meant for free flight is
destroyed. REBOUND's own adaptive binary64 run of the same problem
(`ias15_ref`, 200,000 steps to t = 1,142, 0.9 s) keeps its energy to
3.0e-9 throughout. The honest reading is the row above it: the
binary64 solution is worthless from t = 60, whichever way the steps
are chosen.

**What an encounter costs each format.** The energy error of the
replays against the exact initial energy (tools/oracle.py, which now
reports the angular momentum as an absolute error when the initial
one is exactly zero, as it is for a system started at rest):

    t        binary64      binary128      binary256
    3.8      4.05e-14      1.73e-20       1.73e-20
    11.9     3.17e-13      1.91e-20       1.91e-20
    15.8     3.50e-11      3.62e-20       3.62e-20      <- the first close passage
    19.8     7.15e-11      4.06e-20       4.06e-20
    41.9     7.19e-11      5.36e-20       5.36e-20
    59.8     7.17e-11      1.50e-19       1.50e-19
    73.4     (replay lost) 6.23e-19       6.23e-19

One close passage costs binary64 three orders of magnitude of energy
(4e-14 to 7e-11 in one step); the same passage costs binary128 2e-20
to 4e-20, and the binary128 column is the binary256 column to every
digit printed - the round-off difference between them is 3e-29 - so
at binary128 the energy error of this problem is entirely the
method's truncation at epsilon 1e-9, and the wider format's
arithmetic is invisible in it. The binary256 adaptive run's own
energy history (its unrounded steps) is the same to four digits at
every common step: rounding the sequence to binary64 changed nothing
that matters.

**What this does not measure.** The truncation seed is common to all
three replays and cancels in every difference above. At epsilon 1e-9
it is about 1e-19 in energy, three decades below the binary64
round-off floor, so by the staircase's rate it should be lost about
twelve time units after binary64 is - the horizon of *any* format at
this tolerance would then be near t = 78, and the eighteen decades
binary128 has in hand are only for sale against a tighter step
control (the sweep: epsilon 1e-16 is ten times the steps). A direct
measurement of the truncation-seeded error at one time is in the
entry below.

**The ensemble as the instrument.** tools/make_ensemble.py wrote eight
copies of the problem with body 3's x displaced by 1, 2, 4, ... 64
binary64 ulps (`--geometric m3 x 0x1p-52`; member 0 untouched) and
they were integrated in one run on the same replayed sequence, at
binary64 and at binary256, tools/divergence.py measuring each member
against member 0 at every sample. At binary64 the members reach
1e-12 by t = 7 in no order - m6 (32 ulps) below m1 (1 ulp) at t = 20,
and all seven within a factor of four of each other from t = 7 on:
the perturbation is the size of the round-off noise and the ensemble
measures the noise. The binary256 run of the same file, where the
same perturbations are 2^184 times the round-off, is the measurement
that ensemble was made for; it is recorded in the entry below when it
completes (a 72-lane binary256 run on the software backend is
element-bound and takes hours).
