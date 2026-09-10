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
