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


## 2026-09-10 - the truncation-seeded error of the chaotic problem, measured at one time; a correction to the entry above

The entry above inferred that the truncation seed at epsilon 1e-9,
"about 1e-19 in energy, three decades below the binary64 round-off
floor", would be amplified like the round-off seeds and lost about
twelve time units after binary64 is, near t = 78. That inference was
wrong, and the measurement that shows it is this one.

**Method.** The replay mechanism compares formats at identical
times but cannot compare step controls, whose sequences differ. Two
sequences can, however, be made to end at the same time: an fp256
adaptive run at epsilon 1e-11 (12,000 steps to t = 76, 14.4 corrector
passes a step, none rejected, 41 minutes) recorded its sequence; that
sequence and the epsilon 1e-9 one were each rounded to binary64 and
cut at T* = 60, the final step of each computed exactly (as a
rational) as the remainder and rounded once, so that both end at
t = 60 to within 5e-19 (`T* - 1.7e-19` and `T* + 4.9e-19`); both were
replayed at binary256 (3,650 and 7,049 steps), where the arithmetic's
own contribution is 1e-70 and irrelevant, and the epsilon 1e-9
sequence at binary64 as well. The epsilon 1e-11 run's truncation
error is (1e-11/1e-9)^(15/7) = 1.9e4 times smaller than the epsilon
1e-9 run's, so their difference at T* is the epsilon 1e-9 run's
truncation-seeded error to five parts in 1e5.

**Result**, the largest relative coordinate difference at t = 60.000:

    epsilon 1e-9 at binary256   against epsilon 1e-11 at binary256:   1.333e-13
    epsilon 1e-9 at binary64    against epsilon 1e-11 at binary256:   3.412e-4

(the two final times differ by 6.6e-19; the binary64 number agrees
with the staircase of the entry above, 1.27e-3 at t = 59.78 and
4.39e-4 at t = 60.76.) The truncation-seeded error at REBOUND's
default tolerance is nine decades *below* the binary64 round-off-
seeded error at the same time, not three decades above it.

**Why the inference was wrong.** The round-off seed is random and has
components across the flow, which the encounters amplify; a
fifteenth-order method's truncation error is smooth and lies mostly
along the flow - a time shift on the trajectory, which is why the
energy stays at 1e-19 while the round-off-seeded position error
reaches 1e-3 - and a shift along the trajectory is a neutral
direction that chaos does not amplify. The sweep's "four orders above
binary64's floor" was an energy statement; the divergence is a
transverse one. One time was measured, not a curve, so the growth
law of the truncation-seeded error is not known from this; what is
known is its size at the time binary64 is lost.

**What it changes in the answer.** At REBOUND's default tolerance,
on this problem, the binary128 solution is not merely "binary256's
to 1e-21": it is nine decades closer to the true trajectory than the
binary64 solution at the moment the binary64 solution becomes
worthless, with no change of step control. The eighteen decades of
arithmetic headroom are not all needed to buy that; the first nine
come free. Tightening the step control buys the rest, if the science
needs a trajectory past that point. docs/HORIZON.md is corrected
accordingly, and the entry above stands as written with this entry
beside it.


## 2026-09-10 - high eccentricity: the floors of every format at e = 0.99, on one step sequence

The sweep measured the round-off floors on the e = 1/2 orbit with
fixed dyadic steps. A highly eccentric orbit cannot take a fixed step
(its pericentre passage is 200 times faster than its apocentre), so
the replay mechanism was used instead: an fp256 adaptive run at
REBOUND's default epsilon 1e-9 records its sequence, the sequence is
rounded to binary64, and all three formats take literally the same
steps. Two seeds of 20,000 steps each, `data/problems/kepler_e099.txt`
(e = 0.99, dt0 = 0.001; 17.24 corrector passes a step, 42 at worst at
the pericentre, one step rejected - the first - in 71 minutes) and
the e = 1/2 problem (19.72 passes, none rejected), then six replays.

    e = 0.99, 20,000 steps = 125.8 orbits, 159 steps an orbit:
      fp64  - fp256   max |dE/E| 1.694e-13 (last 1.671e-13)   max |dx|/L 1.228e-7 (last 7.441e-9)
      fp128 - fp256   max |dE/E| 9.719e-32 (last 8.110e-32)   max |dx|/L 3.337e-26 (last 2.318e-27)
      against the closed form: fp64 phase error 1.71e-11 orbits; fp128 and fp256 2.535e-16, identical to every digit,
      energy 2.690e-18 at both, |dx|/a 1.828e-14 at both
    e = 1/2, 20,000 steps = 390.7 orbits, 51.2 steps an orbit:
      fp64  - fp256   max |dE/E| 3.748e-15 (last 2.950e-15)   max |dx|/L 4.407e-12 (last 1.280e-12)
      fp128 - fp256   max |dE/E| 2.325e-33 (last 1.068e-33)   max |dx|/L 9.464e-30 (last 4.356e-30)
      against the closed form: fp64 phase error 1.39e-13 orbits; fp128 and fp256 2.694e-16, identical, energy 9.201e-19

Three things to read off. **The floors scale as the format says at
either eccentricity**: binary64 over binary128 is 1.7e18 = 2^60.6 in
energy and 3.7e18 = 2^61.7 in position at e = 0.99, 1.6e18 = 2^60.4
and 4.7e17 = 2^58.7 at e = 1/2. **Eccentricity multiplies both floors
by the same factor**, about fifty in energy (1.69e-13 against
3.75e-15 at binary64, 9.7e-32 against 2.3e-33 at binary128) with
three times fewer orbits in the eccentric run - the amplification is
the orbit's ((1+e)/(1-e) = 199 between the kinetic and potential
terms at pericentre and the total), not the format's. And **the
method's own error grows with eccentricity as well**: the truncation
ruler of the fp256 seeds is `phase = 1.57e-20 orbits^2` at e = 0.99
against `1.72e-21 orbits^2` at e = 1/2, both t^2.00 to 0.00 dex,
nine times, at three times the steps an orbit. The e = 1/2 seed's
coefficient reproduces the sweep's 1.71e-21 from a different run.

The binary64 phase error of the e = 1/2 replay, 1.39e-13 at 391
orbits, sits seven times under the 64-member ensemble's median law
(9.7e-13 at that count): one draw, on the low side, as the 1e7-orbit
run was. The position floor at e = 0.99 peaks at 1.2e-7 at samples
that land near pericentre and reads 7e-9 at the last, which is the
oracle's reason for reporting the along-track phase error instead of
`|dx|/a` (the 2026-09-10 horizon entry).


## 2026-09-10 - the ensemble as an instrument: eight members one to sixty-four ulps apart, at binary64 and at binary256

The chaotic-horizon entry above left one run unfinished: the eight-
member Pythagorean ensemble (`--geometric m3 x 0x1p-52`: member k
has body 3's x displaced by 2^(k-1) binary64 ulps, member 0 is
untouched) on the replayed 6,000-step sequence at binary256. It took
2 h 25 min on the software backend - 72 lanes at binary256 are
element-bound there, 5.22 system-steps a second, 16.93 corrector
passes a step, none rejected, inexact the only flag - and it is the
measurement the ensemble mode was built for. tools/divergence.py,
each member against member 0, `|dx|/L`, and the ratio of each
member's divergence to member 1's:

    binary256                                                    binary64
    t       d(m1)      d(m2)/d(m1) d(m3) d(m4) d(m5)  d(m6)  d(m7)   d(m1)      ratios
    0.0     7.40e-17   2.00  4.00  8.00  16.00  32.00  64.00      7.40e-17   2.00 4.00 8.00 16.00 32.00 64.00
    9.8     1.89e-16   2.00  4.00  8.00  16.00  32.00  64.00      5.76e-13   0.90 1.97 1.51 0.71 2.91 1.30
    19.8    4.05e-16   2.00  4.00  8.00  16.00  32.00  64.00      4.91e-10   1.34 1.90 1.66 1.46 0.63 2.38
    29.8    1.88e-15   2.00  4.00  8.00  16.00  32.00  64.00      2.11e-8    1.34 1.90 1.65 1.46 0.62 2.39
    41.9    2.19e-13   2.00  4.00  8.00  16.00  32.00  64.00      1.22e-7    1.34 1.90 1.65 1.46 0.63 2.38
    52.1    4.22e-11   2.00  4.00  8.00  16.00  32.00  64.00      2.24e-5    1.33 1.90 1.65 1.46 0.63 2.38
    59.8    3.90e-9    2.00  4.00  8.00  16.00  32.00  64.00      2.20e-3    1.32 1.85 1.62 1.44 0.63 2.27
    63.2    1.03e-8    2.00  4.00  8.00  16.00  31.99  63.97      7.80e-4    1.33 1.88 1.64 1.45 0.63 2.34
    73.4    1.13e-9    2.00  4.00  8.00  16.00  32.00  64.00      (replay lost)

At binary256 every member is its seed times one common factor, to
three digits, at every one of 301 samples through the whole chaotic
evolution, every close encounter included - the linear-response
regime of the dynamics, which is what a Lyapunov measurement assumes
and here is seen directly; the common factor climbs from 1 to 2e8
(1.49e-8 / 7.4e-17 at t = 71) along the same staircase as the format
differences of the entry above. At binary64 the same eight members
are within a factor of four of each other in no order from t = 10 on,
member 6 (32 ulps) below member 1 (1 ulp) throughout, because a
perturbation of one to sixty-four ulps is the size of the round-off
noise and the ensemble is measuring the noise. **At binary64 an
ulp-scale ensemble measures the arithmetic; at binary256 it measures
the dynamics.** That is the scientific case for a wide format in a
chaotic system in one table, and it is the ensemble mode's case for
existing: 72 lanes in one run, every member bit for bit its solo run
(the gate), a perturbation 2^184 above the floor instead of at it.


## 2026-09-10 - the ensemble's throughput on the software backend, and the scalar control made independent of E

The same Kepler integration the tile ran (binary256, `--arith fma
--engine program`, fixed dt = 0.05, 20 steps) at E = 1, 4, 16, 64,
256 and 1,024 members two ulps apart, on the software backend, one
job at a time on an otherwise idle box (the binary256 ulp-ensemble
run of the entry above was still on one core during the first scan):

**The first ensemble build.** Program engine, system-steps a second:
18.3 at E = 1, 22.9 at 4, 29.1 at 16, 45.2 at 64, 54.6 at 256, 55.0 at
1,024. Loop engine: 53.5, 58.4, 58.7, 65.0 (E = 64), flat. Binary64,
program engine: 76.6 at E = 1, 259 at 64, 368 at 1,024. The software
backend is element-bound at binary256 from E = 1 - its loop engine
gains nothing from width - and its program engine gains 3x, which is
its per-run overhead amortising; at binary64 the program engine gains
4.8x. Its asymptote, about 60 binary256 system-steps a second on one
core, is the element rate: 84,000 operations a system-step at 5
million operations a second.

The scan also showed something the design document had predicted and
that the card would have paid for: the library-call count grew with
E, from 73,947 calls for the 20 steps at E = 1 to 4,171,889 at
E = 1,024 - 204 calls per system-step, the divides among them from
5,513 to 399,811. A fixed-step run has no step control, so this was
the corrector's convergence test: two scalar compares per lane per
pass for REBOUND's `max` over the coordinates, kept scalar so that the
maximum was the library's, and one divide per system per pass. On the
software backend that was 3 percent of the time. On the tile it would
have been 220,000 round trips of 35 us per step at E = 1,000 - about
eight seconds - against 0.8 s of element work.

**The change, in two steps, each gated at every format.** The maximum
over a system's lanes and the minimum over its particles became
selections by bit pattern on the host (for non-negative IEEE
encodings, neither NaN, the unsigned order is the numeric order; the
element chosen is the one `CFT_CMPLT` would choose, and it is a
selection, not an arithmetic operation), the convergence quotient and
the two exit tests became one E-wide divide and two E-wide compares
per pass, and the step control became a dozen E-wide calls per step:
the particle timescales in one call with the inputs of the particles
REBOUND skips replaced by 1 (their lanes raise nothing and are never
read), both branches of `dt_new` computed for every system and each
keeping its own, the reject and clamp decisions as E-wide predicates.
Then the exact-time update and the prediction ratio, which were three
adds and a divide per accepted system, became E-wide as well, a
system that did not accept adding +0 (exact, never -0 here). Every
element is still the operation the system would issue alone, and the
four gates said so on each build: check_equivalence 1,576 values,
check_program_engine 960, check_records 576, check_ensemble 27 cases
and 11,124 values, all identical, at every format, twice.

**The scans after each step:**

    E       calls, first build   after step 1   after step 2   system-steps/s (fp256 program): first / final
    1       73,947               71,213         71,213         18.3 / 17.7
    4       87,370               72,969         72,669         22.9 / 24.0
    16      135,583              74,351         72,851         29.1 / 33.2
    64      328,130              79,697         73,397         45.2 / 51.3
    256     1,097,051            99,443         73,943         54.6 / 59.9
    1,024   4,171,889            176,971        74,671         55.0 / 62.3
    fp64, E = 1,024              794,128        123,463        21,163         368 / 405 system-steps/s

The call count of a step is now independent of E to within 0.2 calls
a system-step (the remainder is the sample's energy, E-wide already,
and the per-system copies that are not calls), and the tile's
round-trip count per step - about 3,560 for the element work, plus
five per corrector pass, a dozen for the step control and four for
the time and the prediction - is the E = 1 count. The software
backend's rate hardly moved, as it should not have: the calls it
shed were cheap there. What remains E-proportional on the card is the
convergence test's read-back of 3NE deposits per pass, one bulk
transfer, for which the `CFT_MAX` ask of docs/HARDWARE.md stands.

**What the software backend's curve is for.** It is the baseline the
card must beat, measured: 62 binary256 system-steps a second at
E = 1,024 on one core of this host, 405 at binary64. The card's
projection at E = 1,000 (docs/HARDWARE.md) is 1,300 binary256
system-steps a second per tile, engine-bound; against this baseline
that is 21x, not the 100x the projection quoted against the E = 1
software rate, because the software backend gains 3.5x from width
too. That is the honest number to put beside the card's, when it is
measured.
## 2026-09-10 - Parcel A: the drop-in, and the same bits from REBOUND's side

ROADMAP.md's first parcel: `struct reb_integrator` over the engine in
src/ias15_cft.c, so that a REBOUND user changes one line and keeps
their particles, outputs, callbacks and visualisation.

    cft_ias15_register("ias15_cft");
    struct cft_ias15_state* st = reb_simulation_set_integrator(r, "ias15_cft");
    st->format = CFT_FP128;          // optional; the default is CFT_FP64

New: src/cft_ias15.h (the public header and ROADMAP.md's state struct
verbatim), src/reb_integrator_cft.c (the shim), src/ias15_engine.h and a
library mode in src/ias15_cft.c, src/cft_ias15_fields.c (parcel B's
file, today only the list terminator), tools/check_dropin.c (the gate).

### The gate

`build/check_dropin`: the same program twice inside one process, once
with REBOUND's own `"ias15"` and once with the registered `"ias15_cft"`
at binary64, comparing every particle's nine values (x, y, z, vx, vy,
vz, ax, ay, az) and `r->t`, `r->dt`, `r->dt_last_done` as **bit
patterns**. Thirteen cases, 54.2 s on this host:

    kepler, fixed step        N = 2, dt = 0.05, epsilon = 0,      400 steps
    kepler, adaptive          N = 2, dt = 0.05, epsilon = 1e-9,   400 steps
    kepler, tight tolerance   N = 2, dt = 0.05, epsilon = 1e-12,  200 steps
    pythagorean, adaptive     N = 3, dt = 0.01, epsilon = 1e-9,   400 steps
    five bodies, adaptive     N = 5, dt = 0.5,  epsilon = 1e-9,   300 steps
    five bodies, fixed step   N = 5, dt = 0.5,  epsilon = 0,      200 steps
    kepler, rejected steps    N = 2, dt = 8,    epsilon = 1e-9,   200 steps
    pythagorean, rejected     N = 3, dt = 5,    epsilon = 1e-9,   200 steps
    kepler, integrate to tmax        reb_simulation_integrate to t = 37
    five bodies, to tmax             reb_simulation_integrate to t = 211
    kepler, two integrate calls      to 20, then to 53
    kepler, a coordinate edited      200 steps, vy nudged by 1e-6, 200 more
    kepler, a particle added         150 steps of 2 bodies, a third added, 150 more

All identical, every value. `make check-quick` was run after the change
and the four existing gates still say what they said: check_equivalence
6 cases and 528 values, check_program_engine 3 cases and 320 values,
check_records 4 records and 264 values, check_ensemble 9 cases and 1,728
values, all at binary64, all PASS. The engine's arithmetic is untouched,
which is the point: this parcel had to add a shim without moving a bit.

The last four cases exist because each is a place where a shim can be
wrong without the arithmetic being wrong.

- **exact_finish_time.** ROADMAP.md asks the integrator to honour it. It
  does not have to: `reb_check_exit` in simulation.c shrinks `r->dt`
  itself before the last step and restores `last_full_dt` afterwards, so
  the integrator's job is only to use `r->dt` as it finds it and to put
  back the step it actually took. What the shim must do is notice that
  REBOUND changed `r->dt` under it.
- **Two `reb_simulation_integrate` calls.** `reb_simulation_integrate_raw`
  sets `r->dt_last_done = 0` at the top of every call, and step_try reads
  exactly that to decide whether to predict e and b after a rejected
  first attempt. A shim that kept its own `dt_last_done` would diverge
  at the first rejection of the second call - and only there, which is
  why this is a case and not a comment.
- **A coordinate edited between steps.** REBOUND has no reliable flag
  for it: `r->did_modify_particles` is set by `reb_simulation_add` and
  by three integrators, not by a user assigning to `particles[i].vy`.
- **A particle added.** This is an equivalence case, which was not
  obvious. `reb_integrator_ias15_alloc` only reallocates when the array
  grows past its high-water mark, and `realloc_dp7` then zeroes the
  *whole* array, so REBOUND itself discards g, e, b, csb, er, br and
  zeroes csx, csv on a grow - which is exactly what the shim's
  invalidate-and-reset does. **Removal is not equivalent and is not
  claimed:** REBOUND keeps a polynomial that no longer describes the
  particle set (the array only shrinks logically), and the shim
  deliberately resets. Nothing here tests removal.

### The refusals

**Superseded by entry 31** — eight of these are supported now and the
gate asserts them rather than refusing them; the list below is what was
true when it was written.

Thirteen more cases, each asserting that the step names the feature in a
REBOUND error message, sets `REB_STATUS_GENERIC_ERROR` and does not move
the clock: non-zero softening, `additional_forces`,
`force_is_velocity_dependent`, ghost boxes, `REB_GRAVITY_COMPENSATED`,
the tree code, collision detection, periodic boundaries, test particles
(`N_active`), variational particles, MEGNO, `min_dt` and
`adaptive_mode != PRS23`. `r->map` is refused too and is not in the
gate. The assertion is the clock and the status, not the positions,
because REBOUND's own boundary check and collision search run *after*
the integrator callback returns and may touch particles - which is
REBOUND's doing, not this integrator's.

### The wide path, which the binary64 gate cannot reach

At CFT_FP64 promotion and rounding are `memcpy`: the equivalence gate
above exercises no conversion at all. `build/check_dropin --wide` is a
separate run of the same program at binary128 - separate because the
engine's format is fixed for the process once opened - on Kepler at a
fixed dt = 0.05 for 2,000 steps, 16 orbits:

    t: ias15 99.999999999996461, ias15_cft 100 (the wide clock, rounded)
    largest relative difference in the twelve coordinates: 1.619e-13

The clock is the sharpest thing in it. REBOUND adds fl(0.05) to a double
2,000 times and lands on 99.999999999996461; the shim adds it in
binary128, where 2,000 fl(0.05) is exact, and the rounded view is 100.
That is the wide state surviving the step boundary, which is the whole
claim of the division of state, in one number. The 1.6e-13 in the
coordinates is the binary64 run's own accumulated round-off.

**A failure kept beside it.** The first version of that case reported
`0.000e+00` and "bit-identical to binary64 - the format did not take",
and the format had taken: the fault was the gate's own initial
condition. The Kepler pair had been given `v = 2.4506122933227936` at
r = 0.5 against an escape velocity of 2.0009997501249219, so it was a
hyperbolic flyby, and after 16 units of nearly free flight at x = -71
the two runs rounded to the same double. Replaced by a barycentric,
zero-momentum pericentre state, v_rel = 1.7329166165744962, period
6.2800460687587085 - and the difference appeared. Every binary64 case
had passed on the escaping problem too, which is true but was not the
test it said it was.

### What the design looked like when it met the code

Five things, in the order they cost time.

1. **`create()` cannot see `r`.** `void* (*create)()` takes no
   arguments, so ROADMAP.md's "create/free: allocate and release the
   state for `r->N` particles" is not implementable as written.
   Everything sized by the body count is allocated at the first step
   instead, which is the first moment `r->N` is knowable. `free()` is
   symmetric only for the state struct; see 2.
2. **The engine is one global instance, and the scratch is sized once.**
   Every buffer in ias15_cft.c is a file-scope static - correct for a
   standalone program that runs one problem - and several of the step's
   vectors are allocated lazily at their first call and never resized
   (`static V gk; if (!gk) gk = valloc(N3);` and several dozen more like
   it). Rather than touch every one of those declarations in a file
   whose arithmetic is the thing being preserved, ias15_cft.c gained one
   capacity floor:
   `cap_elems`, which `valloc` and the class arrays' `cbytes` apply to
   every allocation, set once from the reserved body count. It is 0 in
   the standalone program, where nothing resizes, so nothing there
   changed. The consequences the user meets are honest and named:
   `cft_ias15_reserve(N)` before the first step if the simulation will
   grow, and one simulation at a time (a second is refused, a freed one
   releases the engine for the next).
3. **`min_dt` is in the state struct and the engine does not have it.**
   The step control issues `fabs(dt_new) < 0` as a comparison against an
   exact zero, because that is what REBOUND's default does and a real
   floor would change the arithmetic. The field is kept, since it is
   ROADMAP.md's struct and parcel B archives it, and a non-zero value is
   refused.
4. **`dt_last_done` is not missing from the state struct.** It looked
   like a gap - the struct has nowhere to keep it and a restart needs it
   - until REBOUND turned out not to keep it in its own IAS15 state
   either: `struct reb_integrator_ias15_state` has no such field and
   step_try reads `r->dt_last_done`, which is a top-level simulation
   field and is already archived. Parcel B should not add one.
5. **`at` is named in ROADMAP.md's prose and absent from its struct, and
   that is right.** `at` is the acceleration at the current corrector
   substep: written and read within one substep, never live across a
   step boundary. REBOUND archives it because its descriptor list is
   mechanical, not because a restart needs it. The shim exposes it
   through `struct ias15_engine_view` for anyone who wants it and the
   state struct does not carry it.

And one thing that is upstream's, found by trying it:
**REBOUND cannot register two custom integrators.**
`reb_integrator_register` scans with
`while (list[N].name){ N++; if (strcmp(list[N].name, name)==0) ... }`,
reading `list[N].name` after the increment, so the second registration
passes the `{0}` terminator's NULL name to `strcmp`. A two-line scratch
program that registers two no-op integrators prints "registering a
second..." and never returns on this host. It does not affect
"ias15_cft" alone, and it is written down in src/cft_ias15.h because
docs/INTEGRATORS.md ranks WHFast next and a "whfast_cft" would meet it.

### What parcel B has, and what it should not do

`struct cft_ias15_state` in src/cft_ias15.h is ROADMAP.md's, field for
field and in its order; nothing has been inserted among the fields and
nothing appended. `n_elem` (3N), `E` (1), `format`, `cft_abi` (from
`cft_abi_version()`, "0.11") and `constants_digest` are filled in and
the blob pointers are refreshed after every step, so an archive written
from a heartbeat sees live pointers. `constants_digest` is FNV-1a over
the exact bytes of the derived h, rr, c and d at the run's format - what
the run used, not what a table says it used.

`cft_ias15_field_descriptor_list` is declared in the header and defined
in src/cft_ias15_fields.c, which is parcel B's file and today holds only
the terminator - which REBOUND reads as "this integrator adds no
fields", so a cft archive written today is a plain, complete binary64
REBOUND archive. REBOUND prefixes the names itself
(`output_fields_from_list`, binarydata.c:646), so a descriptor named
`cft_x0` reaches the file as `integrator.ias15_cft.cft_x0`.

The one thing that needs a decision and cannot be taken here:
`element_size` is a plain field of every REBOUND descriptor and every
built-in list fills it with a compile-time `sizeof`, while W is 8, 16 or
32 according to `state->format`. The writer computes
`size_data = *(size_t*)(base + offset_N) * element_size`
(binarydata.c, `case REB_POINTER`), so there are two shapes that work
and they are not equivalent:

- **three lists, one per format, chosen per simulation.**
  `set_integrator` copies the whole `struct reb_integrator` into
  `r->integrator.callbacks` (simulation.c:205), so
  `r->integrator.callbacks.field_descriptor_list` is a per-simulation
  pointer and may be repointed after
  `reb_simulation_set_integrator`. `n_elem` stays 3N, which is what
  ROADMAP.md says it is, and the archived `size_data` is right.
- **one list, `element_size = 1`, the blob declared as bytes.** Then
  `offset_N` must name a member holding 3N*W and `n_elem` is no longer
  3N. The struct has no spare size_t for it (`E` is the ensemble count),
  so this shape costs a field, and the field cannot be added without
  moving the offsets parcel B is writing against.

The first is the one that fits the struct as it stands. Neither is
implemented here; it is written down rather than guessed at.

## 2026-09-10 - the Simulationarchive: the wide state as `cft_` fields, and what a reader really does with a name it does not know

Parcel B of four, written against ROADMAP.md's shared state struct.
What was added: `src/cft_ias15_state.h` (that struct, verbatim, until
the integrator shim brings its own), `src/cft_archive.c` and its
header - the field descriptor lists, a probe that reads an archive's
own account of itself from the file, and a load that restores exactly,
promotes, or refuses - and four gate programs in `tests/`, run by
`tools/check_archive.py`. `make check-archive` runs them; `make check`
and `make check-quick` now end with them.

**The fields: 59, every one `cft_`-prefixed.** On disk they are
`integrator.ias15_cft.cft_<name>`, because REBOUND builds an
integrator field's name as `"integrator." + r->integrator.name + "." +
descriptor name` and nothing else is possible. The names mirror
REBOUND's own IAS15 descriptor list in `integrator_ias15.c` so the
mapping can be read line by line:

    48 wide blobs, REB_POINTER, element_size = W (8, 16 or 32):
      cft_x0  cft_v0  cft_a0  cft_csx  cft_csv  cft_csa0
      cft_g0..cft_g6    cft_b0..cft_b6    cft_csb0..cft_csb6
      cft_e0..cft_e6    cft_br0..cft_br6  cft_er0..cft_er6
    11 scalars:
      cft_epsilon (REB_DOUBLE)  cft_min_dt (REB_DOUBLE)
      cft_adaptive_mode  cft_format  cft_max_iter  cft_arith_fma (REB_INT)
      cft_E (REB_SIZE_T)  cft_abi_0  cft_abi_1  cft_constants_digest
      (REB_UINT64)  cft_n_elem (REB_SIZE_T, and it must be last)

Three lists exist, identical but for `element_size`, one per format;
`cft_archive_bind()` points `r->integrator.callbacks.field_descriptor_list`
at the one matching the run. A `cft_archive_selftest()` checks their
shape - 48 blobs, 11 scalars, every name prefixed and unique, every
blob's `element_size` the format's width and its `offset_N` the one
member the reader is allowed to overwrite, `cft_n_elem` last - because
a list that has drifted writes wrong bytes silently.

Three of ROADMAP.md's decisions did not survive contact with
`binarydata.c`. They are the useful part of this entry.

### 1. An unknown field does not "warn and seek past it". It stops the snapshot.

The roadmap quotes the right three lines and the wrong control flow.
The case ends:

    case REB_FIELD_NOT_FOUND:
        *warnings |= REB_BINARYDATA_WARNING_FIELD_UNKNOWN;
        int err = fseek(inf, field.size_data, SEEK_CUR);
        ...
        goto finish_fields;

`finish_fields`, not `next_field`: the reader abandons the rest of the
snapshot. Extra fields are safe for a stock reader only because
REBOUND writes every simulation field - particles, `t`, `dt`, `N`, `G`
- before any integrator field, so everything a stock reader needs has
already arrived by the time it meets the first `cft_` one. What it
loses is what follows: the `functionpointers` flag, written after the
integrator's fields. An archive that put a `cft_` field anywhere in
the main list's namespace would truncate a stock read; none does, and
none can, because the prefix is forced.

### 2. A stock reader does not even reach that path. It says the file is corrupt.

Integrator fields are gated on the integrator NAME before the
descriptor lookup happens:

    if (strncmp("integrator.", name, 11)==0 && (name_sub = strchr(name+11,'.'))){
        if (!r->integrator.name[0] || strncmp(name+11, r->integrator.name, name_sub - name -11)){
            *warnings |= REB_BINARYDATA_WARNING_CORRUPTFILE;
            goto finish_fields;

So a reader without `ias15_cft` registered first prints
`Error! Integrator not found.` when it reads `integrator.name` (from
`reb_simulation_set_integrator`, which leaves the reader's own
integrator in place), and then raises CORRUPTFILE - not FIELD_UNKNOWN
- at the first `integrator.ias15_cft.cft_*` field. `gate_stock`,
which links the pinned upstream `librebound.a` and nothing of
cft-rebound's, said exactly that, on a one-snapshot archive and on the
last snapshot of a three-snapshot one:

    Error! Integrator not found.
    gate 2: a cft archive opened by a stock REBOUND reader
      one snapshot: 1 snapshots indexed; warnings = 0x200: WARNING_CORRUPTFILE
        integrator.name in this reader after the load: "ias15"
        binary64 state recovered: bit for bit as written (3 particles, t = 0x1.1eb851eb851ecp-4)
      last of three: 3 snapshots indexed; warnings = 0x200: WARNING_CORRUPTFILE
        integrator.name in this reader after the load: "ias15"
        binary64 state recovered: bit for bit as written (3 particles, t = 0x1.3333333333333p-3)
    gate 2: PASS

The claim the roadmap makes is still true where it matters: the file
opens, the archive indexes (the blob scan seeks by `size_data` and
never consults a descriptor, so unknown fields cost it nothing), and
the binary64 state comes back bit for bit - every particle
coordinate, mass, `t` and `dt` identical to what the writer recorded.
What is not true is the tone. A collaborator without cft is told
`Error! Integrator not found.` and
`The binary file seems to be corrupted. An attempt has been made to
read the uncorrupted parts of it.` for a file that is neither
erroneous nor corrupt. That is structural: any custom integrator name
produces it, and there is no writer-side choice that avoids it short
of patching REBOUND, which this project does not do. It belongs in the
README next to the compatibility claim, not discovered by a
collaborator.

A second consequence, worth stating plainly: what a stock reader
restarts from is the particles and the time, with REBOUND's own IAS15
series cold (`b`, `e` zero), because the archive's `integrator.ias15.*`
fields do not exist. That is a correct restart - it is what REBOUND
does for a new simulation - but it is not a bit-identical
continuation. Only a reader that knows `cft_` gets one.

### 3. `element_size = W` cannot survive a load, and gate 1 caught it

The roadmap says to archive the wide arrays "as `REB_POINTER` with
`element_size` set to the wide width". That is right for writing:
`size_data = pointer_N * element_size` needs `element_size = W` to
produce the correct byte count. It is wrong for reading, and the
reason is that `element_size` is fixed in the descriptor at compile
time while W is a property of the run. On load the list in force is
whichever one the integrator was REGISTERED with - `set_integrator`
installs it when `integrator.name` is read and nothing re-points it
before the blobs arrive - so `binarydata.c`'s

    if (fd.offset_N!=SIZE_MAX){ *pointer_N = (size_t)field.size_data/fd.element_size; }

writes `size_data/8` into `n_elem` 48 times over.

Putting `cft_n_elem` last in the list repairs that for a whole
snapshot: it is written last and read last, so the true count
overwrites the wrong ones. It does NOT repair an appended snapshot,
because REBOUND stores those as a diff against the first
(`reb_binarydata_diff` in the append path of `reb_simulation_save_to_file`)
and the diff omits `cft_n_elem` precisely because it did not change.
The first run of gate 1 said so:

    gate 1: mid-run archive, restart, bit-identical continuation
      fp64      PASS  ... one-snapshot identical, appended-diff identical
      fp128     FAIL  load(-1): refused: element count mismatch
      fp256     FAIL  load(-1): refused: element count mismatch
    gate 1: FAIL

with, from the refusal message,
`cft_n_elem is 9 and cft_E is 1, but the archive holds 48 of 48 wide
blobs and each is 144 bytes, against 3 N particles at 16 bytes an
element` - 9 x 16 = 144, so the file was consistent and the loaded
state was not. binary64 passed only because 72/8 happens to be 9.

The fix is to stop leaning on the mechanism. `cft_archive_probe`
reads `cft_n_elem`, `cft_format`, `cft_E` and the blob lengths from
the file, checks they agree with each other and with 3N (or 3NE), and
`cft_archive_finish_load` then writes the true count into the state.
The blob CONTENTS were never at risk: the read path `fread`s exactly
`size_data` bytes into a buffer it reallocates to that size, so only
the derived count was ever wrong. The lists keep `element_size = W`,
because the write side needs it, and `cft_n_elem` stays last so a
whole-snapshot load is self-consistent even before the repair runs.
The cost is an API rule: an archive loaded with plain
`reb_simulation_create_from_file` and never handed to
`cft_archive_finish_load` holds a state whose `n_elem` is wrong, and
the next save would write blobs of the wrong length. Both entry
points here do the repair.

### Three smaller collisions with the struct

- ROADMAP.md's mirror list names `at`; the struct it defines has no
  `at` member, so there is no descriptor for it. REBOUND archives its
  own `at`. `at` is written before it is read in every substep of
  `reb_integrator_ias15_step_try`, so a restart does not need it; if
  the integrator shim's state grows one, it is one more line.
- REBOUND keeps each seven-level array as ONE allocation of `7*N3`
  doubles and archives it as one field with `element_size = 7*8`
  (`dpcast()` slices it level-major). The roadmap's struct keeps seven
  separate pointers, which cannot be one `REB_POINTER` field, so each
  level gets its own name: `cft_g0..cft_g6` against REBOUND's `g`.
  The bytes on disk are the same in the same order; only the field
  boundaries differ.
- `char cft_abi[16]` cannot be archived by any REBOUND dtype.
  `REB_STRING` and `REB_POINTER` both dereference the field as a
  pointer, and no simple dtype is 16 bytes wide. It goes out as
  `cft_abi_0` and `cft_abi_1`, two `REB_UINT64` at their own offsets;
  the bytes on disk are the string's bytes in order, so it round-trips
  exactly. There is also no member in which to record "this run began
  from a promoted binary64 state" - the load returns it, and nothing
  persists it. A one-word `provenance` member would be worth having.

### The gates

Quick tests only, seconds each. **They do not run IAS15.** Parcel A
owns the integrator shim and it is being written in parallel, so the
gates drive a stand-in (`tests/cft_shim_stub.c`): a deterministic
mixing over the wide state at a real libcft format, arranged so that
every one of the 48 blobs feeds the next step, which is what an
archive gate needs - a blob missing from the list, or written at the
wrong width, breaks the continuation. Gate 1 also asserts that all 48
blobs hold a non-zero byte, so it cannot pass by comparing zeros.
When the shim lands, the gates should be re-pointed at it; the archive
module does not change.

**Gate 1 - write mid-run, restart, continue bit-identically.** Three
particles, `n_elem = 9`, 20 steps against 10 + archive + 10, at every
format, from a fresh single-snapshot file and from the second snapshot
of a two-snapshot one:

    gate 1: mid-run archive, restart, bit-identical continuation
      (superseded by entry 31: 50 blobs + 44 alias blobs + 15 scalars now)
      descriptor lists: 48 blobs + 11 scalars at each of fp64/fp128/fp256, every name cft_-prefixed and unique, cft_n_elem last
      fp64      PASS  3840 state bytes, 48/48 blobs non-zero, 107 cft_ fields, n_elem 9, 72 B a blob, 2 snapshots; one-snapshot identical, appended-diff identical (restored exactly, restored exactly)
      fp128     PASS  7296 state bytes, 48/48 blobs non-zero, 107 cft_ fields, n_elem 9, 144 B a blob, 2 snapshots; one-snapshot identical, appended-diff identical (restored exactly, restored exactly)
      fp256     PASS  14208 state bytes, 48/48 blobs non-zero, 107 cft_ fields, n_elem 9, 288 B a blob, 2 snapshots; one-snapshot identical, appended-diff identical (restored exactly, restored exactly)
    gate 1: PASS

The compared bytes are all 48 blobs, the particles (with the `ap`,
`sim` and `name` pointers zeroed, which are not state), `t`, `dt`,
`n_elem`, `E`, `format`, `max_iter`, `arith_fma` and `adaptive_mode`.
107 `cft_` fields in the two-snapshot file is 59 + 48: the diff
carried every blob and not one scalar, which is the same fact as the
failure above seen from the other side.

**Gate 2 - stock REBOUND opens it.** Quoted above.

**Gate 3 - a stock archive opens here, and is promoted.** REBOUND's
own IAS15, twelve real steps on a three-body problem, saved; then
loaded through `cft_archive_load` at binary128. The promotion is
checked value by value against the same archive loaded plainly, by
rounding the wide bytes back to binary64:

    gate 3: stock archive integrator "ias15", has_cft = 0, cft_ fields = 0
      load said: promoted from binary64 (1)
      promotion of REBOUND's binary64 IAS15 state into fp128: exact, every value round-trips (epsilon 1e-09, adaptive_mode 2, abi "0.11")

That covers `x0` and `v0` (from the particles, which are the authority
at a snapshot - REBOUND's `ias15->x0` holds the start of the last
completed step), `a0`, `csx`, `csv`, `csa0`, and all seven levels of
`g`, `b`, `csb`, `e`, `br`, `er`: 6 x 9 + 6 x 7 x 9 = 432 values, each
identical after the round trip, as widening must be. `epsilon` and
`adaptive_mode` come across from REBOUND's IAS15 state.

**The refusal.** Same gate. A binary256 archive offered to a binary64
run is refused before anything is loaded, with a NULL simulation and,
on stderr:

    cft-rebound: refusing to load ".../gate3_fp256.bin" snapshot 0: it was written at a different precision
      the archive's cft_format is fp256 (32 bytes an element); this run is configured for fp64 (8 bytes an element).
      The wide state is not reinterpretable across this difference and
      will not be silently narrowed or reshaped. Load the archive with
      the settings it was written with, or start a new run.

and the same file at its own format loads `restored exactly`. The
count mismatch has the same shape and its message names `cft_n_elem`,
`cft_E`, the blob count, the blob length, N and the element width, so
the arithmetic that failed is on the screen. `want_format = -1` adopts
whatever the archive says, for a caller who wants the file to decide.

**`make check-quick`**, the whole of it, on this host: PASS, including
`check_equivalence` (6 cases, 528 values), `check_program_engine`
(3 cases, 320), `check_records` (4 cases, 264), `check_ensemble`
(9 cases, 1,728) and the four archive gates.

### Two upstream notes, in passing

Neither is this project's to fix and neither is triggered by anything
here, but both would bite the next parcel.

- **`reb_integrator_register` cannot register a second custom
  integrator.** `rebound.c`'s duplicate-name loop does `N++` and then
  `strcmp(reb_integrator_configurations_custom[N].name, name)`, which
  on the last pass reads the `{0}` terminator's NULL. Exercised: a
  program that registers two printed `registering first...`,
  `first ok`, `registering second...` and never returned; killed by
  PID. (The same loop never compares entry 0, so it would not detect
  a duplicate of the first name either.) Registering exactly one -
  which is all `ias15_cft` needs - is fine, and is what the gates do.
- **`offset_N = SIZE_MAX` is broken by the address fixup.**
  `reb_binarydata_field_descriptor_for_name` adds the base address to
  `offset_N` unconditionally, so the sentinel becomes `base - 1` and
  the read path's `if (fd.offset_N!=SIZE_MAX)` is true, writing eight
  bytes just below the object. REBOUND's own `display_settings`
  descriptor uses that sentinel. Nothing here does, and the archive's
  blobs deliberately do not, which is why the count is repaired from
  the file instead.

## 2026-09-10 - packaging: the body-count cap raised to 1024, the scope stated, a worked round trip, and `make install`

ROADMAP parcel C. Four things, and two bugs that only appeared because
the example was made to run.

### The cap, and what the ROADMAP got wrong about it

`ias15_cft` refused `N > 64`. The ROADMAP said the reason was
`body_names[64][32]`, "the only fixed object". **That is no longer
true of the code**: `body_names` is already `static char (*)[32]`,
`calloc`d from `NB` in `read_problem`, and so is `sys_names`. There
was nothing behind the 64 at all - `valloc` is `calloc`, every array
is sized from N at run time, and the check was a bare literal.

What does constrain N is the shape of the port. Gravity is an explicit
pair list, `P = E*N(N-1)/2`, with 23 vectors of that length; and
`NMAX = max(3NB, P)`, which above N = 7 is P, is also the length of
about 270 broadcast constant vectors (`KRR[28]`, `KC[21]`, `KD[21]`,
`KRINV[28]`, `KHF[7][7]`, and the rest). Both grow as N^2 and both
carry the format's element size. Peak working set, one system,
allocation plus one energy evaluation (`--steps 0`), measured:

    N       binary64    binary128    binary256
    64        8.3 MB       8.5 MB      23.2 MB
    256      79.1 MB     154.2 MB     283.2 MB
    512     283.4 MB     559.9 MB   1,113.4 MB
    1024  1,117.6 MB   2,222.6 MB   4,432.7 MB
    2048  4,454.3 MB   8,872.0 MB  17,707.3 MB

**The cap is now 1024** (`CFT_MAX_BODIES`, src/ias15_cft.c), chosen as
the largest power of two whose worst case - binary256, 4.4 GB - still
fits an ordinary workstation. 2048 would want 17.7 GB at binary256 and
is not a limit anybody could use.

Time reaches further than memory does and is the real constraint. One
fixed binary64 step, the program's own clock, so setup excluded:

    N = 64    1.130 s
    N = 256  13.598 s
    N = 512  51.067 s
    N = 1024 195.546 s

The diff to src/ias15_cft.c is the `#define`, the one changed check
and its message, and nothing else - parcels A and B are editing that
file. `ref/ias15_ref.c` had its own `struct body bodies[64]` and it
grows by `realloc` now, because a program that is diffed against the
port must not impose a smaller limit than the port has.

### The gate at a few hundred bodies

tools/check_bodycount.py, added to `make check` and `make check-quick`.
It reruns the equivalence gate where it could plausibly break and
nowhere else: `gravity()` sums each particle's partners in ascending
partner order to match REBOUND's `(i, j<i)` loop, and a drift in that
order would show as a last-bit difference at large N only. The
problems come from a new tools/make_nbody.py - a star and N-1 circular
test bodies at golden-angle longitudes, barycentric, deterministic
from N alone - so nothing large is committed.

    check_bodycount: cap = 1024, equivalence at N = 256, 1 step(s)
      ok   nbody256 fixed dt=0.005: 2 samples, 3080 values identical
      ok   nbody256 adaptive eps=1e-9: 2 samples, 3080 values identical
      ok   N = 1024 (the cap): allocated and evaluated
      ok   N = 1025 refused: ias15_cft: N = 1025 bodies per system is
           outside 1..1024; memory and time here grow as E*N^2
    RESULT: PASS

and at the size the ROADMAP claimed had been exercised:

    check_bodycount: cap = 1024, equivalence at N = 512, 1 step(s)
      ok   nbody512 fixed dt=0.005: 2 samples, 6152 values identical
      ok   nbody512 adaptive eps=1e-9: 2 samples, 6152 values identical
    RESULT: PASS

The default is N = 256 and two steps, about a minute; `--quick` is one
step; `--n 512` costs about 50 s a step and is therefore not the
default. The cap is read out of the C source by the test, so the two
cannot drift apart.

`make check-quick` passes whole with the new gate in it: equivalence
528 values, program engine 320, records 264, ensemble 9 cases 1,728
values, body count as above.

### The packaging layer

`include/cft_rebound.h` and `src/cft_rebound_run.c`, built into
`build/libcft_rebound.a`. `make install` puts the header, the library,
the `ias15_cft` program and (because `-lcft` has to resolve) libcft
and its headers under `$(PREFIX)`, with `DESTDIR` honoured; PREFIX is
compiled in as the program's location, DESTDIR is not, so a
distribution can build once and stage anywhere. `BINDIR` is a
prerequisite through a stamp file, so `make install PREFIX=X` after a
build for another prefix rebuilds rather than installing a library
pointing at the old one. The two lines a user adds are echoed by the
install itself and marked in examples/Makefile.

The public surface is `cft_rebound_check()` - the refusal list, which
is README "Scope" in code - and `cft_rebound_steps(r, nsteps, opt,
res)`, the analogue of `reb_simulation_steps()`. Today its inside
writes the particles out as an exact binary64 problem file, runs
`ias15_cft`, and converts the record back with libcft's own
`cft_from_hex_char` at `CFT_FP64` - one correctly rounded conversion
per value, which is ROADMAP's "binary64 view" done by the library that
defines the bits rather than by a second parser. The wide state does
not survive a call, so a run is one call rather than a loop of them.
When parcel A's registration lands, that middle is replaced and the
signature does not move.

### The worked round trip

examples/roundtrip.c: Sun, Jupiter and Saturn from orbital elements,
300 fixed steps of dt = 0.5, at binary64 and binary128 on identical
steps. `make example` stage-installs into `build/stage`, builds the
example against that with examples/Makefile, and runs it, so the
example is also the install's test. Verbatim:

      format     steps        |dE/E|        seconds
      fp64         300       5.073e-16       3.23
      fp128        300       0.000e+00      10.40

        fp128 start -0x1.c0343a985a8e5377b710aab7e36ap-14
        fp128 end   -0x1.c0343a985a8e5377b710aab7e368p-14
        fp64  start -0x1.c0343a985a8e4p-14
        fp64  end   -0x1.c0343a985a8e8p-14

      the two runs took the same steps and ended 8.937e-15 AU apart.

binary128's `|dE/E|` prints as exactly zero because `result.energy` is
a `double` and the drift is under its last bit; that is why
`result.energy_hex` carries the run's own unrounded digits as well,
and they show the drift in the last hex digit. binary64's own initial
energy is already half an ulp off the same initial condition's true
value. The example ends by printing five refusals with their messages,
so a reader meets the limits there rather than mid-run.

### Two bugs the example found, which nothing else would have

**`r->OMEGAZ` is -1, not 0.** The first refusal list tested
`r->OMEGA != 0 || r->OMEGAZ != 0` for the shearing sheet, and
`reb_simulation_create()` sets `OMEGAZ = -1.0` as a sentinel meaning
"use OMEGA". Every default simulation was refused. The check is gone
rather than corrected: OMEGA reaches the dynamics only through the
SHEAR boundary or the SEI integrator, and both are already refused.

**Binary256 was silently running a truncated corrector.**
`ias15_cft`'s `max_iter` defaults to 12 - REBOUND's, and correct at
binary64, where it is part of what makes that run REBOUND's own - at
every format. Twenty binary256 steps of a two-body problem through the
packaging layer with that default: `mean_pc = 12.000, max_pc = 12`,
i.e. the cap on every step. The layer now picks the cap by format when
the caller leaves it at 0: binary64 unchanged at 12, binary128 24,
binary256 60, clear of the 2-3 / 4-9 / 10-20 passes measured in this
ledger. The same run then converges at `mean_pc = 14.450, max_pc = 18`
and the energy differs from the truncated one in the last 12 hex
digits. A converged step exits early, so a cap only costs when it is
reached, and `result.iterations_max_exceeded` reports when it was.

### One call, not a loop of calls - measured

The subprocess form rebuilds the wide state per call, and IAS15 takes
the previous step's b coefficients as the next step's predictor, so a
split run is not the same run. Kepler, binary64, fixed dt = 0.05:
twenty steps in one go against ten steps, the recorded state fed back
through tools/state_to_problem.py (exact at binary64), and ten more.
The final states differ:

    one call of 20   ... 0x1.0f08a10952d22p-10  -0x1.0d452883db3a6p-14
                         -0x1.b654664b664cep-2   0x1.06f58990c416ap-4
    two calls of 10  ... 0x1.0f08a10952d21p-10  -0x1.0d452883db3ap-14
                         -0x1.b654664b664cfp-2   0x1.06f58990c4166p-4

one to six ulps, in four of the twelve coordinates. Neither is wrong -
the second is simply a run that threw its predictor away halfway - but
they are not the same answer, so the header, the example and the
README all now say to ask for a whole run in one call. The
registration form does not have the problem: there the state persists.

### The cost, stated in the README because it decides the question

Outer solar system, N = 6, fixed dt = 40. REBOUND's own IAS15, native
doubles: 20,000 steps in 0.22 s, 10.8 us a step, 92,400 steps a
second. The port on the same problem: 40.2 steps a second at binary64,
11.3 at binary128, 3.6 at binary256. **The first factor of 2,300 is
the software library and buys no accuracy at all** - it buys only the
ability to change format. That number was nowhere in the README and it
is the first thing a reader deciding whether to use this needs.

### Not run here

`make check` at every format, the census replays and anything on the
card: this is parcel work and those belong to integration. The
binary256 path through the packaging layer was checked by hand (the
run quoted above) rather than by a committed gate. Nothing touched
XRT, `--artifact` or `cft://`.

---

## 2026-09-10 - Parcel D: is a custom integrator reachable from REBOUND's Python layer? Yes, for one line of loader

The host note at the top of this file does not hold for this entry.
The mechanism under test is POSIX dynamic linking, so it was measured
in the WSL distro (Ubuntu 22.04, gcc 11.4.0, CPython 3.10.12) as well
as on the Windows host (Miniconda CPython 3.12.9, MSYS2 mingw64 gcc).
Two REBOUND builds were used: the pinned clone at bdfda4bd, built here
as a CPython extension module the way `setup.py` does, and - to be sure
the result is not an artefact of a hand build - the **official PyPI
wheels** for rebound 5.1.1, which are at githash 33549d1d. Every claim
below was reproduced on both. No card, no XRT. (The local build reports
`githash b0c25d43`, setup.py's hard-coded fallback: the sources were
copied out of the verified clone at bdfda4bd into the WSL filesystem
without their `.git`, so setup.py's `git rev-parse` had nothing to read.
The bytes are the pinned ones; only the string is the fallback.)

**The question.** `reb_simulation_set_integrator(r, name)` takes a
string, so a Python user might reach a registered C integrator with
nothing but a `ctypes.CDLL`. Parcel A's integrator did not exist while
this ran, so the experiments used a probe integrator written for the
purpose: a drift-only `step` (`x += v*dt`), a two-field state
(`nsteps`, `scale`) with a `field_descriptor_list`, a process-wide step
counter readable from Python, and registration from the library's
constructor. It is not IAS15 and is not meant to be; the load-bearing
question is reachability.

**The Python side does not stand in the way.** `Simulation.__setattr__`
(rebound/simulation.py) lowercases the string, takes five WHFast/SABA
shortcuts, and otherwise passes it to
`clibrebound.reb_simulation_set_integrator` unmodified. There is no
enum and no list to extend.

**The one obstacle is symbol scope, and only on POSIX.**
`rebound/__init__.py` loads librebound with `cdll.LoadLibrary`, i.e.
`RTLD_LOCAL`, so:

    >>> import rebound, ctypes
    >>> ctypes.CDLL("/tmp/parcelD/libstubint_unlinked.so")
    OSError: /tmp/parcelD/libstubint_unlinked.so: undefined symbol: reb_integrator_register

Re-opening the same file with `RTLD_GLOBAL` promotes the mapping
already in the process rather than making a second copy - the address
of `reb_integrator_configurations_custom` read through the promoted
handle and through `rebound.clibrebound` was the same, `0x7f9e7f53a090`
- and the library then loads and registers:

    ctypes.CDLL(rebound.__libpath__, mode=ctypes.RTLD_GLOBAL)
    ctypes.CDLL(".../libstubint_unlinked.so")      # cftstub_registered = 1
    sim.integrator = "stub_cft"                    # -> 'stub_cft'

Linking the library directly against the extension module works too and
needs no promotion, but records an absolute path in `DT_NEEDED` (the
module carries no SONAME), which ties the build to one install.

**On Windows nothing is needed.** The `.pyd` exports
`reb_integrator_register` and `reb_simulation_set_integrator`
(`objdump -p`; the data symbol `reb_integrator_configurations_custom`
is not exported, and does not need to be, because both the write and
the lookup happen inside the DLL). A MinGW DLL that links against
`librebound.cp312-win_amd64.pyd` binds to the module Python has already
loaded, and a plain `ctypes.CDLL` registers. Full round trip on the
Windows host against the official wheel: `sim.integrator -> stub_cft`,
`t = 1.0 y = 1.0` after four drift steps of `dt = 0.25` with `vy = 1`,
`steps_total = 4`, an archive written (4,869 bytes) and read back with
`t = 2.0 integrator = stub_cft nsteps = 8`.

**Everything else in REBOUND's Python layer then works untouched**,
because that layer drives the integrator through the registered
`struct reb_integrator`. Measured: `sim.integrate()` calls `step` (the
process counter reached 4 for `integrate(1.0)` at `dt = 0.25`);
`sim.integrator.nsteps` read 4; `sim.integrator.scale = 2.0` wrote
through to the C state and the next four steps moved the particle twice
as far (y = 1.0 -> 3.0, as predicted); `repr(sim.integrator)` printed
`name=stub_cft, nsteps=4, scale=1.0`; `sim.integrator.__doc__` was
generated from the C documentation string and the field documentation;
`sim.status()` reported `Selected integrator: stub_cft`. An
unregistered name is a clean `RuntimeError: Integrator not found.` and
leaves the simulation on `ias15`.

**Registration from Python alone also works**, and is worth recording
because it removes any need for the library to resolve a REBOUND symbol
at all: `rebound.clibrebound.reb_integrator_register` with
`argtypes = [rebound.integrator.Integrator, c_char_p]` passes the
64-byte struct by value correctly (`ctypes.sizeof(Integrator)` is 64,
as C computes), and `sim.integrator = "stub_py"` then ran the step.

**What the Simulationarchive does.** With the library loaded, a
five-snapshot archive written by `sim.save_to_file(fn, step=2)` was
read back with `rebound.Simulation(fn)` and with
`rebound.Simulationarchive(fn)`: `integrator = stub_cft` in every
snapshot and `nsteps` 0, 2, 4, 6, 8 across them - the custom
integrator's own fields, restored through the Python API. The file
contains the names `integrator.name`, `stub_cft`,
`integrator.stub_cft.nsteps`, `integrator.stub_cft.scale`.

**A failure, and it contradicts the roadmap's summary.** Reading that
same archive in a process where the integrator is *not* registered:

    RuntimeWarning: The binary file seems to be corrupted. An attempt has been made
    to read the uncorrupted parts of it.
    integrator after load: ias15
    t = 2.0 N = 2 dt = 0.25   p1 = 1.0 2.0 1.0
    process_messages(): RuntimeError : Integrator not found.

The binary64 state is recovered correctly, which is what the roadmap
wants, but not by the "warn and seek past it" path. What happens is
that `reb_simulation_set_integrator` fails, `r->integrator.name` stays
`ias15`, and the next field - prefixed `integrator.stub_cft.` - trips
the prefix check in binarydata.c, which sets
`REB_BINARYDATA_WARNING_CORRUPTFILE` and does `goto finish_fields`,
abandoning the rest of the snapshot. The state survives only because
the integrator's fields are written last. Two consequences worth
carrying into parcel B: the graceful `REB_FIELD_NOT_FOUND` seek applies
to simulation-level names, **not** to `integrator.<name>.*` names of an
integrator the reader does not have; and the `Integrator not found.`
error is left queued in `r->messages`, because the file-load path in
simulation.py only decodes the warning bitmask and never drains the
queue - so it is raised by the *next* call that runs
`process_messages()`. In this run that was a later
`sim.integrator = "stub_cft"` which had itself succeeded:

    RuntimeError: Integrator not found.        <- from the load, one statement later

`REB_BINARYDATA_WARNING_CUSTOM_INTEGRATOR` is declared in binarydata.h
and printed for in `reb_binarydata_process_warnings`, but no line in
REBOUND sets that bit. The friendly message the roadmap quotes is dead
code at bdfda4bd.

Its advice - set the integrator after loading - was tried anyway, and
half works. Draining the stale message and then assigning the name
succeeds, but `set_integrator` calls `create()`, so the integrator's own
state is fresh: `nsteps` read 0 where the archive held 8, with
`t = 2.0` and `p1.y = 2.0` intact. Load the library before opening the
archive.

**A second failure, upstream, and it hangs.** The first probe run of
this campaign never returned. `reb_integrator_register` (src/rebound.c)
scans the existing registrations with

    while(reb_integrator_configurations_custom[N].name){
        N++;
        if (strcmp(reb_integrator_configurations_custom[N].name, name)==0){

which increments before testing and so calls `strcmp` on the `{0}`
terminator's NULL `name`. Reduced to a standalone C program against the
same library, the second `reb_integrator_register` in a process never
returns (`timeout 10` -> exit 124), in C as in Python. The loop shape
compiled on its own says why:

    -O0 : [SIGSEGV -- dereferenced the NULL terminator, as written]
    -O2 : LOOP DID NOT TERMINATE (N reached 5000001)
    -O3 : LOOP DID NOT TERMINATE (N reached 5000001)

gcc takes `strcmp`'s `nonnull` attribute to mean the terminator test
can never be false and deletes the loop's exit. REBOUND's wheels are
built `-O3`. So: **one custom registration per process**, and a stale
second copy of the library on `sys.path` is a hang, not an error.
Registering a name that collides with a built-in is by contrast clean -
`Error! Integrator name must be unique but name already exists.`, and a
following `reb_integrator_register(ig, "ias15_cft")` then registered
and selected normally (`integrator is now: ias15_cft`). Loading the
same file twice is a no-op, as `dlopen`/`LoadLibrary` refcount without
re-running the constructor: checked, `sim.integrator -> stub_cft`, no
hang.

**A third, smaller edge.** `__setattr__` lowercases the name, so a
registered name with an upper-case letter is unreachable from Python:
registering `Stub_CFT` and assigning `"Stub_CFT"` gives
`RuntimeError: Integrator not found.` `ias15_cft` is safe.

**Packaging, which turned out to be the easy part.** The official wheel
installs `librebound.<abi>.so` (or `.pyd`) and a `src/` directory
holding REBOUND's headers side by side in site-packages, so a user who
has run `pip install rebound` already has both halves needed to compile
a library against the exact REBOUND that will call it. No source
checkout. `cft_rebound.include_dir()` is
`os.path.dirname(rebound.__libpath__) + "/src"`, verified present on
both platforms. Every function such a library needs -
`reb_integrator_register`, `reb_simulation_set_integrator`,
`reb_simulation_update_acceleration`, `reb_simulation_error`,
`reb_simulation_warning`, `reb_integrators_registered` - is `REB_API`,
so all of them are exported on Windows as well.

**What was shipped, and what it was run against.**
`python/cft_rebound.py` (the loader: `load`, `registered`,
`include_dir`, `library_path`, and a `__main__` that prints all three)
and `python/example_equivalence.py` (the binary64 equivalence gate from
Python: the same data/problems/kepler.txt initial conditions, the same
fixed step, REBOUND's `ias15` against the registered integrator,
printed as exact hex floats). Both were run on Linux and on Windows
against the official wheel. The example's identical path -
`--integrator ias15`, 20 steps of dt = 0.05 - gave
`IDENTICAL: 13 values, bit for bit`, and the same thirteen hex floats
on both hosts (`p1.vx -0x1.08ae6d431ae13p+0`,
`p1.vy 0x1.06f58990c416ap-4`). Its differing path, against the
drift-only probe, gave `DIFFERS: 8 of 13 values, largest absolute
difference 7.074e-01`, which is what a drift-only integrator should
give and is the branch working, not a result about arithmetic.

**Not tested, and stated so.** Parcel A's `ias15_cft` (it does not
exist yet), so the example has never been run against the integrator it
is written for; macOS; anything on the card; and any archive carrying
parcel B's wide `cft_` fields. This repository also has no
shared-library target yet - `src/ias15_cft.c` is a program with a
`main` - so step 2 of docs/PYTHON.md's install order is owed by
packaging, not by Python.


---

## 25. The checkpoint the parcels left broken between them

**2026-09-10, after the four parcels were merged. Windows, mingw64,
software backend, at every format.**

**What was wrong.** Four archive gates passed and the drop-in gate
passed, and a checkpoint taken through the registered integrator still
did not work. Each gate was right about its own half; nothing exercised
the seam. `tests/cft_shim_stub.c` says so in its own header - "when
parcel A lands, the gates should be re-pointed at the real integrator" -
and it never happened, so every archive test drove a stand-in whose step
is not IAS15, and `tools/check_dropin.c` drove the real integrator and
never wrote a file.

`tests/gate_real.c` is that seam, and it is the shape a checkpoint has
to survive: 30 steps straight against 20 steps, save, free, load, 10
more, compared bit for bit in the particles and in `t`. It found three
defects, none visible from either side alone.

**1. Nothing was read back.** REBOUND resolves
`integrator.ias15_cft.cft_x0` against the simulation's descriptor list,
then the built-ins, then the registered custom integrators
(`binarydata.c`, `reb_binarydata_field_descriptor_for_name`). All three
are whatever registration installed, and `src/cft_ias15_fields.c` was
still parcel A's placeholder holding nothing but its terminator, because
parcel B had put its lists in `src/cft_archive.c` as file-scope statics.
So the write side was correct - the file carried 59 `cft_` fields and
all 48 blobs - and the read side skipped every one of them with "Could
not find field descriptor for name". The state kept its `create()`
defaults, `epsilon` came back 1e-9 instead of the archived 0, and the
restarted run took adaptive steps and reached t = 1.63 where the
straight run was at t = 0.30.

Two files each holding half of one fact - the same duplication as the
state struct at merge time, one layer down. The lists moved to the file
whose stated job is to hold them, `src/cft_ias15_fields.c`, and
`src/cft_ias15_fields.h` now carries the macros. The registered list is
the binary64 one at every format, and that is correct rather than
tolerated: REBOUND allocates `field.size_data` bytes from the file and
reads them, using `element_size` only for
`*pointer_N = size_data/element_size`, which
`cft_archive_finish_load()` already overwrites from the file (its
NOTE 3); and 8 divides 16 and 32, so the "Inconsistent size_data"
warning cannot fire for a wide archive read through the narrow list.

**2. The live coordinates were not in the archive.** The 48 blobs were
`x0`, `v0`, `a0`, the three compensated sums and the six seven-element
arrays. But `step_attempt()` sets `x0` FROM `x` at the top of a step and
leaves the advanced position in `x`, so after a completed step `x0` is
the previous step's start and `x` is where the bodies are. At binary64
nothing was lost, because `r->particles` carry the same bits and the
step re-promotes them - which is why this could sit under a passing
binary64 gate indefinitely. At binary128 and binary256 they do not: a
checkpoint silently truncated position and velocity to binary64, the one
thing a wide run exists to avoid. `x` and `v` joined the list; 48 blobs
became 50.

**3. A loaded state was dropped on the first step.**
`publish_state()` points the state's blobs at the engine's buffers, and
REBOUND's loader had just filled those same pointers with the archive's
bytes, so publishing over them discarded the restored state and the run
continued from an engine `reset_state()` had just zeroed.
`adopt_loaded_state()` copies them in, releases them, and lets
`publish_state` repoint as usual. It needs no new field: `publish_state`
only ever sets these to the engine's own buffers and `bind_engine`
refuses a second simulation while the first owns the engine, so a
non-NULL blob pointer that is not the engine's can only have been
loaded. It also writes `t`, `dt` and `dt_last_done` into the engine and
publishes the binary64 view, so the step's "has the user touched these"
guard does not re-promote `r->particles` over the wide `x` and `v`.

**Two more, found while fixing those.** The probe decided what a blob
was from a hand-written list of names, a third copy of `CFT_FD_BLOBS`;
it reported "48 of 50" and every load was refused for an element count
that did add up. It asks the descriptor list now - a blob is exactly a
`REB_POINTER` field named `cft_<tag>`. And `tests/cft_shim_stub.c` kept
a fourth copy of the index-to-member walker, which returned NULL for the
two new blobs and took the three stub gates down with 0xC0000005;
`cft_archive_state_blob()` is exported and the copy is gone.

**The result.** `gate_real` passes at fp64, fp128 and fp256: 30 values
identical across the checkpoint at each, `restored exactly`, 61 `cft_`
fields and 50 of 50 blobs in the file. `make check-quick` is rc = 0 in
5 min 46 s with every other gate unchanged and passing, including
gate 2, which is a stock REBOUND reader opening a cft archive.

**Two more seams closed in the same pass, both about reaching the
card.** Nothing in the gate suite could open an artifact: the drop-in
had `cft_ias15_set_artifact()` and no gate called it,
`struct cft_rebound_options` had no artifact field at all, so a caller
with a U50C in the machine got the software backend and no indication
of it. All three paths now resolve the same way - an explicit setting
first, then `CFT_REBOUND_ARTIFACT`, then software - so
`CFT_REBOUND_ARTIFACT=... make check` runs the whole existing suite on
a card with no gate modified. And the Makefile could not link an
XRT-enabled `libcft.a` at all: every hardware run this repository has
done passed a hand-written `LIBS="-lm -L/opt/xilinx/xrt/lib
-lxrt_coreutil -lstdc++ -lpthread -luuid"` that lived only in a scratch
script, and a plain `make` against such a library failed with a page of
undefined references to `xrt::bo` naming neither cause nor fix. Sourcing
XRT's `setup.sh` sets `XILINX_XRT`; that is the signal now, `XRT=0`
turns it off, and `LIBS=` still overrides.

**Not tested here.** Any of this on the card - the artifact path is
newly reachable and has not been run through a gate yet; macOS; and a
checkpoint written by one process and read by another, which is a
stronger test than `gate_real`'s single process and is what a real long
run does.

*Later: two of those three were closed the same day. The cross-process
checkpoint is `tools/check_checkpoint.py` - entry 26. The card is
entry 27, which is this entry's own gates run on a U50C quad tile with
output identical to software. macOS still stands. Entry 26 also raised
`gate_real`'s step counts from the 20 + 10 recorded above to 60 + 120,
for the margin reason it gives, so the "30 steps straight" shape
described here is not the gate's shape now.*


---

## 26. `make install` shipped the slow path and not the fast one

**2026-09-10, the same integration pass. Windows, mingw64, staged into
a temporary prefix.**

**What was wrong.** `make install` put `cft_rebound.h` and
`libcft_rebound.a` on the system: the subprocess API, where
`cft_rebound_steps()` writes a problem file, runs the `ias15_cft`
program and reads a record back. The registered integrator - what the
README calls the drop-in, what `check_dropin` and every archive gate
exercise, and the only path that keeps the engine in the caller's own
process - was built as `build/libcftrebound.a`, referenced by nothing
but a `dropin` convenience target, installed nowhere and named in no
document. Its headers were not installable either. A user following the
README could install the slow path and had no way to install the fast
one.

The two names also differed by one underscore, which in a single `lib`
directory is a trap rather than a distinction.

**What was done.** The drop-in library is `libcft_ias15.a`, named for
the integrator it registers. It carries `cft_archive.o` as well, so the
installed library can take a checkpoint - previously the archive was
compiled per-gate and existed nowhere a user could link it. `install`
ships it with `cft_ias15.h`, `cft_ias15_state.h` and `cft_archive.h`;
`uninstall` removes them.

**How it was checked, which is the point.** `examples/dropin.c` is a
program written from outside: it includes only headers from the install
prefix, links only `-lcft_ias15 -lcft` from it plus the caller's own
REBOUND, and mentions no path inside the source tree. It registers the
integrator, runs 40 fixed steps at binary256, binds the archive
descriptors, saves, frees the simulation, reloads and continues. Run
through the repository's own `make example`:

```
40 binary256 steps: t = 0.40000000000000002
  p1.x = 0x1.d7975ec563c7fp-1
reloaded (restored exactly): t = 0.40000000000000002,
  p1.x = 0x1.d7975ec563c7fp-1  -> the same bits it was saved with
ten more: t = 0.5
```

`examples/roundtrip.c` still covers the subprocess API, so `make
example` now builds and runs both ways in, and a change that breaks
either is a build failure rather than something a user finds.

**One thing measured and not fixed.** `make install PREFIX=` with a
Windows path written in backslashes fails to compile:
`CFT_REBOUND_IAS15_DEFAULT` is a `-D` string literal and `C:\Users\...`
becomes `\U`, an incomplete universal character name. Forward slashes
work, which is what the README's Windows section already tells a reader
to use for `TMP` and `TEMP`. Recorded rather than fixed because the fix
is a quoting change in a path that reaches a preprocessor definition
and it deserves its own test, not a fix bundled into an integration
pass.

**Also recorded here: the cross-process checkpoint.** Entry 25's "not
tested" list named a checkpoint written by one process and read by
another. It is `tools/check_checkpoint.py` now - three separate
processes per format, exact hex floats compared by the checker - and it
passes at fp64, fp128 and fp256. It carries its own control: a fourth
process runs the resume with `CFT_REBOUND_NO_ADOPT` set, which discards
the loaded state instead of carrying it, and the checker requires that
run to DIFFER. Without the control the comparison would not be
evidence: at 60 + 120 steps of this problem only 2 of the 14 dumped
lines differ between fp64 and fp256, so agreement on the other 12 says
nothing about whether the wide state was carried at all.

**Still not tested.** macOS. Any of this on a card at the time of
writing - the artifact path is newly reachable and a run is in progress.
The Python shared library, whose requirements are now specified in
docs/PYTHON.md rather than described: two of the three unknowns are
answered there, and the third is that neither host has the `rebound`
wheel installed, so the gate that would prove it cannot run as things
stand.

**A long checkpoint, measured.** The gate uses 60 steps and 120 because
it has to be affordable on a card. The question a long run actually
asks is whether the state still restores exactly after thousands of
steps, where the compensated sums and the seven polynomial arrays have
had time to accumulate. With `--steps-a 2000 --steps-b 4000` at
binary256, three separate processes:

```
2000 steps, checkpoint, 4000 more, against 6000 straight, binary256
  differing lines: 0 of 14
  t = 0x1.ep+5   (exactly 60.0)
  elapsed: 1747 s for all three runs
```

Bit for bit, so nothing about the restore degrades with run length -
which is the property that makes the archive useful for the workloads
this card exists for, and it had not been measured.

*Later: this entry's "still not tested" list is closed except macOS.
The card run in progress at the time of writing is entry 27. The Python
shared library, blocked here on neither host having the `rebound` wheel,
is entry 28 - built, gated, and bit for bit from Python at every format.
Entry 28 also corrects what this entry pointed a reader at:
docs/PYTHON.md's claim that `libcft.a` links into a shared object as it
is was wrong, and `backend_xrt.o` is why.*


---

## 27. The first hardware run of any gate in this repository

**2026-09-10, amd-arc-box, Alveo U50C, `cft_hw_quad.xclbin` (four
tiles), XRT 2022.2. Reached with `CFT_REBOUND_ARTIFACT`, which did not
exist before entry 25.**

Every hardware number in docs/HARDWARE.md until today came from the
standalone `ias15_cft` program, because nothing in the gate suite could
open an artifact: the drop-in had a setter no gate called, the
subprocess API had no field for one, and the Makefile could not even
link against an XRT build of libcft. All three are fixed, so these are
the gates themselves, unmodified, run on the tile.

| gate | what it establishes | card |
|---|---|---|
| `check_dropin` | the registered integrator is REBOUND's own IAS15 at binary64, bit for bit - every case and every refusal | 2,030 s |
| `check_dropin --wide` | the same shim at binary128 | 2,195 s |
| `gate_real --fp64` | a Simulationarchive checkpoint written and restored on the tile, 30 values identical across it | 27 s |
| `gate_real --fp128` | the same at binary128 | 48 s |
| `gate_real --fp256` | the same at binary256 | 157 s |

And the comparison that matters most, run last: the same gate in
software and on the card, output diffed.

```
fp64:  software and card report identical output
fp128: software and card report identical output
fp256: software and card report identical output
```

**What is and is not new here.** The cost is not new. 2,030 s against a
few seconds in software is the per-call cost at N = 2 and N = 5 that
docs/HARDWARE.md predicted from the start, and the ensembles are where
the card earns its keep. What is new is that the *answer* does not
depend on the backend, measured through the gates rather than through
one program, and that a checkpoint round-trips through the tile - the
archive was written by a run on the card, read back, and the
continuation was identical to an uninterrupted card run.

**One number worth reading.** `gate_real` on the same 30 steps took
27 s, 48 s and 157 s at binary64, binary128 and binary256. The ratio
157/27 is 5.8. The corrector's measured pass counts are about 2.4, 6
and 18-22, a ratio of 9.2. So the tile recovers roughly a third of what
the extra passes cost, which is what a format-independent per-call cost
and format-dependent per-pass work predict together. It is a
consistency check on the model, not a new result.

**Two things the box itself revealed.**

`make check` had **never once run on this Linux host**. It dies two
seconds in: `tools/gen_constants.py` imports mpmath and the system
python has never had it, because every hardware session here built
`build/ias15_cft` directly and none ran the suite. Installed into
`~/.local` (PEP 668 refuses a plain `--user`; `python3.12-venv` is not
installed either, so `--break-system-packages --user` was the option
that needed no sudo). The full software suite on Linux is therefore
still owed and is the one leg of this pass that did not run.

And the link flags. Every hardware run this project has ever done
passed `LIBS="-lm -L/opt/xilinx/xrt/lib -lxrt_coreutil -lstdc++
-lpthread -luuid"` on the make command line, and that string existed
only in a scratch script on the Windows host. A plain `make` against an
XRT libcft stopped on a page of undefined references to `xrt::bo`
naming neither the cause nor the fix. The Makefile derives them from
`XILINX_XRT` now, and the log of this run shows them appearing in every
link line without a script's help.

**Still not run.** The full software `make check` on Linux (above); the
cross-process `tools/check_checkpoint.py` on the card, which is a
stronger test than `gate_real` and is what a long run actually needs;
macOS; the Python shared library, whose requirements are specified in
docs/PYTHON.md and whose gate cannot run until a `rebound` wheel is
installed on one of these hosts.

*Later: the last of those is closed - entry 28, on this same host, with
the rebound 5.1.1 wheel installed. The other three stand: the full
software `make check` on Linux, `check_checkpoint.py` on the card, and
macOS.*


---

## 28. Python, closed

**2026-09-10, amd-arc-box, Ubuntu 24.04, gcc 13, the rebound 5.1.1
wheel from PyPI.**

docs/PYTHON.md's step 2 had been owing since parcel D: the loader and
the example existed and there was no library for them to load. It is
built and it passes.

```
loaded build/libcft_ias15.so; registered integrators:
  ias15 whfast sei leapfrog janus mercurius saba eos bs whfast512 trace none ias15_cft
problem kepler.txt: N = 2, G = 0x1.0000000000000p+0, dt = 0x1.999999999999ap-5, 252 steps
  configured through the library: fp64, fixed step
IDENTICAL: 13 values, bit for bit.
```

REBOUND's own `ias15` against `ias15_cft` selected by name from Python,
same problem, same fixed step, at binary64.

**That the format reaches the arithmetic, and is not merely accepted.**
The same 252 steps:

| format | wall |
|---|---|
| fp64 | 5.2 s |
| fp128 | 15.3 s |
| fp256 | 49.4 s |

which is the software backend's known ratio (120, 39 and about 10 steps
a second). At fp128 and fp256 the run differs from REBOUND's by
4.263e-14, and the two wide formats agree with each other, which is what
the Phase 2 precision entry says should happen at a step this size.

**Four things it needed. One of them was a correction to this file.**

Earlier today docs/PYTHON.md recorded that `libcft.a` links into a
shared object as it is, on the evidence that its `device.o` carries no
absolute relocations, and entry 26 pointed a reader at that section as
one of the answered unknowns. **That was wrong, and wrong in a way
worth naming: I sampled one member of an archive and generalised.**
`backend_xrt.o` is C++ and the link stops on `R_X86_64_PC32 against
symbol _ZSt7nothrow@@GLIBCXX_3.4`. The fix needed no new code -
cft-fp256 already builds a real `libcft.so` from its own PIC objects,
and with `XRT=1` it carries the card backend too (329 KB without,
393 KB with).

The other three: PIC objects of the four drop-in sources; REBOUND's
symbols left undefined and resolved because `cft_rebound.load()`
promotes the wheel's librebound to `RTLD_GLOBAL` first, with
registration in a constructor because a Python caller's first contact
is an assignment to `sim.integrator`; and a way to configure it.

**The last of those did not exist and is the interesting one.**
Selecting the integrator by name always worked. Configuring it never
did: REBOUND generates `sim.integrator.epsilon` and its siblings from
structs it knows, and it does not know `struct cft_ias15_state`. So a
Python caller could select `ias15_cft` and then reach binary64 with the
default step control and nothing else - the wide formats, the whole
point, were unreachable from Python and nothing said so.

The example had been showing it all along and it read as an arithmetic
failure: it sets `epsilon = 0` for a fixed reference, could not for
ours, printed a note, and compared a fixed-step run against an adaptive
one. 9 of 13 values differed, t was 31.2 against 12.6, and the
arithmetic was never involved. `cft_ias15_configure()` takes format,
epsilon and max_iter by value - which is what a ctypes caller can call -
and `cft_ias15_format_code()` spells the format so Python need not
hard-code an enumerator belonging to libcft.

**One guard, because the agreement is luck.** The library is compiled
against the pinned REBOUND headers and loaded into a process running the
wheel's REBOUND. If those differ, `struct reb_simulation` is laid out
differently on each side of the call and nothing says so. Today the
5.1.1 wheel ships a `rebound.h` byte for byte identical to the pinned
bdfda4bd; `make python-lib` runs `check-rebound-match` first and refuses
rather than build a library that would be quietly wrong.

**Two more things this pass fixed, both found by trying to compile.**
`CFT_IAS15_INTEGRATOR_NAME` was defined in `src/cft_archive.h` - the
integrator's own name, in the archive's header, because parcel B needed
it and parcel A's header did not have it. And `src/cft_ias15.h` declared
an array of `struct reb_binarydata_field_descriptor` without including
the header that defines it, so it compiled only after whatever had
included `binarydata.h` first. A caller following the header's own usage
comment got an incomplete element type. Both are the same shape as
entry 25: a header that is not self-contained is a fact stated
somewhere else.

**Not done: Windows.** A DLL may not carry undefined symbols, so it must
link against `rebound.__libpath__` rather than leave REBOUND's symbols to
the loader. Different rule, its own change, its own test.


---

## 29. The Linux legs, the card checkpoint, and what a documentation audit found

**2026-09-10, amd-arc-box (Ubuntu 24.04, gcc 13, Alveo U50C) unless
stated.**

### The full software suite on Linux, for the first time

```
software make check: rc=0 in 4454 s
  RESULT: PASS lines: 6
```

This box has had a working card for weeks and had never run the gate
suite once, because `tools/gen_constants.py` imports mpmath and the
system python had never had it — every session here built
`build/ias15_cft` directly. So this is the first confirmation that the
whole suite passes on a second platform and compiler at every format.
It takes 74 minutes there against 48 on the Windows host.

### The cross-process checkpoint, on the card

`tools/check_checkpoint.py` with `CFT_REBOUND_ARTIFACT` pointing at the
quad tile: three separate processes per format, each opening the card
for itself, plus a fourth that discards the loaded state and must
differ.

```
card checkpoint: rc=0 in 1766 s
  fp64   PASS  14 values identical, control differs, archive 10536 B
  fp128  PASS  14 values identical, control differs, archive 15336 B
  fp256  PASS  14 values identical, control differs, archive 24936 B
```

This is the strongest form of the claim: a run on the tile wrote a
Simulationarchive, exited, and a different process opened the tile
again and continued from it bit for bit. Entry 27's `gate_real` runs
were the single-process version of the same test.

### `make example` failed, and the reason was a real user-facing gap

It failed instantly, with a page of undefined references to `std::`
symbols out of `backend_xrt.o`. `examples/Makefile` is a standalone
makefile whose whole purpose is to MODEL a user's own build, and it had
no XRT link flags — the top-level Makefile learned to derive them from
`XILINX_XRT` earlier the same day and the example did not. So the
README's "two lines go into your own build" was incomplete for anybody
whose libcft carries the card backend, which is everybody who has a
card.

It passes on Windows and on any software-only build, which is how it
survived a whole day of testing. It took a build on a machine that
actually has an XRT libcft installed. Fixed, and `make example` on the
box now runs both examples:

```
cft-rebound worked round trip: Sun, Jupiter, Saturn; 300 fixed steps of dt = 0.5
40 binary256 steps: t = 0.40000000000000002
reloaded (restored exactly): p1.x = 0x1.d7975ec563c7fp-1 -> the same bits it was saved with
```

`p1.x` there is bit-identical to the Windows run recorded in entry 26.

### An Opus agent audited every document against the code

It found eleven things wrong in the sense that sends a reader somewhere
that does not exist: a `cft_rebound_register()` reserved in
`include/cft_rebound.h` and never written; a library named
`libias15_cft.so` in four places where the Makefile builds
`libcft_ias15.so`; `make examples` where the target is `make example`
(the plural matches the *directory* and silently does nothing); and —
the one worth naming — `src/cft_ias15.h` still calling `x0/v0/a0`
"position, velocity, acceleration", which is the exact comment
ROADMAP identifies as the root cause of the broken checkpoint, still
sitting there after the fix.

`docs/VALIDATION.md` was left alone as it must be: 21 insertions, 0
deletions, all of them pointer lines on entries 25 to 27 naming the
later entry that closed each "not tested" item.

### Three code defects it found, and one of them needed measuring

**1. The drop-in ran every format at REBOUND's cap of 12.** ROADMAP
says the shim must pick per-format caps; only
`cft_ias15_configure(..., 0)` ever did, so every caller that set
`state->format` directly — `check_dropin`, `gate_real`,
`examples/dropin.c` and the README's own snippet — got 12.
`state->max_iter = 0` now means the default for the format, resolved on
the first step and written back so an archive carries the number
actually used.

The audit's note said the binary256 run "reaches the cap on every step
and the iteration is truncated". **Measured, that is not what happens
on this problem.** Asking the engine directly, 180 steps of the
four-body problem in `examples/dropin.c`:

| cap | step attempts that hit it | library calls | p1.x |
|---|---|---|---|
| 60 (the new default) | 0 | 1,284,801 | `-0x1.d0e9da952b141p-3` |
| 12 (the old one) | 2 | 1,282,743 | `-0x1.d0e9da952b141p-3` |

Two attempts out of 180, and the binary64 output is bit-identical. The
fix is still right — an intent stated in ROADMAP and not implemented is
a defect — but nothing recorded in entry 26 is invalidated, and the
claim of a truncated iteration would have been wrong to publish.

**2. `cft_archive_state_blob()` hardcoded its last two indices** as 48
and 49 rather than deriving them from `CFT_N_BLOBS`. They derive now,
and `cft_archive_selftest()` walks every index requiring each to yield
a distinct non-NULL member and index `CFT_N_BLOBS` to yield none. The
guard is the part that matters: a stale copy of this walker is what
took three gates down with an access violation earlier the same day.

**3. Superseded by entry 31 — the REBOUNDx force is computed now, and
`r->gravity_custom` is refused on both paths.** As it stood:

**3. The drop-in's refusal list is a strict subset** of the subprocess
API's in three places: `pre_`/`post_timestep_modifications`,
`r->gravity_custom` and `r->N_odes`. A REBOUNDx force or an attached ODE
set reaches the drop-in unrefused. **Not fixed** — the timestep hooks
look deliberate (`reb_integrator_cft.c` relies on them) and the other
two need a decision rather than a patch. Recorded so it is not
discovered by someone's wrong answer.

### An open question, closed by running it

The audit corrected my explanation of why `cft_rebound.configure()`
exists. I had written that REBOUND generates a built-in integrator's
settings from structs it knows and does not know `struct
cft_ias15_state`. The audit said that is the 3.x/4.x story and that
5.1.1 resolves `sim.integrator.<field>` against the **registered**
descriptor list, custom integrators included — and honestly marked it
unverified, having no wheel to test on.

Run on the box:

```
type(sim.integrator)       : IntegratorConfiguration
  get epsilon      -> AttributeError: Field 'epsilon' not found in BinarydataFieldDescriptor
  get cft_epsilon  -> 1e-09
  get cft_format   -> 1
  get cft_max_iter -> 12

attribute path vs configure(): IDENTICAL
```

Setting `sim.integrator.cft_format`, `cft_epsilon` and `cft_max_iter`
directly produces a run identical to `cft_rebound.configure()`. So the
audit was right and I was wrong: the fields ARE reachable from Python.
What is not reachable is the spelling every REBOUND example uses,
because our descriptor names carry a `cft_` prefix. `configure()`
remains worth having for the validation and the per-format default, but
it is a convenience, not the only way in, and the documentation now
says so.

The same probe also confirmed `cft_max_iter` reading back as 12, which
is defect 1 seen from Python.

### Still not run

macOS. Windows Python (a DLL may not carry undefined symbols, so it has
to link `rebound.__libpath__`, which ties the build to a Python minor
version — docs/PYTHON.md has the detail). CI, which this repository does
not have. And four documents — `docs/ENSEMBLE.md`, `docs/HORIZON.md`,
`results/horizon/README.md` and `NOTICE` — got only a targeted grep in
the audit and were not checked line by line; the audit says so itself
rather than implying coverage it did not have.

**The audit's own coverage gap, closed.** It said plainly that
`docs/ENSEMBLE.md`, `docs/HORIZON.md`, `results/horizon/README.md` and
`NOTICE` got a targeted grep rather than a line-by-line check, and that
it could not verify their internals. Checked afterwards: every path,
tool, flag and column name in all four resolves - `--dt-file`,
`--member`, `--fit-from`, `pc_lane_efficiency`, the `rel_phase_err`
column `tools/oracle.py` writes and `tools/horizon.py` reads, the
`predict-ens` programs at all three formats, and the ten entries
`results/horizon/README.md` describes against the ten that are there.
`SETACT`, which is not in `cft.h`, is a real sequencer instruction
upstream and the present tense is right. NOTICE attributes REBOUND
under GPL-3.0 and libcft under Apache-2.0 with the one-way combination
rule stated correctly. Nothing to fix in any of the four.

## 30. Collisions, and a defect that was not the one we wrote down

Collision detection and resolution are supported on the drop-in path.
The route there went through a wrong belief that this project had
written into three comments and a document, and the correction is the
substance of this entry.

### What was believed

REBOUND's IAS15 was thought to leave a *stale polynomial* behind when a
particle is removed: the arrays keep describing the old particle set,
so each survivor inherits its neighbour's `b`, `e`, `g`, `br`, `er`,
`csb`. Reproducing that looked like a choice between being faithful and
being right, which is how it was put to the decision and how the
decision was made: match REBOUND, defect and all.

### What is actually there

All seven coefficient levels share **one flat allocation**, and

    static struct reb_dpconst7 dpcast(double* dp, size_t N3){
        .p0 = dp, .p1 = dp+1*N3, ... .p6 = dp+6*N3

`integrator_ias15.c:217`, called at 307-312 with the **current** `N3`,
every step. `reb_integrator_ias15_alloc()` reallocates only when
`N3 > ias15->N_allocated`, so a removal - or a regrow under the mark -
leaves the buffer untouched and re-slices it at a shorter stride. Level
*m* is then read starting `m*(old_N3 - new_N3)` elements early: level 1
begins in the tail of level 0, level 2 in the tail of level 1, and so
on up.

That is a falsifiable statement, unlike "stale", so it was tested. A
probe compared REBOUND's arrays against the port's, member by member,
on the step a merge took three bodies to two (`old_N3 = 9`,
`new_N3 = 6`):

    b1 [0] reb 0x1.81c06a463eabfp+5    cft -0x1.44481639f6f84p-7
    b1 [3] reb -0x1.44481639f6f84p-7   cft -0x1.1ab9834c788ffp-13

REBOUND's `b1[3]` **is** the port's `b1[0]`, bit for bit - a
displacement of exactly 3, which is `old_N3 - new_N3`. Level 0 does not
appear in the diff at all, because offset 0 does not depend on the
stride. Both halves of the prediction, on the nose.

So the polynomial is not stale. It is aliased, and the survivors do not
inherit each other's coefficients - they inherit fragments of the level
above.

### What was done about it

Performed, rather than approximated. `ias15_engine_alias_resize()`
gathers the seven levels into a flat shadow at the old stride and reads
them back at the new one. Through a shadow for two reasons, both
measured rather than reasoned about afterwards:

- **In place is wrong.** At `9 -> 6`, level 2's element 0 is flat index
  12, which is level 1's element 3 - already overwritten by the time
  level 2 is reached, in any ascending order.
- **The tail has to persist.** A shrink strands the part of the buffer
  between `7*new_N3` and `7*old_N3`: no stride reaches it while the
  count is small, and a regrow below the high-water mark reads it
  straight back. Holding the shadow between resizes is what makes that
  exact.

Neither claim is left as reasoning. Both were run as negative
controls, against binaries built for the purpose and then deleted:

- **the re-slice itself**, made a no-op, which is what this port did
  before: `case_collision` and `case_remove_then_add` both fail, the
  merge case in `x`, `y`, `vx` and `vy` of the surviving pair.
- **the shadow, zeroed on every resize** rather than kept between
  them: `case_remove_then_add` fails - `x`, `y`, `z`, `vx`, `t` and
  `dt` - and `case_collision` still passes, because a merge is a pure
  shrink and never reads the stranded tail. So the two cases are not
  redundant, and the second one is the only thing holding the shadow's
  persistence down.

### The results

- `case_collision`: a star and two small bodies that merge, REBOUND's
  own `reb_collision_resolve_merge`, 2000 steps. Bit-identical - two
  particles by nine values, plus the clock - with the merge at step
  1291 and 709 steps after it for a divergence to show in. A
  step-by-step probe on the older setup put the first divergence on the
  very next step after the removal, so 709 is generous rather than
  necessary.
- `case_remove_then_add`: bit-identical, and the negative control above
  says it is testing what it claims.
- `case_wide_collision`: binary128 with `state->accurate = 1` across the
  same merge. Survives, finite, both sides remove the same particle,
  largest relative difference from binary64 afterwards `3.403e-10`, and
  **not** bit-identical - which is the assertion, since matching would
  mean the wide tails had been lost.
- `case_wide_collision_refused`: binary128 with `accurate = 0` refuses
  and the clock does not move.

### The flag, and why the rule is what it is

`state->accurate` is required above binary64 and is refused there
rather than silently approximated. A removal shifts `r->particles`,
which are binary64; at `accurate = 0` - where this port reproduces
REBOUND - the wide state is re-promoted from them and every coordinate
loses its tail. At `accurate = 1` the wide state shifts with the
particles (`ias15_engine_remove_body()`), each survivor keeps its own
polynomial, and the run is deliberately not REBOUND.

At binary64 the flag changes nothing that matters, because the view is
the state.

### One test that was measuring the wrong thing

The merge case first gave the two small bodies `+0.3` and `-0.3` in
`vy`. The remnant of a merge between them carries almost no angular
momentum, so it fell into the star: `3 -> 2` at step 64, `2 -> 1` at
step 94. A one-body IAS15 run has no pairs, so the step controller
divides 0 by 0, and the two sides then disagreed on the **sign of the
quiet NaN** - `fff8000000000000` against `7ff8000000000000`. 754 leaves
the sign of a NaN unspecified and x86's default QNaN is negative where
libcft's is positive, so that difference is real, permitted, and says
nothing whatever about the removal the case exists to test. Straddling
the circular speed instead (`1.3` and `0.7`) keeps the remnant in
orbit. The case now also fails outright if the run collapses to one
particle, so it cannot quietly go back to comparing two NaNs.

Late also matters: bodies that overlap at `t = 0` merge on step 1,
where the polynomial being re-read is still all zeros and re-reading it
proves nothing. Two of the five candidate setups did exactly that.

### What is still open

Across a checkpoint. REBOUND archives its coefficient arrays as one
`REB_POINTER` of `N_allocated` elements of `7*sizeof(double)`, so both
the stranded tail and `N_allocated` survive a save; this port archives
each level as its own blob of `3N` and has no `N_allocated`. A run that
removes, checkpoints, resumes and then regrows below the *original*
mark therefore diverges. Closing it changes the on-disk field list, and
the sequence is narrow enough that it is recorded in ROADMAP.md rather
than fixed on the way past.

## 31. Round 2: eight capabilities, five gates that could not fail

Five parcels in parallel worktrees, two verifiers and a follow-up
parcel, against `docs/PARCEL-ROUNDS.md`'s method. What landed matters
less than what nearly did not get caught, so the gates come first.

### Five gates that could not fail

Each found by the parcel that wrote the code the gate was meant to hold
down, except the last, which a verifier found.

**1. A control that left every case passing.** Parcel P3 disabled the
INDIVIDUAL branch it had just written in `pc_error()` and **every case
still passed**. The corrector's pass count does move - 400 steps of
kepler at `dt = 0.01` give a mean of 2.595 against 2.570 and a peak of 5
against 4 - but on a converged step an extra pass changes no bits. It
swept for a configuration that reaches the state and found one: five
bodies from `dt = 500` at `epsilon = 1e-9`, where the first attempt is
rejected hard. With that case added the control gives **1 failure, 48
values**, on INDIVIDUAL, while GLOBAL correctly still passes - GLOBAL
shares PRS23's corrector error.

**2. A case that triggered before the thing it tested existed.** The
merge case first fired on **step 1**, where the polynomial being
re-read is still all zeros. Re-read zeros, got zeros, passed. It also
collapsed to one particle, where a one-body IAS15 run has no pairs,
divides 0 by 0 in the step controller, and the two sides differ only in
the **sign of the quiet NaN** - `fff8000000000000` against
`7ff8000000000000`, which 754 leaves unspecified and which says nothing
about a removal.

**3. A control that was vacuous at two of three formats.**
`CFT_REBOUND_NO_ALIAS` **passed** at binary128 and binary256: the
particles are the binary64 *view*, and a corrector starting from
different coefficients re-converges to within 1e-34 and rounds to the
same doubles. Every dump now ends with a digest of all 50 wide blobs.
The verifier confirmed the fix bites - 1 of 15 lines differs at both
wide formats, and that line is the digest.

**4. A configuration where the physics made the difference exactly
zero.** With strictly massless test particles, setting `N_active`
changes **no bits**: every pair the skip removes contributes
`prefact * 0 = ±0`, the sum starts at `+0` (a memset) and can never be
`−0`, and the real contributions keep their order because the massive
particles are the low indices. **Star, planets, massless dust - the
obvious case - passes with nothing implemented.** Every divergence
control now runs on *massive* semi-active bodies, and the inertness is
recorded as a named case.

**5. A bit comparison that could become a value comparison.** The
`memcmp` deciding whether a force-touched component keeps its wide
value can be swapped for `!=` and **pass the entire suite** while
breaking signed zero (`az[1]`: REBOUND `-0` = `8000000000000000`, port
`+0` = `0000000000000000`). Found by a verifier; the case that catches
it counts the **40,128** zeros it met so it cannot pass vacuously.

**And a sixth, in shared scaffolding.** With the force-write detection
kept but the step-discard removed, all three `refused()` cases **still
pass**: the refusal fires at the n = 0 call, where the hook has just set
`r->t = t_beginning` and `r->dt = dt`, so `refused()`'s "the clock did
not move" assertion is satisfied while the particles sit at
Gauss-Radau node values. `refused()` is used by every refusal case in
this repository.

### What landed

Eight capabilities came off the refusal list, 21 rows to 14 on the
drop-in.

- **`r->additional_forces` at every Gauss-Radau node.** 22 cases, all
  bit-identical, including the whole **call stream** - `r->t`, `r->dt`,
  `x`, `v` and the gravity handed over, per call, in order: **1041
  calls, 5 values each**, and 929 on pythagorean.
- **Velocity-dependent forces**, which were not free: the port had no
  velocity predictor at all. REBOUND's is separate
  (`integrator_ias15.c:435`) and runs **only** under
  `calculate_megno || (additional_forces && force_is_velocity_dependent)`,
  so a port that predicted always would hand a velocity-*in*dependent
  routine a different `v`.
- **Collisions**, bit-identical across a removal - see entry 30.
- **Test particles and `testparticle_type`**, with a 198-configuration
  sweep against REBOUND's gravity directly (N 1..6 × N_active 0..N and
  unset × testparticle_type × three ignore-terms values), every one bit
  for bit, **155 of them differing from the unrestricted answer**.
- **All four step criteria.** The control: 199 of 200 steps differ from
  PRS23, first at step 1, `0.043995128069094162` against
  `0.015442074500284753`, reaching `t = 24.91` against `t = 8.97`.
- **The checkpoint carries the mark and the stranded tail.**
- **Python on Windows**, by the opposite loading rule.
- **The subprocess API honours the step criterion** it had been
  dropping.

### The force ceiling, which is worse than it was described

The brief said a force enters as an exact binary64 addend. It does not:
REBOUNDx does `ax += f` and the callback returns **only the binary64
total**. The delta form `a_wide += promote(out − in)` is not
bit-identical at binary64 - `fl(u + fl(w−u)) ≠ w`, with `u = 1`,
`w = 2⁻⁶⁰` giving 0 - and flips a `−0` to `+0` even for a do-nothing
routine. What ships is replacement narrowed by `memcmp`, so an untouched
component keeps its wide value exactly and a written one is binary64
from that node on. At binary128 a routine writing `ax` only leaves the
run **1.289e-13** from REBOUND's binary64 running the same routine.

### And a refusal that is not about a setting

A routine that writes anything but `ax/ay/az` is refused, because
REBOUND uses more through three channels: `x`/`v` at the n = 0 call
become the step's initial condition (`:296`, then `:311`), a velocity
written at a node survives into the next while the predictor is off
(`:434`), and a mass is picked up at the next node. Measured with the
detection compiled out: **15 of 21 values differ** in each of the four
channels. The step in progress is **discarded**, restoring the
invariant every other refusal here has.

### Two upstream defects, both reproduced rather than inferred

**`reb_binarydata_diff` writes a field header with no data** for a field
present in the old snapshot and absent from the new
(`binarydata.c:325-333`). End to end: **0xC0000005** in `gate_real`.
And a consequence the parcel did not name, found by the verifier: an
archive written by a build carrying round 2's three new scalars and
appended by one that does not is **unloadable** - REBOUND itself prints
"The binary file seems to be corrupted."

**`N_active > N` is a heap overread** in REBOUND's own gravity loop, so
it is refused rather than reproduced.

### A claim in the code that was not true

An earlier draft said an archive from a run that never shrank is "byte
for byte what this project wrote before these members existed". Measured
by a verifier: **10,590 → 10,783 bytes, 62 → 65 named `cft_` fields**.
The 44 alias blobs are correctly skipped; the three scalars are
unconditional. What holds is that 123 fields are byte-identical and an
older archive resumes bit-identically including the wide digest. Both
comments say that now.

Relatedly: the constants digest is **FNV-style, not FNV-1a** - the prime
is FNV's, the offset basis is FNV-1a's `14695981039346656037` with its
last decimal digit missing. Kept, because the value is archived; renamed,
because the name was the wrong part. It had been hand-copied into a
second file three hours after this repository wrote down its rule
against transcribed constants, and is `CFT_DIGEST_SEED` in
`src/ias15_limits.h` now.

### Two modes that were unusable as shipped

INDIVIDUAL divides by `at[k]` for every coordinate, so a **planar
problem** divides 0 by 0 on every z lane, every step; GLOBAL divides by
`|x|²`, so does a particle at the origin. REBOUND's `isnormal()` rejects
the result and the coordinate drops out of the maximum - its own
estimate showing through. The standalone program aborted on the flag,
which made both modes unusable on most problems. The abort mask now
drops INVALID and DIVBYZERO for those two only; an overflow there is
still a real fault, and `flags_seen` still reports everything.

### The suite

**15 verdict lines, 211 `ok` assertions**, `coverage: all 14 drop-in
rows` and `all 16 subprocess rows`. The new `gate_subprocess` is the
first coverage `make check` has ever had of `libcft_rebound.a` - a
shipped public entry point that nothing in the suite touched, which is
how a flag went unforwarded for a month.

### Still open

The **wide clock is not archived**: a plain binary128 checkpoint differs
by one binary64 ulp in `t` for 2 of 40 step splits tried, and it
reproduces on a build from before this round. A remove → checkpoint →
regrow at `state->accurate = 1` is still inexact, by construction.
Reading a cft archive from Python on **Windows** is ungated - the wheel
is MSVC/UCRT and mingw64 links `msvcrt`, two heaps, and the load path
frees blobs REBOUND's loader allocated.

## 32. The gate cache, and the control that decides whether to trust it

`make check` reached about an hour during round 2 - five case files in
`check_dropin` where there was one, a new `gate_subprocess` leg, and the
archive gates carrying the alias family. That is the wrong instrument
for "I edited a docstring, did I break anything", so `make check-quick`
now skips a leg whose inputs have not changed since it passed
(`tools/gate_cache.py`).

The design decision worth recording is what the key is over: **the
artifacts a leg executes**, its own script, the data it reads and
`$CFT_REBOUND_ARTIFACT` - never source files. A source list must
enumerate every header and flag, and missing one skips a leg that would
have failed. A binary is the closed-over result of all of them. Finding
the inputs by grepping each checker rather than assuming immediately
justified itself: `gate_subprocess` spawns `ias15_cft`, so its key names
**both** binaries, and a key covering only the gate's own executable
would have skipped it after the program it drives had changed.

### Three runs, measured

| run | result |
|---|---|
| cold, cache deleted | every leg executed, **zero** skips, 14 stamps written, about 18 minutes |
| warm, nothing changed | **14 of 14 skipped by name**, `rc = 0`, **5 seconds** |
| invalidated: a comment appended to `src/ias15_cft.c`, rebuilt | the legs that link it **re-ran and re-passed** |

The third is the one that matters, and it is the reason this entry
exists: a cache that never invalidates is worse than no cache, and it is
the same shape as the seven gates-that-cannot-fail entry 31 records - a
mechanism reporting success without doing the work.

`check_dropin`, `check_dropin --wide`, `check_equivalence`,
`check_program_engine`, `check_records` and `check_ensemble` all re-ran
and passed before a `timeout 900` on the harness cut the run short -
`Error 143`, SIGTERM, not a gate. So the invalidation is demonstrated on
six legs rather than twelve; the remainder were never reached.

### Two legs skipped, and both correctly

`constants` keys only on `tools/gen_constants.py`, which was untouched.

`check_archive` is the informative one, because it was the specific
outcome that would have meant an under-inclusive key. It does not
re-run, and it should not: `gate_restart`, `gate_write`, `gate_stock`
and `gate_promote` link `tests/cft_shim_stub.c` plus `src/cft_archive.c`
and `src/cft_ias15_fields.c` - a deliberate **stub** shim, so the
archive is testable without the engine. `src/ias15_cft.c` is not in
those binaries at all, so touching it cannot change their behaviour. The
per-leg key was right rather than merely broad.

## 33. --light, and the banner that overstated it

`make check-light` runs every case with most step counts scaled to a
tenth, floor twenty: **559 s end to end** against the full suite's hour,
and **159 s against 585** on `check_dropin` alone, the same 130 cases,
zero failures. It is a fast signal from everything, and never a
substitute for `make check`.

Three things it had wrong first, every one of them found by measurement
rather than by design.

### Two cases must opt out, not one

The merge cases fire at step 1291 of 2000; at 200 the merge never
happens. Both refused to pass vacuously - `no collision occurred (N is
still 3)` - which is the guard those cases have carried since they were
written, and the third time in one day it earned its place.

The binary64 exemption went in first. The binary128 one was missed
because the probe ran `check_dropin --light` and never `--wide --light`,
and the first end-to-end `make check-light` failed on exactly that, at
`rc = 2`. Three green component measurements did not add up to a green
target: **a suite is only as tested as the composition you actually
ran.** Same shape as entry 30's seam - every part correct, the
combination untested.

### One family is not scaled, and must not be

`dt_sequence` (`tools/cases_modes.c`) drives `reb_simulation_steps(r, 1)`
inside its own `for (i < steps)` loop, and `light_steps(1)` returns 1,
so the step-sequence cases run at **full length** under `--light`. That
is correct: their subject is where two adaptive criteria diverge over
many steps, and a tenth of the steps is a different question. The
integrate-to-a-time cases are likewise untouched - they name a time, not
a count. Both were already truthful, and "fixing" either would have made
it lie in the other direction.

### The banners overstated the work

A banner is the only place in this suite that states a step count
without executing it, so under `--light` it was the only place that
could overstate the work done. Every one of them printed the declared
count:

    kepler, fixed step (N = 2, dt = 0.05, epsilon = 0, 400 steps)
      ok   kepler, fixed step: 2 particles x 9 values and the clock, bit for bit

over a run of 40. The suite exited 0, all 130 cases passed, and the
leg's own 159-against-585 seconds was the only evidence to the contrary.
Fourteen banners now print `light_ran(steps)` - the count that runs.
The three families above are excluded by hand, having been right
already.

**The control.** Without `--light`, `light_ran(n)` is `n`, so the whole
non-light log must be byte for byte what it was before the change.
Captured before, diffed after: `check_dropin` and `check_dropin --wide`,
both **identical**. That is what separates a printf change from a
behaviour change, and it is cheap enough that there was no excuse for
not running it.

## 34. A bitstream without binary256: what the rung cost, what the image refuses, and what IAS15 made of it

**2026-09-13, amd-arc-box (Ubuntu 24.04, 36 threads, 46 GB, Vitis
2022.2, Alveo U50C at `02:00.0` on the GA 6.8 kernel), cft-fp256 at
`ca19fe3` (ABI 0.12, XRT backend), REBOUND at `bdfda4bd`. Unlike every
entry above this one, the numbers here are from a card.** The build
flow is this repository's, in `hw/`; nothing in cft-fp256 was modified,
and every manifest records its `rtl/` tree hash so that claim is
checkable rather than asserted.

The occasion is the precision entries above and docs/HORIZON.md:
binary128 recovers the four orders binary64 loses to round-off, and
binary256 matches binary128 to every printed digit at every step a user
of this integrator would take. So binary256 is the rung this workload
can do without, and it is the most expensive one the tile carries. This
entry is what came of leaving it out - a **performance** artifact for a
port whose correctness was settled on the full image first.

### The packaging gap, and a check that caught its own mistake

cft-fp256's `rtl/cft_krnl.sv` has carried `EN_FP64`, `EN_FP128` and
`EN_FP256` since 2026-08-29, and its `docs/LAYOUTS.md` publishes a
catalogue of narrow layouts built on them. None of them had ever been
built, and the reason is one line of Tcl: `hw/package_kernel.tcl`
removes **every** user parameter from the packaged IP, so the `.xo` it
writes can only ever carry the RTL defaults. The generics were real,
reachable from simulation, and unreachable from a bitstream.

`hw/package_variant.tcl` here sets the value on the packaged IP's HDL
parameter *before* that removal. `hw/verify_variant_xo.tcl` then
instantiates the packaged IP in a throwaway project and prints the
overrides its synthesis wrapper hands to `cft_krnl`, because "the `.xo`
was written" and "the `.xo` carries the generic" are different claims
and only the second one matters:

    HDLPARAM: EN_FP64 = 1  EN_FP128 = 1  EN_FP256 = 0
    WRAPPER_PARAM: .EN_FP256(1'B0),

**That check refused the first correctly packaged `.xo` it was given.**
Vivado 2022.2 writes the sized bit literal with an uppercase `B` and
the normalisation knew only `1'b0`. It cost ten minutes, and it is the
direction a check of this kind should fail in: its false negatives cost
minutes, its false positives cost a two-hour link of the wrong tile
that nothing downstream would contradict.

### The single tile: 150 MHz, and a quarter of the tile gone

One compute unit, `EN_FP256=0`, cft-fp256's card-day recipe unchanged
(`RETIMING=1 PHYS_OPT=1`, default directives). 135 minutes to link.
`hw/verify-image.sh`: **8 of 8 PASS**. sha256 `90457667...`,
36,646,752 bytes, staged at `~/cardday-f128s/cft_hw_f128_1x.xclbin`.

| | full tile (rev-4 single) | f128 |
|---|---|---|
| kernel clock asked | 135 MHz | **150 MHz** |
| kernel WNS | +0.210 ns | +0.027 ns |
| failing endpoints | 0 of 129,804 | **0 of 105,669** |
| implied path delay | 7.197 ns | 6.640 ns |
| CLB LUTs, placed | 252,733 | 220,335 |
| CLB registers | 230,610 | 216,079 |
| DSPs | 311 | 168 |
| Block RAM tiles | 258 | 258 |
| URAM | 12 | 12 |

Net of the 123,897-LUT shell the tile is **96,438 LUTs against
128,836**: the binary256 rung was **32,398 LUTs and 143 DSPs**, a
quarter of the tile, against cft-fp256's own model of 35,333. Block RAM
and UltraRAM do not move at all, which is the expected result stated as
a check: they belong to the engine and the sequencer, which are
`BEAT_BITS`-wide and do not shrink with the rungs.

cft-fp256's `docs/LAYOUTS.md` predicted 150 MHz for this variant and
marked it "target, unmeasured". It is measured now and the prediction
was right to the megahertz.

**+0.027 ns is arrival, not headroom.** Vivado works a path exactly as
hard as the constraint asks and stops, so a closing build always lands
just above zero. The number that transfers between builds is the path
delay: 6.640 ns here against the full tile's 7.197 at its own ask. Read
150 MHz as at or near this tile's single-unit ceiling.

### The tile count, derived from a linked pair rather than a model

`hw/fit.py`. The per-compute-unit cost beyond the tile itself comes
from cft-fp256's **own linked revision-4 pair** - single 252,733 and
quad 640,500 placed, so three extra units cost 387,767 and each carries
**419 LUTs** of crossbar. (cft-fp256's `hw/gen_layouts.py` uses 12,626,
differenced from 2026-09-02 builds; the revision-4 pair is the same RTL
generation and the same tools as this image and is the better
calibration for it.)

| tiles | HBM PCs | model LUT | of device | verdict |
|---|---|---|---|---|
| 4 | 16 | 510,906 | 58.7% | fits |
| 5 | 20 | 607,763 | 69.8% | fits |
| **6** | **24** | **704,620** | **80.9%** | **tight** |
| 7 | 28 | 801,477 | 92.0% | no |

Six is the maximum worth attempting: it sits where the full-tile quad
closed (80.6%) and inside the HBM wall at 24 of 32 pseudo-channels.
Seven is past cft-fp256's 85% practical routing limit.

**The invocation first written into docs/BITSTREAM.md was circular** -
it passed this script's own predicted four-tile total back in as the
measurement from which the per-unit cost is derived, so the table would
have proved only that the arithmetic is self-consistent. `--per-cu`
exists now so that cannot happen silently.

### What it refuses, at three layers, each one measured

The image publishes `CAPS[3:0] = 0x7`. `ias15_cft --probe` on it:
`formats: fp32 fp64 fp128`, contract `0x00000800`, and every
revision-4 sequencer capacity intact - 64 deposits, 16,384
instructions, 512 constants, 256 scratch slots, `seq_features 0xf1f`,
resident buffers.

- **The tile, for itself.** `cft-resident -f fp256` issues `MODE`
  precision 3 straight at the kernel with no capability check. The run
  comes back **`status 0x8`** - `STATUS[3]`, the refusal bit - with the
  hardware's flags clean and its bytes differing from software, because
  the tile computed nothing.
- **libcft.** `cft_run(FMA, fp256)` and `cft_reduce(SUM, fp256)` both
  answer **status 2, "operation or format not available on this
  device", with the output buffer untouched** - not zeroed, which would
  be the worst available shape of a wrong answer. `cft_last_error()` is
  **empty**, which is a real gap and is filed below.
- **This repository.** The standalone program exits 3 and the drop-in
  fails the simulation, both with the same sentence:

      ias15_cft: .../cft_hw_f128_1x.xclbin carries fp32 fp64 fp128
      (CAPS[3:0] = 0x7, 1 tile, contract 0x00000800) and this run asked
      for fp256. The image publishes what it implements and libcft
      refuses the rest before a byte is issued, so nothing was computed.

  Before this session that refusal arrived at the *first operation* as
  a bare status and printed `fma: operation or format not available on
  this device ()`, naming neither the artifact, nor what it carries,
  nor what was asked.

**The control matters more than the refusal.** `hw/refusal_probe.c`
run against the full-tile revision-4 image reports binary256 **running**
there, output written, on the same card in the same session. Without
that, "refused" is equally consistent with a broken probe.

### What it computes

| gate | result |
|---|---|
| `cft-selftest`, cft-fp256's published vectors | 111 sets, **892,548 cases, all matching**, 168 s; binary256 sets skipped by name |
| `device-test`, three modes | 813, 2004, 687 checks, **0 failed** |
| `gate_real --fp64` | **PASS**, 170 s |
| `gate_real --fp128` | **PASS**, 253 s |
| `check_dropin --light` | **every case passed**, all 14 drop-in refusal rows exercised, 3,062 s |
| `check_program_engine fp128` | **PASS**; the FMA form against REBOUND's rounding sequence differs by at most 3.245e-34 |

892,548 rather than the full tile's 1,071,635 because one format of
four is gone. Nothing else in the census changed.

### The rate: the clock reaches the arithmetic

Both images run back to back on the same card in the same session with
the same tool, rather than against a number quoted from cft-fp256's
record. `cft-resident`, fma, n = 1,048,576, 20 reps, one compute unit:

| format | full tile @135 | f128 @150 | ratio |
|---|---|---|---|
| fp32 | 808.92 M/s (101.11 Mbeat/s) | 897.27 M/s (112.16) | **1.109** |
| fp64 | 417.23 M/s (104.31 Mbeat/s) | 462.03 M/s (115.51) | **1.107** |
| fp128 | 211.96 M/s (105.98 Mbeat/s) | 234.95 M/s (117.47) | **1.108** |

The clock ratio is 150/135 = **1.1111**. The measured gain is 99.7 to 99.8% of
it at all three rungs, so the extra megahertz turn into arithmetic and
removing the rung introduced no new wall. Every row: `status 0x0`,
bytes identical to the software backend, identical across repeats and
across units - **and the per-format digests are identical between the
two images** (`0ab39ba8b1e1ef3a`, `b4dd984b36bc7e7c`,
`ca3e223f7f0acad8`), which is the determinism claim holding across a
change of silicon.

Run-to-run spread is about 1.5%: the same f128 rows measured 883.93 /
458.29 / 234.23 M/s an hour earlier in the card series. The back-to-back
pair above is the controlled one and is what the ratio is taken from.

**The staged path does not move and cannot.** `cft-bench` on the same
image: fp32 169,914,652 elem/s at 2,718.6 MB/s, fp64 95,645,632 at
3,060.7, fp128 46,201,374 at 2,956.9 - the same 2.7 to 3.1 GB/s band
every image of this project has measured, because that path stages
operands across PCIe on every call. **That is the path `cft_run`
issues, and therefore the path IAS15 runs on.** Power at the tile was
15.0 W, FPGA 37 C.

### What IAS15 made of it: nothing

Wall clock against the software backend on one core, program engine,
FMA form, `--max-iter 60`. The corrector's mean pass count is beside
every ratio, because bit-identity alone would not prove the problem
posed was the one meant (CLAUDE.md's first way this suite has lied).

| problem | steps | format | software | card | ratio | passes |
|---|---|---|---|---|---|---|
| Kepler, N=2 | 200 | fp64 | 2.61 s | 22.02 s | **0.12x** | 3.43 |
| Kepler, N=2 | 200 | fp128 | 6.67 s | 40.25 s | **0.17x** | 8.56 |
| outer, N=6 | 100 | fp128 | 8.67 s | 24.20 s | **0.36x** | 8.42 |
| nbody, N=64 | 20 | fp64 | 24.19 s | 8.77 s | **2.76x** | 3.20 |
| nbody, N=64 | 20 | fp128 | 69.63 s | 22.49 s | **3.10x** | 7.95 |

Every record byte-identical to the software backend's in every value
and every counter.

That is the same crossover and the same plateau docs/HARDWARE.md
measured on the full tile, within the noise. **Removing a rung did not
move either, and could not have:** the crossover is set by per-call
staging cost and by the scatter in gravity's accumulate half, neither
of which is a function of how many rungs the silicon carries. A tile
that is 1.108x faster at arithmetic it is not waiting on buys nothing.

So the result of this entry is two numbers that must be read together.
**1.108x** where the tile's own rate is the thing measured, and
**1.00x** where IAS15's wall clock is. The image is a real improvement
to the coprocessor and not yet an improvement to this integrator; what
would change that is the library asks docs/HARDWARE.md already ranks -
resident state and a device-side scatter - not more megahertz.

### A comparison that could only ever fail

The card series' first run printed **"records differ" on all five
rows** of the table above while every data line was byte-identical. The
comparison was `cmp` on whole files, and a record carries a header
naming the backend and a trailer carrying elapsed seconds. Neither can
match across two runs, so the check had never been capable of passing -
and it was reporting a divergence in the one claim this whole port
exists to make.

Normalising exactly those two fields and nothing else: **0 lines differ**
at fp64 and at fp128. Both scripts now compare the record rather than
the file, keep the physics counters in the comparison, and report *how
many* lines differ rather than only that they do. Same disease as the
`--light` banner in entry 33: the check measured something adjacent to
the claim, and passed.

### Open at the time of writing

- **The six-tile image is in placement, not a result.** Launched
  2026-09-13 12:43 at 140 MHz - not 150, because the full tile lost
  0.188 ns of path going from one unit to four at a fixed ask, and six
  units at 80.9% occupancy should cost at least as much. If it closes,
  six tiles at 140 MHz against four at 135 is 1.56x the aggregate
  tile-clock product; whether IAS15 collects any of that is exactly the
  question the table above answers "no" to for the clock, and
  `hw/compare-images.sh` is written to answer it for the tiles. A
  six-unit `place_design` on this box is also the step that has been
  OOM-killed before, and an OOM kill reads like a design failure.
- `make check` **fails its three binary256 legs on this image, by
  refusal, and that is correct**: `make check` never skips, and an
  image that cannot do binary256 does not pass a binary256 gate.
  `hw/cardtest-f128.sh` runs the binary64 and binary128 gates and
  demonstrates the refusal instead.
- Reported to cft-fp256, not changed there: the packager strips
  generics; the XRT backend opens compute units by the literal name
  `cft_krnl:{cft_krnl_N}`, so the variant kernel names in its own
  `hw/layouts/*.cfg` would open zero tiles; a format refusal sets no
  `cft_last_error`; no bench simulates a trimmed **full-beat** tile
  (`test_krnl_quarter` is `BEAT_BITS=64`), so this image's RTL
  configuration had never been simulated before it was linked; and
  fp32 has no generic, so the fp64/fp128-only tile this workload would
  actually want cannot be expressed.

## 35. Six tiles do not close at 140 MHz: the area model was right and the timing extrapolation was not

**2026-09-13, amd-arc-box, same toolchain and pins as entry 34, which
this entry closes the open item of.** Entry 34 left the six-tile image
"in placement, not a result". It is a result now, and the result is a
failure: **389 minutes, kernel WNS -0.558 ns, no image produced.**

### What was asked and what happened

Six compute units of the f128 tile at **140 MHz**, cft-fp256's card-day
recipe (`RETIMING=1 PHYS_OPT=1`, default directives). Placement
succeeded and was not memory-killed; routing ran to completion and
missed.

    kernel_wns_ns: -0.558
    kernel_failing_endpoints: 3015 of 632233
    violating_clocks: clk_out1_ulp_clk_wiz_0 -0.558 3015/632233;
                      hbm_aclk -0.005 9/103249

**3,015 failing endpoints of 632,233 is 0.48% - broad and shallow.**
That is the signature of routing congestion rather than one pathological
path: a single bad path fails by a lot at a handful of endpoints, and
this failed by half a nanosecond across three thousand of them. The
`hbm_aclk` violation at -0.005 ns on 9 endpoints is the shell's own and
is the same one Vitis has auto-scaled away on past builds; it is not
this design's problem.

The implied path delay is 7.143 + 0.558 = **7.701 ns**, so six of these
tiles top out at about **130 MHz**, not 140.

### The area model was accurate

Predicted from the linked single and cft-fp256's own revision-4 pair
(entry 34's table), against what actually placed:

| resource | predicted | placed | error |
|---|---|---|---|
| CLB LUTs | 704,620 (80.9%) | 698,188 (80.19%) | 0.9% pessimistic |
| Block RAM | 643.8 | 642.5 (47.81%) | 0.2% |
| URAM | 52 | 52 (8.13%) | exact |
| DSPs | 988 | 988 (16.60%) | exact |

Three of four exact and the fourth within one percent, which also
confirms the tile really did replicate six times rather than something
being shared or dropped. The super-logic-region split was **84.17% in
SLR0 against 75.96% in SLR1**.

### The timing extrapolation was not, and the reason is the interesting part

The estimate that chose 140 MHz came from cft-fp256's full tile, where
going from one compute unit to four cost **0.188 ns** of path at a
fixed ask (+0.210 at one unit, +0.022 at four, both at 135 MHz). Six
units being more crowded than four, that was scaled to an assumed
0.2 to 0.3 ns and taken off the f128 single's 6.640 ns path, giving a
predicted six-unit ceiling of 144 to 146 MHz and 140 as the safe ask.

The actual penalty was about **1.06 ns**, three to five times the
assumption.

**Area is a sum and timing is not.** The full tile's four units sit at
73.6% of the device; these six sit at 80.19%. Six and a half points of
occupancy turned a 0.188 ns penalty into roughly a nanosecond, because
the cost of congestion is not linear in how full the part is - it is
flat until the router has room to spare and then it is not. Every LUT
figure in entry 34 extrapolated correctly because adding a tile adds
its LUTs; the timing did not, because adding a tile does not add a
fixed quantity of delay, it changes how hard every other tile is to
route.

So the transferable rule, and the one this project had already written
down in another form: **a per-unit penalty measured at one occupancy
does not predict timing at a higher one.** cft-fp256's "read the path
delay, not the slack" says the tool works exactly as hard as the ask;
this adds that how hard it *can* work depends on how full the die is,
and that relationship has to be sampled rather than extrapolated. Two
points at different occupancies would have been worth the six and a half
hours this cost, and one point plus an assumption was not.

### What was done about it

Relaunched at **125 MHz**, 19:16, same six tiles and the same recipe.
125 asks for 8.000 ns against the 7.701 ns the router actually
achieved while failing, so there is 0.3 ns of real margin rather than
an extrapolation. 130 MHz was rejected deliberately: it asks for
7.692 ns, nine picoseconds *under* a path the tools have already
demonstrated, and a coin flip is not worth six and a half hours.

If it closes, six tiles at 125 MHz against the shipped four at 135 is
**1.39x** the aggregate tile-clock product rather than the 1.56x
entry 34 projected at 140. Whether IAS15 collects any of it is the same
open question, and entry 34's wall-clock table is the reason to expect
very little: the clock bought 1.108x of arithmetic rate and 1.00x of
IAS15 wall clock, because per-call staging and the scatter in gravity
set that number and neither is a function of the silicon's rung count
or its clock.

The failed build's manifest is kept. It records the timing it reached
and is the only evidence of this tile's six-unit ceiling.
