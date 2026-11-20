# Validation ledger

What ran, against which commits, and exactly what it said - passes and
failures alike, numbers verbatim from the runs. One entry per
campaign, newest last; a correction is a new entry naming the old one.
The style is cft-fp256's docs/VALIDATION.md, on purpose.

Host for the entries that do not say otherwise (every entry from 28 on
names its own host, and from 34 on its own cft-fp256 pin): DESKTOP-class
Windows 11 box, MSYS2 mingw64 gcc 16.1.0, `-std=c99 -O2`, Miniconda
Python 3.12.9 with mpmath 1.3.0 (on gmpy2), libcft from cft-fp256 at
40ce35e3 (ABI 0.11, software backend only), REBOUND at bdfda4bd
(5.1.1 + 6). No card, no XRT, no `cft://` server for those entries:
both were in use, and everything in them is the software backend, which
is the contract's definition of the bits. The pin moved to ca19fe3
(ABI 0.12) on 2026-09-13 (entry 34) and to 922c2ad (ABI 0.18) on
2026-10-30 (entry 38); third_party/MANIFEST is the authority and
`make check-pins` holds the tree to it.

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

## 36. The six-tile image closes, and tile count turns out to be a cost: what the specialised bitstream is actually worth to IAS15

**2026-09-14, amd-arc-box, Alveo U50C, same pins as entries 34 and 35.**
The six-tile image closed at 125 MHz and was measured against the
shipped four-tile full image and against the one-tile images of both
kinds. The headline is not the one this work was started for.

### The image

Six compute units, `EN_FP256=0`, **125 MHz**, 322 minutes, **kernel WNS
+0.061 ns**, `hw/verify-image.sh` 8 of 8 PASS including *24 masters, all
channels HBM, no channel shared*. sha256 `51f00fb6...`, staged at
`~/cardday-f128x6b/cft_hw_f128_6x.xclbin`. Placed at 696,806 LUTs
(80.03%), slightly under the 140 MHz attempt's 698,188, which is what a
relaxed constraint should do. `ias15_cft --probe`: 6 tiles, `CAPS[3:0]
= 0x7`, contract `0x800`; binary256 refused by name, `device-test`
813 checks 0 failed.

### The arithmetic scaled exactly as designed

`cft-resident`, fma, n = 1,048,576, all units engaged:

| format | quad @135, 4 units | f128 @125, 6 units | ratio |
|---|---|---|---|
| fp32 | 3,224.04 M/s (51.6 GB/s) | 4,476.41 M/s (71.6 GB/s) | 1.388 |
| fp64 | 1,664.39 M/s (53.3 GB/s) | 2,308.92 M/s (73.9 GB/s) | 1.387 |
| fp128 | 845.20 M/s (54.1 GB/s) | 1,177.71 M/s (75.4 GB/s) | 1.393 |

Predicted 6x125 / 4x135 = **1.389**; measured 1.387 to 1.393. Per-unit
rate is 93.3 to 98.1 Mbeat/s against the quad's 100.8 to 105.7, exactly
the 125/135 clock ratio. Digests identical to the single-tile images and
to software at every format, so determinism holds across one, four and
six tiles and across two different bitstreams.

**Everything about the hardware worked.** What follows is about the
workload.

### More tiles is slower. Monotonically, everywhere.

`hw/compare-images.sh`, program engine, FMA form, `--max-iter 60`,
ratios against the software backend on one core. Every record on every
row **byte-identical** to software; the corrector's pass count matches
what each format needs.

| problem | steps | format | f128 x1 @150 | quad x4 @135 | f128 x6 @125 |
|---|---|---|---|---|---|
| Kepler N=2 | 100 | fp64 | 0.09x | 0.06x | 0.06x |
| Kepler N=2 | 100 | fp128 | 0.13x | 0.09x | 0.09x |
| outer N=6 | 100 | fp64 | 0.16x | 0.06x | 0.05x |
| outer N=6 | 100 | fp128 | 0.25x | 0.09x | 0.08x |
| nbody N=64 | 20 | fp64 | **1.47x** | 0.63x | 0.48x |
| nbody N=64 | 20 | fp128 | **2.23x** | 1.00x | 0.79x |
| nbody N=256 | 5 | fp64 | **2.99x** | 1.93x | 1.60x |
| nbody N=256 | 5 | fp128 | **3.80x** | 2.61x | 2.22x |
| nbody N=512 | 3 | fp64 | **2.94x** | 2.34x | 2.07x |
| nbody N=512 | 3 | fp128 | **3.91x** | 3.16x | 2.86x |

**One tile beats four beats six, in all ten rows, without a single
exception.** The six-tile image carries 1.39x the aggregate arithmetic
rate of the shipped quad and integrates about 10% *slower* than it.

The mechanism is not mysterious. libcft partitions every `cft_run`
across all the tiles a device presents, so one library call becomes one
kernel launch and one staging round **per tile**. IAS15 issues thousands
of small calls per step - 51,911 for twenty steps at N=64 and binary64 -
and at these vector lengths the per-call cost multiplied by the tile
count swamps whatever the extra arithmetic contributes. Tile count is
not a throughput knob for this workload; it is a divisor on an already
small n and a multiplier on the fixed cost.

### The control: what the specialised image is actually worth

The table above cannot separate "fewer tiles" from "higher clock and
one rung fewer", because its winner differs from the quad in both. So
both **one-compute-unit** images were run head to head, same card, same
session, card runs only:

| problem | steps | format | f128 @150 | full tile @135 | ratio | records |
|---|---|---|---|---|---|---|
| nbody N=64 | 20 | fp64 | 13.47 s | 13.62 s | **1.01** | identical |
| nbody N=64 | 20 | fp128 | 29.21 s | 29.72 s | **1.02** | identical |
| nbody N=256 | 5 | fp64 | 31.98 s | 32.15 s | **1.01** | identical |
| nbody N=256 | 5 | fp128 | 63.61 s | 64.90 s | **1.02** | identical |
| nbody N=512 | 3 | fp64 | 75.41 s | 77.67 s | **1.03** | identical |
| nbody N=512 | 3 | fp128 | 148.09 s | 151.51 s | **1.02** | identical |

**One to three percent.** That is the whole value of this bitstream to
this integrator.

So the 1.24x to 2.67x by which the f128 single beat the shipped quad was
**tile count, essentially all of it, and not the specialised image**.
Loading cft-fp256's existing one-tile full image - which still carries
binary256 - recovers 97 to 99% of it. The correction matters because the
three-image table invites exactly the wrong reading, and I published
that reading before running this control.

Why so little, when the arithmetic is 1.108x faster (entry 34)? Because
the arithmetic is not what the wall clock is made of. An 11% faster
engine moved the integration by 2%, which puts the arithmetic at
roughly a fifth of it; the rest is per-call staging across PCIe and the
scatter in gravity's accumulate half, and neither is a function of the
clock or of how many rungs the silicon carries. Entry 34 measured this
at one tile and called it "the crossover did not move". It is the same
fact, now with the size of the effect on it.

### What this says to do

- **For IAS15 on this card, use a one-tile image.** Either one. The
  six-tile image is the wrong artifact for this workload and the
  four-tile shipped image is worse than one tile at every size tested.
- **A multi-tile image is still right for a workload with long vectors
  or independent work per tile** - the arithmetic scaling above is
  perfect, 1.39x for 1.389x asked. IAS15 is simply not that workload at
  any body count reachable here.
- **The lever is the library, not the silicon.** docs/HARDWARE.md
  already ranks the asks - state resident on the card across a step, and
  a device-side scatter for gravity's accumulate half - and this entry
  is the measurement that says nothing in the bitstream substitutes for
  them. A tile that is 1.39x faster at arithmetic it spends four fifths
  of its time not doing is worth 2%.
- **The cheapest experiment would have predicted all of it.** Running
  the existing one-tile and four-tile images against the workload costs
  an hour and needs no build; it would have shown tile count to be a
  cost before thirteen hours went into two multi-tile links. The signal
  was on the record too: the ensemble measurement in this ledger already
  found many small systems below one large one, which is the same effect
  seen from the other side. **Measure the partitioning before buying
  more of it.**

### What the exercise did buy

Not nothing, and worth separating from the disappointment: cft-fp256's
trim generics now reach a bitstream at all, and its `docs/LAYOUTS.md`
catalogue is no longer entirely placeholders - two of its rows are
built, one of them at a clock the catalogue had guessed and never
measured. The refusal path is real, exercised and controlled at three
layers. The area model is calibrated against two linked multi-tile
images. And the six-unit timing ceiling of this tile, roughly 130 MHz,
is now a measured number rather than an extrapolation. The bitstreams
are staged; what they are worth to IAS15 is 2%, and that is the number
to quote.

## 37. Five benchmarks on a real integration: where the hardware pays, and where a step's time actually goes

**2026-09-14, amd-arc-box, same pins as entries 34-36; the port at
`fbd1532`.** Everything below is an actual IAS15 integration of an
actual problem - not an elementwise microbenchmark - and the baseline
row of every table is REBOUND's own IAS15 in hardware double, the thing
a person runs today. Five datasets, one page drawn from them by
`hw/report.py`, and every card row held byte-identical to its software
counterpart. The tools: `hw/bench-workload.py`, `hw/bench-ensemble.sh`,
the trailer instrumentation in `src/ias15_cft.c`, and `hw/bench-modes.sh`
for the mechanism behind them.

### The workload table: 60 rows, 0 failures

Kepler (2 bodies, 200 steps), the outer solar system (6, 100), Burrau's
Pythagorean problem (3, 100), and n-body at 64, 256 and 512 (20, 5, 3
steps); REBOUND in double, the port in REBOUND's own rounding order at
binary64, the port in the FMA form on the program engine at binary64
and binary128 on software and on the one-tile f128, four-tile full and
six-tile f128 images.

**The port is REBOUND bit for bit, on every problem** - the six
`sw-rebound-form` rows are identical to `ias15_ref` on every value the
equivalence gate compares - **and that costs about 2,000x at
binary64**: 248.7 s against 0.11 s at 512 bodies, 4.9 s against 4 ms
on Kepler. That number is the honest price of the contract in
software, and it is why the case for the card lives at binary128.

| bodies | binary128, software | one f128 tile | four full | six f128 | one tile vs software |
|---|---|---|---|---|---|
| 2 | 7.01 s | 49.97 s | 70.20 s | 75 s | 0.14x |
| 64 | 65.22 s | 29.18 s | 62.65 s | 81.96 s | **2.24x** |
| 256 | 238.28 s | 64.53 s | 90.61 s | 108.09 s | **3.69x** |
| 512 | 577.97 s | 147.64 s | 181.78 s | 201.55 s | **3.91x** |

One tile beats four beats six in every row again, and at binary64 the
one-tile image is the first hardware row in this repository to beat
the software backend on a single system: 13.75 s against 19.94 s at 64
bodies.

### Accuracy for cost

The energy-drift column, once it was computed exactly (the first run
computed it in doubles and printed 0 for every binary128 row - the
accuracy column measured with the instrument whose limitation it
exists to show):

| problem | binary64 (REBOUND and port) | binary128 | binary256 |
|---|---|---|---|
| Kepler, 200 steps | 1.30e-15 | 9.41e-21 | 9.41e-21 |
| outer, 100 | 4.11e-16 | 3.45e-25 | 3.45e-25 |
| Pythagorean, 100 | 1.38e-13 | 1.70e-20 | 1.70e-20 |
| n-body 64, 20 | 1.04e-14 | 2.57e-28 | - |
| n-body 256, 5 | 4.03e-14 | 1.58e-28 | - |
| n-body 512, 3 | 6.78e-14 | 8.40e-29 | - |

Binary128 and binary256 agree to every printed digit on the three
small problems, which is docs/HORIZON.md's result seen again from the
energy: at these step sizes the method's truncation error is what
remains, and binary128 is already below it. So the accuracy the card is
bought for is binary128's, and on the frontier the cheapest binary128
point at every body count from 64 up is the single f128 tile.

### The ensemble: the workload with long vectors

E independent Kepler systems, members one ulp apart on a dyadic offset,
fixed dt 0.05, 20 steps, the ledger's own integration. System-steps a
second:

| E | format | software | one f128 tile | four full | six f128 |
|---|---|---|---|---|---|
| 1 | binary128 | 29.9 | 2.5 | 1.5 | 1.3 |
| 8 | binary128 | 85.1 | 18.8 | 7.3 | 6.7 |
| 64 | binary128 | 108.5 | **117.1** | 53.2 | 40.3 |
| 512 | binary128 | 109.3 | **418.5** | 263.0 | 215.4 |
| 4,096 | binary128 | 103.3 | **562.3** | 493.2 | 454.6 |
| 4,096 | binary64 | 243.2 | **1,405.6** | 1,176.5 | 1,063.5 |

The card crosses the software backend at about **64 members** on one
tile and reaches **5.4x at binary128 and 5.8x at binary64** at 4,096.
And the multi-tile images lose here too - at 24,576 lanes, the widest
vector this integrator can present - by 12% (four) and 19% (six) at
4,096 members. The margin narrows as E grows and never inverts: the
ensemble's call count is independent of E by construction (entry 20),
so each extra tile adds its launch and staging round to every one of
those calls and the extra lanes never pay it back. Lane efficiency
0.77 to 0.97 throughout, so this is not a ratio on idle silicon.

*Verdict column: the first run of this sweep reported every card row
as differing from software in two lines. A diff of one pair showed the
only differing line was the trailer carrying the new wall-clock keys,
and the data-only diff was zero lines; the normalisers were fixed and
the sweep re-run: 30 of 30 card rows identical in every value and every counter, and the three binary256 quad rows likewise. The timings in the tables above are from the re-run.*

### Where a step goes: the measurement that reorders the residency work

The instrumented trailer, one f128 tile, share of wall clock:

| problem | format | program runs | elementwise | div and sqrt | gravity, inclusive |
|---|---|---|---|---|---|
| Kepler | binary64 | 8.1% | 77.3% | 13.0% | 59.8% |
| Kepler | binary128 | 10.5% | 75.9% | 12.2% | 71.7% |
| n-body 64 | binary64 | 3.7% | 63.6% | 30.6% | 82.7% |
| n-body 64 | binary128 | 4.3% | 54.8% | 39.0% | 89.9% |
| n-body 256 | binary64 | 0.9% | 48.5% | 48.0% | 73.5% |
| n-body 256 | binary128 | 1.1% | 37.4% | 59.2% | 85.8% |

**The program runs - where every byte the residency plan would keep on
the card actually moves - are 1 to 10% of the step.** docs/HARDWARE.md's
morning ranking (allocate through `cft_alloc`, keep the scratch block
resident, then the convergence test) was costed correctly and ordered
wrongly, and its own last paragraph, which said to measure before
building any of it, was the part that was right. Two things are the
wall instead:

- **Per-call cost at small sizes.** 172,000 six-element calls in a
  200-step Kepler run at **108 microseconds each** on the card against
  0.73 in software: 77% of the step. A resident six-element buffer still
  pays the launch. This is a call-count problem, and the fix is more of
  the substep as a program.
- **The correctly rounded divide and square root at large ones.** At
  256 bodies and binary128, 674 `cft_div`/`cft_sqrt` calls take 35.1 of
  59.3 seconds: 52 ms a call over 32,640 pairs, **1.6 microseconds an
  element against 4.3 nanoseconds for an FMA on the same tile**; in
  software 7.6 microseconds an element and 71% of the step. That is the
  composed correctly-rounded route (a program core with host prep and
  finish on every operand), not the tile - cft-fp256's own peers tool
  prices its software binary128 square root at ~10,900 ns against
  MPFR's 49 - and it is filed to cft-fp256 as the first ask now.

Bytes are not the wall: 2.4 GB a step at 256 bodies and binary128 is
about 0.8 s of an 11.9 s step at the measured bus rate. **Cut the call
count, then fix the divide, then residency.**

### Two more of my comparisons that could only fail

Recorded because the pattern is now four for four this week and worth
naming as one thing. (1) The workload table's binary64 software row
was held to bit-identity with REBOUND while running the FMA form, which
is documented to differ at round-off; there are two software rows now,
and only the REBOUND-form one is held to REBOUND. (2) Even that one then
"differed" on two trailing fields: the exact-time pair, a double-double
split of the same instant the two programs divide differently between
words, which `tools/check_equivalence.py` has always truncated at the
`|`. The benchmark now uses the gate's definition of agreement rather
than one invented beside it. In every case the disease was the same as
entry 33's banner: **the check measured something adjacent to the
claim, and the adjacent thing could not pass.** A comparison that has
never been seen to pass on known-identical input has not been tested.

### What the page shows

`hw/report.py` draws the five figures - seconds a step against bodies,
accuracy for cost, the ensemble, where a step goes, and the same-bits
matrix - from the CSVs under `bench-modes/`, with a table twin under
each; the series keep one hue per image across every figure and the
software backend is the de-emphasis gray. The digests are identical
across one, four and six tiles, two bitstreams and the software
backend, which is the property that lets a pooled node's answer be
trusted without recomputation, and the matrix is that claim drawn from
the rows rather than asserted.

## 38. The pin moves to cft-fp256 922c2ad, libcft ABI 0.18: the same bits from the same calls, and nothing new called yet

**2026-10-30, this desktop (Windows 11, MSYS2 mingw64 gcc 16.1.0,
`-std=c99 -O2`, Miniconda Python 3.12.15 with mpmath 1.3.0), REBOUND at
bdfda4bd, the software backend only; the port at `039e3c3`.** cft-fp256
moved 1,323 commits and six ABI minors past `ca19fe3` while this
repository sat idle from 2026-09-14, and by its own record built or
answered every ask docs/HARDWARE.md made of it on the way (its
docs/ROUND2.md; ROADMAP.md, "What entry 37 leaves to build", says which
ask became which call). This entry moves the pin and establishes one
thing: the calls this port already makes return the same bits through
the new library. It adopts nothing.

### Which commit, and why not HEAD

`922c2ade`, the newest commit cft-fp256's own gate has passed at:
`verify/run.sh --only docs,golden,generated,seq,cardstress`, rc 0,
VERDICT PASS, on 2026-10-30 (its docs/VALIDATION.md, the lead's entry of
that day; CI green on `77dfc230`, the base). It is parcel M2F2's tip;
its merge onto main, `d6f55748`, has the identical tree (`6dd047b0...`),
and HEAD that evening, `f68c6301`, is that tree plus the three documents
recording the gate.

The part of cft-fp256 this repository builds - `host/src`,
`host/include` and `host/Makefile`, compared by git tree hash - is
identical at `6f7ad1d5`, `0f9ff761`, `9dbfa501`, `922c2ade` and
`f68c6301`, and has not changed since 2026-10-28. So the pinned libcft
is also the one cft-fp256's night 1029 ran on the U50 (at `6f7ad1d5`:
1,071 legs, 0 failed, 2,552,252,429,426 cases, 0 disagreements) and its
night 1030 runs over this weekend. `rtl/` is `44b963dd...` at all five;
entry 34 built its image from `ca19fe3`'s, `0aef8eea...`.

### How it was run

Two worktrees from main at `039e3c3`, each with fresh clones of its own
pins: the baseline at `ca19fe3`, and the re-pin at `922c2ade`, fetched
from GitHub by `tools/fetch_third_party.sh` from the new MANIFEST line
(rc 0, both pins `ok`) - so the line reproduces from the URL it names.
`make check-pins` passes in both. Each then ran `make all programs
check`, the re-pin `make libcft` first, at the same time as the other.

- **The build.** libcft at the new pin builds in 12 s, `all` in 43 s
  (50 s at the old pin), `programs` in 6 s. Both build logs carry the
  same five warnings, all REBOUND's own (`simulationarchive.c`,
  `-Woverflow`); none comes from libcft or from this repository.
- **No source edit.** `cft.h` grew from 2,798 lines to 3,539, and two
  input structs grew with it (`cft_run_args` at 0.14 and again at 0.17,
  refused at an old `struct_size`). The port sets `struct_size =
  sizeof` at all three of its sites (`src/ias15_cft.c:115`, `:148`,
  `:1008`), so the rebuild is the whole migration.

### The suite at both pins

| | old pin `ca19fe3` | new pin `922c2ad` |
|---|---|---|
| `make check` | rc 0, 4,422 s | rc 0, 4,245 s |
| legs run | 14 | 14 |
| `ok` lines | 222 | 222 |
| `FAIL` lines | 0 | 0 |
| skips | gate_subprocess's `body_count` row, deliberately untested and named | the same |

Read from the logs, not the exit codes. With paths, temporary
directories and timings normalised, the two logs differ in exactly
three things: the pin `check-pins` names, the random temp-directory
names, and the archive gates' recorded `abi "0.12"` against
`abi "0.18"` - the archive's own provenance field doing its job. Every
verdict line, case and count is otherwise the same line in both, and
that includes the six wide-state digests `gate_real` prints for its
plain and regrow sequences - fp64 `e08541bc0af9b3ec` and
`56bacf3bf63f19a7`, fp128 `5d7cf4588fad050e` and `cd16e6dc29aa0589`,
fp256 `a1052750e206bebe` and `ba4d447c8e58bd20` - which hold the whole
wide state to the bit at every format across the two libraries, not
only its binary64 view.

### The same integrations through both libraries

The suite holds each pin to REBOUND at binary64 and to itself across
engines, members and checkpoints. Above binary64 it holds the new
library to the old one only through `check_records`' committed samples
and the digests above. So 27 integrations were run through both builds
and compared on what the gates compare: every non-comment line of the
record, byte for byte, and the trailer's counts - steps, rejections,
capped correctors, mean and maximum passes, flags, library calls. The
only comment lines that differ are the header's `impl=libcft-0.12`
against `impl=libcft-0.18`, and the timings.

| problem | formats | variants | identical |
|---|---|---|---|
| Kepler, fixed dt 0.05, 200 steps | fp64, fp128, fp256 | REBOUND's roundings; the FMA form on the loop engine; the FMA form on the program engine; `--cs augmented` | 12 of 12 |
| Kepler, adaptive 1e-9, 200 | fp64, fp128, fp256 | | 3 of 3 |
| outer solar system, fixed dt 40, 40 | fp64, fp128, fp256 | | 3 of 3 |
| Pythagorean, adaptive 1e-9, 100 | fp64, fp128, fp256 | | 3 of 3 |
| Kepler e = 0.99, adaptive 1e-9, 200 | fp64, fp128 | | 2 of 2 |
| 16 bodies (`tools/make_nbody.py`), fixed, 5 | fp64, fp128 | | 2 of 2 |
| a 4-member Kepler ensemble, adaptive | fp128 on the loop engine, fp256 on the program engine | | 2 of 2 |

**27 of 27**, and none of them trivially: the fixed Kepler run's
corrector averages 2.5, 6.2 and 15.0 passes at the three formats, so
each format took its own path; each e = 0.99 run carries a rejected
step and a capped corrector; the library-call counts run from 12,430 to
2,720,328 a run. The 27 program images `make programs` assembles are
byte-identical between the two pins' `cft-asm`.

**The control**, so that "identical" is a measurement: the same
comparison on one build, binary256 Kepler with the corrector capped at
REBOUND's 12 against 60, reports 4 of the 5 sample lines different
(mean passes 12.000 against 15.035) - ROADMAP's "footgun found by
building on it", seen by the comparison that has to see it.

The comparison's first pass needed one correction, and it was the
comparison's: it read `sample` lines only, an ensemble record writes
`system` lines, and both ensemble cases came back "DIFF, 0 samples" on
both sides. The guard against an empty comparison caught it before it
could read as a finding. Re-compared on every non-comment line, both are
identical, 12 `system` lines each.

### A checkpoint across the re-pin

What a user who upgrades actually does is save under one library and
resume under the other. The archive records `cft_abi`, but nothing in
the loader refuses on it (`src/cft_archive.c` reads it as provenance),
so the question is only whether the bits survive. `gate_real`'s
three-process form - `--straight` against `--save` then `--resume`, as
`tools/check_checkpoint.py` drives it - was run with the save and the
resume in different builds: written by the old pin and resumed by the
new, and the reverse, at fp64, fp128 and fp256, plain and regrow.
**12 of 12 identical** to the new build's straight run, all 15 lines of
each dump including its wide-state hash. The control, a resume with the
loaded wide state dropped (`CFT_REBOUND_NO_ADOPT=1`), differs in 6 of
the 15.

### The install path

`make example` at both pins, rc 0: the subprocess round trip asserting
its declared answers, and the drop-in at binary256 through a checkpoint,
"restored exactly". With paths normalised the two outputs are identical
except in one column, the round trip's wall seconds.

### What this does not establish

- **The card.** amd-arc-box's U50 is in cft-fp256's weekend card-stress
  window, Friday 18:00 to Monday 2026-11-02 07:45 (its rule of
  2026-10-24), so nothing here opened an artifact. The new library
  still lists contract `0x800`, the f128 image's (docs/BITSTREAM.md),
  among its `KNOWN_VERSIONS` (`host/src/backend_xrt.cpp:267`), so it is
  believed to open that image. That it does, and returns the software
  backend's bits there, is owed to the next card day.
- **Linux.** Entry 29's open item stands: no `make check` on
  amd-arc-box at either pin.
- **Python.** `make check-python` needs a `rebound` wheel this host's
  Miniconda does not have, and was not run at either pin.
- **The bitstream flow.** `hw/build-variant.sh` packages
  `third_party/cft-fp256/rtl`, which this pin changes. Not run: it needs
  Vitis on amd-arc-box.
- **Anything new.** No call this port makes is one the new ABI added.
  The index tables, `cft_reduce_seg`, `CFT_XRT_TILES` and the identity
  calls are in the pinned library and unused, as `scalar_mask` already
  was at the old pin.

## 39. Round 3's P0: one door to what libcft added for this port, sent back twice and changed in no record line; and entry 38's Linux and Python items

**2026-10-31.** Two hosts:
- this desktop: Windows 11, MSYS2 mingw64 gcc 16.1.0, `-std=c99 -O2`,
  Miniconda Python 3.12.15 with mpmath 1.3.0;
- cft2204, a WSL2 Ubuntu 22.04 distro on the same box: gcc 11.4.0,
  Python 3.10.12, mpmath 1.4.1.

Pins: REBOUND at bdfda4bd and cft-fp256 at 922c2ad (ABI 0.18), on the
software backend only. The port was at `ebaeaa1` for the Linux leg, and
P0's commits on branch `round3-p0` are `971f93d` and `20c05d8`. The fixes
on `round3-p0-v1` are `47f0bc3`, `21579ee`, `90e6f3f`, `d7d964b` and
`733b067`.

ROUND3.md is the round's plan and has P0's design. This entry is what
ran.

### Entry 38's two open items

**Linux, at `ebaeaa1`, in cft2204, from local clones only.**
- `make check`: rc 0, 14 legs, 222 `ok` lines (the same count as
  Windows at that pin), 0 FAIL, in 5,608 s on a loaded box.
- All six wide-state digests are identical to Windows': fp64
  `e08541bc0af9b3ec` and `56bacf3bf63f19a7`; fp128 `5d7cf4588fad050e`
  and `cd16e6dc29aa0589`; fp256 `a1052750e206bebe` and
  `ba4d447c8e58bd20`.
- `make example`: rc 0.

The leg found two warnings that the mingw build hides, both in this
repository's tests:
- `tests/gate_real.c:650` printed a `uint64_t` with `%llu`, which is
  `unsigned long` on LP64.
- `tests/gate_subprocess.c` called `putenv()` with no declaration under
  `-std=c99`. That is the same trap as the trailer's clock of
  2026-09-14, and gcc 14 or later rejects it outright.

Both are fixed in `90e6f3f`. At `733b067` the Linux build has no warning
outside `third_party/`.

This is the first full Linux suite in this ledger, but it was run in
WSL, not on amd-arc-box. So entry 29's item is answered for Linux and
not for that host.

**Python, on Windows.** The owner allowed the download, and rebound
5.1.1 went into the Miniconda env: the `cp312-win_amd64` wheel at
githash `33549d1d`. `check-rebound-match` found the pinned and installed
headers identical.

The first `check-python` failed. `ctypes` could not load
`build/libcft_ias15.dll` "or one of its dependencies". The DLL imported
`libwinpthread-1.dll`, which it had done since `24900ee` on 2026-09-14:
the trailer's `clock_gettime()` comes from winpthreads under mingw-w64,
and nothing had run the Windows Python leg since 2026-09-11. The link
rule printed the import table on every build, and nothing read it.

`47f0bc3` links winpthreads statically and turns the print into a
refusal of any `lib*.dll` import.
- Control: `make python-lib SHSTATIC=` is refused, naming
  `libwinpthread-1.dll`, with rc 2 and no library left behind.
- After the fix, `check-python` gives "IDENTICAL: 13 values, bit for
  bit", and `check-python-control` gives "DIFFERS: 9 of 13 values".

### What P0 is, in one paragraph

ABI 0.12 to 0.14 added three ways of issuing a call that this port had
never used: a scalar operand, an indexed gather and a segmented
reduction. Each sits behind a capability bit, and each has an older
route that the contract makes give the same bits.
- `engine_routes()` decides each family once, at open, in the program
  and the library half alike. The variables are `CFT_REBOUND_ROUTE_*`
  and `CFT_REBOUND_HIDE_FEATURES`, and every misuse is refused by name.
- Three wrappers carry the calls: `vrun_k`, `vrun_ix` and
  `vreduce_seg`. Nothing in the integration calls them yet.
- Every record now has a `# route` line per family and a `# identity`
  line.
- `--door-selftest` and `tools/check_door.py` hold the routes to each
  other, in all three suite tiers.
- The trailer's timing keys print from one table.
- `hw/bench-workload.py` writes `decompose.csv` with every trailer key.
- The gate cache keys on every `CFT_REBOUND_` variable.

### The gate, and the controls that show it can fail

At `733b067`, `--door-selftest` holds 131 cases at each format against
a reference computed outside the wrappers, straight from libcft.
- **Each case runs both routes:** 98 scalar (two with the result
  aliasing the scalar), 15 indexed, and 18 segmented, six of those a
  two-operand dot product.
- **Every wrapper call's accounting is checked:** one operation in
  `calls=`, and exactly the counter meant.
- **The data checks itself** at every format:
  - every nonzero value at binary128 and binary256 carries bits that
    binary64 cannot hold, by a round trip in libcft;
  - every one of those bit positions is set in some value, with the
    positions measured by widening 1 + 2^-52;
  - zeros of both signs occur.

  Measured: 9,423,979 of 9,423,979 values at binary128 and 9,425,290 of
  9,425,290 at binary256. All 60 and all 184 positions were set. There
  were about 314,000 zeros of each sign at every format.

`tools/check_door.py` passes with 56 `ok` and 0 FAIL, on Windows and on
Linux. Besides the cases it checks:
- every refusal by name, on both routes where a route is involved;
- that hidden bits count as not compared;
- that nothing compared is a failure;
- that a forced-off run gives the default's record.

That last comparison is vacuous until a parcel routes a call, and the
gate says so.

**Twenty-one sabotaged copies of the engine** were built from copies of
the source; the tree was never patched. Every one made the gate fail,
at the check meant for it:
- the five the first version caught;
- the five V0 found that version passing, N1-N4 and N7/N8;
- the four its first re-check found passing, F4, F4+N8, F6 and F11;
- the two refusals that re-check added, each removed, F12 and F13;
- the three its second re-check found, F15 (at two formats), F16 and
  F18.

### No record line moved

- **The A/B.** The same 27 integrations as entry 38 ran through the
  pre-P0 binaries (`21579ee`'s, built and kept) and through each P0
  build: `971f93d`, `90e6f3f`, `d7d964b` and `733b067`.
  - Every time, 27 of 27 were identical on every non-comment line and
    every counter. The two ensembles were compared on all their lines,
    for the reason entry 38 gives.
  - The control saw 4 of 5 lines differ every time.
  - The 27 assembled program images were byte-identical every time.
- **The suite.** `make check` at `971f93d` gave rc 0, 15 legs, 250
  `ok` lines and 0 FAIL. At `733b067` it gave rc 0 in 2,828 s, 15 legs,
  278 `ok` lines and 0 FAIL. Both had the same two skip lines as entry
  38's run.
  - The six wide-state digests are identical to entry 38's at both
    commits.
  - The only comment lines that differ from the pre-P0 records are the
    `# identity` line, the `# route` lines and `t_reduce` in the
    trailer.

### The verifier, three times

V0 is a separate agent with its own worktree and builds. Its entries
are in the round's ledger, which is to be archived as
`archive/round3-ledger.zip` at the round's end.

**On `971f93d`.** It confirmed that no record line moved, across 51
integrations of its own. It also confirmed the library path's refusal
and the DLL fix. It sent P0 back on what P0 hands the parcels:
- **Five sabotages the gate passed.** One of them, a probe that reads
  no bit, passed with nothing compared.
- **Two wrappers answering differently by route on aliased calls.**
  libcft's software backend reads a scalar again after writing the
  result (cft-fp256 `host/src/device.c`, 2068-2086). Whether that is
  libcft's intended contract is a question for cft-fp256.
- **Three false claims.** The worst was "the largest shapes the
  integrator presents". An ensemble has no width cap: `CFT_MAX_BODIES`
  bounds the bodies of one system, not the number of systems.
- **Smaller items**, such as the gate cache's key, `staged_bytes`
  depending on the route, and a misspelled variable being silently
  ignored.

**On `90e6f3f`.** It confirmed every fix. It found four of its ten new
sabotages passing, one more route split in `vrun_ix`, and the Windows
limit below.

**On `d7d964b`.** It confirmed every fix, and found that no new refusal
rejects a call the engine makes. It found three ways to fool the data's
check. One of them was a gap in the shipped data: binary128's bits 56-59
were never set.

`733b067` closes those three, **after V0's last cut**. The lead checked
it with V0's own three sabotages, which all fail. V0 has not re-checked
it.

### Stated limits, not closed

- **Arity.** `vrun_ix` does not refuse a table on an operand the opcode
  does not read, because libcft exports no arity. On that call the new
  route refuses by name, and today's route ignores the operand, so its
  result is the opcode's own.
- **Prefix typos.** A typo inside the `CFT_REBOUND_ROUTE_` or
  `CFT_REBOUND_HIDE` prefix makes an unrelated variable, which nothing
  reads.
- **Windows and Python.** In a Python process on Windows, the library
  half's `getenv()` and `_environ` are msvcrt's, and that runtime's copy
  of the environment is taken when the process starts. A variable set
  through `os.environ` is never seen. That holds for the door's
  variables and for `CFT_REBOUND_ARTIFACT`. V0 measured it, and
  docs/PYTHON.md, "The Windows rule", says so.
- **Before P0, found during it.** On the drop-in's refusal path, the
  sentence prints twice and `steps_done` becomes 1. The pre-P0 build did
  the same for a missing artifact. From Python, REBOUND's message store
  cuts the sentence at 254 characters. ROUND3.md files both for the next
  round.

### What this does not establish

- **The card.** No artifact was opened.
  - Nothing here shows that 0.18 opens this repository's f128 image.
  - A tile's own route for any family has not run. On the software
    backend each new route is the same primitive today's route uses,
    so these cases test the wrappers, not a tile.
  - The card day's gravity split (ROADMAP item 4) is still owed.
- **Use.** Nothing in the integration calls the door. The forced-off
  comparison is vacuous, and the gate says so.
- **amd-arc-box's Linux**, as against WSL's.
- **UCRT64**, which docs/PYTHON.md reasons would give the Windows
  library one heap and one environment. It has not been built.

## 40. Round 3's P4: `predict_next_step` as one program run on `--engine program`, in the elementwise route's bits at every format; and the plan's compensated additions, which were never there

**2026-11-02.** Two hosts: this desktop, as in entry 39, on the software
backend; and amd-arc-box's U50, for the card row, with the rev8a
single-tile image. Pins: REBOUND at bdfda4bd and cft-fp256 at 922c2ad
(ABI 0.18).

P4 branched from `af2f400`. Its commits on `round3-p4` are `cae9916` and
the follow-up `0c7ad74`. Merge 1 is `00f89b2`, on `32dc962`, which
merges the card tooling. `0c7ad74` came in with merge 2, as `3159d92`.
The card ran `c3471bc`, which is `0c7ad74` merged with `round3-cardday`.

### What P4 is, in one paragraph

Under `--engine program`, `predict_next_step_vec()` now runs
`nextstep-<fmt>` and returns. That is one program run where there were
69 elementwise calls (34 `vmul`, 28 `vadd`, 7 `vsub`).
- `programs/nextstep-{fp64,fp128,fp256}.cfta` are checked in like the
  other 27: 119 instructions, 11 external constants, scratch 14 in and
  14 out, no deposits.
- Stream a is `rb`, each lane's own ratio. The scratch block holds
  e0..e6 and b0..b6 in, and the predicted e and b out, in place.
- The 21 binomial products use only 11 distinct values. Deduplicating
  them keeps the bank inside the port's 16-slot allocation.
- A new trailer key, `nextstep_calls`, counts these runs.
- The elementwise body is unchanged. `--engine loop` and the library
  half still run it: `ias15_engine_alloc()` sets `engine_program = 0`.

### What the plan had wrong

`predict_next_step` has no "sixteen compensated additions", and neither
has REBOUND's. The compensated additions are the end-of-step update in
`step_attempt()`: 17 `add_cs` a coordinate, 9 into x0/csx and 8 into
v0/csv. ROADMAP.md and ROUND3.md ran the two together, and
docs/HARDWARE.md:147 keeps them apart. P4 read this at `af2f400`, and
the lead checked it there. So the brief's negative control could not
exist, and the lead amended it to one plain addition or multiplication
dropped. The end-of-step update stays elementwise, filed for the next
round.

### The gate, and the controls that show it can fail

`tools/check_program_engine.py` compares `--engine program` with
`--engine loop`, both at `--arith fma`, at fp64, fp128 and fp256. The
brief's trap is that the step ratio is per system. A program that took
it from the bank would be right for one system and wrong for an
ensemble, and no existing case could see that. So the gate gained:
- **the "spread" Kepler ensemble**, E=5, adaptive, dt0 = 0.6, unscaled
  in every tier, which asserts that the systems' ratios differ by more
  than 2^-40 relative and that the first attempt is mixed;
- **on every case**, the step counters and `flags_seen` held to the
  loop's, the ratio stated and asserted (exactly 1 for fixed steps, not
  1 somewhere for adaptive), and `nextstep_calls`: 0 on the loop engine,
  and on the program engine between the most steps any system accepted
  and the total attempts;
- **from `0c7ad74`**, a Kepler e = 0.99 case at epsilon 1e-2 whose 11
  rejections all predict from er and br (V4's finding, below).

The ensemble case reads the same at every format:

```
ok   kepler spread E=5 adaptive from dt0=0.6 (mixed rejection): 20 record lines identical, 5 systems x 4 samples
     counters and flags the loop's; nextstep ran 121 times; the systems' ratios differ at 3 of 3 sampled steps (1.078295 to 1.112737 at the last); a mixed attempt: steps_rejected=[1, 1, 1, 1, 0]
```

| sabotage | property | with it |
|---|---|---|
| N1: level 0's `mul en, q1, s` replaced by a copy, all formats | a dropped multiplication changes the record | FAIL on Kepler adaptive and the ensemble; inert on Kepler fixed and outer fixed, where every ratio is exactly 1 |
| N1a: level 0's `+ be` dropped, fp64 | a dropped addition changes the record | FAIL on Kepler adaptive and the ensemble; on Kepler fixed only `mean_pc_iterations` catches it, 2.455 against 2.460; inert on outer fixed |
| N2: the ratio taken from the bank, lane 0's, all formats | the per-system ratio | the three single-system cases pass; the ensemble FAILs at system 1 |
| N3: the program never runs, fp64 | the program route ran | records identical in all 4 cases; caught only by `nextstep_calls` (0, must be 200 to 200) |
| N4: the bank's order reversed, fp64 | the host's bank order is the generator's | FAIL; on outer fixed only through the counters, 2.283 against 3.933 |
| N5: the copy-in reads `eo`/`bo`, the state's e and b, not `_e`/`_b` | a prediction from er and br | the rejection case FAILs at all three formats, and the older four stay `ok` |

Each passes without its sabotage. N2, N3 and N4 each pass the plan's
existing comparison somewhere, which is why the gate holds the counters
and `nextstep_calls` as well.

### No record line moved

- **P4's A/B against `af2f400`.** 12 problems on both engines at all
  three formats, run twice: 72 of 72 integrations identical in every
  record line, every step counter and `flags_seen`.
- **The accounting is the change.** On `--engine program`, `calls=`
  falls by exactly 68 and `program_calls=` rises by exactly 1 for each
  run of the program. On `--engine loop` neither moves.
- **Merge 1's suite.** The full `make check`, from clean at `00f89b2`,
  read from its log:
  - 15 legs, all in full mode, **281 `ok`, 0 FAIL**;
  - the only skip is the deliberate `body_count`;
  - every `RESULT: PASS`, `check_dropin` every case passed, and
    `gate_real` passing at all three formats;
  - the six wide-state digests unchanged: fp64 `e08541bc0af9b3ec` and
    `56bacf3bf63f19a7`; fp128 `5d7cf4588fad050e` and
    `cd16e6dc29aa0589`; fp256 `a1052750e206bebe` and
    `ba4d447c8e58bd20`;
  - `check_door`'s forced-off comparison still "vacuous", as it must be
    with no routed call merged.
- **`0c7ad74` is not in that suite.** Its diff is one file,
  `tools/check_program_engine.py`, so no record can move. It went into
  merge 2's tree (entry 41).

### The verifier

V4 has its own worktree and builds, as V0 did (entry 39). It found no
blocking defect: READY on `cae9916`, and READY on `0c7ad74`.

**On `cae9916`.**
- **Records against `af2f400`.** 184 pairs, 368 integrations over 32
  configurations, on both engines and at `--arith rebound`. They
  include mid-run rejections, a `--dt-file` with varying steps so that q
  is not 1 in a fixed run, and ensembles of 3 to 64 systems. 180 pairs
  are identical. The other 4, the two loosest ensembles at fp64, stop on
  the flag abort (`mul raised 0x18`) on both builds and both engines,
  with the same rc, stderr and record lines.
- **The program against the elementwise route, directly.** An
  instrumented engine ran the elementwise route after every nextstep
  run, on the same inputs, and compared all 14 x N3 outputs. Over 80
  integrations: 34,338 program runs, 0 elements different. N1's programs
  give 6,408 differing elements over 581 runs, and N4's 12,960 over 60
  runs of outer fixed, whose records stay identical.
- **P4's five controls** reproduce with P4's numbers.
- **13 sabotages of its own, at fp64: 12 fail the gate.** One of the 12,
  the copy-back skipping the last lane, is caught only on outer fixed,
  through `mean_pc_iterations`, 2.283 against 3.067. Both Kepler
  problems are planar, so their last lane is a z whose b is 0. The
  thirteenth, V11, passes (below).
- **The finding.** P4 had stated that no case predicts after a mid-run
  rejection, and that none of seven longer runs reached one. V4 showed
  the path is reachable: Kepler e = 0.99, dt0 1e-4, epsilon 1e-4 gives
  93 rejections in 1,500 steps, and `nextstep_calls` is 1,593 = 1,500 +
  93. Those runs are among the identical pairs.

Both (b) findings went back to P4 as one commit, `0c7ad74`: the
rejection case, and V11 stated as a limit in the gate's header.

**On `0c7ad74`.** Rebuilt from clean, the gate passes at all three
formats. N5 fails the new case and nothing else. The case's own
assertion fails when its premise is removed. No earlier sabotage's
verdict changes.

### On the card

The leg `p4-program` ran `tools/check_program_engine.py` at `c3471bc`,
with `CFT_REBOUND_ARTIFACT` set to the rev8a single: rc 0, **15 ok, 0
FAIL, `RESULT: PASS`**, in 56 minutes. That is `--engine program`
against `--engine loop` at all three formats on a tile, with the
nextstep, corrector and predictor programs and the rejection case.

`c3471bc` is P4's code and the card tooling only, not a commit on main.
Believed by the lead: the merge changes only shared one-liners.

### Stated limits, not closed

- **Program engine only.** A drop-in or Python user always gets the
  elementwise route, and no library-path gate can exercise this program.
- **Fixed-step cases cannot see a fault in q's powers.** The ratio is
  exactly 1 there, and the gate prints that on each fixed case.
- **The "it ran" check is exact only on the single-system cases.** On
  the ensemble it is a range, 120 to 604, holding 121. An ensemble that
  fell back to the elementwise route on some predicting attempts would
  pass, with the same bits.
- **`outer fixed dt=40`'s records cannot see a wrong corrector start.**
  With the bank reversed all 160 record values stay identical. Believed
  by P4: on that problem the corrector reaches its fixed point bit for
  bit from any reasonable start.
- **V11: the nextstep run's flags are not held.** Dropping them passes
  the gate. Believed by V4: inert by construction, since the program
  issues the elementwise route's 69 operations, and on every problem V4
  ran their flags are inexact alone, which other calls raise every step.
- **`staged_bytes` over-counts this program.** Nothing compares it.

### What this does not establish

- **Speed.** P4 made no speed claim, and the card leg compares bits.
- **The merged tree on a card,** as against `c3471bc`.

## 41. Round 3's P2: the corrector's exit test on the device, and a recorded schedule that replays without it; sent back twice, and changed in no record line

**2026-11-02.** Two hosts, as in entry 40. Pins: REBOUND at bdfda4bd and
cft-fp256 at 922c2ad (ABI 0.18). P2 branched from `af2f400`. Its commits
on `round3-p2` are `0b5e17e`, then `faa9f22` and `53ffd5f`, then
`8006820`, then `9539f1d`. Merge 2 is `42b4863`, which merges `8006820`,
after P4's `3159d92` and with the lead's `ce593f7` on top. `9539f1d`
came in with merge 3, as `ea8d690`. The card ran `aa3aa7e`, which is
`53ffd5f` merged with `round3-cardday`.

### What P2 is, in one paragraph

cft-fp256 measured removing the corrector's exit test at 4-9% of a card
step. P2 builds two mechanisms aimed at that share.
- **The recorded schedule.** `--record-iters FILE` writes the pass
  count of each step attempt. `--replay-iters FILE` runs exactly those
  counts, and skips `pc_error()` and the exit test. A schedule carries
  a run key, and one from any other run is refused by name.
- **The device test.** Each maximum in `pc_error()` is now cft-fp256's
  `normalabs-<fmt>` program followed by `CFT_MAXALL` through
  `vreduce_seg()`: `cft_reduce_seg` where the device publishes CAPS2[8],
  and today's host loop as the fallback. All four modes take it. Each
  maximum costs two library calls on either route, and counts as one
  `reduce_seg` operation.
- **The images.** At `af2f400` the library half had no program loader.
  `build/normalabs_images.h` is assembled from the pinned tree by the
  pinned `cft-asm`, and both halves compile it in.
- **`--corrector-selftest`** holds `pc_error()` on both routes against
  a reference built from REBOUND's definition with libcft's own calls:
  11 crafted systems and 53 random ones.

### The gate, and the controls that show it can fail

`tools/check_corrector.py` is a new leg, and `tools/cases_corrector.c`
its library-path case. The brief's trap is a replay that reproduces the
state but not the counters, or keys on the wrong unit. At fp256, from
`0b5e17e`'s full suite:

```
ok   it computed no test: reduce_seg carried nothing; divsqrt_calls 87 fewer (87 passes x 1); calls 495 fewer (wanted 495: 3 a pass, and two compares at each of 117 tests); flags 0x11 within 0x11
ok   fp256: the replays crossed 4 runs with a rejected step, 4 with a capped attempt (each counted as iterations_max_exceeded), 1 ensembles whose systems left at different passes, 1 with some members rejecting and some not
```

Each control ran through the full gate at `0b5e17e`, which gives 243
`ok` and 0 FAIL without one:

| sabotage | property | FAIL lines |
|---|---|---|
| S0: the replay ignores its counts and runs the test | a replay runs the recorded counts | 59 (38 + 10 + 11, one count per format) |
| S1a, S1b: masked lanes become -0, or -inf | no normal coordinate gives +0 | 6 each, at the self-test and the image digest |
| S2: the exponent-zero test dropped | subnormals are masked | 6: `0x1.3028p-1027 where the reference is 0x0p+0` |
| S3, S3b: the 1/1 removed from the device branch, or from INDIVIDUAL's shared lane | a system out of the corrector is given 1/1 | 3 each: `nan where the reference is 0x1p+0` |
| S4: the prototype's bug | a replay counts a recorded cap | 14 |
| S5: keyed on accepted steps | the schedule is per step attempt | 47; V2's copy, keyed on any system accepting, 39 |
| S6: the `normalabs` run counts as a program call | `program_calls` does not move | 50 |

S1a to S3b each left all 237-240 integration-level checks passing.
Records, counters and replays cannot see those properties, which is why
the self-test exists.

**ROUND3.md's replay control, one count moved by one, is inert on most
attempts.** 0 of 122 such moves changed a record line, over the first 8
attempts of 8 configurations. It does bite: 4 of 60 moves at fp128, at
attempts 25 and 28. So P2's first reading, that an attempt's last pass
changes no record bit, was false as a general statement, and the lead
withdrew it. The lead accepted P2's replacements: the first attempt cut
to one pass, a capped count refused by name, and the counters moving by
exactly the change.

**`check_door` holds `program_calls` from merge 2.** S6 passed it until
the lead's `ce593f7`, and then fails, "(DIFFER: program_calls 344
against 0)".

### No record line moved

- **P2's A/B against `af2f400`.** 87 of 87 integrations identical at
  all three formats: every record line, every counter, and each
  system's trailer. Its control saw 2 of 3 lines differ.
- **`check_door`'s forced-off comparison** no longer says vacuous, and
  passes at all three formats.
- **The library path.** All 5 corrector cases in `check_dropin` were
  bit for bit REBOUND's ias15, on the device route with no fallback.
- **P2's own suite** at `0b5e17e`: 16 legs, 531 `ok`, 0 FAIL, the six
  `gate_real` digests equal to P0's.
- **Merge 2's suite.** The full `make check`, from clean at `ce593f7`,
  with P4's `0c7ad74` and P2's `8006820`, read from its log:
  - rc 0, 16 legs, all in full mode, with no cache skips;
  - **571 `ok`, 0 FAIL**, no DIFFER;
  - the six wide-state digests unchanged;
  - `check_door`'s forced-off comparison at every format: "all 9
    counters identical".
- **`9539f1d` is not in that suite.** It changes one file,
  `tools/check_corrector.py`, and went into merge 3's (entry 42).

### The verifier, four times

**On `0b5e17e`: SENT BACK, 2.** V2 held 258 pairs against `af2f400`, all
identical: 129 on the default route and 129 with
`CFT_REBOUND_ROUTE_REDUCE_SEG=off`. Its comparator's controls moved 38
of 43 and 7 of 43 fp64 pairs. P2's seven controls reproduced.
- **(a) M2.** A replay that keeps every system in the corrector to the
  attempt's last count, still counting each system's passes exactly,
  passed at every format: passes past a system's exit land on its fixed
  point. V2's control cuts one system that is not the longest to 1 pass.
  Over 30 cuts the real replay changed only that system's lines in 28,
  was refused by name in 1 and was inert in 1. M2 changed nothing in 30
  of 30. Believed by V2: M2 gives a wrong record where extra passes miss
  a fixed point. It found no ensemble where that shows.
- **(a) M4.** A library half silently on today's route passed
  `check_dropin --light`, reading `off by=device`, a state a correct
  door cannot reach on software with no door variable.
- **(b) M3, M3s.** A run key without epsilon or the state hash passed.

`faa9f22` holds each system's own exit and all 20 run-key fields, and
`53ffd5f` makes the library case derive the state the door must be in.

**On `53ffd5f`: SENT BACK, 1.** M2, M3, M3s and M4 fail where meant.
- **(a) M2b.** A system at its count rides one more pass, uncounted. It
  passed at every format, since a system frozen after 2 passes still
  moves. An exact reference sees it: the member alone with the same cut.
- **(b) M3v, M3m, M3dt.** `state` is one hash, and a hash without the
  velocities, masses or dt-file passed.

`8006820` closes both, in `tools/check_corrector.py`.

**On `8006820`: READY,** with one should-fix, Rboth: a replay that
refuses every uncapped count of 1, by name, passed at every format,
because all three cut controls took their refusal branch. `9539f1d`
requires a cut replay to run.

**On `9539f1d`: READY.** The clean gate gives 277 `ok`, 0 FAIL at all
three formats, and every cut replay runs. Rboth gives 268 `ok`, 9 FAIL.
M2b gives 274 `ok`, 3 FAIL: 12 of 13, 2 of 3 and 2 of 3 lines differ
from the member reference, where the clean tree differs in none. M2
gives 271 `ok`, 6 FAIL. Since `0b5e17e` only `tools/cases_corrector.c`
and `tools/check_corrector.py` have changed.

### On the card

The leg `p2-replay` ran at `aa3aa7e`, on the rev8a single, each run
binary128, `--engine program --arith fma`. `53ffd5f` is not what merged,
but the engine has not changed since `0b5e17e` (V2, above).

| problem | card recording | its replay | change | replay against recording | software recording replayed on the card |
|---|---|---|---|---|---|
| kepler, 200 steps | 55.64 s | 53.95 s | -3.0% | data identical; counters equal | data identical |
| outer, 100 steps | 34.48 s | 32.24 s | -6.5% | data identical; counters equal | data identical |
| nbody64, 20 steps | 25.77 s | 25.45 s | -1.2% | data identical; counters equal | data identical |
| nbody256, 5 steps | 50.46 s | 50.92 s | +0.9% | data identical; counters equal | data identical |

- Measured: the replay's byte-identical record is the proof that
  skipping the exit test is exact, on a tile. A schedule recorded on the
  software backend replays on the card to the card's own recording.
- Believed: the saving, 1 to 6.5% on the small problems and nothing
  measurable at 256 bodies, is the exit test's share of a step. One run
  per row; timing noise is not measured.
- On 2026-11-02 the leg `p2-corrector` was skipped: its estimate, 60
  minutes, would have ended past the card's window.

On 2026-11-06 the leg `p2-corrector` ran at `85fd0aa`, `e47b04e`'s
engine, on the rev8a single, its runs one at a time:
`tools/check_corrector.py` rc 0, **277 ok, 0 FAIL, `RESULT: PASS`**, in
3,467 s. That is the count software gives. So the device's exit test is
held against the host loop on every configuration, on a tile, with the
replays, the cuts and the run key's refusals.

The leg `p2-bench` ran at `e90d5df`. The device's test gave software's
bits and was slower than the host's loop on all 12 rows; entry 45 has
them. So `eb6976a` put the route on the opt-in
list, and from there a default run takes the host's loop even where a
device publishes the bit. "The route" in this entry is now the forced-on
one.

### Stated limits, not closed

- **A replay's `flags_seen` can be narrower than its recording's,** by
  flags only the test raises. The gate asserts within.
- **A device that publishes CAPS2[8] but refuses the image is refused
  at open, not fallen back.** The library half's path has not been run.
- **INDIVIDUAL's data movement stays host work on both routes.**
- **A schedule names its run but not the backend or the library
  build.** That is how a software recording replays on a card.
- **The wide formats.** The gate's route comparison there leaves modes
  1 and 3 to the self-test. V2's A/B covers every mode.
- **The member reference is held at one cut per format.**
- **M1as2, a count of 1 replayed as 2 passes, passes.** No recording
  holds a 1: none in the gate's, and none in six more runs.

### What this does not establish

- **Why the device's test is slower.** Entry 45 measures that it is.
  Nothing here separates the launch from the reduction.
- **What either mechanism saves.** The 4-9% is cft-fp256's, for
  removing the test outright. This port's replay rows are one run each.

## 42. Round 3's P3: a checkpoint resumes bit for bit at binary128 and binary256, adaptive runs included, and the archive says who made it; sent back once; and the seam test

**2026-11-02, this desktop, as in entry 39, on the software backend
only.** Pins: REBOUND at bdfda4bd and cft-fp256 at 922c2ad (ABI 0.18).
P3 branched from `af2f400`. Its commits on `round3-p3` are `d492c32`,
the lead's accessor, a cherry-pick of `9b1744ea`; `a8d7a63`, the gate
alone; `3861a7b`, the fix; `dd51e21` and `d8ae085`, after V3; and
`c957eb2`, after V3's READY. Merge 3 is `4291740`. It holds `ea8d690`
(P2's `9539f1d`), `9103ab6` (P3 at `d8ae085`), `b231398` (the lead's
seam test), and `c957eb2`.

### What P3 is, in one paragraph

One on-disk change, with three items.
- **The wide clock and steps.** The archive now carries the engine's
  `t`, `t_hi`, `t_lo`, `dt` and `dt_last`, one wide element each.
- **The restore field.** A new value, `CFT_PROV_EXACT_CLOCK_STEPS64`,
  says "restored exactly, its clock and steps from binary64".
- **Identity.** The archive carries the record's `# identity` words:
  `libcft_build=commit=922c2ade78fad995a0871bf07cebf8b9fbcb2b9f,tracked=clean,untracked=none image=none`.

A load adopts each wide value only where its binary64 view still equals
REBOUND's field, and promotes the field otherwise, because
`reb_simulation_integrate()` rewrites `r->dt` and zeroes `r->dt_last_done`.

### What the plan had wrong

- **The restore field exists,** as `cft_provenance`, from `7fbbd1f`
  (2026-09-11). The name stays: a rename costs every archive since.
- **The clock alone was not enough above binary64.** At `af2f400`,
  binary128, epsilon 1e-9, 20 steps, a save and 20 more, 4 of 18 dump
  lines differ, because `dt` and `dt_last` are wide too.
- **The 2-of-40 setup was never recorded.** The case is now gate_real's
  four bodies, fixed dt 0.01, a checkpoint after s = 1..40 steps, then
  40 more after each. 180 steps split at 1..40 gives 0 of 40.
- **Entry 38's cross-pin check compared dumps,** which carry no
  identity, so nothing there needed excluding.

### The gate, watched failing and then passing

`tests/gate_real.c` now compares the engine's clock and steps as bytes.
`wide_digest()` is unchanged, so the six recorded digests stay
comparable. Before the fix, at `a8d7a63`:

```
fp64   splits   FAIL  35 of 40 resumes DIFFER (r->t differs on 0 [], the clock or steps on 35, other lines on 0)
fp128  splits   FAIL  34 of 40 resumes DIFFER (r->t differs on 2 [29, 30], the clock or steps on 34, other lines on 0)
fp256  splits   FAIL  34 of 40 resumes DIFFER (r->t differs on 2 [29, 30], the clock or steps on 34, other lines on 0)
fp128  adaptive FAIL  the continued run DIFFERS - line 5: straight 'p1 x  0x1.ef3d04630a0acp-1 0x1.1181227b2e8f4p-2 0x0p+0', resumed '... 0x1.1181227b2e8f5p-2 ...'
```

At binary128 and binary256 only t and t_hi differ, below binary64's
reach: at split 30, binary128, `3ffe6666666666666880000000000000`
against `...6600000000000000`. At binary64 only t_lo differs, because a
resume zeroes it.

With the fix, 40 of 40 resumes are identical at every format, and the
adaptive case gives 17 values identical. **The named control,**
`CFT_REBOUND_NO_CLOCK`, the clock's load removed, brings back 35, 34 and
34 of 40, with r->t on [29, 30] at fp128 and fp256. `CFT_REBOUND_NO_STEP`
gives binary128's wide digest as `85173733207ae8fd` straight against
`7bebb83a206f476b` restarted. Nine more sabotages were all caught.

### An older archive, and an older build

- **An `af2f400` archive here:** 6 of 6 load, and say "restored exactly,
  its clock and steps from binary64".
- **One snapshot of this build's archive in `af2f400`, no function
  pointers set:** 9 of 9 load and say "restored exactly", and continue
  as `af2f400` continues its own. REBOUND prints `Error! Could not find
  field descriptor for name.`, because the old reader stops at `cft_t`.
- **Past that, V3 measured three behaviours, and P3 reproduced them.**
  `af2f400` appending to this build's archive writes a snapshot no build
  can read. A snapshot this build appended after a particle removal is
  refused by `af2f400`. REBOUND's function-pointer warning is lost when
  `af2f400` reads this build's archive with a force or merge set.

`tools/check_archive.py` states these, with "never append to these
archives with a pre-round-3 build". Believed by P3 and V3: any appended
field costs a pre-existing build these three.

### The verifier, three times

**On `3861a7b`, the lead's `d492c32` included: SENT BACK, 2.** Nothing
moved: 312 library-half integrations against `af2f400`, and 9 of 9
records from the standalone program. `d492c32` is tree-identical to
`9b1744ea`.
- **(a) A1: `dt_last`'s restore was ungated.** S4 (the `r->dt_last_done`
  guard dropped), S4b (`dt_last` not restored) and S22
  (`put_dt_last_wide()` a no-op) passed every leg. `dt_last` is read
  only after a rejected attempt, and no case rejected the first attempt
  after a load. S4 and S4b give a resume that is not REBOUND's own,
  where P3's build matches it at fp64.
- **(a) A2: four shipped sentences understated the older-build limit,**
  by the three behaviours above.
- **(b) B1:** no case triggered the clock's "reached memory" refusal (S11).
- **(b) B2:** the 2-of-40 count was printed, not asserted.

`dd51e21` adds sequence 6, a rejected first attempt after a load;
sequence 7, the "reached memory" refusals; and the asserted sweep.
`d8ae085` states the limit by its three behaviours.

**On `d8ae085`: READY.** Its library objects are `3861a7b`'s.
- S4 and S4b fail at all three formats. S22 fails at fp128 only, "wide
  state 5168967c5e8aa99c against df037a7d75397f63": at fp64 `dt_last` is
  binary64, and at fp256 the corrector re-converges. The sequence had to
  lift body 2 out of gate_real's plane (vz 0.001), or S22 passed at
  fp128 too. Believed by P3: a control that only rounds a wide value can
  be inert where the corrector converges past it.
- S11 fails at all three formats. With AFTER = 1, 39 and 41 the sweep
  fails "NOT THE 2-OF-40 CASE".
- Two should-fixes: two phrasings of the limit, and sequence 6's
  rejection counter, which S5 (adopt's `r->dt` guard dropped) fooled.

**On `c957eb2`: READY.** S5 now fails that counter at every format: "the
resumed run took the step with no request - it never saw the request".
The six library objects are byte-identical to `d8ae085`'s, and the
rewritten sentences match V3's A/B, 39 of 39 pairs.

### The seam: `tools/check_seams.py` and `gate_real --routes`

The plan's P2 x P3 seam, "a replay across a checkpoint", cannot exist:
the replay is the program's and the checkpoint the drop-in's. So the
seam saves with `reduce_seg` forced one way, resumes with it forced the
other, and requires the straight run's dump.
- `tools/check_seams.py` imports `check_checkpoint`'s sequences. Its
  control is `CFT_REBOUND_NO_ADOPT`, and it is a leg from merge 3.
- **It needed a physics counter.** Comparing dumps alone, a crossing
  whose forced route carried nothing would compare two runs of today's
  route. `gate_real --routes` prints, on stderr after each phase, each
  family's route from `ias15_engine_route()`. `check_seams` asserts the
  forced route was taken, by=forced, and carried the run: calls and no
  fallback when on, fallbacks and no call when off.
- **Measured on the merged tree:** 12 of 12 crossings at all three
  formats, each carried on both sides. At binary256 the plain sequence
  carried 1,214 calls at the save and 2,400 fallbacks at the resume.
  Every resume gave the straight run's values, and every control differed.
- **What it cannot see.** `reduce_seg` keeps no state from one step to
  the next, so the route has nothing for an archive to drop. The seam
  can catch only a restore that depends on the route that wrote it. No
  such sabotage was built, so that it would show is believed, not measured.

### Merge 3's `make check`

From clean at `4291740`, read from its log:
- rc 0, 17 legs, all `--mode full`, with no cache skips;
- **584 `ok`, 0 FAIL**, no DIFFER and no traceback;
- 9 `RESULT: PASS` lines, and `gate_real: PASS` three times;
- the six digests unchanged: `e08541bc0af9b3ec`/`56bacf3bf63f19a7`,
  `5d7cf4588fad050e`/`cd16e6dc29aa0589`,
  `a1052750e206bebe`/`ba4d447c8e58bd20`;
- `check_door` at every format: "all 9 counters identical";
- `check_seams`: "12 route crossings carried by the route on both sides,
  of 12".

### On the card

On 2026-11-02 the leg `p3-identity` was never chained, because P3 was
sent back at 12:45. The card branch of the identity check, that an
xclbin must yield a 64-digit image digest, was written and had not run.

On 2026-11-06 the leg ran at `85fd0aa`, `e47b04e`'s engine, on the rev8a
single tile, in 3,111 s:
- `gate_real: PASS` at binary64, binary128 and binary256, `identity_ok()`'s
  card branch among what it holds;
- a binary128 save and resume on the card reads `written on
  "libcft_build=commit=922c2ade78fad995a0871bf07cebf8b9fbcb2b9f,tracked=clean,untracked=none
  image_sha256=0639248e9c66f779e7c126f4449062af646d3d7e8abd0173926812790bc80450
  image_bytes=40076366 image_version=0xb00`;
- `sha256sum` of the xclbin gives the same digest.

So a card-written archive names the image it ran on, and the name is the
image file's own digest.

### Stated limits, not closed

- **A particle add or remove truncates the wide clock and steps to
  binary64,** in a straight run and a resumed one alike, so no gate can
  see it. Believed by P3: at `accurate = 1` with an adaptive step above
  binary64, it costs `dt`'s tail at every merge. Filed for the next round.
- **The identity is formatted in two places,** the program's `main()`
  and the shim. Their equality was measured once, not gated.
- **S10, the clock count check removed, passes:** the memory check
  refuses the same files.

### What this does not establish

- **A second image.** The digest is measured on rev8a only. A software
  archive records `image=none`.
- **Python.** Nothing in python/ loads a cft archive.

## 43. Round 3's P1: gravity as a handful of program runs a force evaluation, on two switches of the door, both opt-in because each made a card step slower at 64 to 512 bodies; sent back once, and changed in no record line

**2026-11-02.** Two hosts, as in entry 40. Pins: REBOUND at bdfda4bd and
cft-fp256 at 922c2ad (ABI 0.18). P1 branched from `af2f400`, and its
commits on `round3-p1` are `d4e710b`, `3a00e18` after V1, and `0a0d96d`,
the opt-in. Merge 4 is `e47b04e`: `5c525b4` merges P1 at `3a00e18`,
`2b497a6` is the lead's seams, `e77486e` merges `0a0d96d`, and `e47b04e`
is the lead's fix to `hw/cardday-round3.sh`. The card ran `afc83d8`,
which is `d4e710b` merged with `round3-cardday`.

### What P1 is, in one paragraph

A force evaluation was 22 elementwise calls plus one add per partner. It
is now `gravity-pair`, `cft_sqrt`, `gravity-cube`, `cft_div` and
`gravity-force`, then one fold run per 256 partners.
- **Two switches.** `indexed` is the accumulate fold. `pair_gather`, a
  new row, `CFT_REBOUND_ROUTE_PAIR_GATHER`, is the pair program's six
  inputs. Both read `CFT_SEQ_FEAT_INDEXED`, so the hide hook hides both.
- **Off is not the old gravity.** The off route of each switch is the
  host's gather into the same program, run dense. So `calls=` and
  `program_calls=` are equal on all four combinations.
- **Every mode, both halves.** Neither route needs `--engine program`.
- **The fold.** libcft only gathers. The order and the +0 start are the
  program's: its slot order, from a zero stream.

### On the card, and why both switches are opt-in

The leg `p1-gravity` ran at `afc83d8` on the rev8a single tile.
`check_gravity` there: rc 0, **23 ok, 0 FAIL**. `bench-workload.py` ran
once per route combination, each route forced, and every card row of
all four is **identical to software**: 12, 12, 12 and 6 rows. The both-on
bench dropped the three small problems to fit the card's window. Card
seconds, and the change against both off; one run per row, so noise is
not measured:

| problem | format | both off, s | indexed only | pair_gather only | both on (what auto then picked on rev8a) |
|---|---|---|---|---|---|
| kepler | fp64 | 22.44 | -5.6% | +0.9% | not run |
| kepler | fp128 | 44.06 | -7.4% | -2.2% | not run |
| outer | fp64 | 11.08 | +3.7% | -0.3% | not run |
| outer | fp128 | 22.72 | +3.7% | -0.4% | not run |
| pythagorean | fp64 | 10.92 | +1.3% | +5.6% | not run |
| pythagorean | fp128 | 19.39 | +4.5% | -0.5% | not run |
| nbody64 | fp64 | 6.57 | +22.2% | +22.7% | +47.2% |
| nbody64 | fp128 | 15.93 | +25.2% | +24.5% | +52.1% |
| nbody256 | fp64 | 23.14 | +40.0% | +37.7% | +73.7% |
| nbody256 | fp128 | 47.45 | +32.8% | +31.7% | +62.2% |
| nbody512 | fp64 | 64.86 | +29.6% | +28.2% | +60.5% |
| nbody512 | fp128 | 118.29 | +30.3% | +29.8% | +59.3% |

Measured: each gravity route costs 22-40% at 64 to 512 bodies, and with
both on the costs about add, to 47-74%. On two and three bodies, each is
within -7% to +6%.

`d4e710b` is not what merged. V1 measured that engine the same as
`3a00e18` (below), and measured no integration function changed at
`0a0d96d`. On 2026-11-02 `0a0d96d` had not run on a card.

On a card that published the bit, auto turned both routes on. So P1
merged as it stood would have made a default card run at 64 to 512
bodies 47-74% slower, with identical bits. Given that table, the owner
chose opt-in first. In `0a0d96d`, `route_def[]` gains `opt_in`, set for
`indexed` and `pair_gather`. Under auto both stay off even where the bit
is published, and say `by=opt_in`. `on` and `off` are unchanged, and
`scalar` and `reduce_seg` keep auto as it was.

On 2026-11-06 the leg `door-rev8a` ran at `e90d5df`, `e47b04e`'s engine,
on the rev8a single tile: `check_door` rc 0, **77 ok, 0 FAIL**, in 766 s,
"164 of 164 cases identical to the reference on both routes" at each
format. The opt-in branches, right by reading until then, are measured.
At binary64: "the default leaves indexed and pair_gather off by=opt_in
(fallbacks 1254/1254), forced on they carry calls 1254/1254 and no
fallback". Entry 45 has the same gate on the final tree.

**"Both off" is not the gravity from before P1.** The table above compares
P1's two device gathers with P1's host gather. The 22-call gravity
survives only as the self-test's reference. So the lead's word to the
owner on 2026-11-02, that opt-in keeps a default card step as fast as
before, undersold it and was not measured. Against `af2f400`, the round's
"before", `d4e710b`'s both-off rows are -16% to -35% up to 64 bodies, -8%
to -10% at 256, and +1.5% to +6.3% at 512. Entry 45 has the merged tree's.

### The gate, and the controls that show it can fail

`tools/check_gravity.py` runs every combination the switches allow, at
every format, with route lines asserted. `tools/cases_gravity.c` is the
library half's gate, through `ias15_engine_route()`, and asserts physics
counters: "a merge removed one: 4 particles became 3, on both sides".
`--door-selftest` holds 164 cases at each format.

| control | sabotage | with it |
|---|---|---|
| Order | body 0's 2nd and 3rd addends swapped in the fold table | outer's binary64 record differs: fixed 3 of 230 fields, adaptive 126 of 230, program engine 85 of 184; `check_equivalence` FAILs 3 outer cases, first at step 25 |
| Route | one index given to libcft points at the wrong partner, new route only, on either table | outer fixed: 148 of 230 fields differ when the route is on, 0 when off; `check_door` and `check_gravity` FAIL |
| Trap 1, the fold's zero | the fold started from -0 | the self-test's signed-zero case, 2 of 15 lanes |
| Trap 2, loud | tables never rebuilt | refused by name: "idx_scratch_in[1] = 9 ... is at or past the 6 elements" |
| Trap 2, silent | `set_active` skips the rebuild | `cases_pairs`: 155 of 198 configurations disagree with REBOUND |

Each passes without its sabotage. The order-sensitive problem is the
outer solar system. One force evaluation rarely shows one swapped pair,
so the order is held by the REBOUND gates over integrations.

### No record line moved

- **P1's A/B against `af2f400`.** 34 integrations at all three formats,
  each on the default route and the four forced combinations: 170 of
  170 identical.
- **`calls=` and `program_calls=` do move against `af2f400`.** Kepler
  fixed, fp64: 321,432 becomes 269,956, and 0 becomes 12,112.
- **The opt-in, against `3a00e18`.** Every forced combination identical,
  136 of 136. The default identical to the both-off run, 34 of 34.

### The verifier, three times

**On `d4e710b`: SENT BACK, 1.** No record moved: 108 of 108 program-half
pairs on the default routes, 120 of 120 on the forced combinations, and
56 library-half pairs. P1's six controls reproduced.
- **(a) `cases_gravity.c` held M4 only on unforced software.** M4 sets
  `route_on[]` to 0 for both routes after `engine_routes()`. Under
  `CFT_REBOUND_HIDE_FEATURES=scalar` it passed, reading `off by=device`,
  a state impossible on any device.
- **(b)** Three sabotages passed: a read past a table's source,
  `vprog_ix`'s NULL-source and zero-source refusals removed, and counters
  reporting the new route. P1's "no automated wide-format gate catches a
  swapped pair" was false at binary128.

**On `3a00e18`: READY.** `expected()` derives each route's state, six
`--door-refuse` cases hold `vprog_ix()`'s rules, and the swap limit is
restated, below. The library object is byte-identical to `d4e710b`'s, 61
of the program's 62 shared functions are identical, and the 62nd is
`main()`. 47 of 47 default and 55 of 55 forced pairs are identical. M4
fails all 10 states that then required on.

**On `0a0d96d`: READY.** No integration function changed: 75 of 78
functions identical in the library, 60 of 63 in the program. Records:
115 of 115, 23 of them the default against `3a00e18` forced off.
- O1 (auto ignores opt_in), O2 (by=opt_in while on) and O3 (forced on
  silently ignored) each fail `check_door` and `check_gravity`.
- O4 (a counter that moves only when forced on) fails `check_door`'s new
  part 4b, "(DIFFER: divsqrt_calls 18245 against 20753)".
- M4 now fails only `check_gravity --library`, in the full tier, and
  exactly the 5 of 14 states that force a route on.

### P1's three seams, and merge 4's `make check`

In `tools/check_seams.py` (entry 42), from `2b497a6`:
- **P1 x P3.** `indexed` and `pair_gather` join `reduce_seg` across a
  checkpoint, and the regrow sequence now runs too, so gravity's tables
  are rebuilt on both sides of the save.
- **P1 x P2.** A schedule recorded with both gravity routes forced one
  way is replayed with both forced the other, and must give the
  recording's lines and counters with no exit test computed.
- **P1 x P4.** The program engine, gravity forced on then off, gives the
  loop's record and counters, with `nextstep_calls` equal and above zero.

Measured at all three formats on `2b497a6`: **84 ok, 0 FAIL**. At
binary256 `indexed` carried 4,309 calls at the save and 8,520 fallbacks
at the resume of the plain sequence. The longest merged trailer found
is 455 bytes after the "# ", and the subprocess API keeps 511.

The full `make check`, from clean at `e47b04e`, read from its log:
- rc 0, 18 legs, all full, with no cache skips;
- **715 `ok`, 0 FAIL**, no DIFFER, traceback or make error;
- 10 `RESULT: PASS`, `check_gravity`'s now among them;
- the six digests unchanged: `e08541bc0af9b3ec`/`56bacf3bf63f19a7`,
  `5d7cf4588fad050e`/`cd16e6dc29aa0589`,
  `a1052750e206bebe`/`ba4d447c8e58bd20`;
- `check_seams`: 54 of 54 checkpoint crossings, 18 of 18 replays and 9
  of 9 program-engine pairs carried.

### Stated limits, not closed

- **A route that says on while today's runs.** V1 built it for valid
  tables, and it passes: on software both routes compute the same bits.
- **A swapped pair.** At binary64 `check_records` catches it. At
  binary128 the catches are in the full tier only. At binary256 nothing
  automated catches it: `gate_real`'s digests move, but it exits 0.
- **`make check-light` could not see M4 at `e47b04e`,** and `check_door`'s
  failing lines opened "all 9 counters identical". The fixes, `e90d5df`
  and `f9d2187`, landed with merge 5 (entry 44).
- **One run per force evaluation holds up to 257 bodies.** On a 256-deep
  device, 1,024 bodies takes 4 runs.
- **Gravity's images need SCRATCH and SCRATCH_IO.** A device without
  them is refused by name. Not measured on a real device.

### What this does not establish

- **`check_gravity` on the merged tree, on a tile.** It passed at
  `afc83d8`. Its leg, `p1-gate`, did not fit the window on 2026-11-06, so
  its opt-in branches are right by reading only. `check_door` holds the
  routes on a tile at the final tree.
- **A route that wins.** Neither was measured faster at 64 to 512 bodies.
- **128 bodies,** which ROUND3.md's card day named.
- **The whole-program divide and square root,** priced, not built:
  `divfull` is 212, 214 and 216 instructions, `sqrtfull` 189, 192 and 195.
- **Linux.** The ledger records no full `make check` there at these merges.

## 44. Round 3's P5: a constant goes to the library as one element, on an opt-in route; neutral on a card, and the route comparison made able to see its own trap

**2026-11-06.** Two hosts, as in entry 40. Pins: REBOUND at bdfda4bd and
cft-fp256 at 922c2ad (ABI 0.18). P5 branched from `e47b04e`. Its commits
on `round3-p5` are `2f42423`, `e98dbb1` and `777f160`. Merge 5 is
`c785263`: `c893e5a` merges P5 at `777f160`, and `c785263` is the lead's
`scalar` in `tools/check_seams.py`. The card ran `c785263`.

### What P5 is, in one paragraph

ROADMAP item 2. A constant operand was staged as N3 copies of one value.
- **24 sites** now pass it through P0's `vrun_k()` as one element, where
  the route is on, by `vmulk()`, `vfmak()`, `vkmul()` and `vaddk()`.
- **The route is opt-in from the start,** because a route no card has
  measured does not go on by default.
- **Off, each site is today's call on today's vector.** A default run
  is `e47b04e`'s byte for byte except the `# route scalar=` line: 66
  integrations at three formats, 1,407 stdout lines, whole stdout with
  only the wall-clock fields removed.
- **Left as they are:** 9 divides, because `cft_div` is not an opcode
  and takes no mask; 14 sites in `dtnew_*()` and `choose_timestep()`,
  out of scope; 2 E-wide compares, because `vcmplt()` counts no staged
  bytes.

| run | sites | a pass | an attempt | a prediction | an energy |
|---|---|---|---|---|---|
| `--engine loop --arith rebound` | 11 | 133 | 21 | 21 | 3 |
| `--engine loop --arith fma` | 15 | 98 | 50 | 21 | 3 |
| `--engine program --arith fma` (the card) | 9 | 0 | 50 | 0 | 3 |
| library half, REBOUND's form | 12 | 133 (+105 with a velocity-dependent force) | 21 | 21 | 0 |

So on the card's configuration the bytes fall 0.2-0.8% at binary128 for 2
to 64 bodies, against 2.8-7.4% on `--engine loop --arith fma`.

### What the plan had wrong

- **The route comparison could not see the trap, as first built.** `DTB`,
  `DTDB`, `RB`, `D1B` and `D2B` hold one value a system. Passed as a
  scalar, element 0 is system 0's step, right for one system and wrong
  for an ensemble. But `vrun_k()`'s fallback also broadcasts element 0,
  so both routes gave the same wrong answer: on T1, forced on against
  forced off differed in 0 of 15 lines, while 8 of 15 differed from the
  correct engine. P5 found this, and its gate held each ensemble to its
  members run alone.
- **"That is one flag" is four edits:** the flag, `check_door`'s `OPT_IN`,
  the library case's `expected()`, and the gate's default.
- **ROADMAP's "~300 staged vectors a step"** is 50 an attempt on the
  card's engine, and its 24 sites are not these 24.
- **`KY4` and `KY5`'s only sites are in `dtnew_prs23()`,** out of scope:
  5 sites and 17 calls an attempt, on the card's engine too.

### The change that closed it

With both routes through `vrun_k()`, the default route copied element 0
n times at every site, which `e47b04e` did not. `777f160` makes the off
branch today's call on the prebuilt vector. That gives two things:
- the default path is `e47b04e`'s exactly, so an opt-in route costs
  nothing to have merged;
- the routes can now see the trap, because off reads the whole vector
  and on reads element 0.

So `tools/check_scalar.py` has two oracles. Forced on must equal forced
off, and each ensemble's systems must equal themselves run alone.

### The gate, and the controls that show it can fail

At every format, `check_scalar --library` gives 90 ok, 0 FAIL, and
`check_door` 80 ok, 0 FAIL. Its cases include:

```
ok   kepler e=0.99, rejections that predict: the same record, schedule and counters on all three routes; 40 steps, 11 rejected, 51 attempts, 1470 passes; 11 rejected attempts predicted; the route carried 2559 calls and staged_bytes fell by 408096, both as derived
ok   check_dropin --light, forced onto the new route: rc 0, every case passed; scalar=on by=forced calls=5815433 fallbacks=0
```

| property | sabotage | with it |
|---|---|---|
| The trap | T1: `DTDB` as a scalar at 21 sites, count-neutral | every E = 1 case passes; 6 FAIL, on the two adaptive ensembles; every format 68 ok, 18 FAIL, all naming both oracles |
| The trap | T2 `DTB`, T3 `D1B`, T4 `RB`'s powers | 2, 2 and 4 FAIL, E = 1 all passing |
| `staged_bytes` derived | T5, R1, R2, R3 | 27, 18, 27 and 27 FAIL: "fell by 0, and the sites derive 689664" |
| The route's policy | O1 to O4 | `check_scalar` 27, 27, 27, 28; `check_door` 3, 3, 2, 5 |
| Library half | M4a, M4b | pass unforced `check_dropin`, fail the forced-on leg: 30 ok, 9 FAIL |
| Library sites | L1: a `predict_velocities()` site out | 2 FAIL: "added 4284 calls, and 68 passes of 105 are 7140" |

The gate's `staged_bytes` fall is derived from the recounted sites, their
widths and the run's call counts, not transcribed.

### No record line moved

- **P5's A/B against `e47b04e`, at `777f160`:** 66 integrations at three
  formats, whole stdout, one line differing a run, `# route scalar=`. On
  P5's first build, before `777f160` changed the default path, 93 of 93
  on all three routes, across all four adaptive modes, `--cs augmented`,
  softening and min-dt. V5's 171 rows, below, are at `777f160`.
- **Merge 5's `make check`,** from clean at `c785263`, read from its log:
  rc 0, 19 legs, all full, with no cache skips; **840 ok, 0 FAIL**, no
  DIFFER, traceback or make error; 11 `RESULT: PASS`; the six digests
  unchanged; `check_seams` 72 of 72 crossings, 18 of 18 replays and 9 of
  9 program-engine pairs carried.

### The verifier

**On `777f160`: READY, 0 (a), 4 (b).** V5 read all 24 sites. A scan finds
no write to any constant after it is built. An instrument in every
helper asserted each constant is one value over the whole run: 32 ok, and
its control died, `vmulk: the constant differs at lane 6 of 30`. With the
helpers undone the tip disassembles to `e47b04e`'s code in both halves,
every object section equal. Held to the base: 171 program rows, 4,272
lines, 0 differ, and 15 `gate_real` sequences, 1,284 lines. The whole
`check-light` passes with the route forced on, and `check-python` too.

- **(b) b1.** P5 said the memory half was blocked, every K vector being
  read whole with the route off. V5 measured otherwise. A K vector is
  `NMAX` wide, NMAX = max(N3, P), and only `G` is read past N3. See
  "What this leaves".
- **(b) b2.** A route that says on and is not passes every tier:
  `vrun_k()`'s on branch running today's `cft_run` while counting
  `route_calls`. `staged_bytes` is the door's own tally. Stated below.
- **(b) b3.** The recount's hold on the off branch was thin: H4, H6 and
  H7 passed. `9335a59`, in merge 6, tightens it, and V5 measured H4, H6
  and H7 failing there.
- **(b) b4.** Two sentences went stale at merge 5. From merge 5,
  `check_seams` forces `scalar` on an adaptive ensemble and fails T1, T2
  and T4 too; T3 still has only `check_scalar`.

### On the card

The leg `p5-scalar` ran at `c785263` on the rev8a single tile, where
`reduce_seg` was still on by device. All 24 card rows are identical to
software, and `SCALAR=on` carried 106 calls and no fallback on the route
check.

| problem | format | off, s | on, s | on against off |
|---|---|---|---|---|
| kepler | fp64 | 22.15 | 21.67 | -2.2% |
| kepler | fp128 | 42.68 | 42.33 | -0.8% |
| outer | fp64 | 9.92 | 10.45 | +5.3% |
| outer | fp128 | 21.70 | 22.12 | +1.9% |
| pythagorean | fp64 | 10.93 | 10.65 | -2.6% |
| pythagorean | fp128 | 20.50 | 18.08 | -11.8% |
| nbody64 | fp64 | 6.50 | 6.28 | -3.4% |
| nbody64 | fp128 | 16.21 | 15.89 | -2.0% |
| nbody256 | fp64 | 24.07 | 23.82 | -1.0% |
| nbody256 | fp128 | 48.44 | 47.65 | -1.6% |
| nbody512 | fp64 | 65.58 | 64.62 | -1.5% |
| nbody512 | fp128 | 118.87 | 119.28 | +0.3% |

One run per row; entry 45 gives the noise. Measured: the route is faster
on 9 rows of 12, with a median of -1.5%. Two rows stand past what two runs
of one configuration differ by, pythagorean fp128 at -11.8% and nbody64
fp64 at -3.4%, and one run a row cannot say whether they are real. So it
neither wins nor loses measurably, and it stays opt-in, since a route
leaves the list when a card measures it winning. Of ROADMAP's two reasons
the honesty of having asked stands. The memory needs no route, below.

### Stated limits, not closed

- **A route that says on and is not** passes every tier (b2).
- **What both routes share.** A site given the wrong constant has no
  route or member line; `check_equivalence` catches it directly.
- **The library FMA form** is held by its call counts only, because the
  engine opens once a process.
- **The recount reads only the wrappers it names.**

### What this leaves

The K vectors, measured by V5 on the default route: with them N3 wide, 66
of 66 integrations are the same in every record line and in
`staged_bytes`, and the peak falls 80.6% at N = 256 binary64 and 90-95%
from N = 512 up: 1,206 MB to 97 MB at N = 1024 binary64, and 4,682 MB to
249 MB at binary256. The library half is not built. Filed for the next
round.

### What this does not establish

- **A full `make check` on Linux.** The ledger records none at these
  merges.
- **The memory variant in the library half.**

## 45. Round 3 on the card: the merged tree gives software's bits on a tile and 21-36% off a run up to 64 bodies; every new route lost or tied, and all four are opt-in

**2026-11-02 and 2026-11-06.** amd-arc-box's U50, one tile at a time. Pins:
REBOUND at bdfda4bd and cft-fp256 at 922c2ad (ABI 0.18). The images are
the f128 single (`90457667a76df67b...`, entry 37's), the revision-8
single, "rev8a" (`0639248e9c66f779...`), the quad (`4f9921a25af65444...`)
and the deep single, "rev8d" (`bbf08fb3feb60940...`). Three commits ran
that are not on main, `c3471bc`, `afc83d8` and `aa3aa7e`, each a parcel or
staging commit with `hw/cardday-round3.sh`. Session r3a ran `63a4e4f`,
which is `af2f400` with the card script, and this entry calls it
`af2f400`.

This entry is what the sessions measured across the parcels. Entries 40
to 44 carry each parcel's own card rows.

### The pin on the card

Entry 38's missing row. At `af2f400`, through libcft 0.18:
- **The f128 image opens.** `gate_real` passes at binary64 and binary128,
  and binary256 is refused by name: "carries fp32 fp64 fp128".
- **On rev8a, what ran passes:** `gate_real` at all three formats, and
  `check_dropin`, **131 ok, 0 FAIL, "every case passed"**, in 3 h 34 min.
  That leg's `check_dropin --wide` and `check_ensemble` were left for a
  later window.

### The round's before and after

`hw/bench-workload.py` on the rev8a single tile, default engine for the
card (`--engine program --arith fma`), every card row byte-identical to
software. "Before" is `af2f400`. "After" is the merged engine, at
`e90d5df`, which has `e47b04e`'s, with the device exit test off, which is
the default as the round ends. Both columns are `workload.csv`'s
`seconds`, the bench's clock around the process, one run per row:

| problem | format | before | after | change |
|---|---|---|---|---|
| kepler, 200 steps | fp64 | 26.57 | 20.93 | -21.2% |
| kepler | fp128 | 53.87 | 37.97 | -29.5% |
| outer, 100 steps | fp64 | 13.66 | 9.45 | -30.8% |
| outer | fp128 | 31.78 | 21.38 | -32.7% |
| pythagorean, 100 steps | fp64 | 13.82 | 9.52 | -31.1% |
| pythagorean | fp128 | 24.96 | 18.01 | -27.8% |
| nbody64, 20 steps | fp64 | 9.66 | 6.30 | -34.8% |
| nbody64 | fp128 | 24.45 | 15.55 | -36.4% |
| nbody256, 5 steps | fp64 | 25.57 | 23.38 | -8.6% |
| nbody256 | fp128 | 51.33 | 47.48 | -7.5% |
| nbody512, 3 steps | fp64 | 61.00 | 62.58 | +2.6% |
| nbody512 | fp128 | 116.56 | 116.84 | +0.2% |

**What the noise is.** Two benches ran one configuration in different
sessions, `REDUCE_SEG=on` at `e90d5df` and `SCALAR=off` at `c785263`, whose
off route P5 and V5 measured as `e47b04e`'s calls. Row against row they
differ by -5.8% to +5.7% on kepler, outer and pythagorean, and by -0.6% to
+1.9% on the n-body rows. Add the deep image, which should cost the same,
and a small-problem row spans up to 8.3%, outer fp64. So a single row is
good to about 6-8% and 2% at best. `SCALAR=on` would widen those to 13.4%
and 3.5%, but it is the route under test (entry 44).

Measured, with that noise: a row is 21-36% cheaper up to 64 bodies, about
8% cheaper at 256, and no better at 512, where the rows move +0.2% to
+2.6%. Believed: most of the gain is P1's program form and P4's program.
No bench separates them on one commit. `d4e710b`'s both-off rows, which
have P1 and not P4 or P2, are -16% to -35% up to 64 bodies, -8% to -10% at
256, and +1.5% to +6.3% at 512.

**At 256 and 512 bodies a row is not a step.** Each row evaluates the
energy at its first and last sample, and `energy_all()` sums it with one
width-one call a pair and one a body. That is 65,792 of the row's 68,584
calls at 256 bodies and 262,656 of 264,517 at 512. At those rows' 112-119
us a call it is about 7.4 s of the 256 rows and 31 s of the 512 rows, in
both columns, and half of the nbody512 fp64 row. Estimated with it taken
out, at that cost on both sides, a step moves -9% to -12% at 256 bodies
and +0.4% to +5.5% at 512. That is arithmetic on the rows. No run measured
a step apart from its energy.

**Where the gain goes at 512 bodies.** The program's own buckets, before
and after, in seconds:

| bucket | binary64 | binary128 |
|---|---|---|
| elementwise calls, `t_elem`, about 31 s of each the energy | 39.24 to 31.47 | 54.71 to 31.23 |
| program runs, `t_program` | 0.21 to 6.52 | 0.56 to 19.05 |
| divide and square root, `t_divsqrt` | 19.26 to 22.20 | 56.62 to 62.26 |

So the round moved work from elementwise calls into program runs, and at
512 bodies the runs and the divides take back what the calls saved.
Measured, not explained: the same divide and root calls take 10-19%
longer at 256 and 512 bodies.

### What the software backend paid

The gain is a card's. On this desktop, with nothing else running, `af2f400`
against the tree this entry lands with, on the software backend, the
bench's arguments, best of 5 runs and of 3 from 64 bodies up, parent's
clock:

| problem | REBOUND's form, fp64 | fp128 | `--engine program --arith fma`, fp64 | fp128 |
|---|---|---|---|---|
| kepler, 200 steps | 2.44 to 2.94 s, +20.5% | +19.7% | 0.98 to 1.38 s, +39.7% | +37.7% |
| outer, 100 steps | 3.11 to 3.64 s, +17.0% | +13.1% | 1.47 to 1.91 s, +29.9% | +26.3% |
| pythagorean, 100 steps | 1.56 to 1.86 s, +19.0% | +15.1% | 0.69 to 0.92 s, +34.1% | +32.7% |
| nbody64, 20 steps | 17.52 to 18.70 s, +6.8% | +0.6% | 12.48 to 13.50 s, +8.1% | +2.1% |
| nbody256, 5 steps | 61.36 to 62.91 s, +2.5% | -2.8% | 57.09 to 59.88 s, +4.9% | +2.6% |

Measured: `t_program` is what rises. Kepler in REBOUND's form at binary64
goes from 0 to 0.53 s of it, while `t_elem` falls from 0.41 to 0.38 s on
84,915 fewer calls. So on software a program run costs more than the
calls it replaces, and the gap closes as the pairs grow. REBOUND's form is
the drop-in's, so a user with no card pays 13-21% on a few-body problem.
Gravity's off route is the program form too (entry 43), so no switch
restores the old calls. No entry said this until the docs check at the
round's end found the bench's software rows had moved.

### The gravity split, ROADMAP item 4's table

At `af2f400`. Each figure is a share of the row's wall clock, and the
buckets nest: `t_gravity` is about `t_g_lib` plus `t_g_copy`, and
`t_g_lib` holds `t_g_accum`, `t_g_divsqrt` and the pair arithmetic.

| problem | format | wall s (rev8a / f128) | gravity | host gather | accumulate | divide and sqrt in gravity |
|---|---|---|---|---|---|---|
| kepler | fp64 | 25.9 / 22.4 | 58.3% | 0.0% | 2.2% | 11.6% |
| kepler | fp128 | 53.1 / 46.2 | 69.7% | 0.0% | 2.6% | 14.5% |
| outer | fp64 | 13.0 / 11.3 | 59.5% | 0.0% | 9.7% | 10.9% |
| outer | fp128 | 31.0 / 25.4 | 72.2% | 0.1% | 11.6% | 13.9% |
| pythagorean | fp64 | 13.0 / 10.4 | 57.3% | 0.0% | 4.1% | 11.9% |
| pythagorean | fp128 | 24.3 / 20.1 | 68.8% | 0.0% | 5.0% | 14.0% |
| nbody64 | fp64 | 9.0 / 9.1 | 81.2% | 1.2% | 38.8% | 26.4% |
| nbody64 | fp128 | 23.7 / 24.2 | 88.7% | 1.2% | 38.3% | 33.3% |
| nbody256 | fp64 | 24.8 / 27.0 | 64.7% | 2.1% | 20.0% | 37.6% |
| nbody256 | fp128 | 50.6 / 51.9 | 82.2% | 2.1% | 20.4% | 51.0% |
| nbody512 | fp64 | 60.0 / 69.9 | 46.7% | 2.0% | 10.2% | 31.4% |
| nbody512 | fp128 | 115.4 / 139.2 | 71.0% | 2.8% | 12.9% | 48.6% |

The wall column is the program's own clock, `wall_seconds`, 0.6 to 1.2 s
under a bench's `seconds`, so it is not the "before" above. The
percentages are rev8a's. rev8a is 15-25% slower than the f128 tile on the
three small problems, about even at 64 bodies, and 14-17% faster at 512.
Measured: the host's gather is at most 2.8% of a step, so a device gather
had almost nothing to win, and the device gathers then lost (entry 43).
Inside gravity the divide and square root are the largest cost from 256
bodies up, 31-51% of a row. Estimated without the energy evaluations, they
are 54-67% of a step at 256 and 512 bodies before the round, and 73-74%
as it ends, at both formats.

### Every new route, on a tile

| route | what it carries | against today's route on rev8a | as the round ends |
|---|---|---|---|
| `indexed` | gravity's accumulate fold | +22% to +40% at 64-512 bodies | opt-in |
| `pair_gather` | the pair program's gather | +22% to +38% at 64-512 bodies; with `indexed`, +47% to +74% | opt-in |
| `reduce_seg` | the corrector's exit test | slower on 12 rows of 12: +0.8% to +4.8% on the n-body rows, +2.0% to +11.4% on the small ones | opt-in |
| `scalar` | a constant operand | 9 rows of 12 faster, median -1.5%; neither wins nor loses measurably | opt-in |

All four give software's bits: every card row of every bench is
identical. `reduce_seg`'s rows are:

| problem | format | off, s | on, s | on against off |
|---|---|---|---|---|
| kepler | fp64 | 20.93 | 21.35 | +2.0% |
| kepler | fp128 | 37.97 | 42.16 | +11.0% |
| outer | fp64 | 9.45 | 10.53 | +11.4% |
| outer | fp128 | 21.38 | 22.35 | +4.5% |
| pythagorean | fp64 | 9.52 | 10.37 | +8.9% |
| pythagorean | fp128 | 18.01 | 19.40 | +7.7% |
| nbody64 | fp64 | 6.30 | 6.48 | +2.9% |
| nbody64 | fp128 | 15.55 | 16.30 | +4.8% |
| nbody256 | fp64 | 23.38 | 24.22 | +3.6% |
| nbody256 | fp128 | 47.48 | 47.86 | +0.8% |
| nbody512 | fp64 | 62.58 | 64.36 | +2.8% |
| nbody512 | fp128 | 116.84 | 118.72 | +1.6% |

The direction is measured, on 12 rows of 12. The size is less sure. Outer
fp64's +11.4% is mostly noise: with the test on, the deep image took 9.72 s
and `SCALAR=off` 9.92 s. Kepler fp128's +11.0% may not be: four runs with
the test on took 42.16 to 42.68 s, against one run off at 37.97 s.

**So `auto` takes none of the four.** The door keeps each route and its
switch, and a route leaves the list when a card measures it winning. The
owner chose this for gravity's two on 2026-11-02 and for `reduce_seg` on
2026-11-06, and `scalar` was opt-in from the start. `eb6976a` is the
flip for `reduce_seg`, and V5 verified it:
- 171 rows, 7 runs a row, 0 differ against `c785263`;
- the default equals `c785263` forced off in every line but the by=
  word, `staged_bytes` included;
- all four forced on equals the default in records and counters;
- M4, a library half silently on today's route, passes `check_door` and
  the default drop-in, and fails only the forced-on `check_corrector
  --library` that the flip added.

### The ensemble on a tile, and what idle lanes cost

`hw/bench-ensemble.sh` at `e830230`, on the rev8a single tile: 20 steps of
fixed dt 0.05, the FMA form on the program engine. "Together" is members
one ulp apart. "Apart" is member k's planet faster by k/2000, so they span
eccentricities. All 8 card rows are identical to software.

| members | format | together: card s, against software, lane efficiency, passes | apart: card s, against software, lane efficiency, passes |
|---|---|---|---|
| 64 | fp128 | 6.46, 1.82x, 0.943, 8.062 | 6.32, 1.82x, 0.940, 8.039 |
| 64 | fp256 | 14.56, 2.21x, 0.976, 18.247 | 14.68, 2.17x, 0.975, 18.180 |
| 512 | fp128 | 17.77, 4.74x, 0.910, 8.056 | 17.06, 4.84x, 0.907, 7.887 |
| 512 | fp256 | 50.10, 4.75x, 0.966, 18.258 | 49.25, 4.78x, 0.940, 17.633 |

Measured:
- **A tile beats software here,** 1.8-2.2x at 64 members and 4.7-4.8x at
  512. `docs/bench/ensemble.csv`'s binary128 rows on the f128 tile, from
  2026-09-14, are 10.70 s and 24.66 s, so these are 40% and 28% less.
  They are a different image and engine, and nothing separates the two.
- **Idle lanes are 2.4-9.3% of lane-passes.** Apart moved lane efficiency
  by at most 2.6 points, and card seconds by -4.0% to +0.8%, inside the
  noise above. No row rose past it.

So this does not price the lane mask at a wide spread. HARDWARE's item 4
asks for members needing 10 and 23 passes in one step. This perturbation,
over 20 steps, moved the mean by at most 0.6 of a pass.

### The other images

- **A quad opened as one tile, `CFT_XRT_TILES=1`,** against the rev8a
  single, on the three small problems: 1-10% slower, every row identical
  to software. Believed, not separated: some of that may be the two
  images' kernel clocks.
- **The deep single, rev8d.** `max_scratch: 8192` against the single's
  4,096, `seq_features: 0x7ff1f`. 6 of 6 rows identical to software, and
  5 of them within 2% of rev8a; the sixth, at -7.7%, is inside the span
  above. The port's widest program, gravity's fold, uses 256 slots, well
  under the single's 4,096, so the depth has nothing to buy here, and it
  costs nothing.

### The gates on a tile

| gate | commit | result |
|---|---|---|
| `check_door` | `af2f400` | 53 ok; 131 of 131 cases on both routes at each format |
| `check_door` | `e90d5df` | 77 ok, 0 FAIL; 164 of 164 at each format; the opt-in branches measured |
| `check_door`, f128 tile | `af2f400` | 33 ok; 98 of 131, with `indexed` and `reduce_seg` absent |
| `check_program_engine` | `c3471bc` | 15 ok, 0 FAIL, in 56 min |
| `check_gravity` | `afc83d8` | 23 ok, 0 FAIL, in 59 min |
| `check_corrector` | `85fd0aa` | 277 ok, 0 FAIL, in 58 min |
| `gate_real`, identity | `85fd0aa` | PASS at all three formats |
| `check_ensemble` | `c785263` | 27 ok, 0 FAIL, at all three formats, 3 to 5 members, in 4 h 27 min |
| `check_door` | `e830230` | 77 ok, 0 FAIL, in 780 s; 164 of 164 at each format; all four routes forced on carry their calls with no fallback, and the records and 9 counters equal the default's |
| `check_seams --quick` | `e830230` | 36 ok, 0 FAIL, at binary64, in 30 min; 24 of 24 checkpoint crossings, 6 of 6 replays and 3 of 3 program-engine pairs carried |
| `gate_real --fp64`, all four routes forced on | `e830230` | PASS, in 6 min 46 s; all 8 digest lines software's |

`check_seams` takes each opt-in route through the library half on a tile
at the final tree, one at a time: forced at a save and off at the resume,
then the other way, each crossing giving the straight run's 17 values,
with `CFT_REBOUND_NO_ADOPT` differing. In the plain sequence `scalar`
carried 19,012 calls at the save and `indexed` 928.

The window's last minutes ran `gate_real --fp64` with all four on at once.
All 8 digest lines equal software's default in merge 6's suite,
`e08541bc0af9b3ec` and `56bacf3bf63f19a7` among them, and it printed twice
"the resume is REBOUND's own ias15 across its own checkpoint, 30 values
bit for bit". A straight run under `--routes` carried `scalar` 56,105
calls, `indexed` and `pair_gather` 2,735 each and `reduce_seg` 730, with
no fallback.

Legs ran from well under their estimates to 3 times them: `check_ensemble`
took 4 h 27 min against 90 minutes, and `mixed` 15 minutes against 60.
`seams`, whose 5 minutes was a placeholder, took 30. The Python checkers
open the device once per case, as docs/HARDWARE.md warns.

### Merge 6, the round's last

`eb6976a` is the flip above. `6b2b5f0` merges P5's `9335a59` (entry 44).
`e830230` adds the leg `dropin-forced`, after V5 caught that with all four
routes opt-in an unforced drop-in on a tile takes no device route, and
would have read as a pass for them. `148bb50` takes V5's follow-ups:
`check_corrector` reads a route's state as well as its by=, both gates
hold `nextstep_calls`, so 10 counters, and `check_corrector`'s header
states the limit below. V5 on `6b2b5f0`: READY, 0 (a), 5 (b); on `148bb50`,
four of the five closed, and two comments left. The commit that adds this
entry fixes those two.

The full `make check`, from clean at `148bb50`, read from its log:
- rc 0, 19 legs, all full, with no cache skips, in **2 h 36 min**, the
  first time with little else on the host since the suite grew;
- **841 ok, 0 FAIL**, no DIFFER, traceback or make error;
- 11 `RESULT: PASS`;
- the six digests unchanged: `e08541bc0af9b3ec`/`56bacf3bf63f19a7`,
  `5d7cf4588fad050e`/`cd16e6dc29aa0589`,
  `a1052750e206bebe`/`ba4d447c8e58bd20`;
- each route family forced through `check_dropin --light`, against
  REBOUND, on software: 5 cases carried 574 maxima on the device's exit
  test, `indexed` carried 236,346 calls and `pair_gather` 236,397, and
  `scalar` 5,815,433;
- `check_seams`: 72 of 72 crossings, 18 of 18 replays and 9 of 9
  program-engine pairs carried.

`make check-light` and `make check-quick`, each from clean with its build,
on the tree this entry lands with, docs aside: rc 0, 19 legs, **550 ok, 0
FAIL**, 11 `RESULT: PASS` each, in 31 and 42 minutes, and 13 s warm.

### Stated limits, not closed

- **A route that says on and is not** passes every tier, for all four
  routes. `route_calls` and `staged_bytes` are the door's own tallies, and
  nothing derives `reduce_seg`'s rise in `staged_bytes`.

### What this does not establish

**Not run on a tile.** `check_ensemble` took 4 h 27 min of the day's 9
hours, and a drop-in run cannot stop between cases, so these were left
rather than risk the card's nightly job:
- **`check_dropin` at any merged tree.** It last ran on a tile at
  `af2f400`, so `dropin-forced`, `dropin-light` and `dropin-wide` remain.
  What a tile has against REBOUND's own ias15 at these trees is `gate_real`
  at binary64, with all four forced at `e830230` and with `reduce_seg` on
  at `85fd0aa`: one sequence, not the drop-in's case-by-case comparison.
  The default the round ends with has not met REBOUND on a tile.
- **`check_gravity` at the merged tree,** the leg `p1-gate`.
- **`check_scalar` at all.** The leg `p5-scalar` is two benches and a route
  check.
- **`check_program_engine` at a merged tree.** It ran at `c3471bc` only.
- **`check_seams` above binary64, and `gate_real` at binary128 and
  binary256 after `85fd0aa`.**
- **Anything on the f128 tile at a merged tree,** but r3c's probe.

- **Noise beyond two runs.** One configuration ran twice, and its small
  rows a third time on the deep image. A row's spread is an estimate.
- **What separates P1 from P4 in the gain.** No bench holds one fixed.
- **A step apart from its energy.** The figures for that are estimates.
- **Why the runs and the divides rise at 512 bodies.** The buckets say
  where the gain goes, and not why.
- **A second card, or a second host.** Every card row is one U50's.

## 46. Round 4's P0: a problem's `active` and `zonal` lines, and a seam for the zonal term; sent back once, and changed in no record line

**2026-11-11.** This desktop, on the software backend only, shared with
cft-fp256's agents and at 100% load at 17:22; the suites' logs print
mpmath 1.3.0, and no compiler version. Pins:
REBOUND at bdfda4bd and cft-fp256 at 922c2ad (ABI 0.18). P0 is the
lead's. Its commits on `round4-p0`, from `7a7a279`, are `4b52705` and,
after V0, `9a42b1d`. main fast-forwarded to `9a42b1d` and was pushed on
2026-11-11.

ROUND4.md is the round's plan and has P0's design; docs/RING.md is what
the round left. The round's ledger is archived as
`archive/round4-ledger.zip`. This entry is what ran.

### What P0 is, in one paragraph

Every parcel reads or writes problem files, and three of them write
src/ias15_cft.c, so the shared part landed first.
- **`active n`.** The first n bodies of every system, 1 to N, are
  massive and the rest test particles, through the engine's own
  `eng_active` with `testparticle_type` 0, the path the drop-in already
  takes. The record says so in a header line with the pairs its rows
  fold, `# active 5 of 8 bodies a system: 25 pairs a force evaluation, a
  system`. Only a problem with the line prints it, so no older record
  moves.
- **`zonal k R J2 J4`.** Parsed, stored and keyed, then refused by name
  until P1 built the term.
- **The seam.** `src/ias15_zonal.inc`, empty and unreachable.
  `gravity_body()` calls it after the fold, and `energy_all()` and the
  record's header call it too.
- **The run key** names both lines, only when a problem has them, so no
  older key moves.
- **The energy follows `N_active`,** from `9a42b1d`: `energy_all()` runs
  REBOUND's two loops. With no active line every body is active and the
  sums are the ones they were.

### The gate, and the controls that show it can fail

`tools/check_gravity.py` gains two parts. The problem is outer with two
more massless bodies, made from Pluto by a quarter and a half turn about
z, so every value is exact.
- **2b.** `active 5` against every body active, at every format, on the
  loop engine in both forms and on the program engine. Each of the 9
  rows, here the first from `291c99c`'s suite:

```
ok   fp64, loop engine, rebound form: 5 sample lines identical, 25 pairs a force evaluation against 28 all active, the test particles move, a 2^-40 mass changes nothing under active 5 (energy included) and everything without it, and the two-system ensemble and its member agree
```

- **2c.** Refused by name: an active count of 0 or past N, and a second
  active line. The zonal line was refused here until P1's merge made the
  case an acceptance, and `291c99c` added three refusals (entry 47).
- **The replay key.** A schedule recorded with `active 5` replays with
  it and is refused without it.

The plan's control, a mass of 1e-30, would be inert at binary64, below
its last bit, so the gate's is 2^-40. Four sabotages, each built from a
copy and deleted by the lead:

| sabotage | with it |
|---|---|
| (a) the active line parsed and `eng_active` never set | FAIL at all three engine and form rows: no 25-pair line, and the 2^-40 control does not bite |
| (b) the run key not naming the line | FAIL: the replay without it is accepted, rc 0 |
| (c) `testparticle_type` 1 | FAIL at all three rows |
| (d) the energy counting every body | FAIL at all three rows |

(c) and (d) are the gate's answer to V0's two findings, below.

### The verifier, twice

**On `4b52705`: SENT BACK.**
- **The untouched path held.** 117 runs of data/problems/ and a
  4-member ensemble, at binary64, binary128 and binary256 on both engines
  and in both forms, were byte-identical between `7a7a279` and `4b52705`
  once the timing fields were masked. The schedule key was
  byte-identical in 9 of 9 cases, and the seam is unreachable in all
  three builds.
- **Defect: `energy_all()` summed every body.** REBOUND's
  `reb_simulation_energy` (third_party/rebound/src/tools.c:108) counts
  only the active bodies under `testparticle_type` 0. The gate missed it
  because its test particles were all massless.
- **Defect: `testparticle_type` 0 was unheld.** A copy at type 1 passed
  `check_gravity --quick`.
- **Minor.** `active 0` was accepted and then aborted on a 0/0 at the
  first step, and a second active line silently won.

`9a42b1d` runs REBOUND's two loops in the energy, requires that a test
particle's 2^-40 mass change nothing under `active 5`, the energy column
included, adds a two-system ensemble and its member, and refuses
`active 0` and a second line.

**On `9a42b1d`: READY.** The record's E column equals REBOUND's own
`sim.energy()` (the rebound 5.1.1 wheel, with `N_active` and
`testparticle_type` 0) bit for bit at 7 samples on three problems:
kepler under `active 1`, the star's kinetic energy alone; outer and two
test particles under `active 5`, one of them at mass 2^-40, where the
all-active sum would differ; and kepler under `active 2`. That compared
energies at the same states; no stepped trajectory went through the
wheel. 117 of 117 untouched-path runs identical, and both new controls
bite at all three rows. Believed by V0: `testparticle_type` 0 is held as
a consequence, a test particle's mass changing nothing, not by an
assertion of the type itself.

### No record line moved, and the suites

- **The first full `make check`, on `4b52705`,** ended rc 2 at leg 7 of
  19. `check_corrector`'s forced drop-in run (`check_dropin --light` with
  `CFT_REBOUND_ROUTE_REDUCE_SEG=on`) returned rc 1 with no cases, and the
  log ends in MSYS fork failures: `child_copy`, Win32 error 299, and
  "couldn't create signal pipe", Win32 error 5. Run alone, the same run
  gave rc 0 and 170 ok. The host was at 100% CPU. Believed by the lead:
  environmental.
- **The full `make check` on `9a42b1d`,** from PowerShell, read from its
  log: 18:49:08 to 21:08:53, 2 h 20 min on a shared host;
  - rc 0, 19 legs, each `--mode full` once;
  - **857 `ok`, 0 FAIL**, no fork errors; 11 `RESULT: PASS`,
    `gate_real` PASS at all three formats, `gate_subprocess` PASS,
    `check_dropin` "every case passed", `--wide` passed, `check_archive`
    gates 1-3 PASS, and "78 constants compared, all match";
  - the forced drop-in run that died on `4b52705` passed;
  - `check_gravity` 2b at all 9 rows, and 2c's 4 refusals and the
    replay-key case;
  - the six digests of entry 45 are in the log:
    `e08541bc0af9b3ec`/`56bacf3bf63f19a7`,
    `5d7cf4588fad050e`/`cd16e6dc29aa0589`,
    `a1052750e206bebe`/`ba4d447c8e58bd20`.
- **Counting `ok` lines depends on the pattern** (P4). That log has 857
  lines matching `^  ok ` and 858 matching `^\s+ok\b`, because one
  `check_dropin` line is indented one space. Compare counts by one
  pattern only; the full suites' counts in entries 46 to 51 are by the
  first.

### Stated limits, not closed

- **`tools/state_to_problem.py` carries neither the `active` nor the
  `zonal` line.** READ by study DEBRIS-A (2026-11-12) and by the lead:
  it writes the name, G, N and body lines only. So `tools/reversal.py`,
  run on such a problem, would integrate back on a different problem
  (BELIEVED; not run). No gate runs either tool. P0 changed the grammar
  and no one owned this tool: not fixed in round 4. **Fixed by entry
  52,** which measured the consequence: a zonal problem's return missed
  by 3.95e-3 of its size.
- **`tools/make_ensemble.py` refuses `active` lines** (V0). An ensemble
  with test particles is written whole, as P2's generator does.
- **A problem file through a pipe is refused, not misread.**
  `read_problem()` rewinds between its two passes, and a pipe cannot
  rewind. Measured by P3 on 2026-11-12: kepler.txt through a Windows
  named pipe is refused by name, "read 0 bodies, expected 2", rc 3. The
  dt file's reader had the same shape and misread (entry 49); this one
  refuses, so it stays.
- **Inline comments,** from `291c99c`, are refused by name on the active
  and zonal lines, and ignored after a body line's seventh field, as
  before round 4.
- **P0's parser accepted a second zonal line and a fifth field
  silently.** That was unreachable while the line was refused. V1 found
  it once P1 lifted the refusal, and `291c99c` refuses both (entry 47).

### What this does not establish

- **The card, or Linux.** No artifact was opened and nothing ran in
  WSL, in P0 or anywhere in round 4.

## 47. Round 4's P1: J2 and J4 of one body in libcft's format, held to the golden model, to mpmath's gradient and to theory's rates; and a format difference cannot see a narrowing both formats share

**2026-11-11 and 2026-11-12.** This desktop, on the software backend,
shared with the other three parcels and their verifiers. Pins as in
entry 46. P1 branched from `9a42b1d`. Its commits are `c0c99be`,
`06adb1a`, `b3f664b` and `22ca3de`. `4dd53f3` merges `22ca3de` on
`round4-int`, and `291c99c` is the lead's fix after V1, on top. main
fast-forwarded to `291c99c` and was pushed on 2026-11-12.

### What P1 is, in one paragraph

`src/ias15_zonal.inc` fills P0's seam.
- **The term.** Body k's J2 and J4 on every other body of its system,
  in the run's own format, every operation a correctly rounded library
  call, no FMA, added after gravity's fold. So one order serves both
  `--arith` forms.
- **The reaction on k,** `(G m_j) prefac fac d`: the plan's "term times
  m_j/M", with no divide by M. It is taken only from the bodies that act
  on active ones, so under `active n` a test particle feels the term and
  gives no reaction, whatever its mass, as `testparticle_type` 0 says.
- **The potential** in `energy_all()`, after the pair terms, and a
  `# zonal` record line naming R, J2, J4, the bodies that feel the term
  and those whose reaction k takes.
- **The order** is REBOUNDx's `gravitational_harmonics.c` as recalled,
  not as read. Believed: REBOUNDx is not pinned (decision 4), and the
  include's head says so.
- **`--zonal-selftest`** evaluates the term and the potential at fixed
  points and traces every library operation, its operands and result as
  bits, and every element moved between them.
- **The constants.** `tools/gen_constants.py` derives the 19 integer
  coefficients from P2 and P4, and its check 3 holds them to mpmath's
  derivative of the potential. The header and the JSON gain them,
  additions only.
- **`data/problems/zonal_saturn.txt`:** Saturn with Cassini's J2 and J4,
  a Mimas-like moon inclined 3 degrees and a massless ring particle at the
  Encke gap inclined 2 degrees, under `active 2`. `check_zonal.py
  --write-problem` derives it, and the gate holds the file to that byte
  for byte. Its R = 60,330 km is believed, GPU-E's belief: its J4 table
  names no reference radius.

**What the term costs, measured by V1.** Counted by hand, and equal to
the trace: a body a force evaluation pays 3 sub, 34 mul, 14 add, 1 sqrt
and 10 div, two of the divides by 2 and 8 and so exact. P1 first counted
"21 mul (+13 for a reacting body)". The 13 reaction multiplies are as
wide as the bodies, so every body pays them, test particles included,
and a test particle's are thrown away. GPU-E prices the term at "1 DIV,
21 MUL, 6 FMA, 6 ADD" a body an evaluation (its section 2), an estimate
made before the code existed.

### The gate, and the controls that show it can fail

`tools/check_zonal.py`, one leg in each tier:
1. **Arithmetic.** Every traced operation against cft-fp256's golden
   model, bits and flags, and every operand a declared input or the
   trace's own last result.
2. **Formula.** 14 points, 10 off the plane and 4 on it, and each
   system's reaction and potential, against minus the gradient of the
   potential, taken by `mpmath.diff` at 200 digits. The bound is derived
   from the trace's own dataflow, and the dataflow evaluated exactly must
   equal mpmath's within 1e-150.
3. **Physics,** on zonal_saturn: (a) energy and (c) momentum at every
   format, and (b) the apsidal and nodal rates at binary64 against
   first-order theory in J2 and J4.
4. **Format difference.** binary128 and binary256 on binary128's own
   steps.

Then both engines, an ensemble against its members, the 2^-40 test
particle, the record line, the run key's `zonal=` field and three
refusals. From `291c99c`'s suite, 87 `ok`:

```
ok   126 operations, 1382 elements, each against the golden model's bits and flags (flags seen 0x10), ...: 0 mismatches   [at each format]
ok   fp256: the bound is 1.33e-69 of the term at its widest, against the 1e-68 asked
ok   fp64, 1100 steps adaptive: the energy, Phi's J terms included, moves 1.56e-15 of itself; the same orbits without the term, 1.22e-15 (the integrator's own level; the bound is 10 times it)
ok   the apsidal rate (wanted positive: the apsides advance): 6.61078e-7 rad/s over 30.2 orbits; first-order theory 6.6178e-7 (n - kappa); off by 0.00106, against the second order's 0.00359 (...)
ok   the nodal rate (wanted negative: the node regresses): -6.5874e-7 rad/s over 30.2 orbits; first-order theory -6.6178e-7 (n - nu); off by 0.00459, against the second order's 0.00866 (...)
ok   the nodal rate, J4 resolved: the same theory without J4 gives -6.4278e-7, 0.0248 off the measured rate (wanted past the second order's 0.00866)
ok   121 samples on the same 120 steps: the positions differ by 3.4e-33 of the distance from Saturn at most, against binary128's round-off, u128 * steps^2 = 1.39e-30 (binary64's: 1.6e-12), and must differ at all; the term itself moved them 0.201. Its scope: a precision that differs between the formats, never a narrowing both share (parts 1 and 2 hold that)
```

Seven sabotages, each built by P1 from a copy of src/, run through
`check_zonal --quick` and deleted:

| sabotage | with it |
|---|---|
| the term rounded to binary64 before it is added, in both runs (the plan's control for gate 4) | **gate 4 passes**: 3.371e-33 over 120 steps, the correct build's 3.402e-33. Gate 1 fails at operation 29; gate 2 at 1.4e-17 to 8.3e-17 of the term, against bounds near 1e-32; 3(c)'s momentum at 2.22e-20 against 1.55e-34 |
| the same rounding at binary128 only (gate 4's control, from the ruling below) | gate 4 fails: 4.0e-18 over 120 steps |
| J4's z part dropped, through the trace's own arithmetic | gate 1 passes; gate 2 fails at all 30 checks off the plane (at u = -0.998, 2.37e-2 of the term against a bound of 3.16e-14) and passes all 12 on it |
| the force's z sign flipped, the potential left alone | 3(a) fails at every format, the energy moving 4.35e-5 at binary64 and 4.16e-5 above, against 1.2e-15 and 7e-31; the node advances, +1.32e-6 rad/s |
| the reaction on k dropped | 3(c) fails, momentum moving 5.97e-3 and 5.06e-3 of \|m v\|; gate 2's line for body k fails at 0.998 of the reaction |
| one traced result one ulp off (operation 54, element 5) | gate 1 alone fails, at each format, naming it |
| the reaction taken from test particles too | gate 2's line for body k fails at 1.19e-4 against 1.5e-14, and the 2^-40 case says "THEY DIFFER" |

**P1 found a dead check in its own gate.** The dropped-reaction control
first passed gate 2's line for body k, "|error| 0.998, the bound 0.998":
the bound added |exact dataflow - reference|, which absorbed the missing
reaction. `06adb1a` gives every component a fixed 1e-150 allowance and
requires the exact dataflow to equal the reference within it.

### A format difference cannot see a narrowing both formats share

**Measured by P1, 2026-11-11, in an urgent file** (the ledger's
`urgent/P1-gate4-blind.md`). The plan's first control, the term rounded
to binary64 before it is added, leaves binary128 against binary256 on
the same 120 steps at 3.371e-33 of the distance from Saturn; the correct
build gives 3.402e-33, so gate 4 passes it. Rounded at binary128 only,
the difference is 4.029e-18 and gate 4 fails. The control was caught
instead by the gates that compare with a reference without the defect:
the golden model, gate 1, and mpmath, gate 2. V1's own copy of the
control gave 1.32e-33 against gate 4's 1.39e-30 at 120 steps.

Believed by P1, the mechanism, not measured apart: rounding to binary64
is a deterministic function, and the two runs' terms agree to about
1e-33, so both round to the same binary64 value and the defect cancels
in their difference. The same holds for any binary64 bottleneck two runs
of one build share, the drop-in's binary64 `additional_forces` among
them (README, "A force routine is evaluated at binary64").

**The lead's ruling, 2026-11-12 01:33** (`urgent/lead-gate4-ruling.md`).
Gate 4 stays, and its head states its scope: it can fail on a defect that
differs between the two formats, and cannot fail on a narrowing both
share. Its control is the binary128-only rounding. The plan's control is
kept, shown to fail gates 1 and 2, with its gate 4 number reported as the
measured blind spot. Both were run from copies, by P1 and again by V1, and
deleted: neither is a permanent control in the gate, whose head names
them. P4's format difference prints its scope beside
its number (entry 50).

GPU-E's section 6, item 3(b) proposes this check as "the check that
fails on exactly the defect the README names". It cannot. Believed here,
by message: cft-fp256's lead confirmed it and will correct the study
(entry 51). P1's urgent file reached the ledger 3 h 17 min after it was
written (entry 51).

### The verifier

**On `22ca3de`: READY.** V1 built `22ca3de` and `9a42b1d` from clean.
- **The force, derived by hand,** against `--zonal-selftest` at every
  format, decoded by V1's own decoder in mpmath: worst relative error
  6.7e-16, 2.5e-34 and 1.6e-71, reaction and potential included. Gate
  2's reference matches V1's closed form to 2.8e-201.
- **The count.** 70 `zcall()` sites and no other library call in the
  term; a hand count, 80 + 46 = 126 operations and 1,382 elements,
  equals the trace.
- **No FMA, measured.** An instrumented copy hashed x and a around
  every `zonal_accel()` call. The program engine in the FMA form equals
  the loop engine's at all 160, 349 and 510 calls, by format, and
  REBOUND's form equals the FMA form at all 62 calls whose state matched
  before the term.
- **The constants,** regenerated into scratch, reproduce the header and
  JSON byte for byte, and four copies with a coefficient moved each fail
  check 3.
- **Seven controls rebuilt from copies,** each failing where claimed, by
  9.6 to 150 decades over gate 2's bound at its least favourable point.
- **The untouched path.** 100 configurations identical between the
  commits, timing masked and `calls=` included. `check_records` PASS,
  9 of 9. `make check-light` from a clean export: rc 0, 20 legs, 671 ok,
  0 FAIL, no cache skips.
- **A limit, stated.** At binary64, dz^2/r^2 underflows when |dz| is
  below about 2e-149 km at the Encke gap, and the run aborts, "div raised
  0x18". Gravity's own window is 1.5e-154 km. Neither is reachable.

Its findings, all taken in `291c99c`:
1. `check_zonal.py:82` gave gate 4's bound as "eleven decades" under
   binary64's; u128/u64 = 2^-60, eighteen.
2. The cost: 34 mul, above.
3. **A seam in P0's parser, reachable once P1 lifted the refusal.** A
   second `zonal` line was accepted and replaced the first, so a problem
   written with two oblate bodies ran with one, rc 0, and a fifth field
   was ignored. A second `active` line was already refused.
4. The record's reaction count came from `zonal_ready()`, not from the
   loop that subtracts: with every body's reaction taken, the record
   still said 1.
5. `check_scalar`'s `bodies.get("zonal_ready", "")` would have dropped
   the term's constants from the recount, silently, had the function been
   renamed.

And `src/ias15_zonal.inc:16` and `:240` stated REBOUNDx's order flat;
only `:38` said believed.

### The lead's fix, `291c99c`, and the suites

- `read_problem()` refuses a second zonal line, a fifth zonal field and
  a second active field, by name, and `check_gravity` 2c holds the three,
  with nothing printed. P1's build ran all three files with rc 0, and with
  two zonal lines made body 1, the second line's, the oblate one: the
  cases fail without the fix.
- `zonal_gives()` is the one predicate behind the record's reaction
  count and `zonal_apply()`'s loop.
- The include's two sentences say believed; the selftest's R = 1.5 and
  `check_gravity`'s R = 1 are labelled synthetic; `check_zonal.py:82`
  says eighteen; `check_scalar` requires each `CONST_FUNCS` body by name.
- The refactor moves no bit: `--zonal-selftest` at every format, and
  zonal_saturn at 3 formats on 3 engine and form settings, are identical
  to P1's build with timing masked.

The suites, read from their logs:
- **P1's own, on `b3f664b`:** rc 0, 2 h 32 min on a shared host, 20
  legs in full, **977 `ok`, 0 FAIL**. Seventeen of the 19 existing legs
  had P0's counts; `check_corrector` +24 and `check_scalar` +9, from
  zonal_saturn. `22ca3de` changes only check_zonal.py's head and one
  verdict, and P1 reran that gate in full on it: 87 ok, PASS.
- **The lead's, on `291c99c`,** the tree that landed: 03:47:46 to
  06:34:31, 2 h 47 min on a shared host;
  - rc 0, 20 legs, each `--mode full` once;
  - **980 `ok`, 0 FAIL**, no fork errors, 12 `RESULT: PASS`;
  - leg by leg P0's counts, but `check_corrector` 302 (+24),
    `check_gravity` 44 (+3, the fix's refusals), `check_scalar` 99 (+9)
    and the new `check_zonal` 87;
  - the six digests of entry 45 are in the log.

### Stated limits, not closed

- **The recalled order is REBOUNDx's of 2020, not today's.** READ by study
  DEBRIS-A (2026-11-12, from REBOUNDx's repository on the web; nothing
  run against REBOUNDx): it matches gravitational_harmonics.c at 5406ae5
  (2020-04-13). The file was rewritten on 2024-10-20 (8402022, a general
  spin axis) with another operation order, so an operation-for-operation
  check against a current REBOUNDx would not match bit for bit
  (BELIEVED). The study did MEASURE the port at binary64 bit-identical,
  at the turn and on return, to REBOUND 5.1.1's IAS15 with J2 written in
  Python in that order: low Earth orbit in 100 s steps over a quarter day
  and a day, J4 zero.
- **The order is REBOUNDx's as recalled** (decision 4). The formula and
  the physics hold it either way; only an operation-for-operation binary64
  comparison with REBOUNDx would confirm the order.
- **The rates are measured at binary64 only.**
- **No J2 or J4 on the drop-in.** The library half compiles a stub, and
  ROUND4.md says why.
- **The constants leg does not check that src/ias15_constants.h is
  current.** `gen_constants.py --no-write` returns after its derivation
  checks, so a stale header passes every tier (P1). The lead's to
  decide.
- **A file in data/problems/ binds two gates.** check_scalar fails
  unless every file there is one of its cases, and check_corrector runs
  every one at binary64: on zonal_saturn its adaptive runs from dt 0.01 s
  grow to about 1,200 s, roughly an orbit, and its fixed runs reach t = 3
  s, exercising the replay and nothing physical (V1, measured; P1 had
  believed all of them non-physical).

One limit P1 stated is closed: the usage block at the head of
src/ias15_cft.c did not name `--zonal-selftest`, and the commits that
land these records add it there, beside the README's line.

### What this does not establish

- **That 60,330 km is the radius Cassini's J2 and J4 are referred to.**
- **Agreement with REBOUNDx,** at any format.
- **The card.**

## 48. Round 4's P2: the ring generator and its gate; Pan caps the step, and 10^6 particles in waves of 57 is 1,000,008; sent back once, for moons in Saturn's plane

**2026-11-11 and 2026-11-12.** Host and pins as in entry 47. P2 branched
from `9a42b1d`. Its commits are `94b0daa`, work in progress, `250789e`,
and `8bb49b9` after V2. `9dd4051` merges `8bb49b9` on `round4-int2`,
which with the lead's seam and P4 became `round4-int4`, at `7d50392`.

### What P2 is, in one paragraph

`tools/make_ring.py` writes Saturn, six moons and test particles, as one
system or as an ensemble of 64-body waves, 7 massive bodies then 57
particles: GPU-E's shape a2.
- **The constants** are read in GPU-E's section 4: Saturn's GM and
  radius from NASA's fact sheet, J2 and J4 from Cassini via Fortney et
  al., the moons' GM, a and e from JPL. R = 60,330 km is believed. The
  moons' inclinations, nodes, pericentres and phases are synthetic
  (decision 5), small and not zero.
- **The particles** come from Python's `random.Random(seed)`, inclined
  1e-4 to 1e-3 rad, sorted by radius and cut into waves. A start within 5
  Hill radii of a moon is drawn again, and counted in the file's head.
- **Exactness.** Everything is computed in mpmath at 160 bits, and each
  written value is rounded to binary64 once; the moons by REBOUND's own
  `reb_particle_from_orbit`; the frame is the massive bodies' barycentre.
- **Each file is circular in its own field:** with its
  `zonal 0 R J2 J4` line, at the oblate field's circular speed in the
  equatorial plane; with `--no-zonal`, a point mass's. `active 7` comes
  with every file.

docs/RING.md describes the ring and its generator in full.

### The gate, and the controls that show it can fail

`tools/check_ring.py`, one leg in each tier: the generator's bytes and
its file, one wave at every format on both engines, `active 7` against
every body active, REBOUND's own IAS15 on a wave, the ensemble against
its members, and the physics. From P2's full suite at `250789e`, 37 ok:

```
ok   the same arguments twice give the same 41859 bytes, and seed 2 gives other bytes (differing)
ok   fp64: the loop engine in REBOUND's form runs it, differing from the FMA form (first at sample 2, dt_next), its largest relative difference 1.698e-11 against a bound of 2^-26.5 = 1.1e-08 ...
ok   fp256 loop engine rebound form: `active 7` against all 64 active, 3 sample lines identical, the trailer too: 420 pairs a force evaluation against 2016
ok   fp64: w00001 as one system of 64, every body active: REBOUND's own IAS15 (ias15_ref) against ias15_cft in REBOUND's form at its defaults: 7 samples, 2716 values identical
ok   the moons act: with their masses zeroed the wave differs; after 6 steps a particle (p0000110) is 2.161e-02 km from where the moons put it, 1.6e-07 of its radius
ok   the far wave w00000 (radii to 132941 km, wholly inside pan's pericentre less 5 Hill radii, 133510 km): 35.89 steps an orbit of p0000000, its shortest (47845 s), over steps 12 to 24 at epsilon 1e-9, against GPU-E's MEASURED 36.2: -0.9% (the bound is 10%)
```

Over the gate's 23 runs every particle stays in the annulus and off
Saturn, the largest change of a radius 0.149 km.

Controls, each from a copy, deleted:
- **The plan's.** In the ensemble's copy, p0000143's x one ulp up:
  "FAIL w00002: first difference at sample 0, body 36 x", while w00000,
  w00001 and w00003 still match.
- **A flat ring,** P2's own: every inclination zero. The |z| counter
  fails, and 42 physics lines with it, while all 27 bit comparisons still
  pass: ROUND4.md's silent wrong answer, caught by the counter alone.
- **The moons flat and the particles inclined,** V2's, below.

### What P2 measured that the plan did not have

1. **Pan caps the step of every wave outside its orbit.** On a
   228-particle ring at binary64, REBOUND's form, epsilon 1e-9, waves
   w00002 and w00003, wholly outside Pan's orbit at 133,600 km, settle to
   exactly 1,388.09 s: Pan's period, 49,819 s, over 35.89. Their own
   innermost particles would give 1,423 and 1,471 s. w00003 takes 38.4
   steps per orbit of its innermost particle and 39.3 per its outermost.
   PRS23's minimum runs over every body, so a system's shortest orbit sets
   its step. Computed by P2 from a uniform radius, not run: about 3% more
   steps over the whole ring. GPU-E's cost model does not carry it. The
   far wave's own figure holds, 35.89 steps an orbit, 0.9% under 36.2.
2. **10^6 is not a multiple of 57.** `--particles 1000000` is refused by
   name, rc 1, and no file is written. 17,544 waves hold 1,000,008
   particles: `--particles 1000008` wrote 174,691,521 bytes, the particles
   from 130,000.025 to 141,999.987 km and the largest |z| from Saturn
   141.850 km, with 49 of 1,000,057 draws refused (Pan 2, Prometheus 19,
   Pandora 28). The file was read and deleted. An ensemble's systems all
   have N bodies, so a short last wave cannot be written here.
3. **The plan's "oblate speed" and "the gates use `--no-zonal`" cannot
   hold for one file.** Computed from a file with its zonal line, read in
   a point mass's field: each particle is on an orbit of e 0.0045 to
   0.0053, its apocentre 1,275 to 1,396 km out, the outermost reaching
   143,253 km, past the annulus. So each file is circular in its own
   field, and the zonal line is isolated by deleting it from a file, never
   by comparing with a `--no-zonal` one.
4. **The two forms differ in the step, not at round-off.** On w00001 at
   binary64, `--arith fma` and `--arith rebound` agree in dt_next to 7e-13
   after 2 steps and 1.7e-11 after 3, and positions then differ by 3e-12.
   Believed by P2: PRS23 builds the step from the high b coefficients,
   which amplify round-off. Records are compared like with like.
5. **The brief's third_party clone loop is refused by the worktree
   guard.** Plain, separate `git clone` and `git checkout` worked.

### The verifier, twice

**On `250789e`: SENT BACK.** V2 built `250789e` and `9a42b1d`;
`ias15_cft.exe` is byte-identical between them, and so are the 66
program images.
- **Defect: nothing held the moons out of the plane.** Control C, the
  generator with the six moons' inclinations set to 0 and the particles
  untouched, passed: 31 ok, PASS, rc 0. On such a ring every z term of J2
  and J4 on a moon stays exactly zero, and so does the z of its reaction
  on Saturn; the particles are massless, so the moons carry all of the
  ring's reactions. A wrong reaction z would match bit for bit, here and
  through P4's driver.
- **Notes.** The gate neither checked ias15_ref's `impl=` nor printed its
  build time, so a stand-in at that path would pass; PER_WAVE = 57 was
  typed, not derived; `--exclude-hill 1000` never finished.
- **Confirmed.**
  - Determinism: ten runs, one sha256, under another PYTHONHASHSEED,
    without gmpy and to another path. An independent 400-bit
    recomputation reproduced every value, 1,645 in the 228-particle ring
    and 399,049 in a 57,000-particle ring with its zonal line.
  - Every constant against GPU-E's section 4, and the orbits against
    REBOUND's tools.c:972-992, term for term.
  - The 420 pairs read from the record; 36.22 steps an orbit by the
    study's own method, the mean of its halves.
  - The exclusion at the default refuses 2 of 57,002 draws at 57,000
    particles, the same two V2's own replay of the draws refuses.
  - The untouched path: 72 runs a build, 315 sample lines, identical with
    timing masked. `check_records` PASS. `make check-light` rc 0.

**On `8bb49b9`: READY.** The gate asserts each moon's sin i, from r x v
relative to Saturn in the written file, above zero: not |z|, since Pan
starts at its node with z = 0. ias15_ref's record must say
`impl=rebound`, and both binaries' build times are printed. PER_WAVE is
derived as WAVE - MASSIVE, the file's N and active are held to the
generator's, and the exclusion's redraws are bounded.
- `git diff --stat 250789e 8bb49b9` is the two tools only, and the
  generator's bytes are `250789e`'s on five rings, one of them 57,000
  particles with zonal lines.
- Control C fails by name: 33 ok, 1 FAIL, "the moons are not flat: ...
  0.0000 degrees". V2's own control, Pan flat and the other five
  inclined, fails by name too, so the check is moon by moon.
- A copy of `ias15_cft.exe` at ias15_ref's path passes the bit
  comparison, 1,940 values, and fails the identity line,
  "program=cft impl=libcft-0.18".
- `--exclude-hill 1000` stops after 6,644 refused draws, rc 1, no file.
- The gate in full: 40 ok, PASS; with `--quick`, 34 ok.

### No record line moved, and the suite

- **P2's own full `make check`, on `250789e`:** rc 0, 2 h 36 min on a
  shared host, 20 legs in full, **894 `ok`**, main's 857 plus the leg's
  37, 0 FAIL, no fork errors, and every other verdict line present.
- **The merged tree's suite** is entry 51's.

### Stated limits, not closed

- **One host.** mpmath's two backends give the same bytes here. The same
  file on every host is believed; the check across hosts waits for a WSL
  slot.
- **The wide-format form comparisons are two runs of one build.** At
  binary128 and binary256 the two forms agree to 4.3e-30 and 1.9e-68.
  They share gravity, `pc_error` and the step control, so a binary64
  bottleneck there would cancel, as in entry 47 (believed by P2). They
  show the forms agree, not that a ring runs at full width. At binary64,
  REBOUND's own IAS15 is an independent build, and agrees bit for bit.
- **The ensemble comparison runs at binary64,** on the program engine.
  check_ensemble covers ensembles at every format.
- **No rejected step and no close encounter** in the gate's runs, which
  are 2 to 24 steps.
- **mpmath is not a gate-cache input.**
- **"The engine computes all 2,016 pairs"** with every body active is
  inferred from the code.
- **check_ring runs `--no-zonal` rings.** The ring with its zonal line is
  the seam test's (entry 51).

### What this does not establish

- **The ring at the study's size, integrated.** It was written and read,
  not run.
- **An epoch's moons** (decision 5).

## 49. Round 4's P3: a test particle replayed out of its wave, rejected attempts included; only a rejection before the first accepted step leaves state; sent back once, for a dt file through a pipe

**2026-11-11 and 2026-11-12.** Host and pins as in entry 47. P3 branched
from `9a42b1d`. Its commits are `60457ce`, `569ee4e` and, after V3,
`2638e72`. `7ece3e3` merges `2638e72` on `round4-int5`, after P2 and P4.
A first merge of P3, `3b69eae` on `round4-int3`, was set aside when the
merge order changed (entry 51); `7ece3e3`'s src is identical to it.

### What P3 is, in one paragraph

- **`--attempts-out FILE`** writes one system's every step attempt:
  `# cft-rebound attempts v1`, `# run` and the run's key, a legend, a
  line an attempt (its dt as exact hex, then `a` or `r`), and last the dt
  the run takes next. An ensemble is refused by name.
- **`--dt-file`** reads that form by its first line and replays it. The
  recorded outcome and next dt take `choose_timestep()`'s place, so a
  recorded rejection restores the step's start and predicts from `er`
  and `br` as the run did. An accepted-steps dt file reads and runs as
  before.
- **The run key** hashes the outcomes with the dts under that form, and
  names it `dt_file=attempts`. No older key moves.
- **`--run-key`** prints the key and exits. It reads the problem and any
  dt file, and opens no output file.
- **`tools/extract_particle.py`** writes one test particle alone with its
  wave's G, active and zonal lines and first N_active bodies, verbatim and
  in order, after asserting by the program's own key that the attempts
  record and the schedule are this wave's. Then it re-keys the wave's
  corrector schedule for the sub-system. Not a key that leaves the test
  particles out: every wave of a ring starts its massive bodies from one
  state, so such a key could not tell one wave's schedule from another's.
  No Python hashes anything, so no copy of the key can drift from
  `sched_key()`.

### What the plan had wrong: a rejection after an accepted step leaves no state

ROUND4.md said the replay cannot carry a rejected step, since a
rejection changes the next step's prediction. Measured by P3 on
2026-11-11, then by P4 and V3:
- **The mechanism, read and then run.** After an accepted step `er` and
  `br` are copied from the corrected e and b before that step's
  prediction. A later rejection's prediction, from `er` and `br` at the
  ratio of dt_next to dt_last_done, is the same call on the same inputs
  the accepted step's own prediction would have made. Before the first
  accepted step dt_last_done is 0, nothing is predicted, and the rejected
  attempt's corrected b is the next attempt's first guess. That is the
  only state a rejection leaves.
- **The old accepted-steps replay,** its schedule with the rejections
  after accepted steps deleted and re-keyed: 0 of 7 record lines differ,
  at binary64 and binary128. With one rejection as the first attempt (a
  ring wave from dt0 1e5 s, kepler from 1.0, kepler_e099 from 0.01): 5 of
  5, 7 of 7 and 7 of 7 lines differ.
- **The one token it moves** is the dt_next of a sample printed just
  before the rejected attempt, which carries that attempt's dt. P4
  measured it on its member 8 at binary64 and binary128 (steps 4, 8 and
  10), and V3 on its own wave at binary256 and binary64.
- **A rejected attempt's own dt moves no bit either.** One ulp on each of
  five rejections after accepted steps changed 0 of 48 comparisons, where
  one ulp on an accepted attempt changed 19 and 21.

So the attempts form is needed for waves that reject before their first
accepted step, and for the counters. The plan's fallback, rerunning every
wave that rejects, narrows to those waves.

### The gate, and the controls that show it can fail

`tools/check_replay.py`, one leg in each tier. Three hand-made waves,
Saturn, six moons and eight test particles under `active 7`, every value
one correctly rounded binary64 operation from exact decimals: A rejects
nothing; B starts from a dt0 far too large and rejects its first attempt;
C adds a particle aimed past Janus, which rejects after accepted steps.
Each runs at every format, and B and C on the loop engine in both forms
and on the program engine too. Three particles a wave are replayed
alone, each held to the wave's bits for itself and all seven massive
bodies, every sample's head and exact-time pair, and five counters. The
massive bodies are held because they match only if the particle adds
exactly nothing to them. From P3's full suite, 276 ok at `569ee4e`:

```
ok   a rejection before any step is accepted: 1
ok   body 7 replayed alone: 5 samples, the particle and all 7 massive bodies, 240 values, and every sample's step, t, dt_next, dt_last, E and exact-time pair are the wave's
ok     and its counters are the wave's: steps_done=40 steps_rejected=1 iterations_max_exceeded=0 mean_pc_iterations=3.293 max_pc_iterations=6
ok   attempt 2, a rejection after an accepted step, deleted from both records and re-keyed: every bit is the wave's (0 mismatches) and only the counters see it, steps_rejected 4 against 5
```

Sabotages, each built by P3 from a copy and deleted:

| sabotage | with it |
|---|---|
| one accepted attempt's dt moved one ulp, re-keyed | FAIL at every format, all three of wave A's particles at 26 mismatches each; not re-keyed, refused by name |
| the rejected first attempt deleted | the bits differ, 37 mismatches a particle at every format |
| a rejection after an accepted step deleted | every bit matches, and only the counters fail: steps_rejected 4 against 5, and at binary256 max_pc_iterations 29 against 42 |
| after a recorded rejection, the replay predicts as after an accepted step (a source copy) | FAIL at every format, engine and form: wave B, 34 to 36 mismatches of 40; wave C at binary256, 21 to 22; wave C at binary64 and binary128 in REBOUND's form, only the flyby particle, 5 |

**Bits alone cannot hold the passes:** the sub-system's own exit test
took other passes in up to 33 of 64 attempts with identical bits, so the
gate holds the counters too. Run outside the gate at binary64: a 64-body
wave, 4 particles at 240 values each, 0 mismatches; the flyby wave over
400 steps with 5 rejections, 3 particles at 432 values, 0 mismatches.

The last row matters for the GPU half. After a mid-run rejection at
binary64 and binary128 in REBOUND's form, a wrong prediction showed in
the particle that caused the rejection alone. So a driver sampling a rejecting wave
should take that particle (P3); `compare_run` takes each member's least
closest approach (entry 51).

### The verifier, twice

**On `569ee4e`: SENT BACK.**
- **D1, a defect.** `read_dt_file()` read the first line, then called
  `rewind()` without checking it. A pipe cannot rewind, and the C
  runtime drops its 4,096-byte buffer. Through a Windows named pipe,
  `9a42b1d` read old-form dt files exactly as from a file; `569ee4e`
  refused kepler's ("has 0 steps") and outer's, and ran a file padded to
  32-byte lines with rc 0 from its line 129, all 6 sample lines wrong
  and nothing printed. No gate covered it.
- **F2, a false comment.** The head said a replay with `--replay-iters`
  of the run's own schedule is that run, bit for bit. That schedule is
  refused, rc 3, its key differing in epsilon and `dt_file=attempts`.
- **Confirmed.** On V3's own waves, the 64-body one included, every
  particle alone: 6,528 values, 0 mismatches. 144 of 144 base schedules
  replay unchanged; no older key moves in 72 configurations; `--dt-out`
  is byte-identical in 72; check_corrector's leg is P0's. The
  massive-body comparison is load bearing: with test particles pulling
  on the massive bodies, only Saturn differs at binary64. Both controls
  re-run and fail; the untouched path holds; `make check-light` rc 0.
- **A note.** Wave C's deletion case held only because no sample fell
  just before the deleted attempt.

`2638e72` reads the first line by the same `fgets()` the accepted-steps
loop uses, and starts the loop from it: nothing is read twice or rewound.
A new case feeds wave A's `--dt-out`, padded to 385 lines of 32 bytes,
12,320 bytes, and its attempts record through a named pipe
(`/dev/stdin` on a POSIX host). F2's sentence is corrected and gated:
every wave, at every format and form, replays whole through its attempts
record with no schedule. C's precondition is asserted by name, and where
a rejection does have a sample just before it, only that line's dt_next
may move. B's deletion must move the bodies' state, 32 of 40 pairs, not
only a dt_next.

**On `2638e72`: READY.** Through a named pipe, kepler, outer,
kepler_e099 and a padded outer read identically at `9a42b1d` and
`2638e72`, and an attempts record differs in 0 of 61 lines. From regular
files, 31 configurations are identical, and 62 of 62 base-recorded
schedules replay to the base's run. The pipe case fails on a build with
the old reader: 118 ok, 1 FAIL. Nothing is loosened. The gate in full:
303 ok, PASS; with `--quick`, 119.

### No record line moved, and the suite

- **P3's own full `make check`, on `569ee4e`:** 22:52:32 to 01:29:13,
  2 h 36 min on a shared host; rc 0, 20 legs in full, **1,133 `ok`**,
  P0's 857 plus the leg's 276, 0 FAIL, no fork errors; every existing leg
  at P0's count, and all 21 wide-state digest lines P0's.
- **The merged tree's suite** is entry 51's.

### Stated limits, not closed

- **The gate compares two runs of one build.** A defect in code both
  share, gravity, the corrector or the prediction itself, would cancel,
  as entry 47's format difference does. Agreement with REBOUND is
  check_equivalence's and check_dropin's.
- **The POSIX pipe branch has not run** on this host. Believed.
- **A wider format's record replayed at a narrower one is rounded, not
  refused,** as the older dt file is.
- **The replay's `flags_seen` and `calls` are not compared.** The replay
  runs neither the step control nor the exit test.
- **`ATT_MAGIC` is a hand copy in the tool.** If it drifts, every
  extraction is refused, loudly.
- **The zonal path is the seams'.** P3's waves have no zonal line, and
  extract_particle refuses one whose body is not massive, since the
  sub-system keeps only the indices below N_active (entry 51, part 2).

### What this does not establish

- **A replay from another producer's records.** Every attempts record
  and schedule here is this tree's own; taking a GPU run's is not built
  (entry 51).

## 50. Round 4's P4: `tools/compare_run.py` holds a run's records against this tree, by bits; resumed three times, twice because a driver comparing less passed the gate

**2026-11-11 and 2026-11-12.** Host and pins as in entry 47. P4 branched
from `9a42b1d`. Its commits are `c640f65`; `517e8d5`, after the lead's
ruling on entry 47's finding; and `f27f9c6` and `79bd206`, after V4.
`801e2be` merges `79bd206` on `round4-int4`, and `7d50392` is the lead's
hardening of control a5, from V4's last note.

### What P4 is, in one paragraph

`compare_run.py RECORD --problem FILE --format F --arith A` reads only
the record format, the lines `ias15_cft`'s `main()` writes, never the
program's internals, so an engine that writes the same records is held
with no change to the driver.
- **Refused, rc 2, by name, before anything is compared:** a file that is
  not a record, or whose header lacks a field the rerun is built from; a
  record of another format or arithmetic than `--format` and `--arith`
  name; a problem file that is not the record's problem, held by name, G,
  systems, every body's name and mass, the `active` line and every
  system's sample 0; a record run from a dt file; and a rerun whose own
  header is not the record's run.
- **FAIL, rc 1, by member:** a member missing or incomplete, every
  member checked; a sampled member's first differing line, by the
  record's line number, sample, step and column; its five trailer
  counters; and the members' rejections against the record's total.
- **The sample:** strata of member index, which is radius in P2's waves,
  one member from each by the least SHA-256 of seed and index, so a seed
  picks the same members on every host, and every member whose trailer
  shows a rejection.
- **Beside each verdict:** steps, rejections, passes and cap hits; each
  test particle's closest approach to each active body at the samples;
  and with `--format-diff`, the difference at a second format on the
  same accepted steps, by `tools/compare_formats.py`.
- **Parsing is imported,** `check_ensemble.parse()` and
  `check_records.header_args()`, neither changed.
- **The particles' place** was named and empty, with a tripwire: once
  P3's tools were in the tree the driver said NOT WIRED and its gate
  failed, in whichever order the two merged, until the lead's P3 x P4
  seam filled it (entry 51).

### The gate, and the controls that show it can fail

`tools/check_compare_run.py`, one leg in each tier. Its ensemble is made
by hand, because `make_ensemble.py` refuses `active` lines: 9 waves, each
a primary, a moon and three massless test particles inclined about 1.8
degrees, under `active 2`, every value a chain of correctly rounded
binary64 operations. From dt0 0.3, members 0 and 1 reject their first
attempt and the rest accept it. Member 8, appended at `517e8d5`, has a
retrograde particle that meets the moon head on near t = 1, so its 3
rejections all follow accepted steps; members 0 to 7 kept their bits.
Every control is permanent, asserted to FAIL or be refused by name:

| control | what it plants |
|---|---|
| a, a2 | one byte of a sampled member's line, in an ensemble's record and in a one-system one |
| a3 | the sign of a zero, a sample-0 dt_last made -0x0p+0, binary64 |
| a4 | one hex digit past binary64 in a binary128 record, its binary64 value asserted unmoved |
| a5 | one change in each column the record's own `# columns:` line names: t, dt_next, dt_last, E, the last body's six fields, one field of every other body, and t_hi and t_lo |
| b | one byte of a trailer counter |
| c, d | a member deleted; the record cut short |
| e, f | a record of another format, of another arithmetic, and REBOUND's own ias15_ref record, which shares the first line |
| g | four problem files under the record's name: a test particle's mass 0 made 2^-40, the active line dropped, a coordinate one ulp off, the name changed. The first two move no bit |
| h | the `# active` line with one pair too many |
| i | an empty stand-in for each of P3's two files: at `79bd206` the driver had to say NOT WIRED and exit 1; since the lead's wiring (`e5636a9`, entry 51) it must say P3's replay cannot be called, and FAIL |

P4 ran one-defect copies of the driver through their own gate: 12 at
`c640f65`, 14 at `517e8d5`, 18 at `f27f9c6` and 21 at `79bd206`, each
caught, the unmutated copy passing. With the problem check removed, a
problem with a test particle's mass made 2^-40, or one initial x one ulp
off, in a member the sample does not rerun, passed with rc 0: only the
problem check refuses those. Case 1 at binary64, from `517e8d5`:

```
ok   fp64: the sample [0, 1, 2, 3, 6, 8] holds every rejecting member [0, 1, 8] and one member of each stratum [(0, 2), (3, 5), (6, 8)], and is 6 of 9
ok   fp64: member 8, whose 3 rejected attempts all followed an accepted step, replayed from its dt file at fp64: every value but dt_next is the record's at all 3 samples (dt_next differs at 0), so they leave no state bit, as its note says
```

**The scope beside the number, after the ruling (`517e8d5`).** Every
`--format-diff` line ends with the ruling's words, "the difference shows
the narrower format's round-off on these steps; a narrowing both runs
share passes it", and the gate states the text itself rather than read
it from the driver. The note on a rejecting member follows whether its
first attempt was rejected: the dt file's first step is the record's dt0
exactly when that attempt was accepted. V4 read why that cannot misfire
(src/ias15_cft.c:2545-2553, :2964 and :3001-3002) and measured it against
a one-step run's rejection count: 17 of 17 members at binary64, 17 of 17
at binary128 and 9 of 9 at binary256.

### The verifier, three times

**On `517e8d5`: SENT BACK.**
- **D: the gate could not tell bits from binary64 values.** A copy
  comparing `float(a[0]) != float(b[0])`, binary64 values without their
  sign, passed the whole gate in full, 54 ok. It also passed the binary64
  record whose sample-0 dt_last is -0 and a binary128 record one hex
  digit off below binary64, both of which the shipped driver fails. Every
  byte-change control was on a binary64 record, and none flipped a zero's
  sign.
- **Notes.** NOT WIRED exited 0 with RESULT: PASS. P3's file names were
  not inputs of the leg's cache key, so check-quick and check-light
  would skip the leg when only P3's files changed. A comment's rounding
  count was wrong, and the pair-count refusal named the wrong field. The
  sample trusts the trailers it tests, a limit of the briefed design.
- **Confirmed.** Five planted differences each FAIL by name. A record
  rewritten as another producer's, `program=gpu impl=hip-6.2
  backend=hip:0` with one extra header and one extra trailer line, passes,
  so the driver reads only the record format. The sampling, recomputed
  independently; members 0 to 7 byte-identical to `c640f65`'s; the
  untouched path, 144 of 144; `check_records` PASS; `make check-light`
  rc 0, cold, all 20 legs run.

`f27f9c6` adds a3 and a4, makes NOT WIRED end RESULT: FAIL and exit 1
(control i), puts P3's two files in the leg's key in all three tiers,
prints the sample's limit beside the sample, and fixes the text.

**On `f27f9c6`: SENT BACK.** V4's float() mutant now fails through a3 and
a4, and a sign-only mutant through a3 alone. **Defect:** a copy skipping
the exact-time pair, `zip(rv[:-2], mv[:-2])`, which is what
`check_equivalence.parse_record` does, passed the full gate, 63 ok. The
shipped driver fails a planted t_hi digit, and the mutant passes it with
rc 0. That pair comes from `vaugadd`, the path G1B is expected to refuse
as AUGADD on hip:N. [Corrected in round 5: G1B does not refuse it.
`vaugadd` is the API call `cft_augmented_add`, which at the pin runs on
the host on every device, and every `hip:0` record measured holds the
pair at the software backend's bits (entries 54 and 55).]

`79bd206` adds a5: it reads the record's own `# columns:` line and plants
one change in each of its 12 columns, one column a sampled member, over
2 records. It also makes `--exe` absolute: a relative path written with
forward slashes passed the file check, and Windows' CreateProcess could
not run it.

**On `79bd206`: READY.** a5's columns are derived: a `# columns:` line
with three invented columns gave 15, unedited. The pair-skipping mutant
fails through a5 alone, and so does V4's own, which skips the
velocities. **Not blocking:** a5 planted body fields in the last body
only, and a copy skipping body 0's six values passed, 42 ok. Body 0 is
the primary, whose last bits the zonal reaction can move alone.

The lead's `7d50392` plants one field of every other body as well, the
fields in turn, 16 columns over 3 records. A copy skipping body 0 fails
a5 by name, "NOT: x of body 0 (primary) in member 0: IDENTICAL". The gate
in full at `79bd206`: 64 ok, PASS; with `--quick` at `7d50392`, 42 ok.

### The suite

**P4's own full `make check`, on `c640f65`:** 22:50:46 to 01:27:14, 2 h
36 min on a shared host; rc 0, 20 legs in full, **907 `ok`**, P0's 857
plus the leg's 50, 0 FAIL, each existing leg at P0's count. The merged
tree's suite is entry 51's.

### Stated limits, not closed

- **What a record does not carry.** `# program=` has format, backend,
  problem, E, N, steps, sample, cs, arith, engine, pc_tol_shift,
  max_iter, member and system, and dt_file; `# dt0=` has dt0, epsilon, G
  and pc_tol. Nothing names `--adaptive-mode`, `--softening`, `--min-dt`
  or `--no-flag-abort`, though the corrector schedule's run key names the
  first three (P4, read in `main()`). The rerun takes their defaults, so a
  run made otherwise passes only where the setting moved no bit.
  ias15_ref's record prints adaptive_mode=2; ias15_cft's does not.
- **A record's first line does not name its producer:** ias15_ref
  writes `# cft-rebound ias15 record v1` too.
- **The sample trusts the trailers it tests,** measured by V4: a member
  whose trailer under-reports its rejections, the record's total
  adjusted to match, leaves the sample, and the record passes. An
  unsampled member's lines are never rerun. The driver prints this
  beside the sample.
- **Not built:** sampling by a close approach to a moon, for which no
  threshold was defined, the closest approach staying a counter; steps
  and passes an orbit, since a record carries no period and a wave's step
  is its shortest orbit's (entry 48); both format pairs in one call.
- **The format difference is a counter, not a detector** (entry 47).
- **A record whose sample 0 is not the problem file's** is refused as
  the wrong problem, rc 2, not reported as a FAIL (V4).
- **`check_records.header_args()` returns nothing for ias15_ref's
  unspaced `# dt0=` line;** the driver falls back to its own reading.

### What this does not establish

- **A GPU's record.** Every record the driver has held is this tree's,
  or one of them rewritten by hand as another producer's. The GPU seam
  test, a difference planted in the first GPU record and shown caught, is
  the GPU half's (entry 51).

## 51. Round 4's integration: three seams, the P3 x P4 wiring, a verifier of the lead's own commits that sent them back twice, the merges and suites, and two lapses

**2026-11-11 and 2026-11-12.** Host and pins as in entry 47. The lead's
commits: `4b52705` and `9a42b1d` (entry 46); `291c99c` (entry 47); the
seams, `ecd2b20`, `110834f` and `acf9fff`; the wiring, `e5636a9`;
`7d50392` (entry 50); and, after VL, `7a339d1` and `d9c8c8c`. The merges
are `4dd53f3` (P1), `9dd4051` (P2), `801e2be` (P4) and `7ece3e3` (P3).

### The merges

Logan delegated the round's merges and pushes to the lead on 2026-11-11,
on round 3's bar: ROUND4.md's decision 1. The bar is a parcel's READY,
then a full `make check` on the merged tree, read from its log, before
main moves. Each merge is `--no-ff`, on an integration branch:
- **P1** on `round4-int`: `4dd53f3`, no conflict, then the fix
  `291c99c`. main fast-forwarded `9a42b1d..291c99c` and was pushed after
  its suite (entry 47).
- **P2 and P4 under one suite,** as ROUND4.md allows two parcels that
  share no file beyond a leg each: `9dd4051`, no conflict, and seam part
  1 on `round4-int2`; then `801e2be`, the Makefile's four conflicts each
  kept both ways, and `7d50392` on `round4-int4`. `round4-int`
  fast-forwarded to `7d50392`.
- **P3** after them, on `round4-int5`: `7ece3e3`, one conflict, P1's
  `do_zonal` beside P3's `do_key` in `main()`, both kept; `acf9fff`, seam
  part 2, cherry-picked from the first P3 merge's branch, a Makefile
  conflict resolved; then the wiring and VL's fixes.

The order changed at 05:38 on 2026-11-12: P1, then P2 and P4 under one
suite, then P3 with the seams and the wiring. At 08:56 the records were
given a suite of their own after P3's, so that P3 landed sooner: four
suites after P0's. Each parcel added a leg to all three tiers, and the seams one: 19
at P0, then check_zonal, check_ring, check_ring_seams, check_compare_run
and check_replay, 24 at `d9c8c8c`, counted from the Makefile's `--leg`
lines. VL measured each merge's tiers to be the union of its parents',
the three listing the same legs in the same order.

### The seams, `tools/check_ring_seams.py`

One leg in each tier. Each part holds two parcels together where
neither's gate looks.

**Part 1, the ring in its own field (P1 x P2),** `ecd2b20` and `110834f`.
A ring generated with its zonal line, against the same file with the
line deleted:
- **without the line,** each particle starts too fast by
  delta = 3/2 J2 (R/r)^2 - 15/8 J4 (R/r)^4, so it moves outward on an
  orbit of eccentricity delta: in full mode 1,169 to 1,233 km over 24
  steps, within 2.0% of delta r (1 - cos M), against a bound of 10%;
- **with the line,** P1's field holds each particle's distance from
  Saturn to 0.18 km, 1.5e-4 of that excursion, against a bound of 1e-3
  (3.5e-5 over the 12 steps of `--quick`);
- **the control, permanent:** J4 zeroed in the zonal line moves each
  particle outward 16.6 to 18.1 km, within 1.8% of J4's share of delta,
  18 times the bound, so the bound resolves J4;
- **the nodes,** ROUND4.md's precession counter: over about ten orbits
  each particle's node, from h = r x v relative to Saturn, regresses at
  n - nu, first-order theory in J2 and J4, to 9.0e-4 at worst against
  5e-3, and the same theory without J4 misses by 2.8e-2 at least, 31
  times as far, so the run resolves J4's vertical force. In `--quick`,
  8.3e-3 against 2e-2 over three orbits, J4 not claimed;
- **the energy,** the zonal potential included, moves 1.2e-15 over 360
  steps;
- **the bits,** with the line: the two engines agree at every format,
  `active 7` gives every body active's samples (420 pairs against none, a
  reaction of 6 against 63), a two-wave ensemble gives each wave alone,
  and the record's `# zonal` line is the file's, its counts derived from N
  and active.

**Part 2, a ring particle replayed with its zonal line (P2 x P3, P1 x
P3),** `acf9fff`. P2's waves with P1's line; three particles a wave
replayed alone through P3's `extract()`, imported, from a first step of
100 s and from one of 1e5 s, whose first attempt is rejected, asserted.
Fifteen replays at binary64, binary128 and binary256, 192 to 432 values
each, give the wave's bits for the particle, all seven massive bodies,
every head and exact-time pair, and the counters, with the record naming
the term felt by 7 bodies and 6 reacting. The control, permanent: a
sub-system with its zonal line deleted and its schedule re-keyed by hand
through `--run-key` replays to other bits in 64 of 72 sample-and-body
pairs, so the line is load bearing in the replay.

**Part 3, the driver on the ring's own shape (P2 x P4, P1 x P4),**
`7a339d1`, on VL's finding that no gate ran it. A two-wave ensemble of
P2's ring with its zonal line, recorded from a first step far too large,
each wave rejecting at least once: `compare_run` holds both members
IDENTICAL, and 4 particles, replayed with the zonal line, IDENTICAL. The
control, permanent: one hex digit of a moon's x on wave 0's last line
fails wave 0 and both its particles' replays by name, while wave 1 stays
IDENTICAL.

`check_ring_seams --quick` at `7a339d1`: 21 ok, PASS. VL derived part
1's excursion and its n and nu independently, and its own copies of
`ias15_zonal.inc` each fail part 1 in full: J2's x and y factor times 1.01
fails the circular check, 1.17e-2 against 1e-3; J4's z factor zeroed
passes it and fails the node counter, 3.72e-2 against 5e-3; the z factor
negated makes the nodes advance. The node fit's residual is short-period
noise, 31 times under J4's signal.

### The P3 x P4 wiring, `e5636a9`

P4 left `compare_run`'s `sample_particles()` empty and failing once
P3's tools were in the tree. Now:
- each rerun member's rerun also writes its attempts record and
  schedule, which move no record line (VL: 30 of 30 configurations);
- two test particles of each rerun member are sampled: its least closest
  approach to an active body at the samples, since a rejection's path can
  show in that particle alone (entry 49), and one by the seed;
- each is written alone with the member's massive bodies and its zonal
  line by P3's `extract()`, replayed, and held to the record's own lines,
  not the rerun's: every head field, the particle and every massive body,
  the exact-time pair, and the trailer's counters;
- P3's tools present and not callable FAIL by name.

`check_compare_run` asserts the picks, the least being the body the
gate's own mpmath distances name, and gains control p: a digit planted
in a sampled particle's x fails exactly that particle's replay, and one
in the primary's x fails both. In full: 66 ok, PASS, 12 particles a
format replayed alone and IDENTICAL to the record.

### VL, a verifier of the lead's own commits

No parcel's verifier saw the fix after V1, the seams, the a5 hardening,
the wiring or the merges, so VL was dispatched on them at 06:00 on
2026-11-12. It did not block the merges, whose bar is a parcel's READY
and the suite; its findings were to be fixed before the round closed.

**On `e5636a9`: SENT BACK.**
- **Defect: the particle path was held in two values.** A copy whose
  `particle_compare` checked only the particle's x and the primary's x
  passed `check_compare_run --quick`, 44 ok, while still printing "N
  values and the trailer's counters": control p planted x columns only.
  The hole V4 sent P4 back for, in the lead's code.
- **A combination no gate ran:** `compare_run` on P2's ring with P1's
  zonal line, the shape a GPU record will have. By hand it passed, and a
  digit planted in Mimas's x failed the member and both its particles.
- **Four texts:** check_ring_seams' header gave the `--quick` residual as
  the full one, and credited J4 to the 10% check, where J4's share of
  delta is 1.4%; compare_run's header had two stale lines about P1's
  zonal line; `read_problem()`'s comment did not say an inline comment is
  refused on the active and zonal lines.
- **Confirmed:** `291c99c`'s refusals, with 24 of 24 valid files still
  running in CRLF and with trailing blanks, and `zonal_gives()` moving no
  bit, 12 of 12 and 9 of 9 configurations; part 1's physics and controls,
  above; part 2 by hand at binary128, 406 of 406 tokens, its control
  biting through the bits; a5, a skip-body-1 copy failing; the wiring's
  picks, and its records unmoved; the four merges, each resolution
  keep-both with nothing lost; the untouched path, 60 of 60 identical;
  `make check-light` cold, 24 legs, 890 ok, 0 FAIL.

`7a339d1` adds p2: a planted change in every column `particle_compare`
compares, read from `# columns:`, 24 a particle, and a trailer counter,
each failing the least-closest particle's replay by column, 25 held.
Three mutants of `particle_compare` each fail p2 alone, on exactly what
they drop. And part 3, and the four texts.

**On `7a339d1`: SENT BACK.** p2's plan is derived, 29 columns from a
line with invented ones, and the two-x mutant fails p2 alone. **Defect:**
a `particle_compare` comparing every column on the last sample line only
passed, 45 ok: control p and every p2 plant sat on the last line. It
would have missed a replay that differs only in an intermediate line's
dt_next, the one difference entry 49 measured can occur alone.

`d9c8c8c`: p2's plants cycle over each member's sample lines after sample
0, which the problem check refuses on, and the gate asserts that more
than one line is planted; three non-structural trailer counters are
planted, one a member. The lead's two line mutants, applied: the last
line only fails p2 alone, 15 of 27 held; one middle line only fails p2
and control p.

**On `d9c8c8c`: READY.** The tree: 45 ok, p2 27 held on lines 1 and 2.
VL's last-line mutant fails p2 alone, 15 held, and its new one, the
particle's fields on every line and the rest on the last line only, 18
held. p2 is neither inert nor self-confirming: its planting reads
`# columns:` apart from the driver's reader, each check names its line,
its column and both values, line coverage is enforced, a shifted pick
fails, and the trailer plants are isolated.

Every commit of the round then had a verifier's READY: P0 (V0), P1 (V1),
P2 (V2), P3 (V3), P4 (V4), and the lead's (VL).

### The suites

Read from their logs, each on a shared host:

| commit | tree | legs | result |
|---|---|---|---|
| `4b52705` | P0's first | 19 | rc 2 at leg 7: MSYS fork failures (entry 46) |
| `9a42b1d` | P0 | 19 | rc 0, 857 `ok`, 0 FAIL, 2 h 20 min (entry 46) |
| `291c99c` | P1 and the fix | 20 | rc 0, 980 `ok`, 0 FAIL, 2 h 47 min (entry 47) |
| `7d50392` | P2 and P4 with seam part 1 and the a5 hardening | 23 | rc 0, 1097 `ok`, 0 FAIL, 2 h 34 min (06:35:43 to 09:09:48); each of the 20 legs both ran has `291c99c`'s ok count, and the new check_ring 40, check_compare_run 64, check_ring_seams 13; entry 45's six wide-state digests present; main fast-forwarded 291c99c..7d50392 and pushed |
| `d9c8c8c` | P3, seam parts 2 and 3, the P3 x P4 wiring and VL's fixes | 24 | rc 0, 1424 `ok`, 0 FAIL, 2 h 37 min (09:11:20 to 11:47:58); of the 23 legs both ran, 21 have `7d50392`'s ok count, check_compare_run 67 (+3, the wiring's controls) and check_ring_seams 34 (parts 2 and 3 joined part 1); the new check_replay 303; entry 45's six wide-state digests present; main fast-forwarded 7d50392..d9c8c8c and pushed |

The commit that lands these records ran its own full make check before
its push. A commit cannot record its own suite, so that result is in the
round's closing report; the ledger archived with
these records (`archive/round4-ledger.zip`) ends before it.

### The lead's two lapses

**The drop went unwatched from 21:20 to 01:31.** Agents in isolated
worktrees cannot write `.ledger/`, so their entries go to a drop the lead
places.
At dispatch, 21:20 on 2026-11-11, the lead armed no watcher on it. P1's
urgent file, complete at 22:14, and every drop entry before 01:31
reached the ledger at 01:31: P1's finding 3 h 17 min late, the earliest
entries up to 4 h, P3's and P4's addressed to each other and P1's to the
other three among them. P2's, P3's and P4's full suites had all ended,
01:27 to 01:29, before it was placed, and P4's format difference took its
scope in a follow-up, `517e8d5`. From 01:33 a watcher places every drop
file every 60 s and wakes the lead on any urgent one.

**A mutant run that ran no mutant.** At 08:05 on 2026-11-12, testing the
fix for VL's line defect, the lead's "last line" and "first line"
mutants of `particle_compare` both passed. They had never been mutated:
a heredoc edit of the mutation script put a raw newline in a string, the
script died on a SyntaxError, and the runner went on and ran the
unmutated copies. A stand-in passing for a control, the shape this round
hunted. The script now exits non-zero unless its replacement applies
exactly once, and the runner writes "=== NOT MUTATED" and runs nothing
when it does not. Re-run, both mutants fail (above).

### What went to cft-fp256

At Logan's word, one status message at 05:55 on 2026-11-12, before VL's
report: P0 pushed, P1 to P4 READY; GPU-E's section 6, item 3(b) cannot
fail on a shared binary64 bottleneck (entry 47); the term costs a body 3
sub, 34 mul, 14 add, 1 sqrt and 10 div a force evaluation; 17,544 waves
hold 1,000,008 particles; Pan caps the step; only a first-attempt
rejection leaves state; a record names no non-default `--adaptive-mode`,
`--softening`, `--min-dt` or `--no-flag-abort`; the exact-time pair comes
from `vaugadd`; R = 60,330 km is still believed.

Believed here, by message from cft-fp256's lead: finding 1 is confirmed,
and GPU-E's item 3(b) and its question 4, which leans on it, are to be
corrected in a new commit on cft-fp256's main, credited as measured
here. Its J2 and J4 row, an estimate made before code existed, is to
become P1's operations as read at `291c99c`. Its cost model takes
ceil(N/57), a last wave of 49 particles, which this tree cannot write.
Measured on 2026-11-12: at cft-fp256's `be40005b`, the study's file is
unchanged since `56006987`.

### Stated limits, not closed

- **A producer's own attempts records.** A particle is sampled only in a
  member this tree reran, from the rerun's own records. Taking a GPU
  run's attempts records and schedules, so particles are sampled without
  a rerun, is not built (docs/RING.md).
- **p2 does not plant steps_done or steps_rejected,** which stay with the
  record's structure checks and the member path, and it plants each
  column on one line (VL).
- **Part 3's control finds its token with `compare_run`'s own column
  reader.** `check_compare_run`'s independent reader covers a bug shared
  there (VL).
- **The GPU seam test,** a difference planted in the first GPU record and
  shown caught, belongs to the GPU half's first brief (ROUND4.md).

### What this does not establish

- **Anything on a GPU, a card or WSL.** None was used in round 4.
- **A full suite on Linux** at any of these merges.

## 52. After round 4: `tools/state_to_problem.py` carries the `active` and `zonal` lines, and a reversal is gated; forward with J2 and J4, it had come back without them

**2026-11-12.** This desktop, on the software backend, shared with round
4's last suite; pins as in entry 46. The lead's commits on
`state-to-problem-lines`, from `de376d8`: `68551e2`, the tools, the gate
and its leg; `ac2a152` and `06c4384`, after VS's first and second looks;
and this one, the records. Logan was asked, after study DEBRIS-A found
the gap, whether to fix it with a test once round 4 had landed: "Yes,
after the round". VS, the fix's verifier, said READY on `68551e2`, with a
defect and a gap in the gate that `ac2a152` answers, and READY on
`ac2a152`, with three items that `06c4384` answers. Its look at `06c4384`
and these records, and this commit's own suite, are in the report to
Logan.

### What was wrong

`tools/reversal.py` runs N fixed steps forward, writes the last sample as
a problem through `tools/state_to_problem.py`, runs N steps back from it
and compares the return with the start (the reversal test of
2026-09-10). state_to_problem.py wrote the name, G, N and the bodies.
Round 4's P0 gave the grammar its `active` and `zonal` lines, and no one
owned this tool, so it dropped them (entry 46's limits, READ by
DEBRIS-A). No gate ran either tool.

Measured here, with `de376d8`'s tools, on `zonal_saturn` with `active 1`
for 100 steps of 100 s each way: the return misses the start by 3.952e-3
in max |dx| / L at binary64, binary128 and binary256 alike. That is the
J2 and J4 the backward leg never felt. With `68551e2`'s tools it misses
by 9.265e-17, 8.036e-35 and 1.549e-43.

### The fix, `68551e2` to `06c4384`

- **state_to_problem.py** reads the record's `# active` and `# zonal`
  lines, which the program writes only for a problem that has them, and
  writes them back as the problem's `active` and `zonal` lines. Each value
  is the record's own exact hex, never a float. It refuses an unreadable
  line, a second one, and a count that does not fit the system's N. A
  record without them gives the file the tool always gave.
- **reversal.py** refuses a turn or return problem whose lines, other
  than comments, name, N and bodies, are not the original's as the run
  reads them, value for value. That is how `read_problem()` reads them:
  - hex is rounded to the run's format, ties to even, with its sign, as
    the program rounds it. VS measured `rounded()` against the program in
    87 comparisons, 29 values at each of the three formats, with
    subnormals, ties each way, carries into the next binade, negatives
    and zeros among them. All agree but the two values that overflow
    binary64, at binary64, where the program reads infinity and stops on
    its flags before any comparison;
  - integers compare by value;
  - G, E and a system's name compare by their first field only.

  `68551e2` compared every token exactly and as text. VS found that
  refused a value wider than the run's format, which the program rounds
  and the turn problem carries rounded: a regression for G, which
  `de376d8`'s tool had passed. It refused `active 01`, `zonal 00` and a G
  line with a trailing comment too. `ac2a152` compares as above.

  A line the tool does not carry is now an error by name, not a reversal
  of another problem: an ensemble's `E` and `system` lines are refused,
  as they failed before. reversal.py also refuses a return whose bodies
  are not the start's: it had zipped them, and a short return would have
  compared short.

  Its legs' command became `leg_command()` and its comparison
  `deviation()`; the command is the same list in the same order. Its
  energy, the point-mass energy of every body, is now labelled a
  function of the state, not the energy an `active` or `zonal` run
  conserves.
- **`tools/check_reversal.py`,** a leg in every tier, the 25th: 39
  checks in 71 s in `make check`, and 21 in 31 s in `--quick`, on a host
  shared with a full suite.

### The gate

The problem is `zonal_saturn` with three changes, each there so that a
wrong tool fails:
- **`active 01` for its `active 2`.** Mimas, GM 2.5, becomes a test
  particle with a mass. As committed, the file's one test particle is the
  massless ring particle, and a massless body's reaction is an exact zero
  when every body is active. Measured: the return's three body lines are
  identical with and without the `active` line, so the file as committed
  cannot hold that line. check_ring_seams 1(b) found the same of the
  ring's particles. The leading zero, which the program reads as 1,
  holds a comparison of integers as text.
- **G = 1 + 3 * 2^-53, with a comment after it.** At binary64 it is a
  tie, which the program rounds to even, up, to 1 + 2^-51; the comment
  holds a comparison that reads more of G's line than its first field.
- **R = 60,330 km + 2^-38.** At binary64 it is a tie, which rounds to
  even, down, to 60,330 km.

G and R have 54 significant bits, exact at binary128 and binary256. A
value that binary64 holds cannot show a carry that narrows it, and a tie
each way shows a comparison that truncates or rounds ties away. The
gate came to this in two of VS's looks:
- **On `68551e2`** its values were all binary64. It passed a refusal
  that compared first words only (VS's M2) and a carry through a float
  (M3).
- **On `ac2a152`** G and R were 1 + 2^-60 and 60,330 km + 2^-80, wide
  but not ties. It passed a comparison that truncated (M6), rounded
  ties away (M12), read all of G's line (M8) or compared integers as
  text (M9). Its check that the refusal named the line could not fail,
  because the problem's file name holds both words (M15). And it held
  positions only.

The gate's four parts:
- **0. The carry alone.**
  - (a) The 12 committed records of results/horizon/reversal/ (kepler
    and outer, three formats, forward and back) rebuild their 12
    committed problems, every line but the one naming the record's path.
  - (b) Two steps at binary128 and at binary256: the turn problem
    carries G, `active` and the zonal line exactly, the wide G and R
    included. This runs in `--quick` too.
  - (c) The hex reader reversal.py and the gate share, check_equivalence's
    `hex_to_fraction`, reads known answers right: -0 and +0 apart,
    -0x1.8p+1 as -3, and the widened R. state_to_problem.py parses no
    hex: it copies the record's text.
- **1. The backward leg runs the problem.** At each format (binary64
  only in `--quick`), reversal.py runs (rc 0).
  - Its turn and return problems carry G, `active` and the zonal line,
    value for value as the run's format reads the problem's. The gate
    computes this itself, not through the tools: hex at binary64 by
    CPython's correctly rounded division and wider exactly, an integer
    by its value, and G by its first field.
  - The backward record's `# active` and `# zonal` lines are the forward
    record's, word for word, and name the file's counts: N - 1 = 2
    bodies feel the term, and body 0 takes the reaction of
    active - 1 = 0.
  - The return lies within the bound in positions and in velocities.
    The gate computes both exactly, and reversal.py prints the same:

| format | max abs dx / L | max abs dv / V | bound |
|---|---|---|---|
| binary64 | 9.265e-17 | 1.212e-16 | 1e-14 |
| binary128 | 8.036e-35 | 1.052e-34 | 1e-32 |
| binary256 | 1.549e-43 | 1.339e-43 | 1e-40 |

  At binary64 and binary128 the return is round-off. Its largest
  position error is one unit in the last place of a coordinate in
  [2^16, 2^17), with L = 157,068 km: 2^-36 / L and 2^-96 / L (VS), 0.83 of
  each format's unit roundoff (COMPUTED). At binary256 it is 3.4e28 unit
  roundoffs: IAS15's own reversal error. The method is not
  time-symmetric (reversal.py's head), and the narrower formats'
  round-off hides it. The bounds are about 100 unit roundoffs at
  binary64 and binary128, and 646 times the binary256 return.

  `ac2a152`'s message calls its binary128 return, of 1 + 2^-60 and
  60,330 km + 2^-80, bit-identical. It was not: the positions came back
  exactly, but three velocity components missed (max |dv| / V =
  5.258e-35, VS). The gate now holds velocities too.
- **2. Each line carries weight.** At binary64 the turn problem, with its
  zonal line deleted and then with its active line deleted, runs back as
  reversal.py runs its backward leg. Without the zonal line the return
  misses by 3.952e-3, 4.0e11 times the bound; without the active line,
  by 1.666e-7, 1.7e7 times.
- **3. A dropped or changed line fails.** These are controls from copies
  of the three files, each mutation shown to apply by its diff. Three
  copies of state_to_problem.py:
  - one drops the zonal line;
  - one drops the active line;
  - one writes J4 one unit up in its last hex digit.

  reversal.py refuses each copy's turn problem before the backward leg
  (rc 1), in its own words: "it lacks its zonal line", "it lacks its
  active line", "its zonal line(s) differ". With reversal.py's two
  refusals deleted as well, each copy runs to its end (rc 0), and part 1
  fails exactly the checks the copy should break:
  - a dropped line: that line in both problems, the records' lines, and
    the return, 3.952e-3 (zonal) or 1.666e-7 (active) in positions;
  - the changed J4: the zonal line in both problems and the records'
    lines. The return stays at 9.265e-17, so at binary64 only the value
    checks see a one-unit J4. At binary128 the return would: 6.781e-21
    (VS).

  Each of VS's mutants, M2, M3, M6, M7, M8, M9, M12 and M15, now fails
  the gate in `--quick`, measured by the lead, each in a scratch tree of
  its own.

### The suites

Round 4's records ran their own suites, which entry 51 left to the
round's closing report. They are recorded here, in the first commit
after them; this commit's own suite is in its report to Logan. Read from
their logs:

| commit | tree | legs | result |
|---|---|---|---|
| `37a1ce0` | round 4's records and the ledger's archive | 24 | rc 0, 1424 `ok`, 0 FAIL, 2 h 30 min (11:51:49 to 14:21:52); every leg at `d9c8c8c`'s ok count, and the log `d9c8c8c`'s line for line but for times, temp-directory names and the commit |
| `de376d8` | VD's second-round fixes to the records | 24 | rc 0, 1424 `ok`, 0 FAIL, 3 h 04 min (14:39:42 to 17:44:10), on a host shared with the lead's runs and VS's; every leg at `37a1ce0`'s ok count, and the log `37a1ce0`'s line for line but for times, temp-directory names and the commit; main fast-forwarded `d9c8c8c..de376d8` and pushed, closing round 4 |

### Two texts of round 4's records, corrected here

- **ROUND4.md's status line** put seam part 1 with P3. It landed with P2
  and P4 at `7d50392` (VD's third look, not blocking there).
- **The usage block's `--zonal-selftest` line** said "its own synthetic
  constants". R = 1.5 is synthetic, and J2 and J4 are Saturn's (VD).

### Stated limits, not closed

- **One problem holds it,** with one oblate body, one test particle with
  a mass and one without, at one dt. The binary256 bound is IAS15's
  reversal error at this dt; another dt or problem needs its own.
- **A problem line added later** is dropped by state_to_problem.py until
  it is carried there. reversal.py refuses that by name, but
  state_to_problem.py's other use, `--system` on an ensemble record
  (docs/ENSEMBLE.md), has no such check.
- **reversal.py's reading of each line is a copy of read_problem()'s:**
  which fields G, E and `system` read, and that hex is rounded to
  nearest, ties to even. A change to the grammar must change it too.
- **reversal.py runs one problem.** An ensemble is not reversed.

### What this does not establish

- **That IAS15's reversal error is 1.5e-43 anywhere else.** It is this
  problem's, at this dt, at binary256.
- **Anything on a GPU, a card, WSL or Linux.**

## 53. The pin moves to cft-fp256 `eae74e3`, its main with G1B's HIP backend and G1C's gates: the same bits from the same calls, through a library whose executor moved

**2026-11-13.** This desktop, on the software backend only, booked
shared in the TimeSheet (B0033 and B0041); MSYS2 mingw64 gcc, Miniconda
Python 3.12.15 with mpmath 1.3.0; REBOUND at bdfda4bd. The lead's
commits on `pin-g1c`, from main at `2d52e32`, are:
- `47b2398`, the pin;
- `e520bc8`, the records;
- `58b085c`, after the verifier VP sent them back, with check-pins
  mended;
- `4dd8bd9`, after VP's second look: its two sentences, `58b085c`'s
  suite and check-pins' two edges;
- this one, the line above.

It is Logan's decision 7, relayed by cft-fp256's lead on 2026-11-13:
"Yes, to G1C", the SHA to follow once G1C's merge was pushed. As entry
38 did, this entry moves the pin and establishes one thing: the calls
this port makes return the same bits through the new library. It adopts
nothing, and builds no HIP.

### Which commit

`eae74e3eddc67836d4b7c8922badd40e6e1bdb2f` is cft-fp256's main as pushed
on 2026-11-13. READ, from its lead's message, it holds:
- G1C's merge, `533bbb58`: the stage `hip` in cft-fp256's gate budget,
  and `hip-nvcc` and `hip-source` in its full census, held to the
  software backend's bits on the RTX 5060 Ti and the RX 7600;
- G1B's merge, `9b358e66`: the HIP backend, host/src/backend_hip.cpp,
  and the `hip:N` device;
- docs, a study's merge (GPU-F) and night 20's page data.

Its quick gate passed on amd-arc-box at the G1C merge's tree: 24 ok and
7 skipped by name, the three GPU stages among them, since that box has
no GPU. Its CI on `eae74e3e` had started when the SHA was sent. MEASURED
(`git rev-parse`), `533bbb58`'s host/src, host/include and host/Makefile
trees are `eae74e3e`'s, so that gate covered the pinned library's
sources.

What this repository builds of it is host/src, host/include and
host/Makefile. MEASURED with `git diff --stat 922c2ade eae74e3e` and
`git rev-list`: 13 of those files changed, +3,752 and -907 lines, in 12
of the range's 239 commits (10 of them not merges; 28 touched host/).
- The program executor moved out of program.c (+20, -793) into a new
  seq_exec.h (925 lines), which the device shares. VP measured it the
  same code but for the helpers' parameters and the ALU call behind
  `CFT_SEQ_ALU`, whose default is `cft_sf_compute`. seq_exec.h's own
  header says the moves are not pure.
- softfloat and bigint gained only `CFT_HD` annotations, and their
  objects are identical (VP). transcend.c changed only `do_exp2m1`, and
  this port calls no transcendental function. device.c gained HIP-only
  blocks and one refusal's sentence (VP).
- cft.h changed by +5 and -2 lines, and cft_config.h gained 30.
- backend_hip.cpp (2,321 lines) is new, and not built here: libcft's
  Makefile builds it only under `HIP=1`, which this repository's sub-make
  does not pass, so libcft.a is its 16 C objects (VP).

The ABI is 0.18 at both commits (cft.h's `CFT_ABI_VERSION_MINOR`).

### How it was run

One worktree from main at `2d52e32`, with fresh clones of its pins:
REBOUND from main's own clone at `bdfda4bd`, and cft-fp256 from the
sibling checkout `../cft-fp256` at `eae74e3e`. Both are local; nothing
was fetched from GitHub, by this repository's rule since round 4, so
unlike entry 38 the MANIFEST line's URL was not exercised. `make
check-pins` passes.
- **The build.** `make all programs` took 63 s. It carries the five
  warnings every build has carried since entry 38, all REBOUND's own
  (simulationarchive.c), and none from libcft or this repository.
- **No source edit,** as at entry 38.
- **The baseline** is `2d52e32`'s own suite (entry 52): the same tree at
  the old pin, run the evening before.

### The suite at both pins

| | old pin `922c2ad`, at `2d52e32` | new pin `eae74e3`, at `47b2398` |
|---|---|---|
| `make check` | rc 0, 3 h 02 min | rc 0, 2 h 44 min |
| legs | 25 | 25 |
| lines beginning `  ok ` | 1463 | 1463 |
| FAIL verdicts | 0 | 0 |
| named skips | gate_subprocess's body_count row, deliberately untested | the same |

Read from the logs, not the exit codes. Every leg has the same ok count,
digest count and line count in both. Normalising dates, temporary
directories, the tree's path and check_ring's elapsed seconds leaves 56
differing lines:
- 19 in each log carry the library's own build identity,
  `libcft_build=commit=...`: 18 of them gate_real's and 1 check_door's;
- 2 in each are check-pins' pin lines;
- 2 in each are the TimeSheet's booking lines, and 1 in each is the
  commit line;
- the new log alone has the build lines of six gate programs (7 lines).
  The new tree built them in its suite, against the new library; the
  old tree's had been built before its own. The old log alone has one
  TimeSheet notice.

No verdict, case or count differs. The archive gates' `abi` is "0.18" in
both, where entry 38 saw it move. Entry 45's six wide-state digests are
in both logs at the same counts, so the whole wide state is held to the
bit at every format across the two libraries: fp64 `e08541bc0af9b3ec`
and `56bacf3bf63f19a7`, fp128 `5d7cf4588fad050e` and `cd16e6dc29aa0589`,
fp256 `a1052750e206bebe` and `ba4d447c8e58bd20`.

**The branch tip.** `58b085c`, which carries check_pins.sh's fix, ran its
own full make check: rc 0, 25 legs, 1463 lines of `  ok `, 0 FAIL, from
12:28:31 to 15:18:45. Every leg has `47b2398`'s counts, and the six
digests are present. The commit after it changes two sentences of these
records, so that suite holds it docs aside.

### The same integrations through both libraries

**The lead's,** a script in the session's scratchpad. The program,
built at each pin, ran every problem in data/problems, each at
binary64, binary128 and binary256: on the loop engine in REBOUND's order
and in the FMA form, and on the program engine with that build's own
assembled programs. That is 54 runs through each library. The arguments
were:
- kepler, kepler_e09 and pythagorean: `--dt 0.01 --epsilon 1e-9
  --steps 100 --sample 25`;
- kepler_e099: `--dt 1e-4 --epsilon 1e-2 --steps 40 --sample 20`;
- outer: `--dt 40 --epsilon 0 --steps 60 --sample 20`;
- zonal_saturn: `--dt 100 --epsilon 1e-9 --steps 60 --sample 20`.

All 54 records are the same, line for line, every counter included,
with the identity line set apart and the trailer's timings removed. It
had no control. The identity lines, read afterwards, name
`922c2ade...` and `eae74e3e...`, and the two libcft.a files differ
(sha256 `1c18ac85...` and `e07daf2f...`).

**VP's own,** with other arguments and with controls:
- **95 of 95 runs identical,** every record line byte for byte, with
  stderr and output files, and each identity line checked to name its
  own build.
- **What they covered:**
  - the 54 again, at other dt, epsilon and steps;
  - `--cs augmented`; adaptive modes 0, 1 and 3; softening with a
    minimum dt; `--max-iter 30`; `--pc-tol-shift`;
  - all four routes forced on;
  - 16 bodies, and ensembles at three formats on both engines;
  - the selftests and `--dump-constants`;
  - the `--dt-out`, `--attempts-out` and `--record-iters` files.
- **The controls,** through the same comparator, each differed: a
  corrector cap of 12 against 60 at binary256, dt one ulp apart, one
  planted digit, and two records in swapped slots.
- **A checkpoint across the pins,** as entry 38 did. One build saved and
  the other resumed, both ways, at three formats and in the plain,
  regrow and adaptive sequences: 18 of 18 runs are identical to both
  straight runs, and the NO_ADOPT control differs in all 18.
- **The compiled outputs.** All 66 program images and both embedded
  image headers are byte-identical across the pins, and so are the
  port's own objects. Of libcft's 16 objects, only device.o, transcend.o,
  program.o and build_id.o differ.

**`make example`,** as entry 38 ran it, at both pins: rc 0 at each, and
the same 38 lines of output with the tree's path normalised, but for the
column of seconds. The drop-in example's reload is restored exactly at
both.

### Also in these records

- **HORIZON.md's two tile lines are corrected.** "Nothing on the tile,
  whose rate is format-blind" and "a few thousand system-steps a second
  per tile" stated a belief and a projection as rates. COMPUTED from
  entry 37's measured card rates and pass counts, binary128 costs about
  1.1x binary64 a corrector pass: 4,529 against 4,917 system-passes a
  second on one tile at E = 4,096. COMPUTED from entry 45's run times,
  binary256 runs 88 system-steps a second for 64 Kepler systems on one
  tile and 204 for 512. HARDWARE.md's 1,300 is a projection. They were
  found while answering cft-debris's question for its study DEBRIS-B.
- **Logan's answers to decisions 3, 4, 6 and 7, and to GPU-E's question
  4,** as cft-fp256's lead relayed them, are in ROADMAP.md and
  docs/RING.md.
- **check-pins' README test is mended.** It was a substring test, and
  README's row names earlier pins as its history, so it passed a
  MANIFEST and a clone reverted to `922c2ad` beside today's README (VP).
  tools/check_pins.sh now requires one cft-fp256 row whose commit cell
  begins with the pin.

  Measured on four scratch roots, with the old script and the new:
  - the real tree passes both;
  - the reverted one passes the old and fails the new;
  - a README with two cft-fp256 rows passes the old and fails the new;
  - a row beginning with the old pin fails both.

  Its ok line reads as it did, so suite logs stay comparable.
- HARDWARE.md:502 and ROADMAP's paragraph of 2026-10-30 now give
  `922c2ad`'s time as the pin an end.

### What this does not establish

- **Anything on a GPU.** HIP is not built here; `hip:N` is cft-fp256's,
  gated there (its G1C).
- **The card.** No XRT build was made at this pin.
- **cft-fp256's CI on `eae74e3e`,** which had started when the SHA was
  sent (READ).
- **A full suite on Linux** at this pin.
- **check-pins against every bad README.** It still passes a correct row
  beside a second, stale row written with two spaces after its pipe,
  which it does not count. It refuses a correct row that writes the pin
  in 8 characters or in full (VP).

## 54. Round 5's P0: this tree built with libcft's GPU backend, and the software backend's bits on an RX 7600 through `hip:0`

**2026-11-16.** The desktop for the build rules, and nvidia-box for the
GPU: the RX 7600 (gfx1102), Linux 7.0, ROCm with HIP 7.16, gcc 13.3,
Python 3.12.3. Pins as in entry 53. The GPU window was booked in the
TimeSheet (B0043 and B0044); P0's build and smoke, short runs, went
first, after `check --short` found the window free. The lead's commits
on `round5-p0`, from main at `1354a05`: `d95ac0d`, the Makefile;
`5389ff7`, ROUND5.md; and this one.

### The build

The Makefile's `HIP=1` (`d95ac0d`) passes `HIP_ROOT`, `HIP_PLATFORM`,
`HIP_CUDA` and `HIP_ARCH` through to libcft's own Makefile. It links the
runtime libcft's HIPLIBS names: libamdhip64 and libstdc++ on ROCm; the
CUDA runtime with libstdc++, pthread, dl and rt under HIP's NVIDIA layer.
`HIP=1` with `XRT=1` is refused, as libcft refuses it. `HIP=0`, the
default, adds nothing.

On nvidia-box, ~/cft-rebound-r5 was cloned from git bundles of this
repository at `round5-p0` and of its two pins, made from the desktop's
local clones and carried over SSH. Nothing was fetched from the web.
`check_pins.sh` passes there. `make HIP=1 PYTHON=python3 all programs`
built in 31 s, and `build/ias15_cft` links ROCm's libamdhip64.so.7.

It has 15 warnings:
- 14 are REBOUND's own: `fread` and `vasprintf` results unchecked, in
  binarydata.c and simulationarchive.c, under Ubuntu's fortified glibc.
- 1 is this port's. GCC 13 warns at src/ias15_cft.c:2155 that
  `vcopy`'s source is null. V0 found that it comes only from the
  library object, build/ias15_cft_lib.o (`-DIAS15_CFT_LIBRARY`). There,
  `load_programs()`, which holds `prog_dep`'s only allocation (:1969), is
  compiled out, so `prog_dep` is null outright. It is harmless only
  because the library sets `engine_program = 0` (:4955) and nothing sets
  it to 1, so the call at :2953 is never reached. A library path that
  turned the program engine on would copy from null. The program's own
  object does not warn, and MinGW's GCC 16.1 raises it in neither. The
  lead's first account, that the deposit is allocated before any program
  runs, holds for the program, not for the object that warns.

### The smoke: one binary, the software backend against `hip:0`

kepler, `--dt 0.01 --epsilon 1e-9 --steps 50 --sample 10`, at binary64,
binary128 and binary256, on the loop engine and on the program engine in
the FMA form. Each case ran twice from the same build, once on the
software backend and once with `CFT_REBOUND_ARTIFACT=hip:0`
`CFT_HIP_INFO=1`. MEASURED:
- **The bits.** Every sample line and every counter is the same in all
  6 cases. The physics counters, the same both ways: 50 steps; 3.44
  passes a step at binary64, 8.52 and 8.58 at binary128 (loop and
  program), and 12.00, the corrector's cap, at binary256.
- **The two data lines that differ are routes.** The comparison sets
  aside the header's `backend=` field and the trailer's timings. Of the
  rest, `indexed` and `pair_gather` read `by=absent` on `hip:0` and
  `by=opt_in` on the software backend. The device publishes no INDEXED,
  and both routes are off either way.
- **The device line:** "hip:0 = AMD Radeon RX 7600 (gfx1102, 16
  multiprocessors, bus 0000:03:00.0, 8176 MiB), HIP platform amd,
  runtime 71626385, driver 71626385; of 1 visible; this libcft: HIP_ARCH
  gfx1102, stage 64 MiB, page 256 bytes, kernels 46/46/64 registers and
  5572/5572/13716 bytes local a thread (run/level/seq)".
- **Wall seconds,** the program's own clock, one run each, not ranked:

| format | engine | software | `hip:0` |
|---|---|---|---|
| binary64 | loop | 1.10 | 31.6 |
| binary64 | program | 0.85 | 7.98 |
| binary128 | loop | 3.29 | 85.9 |
| binary128 | program | 2.40 | 18.8 |
| binary256 | loop | 7.13 | 154.0 |
| binary256 | program | 3.78 | 35.1 |

  The default path is bound by its calls on a GPU, as on the card (entry
  37), from 2 to 64 bodies as V0 measured it; nothing above 64 bodies is
  measured. cft-fp256's G1B measured 234 µs a call on this GPU, staged, in its
  first part (READ). V0 measured this port's own, `wall_seconds` over the
  trailer's calls: 286 µs on kepler's loop engine at binary64, and 218
  to 769 µs over its 24 runs (binary64 218 to 364, binary128 336 to 472,
  binary256 364 to 769).

The comparison first called all 6 cases different. It had not expected
the route lines; the data never differed.

### V0's look

V0 verified P0 independently: its own bundles, its own build in
~/cft-rebound-r5-v0 (29.4 s), and its own comparator, booked in the
TimeSheet (B0050). MEASURED:
- **`HIP=0` changes nothing.** `make -n all programs` is the same at
  `1354a05` and `d7b5047`, on the desktop and on nvidia-box. `HIP=1`'s
  link lines match libcft's HIPLIBS for `amd` and for `nvidia`, and
  `HIP=1 XRT=1` is refused.
- **24 of 24 cases identical** in every sample line and counter:
  - kepler, outer, pythagorean, and zonal_saturn with its `active` and
    `zonal` lines;
  - a 4-member ensemble, and a 64-body ring wave on both engines;
  - `--cs augmented`;
  - kepler_e099 rejecting its first two attempts before any accepted
    step, with identical `--attempts-out` and `--record-iters` files;
  - a 3-member rejecting ensemble, and the six cases above.

  Every `hip:0` run carries exactly one device line, and no software run
  has one. Physics beside the bits: 6 to 50 steps, 0 to 6 rejections,
  3.15 to 12.00 passes a step, |dE/E| at most 1.4e-14. Zonal on against
  off, and augmented against Kahan, differ where they should, so neither
  comparison is vacuous.
- **Controls:** eight kinds of planted change, on copies of eight cases'
  GPU files, were all caught: three hex digits, two counters, the
  header's backend, the device line and an attempts file.
- **The routes,** forced on, at binary64 and binary128: on `hip:0`,
  `scalar` and `reduce_seg` run in the software backend's bits, and
  `indexed` and `pair_gather` are refused by name. `--door-selftest` on
  `hip:0`, at binary64, holds 116 cases identical, with 48 not compared
  because their routes are absent.
- **The ratio to software depends on the problem.** A 64-body ring wave
  ran faster on `hip:0` than in software: 6.31 against 7.19 s on the
  loop engine, and 4.85 against 7.26 s on the program engine, one run
  each. Its GPU time is still its calls, at 244 and 271 µs a call. What
  changes at 64 bodies is that software is slower.
- **The record's `backend=hip:0` echoes the artifact string**
  (src/ias15_cft.c:4697-4699). It does not prove the GPU ran; the device
  line does. The trailer's `seconds` and `steps_per_s` are CPU time, 1.5
  to 2.0 times the wall time over V0's 24 `hip:0` runs. That ratio is
  ROCm's. On the RTX 5060 Ti, through HIP's NVIDIA layer, it is 1.02 to
  1.03 in all 12 of P1's cost rows, where P1 measured 1.54 to 1.92 on this
  GPU (entry 55).
- **Found:**
  - `HIP_PLATFORM=cuda` falls into the amd branch. `make -n` shows
    `-lamdhip64` on its link, so it is expected to fail there; not run,
    since nvidia-box has no CUDA.
  - `make HIP=1 example` fails (run): examples/Makefile has no HIP block.
  - In a tree first built without HIP, `make HIP=1 all` leaves libcft.a
    software-only, and `hip:0` is then refused by name: `make clean`
    first.

  The first two are not this round's path. The third is, and ROUND5.md
  lists it among the rules. Each that ran failed loudly.

V0 sent P0 back for the warning's account, the per-call figures' source
(G1B's, not G1C's), "refused by its loader" for index tables (they are
refused at the call), P1's proof that the GPU ran, and RING.md's stale
belief about AUGADD. This entry and ROUND5.md now carry its
corrections.

### What this does not establish

- **The RTX 5060 Ti.** P0 did not build in WSL `cft-gpu`; that is P1's.
- **The suite on `hip:0`.** One problem was run, not the gates: P1.
- **The Python gates on nvidia-box.** Its Python has no mpmath, which
  they need (ROUND5.md, decision 2).
- **This commit's own suite,** which is in the report to Logan.

## 55. Round 5's P1: this tree on `hip:0` as it is, every process a leg starts held to the software backend's record; `make HIP=1 check-gpu`, whose GPU seam went back three times for a copy of the driver comparing less

**2026-11-16 to 2026-11-18.** Two GPUs, each window booked in the
TimeSheet: the RX 7600 on nvidia-box (ROCm with HIP 7.16, gfx1102; strace
and gdb installed), and the RTX 5060 Ti in WSL `cft-gpu` (HIP 7.2.4's
NVIDIA layer, sm_120; no strace and no gdb), where
`CUDA_DEVICE_ORDER=PCI_BUS_ID CUDA_VISIBLE_DEVICES=1` keeps `hip:0` off
the RTX 4070. The desktop, shared, for P1's suite and its dry-run
harness. Pins as in entry 53. mpmath 1.3.0 on both GPU hosts, copied
from the desktop: cft-gpu's own Python has 1.2.1 (P1). P1 branched from
`5b1cac6`. Its commits are 16 to `ab75fbe`; `a068acc` and `6a5e448`,
send-back 1; `fc4ba21`, send-back 2; and `0f324bd`, send-back 3.
`fda3b40` merges `0f324bd` on `round5-int` at `4a8fe99`, and main
fast-forwarded `4a8fe99..fda3b40` and was pushed on 2026-11-18.

### What P1 is, in one paragraph

`make HIP=1 check-gpu` runs `tools/check_gpu.py` on a GPU host. It is in
none of the three tiers and does not use the gate cache: a GPU run is its
own evidence, and a cached pass would carry no device line. On a host
with no `hip:0`, the desktop among them, it refuses by name, rc 2 (V1).
Its variables are `GPU_ARTIFACT` (default `hip:0`), `GPU_EXPECT` (text
the device line must contain, or nothing runs), `GPU_LEGS` (a named
subset) and `GPU_SKIP` (the legs left out). The driver has five parts:
- **0, the device.** `ias15_cft --probe` opens the artifact with
  `CFT_HIP_INFO=1`. A host without it, or a libcft.a without
  backend_hip.o, is refused by name; under the Makefile's `HIP=1` the
  second is named as V0's stale-library trap (entry 54). Every later
  process must print the probe's device line, and the routes the device
  does not publish are read from the probe, not typed.
- **1, the witness.** Each leg runs against a shadow of the build
  directory in which `ias15_cft` and `check_dropin` are a witness: a
  script that runs the real program with the caller's arguments,
  environment and stdin, keeps what it printed, and, when the run named a
  device, runs the same arguments again on the software backend. So for
  every process a leg starts, whatever starts it, the driver holds:
  - its device line: at least one, and every one the probe's, or the
    process is a stand-in and FAILs its leg by name;
  - its stdout and its `--dt-out`, `--attempts-out` and `--record-iters`
    files, byte for byte against the software backend's. Set aside: the
    header's `backend=` value; the door and corrector self-tests' first
    lines, which name the backend they ran on; the trailer's 13 timings,
    each by its name; and the `by=` of a route the device does not
    publish, counted;
  - its physics beside the bits: steps, rejections, passes and flags.

  A process refused by its own arguments before it opens a device, as
  the software backend refuses the same arguments, is counted apart.
  stderr is compared, and a difference printed. The leg's own verdict is
  required too. Where strace is installed, nvidia-box, every exec of a
  witnessed program must be one its witness accounts for.
- **2, the routes,** at every format.
- **3, the GPU seam** (ROUND4.md, "What the integrator keeps"), below.
- **4, the legs,** eleven, ten in their quick form and check_dropin in
  its `--light` form, their commands read from the Makefile's recipes
  and never retyped: check_records, check_equivalence,
  check_program_engine, check_ensemble, check_ring, check_replay,
  check_compare_run, check_zonal, check_reversal and check_corrector as
  check-quick runs them, and check_dropin_light as check-light runs it.

`--reach` runs legs it does not hold through the witness and reports
which processes opened the device: a survey, never a verdict.
`--analyse` re-reads a kept leg and opens no device. `hw/gpu/` holds the
measurements: `where.py`, what runs where; `costs.py`, what a call
costs; and `program_caps.c`, what each program asks of the device's
loader.

### The gate, and the controls that show it can fail

**The witness's controls,** permanent, each a real run through the
witness that must be caught by name. As the GPU seam run on the merged
tree printed them, live on the RTX 5060 Ti (below):

| control | planted | result |
|---|---|---|
| a | one hex digit of body 1's x in the device's record | caught: stdout line 10, x of body 1 |
| a2 | the trailer's `calls=` one more, on the line whose timings are set aside | caught |
| a3 | `t_new=1`, a `t_` key that is not one of today's timings, on that line | caught at line 11, field 13, beside the clean run's 12 timings there, each named and set aside |
| b | an attempt's outcome turned from a to r in the device's attempts record | caught |
| c | the device line's last character changed | caught as not the probe's |
| d | the run with `CFT_REBOUND_ARTIFACT` empty | caught as a stand-in |
| e | the run without `CFT_HIP_INFO` | caught for its missing device line |
| f | check_records' quick form with no artifact, which passes by itself | its 4 processes caught as stand-ins, so the leg FAILs |
| g | a process that runs `ias15_cft` past its witness | caught by strace, on the RX 7600; not run in cft-gpu, which has no strace |

**V1's own plants:** 32, in copies of real `hip:0` witness data, read by
the driver's own `analyse()`. Each plant beside a set-aside rule was
caught by line and column: the tokens either side of `backend=`; the
identity line; eight counters beside the timings; `calls=`,
`fallbacks=` and on or off beside `by=absent`; the corrector banner's
tail and its format; the zonal line, a member's trailer,
`pc_lane_efficiency` and `t_lo`. Not caught, by design: the `backend=`
value and `wall_seconds`.
Six real processes through the witness on the RX 7600: the argument
refusal of check_corrector's `--trace-pc` beside `--replay-iters` was
counted as refused, and `hip:1`, a hidden device, `CFT_HIP_INFO=2`, a
SIGSEGV after the device line and a SIGKILL at once each FAILed. A
failure the device caused cannot be refused alike in software, which
never opens a device.

**The GPU seam, as it landed at `0f324bd`.** tools/make_ring.py's
two-wave ring with its zonal line, recorded on `hip:0` over 8 steps from
a first step of 1e5 s, so that each wave rejects once and both are
sampled. tools/compare_run.py, run with no artifact, holds the record
against this tree's CPU reruns. A particle's replay carries four classes
of value: the schedule (t, dt_next, dt_last, E and the exact-time pair),
the trailer's counters, the wave's massive bodies, and the particle's
own body. Each plant is one digit on one line, one plant a run, and its
class is derived from compare_run.py's own particle lines: the sampled
bodies, 22 and 30 in wave 0 and 57 and 31 in wave 1, and 7 massive
bodies.
- **The schedule:** t, dt_next and the exact-time pair's t_lo.
- **The trailer's counters:** mean_pc_iterations, on the wave's trailer
  line.
- **A massive body:** x of body 1.
- **A sampled particle's own body:** vy of the wave's first sampled
  particle.
- **A body no sampled replay carries:** vz of body 63.

Each plant must FAIL its wave, naming the line and column, and FAIL
those of the wave's sampled particles whose replay carries it: both for
the first three classes, one for the fourth, none for the fifth. The
other wave and its particles must stay IDENTICAL. The sample-line plants
sit on each wave's first and last sampled lines. REBOUND's form makes 14
plants, 7 a wave; the FMA form, the form P2's kernel is held in, makes
8, every class but the fifth in each wave. A last verdict in each form
holds that every class has a plant in each wave.

Mutants of compare_run.py, each one line, landing once or refused: P1's
two member mutants and V1's four. On P1's desktop harness at `0f324bd`
(a software record, the witness stubbed), each gives the same FAIL set,
plant for plant, with this tree's compare_run.py and with P3's:

| compare_run.py that | ok | FAIL |
|---|---|---|
| is the tree's | 28 | 0 |
| never compares member 1, or member 0 | 17 | 11 |
| never compares a particle | 8 | 20 |
| skips a particle's own body | 24 | 4 |
| never compares a member's counters | 24 | 4 |
| never compares a particle's counters | 24 | 4 |

### The routes, what runs where, and what a call costs

**The routes,** at binary64, binary128 and binary256, on both GPUs.
Forced on, `scalar` and `reduce_seg` give the software backend's bits:
every line is the unforced software run's but the route's own, and its
fallbacks became calls one for one (4,341 and 60 at binary64 on the RX
7600). Forced on, `indexed` and `pair_gather` are refused by name, where
the software backend runs the same command. `--door-selftest` on
`hip:0` gives "116 identical, 0 differ, 48 not compared, 0 accounting
errors; controls 2 of 2", what the software backend gives with those two
routes hidden.

**What runs where,** hw/gpu/where.py under gdb on the RX 7600 at
`8ea24bd`, kepler with a rejected attempt at binary64, binary128 and
binary256, alone on the GPU:
- **The clock's augmented sum is the host's, at every attempt.**
  `cft_augmented_add` was entered once a step attempt, rejected ones
  included: 5 of 5, 4 of 4 and 3 of 3. None of the run's 19,351, 25,105
  and 20,941 calls into the GPU backend was made inside it.
- **The divide and square root take the program route.** `cft_div` and
  `cft_sqrt` were entered exactly `divsqrt_calls` times, 3,346, 4,412 and
  3,696, each running one program on the device and no elementwise run,
  with the host's rounding inside (11,523, 15,315 and 12,791 calls of
  `cft_sf_round_pack`). `CFT_DIVSQRT_SEQ=1` gives the same counts; `=0`
  gives no program and 72,917, 105,256 and 95,819 elementwise runs. The
  three records are the same at each format.
- **No reduction reaches the device by default.** The corrector's maxima
  are `cft_class` and a host loop unless `reduce_seg` is forced on
  (src/ias15_cft.c:1835). P1 found no backend reduction call; V1 found no
  entry to `cft_reduce` or `cft_reduce_seg`.
- **The index-table routes** read `by=absent`, with no call.
- **The programs.** All 69 images the port loads, programs/out's 66 and
  the Makefile's three normalabs images, load on both GPUs. The widest,
  gravity-fold256-fp128, uses all 256 scratch slots, the device's whole
  depth: it fits with no margin. The most constants any program
  addresses is 13. A block of 256 slots and a bank of 512 constants
  load; 257 and 513 are refused by name.
- **On the RTX 5060 Ti,** which has no gdb, the route is read from the
  forcing variable and the times alone, and the three records are the
  same each way.
- **What the device line proves** (V1, under gdb on the RX 7600, kepler
  at binary64). 15,089 `cft_run` calls became 15,089 kernel launches,
  and `cft_class`'s 79 calls and `cft_augmented_add`'s 5 reached no
  backend function: 84 of 19,434 counted calls are the host's by design.
  The line is printed once the device is selected, its kernels read and
  its stage allocated (backend_hip.cpp:1743-1800). So it proves that the
  process opened the GPU, not that every operation ran there: a libcft
  that ran everything on the host would pass the witness bit for bit.

**What a call costs,** hw/gpu/costs.py, MEASURED by P1: one run a row,
the trailer's `wall_seconds`, unranked within the noise, and every
`hip:0` row's record its software row's. kepler is P0's smoke (entry
54); the ring is one 64-body wave with its zonal line, 8, 6 and 4 steps
at the three formats. Each cell gives binary64, binary128 and binary256,
separated by semicolons, the loop engine in REBOUND's order and then the
program engine in the FMA form. "Over software" is the `hip:0` run's
wall time over the same host's software run's.

| GPU | kepler, µs a call | kepler, over software | ring wave, µs a call | ring wave, over software |
|---|---|---|---|---|
| RX 7600, at `20dad58`, 2026-11-16 10:42 to 10:54, the GPU idle but for these runs | 288, 333; 346, 476; 449, 690 | 18.0, 5.7; 17.3, 5.6; 14.4, 7.8 | 274, 281; 322, 346; 430, 457 | 0.6, 0.5; 0.5, 0.4; 0.4, 0.3 |
| RTX 5060 Ti, at `ec66a15`, 2026-11-16 15:46 to 15:51, the GPU otherwise idle, its host the shared desktop | 119, 140; 130, 176; 150, 229 | 11.6, 3.4; 10.0, 2.9; 8.2, 3.3 | 141, 157; 161, 187; 197, 230 | 0.5, 0.4; 0.4, 0.4; 0.3, 0.3 |

- A call is `wall_seconds` over the trailer's `calls=`, which counts the
  host's augmented sums too, one an attempt.
- P1's kepler rows on the RX 7600 are P0's runs again: 31.93, 8.00,
  87.02, 19.15, 155.96 and 35.18 s against P0's 31.6, 7.98, 85.9, 18.8,
  154.0 and 35.1, so +0.23% to +1.86%, the largest at binary128 on the
  program engine. P1 first wrote "within 1.5%"; V1 recomputed all six.
- The trailer's CPU `seconds` over its `wall_seconds` is 1.54 to 1.92 on
  the RX 7600, and 1.02 to 1.03 on the RTX 5060 Ti in all 12 rows. V0's
  1.5 to 2.0 (entry 54) is ROCm's.
- A divide or root, from where.py's runs: on the RX 7600, 661, 868 and
  1,341 µs by default, 651, 866 and 1,344 forced to the program route,
  and 4,603, 5,882 and 6,367 forced to the chunk route; on the RTX 5060
  Ti, 242, 299 and 386, then 261, 296 and 389, then 1,818, 2,091 and
  2,456. V1: the 5060 Ti's default and program-route figures at
  binary64, 242 and 261 µs, differ by 7.9%, more than CLAUDE.md's 6%
  between two runs, and that GPU's part is timing only, labelled so.

### The legs on each GPU

Every leg's own verdict PASS, every process on the device with that
GPU's device line, and every record the software backend's. Process
counts are the same on both GPUs. The RX 7600's times, at `8ea24bd`
under strace, are not costs: each includes the witness's software
replay, and from 12:40 a second process, check_dropin_light, shared the
GPU (P1).

| leg | processes | RX 7600 |
|---|---|---|
| check_records | 4 | 30m54s |
| check_equivalence | 6 | 6m12s |
| check_program_engine | 14 | 11m42s |
| check_ensemble | 49 | 15m38s |
| check_ring | 15 | 4m06s |
| check_replay | 141 | 20m30s |
| check_compare_run | 774 | 161m48s; strace: 774 runs of the witnesses and 1,548 of the real program |
| check_zonal | 31 | 31m29s |
| check_reversal | 15 | 15m06s |
| check_dropin_light | 1 | 147m53s, beside the main run; 170 ok in its own stdout |
| check_corrector | 165, 2 of them refused by their own arguments | 192m16s (V1, at `ab75fbe`); strace: 165 runs of the witnesses and 330 of the real programs |

The RTX 5060 Ti ran the same legs at `20dad58`. check_corrector's own
verdict there was PASS, 212 ok, and 165 of its 165 records were the
software backend's, but the witness FAILed it for 2 processes with no
device line: check_corrector's own refusal cases, `--trace-pc` or
`--record-iters` beside `--replay-iters`, which ias15_cft refuses with rc
3 before it opens any device, as it does in software. The rule was the
driver's: `ab75fbe` counts such a process apart, and V1's live runs at
`ab75fbe` passed check_corrector on both GPUs, the two refusals the same
two on each. ROUND5.md had check_corrector "through the registered
integrator"; only its `--library` step is (P1). That step runs
`check_dropin --light` with the exit test forced on. It took 62 min alone
on the RTX 5060 Ti (P1, at `20dad58`, 2026-11-16) and 66 min there (V1,
at `ab75fbe`, 2026-11-18), each time with its host, the desktop, shared.

**Which legs reach a GPU, beyond the eleven** (P1's `--reach` on the RX
7600, and V1's on the RTX 5060 Ti):
- **gate_subprocess** reaches it and passes, 31 ok: its 8 ias15_cft
  children each print the device line, and their 8 records are the
  software backend's.
- **gate_real** `--fp64` reaches it and gives the software suite's
  wide-state digests, `e08541bc0af9b3ec` and `56bacf3bf63f19a7` (entry
  45's), and FAILs only its identity check: `identity_ok()`
  (tests/gate_real.c:326) takes any artifact that is not `cft://` for a
  card image, and a GPU loads none. It is called only at :646, in
  gate_real's self-contained mode.
- **check_checkpoint** runs gate_real in its `--straight`, `--save` and
  `--resume` modes alone, which never call it, and passes on `hip:0` as
  it is: `--quick`, 7m56s, all 18 gate_real processes with the device
  line (V1).
- **check_seams** `--quick` fails there, 12 ok and 24 FAIL. Every FAIL is
  a forced `indexed` or `pair_gather` case, or a count of them, which the
  device refuses by name; its `reduce_seg` and `scalar` cases pass (V1).
- **check_archive's four gates** open no device: three open the software
  backend by design, and gate_stock links no libcft.

P1 first read gate_real's refusal as binding check_checkpoint and
check_seams too. V1 found by reading that it binds gate_real's
self-contained mode alone, and then measured the two above (the lead's
ledger, 2026-11-16 17:40 and 2026-11-18 12:14).

### The verifier, five times

**On `ab75fbe`, interim (2026-11-16, paused at Logan's word).**
Confirmed: no leg bypasses the shadow, by reading, and control g caught
a bypass; the set-aside rules and the refusal rule, above; the FMA form
really run; `make -n -B` of all, programs, check, check-quick and
check-light the same byte for byte with `5b1cac6`'s Makefile and
`ab75fbe`'s, `HIP` unset, where a copy with one word added to a recipe
differs; where.py's counts reproduced exactly at binary64; `cft_reduce`
never called; the clock's sum at every attempt. Found: check_corrector
had never run under strace; on the 5060 Ti nothing catches a bypass;
any `t_*` trailer key was set aside; stderr never fails a leg; a refused
process counts as a device process; each seam plant is one digit on one
line, not two, and every plant sat in wave 0, so a compare_run.py that
never compares member 1 would pass (read, not run); gate_real's refusal
does not bind check_checkpoint or check_seams; "within 1.5%" is 1.9%;
and the 5060 Ti's divide route rests on timing alone.

**On `ab75fbe`: READY** (2026-11-18). `make HIP=1 check-gpu` passed on
both GPUs:
- the RTX 5060 Ti, with `GPU_SKIP` leaving check_records and
  check_corrector: "check_gpu: 36 ok, 0 FAIL, in 123m04s on hip:0", the
  probe, controls a to f, the routes at three formats and the seam in
  both forms among them, control g not run;
- the RX 7600, check_corrector alone under strace: 212 ok, 192m16s, every
  exec accounted for.

Found, MEASURED: every seam plant sat in wave 0, and a one-line
compare_run.py that never compares member 1 passed the seam in both
forms; on a wave-1 plant the tree's FAILs member 1 and the copy calls it
IDENTICAL. Latent: TIMING set aside any `t_*` key, so a planted
`t_new=1` passed; stderr never fails a leg; a refused process counts as
a device process; on the 5060 Ti nothing can see a bypass.

**The lead sent P1 back** on that READY, narrowly, before the merge.
check-gpu is in no tier, so no desktop suite runs check_gpu.py, and a fix
after the merge would cost its own verifier and suite; and the gap is
CLAUDE.md's planted-control shape (the lead's ledger, 2026-11-18 12:14).
`a068acc` plants both waves: five columns on each wave's first and last
sampled lines in REBOUND's form, and x of body 1 in each wave in the FMA
form. `6a5e448` names the 13 timing keys and adds control a3. On the RTX
5060 Ti: the tree 26 ok, 0 FAIL, 7m24s; each member mutant 11 ok and 6
FAIL, its wave's six plants; P3's compare_run.py in the tree's place, 17
ok, 0 FAIL.

**On `6a5e448`: NOT READY.** No plant's verdict asked anything of the
planted wave's own particles. A compare_run.py whose `particle_compare()`
returns IDENTICAL without comparing scored 16 ok, 0 FAIL on P1's harness,
as the tree does. Confirmed: the 13 keys are exactly the trailer's
timings (src/ias15_cft.c:4773, :3281-3295 and :4782), no counter is
named, and a name must be followed by "=".

`fc4ba21`: each plant's verdict holds the planted wave's sampled
particles by what their replays carry, read from compare_run.py's own
lines. On the harness: the tree 16 ok; each member mutant 10 ok, 6 FAIL;
the particle mutant 6 ok, 10 FAIL.

**On `fc4ba21`: NOT READY.** No plant landed on a sampled particle's
own body. A copy whose `particle_compare()` skips it,
`enumerate(bodies[:-1])`, passed at 16 ok, 0 FAIL, comparing 144 values
where the tree compares 162.

`0f324bd`: the own-body plant in each wave, a class per plant with its
FAIL count, the coverage verdict, and a class P1 added unasked, the
trailer's counters, which both comparisons hold and no plant reached.

**On `0f324bd`: READY** (2026-11-18). The classes are derived from
compare_run.py's lines; the coverage verdict cannot pass with a class
missing; no class a replay or the member comparison carries is left
unplanted; V1's four mutants, its own edits, FAIL with the same sets on
both compare_run.py's.

Each look found the next class the plants missed: a wave, the planted
wave's particles, a particle's own body. Inside the classes, dt_last, E,
t_hi, the trailer's other counters and the bodies' other values stay
unplanted, as V1 accepted.

### The merge, the seam on the merged tree, and the suite

**The merge** (the lead's ledger, 2026-11-18 15:46). `round5-int` at
`4a8fe99`, which was main and origin/main, clean; `git merge --no-ff
0f324bd` gave `fda3b40`, with no conflict. `git diff --stat 4a8fe99
fda3b40`: the Makefile, hw/gpu/costs.py, hw/gpu/program_caps.c,
hw/gpu/where.py and tools/check_gpu.py, +2131 and nothing removed.

**The GPU seam on the merged tree,** the lead's, in WSL cft-gpu on the
RTX 5060 Ti, from a bundle at `fda3b40`, `make clean` first, check-pins
ok. Its device line:

```
cft hip: hip:0 = NVIDIA GeForce RTX 5060 Ti (sm_120, 36 multiprocessors, bus 0000:05:00.0, 16283 MiB), HIP platform nvidia, runtime 12090, driver 13040; of 1 visible, CUDA_VISIBLE_DEVICES=1; this libcft: HIP_ARCH 120, stage 64 MiB, page 256 bytes, kernels 40/40/64 registers and 4728/4208/12624 bytes local a thread (run/level/seq)
```

`check_gpu.py --artifact hip:0 --hip 1 --expect "RTX 5060 Ti" --no-legs
--keep --parts controls,seam`, 15:47:31 to 16:19:49: "check_gpu: 38 ok,
0 FAIL, in 32m16s on hip:0", RESULT: PASS. Controls a to f were caught,
and g was not run. In each form both waves rejected once, and the merged
tree's compare_run.py, P3's, held both members and the 4 particles
(bodies 22, 30, 57 and 31) IDENTICAL; all 22 plants FAILed as their
classes require, and both coverage verdicts passed. So check_gpu.py's
patterns read the merged compare_run.py's lines.

**The suite,** read from its log (int-check-fda3b40.log): 16:22:07 to
19:48:06 on 2026-11-18, 3 h 26 min on a shared host, P2's check_wave
runs and mutants beside it under B0171 until 17:35; clean, build and
check each rc 0; 25 legs; 1,504 lines beginning `  ok `; 17 leg RESULT
lines, all PASS; no line beginning FAIL. Every leg has `4a8fe99`'s ok
count (entry 57), check_compare_run
108; the build's 5 compiler warnings, as there; entry 45's six digests
present. check-gpu is in no tier, so 25 legs equal to `4a8fe99`'s is
what P1's diff should give.

P1's own suite, on `a29f602` (the final Makefile but for one comment),
10:34 to 14:23 on 2026-11-16, a shared host: rc 0, 25 legs, 1,463 ok,
every leg at P0's suite's count (`d7b5047`).

### Stated limits, not closed

The merge message's three open notes, and what the verifier found:
- **stderr never fails a leg.** `analyse()` keeps a difference and
  `witness_verdict()` prints it as a note: a device run that differs from
  software in stderr alone passes.
- **A refused process counts as a device process** (check_gpu.py:562 at
  `fda3b40`). "N of M process(es) on the device" subtracts the refused,
  but "165 of 165 record(s) identical" counts them, each as identical,
  and a leg whose processes were all refused passes the driver's "no
  process ran on the device" check; only its own verdict would fail it
  (V1).
- **On the RTX 5060 Ti nothing can see a bypass.** cft-gpu has no
  strace, so control g and the per-leg exec accounting do not run there.
  The evidence is the same code and the same process count a leg as on
  the RX 7600, where strace accounts for every exec. V1's suggestion, not
  built: run the real programs from a copy only the witness knows, so
  that a bypass fails instead of passing.
- **`--analyse` is weaker than a live run:** it checks no rc and takes
  the device line from the processes themselves (V1).
- **`GPU_LEGS` turns the controls, the routes and the seam off;** a
  subset that keeps them is spelled with `GPU_SKIP` (V1).
- **hw/gpu/where.py's `records_same()` sets aside any `t_[a-z_]+` key**
  (V1): a measurement tool, not the witness.
- **`--door-selftest` with nothing hidden** prints 191 lines on `hip:0`
  against software's 206 (V1). No leg runs it through the witness.
- **gate_real's self-contained mode refuses `hip:0`** on its identity
  check; the fix is a tests/ change, not made. check_seams fails on
  `hip:0` on the routes the device does not publish.
- **hw/gpu/costs.py's header** (lines 18-19) still gives V0's 1.5 to 2.0
  as a GPU's ratio of CPU time to wall time, which is ROCm's (above). The
  rows it prints carry the ratio measured beside each. A comment in
  code, not changed.

### What this does not establish

- **That every operation ran on the GPU.** The device line proves the
  process opened it (V1, above).
- **A full `make check` on `hip:0`.** Eleven legs are held; gate_real,
  check_seams and the legs above are not.
- **The final commit's legs on a GPU.** P1's main runs were at `20dad58`
  (the 5060 Ti) and `8ea24bd` (the RX 7600); V1 ran `ab75fbe`'s
  check-gpu with two legs; the seam as it landed ran live once, on the
  merged tree, on the RTX 5060 Ti alone.
- **cft-fp256's GPU round 2.** All of this is at the pin `eae74e3`.
  augmentedAddition on the device, index tables and a deeper scratch
  block reach this tree only through a pin move, which is Logan's call,
  after the round (the lead's ledger, 2026-11-17 09:44 and 2026-11-18
  05:41).
- **The card, and the MI300X.**

## 56. Round 5's P2: the wave kernel, the port's step at the hardware's binary64, its records the port's bit for bit on both GPUs; check_wave sent back seven times, four on V2's word, two on P2's own disclosures and one on the merge's GPU check

**2026-11-16 to 2026-11-20.** The desktop, shared, for the host twin,
check_wave's plants and mutants and P2's suites; the RX 7600 on
nvidia-box (HIP 7.16 with AMD clang, gfx1102) and the RTX 5060 Ti in WSL
`cft-gpu` (HIP 7.2.4 over nvcc 12.9, sm_120), each window booked in the
TimeSheet. Pins as in entry 53. P2 branched from `5b1cac6`. Its commits
are nine to `82b75b8`, then seven send-backs: `8f559c8`, `f3ef701`,
`1e13456`, `ff7221a`, `2d062cd`, `9fa5956` and `3d8ab5b`. `61e45bc`
merges `3d8ab5b` on `round5-int` at `620c98b`, after a first merge of
`9fa5956`, `7444157`, failed its GPU check and never reached main
(entry 58). main fast-forwarded `620c98b..61e45bc` and was pushed on
2026-11-20.

### What P2 is, in one paragraph

`src/wave/` is the port's step under `--arith fma` at binary64, written
once as `__host__ __device__` functions over `double`:
- **wave_arith.h:** add, sub, mul, div, sqrt and fma, each checked for
  the four flags the port aborts on, INVALID among them; the clock's
  augmentedAddition, with libcft's ties toward zero; and on an x86-64
  host, fma and sqrt through their intrinsics, so that every
  optimisation level gives the same bits.
- **wave_step.h:** the predictor and the corrector's passes, gravity with
  the `active` split, J2 and J4 in the port's order, the step control,
  the compensated sums and the energy.
- **wave_gpu.inc:** one system a workgroup, one thread a body. Every
  phase ends in a barrier, thread 0 folds the energy, the reactions and
  the maxima serially, and no wave intrinsic is used (V2, READ). A sample
  block is many short launches, so none meets the display driver's
  watchdog (`7dec985`).
- **wave_main.cpp:** the port's problem file in, and this repository's
  record format out, naming its producer: `program=wave`, `engine=wave`,
  `backend=host` or `hip:<index>:<arch>`, and on a GPU a `# device` line.
  It writes each member's attempts record and corrector schedule under
  the port's `--run-key` (`--attempts-dir`), and a flag stops its system
  by name. It refuses by name what the port has and it has not: any
  format but fp64, REBOUND's form, `--cs augmented`, adaptive modes other
  than PRS23, softening, a step floor, dt files, a system of more bodies
  than it holds, and a non-finite input value (READ, its head).
- **wave_corpus.cpp:** gate 4's corpus program.

`make wave-host` builds the host twin, its two contraction controls and
gate 4's corpus, in every tier. `make HIP=1 wave` builds the kernel,
`build/wave`, with the GPU's controls and corpus: hipcc for an AMD GPU,
or HIP's NVIDIA layer. `tools/check_wave.py` holds the result: a leg in
all three tiers on the twin, and check-gpu's second line on the GPU
(entry 58).

### The gate, and the controls that show it can fail

check_wave's six parts, each asserted (READ, its head):
1. **Gate 4, the arithmetic, measured** against libcft's software
   binary64: add, sub, mul, div, sqrt, fma and augmentedAddition over
   random operands, ties, subnormals and extremes, every result's bits
   and the four flags. At n 20000, 1,064,352 cases: 0 mismatches on the
   desktop, the RX 7600 and the RTX 5060 Ti, with one operand digest,
   `3cd1af43b434860c`, on every build (V2). Among them are 152,735
   augmentedAddition cases, 7,259 of them ties on which round-to-even and
   the clock's ties toward zero differ (int-check-61e45bc.log).
2. **The constants,** against the port's `--dump-constants`.
3. **J2 and J4** at the port's own `--zonal-selftest` points: every
   acceleration and both potentials are the port's bits.
4. **The records,** through compare_run.py `--rerun-engine loop`, every
   member sampled and its particles replayed:
   - A, a ring of 3 waves with its zonal line;
   - B, the same from a first step far too large, every wave rejecting
     its first attempt;
   - C, a particle aimed past Janus, rejecting after accepted steps;
   - D, kepler_e099, zonal_saturn and kepler at fixed steps, with outer
     and pythagorean in full;
   - E, one wave run alone;
   - F, A's ring with every test particle given a mass of 2^-20, which
     must move no bit.

   Each run's attempts records and schedules must be the port's own,
   under its `--run-key`, and on a GPU every record must also be the host
   twin's, line for line. Beside each go the steps, rejections, passes
   and cap hits, and for a ring the physics, asserted: each particle's
   distance from Saturn within 1e-2 of what the field's absence would
   move it, the angle swept within 1e-4 of n t, and the energy within
   1e-13.
5. **The controls,** below.
6. **The flags.** A problem made to raise each of INVALID, DIVBYZERO,
   OVERFLOW and UNDERFLOW stops the port and the program, rc 3, naming
   the flag. With `--no-flag-abort` both run on and report the same
   flags.

**The controls, each asserted to fail by name:**
- **Gate 4 can fail** (V2, desktop, n 4000, every plant named by its
  operation): a divide by a reciprocal, 3,400 mismatches; a square root
  by a reciprocal, 4,729; product subnormals flushed, 8,250; subnormal
  inputs to add read as zero, 8,156; augmentedAddition as TwoSum, 1,501;
  fma as a multiply then an add, 12,981; underflow decided one exponent
  early, 4,868 flags. On the RX 7600, the corpus built with
  `-fdenormal-fp-math=preserve-sign` gives 65,602 mismatches of 228,352,
  in all seven operations.
- **Contraction.** Every shipped object is built with contraction off,
  and `--fmad=false` under nvcc. The control that can fail compiles the
  flag checks out: with contraction off that source passes compare_run,
  and with a*b + c fused it FAILs. V2's counts of fused instructions:

| compiler | shipped | the checked source, contraction on | the unchecked source, off, then on |
|---|---|---|---|
| MinGW g++ 16.1 | 10 vfmadd, all from `_mm_fmadd_sd` | 11, the extra one in the timer | 3, then 44 |
| hipcc, gfx1102 | 632 v_fma_f64 | the same assembly, line for line | 403, then 597 |
| nvcc, sm_120 | 1,169 DFMA | the same SASS, 56,840 lines | 848, then 1,084 |

  So the flags reach every compiler. With the checks in, each product
  feeds its check and none fuses: V2's belief about the mechanism, the
  effect measured. One product-sum fused by hand, gravity's x
  accumulation, FAILs compare_run in member 1 at its last sample (V2).
- **The plants:** one value a plant, in 14 copies of an ensemble record
  (below).
- **The attempts comparison:** a member's attempts record with one dt
  changed must not be the port's.
- **The physics:** wave 0 of A's ring with its zonal line deleted has
  sample lines that are the port's bit for bit, and physics_verdict()
  must refuse it, on the distance and on the angle.

**What the gate found.**
- **The clock is not TwoSum.** libcft's augmentedAddition rounds an exact
  tie toward zero. A twin whose clock is plain TwoSum fails the corpus,
  and writes records unlike the port's in two of the three cases P2
  tried (2026-11-16). ROUND5.md is corrected.
- **The flags.** The port aborts on INVALID too, which ROUND5.md's list
  of flags leaves out. The kernel checks all four (P2).
- **nvcc quiets a signaling NaN in a sign flip.** At `bc9feb3`, gate 4 on
  the 5060 Ti found 3 mismatches in 1,064,352 cases, every one a
  subtract with a signaling NaN as b. nvcc compiles the bit-level
  negation as `DADD R6, -RZ, -R2`, and the abs mask as an add too.
  `81af26a` has the subtract's slow path read b itself, and no finite
  value moves (P2).
- **The corpus depended on the compiler.** C++ leaves the order of a
  call's argument draws to the compiler: g++ draws the last first, clang
  the first. So the RX 7600's first gate 4 ran another corpus, though
  with 0 mismatches. At `82b75b8` every draw is in g++'s order, written
  out, and the program prints an operand digest, which check_wave holds
  GPU build against host build (P2).
- **MinGW's library fma and sqrt are not correctly rounded, and g++
  calls them at -O0.** Built so, the corpus found 20,461 of 178,342 fma
  results and 14 of 100,051 square roots that are not libcft's:
  sqrt(0x1.fffffffffffffp+1023) is 0x1p+512 there, against libcft's
  0x1.fffffffffffffp+511 (P2). V2 measured the library fma differing
  from vfmadd on 18,013 of 199,714 random finite triples. Every Makefile
  build is -O2, and since `1741237` the twin takes both instructions
  through their intrinsics at every level. BELIEVED (P2): REBOUND's own
  C, ias15_ref, built at -O0 here would take the library's root.
- **A transcribed bound.** P2's first full make check, at `7dec985`,
  failed one verdict of check_wave's own: a 1 km physics bound, copied
  from round 4's seam figure for another run. A particle that passed 679
  km from Pandora moved 1.7622 km. `bc9feb3` derives the bounds above:
  1.43e-3 at worst against 1e-2 for the distance, and 9e-6 against 1e-4
  for the angle (P2, 2026-11-16).

**On both GPUs** (MEASURED):
- **The records.** On the RX 7600 at `82b75b8`, check_wave `--device gpu`
  passes, 78 ok in full and 66 in `--quick`; on the RTX 5060 Ti at
  `81af26a`, 77 and 65 (P2). V2 reran the commit of its first look: 78
  ok in full on the RX 7600, 66 in `--quick` on the 5060 Ti. At
  `8f559c8`, `--quick` passes, 91 ok, on the RX 7600 (P2) and on the 5060
  Ti (V2). Every record is compare_run-IDENTICAL member by member, its
  particles replayed, and line for line the host twin's.
- **Across GPUs** (V2). P2's two 4,096-wave records of 362 steps, one
  from each GPU, differ in 0 of 278,535 lines, with `backend=`, the
  identity and device lines and the trailer's five timings set aside. A
  digit planted in a copy shows as one line.
- **Determinism** (V2). A 3-wave ring run twice and a 527-wave ring run
  three times on the RX 7600, and a 2-wave ring run twice on the 5060
  Ti, are identical but for their timings. A build with
  `-mwavefrontsize64` gives the wave32 build's records.
- **Three ways.** One ring wave from a first step of 1e5 s, through a
  first-attempt rejection: the port on the CPU, the port on `hip:0` and
  the kernel give the same lines and counters, on the 5060 Ti and on the
  RX 7600 (P2; V2 again on the RX 7600).
- **The kernel's frame** (P2): on gfx1102, 120 VGPRs, 144 bytes of
  scratch a lane and 22,544 bytes of LDS a workgroup; on sm_120, 182
  registers, a 112-byte stack and 23,568 bytes of shared memory.

### What a system-step costs on each GPU

P2's rates, MEASURED, one run a row, unranked within the noise:
make_ring.py's ring (seed 1) with its zonal line, a first step of 600 s,
epsilon 1e-9, one sample at the end. Two clocks, both from the record's
trailer: the kernel's own timer around each launch through its
synchronize, the first launch left out, and `wall_seconds`, the run
loop with the record's writing. System-steps a second, the kernel clock
then the wall clock:

| waves | steps | RTX 5060 Ti, at `81af26a`, 2026-11-16, its host the shared desktop | RX 7600, at `82b75b8`, 2026-11-18, nvidia-box quiet |
|---|---|---|---|
| 1 | 24 | 857; 828 | 800; 774 |
| 16 | 24 | 12,800; 12,000 | 11,636; 10,971 |
| 128 | 24 | 53,895; 42,667 | 40,960; 33,758 |
| 527 | 24 | not run | 49,600; 40,280 |
| 1,755 | 24 | 63,529; 49,437 | not run |
| 4,096 | 24 | 64,504; 48,473 | 54,222; 42,909 |
| 527 | 362 | not run | 57,016; 56,028 |
| 1,755 | 362 | 67,143; 65,883 | not run |
| 4,096 | 362 | 68,595; 67,288 | 59,315; 58,255 |

- **Beside GPU-E2.** Its estimates, which priced libcft's software
  binary64, were 13 minutes for 3 x 10^4 particles over 10 orbits on the
  RX 7600 and 14 minutes for 10^5 on the 5060 Ti: about 245 and 756
  system-steps a second (COMPUTED by P2, as 527 x 362 / 780 s and 1,755
  x 362 / 840 s). The kernel, in the hardware's binary64, ran those
  rungs, 527 and 1,755 waves of 362 steps, 10 orbits at GPU-E's 36.2
  steps an orbit, in 3.35 s and 9.46 s of kernel time, and 3.41 s and
  9.64 s of wall time.
- **The rungs' records held.** compare_run reran 2 members of each,
  w00070 and w00301 of the 527 and w00301 and w01151 of the 1,755, with
  their particles: IDENTICAL (P2).
- **One row again** (V2): 527 waves of 24 steps on the RX 7600, three
  runs, 0.257, 0.254 and 0.254 s by the kernel clock against P2's 0.255,
  and 0.317, 0.313 and 0.314 s by the wall clock against 0.314.
- **The physics beside the rates,** asserted on none of these runs: the
  energy moved 3.5e-15 and 4.1e-15 of itself over the 5060 Ti's 1,755-
  and 4,096-wave rungs, and 3.68e-15 over the RX 7600's 527. The two
  GPUs' physics agree to every printed digit on 16 waves sampled every
  10 steps, and on the 4,096-wave rung (P2).

### The plants, and the verifier, six times

**The plants as they landed, at `3d8ab5b`,** in an ensemble record of
the program under test: P, the ring of 3 waves at 12 steps, or A itself
in `--quick`. One value a plant, in 14 copies:
- **copies 1 to 11** (send-backs 2 and 3): every class `record_classes()`
  reads from the record, in every member. These are the schedule (t,
  dt_next, dt_last, E, t_hi, t_lo), a massive body's x to vz, a sampled
  particle's own x to vz, the step field, and the trailer's counters,
  the sample-line classes each on the member's first and last sampled
  lines. One plant a member a copy, and each member is left out of one
  copy;
- **copy 12** (send-back 4): the last massive body, 6 (janus), in every
  member;
- **copies 13 and 14** (send-back 5): sample 0 in every member, a
  schedule column (t, E, t_lo), then the step.

Each copy must hold exactly its planned values against the clean record
(`measured_plants()`). Each plant must FAIL its own member, naming its
line and column, and exactly those of the member's particles whose
replay carries it, while every other member stays IDENTICAL. One closing
verdict (`coverage()`), from its own parse of the record, names any
class, member, body or line left unplanted. And every IDENTICAL member
and particle line of every compare_run run must report the record's own
count of that member's sample lines, and that many lines' worth of
values (`hold_counts()`: 390 a member line and 54 a particle line for the
ring).

**On `82b75b8`: NOT READY** (V2, 2026-11-18). The kernel held: its bits
on both GPUs, gate 4 and its plants, contraction, determinism, the
Makefile, and the P2 x P3 and P1 x P2 seams, checked ahead. The planted
controls did not. Every plant sat in case E, one wave run alone, on
sample lines 1 and 2, so compare_run copies that compared less passed
check_wave, 49 ok in `--quick` and 57 in full. They never compared:
- ensemble member 0's lines, member 1's, member 1 at all, or member 2
  (full);
- any sample line past the third;
- any particle;
- ten of the columns;
- the trailer's counters, or the step field.

So did check_wave copies without their physics decision, attempts
comparison or particle parse. With one digit planted in wave 1's last
sample line, the tree's compare_run returns rc 1, and the never-member-1
and first-three-lines copies rc 0. check_compare_run `--quick`, in every
tier, fails the never-member-1 copy, 42 ok and 3 FAIL. Also:
- nothing asserted that a particle was replayed, so 0 of 0 passed;
- the physics control computed its own verdict;
- the identity line was a literal;
- the probe ran two kepler steps before it knew the device.

**Send-back 1, `8f559c8`:**
- the plants in the three-member ring record, 24 a run: every column
  kind, members 0 to 2, first and last lines;
- the particle count held to a plan, 6 for the ring;
- the physics control reading physics_verdict();
- one attempts comparison, with a one-dt control;
- the identity line written from the build's defines;
- a probe, `--device-info`, that launches nothing and refuses an
  unexpected device before anything runs.

The tree: 74 ok in `--quick` and 83 in full. All 27 of V2's mutants FAIL
in both forms (P2).

**On `8f559c8`: NOT READY** (V2, 2026-11-19). Four one-line edits passed
check_wave `--quick`:

| copy that | ok/FAIL | why it passed |
|---|---|---|
| compare_run checking the step in member 0 only | 74/0 | the step field was planted once, in member 0's last line |
| compare_run never checking the step on sample 1 | 74/0 | the same |
| check_wave never planting member 2 | 66/0 | the plan was 24 hand-written tuples |
| check_wave planting no counter | 69/0 | the same |

Everything else held: every mutant on the list FAILs both forms; the
particle plan; physics_verdict on both criteria; the identity line in
all three builds, on both hosts; and the probe. On the 5060 Ti,
`--quick` passes, 91 ok, and `--expect "RX 7600"` is refused on the
probe's line with an empty work directory. The step and steps_done turn
out to be held by compare_run's completeness check alone: without it,
check_wave FAILs 73/1 for each, and without the comparisons it passes
74/0.

**Send-back 2, `f3ef701`:** the classes read from the record
(`record_classes()`), every class in every member on its first and last
lines, 30 plants in 11 copies, each member left out of one copy, and a
closing coverage verdict. The tree: 81/0 and 90/0. All 31 mutants FAIL
in both forms: (i) 77/4, (ii) 78/3, (iii) 70/1 and (iv) 74/1 in
`--quick` (P2).

**On `f3ef701`: NOT READY** (V2). The coverage verdict read the plan,
not the copies, and four one-line edits of check_wave passed:

| check_wave that | --quick | full |
|---|---|---|
| drops the step class from record_classes() | 75/0 | 84/0 |
| drops max_pc_iterations from the counter class | 81/0 | not run |
| writes every last-line plant on sample 1 | 81/0 | 90/0 |
| moves the left-out copies past the end, copy 11 empty | 81/0 | not run |

With mutant (i) beside the first, it is still 75/0: send-back 1's gap,
open again behind one line. Held: the tree, the four mutants, a plants()
that skips the last line (80/1), and the left-out members (75/6, twice).

**Send-back 3, `1e13456`:**
- measured_plants();
- coverage() from its own parse, with record_classes() required to
  agree;
- the first and last lines found there, and required to differ;
- a member counted as left out only where measured;
- the particle class only on replayed bodies (P2's addition).

The tree: 92/0 and 101/0. V2's four FAIL by name in `--quick`: 83/1,
91/1, 87/5 and 90/2.

**Send-back 4, `ff7221a`,** from P2's own disclosure, ruled in scope by
the lead before V2 looked. The massive class sat on bodies 0 to 5, never
on 6, and cr_last_massive_short, whose massive slice stops one short,
passed `1e13456` at 92/0. A twelfth copy plants body 6 in every member,
and the coverage verdict requires every massive body, and the last in
every member. The mutant FAILs, 93/3 and 102/3. The tree: 96/0 and
105/0. The eleven copies are unchanged once the totals line's two
timings are masked.

**Send-back 5, `2d062cd`,** from P2's next disclosure, ruled likewise.
Sample 0 was never planted, and cr_sample0_skipped, both loops starting
one late, passed `ff7221a` at 96/0. P2 first swept every loop compare_run
compares over, naming its two ends and what holds each (READ). The two
open ends were at sample 0, and copies 13 and 14 with one coverage rule
close them: the mutant, and each loop started late on its own, FAIL
101/3 (110/3 in full). The tree: 104/0 and 113/0. Six other shapes were
reported, not fixed. Four stand as open notes below; the totals sum
proved held (V2, next), and the printed value count is send-back 6's.

**On `2d062cd`: NOT READY** (V2), on one point. A compare_run that skips
the middle sample line passed, 104/0 and 113/0, and one that stops after
four lines passed full, 113/0. compare_run prints its counts, 1170 and
162 values against the tree's 1560 and 216, and 1560 on A against 1950,
and check_wave read neither. Everything else held:
- V2's six third-look mutants FAIL by name;
- the four loop-end mutants FAIL, 101/3;
- 12 of V2's 13 mutants of the new code FAIL;
- the copies are unchanged across the four commits.

P2's READ that a totals sum stopping one short passes was wrong: it
FAILs, 103/1 and 112/1, on record B.

**Send-back 6, `9fa5956`,** the count rule:

| run | --quick | full |
|---|---|---|
| the tree | 130/0 | 144/0 |
| check_wave's member width one short | 118/12 | not run |
| compare_run skipping the middle line | 116/14 | 126/18 |
| compare_run stopping after four lines | 128/2 | 134/10 |

`--quick` holds the four-line cap through records C, of 13 lines, and
D1, of 5 (P2).

**On `9fa5956`: READY** (V2, 2026-11-20):
- both of V2's mutants FAIL in both forms, on count verdicts only;
- three narrower ones FAIL: no counts printed (116/14), a middle line
  skipped in the particle loop alone (118/12), and both loops cut to
  four lines with 4 samples reported (128/2);
- the particle width one short fails the tree (118/12);
- three weakenings of the rule each still FAIL both mutants: the counts
  taken from compare_run's output, the particle lines skipped, and the
  plant copies skipped;
- the 14 copies are unchanged.

**Send-back 7, `3d8ab5b`,** after the first merge's GPU check (entry 58).
check_wave's default work directory, `build/wave/<pid>`, was where `make
HIP=1 wave` puts the kernel, as a file. It is now `build/check_wave/<pid>`.
With a plain file at `build/wave`, `9fa5956` fails at `os.makedirs`
(check_wave.py:1088) before running anything, and `3d8ab5b` passes,
130/0 (P2). Every GPU run of check_wave by P2 and by V2 had passed
`--keep` outside build/ (READ).

**On `3d8ab5b`: READY** (V2, 2026-11-20). The same clash, on V2's own
setup, fails `9fa5956` and passes `3d8ab5b`, 130/0, the file untouched.
The tree passes `--quick` at 130/0 and full without `--keep` at 144/0,
leaving no work directory behind. The copies are unchanged.

The copies are always compared with the totals line's timings masked:
it carries the run's measured throughput, which differs between any two
runs (P2; V2 measured two runs of one commit differing there alone).

### The suites

- **P2's own full make check, at `82b75b8`:** 05:59:39 to 09:27:33 on
  2026-11-18, 3 h 28 min on a shared host. rc 0, 26 legs, the 25 at P0's
  suite's counts (`d7b5047`), and check_wave 57 ok, RESULT: PASS (P2; the
  lead's leg-by-leg reading).
- **At `7dec985`,** 11:55:25 to 15:51:55 on 2026-11-16: rc 2. 25 legs
  passed, and check_wave failed its own transcribed bound (above).
- The merge's GPU check and suite are entry 58's.

### Stated limits, not closed

The merge message's open notes:
- **Test-test pairs.** The port computes and discards a ring wave's 1,596
  test-test pairs, and the kernel never computes them, so two
  coincident test particles stop the port with DIVBYZERO but not the
  kernel.
- **The unchecked control holds a different source from the shipped
  one.** Its PASS shows that the checks move no bit in the arithmetic
  the two builds share. It does not show that the shipped build's checks
  would stop a fused a*b + c.
- **GPU_LEGS and GPU_SKIP do not restrict check_wave's line.**
- **A planted step or steps_done** is caught by compare_run's
  completeness checks, hold_structure() and scan(), which mark the
  member incomplete before any comparison. No plant reaches compare_run's
  own comparisons of those fields, or of a line's (sample, step) pair.
- **A plant compare_run refuses outright** (hold_initial(),
  hold_problem(), hold_rerun()) has no verdict in check_wave that
  expects a refusal.
- **The particle-selection loops** (closest(), least()) compare nothing.
  An off-by-one there replays a different particle, which check_wave
  plants just as faithfully.
- **No plant deletes a trailer key.** A trailer missing one still FAILs
  the counter loops (READ).
- **The sample-0 rule edited to look at sample 1 passes** (V2). That
  edit weakens the rule; every edit that drops or moves a sample-0 plant
  FAILs.
- **The copies are compared with the totals line's timings masked.**

And:
- **The kernel opens device 0,** the first the GPU runtime shows it
  (wave_gpu.inc, `hipSetDevice(0)`), whatever `GPU_ARTIFACT` names. Its
  probe holds that device's line to `GPU_EXPECT` before anything runs
  (the Makefile's check-gpu comment).
- **Binary64 only.** Binary128 and binary256 in a kernel need a divide
  and a square root inside it (ROUND5.md).

### What this does not establish

- **The ladder as ROADMAP's GPU half names it.** The rungs ran as rate
  rows, one run each, with two members of each rerun. Their records were
  not compared rung against rung; the comparison across GPUs was of the
  4,096-wave records.
- **Binary128 or binary256 on a GPU, the MI300X, or more than one GPU at
  a time.**
- **Any speed claim beyond the rows above.**

## 57. Round 5's P3: `tools/compare_run.py --producer-files` replays particles from a producer's own attempts records and schedules; they hold its arithmetic, not its decisions, and the driver now prints that; sent back twice

**2026-11-16 to 2026-11-18.** This desktop, on the software backend,
shared, every run booked in the TimeSheet; no GPU. Pins as in entry 53.
P3 branched from `5b1cac6`. Its commits are `2922fa1`; `0bee543`, after
V3's first look; and `e3f6a39`, after its second. `4a8fe99` merges
`e3f6a39` on `round5-int` at main `5b1cac6`, and main fast-forwarded
`5b1cac6..4a8fe99` and was pushed on 2026-11-18.

### What P3 is, in one paragraph

docs/RING.md had it as not built: a GPU run's particles sampled without
a rerun of their wave. compare_run.py takes, beside a producer's record,
each member's attempts record and corrector schedule:
- **`--producer-files DIR`:** member k's `DIR/member<k>.att` and
  `DIR/member<k>.sched`, k its index in the problem file and 0 for a
  single problem, in exactly the forms `--attempts-out` and
  `--record-iters` write for that member run alone. Each is keyed, on its
  `# run` line, with what `ias15_cft --run-key` prints for the arguments
  the member is rerun with: the record's own engine when this tree has
  it, and `--rerun-engine`'s when not. A single problem with no system
  lines keys its files `member=-1` and names them member0's (P2 found it,
  P3 measured it).
- **Before anything runs,** each pair is held to its member's run by that
  key. A file of another wave, format, `active` line or `zonal` line is
  REFUSED, rc 2, naming the member, the file and every key field that
  differs; a file not of its form is refused by its first line.
- **The particles sampled** are the member sample's members' and one
  member's from each of `--particle-strata` ranges of member index
  (default: `--strata`'s), drawn apart under a SHA-256 name of their own.
- **Each is replayed** against the producer's own files, with no rerun of
  its member. A member the member sample holds is still rerun; one it
  does not is NOT RERUN, held through its particles' replays alone, and
  the driver says which, and why. A member whose pair is not both there
  falls back to today's rerun, naming the file that is not.
- **Each particle line** names the files it replayed from, and gives its
  replay's counters and its closest approach at the samples. Its prefix
  and its `held to the record: IDENTICAL|FAIL` are as before. V3's first
  look ran P1's patterns and check_ring_seams part 3's on the new
  output, and they captured what they had. Its second ran those and
  P2's check_wave.py patterns on three cases at both commits, with the
  same result. Its third checked by reading: the driver's printed text
  had changed only in the limit's line, which holds neither IDENTICAL
  nor FAIL, and check_wave.py:361-364 were unchanged.
- **The limit,** printed once beside the plan under `--producer-files`,
  and on no other path, as the driver's `PRODUCER_LIMIT` states it:

  > particles: the limit of the producer's own files: each replay takes
  > the producer's dts, outcomes and passes as given, so only a member
  > that is rerun holds its decisions; a NOT RERUN member's particles
  > that are not sampled are held by nothing; and the run key holds the
  > files to the arguments the record names, not to those the producer
  > ran with, so its epsilon, pc_tol_shift and max_iter, and the settings
  > no record carries, are held only in a member that is rerun, and there
  > only through what they move in its lines

### What a producer's files hold: its arithmetic, not its decisions

MEASURED by P3, and the round's finding against ROUND5.md's P3. Without
a rerun, a replay takes the producer's decisions as given, each
attempt's dt and outcome and its passes, and holds the arithmetic that
reaches its samples and trailer: every sample's head fields and
exact-time pair, the particle, every massive body and the trailer's
counters. Only sample lines and the trailer are compared, and a rejected
attempt after an accepted step moves no state bit (entry 49): its dt
shows only as the dt_next of a sample line printed just before it. Take
a producer whose member 4 took 5 passes on its first attempt where this
tree takes 7, its record that schedule's replay, so that 2 of its 3
sample lines are not this tree's. With that member NOT RERUN, the driver
holds it, both its particles IDENTICAL. GATE 1b's run ends rc 1, since
member 2, rerun, FAILs (check_compare_run.py:836); where nothing else
differs, the run says RESULT: PASS, as V3 measured on a three-wave ring
of make_ring.py's, wave 0 at 2 passes where this tree takes 4. The same
member FAILs once it is rerun. So the decisions are held only where the
driver reruns a member: one member a stratum, and every member with a
rejected step.
- **Pass counts barely move bits.** At `2922fa1`, one pass more or fewer
  moved no bit at the samples in 7 of 7 tries at binary64, and two fewer
  on an accepted first attempt moved none at binary128 or binary256.
  Two fewer on member 0's rejected first attempt, planted since
  `0bee543`, moved 8 tokens at binary256 (q2, below). The trailer
  catches them: mean_pc_iterations, and on a first attempt
  max_pc_iterations too (V3; P3 had first written "only the mean").
  Since `0bee543` the gate derives what each plant moves and prints it.
- **The 3-decimal mean.** COMPUTED by P3: one pass is not certain to show
  in a member of more than about 1,000 attempts. COMPUTED by V3, from
  GPU-E2's 36.2 steps an orbit and docs/RING.md's Pan figures: a ring
  wave takes 359 to 374 steps over 10 orbits, one pass moves its mean by
  about 0.0027, and the limit binds only past about 28 orbits. Today's
  rerun path has the same blind spot where a pass off moves no bit.
- **A rejecting wave is always rerun.** Every member with a rejected step
  stays in the member sample, so a ring started from too large a first
  step is rerun whole, files or not.
- **The run key is a label for any other producer** (V3). One that keys
  its files with `--run-key` on the record's arguments passes the check
  whatever it ran with. Member 4 produced at pc_tol_shift 2 or 8, or at
  max_iter 12, its header rewritten to name the defaults, is IDENTICAL on
  its rerun; at max_iter 6, with two cap hits, or epsilon 1e-8, it FAILs
  (V3, and P3 at pc_tol_shift 2 and epsilon 1e-8, which FAILs at record
  line 13, t).

### The gate, and the controls that show it can fail

`tools/check_compare_run.py`, case 7, in every tier: at `e3f6a39`,
`--quick` 70 ok and full 108. The producer is ias15_cft, each of the 9
members run alone with `--attempts-out` and `--record-iters`, each run's
lines and counters asserted to be the record's member's. As the merged
tree's suite printed it, at every format case 1 runs:
- **From the producer's files,** every member's particles sampled
  (`--particle-strata 9`): PASS, 18 particles replayed from the
  producer's own files, all IDENTICAL. From an empty directory, each
  member falls back to a rerun, naming both files: PASS.
- **GATE 1, the same particles.** Each of the 18, replayed from the
  producer's files and from a rerun's, is IDENTICAL: every sample line,
  1,512 tokens, and the trailer's counters, flags and calls. Each
  replay's record names the attempts record it ran from, 18 of 18 each
  way, and members 0 and 1, whose first attempt was rejected, are among
  them.
- **GATE 2, the physics beside each.** Every particle's replay counters
  are its member's trailer, at binary64 40 steps, 0 to 3 rejected, 3.293
  to 3.875 mean passes, 7 to 9 at most, no cap hit; and its 36 closest
  approaches agree with the gate's mpmath to 3.29e-7 relative at
  binary64 and 3.14e-7 at the other two.
- **GATE 1b, not against itself,** at binary64. The producer above:
  member 4 not rerun, and member 2 rerun with one pass fewer on its first
  attempt, which moves no line and moves mean_pc_iterations, 3.275
  against 3.300, and max_pc_iterations, 7 against 8. From the producer's
  files member 4 holds, NOT RERUN, both particles IDENTICAL, and member 2
  FAILs on its rerun, its particles IDENTICAL (rc 1). Replayed from this
  tree's reruns instead, members 4 and 2 FAIL and all 4 of their
  particles FAIL. So the same-particles gate cannot pass against itself.
- **And:** the waves in their primary's zonal field, 18 particles
  replayed with the line; half a pair, the schedule alone named missing
  and that member falling back; a `--member 4` record; a single problem
  keyed `member=-1`; kepler, every body active, its files not read and no
  limit printed; the record rewritten as another producer's
  (`program=wave impl=hip-7.16 backend=hip:0 engine=wave`), PASS with
  `--rerun-engine loop` and REFUSED without; and the FMA form's
  program-engine record with `--rerun-engine loop`, its files keyed at
  `engine=program`, 4 steps in `--quick` and 20 in full.
- **The limit's line,** held as lines: exactly one, word for word, after
  the plan's line, and no other line opening with its words, in each of
  the 12 runs at binary64 that read a producer's files, 20 in full; and
  no line of it, in any wording, where no files are read: case 1 at each
  format, case 5's two records, and case 7's kepler.

**The controls,** each on a copy of the producer's files, deleted after;
q1 and q2 at every format:

| control | planted | result |
|---|---|---|
| q1 | one attempt's dt, a member a line: member 0's rejected first attempt, member 4's first, member 5's middle, member 7's `next` line; and member 1's outcome turned | each FAILs both its member's particles by name; the turned outcome is refused by the replay |
| q2 | one pass count: member 4's first attempt two fewer; member 5's middle one more; member 2's first one fewer, rerun; member 0's rejected first attempt two fewer, rerun | each FAILs both its member's particles by name. Printed, not asserted, at `e3f6a39` in full: 26, 0, 0 and 9 tokens of their replays moved at binary64, none at binary128, and 0, 0, 0 and 8 at binary256; the mean passes moved every time |
| q3 | the schedule, then the attempts record, from another wave | REFUSED, `member, state` |
| q4 | the zonal line's J2 changed | REFUSED, `zonal` and no other |
| q5 | the schedule recorded at binary128 | REFUSED, `format, pc_tol_shift, dt0, epsilon, state` |
| q6 | the schedule recorded under `active 3` | REFUSED, `active` and no other |
| q7 | the attempts record given as the schedule | REFUSED by its first line |

**Mutants of the driver,** each edit landing exactly once or refused,
against the `--quick` gate at `e3f6a39` (P3's runs; V3's matched for M1
to Mj):

| mutant | the driver | ok | FAIL |
|---|---|---|---|
| M0 | unchanged | 70 | 0 |
| M1 | replays rerun members from the rerun's schedule, naming the producer's | 66 | 4 |
| M2 | M1, and reruns every particle member | 65 | 5 |
| Mf2 (V3's) | the rerun's schedule for rejecting members only | 68 | 2 |
| Mj (V3's) | the files keyed to `--rerun-engine`'s engine | 69 | 1 |
| Ltwice, Lgone, Lword, Lbefore | the limit's line twice and adjacent; removed; one word changed; before the plan | 58 | 12 each |
| Ldefault | the line printed on the default path too | 59 | 11 |
| Grm | (a copy of the gate) one of gate 1's kept replays deleted | 69 | 1 |

### The verifier, three times

**On `2922fa1`: NOT READY, on D1 alone** (2026-11-16; the gate 105 ok in
full, 68 with `--quick`). The driver's bits were right, but its output
never stated the limit: only compare_run.py's source header did. The
output said the opposite, "NOT RERUN: its lines are held by the record
alone, and by its particles' replays", and V3 measured a producer whose
exit test stopped two passes early, in a wave not rerun, get RESULT:
PASS. Confirmed: the scope; the gate; the default path, where only the
particle lines gained their new clauses, on six setups; the key, each
wrong file refused, rc 2, before anything ran; the readers; E and the
exact-time pair held for NOT RERUN members through the replays. Also
found: the pass-count claim was wrong, five of eight plants moving
max_pc_iterations too; V3's mutant Mf2 passed, since q2 planted no
rejected attempt; Mj was caught only by the full tier's program-engine
case; the gate split paths at " and "; and a particle member costs five
whole-file reads by the code, not P3's three.

`0bee543` prints the limit once beside the plan, asserted verbatim;
rewords the NOT RERUN member line, its plan line and RESULT's suffix;
derives and prints what each pass-count plant moves; plants two fewer
passes on member 0's rejected first attempt; runs the program-engine
case in `--quick`, at 4 steps; and matches paths whole.

**On `0bee543`: NOT READY** (2026-11-18). The line printed twice passed
the gate, 69 ok: `out.count("\n" + LIMIT + "\n")` counted two adjacent
copies, which share a newline, as one. The line printed on the default
path passed too: nothing asserted it absent. And the settings clause
said more than a rerun holds: produced at pc_tol_shift 2 or 8 or
max_iter 12 and rewritten to name the defaults, a member came back
IDENTICAL on its rerun. Weak: 4 of the 9 producer-path runs asserted the
line, and `moves()` reported -1 tokens on a missing replay instead of
failing.

`e3f6a39` counts the line as lines, asserts it absent where no files
are read, ends the settings clause "and there only through what they
move in its lines", asserts the line in all 12 producer-path runs at
binary64, and makes `moves()` FAIL by name on a missing replay.

**On `e3f6a39`: READY** (2026-11-18). The two mutants that had passed
FAIL by name: two adjacent copies 58 ok and 12 FAIL, the default path 67
ok and 3 FAIL. The settings clause is true against V3's own runs; 12
runs assert the line at binary64 and 20 in full; the gate 70 and 108 ok,
the 3 over 105 being case 1's absence verdict at three formats. **One
open point, minor:** an indented copy of the line, or a print confined
to runs without producer files that the absence check skips, still
passes. Those runs are `--rerun-engine` without files, which is P2's
check_wave's own call, case 3's `--plan` and case 4's records. No reader
matches the line and no bit is at stake. The fix, `limit_absent`
asserted inside `drive()` for every run without `--producer-files`, was
queued for the lead.

### The merge and the suite

P3's own suites at `2922fa1`, on a shared host, 2026-11-16: `make
check-light`, rc 0, 25 legs, none skipped, 937 ok; and `make check`,
12:39:51 to 16:40:26, rc 0, 25 legs in full, 1,501 ok and 0 FAIL,
check_compare_run 105 ok (its 67 and case 7's 38). The lead read it
against P0's suite on `d7b5047`: every leg's count the same but
check_compare_run's.

**The merge** (the lead's ledger, 2026-11-18 09:42). `round5-int` made
at main `5b1cac6`; `git merge --no-ff e3f6a39` gave `4a8fe99`, with no
conflict. `git diff --stat 5b1cac6 4a8fe99`: the two tools, +1247 and
-58; `git diff e3f6a39 4a8fe99` is empty.

**The suite,** read from its log (int-check-4a8fe99.log): 09:41:33 to
12:16:29 on 2026-11-18, 2 h 35 min on a shared host, with V1's host
side and a cft-debris notice beside it; clean, build and check each
rc 0; 25 legs; 1,504 lines beginning `  ok `; 17 leg RESULT lines,
all PASS; no line beginning FAIL. Against P0's suite on `d7b5047`
(r5-check-d7b5047.log), every leg has the same ok count but
check_compare_run, 67 to 108; the build has the same 5 compiler
warnings; entry 45's six digests are present at the same counts.
check_ring_seams, the round-4 seam whose part 3 reads compare_run.py's
lines, has 34 ok, as at `d7b5047`.

### Stated limits, not closed

- **V3's open point:** the limit's absence is asserted in case 1, case
  5's two runs and case 7's kepler, not in every run without producer
  files (above).
- **A rerun member's schedule is not held pass for pass** against the
  producer's, though the files would allow it: not built, since the plan
  kept the member lines as before (P3).
- **A producer's passes are held by the trailer's counters,** a
  3-decimal mean and a maximum, not by the bits (above).
- **Not built:** sampling by a close approach to a moon, and a record
  naming non-default run settings (entry 50).
- **On the 10^6 ring,** a particle member costs five whole-file reads of
  its 175 MB problem (READ, the code; V3). Not run at that size.
- **`--particle-strata`'s default** is `--strata`'s, P3's choice.

### What this does not establish

- **A producer other than ias15_cft in the gate.** P2's files, and the
  P2 x P3 seam, are entries 56 and 58.
- **Anything on a GPU.**

## 58. Round 5's integration: PL's gate-cache keys and limit check; P2's first merge, whose GPU check failed, and its redo; the GPU checks and the suites

**2026-11-18 to 2026-11-20.** The lead's merges, one at a time, each on
its verifier's READY and a full make check on the desktop, read from its
log before main moved, and P1's, PL's and P2's also on a GPU check on
the RTX 5060 Ti in WSL `cft-gpu`:
- `4a8fe99`, P3 (entry 57), and `fda3b40`, P1 (entry 55), on 2026-11-18;
- `620c98b`, PL, on 2026-11-19;
- `61e45bc`, P2, on 2026-11-20, after a first merge, `7444157`, that never
  reached main.

Logan delegated round 5's merges and pushes on 2026-11-16 ("Yes,
delegate"). Pins as in entry 53.

### PL: the lead's own fixes, done by a parcel

PL made five commits on `fda3b40`, verified by VL:
- **`3012353`, the PROGRAMS key.** PROGRAMS was a wildcard over
  programs/out, expanded when make read the Makefile, so in a fresh clone
  the 10 legs that list the program images keyed on none of them (the
  lead's find, 2026-11-16). It now comes from programs/*.cfta, united
  with whatever programs/out holds. MEASURED through make (VL): with
  programs/out absent at both parses, a program changed in
  tools/gen_programs.py was SKIPPED by `fda3b40`'s check_program_engine
  line, of 10 inputs, and the new line, of 76, ran and FAILed it.
- **`5fa7534`, the limit held absent everywhere** (V3's open point, entry
  57). `limit_lines()` finds the line however it is indented or
  embedded, in limit_once and limit_absent. drive() records every run
  without `--producer-files`, and one closing verdict holds them all.
  Case 6 runs in every tier, 4 steps in `--quick`. Copies of compare_run
  that print the line where it must not appear each FAIL the gate by
  name: PL's M1 to M5, and VL's (a) to (f), five of which pass
  `fda3b40`'s gate. The gate gives 73 ok in `--quick` and 109 in full,
  against 70 and 108.
- **`cc3d87a`, check_ensemble's key.** Its program-engine cases load
  images from programs/out, and its key listed none. MEASURED (PL, at the
  cache): with one byte of predict-ens-fp64.cftp changed, `5fa7534`'s
  line SKIPPED, and the new line ran and FAILed the two adaptive
  program-engine cases. VL's, through make: a change to predict_ens in
  gen_programs.py, SKIPPED by the old line and FAILed by the new.
- **`e8a7d89`, the audit.** PL audited every leg's key statically, in
  all three tiers, and added the 6 inputs missing in 3 legs:
  - tools/check_equivalence.py, imported by check_records,
    check_ensemble and check_bodycount;
  - tools/make_ensemble.py and tools/make_nbody.py, run by check_ensemble
    and check_bodycount;
  - src/ias15_limits.h, read by check_bodycount.

  MEASURED for each kind of input, an import, a run script and a file
  read: after a one-byte change, `cc3d87a`'s line SKIPPED the leg, and
  the new line ran it.
- **`4fcd775`, the constants leg,** VL's find. VL traced a cold
  check-quick, logging every file opened, process started, directory
  listed and module imported. The leg reads
  third_party/rebound/src/integrator_ias15.c, which its key did not list.
  MEASURED (PL, then VL, also through make): with one digit of h[1]
  changed in the clone, `e8a7d89`'s line SKIPPED, and the new line ran
  and FAILed, "derived 0.05626256053692215, REBOUND literal
  0.05626256053693215"; restored, it passed, then skipped warm. No leg
  reads anything else under third_party/ that its key does not list
  (PL's search, audit and trace; VL's trace).

**VL** said READY on `e8a7d89` with one defect, the constants key, for
which PL went back (send-back 3), and READY on `4fcd775`. VL's cache end
to end, MEASURED: cold, all 25 legs ran and passed; warm, all 25 SKIPPED
by name; and a change to predict_ens that moves three images' bytes and
no results ran exactly the 11 legs that list the images, while the other
14 skipped by name.

**The merge, `620c98b`** (2026-11-19; the lead's ledger, 13:39): `git
merge --no-ff round5-pl`, its tree equal to `4fcd775`'s.
`fda3b40..620c98b` changes the Makefile and tools/check_compare_run.py
alone, 120 insertions and 41 deletions.

**Its GPU check** (lead-wsl/plgpu-620c98b-1119-1339; the RTX 5060 Ti, a
fresh HIP=1 build in WSL cft-gpu): `check_gpu.py --legs
check_compare_run`, run directly rather than through make, 13:39:30 to
16:12:17, its software half on the shared desktop: "check_gpu: 2 ok, 0
FAIL, in 152m46s on hip:0", RESULT: PASS. The leg: rc 0, RESULT: PASS,
73 ok, 0 FAIL; 1606 of 1606 processes on the device, and 1606 of 1606
records the software backend's. Beside it: 665 records, 22,347 steps and
435 rejected, 3.61 passes an attempt (at most 20), none capped. VL's
four checks all held:
- the leg's line;
- case 6's "4 steps" header, with its two ok lines;
- the closing verdict, "46 runs (10 with --plan, 1 with --rerun-engine,
  2 at binary128 or binary256, 3 in the FMA form, 3 through a copy of
  the driver)";
- 1606 processes, above the 774 of P1's run of this leg on `5b1cac6`. No
  hip:0 count of it exists at `fda3b40`, since P1's merge seam did not
  run it.

**Its suite** (int-check-620c98b.log): 13:38:35 to 16:55:59 on
2026-11-19, 3 h 17 min on a shared host, beside V2's check_wave runs,
the GPU check's software half and cft-fp256's M3S. Clean, build and
check each rc 0; 25 legs; 1,505 lines beginning `  ok `; 17 leg RESULT
lines, all PASS; no line beginning FAIL; 5 compiler warnings, as before.
Every leg has `fda3b40`'s counts but check_compare_run, 108 to 109 (my
count, and the lead's). With timings and temporary directories masked,
the log differs from `fda3b40`'s in:
- the commit line;
- the gate-cache lines of constants, check_records, check_ensemble and
  check_bodycount, which carry PL's inputs;
- case 6's header and verdicts;
- the new closing verdict, "none of the 48 runs of the driver without
  --producer-files printed a line of it".

Every line carrying a token of 16 or more hex digits, the digests among
them, is the same, 122 of 122 (the lead's). main fast-forwarded
`fda3b40..620c98b` and was pushed.

**PL's open notes** (the merge message):
- check_bodycount and check_gravity reuse the nbody<N>.txt problems they
  wrote earlier in build/, so their presence, not their bytes, steers the
  leg. With make_nbody.py broken, check-quick reran check_bodycount on the
  old problems and passed, then skipped it after build/problems was
  deleted (VL). make check after make clean catches it; check-quick and
  check-light do not.
- A program new in gen_programs.py is keyed from the next make on.
- An edit of a .cfta never reaches an image: make programs runs
  gen_programs.py first, which rewrites it.
- Control i copies the tree's own compare_run.py, so it ignores
  `--driver`.
- The closing verdict asks for a run of each of four kinds. Its counts,
  46 and 48, are printed, not asserted.
- Runs refused after reading a producer's files are not checked for the
  limit line. drive() passes the per-run checks only the last 600
  characters of stderr, so a second copy earlier on a producer run's
  stderr passes (VL's mutant g).
- `--quick`'s case 6 runs 4 steps, in which no member rejects an attempt
  after an accepted step; the full tier's 20 steps hold one.

### P2's first merge, `7444157`, which never reached main

**The merge** (2026-11-20, 00:59; the lead's ledger, 01:00). On
`round5-int` at `620c98b`, `-c merge.conflictstyle=diff3 merge --no-ff
--no-commit 9fa5956` stopped on the Makefile alone, in three conflict
blocks, and the 8 new files merged cleanly. resolve_tiers.py resolved the
three (below), and p2_hand_edit.py made the hand edit: check-gpu's
second recipe line, `tools/check_wave.py --build $(B) --device gpu
--quick --expect '$(GPU_EXPECT)'`, with `wave` among check-gpu's
prerequisites when HIP=1, and its header comment. The Makefile's sha256
was f82c107b, the expected file's. It was committed as `7444157`.

**Its GPU check** (lead-wsl/p2seam-7444157-1120-0100; the RTX 5060 Ti,
01:00:54 to 01:17:26):
- **`make HIP=1 check-gpu GPU_LEGS=check_records`: rc 2.** Its first
  line passed: "check_gpu: 2 ok, 0 FAIL, in 15m19s on hip:0". Its second,
  check_wave `--device gpu --quick`, died before running anything, at
  check_wave.py:1088:

  ```
  NotADirectoryError: [Errno 20] Not a directory: '/root/cft-rebound-r5-lead/tree/build/wave/1549'
  ```

  `make HIP=1 wave` builds the kernel as the file build/wave, and
  check_wave's default work directory had been build/wave/<pid> since
  P2's first commit, `3f7d67e` (READ by the lead). Every GPU run of
  check_wave before the merge had passed `--keep` outside build/ (entry
  56), so this was the first run of check-gpu's own line.
- **The lead's seam script, p1p2_seam.py: 15 ok, 8 FAIL.** All eight were
  the script's. It did not know the kernel record's `# device` line
  (P2's `46ddfed`), so wave against cpu found "1 difference(s)", the
  device line alone, and each of six wave controls found 2. A second
  cause hid behind the first: the script wanted the kernel's backend as
  `gpu`, and the kernel writes `hip:0:sm_120`. compare_run held the
  kernel's record IDENTICAL, and the P2 x P3 part passed.

**make check on `7444157` was stopped** at 01:21, in its sixth leg,
check_door (int-check-7444157.log). Send-back 7 would change the merged
tree, so that run could not be the bar.

### The response

- **Send-back 7** (entry 56): check_wave's default work directory is
  build/check_wave/<pid>. V2 said READY on `3d8ab5b`.
- **The seam script, fixed and verified twice.** At c9e39a38, the kernel
  record's device line is set aside by name and held by a verdict of its
  own: exactly one in the kernel's GPU record, naming `--expect`, and none
  in the port's records. The kernel's backend is wanted as
  hip:0:<arch>, the arch read from that line.
  - MEASURED by the lead: on the failed run's saved records, 16 ok and 0
    FAIL, where the original script FAILs; on the GPU, against the built
    `7444157` tree, 24 ok and 0 FAIL.
  - VL said READY through the script's own main(), its run() stubbed
    from the saved records: all nine planted cases FAIL by name, and
    unknown header lines still FAIL.
  - Three cases today's kernel cannot produce still passed. The lead
    applied VL's three hardenings (f4e475dd): the record's device line
    must equal the probe's `--device-info` line, `--dry-host` requires no
    device line, and the arch may end at ':' as well as ','. VL said
    READY again: the three FAIL by name, and the saved records still
    give 24 ok.
- **The redo.** `git reset --keep 620c98b` on round5-int, `7444157` kept
  in the reflog and never pushed, and the merge redone on `3d8ab5b` with
  the same kit: the same three blocks, the same resolution and hand
  edit, and the Makefile again f82c107b. P2's files are `3d8ab5b`'s and
  every other path is `620c98b`'s. The merged tree differs from
  `7444157`'s in tools/check_wave.py alone, +6 and -2. It was committed
  as `61e45bc`.

### The conflict, resolved by a script and checked against a Makefile built another way

MEASURED by the lead before either merge (git merge-tree, in both
orders): P2 changed each tier's target line, `check: ... wave-host ...`,
and PL's `4fcd775` changed the constants line just under it, so the
three tiers make three adjacent-line conflicts. With `e8a7d89` there
were none.
- **resolve_tiers.py** works on the diff3 blocks and resolves each line
  by line: where one side kept the base's line, it takes the other
  side's. It refuses, and writes nothing, unless:
  - there are exactly 3 blocks;
  - ours, base and theirs are the same length in each;
  - no line was changed by both sides;
  - every resolved line is a tier's target line or a constants line;
  - the result equals an expected file byte for byte.
- **The expected Makefile, by another road:** P2 merged with `e8a7d89`, a
  clean merge to tree 399adaf7, plus `4fcd775`'s three lines by
  substitution; with P2's hand edit, sha256 f82c107b.
- **MEASURED:** both orders, PL then P2 and P2 then PL, resolve equal to
  it. The refusals shown, each leaving the file unchanged: merge-style
  markers, two blocks, and a one-byte-wrong expected file. The lead's own
  first expected file was wrong by a doubled space, and the equality
  check caught it.
- **p2_hand_edit.py** lands each of its 4 edits once and refuses a second
  run. `make -n -p HIP=1` resolves `check-gpu: check-pins all programs
  wave`.

### `61e45bc`'s checks

**The GPU check** (lead-wsl/p2seam-61e45bc-1120-0213; the RTX 5060 Ti, a
fresh clone, 02:13:13 to 02:42:23):
- **`make HIP=1 check-gpu GPU_LEGS=check_records`: rc 0.**
  - check_records: 4 ok, 0 FAIL, in 15m21s, 4 of 4 processes on the
    device.
  - check_wave `--device gpu --quick`, the second line, its command
    without `--keep`: RESULT: PASS, 147 ok, 0 FAIL.
  - Its gate 4 on the GPU: 228,352 cases and 0 mismatches, the corpus
    digest ca7246a62d9533f6, the host build's.
  - Every record's one device line is the GPU's, and the count verdicts
    hold on every record.
- **p1p2_seam.py at f4e475dd: 24 ok, 0 FAIL.**
  - **P1 x P2, three ways:** one ring wave at fp64 in the FMA form, from
    a first step of 1e5 s, run by the port on the CPU, by the port on
    hip:0 and by the kernel. The kernel's backend is hip:0:sm_120, as
    wanted, and its one device line names the 5060 Ti and equals its
    `--device-info` line. hip against cpu and wave against cpu are
    IDENTICAL: 4 sample lines of 394 values, 64 body lines, the header
    lines and the trailer's five counters. Beside them: 12 steps, 1
    rejected, 3.769 mean passes and 6 at most, and all 64 bodies moved.
  - **The controls:** 11, one value planted in a copy of the hip:0
    record or the kernel's, each FAILing by name on its one planted
    value.
  - **compare_run** holds the hip:0 record and the kernel's IDENTICAL to
    this tree's CPU rerun, 2 particles replayed in each.
  - **P2 x P3:** the kernel's own attempts records and schedules for
    three waves, through compare_run `--producer-files`: members 0 and 1
    not rerun, and 6 particles from P2's files IDENTICAL. One dt of
    member 0's attempts record moved one ulp FAILs that member's two
    particles by name.

**The suite** (int-check-61e45bc.log): 02:11:37 to 05:03:15 on
2026-11-20, 2 h 52 min on a shared host, beside the GPU check for its
first 31 minutes. Clean, build and check each rc 0; 26 legs; 1,649
lines beginning `  ok `; 18 leg RESULT lines, all PASS; no line
beginning FAIL; 5 compiler warnings.
- The 25 legs it shares with `620c98b` have that suite's counts (my
  count, and the lead's).
- check_wave, the 26th, after check_reversal, has 144 ok and RESULT:
  PASS. Its count rule holds on A at 5 lines of 390 and of 54.
- Outside check_wave's leg, with timings and temporary directories
  masked, the log differs from `620c98b`'s only in the commit line and
  the four new wave-host build lines. The 122 lines carrying a token of
  16 or more hex digits are all the same (the lead's).

main fast-forwarded `620c98b..61e45bc` and was pushed (the lead's
ledger, 05:05).

### The outage

The API was unreachable from about 12:35 to 12:56 on 2026-11-19. P2's
session ended at 12:38, and VL's after its 12:12 entry. Nothing was lost
(MEASURED by the lead):
- P2's `f3ef701` was committed, its worktree clean, and its mutant
  runners finished alone, every one of its 31 full-form mutants logged;
- VL's item-8 script carried on alone, and restored the files it had
  edited.

Both were resumed by message, and nothing needed re-running.

### Stated limits, not closed

- **PL's open notes** (above), and P2's (entry 56).
- **The pin.** Everything in round 5 was measured at cft-fp256 `eae74e3`.
  cft-fp256's GPU round 2 landed during the round, ending at 9e8ff270,
  with 340cd287 and c3c39c87 docs-only past it (READ by the lead in its
  repository). The move waits on Logan's word, after the round closes:
  "Agreed, finish round 5 on eae74e3e first" (Logan, 2026-11-18).

### What this does not establish

- **check-gpu in full at the merges.** At `620c98b` it ran
  check_compare_run alone, and at `61e45bc` check_records and check_wave's
  line. The other legs last ran on a GPU at P1's commits (entry 55).
- **The RX 7600 at any merge.** Every merge's GPU check ran on the RTX
  5060 Ti.
- **A full suite on Linux** at any of these merges.
- **These records' own merge and suite,** which follow this entry.
