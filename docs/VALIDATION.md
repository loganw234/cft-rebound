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
predictor, caught by the compiler). What the gate does NOT exercise:
the rejected-step path (`dt_new/dt_done < 0.25`) and the `ratio > 20`
branch of `predict_next_step`, because neither occurred in 1,000 steps
of these two problems. They are written to mirror REBOUND and are
untested.

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
