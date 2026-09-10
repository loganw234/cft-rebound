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

    src/ias15_cft.c          the port: every floating-point operation is a cft.h call
    src/ias15_constants.h    GENERATED: the Gauss-Radau constants at every format
    src/hexfloat.h           exact hex-float text for binary64, libc-independent
    ref/ias15_ref.c          REBOUND's own IAS15 on the same problems, plain double
    ref/whfast512_stub.c     why REBOUND's WHFast512 is not built here (MinGW)
    tools/gen_constants.py   derives and CHECKS the constants (mpmath, 130 digits)
    tools/gen_programs.py    the predictor and corrector as sequencer programs
    tools/make_problems.py   the two problems as exact binary64 bit patterns
    tools/check_equivalence.py     the gate: binary64 port == REBOUND, bit for bit
    tools/check_program_engine.py  the gate: programs == host loop, bit for bit
    tools/oracle.py          scores a record from its exact bits (mpmath)
    tools/sweep.py           the precision-versus-error measurement
    tools/tabulate.py        its tables and figure
    programs/*.cfta          GENERATED: 24 orbit-sequencer programs
    data/problems/           kepler.txt, outer.txt
    data/ias15_constants.json
    results/                 the sweep's per-sample tables and summary
    docs/VALIDATION.md       what ran and what it said, failures included
    docs/HARDWARE.md         the resident design, its break-even, what is unverified
    docs/INTEGRATORS.md      which other REBOUND integrators port well, ranked

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
- **Precision.** The tables in docs/VALIDATION.md. The wider formats
  lower the round-off floor by exactly the format's precision, but
  IAS15 collects that only where its own truncation error is already
  below it, which costs a smaller step and, at binary256, about seven
  times the corrector passes per step.
- **The programs.** The predictor and corrector run as 24
  orbit-sequencer programs with the state in the per-lane scratch
  block, bit-identical to the host loop at every format.
- **The hardware.** docs/HARDWARE.md: the design keeps 47 values per
  coordinate resident; the projection is that a single system on the
  tile is 2-3x the software backend and an ensemble of a thousand is
  100-250x, that at binary64 the tile is three times slower than one
  CPU core running REBOUND itself, and that the scratch block - where
  the state lives - is staged every run under today's contract, which
  is the first thing to change. None of that was run on a card.
