# Working on cft-rebound

Short, and only things that are non-obvious AND have already cost hours.
`README.md` says what this is; this says what will bite you.

## Building on Windows

```bash
export PATH="/c/msys64/mingw64/bin:$PATH"
make OS=Windows_NT TMP="$TEMP" TEMP="$TEMP" \
     PYTHON="/c/Users/logan/AppData/Local/Programs/Miniconda3/python.exe" <target>
```

**The PATH line is not optional.** Bare `cc` resolves to
`/c/msys64/mingw32/bin/cc` — the **32-bit** toolchain — because mingw32
is on the user PATH and mingw64 is not. The failure does not look like a
PATH problem: the link dies with

```
build/ias15_cft_lib.o: file format not recognized
```

on an object `file` reports as a perfectly good `Intel amd64 COFF`.
`cc -dumpmachine` settles it in one call. A fresh shell does **not**
inherit the PATH from an earlier successful build in the same session.

## The three suite tiers

| target | what it means | cost |
|---|---|---|
| `make check` | every gate at every format, **never skips** | 2 h 36 min |
| `make check-quick` | cheaper arguments + the gate cache | 42 min cold, 13 s warm |
| `make check-light` | cheaper arguments + cache + **scaled step counts** | 31 min cold |

19 legs each. The times are from clean, build included, on 2026-11-06.

`make check` is the one whose meaning must not change, so it always
executes. The cache (`tools/gate_cache.py`) keys on the **artifacts a leg
executes**, never on source files, and a failure **retires** the stamp.
Every skip prints itself by name; `make gate-cache-report` shows what is
held.

`--light` does **not** scale everything, and three families must stay
unscaled: the integrate-to-a-time cases (they name a time, not a count),
the step-sequence cases (`dt_sequence` steps one at a time and the scaler
leaves a 1 alone), and the two merge cases, which opt out by name through
`(reb_simulation_steps)(...)` — the parentheses suppress the macro.
Scaling a merge case makes it a run where nothing happens.

`CFT_REBOUND_ARTIFACT` is part of every cache key, so a pass on the
software backend can never satisfy a run against a card.

## Three ways this suite has lied

- **Bit-identity is not validity.** Identical results prove the two
  backends agree, not that the problem was the one you meant. Massless
  test particles made a difference exactly zero and every case passed.
  Print a physics counter — "a merge removed one: 3 particles became 2" —
  beside the comparison, and assert it.
- **A suite is only as prebuilt as its least-compiled component.** A
  merge landed mid-run once and the suite read a new Python checker
  against an old binary, producing three failures that were not real.
- **Read the log, not the exit code.** The first green `make check-light`
  printed `400 steps` over a leg that ran 40; the only contrary evidence
  was the leg finishing in 159 s instead of 585.

## On a card, tile count is a cost

Point `CFT_REBOUND_ARTIFACT` at a **one-tile** image, or open a quad with
`CFT_XRT_TILES=1`, which cost 1-10% against the one-tile image on the
small problems (measured 2026-11-02). libcft partitions
every `cft_run` across every compute unit a device presents, so one
library call becomes a kernel launch and a staging round *per tile*, and
this integrator issues thousands of small calls per step. Measured
2026-09-14: one tile beat four beat six in all ten card rows, every body
count and both formats, and a six-tile image with **1.39x** the
aggregate arithmetic integrated **10% slower** than the shipped quad.

The same measurement priced the arithmetic at about **a fifth of card
wall clock** - an 11% faster engine moved the integration 2%. So a
faster or wider bitstream is not the lever; the engine's vectors being
plain `calloc` rather than device-resident `cft_alloc` is. docs/
HARDWARE.md ranks the rest, docs/BITSTREAM.md has the bitstream flow,
and docs/VALIDATION.md entries 34-36 are the numbers.

Two runs of one configuration differ by up to 6% on the small problems
and 2% on the n-body rows, and three span up to 8% (entry 45), so one run
per row cannot rank two routes by less than that.

A bench has three clocks. `workload.csv`'s `seconds` is the parent's,
around the process. `decompose.csv`'s `wall_seconds` is the program's own,
0.6 to 1.2 s less, and its `seconds` is CPU `clock()`, under the same
name. Compare one with itself: round 3's first before-and-after mixed two,
and the headline was 2 points out.

And a row is not a step. `energy_all()` sums with one width-one call a
pair, twice a row: 96% of the calls at 256 bodies and half the nbody512
fp64 row's time. Round 3 first read that as a cliff in the step.

## Every route is opt-in, and off is not the old code

`auto` takes none of the door's four routes. `indexed`, `pair_gather` and
`reduce_seg` measured slower on a tile with identical bits, `scalar`
neutral, so each stays off, `by=opt_in`, even where a card publishes its
bit. A route leaves the list when a card measures it winning. Two things
follow, and both have already bitten:

- **An unforced run sees today's route only.** A library half silently on
  today's route passes every unforced leg. So each route has a forced-on
  `check_dropin --light` in all three tiers, and a new route needs one.
- **Off is not the code from before the round.** Gravity with both
  switches off is still P1's program form, with the host gathering; the
  22-call gravity survives only as the self-test's reference. Route
  against off compares two new things. The "before" is entry 45's.

## Launching a long suite from an agent

Jobs the Bash tool starts run at BelowNormal and starve on this host, so
start one from PowerShell, through Git's bash:

```
Start-Process -FilePath "C:\Program Files\Git\bin\bash.exe" -ArgumentList "<script.sh>" -WindowStyle Hidden -PassThru
```

- The script's path is the only argument. PowerShell 5.1 re-splits `-c
  "..."` on embedded quotes, and drops the make target.
- Inside the script, spell TMP and TEMP out as the Windows directory.
  Started this way `TEMP=/tmp`, which msys64's make hands to native
  children as `C:\msys64\tmp`.
- `make -j` races on `libcft.a` in a fresh tree. Two `gate_real` runs in
  one directory overwrite each other's `gate_real.bin`.

## The trap in the port itself

`dpcast()` ([src/ias15_cft.c](src/ias15_cft.c)) slices **one flat buffer**
into all seven coefficient levels at the **current** `3N`, every step. A
particle removed below the high-water mark therefore **aliases** levels
by `old_3N - new_3N`. REBOUND does this too, and at `accurate = 0` this
port reproduces it exactly — defect included — because bit-identity with
REBOUND is the claim. `accurate = 1` shifts the wide state instead and is
deliberately **not** bit-identical.

Agents: work in a worktree or clone, in the background, and **never push
or merge** — the integrator does both.
