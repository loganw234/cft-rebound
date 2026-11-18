#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of tools/compare_run.py, the driver that holds a run's records
# against this tree (ROUND4.md, P4).
#
# Today every record it is pointed at is this tree's own, so every case
# here compares the CPU with the CPU - METHOD.md's canonical failure, a
# gate driving a stand-in, unless the driver is shown to FAIL. So the
# controls are permanent and each must bite, by name: a driver that cannot
# fail is the failure this parcel exists to prevent.
#
#   python tools/check_compare_run.py [--build build] [--quick]
#
# THE WAVES. A small ring, hand-made here because tools/make_ensemble.py
# refuses `active` lines: E=9 systems of a primary, a moon and three
# massless test particles, `active 2`. Member k's particles start at radii
# 0.5 + k/16 + j/64, so the member index is radius, as in P2's generator.
# Orbits are near circular and inclined about 1.8 degrees, alternately up
# and down, from Pythagorean triples. Every coordinate is a chain of
# correctly rounded binary64 operations - divides, multiplies and one square
# root, and no other library function: a particle's x takes two roundings,
# its y, z and vx four, and its vy and vz six (V4 counted them) - so the
# file is the same on every IEEE host. A first step of 0.3 is rejected by
# the innermost members and accepted by the rest, so the sample's "every
# member with a rejected step" is exercised, and its strata are a proper
# subset.
# Member 8, the last, holds the other kind of rejection: its first particle
# is retrograde in the moon's own plane, 1/128 outside its radius, and meets
# it head on near t = 1, after the first attempt and a few more steps have
# been accepted; so it rejects attempts only after an accepted step. It was
# appended for that (2026-11-12), and the ensemble's members run apart, so
# members 0 to 7 have the bits they had when there were eight.
#
# CASES (all of them in --quick, at binary64, unless marked):
#   1. The driver on the ensemble's record, with --format-diff fp128: every
#      sampled member IDENTICAL, and asserted beside it from this gate's own
#      reading of the record - each member's steps, rejections and passes;
#      members that reject their first attempt, one that rejects only after
#      an accepted step, and members that reject nothing (each rejecting
#      member run alone for one step tells the first two apart); a sample
#      that holds every rejecting member and one member of every stratum,
#      and is not every member; each test particle's closest approach to
#      each active body, recomputed here in mpmath from tools/oracle.py's
#      parse; and each format difference nonzero and below 10^4 units of
#      the narrower format's round-off, with its scope beside its number in
#      the lead's ruling's words, and the note its member's rejections call
#      for. What each note says is held by a replay from the member's own
#      dt file at its own format: a member that rejected nothing is the
#      record's bit for bit, so its nonzero difference is another format
#      really run; one whose first attempt was rejected differs in state;
#      and member 8's state is the record's in every value but dt_next,
#      which differs only on a sample line printed just before a rejected
#      attempt.
#   2. The particles, through P3's replay (the lead's P3 x P4 seam): each
#      sampled member's two particles, its least closest approach at the
#      samples (the one this gate's own mpmath distances name) and one by
#      the seed, each replayed alone and IDENTICAL to the record; and beside
#      the sample, its limit, exactly: it trusts the trailers it is testing.
#   3. The sample: --plan for seeds 1 to 10, each with every rejecting
#      member and one member a stratum; not all ten alike; seed 1's is case
#      1's.
#   4. The same problem's records at binary128 and in the FMA form (short
#      runs): each held at its own format and arithmetic, PASS.
#   5. One system's records: a --member 3 run of the waves, and kepler run
#      as a single problem (every body active): each PASS, IDENTICAL.
#   (full only) Case 1 at binary128 and binary256, the binary128 one with
#      --format-diff fp256; and 6, the FMA form on the program engine, held
#      as recorded and again with --rerun-engine loop.
#   7. A producer's own attempts records and schedules (ROUND5.md, P3), at
#      every format case 1 runs. The producer is ias15_cft, each member run
#      alone with --attempts-out and --record-iters, under the names the
#      driver states (PRODUCER_FILES, imported here, so a driver looking
#      anywhere else finds nothing); each run's lines and counters are
#      asserted to be the record's member's, so the files are the record's
#      run's. Their directory's name holds " and ", which the driver's lines
#      put between the two paths, so no reading of those lines here may
#      split them there (V3, 2026-11-16): each is matched whole, against the
#      paths this gate wrote. The driver on case 1's record with
#      --producer-files and --particle-strata 9 samples every member's
#      particles, the members its member sample leaves NOT RERUN: PASS,
#      every pair held to its member's run by its key, and the limit of the
#      producer's files printed beside the plan, word for word as
#      PRODUCER_LIMIT here states it (V3's D1: only the driver's source
#      said it); and with an empty directory instead, each member falls
#      back to today's rerun, naming both its files. The limit is held as
#      lines, not substrings (two adjacent copies share a newline, and a
#      substring count saw them as one: V3, 2026-11-18): one line, after the
#      plan's, and no other line in any wording, in every run of this case
#      that reads a producer's files, twelve at binary64; and no line of it
#      in any wording where no files are read: case 1 at each format, case
#      5's two records, and kepler here.
#      GATE 1, the same particles: each particle, picked alike in the two
#      runs, replayed from the producer's files is IDENTICAL to the same
#      particle replayed from a rerun's - every sample line and the
#      trailer's counters, flags and calls, read from the replays' own
#      records, which the driver keeps and which each name the attempts
#      record they ran from - with the members whose first attempt was
#      rejected, the one rejection that leaves state, among them.
#      GATE 2, beside each particle: its replay's counters are its member's
#      trailer, read here, and so is each not-rerun member's trailer line;
#      and its closest approaches agree with this gate's mpmath.
#      GATE 1b, at binary64: gate 1's producer is this tree, so a driver
#      reading a rerun's schedule in place of the producer's passes it
#      (measured, with such a mutant of the driver). So a
#      producer whose exit test is not this tree's: two fewer passes on the
#      first attempt of a member that is not rerun, which moves its bits,
#      one fewer on a member that is, which moves no bit at the samples and
#      only its trailer's counters (which ones, the gate derives and
#      prints: V3 found the max moving with the mean where this said "only
#      the mean"), each member's lines that schedule's replay. From the
#      producer's files the first member holds - the limit, measured, as the
#      driver now prints it: without a rerun the producer's passes are
#      taken as given - and the second FAILs on its rerun, its particles
#      IDENTICAL; from a rerun's, all four of their particles FAIL. So the
#      same-particles gate cannot pass against itself. And the waves in
#      their primary's zonal field (`zonal 0 R J2
#      J4`): PASS, every particle replayed with the line. And half a pair,
#      one member's schedule taken away: that member falls back to a rerun,
#      naming the schedule alone, PASS. One system's records: a --member k
#      record takes member k's pair, PASS; a single problem with no system
#      lines and an active line (the waves' w0), run with no --member, keys
#      its files member=-1 and names them member0's, PASS; kepler, every
#      body active, says the files are not read, PASS. And case 1's
#      record rewritten as another producer's (program=wave backend=hip:0
#      engine=wave), held with --rerun-engine loop: its files keyed at the
#      engine it is rerun on, PASS; with no --rerun-engine, REFUSED, naming
#      the engine. And the FMA form's program-engine record, rerun on the
#      loop engine: its files keyed at engine=program, the engine it ran,
#      and its particles replayed there, PASS - 4 steps in --quick, 20 in
#      full, so that a driver keying the files to --rerun-engine's engine
#      (V3's mutant Mj, which only the full gate's case saw) is refused by
#      name in every tier.
#
# NEGATIVE CONTROLS, each asserted to FAIL or be REFUSED by name:
#   a. one byte of one sampled member's last sample line (an x of a test
#      particle, its last hex digit): FAIL, naming that member and line and
#      the column; every other sampled member IDENTICAL; and a2, the same
#      in case 5's one-system record. Each of those moves a binary64 bit,
#      and a driver comparing binary64 values passed them all (V4); so
#      a3, the sign of a zero flipped (a sample-0 dt_last made -0), and a4,
#      one hex digit past binary64 changed in case 4's binary128 record
#      (its binary64 value, asserted here, unmoved): each FAIL, naming the
#      member, the line and the column; and a5, every column the record
#      names, read from its own `# columns:` line - each field before the
#      bodies, each field of the last body (so a comparison stopping short
#      of the end fails) and the exact-time pair - one planted change each,
#      one column a sampled member and so two records for twelve columns:
#      each FAIL, naming its column (a driver skipping the pair passed
#      every control before this one, V4);
#   b. one byte of a sampled member's trailer (mean_pc_iterations): FAIL,
#      naming the member and the counter;
#   c. every line of one member that is not sampled deleted: FAIL, naming
#      it as missing; the sampled ones IDENTICAL;
#   d. the record cut short before its last sample: FAIL, every member
#      incomplete and no trailer;
#   e. a record of another format (case 4's binary128) held as binary64,
#      and of another arithmetic (case 4's FMA form) held as REBOUND's:
#      REFUSED, naming the format and the arithmetic;
#   f. REBOUND's own record (ias15_ref), which shares the first line:
#      REFUSED, naming what it does not carry; a problem file given as the
#      record: REFUSED, not a record;
#   g. a different problem file under the record's name: a test particle's
#      mass 0 made 2^-40, the active line dropped (neither changes a bit:
#      ROUND4.md P0's inert claim), one initial coordinate moved one ulp,
#      and the name changed. Each REFUSED, naming what differs;
#   h. the record's `# active` line with one more pair than this tree's
#      rows: REFUSED, naming the pair count;
#   i. P3's tools in the tree and not callable: a copy of the driver beside
#      an empty stand-in for each of P3's two files in turn says P3's replay
#      cannot be called, ends RESULT: FAIL and exits 1, though its member is
#      IDENTICAL; beside neither, it passes, saying no particle is sampled;
#   p. the particle path is held to the record, not to the rerun: one hex
#      digit of a sampled particle's x in a sampled member's line, and of
#      that member's primary's x: the member FAILs, and so does each
#      sampled particle's replay the change reaches, naming the line and
#      the column (P3 x P4, the lead's); and p2, every column the particle
#      path compares, read from the record's own `# columns:` line - each
#      field before the bodies, every field of each massive body and of the
#      particle, and the exact-time pair - one planted change a member a
#      record, on the member's sample lines after sample 0 in turn, and
#      three trailer counters: each FAILs that member's least-closest
#      particle's replay, naming its line and column (VL: control p alone
#      held two columns, and plants on the last line alone held one line).
#   Case 7's, each on a copy of the producer's files, deleted after it runs:
#   q1. one attempt's dt in an attempts record (its last hex digit), a
#      member a line: the first attempt of a member whose first attempt was
#      rejected, and of one whose was accepted, a middle accepted attempt,
#      and the next line; and the other column, a rejected first attempt's
#      outcome made an acceptance. At every format: each FAILs both its
#      member's particles by name; every other particle, and every member
#      rerun, IDENTICAL.
#   q2. one pass count in a schedule, a member a line: two fewer on a first
#      attempt and one more on a middle one, in members not rerun; one
#      fewer on a first attempt in a member rerun; and two fewer on the
#      rejected first attempt of a member whose first attempt was rejected
#      (V3, 2026-11-16: with no plant on a rejected attempt, a driver that
#      read a rerun's schedule for rejecting members alone passed, its
#      mutant Mf2). At every format: each FAILs both its member's particles
#      by name, by a line and column or a trailer counter, while that
#      member's own rerun is IDENTICAL; and the gate prints, for each, what
#      the plant moves - the bits, or which of the trailer's counters - as
#      measured, not asserted: what a plant moves at a given format is a
#      description (V3, 2026-11-18), and only that it moves something is
#      held. It reads the replays gate 1 kept, and a replay that is not
#      there FAILs by name.
#   q3. a schedule, and an attempts record, from another wave of the same
#      ring: REFUSED by the run key, naming member and state, and nothing
#      compared;
#   q4. a schedule from the zonal waves under another J2: REFUSED, naming
#      zonal and no other field;
#   q5. a schedule recorded at binary128: REFUSED, naming format among the
#      fields;
#   q6. a schedule recorded under `active 3`: REFUSED, naming active and no
#      other field;
#   q7. an attempts record given as the schedule: REFUSED by its form.
#
# Every case runs a fixed seed and fixed bits, so the verdict is the same
# on every run (tools/gate_cache.py). None of its step counts is scaled by
# --light: check-light runs it as check-quick does, like every Python leg.

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import time

import mpmath

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_ensemble import parse as parse_samples          # noqa: E402
from check_equivalence import hex_to_fraction               # noqa: E402
from oracle import parse as oracle_parse, F2m               # noqa: E402

root = os.path.dirname(here)
mpmath.mp.dps = 40

E = 9
NTP = 3
DT0 = "0.3"
PHASES = [(3, 4, 5), (5, 12, 13), (8, 15, 17), (7, 24, 25), (20, 21, 29), (12, 35, 37)]
INC = (4095, 128, 4097)        # cos and sin of about 1.8 degrees, exactly rational
# The last member's first particle meets the moon head on: retrograde, in
# the moon's own orbital plane, 1/128 outside its radius, starting 112.6
# degrees ahead of it (cos, sin = -5/13, 12/13), so they meet near t = 1.
FLYBY, FLYBY_R, FLYBY_AT = E - 1, 1.0 + 1.0 / 128, (-5, 12, 13)
PARTICLES = "particles: not sampled (no particle replay in this tree)"
BROKEN = "particles: FAIL: P3's replay is in this tree and cannot be called: "
PICK = re.compile(r"^  member (\d+): particle (\S+) \(body (\d+), (its least closest approach at the samples|by the "
                  r"seed)\): replayed alone .*?held to the record: (IDENTICAL|FAIL): (.*)$", re.M)
# The sample's limit, as V4 measured it: printed beside the sample.
SAMPLE_LIMIT = ("sample: it trusts the trailers it is testing: a member whose trailer under-reports its rejections "
                "leaves the sample, and an unsampled member's lines are never rerun")
# What --format-diff must say beside its number (the lead's ruling,
# .ledger/urgent/lead-gate4-ruling.md), and the note a rejecting member's
# difference must carry, by whether its first attempt was rejected. Stated
# here as this gate's expectation, not read from the driver.
SCOPE = "the difference shows the narrower format's round-off on these steps; a narrowing both runs share passes it"
FIRST_NOTE = ("; this member's first attempt was rejected, before any accepted step, which leaves state the replay at %s "
              "does not take, so the difference includes that")
LATER_NOTE = ("; this member's %d rejected attempt%s all followed an accepted step, which leaves no state bit, so the "
              "replay at %s follows its path")
FORMATS = ("fp64", "fp128", "fp256")
# 2^-p for a p-bit significand: binary64, binary128 and libcft's binary256
# (237 bits, as ias15_cft's own pc_tol_shift of 237 - 53 says)
UNIT_ROUNDOFF = {"fp64": 2.0 ** -53, "fp128": 2.0 ** -113, "fp256": 2.0 ** -237}


def body(name, m, r, ph, up, gm):
    c, s, n = ph
    cp, sp = c / n, s / n
    ci, si = INC[0] / INC[2], INC[1] / INC[2]
    v = math.sqrt(gm / r)
    vals = (m, r * cp, r * sp * ci, up * r * sp * si, -v * sp, v * cp * ci, up * v * cp * si)
    return "body %s %s" % (name, " ".join(float.hex(x) for x in vals))


def flyby(name, r, at):
    """A massless body on a circular orbit in the moon's orbital plane,
    going the other way round. The plane's axes are the moon's starting
    direction e_a = (0, ci, si) and its velocity's, e_b = (-1, 0, 0); the
    body starts at angle th from e_a toward e_b, cos and sin = c/n and s/n,
    at r (cos th e_a + sin th e_b), with velocity v (sin th e_a - cos th e_b)."""
    c, s, n = at
    ct, st = c / n, s / n
    ci, si = INC[0] / INC[2], INC[1] / INC[2]
    v = math.sqrt(1.0 / r)
    vals = (0.0, -r * st, r * ct * ci, r * ct * si, v * ct, v * st * ci, v * st * si)
    return "body %s %s" % (name, " ".join(float.hex(x) for x in vals))


def write_waves(path, name="cr_waves", tp_mass=None, ulp=None, active=True, zonal=None):
    """The waves. tp_mass=(k, j, m) gives system k's particle j mass m;
    ulp=(k, j) moves system k's particle j's x by one ulp; active=n other
    than True writes `active n`; zonal, the text of a zonal line after the
    word (case 7's waves in their primary's own field)."""
    lines = ["# tools/check_compare_run.py: a small ring of waves, member index = radius",
             "name " + name, "G " + float.hex(1.0), "N %d" % (2 + NTP), "E %d" % E]
    for k in range(E):
        lines.append("system w%d" % k)
        lines.append("body primary %s %s" % (float.hex(1.0), " ".join(["0x0.0p+0"] * 6)))
        lines.append(body("moon", 2.0 ** -12, 1.0, (0, 1, 1), 1, 1.0 + 2.0 ** -12))
        for j in range(NTP):
            m = tp_mass[2] if tp_mass and tp_mass[:2] == (k, j) else 0.0
            if (k, j) == (FLYBY, 0):
                lines.append(flyby("tp0", FLYBY_R, FLYBY_AT))
                continue
            b = body("tp%d" % j, m, 0.5 + k / 16 + j / 64, PHASES[(k + j) % len(PHASES)], 1 - 2 * (j % 2), 1.0)
            if ulp and ulp == (k, j):
                t = b.split()
                t[3] = float.hex(math.nextafter(float.fromhex(t[3]), math.inf))
                b = " ".join(t)
            lines.append(b)
    if active:
        lines.append("active %d" % (2 if active is True else active))
    if zonal:
        lines.append("zonal " + zonal)
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return path


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=root)
    return p.returncode, p.stdout, p.stderr


def columns(text):
    """This gate's own reading of a record's `# columns:` line - not the
    driver's: the fields before the first body, and each body's names."""
    line = next(l for l in text.splitlines() if l.startswith("# columns:"))
    head, _, rest = line.split(":", 1)[1].partition(" then per body ")
    per = rest.split(";")[0].split()
    lead = head.split()
    # "system s" is two fields and so is "sample" (the word and its index)
    nlead = len(lead) + 1
    return nlead, per


def column_plan(text, N, names):
    """Every column this gate plants a change in, read from the record's own
    `# columns:` line: each field before the bodies after the sample's index
    and step, each field of the LAST body's block (so a comparison that
    stops short of the end fails too), the exact-time pair, and one field of
    every other body, the fields in turn (V4's third look, 2026-11-12: a
    driver that skipped body 0 passed, and body 0 is the primary, whose last
    bits the zonal reaction moves). Each comes as (the name the driver must
    give it, its token in an ensemble line). `names` are the bodies' names,
    from the record's `# body` lines."""
    line = next(l for l in text.splitlines() if l.startswith("# columns:"))
    head, _, rest = line.split(":", 1)[1].partition(" then per body ")
    per, _, pair = rest.partition(" ; then exact-time pair ")
    head, per, pair = head.split(), per.split(), pair.split()
    pre = 0
    if head[:1] == ["system"]:
        head, pre = head[2:], 2
    first_body = pre + 1 + len(head)
    bar = first_body + N * len(per)
    plan = [(f, pre + 1 + head.index(f)) for f in head[2:]]
    plan += [("%s of body %d (%s)" % (f, N - 1, names[N - 1]), first_body + (N - 1) * len(per) + per.index(f)) for f in per]
    plan += [("the exact-time pair's %s" % f, bar + 1 + pair.index(f)) for f in pair]
    plan += [("%s of body %d (%s)" % (per[b % len(per)], b, names[b]), first_body + b * len(per) + b % len(per))
             for b in range(N - 1)]
    return plan


def particle_plan(text, N, names, nact, j):
    """Every column the particle path compares for particle j of a member
    (tools/compare_run.py, particle_compare), read from the record's own
    `# columns:` line: each field before the bodies after the sample's index
    and step, every field of each massive body and of the particle, and the
    exact-time pair. Each as (the name the driver gives it, its token in an
    ensemble line). VL, 2026-11-12: control p planted only two of these, and
    a particle_compare that checked only those two passed."""
    line = next(l for l in text.splitlines() if l.startswith("# columns:"))
    head, _, rest = line.split(":", 1)[1].partition(" then per body ")
    per, _, pair = rest.partition(" ; then exact-time pair ")
    head, per, pair = head.split(), per.split(), pair.split()
    pre = 0
    if head[:1] == ["system"]:
        head, pre = head[2:], 2
    first_body = pre + 1 + len(head)
    bar = first_body + N * len(per)
    plan = [(f, pre + 1 + head.index(f)) for f in head[2:]]
    for b in list(range(nact)) + [j]:
        plan += [("%s of body %d (%s)" % (f, b, names[b]), first_body + b * len(per) + per.index(f)) for f in per]
    plan += [("the exact-time pair's %s" % f, bar + 1 + pair.index(f)) for f in pair]
    return plan


def flip_last_digit(tok):
    """A hex float with its mantissa's last hex digit changed (XOR 1)."""
    mant = tok.split("p")[0]
    return mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + tok[len(mant):]


def scalar_tok(text, name):
    """Where a field before the bodies is in a split sample line, read from
    the record's own `# columns:` line: "system s" is two tokens, the word
    sample one more, and each name after it one."""
    line = next(l for l in text.splitlines() if l.startswith("# columns:"))
    head = line.split(":", 1)[1].split(" then per body ")[0].split()
    pre = 0
    if head[:1] == ["system"]:
        head, pre = head[2:], 2
    return pre + 1 + head.index(name)


def value_index(text, name):
    """Where a field is in check_ensemble.parse()'s list of a sample's values
    (every field after the sample's index and step), read from the record's
    own `# columns:` line."""
    line = next(l for l in text.splitlines() if l.startswith("# columns:"))
    head = line.split(":", 1)[1].split(" then per body ")[0].split()
    if head[:1] == ["system"]:
        head = head[2:]
    return head[2:].index(name)


# ------------------------------------------------------------------ case 7
# The producer's own files (ROUND5.md, P3; tools/compare_run.py, THE
# PRODUCER'S OWN FILES). The names a producer writes are stated once, in the
# driver, and imported here to write ias15_cft's under them: so a driver
# that looked anywhere else would find none, fall back, and fail gate 1.
from compare_run import PRODUCER_FILES                       # noqa: E402

# a member's counters, as the record's trailer carries them
TRAILER = ("steps_done", "steps_rejected", "iterations_max_exceeded", "mean_pc_iterations", "max_pc_iterations")
SAID = ("steps_done", "steps_rejected", "mean_pc_iterations", "max_pc_iterations", "iterations_max_exceeded")
PF_LINE = re.compile(r"^  member (\d+): its lines held by (a rerun here|its particles' replays alone, NOT RERUN): (.*)$",
                     re.M)
PF_PART = re.compile(r"^  member (\d+): particle (\S+) \(body (\d+), (its least closest approach at the samples|by "
                     r"the seed)\): (.*)$", re.M)
PF_VERDICT = re.compile(r"held to the record: (IDENTICAL|FAIL): (.*?); its replay's counters: steps (\S+), rejected "
                        r"(\S+), mean passes (\S+), max passes (\S+), cap hits (\S+); its closest approach AT THE "
                        r"SAMPLES: (.*)$")
# The two sources a particle line can name. Paths follow the producer's
# phrase, and are never split out of it (a path may hold " and ", the word
# the driver puts between them: V3, 2026-11-16): the source is the phrase
# that comes first, and a pair is matched whole against the paths written.
OWN, RERUNS = "the producer's own attempts record and schedule, ", "this tree's rerun's attempts record and schedule"
# A member that is not rerun, as the driver must report it: held by its
# particles' replays only, and within the limit printed beside the plan
# (V3's D1: "its lines are held by the record alone, and by its particles'
# replays" claimed more).
PF_ALONE = re.compile(r"^member (\d+) \(\w+\): NOT RERUN: held only by its particles' replays from the producer's own "
                      r"files - its massive bodies, its sampled particles, each sample's head fields and exact-time "
                      r"pair, and its trailer's counters - and its other particles by nothing, within the limit printed "
                      r"above; its trailer: steps (\S+), rejected (\S+), mean passes (\S+), max passes (\S+), cap hits "
                      r"(\S+)$", re.M)
# What the producer's files cannot hold, which the driver must print beside
# its plan, once, and never without --producer-files. Stated here as this
# gate's expectation, not read from the driver, as SAMPLE_LIMIT and SCOPE are
# (V3's D1, 2026-11-16). Its last clause: a rerun holds a setting only
# through what it moves (V3, 2026-11-18: pc_tol_shift 2 and 8 and max_iter
# 12 IDENTICAL on a rerun, epsilon 1e-8 and max_iter 6 FAILing).
PRODUCER_LIMIT = ("particles: the limit of the producer's own files: each replay takes the producer's dts, outcomes "
                  "and passes as given, so only a member that is rerun holds its decisions; a NOT RERUN member's "
                  "particles that are not sampled are held by nothing; and the run key holds the files to the "
                  "arguments the record names, not to those the producer ran with, so its epsilon, pc_tol_shift and "
                  "max_iter, and the settings no record carries, are held only in a member that is rerun, and there "
                  "only through what they move in its lines")
PF_PLAN = "\nparticles: from the producer's own files in "
# the opening words of the limit's line, by which a copy of it in any wording
# is found where none may be
LIMIT_OPENS = "particles: the limit of the producer's own files"
PF_REFUSED = re.compile(r"^RESULT: REFUSED: member (\d+)'s (attempts record|corrector schedule) (.+?) was recorded for "
                        r"another run, and a schedule replays only its own run: its run key differs from member \d+'s "
                        r"in (.+?) \(", re.M)
CLOSE = re.compile(r"(\S+) (\S+) at sample (\d+) \(step \d+\)")


def pf(d, k):
    """Member k's pair in a producer's directory, by the driver's names."""
    return tuple(os.path.join(d, name % k) for _, name in PRODUCER_FILES)


def same_path(a, b):
    return a is not None and b is not None and \
        os.path.normcase(os.path.normpath(a)) == os.path.normcase(os.path.normpath(b))


def pf_members(out):
    """The driver's line for each member it samples particles in: how its
    lines are held, and whether it names a producer's pair (which one is
    names_pair()'s question, asked with the paths this gate wrote)."""
    res = {}
    for m in PF_LINE.finditer(out):
        res[int(m.group(1))] = {"rerun": m.group(2) == "a rerun here", "text": m.group(3),
                                "producer": "; its particles from the producer's own " in m.group(3)}
    return res


def names_pair(entry, pair):
    """Whether a member's plan line or a particle's line names exactly this
    producer's pair, matched whole: no path is split out of the line."""
    if "verdict" in entry:                                   # a particle's line
        return (OWN + "%s and %s, " % pair) in entry["text"]
    return entry["text"].endswith("; its particles from the producer's own %s and %s" % pair)


def pf_particles(out):
    """Each sampled particle's line, by (member, body): its verdict and what
    it says, the source it names (the phrase that comes first), its counters
    and its closest approaches."""
    res = {}
    for m in PF_PART.finditer(out):
        rest = m.group(5)
        v = PF_VERDICT.search(rest)
        at = {src: rest.find(phrase) for src, phrase in (("producer", OWN), ("rerun", RERUNS)) if phrase in rest}
        res[(int(m.group(1)), int(m.group(3)))] = {
            "name": m.group(2), "why": m.group(4), "text": rest,
            "verdict": v.group(1) if v else "FAIL" if rest.startswith("FAIL") else "?",
            "what": v.group(2) if v else rest,
            "counters": dict(zip(SAID, v.groups()[2:7])) if v else {},
            "closest": [(c.group(1), c.group(2), int(c.group(3))) for c in CLOSE.finditer(v.group(8))] if v else [],
            "source": min(at, key=at.get) if at else None,
            "zonal": " massive bodies and its zonal line through " in rest}
    return res


def limit_once(out):
    """The limit of the producer's files printed as one line, word for word,
    after the plan's line, and as no other line. Lines are counted, not
    substrings: two adjacent copies share a newline, and a count of
    "\\n" + line + "\\n" saw them as one (V3, 2026-11-18)."""
    lines = out.splitlines()
    plan = next((i for i, l in enumerate(lines) if l.startswith(PF_PLAN[1:])), None)
    return lines.count(PRODUCER_LIMIT) == 1 and plan is not None and lines.index(PRODUCER_LIMIT) > plan and \
        sum(1 for l in lines if l.startswith(LIMIT_OPENS)) == 1


def limit_absent(out):
    """No line of the limit of a producer's files, in any wording: a run
    without --producer-files reads no producer's files (V3, 2026-11-18: the
    line printed on the default path too passed every case)."""
    return not any(l.startswith(LIMIT_OPENS) for l in out.splitlines())


def replays(keep):
    """The replays' records the driver kept, <keep>/particles/<stem>.rec, by
    (member, body), from tools/extract_particle.py's stem."""
    d = os.path.join(keep, "particles")
    res = {}
    for f in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        m = re.search(r"-body(\d+)-m(\d+)\.rec$", f)
        if m:
            res[(int(m.group(2)), int(m.group(1)))] = os.path.join(d, f)
    return res


def replay_view(path):
    """A replay's record as two of them are compared here: its sample lines,
    its trailer's counters, flags and calls, and the dt file it names."""
    with open(path) as f:
        lines = f.read().splitlines()
    tail = next((l for l in lines if l.startswith("# steps_done=")), "# ")
    kvs = dict(t.partition("=")[::2] for t in tail[2:].split())
    dtf = re.search(r" dt_file=(.*)$", next((l for l in lines if l.startswith("# program=")), ""))
    return ([l for l in lines if l.startswith("sample ")],
            {c: kvs.get(c) for c in TRAILER + ("flags_seen", "calls", "divsqrt_calls")}, dtf.group(1) if dtf else None)


def plant_dt(tok):
    """A hex float with its mantissa's last hex digit changed (XOR 1), a
    fraction digit made first where it has none."""
    mant, _, exp = tok.partition("p")
    if "." not in mant:
        mant += ".0"
    return mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + "p" + exp


def body_lines(lines):
    """The indices of a producer's file's lines that are not comments."""
    return [i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")]


def edit_file(path, fn):
    with open(path) as f:
        lines = f.read().splitlines()
    got = fn(lines)
    with open(path, "w", newline="\n") as f:
        f.write("\n".join(lines) + "\n")
    return got


def case7(clean, fmts, waves, scratch, cft, progs, quick, record, drive, run, verdict, trailers, member_lines):
    """CASE 7, the producer's own attempts records and schedules (ROUND5.md,
    P3): the head of this file says what each part holds."""
    print("check_compare_run: case 7 - the producer's own attempts records and schedules (ROUND5.md, P3)")
    fmts = [f for f in fmts if f in clean]
    if "fp64" not in fmts:
        verdict(False, "case 7: case 1's binary64 record, to take a producer's files beside")
        return
    empty = os.path.join(scratch, "producer_none")
    shutil.rmtree(empty, ignore_errors=True)
    os.makedirs(empty)

    def member_run(problem, fmt, steps, sample, k, extra=()):
        return run([cft, "--format", fmt, "--problem", problem, "--dt", DT0, "--epsilon", "1e-9", "--steps", str(steps),
                    "--sample", str(sample), "--max-iter", "60", "--member", str(k), "--quiet"] + list(extra))

    def produce(pdir, problem, fmt, steps, sample, text, how=()):
        """ias15_cft as the producer: each member run alone (`how`: its
        arithmetic and engine, if not the defaults), writing its own attempts
        record and schedule into pdir by the driver's names. Returns the
        members whose run is not the record's member, line for line and
        counter for counter, so whose files would not be the record's run's."""
        shutil.rmtree(pdir, ignore_errors=True)
        os.makedirs(pdir)
        off = []
        for k in range(E):
            att, sch = pf(pdir, k)
            rc, out, _ = member_run(problem, fmt, steps, sample, k,
                                    list(how) + ["--attempts-out", att, "--record-iters", sch])
            ms, mt = parse_samples(out) if rc == 0 else ([], {})
            es, et = parse_samples(text, k)
            if rc != 0 or not ms or ms != es or any(mt.get(c) != et.get(c) for c in TRAILER):
                off.append(k)
        return off

    def pf_drive(rec, fmt, pdir, keep=None, problem=None):
        extra = ["--producer-files", pdir, "--particle-strata", str(E)]
        if keep:
            shutil.rmtree(keep, ignore_errors=True)
            extra += ["--keep", keep]
        return drive(rec, fmt, extra=extra, problem=problem)

    def mp_closest(rec_path, k, j):
        """Particle j of member k: its least distance to the primary and to the
        moon over the samples, and at which sample, in mpmath from
        tools/oracle.py's parse of the record."""
        _, _, smps, _ = oracle_parse(rec_path, system=k)
        res = []
        for i in (0, 1):
            dists = [mpmath.sqrt(sum((F2m(s["state"][6 * j + c]) - F2m(s["state"][6 * i + c])) ** 2 for c in range(3)))
                     for s in smps]
            best = min(dists)
            res.append((best, dists.index(best)))
        return res

    def refused(label, rec, pdir, k, what, want, exact=True, problem=None):
        """A copy of the producer's files whose member k's file is another
        run's: REFUSED by its key, naming the member, the file and the key
        fields that differ; then the copy is deleted."""
        rc, out = pf_drive(rec, "fp64", pdir, problem=problem)
        m = PF_REFUSED.search(out)
        fields = m.group(4).split(", ") if m else []
        ok = (rc == 2 and m is not None and int(m.group(1)) == k and m.group(2) == what and
              same_path(m.group(3), pf(pdir, k)[0 if what == "attempts record" else 1]) and
              (fields == want if exact else set(want) <= set(fields)) and "IDENTICAL" not in out)
        shutil.rmtree(pdir, ignore_errors=True)
        verdict(ok, "%s: REFUSED by its run key before anything is compared, naming member %d's %s and the key fields "
                "that differ, %s (rc %d); the copy deleted%s"
                % (label, k, what, ", ".join(fields) or "none", rc,
                   "" if ok else ": " + (re.findall(r"^RESULT: .*$", out, re.M) or [out[-400:]])[-1][:500]))

    def fails_of(p, k):
        return [d for x, d in sorted(p.items()) if x[0] == k]

    # ------------------------------------------------ gates 1 and 2, each format
    pdirs, alone_of = {}, {}
    for fmt in fmts:
        c = clean[fmt]
        smp, text, rec = c["smp"], c["text"], c["rec"]
        # " and " in the name, the word the driver's lines put between the
        # pair's two paths: every reading of those lines here must hold
        pdir = os.path.join(scratch, "producer %s and its pairs" % fmt)
        pdirs[fmt] = pdir
        off = produce(pdir, waves, fmt, c["steps"], c["sample"], text)
        verdict(not off, "%s: ias15_cft as the producer, each of the %d members run alone with --attempts-out and "
                "--record-iters, writing %s and %s (the driver's names): each run's lines and counters are the record's "
                "member's, so the files are the record's run's%s"
                % (fmt, E, PRODUCER_FILES[0][1], PRODUCER_FILES[1][1], "; NOT for members %s" % off if off else ""))
        keep_a, keep_b = (os.path.join(scratch, "producer_keep_%s_%s" % (t, fmt)) for t in ("own", "rerun"))
        rca, outa = pf_drive(rec, fmt, pdir, keep_a)
        rcb, outb = pf_drive(rec, fmt, empty, keep_b)
        alone = [k for k in range(E) if k not in smp["members"]]
        alone_of[fmt] = alone
        ma, pa, mb, pb, lb = pf_members(outa), pf_particles(outa), pf_members(outb), pf_particles(outb), member_lines(outb)
        n = 2 * E
        # from the producer's own files
        not_rerun = sorted(int(m.group(1)) for m in PF_ALONE.finditer(outa))
        bad = ["member %d: %s" % (k, ma[k]["text"][:200] if k in ma else "no line") for k in range(E)
               if k not in ma or not names_pair(ma[k], pf(pdir, k)) or ma[k]["rerun"] != (k in smp["members"])]
        bad += ["member %d body %d: %s" % (x + (d["text"][:200],)) for x, d in sorted(pa.items())
                if d["verdict"] != "IDENTICAL" or d["source"] != "producer" or not names_pair(d, pf(pdir, x[0]))]
        summary = ("\nparticles: %d sampled in %d member(s), each replayed alone and held to the record's own lines: %d "
                   "from the producer's own attempts records and schedules, in %d member(s), %d of them not rerun; 0 "
                   "from this tree's reruns', where the producer's pair is not there; %d identical\n"
                   % (n, E, n, E, len(alone), n))
        keyed = "\nparticles: %d producer's pair(s), each keyed to its member's run: " % E
        verdict(rca == 0 and "\nRESULT: PASS: " in outa and not bad and len(pa) == n and not_rerun == alone and alone
                and summary in outa and keyed in outa and limit_once(outa),
                "%s: from the producer's own files, in a directory named %r: PASS (rc %d); all %d members' particles "
                "sampled (--particle-strata %d), each pair held to its member's run by its key and matched whole, %d "
                "particles replayed from the producer's own attempts records and schedules, all IDENTICAL to the "
                "record; members %s NOT RERUN, the rest rerun; beside the plan, once, word for word, the limit: '%s'%s"
                % (fmt, os.path.basename(pdir), rca, E, E, len(pa), not_rerun, PRODUCER_LIMIT,
                   ": " + "; ".join(bad[:3]) if bad else ""))
        # an absent pair: today's rerun, said by name
        badb = ["member %d: %s" % (k, mb[k]["text"][:200] if k in mb else "no line") for k in range(E)
                if k not in mb or not mb[k]["rerun"] or mb[k]["producer"] or
                ("the producer's pair is not there (no %s and no %s)" % pf(empty, k)) not in mb[k]["text"]]
        badb += ["member %d body %d: %s" % (x + (d["text"][:200],)) for x, d in sorted(pb.items())
                 if d["verdict"] != "IDENTICAL" or d["source"] != "rerun"]
        verdict(rcb == 0 and "\nRESULT: PASS: " in outb and not badb and len(pb) == n and limit_once(outb) and
                sorted(lb) == list(range(E)) and all(v[0] == "IDENTICAL" for v in lb.values()),
                "%s: with no producer's pair (an empty directory), each member falls back to today's rerun, naming both "
                "files that are not there: PASS (rc %d), all %d members rerun and IDENTICAL, %d particles replayed from "
                "the reruns' files, IDENTICAL%s" % (fmt, rcb, len(lb), len(pb), ": " + "; ".join(badb[:3]) if badb else ""))
        # GATE 1: the same particles, replay record against replay record
        ra, rb = replays(keep_a), replays(keep_b)
        diffs, own, theirs, ntok = [], 0, 0, 0
        for x in sorted(pa):
            if x not in pb or (pa[x]["name"], pa[x]["why"]) != (pb[x]["name"], pb[x]["why"]):
                diffs.append("member %d body %d: not picked alike in the two runs" % x)
                continue
            if x not in ra or x not in rb:
                diffs.append("member %d body %d: no replay record kept" % x)
                continue
            (sa, ca, da), (sb, cb, db) = replay_view(ra[x]), replay_view(rb[x])
            if not sa or sa != sb:
                diffs.append("member %d body %d: %s" % (x + ("its sample lines differ" if sa else "no sample line",)))
            elif ca != cb:
                diffs.append("member %d body %d: counters %s and %s" % (x + (ca, cb)))
            ntok += sum(len(l.split()) for l in sa)
            own += same_path(da, pf(pdir, x[0])[0])
            theirs += db is not None and same_path(os.path.dirname(db), keep_b)
        first = [k for k in c["first"] if any(x[0] == k for x in pa)]
        verdict(not diffs and len(pa) == n and own == n and theirs == n and first,
                "%s: GATE 1, the same particles: each of the %d, picked alike in both runs, replayed from the producer's "
                "own files is IDENTICAL to the same particle replayed from a rerun's - every sample line, %d tokens, and "
                "the trailer's counters, flags and calls; each replay's record names the attempts record it ran from, "
                "the producer's in %d of %d and the rerun's in %d of %d; members %s rejected their first attempt, the "
                "rejection that leaves state, and their particles are among them%s"
                % (fmt, len(pa), ntok, own, n, theirs, n, first, ": " + "; ".join(diffs[:3]) if diffs else ""))
        # GATE 2: the physics beside each particle
        tr = trailers(text)
        g2, worst, nclose = [], mpmath.mpf(0), 0
        for (k, j), d in sorted(pa.items()):
            if any(d["counters"].get(cn) != tr[k].get(cn) for cn in TRAILER):
                g2.append("member %d body %d: counters %s, its trailer's %s" % (k, j, d["counters"], tr[k]))
            if len(d["closest"]) != 2:
                g2.append("member %d body %d: closest approach %r" % (k, j, d["closest"]))
                continue
            for (nm, val, at), (best, at_mp) in zip(d["closest"], mp_closest(rec, k, j)):
                nclose += 1
                worst = max(worst, abs(mpmath.mpf(val) - best) / best)
                if at != at_mp:
                    g2.append("member %d body %d: closest to %s at sample %d, by mpmath at %d" % (k, j, nm, at, at_mp))
        for m in PF_ALONE.finditer(outa):
            k = int(m.group(1))
            if any(dict(zip(SAID, m.groups()[1:])).get(cn) != tr[k].get(cn) for cn in TRAILER):
                g2.append("member %d, not rerun: its trailer line says %s" % (k, m.groups()[1:]))

        def span(key):
            vals = sorted({d["counters"].get(key) for d in pa.values()}, key=lambda v: float(v or "nan"))
            return vals[0] if len(vals) == 1 else "%s to %s" % (vals[0], vals[-1])
        verdict(not g2 and nclose == 2 * n and worst < 1e-6,
                "%s: GATE 2, the physics beside each: every particle's replay counters are its member's trailer (steps "
                "%s, rejected %s, mean passes %s, max passes %s, cap hits %s), and so is each not-rerun member's trailer "
                "line; its %d closest approaches, each to the primary and the moon at the samples, agree with this "
                "gate's mpmath from tools/oracle.py's parse to %s relative, at the samples they name%s"
                % (fmt, span("steps_done"), span("steps_rejected"), span("mean_pc_iterations"),
                   span("max_pc_iterations"), span("iterations_max_exceeded"), nclose, mpmath.nstr(worst, 3),
                   ": " + "; ".join(g2[:3]) if g2 else ""))

    # ----------------------- GATE 1b: a producer whose exit test is not this tree's
    # Gate 1 compares ias15_cft with itself, so the producer's files are a
    # rerun's to the byte and a driver that read a rerun's schedule instead
    # would pass it. Here the producer's first attempt took other passes:
    # two fewer in a member that is not rerun (it moves bits), one fewer in
    # one that is (it moves no bit at the samples, only trailer counters -
    # which ones is derived and printed below: V3, 2026-11-16, found the max
    # moving with the mean where this once said "only the mean passes"),
    # each member's record lines that schedule's replay. The producer's
    # files must hold its particles to its record, and a rerun's must not.
    c = clean["fp64"]
    smp, text = c["smp"], c["text"]
    alone = alone_of.get("fp64", [])
    u = alone[0] if alone else None
    m = next((k for k in smp["picks"] if k not in smp["rejected"]), None)
    if u is None or m is None:
        verdict(False, "GATE 1b: a member that is not rerun (%s) and a stratum's pick that rejected nothing (%s)" % (u, m))
    else:
        disc = os.path.join(scratch, "producer_disc")
        shutil.rmtree(disc, ignore_errors=True)
        shutil.copytree(pdirs["fp64"], disc)
        rows = text.splitlines()
        moved = {}
        for k, delta in ((u, -2), (m, -1)):
            att, sch = pf(disc, k)

            def fewer(lines, delta=delta):
                i = body_lines(lines)[0]
                was = int(lines[i])
                lines[i] = str(was + delta)
                return was, was + delta
            was, now = edit_file(sch, fewer)
            rc, out, err = member_run(waves, "fp64", c["steps"], c["sample"], k,
                                      ["--replay-iters", sch, "--attempts-out", att])
            ms, mt = parse_samples(out) if rc == 0 else ([], {})
            es, et = parse_samples(text, k)
            # what the change moves in the trailer, read here, counter by counter
            counters = ", ".join("%s %s where this tree's is %s" % (cn, mt.get(cn), et.get(cn)) for cn in TRAILER
                                 if mt.get(cn) != et.get(cn)) or "no counter"
            moved[k] = (was, now, rc, sum(1 for a, b in zip(es, ms) if a != b), len(es), counters)
            mine = iter([l for l in out.splitlines() if l.startswith("sample ")])
            rows = ["system %d %s" % (k, next(mine)) if l.startswith("system %d sample " % k) else
                    "# system %d %s" % (k, " ".join("%s=%s" % (cn, mt.get(cn)) for cn in (
                        "steps_done", "iterations_max_exceeded", "steps_rejected", "mean_pc_iterations",
                        "max_pc_iterations"))) if l.startswith("# system %d steps_done=" % k) else l for l in rows]
        tot = next(i for i, l in enumerate(rows) if l.startswith("# steps_done="))
        for cn in ("steps_rejected", "iterations_max_exceeded"):
            n_all = sum(int(re.search(r"\b%s=(\d+)" % cn, l).group(1)) for l in rows if re.match(r"# system \d+ steps_done=", l))
            rows[tot] = re.sub(r"\b%s=\d+" % cn, "%s=%d" % (cn, n_all), rows[tot])
        discrec = os.path.join(scratch, "producer_disc.rec")
        with open(discrec, "w", newline="\n") as f:
            f.write("\n".join(rows) + "\n")
        (wu, nu, rcu, du, lu, cu), (wm, nm, rcm, dm, lm, cm) = moved[u], moved[m]
        verdict(rcu == 0 and rcm == 0 and du > 0 and (dm > 0 or cm != "no counter"),
                "GATE 1b, not against itself: the producer: member %d's first attempt took %d passes where this tree's "
                "takes %d, and its record lines are that schedule's replay, %d of %d of them not this tree's (and %s); "
                "member %d's took %d where this tree's takes %d, which moves %d of %d lines and %s"
                % (u, nu, wu, du, lu, cu, m, nm, wm, dm, lm, cm))
        rc1, out1 = pf_drive(discrec, "fp64", disc)
        rc2, out2 = pf_drive(discrec, "fp64", empty)
        p1, p2, l1, l2 = pf_particles(out1), pf_particles(out2), member_lines(out1), member_lines(out2)
        alone1 = [int(x.group(1)) for x in PF_ALONE.finditer(out1)]
        ok1 = (rc1 == 1 and u in alone1 and u not in l1 and l1.get(m, ("",))[0] == "FAIL" and len(p1) == 2 * E and
               all(v[0] == "IDENTICAL" for k, v in l1.items() if k != m) and limit_once(out1) and
               all(d["verdict"] == "IDENTICAL" and d["source"] == "producer" for d in p1.values()))
        verdict(ok1, "GATE 1b, from the producer's own files: member %d, NOT RERUN, holds - both its particles IDENTICAL to "
                "the producer's record, whose lines are not this tree's: the limit, measured, and printed beside the "
                "plan word for word (the producer's passes taken as given); member %d FAILs on its rerun (%s), while "
                "both its particles, replayed against the producer's schedule, are IDENTICAL; every other member and "
                "particle IDENTICAL (rc %d)" % (u, m, l1.get(m, ("", "not reported"))[1][:110], rc1))
        bad2 = ["member %d body %d %s" % (x + (d["verdict"],)) for x, d in sorted(p2.items())
                if (d["verdict"] == "FAIL") != (x[0] in (u, m)) or d["source"] != "rerun"]
        why2 = ["member %d: %s" % (k, fails_of(p2, k)[0]["what"][:90]) for k in (u, m) if fails_of(p2, k)]
        verdict(rc2 == 1 and l2.get(u, ("",))[0] == "FAIL" and l2.get(m, ("",))[0] == "FAIL" and len(p2) == 2 * E
                and not bad2 and limit_once(out2),
                "GATE 1b, the same record replayed from this tree's reruns instead (no producer's pair): members %d and %d "
                "FAIL, and all 4 of their particles FAIL (%s): so the replays above read the producer's own schedule, "
                "and the same-particles gate cannot pass against itself (rc %d)%s"
                % (u, m, "; ".join(why2), rc2, ": " + "; ".join(bad2[:3]) if bad2 else ""))
        shutil.rmtree(disc, ignore_errors=True)

    # ------------------------------- half a pair: the fallback names the file
    alone = alone_of.get("fp64", [])
    if alone:
        u = alone[0]
        half = os.path.join(scratch, "producer_half")
        shutil.rmtree(half, ignore_errors=True)
        shutil.copytree(pdirs["fp64"], half)
        os.remove(pf(half, u)[1])
        rch, outh = pf_drive(c["rec"], "fp64", half)
        mh, ph, lh = pf_members(outh), pf_particles(outh), member_lines(outh)
        n = 2 * E
        summary = ("\nparticles: %d sampled in %d member(s), each replayed alone and held to the record's own lines: %d "
                   "from the producer's own attempts records and schedules, in %d member(s), %d of them not rerun; 2 "
                   "from this tree's reruns', where the producer's pair is not there; %d identical\n"
                   % (n, E, n - 2, E - 1, len(alone) - 1, n))
        verdict(rch == 0 and "\nRESULT: PASS: " in outh and u in mh and mh[u]["rerun"] and limit_once(outh) and
                ("the producer's pair is not there (no %s): its particles fall back to today's rerun" % pf(half, u)[1])
                in mh[u]["text"] and lh.get(u, ("",))[0] == "IDENTICAL" and summary in outh and
                all(d["verdict"] == "IDENTICAL" and d["source"] == ("rerun" if x[0] == u else "producer")
                    for x, d in ph.items()) and len(ph) == n,
                "fp64, member %d's schedule taken away, its attempts record left: it falls back to today's rerun, "
                "naming the schedule alone as not there; its lines rerun, IDENTICAL, and its 2 particles from the "
                "rerun's files, the other %d from the producer's, all IDENTICAL: PASS (rc %d)" % (u, n - 2, rch))
        shutil.rmtree(half, ignore_errors=True)

    # ------------- one system's records: a --member k run's, and a single problem
    k1 = alone[0] if alone else 0
    rc, t1, _ = member_run(waves, "fp64", c["steps"], c["sample"], k1)
    rec1 = os.path.join(scratch, "producer_member%d.rec" % k1)
    with open(rec1, "w", newline="\n") as f:
        f.write(t1)
    rc1, out1 = pf_drive(rec1, "fp64", pdirs["fp64"])
    m1, p1 = pf_members(out1), pf_particles(out1)
    verdict(rc == 0 and rc1 == 0 and "\nRESULT: PASS: " in out1 and sorted(m1) == [k1] and m1[k1]["rerun"] and
            names_pair(m1[k1], pf(pdirs["fp64"], k1)) and len(p1) == 2 and limit_once(out1) and
            all(d["verdict"] == "IDENTICAL" and d["source"] == "producer" for d in p1.values()),
            "fp64, a --member %d record of the waves, one system: its pair is member %d's, %s; it is rerun as the "
            "record's one system, and its 2 particles replayed from the producer's files, IDENTICAL: PASS (rc %d)"
            % (k1, k1, " and ".join(os.path.basename(x) for x in pf(pdirs["fp64"], k1)), rc1))
    # a single problem, no system lines, with an active line: its run has no
    # --member, so its key says member=-1, and its pair is member0's (P2,
    # 2026-11-16: `--member 0` is refused on such a file)
    with open(waves) as f:
        wl = f.read().splitlines()
    one = os.path.join(scratch, "cr_single.txt")
    s0 = wl.index("system w0")
    with open(one, "w", newline="\n") as f:
        f.write("\n".join(["# check_compare_run: system w0 of the waves, as a single problem", "name cr_single",
                           "G " + float.hex(1.0)] + wl[s0 + 1:s0 + 3 + NTP] + ["active 2"]) + "\n")
    ps = os.path.join(scratch, "producer_single")
    shutil.rmtree(ps, ignore_errors=True)
    os.makedirs(ps)
    att, sch = pf(ps, 0)
    rc, ts, _ = run([cft, "--format", "fp64", "--problem", one, "--dt", DT0, "--epsilon", "1e-9", "--steps",
                     str(c["steps"]), "--sample", str(c["sample"]), "--max-iter", "60", "--quiet",
                     "--attempts-out", att, "--record-iters", sch])
    recs = os.path.join(scratch, "producer_single.rec")
    with open(recs, "w", newline="\n") as f:
        f.write(ts)
    keys = []
    for path in (att, sch):
        with open(path) as f:
            keys.append(f.read().splitlines()[1] if rc == 0 else "")
    rcs, outs = pf_drive(recs, "fp64", ps, problem=one)
    p1 = pf_particles(outs)
    verdict(rc == 0 and rcs == 0 and "\nRESULT: PASS: " in outs and all(" member=-1 " in k for k in keys) and
            len(p1) == 2 and limit_once(outs) and
            all(d["verdict"] == "IDENTICAL" and d["source"] == "producer" for d in p1.values()),
            "fp64, a single problem with no system lines and `active 2` (the waves' w0): its run has no --member, so "
            "both its files say member=-1 in their run key, and its pair is %s; PASS, its 2 particles replayed from "
            "them, IDENTICAL (rc %d)" % (" and ".join(os.path.basename(x) for x in (att, sch)), rcs))
    shutil.rmtree(ps, ignore_errors=True)
    kep = os.path.join(root, "data", "problems", "kepler.txt")
    rck, tk, _ = run([cft, "--format", "fp64", "--problem", kep, "--dt", "0.05", "--epsilon", "1e-9", "--steps", "20",
                      "--sample", "10", "--max-iter", "60", "--quiet"])
    reck = os.path.join(scratch, "producer_kepler.rec")
    with open(reck, "w", newline="\n") as f:
        f.write(tk)
    rck2, outk = drive(reck, "fp64", extra=["--producer-files", empty], problem=kep)
    verdict(rck == 0 and rck2 == 0 and "\nRESULT: PASS: " in outk and
            ("particles: the producer's files in %s are not read: the record has no `# active` line, so it has no test "
             "particle" % empty) in outk and "member 0 (kepler): IDENTICAL: " in outk and limit_absent(outk),
            "fp64, kepler, a single problem with every body active: the producer's files are not read, and it says so, "
            "printing no limit of files it never read; its one system IDENTICAL: PASS (rc %d)" % rck2)

    # ------------------- a record another producer wrote, naming its own engine
    # P2's kernel writes this tree's record format under its own program and
    # engine: its files are keyed to the engine the caller reruns it on
    # (--rerun-engine), which this tree has. Case 1's record rewritten so.
    other = re.sub(r"^# program=cft impl=\S+ (format=\S+) backend=\S+ ", r"# program=wave impl=hip-7.16 \1 backend=hip:0 ",
                   text, count=1, flags=re.M)
    other = re.sub(r"^(# program=wave .*?) engine=loop ", r"\1 engine=wave ", other, count=1, flags=re.M)
    orec = os.path.join(scratch, "producer_other.rec")
    with open(orec, "w", newline="\n") as f:
        f.write(other)
    rco, outo = drive(orec, "fp64", extra=["--producer-files", pdirs["fp64"], "--particle-strata", str(E),
                                           "--rerun-engine", "loop"])
    po = pf_particles(outo)
    keyed = "\nparticles: %d producer's pair(s), each keyed to its member's run: its run key, both files', is the one " \
            "this tree's ias15_cft --run-key gives the member's arguments, at engine=loop\n" % E
    verdict(other.count(" engine=wave ") == 1 and "program=wave impl=hip-7.16 format=fp64 backend=hip:0" in other and
            rco == 0 and "\nRESULT: PASS: " in outo and keyed in outo and len(po) == 2 * E and limit_once(outo) and
            "record: written by program=wave impl=hip-7.16 backend=hip:0" in outo and
            all(d["verdict"] == "IDENTICAL" and d["source"] == "producer" for d in po.values()),
            "fp64, the record rewritten as another producer's (program=wave impl=hip-7.16 backend=hip:0 engine=wave), "
            "held with --rerun-engine loop: its files keyed at engine=loop, the engine it is rerun on; PASS (rc %d), %d "
            "particles from the producer's own files, all IDENTICAL" % (rco, len(po)))
    rco, outo = drive(orec, "fp64", extra=["--producer-files", pdirs["fp64"], "--particle-strata", str(E)])
    verdict(rco == 2 and "RESULT: REFUSED: the record ran engine=wave, which this tree does not have, so its producer's "
            "files are keyed to the engine its rerun takes, and no --rerun-engine names one" in outo and
            "IDENTICAL" not in outo,
            "fp64, the same with no --rerun-engine: REFUSED, naming the engine its files' key cannot name (rc %d)" % rco)
    # The FMA form's program-engine record, rerun on the loop engine: its
    # producer's files are keyed to the engine it ran, program, and its
    # particles are replayed there, while its members are rerun on loop. In
    # every tier, short in --quick: only this case tells the record's engine
    # from --rerun-engine's, and a driver keying the files to the second
    # (V3's mutant Mj) passed every case --quick ran (V3, 2026-11-16).
    steps, sample = (4, 2) if quick else (20, 10)
    how = ["--arith", "fma", "--engine", "program", "--programs", progs]
    recp, textp = record("producer_fma_program.rec", "fp64", steps, sample, arith="fma", engine="program")
    pp = os.path.join(scratch, "producer_fma_program")
    offp = produce(pp, waves, "fp64", steps, sample, textp, how) if recp else list(range(E))
    rcp, outp = drive(recp, "fp64", "fma", extra=["--producer-files", pp, "--particle-strata", str(E),
                                                  "--rerun-engine", "loop"]) if recp else (-1, "")
    ppp = pf_particles(outp)
    keyed = "\nparticles: %d producer's pair(s), each keyed to its member's run: its run key, both files', is the " \
            "one this tree's ias15_cft --run-key gives the member's arguments, at engine=program\n" % E
    verdict(not offp and rcp == 0 and "\nRESULT: PASS: " in outp and keyed in outp and len(ppp) == 2 * E and
            "rerun: on the loop engine, where the record ran program (--rerun-engine)" in outp and limit_once(outp) and
            all(d["verdict"] == "IDENTICAL" and d["source"] == "producer" for d in ppp.values()),
            "fp64, the FMA form's program-engine record (%d steps) with --rerun-engine loop: its members rerun on the "
            "loop engine, its producer's files keyed at engine=program, the engine it ran, and its %d particles "
            "replayed there, all IDENTICAL: PASS (rc %d)%s"
            % (steps, len(ppp), rcp, "" if rcp == 0 else ": " + (re.findall(r"^RESULT: .*$", outp, re.M) or [""])[-1][:300]))
    shutil.rmtree(pp, ignore_errors=True)

    # ---------------------------------- the waves in their primary's zonal field
    zl = "0 %s %s %s" % (float.hex(0.25), float.hex(2.0 ** -6), float.hex(-2.0 ** -10))
    wz = write_waves(os.path.join(scratch, "cr_waves_z.txt"), name="cr_waves_z", zonal=zl)
    recz, textz = record("waves_z_fp64.rec", "fp64", c["steps"], c["sample"], problem=wz)
    pz = os.path.join(scratch, "producer_z")
    if recz is None:
        verdict(False, "fp64: the zonal waves' record")
    else:
        offz = produce(pz, wz, "fp64", c["steps"], c["sample"], textz)
        rcz, outz = pf_drive(recz, "fp64", pz, problem=wz)
        pzp = pf_particles(outz)
        zrec = re.search(r"^# zonal .*$", textz, re.M)
        verdict(not offz and rcz == 0 and "\nRESULT: PASS: " in outz and len(pzp) == 2 * E and zrec is not None and
                limit_once(outz) and
                all(d["verdict"] == "IDENTICAL" and d["source"] == "producer" and d["zonal"] for d in pzp.values()),
                "fp64, the waves with `zonal %s`, the primary's J2 and J4: from the producer's own files, PASS (rc %d), "
                "%d particles replayed alone with the member's 2 massive bodies and its zonal line, all IDENTICAL; the "
                "record says '%s'" % (zl, rcz, len(pzp), zrec.group(0)[:90] if zrec else "no # zonal line"))

    # --------------------------------------------- the controls on the files
    print("check_compare_run: case 7's controls, each on a copy of the producer's files, which must FAIL or be "
          "refused by name, and is then deleted")
    for fmt in fmts:
        c = clean[fmt]
        smp, rec = c["smp"], c["rec"]
        alone = alone_of[fmt]
        first = c["first"]
        pool = alone + [k for k in smp["picks"] if k not in smp["rejected"] and k not in alone]
        # q1: one attempt's dt, a member a line, the first attempt included;
        # and the other column, a rejected first attempt recorded as accepted
        plants = []
        if first:
            plants.append((first[0], "the first attempt's dt, a rejection before any accepted step", "dt", "first"))
        for (what, col, where), k in zip((("the first attempt's dt, accepted", "dt", "first"),
                                          ("a middle accepted attempt's dt", "dt", "middle"),
                                          ("the last line's dt, the step after its last attempt", "dt", "next")),
                                         [k for k in pool if k not in first]):
            plants.append((k, what, col, where))
        if len(first) > 1:
            plants.append((first[1], "the outcome of its first attempt, a rejection, made an acceptance", "outcome",
                           "first"))
        q1 = os.path.join(scratch, "producer_q1_%s" % fmt)
        shutil.rmtree(q1, ignore_errors=True)
        shutil.copytree(pdirs[fmt], q1)
        said = {}
        for k, what, col, where in plants:
            def change(lines, col=col, where=where):
                body = body_lines(lines)
                tried = body[:-1]
                i = {"first": tried[0], "next": body[-1]}.get(where)
                if i is None:
                    acc = [i for i in tried if lines[i].split()[1] == "a"]
                    i = acc[len(acc) // 2]
                dt, outcome = lines[i].split()
                lines[i] = "%s %s" % ((plant_dt(dt), outcome) if col == "dt" else (dt, "a" if outcome == "r" else "r"))
                return i + 1, dt + " " + outcome, lines[i]
            said[k] = edit_file(pf(q1, k)[0], change)
        rc, out = pf_drive(rec, fmt, q1)
        p, l = pf_particles(out), member_lines(out)
        caught, missed = [], []
        for k, what, col, where in plants:
            got = fails_of(p, k)
            if len(got) == 2 and all(d["verdict"] == "FAIL" and ("record line " in d["what"] or "did not run" in d["what"])
                                     for d in got):
                caught.append("member %d, %s (line %d, '%s' made '%s'): %s" % ((k, what) + said[k] + (got[0]["what"][:100],)))
            else:
                missed.append("member %d, %s: %s" % (k, what, [d["verdict"] + " " + d["what"][:80] for d in got]))
        calm = [x for x, d in p.items() if x[0] not in said and d["verdict"] != "IDENTICAL"]
        shutil.rmtree(q1, ignore_errors=True)
        verdict(rc == 1 and len(plants) >= 4 and not missed and not calm and l and limit_once(out) and
                all(v[0] == "IDENTICAL" for v in l.values()),
                "%s: q1, one change planted in the attempts record, a member a line, %d lines: each FAILs both its "
                "member's particles by name, every other particle and every rerun member IDENTICAL; the copy deleted "
                "(rc %d): %s%s" % (fmt, len(plants), rc, "; ".join(caught),
                                   "; NOT: " + "; ".join(missed + ["member %d body %d not IDENTICAL" % x for x in calm])
                                   if missed or calm else ""))
        # q2: one pass count, a member a line, in members rerun and not
        plants = []
        for (what, where, delta), k in zip((("the first attempt's passes, two fewer", "first", -2),
                                            ("a middle attempt's passes, one more", "middle", 1)), alone):
            plants.append((k, what, where, delta))
        mm = next((k for k in pool if k not in alone), None)
        if mm is not None:
            plants.append((mm, "the first attempt's passes, one fewer, in a member rerun", "first", -1))
        # and a rejected attempt's passes: a rejection before any accepted
        # step leaves its corrected b behind, so its passes reach the state
        # (V3, 2026-11-16: with no plant here, a driver reading a rerun's
        # schedule for rejecting members alone passed q2, its mutant Mf2)
        if first:
            plants.append((first[0], "its rejected first attempt's passes, two fewer, in a member rerun", "first", -2))
        q2 = os.path.join(scratch, "producer_q2_%s" % fmt)
        q2keep = os.path.join(scratch, "producer_q2_keep_%s" % fmt)
        shutil.rmtree(q2, ignore_errors=True)
        shutil.copytree(pdirs[fmt], q2)
        said = {}
        for k, what, where, delta in plants:
            def change(lines, where=where, delta=delta):
                body = body_lines(lines)
                i = body[0] if where == "first" else body[len(body) // 2]
                was = lines[i]
                lines[i] = str(int(was) + delta)
                return i + 1, was, lines[i]
            said[k] = edit_file(pf(q2, k)[1], change)
        rc, out = pf_drive(rec, fmt, q2, keep=q2keep)
        p, l = pf_particles(out), member_lines(out)
        # what each plant moves, derived: its member's replays against gate
        # 1's unplanted ones of the same particles (tokens), and its replay's
        # counters against the record's trailer (V3: "only the mean passes"
        # was said of plants that moved the max too). Printed as measured:
        # the gate asserts only that each plant FAILs its particles and moves
        # something, never which things at which format.
        gate1 = os.path.join(scratch, "producer_keep_own_%s" % fmt)
        rq, ra, tr = replays(q2keep), replays(gate1), trailers(c["text"])

        def moves(k):
            """(tokens, counters, what is missing): a replay this run or gate 1
            kept for one of member k's sampled particles, if one is not there,
            is a FAIL by name, never a count of -1 (V3, 2026-11-18)."""
            js = sorted(x for x in p if x[0] == k)
            lost = ["no replay of member %d body %d kept by %s" % (x + (by,)) for x in js
                    for by, kept in (("this run", rq), ("gate 1, in %s" % gate1, ra)) if x not in kept]
            if len(js) != 2 or lost:
                return None, None, "what the plant moves cannot be read: %s" % (
                    "; ".join(lost) if lost else "%d sampled particles, not 2" % len(js))
            ntok = max(sum(1 for a, b in zip(" ".join(replay_view(rq[x])[0]).split(),
                                             " ".join(replay_view(ra[x])[0]).split()) if a != b) for x in js)
            cq = replay_view(rq[js[0]])[1]
            return ntok, [cn for cn in TRAILER if cq.get(cn) != tr[k].get(cn)], None
        caught, missed = [], []
        for k, what, where, delta in plants:
            got = fails_of(p, k)
            ntok, moved, lost = moves(k)
            if lost:
                missed.append("member %d, %s: %s" % (k, what, lost))
            elif len(got) == 2 and all(d["verdict"] == "FAIL" and ("record line " in d["what"] or "trailer " in d["what"])
                                       for d in got) and (ntok > 0 or moved):
                caught.append("member %d, %s (line %d, %s made %s), which moves, as measured here, %d tokens of its "
                              "replays and %s: %s" % ((k, what) + said[k] + (ntok, ", ".join(moved) or "no counter",
                                                                               got[0]["what"][:100])))
            else:
                missed.append("member %d, %s: %s" % (k, what, [d["verdict"] + " " + d["what"][:80] for d in got]))
        calm = [x for x, d in p.items() if x[0] not in said and d["verdict"] != "IDENTICAL"]
        shutil.rmtree(q2, ignore_errors=True)
        shutil.rmtree(q2keep, ignore_errors=True)
        verdict(rc == 1 and len(plants) == 4 and not missed and not calm and l and limit_once(out) and
                all(v[0] == "IDENTICAL" for v in l.values()),
                "%s: q2, one pass count planted in the schedule, a member a line, %d lines: each FAILs both its member's "
                "particles by name, while that member's own rerun, where it has one, is IDENTICAL; every other particle "
                "IDENTICAL; the copy deleted (rc %d): %s%s" % (fmt, len(plants), rc, "; ".join(caught),
                                                              "; NOT: " + "; ".join(missed) if missed else ""))

    # q3 to q7, refused by the key or the form, at binary64
    c = clean["fp64"]
    u, v = (alone_of["fp64"] + [None, None])[:2]
    if u is None or v is None:
        verdict(False, "q3 to q7: two members that are not rerun, to plant in (%s)" % alone_of["fp64"])
        return

    def copy(tag, src=None):
        d = os.path.join(scratch, "producer_%s" % tag)
        shutil.rmtree(d, ignore_errors=True)
        shutil.copytree(src or pdirs["fp64"], d)
        return d
    d = copy("q3s")
    shutil.copyfile(pf(d, v)[1], pf(d, u)[1])
    refused("q3, member %d's schedule taken from member %d, another wave of the same ring" % (u, v), c["rec"], d, u,
            "corrector schedule", ["member", "state"])
    d = copy("q3a")
    shutil.copyfile(pf(d, v)[0], pf(d, u)[0])
    refused("q3, member %d's attempts record taken from member %d" % (u, v), c["rec"], d, u, "attempts record",
            ["member", "state"])
    d = copy("q5")
    member_run(waves, "fp128", c["steps"], c["sample"], u, ["--record-iters", pf(d, u)[1]])
    refused("q5, member %d's schedule recorded at binary128" % u, c["rec"], d, u, "corrector schedule", ["format"],
            exact=False)
    wa = write_waves(os.path.join(scratch, "cr_waves_active3.txt"), active=3)
    d = copy("q6")
    member_run(wa, "fp64", c["steps"], c["sample"], u, ["--record-iters", pf(d, u)[1]])
    refused("q6, member %d's schedule recorded under `active 3`, the problem otherwise the same" % u, c["rec"], d, u,
            "corrector schedule", ["active"])
    if recz is not None:
        wz2 = write_waves(os.path.join(scratch, "cr_waves_z2.txt"), name="cr_waves_z", zonal="0 %s %s %s" % (
            float.hex(0.25), float.hex(2.0 ** -6 * (1 + 2.0 ** -4)), float.hex(-2.0 ** -10)))
        d = copy("q4", pz)
        member_run(wz2, "fp64", c["steps"], c["sample"], u, ["--record-iters", pf(d, u)[1]])
        refused("q4, member %d's schedule recorded with the zonal line's J2 another, the problem otherwise the same"
                % u, recz, d, u, "corrector schedule", ["zonal"], problem=wz)
    d = copy("q7")
    shutil.copyfile(pf(d, u)[0], pf(d, u)[1])
    rc, out = pf_drive(c["rec"], "fp64", d)
    shutil.rmtree(d, ignore_errors=True)
    want = "RESULT: REFUSED: member %d's corrector schedule %s is not one: its first line is '# cft-rebound attempts v1'" \
           % (u, pf(d, u)[1])
    verdict(rc == 2 and want in out and "IDENTICAL" not in out,
            "q7, member %d's attempts record given as its schedule: REFUSED by its form, naming the file and its first "
            "line (rc %d); the copy deleted" % (u, rc))
    shutil.rmtree(pz, ignore_errors=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true", help="binary64, and short binary128 and FMA records")
    ap.add_argument("--driver", default=os.path.join(here, "compare_run.py"),
                    help="hold this copy of tools/compare_run.py instead (a mutant of it, to show the gate fails)")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, "ias15_cft" + exe)
    ref = os.path.join(args.build, "ias15_ref" + exe)
    progs = os.path.join(root, "programs", "out")
    py = sys.executable
    driver = os.path.abspath(args.driver)
    if driver != os.path.join(here, "compare_run.py"):
        print("check_compare_run: holding %s, not this tree's tools/compare_run.py" % driver)
    scratch = os.path.join(args.build, "check_compare_run")
    os.makedirs(scratch, exist_ok=True)
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= bool(ok)
        print("  %s %s" % ("ok  " if ok else "FAIL", text))
        return ok

    for b in (cft, ref):
        if not os.path.isfile(b):
            print("check_compare_run: %s is not built" % b)
            print("RESULT: FAIL")
            return 1
        print("check_compare_run: %s, built %s" % (b, time.strftime("%Y-%m-%d %H:%M:%S",
                                                                   time.localtime(os.path.getmtime(b)))))

    waves = write_waves(os.path.join(scratch, "cr_waves.txt"))

    def record(name, fmt, steps, sample, arith="rebound", engine="loop", problem=None):
        path = os.path.join(scratch, name)
        cmd = [cft, "--format", fmt, "--problem", problem or waves, "--dt", DT0, "--epsilon", "1e-9", "--steps", str(steps),
               "--sample", str(sample), "--max-iter", "60", "--arith", arith, "--engine", engine, "--quiet"]
        if engine == "program":
            cmd += ["--programs", progs]
        rc, out, err = run(cmd)
        if rc != 0:
            print("  FAIL the ensemble run %s (rc %d): %s" % (name, rc, err.strip()[-300:]))
            return None, None
        with open(path, "w", newline="\n") as f:
            f.write(out)
        return path, out

    def drive(rec, fmt, arith="rebound", extra=(), problem=None, via=None):
        cmd = [py, via or driver, rec, "--problem", problem or waves, "--format", fmt, "--arith", arith,
               "--exe", cft, "--programs", progs, "--strata", "3", "--seed", "1"] + list(extra)
        rc, out, err = run(cmd)
        if err.strip():
            out += "\n[stderr] " + err.strip()[-600:]
        return rc, out

    def trailers(text):
        """Each member's counters, read here from the record's own trailer."""
        return {k: parse_samples(text, k)[1] for k in range(E)}

    def sampled(out):
        m = re.search(r"^sample: seed (\d+), (\d+) strata by member index \(([^)]*)\): ([\d, ]+); with a rejected "
                      r"step: ([\d, ]+|none); (\d+) of (\d+) members$", out, re.M)
        if not m:
            return None
        strata = [tuple(int(x) for x in r.split("-")) for r in m.group(3).split(", ")]
        picks = [int(x) for x in m.group(4).split(", ")]
        rej = [] if m.group(5) == "none" else [int(x) for x in m.group(5).split(", ")]
        members = sorted(set(picks) | set(rej))
        return {"strata": strata, "picks": picks, "rejected": rej, "members": members, "n": int(m.group(6)),
                "E": int(m.group(7))}

    def member_lines(out):
        return {int(m.group(1)): (m.group(3), m.group(4)) for m in
                re.finditer(r"^member (\d+) \((\w+)\): (IDENTICAL|FAIL): (.*)$", out, re.M)}

    def hold_clean(label, rec_path, text, out, rc, fmt, steps, fdiff=None):
        """Case 1's assertions on one clean record. Returns the sample, and
        the rejecting members split by whether their FIRST attempt was."""
        tr = trailers(text)
        rej = sorted(k for k in range(E) if int(tr[k]["steps_rejected"]) > 0)
        # Which of them rejected their first attempt: each run alone for one
        # step, its trailer read here - not the driver's test of the dt file.
        first, later, unread = [], [], []
        for k in rej:
            rc1, o1, _ = run([cft, "--format", fmt, "--problem", waves, "--dt", DT0, "--epsilon", "1e-9", "--steps", "1",
                              "--sample", "1", "--max-iter", "60", "--member", str(k), "--quiet"])
            n1 = parse_samples(o1)[1].get("steps_rejected") if rc1 == 0 else None
            (unread if n1 is None else later if n1 == "0" else first).append(k)
        smp = sampled(out)
        lines = member_lines(out)
        verdict(rc == 0 and "\nRESULT: PASS: " in out, "%s: the driver PASSES (rc %d): %s" % (
            label, rc, (re.findall(r"^RESULT: .*$", out, re.M) or ["no RESULT line"])[-1]))
        if not verdict(smp is not None, "%s: the driver prints its sample" % label):
            print(out[-1500:])
            return None, first, later
        strata_ok = all(any(lo <= m <= hi for m in smp["members"]) for lo, hi in smp["strata"])
        verdict(first and later and not unread and len(rej) < E,
                "%s: members %s reject their first attempt, %s reject only after an accepted step, and %d reject nothing "
                "(this gate's reading of the record's trailers, and of each rejecting member run alone for one step)%s"
                % (label, first, later, E - len(rej), "; no one-step run read for %s" % unread if unread else ""))
        verdict(set(rej) <= set(smp["members"]) and smp["rejected"] == rej and strata_ok and
                len(smp["strata"]) == 3 and smp["n"] == len(smp["members"]) < E and smp["E"] == E,
                "%s: the sample %s holds every rejecting member %s and one member of each stratum %s, and is %d of %d"
                % (label, smp["members"], rej, smp["strata"], smp["n"], E))
        verdict(sorted(lines) == smp["members"] and all(v[0] == "IDENTICAL" for v in lines.values()),
                "%s: every sampled member is IDENTICAL: %s" % (label, ", ".join("%d %s" % (k, lines[k][0])
                                                                                 for k in sorted(lines))))
        bad = []
        for k in smp["members"]:
            t = tr[k]
            want = "steps %s, rejected %s, mean passes %s, max passes %s, cap hits %s" % (
                t["steps_done"], t["steps_rejected"], t["mean_pc_iterations"], t["max_pc_iterations"],
                t["iterations_max_exceeded"])
            if k not in lines or not lines[k][1].endswith(want) or t["steps_done"] != str(steps):
                bad.append("member %d: wanted '%s' for steps=%d" % (k, want, steps))
        verdict(not bad, "%s: each sampled member's counters are its trailer's, %d steps done each%s" % (
            label, steps, ": " + "; ".join(bad) if bad else ""))
        # closest approaches, recomputed from oracle.py's parse, in mpmath
        worst, missing, at_wrong, n = mpmath.mpf(0), [], [], 0
        least = {}                     # (member, body) -> its least distance to any active body
        for k in smp["members"]:
            _, _, smps, _ = oracle_parse(rec_path, system=k)
            for j in range(2, 2 + NTP):
                m = re.search(r"^  member %d: closest approach AT THE SAMPLES \((\d+) of them, not between them\).*\n"
                              r"(?:    .*\n)*?    tp%d \(body %d\): primary (\S+) at sample (\d+) \(step \d+\); "
                              r"moon (\S+) at sample (\d+) \(step \d+\)$" % (k, j - 2, j), out, re.M)
                if not m:
                    missing.append("member %d tp%d" % (k, j - 2))
                    continue
                for i, (val, at) in enumerate(((m.group(2), m.group(3)), (m.group(4), m.group(5)))):
                    dists = []
                    for s in smps:
                        st = s["state"]
                        dists.append(mpmath.sqrt(sum((F2m(st[6 * j + c]) - F2m(st[6 * i + c])) ** 2 for c in range(3))))
                    best = min(dists)
                    least[(k, j)] = min(least.get((k, j), best), best)
                    n += 1
                    worst = max(worst, abs(mpmath.mpf(val) - best) / best)
                    if dists[int(at)] != best:
                        at_wrong.append("member %d tp%d body %d" % (k, j - 2, i))
        verdict(not missing and not at_wrong and n == 2 * NTP * len(smp["members"]) and worst < 1e-6,
                "%s: %d closest approaches (each test particle to each active body, at the samples) agree with this "
                "gate's mpmath from tools/oracle.py's parse to %s relative, at the samples they name%s" % (
                    label, n, mpmath.nstr(worst, 3), "; " + ", ".join(missing + at_wrong) if missing or at_wrong else ""))
        # the particles (the lead's P3 x P4 seam): two a sampled member, its
        # least closest approach (by this gate's mpmath, ties to the lower
        # body) and one by the seed, each IDENTICAL to the record
        picks = {}
        for m in PICK.finditer(out):
            picks.setdefault(int(m.group(1)), []).append((int(m.group(3)), m.group(4), m.group(5), m.group(6)))
        pbad = []
        for k in smp["members"]:
            got = picks.get(k, [])
            near = min(range(2, 2 + NTP), key=lambda j: (least.get((k, j), mpmath.inf), j))
            kinds = sorted(g[1].split()[0] for g in got)
            if kinds != ["by", "its"] or any(g[2] != "IDENTICAL" for g in got):
                pbad.append("member %d: %s" % (k, got))
            elif next(g[0] for g in got if g[1].startswith("its")) != near:
                pbad.append("member %d: its least closest approach is body %d by mpmath, the driver chose %s"
                            % (k, near, [g[0] for g in got]))
            elif len({g[0] for g in got}) != 2:
                pbad.append("member %d: the same particle twice" % k)
        nwant = 2 * len(smp["members"])
        summary = "particles: %d sampled in %d member(s), each replayed alone from its member's attempts record " \
                  "and schedule as this tree reran it, and held to the record's own lines: %d identical" % (
                      nwant, len(smp["members"]), nwant)
        verdict(not pbad and summary in out and PARTICLES not in out,
                "%s: the particles, through P3's replay: %d sampled, two a member (each member's least closest approach "
                "at the samples, the body this gate's mpmath names, and one by the seed), each replayed alone and "
                "IDENTICAL to the record%s" % (label, nwant, ": " + "; ".join(pbad[:3]) if pbad else ""))
        verdict(out.count("\n" + SAMPLE_LIMIT + "\n") == 1 and out.find(SAMPLE_LIMIT) > out.find("\nsample: seed "),
                "%s: beside the sample, its limit, exactly: '%s'" % (label, SAMPLE_LIMIT))
        # The default path, where particles are sampled from this tree's
        # reruns, prints no limit of a producer's files (V3, 2026-11-18)
        verdict(limit_absent(out), "%s: without --producer-files, no line of the limit of a producer's files is "
                "printed, in any wording" % label)
        if fdiff:
            diffs = {int(m.group(1)): m for m in re.finditer(
                r"^  member (\d+): format difference: (\S+) against (\S+) on the same (\d+) accepted steps .*?"
                r"max \|dx\|/L ([^,\s]+), max \|dE/E\| (\S+) over (\d+) common samples \(tools/compare_formats.py\): "
                r"(.*)$", out, re.M)}
            bad, scoped = [], 0
            pair = tuple(sorted((fmt, fdiff), key=FORMATS.index))
            bound = 1e4 * UNIT_ROUNDOFF[pair[0]]
            for k in smp["members"]:
                m = diffs.get(k)
                if not m:
                    bad.append("member %d: no format difference line" % k)
                    continue
                dx = float(m.group(5))
                n = int(tr[k]["steps_rejected"])
                note = (FIRST_NOTE % fdiff if k in first else
                        LATER_NOTE % (n, "" if n == 1 else "s", fdiff) if k in later else "")
                scoped += m.group(8) == SCOPE + note
                if (m.group(2), m.group(3)) != pair or int(m.group(4)) != steps or not 0 < dx < bound or \
                        m.group(8) != SCOPE + note:
                    bad.append("member %d: %s" % (k, m.group(0).strip()))
            verdict(not bad, "%s: each sampled member's format difference, %s against %s on its own steps, is nonzero "
                    "and below 1e4 units of %s's round-off, %.2g (max |dx|/L %s)%s" % (
                        label, pair[0], pair[1], pair[0], bound,
                        ", ".join("%d: %s" % (k, diffs[k].group(5)) for k in sorted(diffs)),
                        ": " + "; ".join(bad) if bad else ""))
            verdict(not bad and scoped == len(smp["members"]) and set(first) | set(later) <= set(smp["members"]),
                    "%s: each of the %d says its scope beside its number, '%s', and the note its rejections call for: "
                    "first attempt rejected for %s, rejected only after an accepted step for %s, none for the rest"
                    % (label, scoped, SCOPE, first, later))
        return smp, first, later

    # ---------------------------------------------------------------- case 1
    fmt_cases = [("fp64", 40, 20, "fp128")] + ([] if args.quick else [("fp128", 20, 10, "fp256"), ("fp256", 20, 10, None)])
    smp64 = None
    text64 = None
    rec64 = None
    out64 = None
    clean = {}                   # fmt -> case 1's record there, for case 7
    for fmt, steps, sample, fdiff in fmt_cases:
        print("check_compare_run: case 1 at %s - the driver on an ensemble of %d waves, %d steps%s" % (
            fmt, E, steps, ", with --format-diff " + fdiff if fdiff else ""))
        rec, text = record("waves_%s.rec" % fmt, fmt, steps, sample)
        if rec is None:
            allok = False
            continue
        rc, out = drive(rec, fmt, extra=["--format-diff", fdiff] if fdiff else [])
        smp, first, later = hold_clean(fmt, rec, text, out, rc, fmt, steps, fdiff)
        if smp is None:
            allok = False
            continue
        clean[fmt] = {"rec": rec, "text": text, "steps": steps, "sample": sample, "smp": smp, "first": first,
                      "later": later}
        if fmt == "fp64":
            smp64, text64, rec64, out64 = smp, text, rec, out
            # One member of each kind replayed from its own dt file at its own
            # format, against the record: what each note above says of it.
            # dt_next is set aside where said, and only there: a sample line
            # printed just before a rejected attempt has that attempt's dt as
            # its dt_next, and the replay never tries that attempt.
            skip = value_index(text, "dt_next")

            def replayed(k):
                dt = os.path.join(scratch, "member%d.dt" % k)
                base = [cft, "--format", "fp64", "--problem", waves, "--steps", str(steps), "--sample", str(sample),
                        "--max-iter", "60", "--member", str(k), "--quiet"]
                rc1, a, _ = run(base + ["--dt", DT0, "--epsilon", "1e-9", "--dt-out", dt])
                rc2, b, _ = run(base + ["--dt-file", dt, "--epsilon", "0"])
                es = parse_samples(text, k)[0]
                if rc1 or rc2:
                    return es, None, None
                return es, parse_samples(a)[0], parse_samples(b)[0]

            def state(samples):
                return [(k_, st, [v for i, v in enumerate(vals) if i != skip]) for k_, st, vals in samples or []]

            quiet = [k for k in smp["members"] if k not in smp["rejected"]]
            if verdict(quiet, "fp64: a sampled member that rejected nothing exists (%s)" % quiet):
                k = quiet[-1]
                es, ra, rb = replayed(k)
                verdict(es and es == ra == rb,
                        "fp64: member %d, which rejected nothing, replayed from its dt file at fp64 is the record's member "
                        "bit for bit (%d samples): so its nonzero difference above is another format" % (k, len(es)))
            if first:
                k = first[0]
                es, ra, rb = replayed(k)
                agree = sum(1 for x, y in zip(state(es), state(rb)) if x == y)
                verdict(es and es == ra and rb and len(rb) == len(es) and state(rb) != state(es),
                        "fp64: member %d, whose first attempt was rejected, replayed from its dt file at fp64: its state "
                        "is the record's at %d of %d samples, so that rejection leaves state, as its note says"
                        % (k, agree, len(es)))
            if later:
                k = later[-1]
                es, ra, rb = replayed(k)
                ndt = sum(1 for x, y in zip(es, rb or []) if x[2][skip] != y[2][skip])
                verdict(es and es == ra and rb and state(rb) == state(es),
                        "fp64: member %d, whose %s rejected attempts all followed an accepted step, replayed from its dt "
                        "file at fp64: every value but dt_next is the record's at all %d samples (dt_next differs at %d), "
                        "so they leave no state bit, as its note says" % (
                            k, parse_samples(text, k)[1].get("steps_rejected"), len(es), ndt))
    if smp64 is None:
        print("RESULT: FAIL")
        return 1

    # ---------------------------------------------------------------- case 3
    print("check_compare_run: case 3 - the sample, seeded and stratified (--plan)")
    plans = {}
    for seed in range(1, 11):
        rc, out = drive(rec64, "fp64", extra=["--plan", "--seed", str(seed)])
        plans[seed] = sampled(out) if rc == 0 else None
    good = [s for s in plans if plans[s] and set(smp64["rejected"]) <= set(plans[s]["members"]) and
            all(any(lo <= m <= hi for m in plans[s]["picks"]) for lo, hi in plans[s]["strata"]) and
            len(plans[s]["picks"]) == 3]
    verdict(len(good) == 10, "every seed's sample holds every rejecting member and one member a stratum (%d of 10 seeds)"
            % len(good))
    verdict(len(set(tuple(plans[s]["picks"]) for s in good)) > 1,
            "the seed chooses: %d different picks over 10 seeds" % len(set(tuple(plans[s]["picks"]) for s in good)))
    verdict(plans.get(1) and plans[1]["members"] == smp64["members"], "seed 1's plan is case 1's sample (%s)"
            % (plans[1]["members"] if plans.get(1) else None))

    # ---------------------------------------------------------------- case 4
    print("check_compare_run: case 4 - the same problem at binary128 and in the FMA form, held as recorded")
    rec128, text128 = record("waves_fp128_short.rec", "fp128", 4, 2)
    recfma, _ = record("waves_fma_short.rec", "fp64", 4, 2, arith="fma")
    smp128 = None
    for label, rec, fmt, arith in (("fp128", rec128, "fp128", "rebound"), ("fp64 FMA form", recfma, "fp64", "fma")):
        rc, out = drive(rec, fmt, arith) if rec else (-1, "no record")
        lines = member_lines(out)
        verdict(rc == 0 and "\nRESULT: PASS: " in out and lines and all(v[0] == "IDENTICAL" for v in lines.values()),
                "%s: held at its own format and arithmetic, PASS, %d members IDENTICAL" % (label, len(lines)))
        if fmt == "fp128":
            smp128 = sampled(out)

    # ---------------------------------------------------------------- case 5
    print("check_compare_run: case 5 - one system's records: a --member run's, and a single problem's")
    kep = os.path.join(root, "data", "problems", "kepler.txt")
    k6 = 3
    rc, t6, err = run([cft, "--format", "fp64", "--problem", waves, "--dt", DT0, "--epsilon", "1e-9", "--steps", "40",
                       "--sample", "20", "--max-iter", "60", "--member", str(k6), "--quiet"])
    rec6 = os.path.join(scratch, "waves_member%d.rec" % k6)
    with open(rec6, "w", newline="\n") as f:
        f.write(t6)
    rc, out = drive(rec6, "fp64")
    lines = member_lines(out)
    t = parse_samples(t6)[1]
    want = "steps %s, rejected %s, mean passes %s, max passes %s, cap hits %s" % (
        t.get("steps_done"), t.get("steps_rejected"), t.get("mean_pc_iterations"), t.get("max_pc_iterations"),
        t.get("iterations_max_exceeded"))
    verdict(rc == 0 and "\nRESULT: PASS: " in out and sorted(lines) == [k6] and lines[k6][0] == "IDENTICAL" and
            lines[k6][1].endswith(want) and "one system of N=%d bodies, member %d of its problem" % (2 + NTP, k6) in out
            and "  member %d: closest approach AT THE SAMPLES (3 of them" % k6 in out and
            "\nparticles: 2 sampled in 1 member(s), " in out and "held to the record's own lines: 2 identical\n" in out
            and len([m for m in PICK.finditer(out) if m.group(5) == "IDENTICAL"]) == 2 and limit_absent(out),
            "a --member %d record of the waves, held as one system: PASS, member %d IDENTICAL with its trailer's "
            "counters (%s), its closest approaches, and its two particles replayed alone, IDENTICAL; and, without "
            "--producer-files, no line of the limit of a producer's files (rc %d)" % (k6, k6, want, rc))
    rck, tk, _ = run([cft, "--format", "fp64", "--problem", kep, "--dt", "0.05", "--epsilon", "1e-9", "--steps", "20",
                      "--sample", "10", "--max-iter", "60", "--quiet"])
    reck = os.path.join(scratch, "kepler.rec")
    with open(reck, "w", newline="\n") as f:
        f.write(tk)
    rc, out = drive(reck, "fp64", problem=kep)
    lines = member_lines(out)
    verdict(rck == 0 and rc == 0 and "\nRESULT: PASS: " in out and sorted(lines) == [0] and
            lines[0][0] == "IDENTICAL" and "member 0 (kepler): IDENTICAL: 3 samples" in out and
            "  member 0: no test particles (no `# active` line: every body is active)" in out and limit_absent(out),
            "a single problem's record (kepler, every body active): PASS, IDENTICAL, and no test particles said; and, "
            "without --producer-files, no line of the limit of a producer's files (rc %d)" % rc)

    if not args.quick:
        print("check_compare_run: case 6 - the FMA form on the program engine, as recorded and on the loop engine")
        recprog, textp = record("waves_fma_program.rec", "fp64", 20, 10, arith="fma", engine="program")
        if recprog:
            for extra, label in (([], "the program engine"), (["--rerun-engine", "loop"], "the loop engine")):
                rc, out = drive(recprog, "fp64", "fma", extra)
                lines = member_lines(out)
                verdict(rc == 0 and "\nRESULT: PASS: " in out and lines and
                        all(v[0] == "IDENTICAL" for v in lines.values()) and
                        (not extra or "rerun: on the loop engine, where the record ran program" in out),
                        "the FMA form's program-engine record, rerun on %s: PASS, %d members IDENTICAL"
                        % (label, len(lines)))
        else:
            allok = False

    # ---------------------------------------------------------------- case 7
    case7(clean, [f for f, *_ in fmt_cases], waves, scratch, cft, progs, args.quick, record, drive, run, verdict,
          trailers, member_lines)

    # ---------------------------------------------------------- the controls
    print("check_compare_run: the negative controls, each of which must FAIL or be refused, by name")
    text = text64
    rows = text.splitlines()
    nlead, per = columns(text)
    lines64 = [i for i, l in enumerate(rows) if not l.startswith("#")]

    def write(name, lines):
        path = os.path.join(scratch, name)
        with open(path, "w", newline="\n") as f:
            f.write("\n".join(lines) + "\n")
        return path

    # a. one byte of a sampled member's last sample line
    k = smp64["picks"][-1]
    i = max(i for i in lines64 if rows[i].startswith("system %d sample " % k))
    toks = rows[i].split()
    b = NTP + 1                                        # the last test particle, body 4
    t = nlead + b * len(per) + per.index("x")
    tok = toks[t]
    mant = tok.split("p")[0]
    toks[t] = mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + tok[len(mant):]
    changed = rows[:i] + [" ".join(toks)] + rows[i + 1:]
    rc, out = drive(write("control_byte.rec", changed), "fp64")
    lines = member_lines(out)
    others = [m for m in lines if m != k]
    verdict(rc == 1 and "\nRESULT: FAIL: " in out and lines.get(k, ("",))[0] == "FAIL" and
            ("first differing line: record line %d (system %d sample " % (i + 1, k)) in lines[k][1] and
            ("x of body %d (tp%d): record %s, here %s" % (b, NTP - 1, toks[t], tok)) in lines[k][1] and
            others and all(lines[m][0] == "IDENTICAL" for m in others),
            "a. one byte of member %d's line %d (tp%d's x, %s for %s): FAIL naming member %d, line %d and the column; "
            "the other %d sampled IDENTICAL (rc %d)%s" % (k, i + 1, NTP - 1, toks[t], tok, k, i + 1, len(others), rc,
                                                          "" if rc == 1 else ": " + out[-800:]))

    # a2. the same in one system's record: the --member 3 run of case 5
    rows6 = t6.splitlines()
    nlead6, per6 = columns(t6)
    i = max(i for i, l in enumerate(rows6) if l.startswith("sample "))
    toks = rows6[i].split()
    t = nlead6 + b * len(per6) + per6.index("x")
    tok = toks[t]
    mant = tok.split("p")[0]
    toks[t] = mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + tok[len(mant):]
    rc, out = drive(write("control_byte_member.rec", rows6[:i] + [" ".join(toks)] + rows6[i + 1:]), "fp64")
    lines = member_lines(out)
    verdict(rc == 1 and "\nRESULT: FAIL: " in out and lines.get(k6, ("",))[0] == "FAIL" and
            ("first differing line: record line %d (sample " % (i + 1)) in lines[k6][1] and
            ("x of body %d (tp%d): record %s, here %s" % (b, NTP - 1, toks[t], tok)) in lines[k6][1],
            "a2. one byte of the --member %d record's line %d: FAIL naming member %d (w%d), line %d and the column (rc %d)"
            % (k6, i + 1, k6, k6, i + 1, rc))

    # The next two are V4's (2026-11-12): a comparison of binary64 values,
    # not of bits, passed every control above. Every one of them moves a
    # binary64 bit; these two move none.
    # a3. the sign of a zero: a sampled member's sample-0 dt_last made -0
    k = smp64["members"][0]
    i = min(i for i in lines64 if rows[i].startswith("system %d sample 0 " % k))
    toks = rows[i].split()
    t = scalar_tok(text, "dt_last")
    tok = toks[t]
    toks[t] = "-" + tok
    rc, out = drive(write("control_signed_zero.rec", rows[:i] + [" ".join(toks)] + rows[i + 1:]), "fp64")
    lines = member_lines(out)
    others = [m for m in lines if m != k]
    verdict(hex_to_fraction(tok) == (0, False) and rc == 1 and "\nRESULT: FAIL: " in out and
            lines.get(k, ("",))[0] == "FAIL" and
            ("first differing line: record line %d (system %d sample 0, step 0): dt_last: record -%s, here %s"
             % (i + 1, k, tok, tok)) in lines[k][1] and others and all(lines[m][0] == "IDENTICAL" for m in others),
            "a3. the sign of a zero, member %d's dt_last %s made -%s on line %d: FAIL naming member %d, line %d and "
            "dt_last; the other %d sampled IDENTICAL (rc %d)" % (k, tok, tok, i + 1, k, i + 1, len(others), rc))

    # a4. one hex digit past binary64, in case 4's binary128 record: tp0's x
    # in a sampled member's last line, its last hex digit, which leaves the
    # value's binary64 rounding where it was
    if smp128 and text128:
        rows128 = text128.splitlines()
        nlead128, per128 = columns(text128)
        k = smp128["picks"][-1]
        i = max(i for i, l in enumerate(rows128) if l.startswith("system %d sample " % k))
        toks = rows128[i].split()
        t = nlead128 + 2 * len(per128) + per128.index("x")
        tok = toks[t]
        mant = tok.split("p")[0]
        toks[t] = mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + tok[len(mant):]
        digits = len(mant.split(".")[1]) if "." in mant else 0
        past64 = digits > 13 and float(hex_to_fraction(toks[t])[0]) == float(hex_to_fraction(tok)[0])
        rc, out = drive(write("control_past64.rec", rows128[:i] + [" ".join(toks)] + rows128[i + 1:]), "fp128")
        lines = member_lines(out)
        others = [m for m in lines if m != k]
        verdict(past64 and rc == 1 and "\nRESULT: FAIL: " in out and lines.get(k, ("",))[0] == "FAIL" and
                ("first differing line: record line %d (system %d sample " % (i + 1, k)) in lines[k][1] and
                ("x of body 2 (tp0): record %s, here %s" % (toks[t], tok)) in lines[k][1] and
                others and all(lines[m][0] == "IDENTICAL" for m in others),
                "a4. one hex digit past binary64 (digit %d of %d, the same binary64 value): member %d's tp0 x %s made %s "
                "on line %d of the binary128 record: FAIL naming member %d, line %d and the column; the other %d "
                "sampled IDENTICAL (rc %d)" % (digits, digits, k, tok, toks[t], i + 1, k, i + 1, len(others), rc))
    else:
        verdict(False, "a4. case 4's binary128 record and its sample, to change a digit past binary64 in")

    # a5. Every column the record names is compared (V4, 2026-11-12: a driver
    # that skipped the exact-time pair passed every control above). From the
    # record's own `# columns:` line, one planted change a column (its last
    # hex digit) in one sampled member's last sample line. The driver names
    # each member's first differing line, so a record carries one column a
    # sampled member, and a few records carry them all.
    Nb = int(re.search(r" N=(\d+)", next(l for l in rows if l.startswith("# program="))).group(1))
    body_names = [l.split()[3] for l in rows if l.startswith("# body ")][:Nb]
    plan = column_plan(text, Nb, body_names)
    members = smp64["members"]
    held, missed, nrec = [], [], 0
    for start in range(0, len(plan), len(members)):
        nrec += 1
        doctored, want = list(rows), {}
        for (name, t), k in zip(plan[start:start + len(members)], members):
            i = max(i for i in lines64 if rows[i].startswith("system %d sample " % k))
            toks = doctored[i].split()
            old = toks[t]
            toks[t] = flip_last_digit(old)
            doctored[i] = " ".join(toks)
            want[k] = (name, i + 1, toks[t], old)
        rc, out = drive(write("control_columns_%d.rec" % nrec, doctored), "fp64")
        lines = member_lines(out)
        if rc != 1:
            missed.append("record %d: rc %d" % (nrec, rc))
        for k in members:
            got = lines.get(k, ("", ""))
            if k not in want:
                if got[0] != "IDENTICAL":
                    missed.append("member %d, with nothing planted, is %s" % (k, got[0] or "not reported"))
                continue
            name, lineno, new, old = want[k]
            if got[0] == "FAIL" and ("first differing line: record line %d (system %d sample " % (lineno, k)) in got[1] \
                    and ("): %s: record %s, here %s" % (name, new, old)) in got[1]:
                held.append(name)
            else:
                missed.append("%s in member %d: %s %s" % (name, k, got[0] or "not reported", got[1][:160]))
    verdict(not missed and len(held) == len(plan),
            "a5. every column the record names, %d of them, each one planted change in one sampled member's last line, "
            "%d to a record over %d records: each FAILs naming its column (%s)%s"
            % (len(plan), len(members), nrec, ", ".join(held), "; NOT: " + "; ".join(missed) if missed else ""))

    # b. one byte of a sampled member's trailer
    k = smp64["picks"][0]
    i = next(i for i, l in enumerate(rows) if l.startswith("# system %d steps_done=" % k))
    old = re.search(r"mean_pc_iterations=(\S+)", rows[i]).group(1)
    new = old[:-1] + str((int(old[-1]) + 1) % 10)
    rc, out = drive(write("control_trailer.rec", rows[:i] + [rows[i].replace("mean_pc_iterations=" + old,
                                                                            "mean_pc_iterations=" + new)] + rows[i + 1:]),
                    "fp64")
    lines = member_lines(out)
    verdict(rc == 1 and lines.get(k, ("",))[0] == "FAIL" and
            ("trailer mean_pc_iterations is %s in the record (line %d), and %s here" % (new, i + 1, old)) in lines[k][1],
            "b. member %d's trailer mean_pc_iterations %s made %s (line %d): FAIL naming member %d and the counter (rc %d)"
            % (k, old, new, i + 1, k, rc))

    # c. every line of a member that is not sampled
    u = max(m for m in range(E) if m not in smp64["members"])
    drop = set()
    for i, l in enumerate(rows):
        if l.startswith("# system %d " % u) and "steps_done=" not in l:
            drop.update(range(i, i + 1 + 2 + NTP))
        elif l.startswith("system %d sample " % u) or l.startswith("# system %d steps_done=" % u):
            drop.add(i)
    rc, out = drive(write("control_missing.rec", [l for i, l in enumerate(rows) if i not in drop]), "fp64")
    lines = member_lines(out)
    sampled_ok = [m for m in smp64["members"] if lines.get(m, ("",))[0] == "IDENTICAL"]
    verdict(rc == 1 and lines.get(u, ("",))[0] == "FAIL" and
            lines[u][1].startswith("missing from the record: no `# system %d` line; 0 of %d `# body` lines; 0 of "
                                   % (u, 2 + NTP)) and sampled_ok == smp64["members"],
            "c. member %d's %d lines deleted, a member not sampled: FAIL naming it missing; the %d sampled IDENTICAL "
            "(rc %d)" % (u, len(drop), len(sampled_ok), rc))

    # d. the record cut short before its last sample
    last = max(int(l.split()[3]) for l in rows if l.startswith("system "))
    cut = next(i for i, l in enumerate(rows) if l.startswith("system 0 sample %d " % last))
    rc, out = drive(write("control_cut.rec", rows[:cut]), "fp64")
    lines = member_lines(out)
    verdict(rc == 1 and sorted(lines) == list(range(E)) and
            all(v[0] == "FAIL" and v[1].startswith("incomplete in the record:") for v in lines.values()) and
            "record: FAIL: no `# steps_done=` trailer: the run did not finish" in out,
            "d. the record cut short at line %d: FAIL, all %d members incomplete, and no trailer (rc %d)"
            % (cut + 1, len(lines), rc))

    # e. another format, another arithmetic
    for label, rec, fmt, arith, needle in (
            ("a binary128 record held as binary64", rec128, "fp64", "rebound",
             "the record is format=fp128, and this comparison holds format=fp64 (--format)"),
            ("an FMA-form record held as REBOUND's arithmetic", recfma, "fp64", "rebound",
             "the record is arith=fma, and this comparison holds arith=rebound (--arith)")):
        rc, out = drive(rec, fmt, arith) if rec else (-1, "")
        verdict(rc == 2 and ("RESULT: REFUSED: " + needle) in out and "IDENTICAL" not in out,
                "e. %s: REFUSED by name (rc %d)" % (label, rc))

    # f. REBOUND's own record, and a file that is no record
    rc, rtext, _ = run([ref, "--problem", kep, "--dt", "0.05", "--steps", "4", "--sample", "2"])
    refrec = write("control_ias15_ref.rec", rtext.splitlines())
    rc, out = drive(refrec, "fp64", problem=kep)
    verdict(rc == 2 and rtext.startswith("# cft-rebound ias15 record v1\n") and
            "names no arith, cs, engine, pc_tol_shift, max_iter in its `# program=` line (it is program=ref" in out,
            "f. REBOUND's own record (ias15_ref, the same first line): REFUSED, naming the fields it lacks (rc %d)" % rc)
    rc, out = drive(waves, "fp64")
    verdict(rc == 2 and "is not a cft-rebound ias15 record: its first line is" in out,
            "f. a problem file given as the record: REFUSED, not a record (rc %d)" % rc)

    # g. a different problem file under the record's name
    for label, path, needle in (
            ("a test particle's mass 0 made 2^-40", write_waves(os.path.join(scratch, "cr_waves_tpmass.txt"),
                                                                tp_mass=(3, 1, 2.0 ** -40)),
             "system 3 (w3) body 3 (tp1) has m 0x1.0000000000000p-40, and the record m=0x0p+0"),
            ("the active line dropped", write_waves(os.path.join(scratch, "cr_waves_noactive.txt"), active=False),
             "has no active line, and the record says active 2"),
            ("one initial coordinate moved one ulp", write_waves(os.path.join(scratch, "cr_waves_ulp.txt"), ulp=(5, 0)),
             "is not the problem the record ran, though both are named cr_waves: system 5 (w5) body 2 (tp0) starts at x"),
            ("the name changed", write_waves(os.path.join(scratch, "cr_waves_name.txt"), name="cr_waves2"),
             "is named cr_waves2, and the record ran problem=cr_waves")):
        rc, out = drive(rec64, "fp64", problem=path)
        verdict(rc == 2 and needle in out and "IDENTICAL" not in out,
                "g. a different problem file under the record's name, %s: REFUSED by name (rc %d)%s"
                % (label, rc, "" if rc == 2 and needle in out else ": " + out[-600:]))

    # h. the record's `# active` line with another pair count than this
    # tree's rows: refused at the first rerun, naming the field that differs
    i = next(i for i, l in enumerate(rows) if l.startswith("# active "))
    m = re.match(r"^(# active \d+ of \d+ bodies a system: )(\d+)( pairs .*)$", rows[i])
    pairs = int(m.group(2))
    rc, out = drive(write("control_pairs.rec", rows[:i] + [m.group(1) + str(pairs + 1) + m.group(3)] + rows[i + 1:]),
                    "fp64")
    verdict(rc == 2 and ("the `# active` line's pairs a force evaluation is %d here, and %d in the record"
                         % (pairs, pairs + 1)) in out and "IDENTICAL" not in out,
            "h. the record's `# active` line saying %d pairs a force evaluation where this tree's rows make %d: REFUSED, "
            "naming the field (rc %d)" % (pairs + 1, pairs, rc))

    # p. The particle path is held to the record, not to the rerun (the lead's
    # P3 x P4 seam). From case 1's fp64 run, member k's least-closest-
    # approach particle j. In a copy of the record, one hex digit of j's x on
    # k's last sample line: the member FAILs, and so does j's replay, naming
    # the line and the column, while k's other sampled particle, whose
    # sub-system has no body j, stays IDENTICAL. Then k's primary's x: the
    # primary is in every sub-system, so both of k's replays FAIL.
    picks64 = {}
    for m in PICK.finditer(out64 or ""):
        picks64.setdefault(int(m.group(1)), []).append((int(m.group(3)), m.group(4)))
    kp = smp64["members"][0]
    jp = next((j for j, why in picks64.get(kp, []) if why.startswith("its")), None)
    if jp is None:
        verdict(False, "p. case 1's fp64 run named no particle of member %d to plant a change in" % kp)
    else:
        ip = max(i for i in lines64 if rows[i].startswith("system %d sample " % kp))
        for body, what in ((jp, "its sampled particle"), (0, "its primary")):
            toks = rows[ip].split()
            t = nlead + body * len(per) + per.index("x")
            toks[t] = flip_last_digit(toks[t])
            rc, out = drive(write("control_particle_body%d.rec" % body, rows[:ip] + [" ".join(toks)] + rows[ip + 1:]),
                            "fp64")
            lines = member_lines(out)
            pm = [m for m in PICK.finditer(out) if int(m.group(1)) == kp]
            hit = [m for m in pm if m.group(5) == "FAIL" and ("record line %d (system %d sample " % (ip + 1, kp))
                   in m.group(6) and ("x of body %d " % body) in m.group(6)]
            calm = [m for m in pm if m.group(5) == "IDENTICAL"]
            want = 1 if body == jp else 2
            verdict(rc == 1 and lines.get(kp, ("",))[0] == "FAIL" and len(pm) == 2 and len(hit) == want
                    and len(calm) == 2 - want,
                    "p. one hex digit of member %d's %s's x (body %d) on line %d: the member FAILs, and %d of its 2 "
                    "sampled particles' replays FAIL naming that line and column (%s IDENTICAL): the replay is held "
                    "to the record, not to the rerun (rc %d)"
                    % (kp, what, body, ip + 1, len(hit), "the other" if want == 1 else "none", rc))

    # p2. Every column the particle path compares is compared (VL, 2026-11-12:
    # control p plants two of them, and a particle_compare that checked only
    # those two passed every control). From the record's own `# columns:`
    # line, each sampled member's least-closest-approach particle j gets its
    # plan (particle_plan): each field before the bodies, every field of each
    # massive body and of j, and the exact-time pair. One planted change a
    # member a record, on one of the member's sample lines after sample 0 in
    # turn - VL again, 2026-11-12: plants on the last line alone passed a
    # particle_compare that read the last line alone; sample 0 is the
    # problem's initial state, which the problem check holds and refuses on.
    # Each must FAIL j's replay, naming its line and column. Then three
    # trailer counters, one a member in one record (mean_pc_iterations,
    # max_pc_iterations, iterations_max_exceeded), each of which must FAIL
    # its member's replay; steps_done and steps_rejected are held by the
    # record's structure and the member path, and are not planted here.
    body_names = [l.split()[3] for l in rows if l.startswith("# body ")][:Nb]
    nact64 = int(re.search(r"^# active (\d+) of ", text, re.M).group(1))
    plans = {}
    for k in smp64["members"]:
        j = next((j for j, why in picks64.get(k, []) if why.startswith("its")), None)
        if j is not None:
            plans[k] = (j, particle_plan(text, Nb, body_names, nact64, j))
    held, missed, nrec, ncols = [], [], 0, 0
    if len(plans) != len(smp64["members"]):
        missed.append("case 1's fp64 run named no least-closest particle for some sampled member: %s" % sorted(plans))
    width = len(next(iter(plans.values()))[1]) if plans else 0
    order = sorted(plans)
    lines_of = {k: sorted(i for i in lines64 if rows[i].startswith("system %d sample " % k)) for k in order}
    positions = set()
    for start in range(0, width, len(order)):
        nrec += 1
        doctored, want = list(rows), {}
        for c, k in zip(range(start, min(width, start + len(order))), order):
            j, plan = plans[k]
            name, t = plan[c]
            later = lines_of[k][1:]
            i = later[c % len(later)]
            positions.add(lines_of[k].index(i))
            toks = doctored[i].split()
            old = toks[t]
            toks[t] = flip_last_digit(old)
            doctored[i] = " ".join(toks)
            want[k] = (j, name, i + 1, toks[t], old)
            ncols += 1
        rc, out = drive(write("control_particle_columns_%d.rec" % nrec, doctored), "fp64")
        found = {(int(m.group(1)), int(m.group(3))): m for m in PICK.finditer(out)
                 if m.group(4).startswith("its")}
        for k, (j, name, lineno, new, old) in want.items():
            m = found.get((k, j))
            if m and m.group(5) == "FAIL" and ("record line %d (system %d sample " % (lineno, k)) in m.group(6) \
                    and ("): %s: record %s, replay %s" % (name, new, old)) in m.group(6):
                held.append(name)
            else:
                missed.append("%s in member %d's particle body %d: %s" % (
                    name, k, j, ("%s %s" % (m.group(5), m.group(6)[:140])) if m else "no particle line"))
    # three trailer counters, one a sampled member, in one record
    planted = ("mean_pc_iterations", "max_pc_iterations", "iterations_max_exceeded")
    doctored, want_t = list(rows), {}
    for counter, k in zip(planted, order):
        i = next(i for i, l in enumerate(rows) if l.startswith("# system %d steps_done=" % k))
        doctored[i] = re.sub(r"\b%s=(\S+)" % counter,
                             lambda mm: "%s=%s" % (counter, mm.group(1) + "1" if "." in mm.group(1)
                                                   else str(int(mm.group(1)) + 1)), doctored[i], count=1)
        want_t[k] = (plans[k][0], counter)
    if want_t:
        rc, out = drive(write("control_particle_trailer.rec", doctored), "fp64")
        for k, (j, counter) in want_t.items():
            m = next((m for m in PICK.finditer(out) if int(m.group(1)) == k and int(m.group(3)) == j), None)
            if m and m.group(5) == "FAIL" and ("trailer %s is " % counter) in m.group(6):
                held.append("the trailer's %s" % counter)
            else:
                missed.append("%s in member %d's trailer: %s" % (
                    counter, k, ("%s %s" % (m.group(5), m.group(6)[:140])) if m else "no particle line"))
    if len(positions) < 2:
        missed.append("every plant fell on one sample line (%s)" % sorted(positions))
    verdict(not missed and len(held) == ncols + len(planted),
            "p2. every column the particle path compares, %d of them a particle (each field before the bodies, every "
            "field of each massive body and of the particle, and the exact-time pair), on sample lines %s of each "
            "member in turn, and %d trailer counters, one planted change a member over %d records and a trailer "
            "record: each FAILs its member's least-closest particle's replay, naming its line and column (%d held)%s"
            % (width, ", ".join(str(p) for p in sorted(positions)), len(planted), nrec, len(held),
               "; NOT: " + "; ".join(missed[:4]) if missed else ""))

    # i. P3's tools in the tree and not callable: a copy of the driver, with
    # what it imports and runs, beside a stand-in for each of P3's two files
    # in turn, and beside neither, on case 5's kepler record. The stand-in is
    # an empty file, with no extract() to call.
    p3 = os.path.join(scratch, "p3_present", "tools")
    shutil.rmtree(os.path.dirname(p3), ignore_errors=True)
    os.makedirs(p3)
    for f in ("compare_run.py", "check_ensemble.py", "check_records.py", "check_equivalence.py", "compare_formats.py",
              "oracle.py"):
        shutil.copy2(os.path.join(here, f), p3)
    for present in ((), ("extract_particle.py",), ("check_replay.py",)):
        for f in ("extract_particle.py", "check_replay.py"):
            path = os.path.join(p3, f)
            if f in present:
                with open(path, "w", newline="\n") as fh:
                    fh.write("# check_compare_run's stand-in for P3's file: only that it is here is tested\n")
            elif os.path.exists(path):
                os.remove(path)
        rc, out = drive(reck, "fp64", problem=kep, via=os.path.join(p3, "compare_run.py"))
        lines = member_lines(out)
        result = (re.findall(r"^RESULT: .*$", out, re.M) or ["no RESULT line"])[-1]
        if not present:
            verdict(rc == 0 and result.startswith("RESULT: PASS: ") and out.count("\n" + PARTICLES + "\n") == 1 and
                    lines.get(0, ("",))[0] == "IDENTICAL",
                    "i. the driver's copy beside neither of P3's files: '%s', PASS, rc %d" % (PARTICLES, rc))
        else:
            why = ("tools/extract_particle.py has no extract() and Refused to call" if present[0] == "extract_particle.py"
                   else "tools/check_replay.py is in this tree without tools/extract_particle.py")
            verdict(rc == 1 and (BROKEN + why) in out and PARTICLES not in out and
                    result.startswith("RESULT: FAIL: ") and "P3's replay is in this tree and cannot be called" in result
                    and lines.get(0, ("",))[0] == "IDENTICAL",
                    "i. the driver's copy beside a stand-in tools/%s: P3's replay cannot be called (%s), RESULT: FAIL "
                    "and rc %d, though member 0 is IDENTICAL" % (present[0], why, rc))
    shutil.rmtree(os.path.dirname(p3), ignore_errors=True)

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
