#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# A run's record against this tree (ROUND4.md, P4): a seeded sample of its
# members rerun here, each reported identical or not - bit for bit, by name
# and first differing line - with its physics counters beside the bits.
#
#   python tools/compare_run.py RECORD --problem FILE --format F --arith A
#          [--exe build/ias15_cft] [--programs programs/out]
#          [--seed 1] [--strata 8] [--format-diff F2] [--plan]
#          [--rerun-engine loop|program] [--keep DIR]
#          [--producer-files DIR [--particle-strata P]]
#
# RECORD is the run's record: an ensemble's, or one system's (a --member k
# run's, or a single problem's). Today it is this tree's own run; the point
# of the driver is a GPU's. So it READS ONLY THE RECORD FORMAT - the lines
# src/ias15_cft.c's main() writes, from `# cft-rebound ias15 [ensemble]
# record v1` to the trailer - and never ias15_cft's internals. An engine
# that writes the same records is held without changing a line here.
#
# REFUSED (exit 2), by name and before anything is compared - each one a
# comparison that would be meaningless rather than a difference:
#   - a file that is not such a record, or whose header lacks a field the
#     rerun is built from (REBOUND's own ias15_ref writes the same first
#     line, with program=ref and no arith=), or one run from a dt file;
#   - a record of another format or arithmetic than --format and --arith,
#     which name the run this comparison holds. A genuine record of
#     another format would be rerun at its own format and agree, so the
#     caller says which run it means, as a test says what it expects;
#   - a problem file that is not the problem the record ran. The record
#     names its problem only by the file's `name` line, so the file is an
#     argument and is held to the record first: its name, G, its systems,
#     every body's name and mass, its `active` line, and every system's
#     sample 0, which is that system's initial state. A test particle given
#     a tiny mass, or the active line dropped from massless ones, changes
#     no bit (testparticle_type 0, ROUND4.md P0) - so only this check sees
#     it;
#   - a rerun whose own header is not the record's run: format, arith, cs,
#     engine, steps, sample, N, pc_tol_shift, max_iter, problem, dt0,
#     epsilon, G, pc_tol, the member's index, name and bodies, the columns,
#     the `# active` line both ways, and every other header line this tree
#     writes, P1's `# zonal` line among them;
#   - with --producer-files, a directory that is not one, or is --keep's;
#     a member's file that is not of its form or has no run key; and a file
#     whose run key is not its member's run (THE PRODUCER'S OWN FILES).
#
# FAIL (exit 1), each by member, with its name:
#   - a member missing from the record or incomplete in it: its system
#     line, its bodies, its sample lines (as many as steps and sample make,
#     at the steps they make) and its trailer line. EVERY member is held to
#     this, sampled or not;
#   - a sampled member whose lines are not this tree's: the first differing
#     line by the record's line number, its sample and step, and the column
#     named from the record's own `# columns:` line;
#   - a sampled member whose trailer is not this tree's. Nothing compares
#     the trailer for you: steps_done, steps_rejected,
#     iterations_max_exceeded, mean_pc_iterations and max_pc_iterations are
#     compared here, and the members' rejections are summed against the
#     record's total.
#
# THE SAMPLE. Seeded and stratified by member index: the members are cut
# into --strata contiguous ranges, and each range gives the member with the
# least SHA-256 of (seed, index), so a seed picks the same members on every
# host and every Python. P2's generator sorts waves by radius, so the index
# is radius. Every member whose trailer shows a rejected step is added.
# The sample trusts the trailers it is testing (V4 measured it): a member
# whose trailer under-reports its rejections leaves the sample, and an
# unsampled member's lines are never rerun. Its trailer is still held to
# the record's total and its lines to their places, but its bits are not
# compared. The driver says so beside the sample.
#
# THE COUNTERS, beside each member's verdict:
#   - steps, rejections, mean and largest passes and cap hits, from the
#     member's trailer (and equal to the rerun's, or the member FAILs);
#   - the closest approach of each test particle to each active body AT
#     THE SAMPLES, not between them, from the record's own positions (its
#     x y z columns, found by name), in binary64;
#   - with --format-diff F2, the format difference: the member rerun here
#     with --dt-out, then at F2 with --dt-file and --epsilon 0, so on the
#     same accepted steps, measured by tools/compare_formats.py. Its scope
#     is printed beside it, by the lead's ruling on P1's measurement: it
#     shows the narrower format's round-off on these steps, and a narrowing
#     both runs share passes it. It is a counter, never evidence that the
#     arithmetic is the format's everywhere.
#     The replay takes no rejected attempt. A rejection before the first
#     accepted step leaves state (that attempt's corrected b is the next
#     one's first guess), so that member's difference includes it. A
#     rejection after an accepted step leaves no state bit, and the replay
#     follows the member's path; only the dt_next of a sample line printed
#     just before a rejected attempt differs. P3 measured both, and
#     tools/check_compare_run.py holds both. The dt file's first step is
#     the record's dt0 exactly when the first attempt was accepted, which
#     is how the note beside the number tells them apart.
#
# PARTICLES, through P3's replay (the lead's P3 x P4 seam, ROUND4.md "What
# the integrator keeps"; sample_particles()). Each rerun member's rerun
# also writes its attempts record and corrector schedule (--attempts-out,
# --record-iters; neither moves a record line). From each such member two
# test particles are sampled: the one whose closest approach to an active
# body, at the samples, is the least (a rejection's path can show in that
# particle alone: P3, measured), and one more by the seed, as the members
# are. Each is written alone with the member's massive bodies, and its
# zonal line if any, by P3's tools/extract_particle.py, replayed here from
# those records, and held to the RECORD's own lines, not to the rerun's:
# every sample's fields before the bodies, the particle and every massive
# body, value for value, the exact-time pair, and the trailer's counters.
# A replay that is not the record's FAILs, naming the member, the particle,
# the line and the column. Beside each particle, its replay's counters and
# its closest approach to each active body at the samples.
# Without --producer-files the attempts records are this tree's reruns', so
# a particle is sampled only in a member that was rerun.
#
# THE PRODUCER'S OWN FILES (ROUND5.md, P3: a GPU run of 10^6 particles
# cannot afford a rerun of every member whose particles are sampled).
# --producer-files DIR names a directory holding, for member k (its index in
# the problem file: system k of an ensemble file, the k of a --member k
# record, 0 for a single problem), the producer's own
#     DIR/member<k>.att    the attempts record, in --attempts-out's form
#     DIR/member<k>.sched  the corrector schedule, in --record-iters's form
# each written by the member's run alone and keyed, on its `# run` line,
# with exactly what `ias15_cft --run-key` prints for that run: the arguments
# the member is rerun with here, its engine the record's own if this tree
# has it and --rerun-engine's if not (PRODUCER_FILES, below; P2's kernel
# writes them, and tools/check_compare_run.py writes its producer's by
# these names, imported from here). Then:
#   - the members whose particles are sampled are the member sample's and,
#     besides them, one from each of --particle-strata ranges of member
#     index (default: as many as --strata), by the least SHA-256 of (seed,
#     index) under a name of their own, so the two samples are drawn apart;
#   - before anything runs, each pair found is held to its member's run by
#     its key, the program's own (--run-key): a file keyed for another wave,
#     format, active line or zonal line is REFUSED, naming the member, the
#     file and every key field that differs, as --replay-iters refuses a
#     schedule recorded from another run;
#   - each sampled particle is replayed against the producer's own
#     attempts record and schedule, with no rerun of its member;
#   - a member's own lines are still held as without the files: rerun here
#     where the member sample holds it, and by the record alone where it
#     does not (and by its particles' replays, which hold every massive
#     body and the exact-time pair). The driver says which, and why;
#   - a member whose pair is not both in DIR falls back to today's rerun,
#     which then holds its lines too and writes the files its particles
#     replay from; the driver names the file that is not there.
# What the replay cannot hold, stated: it takes the producer's decisions
# as given - each attempt's dt and outcome, and its passes - and holds the
# arithmetic of each attempt for the particle, every massive body and the
# exact-time pair. The decisions are held only where the member is rerun,
# by this tree's own step control and exit test: a producer whose exit test
# took another pass, its record that schedule's replay, passes in a member
# that is not rerun and FAILs in one that is (tools/check_compare_run.py
# measures both). And a rejected attempt after an accepted step moves no
# state bit (docs/VALIDATION.md entry 49): its dt shows only as the dt_next
# of a sample line printed just before it.
# Without P3's tools in the tree the driver prints
#     particles: not sampled (no particle replay in this tree)
# and with them present but not callable (an extract() that will not
# import) it FAILs and exits 1, whatever the members did.
#
# LIMITS, stated:
#   - the record does not carry --adaptive-mode, --softening, --min-dt or
#     --no-flag-abort (src/ias15_cft.c main(), read at 9a42b1d), so the
#     rerun takes their defaults. A run made otherwise passes here only if
#     the setting changed no bit of it (a floor that never bound), and
#     otherwise fails as not identical.
#   - a record run from a dt file (dt_file= in its header) is refused: the
#     file is not in the record.
#   - a problem's zonal line is held through the record's `# zonal` line,
#     compared verbatim with the rerun's like every other header line, and
#     through the bits the term moves; the term's own correctness is
#     tools/check_zonal.py's.
#   - the record is read twice, as a stream, so its size is not held in
#     memory; the problem file is.
#   - the sample trusts the trailers it is testing (THE SAMPLE, above).
#
#   exit 0: every sampled member and particle identical; 1: a FAIL (a
#   member, a particle, or P3's replay present and not callable); 2: refused

import argparse
import hashlib
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import time

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import hex_to_fraction               # noqa: E402
from check_ensemble import parse as parse_samples           # noqa: E402
from check_records import header_args                       # noqa: E402

root = os.path.dirname(here)

KINDS = {"# cft-rebound ias15 record v1": "single", "# cft-rebound ias15 ensemble record v1": "ensemble"}
# Everything a rerun is built from. The record must name each: a default
# would be a guess about somebody else's run.
NEED_PROGRAM = ("format", "arith", "cs", "engine", "pc_tol_shift", "max_iter", "steps", "sample", "N", "problem")
NEED_DT0 = ("dt0", "epsilon", "G", "pc_tol")
# The program fields a rerun must report as the record does: the run, not
# who ran it or where. backend, impl and the identity line are the
# producer's own, and are printed, not compared.
SAME_PROGRAM = ("format", "arith", "cs", "pc_tol_shift", "max_iter", "steps", "sample", "N", "problem")
# A member's counters: what an ensemble's trailer carries for each system.
COUNTERS = ("steps_done", "steps_rejected", "iterations_max_exceeded", "mean_pc_iterations", "max_pc_iterations")
FORMATS = ("fp64", "fp128", "fp256")          # narrow to wide
# The problem format's body fields after the name: src/ias15_cft.c,
# read_problem(), "body NAME m x y z vx vy vz". The record's per-body
# columns are matched to these by name.
PROBLEM_FIELDS = ("m", "x", "y", "z", "vx", "vy", "vz")
# The record cuts names to its fields' widths (read_problem: the problem's
# name to 63 characters, a system's and a body's to 31).
NAME_WIDTH, SYS_WIDTH = 63, 31
# P3's tools (ROUND4.md, P3). Either one in the tree means the particle
# place below is due.
P3_TOOLS = ("extract_particle.py", "check_replay.py")
# A producer's own files, one pair a member k (ROUND5.md, P3), in the forms
# ias15_cft's --attempts-out and --record-iters write. The one statement of
# the names: tools/check_compare_run.py imports it to write its producer's.
PRODUCER_FILES = (("attempts record", "member%d.att"), ("corrector schedule", "member%d.sched"))
# The schedule's first line, as src/ias15_cft.c's sched_open() writes it and
# --replay-iters and tools/extract_particle.py hold it; the attempts
# record's is extract_particle's ATT_MAGIC. A copy that drifted would refuse
# every producer's file, loudly, never pass a wrong one.
SCHED_MAGIC = "# cft-rebound corrector schedule v1"
ENGINES = ("loop", "program")                 # this tree's, as --engine takes them
FLOOR = re.compile(r"^FLOOR (\S+) vs (\S+): (\d+) common samples, last step (\d+), max \|dE/E\| (\S+) \(last \S+\), "
                   r"max \|dx\|/L (\S+) \(last \S+\)$", re.M)
# What a format difference can and cannot show, printed beside its number
# (the lead's ruling, .ledger/urgent/lead-gate4-ruling.md, 2026-11-12, on
# P1's measurement that rounding a term to binary64 in BOTH runs leaves
# binary128 against binary256 where the correct build has it).
FORMAT_DIFF_SCOPE = "the difference shows the narrower format's round-off on these steps; a narrowing both runs share passes it"
# What the sample cannot see, printed beside it (V4, 2026-11-12): it is
# chosen from the trailers under test.
SAMPLE_LIMIT = ("sample: it trusts the trailers it is testing: a member whose trailer under-reports its rejections "
                "leaves the sample, and an unsampled member's lines are never rerun")


class Refused(Exception):
    """Not comparable: refused by name, exit 2."""


def same_name(rec, prob, width):
    return rec == prob or (len(prob) > width and rec == prob[:width])


def kv(tokens):
    d = {}
    for t in tokens:
        k, eq, v = t.partition("=")
        if eq:
            d[k] = v
    return d


def exact(tok):
    """A hex float's exact value and sign (check_equivalence.py's), or None
    for a token that is not one."""
    try:
        return hex_to_fraction(tok)
    except (ValueError, IndexError):
        return None


class Columns:
    """A sample line's layout, from the record's own `# columns:` line:

        # columns: [system s] sample step t dt_next dt_last E then per body
                   x y z vx vy vz ; then exact-time pair t_hi t_lo

    "system s" (an ensemble's) is the word system and the system's index;
    "sample" is the word sample, and names the field after it, the sample's
    index; each name up to "then per body" is one field; then each body's
    fields; then "|" and the pair. Nothing here is transcribed but the two
    leading fields check_ensemble.parse() reads by position, and those are
    asserted."""

    def __init__(self, line, N):
        text = line.split(":", 1)[1].strip() if ":" in line else ""
        head, s1, rest = text.partition(" then per body ")
        per, s2, pair = rest.partition(" ; then exact-time pair ")
        if not (s1 and s2):
            raise Refused("the `# columns:` line does not say \"then per body\" and \"; then exact-time pair\": %r"
                          % line)
        head, per, pair = head.split(), per.split(), pair.split()
        self.ensemble = head[:1] == ["system"]
        if self.ensemble:
            head = head[2:]
        if head[:2] != ["sample", "step"]:
            raise Refused("the columns do not begin \"sample step\", the two fields check_ensemble.parse() reads "
                          "first: %r" % line)
        lost = [c for c in ("x", "y", "z") if c not in per]
        if lost:
            raise Refused("the columns name no %s for a body: %r" % (" or ".join(lost), line))
        self.N = N
        self.scalars, self.per, self.pair = head, per, pair
        self.base = 2 if self.ensemble else 0             # "system s"
        self.first_body = self.base + 1 + len(head)      # the word sample, then the scalars
        self.bar = self.first_body + N * len(per)
        self.ntok = self.bar + 1 + len(pair)

    def same_layout(self, other):
        return (self.scalars, self.per, self.pair) == (other.scalars, other.per, other.pair)

    def body_tok(self, b, name):
        return self.first_body + b * len(self.per) + self.per.index(name)

    def value_name(self, j, body_names):
        """Value j of check_ensemble.parse()'s list - every field after the
        sample's index and step, the bar dropped - by its column's name."""
        ns = len(self.scalars) - 2
        if j < ns:
            return self.scalars[2 + j]
        j -= ns
        if j < self.N * len(self.per):
            b, f = divmod(j, len(self.per))
            return "%s of body %d (%s)" % (self.per[f], b, body_names[b] if b < len(body_names) else "?")
        j -= self.N * len(self.per)
        return "the exact-time pair's %s" % self.pair[j] if j < len(self.pair) else "value %d past the columns" % j

    def value_tok(self, j):
        """Where value j is in a split sample line."""
        t = self.base + 3 + j
        return t if t < self.bar else t + 1

    def split(self, toks):
        """(system, sample index, step) of a split sample line, or why it is
        not one of this layout's."""
        if len(toks) != self.ntok or toks[self.bar] != "|":
            return "%d fields, where the columns and N=%d make %d with the bar at field %d" % (
                len(toks), self.N, self.ntok, self.bar + 1)
        if self.ensemble and toks[0] != "system":
            return "it does not begin with \"system\""
        if toks[self.base] != "sample":
            return "\"sample\" is not where the columns put it"
        try:
            return (int(toks[1]) if self.ensemble else 0), int(toks[self.base + 1]), int(toks[self.base + 2])
        except ValueError:
            return "its system, sample or step is not a whole number"


class Record:
    """A record's header, read up to its first sample line, and - after
    scan() - where every later line belongs."""

    def __init__(self, label):
        self.label = label
        self.kind = None
        self.program = {}
        self.dt0 = {}
        self.cols_line = None
        self.cols = None
        self.active = None              # (n, N, pairs), from `# active`
        self.sysname = {}               # s -> name, `# system s NAME`
        self.bodies = {}                # s -> [(name, m)], `# body`
        self.other = []                 # every other header line, verbatim
        self.first = None               # (lineno, line) after the header
        self.nsamples = {}
        self.misplaced = {}             # s -> (lineno, k, step, n, step it should be)
        self.trailer = {}               # s -> (lineno, line, counters)
        self.totals = None              # (lineno, counters)
        self.strays = []                # (lineno, why)

    @property
    def E(self):
        return int(self.program["E"]) if self.kind == "ensemble" else 1

    @property
    def N(self):
        return int(self.program["N"])

    @property
    def member(self):
        """The problem's system a one-system record ran (--member k), or
        None for a record of a whole problem."""
        return int(self.program["member"]) if self.kind == "single" and "member" in self.program else None

    def label_of(self, s):
        """A member by its index in the problem: an ensemble's own, or the
        --member a one-system record ran (0 for a whole single problem)."""
        return s if self.kind == "ensemble" or self.member is None else self.member

    def name_of(self, s, prob=None):
        if self.sysname.get(s):
            return self.sysname[s]
        if self.kind == "single" and self.program.get("system"):
            return self.program["system"]
        ps = prob.system_for(self, s) if prob is not None else None
        return ps[0] if ps is not None and ps[0] else "?"


def numbered(path):
    with open(path, encoding="ascii", errors="replace") as f:
        for n, line in enumerate(f, 1):
            yield n, line.rstrip("\r\n")


def read_header(lines, label):
    """lines: an iterator of (lineno, text), which is left at the line after
    the header. The header is held to what a rerun needs, or refused."""
    rec = Record(label)
    cur = None
    for lineno, line in lines:
        if lineno == 1:
            rec.kind = KINDS.get(line)
            if rec.kind is None:
                raise Refused("%s is not a cft-rebound ias15 record: its first line is %r" % (label, line[:100]))
            continue
        if not line.startswith("#") or "steps_done=" in line:
            rec.first = (lineno, line)
            break
        if line.startswith("# program="):
            rec.program = header_args(line)
        elif line.startswith("# dt0="):
            rec.dt0 = header_args(line) or kv(line[2:].split())    # ias15_ref writes it unspaced
        elif line.startswith("# identity"):
            pass                                                    # the producer's build: its own
        elif line.startswith("# columns:"):
            rec.cols_line = line
        elif line.startswith("# system "):
            t = line.split()
            if len(t) < 3 or not t[2].isdigit():
                raise Refused("%s line %d: a system line with no system index: %r" % (label, lineno, line))
            cur = int(t[2])
            rec.sysname[cur] = t[3] if len(t) > 3 else ""
            rec.bodies.setdefault(cur, [])
        elif line.startswith("# body "):
            t = line.split()
            m = t[5] if len(t) > 5 and t[4] == "m=" else t[4][2:] if len(t) > 4 and t[4].startswith("m=") else None
            s = cur if rec.kind == "ensemble" else 0
            if s is None or m is None:
                raise Refused("%s line %d: a body with no system line before it, or no mass: %r" % (label, lineno, line))
            rec.bodies.setdefault(s, []).append((t[3], m))
        else:
            if line.startswith("# active "):
                m = re.match(r"^# active (\d+) of (\d+) bodies a system: (\d+) pairs", line)
                if not m:
                    raise Refused("%s line %d: an active line this driver cannot read: %r" % (label, lineno, line))
                rec.active = tuple(int(x) for x in m.groups())
            rec.other.append(line)
    if rec.kind is None:
        raise Refused("%s is empty" % label)
    missing = [k for k in NEED_PROGRAM if k not in rec.program] + \
              (["E"] if rec.kind == "ensemble" and "E" not in rec.program else [])
    if missing:
        raise Refused("%s names no %s in its `# program=` line (it is program=%s impl=%s): the rerun is built from "
                      "the record's own arguments, and a default would be a guess"
                      % (label, ", ".join(missing), rec.program.get("program", "?"), rec.program.get("impl", "?")))
    missing = [k for k in NEED_DT0 if k not in rec.dt0]
    if missing:
        raise Refused("%s names no %s in its `# dt0=` line" % (label, ", ".join(missing)))
    if "dt_file" in rec.program:
        raise Refused("%s was run from a dt file (dt_file=%s), and the file is not in the record"
                      % (label, rec.program["dt_file"]))
    if rec.cols_line is None:
        raise Refused("%s has no `# columns:` line" % label)
    try:
        if int(rec.program["sample"]) < 1 or int(rec.program["steps"]) < 0 or rec.N < 1 or rec.E < 1:
            raise ValueError
    except ValueError:
        raise Refused("%s has N=%s E=%s steps=%s sample=%s" % (label, rec.program.get("N"), rec.program.get("E", 1),
                                                               rec.program["steps"], rec.program["sample"]))
    rec.cols = Columns(rec.cols_line, rec.N)
    if rec.cols.ensemble != (rec.kind == "ensemble"):
        raise Refused("%s is a%s record, and its columns line %s a system field" % (
            label, "n ensemble" if rec.kind == "ensemble" else " one-system", "has" if rec.cols.ensemble else "has no"))
    return rec


class Problem:
    """A problem file, as src/ias15_cft.c's read_problem() reads it: `#` and
    blank lines skipped; name, G, N, E, system, body, active and zonal; any
    other line refused by name. Body lines are kept as text."""

    def __init__(self, path):
        self.path = path
        self.name = self.G = self.N = self.active = self.zonal = None
        self.systems = []                # [name, [body lines]]
        loose = []
        if not os.path.isfile(path):
            raise Refused("there is no problem file %s" % path)
        with open(path, encoding="ascii", errors="replace") as f:
            for n, line in enumerate(f, 1):
                if line.startswith("#") or not line.strip():
                    continue
                t = line.split()
                if t[0] == "name" and len(t) > 1:
                    self.name = t[1]
                elif t[0] == "G" and len(t) > 1:
                    self.G = t[1]
                elif t[0] == "N" and len(t) > 1:
                    if not t[1].isdigit():
                        raise Refused("the problem file %s line %d: N %s is not a count" % (path, n, t[1]))
                    self.N = int(t[1])
                elif t[0] == "E":
                    pass                                   # the system lines are counted
                elif t[0] == "system":
                    if loose:
                        raise Refused("the problem file %s has body lines before its first system line" % path)
                    self.systems.append([t[1] if len(t) > 1 else "", []])
                elif t[0] == "body":
                    if len(t) < 9:
                        raise Refused("the problem file %s line %d is a short body line" % (path, n))
                    (self.systems[-1][1] if self.systems else loose).append(line.strip())
                elif t[0] == "active":
                    if self.active is not None:
                        raise Refused("the problem file %s has a second active line" % path)
                    self.active = t[1] if len(t) > 1 else ""
                elif t[0] == "zonal":
                    self.zonal = " ".join(t[1:])
                else:
                    raise Refused("the problem file %s line %d is %r, a line ias15_cft does not read" % (path, n, t[0]))
        self.single = not self.systems
        if self.single:
            self.systems = [[self.name or "", loose]]
        if self.N is None:
            self.N = len(self.systems[0][1])

    def system_for(self, rec, s):
        """The problem's [name, body lines] that the record's system s ran."""
        i = rec.member if rec.member is not None else s
        return self.systems[i] if 0 <= i < len(self.systems) else None


def hold_problem(rec, prob):
    """The problem file against the record's header, before anything runs."""
    where = "the problem file %s" % prob.path
    if prob.name is None or not same_name(rec.program["problem"], prob.name, NAME_WIDTH):
        raise Refused("%s is named %s, and the record ran problem=%s" % (where, prob.name, rec.program["problem"]))
    if prob.N != rec.N:
        raise Refused("%s has N=%d bodies a system, and the record N=%d" % (where, prob.N, rec.N))
    if rec.kind == "ensemble" and (prob.single or len(prob.systems) != rec.E):
        raise Refused("%s has %d system(s), and the record E=%d" % (where, len(prob.systems), rec.E))
    if rec.kind == "single" and rec.member is None and len(prob.systems) != 1:
        raise Refused("%s has %d systems, and the record is one system run without --member" % (where, len(prob.systems)))
    if rec.member is not None and prob.system_for(rec, 0) is None:
        raise Refused("%s has %d system(s), and the record ran --member %d" % (where, len(prob.systems), rec.member))
    if prob.G is None or exact(prob.G) is None or exact(prob.G) != exact(rec.dt0["G"]):
        raise Refused("%s has G %s, and the record G=%s" % (where, prob.G, rec.dt0["G"]))
    ra = rec.active[0] if rec.active else None
    if (str(ra) if ra is not None else None) != prob.active:
        raise Refused("%s has %s, and the record %s" % (
            where, "active %s" % prob.active if prob.active is not None else "no active line",
            "says active %d" % ra if ra is not None else "has no `# active` line"))
    if rec.active and rec.active[1] != rec.N:
        raise Refused("the record's `# active` line counts %d bodies a system, and its N=%d" % (rec.active[1], rec.N))
    if rec.member is not None and "system" in rec.program and \
            not same_name(rec.program["system"], prob.system_for(rec, 0)[0], SYS_WIDTH):
        raise Refused("%s's system %d is named %s, and the record ran system=%s"
                      % (where, rec.member, prob.system_for(rec, 0)[0], rec.program["system"]))
    for s in sorted(rec.bodies):
        ps = prob.system_for(rec, s)
        if ps is None or len(rec.bodies[s]) != len(ps[1]):
            continue                     # the record's own defect: that member is reported FAIL
        if rec.kind == "ensemble" and s in rec.sysname and not same_name(rec.sysname[s], ps[0], SYS_WIDTH):
            raise Refused("%s's system %d is named %s, and the record's %s" % (where, s, ps[0], rec.sysname[s]))
        for b, ((rname, rm), pl) in enumerate(zip(rec.bodies[s], ps[1])):
            pt = pl.split()
            if not same_name(rname, pt[1], SYS_WIDTH):
                raise Refused("%s's system %d body %d is %s, and the record's %s" % (where, s, b, pt[1], rname))
            if exact(pt[2]) is None or exact(pt[2]) != exact(rm):
                raise Refused("%s's system %d (%s) body %d (%s) has m %s, and the record m=%s"
                              % (where, s, ps[0], b, pt[1], pt[2], rm))


def hold_initial(rec, prob, s, toks, lineno):
    """Sample 0 of system s is its initial state: the problem's, exactly."""
    ps = prob.system_for(rec, s)
    if ps is None or len(ps[1]) != rec.N:
        return
    for b, pl in enumerate(ps[1]):
        pt = pl.split()
        for c in rec.cols.per:
            if c not in PROBLEM_FIELDS:
                continue
            want, got = pt[2 + PROBLEM_FIELDS.index(c)], toks[rec.cols.body_tok(b, c)]
            if exact(want) is None or exact(want) != exact(got):
                raise Refused("the problem file %s is not the problem the record ran, though both are named %s: "
                              "system %d (%s) body %d (%s) starts at %s %s in the file, and at %s in the record's "
                              "sample 0 (line %d)" % (prob.path, prob.name, s, ps[0], b, pt[1], c, want, got, lineno))


def scan(rec, lines, prob):
    """Pass 1, after the header: every line's place, without keeping the
    lines, and every system's sample 0 against the problem."""
    steps, sample = int(rec.program["steps"]), int(rec.program["sample"])
    in_trailer = False

    def one(lineno, line):
        nonlocal in_trailer
        if line.startswith("#"):
            if rec.kind == "ensemble" and line.startswith("# system ") and "steps_done=" in line:
                t = line.split()
                try:
                    s = int(t[2])
                except ValueError:
                    rec.strays.append((lineno, "a trailer line with no system index"))
                    return
                if s in rec.trailer:
                    rec.strays.append((lineno, "a second trailer line for system %d" % s))
                rec.trailer[s] = (lineno, line, kv(t[3:]))
                in_trailer = True
            elif line.startswith("# steps_done="):
                rec.totals = (lineno, kv(line[2:].split()))
                if rec.kind == "single":
                    rec.trailer[0] = (lineno, line, rec.totals[1])
                in_trailer = True
            return                          # `# route` and the rest: the producer's own
        if not line.strip():
            return
        if in_trailer:
            rec.strays.append((lineno, "a sample line inside the trailer"))
        toks = line.split()
        r = rec.cols.split(toks)
        if isinstance(r, str):
            rec.strays.append((lineno, "not a sample line of this record: " + r))
            return
        s, k, step = r
        if not 0 <= s < rec.E:
            rec.strays.append((lineno, "a sample line for system %d, past the record's E=%d" % (s, rec.E)))
            return
        n = rec.nsamples.get(s, 0)
        rec.nsamples[s] = n + 1
        if (k != n or step != min(n * sample, steps)) and s not in rec.misplaced:
            rec.misplaced[s] = (lineno, k, step, n, min(n * sample, steps))
        if n == 0 and k == 0:
            hold_initial(rec, prob, s, toks, lineno)

    if rec.first:
        one(*rec.first)
    for lineno, line in lines:
        one(lineno, line)


def hold_structure(rec):
    """Every member present and complete, sampled or not: ({system: [what is
    wrong]}, {system: True when nothing of it is in the record}, [what is
    wrong with the record as a whole])."""
    steps, sample = int(rec.program["steps"]), int(rec.program["sample"])
    want = 1 + (steps + sample - 1) // sample
    bad, gone = {}, {}
    for s in range(rec.E):
        why = []
        if rec.kind == "ensemble" and s not in rec.sysname:
            why.append("no `# system %d` line" % s)
        nb = len(rec.bodies.get(s, []))
        if nb != rec.N:
            why.append("%d of %d `# body` lines" % (nb, rec.N))
        ns = rec.nsamples.get(s, 0)
        if ns != want:
            why.append("%d of %d sample lines" % (ns, want))
        if s in rec.misplaced:
            why.append("line %d is sample %d at step %d, where its place makes it sample %d at step %d" % rec.misplaced[s])
        if s not in rec.trailer:
            why.append("no trailer line")
        else:
            lineno, _, t = rec.trailer[s]
            lost = [c for c in COUNTERS if c not in t]
            if lost:
                why.append("its trailer (line %d) has no %s" % (lineno, ", ".join(lost)))
            elif t["steps_done"] != str(steps):
                why.append("its trailer (line %d) says steps_done=%s, and the run is steps=%d" % (lineno, t["steps_done"], steps))
        if why:
            bad[s] = why
            gone[s] = (nb == 0 and ns == 0 and s not in rec.trailer and
                       (rec.kind != "ensemble" or s not in rec.sysname))
    whole = ["line %d: %s" % x for x in rec.strays]
    if rec.totals is None:
        whole.append("no `# steps_done=` trailer: the run did not finish")
    elif rec.kind == "ensemble" and "steps_rejected" in rec.totals[1]:
        try:
            total = sum(int(rec.trailer[s][2]["steps_rejected"]) for s in rec.trailer)
            if total != int(rec.totals[1]["steps_rejected"]):
                whole.append("the members' trailers reject %d steps, and the record's total (line %d) says "
                             "steps_rejected=%s" % (total, rec.totals[0], rec.totals[1]["steps_rejected"]))
        except (KeyError, ValueError):
            pass                            # a member's trailer lacks it: that member is reported
    return bad, gone, whole


def stratum_pick(seed, lo, hi):
    return min(range(lo, hi), key=lambda m: hashlib.sha256(b"cft-rebound compare_run %d %d" % (seed, m)).digest())


def plan(E, strata, seed, rejected):
    """The strata, each one's pick, and the members to rerun."""
    K = max(1, min(strata, E))
    ranges = [(i * E // K, (i + 1) * E // K) for i in range(K)]
    picks = [stratum_pick(seed, lo, hi) for lo, hi in ranges]
    return ranges, picks, sorted(set(picks) | set(rejected))


def particle_plan(E, strata, seed):
    """With the producer's files, the members whose particles are sampled
    besides the member sample's: one from each of `strata` ranges of member
    index, by the least SHA-256 of (seed, index) under a name of its own, so
    that it is drawn apart from the member sample. None for strata 0."""
    if strata < 1:
        return [], []
    K = min(strata, E)
    ranges = [(i * E // K, (i + 1) * E // K) for i in range(K)]
    picks = [min(range(lo, hi), key=lambda m: hashlib.sha256(b"cft-rebound compare_run particle member %d %d"
                                                           % (seed, m)).digest()) for lo, hi in ranges]
    return ranges, picks


def producer_paths(directory, k):
    """Member k's attempts record and corrector schedule in a producer's
    directory (--producer-files)."""
    return tuple(os.path.join(directory, name % k) for _, name in PRODUCER_FILES)


def run_key(rec, args, s, engine):
    """Member s's run key, as its schedule's and its attempts record's `# run`
    lines hold it: the program's own (ias15_cft --run-key, on the arguments
    the member is rerun with), so that nothing here hashes anything and no
    copy of the key can drift from sched_key()."""
    rc, out, err = run(rerun_cmd(rec, args, s, engine=engine) + ["--run-key"])
    if rc != 0 or len(out.splitlines()) != 1:
        raise Refused("this tree's ias15_cft gave no run key for member %d's run (rc %d): %s"
                      % (rec.label_of(s), rc, (err or out).strip()[-300:]))
    return out.strip()


def file_key(k, what, path, magic):
    """A producer's file's run key, its first two lines held to its form."""
    try:
        with open(path, encoding="ascii", errors="replace") as f:
            head = [f.readline().rstrip("\r\n") for _ in range(2)]
    except OSError as e:
        raise Refused("member %d's %s %s cannot be read: %s" % (k, what, path, e))
    if head[0] != magic:
        raise Refused("member %d's %s %s is not one: its first line is %r, where the form's is %r"
                      % (k, what, path, head[0][:100], magic))
    if not head[1].startswith("# run "):
        raise Refused("member %d's %s %s has no run key: its second line is %r, not `# run <key>`"
                      % (k, what, path, head[1][:100]))
    return head[1][6:]


def hold_key(k, what, path, have, want):
    """A producer's file against its member's run, by the run key: refused by
    name, every field that differs named with both values."""
    if have == want:
        return
    a, b = kv(have.split()), kv(want.split())
    fields = [f for f in b if a.get(f) != b[f]] + [f for f in a if f not in b]
    raise Refused("member %d's %s %s was recorded for another run, and a schedule replays only its own run: its run "
                  "key differs from member %d's in %s (%s)" % (
                      k, what, path, k, ", ".join(fields) or "its text, field for field alike",
                      "; ".join("%s: %s in the file, %s for member %d" % (f, a.get(f, "absent"), b.get(f, "absent"), k)
                                for f in fields)))


def extract(path, rec, members):
    """Pass 2: the chosen members' sample lines, with their line numbers."""
    want = set(members)
    got = {s: [] for s in members}
    for lineno, line in numbered(path):
        if line.startswith("#"):
            continue
        if rec.kind == "ensemble":
            t = line.split(None, 2)
            if len(t) < 2 or t[0] != "system" or not t[1].isdigit():
                continue
            s = int(t[1])
        elif line.startswith("sample "):
            s = 0
        else:
            continue
        if s in want:
            got[s].append((lineno, line))
    return got


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stdout, p.stderr


def rerun_cmd(rec, args, s, fmt=None, dt_file=None, dt_out=None, att_out=None, sched_out=None, engine=None):
    p = rec.program
    engine = engine or args.rerun_engine or p["engine"]
    cmd = [args.exe, "--format", fmt or p["format"], "--problem", args.problem, "--steps", p["steps"],
           "--sample", p["sample"], "--cs", p["cs"], "--arith", p["arith"], "--engine", engine,
           "--max-iter", p["max_iter"], "--quiet"]
    if dt_file:
        cmd += ["--dt-file", dt_file, "--epsilon", "0"]            # at that format's own pc_tol_shift
    else:
        cmd += ["--dt", rec.dt0["dt0"], "--epsilon", rec.dt0["epsilon"], "--pc-tol-shift", p["pc_tol_shift"]]
    if dt_out:
        cmd += ["--dt-out", dt_out]
    if att_out:
        # the attempts record and the corrector schedule P3's replay reads;
        # writing them moves no record line
        cmd += ["--attempts-out", att_out, "--record-iters", sched_out]
    member = s if rec.kind == "ensemble" else rec.member
    if member is not None:
        cmd += ["--member", str(member)]
    if engine == "program":
        cmd += ["--programs", args.programs]
    return cmd


def hold_rerun(rec, s, mine, engine):
    """The rerun's own header against the record's: the same run, or refused."""
    what = "the rerun of member %d is not the run the record describes" % rec.label_of(s)
    if mine.kind != "single":
        raise Refused("%s: it wrote an ensemble record" % what)
    want = s if rec.kind == "ensemble" else rec.member
    if want is not None and mine.program.get("member") != str(want):
        raise Refused("%s: it ran member=%s, and was asked for %d" % (what, mine.program.get("member"), want))
    for k in SAME_PROGRAM:
        if mine.program.get(k) != rec.program.get(k):
            raise Refused("%s: %s=%s here, %s=%s in the record" % (what, k, mine.program.get(k), k, rec.program.get(k)))
    if mine.program.get("engine") != engine:
        raise Refused("%s: engine=%s here, and it was asked for %s" % (what, mine.program.get("engine"), engine))
    for k in NEED_DT0:
        if exact(mine.dt0[k]) != exact(rec.dt0[k]):
            raise Refused("%s: %s=%s here, %s=%s in the record" % (what, k, mine.dt0[k], k, rec.dt0[k]))
    want = rec.sysname.get(s) if rec.kind == "ensemble" else rec.program.get("system")
    if want is not None and mine.program.get("system") != want:
        raise Refused("%s: it ran system=%s, and the record's system %d is %s" % (what, mine.program.get("system"), s, want))
    rb, mb = rec.bodies.get(s, []), mine.bodies.get(0, [])
    if len(rb) != len(mb):
        raise Refused("%s: %d bodies here, %d in the record" % (what, len(mb), len(rb)))
    for b, ((rn, rm), (mn, mm)) in enumerate(zip(rb, mb)):
        if rn != mn or exact(rm) != exact(mm):
            raise Refused("%s: body %d is %s m=%s here, and %s m=%s in the record" % (what, b, mn, mm, rn, rm))
    if not mine.cols.same_layout(rec.cols):
        raise Refused("%s: its columns are %r, the record's %r" % (what, mine.cols_line, rec.cols_line))
    # both ways for the active line: massless test particles run all active
    # give the same bits, so only this line says which run it was
    if mine.active != rec.active:
        if not mine.active or not rec.active:
            raise Refused("%s: %s here, and %s in the record" % (
                what, "active %d" % mine.active[0] if mine.active else "no `# active` line",
                "active %d" % rec.active[0] if rec.active else "no `# active` line"))
        for field, a, b in zip(("active count", "bodies a system", "pairs a force evaluation"), mine.active, rec.active):
            if a != b:
                raise Refused("%s: the `# active` line's %s is %d here, and %d in the record" % (what, field, a, b))
    # the rest one way: every line this tree writes must be the record's
    # too; a line only the record has is its producer's own
    absent = [l for l in mine.other if l not in rec.other]
    if absent:
        raise Refused("%s: this tree writes %r, and the record has no such line" % (what, absent[0]))


def compare(rec, s, lines, mine, mine_text, names):
    """The member's lines and trailer against the rerun's, through
    check_ensemble.parse(). Returns (identical, what, record's counters)."""
    tl, tline, _ = rec.trailer[s]
    text = "\n".join(l for _, l in lines) + "\n" + tline + "\n"
    try:
        rs, rt = parse_samples(text, s if rec.kind == "ensemble" else None)
        ms, mt = parse_samples(mine_text)
    except (ValueError, IndexError) as e:
        return False, "a sample line does not parse as hex floats: %s" % e, {}
    if len(rs) != len(lines):
        return False, "%d of its %d lines parse as its samples" % (len(rs), len(lines)), rt
    if len(rs) != len(ms):
        return False, "%d samples in the record, %d here" % (len(rs), len(ms)), rt
    mine_lines = [l for l in mine_text.splitlines() if l.startswith("sample ")]
    nvals = 0
    for i, ((rk, rstep, rv), (mk, mstep, mv)) in enumerate(zip(rs, ms)):
        lineno, line = lines[i]
        where = "first differing line: record line %d (%ssample %d, step %d)" % (
            lineno, "system %d " % s if rec.kind == "ensemble" else "", rk, rstep)
        if (rk, rstep) != (mk, mstep):
            return False, "%s: it is sample %d at step %d here" % (where, mk, mstep), rt
        if len(rv) != len(mv):
            return False, "%s: %d values in the record, %d here" % (where, len(rv), len(mv)), rt
        for j, (a, b) in enumerate(zip(rv, mv)):
            if a != b:
                return False, "%s: %s: record %s, here %s" % (
                    where, rec.cols.value_name(j, names), line.split()[rec.cols.value_tok(j)],
                    mine_lines[i].split()[mine.cols.value_tok(j)]), rt
            nvals += 1
    for c in COUNTERS:
        if rt.get(c) != mt.get(c):
            return False, "trailer %s is %s in the record (line %d), and %s here" % (c, rt.get(c), tl, mt.get(c)), rt
    return True, "%d samples, %d values, and its trailer" % (len(rs), nvals), rt


def closest_block(rec, k, names, lines):
    """Prints member k's closest approaches, each test particle to each active
    body at the samples (closest()), and returns them."""
    best = closest(rec, lines)
    print("  member %d: closest approach AT THE SAMPLES (%d of them, not between them), each test particle to each "
          "active body:" % (k, len(lines)))
    for j in range(rec.active[0], rec.N):
        print("    %s (body %d): %s" % (names[j], j, "; ".join(
            "%s %.6e at sample %d (step %d)" % ((names[i],) + best[(j, i)]) if (j, i) in best else
            "%s: no sample line to measure" % names[i] for i in range(rec.active[0]))))
    return best


def closest(rec, lines):
    """Each test particle's closest approach to each active body AT THE
    SAMPLES, from the record's own x y z columns, in binary64."""
    nact = rec.active[0]
    best = {}
    for lineno, line in lines:
        toks = line.split()
        r = rec.cols.split(toks)
        if isinstance(r, str):
            continue
        _, k, step = r
        pos = []
        for b in range(rec.N):
            v = [exact(toks[rec.cols.body_tok(b, c)]) for c in ("x", "y", "z")]
            pos.append(tuple(float(x[0]) if x is not None else math.nan for x in v))
        for j in range(nact, rec.N):
            for i in range(nact):
                d = math.dist(pos[j], pos[i])
                if (j, i) not in best or d < best[(j, i)][0]:
                    best[(j, i)] = (d, k, step)
    return best


def p3_replay():
    """Whether P3's replay is in this tree: ("absent", None), ("wired", its
    extract_particle module), or ("broken", why it cannot be called)."""
    present = [t for t in P3_TOOLS if os.path.exists(os.path.join(here, t))]
    if not present:
        return "absent", None
    if "extract_particle.py" not in present:
        return "broken", "tools/%s is in this tree without tools/extract_particle.py" % present[0]
    try:
        import importlib
        xp = importlib.import_module("extract_particle")
    except Exception as e:  # noqa: BLE001 - any failure to import is the finding
        return "broken", "tools/extract_particle.py does not import: %s" % e
    if not callable(getattr(xp, "extract", None)) or not isinstance(getattr(xp, "Refused", None), type):
        return "broken", "tools/extract_particle.py has no extract() and Refused to call"
    return "wired", xp


def particle_compare(rec, s, lines, out, bodies, names):
    """A particle's replay against the RECORD's lines for its member: every
    sample's fields before the bodies, each sub-system body against the
    member's body it is (`bodies`: the massive ones, then the particle),
    and the exact-time pair, value for value; then the trailer's counters.
    Returns (identical, what)."""
    sub_cols = None
    for l in out.splitlines():
        if l.startswith("# columns:"):
            sub_cols = Columns(l, len(bodies))
            break
    if sub_cols is None or not sub_cols.same_layout(rec.cols):
        return False, "the replay's `# columns:` line is not the record's layout"
    try:
        rs, _ = parse_samples("\n".join(l for _, l in lines), s if rec.kind == "ensemble" else None)
        ps, pt = parse_samples(out)
    except (ValueError, IndexError) as e:
        return False, "a sample line does not parse as hex floats: %s" % e
    if len(rs) != len(ps):
        return False, "%d samples in the record, %d in the replay" % (len(rs), len(ps))
    ns, per, npair = len(rec.cols.scalars) - 2, len(rec.cols.per), len(rec.cols.pair)
    pairs = [(j, j) for j in range(ns)]
    for i, b in enumerate(bodies):
        pairs += [(ns + b * per + f, ns + i * per + f) for f in range(per)]
    pairs += [(ns + rec.N * per + q, ns + len(bodies) * per + q) for q in range(npair)]
    sub_lines = [l for l in out.splitlines() if l.startswith("sample ")]
    nvals = 0
    for i, ((rk, rstep, rv), (pk, pstep, pv)) in enumerate(zip(rs, ps)):
        lineno, line = lines[i]
        where = "record line %d (%ssample %d, step %d)" % (lineno, "system %d " % s if rec.kind == "ensemble" else "",
                                                           rk, rstep)
        if (rk, rstep) != (pk, pstep):
            return False, "%s: it is sample %d at step %d in the replay" % (where, pk, pstep)
        if len(pv) != ns + len(bodies) * per + npair:
            return False, "%s: %d values in the replay, where %d bodies make %d" % (
                where, len(pv), len(bodies), ns + len(bodies) * per + npair)
        for rj, pj in pairs:
            if rv[rj] != pv[pj]:
                return False, "%s: %s: record %s, replay %s" % (
                    where, rec.cols.value_name(rj, names), line.split()[rec.cols.value_tok(rj)],
                    sub_lines[i].split()[sub_cols.value_tok(pj)])
            nvals += 1
    rt = rec.trailer[s][2]
    for c in COUNTERS:
        if rt.get(c) != pt.get(c):
            return False, "trailer %s is %s in the record, and %s in the replay" % (c, rt.get(c), pt.get(c))
    return True, "%d samples, %d values and the trailer's counters" % (len(rs), nvals)


def sample_particles(rec, info, args, work, state, xp, producer=None):
    """THE P3 x P4 SEAM - the lead's (ROUND4.md, "What the integrator keeps"),
    and the producer's own files (ROUND5.md, P3).

    Two test particles of each member in `info`: its least closest approach
    at the samples, and one by the seed. Each is written alone with the
    member's massive bodies by P3's extract(), from the member's attempts
    record and schedule - the producer's own where --producer-files gave a
    pair (`producer`), and otherwise this tree's rerun's - replayed on the
    engine those files are keyed to, and held to the record's lines
    (particle_compare). Each replay's record is kept beside its sub-system,
    <work>/particles/<stem>.rec. Returns (whole, failures): whole is False
    when P3's replay is in the tree and cannot be called."""
    if state == "absent":
        print("particles: not sampled (no particle replay in this tree)")
        return True, []
    if state == "broken":
        print("particles: FAIL: P3's replay is in this tree and cannot be called: %s" % xp)
        return False, []
    if not rec.active:
        print("particles: none to sample (no `# active` line: every body is active)")
        return True, []
    nact = rec.active[0]
    fails, n, members_done = [], 0, 0
    n_prod, m_prod, m_alone = 0, 0, 0
    for s, k, names, lines, best, att, sched, src, engine, rerun in info:
        parts = list(range(nact, rec.N))
        if not parts:
            continue
        members_done += 1
        if src == "producer":
            m_prod += 1
            m_alone += not rerun
            source = "the producer's own attempts record and schedule, %s and %s" % (att, sched)
        else:
            source = "this tree's rerun's attempts record and schedule%s" % (
                "" if producer is None else ", the producer's pair not being there")
        def least(j):
            return min((best[(j, i)][0] for i in range(nact) if (j, i) in best), default=math.inf)
        near = min(parts, key=lambda j: (least(j), j))
        rest = [j for j in parts if j != near]
        picks = [(near, "its least closest approach at the samples")]
        if rest:
            seeded = min(rest, key=lambda j: hashlib.sha256(b"cft-rebound compare_run particle %d %d %d"
                                                             % (args.seed, k, j)).digest())
            picks.append((seeded, "by the seed"))
        # the run the files are keyed to: a rerun's own engine, or for the
        # producer's files the engine their key names (compare_run(), below)
        base = rerun_cmd(rec, args, s, engine=engine)[1:]
        run_args = []
        i = 0
        while i < len(base):
            if base[i] in ("--problem", "--member"):
                i += 2
                continue
            run_args.append(base[i])
            i += 1
        member_arg = s if rec.kind == "ensemble" else rec.member
        for j, why in picks:
            n += 1
            n_prod += src == "producer"
            label = "  member %d: particle %s (body %d, %s)" % (k, names[j], j, why)
            try:
                res = xp.extract(args.exe, args.problem, att, sched, j, os.path.join(work, "particles"), run_args,
                                 member=member_arg)
            except xp.Refused as e:
                print("%s: FAIL: P3's extract_particle refused %s: %s" % (label, source, e))
                fails.append((s, "particle %s refused" % names[j]))
                continue
            rc, out, err = run([args.exe] + res["args"])
            with open(os.path.splitext(res["problem"])[0] + ".rec", "w", newline="\n") as f:
                f.write(out)
            if rc != 0:
                print("%s: FAIL: its replay from %s did not run (rc %d): %s" % (label, source, rc,
                                                                            (err or out).strip()[-300:]))
                fails.append((s, "particle %s's replay did not run" % names[j]))
                continue
            ok, what = particle_compare(rec, s, lines, out, res["bodies"], names)
            pt = kv(next((l[2:].split() for l in out.splitlines() if l.startswith("# steps_done=")), []))
            print("%s: replayed alone with the member's %d massive bodies%s through tools/extract_particle.py, from %s, "
                  "and held to the record: %s: %s; its replay's counters: steps %s, rejected %s, mean passes %s, max "
                  "passes %s, cap hits %s; its closest approach AT THE SAMPLES: %s" % (
                      label, nact, " and its zonal line" if " zonal=" in res["sub_key"] else "", source,
                      "IDENTICAL" if ok else "FAIL", what, pt.get("steps_done"), pt.get("steps_rejected"),
                      pt.get("mean_pc_iterations"), pt.get("max_pc_iterations"), pt.get("iterations_max_exceeded"),
                      "; ".join("%s %.6e at sample %d (step %d)" % ((names[i],) + best[(j, i)]) if (j, i) in best
                                else "%s: no sample line to measure" % names[i] for i in range(nact))))
            if not ok:
                fails.append((s, "particle %s" % names[j]))
    if producer is None:
        print("particles: %d sampled in %d member(s), each replayed alone from its member's attempts record and "
              "schedule as this tree reran it, and held to the record's own lines: %d identical"
              % (n, members_done, n - len(fails)))
    else:
        print("particles: %d sampled in %d member(s), each replayed alone and held to the record's own lines: %d from "
              "the producer's own attempts records and schedules, in %d member(s), %d of them not rerun; %d from this "
              "tree's reruns', where the producer's pair is not there; %d identical"
              % (n, members_done, n_prod, m_prod, m_alone, n - n_prod, n - len(fails)))
    return True, fails


def format_diff(rec, args, s, rerun_path, dt_path, work, rejected):
    """The member at a second format, on its own accepted steps."""
    f1, f2 = rec.program["format"], args.format_diff
    rc, out, err = run(rerun_cmd(rec, args, s, fmt=f2, dt_file=dt_path))
    if rc != 0:
        return False, "the replay at %s failed (rc %d): %s" % (f2, rc, (err or out).strip()[-300:])
    p2 = os.path.join(work, "member%d_%s_dtfile.rec" % (s, f2))
    with open(p2, "w", newline="\n") as f:
        f.write(out)
    head = header_args("\n".join(l for l in out.splitlines() if l.startswith("# program=") or l.startswith("# dt0=")))
    trail = kv(next((l[2:].split() for l in out.splitlines() if l.startswith("# steps_done=")), []))
    with open(dt_path) as f:
        dts = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    steps = int(rec.program["steps"])
    eps = exact(head.get("epsilon", "x"))
    if head.get("format") != f2 or "dt_file" not in head or eps is None or eps[0] != 0 or len(dts) != steps + 1 \
            or trail.get("steps_done") != rec.program["steps"] or trail.get("steps_rejected") != "0":
        return False, "the replay at %s is not on the recorded steps: format=%s dt_file=%s epsilon=%s, %d lines in " \
                      "the dt file for %d steps, steps_done=%s steps_rejected=%s" % (
                          f2, head.get("format"), head.get("dt_file"), head.get("epsilon"), len(dts), steps,
                          trail.get("steps_done"), trail.get("steps_rejected"))
    narrow, wide = (rerun_path, p2) if FORMATS.index(f1) < FORMATS.index(f2) else (p2, rerun_path)
    cmd = [sys.executable, os.path.join(here, "compare_formats.py"), narrow, wide, "--quiet"]
    # The dt file's first step is the first ACCEPTED one: the record's dt0
    # exactly when the first attempt was accepted, and smaller when it was
    # rejected. That is what tells the two kinds of rejection apart below.
    first_rejected = exact(dts[0]) != exact(rec.dt0["dt0"])
    if first_rejected:
        # the replay's dt0 is then not the record's; its steps are still
        # the dt file's, so there is no time offset to warn of
        cmd.append("--allow-dt-mismatch")
    rc, out, err = run(cmd)
    m = FLOOR.search(out)
    if rc != 0 or not m:
        return False, "tools/compare_formats.py did not measure it (rc %d): %s" % (rc, (out + err).strip()[-300:])
    note = ""
    if first_rejected:
        note = "; this member's first attempt was rejected, before any accepted step, which leaves state the replay " \
               "at %s does not take, so the difference includes that" % f2
    elif rejected:
        note = "; this member's %d rejected attempt%s all followed an accepted step, which leaves no state bit, so " \
               "the replay at %s follows its path" % (rejected, "" if rejected == 1 else "s", f2)
    return True, "%s against %s on the same %d accepted steps (--dt-out, then --dt-file with --epsilon 0): " \
                 "max |dx|/L %s, max |dE/E| %s over %s common samples (tools/compare_formats.py): %s%s" % (
                     m.group(1), m.group(2), steps, m.group(6), m.group(5), m.group(3), FORMAT_DIFF_SCOPE, note)


def built(path):
    return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(path)))


def producer_plan(rec, args, sample, ranges, picks, rej, state, bad):
    """THE PRODUCER'S OWN FILES (--producer-files): the members whose
    particles are sampled, where each one's pair is, and how and why each
    one's own lines are held. `sample` is the member sample less `bad`, the
    members the record does not hold whole, which are FAILed already and
    sampled for nothing; `ranges` and `picks` are its strata (None for a
    one-system record); `rej` its members with a rejected step.
    Returns ({member: (attempts record, schedule)} for each pair found, the
    members whose particles are sampled, the members that fall back to a
    rerun). Nothing is run."""
    d = args.producer_files
    if state == "absent":
        raise Refused("--producer-files %s: this tree has no particle replay (tools/extract_particle.py) to take a "
                      "producer's files" % d)
    if not os.path.isdir(d):
        raise Refused("--producer-files %s is not a directory" % d)
    if args.keep and os.path.normcase(os.path.abspath(args.keep)) == os.path.normcase(os.path.abspath(d)):
        raise Refused("--keep %s is the --producer-files directory, where this tree's reruns would write their own "
                      "attempts records and schedules over the producer's, under the same names" % args.keep)
    if not rec.active:
        print("particles: the producer's files in %s are not read: the record has no `# active` line, so it has no "
              "test particle" % d)
        return {}, [], []
    if rec.kind == "ensemble":
        strata = args.strata if args.particle_strata is None else args.particle_strata
        pranges, ppicks = particle_plan(rec.E, strata, args.seed)
        members = sorted(set(sample) | set(m for m in ppicks if m not in bad))
        print("particles: from the producer's own files in %s: %s; and the member sample's: %s; %d of %d members"
              % (d, "%d strata by member index (%s), seed %d: %s" % (
                  len(pranges), ", ".join("%d-%d" % (lo, hi - 1) for lo, hi in pranges), args.seed,
                  ", ".join(str(m) for m in ppicks)) if pranges else "no strata of their own (--particle-strata 0)",
                 ", ".join(str(m) for m in sample) or "none", len(members), rec.E))
    else:
        members = list(sample)
        print("particles: from the producer's own files in %s: the record's one system" % d)
    producer, fallback = {}, []
    for s in members:
        k = rec.label_of(s)
        att, sched = producer_paths(d, k)
        lost = [x for x in (att, sched) if not os.path.isfile(x)]
        if s in sample:
            how = "its lines held by a rerun here: the member sample holds it (%s)" % (
                "the record's one system" if ranges is None else "; ".join(
                    ["the pick of stratum %d-%d" % (lo, hi - 1) for (lo, hi), m in zip(ranges, picks) if m == s] +
                    ["its trailer shows %d rejected step(s)" % rec_rejected(rec, s)] * (s in rej)))
        else:
            lo, hi = next((lo, hi) for lo, hi in ranges if lo <= s < hi)
            how = "its lines held by %s: the member sample does not hold it (stratum %d-%d gives member %d, and its " \
                  "trailer shows no rejected step)%s" % (
                      "a rerun here" if lost else "the record alone, not rerun", lo, hi - 1, picks[ranges.index((lo, hi))],
                      ", and its particles need one" if lost else "")
        if lost:
            fallback.append(s)
            src = "the producer's pair is not there (no %s): its particles fall back to today's rerun, which writes " \
                  "the attempts record and schedule they replay from" % " and no ".join(lost)
        else:
            producer[s] = (att, sched)
            src = "its particles from the producer's own %s and %s" % (att, sched)
        print("  member %d: %s; %s" % (k, how, src))
    return producer, members, fallback


def rec_rejected(rec, s):
    t = rec.trailer.get(s)
    return int(t[2]["steps_rejected"]) if t and t[2].get("steps_rejected", "").isdigit() else 0


def compare_run(args):
    print("compare_run: %s against this tree" % args.record)
    if not os.path.isfile(args.exe):
        raise Refused("this tree's ias15_cft is not at %s" % args.exe)
    print("  this tree's ias15_cft: %s, built %s" % (args.exe, built(args.exe)))
    if not os.path.isfile(args.record):
        raise Refused("there is no record %s" % args.record)
    prob = Problem(args.problem)
    lines = numbered(args.record)
    rec = read_header(lines, args.record)
    p = rec.program
    if p["format"] != args.format:
        raise Refused("the record is format=%s, and this comparison holds format=%s (--format): a record of another "
                      "format is not comparable bit for bit" % (p["format"], args.format))
    if p["arith"] != args.arith:
        raise Refused("the record is arith=%s, and this comparison holds arith=%s (--arith): a record of another "
                      "arithmetic is not comparable bit for bit" % (p["arith"], args.arith))
    if args.format_diff == p["format"]:
        raise Refused("--format-diff %s is the record's own format" % args.format_diff)
    hold_problem(rec, prob)
    scan(rec, lines, prob)
    print("record: %s, %s; format=%s arith=%s cs=%s engine=%s pc_tol_shift=%s max_iter=%s steps=%s sample=%s; "
          "problem=%s" % (
              "an ensemble of E=%d systems of N=%d bodies" % (rec.E, rec.N) if rec.kind == "ensemble" else
              "one system of N=%d bodies%s" % (rec.N, ", member %d of its problem" % rec.member
                                               if rec.member is not None else ""),
              "active %d with %d pairs a force evaluation, a system" % (rec.active[0], rec.active[2])
              if rec.active else "every body active",
              p["format"], p["arith"], p["cs"], p["engine"], p["pc_tol_shift"], p["max_iter"], p["steps"],
              p["sample"], p["problem"]))
    print("record: written by program=%s impl=%s backend=%s" % (p.get("program", "?"), p.get("impl", "?"),
                                                                 p.get("backend", "?")))
    print("problem: %s is the record's: its name, G, %d system(s) of %d bodies, every body's name and mass, %s, and "
          "every system's initial state against the record's sample 0" % (
              args.problem, rec.E, rec.N, "active %s" % prob.active if prob.active else "no active line"))
    bad, gone, whole = hold_structure(rec)
    for s in sorted(bad):
        print("member %d (%s): FAIL: %s the record: %s" % (rec.label_of(s), rec.name_of(s, prob),
                                                            "missing from" if gone[s] else "incomplete in",
                                                            "; ".join(bad[s])))
    for w in whole:
        print("record: FAIL: %s" % w)

    rejected = {s: int(t[2]["steps_rejected"]) for s, t in rec.trailer.items()
                if s not in bad and t[2].get("steps_rejected", "").isdigit()}
    ranges = picks = None
    rej = []
    if rec.kind == "ensemble":
        rej = sorted(s for s in rejected if rejected[s] > 0)
        ranges, picks, members = plan(rec.E, args.strata, args.seed, rej)
        print("sample: seed %d, %d strata by member index (%s): %s; with a rejected step: %s; %d of %d members" % (
            args.seed, len(ranges), ", ".join("%d-%d" % (lo, hi - 1) for lo, hi in ranges),
            ", ".join(str(m) for m in picks), ", ".join(str(s) for s in rej) or "none", len(members), rec.E))
        print(SAMPLE_LIMIT)
    else:
        members = [0]
        print("sample: the record's one system%s" % (", member %d of its problem" % rec.member
                                                     if rec.member is not None else ""))
    todo = [s for s in members if s not in bad]
    state, xp = p3_replay()
    # THE PRODUCER'S OWN FILES: which members' particles are sampled, from
    # whose files, and how each one's lines are held (producer_plan())
    producer, pmembers, fallback = None, [], []
    if args.producer_files is not None:
        producer, pmembers, fallback = producer_plan(rec, args, todo, ranges, picks, rej, state, bad)
    extra = [s for s in fallback if s not in todo]                   # rerun only for the files
    alone = [s for s in pmembers if s not in todo and s not in extra]   # the producer's pair, and no rerun
    if args.plan:
        print("RESULT: PLAN: %d member(s) would be rerun%s%s" % (
            len(members) + len(extra), "; %d more sampled for particles alone, from the producer's files" % len(alone)
            if alone else "", "; and the record FAILs, above" if bad or whole else ""))
        return 1 if bad or whole else 0
    # each pair held to its member's run by the program's own key, before
    # anything runs: the engine its key names is the record's when this tree
    # has it, and otherwise the one the caller reruns on
    key_engine = None
    if producer and state == "wired":
        key_engine = p["engine"] if p["engine"] in ENGINES else args.rerun_engine
        if key_engine is None:
            raise Refused("the record ran engine=%s, which this tree does not have, so its producer's files are keyed "
                          "to the engine its rerun takes, and no --rerun-engine names one" % p["engine"])
        for s in sorted(producer):
            k, want = rec.label_of(s), run_key(rec, args, s, key_engine)
            for (what, _), path, magic in zip(PRODUCER_FILES, producer[s], (xp.ATT_MAGIC, SCHED_MAGIC)):
                hold_key(k, what, path, file_key(k, what, path, magic), want)
        print("particles: %d producer's pair(s), each keyed to its member's run: its run key, both files', is the one "
              "this tree's ias15_cft --run-key gives the member's arguments, at engine=%s" % (len(producer), key_engine))

    work = args.keep or tempfile.mkdtemp(prefix="compare_run_")
    os.makedirs(work, exist_ok=True)
    failed = set(bad)
    same = 0
    try:
        got = extract(args.record, rec, todo + extra + alone)
        engine = args.rerun_engine or p["engine"]
        if engine != p["engine"]:
            print("rerun: on the %s engine, where the record ran %s (--rerun-engine)" % (engine, p["engine"]))
        want_particles = state == "wired" and bool(rec.active)
        info = []
        for s in todo + extra:
            name, k = rec.name_of(s, prob), rec.label_of(s)
            # this rerun writes the files its particles replay from, unless
            # the producer's own pair is there for them
            files_here = want_particles and (producer is None or s in fallback)
            dt_path = os.path.join(work, "member%d.dt" % k) if args.format_diff else None
            att = os.path.join(work, "member%d.att" % k) if files_here else None
            sched = os.path.join(work, "member%d.sched" % k) if files_here else None
            rc, out, err = run(rerun_cmd(rec, args, s, dt_out=dt_path, att_out=att, sched_out=sched))
            if rc != 0:
                raise Refused("this tree's ias15_cft refused member %d's rerun with the record's arguments (rc %d): %s"
                              % (k, rc, (err or out).strip()[-400:]))
            mine = read_header(iter(enumerate(out.splitlines(), 1)), "the rerun of member %d" % k)
            hold_rerun(rec, s, mine, engine)
            rerun_path = os.path.join(work, "member%d_%s.rec" % (k, p["format"]))
            with open(rerun_path, "w", newline="\n") as f:
                f.write(out)
            names = [n for n, _ in rec.bodies[s]]
            ok, what, rt = compare(rec, s, got[s], mine, out, names)
            counters = "steps %s, rejected %s, mean passes %s, max passes %s, cap hits %s" % (
                rt.get("steps_done"), rt.get("steps_rejected"), rt.get("mean_pc_iterations"),
                rt.get("max_pc_iterations"), rt.get("iterations_max_exceeded"))
            if ok:
                same += 1
            else:
                failed.add(s)
            print("member %d (%s): %s: %s; %s" % (k, name, "IDENTICAL" if ok else "FAIL", what, counters))
            if not rec.active:
                print("  member %d: no test particles (no `# active` line: every body is active)" % k)
            else:
                best = closest_block(rec, k, names, got[s])
                if want_particles and producer is not None and s in producer:
                    info.append((s, k, names, got[s], best, producer[s][0], producer[s][1], "producer", key_engine, True))
                elif want_particles and (producer is None or s in pmembers):
                    info.append((s, k, names, got[s], best, att, sched, "rerun", engine, True))
            if args.format_diff:
                fok, fmsg = format_diff(rec, args, s, rerun_path, dt_path, work, int(rt.get("steps_rejected", "0")))
                print("  member %d: format difference: %s%s" % (k, "" if fok else "FAIL: ", fmsg))
                if not fok:
                    failed.add(s)
        # the members sampled for their particles alone: not rerun, their
        # lines held by the record and by those particles' replays, which
        # hold every massive body and the exact-time pair
        for s in alone:
            name, k = rec.name_of(s, prob), rec.label_of(s)
            names = [n for n, _ in rec.bodies[s]]
            rt = rec.trailer[s][2]
            print("member %d (%s): NOT RERUN: its lines are held by the record alone, and by its particles' replays from "
                  "the producer's own files; its trailer: steps %s, rejected %s, mean passes %s, max passes %s, cap "
                  "hits %s" % (k, name, rt.get("steps_done"), rt.get("steps_rejected"), rt.get("mean_pc_iterations"),
                               rt.get("max_pc_iterations"), rt.get("iterations_max_exceeded")))
            best = closest_block(rec, k, names, got[s])
            if want_particles:
                info.append((s, k, names, got[s], best, producer[s][0], producer[s][1], "producer", key_engine, False))
        particles_whole, pfails = sample_particles(rec, info, args, work, state, xp, producer)
    finally:
        if not args.keep:
            shutil.rmtree(work, ignore_errors=True)
    nrerun = len(members) + len(extra)
    if failed or whole or not particles_whole or pfails:
        why = []
        if failed:
            why.append("FAIL: " + ", ".join("member %d" % rec.label_of(s) for s in sorted(failed)))
        if whole:
            why.append("FAIL: the record as a whole")
        if pfails:
            why.append("FAIL: " + ", ".join("member %d's %s" % (rec.label_of(s), what) for s, what in pfails))
        if not particles_whole:
            why.append("P3's replay is in this tree and cannot be called")
        print("RESULT: FAIL: %d of %d sampled member(s) identical; %s" % (same, nrerun, "; ".join(why)))
        return 1
    print("RESULT: PASS: %d sampled member(s) of %d, every one identical to this tree's rerun, its counters beside it%s%s"
          % (same, rec.E, "; every sampled particle's replay identical to the record" if info else "",
             "; and %d member(s) not rerun, each held by the record and its particles' replays from the producer's own "
             "files" % len(alone) if alone else ""))
    return 0


def main():
    exe = ".exe" if os.name == "nt" else ""
    ap = argparse.ArgumentParser(description="A run's record against this tree: a sample of its members rerun here, "
                                             "bit for bit, with their physics counters.")
    ap.add_argument("record")
    ap.add_argument("--problem", required=True, help="the problem file the record ran")
    ap.add_argument("--format", required=True, choices=FORMATS, help="the format of the run this comparison holds")
    ap.add_argument("--arith", required=True, choices=("rebound", "fma"), help="the arithmetic of the run it holds")
    ap.add_argument("--exe", default=os.path.join(root, "build", "ias15_cft" + exe), help="this tree's ias15_cft")
    ap.add_argument("--programs", default=os.path.join(root, "programs", "out"), help="for a program-engine rerun")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--strata", type=int, default=8, help="contiguous ranges of member index, one member from each")
    ap.add_argument("--format-diff", choices=FORMATS, help="also rerun each sampled member at this format, on its steps")
    ap.add_argument("--rerun-engine", choices=("loop", "program"),
                    help="rerun on this engine, not the record's (the two give the same bits in the FMA form)")
    ap.add_argument("--plan", action="store_true", help="say what would be rerun, and stop")
    ap.add_argument("--keep", help="keep the reruns' records in this directory")
    ap.add_argument("--producer-files", metavar="DIR",
                    help="the producer's own attempts record and corrector schedule for member k, DIR/member<k>.att "
                         "and DIR/member<k>.sched: particles are replayed from them, with no rerun of their member")
    ap.add_argument("--particle-strata", type=int, metavar="P",
                    help="with --producer-files, sample particles in one more member from each of P ranges of member "
                         "index, besides the member sample's (default: --strata's)")
    args = ap.parse_args()
    if args.particle_strata is not None and (args.particle_strata < 0 or args.producer_files is None):
        ap.error("--particle-strata takes a count of 0 or more, with --producer-files")
    args.problem = os.path.abspath(args.problem)
    # absolute, as Windows' CreateProcess runs no relative path written with
    # forward slashes (a traceback otherwise, measured 2026-11-12)
    args.exe = os.path.abspath(args.exe)
    try:
        return compare_run(args)
    except Refused as e:
        print("RESULT: REFUSED: %s" % e)
        return 2


if __name__ == "__main__":
    sys.exit(main())
