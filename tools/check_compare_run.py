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
#   2. The particles line, exactly: no particle replay is in this tree; and
#      beside the sample, its limit, exactly: it trusts the trailers it is
#      testing.
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
#      member, the line and the column;
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
#   i. P3's tools in the tree and the particle place empty: a copy of the
#      driver beside an empty stand-in for each of P3's two files in turn
#      says NOT WIRED, ends RESULT: FAIL and exits 1, though its member is
#      IDENTICAL; beside neither, it passes with the line in case 2.
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
NOT_WIRED = "particles: NOT SAMPLED, and NOT WIRED: P3's tools/%s is in this tree"
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


def write_waves(path, name="cr_waves", tp_mass=None, ulp=None, active=True):
    """The waves. tp_mass=(k, j, m) gives system k's particle j mass m;
    ulp=(k, j) moves system k's particle j's x by one ulp."""
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
        lines.append("active 2")
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true", help="binary64, and short binary128 and FMA records")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, "ias15_cft" + exe)
    ref = os.path.join(args.build, "ias15_ref" + exe)
    progs = os.path.join(root, "programs", "out")
    py = sys.executable
    driver = os.path.join(here, "compare_run.py")
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

    def record(name, fmt, steps, sample, arith="rebound", engine="loop"):
        path = os.path.join(scratch, name)
        cmd = [cft, "--format", fmt, "--problem", waves, "--dt", DT0, "--epsilon", "1e-9", "--steps", str(steps),
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
                    n += 1
                    worst = max(worst, abs(mpmath.mpf(val) - best) / best)
                    if dists[int(at)] != best:
                        at_wrong.append("member %d tp%d body %d" % (k, j - 2, i))
        verdict(not missing and not at_wrong and n == 2 * NTP * len(smp["members"]) and worst < 1e-6,
                "%s: %d closest approaches (each test particle to each active body, at the samples) agree with this "
                "gate's mpmath from tools/oracle.py's parse to %s relative, at the samples they name%s" % (
                    label, n, mpmath.nstr(worst, 3), "; " + ", ".join(missing + at_wrong) if missing or at_wrong else ""))
        verdict(out.count("\n" + PARTICLES + "\n") == 1, "%s: the particles line, exactly: '%s'" % (label, PARTICLES))
        verdict(out.count("\n" + SAMPLE_LIMIT + "\n") == 1 and out.find(SAMPLE_LIMIT) > out.find("\nsample: seed "),
                "%s: beside the sample, its limit, exactly: '%s'" % (label, SAMPLE_LIMIT))
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
        if fmt == "fp64":
            smp64, text64, rec64 = smp, text, rec
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
            out.count("\n" + PARTICLES + "\n") == 1,
            "a --member %d record of the waves, held as one system: PASS, member %d IDENTICAL with its trailer's "
            "counters (%s), its closest approaches, and the particles line (rc %d)" % (k6, k6, want, rc))
    rck, tk, _ = run([cft, "--format", "fp64", "--problem", kep, "--dt", "0.05", "--epsilon", "1e-9", "--steps", "20",
                      "--sample", "10", "--max-iter", "60", "--quiet"])
    reck = os.path.join(scratch, "kepler.rec")
    with open(reck, "w", newline="\n") as f:
        f.write(tk)
    rc, out = drive(reck, "fp64", problem=kep)
    lines = member_lines(out)
    verdict(rck == 0 and rc == 0 and "\nRESULT: PASS: " in out and sorted(lines) == [0] and
            lines[0][0] == "IDENTICAL" and "member 0 (kepler): IDENTICAL: 3 samples" in out and
            "  member 0: no test particles (no `# active` line: every body is active)" in out,
            "a single problem's record (kepler, every body active): PASS, IDENTICAL, and no test particles said (rc %d)"
            % rc)

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

    # i. P3's tools in the tree and the particle place still empty: a copy of
    # the driver, with what it imports and runs, beside a stand-in for each of
    # P3's two files in turn, and beside neither, on case 5's kepler record.
    # The stand-in is an empty file: only its presence is under test.
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
            verdict(rc == 1 and (NOT_WIRED % present[0]) in out and PARTICLES not in out and
                    result.startswith("RESULT: FAIL: ") and "particles are due and NOT sampled" in result and
                    lines.get(0, ("",))[0] == "IDENTICAL",
                    "i. the driver's copy beside a stand-in tools/%s: NOT WIRED, RESULT: FAIL and rc %d, though member 0 "
                    "is IDENTICAL" % (present[0], rc))
    shutil.rmtree(os.path.dirname(p3), ignore_errors=True)

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
