#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of a particle replayed out of its wave (ROUND4.md, P3): the
# attempts record (--attempts-out), its replay (--dt-file in that form),
# the run key (--run-key) and tools/extract_particle.py.
#
#   python tools/check_replay.py [--build build] [--exe ias15_cft]
#                                [--formats fp64,fp128,fp256] [--quick] [--keep DIR]
#   (`make check-quick` and `make check-light` pass --quick: binary64 only,
#    and shorter waves.)
#
# The waves are hand-made here: Saturn, six moons and eight test particles
# on inclined orbits, `active 7`. Every value is one correctly rounded
# binary64 operation away from exact decimals (the angles by the rational
# parametrisation, never a library's trig), so the file is the same bits on
# every IEEE host. The values are only roughly Saturn's: tools/make_ring.py
# (P2) writes the sourced ones. Three waves at every format:
#   A. the ring from dt0 = 600 s, which rejects nothing;
#   B. the ring from dt0 = 1e5 s, two orbits, far too large: the first
#      attempt is rejected, before any step is accepted;
#   C. the ring and a ninth test particle aimed past Janus, from dt0 = 100 s:
#      the attempts that meet Janus are rejected AFTER accepted steps, the
#      case whose next step is predicted from er and br.
# Each wave runs adaptively with --attempts-out and --record-iters. Then,
# for each of three sampled test particles (the first, the middle one and
# the last, which in C is the flyby), tools/extract_particle.py writes the
# particle alone with the seven massive bodies and re-keys the schedule,
# and the replay runs. All of that is on the loop engine in REBOUND's form.
# B and C, the two kinds of rejection, also run on the loop engine in the
# FMA form and on the program engine (which takes only that form): the
# wave, its rejections counted, and each particle's bits, counters and
# travel, items 1 and 2 and the first two counts of item 3. The step after
# a rejection is restored and predicted through each engine's own code;
# the records, the tool and the refusals are the same on all three.
#
# What it asserts, and what lets each assertion fail:
#   1. Bits. At every sample, the replay's particle AND every massive body,
#      six values each, and the line's step, t, dt_next, dt_last, energy and
#      exact-time pair, are the wave's, as hex. The massive bodies are held
#      because they match only if the particle adds exactly nothing to them:
#      a "test" particle that was not one, or a body out of order, would
#      move them while the particle's own bits could still look right.
#   2. Counters. The replay's steps_done, steps_rejected,
#      iterations_max_exceeded, mean_pc_iterations and max_pc_iterations are
#      the wave's. Bits alone cannot hold the passes: extra passes land on
#      the corrector's fixed point (docs/VALIDATION.md entry 41), and a
#      rejection after an accepted step leaves no state behind at all.
#   3. Physics, each counted and asserted. The wave's rejected attempts:
#      none in A, one or more before the first accepted step in B, one or
#      more after an accepted step in C; the attempts record agrees with the
#      wave's trailer and its schedule, attempt for attempt. The particle
#      moved: the angle it swept about Saturn, summed sample to sample. And
#      the records are load bearing, or the comparison would pass against
#      itself: in each wave, some sampled particle run alone takes other
#      bits on its own steps, and some sampled particle run alone on the
#      recorded steps takes other passes on its own exit test. The whole
#      wave, replayed through its attempts record alone on its own exit
#      test, gives its record lines and counters (every wave and form): the
#      run needs no schedule, and only a sub-system does.
#   4. The comparison sees what it must, on copies of the records:
#      - one accepted attempt's dt moved by one ulp: the replay keyed on the
#        true record is refused by name, and re-keyed it runs and the
#        particle's bits differ;
#      - B's rejected first attempt deleted from both records, re-keyed: the
#        bodies' values differ, the state and not only the dt_next of the
#        sample printed before the attempt;
#      - C's first rejection after an accepted step that has no sample line
#        just before it, deleted from both records and re-keyed: the bits
#        are the same (the step after a rejection is predicted from er and
#        br exactly as the accepted step before it would have predicted it)
#        and steps_rejected falls by one, so item 2 is what catches it. A
#        sample printed just before a rejected attempt carries that
#        attempt's own dt as its dt_next, so the attempt chosen has none, and
#        that precondition is asserted by name: another --sample cannot fail
#        the case with nothing wrong.
#   5. Refusals by name, at the first format: the records replayed for
#      another wave (a particle moved, the sampled one or another); the
#      sub-system with a massive body moved, two swapped, or its active line
#      changed, by the tool and by the program; a "particle" that is a
#      massive body; a wave with no active line; an attempts record with an
#      outcome flipped; the attempts form on an ensemble, in either
#      direction; a malformed line, a record cut short, a line after its
#      end, too few accepted attempts, and no --epsilon 0. And --run-key
#      prints exactly the run line --record-iters writes.
#   6. A wave that is one system of an ensemble, as P2's waves are, at the
#      first format: member 1 of two, recorded with --member 1, gives its
#      particles alone (one of them the body that tells the two waves apart),
#      and its records extracted for member 0 are refused by name.
#   7. A dt file through a pipe, which cannot be rewound, at the first format
#      (V3, 2026-11-12: a rewind there once dropped the C runtime's buffer).
#      Wave A's accepted-steps file (--dt-out), padded to fixed-width lines
#      past 4,096 bytes so that a lost buffer would cut it on a line
#      boundary and run silently, and its attempts record each replay
#      through a pipe to the lines and counters the same file gives as a
#      file. A Windows named pipe here; /dev/stdin on a POSIX host.
#
# The --light rules (CLAUDE.md): this gate has no scaling of its own. Every
# replay runs its wave's own --steps, never a count scaled apart from the
# record it replays, and C needs at least its first dozen attempts to reach
# Janus. --quick runs binary64 only, with shorter waves.
#
# Each run works in <build>/replay/<pid>; a passing run removes it, a
# failing one keeps it.

import argparse
import itertools
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
import extract_particle as xp  # noqa: E402

COUNTERS = ("steps_done", "steps_rejected", "iterations_max_exceeded", "mean_pc_iterations", "max_pc_iterations")
ok_all = True

# km, s, G = 1, masses as GM. Roughly Saturn's; the sourced values are P2's.
GM_SATURN = 37931187.0
MOONS = (("pan", 3.3e-4, 133584.0), ("atlas", 4.4e-4, 137670.0), ("prometheus", 1.074e-2, 139380.0),
         ("pandora", 9.24e-3, 141720.0), ("epimetheus", 3.51e-2, 151422.0), ("janus", 1.266e-1, 151472.0))
T_MOON = (0.45, -1.8, 0.95, 3.1, -0.42, 1.35)        # tan of half the phase
U_MOON = (5e-4, 1e-3, 1.5e-3, 2e-3, 2.5e-3, 3e-3)   # tan of half the inclination
NPART = 8
T_PART = (0.15, 1.2, -0.7, 2.6, -2.2, 0.55, -0.1, 4.0)
U_PART = (1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 1e-4, 2e-4, 3e-4)
FLYBY = (1500.0, 1.0, 20.0)   # km out from Janus, km/s toward it, km off its plane


def verdict(ok, text):
    global ok_all
    ok_all &= bool(ok)
    print("  %s %s" % ("ok  " if ok else "FAIL", text))
    return ok


def circular(r, t, u):
    """A circular orbit about Saturn of radius r, phase 2 atan(t), inclined
    2 atan(u): the rational parametrisation, so nothing but +-*/ and sqrt."""
    c, s = (1.0 - t * t) / (1.0 + t * t), 2.0 * t / (1.0 + t * t)
    ci, si = (1.0 - u * u) / (1.0 + u * u), 2.0 * u / (1.0 + u * u)
    v = math.sqrt(GM_SATURN / r)
    return (r * c, r * s * ci, r * s * si), (-v * s, v * c * ci, v * c * si)


def body_line(name, m, x, v):
    return "body %s %s" % (name, " ".join(float(q).hex() for q in (m,) + tuple(x) + tuple(v)))


def write_wave(path, flyby, name=None, active=True, edit=None):
    lines = ["# tools/check_replay.py: a hand-made ring wave, %s" % ("with a test particle aimed past Janus" if flyby
                                                                    else "Saturn, six moons and eight test particles"),
             "name " + (name or ("ringflyby" if flyby else "ringwave")), "G " + (1.0).hex()]
    if active:
        lines.append("active %d" % (1 + len(MOONS)))
    bodies = [body_line("saturn", GM_SATURN, (0.0, 0.0, 0.0), (0.0, 0.0, 0.0))]
    states = []
    for k, (nm, gm, a) in enumerate(MOONS):
        x, v = circular(a, T_MOON[k], U_MOON[k])
        states.append((x, v))
        bodies.append(body_line(nm, gm, x, v))
    for j in range(NPART):
        x, v = circular(130500.0 + 1500.0 * j, T_PART[j], U_PART[j])
        bodies.append(body_line("p%d" % j, 0.0, x, v))
    if flyby:
        (px, py, pz), (vx, vy, vz) = states[-1]
        d, w, b = FLYBY
        r = math.sqrt(px * px + py * py + pz * pz)
        ux, uy, uz = px / r, py / r, pz / r
        bodies.append(body_line("flyby", 0.0, (px + d * ux, py + d * uy, pz + d * uz + b),
                                (vx - w * ux, vy - w * uy, vz - w * uz)))
    if edit:
        bodies = edit(bodies)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines + bodies) + "\n")
    return len(bodies)


def run(cmd):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    return p.returncode, p.stdout, p.stderr


_pipes = itertools.count()


def run_piped(cmd, data):
    """cmd, then the path of a pipe that carries data and cannot be rewound:
    a Windows named pipe, or /dev/stdin on a POSIX host."""
    if os.name != "nt":
        p = subprocess.run(cmd + ["/dev/stdin"], input=data, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return p.returncode, p.stdout.decode(), p.stderr.decode()
    import _winapi
    name = r"\\.\pipe\cft-check-replay-%d-%d" % (os.getpid(), next(_pipes))
    # outbound (PIPE_ACCESS_OUTBOUND = 2), a byte stream that blocks, one instance
    h = _winapi.CreateNamedPipe(name, 0x00000002, 0, 1, 65536, 65536, 0, _winapi.NULL)

    def serve():
        try:
            try:
                _winapi.ConnectNamedPipe(h, _winapi.NULL)
            except OSError as e:
                if getattr(e, "winerror", None) != 535:        # ERROR_PIPE_CONNECTED: the program opened it first
                    raise
            _winapi.WriteFile(h, data)
        except OSError:
            pass
        finally:
            _winapi.CloseHandle(h)
    t = threading.Thread(target=serve, daemon=True)
    t.start()
    r = run(cmd + [name])
    if t.is_alive():            # the program never opened the pipe: release its server
        try:
            open(name, "rb").close()
        except OSError:
            pass
    t.join(timeout=30)
    return r


def sample_before(outc, k, sample):
    """Whether a record prints a sample line just before attempt k, which
    then carries attempt k's own dt as its dt_next: the first sample is
    printed before attempt 0, and each later one after the accepted attempt
    that brings the accepted count to a multiple of --sample."""
    return k == 0 or (outc[k - 1] == "a" and outc[:k].count("a") % sample == 0)


def kv(text):
    return dict(x.partition("=")[::2] for x in text.split())


def record(out):
    """The sample lines, as token lists, and the trailer's counters."""
    rec = [l.split() for l in out.splitlines() if l.startswith("sample ")]
    tail = next((kv(l[2:]) for l in out.splitlines() if l.startswith("# steps_done=")), {})
    return rec, tail


def outcomes(path):
    with open(path, encoding="utf-8") as f:
        return [l.split()[1] for l in f.read().splitlines() if l.strip() and not l.startswith("#")][:-1]


def passes(path):
    with open(path, encoding="utf-8") as f:
        return [l for l in f.read().splitlines() if l and not l.startswith("#")]


def edit_lines(src, dst, fn):
    with open(src, encoding="utf-8") as f:
        lines = f.read().splitlines()
    with open(dst, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(fn(lines)) + "\n")
    return dst


def body_lines(lines):
    """Indices of a record's attempt or schedule lines (not # lines)."""
    return [i for i, l in enumerate(lines) if l.strip() and not l.startswith("#")]


def ulp_up(h, digits):
    """A positive normal hex float moved up by one unit in the last place
    of a format with `digits` hex digits of fraction."""
    m = re.match(r"^0x1(?:\.([0-9a-f]*))?p([+-]\d+)$", h)
    if not m:
        raise ValueError("not a positive normal hex float: %s" % h)
    frac, exp = (m.group(1) or "").ljust(digits, "0"), int(m.group(2))
    f = int(frac, 16) + 1
    if f >> (4 * digits):
        f, exp = 0, exp + 1
    return "0x1.%0*xp%+d" % (digits, f, exp)


def compare(wave, sub, bodies):
    """Every sample's head and tail and each sub-system body against the
    wave's body it is, value for value. Returns (samples, values, mismatches)."""
    bad, values = [], 0
    if len(wave) != len(sub):
        bad.append("the wave has %d samples and the replay %d" % (len(wave), len(sub)))
    for w, s in zip(wave, sub):
        if w[:7] != s[:7] or w[-3:] != s[-3:]:
            bad.append("sample %s: its head or exact-time pair differs" % w[1])
        for i, b in enumerate(bodies):
            wv, sv = w[7 + 6 * b:13 + 6 * b], s[7 + 6 * i:13 + 6 * i]
            values += 6
            if wv != sv:
                bad.append("sample %s, wave body %d: %s where the wave has %s" % (w[1], b, " ".join(sv), " ".join(wv)))
    return len(sub), values, bad


def particle_differs(wave, sub, j, i):
    return sum(1 for w, s in zip(wave, sub) if w[7 + 6 * j:13 + 6 * j] != s[7 + 6 * i:13 + 6 * i]) + abs(len(wave) - len(sub))


def body_mismatches(wave, sub, bodies):
    """The (sample, body) pairs whose six values differ: the state, without
    the head fields."""
    return sum(1 for w, s in zip(wave, sub) for i, b in enumerate(bodies)
               if w[7 + 6 * b:13 + 6 * b] != s[7 + 6 * i:13 + 6 * i])


def travel(wave, j):
    """The angle body j swept about body 0, summed from each sample to the
    next, and the length of those chords, read back from the record's hex.
    Summed rather than first to last, so that a whole orbit does not read as
    none: a gate wave's samples are under half an orbit apart."""
    def rel(line):
        return [float.fromhex(a) - float.fromhex(b) for a, b in zip(line[7 + 6 * j:10 + 6 * j], line[7:10])]
    swept = chords = 0.0
    for l0, l1 in zip(wave, wave[1:]):
        r0, r1 = rel(l0), rel(l1)
        cosang = sum(p * q for p, q in zip(r0, r1)) / (math.sqrt(sum(p * p for p in r0)) * math.sqrt(sum(q * q for q in r1)))
        swept += math.acos(max(-1.0, min(1.0, cosang)))
        chords += math.sqrt(sum((p - q) ** 2 for p, q in zip(r1, r0)))
    return swept, chords


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(os.path.dirname(here), "build"))
    ap.add_argument("--exe", default="ias15_cft")
    ap.add_argument("--formats", default="fp64,fp128,fp256")
    ap.add_argument("--quick", action="store_true", help="binary64 only, shorter waves")
    ap.add_argument("--keep", default=None, help="directory for the waves and records, kept "
                    "(default <build>/replay/<pid>, removed after a pass)")
    ap.add_argument("--programs", default=os.path.join(os.path.dirname(here), "programs", "out"),
                    help="the program engine's images")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, args.exe + exe)
    formats = ["fp64"] if args.quick else args.formats.split(",")
    work = args.keep or os.path.join(args.build, "replay", str(os.getpid()))
    if not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work, exist_ok=True)
    print("check_replay: %s, built %s" % (cft, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(cft)))))

    waves = {"A": os.path.join(work, "ring.txt"), "C": os.path.join(work, "flyby.txt")}
    nring = write_wave(waves["A"], False)
    nfly = write_wave(waves["C"], True)
    waves["B"] = waves["A"]
    n_active = 1 + len(MOONS)
    cases = {"A": ("the ring from dt0 = 600 s (no rejection)", "600"),
             "B": ("the ring from dt0 = 1e5 s (a first dt far too large)", "1e5"),
             "C": ("the ring and a particle aimed past Janus, from dt0 = 100 s (rejections after accepted steps)", "100")}

    # The engine and form: REBOUND's on the loop engine, which every check
    # runs on, and the FMA form on each engine, where B and C replay too.
    forms = {"loop": ([], "the loop engine, REBOUND's form"),
             "loop-fma": (["--arith", "fma"], "the loop engine, the FMA form"),
             "program": (["--arith", "fma", "--engine", "program", "--programs", args.programs], "the program engine")}

    def run_args(fmt, case, form="loop"):
        steps = 24 if args.quick else (20 if fmt == "fp256" else 40)
        if case == "C":
            steps = 24 if args.quick or fmt == "fp256" else 30
        a = ["--format", fmt, "--dt", cases[case][1], "--epsilon", "1e-9", "--steps", str(steps),
             "--sample", str(steps // 4), "--quiet"] + forms[form][0]
        return a + (["--max-iter", "60"] if fmt == "fp256" else [])   # so its pass counts vary

    def sampled(case):
        """The first test particle, the middle one and the last (in C, the flyby)."""
        n = nfly if case == "C" else nring
        return [n_active, n_active + (n - n_active) // 2, n - 1]

    # One job a format, wave and form: the wave, then each particle's
    # extraction and replay; on the loop engine in REBOUND's form, each
    # particle's runs alone too, and the controls on copies of the records.
    def job(fcf):
        fmt, case, form = fcf
        res = {"fmt": fmt, "case": case, "form": form}
        d = os.path.join(work, "%s-%s-%s" % (fmt, case, form))
        os.makedirs(d, exist_ok=True)
        ra = run_args(fmt, case, form)
        att, sch, dtf = os.path.join(d, "wave.att"), os.path.join(d, "wave.sched"), os.path.join(d, "wave.dt")
        rc, out, err = run([cft, "--problem", waves[case]] + ra + ["--attempts-out", att, "--record-iters", sch]
                           + (["--dt-out", dtf] if case == "A" and form == "loop" else []))
        res["wave_rc"], res["wave_err"] = rc, err
        if rc:
            return res
        wrec, wtail = record(out)
        res.update(wrec=wrec, wtail=wtail, att=att, sch=sch, dtf=dtf, args=ra, outc=outcomes(att), passes=passes(sch),
                   parts={})
        # the whole wave through its attempts record alone, on its own exit test
        rc, out, err = run([cft] + xp.replay_args(ra, waves[case], att))
        res["whole"] = (rc, err) + record(out)
        for j in sampled(case):
            p = {}
            try:
                r = xp.extract(cft, waves[case], att, sch, j, d, ra)
            except xp.Refused as e:
                p["refused"] = str(e)
                res["parts"][j] = p
                continue
            p["x"] = r
            rc, out, err = run([cft] + r["args"])
            p["rc"], p["err"] = rc, err
            p["rec"], p["tail"] = record(out)
            res["parts"][j] = p
            if form != "loop":
                continue
            # alone: on its own steps, and on the recorded steps with its own exit test
            rc2, out2, _ = run([cft, "--problem", r["problem"]] + ra)
            p["own_rc"], (p["own"], _) = rc2, record(out2)
            own_sched = os.path.join(d, "own-%d.sched" % j)
            rc3, out3, _ = run([cft] + r["args"][:-2] + ["--record-iters", own_sched])
            p["ownx_rc"], (p["ownx"], p["ownx_tail"]) = rc3, record(out3)
            p["ownx_passes"] = passes(own_sched) if rc3 == 0 else []
        res["ctl"] = ctl = {}
        if form != "loop":
            return res
        j0 = sampled(case)[0]
        acc = [i for i, o in enumerate(res["outc"]) if o == "a"]
        if case == "A" and j0 in res["parts"] and "x" in res["parts"][j0]:
            # one accepted attempt's dt moved by one ulp, in a copy of the record
            k = acc[len(acc) // 2]
            digits = max(len(t.split("p")[0].split(".")[1]) for line in wrec for t in line
                         if t.startswith("0x1.") and "p" in t)
            moved = os.path.join(d, "moved.att")

            def move(lines):
                i = body_lines(lines)[k]
                dt, what = lines[i].split()
                lines[i] = "%s %s" % (ulp_up(dt, digits), what)
                ctl["moved"] = (k, dt, lines[i].split()[0])
                return lines
            edit_lines(att, moved, move)
            x0 = res["parts"][j0]["x"]
            args_true_key = [moved if a == att else a for a in x0["args"]]
            ctl["moved_refused"] = run([cft] + args_true_key)
            try:
                rk = xp.extract(cft, waves[case], moved, sch, j0, os.path.join(d, "moved"), ra)
                ctl["moved_rekeyed"] = run([cft] + rk["args"])
            except xp.Refused as e:
                ctl["moved_rekeyed"] = (3, "", "extract_particle: refused: %s" % e)
        if case in ("B", "C") and j0 in res["parts"] and "x" in res["parts"][j0]:
            # one rejected attempt deleted from copies of both records: in B
            # the first attempt; in C the first rejection after an accepted
            # step with no sample line just before it, since such a line
            # carries the attempt's own dt as its dt_next
            o = res["outc"]
            sample = int(ra[ra.index("--sample") + 1])
            late = [i for i, w in enumerate(o) if w == "r" and "a" in o[:i]]
            blocked = [i for i in late if sample_before(o, i, sample)]
            if case == "B":
                k = next((i for i, w in enumerate(o) if w == "r"), None)
            else:
                k = next((i for i in late if not sample_before(o, i, sample)), None)
                ctl["deleted_pre"] = (late, k, sample, blocked)

            def delete(kk, tag):
                datt, dsch = os.path.join(d, tag + ".att"), os.path.join(d, tag + ".sched")
                edit_lines(att, datt, lambda L: [l for i, l in enumerate(L) if i != body_lines(L)[kk]])
                edit_lines(sch, dsch, lambda L: [l for i, l in enumerate(L) if i != body_lines(L)[kk]])
                try:
                    rk = xp.extract(cft, waves[case], datt, dsch, j0, os.path.join(d, tag), ra)
                    return kk, run([cft] + rk["args"])
                except xp.Refused as e:
                    return kk, (3, "", "extract_particle: refused: %s" % e)
            if k is not None:
                ctl["deleted"] = delete(k, "deleted")
            # and in C, where one exists, a late rejection WITH a sample line
            # just before it: only that line's dt_next may move
            if case == "C" and blocked:
                with open(att, encoding="utf-8") as f:
                    adts = [l.split()[0] for l in f.read().splitlines() if l.strip() and not l.startswith("#")]
                ctl["deleted_b"] = delete(blocked[0], "deleted-b") + (adts,)
        res["ctl"] = ctl
        return res

    order = [("A", "loop"), ("B", "loop"), ("C", "loop"), ("B", "loop-fma"), ("C", "loop-fma"),
             ("B", "program"), ("C", "program")]
    jobs = [(fmt, case, form) for fmt in reversed(formats) for case, form in order]   # the widest first: the long poles
    with ThreadPoolExecutor(max_workers=1 if os.environ.get("CFT_REBOUND_ARTIFACT") else 4) as pool:
        results = dict(zip(jobs, pool.map(job, jobs)))

    for fmt in formats:
        for case, form in order:
            res = results[(fmt, case, form)]
            print("check_replay at %s: %s, on %s" % (fmt, cases[case][0], forms[form][1]))
            if res["wave_rc"]:
                verdict(False, "the wave did not run (rc %d): %s" % (res["wave_rc"], res["wave_err"].strip()[-300:]))
                continue
            o, wt, wrec = res["outc"], res["wtail"], res["wrec"]
            rej = o.count("r")
            first = o.index("a") if "a" in o else len(o)
            early, late = sum(1 for w in o[:first] if w == "r"), sum(1 for w in o[first:] if w == "r")
            npass = sum(int(p.split()[0].rstrip("c")) for p in res["passes"])
            verdict(str(rej) == wt.get("steps_rejected") and str(len(o) - rej) == wt.get("steps_done")
                    and len(res["passes"]) == len(o),
                    "the wave: %s steps, %d attempts, %d rejected (%d before the first accepted step, %d after one), "
                    "%d passes; the attempts record agrees with its trailer (steps_done=%s steps_rejected=%s) and its "
                    "schedule (%d lines)" % (wt.get("steps_done"), len(o), rej, early, late, npass, wt.get("steps_done"),
                                             wt.get("steps_rejected"), len(res["passes"])))
            want = {"A": rej == 0, "B": early > 0, "C": late > 0}[case]
            verdict(want, {"A": "it rejects nothing (%d rejected)" % rej,
                           "B": "a rejection before any step is accepted: %d" % early,
                           "C": "a rejection after an accepted step, the case predicted from er and br: %d" % late}[case])
            wrc, werr, whole, wholet = res["whole"]
            nwd = sum(a != b for a, b in zip(wrec, whole)) + abs(len(wrec) - len(whole))
            wsame = all(wholet.get(k) == wt.get(k) for k in COUNTERS)
            verdict(wrc == 0 and nwd == 0 and len(whole) > 0 and wsame,
                    "the whole wave through its attempts record alone, on its own exit test and with no schedule: %d of "
                    "%d record lines differ, and its counters are%s the wave's%s"
                    % (nwd, len(wrec), "" if wsame else " NOT", "" if wrc == 0 else " (rc %d: %s)" % (wrc, werr.strip()[-200:])))
            own_bits = own_passes = 0
            for j, p in sorted(res["parts"].items()):
                if "refused" in p:
                    verdict(False, "body %d: extract_particle refused: %s" % (j, p["refused"]))
                    continue
                x = p["x"]
                if p["rc"]:
                    verdict(False, "body %d: the replay did not run (rc %d): %s" % (j, p["rc"], p["err"].strip()[-300:]))
                    continue
                n, values, bad = compare(wrec, p["rec"], x["bodies"])
                for b in bad[:3]:
                    print("       " + b)
                swept, dist = travel(wrec, j)
                same_c = all(p["tail"].get(k) == wt.get(k) for k in COUNTERS)
                verdict(not bad and n >= 2 and n == len(wrec),
                        "body %d replayed alone: %d samples, the particle and all %d massive bodies, %d values, and every "
                        "sample's step, t, dt_next, dt_last, E and exact-time pair are the wave's%s"
                        % (j, n, x["n_active"], values, "" if not bad else " - NOT: %d mismatches" % len(bad)))
                verdict(same_c, "  and its counters are the wave's: %s" % " ".join(
                    "%s=%s%s" % (k, p["tail"].get(k), "" if p["tail"].get(k) == wt.get(k) else "(wave %s)" % wt.get(k))
                    for k in COUNTERS))
                verdict(swept > 0.05 and dist > 1000.0,
                        "  it moved: %.3f rad about Saturn, summed sample to sample, %.0f km of chords, over the run's %.0f s"
                        % (swept, dist, float.fromhex(wrec[-1][3])))
                if form != "loop":
                    continue
                d_own = particle_differs(wrec, p["own"], j, x["n_active"]) if p["own_rc"] == 0 else 0
                d_pass = sum(1 for a, b in zip(res["passes"], p["ownx_passes"]) if a != b) if p["ownx_rc"] == 0 else 0
                own_bits += d_own > 0
                own_passes += d_pass > 0
                print("       alone on its own steps its bits differ at %d of %d samples; alone on the recorded steps its "
                      "own exit test takes other passes in %d of %d attempts (mean_pc_iterations %s, the wave's %s)"
                      % (d_own, len(wrec), d_pass, len(o), p["ownx_tail"].get("mean_pc_iterations"),
                         wt.get("mean_pc_iterations")))
            if form != "loop":
                continue
            verdict(own_bits > 0 and own_passes > 0,
                    "the records are load bearing here: %d sampled particles take other bits on their own steps, and %d "
                    "other passes on their own exit test" % (own_bits, own_passes))
            ctl = res.get("ctl", {})
            if case == "A":
                if "moved" not in ctl:
                    verdict(False, "the moved-dt control did not run")
                else:
                    k, was, now = ctl["moved"]
                    rc, out, err = ctl["moved_refused"]
                    verdict(rc == 3 and "was recorded from a different run" in err,
                            "attempt %d's dt moved by one ulp (%s to %s): the schedule keyed on the true record is refused "
                            "by name (rc %d)" % (k, was, now, rc))
                    rc, out, err = ctl["moved_rekeyed"]
                    rec, _ = record(out)
                    j0 = sampled(case)[0]
                    nd = particle_differs(wrec, rec, j0, n_active) if rc == 0 else 0
                    verdict(rc == 0 and nd > 0, "  and re-keyed for it, the replay runs and body %d's bits differ at %d "
                            "of %d samples%s" % (j0, nd, len(wrec), "" if rc == 0 else " (rc %d: %s)" % (rc, err.strip()[-200:])))
            else:
                if case == "C":
                    # the precondition, by name: a sample line printed just
                    # before the deleted attempt carries that attempt's dt as
                    # its dt_next, and the deletion moves it with nothing wrong
                    lates, kpre, sample, blocked = ctl.get("deleted_pre", ([], None, 0, []))
                    if kpre is None:
                        verdict(False, "the deletion case's precondition fails, by name: each rejection after an accepted "
                                "step (attempts %s) has a sample line just before it at --sample %d, which carries the "
                                "deleted attempt's own dt as its dt_next; the case needs another --sample" % (lates, sample))
                        continue
                    why = ("it follows rejected attempt %d, and no sample is printed after a rejection" % (kpre - 1)
                           if o[kpre - 1] == "r" else
                           "it follows accepted attempt %d, which brings the accepted count to %d, not a multiple of "
                           "--sample %d" % (kpre - 1, o[:kpre].count("a"), sample))
                    verdict(not sample_before(o, kpre, sample),
                            "the deletion case's precondition, by name: attempt %d has no sample line just before it (%s), "
                            "so no sample carries its dt as dt_next%s"
                            % (kpre, why, "; attempts %s would have, and are not chosen" % blocked if blocked else ""))
                if "deleted" not in ctl:
                    verdict(False, "the deletion control did not run")
                else:
                    k, (rc, out, err) = ctl["deleted"]
                    rec, tail = record(out)
                    j0 = sampled(case)[0]
                    if rc:
                        verdict(False, "attempt %d deleted from both records and re-keyed: the replay did not run (rc %d): %s"
                                % (k, rc, err.strip()[-200:]))
                    elif case == "B":
                        bodies = res["parts"][j0]["x"]["bodies"]
                        n, values, bad = compare(wrec, rec, bodies)
                        nb = body_mismatches(wrec, rec, bodies)
                        verdict(nb > 0, "the rejected first attempt deleted from both records and re-keyed: the replay runs, "
                                "and its bodies differ in %d of its %d sample-and-body pairs: the state, not only the "
                                "dt_next of the sample printed before the attempt (%d comparisons in all; steps_rejected "
                                "%s, the wave's %s)" % (nb, n * len(bodies), len(bad), tail.get("steps_rejected"),
                                                        wt.get("steps_rejected")))
                    else:
                        n, values, bad = compare(wrec, rec, res["parts"][j0]["x"]["bodies"])
                        verdict(not bad and tail.get("steps_rejected") == str(int(wt["steps_rejected"]) - 1),
                                "attempt %d, a rejection after an accepted step, deleted from both records and re-keyed: "
                                "every bit is the wave's (%d mismatches) and only the counters see it, steps_rejected %s "
                                "against %s" % (k, len(bad), tail.get("steps_rejected"), wt.get("steps_rejected")))
                if case == "C" and "deleted_b" in ctl:
                    # the precondition's other side, held to the program: with a
                    # sample line just before it, a deleted late rejection moves
                    # that line's dt_next, from its own dt to the next attempt's,
                    # and nothing else
                    kb, (rc, out, err), adts = ctl["deleted_b"]
                    rec, tail = record(out)
                    bodies = res["parts"][sampled(case)[0]]["x"]["bodies"]
                    m = o[:kb].count("a")
                    heads = [(w[2], t) for w, s in zip(wrec, rec) for t in list(range(7)) + [-3, -2, -1] if w[t] != s[t]]
                    at = next((i for i, w in enumerate(wrec) if w[2] == str(m)), None)
                    verdict(rc == 0 and len(rec) == len(wrec) and body_mismatches(wrec, rec, bodies) == 0
                            and heads == [(str(m), 4)] and at is not None
                            and wrec[at][4] == adts[kb] and rec[at][4] == adts[kb + 1],
                            "attempt %d, a rejection after an accepted step WITH a sample line just before it (the step-%d "
                            "line), deleted the same way: only that line's dt_next moves, from the deleted attempt's dt "
                            "%s to the next attempt's %s; every body and every other field is the wave's (%d fields "
                            "differ: %s)" % (kb, m, adts[kb], adts[kb + 1] if kb + 1 < len(adts) else "?", len(heads),
                                             heads[:4] if rc == 0 else "rc %d: %s" % (rc, err.strip()[-150:])))
                elif case == "C":
                    print("       (no rejection after an accepted step has a sample line just before it at --sample %d here, "
                          "so the deletion of one is not exercised at this format)" % ctl.get("deleted_pre", (0, 0, 0))[2])

    # A wave that is a system of an ensemble, as P2's are: recorded with
    # --member, extracted with --member, and refused for the other member.
    fmt0 = formats[0]
    print("check_replay at %s: a wave that is system 1 of a two-wave ensemble (--member 1)" % fmt0)
    member_case(cft, work, run_args(fmt0, "B"), n_active)

    # 7. A dt file through a pipe, which cannot be rewound, in both forms.
    print("check_replay at %s: a dt file through a pipe, which cannot be rewound" % fmt0)
    pipe_case(cft, work, results[(fmt0, "A", "loop")], waves["A"])

    # 5. Refusals by name, at the first format, on wave A's records.
    res = results[(fmt0, "A", "loop")]
    print("check_replay at %s: refusals" % fmt0)
    j0 = sampled("A")[0]
    if res["wave_rc"] or "x" not in res["parts"].get(j0, {}):
        verdict(False, "wave A did not run or its first particle was not extracted, so nothing can be refused")
    else:
        refusals(cft, work, waves["A"], res, results[(fmt0, "B", "loop")], j0, n_active)

    if ok_all and not args.keep:
        shutil.rmtree(work, ignore_errors=True)
    elif not ok_all:
        print("the waves and records of this run are kept in %s" % work)
    print("RESULT: %s" % ("PASS" if ok_all else "FAIL"))
    return 0 if ok_all else 1


def member_case(cft, work, ra, n_active):
    d = os.path.join(work, "member")
    os.makedirs(d, exist_ok=True)
    w0, w1 = os.path.join(d, "w0.txt"), os.path.join(d, "w1.txt")
    write_wave(w0, False)
    n = write_wave(w1, False, edit=lambda b: b[:n_active + 2] + [shift(b[n_active + 2])] + b[n_active + 3:])
    sys_lines = []
    for k, w in enumerate((w0, w1)):
        with open(w, encoding="utf-8") as f:
            sys_lines += ["system w%d" % k] + [l for l in f.read().splitlines() if l.startswith("body ")]
    ens = os.path.join(d, "ensemble.txt")
    with open(ens, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(["# two ring waves, the second with body %d moved 1 km" % (n_active + 2), "name ringens",
                           "G " + (1.0).hex(), "active %d" % n_active, "N %d" % n, "E 2"] + sys_lines) + "\n")
    att, sch = os.path.join(d, "m1.att"), os.path.join(d, "m1.sched")
    rc, out, err = run([cft, "--problem", ens, "--member", "1"] + ra + ["--attempts-out", att, "--record-iters", sch])
    if rc:
        verdict(False, "member 1 did not run (rc %d): %s" % (rc, err.strip()[-300:]))
        return
    wrec, wt = record(out)
    for j in (n_active, n_active + 2):
        try:
            r = xp.extract(cft, ens, att, sch, j, d, ra, member=1)
        except xp.Refused as e:
            verdict(False, "body %d of member 1: extract_particle refused: %s" % (j, e))
            continue
        rc, out, err = run([cft] + r["args"])
        rec, tail = record(out)
        nn, values, bad = compare(wrec, rec, r["bodies"])
        verdict(rc == 0 and not bad and nn == len(wrec) and all(tail.get(k) == wt.get(k) for k in COUNTERS),
                "body %d of member 1%s replayed alone: %d samples, %d values and the counters are member 1's "
                "(steps_rejected=%s)" % (j, " (the body moved)" if j == n_active + 2 else "", nn, values,
                                         tail.get("steps_rejected")))
    try:
        xp.extract(cft, ens, att, sch, n_active, d, ra, member=0)
        why = None
    except xp.Refused as e:
        why = str(e)
    verdict(why is not None and "was recorded from a different run than this wave" in why,
            "member 1's records extracted for member 0: refused by name")


def pipe_case(cft, work, res, wave):
    """Wave A's accepted-steps dt file and its attempts record, each read as a
    file and through a pipe. The accepted-steps file is padded to lines of
    32 or 128 bytes, both of which divide the C runtime's 4,096-byte buffer,
    and to more than three buffers: a reader that lost a buffer would then
    lose whole lines and run on silently, which is the case that hides."""
    if res["wave_rc"] or not os.path.exists(res.get("dtf", "")):
        verdict(False, "wave A or its --dt-out did not run, so no dt file can go through a pipe")
        return
    d = os.path.join(work, "pipe")
    os.makedirs(d, exist_ok=True)
    with open(res["dtf"], encoding="utf-8") as f:
        lines = [l.strip() for l in f.read().splitlines() if l.strip()]
    width = 32 if max(len(l) for l in lines) < 31 else 128
    body = lines + lines[-1:] * max(0, 3 * 4096 // width + 1 - len(lines))
    data = "".join(l.ljust(width - 1) + "\n" for l in body).encode()
    padded = os.path.join(d, "wave-padded.dt")
    with open(padded, "wb") as f:
        f.write(data)
    head = [cft] + xp.replay_args(res["args"], wave, "")[:-1]      # ... --dt-file, and the file or the pipe next
    wrec, wt = res["wrec"], res["wtail"]
    for label, path, payload, against in (
            ("the accepted-steps file (--dt-out), padded to %d lines of %d bytes, %d bytes in all"
             % (len(body), width, len(data)), padded, data, "the wave's"),
            ("the attempts record, %d bytes" % os.path.getsize(res["att"]), res["att"], None, "the wave's")):
        if payload is None:
            with open(path, "rb") as f:
                payload = f.read()
        rc1, out1, err1 = run(head + [path])
        rc2, out2, err2 = run_piped(head, payload)
        (f1, t1), (p2, t2) = record(out1), record(out2)
        nd_file = sum(a != b for a, b in zip(wrec, f1)) + abs(len(wrec) - len(f1))
        nd_pipe = sum(a != b for a, b in zip(f1, p2)) + abs(len(f1) - len(p2))
        same = all(t1.get(k) == t2.get(k) for k in COUNTERS)
        verdict(rc1 == 0 and rc2 == 0 and nd_file == 0 and nd_pipe == 0 and len(f1) > 0 and same,
                "%s: as a file its replay gives %s %d record lines (%d differ), and through a pipe rc %d, %d of them "
                "differ and its counters are%s the file's%s"
                % (label, against, len(wrec), nd_file, rc2, nd_pipe, "" if same else " NOT",
                   "" if rc2 == 0 else ": " + err2.strip()[-160:]))


def shift(line):
    """A body line with its x moved 1 km."""
    t = line.split()
    t[3] = (float.fromhex(t[3]) + 1.0).hex()
    return " ".join(t)


def refusals(cft, work, wave, res, resb, j, n_active):
    d = os.path.join(work, "refuse")
    os.makedirs(d, exist_ok=True)
    x, att, sch, ra = res["parts"][j]["x"], res["att"], res["sch"], res["args"]
    sub = x["problem"]

    def tool(w, a=att, s=sch, particle=j, sub_file=None):
        try:
            xp.extract(cft, w, a, s, particle, d, ra, None, sub_file)
            return None
        except xp.Refused as e:
            return str(e)

    def program(problem=sub, a=att):
        return run([cft] + [problem if t == sub else a if t == att else t for t in x["args"]])

    def nudge(line, field):
        t = line.split()
        t[field] = ulp_up(t[field], 13) if not t[field].startswith("-") else "-" + ulp_up(t[field][1:], 13)
        return " ".join(t)

    def said(rc, err, phrase):
        return rc == 3 and phrase in err

    # the records replayed for another wave: the sampled particle moved, then another one
    for label, b in (("its particle", j), ("another particle, body %d" % (j + 2), j + 2)):
        w2 = os.path.join(d, "other-%d.txt" % b)
        write_wave(w2, False, edit=lambda bodies, b=b: bodies[:b] + [nudge(bodies[b], 3)] + bodies[b + 1:])
        why = tool(w2)
        verdict(why is not None and "was recorded from a different run than this wave" in why,
                "the records replayed for another wave, %s moved by one ulp: extract_particle refuses by name"
                % label)
        if b == j:
            s2 = edit_lines(sub, os.path.join(d, "other-sub.txt"),
                            lambda L: L[:-1] + [next(l for l in open(w2, encoding="utf-8").read().splitlines()
                                                     if l.startswith("body p%d " % (j - n_active)))])
            rc, out, err = program(s2)
            verdict(said(rc, err, "was recorded from a different run"),
                    "  and the program refuses the re-keyed schedule for that wave's particle by name (rc %d)" % rc)
        else:
            print("       (the sub-system of body %d is the same file for both waves, so only the tool can tell them apart)" % j)

    # the sub-system with a massive body changed, two swapped, its active line changed
    with open(sub, encoding="utf-8") as f:
        sl = f.read().splitlines()
    bi = [i for i, l in enumerate(sl) if l.startswith("body ")]
    ai = next(i for i, l in enumerate(sl) if l.startswith("active "))
    for label, fn, phrase in (
            ("a massive body moved (moon 3's x by one ulp)", lambda L: L[:bi[3]] + [nudge(L[bi[3]], 3)] + L[bi[3] + 1:],
             "is not the wave's massive body 3"),
            ("two massive bodies swapped (1 and 2)", lambda L: L[:bi[1]] + [L[bi[2]], L[bi[1]]] + L[bi[2] + 1:],
             "is not the wave's massive body 1"),
            ("its active line made 8, the particle massive", lambda L: L[:ai] + ["active %d" % (n_active + 1)] + L[ai + 1:],
             "active line is not the wave's")):
        s2 = edit_lines(sub, os.path.join(d, "sub-%d.txt" % len(label)), fn)
        why = tool(wave, sub_file=s2)
        verdict(why is not None and phrase in why, "the sub-system with %s: extract_particle refuses by name" % label)
        rc, out, err = program(s2)
        verdict(said(rc, err, "was recorded from a different run"),
                "  and the program refuses its re-keyed schedule by name (rc %d)" % rc)
    why = tool(wave, particle=3)
    verdict(why is not None and "is not a test particle" in why, "a \"particle\" that is massive body 3: refused by name")
    w3 = os.path.join(d, "no-active.txt")
    write_wave(w3, False, active=False)
    why = tool(w3)
    verdict(why is not None and "has no active line" in why, "a wave with no active line: refused by name")
    # a zonal line on a test particle: the tool refuses it before the program
    # is asked anything (the program refuses every zonal line until P1)
    w4 = edit_lines(wave, os.path.join(d, "zonal-on-a-particle.txt"),
                    lambda L: L[:3] + ["zonal %d 0x1p+0 0x0p+0 0x0p+0" % n_active] + L[3:])
    why = tool(w4)
    verdict(why is not None and "is not one of its %d massive bodies" % n_active in why,
            "a wave whose zonal line names a test particle, body %d: refused by name" % n_active)

    # the attempts record with two outcomes swapped, on wave B's: its rejected
    # first attempt made accepted and its last attempt rejected, so the dts
    # and the count of accepted attempts are what they were and only the
    # outcomes differ
    if resb["wave_rc"] or "x" not in resb["parts"].get(j, {}):
        verdict(False, "wave B did not run, so an outcome cannot be flipped in its record")
    else:
        xb, attb = resb["parts"][j]["x"], resb["att"]
        ob = resb["outc"]

        def swap(L):
            i0, i1 = body_lines(L)[0], body_lines(L)[len(ob) - 1]
            L[i0], L[i1] = L[i0].split()[0] + " a", L[i1].split()[0] + " r"
            return L
        flipped = edit_lines(attb, os.path.join(d, "flipped.att"), swap)
        rc, out, err = run([cft] + [flipped if t == attb else t for t in xb["args"]])
        verdict(ob[0] == "r" and ob[-1] == "a" and said(rc, err, "was recorded from a different run"),
                "wave B's attempts record with two outcomes swapped (attempt 0 to a, attempt %d to r), its dts and its "
                "accepted count unchanged: the program refuses its schedule by name, the outcomes being in the key "
                "(rc %d)" % (len(ob) - 1, rc))
    base = [t for t in x["args"]]
    i = base.index("--replay-iters")
    base = base[:i] + base[i + 2:]
    ens = os.path.join(d, "ensemble.txt")
    with open(sub, encoding="utf-8") as f:
        sl = [l for l in f.read().splitlines() if not l.startswith("#")]
    head = [l for l in sl if not l.startswith("body ")]
    bodies = [l for l in sl if l.startswith("body ")]
    with open(ens, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(head + ["N %d" % len(bodies), "E 2", "system a"] + bodies + ["system b"] + bodies) + "\n")
    ens_args = [ens if t == sub else t for t in base]
    steps = int(base[base.index("--steps") + 1])
    cases = [
        ("the attempts form on an ensemble", ens_args, "is an attempts record, which is one system's attempts"),
        ("--attempts-out on an ensemble", ["--problem", ens] + list(ra) + ["--attempts-out", os.path.join(d, "x.att")],
         "--attempts-out records one system's attempts"),
        ("a malformed attempt line", subst(base, att, edit_lines(att, os.path.join(d, "bad.att"),
                                           lambda L: L[:body_lines(L)[1]] + [L[body_lines(L)[1]].split()[0] + " x"]
                                           + L[body_lines(L)[1] + 1:])), "is not an attempt"),
        ("a record cut short (no next line)", subst(base, att, edit_lines(att, os.path.join(d, "cut.att"),
                                                                       lambda L: L[:-1])), "has no next line"),
        ("a line after the record's next line", subst(base, att, edit_lines(att, os.path.join(d, "long.att"),
                                                                           lambda L: L + [L[body_lines(L)[0]]])),
         "follows the attempts record's next line"),
        ("too few accepted attempts for the run", subst(base, str(steps), str(steps + 1), after="--steps"),
         "accepted attempts, the run needs"),
        ("an attempts record without --epsilon 0", subst(base, "0", "1e-9", after="--epsilon"), "it needs --epsilon 0"),
    ]
    for label, a, phrase in cases:
        rc, out, err = run([cft] + a)
        verdict(said(rc, err, phrase) and not [l for l in out.splitlines() if l.startswith("# steps_done=")],
                "%s is refused by name, with no trailer (rc %d)" % (label, rc))
    # --run-key is the run line --record-iters writes, for the replay's own arguments
    rk = os.path.join(d, "own.sched")
    rc1, key, _ = run([cft] + base + ["--run-key"])
    rc2, out, _ = run([cft] + base + ["--record-iters", rk])
    line = next((l[6:] for l in open(rk, encoding="utf-8").read().splitlines() if l.startswith("# run ")), "") if rc2 == 0 else ""
    verdict(rc1 == 0 and key.strip() == line and kv(line).get("dt_file") == "attempts"
            and kv(res_key(cft, ra, wave)).get("dt_file") == "no",
            "--run-key prints exactly the run line --record-iters writes for the replay (rc %d, %d); its dt_file is %s, "
            "and the wave's %s" % (rc1, rc2, kv(line).get("dt_file"), kv(res_key(cft, ra, wave)).get("dt_file")))


def res_key(cft, ra, wave):
    rc, out, _ = run([cft, "--problem", wave] + list(ra) + ["--run-key"])
    return out.strip() if rc == 0 else ""


def subst(args, old, new, after=None):
    """args with the value `old` replaced by `new`; with `after`, only the
    value that follows that option."""
    out = list(args)
    for i, a in enumerate(out):
        if a == old and (after is None or (i > 0 and out[i - 1] == after)):
            out[i] = new
            if after is not None:
                break
    return out


if __name__ == "__main__":
    sys.exit(main())
