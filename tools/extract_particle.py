#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# One test particle of a wave, alone with the wave's massive bodies, and a
# corrector schedule this tree accepts for it (ROUND4.md, P3).
#
#   python tools/extract_particle.py --exe build/ias15_cft --wave WAVE.txt [--member K]
#          --attempts WAVE.att --schedule WAVE.sched --particle J --out DIR [--sub FILE]
#          -- <the wave's run arguments: --format, --dt, --epsilon, --steps, --sample, ...>
#
# The wave ran with --attempts-out WAVE.att and --record-iters WAVE.sched:
# the attempts record holds every step attempt's dt and outcome, and the
# schedule every attempt's corrector passes, a line an attempt in each. This
# writes the sub-system, DIR/<name>.txt: the wave's G, active and zonal
# lines and its first N_active body lines, verbatim and in order, then body
# J, the particle. Then it re-keys the wave's schedule for the sub-system's
# replay, DIR/<name>.sched, and prints the replay's arguments:
#
#   ias15_cft --problem DIR/<name>.txt <the run arguments, --epsilon 0>
#             --dt-file WAVE.att --replay-iters DIR/<name>.sched
#
# Why a re-key, and not a key that leaves the test particles out. Every
# wave of a ring starts its massive bodies from one state, so a key over
# the massive bodies and the settings alone could not tell one wave's
# schedule from another's. The wave's own key is written before its first
# attempt, so it cannot carry the attempts either. So the program's rule
# stays as it was - a schedule replays only the run whose key it carries -
# and the step that crosses from the wave to the sub-system is here, where
# it is asserted before anything is written:
#   - the schedule and the attempts record were written by THIS wave's run:
#     each one's `# run` line is the key the program gives for the wave and
#     these arguments (--run-key), and they hold one attempt a line alike;
#   - J is a test particle: the wave has an active line and N_active <= J;
#   - the sub-system's G, active and zonal lines are the wave's, and its
#     first N_active body lines are the wave's first N_active, line for
#     line and in order, and its last is the wave's body J.
# The re-keyed schedule's key is the program's own for the sub-system's
# replay. It hashes the sub-system's bodies, the particle among them, and
# the attempts record's dts and outcomes, so it is refused by name for any
# other particle, any changed massive body, or any other attempts. Every
# key is the program's: nothing here hashes anything, so no copy of the
# key can drift from sched_key() in src/ias15_cft.c.
#
# --sub FILE asserts and re-keys an existing sub-system file instead of
# writing one. tools/check_replay.py uses it to show the assertion refuses.
#
# A refusal exits 3 with a line that starts "extract_particle: refused:".

import argparse
import os
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))

# the run arguments a wave's are, without what records or replays it: this
# tool supplies the replay's own
RECORD_OPTS = ("--record-iters", "--replay-iters", "--attempts-out", "--dt-out", "--run-key", "--member",
               "--problem")
ATT_MAGIC = "# cft-rebound attempts v1"


class Refused(Exception):
    pass


def refuse(text):
    raise Refused(text)


def read_lines(path):
    with open(path, encoding="utf-8") as f:
        return f.read().splitlines()


def wave_lines(path, member):
    """The wave's G, active and zonal lines and the body lines of its system
    (system `member` of an ensemble file), each exactly as written."""
    head = {"name": None, "G": None, "active": None, "zonal": None}
    systems, cur = [], None
    for line in read_lines(path):
        tok = line.split()
        if not tok or tok[0].startswith("#"):
            continue
        if tok[0] in head:
            if head[tok[0]] is not None and tok[0] != "name":
                refuse("%s has a second %s line" % (path, tok[0]))
            head[tok[0]] = line
        elif tok[0] in ("N", "E"):
            pass
        elif tok[0] == "system":
            cur = []
            systems.append(cur)
        elif tok[0] == "body":
            if cur is None:
                cur = []
                systems.append(cur)
            cur.append(line)
        else:
            refuse("%s has a line '%s' this tool does not know, and would copy wrongly" % (path, tok[0]))
    if member is None:
        if len(systems) != 1:
            refuse("%s holds %d systems: name one with --member" % (path, len(systems)))
        bodies = systems[0]
    else:
        if not 0 <= member < len(systems):
            refuse("--member %d, and %s holds %d systems" % (member, path, len(systems)))
        bodies = systems[member]
    if head["G"] is None:
        refuse("%s has no G line" % path)
    return head, bodies


def program_key(exe, args):
    p = subprocess.run([exe] + args + ["--run-key"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if p.returncode != 0 or len(p.stdout.splitlines()) != 1:
        refuse("the program gave no run key for %s (rc %d): %s" % (" ".join(args), p.returncode, p.stderr.strip()[-300:]))
    return p.stdout.strip()


def run_line(lines, path, magic):
    if not lines or lines[0] != magic:
        refuse("%s is not what it is passed as: its first line is not \"%s\"" % (path, magic))
    if len(lines) < 2 or not lines[1].startswith("# run "):
        refuse("%s has no run line" % path)
    return lines[1][6:]


def replay_args(run_args, sub, attempts):
    """The wave's run arguments as the sub-system's replay takes them: its
    problem, epsilon 0, and the attempts record as the dt file in place of
    any the wave ran from."""
    out, i = [], 0
    while i < len(run_args):
        a = run_args[i]
        if a in RECORD_OPTS:
            refuse("the run arguments hold %s, and this tool supplies the replay's own: pass the wave's run without "
                   "its problem, its records and --member" % a)
        if a in ("--epsilon", "--dt-file"):
            i += 2
            continue
        out.append(a)
        i += 1
    return ["--problem", sub] + out + ["--epsilon", "0", "--dt-file", attempts]


def extract(exe, wave, attempts, schedule, particle, out_dir, run_args, member=None, sub=None):
    """Assert, write and re-key; returns what the replay needs, or raises
    Refused with the reason."""
    head, bodies = wave_lines(wave, member)
    if head["active"] is None:
        refuse("%s has no active line, so it has no test particle: a particle replays alone only if it adds "
               "nothing to the massive bodies" % wave)
    n_active = int(head["active"].split()[1])
    if not n_active <= particle < len(bodies):
        refuse("body %d of %s is not a test particle: they are bodies %d to %d (active %d)"
               % (particle, wave, n_active, len(bodies) - 1, n_active))

    # the records are this wave's: their run lines are the program's key for it
    wave_args = ["--problem", wave] + (["--member", str(member)] if member is not None else []) + list(run_args)
    for a in run_args:
        if a in RECORD_OPTS:
            refuse("the run arguments hold %s, and this tool supplies the replay's own: pass the wave's run without "
                   "its problem, its records and --member" % a)
    wave_key = program_key(exe, wave_args)
    s_lines, a_lines = read_lines(schedule), read_lines(attempts)
    s_key = run_line(s_lines, schedule, "# cft-rebound corrector schedule v1")
    a_key = run_line(a_lines, attempts, ATT_MAGIC)
    for what, path, key in (("schedule", schedule, s_key), ("attempts record", attempts, a_key)):
        if key != wave_key:
            refuse("the %s %s was recorded from a different run than this wave:\n  it was recorded from: %s\n"
                   "  this wave is:          %s" % (what, path, key, wave_key))
    n_sched = sum(1 for l in s_lines if l and not l.startswith("#"))
    n_att = sum(1 for l in a_lines if l.strip() and not l.startswith("#") and l.split()[-1] != "next")
    if n_sched != n_att:
        refuse("the schedule holds %d attempts and the attempts record %d: one run writes one line an attempt in each"
               % (n_sched, n_att))

    # the sub-system: written here, or asserted from the file given
    stem = "%s-body%d" % (os.path.splitext(os.path.basename(wave))[0], particle)
    if member is not None:
        stem += "-m%d" % member
    if sub is None:
        os.makedirs(out_dir, exist_ok=True)
        sub = os.path.join(out_dir, stem + ".txt")
        name = (head["name"].split()[1] if head["name"] else "wave")[:40] + "-b%d" % particle
        text = ["# tools/extract_particle.py: body %d of %s%s, a test particle, alone with the wave's %d massive bodies"
                % (particle, os.path.basename(wave), "" if member is None else ", system %d" % member, n_active),
                "name " + name, head["G"], head["active"]]
        if head["zonal"] is not None:
            text.append(head["zonal"])
        text += bodies[:n_active] + [bodies[particle]]
        with open(sub, "w", encoding="utf-8", newline="\n") as f:
            f.write("\n".join(text) + "\n")
    s_head, s_bodies = wave_lines(sub, None)
    for k in ("G", "active", "zonal"):
        if s_head[k] != head[k]:
            refuse("the sub-system's %s line is not the wave's:\n  wave:       %s\n  sub-system: %s" % (k, head[k], s_head[k]))
    want = bodies[:n_active] + [bodies[particle]]
    if len(s_bodies) != len(want):
        refuse("the sub-system has %d bodies, and the wave's %d massive bodies and one particle are %d"
               % (len(s_bodies), n_active, len(want)))
    for i, (w, s) in enumerate(zip(want, s_bodies)):
        if w != s:
            who = "massive body %d" % i if i < n_active else "particle, body %d of the wave" % particle
            refuse("the sub-system's body %d is not the wave's %s, line for line:\n  wave:       %s\n  sub-system: %s"
                   % (i, who, w, s))

    # the re-key: the program's own key for the replay, on the wave's schedule
    args = replay_args(run_args, sub, attempts)
    sub_key = program_key(exe, args)
    sched = os.path.join(os.path.dirname(sub) if out_dir is None else out_dir, stem + ".sched")
    out = [s_lines[0], "# run " + sub_key,
           "# re-keyed by tools/extract_particle.py for body %d of %s, from the schedule %s, whose run was: %s"
           % (particle, os.path.basename(wave), os.path.basename(schedule), wave_key)] + s_lines[2:]
    with open(sched, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    # sub body i is wave body i below N_active, and the particle after them
    return {"problem": sub, "schedule": sched, "args": args + ["--replay-iters", sched],
            "wave_key": wave_key, "sub_key": sub_key, "bodies": list(range(n_active)) + [particle],
            "n_active": n_active}


def main():
    ap = argparse.ArgumentParser(description="one test particle of a wave, alone with its massive bodies")
    ap.add_argument("--exe", required=True)
    ap.add_argument("--wave", required=True)
    ap.add_argument("--member", type=int, default=None)
    ap.add_argument("--attempts", required=True)
    ap.add_argument("--schedule", required=True)
    ap.add_argument("--particle", type=int, required=True, help="the particle's body index in the wave's system")
    ap.add_argument("--out", default=None, help="directory for the sub-system and its schedule")
    ap.add_argument("--sub", default=None, help="assert and re-key this sub-system file instead of writing one")
    ap.add_argument("run", nargs=argparse.REMAINDER, help="-- then the wave's run arguments")
    a = ap.parse_args()
    run_args = a.run[1:] if a.run[:1] == ["--"] else a.run
    if a.out is None and a.sub is None:
        ap.error("--out DIR, or --sub FILE")
    try:
        r = extract(a.exe, a.wave, a.attempts, a.schedule, a.particle, a.out, run_args, a.member, a.sub)
    except Refused as e:
        print("extract_particle: refused: %s" % e, file=sys.stderr)
        return 3
    print("problem %s" % r["problem"])
    print("schedule %s" % r["schedule"])
    print("bodies %s" % " ".join(str(b) for b in r["bodies"]))
    print("replay %s" % " ".join(r["args"]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
