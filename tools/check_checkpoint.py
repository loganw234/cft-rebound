#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The checkpoint gate, across processes.
#
# tests/gate_real.c already compares a restart against an uninterrupted
# run inside one process, and that is a real test: the engine is reset
# when a new simulation binds to it, so a state that was not carried is
# caught. What it cannot test is the thing a checkpoint is FOR. A long
# run exits. The machine may reboot. Another process picks the archive
# up later, with a freshly opened engine, freshly derived constants and
# nothing whatever inherited from the run that wrote the file.
#
# So this drives three separate processes per format:
#
#   gate_real --fpNNN --straight            -> the reference dump
#   gate_real --fpNNN --save   ckpt.bin     -> stops after --steps-a
#   gate_real --fpNNN --resume ckpt.bin     -> the continued dump
#
# and requires the two dumps to be identical. The dumps are exact hex
# floats (%a), so the comparison is over bits rather than over a decimal
# rendering of them: a last-bit difference shows as a different line.
#
# FOUR SEQUENCES, and THE SWEEP.
#   plain     a checkpoint, 60 steps then 120, at a fixed step.
#   old       plain's checkpoint written with the field list as it was
#             before 2026-11-02 (gate_real --old-list), which is how the
#             suite holds "old archives still load" with no binary in the
#             tree. Its resume must equal plain's archive resumed with the
#             clock's and the steps' loads removed, and differ from the
#             straight run in the clock alone (old_sequence()).
#   regrow    a particle REMOVAL before the checkpoint and one added back
#             below the original high-water mark afterwards. REBOUND keeps
#             the arrays a count change re-reads at N_allocated rather
#             than at the live 3N, so a shrink strands a tail and a
#             regrow under the mark reads it straight back; the engine
#             here is one process-wide instance whose shadow outlives the
#             simulation that filled it, so an in-process load can
#             inherit very nearly the right tail and be right by accident.
#             Three processes inherit nothing.
#   adaptive  plain at epsilon 1e-9, 20 steps then 20. The step control
#             computes dt in the run's format, so above binary64 the
#             engine's dt and dt_last carry bits r->dt and r->dt_last_done
#             cannot; until 2026-11-02 a load rebuilt both from those
#             binary64 views, and this sequence differed at binary128 in
#             a particle coordinate. Nothing had checkpointed an adaptive
#             run before it.
#   splits    THE 2-OF-40 CASE, whose setup no record kept, so it is
#             written here: gate_real's four bodies at a fixed step of
#             0.01, a checkpoint after s = 1..40 steps, then 40 more
#             steps after each, against one straight run dumped after
#             each of its steps 41..80. ROADMAP.md recorded the clock as
#             "1 ulp in t on 2 of 40 step splits", at binary128 and at
#             binary256. On af2f400 this sweep measured exactly that, at
#             splits 29 and 30, at both formats - and 34 of 40 splits in
#             t's WIDE bytes. One run saves after every step
#             (--save-each), one straight run dumps after every step
#             (--dump-from), and each resume is its own process.
#
# THE DUMP'S LAST THREE LINES are the engine's clock (t, t_hi, t_lo), its
# two steps (dt, dt_last), each as wide bytes, and a digest of the WIDE
# state - all there because a control was not sensitive without them.
# The particles are the binary64 view, and above binary64 the view does
# not see everything the state carries: a corrector that starts from
# different coefficients re-converges to within 1e-34 at binary128 and
# rounds to the same doubles, so NO_ALIAS matched the straight run at
# fp128 and fp256 until gate_real hashed the bytes. And r->t is t's
# binary64 view, so a resume that lost t's low bits matched on the `t`
# line in 38 of the sweep's 40 splits. Comparing a view is vacuous
# through the observable; these compare the bytes.
#
# THE NEGATIVE CONTROLS ARE NOT OPTIONAL. Each is the resume again with
# one part of the restore dropped, and each must make its sequence DIFFER:
#
#   CFT_REBOUND_NO_ADOPT   the shim discards the whole loaded state
#   CFT_REBOUND_NO_ALIAS   it adopts the state but drops the high-water
#                          mark and the arrays sized by it - exactly the
#                          world before 2026-09-11
#   CFT_REBOUND_NO_CLOCK   it adopts everything but the wide clock, which
#                          it rebuilds from r->t - the world before
#                          2026-11-02. The sweep's control: with it, the
#                          2-of-40 case must fail again
#   CFT_REBOUND_NO_STEP    it adopts everything but the wide dt and
#                          dt_last, which it rebuilds from r->dt and
#                          r->dt_last_done
#
# A control is required only where it can bite, and the list says where.
# NO_STEP is inert by construction at binary64, where dt and dt_last ARE
# binary64 values, and at a fixed step, where they are 0.01 at every
# format. NO_CLOCK is inert in the regrow sequence, because the add after
# the load rebuilds the clock from r->t in the straight run and the
# resumed one alike.
#
#   python tools/check_checkpoint.py [--build build] [--quick]
#
# --quick is binary64 only, and skips the sweep by name: it is the full
# tier's (make check). An artifact in CFT_REBOUND_ARTIFACT is inherited
# by every process, so this gate runs on a card without a flag of its
# own.

import argparse
import os
import subprocess
import sys
import tempfile

FORMATS = ["fp64", "fp128", "fp256"]
WIDE = ["fp128", "fp256"]

# name, extra gate_real flags, and the controls that must make it differ,
# each with the formats at which it can.
SEQUENCES = [
    ("plain", [],
     [("CFT_REBOUND_NO_ADOPT", FORMATS), ("CFT_REBOUND_NO_CLOCK", FORMATS)]),
    ("regrow", ["--regrow"],
     [("CFT_REBOUND_NO_ADOPT", FORMATS), ("CFT_REBOUND_NO_ALIAS", FORMATS)]),
    ("adaptive", ["--epsilon", "1e-9", "--steps-a", "20", "--steps-b", "20"],
     [("CFT_REBOUND_NO_ADOPT", FORMATS), ("CFT_REBOUND_NO_CLOCK", FORMATS),
      ("CFT_REBOUND_NO_STEP", WIDE)]),
]

# The sweep: a checkpoint after each of s = 1..SPLITS steps, and AFTER
# more steps after each one.
SPLITS = 40
AFTER = 40
SWEEP_CONTROL = "CFT_REBOUND_NO_CLOCK"


def run(exe, args, want_stdout, extra_env=None):
    """Returns (ok, stdout, stderr). A phase writes its dump to stdout
    and everything else to stderr, so the two never mix."""
    env = None
    if extra_env:
        env = dict(os.environ)
        env.update(extra_env)
    try:
        p = subprocess.run([exe] + args, capture_output=True, text=True, env=env)
    except OSError as e:
        return False, "", str(e)
    ok = (p.returncode == 0) and (not want_stdout or p.stdout.strip() != "")
    return ok, p.stdout, p.stderr


def first_difference(a, b):
    la, lb = a.splitlines(), b.splitlines()
    for i in range(min(len(la), len(lb))):
        if la[i] != lb[i]:
            return "line %d: straight %r, resumed %r" % (i + 1, la[i], lb[i])
    if len(la) != len(lb):
        return "%d lines straight, %d resumed" % (len(la), len(lb))
    return "no difference found, which contradicts the comparison"


def one_sequence(exe, fmt, workdir, seq, extra, controls, keep=None):
    ckpt = os.path.join(workdir, "ckpt_%s_%s.bin" % (fmt, seq))
    flag = "--" + fmt

    ok, straight, err = run(exe, [flag] + extra + ["--straight"], True)
    if not ok:
        return False, "the straight run failed: %s" % err.strip().splitlines()[-1:]
    if keep is not None:
        keep["straight"], keep["ckpt"] = straight, ckpt

    ok, _, err = run(exe, [flag] + extra + ["--save", ckpt], False)
    if not ok:
        return False, "the save failed: %s" % err.strip()
    if not os.path.exists(ckpt):
        return False, "no archive at %s" % ckpt
    size = os.path.getsize(ckpt)

    ok, resumed, err2 = run(exe, [flag] + extra + ["--resume", ckpt], True)
    if not ok:
        return False, "the resume failed: %s" % err2.strip()

    # What the two side processes said about themselves, for the record.
    note = " | ".join(s.strip() for s in (err.strip(), err2.strip()) if s.strip())

    if straight != resumed:
        return False, "the continued run DIFFERS - %s" % first_difference(straight, resumed)

    # The physics counter for the adaptive sequence: its step control must
    # have moved dt off the 0.01 every sequence starts from, or it was a
    # fixed-step run passing for the wrong reason.
    if "--epsilon" in extra:
        dts = [l for l in straight.splitlines() if l.startswith("dt ")]
        if not dts or float.fromhex(dts[0].split()[1]) == 0.01:
            return False, "the step control never moved dt (%s): not an adaptive run" % dts

    # The controls. Same three processes, except the resume is told to
    # drop part of what it loaded. If one of THOSE still matches, the
    # comparison above was not testing what it claims to test.
    used = []
    for var, where in controls:
        if fmt not in where:
            continue
        ok, control, _ = run(exe, [flag] + extra + ["--resume", ckpt], True, {var: "1"})
        if not ok:
            return False, "the %s control run would not start" % var
        if control == straight:
            return False, ("NOT SENSITIVE to %s - a resume that drops that part of "
                           "the loaded state still matches, so this comparison "
                           "proves nothing" % var)
        used.append(var)

    return True, "%d values identical, %s differ, archive %d B%s" % (
        len(straight.splitlines()), " and ".join(used), size,
        (" | " + note) if note else "")


def old_sequence(exe, fmt, workdir, straight, new_ckpt):
    """AN ARCHIVE WRITTEN BEFORE 2026-11-02 STILL LOADS, and does what
    every load did before then: it restores the wide state exactly and
    rebuilds the clock and steps from binary64. gate_real --old-list
    writes one, with the field list as it was (field for field af2f400's
    own; tools/check_archive.py's docstring). Its resume must equal the
    plain sequence's archive resumed with the clock's AND the steps'
    loads removed, and must differ from the straight run only in the
    clock: at this fixed step the trajectory does not depend on it."""
    flag = "--" + fmt
    ckpt = os.path.join(workdir, "ckpt_%s_old.bin" % fmt)
    ok, _, err = run(exe, [flag, "--old-list", "--save", ckpt], False)
    if not ok:
        return False, "the save with the old field list failed: %s" % err.strip()
    ok, old, err2 = run(exe, [flag, "--resume", ckpt], True)
    if not ok:
        return False, "the resume of the old archive failed: %s" % err2.strip()
    ok, rebuilt, _ = run(exe, [flag, "--resume", new_ckpt], True,
                         {"CFT_REBOUND_NO_CLOCK": "1", "CFT_REBOUND_NO_STEP": "1"})
    if not ok:
        return False, "the resume with the clock's and steps' loads removed would not start"
    if old != rebuilt:
        return False, ("an old archive does not resume as a new one whose clock and steps "
                       "were rebuilt from binary64 - %s" % first_difference(rebuilt, old))
    kinds = sorted(set(kind(a) for a, b in zip(straight.splitlines(), old.splitlines()) if a != b))
    if not kinds:
        return False, ("the old archive resumed identically to the straight run, clock "
                       "included, so this case cannot show what an old archive concedes")
    if "other" in kinds:
        return False, ("the old archive's resume differs from the straight run beyond its "
                       "clock and steps - %s" % first_difference(straight, old))
    return True, ("resumes as one whose clock and steps were rebuilt from binary64, and "
                  "differs from the straight run in %s lines only | %s" % (
                      " and ".join(kinds), err2.strip()))


def kind(line):
    """Which observable a dump line is: r->t's binary64 view, the
    engine's wide clock or steps, or anything else."""
    if line.startswith("t "):
        return "view"
    if line.startswith(("clock ", "steps ")):
        return "time"
    return "other"


def parse_at(text):
    """--dump-from output: "at <n>" then a dump, repeated. Returns
    {n: [lines]}."""
    out, cur = {}, None
    for line in text.splitlines():
        if line.startswith("at "):
            cur = int(line.split()[1])
            out[cur] = []
        elif cur is not None:
            out[cur].append(line)
    return out


def resume_all(exe, flag, prefix, straight, extra_env=None):
    """Resumes every split's archive for AFTER steps, in its own process.
    Returns (error or None, {kind: [splits]}, [(split, line pairs)])."""
    by = {"view": [], "time": [], "other": []}
    shown = []
    for s in range(1, SPLITS + 1):
        ok, res, err = run(exe, [flag, "--steps-b", str(AFTER),
                                 "--resume", "%s%d.bin" % (prefix, s)], True, extra_env)
        if not ok:
            return "the resume of split %d failed: %s" % (s, err.strip()), by, shown
        want = straight.get(s + AFTER)
        if want is None:
            return "the straight run has no dump at step %d" % (s + AFTER), by, shown
        got = res.splitlines()
        if len(got) != len(want):
            by["other"].append(s)
            continue
        pairs = [(a, b) for a, b in zip(want, got) if a != b]
        for k in sorted(set(kind(a) for a, _ in pairs)):
            by[k].append(s)
        if pairs and len(shown) < 3:
            shown.append((s, pairs))
    return None, by, shown


def splits_sweep(exe, fmt, workdir):
    flag = "--" + fmt
    prefix = os.path.join(workdir, "split_%s_" % fmt)
    ok, out, err = run(exe, [flag, "--steps-a", str(SPLITS), "--steps-b", str(AFTER),
                             "--dump-from", str(AFTER + 1), "--straight"], True)
    if not ok:
        return False, "the straight run failed: %s" % err.strip().splitlines()[-1:]
    straight = parse_at(out)
    if sorted(straight) != list(range(AFTER + 1, SPLITS + AFTER + 1)):
        return False, "the straight run dumped steps %s, not %d..%d" % (
            sorted(straight)[:3], AFTER + 1, SPLITS + AFTER)
    ok, _, err = run(exe, [flag, "--steps-a", str(SPLITS), "--save-each", prefix], False)
    if not ok:
        return False, "the saves failed: %s" % err.strip()
    missing = [s for s in range(1, SPLITS + 1) if not os.path.exists("%s%d.bin" % (prefix, s))]
    if missing:
        return False, "no archive for splits %s" % missing[:5]

    why, by, shown = resume_all(exe, flag, prefix, straight)
    if why:
        return False, why
    bad = sorted(set(by["view"]) | set(by["time"]) | set(by["other"]))
    counts = "r->t differs on %d %s, the clock or steps on %d, other lines on %d" % (
        len(by["view"]), by["view"], len(by["time"]), len(by["other"]))
    if bad:
        detail = "; ".join("split %d: %s" % (s, " / ".join("%s -> %s" % (a, b) for a, b in p))
                           for s, p in shown)
        return False, "%d of %d resumes DIFFER (%s). First: %s" % (len(bad), SPLITS, counts, detail)

    # The control: the same forty resumes with the clock's load removed.
    # The case must fail again, or the sweep was not testing the clock.
    why, cby, _ = resume_all(exe, flag, prefix, straight, {SWEEP_CONTROL: "1"})
    if why:
        return False, "the %s control: %s" % (SWEEP_CONTROL, why)
    cbad = sorted(set(cby["view"]) | set(cby["time"]) | set(cby["other"]))
    if not cbad:
        return False, ("NOT SENSITIVE to %s - forty resumes with the clock's load removed "
                       "all match, so this sweep proves nothing" % SWEEP_CONTROL)
    return True, ("%d of %d resumes identical (%s); with %s, %d of %d differ: r->t on %d %s, "
                  "the clock or steps on %d" % (
                      SPLITS, SPLITS, counts, SWEEP_CONTROL, len(cbad), SPLITS,
                      len(cby["view"]), cby["view"], len(cby["time"])))


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true")
    a = ap.parse_args()

    exe = os.path.join(a.build, "gate_real")
    if os.name == "nt" and not os.path.exists(exe):
        exe += ".exe"
    if not os.path.exists(exe):
        print("check_checkpoint: no %s; run `make` first" % exe)
        print("RESULT: FAIL")
        return 1

    formats = FORMATS[:1] if a.quick else FORMATS
    art = os.environ.get("CFT_REBOUND_ARTIFACT", "")
    print("check_checkpoint: three separate processes per format and sequence, "
          "on %s" % (art if art else "the software backend"))
    print("                  plain = a checkpoint; old = plain's, written as before "
          "2026-11-02; regrow = a removal, a checkpoint and one added back under the "
          "mark; adaptive = plain at epsilon 1e-9")
    print("                  splits = a checkpoint after each of s = 1..%d steps, then %d "
          "more, at a fixed step of 0.01 (the 2-of-40 case)" % (SPLITS, AFTER))
    print("                  each also runs controls that drop part of the loaded "
          "state and must differ")

    allok = True
    with tempfile.TemporaryDirectory(prefix="cft_ckpt_") as workdir:
        for fmt in formats:
            for seq, extra, controls in SEQUENCES:
                keep = {} if seq == "plain" else None
                ok, detail = one_sequence(exe, fmt, workdir, seq, extra, controls, keep)
                print("  %-6s %-8s %s  %s" % (fmt, seq, "PASS" if ok else "FAIL", detail),
                      flush=True)
                allok = allok and ok
                if keep is not None:
                    if ok:
                        ok, detail = old_sequence(exe, fmt, workdir, keep["straight"], keep["ckpt"])
                    else:
                        detail = "not run: it resumes the plain sequence's archive, which failed"
                    print("  %-6s %-8s %s  %s" % (fmt, "old", "PASS" if ok else "FAIL", detail),
                          flush=True)
                    allok = allok and ok
            if a.quick:
                # The sweep is eighty resumes, about two and a half minutes at
                # fp64 on a busy host, and the case it holds is binary128's
                # and binary256's. --quick checks the clock at fp64 through
                # plain and adaptive and their CFT_REBOUND_NO_CLOCK control.
                print("  %-6s %-8s SKIPPED  --quick: the 2-of-40 sweep is the full tier's "
                      "(make check), at every format" % (fmt, "splits"), flush=True)
                continue
            ok, detail = splits_sweep(exe, fmt, workdir)
            print("  %-6s %-8s %s  %s" % (fmt, "splits", "PASS" if ok else "FAIL", detail),
                  flush=True)
            allok = allok and ok

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
