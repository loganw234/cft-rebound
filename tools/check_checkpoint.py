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
# TWO SEQUENCES, and the second is the one that needs three processes
# most. --regrow adds a particle REMOVAL before the checkpoint and adds
# one back below the original high-water mark afterwards. REBOUND keeps
# the arrays a count change re-reads at N_allocated rather than at the
# live 3N, so a shrink strands a tail and a regrow under the mark reads
# it straight back; the engine here is one process-wide instance whose
# shadow outlives the simulation that filled it, so an in-process load
# can inherit very nearly the right tail and be right by accident.
# Three processes inherit nothing.
#
# THE NEGATIVE CONTROLS ARE NOT OPTIONAL, and are a fourth process per
# sequence:
#
#   CFT_REBOUND_NO_ADOPT   the shim discards the whole loaded state
#   CFT_REBOUND_NO_ALIAS   the shim adopts the state but drops the
#                          high-water mark and the arrays sized by it -
#                          which is exactly the world before 2026-09-11
#
# and each required run must DIFFER. Without them the comparison is not
# evidence: at the default step counts only a few of the dumped lines
# differ between binary64 and binary256, so agreement on the rest would
# say nothing about whether the wide state was carried at all.
#
# The dump's last line is a digest of the WIDE state, and it is there
# because the controls were not sensitive without it. The particles are
# the binary64 view, and above binary64 the view does not see everything
# the state carries: a corrector that starts from different coefficients
# re-converges to within 1e-34 at binary128 and rounds to the same
# doubles. NO_ALIAS therefore matched the straight run at fp128 and
# fp256 - a control that could not fail - until gate_real began hashing
# the bytes rather than what they round to.
#
#   python tools/check_checkpoint.py [--build build] [--quick]
#
# --quick is binary64 only. An artifact in CFT_REBOUND_ARTIFACT is
# inherited by all three processes, so this gate runs on a card without
# a flag of its own.

import argparse
import os
import subprocess
import sys
import tempfile

FORMATS = ["fp64", "fp128", "fp256"]

# name, extra gate_real flags, the controls that must make it differ
SEQUENCES = [
    ("plain", [], ["CFT_REBOUND_NO_ADOPT"]),
    ("regrow", ["--regrow"], ["CFT_REBOUND_NO_ADOPT", "CFT_REBOUND_NO_ALIAS"]),
]


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


def one_sequence(exe, fmt, workdir, seq, extra, controls):
    ckpt = os.path.join(workdir, "ckpt_%s_%s.bin" % (fmt, seq))
    flag = "--" + fmt

    ok, straight, err = run(exe, [flag] + extra + ["--straight"], True)
    if not ok:
        return False, "the straight run failed: %s" % err.strip().splitlines()[-1:]

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

    # The controls. Same three processes, except the resume is told to
    # drop part of what it loaded. If one of THOSE still matches, the
    # comparison above was not testing what it claims to test.
    for var in controls:
        ok, control, _ = run(exe, [flag] + extra + ["--resume", ckpt], True, {var: "1"})
        if not ok:
            return False, "the %s control run would not start" % var
        if control == straight:
            return False, ("NOT SENSITIVE to %s - a resume that drops that part of "
                           "the loaded state still matches, so this comparison "
                           "proves nothing" % var)

    return True, "%d values identical, %s differ, archive %d B%s" % (
        len(straight.splitlines()), " and ".join(controls), size,
        (" | " + note) if note else "")


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
    print("                  plain = a checkpoint; regrow = a removal, a checkpoint "
          "and one added back under the mark")
    print("                  each also runs controls that drop part of the loaded "
          "state and must differ")

    allok = True
    with tempfile.TemporaryDirectory(prefix="cft_ckpt_") as workdir:
        for fmt in formats:
            for seq, extra, controls in SEQUENCES:
                ok, detail = one_sequence(exe, fmt, workdir, seq, extra, controls)
                print("  %-6s %-6s %s  %s" % (fmt, seq, "PASS" if ok else "FAIL", detail))
                allok = allok and ok

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
