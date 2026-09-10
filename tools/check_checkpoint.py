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
#   gate_real --fpNNN --save   ckpt.bin     -> stops half way
#   gate_real --fpNNN --resume ckpt.bin     -> the continued dump
#
# and requires the two dumps to be identical. The dumps are exact hex
# floats (%a), so the comparison is over bits rather than over a decimal
# rendering of them: a last-bit difference shows as a different line.
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


def one_format(exe, fmt, workdir):
    ckpt = os.path.join(workdir, "ckpt_%s.bin" % fmt)
    flag = "--" + fmt

    ok, straight, err = run(exe, [flag, "--straight"], True)
    if not ok:
        return False, "the straight run failed: %s" % err.strip().splitlines()[-1:]

    ok, _, err = run(exe, [flag, "--save", ckpt], False)
    if not ok:
        return False, "the save failed: %s" % err.strip()
    if not os.path.exists(ckpt):
        return False, "no archive at %s" % ckpt
    size = os.path.getsize(ckpt)

    ok, resumed, err2 = run(exe, [flag, "--resume", ckpt], True)
    if not ok:
        return False, "the resume failed: %s" % err2.strip()

    # What the two side processes said about themselves, for the record.
    note = " | ".join(s.strip() for s in (err.strip(), err2.strip()) if s.strip())

    if straight != resumed:
        return False, "the continued run DIFFERS - %s" % first_difference(straight, resumed)

    # The control. Same three processes, except the resume is told not
    # to carry the loaded state into the engine. If THAT still matches,
    # the comparison above was not testing what it claims to test.
    ok, control, _ = run(exe, [flag, "--resume", ckpt], True,
                         {"CFT_REBOUND_NO_ADOPT": "1"})
    if not ok:
        return False, "the control run would not start"
    if control == straight:
        return False, ("NOT SENSITIVE - a resume that discards the loaded "
                       "state still matches, so this comparison proves nothing")

    return True, "%d values identical, control differs, archive %d B%s" % (
        len(straight.splitlines()), size, (" | " + note) if note else "")


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
    print("check_checkpoint: three separate processes per format, "
          "on %s" % (art if art else "the software backend"))
    print("                  each format also runs a control that discards the "
          "loaded state and must differ")

    allok = True
    with tempfile.TemporaryDirectory(prefix="cft_ckpt_") as workdir:
        for fmt in formats:
            ok, detail = one_format(exe, fmt, workdir)
            print("  %-6s %s  %s" % (fmt, "PASS" if ok else "FAIL", detail))
            allok = allok and ok

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
