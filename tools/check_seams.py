#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# Round 3's seam tests. Each holds two parcels together, where neither
# parcel's own gate looks (ROUND3.md, "What the integrator keeps").
#
# A ROUTE ACROSS A CHECKPOINT. A route never changes a record line (the
# door's claim, held by each routed parcel's gate). A checkpoint resumes
# bit for bit (P3's claim, held by tools/check_checkpoint.py). Each gate
# holds its claim with the other's variable at its default. This one
# crosses them, with three processes as check_checkpoint's sequences use:
#
#   gate_real --fpNNN --straight                  -> the reference dump
#   gate_real --fpNNN --save   ckpt.bin           (the route forced one way)
#   gate_real --fpNNN --resume ckpt.bin           (the route forced the other)
#
# The resumed dump must still be the straight run's, in both directions.
# Two defects show here and nowhere else: a route's state that an
# archive does not carry, and an archive whose restore depends on the
# route that wrote it.
#
#   P2 x P3   reduce_seg, the corrector's device exit test.
#   P1 x P3   indexed and pair_gather, gravity's routes. They join ROUTES
#             when P1 merges.
#
# ROUND3.md named the P2 x P3 seam "a replay across a checkpoint". The two
# do not meet. The replay is the standalone program's (main()'s parser),
# and the checkpoint is the drop-in's Simulationarchive. The route is what
# they share, so the route is the seam.
#
# THE PHYSICS COUNTER. The save and the resume run with --routes, so each
# says which route every family took in its process and what that route
# carried. The forced route must be the one taken, by=forced, and it must
# have carried the run: calls and no fallback when on, fallbacks and no
# call when off. A crossing whose route carried nothing compares two runs
# of today's route and passes for the wrong reason.
#
# The control: the same resume with CFT_REBOUND_NO_ADOPT, the shim
# discarding what it loaded, must differ from the straight run. If it
# does not, the comparison proves nothing.
#
# The sequences, the process runner and the comparison are
# check_checkpoint's own, imported, never copied. A new sequence there,
# or a fix to how a phase is run, reaches this gate unedited.
#
#   python tools/check_seams.py [--build build] [--quick]
#   (--quick: binary64 only, the light tiers)

import argparse
import os
import re
import sys
import tempfile

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
sys.path.insert(0, here)
import check_checkpoint as CC  # noqa: E402

# The routes a merged parcel carries. Add a name here when its parcel
# merges.
ROUTES = ["reduce_seg"]
SEQUENCES = [s for s in CC.SEQUENCES if s[0] in ("plain", "adaptive")]

# What gate_real --routes prints, one line per family (tests/gate_real.c,
# print_routes()).
ROUTE_LINE = re.compile(r"^gate_real: route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)\s*$", re.M)


def carried(err, route, forced):
    """Whether a phase took the route it was forced to and that route
    carried the run, with a phrase that says what it carried."""
    for m in ROUTE_LINE.finditer(err):
        if m.group(1) != route:
            continue
        on, by, calls, falls = m.group(2) == "on", m.group(3), int(m.group(4)), int(m.group(5))
        if on != (forced == "on") or by != "forced":
            return False, "%s=%s by=%s, though it was forced %s" % (route, m.group(2), by, forced)
        if on and (calls == 0 or falls != 0):
            return False, "forced on, it carried %d calls and %d fallbacks" % (calls, falls)
        if not on and (falls == 0 or calls != 0):
            return False, "forced off, it carried %d calls and %d fallbacks" % (calls, falls)
        return True, ("%d calls on the new route" % calls) if on else ("%d fallbacks to today's" % falls)
    return False, "no route line names %s" % route


def main():
    ap = argparse.ArgumentParser(description="round 3's seam tests: a route across a checkpoint")
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--quick", action="store_true", help="binary64 only")
    a = ap.parse_args()
    exe = os.path.join(a.build, "gate_real")
    formats = ["fp64"] if a.quick else CC.FORMATS
    allok = True
    crossings = 0

    def verdict(ok, text):
        nonlocal allok
        allok &= ok
        print("  %s %s" % ("ok  " if ok else "FAIL", text))

    print("check_seams: a route across a checkpoint, routes %s, formats %s, sequences %s"
          % (", ".join(ROUTES), ", ".join(formats), ", ".join(s[0] for s in SEQUENCES)))
    with tempfile.TemporaryDirectory(prefix="cft_seams_") as wd:
        for fmt in formats:
            flag = "--" + fmt
            for seq, extra, _ in SEQUENCES:
                ok, straight, err = CC.run(exe, [flag] + extra + ["--straight"], True)
                if not ok:
                    verdict(False, "%s %s: the straight run failed: %s" % (fmt, seq, err.strip()[-200:]))
                    continue
                for route in ROUTES:
                    var = "CFT_REBOUND_ROUTE_" + route.upper()
                    for wrote, read in (("on", "off"), ("off", "on")):
                        ckpt = os.path.join(wd, "seam_%s_%s_%s_%s.bin" % (fmt, seq, route, wrote))
                        ok, _, err = CC.run(exe, [flag] + extra + ["--routes", "--save", ckpt], False, {var: wrote})
                        if not ok or not os.path.exists(ckpt):
                            verdict(False, "%s %s: the save with %s=%s failed: %s"
                                    % (fmt, seq, var, wrote, err.strip()[-200:]))
                            continue
                        okw, at_save = carried(err, route, wrote)
                        ok, resumed, err = CC.run(exe, [flag] + extra + ["--routes", "--resume", ckpt], True,
                                                  {var: read})
                        if not ok:
                            verdict(False, "%s %s: the resume with %s=%s failed: %s"
                                    % (fmt, seq, var, read, err.strip()[-200:]))
                            continue
                        okr, at_resume = carried(err, route, read)
                        okc, control, _ = CC.run(exe, [flag] + extra + ["--resume", ckpt], True,
                                                 {var: read, "CFT_REBOUND_NO_ADOPT": "1"})
                        sensitive = okc and control != straight
                        same = resumed == straight
                        if okw and okr:
                            crossings += 1
                        verdict(same and sensitive and okw and okr,
                                "%s %s: %s %s at the save (%s), %s at the resume (%s): %s; %s"
                                % (fmt, seq, route, wrote, at_save, read, at_resume,
                                   "%d values identical to the straight run" % len(straight.splitlines())
                                   if same else "the resume DIFFERS - %s" % CC.first_difference(straight, resumed),
                                   "CFT_REBOUND_NO_ADOPT differs, so the comparison bites" if sensitive
                                   else "NOT SENSITIVE: a resume that discards the archive still matches"))
    # A gate that crossed nothing passes nothing: every crossing's route
    # must have carried both of its halves.
    want = len(formats) * len(SEQUENCES) * len(ROUTES) * 2
    verdict(crossings == want, "%d route crossings carried by the route on both sides, of %d" % (crossings, want))
    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
