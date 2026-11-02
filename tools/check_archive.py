#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The Simulationarchive gates, all short. In order:
#
#   1. gate_restart  an archive written mid-run, restarted from, and the
#                    continuation compared BIT FOR BIT with an
#                    uninterrupted run - at binary64, binary128 and
#                    binary256, from a single-snapshot file and from the
#                    second snapshot of a two-snapshot one (which
#                    REBOUND stores as a diff).
#   2. gate_write    a cft archive at binary128, plus the exact binary64
#                    view a reader must recover.
#      gate_stock    that archive opened by a program linked against the
#                    pinned upstream librebound and NOTHING of
#                    cft-rebound's: no integrator registered, no cft_
#                    descriptor in the process. Its warning bits and the
#                    recovered state are printed verbatim.
#   3. gate_promote  a stock REBOUND archive opened here: promoted, and
#                    the promotion checked to be exact against REBOUND's
#                    own numbers; then a binary256 archive refused by a
#                    binary64 run, loudly.
#
# The stub integrator these four drive keeps no wide clock or steps, so
# its archives carry none of the fields added on 2026-11-02 (the clock
# family and the identity; src/cft_archive.c, NOTE 6), exactly like
# every archive written before that day. Their loads must therefore say
# "restored exactly, its clock and steps from binary64", and gate_restart
# and gate_promote hold them to it. The real integrator's archives, which
# carry those fields, are tests/gate_real.c's and
# tools/check_checkpoint.py's.
#
# COMPATIBILITY ACROSS 2026-11-02, measured on 2026-11-02 with af2f400's
# own gate_real against this tree's, at fp64, fp128 and fp256, in three
# processes each. It is a measurement and not a leg; the legs that hold
# the first half are gate_real's old-archive case and check_checkpoint's
# "old" sequence.
#
#   An archive af2f400 wrote, loaded here: 6 of 6 (plain and regrow at
#     each format) load, and say "restored exactly, its clock and steps
#     from binary64". Each resumes identically to a new archive resumed
#     with CFT_REBOUND_NO_CLOCK and CFT_REBOUND_NO_STEP - the clock and
#     steps rebuilt from binary64, as every load did then - and differs
#     from the straight run in the clock lines alone (plain) or not at
#     all (regrow, whose add rebuilds the clock in both).
#     gate_real --old-list writes such an archive. Field for field it is
#     af2f400's own: the same header, the same 130 fields (174 for
#     regrow) in the same order, the same bytes in all but particles (the
#     struct's three pointers only), walltime, walltime_last_steps and
#     rand_seed - the four fields two af2f400 runs differ in too.
#
#   An archive this tree writes, loaded by af2f400: 9 of 9 (plain, regrow
#     and adaptive at each format) load and continue. That build's reader
#     stops at cft_t, the first field it does not know, and REBOUND says
#     so on stderr:
#         Error! Could not find field descriptor for name.
#         Warning! Unknown field found in binary file.
#     Everything it knows comes before cft_t, and its finish_load repairs
#     cft_n_elem, the one field after, from the file, so it reports
#     "restored exactly" and continues as it continues its own archives:
#     identical to its own resume of its own archive, and to its straight
#     run, on every line it prints. THE STATED LIMIT: such a load prints
#     that Error and that Warning, and gets only what that build ever
#     got - the clock and steps rebuilt from binary64, so an adaptive
#     resume above binary64 diverges there as it always did (measured at
#     fp128 and fp256; at fp64 it does not).
#
#   python tools/check_archive.py [--build build] [--keep]

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="build")
    ap.add_argument("--keep", action="store_true",
                    help="do not delete the archives the gates wrote")
    args = ap.parse_args()

    b = args.build if os.path.isabs(args.build) else os.path.join(root, args.build)
    exe = ".exe" if os.name == "nt" else ""
    work = tempfile.mkdtemp(prefix="cft_archive_")

    gates = [("gate_restart", []), ("gate_write", []), ("gate_stock", []),
             ("gate_promote", [])]
    failed = []
    try:
        for name, extra in gates:
            path = os.path.join(b, name + exe)
            if not os.path.exists(path):
                print("MISSING %s - run `make archive` first" % path)
                failed.append(name)
                continue
            p = subprocess.run([path, work] + extra, cwd=root,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
            sys.stdout.write(p.stdout)
            if p.returncode != 0:
                failed.append("%s (exit %d)" % (name, p.returncode))
    finally:
        if args.keep:
            print("archives kept in %s" % work)
        else:
            shutil.rmtree(work, ignore_errors=True)

    if failed:
        print("check_archive: FAIL - " + ", ".join(failed))
        return 1
    print("check_archive: all archive gates pass")
    return 0


if __name__ == "__main__":
    sys.exit(main())
