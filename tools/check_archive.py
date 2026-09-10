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
