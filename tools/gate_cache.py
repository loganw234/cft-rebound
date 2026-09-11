#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
"""Skip a suite leg whose inputs have not changed since it passed.

The full suite is about forty minutes, which is the wrong tool for
"I edited a docstring, did I break anything". This makes `make
check-quick` skip a leg whose inputs are byte-for-byte what they were
when that leg last passed, and leaves `make check` alone.

WHY IT KEYS ON THE BINARY AND NOT THE SOURCES
---------------------------------------------
The obvious design hashes the source files a leg "regards". It is the
fragile one: it has to enumerate every header, every compiler flag,
every data file and both pinned upstreams, and the failure mode of
missing one is skipping a leg that would have failed - a silent wrong
answer, which is the single outcome this repository is built to
prevent.

So the key is over the ARTIFACTS a leg actually executes, plus its own
script and the data it reads. A binary is the closed-over result of
every source, header and flag that went into it: if `check_dropin.exe`
is byte-identical, nothing that could change its behaviour changed. The
build system already decides when to relink; this trusts that decision
rather than re-deriving it badly.

It works here because every verdict in this suite is bit-identity -
no tolerance, no timing, no seed. A rerun of the same bytes over the
same inputs cannot answer differently. A suite whose assertions were
statistical could not use this.

WHAT IS IN THE KEY
------------------
  - the bytes of every --input (binaries the leg runs, its own script,
    data files it reads);
  - the value of every --env (CFT_REBOUND_ARTIFACT above all: a pass
    from the software backend must NEVER satisfy a card run);
  - the command line itself, so `--fp64` and `--fp256` are different
    legs of the same binary.

Name every binary a leg EXECUTES, not only the one it is. gate_subprocess
spawns ias15_cft, and a key that omitted it would skip after the program
it drives had changed.

WHAT KEEPS IT HONEST
--------------------
  - `make check` never consults this. "The full suite passed" keeps
    meaning what it has always meant.
  - Every skip prints itself, by name and key. A silently skipped leg
    is a gate that cannot fail, and this repository has found seven of
    those in one day.
  - A stamp records whether the pass that wrote it was full or quick,
    and `--report` says so. A quick pass can save you a quick rerun; it
    can never be mistaken for a full one, because nothing reads these
    stamps except the quick run.
  - The cache lives under the build tree, so it cannot travel between
    build/ and build2/, or between machines.

  tools/gate_cache.py --leg NAME --cache-dir DIR --mode quick|full \\
      --input FILE [--input FILE ...] [--env VAR ...] -- CMD [ARGS ...]
  tools/gate_cache.py --report --cache-dir DIR
"""
import argparse
import hashlib
import json
import os
import subprocess
import sys
import time


def key_of(leg, inputs, envs, cmd):
    """SHA-256 over everything that could change this leg's verdict.

    A missing input is fatal rather than skipped: it means the caller
    named something that is not there, and hashing "absent" would make
    two different states share a key."""
    h = hashlib.sha256()
    h.update(b"cft-rebound gate_cache v1\0")
    h.update(leg.encode() + b"\0")
    for path in sorted(inputs):
        if not os.path.isfile(path):
            sys.exit("gate_cache: %s: no such input for leg %s" % (path, leg))
        h.update(path.replace("\\", "/").encode() + b"\0")
        with open(path, "rb") as f:
            for chunk in iter(lambda: f.read(1 << 20), b""):
                h.update(chunk)
        h.update(b"\0")
    for var in sorted(envs):
        h.update(var.encode() + b"=" + os.environ.get(var, "").encode() + b"\0")
    for part in cmd:
        h.update(part.encode() + b"\0")
    return h.hexdigest()


def stamp_path(cache_dir, leg):
    safe = "".join(c if c.isalnum() or c in "-_." else "_" for c in leg)
    return os.path.join(cache_dir, safe + ".json")


def report(cache_dir):
    if not os.path.isdir(cache_dir):
        print("gate cache: empty (no leg has passed in this build tree yet)")
        return 0
    rows = []
    for name in sorted(os.listdir(cache_dir)):
        if not name.endswith(".json"):
            continue
        try:
            with open(os.path.join(cache_dir, name), encoding="utf-8") as f:
                rows.append(json.load(f))
        except (OSError, ValueError):
            continue
    if not rows:
        print("gate cache: empty")
        return 0
    full = sum(1 for r in rows if r.get("mode") == "full")
    print("gate cache: %d leg(s) recorded, %d from a full run, %d from a quick one"
          % (len(rows), full, len(rows) - full))
    for r in sorted(rows, key=lambda r: r.get("leg", "")):
        print("  %-26s %s  %s  %s" % (r.get("leg", "?"), r.get("key", "")[:12],
                                      r.get("mode", "?"), r.get("when", "")))
    return 0


def main():
    ap = argparse.ArgumentParser(add_help=False)
    ap.add_argument("--leg")
    ap.add_argument("--cache-dir", required=True)
    ap.add_argument("--mode", choices=("quick", "full"), default="full")
    ap.add_argument("--input", action="append", default=[])
    ap.add_argument("--env", action="append", default=[])
    ap.add_argument("--report", action="store_true")
    ap.add_argument("cmd", nargs=argparse.REMAINDER)
    a = ap.parse_args()

    if a.report:
        return report(a.cache_dir)

    if not a.leg:
        sys.exit("gate_cache: --leg is required unless --report")
    cmd = a.cmd[1:] if a.cmd and a.cmd[0] == "--" else a.cmd
    if not cmd:
        sys.exit("gate_cache: no command given after --")

    key = key_of(a.leg, a.input, a.env, cmd)
    sp = stamp_path(a.cache_dir, a.leg)

    # Only the quick run reads a stamp. The full run is the thing whose
    # meaning must not change, so it always executes.
    if a.mode == "quick" and os.path.isfile(sp):
        try:
            with open(sp, encoding="utf-8") as f:
                prev = json.load(f)
        except (OSError, ValueError):
            prev = {}
        if prev.get("key") == key:
            print("%s: SKIPPED - inputs unchanged since it passed (%s, %s run)"
                  % (a.leg, key[:12], prev.get("mode", "?")))
            sys.stdout.flush()
            return 0

    rc = subprocess.call(cmd)
    if rc == 0:
        os.makedirs(a.cache_dir, exist_ok=True)
        tmp = sp + ".tmp"
        with open(tmp, "w", encoding="utf-8") as f:
            json.dump({"leg": a.leg, "key": key, "mode": a.mode,
                       "when": time.strftime("%Y-%m-%d %H:%M:%S"),
                       "cmd": cmd, "inputs": sorted(a.input),
                       "env": sorted(a.env)}, f, indent=1)
        os.replace(tmp, sp)      # never a half-written stamp
    else:
        # A failure retires the stamp. Otherwise a leg that failed today
        # could be skipped tomorrow on yesterday's pass.
        try:
            os.remove(sp)
        except OSError:
            pass
    return rc


if __name__ == "__main__":
    sys.exit(main())
