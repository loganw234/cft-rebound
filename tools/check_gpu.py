#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# This tree on a GPU as it is, held to the software backend gate by gate
# (ROUND5.md, P1; ROADMAP's GPU half, step 2). `make HIP=1 check-gpu` runs
# it on a GPU host. It is in none of the desktop's three tiers: the desktop
# has no GPU, and there it refuses by name.
#
#   python tools/check_gpu.py [--build build] [--artifact hip:0] [--hip 1]
#                             [--legs a,b,...] [--skip a,b,...] [--no-legs]
#                             [--keep] [--list]
#
# 0. THE DEVICE. ias15_cft --probe opens the artifact with CFT_HIP_INFO=1. A
#    host with no such device, or a build whose libcft has no GPU backend, is
#    REFUSED by name (rc 2) and nothing else runs. With --hip 1 (the
#    Makefile's HIP) a libcft.a without backend_hip.o is named as V0's
#    stale-library trap. The probe's one device line is the line every later
#    process must print; the routes the device does not publish are read
#    from the probe, never typed.
#
# 1. THE WITNESS. A record's `backend=hip:0` only echoes the artifact string
#    the run was given (src/ias15_cft.c, main(); V0, entry 54), so it proves
#    nothing. Each leg runs against a shadow of the build directory, its
#    top-level files linked, in which ias15_cft and check_dropin are a
#    witness: a script that runs the real program with the caller's
#    arguments, environment and stdin, passes its output back as it comes,
#    and keeps it. When the run was on a device (CFT_REBOUND_ARTIFACT, or
#    --artifact, set) it runs the same arguments once more on the software
#    backend, its output files redirected, and keeps that too. So for every
#    process a leg starts, whatever spawns it (a gate, compare_run.py,
#    extract_particle.py, check_corrector's check_dropin), this driver holds:
#    - its device line: at least one, and every one the probe's. A process
#      that ran the software backend while its leg was pointed at the device
#      is a STAND-IN and FAILs its leg by name, with its arguments;
#    - its record: the software backend's from the same arguments, every
#      line, every counter and every output file (--dt-out, --attempts-out,
#      --record-iters) byte for byte. Set aside, and counted: the header's
#      backend= value, the trailer's timings (seconds, steps_per_s,
#      wall_seconds, the t_ buckets, system_steps_per_s), and the by= of a
#      route the device does not publish (absent on the device, opt_in in
#      software; entry 54). stderr is compared too, device lines and the
#      `done:` line's seconds set aside, and a difference there is printed;
#    - its physics, beside the bits: steps, rejections, passes, flags.
#    The leg's own verdict (rc 0, RESULT: PASS or check_dropin's "every case
#    passed", no FAIL line) is required as well. So a leg passes on the GPU
#    only if it passes, every process it started ran on the device, and every
#    record it made is the software backend's.
#
#    THE CONTROLS, permanent, each a real run through the witness: a record
#    with one hex digit planted after the run must DIFFER by line and column;
#    a run whose environment names no artifact must FAIL as a stand-in; a run
#    without CFT_HIP_INFO must FAIL for its missing device line; and a device
#    line changed by one character must FAIL as not the probe's.
#
# 2. THE ROUTES (V0, entry 54), at every format. Forced on, `scalar` and
#    `reduce_seg` run on the device in the software backend's bits: every
#    line of the forced run is the unforced software run's but the route's
#    own, and the operations the software run carried as fallbacks the
#    device carried as calls, with no fallback. Forced on, every route the
#    device does not publish is refused by name, where the software backend
#    runs the same command. --door-selftest on the device gives the summary
#    the software backend gives with those routes hidden (116 identical and
#    48 not compared at binary64 at this commit, derived, not typed), and
#    every case line the same, `absent` for `hidden`.
#
# 3. THE GPU SEAM (ROUND4.md, "What the integrator keeps"; ROUND5.md).
#    tools/make_ring.py's two-wave ring, its zonal line kept, recorded on the
#    device from a first step far too large, so that both waves reject and
#    are sampled. tools/compare_run.py, run with no artifact, holds that GPU
#    record against this tree's CPU reruns: both members IDENTICAL, and
#    their particles. Then one hex digit of one column on one line of wave 0
#    is planted, a column at a time, on more than one line (the exact-time
#    pair among them, which the host's augmented sum makes): each must FAIL
#    wave 0 naming that line and column, while wave 1 stays IDENTICAL. The
#    record's columns are read here from its own `# columns:` line, apart
#    from compare_run.py's reader.
#
# 4. THE LEGS, each in its quick form as the Makefile's check-quick runs it
#    (check_dropin as check-light runs it, --light), the commands read from
#    the Makefile and never retyped: check_records, check_equivalence and
#    check_program_engine at binary64; check_ensemble, check_ring,
#    check_replay and check_compare_run; check_zonal and check_reversal;
#    check_corrector (with --library, check_dropin --light with the exit
#    test forced on) and check_dropin, the registered integrator's.
#
# Every verdict carries the device line captured from the run. Exit 0 PASS,
# 1 FAIL, 2 REFUSED. The leg logs and the witness's files are kept under
# <build>/check_gpu/<pid>/ after a FAIL or with --keep, and removed after a
# PASS.
#
# Not through the gate cache: a GPU run is its own evidence, and a cached
# pass would carry no device line.

import argparse
import json
import os
import platform
import re
import shutil
import socket
import subprocess
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)
EXE = ".exe" if os.name == "nt" else ""

# The legs: (the Makefile target whose recipe holds the command, the leg's
# name there). The order puts the long ones last.
GPU_LEGS = [
    ("check-quick", "check_records"),
    ("check-quick", "check_equivalence"),
    ("check-quick", "check_program_engine"),
    ("check-quick", "check_ensemble"),
    ("check-quick", "check_ring"),
    ("check-quick", "check_replay"),
    ("check-quick", "check_compare_run"),
    ("check-quick", "check_zonal"),
    ("check-quick", "check_reversal"),
    ("check-quick", "check_corrector"),
    ("check-light", "check_dropin_light"),
]
# The programs that open the engine, and so are witnessed. ias15_ref is
# REBOUND's own and opens no device.
WITNESSED = ("ias15_cft", "check_dropin")
OUT_FLAGS = ("--dt-out", "--attempts-out", "--record-iters")
DEVICE_LINE = "cft hip: "
TIMING = re.compile(r" (?:seconds|steps_per_s|wall_seconds|system_steps_per_s|t_[a-z_]+)=\S+")

FAILS = []
OKS = [0]


def verdict(ok, text, detail=()):
    print("  %s %s" % ("ok  " if ok else "FAIL", text))
    for d in detail:
        print("       %s" % d)
    if ok:
        OKS[0] += 1
    else:
        FAILS.append(text)
    sys.stdout.flush()
    return ok


class Refused(Exception):
    pass


def built(path):
    try:
        return time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(path)))
    except OSError:
        return "MISSING"


def short(text, n=300):
    text = (text or "").strip()
    return text if len(text) <= n else "..." + text[-n:]


# ---------------------------------------------------------------- 0. the device

def probe(real, art, info=True):
    env = {k: v for k, v in os.environ.items() if not k.upper().startswith("CFT_REBOUND_")}
    env["CFT_REBOUND_ARTIFACT"] = art
    if info:
        env["CFT_HIP_INFO"] = "1"
    p = subprocess.run([real, "--probe"], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env,
                       cwd=root, stdin=subprocess.DEVNULL)
    caps = {}
    for line in p.stdout.splitlines():
        k, sep, v = line.partition(": ")
        if sep:
            caps.setdefault(k, v)
    routes = {}
    for m in re.finditer(r"^route (\w+): (\S+) .*?(published|NOT published)", p.stdout, re.M):
        routes[m.group(1)] = (m.group(2), m.group(3) == "published")
    devlines = [l for l in p.stderr.splitlines() if l.startswith(DEVICE_LINE)]
    return p.returncode, p.stdout, p.stderr, caps, routes, devlines


def feature_names(bits):
    """The cft_caps.seq_features bits in a mask, by their names in the pinned
    cft.h, read rather than typed."""
    names = []
    try:
        with open(os.path.join(root, "third_party", "cft-fp256", "host", "include", "cft.h")) as f:
            for m in re.finditer(r"^#define (CFT_(?:SEQ_FEAT|FEAT|ALU_EXT)_\w+)\s+0x([0-9a-fA-F]+)u", f.read(), re.M):
                if bits & int(m.group(2), 16):
                    names.append(m.group(1))
    except OSError:
        pass
    return " and ".join(names) or "?"


def libcft_has_hip(lib):
    """Whether libcft.a carries the GPU backend's object, by its member list."""
    try:
        p = subprocess.run(["ar", "t", lib], stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    except OSError:
        return None
    if p.returncode != 0:
        return None
    return "backend_hip.o" in p.stdout.split()


# ---------------------------------------------------------------- 1. the witness

WRAPPER = r'''#!%(python)s -B
# check_gpu's witness around %(real)s.
# tools/check_gpu.py writes this file into its shadow build directory; it is
# not a tool of its own. It runs the real program with the caller's arguments,
# environment and stdin, passes its output back as it comes, and keeps it; and
# when the run was on a device, it runs the same arguments once more on the
# software backend, output files redirected, and keeps that too.
import json, os, subprocess, sys, threading, time
REAL = %(real)r
SHADOW = %(shadow)r
OUT_FLAGS = %(out_flags)r


def pump(src, sink, keep):
    while True:
        b = src.read1(65536)
        if not b:
            break
        keep.append(b)
        if sink is not None:
            sink.write(b)
            sink.flush()


def run(argv, env, data, live):
    p = subprocess.Popen([REAL] + argv, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, env=env)
    out, err = [], []
    ts = [threading.Thread(target=pump, args=(p.stdout, sys.stdout.buffer if live else None, out)),
          threading.Thread(target=pump, args=(p.stderr, sys.stderr.buffer if live else None, err))]
    for t in ts:
        t.start()
    try:
        if data:
            p.stdin.write(data)
        p.stdin.close()
    except (BrokenPipeError, OSError):
        pass
    for t in ts:
        t.join()
    return p.wait(), b"".join(out), b"".join(err)


def keep(d, name, data):
    with open(os.path.join(d, name), "wb") as f:
        f.write(data)


def main():
    argv = sys.argv[1:]
    data = b""
    try:
        if not os.isatty(0):
            data = sys.stdin.buffer.read()
    except (OSError, ValueError):
        data = b""
    env = dict(os.environ)
    art = None
    for i in range(len(argv) - 1):
        if argv[i] == "--artifact":
            art = argv[i + 1]
    if art is None:
        art = env.get("CFT_REBOUND_ARTIFACT", "")
    t0 = time.time()
    rc, out, err = run(argv, env, data, True)
    t1 = time.time()
    try:
        with open(os.path.join(SHADOW, ".leg")) as f:
            d = f.read().strip()
    except OSError:
        d = ""
    d = d or os.path.join(SHADOW, "unattributed")
    os.makedirs(d, exist_ok=True)
    tag = "inv-%%d-%%d" %% (os.getpid(), time.time_ns())
    rec = {"real": REAL, "argv": argv, "cwd": os.getcwd(), "artifact": art, "hip_info": env.get("CFT_HIP_INFO"),
           "rc": rc, "seconds": t1 - t0, "stdin_bytes": len(data), "started": t0,
           "door_env": sorted("%%s=%%s" %% (k, v) for k, v in env.items()
                              if k.upper().startswith("CFT_REBOUND_") and k != "CFT_REBOUND_ARTIFACT")}
    keep(d, tag + ".dev.out", out)
    keep(d, tag + ".dev.err", err)
    outs = [(argv[i], argv[i + 1]) for i in range(len(argv) - 1) if argv[i] in OUT_FLAGS]
    rec["outputs"] = []
    for flag, path in outs:
        try:
            with open(path, "rb") as f:
                keep(d, tag + ".dev" + flag, f.read())
            rec["outputs"].append(flag)
        except OSError:
            rec["outputs"].append(flag + ":missing")
    if art:
        sargv, i = [], 0
        while i < len(argv):
            a = argv[i]
            if a == "--artifact" and i + 1 < len(argv):
                i += 2
                continue
            if a in OUT_FLAGS and i + 1 < len(argv):
                sargv += [a, os.path.join(d, tag + ".sw" + a)]
                i += 2
                continue
            sargv.append(a)
            i += 1
        senv = dict(env)
        senv["CFT_REBOUND_ARTIFACT"] = ""
        s0 = time.time()
        src, sout, serr = run(sargv, senv, data, False)
        rec.update({"sw_rc": src, "sw_seconds": time.time() - s0, "sw_argv": sargv})
        keep(d, tag + ".sw.out", sout)
        keep(d, tag + ".sw.err", serr)
    with open(os.path.join(d, tag + ".json.tmp"), "w") as f:
        json.dump(rec, f)
    os.replace(os.path.join(d, tag + ".json.tmp"), os.path.join(d, tag + ".json"))
    sys.exit(rc if rc >= 0 else 128 - rc)


main()
'''


class Shadow:
    """The build directory again, its top-level files linked, with a witness
    in place of each program that opens the engine."""

    def __init__(self, build, work):
        self.build = build
        self.dir = os.path.join(work, "build")
        os.makedirs(self.dir)
        for name in sorted(os.listdir(build)):
            src = os.path.join(build, name)
            if os.path.isfile(src) and os.path.splitext(name)[0] not in WITNESSED:
                os.symlink(src, os.path.join(self.dir, name))
        for prog in WITNESSED:
            real = os.path.join(build, prog + EXE)
            path = os.path.join(self.dir, prog + EXE)
            with open(path, "w") as f:
                f.write(WRAPPER % {"python": sys.executable, "real": real, "shadow": self.dir,
                                   "out_flags": OUT_FLAGS})
            os.chmod(path, 0o755)
            st = os.stat(real)
            os.utime(path, (st.st_atime, st.st_mtime))   # a leg that prints a binary's build time prints the real one's
        self.legs = os.path.join(work, "witness")

    def start(self, leg):
        d = os.path.join(self.legs, leg)
        shutil.rmtree(d, ignore_errors=True)
        os.makedirs(d)
        with open(os.path.join(self.dir, ".leg.tmp"), "w") as f:
            f.write(d + "\n")
        os.replace(os.path.join(self.dir, ".leg.tmp"), os.path.join(self.dir, ".leg"))
        return d


def read(path, binary=False):
    try:
        with open(path, "rb") as f:
            b = f.read()
    except OSError:
        return None
    return b if binary else b.decode("utf-8", "replace")


def norm_stdout(text, unpublished, side):
    """The lines compared, with what may differ set aside. Each line keeps
    its number, so a difference is named where it is."""
    keep = []
    rt = re.compile(r"\b(%s)=(on|off) by=(%s)\b" % ("|".join(map(re.escape, unpublished)) or "(?!)",
                                                    "absent" if side == "dev" else "opt_in"))
    routed = 0
    for line in text.splitlines():
        if line.startswith("# program="):
            line = re.sub(r" backend=\S+", " backend=*", line)
        elif line.startswith("# steps_done="):
            line = TIMING.sub("", line)
        line, n = rt.subn(r"\1=\2 by=~", line)
        routed += n
        keep.append(line)
    return keep, routed


def norm_stderr(text):
    out = []
    for line in text.splitlines():
        if line.startswith(DEVICE_LINE):
            continue
        out.append(re.sub(r"^(done: .*?), [0-9.]+ s, ", r"\1, ", line))
    return out


def first_diff(a, b):
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            return i, x, y
    if len(a) != len(b):
        i = min(len(a), len(b))
        return i, a[i] if i < len(a) else "<end>", b[i] if i < len(b) else "<end>"
    return None


def name_column(x, y):
    """Where two lines first differ: the field, by the record's column layout
    on a sample line (system s sample k step, t dt_next dt_last E, six a body,
    then | t_hi t_lo), and its two values."""
    xs, ys = x.split(), y.split()
    j = next((k for k, (p, q) in enumerate(zip(xs, ys)) if p != q), min(len(xs), len(ys)))
    vals = (xs[j] if j < len(xs) else "<end>", ys[j] if j < len(ys) else "<end>")
    off = 2 if xs[:1] == ["system"] else 0
    if xs[off:off + 1] != ["sample"]:
        return "field %d" % (j + 1), vals
    k = j - off - 3
    if "|" in xs and j > xs.index("|"):
        return ("t_hi", "t_lo")[min(j - xs.index("|") - 1, 1)], vals
    if k < 0:
        return ("sample", "step")[min(j - off - 1, 1)] if j > off else "sample", vals
    if k < 4:
        return ("t", "dt_next", "dt_last", "E")[k], vals
    b, c = divmod(k - 4, 6)
    return "%s of body %d" % (("x", "y", "z", "vx", "vy", "vz")[c], b), vals


def physics(text):
    """The counters of one record, from its trailer."""
    m = re.search(r"^# steps_done=.*$", text, re.M)
    if not m:
        return None
    kv = dict(t.partition("=")[::2] for t in m.group(0)[2:].split())
    try:
        e = int(kv.get("E", "1"))
        return {"steps": int(kv["steps_done"]), "rejected": int(kv["steps_rejected"]),
                "attempts": int(kv["steps_done"]) + int(kv["steps_rejected"]),
                "mean_passes": float(kv["mean_pc_iterations"]), "max_passes": int(kv["max_pc_iterations"]),
                "capped": int(kv["iterations_max_exceeded"]), "flags": int(kv["flags_seen"], 16),
                "calls": int(kv["calls"]), "wall": float(kv.get("wall_seconds", "0")), "E": e}
    except (KeyError, ValueError):
        return None


def analyse(d, probe_line, unpublished):
    """Every process a leg started, from what the witness kept."""
    invs = []
    for name in sorted(os.listdir(d)) if os.path.isdir(d) else []:
        if name.endswith(".json"):
            with open(os.path.join(d, name)) as f:
                rec = json.load(f)
            rec["tag"] = name[:-5]
            invs.append(rec)
    invs.sort(key=lambda r: r.get("started", 0))
    res = {"n": len(invs), "by_prog": {}, "device": 0, "software": [], "no_line": [], "bad_line": [],
           "lines": set(), "identical": 0, "routed": 0, "differ": [], "stderr": [], "phys": [], "dropin": []}
    for r in invs:
        prog = os.path.basename(r["real"])
        res["by_prog"][prog] = res["by_prog"].get(prog, 0) + 1
        dev_out = read(os.path.join(d, r["tag"] + ".dev.out")) or ""
        dev_err = read(os.path.join(d, r["tag"] + ".dev.err")) or ""
        what = "%s %s" % (prog, " ".join(r["argv"]))
        if not r["artifact"]:
            res["software"].append(what)
            continue
        res["device"] += 1
        lines = [l for l in dev_err.splitlines() if l.startswith(DEVICE_LINE)]
        res["lines"].update(lines)
        if not lines:
            res["no_line"].append("%s (CFT_HIP_INFO=%s, rc %d)" % (what, r.get("hip_info"), r["rc"]))
        bad = [l for l in lines if l != probe_line]
        if bad:
            res["bad_line"].append("%s: %s" % (what, bad[0]))
        sw_out = read(os.path.join(d, r["tag"] + ".sw.out")) or ""
        sw_err = read(os.path.join(d, r["tag"] + ".sw.err")) or ""
        a, ra = norm_stdout(dev_out, unpublished, "dev")
        b, rb = norm_stdout(sw_out, unpublished, "sw")
        why = []
        if r["rc"] != r.get("sw_rc"):
            why.append("rc %d on the device, %s in software" % (r["rc"], r.get("sw_rc")))
        fd = first_diff(a, b)
        if fd:
            i, x, y = fd
            col, (vx, vy) = name_column(x, y)
            why.append("stdout line %d, %s: device %s, software %s" % (i + 1, col, short(vx, 80), short(vy, 80)))
        for flag in r.get("outputs", []):
            if flag.endswith(":missing"):
                continue
            x = read(os.path.join(d, r["tag"] + ".dev" + flag), True)
            y = read(os.path.join(d, r["tag"] + ".sw" + flag), True)
            if x != y:
                fx = first_diff((x or b"").decode("utf-8", "replace").splitlines(),
                                (y or b"").decode("utf-8", "replace").splitlines())
                why.append("its %s file, line %s" % (flag, fx[0] + 1 if fx else "?"))
        if why:
            res["differ"].append("%s: %s" % (what, "; ".join(why)))
        else:
            res["identical"] += 1
            res["routed"] += ra
        if norm_stderr(dev_err) != norm_stderr(sw_err):
            fe = first_diff(norm_stderr(dev_err), norm_stderr(sw_err))
            res["stderr"].append("%s: stderr line %d: device %r, software %r"
                                 % (what, fe[0] + 1, short(fe[1], 120), short(fe[2], 120)))
        ph = physics(dev_out)
        if ph:
            res["phys"].append(ph)
        if prog.startswith("check_dropin"):
            res["dropin"].append((len(re.findall(r"^\s*ok\b", dev_out, re.M)),
                                  len(re.findall(r"^\s*FAIL\b", dev_out, re.M)),
                                  "every case passed" in dev_out))
    return res


def physics_text(res):
    ph = res["phys"]
    parts = []
    if ph:
        att = sum(p["attempts"] for p in ph)
        passes = sum(p["mean_passes"] * p["attempts"] for p in ph)
        flags = 0
        for p in ph:
            flags |= p["flags"]
        parts.append("%d records, %d steps and %d rejected (systems' first), %.2f passes an attempt (at most %d), "
                     "%d capped, flags 0x%02x, %d library calls in %.1f s on the device's own clock"
                     % (len(ph), sum(p["steps"] for p in ph), sum(p["rejected"] for p in ph),
                        passes / att if att else 0.0, max(p["max_passes"] for p in ph), sum(p["capped"] for p in ph),
                        flags, sum(p["calls"] for p in ph), sum(p["wall"] for p in ph)))
    for ok, fail, every in res["dropin"]:
        parts.append("check_dropin: %d ok, %d FAIL, %s" % (ok, fail, "every case passed" if every else
                                                            "NOT every case passed"))
    return "; ".join(parts) or "no record"


def witness_verdict(label, res, need_device=True):
    """The witness's half of a leg's verdict: (ok, text, detail)."""
    detail = []
    ok = True
    if need_device and res["device"] == 0:
        ok = False
        detail.append("no process of this leg ran on the device")
    if res["software"]:
        ok = False
        detail.append("%d STAND-IN process(es) ran the software backend, the leg's environment naming no artifact: %s"
                      % (len(res["software"]), short(res["software"][0], 240)))
    if res["no_line"]:
        ok = False
        detail.append("%d device process(es) printed no device line: %s" % (len(res["no_line"]),
                                                                          short(res["no_line"][0], 240)))
    if res["bad_line"]:
        ok = False
        detail.append("%d device line(s) not the probe's: %s" % (len(res["bad_line"]), short(res["bad_line"][0], 300)))
    if res["differ"]:
        ok = False
        detail.append("%d of %d device record(s) NOT the software backend's: %s"
                      % (len(res["differ"]), res["device"], short(res["differ"][0], 600)))
        for extra in res["differ"][1:4]:
            detail.append("  and %s" % short(extra, 300))
    progs = ", ".join("%s %d" % kv for kv in sorted(res["by_prog"].items()))
    text = ("%s: %d of %d process(es) on the device (%s), a device line each and every line the probe's; "
            "%d of %d record(s) identical to the software backend's from the same arguments%s"
            % (label, res["device"] - len(res["no_line"]), res["n"], progs or "none", res["identical"], res["device"],
               ", %d route line(s) by=absent against opt_in" % res["routed"] if res["routed"] else ""))
    for s in res["stderr"][:2]:
        detail.append("note: stderr differs, %s" % short(s, 300))
    return ok, text, detail


def run_witnessed(shadow, leg, cmd, env, log, errlog=None):
    """A leg's command, its stdout to log and its stderr with it unless
    errlog names a file of its own."""
    d = shadow.start(leg)
    t0 = time.time()
    with open(log, "w") as f, open(errlog or os.devnull, "w") as g:
        p = subprocess.run(cmd, stdout=f, stderr=g if errlog else subprocess.STDOUT, stdin=subprocess.DEVNULL,
                           env=env, cwd=root)
    return p.returncode, time.time() - t0, d


def controls(shadow, work, art, probe_line, unpublished, env_dev):
    """The witness's own controls: each must FAIL what it plants."""
    print("1. the witness's controls: each plant must be caught by name")
    kepler = os.path.join(root, "data", "problems", "kepler.txt")
    att = os.path.join(work, "control.att")
    cmd = [os.path.join(shadow.dir, "ias15_cft" + EXE), "--format", "fp64", "--problem", kepler, "--dt", "0.01",
           "--epsilon", "1e-9", "--steps", "6", "--sample", "3", "--quiet", "--attempts-out", att]
    rc, secs, d = run_witnessed(shadow, "control-clean", cmd, env_dev, os.path.join(work, "control-clean.log"))
    res = analyse(d, probe_line, unpublished)
    ok, text, detail = witness_verdict("the clean run", res)
    verdict(ok and rc == 0 and res["identical"] == 1, text + " (rc %d)" % rc, detail)
    if not res["identical"]:
        return
    # (a) one hex digit of the device's record, planted after the run
    pd = os.path.join(shadow.legs, "control-planted")
    shutil.copytree(d, pd)
    tag = next(n[:-5] for n in os.listdir(pd) if n.endswith(".json"))
    out = read(os.path.join(pd, tag + ".dev.out")).splitlines()
    i = max(k for k, l in enumerate(out) if l.startswith("sample "))
    toks = out[i].split()
    j = 3 + 4 + 6                                  # sample k step, t dt_next dt_last E, body 0's six, then body 1's x
    mant = toks[j].split("p")[0]
    toks[j] = mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + toks[j][len(mant):]
    out[i] = " ".join(toks)
    with open(os.path.join(pd, tag + ".dev.out"), "w", newline="\n") as f:
        f.write("\n".join(out) + "\n")
    pres = analyse(pd, probe_line, unpublished)
    hit = pres["differ"] and ("stdout line %d, x of body 1: device %s" % (i + 1, toks[j])) in pres["differ"][0]
    verdict(hit and not witness_verdict("", pres)[0],
            "control a: one hex digit of body 1's x on the record's line %d, planted in the device's output: %s"
            % (i + 1, "caught: " + pres["differ"][0].split(": ", 1)[1] if pres["differ"] else "NOT CAUGHT"))
    # (b) one byte of the device's attempts record
    with open(os.path.join(pd, tag + ".dev.out"), "w", newline="\n") as f:
        f.write(read(os.path.join(d, tag + ".dev.out")))
    a = read(os.path.join(pd, tag + ".dev--attempts-out")).splitlines()
    k = max(n for n, l in enumerate(a) if not l.startswith("#") and l.endswith(" a"))
    a[k] = a[k][:-1] + "r"
    with open(os.path.join(pd, tag + ".dev--attempts-out"), "w", newline="\n") as f:
        f.write("\n".join(a) + "\n")
    pres = analyse(pd, probe_line, unpublished)
    hit = pres["differ"] and ("its --attempts-out file, line %d" % (k + 1)) in pres["differ"][0]
    verdict(hit, "control b: an attempt's outcome turned from a to r on line %d of the device's attempts record: %s"
            % (k + 1, "caught" if hit else "NOT CAUGHT"))
    # (c) a device line one character off
    with open(os.path.join(pd, tag + ".dev--attempts-out"), "w", newline="\n") as f:
        f.write("\n".join(read(os.path.join(d, tag + ".dev--attempts-out")).splitlines()) + "\n")
    err = read(os.path.join(pd, tag + ".dev.err"))
    bad = err.replace(probe_line, probe_line[:-1] + chr(ord(probe_line[-1]) ^ 1), 1)
    with open(os.path.join(pd, tag + ".dev.err"), "w", newline="\n") as f:
        f.write(bad)
    pres = analyse(pd, probe_line, unpublished)
    verdict(len(pres["bad_line"]) == 1 and not witness_verdict("", pres)[0],
            "control c: the device line's last character changed: %s"
            % ("caught as not the probe's" if pres["bad_line"] else "NOT CAUGHT"))
    # (d) a stand-in: the same run with no artifact in its environment
    env_sw = dict(env_dev)
    env_sw["CFT_REBOUND_ARTIFACT"] = ""
    rc, secs, d2 = run_witnessed(shadow, "control-standin", cmd, env_sw, os.path.join(work, "control-standin.log"))
    sres = analyse(d2, probe_line, unpublished)
    sok = witness_verdict("", sres)[0]
    verdict(rc == 0 and len(sres["software"]) == 1 and not sok,
            "control d: the run with CFT_REBOUND_ARTIFACT empty, a leg that quietly ran the software backend: rc %d, "
            "%s" % (rc, "caught as a stand-in" if sres["software"] and not sok else "NOT CAUGHT"))
    # (e) a device run with no device line: CFT_HIP_INFO unset
    env_ni = {k: v for k, v in env_dev.items() if k != "CFT_HIP_INFO"}
    rc, secs, d3 = run_witnessed(shadow, "control-noinfo", cmd, env_ni, os.path.join(work, "control-noinfo.log"))
    nres = analyse(d3, probe_line, unpublished)
    verdict(rc == 0 and len(nres["no_line"]) == 1 and not witness_verdict("", nres)[0],
            "control e: the run without CFT_HIP_INFO, on the device with no evidence of it: rc %d, %s"
            % (rc, "caught for its missing device line" if nres["no_line"] else "NOT CAUGHT"))


# ---------------------------------------------------------------- 2. the routes

def kvs(line):
    return dict(t.partition("=")[::2] for t in line.split())


def route_lines(text):
    out = {}
    for m in re.finditer(r"^# route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)$", text, re.M):
        out[m.group(1)] = (m.group(2), m.group(3), int(m.group(4)), int(m.group(5)))
    return out


def run_env(cmd, env):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env, cwd=root,
                       stdin=subprocess.DEVNULL)
    return p.returncode, p.stdout, p.stderr


def routes(real, art, probe_line, unpublished, published, envs, formats):
    env_dev, env_sw = envs
    print("2. the routes on %s, at every format (V0, entry 54)" % art)
    # The door's own names and variables, from the program's self-test lines.
    rc, out, err = run_env([real, "--door-selftest", "--format", "fp64"], env_sw)
    door = {m.group(1): (m.group(2), m.group(3)) for m in
            re.finditer(r"^door route (\w+)=\w+ by=\w+ env=(\S+) bit=(\S+)$", out, re.M)}
    if not door:
        verdict(False, "--door-selftest named no route (rc %d): %s" % (rc, short(err)))
        return
    probs = {"fp64": ("kepler.txt", "0.01", "8"), "fp128": ("kepler.txt", "0.01", "6"), "fp256": ("kepler.txt", "0.01", "4")}
    for fmt in formats:
        prob, dt, steps = probs[fmt]
        base = [real, "--format", fmt, "--problem", os.path.join(root, "data", "problems", prob), "--dt", dt,
                "--epsilon", "1e-9", "--steps", steps, "--sample", "2", "--quiet"]
        rc0, sw, e0 = run_env(base, env_sw)
        if rc0:
            verdict(False, "%s: the software backend's unforced run failed (rc %d): %s" % (fmt, rc0, short(e0)))
            continue
        swr = route_lines(sw)
        for name in sorted(door):
            var = door[name][0]
            if name in unpublished:
                rc, out, err = run_env(base, dict(env_dev, **{var: "on"}))
                lines = [l for l in err.splitlines() if l.startswith(DEVICE_LINE)]
                want = "%s=on asks for the %s route, and %s does not publish %s" % (var, name, art, door[name][1])
                rcs, outs, errs = run_env(base, dict(env_sw, **{var: "on"}))
                verdict(rc != 0 and want in err and lines == [probe_line] and rcs == 0,
                        "%s: %s=on on %s is REFUSED by name (rc %d): \"%s ...\"; the software backend runs the same "
                        "command (rc %d), so the refusal is the device's" % (fmt, var, art, rc, want, rcs),
                        [] if want in err else ["its stderr: %s" % short(err, 300)])
                continue
            rc, out, err = run_env(base, dict(env_dev, **{var: "on"}))
            lines = [l for l in err.splitlines() if l.startswith(DEVICE_LINE)]
            a, _ = norm_stdout(out, unpublished, "dev")
            b, _ = norm_stdout(sw, unpublished, "sw")
            dr, srr = route_lines(out).get(name), swr.get(name)
            # the route's own line is held below; staged_bytes is the door's
            # tally of what a route presents, which a route moves (entry 45)
            a = [re.sub(r" staged_bytes=\d+", "", l) for l in a if not l.startswith("# route %s=" % name)]
            b = [re.sub(r" staged_bytes=\d+", "", l) for l in b if not l.startswith("# route %s=" % name)]
            fd = first_diff(a, b)
            moved = dr and srr and dr[0] == "on" and dr[1] == "forced" and dr[3] == 0 and dr[2] > 0 \
                and srr[0] == "off" and srr[2] == 0 and dr[2] == srr[3]
            ph = physics(out) or {}
            verdict(rc == 0 and not fd and moved and lines == [probe_line],
                    "%s: %s=on on %s runs in the software backend's bits: %d lines identical to the unforced software "
                    "run but the route's own; the route on by=forced carried %s operation(s) the software run carried "
                    "as %s fallback(s), no fallback; %s steps, %s rejected, %s passes a step"
                    % (fmt, var, art, len(a), dr[2] if dr else "?", srr[3] if srr else "?", ph.get("steps"),
                       ph.get("rejected"), ph.get("mean_passes")),
                    ([] if not fd else ["line %d: device %r, software %r" % (fd[0] + 1, short(fd[1], 160),
                                                                             short(fd[2], 160))])
                    + ([] if rc == 0 else ["rc %d: %s" % (rc, short(err))]))
        # the door's self-test, against the software backend with the same routes hidden
        hide = ",".join(sorted(unpublished))
        rc, out, err = run_env([real, "--door-selftest", "--format", fmt], env_dev)
        rcs, outs, errs = run_env([real, "--door-selftest", "--format", fmt],
                                  dict(env_sw, CFT_REBOUND_HIDE_FEATURES=hide) if hide else env_sw)
        summ = re.compile(r"^door: (\d+) identical, (\d+) differ, (\d+) not compared, (\d+) accounting errors; "
                          r"controls (\d+) of (\d+)$", re.M)
        m, ms = summ.search(out), summ.search(outs)
        lines = [l for l in err.splitlines() if l.startswith(DEVICE_LINE)]
        # every line but the first, which names the backend, and the summary,
        # compared above; a route the device lacks reads absent, and hidden
        # on the software backend
        a = [re.sub(r"\babsent\b", "~", l) for l in out.splitlines() if not l.startswith("door: ")]
        b = [re.sub(r"\bhidden\b", "~", l) for l in outs.splitlines() if not l.startswith("door: ")]
        fd = first_diff(a, b)
        verdict(rc == 0 and m and ms and m.groups() == ms.groups() and m.group(2) == "0" and m.group(4) == "0"
                and m.group(5) == m.group(6) and not fd and lines == [probe_line],
                "%s: --door-selftest on %s: \"%s\", the software backend's with %s hidden: \"%s\"; %d case lines the "
                "same, absent for hidden" % (fmt, art, m.group(0)[6:] if m else "no summary (rc %d)" % rc,
                                             hide or "nothing", ms.group(0)[6:] if ms else "no summary", len(a)),
                [] if not fd else ["line %d: device %r, software %r" % (fd[0] + 1, short(fd[1], 160), short(fd[2], 160))])


# ---------------------------------------------------------------- 3. the seam

class Cols:
    """The record's own `# columns:` line, read here, apart from
    compare_run.py: system s sample k step, then the head columns, then six
    a body, then the exact-time pair after a `|`."""

    def __init__(self, line, n):
        words = line.split(":", 1)[1].split()
        step = words.index("step")
        self.head = words[step + 1:words.index("then")]
        self.body = words[words.index("body") + 1:words.index(";")]
        self.pair = words[words.index("pair") + 1:]
        # in a sample line the word sample is followed by its index, which
        # the columns line does not name: one token more than its words
        self.lead, self.n = step + 2, n

    def token(self, col, body=None):
        if body is None and col in self.head:
            return self.lead + self.head.index(col)
        if body is None and col in self.pair:
            return self.lead + len(self.head) + self.n * len(self.body) + 1 + self.pair.index(col)
        return self.lead + len(self.head) + body * len(self.body) + self.body.index(col)


def seam(shadow, work, art, probe_line, unpublished, env_dev, env_sw, real):
    print("3. the GPU seam: tools/compare_run.py holds a %s record of the ring against this tree's CPU reruns" % art)
    waves = os.path.join(work, "ring2w_z.txt")
    rc, out, err = run_env([sys.executable, "-B", os.path.join(here, "make_ring.py"), waves, "--seed", "1",
                            "--particles", "114"], env_sw)
    if rc:
        verdict(False, "make_ring.py failed (rc %d): %s" % (rc, short(err or out)))
        return
    zl = [l for l in open(waves) if l.split()[:1] == ["zonal"]]
    run_args = ["--format", "fp64", "--problem", waves, "--dt", "1e5", "--epsilon", "1e-9", "--max-iter", "60",
                "--steps", "8", "--sample", "4", "--quiet"]
    recpath = os.path.join(work, "ring2w_z.%s.rec" % art.replace(":", ""))
    rc, secs, d = run_witnessed(shadow, "seam", [os.path.join(shadow.dir, "ias15_cft" + EXE)] + run_args, env_dev,
                                recpath, recpath + ".err")
    rec_text = read(recpath)
    res = analyse(d, probe_line, unpublished)
    ok, text, detail = witness_verdict("the ring's record", res)
    E = len(re.findall(r"^# system \d+ steps_done=", rec_text, re.M))
    rejected = [int(x) for x in re.findall(r"^# system \d+ steps_done=\S+ \S+ steps_rejected=(\d+)", rec_text, re.M)]
    verdict(ok and rc == 0 and len(zl) == 1 and "# zonal " in rec_text and E == 2 and all(rejected),
            "the record: make_ring.py's %d waves of 64 bodies with the zonal line, 8 steps from a first step of 1e5 s "
            "on %s (rejections %s, so both waves are sampled); %s" % (E, art, rejected, text), detail)
    if not (rc == 0 and E == 2):
        return
    print("       device: %s" % (sorted(res["lines"])[0] if res["lines"] else "NONE"))
    drv = [sys.executable, "-B", os.path.join(here, "compare_run.py"), None, "--problem", waves, "--format", "fp64",
           "--arith", "rebound", "--exe", real, "--strata", "2", "--seed", "1"]

    def drive(path, tag):
        cmd = list(drv)
        cmd[3] = path
        rc, out, err = run_env(cmd, env_sw)
        with open(os.path.join(work, "compare_run.%s.log" % tag), "w") as f:
            f.write(out + err)
        mem = {int(m.group(1)): m.group(2) for m in re.finditer(r"^member (\d+) \(\S+\): (IDENTICAL|FAIL): ", out, re.M)}
        parts = re.findall(r"^  member (\d+): particle .*?held to the record: (IDENTICAL|FAIL)", out, re.M)
        return rc, out, err, mem, parts

    rc, out, err, mem, parts = drive(recpath, "clean")
    verdict(rc == 0 and mem == {0: "IDENTICAL", 1: "IDENTICAL"} and parts and all(p[1] == "IDENTICAL" for p in parts)
            and "\nRESULT: PASS: " in out and "backend=%s" % art in out,
            "compare_run.py, no artifact, on the %s record: %s; %d particle(s) replayed alone with the zonal line, %s "
            "(rc %d)" % (art, ", ".join("member %d %s" % kv for kv in sorted(mem.items())) or "no member", len(parts),
                         "all IDENTICAL" if parts and all(p[1] == "IDENTICAL" for p in parts) else "NOT all identical",
                         rc),
            [] if rc == 0 else [short(out + err, 400)])
    rows = rec_text.splitlines()
    n = int(re.search(r"\bN=(\d+)", rec_text).group(1))
    cols = Cols(next(l for l in rows if l.startswith("# columns:")), n)
    lines0 = [k for k, l in enumerate(rows) if l.startswith("system 0 sample ")][1:]   # sample 0 is the problem's
    # each label is the column as compare_run.py names it in its verdict
    plants = [("x", 1, lines0[-1], "x of body 1"), ("dt_next", None, lines0[0], "dt_next"),
              ("t_lo", None, lines0[0], "the exact-time pair's t_lo"), ("vz", n - 1, lines0[-1], "vz of body %d" % (n - 1)),
              ("t", None, lines0[-1], "t")]
    for col, body, i, label in plants:
        toks = rows[i].split()
        t = cols.token(col, body)
        old = toks[t]
        mant = old.split("p")[0]
        toks[t] = mant[:-1] + "%x" % (int(mant[-1], 16) ^ 1) + old[len(mant):]
        slug = re.sub(r"\W+", "_", label)
        bad = os.path.join(work, "ring2w_z.planted-%s.rec" % slug)
        with open(bad, "w", newline="\n") as f:
            f.write("\n".join(rows[:i] + [" ".join(toks)] + rows[i + 1:]) + "\n")
        rc, out, err, mem, parts = drive(bad, "planted-" + slug)
        named = re.search(r"^member 0 \(\S+\): FAIL: first differing line: record line %d \([^)]*\): %s(?: \([^)]*\))?: "
                          r"record " % (i + 1, re.escape(label)), out, re.M)
        verdict(rc == 1 and mem.get(0) == "FAIL" and mem.get(1) == "IDENTICAL" and named is not None,
                "the plant: one hex digit of %s on wave 0's record line %d (%s to %s): wave 0 FAILs naming that line "
                "and column, and wave 1 stays IDENTICAL (rc %d; particles %s)"
                % (label, i + 1, old, toks[t], rc, ", ".join("%s %s" % p for p in parts) or "none"),
                [] if named else [short(next((l for l in out.splitlines() if l.startswith("member 0 ")), out), 400)])


# ---------------------------------------------------------------- 4. the legs

def makefile_legs():
    """{(tier, leg): command} from the Makefile's quick tiers."""
    with open(os.path.join(root, "Makefile")) as f:
        text = f.read()
    legs = {}
    for tier in ("check-quick", "check-light"):
        m = re.search(r"^%s:[^\n]*\n((?:\t[^\n]*\n)+)" % re.escape(tier), text, re.M)
        if not m:
            raise Refused("the Makefile has no %s recipe to read the legs from" % tier)
        for line in m.group(1).splitlines():
            g = re.match(r"\t\$\(GATE_CACHE\) --leg (\S+) --mode \S+ .*? -- (.*)$", line)
            if g:
                legs[(tier, g.group(1))] = g.group(2)
    return legs


def expand(cmd, shadow_dir):
    words = []
    for w in cmd.split():
        w = w.replace("$(PYTHON)", sys.executable).replace("$(B)", shadow_dir).replace("$(EXE)", EXE)
        if "$(" in w:
            raise Refused("the leg's command has a variable this driver does not expand: %s" % cmd)
        words.append(w)
    return words


def leg_result(text):
    oks = len(re.findall(r"^\s*ok\b", text, re.M))
    fails = [l.strip() for l in text.splitlines() if re.match(r"^\s*FAIL\b", l)]
    m = re.findall(r"^RESULT: (PASS|FAIL)\b.*$", text, re.M)
    if m:
        return m[-1] == "PASS", oks, fails, re.findall(r"^RESULT: .*$", text, re.M)[-1]
    if "check_dropin: every case passed" in text:
        return True, oks, fails, "check_dropin: every case passed"
    g = re.search(r"^check_dropin: \d+ FAILURES$", text, re.M)
    return False, oks, fails, g.group(0) if g else "no RESULT line"


def main():
    ap = argparse.ArgumentParser(description="this tree on a GPU, held to the software backend gate by gate")
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--artifact", default="hip:0")
    ap.add_argument("--hip", default="", help="the Makefile's HIP, so a stale libcft.a can be named")
    ap.add_argument("--legs", default="", help="only these legs, by name, comma separated")
    ap.add_argument("--skip", default="", help="every leg but these")
    ap.add_argument("--no-legs", action="store_true", help="the device, the controls, the routes and the seam only")
    ap.add_argument("--parts", default="controls,routes,seam",
                    help="which of the controls, the routes and the seam run before the legs (all three by default)")
    ap.add_argument("--formats", default="fp64,fp128,fp256", help="the routes' formats")
    ap.add_argument("--keep", action="store_true")
    ap.add_argument("--list", action="store_true", help="print each leg's command and exit")
    args = ap.parse_args()
    build = os.path.abspath(args.build)
    art = args.artifact
    real = os.path.join(build, "ias15_cft" + EXE)
    print("check_gpu: this tree on %s, held to the software backend gate by gate (ROUND5.md, P1)" % art)
    try:
        legs = makefile_legs()
        names = [leg for _, leg in GPU_LEGS]
        for tier, leg in GPU_LEGS:
            if (tier, leg) not in legs:
                raise Refused("the Makefile's %s has no leg named %s" % (tier, leg))
        want = [x for x in args.legs.split(",") if x]
        skip = [x for x in args.skip.split(",") if x]
        unknown = [x for x in want + skip if x not in names]
        if unknown:
            raise Refused("no such leg: %s; the legs are %s" % (", ".join(unknown), ", ".join(names)))
        chosen = [(t, l) for t, l in GPU_LEGS if (not want or l in want) and l not in skip and not args.no_legs]
        if args.list:
            for tier, leg in chosen:
                print("  %s (%s): %s" % (leg, tier, " ".join(expand(legs[(tier, leg)], "<shadow>"))))
            return 0
        print("  host %s, %s %s; Python %s" % (socket.gethostname(), platform.system(), platform.release(),
                                               platform.python_version()))
        for prog in ("ias15_cft", "check_dropin", "ias15_ref"):
            print("  %s, built %s" % (os.path.join(build, prog + EXE), built(os.path.join(build, prog + EXE))))
        lib = os.path.join(root, "third_party", "cft-fp256", "host", "libcft.a")
        has_hip = libcft_has_hip(lib)
        print("  %s, built %s; backend_hip.o in it: %s" % (lib, built(lib), {True: "yes", False: "NO",
                                                                             None: "unknown (no ar)"}[has_hip]))
        if not os.path.isfile(real):
            raise Refused("%s is not built: make HIP=1 all programs" % real)

        # 0. the device
        print("0. the device")
        rc, out, err, caps, rts, lines = probe(real, art)
        if rc != 0:
            why = short(err, 600)
            if args.hip == "1" and has_hip is False:
                why += (" -- and this make asked for HIP=1 while libcft.a has no backend_hip.o: the stale-library "
                        "trap (V0, entry 54). make clean, then make HIP=1 all programs")
            raise Refused("this host has no %s for this build: %s" % (art, why))
        if art.startswith("hip:") and caps.get("backend") != "hip":
            raise Refused("%s opened, and its backend is %r, not hip" % (art, caps.get("backend")))
        if len(lines) != 1 or not lines[0].startswith("%s%s = " % (DEVICE_LINE, art)):
            raise Refused("%s opened with CFT_HIP_INFO=1 and printed %d device line(s), not one naming it: %s"
                          % (art, len(lines), short(err, 400)))
        probe_line = lines[0]
        rcs, outs, errs, scaps, srts, slines = probe(real, "")
        unpublished = sorted(r for r, (bit, pub) in rts.items() if not pub)
        published = sorted(r for r, (bit, pub) in rts.items() if pub)
        print("  device: %s" % probe_line)
        try:
            diff = int(scaps.get("seq_features", "0"), 16) & ~int(caps.get("seq_features", "0"), 16)
        except ValueError:
            diff = -1
        verdict(rcs == 0 and not slines and caps.get("formats") and caps.get("max_scratch") and caps.get("max_consts"),
                "%s --probe: backend %s, formats %s, seq_features %s (the software backend's %s less 0x%x, %s), "
                "max_scratch %s, max_consts %s, buffers_resident %s; routes published: %s; NOT published: %s"
                % (art, caps.get("backend"), caps.get("formats"), caps.get("seq_features"),
                   scaps.get("seq_features"), diff, feature_names(diff) if diff > 0 else "none",
                   caps.get("max_scratch"), caps.get("max_consts"), caps.get("buffers_resident"),
                   ", ".join(published) or "none", ", ".join(unpublished) or "none"))
        if os.name == "nt":
            raise Refused("the witness is a POSIX script; run this on the GPU host's Linux")
    except Refused as e:
        print("REFUSED: %s" % e)
        print("RESULT: REFUSED")
        return 2

    work = os.path.join(build, "check_gpu", str(os.getpid()))
    shutil.rmtree(work, ignore_errors=True)
    os.makedirs(work)
    shadow = Shadow(build, work)
    base_env = {k: v for k, v in os.environ.items()
                if not k.upper().startswith("CFT_REBOUND_ROUTE_") and not k.upper().startswith("CFT_REBOUND_HIDE")}
    dropped = sorted(set(os.environ) - set(base_env))
    if dropped:
        print("  note: %s left out of every leg's environment, so it cannot choose a route behind the gate's back"
              % ", ".join(dropped))
    env_dev = dict(base_env, CFT_REBOUND_ARTIFACT=art, CFT_HIP_INFO="1")
    env_sw = {k: v for k, v in base_env.items() if k not in ("CFT_REBOUND_ARTIFACT", "CFT_HIP_INFO")}
    print("  work: %s; every witnessed program's kept output is under %s" % (work, shadow.legs))

    t_start = time.time()
    parts = [p for p in args.parts.split(",") if p]
    if [p for p in parts if p not in ("controls", "routes", "seam")]:
        print("REFUSED: --parts takes controls, routes and seam")
        return 2
    if not want or args.no_legs:
        if "controls" in parts:
            controls(shadow, work, art, probe_line, unpublished, env_dev)
        if "routes" in parts:
            routes(real, art, probe_line, unpublished, published, (env_dev, env_sw),
                   [f for f in args.formats.split(",") if f])
        if "seam" in parts:
            seam(shadow, work, art, probe_line, unpublished, env_dev, env_sw, real)
        skipped = [p for p in ("controls", "routes", "seam") if p not in parts]
        if skipped:
            print("  note: NOT RUN, by --parts: %s" % ", ".join(skipped))
    else:
        print("  note: the controls, the routes and the seam NOT RUN: --legs names a subset")

    if chosen:
        print("4. the legs, each in its quick form, as the Makefile runs it, on %s through the witness" % art)
    summary = []
    for tier, leg in chosen:
        cmd = expand(legs[(tier, leg)], shadow.dir)
        log = os.path.join(work, leg + ".log")
        print("  ... %s: %s (log %s), from %s" % (leg, " ".join(cmd[1:] if cmd[0] == sys.executable else cmd), log,
                                                  time.strftime("%H:%M:%S")))
        sys.stdout.flush()
        rc, secs, d = run_witnessed(shadow, leg, cmd, env_dev, log)
        text = read(log) or ""
        passed, oks, fails, result = leg_result(text)
        res = analyse(d, probe_line, unpublished)
        wok, wtext, wdetail = witness_verdict("witness", res)
        ok = rc == 0 and passed and not fails and wok
        verdict(ok, "%s (%s form): rc %d, %s, %d ok, %d FAIL, %dm%02ds; %s"
                % (leg, "--light" if tier == "check-light" else "quick", rc, result, oks, len(fails),
                   int(secs) // 60, int(secs) % 60, wtext),
                (["leg FAIL: %s" % short(f, 300) for f in fails[:5]] + wdetail
                 + ["physics: %s" % physics_text(res)]
                 + ["device: %s" % l for l in sorted(res["lines"])]))
        summary.append((leg, ok, secs))

    print()
    print("check_gpu: %d ok, %d FAIL, in %dm%02ds on %s" % (OKS[0], len(FAILS), int(time.time() - t_start) // 60,
                                                          int(time.time() - t_start) % 60, art))
    print("  device: %s" % probe_line)
    for leg, ok, secs in summary:
        print("  %-22s %s  %dm%02ds" % (leg, "PASS" if ok else "FAIL", int(secs) // 60, int(secs) % 60))
    if FAILS:
        print("  kept: %s" % work)
        for f in FAILS:
            print("  FAILED: %s" % short(f, 200))
    elif args.keep:
        print("  kept: %s" % work)
    else:
        shutil.rmtree(work, ignore_errors=True)
    print("RESULT: %s" % ("FAIL" if FAILS else "PASS"))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
