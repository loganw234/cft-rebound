#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# What runs where when this port runs on a GPU, measured rather than read
# (round 5, P1; ROUND5.md, "What hip:N does with this port's calls"). Run on
# a GPU host, in a tree built with HIP=1:
#
#   python3 hw/gpu/where.py [--build build] [--artifact hip:0]
#                           [--make-args "HIP=1"] [--formats fp64,fp128,fp256]
#
# A. THE CLOCK'S AUGMENTED SUM, on the host at every step. ias15_cft runs
#    under gdb with a breakpoint on cft_augmented_add and on each of the GPU
#    backend's four entry points (cfth_run, cfth_reduce, cfth_reduce_seg,
#    cfth_program_run: every operation libcft sends to the device goes
#    through one). Each backend call is attributed by its stack to the libcft
#    call it serves. Measured: how many times cft_augmented_add was entered,
#    against the record's own step attempts (steps_done + steps_rejected);
#    and how many backend calls were made from inside it.
#
# B. THE DIVIDE AND SQUARE ROOT, by libcft's program route. The same runs
#    attribute every backend call made inside cft_div and cft_sqrt, and the
#    host's cft_sf_round_pack calls inside them, the per-element finish the
#    program route does on the host. Then CFT_DIVSQRT_SEQ forces each route:
#    =1 the program route, =0 the chunk route. The default must attribute
#    as =1 does, and =0 must not; and the records of all three must be the
#    same, with the time a divide takes on each, from the trailer.
#
# C. THE INDEX-TABLE ROUTES, by=absent. Each default record's `# route`
#    lines: indexed and pair_gather off by=absent, no call, their operations
#    counted as fallbacks. (Forced on, check_gpu.py holds their refusal.)
#
# D. THE PROGRAMS, within the device's scratch slots and bank. hw/gpu/
#    program_caps.c, built here with ias15_cft's own link line, loads every
#    image the port runs - programs/out/*.cftp, gravity's among them, and the
#    corrector's normalabs images from the pinned tree - on the device, and
#    prints what cft_program_get_info() says each uses against what the
#    device publishes. The boundary controls rewrite one header word: a
#    scratch block of the device's own depth must load, and one slot past it
#    must be refused by name; likewise the bank's 512 constants and 513.
#
# Everything is printed as measured, with the commands; nothing is asserted
# beyond what each part says it holds, and the summary says which held.

import argparse
import glob
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
import time

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(os.path.dirname(here))
EXE = ".exe" if os.name == "nt" else ""

GDB_SCRIPT = r'''
import gdb
TARGETS = ("cft_augmented_add", "cft_div", "cft_sqrt")
BACKEND = ("cfth_run", "cfth_reduce", "cfth_reduce_seg", "cfth_program_run")
HOST = ("cft_sf_round_pack",)
count = {}
inside = {}


class Counter(gdb.Breakpoint):
    def __init__(self, name):
        super().__init__(name, internal=True)
        self.fn = name

    def stop(self):
        count[self.fn] = count.get(self.fn, 0) + 1
        if self.fn in TARGETS:
            return False
        f = gdb.newest_frame().older()
        depth = 0
        while f is not None and depth < 40:
            n = f.name()
            if n in TARGETS:
                inside[(n, self.fn)] = inside.get((n, self.fn), 0) + 1
                break
            f = f.older()
            depth += 1
        return False


gdb.execute("set pagination off")
gdb.execute("set confirm off")
for sig in ("SIGUSR1", "SIGUSR2", "SIGPIPE", "SIG32", "SIG33", "SIG34", "SIG35"):
    try:
        gdb.execute("handle %s nostop noprint pass" % sig)
    except gdb.error:
        pass
for n in TARGETS + BACKEND + HOST:
    Counter(n)
gdb.execute("run")
for k in sorted(count):
    print("WHERE count %s %d" % (k, count[k]))
for (t, b) in sorted(inside):
    print("WHERE inside %s %s %d" % (t, b, inside[(t, b)]))
'''


def run(cmd, env=None, cwd=None):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=env, cwd=cwd or root,
                       stdin=subprocess.DEVNULL)
    return p.returncode, p.stdout, p.stderr


def trailer(text):
    m = re.search(r"^# steps_done=.*$", text, re.M)
    return dict(t.partition("=")[::2] for t in m.group(0)[2:].split()) if m else {}


def routes_of(text):
    return {m.group(1): m.group(0)[2:] for m in re.finditer(r"^# route (\w+)=.*$", text, re.M)}


def records_same(a, b):
    def norm(t):
        out = []
        for l in t.splitlines():
            if l.startswith("# steps_done="):
                l = re.sub(r" (?:seconds|steps_per_s|wall_seconds|t_[a-z_]+)=\S+", "", l)
            out.append(l)
        return out
    return norm(a) == norm(b)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--artifact", default="hip:0")
    ap.add_argument("--make-args", default="HIP=1", help="the make variables the tree was built with")
    ap.add_argument("--formats", default="fp64,fp128,fp256")
    ap.add_argument("--work", default=None)
    ap.add_argument("--parts", default="A,B,C,D")
    args = ap.parse_args()
    build = os.path.abspath(args.build)
    real = os.path.join(build, "ias15_cft" + EXE)
    art = args.artifact
    work = args.work or tempfile.mkdtemp(prefix="where-")
    os.makedirs(work, exist_ok=True)
    parts = args.parts.split(",")
    formats = args.formats.split(",")
    held = []
    env_dev = {k: v for k, v in os.environ.items() if not k.upper().startswith("CFT_")}
    env_dev.update({"CFT_REBOUND_ARTIFACT": art, "CFT_HIP_INFO": "1"})
    env_sw = {k: v for k, v in os.environ.items() if not k.upper().startswith("CFT_")}
    print("where: what runs where on %s, measured; %s, built %s; work %s"
          % (art, real, time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(real))), work))
    kepler = os.path.join(root, "data", "problems", "kepler.txt")
    rc, o, e = run([real, "--probe"], env_dev)
    dev = [l for l in e.splitlines() if l.startswith("cft hip: ")]
    print("  device: %s" % (dev[0] if dev else "NONE (rc %d): %s" % (rc, e.strip()[-300:])))
    if rc or not dev:
        print("RESULT: REFUSED: %s does not open" % art)
        return 2

    have_gdb = shutil.which("gdb") is not None
    if not have_gdb:
        print("  note: no gdb on this host, so A and B's attribution are NOT MEASURED here; B's routes are held by "
              "CFT_DIVSQRT_SEQ, their records and their times alone")
    if "A" in parts or "B" in parts or "C" in parts:
        gdbpy = os.path.join(work, "where_gdb.py")
        with open(gdbpy, "w") as f:
            f.write(GDB_SCRIPT)
        # A rejected first attempt and accepted ones: dt0 far too large.
        steps = {"fp64": "4", "fp128": "3", "fp256": "2"}
        for fmt in formats:
            base = ["--format", fmt, "--problem", kepler, "--dt", "1.0", "--epsilon", "1e-9", "--steps", steps[fmt],
                    "--sample", steps[fmt], "--quiet"]
            res = {}
            for label, extra in (("default", {}), ("CFT_DIVSQRT_SEQ=1", {"CFT_DIVSQRT_SEQ": "1"}),
                                 ("CFT_DIVSQRT_SEQ=0", {"CFT_DIVSQRT_SEQ": "0"})):
                env = dict(env_dev, **extra)
                rc, out, err = run([real] + base, env)              # timed and recorded, no debugger
                if have_gdb:
                    g_rc, g_out, g_err = run(["gdb", "-q", "-batch", "-x", gdbpy, "--args", real] + base, env)
                else:
                    g_rc, g_out, g_err = -1, "", "no gdb"
                count = {m.group(1): int(m.group(2)) for m in re.finditer(r"^WHERE count (\S+) (\d+)$", g_out, re.M)}
                inside = {(m.group(1), m.group(2)): int(m.group(3))
                          for m in re.finditer(r"^WHERE inside (\S+) (\S+) (\d+)$", g_out, re.M)}
                rec = "\n".join(l for l in g_out.splitlines() if not l.startswith("WHERE ") and
                                (l.startswith("#") or l.startswith("sample ")))
                res[label] = (rc, out, err, count, inside, rec, g_rc, g_err)
                t = trailer(out)
                print("  %s %s: rc %d; %s steps, %s rejected, %s calls, %s divsqrt_calls, wall_seconds %s, t_divsqrt %s"
                      % (fmt, label, rc, t.get("steps_done"), t.get("steps_rejected"), t.get("calls"),
                         t.get("divsqrt_calls"), t.get("wall_seconds"), t.get("t_divsqrt")))
                if have_gdb:
                    print("    under gdb (rc %d): entered %s" % (g_rc, ", ".join("%s %d" % kv for kv in sorted(count.items()))
                                                                 or "NOTHING: %s" % g_err.strip()[-300:]))
                for (tg, b), n in sorted(inside.items()):
                    print("    inside %s: %s %d" % (tg, b, n))
                if have_gdb and not records_same(rec, out):
                    print("    note: the record under gdb is not the record without it")
            rc, out, err, count, inside, rec, g_rc, g_err = res["default"]
            t = trailer(out)
            attempts = int(t.get("steps_done", 0)) + int(t.get("steps_rejected", 0))
            if "A" in parts and not have_gdb:
                print("  A %s: NOT MEASURED, no gdb on this host" % fmt)
            elif "A" in parts:
                n_aug = count.get("cft_augmented_add", 0)
                from_aug = sum(n for (tg, b), n in inside.items() if tg == "cft_augmented_add" and b.startswith("cfth_"))
                total = sum(count.get(b, 0) for b in ("cfth_run", "cfth_reduce", "cfth_reduce_seg", "cfth_program_run"))
                ok = n_aug == attempts and from_aug == 0 and total > 0 and int(t.get("steps_rejected", 0)) > 0
                held.append(("A %s" % fmt, ok))
                print("  A %s: cft_augmented_add entered %d times in %d step attempts (%s accepted, %s rejected); "
                      "%d of the run's %d backend calls were made inside it: %s"
                      % (fmt, n_aug, attempts, t.get("steps_done"), t.get("steps_rejected"), from_aug, total,
                         "HELD, on the host at every attempt" if ok else "NOT HELD"))
            if "B" in parts:
                def split(lab):
                    ins = res[lab][4]
                    return {b: sum(n for (tg, bb), n in ins.items() if tg in ("cft_div", "cft_sqrt") and bb == b)
                            for b in ("cfth_run", "cfth_program_run", "cfth_reduce", "cft_sf_round_pack")}
                d, s1, s0 = split("default"), split("CFT_DIVSQRT_SEQ=1"), split("CFT_DIVSQRT_SEQ=0")
                n_ds = res["default"][3].get("cft_div", 0) + res["default"][3].get("cft_sqrt", 0)
                same = all(records_same(res[l][1], out) for l in res) and all(res[l][0] == 0 for l in res)

                def per(lab):
                    tt = trailer(res[lab][1])
                    return 1e6 * float(tt.get("t_divsqrt", "0")) / max(1, int(tt.get("divsqrt_calls", "0")))
                if have_gdb:
                    ok = (d == s1 and d["cfth_program_run"] > 0 and d["cfth_run"] == 0 and d["cft_sf_round_pack"] > 0
                          and s0["cfth_program_run"] == 0 and s0["cfth_run"] > 0 and same
                          and n_ds == int(t.get("divsqrt_calls", -1)))
                else:
                    # by time alone: the default within 25% of the program route, and the chunk route at least
                    # twice either; one run each, so this holds a difference of several times, never a few percent
                    pd_, p1_, p0_ = per("default"), per("CFT_DIVSQRT_SEQ=1"), per("CFT_DIVSQRT_SEQ=0")
                    ok = same and p1_ > 0 and abs(pd_ - p1_) <= 0.25 * p1_ and p0_ >= 2 * max(pd_, p1_)
                held.append(("B %s" % fmt, ok))
                if not have_gdb:
                    print("  B %s: the three records %s; a divide or root %.0f us by default, %.0f us forced to the "
                          "program route, %.0f us forced to the chunk route (t_divsqrt over divsqrt_calls, one run "
                          "each; no gdb here, so by the forcing variable and the time alone): %s"
                          % (fmt, "IDENTICAL" if same else "DIFFER", per("default"), per("CFT_DIVSQRT_SEQ=1"),
                             per("CFT_DIVSQRT_SEQ=0"), "HELD, the program route" if ok else "NOT HELD"))
                else:
                    print("  B %s: cft_div and cft_sqrt entered %d times (the trailer's divsqrt_calls %s). Inside them, "
                          "by default: %d program runs and %d elementwise runs on the device, %d host round_packs; "
                          "CFT_DIVSQRT_SEQ=1: %d, %d, %d; CFT_DIVSQRT_SEQ=0: %d, %d, %d. The three records %s. "
                          "A divide or root: %.0f us by default, %.0f us forced to the program route, %.0f us forced to "
                          "the chunk route (t_divsqrt over divsqrt_calls, one run each): %s"
                          % (fmt, n_ds, t.get("divsqrt_calls"), d["cfth_program_run"], d["cfth_run"], d["cft_sf_round_pack"],
                             s1["cfth_program_run"], s1["cfth_run"], s1["cft_sf_round_pack"], s0["cfth_program_run"],
                             s0["cfth_run"], s0["cft_sf_round_pack"], "IDENTICAL" if same else "DIFFER",
                             per("default"), per("CFT_DIVSQRT_SEQ=1"), per("CFT_DIVSQRT_SEQ=0"),
                             "HELD, the program route" if ok else "NOT HELD"))
            if "C" in parts:
                r = routes_of(out)
                absent = sorted(k for k, v in r.items() if " by=absent " in v + " ")
                ok = bool(absent) and all(re.search(r"=off by=absent calls=0 fallbacks=\d+", r[k]) for k in absent)
                held.append(("C %s" % fmt, ok))
                print("  C %s: the record's route lines: %s: %s" % (fmt, "; ".join(r[k] for k in sorted(r)),
                                                                  "HELD, %s absent with no call" % " and ".join(absent)
                                                                  if ok else "NOT HELD"))

    if "D" in parts:
        print("  D: the programs against the device's sequencer capacities")
        rc, out, err = run(["make", "-n", "-W", "src/ias15_cft.c"] + shlex.split(args.make_args)
                           + [os.path.relpath(real, root)])
        link = next((l for l in out.splitlines() if " -o %s " % os.path.relpath(real, root) in l + " "), None)
        if not link:
            print("    NOT RUN: make -n printed no link line for %s: %s" % (real, (out + err)[-300:]))
            held.append(("D", False))
        else:
            words = shlex.split(link)
            words = [w for w in words if w != "src/ias15_cft.c"]
            o = words.index("-o")
            probe = os.path.join(work, "program_caps" + EXE)
            words[o + 1] = probe
            words.insert(o + 2, os.path.join(here, "program_caps.c"))
            rc, out, err = run(words)
            print("    built with ias15_cft's link line (rc %d): %s" % (rc, " ".join(words)))
            asm = os.path.join(root, "third_party", "cft-fp256", "host", "cft-asm" + EXE)
            images = sorted(glob.glob(os.path.join(root, "programs", "out", "*.cftp")))
            # the corrector's images, at the formats the Makefile embeds them (NORMALABS_SRC)
            with open(os.path.join(root, "Makefile")) as f:
                m = re.search(r"^NORMALABS_SRC\s*:=\s*\$\(foreach f,([^,]+),", f.read(), re.M)
            fmts = m.group(1).split() if m else []
            if not fmts:
                print("    note: the Makefile's NORMALABS_SRC was not read, so no normalabs image is loaded")
            for src in [os.path.join(root, "third_party", "cft-fp256", "programs", "normalabs-%s.cfta" % f) for f in fmts]:
                dst = os.path.join(work, os.path.basename(src)[:-5] + ".cftp")
                r2, o2, e2 = run([asm, src, "-o", dst])
                if r2 == 0:
                    images.append(dst)
                else:
                    print("    %s did not assemble: %s" % (src, (o2 + e2)[-200:]))
            for label, env in (("device", env_dev), ("software", env_sw)):
                a = art if label == "device" else ""
                rc, out, err = run([probe, a] + images, env)
                info = [m.groupdict() for m in re.finditer(
                    r"^image (?P<name>\S+): loaded n_insns=(?P<insns>\d+) n_consts=(?P<consts>\d+) max_deposits=(?P<dep>\d+) "
                    r"flags=(?P<flags>\S+) n_scratch_in=(?P<sin>\d+) n_scratch_out=(?P<sout>\d+) "
                    r"scratch_used=(?P<used>\d+)$", out, re.M)]
                refused = re.findall(r"^image (\S+): refused: (.*)$", out, re.M)
                caps = re.search(r"^device \S+ backend=(\S+) max_scratch=(\d+) max_consts=(\d+) max_deposits=(\d+) "
                                 r"max_insns=(\d+) seq_features=(\S+)$", out, re.M)
                if not caps:
                    print("    %s: the probe did not open (rc %d): %s" % (label, rc, (out + err)[-300:]))
                    held.append(("D", False))
                    continue
                ms, mc = int(caps.group(2)), int(caps.group(3))
                top = lambda k: max(info, key=lambda i: int(i[k])) if info else {"name": "-", k: "0"}
                used, cons, sin, sout = top("used"), top("consts"), top("sin"), top("sout")
                ok = (len(info) == len(images) and not refused and int(used["used"]) <= ms and int(cons["consts"]) <= mc
                      and int(sin["sin"]) <= ms and int(sout["sout"]) <= ms)
                print("    %s (backend %s, max_scratch %d, max_consts %d, seq_features %s): %d of %d images loaded, %d "
                      "refused; the most scratch an image's instructions touch %s (%s), its block preloads %s (%s) and "
                      "reads back %s (%s); the most constants one addresses %s (%s): %s"
                      % (label, caps.group(1), ms, mc, caps.group(6), len(info), len(images), len(refused),
                         used["used"], used["name"], sin["sin"], sin["name"], sout["sout"], sout["name"], cons["consts"],
                         cons["name"], "HELD, every program within the device's depth and bank" if ok else "NOT HELD"))
                for name, why in refused:
                    print("      refused %s: %s" % (name, why))
                if label == "device":
                    held.append(("D images", ok))
                    if not ok:
                        continue
                    scr_io = [i for i in info if int(i["flags"], 16) & 2]
                    bank_ext = [i for i in info if int(i["flags"], 16) & 1]
                    for flag, n, what in (("--scratch-in", ms, "a scratch block of the device's depth"),
                                          ("--scratch-in", ms + 1, "one slot past it"),
                                          ("--consts", mc, "a bank of the device's 512"),
                                          ("--consts", mc + 1, "one constant past it")):
                        rc, o2, e2 = run([probe, a, flag, str(n)] + images, env)
                        ld = re.findall(r"^image (\S+): loaded", o2, re.M)
                        rf = re.findall(r"^image (\S+): refused: (.*)$", o2, re.M)
                        want_load = n in (ms, mc)
                        okb = (ld and not rf) if want_load else (rf and not ld)
                        held.append(("D %s %d" % (flag, n), bool(okb)))
                        print("    control %s %d, %s, on the %d images whose flags allow it: %d loaded, %d refused%s: %s"
                              % (flag, n, what, len(scr_io) if flag == "--scratch-in" else len(bank_ext), len(ld),
                                 len(rf), (" - \"%s\"" % rf[0][1][:220]) if rf else "",
                                 "HELD" if okb else "NOT HELD"))

    print("where: %s" % ", ".join("%s %s" % (k, "held" if v else "NOT HELD") for k, v in held))
    print("RESULT: %s" % ("HELD" if held and all(v for _, v in held) else "NOT HELD"))
    return 0 if held and all(v for _, v in held) else 1


if __name__ == "__main__":
    sys.exit(main())
