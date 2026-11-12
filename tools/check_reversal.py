#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The reversal runs the problem it was given (VALIDATION entry 52).
#
# tools/reversal.py runs N fixed steps forward, writes the last sample as a
# problem through tools/state_to_problem.py, runs N steps back from it, and
# holds the return to the start. Until entry 52 the turn problem kept G and
# the bodies and dropped the problem's `active` and `zonal` lines (ROUND4.md,
# P0): a zonal problem went forward with J2 and J4 and came back without
# them, and missed by 4e-3 of its size where binary64 returns to 1e-16.
# No gate ran either tool; study DEBRIS-A found the gap by reading them.
#
# THE PROBLEM: data/problems/zonal_saturn.txt with `active 1` for its
# `active 2`, so that Mimas (GM 2.5) is a test particle with a mass. As
# committed, the file's one test particle is the massless ring particle,
# whose reaction is an exact zero when every body is active: its active
# line changes no bit of a return (measured for entry 52: the return's
# three body lines identical with and without it), so it could not be
# held here.
#
# 0. A RECORD WITHOUT THE LINES GIVES THE FILE IT ALWAYS GAVE. The committed
#    records of results/horizon/reversal/ (kepler and outer, three formats,
#    forward and back) rebuild their committed turn and return problems,
#    every line but the one that names the record's path.
# 1. THE BACKWARD LEG RUNS THE PROBLEM. At each format (binary64 only in
#    --quick), tools/reversal.py runs 100 steps of 100 s each way (rc 0), and:
#    - its turn and return problems carry `active 1` and the zonal line,
#      value for value as the problem gives them (read here, not through
#      the tools);
#    - the backward record's `# active` and `# zonal` lines are the forward
#      record's, word for word, and name what the file names: N - 1 = 2
#      bodies feel the term, and body 0 takes the reaction of active - 1 = 0;
#    - the return lies within B of the start, in max |dx| / L as computed
#      here, and reversal.py prints the same figure: B is 1e-14 at binary64
#      and 1e-32 at binary128, about 100 unit roundoffs, where the arithmetic
#      dominates; 1e-40 at binary256, where IAS15's own reversal error does
#      (it is not time-symmetric: tools/reversal.py's head). Measured for
#      entry 52: 9.3e-17, 8.0e-35 and 1.5e-43.
# 2. EACH LINE CARRIES WEIGHT, so the bound resolves it. At binary64 the turn
#    problem with its zonal line deleted, and with its active line deleted,
#    runs back as reversal.py runs its backward leg (its leg_command()), and
#    each return must miss by more than 1000 B. Measured: 4.0e-3 and 1.7e-7.
# 3. A DROPPED LINE FAILS: controls from copies of the tools, each mutation
#    shown to apply. A copy of state_to_problem.py that drops the zonal line,
#    and one that drops the active line: reversal.py refuses the turn
#    problem, naming the line (rc 1). With reversal.py's own refusal deleted
#    as well, the copy runs to its end (rc 0) and part 1's checks FAIL on its
#    run: the line missing from both problems, the records' lines apart, and
#    the return outside B.
#
#   python tools/check_reversal.py [--build build] [--quick] [--keep]

import argparse
import difflib
import os
import re
import shutil
import subprocess
import sys
import time

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import hex_to_fraction  # noqa: E402
import reversal  # noqa: E402  (leg_command() only: part 2 runs the backward leg as it does)

root = os.path.dirname(here)
PY = sys.executable
STEPS, DT = 100, "100"
BOUND = {"fp64": 1e-14, "fp128": 1e-32, "fp256": 1e-40}
TOOLS = ("reversal.py", "state_to_problem.py", "check_equivalence.py")

FAILS = []


def verdict(ok, text):
    print("  %s %s" % ("ok  " if ok else "FAIL", text))
    if not ok:
        FAILS.append(text)
    return ok


def run(cmd):
    p = subprocess.run(cmd, cwd=root, capture_output=True, text=True)
    return p.returncode, p.stdout, p.stderr


def words(path, word):
    """The lines of a problem file that start with `word`, split."""
    return [l.split() for l in open(path) if l.split()[:1] == [word]]


def exact(tokens):
    """A line's values: hex as its exact value and sign, the rest as text."""
    return [hex_to_fraction(t) if t.lstrip("+-")[:2].lower() == "0x" else t for t in tokens]


def bodies(path):
    return [(t[1], [hex_to_fraction(v)[0] for v in t[2:9]]) for t in words(path, "body")]


def miss(start, ret):
    """max |dx| / L of the return against the start, L the start's largest
    coordinate: computed here, exactly, and only then rounded."""
    b0, b1 = bodies(start), bodies(ret)
    if [n for n, _ in b0] != [n for n, _ in b1]:
        return None
    L = max(abs(q) for _, u in b0 for q in u[1:4])
    dx = max(abs(u[c] - v[c]) for (_, u), (_, v) in zip(b0, b1) for c in (1, 2, 3))
    return float(dx / L)


def header(rec):
    return [l.rstrip("\n") for l in open(rec) if l.startswith(("# active ", "# zonal "))]


def one(keep, suffix):
    hits = [f for f in os.listdir(keep) if f.endswith(suffix)]
    return os.path.join(keep, hits[0]) if len(hits) == 1 else None


def reverse(tools, problem, fmt, keep, build):
    """tools/reversal.py from the directory `tools` on the problem; its
    rc, its output and the four files it wrote (None where absent)."""
    shutil.rmtree(keep, ignore_errors=True)
    rc, out, err = run([PY, os.path.join(tools, "reversal.py"), "--problem", problem, "--format", fmt,
                        "--dt", DT, "--steps", str(STEPS), "--build", build, "--keep", keep])
    files = {s: one(keep, s) if os.path.isdir(keep) else None
             for s in ("_fwd.rec", "_turn.txt", "_back.rec", "_return.txt")}
    return rc, out + err, files


def checks(problem, fmt, out, files):
    """Part 1's checks on one run: a list of (ok, text), and the miss."""
    res = []
    for what in ("_turn.txt", "_return.txt"):
        path = files[what]
        for word in ("active", "zonal"):
            want = [exact(t[1:]) for t in words(problem, word)]
            got = [exact(t[1:]) for t in words(path, word)] if path else None
            res.append((got == want and len(want) == 1,
                        "the %s problem's %s line is the problem's, value for value: %s"
                        % (what[1:-4], word, " ".join(words(path, word)[0]) if got else "MISSING")))
    hf = header(files["_fwd.rec"]) if files["_fwd.rec"] else []
    hb = header(files["_back.rec"]) if files["_back.rec"] else []
    N = len(words(problem, "body"))
    n_act = int(words(problem, "active")[0][1])
    k = int(words(problem, "zonal")[0][1])
    named = (len(hf) == 2 and hf[0].startswith("# active %d of %d bodies a system:" % (n_act, N))
             and re.search(r"it acts on %d bodies a system, and body %d takes the reaction of %d$" % (N - 1, k, n_act - 1), hf[1]))
    res.append((hf == hb and bool(named),
                "the backward record's # active and # zonal lines are the forward record's, and name the file's "
                "counts (N - 1 = %d feel the term, body %d takes the reaction of active - 1 = %d): %s"
                % (N - 1, k, n_act - 1, "the same" if hf == hb else "%r against %r" % (hb, hf))))
    m = miss(problem, files["_return.txt"]) if files["_return.txt"] else None
    printed = re.search(r"max\|dx\|/L=(\S+)", out)
    agree = m is not None and printed and abs(float(printed.group(1)) - m) <= 1e-3 * m
    res.append((m is not None and m <= BOUND[fmt] and bool(agree),
                "the return lies %s of the start in max |dx| / L (the bound is %.0e), and reversal.py prints %s"
                % ("%.3e" % m if m is not None else "nowhere", BOUND[fmt], printed.group(1) if printed else "nothing")))
    return res, m


def mutate(src_dir, dst_dir, edits):
    """A copy of the tools with each (file, old, new) applied exactly once;
    the number of lines each changed, or SystemExit if one does not apply."""
    shutil.rmtree(dst_dir, ignore_errors=True)
    os.makedirs(dst_dir)
    for f in TOOLS:
        shutil.copyfile(os.path.join(src_dir, f), os.path.join(dst_dir, f))
    changed = []
    for f, old, new in edits:
        p = os.path.join(dst_dir, f)
        text = open(p).read()
        if text.count(old) != 1:
            raise SystemExit("check_reversal: the mutation of %s does not apply: %r occurs %d times"
                             % (f, old, text.count(old)))
        after = text.replace(old, new)
        open(p, "w", newline="\n").write(after)
        d = [l for l in difflib.unified_diff(text.splitlines(), after.splitlines(), lineterm="", n=0)
             if l[:1] in "+-" and l[:3] not in ("+++", "---")]
        changed.append("%s -%d +%d" % (f, sum(l[0] == "-" for l in d), sum(l[0] == "+" for l in d)))
    return ", ".join(changed)


DROP = {   # state_to_problem.py's write of each line, made never to happen
    "zonal": ("state_to_problem.py", "        if zonal is not None:\n            f.write(\"zonal",
              "        if False:\n            f.write(\"zonal"),
    "active": ("state_to_problem.py", "        if active is not None:\n            f.write(\"active",
               "        if False:\n            f.write(\"active"),
}
NO_REFUSAL = [("reversal.py", "    same_problem(args.problem, turn, \"turn\")\n", ""),
              ("reversal.py", "    same_problem(args.problem, ret, \"return\")\n", "")]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default="build")
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--keep", action="store_true")
    args = ap.parse_args()
    build = os.path.join(root, args.build)
    exe = os.path.join(build, "ias15_cft.exe" if os.name == "nt" else "ias15_cft")
    if not os.path.exists(exe):
        exe = os.path.join(build, "ias15_cft")
    wd = os.path.join(build, "check_reversal")
    shutil.rmtree(wd, ignore_errors=True)
    os.makedirs(wd)
    print("check_reversal: %s, built %s" % (os.path.relpath(exe, root),
          time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(os.path.getmtime(exe)))))

    # ---------------------------------------------------------------- 0
    print("check_reversal 0: a record without the lines gives the file it always gave")
    hz = os.path.join(root, "results", "horizon", "reversal")
    pairs = [(r, r[:-len("_fwd.rec")] + "_turn.txt") for r in sorted(os.listdir(hz)) if r.endswith("_fwd.rec")]
    pairs += [(r, r[:-len("_back.rec")] + "_return.txt") for r in sorted(os.listdir(hz)) if r.endswith("_back.rec")]
    same, differ = 0, []
    for rec, prob in pairs:
        out = os.path.join(wd, "old_" + prob)
        rc, o, e = run([PY, os.path.join(here, "state_to_problem.py"), os.path.join(hz, rec), "--out", out, "--last"])
        a = open(os.path.join(hz, prob)).read().splitlines()
        b = open(out).read().splitlines() if rc == 0 else []
        if rc == 0 and len(a) == len(b) and a[:1] + a[2:] == b[:1] + b[2:]:
            same += 1
        else:
            differ.append(prob)
    verdict(len(pairs) == 12 and not differ,
            "%d of the %d committed records of results/horizon/reversal/ rebuild their committed problem, every "
            "line but the one naming the record's path%s" % (same, len(pairs), "; not: " + ", ".join(differ) if differ else ""))

    # the problem: zonal_saturn with Mimas a test particle with a mass
    text = open(os.path.join(root, "data", "problems", "zonal_saturn.txt")).read()
    if text.count("\nactive 2\n") != 1:
        raise SystemExit("check_reversal: data/problems/zonal_saturn.txt has no single `active 2` line")
    prob = os.path.join(wd, "zonal_saturn_active1.txt")
    open(prob, "w", newline="\n").write(text.replace("\nactive 2\n", "\nactive 1\n"))

    # ---------------------------------------------------------------- 1
    print("check_reversal 1: the backward leg runs the problem (zonal_saturn with `active 1`, %d steps of %s s "
          "each way)" % (STEPS, DT))
    turn64 = None
    for fmt in (("fp64",) if args.quick else ("fp64", "fp128", "fp256")):
        rc, out, files = reverse(here, prob, fmt, os.path.join(wd, "real_" + fmt), build)
        if not verdict(rc == 0 and all(files.values()),
                       "%s: tools/reversal.py runs (rc %d), and writes its four files" % (fmt, rc)):
            print(out[-1500:])
            continue
        res, m = checks(prob, fmt, out, files)
        for ok, t in res:
            verdict(ok, "%s: %s" % (fmt, t))
        if fmt == "fp64":
            turn64 = files["_turn.txt"]

    # ---------------------------------------------------------------- 2
    print("check_reversal 2: each line carries weight, at fp64: the turn problem with one line deleted, run back")
    if turn64:
        lines = open(turn64).read().splitlines(True)
        for word in ("zonal", "active"):
            kept = [l for l in lines if l.split()[:1] != [word]]
            t = os.path.join(wd, "turn_without_%s.txt" % word)
            open(t, "w", newline="\n").writelines(kept)
            rc, o, e = run(reversal.leg_command(exe, "fp64", STEPS, [], t, "-" + DT))
            back = os.path.join(wd, "back_without_%s.rec" % word)
            open(back, "w", newline="\n").write(o)
            ret = os.path.join(wd, "return_without_%s.txt" % word)
            rc2, o2, e2 = run([PY, os.path.join(here, "state_to_problem.py"), back, "--out", ret, "--last"])
            m = miss(prob, ret) if rc == 0 and rc2 == 0 else None
            verdict(len(kept) == len(lines) - 1 and m is not None and m > 1000 * BOUND["fp64"],
                    "without its %s line (%d line deleted) the backward leg (rc %d) returns %s of the start, "
                    "%s times the bound" % (word, len(lines) - len(kept), rc, "%.3e" % m if m is not None else "nowhere",
                                            "%.1e" % (m / BOUND["fp64"]) if m is not None else "no"))
    else:
        verdict(False, "part 1 wrote no fp64 turn problem to delete lines from")

    # ---------------------------------------------------------------- 3
    print("check_reversal 3: a dropped line fails (controls from copies of the tools, at fp64)")
    for word in ("zonal", "active"):
        cp = os.path.join(wd, "copy_drops_" + word)
        how = mutate(here, cp, [DROP[word]])
        rc, out, files = reverse(cp, prob, "fp64", os.path.join(wd, "run_drops_" + word), build)
        verdict(rc == 1 and ("lacks its %s line" % word) in out and files["_back.rec"] is None,
                "a copy of state_to_problem.py that drops the %s line (the mutation applied: %s): reversal.py "
                "refuses its turn problem, naming the line, before the backward leg (rc %d)" % (word, how, rc))
        cp2 = os.path.join(wd, "copy_drops_%s_unrefused" % word)
        how2 = mutate(here, cp2, [DROP[word]] + NO_REFUSAL)
        rc, out, files = reverse(cp2, prob, "fp64", os.path.join(wd, "run_drops_%s_unrefused" % word), build)
        if rc == 0 and all(files.values()):
            res, m = checks(prob, "fp64", out, files)
            failed = [t for ok, t in res if not ok]
            line_failed = sum(1 for t in failed if (" %s line " % word) in t)
            verdict(line_failed == 2 and any("backward record" in t for t in failed) and any("return lies" in t for t in failed),
                    "the same with reversal.py's refusal deleted too (%s): it runs (rc 0), and part 1 FAILs %d of its "
                    "%d checks: the %s line missing from both problems, the records' lines apart, and the return %s "
                    "of the start" % (how2, len(failed), len(res), word, "%.3e" % m if m is not None else "nowhere"))
        else:
            verdict(False, "the copy without reversal.py's refusal did not run to its end (rc %d)" % rc)

    if not args.keep:
        shutil.rmtree(wd, ignore_errors=True)
    print("RESULT: %s" % ("FAIL" if FAILS else "PASS"))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
