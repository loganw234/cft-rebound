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
# THE PROBLEM: data/problems/zonal_saturn.txt, with three changes.
# - `active 01` for its `active 2`, so that Mimas (GM 2.5) is a test particle
#   with a mass, written with a leading zero, which the program reads as 1.
#   As committed, the file's one test particle is the massless ring particle,
#   whose reaction is an exact zero when every body is active: its active
#   line changes no bit of a return (measured for entry 52: the return's
#   three body lines identical with and without it), so it could not be held.
# - G = 1 + 3 * 2^-53, with a comment after it on its line: at binary64 a
#   tie, which the program rounds to even, up, to 1 + 2^-51.
# - R = 60,330 km + 2^-38: at binary64 a tie, which rounds to even, down, to
#   60,330 km.
# Each has 54 significant bits, exact at binary128 and binary256. A value
# that binary64 holds cannot show a tool that narrows it; a tie each way
# shows a comparison that does not round as the program does; and the
# leading zero and the comment show one that reads a line otherwise.
#
# 0. THE CARRY ALONE.
#    (a) A record without the lines gives the file it always gave: the
#        committed records of results/horizon/reversal/ (kepler and outer,
#        three formats, forward and back) rebuild their committed turn and
#        return problems, every line but the one that names the record's
#        path.
#    (b) Two steps of the problem at binary128 and at binary256: the turn
#        problem carries G, `active` and the zonal line exactly, the wide G
#        and R included. In --quick too, so a narrowing carry fails there.
#    (c) The hex reader reversal.py and this gate use, check_equivalence's
#        hex_to_fraction (state_to_problem.py parses no hex: it copies the
#        record's text), on known answers: -0 and +0 apart, a negative
#        value, and the widened R.
# 1. THE BACKWARD LEG RUNS THE PROBLEM. At each format (binary64 only in
#    --quick), tools/reversal.py runs 100 steps of 100 s each way (rc 0), and:
#    - its turn and return problems carry G, `active` and the zonal line,
#      value for value as the run's format reads the problem's, read here
#      and not through the tools: hex at binary64 by CPython's correctly
#      rounded division, at binary128 and binary256 exactly; an integer by
#      its value; G by its first field;
#    - the backward record's `# active` and `# zonal` lines are the forward
#      record's, word for word, and name what the file names: N - 1 = 2
#      bodies feel the term, and body 0 takes the reaction of active - 1 = 0;
#    - the return lies within B of the start in both max |dx| / L and
#      max |dv| / V as computed here, and reversal.py prints the same
#      figures: B is 1e-14 at binary64 and 1e-32 at binary128, about 100
#      unit roundoffs, where the arithmetic dominates; 1e-40 at binary256,
#      where IAS15's own reversal error does (it is not time-symmetric:
#      tools/reversal.py's head).
# 2. EACH LINE CARRIES WEIGHT, so the bound resolves it. At binary64 the turn
#    problem with its zonal line deleted, and with its active line deleted,
#    runs back as reversal.py runs its backward leg (its leg_command()), and
#    each return must miss by more than 1000 B in max |dx| / L.
# 3. A DROPPED OR CHANGED LINE FAILS: controls from copies of the tools, each
#    mutation shown to apply. Copies of state_to_problem.py that drop the
#    zonal line, drop the active line, and write J4 one unit up in its last
#    hex digit: reversal.py refuses each turn problem with its own words for
#    that line (rc 1). With reversal.py's two refusals deleted as well, each
#    copy runs to its end (rc 0) and part 1's checks FAIL on its run exactly
#    where the line is wrong: a dropped line missing from both problems, the
#    records' lines apart and the return outside B; a changed J4, both
#    problems' zonal lines and the records' lines.
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
from fractions import Fraction

here = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, here)
from check_equivalence import hex_to_fraction  # noqa: E402
import reversal  # noqa: E402  (leg_command() only: part 2 runs the backward leg as it does)

root = os.path.dirname(here)
PY = sys.executable
STEPS, DT = 100, "100"
BOUND = {"fp64": 1e-14, "fp128": 1e-32, "fp256": 1e-40}
TOOLS = ("reversal.py", "state_to_problem.py", "check_equivalence.py")
WORDS = ("G", "active", "zonal")       # the lines the gate holds, after the bodies
FIRST_ONLY = ("G", "E", "system")      # read_problem() reads only their first field

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


def as_read(word, tokens, fmt):
    """A line's values as a run at fmt reads them: G, E and a system's name
    by their first field; hex by value and sign, at binary64 through
    CPython's correctly rounded int division (Fraction to float), wider
    exactly (every value this gate writes fits binary128, asserted where it
    is made); a decimal integer by its value; anything else as text."""
    out = []
    for t in (tokens[:1] if word in FIRST_ONLY else tokens):
        if t.lstrip("+-")[:2].lower() == "0x":
            v, neg = hex_to_fraction(t)
            out.append((Fraction(float(v)) if fmt == "fp64" else v, neg))
        elif re.fullmatch(r"[+-]?[0-9]+", t):
            out.append(int(t))
        else:
            out.append(t)
    return out


def hexs(v):
    """A positive dyadic Fraction as exact canonical hex."""
    e = v.numerator.bit_length() - v.denominator.bit_length()
    if Fraction(2) ** e > v:
        e -= 1
    f, digits = v / Fraction(2) ** e - 1, ""
    while f:
        f *= 16
        digits += "%x" % int(f)
        f -= int(f)
    return "0x1%sp%+d" % ("." + digits if digits else "", e)


def bits(v):
    """The significant bits of a positive dyadic Fraction."""
    n = v.numerator
    return (n >> ((n & -n).bit_length() - 1)).bit_length()


def bodies(path):
    return [(t[1], [hex_to_fraction(v)[0] for v in t[2:9]]) for t in words(path, "body")]


def miss(start, ret):
    """The return against the start: max |dx| / L and max |dv| / V, L and V
    the start's largest coordinate and velocity, computed here exactly and
    only then rounded; None if the bodies are not the start's."""
    b0, b1 = bodies(start), bodies(ret)
    if [n for n, _ in b0] != [n for n, _ in b1]:
        return None
    L = max(abs(q) for _, u in b0 for q in u[1:4])
    V = max(abs(q) for _, u in b0 for q in u[4:7]) or Fraction(1)
    dx = max(abs(u[c] - v[c]) for (_, u), (_, v) in zip(b0, b1) for c in (1, 2, 3))
    dv = max(abs(u[c] - v[c]) for (_, u), (_, v) in zip(b0, b1) for c in (4, 5, 6))
    return float(dx / L), float(dv / V)


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


def carried(problem, path, fmt, word):
    """Whether path carries the problem's `word` line as a run at fmt
    reads it, and the line path has."""
    want = [as_read(word, t[1:], fmt) for t in words(problem, word)]
    got = [as_read(word, t[1:], fmt) for t in words(path, word)] if path else []
    line = " ".join(words(path, word)[0]) if path and words(path, word) else "MISSING"
    return got == want and len(want) == 1, line


def agrees(printed, m):
    """reversal.py's printed figure is the gate's, to its four digits."""
    return printed is not None and abs(float(printed) - m) <= 1e-3 * m


def checks(problem, fmt, out, files):
    """Part 1's checks on one run: {key: (ok, text)}, and the miss."""
    res = {}
    for what in ("turn", "return"):
        for word in WORDS:
            ok, line = carried(problem, files["_%s.txt" % what], fmt, word)
            res["%s %s" % (what, word)] = (ok, "the %s problem's %s line is the problem's as %s reads it: %s"
                                           % (what, word, fmt, line))
    hf = header(files["_fwd.rec"]) if files["_fwd.rec"] else []
    hb = header(files["_back.rec"]) if files["_back.rec"] else []
    N = len(words(problem, "body"))
    n_act = int(words(problem, "active")[0][1])
    k = int(words(problem, "zonal")[0][1])
    named = (len(hf) == 2 and hf[0].startswith("# active %d of %d bodies a system:" % (n_act, N))
             and re.search(r"it acts on %d bodies a system, and body %d takes the reaction of %d$"
                           % (N - 1, k, n_act - 1), hf[1]))
    res["records"] = (hf == hb and bool(named),
                      "the backward record's # active and # zonal lines are the forward record's, and name the file's "
                      "counts (N - 1 = %d feel the term, body %d takes the reaction of active - 1 = %d): %s"
                      % (N - 1, k, n_act - 1, "the same" if hf == hb else "%r against %r" % (hb, hf)))
    m = miss(problem, files["_return.txt"]) if files["_return.txt"] else None
    px = re.search(r"max\|dx\|/L=(\S+)", out)
    pv = re.search(r"max\|dv\|/V=(\S+)", out)
    ok = (m is not None and max(m) <= BOUND[fmt]
          and agrees(px and px.group(1), m[0]) and agrees(pv and pv.group(1), m[1]))
    res["return"] = (ok, "the return lies %s of the start in max |dx| / L and %s in max |dv| / V (the bound is "
                     "%.0e), and reversal.py prints %s and %s"
                     % ("%.3e" % m[0] if m else "nowhere", "%.3e" % m[1] if m else "nowhere", BOUND[fmt],
                        px.group(1) if px else "nothing", pv.group(1) if pv else "nothing"))
    return res, m


def mutate(src_dir, dst_dir, edits):
    """A copy of the tools with each (file, old, new) applied exactly once;
    each edit's lines removed and added, or SystemExit if one does not apply."""
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


# state_to_problem.py's writes, made to drop a line or to change a value;
# each control's refusal, in reversal.py's own words, and the part 1 checks
# its run must FAIL once reversal.py's refusals are gone (the others hold)
WRITE_ZONAL = "            f.write(\"zonal %s %s %s %s\\n\" % zonal)\n"
CONTROLS = [
    ("drops the zonal line", "it lacks its zonal line",
     ("state_to_problem.py", "        if zonal is not None:\n            f.write(\"zonal",
      "        if False:\n            f.write(\"zonal"),
     {"turn zonal", "return zonal", "records", "return"}),
    ("drops the active line", "it lacks its active line",
     ("state_to_problem.py", "        if active is not None:\n            f.write(\"active",
      "        if False:\n            f.write(\"active"),
     {"turn active", "return active", "records", "return"}),
    ("writes J4 one unit up in its last hex digit", "its zonal line(s) differ",
     ("state_to_problem.py", WRITE_ZONAL,
      "            j4 = zonal[3]; i = j4.index(\"p\") - 1\n"
      "            f.write(\"zonal %s %s %s %s\\n\" % (zonal[:3] + (j4[:i] + \"%x\" % ((int(j4[i], 16) + 1) % 16) + j4[i + 1:],)))\n"),
     {"turn zonal", "return zonal", "records"}),
]
NO_REFUSAL = [("reversal.py", "    same_problem(args.problem, turn, \"turn\", args.format)\n", ""),
              ("reversal.py", "    same_problem(args.problem, ret, \"return\", args.format)\n", "")]


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

    # the problem: zonal_saturn with Mimas a test particle with a mass, its
    # active count written `01`, and G and R binary64 ties, each way, within
    # binary128, G with a comment after it
    lines = open(os.path.join(root, "data", "problems", "zonal_saturn.txt")).read().splitlines(True)
    ga = [i for i, l in enumerate(lines) if l.split()[:1] == ["G"]]
    za = [i for i, l in enumerate(lines) if l.split()[:1] == ["zonal"]]
    aa = [i for i, l in enumerate(lines) if l.split() == ["active", "2"]]
    if not (len(ga) == len(za) == len(aa) == 1):
        raise SystemExit("check_reversal: data/problems/zonal_saturn.txt has not one G, zonal and `active 2` line")
    g0 = hex_to_fraction(lines[ga[0]].split()[1])[0]
    z = lines[za[0]].split()
    r0 = hex_to_fraction(z[2])[0]
    wg = g0 + 3 * Fraction(1, 2 ** 53) * g0                 # 1.5 binary64 units above G: a tie, to even, up
    wr = r0 + Fraction(1, 2 ** 38)                          # 0.5 binary64 unit above R: a tie, to even, down
    up, down = Fraction(float(wg)), Fraction(float(wr))
    if not (53 < bits(wg) <= 113 and 53 < bits(wr) <= 113 and up == g0 + 4 * Fraction(1, 2 ** 53) * g0
            and down == r0):
        raise SystemExit("check_reversal: G (%d bits) and R (%d bits) are not binary64 ties, up and down, within "
                         "binary128" % (bits(wg), bits(wr)))
    lines[ga[0]] = "G %s # 1 + 3 * 2^-53: a tie at binary64, to even, up\n" % hexs(wg)
    lines[za[0]] = "zonal %s %s %s %s\n" % (z[1], hexs(wr), z[3], z[4])
    lines[aa[0]] = "active 01\n"
    prob = os.path.join(wd, "zonal_saturn_active1_wide.txt")
    open(prob, "w", newline="\n").writelines(lines)

    # ---------------------------------------------------------------- 0
    print("check_reversal 0: the carry alone")
    hz = os.path.join(root, "results", "horizon", "reversal")
    pairs = [(r, r[:-len("_fwd.rec")] + "_turn.txt") for r in sorted(os.listdir(hz)) if r.endswith("_fwd.rec")]
    pairs += [(r, r[:-len("_back.rec")] + "_return.txt") for r in sorted(os.listdir(hz)) if r.endswith("_back.rec")]
    same, differ = 0, []
    for rec, p in pairs:
        out = os.path.join(wd, "old_" + p)
        rc, o, e = run([PY, os.path.join(here, "state_to_problem.py"), os.path.join(hz, rec), "--out", out, "--last"])
        a = open(os.path.join(hz, p)).read().splitlines()
        b = open(out).read().splitlines() if rc == 0 else []
        if rc == 0 and len(a) == len(b) and a[:1] + a[2:] == b[:1] + b[2:]:
            same += 1
        else:
            differ.append(p)
    verdict(len(pairs) == 12 and not differ,
            "(a) %d of the %d committed records of results/horizon/reversal/ rebuild their committed problem, every "
            "line but the one naming the record's path%s" % (same, len(pairs), "; not: " + ", ".join(differ) if differ else ""))
    for fmt in ("fp128", "fp256"):
        rc, o, e = run(reversal.leg_command(exe, fmt, 2, [], prob, DT))
        rec = os.path.join(wd, "two_steps_%s.rec" % fmt)
        open(rec, "w", newline="\n").write(o)
        turn = os.path.join(wd, "two_steps_%s_turn.txt" % fmt)
        rc2, o2, e2 = run([PY, os.path.join(here, "state_to_problem.py"), rec, "--out", turn, "--last"])
        held = [carried(prob, turn if rc == 0 and rc2 == 0 else None, fmt, w) for w in WORDS]
        verdict(all(ok for ok, _ in held),
                "(b) %s, two steps (rc %d, %d): the turn problem carries G, `active` and the zonal line exactly: %s"
                % (fmt, rc, rc2, "; ".join(line for _, line in held)))
    known = [("-0x0p+0", (0, True)), ("0x0p+0", (0, False)), ("-0x1.8p+1", (-3, True)), (hexs(wr), (wr, False))]
    verdict(all(hex_to_fraction(s) == v for s, v in known),
            "(c) hex_to_fraction reads -0x0p+0 as -0, 0x0p+0 as +0, -0x1.8p+1 as -3 and %s as 60330 + 2^-38"
            % hexs(wr))

    # ---------------------------------------------------------------- 1
    print("check_reversal 1: the backward leg runs the problem (zonal_saturn with `active 01`, G and R binary64 "
          "ties, %d steps of %s s each way)" % (STEPS, DT))
    turn64 = None
    for fmt in (("fp64",) if args.quick else ("fp64", "fp128", "fp256")):
        rc, out, files = reverse(here, prob, fmt, os.path.join(wd, "real_" + fmt), build)
        if not verdict(rc == 0 and all(files.values()),
                       "%s: tools/reversal.py runs (rc %d), and writes its four files" % (fmt, rc)):
            print(out[-1500:])
            continue
        res, m = checks(prob, fmt, out, files)
        for key in res:
            verdict(res[key][0], "%s: %s" % (fmt, res[key][1]))
        if fmt == "fp64":
            turn64 = files["_turn.txt"]

    # ---------------------------------------------------------------- 2
    print("check_reversal 2: each line carries weight, at fp64: the turn problem with one line deleted, run back")
    if turn64:
        tl = open(turn64).read().splitlines(True)
        for word in ("zonal", "active"):
            kept = [l for l in tl if l.split()[:1] != [word]]
            t = os.path.join(wd, "turn_without_%s.txt" % word)
            open(t, "w", newline="\n").writelines(kept)
            rc, o, e = run(reversal.leg_command(exe, "fp64", STEPS, [], t, "-" + DT))
            back = os.path.join(wd, "back_without_%s.rec" % word)
            open(back, "w", newline="\n").write(o)
            ret = os.path.join(wd, "return_without_%s.txt" % word)
            rc2, o2, e2 = run([PY, os.path.join(here, "state_to_problem.py"), back, "--out", ret, "--last"])
            m = miss(prob, ret) if rc == 0 and rc2 == 0 else None
            verdict(len(kept) == len(tl) - 1 and m is not None and m[0] > 1000 * BOUND["fp64"],
                    "without its %s line (%d line deleted) the backward leg (rc %d) returns %s of the start in "
                    "max |dx| / L, %s times the bound" % (word, len(tl) - len(kept), rc,
                                                         "%.3e" % m[0] if m else "nowhere",
                                                         "%.1e" % (m[0] / BOUND["fp64"]) if m else "no"))
    else:
        verdict(False, "part 1 wrote no fp64 turn problem to delete lines from")

    # ---------------------------------------------------------------- 3
    print("check_reversal 3: a dropped or changed line fails (controls from copies of the tools, at fp64)")
    for n, (what, said, edit, must_fail) in enumerate(CONTROLS):
        cp = os.path.join(wd, "copy_%d" % n)
        how = mutate(here, cp, [edit])
        rc, out, files = reverse(cp, prob, "fp64", os.path.join(wd, "run_%d" % n), build)
        verdict(rc == 1 and said in out and files["_back.rec"] is None,
                "a copy of state_to_problem.py that %s (the mutation applied: %s): reversal.py refuses its turn "
                "problem before the backward leg (rc %d), saying \"%s\"" % (what, how, rc, said))
        cp2 = os.path.join(wd, "copy_%d_unrefused" % n)
        how2 = mutate(here, cp2, [edit] + NO_REFUSAL)
        rc, out, files = reverse(cp2, prob, "fp64", os.path.join(wd, "run_%d_unrefused" % n), build)
        if rc == 0 and all(files.values()):
            res, m = checks(prob, "fp64", out, files)
            failed = {k for k, (ok, _) in res.items() if not ok}
            verdict(failed == must_fail,
                    "the same with reversal.py's two refusals deleted too (%s): it runs (rc 0), and part 1 FAILs "
                    "exactly %s of its %d checks (the return %s of the start in max |dx| / L)%s"
                    % (how2, ", ".join(sorted(must_fail)), len(res), "%.3e" % m[0] if m else "nowhere",
                       "" if failed == must_fail else "; it failed " + ", ".join(sorted(failed))))
        else:
            verdict(False, "the copy without reversal.py's refusals did not run to its end (rc %d)" % rc)

    if not args.keep:
        shutil.rmtree(wd, ignore_errors=True)
    print("RESULT: %s" % ("FAIL" if FAILS else "PASS"))
    return 1 if FAILS else 0


if __name__ == "__main__":
    sys.exit(main())
