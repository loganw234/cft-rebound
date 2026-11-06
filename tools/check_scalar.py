#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
# Copyright 2026 the cft-rebound contributors.
#
# The gate of the scalar route (ROUND3.md, P5).
#
# A constant operand of an elementwise call is a vector of N3 copies of one
# value. With the scalar route on it goes to the library as ONE element,
# through the door's vrun_k(). The route is opt-in (src/ias15_cft.c,
# route_def[]): unforced it is off, by=opt_in, even where the device publishes
# the bit, because no card has measured it, and each site is then today's
# call on the whole vector (vmulk() and its kin). So nothing but a forced run
# takes the new route, and every run here says which one it wants.
#
#   python tools/check_scalar.py [--build build] [--quick] [--library]
#                                [--exe ias15_cft] [--source FILE]
#   (every tier passes --library; check-quick and check-light pass --quick,
#   which is binary64 only. No case's step count is scaled: the rejections
#   and the mixed attempt are what those counts are there for)
#
# What it asserts, and what lets each assertion fail:
#   1. THE RECOUNT, from the source. Every site that hands the door a
#      constant, by function, constant and width, must be a row of SITES
#      below, and each row must be as many sites as it says. A site's scalar
#      must be a constant: a name make_constants() fills or vbcast() broadcasts, and
#      never one bcast_sys() fills, which holds one value a SYSTEM. Its
#      result must not be its scalar, which is the door's read-once rule
#      made moot. A constant left on a plain call must be one of the kinds
#      this parcel leaves, each of which is printed: a divide, a compare, or
#      a site in dtnew_*() and choose_timestep(), which are nobody's this
#      round. Anything else is a site the conversion missed.
#      And each helper must be exactly its two statements: today's call on its
#      own arguments with the route off, counting the fallback, then vrun_k(),
#      with no site calling vrun_k() itself. Its fallback copies element 0 into
#      scratch. So a site through it, or a helper that copied first, would
#      cost the default route that copy and blind the comparison below, and no
#      run could tell: the bits, the counts and the bytes are the same either
#      way. That is why the text is held whole (V5's H4, H6 and H7, 2026-11-06).
#   2. THE MATRIX. Every problem, at every format, on both engines and both
#      forms of the arithmetic, on three routes: default, forced off and
#      forced on.
#        - The records, the schedule of corrector passes, and every counter
#          but staged_bytes are the same on all three. Timings are not
#          counters. The trailer's keys are read, never listed.
#        - The `# route scalar=` line: calls above zero and no fallback
#          when on, the reverse when off, and off by=opt_in in the default,
#          which only a published bit can say. The other routes say the
#          same in all three.
#        - calls= is DERIVED. SITES gives each site's calls a pass, an
#          attempt, a prediction or an energy evaluation from the
#          integrator's shape, seven nodes and seven levels. The run's own
#          schedule (--record-iters) gives its attempts and passes, and a
#          one-step run gives the attempts before any system was accepted,
#          which predict nothing. On the program engine that prediction
#          count must be the trailer's nextstep_calls: the engine's own word
#          on the same number.
#        - staged_bytes falls by exactly those calls times ESZ * (width - 1),
#          and only when the route is on. A site that fell out of the door
#          moves both.
#        - Beside each, the physics that case is there for (below).
#   3. THE TRAP. DTB, DTDB, RB, D1B and D2B look exactly like scalars with
#      one system, and with systems that happen to share a step. In any
#      other ensemble element 0 is system 0's step, applied to every system.
#      Two oracles see it, and every ensemble here is held to both:
#        - The routes, part 2's comparison. With the route off a site reads
#          the whole vector, and with it on element 0 alone, so a per-system
#          vector gives different records on the two.
#        - The members. Forced on, each system must equal itself run alone
#          (--member k), where element 0 is its own. This one also holds
#          what the two routes share.
#      Only a run with the route forced on can see it, because today's route
#      stays right under this defect. check_seams forces it on an adaptive
#      ensemble and fails DTDB, DTB and RB's powers. It never runs --engine
#      loop --arith fma with the route on, so D1B and D2B have this gate alone.
#      The case that exists for the trap must SHOW its systems' steps
#      differing, and a mixed attempt. The fixed-step ensemble must show them
#      equal: it is the inert case (METHOD.md 5), there so that it is not
#      taken for a test.
#   4. --library. The library half takes the route from the same door, and
#      the route being opt-in, the suite's own check_dropin leg runs today's.
#      So a library half that ignored the forced route passes that leg.
#      Here: check_dropin --light with the route forced on passes against
#      REBOUND, and its own scalar topic reads on by=forced through
#      ias15_engine_route(), with calls and no fallback. Then gate_real's
#      uninterrupted run at each format, forced on against forced off: the
#      same dump, wide state included, each side's route read.
#
# Neither oracle sees a site given the wrong constant, which both routes and
# E = 1 share. The suite's REBOUND gates do, because they run the same sites
# on today's route.
#
# Stated, not closed: a route that says it is on and is not passes every tier.
# A vrun_k() whose on branch ran today's cft_run, still counting route_calls and
# tallying one element, gives the same records, the same `# route` line and the
# same fall in staged_bytes, in either half or both. On software both routes
# compute the same bits, and calls= and staged_bytes are the door's own tallies,
# not anything libcft reports, so nothing here reads whether cft_run_ex ran.
# (V5's m4_lie_both, measured on 2026-11-06: check_scalar --quick --library 32 ok,
# check_door 50 ok. It is P1's m4_lie for this route. The M4 forms this gate
# does fail leave a counter or the route line honest.)
#
# A device without the bit refuses the route forced on, here as in
# check_gravity.py, so this gate needs one that publishes it.
#
# The route names and their variables are read from the program's own
# output, never typed here, except the one name this gate exists for.

import argparse
import concurrent.futures
import collections
import os
import re
import subprocess
import sys

here = os.path.dirname(os.path.abspath(__file__))
root = os.path.dirname(here)

ROUTE = "scalar"
ROUTE_LINE = re.compile(r"^door route (\w+)=(on|off) by=(\w+) env=(\S+) bit=(\S+)$")
TRAILER_ROUTE = re.compile(r"^# route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)$")
GATE_ROUTE = re.compile(r"^gate_real: route (\w+)=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)\s*$", re.M)
HEADER = re.compile(r"^# program=.*?(?: E=(\d+))? N=(\d+) ", re.M)

# REBOUND's seven Gauss-Radau nodes past the start, and the seven levels of
# its b polynomial. Every count below is derived from these two.
NODES, LEVELS = 7, 7
NODE = range(1, NODES + 1)
TRIANGLE = sum(LEVELS - 1 - m for m in range(LEVELS))     # j > m over the levels: KD's and KBIN's pairs

# Every site, by the function it is in and the constant it passes:
#   (function, constant, width, arithmetic, engine, once per, calls then, sites)
# The last is how many call sites in the source make those calls. The arithmetic is rebound, fma or both. The engine is loop, or both where
# the code is outside the programs. library is the library half alone, which
# has no program engine, and standalone is the standalone program alone, on
# either engine.
#   per pass       a corrector pass: every node's predictor and corrector
#   per attempt    a step attempt, accepted or not
#   per predict    an attempt that predicts the next step's e and b
#   per energy     an energy evaluation, one per sample
SITES = (
    ("predict_positions", "KHF", "N3", "fma", "loop", "pass", NODES * LEVELS, 2),
    ("predict_positions", "KNUM", "N3", "rebound", "loop", "pass", NODES * LEVELS, 2),
    ("predict_positions", "KH", "N3", "rebound", "loop", "pass", NODES * (LEVELS + 2), 3),      # each level, dt*h/2 and dt*h
    ("correct", "KRINV", "N3", "fma", "loop", "pass", sum(NODE), 2),                            # n of them at node n
    ("correct", "KC", "N3", "both", "loop", "pass", sum(n - 1 for n in NODE), 1),               # n - 1
    ("predict_next_step_vec", "KBIN", "N3", "both", "loop", "predict", TRIANGLE, 1),
    ("predict_velocities", "KVNUM", "N3", "rebound", "library", "pass", NODES * LEVELS, 2),
    ("predict_velocities", "KH", "N3", "rebound", "library", "pass", NODES * (LEVELS + 1), 2),  # each level and dt*h
    ("step_attempt", "KD", "N3", "both", "both", "attempt", TRIANGLE, 1),
    ("step_attempt", "KH", "N3", "fma", "both", "attempt", NODES, 1),                           # fl(dt h_n)
    ("step_attempt", "KHALF", "N3", "fma", "both", "attempt", NODES, 1),                        # and its half
    ("step_attempt", "KPOSR", "N3", "fma", "both", "attempt", LEVELS + 1, 2),                   # with a0's
    ("step_attempt", "KVELR", "N3", "fma", "both", "attempt", LEVELS, 1),
    ("energy_all", "halfb", "NB", "both", "standalone", "energy", 1, 1),
    ("energy_all", "G", "P", "both", "standalone", "energy", 1, 1),
    ("energy_all", "K0", "E", "both", "standalone", "energy", 1, 1),
)
ENGINES = (("loop", "rebound"), ("loop", "fma"), ("program", "fma"))
WORKERS = 1 if os.environ.get("CFT_REBOUND_ARTIFACT") else 4      # on a card the runs go one at a time

# The parts of the source this recount does not read: the door's own
# wrappers and their self-test, and the reference gravity that every record
# before P1 was made with, which --door-selftest holds gravity to.
# The helpers, whole: each one's parameters, today's call with the route off,
# vrun_k()'s arguments with it on, and which argument at a site is the constant.
HELPERS = {
    "vmulk": ("V d, const V a, const V k, size_t n", "vmul(d, a, k, n)", 'CFT_MUL, d, a, k, NULL, n, 2u, "mul"', 2),
    "vfmak": ("V d, const V a, const V k, const V c, size_t n", "vfma(d, a, k, c, n)", 'CFT_FMA, d, a, k, c, n, 2u, "fma"', 2),
    "vkmul": ("V d, const V k, const V b, size_t n", "vmul(d, k, b, n)", 'CFT_MUL, d, k, b, NULL, n, 1u, "mul"', 1),
    "vaddk": ("V d, const V a, const V k, size_t n", "vadd(d, a, k, n)", 'CFT_ADD, d, a, NULL, k, n, 4u, "add"', 2),
}
NOT_SITES = {"vrun_k", "door_k", "door_selftest", "gravity_reference"} | set(HELPERS)
# energy_all()'s 1/2, a quotient of two constants broadcast one line above its site.
LOCAL_CONSTANTS = {"energy_all": {"halfb"}}
CALL = re.compile(r"\b(vmulk|vfmak|vkmul|vaddk|vrun_k|vmul|vfma|vadd|vsub|vcopysign|vdiv|vcmplt|vcmple)\(")
COMPARES = ("vcmplt", "vcmple")


def strip_comments(src):
    """The source with comments blanked, newlines kept, strings left alone."""
    out, i, n = [], 0, len(src)
    while i < n:
        c = src[i]
        if src.startswith("/*", i):
            j = src.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append("".join(ch if ch == "\n" else " " for ch in src[i:j]))
            i = j
        elif src.startswith("//", i):
            j = src.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c == '"':
            j = i + 1
            while j < n and src[j] != '"':
                j += 2 if src[j] == "\\" else 1
            out.append(src[i:j + 1])
            i = j + 1
        elif c == "'":
            j = i + 1
            while j < n and src[j] != "'":
                j += 2 if src[j] == "\\" else 1
            out.append(src[i:j + 1])
            i = j + 1
        else:
            out.append(c)
            i += 1
    return "".join(out)


def functions(src):
    """Each top-level function: (name, offset of its body, its body)."""
    found, depth, start, name = [], 0, 0, None
    for i, c in enumerate(src):
        if c == "{":
            if depth == 0:
                m = re.search(r"\b(\w+)\s*\([^;{}()]*(?:\([^()]*\)[^;{}()]*)*\)\s*$", src[max(0, i - 400):i])
                name, start = (m.group(1) if m else None), i
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0 and name:
                found.append((name, start, src[start:i + 1]))
                name = None
    return found


def arguments(text, at):
    """The top-level comma-separated arguments of the call open at text[at - 1]."""
    args, depth, cur, i = [], 0, [], at
    while i < len(text):
        c = text[i]
        if c == '"':
            j = i + 1
            while text[j] != '"':
                j += 2 if text[j] == "\\" else 1
            cur.append(text[i:j + 1])
            i = j + 1
            continue
        if c in "([{":
            depth += 1
        elif c in ")]}":
            if depth == 0:
                args.append("".join(cur).strip())
                return args
            depth -= 1
        if c == "," and depth == 0:
            args.append("".join(cur).strip())
            cur = []
        else:
            cur.append(c)
        i += 1
    raise SystemExit("check_scalar: an unterminated call in the source")


def base(token):
    m = re.match(r"\s*(\w+)", token)
    return m.group(1) if m else ""


def squeezed(text):
    return re.sub(r"\s+", "", text)


def recount(path):
    """Everything the source says about the sites."""
    with open(path, encoding="utf-8") as f:
        raw = f.read()
    src = strip_comments(raw)
    funcs = functions(src)
    bodies = {n: b for n, _, b in funcs}
    if "make_constants" not in bodies or "force_constants" not in bodies:
        raise SystemExit("check_scalar: no make_constants() or force_constants() in %s" % path)
    # What is a constant is read from how it is filled: by the two functions
    # that derive them, or by vbcast(), one value into every lane.
    consts = set()
    for fn in ("make_constants", "force_constants"):
        consts |= set(re.findall(r"\b([A-Z]\w*)(?:\[[^\]]*\])*\s*=[^=;][^;]*?\b(?:kderived|kint|kdec|valloc)\(", bodies[fn]))
    consts |= {n for n in re.findall(r"\bvbcast\(\s*(\w+)\s*[,\[]", src) if n[0].isupper()}
    per_system = set(re.findall(r"\bbcast_sys\(\s*(\w+)\s*,", src))
    sites, left, bad = [], [], []
    aliased = 0
    # Each helper, whole, whitespace aside: its parameters and two statements.
    for name, (params, plain, door, _) in HELPERS.items():
        head = re.search(r"\bstatic\s+void\s+%s\s*\(([^)]*)\)\s*\{" % name, src)
        want = "{ if (!route_on[ROUTE_SCALAR]){ route_fallbacks[ROUTE_SCALAR]++; %s; return; } vrun_k(%s); }" % (plain, door)
        if not head or squeezed(head.group(1)) != squeezed(params) or squeezed(bodies.get(name, "")) != squeezed(want):
            bad.append("%s() is not (%s) %s: anything else can copy on the default route, or read element 0 on both, "
                       "and no run would show it" % (name, params, want))
    for name, off, body in funcs:
        if name in NOT_SITES:
            continue
        alias = dict(re.findall(r"\bV\s+(\w+)\s*=\s*(K\w+)\[", body))
        known = consts | set(alias) | LOCAL_CONSTANTS.get(name, set())
        for m in CALL.finditer(body):
            line = src.count("\n", 0, off + m.start()) + 1
            call = m.group(1)
            args = arguments(body, m.end())
            if call == "vrun_k":
                bad.append("%s:%d: a site calls vrun_k() itself, whose fallback copies into scratch: the default route "
                           "must be today's call, through a helper" % (name, line))
                continue
            if call in HELPERS:
                pos = HELPERS[call][3]
                k, d, width = args[pos], args[0], args[-1]
                b = base(k)
                konst = alias.get(b, b)
                if b not in known:
                    why = " It is filled by bcast_sys(), one value a SYSTEM." if b in per_system else ""
                    bad.append("%s:%d: %s passes %s as a scalar, and it is not a constant.%s" % (name, line, call, b, why))
                if base(d) == b:
                    bad.append("%s:%d: %s's result is its own scalar, %s" % (name, line, call, b))
                if any(base(o) == base(d) for i, o in enumerate(args[1:-1], 1) if i != pos):
                    aliased += 1
                sites.append((name, konst, width, line))
                continue
            width = args[-1]
            if width == "1":
                continue
            hit = [base(o) for o in args[1:-1] if base(o) in known]
            if not hit:
                continue
            if name.startswith("dtnew_") or name == "choose_timestep":
                left.append(("dtnew_*() and choose_timestep(), which are nobody's this round", name, hit[0], width, call, line))
            elif call == "vdiv":
                left.append(("a divide: cft_div has no scalar form, and is not an opcode cft_run_ex can take", name, hit[0], width, call, line))
            elif call in COMPARES:
                left.append(("a compare: vcmplt() and vcmple() count no staged bytes, so vrun_k() would add to staged_bytes and not take from it", name, hit[0], width, call, line))
            else:
                bad.append("%s:%d: %s(%s) passes the constant %s as a vector: a site the conversion missed"
                           % (name, line, call, ", ".join(args), hit[0]))
    return sites, left, bad, consts, per_system, aliased


def rows_for(engine, arith):
    """The sites the standalone program reaches on this engine and form."""
    return [r for r in SITES if (r[3] == "both" or r[3] == arith) and
            (r[4] in ("both", "standalone") or r[4] == engine)]


def derive(engine, arith, st):
    """calls= and the fall in staged_bytes this run's shape requires."""
    calls = fall = 0
    width = {"N3": st["N3"], "NB": st["NB"], "P": st["P"], "E": st["E"]}
    for r in rows_for(engine, arith):
        w = width[r[2]]
        if not w:                 # energy_all()'s pair site stands behind `if (P)`
            continue
        n = r[6] * st[r[5]]
        calls += n
        fall += n * (w - 1)
    return calls, fall * st["esz"]


def clean_env(extra=None):
    """No door variable inherited: the environment may not choose a route
    behind the gate's back. CFT_REBOUND_ARTIFACT stays."""
    env = {k: v for k, v in os.environ.items()
           if not k.upper().startswith("CFT_REBOUND_ROUTE_") and not k.upper().startswith("CFT_REBOUND_HIDE")}
    env.update(extra or {})
    return env


def run(cmd, extra=None):
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, env=clean_env(extra), cwd=root)
    return p.returncode, p.stdout, p.stderr


class Record:
    """One integration, read."""

    def __init__(self, cft, fmt, prob, args, env, sched):
        cmd = [cft, "--format", fmt, "--problem", prob, "--quiet"] + args
        if sched:
            cmd += ["--record-iters", sched]
        self.rc, self.out, self.err = run(cmd, env)
        self.recs = [l for l in self.out.splitlines() if l and not l.startswith("#")]
        self.tail, self.systems, self.routes = {}, {}, {}
        for line in self.out.splitlines():
            if line.startswith("# steps_done="):
                for kv in line[2:].split():
                    k, _, v = kv.partition("=")
                    self.tail[k] = v
            m = re.match(r"^# system (\d+) (steps_done=.*)$", line)
            if m:
                self.systems[int(m.group(1))] = dict(kv.split("=") for kv in m.group(2).split())
            m = TRAILER_ROUTE.match(line)
            if m:
                self.routes[m.group(1)] = (m.group(2), m.group(3), int(m.group(4)), int(m.group(5)))
        m = HEADER.search(self.out)
        self.E, self.N = (int(m.group(1) or 1), int(m.group(2))) if m else (0, 0)
        self.sched = []
        if sched and self.rc == 0:
            with open(sched) as f:
                self.sched = [l.split() for l in f if not l.startswith("#")]

    def counters(self):
        """Every trailer key but the timings and staged_bytes."""
        return {k: v for k, v in self.tail.items()
                if k != "staged_bytes" and not k.startswith("t_") and not k.endswith(("seconds", "_per_s"))}


def passes(sched):
    """The shared loop's passes: an attempt runs until its last system leaves."""
    return sum(max(int(f.rstrip("c")) for f in fields if f != "-") for fields in sched)


def before_first_accepted(sched):
    """From a one-step run's schedule: the attempts before any system was
    accepted. A system accepted at attempt i is out of attempt i + 1."""
    for i in range(len(sched) - 1):
        if "-" in sched[i + 1]:
            return i
    return len(sched) - 1


def spread_of(rec):
    """The samples past the first at which the systems' dt_next differ by
    more than any rounding could hide, how many samples there are, and the
    widest pair. Read at binary64's reach, which is plenty: float.fromhex
    rounds a wide format's digits as it reads them. Sample 0 is dt0, the same
    for every system by construction. A line is `system s sample k step t
    dt_next dt_last E ...`."""
    by = collections.defaultdict(list)
    for line in rec.recs:
        t = line.split()
        by[int(t[3])].append(float.fromhex(t[6]))
    n, total, widest = 0, 0, (1.0, 1.0)
    for k, vals in sorted(by.items()):
        if k == 0:
            continue
        total += 1
        lo, hi = min(vals), max(vals)
        if hi - lo > lo * 2.0 ** -40:
            n += 1
            if hi / lo > widest[1] / widest[0]:
                widest = (lo, hi)
    return n, total, widest


def per_system_rejected(rec):
    return [int(rec.systems[s]["steps_rejected"]) for s in sorted(rec.systems)]


class Case:
    def __init__(self, name, prob, kind, args):
        self.name, self.prob, self.kind, self.args = name, prob, kind, args


def check_case(cft, fmt, engine, arith, case, env_name, want_default, work, progs):
    """One case on one engine at one format: a list of failures and a line
    saying what it held."""
    bad = []
    extra = list(case.args) + ["--arith", arith, "--engine", engine]
    if engine == "program":
        extra += ["--programs", progs]
    tag = "%s-%s-%s-%s" % (re.sub(r"\W+", "_", case.name), fmt, engine, arith)
    runs = {}
    for label, env in (("default", {}), ("off", {env_name: "off"}), ("on", {env_name: "on"})):
        runs[label] = Record(cft, fmt, case.prob, extra, env, os.path.join(work, "%s-%s.sched" % (tag, label)))
        if runs[label].rc:
            bad.append("%s: rc %d: %s" % (label, runs[label].rc, runs[label].err.strip()[-240:]))
    if bad:
        return bad, ""
    d, off, on = runs["default"], runs["off"], runs["on"]

    # 2. The records, the schedule, every counter.
    for label, r in (("default", d), ("on", on)):
        if r.recs != off.recs or not r.recs:
            bad.append("the record %s differs from forced off's in %d of %d lines%s"
                       % (label, sum(a != b for a, b in zip(r.recs, off.recs)) + abs(len(r.recs) - len(off.recs)), len(off.recs),
                          ": today's route reads each lane, and the new one element 0" if label == "on" else ""))
        if r.sched != off.sched or not r.sched:
            bad.append("the corrector's schedule %s differs from forced off's" % label)
        if r.counters() != off.counters():
            bad.append("the counters %s differ from forced off's: %s" % (label, ", ".join(
                "%s %s against %s" % (k, r.counters().get(k), off.counters().get(k))
                for k in sorted(set(r.counters()) | set(off.counters()))
                if r.counters().get(k) != off.counters().get(k))))
        if r.systems != off.systems:
            bad.append("a system's trailer line %s differs from forced off's" % label)

    # The route lines.
    for label, r, want in (("default", d, want_default), ("off", off, ("off", ("forced",))), ("on", on, ("on", ("forced",)))):
        if ROUTE not in r.routes:
            bad.append("%s: no # route %s line" % (label, ROUTE))
            continue
        state, by, calls, falls = r.routes[ROUTE]
        ran = (calls > 0 and falls == 0) if want[0] == "on" else (falls > 0 and calls == 0)
        if state != want[0] or by not in want[1] or not ran:
            bad.append("%s: # route %s=%s by=%s calls=%d fallbacks=%d (wanted %s by=%s, and the counter of the route "
                       "that ran above zero, the other zero)" % (label, ROUTE, state, by, calls, falls, want[0], "/".join(want[1])))
    others = {k: v for k, v in off.routes.items() if k != ROUTE}
    if any({k: v for k, v in r.routes.items() if k != ROUTE} != others for r in (d, on)):
        bad.append("another route's line moved with this one")

    # The derivation.
    one = Record(cft, fmt, case.prob, steps_once(extra), {env_name: "off"}, os.path.join(work, tag + "-one.sched"))
    if one.rc or not one.sched:
        bad.append("the one-step run: rc %d: %s" % (one.rc, one.err.strip()[-200:]))
        return bad, ""
    E, N = off.E, off.N
    attempts, npass = len(off.sched), passes(off.sched)
    skipped = before_first_accepted(one.sched)
    st = {"esz": int(fmt[2:]) // 8, "N3": 3 * N * E, "NB": N * E, "P": E * N * (N - 1) // 2, "E": E,
          "pass": npass, "attempt": attempts, "predict": attempts - skipped,
          "energy": sum(1 for l in off.recs if E == 1 or l.startswith("system 0 "))}
    want_calls, want_fall = derive(engine, arith, st)
    calls = on.routes.get(ROUTE, (0, 0, 0, 0))[2]
    falls = off.routes.get(ROUTE, (0, 0, 0, 0))[3]
    if calls != want_calls or falls != want_calls:
        bad.append("calls through the door: %d on, %d off, and the sites derive %d from %d passes, %d attempts, %d "
                   "predictions and %d energy evaluations" % (calls, falls, want_calls, npass, attempts, st["predict"], st["energy"]))
    s_d, s_off, s_on = int(d.tail.get("staged_bytes", -1)), int(off.tail.get("staged_bytes", -1)), int(on.tail.get("staged_bytes", -1))
    if s_d != s_off:
        bad.append("staged_bytes moved with the route off: %d default, %d forced off" % (s_d, s_off))
    if s_off - s_on != want_fall or want_fall <= 0:
        bad.append("staged_bytes fell by %d, and the sites derive %d" % (s_off - s_on, want_fall))
    if engine == "program" and int(off.tail.get("nextstep_calls", -1)) != st["predict"]:
        bad.append("nextstep_calls=%s, and the schedule derives %d predictions" % (off.tail.get("nextstep_calls"), st["predict"]))
    if engine == "loop" and int(off.tail.get("nextstep_calls", -1)) != 0:
        bad.append("nextstep_calls=%s on the loop engine" % off.tail.get("nextstep_calls"))

    # What this case is here for.
    rejected = int(off.tail.get("steps_rejected", -1))
    said = ["%s steps, %d rejected, %d attempts, %d passes" % (off.tail.get("steps_done"), rejected, attempts, npass)]
    if case.kind == "rejection":
        if rejected <= 0 or st["predict"] <= int(off.tail["steps_done"]):
            bad.append("the physics: %d rejected and %d predictions in %s steps, and this case is here for rejected "
                       "attempts that predict" % (rejected, st["predict"], off.tail["steps_done"]))
        said.append("%d rejected attempts predicted" % (st["predict"] - int(off.tail["steps_done"])))
    if case.kind == "first rejected":
        if skipped < 1 or rejected < skipped:
            bad.append("the physics: %d attempts before the first accepted, and this case is here for a first step "
                       "rejected, which predicts nothing" % skipped)
        said.append("%d attempts before the first accepted predicted nothing" % skipped)
    if E > 1:
        n, total, wide = spread_of(off)
        if case.kind == "systems":
            rej = per_system_rejected(off)
            if n == 0:
                bad.append("the physics: at no sampled step do the %d systems' steps differ, and this case is the trap's" % E)
            if case.mixed and (not any(rej) or all(rej)):
                bad.append("the physics: steps_rejected=%s, and this case is here for a mixed attempt" % rej)
            said.append("the systems' steps differ at %d of %d samples (%.6g to %.6g at the widest), rejected %s"
                        % (n, total, wide[0], wide[1], rej))
        elif case.kind == "shared":
            if n:
                bad.append("the physics: the systems' steps differ at %d of %d samples, and this case is the one where "
                           "they share a step" % (n, total))
            said.append("the systems share a step, so a per-system vector is inert here")
        # 3. The ensemble against its members, forced on.
        same = 0
        for k in range(E):
            mem = Record(cft, fmt, case.prob, extra + ["--member", str(k)], {env_name: "on"}, None)
            if mem.rc:
                bad.append("member %d alone: rc %d: %s" % (k, mem.rc, mem.err.strip()[-200:]))
                continue
            mine = [l.split(" ", 2)[2] for l in on.recs if l.startswith("system %d " % k)]
            if mine != mem.recs or not mine:
                wrong = sum(1 for a, b in zip(mine, mem.recs) if a != b)
                bad.append("forced on, system %d differs from itself run alone (%d of %d lines): each lane must keep its "
                           "own system's values" % (k, wrong, len(mine)))
                continue
            theirs = on.systems.get(k, {})
            if not theirs or any(theirs[key] != mem.tail.get(key) for key in theirs):
                bad.append("forced on, system %d's trailer is not its own: %s, and alone %s"
                           % (k, theirs, {key: mem.tail.get(key) for key in theirs}))
                continue
            same += 1
        said.append("the two oracles: forced on %s forced off, and %d of %d systems equal themselves run alone, forced on"
                    % ("equals" if on.recs == off.recs else "DIFFERS FROM", same, E))
    # A failing line says what the run did, never the claim it failed.
    derived = calls == want_calls and falls == want_calls and s_off - s_on == want_fall
    said.append("the route carried %d calls and staged_bytes fell by %d%s" % (calls, s_off - s_on, ", both as derived" if derived else ""))
    return bad, "; ".join(said)


def steps_once(extra):
    """The same run, stopped at its first accepted step."""
    out, i = [], 0
    while i < len(extra):
        if extra[i] in ("--steps", "--sample"):
            out += [extra[i], "1"]
            i += 2
        else:
            out.append(extra[i])
            i += 1
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--build", default=os.path.join(root, "build"))
    ap.add_argument("--exe", default="ias15_cft")
    ap.add_argument("--source", default=os.path.join(root, "src", "ias15_cft.c"))
    ap.add_argument("--quick", action="store_true", help="binary64 only")
    ap.add_argument("--library", action="store_true", help="also the drop-in's suite and gate_real, the route forced on")
    args = ap.parse_args()
    exe = ".exe" if os.name == "nt" else ""
    cft = os.path.join(args.build, args.exe + exe)
    formats = ["fp64"] if args.quick else ["fp64", "fp128", "fp256"]
    software = not os.environ.get("CFT_REBOUND_ARTIFACT")
    progs = os.path.join(root, "programs", "out")
    allok = True

    def verdict(ok, text):
        nonlocal allok
        allok &= ok
        print("  %s %s" % ("ok  " if ok else "FAIL", text))

    # The route, from the program. A self-test with every bit hidden runs
    # no case at all, and still prints each route's variable.
    rc, out, err = run([cft, "--probe", "--format", "fp64"])
    names = re.findall(r"^route (\w+): ", out, re.M)
    rc, out, err = run([cft, "--door-selftest", "--format", "fp64"], {"CFT_REBOUND_HIDE_FEATURES": ",".join(sorted(names))})
    routes = {}
    for line in out.splitlines():
        m = ROUTE_LINE.match(line)
        if m:
            routes[m.group(1)] = {"env": m.group(4), "bit": m.group(5)}
    if ROUTE not in routes or sorted(routes) != sorted(names):
        print("check_scalar: the program does not name the route %s (--probe names %s, the self-test %s; rc %d): %s"
              % (ROUTE, ", ".join(names) or "none", ", ".join(sorted(routes)) or "none", rc, (err or out)[-300:].strip()))
        print("RESULT: FAIL")
        return 1
    env_name = routes[ROUTE]["env"]
    want_default = ("off", ("opt_in",) if software else ("opt_in", "absent"))
    print("check_scalar: the route %s (%s, %s)" % (ROUTE, env_name, routes[ROUTE]["bit"]))

    # 1. The recount.
    print("check_scalar: the recount, from %s" % os.path.relpath(args.source, root))
    sites, left, bad, consts, per_system, aliased = recount(args.source)
    have = collections.Counter((s[0], s[1], s[2]) for s in sites)
    model = {(r[0], r[1], r[2]): r[7] for r in SITES}
    missing = sorted(set(model) - set(have))
    strange = sorted(set(have) - set(model))
    uneven = sorted(k for k in set(model) & set(have) if have[k] != model[k])
    held = not bad and not missing and not strange and not uneven
    verdict(held, ("%d sites: each row of the gate's table as many as it says, each scalar a constant and never its own result, "
                   "each helper its two statements, and no constant left on a call that could take it" if held
                   else "%d sites, and the recount found:")
            % len(sites) + " (%d constants; per system, so never a scalar: %s; %d sites' results alias their vector operand, "
            "which cft_run allows)%s%s%s%s"
            % (len(consts), ", ".join(sorted(per_system)) or "none", aliased,
               "".join("\n         " + b for b in bad),
               "".join("\n         the gate's table has %s in %s at %s, and the source does not" % (k, f, w) for f, k, w in missing),
               "".join("\n         the source has %s in %s at %s, and the gate's table does not" % (k, f, w) for f, k, w in strange),
               "".join("\n         the source has %d sites of %s in %s at %s, and the gate's table says %d" % (have[(f, k, w)], k, f, w, model[(f, k, w)])
                       for f, k, w in uneven)))
    for engine, arith in ENGINES:
        rows = rows_for(engine, arith)
        n = sum(have[(r[0], r[1], r[2])] for r in rows)
        per = {p: sum(r[6] for r in rows if r[5] == p) for p in ("pass", "attempt", "predict", "energy")}
        print("       --engine %s --arith %s: %d sites, in calls %d a pass, %d an attempt, %d a prediction, %d an energy evaluation"
              % (engine, arith, n, per["pass"], per["attempt"], per["predict"], per["energy"]))
    lib = [r for r in SITES if r[4] != "standalone" and r[3] != "fma"]
    print("       the library half, REBOUND's form: %d sites, four of them predict_velocities()'s, which only a "
          "velocity-dependent routine reaches" % sum(have[(r[0], r[1], r[2])] for r in lib))
    kinds = collections.OrderedDict()
    for why, fn, k, width, call, line in left:
        kinds.setdefault(why, []).append("%s %s(%s, %s)" % (fn, call, k, width))
    for why, what in kinds.items():
        print("       left, %d sites: %s: %s" % (len(what), why, "; ".join(what)))

    # 2 and 3. The matrix.
    work = os.path.join(args.build, "check_scalar")
    os.makedirs(work, exist_ok=True)
    py = sys.executable
    probs = os.path.join(root, "data", "problems")
    spread = os.path.join(work, "kepler_spread_E5.txt")
    outer3 = os.path.join(work, "outer_jup_E3.txt")
    for cmd in ([py, os.path.join(here, "make_ensemble.py"), os.path.join(probs, "kepler.txt"), "--members", "5", "--out", spread,
                 "--offset", "planet", "x", "0.4"],
                [py, os.path.join(here, "make_ensemble.py"), os.path.join(probs, "outer.txt"), "--members", "3", "--out", outer3,
                 "--offset", "jupiter", "x", "0x1p-20"]):
        subprocess.run(cmd, check=True, cwd=root, stdout=subprocess.DEVNULL)

    def p(name):
        return os.path.join(probs, name)

    cases = [
        Case("kepler, fixed", p("kepler.txt"), "fixed", ["--dt", "0.05", "--epsilon", "0", "--steps", "40", "--sample", "20"]),
        Case("kepler, adaptive", p("kepler.txt"), "adaptive", ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "40", "--sample", "20"]),
        Case("kepler e=0.9, a first step rejected", p("kepler_e09.txt"), "first rejected",
             ["--dt", "1", "--epsilon", "1e-9", "--steps", "30", "--sample", "10"]),
        Case("kepler e=0.99, rejections that predict", p("kepler_e099.txt"), "rejection",
             ["--dt", "1e-4", "--epsilon", "1e-2", "--steps", "40", "--sample", "20", "--max-iter", "60"]),
        Case("pythagorean, three bodies", p("pythagorean.txt"), "adaptive", ["--dt", "0.01", "--epsilon", "1e-9", "--steps", "40", "--sample", "20"]),
        Case("outer, six bodies, pluto massless", p("outer.txt"), "fixed", ["--dt", "40", "--epsilon", "0", "--steps", "20", "--sample", "10"]),
        Case("kepler spread x5, steps that differ and a mixed attempt", spread, "systems",
             ["--dt", "0.6", "--epsilon", "1e-9", "--steps", "60", "--sample", "20", "--max-iter", "60"]),
        Case("outer x3, adaptive", outer3, "systems", ["--dt", "40", "--epsilon", "1e-9", "--steps", "12", "--sample", "6", "--max-iter", "60"]),
        Case("outer x3, fixed: systems that share a step", outer3, "shared", ["--dt", "40", "--epsilon", "0", "--steps", "12", "--sample", "6"]),
    ]
    for c in cases:
        c.mixed = "mixed" in c.name
    used = {os.path.basename(c.prob) for c in cases} | {"kepler.txt", "outer.txt"}      # the two the ensembles are built from
    unseen = sorted(set(os.listdir(probs)) - used)
    verdict(not unseen, "every problem in data/problems is a case here%s" % (": not %s" % ", ".join(unseen) if unseen else ""))

    with concurrent.futures.ThreadPoolExecutor(max_workers=WORKERS) as pool:
        jobs = [(fmt, engine, arith, c, pool.submit(check_case, cft, fmt, engine, arith, c, env_name, want_default, work, progs))
                for fmt in formats for engine, arith in ENGINES for c in cases]
        last = None
        for fmt, engine, arith, c, job in jobs:
            if (fmt, engine, arith) != last:
                print("check_scalar at %s, --engine %s --arith %s: default, forced off and forced on" % (fmt, engine, arith))
                last = (fmt, engine, arith)
            bad, said = job.result()
            if not bad:
                verdict(True, "%s: the same record, schedule and counters on all three routes; %s" % (c.name, said))
            else:
                verdict(False, "%s: %d of the checks failed%s:%s" % (c.name, len(bad), " (%s)" % said if said else "",
                                                                     "".join("\n         " + b for b in bad[:8])))

    # The route hidden: today's route, and it says why.
    hide = {"CFT_REBOUND_HIDE_FEATURES": ROUTE}
    print("check_scalar: the route hidden")
    for fmt in formats:
        c = cases[1]
        a = Record(cft, fmt, c.prob, c.args, hide, None)
        b = Record(cft, fmt, c.prob, c.args, {env_name: "off"}, None)
        st = a.routes.get(ROUTE)
        verdict(a.rc == 0 and a.recs == b.recs and bool(a.recs) and st is not None and st[:2] == ("off", "hidden") and st[2] == 0 and st[3] > 0
                and a.tail.get("staged_bytes") == b.tail.get("staged_bytes"),
                "%s: %s: the record and staged_bytes of forced off, # route %s=%s (rc %d)"
                % (fmt, c.name, ROUTE, "%s by=%s calls=%d fallbacks=%d" % st if st else "unread", a.rc))

    # 4. The library half on the new route, forced on: the route is opt-in,
    # so the suite's own check_dropin leg runs today's, and this is the one
    # run of the suite that takes the new one.
    if args.library:
        print("check_scalar: the drop-in's suite with the route forced on")
        dropin = os.path.join(args.build, "check_dropin" + exe)
        rc, out, err = run([dropin, "--light"], {env_name: "on"})
        seen = None
        for line in out.splitlines():
            m = re.match(r"^\s+(ok  |FAIL) route scalar=(on|off) by=(\w+) calls=(\d+) fallbacks=(\d+)", line)
            if m and seen is None:
                seen = (m.group(1).strip(), m.group(2), m.group(3), int(m.group(4)), int(m.group(5)))
        right = seen is not None and seen[:3] == ("ok", "on", "forced") and seen[3] > 0 and seen[4] == 0
        verdict(rc == 0 and "check_dropin: every case passed" in out and right,
                "check_dropin --light, forced onto the new route: rc %d, %s; %s"
                % (rc, "every case passed" if "check_dropin: every case passed" in out else "NOT every case passed",
                   "%s=%s by=%s calls=%d fallbacks=%d" % ((ROUTE,) + seen[1:]) if seen else ROUTE + " unread"))
        if rc != 0:
            for line in out.splitlines():
                if "FAIL" in line:
                    print("       " + line.strip())
        print("check_scalar: gate_real's run, forced on against forced off, the wide state dumped")
        real = os.path.join(args.build, "gate_real" + exe)
        for fmt in formats:
            flag = "--" + fmt
            common = [real, flag, "--straight", "--routes", "--steps-a", "20", "--steps-b", "20"]
            rc1, d1, e1 = run(common, {env_name: "on"})
            rc0, d0, e0 = run(common, {env_name: "off"})
            r1 = {m[0]: m[1:] for m in GATE_ROUTE.findall(e1)}.get(ROUTE)
            r0 = {m[0]: m[1:] for m in GATE_ROUTE.findall(e0)}.get(ROUTE)
            on_ok = r1 is not None and r1[:2] == ("on", "forced") and int(r1[2]) > 0 and int(r1[3]) == 0
            off_ok = r0 is not None and r0[:2] == ("off", "forced") and int(r0[2]) == 0 and int(r0[3]) > 0
            wide = [l for l in d1.splitlines() if l.startswith("wide ")]
            verdict(rc1 == 0 and rc0 == 0 and d1 == d0 and bool(wide) and on_ok and off_ok,
                    "%s: %d dump lines %s, %s; on calls=%s fallbacks=%s, off calls=%s fallbacks=%s (rc %d/%d)"
                    % (fmt, len(d1.splitlines()), "identical" if d1 == d0 else "DIFFER", wide[0] if wide else "no wide digest",
                       r1[2] if r1 else "?", r1[3] if r1 else "?", r0[2] if r0 else "?", r0[3] if r0 else "?", rc1, rc0))

    print("RESULT: %s" % ("PASS" if allok else "FAIL"))
    return 0 if allok else 1


if __name__ == "__main__":
    sys.exit(main())
