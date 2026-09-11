# Parcel rounds

How work gets split across several agents at once without the pieces
failing to meet — the method, not this project's instance of it. The
worked instance is `PARCELS.md`; the post-mortem that produced most of
these rules is the last section of `ROADMAP.md`.

Written to be copied into other repositories. Where it cites this one,
that is evidence rather than context — every rule below is here because
something broke without it.

**This is the working copy, and it is no longer the canonical one.**
The method was lifted into a standalone repository - `ParcelRound`,
with the project-specific parts turned into templates - so it could be
reused and handed to other people. That leaves one document in two
places, which is exactly what "exactly one file owns each shared fact"
forbids, and saying so is better than pretending the copies will stay
in step. When `ParcelRound` has a published home this file becomes a
pointer to it; until then, where the two differ, the standalone copy
is right.

---

## When this is worth it

Use a parcel round when there are **three or more pieces of work that
are each session-sized, mostly independent, and all land in one
codebase**. The gain is wall-clock; the cost is that you now own a merge
problem and a verification problem that you did not have before.

Do not use it for one task, for exploratory work where the split is not
yet obvious, or for anything whose pieces cannot be tested separately.
Two agents on a two-way split is usually slower than doing it yourself,
because the brief costs more than the work.

**Agent count is not the bottleneck.** The lead's capacity to merge and
verify is. Four parcels is comfortable; eight is a queue.

---

## The one failure mode

Not "an agent does bad work". Agents mostly do the work in front of
them. The failure is that **the work between the parcels belongs to
nobody**, and it is invisible because every parcel's own gate is green.

The instance: four parcels built a feature, every gate passed, and the
feature did not work. The archive gates drove a stand-in integrator (a
stub whose header said to re-point them once the real one landed; it
landed and nobody did), and the other parcel's gate never wrote a file.
So the whole path had never been executed, and the load routine reported
"restored exactly" for a restore that had not happened, because it
validated the *file* rather than what reached memory. Nobody was wrong.
The seam belonged to no parcel.

Everything below is a defence against that.

### Two rules that prevent most of it

**1. Exactly one file owns each shared fact; everyone else includes it.**
When the same fact is stated in two places, the copies drift — not
might, do. In the instance above: the state struct was defined in two
headers, the descriptor lists in two files, the walker over them in two
files, a count written out by hand in four places and a name list in a
fifth. Growing that list broke all six in different ways: a refused
load, three crashed gates, and a passing gate that printed "only 50
of 48".

**2. Derive counts and name lists; never transcribe them.** A number
kept in step by a comment is a number that will drift. If a test needs a
constant the code defines, have it *read* the definition.

There are two grades of this and it is worth being precise about which
you are asking for, because a brief that conflates them is wrong about
the thing it is trying to protect. **Derived**: the consumer reads the
definition, so there is one value and no way to disagree — a test that
greps the header for a `#define`. **Declared once and checked**: the
count is still a literal, but it sits beside the list it counts and a
selftest asserts it against the list, so a disagreement is caught on the
next run rather than at the next release. The second is often all that a
macro-generated list allows. Ask for the first where it is possible and
the second where it is not, and do not call the second the first.

That second rule pays immediately and visibly. When a cap moved from one
header to another in this repo, the test that reads it stopped with
`no #define CFT_MAX_BODIES in src/ias15_cft.c` — naming the file it
could not find it in. A transcribed cap would have gone on silently
testing a stale number. **A test that breaks when you move a definition
is working.**

---

## P0: make the shared thing shared, before you split

Before dispatching anything, list every file each parcel would touch.
**Anything appearing in three or more columns is not a conflict to
manage — it is the first parcel, and it is the lead's.**

The diagnostic is mechanical. In this round, five parcels each removed
one entry from a list that was stated in four places: two functions, a
test, and the README table. Five agents editing four files each is
twenty edits and a guaranteed drift. So P0 turned the list into one
table both callers walk, and **unlocking a capability became deleting
one row**.

Same move on the test file: one 718-line file with a single `main()`
would have been five agents editing three regions. It became a file per
topic plus a header saying how to add one, so a parcel writes its own
file and makes four small edits outside it.

**Aim for a glob, not a list.** Every parcel in this round still had
to add its file to a build variable, a declaration to a header and a
call to a `main()` — four one-line edits in shared files, which meant
four trivial conflicts per merge. Trivial is a P0 win over what it
would otherwise have been, but zero was available: had the build
globbed the topic directory and the entry points been discovered
rather than declared, the parcels would have shared **no** file at
all. When you design the P0 refactor, ask what would make the seam
disappear rather than what would make it small.

Two conditions on P0:

- **It must preserve behaviour**, and you prove that by running the full
  suite before and after. A P0 that also fixes things is a P0 whose
  green suite means nothing.
- **It must land and be pushed before any parcel starts**, so every
  worktree branches from it. A P0 landing mid-round is worse than no P0.

### Make the registry self-enforcing

When P0 turns scattered statements into a table, make the test **walk
that table and fail by name, in both directions**:

- a row no test exercised → fail, naming the row;
- a test naming a row that no longer exists → fail, naming the test.

The first stops a capability being added untested. The second catches a
test left behind when a capability is removed, which is the one a round
of deletions will actually hit.

This earned its place on the first run here: four rows had **never been
exercised by anything**. The table did not create that gap, it revealed
it.

Then prove the check can fail. Add a dummy row, confirm the failure
names it, remove the row.

---

## Anatomy of a brief

Every section below is here because leaving it out cost something.

**Start here** — what to read, in what order, with file paths. An agent
that has read the post-mortem does not repeat it.

**Your job** — one paragraph. If it needs three, the parcel is two
parcels.

**Files you own** — by path, and by *function* where a file is shared.
"`choose_timestep()` in `src/ias15_cft.c`, and nothing else in that
file" is a workable boundary; "the step control" is not.

**Files you must NOT touch, and who owns each.** Naming the owner
matters: a boundary with a reason behind it gets respected, and an
arbitrary one gets litigated in the report.

**The small edits outside your files that ARE expected**, enumerated.
Without this an agent either avoids a necessary edit and delivers
something that does not build, or takes the absence of a rule as
permission and sprawls. Say "one declaration here, one line in the
Makefile, one call in `main()`" and you get exactly that.

**Host build traps, as verbatim commands.** Not "build it" — the exact
invocation, including whatever is non-obvious on this machine. Here that
is a default compiler of the wrong architecture and a build system that
drops certain environment variables unless they are passed as
variables. An agent will rediscover these in forty minutes; you can
spend three lines instead.

**How to make the tree buildable**, if that takes a step. Three
parcels in one round each discovered independently that a directory
the brief told them to read did not exist in a fresh worktree, and
each solved it from scratch. One line in the brief, or one line in
the ledger from whoever hit it first.

**The base commit, and an instruction to verify it.** A worktree does
not necessarily branch from where you think: one parcel found itself
a merge behind the commit its brief named, noticed, and reset. Say
the SHA and say "check you are on it".

**Working rules the repository learned the hard way** — the two or three
environment-specific traps that have actually cost hours. Keep this list
short and true. A long one gets skimmed.

**The negative control, named specifically.** See below.

**Do not** — push, merge, rebase, or commit to the main branch; weaken
an assertion to make something pass; fix anything outside scope.

**Report format**, and always including: *"anything you found that the
brief got wrong."*

That last line is the highest-yield sentence in the whole brief. In this
round the first parcel back reported two brief errors and a design fact
the brief had asserted backwards. Without the invitation, an agent
routes around a wrong brief silently.

### Expect boundary violations, and judge them on disclosure

A seam drawn by function is a hypothesis like any other, and one round
falsified three of them: a parcel told to own the step control found
that the setting it was implementing also changes the *corrector*, in a
function outside its seam; another needed one call in a function its
brief had not mentioned; a third put its plumbing in an existing file
rather than the new one the brief named, because a new object would have
meant a larger edit to the build than the code was to the file.

All three were right, and all three said so unprompted. That is the
standard: **judge a boundary crossing on whether it was disclosed and
reasoned, not on whether it happened.** A parcel that silently stays
inside a wrong boundary ships something half-implemented; a parcel that
crosses it quietly is the thing the ownership list exists to prevent.
Ask for the disclosure explicitly and you get it.

### Verify the constraints you write down

A forbidden-files list carried from another project is noise at best.
Four of the six files named "never edit" in this round's briefs **do not
exist in this repository** — they belong to a sibling repo, and the list
came across unchecked. Harmless here, but the same carelessness in an
*owned*-files list sends an agent to edit the wrong thing.

Check that every path in a brief exists before you send it.

---

## The ledger

A brief is written once, at dispatch, and cannot be updated. Everything
learned while the parcels are running reaches them **never** — it
arrives in a final report, after the sibling who needed it has already
finished. One round of evidence:

- three parcels independently discovered the same setup problem and each
  solved it from scratch;
- one parcel ended its report with a paragraph headed "for P3" — a fact
  about which velocities the step controller reads — which reached P3
  hours after P3 had started and made its own discovery of the same
  thing;
- the lead learned at the second merge that a refusal row was refusing a
  run that would have been correct, and had no way to tell the three
  agents still running.

So: **an append-only ledger the agents read and write while they work.**

**Where it lives matters.** Agents in worktrees cannot see each other's
files — that is the point of a worktree. The ledger has to sit outside
every worktree, at a fixed absolute path each brief states, and be
**ignored by version control** so it never conflicts and never lands in
history. A directory beside the repository is the simple answer.

**One file per author, append only.** `ledger/P1.md`, `ledger/lead.md`,
`ledger/verifier-P2.md`. Not one shared file: concurrent appends to one
file lose writes, and a per-author file makes locking unnecessary by
construction. Everyone reads the whole directory; everyone writes only
their own.

**When to read**, stated as moments rather than "periodically", which
means never:

- once before starting anything;
- again before designing anything that touches a file the brief marked
  as shared or forbidden;
- again before writing the final report.

**What goes in**, with a bar: *would this have changed another parcel's
work, or the lead's?* Three things clear it.

1. **Environment and setup** — what the tree needs before it builds,
   what the base commit actually is, a tool that is not where it looks.
   The cheapest entries and the ones that save the most.
2. **A brief that turned out wrong.** The moment you find it, not in the
   report. A sibling may be acting on the same wrong premise right now.
3. **A finding about shared code** — behaviour in a function another
   parcel owns, an upstream defect, a constraint that will bind someone
   else.

What does **not** go in: progress, plans, anything only your own parcel
cares about. A ledger nobody reads because it is full of status is worse
than none.

**Mark what you measured.** Every entry says which of its claims were
*measured* and which are *believed*, to the same standard as a report.
Agents build on each other's entries, and an unverified claim propagates
faster than a verified one.

**The lead writes to it too**, and this is half the value: it is the
only channel for correcting a brief after dispatch. Re-briefs, merge
findings, "the row you were told to delete turned out to be load
bearing" — all of it goes there, and the read-before-you-design rule is
what makes agents see it.

At the end of the round the lead folds anything durable into the
repository's own records and throws the ledger away. It is scaffolding,
not history.

---

## Gates and negative controls

**A gate that cannot fail is not a gate**, and the ways one dies are
specific enough to check for:

- **Passing against itself.** A mode the port silently ignored would
  have left every case comparing the default against the default. The
  control is to prove the two configurations *diverge* before asserting
  they agree.
- **Proving nothing.** A case that exercised a polynomial re-read
  triggered on step 1, when the polynomial was still all zeros. Re-read
  zeros, got zeros, passed. The control is to prove the thing being
  tested was non-trivial when the test ran.
- **Validating the wrong artifact.** The load routine checked the file
  on disk rather than what reached memory, and returned success for a
  restore that never happened.
- **Measuring an unspecified difference.** One case collapsed to a
  configuration where the answer was `0/0`, and the two sides differed
  only in the *sign of a quiet NaN* — which the standard leaves
  unspecified. Real difference, permitted, and nothing to do with what
  the case was for.

So: **name each parcel's negative control in its brief.** Not "test it
properly" — the specific control. "A do-nothing routine must leave the
run bit-identical to no routine at all, *and* the active case must
differ from the inactive one." Require the parcel to build the control,
run it, confirm it fails, delete the control artifact, and **report both
results**.

A control described but not run is worth nothing, so ask for its output.

---

## The verifier

For anything hard to check by reading, put a **second agent between the
parcel and the merge**, whose job is to disconfirm.

**Inputs:** the brief, the diff, and the parcel's claimed results.

**It must:**

- **re-run the gate itself**, from a clean build, rather than trusting
  pasted output;
- **re-run the negative control**, or build one if the parcel did not;
- **diff the changes against the ownership list** — anything outside it
  is a finding, even if it looks correct;
- look for the specific cheat shapes: an assertion loosened, a tolerance
  introduced where exactness was required, a skipped case, a constant
  transcribed rather than derived, a `TODO` where work was claimed;
- **verify the "nothing else regressed" claim by running the rest**, not
  by reading it.

**It must NOT fix anything.** A verifier that edits is a second author
with none of the first one's context, and you lose the independence that
made it worth running. It reports; the parcel or the lead acts.

**When it is worth the extra agent:**

- the gate is slow, so there is a real chance it was not run;
- the change touches a numerical invariant, a hot path, or anything
  where "looks right" and "is right" come apart;
- the report claims green without pasting output;
- the parcel is the long pole, where a late-discovered problem is most
  expensive.

**Its own failure mode is rubber-stamping.** Defend against it by
requiring the verifier to state *what it actually ran*, command by
command, and by making "found nothing" an acceptable, unpenalised
answer. A verifier under pressure to produce findings invents them, and
that is worse than no verifier.

### What the first run of this taught

**Give it a numbered list to attack, and expect its best work to be
outside the list.** The first verifier run here confirmed all eight
items it was handed and found two genuine divergences nobody had
thought to ask about — the more serious being that a user callback
writing anything other than the field the port reads back was silently
discarded, where the original would have used it. Eight confirmations
were worth the run; the two findings were worth more. So: name what you
suspect, then ask explicitly what *else* is there.

**One technique worth copying.** To prove a change did not touch a build
it was not supposed to, diff the **preprocessed translation unit** at
both commits, not the source. Here the same file compiles two ways and
only one was in scope; the preprocessed output at the two commits was
3610 lines each and differed by a single line, which settles the
question in a way that reading a diff cannot.

**Distinguish "the shipped code is right" from "the gate would catch it
if it weren't."** Four of that run's findings were of the second kind —
correct behaviour with nothing holding it down, including a bit
comparison that could be swapped for a value comparison and still pass
the entire suite while breaking signed zero. Those are worth as much as
defects, because they are the defects of the next change, and a
verifier is the only role positioned to notice them.

---

## What the lead keeps

- **P0.**
- **The seam tests.** They belong to no parcel, which is exactly why
  they get skipped. Each merge should get a test that exercises *two*
  parcels together. The test that would have caught round 1's failure
  here is nine lines and under a second: N steps straight against a save
  and a resume that add to N, compared bit for bit.
- **Files that are the lead's alone** — the README, published docs, CI
  config. Not because agents cannot write prose, but because those files
  state the project's claims and the claims must be one voice.
- **Every merge, and the full suite after each one.**
- **The ledger** — whatever records what was measured and when.

### Two habits

**Never merge on the strength of a report.** Run the suite yourself.
Agents report in good faith and are sometimes wrong about what their own
change did.

**Resolve conflicts from what git says is conflicted**, never from
what the merge output happened to print, and assert that no marker
survives before you stage. The lead in this round resolved the two
files the merge message named, ran `git add -A`, and committed a
third file still holding its conflict — along with eight agent
worktrees as embedded repositories. Both are one command to check.

**Never merge while a suite is running.** This looks safe — the suite
built its binaries at the start, so a source edit cannot reach it — and
it is wrong for any part of the suite that is *interpreted*. The lead
merged a parcel mid-run on exactly that reasoning; the run then reached
a Python checker, read the merged parcel's **new** checker off disk, and
drove it against the **old** compiled gate, producing three failures
that were entirely self-inflicted. A suite is only as prebuilt as its
least-compiled component, and scripts are never prebuilt.

Two corollaries. A failure you cannot immediately attribute deserves a
timestamp check before a diagnosis: comparing the build time of the
binary against the mtime of the script settled this one in a single
command. And never report such a failure as a defect, or as a false
alarm, until a clean re-run says which it was.

**Read the log, not the exit code.** A background wrapper in this
session reported exit 0 for a run whose log said `Error 1`, because a
trailing `echo` succeeded. If a claim of green rests on a status code
you did not watch produce, it is not a claim yet.

---

## Sequencing

- **Dispatch the long pole first.** It sets the merge order and it is
  where a late problem hurts most.
- **Hold any parcel that shares a file intimately with another** until
  the first lands. Git will merge two disjoint functions in one file
  without complaint; what it cannot merge is one parcel reshaping the
  thing another is hooking into.
- **Parcels finishing early is fine.** Merging is serial and that is the
  real constraint — and the constraint is not the merge, it is the
  **suite run after it**. Budget that explicitly: at twenty-five
  minutes a run, five parcels is two hours of waiting that nobody
  plans for. Where two parcels share no file, merging both and
  running one suite is a defensible trade if you say you made it;
  a failure then costs a bisection over two candidates, which is
  cheap when each arrived green on its own.
- **Re-brief from what comes back.** The first report usually corrects
  the plan. Fold it into the remaining briefs before sending them.

---

## Observed, one round

Five parcels plus a P0 and two verifiers, on a codebase of about
fifteen thousand lines with an exacting bit-identity gate.

**What the system caught that would otherwise have shipped:**

- **Twelve brief errors**, from all five parcels. Every single parcel
  corrected its brief, including on a point the lead had asserted
  backwards: a mechanism described as an exact addition was in fact a
  replacement, which changed what the feature costs at every format
  above the baseline.
- **Four gates that could not fail**, each found by the parcel that
  wrote the code the gate was meant to hold down: a control that left
  every case passing; a control that was vacuous at two of three
  formats because the observable was a rounded view; a configuration
  where the physics made the difference exactly zero, so the obvious
  test would have passed with nothing implemented; and a bit comparison
  that could be swapped for a value comparison and still pass the whole
  suite.
- **Two divergences and four coverage gaps** from one verifier on the
  long pole, none of them on the list it was given.
- **Two previously unrecorded defects**, one of them upstream.

**What it cost:** every merge conflicted, in the same four places every
time, and every one was trivial — a P0 outcome, not luck. Suite time
dominated: about twenty-five minutes per merge, serially, which is more
than the merges themselves and more than anyone budgets for.

**The single highest-yield line in the whole system** is the one asking
what the brief got wrong. It cost nine words and returned twelve
corrections.

---

## Checklist

Before dispatch:

- [ ] Every file each parcel would touch is listed; anything in 3+
      columns became P0.
- [ ] P0 landed, pushed, behaviour-preserving, full suite green both
      sides.
- [ ] Every shared fact has exactly one owning file.
- [ ] Every count is derived from the definition, not transcribed.
- [ ] Registries walk themselves and fail by name, both directions —
      proven with a dummy entry.
- [ ] Each brief names its owned files, its forbidden files *with
      owners*, and the small outside edits that are expected.
- [ ] Every path named in a brief exists.
- [ ] Each brief names its specific negative control.
- [ ] Each brief asks what the brief got wrong.
- [ ] Host build traps are in the brief as verbatim commands.
- [ ] The brief says how to make the tree buildable, and names the
      base commit with an instruction to verify it.
- [ ] The ledger exists, its path is in every brief, and the read
      moments are stated as moments.

At each merge:

- [ ] The control was run and its output reported.
- [ ] The diff stayed inside the ownership boundary.
- [ ] A seam test exercises this parcel against an already-merged one.
- [ ] Conflicts were enumerated from git, and no marker survived
      into the commit.
- [ ] The full suite was run by the lead, on a tree nothing was
      merged into while it ran, and the log read.
- [ ] Anything the merge taught the lead went into the ledger, for
      the parcels still running.
