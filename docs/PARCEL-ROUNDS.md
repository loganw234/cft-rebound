# Parcel rounds

How work gets split across several agents at once without the pieces
failing to meet — the method, not this project's instance of it. The
worked instance is `PARCELS.md`; the post-mortem that produced most of
these rules is the last section of `ROADMAP.md`.

Written to be copied into other repositories. Where it cites this one,
that is evidence rather than context — every rule below is here because
something broke without it.

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

### Verify the constraints you write down

A forbidden-files list carried from another project is noise at best.
Four of the six files named "never edit" in this round's briefs **do not
exist in this repository** — they belong to a sibling repo, and the list
came across unchecked. Harmless here, but the same carelessness in an
*owned*-files list sends an agent to edit the wrong thing.

Check that every path in a brief exists before you send it.

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
  real constraint.
- **Re-brief from what comes back.** The first report usually corrects
  the plan. Fold it into the remaining briefs before sending them.

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

At each merge:

- [ ] The control was run and its output reported.
- [ ] The diff stayed inside the ownership boundary.
- [ ] A seam test exercises this parcel against an already-merged one.
- [ ] The full suite was run by the lead, and the log read.
