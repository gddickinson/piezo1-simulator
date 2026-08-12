# Testing your own central claim: what this project's machinery does

*A methods note for someone building a structural-biology tool around a
prediction they believe in.*

This project set out to predict whether a PIEZO1 variant causes gain or loss of
function from structure. **It failed**, across five pre-registered tests and
five predictor families, and then showed the claim cannot be settled with data
that could exist. `docs/CONCLUSION.md` has the numbers.

That is the context for this note, and the reason to read it. The predictor is
not reusable. The apparatus that established the predictor could not be
validated is, and most of it is not PIEZO1-specific.

**Why a failure is the right advertisement for a method.** A pipeline that only
ever confirmed things would tell you nothing about whether its safeguards work.
Every mechanism below is here because it caught something — usually something I
believed, and twice something I had already written down as a result.

---

## 1. Pre-register, with a decision rule that has more than one clause

Write the hypothesis, the frozen predictor, the inclusion criteria, one primary
endpoint, and the decision rule into a file. Commit it **alone**, before running
the comparison. The design numbers — counts, power, exclusions — are properties
of the design and do not depend on the association, so there is no excuse for
computing them afterwards.

**The clause that earned its place.** Round 41's decision rule required three
things: p < 0.05, the effect in the hypothesised direction, *and* its bootstrap
interval excluding zero. The result was p = 0.0477 with an interval spanning
zero. A rule written as "p < 0.05" alone would have produced this project's only
positive result, on its weakest evidence, and it would have been reported.

**State the ceiling before the hypothesis, not after.** Round 48 tested whether
wild-type structural features separate the classes. A feature computed on the
wild-type structure has *exactly zero* within-position variance — the same value
for every substitution at a residue — so it could never assign a direction to a
variant however significant it came out. Putting that in §2 rather than the
closing caveats meant a positive could not have been over-read.

*Implementation:* `docs/NEGATIVE_RESULT_PROTOCOL.md`,
`docs/PREREGISTRATION_ROUND*.md`.

---

## 2. Put a negative control in every test

Choose a predictor that *should* be meaningless and run it alongside the real
one, pre-registered as a control.

In Round 41 the control — raw per-residue variant counts, expected to be shot
noise — gave δ = −0.231 against the real predictor's −0.269. In Round 48 the
control (distance from the symmetry axis) had a **larger** effect than every
mechanistic endpoint. Neither result is visible from a p-value. Both are
stronger grounds for disbelief than any interval, because they show the spread
across endpoints is noise at that sample size.

---

## 3. Cost the design before running another test

When a test returns null the reflex is "we need more data". Check it.

Round 47 asked what sample size the observed effect would need: **134**
directional variants, against **34** available and an optimistic ceiling of
**59**. Round 61 asked the same of the alternative paired design: **8** shared
positions at an implausibly good predictor, against **1**.

That converts "we need more data" into **"the data that could exist is not
enough"** — a different statement, and the only one that is actionable. It says
a sixth test should not be run *whatever predictor goes into it*, which saves
the effort that would otherwise go into a seventh predictor.

*Implementation:* `analysis/feasibility.py`, `analysis/data_routes.py`.

---

## 4. Calibrate every checking instrument before believing it

This is the rule that caught the most, and it is the one most likely to be
skipped.

A cross-check, audit or probe is a *measuring instrument*. Run it first on a
case whose answer is known independently — an analytic shape, a planted signal,
a true null, an enumerable exact value, a deliberately inert input. The
calibration must be able to fail: if no input makes it say "no", it is not a
calibration.

**Nine times in this project an instrument built to check something was itself
wrong**, and every time it returned a plausible number rather than an error:

- a spheroid fitter that would have reported 89% model error;
- a document checker that could not read the Unicode minus its own documents use;
- two dead-code detectors that would have deleted the command-line interface;
- an audit that missed calibrations named in test names rather than test bodies;
- a graft anchored on a whole flexible arm, reporting 19 Å where the local fit
  is 2.4 Å.

**When a checker disagrees with the pipeline, suspect the checker first.**
Historically it has been wrong more often than the thing it was checking.

**Once, it was the other way round, and that is worth stating too.** Round 81
calibrated a new selectivity measurement on a pore with no charge in it, where
the answer is the two ions' mobility ratio. It came back inverted. Suspecting
the instrument was right as a *first* move and wrong as a conclusion: the
inversion was a sign error in the drift term of the transport solver, five
rounds older than the checker, which had made cations drift *up* the potential
gradient since the solver was written. It had survived two independent
cross-checks because both compared magnitudes, and every current the project
had computed was between identical baths — where reversing the field only
reverses the sign, and the sign was then discarded. The rule that saved it is
the second half of the calibration rule rather than the first: the known answer
has to be one the instrument can get *wrong*, and "which way does a cation go"
is such an answer where "how much current flows" is not.

*Implementation:* `tests/test_calibration.py` holds a register mapping every
public checking callable to the test that calibrates it, and fails if one is
added without.

---

## 5. Make every number traceable, and measure the trace

Registering parameters is not enough. Round 49 measured which registered
parameters the code actually *reads* and found **26 of 101 were read by
nothing** — they appeared in the UI with units and citations, overrides on them
were recorded, reports carried a non-default banner because of them, and the
numbers did not move.

That is worse than an unregistered constant, which is at least honestly
invisible. A declaration is not a wire, and only running the code shows the
difference.

*Implementation:* `analysis/provenance_chain.py` records which parameters and
files a computation touches while it runs; `analysis/parameter_effect.py`
overrides each one and checks the answer moves *and* that reset restores it.

---

## 6. Guard the prose, because prose is what goes stale

Twice this project shipped documentation that had quietly stopped being true:
the guided tour still said the claim had been "tested twice" after five tests,
and a progress review counted 40 usable variant positions where there is one.
Both survived because each surface was written by hand at a different time.

Link them by a test. Extract every number from the summary documents and require
each to come from a registry, a frozen record, or an allowlist whose entries
each state why they are exempt — and calibrate that guard too, so an invented
number fails it.

*Implementation:* `analysis/claims.py`, `tests/test_conclusion.py`.

---

## 7. Record a decision not to test

When the feasibility work says a test should not be run, write that down with
its numbers, in the same shape as a result. Otherwise the next person to look at
the data simply runs the comparison.

Add a ratchet: a test that fails if the constraint that forced the refusal ever
lifts. The question then reopens by itself rather than depending on someone
finding the document.

*Implementation:* `docs/NOT_PREREGISTERED_ROUND64.md` and its test, which
ratchets the count of usable positions.

---

## 8. What this cost, and what it bought

Roughly 65 development rounds, 103 registered parameters, 21 guarded claims, a
suite of about 930 tests, and five pre-registered tests.

It bought a negative result that is worth stating: PIEZO1 gain- versus
loss-of-function is not predictable from these structural descriptors, and the
data that would decide it does not exist. Without the apparatus that would have
been "our predictor did not reach significance" — an outcome indistinguishable
from a bug, an underpowered design, or a predictor that was never given a fair
test.

The apparatus is what makes a null informative. That is the case for building
it before you need it, rather than after a result you did not like.

## A null needs a positive control

Every mechanism above is aimed at not over-claiming: calibrate the instrument
before believing it, suspect the checker first, report a null as a null. They
work. Between them they caught an electrostatics constant 10¹⁰ too large, a
spheroid fitter that would have reported 89% model error, a helix detector that
passed 41% of a random walk, and a numbering test that read one protein as
another.

Not one of them fires at a negative result.

A negative from a broken instrument passes every guard in this repository,
because a guard that asks *are you sure?* only ever fires at confidence. That is
not hypothetical here. One round reported that a conduction pathway "does not
separate open from closed", recorded it as an honest null, and **pinned it with
a test**. It was an artefact of how that round composed two criteria: the
discriminating one was a sum over a truncated profile and had collapsed to zero
on every input. Nobody would have looked again, because nothing in the project
looks at a "no".

    A null needs a positive control — an input on which the instrument MUST
    return "yes" — exactly as a checking instrument needs an input on which it
    must say "no". A null from an instrument that has never returned a positive
    is not a result; it is an untested instrument.

Pre-registered statistical nulls tend to satisfy this already without naming it:
a power curve states what effect the test *would* have detected, which is
precisely a positive control. Engineering negatives — "this pathway does not
separate", "this signal is not there", "no entry passes" — usually do not, and
they are exactly the ones nobody re-examines.

The cheap version is a planted signal. Put a known answer into the instrument
and require it to come back out, in the same run that reports the null. If it
cannot, the null says nothing about the world.

