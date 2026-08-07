# SESSION LOG

Running record of what was done and — more importantly — *why*. Newest first.

---

## Round 42 — the premise was wrong, and that is the result

### What the round assumed
"Two sources make this cheap." The idea was sound: comparing our geometrically
found lipid sites against occupancies from simulations other people have already
run would be an independent check, and an independent method agreeing is worth
more than a better version of ours.

The data is not there.

### What was measured
**MemProtMD holds 1 of 21** catalogued PIEZO entries — only 3JAC, the 2015
structure. Absent are 7WLT, 7WLU, 6B3R and 8YEZ: every structure this project
actually uses.

That absence is only evidence because the probe has a **control**. 2RH1 and
1M0L return 200 on the identical request, so "PIEZO is not there" is
distinguishable from "the URL is wrong". This project has been caught by that
class of error more than once — a PMID resolving to an unrelated paper, a frame
mismatch reporting a closed channel as conducting — and an absence without a
control is not a finding.

**And the one entry cannot answer the question.** 3JAC resolves 918 of 2,547
residues (36%), and of the 15 curated lipid-associated residues it resolves 4:
the PIP2 cluster in full, and none of the three blade basic clusters. A
simulation of a model that omits the lipid-binding residues cannot report their
contacts.

Zenodo's PIEZO1 records turned out to be microscopy TIFFs and PDFs rather than
trajectories. GPCRmd is GPCR-specific and PIEZO1 is not a GPCR.

### Where I stopped
MemProtMD's site is a single-page app and its analysis is browsable but not
fetchable: ten candidate API paths returned 404 against page URLs that return
200. I stopped guessing endpoints at that point rather than continuing, because
the coverage result already answered the round — the API being inaccessible is a
second obstacle behind a first one that is decisive on its own.

### What was built instead
`analysis/external_md.py` implements the **check**, not the comparison. The
conclusion is therefore reproducible, and it will change by itself: if MemProtMD
ever ingests a modern PIEZO structure the coverage function says so and the test
guarding this null fails.

One test exists specifically for the way this round could have fooled itself: an
offline run must return "not checked", never an empty coverage that reads as a
measured absence. A network failure would otherwise manufacture exactly the
conclusion the round reached.

### The pattern
This is the fourth round in a row where the answer was about data rather than
method — Round 34 (one informative variant structure), Round 36 (34 variants
against the 130 needed), Round 41 (an unconstrained gene), and now this. The
project's models are not the limiting factor and have not been for some time.

674 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 41 — the fourth null, and the clause that caught it

### The gene-level answer arrived before the test
gnomAD's constraint metrics for PIEZO1 came back in the first query, and they
are discouraging: **LOEUF 1.10**, **pLI 1.5e-100**, `oe_mis` **1.45**, `mis_z`
**−11.3**. The gene is not merely unconstrained — missense variation is
*enriched* relative to the mutational model. Loss-of-function tolerance is
biologically unsurprising (E756del is a common African allele), but the missense
enrichment means selection is not removing missense variation from PIEZO1 at
all.

That went into §2 of the pre-registration, before the per-position test, because
it changes what a null means. Looking for regional depletion inside a gene with
no global depletion was always a long shot, and saying so beforehand is the
difference between a predicted null and an excuse.

### The result, and the clause that decided it
Primary Cliff's δ **−0.269**, CI **[−0.595, +0.074]**, **p = 0.0477**, on 18 LoF
versus 24 GoF positions. Medians differ in the hypothesised direction — LoF
2.275 against GoF 2.520 missense per residue, about 10% lower.

**p = 0.0477 is below the conventional threshold.** The pre-registered rule
required three things: p < 0.05, δ negative, *and* the bootstrap interval
excluding zero. The first two pass; the interval contains zero. **Fail to
reject.**

Had the rule been written as "p < 0.05" — which is how most of these are written
— this would have been the project's first positive result, on a p-value 0.0023
under the line, with an effect estimate whose interval comfortably contains no
effect. The conjunction was fixed in advance for precisely this case, and this
is the first time it has been the deciding clause rather than a formality.

### The negative control did more work than the primary
§7.4 pre-registered the raw per-residue count as a control that "should be
dominated by shot noise and show nothing". It gives δ **−0.231, p 0.078** —
statistically indistinguishable from the ±25 smoothed predictor, and from both
the ±10 and ±50 alternatives.

So the windowing that makes this a *regional* constraint estimate contributes
nothing. Whatever weak tendency is there is present in the unsmoothed counts,
and the predictor cannot be told apart from its own negative control. That is a
stronger reason to disbelieve the primary than the interval alone, and it exists
only because the control was pre-registered rather than added afterwards to
explain a result.

The summed allele count ran the *opposite* direction (+0.250): how often a
position is hit says the reverse of how many ways it is hit. Uninformative at
q = 0.716, and reported because the family is reported whole.

### R2456 again
Excluded in writing beforehand as the only position carrying both directions —
R2456H/K/P gain, R2456C loss. A position-level predictor is structurally
incapable of separating them, so the exclusion is both necessary and the
cleanest available statement of this predictor's ceiling. The same residue was
Round 7's diagnostic example for the mechanical predictor. Two different
predictors, the same residue, the same reason.

### Where four nulls leave the project
Four pre-registered tests, four nulls, four different predictors: elastic-network
ΔΔG, FoldX, substitution-aware ΔΔG, and now population constraint. Every
diagnosis has ended in data rather than method.

`docs/SCIENCE.md` §8b now says four rather than three. The honest summary is that
this project has found **no signal that separates gain- from loss-of-function**,
from structure or from population genetics.

I named the biggest methodological weakness for whoever picks this up: observed
counts conflate "few variants seen" with "few expected". gnomAD's own regional
missense constraint model fits expected counts per region and would separate
them — that is a better instrument than the one used here, and it is the fair
next test rather than another window width.

668 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 40 — reproducing a published figure, and getting one of two

### Choosing what to reproduce
Haselwandter & MacKinnon's footprint would largely have tested our own linear
solver against our own parameters. Young et al. 2023 was the better choice
because their **full rate set is published and registered**, so the output can
be checked against papers the model was not built from. That is the difference
between an integration test and a restatement.

### One agreement, one disagreement
**T₅₀ = 2.711 mN/m** against Lewis et al. 2015's measured **2.7 ± 0.1** — a
0.4% difference. Young's rates, this project's solver, a third group's
measurement: three independent things landing on the same number. That is the
strongest cross-source agreement in the project.

**τ_inact = 73.3 ms** against Bae et al. 2013's **8.6 ± 0.4 ms** — **8.5×**
slower. Checked that this is not a fitting artefact: the decay is cleanly
mono-exponential, and a bi-exponential fit sends its second time constant to
2.2e6 ms, i.e. a constant. k₂ carries the timescale — at 8 s⁻¹ it sets ~125 ms
before the rest of the system pulls it to 73 — and reaching 8.6 ms needs a
12.8-fold increase to ~103 s⁻¹.

The agreement is what makes the disagreement informative. If both had failed the
natural conclusion would be that our solver is wrong.

### It justifies a policy that already existed
`kinetics.wt_tau_ms` already carried the note that mutants are calibrated by
**fold change** against the wild-type τ, "never by absolute τ across
preparations". That was written as caution, before anyone measured how far apart
the preparations are. This round supplies the number: 8.5×. The policy was right
and is now quantified rather than merely prudent.

### An API I misread, and why it matters here specifically
`with_modification` takes **fold changes, not absolute rates**. Passing
``k2=8.0`` to a model whose k2 is already 8 gives **64**. I read it as a setter,
measured τ = 13 ms, and briefly had two irreconcilable numbers for the same
model — 13 ms and 73 ms.

What makes this worth recording is that the misreading was wrong by roughly the
same factor the round was measuring. Had I not already had the 73 ms figure from
the default model, 13 ms would have looked like a modest disagreement with Bae
rather than a large one, and the round's conclusion would have flipped. The
reconciliation came from diffing the two model objects rather than from
re-reading the docstring, which is the faster route when two numbers disagree
and both look plausible.

`calibrate_k2_for_tau` came out well: it **refuses** an out-of-reach target and
states the reachable range instead of returning its search bound. Clipping would
have looked like an answer.

### Notes
- The figure is `docs/img/young2023_response.png`: tension steps, the T₅₀ curve
  with Lewis's band, and τ against tension on a log axis with Bae's band, so
  both the agreement and the 8.5× gap are visible at a glance.
- No parameters changed; this round consumed the registry rather than adding to
  it, which is what an integration test should do.

660 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 39 — putting the record where the user meets the claim

### The problem
The mechanical ΔΔG is this project's central claim and it has failed three
pre-registered tests. All of that lived in `docs/` and the CLI. A GUI user could
select a variant, see a classification and a score, and have no way of learning
that the score has never predicted anything.

That is the worst possible arrangement: the machinery is visible and the
evidence about the machinery is not.

### What was built
`analysis/prediction_record.py` holds the record as **data, Qt-free**, so the
GUI, the CLI and the tests read the same numbers and cannot drift. Following
`tour.py`, which solved the same problem for the guided tour.

Three things now reach the user:

**Evidence level, beside every variant.** A classification alone reads as a
fact. For 20 of the 46 directional variants the direction is *inferred* from
which disease the variant causes rather than measured, and the panel now says
which, colour-coded, with the sentence explaining the difference.

**The conflicting variant, flagged where it is met.** V598M is the one variant
where the curated and ClinVar-inferred directions disagree. Round 36 excluded it
in writing beforehand; the panel now shows "sources disagree: curated says GoF,
the disease mechanism implies LoF — this project reports the disagreement rather
than resolving it".

**The record itself**, under *Analysis → Variant prediction record…*: all three
tests, the power statement, and five standing caveats.

### The caveat that took the most thought
Four of the five caveats say what the score cannot do. The fifth says what it
still can: *use it to ask which residues sit in mechanically coupled positions;
do not use it to assign a direction to a variant.*

That distinction is real and worth stating. Round 7's diagnostic was that the
predictor reports **where a residue sits** rather than which substitution
occurred — which is a failure for direction prediction and is exactly what a
coupling map should do. Leaving it out would have made the interface
discouraging rather than accurate, and a caveat nobody believes is a caveat
nobody reads.

### Keeping it honest
`verify_record()` re-reads the stored Round 36 run and fails if the frozen
constants drift; measured agreement is 1.2e-4. Same discipline `analysis.claims`
applies to the documented numbers — a frozen record that drifts from its own run
is worse than no record, because it looks like provenance.

The GUI-reachability guard from Round 34 caught the new analysis the moment it
entered the registry without a menu entry. Second time it has paid for itself,
and both times within a round of being written.

### Block L appended
Five rounds since the last review. The pattern worth naming: **three times in
five rounds an alternative built to check the pipeline was itself defective**,
and each time it produced a plausible number rather than an error — a 15–16 nm
tag distance, an 89% dome model error, 32 pA through a closed channel. What
caught all three was calibrating the alternative against a known answer first.
Round 51 proposes making that a written rule rather than a habit.

651 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 38 — the error the intervals never contained

### What the round is for
Round 29 attached three kinds of spread to the headline numbers and stated on
each that **model error is not captured**. Bootstrapping a sphere fit says how
well a sphere is determined; it cannot say whether a sphere was the right shape.
This estimates that missing term where a second defensible model exists.

`ModelError` is deliberately a separate type from `Sensitivity`. Changing a
spring exponent is a knob within one model; changing a sphere into a spheroid is
a different claim about the object, and letting a "sensitivity range" stand in
for "we do not know the shape" would be the wrong kind of tidy.

### The headline: the dome interval is too narrow by six times
A sphere fitted to the transmembrane surface gives **9.454 nm** with a geometric
rmse of 6.180 Å. An oblate spheroid gives flattening **+0.431** and an apex
curvature of **14.991 nm** with rmse **5.243 Å** — it fits *better*, as it must
with an extra parameter, and the surface is plainly not spherical.

Those two radii differ by **5.54 nm (59%)**. The bootstrap interval on the
sphere is **0.92 nm**. So model error is **6.0×** sampling error, and the
interval this project has been quoting measures how well a sphere is determined
rather than whether a sphere was right.

The published 10.2 nm is itself a sphere-based number, so the sphere stays the
right comparator for the literature. What changes is the honest width.

### The fitter was wrong, and would have been reported as science
The first `fit_spheroid` alternated between the centre and the semi-axes with a
hand-rolled gradient step. On the dome it produced an apex curvature of 18.4 nm
and an "89% model error" — a number I nearly wrote down.

Testing it on a *known* spheroid first is what stopped that. Given a full
surface with true (a, c) = (100, 60) it returned (163, 98) and put the centre
89 Å away. Both axes inflated by the same 1.63×, which is the signature of a
drifting centre rather than a shape error — the ratio was right, so any check on
shape alone would have passed it.

Replaced with the exact linear solution: in the axis frame the implicit
quadric is linear in six coefficients, the null vector of the design matrix
gives them up to scale, and completing the square recovers centre and axes with
no iteration. It now recovers known spheroids to **0.01 Å**.

The lesson is the ordering. An alternative model is an instrument, and an
instrument has to be calibrated on something with a known answer before it is
pointed at the thing you do not know. This is the second time in two rounds that
the alternative route turned out to be the faulty one.

### Two smaller results
**Springs, 5.2%.** Cumulative gating overlap is 0.890 / 0.912 / 0.937 across the
three published spring models. All find the transition; the elastic network is
far less model-sensitive than the dome geometry.

**Pore, exactly zero — with a mechanism.** Apollonius and a uniform probe agree
*exactly* at 1.70 Å, because 7WLT's bottleneck lining is carbon and 1.70 Å is
carbon's radius. Moving the probe off 1.70 shifts the answer by precisely the
offset (0.300 Å at both 1.40 and 2.00), which proves the check is live rather
than silently returning one number. So the per-atom refinement buys nothing at a
carbon-lined constriction — a real null, with a reason.

### What is claimed
Every result says "**lower bound**" on its face. Two models disagreeing bounds
model error from below. Two agreeing does not bound it from above, because both
may be wrong in the same direction — and for a dome fitted only over a cap, that
is a live possibility rather than a formality.

640 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 37 — cross-checking four methods, and being wrong about why one disagreed

### What was built
Round 30 re-derived three headline *physics results* by independent routes. This
does the same for four *algorithms*, each swapping out the part most likely to
be wrong:

* **pore radius** — the pipeline runs a coarse polar grid then a shrinking
  pattern search, a local optimiser on a piecewise-smooth surface whose failure
  mode is a local maximum. The alternative maximises the same clearance by 20k
  uniform random probe centres, which has no basin structure to get stuck in.
* **SASA** — Shrake–Rupley with a fixed golden-spiral point set, against
  Monte-Carlo with independent random directions. A defect in the point set is
  invisible to a test that uses the same point set.
* **conservation** — Needleman–Wunsch anchoring, against exact k-mer seeds with
  no dynamic programming and no gap penalties at all.
* **PCA** — the SVD of the centred coordinate matrix, against power iteration on
  XᵀX, which touches no library eigensolver and never forms the 7389×7389
  covariance.

Each alternative is tested first against a case with an analytic answer — an
isolated sphere, a planted component, a slab with a known gap — because two
routes agreeing means nothing if neither is right.

### Results
PCA agrees **exactly** (0.0%, |cos| = 1.000000). SASA to **0.1%**.

The pore check agrees to 5.2%, and the *sign* is the useful part: brute force
finds **0.9783 Å** where the pattern search finds **0.9300 Å**. A brute force can
only match or beat a local optimiser, so the larger value means the pattern
search stops slightly short. That is under-convergence, not a wrong answer, and
a test now bounds the gap so a real regression would widen it.

### The one I got wrong
Conservation correlates at **0.817**, and I wrote in the module docstring —
before measuring — that the residual would concentrate near indels, since that
is where an alignment does work a seed-based anchor cannot. That was a
plausible story and it was wrong. The eight worst-disagreeing positions all have
coverage 1.00: perfectly aligned.

The real cause is a bias in *my* alternative. Anchoring by **maximum exact
matches** is a selection — given a choice of offsets it prefers the one where
residues agree — so it inflates conservation exactly where a position is
variable. The evidence is unambiguous: the k-mer profile's floor is **0.36
rather than 0**, it agrees at invariant positions (0.993 where the pipeline says
1.00), and it reads **0.653** where the pipeline says below 0.50. A test now
isolates the bias on random unrelated sequences, where it reads them as partly
conserved.

So the k-mer route is the weaker instrument, and the check confirms the pipeline
rather than indicting it — the same shape of conclusion Round 30 reached about
the parabola. The docstring was corrected rather than left with the tidier
explanation.

### Why that matters beyond this round
A cross-check whose disagreement is explained by a guess is not a cross-check.
The value came from measuring *where* the two routes part company, which took
about as long as writing the alternative did, and overturned the explanation I
had already committed to prose.

### Notes
- Three cross-check arguments map onto the pipeline's own registered parameters
  (`pore.leash`, `pore.search`, `sasa.probe_radius`); the two k-mer arguments are
  exempted with the reason that they are properties of the alternative
  instrument, not of PIEZO1.
- Six background research agents were restarted on request; all hit an API
  session limit and none completed.

627 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 36 — the third null, pre-registered first

### Why a third test at all
Two nulls already stand. Re-testing the same hypothesis with the same predictor
on the same data is precisely the drift `NEGATIVE_RESULT_PROTOCOL.md` exists to
prevent. What justified a third was that **both sides of the question had
changed, and neither change was chosen after seeing an outcome**: Round 26 took
the predictor's within-position variance from 4.9% to 52.5% (and deliberately
contained no phenotype comparison), and Round 27 grew the directional set from
26 to 46.

The pre-registration was written and committed **alone**, in `af37a82`, before
anything ran. That separation is the whole mechanism: with the document in
history, the endpoints, the decision rule, the exclusions and the caveats cannot
be adjusted to suit what came back.

### The result
Primary Cliff's δ = **−0.249**, CI **[−0.628, +0.151]**, p = **0.405**,
AUROC 0.625, on 19 GoF vs 15 LoF. The direction is as hypothesised; the
significance is not. **Fail to reject.** Nothing in the six-test secondary family
survives Benjamini–Hochberg (min q = 0.591).

Per the protocol, the predictor is **not adjusted and re-run**, and the record
is written in the prescribed order in `docs/VALIDATION_ROUND36.md`.

### What is genuinely interesting, and what it is not
The substitution-aware predictor gives δ = −0.249 where the volume-only
predictor Round 7 used gives −0.025 on the *same 34 variants*. Tenfold larger,
in the predicted direction. Round 26's improvement is real and shows up here.

And across the three tests the effect has grown monotonically: −0.083, −0.211,
−0.249. It is very tempting to read that as a signal emerging as the method
improves. It is **not evidence**, and the write-up says so in those words. Three
point estimates drifting the right way across designs with 13%, low and 50%
power is exactly what either a real medium effect *or* chance would produce.
Separating them needs about 130 variants; there are 34.

That is the honest reading, and writing it down was the point of having
pre-committed the power section.

### Two defects, both found by doing the round properly
**`design.sample_size_for` returned `max_n` for any positive effect size.**
`power_curve` defaults to `alternative="less"`, so a positive δ injected the
effect *against* the alternative, gave ~zero power at every size, and the
bisection walked to its ceiling — which the docstring described as "not reached
within the search range", i.e. as a statement about the design. Anyone asking
"how many variants would a large effect need?" got a confident wrong answer
dressed as a finding. Fixed to use the magnitude and take the sign from
`alternative`; it now reproduces the protocol's own recorded table.

That surfaced while *writing* a power section rather than reading one, which is
the argument for the protocol requiring the section at all.

**The AUROC came back `nan`.** `auroc()` takes scores and a boolean mask; I
passed it two groups, so the second was cast to all-True, there were no
negatives, and it returned nan. My misuse, not a library defect — caught only
because a nan in the output is loud. Had the arrays been ordinary integers it
would have returned a plausible number instead.

### One endpoint could not be run
FoldX ΔΔG is absent from the offline ProtVar cache for **0 of 34** variants, so
a pre-registered secondary endpoint was untestable. It is recorded as such
rather than dropped: §7 forbids removing a test from the family after the fact,
and "could not be run" is a different statement from "was not significant".

### Where this leaves the destination
Three pre-registered tests, three nulls, three different predictors. Round 34
established the structural side cannot supply more data — one informative
variant structure, all gain-of-function. The remaining route that needs no new
experiments is Block K §45: harvesting T50 and inactivation constants from
published supplementary tables, behind the same wild-type gate the curated set
uses.

Until that exists the project should state the position plainly rather than keep
testing, and `docs/SCIENCE.md` §8b now does.

### Notes
- Also fixed a numbering collision I introduced in Round 34: Block J claimed
  Rounds 36–40, which Block I already had. Block J renumbered to 46–50.
- The Round 36 record is pinned by `tests/test_validation_round36.py`, including
  a guard that the write-up still states what the null does and does not exclude.

618 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 35 — the calcium nanodomain, and a prediction that held

### The result
Unusually for this project, the round's prediction survived. The screened
Green's function, fed the two numbers Rounds 31 and 33 already produce, gives
**113.8 µM at the tag** against a **0.2 µM** sensor Kd — **99.82% occupancy**.
A BAPTA-based sensor on a C-terminal HaloTag is therefore saturated whenever its
own channel opens, and reports opening as a binary event.

The roadmap expected ~200 µM at 4–6 nm. The measured value is half that, because
Round 31 moved the tag to 3.95 nm rather than the 4–6 nm assumed when the
prediction was written. Same order, and far above the Kd either way, so the
conclusion the number existed to support is unchanged — which is worth saying
explicitly rather than quietly reporting a different figure.

λ comes out at **148 nm**, far larger than the tag distance, so the exponential
is essentially 1 where it matters: this nanodomain is set by geometry, not by
buffering. That is worth knowing because it means the two buffer parameters —
which are not separately identifiable anyway, only their product enters — barely
influence the answer.

### Making the claim falsifiable rather than merely robust
A prediction that survives every parameter in its own model is only interesting
if the range swept could have broken it, so the sweep spans 80 combinations of
tag distance, calcium fraction and buffering, and a test asserts that **at least
one combination desaturates the sensor**. Two do — and each needs 20 nm *and*
0.5% calcium *and* ≥1 mM buffer simultaneously.

The falsifiers are stated as numbers someone could go and measure: the tag would
have to sit at 373 nm, or calcium carry 4.4e-5 of the current, or free buffer
reach 0.14 M. All are two to three orders of magnitude from reality.

### Two things found on the way
**The sensor has a floor.** Resting calcium at 100 nM against a 0.2 µM Kd
already holds it 33% occupied, so its dynamic range is 33–100%, not 0–100%.
Asking for an occupancy below that floor has no answer at any distance, and the
solver now returns infinity rather than its search bound. This surfaced as two
test failures that looked like numerical noise and were not.

**A silent frame bug.** The report entry detected the C3 axis on the *unframed*
structure and applied it to the *framed* one. The axis and the coordinates were
then in different frames, so the pore was measured along a line that misses the
pore — and the **closed** 8YEZ came back carrying **32 pA** and making a 1.5 mM
nanodomain. Every number downstream stayed finite and plausible; the only tell
was that a closed structure was conducting. Caught by comparing the CLI output
against the standalone calculation, fixed, and pinned by a test asserting a
closed structure never reports its own current.

That is the second time in three rounds that mixing framed and unframed
coordinates has produced a confident wrong answer. Both times the symptom was a
plausible number rather than an error.

### What it joins up with
Round 32 found that a saturating labelling protocol puts a dye on all three
tags, so a 1:2:3 brightness mixture cannot come from a short incubation. Round 35
finds the sensor saturated whenever the channel opens, so brightness cannot
track calcium amplitude either. Together they say published puncta-brightness
heterogeneity points at **unreactive tags and open probability** — two things
that are measurable — rather than at dye kinetics or graded calcium.

### Notes
- Four references added through the title gate: `stern1992`, `naraghi1997`,
  `allbritton1992`, `tsien1980bapta`. 72 resolve.
- 6 parameters registered; the calcium share of the current is `unverified` and
  is the one the answer scales linearly with, so it is swept rather than trusted.
- The GUI-reachability test written last round earned its keep immediately: it
  failed the moment `nanodomain` entered the registry without a menu entry.

607 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 34 — the comparison that was not available, and a GUI audit

### What the round wanted, and why it could not happen
Compare permeation across the four deposited variant structures and read a
direction of change against the measured phenotype. Three things stopped it, and
each is measured rather than asserted:

**Every deposited human PIEZO1 structure is closed.** Bottleneck radii 0.67–0.93
Å against a 1.38 Å cation, so every conductance is exactly zero. A *difference*
in conductance between them does not exist to be measured.

**Three of the four variant entries do not contain their variant.** A1988 is
unmodelled in both entries named for A1988V; E756 is unmodelled in the E756del
entry. Only 8YFG (R2456H) shows its mutation — histidine there, arginine in
every other entry, which is the control that makes it meaningful.

**Three of the entries share one model.** 8ZU3, 8YFC and 9VMX have byte-identical
protein coordinates: 31,839 atoms, the same hash, 0.000 Å RMSD.

That last one needed ruling out before it could be reported. Far the likeliest
explanation was our own fetch cache writing one file under three names, which
would have made it a fact about this project rather than about the PDB. The
files have different sizes, different md5s, and each identifies itself with its
own `data_` block and its own title, so the identity is in the depositions. A
test pins that check alongside the finding, because the finding is worthless
without it.

### What follows
Four deposited variant entries → one resolves its own mutation → one informative,
against 68 curated variants. And all four are gain-of-function: there is no
deposited loss-of-function structure at all, so this route cannot discriminate
direction even in principle. Round 22 found too few phenotyped variants; this
finds too few structures. Both ends are data-limited.

The null is pinned by tests written to **fail if the situation improves** — if a
human structure ever conducts, or a second variant becomes informative, the
suite says so rather than letting a stale limitation stand.

### Out of band: is everything reachable from the GUI?
Asked to check, and the answer was no. `permeation` and `interactions` were in
the shared `ANALYSES` registry and wired into the CLI but absent from every
menu — computable, and invisible to anyone not using a terminal. The pattern is
easy to fall into: an analysis with no picture to draw has no obvious home in a
3-D viewer.

`ui/result_dialog.py` renders any registry result as formatted text, and
`ui/tabular_analyses.py` adds the menu entries. Deliberately generic over the
result dict: a bespoke panel per analysis would look better and would be one
more thing to drift from the function it displays. Each entry carries the caveat
its numbers need *above* them rather than in a docstring.

The guide had stopped at Round 30 and therefore misdescribed the application —
nothing on the tags, the labelling, the ion current, the canonical framing or
multiple structures. Two topics added, including what those models cannot do.

Three tests keep it from decaying: every registry analysis must have a GUI entry
point, every menu action must carry a tooltip, and the guide must mention both
what was added and its limits. The last one caught the hand-wrapped-HTML trap
this project has hit before — a phrase split across a line break fails a naive
substring test for reasons unrelated to content, so the assertion normalises
whitespace.

### Notes
- The audit exemption list grew by two: a zero-initialised counter and a pore
  step already mapped elsewhere. Both mechanical.
- No new parameters: this round measured structures rather than introducing
  physics.

592 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 33 — ion permeation, and where continuum theory runs out

### What was built
`physics/permeation.py` turns the measured pore radius profile into a current:
steady-state drift-diffusion per species, Scharfetter–Gummel discretised so it
stays stable where drift dominates, with Hall access resistance in series at
both mouths — for a pore this short the mouths are a real share of the total
resistance, not a correction.

### The Poisson half does not converge, and that is the finding
Feeding the solved potential back through Poisson diverged immediately: −0.37 V,
then −171 V, then −2×10¹⁶ V. The obvious diagnosis is a missing restoring term,
so I added the one a proper Newton step needs, ∂ρ/∂φ = −F²/(RT)·Σz²c, which is
negative and makes the operator negative-definite. It still did not converge —
the update plateaued and the potential swung ±1.5 V at every damping I tried.

The reason turned out to be physics rather than arithmetic. In 150 mM the Debye
length is **5.7–8.1 Å**; PIEZO1's open bottleneck radius is **3.3 Å**. The
double layers from opposite walls overlap completely, so the pore has no
electroneutral core for a Gummel map to relax onto. That is not a solver defect
to be worked around; it is the statement that a continuum treatment of this pore
is at the edge of its validity.

So the potential is solved in the **electroneutral limit** — current continuity
with the local conductivity — which converges, and agrees with the independent
closed-form series-resistance formula to **1.5%** (41.0 vs 40.4 pS). The Poisson
machinery is kept and `debye_length()` is reported on every result, because the
honest record is that it was attempted and why it was set aside.

### Validation
**41.0 pS for the open 11ZC against a published 25–30 pS** — high by about half.
Closed structures give exactly zero.

The roadmap asked specifically *which* mechanism blocks which, and getting that
right needed a change: my first `_blocked_reason` returned on the first match,
which made 8YEZ and 7WLU look identical when the whole point is that they are
not. `blocking_mechanisms()` now returns every reason. 8YEZ is shut by **two**
(0.95 Å bottleneck *and* a hydrophobic gate at wetting 0.82); 7WLU by **one**
(0.98 Å, wetting 0.11 — no gate).

### Why the agreement is not a prediction
Sweeping the two unmeasured confinement parameters over plausible ranges —
in-pore diffusivity 0.25–1.0 of bulk, ion radius 1.0–2.0 Å — moves the answer
across **16–94 pS**, a 5.8× span straddling the measurement. The model can be
made to agree by choosing values nobody has measured. Both are registered
`unverified`, and a test pins the span so the point cannot be quietly forgotten.

### Two numerical faults, both of which produced a number before an error
The interior rows of the discretised system are built from `D·A/h` with `A` a
few square angstroms, so they sit around 1e-18 while a Dirichlet row is 1.
LAPACK called that singular. It is not singular, it is badly scaled, and since
the interior rows are homogeneous, dividing each by its own largest coefficient
is free. The Bernoulli function also overflowed at large argument, propagating
inf into the matrix as nan; its asymptotes are now explicit. Both live in
`_pnp_kernels.py` — pure linear algebra with no ions in it, testable alone.

### Out of band: seeing the HaloTag fusion
Asked mid-round how to view the tags, the honest answer was that you could not:
Rounds 31 and 32 computed the fusion and the labelling but drew nothing, and the
numbers were reachable only through the CLI. `ui/fusion_controller.py` now draws
them under **View → HaloTag fusion** — tag bodies, linker seams in a colour the
channel never uses, the accessible-volume cloud, and the bound dyes. Everything
is styled so it reads as a model: sphere-of-gyration bodies rather than the
fold, straight seams rather than a conformation, and the envelope shown
precisely so a single sphere is not mistaken for a determined position.

While doing it I found I had **misreported** where the Round 31/32 options live.
I described them as *View*; they had in fact been added to *Options*. Moved to
*View*, which is both where they belong and what the documentation says.

### Notes
- Three references added through the title gate: `coste2010piezo`,
  `gnanasambandam2015`, `hall1975access`. The gate earned its keep again —
  the Coste entry was first rejected because my `expect` word did not appear in
  the resolved title, and a PMID lookup I tried resolved to an unrelated 2025
  paper.
- 12 parameters registered. Several tabulated diffusivities and Shannon radii
  had to move from `physical` to `convention`: the audit is right that a
  `physical` kind promises a citation, and a community-standard tabulated value
  does not have one in this bibliography.
- `report.py` and `permeation.py` both passed 500 lines and were split at real
  seams — the newer analyses into `report_tags.py`, the linear algebra into
  `_pnp_kernels.py`.
- **Not done:** the particle animation driven by the computed current. The
  physics is in place and the current is available to drive it.

580 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 32 — HaloTag labelling: an import, and what it then says

### The criterion was "exactly", so that is what was tested
The kinetics were not derived here. `analysis/labelling.py` brings across the
three equations from `halotag_binding_sim` — exposure, per-site `p(t)`, and
Binomial(3, p) over the trimer — and `compare_with_source()` re-runs the
original functions to check. Over 241 time points the maximum absolute
difference is **0.0** for `p_site`, `p³`, the occupancy distribution and the
sampled histogram, and the Monte-Carlo populations are identical channel for
channel.

That last part drove a design decision. The sampler reproduces the source's two
`rng.random((n_channels, 3))` draws **in the same order**, reactivity first.
A different order gives a statistically identical population and a numerically
different one — which would have satisfied a "close enough" test while hiding a
real divergence behind sampling noise. Reusing a single uniform draw across
snapshots is also what makes an individual channel's dye count monotonic; the
covalent bond does not reverse, and resampling per snapshot would let a
three-dye channel drop back to two while still giving the right marginals.

The equations are vendored rather than imported at runtime, because the sibling
project is not on a fresh clone's path and A5 says a fresh clone must reproduce
the working state. `compare_with_source()` returns `{"available": False}` there
and the test skips; where the source *is* present it must agree exactly.

### What the model then says, which was not what the round expected
At the standard protocol — 200 nM JF646, 30 min, live cell — labelling is
complete in **54 s** to 99%. Per-site p = 1.0000 and the population is **100%
three-dye**. There is no kinetic dye mixture at any realistic concentration;
producing one needs sub-nanomolar ligand or an incubation under a minute.

This bears on why puncta brightness is heterogeneous, where sub-saturation
labelling is a named candidate. The model separates two things that share that
name. The kinetic route is closed at a saturating concentration. The other route
is not: a population of chemically unreactive tags leaves a mixture at **every**
time, because the ceiling is `active_fraction³`. At 90% reactive that is 72.9%
three-dye and 24.3% two-dye, and no incubation removes it. So an observed 1:2:3
mixture under a saturating protocol argues for unreactive tags, not for a short
incubation.

That conclusion rests on two registered-`unverified` numbers — `k_perm_live`
and `active_fraction` — and both say so in their `source_note`. The first is a
transport estimate, the second an assumption; neither is measured, in this
project or the source.

### What the structure adds
`label_sites()` puts the statistics on the three tag centres from Round 31's
fusion model, so occupancy is a set of places rather than a count. The sites are
treated as **equivalent**, and that is a claim worth stating: the trimer is
C3-symmetric, the three C-termini sit at the same height and radius, and the
ligand reaches them from one cytosolic pool. Geometry therefore decides where a
dye is drawn, not whether it binds.

### Notes
- Three references added through `build_references.py`'s title-verification
  gate rather than by remembered PMID: `los2008halotag`, `grimm2015jf`,
  `bertaccini2025piezo1`. 65 references resolve.
- 8 parameters registered. Two of them are `unverified`, and two more had to be
  reclassified from `empirical` to `method` — the audit correctly refused to let
  an unmeasured model estimate sit in a category that promises a citation.
- `parameter_table.py` passed 500 lines; the tag parameters moved to
  `parameter_table_tags.py`, a real seam rather than an arbitrary cut.
- **Not done:** the brightness *animation*. The histogram, the per-site
  occupancy and the figure are there; nothing is rendered on the trimer over
  time. Recorded on the roadmap rather than quietly dropped.

566 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Round 31 — HaloTag fusion geometry, and a ghost structure

### What was built
`structure/fusion.py` places a HaloTag at each of PIEZO1's three C-termini.
There is no structure of the fusion, so the module produces an **accessible
volume** — the region the tag centre can occupy without clashing — rather than a
pose, in the manner used for FRET position modelling. Measured inputs from 6U32
(1.8 Å, TMR ligand bound): Rg **17.6 Å**, N-terminus **19.9 Å** from the centre,
ligand 21.8 Å from that N-terminus. A C-terminal fusion attaches to the tag's
*N*-terminus, so it is that offset, not the radius of gyration, that decides
where the body sits.

### The validation: two of three
**C3 symmetry — pass.** 0.0000 Å on all 20 entries. Exact by construction, since
one envelope is solved and rotated, but measured anyway: a placement that
quietly broke the symmetry would look right in a picture and be wrong in every
calculation after it.

**No steric clash — pass on 18 of 20.** 21.5 Å clearance against the tag's
17.6 Å radius on 8YEZ. 3JAC and 11ZC are marginal at 17.6 and 15.7 Å.

**Tag centre 4–6 nm from the pore exit — misses.** 3.95 nm on 8YEZ, and
3.27–4.21 nm across all 20 structures. The interesting part is that this is not
the unverified linker's fault: sweeping it from 1 to 30 residues — a 30× change
in accessible volume — moves the answer only between 3.0 and 4.0 nm, and moves
it *down*, because a longer tether wraps further around the channel and pulls
the centroid back. So the miss is structural, and it is explained: the 4–6 nm
estimate added the tag's ~2 nm anchor-to-centre offset to the anchor's 2.6 nm
from the pore exit, which assumes the tag points straight away from the channel.
Averaged over accessible directions, many of which run sideways along the
membrane, the mean is pulled in. The band is not unreachable — the envelope
spans 1.7–7.9 nm and 51% of it lies inside 4–6 nm — so the honest statement is
that **the window describes an achievable position, not the ensemble mean**.
Round 35 should take its nanodomain distance from the envelope, not the centroid.

### Two sign faults, both of which gave confident wrong answers
The sweep across all 20 structures is what exposed them; a single structure
would have passed silently.

**The pore exit is not the lowest atom.** In a curved trimer that is a distal
blade tip. Restricting to atoms within a registered `fusion.pore_mouth_radius`
of the axis finds the CTD bundle, which is what an ion actually leaves through.

**`SymmetryAxis.direction` has no fixed sign.** It comes from the rotation
operator relating two protomers, and re-detecting the axis on an already-framed
structure returns −z as readily as +z — it does for 7WLT and 8YFG but not for
8YEZ. Trusting it put those two structures' pore exit at the *extracellular*
end. Together the two faults reported the tag 15–16 nm from the pore exit for
7WLT and 8YFG against 3.9 nm for the same construct on 8YEZ. The cytosolic
direction is now taken from the C-terminal anchors, which are intracellular by
topology.

### The ghost structure
Reported mid-round: loading a new structure left the previous one drawn.
`OverlayController.clear()` finishes by rebuilding the primary view to undo
deviation colouring, and `load_structure` called it *after* clearing the old
view — so the rebuild put the old batches straight back. The bug predates this
round; it was survivable only while deposited frames sat 100 Å apart, and became
obvious as soon as canonical framing made structures superimpose. Fixed by
clearing the overlay first and nulling the view before anything can rebuild it.

**And the feature it suggested.** `ui/companions.py` adds deliberate
multi-structure display: off by default, opt-in under *View → Show multiple
structures at once*, each extra structure drawn in its own colour in the shared
frame, with the Structure panel naming what is on screen. Loading with the
option on demotes the previous primary to a companion, so the feature needs no
second way of opening a file. Turning the option off drops the extras, because a
setting reading "one structure" while three are drawn is worse than either
state. It is deliberately **not** `overlay_controller`, which superposes one
nominated structure and reports an RMSD — that measures; this displays. Every
analysis still runs on the primary, whatever else is drawn.

### Notes
- `fusion.linker_residues` is registered with the `unverified` sentinel. No
  source for the construct states a linker, so it is the one assumed input, and
  the module says to vary it rather than trust it. The sweep above is why the
  result survives that.
- `build_parameters.py` passed 500 lines, so the table moved to
  `scripts/parameter_table.py` and the gate stayed behind. 75 parameters now.
- The GL-dependent companion tests skip on the offscreen Qt platform the rest of
  the UI suite forces, so the ordering fault is *also* pinned by a test that
  needs no GL at all.

550 tests pass, 10 skipped; `parameter_audit` clean; no file over 500 lines;
`screenshot_app.py --structure 8YEZ` completes.

---

## Structure framing — two faults behind one complaint

### What was reported
Reloading 8YEZ after visiting another entry "looked different", and different
structures loaded in orientations that did not overlap.

### Fault 1 — the camera drifted on every load
`_reset_camera` called `Camera.orbit(0.0, -0.42)`. `orbit` **composes** with the
current rotation, which is right for a mouse drag and wrong for restoring a
standard view. Every structure load therefore added another 0.42 rad of pitch,
so the same entry came back at a different angle on each visit. Added
`Camera.set_orientation`, which sets the rotation absolutely, and used it for the
reset. `test_orbit_accumulates_but_set_orientation_does_not` pins both halves —
including a guard that fails if `orbit` ever stops accumulating, so the test
cannot silently stop testing anything.

### Fault 2 — deposited frames are unrelated
Nothing about a PDB frame is canonical. Measured across the 20 downloaded
entries, pairs sit **29–147 Å** apart before any alignment. New
`piezo1/structure/frame.py` offers three modes, exposed under *View → Structure
alignment*: `deposited` (unchanged), `canonical` (the structure's own C3 axis on
+z, cytosolic side at −z, centred) and `reference` (least-squares onto the first
structure loaded).

**Result.** Canonical framing takes those pairs to 0.9–25 Å, in every case
within about an angstrom of what an unconstrained least-squares superposition
achieves — 19.77 vs 19.70 Å for 7WLU on 7WLT, 0.91 vs 0.91 Å for 8ZU8 on 8YEZ.
The residual is real conformational difference, not framing. The one large
value, 6KG7 at 57.8 Å, is PIEZO2 rather than PIEZO1 and its free superposition
is no better at 55.4 Å.

### Three things that had to be got right, each found by measuring

**The z sign straddled the cap.** The rule "the top 10% of residues by number is
the cytosolic end" is wrong: PIEZO's extracellular cap runs to ~2457 of 2547, so
the top tenth mixes cap and CTD and its mean z depends on which is better
resolved. 7WLU and 11ZC loaded **upside down** while reporting a perfect C3 fit
— the failure was silent, which is the dangerous kind. Checked 0.02 / 0.05 /
0.10 / 0.20 against the last 15 residues over all 20 entries: the first two
agree everywhere, 0.10 fails on two, 0.20 on nineteen. `CTERM_FRACTION = 0.02`,
and the rule now *measures where the C-terminus lands* rather than predicting it.

**Chain labels can run either way round the ring.** 8YFG and 8ZU3 both present
chains A, B, D, but numbered in opposite rotational senses. Taking the labels at
face value costs 60 Å of apparent RMSD against 8YEZ — a bookkeeping error that
reads as a conformational change. All six protomer permutations collapse into
two classes of three (cyclic 72.88 Å, reversed 12.45 Å, identical within each),
so `PERMUTATIONS` holds one representative of each and the search is complete.

**Residue number alone is not a correspondence.** A residue number occurs once
*per protomer*, so `dict(zip(res_seq, xyz))` silently keeps whichever chain was
read last and discards two thirds of the atoms — and if the two structures order
their chains differently, the survivors are not equivalent. Correspondence is now
built on protomer blocks over the shared residue basis.

### Why the roll is solved for, not chosen from three
C3 symmetry leaves the roll about z free, so it can only be pinned against a
reference. The starting rule (protomer 0's centroid on +x) is coverage
dependent: 7WLU resolves from residue 576 and 7WLT from 784, and the extra blade
density swings the azimuth by 47.7° — not a multiple of 120°, so no symmetry-
equivalent choice undoes it. Given a reference, the roll is now the closed-form
one-parameter Procrustes optimum. This is a constrained fit and is documented as
one: the axis and the z sign still come from the molecule, only the genuinely
free degree of freedom is set by the reference.

### Consequences elsewhere
- The load path now hands **transformed** coordinates to the physics, so
  `test_measured_geometry_is_unchanged_by_reframing` pins that curvature, depth
  and areas are invariant (they agree to 2e-8 relative). `pore_profile` shifts by
  ≤0.22 Å because it samples on a fixed step whose origin moves with the frame —
  discretisation, not error.
- `ui/model_utils.py` moved to `structure/protomers.py`. `structure/frame.py`
  needed it, and importing from `ui` would have pointed the dependency arrow
  backwards; `analysis/claims.py` was **already** reaching into `ui` for it,
  which quietly made a headless analysis depend on PyQt. `ui/model_utils.py`
  remains as a re-export shim.
- `_open_file` referenced an undefined `rec` on four lines — a guaranteed
  `NameError` for anyone opening a file outside the catalogue. It also skipped
  the per-structure resets (modes, measurements, overlay) that `load_structure`
  does, so stale state carried over. Fixed both.
- `main_window.py` reached 574 lines; the framing, file-open and camera-reset
  code moved to `ui/alignment.py` (`AlignmentMixin`).
- `CTERM_FRACTION` and `MIN_CA_PER_PROTOMER` added to the audit's `EXEMPT` with
  reasons: both are criteria for reading a deposited file, not measured
  quantities with a literature source.

534 tests pass; `parameter_audit` clean; `screenshot_app.py --structure 8YEZ`
completes and now reports the frame in the status line.

---

## Session 1 — 2026-08-05

### Goal
Establish the project: research PIEZO1 thoroughly, acquire all structural and
sequence data, choose the technology stack, and build the foundation layers.

### Environment decisions

**Conda environment `piezo1` (Python 3.11.15).** Created by
`scripts/create_env.sh`. Scientific core from conda-forge (numpy 2.4.6,
scipy 1.17.1, numba 0.66, MDAnalysis 2.10, mdtraj 1.11.1, OpenMM 8.5.2,
pdbfixer, biotite 1.4, Biopython 1.87, scikit-image 0.26, networkx 3.6);
GUI/GL and structural-bioinformatics layer from pip (PyQt6, moderngl 5.12,
PyOpenGL, pyqtgraph 0.14, ProDy 2.6.1, pydssp, freesasa). All 21 imports
verified working on macOS ARM.

**Renderer: moderngl + QOpenGLWidget on PyQt6, OpenGL 4.1 core.**
This was the highest-risk decision, so it was de-risked first with a probe
before any renderer code was written. Result on this machine (Apple M1 Max):
context reports `4.1 Metal - 90.5`, and a fragment shader writing
`gl_FragDepth` compiles. That last point is what matters — it means ray-cast
**impostor** rendering is available, so spheres and cylinders can be drawn as
screen-space quads with per-pixel-exact geometry instead of tessellated meshes.
That is how PyMOL, VMD and ChimeraX get their speed, and it is the difference
between 120k atoms being interactive and being a slideshow.
Rejected: VTK (heavy, harder to style, awkward instancing), pyqtgraph.opengl
(too limited for custom shaders), Qt3D (immature Python bindings), embedding a
web renderer (loses direct access to our numpy arrays).

### Data acquired

- **UniProt**: human Q92508 (2521 aa, verified) and mouse E2JF22 (2547 aa),
  FASTA + full JSON. 38 transmembrane segments annotated per protomer, plus
  topology, PTMs, disulfide C2411–C2415, coiled coil 1339–1368, and 26 natural
  variants with disease annotations.
- **RCSB**: a targeted query (by UniProt accession *and* by entity description,
  rather than free text, which returned 236 mostly irrelevant hits) found
  **28 PIEZO entries**. 21 downloaded as mmCIF. The find that matters most:
  **human PIEZO1 structures now exist** — 8YEZ (3.3 Å apo), 8ZU3 (3.1 Å
  PIEZO1–MDFIC), and, remarkably, **cryo-EM structures of three disease
  variants**: 8YFG (R2456H), 8ZU8/8YFC (A1988V) and 9VMX (E756del). Earlier
  work in this field was almost entirely mouse.
- **Curved/flattened pairs** suitable as morphing endpoints: 7WLT/7WLU (mouse in
  bilayer, 2022) and 11YE/11ZC (mouse in plasma-membrane vesicles, 2026).
- **AlphaFold DB**: `AF-Q92508-F1-model_v6`. Note the v4 URLs 404 — the current
  version is **v6** (created 2025-08-01) and must be discovered via the API
  endpoint, not guessed. This model covers residues 1–2521, which matters
  because no experimental structure resolves the distal blade below ~570.

### Code written

`config.py`, `io/cif_reader.py`, `core/structure.py`,
`structure/superpose.py`, `structure/geometry.py`,
`scripts/build_uniprot_annotations.py`.

**Why a custom mmCIF reader rather than Biopython.** Biopython allocates a
Python object per atom. For a 34k-atom PIEZO1 trimer — and we will load several
at once, plus morph trajectories — that is the wrong shape entirely. The custom
reader walks the file once into contiguous numpy arrays and takes ~0.6 s.
Selections then become boolean masks and coordinate maths becomes vectorised,
which is what both the renderer and the elastic-network model want anyway.

*Bug worth remembering:* the first version of the tokenizer treated only space
and tab as whitespace, so the trailing newline of each `_atom_site` row became a
22nd token and silently shifted every subsequent row by one column. It surfaced
as `invalid literal for int(): 'ATOM'`. Whitespace handling in hand-written
parsers deserves a test — one is planned in `tests/test_cif_reader.py`.

### Scientific validation achieved

The dome-geometry pipeline was checked against the literature rather than
assumed correct. Taking the mid-point of each of the 38 transmembrane helices in
each protomer as a sample of the mid-membrane surface, recovering the three-fold
axis, and fitting a sphere gives:

| Structure | State | Radius of curvature | Dome depth | Excess area |
|---|---|---|---|---|
| 7WLT | curved, bilayer | **9.7 nm** | 4.9 nm | 256 nm² |
| 11YE | curved, PM vesicle | **10.4 nm** | 4.6 nm | 293 nm² |
| 8YEZ | human apo | 12.0 nm | 5.8 nm | 279 nm² |
| 8ZU3 | human + MDFIC | 12.5 nm | 5.4 nm | 270 nm² |
| 7WLU | flattened | 18.4 nm | 2.5 nm | 379 nm² |
| 11ZC | flat, PM vesicle | 21.6 nm | 3.5 nm | 335 nm² |

Published value for the closed state is **10.2 nm** (Haselwandter & MacKinnon
2018, eLife), and ~11.8 nm outside-in (Vaisey & MacKinnon 2026). Our 9.7 and
10.4 nm land squarely there, and the curved→flat contrast is unambiguous. C3
axis recovery is exact (120.00°, 0.00 Å RMSD). **This is now the standing
regression test for the geometry pipeline.**

*Known caveat:* the curved and flattened entries resolve different residue
ranges (7WLT 784–2547, 7WLU 576–2547), so their footprint radii are not directly
comparable. A fair comparison must restrict both to the commonly resolved
residues — to be implemented alongside the morphing code.

### Research

Six parallel literature-research agents were dispatched, writing dossiers to
`ref/research/`. The first wave hit an API session limit and was relaunched
against REST endpoints (Europe PMC, UniProt, RCSB, PubChem) instead of web
search, which had also been exhausted. The lipid-modulation dossier completed
and is a strong result: it supplies the membrane parameters (κ = 20–25 k_BT,
footprint decay λ = 14 nm, ΔA values), ligand potencies, and the important
2026 finding from Vaisey & MacKinnon that mechanical force alone is *not*
sufficient to gate PIEZO1 — a specific lipid cofactor is also required.

### Next

Sequence numbering map (human↔mouse), the annotation layer, then the elastic
network model, then the renderer.

---

## Variants & disease dossier — `ref/research/04_variants_disease.md`

Written 2026-08-05. 68-entry curated JSON variant table (22 GoF, 17 LoF, 8 VUS,
6 blood-group, 15 engineered) plus narrative. Every wild-type residue in the
table was validated programmatically against `ref/sequences/Q92508_human_PIEZO1.fasta`
— zero mismatches.

**The human↔mouse numbering map now exists** (Biopython global alignment,
BLOSUM62, 82.47% identity) and is tabulated in §8.2. The headline result is that
**the offset is not constant**: −6 at E756, −5 through the beam, **+16** across
THU9/anchor/outer-helix, and **+26** from the cap onward. Validated against the
known anchor mouse R2482 ≡ human R2456 (Ma et al. 2018). Traps to remember:
human E2496 ≠ mouse E2496 (mouse E2496 = human **E2470**); human S2446 = mouse
S2472; human E2117 = mouse E2133; human S1335 = mouse S1330.

**Two literature errors found and documented:** the Open Biology 2025 review
gives T2127's mouse equivalent as 2142 (correct: **2143**; the local alignment is
gapless), and Albuisson 2013 gives R1358P as c.4072G>C (ClinVar: **c.4073G>C**).
Also note a WebFetch summarisation hazard — the summariser fabricated a mouse
numbering column for that review; the raw XML had to be parsed directly to catch
it. Prefer parsing tables from XML over summarised fetches for numeric data.

**Blocking finding for the renderer/annotation layer:** the mutated residue is
often *absent* from its own structure. Verified against the local CIFs — all six
human PIEZO1 entries model chain A from residue **570** only, and **E756 is not
modelled in 9VMX**, nor **A1988 in 8ZU8 or 8YFC**. Only R2456 (8YFG) is present.
Residues 2060–2521 are continuously modelled, so the whole pore module renders
cleanly; G253, D669, G718, E756, C1064, K1877 and A1988 need AlphaFold and must
be flagged model-only in the UI.

Population genetics: gnomAD v4 gives pLI ≈ 0 and LOEUF 1.097 — PIEZO1 is
unconstrained for heterozygous LoF, exactly as expected for recessive LMPHM6.
The missense Z of −11.3 is *not* usable: synonymous Z is −14.1, so the mutational
model fails at this locus. ClinVar (2858 records) yields 133 sequence-level P/LP;
truncating alleles are spread uniformly while missense P/LP density is 6–8×
enriched in the pore module (inner helix 5.88 per 100 aa vs blade 0.69).

---

## Session 1, part 2 — renderer, GUI, validation

### The result that matters

An elastic network model built from the **closed** structure alone reproduces
the experimentally observed gating transition. Comparing curved 7WLT with
flattened 7WLU over the 1274 residues common to all six protomers (19.7 Å
trimer RMSD):

| | |
|---|---|
| Best single mode overlap | **0.705** (mode 3, symmetry A, collectivity 0.610) |
| Cumulative overlap over 40 modes | **0.964** |
| Best E-mode overlap | **0.0011** |
| Fraction of overlap² in A modes | **100.00%** |

The symmetry result is the part worth dwelling on. PIEZO1 is a C3 trimer, so
every mode carries an irreducible-representation label. Isotropic membrane
tension is itself C3-symmetric, so only A modes can couple to it at first
order. Every E mode scores essentially zero overlap with the real transition —
the analysis recovers the selection rule without being told about it. That is a
strong internal consistency check, and it means the app can tell a user which
modes are candidate gating coordinates on principle rather than by eye.

### Traps found by checking rather than trusting

**Protomer labels lie.** 7WLT and 7WLU label their three chains in opposite
rotational order around the symmetry axis. Superposing by chain label gave
71.2 Å RMSD instead of 19.7 Å, and the difference vector built from it was
meaningless — the first overlap calculation returned a misleading 0.213 with no
obvious sign of error. `match_protomers()` now always determines correspondence
by superposition, and a test pins it.

**Secondary structure came out 100% coil.** The distance criteria were fine;
the C-alpha pseudo-torsion had an inverted sign, giving −51° for the PIEZO1
inner helix where the IUPAC convention gives +51°, so the helix test never
fired. Fixed by correcting the cross-product order. PIEZO1 now assigns as 77%
helix / 10% strand / 13% coil, with OH, IH and the beam coiled coil all 100%
helix and the cap correctly β-rich.

**Camera framing by bounding sphere** left the molecule filling ~55% of the
viewport. PIEZO1 is a wide flat propeller, so the bounding-sphere radius badly
overestimates what needs to fit. `frame()` now projects into the current camera
orientation and solves for the exact containing distance.

**Disconnected networks.** Writing the ANM tests surfaced that a contact network
in several pieces contributes six rigid-body modes *per piece*. A model with a
detached fragment would have returned rigid-body motions as its lowest
functional modes — silent and thoroughly misleading. `ANM.n_components()` now
detects this and `calc_modes` discards `6 × n_components` by default.

### Research integration

All six literature agents completed on the second attempt (the first wave hit an
API session limit; relaunching them against REST endpoints rather than web
search worked). The most valuable outcome was **independent convergence**: an
agent's alignment of human and mouse PIEZO1 reproduced this project's offset map
block-for-block, and its Yoda1 pocket conversion (human A1718/A2075/A2078) and
selectivity residues (E2117/E2470) matched the values computed here exactly.

Two corrections were applied as a result. The "clasp" and "latch" domains were
**removed** — neither survived verification against primary sources, and
inventing boundaries would colour residues with false confidence. They were
replaced with two verifiable elements: the UniProt-annotated beam coiled coil
and the Piezo1.1 spliced segment. The anchor description was corrected:
P2113/F2114 is primarily the anchor apex brake on the inner helix, with the
cholesterol context secondary.

### Coverage honesty

Building the variant resource revealed something that shapes the whole UI: all
six human PIEZO1 structures model from residue 570 only, and **14 of 68 curated
variants — including the E756del malaria-associated allele — are resolved in no
human structure at all.** Only R2456 appears in its own structure (8YFG). The
viewer now greys those out, states the count, and warns on selection, rather
than highlighting nothing and letting the user assume it worked.

### Renderer

moderngl + QOpenGLWidget, OpenGL 4.1 core, ray-cast impostors. 31 599 atoms and
275k ribbon vertices render in 14–20 ms on an M1 Max. One Qt-specific trap worth
recording: `QOpenGLWidget` does **not** render to framebuffer 0, so moderngl
must be pointed at `defaultFramebufferObject()` every frame or nothing appears.

### Verification

43 tests, ~10 s. They pin the tokenizer's whitespace handling, the reversed-
handedness detection, dome curvature against the published 10.2 nm, the ANM
symmetry characters, the gating-overlap result itself, ten cross-species residue
equivalences (each also checked for matching amino-acid identity), and an
assertion that the numbering offset is *not* constant.

A scripted GUI smoke test (`scripts/screenshot_app.py`) drives the real
application and checks its outputs, so "the app still starts and computes" is a
test rather than a hope.

### Next

Helfrich membrane-footprint solver, tension-dependent Markov gating kinetics
(the Young et al. 2023 PNAS four-state model is fully parameterised in the
research dossier and is the one to implement), conformational morphing between
curved and flat endpoints, and pore-radius profiling.

---

## Session 1, part 3 — morphing, and the repo

`piezo1/structure/morph.py` interpolates between two experimental endpoints.
Three methods, each reporting its own geometric error so the cost of the
interpolation is visible rather than hidden:

| method | worst C-alpha bond error | note |
|---|---|---|
| linear | 2.94 A | the chord artefact |
| restrained | **0.00 A** | distances restored to interpolated targets |
| modal | 1.60 A | 30 modes capture 95.2% of the change |

The chord artefact is worth naming because it is not obvious: under
straight-line interpolation atoms cut chords through space, so C-alpha-C-alpha
distances contract wherever the local motion is rotational. PIEZO1's blades
swing through large arcs, so mid-path frames are measurably wrong. For
comparison, ProDy's `calcAdaptiveANM` was benchmarked on this same 3822-site
trimer at 18 minutes and still left bonds stretched to 5.39 A.

Dome geometry tracks along the restrained path — radius of curvature
9.2 -> 13.1 nm, dome depth 4.6 -> 2.7 nm, excess area 278 -> 316 nm2 — so the
physics follows the morph rather than being imposed on it.

**Refactor.** `main_window.py` reached 657 lines, past the project's 500-line
rule, so dome/mode handling moved to `ui/physics_controller.py` and morph
handling to `ui/morph_controller.py`; the window is now 364 lines. The split
introduced two bugs, both caught by re-running the scripted GUI test rather
than by assuming the move was safe: a guard still testing `hasattr` on an
attribute now initialised to `None`, and a method left with its pre-split name.
Mechanical refactors of GUI code need the smoke test run afterwards, every time.

**Final research agent** returned benchmarks that independently validated the
architecture: sparse ANM at N=7500 takes 4.4 s versus 311 s dense (71x); full
C3 block-diagonalisation would give only 1.76x on top of that, confirming the
decision to do symmetry *labelling* rather than symmetry-adapted solving;
ProDy's `imANM` is an O(N^2) Python loop and unusable at this scale;
MDAnalysis's `hole2` is an empty stub in 2.10 and HOLE has no arm64 build, so a
native pore profiler is required.

**Repository** published to github.com/gddickinson/piezo1-simulator (private),
eight commits, 86 tracked files. `ref/` and `data/` are git-ignored and fully
regenerable with `python -m piezo1.io.fetch`.

### State at end of session

Working: data layer, structure model, cross-species numbering, dome geometry,
elastic network models with C3 symmetry labelling, conformational morphing, the
renderer, and the GUI. 51 tests, ~15 s.

Next, in order of value: the tension-dependent Markov gating model (Young et al.
2023 PNAS four-state, fully parameterised in `ref/research/03a_kinetic_models.md`
— sigma_50 = 1.4, b = 0.8 mN/m, all rates given); the Helfrich membrane
footprint solver (1-D radial validates to 3.8e-3 against the exact K_0 in
0.9 ms, then revolve; lambda = 14 nm with kappa = 20 kT implies gamma =
0.42 mN/m); a native pore-radius profiler (a leash is mandatory — unconstrained
the probe sphere escapes to R = 6188 A); and the hybrid full-length model,
noting that AlphaFold **cannot** place residues 1-570 relative to the core
(PAE 25-29 A against a 31.75 A maximum), so PIEZO2 6KG7 is the better guide.

---

## Session 2 — 2026-08-05/06 — rounds, references, measurement, animation

Set up a ROADMAP.md organised into rounds, and a recurring 20-minute loop that
works the next unchecked one.

### Rounds completed

**Round 1 — pore geometry.** Closed human 8YEZ has a 0.76 Å bottleneck and is
non-conductive; flat 11ZC has 3.25 Å and is conductive. The profiler
rediscovered all three curated constrictions from coordinates alone: the V2450
hydrophobic gate (3.0 Å) and the CTD constrictions at M2467 (1.2 Å) and P2510
(1.4 Å). *The leash is a correctness requirement, not a convenience* — the
clearance function has no interior maximum, and an untethered probe escapes to
R ≈ 6188 Å, which is a true maximum and a useless answer.

**Round 2 — gating kinetics.** Young et al. 2023 four-state tension model.
Emergent half-activation **2.71 mN/m** against a measured cell-attached T50 of
**2.7 ± 0.1 mN/m**. Detailed balance exact to 1e-16, because C→I₁ is fixed by
microscopic reversibility rather than fitted.

**Rounds 11–12 — measurement and interactions.** The C2411–C2415 disulfide is
recovered in all three protomers at 2.04 Å. TM38 is the least tilted
pore-proximal helix at 6.9°. And **R2456 salt-bridges to E2117 of the
neighbouring protomer** in all three copies — the archetypal gain-of-function
residue pairing with the selectivity glutamate, a concrete structural route
from mutation to phenotype.

**Rounds 14–15 — animation.** Seven animations, offscreen-rendered.

### Errors caught, and what they teach

**Mutant calibration was inverted.** Presets originally solved for the rate that
reproduced an *absolute* measured inactivation τ. But Young's parameterisation
gives a wild-type τ of 35–80 ms while Bae's whole-cell measurement gives 8.6 ms
— different preparations. Calibrating R2456H to 22.2 ms therefore made it
*faster* than the model's own wild type, exactly inverting the phenotype of the
best-known gain-of-function variant. Fold changes transfer between
preparations; absolute time constants do not.

**Six citations were confidently wrong.** PMIDs entered from memory resolved
cleanly to unrelated papers — a Piezo1 structure citation came back as a
bone-marrow transplantation study, another as stem-cell reprogramming. Europe
PMC returns whatever the identifier points at. Every seed entry now carries an
`expect` keyword checked against the resolved title, and anything failing is
reported rather than written out. 51/51 now verify.

**Interaction criteria were too loose.** PLIP's 4.1 Å hydrogen-bond cutoff is
only valid *with* hydrogens and an angle test; on heavy atoms alone it produced
8005 "bonds" per trimer including donor–donor N···N pairs. Tightened to 3.5 Å
with N···N excluded except for histidine.

**A test claim was stronger than the data.** "TM38 is the least tilted helix"
is false — blade helices 50–60 Å out are also near-vertical. Narrowed to "least
tilted of the pore module", which is what the coordinates support.

The pattern across all four: the failure mode is never a crash, it is a
confident wrong number. Every one was caught by checking the output against
something independent rather than by the code raising an error.

### Bibliography

51 references resolved from Europe PMC into a committed JSON plus a generated
`docs/REFERENCES.md`, each recording what the project uses it *for*. 29
open-access full texts downloaded to `ref/papers/` (git-ignored — other
people's copyright, and the bibliography suffices to retrieve them).

### Operational note

The 20-minute cron cadence is shorter than a round actually takes, so
invocations stack up. Harmless — each fire simply picks up the next unchecked
item — but worth knowing when reading the history.

---

## Round 3 — membrane mechanics (2026-08-06)

Implemented the Helfrich footprint solver and the dome energetics. Four
published numbers reproduced: γ = 0.420 mN/m from λ = 14 nm and κ = 20 k_BT;
4.116 mN/m per k_BT/nm²; λ = 13.998 nm recovered from the solver's own output;
and T₅₀ = 4.99 mN/m from Cox's ΔG₀ and ΔA against their measured 5.1 ± 0.2.

The footprint result is the one that matters scientifically: around the
measured 7WLT dome it stores **622 nm² of excess area against the dome's own
256 nm²**. Haselwandter & MacKinnon argued the footprint dominates tension
sensitivity; this puts a number on it.

### Two errors that produced plausible numbers

**`L @ L` is not how you build a biharmonic operator.** Squaring the discrete
Laplacian squares its condition number. The solver converged — to a profile
with a 47 nm decay length where the exact answer is 14 nm, and a 59% energy
error that stayed at 59% under grid refinement. That last detail is the tell:
an error that does not shrink when the grid does is not a discretisation error,
it is convergence to the wrong problem. Rewritten as a coupled second-order
system, now second-order convergent.

**The closed-form energy had the Bessel ratio upside down.** K₁/K₀ instead of
K₀/K₁, which is 2.5× too large at PIEZO1's r₀/λ. Both the formula and the
solver looked reasonable in isolation; the disagreement only surfaced when the
functional was integrated over the *exact* analytic profile — a third,
independent route to the same number. Two implementations agreeing is weak
evidence when they share an author; three disagreeing is what localises the
fault.

### One thing flagged rather than fixed

PIEZO1's dome meets the bilayer at a contact slope near 2.0 — about 63°. The
Monge gauge assumes |∇h| ≪ 1 and drops terms of order |∇h|², so at that slope
the neglected terms exceed the retained ones. The code solves it anyway,
because the trend and scale are still informative, but
`FootprintSolution.validity_note()` states plainly that the numbers are not
quantitative and that a nonlinear Helfrich or Euler–elastica treatment is what
the problem actually needs. Reporting a number with a caveat is better than
either silently reporting it or refusing to compute it.

---

## Round 4 — experimental conformational space (2026-08-06)

PCA over the deposited structures, compared with the elastic-network modes.
This is the strongest validation the project has produced.

**PC1 = 90.0% of variance, and it is the gating coordinate.** The PCA sees only
coordinates — no state labels — yet PC1 orders every structure correctly:
seven curved entries negative, the 8IXO intermediate at +334, flattened 7WLU at
+678, flat 11ZC at +1045. It overlaps **0.804** with ANM mode 6, cumulative
0.960 over 30 modes, RWSIP 0.555 against a random control of 0.001.

And the top three principal components all match **A**-symmetric modes, even
though E modes outnumber A two to one in the mode set. The symmetry selection
rule — only C3-symmetric modes can couple to isotropic tension — now shows up
in the deposited structural record, not merely in one pairwise transition.

### Four traps, all of which return a number rather than an error

1. **Species.** Human entries are numbered by Q92508 and mouse by E2JF22, with
   a non-constant offset. Everything converts to human numbering first.
2. **Coverage.** All 20 usable PIEZO1 entries share only 325 residues, because
   a couple of poorly-ordered structures drag the intersection down. Entries
   are dropped worst-first with the cost recorded.
3. **Protomer correspondence.** Four entries label their protomers in reversed
   rotational order. A test now asserts that at least one is detected — if that
   ever returns none, matching has silently broken.
4. **Paralogues.** 6KG7 is *PIEZO2*. Putting a 40%-identity paralogue into an
   ensemble meant to describe one protein's motion would be a category error,
   so it is excluded by default.

### The exclusion that changed the answer

6LQI is the Piezo1.1 splice isoform, missing residues 1382–1405. What sets it
apart from the rest is a **sequence** difference, not a conformational one, but
PCA cannot tell those apart — it sees only coordinate variance. Included, it
dominates an entire component by itself and splits the gating coordinate across
PC1 (58%) and PC2 (36%). Excluded, PC1 is a single clean 90%.

The general lesson: an ensemble method will happily report the largest axis of
variation without caring whether that variation is the biology you were asking
about. Deciding what belongs in the ensemble is part of the analysis, not a
preliminary to it, and every exclusion here carries its reason in the code.

---

## Round 5 — allostery and force transmission (2026-08-06)

Completes Block A. Perturbation response scanning, dynamic cross-correlation
and correlation-weighted allosteric pathways, all from the elastic network's
covariance — which is never formed in full, since for a PIEZO1 trimer it would
be an 11466² matrix of about a gigabyte.

**The anchor is the transmission hub.** Forcing the blade→gate path through it
costs a detour penalty of −0.000: it is already on the optimal route. It ranks
second by betweenness (5.19) behind only the CTD (7.67), and the cap is clearly
not a transmission route (+0.055).

**The beam result is softer than the prediction, and is reported that way.**
The lever model says the beam carries blade motion to the pore. It does not
appear on the single shortest path; but forcing the route through it costs only
+0.010, so it is a near-degenerate parallel channel rather than an excluded
one, with low but real betweenness (1.30). The honest statement is "viable
parallel route, not the dominant one" — not "confirms the lever model" and not
"refutes it".

Worth noting the beam is fully resolved (all 66 residues) in 8YEZ, so this is a
genuine negative rather than an artefact of missing density. That was checked
before drawing the conclusion.

### The error: a detour cheaper than the direct path

Asking whether the signal passes through a region by computing source→region
and region→target separately and adding their costs is wrong. Each leg
independently picks its best endpoints, and on a C3 trimer those can be in
*different protomers* — so the two legs never join into a path at all. Done
that way, routing "via the beam" came out at 0.101 against a direct path of
0.223: cheaper than the shortest path, which is impossible by definition.

The tell was the impossibility itself. A constrained optimum can never beat an
unconstrained one over the same feasible set, so a negative penalty is not a
surprising result, it is proof of a bug. `detour_cost()` now minimises
`d(source→v) + d(v→target)` over a shared via-point `v`, and an invariant test
asserts the constrained path is never cheaper.

Single shortest paths also turned out to be fragile — one marginally better
edge reroutes the whole thing — so `path_betweenness()` aggregates over many
source/target pairs, which is the standard dynamical-network-analysis answer
and far more stable.

### Block A review

The physics chain is closed and every link is validated against a published
number. The two results strong enough to build on: PC1 of the experimental
ensemble *is* the gating coordinate and matches an A-symmetric mode at 0.804,
and the anchor is the dominant force-transmission hub. Block E has been added
to the roadmap: allostery-derived per-variant features, licence-clean external
predictors via the ProtVar API, nonlinear membrane mechanics to fix the flagged
small-slope violation, pore wetting prediction, and — written *before* the
blind test rather than after — a statistical protocol under which "the
predictor does not separate GoF from LoF" is a recordable outcome rather than a
prompt to keep tuning.

---

## Rounds 6-7 — the blind test, and a null result (2026-08-06)

Round 6 built the predictor: ΔΔG_gating = ½dᵀ(H_mut − H_wt)d, the change in
elastic cost of the *observed* gating motion. Exact to 7e-16 against an
explicitly rebuilt Hessian, and cheap, because H_mut − H_wt is non-zero only at
the mutated residue's contacts.

**No phenotype comparison was made in Round 6, deliberately**, and
`docs/PREREGISTRATION.md` was written before Round 7 ran — pulling Round 20's
protocol forward, because a blind test is only blind if the rule was fixed
first.

### The result

**H0 not rejected.** 25 variants (16 GoF, 9 LoF) after pre-registered
inclusion. Permutation p = 0.234, Cliff's delta −0.083 with CI spanning zero,
AUROC 0.542. The mean difference points the predicted way; the effect is
negligible.

Writing the pre-registration turned out to matter more than expected. Two
numbers were close enough to tempt: the secondary analysis on the normalised
score came in at **p = 0.054**, and the primary difference had the right sign.
Without a decision rule fixed in advance, "the normalised version is basically
significant" is exactly the sentence one writes. The rule said p ≥ 0.05 is a
null result, and it is reported as one.

### Why it fails, and why that is the useful part

Post-hoc, labelled as such: **99.8% of the ΔΔG variance is between-position,
0.2% within-position.** The predictor reports *where a residue is*, not *what
happened to it*. That falls straight out of the construction — ΔΔG scales with
local gating strain and contact count, both properties of the position, while
the substitution enters through a single scalar spring multiplier, a far weaker
lever.

The R2456 series makes it concrete: four substitutions at one position, GoF and
LoF among them, all predicted to soften, with the *largest* softening belonging
to the loss-of-function variant R2456C. Positions do not have phenotypes;
variants do.

### What this does not mean

It does not touch the physics chain. Dome curvature 9.7 nm against a published
10.2; T₅₀ 4.99 mN/m against 5.1 ± 0.2; footprint decay 13.998 nm against 14.0;
PC1 of the experimental ensemble matching an A-symmetric mode at 0.804. Those
stand. What fails is a single scalar as a phenotype call.

The right reading is that the elastic network models the *machine* well and the
*substitution* badly. Round 17's sequence-based predictors (AlphaMissense, EVE,
ESM-1b via ProtVar, all substitution-aware by construction) are the natural
complement — they lack the mechanism, which is exactly what this supplies.
`docs/VALIDATION.md` §8 records what a fair next test would look like, written
now so that whatever is tried next is a stated new hypothesis rather than a
retrofit.

---

## Round 8 — pockets and ligands (2026-08-06)

Delaunay alpha spheres, the fpocket construction, reimplemented in numpy. Every
tetrahedron of the Delaunay triangulation has an empty circumsphere; its radius
says whether it sits in packed interior, in bulk solvent, or in a cavity of
small-molecule size.

**Two annotated sites recovered from geometry alone:** the transmembrane
hydrophobic gate (2/3 residues) and the anchor-domain apex brake (2/2). Neither
the detector nor the annotation knows about the other.

### The percolation trap

The first run produced a top "pocket" of **408 000 Å³ with 601 lining
residues** — which is the outside of the protein. On a large open structure a
radius filter alone is not enough: PIEZO1 is a curved propeller with enormous
grooves between its blades, and single-linkage clustering happily percolates
the whole exterior into one object.

Requiring each alpha sphere to have at least 30 atoms within 8 Å discards the
surface spheres and stops the merge; the largest pocket becomes 6 691 Å³ with
63 residues. **The parameters were chosen on pocket-size plausibility — a
druggable cavity is hundreds to a couple of thousand cubic Ångström — and
fixed before any site recovery was checked.** Tuning them until Yoda1 appeared
would have made the recovery meaningless.

### The Yoda1 result, and why the negative is the interesting part

The detector does *not* recover the Yoda1 site as an enclosed cavity: at most
one of its three residues, in either the human apo structure or the mouse
lipid-bound one. Allowing surface grooves recovers two.

Rather than report that as a failure, it is worth asking what it means, and
three independent facts line up. Yoda1 is proposed to act as a **molecular
wedge from the lipid phase**, which is a description of an interfacial site
rather than a pocket. A **PLX lipid occupies part of the site** in 7WLT — the
contact mapping here confirms it touches A2091. And the site has never been
seen in a co-structure: every PIEZO entry in the PDB contains only lipids, so
the mapping rests on mutagenesis and docking, which is exactly why this
project's annotation labels its evidence as *predicted*.

So the honest statement is that the Yoda1 site is interfacial, not enclosed,
and a cavity detector is the wrong instrument for it. That is a more useful
sentence than either "not recovered" or a tuned-up "recovered".

### Also fixed

`np.linalg.solve` in numpy 2 treats a 2-D right-hand side as a single matrix
rather than a batch of vectors, so the batched circumsphere solve needed an
explicit trailing axis. And the pocket tests originally took 137 s because each
one recomputed the pockets; module-scoped fixtures and a smaller Monte-Carlo
sample bring that to 16 s.

---

## Round 9 — conservation and constraint (2026-08-06)

62 vertebrate orthologs from UniProt, deduplicated to one per species,
reference-anchored to human numbering.

**The result worth reporting is a convergence.** Ranked by mean conservation,
the **anchor domain is the most constrained region of PIEZO1 (0.987)**. Round 5
had already identified the anchor as the force-transmission hub — the only
region with a zero detour penalty on the blade-to-gate path. Those are entirely
independent lines of evidence: one is elastic mechanics on a single structure,
the other is 62 genomes' worth of selection. They agree.

The gradient across the whole protein is consistent with that picture: pore
module and anchor at 0.95–0.99, the cap and beam in the middle, and the distal
blade THU1 at 0.719 — least constrained, and also the region no experimental
structure resolves.

### Conservation alone is not a hypothesis

426 positions are invariant, carry no reported variant, and are structurally
resolved. That is about a quarter of the modelled protein, and as a "candidate
functional site" list it is useless — PIEZO1 is simply very conserved.

What makes it a hypothesis is the intersection with mechanics, which is
something this project can do and a sequence-only method cannot. Crossing
conservation with the Round 5 perturbation-response and path-betweenness
profiles narrows it sharply, and the top distal candidates are dominated by the
anchor (20 of 40). Two residues, **2021 and 2034**, are invariant across all 62
species, have never been reported as variants, and lie *on the blade-to-gate
allosteric path* — they appear literally in that path's residue list. Those are
well-motivated mutagenesis targets.

Also worth noting as a negative control that behaved: the **Yoda1 pocket is the
least conserved annotated site** (0.859, with A2075 at 0.63). A synthetic
agonist's binding site has no reason to be under selection, and this is the
third independent observation pointing at that site being unusual — the pocket
detector found it interfacial rather than enclosed, the PDB contains no
agonist co-structure, and now evolution is indifferent to it.

### The bug

The reference sequence was skipped by testing `member.sequence == ref`, which
discards *every* sequence identical to the reference rather than the reference
entry itself. Two closely related species with identical sequences would both
vanish, and on a small set that empties the alignment completely — the unit
test with three identical sequences returned conservation 0 everywhere. Fixed
to skip at most one exact match. Real-data numbers are unchanged, which is the
point: the bug was only reachable with near-identical inputs, which is exactly
what a synthetic test provides and real orthologs do not.

---

## Round 10 — research workflow (2026-08-06) — Block B complete

A headless CLI, provenance-stamped reports, session save/load and a documented
notebook API.

The design decision worth recording is that **the CLI and the report share one
analysis registry**. Adding an analysis to `ANALYSES` makes it available in
both at once, and a test asserts every registered analysis is reachable from
the command line. Two parallel dispatch tables would have drifted apart within
a round or two.

**Sessions store what you were looking at, never the data.** Structure, style,
camera, selection, which analyses had been run and with what parameters — but
no coordinates and no results. A session carrying its own copy of the numbers
would let a saved file drift silently out of step with the code that produced
them, which is precisely backwards for a reproducibility feature.

### What the batch run turned up

`python -m piezo1.cli batch` over all 20 structures reproduces the whole gating
series in one command: curved entries clustered at R_c 9.3–12.5 nm against a
published 10.2, the 8IXO intermediate at 16.5, and flat 11ZC at 21.6 as the
only entry called conductive.

It also **independently flagged 3JAC** — R_c 5.3 nm and spuriously conductive.
That is the same entry the Round 4 ensemble excluded for having poly-UNK
regions with arbitrary residue numbering, found again by a completely different
route. A batch mode earns its place partly by surfacing this kind of thing
without being asked.

### The argparse trap

`--json` was a top-level flag, so `cli dome 8YEZ --json` — which is how anyone
would actually type it — failed with "unrecognized arguments". Adding it to
each subparser via a shared parent fixed that but broke the other order:
a subparser writes its own default over whatever the parent already parsed, so
`cli --json dome 8YEZ` silently came back as `json=False`. `default=SUPPRESS`
plus a default supplied in `main()` makes both orders work, and a test pins
both.

### Block B review

The blind test came back null and its diagnostic was precise: the mechanical
predictor reports *where* a residue sits, not *which* substitution occurred.
Round 9 then showed conservation crossed with mechanics is sharp where either
alone is blunt, and Rounds 5 and 9 converged independently on the anchor.

The clearest gap now is that **the engine has raced ahead of the interface**.
The GUI can show a structure, a dome measurement and normal modes; it cannot
reach the pore profile, pockets, conservation, allostery or any reporting.
Block F leads with that, then a *new* pre-registration for a second variant
hypothesis (the Round 7 result stands as recorded and is not to be revised),
packaging, performance, and the teaching layer — which project aim A1 asks for
and which has had the least attention of anything.

---

## Round 16 — the feature table (2026-08-06)

`analysis/features.py` assembles everything the project computes into one row
per residue: PRS gate response and coupling, path betweenness, correlation to
the gate, amplitude along the lowest symmetric mode, fluctuation, relative
SASA, conservation, geometry and domain. 1279 residues × 11 features for 8YEZ
in 9 s.

**No phenotype comparison was run, deliberately.** The Round 7 blind test
returned a null result which stands as recorded; re-testing these features
against the variant labels needs the new pre-registration scheduled for Round
22. Assembling predictors and evaluating them in the same breath is exactly how
a blind test stops being blind, so the validation here is entirely structural:
gate response falls off with distance to the gate (r = −0.55), relative SASA
lies in [0, 1], the mode used as the gating coordinate is confirmed to be
A-symmetric, and the conservation join reproduces Round 9's domain ranking
exactly — which it would not if the residue join were off by a single position.

### The redundancy the tests caught

I wrote a test asserting that no two feature columns correlate above 0.99, on
the general principle that a near-perfect correlation means a duplicated
column. It fired immediately, and what it exposed is worth stating.

**The PRS response matrix is symmetric.** Its entries are Frobenius norms of
covariance blocks, and the Frobenius norm is invariant under transpose, so
‖C_ij‖ = ‖C_ji‖. Row and column means are therefore the *same numbers*:
"effectiveness" and "sensitivity" — effector and sensor — are not two
populations to compare, they are one quantity written twice.

Row normalisation appears to break the symmetry, and I initially "fixed" the
problem that way. It does not work either. Normalisation forces every row mean
to 1, so effectiveness becomes near-constant (spread/mean 0.0018, a range of
0.9655–0.9999), and the normalised *column* mean still correlates with the raw
row mean at **0.998**. Two columns, one quantity, wearing different scales.

Shipping both would have looked like two independent lines of mechanical
evidence and been one — the sort of thing that quietly inflates a combined
score. The table now carries a single, honestly named `prs_coupling`, and the
symmetry is documented on `PRSResult.is_symmetric` so the next person does not
re-derive it.

Highest remaining inter-feature correlation is 0.953 (`prs_coupling` against
`gating_amplitude`), which is high but genuinely two things: how strongly a
residue couples to the whole protein, and how far it moves along the specific
gating coordinate.

## Round 17 — external predictors, and a payload that lies by omission (2026-08-06)

The mechanical ΔΔG cannot see the substitution (Round 7). The obvious response
is to bring in predictors that can — AlphaMissense, EVE, ESM-1b, FoldX. The
obstacle was never the science, it was the licences.

### Why an API rather than local tools

FoldX is not redistributable at all: academic use needs a signed agreement, and
the community Python wrapper carries no licence file, which means all rights
reserved rather than permissive. SIFT4G is GPL-3.0, so calling it would drag
this whole PyQt application under copyleft. Everything on biosig.lab.uq.edu.au
— mCSM, DynaMut2 and relatives — publishes no licence text whatsoever, so
nothing is granted and using it would be a bet, not a permission. VarSite and
VarMap are both retired.

ProtVar (EMBL-EBI) serves all four predictors plus per-position conservation
from one endpoint, and I confirmed **CC BY 4.0** by reading `info.license` in
the service's own OpenAPI document rather than assuming it from the EBI's usual
terms. One source, one licence, no local models, attribution recorded on every
cached response and in `docs/REFERENCES.md`.

### The trap

`/score/{acc}/{pos}` returns nineteen entries per predictor — one per possible
substitution — and **no field anywhere in the payload says which entry is which
mutation**. Nothing errors. You get a well-formed JSON array of plausible
pathogenicity scores, and if you read them in array order, or assume they are
alphabetical by residue, you attribute the wrong score to every variant in the
study and never find out.

The fix is an undocumented `mt=` query parameter, which returns the single
score for the substitution asked for. A position-only query now keeps *only*
the conservation value, which is genuinely position-level; the missense scores
are dropped rather than guessed. `/prediction/foldx/` needed no such care — it
labels every entry with `mutatedType`, which is how the endpoint should behave
and how I noticed the other one did not. (The documented
`/prediction/interaction/` endpoint 404s; not needed here.)

This is the second time in this project that a silent field-alignment error has
been the real danger — the first was the mmCIF tokenizer shifting every column
by one because a trailing newline was not treated as whitespace. Both produce
confident numbers. Neither raises.

### What came back

64 of 65 single substitutions annotated in 77 s: conservation 64,
AlphaMissense 51, EVE 51, ESM-1b 51, FoldX ΔΔG 50. The 13 without missense
scores are nonsense and frameshift variants, where a missense predictor
correctly has nothing to say — absence there is right behaviour, and there is a
test asserting it stays that way.

### An unplanned validation worth more than the scores

ProtVar reports the wild-type residue it holds at each position, so annotating
the variants cross-checked our numbering against Q92508 from outside the
project. **Zero mismatches out of 64.** Rounds 1–2 established that the human
and mouse numbering differ by a non-constant offset across twelve blocks, and
that work has until now been checked only against itself. This is the first
independent confirmation that every curated variant sits on the residue we say
it does, and it is now a test.

### What these predictors still cannot do

All three missense predictors emit a single *pathogenicity* axis, benign to
damaging. That axis has no dimension in which to express **direction**. The
demonstration is R2456: all four substitutions score PATHOGENIC, but R2456H,
R2456K and R2456P are gain-of-function and R2456C is loss-of-function. A
predictor that calls all four damaging is not wrong — it is answering a
different question from the one this project is aiming at.

So the position after this round is that we hold two families of features that
fail in *opposite* directions: mechanical ones that see the position but not
the substitution, and sequence ones that see the substitution but not the
direction. Whether combining them recovers direction is a real hypothesis and
the only reason to run Round 22.

**No phenotype comparison was run.** Round 7's null result stands as recorded.
Touching the 68 labels again requires the new pre-registration to be written
first — that is the whole point of pre-registering, and the temptation to
"just peek" at 51 fresh AlphaMissense scores is exactly what the discipline is
for. Suite 243 → 254 passing, all of the new tests offline from cache.

## Round 18 — the nonlinear footprint, and a Round 3 result overturned (2026-08-06)

Round 3 built the linearised Helfrich footprint solver, got it right against the
exact Bessel solution, and then *flagged* that PIEZO1's contact slope of ~2.0
(63°) is far outside the regime where that theory holds. Having flagged it, it
quoted the numbers anyway: 622 nm² of footprint excess area against the dome's
256 nm², "about 2.4× as much deformable area as the dome". This round did the
calculation the caveat was standing in for, and the caveat was not strong
enough. The linear number is **3.5× too large**, and the conclusion reverses.

### The formulation

The Monge gauge writes the surface as a height field h(r) and expands in |∇h|.
There is no fixing that at 63°; the terms it drops are larger than the ones it
keeps. So parametrise the meridian by arc length instead, with ψ(s) the tangent
angle, which makes the principal curvatures exact with no expansion at all:
c₁ = ψ̇ and c₂ = sin ψ / r. Minimising the Helfrich energy subject to ṙ = cos ψ
with a multiplier η gives a first-order system in (r, z, ψ, M, η), solved as a
boundary-value problem.

The nice part is free. The Lagrangian has no explicit s dependence, so its
Hamiltonian is conserved, and that Hamiltonian *is* the axial force transmitted
through the membrane — zero for an inclusion nobody is pulling on. Imposing
H = 0 as a boundary condition and then measuring how far H drifts along the
solved profile gives an error estimate that costs nothing and depends on
nothing I derived. It runs at 7e-11.

### Checking a derivation I did by hand

A BVP solver will converge happily onto the wrong equations. So the checks that
count are the ones that do not reuse the derivation:

1. **Small-slope agreement.** The roadmap's own criterion. The relative
   discrepancy divided by slope² converges to a constant, 0.746 — not merely
   shrinking, but shrinking at exactly the order the Monge expansion discards.
2. **The exact functional in another gauge.** Re-evaluate the energy of the
   solved shape in Monge form with the unexpanded expressions
   (dA = √(1+h'²)·2πr·dr, c₁+c₂ = h''/(1+h'²)^{3/2} + h'/(r√(1+h'²))). Agrees
   to 1.3e-3, which is finite-difference noise.
3. **No nearby shape is cheaper.** Perturb the profile with bumps that respect
   both boundary conditions; the energy may only rise.

Check 2 is where I nearly went wrong in an instructive way. My first version of
it minimised the exact functional from scratch, starting from a linear guess,
and it *disagreed with the solver by 142%* at slope 2. That looks like a
refutation. It was not: the minimiser returned a **higher** energy than the BVP,
and a minimiser that has not converged always errs high. Starting it from the
BVP solution, it could not improve on it in 151 iterations. The independent
check was the weaker instrument, and the tell was the sign of the disagreement.

### The result

At the measured 7WLT geometry (inclusion radius 8.69 nm, contact slope 1.99):

| | linear | nonlinear |
|---|---|---|
| footprint energy | 92.2 k_BT | **25.3 k_BT** |
| footprint excess area | 622 nm² | **179 nm²** |

Invariant to domain truncation from 8λ to 40λ at six significant figures, to
grid, to solver tolerance and to whether the slope is walked up by continuation
or hit directly. The correction factor is 3.46–3.67× across κ = 20–25 k_BT and
γ = 0.42–3.0 mN/m, so it is not an artefact of one parameter choice.

### What was actually wrong with Round 3

Two things, and the second is the one worth remembering.

The first is the obvious one: linear theory used outside its range.

The second is that the comparison was never like for like. The dome's 256 nm² is
an **exact** area difference, measured from the fitted spherical cap. The
footprint's 622 nm² was **linearised**. Putting them in a ratio compares a
quantity to an approximation of a different quantity. Measured consistently the
footprint holds 179 nm², which is **0.70× the dome — less than the dome, not
2.4× more**.

I have left the Round 3 entry in ROADMAP.md struck through rather than edited
away, with a pointer to this round. The caveat was recorded honestly at the
time; the lesson is that recording a caveat is not a substitute for doing the
calculation, and deleting the evidence of that would remove the lesson.

### What this does *not* say

It does not refute Haselwandter & MacKinnon. Their argument concerns the
footprint's contribution to tension *sensitivity* — the area released between
closed and open states — and absolute stored area was never a test of it. Round
3's error was rhetorical as much as numerical: it presented our absolute-area
ratio as "the quantitative form of" their claim, which it is not. The defensible
statement is narrower: at PIEZO1's contact slope the linearised footprint is not
quantitatively usable, and corrected, dome and footprint hold comparable excess
areas with the dome slightly larger.

Suite 254 → 272 passing; GUI smoke test clean.

## Round 19 — hydrophobic gating, and what a heuristic is actually for (2026-08-06)

Round 1 built the pore profiler and it worked: 8YEZ's bottleneck is 0.095 nm,
11ZC's is 0.330 nm, and the profiler rediscovered the curated gate and CTD
constrictions from coordinates alone. But radius is a weak predictor of
conduction — Rao et al. put its AUROC at **0.59**, barely better than a coin —
because a pore can be wide enough for a hydrated ion and still block when a
hydrophobic neck expels liquid water. Their heuristic, combining radius with
local hydrophobicity, reaches **0.91**. This round implements it.

### Getting the boundary rather than redrawing it

The paper gives the construction, the 1 RT = 2.6 kJ/mol contour and the
Σd > 0.55 cutoff, but not the classification line itself — that lives in a
figure. Digitising a figure by eye would have been precisely the kind of silent
correctness bug Round 17 was about.

It turned out not to be necessary. CHAP is **MIT licensed** and its repository
ships `heuristic_grid.json`: the actual 100×100 water free-energy landscape over
(hydrophobicity, radius), built from ~600 MD simulations. It also ships the
exact normalised Wimley–White scale the landscape is indexed by, and the default
kernel bandwidth (0.35 nm). So the published artefact is used directly, as a
download rather than a commit, and analyses degrade to "unavailable" without it.

An independent check that we read it correctly: our extracted 1 RT contour gives
a critical radius rising from **0.10 nm** at the hydrophilic end to **0.43 nm**
at the hydrophobic end, against the paper's prose "hydrophilic pores wet below
0.2 nm; hydrophobic ones can hold a barrier out to ~0.4 nm". Nothing in our code
was fitted to those numbers.

### The bug that returned confident nonsense

My first hydrophobicity profile averaged residues in a sphere around each probe
centre. It ran, produced smooth plausible values, and gave 8YEZ a score of 0.45
— just under the cutoff, so the closed structure was called **conductive**.

The tell was the range. CHAP smooths **along the pore coordinate** over
pore-facing residues, and the published grid spans −0.45 to +0.30. My 3-D
neighbourhood was 1.85 nm wide, which pulls in the entire shell of residues
surrounding the lumen, and the profile collapsed into a band from −0.12 to
+0.02. Every energy was then read out of the landscape at a coordinate the
landscape was never built on. Nothing errored; the numbers merely meant nothing.

Rewritten as a proper Nadaraya–Watson average along the axis, using side-chain
centroids because it is the side chain that faces the lumen, the range opens to
−0.635…+0.229 and 8YEZ scores 0.82 — non-conductive. There is now a test
asserting the profile uses most of the grid's range, because "plausible but
compressed" is what this failure looks like from the outside.

### Right answer, right reason

The roadmap asked for the right answer *for the right reason*, which is a
demand for a control rather than an assertion. So: hold every radius fixed and
replace the hydrophobicity scale with a uniform hydrophilic value. If the closed
verdict were a radius threshold wearing a disguise, the score would not move.
It goes from **0.82 to 0.00** — conductive. The verdict is chemistry.

The sharpest single fact: 8YEZ's F2451 and V2454 sit at **0.325 nm** and are
called dewetted; 11ZC's *bottleneck* sits at **0.330 nm** and is called wet.
Same radius, opposite verdict. And the flagged set — F2451, V2454, R2467,
F2468 — is the curated hydrophobic gate and cytoplasmic constrictions, which the
heuristic never sees.

### The limitation I would have missed by testing only what was asked

The round specified two structures. Running five exposed something the two
would not have: **7WLU and 8IXO have 0.098 nm bottlenecks — less than a water
molecule's 0.15 nm radius — yet score 0.11 and 0.30, i.e. open.**

That is not a bug in the heuristic; it is what the heuristic is *for*. It
answers "would water dewet here?", not "does water fit here?". Rao et al. built
it to find hydrophobic gates, which are by definition blockages *without* steric
occlusion. A pore too narrow for water never poses the wetting question.

Merging the two into one verdict would have hidden this, and would have made the
heuristic look like a general conduction predictor it does not claim to be. So
`WettingPrediction` exposes `hydrophobic_gate` and `sterically_occluded`
separately and `conductive` requires neither. With both, all five states come
out right: 8YEZ and 7WLT shut on both counts, 7WLU and 8IXO shut on sterics
alone, 11ZC open.

Suite 272 → 288 passing; GUI smoke test clean.

## Round 20 — what the null was entitled to claim (2026-08-06)

Most of this round's list was already done. Rounds 6–7 delivered the permutation
test, the bootstrap effect size and the pre-registered decision rule; Round 7
recorded the null. What was missing were the two questions asked *around* a
result rather than by it — could the design have detected the effect, and did we
look in too many places — plus cross-validation. Those went into a new
`analysis/design.py`, kept separate from `validation.py` because "did it work?"
and "could it have worked?" are different questions and conflating them is how
underpowered nulls get overclaimed.

### The finding

Simulating the pre-registered test at Round 7's actual group sizes, 16
gain-of-function against 9 loss-of-function, one-sided at α = 0.05:

| effect | Cliff's δ | power |
|---|---|---|
| **observed in Round 7** | −0.083 | **0.13** |
| small | −0.11 | 0.16 |
| medium | −0.28 | 0.35 |
| large | −0.43 | 0.60 |

**80% power is reached only at |δ| ≥ 0.55**, past 'large' on the usual
thresholds. I ran it twice — once from a normal model, once resampling the
observed heavy-tailed ΔΔG values, because a difference-in-means test loses power
on heavy tails and the normal model would flatter the design. They agree: 0.55
and 0.56.

So the Round 7 null **excludes a large mechanical effect and is close to
uninformative about a small or medium one.** That is a real qualification and it
was not stated at the time. I have added it as `VALIDATION.md` §6b, explicitly
marked as not amending §§1–3: the result stands, its scope is now bounded. The
existing §6 diagnostic — 99.8% of ΔΔG variance is between-position — remains the
mechanistic explanation and is independent of this. Both are true. The power
limit is simply the one nobody wrote down.

### The constraint that binds the rest of the project

At 80% power with equal groups: **42 variants for a large effect, 98 for a
medium one, 600+ for a small one.** Twenty-five survive Round 7's inclusion
criteria, and relaxing every criterion cannot reach forty-five.

So a *confirmatory* test of anything below a large effect is not available from
this variant set, however good the predictor gets. Round 22 is now explicitly
bound to declare in advance whether it is confirmatory-for-a-large-effect or
exploratory, and I have added Round 27 — expand the phenotyped set — because the
binding constraint on this project's central claim turns out to be **data, not
method**. That is worth knowing at round 20 rather than round 40.

### Why the protocol had to be written now

`docs/NEGATIVE_RESULT_PROTOCOL.md` is deliberately written *between* Round 7's
null and Round 22's hypothesis. Written after Round 22 it would be a rule fitted
to a result already seen, which is the thing it exists to prevent.

The failure mode it guards against is not dishonesty, it is drift: a null comes
back, and there is always one more defensible adjustment — drop the noisy
variants, use the normalised score, go two-sided. Each is reasonable alone.
Together they are an unrecorded search, and the p-value at the end means nothing
like what it appears to. With 68 curated variants and one primary claim, that
resource is finite and non-renewable.

The multiplicity worked example makes the point better than the argument does.
Six candidate predictors with plausible p-values: three clear 0.05, and after
Benjamini–Hochberg **none survives**. Reporting AlphaMissense at p = 0.012
without its family would be a false discovery manufactured by looking six times.

### A sign error caught by its own diagnostic

`power_curve` injected the effect into the wrong group. `shift_for_delta` is
defined with a ~ N(0,1) and b ~ N(shift,1), so the displacement belongs on group
b; I added it to a. The resulting power curve rose monotonically with effect
size, saturated at 1.0, and was entirely plausible — it was simply the power to
detect the opposite direction.

It was caught because `power_curve` measures the Cliff's delta it actually
achieved rather than assuming it got what it asked for. That diagnostic existed
only because Round 18 had just taught the same lesson: the useful check is the
one that does not reuse the thing being checked. It is now a test.

### Cross-validation

Leave-one-out on the Round 7 predictors gives AUROC 0.535 out-of-sample against
0.542 in-sample, optimism +0.007. Small because the default combination is
unsupervised, and reassuring: there was no hidden overfitting inflating the
original number. The machinery matters more for Round 22, where a *fitted*
combination of six predictors on 25 variants would otherwise measure how well 25
points can be fitted.

Suite 288 → 311 passing; GUI smoke test clean.

## Round 21 — the engine reaches the interface (2026-08-06)

Twenty rounds built an engine the GUI could not reach. It could show a
structure, measure a dome and animate normal modes; the pore profiler, the
pockets, conservation, allostery and every report were CLI-only. This round
closes that gap.

### The thing worth guarding against

A GUI that recomputes what the engine already computes becomes a second
implementation, and second implementations diverge. So the workers call the
same functions the CLI calls, and a test asserts the pore worker reproduces
`pore_profile` and `predict_wetting` *exactly* rather than producing something
plausible. The smoke test then reads the panel's own label back: bottleneck
0.95 Å, Rao score 0.82, "non-conductive (sterically occluded + hydrophobic
gate)" — identical to `piezo1 hydration 8YEZ`.

### Two axes, because that is the whole point

Round 19's finding is that radius alone predicts conduction at AUROC 0.59 and
radius-with-hydrophobicity at 0.91. A plot that puts both on one y-scale makes
the hydrophobicity trace a flat line at the bottom of a 0–10 Å axis, which
destroys exactly the comparison the plot exists to show. So the widget has
independent left and right axes, with the left anchored at zero (a radius of
zero is meaningful) and the right not (hydrophobicity is signed, and anchoring
it at zero squashes the range carrying the signal).

I wrote it as a QPainter widget rather than adding matplotlib or pyqtgraph.
One plot type, has to repaint inside a dock at interactive rates, has to match
the dark theme, and both libraries drag in their own event-loop integration —
more to configure than to write.

### Where the colour schemes went, and why not where the roadmap said

The roadmap asked for conservation and PRS as "first-class colour schemes",
which I read as the main colour dropdown next to Domain and Chain. I put them
in the Analysis dock instead, appearing once computed.

The reason is that they are not the same kind of thing. Domain and chain are
properties of a loaded structure and are available the instant it loads.
Conservation needs 61 orthologs fetched and aligned; PRS needs normal modes
computed first. An entry in the main dropdown that silently does nothing until
some other panel has been used is a dead control, and a user who selects it and
sees no change has been told the feature is broken.

They do go through the existing `ColorBy.VALUE` path rather than a new one —
the same mechanism mode-displacement colouring already used.

### Two traps

**Unmeasured residues take the floor, not zero.** Conservation runs about
0.6–1.0. Giving an unmeasured residue a zero stretches the colour scale across
a range containing no data and washes out every real difference — while looking
like a perfectly valid map. The value map now fills with the minimum of what it
was actually handed.

**Restoring a session rebuilt the scene four times.** Style, colouring, ligands
and atom size each emit a signal. Setting all four from a saved session fires
four rebuilds of a 120k-atom trimer, which is visible as flicker through three
intermediate views nobody asked for. `set_state` blocks signals, sets the
widgets, then emits once.

Also worth recording: conservation drops positions with ortholog coverage below
0.7. Those values measure how well the alignment covered the position, not
selection pressure on it, and left in they paint a confident band across the
unresolved distal blade. 2484 residues survive, mean 0.770, 901 above 0.95.

### Testing a GUI

Qt refactors have silently broken this application twice, so the tests run real
widgets on the offscreen platform rather than mocks. Mocks would not have caught
either breakage — both were real widgets behaving differently from how the code
assumed, which is precisely what a mock encodes rather than tests.
`scripts/screenshot_app.py --analysis` now drives the new panel end to end and
round-trips a session through the controller.

Suite 311 → 330 passing.

## Interlude — the interface, driven by using it (2026-08-06)

Round 21 put the engine behind panels; actually running the application turned
up four things that no amount of testing from inside had shown. All four came
from the user opening it on their own machine, which is the only way some of
these surface.

### The window did not fit the screen

`resize(1680, 1000)` was hard-coded. On a laptop or a 1080p display that puts
the title bar above the top of the screen, and on some window managers the
window then cannot be dragged back or resized. It now opens at the smaller of
the preferred size and 95% of the *available* geometry — which already excludes
the menu bar and dock — and centres itself.

That was necessary but not sufficient. Qt will not shrink a dock below its
content's minimum, so the tallest panel still set a floor on the whole window.
Every panel is now wrapped in a `QScrollArea`, and the window's minimum size
hint dropped to 535×258. There is also `--geometry WxH` for when screen
detection misjudges a multi-monitor setup.

### The model jumped whenever anything was selected

Selecting a domain, site or variant called `_focus_residues`, which moved the
camera pivot to the selection's centroid. On a 2500-residue trimer that means
the whole structure lurches every time you click a list entry.

The honest fix was not to remove it — on a structure this size a highlight can
easily be behind the protein, and centring is genuinely useful — but to make it
a choice. Options → *When something is selected* now offers **keep the view
still** (the new default), **centre on the selection**, or **centre and zoom**.
The zoom variant deliberately preserves the current *orientation*: reframing
rotation as well would throw away the view the user had set up, which is the
same complaint arriving by another route.

### Panels could not be rearranged

Docks were created with Qt's defaults, which confine each one to the two areas
it was born in and give no way to float or close it. All five are now movable to
any edge, floatable into their own windows, closable, and tabbable.

Making them closable creates the obvious problem: a user who closes two panels
and floats two more has no route back. So the shipped arrangement is captured
once at startup — *before* any remembered layout is applied, so Reset is always
the application's own layout rather than whatever was last left — and
View → Reset layout (Ctrl+R) restores it. The trap here is `objectName`:
`saveState`/`restoreState` key on it, and without one Qt silently declines to
restore that dock, which looks like a corrupt settings file rather than a
missing property.

Layouts persist between runs, with the geometry clamped to the current screen on
restore — otherwise a layout saved on a large external monitor reopens mostly
off-screen, which is the first bug again by a different path.

### Menus, help and tooltips

Options and Help menus now exist. The help window is non-modal, because the
point of a feature guide is to be read *while* driving the application.

The guide has seven topics, and the one that matters most is **"What this
application will not do"** — the pre-registered blind test that returned
p = 0.234 and AUROC 0.542, the diagnostic that 99.8% of the predictor's variance
is between-position, the power limit that lets the null exclude only a large
effect, and the footprint number that turned out 3.5× too large. A test asserts
those figures are present, because a help file is exactly where an inconvenient
result would quietly stop being mentioned.

Tooltips now carry provenance rather than restating labels: the dome button says
what number it should produce and against which publication, the hydrophobicity
checkbox gives 0.59 against 0.91, the mode selector explains why only A modes
can couple to isotropic tension.

Suite 330 → 342 passing, and the GUI smoke test now checks the two behaviours
that were reported — that focus-off leaves the camera still and focus-centre
moves it, and that reset restores hidden and floated panels.

## Round 22 — the second null, and why it was foreseen (2026-08-06)

Round 7 tested a mechanical predictor and returned a null whose diagnostic was
precise: 99.8% of its variance was between-position, so it reported *where a
residue sits* rather than *which substitution occurred*. Round 17 brought in
four substitution-aware predictors. This round asked whether they do better.

The pre-registration was written and **committed in its own commit before the
test ran**. That is not ceremony. A rule written after seeing a p-value is a
rule fitted to it, and with 68 curated variants and one primary claim, the
resource being spent is finite.

### The design was worse than Round 7, and I knew that before testing

Of 39 variants carrying a directional label, only 26 are single-residue
substitutions. **Eleven of the thirteen dropped are loss-of-function** — stop
codons, frameshifts, deletions.

That is not a curation artefact. It is how loss of function actually happens:
you break the protein. And no missense predictor can score a stop codon.

So the usable design is 20 GoF against 6 LoF, and 80% power needs |δ| ≥ 0.61
where Round 7 needed 0.55. Power at a conventionally large effect is 0.52 — a
coin toss. Under the Round 20 protocol that forces a declaration, and I declared
the round **exploratory**: at this n a "confirmatory null" would exclude only
effects beyond large, and calling that confirmation of anything would overstate
it.

### The hypothesis, and the objection recorded alongside it

Primary endpoint: FoldX ΔΔG. The mechanism is that loss of function can be
achieved by breaking the protein, while gain of function cannot — a channel that
opens too readily must still fold, traffic and gate, so gain-of-function
variants are constrained to be structurally tolerable. Prediction: LoF more
destabilising.

I also recorded, in the same document, the objection: excluding the truncating
LoF variants removes precisely the "break the protein" mechanism, leaving a LoF
subset selected for **not** being truncating. Writing that down first is what
stops it later being a post-hoc rescue.

### The result

**Cliff's δ = −0.211, CI [−0.684, +0.298], AUROC 0.395.** Mean ΔΔG 0.767 for
LoF against 1.309 for GoF — the point estimate runs *opposite* to the
hypothesis. Mean difference, Cliff's delta and AUROC all agree on the direction,
so it is not a sign error. The interval spans zero, so it is uninformative
rather than a reversal.

The secondary family: nothing significant, nothing close, smallest adjusted
q = 0.448. The pre-registered expectation that AlphaMissense, EVE and ESM-1b
would show *no* separation held — a single pathogenicity axis has no room for
direction when both classes are pathogenic. Combining all five by equal weights
gives leave-one-out AUROC 0.535 with zero optimism: nothing recovered.

### What the foreseen objection now looks like

It looks like the explanation. A destabilisation predictor was asked about the
one subset of loss-of-function variants that does not act by destabilisation,
because the subset that does was excluded for being unscoreable. That it was
written down before the run is the only reason it can be offered at all.

### Where this leaves the project

Two nulls, two different predictors, both pre-registered, neither revised. The
central claim is not refuted — it is **untested at adequate power**, and after
two attempts the constraint is clearly data rather than method. Fifty-three
variants would be needed for a large effect at this class ratio, 148 for a
medium one, against 26 available. Round 27 is now the round that matters.

There is also a sharper question hiding in the counts. Eleven of seventeen
loss-of-function variants are truncating, and predicting *those* from structure
is close to trivial. The hard and interesting question — distinguishing "opens
too easily" from "does not open" among missense variants — has a sample size of
six. That is the honest scope of what this project can currently ask.

Suite 342 → 352 passing.

## Interlude II — measurement furniture, sequences and overlays (2026-08-06)

Four more requests from actually using the application. Each turned out to have
a correctness question hiding inside what sounded like a display feature.

### Scale bars and a clock

A structure screenshot without a scale bar states nothing quantitative; a reader
cannot recover the scale of a perspective projection by eye. The bar rounds to a
tidy length — 1, 2, 5, 10, 20 … Å — because a bar labelled 47 Å is unreadable.
It is exact in the plane through the camera pivot, which is where the molecule
is, and the docstring says so rather than implying a perspective bar is exact
everywhere.

The animation clock raised the more interesting question. For a **mode sweep**
it reports elapsed seconds and the fraction of a display cycle, explicitly
labelled *not a physical period*: the frequency of an elastic-network mode is
not time. For a **morph** it refuses to report seconds at all and gives percent
along the path, because a morph is an interpolation between two endpoints and a
seconds axis would imply kinetics the model does not contain.

### Presentation mode

Full screen hides the panels and menu. The detail worth recording is that
leaving it restores each panel's **previous** visibility rather than showing
them all — someone who had closed a panel before presenting should not find it
reopened afterwards. And Escape leaves, because with the menu bar hidden there
is otherwise no visible way out.

### The sequence window

Three kinds of sequence are involved and they are **not** interchangeable:
UniProt (the reference numbering variants are quoted in), structure (only
resolved residues, with gaps, starting at 570 for 8YEZ), and translated CDS.
`NamedSequence` carries its own positions for exactly this reason — a viewer
that showed "the sequence" without saying which one is a numbering bug waiting
to happen, and this project has already been bitten by three of those.

For DNA I fetched the real Ensembl transcript rather than back-translating. A
back-translation would look like a gene, invent every codon choice, and make
silent variants unrepresentable. The check that the right transcript was
fetched: **the human CDS translates to Q92508 at 100% identity over 2521
residues.** Mouse gives 99.9% — three differences at 147, 229 and 1572, being
the reference genome against the TrEMBL submission. That is recorded in a test
rather than smoothed over.

I also guessed the mouse transcript ID from memory and it 400'd, which is the
same failure mode as the six wrong PMIDs in Round 8. Looked it up properly.

Comparison offers global alignment or pairing by residue number, and the second
is not a lesser option: when two sequences already share a numbering, an aligner
can slide a run of residues to buy score and manufacture differences that are
alignment artefacts. Pairing by number is refused across species, where it would
be meaningless.

### The structure overlay

The comparison this project needs constantly — curved against flattened, wild
type against variant — is two structures in one frame.

The correctness question is protomer correspondence. Round 4 found four
deposited entries labelled in the reverse rotational order. Overlaying 7WLU on
7WLT, correspondence search rematches to (2, 1, 0) and gives 12.3 Å; taken at
chain-label face value the two sit **90.7 Å** apart. A viewer that trusted
labels would display an enormous conformational change that does not exist —
and it is exactly the conclusion someone would want to draw. Both numbers are
reported.

A second guard: cross-species overlays are refused outright. Residue numbers are
the join key, and human and mouse numbering do not correspond, so the fit would
be confidently wrong rather than obviously wrong.

Building the protomer blocks per structure also failed, because two entries
rarely resolve the same residues and `match_protomers` cannot compare blocks of
different length. The basis has to be shared across both.

Suite 352 → 376 passing.

## Round 23 — one command, and a guard against the docs going stale (2026-08-06)

Aim A5 says a fresh clone plus an environment plus a fetch should reproduce the
whole working state. The packaging half of that is routine: `pyproject.toml`,
two pinned lock files, a `Makefile` whose targets each wrap `conda run` so they
work from a bare shell. One decision worth recording — **the GUI dependencies
are an optional extra, not a requirement**. Everything below `render` runs
headless, and that is precisely what makes the science testable without a
display; making PyQt mandatory would quietly give that up.

### The half that actually matters

This project states a great many specific numbers in prose. A dome radius of
9.7 nm. A mode overlap of 0.705. A half-activation tension of 2.71 mN/m. A
footprint of 179 nm² that replaced one of 622. Two null results with their
p-values and effect sizes.

**Prose does not fail a test suite.** A solver rewrite, a re-fetched structure
or a changed default can leave `docs/SCIENCE.md` confidently asserting a value
the code no longer produces, and nothing in the repository would notice. Given
how many recorded numbers have already needed correcting here — the footprint
was wrong by 3.5×, the Bessel ratio by 2.5×, the decay length by a factor of
three — the risk is not hypothetical.

So `piezo1/analysis/claims.py` holds seventeen claims. Each names the documented
value, its tolerance, the document it appears in, and a callable that recomputes
it from scratch. `make verify` runs them in about ten seconds. All seventeen
currently reproduce.

A claim is not a test, and the distinction is the point. Tests assert the code
behaves; claims assert **the documentation is still true of the code**. They
fail for different reasons and both are worth having.

### Two design decisions inside that

**Frozen claims.** The four recorded validation numbers — Round 7's p-value and
AUROC, Round 22's effect size and AUROC — are marked `frozen`. If one drifts,
the report prints an explicit instruction *not* to edit the document to match,
but to work out why the computation changed. The obvious way to resolve a
failing claim is to update the prose, and for a recorded null result that would
be exactly the wrong move.

**Skipped is not drifted.** A claim that cannot run because a structure has not
been downloaded reports as skipped. Conflating that with drift would make every
fresh clone look like it had a broken scientific record, and the person seeing
it would learn to ignore the report.

### Testing the detector

A drift detector that has never detected drift is decoration. There is a test
that feeds it a deliberately wrong claim and requires it to complain, and
another that feeds it a claim raising `FileNotFoundError` and requires *skipped*
rather than *drifted*. Both would have passed vacuously if the registry simply
returned success.

Suite 376 → 390 passing.

## Parameter registry — every number gets a citation (2026-08-06)

An audit of `physics/`, `structure/` and `analysis/` found **203 numeric
literals**. Some were physical constants from named papers; some were
convergence tolerances; and from the outside there was no way to tell which was
which. That is the problem: a constant written into a function default is
invisible. You cannot list them, show them to a user, or trace one to a paper.

This project has already had to correct several numbers that were invisible in
exactly that way — the footprint area wrong by 3.5×, the Bessel ratio by 2.5×,
the biharmonic decay length by a factor of three.

### What was built

`piezo1/resources/parameters.json` holds **61 parameters** across 12
categories, each with a unit, bounds, a kind (physical / empirical / method /
convention), a description and a citation. It is authored by
`scripts/build_parameters.py`, in the same shape as `build_variants.py`:
authored content, validated on the way out, committed as a resource. Keeping it
as data means the whole parameter set can be read and diffed without opening a
module.

**The gate is the point.** A citation must resolve to a key in
`references.json`, or be one of five sentinels — `method_choice`,
`measured_here`, `derived`, `convention`, `unverified` — each of which *obliges
the entry to say why* in `source_note`. The build refuses to write otherwise. It
caught eight entries where I had left the note blank, and I filled them in
rather than weakening the gate. Two references had to be added first (fpocket
and Shrake–Rupley); the reference builder's own `expect` gate then rejected my
first attempt at Shrake–Rupley because I picked a word that is not in the title.

31 of the 61 cite a published paper. The other 30 are method choices, and saying
so explicitly is more honest than dressing them up.

### Overrides are tracked, not silent

Modules resolve through the registry **at call time**, so a change takes effect
on the next call. That makes the values genuinely editable — and immediately
raises the real problem: a number computed with a changed parameter is not the
number in `docs/SCIENCE.md`, and nothing would have said so.

So three things enforce it:

- `verify_claims` **refuses to run** against a modified registry. Every
  documented number was produced at the defaults; recomputing with a changed
  value would report drift the *user* caused, and the obvious reading of that
  report is that the code is broken.
- Reports carry a warning banner **at the top**, not in the provenance footer,
  and `Provenance` gained a `parameter_overrides` field.
- The GUI keeps a persistent amber indicator in the status bar.

### The audit is what makes it a rule

A convention nobody checks decays at the first hurry. `piezo1/parameter_audit.py`
scans the three scientific packages and fails on any numeric literal that is
neither registered nor listed in `EXEMPT` **with a stated reason**. It started
at 18 unaccounted; 11 were registered-but-unmigrated and 7 were genuinely
unregistered — those became new parameters rather than new exemptions.

The exemption categories matter as much as the registrations. A convergence
tolerance, an iteration cap, a random seed and a zero-initialised dataclass
field are implementation details, and pretending otherwise would bury the real
parameters in noise. What the audit insists on is that the exemption is a
decision someone made and wrote down.

There is a test that invents a module containing `binding_energy = 7.25` and
requires the audit to catch it — a detector that has never detected is
decoration.

`make audit` and `make params` are build steps, and both run inside
`make reproduce`. The rule is written into `CLAUDE.md`. Suite 390 → 408 passing,
and all 17 documented numbers still reproduce, which is the check that the
migration changed nothing scientific.

## Round 24 — performance, and profiling before optimising (2026-08-06)

The roadmap named three slow paths. Profiling found that two of them were not
the problem and the real one was not on the list.

- **PRS**: 0.52 s. Never a problem.
- **Pocket detection**: 4.2 s, not the ~10 s the roadmap assumed.
- **Ensemble PCA**: 12.4 s — but **99% of that was mmCIF parsing**, not the
  PCA at all.
- **SASA**: 7.5 s, and the roadmap did not mention it.

Recording that is worth more than the speedups. A round spent optimising PRS
because a note said it was slow would have been a round wasted.

### The rule for this round

An optimisation that alters a number is a bug, not a speedup. Every change here
is a reformulation with the same value, and the tests assert **identity**, not
closeness.

### What was actually slow

**SASA, 7.54 → 1.27 s.** The inner loop built a (256 × neighbours × 3) array per
atom and took a square root. Two observations remove both: `d ≥ r` and `d² ≥ r²`
decide the same way for non-negative values, so the root was never needed; and
expanding `|t − x_j|² = |v_j|² + r_i² + 2r_i(p_k · v_j)` turns the 3-D broadcast
into a single BLAS product. Bit-identical across all 31,599 atoms —
`np.array_equal`, not `allclose`.

**The mmCIF tokenizer, and through it the ensemble, 12.35 → 2.05 s.** The
tokenizer walked characters in a Python loop, 321,913 times per ensemble load.
But 99.5% of mmCIF lines contain no quote and no comment, and for those
`str.split()` is exactly equivalent — its default whitespace set is the same one
the function uses, and it discards empty fields the same way.

I checked that claim rather than asserting it: 245,528 lines of deposited
structure, **zero mismatches**. This is the function whose whitespace handling
once shifted every column by one and produced `invalid literal for int(): 'ATOM'`,
so the careful path is left exactly as it was and merely bypassed.

**Conservation, 3.67 s → 0.003 s.** Aligning 61 orthologs dominated anything
needing conservation, and the result depends only on the sequences. The cache is
keyed on a **content hash** of the reference and the ortholog sequences rather
than on a filename or timestamp, so re-fetching or changing the reference
invalidates it automatically. A cache that can go stale is worse than no cache:
it would report last week's conservation against this week's alignment, and
nothing would look wrong.

### The check that matters

Overall 33.8 → 12.5 s, and every number unchanged: SASA total 197490.5582 Å²,
top pocket 6593.6 Å³, ensemble PC1 0.9000, pore bottleneck 0.9518 Å. All 17
documented numbers still reproduce, which is the guard Round 23 built for
exactly this kind of change. The suite itself went 118 → 97 s.

Timing assertions in the tests are loose ceilings only. Pinning a runtime would
fail on a slower machine for no scientific reason; what is worth pinning is that
the fast path and the careful path agree.

## Round 25 — the teaching layer, with nothing narrated (2026-08-06)

Aim A1 says this should be a learning instrument, and of the six aims it had
received the least attention. The tour walks the mechanism in eleven steps:
trimer, blades, dome, footprint, lever, gate, open state, normal modes, gating
energetics, a variant, and what the project cannot do.

### The rule that shapes it

**Every number a step states is computed when the step runs.**

A tour is prose, and prose is where numbers go to rot. Writing "the dome radius
is 9.7 nm" into a tour step would create a fourth place for that value to live —
beside the code, `docs/SCIENCE.md` and the claims registry — and the fourth copy
is the one nobody remembers to update. So each step carries a callable that
reads whatever the application has actually computed, and where a step quotes a
published comparison it reads it from the **parameter registry**. There is a
test that changing `dome.published_radius_closed` changes the tour text, which
is the proof it is not a literal.

The controller calls the same controllers the panels use. A teaching tool that
quietly disagreed with the application it was teaching would be worse than none.

### Where it ends

On the failures. The last step states both null results — p = 0.234, AUROC
0.542, and the second test's −0.211 with an interval spanning zero — and the
power limit that means the first excludes a large effect and little else. There
is a test asserting those figures are present, because a tour is exactly the
sort of document where an inconvenient result quietly stops being mentioned.

A learning instrument that only shows its successes teaches the wrong lesson,
and the lesson worth teaching here is that the interesting question about R2456
is why four substitutions at one residue do not all do the same thing.

### Robustness

Two tests exist because a tour must never take the application down: every step
must degrade to a readable message before anything has been computed, and none
may raise when handed junk. The GUI smoke test walks all eleven steps.

Suite 419 → 431 passing.

## Structure composition — what is actually in these files (2026-08-06)

Asked to check that loading is consistent across structures and that bound
elements are handled correctly. All 23 entries load, but auditing what is *in*
them turned up something worth fixing.

### Seven entries contain a second protein

Six carry three copies of **MDFIC** — a 21-residue cysteine-rich peptide
(`CCESSDCLEICMECCGICFPS`) that is a genuine auxiliary subunit inserting into the
pore (Zhou et al. 2023, already in our bibliography). 6B3R carries three
16-residue poly-UNK chains.

These are *protein*, so a protein mask includes them. Worse, MDFIC is numbered
**226–247**, which sits inside PIEZO1's own numbering. A selection keyed on
residue number alone would pool the two.

Nothing has actually been wrong: PIEZO1 resolves from residue 570 upward in
exactly the entries that carry MDFIC, so the ranges happen not to overlap. But
that is luck, not design, and the `> 300 C-alpha` rule that kept them out of
protomer blocks was a coincidence of scale rather than a statement about what
they are. There is now a test asserting no auxiliary residue number reaches a
protomer basis, so if a future entry resolves further into the N-terminus it
fails loudly instead of quietly averaging MDFIC into a mode.

`core/entities.py` classifies every atom. The principal-chain rule is
**relative to the largest chain**, not an absolute threshold — 4RAX is a lone
227-residue domain and is the entire structure, while a 21-residue peptide
beside a 1,280-residue protomer is not a protomer. Across all 23 entries the
classifier returns 1 or 3 protomers, never anything else.

The heterogens were looked up in the PDB chemical component dictionary rather
than inferred: PLX, PEE, P5S and L9Q are phospholipids, D12 is dodecane, NAG is
a glycan on PIEZO2.

### Display is now a choice, and only display

The Model panel lists the categories actually present with their atom counts,
each independently switchable. Only present categories get a control — a
permanent list mostly greyed out would say nothing about what you are looking
at, and what is in *this* file is the thing worth surfacing.

There is a test that hiding a category cannot change an analysis. Display and
computation are separate questions, and the analyses always use the channel
protomers whatever is drawn.

### Three structures that legitimately cannot answer

Extending the GUI smoke test beyond 8YEZ exposed that it assumed human
numbering and full side chains. Fixed, and the three failures turned out to be
correct behaviour needing to be *reported* rather than *fixed*:

- **11ZC** is a backbone-only model — C, CA, N, O and nothing else — so there
  are no sulfur atoms to measure.
- **6KG7** is PIEZO2. PIEZO1 numbering does not transfer to a paralogue; the
  second position is isoleucine.
- **4RAX** is one isolated domain, so there is no trimer to fit a dome to.

Each is now skipped with its reason. Conflating "cannot answer" with "answered
wrongly" is the same mistake the claims registry avoids, and it makes correct
behaviour look broken.

### An observation about the disulfide

Converting the curated C2411–C2415 disulfide into mouse numbering (2437/2441,
via the sanctioned path rather than a constant) shows both residues are cysteine
in every mouse entry — a 26-residue offset landing on cysteine twice is not
chance, so the numbering map is doing its job.

But the **bond is only modelled as formed in some of them**. Human 8YEZ and
8YFC give 2.04 Å, and so does mouse 4RAX — the 1.45 Å X-ray structure of the
isolated CED domain. The full-length mouse cryo-EM models put the same sulfurs
**5.2 Å apart** (7WLT, 7WLU) or 6.7 Å (6B3R). That is a resolution and
modelling difference between depositions rather than a biological one, and the
smoke test now records it instead of failing on it.

Suite 431 → 445 passing.

## Round 26 — the substitution finally enters the mechanics (2026-08-06)

Round 7's blind test failed and its diagnostic was precise: 99.8% of the
mechanical ΔΔG's variance was between position, not between substitutions. Two
rounds since have worked around that by adding other predictors. This one
attacks the cause.

### The cause was algebra, not statistics

The old model scaled **every contact of the mutated residue by one number**. So

    ΔΔG = (s − 1) · Q(position)

which is a rank-one product: the substitution enters only as a multiplicative
scalar, and the positional factor Q — contact count times local strain — varies
enormously between positions and not at all within one. Four substitutions at
R2456 could therefore differ only by a factor and had to rank every position
identically.

Seeing that written out matters, because it says no amount of refining `s` could
ever have helped. A better volume term, a substitution matrix, a fitted
sensitivity — all of them are still a scalar. The separability itself had to go.

### The repair

Scale each contact individually, by properties of the new residue *and of the
partner it touches*:

- **packing** — as before, but weighted by how close the contact is;
- **charge** — felt only at contacts with charged partners, so losing R2456's
  salt bridge softens the contact to an aspartate and not the one to a leucine;
- **hydrogen bonding** — donor against acceptor, both ways;
- **proline** — stiffens sequence-local contacts specifically, because that is
  where a backbone restraint is felt;
- **glycine** — softens whatever the side chain was mediating, since there is no
  longer a side chain.

Different substitutions now perturb *different subsets* of contacts. Measured
directly rather than assumed: the per-contact patterns for R2456H/K/P/C
correlate 0.62–0.98, where under the old model they were 1.00 by construction.

An elastic network is a mechanical model, and its springs are an effective
stiffness standing in for packing, hydrogen bonds and ion pairs together.
Letting that stiffness depend on charge is a statement about what the springs
represent, not a claim to have added electrostatics, and the module says so
rather than dressing it up.

### The criterion, and reporting both numbers

The success criterion was fixed in the roadmap before the work: within-position
variance above 20%, measured on the multiply-substituted positions. **It is met:
4.9% → 52.5%.**

Across *all* 35 substituted positions the figure is 2.4%, up from 0.8%. Both are
reported because the difference is instructive rather than embarrassing: 29 of
those positions carry a single substitution and contribute exactly zero
within-variance by construction, so including them pushes the statistic down for
arithmetic reasons that say nothing about the model. The criterion named the
multiply-substituted subset for exactly this reason.

### What this does not license

**No phenotype comparison was run.** The criterion was about variance
decomposition; whether the new distinctions point in the *right* direction is a
different question, it has not been asked, and asking it needs a new
pre-registration under `NEGATIVE_RESULT_PROTOCOL.md` §7. The temptation to
"just check" now that the model can tell substitutions apart is precisely what
that document exists to resist.

Round 7's recorded result is untouched: its script passes no sequence, so the
model falls back to the uniform scale, and there is a test asserting the two
agree to 1e-12.

Suite 445 → 462 passing; 18 documented numbers reproduce.

## Round 27 — more variants, and what kind of evidence they are (2026-08-06)

Round 20 measured the constraint: 42 variants for a large effect, 98 for a
medium one, against 25 usable. Round 22 then ran into it from the other side,
with six loss-of-function missense variants. This round went looking for more.

### ClinVar does not answer the question directly

It reports **pathogenicity**, not **direction**. A gain-of-function and a
loss-of-function variant are both "Pathogenic".

What makes direction recoverable for PIEZO1 specifically is that its two
diseases have opposite, well-established mechanisms — dehydrated hereditary
stomatocytosis is dominant and acts by slowed inactivation; generalised
lymphatic dysplasia is recessive and acts by loss of function. So the condition
implies the direction.

That is **weaker evidence than measuring the current**, because it assumes the
variant acts by the usual mechanism for its disease. So it is recorded as a
separate evidence level rather than pooled, and `build_analysis_set` defaults to
the measured set: a caller who does not think about evidence strength gets the
smaller, better-supported answer rather than the larger one.

### The ambiguity is bigger than expected

**Eleven of 63 directed records are reported under both diseases at once.**
ClinVar submitters routinely attach the whole gene's disease list to a variant,
so a single entry carries both "Dehydrated hereditary stomatocytosis" and
"Lymphatic malformation 6". Those get no direction. Resolving them by preferring
one disease would have manufactured 11 labels out of nothing.

The wild-type gate also did its job: three records disagree with Q92508 at their
stated position and were dropped, and 117 had a protein change that would not
parse unambiguously.

### An independent check that came free

Nine of the ClinVar variants are already in the curated set from
electrophysiology. The condition-based inference agrees with the measured
direction **8 times out of 9** — which is the only evidence available that the
inference is worth anything at all.

The one disagreement is V598M: we have it as gain-of-function, the inference
says loss. Our own record reads *"increased opening (one report); no change in
another"*, so the literature is mixed and neither side is simply wrong. It is
reported rather than resolved, because picking a side by fiat would hide that
from any test built on the set.

### What it bought

The directional missense set goes from **26 (20 GoF, 6 LoF) to 46 (27 GoF,
19 LoF)**. The loss-of-function class more than triples, which is precisely what
Round 22 was starved of, and the design becomes far better balanced.

Minimum detectable effect: **0.61 → 0.41**. Power at a large effect rises from
0.50 to 0.83.

### What it did not buy, stated because the roadmap asked

**A medium effect is still out of reach.** Power at δ = 0.28 is **0.49**, not
0.80, and reaching it at this ratio would need **104** variants against the 46
now available. The constraint has loosened by a useful amount; it has not
lifted.

The original `variants.json` is untouched. Round 7 and Round 22 reference it,
and growing it underneath a recorded result would invalidate that result while
nothing appeared to change.

Suite 462 → 477 passing.

## Round 28 — the footprint in the area change (2026-08-06)

Round 18 built the nonlinear footprint solver and showed the linearised version
is 3.5× too large at PIEZO1's 63° contact slope. Only `DomeModel` consumed it;
the gating energetics still ran on linear numbers.

### The quantity was wrong before the model was

ΔA is a **change between states**, not an absolute area. Round 3's framing —
"the footprint stores 622 nm² against the dome's 256" — invites treating a
stored area as the gating area, and they are different things. What tension does
work on is the *increase* in projected in-plane area on opening.

So both endpoints had to be measured: 7WLT closed at R_c 9.72 nm and contact
slope 1.992 (63.3°), 7WLU flattened at R_c 18.38 nm and slope 0.839 (40.0°).

### The correction is bigger on the difference than on either endpoint

This is the part I did not anticipate. The closed state sits at 63°, where the
small-slope expansion fails badly; the open state at 40°, where it is much less
wrong. Taking a difference between one badly overestimated number and one mildly
overestimated one amplifies the error rather than cancelling it:

- footprint stored closed: 622 → 179 nm² (3.5×, as Round 18 found)
- footprint stored open: 159 → 108 nm² (1.5×)
- **footprint released on opening: 463 → 71 nm² (6.5×)**

Total ΔA falls from 664 to 272 nm².

### The roadmap's question, answered

It asked to report the change *including if the linear version happened to agree
better*, since a wrong model can fit a right number. It did not: T₅₀ moves from
0.060 to 0.147 mN/m, toward the measured 2.7 ± 0.1 rather than away. There is a
test asserting that direction, because the opposite was a live possibility and
would have been the more interesting result.

### What the round actually established

**The correction does not close the gap, and that is worth more than the
correction.** Improving the membrane physics by a factor of six moved T₅₀ by a
factor of 2.4 and left it about **eighteen times below measurement**. The
structural ΔA is still 34× the functional 8 nm².

So the structural-versus-functional discrepancy that `docs/SCIENCE.md` has
carried since Round 3 is **not a membrane-modelling error**. No further
refinement of the footprint will fix it. The two numbers measure different
things: the functional ΔA is the area change along the gating reaction
coordinate, the structural one is the whole protein-plus-footprint deformation,
and only part of that is coupled to the gate.

Knowing that a candidate explanation has been ruled out is a real result, and it
took building the better model to rule it out.

Suite 477 → 489 passing; 19 documented numbers reproduce.

## Round 29 — intervals, and three different things people call one (2026-08-06)

Rounds 18 to 28 kept finding the same failure: a number recorded with more
confidence than it had earned. The footprint area wrong by 3.5×. A null result
that could only ever have excluded a large effect. A T₅₀ eighteen times off. A
point estimate invites that, so this round attaches a spread.

### The distinction that shapes the module

Three things get called "the uncertainty" and they mean different things:

- **Bootstrap** — resample the data. This is a genuine confidence interval for
  sampling variability and the only one that deserves the name.
- **Sensitivity** — vary a method choice, such as the elastic-network cutoff.
  There is no sampling distribution here; a cutoff is not a random variable.
  Quoting its spread as a confidence interval would be a second kind of
  overconfidence dressed as a cure for the first.
- **Parameter range** — vary a registered input over its published values.
  Propagated uncertainty from someone else's measurement, not a statement about
  this dataset.

They are separate classes with different names, and the printed summary of a
sensitivity range says *"sensitivity, not a confidence interval"* in so many
words.

### Two headline numbers change meaning

**The dome radius.** 9.73 nm with a 95% CI of **[8.83, 10.34]**, bootstrapped
over 66 transmembrane surface points. The published 10.2 nm is *inside* that
interval.

The project has been reporting this as agreement-with-a-small-gap — "our 9.7
against their 10.2". The honest statement is that the two are **statistically
indistinguishable**: a stronger claim of consistency, and a weaker claim of
precision than the two-decimal figure implied.

**The gating overlap.** 0.705 at the default cutoff, but **0.554–0.723** across
cutoffs from 10 to 20 Å, non-monotonically. The qualitative conclusion — a
substantial overlap carried entirely by A-symmetric modes — holds at every
cutoff, and the E-mode contribution stays negligible throughout. The third digit
was never meaningful, and it has been quoted as though it were since Round 4.

Ensemble PC1 is 0.900 [0.796, 0.972] over ten structures; dropping 11ZC alone
moves it to 0.832, which is the largest single influence but not a dominating
one.

### Two details worth recording

The confidence level is **derived from `stats.alpha`** rather than written as
0.95. Two copies of a significance level drift, and this project already
registers one — setting α to 0.01 now widens every interval automatically.

The parameter sweep restores the registry in a `finally`, including when the
statistic raises. Leaving it modified would make every later number in the
session incomparable with the documentation, which is exactly what the override
tracking from the parameter round exists to prevent.

### What none of it does

**It does not capture model error.** Bootstrapping a sphere fit tells you how
well the sphere is determined; it says nothing about whether a sphere was the
right shape for the dome, or whether springs are the right physics for a
2500-residue trimer. That is stated on every result rather than left for the
reader to remember, because it is precisely the assumption an interval invites
you to stop questioning.

The parameter audit caught the new module on its first run — seven unregistered
literals — which is the rule from the parameter round working as intended.

Suite 489 → 507 passing.

## Round 30 — checking without reusing the derivation (2026-08-06)

Round 18 is the reason this round exists. The elastica solver integrated
Euler–Lagrange equations derived by hand, and a boundary-value solver converges
happily onto wrong equations. What caught the error was evaluating the exact
functional in a different gauge — a check that reused none of the derivation
being checked. A test written from the same understanding as the code shares its
blind spots; if the author misread the physics, the test encodes the misreading.

Three headline results, each re-derived by a route sharing no machinery with the
original. Two of the three taught me something.

### The dome: a disagreement, and the checking route was the wrong one

A parabola through the radial height profile gave **8.12 nm against the sphere
fit's 9.72** — 16.6% apart. That looks like a finding.

It was, but not the one it appeared to be. On synthetic caps of *known* radius,
the parabola is accurate to 0.6% at an 8.6° contact angle and **25.8% low at
63.4°** — and PIEZO1's dome sits at 63.3°. The parabola `h = h₀ − r²/2R` is a
shallow-cap approximation, which is the *same small-slope failure* Round 18
found in the membrane theory, now turning up in a geometry method. The sphere
fit is exact at every slope on the same synthetics.

So I built a route that is actually valid: the exact cap relation
`R = −(h² + r²)/2h`, inverted per point, no fitting and no expansion. It gives
**10.17 nm against 9.72** — 4.5% apart, both inside Round 29's bootstrap
interval and both consistent with the published 10.2.

Both alternative routes are kept. The invalid one is worth keeping precisely
because knowing *why* a check disagrees is the useful part, and a future reader
who reaches for a parabolic fit should find the answer already written down.

### The overlap: superposition-free, and it holds

Distances between sites are invariant to rotation and translation. Comparing the
observed transition and each normal mode in *pairwise distance changes* uses no
Kabsch fit and no protomer matching — so neither can manufacture the result.
**0.641 against 0.705, 9.0% apart.** Given this project has been bitten twice by
protomer correspondence, that is a check worth having.

### T₅₀: a disagreement that was mine

The first alternative used the analytic steady state and returned zero. I nearly
wrote that up as a discrepancy.

The diagnosis: at equilibrium this channel sits **~96% inactivated at every
tension**, so steady-state open occupancy runs only 0.030 to 0.036 and never has
a half-maximum. The route was computing a *different quantity*. T₅₀ is
necessarily a property of the peak transient — which is also what a patch-clamp
measures, so the pipeline's definition is the experimentally right one.

Replaced with adaptive Runge–Kutta integration of the same master equation: same
quantity, different numerics, no matrix exponential. **2.727 against 2.711 —
0.6%.** The steady-state function is kept, with a test asserting it must *not*
reproduce the peak-based number: if it ever does, one of the two has stopped
computing what it claims.

### A small thing the audit prompted

The parameter audit flagged the new module twice. The second flag was
`OPEN_STATE = 1`, an array index. Rather than exempt it, it now derives from
`STATE_NAMES.index("O")` — a hardcoded index would silently read the wrong
occupancy if the state order ever changed, and that is exactly the class of
error this round is about.

Suite 507 → 518 passing. Block I appended after the thirty-round review,
numbered from 36 since Block H already claims 31–35.
