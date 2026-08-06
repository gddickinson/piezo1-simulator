# Negative-result protocol

**Status: standing policy. Written 2026-08-06 (Round 20), before the Round 22
test it governs.**

This document exists so that "the predictor does not separate gain- from
loss-of-function" is a *finding this project reports*, not a signal to keep
adjusting until it separates them. It was written deliberately after Round 7
returned a null and before Round 22 proposes a new hypothesis, so that the rule
cannot be tuned to a result already seen.

---

## 1. Why this is needed here specifically

The failure mode is not dishonesty; it is drift. A predictor is built, tested,
comes back null, and there is always one more defensible-sounding adjustment —
drop the noisier variants, try the normalised score, use a two-sided test,
restrict to the resolved domains. Each step is individually reasonable. Together
they are a search over analyses, and the p-value at the end of that search means
nothing like what it appears to mean.

The project has one primary scientific claim available to it and 68 curated
variants to test it against. That is a small, finite, non-renewable resource.
Spending it on an unrecorded search would leave nothing to validate against.

---

## 2. What must exist before any test is run

1. **A pre-registration file** in `docs/`, committed *before* the comparison is
   executed, containing: the hypothesis, the frozen predictor, the inclusion
   criteria, the primary statistic, the decision rule, and the pre-committed
   caveats. `docs/PREREGISTRATION.md` is the Round 7 template.
2. **A named primary endpoint.** Exactly one. Everything else is secondary and
   is corrected (§4).
3. **A power statement** (§3). If the design cannot detect the effect being
   hypothesised, that must be known and written down beforehand, because it
   determines what a null will be allowed to mean.
4. **The decision rule fixed in advance**, including what happens on a null.

A test run without these is *exploratory*. Exploratory results may be reported
but must be labelled as such and may never be described as validation.

---

## 3. Power, and what a null is allowed to mean

A null result from an underpowered design excludes nothing. Before testing,
compute the minimum detectable effect with
`piezo1.analysis.design.minimum_detectable_effect`, using the real score
distribution as the resampling pool when one is available — the mechanical ΔΔG
values are heavy tailed, and a difference-in-means test loses power on heavy
tails, so a normal model flatters the design.

**The measured position for this project's available data (Round 20):**

| Design | 80% power reached at |
|---|---|
| 16 GoF vs 9 LoF (Round 7, as run) | **\|δ\| ≥ 0.55** |

Power at conventional effect sizes for that same design:

| Effect | Cliff's δ | Power |
|---|---|---|
| observed in Round 7 | −0.083 | **0.13** |
| small | −0.11 | 0.16 |
| medium | −0.28 | 0.35 |
| large | −0.43 | 0.60 |

And the sample sizes that would be needed, at 80% power with equal groups:

| Effect | Variants needed (total) |
|---|---|
| large (δ = 0.43) | **42** |
| medium (δ = 0.28) | **98** |
| small (δ = 0.11) | ≥ 600 |

**Consequence, and it binds Round 22.** Only 25 of the 68 curated variants
survive the Round 7 inclusion criteria, and even relaxing those cannot reach
~45. So a *confirmatory* test on this variant set is adequately powered for a
large effect and nothing smaller. Round 22 must therefore either

- pre-register a large-effect hypothesis and accept that a null excludes only
  large effects, or
- declare itself exploratory, report effect sizes with intervals rather than
  decisions, and defer confirmation until more phenotyped variants exist.

Choosing between those *after* seeing the p-value is prohibited by this
document.

---

## 4. Multiple comparisons

Round 22 has at least six candidate predictors — mechanical coupling,
conservation, AlphaMissense, EVE, ESM-1b, FoldX — plus combinations. Testing
all and reporting the best is how a null becomes a false positive.

- The **primary endpoint** is tested at its uncorrected p-value and is the only
  input to the decision rule.
- All **secondary** predictors form a family corrected by Benjamini–Hochberg
  (`piezo1.analysis.design.benjamini_hochberg`). FDR rather than Bonferroni
  because the sequence-based predictors are strongly correlated — they read the
  same evolutionary signal — and family-wise correction across correlated tests
  is conservative enough to guarantee another null.
- The correction is applied to **every** test performed, including ones that
  came out uninteresting. Tests are not removed from the family after the fact.

Worked illustration of why this matters, using plausible p-values:

| Predictor | p | q (BH) | significant |
|---|---|---|---|
| mechanical (primary) | 0.234 | 0.281 | no |
| conservation | 0.041 | 0.098 | no |
| AlphaMissense | 0.012 | 0.072 | no |
| EVE | 0.180 | 0.270 | no |
| ESM-1b | 0.049 | 0.098 | no |
| FoldX | 0.310 | 0.310 | no |

Three predictors clear p < 0.05. **None survives correction.** Reporting
AlphaMissense at p = 0.012 without the family would be a false discovery
manufactured by the act of looking six times.

---

## 5. Any fitted combination is cross-validated

A combined score with weights fitted on 25 variants and evaluated on those same
25 measures how well 25 points can be fitted. Report leave-one-out performance
(`piezo1.analysis.design.leave_one_out`) alongside in-sample, and report the
**optimism** — the gap between them — explicitly.

Every step that consults the labels goes inside the fold, including
standardisation and feature selection. Standardising on the full dataset before
cross-validating is the commonest way to leak, and it leaks quietly.

---

## 6. What is written when the result is null

A null is written up **in the pre-registered order**, in `docs/VALIDATION.md`:

1. counts, including every exclusion and why;
2. the primary test result;
3. the decision under the fixed rule;
4. the pre-committed caveats, unedited;
5. secondary analyses, labelled as such;
6. **a power statement** — what effect the design could have detected — so the
   null's scope is explicit;
7. a post-hoc diagnostic of *why*, clearly labelled post-hoc;
8. what would constitute a fair next test.

A recorded result is never revised. If a later round shows the analysis was
wrong, the correction is a new entry that says so and points back — as Round 18
did to Round 3 — leaving the original visible. Deleting a superseded result
removes the evidence of how the conclusion was reached, which is the part worth
keeping.

---

## 7. What counts as a new test

Re-running the same hypothesis on the same variants after changing the
predictor is a **new test** and needs a **new pre-registration**. It does not
amend the old one, and the old result stands as recorded.

This is why Round 7's null is not revisited by Rounds 16–19, all of which added
predictors that could have been thrown at it.

---

## 8. Instruments

| Question | Function |
|---|---|
| Did the groups separate? | `validation.permutation_test` |
| How large is the effect, with an interval? | `validation.bootstrap_cliffs_delta` |
| How well does it rank? | `validation.auroc` |
| Could this design have detected it? | `design.power_curve`, `design.minimum_detectable_effect` |
| How many variants would we need? | `design.sample_size_for` |
| Does it survive having looked N times? | `design.benjamini_hochberg` |
| Does it generalise? | `design.leave_one_out` |
