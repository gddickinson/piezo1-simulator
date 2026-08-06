# Round 36 — the third null

**Executed against `docs/PREREGISTRATION_ROUND36.md`, which was committed first,
in its own commit (`af37a82`), before any comparison was run.**

Written in the order §6 of `NEGATIVE_RESULT_PROTOCOL.md` prescribes. This is a
record, not an argument; it is not revised if a later round disagrees with it.

---

## 1. Counts, including every exclusion

| Stage | n |
|---|---|
| Directional missense, both evidence levels | 46 |
| − V598M, pre-excluded (curated and inferred directions disagree) | 45 |
| − not modelled in the mouse reference structure | **34** |
| **Analysis set** | **34** (19 GoF, 15 LoF) |

The eleven dropped variants sit at positions the curved/flat mouse pair does not
resolve. A variant at an unmodelled position has no contacts to perturb, so the
predictor returns nothing for it. This criterion was fixed in §4.3 before the
run.

By evidence level within the analysis set: 22 measured, 12 disease-mechanism.

---

## 2. Primary test result

**Substitution-aware mechanical ΔΔG, GoF versus LoF, combined 34-variant set,
one-sided in the hypothesised direction, permutation test.**

| Statistic | Value |
|---|---|
| Cliff's δ | **−0.249** |
| 95% bootstrap CI | **[−0.628, +0.151]** |
| p (permutation, one-sided) | **0.405** |
| AUROC | 0.625 |
| n | 19 GoF vs 15 LoF |

---

## 3. Decision under the fixed rule

§8 required p < 0.05 **and** δ < 0 **and** the CI excluding zero.

δ is negative — the point estimate is in the hypothesised direction — but
p = 0.405 and the interval spans zero comfortably.

> **FAIL TO REJECT H₀.**

Per §8 and protocol §6, the predictor is **not adjusted and re-run**.

---

## 4. Pre-committed caveats, unedited

Reproduced verbatim from §9 of the pre-registration, written before the result:

- **Twenty of the 46 directions are inferred**, not measured. A positive result
  would need replication on measured labels before it could be called validated.
- **The predictor sees one conformational pair.** `d` is the curved→flat
  displacement of the mouse structures; a variant acting through a motion not in
  that vector is invisible to it, whatever its phenotype.
- **Position dominates.** Even after Round 26, 47.5% of the variance is
  between-position, so the predictor still partly reports *where* a residue sits.
- **All deposited variant structures are gain-of-function** (Round 34), so
  nothing in this test is anchored on an experimentally observed LoF structure.
- **A null here is the third.** Three nulls under three designs would be
  substantive evidence that the mechanical route does not carry the signal — a
  conclusion this project should then state plainly rather than continue to
  work around.

---

## 5. Secondary analyses

All corrected together by Benjamini–Hochberg at q = 0.05, including the ones
that came out uninteresting. None was promoted, removed, or added after the run.

| Endpoint | n | Cliff's δ | p | q | significant |
|---|---|---|---|---|---|
| mechanical, measured labels only | 16/6 | −0.188 | 0.515 | 0.669 | no |
| mechanical, volume-only (negative control) | 19/15 | −0.025 | 0.558 | 0.669 | no |
| AlphaMissense | 16/6 | −0.062 | 0.289 | 0.591 | no |
| EVE | 16/6 | −0.042 | 0.295 | 0.591 | no |
| ESM-1b | 16/6 | **+0.250** | 0.826 | 0.826 | no |
| conservation | 16/6 | −0.062 | 0.197 | 0.591 | no |

**Nothing survives correction.** The smallest q is 0.591.

**One pre-registered endpoint could not be run.** FoldX ΔΔG is absent from the
offline ProtVar cache for **0 of 34** variants. It is recorded here as
untestable rather than dropped, because "could not be run" and "was not
significant" are different statements and §7 forbids removing a test from the
family after the fact.

The external predictors are available only for the 22 measured-label variants;
the ClinVar-derived additions are not in the cache, so those rows are 16/6.

---

## 6. Power statement — what this null is entitled to claim

Pre-registered in §6 of the pre-registration, before the run:

| Effect | Cliff's δ | Power (27 vs 19 design) |
|---|---|---|
| small | −0.11 | 0.17 |
| medium | −0.28 | 0.50 |
| **large** | **−0.43** | **0.84** |
| very large | −0.55 | 0.97 |

The design reaches 80% power at **\|δ\| ≥ 0.41**.

> **This null excludes a large effect. It does not exclude a medium one**, where
> power was 50%, and says essentially nothing about a small one. The observed
> δ = −0.249 sits just below "medium", exactly in the region this design cannot
> resolve.

The realised set was 19 vs 15 rather than the 27 vs 19 the power was computed
for, so the achieved power is somewhat **lower** than the table above.

---

## 7. Post-hoc diagnostic — clearly labelled post-hoc

Nothing in this section was pre-registered. It is a description of the result,
not a test of it.

**The predictor improved, and the improvement is visible.** The
substitution-aware ΔΔG gives δ = −0.249; the volume-only predictor Round 7 used
gives δ = −0.025 on the same 34 variants. That is a tenfold larger effect in the
hypothesised direction, consistent with Round 26 having genuinely made the
predictor sensitive to *which* substitution occurred. Neither is significant.

**Across the three tests the effect has grown monotonically:**

| Round | predictor | n | Cliff's δ | p |
|---|---|---|---|---|
| 7 | elastic-network ΔΔG (volume) | 16/9 | −0.083 | 0.234 |
| 22 | FoldX ΔΔG | 20/6 | −0.211 | — |
| 36 | substitution-aware ΔΔG | 19/15 | **−0.249** | 0.405 |

That is suggestive and it is **not evidence**. Three point estimates drifting in
the predicted direction across designs with 13%, low, and 50% power at the
relevant effect size is what one would expect from either a real medium effect
or from chance. Distinguishing them needs data, not another analysis: at
δ = −0.25, roughly **130 variants** would be required for 80% power, against the
34 available.

**ESM-1b came out with the wrong sign** (+0.250). With q = 0.826 this is
uninformative, and is reported only because §7 requires every test in the family
to appear.

---

## 8. What would constitute a fair next test

Not another predictor on these variants. The three nulls differ in predictor and
agree in outcome, and the diagnostic above says the binding constraint is n.

A fair next test needs **more phenotyped variants with measured directions**.
Round 34 established the structural side cannot supply them — one informative
variant structure, all gain-of-function. Block K §45 proposes the only route
that does not require new experiments: harvesting T50 and inactivation constants
from the supplementary tables of the published mutagenesis literature, behind
the same wild-type-residue gate the curated set already uses.

Until that exists, this project should state the position plainly rather than
keep testing: **the mechanical elastic-network route does not demonstrably
predict gain- versus loss-of-function, and the available data cannot resolve
whether a medium-sized effect is there.**
