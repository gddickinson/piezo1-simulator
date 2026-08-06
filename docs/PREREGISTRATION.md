# Pre-registration: does the elastic-network predictor separate GoF from LoF?

**Written 2026-08-06, at the end of Round 6 and before any comparison between
predicted ΔΔG and variant phenotype has been performed.**

This document exists because the temptation it guards against is real and
strong. A predictor that is adjusted until it separates the classes is not a
prediction, it is a fit; and the difference is invisible in the final number.
Fixing the hypothesis, the statistic and the decision rule in advance is the
only thing that makes the Round 7 result mean anything.

Nothing below may be changed after Round 7 runs. If a different model is tried
later it is a **new hypothesis**, reported as such, with the original result
left standing.

---

## 1. Hypothesis

**H1.** Variants that soften the gating coordinate (ΔΔG_gating < 0, the motion
becomes cheaper) are enriched for **gain of function**; variants that stiffen it
(ΔΔG_gating > 0) are enriched for **loss of function**.

**H0.** Predicted ΔΔG_gating is unrelated to phenotype class.

The direction is not arbitrary and is fixed by the physics: PIEZO1 opens by
flattening, so anything that lowers the elastic cost of that motion should
lower the tension threshold.

## 2. The predictor, frozen

- **Model:** `piezo1.analysis.variant_impact.VariantImpactModel`
- **Statistic:** ΔΔG_gating = ½·dᵀ(H_mut − H_wt)·d, in elastic-network energy
  units.
- **Gating coordinate d:** PC1 of the mouse experimental ensemble
  (`build_ensemble(species="mouse", min_common=900)`), which carries 90.0% of
  the variance and matches an A-symmetric elastic mode at 0.804 overlap. Chosen
  in Round 4 on structural grounds, before any variant work.
- **Perturbation:** spring constants at the mutated residue scaled by
  `1 + (V_mut − V_wt)/V_wt`, with Zamyatnin side-chain volumes, sensitivity 1.0.
- **Network:** cutoff 15 Å, `inverse_square` springs, d0 = 7.5 Å, γ = 1.0 —
  the same parameters used for every elastic-network result in this project.
- **Mutation applied to all three protomers.**

No parameter above may be adjusted in Round 7.

## 3. Inclusion criteria

From the 68 curated variants in `piezo1/resources/variants.json`:

- **Include** single-residue substitutions whose residue is resolved in the
  ensemble's shared basis, and whose classification is `GoF` or `LoF`.
- **Exclude** deletions, insertions and multi-residue changes (the volume model
  does not define a spring scale for them).
- **Exclude** `VUS`, `blood-group` and `engineered` classes from the primary
  test. Blood-group variants in particular are documented as antigenic but
  electrophysiologically wild-type, so they are not informative about gating.
- **Exclude** variants whose residue is not in the network, and report how many.

The counts must be reported before the test statistic.

## 4. Primary test

- **Statistic:** difference in mean ΔΔG_gating between the GoF and LoF groups,
  in the direction GoF < LoF (softening for GoF).
- **Test:** one-sided **permutation test**, 10 000 label shuffles, since the
  group sizes are small and unequal and ΔΔG is not expected to be normal.
- **Effect size:** Cliff's delta, with a bootstrap 95% confidence interval
  (10 000 resamples). Cliff's delta rather than Cohen's d because it is
  non-parametric and robust to the outliers a mechanical model will produce.
- **Significance threshold:** p < 0.05, one-sided.
- **Also reported regardless of outcome:** AUROC for ranking GoF above LoF, and
  the same analysis on `ddg_normalised`.

## 5. Decision rule, fixed in advance

| Outcome | Interpretation to be reported |
|---|---|
| p < 0.05 **and** Cliff's delta CI excludes 0, in the predicted direction | H1 supported. Report as a *structure-derived* signal, still requiring prospective testing. |
| p ≥ 0.05 | **H0 not rejected.** Report as a null result. |
| p < 0.05 in the **opposite** direction | Report as a significant result contradicting the physical hypothesis, and say so plainly. |

## 6. Pre-committed caveats

These will be reported whatever the outcome, because they limit the result in
either direction:

1. **Coverage.** Roughly 17 of the 68 variants sit outside the resolved range,
   including E756del. The test speaks only to the subset that is modelled.
2. **Class imbalance and provenance.** GoF variants are mostly xerocytosis
   alleles from a handful of papers; LoF are mostly lymphatic-dysplasia
   alleles. Any separation could reflect *where in the protein each literature
   looked* rather than mechanics. A domain-stratified secondary analysis will
   be reported for this reason.
3. **The volume model is crude.** It captures packing, not charge, hydrogen
   bonding, proline backbone effects or folding stability. A null result is
   therefore weak evidence against the *approach*, and a positive result is not
   evidence that volume is the operative variable.
4. **ΔΔG is in arbitrary elastic units.** Only its sign and relative magnitude
   are interpretable; it is not a kcal/mol.
5. **Not independent of the ensemble.** The gating coordinate comes from the
   same structures the variants are mapped onto. This is a test of internal
   mechanical consistency, not an out-of-sample prediction.

## 7. Secondary analyses, also pre-specified

Reported as exploratory, and **not** used to rescue a null primary result:

- Leave-one-out cross-validation of a single-feature classifier.
- Stratification by domain, to test caveat 2.
- The same test using PRS gate-response and path betweenness in place of ΔΔG.
- Whether the archetypal gain-of-function variant **R2456H** is predicted to
  soften. This is a single case and cannot support a conclusion on its own; it
  is recorded because it is the variant with a solved structure and the
  best-characterised phenotype.

## 8. What will be written

`docs/VALIDATION.md`, containing the counts, the primary statistic, the effect
size with its interval, the decision under §5, every caveat in §6, and the
secondary analyses in §7 — **in that order, whatever the numbers turn out to
be.**
