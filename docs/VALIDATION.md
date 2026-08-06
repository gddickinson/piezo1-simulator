# Blind validation: does the elastic-network predictor separate GoF from LoF?

**Result: it does not. H0 is not rejected.**

Run 2026-08-06 following `docs/PREREGISTRATION.md`, which fixed the hypothesis,
the predictor, the statistic and the decision rule before any comparison was
made. Reported in the order that document specifies, whatever the numbers.

---

## 1. Counts

Gating coordinate: PC1 of the mouse experimental ensemble, 90.0% of variance,
over 10 structures and 1091 residues per protomer.

| | n |
|---|---|
| Curated variants | 68 |
| Dropped — class `engineered` | 15 |
| Dropped — residue not in the resolved network | 12 |
| Dropped — class `VUS` | 8 |
| Dropped — class `blood-group` | 6 |
| Dropped — not a single-residue substitution | 2 |
| **Included** | **25** (16 GoF, 9 LoF) |

## 2. Primary test

| | GoF | LoF |
|---|---|---|
| n | 16 | 9 |
| mean ΔΔG_gating | −4.11 × 10⁻⁶ | +2.20 × 10⁻⁵ |
| median ΔΔG_gating | see `data/derived/validation_round7.json` | |

- **Difference in means:** −2.61 × 10⁻⁵ (GoF lower, i.e. *in the predicted
  direction*)
- **One-sided permutation test**, 10 000 shuffles: **p = 0.234**
- **Cliff's delta:** **−0.083**, bootstrap 95% CI **[−0.528, +0.403]** —
  *negligible*, and the interval spans zero
- **AUROC:** **0.542** (0.5 is chance)

## 3. Decision

Under §5 of the pre-registration: **p ≥ 0.05 → H0 not rejected.**

The mean difference points the right way, but the effect is negligible, its
confidence interval comfortably includes zero, and the ranking performance is
barely distinguishable from a coin toss. This is a null result and is reported
as one.

## 4. Pre-committed caveats

All four were written down before the test and all four apply.

1. **Coverage.** 12 of the 68 variants sit outside the resolved range,
   including E756del. The test speaks only to the 25 that survived inclusion.
2. **Class imbalance and provenance.** 16 versus 9 is a small, unequal split.
   The domain distributions do overlap — GoF spread across anchor (4), CTD (4),
   THU9 (3), cap (2); LoF across cap (3), THU9 (2), CTD (2), anchor (1) — so
   this is not a case of the two literatures having looked in entirely
   different places, but the numbers per domain are far too small to stratify
   meaningfully.
3. **The volume model is crude.** It captures packing and nothing else: no
   charge, no hydrogen bonding, no proline backbone effect, no folding
   stability. A null result is therefore weak evidence against the *approach*.
4. **ΔΔG is in arbitrary elastic units.** Only sign and relative magnitude are
   interpretable.
5. **Not out-of-sample.** The gating coordinate derives from the same
   structures the variants are mapped onto.

## 5. Secondary analyses

Pre-specified, and — as the pre-registration requires — **not** used to rescue
the primary result.

- **On `ddg_normalised`:** p = 0.054, Cliff's delta −0.111 (negligible),
  AUROC 0.556. Closer to threshold but still not significant, and the effect
  remains negligible.
- **R2456, the archetypal position.** Four included variants sit here, and the
  model predicts *softening* for all four — including **R2456C, which is
  loss-of-function**:

  | Variant | Class | ΔΔG | Predicted |
  |---|---|---|---|
  | R2456H | GoF | −2.61 × 10⁻⁶ | softening |
  | R2456K | GoF | −6.20 × 10⁻⁷ | softening |
  | R2456P | GoF | −7.84 × 10⁻⁶ | softening |
  | R2456C | **LoF** | −8.38 × 10⁻⁶ | softening |

  The largest predicted softening belongs to the one loss-of-function variant.

## 6. Why it fails — post-hoc diagnostic

**Labelled post hoc. It played no part in the decision above and is recorded
because it is the useful part of a null result.**

Partitioning the ΔΔG variance:

| | |
|---|---|
| Total variance | 5.42 × 10⁻⁹ |
| Within-position variance | 1.11 × 10⁻¹¹ |
| **Within-position share** | **0.2%** |

**99.8% of the signal is between-position.** The predictor is, in effect,
reporting *where a residue sits in the structure* rather than *which
substitution occurred there*. The R2456 series makes this concrete: four
different substitutions at one position, with phenotypes spanning GoF and LoF,
all receive nearly the same answer.

That is a structural property of the method, not a bug. ΔΔG = ½dᵀ(H_mut−H_wt)d
scales with the local strain of the gating coordinate at that residue and with
its contact count, both of which are properties of the *position*. The
substitution enters only through a single scalar spring multiplier, which is a
far weaker lever than position. Any predictor built this way will be dominated
by position, and positions do not have phenotypes — variants do.

The GoF and LoF ΔΔG ranges overlap completely (GoF [−1.95 × 10⁻⁴, +2.20 × 10⁻⁴],
LoF [−8.38 × 10⁻⁶, +1.70 × 10⁻⁴]).

## 6b. Power — added Round 20, qualifying the scope of this null

**This section does not revise the result above.** H0 was not rejected and that
stands. It records what the design was capable of detecting, which determines
what the null is entitled to claim.

Simulating the pre-registered test at the design's actual group sizes —
16 versus 9, one-sided, α = 0.05 — using both a normal model and resampling
from the observed heavy-tailed ΔΔG pool, which agree:

| Effect | Cliff's δ | Power |
|---|---|---|
| **observed here** | −0.083 | **0.13** |
| small | −0.11 | 0.16 |
| medium | −0.28 | 0.35 |
| large | −0.43 | 0.60 |

**80% power is reached only at \|δ\| ≥ 0.55**, which is past 'large' on the
usual thresholds.

So this null **excludes a large mechanical effect and is close to uninformative
about a small or medium one**. §7 below should be read with that scope: the
correct statement is not "there is no effect" but "there is no large effect,
and this design could not have seen a smaller one". The §6 diagnostic — that
99.8% of the ΔΔG variance is between-position — remains the mechanistic
explanation, and it is independent of the power argument. Both are true, and
the power limit was the one not stated at the time.

Sample sizes that would be needed at 80% power, equal groups: **42** variants
total for a large effect, **98** for medium, ≥600 for small. Only 25 of the 68
curated variants survive the inclusion criteria, so a confirmatory test of
anything below a large effect is not available from this variant set. That
constraint is now binding policy in `docs/NEGATIVE_RESULT_PROTOCOL.md`.

Leave-one-out cross-validation of the two Round 7 predictors gives AUROC 0.535
out-of-sample against 0.542 in-sample (optimism +0.007) — consistent with the
null rather than additional evidence for it, and confirming there was no
hidden overfitting inflating the original number.

## 7. What this does and does not mean

**It does not** invalidate the physics chain. Every link in it was validated
against an independent published number: dome curvature 9.7 nm against 10.2;
T₅₀ 4.99 mN/m against 5.1 ± 0.2; footprint decay length recovered as 13.998 nm
against 14.0; PC1 of the experimental ensemble matching an A-symmetric elastic
mode at 0.804 overlap. Those results stand.

**It does** mean that a single scalar derived from a volume-scaled elastic
network is not sufficient to call a variant's phenotype, on this dataset at
this sample size. Stated plainly: the mechanical model knows about the
*machine*, and this is the wrong instrument for asking about the *substitution*.

## 8. What would be a fair next test

Recorded now so that whatever is tried next is a stated new hypothesis rather
than a retrofit:

1. **A substitution-aware perturbation.** Charge, hydrogen-bonding capacity and
   proline backbone effects, so that different mutations at one position can
   differ. This is the direct fix for the diagnostic in §6.
2. **Combine with sequence-based predictors.** Round 17 brings in
   AlphaMissense, EVE, ESM-1b and precomputed FoldX ΔΔG through the ProtVar API
   under CC BY 4.0. Those are substitution-aware by construction; the
   mechanical score would then contribute what they lack — the *mechanism*.
3. **Ask a better-posed question.** "Which residues are mechanically coupled to
   the gate" is a question the elastic network can answer well (Round 5 showed
   the anchor is the transmission hub). "Is this particular amino-acid swap
   GoF or LoF" may simply be the wrong question to put to it.
4. **More variants.** n = 25 gives little power for a small effect — Round 20
   now quantifies that: 0.13 at the observed effect, and 98 variants would be
   needed for a medium one. A pre-registered replication on an expanded,
   independently curated set would be worth more than any reanalysis of these
   25, and no reanalysis of these 25 can settle the question either way.
