# Round 22: does substitution-aware evidence separate GoF from LoF?

**Exploratory. Not a validation.** Run 2026-08-06 against
`docs/PREREGISTRATION_ROUND22.md`, which was written and committed before any
comparison with the labels. Reported in the order that document fixed.

**Result in one line: nothing separates the classes, and the primary endpoint's
point estimate runs opposite to the hypothesis with an interval spanning zero.**

---

## 1. Counts achieved

| | n |
|---|---|
| Curated variants | 68 |
| With a directional label | 39 |
| Dropped — not a single-residue substitution | 13 (**11 of them LoF**) |
| **Included** | **26** — 20 GoF, 6 LoF |

Dropped: `R49*`, `Q422*`, `E679*`, `E755*`, `E756-`, `W1069*`, `S1153fs`,
`E1630*`, `K1877-`, `Q2228*`, `E2322*`, `E2496ELE`, `E2498-`.

FoldX covers 25 of the 26; the mechanical ΔΔG joins for 22.

## 2. Primary endpoint — FoldX ΔΔG, hypothesised LoF > GoF

| | value |
|---|---|
| n | 6 LoF vs 19 GoF |
| mean ΔΔG, LoF | **0.767 kcal/mol** |
| mean ΔΔG, GoF | **1.309 kcal/mol** |
| Cliff's delta | **−0.211** (small) |
| bootstrap 95% CI | **[−0.684, +0.298]** — spans zero |
| AUROC (LoF positive) | **0.395** |
| one-sided p, in the hypothesised direction | 0.747 |

**The point estimate is in the opposite direction to H1.** Loss-of-function
missense variants are, in this sample, *less* structurally destabilising than
gain-of-function ones, not more. The confidence interval is wide and includes
zero, so this is an uninformative result rather than a reversal — but it
contains no support whatever for the hypothesis.

Three independent statistics agree on the direction (mean difference, Cliff's
delta and AUROC all point the same way), so this is not a sign error.

## 3. This was exploratory; no decision follows

Per §5 of the pre-registration, **no reject/do-not-reject decision is made and
none may be inferred.** The p-value above is reported for completeness and was
not a decision input. Nothing here may be described as validation, and it is not
evidence that the pipeline works.

## 4. The pre-registered objection, restated against the result

§2 of the pre-registration recorded, *before testing*, that excluding the 11
truncating LoF variants removes precisely the "break the protein" mechanism that
H1 leans on, leaving a LoF subset selected for **not** being truncating.

That objection now looks like the most likely explanation of what was found. The
LoF variants that a missense predictor can score are the structurally mild ones
by construction; the ones that destroy the protein are stop codons and
frameshifts, and they were excluded because no missense predictor can score
them. A predictor of destabilisation was therefore asked about the one subset of
loss-of-function variants that does not act by destabilisation.

Recording this in advance is what stops it being a post-hoc rescue. It was
foreseen, it is consistent with the data, and it does not make the result
positive.

## 5. Secondary family

Corrected together by Benjamini–Hochberg at q < 0.05.

| Predictor | n | Cliff's δ | 95% CI | AUROC | p | q |
|---|---|---|---|---|---|---|
| AlphaMissense | 6+20 | +0.183 | [−0.333, +0.667] | 0.592 | 0.375 | 0.469 |
| EVE | 6+20 | +0.133 | [−0.367, +0.600] | 0.567 | 0.502 | 0.502 |
| ESM-1b | 6+20 | −0.350 | [−0.733, +0.067] | 0.325 | 0.226 | 0.448 |
| Conservation | 6+20 | +0.200 | [−0.250, +0.600] | 0.600 | 0.269 | 0.448 |
| Mechanical ΔΔG | 6+16 | +0.125 | [−0.458, +0.667] | 0.562 | 0.197 | 0.448 |

**Nothing is significant, and nothing is close.** Every interval spans zero; the
smallest adjusted value is q = 0.448.

The pre-registration recorded an expectation of **no separation** for
AlphaMissense, EVE and ESM-1b, on the grounds that a single pathogenicity axis
cannot express direction when both classes are pathogenic. That expectation
held. ESM-1b's point estimate runs opposite to the others, which is what noise
at this sample size looks like.

The mechanical ΔΔG is null again, consistent with Round 7 and with no
information added.

## 6. The combined score

Equal weights, standardised, no fitting, five features, 25 complete rows.

| | AUROC |
|---|---|
| in sample | 0.535 |
| **leave-one-out** | **0.535** |
| optimism | **0.000** |

Combining them recovers nothing that the parts lacked. The zero optimism is
expected rather than reassuring: an equal-weight sum has no fitted parameters,
so there is nothing for cross-validation to penalise. It does confirm that the
0.535 is not an artefact of fitting.

## 7. Power — what this result can and cannot mean

At the achieved 20 versus 6:

| Cliff's δ | Power |
|---|---|
| 0.20 | 0.18 |
| 0.30 | 0.29 |
| 0.43 (large) | **0.52** |
| 0.61 | 0.81 |

The design is a coin toss at a conventionally large effect. **This result rules
out essentially nothing.** It is not evidence of absence; it is close to an
absence of evidence.

## 8. What a fair confirmatory test would require

1. **More phenotyped variants.** 53 at this 20:6 ratio for a large effect, 148
   for a medium one. This is Round 27 and it is the binding constraint on the
   project's central claim — the limitation is data, not method.
2. **A better-posed question for the truncating variants.** Eleven of seventeen
   LoF variants are stop codons or frameshifts. Predicting *those* from
   structure is nearly trivial and probably not worth doing; the interesting and
   hard question is confined to missense, where the sample is six.
3. **A directional predictor, which none of these is.** Four of the five
   secondary features emit a pathogenicity axis with no room for direction, and
   the fifth reports position rather than substitution (Round 7). A predictor
   that could distinguish "opens too easily" from "does not open" would need to
   model the gating equilibrium itself, not the folding free energy.
4. **Class balance.** Six variants in the smaller class means the permutation
   null has few distinct arrangements and every interval is wide, whatever the
   effect.

## 9. Standing record

Round 7's null stands. Round 22 adds a second null, from a different predictor,
with a pre-registered objection that explains it and was recorded in advance.
Neither is revised. The project's central claim remains **untested at adequate
power**, and saying so is more useful than any reanalysis of these 26 variants.
