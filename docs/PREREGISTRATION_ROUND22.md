# Pre-registration II: does substitution-aware evidence separate GoF from LoF?

**Written and committed 2026-08-06, before any comparison against the phenotype
labels was run.** Round 7's null result stands as recorded and is not revised by
anything here; this is a *new* test of a *new* predictor, as required by
`docs/NEGATIVE_RESULT_PROTOCOL.md` §7.

Everything below was fixed while blind to outcome. The counts and coverage in §3
were obtained by asking which variants have data, which does not involve the
labels. The power analysis in §5 depends only on group sizes. No score has been
compared with any phenotype at the time of writing.

---

## 1. Why a second test at all

Round 7 tested an elastic-network ΔΔG and returned a null. Its diagnostic was
precise and is the reason to try again rather than to stop: **99.8% of that
predictor's variance was between-position**, meaning it reported *where a residue
sits* rather than *which substitution occurred*. Four substitutions at R2456 all
scored "softening" although three are gain-of-function and one is
loss-of-function.

The obvious repair is a substitution-aware predictor. Round 17 brought in four,
through the ProtVar API under CC BY 4.0: AlphaMissense, EVE, ESM-1b and
precomputed FoldX ΔΔG.

But Round 17 also recorded why three of those four cannot, by construction,
answer this question. AlphaMissense, EVE and ESM-1b each emit a **single
pathogenicity axis**, benign to damaging. A variant that opens a channel too
easily and one that stops it opening are both damaging. Observed concretely: all
four R2456 substitutions score PATHOGENIC.

That leaves one predictor with a *directional* mechanism available to it.

## 2. Hypothesis

**H1 — Loss-of-function missense variants are more structurally destabilising
than gain-of-function ones.**

The reasoning, stated before testing: loss of function can be achieved by
breaking the protein — misfolding, mistrafficking, gross destabilisation.
Gain of function cannot. A channel that opens too readily must still fold,
reach the membrane and gate; it is a *retuning*, and retunings are constrained
to be structurally tolerable. So the distribution of destabilisation should sit
higher for LoF than for GoF.

**H0** — FoldX ΔΔG is distributed identically in the two classes.

### The objection to H1, recorded now rather than after the result

The 11 truncating LoF variants (nonsense, frameshift, in-frame deletion) are
excluded by §3 because no missense predictor can score them. Those are precisely
the "break the protein" mechanism. The 6 LoF variants that survive are, *by
selection*, the ones that are **not** truncating — so the analysed subset is
depleted of exactly the cases H1 predicts most strongly.

This is a real weakness and it is not fixable within this dataset. It is
recorded here so that a null cannot later be presented as a clean refutation of
H1, nor a positive result as a stronger confirmation than the design allows.

## 3. Inclusion criteria — fixed, and already applied blind

From the 68 curated variants:

| | n |
|---|---|
| Curated variants | 68 |
| Have a directional label (GoF or LoF) | 39 |
| — dropped: not a single-residue substitution | 13 |
| **Included: single missense with a directional label** | **26** (20 GoF, 6 LoF) |

The 13 dropped are `R49*`, `Q422*`, `E679*`, `E755*`, `E756-`, `W1069*`,
`S1153fs`, `E1630*`, `K1877-`, `Q2228*`, `E2322*`, `E2496ELE`, `E2498-` —
**11 of them loss-of-function**. Engineered mutations, variants of uncertain
significance and blood-group antigens are excluded as in Round 7, since none
carries a directional phenotype.

Coverage within those 26, from the cached ProtVar responses:

| Predictor | GoF | LoF | total |
|---|---|---|---|
| FoldX ΔΔG | 19 | 6 | **25** |
| AlphaMissense | 20 | 6 | 26 |
| EVE | 20 | 6 | 26 |
| ESM-1b | 20 | 6 | 26 |
| Conservation | 20 | 6 | 26 |

A variant missing the predictor under test is dropped from that test only, and
the achieved n is reported with every result.

## 4. The primary endpoint — one, named in advance

**FoldX ΔΔG**, as served by ProtVar, keyed by `mutatedType`.

Chosen over the alternatives for reasons that are all pre-specified:

- it is **substitution-aware**, which is the specific defect Round 7 diagnosed;
- it is **directional by mechanism** (§2), where a pathogenicity axis is not;
- it involves **no fitting**, so there are no researcher degrees of freedom
  between the data and the number;
- it is **external** — computed by someone else, from a structure, with no
  input from this project's models.

Statistic: **difference in means, LoF minus GoF**, tested by a one-sided
permutation test with 10 000 shuffles under the `(r+1)/(n+1)` convention.
Direction: **LoF > GoF**. Effect size: **Cliff's delta** with a percentile
bootstrap 95% CI. Ranking: **AUROC** with LoF as the positive class.

## 5. Power, and the consequence for what may be concluded

Simulated at the achieved group sizes, 20 versus 6, one-sided at α = 0.05:

| Effect | Cliff's δ | Power |
|---|---|---|
| small | 0.20 | 0.19 |
| medium | 0.30 | 0.32 |
| large | 0.43 | **0.52** |
| very large | 0.55 | 0.72 |
| — | 0.61 | 0.80 |

**80% power is reached only at \|δ\| ≥ 0.61**, which is *worse* than Round 7's
0.55 despite one more variant, because the smaller group has shrunk from 9 to 6.
Power is barely a coin toss even at a conventionally large effect.

### Declaration required by `NEGATIVE_RESULT_PROTOCOL.md` §3

**This round is EXPLORATORY.** It is not a confirmatory test and its result will
not be reported as a decision.

The protocol offers a confirmatory alternative — pre-register a large-effect
hypothesis and accept that a null excludes only large effects. That is declined
here, because at this n a "confirmatory null" would exclude only effects beyond
large, and describing that as confirmation of anything would overstate it.

Consequences, binding:

- **Effect sizes with confidence intervals are the reported output**, not
  reject/do-not-reject.
- p-values are reported for completeness and are **not** decision inputs.
- Whatever the result, it may **not** be described as validation, and may not be
  cited as evidence that the pipeline works.
- A positive result here licenses one thing only: a properly powered
  confirmatory test on an expanded variant set (Round 27), pre-registered
  separately.

To reach 80% power at a large effect at this 20:6 ratio would need **53**
variants; at a medium effect, **148**.

## 6. Secondary analyses, also fixed in advance

A family of five, corrected together by Benjamini–Hochberg at q < 0.05. The
correction covers **every** test listed, including those expected to be null;
none may be removed from the family after the fact.

| Predictor | Direction tested | Prior expectation, stated now |
|---|---|---|
| AlphaMissense | two-sided | **no separation** — single pathogenicity axis |
| EVE | two-sided | **no separation** — same reason |
| ESM-1b | two-sided | **no separation** — same reason |
| Conservation (ProtVar) | two-sided | no directional prediction |
| Mechanical ΔΔG (Round 7) | one-sided, GoF lower | null, as already recorded |

Recording the expectation of *no separation* for three of these matters: if they
do separate, that is a surprise requiring explanation, not a success.

## 7. The combined score

The roadmap asks whether substitution-aware evidence combined with this
project's mechanical and conservation features does better than either alone.

Pre-specified: features are standardised and summed with **equal weights and no
fitting**, sign-oriented so that each is expected positive for LoF under §2.
Evaluated by **leave-one-out AUROC**, with in-sample AUROC and the optimism
between them both reported.

An equal-weight combination is chosen deliberately over a fitted one. Fitting
six weights on 26 points measures how well 26 points can be fitted; with 6
variants in the smaller class, a fitted model has roughly four points per
parameter and its in-sample performance is uninterpretable.

## 8. What will be written, regardless of outcome

In `docs/VALIDATION_ROUND22.md`, in this order:

1. counts achieved, with every exclusion;
2. the primary effect size and its interval, with the p-value alongside;
3. the explicit statement that this was exploratory and no decision follows;
4. the §2 objection restated against whatever was found;
5. the secondary family with BH-adjusted values, including the null ones;
6. the combined score in and out of sample, with the optimism;
7. a power statement bounding what the result can mean;
8. what a fair confirmatory test would require.

No result will be removed, and no analysis not listed here will be added without
being labelled post-hoc.
