# Round 41 — the fourth null, and the clause that caught it

**Executed against `docs/PREREGISTRATION_ROUND41.md`, committed first in its own
commit (`eb9ae1f`) before any comparison was run.**

Written in the order §6 of `NEGATIVE_RESULT_PROTOCOL.md` prescribes.

---

## 1. Counts, including every exclusion

| Stage | n |
|---|---|
| Directional missense variants, both evidence levels | 46 |
| Distinct positions | 43 |
| − position 2456, pre-excluded (carries **both** directions) | 42 |
| **Analysis set** | **42 positions** (24 GoF, 18 LoF) |

Positions are deduplicated because the predictor is a property of the position:
several variants at one residue would have entered identical values and inflated
the effective sample size. Both criteria were fixed in §5 before the run.

---

## 2. Primary test result

**Local missense rate (±25 residues), loss-of-function versus gain-of-function
positions, one-sided in the hypothesised direction, permutation test.**

| Statistic | Value |
|---|---|
| Median missense rate, LoF positions | 2.275 /residue |
| Median missense rate, GoF positions | 2.520 /residue |
| Cliff's δ | **−0.269** |
| 95% bootstrap CI | **[−0.595, +0.074]** |
| p (permutation, one-sided) | **0.0477** |
| AUROC | 0.634 |
| n | 18 LoF vs 24 GoF |

---

## 3. Decision under the fixed rule

§8 required **all three**: p < 0.05, δ < 0, and the interval excluding zero.

- p = 0.0477 — **passes**
- δ = −0.269, negative as hypothesised — **passes**
- CI [−0.595, **+0.074**] — **fails**: it contains zero

> **FAIL TO REJECT H₀.**

**This is the clause earning its place.** A rule written as "p < 0.05" alone
would have returned this project's first positive result, on a p-value 0.0023
below threshold, with an effect estimate whose interval comfortably contains no
effect at all. The three-part rule was fixed in advance precisely so that a
marginal p could not be promoted by the act of looking at it.

---

## 4. Pre-committed caveats, unedited

Reproduced from §9 of the pre-registration, written before the result:

- **Constraint is not direction.** Even a positive result would say the two
  classes sit in differently-constrained regions, not that constraint predicts
  direction for an unseen variant.
- **The gene shows no global constraint**, so any regional signal is being
  sought against an unpromising background.
- **This is a position-level test.** The predictor can never distinguish two
  variants at one residue, as R2456 demonstrates.
- **Ascertainment.** gnomAD is not a random sample of humanity, and PIEZO1
  carries alleles at high frequency in specific populations for malaria-related
  reasons.
- **Eighteen of the 42 directions are inferred**, not measured.
- A null here is the **fourth**.

---

## 5. Secondary analyses

All corrected together by Benjamini–Hochberg at q = 0.05.

| Endpoint | n (LoF/GoF) | Cliff's δ | p | q | significant |
|---|---|---|---|---|---|
| same, measured labels only | 5/17 | −0.388 | 0.128 | 0.160 | no |
| summed allele count, ±25 | 18/24 | **+0.250** | 0.716 | 0.716 | no |
| local missense rate, ±10 | 18/24 | −0.231 | 0.0625 | 0.140 | no |
| local missense rate, ±50 | 18/24 | −0.220 | 0.0837 | 0.140 | no |
| **raw per-residue count (negative control)** | 18/24 | **−0.231** | 0.0776 | 0.140 | no |

**Nothing survives correction.** The smallest q is 0.140.

**The negative control is the informative row.** §7.4 pre-registered the raw
per-residue count as a control that "should be dominated by shot noise and show
nothing". It shows **δ = −0.231, p = 0.078** — statistically indistinguishable
from the ±25 predictor the primary endpoint used, and from both alternative
windows.

So the windowing that makes this a *regional* constraint estimate contributes
nothing: the same weak signal is present in the unsmoothed counts. The predictor
cannot be distinguished from its own negative control, which is a stronger reason
to disbelieve the primary than the interval alone.

The summed allele count runs the **opposite** direction (+0.250): how *often* a
position is hit says the reverse of how many *ways* it is hit. With q = 0.716
this is uninformative, and it is reported because §7 requires the whole family.

---

## 6. Power statement — what this null is entitled to claim

Pre-registered in §6, before the run.

| Design | 80% power at |
|---|---|
| measured only, 17 GoF vs 5 LoF | \|δ\| ≥ 0.64 |
| combined, 24 GoF vs 18 LoF | **\|δ\| ≥ 0.43** |

The observed δ = −0.269 sits below that, in the region this design cannot
resolve. **The null excludes a large effect and does not exclude a medium one.**

---

## 7. Post-hoc diagnostic — clearly labelled post-hoc

Nothing here was pre-registered.

**The gene-level result predicted this.** §2 of the pre-registration recorded,
before testing, that PIEZO1 has LOEUF 1.10, pLI ≈ 0 and `oe_mis` 1.45 with
`mis_z` −11.3 — not merely unconstrained but missense-*enriched*. Looking for
regional depletion inside a gene where selection is not removing missense
variation was always a long shot, and the outcome is consistent with that.

**The direction of the point estimate is as hypothesised**, and the medians
differ in the predicted way (LoF 2.275 versus GoF 2.520 missense per residue,
about 10% lower). Whether that is a real 10% or noise is exactly what the
interval cannot say at n = 42.

---

## 8. What would constitute a fair next test

Not a different window; §5 froze it and §7 already showed the answer does not
depend on it.

Two things would change the situation:

1. **More measured directions.** The measured-only subset has 5 LoF positions.
   Every route the project has tried ends here.
2. **A constraint metric with more resolution than observed counts** — gnomAD's
   own regional missense constraint model, which fits expected counts per
   region rather than smoothing observations, would separate "few variants seen"
   from "few variants expected". This analysis conflates them, and that is its
   biggest methodological weakness.

Four pre-registered tests, four nulls, four different predictors — mechanical,
FoldX, substitution-aware mechanical, and now population constraint. The honest
summary is that **this project has found no signal that separates gain- from
loss-of-function**, from structure or from population genetics, and that the
binding constraint has been data every time.
