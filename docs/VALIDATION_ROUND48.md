# Round 48 — the fifth null, and the flattest of them

**Executed against `docs/PREREGISTRATION_ROUND48.md`, committed first in its own
commit (`7ffb008`) before any comparison was run.**

Written in the order §6 of `NEGATIVE_RESULT_PROTOCOL.md` prescribes.
Reproduce with `python scripts/run_validation_round48.py`.

---

## 1. Counts, including every exclusion

| Stage | n |
|---|---|
| Directional missense variants, both evidence levels | 46 |
| Distinct positions | 43 |
| − position 2456, pre-excluded (carries **both** directions) | 42 |
| − not modelled in 7WLT | −12 |
| **Analysis set** | **30 positions** (16 GoF, 14 LoF) |

The twelve unmodelled positions are 46, 214, 276, 338, 445, 549, 560, 598, 669,
737, 782 and 1988 — almost all in the distal blade the cryo-EM model does not
resolve, which is the same coverage limit Round 36 hit.

---

## 2. Primary test result

**Relative SASA, loss-of-function versus gain-of-function positions, one-sided
in the hypothesised direction (LoF more buried), permutation test.**

| Statistic | Value |
|---|---|
| Median relative SASA, LoF positions | 0.2190 |
| Median relative SASA, GoF positions | 0.1431 |
| Cliff's δ | **+0.036** |
| 95% bootstrap CI | **[−0.384, +0.473]** |
| p (permutation, one-sided) | **0.5092** |
| AUROC | 0.482 |
| n | 14 LoF vs 16 GoF |

---

## 3. Decision under the fixed rule

§8 required **all three**: p < 0.05, δ < 0, and the interval excluding zero.

- p = 0.5092 — **fails**
- δ = +0.036, *positive*, opposite to the hypothesis — **fails**
- CI [−0.384, +0.473] — **fails**: it contains zero

> **FAIL TO REJECT H₀.**

And under §6, \|δ\| = 0.036 against the 0.495 confirmatory threshold:

> **EXPLORATORY.** This result is not evidence in either direction.

**A note on reading the medians rather than the statistic.** The median LoF
position is *more* exposed than the median GoF position (0.219 vs 0.143), which
looks like a reversal of the hypothesis. Cliff's δ is +0.036 — a rank statistic
counting pairwise dominance — so the two distributions overlap almost
completely and the median gap is carried by a few positions. Anyone reporting
the medians alone would describe a reversal that the effect size says is not
there. That is why the pre-registered statistic was fixed in advance.

---

## 4. The pre-committed ceiling, now measured

§2 stated before the run that any wild-type positional feature has zero
within-position variance. Measured over all 34 modelled directional variants at
31 positions:

| Quantity | Value |
|---|---|
| Between-position variance share | **1.000000** |
| Within-position variance share | **0.000000** |

The demonstration is position 2456, which carries four curated variants:

| Variant | Direction | Wild-type relative SASA |
|---|---|---|
| R2456H | GoF | 0.127326 |
| R2456K | GoF | 0.127326 |
| R2456C | **LoF** | 0.127326 |
| R2456P | GoF | 0.127326 |

Identical to every digit, because the feature never saw the substitution. Set
against the project's other predictors:

| Predictor | Within-position variance |
|---|---|
| Round 7, volume-based ΔΔG | 4.9% |
| Round 26, substitution-aware ΔΔG | 52.5% |
| **Round 48, wild-type positional features** | **0.0%** |

So a positive result here could never have become a variant-direction
predictor. That was stated in the pre-registration, and it is now a measured
number rather than an argument.

---

## 5. Secondary analyses

All corrected together by Benjamini–Hochberg at q = 0.05.

| Endpoint | n (LoF/GoF) | Cliff's δ | p | q | significant |
|---|---|---|---|---|---|
| conservation | 13/16 | −0.062 | 0.7532 | 0.9300 | no |
| PRS gate response | 14/16 | −0.062 | 0.6929 | 0.9300 | no |
| gating-mode amplitude | 14/16 | +0.152 | 0.9300 | 0.9300 | no |
| distance to gate | 14/16 | **+0.000** | 0.5890 | 0.9300 | no |
| **distance to axis (negative control)** | 14/16 | **+0.268** | 0.2259 | 0.9300 | no |
| relative SASA, measured labels only | 5/13 | +0.292 | 0.8410 | 0.9300 | no |

**Nothing survives correction.** The smallest q is 0.9300 — the least
significant secondary family this project has produced.

**The negative control is again the informative row.** §7.6 pre-registered
distance from the three-fold axis as a bulk geometric coordinate with no
proposed mechanism, which "should show nothing". At δ = +0.268 it has a **larger
effect than every pre-registered real endpoint**. The only endpoint above it is
the measured-only subgroup at +0.292, which rests on 5 LoF positions.

The reading is the same one Round 41 forced: at this sample size the spread of
effect sizes across endpoints is dominated by noise, and a feature chosen
*because* it should be meaningless lands in the middle of it. Any single
endpoint that had come out large would not have been distinguishable from this.

**Distance to gate is exactly zero.** δ = +0.000 to three decimals — the two
classes are perfectly interleaved by how far they sit from the pore gate. Of all
the endpoints this is the one with the clearest mechanical story, and it is the
one that separates least.

---

## 6. Power statement — what this null is entitled to claim

Pre-registered in §6, before the run.

| Design | 80% power at |
|---|---|
| measured only, 13 GoF vs 5 LoF | \|δ\| ≥ 0.665 |
| combined, 16 GoF vs 14 LoF | **\|δ\| ≥ 0.495** |

The observed δ = +0.036 is far below that. **This null excludes a very large
effect and excludes nothing else.** It is consistent with a true effect of
δ = 0.3 that the design cannot see — and equally consistent with no effect.

Round 47 measured that no reachable dataset resolves δ ≈ 0.25. This design needs
0.495, so the same conclusion applies with more room to spare.

---

## 7. Post-hoc diagnostic — clearly labelled post-hoc

Nothing here was pre-registered.

**The features are not weak; they are unrelated to direction.** These same
columns do real work elsewhere in the project — conservation recovers the
pore-versus-blade gradient, PRS gate response falls off with distance to the
gate, and `test_features.py` guards both. They are informative about the
structure. They carry no information about which way a variant breaks it.

**The ascertainment worry runs the wrong way to rescue this.** If GoF variants
had been studied preferentially at, say, gate-adjacent positions, that bias
would *create* separation rather than hide it. The endpoints closest to that
story — distance to gate and PRS gate response — are the two flattest.

---

## 8. What would constitute a fair next test

Not another feature on the same 30 positions. Six endpoints have now been
spent, the negative control matched them all, and §8 forbids swapping the
predictor and re-running.

Two things would change the situation, and neither is a modelling change:

1. **Variants at shared positions.** The entire curated set contains one
   position with more than one variant. A design matched within position is the
   only one that escapes the 0% ceiling measured in §4, and it needs curation,
   not code.
2. **Loss-of-function structures**, which Round 34 established do not exist —
   one informative variant entry, gain-of-function.

**Five pre-registered tests, five nulls**, across five distinct predictor
families: elastic-network ΔΔG, FoldX ΔΔG, substitution-aware ΔΔG, population
constraint, and now wild-type structural context. The honest summary is
unchanged and now better supported: **this project has found no signal that
separates gain- from loss-of-function**, and Round 47 established that the data
which could exist would not settle it either.
