# Pre-registration — Round 36, the third variant test

**Status: committed before the comparison was executed.** Written under
`docs/NEGATIVE_RESULT_PROTOCOL.md` §2, which requires this file to exist, in its
own commit, before any predictor is compared against any phenotype.

Nothing in this document was written after seeing a result. The design numbers
below (counts, power, minimum detectable effect) are properties of the *design*
and depend on group sizes and the score distribution, not on the association
between scores and labels. The protocol explicitly permits computing them
beforehand (§3) — indeed it requires it.

---

## 1. Why a third test, when two have returned null

Round 7 and Round 22 both returned nulls and both stand unrevised. Re-testing
the same hypothesis with the same predictor on the same data would be the drift
this project's protocol exists to prevent.

**Two things changed, on opposite sides of the question, and neither was chosen
after seeing an outcome:**

- **Round 26 changed the predictor.** The old elastic-network ΔΔG was
  effectively blind to *which* substitution occurred: 99.8% of its variance was
  between-position, so all four R2456 substitutions scored alike although three
  are gain-of-function and one is loss. The substitution-aware model scales each
  contact by properties of the new residue and its partner, raising
  within-position variance from **4.9% to 52.5%**. That was validated in Round
  26 against a pre-registered criterion of 20% and **deliberately contained no
  phenotype comparison**.
- **Round 27 changed the data.** The directional missense set grew from 26 to
  **46** by adding ClinVar variants whose direction is *inferred from the
  disease mechanism*. That is weaker evidence than electrophysiology, and it is
  recorded per variant rather than pooled away.

A predictor that can now distinguish substitutions, against a set that is now
large enough to detect a large effect, is a genuinely new test. It is also the
last one this project can afford: the curated set is finite and non-renewable.

---

## 2. The hypothesis

> **H₁.** The substitution-aware mechanical ΔΔG of gating
> (`analysis.variant_impact` with `analysis.substitution` contact scaling)
> is **lower for gain-of-function than for loss-of-function** missense variants.

The direction is not arbitrary. A gain-of-function variant is hypothesised to
*soften* the elastic cost of the gating motion — making opening cheaper — while
a loss-of-function variant stiffens or disrupts it. That is the same directional
claim Round 7 made, and it is retained precisely so that the two tests remain
comparable.

**H₀** is that the two distributions do not differ in that direction.

---

## 3. The frozen predictor

`VariantImpactModel` with substitution-aware contact scaling, evaluated as
½·dᵀ(H_mut − H_wt)·d — the change in elastic cost of the *observed* gating
motion, with `d` the curved→flat displacement.

Frozen means: the predictor is not adjusted after this file is committed. No
re-tuning of the spring model, the cutoff, the contact scaling, or the
normalisation. If the predictor is changed for any reason, this pre-registration
is void and a new one is required.

---

## 4. Inclusion criteria, fixed now

1. Missense only. Nonsense and frameshift variants are excluded: a truncation
   removes the structure the predictor operates on, so its score is not
   comparable.
2. A direction (`GoF` or `LoF`) must be recorded. Variants of uncertain
   significance, engineered controls and blood-group antigens are excluded.
3. The wild-type residue must be resolved in the reference structure. A variant
   at an unmodelled position has no contacts to perturb.
4. **V598M is excluded.** It is the one variant where the curated (measured)
   direction and the ClinVar-inferred direction disagree — GoF versus LoF. It is
   removed *before* the test rather than resolved by fiat, and this is recorded
   here so the exclusion cannot later look like a response to its score.

---

## 5. The two evidence levels, and how they are handled

The protocol forbids pooling them silently. They are not pooled.

| Level | n | GoF | LoF | What "direction" means |
|---|---|---|---|---|
| `measured` | 26 | 20 | 6 | determined by electrophysiology |
| `measured` + `disease_mechanism` | 46 | 27 | 19 | the second group's direction is *inferred* from which disease the variant causes |

**The primary endpoint uses the combined set (46).** Stated now, with the
reason: the measured-only design is powered only at |δ| ≥ 0.61, which is beyond
"large" and would make a null uninformative. The combined set reaches |δ| ≥ 0.41.

The cost is that 20 of 46 labels are inferred rather than measured, so the
primary endpoint is **exploratory-strength evidence about a measured question**.
That trade is made here, before the result, rather than chosen afterwards to
suit it.

The measured-only set is a **pre-declared secondary** endpoint (§7), so that a
reader can see whether the weaker labels carried the result.

---

## 6. Power, and what a null will be allowed to mean

Computed with `piezo1.analysis.design`, using the real score distribution as the
resampling pool.

| Design | 80% power at | small (0.11) | medium (0.28) | large (0.43) |
|---|---|---|---|---|
| measured only, 20 vs 6 | **\|δ\| ≥ 0.61** | 0.11 | 0.25 | 0.51 |
| combined, 27 vs 19 | **\|δ\| ≥ 0.41** | 0.17 | 0.50 | **0.84** |

Total variants needed at 80% power with equal groups: **40** for a large effect,
**108** for a medium one.

**This test is therefore declared CONFIRMATORY FOR A LARGE EFFECT AND
EXPLORATORY BELOW IT.** A null excludes a large effect (δ ≤ −0.43) at 84% power.
It excludes a medium effect at only 50% power and says essentially nothing about
a small one. That sentence is written here so it cannot be softened later.

---

## 7. Endpoints

**Primary — exactly one.** Cliff's δ for substitution-aware ΔΔG, GoF versus LoF,
on the combined 46-variant set, one-sided in the hypothesised direction
(`alternative="less"`), by permutation test with the (r+1)/(n+1) convention,
α = 0.05, tested at its uncorrected p-value.

**Secondary family**, all corrected together by Benjamini–Hochberg at q = 0.05,
including any that come out uninteresting:

1. The same statistic on the measured-only set (26).
2. AlphaMissense, EVE, ESM-1b and precomputed FoldX ΔΔG from ProtVar.
3. Per-residue conservation.
4. The unmodified (volume-only) ΔΔG from Round 7, as a negative control — it
   returned null before and should again; if it now separates, the difference is
   the data rather than the predictor.

No secondary endpoint may be promoted to primary. No test is removed from the
family after the fact.

---

## 8. Decision rule, fixed in advance

- **Reject H₀** only if the primary p < 0.05 *and* the primary Cliff's δ is
  negative *and* its bootstrap 95% confidence interval excludes zero.
- **Fail to reject** otherwise. On a null, the result is recorded as a third
  null, `docs/VALIDATION_ROUND36.md` is written in the order prescribed by
  §6 of the protocol, and **the predictor is not adjusted and re-run.**
- Any fitted combination of predictors is reported with leave-one-out
  performance and its optimism, with every label-consuming step inside the fold.

---

## 9. Pre-committed caveats

Written before the result so they cannot be selected to suit it.

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

## 10. What was found while writing this

`design.sample_size_for` returned `max_n` for any **positive** effect size.
`power_curve` defaults to `alternative="less"`, so a positive δ injected the
effect against the alternative, gave ~zero power at every size, and the
bisection walked to its ceiling — which the docstring presented as "not reached
within the search range", i.e. as a statement about the design. A caller asking
the obvious question, "how many variants would a large effect need?", got a
confident wrong answer that looked like a finding.

Fixed to use the magnitude and take the sign from `alternative`; the answers are
now symmetric in the sign and reproduce the protocol's recorded table (40 for a
large effect against the recorded 42, 108 for medium against 98, the differences
being simulation noise).

That this surfaced while *writing a power section* rather than while reading one
is the argument for the protocol requiring the section at all.
