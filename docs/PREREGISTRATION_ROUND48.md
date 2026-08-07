# Pre-registration — Round 48, wild-type structural features

**Status: committed before the comparison was executed.** Written under
`docs/NEGATIVE_RESULT_PROTOCOL.md` §2. Nothing below was written after seeing a
result.

The design numbers — counts, exclusions, power — are properties of the design
and do not depend on the association being tested. §3 of the protocol requires
them beforehand, and they are in §5 and §6.

---

## 1. The question, and why it is not the previous four

Rounds 7, 22 and 36 tested *per-variant* predictors: an energy computed for the
specific substitution. Round 41 tested a *positional* predictor from population
genetics. All four returned null.

Round 34 established that the structures cannot supply loss-of-function
examples — one informative variant structure exists and it is gain-of-function.
The roadmap's proposal is to stop waiting for LoF structures and ask whether LoF
positions are distinguishable **in the wild type**:

> Do loss-of-function and gain-of-function variants occur at positions that
> differ in burial, conservation, or coupling to the gate?

This is a different signal from any of the four. It asks nothing about the
substitution and everything about where the position sits in the structure.

---

## 2. The ceiling on what this can ever be, stated first

**Any wild-type positional feature has exactly 0% within-position variance.**
Every variant at a residue receives the same value, by construction.

This is the confound that killed Round 7, in its most extreme form. Round 7's
predictor had 4.9% within-position variance; Round 26 raised it to 52.5%. A
feature computed on the wild-type structure has none at all.

The consequence is concrete and is **not** a caveat to be added afterwards:

- The curated set contains **exactly one** position carrying more than one
  variant — R2456 (H, K, P gain-of-function; C loss-of-function).
- It is therefore the only position that could ever demonstrate a positional
  feature discriminating direction, and it is the one position where such a
  feature is guaranteed to fail.
- It is **excluded** under §5, as it was in Round 41, and its exclusion is
  recorded here so it cannot later look like a response to its value.

So even a clearly positive result would license the statement *"GoF and LoF
variants occur at structurally different positions"* and **not** *"this predicts
the direction of a given substitution"*. The second is the project's actual
goal. This round cannot reach it, and is being run because the first statement
would still be worth knowing and is cheap to test.

---

## 3. The hypotheses

Reasoned rather than free. A loss-of-function variant breaks something the
protein needs; a gain-of-function variant alters how it is regulated.

> **H₁ (primary).** Loss-of-function positions are **more buried** — lower
> relative SASA — than gain-of-function positions.

Buried positions are where a substitution disrupts folding or packing, the
classic route to loss of function. Gain-of-function variants are expected at
regulatory, interfacial or gate-coupled positions, which are more exposed.

**H₀** is that the two classes do not differ in that direction.

---

## 4. The predictor, frozen

`analysis.features.build_feature_table` on **7WLT** (mouse, curved/closed, the
project's standard reference), at its documented defaults: `n_modes=30`,
`cutoff=15.0 Å`, `gate_group="hydrophobic_gate"`, conservation and SASA
included. Values are averaged over the three protomers.

The primary endpoint is the `relative_sasa` column: Shrake–Rupley solvent
accessibility as a fraction of the Gly-X-Gly maximum. Low means buried.

Nothing about this table is tuned for this round; it is the same table Round 27
built and `test_features.py` already guards.

**Numbering.** Variants are curated in human numbering and 7WLT is mouse, so
each position is converted through `core.sequence.human_to_mouse` — never a
constant offset.

---

## 5. Inclusion criteria, fixed now

1. Missense variants with a recorded direction (`GoF` or `LoF`), at evidence
   levels `measured` + `disease_mechanism`.
2. **One value per position, not per variant.** The predictor is a property of
   the position; several variants at one residue would enter identical values
   and inflate the effective sample size.
3. **Position 2456 is excluded** — the only position carrying both directions,
   for the reason given in §2.
4. The position must be **modelled in 7WLT**. Twelve are not, almost all in the
   distal blade the cryo-EM model does not resolve. They are excluded, and the
   count is reported rather than quietly dropped.

---

## 6. Counts and power, fixed now

| Stage | n |
|---|---|
| Directional missense variants | 46 |
| Distinct positions | 43 |
| − position 2456 (both directions) | 42 |
| − not modelled in 7WLT | **30** |

| Level | Positions | GoF | LoF | 80% power at |
|---|---|---|---|---|
| `measured` | 18 | 13 | 5 | \|δ\| ≥ 0.665 |
| `measured` + `disease_mechanism` | **30** | **16** | **14** | **\|δ\| ≥ 0.495** |

**This design is confirmatory only for a very large effect.** Round 47 measured
that no reachable dataset can resolve an effect of the size the mechanical
predictor produces (δ ≈ 0.25); this design needs 0.495, which is larger still.

Accordingly, and following the precedent of Round 22:

> **This round is declared EXPLORATORY for any \|δ\| < 0.495.** A result below
> that threshold, in either direction, is not evidence and will not be reported
> as though it were.

---

## 7. Endpoints

**Primary — exactly one.** Cliff's δ for `relative_sasa`, LoF versus GoF
positions, on the combined 30-position set, one-sided in the hypothesised
direction (LoF lower), permutation test with the (r+1)/(n+1) convention,
α = 0.05, uncorrected, with a bootstrap 95% CI.

**Secondary family**, corrected together by Benjamini–Hochberg at q = 0.05:

1. `conservation` — LoF positions more conserved.
2. `prs_gate_response` — GoF positions more strongly coupled to the gate.
3. `gating_amplitude` — GoF positions move more along the gating mode.
4. `distance_to_gate` — GoF positions closer to the gate.
5. The primary statistic on the `measured`-only positions (13 GoF vs 5 LoF).
6. **`distance_to_axis` as a negative control.** Perpendicular distance from the
   three-fold axis is a bulk geometric coordinate with no directional
   mechanism proposed for it. It should show nothing. If it matches the primary,
   the primary is measuring where variants sit in the molecule rather than
   anything about direction — which is what happened to Round 41's control.

---

## 8. Decision rule, fixed in advance

- **Reject H₀** only if the primary p < 0.05 **and** Cliff's δ is negative
  **and** its bootstrap 95% interval excludes zero.
- **Fail to reject** otherwise. The predictor is **not** swapped for another
  feature and re-run; that is what the secondary family is for.
- Any \|δ\| < 0.495 is reported as exploratory regardless of p (§6).
- A null is written up in the order of protocol §6.

The three-part rule is retained deliberately: in Round 41 it was the interval
clause that caught a p of 0.0477 whose effect interval spanned zero.

---

## 9. Pre-committed caveats

- **Zero within-position variance** (§2). Even a positive result is a statement
  about positions, not about substitutions, and cannot become the project's
  variant-direction predictor.
- **The measured/inferred mix.** Sixteen of the 30 directions are *inferred*
  from which disease the variant causes, not measured by electrophysiology.
- **One structure, one state.** Features come from 7WLT, a closed/curved mouse
  structure. A position's burial or gate coupling may differ in the open state.
- **Ascertainment.** Which positions have curated variants at all reflects what
  has been studied, and gain-of-function PIEZO1 variants have been studied more.
- **Feature correlation.** Burial, contact number and conservation are not
  independent; the BH family is corrected for multiplicity but the endpoints are
  not orthogonal, so the family-wise picture is optimistic.
- **This would be the fifth null.** Round 47 has already recorded that the data
  cannot support the central claim; a null here adds a further signal that does
  not separate the classes, and should not be read as a new diagnosis.
