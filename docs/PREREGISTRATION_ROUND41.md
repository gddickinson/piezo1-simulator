# Pre-registration — Round 41, population constraint

**Status: committed before the comparison was executed.** Written under
`docs/NEGATIVE_RESULT_PROTOCOL.md` §2. Nothing below was written after seeing a
result.

The design numbers — counts, power, exclusions — are properties of the design
and do not depend on the association between constraint and phenotype. §3 of the
protocol requires them beforehand.

---

## 1. Why look at gnomAD at all

Three pre-registered tests have now returned null, and every diagnosis has ended
in the same place: **not enough phenotyped variants**. Round 34 added that the
structural side cannot supply more — one informative variant structure, all
gain-of-function.

Population constraint is the one signal that exists at **every residue** rather
than at 39. It is not a direction, and cannot become one: a constrained position
is one where variation is selected against, which is consistent with either
gain- or loss-of-function being harmful there. So the question this round can
ask is narrower than the project's central claim:

> Do curated gain-of-function and loss-of-function variants sit in
> **differently-constrained positions**?

That is worth asking because if constraint separates them, it is a feature that
exists genome-wide and could be combined with the mechanical score.

---

## 2. The gene-level result, which is already known and is discouraging

This was measured before writing this file and is recorded here rather than
buried, because it changes what a null would mean.

| Metric | PIEZO1 | Interpretation |
|---|---|---|
| LOEUF (`oe_lof_upper`) | **1.10** | LoF variants **not depleted at all** (constrained is < 0.35) |
| pLI | **1.5 × 10⁻¹⁰⁰** | no evidence of haploinsufficiency |
| `oe_mis` | **1.45** | 45% **more** missense than expected |
| `mis_z` | **−11.3** | missense **enriched**, strongly |

**PIEZO1 is not a constrained gene.** The LoF tolerance is biologically
unsurprising — E756del is a common African allele — but the missense *enrichment*
is the part that matters here: gene-wide, selection is not removing missense
variation from PIEZO1, it is over-represented relative to the mutational model.

A regional signal can still exist inside a gene with no global one. But the
prior is worse than the roadmap assumed when it proposed this round, and that is
stated now rather than after the fact.

---

## 3. The hypothesis

> **H₁.** The local missense rate is **lower** (more constrained) at
> loss-of-function positions than at gain-of-function positions.

The direction is reasoned, not free: a loss-of-function variant disrupts a
residue the protein needs, so such positions should be under stronger purifying
selection. Gain-of-function variants act at regulatory or allosteric positions,
which need not be depleted.

**H₀** is that the two do not differ in that direction.

---

## 4. The predictor, frozen

Local missense rate from `analysis.gnomad.MissenseDensity.local_rate()`: the
count of distinct observed missense alleles per residue, averaged in a
**±25-residue window** centred on the position.

The window is fixed now and is not tuned. A single residue carries a handful of
observed variants at most, so a per-residue count is mostly shot noise; the
window is what makes this a *regional* constraint estimate, which is the form the
constraint literature reports as informative.

Data: gnomAD v4 exomes and genomes via the public GraphQL API, cached. 6,708
missense variants placed on 2,521 residues, none unplaced.

---

## 5. Inclusion criteria, fixed now

1. Missense variants with a recorded direction (`GoF` or `LoF`).
2. **One value per position, not per variant.** The predictor is a property of
   the position, so several variants at one residue would enter identical values
   and inflate the effective sample size. Positions are deduplicated.
3. **Position 2456 is excluded.** It is the only position carrying variants of
   *both* directions — R2456H/K/P are gain-of-function, R2456C is loss — and a
   position-level predictor is structurally incapable of separating them. Its
   exclusion is recorded here so it cannot later look like a response to its
   value. It is also the cleanest possible demonstration of this predictor's
   ceiling.
4. Positions outside 1–2521 are excluded; there are none expected.

---

## 6. Evidence levels

Not pooled silently. The primary uses the combined set for power, exactly as
Round 36 did and for the same reason.

| Level | Positions | GoF | LoF | 80% power at |
|---|---|---|---|---|
| `measured` | 22 | 17 | 5 | \|δ\| ≥ 0.64 |
| `measured` + `disease_mechanism` | 42 | 24 | 18 | **\|δ\| ≥ 0.43** |

The measured-only set is a pre-declared secondary.

---

## 7. Endpoints

**Primary — exactly one.** Cliff's δ for local missense rate, LoF versus GoF
positions, on the combined 42-position set, one-sided in the hypothesised
direction (LoF lower), permutation test with the (r+1)/(n+1) convention,
α = 0.05, uncorrected.

**Secondary family**, corrected together by Benjamini–Hochberg at q = 0.05:

1. The same statistic on the measured-only positions.
2. Summed allele count per position rather than distinct alleles — how *often* a
   position is hit rather than how many ways.
3. Window ±10 and ±50, to show whether any result depends on the window this
   pre-registration froze at ±25.
4. Raw per-residue count with no window, as a negative control: it should be
   dominated by shot noise and show nothing.

---

## 8. Decision rule, fixed in advance

- **Reject H₀** only if the primary p < 0.05 **and** Cliff's δ is negative
  **and** its bootstrap 95% interval excludes zero.
- **Fail to reject** otherwise, and the predictor is **not** adjusted and re-run.
  The window is not re-tuned; that is what the secondary window entries are for.
- A null is written up in the order of protocol §6.

---

## 9. Pre-committed caveats

- **Constraint is not direction.** Even a positive result would say the two
  classes sit in differently-constrained regions, not that constraint predicts
  direction for an unseen variant.
- **The gene shows no global constraint** (§2), so any regional signal is being
  sought against an unpromising background.
- **This is a position-level test.** Round 7 died of the position confound, but
  here it is the hypothesis rather than a confound — with the consequence that
  the predictor can never distinguish two variants at one residue, as R2456
  demonstrates.
- **Ascertainment.** gnomAD is not a random sample of humanity; population
  composition affects which variants are observed, and PIEZO1 carries alleles at
  high frequency in specific populations for malaria-related reasons.
- **Eighteen of the 42 directions are inferred**, not measured.
- A null here would be the **fourth**. At that point the honest summary is that
  this project has not found any signal that separates gain- from
  loss-of-function, from structure or from population genetics.
