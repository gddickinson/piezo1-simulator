# Round 47 — what the substitution-aware predictor would need, and what exists

**This is a design analysis, not a test.** No comparison between prediction and
phenotype is run here. Every effect size below is read from
`analysis.prediction_record.VALIDATION_RECORD`, which holds the results already
recorded in `VALIDATION.md`, `VALIDATION_ROUND22.md` and `VALIDATION_ROUND36.md`.
Running a fresh comparison would require its own pre-registration under
`NEGATIVE_RESULT_PROTOCOL.md` §2, and the roadmap item that commissioned this
round said so explicitly.

Reproduce with `python -m piezo1.analysis.feasibility` or `assess()`.

---

## 1. The question

Round 26 made the mechanical predictor sensitive to *which* substitution
occurred — within-position variance rose from 4.9% to 52.5%, the pre-registered
criterion it was built to meet. Round 36 then tested it and got Cliff's
δ = −0.249, p = 0.405: the fourth in a row of nulls.

The question that follows is not "what should we try next" but:

> Is an effect of that size detectable with the variants this project has, or
> with the most it could ever have?

---

## 2. What each recorded effect would require

Sample size for 80% power at the effect each predictor actually produced,
by simulation through `design.sample_size_for` (equal groups):

| Round | Predictor | δ | n had | n needed |
|---|---|---|---|---|
| 7 | elastic-network ΔΔG | −0.083 | 25 | **> 800** |
| 22 | FoldX ΔΔG | −0.211 | 26 | 176 |
| 36 | substitution-aware ΔΔG | **−0.249** | **34** | **134** |

Round 26's improvement did help: the effect grew from −0.083 to −0.249, and the
requirement fell from over 800 variants to 134. The predictor got roughly six
times cheaper to validate. It is still four times more expensive than the data
allows.

---

## 3. What could ever exist

| Source | n |
|---|---|
| Curated variants | 68 |
| — with a direction, missense | 46 |
| — surviving Round 36's modelling gate (74%) | **34** |
| Round 45's literature harvest, fresh candidates | +35 |
| **Optimistic ceiling** | **59** |

The ceiling is generous on purpose. It assumes **every** one of the 35 harvest
candidates is hand-curated with a direction and survives the modelling gate —
where today *none* of them carries a direction at all, and Round 45 found only
two with any electrophysiological measurement behind them.

---

## 4. The answer

Minimum detectable effect at 80% power, holding Round 36's 19:15
gain-to-loss imbalance:

| Scenario | n | MDE | power at \|δ\| = 0.249 |
|---|---|---|---|
| today | 34 (19/15) | 0.465 | 0.32 |
| every directional variant modelled | 46 (26/20) | 0.405 | 0.45 |
| **optimistic ceiling** | **59 (33/26)** | **0.356** | **0.51** |
| required | 134 (75/59) | 0.247 | 0.80 |

> **No reachable dataset can detect the effect this predictor produces.**
> At the ceiling the minimum detectable effect is 0.356 against an observed
> 0.249, and power is 0.51 — a coin flip.

---

## 5. Why this is worth stating

"We need more data" and "the data that could exist is not enough" sound alike
and are not. The first invites another round of curation; the second says that a
fifth pre-registered test on this variant set should not be run **whatever
predictor is put into it**, because a null would be uninformative by
construction and a positive would be under-powered enough to distrust.

That is a claim about the design, and it is falsifiable in a specific way: it
fails the moment a source of directional variants larger than the harvest
appears. The test guarding this document asserts the conclusion still follows
from the numbers, so it will fail rather than quietly become untrue.

---

## 6. Discrepancy with the roadmap item

The roadmap asked this question for "39 directioned variants". The current count
is **46** directional missense variants, of which **34** survived Round 36's
modelling gate. The 39 appears to date from before the ClinVar expansion. The
analysis uses the measured counts, and the conclusion does not depend on which
of the three numbers is used — all of them are far below 134.

---

## 7. What would change it

Not a better predictor. An effect twice as large would still need n = 34 to
reach an MDE of 0.465, and no mechanical predictor this project has built comes
close to δ = 0.47.

Two things would:

1. **A directional variant set several times larger.** Deep mutational scanning
   of PIEZO1 with a gain/loss readout would supply hundreds. None exists.
2. **A within-position design.** Every test so far compares variants *across*
   positions, where Round 7 measured 99.8% of the predictor's variance as
   between-position. R2456 — gain-of-function for H/K/P, loss for C — is a
   position that could discriminate on its own, and Round 41 had to exclude it
   for exactly that reason. A design matched *within* position removes the
   between-position variance rather than fighting it, so it needs far fewer
   variants. The curated set currently contains one such position, which is why
   this is not offered as the next round: it is a curation problem, not a
   modelling one, and it is the specific curation problem worth solving.

   (`analysis.paired_variant` is a different pairing — one variant *structure*
   against the wild-type structures, which is Round 34's n = 1 problem. It does
   not address this.)
