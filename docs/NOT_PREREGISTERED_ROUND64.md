# Round 64 — a decision not to run the within-position test

**This is not a pre-registration. It is a record of declining to write one**,
and of why, so that the question is not reopened by someone who has not seen
the numbers.

Round 64 was written conditionally: *"Only if Rounds 61–63 leave a design with
adequate power."* They do not. This file exists because "we decided not to test"
is itself a result that should be recorded with its reasoning, in the same way
a null is — otherwise the next person to look at the one discriminating position
simply runs the comparison and reports whatever comes out.

---

## 1. The design that was proposed

Compare two variants at the **same** residue: does the mechanical coupling score
rank the gain-of-function one above the loss-of-function one? The statistic is a
sign test, so nothing has to be assumed about a distribution there is no sample
to estimate.

The motivation was sound and remains sound. Round 7 measured that **99.8%** of
the predictor's variance is between positions, so an across-position comparison
spends almost all of its variance on which residue is being looked at rather
than on which substitution occurred. Pairing removes exactly that term.

---

## 2. Why it will not be run

**Round 61 — what the design needs.**

| Paired δ | Shared positions required |
|---|---|
| 0.25 (the across-position effect) | 102 |
| 0.50 | 26 |
| 0.80 — the predictor ordering nine pairs in ten correctly | **8** |

**Round 62 — what exists.**

| Evidence level | Usable now | Reachable by curating three named variants |
|---|---|---|
| `measured` (electrophysiology) | **1** | 3 |
| `measured` + `disease_mechanism` | 1 | 4 |

**Round 63 — what the held-back evidence adds.** Zero. Five of the fifteen
engineered variants are on the right axis, but none sits beside a directional
curated variant, and the only engineered pair (S1335A/S1335V) is
same-direction.

So the best case is **8 required against 1 available**, or 3–4 after curation
that Round 54 costed and whose yield is an upper bound. There is no effect size
at which this design is powered by the data that exists.

---

## 3. Why not run it anyway and label it exploratory

§2 of `NEGATIVE_RESULT_PROTOCOL.md` permits exploratory work, provided it is
labelled and never described as validation. That permission is not a licence
here, for three reasons.

**A sign test on n = 1 has no useful outcome.** The minimum one-sided p from a
single correctly-ordered pair is 0.5. The test cannot reject at any α a reader
would accept, so the only possible results are "not significant" and "not
significant".

**The one available position is the one the project has already looked at.**
R2456 has been named in every round since Round 7 as the example that breaks the
predictor — four substitutions, three gain-of-function and one loss, scoring
nearly alike. Testing on it is not a blind test; the answer is already known
informally, and a pre-registration written afterwards would be a formality
dressed as a safeguard.

**It would produce a number that outlives its caveat.** Five rounds of this
project's history are numbers that travelled further than the sentence attached
to them. A δ from one pair would be quoted.

---

## 4. What would change this

Not a better predictor, and not a different statistic. Only more discriminating
positions — residues carrying two variants of *opposite* direction, both
missense, both at the same evidence level.

Eight is the threshold at an implausibly good predictor; twenty-six at a
plausible one. The route reopens if a deep mutational scan with a directional
readout is published, or if targeted electrophysiology on the three variants
Round 54 named (M870V, R1358C, A2020V) is done and the results happen to
disagree with their partners — which is two positions at best, since M870's
partner is inferred rather than measured.

---

## 5. Standing instruction

**Do not run a within-position comparison on this variant set.** If a later
round wants to, it must first show that the number of discriminating positions
has reached the requirement in §2 — and `test_not_preregistered_round64.py`
fails if that count grows, so the question resurfaces automatically rather than
depending on someone remembering this file.
