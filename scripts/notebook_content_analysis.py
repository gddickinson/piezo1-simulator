"""Cells for the last two notebooks: the conducting pathway, and the variants.

Split from `notebook_content.py` at the seam between "how the channel works"
and "what we tried to predict with it" — and to keep both files well under the
project's length limit.

The fourth notebook is the awkward one and is meant to be. It walks the variant
workflow to the point where it stops working, because a reader who takes the
machinery and not the result would be repeating a mistake this project already
made five times.
"""

from __future__ import annotations

__all__ = ["NOTEBOOKS"]


_PORE = [
    ("markdown", """
# 3 · From the pore to a current

Whether a channel conducts is three questions, not one, and they are kept apart
here because a pore can be shut in more than one way:

1. **Is it wide enough?**
2. **Would water stay in it?** A pore can be geometrically open and still block,
   because a hydrophobic neck expels liquid water.
3. **What current would flow?**

Run `python -m piezo1.io.fetch` first — question 2 needs the CHAP hydration
grid as well as the structures.
"""),
    ("code", """
import numpy as np

from piezo1.config import STRUCTURE_DIR
from piezo1.core.structure import Structure
from piezo1.structure.pore import pore_profile
from piezo1.structure.protomers import protomer_blocks
from piezo1.structure.superpose import detect_c3_axis

st = Structure.from_file(STRUCTURE_DIR / "11ZC.cif")
blocks, _ = protomer_blocks(st)
profile = pore_profile(st, detect_c3_axis(blocks), step=1.0)

print(f"bottleneck radius : {profile.bottleneck_radius:.2f} A")
print(f"at z              : {profile.bottleneck_z:.1f} A")
print(f"constrictions     : {len(profile.constrictions())}")
"""),
    ("markdown", """
## Why the probe is on a leash

The radius at each height is the largest sphere that fits — but "largest sphere
that fits" has no interior maximum in an open system. Without a constraint
tethering the probe near the conduction axis, the optimiser walks out of the
protein entirely and reports a radius of about 6,000 Å. That is a true maximum,
and completely useless.

`pore.leash` is a registered parameter, so you can see it, change it, and see
what it does to the answer.
"""),
    ("code", """
from piezo1.parameters import PARAMETERS, reset, set_value

print("leash  :", PARAMETERS.value("pore.leash"), "A")
print("step   :", PARAMETERS.value("pore.step"), "A")

set_value("pore.step", 2.0)
coarse = pore_profile(st, detect_c3_axis(blocks))
print(f"\\nbottleneck at 2.0 A slices: {coarse.bottleneck_radius:.2f} A")
reset()                       # always put the registry back
print("restored:", PARAMETERS.value("pore.step"), "A")
"""),
    ("markdown", """
## Would water stay in it?

Radius alone predicts the conducting state at AUROC 0.59. Radius *combined
with* the hydrophobicity of the lining reaches 0.91 (Rao et al. 2019). The
heuristic scores the pore against an MD-derived free-energy grid; above 0.55
the lining would dewet and the channel is shut whatever its radius.

The two ways of being shut are reported **separately**, because PIEZO1 has
structures that are one and not the other.
"""),
    ("code", """
from piezo1.analysis.hydration import load_grid, predict_wetting

for pdb in ("8YEZ", "11ZC"):
    entry = Structure.from_file(STRUCTURE_DIR / f"{pdb}.cif")
    eb, _ = protomer_blocks(entry)
    prof = pore_profile(entry, detect_c3_axis(eb), step=1.0)
    wet = predict_wetting(entry, prof, load_grid())
    print(f"{pdb}: score {wet.score:5.2f}  "
          f"hydrophobic gate {str(wet.hydrophobic_gate):5s}  "
          f"sterically occluded {str(wet.sterically_occluded):5s}")
    print(f"      {wet.verdict}")
"""),
    ("markdown", """
## What current would flow?

A one-dimensional drift-diffusion calculation over the measured pore, with the
access resistance at each mouth. The potential is solved in the electroneutral
limit, because the Debye length here (5.7–8.1 Å) is larger than the pore radius
(3.3 Å) — the usual Poisson–Nernst–Planck assumption of a well-screened channel
simply does not hold, and the module says so rather than proceeding quietly.
"""),
    ("code", """
from piezo1.physics.permeation import (blocking_mechanisms, default_species,
                                       solve_pnp)

blocks_11zc, _ = protomer_blocks(st)
prof = pore_profile(st, detect_c3_axis(blocks_11zc), step=1.0)
wet = predict_wetting(st, prof, load_grid())

result = solve_pnp(prof, default_species())
print(f"current      : {result.current * 1e12:6.2f} pA")
print(f"conductance  : {result.conductance * 1e12:6.1f} pS   (published 25-30)")
print(f"converged    : {result.converged}")
print(f"blocked by   : {result.blocked_by}")
print(f"mechanisms   : {blocking_mechanisms(wet, prof.radius, default_species())}")
"""),
    ("markdown", """
## The number that does not agree, and why that is reported

40.7 pS against a published 25–30 pS. This project does **not** tune it into
agreement, because two of the inputs have never been measured for PIEZO1: the
in-pore diffusivity and the effective ion radius. Across the plausible ranges
of **both**, the answer spans 16–94 pS, which contains the published value
comfortably. The cell below varies only the first of the two, so it shows a
narrower spread than that.

Agreement reached by choosing values inside that range would be fitting, not
prediction. So the disagreement stands, with the reason attached.
"""),
    ("code", """
from piezo1.analysis.uncertainty import parameter_range

def conductance_at(_scale):
    # parameter_range sets the registry key before each call and restores it
    # afterwards, even if the statistic raises - so the argument is the value
    # already in force rather than something to apply by hand.
    return solve_pnp(prof, default_species()).conductance * 1e12


spread = parameter_range(conductance_at, "permeation.diffusion_scale",
                         [0.4, 0.7, 1.0, 1.5, 2.0],
                         what="single-channel conductance, pS")
print(spread.summary())
print("\\nThis is a PARAMETER range propagated from an unmeasured input.")
print("It is not a confidence interval and the class will not call it one.")
"""),
    ("markdown", """
## Seeing it move

`piezo1.render.flux` turns the computed current into an animation time base. A
channel passes about 10⁷ ions per second, so any watchable stream runs roughly
a millionfold slow — and the number is computed from the solver's own output
rather than chosen to look good, so the display can state it.
"""),
    ("code", """
from piezo1.render.flux import ion_rate, timebase

current_pA = result.current * 1e12
tb = timebase(current_pA)
print(f"{tb.ions_per_second:.3e} ions per second")
print(f"slowdown for a watchable stream: {tb.slowdown:,.0f}x")
print(f"\\n{tb.statement()}")
"""),
]


_VARIANTS = [
    ("markdown", """
# 4 · Variants, and the result that did not work

This notebook is the awkward one, and it is meant to be.

The project was built to predict whether a PIEZO1 variant causes **gain** or
**loss** of function from structure. That is the whole point of the machinery
in the other three notebooks. It does not work — five pre-registered tests,
five different predictor families, five nulls — and, more usefully, the data
that would settle it does not exist and cannot be assembled.

If you take the machinery and skip this, you will repeat a mistake this project
already made five times.
"""),
    ("code", """
from piezo1.core.annotations import load_annotations

ann = load_annotations("human")
by_class = {}
for v in ann.variants:
    by_class.setdefault(v.classification, []).append(v)

print(f"{len(ann.variants)} curated variants, every wild-type residue verified")
for name, group in sorted(by_class.items(), key=lambda kv: -len(kv[1])):
    print(f"  {name:24s} {len(group):3d}")
"""),
    ("markdown", """
## Coverage, reported rather than hidden

Most variants of interest are not resolved in any deposited human structure.
All six human entries model from residue 570 only. A viewer that highlighted
nothing would be indistinguishable from one that had found nothing, so the
count is stated.
"""),
    ("code", """
unresolved = [v for v in ann.variants if not v.modelled_in]
print(f"{len(unresolved)} of {len(ann.variants)} variants are resolved in NO "
      f"human structure")
print("among them:", ", ".join(sorted(v.label for v in unresolved)[:6]), "...")
"""),
    ("markdown", """
## What a prediction is worth

`prediction_record` holds the record as data rather than prose, so the GUI, the
command line, the tour and this notebook cannot drift apart. Each entry is a
test that was **pre-registered** — hypothesis, statistic and decision rule
fixed and committed *before* the comparison was run.
"""),
    ("code", """
from piezo1.analysis.prediction_record import ALL_PREREGISTERED, headline

print(headline())
print()
for e in ALL_PREREGISTERED:
    p = "  n/a " if e.p_value is None else f"{e.p_value:6.3f}"
    print(f"  Round {e.round:2d}  delta {e.cliffs_delta:+.3f}  p {p}  "
          f"n {e.n_gof:2d}/{e.n_lof:2d}  {e.predictor}")
"""),
    ("markdown", """
## Why five nulls are not "nearly something"

One of them, Round 41, returned p = 0.0477 — below the conventional 0.05. It is
still a null, and it was declared one, because the decision rule fixed in
advance had **three** clauses: p below threshold, effect in the predicted
direction, *and* a confidence interval excluding zero. The interval spanned
zero, so the rule was not met.

Deciding that after seeing the number is how a null becomes a finding.
"""),
    ("code", """
from piezo1.analysis.prediction_record import what_it_means

for line in what_it_means():
    print("*", line)
"""),
    ("markdown", """
## The part that makes this more than a list of failures

A null result normally means *get more data*. Here the amount of data that
would be needed was costed, and it exceeds what could ever exist.

**Across positions:** the effect the best predictor produces would need 134
directional variants. The ceiling — every curated variant, every ClinVar entry
with an inferable direction, plus everything the literature harvest could add —
is 59.
"""),
    ("code", """
from piezo1.analysis.feasibility import assess

report = assess(n_simulations=300)
print(report.summary())
"""),
    ("markdown", """
**Within positions:** comparing two variants at the *same* residue removes the
between-position variance, which consumed 99.8% of the first predictor's
signal. It is the obvious way out, and it is closed too: the design needs 8
shared positions even at an implausibly good predictor, and the curated and
ClinVar sets together contain exactly **one**.
"""),
    ("code", """
from piezo1.analysis.data_routes import evidence_summary

for key, value in evidence_summary().items():
    print(f"  {key:34s} {value}")
"""),
    ("markdown", """
## So what is reusable?

Not the predictor. The apparatus that established it could not be validated:

* **pre-registration** with a decision rule fixed in advance, so a marginal p
  cannot be promoted after the fact;
* a **negative control** in every test — in Round 48 the control out-performed
  every mechanistic endpoint, which is what a null looks like from the inside;
* **feasibility costed before another attempt**, which is what turned "we need
  more data" into "the data that could exist is not enough";
* every **checking instrument calibrated** against a known answer before its
  disagreement is believed. Six times in this project, the instrument built to
  check the pipeline was itself the thing at fault.

`docs/METHODS_NOTE.md` writes this up for someone else's project.
"""),
    ("code", """
# The coupling score still exists and still computes. What was removed in
# Round 58 is the claim that its sign means gain or loss of function - the
# reading five pre-registered tests failed to support.
from piezo1.analysis.variant_impact import CouplingScore

print("attributes:", [f for f in CouplingScore.__dataclass_fields__])
print()
print("There is deliberately no `.direction` here, and no alias for it.")
print("A test enforces that, because the name was the misleading part.")
"""),
]


NOTEBOOKS = {
    "03_pore_to_current": {
        "title": "From the pore to a current",
        "cells": _PORE,
    },
    "04_variants_and_the_null": {
        "title": "Variants, and the result that did not work",
        "cells": _VARIANTS,
    },
}
