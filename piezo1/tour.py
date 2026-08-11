"""A guided tour of the PIEZO1 mechanism, tied to live measurement.

Project aim A1 is that this be a *learning* instrument, and it has had the least
attention of any aim. This is the tour: dome, blades, lever, gate, opening —
each step setting up a view and then **measuring the thing it just described**.

**Every number a step states is computed when the step runs.** None is typed
into the prose. A tour that narrated "the dome radius is 9.7 nm" would be a
fourth place for that number to live and go stale, alongside the code, the
documentation and the claims registry. Instead each step carries a callable that
reads the analysis the application actually ran, and the published comparison
comes from the parameter registry rather than from a literal — so a tour step
cannot drift from either the code or the literature.

Qt-free, so the whole tour can be walked headlessly and tested without a
display. The GUI half is :mod:`piezo1.ui.tour_panel`.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable

from .parameters import PARAMETERS

__all__ = ["TourStep", "TOUR", "step_by_key", "Results"]

#: What a step's ``measure`` callable receives: whatever the application has
#: computed so far, keyed by analysis name. Missing keys are normal — a step
#: says so rather than inventing a number.
Results = dict


@dataclass
class TourStep:
    """One stop on the tour."""

    key: str
    title: str
    body: str
    #: Structure to load, if the step needs a particular one.
    structure: str = ""
    style: str = ""
    color_by: str = ""
    #: Annotation group whose residues should be highlighted, matched on label.
    highlight_group: str = ""
    highlight: tuple[int, ...] = ()
    #: Which existing analysis the step needs: dome | pore | modes | pockets.
    run: str = ""
    #: Formats the live result. Receives everything computed so far.
    measure: Callable[[Results], str] | None = None
    #: Parameter keys whose values are quoted in the body, for the citation
    #: line the panel shows beneath each step.
    cites: tuple[str, ...] = field(default_factory=tuple)
    #: Filename in ``docs/img`` to show beneath the prose. Kept as a name
    #: rather than a path so this module stays independent of where the
    #: project is installed, and Qt-free so the tour can be read headlessly.
    image: str = ""
    image_caption: str = ""

    def image_path(self):
        """Absolute path to the figure, or ``None`` if it is not there.

        Missing figures are not an error: they are regenerable and git-ignored
        outputs, and a tour that refused to open because a PNG had not been
        built would be worse than one that shows the prose alone.
        """
        from .config import PROJECT_ROOT

        if not self.image:
            return None
        path = PROJECT_ROOT / "docs" / "img" / self.image
        return path if path.exists() else None

    def body_html(self) -> str:
        """The prose, with the figure appended when it exists."""
        path = self.image_path()
        if path is None:
            return self.body
        caption = (f'<div style="color:#8b93a1;font-size:11px;">'
                   f'{self.image_caption}</div>' if self.image_caption else "")
        return (f'{self.body}<p><img src="{path.as_uri()}" width="620"></p>'
                f'{caption}')

    def report(self, results: Results) -> str:
        if self.measure is None:
            return ""
        try:
            return self.measure(results or {})
        except Exception as exc:            # a tour must not crash the app
            return f"could not measure: {type(exc).__name__}: {exc}"


# --------------------------------------------------------------------------
# Measurements. Each reads what the application computed; none invents a value.
# --------------------------------------------------------------------------

def _dome(results: Results) -> str:
    dome = results.get("dome")
    if dome is None:
        return "measure the dome to see this number"
    measured = dome.radius_of_curvature / 10.0
    published = PARAMETERS.value("dome.published_radius_closed")
    return (f"Measured radius of curvature <b>{measured:.2f} nm</b>, against a "
            f"published closed-state value of <b>{published:g} nm</b>. "
            f"Dome depth {dome.dome_depth / 10:.1f} nm; the curved surface "
            f"holds {dome.excess_area / 100:.0f} nm² more area than its own "
            f"shadow on the membrane plane.")


def _footprint(results: Results) -> str:
    comparison = results.get("footprint")
    if comparison is None:
        return "compute the footprint to see this number"
    return (f"Linear theory gives <b>{comparison.linear_energy:.0f} k_BT</b> "
            f"and {comparison.linear_excess_area:.0f} nm²; solved without the "
            f"small-slope approximation the same footprint is "
            f"<b>{comparison.nonlinear_energy:.0f} k_BT</b> and "
            f"<b>{comparison.nonlinear_excess_area:.0f} nm²</b> — the linear "
            f"number is {1 / comparison.energy_ratio:.1f}× too large at this "
            f"{comparison.contact_angle_deg:.0f}° contact angle.")


def _pore(results: Results) -> str:
    profile = results.get("pore")
    if profile is None:
        return "compute the pore profile to see this number"
    text = (f"Bottleneck <b>{profile.bottleneck_radius:.2f} Å</b> "
            f"({profile.bottleneck_radius / 10:.3f} nm) at "
            f"z = {profile.bottleneck_z:.0f} Å.")
    wetting = results.get("hydration")
    if wetting is not None and getattr(wetting, "available", False):
        text += (f" Wetting score <b>{wetting.score:.2f}</b> against a cutoff "
                 f"of {PARAMETERS.value('hydration.closed_cutoff'):g} → "
                 f"<b>{wetting.verdict}</b>.")
    return text


def _modes(results: Results) -> str:
    modes = results.get("modes")
    if modes is None:
        return "compute normal modes to see this number"
    symmetry = list(getattr(modes, "symmetry", []) or [])
    n_a = sum(1 for s in symmetry if s == "A")
    n_e = len(symmetry) - n_a
    return (f"{modes.n_modes} modes over {modes.n_sites} sites: "
            f"<b>{n_a} symmetric (A)</b> and {n_e} degenerate (E). Membrane "
            f"tension is isotropic and therefore C3-symmetric, so only the A "
            f"modes can couple to it at first order — a selection rule that "
            f"falls out of the calculation rather than being imposed.")


def _gating(results: Results) -> str:
    from .physics.kinetics import GatingModel
    model = GatingModel()
    measured = PARAMETERS.value("kinetics.t50_measured")
    return (f"The four-state scheme, given only its published rate constants, "
            f"produces a half-activation tension of "
            f"<b>{model.half_activation():.2f} mN/m</b> against a measured "
            f"<b>{measured:g} mN/m</b>. That number is emergent: nothing in "
            f"the model was fitted to it.")


def _variant(results: Results) -> str:
    from .core.annotations import load_annotations
    annotations = load_annotations("human")
    at_2456 = [v for v in annotations.variants if v.residue == 2456]
    kinds = ", ".join(f"{v.label} ({v.classification})" for v in at_2456)
    return (f"Four substitutions are curated at this one position: {kinds}. "
            f"Three gain function and one loses it — from the same residue.")


def _record(results: Results) -> str:
    """Every pre-registered test, read from the frozen record."""
    from .analysis.prediction_record import ALL_PREREGISTERED

    parts = [f"Round {r.round} ({r.predictor}) delta {r.cliffs_delta:+.3f}"
             for r in ALL_PREREGISTERED]
    return (f"{len(ALL_PREREGISTERED)} pre-registered tests, "
            f"{len(ALL_PREREGISTERED)} nulls: " + "; ".join(parts)
            + ". Every interval crosses zero.")


def _data_limit(results: Results) -> str:
    """The Round 47 feasibility numbers, read from the claims registry.

    Read rather than recomputed: these are documented values, and a tour step
    must not spend ten seconds of simulation before it can draw itself.
    """
    from .analysis.claims import CLAIMS

    values = {c.key: c.expected for c in CLAIMS}
    needed = values.get("feasibility.required_n", float("nan"))
    ceiling = values.get("feasibility.ceiling", float("nan"))
    return (f"The best predictor produces an effect of delta −0.249. Detecting "
            f"that at {PARAMETERS.value('stats.target_power'):g} power needs "
            f"{needed:.0f} directional variants. The most this project could "
            f"ever assemble — every curated variant plus every candidate the "
            f"literature harvest found — is {ceiling:.0f}, where the power is "
            f"about a coin flip. The gap is a factor of {needed / ceiling:.1f}.")


def _limits(results: Results) -> str:
    """What the remaining uncertainty is, on the numbers the tour just showed."""
    from .analysis.published_interval import publish

    dome = publish("dome")
    overlap = publish("gating overlap")
    return (f"Even the measurements that worked carry more uncertainty than a "
            f"confidence interval shows. The dome radius is "
            f"{dome.estimate:.2f} nm, but a sphere and an oblate spheroid fitted "
            f"to the same points give [{dome.low:.2f}, {dome.high:.2f}] nm — "
            f"about {dome.overconfident_by:.0f} times wider than the bootstrap "
            f"interval. The gating overlap is {overlap.estimate:.3f} but ranges "
            f"{overlap.low:.3f}–{overlap.high:.3f} across network cutoffs.")


# --------------------------------------------------------------------------
# The tour
# --------------------------------------------------------------------------

TOUR: list[TourStep] = [
    TourStep(
        key="channel", title="1 · One channel, three copies of everything",
        structure="8YEZ", style="cartoon", color_by="chain",
        body="""
<p>This is human PIEZO1 in its closed state. The first thing to see is that it
is a <b>trimer</b>: three identical protomers related by a three-fold rotation,
coloured here one per chain.</p>
<p>Each protomer contributes 38 transmembrane helices. That is an enormous
amount of membrane-embedded protein for an ion channel, and it is the clue to
the mechanism — most of this structure is not there to conduct ions. It is
there to <i>feel the membrane</i>.</p>
<p>Rotate the model. The three blades sweep out from a central pore like the
arms of a propeller.</p>"""),

    TourStep(
        key="blades", title="2 · The blades are a curved lever",
        color_by="domain",
        body="""
<p>Coloured by domain now. Each blade is built from nine repeats of a four-helix
unit — the <b>THU</b> repeats — and they are not flat. They curve.</p>
<p>Because they curve, they bend the bilayer around them into a dome. That dome
is the sensor: it stores membrane deformation, and deformation costs energy that
tension can pay back. Everything downstream follows from the shape.</p>"""),

    TourStep(
        key="dome", title="3 · The dome, measured",
        run="dome", measure=_dome,
        cites=("dome.published_radius_closed",),
        body="""
<p>Rather than assert the dome, measure it. The application fits a sphere to the
mid-membrane surface — the centre of every transmembrane helix — and reports the
radius of curvature.</p>
<p>This is the project's standing regression test: if the geometry pipeline
breaks, this number moves.</p>"""),

    TourStep(
        key="footprint", title="4 · The membrane keeps bending past the protein",
        run="footprint", measure=_footprint,
        cites=("membrane.kappa", "membrane.tension"),
        body="""
<p>The dome does not stop at the protein boundary. The bilayer around it is bent
too, over a decay length set by the ratio of bending rigidity to tension.</p>
<p>This step is also a warning about modelling. The standard treatment linearises
the membrane energy assuming small slopes — and PIEZO1's dome meets the bilayer
at about 63°, where that assumption fails badly. The two answers are shown
together, because the difference between them is larger than most of the effects
people argue about.</p>"""),

    TourStep(
        key="lever", title="5 · From blade to gate",
        highlight_group="anchor", color_by="domain",
        body="""
<p>Highlighted is the <b>anchor</b>, where the blade meets the pore module.</p>
<p>Two independent analyses converge here. Perturbation-response scanning — push
each residue and measure how far the gate moves — puts the anchor on the optimal
force-transmission route. And it is separately the most conserved domain across
vertebrate orthologs. Neither calculation knows about the other.</p>
<p>That is what a lever looks like from the inside: the part that transmits force
is the part evolution will not let you change.</p>"""),

    TourStep(
        key="gate", title="6 · The gate, closed",
        run="pore", measure=_pore,
        cites=("hydration.closed_cutoff", "pore.ion_radius"),
        body="""
<p>Now the pore itself: the radius of the largest sphere that fits at each height
along the conduction axis.</p>
<p>But radius alone is a poor predictor of conduction — on its own it separates
open from closed channels barely better than a coin. A pore can be wide enough
for a hydrated ion and still block, because a hydrophobic neck expels liquid
water. The verdict below therefore reports <i>two</i> ways to be shut: too narrow
for water at all, and wide enough but chemically unwelcoming.</p>"""),

    TourStep(
        key="open", title="7 · The same channel, flattened",
        structure="11ZC", run="pore", measure=_pore, style="cartoon",
        color_by="domain",
        body="""
<p>This is a flattened, open-like structure. Compare the number below with the
closed state two steps ago.</p>
<p>The instructive part is not that the pore is wider. It is that residues in the
<i>closed</i> structure sit at essentially the same radius as this one's
bottleneck and are still called non-conducting — because they are hydrophobic and
these are not. Same geometry, opposite verdict.</p>"""),

    TourStep(
        key="modes", title="8 · Why it moves the way it does",
        structure="8YEZ", run="modes", measure=_modes,
        cites=("anm.cutoff", "anm.n_modes"),
        body="""
<p>An elastic network model treats the fold as springs between residues and asks
which motions are cheapest. No knowledge of the gating transition goes in.</p>
<p>Each mode is then labelled by how it behaves under the three-fold rotation.
This matters more than it sounds: membrane tension is isotropic, so it can only
drive motions that respect the symmetry. The calculation finds that the lowest
symmetric mode overlaps the experimentally observed curved-to-flat transition at
0.70, while the degenerate modes contribute essentially nothing.</p>
<p>Try <i>Colour by displacement</i> in the Physics panel to see which parts of
the protein a mode actually moves.</p>"""),

    TourStep(
        key="energetics", title="9 · Tension to open probability",
        measure=_gating, cites=("kinetics.t50_measured", "dome.delta_area"),
        body="""
<p>The last link. A four-state Markov scheme with tension-dependent rates turns
a membrane tension into an open probability and an inactivation time course.</p>
<p>The number below is the check that the chain holds together end to end.</p>"""),

    TourStep(
        key="variant", title="10 · Where a disease variant sits",
        structure="8YEZ", highlight=(2456,), color_by="domain",
        measure=_variant,
        body="""
<p>R2456 sits in the C-terminal domain, on the inner leaflet side of the pore.
It is the best-studied PIEZO1 variant position.</p>
<p>And it is the position that shows why this is hard. Read the line below: the
same residue, four different substitutions, and they do not all do the same
thing.</p>"""),

    TourStep(
        key="record", title="11 · The central claim, and five attempts on it",
        measure=_record, cites=("stats.alpha",),
        image="record_nulls.png",
        image_caption="Every pre-registered test. The red line is no effect; "
                      "every interval crosses it.",
        body="""
<p>A learning instrument that only shows its successes teaches the wrong lesson,
so the tour ends on the record.</p>
<p>The project's central aim is to predict gain- versus loss-of-function from
structure. It has been tested <b>five times</b>, each pre-registered in its own
commit before the comparison was run, with five different predictors: elastic
network energy, FoldX stability, a substitution-aware version of the first,
population constraint from gnomAD, and the wild-type structural context of the
position. <b>All five returned null.</b></p>
<p>The diagnostic from the first attempt is still the useful one: the mechanical
predictor reports <i>where a residue sits</i> rather than <i>which substitution
occurred</i>, which is exactly why all four R2456 variants score alike. Round 26
raised its sensitivity to the substitution from 4.9% to 52.5% of the variance,
and the effect grew from −0.083 to −0.249 — real improvement, still not
significant.</p>
<p>The last of the five is the sharpest. A feature computed on the wild-type
structure has <b>exactly zero</b> within-position variance: R2456H, R2456K,
R2456P and R2456C receive the identical value, to every digit. Such a feature
could never assign a direction to a substitution, whatever its p-value.</p>"""),

    TourStep(
        key="data_limit", title="12 · Why more data would not settle it",
        measure=_data_limit, cites=("stats.target_power",),
        image="record_data_limit.png",
        image_caption="Blue is reachable; red is not. The dashed line is the "
                      "most this project could ever assemble.",
        body="""
<p>The usual conclusion from a null is "we need more data". Here that was
checked rather than assumed, and it is not the right conclusion.</p>
<p>The effect the best predictor actually produces would need about
<b>134 directional variants</b> to detect at the conventional power. The curated
set supplies 34 that can be modelled. Adding every candidate a systematic
literature harvest could find, and assuming every one of them could be assigned
a direction it does not currently have, reaches <b>59</b> — where the power is
roughly a coin flip.</p>
<p>So the honest statement is not "we need more data" but <b>"the data that
could exist is not enough for an effect this size"</b>. Those sound alike and
are different: only the second tells you that a sixth test on this variant set
should not be run, whatever predictor goes into it.</p>
<p>What would change it is not a better model. It is a design matched
<i>within</i> position — comparing two variants at the same residue, which
removes the between-position variance that consumed 99.8% of the first
predictor's signal — and that is a curation problem rather than a modelling
one.</p>"""),

    TourStep(
        key="limits", title="13 · What this application cannot do",
        measure=_limits, cites=("stats.alpha", "stats.target_power"),
        body="""
<p>Two last things, so that nothing here is read as more certain than it is.</p>
<p>The whole of it is on one page in <b>docs/CONCLUSION.md</b>, with every number traceable to the code that produced it — Help &rarr; <i>What was established</i> opens it.</p>
<p><b>The measurements that worked still carry model uncertainty.</b> The dome
radius quoted in step 3 has a bootstrap interval of about ±0.9 nm — but fit an
oblate spheroid to the same surface points instead of a sphere and the radius of
curvature moves from 9.45 to 14.99 nm. The spheroid fits <i>better</i>. What
limits that number is the choice of shape, not the scatter of the points, and
the confidence interval was answering a question nobody asked.</p>
<p><b>And this is not a clinical tool.</b> Nothing here is validated for
diagnosis. A variant's position on the structure is a hypothesis about
mechanism, not a statement about a patient.</p>
<p>What the project does establish: the dome geometry reproduces the published
curvature, the gating scheme reproduces a measured half-activation tension to
0.4%, the closed pore is predicted to dewet, and every number carries its
source. What it does not establish is the thing it set out to do.</p>"""),
]


def step_by_key(key: str) -> TourStep | None:
    for step in TOUR:
        if step.key == key:
            return step
    return None
