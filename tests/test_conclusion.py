"""Every number in the conclusion must come from the code, not from prose.

`docs/CONCLUSION.md` is the one page a reader is meant to trust without reading
the working. This project has twice shipped prose that went stale — the tour
still said "tested twice" after five tests, and a review counted 40 usable
positions where there was one — so the document is guarded rather than
proof-read.

The guard extracts every number from the document and requires each to be
supported by `CLAIMS`, `ALL_PREREGISTERED`, `HEADLINE`, or a small allowlist of
structural constants each of which states why it is exempt.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

from piezo1.analysis.claims import CLAIMS
from piezo1.analysis.prediction_record import ALL_PREREGISTERED
from piezo1.analysis.published_interval import HEADLINE

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "CONCLUSION.md"
README = ROOT / "README.md"

#: Numbers that are structural facts or counts rather than computed results.
#: Each must say where it comes from; an unexplained entry is how a guard like
#: this decays into a rubber stamp.
ALLOWED = {
    68: "curated variants in resources/variants.json",
    232: "ClinVar variants in resources/variants_clinvar.json",
    34: "variants Round 36 could model, recorded in VALIDATION_ROUND36.md",
    2456: "residue number, the one usable within-position site",
    870: "residue number (M870V)",
    1358: "residue number (R1358C)",
    2020: "residue number (A2020V)",
    7: "a roadmap round number", 22: "a roadmap round number",
    26: "a roadmap round number", 36: "a roadmap round number",
    41: "a roadmap round number", 47: "a roadmap round number",
    48: "a roadmap round number", 54: "a roadmap round number",
    50: "a roadmap round number",
    0.4: "derived: |2.71 - 2.7| / 2.7, the T50 agreement stated in SCIENCE.md",
    0.01: "derived: |13.998 - 14.0| / 14.0, the decay-length round trip",
    5: "the number of pre-registered tests",
    6: "the number of times a checking instrument was itself wrong",
    3: "the number of variants one label away from a usable pair",
    1: "the number of usable within-position positions (R2456)",
    2: "small counts in prose (two routes, two of the three VUS)",
    0.05: "the conventional significance threshold",
    0.2: "the Round 26 pre-registered within-position criterion",
    0.55: "Rao 2019 wetting cutoff, a published threshold",
    10.2: "published dome radius (Haselwandter & MacKinnon 2018)",
    2.7: "published T50 (Lewis & Grandl 2015)",
    0.1: "the uncertainty on the published T50",
    14.0: "the decay length imposed on the solver by construction",
    20: "the bending modulus kappa in k_BT, a published range endpoint",
    463: "the linear route's footprint area, recorded in SCIENCE.md",
    63: "the dome contact slope in degrees, recorded in SCIENCE.md",
    99.8: "Round 7's between-position variance share, in VALIDATION.md",
    25: "published conductance range low (Coste 2010)",
    30: "published conductance range high (Shi 2020)",
    41: "the model's conductance, pinned in test_permeation.py",
    0.9: "the dome bootstrap half-width, from published_interval",
    2019: "a citation year (Rao)", 2018: "a citation year (Haselwandter)",
    2015: "a citation year (Lewis & Grandl)", 2010: "a citation year (Coste)",
}

_NUMBER = re.compile(r"(?<![\w/])[-+−]?\d+(?:[.,]\d+)?(?![\w])")


def _numbers_in(text: str) -> set:
    out = set()
    for match in _NUMBER.finditer(text):
        raw = match.group().replace("−", "-").replace(",", "")
        try:
            value = float(raw)
        except ValueError:
            continue
        out.add(value)
    return out


def _supported() -> set:
    values = {abs(c.expected) for c in CLAIMS} | {c.expected for c in CLAIMS}
    for entry in ALL_PREREGISTERED:
        values |= {entry.cliffs_delta, abs(entry.cliffs_delta),
                   float(entry.n_gof), float(entry.n_lof), float(entry.round)}
        if entry.p_value is not None:
            values.add(entry.p_value)
    for headline in HEADLINE:
        values |= {headline.estimate, headline.low, headline.high}
        values.add(round(headline.overconfident_by))
    return values


@pytest.fixture(scope="module")
def text():
    assert DOC.exists(), "docs/CONCLUSION.md is missing"
    return DOC.read_text()


@pytest.fixture(scope="module")
def flowed(text):
    """The prose with wrapping removed, so a phrase check is not a line check."""
    return " ".join(text.split())


# ------------------------------------------------- the guard

def test_every_number_is_traceable(text):
    """The load-bearing test: no number may be here that the code cannot make."""
    supported = _supported()
    unexplained = []
    for value in sorted(_numbers_in(text)):
        if value in ALLOWED or abs(value) in ALLOWED:
            continue
        if any(abs(value - s) <= max(abs(s) * 0.005, 1e-9) for s in supported):
            continue
        unexplained.append(value)
    assert not unexplained, (
        "these numbers are in CONCLUSION.md but come from neither the claims "
        f"registry, the validation record, nor ALLOWED: {unexplained}")


def test_the_guard_would_catch_an_invented_number(text):
    """Calibration: a guard that cannot fail asserts nothing."""
    supported = _supported()
    invented = 4271.9
    assert invented not in ALLOWED
    assert not any(abs(invented - s) <= max(abs(s) * 0.005, 1e-9)
                   for s in supported)


def test_every_allowlist_entry_states_a_reason():
    for value, reason in ALLOWED.items():
        assert len(reason) > 12, f"{value} is allowed without a real reason"


# ------------------------------------------- what it must not soften

def test_it_states_all_five_nulls(text):
    for entry in ALL_PREREGISTERED:
        assert f"{entry.cliffs_delta:+.3f}".replace("+", "").replace("-", "") in \
            text.replace("−", "-").replace("+", ""), (
            f"Round {entry.round}'s effect size is missing")
    assert "five nulls" in text.lower()


def test_it_states_the_unprovability_result_not_just_the_nulls(flowed):
    """Rounds 47 and 54 are what make this more than a list of failures."""
    assert "cannot be settled with data that could exist" in flowed
    assert "134" in flowed and "59" in flowed
    assert "should not be run" in flowed


def test_it_does_not_soften_the_null(text):
    lowered = text.lower()
    for banned in ("trend toward", "approaching significance",
                   "marginally significant", "promising", "encouraging"):
        assert banned not in lowered, f"softening language: {banned!r}"


def test_it_states_what_is_not_claimed(flowed):
    assert "not a clinical tool" in flowed.lower()
    assert "reported, not tuned" in flowed.lower()
    assert "14.99" in flowed, "the dome's model spread must appear"


def test_it_links_the_figures_the_tour_uses(text):
    for figure in ("record_nulls.png", "record_data_limit.png"):
        assert figure in text


def test_it_is_reachable_from_the_readme_and_the_help():
    readme = README.read_text()
    assert "CONCLUSION.md" in readme, "the conclusion must be linked from README"
    from piezo1.ui.help_content import DOC_LINKS

    assert any("CONCLUSION" in str(v) for v in
               (DOC_LINKS.values() if isinstance(DOC_LINKS, dict) else DOC_LINKS))


# --------------------------- the README must end where the science does

def _readme_summary() -> str:
    """The closing section, which Round 59 required to match this page."""
    text = README.read_text()
    start = text.index("## What this established, and what it did not")
    return text[start:text.index("## Data sources", start)]


def test_the_readme_ends_on_the_record_not_on_the_machinery():
    summary = _readme_summary()
    assert "five nulls" in summary.lower()
    assert "does not work" in summary
    assert "CONCLUSION.md" in summary


def test_the_readme_summary_states_the_unprovability_result():
    """A list of nulls is not the finding; the closed routes are."""
    flowed = " ".join(_readme_summary().split())
    assert "134" in flowed and "59" in flowed, "the across-position limit"
    assert "one" in flowed and "usable site" in flowed, "the within-position limit"
    assert "should not be run" in flowed


def test_every_number_in_the_readme_summary_is_traceable():
    """The same guard as the conclusion page, on the section that mirrors it."""
    supported = _supported()
    unexplained = []
    for value in sorted(_numbers_in(_readme_summary())):
        if value in ALLOWED or abs(value) in ALLOWED:
            continue
        if any(abs(value - s) <= max(abs(s) * 0.005, 1e-9) for s in supported):
            continue
        unexplained.append(value)
    assert not unexplained, (
        f"untraceable numbers in the README summary: {unexplained}")


def test_the_readme_summary_does_not_soften_the_null():
    lowered = _readme_summary().lower()
    for banned in ("trend toward", "approaching significance", "promising",
                   "encouraging", "preliminary evidence"):
        assert banned not in lowered, f"softening language: {banned!r}"


# ------------------- the three places must not drift from one another

def test_the_tour_the_readme_and_the_conclusion_agree_on_the_count():
    """Three surfaces state the record; a sixth test must update all three.

    They went out of step before: the tour said "tested twice" after five tests
    had run, and it survived because each surface was written by hand at a
    different time. The tour now reads the record, the conclusion is guarded,
    and this asserts the README cannot be the one left behind.
    """
    from piezo1.tour import step_by_key

    n = len(ALL_PREREGISTERED)
    words = {2: "two", 3: "three", 4: "four", 5: "five", 6: "six"}
    word = words[n]

    tour_text = step_by_key("record").report({})
    assert f"{n} pre-registered tests" in tour_text

    for name, body in (("CONCLUSION.md", DOC.read_text()),
                       ("README.md", _readme_summary())):
        lowered = " ".join(body.split()).lower()
        assert f"{word} pre-registered" in lowered or f"{word} nulls" in lowered, (
            f"{name} does not state {word} pre-registered tests")


def test_no_surface_still_says_the_project_is_pursuing_the_claim():
    """The framing Round 59 was written to remove.

    Not a search for the word "predict" — the project legitimately predicts a
    gating motion and a wetting verdict. This looks for the claim being
    presented as open.
    """
    from piezo1.tour import TOUR

    surfaces = {"README.md": README.read_text(),
                "CONCLUSION.md": DOC.read_text()}
    surfaces.update({f"tour:{s.key}": s.body for s in TOUR})
    for name, body in surfaces.items():
        lowered = " ".join(body.split()).lower()
        for phrase in ("aims to predict", "will predict whether",
                       "sets out to predict gain", "hopes to"):
            assert phrase not in lowered, f"{name} still presents the claim as open"
