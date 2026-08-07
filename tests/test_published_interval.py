"""What gets published beside each headline number.

Round 38 measured that the dome radius's model spread is six times its
bootstrap interval. Round 52 acts on it: the interval quoted must be the widest
term, named for its kind. These tests guard the rule and the four numbers it
was applied to.

The point estimates do not move, so `verify_claims` is untouched — a claim
tolerance detects code drift and is a different question from what the science
section is entitled to say.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from piezo1.analysis.published_interval import (HEADLINE, KINDS,
                                                PublishedInterval, Term,
                                                publish)

DOC = Path(__file__).resolve().parents[1] / "docs" / "SCIENCE.md"


# ------------------------------------------------------------ the rule

def test_only_a_bootstrap_may_be_called_a_confidence_interval():
    """A cutoff has no sampling distribution; nor does a choice of shape."""
    assert KINDS["bootstrap"] is True
    for kind in ("sensitivity", "parameter", "model"):
        assert KINDS[kind] is False
        assert not Term(kind, 1.0, 2.0, "x").is_confidence_interval


def test_a_non_bootstrap_term_says_it_is_not_a_confidence_interval():
    for kind in ("sensitivity", "parameter", "model"):
        assert "NOT a CI" in Term(kind, 1.0, 2.0, "x").describe()
    assert "95% CI" in Term("bootstrap", 1.0, 2.0, "x").describe()


def test_a_model_term_declares_itself_a_lower_bound():
    """Two models agreeing does not bound the error from above."""
    assert "LOWER BOUND" in Term("model", 1.0, 2.0, "x").describe()


def test_the_widest_term_is_the_one_published():
    interval = PublishedInterval(
        quantity="q", estimate=1.0, unit="u",
        terms=[Term("bootstrap", 0.9, 1.1, "narrow"),
               Term("model", 0.5, 2.0, "wide")])
    assert interval.dominant.kind == "model"
    assert (interval.low, interval.high) == (0.5, 2.0)


def test_an_unknown_kind_is_rejected():
    with pytest.raises(ValueError, match="unknown kind"):
        Term("vibes", 1.0, 2.0, "x")


def test_an_inverted_interval_is_rejected():
    with pytest.raises(ValueError, match="below low"):
        Term("bootstrap", 2.0, 1.0, "x")


# ------------------------------------------------- the four quantities

def test_every_headline_number_has_more_than_one_term():
    """A single term cannot establish which kind dominates."""
    for entry in HEADLINE:
        assert entry.terms, f"{entry.quantity} has no terms"
        if len(entry.terms) == 1:
            # Legitimate only when no other kind applies; say so explicitly.
            assert entry.terms[0].note, (
                f"{entry.quantity} has one term and no note explaining why")


def test_the_dome_radius_publishes_the_model_spread_not_the_bootstrap():
    """The finding Round 38 made and Round 52 acts on."""
    dome = publish("dome")
    assert dome.dominant.kind == "model"
    assert dome.overconfident_by > 10, (
        "quoting the narrowest term would understate the spread by a lot")
    bootstrap = next(t for t in dome.terms if t.kind == "bootstrap")
    assert dome.dominant.width > 3 * bootstrap.width


def test_the_dome_records_that_the_two_fits_were_not_like_for_like():
    """Found while writing this round, and kept rather than quietly fixed.

    The model comparison is anchored on the untrimmed sphere fit (9.45 nm)
    while the published number is trimmed (9.72 nm). The 0.27 nm gap does not
    change the conclusion against a 5.54 nm model spread, but the mismatch was
    real and unnoticed.
    """
    dome = publish("dome")
    trim = next(t for t in dome.terms if "trim" in t.source)
    assert "UNTRIMMED" in trim.note
    assert trim.width == pytest.approx(0.30, abs=0.02)
    assert trim.width < dome.dominant.width / 10, (
        "the trim choice must remain small against the shape choice")


def test_the_gating_overlap_publishes_the_cutoff_range():
    overlap = publish("gating overlap")
    assert overlap.dominant.kind == "sensitivity"
    assert "cutoff" in overlap.dominant.source
    # The spring-model spread is real but smaller; it must not be what is quoted.
    springs = next(t for t in overlap.terms if t.kind == "model")
    assert springs.width < overlap.dominant.width


def test_t50_is_limited_by_its_input_rates_not_by_the_solver():
    t50 = publish("T50")
    assert t50.dominant.kind == "parameter"
    solver = next(t for t in t50.terms if t.kind == "sensitivity")
    assert solver.width < 0.05, "the two solvers agree closely"
    assert t50.dominant.width > 10 * solver.width


def test_the_published_t50_measurement_lies_inside_the_published_range():
    """The agreement with Lewis 2015 must survive the input uncertainty.

    If 2.7 fell outside the range implied by the rate constants, the reported
    agreement would be a coincidence of the chosen values.
    """
    t50 = publish("T50")
    assert t50.low <= 2.7 <= t50.high
    assert t50.low <= 2.711 <= t50.high


def test_the_footprint_is_limited_by_kappa():
    footprint = publish("footprint")
    assert footprint.dominant.kind == "parameter"
    assert "kappa" in footprint.dominant.source


# --------------------------------------------------- the documentation

def test_the_science_section_publishes_the_dominant_term():
    """The rule is worthless if the document still quotes the narrow one."""
    text = DOC.read_text()
    assert "9.45" in text and "14.99" in text, (
        "SCIENCE.md must state the dome's model range")
    assert "0.554" in text and "0.723" in text, (
        "SCIENCE.md must state the gating overlap's cutoff range")
    assert "2.584" in text or "2.58" in text, (
        "SCIENCE.md must state T50's input-propagated range")


def test_the_science_section_does_not_call_a_model_spread_a_confidence_interval():
    text = DOC.read_text()
    for line in text.splitlines():
        low = line.lower()
        if "confidence interval" in low and "9.45" in line:
            pytest.fail(f"model spread described as a confidence interval: {line}")


def test_every_headline_statement_names_its_kind():
    for entry in HEADLINE:
        statement = entry.statement()
        assert any(word in statement
                   for word in ("95% CI", "NOT a CI")), statement
