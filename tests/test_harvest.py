"""Harvesting mutant phenotypes from the open-access corpus.

The round's premise was that published tables are the route to the ~130
phenotyped variants the blind tests need. The harvest yields 2, and neither
carries a direction. These tests hold the funnel, and — more importantly — hold
the gate that makes the count meaningful.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.harvest import (MEASUREMENT_PATTERNS, SUBSTITUTION,
                                     Candidate, harvest)


@pytest.fixture(scope="module")
def report():
    result = harvest()
    if not result.candidates:
        pytest.skip("no open-access papers downloaded — "
                    "run python scripts/build_references.py --download")
    return result


# ------------------------------------------------------------- the gate

def test_the_wild_type_gate_rejects_a_quarter_of_raw_hits(report):
    """23% of regex hits are not PIEZO1 protein substitutions at all.

    cDNA changes are written in the same shape (C7366T), and unrelated text
    matches too. Without the gate those would enter the set as variants, which
    is how a curated resource quietly acquires nonsense.
    """
    total = len(report.candidates)
    passing = len(report.passing())
    assert total > 50
    assert 0.15 < (total - passing) / total < 0.35


def test_the_gate_catches_a_cdna_change_written_like_a_protein_one():
    """The specific false positive, in isolation."""
    from piezo1.analysis.harvest import _numbering_of

    # C7366T is a cDNA change; position 7366 is past the end of either protein.
    assert _numbering_of("C", 7366) == "none"
    # R2456H is real and human.
    assert _numbering_of("R", 2456) in ("human", "both")


def test_most_of_the_literature_is_mouse_numbered(report):
    """The project's standing trap, measured.

    Most functional literature uses mouse and most disease variants use human.
    A harvest that assumed one numbering would silently mis-assign the majority.
    """
    counts = report.by_numbering()
    assert counts.get("mouse", 0) > counts.get("human", 0)
    assert counts.get("mouse", 0) > 25


def test_conversion_goes_through_the_alignment_not_arithmetic(report):
    """A constant offset would be wrong; the map is not linear."""
    from piezo1.core.sequence import mouse_to_human

    mouse = [c for c in report.candidates if c.numbering == "mouse"
             and c.human_label]
    assert len(mouse) > 10
    offsets = {c.position - int(c.human_label[1:-1]) for c in mouse}
    assert len(offsets) > 1, (
        "if every offset were the same, a constant would have worked and this "
        "harvest is not exercising the numbering map")

    for candidate in mouse[:5]:
        assert mouse_to_human(candidate.position) == int(
            candidate.human_label[1:-1])


# ----------------------------------------------------------- the funnel

def test_the_curated_set_already_holds_most_of_what_is_findable(report):
    """31 of 66 gated candidates are already curated — the curation was good.

    Worth pinning because it bounds how much a harvest could ever have added.
    """
    gated = report.passing()
    already = [c for c in gated if c.already_curated]
    assert len(already) > 20
    assert len(already) / len(gated) > 0.35


def test_almost_nothing_carries_an_extractable_measurement(report):
    """The bottleneck, and it is not the gate.

    Of the fresh candidates, nearly all appear only in prose. The numbers this
    project needs are in sentences and in non-open-access supplements, not in
    the machine-readable tables the round assumed.
    """
    fresh = [c for c in report.passing()
             if c.human_label and not c.already_curated]
    with_measurement = [c for c in fresh if c.measurements]
    assert len(fresh) > 20
    assert len(with_measurement) < 5, (
        f"{len(with_measurement)} of {len(fresh)} carry a number; if this has "
        f"risen substantially the harvest is worth redoing")
    assert len(report.new_usable()) == len(with_measurement)


def test_no_direction_is_ever_assigned(report):
    """The line this module must not cross.

    Reading "slowed inactivation" out of prose and calling it gain-of-function
    would put unreviewed labels into the set the blind tests depend on. The
    Candidate type has no direction field at all, which is the enforcement.
    """
    assert not hasattr(Candidate, "classification")
    assert not hasattr(Candidate, "direction")
    for candidate in report.candidates[:20]:
        assert not hasattr(candidate, "direction")
    assert "NOT assigned" in report.summary()


# ------------------------------------------------------------ the parsing

def test_measurement_patterns_find_what_they_claim():
    text = ("The inactivation tau was 22.2 ms for the mutant, with a unitary "
            "conductance of 27 pS and a current density of -45.3 pA/pF.")
    assert MEASUREMENT_PATTERNS["inactivation_tau"].search(text).group(1) == "22.2"
    assert MEASUREMENT_PATTERNS["conductance"].search(text).group(1) == "27"
    assert MEASUREMENT_PATTERNS["current_density"].search(text).group(1) == "-45.3"


def test_the_substitution_pattern_is_loose_on_purpose():
    """It must admit cDNA changes, because tightening it would drop real ones."""
    assert SUBSTITUTION.search("R2456H")
    assert SUBSTITUTION.search("C7366T"), (
        "the pattern is deliberately loose; the wild-type gate is what "
        "rejects this, not the regex")
    assert not SUBSTITUTION.search("R24H"), "two digits is too few to be a variant"


def test_context_is_captured_for_human_review(report):
    """A candidate without its sentence cannot be curated by anyone."""
    for candidate in report.passing()[:10]:
        assert candidate.context
        assert candidate.source
