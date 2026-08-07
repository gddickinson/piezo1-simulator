"""Reading all 35 harvested candidates by hand, and what it was worth.

Round 45 found 35 substitutions the curated set does not have. Round 57 asked
the question that decides whether another test is possible: how many carry a
direction a human can recover from the sentence alone?

Five. All five say "non-functional", all five are alanine-scanning mutants, and
none sits at a position that has any other variant — so they add nothing to
either the across-position or the within-position design.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.harvest import harvest
from piezo1.analysis.harvest_curation import (CATEGORIES, CURATION, Verdict,
                                              by_category, directional,
                                              summary)


# ------------------------------------------------ the curation is complete

def test_every_harvested_candidate_was_read():
    """A hand-curation that skipped rows would report a yield it did not earn."""
    fresh = {c.human_label for c in harvest().passing()
             if c.human_label and not c.already_curated}
    curated = {v.label for v in CURATION}
    assert not fresh - curated, f"never read: {sorted(fresh - curated)}"
    assert not curated - fresh, f"curated but not a candidate: {sorted(curated - fresh)}"


def test_every_verdict_states_the_phrase_it_rests_on():
    for verdict in CURATION:
        assert verdict.category in CATEGORIES
        assert len(verdict.basis) > 25, (
            f"{verdict.label}'s verdict has no quotable basis")


def test_a_direction_may_only_come_from_a_stated_outcome():
    """The type refuses to attach a direction to a construct list."""
    with pytest.raises(ValueError, match="functional outcome"):
        Verdict("X1A", "no_phenotype", "LoF", "a" * 30)
    with pytest.raises(ValueError, match="unknown category"):
        Verdict("X1A", "vibes", None, "a" * 30)


# ---------------------------------------------------------- the answer

def test_five_of_thirty_five_carry_a_direction():
    """The number Round 57 exists to produce."""
    report = summary()
    assert report["total"] == 35
    assert report["directional"] == 5
    assert report["distinct_directional"] == 5
    assert {v.label for v in directional()} == {
        "D1959A", "D2018A", "L2131A", "R2119A", "W2124A"}


def test_all_five_are_loss_of_channel_function_not_disease_direction():
    """They say "non-functional", which is not the curated set's LoF.

    An alanine-scanning mutant that does not conduct is a loss of channel
    function in a screen. The curated LoF label means a variant that causes a
    loss-of-function disease. Pooling them would be the evidence-level mistake
    `variant_sets` exists to prevent.
    """
    for verdict in directional():
        assert verdict.direction == "LoF"
        assert "non-functional" in verdict.basis
        assert verdict.label.endswith("A"), (
            "all five are alanine substitutions from a mutagenesis screen")


def test_the_five_unlock_no_within_position_pair():
    """The decisive measurement: hand-curation adds nothing to either design.

    None of the five sits at a position carrying any other variant, so they
    create no pair — and Round 54's count of one usable position is unchanged.
    """
    from piezo1.analysis.data_routes import _variants_by_position

    by_position, _ = _variants_by_position()
    for verdict in directional():
        position = int("".join(c for c in verdict.label[1:-1] if c.isdigit()))
        assert not by_position.get(position), (
            f"{verdict.label} now has a partner at {position}; the Round 54 "
            f"and Round 57 conclusions should be revisited")


# -------------------------------------------- the failure modes found

def test_the_harvest_extracted_a_substitution_from_a_different_protein():
    """V190P is STOML3, and no residue-identity gate could have caught it.

    The wild-type gate rejects 23% of raw hits and is the reason to trust the
    rest. It passed this one because position 190 is valine in PIEZO1 too. The
    class of error is recorded rather than filtered away, because a filter
    would hide how the gate can fail.
    """
    wrong = by_category("wrong_protein")
    assert len(wrong) == 1
    assert wrong[0].label == "V190P"
    assert "STOML3" in wrong[0].basis

    candidate = next(c for c in harvest().passing() if c.human_label == "V190P")
    assert "STOML3" in candidate.context
    assert candidate.wt == "V", "the gate passed it on a genuine wild-type match"


def test_clone_sequence_notes_are_not_tested_mutants():
    labels = {v.label for v in by_category("sequence_variant")}
    assert labels == {"V250A", "V394L", "R407G"}
    for verdict in by_category("sequence_variant"):
        assert "clone" in verdict.basis or "sequencing" in verdict.basis


def test_conductance_changes_are_held_not_admitted():
    """The same question Round 54 marked blocked, kept blocked here."""
    conductance = by_category("conductance_only")
    assert len(conductance) == 5
    assert all(v.direction is None for v in conductance)


def test_two_parsed_measurements_are_truncation_artefacts():
    """`measurements={'conductance': 7.0}` comes from '...56.7 pS' split badly.

    Round 45 reported two candidates carrying a measurement. Read in context
    both are fragments of a conductance list, so the count of candidates with
    usable measurements is zero, not two.
    """
    labels = {v.label for v in by_category("conductance_only")}
    assert {"V2116A", "L2118A"} <= labels
    for label in ("V2116A", "L2118A"):
        verdict = next(v for v in CURATION if v.label == label)
        assert "artefact" in verdict.basis or "truncat" in verdict.basis

    measured = [c for c in harvest().passing()
                if c.measurements and not c.already_curated]
    assert len(measured) == 2, "the harvest still reports two; both are artefacts"


def test_the_categories_account_for_every_candidate():
    report = summary()
    assert sum(report["by_category"].values()) == report["total"]
    for name, reason in CATEGORIES.items():
        assert len(reason) > 20, f"{name} is a category without an explanation"
