"""Whether a conductance change may stand for a direction. It may not.

Round 54 marked the fifteen engineered variants `blocked` on a scientific
question rather than a data one. Round 63 answers it, and the answer rests on
this project's own curated annotations rather than on caution: two of the
fifteen explicitly demonstrate that the axes dissociate.

The round is then settled by a measurement rather than by the reasoning —
admitting the five that *are* admissible adds zero discriminating positions.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.engineered import (AXES, ENGINEERED, EngineeredVerdict,
                                        admitting_would_add, by_axis, decision,
                                        mechanosensitivity_variants)


# ------------------------------------- the verdicts match the curated record

def test_every_engineered_variant_is_classified():
    """A verdict per variant, none invented and none skipped."""
    from piezo1.core.annotations import load_annotations

    curated = {v.label for v in load_annotations("human").variants
               if v.classification == "engineered"}
    judged = {v.label for v in ENGINEERED}
    assert judged == curated, (
        f"unjudged: {sorted(curated - judged)}; invented: {sorted(judged - curated)}")


def test_each_verdict_rests_on_the_curated_functional_effect():
    """The basis must be a real phrase from the annotation, not a summary.

    A classification that paraphrases loosely could put a variant on whichever
    axis the argument needed.
    """
    from piezo1.core.annotations import load_annotations

    curated = {v.label: (v.functional_effect or "").lower()
               for v in load_annotations("human").variants}
    for verdict in ENGINEERED:
        effect = curated[verdict.label]
        head = verdict.basis.lower().split(",")[0].split(" while ")[0][:24]
        assert head in effect, (
            f"{verdict.label}: basis {verdict.basis!r} is not in the curated "
            f"effect {effect!r}")


def test_an_unknown_axis_is_refused():
    with pytest.raises(ValueError, match="unknown axis"):
        EngineeredVerdict("X1A", "vibes", "something")


def test_only_one_axis_defines_a_direction():
    directional = [name for name, (_, is_dir) in AXES.items() if is_dir]
    assert directional == ["mechanosensitivity"]


# ------------------------------------------- the evidence for the refusal

def test_the_dissociation_is_shown_by_the_projects_own_annotations():
    """Two variants make the argument; without them this would be an opinion.

    A2078W: agonist response lost, stretch retained.
    KKKK2166-: inactivation removed, mechanics unchanged.
    """
    agonist = next(v for v in ENGINEERED if v.label == "A2078W")
    assert "retained" in agonist.basis.lower()
    assert not agonist.is_a_direction

    inactivation = next(v for v in ENGINEERED if v.label == "KKKK2166-")
    assert "without changing mechanical" in inactivation.basis.lower()
    assert not inactivation.is_a_direction


def test_permeation_changes_are_not_directions():
    permeation = by_axis("permeation_only")
    assert len(permeation) == 5
    assert all(not v.is_a_direction for v in permeation)
    assert {v.label for v in permeation} >= {"E2117A", "E2117D", "E2117Q"}


def test_structural_probes_are_not_phenotypes():
    probes = by_axis("structural_probe")
    assert {v.label for v in probes} == {"Y2438C", "M2441C"}
    assert all(not v.is_a_direction for v in probes)


def test_the_decision_is_recorded_with_its_reason():
    verdict = decision()
    assert verdict["may_permeation_stand_for_direction"] is False
    assert "dissociate" in verdict["why"]
    assert "A2078W" in verdict["why"] and "KKKK2166-" in verdict["why"]


# ---------------------------------- the measurement that settles the round

def test_five_are_admissible_in_principle():
    """Refusing the permeation ones is not refusing all fifteen."""
    admissible = mechanosensitivity_variants()
    assert {v.label for v in admissible} == {
        "S1335A", "S1335V", "A1718W", "P2113A", "S2446E"}
    assert all(v.is_a_direction for v in admissible)


def test_admitting_them_adds_no_discriminating_position():
    """The reason the question is moot for the design, measured not argued."""
    added = admitting_would_add()
    assert added["discriminating_positions_gained"] == []
    assert decision()["would_add_discriminating_positions"] == 0


def test_position_1335_is_a_pair_but_not_a_discriminating_one():
    """Both S1335A and S1335V raise the threshold — same direction.

    This is the near miss worth pinning: it is the only engineered pair, and a
    within-position test needs the two variants to disagree.
    """
    added = admitting_would_add()
    assert 1335 in added["same_direction_pairs"]

    at_1335 = [v for v in mechanosensitivity_variants()
               if v.label.startswith("S1335")]
    assert len(at_1335) == 2
    for verdict in at_1335:
        assert "threshold" in verdict.basis or "force" in verdict.basis


def test_admission_would_not_pool_evidence_levels():
    """The Round 63 validation condition: `variant_sets` still refuses."""
    from piezo1.analysis.variant_sets import EVIDENCE_LEVELS

    assert "engineered" not in EVIDENCE_LEVELS, (
        "admitting them must add a level, never widen an existing one")
    assert decision()["note"].count("evidence level") >= 1
