"""What a within-position design would need, against what exists.

Round 47 closed the across-position route: 134 variants needed, 59 reachable.
The within-position route was the one it left open, because pairing removes the
between-position variance that consumed 99.8% of Round 7's predictor.

Round 61 costs it, and it is closed too. Even at an implausibly large paired
effect the design needs **8** shared positions; Round 54 measured that the
curated and ClinVar sets together contain **one**.

**No comparison is run here**, as in Round 47. The paired effect size cannot be
measured without running the test the pre-registration protocol forbids until a
design is registered, so the requirement is reported across a range of
hypothetical effects and the reader compares it with what exists.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.feasibility import (PairedRequirement,
                                         paired_feasibility,
                                         paired_positions_required,
                                         shared_positions_available)


# ------------------------------------------------- calibrate the instrument

def test_a_coin_flip_predictor_is_never_detectable():
    """delta = 0 must return None, not the search bound.

    Returning n_max would look like an answer — the same failure the kinetics
    calibration was fixed for in an earlier round.
    """
    result = paired_positions_required(0.0, n_max=60, n_simulations=800)
    assert result.positions is None
    assert "not detectable" in result.summary()


def test_a_perfect_predictor_needs_very_few_positions():
    """The other end of the scale, where the answer is known by inspection.

    If the predictor orders every pair correctly, a handful of positions is
    enough — this is what a sign test can do at p = 1.
    """
    result = paired_positions_required(1.0, n_simulations=800)
    assert result.positions is not None
    assert result.positions <= 6


def test_the_requirement_falls_as_the_effect_grows():
    sizes = [paired_positions_required(d, n_simulations=800).positions
             for d in (0.3, 0.5, 0.7, 0.9)]
    assert all(s is not None for s in sizes)
    assert sizes == sorted(sizes, reverse=True), sizes


def test_an_impossible_effect_is_refused():
    with pytest.raises(ValueError, match="outside"):
        paired_positions_required(-1.5)


# --------------------------------------------------------- the answer

def test_even_an_implausible_effect_needs_more_positions_than_exist():
    """The finding. The route Round 47 left open is closed too."""
    for requirement in paired_feasibility(n_simulations=1500):
        assert not requirement.reachable, (
            f"{requirement.summary()} — if this is now reachable the Round 61 "
            f"conclusion must be revisited")
        assert requirement.positions > shared_positions_available()


def test_the_best_case_still_needs_several_positions():
    """At delta 0.8 — 90% correct ordering — it is 8 against 1."""
    result = paired_positions_required(0.8, n_simulations=1500)
    assert 5 <= result.positions <= 12, result.summary()


def test_matching_the_across_position_effect_is_no_cheaper():
    """Pairing only helps if it *raises* the effect; at the same delta it does not.

    At the across-position δ = 0.249 the paired design needs about 100
    positions — comparable to the 134 variants Round 47 costed. The argument
    for pairing was always that it would enlarge the effect, not that it needs
    fewer observations at the same one.
    """
    result = paired_positions_required(0.249, n_simulations=1500)
    assert 80 <= result.positions <= 130, result.summary()


def test_the_available_count_matches_what_data_routes_measured():
    """One source of truth for how many shared positions exist."""
    from piezo1.analysis.data_routes import discriminating_positions

    assert shared_positions_available() == len(discriminating_positions())


def test_no_comparison_is_run():
    """The Round 47 constraint, re-applied.

    This module may reason about sample size; it may not measure a phenotype
    effect, because that test is not pre-registered.
    """
    import inspect

    from piezo1.analysis import feasibility

    source = inspect.getsource(feasibility)
    for banned in ("cliffs_delta(", "permutation_test(", "auroc(",
                   "VariantImpactModel", "CouplingScore"):
        assert banned not in source, (
            f"{banned} appears in feasibility.py — this is a design analysis")


def test_the_requirement_states_both_numbers():
    result = PairedRequirement(delta=0.5, positions=23)
    text = result.summary()
    assert "23" in text and str(shared_positions_available()) in text
