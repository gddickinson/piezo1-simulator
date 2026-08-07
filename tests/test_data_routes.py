"""What each remaining route to more phenotyped variants would yield.

The load-bearing test is the one that pins the *correction*: a within-position
comparison has one usable position, not the forty the Round 50 review claimed.
That review counted positions carrying more than one variant of any kind, which
is not what the design needs.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.data_routes import (DIRECTIONS, CurationTarget, Route,
                                         assess_routes,
                                         discriminating_positions,
                                         one_label_away, usable_positions)


# ------------------------------------------------- the correction

def test_only_one_position_supports_a_within_position_comparison():
    """Round 48's count was right; the Round 50 review's was not.

    Forty positions carry more than one variant. One has two or more missense
    variants that each carry a direction with no source conflict.
    """
    usable = usable_positions()
    assert list(usable) == [2456], (
        f"expected only R2456, got {sorted(usable)}")
    assert len(usable[2456]) == 4
    assert set(usable[2456].values()) == {"GoF", "LoF"}


def test_the_one_usable_position_is_also_discriminating():
    assert list(discriminating_positions()) == [2456]


def test_nonsense_and_insertion_variants_are_not_counted_as_substitutions():
    """The three filters that took forty down to one, checked individually.

    A stop codon truncates the protein and an insertion changes its length;
    neither can be scored against a substitution at the same position. Both
    were being counted.
    """
    from piezo1.analysis.data_routes import _missense

    assert _missense("R2456C", "R", "C")
    assert not _missense("Q1009*", "Q", "*"), "nonsense must not count"
    assert not _missense("E2496ELE", "E", "ELE"), "insertion must not count"
    assert not _missense("E756-", "E", "-"), "deletion must not count"
    assert not _missense("Q1009L", "Q", "L", kind="nonsense"), \
        "an explicit kind must override the label"


def test_a_position_with_conflicting_sources_is_excluded():
    """V598M is GoF in the curated set and LoF in ClinVar.

    An unresolved disagreement is not evidence, so the position is dropped
    rather than resolved by picking a side. `variant_sets.disagreements()`
    surfaces it for the reader, so it is excluded here, not hidden.
    """
    from piezo1.analysis.variant_sets import disagreements

    reported = disagreements()
    assert any(d["label"] == "V598M" for d in reported)
    assert 598 not in usable_positions()
    assert all(t.position != 598 for t in one_label_away())


# ---------------------------------------------------- the actionable list

def test_the_curation_list_is_specific_and_small():
    """"Not enough variants" becomes three named ones."""
    targets = one_label_away()
    assert targets, "the route must name what to curate"
    assert all(isinstance(t, CurationTarget) for t in targets)
    needed = {t.needs for t in targets}
    assert needed == {"M870V", "R1358C", "A2020V"}, needed
    for target in targets:
        assert target.existing_direction in DIRECTIONS
        assert target.needs != target.unlocks_with


def test_curating_the_whole_list_still_does_not_make_a_design():
    """The honest conclusion, pinned so it cannot be quietly outgrown.

    One position now, at most four if every target resolves — and two of the
    three are VUS precisely because the evidence was not found.
    """
    reachable = len(usable_positions()) + len(one_label_away())
    assert reachable <= 5, (
        "the within-position route has grown; Round 54's conclusion should be "
        "revisited rather than left as recorded")


# ------------------------------------------------------------- the routes

def test_the_two_routes_the_roadmap_named_are_already_spent():
    routes = {r.name: r for r in assess_routes()}
    for name in ("Population constraint (gnomAD)",
                 "Published supplementary tables (literature harvest)"):
        assert routes[name].status == "done"
        assert routes[name].yields == 0


def test_every_route_states_a_yield_a_cost_and_a_status():
    for route in assess_routes():
        assert isinstance(route, Route)
        assert route.status in {"done", "open", "blocked"}
        assert route.cost, f"{route.name} states no cost"
        assert len(route.note) > 40, f"{route.name}'s note is too vague"
        assert route.yields >= 0


def test_the_engineered_route_is_blocked_on_a_scientific_question():
    """Fifteen measured effects that are not measurements of direction."""
    route = next(r for r in assess_routes() if "engineered" in r.name)
    assert route.status == "blocked"
    assert route.yields == 15
    assert "conductance" in route.note or "selectivity" in route.note


def test_the_cheapest_open_route_declares_its_yield_an_upper_bound():
    route = next(r for r in assess_routes() if "one label" in r.name.lower())
    assert route.status == "open"
    assert "UPPER BOUND" in route.note
    assert route.yields == len(one_label_away())


def test_no_route_claims_to_solve_the_problem():
    """Round 47 measured 134 variants needed. Nothing here approaches it."""
    total = sum(r.yields for r in assess_routes() if r.status != "done")
    assert total < 50, (
        f"routes now total {total}; if this approaches 134 the feasibility "
        f"conclusion must be recomputed")


# ------------------------------- Round 62: the same count, per evidence level

def test_the_one_usable_position_is_at_the_strongest_evidence_level():
    """R2456's four variants are all electrophysiology, not inference.

    Worth knowing before concluding: the single discriminating position is not
    a weak one that a stricter design would discard. It is the project's
    best-evidenced position. The problem is that there is one of it.
    """
    from piezo1.analysis.data_routes import positions_by_evidence

    levels = positions_by_evidence()
    assert levels["measured"]["usable"] == [2456]
    assert levels["measured+disease_mechanism"]["usable"] == [2456], (
        "pooling in the inferred level adds no usable position")


def test_a_pair_is_no_stronger_than_its_weaker_half():
    """M870 can only ever make an inferred-level pair, so it is not counted
    at `measured` — the ceiling is 3 there, not 4."""
    from piezo1.analysis.data_routes import positions_by_evidence

    levels = positions_by_evidence()
    assert 870 not in levels["measured"]["reachable"], (
        "M870I is disease_mechanism; curating M870V cannot make it measured")
    assert 870 in levels["measured+disease_mechanism"]["reachable"]
    assert len(levels["measured"]["reachable"]) == 3
    assert len(levels["measured+disease_mechanism"]["reachable"]) == 4


def test_no_evidence_level_reaches_the_requirement():
    """The Round 62 answer: splitting by evidence makes the ceiling worse.

    Round 61 measured that even a predictor ordering nine pairs in ten
    correctly needs 8 positions. The best reachable count at any level is 4,
    and 3 if the test is to be confirmatory on measured evidence alone.
    """
    from piezo1.analysis.data_routes import evidence_summary

    report = evidence_summary()
    assert report["best_case_requirement"] >= 5
    assert report["reachable_anywhere"] < report["best_case_requirement"], (
        "a route has opened; Rounds 61 and 62 must be revisited")
    for level, counts in report["by_level"].items():
        assert counts["reachable"] < report["best_case_requirement"], level


def test_the_evidence_levels_are_not_pooled_by_accident():
    """The `measured` view must be a strict subset of the pooled one."""
    from piezo1.analysis.data_routes import positions_by_evidence

    levels = positions_by_evidence()
    strict = set(levels["measured"]["reachable"])
    pooled = set(levels["measured+disease_mechanism"]["reachable"])
    assert strict <= pooled
    assert strict != pooled, "pooling must add something, or the split is moot"
