"""The design analysis behind Round 47, and the guards that keep it honest.

Two things are being protected. First, that this module **runs no comparison** —
it reads a recorded effect size and asks a question about sample size, and a
test that let a phenotype comparison creep in here would let the project run a
fifth test without a pre-registration. Second, that the conclusion is a real
consequence of the numbers rather than a sentence someone wrote.
"""

from __future__ import annotations

import inspect

import pytest

from piezo1.analysis import feasibility
from piezo1.analysis.design import minimum_detectable_effect
from piezo1.analysis.feasibility import (assess, modelling_survival,
                                         observed_effect, observed_split,
                                         recorded_round)


@pytest.fixture(scope="module")
def report():
    return assess(n_simulations=600)


# --------------------------------------------------- no comparison happens

def test_this_module_runs_no_phenotype_comparison():
    """The load-bearing guard on the whole round.

    The roadmap said explicitly: do not run the comparison unless the
    pre-registration protocol is followed first. A design analysis is allowed
    because it uses only the *recorded* effect. If someone later imports a
    statistic here and computes a fresh delta, that would be a fifth test run
    without a pre-registration, and it would arrive wearing the word
    'feasibility'.
    """
    source = inspect.getsource(feasibility)
    # Reading ``record.cliffs_delta`` off the recorded result is exactly what
    # this module should do; *importing* a statistic to compute a fresh one is
    # what must not happen. So the guard is on imports, not on words.
    for banned in ("analysis.validation", "from .validation",
                   "from .substitution", "from .variant_impact",
                   "cliffs_delta(", "permutation_test(", "auroc("):
        assert banned not in source, (
            f"{banned!r} appears in feasibility.py — this module must not run "
            f"a comparison, only reason about sample size")
    # And it must genuinely read the recorded result rather than restate it.
    assert "VALIDATION_RECORD" in source


def test_the_effect_size_is_the_recorded_one_not_a_fresh_measurement():
    """It must match what docs/VALIDATION_ROUND36.md records.

    Derived from VALIDATION_RECORD rather than copied, so the two cannot drift
    apart; this checks the document agrees as well.
    """
    from pathlib import Path

    assert observed_effect() == pytest.approx(-0.249)
    assert recorded_round().round == 36
    doc = Path(__file__).resolve().parents[1] / "docs" / "VALIDATION_ROUND36.md"
    if doc.exists():
        assert "0.249" in doc.read_text(), (
            "the recorded effect and the constant here have diverged")


# ------------------------------------------------------------ the finding

def test_the_observed_effect_is_not_detectable_at_any_reachable_size(report):
    """Round 47's answer: no, and not by a small margin.

    Every reachable scenario has a minimum detectable effect larger than the
    effect actually observed. This is the difference between 'we need more
    data' and 'the data that could exist is not enough', and only the second
    tells you not to run the test.
    """
    assert not report.achievable
    for scenario in report.scenarios:
        if scenario.reachable:
            assert scenario.minimum_detectable > abs(observed_effect()), (
                f"{scenario.label} would now detect the observed effect — "
                f"the conclusion of Round 47 has changed and the docs must too")


def test_even_the_optimistic_ceiling_is_a_coin_flip(report):
    """Power at the ceiling is about 0.5, against the 0.8 a test should have."""
    ceiling = report.get("optimistic ceiling")
    assert ceiling is not None
    assert 0.40 < ceiling.power_at_observed < 0.62, ceiling.summary()
    assert ceiling.n < report.required_n / 1.8, (
        "the ceiling has moved within reach of the requirement")


def test_the_requirement_is_about_four_times_what_exists(report):
    """134-ish against 34. Pinned loosely, since it is simulation-based."""
    assert 110 <= report.required_n <= 165
    today = report.get("today")
    assert today.n == 34
    assert report.required_n > 3 * today.n


def test_required_n_actually_delivers_the_power_it_claims(report):
    """The requirement is checked against the power curve, not just asserted."""
    required = report.get("required for 80% power")
    assert required.power_at_observed == pytest.approx(0.80, abs=0.08)
    assert not required.reachable


# ------------------------------------------------- the ceiling is a ceiling

def test_the_ceiling_is_optimistic_in_every_direction(report):
    """It assumes things that are not true today, and must keep assuming them.

    If this ever became a *realistic* estimate rather than a generous one, the
    conclusion would be weaker than stated. The ceiling counts every harvest
    candidate as though it already carried a measured direction — Round 45
    found that none of them do.
    """
    meta = report.meta
    assert meta["fresh_harvest_candidates"] > 0
    assert report.ceiling_n > meta["directional_variants"], (
        "the ceiling must exceed what exists, or it is not a ceiling")
    projected = ((meta["directional_variants"] + meta["fresh_harvest_candidates"])
                 * modelling_survival())
    assert report.ceiling_n == int(projected)
    # And the survival rate is applied, not ignored.
    assert modelling_survival() < 1.0
    assert report.ceiling_n < (meta["directional_variants"]
                               + meta["fresh_harvest_candidates"])


def test_more_variants_always_helps_but_not_enough(report):
    """Monotonicity — a sanity check on the design functions themselves."""
    ordered = sorted(report.scenarios, key=lambda s: s.n)
    mdes = [s.minimum_detectable for s in ordered]
    powers = [s.power_at_observed for s in ordered]
    assert mdes == sorted(mdes, reverse=True), "MDE must fall as n rises"
    assert powers == sorted(powers), "power must rise with n"


def test_scenarios_agree_with_the_design_module(report):
    """No re-implementation: the numbers come from analysis.design."""
    for scenario in report.scenarios:
        assert scenario.minimum_detectable == pytest.approx(
            minimum_detectable_effect(scenario.n_a, scenario.n_b), abs=1e-9)
        assert scenario.n_a + scenario.n_b == scenario.n


def test_the_split_is_the_one_the_round_actually_had(report):
    """Not a round number someone chose — 19 GoF of 34."""
    record = recorded_round()
    assert observed_split() == pytest.approx(
        record.n_gof / (record.n_gof + record.n_lof))
    assert report.meta["split"] == pytest.approx(observed_split())
    today = report.get("today")
    assert (today.n_a, today.n_b) == (record.n_gof, record.n_lof)


def test_a_skewed_split_is_never_treated_as_better(report):
    """The projection assumes the imbalance persists, which is the harder case."""
    balanced = minimum_detectable_effect(30, 30)
    skewed = minimum_detectable_effect(45, 15)
    assert skewed > balanced, (
        "an imbalanced set must be harder, or the ceiling is flattered")


def test_summary_states_the_answer_in_words(report):
    text = report.summary()
    assert "0.249" in text and str(report.ceiling_n) in text
    assert "False" in text
