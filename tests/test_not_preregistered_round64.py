"""The decision not to run the within-position test, and its enforcement.

Round 64 was conditional: pre-register *only if* Rounds 61–63 leave a design
with adequate power. They do not — 8 positions required at an implausibly good
predictor, 1 available, 3–4 reachable. So the round records a refusal instead.

A refusal is only worth writing down if something makes it stick. These tests
are that: the arithmetic must still hold, no within-position comparison may
exist in the codebase, and the count of discriminating positions is ratcheted so
the question **resurfaces automatically** if the data ever changes rather than
depending on someone finding the document.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from piezo1.analysis.data_routes import (discriminating_positions,
                                         positions_by_evidence)
from piezo1.analysis.feasibility import paired_positions_required

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "NOT_PREREGISTERED_ROUND64.md"


# ------------------------------------------------- the arithmetic still holds

def test_the_requirement_still_exceeds_what_exists():
    """The condition the round was gated on, re-checked rather than assumed."""
    best_case = paired_positions_required(0.8, n_simulations=1500).positions
    available = len(discriminating_positions())
    assert available < best_case, (
        f"{available} discriminating positions against a best-case requirement "
        f"of {best_case} — if this has flipped, Round 64's refusal must be "
        f"revisited and the test pre-registered")


def test_no_evidence_level_reaches_the_requirement():
    best_case = paired_positions_required(0.8, n_simulations=1500).positions
    for level, counts in positions_by_evidence().items():
        assert len(counts["reachable"]) < best_case, level


def test_a_sign_test_on_one_pair_cannot_reject():
    """Section 3's claim, checked rather than asserted.

    With a single correctly-ordered pair the one-sided p is 0.5 exactly, so
    running it "exploratorily" has only one possible outcome.
    """
    from scipy import stats

    assert 1 - stats.binom.cdf(0, 1, 0.5) == pytest.approx(0.5)
    # Even four perfect pairs cannot reach the conventional threshold.
    assert 1 - stats.binom.cdf(3, 4, 0.5) > 0.05


# --------------------------------------------- the ratchet that reopens it

def test_the_discriminating_count_has_not_grown():
    """Fails if a new discriminating position appears — by design.

    This is the mechanism that makes the refusal self-revoking. A later round
    does not have to remember `NOT_PREREGISTERED_ROUND64.md`; the suite tells
    it the situation changed.
    """
    positions = discriminating_positions()
    assert set(positions) == {2456}, (
        f"discriminating positions are now {sorted(positions)} — Round 64's "
        f"refusal was based on there being one; re-read "
        f"docs/NOT_PREREGISTERED_ROUND64.md before proceeding")


# ------------------------------------- no such comparison exists in the code

def test_no_within_position_comparison_is_implemented():
    """The standing instruction, enforced on the codebase rather than trusted.

    The temptation is concrete: one position, four variants, and the ordering
    is a single line of code. A test is cheaper than discipline.
    """
    # Calibration first: `feasibility.py` simulates a sign test to compute the
    # REQUIRED sample size and imports the discriminating positions to compare
    # against. That is the design analysis, not the comparison, and a detector
    # keyed on those two things alone flags it — which it did on the first run.
    # What a real comparison additionally needs is variant SCORES.
    scoring = ("CouplingScore", "VariantImpactModel", "gating_cost_change")
    offenders = []
    for path in list((ROOT / "piezo1").rglob("*.py")) + \
                list((ROOT / "scripts").rglob("*.py")):
        text = path.read_text()
        statistic = any(k in text for k in ("binom", "sign_test", "wilcoxon"))
        if ("discriminating_positions" in text and statistic
                and any(k in text for k in scoring)):
            offenders.append(str(path.relative_to(ROOT)))
    assert not offenders, (
        f"a within-position comparison may be implemented in {offenders}; "
        f"Round 64 declined to pre-register one")


def test_the_feasibility_module_still_runs_no_comparison():
    import inspect

    from piezo1.analysis import feasibility

    source = inspect.getsource(feasibility)
    for banned in ("cliffs_delta(", "permutation_test(", "auroc(",
                   "CouplingScore", "VariantImpactModel"):
        assert banned not in source


# ------------------------------------------------------------ the document

def test_the_document_records_a_refusal_not_a_hypothesis():
    text = DOC.read_text()
    assert "not a pre-registration" in text.lower()
    assert "decision not to run" in text.lower()
    for word in ("hypothesis", "H₁", "H1"):
        assert f"## {word}" not in text, (
            "a refusal must not be shaped like a pre-registration")


def test_the_document_states_the_numbers_that_forced_it():
    flowed = " ".join(DOC.read_text().split())
    assert "8" in flowed and "102" in flowed
    assert "99.8%" in flowed, "the reason pairing was attractive at all"
    assert "R2456" in flowed


def test_the_document_says_why_exploratory_is_not_a_way_round_it():
    flowed = " ".join(DOC.read_text().split()).lower()
    assert "exploratory" in flowed
    assert "0.5" in flowed
    assert "not a blind test" in flowed


def test_the_document_carries_a_standing_instruction_and_names_its_guard():
    text = DOC.read_text()
    assert "Do not run a within-position comparison" in text
    assert "test_not_preregistered_round64" in text, (
        "the document must name the test that makes it self-revoking")
