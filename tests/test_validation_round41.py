"""The Round 41 record: the fourth null, and the clause that caught it.

The primary p-value is 0.0477 — below the conventional threshold — and the
result is still a null, because the pre-registered rule required the confidence
interval to exclude zero as well. These tests hold that outcome, since it is
exactly the kind that later gets remembered as "we found something".
"""

from __future__ import annotations

import json
import re

import pytest

from piezo1.config import DERIVED_DIR, PROJECT_ROOT

RESULT = DERIVED_DIR / "validation_round41.json"
DOCUMENT = PROJECT_ROOT / "docs" / "VALIDATION_ROUND41.md"


@pytest.fixture(scope="module")
def result():
    if not RESULT.exists():
        pytest.skip("run scripts/run_validation_round41.py first")
    return json.loads(RESULT.read_text())


def test_the_primary_is_significant_by_p_and_still_a_null(result):
    """The whole point of the round's decision rule.

    p < 0.05 and delta < 0, but the interval contains zero, so the pre-registered
    conjunction fails. A rule written as "p < 0.05" alone would have produced
    this project's first positive on a p-value 0.0023 below threshold.
    """
    primary = result["primary"]
    assert primary["p_value"] < 0.05
    assert primary["cliffs_delta"] < 0
    assert primary["ci_low"] < 0 < primary["ci_high"], "the interval spans zero"
    assert result["decision"] == "fail_to_reject"


def test_the_primary_numbers_are_what_was_recorded(result):
    primary = result["primary"]
    assert primary["n_lof"] == 18 and primary["n_gof"] == 24
    assert primary["cliffs_delta"] == pytest.approx(-0.269, abs=0.01)
    assert primary["p_value"] == pytest.approx(0.048, abs=0.01)
    assert primary["median_lof"] < primary["median_gof"]


def test_the_negative_control_matches_the_predictor(result):
    """Pre-registered as showing nothing; it shows the same as the predictor.

    So the windowing that makes this a *regional* estimate contributes nothing,
    and the predictor cannot be distinguished from its own control — a stronger
    reason to disbelieve the primary than the interval alone.
    """
    control = next(f for f in result["secondary"]
                   if "negative control" in f["endpoint"])
    primary = result["primary"]
    assert abs(control["cliffs_delta"] - primary["cliffs_delta"]) < 0.10
    assert control["p_value"] < 0.2


def test_the_answer_does_not_depend_on_the_frozen_window(result):
    """+-10, +-25 and +-50 agree, so no result rests on the pre-registered 25."""
    deltas = [f["cliffs_delta"] for f in result["secondary"]
              if "local missense rate" in f["endpoint"]]
    deltas.append(result["primary"]["cliffs_delta"])
    assert max(deltas) - min(deltas) < 0.10


def test_nothing_survives_correction(result):
    assert not result["any_secondary_significant"]
    assert len(result["secondary"]) == 5
    assert min(f["q_value"] for f in result["secondary"]) > 0.1


def test_the_gene_is_not_constrained(result):
    """Recorded in the pre-registration before testing, because it sets the prior."""
    gene = result["gene_constraint"]
    assert gene["loeuf"] > 1.0 and not gene["lof_intolerant"]
    assert gene["oe_mis"] > 1.0 and not gene["missense_depleted"]
    assert gene["mis_z"] < -5


def test_the_ambiguous_position_was_excluded(result):
    """R2456 carries both directions; a position-level predictor cannot separate them."""
    assert result["excluded_positions"] == [2456]
    assert result["n_positions"] == 42


def test_the_write_up_does_not_promote_the_marginal_p():
    text = re.sub(r"\s+", " ", DOCUMENT.read_text())
    assert "FAIL TO REJECT" in text
    assert "fourth null" in text.lower()
    assert "contains zero" in text
    # It must name the negative control as the informative row.
    assert "negative control is the informative row" in text
    # And state what the null does and does not exclude.
    assert "does not exclude a medium one" in text
