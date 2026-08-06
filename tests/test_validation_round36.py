"""The Round 36 record, pinned.

A recorded result is never revised (protocol §6), so these tests exist to make a
silent change loud. They check the numbers the document states, that the
decision follows the rule that was fixed in advance, and that the write-up still
says what a null of this size is and is not entitled to claim.
"""

from __future__ import annotations

import json
import re

import pytest

from piezo1.config import DERIVED_DIR, PROJECT_ROOT

RESULT = DERIVED_DIR / "validation_round36.json"
DOCUMENT = PROJECT_ROOT / "docs" / "VALIDATION_ROUND36.md"
PREREGISTRATION = PROJECT_ROOT / "docs" / "PREREGISTRATION_ROUND36.md"


@pytest.fixture(scope="module")
def result():
    if not RESULT.exists():
        pytest.skip("run scripts/run_validation_round36.py first")
    return json.loads(RESULT.read_text())


def test_the_primary_result_is_what_was_recorded(result):
    primary = result["primary"]
    assert primary["n_gof"] == 19 and primary["n_lof"] == 15
    assert primary["cliffs_delta"] == pytest.approx(-0.249, abs=0.01)
    assert primary["p_value"] == pytest.approx(0.405, abs=0.02)
    assert primary["ci_low"] < 0 < primary["ci_high"], "the interval spans zero"
    assert primary["auroc"] == pytest.approx(0.625, abs=0.02)


def test_the_decision_follows_the_rule_fixed_in_advance(result):
    """p < 0.05 AND delta < 0 AND the CI excluding zero. Only one of three holds."""
    primary = result["primary"]
    assert result["decision"] == "fail_to_reject"
    assert primary["cliffs_delta"] < 0            # direction, yes
    assert primary["p_value"] >= 0.05             # significance, no
    assert primary["ci_high"] > 0                 # interval, no


def test_nothing_in_the_secondary_family_survives_correction(result):
    assert not result["any_secondary_significant"]
    assert len(result["secondary"]) == 6
    assert min(f["q_value"] for f in result["secondary"]) > 0.5


def test_an_endpoint_that_could_not_be_run_is_recorded_not_dropped(result):
    """"Could not be run" and "was not significant" are different statements."""
    untestable = {e["feature"] for e in result["untestable_endpoints"]}
    assert "foldx_ddg" in untestable
    for entry in result["untestable_endpoints"]:
        assert entry["values_available"] == 0
        assert entry["reason"]


def test_the_substitution_aware_predictor_beats_its_own_control(result):
    """Round 26's improvement is visible even though neither is significant.

    The volume-only predictor is the negative control pre-registered in §7. If
    the substitution-aware score ever stops beating it, Round 26's result has
    regressed and this test should say so.
    """
    primary = result["primary"]["cliffs_delta"]
    control = next(f for f in result["secondary"]
                   if f["feature"] == "mechanical_volume_only")["cliffs_delta"]
    assert primary < control, "the aware predictor should show the larger effect"
    assert abs(primary) > 5 * abs(control)


def test_every_exclusion_is_accounted_for(result):
    assert result["n_variants"] == 34
    assert result["dropped"]["excluded_by_name"] == 1      # V598M
    assert result["dropped"]["not_modelled"] == 11
    assert 34 + 1 + 11 == 46


def test_the_write_up_does_not_soften_the_null():
    text = re.sub(r"\s+", " ", DOCUMENT.read_text())
    assert "FAIL TO REJECT" in text
    assert "third null" in text.lower()
    # It must state the scope of the null, not merely report a p-value.
    assert "excludes a large effect" in text
    assert "does not exclude a medium one" in text
    # And it must not claim the trend as evidence.
    assert "is **not evidence**" in text or "not evidence" in text


def test_the_preregistration_came_first_and_says_what_it_committed_to():
    text = re.sub(r"\s+", " ", PREREGISTRATION.read_text())
    assert "committed before the comparison was executed" in text
    assert "CONFIRMATORY FOR A LARGE EFFECT AND EXPLORATORY BELOW IT" in text
    assert "V598M is excluded" in text
    # The two evidence levels must be handled explicitly, not pooled silently.
    assert "not pooled" in text.lower()
