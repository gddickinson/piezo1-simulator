"""Pins the Round 48 record so it cannot decay or be softened.

The load-bearing test here is not the null itself — it is §4's measurement that
a wild-type positional feature has **zero** within-position variance. That is
the reason a positive result would not have been the predictor the project
wants, and it was stated in the pre-registration before the numbers existed.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from piezo1.analysis.substitution import variance_decomposition

ROOT = Path(__file__).resolve().parents[1]
RECORD = ROOT / "data" / "derived" / "validation_round48.json"
DOCUMENT = ROOT / "docs" / "VALIDATION_ROUND48.md"
PREREGISTRATION = ROOT / "docs" / "PREREGISTRATION_ROUND48.md"


@pytest.fixture(scope="module")
def record():
    if not RECORD.exists():
        pytest.skip("Round 48 has not been run; "
                    "python scripts/run_validation_round48.py")
    return json.loads(RECORD.read_text())


# ------------------------------------------------------- the numbers

def test_the_primary_numbers_are_what_the_document_states(record):
    primary = record["primary"]
    assert primary["endpoint"] == "relative_sasa"
    assert primary["n_lof"] == 14 and primary["n_gof"] == 16
    assert primary["cliffs_delta"] == pytest.approx(0.036, abs=0.005)
    assert primary["p_value"] == pytest.approx(0.509, abs=0.02)
    assert primary["ci_low"] < 0 < primary["ci_high"]


def test_the_decision_follows_the_rule_fixed_in_advance(record):
    """All three clauses fail here, so the null is not a close call."""
    clauses = record["decision_clauses"]
    assert clauses["p < 0.05"] is False
    assert clauses["delta negative"] is False
    assert clauses["CI excludes zero"] is False
    assert record["decision"] == "fail_to_reject"


def test_the_round_is_exploratory_by_its_own_threshold(record):
    """|delta| 0.036 against the 0.495 fixed in section 6."""
    assert record["status"] == "exploratory"
    assert abs(record["primary"]["cliffs_delta"]) < record["confirmatory_threshold"]


def test_nothing_survives_correction(record):
    assert record["any_secondary_significant"] is False
    assert min(f["q_value"] for f in record["secondary"]) > 0.5


# ------------------------------------------- the ceiling, which is the point

def test_within_position_variance_is_exactly_zero(record):
    """§2's pre-committed ceiling, measured.

    Not 'small' — zero, because the feature never sees the substitution. This
    is what makes a positive result unusable as a variant-direction predictor,
    and it was written down before the test ran.
    """
    share = record["variance_share"]
    assert share["within_fraction"] == pytest.approx(0.0, abs=1e-12)
    assert share["between_fraction"] == pytest.approx(1.0, abs=1e-12)


def test_the_one_shared_position_gives_all_four_variants_one_value(record):
    """R2456: three GoF and one LoF, indistinguishable by construction."""
    multi = record["variance_share"]["multi_variant_positions"]
    assert list(multi) == ["2456"], (
        "the curated set has gained a second multi-variant position — the "
        "within-position design Round 48 could not run may now be possible")
    group = multi["2456"]
    assert len(group) == 4
    assert {d for _, d, _ in group} == {"GoF", "LoF"}
    values = {round(v, 12) for _, _, v in group}
    assert len(values) == 1, "the feature must be identical at one position"


def test_the_decomposition_would_detect_within_position_variance():
    """The measurement is capable of a non-zero answer — a control on §4.

    Zero is only meaningful if the instrument could have said otherwise.
    """
    positions = np.array([1, 1, 2, 2, 3, 3])
    identical = np.array([0.5, 0.5, 0.2, 0.2, 0.9, 0.9])
    varying = np.array([0.5, 0.1, 0.2, 0.8, 0.9, 0.3])
    assert variance_decomposition(positions, identical).within_fraction == \
        pytest.approx(0.0, abs=1e-12)
    assert variance_decomposition(positions, varying).within_fraction > 0.3


# ------------------------------------------------- the negative control

def test_the_negative_control_beats_every_real_endpoint(record):
    """The diagnostic that decides how to read this round.

    A feature pre-registered *because* it should be meaningless has a larger
    effect than every endpoint with a mechanism behind it. So the spread across
    endpoints here is noise, and any single large result would not have been
    distinguishable from it.
    """
    by_name = {f["endpoint"]: f for f in record["secondary"]}
    control = abs(by_name["distance_to_axis"]["cliffs_delta"])
    assert control == pytest.approx(0.268, abs=0.01)
    mechanistic = ["conservation", "prs_gate_response", "gating_amplitude",
                   "distance_to_gate"]
    for name in mechanistic:
        assert abs(by_name[name]["cliffs_delta"]) < control, (
            f"{name} now exceeds the negative control; the round's "
            f"interpretation in section 5 no longer holds")
    assert abs(record["primary"]["cliffs_delta"]) < control


def test_distance_to_gate_separates_nothing(record):
    """The endpoint with the clearest mechanical story is the flattest."""
    by_name = {f["endpoint"]: f for f in record["secondary"]}
    assert abs(by_name["distance_to_gate"]["cliffs_delta"]) < 0.01


# -------------------------------------------------------- the write-up

def test_the_preregistration_committed_the_ceiling_before_the_run():
    text = PREREGISTRATION.read_text()
    assert "0% within-position variance" in text
    assert "2456" in text
    # The exclusion and the power threshold were both fixed in advance.
    assert "0.495" in text
    assert "EXPLORATORY" in text


def test_the_writeup_does_not_soften_the_null():
    text = DOCUMENT.read_text()
    assert "FAIL TO REJECT" in text
    assert "EXPLORATORY" in text
    assert "1.000000" in text and "0.000000" in text
    assert "five nulls" in text.lower()
    for banned in ("trend toward", "trending", "approaching significance",
                   "marginally significant", "suggests that LoF"):
        assert banned not in text.lower(), f"softening language: {banned!r}"


def test_the_writeup_warns_against_reading_the_medians():
    """The medians reverse; the effect size says they do not mean it."""
    text = DOCUMENT.read_text()
    assert "0.219" in text and "0.143" in text
    assert "rank statistic" in text
