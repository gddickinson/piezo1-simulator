"""The Round 22 exploratory test.

Pins the recorded numbers so a change to any predictor, cache or statistic
cannot move a published result silently — the same guard `test_validation.py`
puts on Round 7.

Nothing here re-runs the hypothesis or re-tests the labels. The result is
recorded; these assert the record is reproducible.
"""

import json
import pathlib

import numpy as np
import pytest

RESULT = pathlib.Path("data/derived/validation_round22.json")


@pytest.fixture(scope="module")
def result():
    if not RESULT.exists():
        pytest.skip("run scripts/run_validation_round22.py first")
    return json.loads(RESULT.read_text())


def test_design_is_declared_exploratory(result):
    """The declaration is load-bearing: it forbids reporting a decision."""
    assert "exploratory" in result["design"].lower()
    assert result["preregistration"] == "docs/PREREGISTRATION_ROUND22.md"


def test_counts_match_the_preregistration(result):
    counts = result["counts"]
    assert counts["curated"] == 68
    assert counts["directional"] == 39
    assert counts["included"] == 26
    assert counts["included_gof"] == 20
    assert counts["included_lof"] == 6
    assert counts["dropped_not_single_substitution"] == 13


def test_most_dropped_variants_are_loss_of_function(result):
    """The design's central problem, and it is biology not curation.

    Loss of function is commonly achieved by truncation, and no missense
    predictor can score a stop codon — so the exclusion removes 11 of 17 LoF
    variants and leaves a subset selected for *not* being truncating.
    """
    assert result["counts"]["dropped_lof"] == 11


def test_primary_endpoint_found_nothing(result):
    """FoldX ddG: point estimate opposite to the hypothesis, interval spans zero."""
    primary = result["primary"]
    assert primary["feature"] == "foldx_ddg"
    assert primary["cliffs_delta"] == pytest.approx(-0.211, abs=0.01)
    assert primary["ci_low"] < 0 < primary["ci_high"], "interval must span zero"
    assert primary["auroc"] == pytest.approx(0.395, abs=0.01)
    assert primary["mean_lof"] < primary["mean_gof"], (
        "recorded direction: LoF less destabilising than GoF")


def test_three_statistics_agree_on_direction(result):
    """Guards against a sign error in the one place it would matter."""
    primary = result["primary"]
    assert primary["cliffs_delta"] < 0
    assert primary["auroc"] < 0.5
    assert primary["mean_lof"] - primary["mean_gof"] < 0


def test_no_secondary_predictor_survives_correction(result):
    for row in result["bh"]:
        assert not row["significant"], f"{row['name']} unexpectedly significant"
    assert min(row["q"] for row in result["bh"]) > 0.4


def test_pathogenicity_predictors_behaved_as_predicted(result):
    """The pre-registration expected NO separation from these three.

    A single pathogenicity axis cannot express direction when both classes are
    pathogenic. Every interval spanning zero is the expectation holding.
    """
    by_name = {r["feature"]: r for r in result["secondary"] if "error" not in r}
    for name in ("alphamissense", "eve", "esm1b"):
        row = by_name[name]
        assert row["ci_low"] < 0 < row["ci_high"], f"{name} separated the classes"


def test_combining_features_added_nothing(result):
    combined = result["combined"]
    assert combined["auroc_leave_one_out"] == pytest.approx(0.535, abs=0.02)
    # Equal weights have no fitted parameters, so there is nothing to overfit.
    assert combined["optimism"] == pytest.approx(0.0, abs=1e-9)


def test_the_design_was_underpowered_and_says_so(result):
    power = result["power"]
    assert power["mde_80"] > 0.55, "80% power should need a very large effect"
    at_large = dict(zip([abs(d) for d in power["deltas"]], power["power"]))
    assert at_large[0.43] < 0.6, "power at a large effect was a coin toss"


def test_the_written_record_states_the_null(result):
    """The document must not soften what the numbers say."""
    import re
    raw = pathlib.Path("docs/VALIDATION_ROUND22.md").read_text()
    # Markdown wraps at 79 columns, so match against whitespace-normalised
    # text rather than requiring a phrase to sit on one line.
    text = re.sub(r"\s+", " ", raw).lower()
    assert "exploratory. not a validation." in text
    assert "opposite direction" in text
    assert "rules out essentially nothing" in text
    assert "no decision" in text
