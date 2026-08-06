"""Statistics for the blind variant test.

These check the instruments, not the outcome. The Round 7 result itself is a
recorded fact in ``docs/VALIDATION.md`` and is pinned here only so that a later
change to the predictor cannot silently alter a published number without a test
going red.
"""

import json
from pathlib import Path

import numpy as np
import pytest

from piezo1.analysis.validation import (auroc, bootstrap_cliffs_delta,
                                        cliffs_delta, interpret_delta,
                                        permutation_test)


# --------------------------------------------------------------------------
# Cliff's delta
# --------------------------------------------------------------------------

def test_cliffs_delta_extremes():
    assert cliffs_delta([1, 2, 3], [4, 5, 6]) == pytest.approx(-1.0)
    assert cliffs_delta([4, 5, 6], [1, 2, 3]) == pytest.approx(1.0)
    assert cliffs_delta([1, 2, 3], [1, 2, 3]) == pytest.approx(0.0)


def test_cliffs_delta_ties_contribute_nothing():
    assert cliffs_delta([1, 1, 1], [1, 1, 1]) == pytest.approx(0.0)
    # Half above, half equal: only the strict comparisons count.
    assert cliffs_delta([1, 2], [1, 1]) == pytest.approx(0.5)


def test_delta_interpretation_thresholds():
    assert interpret_delta(0.05) == "negligible"
    assert interpret_delta(0.2) == "small"
    assert interpret_delta(0.4) == "medium"
    assert interpret_delta(0.8) == "large"
    assert interpret_delta(-0.8) == "large"      # magnitude, not sign


# --------------------------------------------------------------------------
# AUROC
# --------------------------------------------------------------------------

def test_auroc_known_cases():
    assert auroc([1, 2, 3, 4], [False, False, True, True]) == pytest.approx(1.0)
    assert auroc([4, 3, 2, 1], [False, False, True, True]) == pytest.approx(0.0)
    # Positives {1, 3} against negatives {2, 4}: only the (3, 2) pair ranks
    # correctly, so 1 of 4 pairs -> 0.25, not 0.5. Symmetric-looking labels are
    # not the same as chance performance.
    assert auroc([1, 2, 3, 4], [True, False, True, False]) == pytest.approx(0.25)
    # Genuine chance: positives {1, 4} against negatives {2, 3} -> 2 of 4.
    assert auroc([1, 2, 3, 4], [True, False, False, True]) == pytest.approx(0.5)


def test_auroc_of_a_constant_predictor_is_exactly_half():
    """Ties must be averaged, or a useless predictor scores 0 or 1."""
    assert auroc([7, 7, 7, 7], [True, True, False, False]) == pytest.approx(0.5)


def test_auroc_undefined_without_both_classes():
    assert np.isnan(auroc([1, 2, 3], [True, True, True]))


# --------------------------------------------------------------------------
# Permutation test
# --------------------------------------------------------------------------

def test_permutation_p_is_never_zero():
    """The (r+1)/(n+1) convention: a finite shuffle count cannot give p = 0."""
    r = permutation_test(np.full(10, -100.0), np.full(10, 100.0),
                         n_permutations=200, alternative="less")
    assert r.p_value > 0
    assert r.p_value == pytest.approx(1 / 201, rel=0.01)


def test_permutation_on_null_data_is_not_significant():
    rng = np.random.default_rng(0)
    r = permutation_test(rng.normal(size=50), rng.normal(size=50),
                         n_permutations=2000, alternative="less", seed=1)
    assert 0.02 < r.p_value < 0.98


def test_permutation_detects_a_real_shift():
    rng = np.random.default_rng(2)
    r = permutation_test(rng.normal(-2, 1, 40), rng.normal(0, 1, 40),
                         n_permutations=2000, alternative="less", seed=3)
    assert r.p_value < 0.01
    assert r.observed < 0


def test_permutation_direction_matters():
    rng = np.random.default_rng(4)
    a, b = rng.normal(-2, 1, 30), rng.normal(0, 1, 30)
    less = permutation_test(a, b, 2000, alternative="less", seed=5)
    greater = permutation_test(a, b, 2000, alternative="greater", seed=5)
    assert less.p_value < 0.05 < greater.p_value


def test_permutation_rejects_tiny_groups():
    with pytest.raises(ValueError, match="at least two"):
        permutation_test([1.0], [2.0, 3.0])


def test_bootstrap_interval_brackets_the_point_estimate():
    rng = np.random.default_rng(6)
    a, b = rng.normal(0, 1, 30), rng.normal(1, 1, 30)
    e = bootstrap_cliffs_delta(a, b, n_bootstrap=500, seed=7)
    assert e.ci_low <= e.delta <= e.ci_high
    assert e.interpretation == interpret_delta(e.delta)


def test_bootstrap_interval_excludes_zero_for_a_clear_effect():
    rng = np.random.default_rng(8)
    e = bootstrap_cliffs_delta(rng.normal(-3, 0.5, 40), rng.normal(3, 0.5, 40),
                               n_bootstrap=500, seed=9)
    assert e.excludes_zero
    assert e.delta < 0


# --------------------------------------------------------------------------
# The recorded Round 7 result
# --------------------------------------------------------------------------

def test_round7_result_is_reproducible():
    """Pin the published null result.

    If a change to the predictor moves these numbers, this test goes red and
    the change has to be reported as a new hypothesis rather than quietly
    replacing a published figure.
    """
    path = Path("data/derived/validation_round7.json")
    if not path.exists():
        pytest.skip("run scripts/run_validation.py first")
    d = json.loads(path.read_text())

    assert d["counts"]["included"] == 25
    assert d["primary"]["n_gof"] == 16 and d["primary"]["n_lof"] == 9
    assert d["decision"] == "H0 not rejected"
    assert d["primary"]["permutation"]["p_value"] > 0.05
    assert d["primary"]["effect"]["interpretation"] == "negligible"
    # The confidence interval must span zero for a null result.
    assert d["primary"]["effect"]["ci_low"] < 0 < d["primary"]["effect"]["ci_high"]
    assert 0.4 < d["primary"]["auroc"] < 0.7
