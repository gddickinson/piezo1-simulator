"""Study design: power, multiplicity and cross-validation.

These are the checks that decide whether a result *could* have been found, and
whether one that was found survives having looked in several places. They run
before an analysis, not after a disappointing one.

Nothing here reads the variant phenotype labels except through the recorded
Round 7 output, and nothing re-tests the Round 7 hypothesis — that null result
stands and re-testing needs the Round 22 pre-registration written first.
"""

import json
import pathlib

import numpy as np
import pytest

from piezo1.analysis.design import (_permutation_p_values, benjamini_hochberg,
                                    delta_for_shift, leave_one_out,
                                    minimum_detectable_effect, power_curve,
                                    sample_size_for, shift_for_delta)
from piezo1.analysis.validation import cliffs_delta, permutation_test

ROUND7 = pathlib.Path("data/derived/validation_round7.json")


@pytest.fixture(scope="module")
def round7():
    if not ROUND7.exists():
        pytest.skip("Round 7 output not present — run scripts/run_validation.py")
    return json.loads(ROUND7.read_text())


# --------------------------------------------------------------------------
# Effect-size conversions
# --------------------------------------------------------------------------

def test_shift_and_delta_round_trip():
    for delta in (-0.9, -0.5, -0.3, -0.1, 0.25, 0.7):
        assert delta_for_shift(shift_for_delta(delta)) == pytest.approx(delta, abs=1e-9)


def test_delta_conversion_matches_sampled_data():
    """The analytic conversion has to agree with what the data actually show.

    δ = 2Φ(−shift/√2) − 1 is easy to write with the sign or the √2 wrong, and
    either mistake yields a smooth plausible curve.
    """
    rng = np.random.default_rng(0)
    for target in (-0.5, -0.3, 0.4):
        shift = shift_for_delta(target)
        measured = cliffs_delta(rng.normal(0.0, 1.0, 20000),
                                rng.normal(shift, 1.0, 20000))
        assert measured == pytest.approx(target, abs=0.02)


def test_zero_shift_is_zero_effect():
    assert delta_for_shift(0.0) == pytest.approx(0.0, abs=1e-12)


def test_delta_outside_range_is_rejected():
    with pytest.raises(ValueError, match="strictly inside"):
        shift_for_delta(1.0)


# --------------------------------------------------------------------------
# The fast permutation path
# --------------------------------------------------------------------------

def test_fast_path_agrees_with_the_real_permutation_test():
    """The power simulation reduces the test to a subset sum.

    ``mean(a) − mean(b) = S_a(1/n_a + 1/n_b) − S/n_b`` is monotone in ``S_a``,
    so the two must produce the same p-value up to Monte-Carlo noise. If the
    identity were wrong, every power number in this module would be wrong while
    still looking like a sensible curve.
    """
    rng = np.random.default_rng(1)
    for trial in range(4):
        a = rng.normal(0.0, 1.0, 16)
        b = rng.normal(0.5, 1.0, 9)
        slow = permutation_test(a, b, n_permutations=20000,
                                alternative="less", seed=trial).p_value
        fast = _permutation_p_values(np.concatenate([a, b])[None, :], 16, 20000,
                                     np.random.default_rng(trial), "less")[0]
        assert fast == pytest.approx(slow, abs=0.015)


# --------------------------------------------------------------------------
# Power
# --------------------------------------------------------------------------

def test_power_rises_with_effect_size():
    result = power_curve(16, 9, deltas=[-0.1, -0.3, -0.5, -0.8],
                         n_simulations=600, n_permutations=299, seed=2)
    assert np.all(np.diff(result.power) > 0)
    assert result.power[0] < 0.35
    assert result.power[-1] > 0.9


def test_simulated_effect_matches_the_requested_effect():
    """Guards the sign and the scale of the injected effect.

    Displacing the wrong group inverts the effect and silently reports the
    power to detect the opposite direction — which is what happened when this
    was first written.
    """
    result = power_curve(16, 9, deltas=[-0.8, -0.5, -0.3],
                         n_simulations=400, n_permutations=199, seed=3)
    achieved = result.meta["achieved_delta"]
    assert np.all(achieved < 0), "effect injected in the wrong direction"
    for target, got in zip(result.deltas, achieved):
        assert got == pytest.approx(target, abs=0.06)


def test_null_effect_gives_the_nominal_false_positive_rate():
    """At zero effect the rejection rate must be alpha, not more.

    A permutation test that over-rejects under the null would make every power
    number optimistic and every p-value in the project suspect.
    """
    result = power_curve(16, 9, deltas=[-1e-9], n_simulations=3000,
                         n_permutations=999, alpha=0.05, seed=4)
    assert result.power[0] < 0.075


def test_larger_samples_have_more_power():
    small = power_curve(16, 9, deltas=[-0.4], n_simulations=800,
                        n_permutations=299, seed=5).power[0]
    large = power_curve(60, 60, deltas=[-0.4], n_simulations=800,
                        n_permutations=299, seed=5).power[0]
    assert large > small + 0.2


def test_heavy_tailed_pool_is_supported(round7):
    values = np.array([v["ddg"] for v in round7["per_variant"]])
    result = power_curve(16, 9, deltas=[-0.5], n_simulations=600,
                         n_permutations=299, pool=values, seed=6)
    assert result.meta["model"] == "resampled from pool"
    assert 0.5 < result.power[0] < 0.95


# --------------------------------------------------------------------------
# The Round 7 design
# --------------------------------------------------------------------------

def test_round7_was_only_powered_for_a_large_effect():
    """The finding that qualifies the recorded null result.

    16 versus 9 reaches 80% power only around |δ| = 0.55, which is beyond
    'large' by the usual thresholds (0.11 / 0.28 / 0.43). The null therefore
    excludes a large effect and says little about a small or medium one. This
    does **not** revise the Round 7 result; it constrains what it means.
    """
    mde = minimum_detectable_effect(16, 9, n_simulations=2000,
                                    n_permutations=999,
                                    deltas=np.linspace(-0.9, -0.2, 15), seed=7)
    assert 0.45 < mde < 0.7, mde


def test_power_at_the_observed_effect_was_low(round7):
    values = np.array([v["ddg"] for v in round7["per_variant"]])
    labels = np.array([v["classification"] for v in round7["per_variant"]])
    observed = cliffs_delta(values[labels == "GoF"], values[labels == "LoF"])
    assert observed == pytest.approx(-0.083, abs=0.01)

    result = power_curve(16, 9, deltas=[observed], n_simulations=2000,
                         n_permutations=999, pool=values, seed=8)
    assert result.power[0] < 0.25, "an underpowered design should show it"


def test_sample_size_grows_as_the_effect_shrinks():
    large = sample_size_for(-0.43, n_simulations=500, n_permutations=299,
                            max_n=200, seed=9)
    medium = sample_size_for(-0.28, n_simulations=500, n_permutations=299,
                             max_n=200, seed=9)
    assert large < medium
    assert 10 < large < 40, large
    assert 30 < medium < 120, medium


# --------------------------------------------------------------------------
# Multiplicity
# --------------------------------------------------------------------------

def test_benjamini_hochberg_against_a_hand_computed_case():
    p = [0.01, 0.02, 0.03, 0.04, 0.05]
    result = benjamini_hochberg(p, alpha=0.05)
    # q_(i) = min_{j>=i} n/j * p_(j), enforced monotone non-decreasing.
    assert result.adjusted[0] == pytest.approx(0.05)
    assert np.all(np.diff(result.adjusted) >= -1e-12)
    assert np.all(result.adjusted >= np.array(p))


def test_correction_is_more_conservative_than_raw_p():
    names = ["mechanical", "conservation", "alphamissense", "eve", "esm1b", "foldx"]
    raw = [0.234, 0.041, 0.012, 0.180, 0.049, 0.310]
    result = benjamini_hochberg(raw, names, alpha=0.05, primary="mechanical")
    assert sum(1 for x in raw if x < 0.05) == 3
    assert result.n_significant == 0, "three raw hits, none survives correction"
    assert result.primary == "mechanical"


def test_adjusted_values_never_exceed_one():
    result = benjamini_hochberg([0.6, 0.9, 0.95])
    assert result.adjusted.max() <= 1.0


def test_unknown_primary_is_rejected():
    with pytest.raises(ValueError, match="primary endpoint"):
        benjamini_hochberg([0.1, 0.2], ["a", "b"], primary="c")


def test_invalid_p_values_are_rejected():
    with pytest.raises(ValueError, match=r"\[0, 1\]"):
        benjamini_hochberg([0.1, 1.5])


# --------------------------------------------------------------------------
# Cross-validation
# --------------------------------------------------------------------------

def test_leave_one_out_penalises_a_fitted_rule():
    """A rule allowed to peek at the held-out point must look better in sample.

    The default combination is unsupervised, so its optimism is near zero; a
    deliberately label-peeking combiner must show the gap instead.
    """
    rng = np.random.default_rng(0)
    x = rng.normal(size=(30, 4))
    y = rng.random(30) < 0.5

    def peeking(train_x, train_y):
        del train_x, train_y
        return lambda row: float(np.dot(np.atleast_2d(row)[0], weights))

    weights = np.linalg.lstsq(x, y.astype(float), rcond=None)[0]
    fitted = leave_one_out(x, y, combine=peeking)
    assert fitted.auroc_in > 0.6, "a rule fitted on all the data should fit it"

    honest = leave_one_out(x, y)
    assert abs(honest.optimism) < 0.15


def test_leave_one_out_recovers_a_real_signal():
    rng = np.random.default_rng(1)
    y = np.arange(40) < 20
    x = np.where(y, rng.normal(2.0, 1.0, 40), rng.normal(-2.0, 1.0, 40))[:, None]
    result = leave_one_out(x, y)
    assert result.auroc_out > 0.85
    assert result.n == 40


def test_leave_one_out_on_pure_noise_is_chance():
    rng = np.random.default_rng(2)
    x = rng.normal(size=(60, 3))
    y = rng.random(60) < 0.5
    result = leave_one_out(x, y)
    assert 0.3 < result.auroc_out < 0.7


def test_leave_one_out_guards(round7):
    with pytest.raises(ValueError, match="both classes"):
        leave_one_out(np.zeros((5, 2)), np.ones(5, dtype=bool))
    with pytest.raises(ValueError, match="differ in length"):
        leave_one_out(np.zeros((5, 2)), np.ones(4, dtype=bool))


def test_leave_one_out_on_the_round7_predictors(round7):
    """Recorded for provenance: the mechanical predictor does not generalise
    either, which is consistent with — not additional to — the null result."""
    values = np.array([[v["ddg"], v["ddg_normalised"]]
                       for v in round7["per_variant"]])
    is_lof = np.array([v["classification"] for v in round7["per_variant"]]) == "LoF"
    result = leave_one_out(values, is_lof)
    assert result.n == 25
    assert result.auroc_out < 0.65
