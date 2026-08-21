"""Does mechanics predict what evolution protected? — the instrument, first.

This module joins an external evolutionary result to this project's elastic
network, and it is the single most dangerous thing in the family subsystem,
because the answer is a correlation and correlations are easy to manufacture:

* both series are strongly **autocorrelated along the chain**, so a permutation
  null is far too easy to beat and would make almost anything significant;
* **burial** explains most conservation in most proteins and correlates with
  nearly every mechanical quantity, so a raw correlation says little;
* **eight features** are tested against one track, so one of them clearing a
  threshold is not a result.

Every test below drives the instrument against an answer known in advance
before any PIEZO1 number is read. The order matters: a planted signal must be
recovered, a true null must not be, and the burial control must be shown to be
capable of destroying a correlation rather than merely being applied to one.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.constraint_mechanics import (BURIAL_FEATURES,
                                                  MECHANICAL_FEATURES,
                                                  align_track_to_features,
                                                  circular_shift_null, couple,
                                                  partial_spearman)
from piezo1.analysis.features import build_feature_table
from piezo1.analysis.fluctuations import spearman
from piezo1.core.family import load_constraint


# ------------------------------------------------------- calibration: the null

def test_the_shift_null_is_centred_on_zero_for_an_unrelated_feature():
    """A random feature must sit inside the null it is compared against."""
    rng = np.random.default_rng(0)
    values = np.cumsum(rng.normal(size=500))       # autocorrelated, like the track
    feature = np.cumsum(rng.normal(size=500))      # and unrelated
    null = circular_shift_null(values, feature, n=200, seed=1)
    observed = spearman(values, feature)
    z = (observed - null.mean()) / null.std()
    assert abs(z) < 3.0, f"an unrelated pair reached z = {z:.1f}"


def test_a_permutation_null_would_have_been_far_too_easy_to_beat():
    """The measured reason the null is a circular shift and not a shuffle.

    Both series are autocorrelated. Destroying that with a permutation shrinks
    the null's spread, so the same observed correlation scores a much larger z
    — which is how a comparison manufactures significance from structure that
    was in both series to begin with.
    """
    rng = np.random.default_rng(2)
    values = np.cumsum(rng.normal(size=500))
    feature = np.cumsum(rng.normal(size=500))
    shift_sd = circular_shift_null(values, feature, n=300, seed=3).std()
    permuted = np.array([spearman(rng.permutation(values), feature)
                         for _ in range(300)])
    assert shift_sd > 3 * permuted.std(), (
        f"shift null sd {shift_sd:.3f} against permutation {permuted.std():.3f} "
        f"— the shift is supposed to be the wider, harder null")


def test_the_null_recovers_a_planted_correlation_at_the_real_autocorrelation():
    """It must be able to say yes as well as no — planted on the real track.

    The planting is done on the *actual* constraint values rather than on a
    synthetic series, because the null's difficulty is set by how autocorrelated
    the track is and a synthetic one would set it somewhere else. A feature that
    is the real track plus noise must clear the null comfortably.
    """
    track = load_constraint("PIEZO1")
    values = track.values[~np.isnan(track.values)]
    rng = np.random.default_rng(4)
    feature = values + rng.normal(scale=0.05, size=values.size)
    null = circular_shift_null(values, feature, n=200, seed=5)
    z = (spearman(values, feature) - null.mean()) / null.std()
    assert z > 5, f"a planted correlation only reached z = {z:.1f}"


def test_the_shift_null_is_conservative_on_a_random_walk_and_that_is_correct():
    """Pinned because it looks like a failure and is the point of the design.

    Two independent random walks can correlate at |rho| near 1 by chance, so a
    shift null built from one of them is very wide and even a near-perfect
    planted correlation reaches only z ~ 2. A permutation null would call the
    same pair overwhelmingly significant. The real constraint track is nowhere
    near random-walk autocorrelation — the test above clears z = 5 on it — so
    this bounds the instrument rather than describing it.
    """
    rng = np.random.default_rng(11)
    walk = np.cumsum(rng.normal(size=500))
    feature = walk + rng.normal(scale=0.2, size=500)
    null = circular_shift_null(walk, feature, n=300, seed=12)
    z = (spearman(walk, feature) - null.mean()) / null.std()
    assert spearman(walk, feature) > 0.95
    assert z < 5, "a random-walk null is supposed to be hard to beat"


# -------------------------------------------- calibration: the burial control

def test_the_partial_correlation_removes_a_correlation_that_is_only_burial():
    """The control must be capable of destroying a signal, or it asserts nothing.

    Constructed so the feature and the track are related *only* through a
    shared confounder: partialling the confounder out must collapse it.
    """
    rng = np.random.default_rng(6)
    burial = rng.normal(size=400)
    feature = burial + rng.normal(scale=0.1, size=400)
    values = burial + rng.normal(scale=0.1, size=400)
    raw = spearman(values, feature)
    partial = partial_spearman(values, feature, burial[:, None])
    assert raw > 0.8, "the confounded correlation should be strong"
    assert abs(partial) < 0.2, f"burial was not removed: {partial:.2f}"


def test_the_partial_correlation_keeps_a_correlation_that_is_not_burial():
    """And it must not destroy a real one, or every answer would be negative."""
    rng = np.random.default_rng(7)
    burial = rng.normal(size=400)
    feature = rng.normal(size=400)
    values = feature + 0.5 * burial + rng.normal(scale=0.1, size=400)
    partial = partial_spearman(values, feature, burial[:, None])
    assert partial > 0.7, f"a genuine correlation was destroyed: {partial:.2f}"


def test_ranking_before_regressing_removes_more_of_a_skewed_confounder():
    """Why this is a partial *Spearman* rather than a partial Pearson.

    The case that decides it is a confounder acting **monotonically but not
    linearly** on the two series — which is what burial does: constraint rises
    with burial and then saturates, and relative SASA is bounded at zero. Here
    the feature depends on burial as a cube and the values as a log; both are
    driven entirely by burial and nothing else, so a control that works must
    reduce the correlation to nothing.

    A linear regression on the raw values cannot: it can only subtract a
    straight line from a curve. Ranking first turns both monotone dependences
    into the same linear one, which is exactly the situation linear regression
    handles. Measured against the alternative rather than asserted.
    """
    rng = np.random.default_rng(8)
    burial = rng.uniform(0.02, 1.0, size=400)
    feature = burial ** 3 + rng.normal(scale=1e-3, size=400)
    values = np.log(burial) + rng.normal(scale=1e-3, size=400)

    design = np.column_stack([np.ones(400), burial])

    def raw_partial(x, y):
        rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
        ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
        return abs(float(np.corrcoef(rx, ry)[0, 1]))

    ranked = abs(partial_spearman(values, feature, burial[:, None]))
    unranked = raw_partial(values.copy(), feature.copy())
    assert unranked > 0.3, ("the linear control is supposed to fail here, or "
                            "the test has nothing to show")
    assert ranked < 0.1, f"ranked control left {ranked:.2f}"
    assert ranked < unranked


# ------------------------------------------------------- on the real structure

@pytest.fixture(scope="module")
def coupled(curved_structure):
    features = build_feature_table(curved_structure)
    return couple(curved_structure, features=features, structure_id="7WLT"), features


def test_the_track_lines_up_with_the_feature_table_through_the_alignment(
        curved_structure, coupled):
    _, features = coupled
    values, mask = align_track_to_features(features, load_constraint("PIEZO1"),
                                           "mouse")
    assert mask.sum() > 1000
    assert np.isfinite(values[mask]).all()


def test_burial_is_reported_beside_every_mechanical_feature(coupled):
    """Without it the mechanical numbers are unreadable: burial explains most
    conservation in most proteins."""
    result, _ = coupled
    assert set(result.burial_alone) == set(BURIAL_FEATURES)
    assert max(abs(v) for v in result.burial_alone.values()) > 0.3


def test_every_feature_carries_a_corrected_p_because_eight_are_tested(coupled):
    result, _ = coupled
    assert len(result.features) >= 6
    for feature in result.features:
        assert 0.0 < feature.p_empirical <= 1.0
        assert feature.q_value >= feature.p_empirical - 1e-12


def test_coupling_to_the_gate_survives_both_controls(coupled):
    """The joint result neither project could reach alone.

    Not asserted as a large effect — it is not. What is asserted is that a
    gate-coupling feature is among those surviving the shift null, the
    multiplicity correction and the burial control, and that the verdict says
    how many of the eight did.
    """
    result, _ = coupled
    survivors = [f.feature for f in result.features
                 if f.survives_correction and f.survives_burial]
    assert survivors, result.verdict
    assert any("gate" in name for name in survivors), survivors
    assert str(len(survivors)) in result.verdict


def test_mobility_predicts_the_opposite_sign_to_gate_coupling(coupled):
    """The two halves of the census's core-periphery picture, in mechanics:
    residues that move are free to change, residues coupled to the gate are
    not."""
    result, _ = coupled
    msf = result.by_name("msf")
    gate = result.by_name("prs_gate_response")
    assert msf.spearman < 0 < gate.spearman


def test_an_entry_in_another_numbering_is_refused_rather_than_correlated(
        structure_by_id):
    piezo2 = structure_by_id("6KG7")
    if piezo2 is None:
        pytest.skip("6KG7 is not downloaded")
    result = couple(piezo2, structure_id="6KG7")
    assert result.n_residues == 0 and not result.features
    assert "refused" in result.note


def test_the_verdict_says_no_when_nothing_survives(coupled):
    """A verdict that can only say yes is not a verdict. Driven by raising the
    burial-retention bar past what any feature achieves."""
    from piezo1.parameters import PARAMETERS as _P

    result, _ = coupled
    _P.set_value("family.burial_retention", 1.0)
    try:
        assert "no mechanical quantity predicts constraint beyond burial" in \
            result.verdict
    finally:
        _P.reset("family.burial_retention")
    assert "no mechanical quantity" not in result.verdict


def test_every_named_mechanical_feature_exists_in_the_feature_table(coupled):
    """A feature silently missing would shrink the family being corrected for."""
    _, features = coupled
    for name in MECHANICAL_FEATURES + BURIAL_FEATURES:
        assert name in features.columns, f"{name} is not a feature column"
