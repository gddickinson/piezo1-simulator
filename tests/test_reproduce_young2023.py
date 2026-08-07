"""Reproducing Young et al. 2023 end to end — one agreement and one disagreement.

Their rate constants, this project's solver and time-constant extraction,
checked against **two other papers**. That is what makes it an integration test
rather than a restatement: nothing here compares the model to the numbers it was
built from.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.parameters import PARAMETERS
from piezo1.physics.kinetics import GatingModel


@pytest.fixture(scope="module")
def model():
    return GatingModel()


def test_the_rates_are_youngs_as_published(model):
    """The inputs, before anything is concluded from them."""
    assert model.k1_0 == pytest.approx(5.1)
    assert model.km1_0 == pytest.approx(5.0)
    assert model.k2 == pytest.approx(8.0)
    assert model.km2 == pytest.approx(0.4)
    assert model.k3_0 == pytest.approx(34.6)
    assert model.k4 == pytest.approx(4.0)
    assert "Young" in model.provenance and "2023" in model.provenance


def test_microscopic_reversibility_holds(model):
    """A four-state cycle constrains one rate; if it drifts the model is not theirs."""
    for tension in (0.0, 1.0, 3.0, 8.0):
        assert abs(model.detailed_balance_residual(tension)) < 1e-9


# ------------------------------------------------------- what agrees

def test_half_activation_matches_lewis_2015(model):
    """Young's rates, our solver, Lewis's measurement — three separate things.

    2.711 mN/m against a measured 2.7 +/- 0.1. This is the strong result of the
    reproduction and the reason the disagreement below is informative rather
    than a sign the pipeline is broken.
    """
    measured = PARAMETERS.value("kinetics.t50_measured")
    t50 = model.half_activation()
    assert t50 == pytest.approx(measured, abs=0.15)
    assert abs(t50 - measured) / measured < 0.02


# ---------------------------------------------------- what does not

def test_inactivation_is_eight_times_slower_than_bae_2013(model):
    """The finding, pinned so it cannot drift into or out of agreement.

    Young's published rates give tau ~73 ms at 5 mN/m; Bae measured 8.6 +/- 0.4
    ms. The two papers used different preparations, and k2 -- the O->I1 rate --
    is what carries the timescale: at 8/s it sets ~125 ms before the rest of the
    system pulls it to 73.

    This is why the project calibrates mutants by FOLD CHANGE against the
    wild-type tau rather than by absolute tau across preparations. That policy
    was written before this measurement; the measurement is what justifies it.
    """
    measured = PARAMETERS.value("kinetics.wt_tau_ms")
    tau_ms = model.step(5.0, duration=1.0, n_points=8000).inactivation_tau() * 1e3

    assert tau_ms == pytest.approx(73.3, rel=0.05)
    ratio = tau_ms / measured
    assert 7.0 < ratio < 10.0, (
        f"the disagreement has moved: tau {tau_ms:.1f} ms is {ratio:.1f}x the "
        f"measured {measured} ms")


def test_reaching_the_measured_tau_needs_a_large_k2_increase(model):
    """How far the model has to be pushed to match the other paper."""
    measured = PARAMETERS.value("kinetics.wt_tau_ms")
    scale = model.calibrate_k2_for_tau(measured / 1000.0, hi=40.0)
    assert 10.0 < scale < 16.0

    tuned = model.with_modification("calibrated", k2=scale)
    achieved = tuned.step(5.0, duration=1.0,
                          n_points=8000).inactivation_tau() * 1e3
    assert achieved == pytest.approx(measured, rel=0.1)


def test_the_target_is_refused_when_out_of_reach_rather_than_clipped(model):
    """The default search cannot reach 8.6 ms, and says so.

    Returning the search bound would look like an answer. The reachable range
    appears in the message, which is what a caller needs to widen it.
    """
    with pytest.raises(ValueError, match="outside the reachable range"):
        model.calibrate_k2_for_tau(PARAMETERS.value("kinetics.wt_tau_ms") / 1000.0)


# ------------------------------------------------------- the API I misread

def test_with_modification_scales_rather_than_sets():
    """It takes fold changes, not absolute rates, and the difference is silent.

    Passing ``k2=8.0`` to a model whose k2 is already 8 gives **64**, not 8.
    Reading it as a setter produced a tau of 13 ms where the model really gives
    73 — a plausible number, and wrong by the factor being tested. Pinned
    because the two readings differ by exactly the quantity this round measures.
    """
    base = GatingModel()
    assert base.k2 == pytest.approx(8.0)
    doubled = base.with_modification("probe", k2=2.0)
    assert doubled.k2 == pytest.approx(16.0)
    unchanged = base.with_modification("probe", k2=1.0)
    assert unchanged.k2 == pytest.approx(base.k2)

    with pytest.raises(ValueError, match="unknown rate"):
        base.with_modification("bad", not_a_rate=2.0)


# ------------------------------------------------------------ the shape

def test_the_response_activates_then_inactivates(model):
    """The qualitative claim of the four-state scheme, at every tension tested."""
    for tension in (1.0, 3.0, 5.0, 8.0):
        result = model.step(tension, duration=1.0, n_points=4000)
        open_p = result.occupancy[:, 1]
        peak = int(np.argmax(open_p))
        assert 0 < peak < len(open_p) - 1, f"no interior peak at {tension} mN/m"
        assert open_p[-1] < open_p[peak], "must decay after the peak"
        assert np.all(open_p >= -1e-12)
        # Occupancies are a distribution at every instant.
        assert np.allclose(result.occupancy.sum(axis=1), 1.0, atol=1e-8)


def test_peak_open_probability_rises_and_saturates(model):
    tensions = np.array([0.5, 1.0, 2.0, 3.0, 5.0, 8.0, 12.0])
    peak = np.array([float(np.ravel(model.peak_open_probability(t))[0])
                     for t in tensions])
    assert np.all(np.diff(peak) > -1e-9), "peak Po must not fall with tension"
    # Saturating: the last step adds far less than the first.
    assert (peak[-1] - peak[-2]) < 0.2 * (peak[2] - peak[1])
