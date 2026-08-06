"""Tension-dependent Markov gating kinetics."""

import numpy as np
import pytest

from piezo1.physics.kinetics import (MUTANT_PRESETS, STATE_NAMES, GatingModel)


@pytest.fixture(scope="module")
def model():
    return GatingModel()


def test_published_rate_values(model):
    """Rates must match the Young et al. 2023 specification exactly."""
    r0 = model.rates(0.0)
    assert r0["k1"] == pytest.approx(5.1)
    assert r0["k2"] == pytest.approx(8.0)
    assert r0["km2"] == pytest.approx(0.4)
    assert r0["k3"] == pytest.approx(34.6)
    assert r0["k4"] == pytest.approx(4.0)
    assert r0["km4"] == pytest.approx(0.6)
    # k-1 = 5 * exp(sigma_50 / b) = 28.8 s^-1, stated in the paper.
    assert r0["km1"] == pytest.approx(28.8, abs=0.1)


def test_tension_dependence_directions(model):
    """Activation must accelerate with tension and recovery must slow."""
    lo, hi = model.rates(0.0), model.rates(4.0)
    assert hi["k1"] > lo["k1"]         # C -> O
    assert hi["k3"] < lo["k3"]         # I1 -> C
    assert hi["km1"] == pytest.approx(lo["km1"])   # O -> C is tension-free


def test_microscopic_reversibility_holds(model):
    """C→I₁ is fixed by detailed balance, not free.

    If this ever fails the model is generating energy from nothing.
    """
    for sigma in (0.0, 1.4, 3.0, 8.0):
        assert model.detailed_balance_residual(sigma) < 1e-12


def test_generator_matrix_is_valid(model):
    q = model.rate_matrix(2.0)
    assert q.shape == (4, 4)
    assert np.allclose(q.sum(axis=1), 0.0)          # rows sum to zero
    off = q - np.diag(np.diag(q))
    assert (off >= 0).all()                          # off-diagonals are rates
    assert (np.diag(q) <= 0).all()


def test_steady_state_is_a_probability_vector(model):
    for sigma in (0.0, 2.0, 10.0):
        p = model.steady_state(sigma)
        assert p.shape == (4,)
        assert p.sum() == pytest.approx(1.0)
        assert (p >= 0).all()


def test_steady_state_is_a_true_null_vector(model):
    p = model.steady_state(3.0)
    assert np.abs(model.rate_matrix(3.0).T @ p).max() < 1e-9


def test_half_activation_matches_measurement(model):
    """Emergent half-activation should land on the measured T50.

    σ₅₀ = 1.4 mN/m is a *parameter inside a rate*, not the model's output. The
    quantity to compare with experiment is the peak open probability at which
    the channel is half-maximally activated, and the measured cell-attached
    value is 2.7 ± 0.1 mN/m (Lewis & Grandl 2015).
    """
    t50 = model.half_activation(peak=True)
    assert 2.2 < t50 < 3.3, f"half-activation {t50:.2f} mN/m"


def test_step_response_activates_then_inactivates(model):
    res = model.step(5.0, duration=1.0, n_points=2000)
    assert res.occupancy.shape == (2000, 4)
    assert np.allclose(res.occupancy.sum(axis=1), 1.0)
    # Rises to a peak, then decays: the signature of an inactivating channel.
    assert res.peak_open() > res.occupancy[-1, 1]
    assert res.peak_time() > 0
    tau = res.inactivation_tau()
    assert tau is not None and 0.001 < tau < 1.0


def test_current_sign_and_scale(model):
    res = model.step(5.0, duration=0.2, n_channels=100)
    # Inward current at -80 mV must be negative.
    assert res.current.min() < 0
    unitary = model.conductance_pS * model.holding_mV * 1e-3
    assert res.current.min() >= 100 * unitary - 1e-9


def test_mutants_move_in_the_measured_direction(model):
    """Gain-of-function mutants must inactivate more slowly than wild type.

    This is the check that caught a real bug: calibrating to an absolute time
    constant measured in a different preparation made R2456H come out *faster*
    than wild type. Fold changes transfer between preparations; absolute time
    constants do not.
    """
    wt = model.step(5.0, duration=3.0, n_points=4000).inactivation_tau()
    for name, preset in MUTANT_PRESETS.items():
        mut = model.mutant(name)
        tau = mut.step(5.0, duration=3.0, n_points=4000).inactivation_tau()
        ratio = tau / wt
        assert ratio == pytest.approx(preset["tau_ratio"], rel=0.05), name
        if preset["classification"] == "GoF":
            assert ratio > 1.0, f"{name} should inactivate more slowly"


def test_mutant_carries_its_provenance(model):
    mut = model.mutant("R2456H")
    assert "phenomenological" in mut.provenance
    assert "R2456H" in mut.label


def test_single_channel_matches_the_analytic_steady_state(model):
    """A long Gillespie run must reproduce the equilibrium occupancies."""
    sigma = 3.0
    dwell_total = np.zeros(4)
    for seed in range(8):
        t, s = model.simulate_single_channel(sigma, duration=200.0, seed=seed)
        dwell = np.diff(np.append(t, 200.0))
        for i in range(4):
            dwell_total[i] += dwell[s == i].sum()
    empirical = dwell_total / dwell_total.sum()
    analytic = model.steady_state(sigma)
    assert np.abs(empirical - analytic).max() < 0.08, (
        f"empirical {empirical} vs analytic {analytic}")


def test_mean_open_time_is_reciprocal_of_exit_rate(model):
    q = model.rate_matrix(2.0)
    assert model.mean_open_time(2.0) == pytest.approx(-1.0 / q[1, 1])


def test_alternative_i2_topology_is_available(model):
    alt = GatingModel(i2_downstream_of_i1=True)
    q = alt.rate_matrix(2.0)
    assert q[2, 3] > 0 and q[1, 3] == 0      # I1 -> I2, not O -> I2
    p = alt.steady_state(2.0)
    assert p.sum() == pytest.approx(1.0)


def test_ramp_protocol_runs(model):
    res = model.ramp(8.0, duration=0.5, n_points=500)
    assert res.tension[-1] == pytest.approx(8.0)
    assert np.allclose(res.occupancy.sum(axis=1), 1.0)
    assert len(STATE_NAMES) == 4
