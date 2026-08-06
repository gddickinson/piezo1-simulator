"""Helfrich membrane mechanics and dome energetics."""

import numpy as np
import pytest

from piezo1.physics.dome import (PUBLISHED_AREA_ESTIMATES, DomeGeometrySummary,
                                 DomeModel, half_activation_tension,
                                 open_probability)
from piezo1.physics.membrane import (KT_PER_NM2_IN_MNM, FootprintSolution,
                                     MembraneParameters, analytic_energy,
                                     analytic_profile, decay_length,
                                     footprint_energy, kt_per_nm2_to_mnm,
                                     mnm_to_kt_per_nm2, solve_footprint)

KAPPA = 20.0
LAMBDA = 14.0
GAMMA = KAPPA / LAMBDA ** 2


@pytest.fixture(scope="module")
def params():
    return MembraneParameters(kappa=KAPPA, tension=GAMMA)


# --------------------------------------------------------------------------
# Units and the decay length
# --------------------------------------------------------------------------

def test_tension_unit_conversion_matches_the_literature():
    """1 k_BT/nm² = 4.114 mN/m at 298 K. Getting this wrong is a 4x error."""
    assert KT_PER_NM2_IN_MNM == pytest.approx(4.114, abs=0.005)
    assert mnm_to_kt_per_nm2(kt_per_nm2_to_mnm(0.37)) == pytest.approx(0.37)


def test_published_lambda_kappa_gamma_triple(params):
    """λ = 14 nm with κ = 20 k_BT implies γ = 0.42 mN/m.

    The three numbers are quoted together in Haselwandter & MacKinnon 2018 and
    must be mutually consistent through λ = √(κ/γ).
    """
    assert params.decay_length == pytest.approx(LAMBDA, rel=1e-9)
    assert params.tension_mnm == pytest.approx(0.42, abs=0.01)
    assert decay_length(KAPPA, mnm_to_kt_per_nm2(0.42)) == pytest.approx(14.0, rel=1e-2)


def test_zero_tension_gives_an_infinite_footprint():
    assert decay_length(20.0, 0.0) == float("inf")
    with pytest.raises(ValueError, match="infinite footprint"):
        solve_footprint(10.0, 0.2, MembraneParameters(20.0, 0.0))


# --------------------------------------------------------------------------
# The exact solution
# --------------------------------------------------------------------------

def test_analytic_profile_satisfies_its_boundary_slope():
    r = np.linspace(10.0, 10.0 + 12 * LAMBDA, 4000)
    h = analytic_profile(r, 10.0, 0.35, LAMBDA)
    slope = np.gradient(h, r, edge_order=2)[0]
    assert slope == pytest.approx(-0.35, rel=2e-3)
    assert abs(h[-1]) < abs(h[0]) * 1e-4        # decays


def test_closed_form_energy_matches_direct_integration(params):
    """The boundary-term formula must equal the integrated functional.

    This is the check that caught the Bessel ratio being inverted: the closed
    form was written with K₁/K₀ instead of K₀/K₁, which is 2.5x too large at
    the ratio r₀/λ where PIEZO1 sits, and nothing else would have revealed it.
    """
    r0, slope = 10.0, 0.35
    r = np.linspace(r0, r0 + 30 * LAMBDA, 20000)
    ref = FootprintSolution(r=r, h=analytic_profile(r, r0, slope, LAMBDA),
                            params=params, r0=r0, slope=slope)
    assert footprint_energy(ref) == pytest.approx(
        analytic_energy(r0, slope, KAPPA, GAMMA), rel=1e-4)


def test_energy_scales_as_slope_squared():
    e1 = analytic_energy(10.0, 0.2, KAPPA, GAMMA)
    e2 = analytic_energy(10.0, 0.4, KAPPA, GAMMA)
    assert e2 / e1 == pytest.approx(4.0, rel=1e-9)


# --------------------------------------------------------------------------
# The numerical solver
# --------------------------------------------------------------------------

def test_solver_reproduces_the_exact_profile(params):
    sol = solve_footprint(10.0, 0.35, params, n=1600)
    exact = analytic_profile(sol.r, 10.0, 0.35, LAMBDA)
    assert np.abs(sol.h - exact).max() / np.abs(exact).max() < 1e-3


def test_solver_energy_matches_the_closed_form(params):
    sol = solve_footprint(10.0, 0.35, params, n=1600)
    assert sol.energy == pytest.approx(
        analytic_energy(10.0, 0.35, KAPPA, GAMMA), rel=1e-3)


def test_solver_is_second_order_convergent(params):
    """Error must fall ~4x per grid doubling.

    Building the biharmonic operator as L@L gave first order at best and an
    energy that converged to the wrong answer entirely; the coupled
    second-order formulation plus second-order boundary stencils gives this.
    """
    exact_e = analytic_energy(10.0, 0.35, KAPPA, GAMMA)
    errors = []
    for n in (400, 800, 1600):
        sol = solve_footprint(10.0, 0.35, params, n=n)
        errors.append(abs(sol.energy - exact_e) / exact_e)
    for coarse, fine in zip(errors, errors[1:]):
        assert coarse / fine > 3.0, f"convergence ratio {coarse / fine:.1f}"


def test_decay_length_is_recovered_from_the_profile(params):
    sol = solve_footprint(10.0, 0.35, params, n=1600)
    assert sol.fitted_decay_length() == pytest.approx(LAMBDA, rel=0.02)


def test_footprint_shrinks_as_tension_rises():
    """Higher tension pulls the membrane flat, so λ and stored area fall."""
    lam, area = [], []
    for t_mnm in (0.42, 2.0, 10.0):
        p = MembraneParameters.from_mnm(KAPPA, t_mnm)
        sol = solve_footprint(10.0, 0.3, p, n=1000)
        lam.append(sol.decay_length)
        area.append(sol.excess_area())
    assert lam[0] > lam[1] > lam[2]
    assert area[0] > area[1] > area[2]


def test_small_slope_validity_is_reported(params):
    """The linear theory must declare when it is being used out of range."""
    ok = solve_footprint(10.0, 0.2, params, n=800)
    assert ok.within_linear_regime
    assert "within the small-slope regime" in ok.validity_note()

    bad = solve_footprint(10.0, 2.0, params, n=800)
    assert not bad.within_linear_regime
    assert "EXCEEDS" in bad.validity_note()


def test_surface_revolution_shape(params):
    sol = solve_footprint(10.0, 0.2, params, n=200)
    x, y, z = sol.surface(n_angular=32)
    assert x.shape == y.shape == z.shape == (200, 32)
    assert np.allclose(np.sqrt(x[0] ** 2 + y[0] ** 2), sol.r0, atol=1e-6)


# --------------------------------------------------------------------------
# Dome energetics
# --------------------------------------------------------------------------

def test_cox_parameters_reproduce_the_measured_t50():
    """ΔG₀ = 9.7 k_BT and ΔA = 8 nm² must give back T₅₀ ≈ 5.1 mN/m.

    All three numbers are from Cox et al. 2016, so this is an internal
    consistency check on the two-state dome model: fed that paper's energetics
    it must reproduce that paper's half-activation tension.
    """
    m = DomeModel()
    assert m.half_activation_mnm == pytest.approx(5.1, abs=0.3)


def test_open_probability_is_a_sigmoid_through_one_half():
    m = DomeModel()
    assert float(m.open_probability(0.0)) < 0.01
    assert float(m.open_probability(m.half_activation_mnm)) == pytest.approx(0.5, abs=1e-6)
    assert float(m.open_probability(20.0)) > 0.99
    t = np.linspace(0, 15, 60)
    p = m.open_probability(t)
    assert np.all(np.diff(p) >= -1e-12)          # monotonic


def test_half_activation_is_delta_g_over_delta_area():
    assert half_activation_tension(8.0, 9.7) == pytest.approx(9.7 / 8.0)
    assert half_activation_tension(0.0, 9.7) == float("inf")


def test_structural_areas_predict_too_low_a_threshold():
    """The documented discrepancy, made quantitative.

    Structural ΔA values, taken at face value as gating areas, imply
    half-activation tensions far below anything measured. They are measuring a
    different quantity, and the model should show that rather than hide it.
    """
    rows = {r["key"]: r for r in DomeModel().compare_area_estimates()}
    assert rows["cox2016_functional"]["t50_mnm"] == pytest.approx(5.0, abs=0.3)
    for key in ("dixit2025_nanodome", "guo2017_projected", "yang2022_inplane"):
        assert rows[key]["t50_mnm"] < 2.0
        assert rows[key]["kind"] == "structural"
    assert len(PUBLISHED_AREA_ESTIMATES) >= 4


def test_geometry_summary_converts_angstrom_to_nanometres():
    class FakeDome:
        radius_of_curvature = 97.0     # Angstrom
        dome_area = 49300.0            # Angstrom^2
        projected_area = 23700.0
        dome_depth = 49.0
        footprint_radius = 87.0

    g = DomeGeometrySummary.from_measurement(FakeDome())
    assert g.radius_of_curvature == pytest.approx(9.7)
    assert g.dome_area == pytest.approx(493.0)
    assert g.excess_area == pytest.approx(256.0)


def test_model_reports_both_area_families():
    text = DomeModel().report()
    assert "functional" in text and "structural" in text
    assert "Measured T50" in text


def test_footprint_requires_geometry():
    with pytest.raises(ValueError, match="no geometry"):
        DomeModel().footprint()
