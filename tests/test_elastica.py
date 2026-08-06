"""Nonlinear axisymmetric membrane mechanics.

The shape equations in :mod:`piezo1.physics.elastica` were derived by hand.
A BVP solver will happily converge on a wrong derivation, so the checks that
matter here are the ones that do not reuse the derivation: agreement with the
linear theory where both are valid, conservation of the axial force, and
evaluation of the *exact* Helfrich functional in a different gauge.
"""

import numpy as np
import pytest

from piezo1.physics.dome import DomeGeometrySummary
from piezo1.physics.elastica import (axial_force, compare_with_linear,
                                     solve_elastica)
from piezo1.physics.membrane import MembraneParameters

#: The measured 7WLT dome: projected radius and the contact slope its
#: spherical-cap geometry implies. 63 degrees.
PIEZO1_R0 = 8.691
PIEZO1_SLOPE = 1.992


@pytest.fixture(scope="module")
def params():
    return MembraneParameters()


@pytest.fixture(scope="module")
def piezo1(params):
    return solve_elastica(PIEZO1_R0, PIEZO1_SLOPE, params)


# --------------------------------------------------------------------------
# Agreement with the linear theory where the linear theory is valid
# --------------------------------------------------------------------------

def test_converges_to_linear_at_small_slope(params):
    """The whole point of the round: the two must agree as slope -> 0."""
    c = compare_with_linear(PIEZO1_R0, 0.02, params)
    assert c.linear_error < 1e-3
    assert c.energy_ratio == pytest.approx(1.0, abs=1e-3)


def test_linear_error_is_second_order_in_slope(params):
    """Error/slope² must approach a constant, not merely shrink.

    The Monge expansion drops terms of order |∇h|², so this exponent is a
    statement about which terms were dropped — a shrinking error alone would
    be consistent with a solver that is simply converging to the wrong thing.
    """
    ratios = []
    for slope in (0.05, 0.1, 0.2, 0.4):
        c = compare_with_linear(PIEZO1_R0, slope, params)
        ratios.append(c.linear_error / slope ** 2)
    assert np.allclose(ratios, ratios[-1], rtol=0.10), ratios
    assert 0.5 < ratios[-1] < 1.0


def test_excess_area_agrees_at_small_slope(params):
    c = compare_with_linear(PIEZO1_R0, 0.05, params)
    assert c.area_ratio == pytest.approx(1.0, abs=0.01)


# --------------------------------------------------------------------------
# Checks that do not reuse the derivation
# --------------------------------------------------------------------------

def test_axial_force_is_conserved(piezo1):
    """H is a constant of the motion and must be zero for a free inclusion.

    This is a free accuracy diagnostic: nothing in the integration enforces H
    away from the boundary, so drift measures the error directly.
    """
    assert piezo1.force_residual < 1e-8


def test_exact_functional_in_another_gauge_agrees(piezo1, params):
    """Re-evaluate the energy in the Monge gauge with no expansion at all.

    Exact axisymmetric forms: dA = √(1+h'²)·2πr·dr and
    c₁+c₂ = h''/(1+h'²)^{3/2} + h'/(r√(1+h'²)). If the hand-derived
    arc-length equations minimise a different functional from the one intended,
    these two numbers separate.
    """
    r, h = piezo1.height_profile(n=1500)
    dr = r[1] - r[0]
    hp = np.gradient(h, dr, edge_order=2)
    hpp = np.gradient(hp, dr, edge_order=2)
    g = 1.0 + hp ** 2
    curvature = hpp / g ** 1.5 + hp / (r * np.sqrt(g))
    da = np.sqrt(g) * 2 * np.pi * r

    bending = np.trapezoid(0.5 * params.kappa * curvature ** 2 * da, r)
    excess = np.trapezoid(da, r) - np.pi * (r[-1] ** 2 - r[0] ** 2)
    energy = bending + params.tension * excess

    assert energy == pytest.approx(piezo1.energy, rel=5e-3)
    assert excess == pytest.approx(piezo1.excess_area, rel=5e-3)


def test_no_nearby_shape_has_lower_energy(piezo1, params):
    """A solution of the shape equations must be a minimum, not a saddle.

    Perturbing the profile by shape-preserving bumps that respect both boundary
    conditions may only raise the energy. This catches a sign error in the
    Euler-Lagrange derivation, which would leave the solver converging happily
    onto a stationary point that is not the minimum.
    """
    r, h = piezo1.height_profile(n=800)
    dr = r[1] - r[0]

    def energy_of(hh):
        hp = np.gradient(hh, dr, edge_order=2)
        hpp = np.gradient(hp, dr, edge_order=2)
        g = 1.0 + hp ** 2
        curvature = hpp / g ** 1.5 + hp / (r * np.sqrt(g))
        da = np.sqrt(g) * 2 * np.pi * r
        return (np.trapezoid(0.5 * params.kappa * curvature ** 2 * da, r)
                + params.tension * (np.trapezoid(da, r)
                                    - np.pi * (r[-1] ** 2 - r[0] ** 2)))

    base = energy_of(h)
    span = r[-1] - r[0]
    rng = np.random.default_rng(0)
    for _ in range(12):
        centre = r[0] + span * rng.uniform(0.05, 0.5)
        width = span * rng.uniform(0.01, 0.1)
        amp = rng.uniform(-1.0, 1.0) * 0.05
        # Vanishes at both ends, and its derivative vanishes at r0, so the
        # contact slope and both boundary conditions survive the perturbation.
        bump = amp * ((r - r[0]) / span) ** 2 * np.exp(-((r - centre) / width) ** 2)
        assert energy_of(h + bump) >= base - 1e-9


# --------------------------------------------------------------------------
# Numerical robustness
# --------------------------------------------------------------------------

def test_result_is_independent_of_domain_truncation(params):
    """A footprint energy that depends on where you stopped integrating is a
    truncation artefact rather than a physical quantity."""
    values = [solve_elastica(PIEZO1_R0, PIEZO1_SLOPE, params,
                             arc_span=mult * params.decay_length).energy
              for mult in (8, 20, 40)]
    # 1e-5 is the solver's own adaptive-mesh noise; a genuine truncation
    # artefact would grow with the domain, not scatter at the sixth digit.
    assert np.allclose(values, values[0], rtol=1e-5), values


def test_continuation_and_direct_solve_agree(params):
    direct = solve_elastica(PIEZO1_R0, PIEZO1_SLOPE, params, continuation=False)
    walked = solve_elastica(PIEZO1_R0, PIEZO1_SLOPE, params, continuation=True)
    assert direct.energy == pytest.approx(walked.energy, rel=1e-5)
    assert walked.meta["continuation_steps"] > 1


def test_grid_refinement_is_stable(params):
    coarse = solve_elastica(PIEZO1_R0, PIEZO1_SLOPE, params, n=200)
    fine = solve_elastica(PIEZO1_R0, PIEZO1_SLOPE, params, n=1200)
    assert coarse.energy == pytest.approx(fine.energy, rel=1e-5)


# --------------------------------------------------------------------------
# Geometry and guards
# --------------------------------------------------------------------------

def test_geometric_identities(piezo1):
    assert piezo1.converged
    assert piezo1.area > piezo1.projected_area          # deformed area is larger
    assert piezo1.excess_area > 0
    assert piezo1.s[-1] >= piezo1.r[-1] - piezo1.r0     # arc length >= chord
    assert piezo1.contact_angle_deg == pytest.approx(63.3, abs=0.5)
    assert piezo1.z[0] > 0 and piezo1.z[-1] == pytest.approx(0.0, abs=1e-6)
    assert np.all(np.diff(piezo1.r) > 0)                # no overhang


def test_profile_decays_monotonically(piezo1):
    assert np.all(np.diff(piezo1.z) <= 1e-9)
    assert abs(piezo1.z[-1]) < 1e-6 * piezo1.z[0]


def test_negative_slope_is_rejected(params):
    with pytest.raises(ValueError, match="magnitude"):
        solve_elastica(PIEZO1_R0, -1.0, params)


def test_zero_tension_is_rejected():
    with pytest.raises(ValueError, match="infinite footprint"):
        solve_elastica(10.0, 1.0, MembraneParameters(kappa=20.0, tension=0.0))


def test_axial_force_helper_is_zero_on_the_solution(piezo1, params):
    y = np.vstack([piezo1.r, piezo1.z, piezo1.psi, piezo1.curvature, piezo1.eta])
    assert np.abs(axial_force(y, params.kappa, params.tension)).max() < 1e-8


# --------------------------------------------------------------------------
# The PIEZO1 result
# --------------------------------------------------------------------------

def test_linear_theory_overestimates_piezo1_footprint(params):
    """Pin the measured correction so a solver change cannot move it silently.

    At the 7WLT dome's 63 degree contact slope the linearised Helfrich energy
    is 3.65x the true one and the excess area 3.48x. This is why the linear
    number could not support the conclusion Round 3 drew from it.
    """
    c = compare_with_linear(PIEZO1_R0, PIEZO1_SLOPE, params)
    assert c.linear_energy == pytest.approx(92.2, rel=0.02)
    assert c.nonlinear_energy == pytest.approx(25.3, rel=0.02)
    assert c.nonlinear_excess_area == pytest.approx(178.9, rel=0.02)
    assert 1.0 / c.energy_ratio == pytest.approx(3.65, rel=0.03)


def test_footprint_holds_less_area_than_the_dome(params):
    """The corrected comparison, and it reverses the Round 3 ordering.

    Round 3 compared the dome's *exact* excess area against the footprint's
    *linearised* one and concluded the footprint held 2.4x more. Measured
    consistently, the footprint holds 0.70x the dome's 255.9 nm².
    """
    c = compare_with_linear(PIEZO1_R0, PIEZO1_SLOPE, params)
    dome_excess = 255.87
    assert c.nonlinear_excess_area / dome_excess == pytest.approx(0.70, abs=0.03)
    assert c.linear_excess_area / dome_excess == pytest.approx(2.43, abs=0.05)


def test_correction_survives_the_published_parameter_range(params):
    """Not an artefact of one choice of κ and γ."""
    factors = []
    for kappa in (20.0, 25.0):
        for tension_mnm in (0.42, 1.0, 3.0):
            c = compare_with_linear(
                PIEZO1_R0, PIEZO1_SLOPE,
                MembraneParameters.from_mnm(kappa, tension_mnm))
            factors.append(1.0 / c.energy_ratio)
    assert min(factors) > 3.4 and max(factors) < 3.7, factors


def test_dome_model_comparison_on_measured_coordinates(curved_structure):
    """End to end from real coordinates, not the pinned constants above.

    Everything else in this file starts from the recorded 7WLT numbers. This
    one measures the dome from the deposited structure and lets the model
    derive its own contact slope, so a change in the geometry pipeline shows
    up here rather than silently invalidating the constants.
    """
    from piezo1.physics.dome import DomeModel
    from piezo1.structure.geometry import measure_dome

    from conftest import protomer_blocks
    from test_geometry import _tm_surface

    blocks, _ = protomer_blocks(curved_structure)
    dome = measure_dome(blocks, _tm_surface(curved_structure, "mouse"))
    model = DomeModel(geometry=DomeGeometrySummary.from_measurement(dome))

    assert model.contact_slope() == pytest.approx(2.0, abs=0.35)
    report = model.compare_footprint_theories()
    assert report["contact_angle_deg"] > 55.0
    assert report["linear_overestimate_energy"] > 3.0
    assert report["nonlinear_excess_area_nm2"] < report["dome_excess_area_nm2"]
    assert report["axial_force_residual"] < 1e-8
