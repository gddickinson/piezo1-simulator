"""1-D PNP over the measured pore.

The load-bearing checks are the ones that do not reuse the solver's own
machinery: an analytic cylinder, the independent series-resistance formula, and
the measured conductance against the published 25-30 pS. The numerical faults
that had to be fixed to get here — a singular system, a diverging Gummel loop —
are pinned too, because each of them produced a confident number before it
produced an error.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.parameters import PARAMETERS
from piezo1.physics.permeation import (F_FARADAY, R_GAS, IonSpecies,
                                       access_resistance, blocking_mechanisms,
                                       debye_length, default_species,
                                       series_conductance, solve_pnp)


class _Profile:
    """The minimum a PoreProfile needs to expose for the solver."""

    def __init__(self, z, radius):
        self.z = np.asarray(z, dtype=float)
        self.radius = np.asarray(radius, dtype=float)


def _cylinder(radius_A: float, length_A: float, step: float = 1.0):
    z = np.arange(0.0, length_A + step, step)
    return _Profile(z, np.full_like(z, radius_A))


# ------------------------------------------------------- against known answers

def test_uniform_cylinder_matches_the_textbook_resistance():
    """R = L / (sigma A) for a straight pipe, with no solver involved.

    Derived from the bulk conductivity and the accessible area directly, so a
    sign error or a factor of two anywhere in the discretisation shows up here
    rather than being absorbed into a plausible-looking conductance.
    """
    radius_A, length_A = 6.0, 60.0
    profile = _cylinder(radius_A, length_A)
    species = default_species()
    temperature = PARAMETERS.value("permeation.temperature")

    result = solve_pnp(profile)
    assert result.is_conducting

    # Expected pore resistance, computed independently.
    per_length = 0.0
    for s in species:
        usable = (radius_A - s.radius) * 1e-10
        per_length += ((s.valence ** 2 * F_FARADAY ** 2 * s.diffusivity
                        * s.concentration * 1000.0) / (R_GAS * temperature)
                       ) * np.pi * usable ** 2
    expected_ohm = (length_A * 1e-10) / per_length
    assert result.pore_ohm == pytest.approx(expected_ohm, rel=0.02)


def test_pnp_agrees_with_the_independent_series_formula(open_profile):
    """Two routes to the same number, sharing no code.

    `series_conductance` integrates 1/(sigma A) in closed form; `solve_pnp`
    discretises drift-diffusion and iterates a potential. At low voltage they
    must coincide, and on 11ZC they do to 1.5%.
    """
    result = solve_pnp(open_profile)
    series = series_conductance(np.asarray(open_profile.z) * 1e-10,
                                np.asarray(open_profile.radius) * 1e-10)
    assert result.converged
    assert result.conductance == pytest.approx(series["conductance"], rel=0.05)


def test_conductance_scales_the_way_geometry_says_it_should():
    """Twice as long, half the conductance; twice as wide, more than four times.

    The second is not a typo: widening also lowers the access resistance, so the
    gain exceeds the area ratio. Checking the direction and rough magnitude of
    both catches an area/radius confusion that a single-geometry test cannot.
    """
    short = solve_pnp(_cylinder(6.0, 40.0)).conductance
    long = solve_pnp(_cylinder(6.0, 80.0)).conductance
    assert long < short
    # Access resistance is shared, so the ratio is close to but above 0.5.
    assert 0.5 <= long / short <= 0.75

    narrow = solve_pnp(_cylinder(4.0, 60.0)).conductance
    wide = solve_pnp(_cylinder(8.0, 60.0)).conductance
    assert wide / narrow > 4.0


def test_access_resistance_is_hall_s_formula():
    assert access_resistance(1e-9, 1.0) == pytest.approx(1.0 / (4.0 * 1e-9))
    assert access_resistance(0.0, 1.0) == np.inf
    # It is not negligible: for a short wide pore it is a real share of the total.
    result = solve_pnp(_cylinder(10.0, 20.0))
    assert result.access_ohm > 0.05 * result.pore_ohm


# ------------------------------------------------------- the numerical faults

def test_row_scaling_fixes_the_conditioning_without_changing_the_answer():
    """The interior rows are ~1e-18 while a Dirichlet row is 1.

    That is eighteen orders of magnitude between two kinds of row, and on the
    real pore LAPACK reported the matrix as singular because of it. The rows are
    homogeneous, so scaling each by its own largest coefficient is free — this
    checks that it collapses the condition number and still returns the exact
    solution of a problem whose answer is known.
    """
    from piezo1.physics.permeation import _solve_with_dirichlet

    n = 40
    matrix = np.zeros((n, n))
    for k in range(1, n - 1):
        matrix[k, k - 1], matrix[k, k], matrix[k, k + 1] = 1e-18, -2e-18, 1e-18
    matrix[0, 0] = matrix[-1, -1] = 1.0
    rhs = np.zeros(n)
    rhs[0], rhs[-1] = 0.0, 1.0

    scale = np.max(np.abs(matrix), axis=1)
    assert np.linalg.cond(matrix) / np.linalg.cond(matrix / scale[:, None]) > 1e12

    # A Laplace problem with these boundaries has the linear solution.
    solution = _solve_with_dirichlet(matrix, rhs)
    assert np.allclose(solution, np.linspace(0.0, 1.0, n), atol=1e-9)


def test_bernoulli_is_continuous_and_does_not_overflow():
    from piezo1.physics.permeation import _bernoulli

    x = np.array([-1e6, -800.0, -1e-12, 0.0, 1e-12, 800.0, 1e6])
    b = _bernoulli(x)
    assert np.all(np.isfinite(b))
    assert b[3] == pytest.approx(1.0)
    # B(x) - B(-x) = -x, the identity the Scharfetter-Gummel flux relies on.
    y = np.array([-5.0, -0.3, 0.3, 5.0])
    assert np.allclose(_bernoulli(y) - _bernoulli(-y), -y)
    assert b[-1] == pytest.approx(0.0, abs=1e-300)


def test_the_solver_converges_on_a_real_pore(open_profile):
    result = solve_pnp(open_profile)
    assert result.converged, (
        f"did not converge in {result.iterations} iterations — the Gummel loop "
        f"has regressed")
    assert np.all(np.isfinite(result.potential))
    # The potential stays between the two bath values; excursions of volts were
    # the signature of the diverging Poisson coupling.
    assert result.potential.min() >= -1e-6
    assert result.potential.max() <= result.voltage + 1e-6


def test_double_layers_overlap_in_this_pore(open_profile):
    """Why the electroneutral limit is used, stated as a number.

    If this ever becomes false the full Poisson coupling is worth revisiting;
    while it is true, a converged PNP here would still not be trustworthy.
    """
    species = default_species()
    length = debye_length(species, PARAMETERS.value("permeation.temperature"),
                          PARAMETERS.value("permeation.permittivity_pore"))
    assert 4e-10 < length < 9e-10
    assert length > float(np.min(open_profile.radius)) * 1e-10
    assert solve_pnp(open_profile).meta["double_layers_overlap"]


# -------------------------------------------------------------- the gating

def test_every_blocking_mechanism_is_reported_not_just_the_first():
    """The distinction the round exists to make.

    8YEZ is shut sterically *and* by dewetting; 7WLU only sterically. Returning
    on the first match makes the two look identical.
    """
    class _Wetting:
        available = True

        def __init__(self, gate, steric, score):
            self.hydrophobic_gate, self.sterically_occluded = gate, steric
            self.score = score

    species = default_species()
    narrow = np.full(5, 0.9e-10)
    wide = np.full(5, 4.0e-10)

    both = blocking_mechanisms(_Wetting(True, True, 0.99), narrow, species)
    assert len(both) == 2
    assert any("steric" in m for m in both) and any("hydrophobic" in m for m in both)

    steric_only = blocking_mechanisms(_Wetting(False, True, 0.11), narrow, species)
    assert len(steric_only) == 1 and "steric" in steric_only[0]

    # Wide enough to pass an ion, but it would dewet: gated by chemistry alone.
    gate_only = blocking_mechanisms(_Wetting(True, False, 0.82), wide, species)
    assert len(gate_only) == 1 and "hydrophobic" in gate_only[0]

    assert blocking_mechanisms(_Wetting(False, False, 0.0), wide, species) == []


def test_a_dewetted_pore_carries_no_current_even_though_it_is_wide():
    """The gate is chemistry, and the continuum equations cannot see it.

    Nothing in drift-diffusion knows about a hydration shell, so a wide dewetted
    pore would otherwise be reported as freely conducting.
    """
    class _Wetting:
        available, hydrophobic_gate, sterically_occluded, score = True, True, False, 0.82

    profile = _cylinder(6.0, 60.0)
    assert solve_pnp(profile).is_conducting
    blocked = solve_pnp(profile, _Wetting())
    assert not blocked.is_conducting
    assert blocked.current == 0.0
    assert "hydrophobic" in blocked.blocked_by


# --------------------------------------------------- against the measurement

def test_open_structure_against_the_published_conductance(open_profile):
    """Validation: measured 25-30 pS, model 41 pS — reported, not tuned.

    The model is high by about half. What decides whether that matters is the
    sensitivity below: the two confinement parameters are unmeasured, and moving
    them across plausible ranges straddles the measurement. So this number is
    the right order and is not an independent prediction, and the test pins both
    halves of that so neither can drift quietly.
    """
    published = PARAMETERS.value("permeation.published_conductance")
    measured = solve_pnp(open_profile).conductance_pS

    assert 30.0 < measured < 60.0, "the recorded 41 pS has moved"
    assert 1.2 < measured / published < 2.0


def test_the_answer_is_dominated_by_two_unmeasured_parameters(open_profile):
    """Which is why the agreement above cannot be claimed as a prediction."""
    values = []
    try:
        for scale in (0.25, 1.0):
            PARAMETERS.set_value("permeation.diffusion_scale", scale)
            for radius in (1.0, 2.0):
                PARAMETERS.set_value("permeation.radius_cation", radius)
                values.append(solve_pnp(open_profile).conductance_pS)
    finally:
        PARAMETERS.reset("permeation.diffusion_scale")
        PARAMETERS.reset("permeation.radius_cation")

    assert max(values) / min(values) > 4.0, "the sweep should span a wide range"
    assert min(values) < 30.0 < max(values), (
        "the published range should sit inside the span, which is what makes "
        "agreement tuning rather than prediction")


def test_calcium_permeates_but_carries_little_of_the_current(open_profile):
    """PIEZO1 is weakly selective, so calcium is a minority carrier at 2 mM.

    The species are named for the ions their registered constants describe —
    K+, Cl-, Ca2+ — rather than "cation" and "anion", which hid that the model
    was never told which ion it was solving for.
    """
    result = solve_pnp(open_profile, species=default_species(calcium=0.002))
    assert "Ca2+" in result.fluxes
    assert abs(result.fluxes["Ca2+"]) > 0
    carried = 2 * abs(result.fluxes["Ca2+"])
    total = sum(abs(s.valence * result.fluxes[s.name])
                for s in default_species(calcium=0.002))
    assert carried / total < 0.05


def test_parameters_are_registered_and_the_unmeasured_ones_say_so():
    for key in ("permeation.temperature", "permeation.bath_concentration",
                "permeation.test_voltage", "permeation.diffusion_cation",
                "permeation.diffusion_anion", "permeation.diffusion_calcium",
                "permeation.radius_cation", "permeation.diffusion_scale",
                "permeation.permittivity_pore",
                "permeation.published_conductance"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, f"{key} is not registered"
        assert parameter.citation and parameter.description, key

    for key in ("permeation.diffusion_scale", "permeation.permittivity_pore"):
        assert PARAMETERS.get(key).citation == "unverified", key
    assert PARAMETERS.get("permeation.published_conductance").citation == "coste2010piezo"
