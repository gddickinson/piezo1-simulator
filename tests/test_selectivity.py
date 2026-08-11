"""Selectivity from fixed charge, and the sign error that finding it exposed.

The order here is the point. Nothing about PIEZO1 is measured until three
things with known answers have been checked: that a cation drifts *downhill*,
that the GHK inversion returns its two analytic limits, and that an uncharged
pore returns the ratio its ions' mobilities and sizes demand rather than one.
The first of those is not decoration — it failed. The Scharfetter-Gummel drift
term had its two Bernoulli factors on the wrong nodes, so cations climbed the
field, and no test in the suite could see it because every current the project
had ever computed was between identical baths with its sign discarded.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.parameters import PARAMETERS
from piezo1.physics._pnp_kernels import F_FARADAY, R_GAS, _nernst_planck
from piezo1.physics.permeation import solve_pnp, sodium_species
from piezo1.physics.selectivity import (REVERSAL_BRACKET_V,
                                        ghk_permeability_ratio,
                                        measure_selectivity,
                                        reversal_potential)


class _Profile:
    def __init__(self, z, radius):
        self.z = np.asarray(z, dtype=float)
        self.radius = np.asarray(radius, dtype=float)


def _cylinder(radius_A: float, length_A: float, step: float = 1.0):
    z = np.arange(0.0, length_A + step, step)
    return _Profile(z, np.full_like(z, radius_A))


def _thermal() -> float:
    return R_GAS * PARAMETERS.value("permeation.temperature") / F_FARADAY


# ------------------------------------------- calibration: the drift direction

def test_a_cation_drifts_towards_the_lower_potential():
    """The check the sign error failed, on a case with no discretisation in it.

    At a field weak enough that Scharfetter-Gummel is indistinguishable from a
    centred difference, the flux of a uniform concentration is exactly
    ``-D A z c phi' / phi_T`` — magnitude *and* sign, with no Bernoulli
    subtlety anywhere. Before Round 81 this returned the right magnitude with
    the wrong sign, which is invisible to every conductance test in the suite.
    """
    thermal = _thermal()
    z = np.arange(0.0, 61.0, 1.0) * 1e-10
    area = np.full_like(z, np.pi * (5e-10) ** 2)
    slope = 1e-3 / (60e-10)                      # 1 mV rising over the pore
    potential = np.linspace(0.0, 1e-3, len(z))

    _, flux = _nernst_planck(z, area, potential, +1, 1e-9, thermal, 150.0, 150.0)
    expected = -1e-9 * area[0] * 150.0 * slope / thermal
    assert flux < 0.0, "a cation must move away from the high-potential end"
    assert flux == pytest.approx(expected, rel=1e-6)

    _, anion = _nernst_planck(z, area, potential, -1, 1e-9, thermal, 150.0, 150.0)
    assert anion == pytest.approx(-expected, rel=1e-6)


# --------------------------------------------- calibration: the GHK inversion

def test_ghk_inversion_returns_its_two_analytic_limits():
    """Perfect cation selectivity gives the Nernst potential; none gives zero."""
    thermal = _thermal()
    inside, outside = 0.15, 0.03

    nernst = thermal * np.log(outside / inside)
    assert nernst < 0.0
    assert ghk_permeability_ratio(nernst, inside, outside) == np.inf
    assert ghk_permeability_ratio(0.0, inside, outside) == pytest.approx(1.0)

    # An entirely anion-selective pore reverses the other way, by as much.
    assert ghk_permeability_ratio(-nernst, inside, outside) == pytest.approx(
        0.0, abs=1e-12)


def test_ghk_inversion_round_trips_a_known_ratio():
    """Put a ratio in through the forward equation, get it back out."""
    thermal = _thermal()
    inside, outside = 0.15, 0.03
    for ratio in (0.05, 0.5, 1.0, 3.0, 20.0):
        forward = thermal * np.log((ratio * outside + inside)
                                   / (ratio * inside + outside))
        assert ghk_permeability_ratio(forward, inside, outside) == pytest.approx(
            ratio, rel=1e-9)


# ----------------------------------- calibration: what an uncharged pore does

def test_an_uncharged_pore_returns_its_mobility_ratio_not_one():
    """A neutral pore is not unselective, and the answer is calculable.

    With no fixed charge the dilution potential is the liquid-junction
    potential, set by how fast each ion moves and how much of the cross-section
    it can reach. In a wide pore Cl- is the faster ion and the size difference
    hardly matters, so the pore must come out slightly *anion*-favouring —
    ``P_Cl/P_Na`` a little above one. A result of exactly one would mean the
    diffusion potential was missing altogether, which is what the old ohmic
    closure would have given.
    """
    result = measure_selectivity(_cylinder(12.0, 60.0))
    mobility = (PARAMETERS.value("permeation.diffusion_anion")
                / PARAMETERS.value("permeation.diffusion_sodium"))
    assert mobility > 1.0
    assert 1.02 < result.p_anion_over_cation < mobility
    assert result.reversal_mV > 0.0


def test_a_narrow_uncharged_pore_is_already_cation_selective_from_size():
    """The baseline a charged answer has to beat.

    Chloride's crystal radius is nearly twice sodium's, so at a 3 A bottleneck
    the anion loses more cross-section than it gains in mobility and the pore
    prefers cations before any charge is added. Reporting a cation-selective
    charged pore without this number beside it would credit the charge with
    something geometry had already done.
    """
    published = PARAMETERS.value("permeation.published_pcl_pna")
    neutral = measure_selectivity(_cylinder(3.0, 40.0))
    assert neutral.cation_selective
    # ...but nowhere near as selective as the channel measures. Size alone
    # accounts for part of PIEZO1's preference for cations and not for most
    # of it, which is what leaves the charge something to do.
    assert neutral.p_anion_over_cation > 2 * published


# ------------------------------------------- calibration: charge, both signs

def test_fixed_charge_moves_selectivity_the_way_its_sign_says():
    """Negative charge prefers cations, positive prefers anions, monotonically.

    The instrument has to be able to say "no": if it reported cation
    selectivity for both signs it would be measuring something other than the
    charge. It is also the test that fails, loudly and in the right direction,
    when the drift term's sign is wrong.
    """
    profile = _cylinder(6.0, 60.0)
    n = len(profile.z)
    ratios = [measure_selectivity(profile,
                                  fixed_charge=np.full(n, x)).p_anion_over_cation
              for x in (-500.0, -50.0, 0.0, 50.0, 500.0)]
    assert ratios == sorted(ratios), f"not monotone in the charge: {ratios}"
    assert ratios[0] < 0.05, "half-molar negative charge should nearly exclude Cl-"
    assert ratios[-1] > 10.0, "and half-molar positive charge should nearly exclude Na+"


def test_a_strongly_charged_pore_saturates_at_the_nernst_limit():
    """A perfect ion exchanger cannot reverse past the counterion's Nernst potential."""
    profile = _cylinder(6.0, 60.0)
    n = len(profile.z)
    thermal = _thermal()
    high = PARAMETERS.value("permeation.dilution_high")
    low = PARAMETERS.value("permeation.dilution_low")
    limit = thermal * np.log(low / high)

    result = measure_selectivity(profile, fixed_charge=np.full(n, -5000.0))
    assert limit < result.reversal_V < 0.95 * limit
    assert abs(result.reversal_V) < REVERSAL_BRACKET_V


# ------------------------------------------------- the equations, not the fit

def test_zero_fixed_charge_reproduces_the_neutral_pore_exactly():
    """Round 81's validation clause, both halves.

    An explicit array of zeros must give the untouched arithmetic bit for bit,
    so wiring the charge in cannot have moved the recorded conductance. And a
    charge nine orders of magnitude below the bath must reproduce it through
    the *other* code path — which is what shows the branch is a statement about
    the equations rather than a way of skipping them.
    """
    for profile in (_cylinder(6.0, 60.0), _cylinder(4.0, 25.0)):
        n = len(profile.z)
        neutral = solve_pnp(profile)
        zeros = solve_pnp(profile, fixed_charge=np.zeros(n))
        assert zeros.conductance == neutral.conductance
        assert np.array_equal(zeros.potential, neutral.potential)
        assert not zeros.meta["fixed_charge"]

        tiny = solve_pnp(profile, fixed_charge=np.full(n, 1e-9))
        assert tiny.meta["fixed_charge"], "the charged path should have been used"
        assert tiny.conductance == pytest.approx(neutral.conductance, rel=1e-9)


def test_a_charged_pore_at_equilibrium_carries_no_current_and_partitions_exactly():
    """Donnan equilibrium is the one case with a closed-form answer.

    At zero applied voltage between identical baths, the counterion
    concentration must be the Boltzmann factor of the local Donnan potential
    and every species flux must be zero. Both are exact, so this is a hard
    check rather than a tolerance.
    """
    profile = _cylinder(6.0, 60.0)
    n = len(profile.z)
    charge = np.zeros(n)
    charge[28:33] = -13600.0

    result = solve_pnp(profile, voltage=0.0, fixed_charge=charge)
    assert result.converged
    for flux in result.fluxes.values():
        assert abs(flux) < 1e-28
    # Electroneutrality holds where the charge is, which is what makes this a
    # space-charge model rather than a decoration on a neutral one.
    assert result.meta["electroneutrality_residual"] < 1e-8
    cation = result.concentrations["K+"]
    assert cation[30] == pytest.approx(13600.0 + 150.0 ** 2 / 13600.0, rel=1e-6)


def test_the_packing_ceiling_flags_rather_than_clips():
    """A density no solution could reach is reported, and nothing is changed."""
    profile = _cylinder(6.0, 60.0)
    n = len(profile.z)
    ceiling = PARAMETERS.value("pore_charge.max_concentration")

    mild = solve_pnp(profile, fixed_charge=np.full(n, -100.0))
    assert not mild.meta["exceeds_packing_limit"]
    assert mild.meta["peak_in_pore_M"] < ceiling

    absurd = solve_pnp(profile, fixed_charge=np.full(n, -50000.0))
    assert absurd.meta["exceeds_packing_limit"]
    # Flagged, not clipped: the concentration really is what the model says.
    assert absurd.meta["peak_in_pore_M"] > ceiling


def test_a_blocked_pore_reports_no_selectivity_rather_than_a_number():
    class _Wetting:
        available, hydrophobic_gate, sterically_occluded, score = True, True, False, 0.9

    with pytest.raises(ValueError, match="no current to reverse"):
        reversal_potential(_cylinder(6.0, 60.0), wetting=_Wetting())


def test_getting_the_orientation_wrong_inverts_the_answer_plausibly():
    """Why the cytosolic end is measured rather than assumed.

    Describing the *same* experiment from either end must give the same
    membrane potential — the concentrated bath is on the cytosolic side either
    way. But telling the model the wrong end while leaving the baths where they
    are turns a cation-selective pore into an anion-selective one, and the
    number that comes out is entirely reasonable-looking. That is the failure
    a sign convention produces if nobody checks it.
    """
    profile = _cylinder(6.0, 60.0)
    n = len(profile.z)
    charge = np.full(n, -500.0)
    high = PARAMETERS.value("permeation.dilution_high")
    low = PARAMETERS.value("permeation.dilution_low")

    first = measure_selectivity(profile, fixed_charge=charge, cytosolic_index=0)
    last = measure_selectivity(profile, fixed_charge=charge, cytosolic_index=-1)
    assert first.cation_selective and last.cation_selective
    assert first.reversal_V == pytest.approx(last.reversal_V, rel=1e-3)

    # Same physical arrangement, wrong label on the ends.
    fixed_baths = sodium_species(high, low)
    right_way = measure_selectivity(profile, fixed_charge=charge,
                                    cytosolic_index=0, species=fixed_baths)
    wrong_way = measure_selectivity(profile, fixed_charge=charge,
                                    cytosolic_index=-1, species=fixed_baths)
    assert right_way.cation_selective
    assert not wrong_way.cation_selective
    assert right_way.reversal_V == pytest.approx(-wrong_way.reversal_V, rel=1e-3)


# ------------------------------------------------------ against the structure

def test_the_charged_pore_is_cation_selective_on_the_real_structure(open_profile):
    """Round 81's result: the direction is right, the number is not.

    Both routes make the model prefer cations, which is the direction the
    measurement has and which the neutral pore only weakly shows. Neither
    reproduces the published 0.14 — they bracket it about sevenfold apart — and
    the curated route only gets there at an in-pore concentration a continuum
    model has no business quoting. Pinned so none of the three can drift.
    """
    from piezo1.core import Structure
    from piezo1.config import STRUCTURE_DIR
    from piezo1.physics.pore_charge import pore_charge
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    st = Structure.from_file(STRUCTURE_DIR / "11ZC.cif")
    st = apply_frame(st, canonical_transform(st))
    blocks, _ = protomer_blocks(st)
    axis = detect_c3_axis(blocks)
    published = PARAMETERS.value("permeation.published_pcl_pna")

    neutral = measure_selectivity(open_profile, cytosolic_index=0)
    assert neutral.p_anion_over_cation > 5 * published

    ratios = {}
    for mode in ("curated", "lining"):
        charge = pore_charge(st, open_profile, axis, mode=mode, species="mouse")
        assert charge.n_groups > 0
        result = measure_selectivity(open_profile, fixed_charge=charge.density,
                                     cytosolic_index=0)
        assert result.cation_selective, mode
        assert result.p_anion_over_cation < neutral.p_anion_over_cation
        ratios[mode] = result.p_anion_over_cation

    assert ratios["curated"] < published < ratios["lining"]
    assert ratios["lining"] / ratios["curated"] > 5.0, (
        "the two routes are supposed to disagree; if they have converged the "
        "reported bracket is no longer what it says")


def test_the_curated_route_is_outside_the_models_validity(open_profile):
    """The result that stops the agreement being claimed.

    Six carboxylates in a 3 A lumen demand a counterion concentration above ten
    molar. That is not a solution, so the ratio the curated route reports is
    not a prediction — and this is checked rather than left to prose.
    """
    from piezo1.core import Structure
    from piezo1.config import STRUCTURE_DIR
    from piezo1.physics.pore_charge import pore_charge
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    st = Structure.from_file(STRUCTURE_DIR / "11ZC.cif")
    st = apply_frame(st, canonical_transform(st))
    blocks, _ = protomer_blocks(st)
    charge = pore_charge(st, open_profile, detect_c3_axis(blocks),
                         mode="curated", species="mouse")
    result = solve_pnp(open_profile, fixed_charge=charge.density)
    assert result.converged
    assert result.meta["exceeds_packing_limit"]
    assert result.meta["peak_in_pore_M"] > 10.0
    assert result.meta["electroneutrality_residual"] < 1e-6


def test_the_selectivity_parameters_are_registered_with_their_source():
    for key in ("permeation.published_pcl_pna", "permeation.dilution_high",
                "permeation.dilution_low", "permeation.diffusion_sodium",
                "permeation.radius_sodium"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, f"{key} is not registered"
        assert parameter.description
    for key in ("permeation.published_pcl_pna", "permeation.dilution_high",
                "permeation.dilution_low"):
        assert PARAMETERS.get(key).citation == "coste2015", key
