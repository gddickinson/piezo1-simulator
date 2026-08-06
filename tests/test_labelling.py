"""HaloTag labelling — an import, tested as one.

The kinetics were not derived here. They come from ``halotag_binding_sim``, and
the criterion for this round was that they be reproduced **exactly**: any
divergence means the import is wrong, not that this project found something. So
the load-bearing test is a numerical comparison against the original functions,
and the rest check that the vendored copy behaves like the mathematics it claims
to implement — because a fresh clone will not have the source project on its
path and would otherwise be running unchecked code.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.labelling import (LabellingConditions, compare_with_source,
                                       cumulative_exposure, detectable_fraction,
                                       fully_labelled_fraction,
                                       occupancy_distribution,
                                       population_summary, predicted_brightness,
                                       simulate_population,
                                       site_labelled_fraction, time_to_fraction)
from piezo1.parameters import PARAMETERS


# ------------------------------------------------------- the import criterion

def test_reproduces_the_source_project_exactly():
    """Round 32's stated pass condition.

    Not "close": every quantity must agree to the last bit, and the stochastic
    populations must be identical channel for channel. That last part is why
    the sampler reproduces the source's two `rng.random((n, 3))` draws in the
    same order — a different order gives a statistically identical population
    and a numerically different one, and would hide a real divergence behind
    sampling noise.
    """
    result = compare_with_source()
    if not result["available"]:
        pytest.skip(result["reason"])

    assert result["max_abs_diff_p_site"] == 0.0
    assert result["max_abs_diff_fully_labelled"] == 0.0
    assert result["max_abs_diff_occupancy"] == 0.0
    assert result["max_abs_diff_hist"] == 0.0
    assert result["dye_counts_identical"], (
        "the analytical curves agree but the sampled populations do not; the "
        "draw order has diverged from the source")


# ----------------------------------------------------------- the mathematics

def test_exposure_is_the_integral_of_intracellular_concentration():
    """E(t) must be the integral of partition*[L]*(1 - exp(-k_perm t)).

    Checked by quadrature rather than by re-deriving the closed form, so the
    test cannot agree with the implementation by sharing its algebra.
    """
    conditions = LabellingConditions()
    t_end = 600.0
    fine = np.linspace(0.0, t_end, 200001)
    concentration = (conditions.partition * conditions.concentration
                     * (1.0 - np.exp(-conditions.k_perm * fine)))
    assert float(cumulative_exposure(t_end, conditions)) == pytest.approx(
        float(np.trapezoid(concentration, fine)), rel=1e-9)


def test_exposure_becomes_linear_once_access_saturates():
    """After many time constants the lag is a fixed offset, not a slope."""
    conditions = LabellingConditions()
    tau = 1.0 / conditions.k_perm
    steady = conditions.partition * conditions.concentration
    late = np.array([50 * tau, 60 * tau])
    exposure = cumulative_exposure(late, conditions)
    slope = np.diff(exposure)[0] / np.diff(late)[0]
    assert slope == pytest.approx(steady, rel=1e-9)
    # The offset is exactly one time constant's worth of ligand.
    assert float(exposure[0]) == pytest.approx(steady * (late[0] - tau), rel=1e-9)


def test_labelling_is_monotonic_and_bounded_by_the_active_fraction():
    for active in (1.0, 0.9, 0.5):
        conditions = LabellingConditions(active_fraction=active)
        t = np.linspace(0.0, 4 * 3600.0, 2000)
        p = site_labelled_fraction(t, conditions)
        assert np.all(np.diff(p) >= -1e-15), "labelling cannot go backwards"
        assert p[0] == pytest.approx(0.0, abs=1e-15)
        assert np.all(p <= active + 1e-12)
        assert float(p[-1]) == pytest.approx(active, rel=1e-6)


def test_occupancy_is_a_distribution_and_matches_its_moments():
    p = np.linspace(0.0, 1.0, 51)
    dist = occupancy_distribution(p, 3)
    assert dist.shape == (51, 4)
    assert np.allclose(dist.sum(axis=1), 1.0)
    assert np.all(dist >= -1e-15)

    ks = np.arange(4)
    assert np.allclose((dist * ks).sum(axis=1), 3 * p)          # mean = np
    variance = (dist * (ks - 3 * p[:, None]) ** 2).sum(axis=1)
    assert np.allclose(variance, 3 * p * (1 - p))               # var = np(1-p)

    # The two headline summaries are the tails of the same distribution.
    assert np.allclose(dist[:, 3], fully_labelled_fraction(p, 3))
    assert np.allclose(1.0 - dist[:, 0], detectable_fraction(p, 3))


def test_the_p_cubed_amplification():
    """The crux: every site must bind, so a per-site shortfall is cubed."""
    for p, expected in ((0.9, 0.729), (0.8, 0.512), (0.7, 0.343)):
        assert float(fully_labelled_fraction(p, 3)) == pytest.approx(expected,
                                                                     abs=1e-12)
    # A 10% per-site shortfall becomes a 27% channel-level shortfall.
    assert 1.0 - 0.9 ** 3 == pytest.approx(0.271, abs=1e-3)


# ------------------------------------------------------------- the sampler

def test_sampled_population_converges_on_the_analytical_law():
    conditions = LabellingConditions(concentration=1e-9)   # a partial regime
    times = [60.0, 300.0, 900.0]
    population = simulate_population(times, conditions, n_channels=200000, seed=3)
    analytical = population_summary(np.array(times), conditions)

    assert np.allclose(population["fully_labeled"], analytical.fully_labelled,
                       atol=4e-3)
    assert np.allclose(population["hist"], analytical.occupancy, atol=4e-3)
    assert np.allclose(population["mean_dyes"], analytical.mean_dyes, atol=1e-2)


def test_a_channel_never_loses_a_dye():
    """The labelling is covalent, so counts must be non-decreasing in time.

    This is what forces the uniform draw to happen once and be reused across
    snapshots; resampling per snapshot gives the right marginals and lets an
    individual channel go from three dyes back to two.
    """
    conditions = LabellingConditions(concentration=1e-9)
    times = np.linspace(0.0, 1800.0, 40)
    counts = simulate_population(times, conditions, n_channels=500, seed=7)["dye_counts"]
    assert np.all(np.diff(counts, axis=0) >= 0)


def test_unreactive_tags_cap_the_population_too():
    conditions = LabellingConditions(active_fraction=0.8)
    population = simulate_population([6 * 3600.0], conditions,
                                     n_channels=100000, seed=1)
    assert population["fully_labeled"][0] == pytest.approx(0.8 ** 3, abs=5e-3)


def test_brightness_keeps_the_dye_levels_separable():
    """The histogram is only informative if 1, 2 and 3 dyes stay distinct."""
    counts = np.repeat([1, 2, 3], 20000)
    brightness = predicted_brightness(counts, per_dye_intensity=1.0, seed=2)
    means = [brightness[counts == k].mean() for k in (1, 2, 3)]
    assert means == pytest.approx([1.0, 2.0, 3.0], abs=0.05)
    # Spread scales with the signal, so the levels do not smear into each other.
    spreads = [brightness[counts == k].std() for k in (1, 2, 3)]
    cv = PARAMETERS.value("labelling.brightness_noise_cv")
    assert spreads == pytest.approx([cv * 1, cv * 2, cv * 3], rel=0.1)


# ---------------------------------------------------------- what it predicts

def test_the_standard_protocol_saturates_long_before_it_ends():
    """Measured result: 200 nM for 30 min is complete in about a minute.

    Which matters for what the model can be used to claim — see the mixture
    test below.
    """
    conditions = LabellingConditions()
    assert time_to_fraction(0.99, conditions) < 120.0
    protocol = population_summary(
        np.array([PARAMETERS.value("labelling.incubation_time")]), conditions)
    assert float(protocol.fully_labelled[0]) == pytest.approx(1.0, abs=1e-6)


def test_a_dye_mixture_needs_unreactive_tags_not_a_shorter_incubation():
    """The distinction the model actually draws.

    Two different things get called "sub-saturation labelling". At standard
    concentrations the kinetic route is not available: labelling is complete in
    under a minute, so no plausible incubation leaves a mixture. A population of
    chemically unreactive tags does leave one, at every time, because the
    ceiling is active_fraction cubed. So an observed 1:2:3 brightness mixture
    under a saturating protocol argues for unreactive tags rather than for a
    short incubation.
    """
    saturating = LabellingConditions()
    at_protocol = occupancy_distribution(
        float(site_labelled_fraction(1800.0, saturating)), 3)
    assert at_protocol[3] == pytest.approx(1.0, abs=1e-6)
    assert at_protocol[:3].sum() < 1e-6, "no kinetic mixture at 200 nM"

    # Even a tenth of the protocol time is still saturating.
    assert float(occupancy_distribution(
        float(site_labelled_fraction(180.0, saturating)), 3)[3]) > 0.999

    partial = LabellingConditions(active_fraction=0.9)
    ceiling = occupancy_distribution(
        float(site_labelled_fraction(6 * 3600.0, partial)), 3)
    assert ceiling[3] == pytest.approx(0.729, abs=1e-3)
    assert ceiling[2] == pytest.approx(0.243, abs=1e-3)
    assert ceiling[:3].sum() > 0.25, "unreactive tags do leave a mixture"


def test_an_unreachable_target_is_reported_as_unreachable():
    """A ceiling below the target is a real answer, not a failure to converge."""
    conditions = LabellingConditions(active_fraction=0.9)   # ceiling 0.729
    assert not np.isfinite(time_to_fraction(0.8, conditions))
    assert np.isfinite(time_to_fraction(0.7, conditions))


# ------------------------------------------------------------- the structure

def test_labelling_maps_onto_the_modelled_tag_positions(structure_by_id):
    """What the structure adds: the sites are places, not slots."""
    from piezo1.config import STRUCTURE_DIR
    from piezo1.analysis.labelling import label_sites
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.fusion import HALOTAG_PDB, build_fusion, load_halotag

    st = structure_by_id("8YEZ")
    if st is None or not (STRUCTURE_DIR / f"{HALOTAG_PDB}.cif").exists():
        pytest.skip("8YEZ or 6U32 not downloaded")

    framed = apply_frame(st, canonical_transform(st))
    model = build_fusion(framed, load_halotag())
    result = label_sites(model, t=1800.0)

    assert result["tag_centres"].shape == (3, 3)
    assert result["anchor_residues"] == [2521, 2521, 2521]
    assert len(result["occupied"]) == 3
    assert 0 <= result["n_dyes"] <= 3
    assert sum(result["hist"]) == pytest.approx(1.0)
    # Saturating protocol, so every site is taken.
    assert result["n_dyes"] == 3


def test_parameters_are_registered_with_sources():
    for key in ("labelling.k_on", "labelling.n_sites", "labelling.k_perm_live",
                "labelling.partition_live", "labelling.active_fraction",
                "labelling.concentration", "labelling.incubation_time",
                "labelling.brightness_noise_cv"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, f"{key} is not registered"
        assert parameter.citation and parameter.description, key

    # The rate is the one measured quantity in the set and must cite a paper.
    assert PARAMETERS.get("labelling.k_on").citation == "los2008halotag"
    # The transport rate is not measured, and must not pretend to be.
    assert PARAMETERS.get("labelling.k_perm_live").citation == "unverified"


def test_conditions_follow_a_parameter_override():
    """Resolved at construction, so the dialog takes effect on the next call."""
    default = LabellingConditions().k_on
    PARAMETERS.set_value("labelling.k_on", 1.0e6)
    try:
        assert LabellingConditions().k_on == pytest.approx(1.0e6)
    finally:
        PARAMETERS.reset("labelling.k_on")
    assert LabellingConditions().k_on == pytest.approx(default)
