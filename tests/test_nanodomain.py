"""The calcium nanodomain, and the one prediction this project can actually make.

Unusually for this project the round's prediction **held**, so the tests are
shaped differently from the usual null-guarding: they check the Green's function
against analytic limits it must obey, then pin the prediction *and the sweep
that makes it non-trivial*. A claim that survives every parameter in its own
model is only interesting if the range swept was wide enough to have broken it.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.parameters import PARAMETERS
from piezo1.physics.nanodomain import (Nanodomain, calcium_at,
                                       distance_for_occupancy, saturation,
                                       screening_length, sweep)
from piezo1.physics.permeation import F_FARADAY

#: The unitary current Round 33 measured on the one open-like structure.
CURRENT_A = 2.46e-12


# --------------------------------------------------------- the Green's function

def test_profile_matches_the_closed_form():
    """Computed independently of the implementation, from the equation itself."""
    current = 1e-13
    distance = 5e-9
    diffusivity = PARAMETERS.value("nanodomain.d_calcium")
    length = screening_length()
    expected = (current / (4.0 * np.pi * F_FARADAY * diffusivity * distance)
                * np.exp(-distance / length) / 1000.0)
    expected += PARAMETERS.value("nanodomain.resting_calcium")
    assert float(calcium_at(distance, current)) == pytest.approx(expected, rel=1e-12)


def test_it_falls_as_one_over_r_where_buffering_is_weak():
    """λ is ~150 nm, so inside a few nanometres the exponential is ~1.

    Which is worth checking rather than assuming: it means the screening term
    plays almost no role at the distance the tag actually sits, and the answer
    is set by geometry rather than by the buffer.
    """
    resting = PARAMETERS.value("nanodomain.resting_calcium")
    near = float(calcium_at(2e-9, CURRENT_A)) - resting
    far = float(calcium_at(4e-9, CURRENT_A)) - resting
    assert far == pytest.approx(near / 2.0, rel=0.02)


def test_screening_length_depends_only_on_the_buffer_product():
    """k_on and [B] are not separately identifiable from a nanodomain."""
    a = screening_length(buffer_kon=1e8, buffer_concentration=1e-4)
    b = screening_length(buffer_kon=1e7, buffer_concentration=1e-3)
    assert a == pytest.approx(b, rel=1e-12)
    # And it shortens as the square root of buffering, not linearly.
    assert screening_length(buffer_concentration=4e-4) == pytest.approx(a / 2.0,
                                                                        rel=1e-12)


def test_the_profile_is_monotonic_and_settles_on_the_resting_level():
    r = np.geomspace(1e-9, 1e-3, 400)
    c = calcium_at(r, CURRENT_A)
    assert np.all(np.diff(c) <= 0), "concentration must never rise with distance"
    # Strictly falling while the nanodomain is still above the resting level.
    resting = PARAMETERS.value("nanodomain.resting_calcium")
    near = c[c > 1.01 * resting]
    assert len(near) > 50 and np.all(np.diff(near) < 0)
    assert float(c[-1]) == pytest.approx(resting, rel=1e-6)


def test_resting_calcium_sets_a_floor_under_the_sensor():
    """The sensor is already a third occupied with no channel open.

    100 nM resting against a 0.2 uM Kd gives 33%, so its dynamic range is
    33-100%, not 0-100%. A target below the floor is unreachable at any
    distance, and the solver says so rather than returning its search bound.
    """
    nano = Nanodomain(current_A=CURRENT_A, calcium_fraction=0.05,
                      distance_m=3.95e-9)
    assert nano.resting_occupancy == pytest.approx(1 / 3, abs=0.01)
    assert distance_for_occupancy(0.2, CURRENT_A * 0.05) == float("inf")


def test_saturation_is_the_binding_isotherm():
    kd = PARAMETERS.value("nanodomain.sensor_kd")
    assert float(saturation(kd)) == pytest.approx(0.5)
    assert float(saturation(0.0)) == 0.0
    assert float(saturation(1e6 * kd)) > 0.999


# ------------------------------------------------------------- the prediction

def test_the_sensor_is_saturated_at_the_modelled_tag_distance():
    """Round 35's claim, at the numbers Rounds 31 and 33 actually produced.

    The roadmap expected ~200 uM at 4-6 nm. The tag centroid moved to 3.95 nm in
    Round 31 and the concentration comes out at ~114 uM — half the expectation,
    same order, and far above the sensor's 0.2 uM Kd either way. The conclusion
    the number was for is unchanged.
    """
    nano = Nanodomain(current_A=CURRENT_A,
                      calcium_fraction=PARAMETERS.value(
                          "nanodomain.calcium_current_fraction"),
                      distance_m=3.95e-9)
    assert nano.concentration_M * 1e6 == pytest.approx(114.0, rel=0.05)
    assert nano.occupancy > 0.99
    assert nano.saturated


def test_it_stays_saturated_across_the_whole_accessible_envelope():
    """The tag samples a region, so a point estimate would overstate precision.

    Round 31's envelope runs 1.74-7.89 nm on 8YEZ. Both ends must give the same
    verdict or the claim depends on where in the envelope you look.
    """
    nano = Nanodomain(current_A=CURRENT_A, calcium_fraction=0.05,
                      distance_m=3.95e-9, envelope_m=(1.74e-9, 7.89e-9))
    far_c, near_c, far_occ, near_occ = nano.envelope_range()
    assert near_c > far_c
    assert far_occ > 0.99 and near_occ > 0.99
    assert 50e-6 < far_c < 300e-6 and 50e-6 < near_c < 300e-6


def test_it_survives_the_conductance_being_uncertain():
    """Round 33's conductance spans 16-94 pS across unmeasured parameters."""
    for scale in (16 / 41.0, 1.0, 94 / 41.0):
        nano = Nanodomain(current_A=CURRENT_A * scale, calcium_fraction=0.05,
                          distance_m=3.95e-9)
        assert nano.saturated, f"not saturated at {scale * 41:.0f} pS"
        assert nano.occupancy > 0.99


def test_the_sweep_is_wide_enough_to_have_broken_the_claim():
    """A prediction nothing could falsify is not a prediction.

    The sweep must span combinations that *do* desaturate the sensor, so that
    the overwhelming majority which do not is informative rather than an
    artefact of a timid range.
    """
    rows = sweep(CURRENT_A)
    occupancies = [r["occupancy"] for r in rows]
    assert len(rows) >= 60
    assert min(occupancies) < 0.9, (
        "no swept combination desaturates the sensor; widen the range or the "
        "robustness claim means nothing")
    assert sum(o > 0.9 for o in occupancies) / len(rows) > 0.9

    # The exceptions are extreme in every direction at once.
    for row in (r for r in rows if r["occupancy"] < 0.9):
        assert row["distance_nm"] >= 20.0
        assert row["calcium_fraction"] <= 0.005
        assert row["buffer_M"] >= 1e-3


def test_the_falsifiers_are_quantitative_and_far_from_reality():
    """What would have to be true for the claim to fail, as numbers."""
    nano = Nanodomain(current_A=CURRENT_A, calcium_fraction=0.05,
                      distance_m=3.95e-9)
    f = nano.falsifiers()

    # The tag would have to sit ~100x further away than it is modelled to.
    assert f["distance_for_half_occupancy_m"] > 50 * nano.distance_m
    # Or calcium would have to carry ~1000x less of the current.
    assert f["calcium_fraction_for_half_occupancy"] < 1e-3
    # Or free buffer would have to be far above any cytosolic value.
    assert f["buffer_concentration_needed_M"] > 0.01


def test_distance_for_occupancy_inverts_the_profile():
    current = CURRENT_A * 0.05
    for target in (0.5, 0.8, 0.95):
        r = distance_for_occupancy(target, current)
        assert np.isfinite(r)
        assert float(saturation(calcium_at(r, current))) == pytest.approx(
            target, rel=1e-3)

    with pytest.raises(ValueError):
        distance_for_occupancy(1.0, current)


# ------------------------------------------------------- against the structures

def test_it_uses_the_numbers_the_earlier_rounds_produced(structure_by_id):
    """End to end: pore -> current -> tag geometry -> nanodomain.

    The point of the round was that this needed no new measurements, only the
    two numbers Rounds 31 and 33 already produce.
    """
    from piezo1.analysis.hydration import load_grid, predict_wetting
    from piezo1.config import STRUCTURE_DIR
    from piezo1.physics.permeation import solve_pnp
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.fusion import HALOTAG_PDB, build_fusion, load_halotag
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    open_st = structure_by_id("11ZC")
    human = structure_by_id("8YEZ")
    if (open_st is None or human is None
            or not (STRUCTURE_DIR / f"{HALOTAG_PDB}.cif").exists()):
        pytest.skip("11ZC, 8YEZ or 6U32 not downloaded")

    framed = apply_frame(open_st, canonical_transform(open_st))
    blocks, _ = protomer_blocks(framed)
    profile = pore_profile(framed, detect_c3_axis(blocks), step=1.0)
    grid = load_grid()
    wetting = predict_wetting(framed, profile, grid) if grid.available else None
    current = abs(solve_pnp(profile, wetting).current)
    assert current > 0, "the open structure should conduct"

    tagged = apply_frame(human, canonical_transform(human))
    model = build_fusion(tagged, load_halotag())
    distance = float(model.pore_exit_distances()[0]) * 1e-9

    nano = Nanodomain(current_A=current, calcium_fraction=0.05,
                      distance_m=distance)
    assert nano.saturated
    assert 20e-6 < nano.concentration_M < 500e-6


def test_parameters_are_registered_and_the_assumed_one_says_so():
    for key in ("nanodomain.d_calcium", "nanodomain.buffer_kon",
                "nanodomain.buffer_concentration", "nanodomain.sensor_kd",
                "nanodomain.calcium_current_fraction",
                "nanodomain.resting_calcium"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, f"{key} is not registered"
        assert parameter.citation and parameter.description, key

    assert PARAMETERS.get(
        "nanodomain.calcium_current_fraction").citation == "unverified"
    assert PARAMETERS.get("nanodomain.sensor_kd").citation == "tsien1980bapta"


def test_a_closed_structure_does_not_report_its_own_current(structure_by_id):
    """The frame trap, pinned.

    The C3 axis has to come from the same frame as the coordinates it indexes.
    Detecting it on the unframed structure and applying it to the framed one
    measures the pore along a line that misses the pore — and reported the
    *closed* 8YEZ as carrying 32 pA and making a 1.5 mM nanodomain. The failure
    is silent: every number downstream stays finite and plausible.
    """
    from piezo1.analysis.report import ANALYSES
    from piezo1.config import STRUCTURE_DIR
    from piezo1.structure.fusion import HALOTAG_PDB

    st = structure_by_id("8YEZ")
    if st is None or not (STRUCTURE_DIR / f"{HALOTAG_PDB}.cif").exists():
        pytest.skip("8YEZ or 6U32 not downloaded")
    if not (STRUCTURE_DIR / "11ZC.cif").exists():
        pytest.skip("11ZC not downloaded to supply an open-state current")

    result = ANALYSES["nanodomain"](st, "human")
    assert "error" not in result, result

    # 8YEZ is closed, so the current must be borrowed and labelled as borrowed.
    assert result["current_source"] != "8YEZ"
    assert "11ZC" in result["current_source"]
    assert result["unitary_current_pA"] == pytest.approx(2.46, rel=0.1)
    assert result["calcium_uM"] == pytest.approx(114.0, rel=0.05)

    # And the same analysis run directly must agree that 8YEZ carries nothing.
    assert ANALYSES["permeation"](st, "human")["conductance_pS"] == 0.0
