"""Figure 4a: a protomer fits a plane, the trimer does not.

The instrument is calibrated on shapes whose planarity is known by
construction before it is pointed at a structure, and the decomposition it
reports is checked for closure — a residual split that does not add up has
missed a term and would attribute the difference to whichever half was named
first.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.structure.planarity import (PROXIMAL_FIRST_HELIX, beam_angle,
                                        blade_dependence, fit_plane, planarity)


# --------------------------------------------------------------------------
# The plane fitter, on shapes with known answers
# --------------------------------------------------------------------------

def test_a_plane_is_fitted_exactly_to_coplanar_points():
    rng = np.random.default_rng(0)
    xy = rng.normal(size=(200, 2)) * 10.0
    points = np.column_stack([xy, np.zeros(len(xy))])
    fit = fit_plane(points)
    assert fit.rmsd == pytest.approx(0.0, abs=1e-9)
    assert abs(abs(float(fit.normal[2])) - 1.0) < 1e-9


def test_the_fit_survives_a_plane_containing_the_z_axis():
    """A regression of z on (x, y) has no answer here; total least squares does.

    PIEZO1's blades approach the vertical, so this is not a hypothetical.
    """
    rng = np.random.default_rng(1)
    points = np.column_stack([np.zeros(200), rng.normal(size=200) * 10.0,
                              rng.normal(size=200) * 10.0])
    fit = fit_plane(points)
    assert fit.rmsd == pytest.approx(0.0, abs=1e-9)
    assert abs(abs(float(fit.normal[0])) - 1.0) < 1e-9


def test_a_known_out_of_plane_displacement_is_recovered():
    """Points on a cone of known half-angle give a known RMS."""
    angles = np.linspace(0, 2 * np.pi, 361)[:-1]
    radius, rise = 10.0, 3.0
    points = np.column_stack([radius * np.cos(angles), radius * np.sin(angles),
                              rise * np.cos(3 * angles)])
    fit = fit_plane(points)
    assert fit.rmsd == pytest.approx(rise / np.sqrt(2.0), rel=0.02)
    assert fit.max_deviation == pytest.approx(rise, rel=0.02)


def test_fewer_than_three_points_is_refused():
    with pytest.raises(ValueError, match="three points"):
        fit_plane(np.zeros((2, 3)))


# --------------------------------------------------------------------------
# The decomposition
# --------------------------------------------------------------------------

def test_the_decomposition_closes(structure_6b3r):
    """within^2 + arrangement^2 must account for the trimer residual.

    If it does not, a third contribution is being silently attributed to
    whichever of the two is reported first.
    """
    result = planarity(structure_6b3r, "mouse")
    assert abs(result.decomposition_residual) < 0.1 * result.trimer.rmsd


def test_a_curved_trimer_puts_the_non_planarity_in_the_arrangement(
        structure_6b3r):
    """6B3R: Figure 4a's claim, as residuals."""
    result = planarity(structure_6b3r, "mouse")
    assert result.supports_paper
    assert result.arrangement_rmsd > 2 * result.protomer_rmsd
    assert result.ratio > 2.0
    # "approximately 30 degrees out of the plane defined by the pore"
    assert 25.0 <= result.mean_tilt_deg <= 45.0


def test_the_flattened_structure_is_the_control_that_makes_this_mean_something(
        flat_structure):
    """7WLU is flattened, and the claim must fail on it.

    Without a structure the test says no to, "the trimer is less planar than
    a protomer" would be a property of trimers rather than of curvature.
    """
    result = planarity(flat_structure, "mouse")
    assert not result.supports_paper, (
        "the flattened structure must not read as curved")
    assert result.ratio < 1.5
    assert result.mean_tilt_deg < 20.0


def test_coverage_decides_the_answer_and_the_module_says_so(structure_by_id):
    """The trap: 6B3R and 6BPZ disagree only about what they resolve.

    Coverage-matched they agree to within an Angstrom. This is the same
    failure ``analysis/paralogue.py`` was written after.
    """
    a, b = structure_by_id("6B3R"), structure_by_id("6BPZ")
    if a is None or b is None:
        pytest.skip("6B3R and 6BPZ not both downloaded")

    from piezo1.structure.geometry import tm_surface_by_chain
    _, resolved_a = tm_surface_by_chain(a, "mouse")
    _, resolved_b = tm_surface_by_chain(b, "mouse")
    shared = resolved_a & resolved_b
    assert len(resolved_a) > len(shared), "6B3R should resolve more blade"

    naive_a = planarity(a, "mouse").arrangement_rmsd
    naive_b = planarity(b, "mouse").arrangement_rmsd
    assert naive_a > 3 * naive_b, "naively they look like different proteins"

    matched_a = planarity(a, "mouse", keep=shared).arrangement_rmsd
    matched_b = planarity(b, "mouse", keep=shared).arrangement_rmsd
    assert abs(matched_a - matched_b) < 1.0, (
        "coverage-matched, the difference must disappear")


def test_the_blade_carries_the_non_planarity(structure_6b3r):
    """And the split is per-structure, so no second entry can confound it."""
    split = blade_dependence(structure_6b3r, "mouse")
    assert split.split_at == PROXIMAL_FIRST_HELIX
    assert len(split.distal_helices) >= 8
    assert 0.5 < split.blade_share < 1.0
    assert split.proximal.arrangement_rmsd < split.full.arrangement_rmsd


def test_an_entry_with_no_distal_blade_refuses_rather_than_returning_zero(
        structure_by_id):
    """6BPZ resolves nothing distal, so the blade's share is unmeasurable."""
    entry = structure_by_id("6BPZ")
    if entry is None:
        pytest.skip("6BPZ not downloaded")
    split = blade_dependence(entry, "mouse")
    assert split.distal_helices == ()
    assert np.isnan(split.blade_share)
    assert "cannot be measured" in split.summary()


# --------------------------------------------------------------------------
# The beam
# --------------------------------------------------------------------------

def test_the_beam_sits_at_the_angle_the_paper_states(structure_6b3r):
    """'about 60 degrees instead of 90', and the arms 30 out of plane."""
    beam = beam_angle(structure_6b3r)
    assert len(beam.angle_deg) == 3
    assert beam.mean_deg == pytest.approx(60.0, abs=8.0)
    assert beam.out_of_plane_deg == pytest.approx(30.0, abs=8.0)
    assert beam.numbering == "mouse"


def test_the_beam_angle_opens_towards_90_when_the_channel_flattens(
        flat_structure, structure_6b3r):
    """The mechanism, as a comparison: flattening takes the beam towards 90.

    This is what makes the 60 degrees a measurement of curvature rather than
    of where the beam happens to sit.
    """
    curved = beam_angle(structure_6b3r).mean_deg
    flat = beam_angle(flat_structure).mean_deg
    assert flat > curved + 10.0


def test_a_range_in_the_wrong_numbering_raises_rather_than_answering(
        structure_6b3r):
    """A real but different set of residues would give a plausible angle."""
    with pytest.raises(ValueError, match="numbering"):
        beam_angle(structure_6b3r, residue_range=(10, 60))
