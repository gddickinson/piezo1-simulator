"""The idealised dome, against the arithmetic Guo & MacKinnon published.

Every number in Figure 7 and Figure 7-figure supplement 1 follows from two
lengths by closed-form spherical-cap geometry, so this file is a calibration
before it is a test: the answers are known independently of any code here.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.physics.dome_idealised import (PUBLISHED_FIGURE7, IdealisedDome,
                                           compare_with_measured,
                                           flattening_series, guo2017_dome,
                                           spherical_cap_from_measurement)
from piezo1.parameters import PARAMETERS


def test_every_published_figure7_number_is_reproduced():
    """The whole of Figure 7-figure supplement 1, from R and c alone."""
    dome = guo2017_dome()
    values = dome.as_dict()
    values["stabilisation_kT"] = dome.stabilisation(
        PARAMETERS.value("dome.lytic_tension") / 10.0)
    values["lytic_tension_kT_per_nm2"] = PARAMETERS.value("dome.lytic_tension")

    checked = 0
    for key, (published, tolerance, source) in PUBLISHED_FIGURE7.items():
        if key not in values:
            continue
        checked += 1
        assert abs(values[key] - published) <= tolerance, (
            f"{key}: {values[key]:.3f} against the published {published} "
            f"({source})")
    assert checked >= 10, "the published table should be substantially covered"


def test_the_cap_relations_are_exact_on_a_known_hemisphere():
    """A hemisphere: centre on the plane, area 2 pi R^2, projection pi R^2.

    The one cap whose numbers can be written down without trusting anything.
    """
    dome = IdealisedDome(radius=5.0, center_height=0.0)
    assert dome.rim_radius == pytest.approx(5.0)
    assert dome.depth == pytest.approx(5.0)
    assert dome.area == pytest.approx(2 * np.pi * 25.0)
    assert dome.projected_area == pytest.approx(np.pi * 25.0)
    assert dome.excess_area == pytest.approx(np.pi * 25.0)
    assert np.degrees(dome.polar_angle) == pytest.approx(90.0)
    assert dome.contact_slope == float("inf")


def test_bending_energy_is_the_helfrich_sphere_result():
    """E = 2 kappa A / R^2, checked against a whole sphere's 8 pi kappa.

    A complete sphere has A = 4 pi R^2, so E = 8 pi kappa regardless of
    radius — the textbook scale-invariance of Helfrich bending, and a check
    that cannot be satisfied by an expression with the wrong power of R.
    """
    kappa = 20.0
    for radius in (5.0, 10.2, 40.0):
        # Two hemispheres make a sphere.
        hemisphere = IdealisedDome(radius=radius, center_height=0.0,
                                   kappa=kappa)
        assert 2 * hemisphere.bending_energy == pytest.approx(8 * np.pi * kappa)


def test_a_sphere_that_never_reaches_the_plane_is_refused():
    """No cap exists, and a clamped one would report a plausible area."""
    with pytest.raises(ValueError, match="never meets the plane"):
        IdealisedDome(radius=10.0, center_height=10.0)
    with pytest.raises(ValueError, match="never meets the plane"):
        IdealisedDome(radius=10.0, center_height=12.0)


def test_flattening_conserves_membrane_area_and_releases_bending_energy():
    """The constraint Figure 7c draws, checked at every point of the sweep."""
    dome = guo2017_dome()
    series = flattening_series(dome, n=25)

    for point in series:
        area = 2 * np.pi * point.radius_nm * point.depth_nm
        assert area == pytest.approx(dome.area, rel=1e-6), (
            "flattening must transfer area, not create it")

    angles = [p.polar_angle_deg for p in series]
    assert angles == sorted(angles, reverse=True)
    projected = [p.projected_area_nm2 for p in series]
    assert projected == sorted(projected), "flattening must grow the projection"
    bending = [p.bending_energy_kT for p in series]
    assert bending == sorted(bending, reverse=True), (
        "flattening must release bending energy")

    # In the flat limit the whole excess area is recovered and the whole
    # bending energy is released.
    assert series[-1].delta_projected_nm2 == pytest.approx(
        dome.excess_area, rel=1e-3)
    assert series[-1].delta_bending_kT == pytest.approx(
        -dome.bending_energy, rel=1e-3)


def test_flattening_free_energy_uses_both_terms():
    """Equation 3 at a point on the sweep, assembled by hand."""
    dome = guo2017_dome()
    point = flattening_series(dome, n=9)[4]
    tension, delta_g_prot = 0.35, 30.0
    expected = (delta_g_prot + point.delta_bending_kT
                - tension * point.delta_projected_nm2)
    assert point.free_energy(tension, delta_g_prot) == pytest.approx(expected)


def test_the_measured_dome_does_not_match_the_idealised_one(structure_6b3r):
    """And the comparison says so rather than reconciling them.

    The idealisation is a shape chosen for tractability. Ours integrates the
    real radial profile. Both are reported; if this ever starts agreeing, the
    measurement has been fitted to the paper rather than to the coordinates.
    """
    from piezo1.structure.geometry import measure_dome, tm_surface_points
    from piezo1.structure.protomers import protomer_blocks

    blocks, _ = protomer_blocks(structure_6b3r)
    points, _ = tm_surface_points(structure_6b3r, "mouse")
    comparison = compare_with_measured(measure_dome(blocks, points))

    assert comparison["idealised"]["dome_area_nm2"] == pytest.approx(397.3,
                                                                    abs=1.0)
    assert comparison["measured"]["dome_area_nm2"] > 450.0
    assert "idealisation" in comparison["caveat"]
    # The radius is the one thing that should be close: it is what the
    # idealisation was chosen to match.
    assert 0.8 < comparison["measured_over_idealised"]["radius_nm"] < 1.3


def test_a_measurement_that_is_not_a_cap_is_refused_not_reconciled():
    """A footprint wider than the sphere has no cap form."""
    with pytest.raises(ValueError, match="no spherical cap"):
        spherical_cap_from_measurement(radius_nm=10.0, rim_radius_nm=12.0)


def test_the_paper_numbers_follow_the_registry():
    """Overriding the radius must move every derived number with it."""
    try:
        PARAMETERS.set_value("dome.published_radius_closed", 12.0)
        moved = guo2017_dome()
        assert moved.radius == pytest.approx(12.0)
        assert moved.area > 397.3
    finally:
        PARAMETERS.reset()
    assert guo2017_dome().radius == pytest.approx(10.2)
