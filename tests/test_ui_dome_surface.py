"""The drawn dome, and the surface that was deliberately not drawn.

A number cannot show whether a sphere was fitted to the right atoms. A picture
can, so the check that matters here is geometric: the drawn cap has to pass
through the transmembrane mid-points it was fitted to, and if it ever stops
doing that the measurement has drifted from what it claims to measure.

The second half is about restraint. The first version also drew the far-field
Helfrich footprint, which at PIEZO1's contact slope plunges 158 Å over a 526 Å
skirt and overestimates the real thing 3.65x. It is not drawn now, and this
records why so nobody adds it back as an improvement.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.structure.frame import apply_frame, canonical_transform
from piezo1.structure.geometry import measure_dome, tm_surface_points
from piezo1.structure.protomers import protomer_blocks


def _dome(pdb: str, reference: str):
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    st = Structure.from_file(path)
    st = apply_frame(st, canonical_transform(st))
    blocks, _ = protomer_blocks(st)
    points, resolved = tm_surface_points(st, reference)
    return measure_dome(blocks, points), points, resolved


def _drawn_height(geometry):
    """The z(r) the controller draws, reproduced without a GL context."""
    from piezo1.ui.dome_controller import DomeController

    controller = DomeController.__new__(DomeController)
    return DomeController._height(controller, geometry)


@pytest.mark.parametrize("pdb,reference", [("8YEZ", "human"), ("7WLT", "mouse")])
def test_the_drawn_cap_passes_through_the_atoms_it_was_fitted_to(pdb, reference):
    """The check a picture makes and a radius does not.

    A sphere fitted to the wrong atoms returns a perfectly reasonable radius of
    curvature. What it cannot do is pass through the transmembrane helices, so
    that is what is asserted — against the fit's own RMSE, not a round number.
    """
    geometry, points, _ = _dome(pdb, reference)
    height, sign, _centre = _drawn_height(geometry)

    axis = geometry.axis
    direction = axis.direction / np.linalg.norm(axis.direction)
    along = (points - axis.point) @ direction
    radial = np.linalg.norm((points - axis.point)
                            - np.outer(along, direction), axis=1)

    miss = np.abs(height(radial) - along)
    assert miss.mean() < 3.0 * geometry.notes["sphere_rmse"], (
        f"the drawn cap sits {miss.mean():.1f} A from the helices it was "
        f"fitted to, against a fit RMSE of {geometry.notes['sphere_rmse']:.1f}")


def test_the_side_of_the_sphere_is_taken_from_the_data():
    """A sphere has two caps and only one of them is the membrane.

    Choosing by the frame's convention would be an assumption; the controller
    takes it from the measured radial profile. The wrong choice misses by two
    orders of magnitude more, which is what makes this a real decision.
    """
    geometry, points, _ = _dome("8YEZ", "human")
    height, sign, centre = _drawn_height(geometry)

    axis = geometry.axis
    direction = axis.direction / np.linalg.norm(axis.direction)
    along = (points - axis.point) @ direction
    radial = np.linalg.norm((points - axis.point)
                            - np.outer(along, direction), axis=1)

    good = np.abs(height(radial) - along).mean()
    inside = np.clip(geometry.sphere.radius ** 2 - radial ** 2, 0.0, None)
    wrong = np.abs(centre - sign * np.sqrt(inside) - along).mean()
    assert good < wrong / 10.0, (good, wrong)


def test_the_projection_disc_is_the_surface_the_excess_area_is_measured_against():
    """The grey disc is not decoration: it is the denominator of the claim."""
    geometry, _points, _ = _dome("8YEZ", "human")
    height, _sign, _centre = _drawn_height(geometry)

    rim_radius = min(float(geometry.sphere.radius),
                     float(geometry.footprint_radius))
    rim = float(height(rim_radius))
    apex = float(height(0.0))

    # The cap rises above its own projection by the measured dome depth.
    assert abs(abs(apex - rim) - geometry.dome_depth) < 0.15 * geometry.dome_depth
    assert geometry.excess_area > 0


def test_the_far_field_footprint_is_not_drawn_and_the_reason_is_recorded():
    """Restraint, pinned so it is not undone as an improvement.

    The linearised Helfrich solution is a small-slope expansion. PIEZO1's cap
    meets the membrane at a slope near 1.9, which is not small, and Round 18
    measured the resulting 3.65x overestimate.
    """
    from piezo1.ui import dome_controller

    geometry, _points, _ = _dome("8YEZ", "human")
    height, _sign, _centre = _drawn_height(geometry)
    rim_radius = min(float(geometry.sphere.radius),
                     float(geometry.footprint_radius))
    step = rim_radius * 1e-3
    slope = abs((height(rim_radius) - height(rim_radius - step)) / step)
    assert slope > 1.0, (
        f"contact slope {slope:.2f}; if this has become small the linear "
        f"theory would be valid and the skirt could be drawn after all")

    assert not hasattr(dome_controller, "FOOTPRINT_COLOR")
    assert "3.65" in dome_controller.__doc__
    assert "PROJECTION_COLOR" in dome_controller.__all__
