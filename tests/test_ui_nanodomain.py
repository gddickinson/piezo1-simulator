"""The drawn calcium field, and the surface that is deliberately not drawn.

Round 32's conclusion is a distance: the sensor is saturated whenever its own
channel opens, so puncta brightness cannot report calcium amplitude. The
surfaces that carry that conclusion are 119 nm and 372 nm across a channel
reaching about 15 nm, so drawing them would leave a viewport of shell and a
speck of protein. They are reported as numbers instead — the same rule
`dome_controller` applies to the far-field membrane footprint.

That makes the *budget* the load-bearing part of this feature, so it is
measured here rather than trusted: the near shells must fall inside it and the
occupancy surfaces must fall outside it, on the one open-like entry the project
has. A budget that admitted everything, or nothing, would produce a picture
that looks the same either way.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.parameters import PARAMETERS as _P  # noqa: E402
from piezo1.physics.nanodomain import (Nanodomain, calcium_at,  # noqa: E402
                                       saturation)
from piezo1.render.geometry_builders import build_sphere  # noqa: E402
from piezo1.ui.nanodomain_controller import (KD_OCCUPANCY,  # noqa: E402
                                             MIN_SHELL_A, NAME,
                                             NanodomainController,
                                             SATURATION_OCCUPANCY,
                                             SCENE_MULTIPLE,
                                             SHELL_CONCENTRATIONS)


# --------------------------------------------------------------------------
# The primitive, on a shape whose answer is arithmetic
# --------------------------------------------------------------------------

def test_the_sphere_mesh_is_a_sphere():
    centre = np.array([3.0, -4.0, 12.0])
    mesh = build_sphere(centre, 7.5)
    radii = np.linalg.norm(mesh.positions - centre, axis=1)
    assert np.allclose(radii, 7.5, atol=1e-4)
    assert np.allclose(np.linalg.norm(mesh.normals, axis=1), 1.0, atol=1e-5)
    # Outward normals: a mesh lit from inside looks like a hole.
    outward = (mesh.positions - centre)
    outward /= np.linalg.norm(outward, axis=1, keepdims=True)
    assert np.allclose((outward * mesh.normals).sum(axis=1), 1.0, atol=1e-4)
    assert mesh.n_triangles > 100
    assert mesh.indices.max() < mesh.n_vertices


# --------------------------------------------------------------------------
# The window stubs
# --------------------------------------------------------------------------

class _Batch:
    def __init__(self):
        self.args = None

    def upload(self, *args):
        self.args = args


class _Scene:
    """Records what was uploaded, so drawing can be checked without a GPU."""

    def __init__(self):
        self.batches = {}

    def mesh(self, name, **kw):
        return self.batches.setdefault(name, _Batch())

    def spheres(self, name):
        return self.batches.setdefault(name, _Batch())

    def remove(self, name):
        self.batches.pop(name, None)


class _FakeWindow:
    def __init__(self, structure, profile, hydration=None):
        self.structure = structure
        self.analysis = type("A", (), {"pore": profile, "hydration": hydration,
                                       "compute_pore": lambda self: None})()
        self.fusion = type("F", (), {"model": None})()
        self.viewport = type("V", (), {"scene": _Scene(),
                                       "update": lambda self: None})()
        self.status = ""

    def _set_status(self, text):
        self.status = text


def _framed(pdb):
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    from piezo1.structure.frame import apply_frame, canonical_transform

    st = Structure.from_file(path)
    return apply_frame(st, canonical_transform(st))


@pytest.fixture(scope="module")
def open_controller(open_profile):
    controller = NanodomainController(_FakeWindow(_framed("11ZC"), open_profile))
    controller.show(True)
    if controller.model is None:
        pytest.skip(f"11ZC gave no current: {controller.win.status}")
    return controller


# --------------------------------------------------------------------------
# The shells are the inverse of the field they claim to be
# --------------------------------------------------------------------------

def test_each_shell_sits_where_the_field_reaches_its_concentration(open_controller):
    """Round-trip through `calcium_at`, which the shells never call.

    The radii come out of `distance_for_occupancy`, an inverter. Evaluating
    the forward function at the radius it returns is the check that the two
    are inverses of each other rather than of something else — and the
    occupancy detour is where that could go wrong, since the shells are
    specified in concentration and solved for in occupancy.
    """
    calcium = open_controller.model.calcium_current_A
    assert open_controller.shells, "no shell radius was computed at all"
    for concentration, radius_A in open_controller.shells.items():
        recovered = float(calcium_at(radius_A * 1e-10, calcium))
        assert recovered == pytest.approx(concentration, rel=1e-3), (
            f"the {concentration * 1e6:.0f} uM shell was placed where the "
            f"field is {recovered * 1e6:.1f} uM")


def test_the_shells_fall_off_with_distance(open_controller):
    radii = [open_controller.shells[c] for c in
             sorted(open_controller.shells, reverse=True)]
    assert radii == sorted(radii), \
        "a higher concentration was placed further from the source"


# --------------------------------------------------------------------------
# The budget: the part that decides what the picture is of
# --------------------------------------------------------------------------

def test_the_occupancy_surfaces_really_are_too_big_to_draw(open_controller):
    """The measurement behind "not drawn". If they ever fit, draw them."""
    origin = open_controller.source_point()
    limit = SCENE_MULTIPLE * open_controller.scene_radius(origin)
    for occupancy in (SATURATION_OCCUPANCY, KD_OCCUPANCY):
        radius = open_controller.occupancy_radii[occupancy]
        assert radius > limit, (
            f"the {occupancy:.0%} surface is {radius / 10:.0f} nm against a "
            f"{limit / 10:.0f} nm budget — it now fits and should be drawn "
            f"rather than reported")
    assert open_controller.occupancy_radii[KD_OCCUPANCY] > \
        open_controller.occupancy_radii[SATURATION_OCCUPANCY]


def test_the_budget_still_admits_something(open_controller):
    """A filter that excludes everything is indistinguishable from a bug."""
    drawn = open_controller.drawable()
    assert 2 <= len(drawn) <= len(SHELL_CONCENTRATIONS)
    origin = open_controller.source_point()
    limit = SCENE_MULTIPLE * open_controller.scene_radius(origin)
    for _concentration, radius in drawn:
        assert MIN_SHELL_A <= radius <= limit


def test_the_drawn_shells_are_the_ones_uploaded(open_controller):
    scene = open_controller.win.viewport.scene
    names = {f"{NAME}:{c:g}" for c, _ in open_controller.drawable()}
    assert names <= set(scene.batches)
    assert f"{NAME}:source" in scene.batches
    for name in names:
        positions, _normals, _colors, _indices, alpha = scene.batches[name].args
        assert 0.0 < alpha < 0.5, \
            "an iso-surface was uploaded opaque and would hide the protein"
        assert len(positions) > 100


def test_the_source_marker_is_the_cytosolic_mouth_not_the_other_end(open_controller):
    """The sign is the whole answer, so it is measured — and measurably so.

    Flipping the structure must move the marker to the other end of the
    profile. If it did not, the marker would be the frame's convention rather
    than a property of the protein.
    """
    from piezo1.physics.pore_charge import cytosolic_end

    profile = open_controller.win.analysis.pore
    structure = open_controller.win.structure
    end = cytosolic_end(structure, profile.axis)
    assert np.allclose(open_controller.source_point(), profile.centers[end])

    flipped = structure.copy_with_coords(
        structure.xyz * np.array([1.0, 1.0, -1.0], dtype=structure.xyz.dtype))
    assert cytosolic_end(flipped, profile.axis) != end, (
        "turning the structure upside down did not move the cytosolic end, "
        "so the source marker is a frame convention rather than a measurement")


# --------------------------------------------------------------------------
# The threshold the picture shares with the boolean
# --------------------------------------------------------------------------

def test_the_saturation_surface_is_where_saturated_stops_being_true(open_controller):
    """`Nanodomain.saturated` is a literal 0.9 in the physics module and
    `SATURATION_OCCUPANCY` is a literal 0.9 here. Two literals agreeing today
    is not the same as agreeing tomorrow, so the boolean is driven across the
    drawn surface rather than the numbers compared."""
    model = open_controller.model
    radius_m = open_controller.occupancy_radii[SATURATION_OCCUPANCY] * 1e-10
    inside = Nanodomain(current_A=model.current_A,
                        calcium_fraction=model.calcium_fraction,
                        distance_m=radius_m * 0.9)
    outside = Nanodomain(current_A=model.current_A,
                         calcium_fraction=model.calcium_fraction,
                         distance_m=radius_m * 1.1)
    assert inside.saturated
    assert not outside.saturated


def test_the_kd_surface_is_the_registered_kd(open_controller):
    calcium = open_controller.model.calcium_current_A
    radius_m = open_controller.occupancy_radii[KD_OCCUPANCY] * 1e-10
    here = float(calcium_at(radius_m, calcium))
    assert here == pytest.approx(_P.value("nanodomain.sensor_kd"), rel=1e-2)
    assert float(saturation(here)) == pytest.approx(0.5, abs=1e-3)


# --------------------------------------------------------------------------
# The shut structure, and the current that is not borrowed
# --------------------------------------------------------------------------

def test_a_shut_structure_draws_nothing_and_borrows_nothing(human_structure):
    """Round 34's null, in the one place borrowing would be invisible.

    `report_tags` substitutes 11ZC's current when the loaded entry is shut and
    labels the substitution. A picture cannot carry that label convincingly, so
    here the answer is an empty screen and a sentence.
    """
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.superpose import detect_c3_axis

    st = apply_frame(human_structure, canonical_transform(human_structure))
    blocks, _ = protomer_blocks(st)
    profile = pore_profile(st, detect_c3_axis(blocks), step=1.0)

    controller = NanodomainController(_FakeWindow(st, profile))
    controller.show(True)
    assert not controller.visible
    assert controller.model is None
    assert not controller.win.viewport.scene.batches
    status = controller.win.status
    assert "no current" in status and "11ZC" in status
    assert "NOT borrowed" in status


def test_the_status_line_says_the_model_does_not_know_the_protein_is_there(
        open_controller):
    line = open_controller.status_line()
    assert "point source in free solution" in line
    assert "NEITHER is drawn" in line
    assert "pA" in line


def test_nothing_computed_says_so_rather_than_raising():
    controller = NanodomainController(_FakeWindow(None, None))
    assert controller.status_line() == "no calcium nanodomain"
    assert controller.drawable() == []
