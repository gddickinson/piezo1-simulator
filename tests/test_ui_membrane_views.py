"""The three Figure 4 views: micelle, planar membrane, electrostatic colouring.

Each of the three is a picture that can be mistaken for a measurement it is
not, and each is tested for that specifically:

* the **micelle** is a construction, not the density map Figure 4b shows;
* the **planar membrane** is a least-squares fit, and every point set has one;
* the **electrostatic colouring** is not APBS, and its scale is fixed so that
  an almost-neutral protein cannot be painted as violently charged.

The micelle's geometry is calibrated on shapes whose envelope is known in
closed form — a point gives a sphere of the offset radius, a line gives a
capsule — before it is run on a structure.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.structure.micelle import (APOLAR_RESIDUES, belt_atoms,
                                      build_micelle, distance_field)


def _contour(points, offset, spacing=0.5, pad=None):
    from skimage import measure

    field, low = distance_field(np.asarray(points, dtype=float), spacing,
                                pad=(offset + 4.0) if pad is None else pad)
    vertices, faces, _, _ = measure.marching_cubes(
        field, level=float(offset), spacing=(spacing,) * 3)
    return vertices + low, faces


def _area(vertices, faces) -> float:
    a, b, c = (vertices[faces[:, i]] for i in range(3))
    return float(0.5 * np.linalg.norm(np.cross(b - a, c - a), axis=1).sum())


# --------------------------------------------------------------------------
# The envelope, on shapes with a closed-form answer
# --------------------------------------------------------------------------

def test_the_envelope_of_one_point_is_a_sphere_of_the_offset_radius():
    """The calibration everything else rests on."""
    for offset in (5.0, 8.0, 12.0):
        vertices, faces = _contour(np.zeros((1, 3)), offset)
        radii = np.linalg.norm(vertices, axis=1)
        assert radii.mean() == pytest.approx(offset, abs=0.02)
        assert radii.std() < 0.05
        assert _area(vertices, faces) == pytest.approx(
            4 * np.pi * offset ** 2, rel=0.01)


def test_the_envelope_of_a_line_is_a_capsule():
    """Area = cylinder + two hemispherical caps, known exactly."""
    length, offset = 30.0, 6.0
    line = np.column_stack([np.linspace(0, length, 400),
                            np.zeros(400), np.zeros(400)])
    vertices, faces = _contour(line, offset, spacing=0.4)
    expected = 2 * np.pi * offset * length + 4 * np.pi * offset ** 2
    assert _area(vertices, faces) == pytest.approx(expected, rel=0.02)


def test_the_offset_adds_to_a_sphere_rather_than_reshaping_it():
    """The reason the curvature is reported and the thickness is not.

    An offset surface around a sphere is a sphere with the radius increased by
    exactly the offset, so the shell thickness carries no shape information —
    which is what lets the module call the curvature a measurement and the
    thickness a parameter.
    """
    rng = np.random.default_rng(0)
    direction = rng.normal(size=(3000, 3))
    direction /= np.linalg.norm(direction, axis=1, keepdims=True)
    shell = 20.0 * direction
    for offset in (5.0, 10.0):
        vertices, _ = _contour(shell, offset, spacing=0.6)
        outer = np.linalg.norm(vertices, axis=1)
        # Marching cubes finds both the outer and the inner surface here.
        assert outer.max() == pytest.approx(20.0 + offset, abs=0.3)


def test_an_empty_point_set_is_refused():
    with pytest.raises(ValueError, match="no points"):
        distance_field(np.zeros((0, 3)), 1.0, pad=2.0)


# --------------------------------------------------------------------------
# On a structure
# --------------------------------------------------------------------------

def test_the_belt_is_apolar_transmembrane_side_chains(structure_6b3r):
    mask = belt_atoms(structure_6b3r, "mouse")
    assert mask.sum() > 1000
    assert set(np.unique(structure_6b3r.res_name[mask])) <= APOLAR_RESIDUES
    assert not np.isin(structure_6b3r.atom_name[mask], ("N", "C", "O")).any(), (
        "the backbone runs through the helix core and is not the surface "
        "detergent packs against")


def test_the_belt_moves_with_the_numbering_system(structure_6b3r):
    """A human entry read with mouse ranges would give a plausible band."""
    mouse = belt_atoms(structure_6b3r, "mouse")
    human = belt_atoms(structure_6b3r, "human")
    assert mouse.sum() != human.sum() or not np.array_equal(mouse, human), (
        "if the two references gave identical belts, the numbering guard "
        "would be untested")


def test_the_micelle_encloses_the_belt_and_says_it_is_a_model(structure_6b3r):
    envelope = build_micelle(structure_6b3r, "mouse")
    assert envelope.is_observed is False
    assert "not the observed density" in envelope.caveat.lower()
    assert envelope.n_vertices > 1000 and len(envelope.faces) > 1000

    # Every belt atom must be inside: the envelope is the offset surface, so
    # nothing it was built around can be outside it.
    from scipy.spatial import cKDTree
    belt = structure_6b3r.xyz[belt_atoms(structure_6b3r, "mouse")]
    tree = cKDTree(belt)
    distance, _ = tree.query(envelope.vertices, k=1)
    assert distance.min() > envelope.offset - 1.0
    assert distance.max() < envelope.offset + 1.0, (
        "every vertex must sit at the offset distance; a vertex further out "
        "means the contour caught a second component")


def test_the_micelle_curvature_is_close_to_the_published_idealisation(
        structure_6b3r):
    """The one number in the picture that is a measurement of the protein."""
    envelope = build_micelle(structure_6b3r, "mouse")
    assert envelope.sphere is not None
    radius_nm = envelope.sphere.radius / 10.0
    assert 8.0 < radius_nm < 13.0, (
        f"belt curvature {radius_nm:.1f} nm against the paper's 10.2 nm "
        f"idealisation")


def test_the_thickness_does_not_move_the_curvature(structure_6b3r):
    """The claim the status line makes, checked rather than asserted."""
    a = build_micelle(structure_6b3r, "mouse", offset=7.0)
    b = build_micelle(structure_6b3r, "mouse", offset=13.0)
    assert a.sphere.radius == pytest.approx(b.sphere.radius, rel=1e-9)
    assert b.enclosed_volume() > a.enclosed_volume(), (
        "a thicker shell must enclose more")


def test_a_structure_without_a_resolved_belt_is_refused(structure_6b3r):
    """Rather than returning a handful of disconnected blobs that still draw."""
    nothing = np.zeros(structure_6b3r.n_atoms, dtype=bool)
    with pytest.raises(ValueError, match="does not resolve enough"):
        build_micelle(structure_6b3r, "mouse", mask=nothing)


# --------------------------------------------------------------------------
# The potential colouring
# --------------------------------------------------------------------------

def test_the_potential_scale_is_fixed_not_auto_ranged():
    """The whole reason ColorBy.POTENTIAL exists beside ColorBy.VALUE.

    An auto-ranged map paints an almost-neutral protein in full saturation and
    cannot be compared with a published surface or with another state.
    """
    from piezo1.render.colormaps import potential_colors, value_colors

    tiny = np.array([-0.02, 0.0, 0.02])
    fixed = potential_colors(tiny, scale=5.0)
    # Nearly neutral input must give nearly white output.
    assert np.all(fixed > 0.97)
    # The auto-ranging map does the opposite, which is what is being avoided.
    auto = value_colors(tiny)
    assert np.ptp(auto, axis=0).max() > 0.3


def test_the_potential_colours_follow_the_published_convention():
    from piezo1.render.colormaps import potential_colors

    colours = potential_colors(np.array([-5.0, 0.0, 5.0]), scale=5.0)
    negative, neutral, positive = colours
    assert negative[0] > negative[2], "negative must be red"
    assert positive[2] > positive[0], "positive must be blue"
    assert np.all(neutral > 0.95), "zero must be white"
    # Beyond the scale it saturates rather than continuing.
    beyond = potential_colors(np.array([-50.0, 50.0]), scale=5.0)
    assert np.allclose(beyond[0], negative) and np.allclose(beyond[1], positive)


def test_an_unmeasured_atom_is_neutral_rather_than_extreme():
    from piezo1.render.colormaps import potential_colors

    colours = potential_colors(np.array([np.nan]), scale=5.0)
    assert np.all(colours > 0.95)


def test_the_potential_mode_is_dispatched_separately_from_value(structure_6b3r):
    """Both modes read `view.values`, so the dispatch has to tell them apart.

    Driven through `atom_colors` rather than read off the source: an earlier
    version of this test inspected the bytecode with an `or True` in the
    assertion, which could not fail.
    """
    from piezo1.render.colormaps import potential_colors, value_colors
    from piezo1.render.representations import ColorBy, MolecularView

    # MolecularView needs a scene it never touches for this call.
    view = MolecularView(scene=None, structure=structure_6b3r)
    values = np.linspace(-8.0, 8.0, structure_6b3r.n_atoms)
    view.values = values

    view.color_by = ColorBy.POTENTIAL
    assert np.allclose(view.atom_colors(), potential_colors(values))

    view.color_by = ColorBy.VALUE
    assert np.allclose(view.atom_colors(), value_colors(values))

    assert not np.allclose(potential_colors(values), value_colors(values)), (
        "if the two maps agreed, the dispatch would be untested")


# --------------------------------------------------------------------------
# The controllers, without a GL context
# --------------------------------------------------------------------------

class _Window:
    def __init__(self, structure=None):
        self.structure = structure
        self.record = type("R", (), {"numbering_species": "mouse"})()
        self.viewport = type("V", (), {"scene": None,
                                       "update": lambda self: None})()
        self.view = None
        self.status = ""

    def _set_status(self, text):
        self.status = text

    def _current_color(self):
        from piezo1.render.representations import ColorBy
        return ColorBy.CHAIN


def test_the_micelle_status_line_cannot_omit_that_it_is_modelled(
        structure_6b3r):
    from piezo1.ui.micelle_controller import MicelleController

    controller = MicelleController(_Window(structure_6b3r))
    controller.envelope = build_micelle(structure_6b3r, "mouse")
    line = controller.status_line()
    assert "MODELLED" in line
    assert "NOT THE OBSERVED DENSITY" in line
    assert "parameter" in line and "curvature" in line


def test_the_planar_membrane_status_line_reports_the_residual(structure_6b3r):
    """Not the lines — every point set has a best-fit plane."""
    from piezo1.structure.planarity import planarity
    from piezo1.ui.planar_membrane_controller import PlanarMembraneController

    controller = PlanarMembraneController(_Window(structure_6b3r))
    controller.comparison = planarity(structure_6b3r, "mouse")
    controller.chain = sorted(controller.comparison.per_protomer)[0]
    line = controller.status_line()
    assert "RMS departure" in line
    assert "read the residual, not the lines" in line
    assert "slab" in line


def test_the_potential_status_line_leads_with_not_apbs(structure_6b3r):
    from piezo1.physics.electrostatics import surface_potential
    from piezo1.ui.potential_controller import ElectrostaticColourController

    controller = ElectrostaticColourController(_Window(structure_6b3r))
    controller.result = surface_potential(structure_6b3r)
    line = controller.status_line()
    assert line.startswith("NOT APBS")
    assert "UNIFORM dielectric" in line
    assert "no accessible surface" in line


def test_no_controller_draws_without_a_structure():
    from piezo1.ui.micelle_controller import MicelleController
    from piezo1.ui.planar_membrane_controller import PlanarMembraneController
    from piezo1.ui.potential_controller import ElectrostaticColourController

    for factory in (MicelleController, PlanarMembraneController,
                    ElectrostaticColourController):
        window = _Window(None)
        controller = factory(window)
        controller.show(True)
        assert not controller.visible
        assert "load a structure" in window.status
