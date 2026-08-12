"""Does the morph actually end on the structure it says it ends on?

``tests/test_morph.py`` covers the interpolation itself, and it passes: the
trajectory hits both endpoints to 1e-9. None of it touches the controller, and
the controller is where the morph became wrong — because a morph is not drawn
by replacing coordinates but by *adding a displacement field* to the ones
already on screen, and a displacement field only means anything in the frame it
was measured in and against the sites it was indexed by.

Two independent faults, either sufficient on its own to make the endpoint not
be the endpoint:

* the path was built from the deposited file while the viewport shows the
  canonical frame — 180 degrees apart on 7WLT, so every displacement was
  applied backwards and the flattened endpoint landed further from 7WLU than
  the whole conformational change;
* atoms at a residue outside the shared basis were assigned a site by *residue
  number*, which tied every bound lipid to a C-alpha most of the length of the
  molecule away.

So the assertions here are made in the two ways the fault could not survive:
the drawn endpoint is compared with the deposited flat structure **as a shape**
(optimal superposition, so no assumption about frame enters), and the path's
first frame is compared with the displayed coordinates **by value**. Each is
paired with a case that must fail, because a check on a 19.7 A change that
cannot distinguish the two endpoints is asserting nothing.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
moderngl = pytest.importorskip("moderngl")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.config import RenderSettings, STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.structure.superpose import match_protomers, superpose  # noqa: E402

PAIR = ("7WLT", "7WLU")          # curved/closed → flattened
#: 7WLT against 7WLU once the protomer correspondence is fixed. The scale every
#: tolerance here is judged against: a check that cannot resolve this is not
#: measuring the endpoint.
CHANGE_RMSD = 19.7


@pytest.fixture(scope="module")
def app():
    instance = QApplication.instance()
    if instance is None:
        try:
            instance = QApplication([])
        except Exception as exc:                       # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return instance


@pytest.fixture(scope="module")
def window(app):
    """A real MainWindow with a real GL context, as the user runs it."""
    for pdb in PAIR:
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb}.cif not downloaded — "
                        f"run python -m piezo1.io.fetch")
    from piezo1.render.scene import Scene
    from piezo1.ui.gl_widget import configure_surface_format
    from piezo1.ui.main_window import MainWindow

    configure_surface_format()
    win = MainWindow()
    win.resize(900, 700)
    win.show()
    app.processEvents()
    try:
        ctx = moderngl.create_standalone_context(require=410)
    except Exception as exc:                           # pragma: no cover
        pytest.skip(f"no OpenGL 4.1 context available: {exc}")
    scene = Scene(ctx, RenderSettings(samples=1))
    scene.resize(900, 700)
    win.viewport.scene = scene
    win._on_scene_ready(scene)
    app.processEvents()
    if win.structure is None:
        pytest.skip("no default structure could be loaded")
    yield win

    for name in ("analysis", "physics", "overlay"):
        cleanup = getattr(getattr(win, name, None), "cleanup", None)
        if cleanup is not None:
            cleanup()
    app.processEvents()


@pytest.fixture(scope="module")
def built(window, app):
    """The morph the Physics panel builds, on the window that displays it."""
    window.morph_controller.build(
        {"start": PAIR[0], "end": PAIR[1], "method": "restrained"})
    app.processEvents()
    if window.morph_controller.trajectory is None:
        pytest.fail(f"morph not built: {window.status_label.text()}")
    return window.morph_controller


# --------------------------------------------------------------------------
# Helpers — deliberately built from `structure.superpose` and the raw file, so
# nothing the controller computes is used to judge what the controller drew.
# --------------------------------------------------------------------------

def ca_blocks(st: Structure, residues: np.ndarray) -> list[np.ndarray]:
    """C-alpha coordinates at ``residues``, one block per protomer chain."""
    ca = st.mask_ca()
    chains = [ch for ch in st.chains if (ca & (st.chain == ch)).sum() > 300][:3]
    out = []
    for ch in chains:
        sel = ca & (st.chain == ch)
        seq, xyz = st.res_seq[sel], st.xyz[sel]
        out.append(xyz[np.searchsorted(seq, residues)].astype(np.float64))
    return out


def flat_reference(ctl) -> np.ndarray:
    """The deposited flat structure on the morph's own site basis and order."""
    flat = Structure.from_file(STRUCTURE_DIR / f"{PAIR[1]}.cif")
    mobile = ca_blocks(flat, ctl.residues)
    target = ca_blocks(ctl.win.structure, ctl.residues)
    order = match_protomers(mobile, target).order
    return np.vstack([mobile[i] for i in order])


def drawn_sites(ctl, fraction: float) -> np.ndarray:
    """Where the C-alpha of every morph site ends up on screen.

    Read off the coordinates the view is holding — the ones a rendered pixel
    comes from — rather than off anything the controller reports.
    """
    ctl.show_frame(fraction)
    return np.vstack(ca_blocks(ctl.win.view.structure, ctl.residues))


def shape_difference(a: np.ndarray, b: np.ndarray) -> float:
    """RMSD after optimal superposition: a statement about shape, not place."""
    return superpose(a, b)[1]


# --------------------------------------------------------------------------
# The path must live in the frame the viewport is in
# --------------------------------------------------------------------------

def test_the_path_starts_from_the_coordinates_on_screen(built):
    """`frames[0]` must be the displayed C-alphas, by value.

    This is the invariant the whole drawing method rests on. It held against
    the *file* and failed against the viewport: 7WLT's canonical frame is a
    180-degree rotation of its deposited one, and the two centroids are 381 A
    apart, so the field was added to coordinates it had never been measured on.
    """
    ctl = built
    displayed = np.vstack(ca_blocks(ctl.win.structure, ctl.residues))
    assert np.allclose(ctl.trajectory.frames[0], displayed, atol=1e-4)


def test_a_path_from_the_deposited_file_would_not_pass_that(built):
    """The calibration: the old composition must fail the test above.

    Without this, `test_the_path_starts_from_the_coordinates_on_screen` would
    be equally satisfied by a check that cannot tell two frames apart.
    """
    ctl = built
    deposited = Structure.from_file(STRUCTURE_DIR / f"{PAIR[0]}.cif")
    from_file = np.vstack(ca_blocks(deposited, ctl.residues))
    displayed = np.vstack(ca_blocks(ctl.win.structure, ctl.residues))
    offset = np.linalg.norm(from_file.mean(0) - displayed.mean(0))
    assert offset > 100.0, "the frames must really differ for this to be a test"
    assert not np.allclose(from_file, displayed, atol=1e-4)


# --------------------------------------------------------------------------
# The endpoints must be the endpoints
# --------------------------------------------------------------------------

def test_the_last_frame_is_the_flattened_structure(built):
    """At the far end of the slider, what is drawn *is* 7WLU — as a shape.

    Compared after optimal superposition, so the assertion says nothing about
    where the model sits and everything about whether it got there.
    """
    ctl = built
    assert shape_difference(drawn_sites(ctl, 1.0), flat_reference(ctl)) < 0.05


def test_the_first_frame_is_the_curved_structure(built):
    ctl = built
    start = np.vstack(ca_blocks(ctl.win.structure, ctl.residues))
    assert shape_difference(drawn_sites(ctl, 0.0), start) < 0.05


def test_the_two_ends_are_not_the_same_shape(built):
    """The calibration for both of the above.

    A comparison that called the drawn endpoint 7WLU would be worthless if it
    also called the drawn *start* 7WLU. The gap it has to resolve is the
    conformational change itself.
    """
    ctl = built
    flat = flat_reference(ctl)
    assert shape_difference(drawn_sites(ctl, 0.0), flat) > 0.5 * CHANGE_RMSD
    assert shape_difference(drawn_sites(ctl, 1.0), flat) < 0.05


def test_a_field_applied_in_the_wrong_frame_lands_nowhere_near(built):
    """Pins the defect itself, in the units a user would have seen it in.

    The displacement field measured in the deposited frame, added to the
    canonically framed coordinates — exactly what the controller used to do,
    and it reproduces the failure to 0.01 A. The endpoint misses 7WLU by
    **36 A**, nearly twice the 19.7 A change it was interpolating: the motion
    was not merely inaccurate, it ran the wrong way.
    """
    ctl = built
    from piezo1.structure.frame import canonical_transform
    deposited = Structure.from_file(STRUCTURE_DIR / f"{PAIR[0]}.cif")
    rot = canonical_transform(deposited).rotation
    field = ctl.trajectory.frames[-1] - ctl.trajectory.frames[0]
    # Orthogonal rotation, so a field measured in the canonical frame is the
    # deposited-frame field right-multiplied by it.
    wrong = np.vstack(ca_blocks(ctl.win.structure, ctl.residues)) + field @ rot
    assert shape_difference(wrong, flat_reference(ctl)) == pytest.approx(36.0,
                                                                        abs=0.5)


def test_the_status_line_says_what_the_far_end_of_the_slider_is(window, app):
    """The endpoint claim has to follow the method, not the wording.

    `restrained` and `linear` land on the target's C-alpha positions exactly;
    `modal` is confined to the elastic-network subspace, captures 95% of the
    change on this pair and stops about 6 A short. One sentence covering both
    would be false for one of them — and it is the *modal* case a user would
    otherwise report as this same bug.
    """
    ctl = window.morph_controller
    ctl.build({"start": PAIR[0], "end": PAIR[1], "method": "restrained"})
    app.processEvents()
    text = window.status_label.text()
    assert f"ends on {PAIR[1]}'s C-alpha positions" in text
    assert "only c-alphas are interpolated" in text.lower()
    assert f"not the deposited {PAIR[1]}" in text

    ctl.build({"start": PAIR[0], "end": PAIR[1], "method": "modal"})
    app.processEvents()
    text = window.status_label.text()
    assert "stops short" in text
    assert f"ends on {PAIR[1]}" not in text
    captured = ctl.trajectory.meta["fraction_captured_by_modes"]
    assert 0.85 < captured < 1.0
    # And the shortfall is real, so the sentence is describing something.
    assert shape_difference(drawn_sites(ctl, 1.0), flat_reference(ctl)) > 1.0

    ctl.build({"start": PAIR[0], "end": PAIR[1], "method": "restrained"})
    app.processEvents()


# --------------------------------------------------------------------------
# What the C-alpha path does not describe has to be carried sensibly
# --------------------------------------------------------------------------

def test_lipids_are_carried_by_the_helix_they_touch(built):
    """Every atom without a site of its own must be tied to a near one.

    7WLT's lipids are numbered 2601-2609 and the last shared residue is 2546,
    so assignment by residue number tied all 1,407 of them to the C-terminal
    C-alpha — a median of 64.8 A away — and they travelled with the CTD tip.
    """
    ctl = built
    st = ctl.win.structure
    sites = ctl.trajectory.frames[0]
    in_basis = np.isin(st.res_seq, ctl.residues)
    loose = ~in_basis
    assert loose.sum() > 1000, "7WLT carries lipids outside the shared basis"

    tied = np.linalg.norm(st.xyz[loose] - sites[ctl._map[loose]], axis=1)
    assert np.median(tied) < 10.0
    assert tied.max() < 25.0

    # The calibration: the rule that was there put them an order of magnitude
    # further off, so this is a measurement rather than a tautology.
    per = len(ctl.residues)
    by_number = np.clip(np.searchsorted(ctl.residues, st.res_seq[loose]),
                        0, per - 1)
    old = np.linalg.norm(st.xyz[loose] - sites[by_number], axis=1)
    assert np.median(old) > 50.0


def test_a_residue_in_the_basis_keeps_its_own_site(built):
    """Whole residues move together; nothing is re-tied by proximity."""
    ctl = built
    st = ctl.win.structure
    per = len(ctl.residues)
    chains = [ch for ch in st.chains
              if (st.mask_ca() & (st.chain == ch)).sum() > 300][:3]
    for p, ch in enumerate(chains):
        sel = (st.chain == ch) & np.isin(st.res_seq, ctl.residues)
        expected = p * per + np.searchsorted(ctl.residues, st.res_seq[sel])
        assert np.array_equal(ctl._map[sel], expected)


# --------------------------------------------------------------------------
# A path belongs to the structure it was built on
# --------------------------------------------------------------------------

def test_loading_another_structure_drops_the_path(window, app, built):
    """The slider must not push 7WLT's motion onto a different entry.

    Every other overlay is cleared on load; the morph was not, and it is the
    one holding a base coordinate array of its own.
    """
    other = "8YEZ" if (STRUCTURE_DIR / "8YEZ.cif").exists() else None
    if other is None:
        pytest.skip("8YEZ.cif not downloaded — run python -m piezo1.io.fetch")
    assert window.morph_controller.trajectory is not None
    window.load_structure(other)
    app.processEvents()
    ctl = window.morph_controller
    assert ctl.trajectory is None
    assert ctl._base is None and ctl._map is None
    assert not window.physics_panel.morph_slider.isEnabled()
    assert not window.physics_panel.morph_play.isEnabled()
    # And a stale slider drag is then a no-op rather than a wrong picture.
    ctl.show_frame(1.0)
