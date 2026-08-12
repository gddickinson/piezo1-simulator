"""The particle animation, and the two ways it could lie.

Round 33 built the permeation physics and deferred the animation. The physics
is the easy part: a channel passes ~10^7 ions per second, so any watchable
stream is about a millionfold slow, and an unlabelled stream would misstate the
rate by six orders of magnitude.

The second way it could lie is worse. Every deposited human structure is closed,
and the permeation result is gated by the wetting verdict — so a trickle of ions
through a shut gate would contradict the project's own structural result while
looking like a demonstration of it.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.render.flux import (DISPLAY_PARTICLES_PER_SECOND,
                                ELEMENTARY_CHARGE, FluxTimebase, ion_rate,
                                timebase)


# ------------------------------------------------------------ the rate

def test_the_ion_rate_matches_the_definition_of_current():
    """Calibrated against the arithmetic, not against itself.

    One ampere is one coulomb per second, so 1 pA is 1e-12 / e ions per second
    for a monovalent ion. Checked at a round number before the real current is
    trusted through it.
    """
    assert ion_rate(1.0) == pytest.approx(1e-12 / ELEMENTARY_CHARGE, rel=1e-12)
    assert ion_rate(1.0) == pytest.approx(6.2415e6, rel=1e-4)


def test_a_divalent_ion_carries_the_same_current_with_half_the_ions():
    """Matters because PIEZO1 is calcium-permeable."""
    assert ion_rate(2.0, valence=2) == pytest.approx(ion_rate(2.0, valence=1) / 2)


def test_a_nonsensical_valence_is_refused():
    with pytest.raises(ValueError, match="valence"):
        ion_rate(1.0, valence=0)


def test_the_current_a_channel_actually_passes_needs_a_millionfold_slowdown():
    """The measured fact that makes the label non-optional."""
    base = timebase(1.65)
    assert base.ions_per_second == pytest.approx(1.03e7, rel=0.02)
    assert 1e5 < base.slowdown < 1e7, base.slowdown
    assert base.particles_per_second == DISPLAY_PARTICLES_PER_SECOND


# ----------------------------------------------------- the honest label

def test_the_statement_names_the_slowdown_and_the_particle_meaning():
    text = timebase(1.65).statement()
    assert "ions/s" in text
    assert "slower than real time" in text
    assert "1 particle = 1 ion" in text


def test_one_particle_is_one_ion_by_construction():
    """The alternative convention would make a viewer count the wrong thing."""
    assert timebase(1.65).ions_per_particle == 1.0


def test_the_slowdown_scales_with_the_current():
    """Twice the current, twice as slow — the label must track the physics."""
    slow_1 = timebase(1.0).slowdown
    slow_2 = timebase(2.0).slowdown
    assert slow_2 == pytest.approx(2 * slow_1, rel=1e-9)


# --------------------------------------------- a shut pore animates nothing

def test_a_non_conducting_pore_produces_no_particles():
    base = timebase(3.0, conducting=False, reason="hydrophobic gate")
    assert not base.conducting
    assert base.particles_per_second == 0.0
    assert base.ions_per_second == 0.0
    assert base.current_pA == 0.0, (
        "a non-conducting time base must not carry a current a caller could draw")


def test_a_shut_pore_states_why_instead_of_showing_a_rate():
    base = timebase(3.0, conducting=False, reason="sum d = 0.90 > 0.55")
    text = base.statement()
    assert "no conduction" in text
    assert "0.90" in text
    assert "slower than real time" not in text, (
        "a shut pore must not advertise a time base it is not using")


def test_zero_current_is_treated_as_not_conducting_even_if_claimed_open():
    """Defensive: a solver returning zero must not divide into a slowdown."""
    base = timebase(0.0, conducting=True)
    assert not base.conducting
    assert base.slowdown == 1.0


def test_a_default_reason_is_supplied_rather_than_left_blank():
    assert timebase(0.0, conducting=False).reason


# ------------------------------------------------------------ the controller

def test_the_controller_uses_the_shared_gated_computation():
    """The GUI must not reach its own conclusion about whether this conducts.

    This began as a test that read the controller's source for
    ``predict_wetting``. That passed while the units were wrong by 10^12, and
    then broke the moment the gating moved somewhere testable — which is the
    argument against source-inspection tests. It now asserts that the
    controller calls the shared function, and the *behaviour* of that function
    is checked below on real coordinates.
    """
    import inspect

    from piezo1.ui import ion_flux_controller

    source = inspect.getsource(ion_flux_controller.IonFluxController._prepare)
    assert "timebase_for_structure" in source, (
        "the controller must use the shared gated computation, not its own")


def test_the_controller_is_reachable_from_the_view_menu():
    from pathlib import Path

    menus = Path(__file__).resolve().parents[1] / "piezo1" / "ui" / "menus.py"
    text = menus.read_text()
    assert "ion_flux.show" in text
    assert "MILLIONFOLD" in text.upper(), (
        "the menu tooltip must warn about the time base before it is switched on")


def test_the_particles_are_not_coloured_like_any_element():
    """They represent a rate, not resolved atoms."""
    from piezo1.core.structure import ELEMENT_COLORS
    from piezo1.ui.ion_flux_controller import ION_COLOR

    for colour in ELEMENT_COLORS.values():
        assert not np.allclose(np.asarray(colour, float)[:3],
                               np.asarray(ION_COLOR, float), atol=0.02)


# ------------------------------- the real structures, on real coordinates

def _structure(pdb):
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure

    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb} not downloaded")
    return Structure.from_file(path)


def test_the_closed_human_structure_animates_nothing():
    """8YEZ is shut on both criteria; the animation must agree with the pipeline."""
    from piezo1.render.flux import timebase_for_structure

    base = timebase_for_structure(_structure("8YEZ"))
    assert not base.conducting
    assert base.particles_per_second == 0.0
    assert "non-conductive" in base.reason or "occluded" in base.reason


def test_the_open_like_structure_gives_a_watchable_rate():
    """11ZC conducts, and the units are the ones that bit on the first attempt.

    `PermeationResult.current` is in AMPERES. The controller's first version
    read a `current_pA` field that does not exist; had it existed and been
    misread, the rate would have been wrong by 10^12 and the animation would
    have looked plausible either way.
    """
    from piezo1.render.flux import timebase_for_structure

    base = timebase_for_structure(_structure("11ZC"))
    assert base.conducting
    assert 1.0 < base.current_pA < 10.0, (
        f"{base.current_pA} pA is not a single-channel current — check units")
    assert 1e6 < base.ions_per_second < 1e8
    assert 1e4 < base.slowdown < 1e8


def test_the_rate_declares_the_conductance_disagreement_it_inherits():
    """The solver records 41 pS against a published 25-30, reported not tuned.

    An animation calibrated to a number 1.5x too large, shown without saying
    so, is the confident-wrong-picture failure Round 50 audited for.
    """
    from piezo1.parameters import PARAMETERS
    from piezo1.render.flux import timebase_for_structure

    base = timebase_for_structure(_structure("11ZC"))
    assert base.caveat, "the inherited overestimate must be stated"
    assert "pS" in base.caveat
    published = PARAMETERS.value("permeation.published_conductance")
    assert f"{published:.0f}" in base.caveat
    assert base.caveat in base.statement()


def test_a_rate_that_matched_the_measurement_would_carry_no_caveat():
    """Calibration: the caveat must be able to be absent, or it says nothing."""
    from piezo1.render.flux import timebase

    assert timebase(1.65).caveat == ""


# ---------------------------------------- the refusal has to say which one

def test_the_refusal_names_the_constriction_and_the_gate_beside_it():
    """8IXO is the case that makes this non-optional.

    It is the intermediate-open S2472E structure of Liu et al. 2025, its gate
    is 3.52 A, its lining clears the Rao cutoff — and it is refused. A message
    saying only "sterically occluded" reads as a shut gate, which is the one
    thing that is not true of it.
    """
    from piezo1.render.flux import timebase_for_structure

    base = timebase_for_structure(_structure("8IXO"))
    assert not base.conducting
    assert "transmembrane gate" in base.reason, base.reason
    assert "cytosol" in base.reason, base.reason
    assert "2537" in base.reason, base.reason


def test_the_refusal_still_carries_the_verdict_it_came_from():
    """The location is added to the wetting summary, never in place of it."""
    from piezo1.analysis.hydration import load_grid, predict_wetting
    from piezo1.render.flux import timebase_for_structure
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    grid = load_grid()
    if not grid.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")

    st = _structure("7WLT")
    blocks, _ = protomer_blocks(st)
    profile = pore_profile(st, detect_c3_axis(blocks))
    verdict = predict_wetting(st, profile, grid=grid)
    assert verdict.summary() in timebase_for_structure(st, profile=profile).reason


# --------------------------------------------- the controller, against real GL

class _Viewport:
    """Records what the controller registers, and owns a real scene."""

    def __init__(self, scene):
        self.scene = scene
        self.hud = None
        self.animations = []

    def add_animation(self, callback):
        self.animations.append(callback)

    def update(self):
        pass


class _Window:
    def __init__(self, structure, scene):
        self.structure = structure
        self.viewport = _Viewport(scene)
        self.status = ""

    def _set_status(self, text):
        self.status = text


@pytest.fixture(scope="module")
def gl():
    moderngl = pytest.importorskip("moderngl")
    try:
        return moderngl.create_standalone_context(require=410)
    except Exception as exc:                              # pragma: no cover
        pytest.skip(f"no OpenGL 4.1 context available: {exc}")


@pytest.fixture
def flowing(gl):
    """The controller, switched on, on the one entry that conducts."""
    from piezo1.config import RenderSettings
    from piezo1.render.scene import Scene
    from piezo1.ui.ion_flux_controller import IonFluxController

    scene = Scene(gl, RenderSettings(samples=1))
    scene.resize(320, 240)
    window = _Window(_structure("11ZC"), scene)
    controller = IonFluxController(window)
    controller.show(True)
    if controller._timebase is None or not controller._timebase.conducting:
        pytest.skip("11ZC did not come back conducting")
    return controller, window, scene


def test_the_stream_survives_its_empty_first_frame(flowing):
    """The bug: the animation died on frame one and never drew an ion.

    12 particles enter per displayed second, so at 60 fps the first four frames
    hold nothing. The empty upload raised, `_on_tick` caught it and
    unregistered the callback, and the only conducting structure in the
    catalogue looked exactly like the 17 that are refused.
    """
    controller, window, _ = flowing
    assert window.viewport.animations, "the animation was never registered"

    assert len(controller._positions) == 0, (
        "the premise: the first frame really is empty")
    assert controller._step(1.0 / 60.0) is not False, (
        "the first frame unregistered the animation")

    for _ in range(90):
        assert controller._step(1.0 / 60.0) is not False
    assert len(controller._positions) > 0, "no ions ever entered the pore"


def test_the_ions_reach_the_screen_and_leave_it_again(gl, flowing):
    """Counted in pixels, because a built batch proves nothing.

    Every earlier test of this controller passed while it drew no ion at all.
    """
    import numpy as np

    controller, _, scene = flowing
    scene.camera.frame(controller._path.astype(np.float32))

    def lit():
        fbo = gl.simple_framebuffer((320, 240))
        fbo.use()
        fbo.clear(0.05, 0.05, 0.07, 1.0)
        scene.render()
        pixels = np.frombuffer(fbo.read(components=3), np.uint8)
        return int((pixels.reshape(-1, 3).astype(int).sum(axis=1) > 60).sum())

    empty = lit()
    for _ in range(90):
        controller._step(1.0 / 60.0)
    assert lit() > empty + 100, "the ions are not reaching the screen"

    controller.clear()
    assert lit() == empty


def test_the_stream_follows_the_measured_pore_and_not_the_symmetry_axis(flowing):
    """The probe centre is leashed to the axis, not pinned to it.

    On 11ZC the fitted centre sits a median 0.56 A off the axis and up to the
    full 8 A leash, and at 11 of 125 heights the axis line falls *outside* the
    sphere fitted there — so ions drawn on the straight axis cross the wall of
    the pore that was measured. The calibration is that assertion, on this
    entry's own numbers: if the two ever coincide, this test should say so
    rather than passing vacuously.
    """
    import numpy as np

    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    controller, window, _ = flowing
    blocks, _residues = protomer_blocks(window.structure)
    axis = detect_c3_axis(blocks)
    profile = pore_profile(window.structure, axis)

    direction = axis.direction / np.linalg.norm(axis.direction)
    straight = axis.point[None, :] + profile.z[:, None] * direction[None, :]
    offset = np.linalg.norm(profile.centers - straight, axis=1)
    assert (offset > profile.radius).sum() >= 5, (
        "the axis no longer leaves the fitted probe anywhere, so this test "
        "cannot distinguish the two routes")

    # The stream visits every measured probe centre, in order.
    assert np.allclose(controller.points_at(controller._arc), controller._path)
    assert sorted(map(tuple, np.round(controller._path, 6))) == sorted(
        map(tuple, np.round(profile.centers, 6)))


def test_the_stream_runs_towards_the_cytosolic_end(flowing):
    """Which way is out is measured, not left to the sign of the axis.

    `detect_c3_axis` fixes a line, not a direction, so without this half the
    catalogue would show an inward current flowing out of the cell — a picture
    that looks entirely reasonable and is backwards.
    """
    import numpy as np

    from piezo1.physics.pore_charge import cytosolic_end
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    controller, window, _ = flowing
    blocks, _residues = protomer_blocks(window.structure)
    axis = detect_c3_axis(blocks)
    profile = pore_profile(window.structure, axis)

    entry, exit_ = controller.points_at([0.0, controller.length])
    cytosolic_z = (profile.z.min() if cytosolic_end(window.structure, axis) == 0
                   else profile.z.max())
    assert (abs(axis.project(exit_[None, :])[0] - cytosolic_z)
            < abs(axis.project(entry[None, :])[0] - cytosolic_z)), (
        "the stream ends further from the cytosolic mouth than it began")
