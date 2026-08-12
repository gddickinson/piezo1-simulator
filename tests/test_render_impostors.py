"""Do the impostors actually put pixels on the screen?

Every other renderer test in this project checks that batches are *created*.
None checked that anything is *drawn*, and cylinders were invisible for the
whole life of the project because of it: ball-and-stick drew balls only, the
HaloTag seam drew nothing, and the batch existed with the right instance count
the entire time. Two independent bugs, either of which was sufficient.

So these render to a framebuffer and count lit pixels. That is the only
assertion that would have caught it.
"""

from __future__ import annotations

import numpy as np
import pytest

moderngl = pytest.importorskip("moderngl")

from piezo1.config import RenderSettings  # noqa: E402


@pytest.fixture(scope="module")
def context():
    try:
        return moderngl.create_standalone_context(require=410)
    except Exception as exc:                              # pragma: no cover
        pytest.skip(f"no OpenGL 4.1 context available: {exc}")


def _scene(context, size=(320, 240)):
    from piezo1.render.scene import Scene

    scene = Scene(context, RenderSettings(samples=1))
    scene.resize(*size)
    return scene


def _lit(context, scene, size=(320, 240)) -> int:
    fbo = context.simple_framebuffer(size)
    fbo.use()
    fbo.clear(0.05, 0.05, 0.07, 1.0)
    scene.render()
    pixels = np.frombuffer(fbo.read(components=3), np.uint8)
    return int((pixels.reshape(-1, 3).astype(int).sum(axis=1) > 60).sum())


def _one_cylinder(scene):
    batch = scene.cylinders("bond")
    batch.upload(np.array([[-10.0, 0.0, 0.0]], np.float32),
                 np.array([[10.0, 0.0, 0.0]], np.float32),
                 np.array([3.0], np.float32),
                 np.array([[1.0, 0.2, 0.2]], np.float32))
    scene.camera.frame(np.array([[-12.0, -4.0, -4.0],
                                 [12.0, 4.0, 4.0]], np.float32))
    return batch


def test_a_cylinder_is_visible(context):
    """The whole regression, in one assertion."""
    scene = _scene(context)
    _one_cylinder(scene)
    assert _lit(context, scene) > 1000, (
        "a cylinder spanning most of the view drew nothing; impostor "
        "cylinders are invisible again")


def test_cylinder_batches_declare_that_they_must_not_be_culled(context):
    """The first of the two bugs, pinned as the property that fixes it.

    The bounding quad is oriented by the cylinder's own axis, so whether it
    faces the camera depends on the direction of the bond. With back-face
    culling on it is discarded before the fragment stage runs.
    """
    from piezo1.render.primitives import CylinderBatch, SphereBatch

    assert CylinderBatch.cull is False
    # Spheres build a screen-aligned quad, so their winding never flips and
    # they are unaffected — which is why only cylinders were invisible.
    assert getattr(SphereBatch, "cull", True) is True

    scene = _scene(context)
    _one_cylinder(scene)
    lit = _lit(context, scene)
    context.enable(moderngl.CULL_FACE)
    scene.batches["bond"].cull = True                 # put the bug back
    assert _lit(context, scene) == 0
    scene.batches["bond"].cull = False
    assert _lit(context, scene) == lit


def test_the_ray_cylinder_intersection_returns_a_hit_in_front_of_the_eye():
    """The second bug, checked in arithmetic rather than in pixels.

    The shader built the perpendicular offset from ``-oc`` instead of ``oc``,
    which negates B and therefore negates both roots — so the near hit came
    out behind the eye and was discarded. Reproduced here because a shader
    cannot be unit-tested and this is the line that mattered.
    """
    base = np.array([0.0, 0.0, -50.0])
    end = np.array([1.5, 0.0, -50.0])
    radius = 0.20
    ray = (0.5 * (base + end))
    ray = ray / np.linalg.norm(ray)
    axis = end - base
    direction = axis / np.linalg.norm(axis)
    perpendicular_ray = ray - direction * np.dot(ray, direction)

    def near_root(offset):
        a = np.dot(perpendicular_ray, perpendicular_ray)
        b = 2.0 * np.dot(perpendicular_ray, offset)
        c = np.dot(offset, offset) - radius ** 2
        disc = b * b - 4 * a * c
        return (-b - np.sqrt(disc)) / (2 * a) if disc >= 0 else float("nan")

    oc = -base
    correct = oc - direction * np.dot(oc, direction)
    broken = -oc - direction * np.dot(-oc, direction)
    assert near_root(correct) > 0.0
    assert near_root(broken) < 0.0
    assert near_root(correct) == pytest.approx(-near_root(broken), rel=0.02)


def test_ball_and_stick_draws_more_than_balls(context, structure_by_id):
    """The user-visible symptom, asserted against the style beside it."""
    from piezo1.render.representations import ColorBy, MolecularView, Style

    structure = structure_by_id("4RAX")
    if structure is None:
        pytest.skip("4RAX not downloaded — run python -m piezo1.io.fetch")

    counts = {}
    for style in (Style.BALLS, Style.BALL_AND_STICK, Style.STICKS):
        scene = _scene(context)
        view = MolecularView(scene, structure, name="m")
        view.style, view.color_by = style, ColorBy.ELEMENT
        view.rebuild()
        scene.camera.frame(structure.xyz)
        counts[style.value] = _lit(context, scene)

    assert counts["ball_and_stick"] > 1.05 * counts["balls"], (
        f"ball-and-stick lit {counts['ball_and_stick']} pixels against "
        f"{counts['balls']} for balls alone — the sticks are missing again")
    assert counts["sticks"] > 0.5 * counts["balls"]


def test_balls_and_sticks_are_separate_styles():
    """Both were asked for, and neither is the other."""
    from piezo1.render.representations import Style
    from piezo1.ui.panels.structure_panel import STYLE_LABELS

    offered = {value for _label, value in STYLE_LABELS}
    assert {Style.BALLS, Style.STICKS, Style.BALL_AND_STICK} <= offered


# ------------------------------------------------- an empty batch draws nothing

def test_an_empty_upload_does_not_raise(context):
    """Round 84c's regression, and the reason the ion stream never started.

    ``ctx.buffer`` refuses a zero-length payload — "the buffer cannot be
    empty". An animation whose first frame has no particles is the normal
    case, not an error, and the exception escaped the frame callback into
    ``ViewportWidget._on_tick``, which responds by *unregistering the
    animation*. So the one conducting structure in the catalogue died on frame
    one and looked exactly like the 17 that are refused.
    """
    scene = _scene(context)
    empty3 = np.zeros((0, 3), np.float32)

    scene.spheres("s").upload(empty3, np.zeros(0, np.float32), empty3)
    scene.cylinders("c").upload(empty3, empty3, np.zeros(0, np.float32), empty3)
    scene.mesh("m").upload(empty3, empty3, empty3, np.zeros((0, 3), np.int32))


def test_an_empty_batch_lights_no_pixels_and_then_a_full_one_does(context):
    """Calibration: the empty case must be *empty*, not merely quiet.

    A guard that swallowed the upload entirely would also pass the test above,
    so the same batch is filled afterwards and must reach the screen — and
    emptied again, and must leave it.
    """
    scene = _scene(context)
    batch = scene.spheres("ions")
    scene.camera.frame(np.array([[-12.0, -12.0, -12.0],
                                 [12.0, 12.0, 12.0]], np.float32))
    empty3 = np.zeros((0, 3), np.float32)

    batch.upload(empty3, np.zeros(0, np.float32), empty3)
    blank = _lit(context, scene)
    assert batch.count == 0

    batch.upload(np.zeros((1, 3), np.float32), np.array([6.0], np.float32),
                 np.array([[0.35, 0.85, 1.0]], np.float32))
    assert _lit(context, scene) > blank + 500, (
        "the guard swallowed a non-empty upload as well")

    batch.upload(empty3, np.zeros(0, np.float32), empty3)
    assert _lit(context, scene) == blank
