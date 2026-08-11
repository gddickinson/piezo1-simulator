"""Do the new overlays put pixels on the screen?

Every controller test in this project checks that a batch was *built* with the
right contents. All of those passed for the entire life of the renderer while
cylinders drew nothing at all — the batch existed, the instance count was
right, and the fragment shader discarded every fragment. Two independent bugs,
either sufficient on its own, and the only assertion that would have caught
either is this one.

So each overlay is driven through the controller's own draw path against a
real GL scene, and the lit pixels are counted. Each is also checked with the
overlay cleared, because a count that is high with *and* without the batch is
counting the background.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

moderngl = pytest.importorskip("moderngl")

from piezo1.config import RenderSettings, STRUCTURE_DIR  # noqa: E402

SIZE = (320, 240)


@pytest.fixture(scope="module")
def context():
    try:
        return moderngl.create_standalone_context(require=410)
    except Exception as exc:                              # pragma: no cover
        pytest.skip(f"no OpenGL 4.1 context available: {exc}")


@pytest.fixture
def scene(context):
    from piezo1.render.scene import Scene

    scene = Scene(context, RenderSettings(samples=1))
    scene.resize(*SIZE)
    return scene


def _lit(context, scene) -> int:
    fbo = context.simple_framebuffer(SIZE)
    fbo.use()
    fbo.clear(0.05, 0.05, 0.07, 1.0)
    scene.render()
    pixels = np.frombuffer(fbo.read(components=3), np.uint8)
    return int((pixels.reshape(-1, 3).astype(int).sum(axis=1) > 60).sum())


class _Window:
    """A window stub carrying a real scene, so drawing really happens."""

    def __init__(self, scene, structure=None, pore=None, pockets=(),
                 hydration=None, numbering="human"):
        self.structure = structure
        self.record = type("R", (), {"protein": "PIEZO1",
                                     "numbering_species": numbering})()
        self.modes = object()
        self.analysis = type("A", (), {
            "pore": pore, "hydration": hydration, "pockets": list(pockets),
            "compute_pore": lambda self: None,
            "compute_pockets": lambda self: None})()
        self.fusion = type("F", (), {"model": None})()
        self.viewport = type("V", (), {"scene": scene,
                                       "update": lambda self: None})()
        self.status = ""

    def _set_status(self, text):
        self.status = text


def _frame_on(scene, points, pad=8.0):
    points = np.asarray(points, dtype=np.float32)
    lo, hi = points.min(axis=0) - pad, points.max(axis=0) + pad
    scene.camera.frame(np.stack([lo, hi]).astype(np.float32))


@pytest.fixture(scope="module")
def framed_11zc():
    path = STRUCTURE_DIR / "11ZC.cif"
    if not path.exists():
        pytest.skip("11ZC.cif not downloaded — run python -m piezo1.io.fetch")
    from piezo1.core import Structure
    from piezo1.structure.frame import apply_frame, canonical_transform

    st = Structure.from_file(path)
    return apply_frame(st, canonical_transform(st))


# --------------------------------------------------------------------------

def test_the_pore_surface_is_visible(context, scene, framed_11zc, open_profile):
    from piezo1.ui.pore_controller import PoreSurfaceController

    controller = PoreSurfaceController(
        _Window(scene, framed_11zc, pore=open_profile))
    _frame_on(scene, open_profile.centers)
    empty = _lit(context, scene)
    controller._draw()
    drawn = _lit(context, scene)
    assert drawn > empty + 2000, (
        f"the pore drew {drawn - empty} extra lit pixels; the probe spheres "
        f"are not reaching the screen")

    controller.clear()
    assert _lit(context, scene) == empty


def test_the_pockets_are_visible(context, scene):
    """Synthetic clusters: this asks whether alpha spheres draw, not whether
    the detector works, which `test_ui_pockets` covers without a GPU."""
    from piezo1.analysis.pockets import Pocket
    from piezo1.ui.pocket_controller import PocketController

    rng = np.random.default_rng(3)
    pockets = [Pocket(index=i + 1,
                      centers=rng.normal(i * 25.0, 4.0, size=(60, 3)),
                      radii=np.full(60, 3.0),
                      meta={"volume": 100.0})
               for i in range(3)]
    controller = PocketController(_Window(scene, object(), pockets=pockets))
    _frame_on(scene, np.vstack([p.centers for p in pockets]))
    empty = _lit(context, scene)
    controller._draw()
    drawn = _lit(context, scene)
    assert drawn > empty + 2000, f"pockets drew only {drawn - empty} pixels"

    controller.clear()
    assert _lit(context, scene) == empty


def test_the_allosteric_route_draws_its_tube_and_not_only_its_nodes(context,
                                                                    scene):
    """The tube is the part that carries the correlation colouring, and it is
    thinner than the node markers — so "the nodes drew" is not evidence that
    the tube did. Drawn separately, and the tube must contribute on its own.
    """
    from piezo1.ui.path_controller import AllostericPathController

    coords = np.stack([np.linspace(-60, 60, 24),
                       np.zeros(24),
                       10.0 * np.sin(np.linspace(0, 3, 24))], axis=1)
    controller = AllostericPathController(_Window(scene, object()))
    controller.result = {"residues": list(range(600, 624)),
                         "sites": list(range(24)), "cost": 0.2,
                         "correlations": list(np.linspace(0.95, 0.999, 23)),
                         "coords": coords, "alternative_cost": 0.21,
                         "source_name": "THU4 (TM13-TM16)"}
    _frame_on(scene, coords)
    empty = _lit(context, scene)
    controller._draw()
    both = _lit(context, scene)
    assert both > empty + 1000, f"the route drew only {both - empty} pixels"

    # Now hide the node markers and check the tube alone still draws. If it
    # does not, the picture has been the beads all along and the colouring
    # that carries the weakest-link information was never visible.
    scene.set_visible(f"{controller._names[1]}", False)
    tube_only = _lit(context, scene)
    scene.set_visible(f"{controller._names[1]}", True)
    assert tube_only > empty + 500, (
        f"with the nodes hidden the tube drew {tube_only - empty} pixels — "
        f"the connecting tube is not reaching the screen")


def test_the_nanodomain_shells_are_visible_and_translucent(context, scene,
                                                           framed_11zc):
    """Translucency is the property that matters: an opaque shell around the
    protein hides everything the picture exists to put it in context with."""
    from piezo1.physics.nanodomain import Nanodomain
    from piezo1.ui.nanodomain_controller import NanodomainController

    window = _Window(scene, framed_11zc)
    controller = NanodomainController(window)
    controller.model = Nanodomain(current_A=2.4e-12, calcium_fraction=0.05,
                                  distance_m=4e-9)
    origin = framed_11zc.xyz.mean(axis=0)
    controller.shells = {3e-4: 40.0, 1e-4: 90.0}
    controller.occupancy_radii = {0.9: 1.2e3, 0.5: 3.7e3}
    controller.source_point = lambda: origin

    _frame_on(scene, framed_11zc.xyz)
    empty = _lit(context, scene)
    controller._draw()
    shells = _lit(context, scene)
    assert shells > empty + 2000, f"the shells drew {shells - empty} pixels"

    # The protein behind the shells must still be reachable: draw it too and
    # check the shells did not simply replace it.
    from piezo1.render.representations import MolecularView

    view = MolecularView(scene, framed_11zc, name="protein")
    view.rebuild()
    with_protein = _lit(context, scene)
    controller.clear()
    protein_only = _lit(context, scene)
    assert with_protein >= protein_only, \
        "adding translucent shells removed lit pixels from the protein"


# --------------------------------------------------------------------------
# Round 84b's Figure 4 overlays
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def framed_6b3r():
    path = STRUCTURE_DIR / "6B3R.cif"
    if not path.exists():
        pytest.skip("6B3R.cif not downloaded — run python -m piezo1.io.fetch")
    from piezo1.core import Structure
    from piezo1.structure.frame import apply_frame, canonical_transform

    st = Structure.from_file(path)
    return apply_frame(st, canonical_transform(st))


def test_the_micelle_envelope_is_visible(context, scene, framed_6b3r):
    """A translucent closed surface is the easiest thing to upload and not see.

    It encloses the protein, so a back-face cull would hide the near half and
    a front-face cull the far half; either leaves a batch with the right
    triangle count and a picture that is wrong.
    """
    from piezo1.ui.micelle_controller import MicelleController

    controller = MicelleController(
        _Window(scene, framed_6b3r, numbering="mouse"))
    _frame_on(scene, framed_6b3r.xyz[framed_6b3r.mask_ca()], pad=25.0)
    empty = _lit(context, scene)
    controller._build()
    assert controller.envelope is not None, controller.win.status
    drawn = _lit(context, scene)
    assert drawn > empty + 2000, (
        f"the micelle drew {drawn - empty} extra lit pixels")

    controller.clear()
    assert _lit(context, scene) == empty


def test_the_planar_membrane_is_visible(context, scene, framed_6b3r):
    from piezo1.ui.planar_membrane_controller import PlanarMembraneController

    controller = PlanarMembraneController(
        _Window(scene, framed_6b3r, numbering="mouse"))
    _frame_on(scene, framed_6b3r.xyz[framed_6b3r.mask_ca()], pad=25.0)
    empty = _lit(context, scene)
    controller._build()
    assert controller.comparison is not None, controller.win.status
    drawn = _lit(context, scene)
    assert drawn > empty + 2000, (
        f"the membrane planes drew {drawn - empty} extra lit pixels")

    controller.clear()
    assert _lit(context, scene) == empty


def test_fitting_the_planes_to_the_trimer_moves_them(context, scene,
                                                     framed_6b3r):
    """The control that makes Figure 4a a claim: the same construction on the
    assembly, where the paper says it does not work."""
    from piezo1.ui.planar_membrane_controller import PlanarMembraneController

    controller = PlanarMembraneController(
        _Window(scene, framed_6b3r, numbering="mouse"))
    controller._build()
    protomer = controller.status_line()
    controller.use_trimer(True)
    trimer = controller.status_line()
    assert protomer != trimer
    assert "all three protomers" in trimer
    controller.clear()
