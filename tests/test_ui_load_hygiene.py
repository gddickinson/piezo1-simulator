"""What survives a structure change: nothing that describes the old one.

The history of `load_structure`'s clear list is a history of omissions — the
morph was missing until Round 87, the full-length overlay until Round 90c,
and the micelle, the planar membrane, the ion stream and the potential
colouring until Round 91. Each looked exactly like the others when it was
found: a picture of the previous entry left drawn (or, for the ion stream,
*animating*) over the new one, which nothing on screen contradicts. And the
open-a-file path had a two-entry copy of the list, so every overlay it
lacked survived opening a file.

So the guard here is not the list — lists decay, that is the finding — but
an **equivalence**: loading B after using A must leave the same scene as
loading B fresh, whatever was drawn in between. Any overlay any future round
adds is covered the moment its controller is discovered (has ``show`` and
``clear``/``reset``), with a ratchet on the discovery count so the sweep
cannot quietly find nothing.
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
moderngl = pytest.importorskip("moderngl")

from PyQt6.QtWidgets import QApplication, QWidget  # noqa: E402

from piezo1.config import STRUCTURE_DIR, RenderSettings  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402

#: The conducting entry: the one whose ion stream actually runs, which is
#: what makes "the stream stops on load" a test rather than a tautology.
CONDUCTING = "11ZC"
OTHER = "8YEZ"


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
    for pdb in (CONDUCTING, OTHER):
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb} not downloaded; run python -m piezo1.io.fetch")
    from piezo1.render.scene import Scene
    from piezo1.ui.gl_widget import configure_surface_format
    from piezo1.ui.main_window import MainWindow

    configure_surface_format()
    win = MainWindow()
    win.resize(1200, 800)
    win.show()
    app.processEvents()
    try:
        ctx = moderngl.create_standalone_context(require=410)
    except Exception as exc:                           # pragma: no cover
        pytest.skip(f"no OpenGL 4.1 context available: {exc}")
    scene = Scene(ctx, RenderSettings(samples=1))
    scene.resize(1200, 800)
    win.viewport.scene = scene
    win._on_scene_ready(scene)
    app.processEvents()
    if win.structure is None:
        pytest.skip("no default structure could be loaded")
    yield win
    for controller in ("analysis", "physics", "overlay"):
        cleanup = getattr(getattr(win, controller, None), "cleanup", None)
        if cleanup is not None:
            cleanup()
    app.processEvents()


def discovered_overlays(win) -> list:
    """Every controller that draws something: has show, and clear or reset.

    Discovery rather than a list, because the lists are what decayed. QWidgets
    are excluded — every widget has ``show()``, and none of them is an
    overlay.
    """
    found = []
    for value in vars(win).values():
        if isinstance(value, QWidget):
            continue
        if callable(getattr(value, "show", None)) and (
                callable(getattr(value, "clear", None))
                or callable(getattr(value, "reset", None))):
            found.append(value)
    return found


def _switch_on_everything(win, app) -> None:
    for controller in discovered_overlays(win):
        try:
            controller.show(True)
        except Exception:
            # A controller may refuse (no modes yet, no pore profile) — that
            # is its documented behaviour, and refusals draw nothing.
            pass
        app.processEvents()


def test_the_discovery_finds_the_overlays(window):
    """The ratchet. A refactor that renamed show/clear would make the sweep
    silently cover nothing, and this is the assertion that notices."""
    found = discovered_overlays(window)
    assert len(found) >= 12, [type(c).__name__ for c in found]


def test_loading_after_every_overlay_equals_loading_fresh(window, app):
    win = window
    win.load_structure(OTHER)
    app.processEvents()
    baseline = sorted(win.viewport.scene.batches)

    win.load_structure(CONDUCTING)
    app.processEvents()
    plain = len(win.viewport.scene.batches)
    _switch_on_everything(win, app)
    assert len(win.viewport.scene.batches) > plain, (
        "nothing was drawn, so the equivalence below would be vacuous")

    win.load_structure(OTHER)
    app.processEvents()
    after = sorted(win.viewport.scene.batches)
    assert after == baseline, (
        "loading after overlays left a different scene than loading fresh — "
        "something survived that describes the previous entry")


def test_the_ion_stream_does_not_survive_a_load(window, app):
    """The worst single case: left running, the stream keeps animating the
    old entry's ions along the old entry's path over the new structure —
    over 8YEZ, the very entry the wetting model refuses to animate."""
    win = window
    win.load_structure(CONDUCTING)
    app.processEvents()
    win.ion_flux.show(True)
    app.processEvents()
    if not win.ion_flux._running:
        pytest.skip("the conducting entry did not animate; "
                    "the stop-on-load test would be vacuous")

    win.load_structure(OTHER)
    app.processEvents()
    assert win.ion_flux._running is False
    assert len(win.ion_flux._path) == 0
    win.viewport._on_tick()
    assert win.ion_flux._step not in win.viewport._animations, (
        "the animation callback outlived the structure it animates")


def test_opening_a_file_clears_like_loading_one(window, app):
    """The open-from-disk path had its own two-entry copy of the clear list;
    it now shares `load_structure`'s, and this holds it there."""
    win = window
    opened = Structure.from_file(STRUCTURE_DIR / f"{CONDUCTING}.cif")
    win.show_opened_structure(opened)
    app.processEvents()
    baseline = sorted(win.viewport.scene.batches)

    win.load_structure(OTHER)
    app.processEvents()
    _switch_on_everything(win, app)
    # Plant stale spliced-model state: an opened file is deposited
    # coordinates, and any banner state left from the previous entry would
    # put PART PREDICTED over a file it does not describe.
    win.full_length = object()

    win.show_opened_structure(
        Structure.from_file(STRUCTURE_DIR / f"{CONDUCTING}.cif"))
    app.processEvents()
    assert sorted(win.viewport.scene.batches) == baseline, (
        "opening a file left overlays of the previous entry drawn")
    assert win.full_length is None and win.assembly is None
