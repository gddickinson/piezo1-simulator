"""Showing one structure, and showing several on purpose.

The fault these cover was invisible to any test of the pieces in isolation:
`MolecularView.clear()` worked, `Scene.remove()` worked, and the structure was
still on screen — because `overlay.clear()` rebuilt the old view in between.

Two layers, because the offscreen Qt platform the rest of the UI suite runs on
cannot create an OpenGL context:

* :func:`test_clearing_the_overlay_cannot_resurrect_a_cleared_view` needs no GL
  at all and encodes the fault directly, so it runs everywhere.
* the rest drive the real window and **skip without a display**. They are the
  ones that prove the whole path, and they are why the ordering is what it is.
"""

from __future__ import annotations

import pytest

pytest.importorskip("PyQt6")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.ui.companions import COMPANION_COLORS  # noqa: E402


# ------------------------------------------------------- the fault, without GL

class _RecordingScene:
    """Just enough Scene to record what a view puts in and takes out."""

    def __init__(self):
        self.batches: dict[str, str] = {}

    def remove(self, name):
        self.batches.pop(name, None)

    def _add(self, name):
        self.batches[name] = "batch"
        return self

    mesh = spheres = cylinders = lambda self, name, **kw: self._add(name)


class _RecordingView:
    """A stand-in for MolecularView with the same clear/rebuild contract."""

    def __init__(self, scene, name):
        self.scene, self.name = scene, name
        self.values = None
        self.color_by = None

    def rebuild(self):
        self.scene.mesh(f"{self.name}:ribbon")

    def clear(self):
        for suffix in ("ribbon", "atoms", "bonds", "ligands"):
            self.scene.remove(f"{self.name}:{suffix}")


class _FakeWindow:
    def __init__(self):
        self.viewport = type("V", (), {"scene": _RecordingScene(),
                                       "update": lambda self: None})()
        self.view = _RecordingView(self.viewport.scene, "8YEZ")
        self.view.rebuild()

    def _current_color(self):
        return None


def test_clearing_the_overlay_cannot_resurrect_a_cleared_view():
    """`OverlayController.clear()` ends by rebuilding the primary view.

    That is right on its own — it undoes per-residue deviation colouring — but
    it means the caller must not have cleared that view first. `load_structure`
    did exactly that, so the rebuild put the old structure's batches straight
    back and it stayed on screen for good. Setting ``view = None`` before
    anything can rebuild is what makes the order safe, and this pins it.
    """
    from piezo1.ui.overlay_controller import OverlayController

    def clear_overlay(win):
        OverlayController.clear(type("C", (), {
            "cleanup": lambda self: None, "win": win,
            "view": None, "result": None})())

    # With a live view, rebuilding is the intended behaviour: it is how the
    # deviation colouring gets undone.
    win = _FakeWindow()
    win.view.clear()
    clear_overlay(win)
    assert "8YEZ:ribbon" in win.viewport.scene.batches, (
        "the rebuild is deliberate; if it stops happening, deviation colouring "
        "will no longer be undone and this test is guarding nothing")

    # So the caller has to hand over a window with no view — which is what
    # load_structure now does before it swaps structures.
    win = _FakeWindow()
    win.view.clear()
    win.view = None
    clear_overlay(win)
    assert win.viewport.scene.batches == {}, (
        "clearing the overlay resurrected a view that had already been "
        "cleared; load_structure must null the view before calling it")


def _require(*pdbs):
    for pdb in pdbs:
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb} not downloaded")


@pytest.fixture(scope="module")
def window(qt_app):
    from piezo1.ui.main_window import MainWindow
    win = MainWindow()
    win.show()
    qt_app.processEvents()
    if win.viewport.scene is None:
        pytest.skip("no GL scene available")
    yield win
    win.close()


def _batches(win):
    return sorted(win.viewport.scene.batches)


def _structures_drawn(win):
    """Distinct structures with batches in the scene, however they are named."""
    names = set()
    for key in win.viewport.scene.batches:
        head = key.rsplit(":", 1)[0]
        names.add(head[len("extra:"):] if head.startswith("extra:") else head)
    return names


# ------------------------------------------------------- the default: one

def test_loading_replaces_rather_than_accumulates(window, qt_app):
    """The regression.

    `overlay.clear()` ends by rebuilding the primary view to undo deviation
    colouring. `load_structure` called it *after* clearing the old view, so the
    rebuild put the old batches straight back and the previous structure stayed
    on screen permanently. It went unnoticed only because deposited frames sat
    100 A apart; once every structure is framed canonically the two superimpose.
    """
    _require("8YEZ", "7WLT")
    window.set_multi_structure(False)

    for pdb in ("8YEZ", "7WLT", "8YEZ"):
        window.load_structure(pdb)
        qt_app.processEvents()
        assert _structures_drawn(window) == {pdb}, (
            f"after loading {pdb} the scene holds {_batches(window)}")
        assert window.displayed_structures() == [pdb]


def test_the_indicator_stays_hidden_for_a_single_structure(window, qt_app):
    _require("8YEZ")
    window.set_multi_structure(False)
    window.load_structure("8YEZ")
    qt_app.processEvents()
    assert not window.structure_panel.displayed_label.isVisible()


# --------------------------------------------------- the opt-in: several

def test_multi_structure_keeps_the_previous_one(window, qt_app):
    _require("8YEZ", "7WLT", "8ZU8")
    window.set_multi_structure(False)
    window.load_structure("8YEZ")
    qt_app.processEvents()

    window.set_multi_structure(True)
    for pdb in ("7WLT", "8ZU8"):
        window.load_structure(pdb)
        qt_app.processEvents()

    # Primary first, then companions oldest-first in the order they were added.
    assert window.displayed_structures() == ["8ZU8", "8YEZ", "7WLT"]
    assert _structures_drawn(window) == {"8YEZ", "7WLT", "8ZU8"}

    # The most recently loaded one is the primary, and it is the only one the
    # rest of the application knows about.
    assert window.record.pdb == "8ZU8"
    assert window.view.name == "8ZU8"


def test_companions_are_given_distinct_colours(window, qt_app):
    _require("8YEZ", "7WLT", "8ZU8")
    window.set_multi_structure(True)
    window.clear_companions()
    window.load_structure("8YEZ")
    qt_app.processEvents()
    window.add_companion("7WLT")
    window.add_companion("8ZU8")
    qt_app.processEvents()

    colors = [c.color for c in window._companions().values()]
    assert len(set(colors)) == len(colors), "companions must be distinguishable"
    assert all(c in COMPANION_COLORS for c in colors)


def test_the_indicator_names_every_displayed_structure(window, qt_app):
    _require("8YEZ", "7WLT")
    window.set_multi_structure(True)
    window.clear_companions()
    window.load_structure("8YEZ")
    qt_app.processEvents()
    window.add_companion("7WLT")
    qt_app.processEvents()

    label = window.structure_panel.displayed_label
    assert label.isVisible()
    assert "8YEZ" in label.text() and "7WLT" in label.text()
    assert "primary" in label.text()


def test_turning_the_option_off_drops_the_extras(window, qt_app):
    """A setting saying "one structure" while three are drawn is worse than either."""
    _require("8YEZ", "7WLT")
    window.set_multi_structure(True)
    window.clear_companions()
    window.load_structure("8YEZ")
    qt_app.processEvents()
    window.add_companion("7WLT")
    qt_app.processEvents()
    assert len(window.displayed_structures()) == 2

    window.set_multi_structure(False)
    qt_app.processEvents()
    assert window.displayed_structures() == ["8YEZ"]
    assert _structures_drawn(window) == {"8YEZ"}


def test_a_structure_is_never_displayed_twice(window, qt_app):
    _require("8YEZ", "7WLT")
    window.set_multi_structure(True)
    window.clear_companions()
    window.load_structure("8YEZ")
    qt_app.processEvents()

    window.add_companion("8YEZ")             # already primary
    window.add_companion("7WLT")
    window.add_companion("7WLT")             # already a companion
    qt_app.processEvents()

    shown = window.displayed_structures()
    assert shown == ["8YEZ", "7WLT"]
    assert len(shown) == len(set(shown))


def test_removing_one_companion_leaves_the_others(window, qt_app):
    _require("8YEZ", "7WLT", "8ZU8")
    window.set_multi_structure(True)
    window.clear_companions()
    window.load_structure("8YEZ")
    qt_app.processEvents()
    window.add_companion("7WLT")
    window.add_companion("8ZU8")
    qt_app.processEvents()

    window.remove_companion("7WLT")
    qt_app.processEvents()
    assert window.displayed_structures() == ["8YEZ", "8ZU8"]
    assert _structures_drawn(window) == {"8YEZ", "8ZU8"}


def test_what_is_drawn_does_not_change_what_is_analysed(window, qt_app):
    """The standing rule, applied to a new way of putting things on screen.

    Hiding a lipid cannot change a pore profile, and neither can drawing a
    second structure. Analyses run on the primary, whatever else is displayed.
    """
    _require("8YEZ", "7WLT")
    from piezo1.structure.geometry import measure_dome
    from piezo1.structure.protomers import protomer_blocks

    def dome(win):
        blocks, _ = protomer_blocks(win.structure)
        import numpy as np
        return measure_dome(blocks, np.concatenate(blocks)).radius_of_curvature

    window.set_multi_structure(False)
    window.load_structure("8YEZ")
    qt_app.processEvents()
    alone = dome(window)

    window.set_multi_structure(True)
    window.add_companion("7WLT")
    qt_app.processEvents()
    assert window.record.pdb == "8YEZ"
    assert dome(window) == pytest.approx(alone, rel=1e-12)

    window.set_multi_structure(False)
    qt_app.processEvents()
