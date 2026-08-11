"""Selecting atoms: the feedback loop between a click and the Measure panel.

`MeasurementSet` had a thorough suite and the panel that drives it had none, so
three gaps survived in the layer between them — all reported by a user who
could not work out how to select anything:

* a click marked nothing on the model, so there was no way to tell which atom
  had been hit, or that the click had registered at all;
* a pick appeared in the panel only when the *last* atom of a measurement was
  clicked, so selecting one atom of two left the Selection table empty;
* nothing said that picking has to be armed first, and the hint under the
  button described the goal ("Pick 2 atoms for a distance") in words that read
  as though clicking would already work.

These are the tests for that layer. They use the real widgets on the offscreen
platform, because the bug was in the wiring rather than in either end of it.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402


class _Batch:
    def upload(self, *args, **kwargs):
        self.args = args


class _Scene:
    def __init__(self):
        self.batches = {}

    def spheres(self, name):
        return self.batches.setdefault(name, _Batch())

    cylinders = spheres

    def mesh(self, name, **kwargs):
        return self.batches.setdefault(name, _Batch())

    def remove(self, name):
        self.batches.pop(name, None)

    def get(self, name):
        return self.batches.get(name)


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance()
    if app is None:
        try:
            app = QApplication([])
        except Exception as exc:                       # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return app


@pytest.fixture(scope="module")
def structure():
    path = STRUCTURE_DIR / "8YEZ.cif"
    if not path.exists():
        pytest.skip("8YEZ not downloaded; run python -m piezo1.io.fetch")
    return Structure.from_file(path)


@pytest.fixture
def win(qapp, structure):
    """A real MainWindow with a recording scene in place of a GL context."""
    from piezo1.render.representations import ColorBy, MolecularView, Style
    from piezo1.ui.gl_widget import configure_surface_format
    from piezo1.ui.main_window import MainWindow

    configure_surface_format()
    window = MainWindow()
    window.structure = structure
    window.viewport.scene = _Scene()
    window.view = MolecularView(window.viewport.scene, structure, name="8YEZ",
                                style=Style.CARTOON, color_by=ColorBy.CHAIN)
    return window


@pytest.fixture
def atoms(structure):
    return np.flatnonzero(structure.mask_ca())


# ------------------------------------------------- a click has to show itself

def test_an_unarmed_click_marks_the_residue_on_the_model(win, structure, atoms):
    """The missing feedback. Without it a click changed only the status bar."""
    win._on_pick(int(atoms[0]))

    assert win.view.highlight is not None, "clicking highlighted nothing"
    assert win.view.highlight.sum() > 0
    assert "8YEZ:selection" in win.viewport.scene.batches, \
        "the selection is not drawn, so it cannot be seen"


def test_a_click_selects_the_copy_clicked_not_all_three(win, structure, atoms):
    """A residue *number* means three copies in a trimer; a click means one.

    Annotations select by number and want all three. A click points at one
    protomer, and highlighting the other two would claim the user meant them.
    """
    index = int(atoms[0])
    win._on_pick(index)

    assert win.view.highlight is not None, "clicking highlighted nothing"
    chains = set(structure.chain[win.view.highlight].tolist())
    assert chains == {str(structure.chain[index])}
    assert len(structure.chains) > 1, "8YEZ is a trimer; the check needs one"


def test_the_status_bar_still_identifies_the_residue(win, atoms):
    """The highlight is added to inspection, not in place of it."""
    win._on_pick(int(atoms[0]))
    assert "GLU570" in win.status_label.text()
    assert "chain A" in win.status_label.text()


def test_an_unarmed_click_measures_nothing(win, atoms):
    """Clicks must not be silently consumed by a tool that is switched off."""
    win._on_pick(int(atoms[0]))
    win._on_pick(int(atoms[1]))
    assert win.measure_panel.set.measurements == []
    assert win.measure_panel.set.pending == []


def test_nothing_under_the_cursor_is_reported_rather_than_ignored(win):
    win._on_pick(-1)
    assert "nothing under the cursor" in win.status_label.text()


# ------------------------------------------------------- finding the feature

def test_the_hint_names_the_step_not_the_goal(win):
    """The old text implied clicking already worked."""
    hint = win.measure_panel.hint.text()
    assert "Start picking" in hint
    assert win.measure_panel.arm_button.text() == "Start picking"

    win.measure_panel.arm_button.setChecked(True)
    armed = win.measure_panel.hint.text()
    assert "Click atoms" in armed and "0 of 2" in armed
    assert win.measure_panel.arm_button.text() != "Start picking", \
        "the button must say what it is doing, not what it would do"


def test_the_first_few_clicks_say_where_the_measuring_tool_is(win, atoms):
    """Discoverability, bounded so it does not become a permanent nag.

    The status bar's job is to identify the residue; crowding that out forever
    would trade one complaint for another.
    """
    from piezo1.ui.main_window import PICK_HINTS

    seen = []
    for i in range(PICK_HINTS + 2):
        win._on_pick(int(atoms[i]))
        seen.append("Start picking" in win.status_label.text())

    assert all(seen[:PICK_HINTS]), "the hint never appeared"
    assert not any(seen[PICK_HINTS:]), "the hint never stopped"


def test_the_hint_is_not_shown_once_picking_is_armed(win, atoms):
    """It would be advice to do the thing already being done."""
    win.measure_panel.arm_button.setChecked(True)
    win._on_pick(int(atoms[0]))
    assert "Start picking" not in win.status_label.text()


# -------------------------------------------- the selection has to be visible

def test_a_single_pick_appears_in_the_table_immediately(win, atoms):
    """The reported bug: the panel stayed empty until the last atom."""
    panel = win.measure_panel
    panel.arm_button.setChecked(True)
    assert panel.table.rowCount() == 0

    win._on_pick(int(atoms[0]))

    assert panel.table.rowCount() == 1, \
        "one atom of two picked and the Selection table is still empty"
    assert "GLU570" in panel.table.item(0, 0).text()


def test_a_pending_row_carries_no_value_and_says_what_is_missing(win, atoms):
    """It must not be mistaken for a result."""
    from piezo1.ui.panels.measure_panel import PENDING_COLOR

    panel = win.measure_panel
    panel.arm_button.setChecked(True)
    win._on_pick(int(atoms[0]))

    assert panel.table.item(0, 1).text() == "…"
    assert "1 more" in panel.table.item(0, 2).text()
    assert panel.table.item(0, 0).foreground().color() == PENDING_COLOR


def test_completing_a_measurement_replaces_the_pending_row(win, atoms):
    panel = win.measure_panel
    panel.arm_button.setChecked(True)
    win._on_pick(int(atoms[0]))
    win._on_pick(int(atoms[1]))

    assert panel.table.rowCount() == 1
    assert panel.table.item(0, 1).text() == "3.79"
    assert panel.table.item(0, 2).text() == "A"
    assert "0 of 2" in panel.hint.text()


def test_an_angle_shows_two_picks_pending_then_resolves(win, atoms):
    panel = win.measure_panel
    panel.kind_combo.setCurrentText("angle")
    panel.arm_button.setChecked(True)

    win._on_pick(int(atoms[0]))
    win._on_pick(int(atoms[1]))
    assert panel.table.rowCount() == 1
    assert "2 more" not in panel.table.item(0, 2).text()
    assert "1 more for an angle" == panel.table.item(0, 2).text()

    win._on_pick(int(atoms[2]))
    assert panel.table.item(0, 1).text() != "…"
    assert panel.table.item(0, 2).text() == "deg"


def test_the_picked_atoms_are_drawn_in_the_viewport(win, atoms):
    panel = win.measure_panel
    panel.arm_button.setChecked(True)
    win._on_pick(int(atoms[0]))

    batch = win.viewport.scene.get("measure:picks")
    assert batch is not None, "a pick is listed but not marked on the model"
    assert len(batch.args[0]) == 1


def test_disarming_abandons_the_partial_selection(win, atoms):
    panel = win.measure_panel
    panel.arm_button.setChecked(True)
    win._on_pick(int(atoms[0]))
    assert panel.table.rowCount() == 1

    panel.arm_button.setChecked(False)
    assert panel.set.pending == []
    assert panel.table.rowCount() == 0, "the abandoned pick is still listed"


def test_deleting_the_pending_row_does_not_delete_a_measurement(win, atoms):
    """The pending row sits past the end of `measurements`; a naive index
    would have removed the wrong thing, or raised."""
    panel = win.measure_panel
    panel.arm_button.setChecked(True)
    win._on_pick(int(atoms[0]))
    win._on_pick(int(atoms[1]))          # one completed measurement
    win._on_pick(int(atoms[2]))          # and one pending pick
    assert panel.table.rowCount() == 2

    panel.table.selectRow(1)
    panel._delete()

    assert len(panel.set.measurements) == 1, "the completed measurement went"
    assert panel.set.pending == []
    assert panel.table.rowCount() == 1
