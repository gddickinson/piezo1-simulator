"""The right-click menu, and the two ways it could go wrong.

It could **do things twice**: a Representation entry that called `_set_style`
directly would change the model while the Model panel went on reporting the old
value, and the two would drift apart the first time anyone used the menu. So
the appearance entries drive the panel's own combo boxes, and that is tested by
watching the combo move rather than by reading the code.

And it could **act by being opened**: identifying the residue under the cursor
so the entries can name it must not count as selecting it, or dismissing the
menu would leave the model changed.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
moderngl = pytest.importorskip("moderngl")

from PyQt6.QtCore import Qt, QPointF  # noqa: E402
from PyQt6.QtGui import QMouseEvent  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.config import RenderSettings, STRUCTURE_DIR  # noqa: E402
from piezo1.ui.context_menu import build_context_menu  # noqa: E402


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
    if not (STRUCTURE_DIR / "8YEZ.cif").exists():
        pytest.skip("no structures downloaded; run python -m piezo1.io.fetch")
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
    return win


@pytest.fixture
def atom(window):
    return int(np.flatnonzero(window.structure.mask_ca())[0])


def labels(menu) -> list[str]:
    return [a.text() for a in menu.actions() if not a.isSeparator()]


def find(menu, fragment):
    for action in menu.actions():
        if fragment in action.text():
            return action
    raise AssertionError(f"no entry containing {fragment!r} in {labels(menu)}")


# ------------------------------------------------------- what is in the menu

def test_the_menu_names_what_was_clicked(window, atom):
    menu = build_context_menu(window, atom)
    caption = menu.actions()[0]
    assert not caption.isEnabled(), "the caption must not be clickable"
    st = window.structure
    assert str(st.res_name[atom]) in caption.text()
    assert f"chain {st.chain[atom]}" in caption.text()


def test_clicking_empty_space_still_gives_a_usable_menu(window):
    """A control that sometimes does nothing reads as broken."""
    menu = build_context_menu(window, -1)
    assert labels(menu), "right-clicking the background gave an empty menu"
    find(menu, "Representation")
    find(menu, "Reset the view")
    assert not any("Select" in t for t in labels(menu)), \
        "residue entries appeared with no residue under the cursor"


def test_the_residue_entries_appear_only_with_a_residue(window, atom):
    on_atom = labels(build_context_menu(window, atom))
    assert any("Add to measurement" in t for t in on_atom)
    assert any("Centre the view here" in t for t in on_atom)
    assert any("Copy" in t for t in on_atom)


def test_a_variant_at_the_clicked_residue_is_named(window):
    """R2456H is the project's one paired variant; it should be mentioned."""
    st = window.structure
    mask = st.mask_ca() & (st.res_seq == 2456)
    if not mask.any():
        pytest.skip("residue 2456 not resolved in this entry")
    menu = build_context_menu(window, int(np.flatnonzero(mask)[0]))
    assert any("Variant here" in t for t in labels(menu))


# --------------------------------------------- opening it must not change it

def test_building_the_menu_selects_nothing(window, atom):
    window._highlight([], "")
    before = list(window.selected_residues)
    build_context_menu(window, atom)
    assert list(window.selected_residues) == before
    assert window.view.highlight is None


def test_the_viewport_can_identify_an_atom_without_announcing_it(window):
    """`atom_at` is the half of picking the menu needs and selection does not."""
    from PyQt6.QtCore import QPoint

    seen = []
    connection = window.viewport.atom_picked.connect(seen.append)
    window.viewport.atom_at(QPoint(600, 400))
    window.viewport.atom_picked.disconnect(connection)
    assert not seen, "identifying an atom emitted a selection"


# --------------------------------------------------------- the entries act

def test_select_marks_only_the_chain_clicked(window, atom):
    st = window.structure
    menu = build_context_menu(window, atom)
    find(menu, f"(chain {st.chain[atom]})").trigger()

    assert window.view.highlight.sum() > 0
    assert set(st.chain[window.view.highlight]) == {str(st.chain[atom])}


def test_select_in_all_chains_marks_every_protomer(window, atom):
    st = window.structure
    menu = build_context_menu(window, atom)
    find(menu, "in all").trigger()
    assert len(set(st.chain[window.view.highlight])) == len(st.chains)


def test_selecting_a_whole_chain_marks_far_more_than_one_residue(window, atom):
    st = window.structure
    menu = build_context_menu(window, atom)
    find(menu, "whole of chain").trigger()
    marked = window.view.highlight
    assert set(st.chain[marked]) == {str(st.chain[atom])}
    assert marked.sum() > 1000, "a whole protomer is thousands of atoms"


def test_clear_selection_is_offered_only_when_there_is_one(window, atom):
    window._highlight([], "")
    assert not find(build_context_menu(window, atom), "Clear selection").isEnabled()

    build_context_menu(window, atom)
    find(build_context_menu(window, atom), "(chain").trigger()
    action = find(build_context_menu(window, atom), "Clear selection")
    assert action.isEnabled()
    action.trigger()
    assert window.view.highlight is None


def test_add_to_measurement_arms_picking_rather_than_refusing(window, atom):
    """This menu is where a user who never found the Measure button ends up."""
    panel = window.measure_panel
    panel.arm_button.setChecked(False)
    panel._clear()
    assert not panel.armed

    find(build_context_menu(window, atom), "Add to measurement").trigger()

    assert panel.armed, "the entry refused instead of arming"
    assert panel.set.pending == [atom]
    assert panel.table.rowCount() == 1
    panel.arm_button.setChecked(False)


def test_centre_the_view_moves_the_pivot_and_not_the_rotation(window, atom):
    camera = window.viewport.scene.camera
    pivot, rotation = camera.pivot.copy(), camera.rotation.copy()
    find(build_context_menu(window, atom), "Centre the view here").trigger()
    assert not np.allclose(pivot, camera.pivot)
    assert np.allclose(rotation, camera.rotation), \
        "centring threw away the orientation the user had set"


def test_copy_puts_the_residue_label_on_the_clipboard(window, atom):
    find(build_context_menu(window, atom), "Copy").trigger()
    st = window.structure
    expected = f"{st.res_name[atom]}{int(st.res_seq[atom])}{st.chain[atom]}"
    assert QApplication.clipboard().text() == expected


# ------------------------------------- appearance goes through the one panel

def test_representation_drives_the_model_panel_not_the_view_directly(window, atom):
    """The guard against a second implementation.

    Calling `_set_style` here would change the model and leave the Model
    panel's combo showing the old value — two sources of truth for one setting.
    """
    panel = window.structure_panel
    panel.style_combo.setCurrentIndex(0)
    menu = build_context_menu(window, atom)
    styles = find(menu, "Representation").menu()

    assert styles.actions()[0].isChecked(), "the current style is not ticked"
    styles.actions()[3].trigger()          # Spheres
    assert panel.style_combo.currentIndex() == 3, \
        "the menu changed the model without telling the panel"
    from piezo1.render.representations import Style
    assert window.view.style is Style.SPHERES
    panel.style_combo.setCurrentIndex(0)


def test_colour_by_drives_the_model_panel_too(window, atom):
    panel = window.structure_panel
    panel.color_combo.setCurrentIndex(0)
    colours = find(build_context_menu(window, atom), "Colour by").menu()
    assert colours.actions()[0].isChecked()
    colours.actions()[1].trigger()         # Chain / protomer
    assert panel.color_combo.currentIndex() == 1
    from piezo1.render.representations import ColorBy
    assert window.view.color_by is ColorBy.CHAIN
    panel.color_combo.setCurrentIndex(0)


def test_the_view_toggles_report_the_state_they_will_change(window, atom):
    camera = window.viewport.scene.camera
    camera.orthographic = False
    assert not find(build_context_menu(window, atom),
                    "Orthographic").isChecked()
    camera.orthographic = True
    assert find(build_context_menu(window, atom), "Orthographic").isChecked()
    camera.orthographic = False

    window.viewport.set_spin(0.0)
    spin = find(build_context_menu(window, atom), "Spin")
    assert not spin.isChecked()
    spin.trigger()
    assert window.viewport._spin_speed != 0.0
    window.viewport.set_spin(0.0)


# ------------------------------------------------------------- the binding

def _press(view, x, y, button):
    view.mousePressEvent(QMouseEvent(
        QMouseEvent.Type.MouseButtonPress, QPointF(x, y), button, button,
        Qt.KeyboardModifier.NoModifier))


def _move(view, x, y, button):
    view.mouseMoveEvent(QMouseEvent(
        QMouseEvent.Type.MouseMove, QPointF(x, y), Qt.MouseButton.NoButton,
        button, Qt.KeyboardModifier.NoModifier))


def _release(view, x, y, button):
    view.mouseReleaseEvent(QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease, QPointF(x, y), button,
        Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier))


def test_a_right_click_opens_the_menu(window, app):
    right = Qt.MouseButton.RightButton
    window._context_menu = None
    _press(window.viewport, 600, 400, right)
    _release(window.viewport, 600, 400, right)
    app.processEvents()
    assert window._context_menu is not None, "right-click opened nothing"
    assert labels(window._context_menu)
    window._context_menu.close()


def test_a_right_drag_zooms_and_opens_nothing(window, app):
    """Right-drag was the zoom before the menu existed and still is."""
    right = Qt.MouseButton.RightButton
    window._context_menu = None
    camera = window.viewport.scene.camera
    before = camera.distance

    _press(window.viewport, 600, 300, right)
    for y in range(310, 500, 20):
        _move(window.viewport, 600, y, right)
    _release(window.viewport, 600, 490, right)
    app.processEvents()

    assert abs(camera.distance - before) > 1e-9, "right-drag stopped zooming"
    assert window._context_menu is None, "a zoom drag popped up a menu"


def test_the_menu_is_shown_without_blocking(window):
    """`exec()` spins its own event loop and never returns until the user
    chooses; the whole path would be untestable and a hung test would be the
    first anyone knew."""
    import inspect

    from piezo1.ui.context_menu import ContextMenuMixin

    source = inspect.getsource(ContextMenuMixin._show_context_menu)
    assert ".popup(" in source
    assert ".exec(" not in source
