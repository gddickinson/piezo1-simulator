"""The right-click menu: what to do with the thing under the cursor.

Two rules shape what is in here.

**Nothing is implemented twice.** Every entry drives the same widget or
controller the rest of the interface uses — the Model panel's own combo boxes
for style and colouring, `_highlight` for selection, the Measure panel for
picking. Setting the combo rather than calling `_set_style` directly matters:
otherwise the panel would go on saying "Cartoon" while the model showed
spheres, and the two would drift apart the first time anyone used this menu.

**Opening the menu is not a selection.** A right-click identifies what is under
the cursor so the entries can name it, but changes nothing; dismissing the menu
leaves the model exactly as it was. That is why the viewport has `atom_at`
alongside `_pick_at` — the first answers the question, the second announces it.

The menu is context-sensitive in one direction only: with an atom under the
cursor it gains the entries that act on that residue, and without one it keeps
everything that acts on the view. It is never empty, because a menu that
sometimes does not appear reads as a broken control.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtGui import QAction, QActionGroup
from PyQt6.QtWidgets import QMenu

from .panels.structure_panel import COLOR_LABELS, STYLE_LABELS

__all__ = ["ContextMenuMixin", "build_context_menu"]


def _header(menu: QMenu, text: str) -> None:
    """A non-clickable caption naming what the menu was opened on."""
    action = QAction(text, menu)
    action.setEnabled(False)
    menu.addAction(action)
    menu.addSeparator()


def build_context_menu(win, index: int) -> QMenu:
    """The menu for a right-click on atom ``index`` (-1 for empty space)."""
    menu = QMenu(win)
    menu.setToolTipsVisible(True)

    structure = win.structure
    if structure is not None and 0 <= index < structure.n_atoms:
        _residue_entries(win, menu, index)

    _selection_entries(win, menu)
    _appearance_entries(win, menu)
    _view_entries(win, menu)
    return menu


# --------------------------------------------------------- the clicked atom

def _residue_entries(win, menu: QMenu, index: int) -> None:
    st = win.structure
    residue = int(st.res_seq[index])
    chain = str(st.chain[index])
    name = str(st.res_name[index])
    label = f"{name}{residue}{chain}"

    info = win.annotations.annotate_residue(residue)
    caption = f"{name}{residue} · chain {chain} · {st.atom_name[index]}"
    if info["domain"]:
        caption += f" · {info['domain']}"
    _header(menu, caption)

    act = menu.addAction(f"Select {name}{residue} (chain {chain})")
    act.setToolTip("Mark this residue on the model — the copy you clicked")
    act.triggered.connect(
        lambda: win._highlight([residue], label, chains=[chain]))

    n_chains = len(st.chains)
    act = menu.addAction(f"Select {name}{residue} in all {n_chains} chains")
    act.setToolTip("The same residue number in every protomer, which is what "
                   "an annotation means by a residue")
    act.triggered.connect(lambda: win._highlight([residue], f"{name}{residue}"))

    act = menu.addAction(f"Select the whole of chain {chain}")
    act.setToolTip("Every modelled residue of this protomer")
    act.triggered.connect(lambda: _select_chain(win, chain))

    menu.addSeparator()

    act = menu.addAction("Add to measurement")
    act.setToolTip("Start picking if it is not already on, and use this atom "
                   "as the next point of the measurement")
    act.triggered.connect(lambda: _add_to_measurement(win, index, label))

    act = menu.addAction("Centre the view here")
    act.setToolTip("Move the camera pivot to this residue without reframing")
    act.triggered.connect(lambda: _centre_on(win, residue, chain))

    act = menu.addAction(f"Copy “{label}”")
    act.setToolTip("Put the residue label on the clipboard")
    act.triggered.connect(lambda: _copy(win, label))

    for variant in info["variants"]:
        act = menu.addAction(
            f"Variant here: {variant['label']} ({variant['classification']})")
        act.setEnabled(False)

    menu.addSeparator()


def _select_chain(win, chain: str) -> None:
    st = win.structure
    residues = sorted({int(r) for r in st.res_seq[st.chain == chain]})
    win._highlight(residues, f"chain {chain}", chains=[chain])


def _add_to_measurement(win, index: int, label: str) -> None:
    """Route this atom into the Measure panel, arming it if need be.

    Arming here rather than refusing is the point: this menu is where a user
    who has not found the Measure button ends up, and telling them to go and
    press it first would repeat the problem the menu is meant to solve.
    """
    panel = win.measure_panel
    if not panel.armed:
        panel.arm_button.setChecked(True)
    st = win.structure
    panel.add_pick(int(index), st.xyz[index].astype(float), label)


def _centre_on(win, residue: int, chain: str) -> None:
    st = win.structure
    mask = (st.res_seq == residue) & (st.chain == chain)
    if not mask.any() or win.viewport.scene is None:
        return
    win.viewport.scene.camera.pivot = st.xyz[mask].mean(axis=0)
    win.viewport.update()


def _copy(win, text: str) -> None:
    from PyQt6.QtWidgets import QApplication

    QApplication.clipboard().setText(text)
    win._set_status(f"copied “{text}”")


# ------------------------------------------------------------- the selection

def _selection_entries(win, menu: QMenu) -> None:
    act = menu.addAction("Clear selection")
    act.setEnabled(bool(getattr(win, "selected_residues", None)))
    act.setToolTip("Unmark whatever is currently highlighted")
    act.triggered.connect(lambda: win._highlight([], ""))
    menu.addSeparator()


# ------------------------------------------------------------- how it looks

def _appearance_entries(win, menu: QMenu) -> None:
    """Style and colouring, driven through the Model panel's own combos.

    Calling `_set_style` directly would leave the panel showing the old value.
    """
    panel = win.structure_panel

    styles = menu.addMenu("Representation")
    styles.setToolTipsVisible(True)
    group = QActionGroup(styles)
    for i, (text, _style) in enumerate(STYLE_LABELS):
        act = styles.addAction(text)
        act.setCheckable(True)
        act.setChecked(panel.style_combo.currentIndex() == i)
        group.addAction(act)
        act.triggered.connect(
            lambda _checked=False, i=i: panel.style_combo.setCurrentIndex(i))

    colours = menu.addMenu("Colour by")
    colours.setToolTipsVisible(True)
    group = QActionGroup(colours)
    for i, (text, _color) in enumerate(COLOR_LABELS):
        act = colours.addAction(text)
        act.setCheckable(True)
        act.setChecked(panel.color_combo.currentIndex() == i)
        group.addAction(act)
        act.triggered.connect(
            lambda _checked=False, i=i: panel.color_combo.setCurrentIndex(i))
    menu.addSeparator()


# ------------------------------------------------------------- the camera

def _view_entries(win, menu: QMenu) -> None:
    scene = win.viewport.scene

    act = menu.addAction("Reset the view\tR")
    act.setToolTip("Frame the whole model again")
    act.triggered.connect(win._reset_camera)

    act = menu.addAction("Spin\tSpace")
    act.setCheckable(True)
    act.setChecked(bool(win.viewport._spin_speed))
    act.toggled.connect(
        lambda on: win.viewport.set_spin(30.0 if on else 0.0))

    act = menu.addAction("Orthographic projection\tO")
    act.setCheckable(True)
    act.setChecked(bool(scene is not None and scene.camera.orthographic))
    act.setToolTip("No perspective foreshortening — parallel lines stay "
                   "parallel, which is what a measured figure wants")
    act.toggled.connect(lambda on: _set_orthographic(win, on))


def _set_orthographic(win, on: bool) -> None:
    if win.viewport.scene is None:
        return
    win.viewport.scene.camera.orthographic = bool(on)
    win.viewport.update()


# --------------------------------------------------------------- the mixin

class ContextMenuMixin:
    """Shows the right-click menu. Mixed into MainWindow."""

    #: The menu currently on screen. Held as an attribute for two reasons: a
    #: menu shown with ``popup()`` is garbage-collected the moment the last
    #: reference goes, and it is the only way a test can see what was built.
    _context_menu = None

    def _show_context_menu(self, pos, index: int) -> None:
        """Pop the menu up at ``pos``, non-modally.

        ``popup()`` rather than ``exec()`` because ``exec()`` spins its own
        event loop and does not return until the user chooses something — which
        is fine for a person and impossible for a test, so the whole path would
        have gone uncovered.
        """
        if self.structure is None:
            return
        self._context_menu = build_context_menu(self, int(index))
        self._context_menu.popup(self.viewport.mapToGlobal(pos))
