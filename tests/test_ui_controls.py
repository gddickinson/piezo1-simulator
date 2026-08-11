"""Every control the application offers, exercised rather than assumed.

Asked to confirm that the mouse, keyboard and menu controls all work, the only
honest way to answer was to fire them. Reading the handlers is not enough — the
bug this file was written around looked completely correct in source:
``mouseReleaseEvent`` compared the release position against ``_last_pos`` to
tell a click from a drag, and ``_last_pos`` is updated by every move event, so
the distance was always zero and **every rotation ended in a pick**.

The window runs against a real standalone GL context rather than a mocked
scene, because a mock would have accepted the same wrong answer.
"""

from __future__ import annotations

import os
import re
from pathlib import Path

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
moderngl = pytest.importorskip("moderngl")

from PyQt6.QtCore import Qt, QPoint, QPointF  # noqa: E402
from PyQt6.QtGui import QKeyEvent, QMouseEvent, QWheelEvent  # noqa: E402
from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.config import RenderSettings, STRUCTURE_DIR  # noqa: E402
from piezo1.ui.help_content import SHORTCUTS  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


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
    """A real MainWindow with a real GL context and a structure loaded."""
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
    yield win

    # Firing every menu action starts background analyses, and nothing here
    # waits for them. A worker still running a pore profile or an eigensolve
    # while the next module does its own numpy work segfaults the interpreter —
    # two threads inside the same LAPACK. The window already knows how to stop
    # them; the fixture just never asked.
    for controller in ("analysis", "physics", "overlay"):
        owner = getattr(win, controller, None)
        cleanup = getattr(owner, "cleanup", None)
        if cleanup is not None:
            cleanup()
    app.processEvents()


def menu_actions(win) -> list:
    """Every leaf action in the menu bar, with its full path."""
    found = []

    def walk(menu, path):
        for action in menu.actions():
            if action.isSeparator():
                continue
            name = f"{path}{action.text().replace('&', '')}"
            if action.menu():
                walk(action.menu(), name + " > ")
            else:
                found.append((name, action))

    for action in win.menuBar().actions():
        if action.menu():
            walk(action.menu(), action.text().replace("&", "") + " > ")
    return found


# ----------------------------------------------------------------- the menus

def test_every_menu_action_is_reachable_and_explained(window):
    actions = menu_actions(window)
    assert len(actions) > 50, f"only {len(actions)} menu actions found"
    for name, action in actions:
        assert action.text().strip(), "an action with no label"
        assert action.toolTip() and action.toolTip() != action.text(), \
            f"{name} has no tooltip of its own"
        assert action.isEnabled(), f"{name} is disabled with a structure loaded"


def test_no_two_actions_claim_the_same_shortcut(window):
    claimed = {}
    for name, action in menu_actions(window):
        key = action.shortcut().toString()
        if key:
            claimed.setdefault(key, []).append(name)
    clashes = {k: v for k, v in claimed.items() if len(v) > 1}
    assert not clashes, f"shortcut collisions: {clashes}"


def test_every_action_runs_without_raising(window, app):
    """Fired for real, with a structure loaded, and put back if it toggles.

    Dialog actions are excluded — a modal would block the run — but their
    handlers are covered by the panel suites. Everything else executes.
    """
    ran, failures = 0, []
    for name, action in menu_actions(window):
        if name.endswith("…") or "Quit" in name:
            continue
        try:
            action.trigger()
            app.processEvents()
            if action.isCheckable():
                action.trigger()
                app.processEvents()
            ran += 1
        except Exception as exc:
            failures.append(f"{name} -> {type(exc).__name__}: {exc}")
    assert not failures, "\n".join(failures)
    assert ran > 40, f"only {ran} actions were actually fired"


# ----------------------------------------------------------------- the mouse

def _press(view, x, y, button=Qt.MouseButton.LeftButton,
           mods=Qt.KeyboardModifier.NoModifier):
    view.mousePressEvent(QMouseEvent(
        QMouseEvent.Type.MouseButtonPress, QPointF(x, y), button, button, mods))


def _move(view, x, y, button, mods=Qt.KeyboardModifier.NoModifier):
    view.mouseMoveEvent(QMouseEvent(
        QMouseEvent.Type.MouseMove, QPointF(x, y), Qt.MouseButton.NoButton,
        button, mods))


def _release(view, x, y, button=Qt.MouseButton.LeftButton):
    view.mouseReleaseEvent(QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease, QPointF(x, y), button,
        Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier))


def _drag(view, start, end, button=Qt.MouseButton.LeftButton,
          mods=Qt.KeyboardModifier.NoModifier, steps=20):
    """A drag as the window manager delivers one: many small move events."""
    _press(view, *start, button=button, mods=mods)
    for i in range(1, steps + 1):
        _move(view, start[0] + (end[0] - start[0]) * i / steps,
              start[1] + (end[1] - start[1]) * i / steps, button, mods)
    _release(view, *end, button=button)


def test_left_drag_rotates(window):
    camera = window.viewport.scene.camera
    before = camera.rotation.copy()
    _drag(window.viewport, (300, 300), (700, 400))
    assert not np.allclose(before, camera.rotation)


def test_shift_drag_and_middle_drag_both_pan(window):
    camera = window.viewport.scene.camera
    before = camera.pan.copy()
    _drag(window.viewport, (300, 300), (400, 350),
          mods=Qt.KeyboardModifier.ShiftModifier)
    assert not np.allclose(before, camera.pan), "shift+drag did not pan"

    before = camera.pan.copy()
    _drag(window.viewport, (300, 300), (400, 350),
          button=Qt.MouseButton.MiddleButton)
    assert not np.allclose(before, camera.pan), "middle-drag did not pan"


def test_right_drag_and_wheel_both_zoom(window):
    camera = window.viewport.scene.camera
    before = camera.distance
    _drag(window.viewport, (600, 300), (600, 500),
          button=Qt.MouseButton.RightButton)
    assert abs(camera.distance - before) > 1e-9, "right-drag did not zoom"

    before = camera.distance
    window.viewport.wheelEvent(QWheelEvent(
        QPointF(600, 400), QPointF(600, 400), QPoint(0, 0), QPoint(0, 120),
        Qt.MouseButton.NoButton, Qt.KeyboardModifier.NoModifier,
        Qt.ScrollPhase.NoScrollPhase, False))
    assert abs(camera.distance - before) > 1e-9, "the wheel did not zoom"


# ----------------------------------------- the click that a drag must not be

@pytest.fixture
def picks(window):
    seen = []
    connection = window.viewport.atom_picked.connect(seen.append)
    yield seen
    window.viewport.atom_picked.disconnect(connection)


def test_a_click_picks(window, picks):
    _press(window.viewport, 600, 400)
    _release(window.viewport, 600, 400)
    assert picks, "a click did not pick at all"


def test_a_little_tremor_is_still_a_click(window, picks):
    """Trackpads move a pixel or two under a deliberate tap."""
    from piezo1.ui.gl_widget import CLICK_SLOP

    _press(window.viewport, 600, 400)
    _move(window.viewport, 601, 400, Qt.MouseButton.LeftButton)
    _release(window.viewport, 601, 400)
    assert picks, f"a 1 px wobble was rejected; CLICK_SLOP is {CLICK_SLOP}"


def test_a_rotation_drag_never_picks(window, picks):
    """The regression this file exists for.

    ``mouseReleaseEvent`` measured the distance from ``_last_pos``, which every
    move event overwrites, so at release it *was* the release position and the
    distance was always zero. A 390-pixel rotation picked an atom. Harmless
    while a pick only rewrote the status bar; once a pick highlights, every
    turn of the structure repainted the selection.
    """
    _drag(window.viewport, (300, 300), (690, 300))
    assert not picks, "rotating the structure picked an atom"


def test_a_drag_that_ends_where_it_ended_still_does_not_pick(window, picks):
    """The exact shape of the old bug: no move event between the last step and
    the release, so press-vs-release is the only comparison that can work."""
    _press(window.viewport, 300, 300)
    _move(window.viewport, 500, 300, Qt.MouseButton.LeftButton)
    _release(window.viewport, 500, 300)
    assert not picks


def test_only_the_left_button_picks(window, picks):
    """Right-drag zooms; releasing it must not also select something."""
    _press(window.viewport, 600, 400, button=Qt.MouseButton.RightButton)
    _release(window.viewport, 600, 400, button=Qt.MouseButton.RightButton)
    assert not picks


# -------------------------------------------------------------- the keyboard

def test_the_viewport_can_receive_keys_at_all(window):
    """Every viewport key is dead if the widget cannot take focus."""
    assert window.viewport.focusPolicy() in (
        Qt.FocusPolicy.StrongFocus, Qt.FocusPolicy.WheelFocus)


def _key(window, key):
    window.viewport.keyPressEvent(QKeyEvent(
        QKeyEvent.Type.KeyPress, key, Qt.KeyboardModifier.NoModifier))


def test_r_resets_the_camera(window):
    camera = window.viewport.scene.camera
    camera.distance = 999.0
    _key(window, Qt.Key.Key_R)
    assert abs(camera.distance - 999.0) > 1e-9


def test_o_toggles_the_projection(window):
    camera = window.viewport.scene.camera
    before = camera.orthographic
    _key(window, Qt.Key.Key_O)
    assert camera.orthographic != before
    _key(window, Qt.Key.Key_O)
    assert camera.orthographic == before


def test_space_toggles_spin_both_ways(window):
    _key(window, Qt.Key.Key_Space)
    assert window.viewport._spin_speed != 0.0
    _key(window, Qt.Key.Key_Space)
    assert window.viewport._spin_speed == 0.0


def test_plus_and_minus_resize_the_atoms(window):
    scene = window.viewport.scene
    before = scene.radius_scale
    _key(window, Qt.Key.Key_Plus)
    assert scene.radius_scale > before
    _key(window, Qt.Key.Key_Minus)
    assert scene.radius_scale == pytest.approx(before)


# ------------------------------------------- the guide against the bindings

def test_every_menu_shortcut_is_documented(window):
    """A binding nobody can discover is not a feature."""
    documented = " ".join(k for k, _ in SHORTCUTS)
    missing = []
    for name, action in menu_actions(window):
        key = action.shortcut().toString()
        if key and key not in documented:
            missing.append(f"{key} ({name})")
    assert not missing, f"shortcuts absent from the help: {missing}"


def test_every_viewport_key_is_documented():
    """Read out of the handler, so a new key cannot be added unannounced.

    Four bindings — middle-drag, right-drag, O and +/- — were undocumented
    until this check was written.
    """
    source = (ROOT / "piezo1" / "ui" / "gl_widget.py").read_text()
    body = source[source.index("def keyPressEvent"):]
    body = body[:body.index("\n    def ", 1)]
    handled = set(re.findall(r"Qt\.Key\.Key_(\w+)", body))

    documented = " ".join(k for k, _ in SHORTCUTS).lower()
    aliases = {"plus": "+", "equal": "+", "minus": "-", "space": "space"}
    missing = [k for k in handled
               if aliases.get(k.lower(), k.lower()) not in documented]
    assert not missing, f"viewport keys absent from the help: {sorted(missing)}"


def test_the_mouse_bindings_are_documented():
    documented = " ".join(k for k, _ in SHORTCUTS).lower()
    for binding in ("drag", "shift + drag", "middle-drag", "right-drag",
                    "right-click", "wheel", "click"):
        assert binding in documented, f"{binding} is not in the shortcut list"
