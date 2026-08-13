"""The Options menu: every entry observable, and every preference in it.

Two claims, both testable. **Consolidation**: Options holds what is
remembered across sessions — the four persisted preferences that lived under
View (alignment, multi-structure, companion style, display options) moved
here, and the two deliberate exceptions stay put and are pinned as such,
because the flux pathway and voltage change what is computed, not what is
preferred. **Effect**: an option that changes nothing observable is a lie in
a menu, so each one here is fired and its effect measured — the background
down to the rendered pixel, the theme down to the application palette.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
moderngl = pytest.importorskip("moderngl")

from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.config import STRUCTURE_DIR, RenderSettings  # noqa: E402


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

    # These tests change persisted settings and the application-wide theme;
    # both are put back exactly, so a test run does not restyle the user's
    # next session.
    saved = {key: win.settings.value(key) for key in win.settings.allKeys()}
    style_before = app.style().objectName()
    sheet_before = app.styleSheet()
    palette_before = app.palette()

    win.resize(1200, 800)
    win.show()
    app.processEvents()
    try:
        ctx = moderngl.create_standalone_context(require=410)
    except Exception as exc:                           # pragma: no cover
        pytest.skip(f"no OpenGL 4.1 context available: {exc}")
    # The window's own settings object, exactly as the real application
    # wires it — otherwise the background option would be asserted against
    # a copy the scene never reads.
    scene = Scene(ctx, win.viewport.settings)
    scene.resize(320, 240)
    win.viewport.scene = scene
    win._on_scene_ready(scene)
    app.processEvents()
    if win.structure is None:
        pytest.skip("no default structure could be loaded")
    win._fbo = ctx.simple_framebuffer((320, 240))
    win._scene_ctx = ctx
    yield win

    for controller in ("analysis", "physics", "overlay"):
        cleanup = getattr(getattr(win, controller, None), "cleanup", None)
        if cleanup is not None:
            cleanup()
    win.settings.clear()
    for key, value in saved.items():
        win.settings.setValue(key, value)
    # The viewport settings object is the process-wide SETTINGS.render;
    # leave it exactly as a fresh start would have it.
    win.viewport.settings.background = RenderSettings().background
    app.setStyle(style_before)
    app.setStyleSheet(sheet_before)
    app.setPalette(palette_before)
    app.processEvents()


def menu_titles(win) -> dict:
    """Top-level menu title -> set of every action/submenu label under it."""
    def walk(menu, labels):
        for action in menu.actions():
            if action.isSeparator():
                continue
            labels.add(action.text().replace("&", ""))
            if action.menu():
                walk(action.menu(), labels)

    out = {}
    for top in win.menuBar().actions():
        labels: set[str] = set()
        if top.menu():
            walk(top.menu(), labels)
        out[top.text().replace("&", "")] = labels
    return out


def find_action(win, menu_title: str, label: str):
    for top in win.menuBar().actions():
        if top.text().replace("&", "") != menu_title or not top.menu():
            continue

        def walk(menu):
            for action in menu.actions():
                if action.text().replace("&", "") == label:
                    return action
                if action.menu():
                    found = walk(action.menu())
                    if found is not None:
                        return found
            return None
        return walk(top.menu())
    return None


def rendered_corner(win) -> np.ndarray:
    """The bottom-left rendered pixel — pure background on every entry."""
    win._fbo.use()
    win.viewport.scene.resize(320, 240)
    win.viewport.scene.render()
    data = np.frombuffer(win._fbo.read(components=3), dtype=np.uint8)
    return data.reshape(240, 320, 3)[0, 0]


# ------------------------------------------------------------- consolidation

def test_every_persisted_preference_lives_under_options(window):
    menus = menu_titles(window)
    for label in ("Interface theme", "Viewport background", "Display options…",
                  "Structure alignment", "Show multiple structures at once",
                  "Extra structures style", "Spin speed",
                  "When something is selected", "Remember layout on exit"):
        assert label in menus["Options"], f"{label!r} is not under Options"
        assert label not in menus["View"], f"{label!r} is still under View"


def test_the_two_deliberate_exceptions_stay_by_their_features(window):
    """The flux pathway and voltage change what is COMPUTED — each moves the
    physics of the number on the status line — so they stay beside the
    animation they parameterise, exactly as `menus_flux.py` argues. If this
    fails because they moved, move the reasoning too, not just the menu."""
    menus = menu_titles(window)
    assert "Ion flux pathway" in menus["View"]
    assert "Ion flux voltage" in menus["View"]
    assert "Ion flux pathway" not in menus["Options"]


# ------------------------------------------------------- the two new options

def test_the_background_option_changes_the_rendered_pixels(window):
    dark = rendered_corner(window)

    find_action(window, "Options", "White").setChecked(True)
    assert window.viewport.settings.background == (1.0, 1.0, 1.0, 1.0)
    assert window.viewport.scene.settings.background == (1.0, 1.0, 1.0, 1.0), \
        "the scene is reading a different settings object than the option set"
    white = rendered_corner(window)
    assert white.min() > 240, f"white background rendered as {white}"
    assert dark.max() < 40, f"the default background rendered as {dark}"

    find_action(window, "Options", "Midnight (default)").setChecked(True)
    assert (window.viewport.settings.background
            == RenderSettings().background)
    assert rendered_corner(window).max() < 40


def test_the_background_is_also_the_fog_colour(window):
    """The depth cue fades into the background; a background change that
    left the fog behind would haze every atom toward the old colour."""
    find_action(window, "Options", "White").setChecked(True)
    assert window.viewport.scene.uniforms()["u_fog_color"] == (1.0, 1.0, 1.0)
    find_action(window, "Options", "Midnight (default)").setChecked(True)
    expected = RenderSettings().background[:3]
    assert window.viewport.scene.uniforms()["u_fog_color"] == expected


def test_the_theme_option_restyles_the_application(window, app):
    from PyQt6.QtGui import QPalette

    def window_lightness() -> float:
        return app.palette().color(QPalette.ColorRole.Window).lightnessF()

    find_action(window, "Options", "Light").setChecked(True)
    assert window.ui_theme() == "light"
    assert window_lightness() > 0.5, "the light theme left a dark palette"
    assert app.styleSheet(), "the light theme cleared the stylesheet"

    find_action(window, "Options", "System").setChecked(True)
    assert app.styleSheet() == "", "System should hand styling to the platform"

    find_action(window, "Options", "Dark (default)").setChecked(True)
    assert window.ui_theme() == "dark"
    assert window_lightness() < 0.35, "the dark theme left a light palette"


# ------------------------------------------------- the moved options still act

def test_the_alignment_option_still_changes_the_frame(window, app):
    find_action(window, "Options", "As deposited").setChecked(True)
    app.processEvents()
    assert window.alignment_mode == "deposited"
    find_action(window, "Options",
                "Canonical (three-fold axis on z)").setChecked(True)
    app.processEvents()
    assert window.alignment_mode == "canonical"


def test_the_multi_structure_cluster_still_acts(window):
    from piezo1.render.representations import Style

    action = find_action(window, "Options", "Show multiple structures at once")
    was = window.multi_structure
    action.setChecked(not was)
    assert window.multi_structure is (not was)
    action.setChecked(was)
    assert window.multi_structure is was

    find_action(window, "Options", "Cartoon").setChecked(True)
    assert window.companion_style() is Style.CARTOON
    find_action(window, "Options", "Backbone trace").setChecked(True)
    assert window.companion_style() is Style.BACKBONE


def test_every_simple_option_has_an_observable_effect(window):
    find_action(window, "Options", "Show status-bar hints").setChecked(False)
    assert not window.hint_label.isVisible()
    find_action(window, "Options", "Show status-bar hints").setChecked(True)
    assert window.hint_label.isVisible()

    find_action(window, "Options", "Fast").setChecked(True)
    assert window._spin_speed() == 60.0
    find_action(window, "Options", "Normal").setChecked(True)
    assert window._spin_speed() == 28.0

    find_action(window, "Options", "Centre on the selection").setChecked(True)
    assert window.focus_mode() == "centre"
    find_action(window, "Options", "Keep the view still").setChecked(True)
    assert window.focus_mode() == "none"

    action = find_action(window, "Options", "Remember layout on exit")
    action.setChecked(False)
    assert window.settings.value("options/remember_layout", type=bool) is False
    action.setChecked(True)
    assert window.settings.value("options/remember_layout", type=bool) is True


def test_reset_restores_the_appearance_immediately(window, app):
    """A reset that left a white viewport or a light chrome behind would not
    be a reset — the two appearance options act on choice, so they must act
    on forgetting too."""
    from PyQt6.QtGui import QPalette

    find_action(window, "Options", "White").setChecked(True)
    find_action(window, "Options", "Light").setChecked(True)
    window._reset_options()
    assert window.background_key() == "midnight"
    assert window.viewport.settings.background == RenderSettings().background
    assert window.ui_theme() == "dark"
    assert app.palette().color(
        QPalette.ColorRole.Window).lightnessF() < 0.35
