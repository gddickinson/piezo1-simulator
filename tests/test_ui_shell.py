"""The application shell: docks, menus, preferences and help.

Split from `test_ui_analysis.py` when that file passed the project's 500-line
limit. Runs real widgets on the **offscreen** Qt platform rather than mocks —
both of this project's silent GUI breakages were real widgets behaving
differently from how the code assumed, which is precisely what a mock encodes
rather than tests.
"""

import os

import numpy as np
import re

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
from PyQt6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance()
    if app is None:
        try:
            app = QApplication([])
        except Exception as exc:                       # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return app


# --------------------------------------------------------------------------
# Docking, menus and preferences
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def window(qapp):
    """A real MainWindow. Skips where no GL context can be made."""
    from piezo1.config import SETTINGS
    from piezo1.ui.gl_widget import configure_surface_format
    configure_surface_format(SETTINGS.render)
    try:
        from piezo1.ui.main_window import MainWindow
        win = MainWindow()
    except Exception as exc:                           # pragma: no cover
        pytest.skip(f"cannot construct MainWindow: {exc}")
    win.show()
    qapp.processEvents()
    yield win
    win.settings.clear()
    win.close()


def test_every_panel_is_a_full_dock(window):
    """Movable, floatable and closable, and allowed in all four areas.

    Qt's defaults leave a dock confined to the areas it was created for, so
    these have to be set explicitly; forgetting is invisible until a user tries
    to drag a panel somewhere and it refuses.
    """
    from PyQt6.QtCore import Qt
    from PyQt6.QtWidgets import QDockWidget

    assert set(window.docks.docks) == {"model", "physics", "annotation",
                                       "analysis", "measure", "overlay",
                                       "tour"}
    for key, dock in window.docks.docks.items():
        features = dock.features()
        for wanted in (QDockWidget.DockWidgetFeature.DockWidgetMovable,
                       QDockWidget.DockWidgetFeature.DockWidgetFloatable,
                       QDockWidget.DockWidgetFeature.DockWidgetClosable):
            assert features & wanted, f"{key} is missing {wanted}"
        for area in (Qt.DockWidgetArea.LeftDockWidgetArea,
                     Qt.DockWidgetArea.RightDockWidgetArea,
                     Qt.DockWidgetArea.TopDockWidgetArea,
                     Qt.DockWidgetArea.BottomDockWidgetArea):
            assert dock.isAreaAllowed(area), f"{key} rejects {area}"
        assert dock.objectName(), f"{key} has no objectName; saveState needs it"


def test_reset_layout_undoes_hiding_and_floating(window, qapp):
    docks = window.docks.docks
    docks["measure"].hide()
    docks["model"].setFloating(True)
    qapp.processEvents()
    assert not docks["measure"].isVisible()
    assert docks["model"].isFloating()

    window._reset_layout()
    qapp.processEvents()
    assert docks["measure"].isVisible(), "reset did not restore a closed panel"
    assert not docks["model"].isFloating(), "reset did not re-dock a panel"


def test_float_and_dock_all(window, qapp):
    window.docks.float_all(True)
    qapp.processEvents()
    assert all(d.isFloating() for d in window.docks.docks.values())
    window.docks.float_all(False)
    qapp.processEvents()
    assert not any(d.isFloating() for d in window.docks.docks.values())
    window._reset_layout()


def test_menus_exist_with_their_actions(window):
    menus = {a.text().replace("&", ""): a.menu()
             for a in window.menuBar().actions() if a.menu()}
    assert set(menus) == {"File", "View", "Analysis", "Options", "Help"}
    labels = {name: [x.text().replace("&", "") for x in menu.actions()]
              for name, menu in menus.items()}
    assert "Reset layout" in labels["View"]
    assert "Export analysis report…" in labels["File"]
    assert any("not do" in x for x in labels["Help"])
    assert "Panels" in labels["View"]


def test_view_menu_has_one_toggle_per_panel(window):
    actions = dict(window.docks.view_actions())
    assert set(actions) == set(window.docks.docks)
    for action in actions.values():
        assert action.isCheckable()


def test_window_fits_on_its_screen(window):
    """The complaint that started this: a hard 1680x1000 puts the title bar
    off the top of a laptop display, and some window managers then make the
    window unmovable."""
    screen = window.screen()
    if screen is None:
        pytest.skip("no screen")
    available = screen.availableGeometry()
    assert window.width() <= available.width()
    assert window.height() <= available.height()


def test_window_can_shrink_below_its_panels(window):
    """Docks wrap their panels in scroll areas so the window is not held open
    by the tallest one."""
    hint = window.minimumSizeHint()
    assert hint.height() < 700, hint.height()
    assert hint.width() < 900, hint.width()


def test_focus_mode_defaults_to_leaving_the_view_alone(window):
    assert window.focus_mode() == "none"


def test_focus_mode_controls_whether_the_camera_moves(window, qapp):
    """The reported problem: clicking a list entry moved the whole model."""
    import numpy as np
    window.load_structure("8YEZ")
    qapp.processEvents()
    if window.structure is None or window.viewport.scene is None:
        pytest.skip("structure or GL scene unavailable")

    camera = window.viewport.scene.camera
    residues = [int(window.structure.res_seq[100])]

    window._set_focus_mode("none")
    before = np.array(camera.pivot, dtype=float).copy()
    window._focus_residues(residues)
    assert np.allclose(camera.pivot, before), "view moved with focus off"

    window._set_focus_mode("centre")
    window._focus_residues(residues)
    assert not np.allclose(camera.pivot, before), "view did not move with focus on"
    window._set_focus_mode("none")


def test_help_dialog_covers_every_panel():
    from piezo1.ui.help_content import DOC_LINKS, SHORTCUTS, TOPICS
    titles = [t for t, _ in TOPICS]
    for panel in ("Model", "Annotation", "Physics", "Analysis", "Measure"):
        assert any(panel in t for t in titles), f"no help topic for {panel}"
    assert any("Limits" in t for t in titles)
    assert len(SHORTCUTS) > 8
    assert all(len(entry) == 3 for entry in DOC_LINKS)


def test_help_records_the_null_result():
    """The guide must not quietly omit that the central claim failed."""
    from piezo1.ui.help_content import TOPICS
    body = dict(TOPICS)["Limits and honesty"]
    assert "0.234" in body and "0.542" in body
    assert "99.8%" in body
    assert "179 nm" in body, "the corrected footprint area should be stated"


def test_shipped_documents_referenced_by_help_exist():
    from piezo1.config import PROJECT_ROOT
    from piezo1.ui.help_content import DOC_LINKS
    missing = [path for _title, path, _d in DOC_LINKS
               if not (PROJECT_ROOT / path).exists()]
    assert not missing, f"help links to documents that do not exist: {missing}"


def test_controls_carry_tooltips(qapp):
    """A tooltip is where the provenance of a number lives in a GUI."""
    from piezo1.ui.panels.analysis_panel import AnalysisPanel
    from piezo1.ui.panels.physics_panel import PhysicsPanel
    from piezo1.ui.panels.structure_panel import StructurePanel

    analysis = AnalysisPanel()
    assert "0.59" in analysis.hydro_check.toolTip()
    assert analysis.pore_button.toolTip()
    assert analysis.pockets_button.toolTip()

    physics = PhysicsPanel()
    assert "9.7" in physics.measure_button.toolTip()
    assert physics.compute_button.toolTip()

    structure = StructurePanel()
    assert structure.structure_combo.toolTip()
    assert structure.ligand_check.toolTip()


# --------------------------------------------------------------------------
# Everything computable must be reachable from the GUI, and explained there
# --------------------------------------------------------------------------

def test_every_shared_analysis_is_reachable_from_the_gui(qapp):
    """A function only the command line can run is invisible to most users.

    `permeation` and `interactions` were exactly that until Round 34: present in
    the shared ANALYSES registry, wired into the CLI, and absent from every
    menu. This checks the window exposes a way to run each one.
    """
    from piezo1.analysis.report import ANALYSES
    from piezo1.ui.main_window import MainWindow

    # Analyses drawn on the model rather than tabulated, and where they live.
    drawn = {"dome": "physics.measure_dome", "pore": "analysis.compute_pore",
             "hydration": "analysis.compute_pore", "modes": "physics.compute_modes",
             "pockets": "analysis.compute_pockets"}
    tabular = {"permeation": "show_permeation", "interactions": "show_interactions",
               "labelling": "show_labelling", "fusion": "show_fusion_numbers",
               "nanodomain": "show_nanodomain",
               "prediction_record": "show_prediction_record",
               "ligands": "show_ligands",
               "paired_variant": "show_paired_variant",
               "hybrid": "show_hybrid",
               "fluctuations": "show_fluctuations"}

    for name in ANALYSES:
        assert name in drawn or name in tabular, (
            f"{name!r} is in the ANALYSES registry but has no GUI entry point; "
            f"add one or record why it is command-line only")

    for attribute in tabular.values():
        assert callable(getattr(MainWindow, attribute, None)), attribute


def test_every_menu_action_has_a_tooltip(qapp):
    """A menu entry whose meaning is not stated is a menu entry nobody uses."""
    from piezo1.ui.main_window import MainWindow

    window = MainWindow()
    try:
        missing = []

        def walk(menu):
            for action in menu.actions():
                if action.isSeparator():
                    continue
                if action.menu() is not None:
                    walk(action.menu())
                elif (not action.toolTip()
                      or action.toolTip() == action.text().replace("&", "")):
                    missing.append(action.text())

        for action in window.menuBar().actions():
            if action.menu() is not None:
                walk(action.menu())
        assert not missing, f"menu actions without a tooltip: {missing}"
    finally:
        window.close()


def test_the_guide_covers_what_the_recent_rounds_added(qapp):
    """Help that stops at Round 30 is help that misdescribes the application."""
    from piezo1.ui.help_content import TOPICS

    # Whitespace-normalised: the guide is hand-wrapped HTML, so a phrase can
    # sit across a line break and a naive substring test then silently passes
    # or silently fails for a reason that has nothing to do with the content.
    text = re.sub(r"\s+", " ", " ".join(body for _t, body in TOPICS)).lower()
    for phrase in ("halotag", "accessible volume", "permeation", "pore mouths",
                   "canonical", "multiple structures", "unreactive tags"):
        assert phrase in text, f"the guide never mentions {phrase!r}"

    # And it must keep saying what the model cannot do.
    assert "16-94 ps" in text.replace("\u2013", "-")
    assert "no deposited loss-of-function structure" in text
