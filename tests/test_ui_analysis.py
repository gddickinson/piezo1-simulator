"""The Analysis dock, the profile plot, and session round-tripping.

Qt refactors have silently broken the GUI twice in this project, so these run
the real widgets on the **offscreen** platform rather than mocking them. They
skip if a QApplication cannot be created at all.

What is deliberately *not* tested here is scientific content — the pore profile
and the wetting verdict have their own suites. What is tested is that the panel
receives, transforms and hands back what the engine produced, since that is the
layer where a refactor breaks things without any exception being raised.
"""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.io.session import Session, load_session, save_session  # noqa: E402


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
# The profile plot
# --------------------------------------------------------------------------

def test_two_axes_scale_independently(qapp):
    """Radius and hydrophobicity share an x-axis and nothing else.

    Putting them on one y-scale is the obvious shortcut and makes the
    hydrophobicity trace a flat line at the bottom of a 0-10 Angstrom axis,
    which is exactly the comparison the plot exists to show.
    """
    from piezo1.ui.profile_plot import ProfilePlot, Trace
    plot = ProfilePlot()
    z = np.linspace(-60.0, 20.0, 80)
    plot.set_data([Trace("radius", z, np.abs(z) / 8.0, axis=0),
                   Trace("hydrophobicity", z, np.sin(z / 10.0) * 0.4, axis=1)])
    assert plot.left.hi > 5.0
    assert plot.right.hi < 1.0
    assert plot.right.lo < 0.0


def test_left_axis_anchors_at_zero_but_right_does_not(qapp):
    """A radius of zero is meaningful; a hydrophobicity of zero is not a floor."""
    from piezo1.ui.profile_plot import ProfilePlot, Trace
    plot = ProfilePlot()
    x = np.arange(20.0)
    plot.set_data([Trace("radius", x, x / 4.0 + 3.0, axis=0),
                   Trace("hydro", x, np.full(20, -0.4), axis=1)])
    assert plot.left.lo <= 0.0
    assert plot.right.lo > -1.0


def test_non_finite_values_do_not_break_scaling(qapp):
    """Hydrophobicity is NaN where no residue lines the pore."""
    from piezo1.ui.profile_plot import ProfilePlot, Trace
    plot = ProfilePlot()
    y = np.array([1.0, 2.0, np.nan, 4.0, np.nan])
    plot.set_data([Trace("t", np.arange(5.0), y)])
    assert np.isfinite(plot.left.lo) and np.isfinite(plot.left.hi)
    assert plot.left.hi >= 4.0


def test_empty_traces_are_dropped(qapp):
    from piezo1.ui.profile_plot import ProfilePlot, Trace
    plot = ProfilePlot()
    plot.set_data([Trace("empty", np.array([]), np.array([])),
                   Trace("real", np.arange(4.0), np.arange(4.0))])
    assert len(plot.traces) == 1


def test_click_reports_a_position_in_data_units(qapp):
    from piezo1.ui.profile_plot import ProfilePlot, Trace
    from PyQt6.QtCore import QPointF, Qt
    from PyQt6.QtGui import QMouseEvent

    plot = ProfilePlot()
    plot.resize(400, 240)
    plot.set_data([Trace("t", np.linspace(-50.0, 50.0, 50),
                         np.ones(50))])
    seen = []
    plot.position_clicked.connect(seen.append)

    rect = plot._plot_rect()
    centre = QPointF(rect.center().x(), rect.center().y())
    plot.mousePressEvent(QMouseEvent(
        QMouseEvent.Type.MouseButtonPress, centre, centre,
        Qt.MouseButton.LeftButton, Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier))
    assert seen and abs(seen[0]) < 5.0, seen


# --------------------------------------------------------------------------
# Panel plumbing
# --------------------------------------------------------------------------

class _FakeSlice:
    def __init__(self, z, radius, lining):
        self.z, self.radius = z, radius
        self.lining = tuple(lining)
        self.lining_names = tuple("ALA" for _ in lining)


class _FakeProfile:
    def __init__(self):
        self.z = np.linspace(-60.0, 10.0, 71)
        self.radius = 3.0 + np.cos(self.z / 12.0)
        self.slices = [_FakeSlice(z, r, [2440 + i])
                       for i, (z, r) in enumerate(zip(self.z, self.radius))]
        self.bottleneck_radius = float(self.radius.min())
        self.bottleneck_z = float(self.z[int(np.argmin(self.radius))])

    def constrictions(self, **kw):
        return self.slices[:3]


def test_panel_shows_the_pore_and_its_constrictions(qapp):
    from piezo1.ui.panels.analysis_panel import AnalysisPanel
    panel = AnalysisPanel()
    profile = _FakeProfile()
    panel.set_pore(profile, None, np.linspace(-0.4, 0.3, 71))
    assert "bottleneck" in panel.pore_label.text()
    assert panel.constriction_table.rowCount() == 3
    assert len(panel.pore_plot.traces) == 2, "hydrophobicity trace missing"


def test_hydrophobicity_can_be_switched_off(qapp):
    from piezo1.ui.panels.analysis_panel import AnalysisPanel
    panel = AnalysisPanel()
    panel.set_pore(_FakeProfile(), None, np.linspace(-0.4, 0.3, 71))
    assert len(panel.pore_plot.traces) == 2
    panel.hydro_check.setChecked(False)
    assert len(panel.pore_plot.traces) == 1


def test_selecting_a_constriction_emits_its_lining(qapp):
    from piezo1.ui.panels.analysis_panel import AnalysisPanel
    panel = AnalysisPanel()
    panel.set_pore(_FakeProfile(), None, None)
    seen = []
    panel.residues_selected.connect(lambda r, label: seen.append((list(r), label)))
    panel.constriction_table.selectRow(1)
    assert seen and seen[0][0] == [2441], seen


def test_scalar_menu_grows_as_analyses_finish(qapp):
    from piezo1.ui.panels.analysis_panel import AnalysisPanel
    panel = AnalysisPanel()
    assert panel.current_scalar() == ""
    panel.add_scalar("conservation", "Conservation")
    assert panel.current_scalar() == "conservation"
    panel.add_scalar("allostery", "Coupling")
    assert panel.current_scalar() == "allostery"
    panel.add_scalar("conservation", "Conservation")     # no duplicate
    assert panel.scalar_combo.count() == 3               # off + two scalars


def test_busy_disables_every_launch_button(qapp):
    from piezo1.ui.panels.analysis_panel import AnalysisPanel
    panel = AnalysisPanel()
    panel.set_busy(True, "pore")
    assert not panel.pore_button.isEnabled()
    assert not panel.pockets_button.isEnabled()
    panel.set_busy(False)
    assert panel.pore_button.isEnabled()


# --------------------------------------------------------------------------
# Residue values onto atoms
# --------------------------------------------------------------------------

def test_unmeasured_residues_take_the_floor_not_zero(human_structure, qapp):
    """Injecting zero for unmeasured residues rescales the whole legend.

    Conservation runs about 0.6-1.0. A zero for a residue that was never
    measured would stretch the colour scale over a range that contains no data,
    washing out every real difference — and it would look like a valid map.
    """
    from piezo1.ui.analysis_controller import AnalysisController

    class _Win:
        structure = human_structure
    controller = AnalysisController.__new__(AnalysisController)
    controller.win = _Win()

    values = {int(r): 0.6 + 0.4 * (i % 3) / 2.0
              for i, r in enumerate(np.unique(human_structure.res_seq)[:200])}
    atoms = controller.residue_values_to_atoms(values)
    assert atoms.shape == (human_structure.n_atoms,)
    assert atoms.min() >= min(values.values()) - 1e-6
    assert atoms.max() <= max(values.values()) + 1e-6


def test_residue_values_land_on_the_right_atoms(human_structure, qapp):
    from piezo1.ui.analysis_controller import AnalysisController

    class _Win:
        structure = human_structure
    controller = AnalysisController.__new__(AnalysisController)
    controller.win = _Win()

    target = int(np.unique(human_structure.res_seq)[50])
    atoms = controller.residue_values_to_atoms({target: 9.0, target + 1: 1.0})
    assert np.all(atoms[human_structure.res_seq == target] == 9.0)
    assert not np.any(atoms[human_structure.res_seq == target + 2] == 9.0)


# --------------------------------------------------------------------------
# Sessions
# --------------------------------------------------------------------------

def test_session_round_trip_preserves_view_state(tmp_path):
    session = Session(structure="8YEZ", species="human", style="cartoon",
                      color_by="domain", show_ligands=False, radius_scale=1.4,
                      camera_rotation=[0.5, 0.5, 0.5, 0.5],
                      camera_pivot=[1.0, 2.0, 3.0], camera_distance=275.0,
                      selected_residues=[2456, 2457],
                      selection_label="R2456H",
                      analyses={"pore": {"bottleneck_radius_A": 0.95}})
    path = tmp_path / "s.json"
    save_session(session, path)
    back = load_session(path)
    assert back.structure == "8YEZ"
    assert back.selected_residues == [2456, 2457]
    assert back.camera_distance == pytest.approx(275.0)
    assert back.radius_scale == pytest.approx(1.4)
    assert back.analyses["pore"]["bottleneck_radius_A"] == pytest.approx(0.95)


def test_session_records_parameters_not_results(tmp_path):
    """A session must not embed a result.

    Embedded numbers go stale the moment a parameter changes and would then be
    presented against a structure that never produced them.
    """
    session = Session(structure="8YEZ",
                      analyses={"modes": {"n_modes": 20}})
    path = tmp_path / "s.json"
    save_session(session, path)
    text = path.read_text()
    assert "eigenvalues" not in text
    assert "xyz" not in text and "coordinates" not in text


def test_structure_panel_restores_appearance(qapp):
    from piezo1.ui.panels.structure_panel import StructurePanel
    panel = StructurePanel()
    emitted = []
    panel.style_changed.connect(lambda v: emitted.append(("style", v)))
    panel.color_changed.connect(lambda v: emitted.append(("color", v)))

    panel.set_state(style="spheres", color_by="chain", ligands=False,
                    radius_scale=1.5)
    assert panel.current_style().value == "spheres"
    assert panel.current_color().value == "chain"
    assert not panel.ligand_check.isChecked()
    assert panel.radius_slider.value() == 150
    # One rebuild per property, not one per widget touched.
    assert len([e for e in emitted if e[0] == "style"]) == 1


# --------------------------------------------------------------------------
# The workers, run for real
# --------------------------------------------------------------------------

def test_pore_worker_matches_the_headless_result(human_structure, qapp):
    """The GUI must not be a second, divergent implementation.

    Everything the panel shows has to come from the same functions the CLI
    calls, so this asserts the worker reproduces them rather than that it
    produces something plausible.
    """
    from conftest import protomer_blocks
    from piezo1.analysis.hydration import load_grid, predict_wetting
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.superpose import detect_c3_axis
    from piezo1.ui.analysis_controller import AnalysisWorker

    blocks, _ = protomer_blocks(human_structure)
    result = AnalysisWorker("pore", {"structure": human_structure,
                                     "blocks": blocks})._pore()

    axis = detect_c3_axis(blocks)
    expected = pore_profile(human_structure, axis, step=1.0)
    assert result["profile"].bottleneck_radius == pytest.approx(
        expected.bottleneck_radius, rel=1e-9)

    grid = load_grid()
    if grid.available:
        direct = predict_wetting(human_structure, expected, grid)
        assert result["hydration"].score == pytest.approx(direct.score, rel=1e-9)
        assert result["hydrophobicity"] is not None


def test_allostery_worker_folds_the_three_protomers(human_structure, qapp):
    """One value per residue, not three.

    The trimer is C3 symmetric, so per-protomer differences are numerical
    noise. Leaving them unfolded colours one protomer differently from its
    mates, which reads as a rendering bug.
    """
    from conftest import protomer_blocks
    from piezo1.physics.anm import ANM
    from piezo1.ui.analysis_controller import AnalysisWorker

    blocks, residues = protomer_blocks(human_structure)
    modes = ANM.from_trimer(blocks, cutoff=15.0).build().calc_modes(n_modes=12)
    result = AnalysisWorker("allostery", {"modes": modes,
                                          "residues": residues})._allostery()
    assert len(result["values"]) == len(residues)
    assert all(np.isfinite(v) for v in result["values"].values())


def test_conservation_worker_drops_poorly_covered_positions(qapp):
    """Coverage below 0.7 measures the alignment, not selection pressure."""
    from piezo1.ui.analysis_controller import AnalysisWorker
    try:
        result = AnalysisWorker("conservation", {})._conservation()
    except Exception as exc:                     # no ortholog cache
        pytest.skip(f"conservation unavailable: {exc}")
    profile = result["profile"]
    assert len(result["values"]) < len(profile.residues)
    for residue in result["values"]:
        i = int(np.flatnonzero(profile.residues == residue)[0])
        assert profile.coverage[i] >= 0.7


def test_pockets_worker_honours_the_limit(human_structure, qapp):
    from piezo1.ui.analysis_controller import AnalysisWorker
    result = AnalysisWorker("pockets", {"structure": human_structure,
                                        "top": 4})._pockets()
    assert len(result["pockets"]) == 4
    volumes = [p.volume for p in result["pockets"]]
    assert volumes == sorted(volumes, reverse=True), "pockets must be ranked"


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
                                       "analysis", "measure"}
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
