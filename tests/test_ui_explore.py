"""The Explore button, and the window behind it.

Every analysis window in this application is a table, and the reasoning behind
the table — the figure, the model the number came out of, the same result drawn
on the structure — used to live in ``docs/img`` and in scripts. These check the
button that closes that gap, on the offscreen Qt platform with the real
widgets, because the two Qt refactors that silently broke this GUI both looked
correct in source.

What is deliberately checked **in pixels** is that a chart reaches the screen.
Everything in this project that draws has at some point passed every test while
drawing nothing at all — the cylinder impostors for the whole life of the
renderer, the ion stream until Round 84c — and a chart is no different: it
uploads nothing and raises nothing when it paints an empty picture. Writing
these found exactly that, in bars on a logarithmic axis.
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")

from PyQt6.QtCore import QSize  # noqa: E402
from PyQt6.QtGui import QColor, QImage, QPainter  # noqa: E402
from PyQt6.QtWidgets import QApplication, QLabel  # noqa: E402

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.parameters import PARAMETERS  # noqa: E402
from piezo1.ui.exhibits import ChartData, Series, empty_chart, registry  # noqa: E402


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
def window(qapp):
    """A main window with menus and panels, but no GL scene.

    Enough for everything except actually drawing on the model: the exhibit
    list, the panes and the controls a model action presses all exist without
    a context.
    """
    from piezo1.ui.main_window import MainWindow

    win = MainWindow()
    yield win
    for controller in ("analysis", "physics", "overlay"):
        cleanup = getattr(getattr(win, controller, None), "cleanup", None)
        if cleanup is not None:
            cleanup()
    qapp.processEvents()


def lit_pixels(widget, colour: str, size=(560, 360), tolerance: int = 40) -> int:
    """How many pixels of a series colour the widget actually paints."""
    want = QColor(colour)
    widget.resize(*size)
    image = QImage(QSize(*size), QImage.Format.Format_ARGB32)
    image.fill(0)
    painter = QPainter(image)
    widget.render(painter)
    painter.end()
    count = 0
    for y in range(size[1]):
        for x in range(size[0]):
            pixel = image.pixelColor(x, y)
            if (pixel.alpha() > 200
                    and abs(pixel.red() - want.red()) < tolerance
                    and abs(pixel.green() - want.green()) < tolerance
                    and abs(pixel.blue() - want.blue()) < tolerance):
                count += 1
    return count


def labels(widget) -> list[str]:
    return [child.text() for child in widget.findChildren(QLabel)]


# --------------------------------------------------------------------------
# The button on every pop-up
# --------------------------------------------------------------------------

def test_every_analysis_window_offers_something_to_explore(window):
    """The requirement, stated as a guard: a result window with no exhibits
    would show a dead button, and a new analysis added without any would be
    invisible until someone clicked it."""
    from piezo1.ui.tabular_analyses import CAVEATS

    missing = [key for key in CAVEATS if not registry().get(key)]
    assert not missing, f"result windows with nothing to explore: {missing}"


def test_the_button_is_wired_on_a_window_the_application_opens(window, qapp):
    """Driven through the real path rather than by constructing a dialog:
    ``_show_result`` is where the callback and the exhibit count are bound."""
    window._show_result("interactions", "Interactions",
                        {"counts": {"salt_bridge": 4}})
    qapp.processEvents()
    dialog = window._result_dialogs[-1]
    assert dialog.explore_button.isEnabled()
    dialog.close()


def test_a_result_with_nothing_registered_disables_the_button(qapp):
    """The calibration: the button can be off, so its being on means
    something."""
    from piezo1.ui.result_dialog import ResultDialog

    dialog = ResultDialog("nothing", {"a": 1}, explore=None, n_exhibits=0)
    assert not dialog.explore_button.isEnabled()
    assert "Nothing is registered" in dialog.explore_button.toolTip()
    dialog.close()


def test_the_button_opens_a_window_for_that_analysis(window, qapp):
    window._show_result("ligands", "Modulators", {})
    qapp.processEvents()
    dialog = window._result_dialogs[-1]
    dialog.explore_button.click()
    qapp.processEvents()
    opened = window._explore_windows[-1]
    assert opened.analysis == "ligands"
    assert opened.list.count() == len(registry()["ligands"])
    opened.close()
    dialog.close()


def test_the_provenance_stamp_travels_to_the_explore_window(window, qapp):
    """Both windows are non-modal and the registry can move under either, so
    the stamp is passed across rather than read again."""
    window._show_result("interactions", "Interactions", {"counts": {"x": 1}})
    qapp.processEvents()
    dialog = window._result_dialogs[-1]
    dialog.explore_button.click()
    qapp.processEvents()
    opened = window._explore_windows[-1]
    assert dialog.provenance in labels(opened)
    opened.close()
    dialog.close()


# --------------------------------------------------------------------------
# Every exhibit, selected
# --------------------------------------------------------------------------

def test_every_exhibit_selects_and_carries_its_warning(window, qapp):
    """Each exhibit is opened with an empty result — the case a user meets
    when an analysis returned little — and each must still show its pane and
    the line saying what it is not."""
    from piezo1.ui.explore_window import ExploreWindow

    for analysis, items in registry().items():
        opened = ExploreWindow(analysis, analysis, {}, window=window,
                               provenance="stamp")
        opened.show()
        for row in range(len(items)):
            opened.list.setCurrentRow(row)
            qapp.processEvents()
            warnings = [t for t in labels(opened) if t.startswith("NOT THIS:")]
            assert warnings, f"{analysis}/{items[row].title}: no warning shown"
            assert items[row].not_this in warnings[0]
            assert opened._pane is not None and opened._pane.isVisible()
        opened.close()


def test_switching_exhibits_leaves_nothing_of_the_last_one_behind(window,
                                                                  qapp):
    """The defect this was written after: the panel is rebuilt on every
    selection, and a *nested layout* taken out of a layout leaves its widgets
    parented and painting. What the user saw was the previous exhibit's chip
    hanging over a blank panel — while every test passed, because
    ``findChildren`` finds a widget wherever it is.

    Checked in both directions: the new exhibit's text present, the old
    exhibit's text gone.
    """
    from piezo1.ui.explore_window import ExploreWindow

    items = registry()["permeation"]
    assert len(items) >= 2, "needs two exhibits to switch between"
    opened = ExploreWindow("permeation", "Permeation", {}, window=window,
                           provenance="stamp")
    opened.show()
    opened.list.setCurrentRow(0)
    qapp.processEvents()
    assert items[0].what in labels(opened)

    opened.list.setCurrentRow(1)
    qapp.processEvents()
    shown = labels(opened)
    assert items[1].what in shown
    assert items[0].what not in shown, "the previous exhibit is still there"
    assert items[0].not_this not in " ".join(shown)
    opened.close()


# --------------------------------------------------------------------------
# The chart actually paints
# --------------------------------------------------------------------------

def test_a_rebuilt_panel_is_laid_out_and_not_just_added(window, qapp):
    """Selecting a second exhibit rebuilds the panel of a visible window.

    The widgets must end up laid out down the panel rather than stacked at the
    origin. **This guard is weaker than it looks and the next test says why**:
    the offscreen platform lays a rebuilt panel out by itself, so it passes
    either way. It is kept for gross breakage, not as the calibration.
    """
    from piezo1.ui.explore_window import ExploreWindow

    opened = ExploreWindow("permeation", "Permeation", {}, window=window,
                           provenance="stamp")
    opened.resize(900, 640)
    opened.show()
    opened.list.setCurrentRow(1)
    qapp.processEvents()

    tops = []
    for index in range(opened.content_layout.count()):
        widget = opened.content_layout.itemAt(index).widget()
        assert widget is not None and widget.isVisible()
        tops.append(widget.geometry().top())
    assert len(tops) >= 4
    assert tops == sorted(tops) and len(set(tops)) == len(tops), (
        f"the panel was never laid out: {tops}")
    opened.close()


def test_the_rebuild_shows_what_it_adds(window):
    """A discipline checked in the source, because no test on this platform
    can check it in behaviour.

    On macOS the rebuilt panel came up **blank** in the running application:
    the new widgets were in the layout, unshown, each at its default 640x480
    at the origin. Offscreen — where this suite runs — Qt lays the same panel
    out correctly, so the behavioural guard above passes on the broken code.
    Counted in the source instead, the way ``test_martini`` counts the one
    place a run may be constructed: it is a discipline, not a type.
    """
    import inspect

    from piezo1.ui.explore_window import ExploreWindow

    source = inspect.getsource(ExploreWindow._select)
    assert ".show()" in source, "the rebuild does not show what it adds"
    assert "activate()" in source, "the rebuild never activates the layout"


def test_the_panel_clears_a_nested_layout_too(window, qapp):
    """The calibration for the test above, with the bug planted.

    A widget inside a *nested* layout is the case that survived: taking the
    layout out of its parent leaves the widget parented to the panel, still
    painting where it was. Planting one and watching it disappear is what
    makes the guard above evidence rather than a coincidence of the current
    layout being flat.
    """
    from PyQt6.QtWidgets import QHBoxLayout
    from piezo1.ui.explore_window import ExploreWindow

    opened = ExploreWindow("permeation", "Permeation", {}, window=window,
                           provenance="stamp")
    opened.show()
    opened.list.setCurrentRow(0)
    qapp.processEvents()

    nested = QHBoxLayout()
    nested.addWidget(QLabel("STRAY LEFTOVER"))
    opened.content_layout.addLayout(nested)
    qapp.processEvents()
    assert "STRAY LEFTOVER" in labels(opened), "the planted widget never showed"

    opened.list.setCurrentRow(1)
    qapp.processEvents()
    assert "STRAY LEFTOVER" not in labels(opened), (
        "a widget in a nested layout survived the rebuild")
    opened.close()


def test_a_chart_reaches_the_screen(qapp):
    from piezo1.ui.exhibit_chart import ChartView

    view = ChartView()
    view.set_chart(empty_chart("nothing to draw"))
    assert lit_pixels(view, "#6fb1ff") == 0, "an empty chart painted a series"

    view.set_chart(ChartData(
        y_label="y", series=[Series("line", list(range(40)),
                                    [float(i * i) for i in range(40)])]))
    assert lit_pixels(view, "#6fb1ff") > 100, "the line never reached the screen"


def test_bars_on_a_logarithmic_axis_are_drawn(qapp):
    """Pinned because it was broken, and invisibly: a bar is read from zero,
    log(0) is not a number, and the whole series painted nothing while the
    axis, the grid and the labels all came out perfectly."""
    from piezo1.ui.exhibit_chart import ChartView

    view = ChartView()
    chart = ChartData(y_label="count", categories=["a", "b", "c"], log_y=True,
                      series=[Series("bars", [0, 1, 2], [3.0, 900.0, 27000.0],
                                     kind="bar", color="#6fb1ff")])
    view.set_chart(chart)
    assert lit_pixels(view, "#6fb1ff") > 500

    linear = ChartData(y_label="count", categories=["a", "b", "c"],
                       series=[Series("bars", [0, 1, 2], [3.0, 9.0, 27.0],
                                      kind="bar", color="#6fb1ff")])
    view.set_chart(linear)
    assert lit_pixels(view, "#6fb1ff") > 500


def test_a_reference_band_is_drawn_over_the_data_not_under_it(qapp):
    """A single bar covers the whole plot area behind it, and the band is the
    thing the bar has to be read against."""
    from piezo1.ui.exhibit_chart import ChartView
    from piezo1.ui.exhibits import Reference

    view = ChartView()
    view.set_chart(ChartData(
        y_label="A", categories=["8YFG"],
        series=[Series("variant", [0], [0.81], kind="bar", color="#e06c75")],
        references=[Reference(0.67, "wild-type range", high=0.93)]))
    # The dashed reference ink survives on top of the bar.
    assert lit_pixels(view, "#8a919e", tolerance=60) > 30


# --------------------------------------------------------------------------
# Charts built from results the analyses really produce
# --------------------------------------------------------------------------

CHEAP = {"interactions": None, "labelling": None, "fusion": None,
         "hybrid": None, "ligands": None, "prediction_record": None,
         "family": None, "constraint": None, "disease": None,
         "nanodomain": None, "fluctuations": None}


def test_charts_built_from_real_results_are_not_empty(qapp):
    """Synthetic dicts drift from what the analyses return; these are the real
    thing.

    Only the analyses that run in a second or two are here. The expensive ones
    — homology at a minute, liu2025 at half of one, piezo3, paired_variant,
    paralogue, permeation — are named rather than silently skipped, and their
    builders are exercised on empty input in ``test_exhibits.py``.
    """
    from piezo1.analysis.report import ANALYSES
    from piezo1.core import Structure
    from piezo1.ui.exhibit_plots import build_chart

    path = STRUCTURE_DIR / "7WLT.cif"
    if not path.exists():
        pytest.skip("7WLT.cif not downloaded — run python -m piezo1.io.fetch")
    structure = Structure.from_file(path)

    checked = 0
    for analysis in CHEAP:
        data = ANALYSES[analysis](structure, "mouse")
        for item in registry()[analysis]:
            if item.kind != "chart":
                continue
            chart = build_chart(item.plot, data)
            assert not chart.empty, f"{analysis}/{item.plot}: {chart.note}"
            checked += 1
    assert checked >= 10, checked


# --------------------------------------------------------------------------
# Simulations in the window
# --------------------------------------------------------------------------

def _pane_for(simulation_key: str):
    from piezo1.ui.exhibit_models import SIMULATIONS, Context
    from piezo1.ui.exhibit_panes import SimulationPane

    return SimulationPane(SIMULATIONS[simulation_key], Context())


def test_moving_a_slider_recomputes_the_curve(qapp):
    from PyQt6.QtTest import QTest

    pane = _pane_for("dome_activation")
    before = list(pane.view.chart.series[0].y)
    slider = pane._sliders["delta_area"]
    slider.setValue(min(slider.value() + 120, slider.maximum()))
    QTest.qWait(pane.DEBOUNCE_MS * 4)
    qapp.processEvents()
    after = list(pane.view.chart.series[0].y)
    assert before != after, "the curve did not follow the slider"
    assert pane.view.chart.note, "the model said nothing about what it drew"


def test_the_pane_marks_which_controls_are_at_their_registered_default(qapp):
    pane = _pane_for("dome_activation")
    assert any("default" in label.text()
               for label in pane._readouts.values()), \
        "nothing says the sliders start at the registry's own values"


def test_driving_the_pane_leaves_the_registry_alone(qapp):
    """The same rule as the headless test, checked through the widget — which
    is the path a user takes."""
    from PyQt6.QtTest import QTest

    assert not PARAMETERS.modified
    pane = _pane_for("labelling_timecourse")
    for slider in pane._sliders.values():
        slider.setValue(slider.maximum())
        QTest.qWait(pane.DEBOUNCE_MS * 2)
    qapp.processEvents()
    assert not PARAMETERS.modified, PARAMETERS.override_summary()


def test_the_simulation_pane_cannot_be_shown_without_its_caveat(qapp):
    """A curve the user has just produced looks exactly like a measured one."""
    pane = _pane_for("dose_response")
    text = " ".join(labels(pane))
    assert "SENSITIVITY, NOT A MEASUREMENT" in text
    assert "registry is not written to" in text


# --------------------------------------------------------------------------
# Figures and model actions
# --------------------------------------------------------------------------

def test_a_missing_figure_shows_the_command_that_builds_it(qapp):
    from piezo1.ui.exhibit_panes import FigurePane
    from piezo1.ui.exhibits import Exhibit

    exhibit = Exhibit(analysis="x", kind="figure", title="t", what="w" * 40,
                      basis="measured", not_this="n" * 40,
                      figure="not_generated_here.png",
                      rebuild="python scripts/make_figures.py")
    pane = FigurePane(exhibit)
    text = " ".join(labels(pane))
    assert "python scripts/make_figures.py" in text
    assert "not been generated" in text


def test_every_model_action_presses_a_control_the_user_could_press(window):
    """The one-control rule. If this table drifts — a menu entry renamed, a
    panel button moved — the exhibit would silently do nothing, and the only
    symptom would be a button that appears to work."""
    from piezo1.ui.explore_window import MODEL_ACTIONS, BoundAction

    for key, spec in MODEL_ACTIONS.items():
        bound = BoundAction(spec, window)
        assert bound.resolved, f"{key} resolves to no control ({spec})"


def test_a_model_action_does_not_turn_an_overlay_back_off(window, qapp):
    """Pressing "draw this" twice must not undo it, which a naive trigger()
    on a checkable action does."""
    from piezo1.ui.explore_window import MODEL_ACTIONS, BoundAction

    bound = BoundAction(MODEL_ACTIONS["contacts"], window)
    assert bound.action is not None and bound.action.isCheckable()
    bound.action.setChecked(True)
    message = bound.run()
    assert bound.action.isChecked(), "the second press switched it off"
    assert "already drawn" in message
    bound.action.setChecked(False)


def test_an_unavailable_control_says_so_rather_than_doing_nothing(window):
    from piezo1.ui.explore_window import MODEL_ACTIONS, BoundAction

    bound = BoundAction(MODEL_ACTIONS["colour_fluctuation"], window)
    assert bound.button is not None
    bound.button.setEnabled(False)
    assert "not available" in bound.run()
