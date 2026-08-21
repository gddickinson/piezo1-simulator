"""The exploration catalogue: what each result window offers, and whether it
can actually produce it.

Qt-free by design — the catalogue, the chart builders and the simulations are
all plain data and plain functions, so the half of this feature that could be
quietly wrong is testable without a display. The window itself is in
``test_ui_explore.py``.

The tests are shaped around the three ways this could mislead rather than fail:

* a **figure** exhibit naming a command that does not build it, so a reader
  told to regenerate it never gets the picture back;
* a **simulation** writing to the parameter registry, which would leave the
  application quoting non-default numbers with nothing on screen to say so;
* a **chart** that raises on an ordinary result — a shut pore, an entry with no
  partner — instead of saying why there is nothing to draw.
"""

from __future__ import annotations

import math
from pathlib import Path

import pytest

from piezo1.parameters import PARAMETERS
from piezo1.ui.exhibit_models import SIMULATIONS, Context, run_simulation
from piezo1.ui.exhibit_plots import BUILDERS, build_chart
from piezo1.ui.exhibit_plots_family import FAMILY_BUILDERS
from piezo1.ui.exhibits import BASES, KINDS, ChartData, exhibits_for, registry

ROOT = Path(__file__).resolve().parents[1]

ALL_BUILDERS = dict(BUILDERS) | dict(FAMILY_BUILDERS)


def _all_exhibits():
    return [item for items in registry().values() for item in items]


# --------------------------------------------------------------------------
# The catalogue itself
# --------------------------------------------------------------------------

def test_the_catalogue_is_not_empty_and_covers_several_analyses():
    """A ratchet, so a refactor that quietly emptied the registry would fail
    rather than making every window offer nothing."""
    assert len(registry()) >= 20, sorted(registry())
    assert len(_all_exhibits()) >= 50


def test_every_exhibit_states_what_it_is_and_what_it_is_not():
    """The ``not_this`` line is the whole reason a picture is safe to show.

    A figure is more persuasive than the number behind it — that is why every
    overlay controller in this project carries a caveat — and it applies twice
    over to a curve the user has just produced by moving a slider.
    """
    for item in _all_exhibits():
        assert item.kind in KINDS
        assert item.basis in BASES
        assert len(item.what) > 30, f"{item.title}: no description"
        assert len(item.not_this) > 30, f"{item.title}: no 'not this' line"
        assert item.title.strip()


def test_each_kind_names_a_source_it_can_actually_reach():
    """A chart naming a builder that does not exist opens an empty panel, and
    nothing before the click would have said so."""
    for item in _all_exhibits():
        if item.kind == "chart":
            assert item.plot in ALL_BUILDERS, item.plot
        elif item.kind == "simulation":
            assert item.simulation in SIMULATIONS, item.simulation
        elif item.kind == "model":
            from piezo1.ui.explore_window import MODEL_ACTIONS

            assert item.action in MODEL_ACTIONS, item.action
        else:
            assert item.figure and item.rebuild


def test_no_structural_display_is_offered_by_nothing():
    """A display in the table that no result offers is a row nobody can
    reach — scaffolding, which is what ``dead_code.py`` exists to find."""
    from piezo1.ui.explore_window import MODEL_ACTIONS

    used = {item.action for item in _all_exhibits() if item.kind == "model"}
    assert set(MODEL_ACTIONS) == used, set(MODEL_ACTIONS) ^ used


def test_no_builder_or_simulation_is_orphaned():
    """The other direction: code nothing can reach is scaffolding, and
    ``dead_code.py`` exists because this project keeps growing some."""
    used_plots = {i.plot for i in _all_exhibits() if i.kind == "chart"}
    used_models = {i.simulation for i in _all_exhibits()
                   if i.kind == "simulation"}
    assert set(ALL_BUILDERS) == used_plots, set(ALL_BUILDERS) ^ used_plots
    assert set(SIMULATIONS) == used_models, set(SIMULATIONS) ^ used_models


# --------------------------------------------------------------------------
# Figures
# --------------------------------------------------------------------------

def test_every_figure_names_a_script_that_writes_it():
    """The rebuild command has to be true, not plausible.

    Checked by reading the named script and looking for the figure's own name
    in it — the cheap version of running it. Several figures in ``docs/img``
    have no producing script at all any more, and one of those in a catalogue
    would send a reader after a command that cannot bring it back.
    """
    figures = [i for i in _all_exhibits() if i.kind == "figure"]
    assert figures, "no figure exhibits to check"
    for item in figures:
        script = next((part for part in item.rebuild.split()
                       if part.endswith(".py")), "")
        assert script, f"{item.title}: no script in {item.rebuild!r}"
        path = ROOT / script
        assert path.exists(), f"{item.title}: {script} does not exist"
        stem = Path(item.figure).stem
        assert stem in path.read_text(), (
            f"{item.title}: {script} never mentions {stem}")


def test_a_figure_that_has_not_been_generated_resolves_to_nothing():
    """Missing figures are ordinary — ``docs/img`` is regenerable output — so
    the pane shows the command instead. The control is that a figure which
    *is* on disk resolves, or this would pass with a broken resolver."""
    from piezo1.ui.exhibits import Exhibit

    absent = Exhibit(analysis="x", kind="figure", title="t",
                     what="w" * 40, basis="measured", not_this="n" * 40,
                     figure="no_such_figure_at_all.png", rebuild="make it")
    assert absent.figure_file() is None

    present = [i for i in _all_exhibits()
               if i.kind == "figure" and i.figure_file() is not None]
    assert present, "no catalogued figure is on disk; the check is vacuous"


# --------------------------------------------------------------------------
# Charts
# --------------------------------------------------------------------------

#: The one builder that does not read the result at all: the modulator
#: potencies are formatted as prose in the result ("EC50 26.6 uM (...)"), and
#: parsing a sentence back into a number is how a plot ends up disagreeing
#: with the table beside it. It reads the curated resource both come from, so
#: it draws whatever it is handed.
FROM_A_RESOURCE = {"ligand_potency"}


def test_no_chart_builder_raises_on_an_empty_result():
    """Every builder must survive a result it cannot use, and say why.

    This is not hypothetical: ``permeation`` returns no conductance for the 17
    deposited entries the model calls shut, and that is the common case.
    """
    for name in ALL_BUILDERS:
        chart = build_chart(name, {})
        assert isinstance(chart, ChartData)
        if name in FROM_A_RESOURCE:
            assert not chart.empty, f"{name} reads a resource and drew nothing"
            continue
        assert chart.empty, f"{name} drew something out of an empty result"
        assert len(chart.note) > 20, f"{name} gave no reason: {chart.note!r}"


def test_an_unknown_chart_is_reported_rather_than_raised():
    chart = build_chart("no_such_chart", {})
    assert chart.empty and "no_such_chart" in chart.note


def test_a_builder_that_raises_is_caught_and_named():
    """The dispatcher's promise, checked by making one fail on purpose."""
    chart = build_chart("interaction_counts", {"counts": {"a": "not a number"}})
    assert chart.empty
    assert "could not draw" in chart.note or chart.note


def test_a_chart_built_from_a_result_reports_that_result():
    """Nothing is recomputed: the numbers on the chart are the ones handed in."""
    data = {"counts": {"hydrogen_bond": 5619, "disulfide": 3}}
    chart = build_chart("interaction_counts", data)
    assert not chart.empty
    assert chart.series[0].y == [5619.0, 3.0]
    assert chart.log_y, "three orders of magnitude on a linear axis is one bar"


# --------------------------------------------------------------------------
# Simulations
# --------------------------------------------------------------------------

def _headless_simulations():
    return [s for s in SIMULATIONS.values() if not s.needs_structure]


def test_every_simulation_runs_across_the_whole_range_of_its_controls():
    """Moved to either end, every control must still produce a finite curve.

    The corners are where a model breaks — a buffer that screens to nothing, a
    Hill coefficient of a half — and the user is entitled to walk into one.
    """
    for simulation in _headless_simulations():
        base = {c.key: c.start() for c in simulation.controls}
        for control in simulation.controls:
            for value in (control.low, control.start(), control.high):
                chart = run_simulation(simulation.key,
                                       dict(base, **{control.key: value}),
                                       Context())
                assert not chart.empty, (
                    f"{simulation.key} at {control.key}={value}: {chart.note}")
                for series in chart.series:
                    finite = [v for v in series.y
                              if v is not None and math.isfinite(v)]
                    assert finite, (f"{simulation.key} at {control.key}="
                                    f"{value}: {series.name} is all non-finite")


def test_a_simulation_moves_when_its_control_moves():
    """A slider that changes nothing is a lie in a window.

    The calibration is the other half: the same inputs twice must give the
    identical curve, or "it moved" would be no evidence at all.
    """
    for simulation in _headless_simulations():
        base = {c.key: c.start() for c in simulation.controls}
        first = run_simulation(simulation.key, base, Context())
        again = run_simulation(simulation.key, dict(base), Context())
        assert first.series[0].y == again.series[0].y, simulation.key

        moved = False
        for control in simulation.controls:
            other = dict(base, **{control.key: control.high})
            changed = run_simulation(simulation.key, other, Context())
            if any(a.y != b.y for a, b in zip(first.series, changed.series)):
                moved = True
        assert moved, f"no control changes {simulation.key}"


def test_no_simulation_writes_to_the_parameter_registry():
    """The rule this feature would do real damage by breaking.

    An override survives the window: reports carry the amber banner and
    ``verify_claims`` refuses to run against a modified registry, which is
    exactly right and would be baffling if a slider had caused it.

    Calibrated: the check is shown to notice an override before it is trusted.
    """
    assert not PARAMETERS.modified, "the registry was already modified"
    PARAMETERS.set_value("dome.delta_area", 12.0)
    assert PARAMETERS.modified, "the guard cannot see an override at all"
    PARAMETERS.reset()
    assert not PARAMETERS.modified

    for simulation in _headless_simulations():
        base = {c.key: c.start() for c in simulation.controls}
        for control in simulation.controls:
            run_simulation(simulation.key,
                           dict(base, **{control.key: control.high}), Context())
    assert not PARAMETERS.modified, PARAMETERS.override_summary()


def test_a_simulation_needing_coordinates_says_so_rather_than_guessing():
    needs = [s for s in SIMULATIONS.values() if s.needs_structure]
    assert needs, "no simulation reads coordinates; the guard is vacuous"
    for simulation in needs:
        values = {c.key: c.start() for c in simulation.controls}
        chart = run_simulation(simulation.key, values, Context(structure=None))
        assert chart.empty and "structure" in chart.note


def test_an_unknown_simulation_is_reported_rather_than_raised():
    chart = run_simulation("no_such_model", {}, Context())
    assert chart.empty and "no_such_model" in chart.note


# --------------------------------------------------------------------------
# Controls that are registered parameters
# --------------------------------------------------------------------------

def test_a_registered_control_starts_at_the_registry_value():
    controls = [c for s in SIMULATIONS.values() for c in s.controls
                if c.parameter]
    assert controls, "no control is tied to a registered parameter"
    for control in controls:
        expected = PARAMETERS.value(control.parameter) * control.scale
        assert control.start() == pytest.approx(expected), control.key


def test_a_registered_control_follows_an_override_at_call_time():
    """Resolved when the panel is built, not when the module is imported —
    the rule every dataclass default in this project follows, and the reason
    the unit conversion is a factor rather than a literal default."""
    control = next(c for s in SIMULATIONS.values() for c in s.controls
                   if c.parameter == "dome.delta_area")
    before = control.start()
    try:
        PARAMETERS.set_value("dome.delta_area", before + 3.0)
        assert control.start() == pytest.approx(before + 3.0)
    finally:
        PARAMETERS.reset()
    assert control.start() == pytest.approx(before)


def test_a_control_range_lies_inside_the_registered_bounds():
    """A slider that can leave the parameter's declared bounds would let the
    user explore a value the registry itself would refuse."""
    for simulation in SIMULATIONS.values():
        for control in simulation.controls:
            if not control.parameter:
                continue
            parameter = PARAMETERS.parameters[control.parameter]
            assert control.low * (1 / control.scale) >= parameter.minimum - 1e-12, control.key
            assert control.high * (1 / control.scale) <= parameter.maximum + 1e-12, control.key


# --------------------------------------------------------------------------
# Coverage of the windows themselves is in test_ui_explore.py, which needs Qt
# to read the caveat table. What can be checked here is that the catalogue
# names analyses the shared registry knows about.
# --------------------------------------------------------------------------

def test_every_analysis_explored_is_one_the_project_can_run():
    from piezo1.analysis.report import ANALYSES

    extra = {"variant_structures"}          # built in the GUI, not in ANALYSES
    unknown = set(registry()) - set(ANALYSES) - extra
    assert not unknown, f"exhibits for analyses that do not exist: {unknown}"


def test_the_family_analyses_are_explored_as_imported_results():
    """The census entries must not present as measurements. Their exhibits
    carry the imported basis, which is what the window prints beside them."""
    for key in ("family", "constraint"):
        bases = {item.basis for item in exhibits_for(key)}
        assert "imported" in bases, (key, bases)
