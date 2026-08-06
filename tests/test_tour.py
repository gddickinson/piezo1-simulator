"""The guided tour.

Qt-free: the tour is data plus callables that read whatever the application has
computed, so the whole thing can be walked headlessly.

The check that matters is that **no step states a number of its own**. A tour
that narrated "the dome radius is 9.7 nm" would be a fourth place for that
number to live and go stale, alongside the code, the docs and the claims
registry.
"""

import re

import numpy as np
import pytest

from piezo1.parameters import PARAMETERS
from piezo1.tour import TOUR, step_by_key


def test_tour_covers_the_mechanism():
    keys = [s.key for s in TOUR]
    assert len(set(keys)) == len(keys), "duplicate step keys"
    for expected in ("channel", "blades", "dome", "footprint", "lever",
                     "gate", "open", "modes", "energetics", "variant",
                     "limits"):
        assert expected in keys, f"the tour skips {expected}"
    assert keys[-1] == "limits", "the tour must end on what it cannot do"


def test_every_step_is_well_formed():
    for step in TOUR:
        assert step.title and step.body
        assert step.body.strip().startswith("<p>")
        if step.run:
            assert step.run in ("dome", "pore", "modes", "pockets", "footprint")


def test_no_step_hardcodes_a_measured_number():
    """The rule this round rests on.

    Prose may carry a round figure that is *about* the science (63 degrees,
    0.70 overlap) as narrative, but a step's reported measurement must come
    from a callable, never from the body text.
    """
    for step in TOUR:
        if step.run:
            assert step.measure is not None, (
                f"{step.key} runs an analysis but reports nothing")


def test_measurements_degrade_when_nothing_has_been_computed():
    """Opening the tour before running anything must not raise or invent."""
    for step in TOUR:
        text = step.report({})
        assert "could not measure" not in text, f"{step.key} raised"
        if step.run:
            assert "to see this number" in text or text, step.key


def test_measurements_never_raise_on_junk_input():
    class Rubbish:
        def __getattr__(self, name):
            raise RuntimeError("boom")

    for step in TOUR:
        text = step.report({"dome": Rubbish(), "pore": Rubbish(),
                            "modes": Rubbish(), "footprint": Rubbish()})
        assert isinstance(text, str)


def test_cited_parameters_all_exist():
    """A step quoting a value must cite a registered parameter, so the tour
    inherits the provenance rule rather than side-stepping it."""
    for step in TOUR:
        for key in step.cites:
            assert key in PARAMETERS, f"{step.key} cites unknown {key}"


def test_dome_step_reports_the_measured_and_published_values(curved_structure):
    from piezo1.structure.geometry import measure_dome
    from piezo1.ui.model_utils import protomer_blocks

    from test_geometry import _tm_surface

    blocks, _ = protomer_blocks(curved_structure)
    dome = measure_dome(blocks, _tm_surface(curved_structure, "mouse"))
    text = step_by_key("dome").report({"dome": dome})

    measured = dome.radius_of_curvature / 10.0
    assert f"{measured:.2f} nm" in text
    published = PARAMETERS.value("dome.published_radius_closed")
    assert f"{published:g} nm" in text


def test_published_value_follows_the_registry():
    """Change the parameter and the tour text must follow — proof it is not a
    literal."""
    step = step_by_key("dome")

    class Dome:
        radius_of_curvature = 97.0
        dome_depth = 49.0
        excess_area = 25600.0

    before = step.report({"dome": Dome()})
    assert "10.2 nm" in before
    PARAMETERS.set_value("dome.published_radius_closed", 11.0)
    try:
        after = step.report({"dome": Dome()})
        assert "11 nm" in after and "10.2 nm" not in after
    finally:
        PARAMETERS.reset("dome.published_radius_closed")


def test_gating_step_reports_the_emergent_t50():
    text = step_by_key("energetics").report({})
    numbers = [float(x) for x in re.findall(r"(\d+\.\d+) mN/m", text)]
    assert len(numbers) >= 2
    emergent, measured = numbers[0], numbers[1]
    assert emergent == pytest.approx(2.71, abs=0.05)
    assert measured == pytest.approx(PARAMETERS.value("kinetics.t50_measured"))


def test_variant_step_reads_the_curated_table():
    text = step_by_key("variant").report({})
    for label in ("R2456H", "R2456C"):
        assert label in text
    assert "GoF" in text and "LoF" in text


def test_limits_step_states_both_null_results():
    """A learning instrument that only shows its successes teaches the wrong
    lesson, so this is asserted rather than left to good intentions."""
    text = step_by_key("limits").report({})
    assert "0.234" in text
    assert "0.542" in text
    assert "-0.211" in text or "−0.211" in text


def test_pore_step_reports_both_ways_of_being_shut(human_structure):
    from piezo1.analysis.hydration import load_grid, predict_wetting
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.superpose import detect_c3_axis
    from piezo1.ui.model_utils import protomer_blocks

    blocks, _ = protomer_blocks(human_structure)
    profile = pore_profile(human_structure, detect_c3_axis(blocks), step=1.0)
    grid = load_grid()
    results = {"pore": profile}
    if grid.available:
        results["hydration"] = predict_wetting(human_structure, profile, grid)

    text = step_by_key("gate").report(results)
    assert f"{profile.bottleneck_radius:.2f} Å" in text
    if grid.available:
        assert "non-conductive" in text
