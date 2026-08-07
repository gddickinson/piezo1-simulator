"""The scaffolding audit, and the two wrong versions of it that came first.

Round 55's job was to delete what does not earn its place. The hard part was
not the deleting — it was building a detector whose output could be acted on.
Two earlier versions would each have caused real damage:

- a grep over ``__all__`` reported **102** unused public names, including
  ``format_result`` (used inside its own module) and every return-type
  dataclass, which is constructed but never named elsewhere;
- an AST version that collapsed same-file references into a set reported
  **129** dead functions, including ``fetch_pdb``, ``cmd_list`` and
  ``_optimise_slice``.

Acting on either would have deleted the CLI. So every test here runs the
calibration first, per the rule in ``CLAUDE.md``.
"""

from __future__ import annotations

import pytest

from piezo1.dead_code import EXEMPT, DeadName, audit, calibration, reference_counts


@pytest.fixture(scope="module")
def calibrated():
    result = calibration()
    assert not result["false_positives"], (
        f"the detector flags known-used names {result['false_positives']}; "
        f"its output must not be acted on")
    assert result["detects_planted"], (
        "the detector cannot detect an unreferenced name at all")
    return result


def test_the_detector_is_calibrated(calibrated):
    assert calibrated["false_positives"] == []
    assert calibrated["detects_planted"] is True


def test_nothing_in_the_package_is_unreferenced(calibrated):
    """The ratchet. Round 55 took this to zero and it must stay there."""
    report = audit()
    assert report.n_definitions > 400, "the audit is not seeing the codebase"
    assert not report.dead, "\n".join(d.summary() for d in report.dead)


def test_same_file_references_count(calibrated):
    """The bug that made the second attempt report 129 dead functions.

    A helper called only from within its own module is live. Counting distinct
    files instead of occurrences hid exactly those.
    """
    counts = reference_counts()
    # `_optimise_slice` is called once, from inside structure/pore.py.
    internal = counts.get("_optimise_slice")
    assert internal is not None and sum(internal.values()) > 0


def test_names_reached_only_through_a_string_registry_count(calibrated):
    """The CLI dispatches by string; those functions are live.

    A detector that ignored string literals would report every `cmd_*` and
    every registry entry as dead.
    """
    counts = reference_counts()
    for name in ("permeation", "labelling", "interactions"):
        assert sum(counts.get(name, {}).values()) > 0, (
            f"{name} is an ANALYSES key and must count as referenced")


def test_entry_points_are_exempt_with_a_stated_reason(calibrated):
    """`main` is invoked by name from packaging metadata, not from code."""
    assert "main" in EXEMPT
    for name, reason in EXEMPT.items():
        assert len(reason) > 15, f"{name} is exempt without a real reason"


def test_the_report_names_what_it_found():
    dead = DeadName(name="f", kind="function", path="piezo1/x.py")
    text = dead.summary()
    assert "f" in text and "function" in text and "piezo1/x.py" in text


def test_deleted_names_are_really_gone():
    """The five things Round 55 removed, pinned so they cannot creep back."""
    import piezo1.analysis.design as design
    import piezo1.analysis.prediction_record as record
    import piezo1.physics.permeation as permeation
    import piezo1.render.colormaps as colormaps
    import piezo1.tour as tour

    assert not hasattr(design, "permutation_p")
    assert not hasattr(record, "PredictionContext")
    assert not hasattr(permeation, "_poisson_newton_step")
    assert not hasattr(colormaps, "distinct_colors")
    assert not hasattr(tour, "_published")


def test_the_permeation_finding_survived_its_functions_deletion():
    """`_poisson_newton_step` carried a real numerical result in its docstring.

    It was deleted only after checking the divergence it documented is recorded
    in the module, the test and SCIENCE.md. Deleting the last copy of a finding
    would be worse than keeping dead code.
    """
    from pathlib import Path

    root = Path(__file__).resolve().parents[1]
    for path in ("piezo1/physics/permeation.py", "tests/test_permeation.py",
                 "docs/SCIENCE.md"):
        assert "Gummel" in (root / path).read_text(), (
            f"the diverging-Gummel finding is no longer recorded in {path}")
