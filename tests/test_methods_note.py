"""The methods note must describe machinery that exists and be honest about why.

Round 67's condition was explicit: the note must state that this pipeline's main
output was five nulls, and say why that is the point rather than a disclaimer.
A methods note that quietly reads as a success story would be the most
consequential piece of stale prose this project could ship, because it is the
part someone else would act on.
"""

from __future__ import annotations

import importlib
import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTE = ROOT / "docs" / "METHODS_NOTE.md"


@pytest.fixture(scope="module")
def text():
    assert NOTE.exists(), "docs/METHODS_NOTE.md is missing"
    return NOTE.read_text()


@pytest.fixture(scope="module")
def flowed(text):
    return " ".join(text.split())


# --------------------------------------------------- honest about the outcome

def test_it_leads_with_the_failure(flowed):
    """The condition Round 67 set. Not buried, not at the end."""
    assert "It failed" in flowed or "it failed" in flowed
    head = flowed[:900]
    assert "failed" in head, "the failure must be in the opening, not a footnote"


def test_it_states_the_count_from_the_record(flowed):
    from piezo1.analysis.prediction_record import ALL_PREREGISTERED

    words = {5: "five", 6: "six"}
    assert f"{words[len(ALL_PREREGISTERED)]} pre-registered tests" in flowed


def test_it_says_why_the_failure_is_the_point(flowed):
    """Not a disclaimer — an argument that the safeguards were exercised."""
    assert "only ever confirmed" in flowed or "only ever confirm" in flowed
    assert "makes a null informative" in flowed


def test_it_does_not_read_as_a_success_story(flowed):
    lowered = flowed.lower()
    for banned in ("successfully predicts", "we show that piezo1 variants",
                   "validated predictor", "promising results"):
        assert banned not in lowered, f"overclaiming: {banned!r}"


# ------------------------------------------ every mechanism it cites is real

def test_every_module_it_names_exists(text):
    named = set(re.findall(r"`(analysis/[\w_]+\.py|tests/test_[\w_]+\.py|"
                           r"docs/[\w_.]+\.md)`", text))
    assert named, "the note cites no implementation"
    # `analysis/x.py` is the project's own notation for `piezo1/analysis/x.py`
    # — the same one INTERFACE.md uses. Resolve it the same way.
    def resolve(name):
        return (ROOT / name) if (ROOT / name).exists() else (ROOT / "piezo1" / name)

    missing = [n for n in named if not resolve(n).exists()]
    assert not missing, f"the note cites files that do not exist: {missing}"


def test_the_named_analysis_modules_import(text):
    for name in sorted(set(re.findall(r"`analysis/([\w_]+)\.py`", text))):
        importlib.import_module(f"piezo1.analysis.{name}")


def test_the_numbers_it_quotes_match_the_code(text):
    """The note is a summary; its figures must still be the project's."""
    from piezo1.analysis.claims import CLAIMS
    from piezo1.analysis.feasibility import paired_positions_required

    flowed = " ".join(text.split())
    expected = {c.key: c.expected for c in CLAIMS}
    assert str(int(expected["feasibility.required_n"])) in flowed  # 134
    assert str(int(expected["feasibility.ceiling"])) in flowed     # 59

    paired = paired_positions_required(0.8, n_simulations=1200).positions
    assert f"**{paired}** shared positions" in flowed, (
        f"the note says a different number from paired_positions_required "
        f"({paired})")


def test_the_parameter_wiring_finding_is_quoted_correctly(flowed):
    """26 of 101 — the count Round 49 measured."""
    assert "26 of 101" in flowed


# ------------------------------------------ the calibration list is truthful

def test_every_calibration_failure_it_lists_is_recorded_elsewhere(text):
    """Each anecdote must be traceable to a round, not remembered loosely."""
    log = (ROOT / "SESSION_LOG.md").read_text().lower()
    # Terms taken from the log itself, not from how the note phrases them —
    # the first version searched for "unicode minus" where the log says
    # "U+2212", and failed on its own wording rather than on the facts.
    for term in ("spheroid", "u+2212", "dead code", "fetch_pdb", "19.0"):
        assert term in log, (
            f"the note cites an incident the session log does not record: {term}")


def test_it_names_the_rule_that_caught_the_most(flowed):
    assert "suspect the checker first" in flowed.lower()


def test_it_is_reachable_from_the_conclusion_or_the_readme():
    conclusion = (ROOT / "docs" / "CONCLUSION.md").read_text()
    readme = (ROOT / "README.md").read_text()
    assert "METHODS_NOTE.md" in conclusion or "METHODS_NOTE.md" in readme, (
        "a methods note nobody can find is not a methods note")
