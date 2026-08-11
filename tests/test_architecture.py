"""`ARCHITECTURE.md` must explain *why*, not restate `INTERFACE.md`.

Round 79's validation clause, and the reason for it: the row sat at 📋 for the
whole project, and the easy way to close it would have been to write a summary
of the module layout. That document already exists. A second copy of it would
go stale in a different way from the first and settle nothing.

So this checks the document is about **constraints and the incidents that
forced them** — and, more usefully, that the architecture it describes is still
the one the code has. A `why` document that describes a rule the code stopped
following is worse than none.
"""

from __future__ import annotations

import re
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
DOC = ROOT / "docs" / "ARCHITECTURE.md"


@pytest.fixture(scope="module")
def text() -> str:
    assert DOC.exists(), "ARCHITECTURE.md is promised in INTERFACE.md"
    return DOC.read_text()


# ------------------------------------------------- it is the why, not the what

def test_it_explains_rather_than_lists(text):
    """Every section has to give a reason, not an inventory."""
    lowered = text.lower()
    for word in ("constraint", "why", "costs"):
        assert word in lowered, f"no {word!r} anywhere; this reads as a summary"
    assert lowered.count("constraint") >= 5, (
        "fewer than five sections state what forced the decision")


def test_it_does_not_restate_the_navigation_map(text):
    """The specific failure mode: a second INTERFACE.md.

    A per-module table is the giveaway, so the count of module-shaped rows is
    capped well below the number of modules the project has.
    """
    rows = re.findall(r"^\|\s*`[\w./]+\.py`", text, re.M)
    assert len(rows) < 5, (
        f"{len(rows)} module rows: this is turning into INTERFACE.md, which "
        f"already exists and is generated from the same facts")


def test_it_names_the_incidents_rather_than_asserting_the_rules(text):
    """A rule whose reason is forgotten is the next thing to be undone.

    Each of these is a real, recorded failure; citing them is what makes the
    document worth keeping.
    """
    lowered = " ".join(text.split()).lower()
    for incident in ("model_utils", "127-byte", "26 registered parameters",
                     "non-modal", "spheroid"):
        assert incident in lowered, f"the {incident!r} incident is not cited"


# ----------------------------------- the architecture it describes still holds

def test_the_dependency_arrow_it_claims_is_the_one_the_code_has():
    """The central claim, checked against imports rather than trusted."""
    offenders = []
    for layer in ("core", "structure", "physics", "analysis", "io"):
        for path in (ROOT / "piezo1" / layer).rglob("*.py"):
            source = path.read_text()
            for banned in ("from ..ui", "from ..render", "import piezo1.ui",
                           "from piezo1.ui"):
                if banned in source:
                    offenders.append(f"{path.relative_to(ROOT)}: {banned}")
    assert not offenders, (
        "ARCHITECTURE.md claims the arrow points one way, but: " + "; ".join(offenders))


def test_the_structure_of_arrays_claim_still_holds():
    """`Structure` must still be parallel arrays, not a list of atoms."""
    import numpy as np

    from piezo1.core.structure import Structure

    fields = Structure.__dataclass_fields__
    for name in ("xyz", "element", "res_seq", "chain"):
        assert name in fields, f"Structure lost its {name} array"
    from piezo1.config import STRUCTURE_DIR
    path = STRUCTURE_DIR / "8YEZ.cif"
    if not path.exists():
        pytest.skip("8YEZ not downloaded; run python -m piezo1.io.fetch")
    st = Structure.from_file(path)
    assert isinstance(st.xyz, np.ndarray) and st.xyz.ndim == 2
    assert isinstance(st.mask_ca(), np.ndarray)
    assert st.mask_ca().dtype == bool, "selections must be boolean masks"


def test_the_calibration_rule_it_describes_is_enforced_somewhere():
    """The document calls this the most expensive lesson; a document is not
    enforcement, so point at the thing that is."""
    from tests.test_calibration import CALIBRATED

    assert CALIBRATED, "the calibration register is empty"
    assert (ROOT / "tests" / "test_calibration.py").exists()


def test_every_document_it_points_at_exists(text):
    missing = []
    for target in re.findall(r"\[[^\]]+\]\(([^)]+)\)", text):
        if target.startswith(("http", "#")):
            continue
        if not (DOC.parent / target).exists() and not (ROOT / target).exists():
            missing.append(target)
    assert not missing, f"ARCHITECTURE.md links to nothing: {missing}"


def test_interface_no_longer_promises_it_as_planned():
    """The row was 📋 for the whole project. Either written or gone."""
    interface = (ROOT / "INTERFACE.md").read_text()
    # Match the table *row* for the file, not any line that mentions it — the
    # row describing this very test file names ARCHITECTURE.md and legitimately
    # contains a 📋 character, which the looser check tripped over.
    rows = [line for line in interface.splitlines()
            if line.startswith("| `ARCHITECTURE.md`")]
    assert rows, "INTERFACE.md no longer has a row for ARCHITECTURE.md"
    for row in rows:
        assert "📋" not in row, (
            "ARCHITECTURE.md exists but INTERFACE still marks it planned")


def test_no_planned_documents_remain_unwritten():
    """The ratchet Round 79 closes: 📋 rows on the docs table.

    Round 65 deleted four planned module rows rather than let them sit. This
    keeps the documents table honest the same way.
    """
    interface = (ROOT / "INTERFACE.md").read_text()
    section = interface[interface.index("## `docs/`"):]
    planned = [line for line in section.splitlines()
               if line.startswith("|") and line.rstrip().endswith("📋 |")]
    assert not planned, f"documents still promised but not written: {planned}"
