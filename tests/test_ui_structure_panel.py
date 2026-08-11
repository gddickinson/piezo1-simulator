"""The structure chooser's filters, on the real widgets.

Twenty-one entries in one combo is a list to scroll rather than a choice to
make, and the one that is a different *molecule* — 6KG7, PIEZO2 — was reachable
only by knowing which of the fifteen mouse entries it was. These run offscreen
with real Qt widgets, because the two ways this feature goes wrong are both
behavioural: options that go stale against the data, and a filter that hides an
entry something else then asks for by name.
"""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.io.registry import load_registry  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance()
    if app is None:
        try:
            app = QApplication([])
        except Exception as exc:                          # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return app


@pytest.fixture
def panel(qapp):
    from piezo1.ui.panels.structure_panel import StructurePanel

    made = StructurePanel()
    if not made._records:
        pytest.skip("no structures downloaded — run python -m piezo1.io.fetch")
    return made


def _options(panel, field):
    combo = panel.filter_combos[field]
    return [combo.itemText(i) for i in range(combo.count())]


# ------------------------------------------------- the options come from data

def test_every_filters_options_are_read_from_the_records(panel):
    """Hard-coded options are how a box goes stale.

    The species box listed ``["all", "human", "mouse"]`` in the source. That
    happened to be right, and it is right for exactly as long as nobody adds an
    entry of a third species — so the options are now the values the field
    actually takes.
    """
    from piezo1.ui.panels.structure_panel import FILTERS

    available = load_registry().available()
    for field, _title, tip in FILTERS:
        options = _options(panel, field)
        assert options[0] == "all", field
        assert options[1:] == sorted({str(getattr(r, field)) for r in available}), field
        assert panel.filter_combos[field].toolTip() == tip


def test_the_protein_filter_separates_piezo2_from_piezo1(panel):
    """The one the request was about."""
    assert "PIEZO2" in _options(panel, "protein")

    panel.filter_combos["protein"].setCurrentText("PIEZO2")
    assert [r.pdb for r in panel._records] == ["6KG7"]
    assert panel.current_record().pdb == "6KG7"

    panel.filter_combos["protein"].setCurrentText("PIEZO1")
    assert "6KG7" not in [r.pdb for r in panel._records]
    assert len(panel._records) > 10


def test_the_filters_combine(panel):
    panel.filter_combos["protein"].setCurrentText("PIEZO1")
    panel.filter_combos["gating"].setCurrentText("open-like")
    assert [r.pdb for r in panel._records] == ["11ZC"], (
        "11ZC is the only open-like entry in the catalogue")


def test_the_count_says_how_much_is_hidden(panel):
    total = len(load_registry().available())
    assert f"of {total}" in panel.count.text()
    assert panel.count.text().startswith(f"{total} of")

    panel.filter_combos["protein"].setCurrentText("PIEZO2")
    assert panel.count.text().startswith("1 of")


def test_an_empty_combination_says_so_rather_than_leaving_stale_detail(panel):
    """The previous entry's details under an empty chooser read as current."""
    panel.filter_combos["protein"].setCurrentText("PIEZO2")
    before = panel.detail.text()
    assert "6KG7" in before or "PIEZO2" in before or before

    panel.filter_combos["species"].setCurrentText("human")
    assert panel._records == []
    assert panel.current_record() is None
    assert "No entry matches" in panel.detail.text()
    assert panel.count.text().startswith("0 of")


# ---------------------------------------------- the trap: select past a filter

def test_select_clears_a_filter_that_would_hide_the_entry(panel):
    """A silent no-op here would leave the wrong structure on screen.

    ``--structure 6KG7``, a restored session and the morph controller all
    arrive through ``select``. With a filter set, the old code found no match
    and returned, having appeared to honour the request.
    """
    panel.filter_combos["protein"].setCurrentText("PIEZO1")
    assert "6KG7" not in [r.pdb for r in panel._records]

    panel.select("6KG7")
    assert panel.current_record() is not None
    assert panel.current_record().pdb == "6KG7"
    assert panel.filter_combos["protein"].currentText() == "all"


def test_select_leaves_the_filters_alone_when_it_does_not_need_them(panel):
    """Clearing filters on every select would be its own annoyance."""
    panel.filter_combos["protein"].setCurrentText("PIEZO1")
    target = panel._records[-1].pdb
    panel.select(target)
    assert panel.current_record().pdb == target
    assert panel.filter_combos["protein"].currentText() == "PIEZO1"


def test_select_of_something_not_downloaded_changes_nothing(panel):
    panel.filter_combos["protein"].setCurrentText("PIEZO1")
    before = panel.current_record().pdb
    panel.select("1ABC")
    assert panel.current_record().pdb == before
    assert panel.filter_combos["protein"].currentText() == "PIEZO1"


def test_the_species_combo_alias_still_points_at_the_species_filter(panel):
    """`species_combo` was the only filter and other code reads it by name."""
    assert panel.species_combo is panel.filter_combos["species"]


# ------------------------------------------------------- the registry field

def test_the_protein_field_is_measured_and_present_on_every_entry():
    """A stale resource must be visibly stale, not silently PIEZO1.

    The default is "unknown" for exactly that reason: the species field said
    "mouse" for PIEZO2 and that is how it slipped past the overlay guard.
    """
    from piezo1.core.numbering_check import identify_numbering
    from piezo1.core import Structure

    records = load_registry().available()
    assert records, "no structures downloaded"
    for record in records:
        assert record.protein in ("PIEZO1", "PIEZO2"), record.pdb
        measured = identify_numbering(Structure.from_file(record.path))
        assert record.is_piezo2 == measured.is_piezo2, (
            f"{record.pdb}: registry says {record.protein}, the coordinates "
            f"say {measured.reference}")
    assert sum(r.is_piezo2 for r in records) == 1
