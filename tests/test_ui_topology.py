"""The topology diagram — Guo & MacKinnon 2017 Figure 3, for a loaded entry.

Two failure modes are specific to a topology figure and both are checked here.

The first is **silent renumbering**: if an unresolved helix were dropped rather
than drawn dashed, TM13 would appear where TM1 belongs and every label after it
would be wrong, on a picture that looks completely reasonable. So the diagram
must always carry all 38.

The second is **the wrong numbering system**. The diagram places helices by
residue number; reading a human entry with the mouse transmembrane table shifts
every one of them by up to 26 residues and again produces a plausible picture.
The window names the numbering, and the ranges it hands to the 3-D model are
in it.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.topology import (MEMBRANE_HALF, build_topology,
                                      unit_extent)


# --------------------------------------------------------------------------
# The data, without Qt
# --------------------------------------------------------------------------

def test_the_diagram_carries_all_thirty_eight_helices_in_nine_units():
    topology = build_topology("mouse")
    helices = topology.of_kind("tm_helix")
    assert len(helices) == 38
    assert topology.n_units == 9
    assert [e.helix for e in helices] == list(range(1, 39))


def test_the_units_are_fours_and_match_the_papers_own_numbering():
    """Guo & MacKinnon's '4-TM unit 6, consisting of TM 21 to 24'."""
    topology = build_topology("mouse")
    by_unit: dict[int, list[int]] = {}
    for element in topology.of_kind("tm_helix"):
        if element.unit is not None:
            by_unit.setdefault(element.unit, []).append(element.helix)
    assert len(by_unit) == 9
    for unit, helices in by_unit.items():
        assert len(helices) == 4, f"THU{unit} has {len(helices)} helices"
    assert by_unit[6] == [21, 22, 23, 24]
    # TM37 and TM38 are the pore module and belong to no unit.
    pore = [e.helix for e in topology.of_kind("tm_helix") if e.unit is None]
    assert pore == [37, 38]


def test_helices_alternate_sides_of_the_membrane():
    topology = build_topology("mouse")
    loops = topology.of_kind("loop")
    sides = [e.side for e in loops]
    assert set(sides) == {"extracellular", "cytoplasmic"}
    for a, b in zip(sides, sides[1:]):
        assert a != b, "a chain crossing the membrane must alternate sides"


def test_every_helix_spans_the_membrane_and_nothing_else_does():
    topology = build_topology("mouse")
    for element in topology.of_kind("tm_helix"):
        assert element.y0 == pytest.approx(-MEMBRANE_HALF)
        assert element.y1 == pytest.approx(MEMBRANE_HALF)
    for element in topology.of_kind("box"):
        assert element.y0 > MEMBRANE_HALF, "the cap is extracellular"
    for element in topology.of_kind("bar") + topology.of_kind("cuff"):
        assert element.y1 < -MEMBRANE_HALF, "the beam and cuff are cytoplasmic"


def test_no_two_elements_of_the_lower_rows_overlap():
    """The defect the layout was fixed for: 'bo|as' and 'PE|irp'.

    The four cuff elements fall in two tight pairs by residue number, and at
    the drawn scale their labels sat on top of one another.
    """
    topology = build_topology("mouse")
    rows: dict[float, list] = {}
    for element in topology.of_kind("cuff"):
        rows.setdefault(round(element.y0, 3), []).append(element)
    assert rows, "there should be cuff elements"
    for row in rows.values():
        row.sort(key=lambda e: e.x0)
        for a, b in zip(row, row[1:]):
            assert b.x0 >= a.x1, f"{a.label} and {b.label} overlap"


def test_the_layout_extent_contains_every_element():
    """A range stopping at the last helix clipped the hairpin off the edge."""
    topology = build_topology("mouse")
    lo, hi = topology.meta["x_range"]
    for element in topology.elements:
        assert element.x0 >= lo - 1e-9 and element.x1 <= hi + 1e-9, element.label


def test_unresolved_helices_are_marked_not_dropped(structure_6b3r):
    """Figure 3a greys out TM1-12 for this entry; ours must do it from the
    coordinates rather than from a hardcoded list."""
    topology = build_topology("mouse", structure_6b3r)
    assert len(topology.of_kind("tm_helix")) == 38, (
        "dropping an unresolved helix would renumber every one after it")
    assert set(topology.unresolved_helices) == set(range(1, 13))
    assert topology.structure == "6B3R"


def test_a_different_entry_greys_out_a_different_set(curved_structure):
    """7WLT resolves less blade than 6B3R, and the diagram follows."""
    topology = build_topology("mouse", curved_structure)
    assert set(topology.unresolved_helices) == set(range(1, 17))


def test_without_a_structure_nothing_claims_to_be_resolved():
    topology = build_topology("mouse")
    assert topology.structure is None
    assert topology.unresolved_helices == ()
    for element in topology.elements:
        assert element.resolved is None


def test_unit_extents_cover_their_helices():
    topology = build_topology("mouse")
    extents = unit_extent(topology)
    assert set(extents) == set(range(1, 10))
    for element in topology.of_kind("tm_helix"):
        if element.unit is None:
            continue
        lo, hi = extents[element.unit]
        assert lo <= element.x0 and element.x1 <= hi


def test_the_human_reference_gives_human_ranges():
    """The numbering is named, and the ranges are in it."""
    mouse = build_topology("mouse")
    human = build_topology("human")
    assert mouse.numbering == "mouse" and human.numbering == "human"
    assert mouse.sequence_length == 2547 and human.sequence_length == 2521
    mouse_beam = next(e for e in mouse.elements if e.label == "Beam")
    human_beam = next(e for e in human.elements if e.label == "Beam")
    assert (mouse_beam.start, mouse_beam.end) == (1300, 1365)
    assert (human_beam.start, human_beam.end) == (1305, 1370)


# --------------------------------------------------------------------------
# The widget
# --------------------------------------------------------------------------

@pytest.fixture
def window(qt_app, structure_6b3r):
    from piezo1.ui.topology_window import TopologyWindow

    win = TopologyWindow()
    win.resize(1200, 460)
    win.set_structure(structure_6b3r, "mouse")
    return win


def test_the_window_offers_a_box_for_every_unit(window):
    assert sorted(window._unit_boxes) == list(range(1, 10))
    assert window.view.available_units() == list(range(1, 10))


def test_ticking_a_unit_boxes_it_and_reports_its_residues(window):
    """Figure 3b's red boxes, as a selection."""
    seen = []
    window.residues_selected.connect(lambda lo, hi, n: seen.append((lo, hi, n)))
    window._unit_boxes[6].setChecked(True)
    assert window.view.boxed == {6}
    # THU6 is TM21-24, which the paper draws as its inset.
    assert seen and seen[-1] == (990, 1113, "mouse")


def test_the_selection_says_when_it_is_not_contiguous(window):
    """Boxing THU1 and THU9 highlights everything between them on the model,
    and the status line must not let that pass silently."""
    window._unit_boxes[1].setChecked(True)
    window._unit_boxes[9].setChecked(True)
    message = window.statusBar().currentMessage()
    assert "not contiguous" in message


def test_select_all_and_clear(window):
    window.select_all_units()
    assert window.view.boxed == set(range(1, 10))
    window.clear_units()
    assert window.view.boxed == set()


def test_hiding_a_kind_removes_it_from_the_drawing(window):
    assert "cuff" in window.view.shown_kinds
    window._kind_boxes["cuff"].setChecked(False)
    assert "cuff" not in window.view.shown_kinds
    assert "tm_helix" in window.view.shown_kinds


def test_the_caption_states_the_numbering_and_what_is_not_modelled(window):
    text = window.caption.text()
    assert "mouse" in text
    assert "12 helices are drawn dashed" in text or "dashed" in text
    assert "6B3R" in text


def test_the_diagram_actually_paints_something(window, tmp_path):
    """Renders to an image and counts non-background pixels.

    The same assertion the impostor tests make, for the same reason: a widget
    that uploads its geometry and draws nothing produces a clean run.
    """
    from PyQt6.QtGui import QImage, QPainter

    window.show()
    image = QImage(window.view.width(), window.view.height(),
                   QImage.Format.Format_ARGB32)
    image.fill(0)
    painter = QPainter(image)
    window.view.render(painter)
    painter.end()

    background = image.pixel(1, 1)
    lit = sum(1 for y in range(0, image.height(), 3)
              for x in range(0, image.width(), 3)
              if image.pixel(x, y) != background)
    assert lit > 200, f"only {lit} sampled pixels differ from the background"


def test_export_writes_a_png(window, tmp_path):
    window.show()
    path = window.export_png(path=str(tmp_path / "topology.png"))
    assert path and (tmp_path / "topology.png").exists()
    assert (tmp_path / "topology.png").stat().st_size > 2000


def test_an_empty_view_says_so_rather_than_crashing(qt_app):
    from PyQt6.QtGui import QImage, QPainter

    from piezo1.ui.topology_view import TopologyView

    view = TopologyView()
    view.resize(400, 200)
    assert view.topology is None
    assert view.available_units() == []
    assert view.boxed_range() == (0, 0, "")
    image = QImage(400, 200, QImage.Format.Format_ARGB32)
    image.fill(0)
    painter = QPainter(image)
    view.render(painter)          # must not raise
    painter.end()
