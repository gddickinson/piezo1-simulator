"""The family subsystem's two GUI surfaces: the colouring, and the menu entries.

A coloured trimer is the most persuasive thing this application can put on a
screen, and these values were not measured here. So the tests that matter are
not that the colouring works — they are that it **cannot reach the screen
without saying whose numbers they are**, that it **refuses an entry it cannot
read**, and that an unscored residue is **distinguishable from an unconstrained
one**, which is the difference between a coverage hole and a finding.
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("PyQt6")

from piezo1.render.colormaps import constraint_colors  # noqa: E402
from piezo1.render.representations import ColorBy  # noqa: E402


# ------------------------------------------------------- the colour map itself

def test_the_scale_is_fixed_so_two_entries_stay_comparable():
    """The whole reason this is not ColorBy.VALUE.

    An auto-ranged map repaints the same protein differently depending on how
    much blade the entry resolved. Here the colour of 0.8 must not depend on
    what else is in the array.
    """
    alone = constraint_colors(np.array([0.8]))
    among_low = constraint_colors(np.array([0.1, 0.15, 0.8]))[2]
    among_high = constraint_colors(np.array([0.8, 0.95, 0.99]))[0]
    assert np.allclose(alone[0], among_low)
    assert np.allclose(alone[0], among_high)


def test_an_unscored_residue_is_not_painted_as_an_unconstrained_one():
    """The blade tips are where coverage is worst and where low constraint is
    the claim. The two must not share a colour."""
    colours = constraint_colors(np.array([0.0, np.nan]))
    assert not np.allclose(colours[0], colours[1])
    assert np.allclose(colours[1], [0.55, 0.55, 0.55])


def test_the_map_saturates_rather_than_wrapping():
    """A value outside [0, 1] must clamp, not roll round to the other end."""
    assert np.allclose(constraint_colors(np.array([-1.0]))[0],
                       constraint_colors(np.array([0.0]))[0])
    assert np.allclose(constraint_colors(np.array([2.0]))[0],
                       constraint_colors(np.array([1.0]))[0])


# --------------------------------------------------------------- the controller

@pytest.fixture
def window(qt_app, curved_structure):
    from piezo1.ui.constraint_controller import ConstraintColourController

    class FakeView:
        def __init__(self):
            self.values = None
            self.color_by = ColorBy.DOMAIN
            self.rebuilt = 0

        def rebuild(self):
            self.rebuilt += 1

    class FakeViewport:
        def update(self):
            pass

    class FakeWindow:
        def __init__(self, structure):
            self.structure = structure
            self.view = FakeView()
            self.viewport = FakeViewport()
            self.status = []

        def _set_status(self, text):
            self.status.append(text)

        def _current_color(self):
            return ColorBy.DOMAIN

    win = FakeWindow(curved_structure)
    win.constraint_colour = ConstraintColourController(win)
    return win


def test_the_colouring_cannot_reach_the_screen_without_saying_whose_numbers(window):
    window.constraint_colour.show(True)
    assert window.constraint_colour.visible
    assert window.view.color_by is ColorBy.CONSTRAINT
    status = window.status[-1]
    assert "NOT MEASURED HERE" in status
    assert "piezo_genes" in status
    assert "FIXED" in status
    assert "UNSCORED" in status


def test_it_paints_the_track_and_holds_unscored_atoms_out(window):
    window.constraint_colour.show(True)
    values = np.asarray(window.view.values)
    assert values.size == window.structure.n_atoms
    assert np.isnan(values).any(), "a real entry has unresolved residues"
    assert np.nanmax(values) <= 1.0 and np.nanmin(values) >= 0.0


def test_turning_it_off_restores_the_previous_colouring(window):
    window.constraint_colour.show(True)
    window.constraint_colour.show(False)
    assert not window.constraint_colour.visible
    assert window.view.color_by is ColorBy.DOMAIN
    assert window.view.values is None


def test_an_entry_in_another_numbering_is_refused_out_loud(qt_app, structure_by_id):
    """A PIEZO2 entry coloured by PIEZO1's constraint would be wrong at every
    residue and would look exactly like a result."""
    from piezo1.ui.constraint_controller import ConstraintColourController

    piezo2 = structure_by_id("6KG7")
    if piezo2 is None:
        pytest.skip("6KG7 is not downloaded")

    class FakeWindow:
        def __init__(self):
            self.structure = piezo2
            self.view = type("V", (), {"values": None,
                                       "color_by": ColorBy.DOMAIN,
                                       "rebuild": lambda self: None})()
            self.viewport = type("P", (), {"update": lambda self: None})()
            self.status = []

        def _set_status(self, text):
            self.status.append(text)

        def _current_color(self):
            return ColorBy.DOMAIN

    win = FakeWindow()
    controller = ConstraintColourController(win)
    controller.show(True)
    assert not controller.visible
    assert win.view.color_by is ColorBy.DOMAIN
    assert "REFUSED" in win.status[-1]


def test_reset_forgets_the_old_entry_rather_than_repainting(window):
    """For the structure-replacement path: a stale result is a status line
    quoting one structure's coverage over another's picture."""
    window.constraint_colour.show(True)
    before = window.view.rebuilt
    window.constraint_colour.reset()
    assert not window.constraint_colour.visible
    assert window.constraint_colour.result is None
    assert window.view.rebuilt == before, "reset must not repaint"


# ------------------------------------------------------------- the menu entries

def test_every_family_analysis_has_a_gui_entry_and_a_caveat():
    from piezo1.analysis.report import ANALYSES
    from piezo1.ui.tabular_analyses import CAVEATS, TabularAnalysisMixin

    for key, method in (("family", "show_family"),
                        ("constraint", "show_constraint"),
                        ("disease", "show_disease_geography"),
                        ("coreperiphery", "show_core_periphery"),
                        ("piezo3", "show_piezo3")):
        assert key in ANALYSES
        assert hasattr(TabularAnalysisMixin, method), method
        assert key in CAVEATS, f"{key} would be shown with no caveat"
        assert len(CAVEATS[key]) > 200, f"{key}'s caveat is too thin to help"


def test_the_caveats_lead_with_what_the_reader_would_get_wrong():
    """Each of these results has one specific misreading, and the caveat has to
    block it in its first sentence rather than in a footnote."""
    from piezo1.ui.tabular_analyses import CAVEATS

    assert CAVEATS["family"].startswith("AN EXTERNAL PROJECT'S RESULTS")
    assert "THE PER-RESIDUE VALUES ARE THE CENSUS'S" in CAVEATS["constraint"]
    assert "BOUNDARY-DEPENDENT" in CAVEATS["constraint"]
    assert "not a reproduction" in CAVEATS["disease"]
    assert "SAME PROTEIN" in CAVEATS["coreperiphery"]
    assert "no current has ever been recorded" in CAVEATS["piezo3"]
