"""The drawn pockets, and the arithmetic the picture invites you to get wrong.

A pocket is a *cluster of overlapping spheres*. Drawn, it looks like a pile of
beads, and the obvious reading of a pile of beads is that its volume is the sum
of the beads. That reading is wrong by several fold, which is why
`Pocket.volume` integrates a Monte-Carlo union — and the size of the error is
measured here rather than asserted, so the caveat on the status line is backed
by a number from this structure.

The rest is the same join the drawn contacts got wrong: the picture must be the
detector's own output, in the detector's own ranking, or the table and the
picture are two different answers with nothing to say which is right.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.ui.pocket_controller import (DEFAULT_TOP,  # noqa: E402
                                         POCKET_COLORS, PocketController)


@pytest.fixture(scope="module")
def pockets():
    path = STRUCTURE_DIR / "8YEZ.cif"
    if not path.exists():
        pytest.skip("8YEZ.cif not downloaded — run python -m piezo1.io.fetch")
    from piezo1.analysis.pockets import find_pockets

    return find_pockets(Structure.from_file(path))


class _FakeWindow:
    def __init__(self, pockets):
        self.structure = object()
        self.analysis = type("A", (), {"pockets": pockets})()
        self.status = ""

    def _set_status(self, text):
        self.status = text


@pytest.fixture(scope="module")
def controller(pockets):
    return PocketController(_FakeWindow(pockets))


def test_the_controller_draws_the_detector_ranking_not_its_own(controller,
                                                               pockets):
    """Identity and order. A second `find_pockets` call would also be a second
    ranking, and two rankings that disagree is a defect with no symptom."""
    assert controller.pockets is pockets
    drawn = controller.drawn()
    assert drawn == pockets[:controller.top]
    assert [p.index for p in drawn] == list(range(1, len(drawn) + 1)), \
        "the drawn pockets are not the top-ranked ones in rank order"


def test_only_the_top_few_are_drawn_and_the_rest_are_counted(controller):
    """Not dropped silently — the contact drawing's own rule.

    The detector returns up to 30. Drawing all of them fills the protein and
    the ranking stops meaning anything, but a picture of five that reads as a
    picture of all is worse than either.
    """
    assert controller.top == DEFAULT_TOP < 30
    if len(controller.pockets) > controller.top:
        assert "not drawn" in controller.status_line()


def test_every_drawn_pocket_gets_a_colour_of_its_own(controller):
    """Two adjacent pockets in one colour are one pocket on screen."""
    assert len(POCKET_COLORS) >= DEFAULT_TOP
    used = [POCKET_COLORS[i % len(POCKET_COLORS)]
            for i in range(len(controller.drawn()))]
    assert len(set(used)) == len(used)


def test_summing_the_drawn_spheres_overcounts_the_volume(controller):
    """The measurement behind the caveat, on this structure's own pockets.

    If the spheres ever stopped overlapping, counting them *would* be a
    volume and the warning would be noise. They do overlap, by a lot, and this
    says by how much.
    """
    worst = 0.0
    for pocket in controller.drawn():
        naive = float((4.0 / 3.0 * np.pi * np.asarray(pocket.radii) ** 3).sum())
        union = pocket.volume
        assert union > 0.0, "a drawn pocket reports no volume at all"
        assert naive > union, (
            "summing the alpha spheres did not overcount, so they are not "
            "overlapping and the status line's warning is wrong")
        worst = max(worst, naive / union)
    assert worst > 2.0, (
        f"the worst overcount is only {worst:.1f}x; the caveat claims several "
        f"fold and should be reworded if this is now the truth")


def test_the_status_line_keeps_a_cavity_from_reading_as_a_binding_site(controller):
    """Three claims the picture would otherwise make on its own.

    That the beads are a volume; that a cavity is a site; and that an empty
    pocket is empty — when detection removed the lipid sitting in it first.
    """
    line = controller.status_line().lower()
    assert "overlap" in line
    assert "not a binding site" in line
    assert "lipid" in line


def test_nothing_detected_says_so_rather_than_drawing_an_empty_picture():
    controller = PocketController(_FakeWindow([]))
    assert controller.drawn() == []
    assert "no pockets" in controller.status_line()
