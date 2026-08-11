"""Colouring by predicted fluctuation, beside the mode-displacement colouring.

Two buttons, one slot. Both drive `ColorBy.VALUE` through `view.values`, so the
failure mode is not a crash — it is a lit button describing a colour that is
not on screen, which reads as a finding about the protein.

The other thing worth pinning is that the new button is not the old one under a
different name. A single mode's displacement and the mode set's mean-square
fluctuation are different quantities and are routinely confused; if they
happened to be proportional on real coordinates, the button would be a second
way to see one thing and should be removed rather than kept.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")

from piezo1.analysis.fluctuations import predicted_msf  # noqa: E402
from piezo1.physics.anm import ANM  # noqa: E402
from piezo1.render.representations import ColorBy  # noqa: E402
from piezo1.structure.frame import apply_frame, canonical_transform  # noqa: E402
from piezo1.structure.protomers import protomer_blocks  # noqa: E402
from piezo1.ui.physics_controller import PhysicsController  # noqa: E402


@pytest.fixture(scope="module")
def modes(human_structure):
    st = apply_frame(human_structure, canonical_transform(human_structure))
    blocks, _residues = protomer_blocks(st)
    anm = ANM.from_trimer(blocks, cutoff=15.0).build()
    return anm.calc_modes(n_modes=20)


class _View:
    def __init__(self):
        self.values = None
        self.color_by = ColorBy.CHAIN
        self.rebuilt = 0

    def rebuild(self):
        self.rebuilt += 1


class _Button:
    def __init__(self, checked=False):
        self._checked = checked
        self.blocked = False

    def isChecked(self):
        return self._checked

    def setChecked(self, value):
        self._checked = bool(value)

    def blockSignals(self, value):
        self.blocked = bool(value)


class _FakeWindow:
    def __init__(self, structure, modes):
        self.structure = structure
        self.modes = modes
        self.view = _View()
        self.viewport = type("V", (), {"update": lambda self: None})()
        self.physics_panel = type("P", (), {})()
        self.physics_panel.color_button = _Button()
        self.physics_panel.fluctuation_button = _Button()
        self._mode_blocks = []
        self.status = ""

    def _set_status(self, text):
        self.status = text

    def _current_color(self):
        return ColorBy.CHAIN


@pytest.fixture(scope="module")
def controller(human_structure, modes):
    st = apply_frame(human_structure, canonical_transform(human_structure))
    blocks, _residues = protomer_blocks(st)
    win = _FakeWindow(st, modes)
    win._mode_blocks = blocks
    return PhysicsController(win)


# --------------------------------------------------------------------------
# It is the validated quantity, not a lookalike
# --------------------------------------------------------------------------

def test_the_colours_are_the_array_round_82_validated(controller, modes):
    """`predicted_msf` correlates `modes.msf()` against the B-factors. If the
    picture used anything else, it would be showing a quantity nothing tested.
    """
    controller.color_by_fluctuation(True)
    values = controller.win.view.values
    assert values is not None
    assert controller.win.view.color_by is ColorBy.VALUE

    # Every atom carries the value of the site it was mapped to.
    msf = np.asarray(modes.msf(), dtype=float)
    assert np.allclose(values, msf[controller._ca_map])

    # And that array is the one the validation averages over protomers.
    folded = predicted_msf(modes)
    per = len(msf) // 3
    assert np.allclose(folded, msf[:per * 3].reshape(3, per).mean(axis=0))


def test_the_three_protomers_are_not_folded_together(controller, modes):
    """`predicted_msf` averages them because the *observation* is per chain.
    Here the three copies are on screen separately, so averaging would paint an
    agreement the model does not have."""
    msf = np.asarray(modes.msf(), dtype=float)
    per = len(msf) // 3
    spread = np.abs(msf[:per] - msf[per:2 * per]).max()
    assert spread >= 0.0
    controller.color_by_fluctuation(True)
    assert len(np.unique(controller.win.view.values)) > 1


def test_fluctuation_is_not_the_same_picture_as_one_mode(controller, modes):
    """If it were, one of the two buttons should go."""
    controller.select_mode(0)
    controller.color_by_mode(True)
    displacement = np.array(controller.win.view.values, dtype=float, copy=True)
    controller.color_by_fluctuation(True)
    fluctuation = np.array(controller.win.view.values, dtype=float, copy=True)

    a = displacement - displacement.mean()
    b = fluctuation - fluctuation.mean()
    r = float((a * b).sum() / np.sqrt((a * a).sum() * (b * b).sum()))
    assert abs(r) < 0.95, (
        f"the two colourings correlate at r = {r:.3f}; they are one picture "
        f"with two buttons")


# --------------------------------------------------------------------------
# One slot, two buttons
# --------------------------------------------------------------------------

def test_turning_one_on_turns_the_other_off(controller):
    panel = controller.win.physics_panel
    panel.color_button.setChecked(True)
    controller.color_by_fluctuation(True)
    assert not panel.color_button.isChecked(), \
        "the displacement button stayed lit while fluctuation was painted"

    panel.fluctuation_button.setChecked(True)
    controller.color_by_mode(True)
    assert not panel.fluctuation_button.isChecked()


def test_unticking_the_other_button_does_not_re_enter_its_handler(controller):
    """Signals are blocked while it is cleared. Without that, unchecking one
    calls its handler with False, which resets the colouring the other just
    set — the button works and the model stays grey."""
    panel = controller.win.physics_panel
    panel.color_button.setChecked(True)
    controller.color_by_fluctuation(True)
    assert panel.color_button.blocked is False, "signals were left blocked"
    assert controller.win.view.color_by is ColorBy.VALUE, \
        "the colouring was undone while the other button was cleared"


def test_turning_it_off_restores_the_ordinary_colouring(controller):
    controller.color_by_fluctuation(True)
    controller.color_by_fluctuation(False)
    assert controller.win.view.values is None
    assert controller.win.view.color_by is ColorBy.CHAIN


# --------------------------------------------------------------------------
# What it says about itself
# --------------------------------------------------------------------------

def test_the_status_line_says_the_scale_is_arbitrary(controller):
    """A colour ramp implies a quantity. This one has no fitted spring
    constant, so only the ordering means anything — and whether the ordering
    is right is measured somewhere else, which the line has to point at."""
    controller.color_by_fluctuation(True)
    line = controller.win.status
    assert "ARBITRARY" in line
    assert "Fluctuation vs B-factor" in line
    assert "control" in line


def test_the_quoted_numbers_match_the_measured_survey(controller):
    """The line quotes Round 82's medians. They are the survey's, so if the
    survey moves the sentence is wrong — pinned here rather than in prose."""
    line = controller.fluctuation_line()
    assert "0.74" in line and "0.32" in line
    assert "0.48" in line and "0.39" in line, \
        "the Pearson half of the result was dropped, which is the half where " \
        "the network barely beats the control"


def test_it_does_nothing_without_modes(human_structure):
    win = _FakeWindow(human_structure, None)
    controller = PhysicsController(win)
    controller.color_by_fluctuation(True)
    assert win.view.values is None
    assert win.view.color_by is ColorBy.CHAIN
