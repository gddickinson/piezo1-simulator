"""HUD, presentation mode and the structure overlay.

Offscreen widgets and Qt-free logic where possible. The scale bar and the
overlay both carry a claim about the science, so those are what is asserted:
a bar of stated length, and a superposition that does not trust chain labels.
"""

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PyQt6.QtWidgets")
from PyQt6.QtWidgets import QApplication  # noqa: E402

from piezo1.ui.hud import HudSettings, nice_scale_length  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance()
    if app is None:
        try:
            app = QApplication([])
        except Exception as exc:                    # pragma: no cover
            pytest.skip(f"no Qt platform: {exc}")
    return app


# --------------------------------------------------------------------------
# Scale bar
# --------------------------------------------------------------------------

def test_scale_bar_lengths_are_round_numbers():
    """A bar labelled 47 Å is unreadable; the point is a number you can hold."""
    for target in (3.4, 17.0, 46.0, 230.0, 1400.0):
        length = nice_scale_length(target)
        assert length <= target
        assert length in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)


def test_scale_bar_never_exceeds_the_space_it_was_given():
    for target in np.linspace(1.0, 900.0, 40):
        assert nice_scale_length(float(target)) <= target


def test_tiny_target_still_returns_a_usable_length():
    assert nice_scale_length(0.2) == 1.0


def test_hud_settings_round_trip():
    settings = HudSettings(scale_bar=False, clock=True, font_scale=1.6,
                           fields=["pore", "dome"])
    restored = HudSettings.from_dict(settings.as_dict())
    assert restored.clock and not restored.scale_bar
    assert restored.font_scale == pytest.approx(1.6)
    assert restored.fields == ["pore", "dome"]


def test_unknown_settings_keys_are_ignored():
    """A settings file from an older build must not stop the app starting."""
    restored = HudSettings.from_dict({"scale_bar": False, "invented_key": 7})
    assert restored.scale_bar is False


def test_hud_readouts_can_be_set_and_cleared(qapp):
    from piezo1.ui.hud import HudOverlay
    from PyQt6.QtWidgets import QWidget

    parent = QWidget()
    hud = HudOverlay(parent)
    hud.set_readout("pore", "bottleneck 0.95 Å")
    assert hud.readouts["pore"].startswith("bottleneck")
    hud.set_readout("pore", "")
    assert "pore" not in hud.readouts
    hud.set_clock("1.20 s", "note")
    assert hud.clock_text == "1.20 s"
    hud.clear_readouts()
    assert not hud.readouts


def test_world_per_pixel_is_zero_without_a_scene(qapp):
    from piezo1.ui.hud import HudOverlay
    from PyQt6.QtWidgets import QWidget
    hud = HudOverlay(QWidget())
    assert hud.world_per_pixel() == 0.0


# --------------------------------------------------------------------------
# Readout catalogue
# --------------------------------------------------------------------------

def test_every_readout_has_a_label():
    from piezo1.ui.presentation import READOUTS
    keys = [k for k, _ in READOUTS]
    assert len(set(keys)) == len(keys), "duplicate readout keys"
    for key, label in READOUTS:
        assert key and label and label[0].isupper()


# --------------------------------------------------------------------------
# Structure overlay
# --------------------------------------------------------------------------

def test_overlay_rematches_protomers_rather_than_trusting_labels(
        curved_structure, flat_structure):
    """The trap: deposited entries label protomers in either rotational order.

    Taken at chain-label face value 7WLT and 7WLU sit ~90 Å apart; matched
    first they are ~12 Å. A viewer that trusted labels would show an enormous
    conformational change that is not there.
    """
    from piezo1.ui.overlay_controller import OverlayWorker

    result = OverlayWorker(curved_structure, flat_structure,
                           "protomer")._superpose()
    assert result.reordered, "7WLT/7WLU are known to disagree on protomer order"
    assert result.rmsd < 20.0
    assert result.rmsd_by_label > 3 * result.rmsd
    assert result.n_common > 1000


def test_overlay_reports_per_residue_deviation(curved_structure,
                                               flat_structure):
    from piezo1.ui.overlay_controller import OverlayWorker
    result = OverlayWorker(curved_structure, flat_structure,
                           "protomer")._superpose()
    assert len(result.per_residue) == result.n_common
    assert all(v >= 0 for v in result.per_residue.values())
    assert result.meta["max_deviation"] > result.rmsd


def test_superposition_is_a_rigid_motion(curved_structure, flat_structure):
    """Kabsch must not scale or reflect — either would fake a better fit."""
    from piezo1.ui.overlay_controller import OverlayWorker
    result = OverlayWorker(curved_structure, flat_structure,
                           "protomer")._superpose()
    rotation = result.rotation
    assert np.allclose(rotation @ rotation.T, np.eye(3), atol=1e-8)
    assert np.linalg.det(rotation) == pytest.approx(1.0, abs=1e-8)


def test_overlay_needs_shared_residues(human_structure):
    """Too few residues in common must raise, not fit noise."""
    import dataclasses
    from piezo1.ui.overlay_controller import OverlayWorker

    shifted = dataclasses.replace(
        human_structure, res_seq=human_structure.res_seq + 100000)
    with pytest.raises(ValueError, match="in common"):
        OverlayWorker(human_structure, shifted, "chain")._superpose()
