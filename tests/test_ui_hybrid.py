"""The full-length model on screen, and the confidence colouring it depends on.

Round 76's validation clause is not "it draws": it is that **the seam must be
visibly rendered**, because a complete-looking PIEZO1 trimer whose 569 distal
residues are AlphaFold is precisely the confident-wrong-picture failure the
Round 50 hazard audit exists for.

Building it found that the signal it relies on did not work. `plddt_colors`
applied AlphaFold's bands highest-threshold-first, so the final pass at
``>= 0.0`` overwrote every atom: **Colour by → AlphaFold pLDDT painted the whole
model one flat orange**, and had since the feature was written. A confidence
colouring showing no variation is worse than none — it reads as uniformly bad.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.render.colormaps import (PLDDT_BANDS, plddt_band_colors,  # noqa: E402
                                     plddt_colors)

PREDICTED_MODEL = "AF-Q92508-F1-model_v6.cif"


# ------------------------------------------------- the banding, on its own

def test_each_confidence_band_gets_its_own_colour():
    """The bug, pinned on values chosen to land one in each band."""
    values = np.array([95.0, 80.0, 60.0, 30.0])
    colors = plddt_band_colors(values)
    assert len({tuple(c) for c in colors}) == 4, (
        "the pLDDT bands collapse to one colour; they are being applied in "
        "descending order, so >= 0.0 overwrites everything")


def test_a_value_takes_the_colour_of_the_highest_band_it_clears():
    bands = dict(sorted(PLDDT_BANDS))
    for value, threshold in ((95.0, 90.0), (90.0, 90.0), (89.9, 70.0),
                             (70.0, 70.0), (50.0, 50.0), (0.0, 0.0)):
        got = tuple(plddt_band_colors(np.array([value]))[0])
        assert got == pytest.approx(bands[threshold]), (
            f"pLDDT {value} should take the {threshold} band")


def test_the_real_predicted_model_shows_every_band():
    path = STRUCTURE_DIR / PREDICTED_MODEL
    if not path.exists():
        pytest.skip("AlphaFold model not downloaded; run python -m piezo1.io.fetch")
    predicted = Structure.from_file(path)
    colors = plddt_colors(predicted)
    assert len({tuple(c) for c in colors}) == 4, (
        "PIEZO1's prediction spans 22-95 pLDDT; a single colour means the "
        "banding is broken again")


# ----------------------------------------------------- the drawn model

@pytest.fixture(scope="module")
def controller():
    for name in ("8YEZ.cif", PREDICTED_MODEL):
        if not (STRUCTURE_DIR / name).exists():
            pytest.skip(f"{name} not downloaded; run python -m piezo1.io.fetch")

    class _Batch:
        def upload(self, *args, **kwargs):
            self.args = args

    class _Scene:
        def __init__(self):
            self.batches = {}

        def spheres(self, name):
            return self.batches.setdefault(name, _Batch())

        cylinders = spheres

        def remove(self, name):
            self.batches.pop(name, None)

    class _Viewport:
        def __init__(self):
            self.scene = _Scene()

        def update(self):
            pass

    class _Window:
        def __init__(self, structure):
            self.structure = structure
            self.viewport = _Viewport()
            self.status = ""

        def _set_status(self, text):
            self.status = text

    from piezo1.ui.hybrid_controller import HybridController

    window = _Window(Structure.from_file(STRUCTURE_DIR / "8YEZ.cif"))
    return HybridController(window)


@pytest.fixture
def drawn(controller):
    controller.clear()
    controller.show(True)
    if controller.model is None:
        pytest.skip(f"could not build the model: {controller.win.status}")
    return controller


def test_the_experimental_and_predicted_parts_are_not_the_same_colour(drawn):
    """The whole point. One colour would be the failure this guards against."""
    from piezo1.ui.hybrid_controller import EXPERIMENTAL_COLOR

    colors = drawn._colors()
    model = drawn.model
    grey = (colors == np.float32(EXPERIMENTAL_COLOR)).all(axis=1)

    assert grey.sum() == int(model.experimental_only.sum())
    assert not grey[model.predicted].any(), \
        "a predicted atom is drawn in the experimental colour"
    assert len({tuple(c) for c in colors[model.predicted]}) > 1, \
        "the graft is one flat colour; the confidence banding is not reaching it"


def test_the_seam_is_marked(drawn):
    """It is where the model is least trustworthy and otherwise invisible."""
    assert "hybrid:seam" in drawn.win.viewport.scene.batches
    point = drawn._seam_point()
    assert point is not None and point.shape == (3,)
    at_seam = drawn.model.res_seq == drawn.model.seam_residue
    assert at_seam.any(), "the marker is not on a residue the model contains"


def test_the_graft_covers_what_the_experiment_does_not_resolve(drawn):
    model = drawn.model
    residues = model.res_seq[model.predicted]
    assert int(residues.min()) == 1
    assert int(residues.max()) == model.seam_residue - 1
    assert np.unique(residues).size == 569


def test_it_cannot_be_drawn_without_saying_the_blade_is_predicted(drawn):
    """The guard. A caption is the only thing between this picture and a
    reader who takes the distal blade for experimental structure."""
    status = drawn.win.status
    assert "PREDICT" in status.upper()
    assert "1" in status and "569" in status, "the range must be stated"
    assert "NOT validated" in status


def test_the_status_carries_the_number_a_good_local_fit_hides(drawn):
    """2.4 Å at the seam looks like a validated join. It is not: the two
    models differ by 75 Å over the region they share."""
    assert "75" in drawn.win.status
    assert drawn.model.global_rmsd > 50.0
    assert drawn.model.overlap_rmsd < 5.0


def test_less_than_half_the_graft_is_confident_and_that_is_reported(drawn):
    model = drawn.model
    fraction = float(model.confident_prediction.sum()) / int(model.predicted.sum())
    assert 0.4 < fraction < 0.6, f"confident fraction moved to {fraction:.2f}"
    assert f"{fraction:.0%}" in drawn.win.status


def test_turning_it_off_removes_every_batch(drawn):
    drawn.show(False)
    assert not [k for k in drawn.win.viewport.scene.batches
                if k.startswith("hybrid:")]
    assert drawn.model is None


def test_the_analysis_and_the_drawing_agree(drawn):
    """One model, two surfaces; they must not diverge."""
    from piezo1.analysis.report import ANALYSES

    result = ANALYSES["hybrid"](drawn.win.structure, "human")
    assert result["seam_residue"] == drawn.model.seam_residue
    assert result["predicted_atoms"] == int(drawn.model.predicted.sum())
    assert result["global_rmsd_A"] == pytest.approx(drawn.model.global_rmsd)


def test_the_analysis_is_reachable_from_the_registry_with_a_caveat():
    from piezo1.analysis.report import ANALYSES
    from piezo1.ui.tabular_analyses import CAVEATS

    assert "hybrid" in ANALYSES
    assert "hybrid" in CAVEATS
    assert "PREDICTION" in CAVEATS["hybrid"]
