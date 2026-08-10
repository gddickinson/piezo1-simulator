"""Drawing the HaloTag: the sphere, the fold, and the caption that must survive.

The controller had no tests at all before this. That matters more here than for
most panels, because the fold is the one thing this project draws that *looks
like* an experimental result and is not one: the tag's coordinates are
deposited, its position is modelled, and its orientation is nothing at all.

So the load-bearing test is not that the atoms appear. It is that they cannot
appear **without** the statement that the orientation is undetermined — the
same shape of guard the project uses on its documents, applied to a status bar.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.structure.fusion import build_fusion, load_halotag  # noqa: E402
from piezo1.ui.fusion_controller import (CONTACT_COLOR, DYE_COLOR, TAG_COLOR,
                                         FusionController)  # noqa: E402


class _Batch:
    def __init__(self):
        self.args = ()

    def upload(self, *args, **kwargs):
        self.args = args


class _Scene:
    """Enough of the scene to record what a controller drew, without a GL context."""

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
    """The three attributes the controller touches, and nothing else."""

    def __init__(self, structure):
        self.structure = structure
        self.viewport = _Viewport()
        self.status = ""

    def _set_status(self, text):
        self.status = text


@pytest.fixture(scope="module")
def controller():
    for pdb in ("8YEZ", "6U32"):
        if not (STRUCTURE_DIR / f"{pdb}.cif").exists():
            pytest.skip(f"{pdb} not downloaded; run python -m piezo1.io.fetch")
    window = _Window(Structure.from_file(STRUCTURE_DIR / "8YEZ.cif"))
    return FusionController(window)


@pytest.fixture
def drawn(controller):
    controller.clear()
    controller.show_atoms = False
    controller.show_dyes = False
    controller.show_envelope = False
    controller.spin = None
    controller.show(True)
    return controller


# ------------------------------------------------------------ the two bodies

def test_the_sphere_is_what_is_drawn_by_default(drawn):
    """The shape that claims exactly what the model determined."""
    assert "halotag:tags" in drawn.win.viewport.scene.batches
    assert "halotag:fold" not in drawn.win.viewport.scene.batches
    assert drawn.pose is None
    radii = drawn.win.viewport.scene.batches["halotag:tags"].args[1]
    assert radii == pytest.approx(drawn.model.meta["tag_radius"])


def test_the_fold_replaces_the_sphere_rather_than_joining_it(drawn):
    """Both at once would read as a tag inside a bubble."""
    drawn.set_atoms(True)
    batches = drawn.win.viewport.scene.batches
    assert "halotag:fold" in batches
    assert "halotag:tags" not in batches

    coords, radii, colours, _ = batches["halotag:fold"].args
    assert len(coords) == drawn.pose.n_atoms * drawn.pose.n_tags
    assert len(radii) == len(coords) and len(colours) == len(coords)
    assert radii.max() < 1.0, "drawn well under van der Waals, or it is a blob"


def test_the_atoms_inside_the_channel_are_coloured_as_such(drawn):
    """The reported contact count has to be visible, not only stated."""
    from piezo1.structure.fusion_pose import SPIN_SAMPLES, spin_scan

    counts = spin_scan(drawn.win.structure, drawn.model)
    drawn.set_atoms(True)
    drawn.spin = float(np.argmax(counts) * 2 * np.pi / SPIN_SAMPLES)
    drawn._draw()

    colours = drawn.win.viewport.scene.batches["halotag:fold"].args[2]
    red = (colours == np.float32(CONTACT_COLOR)).all(axis=1)
    assert red.sum() > 0, "the worst orientation is well inside the channel"
    assert red.sum() == drawn.pose.body_contacts * drawn.pose.n_tags
    assert (colours == np.float32(TAG_COLOR)).all(axis=1).sum() > 1000


def test_the_seam_ends_on_a_real_atom_once_the_fold_is_shown(drawn):
    """With a sphere the linker runs to a notional middle; with the fold it does not."""
    drawn.set_atoms(True)
    starts, ends = drawn.win.viewport.scene.batches["halotag:seam"].args[:2]
    assert ends == pytest.approx(np.asarray(drawn.pose.anchors, np.float32),
                                 abs=1e-3)
    assert starts == pytest.approx(np.asarray(drawn.model.anchors, np.float32),
                                   abs=1e-3)


# --------------------------------------------------- the undetermined angle

def test_the_fold_can_never_be_drawn_without_saying_the_spin_is_unknown(drawn):
    """The guard this file exists for.

    A drawn fold reads as a determined pose. Every path that puts one on
    screen — first draw, a turn, and a turn with dyes also on — must leave the
    caveat in the status line.
    """
    drawn.set_atoms(True)
    assert "UNDETERMINED" in drawn.win.status

    drawn.rotate_tags()
    assert "UNDETERMINED" in drawn.win.status

    drawn.set_dyes(True)
    assert "UNDETERMINED" in drawn.win.status, (
        "the dye count replaced the caveat instead of prefixing it")
    assert "labelled" in drawn.win.status


def test_colour_is_not_what_distinguishes_a_modelled_tag(drawn):
    """Measured, because the comment used to claim otherwise.

    The controller said its colours were "deliberately unlike any colouring the
    channel uses". They are not: the tag's orange sits 0.10 from the chain
    palette's orange and the dye's red 0.10 from its red, and there is nowhere
    to move to — every colour far from those eight hues is too dark to see on
    the background. That is tolerable while the guard is the status line, so
    this pins both halves together: if someone finds a free colour and drops
    the caveat believing colour now carries it, the second half fails.
    """
    from piezo1.render.colormaps import SEQUENCE_COLORS

    palette = np.asarray(SEQUENCE_COLORS)
    nearest = np.linalg.norm(palette - np.asarray(TAG_COLOR), axis=1).min()
    assert nearest < 0.2, (
        f"TAG_COLOR is now {nearest:.2f} from every chain colour; if it is "
        f"genuinely distinct the comment in fusion_controller should say so")

    drawn.set_atoms(True)
    assert "UNDETERMINED" in drawn.win.status, "the guard that actually works"


def test_turning_advances_the_spin_and_redraws(drawn):
    from piezo1.structure.fusion_pose import SPIN_SAMPLES

    drawn.set_atoms(True)
    before = drawn.pose.spin
    coords = drawn.pose.coords.copy()
    drawn.rotate_tags()
    assert drawn.pose.spin == pytest.approx(
        (before + 2 * np.pi / SPIN_SAMPLES) % (2 * np.pi))
    assert not np.allclose(coords, drawn.pose.coords)


def test_turning_a_full_circle_returns_to_the_start(drawn):
    """A turn that drifted would make the free angle look like a trajectory."""
    from piezo1.structure.fusion_pose import SPIN_SAMPLES

    drawn.set_atoms(True)
    drawn.spin = 0.0
    drawn._draw()
    start = drawn.pose.coords.copy()
    for _ in range(SPIN_SAMPLES):
        drawn.rotate_tags()
    assert drawn.pose.coords == pytest.approx(start, abs=1e-6)


def test_turning_without_the_fold_shown_is_refused(drawn):
    """Spinning an isotropic sphere would silently do nothing at all."""
    drawn.rotate_tags()
    assert "show the tag structure first" in drawn.win.status
    assert drawn.pose is None


def test_switching_back_to_the_sphere_forgets_the_chosen_angle(drawn):
    """Otherwise a user returns to the fold at an orientation they cannot see."""
    drawn.set_atoms(True)
    drawn.rotate_tags()
    drawn.set_atoms(False)
    assert drawn.spin is None
    assert drawn.pose is None


# --------------------------------------------------------------- the dye

def test_the_dye_marker_moves_to_the_resolved_ligand_with_the_fold(drawn):
    """6U32 resolves the conjugate; with the fold drawn, that is where it goes."""
    drawn.set_dyes(True)
    at_centre = drawn.win.viewport.scene.batches["halotag:dyes"].args[0].copy()
    drawn.set_atoms(True)
    on_ligand = drawn.win.viewport.scene.batches["halotag:dyes"].args[0]

    assert len(on_ligand) == len(at_centre)
    offsets = np.linalg.norm(on_ligand - at_centre, axis=1)
    assert offsets.min() > 5.0, "the dye is not at the tag's centroid"
    assert offsets.max() < drawn.model.meta["tag_radius"] + 10.0, \
        "nor is it somewhere unrelated to the tag"
    assert (drawn.pose.ligand.sum() > 0)
    colours = drawn.win.viewport.scene.batches["halotag:fold"].args[2]
    assert (colours == np.float32(DYE_COLOR)).all(axis=1).sum() > 0


# ------------------------------------------------------------- housekeeping

def test_clearing_removes_every_batch_it_added(drawn):
    drawn.set_atoms(True)
    drawn.set_dyes(True)
    drawn.set_envelope(True)
    assert len(drawn.win.viewport.scene.batches) >= 4
    drawn.clear()
    assert drawn.win.viewport.scene.batches == {}
    assert drawn.model is None and drawn.pose is None
