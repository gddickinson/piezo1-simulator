"""The per-feature rendering styles, and the caveats they must not move.

The primary structure has had a style selector since the beginning and the
superposition overlay since the Overlay panel existed; the HaloTag fold, the
full-length graft, the extra structures, the component highlight and the
resolved ligands were each drawn in exactly one hard-coded representation.
These tests cover the controls that opened them up — and, more importantly,
the property every one of them must preserve: **restyling is presentation
only**. The fold keeps its UNDETERMINED status line in every style, the graft
keeps its confidence colouring and its seam, the highlight keeps its gold and
its residues, and nothing computed moves at all.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402


class _Batch:
    def __init__(self):
        self.args = ()

    def upload(self, *args, **kwargs):
        self.args = args


class _Scene:
    """Records what was drawn, without a GL context. Includes ``mesh`` —
    the ribbon styles are what several of these tests are about."""

    def __init__(self):
        self.batches = {}

    def spheres(self, name, **kwargs):
        return self.batches.setdefault(name, _Batch())

    def cylinders(self, name, **kwargs):
        return self.batches.setdefault(name, _Batch())

    def mesh(self, name, **kwargs):
        return self.batches.setdefault(name, _Batch())

    def remove(self, name):
        self.batches.pop(name, None)

    def get(self, name):
        return self.batches.get(name)

    def set_visible(self, name, visible):
        pass


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


def _structure(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb} not downloaded; run python -m piezo1.io.fetch")
    return Structure.from_file(path)


# ---------------------------------------------------------- the HaloTag fold

@pytest.fixture(scope="module")
def fusion():
    from piezo1.ui.fusion_controller import FusionController

    if not (STRUCTURE_DIR / "6U32.cif").exists():
        pytest.skip("6U32 not downloaded; run python -m piezo1.io.fetch")
    window = _Window(_structure("8YEZ"))
    return FusionController(window)


@pytest.fixture
def fold(fusion):
    fusion.clear()
    fusion.show_atoms = False
    fusion.show_dyes = False
    fusion.show_envelope = False
    fusion.spin = None
    fusion.fold_style = "atoms"
    fusion.show(True)
    fusion.set_atoms(True)
    return fusion


def test_the_default_fold_style_is_the_pinned_sphere_cloud(fold):
    """`test_ui_fusion.py` pins the sphere cloud's batch and radii; the style
    machinery must leave that default bit-for-bit alone."""
    batches = fold.win.viewport.scene.batches
    assert fold.fold_style == "atoms"
    assert "halotag:fold" in batches
    assert not any(k.startswith("halotag:fold:") for k in batches)


def test_a_cartoon_fold_is_a_ribbon_and_not_the_sphere_cloud(fold):
    fold.set_fold_style("cartoon")
    batches = fold.win.viewport.scene.batches
    assert "halotag:fold:ribbon" in batches
    assert "halotag:fold" not in batches
    positions = batches["halotag:fold:ribbon"].args[0]
    assert len(positions) > 100, "three tags' worth of ribbon geometry"
    # The dye would vanish with the side chains; the ligand pass keeps it.
    assert "halotag:fold:ligands" in batches


def test_every_fold_style_keeps_the_undetermined_caveat(fold):
    """The guard the fold exists behind, checked across the whole style list.

    A cartoon looks even more like a determined pose than the sphere cloud
    does, so the caveat matters more after restyling, not less.
    """
    from piezo1.ui.fusion_controller import FOLD_STYLES

    for key, _label in FOLD_STYLES:
        fold.set_fold_style(key)
        assert "UNDETERMINED" in fold.win.status, key
        assert "halotag:seam" in fold.win.viewport.scene.batches, key


def test_the_fold_colours_survive_restyling(fold):
    """TAG, DYE and CONTACT colours carry the reported numbers; a style change
    must not hand them to a palette that knows nothing of them."""
    from piezo1.structure.fusion_pose import SPIN_SAMPLES, spin_scan
    from piezo1.ui.fusion_controller import CONTACT_COLOR, DYE_COLOR, TAG_COLOR

    counts = spin_scan(fold.win.structure, fold.model)
    fold.set_fold_style("ball_and_stick")
    fold.spin = float(np.argmax(counts) * 2 * np.pi / SPIN_SAMPLES)
    fold._draw()

    batches = fold.win.viewport.scene.batches
    assert "halotag:fold:atoms" in batches and "halotag:fold:bonds" in batches
    colours = batches["halotag:fold:atoms"].args[2]
    assert (colours == np.float32(TAG_COLOR)).all(axis=1).sum() > 1000
    assert (colours == np.float32(DYE_COLOR)).all(axis=1).sum() > 0
    red = (colours == np.float32(CONTACT_COLOR)).all(axis=1)
    assert red.sum() == fold.pose.body_contacts * fold.pose.n_tags


def test_the_restyled_fold_has_every_placed_atom(fold):
    fold.set_fold_style("sticks")
    coords = fold.win.viewport.scene.batches["halotag:fold:atoms"].args[0]
    assert len(coords) == fold.pose.n_atoms * fold.pose.n_tags
    assert coords == pytest.approx(
        fold.pose.coords.reshape(-1, 3).astype(np.float32))


def test_switching_back_to_atom_spheres_restores_the_pinned_batch(fold):
    fold.set_fold_style("cartoon")
    fold.set_fold_style("atoms")
    batches = fold.win.viewport.scene.batches
    assert "halotag:fold" in batches
    assert "halotag:fold:ribbon" not in batches


def test_an_unknown_fold_style_is_refused(fold):
    fold.set_fold_style("wireframe")
    assert fold.fold_style == "atoms"


def test_a_fold_style_chosen_before_the_fold_is_shown_is_stored(fusion):
    fusion.clear()
    fusion.show_atoms = False
    fusion.fold_style = "atoms"
    fusion.show(True)
    fusion.set_fold_style("tube")
    assert fusion.fold_style == "tube"
    assert "halotag:tags" in fusion.win.viewport.scene.batches, \
        "the radius-of-gyration sphere is not a style and must stay a sphere"
    fusion.set_atoms(True)
    assert "halotag:fold:ribbon" in fusion.win.viewport.scene.batches
    fusion.clear()


# ------------------------------------------------------ the full-length model

@pytest.fixture(scope="module")
def hybrid():
    from piezo1.ui.hybrid_controller import HybridController

    if not (STRUCTURE_DIR / "AF-Q92508-F1-model_v6.cif").exists():
        pytest.skip("AlphaFold model not downloaded; run python -m piezo1.io.fetch")
    window = _Window(_structure("8YEZ"))
    return HybridController(window)


@pytest.fixture
def drawn_hybrid(hybrid):
    hybrid.clear()
    hybrid.style = "spheres"
    hybrid.show(True)
    if hybrid.model is None:
        pytest.skip(f"could not build the model: {hybrid.win.status}")
    return hybrid


def test_the_hybrid_model_records_its_c_alphas(drawn_hybrid):
    model = drawn_hybrid.model
    assert model.ca is not None
    assert 1000 < int(model.ca.sum()) < len(model.xyz)
    assert model.ca[model.predicted].any() and model.ca[~model.predicted].any()


def test_the_hybrid_tube_keeps_the_seam_and_the_status_line(drawn_hybrid):
    drawn_hybrid.set_style("tube")
    batches = drawn_hybrid.win.viewport.scene.batches
    assert "hybrid:ribbon" in batches
    assert "hybrid:atoms" not in batches
    assert "hybrid:seam" in batches
    assert "PREDICTED" in drawn_hybrid.win.status


def test_the_hybrid_ribbon_keeps_the_confidence_colouring(drawn_hybrid):
    """Grey stays grey and the bands stay bands: restyling must not soften
    the one signal that says most of the blade is a prediction."""
    from piezo1.ui.hybrid_controller import EXPERIMENTAL_COLOR

    drawn_hybrid.set_style("backbone")
    colours = drawn_hybrid.win.viewport.scene.batches["hybrid:ribbon"].args[2]
    grey = np.isclose(colours, np.float32(EXPERIMENTAL_COLOR),
                      atol=1e-4).all(axis=1)
    assert grey.any(), "the experimental part lost its flat grey"
    assert not grey.all(), "the graft lost its pLDDT bands"
    assert len({tuple(np.round(c, 3)) for c in colours[~grey]}) > 1


def test_the_hybrid_default_is_still_the_sphere_cloud(drawn_hybrid):
    drawn_hybrid.set_style("tube")
    drawn_hybrid.set_style("spheres")
    batches = drawn_hybrid.win.viewport.scene.batches
    assert "hybrid:atoms" in batches
    assert "hybrid:ribbon" not in batches


def test_an_unknown_hybrid_style_is_refused(drawn_hybrid):
    drawn_hybrid.set_style("wireframe")
    assert drawn_hybrid.style in {"spheres", "tube", "backbone"}


# ------------------------------------------------------- the extra structures

class _FakeSettings:
    """QSettings' value/setValue signature, backed by a dict."""

    def __init__(self):
        self.stored = {}

    def value(self, key, default=None, type=None):  # noqa: A002
        return self.stored.get(key, default)

    def setValue(self, key, value):  # noqa: N802
        self.stored[key] = value


class _FakeCompanionView:
    def __init__(self):
        self.style = None
        self.rebuilt = 0

    def rebuild(self):
        self.rebuilt += 1


def _companion_host():
    from piezo1.ui.companions import CompanionMixin

    class Host(CompanionMixin):
        def __init__(self):
            self.settings = _FakeSettings()
            self.viewport = _Viewport()
            self.status = ""

        def _set_status(self, text):
            self.status = text

    return Host()


def test_companions_default_to_backbone():
    from piezo1.render.representations import Style

    host = _companion_host()
    assert host.companion_style() is Style.BACKBONE


def test_setting_a_companion_style_restyles_every_companion():
    from piezo1.render.representations import Style
    from piezo1.ui.companions import Companion

    host = _companion_host()
    views = [_FakeCompanionView(), _FakeCompanionView()]
    for pdb, view in zip(("7WLT", "7WLU"), views, strict=True):
        host._companions()[pdb] = Companion(
            pdb=pdb, structure=None, view=view, color=(1, 0, 0),
            species="mouse")

    host.set_companion_style("cartoon")
    assert host.companion_style() is Style.CARTOON
    for view in views:
        assert view.style is Style.CARTOON and view.rebuilt == 1

    host.set_companion_style("not-a-style")
    assert host.companion_style() is Style.CARTOON, \
        "an unknown key must not clobber the stored choice"


# ---------------------------------------------------- the component highlight

@pytest.fixture(scope="module")
def component():
    from piezo1.ui.component_controller import ComponentController

    window = _Window(_structure("8YEZ"))
    controller = ComponentController(window)
    controller.show("pore_module")
    if controller.selection is None or not controller.selection.highlight.any():
        pytest.skip("no highlighted residues on this entry")
    return controller


def test_the_highlight_styles_change_radii_and_nothing_else(component):
    scene = component.win.viewport.scene

    component.set_style("ball_and_stick")
    from piezo1.ui.component_controller import (
        BOND_RADIUS,
        HIGHLIGHT_COLOR,
        HIGHLIGHT_RADIUS,
    )
    coords_bs = scene.batches["component:atoms"].args[0].copy()
    assert scene.batches["component:atoms"].args[1][0] == HIGHLIGHT_RADIUS
    assert "component:bonds" in scene.batches

    component.set_style("sticks")
    assert scene.batches["component:atoms"].args[1][0] == BOND_RADIUS
    assert "component:bonds" in scene.batches

    component.set_style("spheres")
    radii = scene.batches["component:atoms"].args[1]
    assert radii.max() > 1.0, "van der Waals, not ball-sized"
    assert "component:bonds" not in scene.batches, \
        "spheres with sticks poking out would read as a fourth style"

    # Which residues are highlighted is curated, not styled.
    assert scene.batches["component:atoms"].args[0] == pytest.approx(coords_bs)
    colours = scene.batches["component:atoms"].args[2]
    assert (colours == np.float32(HIGHLIGHT_COLOR)).all()


def test_an_unknown_highlight_style_is_refused(component):
    component.set_style("ball_and_stick")
    component.set_style("wireframe")
    assert component.highlight_style == "ball_and_stick"


# ------------------------------------------------------- the resolved ligands

@pytest.fixture(scope="module")
def ligand_view():
    from piezo1.render.representations import MolecularView, Style

    st = _structure("7WLT")
    if not st.mask_ligands().any():
        pytest.skip("no ligands in this entry")
    scene = _Scene()
    view = MolecularView(scene, st, name="m", style=Style.CARTOON)
    return view


def test_ligands_default_to_van_der_waals_spheres(ligand_view):
    ligand_view.ligand_style = "spheres"
    ligand_view.rebuild()
    radii = ligand_view.scene.batches["m:ligands"].args[1]
    assert radii.max() > 1.2
    assert "m:ligbonds" not in ligand_view.scene.batches


def test_ligand_ball_and_stick_draws_bonds_in_their_own_batch(ligand_view):
    from piezo1.render.representations import BALL_RADIUS

    ligand_view.ligand_style = "ball_and_stick"
    ligand_view.rebuild()
    batches = ligand_view.scene.batches
    assert (batches["m:ligands"].args[1] == np.float32(BALL_RADIUS)).all()
    assert "m:ligbonds" in batches
    starts = batches["m:ligbonds"].args[0]
    assert len(starts) > 0, "lipids have bonds; none were found"

    ligand_view.clear()
    assert "m:ligbonds" not in batches, "clear() left the ligand bonds behind"


def test_an_unknown_ligand_style_falls_back_to_spheres(ligand_view):
    ligand_view.ligand_style = "wireframe"
    ligand_view.rebuild()
    radii = ligand_view.scene.batches["m:ligands"].args[1]
    assert radii.max() > 1.2


def test_a_colour_override_wins_in_every_representation(ligand_view):
    """The mechanism the fold and the graft rely on, checked at the layer
    that implements it."""
    st = ligand_view.structure
    override = np.tile(np.float32((0.1, 0.9, 0.2)), (st.n_atoms, 1))
    ligand_view.color_override = override
    try:
        assert ligand_view.atom_colors() is not None
        assert (ligand_view.atom_colors() == override).all()
        ligand_view.ligand_style = "spheres"
        ligand_view.rebuild()
        colours = ligand_view.scene.batches["m:ligands"].args[2]
        assert (colours == np.float32((0.1, 0.9, 0.2))).all()
    finally:
        ligand_view.color_override = None
