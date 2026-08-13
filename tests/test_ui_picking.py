"""Clicking anything drawn, and what the click is entitled to say.

Two rules carry the feature. **Nearest visible wins, whatever drew it** —
"what did I click" has one honest answer, the thing in front that is on
screen — tested on rays whose answer is known by construction. And **a
feature identifies as what it is**: a modelled tag says MODELLED, a graft
atom says PREDICTED with its pLDDT, an extra structure names itself and says
the analyses do not run on it.

Also pinned: the HETATM fix (a lipid's author number lands inside the
protein's range, and the old click confidently named a domain for it), the
visibility mask (a hidden category must stop answering clicks), and the
right-click routing through every source.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.ui.gl_widget import PRIMARY_SOURCE, nearest_hit  # noqa: E402


def _structure(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb} not downloaded; run python -m piezo1.io.fetch")
    return Structure.from_file(path)


# ------------------------------------------------------------- the geometry

ORIGIN = np.zeros(3)
FORWARD = np.array([0.0, 0.0, 1.0])


def test_the_nearer_feature_atom_beats_the_primary():
    sources = {
        PRIMARY_SOURCE: np.array([[0.0, 0.0, 10.0]]),
        "halotag": np.array([[0.5, 0.0, 5.0]]),
    }
    assert nearest_hit(sources, ORIGIN, FORWARD) == ("halotag", 0)


def test_the_nearer_primary_atom_beats_the_feature():
    sources = {
        PRIMARY_SOURCE: np.array([[0.0, 0.0, 4.0]]),
        "halotag": np.array([[0.0, 0.5, 9.0]]),
    }
    assert nearest_hit(sources, ORIGIN, FORWARD) == (PRIMARY_SOURCE, 0)


def test_an_atom_behind_the_camera_cannot_be_hit():
    sources = {"halotag": np.array([[0.0, 0.0, -5.0]])}
    assert nearest_hit(sources, ORIGIN, FORWARD) is None


def test_an_atom_off_the_ray_is_not_hit():
    sources = {"halotag": np.array([[8.0, 0.0, 5.0]])}
    assert nearest_hit(sources, ORIGIN, FORWARD) is None


def test_within_one_source_the_nearest_along_the_ray_wins():
    sources = {"x": np.array([[0.0, 0.0, 9.0], [0.1, 0.0, 3.0]])}
    assert nearest_hit(sources, ORIGIN, FORWARD) == ("x", 1)


def test_empty_sources_hit_nothing():
    assert nearest_hit({}, ORIGIN, FORWARD) is None
    assert nearest_hit({"x": np.zeros((0, 3))}, ORIGIN, FORWARD) is None


def test_a_hidden_atom_in_front_cannot_be_hit():
    """The invisible-lipid case: hiding must remove an atom from picking, or
    a click identifies something not on screen."""
    sources = {PRIMARY_SOURCE: np.array([[0.0, 0.0, 5.0], [0.3, 0.0, 9.0]])}
    masks = {PRIMARY_SOURCE: np.array([False, True])}
    assert nearest_hit(sources, ORIGIN, FORWARD, masks=masks) == \
        (PRIMARY_SOURCE, 1)
    # And the same click without the mask hits the front atom — the mask is
    # doing the work, not the geometry.
    assert nearest_hit(sources, ORIGIN, FORWARD) == (PRIMARY_SOURCE, 0)


def test_a_mask_on_one_source_does_not_shadow_another():
    sources = {
        PRIMARY_SOURCE: np.array([[0.0, 0.0, 4.0]]),
        "halotag": np.array([[0.2, 0.0, 8.0]]),
    }
    masks = {PRIMARY_SOURCE: np.array([False])}
    assert nearest_hit(sources, ORIGIN, FORWARD, masks=masks) == ("halotag", 0)


# --------------------------------------------------- what is drawn is pickable

@pytest.fixture(scope="module")
def lipid_view():
    from piezo1.render.representations import MolecularView

    st = _structure("7WLT")
    if not st.mask_ligands().any():
        pytest.skip("no ligands in this entry")
    return MolecularView(_Scene(), st, name="m")


def test_everything_is_pickable_until_something_is_hidden(lipid_view):
    assert lipid_view.pickable_mask().all()


def test_a_hidden_entity_category_stops_answering_clicks(lipid_view):
    entities = lipid_view.entity_map()
    if "lipid" not in entities.present():
        pytest.skip("no lipid category in this entry")
    lipid = entities.mask("lipid")
    lipid_view.visible_entities = frozenset(
        k for k in entities.present() if k != "lipid")
    try:
        mask = lipid_view.pickable_mask()
        assert not mask[lipid].any(), "hidden lipids still answer clicks"
        assert mask[~lipid].all(), "hiding lipids silenced something else"
    finally:
        lipid_view.visible_entities = frozenset()


def test_a_component_selection_restricts_picking_to_its_residues(lipid_view):
    st = lipid_view.structure
    residue = int(st.res_seq[st.mask_ca()][0])
    lipid_view.visible_residues = frozenset([residue])
    try:
        mask = lipid_view.pickable_mask()
        chosen = st.res_seq == residue
        assert mask[chosen].all()
        assert not mask[~chosen].any(), \
            "atoms a component hides still answer clicks"
    finally:
        lipid_view.visible_residues = None


def test_ligands_toggled_off_in_a_ribbon_style_are_unpickable(lipid_view):
    from piezo1.render.representations import Style

    st = lipid_view.structure
    lipid_view.ligands_as_spheres = False
    try:
        lipid_view.style = Style.CARTOON
        assert not lipid_view.pickable_mask()[st.mask_ligands()].any(), (
            "a ribbon style draws ligands only through the ligand pass; "
            "toggled off they are invisible and must not answer clicks")
        # In an atom style every atom is in the atoms batch regardless of the
        # ligand toggle, so the same toggle must NOT silence them there.
        lipid_view.style = Style.SPHERES
        assert lipid_view.pickable_mask()[st.mask_ligands()].all()
    finally:
        lipid_view.ligands_as_spheres = True
        lipid_view.style = Style.CARTOON


# ------------------------------------------------------- right-click routing

def test_a_right_click_through_a_feature_does_not_name_the_atom_behind(qapp_):
    """The context menu's residue entries are annotation read by primary atom
    index; when a drawn feature is nearest, the menu must open generic
    rather than describing the occluded residue behind the tag."""
    from PyQt6.QtCore import QPointF, Qt
    from PyQt6.QtGui import QMouseEvent

    from piezo1.ui.gl_widget import ViewportWidget

    view = ViewportWidget()
    got = []
    view.context_requested.connect(lambda pos, index: got.append(index))

    def click(hit):
        view.hit_at = lambda pos: hit
        press = QMouseEvent(QMouseEvent.Type.MouseButtonPress,
                            QPointF(50, 50), Qt.MouseButton.RightButton,
                            Qt.MouseButton.RightButton,
                            Qt.KeyboardModifier.NoModifier)
        release = QMouseEvent(QMouseEvent.Type.MouseButtonRelease,
                              QPointF(50, 50), Qt.MouseButton.RightButton,
                              Qt.MouseButton.NoButton,
                              Qt.KeyboardModifier.NoModifier)
        view.mousePressEvent(press)
        view.mouseReleaseEvent(release)

    click((PRIMARY_SOURCE, 7))
    click(("halotag", 3))
    click(None)
    assert got == [7, -1, -1]


@pytest.fixture(scope="module")
def qapp_():
    from PyQt6.QtWidgets import QApplication

    instance = QApplication.instance()
    if instance is None:
        try:
            instance = QApplication([])
        except Exception as exc:                       # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return instance


# --------------------------------------------------- the selection behaviour

class _Batch:
    def upload(self, *args, **kwargs):
        self.args = args


class _Scene:
    def __init__(self):
        self.batches = {}

    def spheres(self, name, **kwargs):
        return self.batches.setdefault(name, _Batch())

    cylinders = spheres

    def mesh(self, name, **kwargs):
        return self.batches.setdefault(name, _Batch())

    def remove(self, name):
        self.batches.pop(name, None)


class _Viewport:
    def __init__(self):
        self.scene = _Scene()
        self.feature_sources = {}

    def set_feature_pick_source(self, name, coords):
        if coords is None:
            self.feature_sources.pop(name, None)
        else:
            self.feature_sources[name] = coords

    def update(self):
        pass


class _MeasurePanel:
    def __init__(self):
        self.armed = False
        self.picks = []

    def add_pick(self, index, xyz, label):
        self.picks.append((index, label))


def _selection_host(structure=None, species="human"):
    from piezo1.core.annotations import load_annotations
    from piezo1.ui.selection import SelectionMixin

    class _View:
        name = "m"
        highlight = None

    class Host(SelectionMixin):
        def __init__(self):
            self.viewport = _Viewport()
            self.structure = structure
            self.view = _View() if structure is not None else None
            self.record = None
            self.measure_panel = _MeasurePanel()
            self.annotations = load_annotations(species)
            self.selected_residues = []
            self.selection_label = ""
            self._pick_hints = 0
            self.status = ""

        def _set_status(self, text):
            self.status = text

        def focus_mode(self):
            return "none"

    return Host()


def test_a_feature_pick_marks_the_atom_and_says_what_it_is():
    host = _selection_host()
    coords = np.array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
    host.register_pick_feature("halotag", coords,
                               lambda i: f"tag atom {i} — MODELLED")
    assert "halotag" in host.viewport.feature_sources

    host._on_feature_pick("halotag", 1)
    assert host.status == "tag atom 1 — MODELLED"
    marker = host.viewport.scene.batches["feature:selection"]
    assert marker.args[0] == pytest.approx(coords[[1]])


def test_an_armed_measure_click_on_a_feature_refuses_out_loud():
    """Swallowing the click silently would break inspection; measuring a
    modelled position would be a measurement of a guess. So: neither."""
    host = _selection_host()
    host.register_pick_feature("halotag", np.zeros((1, 3)),
                               lambda i: "tag atom — MODELLED")
    host.measure_panel.armed = True
    host._on_feature_pick("halotag", 0)
    assert host.measure_panel.picks == [], "a modelled position was measured"
    assert "primary structure" in host.status
    assert "MODELLED" in host.status
    assert "feature:selection" not in host.viewport.scene.batches


def test_unregistering_drops_the_source_and_its_marker():
    host = _selection_host()
    host.register_pick_feature("halotag", np.zeros((1, 3)), lambda i: "x")
    host._on_feature_pick("halotag", 0)
    assert "feature:selection" in host.viewport.scene.batches

    host.unregister_pick_feature("halotag")
    assert "halotag" not in host.viewport.feature_sources
    assert "feature:selection" not in host.viewport.scene.batches, \
        "the marker outlived the thing it marked"
    host._on_feature_pick("halotag", 0)          # stale signal: must not raise


def test_a_pick_out_of_range_is_ignored():
    host = _selection_host()
    host.register_pick_feature("halotag", np.zeros((2, 3)), lambda i: "x")
    host._on_feature_pick("halotag", 7)
    assert host.status == ""


# ------------------------------------------------------------ the HETATM fix

@pytest.fixture(scope="module")
def lipid_host():
    st = _structure("7WLT")
    if not st.hetero.any():
        pytest.skip("no HETATM records in this entry")
    return _selection_host(structure=st, species="mouse")


def test_a_lipid_click_does_not_claim_a_protein_domain(lipid_host):
    """The residue-number collision, pinned. Mouse numbering runs to 2547 and
    the lipids' author numbers land inside it; the old code looked them up."""
    st = lipid_host.structure
    hetero = np.flatnonzero(st.hetero)
    # Prefer one whose number a protein-annotation lookup would answer for.
    inside = hetero[(st.res_seq[hetero] > 0) & (st.res_seq[hetero] <= 2547)]
    index = int(inside[0] if len(inside) else hetero[0])

    lipid_host._on_pick(index)
    assert "HETATM" in lipid_host.status
    assert "not a docked pose" in lipid_host.status
    assert "domain:" not in lipid_host.status
    assert "variant" not in lipid_host.status


def test_a_protein_click_still_annotates(lipid_host):
    st = lipid_host.structure
    index = int(np.flatnonzero(~st.hetero & (st.atom_name == "CA"))[0])
    lipid_host._on_pick(index)
    assert "HETATM" not in lipid_host.status
    assert st.res_name[index] in lipid_host.status


# ------------------------------------------------------------ the controllers

class _Window:
    """A window stub that can also receive pick registrations."""

    def __init__(self, structure):
        self.structure = structure
        self.viewport = _Viewport()
        self.status = ""
        self.registered = {}

    def _set_status(self, text):
        self.status = text

    def register_pick_feature(self, name, coords, describe):
        self.registered[name] = (np.asarray(coords), describe)
        self.viewport.set_feature_pick_source(name, coords)

    def unregister_pick_feature(self, name):
        self.registered.pop(name, None)
        self.viewport.set_feature_pick_source(name, None)


@pytest.fixture(scope="module")
def fusion_window():
    from piezo1.ui.fusion_controller import FusionController

    if not (STRUCTURE_DIR / "6U32.cif").exists():
        pytest.skip("6U32 not downloaded; run python -m piezo1.io.fetch")
    window = _Window(_structure("8YEZ"))
    return window, FusionController(window)


def test_the_tag_registers_its_atoms_and_identifies_as_a_model(fusion_window):
    window, fusion = fusion_window
    fusion.clear()
    fusion.show_atoms = False
    fusion.show(True)

    coords, describe = window.registered["halotag"]
    assert len(coords) == fusion.model.n_tags
    assert "MODELLED" in describe(0)

    fusion.set_atoms(True)
    coords, describe = window.registered["halotag"]
    assert len(coords) == fusion.pose.n_atoms * fusion.pose.n_tags
    assert coords == pytest.approx(
        fusion.pose.coords.reshape(-1, 3).astype(np.float64), abs=1e-4)
    text = describe(fusion.pose.n_atoms + 3)     # an atom of the second tag
    assert "MODELLED" in text and "UNDETERMINED" in text
    assert "tag 2" in text

    fusion.clear()
    assert "halotag" not in window.registered


@pytest.fixture(scope="module")
def hybrid_window():
    from piezo1.ui.hybrid_controller import HybridController

    if not (STRUCTURE_DIR / "AF-Q92508-F1-model_v6.cif").exists():
        pytest.skip("AlphaFold model not downloaded; run python -m piezo1.io.fetch")
    window = _Window(_structure("8YEZ"))
    return window, HybridController(window)


def test_the_graft_identifies_each_half_for_what_it_is(hybrid_window):
    window, hybrid = hybrid_window
    hybrid.clear()
    hybrid.show(True)
    if hybrid.model is None:
        pytest.skip(f"could not build the model: {window.status}")

    coords, describe = window.registered["hybrid"]
    assert len(coords) == len(hybrid.model.xyz)
    predicted = int(np.flatnonzero(hybrid.model.predicted)[0])
    experimental = int(np.flatnonzero(~hybrid.model.predicted)[0])
    assert "PREDICTED" in describe(predicted)
    assert "pLDDT" in describe(predicted)
    assert "experimental" in describe(experimental)
    assert "PREDICTED" not in describe(experimental)

    hybrid.clear()
    assert "hybrid" not in window.registered


def test_a_companion_registers_and_says_who_it_is():
    from piezo1.render.representations import Style
    from piezo1.ui.companions import CompanionMixin
    from piezo1.ui.selection import SelectionMixin

    st = _structure("7WLT")

    class _FakeSettings:
        def value(self, key, default=None, type=None):  # noqa: A002
            return default

        def setValue(self, key, value):  # noqa: N802
            pass

    class _FakeRegistry:
        def get(self, pdb):
            return None

    class Host(CompanionMixin, SelectionMixin):
        def __init__(self):
            self.settings = _FakeSettings()
            self.viewport = _Viewport()
            self.record = None
            self.registry = _FakeRegistry()
            self.status = ""

        def _set_status(self, text):
            self.status = text

    host = Host()
    assert host.companion_style() is Style.BACKBONE
    host.add_companion("7WLT", structure=st, species="mouse")

    coords, describe = host._feature_registry()["extra:7WLT"]
    assert len(coords) == st.n_atoms
    text = describe(0)
    assert "extra structure 7WLT" in text
    assert "analyses run on the primary" in text

    host.remove_companion("7WLT")
    assert "extra:7WLT" not in host._feature_registry()
    assert "extra:7WLT" not in host.viewport.feature_sources
