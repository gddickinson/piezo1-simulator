"""Clicking anything drawn, and what the click is entitled to say.

Before this, the pick source was the primary structure's atoms and nothing
else: the HaloTag, the extra structures and the full-length graft were mute,
and a click aimed at them identified whatever primary atom lay behind. Two
rules carry the feature. **Nearest wins, whatever drew it** — "what did I
click" has one honest answer, the thing in front — tested here on rays whose
answer is known by construction. And **a feature identifies as what it is**:
a modelled tag says MODELLED, a graft atom says PREDICTED with its pLDDT, an
extra structure names itself and says the analyses do not run on it. A tag
atom identified like a deposited one would be the confident wrong answer the
rest of this project spends its guards on.

Also pinned: the HETATM fix. A lipid's author-assigned residue number can
land inside the protein's range, and the old click looked it up in the
curated annotation — confidently naming a domain the lipid is not part of.
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
