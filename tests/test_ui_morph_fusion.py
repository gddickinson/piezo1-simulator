"""Does the HaloTag travel with the morph, and is it still a solved position?

The tag is anchored to a C-alpha of the channel, so a morph that leaves it
behind draws a fusion protein detached from the protein it is fused to. Making
it follow is easy; making it follow *honestly* is the part with tests.

The placement is not a pose. It is the centroid of the region the tag centre can
occupy without clashing, and flattening takes a third of that region away —
242 → 177 nm³ on 7WLT → 7WLU. So carrying the tag rigidly with its anchor is
wrong by about 7 Å at the far end, and the model is re-solved on every frame
instead. These tests exist to hold that distinction: that what is drawn at a
frame was solved *against that frame*, that it agrees with the entry's own model
where the two should agree, and that the far end is not quoted as the deposited
structure's answer when it is not.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")
moderngl = pytest.importorskip("moderngl")

from piezo1.config import STRUCTURE_DIR  # noqa: E402

from test_ui_morph import PAIR, app, window  # noqa: E402,F401


@pytest.fixture(scope="module")
def tagged(window, app):
    """A morph built with the HaloTag switched on, as the user would."""
    if not (STRUCTURE_DIR / "6U32.cif").exists():
        pytest.skip("6U32.cif not downloaded — run python -m piezo1.io.fetch")
    window.fusion.show(True)
    app.processEvents()
    if not window.fusion.visible:
        pytest.skip(f"fusion unavailable: {window.status_label.text()}")
    window.morph_controller.build(
        {"start": PAIR[0], "end": PAIR[1], "method": "restrained"})
    app.processEvents()
    ctl = window.morph_controller
    if ctl.trajectory is None or not ctl._fusion_frames:
        pytest.fail(f"no tagged morph: {window.status_label.text()}")
    yield ctl
    window.fusion.clear()
    app.processEvents()


def anchor_now(window):
    """The C-alpha the tag is attached to, read off the drawn coordinates."""
    st = window.view.structure
    residue = window.fusion.model.anchor_residues[0]
    chain = [c for c in st.chains
             if (st.mask_ca() & (st.chain == c)).sum() > 300][0]
    mask = st.mask_ca() & (st.chain == chain) & (st.res_seq == residue)
    return st.xyz[mask][0].astype(np.float64)


def clearance_against(centres, st) -> float:
    """Closest approach of any of the three tags to a set of coordinates.

    All three, because that is what ``min_clearance`` reports — and the three
    are not exactly C3-equivalent, since the channel's own atoms are not.
    Scoring only the first differs from it by about 0.005 A, which is small
    enough to look like a rounding problem and is not one.
    """
    from scipy.spatial import cKDTree
    protein = st.mask_protein() & ~st.hetero
    tree = cKDTree(st.xyz[protein].astype(np.float64))
    return float(min(tree.query(c)[0] for c in np.atleast_2d(centres)))


# --------------------------------------------------------------------------
# It travels, and it travels with the atom it is attached to
# --------------------------------------------------------------------------

def test_the_tag_is_drawn_at_the_anchor_it_claims(tagged, window, app):
    """The reported anchor must be the drawn C-alpha, at both ends.

    Not "near" it. The anchor is a real atom of the channel, and the whole
    fusion model hangs off it — if this drifts, every distance the tag work
    reports is measured from a point that is not on the molecule.
    """
    for fraction in (0.0, 0.5, 1.0):
        tagged.show_frame(fraction)
        app.processEvents()
        assert np.allclose(window.fusion.model.anchors[0], anchor_now(window),
                           atol=1e-3), fraction


def test_the_anchor_really_moves_along_the_path(tagged, window, app):
    """The calibration for the test above: a still anchor would pass it too."""
    tagged.show_frame(0.0)
    start = anchor_now(window)
    tagged.show_frame(1.0)
    end = anchor_now(window)
    # The C-terminus rises ~15 A towards the membrane as the dome flattens.
    assert np.linalg.norm(end - start) > 10.0


def test_the_first_frame_reproduces_the_entrys_own_model(tagged, window, app):
    """Frame 0 is the loaded structure, so its tag model must be that one's.

    Built by a different route — `build_fusion` straight off the loaded
    structure — so this compares the morph's per-frame solve against the
    ordinary one rather than against itself.
    """
    from piezo1.structure.fusion import build_fusion, load_halotag

    tagged.show_frame(0.0)
    app.processEvents()
    direct = build_fusion(window.structure, load_halotag())
    carried = window.fusion.model
    assert np.allclose(carried.tag_centres, direct.tag_centres, atol=1e-6)
    assert np.allclose(carried.anchors, direct.anchors, atol=1e-6)
    assert carried.volume.volume == pytest.approx(direct.volume.volume)


# --------------------------------------------------------------------------
# It is re-solved, not carried
# --------------------------------------------------------------------------

def test_the_tag_is_solved_against_the_frame_it_is_drawn_on(tagged, window, app):
    """The far end's tag must clear the *flattened* channel, not the curved one.

    Measured both ways round, which is what makes it a test: the same centre is
    scored against the frame it belongs to and against the start frame, and the
    two must disagree. A tag merely carried along would score identically
    against both, since nothing about it would have been recomputed.
    """
    tagged.show_frame(0.0)
    app.processEvents()
    start_structure = window.view.structure.copy_with_coords(
        window.view.structure.xyz.copy())

    tagged.show_frame(1.0)
    app.processEvents()
    centres = window.fusion.model.tag_centres
    own = clearance_against(centres, window.view.structure)
    against_start = clearance_against(centres, start_structure)

    assert own == pytest.approx(window.fusion.model.meta["min_clearance"],
                                abs=1e-6)
    assert abs(own - against_start) > 1.0


def test_a_rigid_carry_would_have_been_several_angstrom_out(tagged):
    """Why the envelope is re-solved rather than translated with the anchor.

    The anchor-to-centre offset is not constant along the path: the flattened
    channel occludes a different part of the tether's reach, so the centroid
    moves relative to its own anchor. Rigidly carrying the tag would put it
    that far from where the model says it goes.
    """
    frames = [f for f in tagged._fusion_frames if f is not None]
    first, last = frames[0], frames[-1]
    offset_start = first.tag_centres[0] - first.anchors[0]
    offset_end = last.tag_centres[0] - last.anchors[0]
    assert np.linalg.norm(offset_end - offset_start) > 3.0


def test_the_accessible_volume_shrinks_as_the_dome_flattens(tagged):
    """The measurement the re-solve exists to capture, pinned.

    242 -> 177 nm3 on this pair. Recorded as a direction and a floor rather
    than an exact value, so a change to the linker parameters moves it without
    failing, but a tag that stopped being re-solved does fail.
    """
    volumes = [f.volume.volume for f in tagged._fusion_frames if f is not None]
    assert len(volumes) == len(tagged.trajectory)
    assert volumes[0] > volumes[-1] * 1.2
    assert len(set(np.round(volumes, 3))) > len(volumes) // 2


def test_every_stored_frame_has_a_solved_model(tagged):
    assert all(f is not None for f in tagged._fusion_frames)
    assert len(tagged._fusion_frames) == len(tagged.trajectory)


# --------------------------------------------------------------------------
# What the far end is not
# --------------------------------------------------------------------------

def test_the_far_end_is_not_quoted_as_the_deposited_entrys_answer(
        tagged, window):
    """The status line must carry both numbers and say which one to use.

    The path ends on 7WLU's C-alpha positions but carries 7WLT's side chains
    and unshared residues, and the pore exit is set by whichever atom reaches
    furthest down the axis. So the tag-to-pore-exit distance there is 3.92 nm
    where the deposited entry's own model gives 3.59 — an 8% gap in the number
    the calcium work uses, which a picture would not reveal.
    """
    text = window.status_label.text()
    assert "HaloTag carried on its anchor" in text
    assert "quote the deposited value" in text
    assert tagged._fusion_endpoint is not None
    carried = tagged._fusion_frames[-1].pore_exit_distances()[0]
    # The gap is real, which is the only reason the sentence is there.
    assert abs(carried - tagged._fusion_endpoint) > 0.1
    assert f"{tagged._fusion_endpoint:.2f}" in text


def test_the_carried_distance_is_not_taken_from_the_end_entry(tagged, window):
    """The calibration: the deposited value must come from the end entry.

    Read here off 7WLU itself, by a route the controller does not use, so a
    controller that quietly reported the start entry's number would fail.
    """
    from piezo1.core.structure import Structure
    from piezo1.structure.fusion import build_fusion, load_halotag

    flat = Structure.from_file(STRUCTURE_DIR / f"{PAIR[1]}.cif")
    own = build_fusion(flat, load_halotag()).pore_exit_distances()[0]
    assert tagged._fusion_endpoint == pytest.approx(own, abs=1e-6)


def test_building_a_morph_cannot_lose_the_folds_caveat(tagged, window, app):
    """A drawn fold must never be on screen without `UNDETERMINED` said.

    The guard `test_ui_fusion` enforces everywhere else, and the morph is a new
    way to break it: the build's own status line is written last and replaces
    whatever the fusion controller had put there. So the statement has to travel
    with it.
    """
    window.fusion.set_atoms(True)
    app.processEvents()
    assert "UNDETERMINED" in window.status_label.text()

    window.morph_controller.build(
        {"start": PAIR[0], "end": PAIR[1], "method": "restrained"})
    app.processEvents()
    text = window.status_label.text()
    assert "HaloTag carried on its anchor" in text, "the morph line is the one"
    assert "UNDETERMINED" in text

    window.fusion.set_atoms(False)
    app.processEvents()
    # And the calibration: with only the sphere there is no spin to caveat, so
    # the statement must be absent rather than boilerplate.
    window.morph_controller.build(
        {"start": PAIR[0], "end": PAIR[1], "method": "restrained"})
    app.processEvents()
    assert "UNDETERMINED" not in window.status_label.text()


# --------------------------------------------------------------------------
# It reaches the screen, and it moves
# --------------------------------------------------------------------------

def test_the_tag_is_drawn_and_moves_between_the_two_ends(tagged, window, app):
    """Counted in pixels, because a batch with the right contents can draw
    nothing at all — which is how this project lost every cylinder it ever
    uploaded. The tag's own pixels are isolated by clearing it and diffing.
    """
    scene = window.viewport.scene
    ctx = scene.ctx
    size = (360, 300)
    scene.resize(*size)
    fbo = ctx.simple_framebuffer(size)

    def shot() -> np.ndarray:
        fbo.use()
        fbo.clear(0.05, 0.05, 0.07, 1.0)
        scene.render()
        px = np.frombuffer(fbo.read(components=3), np.uint8)
        return px.reshape(size[1], size[0], 3).astype(int).sum(2) > 60

    centroids = []
    for fraction in (0.0, 1.0):
        tagged.show_frame(fraction)
        app.processEvents()
        with_tag = shot()
        keys = [k for k in list(scene.batches) if k.startswith("halotag:")]
        assert keys, "nothing was uploaded for the tag"
        saved = {k: scene.batches[k] for k in keys}
        for k in keys:
            scene.remove(k)
        without = shot()
        for k, batch in saved.items():
            scene.batches[k] = batch

        only_tag = with_tag & ~without
        assert only_tag.sum() > 50, f"the tag drew {only_tag.sum()} pixels"
        ys, xs = np.nonzero(only_tag)
        centroids.append(np.array([xs.mean(), ys.mean()]))

    # It has to go somewhere: a tag left behind would land on the same pixels.
    assert np.linalg.norm(centroids[0] - centroids[1]) > 3.0


# --------------------------------------------------------------------------
# A path's tags belong to that path
# --------------------------------------------------------------------------

def test_turning_the_tag_off_drops_the_per_frame_models(tagged, window, app):
    assert tagged._fusion_frames is not None
    window.fusion.clear()
    app.processEvents()
    assert tagged._fusion_frames is None
    # And the slider still works, drawing the channel alone.
    tagged.show_frame(0.5)
    assert not window.fusion.visible

    window.fusion.show(True)
    app.processEvents()
    assert tagged._fusion_frames is not None
    assert len(tagged._fusion_frames) == len(tagged.trajectory)


def test_loading_another_structure_drops_them_too(tagged, window, app):
    other = "8YEZ" if (STRUCTURE_DIR / "8YEZ.cif").exists() else None
    if other is None:
        pytest.skip("8YEZ.cif not downloaded — run python -m piezo1.io.fetch")
    window.load_structure(other)
    app.processEvents()
    ctl = window.morph_controller
    assert ctl.trajectory is None
    assert ctl._fusion_frames is None
    assert ctl._fusion_endpoint is None
    assert window.fusion.host is None
