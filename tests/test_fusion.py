"""HaloTag fusion geometry — a model, tested as one.

There is no structure of the PIEZO1–HaloTag fusion, so nothing here can be
checked against a measurement of the fusion itself. What *can* be checked is
that the construction is internally sound: that the symmetry it claims is real,
that the envelope is an envelope rather than a pose, that the geometry it
reports does not depend on quantities it should not depend on, and that the two
sign conventions which produced confident wrong answers stay fixed.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core.structure import Structure
from piezo1.parameters import PARAMETERS
from piezo1.structure.frame import apply_frame, canonical_transform
from piezo1.structure.fusion import (HALOTAG_PDB, accessible_volume,
                                     build_fusion, cterm_anchors, load_halotag)

PORE_MOUTH_RADIUS = PARAMETERS.value("fusion.pore_mouth_radius")


@pytest.fixture(scope="module")
def tag():
    if not (STRUCTURE_DIR / f"{HALOTAG_PDB}.cif").exists():
        pytest.skip(f"{HALOTAG_PDB} not downloaded")
    return load_halotag()


@pytest.fixture(scope="module")
def framed(human_structure):
    return apply_frame(human_structure, canonical_transform(human_structure))


# ------------------------------------------------------------------ the tag

def test_halotag_geometry_matches_the_deposited_structure(tag):
    """The measured inputs the whole model rests on.

    6U32 is 1.8 Å with the TMR-HaloTag ligand covalently bound. These are
    measurements, not choices, so they are pinned tightly — a change here means
    the file changed or the selection did.
    """
    assert tag.radius_of_gyration == pytest.approx(17.6, abs=0.3)
    assert tag.max_extent == pytest.approx(30.0, abs=1.0)
    # The fusion joins PIEZO1's C-terminus to the tag's N-terminus, so it is
    # this distance — not the radius of gyration — that sets where the body sits.
    assert tag.anchor_to_centre == pytest.approx(19.9, abs=0.5)
    assert tag.anchor_to_ligand == pytest.approx(21.8, abs=0.5)
    assert tag.ligand is not None, "6U32's TMR ligand should be resolved"


def test_anchor_is_the_n_terminus_not_the_c_terminus(tag):
    """A C-terminal fusion attaches to the tag's N-terminus.

    Anchoring at the wrong end would move the tag body by the difference
    between the two offsets and still look entirely plausible.
    """
    ca = tag.structure.mask_ca() & ~tag.structure.hetero
    numbers = tag.structure.res_seq[ca]
    n_term = tag.structure.xyz[ca & (tag.structure.res_seq == numbers.min())][0]
    assert np.allclose(tag.anchor, n_term)


# --------------------------------------------------------------- the anchors

def test_anchors_are_the_last_residue_of_each_protomer(framed):
    anchors, residues = cterm_anchors(framed)
    assert len(anchors) == 3
    assert residues == (2521, 2521, 2521), "8YEZ resolves to human 2521"

    # C3-related, so all three sit at the same height and the same radius.
    assert np.ptp(anchors[:, 2]) < 1e-3
    radial = np.hypot(anchors[:, 0], anchors[:, 1])
    assert np.ptp(radial) < 1e-3


# ------------------------------------------------ the envelope, not a pose

def test_envelope_is_a_region_not_a_point(framed, tag):
    volume = accessible_volume(framed, cterm_anchors(framed)[0][0], tag)
    assert len(volume) > 1000
    assert volume.volume > 10.0
    # It must actually be three-dimensional; a degenerate result would still
    # have a centroid and would still look like an answer.
    assert np.ptp(volume.points, axis=0).min() > 5.0


def test_every_accessible_point_clears_the_channel(framed, tag):
    """The clash filter must be exact, not approximate."""
    from scipy.spatial import cKDTree

    anchor = cterm_anchors(framed)[0][0]
    volume = accessible_volume(framed, anchor, tag, linker_residues=5)
    protein = framed.mask_protein() & ~framed.hetero
    tree = cKDTree(framed.xyz[protein].astype(np.float64))

    required = volume.tag_radius + PARAMETERS.value("fusion.clash_clearance")
    nearest = tree.query(volume.points)[0]
    assert nearest.min() >= required - 1e-6

    # ...and within the tether.
    assert np.linalg.norm(volume.points - anchor, axis=1).max() <= volume.reach + 1e-6


def test_a_zero_length_linker_cannot_place_the_tag(framed, tag):
    """A 33 kDa tag welded straight onto the C-terminus does not fit.

    The failure is raised rather than returned as an empty envelope, because an
    empty envelope has a centroid and would flow downstream as a position.
    """
    with pytest.raises(RuntimeError, match="no clash-free position"):
        build_fusion(framed, tag, linker_residues=0)


# --------------------------------------------------- the roadmap's criteria

def test_c3_symmetry_is_preserved_across_the_three_tags(framed, tag):
    """Validation criterion 1.

    Exact by construction — one envelope is solved and rotated — but measured
    anyway, since a placement that broke the symmetry would look right in a
    picture and be wrong in every calculation downstream.
    """
    model = build_fusion(framed, tag)
    assert model.n_tags == 3
    assert model.c3_deviation() < 1e-6

    # Same height, same radius, 120 deg apart.
    assert np.ptp(model.tag_centres[:, 2]) < 1e-6
    angles = np.degrees(np.arctan2(model.tag_centres[:, 1],
                                   model.tag_centres[:, 0]))
    gaps = np.diff(np.sort(angles))
    assert np.allclose(gaps, 120.0, atol=1e-3)


def test_no_steric_clash_with_the_channel(framed, tag):
    """Validation criterion 3."""
    model = build_fusion(framed, tag)
    assert not model.meta["clashes"]
    assert model.meta["min_clearance"] > model.meta["tag_radius"]


def test_tag_sits_below_the_predicted_window_but_the_band_is_reachable(
        framed, tag):
    """Validation criterion 2 — the one that does **not** pass.

    The roadmap predicted a tag centre 4–6 nm from the pore exit. The ensemble
    mean is 3.9 nm on 8YEZ and 3.3–4.2 nm across all twenty entries, so the
    centroid sits below the window. That prediction came from adding the tag's
    ~2 nm anchor-to-centre offset to the anchor's 2.6 nm from the pore exit,
    which implicitly assumes the tag points straight away from the channel;
    averaged over the accessible directions, many of which run sideways along
    the membrane, the mean is pulled in.

    The band is not unreachable — about half the envelope lies inside it — so
    the recorded result is "the window describes an achievable position, not
    the ensemble mean", and this test pins both halves of that.
    """
    model = build_fusion(framed, tag)
    centre = model.pore_exit_distances()
    assert np.allclose(centre, centre[0]), "C3 mates must be equidistant"
    assert centre[0] == pytest.approx(3.95, abs=0.25)
    assert centre[0] < 4.0, "if this passes the window, update the record"

    reachable = model.volume.distances_from(model.pore_exit) / 10.0
    assert reachable.max() > 6.0
    in_band = float(((reachable >= 4.0) & (reachable <= 6.0)).mean())
    assert in_band > 0.3, "the predicted window should be well populated"


# ------------------------------------------------- the two sign conventions

def test_pore_exit_is_on_the_cytosolic_side(framed, tag):
    """`SymmetryAxis.direction` has no fixed sign, and trusting it was wrong.

    Re-detecting the axis on an already-framed structure returns −z as readily
    as +z; it does for 7WLT and 8YFG but not for 8YEZ. That put their pore exit
    at the *extracellular* end and the tag 15–16 nm from it, against 3.9 nm for
    the same construct on 8YEZ.
    """
    model = build_fusion(framed, tag)
    # In the canonical frame the cytosolic side is −z, and the pore exit must be
    # below the C-terminal anchors it opens beside.
    assert model.pore_exit[2] < 0.0
    assert model.pore_exit[2] < model.anchors[:, 2].max()
    # It is a pore exit, not a blade tip: it lies on the conduction axis.
    assert np.hypot(model.pore_exit[0], model.pore_exit[1]) < 1e-6


def test_pore_exit_ignores_off_axis_blade_tips(framed, tag):
    """The other trap: the lowest atom overall is a blade tip, not the mouth."""
    model = build_fusion(framed, tag)
    protein = framed.mask_protein() & ~framed.hetero
    xyz = framed.xyz[protein]

    lowest_anywhere = xyz[:, 2].min()
    near_axis = np.hypot(xyz[:, 0], xyz[:, 1]) <= PORE_MOUTH_RADIUS
    lowest_near_axis = xyz[near_axis][:, 2].min()

    assert model.pore_exit[2] == pytest.approx(lowest_near_axis, abs=1e-3)
    assert lowest_anywhere <= lowest_near_axis


def test_geometry_is_consistent_across_every_downloaded_trimer(structure_by_id,
                                                              tag):
    """The same construct on twenty entries should not give twenty answers.

    Both sign faults showed up exactly here — as two structures reporting four
    times the distance the others did — so the sweep is the regression test.
    """
    from piezo1.io.registry import load_registry

    distances, checked = [], 0
    for entry in load_registry().entries:
        st = structure_by_id(entry.pdb)
        if st is None:
            continue
        st = apply_frame(st, canonical_transform(st))
        try:
            model = build_fusion(st, tag)
        except (ValueError, RuntimeError):
            continue                     # not a trimer, or no room for the tag
        assert model.c3_deviation() < 1e-6, entry.pdb
        distances.append(model.pore_exit_distances()[0])
        checked += 1

    if checked < 5:
        pytest.skip("need several downloaded trimers")
    distances = np.array(distances)
    assert distances.min() > 2.5 and distances.max() < 5.0, (
        f"spread {distances.min():.2f}-{distances.max():.2f} nm suggests a "
        f"frame or sign fault, not real structural variation")


# ---------------------------------------------------------- the assumption

def test_reported_distance_is_robust_to_the_unverified_linker(framed, tag):
    """The linker length is the one assumed input, so its influence is bounded.

    A thirtyfold change in accessible volume moves the reported centre by under
    a nanometre — and *downwards*, because a longer tether wraps further around
    the channel. So the miss on the 4–6 nm criterion cannot be explained away by
    the assumption, which is what makes it worth recording.
    """
    centres, volumes = [], []
    for n in (2, 5, 10, 20, 30):
        model = build_fusion(framed, tag, linker_residues=n)
        centres.append(model.pore_exit_distances()[0])
        volumes.append(model.volume.volume)

    assert volumes[-1] / volumes[0] > 30.0, "the sweep should span a wide range"
    assert max(centres) - min(centres) < 1.0
    assert centres[-1] < centres[0], "a longer tether should not reach further out"
    assert all(c < 4.3 for c in centres)


def test_parameters_are_registered_not_hard_coded():
    """Every number the model uses must be listed, with a stated source."""
    for key in ("fusion.linker_residues", "fusion.residue_extension",
                "fusion.grid_spacing", "fusion.clash_clearance",
                "fusion.pore_mouth_radius"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, f"{key} is not registered"
        assert parameter.citation, key
        assert parameter.description, key

    # The linker is the assumption, and must say so rather than look measured.
    linker = PARAMETERS.get("fusion.linker_residues")
    assert linker.citation == "unverified"
    assert linker.source_note


def test_overriding_the_linker_takes_effect_at_call_time(framed, tag):
    """Resolved in the body, not baked into a default at import."""
    default = build_fusion(framed, tag).volume.reach
    PARAMETERS.set_value("fusion.linker_residues", 20)
    try:
        overridden = build_fusion(framed, tag).volume.reach
    finally:
        PARAMETERS.reset("fusion.linker_residues")
    assert overridden > default
    assert build_fusion(framed, tag).volume.reach == pytest.approx(default)
