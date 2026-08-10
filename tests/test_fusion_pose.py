"""Drawing the real HaloTag fold without claiming a pose.

Three things have to hold for this to be an improvement on the sphere rather
than a prettier lie:

* the tag must be **moved rigidly**, so what is drawn is still the deposited
  structure and not a model of one;
* the placement must determine what it claims to determine — position and the
  seam direction — and **nothing more**, so the spin has to be demonstrably
  free rather than described as free;
* the contact counter must be **calibrated**, because it is the instrument that
  decides whether the fold fits. Its first version was wrong in exactly the way
  this project's checkers keep being wrong: it reported the placement rule's
  own construction as a finding.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core.structure import Structure
from piezo1.structure.fusion import build_fusion, load_halotag
from piezo1.structure.fusion_pose import (SPIN_SAMPLES, align_rotation,
                                          best_spin, drawable_mask, place_tag,
                                          pose_for_display, spin_scan)


@pytest.fixture(scope="module")
def tag():
    if not (STRUCTURE_DIR / "6U32.cif").exists():
        pytest.skip("6U32 not downloaded; run python -m piezo1.io.fetch")
    return load_halotag()


@pytest.fixture(scope="module")
def host():
    if not (STRUCTURE_DIR / "8YEZ.cif").exists():
        pytest.skip("8YEZ not downloaded; run python -m piezo1.io.fetch")
    return Structure.from_file(STRUCTURE_DIR / "8YEZ.cif")


@pytest.fixture(scope="module")
def model(host, tag):
    return build_fusion(host, tag)


# ------------------------------------------------- the rotation, on its own

def test_align_rotation_takes_one_vector_onto_another():
    for _ in range(20):
        rng = np.random.default_rng(_)
        a, b = rng.normal(size=3), rng.normal(size=3)
        turned = align_rotation(a, b) @ (a / np.linalg.norm(a))
        assert turned == pytest.approx(b / np.linalg.norm(b), abs=1e-12)


def test_align_rotation_is_a_rotation_not_a_reflection():
    """A reflection would also map the vector across and mirror the fold."""
    rng = np.random.default_rng(0)
    for _ in range(10):
        matrix = align_rotation(rng.normal(size=3), rng.normal(size=3))
        assert np.linalg.det(matrix) == pytest.approx(1.0, abs=1e-12)
        assert matrix @ matrix.T == pytest.approx(np.eye(3), abs=1e-12)


def test_the_antiparallel_case_has_no_shortest_arc_and_is_handled():
    """Left to the general formula this divides by a vanishing sine."""
    a = np.array([0.0, 0.0, 1.0])
    turned = align_rotation(a, -a) @ a
    assert turned == pytest.approx(-a, abs=1e-12)
    assert np.linalg.det(align_rotation(a, -a)) == pytest.approx(1.0, abs=1e-12)
    assert align_rotation(a, a) == pytest.approx(np.eye(3), abs=1e-12)


# --------------------------------------------------- what the placement does

def test_the_tag_is_moved_rigidly(host, model, tag):
    """It must remain 6U32, not a distorted picture of it.

    Every internal distance is preserved, so nothing downstream can read the
    drawn fold as a refined or relaxed structure.
    """
    keep = drawable_mask(tag.structure)
    original = tag.structure.xyz[keep].astype(float)
    placed = place_tag(host, model, tag, spin=1.1).coords[0]

    rng = np.random.default_rng(3)
    i, j = rng.integers(0, len(original), size=(2, 400))
    before = np.linalg.norm(original[i] - original[j], axis=1)
    after = np.linalg.norm(placed[i] - placed[j], axis=1)
    assert np.abs(before - after).max() < 1e-9


def test_the_tag_centre_lands_on_the_modelled_centre(host, model, tag):
    """The drawn fold must occupy the position every reported number describes.

    If the fold sat anywhere else, the tag-to-pore-exit distance in the status
    bar would be describing something other than the thing on screen.
    """
    pose = place_tag(host, model, tag, spin=0.7)
    # The same atom set `load_halotag` took the centroid of, or the two
    # centres would differ by the dye's offset and the check would be vacuous.
    keep = drawable_mask(tag.structure)
    protein = (tag.structure.mask_protein() & ~tag.structure.hetero)[keep]
    for i, centre in enumerate(model.tag_centres):
        drawn = pose.coords[i][protein].mean(axis=0)
        assert drawn == pytest.approx(np.asarray(centre), abs=1e-9)


def test_the_tag_n_terminus_faces_the_channel_c_terminus(host, model, tag):
    """The two rotations the geometry does fix.

    The placed N-terminus must sit on the line from the tag centre to PIEZO1's
    C-terminus, at the tag's own anchor-to-centre distance — any other
    orientation stretches the linker further than it needs to go.
    """
    pose = place_tag(host, model, tag, spin=2.0)
    gap = np.linalg.norm(pose.anchors[0] - model.anchors[0])
    assert gap == pytest.approx(pose.linker_gap, abs=1e-6)
    assert gap == pytest.approx(
        np.linalg.norm(model.tag_centres[0] - model.anchors[0])
        - tag.anchor_to_centre, abs=1e-6)


def test_the_spin_is_free_and_moves_nothing_it_should_not(host, model, tag):
    """The claim that the third rotation is undetermined, made checkable.

    Spinning must leave the anchor, the centre and the linker length exactly
    where they were — otherwise it is not a free angle, it is a second
    placement rule with a preference.
    """
    a = place_tag(host, model, tag, spin=0.0)
    b = place_tag(host, model, tag, spin=1.4)
    assert b.anchors[0] == pytest.approx(a.anchors[0], abs=1e-9)
    assert b.linker_gap == pytest.approx(a.linker_gap, abs=1e-9)
    assert not np.allclose(a.coords[0], b.coords[0]), "the spin did nothing"


def test_the_other_two_tags_are_c3_images_of_the_first(host, model, tag):
    """A homotrimer carries three identical copies; grid noise must not differ."""
    pose = place_tag(host, model, tag, spin=0.3)
    assert pose.n_tags == 3
    shape = [np.linalg.norm(pose.coords[i] - pose.coords[i].mean(axis=0), axis=1)
             for i in range(3)]
    assert shape[1] == pytest.approx(shape[0], abs=1e-9)
    assert shape[2] == pytest.approx(shape[0], abs=1e-9)


def test_crystallisation_solvent_is_not_drawn(tag):
    """A buffer chloride beside PIEZO1 would read as a permeating ion."""
    keep = drawable_mask(tag.structure)
    assert "HOH" not in set(tag.structure.res_name[keep])
    assert "CL" not in set(tag.structure.res_name[keep])
    assert keep.sum() > 2000, "the fold itself must survive the filter"


def test_the_bound_dye_is_kept_and_marked(host, model, tag):
    """6U32 was chosen for its covalent ligand; it is the dye position."""
    pose = place_tag(host, model, tag)
    assert pose.ligand.sum() > 10, "the TMR conjugate should be drawn"
    assert pose.ligand.sum() < 100, "only the dye, not every het atom"


# ------------------------------------------- the contact counter, calibrated

def test_a_tag_moved_far_away_touches_nothing(host, model, tag):
    """The instrument must be able to say no.

    A contact counter that never returns zero would call every orientation a
    clash and the measurement would mean nothing.
    """
    from dataclasses import replace

    far = replace(model, tag_centres=np.asarray(model.tag_centres) +
                  np.array([0.0, 0.0, 500.0]))
    assert place_tag(host, far, tag).contacts == 0


def test_a_tag_placed_inside_the_channel_touches_a_great_deal(host, model, tag):
    """And the instrument must be able to say yes, loudly."""
    from dataclasses import replace

    centre = np.asarray(host.xyz[host.mask_ca()]).mean(axis=0)
    buried = replace(model, tag_centres=np.tile(centre, (3, 1)))
    assert place_tag(host, buried, tag).body_contacts > 100


def test_the_counter_responds_to_its_threshold(host, model, tag):
    """A count that ignored the distance would be reading something else."""
    close = place_tag(host, model, tag, spin=1.0, contact=2.0).contacts
    wide = place_tag(host, model, tag, spin=1.0, contact=6.0).contacts
    assert wide > close


def test_the_attachment_residue_is_excluded_and_that_is_what_mattered():
    """The defect that made the first version manufacture a finding.

    The placement points the tag's N-terminal residue at PIEZO1's C-terminus,
    so of course they touch. Counting that contact said 0 of 36 orientations
    clear on every structure — the fold appearing to contradict the sphere
    envelope, when it was only reporting the rule's own construction.

    This is checked on synthetic coordinates rather than on the real pair, so
    it fails if the exclusion is dropped even should the deposited structures
    change.
    """
    from piezo1.structure.fusion_pose import TagPose

    touching = np.array([True, True, False, False])
    body = np.array([False, False, True, True])
    pose = TagPose(coords=np.zeros((1, 4, 3)), radii=np.ones(4),
                   ligand=np.zeros(4, bool), body=body, anchors=np.zeros((1, 3)),
                   seams=np.zeros((1, 2, 3)), spin=0.0, linker_gap=1.0,
                   touching=touching, contact_distance=3.4)
    assert pose.contacts == 2
    assert pose.attachment_contacts == 2
    assert pose.body_contacts == 0
    assert pose.clears, "only the anchor residue touches; the fold body is clear"


def test_the_anchor_residue_the_exclusion_drops_is_the_one_the_model_anchors_on(tag):
    """The exclusion must be exactly one residue, and the right one."""
    keep = drawable_mask(tag.structure)
    residues = tag.structure.res_seq[keep]
    protein = ~tag.structure.hetero[keep]
    excluded = set(residues[residues == residues[protein].min()].tolist())
    assert len(excluded) == 1
    anchor_residue = tag.structure.res_seq[
        tag.structure.mask_ca() & ~tag.structure.hetero].min()
    assert excluded == {anchor_residue}


# ------------------------------------------------------- the measured result

def test_the_scan_covers_the_whole_turn(host, model, tag):
    counts = spin_scan(host, model, tag)
    assert len(counts) == SPIN_SAMPLES
    angle, same = best_spin(host, model, tag)
    assert same == pytest.approx(counts)
    assert counts[int(round(angle / (2 * np.pi / SPIN_SAMPLES)))] == counts.min()


def test_the_fold_needs_more_room_than_the_sphere_allows(host, model, tag):
    """The number this module exists to produce.

    `accessible_volume` treats the tag as a sphere of its radius of gyration
    and states in its own docstring that the real fold, reaching 30 A, will
    clash where the sphere does not. Measured on 8YEZ: the sphere sits 21.5 A
    clear, and just 1 of 36 orientations of the real fold avoids the channel.
    """
    assert tag.max_extent > tag.radius_of_gyration * 1.6
    assert not model.meta["clashes"], "the sphere model reports 8YEZ as clear"
    counts = spin_scan(host, model, tag)
    assert (counts == 0).sum() == 1, (
        f"8YEZ: expected exactly one clear orientation, got {(counts == 0).sum()}")
    assert counts.max() >= 20, "and the worst orientation is well inside"


def test_the_two_models_agree_on_which_structures_admit_a_tag():
    """The cross-check that makes the fold worth drawing.

    Sphere and fold are independent statements about the same question. The
    one entry whose sphere clearance falls below the radius of gyration is the
    one where no orientation of the real fold clears — and the three where the
    sphere fits all admit at least one. Were this to fail, one of the two is
    wrong and the drawing would be the less trustworthy of them.
    """
    tag = load_halotag()
    verdicts = {}
    for pdb in ("7WLT", "8YFG", "8YEZ", "11ZC"):
        path = STRUCTURE_DIR / f"{pdb}.cif"
        if not path.exists():
            pytest.skip(f"{pdb} not downloaded; run python -m piezo1.io.fetch")
        host = Structure.from_file(path)
        model = build_fusion(host, tag)
        verdicts[pdb] = (bool(model.meta["clashes"]),
                         int((spin_scan(host, model, tag) == 0).sum()))

    for pdb, (sphere_clashes, clear) in verdicts.items():
        assert sphere_clashes == (clear == 0), (
            f"{pdb}: sphere says clash={sphere_clashes} but {clear} of "
            f"{SPIN_SAMPLES} fold orientations clear the channel")
    assert verdicts["11ZC"] == (True, 0)
    assert verdicts["7WLT"][1] > verdicts["8YEZ"][1], (
        "the sphere is generous about how much room there is")


def test_the_display_pose_reports_how_much_choice_there_was(host, model, tag):
    """A drawn fold that happens to fit must not imply that it had to."""
    pose = pose_for_display(host, model, tag)
    assert pose.clears, "the least-contacting draw is the one to show"
    assert pose.meta["clear_spins"] == 1
    assert pose.meta["spins_sampled"] == SPIN_SAMPLES
    assert pose.meta["most_contacts"] > pose.meta["fewest_contacts"]
    assert "UNDETERMINED" in pose.summary()


def test_an_explicit_spin_is_honoured_even_when_it_is_a_bad_one(host, model, tag):
    """Turning the tag must show what the user asked for, not the best draw."""
    worst = int(np.argmax(spin_scan(host, model, tag)))
    angle = worst * 2 * np.pi / SPIN_SAMPLES
    pose = pose_for_display(host, model, tag, spin=angle)
    assert pose.spin == pytest.approx(angle)
    assert not pose.clears
    assert pose.meta["clear_spins"] == 1, "the scan is still reported honestly"
