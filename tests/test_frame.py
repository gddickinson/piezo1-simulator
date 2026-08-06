"""Standardised framing, and the camera drift it was confused with.

Two separate faults produced the same complaint — "the structure looks
different" — and the tests keep them apart:

* the camera composed a relative orbit every time a structure loaded, so the
  *same* entry came back at a different angle on each visit;
* deposited entries sit in unrelated frames, so *different* entries never
  overlapped.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.render.camera import Camera
from piezo1.structure.frame import (ALIGNMENT_MODES, PERMUTATIONS,
                                    apply_frame, canonical_transform,
                                    reference_transform, standardise)


# --------------------------------------------------------------- the camera

def test_orbit_accumulates_but_set_orientation_does_not():
    """The bug and its fix, side by side.

    ``orbit`` composing is correct for a mouse drag and wrong for restoring a
    standard view. Loading a structure four times must not leave the camera
    somewhere new.
    """
    drifting = Camera()
    for _ in range(4):
        drifting.orbit(0.0, -0.42)

    steady = Camera()
    steady.set_orientation(0.0, -0.42)
    first = steady.rotation.copy()
    for _ in range(3):
        steady.set_orientation(0.0, -0.42)

    assert np.allclose(steady.rotation, first), "reset must be idempotent"
    assert not np.allclose(drifting.rotation, first), (
        "if orbit stopped accumulating this test is guarding nothing")


def test_set_orientation_matches_a_single_orbit_from_rest():
    a, b = Camera(), Camera()
    a.orbit(0.0, -0.42)
    b.set_orientation(0.0, -0.42)
    assert np.allclose(a.rotation, b.rotation)


# ------------------------------------------------------------- the canonical frame

def _overlap(a, b):
    """Best shared-C-alpha RMSD over both protomer correspondence classes."""
    from piezo1.structure.frame import _corresponding_ca
    best = np.inf
    for perm in PERMUTATIONS:
        moving, target, n = _corresponding_ca(a, b, perm)
        if n < 3:
            return np.nan
        best = min(best, float(np.sqrt(((moving - target) ** 2).sum(1).mean())))
    return best


def test_canonical_frame_puts_the_axis_on_z(curved_structure):
    frame = canonical_transform(curved_structure)
    placed = apply_frame(curved_structure, frame)

    ca = placed.xyz[placed.mask_ca()]
    assert np.allclose(ca.mean(axis=0)[:2], 0.0, atol=1.0), \
        "the three-fold axis should pass through the origin in x and y"

    # Recovering the axis from the placed coordinates must give +z back.
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis
    blocks, _ = protomer_blocks(placed)
    axis = detect_c3_axis(blocks)
    assert abs(abs(float(axis.direction[2])) - 1.0) < 1e-3


def test_cytosolic_end_is_placed_at_negative_z(curved_structure):
    """The sign rule, which decides whether a structure loads upside down.

    PIEZO's C terminus is intracellular by topology, so it must land at
    negative z whichever way the deposited axis happened to point.
    """
    placed = apply_frame(curved_structure,
                         canonical_transform(curved_structure))
    mask = placed.mask_ca() & ~placed.hetero
    seq, xyz = placed.res_seq[mask], placed.xyz[mask]
    cterm = xyz[seq >= seq.max() - 40][:, 2].mean()
    nterm = xyz[seq <= seq.min() + 40][:, 2].mean()
    assert cterm < 0.0, "the C-terminal (cytosolic) end must be at negative z"
    assert cterm < nterm


def test_no_downloaded_structure_loads_upside_down(structure_by_id):
    """The sign rule, over every entry rather than one.

    This failed silently on 7WLU and 11ZC while reporting a perfect C3 fit:
    the cytosolic test slice was the top 10% of residues by number, which
    straddles the extracellular cap (to ~2457 of 2547) as well as the CTD, so
    its mean z depended on which of the two happened to be better resolved.
    A structure upside down in the viewport is obvious; one upside down inside
    a calculation is not.
    """
    from piezo1.io.registry import load_registry

    checked = 0
    for entry in load_registry().entries:
        st = structure_by_id(entry.pdb)
        if st is None:
            continue
        frame = canonical_transform(st)
        if frame.mode == "deposited":       # not a trimer; nothing to orient
            continue
        placed = apply_frame(st, frame)
        mask = placed.mask_ca() & ~placed.hetero
        seq, xyz = placed.res_seq[mask], placed.xyz[mask]
        z = xyz[seq >= seq.max() - 15][:, 2].mean()
        assert z < 0.0, f"{entry.pdb} loaded upside down: C-terminus at z={z:.1f}"
        checked += 1

    if checked < 2:
        pytest.skip("need at least two downloaded trimers")


def test_canonical_frame_is_idempotent(curved_structure):
    """Framing an already-framed structure must not move it again."""
    once = apply_frame(curved_structure, canonical_transform(curved_structure))
    twice = apply_frame(once, canonical_transform(once))
    assert _overlap(once, twice) < 1e-6


def test_canonical_frame_makes_independent_structures_overlap(
        curved_structure, flat_structure):
    """The point of the feature, measured rather than asserted qualitatively.

    Two entries refined in unrelated frames start hundreds of angstroms apart.
    The canonical frame must bring them to within the least-squares optimum,
    since anything left over is real conformational difference.
    """
    before = _overlap(curved_structure, flat_structure)

    ref = apply_frame(curved_structure, canonical_transform(curved_structure))
    mob = apply_frame(flat_structure,
                      canonical_transform(flat_structure, reference=ref))
    after = _overlap(mob, ref)

    best = reference_transform(flat_structure, curved_structure).rmsd

    assert before > 50.0, "the deposited frames were already aligned?"
    assert after < before

    # The honest measure is how much of the *achievable* improvement it
    # recovers. An absolute threshold would just encode how different these two
    # conformations happen to be — 7WLT and 7WLU are the curved and flattened
    # states, so ~20 Å of the residual is real biology that no framing removes.
    recovered = (before - after) / (before - best)
    assert recovered > 0.95, (
        f"canonical framing recovered {recovered:.1%} of the achievable "
        f"improvement ({before:.1f} → {after:.1f} Å, optimum {best:.1f} Å)")


def test_a_non_trimer_is_refused_rather_than_mangled(structure_by_id):
    """4RAX is a 227-residue domain, not a trimer: there is no C3 axis to use."""
    st = structure_by_id("4RAX")
    if st is None:
        pytest.skip("4RAX not downloaded")
    frame = canonical_transform(st)
    assert frame.mode == "deposited"
    assert frame.is_identity
    assert "trimer" in frame.note


# ------------------------------------------------------- protomer correspondence

def test_reversed_protomer_labelling_is_detected(structure_by_id):
    """Deposited chain labels can run round the ring either way.

    8YFG and 8ZU3 both present chains A, B and D, but numbered in opposite
    rotational senses. Taking the labels at face value costs about 60 A of
    apparent RMSD against 8YEZ — a bookkeeping error that reads as a
    conformational change.
    """
    ref = structure_by_id("8YEZ")
    mob = structure_by_id("8YFG")
    if ref is None or mob is None:
        pytest.skip("8YEZ/8YFG not downloaded")

    from piezo1.structure.frame import _corresponding_ca
    from piezo1.structure.superpose import kabsch, rmsd

    scores = {}
    for perm in PERMUTATIONS:
        moving, target, n = _corresponding_ca(mob, ref, perm)
        rot, trans, centroid = kabsch(moving, target)
        scores[perm] = float(rmsd((moving - centroid) @ rot.T + trans, target))

    assert scores[(0, 2, 1)] < scores[(0, 1, 2)] / 3.0, \
        "8YFG should match 8YEZ only under the reversed correspondence"

    frame = reference_transform(mob, ref)
    assert frame.reordered
    assert frame.rmsd == pytest.approx(min(scores.values()), abs=1e-6)


def test_search_covers_both_correspondence_classes():
    """One representative per class is a complete search, and there are two."""
    assert len(PERMUTATIONS) == 2
    parity = [sum(1 for i in range(3) for j in range(i + 1, 3)
                  if perm[i] > perm[j]) % 2 for perm in PERMUTATIONS]
    assert set(parity) == {0, 1}, "the two must be of opposite parity"


# --------------------------------------------------------------- refusals

def test_reference_alignment_refuses_across_species(curved_structure,
                                                    flat_structure):
    """Residue 1000 is a different residue in human and mouse.

    Matching by number across species would superpose non-equivalent atoms and
    report a confident RMSD for it, so the request is refused and the caller
    falls back to a frame that does not depend on numbering.
    """
    frame = reference_transform(flat_structure, curved_structure,
                                same_numbering=False)
    assert frame.mode == "deposited"
    assert "numbering" in frame.note

    _placed, used = standardise(flat_structure, mode="reference",
                                reference=curved_structure,
                                same_numbering=False)
    assert used.mode == "canonical"
    assert "fell back" in used.note


def test_deposited_mode_changes_nothing(curved_structure):
    placed, frame = standardise(curved_structure, mode="deposited")
    assert frame.is_identity
    assert np.array_equal(placed.xyz, curved_structure.xyz)


def test_unknown_mode_is_rejected(curved_structure):
    with pytest.raises(ValueError, match="unknown alignment mode"):
        standardise(curved_structure, mode="whatever")
    assert "canonical" in ALIGNMENT_MODES


def test_every_atom_travels_including_ligands(structure_by_id):
    """Ligands, lipids and ions must move with the protein they are bound to."""
    st = structure_by_id("8YEZ")
    if st is None:
        pytest.skip("8YEZ not downloaded")
    frame = canonical_transform(st)
    placed = apply_frame(st, frame)
    assert placed.n_atoms == st.n_atoms

    # Rigid: every interatomic distance is preserved, hetero atoms included.
    rng = np.random.default_rng(0)
    idx = rng.choice(st.n_atoms, size=200, replace=False)
    before = np.linalg.norm(st.xyz[idx][:, None] - st.xyz[idx][None], axis=-1)
    after = np.linalg.norm(placed.xyz[idx][:, None] - placed.xyz[idx][None],
                           axis=-1)
    assert np.allclose(before, after, atol=1e-3)


def test_measured_geometry_is_unchanged_by_reframing(curved_structure):
    """Reframing must not move a single reported number.

    The window now standardises coordinates *before* handing them to the
    physics, so anything that framing could perturb would silently change every
    result the GUI shows relative to the CLI. A rigid transform cannot alter a
    curvature or an area, and this checks that the pipeline actually treats it
    that way — recovering the axis from the placed coordinates rather than
    carrying a stale one.
    """
    from piezo1.structure.geometry import measure_dome
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    def measure(st):
        blocks, _ = protomer_blocks(st)
        detect_c3_axis(blocks)
        dome = measure_dome(blocks, np.concatenate(blocks))
        return (dome.radius_of_curvature, dome.dome_depth, dome.dome_area,
                dome.projected_area, dome.footprint_radius)

    before = measure(curved_structure)
    after = measure(apply_frame(curved_structure,
                                canonical_transform(curved_structure)))
    for name, a, b in zip(("radius_of_curvature", "dome_depth", "dome_area",
                           "projected_area", "footprint_radius"), before, after):
        assert b == pytest.approx(a, rel=1e-6), f"{name} moved on reframing"


def test_transform_is_a_proper_rotation(curved_structure):
    """No reflections: a mirrored protein is a different molecule."""
    for frame in (canonical_transform(curved_structure),
                  reference_transform(curved_structure, curved_structure)):
        assert np.allclose(frame.rotation @ frame.rotation.T, np.eye(3), atol=1e-9)
        assert float(np.linalg.det(frame.rotation)) == pytest.approx(1.0, abs=1e-9)
