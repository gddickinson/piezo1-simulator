"""Conformational morphing between experimental endpoints."""

import numpy as np
import pytest

from piezo1.physics.anm import ANM
from piezo1.structure.morph import (interpolate_linear, morph,
                                    prepare_endpoints)
from conftest import protomer_blocks


@pytest.fixture(scope="module")
def endpoints(curved_structure, flat_structure):
    a_blocks, a_res = protomer_blocks(curved_structure)
    b_blocks, b_res = protomer_blocks(flat_structure)
    return prepare_endpoints(a_blocks, a_res, b_blocks, b_res)


def test_prepare_endpoints_matches_protomers(endpoints):
    start, end, common, info = endpoints
    assert start.shape == end.shape
    assert len(common) > 1000
    assert info["handedness_flipped"] is True     # 7WLT vs 7WLU
    assert 15.0 < info["endpoint_rmsd"] < 25.0


def test_linear_interpolation_hits_both_endpoints(endpoints):
    start, end, _, _ = endpoints
    frames = interpolate_linear(start, end, 11)
    assert np.allclose(frames[0], start)
    assert np.allclose(frames[-1], end)
    assert len(frames) == 11


def test_restrained_morph_removes_the_chord_artefact(endpoints):
    """Linear interpolation contracts bonds; the restrained method must not.

    Blades swinging through large arcs make atoms take chords through space,
    so mid-path C-alpha-C-alpha distances shrink. That is the artefact the
    restrained method exists to remove.
    """
    start, end, _, _ = endpoints
    linear = morph(start, end, n_frames=21, method="linear")
    restrained = morph(start, end, n_frames=21, method="restrained")

    assert linear.bond_error.max() > 1.0
    assert restrained.bond_error.max() < 0.1
    assert restrained.bond_error.max() < linear.bond_error.max() / 10


def test_morph_preserves_experimental_endpoints(endpoints):
    start, end, _, _ = endpoints
    tr = morph(start, end, n_frames=15, method="restrained")
    assert np.allclose(tr.frames[0], start, atol=1e-9)
    assert np.allclose(tr.frames[-1], end, atol=1e-9)


def test_modal_morph_captures_most_of_the_change(endpoints):
    start, end, common, _ = endpoints
    per = len(common)
    anm = ANM.from_trimer([start[p * per:(p + 1) * per] for p in range(3)],
                          cutoff=15.0).build()
    modes = anm.calc_modes(n_modes=30)
    tr = morph(start, end, n_frames=15, method="modal", modes=modes)
    captured = tr.meta["fraction_captured_by_modes"]
    assert 0.85 < captured <= 1.0


def test_modal_morph_rejects_a_mismatched_mode_set(endpoints, curved_structure):
    """A mode set built on a different residue subset must fail loudly."""
    start, end, _, _ = endpoints
    blocks, _ = protomer_blocks(curved_structure)      # full set, not reduced
    anm = ANM.from_trimer(blocks, cutoff=15.0).build()
    modes = anm.calc_modes(n_modes=8)
    with pytest.raises(ValueError, match="reduced residue set"):
        morph(start, end, n_frames=5, method="modal", modes=modes)


def test_trajectory_at_interpolates(endpoints):
    start, end, _, _ = endpoints
    tr = morph(start, end, n_frames=11, method="restrained")
    assert np.allclose(tr.at(0.0), start)
    assert np.allclose(tr.at(1.0), end)
    mid = tr.at(0.5)
    assert mid.shape == start.shape
    assert not np.allclose(mid, start)


def test_morph_carries_its_caveat(endpoints):
    start, end, _, _ = endpoints
    tr = morph(start, end, n_frames=5)
    assert "not a simulated trajectory" in tr.meta["note"]
