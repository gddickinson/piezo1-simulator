"""Superposition, C3 symmetry recovery, and protomer correspondence."""

import numpy as np
import pytest

from piezo1.structure.superpose import (detect_c3_axis, kabsch, match_protomers,
                                        rmsd, rotation_matrix, superpose)
from conftest import protomer_blocks


def test_kabsch_recovers_a_known_transform():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(120, 3)) * 10
    rot = rotation_matrix(np.array([0.3, -0.5, 0.8]), 0.7)
    t = np.array([4.0, -2.0, 9.0])
    y = x @ rot.T + t
    fitted, err = superpose(x, y)
    assert err < 1e-8
    assert np.allclose(fitted, y, atol=1e-6)


def test_kabsch_excludes_reflections():
    rng = np.random.default_rng(1)
    x = rng.normal(size=(50, 3))
    y = x * np.array([1.0, 1.0, -1.0])      # a reflection, not a rotation
    r, _, _ = kabsch(x, y)
    assert np.linalg.det(r) > 0


def test_rmsd_is_zero_for_identical_sets():
    x = np.arange(30, dtype=float).reshape(10, 3)
    assert rmsd(x, x) == pytest.approx(0.0)


def test_c3_axis_is_exact_on_a_real_trimer(curved_structure):
    blocks, _ = protomer_blocks(curved_structure)
    assert blocks is not None
    axis = detect_c3_axis(blocks)
    assert axis.angle_deg == pytest.approx(120.0, abs=0.05)
    assert axis.rmsd < 0.5
    assert np.linalg.norm(axis.direction) == pytest.approx(1.0)


def test_c3_axis_on_synthetic_trimer():
    rng = np.random.default_rng(2)
    proto = rng.normal(size=(80, 3)) * 5 + np.array([12.0, 0.0, 0.0])
    axis_dir = np.array([0.0, 0.0, 1.0])
    blocks = [proto @ rotation_matrix(axis_dir, k * 2 * np.pi / 3).T
              for k in range(3)]
    axis = detect_c3_axis(blocks)
    assert abs(abs(float(np.dot(axis.direction, axis_dir))) - 1.0) < 1e-6
    assert axis.angle_deg == pytest.approx(120.0, abs=1e-3)


def test_match_protomers_detects_reversed_handedness(curved_structure,
                                                     flat_structure):
    """7WLT and 7WLU label their protomers in opposite rotational order.

    Superposing by chain label gives ~71 A RMSD; the correct correspondence
    gives ~19.7 A. This is a real trap, so it is pinned by a test.
    """
    a_blocks, a_res = protomer_blocks(curved_structure)
    b_blocks, b_res = protomer_blocks(flat_structure)
    common = np.array(sorted(set(a_res.tolist()) & set(b_res.tolist())))
    a = [blk[np.searchsorted(a_res, common)] for blk in a_blocks]
    b = [blk[np.searchsorted(b_res, common)] for blk in b_blocks]

    match = match_protomers(b, a)
    assert match.handedness_flipped is True
    assert match.rmsd < 25.0
    naive = match.all_rmsd[(0, 1, 2)]
    assert naive > 60.0
    assert naive / match.rmsd > 3.0
