"""Elastic network model: correctness, symmetry, and the gating-overlap result."""

import numpy as np
import pytest

from piezo1.physics.anm import ANM, build_hessian
from piezo1.structure.superpose import match_protomers, superpose
from conftest import protomer_blocks


def _connected_cloud(n=120, spacing=6.0, seed=5):
    """A compact, definitely-connected bead network."""
    rng = np.random.default_rng(seed)
    side = int(np.ceil(n ** (1 / 3)))
    grid = np.stack(np.meshgrid(*[np.arange(side)] * 3, indexing="ij"), -1)
    pts = grid.reshape(-1, 3)[:n].astype(float) * spacing
    return pts + rng.normal(scale=0.3, size=pts.shape)


def test_hessian_is_symmetric_and_has_six_zero_modes():
    coords = _connected_cloud()
    h = build_hessian(coords, cutoff=15.0)
    assert h.shape == (3 * len(coords), 3 * len(coords))
    dense = h.toarray()
    assert np.allclose(dense, dense.T, atol=1e-10)
    # A connected network has exactly six rigid-body modes at zero.
    vals = np.linalg.eigvalsh(dense)
    assert np.abs(vals[:6]).max() < 1e-8
    assert vals[6] > 1e-6


def test_disconnected_network_is_detected():
    """Two well-separated clusters give twelve zero modes, not six.

    If this were not accounted for, six rigid-body motions would be returned as
    though they were the lowest functional modes.
    """
    a = _connected_cloud(n=60, seed=11)
    b = _connected_cloud(n=60, seed=12) + np.array([500.0, 0.0, 0.0])
    anm = ANM(np.vstack([a, b]), cutoff=15.0).build()
    assert anm.n_components() == 2
    modes = anm.calc_modes(n_modes=6)
    assert modes.meta["zero_modes_dropped"] == 12
    assert modes.eigenvalues.min() > 1e-8


def test_translation_is_an_exact_null_vector():
    coords = _connected_cloud(n=60, seed=6)
    h = build_hessian(coords, cutoff=14.0).toarray()
    for axis in range(3):
        t = np.zeros(3 * len(coords))
        t[axis::3] = 1.0
        assert np.abs(h @ t).max() < 1e-9


def test_spring_models_all_build():
    coords = _connected_cloud(n=80, seed=7)
    for spring in ("uniform", "inverse_square", "inverse_sixth"):
        h = build_hessian(coords, cutoff=15.0, spring=spring)
        assert h.nnz > 0


def test_modes_on_real_trimer_are_symmetry_labelled(curved_structure):
    blocks, _ = protomer_blocks(curved_structure)
    anm = ANM.from_trimer(blocks, cutoff=15.0, spring="inverse_square").build()
    modes = anm.calc_modes(n_modes=12)
    anm.label_symmetry(modes)

    assert modes.n_modes == 12
    assert (modes.eigenvalues > 0).all()
    assert np.all(np.diff(modes.eigenvalues) >= -1e-12)   # sorted

    # Characters must land on the group-theoretic values.
    for ch, sym in zip(modes.character, modes.symmetry):
        expected = 1.0 if sym == "A" else -0.5
        assert ch == pytest.approx(expected, abs=0.02)

    # E modes come in degenerate pairs.
    e_idx = [i for i, s in enumerate(modes.symmetry) if s == "E"]
    assert len(e_idx) >= 2
    for i, j in zip(e_idx[::2], e_idx[1::2]):
        assert modes.eigenvalues[i] == pytest.approx(modes.eigenvalues[j], rel=1e-3)


def test_lowest_symmetric_mode_predicts_the_gating_transition(
        curved_structure, flat_structure):
    """The central scientific claim, pinned as a regression test.

    An ANM built from the closed structure alone should reproduce the observed
    curved-to-flattened change, and by symmetry the overlap must live entirely
    in the A (three-fold symmetric) modes.
    """
    a_blocks, a_res = protomer_blocks(curved_structure)
    b_blocks, b_res = protomer_blocks(flat_structure)
    common = np.array(sorted(set(a_res.tolist()) & set(b_res.tolist())))
    a = [blk[np.searchsorted(a_res, common)] for blk in a_blocks]
    b = [blk[np.searchsorted(b_res, common)] for blk in b_blocks]

    match = match_protomers(b, a)
    curved = np.vstack(a)
    flat = np.vstack([b[i] for i in match.order])
    fitted, trimer_rmsd = superpose(flat, curved)
    disp = fitted - curved
    assert 15.0 < trimer_rmsd < 25.0

    anm = ANM.from_trimer(a, cutoff=15.0, spring="inverse_square").build()
    modes = anm.calc_modes(n_modes=40)
    anm.label_symmetry(modes)

    ov = modes.overlap(disp)
    is_a = modes.symmetry == "A"

    assert ov.max() > 0.6, f"best overlap only {ov.max():.3f}"
    assert modes.symmetry[int(np.argmax(ov))] == "A"
    assert modes.cumulative_overlap(disp)[-1] > 0.9
    # Symmetry selection rule: E modes cannot describe a C3-symmetric change.
    assert ov[~is_a].max() < 0.02
    assert (ov[is_a] ** 2).sum() / (ov ** 2).sum() > 0.99


def test_msf_and_collectivity_are_sane(curved_structure):
    blocks, _ = protomer_blocks(curved_structure)
    anm = ANM.from_trimer(blocks, cutoff=15.0).build()
    modes = anm.calc_modes(n_modes=10)
    msf = modes.msf()
    assert msf.shape == (anm.n_sites,)
    assert (msf > 0).all()
    for i in range(modes.n_modes):
        k = modes.collectivity(i)
        assert 0.0 < k <= 1.0
