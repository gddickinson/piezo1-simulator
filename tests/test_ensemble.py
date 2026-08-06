"""PCA over the experimental structure ensemble, and its agreement with the ANM."""

import numpy as np
import pytest

from piezo1.analysis.ensemble import (DEFAULT_EXCLUSIONS, build_ensemble,
                                      rwsip, subspace_overlap)
from piezo1.physics.anm import ANM


@pytest.fixture(scope="module")
def ensemble():
    try:
        return build_ensemble(species="mouse", min_common=900)
    except ValueError as exc:
        pytest.skip(f"ensemble unavailable: {exc}")


@pytest.fixture(scope="module")
def pca(ensemble):
    return ensemble.pca()


@pytest.fixture(scope="module")
def modes(ensemble):
    """An elastic network built on the ensemble's own residue basis.

    It must be the *same* basis, or the components and modes live in different
    spaces and the overlap is meaningless. The code raises rather than
    broadcasting if they disagree.
    """
    ref = ensemble.members[0].coords
    per = len(ensemble.residues)
    anm = ANM.from_trimer([ref[i * per:(i + 1) * per] for i in range(3)],
                          cutoff=15.0).build()
    m = anm.calc_modes(n_modes=30)
    anm.label_symmetry(m)
    return m


# --------------------------------------------------------------------------
# Construction
# --------------------------------------------------------------------------

def test_ensemble_is_built_on_a_shared_basis(ensemble):
    assert len(ensemble) >= 8
    assert len(ensemble.residues) > 800
    assert ensemble.n_sites == len(ensemble.residues) * 3
    assert ensemble.coords.shape == (len(ensemble), ensemble.n_sites, 3)
    assert ensemble.meta["numbering"] == "human Q92508"


def test_paralogue_and_known_bad_entries_are_excluded(ensemble):
    """6KG7 is PIEZO2 and must never enter a PIEZO1 ensemble."""
    assert "6KG7" in ensemble.excluded
    assert "PIEZO2" in ensemble.excluded["6KG7"]
    assert "6KG7" not in ensemble.labels
    for pdb in DEFAULT_EXCLUSIONS:
        assert pdb not in ensemble.labels


def test_reversed_protomer_order_is_detected(ensemble):
    """Several deposited entries label protomers in the opposite sense.

    If this ever returns none, correspondence matching has silently stopped
    working and every downstream number is suspect.
    """
    flipped = [m.pdb for m in ensemble.members if m.handedness_flipped]
    assert flipped, "expected at least one entry with reversed protomer order"


def test_superposition_leaves_sane_rmsds(ensemble):
    rmsds = [m.rmsd_to_reference for m in ensemble.members]
    assert min(rmsds) == 0.0                    # the reference itself
    assert max(rmsds) < 40.0
    curved = [m.rmsd_to_reference for m in ensemble.members if m.state == "curved"]
    flat = [m.rmsd_to_reference for m in ensemble.members
            if m.state in ("flat", "flattened")]
    if curved and flat:
        assert max(curved) < min(flat), "flat states must differ more from curved"


# --------------------------------------------------------------------------
# The principal components
# --------------------------------------------------------------------------

def test_pc1_dominates(pca):
    assert pca.variance_explained[0] > 0.7
    assert pca.cumulative_variance()[-1] == pytest.approx(1.0, abs=1e-9)
    assert np.all(np.diff(pca.eigenvalues) <= 1e-9)      # sorted


def test_pc1_is_the_gating_coordinate(pca):
    """PC1 must order the structures by gating state.

    No state labels go into the PCA — it sees only coordinates. That the first
    principal component nevertheless separates curved from intermediate from
    flat, in that order, is the result: the dominant axis of experimental
    variation *is* the gating transition.
    """
    proj = pca.projections[:, 0]
    by_state: dict[str, list[float]] = {}
    for state, p in zip(pca.states, proj):
        by_state.setdefault(state, []).append(float(p))

    assert "curved" in by_state
    curved = np.mean(by_state["curved"])
    for flat_state in ("flat", "flattened"):
        if flat_state in by_state:
            assert np.mean(by_state[flat_state]) > curved
    if "intermediate" in by_state:
        mid = np.mean(by_state["intermediate"])
        flat = np.mean([v for s in ("flat", "flattened")
                        for v in by_state.get(s, [])] or [mid + 1])
        assert curved < mid < flat, "the intermediate state should sit between"


def test_per_site_displacement_is_the_interpretable_number(pca):
    """The 3N-space norm is a big number that means nothing per atom."""
    assert pca.amplitude(0) > pca.rms_displacement(0)
    assert 1.0 < pca.rms_displacement(0) < 30.0


def test_pc1_is_collective(pca):
    assert pca.collectivity(0) > 0.3


# --------------------------------------------------------------------------
# Agreement with the elastic network
# --------------------------------------------------------------------------

def test_pc1_matches_a_low_frequency_mode(pca, modes):
    """The headline comparison: does the ANM predict the observed variation?"""
    j, overlap = pca.best_mode_for(modes, pc=0)
    assert overlap > 0.6, f"PC1 best overlap only {overlap:.3f}"
    assert j < 12, "should match a low-frequency mode, not an obscure one"
    assert pca.cumulative_overlap_with_modes(modes, pc=0)[-1] > 0.85


def test_leading_pcs_match_symmetric_modes(pca, modes):
    """A modes must win despite E modes outnumbering them two to one.

    Isotropic tension is C3-symmetric, so only A modes can couple to it. The
    experimental ensemble varies along exactly those directions, which is the
    symmetry selection rule showing up in the deposited structures rather than
    in a single pairwise transition.
    """
    assert modes.symmetry is not None
    n_a = int((modes.symmetry == "A").sum())
    n_e = int((modes.symmetry == "E").sum())
    assert n_e > n_a, "E modes should be the majority, making this a real test"

    best = [pca.best_mode_for(modes, pc=i)[0] for i in range(3)]
    labels = [modes.symmetry[j] for j in best]
    assert labels.count("A") >= 2, f"got {labels}"


def test_overlap_matrix_shape_and_bounds(pca, modes):
    ov = pca.overlap_with_modes(modes, n_pcs=4, n_modes=10)
    assert ov.shape == (4, 10)
    assert (ov >= 0).all() and (ov <= 1.0 + 1e-9).all()


def test_mismatched_basis_raises_rather_than_broadcasting(pca, modes):
    class Fake:
        n_modes = 2
        vectors = np.zeros((2, 5, 3))
    with pytest.raises(ValueError, match="residue basis"):
        pca.overlap_with_modes(Fake())


def test_subspace_overlap_bounds_and_control(pca, modes):
    same = subspace_overlap(pca.components[:5], pca.components[:5])
    assert same == pytest.approx(1.0, abs=1e-9)
    real = subspace_overlap(pca.components[:10], modes.vectors[:10])
    rng = np.random.default_rng(0)
    control = subspace_overlap(rng.normal(size=(10, pca.components.shape[1], 3)),
                               modes.vectors[:10])
    assert real > 20 * control, f"real {real:.3f} vs random {control:.4f}"


def test_rwsip_is_amplitude_weighted(pca, modes):
    w_pca = pca.eigenvalues[:8]
    # 1/lambda, because fluctuation amplitude scales that way. Passing raw
    # eigenvalues would weight the stiffest modes most, which is backwards.
    w_anm = 1.0 / np.maximum(modes.eigenvalues[:8], 1e-12)
    score = rwsip(pca.components, w_pca, modes.vectors, w_anm, n=8)
    assert 0.0 <= score <= 1.0
    assert score > 0.3
    identical = rwsip(pca.components, w_pca, pca.components, w_pca, n=8)
    assert identical == pytest.approx(1.0, abs=1e-9)
