"""The four method-level cross-checks.

Each alternative route is tested first on a case whose answer is known
analytically — otherwise a cross-check that agrees only tells you two routes
share a bug — and then against the pipeline on real data.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.crosscheck_methods import (conservation_by_kmer_anchoring,
                                                pc1_by_power_iteration,
                                                pore_radius_by_random_search,
                                                sasa_by_monte_carlo)


def _synthetic(xyz: np.ndarray, name: str = "synthetic"):
    """Minimal Structure for an analytic test case.

    Built here rather than loaded, so the expected answer is arithmetic rather
    than a previous run of the same code.
    """
    from piezo1.core.structure import Structure

    xyz = np.asarray(xyz, dtype=np.float32)
    n = len(xyz)
    return Structure(
        xyz=xyz, element=np.array(["C"] * n),
        res_name=np.array(["ALA"] * n), res_seq=np.arange(n, dtype=np.int64),
        chain=np.array(["A"] * n), atom_name=np.array(["CA"] * n),
        hetero=np.zeros(n, dtype=bool),
        b_factor=np.zeros(n, dtype=np.float32),
        occupancy=np.ones(n, dtype=np.float32),
        alt_loc=np.array([""] * n), entity=np.array(["1"] * n), name=name)


# ----------------------------------------------- the alternatives, on knowns

def test_monte_carlo_sasa_recovers_an_isolated_sphere():
    """One atom has no neighbours, so its area is exactly 4*pi*(r+probe)^2."""
    st = _synthetic(np.zeros((1, 3)), name="one")
    radius = float(st.vdw_radii()[0]) + 1.4
    assert sasa_by_monte_carlo(st, n_samples=2000) == pytest.approx(
        4.0 * np.pi * radius ** 2, rel=1e-9)


def test_power_iteration_finds_a_planted_component():
    """A matrix with one dominant direction has a known leading eigenvector."""
    rng = np.random.default_rng(0)
    direction = rng.normal(size=40)
    direction /= np.linalg.norm(direction)
    amplitudes = rng.normal(scale=8.0, size=60)
    matrix = (np.outer(amplitudes, direction)
              + rng.normal(scale=0.05, size=(60, 40)))

    vector, eigenvalue = pc1_by_power_iteration(matrix)
    assert abs(float(np.dot(vector, direction))) > 0.999
    # And it must agree with numpy's SVD, which it never calls.
    centred = matrix - matrix.mean(axis=0)
    singular = np.linalg.svd(centred, compute_uv=False)
    assert eigenvalue == pytest.approx(singular[0] ** 2 / (len(matrix) - 1),
                                       rel=1e-6)


def test_random_pore_search_finds_a_known_gap():
    """Two parallel walls with a gap: the largest fitting sphere is analytic."""
    from piezo1.structure.superpose import SymmetryAxis

    # Atoms on two planes at x = +/- 6 A, so a probe on the axis clears
    # 6 A minus the atomic radius.
    grid = np.arange(-20.0, 20.1, 2.0)
    points = [[sx * 6.0, y, z] for sx in (-1, 1) for y in grid for z in grid]
    st = _synthetic(np.array(points), name="slab")
    axis = SymmetryAxis(point=np.zeros(3), direction=np.array([0.0, 0.0, 1.0]))

    expected = 6.0 - float(st.vdw_radii()[0])
    found = pore_radius_by_random_search(st, axis, z=0.0, leash=4.0,
                                         n_samples=8000)
    assert found == pytest.approx(expected, abs=0.15)


def test_kmer_anchoring_is_exact_on_identical_sequences():
    reference = "MEPHVLGAVLYWLLLPCALLAACLLRFSGLSLVYLLFLLLLPWFPGPTRC" * 3
    values = conservation_by_kmer_anchoring(reference, [reference] * 4)
    assert np.nanmin(values) == pytest.approx(1.0)


def test_kmer_anchoring_inflates_conservation_at_variable_positions():
    """The bias that explains the whole conservation residual, in isolation.

    Anchoring by *maximum exact matches* is a selection: given a choice of
    offsets it prefers the one where residues agree. At a genuinely variable
    position that biases the reading upwards, which is why the k-mer profile has
    a floor well above zero on real data.
    """
    rng = np.random.default_rng(0)
    alphabet = "ACDEFGHIKLMNPQRSTVWY"
    reference = "".join(rng.choice(list(alphabet), 300))
    # Orthologs that are independent random sequences: true conservation is
    # near the 1/20 chance level, i.e. essentially zero on this scale.
    others = ["".join(rng.choice(list(alphabet), 300)) for _ in range(6)]
    values = conservation_by_kmer_anchoring(reference, others)
    assert np.nanmean(values) > 0.2, (
        "the bias should be visible: unrelated sequences read as partly "
        "conserved")


# --------------------------------------------------- against the pipeline

def test_sasa_routes_agree_on_a_real_domain(structure_by_id):
    from piezo1.analysis.crosscheck_methods import check_sasa
    from piezo1.analysis.measure import sasa

    st = structure_by_id("4RAX")
    if st is None:
        pytest.skip("4RAX not downloaded")
    mask = st.mask_protein() & ~st.hetero
    check = check_sasa(st, sasa(st, mask=mask).total, mask=mask)
    assert check.agrees
    assert check.relative < 0.01


def test_pca_routes_agree_exactly():
    """Power iteration must reproduce the SVD to numerical precision.

    Unlike the other three this has no approximation in it, so anything but
    exact agreement is a bug rather than a diagnosed difference.
    """
    from piezo1.analysis.crosscheck_methods import check_pca
    from piezo1.analysis.ensemble import build_ensemble

    try:
        ensemble = build_ensemble(species="mouse")
    except Exception:
        pytest.skip("ensemble unavailable")
    if len(ensemble.members) < 3:
        pytest.skip("need at least three ensemble members")

    pca = ensemble.pca()
    matrix = np.array([m.coords.ravel() for m in ensemble.members])
    check = check_pca(pca.eigenvalues[0], pca.components[0], matrix)
    assert check.relative < 1e-6
    assert "1.000000" in check.note, "the two PC1 directions must coincide"


def test_brute_force_never_beats_the_optimiser_by_much(structure_by_id):
    """The informative direction: brute force can only match or beat a local search.

    A random search finding a *larger* clearance means the pattern search stopped
    short. On 8YEZ it finds 0.978 A against the pipeline's 0.930 A — a 5.2% gap,
    which is under-convergence rather than a wrong answer. The test bounds it, so
    a real regression in the optimiser would show up as the gap widening.
    """
    from piezo1.analysis.crosscheck_methods import check_pore_radius
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    st = structure_by_id("8YEZ")
    if st is None:
        pytest.skip("8YEZ not downloaded")
    framed = apply_frame(st, canonical_transform(st))
    blocks, _ = protomer_blocks(framed)
    axis = detect_c3_axis(blocks)
    profile = pore_profile(framed, axis, step=1.0)

    check = check_pore_radius(framed, axis, profile)
    assert check.agrees
    # Brute force should be the larger of the two, and not by much.
    assert check.alternative >= check.primary - 0.02
    assert check.relative < 0.10


def test_conservation_routes_agree_in_shape(structure_by_id):
    """0.817 correlation between an alignment and a DP-free anchor.

    The residual is a bias in the *alternative* — see the module docstring — so
    this asserts the shape agrees rather than the values.
    """
    from piezo1.analysis.conservation import (OrthologSet, conservation_profile,
                                              load_orthologs)
    from piezo1.analysis.crosscheck_methods import check_conservation

    try:
        orthologs = load_orthologs()
    except Exception:
        pytest.skip("orthologs not downloaded")
    reference = orthologs.members[0]
    close = [m for m in orthologs.members[1:]
             if abs(m.length - len(reference.sequence)) <= 60][:8]
    if len(close) < 4:
        pytest.skip("too few length-matched orthologs")

    subset = OrthologSet(members=[reference] + close, meta=dict(orthologs.meta))
    primary = conservation_profile(subset, use_cache=False)
    alternative = conservation_by_kmer_anchoring(
        reference.sequence, [m.sequence for m in close], window=120)

    check = check_conservation(primary.identity, alternative)
    assert check.agrees
    assert check.alternative > 0.75

    # The bias is upward and concentrated at variable positions.
    variable = primary.identity < 0.5
    assert np.nanmean(alternative[variable]) > np.nanmean(primary.identity[variable])
    invariant = primary.identity >= 0.999
    assert np.nanmean(alternative[invariant]) > 0.98
