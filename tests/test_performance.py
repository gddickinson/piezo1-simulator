"""Performance work, and the guard that it changed no result.

The discipline for this round: **an optimisation that alters a number is a bug,
not a speedup.** Every change here was a reformulation with the same value —
squared distances instead of a square root, one BLAS product instead of a 3-D
broadcast, `str.split()` where it is provably equivalent — so the tests assert
*identity*, not closeness.

Timings are asserted only as loose ceilings. A test that pins a runtime fails on
a slower machine for no scientific reason; what is worth pinning is that the
fast path and the careful path agree.
"""

import time

import numpy as np
import pytest


# --------------------------------------------------------------------------
# The mmCIF tokenizer fast path
# --------------------------------------------------------------------------

def test_fast_tokenizer_path_matches_the_careful_one():
    """99.5% of mmCIF lines take the `str.split()` path. It must be exact.

    This is the function whose whitespace handling once shifted every column by
    one, so the careful path is untouched and only bypassed — and the bypass is
    checked against it here, on lines chosen to exercise the edges.
    """
    from piezo1.io.cif_reader import _tokenize

    cases = [
        "ATOM 1 N N . MET A 1 1 ? 12.3 4.5 -6.7 1.00 30.0 ? 1 MET A N 1",
        "ATOM   1  CA  CA . HOH B 2 . ? 0.0 0.0 0.0",
        "   leading and trailing whitespace   ",
        "tabs\tbetween\tfields",
        "trailing newline\n",
        "carriage return\r\n",
        "",
        "   ",
        "single",
    ]
    for line in cases:
        assert _tokenize(line) == line.split(), repr(line)


def test_tokenizer_still_handles_quotes_and_comments():
    from piezo1.io.cif_reader import _tokenize

    assert _tokenize("'quoted value' next") == ["quoted value", "next"]
    assert _tokenize('"double quoted" tail') == ["double quoted", "tail"]
    assert _tokenize("value # comment here") == ["value"]
    # A quote only terminates when followed by whitespace or end of line.
    assert _tokenize("'it's fine' after") == ["it's fine", "after"]


def test_tokenizer_agrees_with_itself_over_a_real_file(curved_structure):
    """Exhaustive rather than illustrative: every line of a deposited entry."""
    from piezo1.config import STRUCTURE_DIR
    from piezo1.io.cif_reader import _tokenize

    path = STRUCTURE_DIR / "7WLT.cif"
    if not path.exists():
        pytest.skip("7WLT not downloaded")
    lines = path.read_text(errors="replace").splitlines()
    plain = [line for line in lines
             if "'" not in line and '"' not in line and "#" not in line]
    assert len(plain) > 1000
    for line in plain:
        assert _tokenize(line) == line.split()


# --------------------------------------------------------------------------
# SASA
# --------------------------------------------------------------------------

def test_sasa_matches_the_direct_formulation(human_structure):
    """The BLAS reformulation must be arithmetically identical.

    |t - x|² is computed as |v|² + r² + 2r(p·v) rather than by forming the
    difference and taking a norm. That is exact algebra, but floating point is
    not algebra, so it is checked rather than assumed.
    """
    from scipy.spatial import cKDTree

    from piezo1.analysis.measure import _sphere_points, sasa

    subset = np.zeros(human_structure.n_atoms, bool)
    subset[:1500] = True
    result = sasa(human_structure, mask=subset)

    xyz = human_structure.xyz[subset].astype(np.float64)
    radii = human_structure.vdw_radii()[subset].astype(np.float64) + 1.4
    points = _sphere_points(256)
    tree = cKDTree(xyz)
    max_r = radii.max()

    reference = np.zeros(len(xyz))
    for i in range(len(xyz)):
        nb = np.asarray(tree.query_ball_point(xyz[i], radii[i] + max_r))
        nb = nb[nb != i]
        test = xyz[i] + points * radii[i]
        if len(nb):
            d = np.linalg.norm(test[:, None, :] - xyz[nb][None, :, :], axis=2)
            accessible = (d >= radii[nb][None, :]).all(axis=1)
        else:
            accessible = np.ones(len(test), bool)
        reference[i] = 4.0 * np.pi * radii[i] ** 2 * accessible.mean()

    assert np.array_equal(result.atom, reference), (
        f"max difference {np.abs(result.atom - reference).max():.3e}")


def test_isolated_atom_takes_the_no_neighbour_branch(human_structure):
    """The branch the optimisation rewrote: an atom with nothing near it now
    short-circuits to 4-pi-r-squared instead of building an empty comparison."""
    from piezo1.analysis.measure import sasa

    lone = human_structure.subset(np.arange(1))
    result = sasa(lone, probe=1.4, n_points=512)
    # In float64, matching what sasa() works in — vdw_radii() is float32, and
    # comparing against a float32-derived expectation only pins ~1e-7.
    radius = np.float64(lone.vdw_radii()[0]) + 1.4
    assert result.atom[0] == pytest.approx(4.0 * np.pi * radius ** 2, rel=1e-12)


# --------------------------------------------------------------------------
# Pockets
# --------------------------------------------------------------------------

def test_monte_carlo_volume_unchanged_by_the_early_exit():
    """Skipping points already inside cannot change the answer: `inside` is
    monotone, so a point once inside stays inside."""
    from piezo1.analysis.pockets import _monte_carlo_volume

    rng = np.random.default_rng(3)
    centers = rng.uniform(-20, 20, size=(600, 3))
    radii = rng.uniform(3.0, 5.5, 600)

    lo = (centers - radii[:, None]).min(axis=0)
    hi = (centers + radii[:, None]).max(axis=0)
    box = float(np.prod(hi - lo))
    pts = np.random.default_rng(0).uniform(lo, hi, size=(6000, 3))
    inside = np.zeros(6000, bool)
    for start in range(0, len(centers), 400):
        block = centers[start:start + 400]
        d = np.linalg.norm(pts[:, None, :] - block[None, :, :], axis=2)
        inside |= (d <= radii[start:start + 400][None, :]).any(axis=1)
    reference = box * float(inside.mean())

    assert _monte_carlo_volume(centers, radii) == reference


def test_union_volume_is_still_not_a_sum(human_structure):
    """Two overlapping spheres must not count their overlap twice."""
    from piezo1.analysis.pockets import _monte_carlo_volume

    centers = np.array([[0.0, 0, 0], [1.0, 0, 0]])
    radii = np.array([3.0, 3.0])
    union = _monte_carlo_volume(centers, radii, n_samples=200000)
    single = 4.0 / 3.0 * np.pi * 27.0
    assert single < union < 2 * single


# --------------------------------------------------------------------------
# Conservation cache
# --------------------------------------------------------------------------

def test_conservation_cache_returns_an_identical_profile():
    from piezo1.analysis.conservation import conservation_profile, load_orthologs

    try:
        orthologs = load_orthologs()
    except Exception as exc:
        pytest.skip(f"orthologs unavailable: {exc}")

    fresh = conservation_profile(orthologs, use_cache=False)
    cached = conservation_profile(orthologs)
    assert np.array_equal(fresh.entropy, cached.entropy)
    assert np.array_equal(fresh.coverage, cached.coverage)
    assert np.array_equal(fresh.identity, cached.identity)
    assert fresh.n_orthologs == cached.n_orthologs
    assert fresh.organisms == cached.organisms


def test_conservation_cache_key_tracks_content():
    """A cache that can go stale is worse than no cache — it would report last
    week's conservation against this week's alignment."""
    from piezo1.analysis.conservation import _profile_cache_key, load_orthologs

    try:
        orthologs = load_orthologs()
    except Exception as exc:
        pytest.skip(f"orthologs unavailable: {exc}")

    base = _profile_cache_key("MAAA", orthologs, 0.5)
    assert base == _profile_cache_key("MAAA", orthologs, 0.5)
    assert base != _profile_cache_key("MAAB", orthologs, 0.5)
    assert base != _profile_cache_key("MAAA", orthologs, 0.7)

    mutated = type(orthologs)(members=list(orthologs.members),
                              meta=dict(orthologs.meta))
    mutated.members[0] = type(mutated.members[0])(
        **{**mutated.members[0].__dict__,
           "sequence": mutated.members[0].sequence + "A"})
    assert base != _profile_cache_key("MAAA", mutated, 0.5)


# --------------------------------------------------------------------------
# Loose ceilings — a machine-independent smoke check, not a pinned runtime
# --------------------------------------------------------------------------

def test_structure_loading_is_not_pathological(human_structure):
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure

    start = time.time()
    Structure.from_file(STRUCTURE_DIR / "8YEZ.cif")
    assert time.time() - start < 5.0, "31k-atom load should take well under 5 s"


def test_sasa_is_not_pathological(human_structure):
    from piezo1.analysis.measure import sasa
    start = time.time()
    sasa(human_structure)
    assert time.time() - start < 20.0
