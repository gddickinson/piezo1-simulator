"""Every wired parameter, proved to reach a number by moving it.

Round 49 showed that reading the source is not enough: `pore.step` was read,
declared and displayed, and changing it did nothing. Round 49b wired the
remaining 21 parameters, and this is the proof obligation that goes with the
wiring — override it, show the answer moves, show reset restores it exactly.

**A probe has to be calibrated before a zero means anything.** Two of these
first reported "no effect" for reasons that had nothing to do with wiring: the
pockets probe used coordinates too diffuse to make any alpha sphere at all, and
`pockets.r_max` was pushed in the direction that cannot change the answer. Both
are noted where they are set, because a badly aimed probe produces exactly the
result a broken wire would.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.parameter_effect import measure_effect, probe_effects
from piezo1.parameters import PARAMETERS


@pytest.fixture(scope="module")
def probes():
    from piezo1.analysis.allostery import build_network
    from piezo1.analysis.conservation import ConservationProfile
    from piezo1.analysis.measure import sasa
    from piezo1.analysis.pockets import alpha_spheres
    from piezo1.analysis.validation import (bootstrap_cliffs_delta,
                                            permutation_test)
    from piezo1.core.structure import Structure

    rng = np.random.default_rng(0)
    coords = rng.normal(scale=12.0, size=(90, 3))
    dcc = np.abs(np.corrcoef(rng.normal(size=(90, 40))))

    # Dense enough that the DEFAULT r_min/r_max and burial filter already find
    # spheres. Diffuse coordinates find none, and every pockets probe then
    # reads "no effect" whether or not the parameter is wired.
    blob = rng.uniform(-7.0, 7.0, size=(200, 3))

    a = np.array([1.0, 2, 3, 4, 5, 6, 7, 8.5])
    b = np.array([2.0, 3, 4, 5, 6, 7, 8, 9.5])

    n = 40
    atoms = Structure(
        name="probe", xyz=rng.normal(scale=8.0, size=(n, 3)),
        element=np.array(["C"] * n), atom_name=np.array(["CA"] * n),
        res_name=np.array(["ALA"] * n), res_seq=np.arange(1, n + 1),
        chain=np.array(["A"] * n), hetero=np.zeros(n, bool),
        b_factor=np.zeros(n), occupancy=np.ones(n),
        alt_loc=np.array([""] * n), entity=np.zeros(n, int))

    profile = ConservationProfile(
        residues=np.arange(1, 51), entropy=rng.random(50),
        identity=rng.random(50), coverage=np.linspace(0.4, 1.0, 50),
        n_orthologs=9)

    return {
        "allostery.contact_cutoff":
            (lambda: float(build_network(coords, dcc).nnz), 4.0),
        "allostery.min_correlation":
            (lambda: float(build_network(coords, dcc).sum()), 0.2),
        "pockets.r_min":
            (lambda: float(alpha_spheres(blob).radii.size), 4.6),
        # 4.0 rather than a larger value: the widest sphere this geometry makes
        # is ~5.0 A, so RAISING r_max above 5.5 admits nothing and the probe
        # would report no effect on a parameter that is correctly wired.
        "pockets.r_max":
            (lambda: float(alpha_spheres(blob).radii.size), 4.0),
        "pockets.min_neighbours":
            (lambda: float(alpha_spheres(blob).radii.size), 2),
        "pockets.neighbour_radius":
            (lambda: float(alpha_spheres(blob).radii.size), 30.0),
        "sasa.probe_radius": (lambda: float(sasa(atoms).total), 3.0),
        "sasa.n_points": (lambda: float(sasa(atoms).total), 64),
        "stats.n_permutations":
            (lambda: float(permutation_test(a, b).p_value), 500),
        "stats.n_bootstrap":
            (lambda: float(bootstrap_cliffs_delta(a, b).ci_low), 400),
        "conservation.min_coverage":
            (lambda: float(len(profile.top_conserved(n=10))), 0.99),
    }


def test_every_probed_parameter_reaches_the_number(probes):
    """The proof obligation Round 49b was set: it moves, and it restores."""
    failures = [e.summary() for e in probe_effects(probes) if not e.ok]
    assert not failures, "\n".join(failures)


def test_no_registered_parameter_is_read_by_nothing():
    """Round 49 measured 26 dead parameters; 49b wired them all."""
    from piezo1.analysis.provenance_chain import unwired_parameters

    dead = unwired_parameters()
    assert dead == [], f"{len(dead)} parameters are read by no code: {dead}"


def test_integer_parameters_are_converted_not_passed_as_floats():
    """A defect the wiring introduced, and the probes caught.

    ``value()`` returns a float, so ``n_permutations`` arrived as ``10000.0``
    and numpy rejected it. Counts must be cast where they are resolved.
    """
    from piezo1.analysis.validation import permutation_test

    result = permutation_test(np.array([1.0, 2, 3, 4]), np.array([2.0, 3, 4, 5]))
    assert 0.0 <= result.p_value <= 1.0

    from piezo1.analysis.measure import sasa
    from piezo1.core.structure import Structure

    n = 12
    st = Structure(
        name="t", xyz=np.random.default_rng(1).normal(scale=4.0, size=(n, 3)),
        element=np.array(["C"] * n), atom_name=np.array(["CA"] * n),
        res_name=np.array(["ALA"] * n), res_seq=np.arange(1, n + 1),
        chain=np.array(["A"] * n), hetero=np.zeros(n, bool),
        b_factor=np.zeros(n), occupancy=np.ones(n),
        alt_loc=np.array([""] * n), entity=np.zeros(n, int))
    assert sasa(st).total > 0


# --------------------------------------------------- the instrument itself

def test_a_parameter_with_no_effect_is_reported_as_such():
    """Calibration: the measurement must be able to say "no".

    Every "moves" above is worthless unless a genuinely inert parameter comes
    back as inert, so one is constructed here.
    """
    effect = measure_effect("pore.step", lambda: 42.0, override=0.25)
    assert not effect.moved
    assert effect.restored
    assert not effect.ok


def test_a_probe_that_raises_is_reported_not_swallowed():
    def broken():
        raise ValueError("probe failed")

    effect = measure_effect("pore.step", broken, override=0.25)
    assert "probe failed" in effect.error
    assert not effect.ok


def test_the_registry_is_restored_even_when_the_probe_raises():
    """A leaked override would corrupt every later measurement in the process."""
    before = PARAMETERS.value("pore.step")

    def broken():
        raise RuntimeError("boom")

    measure_effect("pore.step", broken, override=0.25)
    assert PARAMETERS.value("pore.step") == before
    assert not PARAMETERS.overrides()


def test_an_existing_override_is_preserved():
    """Measuring one parameter must not silently reset a user's other overrides."""
    from piezo1.parameters import reset, set_value

    try:
        set_value("anm.cutoff", 13.0)
        measure_effect("pore.step", lambda: 1.0, override=0.25)
        assert PARAMETERS.value("anm.cutoff") == pytest.approx(13.0)
    finally:
        reset()


# ------------------------------------------- the bug the proof uncovered

def test_top_conserved_excludes_residues_that_fail_the_coverage_filter():
    """Found by probing ``conservation.min_coverage``.

    Sorting a failing residue to the bottom is not excluding it: when fewer
    than ``n`` residues met the requirement, the rest were returned anyway,
    carrying their real conservation value with nothing to mark them. This is
    reachable from the CLI as ``conservation --top``.
    """
    from piezo1.analysis.conservation import ConservationProfile

    rng = np.random.default_rng(0)
    profile = ConservationProfile(
        residues=np.arange(1, 51), entropy=rng.random(50),
        identity=rng.random(50), coverage=np.linspace(0.4, 1.0, 50),
        n_orthologs=9)

    passing = int((profile.coverage >= 0.99).sum())
    got = profile.top_conserved(n=10, min_coverage=0.99)
    assert len(got) == passing == 1, (
        "top_conserved returned residues that fail the coverage requirement")

    coverage = dict(zip(profile.residues.tolist(), profile.coverage.tolist()))
    assert all(coverage[r] >= 0.99 for r, _ in got)

    # And it still returns n when n residues genuinely qualify.
    assert len(profile.top_conserved(n=10, min_coverage=0.0)) == 10
