"""Independent re-derivations of the headline results.

Round 18's lesson generalised: the useful check is the one that does not reuse
the derivation being checked. A test written from the same understanding as the
code shares its blind spots — if the author misread the physics, the test
encodes the misreading.

Every route here answers a question the main pipeline already answers, by
different means. Two of them found something.
"""

import numpy as np
import pytest
from dataclasses import dataclass

from piezo1.analysis.crosscheck import (CrossCheck,
                                        dome_curvature_by_cap_geometry,
                                        dome_curvature_by_parabola,
                                        gating_overlap_by_distances,
                                        t50_by_ode_integration,
                                        t50_by_steady_state)
from piezo1.physics.kinetics import GatingModel
from piezo1.structure.geometry import fit_sphere


@dataclass
class _Axis:
    """Minimal axis for synthetic tests, matching SymmetryAxis's interface."""

    point: np.ndarray
    direction: np.ndarray

    def project(self, xyz):
        return (np.asarray(xyz) - self.point) @ self.direction

    def radial(self, xyz):
        delta = np.asarray(xyz) - self.point
        along = np.outer(delta @ self.direction, self.direction)
        return np.linalg.norm(delta - along, axis=1)


def _cap(radius: float, fraction: float, n: int = 900, seed: int = 0):
    """A spherical cap of known radius, apex at the origin."""
    rng = np.random.default_rng(seed)
    a = fraction * radius
    r = np.sqrt(rng.uniform(0.0, a * a, n))
    theta = rng.uniform(0.0, 2 * np.pi, n)
    z = np.sqrt(radius * radius - r * r) - radius
    points = np.column_stack([r * np.cos(theta), r * np.sin(theta), z])
    slope = a / np.sqrt(radius * radius - a * a)
    return points, slope


AXIS = _Axis(np.zeros(3), np.array([0.0, 0.0, 1.0]))


# --------------------------------------------------------------------------
# Dome curvature — three routes on synthetic caps of known radius
# --------------------------------------------------------------------------

def test_sphere_fit_is_exact_at_every_slope():
    """The route the pipeline uses, validated against ground truth.

    This is the point of the synthetic caps: the real dome has no known answer,
    so exactness has to be shown where one exists.
    """
    for fraction in (0.15, 0.5, 0.894):
        points, _ = _cap(100.0, fraction)
        assert fit_sphere(points, iterations=4).radius == pytest.approx(
            100.0, rel=1e-3)


def test_the_exact_cap_route_is_also_exact_at_every_slope():
    for fraction in (0.15, 0.5, 0.894):
        points, _ = _cap(100.0, fraction)
        assert dome_curvature_by_cap_geometry(points, AXIS) * 10.0 == \
            pytest.approx(100.0, rel=2e-3)


def test_the_parabola_route_degrades_with_slope():
    """The finding. The parabola is a shallow-cap approximation, and its error
    grows exactly where PIEZO1 sits."""
    errors = {}
    for fraction in (0.15, 0.5, 0.894):
        points, slope = _cap(100.0, fraction)
        estimate = dome_curvature_by_parabola(points, AXIS) * 10.0
        errors[round(slope, 2)] = abs(estimate - 100.0) / 100.0
    slopes = sorted(errors)
    assert errors[slopes[0]] < 0.02, "should be accurate for a shallow cap"
    assert errors[slopes[-1]] > 0.20, "should fail badly at 63 degrees"
    assert errors[slopes[0]] < errors[slopes[1]] < errors[slopes[-1]]


def test_the_parabola_always_underestimates():
    """A systematic bias, not noise — the signature of a dropped term."""
    for fraction in (0.3, 0.5, 0.7, 0.894):
        points, _ = _cap(100.0, fraction)
        assert dome_curvature_by_parabola(points, AXIS) * 10.0 < 100.0


def test_the_two_valid_routes_agree_on_the_real_dome(curved_structure):
    """Sphere fit against the exact cap relation, on 7WLT."""
    from piezo1.structure.geometry import measure_dome
    from piezo1.structure.superpose import detect_c3_axis
    from piezo1.ui.model_utils import protomer_blocks

    from test_geometry import _tm_surface

    blocks, _ = protomer_blocks(curved_structure)
    surface = _tm_surface(curved_structure, "mouse")
    axis = detect_c3_axis(blocks)

    primary = measure_dome(blocks, surface).radius_of_curvature / 10.0
    alternative = dome_curvature_by_cap_geometry(surface, axis)
    assert abs(alternative - primary) / primary < 0.10, (
        f"{primary:.3f} vs {alternative:.3f}")
    # Both land inside the Round 29 bootstrap interval.
    assert 8.8 < primary < 10.4 and 8.8 < alternative < 10.4


# --------------------------------------------------------------------------
# Mode overlap without superposition
# --------------------------------------------------------------------------

def test_distance_route_is_invariant_to_rigid_motion():
    """The property that makes it an independent check: rotating one structure
    must not change the answer, because a superposition never happens."""
    rng = np.random.default_rng(0)
    closed = rng.normal(size=(120, 3)) * 10.0
    displacement = rng.normal(size=(120, 3)) * 0.5
    open_ = closed + displacement

    class Modes:
        n_modes = 1
        symmetry = np.array(["A"])
        vectors = displacement[None, :, :]

    direct = gating_overlap_by_distances(closed, open_, Modes(), n_pairs=800)

    angle = 0.7
    rotation = np.array([[np.cos(angle), -np.sin(angle), 0],
                         [np.sin(angle), np.cos(angle), 0], [0, 0, 1.0]])
    moved = open_ @ rotation.T + np.array([31.0, -12.0, 7.0])
    rotated = gating_overlap_by_distances(closed, moved, Modes(), n_pairs=800)
    assert rotated == pytest.approx(direct, abs=1e-9)


def test_distance_route_recovers_a_planted_mode():
    rng = np.random.default_rng(1)
    closed = rng.normal(size=(150, 3)) * 10.0
    true_mode = rng.normal(size=(150, 3))
    open_ = closed + 0.4 * true_mode
    decoy = rng.normal(size=(150, 3))

    class Modes:
        n_modes = 2
        symmetry = np.array(["A", "A"])
        vectors = np.stack([decoy, true_mode])

    assert gating_overlap_by_distances(closed, open_, Modes(),
                                       n_pairs=1500) > 0.9
    assert gating_overlap_by_distances(closed, open_, Modes(), mode_index=0,
                                       n_pairs=1500) < 0.4


def test_distance_route_confirms_the_gating_overlap(curved_structure,
                                                    flat_structure):
    """The check that matters: a superposition-free route agrees, so the
    Kabsch fit and the protomer matching are not manufacturing the result."""
    from piezo1.physics.anm import ANM
    from piezo1.structure.superpose import kabsch, match_protomers
    from piezo1.ui.model_utils import protomer_blocks

    _c, cr = protomer_blocks(curved_structure)
    _f, fr = protomer_blocks(flat_structure)
    common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))

    def blocks(structure):
        out = []
        for chain in structure.chains:
            mask = structure.mask_ca() & (structure.chain == chain)
            if mask.sum() < 300:
                continue
            index = {int(r): i for i, r in enumerate(structure.res_seq[mask])}
            xyz = structure.xyz[mask]
            if all(r in index for r in common):
                out.append(np.array([xyz[index[r]] for r in common], float))
        return out[:3]

    cb, fb = blocks(curved_structure), blocks(flat_structure)
    fb = [fb[i] for i in match_protomers(cb, fb).order]
    anm = ANM.from_trimer(cb, cutoff=15.0, spring="inverse_square").build()
    modes = anm.calc_modes(n_modes=30)
    anm.label_symmetry(modes)

    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    displacement = ((((np.vstack(fb) - centroid) @ rotation.T + translation)
                     - np.vstack(cb))).ravel()
    values = np.abs(np.asarray(modes.overlap(displacement), float))
    symmetric = np.array([s == "A" for s in modes.symmetry])
    primary = float(values[symmetric].max())

    alternative = gating_overlap_by_distances(np.vstack(cb), np.vstack(fb),
                                              modes)
    assert primary == pytest.approx(0.705, abs=0.02)
    assert abs(alternative - primary) / primary < 0.20, (
        f"{primary:.3f} vs {alternative:.3f}")
    assert alternative > 0.5


# --------------------------------------------------------------------------
# T50
# --------------------------------------------------------------------------

def test_ode_integration_reproduces_the_matrix_exponential():
    """Same quantity, different numerics: no expm, no fixed time grid."""
    model = GatingModel()
    assert t50_by_ode_integration(model) == pytest.approx(
        model.half_activation(), rel=0.05)


def test_the_steady_state_is_a_different_quantity():
    """Recorded because it looked like a disagreement and was not.

    At equilibrium this channel sits ~96% inactivated at every tension, so
    steady-state open occupancy runs only 0.030 to 0.036 and has no
    half-maximum. T₅₀ is necessarily a property of the peak transient — which
    is also what a patch-clamp measures.
    """
    model = GatingModel()
    generator = np.asarray(model.rate_matrix(0.0), float)
    size = generator.shape[0]
    system = np.vstack([generator.T[:-1], np.ones(size)])
    target = np.zeros(size)
    target[-1] = 1.0
    resting = np.linalg.lstsq(system, target, rcond=None)[0]
    assert resting[2] + resting[3] > 0.75, "inactivated states dominate"

    result = t50_by_steady_state(model)
    assert not (2.0 < result < 4.0), (
        "the steady-state route should NOT reproduce the peak-based T50; if it "
        "does, one of the two is not computing what it claims")


# --------------------------------------------------------------------------
# The record
# --------------------------------------------------------------------------

def test_crosscheck_reports_disagreement_as_such():
    agreeing = CrossCheck("q", 1.0, 1.02, tolerance=0.05)
    disagreeing = CrossCheck("q", 1.0, 1.5, tolerance=0.05)
    assert agreeing.agrees and "agree" in agreeing.summary()
    assert not disagreeing.agrees and "DISAGREE" in disagreeing.summary()
    assert disagreeing.relative == pytest.approx(0.5)
