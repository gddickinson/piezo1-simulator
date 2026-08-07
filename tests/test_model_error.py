"""Model-form error: the term the confidence intervals do not contain.

The fitter is tested against spheroids whose answer is known before it is
allowed anywhere near the dome, because the first version of it was wrong in a
way that would have been reported as science: it inflated both semi-axes by the
same 1.63x factor on a *known* spheroid, and the resulting "model error" would
have been a fitting bug wearing a result's clothes.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.model_error import (ModelError, compare_with_sampling,
                                         fit_spheroid)


# ------------------------------------------------------- the fitter, on knowns

@pytest.mark.parametrize("equatorial,polar", [(100.0, 60.0), (80.0, 80.0),
                                              (120.0, 40.0), (50.0, 45.0)])
def test_spheroid_fitter_recovers_a_known_shape(equatorial, polar):
    """Including the sphere as a special case, which must come back a sphere."""
    rng = np.random.default_rng(0)
    u = rng.uniform(0, 2 * np.pi, 4000)
    v = np.arccos(rng.uniform(-1, 1, 4000))
    points = np.column_stack([equatorial * np.sin(v) * np.cos(u),
                              equatorial * np.sin(v) * np.sin(u),
                              polar * np.cos(v)])
    offset = np.array([5.0, -3.0, 2.0])
    points = points + offset + rng.normal(scale=0.3, size=points.shape)

    fit = fit_spheroid(points)
    assert fit.equatorial == pytest.approx(equatorial, rel=0.01)
    assert fit.polar == pytest.approx(polar, rel=0.01)
    assert np.linalg.norm(fit.center - offset) < 0.5
    assert fit.rmse < 1.0


def test_the_broken_fitter_would_have_been_caught():
    """The specific failure the first implementation had.

    It returned a = 163 for a true 100 and c = 98 for a true 60 — both inflated
    by the same factor, the signature of a drifting centre. This asserts the
    ratio is right *and* the scale is, since a shape-only check would have
    passed the broken version.
    """
    rng = np.random.default_rng(1)
    u = rng.uniform(0, 2 * np.pi, 3000)
    v = np.arccos(rng.uniform(-1, 1, 3000))
    points = np.column_stack([100 * np.sin(v) * np.cos(u),
                              100 * np.sin(v) * np.sin(u), 60 * np.cos(v)])
    fit = fit_spheroid(points)
    assert fit.equatorial / fit.polar == pytest.approx(100 / 60, rel=0.01)
    assert fit.equatorial == pytest.approx(100.0, rel=0.01)


def test_apex_curvature_is_not_a_semi_axis():
    """For a spheroid the polar radius of curvature is a^2/c, not c.

    Quoting c would understate the curvature by the flattening, which is the
    obvious mistake and is silent.
    """
    rng = np.random.default_rng(2)
    u = rng.uniform(0, 2 * np.pi, 3000)
    v = np.arccos(rng.uniform(-1, 1, 3000))
    points = np.column_stack([100 * np.sin(v) * np.cos(u),
                              100 * np.sin(v) * np.sin(u), 50 * np.cos(v)])
    fit = fit_spheroid(points)
    assert fit.apex_curvature == pytest.approx(100.0 ** 2 / 50.0, rel=0.02)
    assert fit.apex_curvature > fit.equatorial > fit.polar
    assert fit.flattening == pytest.approx(0.5, abs=0.02)


def test_a_sphere_has_zero_flattening():
    rng = np.random.default_rng(3)
    points = rng.normal(size=(3000, 3))
    points = 70.0 * points / np.linalg.norm(points, axis=1, keepdims=True)
    fit = fit_spheroid(points)
    assert abs(fit.flattening) < 0.01
    assert fit.apex_curvature == pytest.approx(70.0, rel=0.02)


# ------------------------------------------------------------- the container

def test_model_error_reports_a_lower_bound_not_an_interval():
    """The claim has to be stated on the object, not left to be remembered."""
    error = ModelError(quantity="q", values={"a": 1.0, "b": 3.0},
                       reference="a")
    assert error.estimate == 1.0          # the reference model, not the median
    assert error.spread == 2.0
    assert "LOWER BOUND" in error.summary()
    assert error.kind == "model error"


def test_compare_names_which_error_dominates():
    from piezo1.analysis.uncertainty import Bootstrap

    model = ModelError(quantity="q", values={"a": 9.45, "b": 14.99},
                       reference="a")
    sampling = Bootstrap(estimate=9.45, low=9.05, high=9.97, n_resamples=400)
    verdict = compare_with_sampling(model, sampling)
    assert verdict["dominant"] == "model"
    assert verdict["ratio"] > 5
    assert "not whether it is the right model" in verdict["verdict"]

    # And the other way round, so the function is not just always saying "model".
    tight = ModelError(quantity="q", values={"a": 9.45, "b": 9.46},
                       reference="a")
    assert compare_with_sampling(tight, sampling)["dominant"] == "sampling"


# ------------------------------------------------------ against the structures

def test_model_error_dominates_the_dome_radius(curved_structure):
    """The round's headline, measured.

    A sphere and an oblate spheroid both fit the transmembrane surface — the
    spheroid slightly better, as it must with an extra parameter — and give
    radii of curvature differing by 59%. The bootstrap interval on the sphere is
    six times narrower than that, so it measures how well a sphere is determined
    rather than whether a sphere was the right shape.
    """
    from piezo1.analysis.claims import _tm_surface
    from piezo1.analysis.model_error import dome_model_error
    from piezo1.analysis.uncertainty import bootstrap
    from piezo1.structure.geometry import fit_sphere
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    blocks, _ = protomer_blocks(curved_structure)
    axis = detect_c3_axis(blocks)
    surface = _tm_surface(curved_structure, "mouse")

    sphere = fit_sphere(surface)
    spheroid_error = dome_model_error(surface, axis, sphere.radius / 10.0)
    assert spheroid_error.values["sphere"] == pytest.approx(9.45, abs=0.2)
    assert spheroid_error.values["spheroid (apex)"] == pytest.approx(15.0, abs=1.5)
    assert spheroid_error.relative_spread > 0.4

    sampling = bootstrap(lambda idx: fit_sphere(surface[idx]).radius / 10.0,
                         surface, n_resamples=200)
    verdict = compare_with_sampling(spheroid_error, sampling)
    assert verdict["dominant"] == "model"
    assert verdict["ratio"] > 3.0


def test_the_spheroid_is_a_competitive_fit_not_a_failed_one(curved_structure):
    """Otherwise the "model error" would just be a bad fit disagreeing.

    A model that fits far worse is rejected, not an alternative. The spheroid
    has one more parameter than the sphere so it must fit at least as well; if
    it ever fits *worse*, the fitter has regressed.
    """
    from piezo1.analysis.claims import _tm_surface
    from piezo1.structure.geometry import fit_sphere
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    blocks, _ = protomer_blocks(curved_structure)
    axis = detect_c3_axis(blocks)
    surface = _tm_surface(curved_structure, "mouse")

    sphere = fit_sphere(surface)
    sphere_residual = np.sqrt(np.mean(
        (np.linalg.norm(surface - sphere.center, axis=1) - sphere.radius) ** 2))
    spheroid = fit_spheroid(surface, axis.direction)
    assert spheroid.rmse <= sphere_residual + 1e-6, (
        f"spheroid rmse {spheroid.rmse:.3f} A is worse than the sphere's "
        f"{sphere_residual:.3f} A; the fitter has regressed")


def test_spring_model_spread_is_modest(curved_structure, flat_structure):
    """The elastic network is less model-sensitive than the dome geometry."""
    from piezo1.analysis.model_error import spring_model_error
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import kabsch, match_protomers

    curved_blocks, curved_res = protomer_blocks(curved_structure)
    flat_blocks, flat_res = protomer_blocks(flat_structure)
    common = np.array(sorted(set(curved_res.tolist()) & set(flat_res.tolist())))
    a = [b[np.isin(curved_res, common)] for b in curved_blocks]
    b = [b[np.isin(flat_res, common)] for b in flat_blocks]
    b = [b[i] for i in match_protomers(a, b).order]
    rotation, translation, centroid = kabsch(np.vstack(b), np.vstack(a))
    displacement = (((np.vstack(b) - centroid) @ rotation.T + translation)
                    - np.vstack(a))

    error = spring_model_error(a, displacement, n_modes=20)
    assert set(error.values) == {"uniform", "inverse_square", "inverse_sixth"}
    assert error.relative_spread < 0.15
    # Every spring model must still find the transition; none is degenerate.
    assert min(error.values.values()) > 0.8


def test_pore_conventions_coincide_at_the_carbon_radius(curved_structure):
    """A null with a mechanism, not a null from a broken check.

    7WLT's bottleneck is carbon-lined, so a uniform probe at carbon's 1.7 A
    reproduces the Apollonius answer exactly. Moving the probe off 1.7 shifts
    the answer by exactly the offset, which proves the check is live rather than
    silently returning the same number.
    """
    from piezo1.analysis.model_error import pore_convention_error
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    axis = detect_c3_axis(protomer_blocks(curved_structure)[0])
    same = pore_convention_error(curved_structure, axis, uniform_radius=1.7)
    assert same.spread == pytest.approx(0.0, abs=1e-6)

    for radius in (1.4, 2.0):
        shifted = pore_convention_error(curved_structure, axis,
                                        uniform_radius=radius)
        assert shifted.spread == pytest.approx(abs(1.7 - radius), abs=0.02)
