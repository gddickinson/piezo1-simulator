"""Variant impact prediction.

Everything here tests **mechanics**, not whether the predictor separates
gain-of-function from loss-of-function. That comparison belongs to Round 7 and
is deliberately kept out of this file: a predictor tuned until its tests pass
against the phenotype labels is no longer blind, and the blind test is the
point.
"""

import numpy as np
import pytest
from scipy.spatial import cKDTree

from piezo1.analysis.variant_impact import (RESIDUE_VOLUME, VariantImpactModel,
                                            spring_scale_from_volume)
from piezo1.physics.anm import SPRING_MODELS


def _lattice(side: int = 5, spacing: float = 6.0, seed: int = 0):
    rng = np.random.default_rng(seed)
    grid = np.stack(np.meshgrid(*[np.arange(side)] * 3, indexing="ij"), -1)
    coords = grid.reshape(-1, 3).astype(float) * spacing
    return coords + rng.normal(scale=0.3, size=coords.shape)


def _explicit_hessian(coords, site, scale, cutoff=15.0):
    """Build the full Hessian with one residue's springs scaled."""
    n = len(coords)
    h = np.zeros((3 * n, 3 * n))
    pairs = np.asarray(list(cKDTree(coords).query_pairs(cutoff)))
    diff = coords[pairs[:, 1]] - coords[pairs[:, 0]]
    dist = np.linalg.norm(diff, axis=1)
    unit = diff / dist[:, None]
    k = SPRING_MODELS["inverse_square"](dist, 1.0, d0=7.5)
    k = k * np.where((pairs[:, 0] == site) | (pairs[:, 1] == site), scale, 1.0)
    for (i, j), kk, u in zip(pairs, k, unit):
        blk = -kk * np.outer(u, u)
        h[3 * i:3 * i + 3, 3 * j:3 * j + 3] += blk
        h[3 * j:3 * j + 3, 3 * i:3 * i + 3] += blk.T
        h[3 * i:3 * i + 3, 3 * i:3 * i + 3] -= blk
        h[3 * j:3 * j + 3, 3 * j:3 * j + 3] -= blk.T
    return h


@pytest.fixture(scope="module")
def toy():
    coords = _lattice()
    rng = np.random.default_rng(1)
    d = rng.normal(size=coords.shape)
    d /= np.linalg.norm(d)
    model = VariantImpactModel(coords=coords, residues=np.arange(len(coords)),
                               gating_vector=d)
    return model, coords, d


# --------------------------------------------------------------------------
# The core identity
# --------------------------------------------------------------------------

def test_quadratic_form_equals_the_explicit_hessian_difference(toy):
    """½dᵀ(H_mut − H_wt)d computed cheaply must equal the full rebuild.

    The whole method rests on this identity. If it drifts, every predicted ΔΔG
    is wrong in a way no downstream test would notice.
    """
    model, coords, d = toy
    dv = d.ravel()
    for site, scale in ((40, 1.35), (7, 0.6), (99, 2.0)):
        exact = 0.5 * dv @ (_explicit_hessian(coords, site, scale)
                            - _explicit_hessian(coords, site, 1.0)) @ dv
        fast = model.quadratic_form_at(site, scale)
        assert fast == pytest.approx(exact, rel=1e-10, abs=1e-14)


def test_unit_scale_costs_nothing(toy):
    model, _, _ = toy
    assert model.quadratic_form_at(40, 1.0) == pytest.approx(0.0, abs=1e-15)


def test_energy_change_is_linear_in_the_scale(toy):
    """ΔΔG is exactly linear in (s − 1), by construction."""
    model, _, _ = toy
    a = model.quadratic_form_at(40, 1.5)
    b = model.quadratic_form_at(40, 2.0)
    assert b / a == pytest.approx(1.0 / 0.5, rel=1e-12)


# --------------------------------------------------------------------------
# The perturbation model
# --------------------------------------------------------------------------

def test_volume_table_is_complete():
    assert len(RESIDUE_VOLUME) == 20
    assert RESIDUE_VOLUME["G"] < RESIDUE_VOLUME["A"] < RESIDUE_VOLUME["W"]


def test_spring_scale_direction():
    assert spring_scale_from_volume("G", "W") > 1.0     # bigger packs harder
    assert spring_scale_from_volume("W", "G") < 1.0     # smaller leaves a cavity
    assert spring_scale_from_volume("L", "L") == pytest.approx(1.0)
    assert spring_scale_from_volume("A", "X") == 1.0    # unknown -> no change


def test_spring_scale_never_goes_negative():
    for wt in RESIDUE_VOLUME:
        for mut in RESIDUE_VOLUME:
            assert spring_scale_from_volume(wt, mut, sensitivity=5.0) > 0.0


def test_prediction_sign_convention(toy):
    """Stiffening is positive (harder to gate), softening negative."""
    model, _, _ = toy
    bigger = model.predict(40, "G", "W")
    smaller = model.predict(40, "W", "G")
    assert bigger.gating_cost_change > 0 and bigger.sign == "stiffening"
    assert smaller.gating_cost_change < 0 and smaller.sign == "softening"
    same = model.predict(40, "A", "A")
    assert same.gating_cost_change == pytest.approx(0.0, abs=1e-15)
    assert same.sign == "neutral"


# --------------------------------------------------------------------------
# Bookkeeping
# --------------------------------------------------------------------------

def test_all_protomer_copies_are_mutated():
    """A homotrimer carries three copies; scoring one would underestimate it.

    Worse, it would break the C3 symmetry the gating coordinate depends on.
    """
    coords = np.vstack([_lattice(seed=s) + np.array([200.0 * s, 0, 0])
                        for s in range(3)])
    per = len(coords) // 3
    residues = np.tile(np.arange(per), 3)
    rng = np.random.default_rng(2)
    d = rng.normal(size=coords.shape)
    d /= np.linalg.norm(d)
    model = VariantImpactModel(coords=coords, residues=residues, gating_vector=d)

    assert len(model.sites_for(10)) == 3
    trimer = model.predict(10, "A", "W")
    single = model.quadratic_form_at(10, trimer.spring_scale)
    assert trimer.gating_cost_change == pytest.approx(
        single + model.quadratic_form_at(10 + per, trimer.spring_scale)
        + model.quadratic_form_at(10 + 2 * per, trimer.spring_scale), rel=1e-12)


def test_unresolved_residue_is_reported_not_guessed(toy):
    model, _, _ = toy
    p = model.predict(99999, "A", "V")
    assert p.modelled is False
    assert p.gating_cost_change == 0.0
    assert "not resolved" in p.note


def test_gating_vector_shape_is_checked():
    coords = _lattice(side=3)
    with pytest.raises(ValueError, match="gating vector"):
        VariantImpactModel(coords=coords, residues=np.arange(len(coords)),
                           gating_vector=np.zeros((5, 3)))


def test_gating_vector_is_normalised(toy):
    model, _, _ = toy
    assert np.linalg.norm(model.gating_vector) == pytest.approx(1.0)


def test_normalised_score_separates_mild_from_rigid(toy):
    """A residue that barely moves cannot score high however drastic the swap.

    cost_change_normalised divides that out, so "mechanically mild mutation" and
    "position in a rigid region" are distinguishable rather than conflated.
    """
    model, coords, d = toy
    amplitudes = np.linalg.norm(d, axis=1)
    mobile = int(np.argmax(amplitudes))
    rigid = int(np.argmin(amplitudes))
    a = model.predict(mobile, "G", "W")
    b = model.predict(rigid, "G", "W")
    assert abs(a.gating_cost_change) > abs(b.gating_cost_change)
    assert a.local_strain > b.local_strain


def test_predict_all_handles_the_real_variant_table():
    """Coverage must be reported honestly, not silently dropped."""
    from piezo1.analysis.ensemble import build_ensemble
    from piezo1.core.annotations import load_annotations
    try:
        ens = build_ensemble(species="mouse", min_common=900)
    except ValueError as exc:
        pytest.skip(f"ensemble unavailable: {exc}")

    pca = ens.pca()
    model = VariantImpactModel(coords=ens.members[0].coords,
                               residues=np.tile(ens.residues, 3),
                               gating_vector=pca.components[0])
    ann = load_annotations("human")
    preds = model.predict_all(ann.variants, annotations=ann)

    assert len(preds) >= 60
    scored = [p for p in preds if p.modelled and p.wt_aa and p.mut_aa]
    unresolved = [p for p in preds if not p.modelled]
    assert len(scored) >= 30
    assert unresolved, "some variants sit outside the resolved range; say so"
    for p in unresolved:
        assert p.note
    # Domains should be attached where the residue is known.
    assert any(p.domain for p in scored)
