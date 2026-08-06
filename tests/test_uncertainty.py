"""Intervals, and keeping the three kinds of spread distinct.

Rounds 18–28 repeatedly found a recorded number stated with more confidence than
it had earned. This module attaches a spread — but a bootstrap over data and a
sweep over a method parameter mean different things, and conflating them would
be a second kind of overconfidence rather than a cure for the first.
"""

import numpy as np
import pytest

from piezo1.analysis.uncertainty import (Bootstrap, ParameterRange,
                                         Sensitivity, bootstrap,
                                         format_with_interval, parameter_range,
                                         sensitivity)
from piezo1.parameters import PARAMETERS


# --------------------------------------------------------------------------
# The bootstrap
# --------------------------------------------------------------------------

def test_bootstrap_recovers_a_known_interval():
    """For a mean of standard normals the 95% interval is about ±1.96/√n."""
    rng = np.random.default_rng(0)
    data = rng.normal(size=500)
    result = bootstrap(lambda idx: float(data[idx].mean()), data,
                       n_resamples=800, seed=1)
    expected_half_width = 1.96 / np.sqrt(len(data))
    assert result.width == pytest.approx(2 * expected_half_width, rel=0.25)
    assert result.contains(0.0)


def test_a_tighter_dataset_gives_a_tighter_interval():
    rng = np.random.default_rng(2)
    wide = rng.normal(scale=5.0, size=300)
    tight = rng.normal(scale=0.5, size=300)
    a = bootstrap(lambda i: float(wide[i].mean()), wide, n_resamples=300)
    b = bootstrap(lambda i: float(tight[i].mean()), tight, n_resamples=300)
    assert b.width < a.width


def test_bootstrap_is_deterministic_for_a_given_seed():
    data = np.arange(100.0)
    a = bootstrap(lambda i: float(data[i].mean()), data, n_resamples=200, seed=7)
    b = bootstrap(lambda i: float(data[i].mean()), data, n_resamples=200, seed=7)
    assert (a.low, a.high) == (b.low, b.high)


def test_a_fragile_statistic_raises_rather_than_narrowing():
    """A resample that fails must not be replaced by the point estimate.

    Silently substituting would pull the distribution toward the centre and
    report an interval far tighter than the data support. Here the statistic
    works on the full data — so there is a point estimate — but fails on almost
    every resample.
    """
    data = np.arange(50.0)

    def fragile(index):
        if len(np.unique(index)) < len(index):     # any resample has repeats
            raise ValueError("cannot handle duplicates")
        return float(data[index].mean())

    with pytest.raises(RuntimeError, match="too fragile"):
        bootstrap(fragile, data, n_resamples=100)


def test_a_statistic_that_fails_on_the_full_data_raises_its_own_error():
    """Not the bootstrap's problem to disguise."""
    with pytest.raises(ValueError, match="no"):
        bootstrap(lambda i: (_ for _ in ()).throw(ValueError("no")),
                  np.arange(50.0), n_resamples=10)


def test_too_few_observations_is_refused():
    with pytest.raises(ValueError, match="at least three"):
        bootstrap(lambda i: 1.0, np.zeros(2))


def test_interval_level_is_honoured():
    rng = np.random.default_rng(3)
    data = rng.normal(size=400)
    narrow = bootstrap(lambda i: float(data[i].mean()), data,
                       n_resamples=400, level=0.50, seed=4)
    wide = bootstrap(lambda i: float(data[i].mean()), data,
                     n_resamples=400, level=0.99, seed=4)
    assert narrow.width < wide.width


# --------------------------------------------------------------------------
# Sensitivity is not a confidence interval
# --------------------------------------------------------------------------

def test_the_three_kinds_are_named_differently():
    """The distinction is the point of the module, so it is asserted."""
    assert Bootstrap.kind == "confidence interval"
    assert Sensitivity.kind == "sensitivity range"
    assert ParameterRange.kind == "parameter range"
    assert Bootstrap.kind != Sensitivity.kind != ParameterRange.kind


def test_sensitivity_summary_refuses_the_words_confidence_interval():
    result = sensitivity(lambda s: float(s) ** 2, [1.0, 2.0, 3.0],
                         reference=2.0, knob="setting")
    text = result.summary()
    assert "sensitivity" in text
    assert "confidence interval" in text and "not a confidence interval" in text


def test_sensitivity_reports_the_reference_as_the_estimate():
    result = sensitivity(lambda s: float(s), [1.0, 5.0, 9.0], reference=5.0,
                         knob="k")
    assert result.estimate == 5.0
    assert (result.low, result.high) == (1.0, 9.0)
    assert result.spread == 8.0


def test_settings_that_fail_are_dropped_not_counted():
    def fussy(setting):
        if setting == 2.0:
            raise ValueError("bad setting")
        return float(setting)

    result = sensitivity(fussy, [1.0, 2.0, 3.0], reference=1.0, knob="k")
    assert list(result.settings) == [1.0, 3.0]
    assert len(result.values) == 2


def test_all_settings_failing_raises():
    with pytest.raises(RuntimeError, match="no setting"):
        sensitivity(lambda s: 1 / 0, [1.0, 2.0], knob="k")


# --------------------------------------------------------------------------
# Parameter propagation restores the registry
# --------------------------------------------------------------------------

def test_parameter_range_restores_the_registry():
    """Leaving the registry modified would make every later number in the
    session incomparable with the documentation."""
    assert PARAMETERS.is_default("membrane.kappa")
    parameter_range(lambda k: PARAMETERS.value("membrane.kappa"),
                    "membrane.kappa", [20.0, 25.0])
    assert PARAMETERS.is_default("membrane.kappa")
    assert PARAMETERS.value("membrane.kappa") == 20.0


def test_parameter_range_restores_even_when_the_statistic_raises():
    """The exception propagates — that is correct — but the registry must not
    be left modified, or every later number in the session becomes
    incomparable with the documentation."""
    def explodes(_value):
        raise RuntimeError("boom")

    with pytest.raises(RuntimeError, match="boom"):
        parameter_range(explodes, "membrane.kappa", [25.0])
    assert PARAMETERS.is_default("membrane.kappa")
    assert PARAMETERS.value("membrane.kappa") == 20.0


def test_parameter_range_actually_varies_the_parameter():
    result = parameter_range(lambda k: PARAMETERS.value("membrane.kappa"),
                             "membrane.kappa", [20.0, 22.5, 25.0])
    assert sorted(result.values) == [20.0, 22.5, 25.0]
    assert result.estimate == 20.0


# --------------------------------------------------------------------------
# On the real numbers
# --------------------------------------------------------------------------

def test_the_dome_radius_interval_contains_the_published_value(
        curved_structure):
    """The reframing this round produced.

    The project reported 9.7 nm against a published 10.2 nm as a near-miss. The
    interval says the two are statistically indistinguishable given 66 surface
    points — a stronger claim of consistency and a weaker claim of precision.
    """
    from piezo1.structure.geometry import fit_sphere

    from test_geometry import _tm_surface

    surface = _tm_surface(curved_structure, "mouse")
    assert len(surface) > 40
    result = bootstrap(
        lambda idx: fit_sphere(surface[idx], iterations=4, trim=0.15).radius / 10.0,
        surface, n_resamples=200, seed=0)
    assert result.estimate == pytest.approx(9.72, abs=0.1)
    assert result.contains(PARAMETERS.value("dome.published_radius_closed"))
    assert result.relative_width > 0.05, (
        "66 points cannot determine a radius to better than a few percent")


def test_the_gating_overlap_is_cutoff_dependent(curved_structure,
                                                flat_structure):
    """0.705 is not robust to three digits.

    The qualitative result — a substantial overlap carried entirely by
    A-symmetric modes — survives every cutoff. The specific value does not, and
    reporting it without that caveat was overconfident.
    """
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
    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    displacement = ((((np.vstack(fb) - centroid) @ rotation.T + translation)
                     - np.vstack(cb))).ravel()

    def overlap_at(cutoff):
        anm = ANM.from_trimer(cb, cutoff=cutoff, spring="inverse_square").build()
        modes = anm.calc_modes(n_modes=30)
        anm.label_symmetry(modes)
        values = np.abs(np.asarray(modes.overlap(displacement), float))
        symmetric = np.array([s == "A" for s in modes.symmetry])
        return float(values[symmetric].max())

    result = sensitivity(overlap_at, [12.0, 15.0, 18.0], reference=15.0,
                         knob="anm.cutoff")
    assert result.estimate == pytest.approx(0.705, abs=0.02)
    assert result.relative_spread > 0.10, "the spread is real and worth stating"
    assert result.low > 0.4, "the qualitative conclusion survives every cutoff"


def test_format_names_the_kind_of_spread():
    result = sensitivity(lambda s: float(s), [1.0, 2.0], reference=1.0, knob="k")
    assert "not a confidence interval" in format_with_interval(1.0, result)
    assert "no interval" in format_with_interval(1.0, None)
