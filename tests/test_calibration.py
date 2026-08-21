"""Every checking instrument, calibrated against a case with a known answer.

This project's most expensive mistakes have all been the same shape: an
*alternative* route, built to check the main one, was itself wrong and produced
a plausible number rather than an error. A spheroid fitter that would have
reported 89% model error. A document checker that could not read its own
documents' minus sign. A probe whose "no effect" came from badly chosen
coordinates. In each case the checker disagreed with the pipeline and the
disagreement looked like a finding.

So a checking instrument is a measuring instrument, and an uncalibrated one is
worse than none — it manufactures findings. This file holds the register of
instruments and the calibration each one is subject to, and a guard that no
instrument can be added without one.

The register was assembled by auditing the modules rather than by trusting
their docstrings. Four instruments turned out to have no known-answer case and
are calibrated here for the first time.
"""

from __future__ import annotations

import importlib
import itertools
import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent

#: Modules whose public callables exist to *check* something. Every lowercase
#: name in their ``__all__`` must appear in :data:`CALIBRATED` below.
from calibration_register import CALIBRATED  # noqa: E402

CHECKING_MODULES = (
    "analysis.crosscheck",
    "analysis.hydropathy",
    "structure.planarity",
    "structure.architecture",
    "structure.micelle",
    "analysis.crosscheck_methods",
    "analysis.model_error",
    "analysis.uncertainty",
    "analysis.validation",
    "analysis.design",
    "analysis.provenance_chain",
    "analysis.parameter_effect",
    "analysis.fluctuations",
    # Round 89. Both are instruments in the strict sense: neither measures a
    # property of PIEZO1, both exist to decide whether another number may be
    # believed. ``alignment_windows`` earns its place twice over — it was
    # wrong twice before it was calibrated, once in its statistic and once in
    # its width.
    "analysis.homology",
    "analysis.alignment_windows",
    "analysis.homology_structure",
    "structure.assembly",
    "structure.clashes",
    # Round 93. The family subsystem joins an *external project's* results to
    # this one's coordinates, and every module in it exists to decide whether
    # something may be believed rather than to measure a property of PIEZO1.
    # Two of them earn the register outright: ``constraint_mechanics`` reports
    # a correlation, which is the easiest thing in science to manufacture, and
    # ``core_periphery`` reports a ratio whose denominator can be anything.
    "analysis.family_constraint",
    "analysis.constraint_mechanics",
    "analysis.core_periphery",
    "analysis.equivalent_positions",
    "analysis.disease_geography",
    "analysis.family_motifs",
    "analysis.piezo3",
    "analysis.piezo3_channel",
    "dead_code",
)

#: instrument -> the test that drives it against a known answer.
# ---------------------------------------------------- the register guard

def test_every_checking_instrument_has_a_calibration():
    """No instrument may be added without a known-answer case.

    This is the mechanism that keeps the rule from decaying into an
    aspiration, in the same way `parameter_audit` does for registered numbers.
    """
    missing = []
    for module in CHECKING_MODULES:
        mod = importlib.import_module(f"piezo1.{module}")
        short = module.split(".")[-1]
        for name in getattr(mod, "__all__", []):
            if name[:1].isupper() or name.startswith("_"):
                continue                      # dataclasses carry no answer
            if f"{short}.{name}" not in CALIBRATED:
                missing.append(f"{short}.{name}")
    assert not missing, (
        "these checking instruments have no calibration registered: "
        + ", ".join(sorted(missing)))


def test_named_calibrating_tests_exist():
    """A register pointing at a test that does not exist proves nothing."""
    missing = []
    for instrument, where in CALIBRATED.items():
        filename = where.split("::")[0].split(" ")[0]
        path = ROOT / filename
        if not path.exists():
            missing.append(f"{instrument} -> {filename} (no such file)")
            continue
        if "::" in where:
            testname = where.split("::")[1].split(" ")[0]
            if not re.search(rf"def {re.escape(testname)}\b", path.read_text()):
                missing.append(f"{instrument} -> {testname} (not in {filename})")
    assert not missing, "\n".join(missing)


# ------------------------------------- the four calibrations added here

def test_spring_model_error_recovers_a_planted_mode():
    """Plant a mode from one spring model; that model must score ~1.

    Previously this instrument was only checked by asserting its spread on real
    structures was "modest" — which is a statement about the answer, not about
    the instrument. If it had been computing overlaps wrongly, a modest spread
    would have looked exactly the same.
    """
    from piezo1.analysis.model_error import spring_model_error
    from piezo1.physics.anm import ANM

    rng = np.random.default_rng(3)
    blocks = [rng.normal(scale=9.0, size=(60, 3)) for _ in range(3)]
    anm = ANM.from_trimer(blocks, cutoff=15.0, spring="inverse_square").build()
    planted = anm.calc_modes(n_modes=12).vectors[0]

    error = spring_model_error(blocks, planted, n_modes=12, cutoff=15.0)
    recovered = error.values["inverse_square"]
    others = [v for k, v in error.values.items() if k != "inverse_square"]

    # The calibration is the SEPARATION, not an absolute threshold. On these
    # random coordinates the low modes are near-degenerate, and ARPACK starts
    # from an unseeded random vector, so the recovered overlap varies run to
    # run between about 0.95 and 1.00 with identical inputs. That is a property
    # of degenerate geometry, not of the instrument: the real structures have
    # well-separated low modes and `anm.gating_overlap` reproduces to every
    # digit. Asserting a tight absolute bound here would make this test flaky
    # for a reason that has nothing to do with what it is checking.
    assert recovered > 0.9, (
        f"the model that generated the displacement must recover it, got "
        f"{recovered:.4f}")
    assert recovered > 1.5 * max(others), (
        f"the generating model must stand clear of the others: "
        f"{recovered:.4f} vs {max(others):.4f}")
    assert error.spread > 0.3, "the instrument must be able to say they differ"


def test_minimum_detectable_effect_delivers_its_target_power():
    """The effect it returns must actually be detectable at the stated power.

    Until now this was only checked against the project's own recorded results,
    which is circular: the same function supplies the power statements those
    rounds are judged by. This closes the loop against the simulator instead.
    """
    from piezo1.analysis.design import minimum_detectable_effect, power_curve

    for n_a, n_b in ((20, 20), (16, 9)):
        delta = minimum_detectable_effect(n_a, n_b, target_power=0.8)
        achieved = power_curve(n_a, n_b, deltas=[-delta], n_simulations=3000,
                               n_permutations=999, seed=7).power[0]
        assert achieved == pytest.approx(0.80, abs=0.06), (
            f"MDE {delta:.3f} at n={n_a}/{n_b} gives power {achieved:.3f}")


def test_sensitivity_range_is_exactly_the_known_extremes():
    """Over known outputs the range must be their min and max, exactly.

    The existing tests check that it refuses to call itself a confidence
    interval — which is about the wording, not the arithmetic.
    """
    from piezo1.analysis.uncertainty import sensitivity

    values = [2.0, 5.0, 3.0]
    result = sensitivity(lambda i: values[i], settings=[0, 1, 2],
                         reference=1, knob="i", what="calibration")
    assert result.low == pytest.approx(2.0)
    assert result.high == pytest.approx(5.0)
    assert result.estimate == pytest.approx(5.0), "estimate is the reference"

    # A single setting has no spread — the degenerate case must not widen it.
    flat = sensitivity(lambda i: 4.0, settings=[0], reference=0)
    assert flat.low == flat.high == pytest.approx(4.0)


def test_permutation_p_matches_exhaustive_enumeration():
    """For samples small enough to enumerate, the sampler must agree.

    Eight values give C(8,4) = 70 partitions, so the exact one-sided p is
    computable directly. The sampled value is very slightly larger because of
    the (r+1)/(n+1) convention, which is deliberately conservative — that
    direction is asserted rather than merely tolerated.
    """
    from piezo1.analysis.validation import permutation_test

    a = np.array([1.0, 2, 3, 4])
    b = np.array([5.0, 6, 7, 8])
    pool = np.concatenate([a, b])
    observed = a.mean() - b.mean()

    hits = total = 0
    for index in itertools.combinations(range(8), 4):
        mask = np.zeros(8, bool)
        mask[list(index)] = True
        if pool[mask].mean() - pool[~mask].mean() <= observed:
            hits += 1
        total += 1
    exact = hits / total

    sampled = permutation_test(a, b, n_permutations=20000, seed=1).p_value
    assert sampled == pytest.approx(exact, abs=0.005)
    assert sampled >= exact, "the (r+1)/(n+1) convention must not anti-conserve"


def test_the_calibrations_added_here_would_fail_on_a_broken_instrument():
    """Each new calibration must be able to reject, or it asserts nothing."""
    from piezo1.analysis.uncertainty import sensitivity

    # A 'sensitivity' that returned the reference for both ends would pass any
    # wording test and fail this one.
    values = [2.0, 5.0, 3.0]
    result = sensitivity(lambda i: values[i], settings=[0, 1, 2], reference=1)
    assert result.low != result.high

    from piezo1.analysis.validation import permutation_test
    identical = np.arange(4, dtype=float)
    assert permutation_test(identical, identical.copy()).p_value > 0.3, (
        "identical samples must not look significant")
