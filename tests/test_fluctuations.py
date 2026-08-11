"""The elastic network against the deposited B-factors, and its calibration.

The comparison is a checking instrument — it is the standard test of whether an
elastic network describes the molecule it was built from — so nothing about
PIEZO1 is believed here until three things with known answers have been
checked. That the two correlation statistics return their analytic limits. That
a network built from the *right* coordinates recovers a planted fluctuation
while one built from shuffled coordinates does not. And that the three ways a
B-factor column can be uninterpretable are each refused rather than averaged
into a number.

The last of those is the one that matters most in practice: an AlphaFold model
carries pLDDT in that field, and comparing against it returns a confident
negative correlation that means nothing at all. The gate is provenance-based,
so the test that it points the right way has to be a measurement — and it is:
the predicted models anti-correlate where the deposited ones do not.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.parameters import PARAMETERS
from piezo1.physics.anm import ANM
from piezo1.analysis.fluctuations import (BFactorQuality, assess_b_factors,
                                          compare_fluctuations,
                                          contact_number, observed_b_factors,
                                          pearson, predicted_msf, spearman,
                                          survey_fluctuations)


def _blob(n_per: int, seed: int = 0):
    """Three identical compact protomer blocks — a trimer with no biology in it."""
    rng = np.random.default_rng(seed)
    core = rng.normal(scale=8.0, size=(n_per, 3))
    blocks = []
    for k in range(3):
        angle = 2.0 * np.pi * k / 3.0
        rotation = np.array([[np.cos(angle), -np.sin(angle), 0.0],
                             [np.sin(angle), np.cos(angle), 0.0],
                             [0.0, 0.0, 1.0]])
        blocks.append((core + np.array([25.0, 0.0, 0.0])) @ rotation.T)
    return blocks


def _require(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


# --------------------------------------------- calibration: the two statistics

def test_the_correlations_return_their_analytic_limits():
    x = np.array([1.0, 2.0, 3.0, 4.0, 5.0])
    assert pearson(x, 2.0 * x + 7.0) == pytest.approx(1.0)
    assert pearson(x, -x) == pytest.approx(-1.0)
    assert spearman(x, np.exp(x)) == pytest.approx(1.0), (
        "a monotone transform must leave the rank correlation at one")
    assert pearson(x, np.exp(x)) < 0.95, (
        "...and must not leave Pearson there, which is why both are reported")
    # A constant has no correlation with anything; nan rather than a number.
    assert np.isnan(pearson(x, np.ones(5)))
    assert np.isnan(pearson(np.array([1.0]), np.array([1.0])))


def test_ranks_average_ties_rather_than_inventing_an_order():
    """A grouped B-factor column is mostly ties, so this decides those entries."""
    from piezo1.analysis.fluctuations import _ranks

    assert list(_ranks(np.array([5.0, 5.0, 1.0]))) == [1.5, 1.5, 0.0]
    # Two columns that are entirely tied cannot correlate with anything.
    assert np.isnan(spearman(np.ones(6), np.arange(6.0)))
    # And a perfect correlation survives ties being present.
    a = np.array([1.0, 1.0, 2.0, 3.0])
    assert spearman(a, 10.0 * a) == pytest.approx(1.0)


# ------------------------------ calibration: good network versus a bad network

def test_the_comparison_recovers_a_planted_fluctuation_and_a_bad_network_does_not():
    """The requirement the round was set: can it tell a good network from a bad one?

    The observation is planted *from* the network built on the real
    coordinates, so the right answer is exactly one. A network built on the
    same points shuffled has the same beads, the same cutoff and the same
    spring model, and differs only in which of them are neighbours — so if the
    comparison could not tell them apart it would not be measuring the network
    at all.
    """
    blocks = _blob(60)
    good = ANM.from_trimer(blocks).build().calc_modes(n_modes=20)
    planted = predicted_msf(good)

    assert pearson(planted, planted) == pytest.approx(1.0)
    assert spearman(planted, planted) == pytest.approx(1.0)

    rng = np.random.default_rng(3)
    scrambled = [b[rng.permutation(len(b))] for b in blocks]
    bad = ANM.from_trimer(scrambled).build().calc_modes(n_modes=20)
    bad_msf = predicted_msf(bad)

    assert abs(spearman(bad_msf, planted)) < 0.6, (
        "a network on shuffled coordinates should not reproduce the planted "
        "fluctuation; if it does, the comparison is measuring the geometry "
        "rather than the network")
    assert spearman(planted, planted) > abs(spearman(bad_msf, planted))


def test_a_shuffled_observation_correlates_with_nothing():
    """The true null. An instrument that cannot say "no" asserts nothing.

    The bound is the null's own standard deviation, ``1/sqrt(n-1)``, rather
    than a round number: with sixty residues a correlation of 0.25 is two
    sigma and turns up about once in twenty draws, so a fixed threshold would
    either be vacuous or flake. Five sigma is a statement about the
    distribution; 0.25 was a statement about one seed.
    """
    blocks = _blob(60, seed=5)
    modes = ANM.from_trimer(blocks).build().calc_modes(n_modes=20)
    planted = predicted_msf(modes)
    sigma = 1.0 / np.sqrt(len(planted) - 1)
    rng = np.random.default_rng(11)
    values = np.array([abs(spearman(planted, rng.permutation(planted)))
                       for _ in range(30)])
    assert values.max() < 5.0 * sigma, f"reached {values.max():.2f}"
    assert np.median(values) < 1.5 * sigma
    # ...and the planted signal is nowhere near that band.
    assert spearman(planted, planted) - values.max() > 0.5


def test_contact_number_is_a_control_and_not_a_copy_of_the_prediction():
    """The control has to be able to disagree, or it is not a control."""
    blocks = _blob(60, seed=7)
    modes = ANM.from_trimer(blocks).build().calc_modes(n_modes=20)
    predicted = predicted_msf(modes)
    control = contact_number(np.vstack(blocks))

    assert len(control) == len(predicted)
    # Related — burial is genuinely part of why a residue is rigid — but not
    # the same quantity, which is what makes beating it meaningful.
    assert 0.1 < spearman(control, predicted) < 0.98


# ----------------------------------------- calibration: the three refusals

def test_a_uniform_column_is_refused_rather_than_correlated():
    import dataclasses

    st = _require("7WLT")
    flat = dataclasses.replace(st, b_factor=np.full_like(st.b_factor, 42.0))
    quality = assess_b_factors(flat)
    assert not quality.usable and "uniform" in quality.reason
    assert not compare_fluctuations(flat).available


def test_a_grouped_column_is_refused_with_the_number_that_refused_it():
    """3JAC and 6BPZ are the real cases; this pins the rule on a made one."""
    import dataclasses

    st = _require("7WLT")
    floor = PARAMETERS.value("fluctuation.min_distinct_fraction")
    # One distinct value per fifty residues: far below the floor.
    grouped = np.round(st.b_factor / 50.0) * 50.0
    coarse = dataclasses.replace(st, b_factor=grouped.astype(np.float32))
    quality = assess_b_factors(coarse)
    assert quality.distinct_fraction < floor
    assert not quality.usable and "grouped" in quality.reason

    # ...and the real column of the same entry passes, so the rule is not
    # simply refusing everything.
    assert assess_b_factors(st).usable


def test_a_predicted_model_is_refused_and_the_gate_points_the_right_way():
    """pLDDT is a confidence, so it must run *against* the fluctuation.

    The refusal is decided by provenance, which is a decision rather than a
    measurement — so the measurement is here instead: build the network on the
    AlphaFold monomer and correlate its own column. A B-factor would agree with
    the fluctuation; pLDDT disagrees with it, which is what the gate exists to
    stop being reported as a result.
    """
    st = _require("AF-Q92508-F1-model_v6")
    quality = assess_b_factors(st)
    assert quality.is_confidence
    assert not quality.usable and "pLDDT" in quality.reason
    assert not compare_fluctuations(st).available

    mask = st.mask_ca()
    modes = ANM(st.xyz[mask].astype(float)).build().calc_modes(n_modes=20)
    plddt = st.b_factor[mask].astype(float)
    assert spearman(modes.msf(), plddt) < -0.2, (
        "pLDDT should anti-correlate with fluctuation; if it does not, the "
        "gate is refusing these entries for the wrong reason")


def test_quality_is_measured_before_any_network_is_built():
    """A refusal must not be able to depend on the answer it would have given."""
    st = _require("3JAC")
    quality = assess_b_factors(st)
    assert not quality.usable and "grouped" in quality.reason
    assert quality.n_distinct == 212 and quality.n_residues == 2754


# ------------------------------------------------------ the measured result

def test_the_survey_of_every_downloaded_entry(structure_by_id):
    """Round 82's numbers, including the ones that are poor.

    Eighteen of twenty-one entries can answer; three cannot and say why. The
    network's median rank correlation is about 0.70 against a contact-number
    control's 0.25, and it beats that control on all but two entries whose
    column behaves like a mobility at all.
    """
    if structure_by_id("7WLT") is None:
        pytest.skip("structures not downloaded — run python -m piezo1.io.fetch")
    rows = survey_fluctuations()
    usable = [r for r in rows if r["usable"]]
    refused = [r for r in rows if not r["usable"]]

    assert len(rows) >= 20
    assert len(usable) >= 17, f"only {len(usable)} entries could answer"
    assert {r["pdb"] for r in refused} == {"3JAC", "6BPZ", "4RAX"}
    for row in refused:
        assert row["reason"], f"{row['pdb']} was refused without a reason"

    spearmans = np.array([r["spearman"] for r in usable])
    controls = np.array([r["control_spearman"] for r in usable])
    assert 0.60 < np.median(spearmans) < 0.80
    assert np.median(controls) < np.median(spearmans)

    # Entries whose control is negative have a column that disagrees with
    # burial, so beating it there is not evidence; they are excluded.
    honest = [r for r in usable if not r["control_inverted"]]
    assert len(honest) >= 14
    wins = sum(r["beats_control"] for r in honest)
    assert wins >= len(honest) - 3, f"network beat the control only {wins} times"

    # The two it loses to are named, so a change in either reopens the question.
    losers = {r["pdb"] for r in honest if not r["beats_control"]}
    assert losers <= {"6KG7", "8IXN"}, f"a new entry now fails: {losers}"

    # The nuance that stops this being read as a clean validation: on the
    # *linear* statistic the network barely beats burial, and on more than a
    # third of entries it loses. If this ever stops being true the headline
    # has changed and the prose has to change with it.
    pearson_wins = sum(r["pearson"] > r["control_pearson"] for r in honest)
    rank_wins = sum(r["spearman"] > r["control_spearman"] for r in honest)
    assert pearson_wins < rank_wins, (
        "the rank correlation is supposed to be where the network wins; if "
        "Pearson has caught up, the reported asymmetry is gone")
    assert pearson_wins <= 0.75 * len(honest)


def test_the_headline_entry_and_its_control(structure_by_id):
    """7WLT, where the two statistics disagree most instructively."""
    st = structure_by_id("7WLT")
    if st is None:
        pytest.skip("7WLT not downloaded — run python -m piezo1.io.fetch")
    result = compare_fluctuations(st)
    assert result.available
    assert 0.3 < result.pearson_r < 0.6
    assert 0.6 < result.spearman_r < 0.85
    assert result.spearman_r > result.pearson_r + 0.2, (
        "the relationship is monotone but not linear, which is the reason "
        "both statistics are reported")
    assert result.beats_control and not result.control_inverted
    # Truncation matters but does not decide the answer.
    assert len(result.by_mode_count) >= 2
    assert all(0.2 < v < 0.7 for v in result.by_mode_count.values())


def test_the_two_entries_whose_column_disagrees_with_burial(structure_by_id):
    """8YEZ and 8ZU8 report a *negative* control, which is a verdict on them.

    More neighbours means less mobile in any packed solid, so a column that
    rises with burial is not reporting mobility. Both are the pair the registry
    calls the highest-resolution apo human entries, and both give the network
    essentially nothing (rank correlation 0.10) — the honest reading is that
    the column is not a temperature factor, not that the network failed.
    """
    st = structure_by_id("8YEZ")
    if st is None:
        pytest.skip("8YEZ not downloaded — run python -m piezo1.io.fetch")
    result = compare_fluctuations(st)
    assert result.available, "the column passes the quality gate"
    assert result.control_inverted
    assert result.control_spearman < -0.4
    assert abs(result.spearman_r) < 0.25


def test_observed_b_factors_are_read_per_residue_not_per_atom(structure_by_id):
    st = structure_by_id("7WLT")
    if st is None:
        pytest.skip("7WLT not downloaded")
    residues, values = observed_b_factors(st)
    assert len(residues) == len(values) == int(st.mask_ca().sum())
    assert values.min() > 0.0


def test_the_resolution_floor_is_registered_with_its_reason():
    parameter = PARAMETERS.get("fluctuation.min_distinct_fraction")
    assert parameter is not None
    assert parameter.citation == "method_choice" and parameter.source_note
    assert 0.0 < parameter.default < 0.5


def test_quality_carries_its_own_arithmetic():
    quality = BFactorQuality(n_residues=100, n_distinct=25, minimum=10.0,
                             maximum=60.0)
    assert quality.span == 50.0
    assert quality.distinct_fraction == 0.25
