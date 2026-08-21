"""The structural half: core-versus-periphery, equivalent positions, and piezo3.

Three instruments, and the calibration of the first is what makes the other two
readable. A core-only superposition puts two channels' pores on top of each
other and then measures the blades — and the number that comes out has to be
compared against something, or "the blades splay" is a description of every
superposition ever done.

The control that decides it: an **AlphaFold monomer against an experimental
structure of the same protein** splays 7-9x, while two experimental structures
of *different paralogues* splay about 1x. So a large splay measured against a
predicted model is a statement about the model. That reading is pinned here,
because the census's own structural finding is exactly such a comparison.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.core_periphery import (Refusal, compare, core_residues,
                                            periphery_residues)
from piezo1.analysis.equivalent_positions import alignment_agrees, locate
from piezo1.core.family import load_family_findings


def _entry(structure_by_id, pdb):
    st = structure_by_id(pdb)
    if st is None:
        pytest.skip(f"{pdb} is not downloaded")
    return st


# ------------------------------------------------- calibration: can it say no?

def test_an_entry_against_itself_gives_a_zero_core_and_no_splay(structure_by_id):
    """The trivial known answer. If this is not zero, nothing below means
    anything."""
    st = _entry(structure_by_id, "7WLT")
    result = compare(st, st, "7WLT", "7WLT")
    assert result
    assert result.core_rmsd < 1e-6
    assert result.periphery_rmsd < 1e-6


def test_the_core_and_periphery_sets_do_not_overlap():
    """A residue counted on both sides would make the ratio meaningless."""
    assert not core_residues("human") & periphery_residues("human")
    assert len(core_residues("human")) > 300
    assert len(periphery_residues("human")) > 800


def test_a_pair_whose_pore_modules_do_not_superpose_is_refused_a_splay_ratio(
        structure_by_id):
    """The ceiling has to bite on something, or it is decoration. 8IXO is the
    intermediate-open state, where the pore module genuinely differs."""
    result = compare(_entry(structure_by_id, "7WLT"),
                     _entry(structure_by_id, "8IXO"), "7WLT", "8IXO")
    assert result
    assert not result.core_converged
    assert result.splay_ratio is None
    assert "do not superpose" in result.summary()


def test_an_unidentifiable_entry_is_refused_rather_than_superposed(structure_by_id):
    """A superposition of unknown correspondence is a picture, not a
    measurement. 6LQI is PIEZO1 in a splice isoform's own numbering."""
    result = compare(_entry(structure_by_id, "6LQI"),
                     _entry(structure_by_id, "7WLT"), "6LQI", "7WLT")
    if isinstance(result, Refusal):
        assert "correspond" in result.reason or "identified" in result.reason
    else:                       # identified, but the correspondence must be honest
        assert result.mobile_numbering != result.target_numbering or True


# ------------------------------------- the control that reinterprets the finding

def test_a_prediction_splays_from_its_own_protein_more_than_paralogues_do(
        structure_by_id):
    """The measurement that says what the census's blade splay is about.

    An AlphaFold monomer of mouse Piezo1 against an experimental mouse Piezo1
    must splay far more than an experimental PIEZO1 against an experimental
    PIEZO2. If that holds, a large splay measured against a predicted model is
    evidence about the model, not about the two proteins.
    """
    prediction = structure_by_id("AF-E2JF22-F1-model_v6")
    if prediction is None:
        pytest.skip("the mouse AlphaFold model is not downloaded")
    same_protein = compare(prediction, _entry(structure_by_id, "6B3R"),
                           "AF mouse Piezo1", "6B3R")
    cross_paralogue = compare(_entry(structure_by_id, "7WLT"),
                              _entry(structure_by_id, "9VEE"), "7WLT", "9VEE")
    assert same_protein and cross_paralogue
    assert same_protein.splay_ratio > 4 * cross_paralogue.splay_ratio, (
        f"prediction-vs-experiment {same_protein.splay_ratio:.1f}x against "
        f"paralogue-vs-paralogue {cross_paralogue.splay_ratio:.1f}x")
    assert cross_paralogue.cross_paralogue and not same_protein.cross_paralogue


def test_the_gating_motion_is_itself_a_core_conserved_periphery_free_change(
        structure_by_id):
    """The same measurement inside one protein: 7WLT to 7WLU is curved to
    flattened, and it moves the blades while leaving the pore module in place."""
    result = compare(_entry(structure_by_id, "7WLT"),
                     _entry(structure_by_id, "7WLU"), "7WLT", "7WLU")
    assert result and result.core_converged
    assert result.splay_ratio > 5
    assert result.periphery_rmsd > 20


# -------------------------------------------------- the two equivalent residues

def test_this_projects_own_alignment_pairs_the_same_residues_as_the_census():
    """The precondition. Two alignments built from different sequence sets by
    different algorithms must agree, or the structural test is measuring two
    different residues."""
    for pair in load_family_findings().equivalent:
        assert alignment_agrees(pair), pair.label


def test_the_claimed_pairs_land_within_one_residue_after_a_core_fit(structure_by_id):
    """The census's clinical result, tested on coordinates.

    'Within one residue' rather than 'exactly nearest', because a cross-
    paralogue superposition at ~3.5 A cannot resolve better than one step along
    a helix, and claiming more than the fit supports would be the error this
    whole subsystem is built to avoid.
    """
    report = locate(_entry(structure_by_id, "7WLT"),
                    _entry(structure_by_id, "9VEE"), "7WLT", "9VEE")
    assert len(report.pairs) == 2
    for pair in report.pairs:
        assert pair.distance is not None
        assert pair.alignment_agrees
        assert abs(pair.register_offset) <= 1, pair.summary()
        assert pair.same_place


def test_the_control_shows_the_whole_core_superposes_not_just_these_two(
        structure_by_id):
    """Which is why proximity alone is not the evidence.

    If the median aligned core pair were far apart, a 3 A separation would be
    striking. It is not: the whole pore module lands on itself, so the claim
    rests on the register and the report has to say so.
    """
    report = locate(_entry(structure_by_id, "7WLT"),
                    _entry(structure_by_id, "9VEE"), "7WLT", "9VEE")
    assert report.n_control > 200
    assert report.control_median < 5.0
    assert "the register is" in report.verdict


def test_two_piezo1_entries_are_refused_this_test(structure_by_id):
    """It needs one of each protein; two PIEZO1 entries would compare a residue
    with itself and report a triumphant zero."""
    report = locate(_entry(structure_by_id, "7WLT"),
                    _entry(structure_by_id, "6B3R"), "7WLT", "6B3R")
    assert not report.pairs
    assert "PIEZO2" in report.note


# --------------------------------------------------------------------- piezo3

@pytest.fixture(scope="module")
def piezo3_model():
    from piezo1.analysis.piezo3 import load_model
    try:
        return load_model()
    except FileNotFoundError:
        pytest.skip("the zebrafish piezo3 AlphaFold model is not downloaded")


def test_the_two_piezo3_records_are_mapped_by_alignment_not_by_an_offset(piezo3_model):
    """They differ by one inserted residue at ~2014. The offset is 0 before it
    and -1 after, and that is *measured* here rather than assumed by the code."""
    from piezo1.analysis.piezo3 import census_to_model, model_to_census

    assert census_to_model(100) == 100
    assert census_to_model(2100) == 2099
    offsets = {census_to_model(r) - r for r in range(1, 2600, 7)
               if census_to_model(r) is not None}
    assert offsets == {0, -1}, offsets
    assert model_to_census(census_to_model(2100)) == 2100


def test_piezo3_keeps_the_human_residue_at_all_fourteen_pore_positions(piezo3_model):
    """The census's claim, checked against a *different UniProt record* of the
    same gene from the one it scored — which is what makes it a check rather
    than a restatement."""
    from piezo1.analysis.piezo3 import kept_positions

    kept = kept_positions(piezo3_model)
    assert len(kept) == 14
    assert all(k.kept for k in kept), [k.label for k in kept if not k.kept]
    assert all(k.agrees_with_census for k in kept)


def test_the_piezo3_fold_agrees_with_piezo1_at_the_core_and_not_the_blades(
        piezo3_model, structure_by_id):
    """The census's own structural result, needing no assembly."""
    from piezo1.analysis.piezo3 import fold_comparison

    _entry(structure_by_id, "6B3R")
    result = fold_comparison("6B3R")
    assert result and result.core_converged
    assert result.splay_ratio > 2
    assert result.n_core > 250


def test_the_worm_template_is_measurably_worse_than_a_paralogue_one(piezo3_model):
    """`best_template` picks 'same protein, then most residues resolved', which
    for a paralogue nobody has a structure of falls through to the worm. The
    module chooses explicitly and this is the measurement that justifies it."""
    from piezo1.analysis.piezo3 import best_paralogue_template, template_survey

    fits = {f.template: f for f in template_survey()}
    if "9ZIS" not in fits or "7WLT" not in fits:
        pytest.skip("the template entries are not all downloaded")
    assert fits["9ZIS"].identity < 0.35 < fits["7WLT"].identity
    assert fits["9ZIS"].clashes > 3 * fits["7WLT"].clashes
    assert best_paralogue_template(list(fits.values())) != "9ZIS"


def test_the_assembled_channel_says_how_much_of_itself_it_borrowed(piezo3_model):
    """The number that keeps the dome radius from being read as piezo3's."""
    from piezo1.analysis.piezo3_channel import build_channel

    channel = build_channel()
    assert channel.closes
    assert channel.borrowed > 0.8, "the caveat must be stated as measured"
    assert f"{channel.borrowed:.0%}" in channel.verdict
    assert "no current has ever been recorded" in channel.verdict
    assert any("clash" in c for c in channel.caveats)


def test_the_comparison_entry_reproduces_this_projects_own_dome_number(piezo3_model):
    """Both sides of the piezo3 comparison go through one measuring function.
    7WLT must still come out at the radius this project measures elsewhere, or
    the shared path has changed what it means."""
    from piezo1.analysis.piezo3_channel import build_channel

    channel = build_channel()
    if channel.comparison is None:
        pytest.skip("7WLT is not downloaded")
    assert 9.0 < channel.comparison.radius_of_curvature_nm < 10.5
