"""The one variant-versus-wild-type pair, and the control that interprets it.

Two structures always differ. What makes a single pair mean anything is whether
they differ by more than two wild-type structures do — and that control is only
honest if the wild-type set contains no coordinate duplicates.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.paired_variant import (VARIANT_ENTRY,
                                            WILD_TYPE_CANDIDATES,
                                            PairedComparison,
                                            StructuralMetrics, compare)


@pytest.fixture(scope="module")
def comparison():
    result = compare()
    if result is None or len(result.wild_type) < 2:
        pytest.skip("structures or CHAP grid not downloaded")
    return result


def test_duplicates_are_excluded_from_the_control(comparison):
    """8YFC and 9VMX share coordinates with 8ZU3.

    Including them would add two zero-difference pairs, narrowing the wild-type
    spread and making the variant look more distinct than it is. They are
    removed by fingerprint rather than by name, so a future duplicate is caught
    without anyone updating a list.
    """
    excluded = {pdb for pdb, _ in comparison.excluded_duplicates}
    assert excluded == {"8YFC", "9VMX"}
    kept = {m.pdb for m in comparison.wild_type}
    assert kept == {"8YEZ", "8ZU3", "8ZU8"}
    assert len(kept) + len(excluded) == len(WILD_TYPE_CANDIDATES)

    # And the fingerprints really are distinct among those kept.
    assert len({m.fingerprint for m in comparison.wild_type}) == len(kept)


def test_the_variant_falls_inside_the_wild_type_range(comparison):
    """The measured answer, on both metrics."""
    report = comparison.report()
    for name, values in report.items():
        assert values["within_wild_type_range"], name
        assert not values["exceeds_wild_type_spread"], name
        assert values["largest_variant_difference"] < values["wild_type_spread"]


def test_r2456h_is_not_structurally_distinguishable(comparison):
    """Stated generously: *any* measure exceeding the spread would count.

    It still does not. Bottleneck 0.808 A against a wild-type 0.673-0.930;
    wetting 0.904 against 0.457-0.986.
    """
    assert comparison.variant.pdb == VARIANT_ENTRY
    assert not comparison.distinguishable
    assert comparison.variant.bottleneck_A == pytest.approx(0.808, abs=0.02)
    assert comparison.variant.wetting_score == pytest.approx(0.904, abs=0.02)


def test_every_entry_is_closed_so_the_comparison_is_of_closed_states(comparison):
    """Which is why a gating variant need not show up here.

    R2456H's phenotype is slowed inactivation. A closed structure has no
    obligation to reveal it, and saying so is part of the result rather than an
    excuse for it.
    """
    for metrics in [comparison.variant] + comparison.wild_type:
        assert metrics.sterically_occluded, metrics.pdb
        assert metrics.bottleneck_A < 1.4, metrics.pdb


def test_the_summary_states_that_n_equals_one(comparison):
    text = comparison.summary()
    assert "n = 1" in text
    assert "supports no inference" in text
    assert "Distinguishable: False" in text


# ------------------------------------------------- the control, on knowns

def test_the_control_would_detect_a_variant_that_did_differ():
    """Otherwise "not distinguishable" and "the test cannot tell" look alike."""
    wild = [StructuralMetrics(pdb=f"WT{i}", bottleneck_A=0.90 + 0.01 * i,
                              wetting_score=0.80 + 0.01 * i,
                              hydrophobic_gate=True, sterically_occluded=True,
                              fingerprint=f"f{i}") for i in range(3)]
    far = StructuralMetrics(pdb="FAR", bottleneck_A=3.5, wetting_score=0.05,
                            hydrophobic_gate=False, sterically_occluded=False,
                            fingerprint="fx")
    outside = PairedComparison(variant=far, wild_type=wild)
    assert outside.distinguishable
    for values in outside.report().values():
        assert not values["within_wild_type_range"]
        assert values["exceeds_wild_type_spread"]

    near = StructuralMetrics(pdb="NEAR", bottleneck_A=0.905, wetting_score=0.805,
                             hydrophobic_gate=True, sterically_occluded=True,
                             fingerprint="fy")
    assert not PairedComparison(variant=near, wild_type=wild).distinguishable


def test_a_narrow_control_makes_anything_look_distinguishable():
    """Why the duplicate exclusion matters, made concrete.

    With a wild-type set of identical structures the spread is zero, and any
    difference at all clears it. That is exactly what including 8YFC and 9VMX
    would have moved the result towards.
    """
    identical = [StructuralMetrics(pdb=f"WT{i}", bottleneck_A=0.90,
                                   wetting_score=0.80, hydrophobic_gate=True,
                                   sterically_occluded=True, fingerprint="same")
                 for i in range(3)]
    barely = StructuralMetrics(pdb="V", bottleneck_A=0.901, wetting_score=0.801,
                               hydrophobic_gate=True, sterically_occluded=True,
                               fingerprint="v")
    assert PairedComparison(variant=barely, wild_type=identical).distinguishable
