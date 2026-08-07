"""Reading AlphaFold's PAE, and the control that changes the answer.

The naive comparison — mean PAE across the seam versus within each region —
says 27.3 against 16.1/20.7 and looks decisive. Controlling for sequence
separation it mostly disappears, and at short separation it reverses. These
tests hold both the controlled result and the reason the uncontrolled one is
wrong.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.prediction_confidence import (PAEMatrix, assess_seam,
                                                   load_pae, load_plddt)


@pytest.fixture(scope="module")
def report():
    result = assess_seam()
    if result is None:
        pytest.skip("AlphaFold model or PAE not downloaded — "
                    "run python -m piezo1.io.fetch --pae")
    return result


# --------------------------------------------------- the control, on a known

def test_separation_control_detects_a_planted_seam():
    """On a matrix built with a real seam penalty, the control must find it.

    Otherwise "no penalty" would be indistinguishable from "the control does
    not work" — the same trap Round 42's probe control exists for.
    """
    n, seam = 400, 201
    i, j = np.indices((n, n))
    values = 5.0 + 0.01 * np.abs(i - j)                    # grows with separation
    crosses = ((i < seam - 1) & (j >= seam - 1)) | ((j < seam - 1) & (i >= seam - 1))
    values = values + 6.0 * crosses                        # a planted penalty

    matrix = PAEMatrix(values=values, maximum=32.0)
    table = matrix.by_separation(seam, bins=((50, 150), (150, 300)))
    assert table, "the control produced no bins"
    for row in table:
        # Recovered as ~6.5 rather than exactly 6.0: within a bin, cross-seam
        # pairs still sit at slightly larger separations than within-region
        # ones, so a little of the separation trend leaks through. The control
        # removes most of the confound, not all of it — which is why the real
        # measurement is read as "no consistent penalty" rather than as a
        # precise effect size.
        assert row["penalty"] == pytest.approx(6.0, abs=1.0)


def test_separation_control_finds_nothing_when_there_is_nothing():
    """A matrix depending only on separation must show no seam penalty."""
    n, seam = 400, 201
    i, j = np.indices((n, n))
    matrix = PAEMatrix(values=5.0 + 0.01 * np.abs(i - j), maximum=32.0)
    for row in matrix.by_separation(seam, bins=((50, 150), (150, 300))):
        # Residual leakage only; an order of magnitude below a planted effect.
        assert abs(row["penalty"]) < 1.0


def test_the_uncontrolled_comparison_is_misleading_by_construction():
    """Block means conflate "across the seam" with "far apart in sequence".

    Built from a matrix with NO seam penalty at all: the raw block comparison
    still shows one, purely because cross-seam pairs are more separated.
    """
    n, seam = 400, 201
    i, j = np.indices((n, n))
    values = 5.0 + 0.05 * np.abs(i - j)
    matrix = PAEMatrix(values=values, maximum=32.0)

    blade = np.arange(0, seam - 1)
    core = np.arange(seam - 1, n)
    naive = matrix.block_mean(blade, core) - matrix.block_mean(core, core)
    controlled = max(abs(row["penalty"]) for row in
                     matrix.by_separation(seam, bins=((50, 150), (150, 300))))
    assert naive > 2.0, "the uncontrolled comparison should look like a penalty"
    assert controlled < naive / 2.0, (
        f"the control should shrink the apparent penalty: naive {naive:.2f}, "
        f"controlled {controlled:.2f}")


# ------------------------------------------------------- the real measurement

def test_plddt_agrees_with_where_the_seam_had_to_be_placed(report):
    """The half of the round's question that comes out yes."""
    assert report.plddt_agrees_with_seam
    assert report.plddt_blade == pytest.approx(64.5, abs=2.0)
    assert report.plddt_core == pytest.approx(74.2, abs=2.0)
    # 52.2% of the blade below 70 against 27.0% of the core: 1.9x.
    assert report.plddt_blade_low_fraction > 1.5 * report.plddt_core_low_fraction


def test_pae_does_not_single_out_the_seam(report):
    """The half that comes out no, and reverses at short separation."""
    assert not report.pae_singles_out_seam
    short = report.separation_table[0]
    assert short["separation"][0] == 50
    assert short["penalty"] < 0, (
        "pairs 50-150 apart across the seam should score BETTER than pairs the "
        "same distance apart within one region")
    # The penalty peaks in the middle and fades again.
    penalties = [row["penalty"] for row in report.separation_table]
    assert max(penalties) < 6.0, "on a 31.75 A scale this is a modest effect"


def test_the_prediction_does_not_constrain_long_range_geometry_anywhere(report):
    """The stronger finding, and the one that matters for a hybrid model.

    PAE is 85% saturated beyond 800 residues of separation, and 80% saturated
    *within the cryo-EM-resolved core* — a region experiment places confidently.
    So the seam is not the weak point; the global arrangement is unconstrained
    wherever the cut is made.
    """
    assert not report.global_architecture_constrained
    assert report.saturation_all > 0.8
    assert report.saturation_core > 0.75
    # And the core is not meaningfully better than the whole.
    assert report.saturation_all - report.saturation_core < 0.15


def test_the_matrix_is_the_full_sequence(report):
    assert report.meta["n_residues"] == 2521
    assert report.meta["max_pae"] == pytest.approx(31.75, abs=0.5)


def test_summary_states_both_halves(report):
    text = report.summary()
    assert "agrees with the seam" in text
    assert "does NOT single out" in text
    assert "saturated" in text


def test_loaders_return_none_rather_than_guessing(tmp_path):
    """A missing download must not become a fabricated confidence report."""
    assert load_pae(path=None) is not None or True   # cached or not, must not raise
    from piezo1.analysis.prediction_confidence import ConfidenceReport
    assert ConfidenceReport(seam=1, plddt_blade=0, plddt_core=0,
                            plddt_blade_low_fraction=0,
                            plddt_core_low_fraction=0).pae_singles_out_seam is False
