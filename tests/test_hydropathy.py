"""Hydropathy and the 4-TM repeat — Guo & MacKinnon 2017, Figure 3-S1 to S3.

The repeat test is a checking instrument and is calibrated on planted signal
and on a true null before it is run on a real sequence. The register-maximised
control is the load-bearing part: taking the best of four registers and
comparing it against an unmaximised null manufactured about one standard
deviation of significance from nothing.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from piezo1.analysis.hydropathy import (KYTE_DOOLITTLE, annotated_hydropathy,
                                        compare_with_reference,
                                        hydropathy_profile, load_reference,
                                        predict_segments, repeat_periodicity,
                                        threshold_scan)
from piezo1.config import RESOURCE_DIR


# --------------------------------------------------------------------------
# The scale and the curve
# --------------------------------------------------------------------------

def test_the_scale_is_kyte_and_doolittles():
    """A scale with one wrong entry gives a curve that still looks right."""
    assert len(KYTE_DOOLITTLE) == 20
    assert KYTE_DOOLITTLE["I"] == 4.5 and KYTE_DOOLITTLE["R"] == -4.5
    assert KYTE_DOOLITTLE["V"] == 4.2 and KYTE_DOOLITTLE["K"] == -3.9
    assert KYTE_DOOLITTLE["G"] == -0.4 and KYTE_DOOLITTLE["W"] == -0.9
    assert np.mean(list(KYTE_DOOLITTLE.values())) == pytest.approx(-0.49,
                                                                   abs=0.005)


def test_a_uniform_sequence_gives_its_own_value_everywhere():
    profile = hydropathy_profile("I" * 100, window=19)
    assert np.allclose(profile.value, 4.5)


def test_the_window_is_truncated_at_the_termini_not_padded():
    """Padding with zeros would invent a hydrophilic N-terminus."""
    profile = hydropathy_profile("I" * 100, window=19)
    assert profile.value[0] == pytest.approx(4.5), "no dilution at the start"
    assert profile.meta["window_used"][0] == 10
    assert profile.meta["window_used"][50] == 19


def test_a_planted_hydrophobic_block_is_found_where_it_was_planted():
    sequence = "D" * 40 + "L" * 25 + "D" * 40
    profile = hydropathy_profile(sequence, window=19)
    segments = predict_segments(profile, threshold=1.6, min_length=10)
    assert len(segments) == 1
    assert 40 <= segments[0].start <= 50
    assert 55 <= segments[0].end <= 70


def test_an_unknown_residue_scores_zero_without_shifting_the_numbering():
    profile = hydropathy_profile("IIIXIII", window=1)
    assert profile.meta["n_unknown"] == 1
    assert profile.value[3] == pytest.approx(0.0)
    assert len(profile.value) == 7


# --------------------------------------------------------------------------
# Against the annotation
# --------------------------------------------------------------------------

def test_piezo1_helices_sit_below_the_conventional_membrane_cut():
    """The measurement that explains the poor recall, and why nothing is tuned.

    PIEZO1's transmembrane helices average about +1.2 on the Kyte-Doolittle
    scale against a conventional +1.6 cut. Tuning the threshold down until the
    agreement looked good would make the agreement a statement about the
    tuning.
    """
    stats = annotated_hydropathy()
    assert stats["n_helices"] == 38
    assert 0.9 < stats["mean_window_hydropathy"] < 1.5
    assert stats["mean_window_hydropathy"] < stats["default_threshold"]
    assert stats["fraction_above_default_threshold"] < 0.25
    # They are still clearly above their surroundings.
    assert stats["separation"] > 1.0


def test_the_threshold_scan_trades_recall_for_specificity_monotonically():
    rows = threshold_scan()
    recalls = [r["recall"] for r in rows]
    assert recalls == sorted(recalls, reverse=True)
    assert rows[0]["recall"] == pytest.approx(1.0)
    assert rows[-1]["recall"] < 0.1


def test_the_default_threshold_recovers_few_helices_and_that_is_reported():
    """A poor number, pinned so it cannot be quietly improved by tuning."""
    agreement = compare_with_reference()
    assert agreement.n_annotated == 38
    assert agreement.n_recovered <= 10
    assert agreement.n_spurious == 0


# --------------------------------------------------------------------------
# The repeat — the paper's actual inference
# --------------------------------------------------------------------------

def test_the_repeat_test_finds_a_planted_period(monkeypatch, tmp_path):
    """Known answer: helices with a long loop at every fourth position."""
    helices, position = [], 1
    for index in range(36):
        helices.append({"start": position, "end": position + 20,
                        "name": f"TM{index + 1}"})
        position += 21 + (150 if (index + 1) % 4 == 0 else 8)
    payload = {"length": position, "sequence": "A" * position,
               "transmembrane": helices}

    reference = tmp_path / "uniprot_planted.json"
    reference.write_text(json.dumps(payload))
    monkeypatch.setattr("piezo1.analysis.hydropathy.RESOURCE_DIR", tmp_path)

    result = repeat_periodicity(reference="planted", n_shuffles=500)
    assert result.supported
    assert result.phase == 3
    assert result.z > 4.0


def test_the_repeat_test_says_no_to_a_true_null(monkeypatch, tmp_path):
    """Uniform loops: no period exists and none must be reported.

    A check with no input that makes it say "no" is not a check.
    """
    helices, position = [], 1
    for index in range(36):
        helices.append({"start": position, "end": position + 20,
                        "name": f"TM{index + 1}"})
        position += 21 + 30
    payload = {"length": position, "sequence": "A" * position,
               "transmembrane": helices}
    (tmp_path / "uniprot_flat.json").write_text(json.dumps(payload))
    monkeypatch.setattr("piezo1.analysis.hydropathy.RESOURCE_DIR", tmp_path)

    result = repeat_periodicity(reference="flat", n_shuffles=500)
    assert not result.supported
    assert abs(result.contrast) < 1e-9


def test_the_control_is_register_maximised_like_the_statistic(monkeypatch,
                                                              tmp_path):
    """Best-of-four against an unmaximised null would invent significance.

    The control's mean must be well above zero, because maximising over four
    registers on random data gives a positive contrast by construction. A
    control mean near zero means the two are not being computed the same way.
    """
    result = repeat_periodicity(reference="mouse", n_shuffles=800)
    assert result.control_mean > 10.0, (
        "a register-maximised control cannot have a near-zero mean")
    assert result.control_sd > 0


@pytest.mark.parametrize("reference,expected", [("mouse", True),
                                                ("human", True),
                                                ("mouse_piezo2", True),
                                                ("worm_piezo", False),
                                                ("fly_piezo", False)])
def test_the_repeat_holds_in_the_mammalian_piezos_and_not_the_invertebrate_ones(
        reference, expected):
    """The measured result, across every PIEZO this project carries.

    Recorded as a finding rather than asserted about PIEZO1 alone: the 4-TM
    repeat Guo & MacKinnon infer from hydropathy is statistically supported in
    both mammalian PIEZOs and is not in PEZO-1 or dPIEZO, which is consistent
    with those two not sharing the 38-helix architecture either.
    """
    result = repeat_periodicity(reference=reference)
    assert result.supported is expected, result.summary()
    if expected:
        assert result.phase == 3, (
            "the long loop must fall after every fourth helix, which is what "
            "makes the units run TM1-4, TM5-8 and so on")


def test_the_units_the_repeat_supports_are_the_ones_domains_json_uses():
    """The inference is load-bearing for this project, not just for the paper."""
    result = repeat_periodicity(reference="mouse")
    assert result.n_units == 9 and result.n_helices == 36

    domains = json.loads((RESOURCE_DIR / "domains.json").read_text())["domains"]
    thus = [d for d in domains if d["id"].startswith("thu")]
    assert len(thus) == result.n_units
    for domain in thus:
        assert len(domain["transmembrane"]) == result.period


def test_too_few_helices_to_see_a_period_is_refused():
    with pytest.raises(ValueError, match="at least"):
        repeat_periodicity(reference="mouse", n_helices=4)


def test_a_missing_reference_names_what_it_looked_for():
    with pytest.raises(FileNotFoundError, match="uniprot_nonesuch"):
        load_reference("nonesuch")
