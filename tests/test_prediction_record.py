"""What a variant prediction is worth, and that the application says so.

The project's central claim has failed three pre-registered tests. That fact
lived in `docs/` and the CLI, where a GUI user would never meet it. These tests
hold the record to the runs that produced it and check the interface still
carries the caveat next to the number.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.prediction_record import (VALIDATION_RECORD, headline,
                                               variant_evidence, verify_record,
                                               what_it_means)


def test_the_record_has_all_three_tests_in_order():
    rounds = [e.round for e in VALIDATION_RECORD]
    assert rounds == [7, 22, 36]
    assert all(not e.rejected for e in VALIDATION_RECORD), (
        "every test failed to reject; if one ever does, this module and the "
        "headline both need rewriting rather than patching")


def test_the_record_matches_the_run_that_produced_it():
    """A frozen record that drifts from its own run is worse than none.

    Same discipline `analysis.claims` applies to the documented numbers: the
    stored Round 36 result is the authority, and these constants must track it.
    """
    check = verify_record()
    if not check["checked"]:
        pytest.skip(check["reason"])
    assert check["agrees"], f"the record has drifted: {check['drift']}"
    assert check["drift"]["n_gof"] == 0 and check["drift"]["n_lof"] == 0


def test_the_effect_grew_across_the_tests_and_none_reached_significance():
    """The pattern that is tempting to over-read, pinned so it stays labelled."""
    deltas = [e.cliffs_delta for e in VALIDATION_RECORD]
    assert deltas == sorted(deltas, reverse=True), (
        "the point estimate should be monotonically more negative")
    assert all(d < 0 for d in deltas)
    for entry in VALIDATION_RECORD:
        if entry.p_value is not None:
            assert entry.p_value >= 0.05


def test_the_headline_states_the_conclusion_plainly():
    text = headline()
    assert "three nulls" in text
    assert "does not predict" in text
    # It must carry the numbers, not just an adjective.
    assert "-0.249" in text or "−0.249" in text


def test_the_caveats_say_what_the_score_may_still_be_used_for():
    caveats = what_it_means()
    joined = " ".join(caveats).lower()
    assert len(caveats) >= 4
    assert "not evidence" in joined            # about the growing trend
    assert "power" in joined                   # scope of the null
    assert "data, not method" in joined        # where the constraint is
    # And it must not end on despair: there is a legitimate use.
    assert "do not use it to assign a direction" in joined
    assert "mechanically coupled" in joined


# ------------------------------------------------------------ per variant

def test_evidence_level_is_reported_per_variant():
    measured = variant_evidence("R2456H")
    assert measured["in_analysis_set"]
    assert measured["direction"] == "GoF"
    assert measured["evidence"] == "measured"
    assert "electrophysiology" in measured["evidence_note"]


def test_a_variant_outside_the_analysis_set_says_so():
    """E756del is a deletion, so it carries no direction and no prediction."""
    absent = variant_evidence("E756del")
    assert not absent["in_analysis_set"]
    assert absent["direction"] is None
    assert absent["evidence"] == "none"


def test_the_one_conflicting_variant_is_flagged():
    """V598M: curated says GoF, the disease mechanism implies LoF.

    The project reports the disagreement rather than resolving it, and Round 36
    excluded the variant in writing beforehand. A user meeting it in the
    interface should see the conflict, not a confident label.
    """
    conflict = variant_evidence("V598M")["conflict"]
    assert conflict is not None
    assert {conflict["curated"], conflict["inferred"]} == {"GoF", "LoF"}


# -------------------------------------------------------------- the GUI

def test_the_panel_shows_the_evidence_level_beside_the_variant(qt_app):
    """A classification alone reads as a fact, and for 20 of 46 it is inferred."""
    from piezo1.ui.panels.annotation_panel import AnnotationPanel

    panel = AnnotationPanel("human")
    try:
        html = panel._evidence_html("R2456H")
        assert "measured" in html
        assert "electrophysiology" in html

        conflicted = panel._evidence_html("V598M")
        assert "sources disagree" in conflicted
        assert "GoF" in conflicted and "LoF" in conflicted

        # A variant with no direction must add nothing rather than guess.
        assert panel._evidence_html("E756del") == ""
    finally:
        panel.deleteLater()


def test_the_record_is_reachable_from_the_gui_and_the_cli():
    """It is the central claim; a user must be able to meet its record."""
    from piezo1.analysis.report import ANALYSES
    from piezo1.ui.main_window import MainWindow

    assert "prediction_record" in ANALYSES
    assert callable(getattr(MainWindow, "show_prediction_record", None))

    result = ANALYSES["prediction_record"](None, "human")
    assert len(result["tests"]) == 3
    assert "does not predict" in result["verdict"]
    assert result["what_this_means"]


def test_the_dialog_caveat_does_not_soften_the_result():
    from piezo1.ui.tabular_analyses import CAVEATS

    caveat = CAVEATS["prediction_record"]
    assert "FAILED three" in caveat
    assert "less than it looks" in caveat
