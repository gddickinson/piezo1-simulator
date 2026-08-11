"""Every way to get a confident wrong number, driven until the guard fires.

The register in `piezo1.ui.hazards` is prose until something makes each guard
actually refuse. So each test here **constructs the dangerous situation** and
asserts the guard responds — a guard nobody has watched fail is not evidence
that it works.

Runs on the offscreen Qt platform with real widgets, following
`test_ui_analysis.py`, because the hazards Round 50 audits are hazards *of the
interface* and mocks would prove nothing about it.
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.parameters import PARAMETERS, reset, set_value  # noqa: E402
from piezo1.ui.hazards import HAZARDS, Hazard, by_status, unguarded  # noqa: E402

pytest.importorskip("PyQt6.QtWidgets")
from PyQt6.QtWidgets import QApplication  # noqa: E402


@pytest.fixture(scope="module")
def qapp():
    app = QApplication.instance()
    if app is None:
        try:
            app = QApplication([])
        except Exception as exc:                       # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return app


@pytest.fixture(autouse=True)
def clean_registry():
    reset()
    yield
    reset()


# ------------------------------------------------------- the register

def test_every_hazard_is_completely_specified():
    assert HAZARDS, "the register must not be empty"
    for hazard in HAZARDS:
        assert isinstance(hazard, Hazard)
        assert hazard.status in {"guarded", "by-design", "accepted"}
        for field in ("scenario", "wrong", "guard"):
            assert len(getattr(hazard, field)) > 30, (
                f"{hazard.key}.{field} is too vague to check")
        assert hazard.where, f"{hazard.key} does not say where its guard lives"


def test_an_accepted_hazard_says_what_to_do_instead():
    """Recording a hazard is not enough if a user is left with no recourse."""
    for hazard in unguarded():
        assert "instead" in hazard.guard.lower() or "use " in hazard.guard.lower(), (
            f"{hazard.key} is accepted but offers no alternative")


def test_the_guard_modules_named_in_the_register_exist():
    """A register pointing at a deleted module is worse than no register."""
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "piezo1"
    for hazard in HAZARDS:
        for path in (p.strip() for p in hazard.where.split(",")):
            assert (root / path).exists(), f"{hazard.key} names missing {path}"


# --------------------------------------- modified registry is never silent

def test_a_result_window_states_the_parameter_set_it_used():
    from piezo1.ui.result_dialog import provenance_line

    assert "documented defaults" in provenance_line("8YEZ", "human")
    set_value("pore.step", 0.5)
    warned = provenance_line("8YEZ", "human")
    assert "NON-DEFAULT" in warned
    assert "pore.step" in warned
    assert "docs/SCIENCE.md" in warned, (
        "the warning must say what the number cannot be compared with")


def test_the_stamp_records_the_state_at_compute_time_not_at_read_time(qapp):
    """The window is non-modal and can outlive the registry state.

    A stamp recomputed on read would quietly agree with whatever the registry
    says later, which is exactly the failure — the numbers in the window were
    computed under the old set.
    """
    from piezo1.ui.result_dialog import ResultDialog

    dialog = ResultDialog("t", {"conductance_pS": 27.4},
                          structure_name="8YEZ", species="human")
    assert "documented defaults" in dialog.provenance

    set_value("pore.step", 0.5)
    assert "documented defaults" in dialog.provenance, (
        "the stamp changed after the fact; it must describe the computation")


def test_a_number_shown_while_overridden_carries_the_warning(qapp):
    from piezo1.ui.result_dialog import ResultDialog

    set_value("pore.step", 0.5)
    dialog = ResultDialog("t", {"bottleneck_A": 0.76}, structure_name="8YEZ")
    assert "NON-DEFAULT" in dialog.provenance


def test_exported_reports_still_carry_their_own_banner():
    """The report path is independent of the dialog path; both must warn."""
    from piezo1.analysis.report import collect_provenance

    set_value("pore.step", 0.5)
    prov = collect_provenance()
    assert prov.parameter_overrides
    assert any("NOT comparable" in w for w in prov.warnings)


def test_verify_claims_refuses_against_a_modified_registry():
    """It refuses outright rather than returning results marked incomparable.

    The stronger behaviour: a caller who ignores a flag still cannot get a
    drift report computed under someone else's parameters, and the message
    names the way to do it on purpose.
    """
    from piezo1.analysis.claims import CLAIMS, verify_claims

    set_value("kinetics.t50_measured", 3.5)
    claim = [c for c in CLAIMS if c.key == "kinetics.t50"]
    with pytest.raises(RuntimeError, match="cannot verify documented numbers"):
        verify_claims(claims=claim, verbose=False)

    # And it can still be done deliberately, which is what makes the refusal a
    # guard rather than an obstruction.
    results = verify_claims(claims=claim, verbose=False, allow_overrides=True)
    assert results and not results[0].comparable


# ------------------------------- companion structures and the primary

def test_a_result_window_names_the_structure_it_used(qapp):
    """With several structures displayed, the number must say which one."""
    from piezo1.ui.result_dialog import ResultDialog

    dialog = ResultDialog("Ion permeation", {"g_pS": 27.4},
                          structure_name="11ZC", species="human")
    assert "11ZC" in dialog.provenance
    assert "human numbering" in dialog.provenance


def test_an_unnamed_structure_is_reported_as_unknown_not_omitted():
    from piezo1.ui.result_dialog import provenance_line

    assert "unknown structure" in provenance_line("", "")


def test_analyses_use_the_primary_structure_not_a_companion():
    """Read off the code path, which is the thing that must stay true."""
    import inspect

    from piezo1.ui import tabular_analyses

    source = inspect.getsource(tabular_analyses.TabularAnalysisMixin
                               ._run_registry_analysis)
    assert "self.structure" in source
    assert "companion" not in source.lower(), (
        "the analysis path must not consult the companion map")


# ----------------------------------------- cross-species numbering

def test_the_overlay_refuses_a_pair_that_shares_no_numbering():
    """Driven rather than grepped: the guard is watched refusing each case.

    Reading the source for a variable name asserted that a particular
    *implementation* was present, which is how it went stale — the guard was
    rewritten in Round 83 to identify the protein rather than trust the
    registry's species label, and the grep failed while the hazard was more
    strongly guarded than before. Three pairs must be refused and one must go
    through.
    """
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure
    from piezo1.ui.overlay_controller import OverlayController

    def load(pdb):
        path = STRUCTURE_DIR / f"{pdb}.cif"
        if not path.exists():
            pytest.skip(f"{pdb} not downloaded — run python -m piezo1.io.fetch")
        return Structure.from_file(path)

    class _Window:
        def __init__(self, structure):
            self.structure, self.record = structure, None

        def _set_status(self, message):
            self.status = message

    reference = load("7WLT")
    controller = OverlayController(_Window(reference))

    for pdb, expected in (("8YEZ", "different species"),
                          ("6KG7", "different proteins"),
                          ("6LQI", "canonical numbering")):
        refusal = controller._numbering_refusal(pdb, load(pdb))
        assert refusal and expected in refusal, f"{pdb}: {refusal!r}"

    # ...and the pair the overlay exists for still goes through.
    assert controller._numbering_refusal("7WLU", load("7WLU")) == ""


def test_no_constant_offset_exists_between_the_two_numbering_systems():
    """The reason the overlay must refuse, measured rather than asserted."""
    from piezo1.core.sequence import human_to_mouse

    offsets = {human_to_mouse(r) - r for r in (600, 1200, 1800, 2400)
               if human_to_mouse(r) is not None}
    assert len(offsets) > 1, (
        "a single offset would make cross-species alignment safe; it is not")


# ------------------------------------------- caveats and the record

def test_no_tabular_analysis_is_shown_without_a_caveat():
    """`interactions` was an empty string until Round 50."""
    from piezo1.ui.tabular_analyses import CAVEATS

    for key, text in CAVEATS.items():
        assert text.strip(), f"{key} is displayed with no caveat"
        assert len(text) > 60, f"{key}'s caveat is too short to be a warning"


def test_the_prediction_caveat_still_refuses_to_soften():
    from piezo1.analysis.prediction_record import headline, what_it_means

    assert "not predict" in headline().lower()
    joined = " ".join(what_it_means()).lower()
    assert "five tests, five" in joined or "five nulls" in joined


def test_hidden_categories_cannot_change_what_is_measured():
    """Visibility is a render property; the entity map decides what is used."""
    import inspect

    from piezo1.ui import appearance

    source = inspect.getsource(appearance)
    for banned in ("ANALYSES[", "pore_profile(", "detect_interactions("):
        assert banned not in source, (
            f"appearance.py calls {banned} — visibility must not reach an analysis")


def test_sessions_record_what_was_viewed_and_never_a_result():
    import inspect

    from piezo1.io import session

    source = inspect.getsource(session)
    for banned in ("bottleneck", "conductance", "cliffs_delta", "wetting"):
        assert banned not in source, (
            f"session.py mentions {banned}; sessions must not carry results")


def test_the_register_covers_the_three_hazards_the_round_named():
    keys = {h.key for h in HAZARDS}
    assert {"cross_species_overlay", "companion_mistaken_for_primary",
            "modified_registry_unmarked"} <= keys
    assert not by_status("accepted"), (
        "an accepted hazard must be listed in ROADMAP.md as well")
