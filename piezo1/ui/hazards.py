"""Ways a user could get a confident wrong number, and what stops each one.

Round 33 audited the menus for *reachability* — could a user get to every
analysis. This is the same exercise pointed at *correctness*: given that they
can reach it, can they be handed a number that is wrong without anything saying
so. A wrong number that announces itself is a bug; a wrong number that looks
exactly like a right one is the failure this project cares about.

Each entry names the scenario, what the user would wrongly conclude, and the
specific guard. ``tests/test_hazards.py`` exercises every guarded entry with a
**positive control**: it constructs the dangerous situation and asserts the
guard actually fires, because a guard nobody has seen fail is not evidence.

Kept as data, and Qt-free, so the help window, the tests and this audit cannot
drift apart — the same reason `prediction_record` and `claims` are data.

Status values:

``guarded``
    A mechanism exists and a test drives it.
``by-design``
    The situation cannot arise, because of how the code is shaped rather than
    because something checks at runtime.
``accepted``
    Real, not fixed, and written down. An accepted hazard must say what a user
    should do instead.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["Hazard", "HAZARDS", "by_status", "unguarded"]


@dataclass(frozen=True)
class Hazard:
    """One way to be handed a confident wrong number."""

    key: str
    scenario: str          # what the user does
    wrong: str             # what they would wrongly conclude
    guard: str             # what prevents or marks it
    status: str            # guarded | by-design | accepted
    where: str = ""        # the module the guard lives in

    def summary(self) -> str:
        return f"[{self.status}] {self.key}: {self.scenario} — {self.guard}"


HAZARDS: tuple = (
    Hazard(
        key="cross_species_overlay",
        scenario="Superpose a mouse structure on a human one and read the "
                 "per-residue deviation, where residue 2456 means different "
                 "things in the two files.",
        wrong="A deviation profile aligned by residue number across species, "
              "which is off by a non-constant offset.",
        guard="OverlayController refuses a pair whose numbering species "
              "differ, naming both, rather than silently aligning them.",
        status="guarded",
        where="ui/overlay_controller.py"),

    Hazard(
        key="companion_mistaken_for_primary",
        scenario="Display three structures at once, then run an analysis and "
                 "read the number as belonging to whichever is being looked at.",
        wrong="A measurement attributed to the wrong structure — the most "
              "likely way to compare two states and get the comparison "
              "backwards.",
        guard="Analyses always run on the primary, and every result window "
              "now states which structure and numbering it used.",
        status="guarded",
        where="ui/tabular_analyses.py, ui/result_dialog.py"),

    Hazard(
        key="modified_registry_unmarked",
        scenario="Override a parameter, close the dialog, run an analysis, "
                 "and read the result as the documented one.",
        wrong="A number that cannot be compared with docs/SCIENCE.md, taken "
              "for one that can.",
        guard="The status bar carries a persistent warning, exported reports "
              "carry a banner, verify_claims refuses to run, and every result "
              "window states the parameter set it was computed under.",
        status="guarded",
        where="ui/preferences.py, ui/result_dialog.py, analysis/report.py"),

    Hazard(
        key="stale_result_after_parameter_change",
        scenario="Open a result window, then change a parameter. The window is "
                 "non-modal and keeps showing the old numbers.",
        wrong="Numbers read as current when the registry has moved under them.",
        guard="The window records the parameter set at the moment it computed, "
              "so a stale window disagrees visibly with the status bar rather "
              "than matching it.",
        status="guarded",
        where="ui/result_dialog.py"),

    Hazard(
        key="hidden_category_changes_an_analysis",
        scenario="Hide lipids or ions for a clearer picture, then run an "
                 "analysis and assume the drawing is what was measured.",
        wrong="A pore radius or interaction count that silently depends on "
              "what happens to be visible.",
        guard="Visibility is a render property only; analyses always select "
              "the channel protomers from the entity map.",
        status="by-design",
        where="core/entities.py, ui/appearance.py"),

    Hazard(
        key="session_replays_a_stale_result",
        scenario="Save a session, change the code or the parameters, reload, "
                 "and see the earlier numbers.",
        wrong="A result attributed to the current code that predates it.",
        guard="Sessions record what was being *viewed* — structure, style, "
              "camera, selection — and never a computed value.",
        status="by-design",
        where="io/session.py"),

    Hazard(
        key="variant_numbering_on_a_mouse_structure",
        scenario="Map human-numbered curated variants onto a mouse structure.",
        wrong="Variants placed at the wrong residues, by a non-constant "
              "offset that looks plausible everywhere.",
        guard="Every conversion goes through core.sequence's alignment map; "
              "no constant offset exists anywhere in the codebase.",
        status="by-design",
        where="core/sequence.py"),

    Hazard(
        key="unresolved_residues_read_as_measured",
        scenario="Run a per-residue analysis on a structure that does not "
                 "resolve the residues of interest.",
        wrong="An absence of signal read as a measurement of zero.",
        guard="Unmeasured residues take the colour map's floor rather than "
              "zero, and coverage is reported alongside variant predictions.",
        status="guarded",
        where="ui/analysis_controller.py, analysis/variant_impact.py"),

    Hazard(
        key="prediction_read_as_validated",
        scenario="Select a variant, see a mechanical ΔΔG, and take it as a "
                 "prediction of gain versus loss of function.",
        wrong="The project's central claim, treated as established when five "
              "pre-registered tests have failed to support it.",
        guard="prediction_record supplies the headline and caveats to the GUI, "
              "the CLI and the tests from one frozen record.",
        status="guarded",
        where="analysis/prediction_record.py"),

    Hazard(
        key="closed_state_interactions_read_as_general",
        scenario="Run the interaction inventory on a closed structure and take "
                 "the salt bridges as properties of the channel.",
        wrong="State-specific contacts read as state-independent ones.",
        guard="Stated in the analysis caveat, which previously was the one "
              "empty string in the caveat table.",
        status="guarded",
        where="ui/tabular_analyses.py"),
)


def by_status(status: str) -> list:
    return [h for h in HAZARDS if h.status == status]


def unguarded() -> list:
    """Hazards that are real, known and not prevented."""
    return by_status("accepted")
