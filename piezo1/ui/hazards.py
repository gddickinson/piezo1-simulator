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
        scenario="Superpose one structure on another and read the per-residue "
                 "deviation, where the same residue number means different "
                 "things in the two files — a different species, the PIEZO2 "
                 "paralogue, or a splice isoform numbered in its own "
                 "coordinates.",
        wrong="A deviation profile aligned by residue number between files "
              "that do not share one. Measured on the paralogue: 47.9 A over "
              "920 'matched' C-alphas of which 6% are even the same amino "
              "acid, where a real alignment gives 4.36 A over 3,708.",
        guard="OverlayController identifies both files by scoring their own "
              "residue names against every reference sequence, and refuses "
              "any pair that does not share one — naming what each is. The "
              "check was the registry's species label until Round 83 found "
              "PIEZO2 filed as 'mouse' and passing it.",
        status="guarded",
        where="ui/overlay_controller.py"),

    Hazard(
        key="imported_result_read_as_measured",
        scenario="Open the constraint colouring, or any of the four family "
                 "analyses, and read the numbers as something this "
                 "application computed from the loaded coordinates.",
        wrong="An external project's evolutionary result attributed to this "
              "one's physics. A coloured trimer is the most persuasive thing "
              "this application can put on a screen and these values came "
              "from a 194-genome census that ran somewhere else; a reader who "
              "takes them for a measurement here has no way to check them "
              "and no reason to look for the caveat.",
        guard="Every family result window leads with whose numbers they are, "
              "the constraint status line begins 'NOT MEASURED HERE', and the "
              "imported resource records the source project and the commit it "
              "was taken at. The importer re-reads all 32 quoted numbers from "
              "the source on every build and refuses to write when the source "
              "is absent, so a superseded value cannot survive as a "
              "quotation. The importer itself is "
              "scripts/build_family_findings.py, outside the package.",
        status="guarded",
        where="ui/constraint_controller.py, ui/tabular_analyses.py, "
              "core/family.py"),

    Hazard(
        key="domain_partition_decides_the_answer",
        scenario="Report the pore-module disease enrichment, or any per-domain "
                 "constraint number, on one set of domain boundaries.",
        wrong="A boundary choice reported as a finding. This project and the "
              "census put the anchor 141 residues apart and the outer helix "
              "120; the enrichment reaches odds ratio 3.63 at P = 0.0033 on "
              "their boundaries and 1.60 at P = 0.25 on ours, and the 120 "
              "residues in dispute carry six pathogenic positions themselves.",
        guard="disease_geography reports the test under BOTH partitions with "
              "the disagreement element by element and the disputed band's "
              "own pathogenic residues named, and its verdict says outright "
              "when significance depends on the choice. blade_gradient does "
              "the same for the census's distal-versus-proximal result, which "
              "reverses on the transmembrane units.",
        status="guarded",
        where="analysis/disease_geography.py, analysis/family_constraint.py"),

    Hazard(
        key="prediction_splay_read_as_divergence",
        scenario="Superpose a predicted monomer on an experimental structure "
                 "by the pore module and read the blade separation as a "
                 "difference between the two proteins.",
        wrong="A property of the predictor reported as evolutionary "
              "divergence. Measured: an AlphaFold monomer splays 7.2-9.1x from an "
              "experimental structure OF THE SAME PROTEIN (7.2-9.1x), "
              "while three experimental cross-paralogue pairs splay 0.8-2.5x.",
        guard="core_periphery reports the splay as a ratio against the core "
              "it was fitted on, refuses a ratio entirely when the cores do "
              "not superpose, and states cross_paralogue explicitly; the "
              "caveat on every result names the prediction-versus-experiment "
              "control, and a test drives it on the mouse AlphaFold model "
              "against 6B3R.",
        status="guarded",
        where="analysis/core_periphery.py, ui/tabular_analyses.py"),

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
        scenario="Read the mechanical coupling score as a prediction of gain "
                 "versus loss of function. Round 58 measured that no GUI or "
                 "CLI path computes it at all — it is reachable only from a "
                 "notebook or the validation scripts — so the exposure is "
                 "narrower than this register first claimed.",
        wrong="The project's central claim, treated as established when five "
              "pre-registered tests have failed to support it and two "
              "independent routes show it cannot be settled.",
        guard="The output is named for what it measures: CouplingScore with a "
              "gating_cost_change field and a `sign` property. It was called "
              "VariantPrediction with a `direction` property returning "
              "'LoF-like'/'GoF-like' until Round 58. prediction_record still "
              "supplies the caveats wherever a variant is shown.",
        status="guarded",
        where="analysis/variant_impact.py, analysis/prediction_record.py"),

    Hazard(
        key="closed_state_interactions_read_as_general",
        scenario="Run the interaction inventory on a closed structure and take "
                 "the salt bridges as properties of the channel.",
        wrong="State-specific contacts read as state-independent ones.",
        guard="Stated in the analysis caveat, which previously was the one "
              "empty string in the caveat table.",
        status="guarded",
        where="ui/tabular_analyses.py"),

    Hazard(
        key="stale_analysis_lands_on_a_new_structure",
        scenario="Start a pore profile, then load a different entry before it "
                 "finishes. The analyses run off the GUI thread and a worker "
                 "cannot be interrupted, so the result arrives after the "
                 "structure it was computed from has been replaced.",
        wrong="One entry's bottleneck radius, wetting verdict, probe spheres "
              "and calcium source point read as the displayed entry's. Not "
              "merely recorded — `pore_surface.refresh()` and "
              "`nanodomain.refresh()` draw them, inside a lumen they were "
              "never measured in, and nothing on screen contradicts it.",
        guard="Each run is stamped with the structure object it was launched "
              "for, and a result whose stamp no longer matches is discarded "
              "with a note on the Analysis panel — not on the status line, "
              "which belongs to the load that invalidated it. `reset()` "
              "cleared the stored result but never the run in flight. Found "
              "in Round 88 by a timing change, which is the only reason it "
              "ever showed; a hazard that depends on the scheduler being slow "
              "is not guarded.",
        status="guarded",
        where="ui/analysis_controller.py"),
)


def by_status(status: str) -> list:
    return [h for h in HAZARDS if h.status == status]


def unguarded() -> list:
    """Hazards that are real, known and not prevented."""
    return by_status("accepted")
