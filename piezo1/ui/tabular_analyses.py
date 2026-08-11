"""Analyses that produce a table rather than something drawn on the model.

Each of these already existed in the shared ``ANALYSES`` registry and was
reachable from the command line only. That is a real gap rather than a cosmetic
one: a GUI user could not get at the ion current, the interaction inventory or
the variant-structure survey at all.

They run on the **primary** structure, whatever else is displayed, and each
carries the caveat its result needs rather than leaving it in a docstring the
user will not read.
"""

from __future__ import annotations

from .result_dialog import ResultDialog

__all__ = ["TabularAnalysisMixin"]

#: What each result must say about itself when shown. Kept beside the call
#: rather than inside the analysis, because the analysis returns data and this
#: is the sentence a reader needs in front of it.
CAVEATS = {
    "permeation": (
        "Continuum model of an atomic-scale pore. The in-pore diffusivity and "
        "the ion radius are UNMEASURED, and the computed conductance spans "
        "16-94 pS across their plausible ranges, so agreement with the "
        "published 25-30 pS would be tuning rather than prediction. The "
        "SELECTIVITY block reports two routes to the pore's fixed charge that "
        "bracket the measured P_Cl/P_Na tenfold apart rather than reproducing "
        "it, from an uncharged baseline that is already cation-selective from "
        "size alone; the curated route reaches an in-pore concentration no "
        "solution could hold, which is flagged on the result."),
    "labelling": (
        "Kinetics imported unchanged from the halotag_binding_sim project and "
        "reproduced to machine precision. The linker length and the reactive "
        "fraction are UNVERIFIED assumptions."),
    "hybrid": (
        "TWO POPULATIONS, NOT ONE STRUCTURE. Roughly 570-2521 is experimental; "
        "the remaining 569 residues are an AlphaFold PREDICTION, only 48% of "
        "which clears pLDDT 70. The seam fits to 2.4 A but the two models "
        "differ by 75 A overall, so a good local fit says nothing about the "
        "rest. Never average across the join."),
    "fusion": (
        "PLACEMENT IS A MODEL. There is no structure of the PIEZO1-HaloTag "
        "fusion; the tag body is a sphere of its radius of gyration and the "
        "linker length is unverified."),
    "interactions": (
        "Contacts are those of THIS structure in THIS state. A closed-state "
        "entry does not show the open-state salt bridges, and unresolved side "
        "chains cannot contribute a bond. Geometric criteria are heavy-atom "
        "based because deposited entries carry no hydrogens."),
    "nanodomain": (
        "The tag distance is MODELLED, not measured, and the calcium share of "
        "the current is unverified. Every deposited human structure is closed, "
        "so the current is borrowed from the open-like 11ZC and labelled as "
        "such."),
    "prediction_record": (
        "This is the project's central claim and it has FAILED three "
        "pre-registered tests. The record below is what a variant score from "
        "this application is entitled to claim, which is less than it looks."),
    "ligands": (
        "Every binding site here is INFERRED from mutagenesis, docking or "
        "geometry. No PIEZO structure with a bound small-molecule modulator "
        "has ever been deposited, so none of these pockets has been observed."),
    "paired_variant": (
        "n = 1. Only one deposited variant entry resolves its own mutation, so "
        "this is the single structural comparison available. It says what the "
        "structures show, not what R2456H does."),
    "variant_structures": (
        "A null result, reported rather than worked around: every deposited "
        "human PIEZO1 structure is closed, so no difference in conductance can "
        "be measured between variants."),
}


class TabularAnalysisMixin:
    """Menu-driven analyses whose output is a table."""

    def _show_result(self, key: str, title: str, data) -> None:
        # Name the structure the numbers came from. With companions displayed
        # there is otherwise nothing on the window saying which one it is, and
        # analyses always use the primary whatever else is drawn.
        name = getattr(self.structure, "name", "") if self.structure else ""
        species = self.record.numbering_species if self.record else ""
        dialog = ResultDialog(title, data, CAVEATS.get(key, ""), self,
                              structure_name=name, species=species)
        # Held on the window so it is not garbage-collected while open.
        if not hasattr(self, "_result_dialogs"):
            self._result_dialogs = []
        self._result_dialogs = [d for d in self._result_dialogs if d.isVisible()]
        self._result_dialogs.append(dialog)
        dialog.show()

    def _run_registry_analysis(self, key: str, title: str) -> None:
        from ..analysis.report import ANALYSES

        if self.structure is None:
            self._set_status("load a structure first")
            return
        species = self.record.numbering_species if self.record else "human"
        self._set_status(f"running {key}…")
        try:
            data = ANALYSES[key](self.structure, species)
        except Exception as exc:
            self._set_status(f"{key} failed: {exc}")
            return
        self._show_result(key, title, data)
        self._set_status(f"{key}: done")

    # --------------------------------------------------------------- entries

    def show_permeation(self) -> None:
        self._run_registry_analysis(
            "permeation", "Ion permeation through the pore")

    def show_interactions(self) -> None:
        self._run_registry_analysis(
            "interactions", "Hydrogen bonds, salt bridges and disulfides")

    def show_labelling(self) -> None:
        self._run_registry_analysis("labelling", "HaloTag labelling kinetics")

    def show_fusion_numbers(self) -> None:
        self._run_registry_analysis("fusion", "HaloTag fusion geometry")

    def show_hybrid(self) -> None:
        """The numbers behind View -> Full-length model.

        Drawn *and* tabulated: the picture shows where the graft is, and only
        the table gives the confident fraction and the 75 A disagreement that
        say how far to trust it.
        """
        self._run_registry_analysis(
            "hybrid", "Full-length model — experimental core plus prediction")

    def show_nanodomain(self) -> None:
        self._run_registry_analysis(
            "nanodomain", "Calcium nanodomain at the HaloTag")

    def show_ligands(self) -> None:
        self._run_registry_analysis(
            "ligands", "Modulators — and how much is known about their sites")

    def show_prediction_record(self) -> None:
        self._run_registry_analysis(
            "prediction_record",
            "Can this predict gain- vs loss-of-function? — the record")

    def show_paired_variant(self) -> None:
        self._run_registry_analysis(
            "paired_variant", "R2456H against wild type — the one pair")

    def show_variant_structures(self) -> None:
        """What the deposited variant structures can actually support."""
        from ..analysis.variant_structures import survey_variant_structures

        self._set_status("surveying the deposited variant structures…")
        try:
            survey = survey_variant_structures()
        except Exception as exc:
            self._set_status(f"variant survey failed: {exc}")
            return
        data = {
            "summary": survey.summary(),
            "coverage": survey.coverage(),
            "identical coordinate groups": survey.duplicate_groups(),
            "entries": [
                {"pdb": e.pdb, "variant": e.variant or "wild type",
                 "resolves its mutation": e.mutation_resolved,
                 "shows the mutant residue": e.shows_mutation,
                 "informative": e.informative,
                 "bottleneck_A": e.bottleneck_A,
                 "wetting_score": e.wetting_score,
                 "conductance_pS": e.conductance_pS,
                 "duplicates": list(e.duplicates)}
                for e in survey.entries],
        }
        self._show_result("variant_structures",
                          "Deposited variant structures — what they support",
                          data)
        self._set_status(survey.summary())
