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
    "fluctuations": (
        "A VALIDATION OF THE NETWORK, and mostly a verdict on the entry. An "
        "observed B-factor in a cryo-EM map absorbs local resolution, "
        "sharpening and the refinement's own restraints, so read the quality "
        "block first. Spearman is the number to use — the relationship is "
        "monotone but not linear. The contact-number CONTROL needs no network "
        "at all; if it wins, the agreement is burial rather than mechanism, "
        "and if it comes out NEGATIVE this entry's column rises with burial "
        "and is not reporting mobility."),
    "paralogue": (
        "THE ONLY GENERALITY CONTROL AVAILABLE. Read the two dome blocks "
        "together: measured naively the two proteins look very different, and "
        "that is because 6KG7 resolves all 38 transmembrane helices where a "
        "PIEZO1 entry resolves 22-26. Coverage-matched they are "
        "indistinguishable. n = 1 PIEZO2 structure, so this says the fold "
        "admits the mechanism, not that every PIEZO uses it."),
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
    "liu2025": (
        "A REPLICATION AUDIT of Liu et al. 2025, not a result about PIEZO1. "
        "Six panels reproduce from coordinates, seven have an ANALOGUE that "
        "is a different quantity, and eleven need patch clamp, a cryo-EM map "
        "or a molecular-dynamics trajectory this project does not hold. The "
        "curvature panel DISAGREES with the paper and says so: our sphere fit "
        "saturates on a nearly flat surface, giving 18 nm where they report "
        "117. Read the reason beside any panel not marked replicated."),
    "guo2017": (
        "A REPLICATION AUDIT, not a result about PIEZO1. It reports how much "
        "of the paper this project can reproduce from coordinates and, just "
        "as importantly, what it cannot: twelve panels need the cryo-EM map, "
        "the micrographs or two structures from other channel families, and "
        "three more have only an ANALOGUE that is a different quantity — a "
        "projection of a model is not a 2D class average, and a "
        "screened-Coulomb surface is not APBS. Read the reason beside any "
        "panel that is not marked replicated before putting it beside the "
        "original."),
    "homology": (
        "THE WIDEST GENERALITY CONTROL THERE IS, and the column to read is the "
        "OVERLAP, not the identity. Sequence identity to the invertebrate "
        "PIEZOs is ~0.30, inside Rost's twilight zone, where a percentage is "
        "barely distinguishable from what a shuffled sequence of the same "
        "composition gives — see docs/HOMOLOGY_SEARCH.md. The gating-mode "
        "overlap does not depend on it and carries a shuffled-correspondence "
        "control. Every partner here is a separate deposited structure with "
        "its own coverage, so the overlaps are not directly comparable with "
        "each other; and the family's two non-animal members (plant, amoeba) "
        "are ABSENT, because the only structural representation either has is "
        "an AlphaFold MONOMER and the elastic network needs three protomers."),
    "family": (
        "AN EXTERNAL PROJECT'S RESULTS, not this one's. Every statement here "
        "was produced by the piezo_genes census — a 194-genome sweep this "
        "application does not perform — and re-verified against its own source "
        "files when the resource was built. Nothing in this window was "
        "measured on the loaded structure. The four entries beside it "
        "(Constraint, Disease geography, Core and periphery, piezo3) are the "
        "ones that measure."),
    "constraint": (
        "THE PER-RESIDUE VALUES ARE THE CENSUS'S; the partition is ours. Two "
        "things follow. A mouse entry is read through the alignment map, so "
        "coverage below 1.0 is residues the map could not carry rather than "
        "residues nobody scored. And the blade gradient is BOUNDARY-DEPENDENT: "
        "the census's finding that the distal blade is more conserved than the "
        "proximal one holds on its chain-cut bands and REVERSES on the "
        "transmembrane units, because its proximal band is 77% inter-unit "
        "linker against the distal band's 29% and linker scores the same "
        "either side."),
    "disease": (
        "AN INDEPENDENT RE-TEST, not a reproduction, and it does not settle. "
        "PIEZO1 only where the census pooled two genes, and against gnomAD "
        "population missense rather than ClinVar benign labels. The enrichment "
        "reaches significance on the census's pore-module boundaries and not "
        "on ours, and the 120 residues the two disagree about (2057-2176) "
        "themselves carry six pathogenic positions. Read both rows."),
    "coreperiphery": (
        "A CORE-ONLY FIT IS DIRECTIONAL and can fail. It asks where the blades "
        "land given that the pore modules are superposed. The control that "
        "makes the answer readable: an AlphaFold MONOMER splays 7-9x from an "
        "experimental structure OF THE SAME PROTEIN, while two experimental "
        "structures of DIFFERENT paralogues splay ~1x. A large splay against a "
        "predicted model is therefore a statement about the model."),
    "piezo3": (
        "A PREDICTED MONOMER ASSEMBLED ON SOMEBODY ELSE'S TRIMER. 96% of the "
        "resulting departure from planarity is the template's arrangement, so "
        "the dome radius is not a measurement of piezo3. What the numbers can "
        "do is FAIL, and they did not: the protomer arranges into a closed "
        "trimer with an axis and a continuous lumen. That is a negative that "
        "survived, not a positive demonstrated — no current has ever been "
        "recorded from any piezo3, and human piezo3 has been a pseudogene "
        "since before the primate radiation."),
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
        # A result computed on a spliced model must say so *before* its own
        # caveat, because the reader's first question is what the numbers were
        # measured on. The structure name already carries the suffix, but a
        # name is a label and this is a sentence.
        caveat = CAVEATS.get(key, "")
        model = getattr(self, "full_length", None)
        if model is not None:
            caveat = (f"COMPUTED ON A PART-PREDICTED MODEL. "
                      f"{model.n_predicted_residues} residues here are "
                      f"AlphaFold, not experiment "
                      f"({model.confident_fraction:.0%} above pLDDT 70). "
                      f"Switch Completeness to 'Deposited only' in the Model "
                      f"panel for experiment alone.\n\n") + caveat
        dialog = ResultDialog(title, data, caveat, self,
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

    def show_paralogue(self) -> None:
        self._run_registry_analysis(
            "paralogue", "PIEZO1 against PIEZO2 — is this the fold?")

    def show_homology(self) -> None:
        self._run_registry_analysis(
            "homology", "Against every PIEZO in the catalogue")

    def show_fluctuations(self) -> None:
        self._run_registry_analysis(
            "fluctuations", "Predicted fluctuation against the B-factors")

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

    def show_family(self) -> None:
        """The imported census. Runs on no structure, and says so."""
        self._run_registry_analysis(
            "family", "The PIEZO family census — what was imported")

    def show_constraint(self) -> None:
        self._run_registry_analysis(
            "constraint", "Evolutionary constraint, on our own domains")

    def show_disease_geography(self) -> None:
        self._run_registry_analysis(
            "disease", "Where human disease sits — re-tested here")

    def show_core_periphery(self) -> None:
        self._run_registry_analysis(
            "coreperiphery", "Core and periphery — fitted on the pore module")

    def show_piezo3(self) -> None:
        self._run_registry_analysis(
            "piezo3", "piezo3 — the third vertebrate PIEZO")

    def show_paired_variant(self) -> None:
        self._run_registry_analysis(
            "paired_variant", "R2456H against wild type — the one pair")

    def show_liu2025(self) -> None:
        """Every panel of the paper the intermediate-open structure comes from."""
        self._run_registry_analysis(
            "liu2025", "Liu et al. 2025 — what reproduces, and what cannot")

    def show_guo2017(self) -> None:
        """Every panel of the paper the dome model comes from."""
        self._run_registry_analysis(
            "guo2017", "Guo & MacKinnon 2017 — what reproduces, and what cannot")

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
