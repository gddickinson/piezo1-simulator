"""What exploring means for the analyses this project measures itself.

Split from :mod:`piezo1.ui.exhibits` at the seam ``report_tags`` and
``report_family`` already use, and for the same reason: everything here is a
question asked of coordinates this project holds, and the companion file is
questions asked about somebody else's result. The two fail differently. A chart
here can be wrong because the measurement is wrong; a chart there can be wrong
because it is quoting a number that has since been superseded, and no amount of
re-running would tell you.

Each entry is data — no Qt, no computation — so the catalogue is readable and
diffable on its own, which is the half a reviewer has to check.

Everything that puts something on the 3-D view moved to
:mod:`piezo1.ui.exhibit_catalogue_structure` in Round 95 — those share a
failure mode of their own (a picture is over-read where a table is not) and
they are the half this application exists for, so they are read together.
"""

from __future__ import annotations

from .exhibits import Exhibit

__all__ = ["MEASURED"]

MEASURED: tuple[Exhibit, ...] = (

    # ----------------------------------------------------- ion permeation --
    Exhibit(
        analysis="permeation", kind="simulation",
        title="How far would the constriction have to open?",
        what="The entry's own measured radius profile with every slice "
             "narrower than R opened to R, through the closed-form "
             "series-resistor formula — the same independent check the "
             "drift-diffusion solver is held to. Two of the three controls are "
             "the inputs nothing has measured.",
        basis="modelled", simulation="pore_conductance",
        not_this="Not a state of the channel. Opening the constriction is a "
                 "construction made to show what the radius decides; this "
                 "entry's own bottleneck is marked, and it is where the pore "
                 "actually is."),
    Exhibit(
        analysis="permeation", kind="chart",
        title="This entry against the published conductance",
        what="What the solver returned, beside the closed-form check and the "
             "25-30 pS single-channel band from the literature.",
        basis="measured", plot="permeation_conductance",
        not_this="Landing in the published band would be tuning rather than "
                 "prediction: the in-pore diffusivity and the ion radius are "
                 "unmeasured and span 16-94 pS between them."),

    # ------------------------------------------- fluctuation vs B-factors --
    Exhibit(
        analysis="fluctuations", kind="chart",
        title="The network against a control that uses no network",
        what="Both correlations with the deposited B-factors, beside the "
             "contact-number control — a count of neighbours, no modes, no "
             "Hessian.",
        basis="measured", plot="fluctuation_controls",
        not_this="Not a validation of the modes unless the network beats the "
                 "control; where it does not, the agreement is burial."),

    # ---------------------------------------------------------- paralogue --
    Exhibit(
        analysis="paralogue", kind="chart",
        title="Gating-mode overlap, against its shuffled control",
        what="How much PIEZO1's gating mode overlaps PIEZO2's symmetric "
             "subspace, beside the same statistic on a shuffled "
             "correspondence.",
        basis="measured", plot="paralogue_overlap",
        not_this="n = 1 PIEZO2 structure. It says the fold admits the "
                 "mechanism, not that every PIEZO uses it."),
    Exhibit(
        analysis="paralogue", kind="chart",
        title="Why the naive dome comparison is a coverage artefact",
        what="Transmembrane helices resolved by each entry, and the shared "
             "set the coverage-matched dome is measured on.",
        basis="measured", plot="paralogue_coverage",
        not_this="The difference between the two dome rows is this bar chart, "
                 "not a difference between the two proteins."),

    # ------------------------------------------------------ whole family ---
    Exhibit(
        analysis="homology", kind="chart",
        title="The overlap is a range, not a number",
        what="Lowest and highest gating-mode overlap over every available "
             "entry pair per partner protein, with the shuffled control each "
             "has to beat.",
        basis="measured", plot="homology_range",
        not_this="Not comparable between proteins: each partner is a separate "
                 "deposited entry with its own coverage, and one pair reaching "
                 "0.98 is a cherry-pick by a factor of five."),

    # ------------------------------------------------------- interactions --
    Exhibit(
        analysis="interactions", kind="chart",
        title="What kind of contact, and how many",
        what="The inventory by kind on a log scale, because hydrogen bonds "
             "outnumber disulfides by three orders of magnitude and a linear "
             "axis would show one bar.",
        basis="measured", plot="interaction_counts",
        not_this="Not the same evidence bar to bar: a disulfide is resolved "
                 "density, a hydrogen bond is heavy-atom geometry with no "
                 "proton anywhere in the model."),

    # --------------------------------------------------------- labelling ---
    Exhibit(
        analysis="labelling", kind="simulation",
        title="Labelling in time, and the p³ that follows from it",
        what="Per-site labelling against exposure, and the fully-labelled "
             "channel fraction it implies for three sites. The protocol's own "
             "incubation time is marked.",
        basis="modelled", simulation="labelling_timecourse",
        not_this="The ceiling is set by the reactive fraction, not by "
                 "patience: at 90% reactive tags, 73% of channels is the most "
                 "any incubation can reach."),
    Exhibit(
        analysis="labelling", kind="chart",
        title="How many dyes a channel ends up carrying",
        what="The Binomial(3, p) occupancy at the protocol, beside the same "
             "distribution when a tenth of tags are unreactive.",
        basis="modelled", plot="labelling_histogram",
        not_this="A distribution over channels, not over cells: single "
                 "molecules, not an ensemble average."),
    Exhibit(
        analysis="labelling", kind="figure", figure="labelling.png",
        rebuild="python scripts/make_labelling_figure.py",
        title="The published labelling figure",
        what="Saturation, the p³ amplification and the 1:2:3-dye histogram, "
             "as the project's own figure.",
        basis="modelled",
        not_this="Kinetics imported from halotag_binding_sim; the linker "
                 "length and reactive fraction behind them are unverified."),

    # ---------------------------------------------------- HaloTag fusion ---
    Exhibit(
        analysis="fusion", kind="figure", figure="halotag_fold.png",
        rebuild="python scripts/make_model_figures.py --only halotag",
        title="The tag as its real fold, in one arbitrary orientation",
        what="6U32's coordinates placed at the modelled position, with the "
             "dye it carries.",
        basis="modelled",
        not_this="One draw of many: the position is modelled and the spin "
                 "about the linker is undetermined."),
    Exhibit(
        analysis="fusion", kind="chart",
        title="Where the tag can be, in nanometres",
        what="The accessible-volume envelope's range and median beside the "
             "modelled tag-to-pore-exit distance and the tag's own radius.",
        basis="modelled", plot="fusion_geometry",
        not_this="A region, not a pose. Nothing determines where in the "
                 "envelope the tag sits."),

    # ------------------------------------------------- full-length model ---
    Exhibit(
        analysis="hybrid", kind="figure", figure="hybrid_model.png",
        rebuild="python scripts/make_model_figures.py --only hybrid",
        title="Measured core, predicted blade",
        what="One protomer with the experimental core in grey and the grafted "
             "residues in AlphaFold's own confidence bands.",
        basis="modelled",
        not_this="The colours are the point: the two halves are not the same "
                 "kind of thing and must never be averaged across."),
    Exhibit(
        analysis="hybrid", kind="chart",
        title="How much is prediction, and how well does the seam fit",
        what="Experimental against predicted atoms, the fraction clearing "
             "pLDDT 70, and the seam fit beside the two models' global "
             "disagreement.",
        basis="modelled", plot="hybrid_split",
        not_this="A good local fit at the seam says nothing about the rest: "
                 "the two models differ by tens of Angstrom away from it."),

    # ------------------------------------------------ calcium nanodomain ---
    Exhibit(
        analysis="nanodomain", kind="simulation",
        title="Calcium against distance from the pore exit",
        what="The buffered-diffusion Green's function with the sensor "
             "occupancy it implies. The modelled tag distance, the sensor's "
             "Kd and the 90%-occupancy distance are all marked.",
        basis="modelled", simulation="calcium_profile",
        not_this="The distance is to the MODELLED tag position, not a "
                 "measured one — and moving the tag does not escape the "
                 "prediction, which is what the 90%-occupancy distance says: "
                 "it is an order of magnitude beyond the whole channel."),

    # ---------------------------------------------------------- ligands ----
    Exhibit(
        analysis="ligands", kind="simulation",
        title="Dose-response from the measured potencies",
        what="A one-site Hill curve through each modulator's measured EC50 or "
             "K_D, with a concentration you can move.",
        basis="curated", simulation="dose_response",
        not_this="The potency is measured; the curve is assumed. No "
                 "cooperativity was fitted, and a real PIEZO1 dose-response is "
                 "not a one-site binding isotherm."),
    Exhibit(
        analysis="ligands", kind="chart",
        title="Potency, and how much is known about the site",
        what="Each modulator's measured potency on a log axis, labelled with "
             "the evidence level its binding site carries.",
        basis="curated", plot="ligand_potency",
        not_this="No PIEZO structure with a bound modulator has ever been "
                 "deposited. Every site here is inferred."),

    # ------------------------------------------- the variant structures ----
    Exhibit(
        analysis="variant_structures", kind="chart",
        title="Every deposited variant entry, measured the same way",
        what="Bottleneck radius and wetting score per entry, with the ones "
             "sharing identical coordinates named in the note.",
        basis="measured", plot="variant_survey",
        not_this="Not a comparison of variants: all four entries are closed "
                 "and three of them are one set of coordinates."),
    Exhibit(
        analysis="paired_variant", kind="chart",
        title="R2456H's bottleneck, inside the wild-type range",
        what="The variant's bottleneck radius against the band three "
             "independent wild-type entries span.",
        basis="measured", plot="paired_bottleneck",
        not_this="Falling inside the band is not evidence the variant does "
                 "nothing — n = 1, and the band is three entries wide."),
    Exhibit(
        analysis="paired_variant", kind="chart",
        title="R2456H's wetting score, inside the wild-type range",
        what="The same comparison on the hydrophobic-gating score, which "
             "spreads more than twice as far among wild-type entries.",
        basis="measured", plot="paired_wetting",
        not_this="The spread is what makes a single pair interpretable, and "
                 "it is wide enough to hide a real difference."),
)
