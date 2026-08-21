"""What exploring means for results this project did not produce.

The companion to :mod:`piezo1.ui.exhibit_catalogue`, split along the seam that
matters most for a picture: everything here is either an **imported** census
result, a **frozen record**, or a **published figure** being reproduced. All
three look exactly like a measurement once they are drawn, which is why every
entry names its basis and every ``not_this`` line says whose number it is.

The census entries are the sharpest case. A trimer coloured by evolutionary
constraint is the most persuasive thing this application can put on a screen,
and not one of those values was measured here.

The structural displays for these results are in
:mod:`piezo1.ui.exhibit_catalogue_structure`, with the rest of them.
"""

from __future__ import annotations

from .exhibits import Exhibit

__all__ = ["IMPORTED"]

IMPORTED: tuple[Exhibit, ...] = (

    # ------------------------------------------------ the imported census --
    Exhibit(
        analysis="family", kind="chart",
        title="What the census actually swept",
        what="PIEZO proteins found per kingdom across 194 genomes, on a log "
             "axis, with the count that does not reconcile stated in the note "
             "rather than quietly dropped.",
        basis="imported", plot="census_kingdoms",
        not_this="Not a measurement on the loaded structure. Nothing in this "
                 "window touched coordinates at all."),
    Exhibit(
        analysis="family", kind="figure", figure="family_mechanics.png",
        rebuild="python scripts/make_family_figures.py",
        title="Does the elastic network explain what evolution kept?",
        what="Every mechanical feature against the census's constraint track, "
             "raw and with burial held fixed, with burial's own correlation "
             "drawn as the line each has to beat.",
        basis="imported",
        not_this="Burial is the reason most of the raw correlations exist. "
                 "Read the partialled column, not the raw one."),

    # ------------------------------------------------------- constraint ----
    Exhibit(
        analysis="constraint", kind="chart",
        title="Constraint domain by domain, on our boundaries",
        what="The most and least constrained domains against the "
             "whole-protein mean, with the pore machinery at one end and the "
             "blades at the other.",
        basis="imported", plot="constraint_domains",
        not_this="The values are the census's; only the partition is ours. A "
                 "domain scoring low may be a domain the alignment could not "
                 "carry, which is what the coverage figure is for."),
    Exhibit(
        analysis="constraint", kind="figure",
        figure="family_constraint_by_domain.png",
        rebuild="python scripts/make_family_figures.py",
        title="The same result as the census's own picture",
        what="Constraint per domain on this project's boundaries, which sit "
             "up to 141 residues from the census's own.",
        basis="imported",
        not_this="Agreement here is a replication of an ordering, not of a "
                 "number: the two partitions disagree about which residues "
                 "belong to which domain."),
    Exhibit(
        analysis="constraint", kind="figure", figure="family_blade_gradient.png",
        rebuild="python scripts/make_family_figures.py",
        title="Where the blade gradient reverses, and why",
        what="The census's distal-beats-proximal result on its own chain-cut "
             "bands, the reversal on the transmembrane units, and what each "
             "band is made of.",
        basis="imported",
        not_this="Not a refutation of the census's arithmetic. The bands "
                 "differ in how much inter-unit linker they contain, and "
                 "linker scores the same either side."),

    # ---------------------------------------------------- disease geography
    Exhibit(
        analysis="disease", kind="chart",
        title="The enrichment under both partitions",
        what="Odds ratio for pathogenic missense in the pore module, on the "
             "census's boundaries and on ours, beside the census's own "
             "published value, against the line of no enrichment.",
        basis="imported", plot="disease_odds",
        not_this="Two bars of one result, not two results. The comparator "
                 "differs too — gnomAD population missense here, ClinVar "
                 "benign labels there."),
    Exhibit(
        analysis="disease", kind="chart",
        title="The 120 residues the two partitions disagree about",
        what="How far each element's boundaries move between the census's "
             "partition and ours, in residues.",
        basis="imported", plot="disease_boundaries",
        not_this="Not a bookkeeping detail: the disputed band carries six "
                 "pathogenic positions, which is what moves the odds ratio."),

    # ------------------------------------------------- core and periphery --
    Exhibit(
        analysis="coreperiphery", kind="chart",
        title="Splay, against the control that reinterprets it",
        what="This pair's core and periphery fit and the ratio between them, "
             "beside the band experimental cross-paralogue pairs occupy and "
             "the band an AlphaFold monomer of the same protein reaches.",
        basis="measured", plot="splay_ratio",
        not_this="A large splay against a predicted model is a statement "
                 "about the model, not about paralogue divergence."),
    Exhibit(
        analysis="coreperiphery", kind="figure", figure="family_splay.png",
        rebuild="python scripts/make_family_figures.py",
        title="Every pair measured the same way",
        what="The splay ratios with the prediction-versus-experiment control "
             "that reinterprets the census's structural finding.",
        basis="measured",
        not_this="The control is the finding here. Without it, 7-9x reads as "
                 "divergence rather than as prediction error."),
    Exhibit(
        analysis="coreperiphery", kind="figure",
        figure="gating_morph_small.gif",
        rebuild="python scripts/make_animations.py --only readme",
        title="What core-conserved and periphery-free looks like",
        what="PIEZO1's own curved-to-flat transition, 7WLT to 7WLU — the pair "
             "that splays 19x by this same measurement.",
        basis="modelled",
        not_this="An interpolation between two deposited endpoints, not a "
                 "trajectory: nothing here says how the channel gets from one "
                 "to the other."),

    # --------------------------------------------------------- piezo3 ------
    Exhibit(
        analysis="piezo3", kind="chart",
        title="The fourteen pathogenic pore positions in piezo3",
        what="How many of the positions human disease strikes in PIEZO1 and "
             "PIEZO2 carry the identical residue in the zebrafish paralogue.",
        basis="imported", plot="piezo3_positions",
        not_this="Identity at a position is not function. No current has ever "
                 "been recorded from any piezo3."),
    Exhibit(
        analysis="piezo3", kind="chart",
        title="How much of the assembled trimer is the template's",
        what="The borrowed fraction of the assembly's departure from "
             "planarity, with the template's identity and the clash count "
             "beside it.",
        basis="modelled", plot="piezo3_template",
        not_this="At this borrowed fraction a dome radius measured on the "
                 "assembly is mostly a measurement of the template."),

    # ----------------------------------------------- the prediction record -
    Exhibit(
        analysis="prediction_record", kind="chart",
        title="Five pre-registered tests, against the line of no effect",
        what="Each round's Cliff's delta and p, drawn against zero — the "
             "central claim's whole record in one picture.",
        basis="record", plot="record_effects",
        not_this="Not a result awaiting more data. The feasibility analysis "
                 "found no dataset that could exist would settle it."),
    Exhibit(
        analysis="prediction_record", kind="figure", figure="record_nulls.png",
        rebuild="python scripts/make_record_figure.py",
        title="The forest plot as published",
        what="Every pre-registered test with its interval, against the line "
             "of no effect.",
        basis="record",
        not_this="An interval spanning zero is a null however small the p "
                 "beside it is — which is the clause that caught Round 41."),
    Exhibit(
        analysis="prediction_record", kind="figure",
        figure="record_data_limit.png",
        rebuild="python scripts/make_record_figure.py",
        title="What could be reached, against what would be required",
        what="Round 47's data limit: the variants that could ever exist "
             "beside the number the observed effect would need.",
        basis="record",
        not_this="Not 'we need more data'. The reachable ceiling is below the "
                 "requirement, which is a different statement."),

    # ------------------------------------------------ Guo & MacKinnon 2017 -
    Exhibit(
        analysis="guo2017", kind="simulation",
        title="Figure 7c as arithmetic: flattening the idealised dome",
        what="The dome flattened at constant membrane area, with the "
             "projected area it releases, the bending energy it gives up and "
             "the free energy at the tension you choose.",
        basis="published", simulation="dome_flattening",
        not_this="The paper's idealisation, not a measurement: two lengths "
                 "and closed-form spherical-cap geometry, with no structure "
                 "in it anywhere."),
    Exhibit(
        analysis="guo2017", kind="simulation",
        title="Figure 7d: tension against open probability",
        what="The two-state Boltzmann with the area change and the intrinsic "
             "bias as sliders, and the measured half-activation band marked.",
        basis="published", simulation="dome_activation",
        not_this="Passing through the measured T50 is what the parameters "
                 "were fitted to do. It is a consistency check, not a "
                 "prediction."),
    Exhibit(
        analysis="guo2017", kind="chart",
        title="How much of the paper reproduces",
        what="Panels that reproduce from coordinates, panels with only an "
             "analogue, and panels needing data this project does not hold.",
        basis="published", plot="paper_coverage",
        not_this="An analogue is a different quantity, not a weaker "
                 "reproduction — a projection of a model is not a 2-D class "
                 "average."),
    Exhibit(
        analysis="guo2017", kind="figure",
        figure="guo2017/figure_7c_flattening.png",
        rebuild="python scripts/make_guo2017_figures.py",
        title="Figure 7c, regenerated",
        what="The flattening series as the paper draws it, from this "
             "project's own arithmetic.",
        basis="published",
        not_this="Every number here follows from the paper's two lengths. "
                 "Agreement is arithmetic, not independent confirmation."),
    Exhibit(
        analysis="guo2017", kind="figure",
        figure="guo2017/figure_6b_pore.png",
        rebuild="python scripts/make_guo2017_figures.py",
        title="Figure 6b: the pore profile",
        what="The radius along the conduction axis, measured here on the "
             "paper's own entry.",
        basis="published",
        not_this="Ours is an Apollonius clearance profile, theirs is HOLE. "
                 "The offset between the two is pinned in both directions."),
    Exhibit(
        analysis="guo2017", kind="figure",
        figure="guo2017/figure_2ab_projection.png",
        rebuild="python scripts/make_guo2017_figures.py",
        title="Figure 2a,b: the analogue, and what it is missing",
        what="Simulated projections of the model at the paper's own pixel "
             "size.",
        basis="published",
        not_this="Not a 2-D class average. No CTF, no solvent and no "
                 "detergent micelle — which is most of what the published "
                 "side view shows."),

    # ------------------------------------------------------- Liu et al 2025
    Exhibit(
        analysis="liu2025", kind="chart",
        title="How much of the paper reproduces",
        what="Six panels reproduce, seven have an analogue, eleven need "
             "electrophysiology, a map or a trajectory this project does not "
             "hold.",
        basis="published", plot="paper_coverage",
        not_this="The refusals are the informative half. A tool that quietly "
                 "covered the tractable panels would leave a reader assuming "
                 "the rest."),
    Exhibit(
        analysis="liu2025", kind="figure", figure="liu2025_5e_iv.png",
        rebuild="python scripts/make_liu2025_figures.py",
        title="Figure 5E: the current-voltage relation",
        what="The drift-diffusion solver swept over their four voltages, with "
             "the slope read off the way they read it, against their 20 pS.",
        basis="published",
        not_this="A 2x overestimate, consistent with the 1.4x this model "
                 "already carries. Not a reproduction of their number."),
    Exhibit(
        analysis="liu2025", kind="figure", figure="liu2025_2d_pore.png",
        rebuild="python scripts/make_liu2025_figures.py",
        title="Figure 2D: the pore through all four states",
        what="Radius along the axis for the closed, intermediate, open and "
             "flattened entries, measured here by one route.",
        basis="published",
        not_this="A radius is not a verdict. The intermediate-open entry is "
                 "refused by an axial model on a neck their paper says the "
                 "lateral portals bypass."),
    Exhibit(
        analysis="liu2025", kind="figure", figure="liu2025_distances.png",
        rebuild="python scripts/make_liu2025_figures.py",
        title="Seven published distances, measured here",
        what="The gate dilation, the cap-gate loops, the spring linker and "
             "the axis shortening, ours against theirs.",
        basis="published",
        not_this="Side-chain distances use side-chain atoms. Where the atoms "
                 "are unresolved the answer is None, never a backbone "
                 "distance standing in."),
    Exhibit(
        analysis="liu2025", kind="figure", figure="liu2025_6_curvature.png",
        rebuild="python scripts/make_liu2025_figures.py",
        title="Where we disagree with the paper",
        what="The curvature comparison on a log axis: our sphere fit "
             "saturates on a nearly flat surface where they report 117 nm.",
        basis="published",
        not_this="Not evidence the paper is wrong. It is evidence a sphere "
                 "fit stops meaning anything as the surface flattens."),
)
