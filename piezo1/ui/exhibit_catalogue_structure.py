"""Every exhibit that puts something on the 3-D view.

Split from the other two catalogues along the seam this application is *for*:
the point of the project is that PIEZO1's shape is its mechanism, so a result
that can only be read as a table is doing half its job. These are the
structural displays — load the entry a result is about, draw a second one
beside it, superpose one on another, show one part of the assembly up close,
mark the residues a finding is made of, recolour by the quantity being
discussed, or build the gating motion.

Each names an action in :mod:`piezo1.ui.model_actions`, which presses the
control a user would press — the Model panel's entry selector, the Overlay
panel's button, the View menu's component group. Nothing here draws anything
itself.

The ``not_this`` lines carry more weight in this file than anywhere else. A
superposition is the most persuasive picture this application can produce and
the easiest to over-read: two structures on top of each other look like a
measurement of how similar they are, when what they are is a fit chosen from
several possible fits.
"""

from __future__ import annotations

from .exhibits import Exhibit

__all__ = ["STRUCTURAL"]


def _e(analysis: str, action: str, title: str, what: str, basis: str,
       not_this: str) -> Exhibit:
    return Exhibit(analysis=analysis, kind="model", action=action, title=title,
                   what=what, basis=basis, not_this=not_this)


STRUCTURAL: tuple[Exhibit, ...] = (

    # ----------------------------------------------------- ion permeation --
    _e("permeation", "pore_surface",
       "Draw the pore as the spheres it was measured with",
       "At each height, the largest probe sphere that fits without touching an "
       "atom — the geometry the conductance was computed over, narrowest slice "
       "marked.", "measured",
       "A probe sphere is the space left over, not the pore wall — and a "
       "radius does not settle whether the pore conducts."),
    _e("permeation", "component_pore",
       "Show the pore module on its own",
       "Outer helix, cap, spring linker, inner helix and CTD with all four "
       "gates — about a quarter of the atoms, and the part the current has to "
       "cross.", "curated",
       "It HIDES rather than subsets: every analysis still runs on the whole "
       "assembly, so the numbers do not change when you switch."),
    _e("permeation", "highlight_gate",
       "Mark the hydrophobic gate",
       "The three curated gate residues, marked on all three protomers in this "
       "entry's own numbering.", "curated",
       "The gate is 2.4-4.7 A across the catalogue — at or above the water "
       "radius everywhere. It is chemistry that shuts it, not width."),
    _e("permeation", "highlight_selectivity",
       "Mark the selectivity glutamates",
       "The curated acidic residues said to set ion selectivity.", "curated",
       "Three of the four are measured NOT to reach the lumen — which is what "
       "Coste et al. concluded about E2117 from function alone."),
    _e("permeation", "load_open_like",
       "Load the one entry that conducts (11ZC)",
       "The open-like structure: on the axial route it is the only entry in "
       "the catalogue the model passes current through.", "measured",
       "Deposited without side chains, so anything needing them falls back to "
       "C-alpha — the pore charge map says so where it does."),
    _e("permeation", "ion_flux",
       "Run ions through at the computed rate",
       "The ion stream, timed from this structure's own current. A shut pore "
       "draws nothing and says why.", "modelled",
       "The stream runs about a millionfold slow — a real channel passes "
       "~10^7 ions/s. Read the HUD's factor, not the pace."),

    # ------------------------------------------- fluctuation vs B-factors --
    _e("fluctuations", "colour_fluctuation",
       "Colour by the network's predicted fluctuation",
       "The array this comparison scores, painted on the structure and read "
       "from the mode set rather than recomputed.", "measured",
       "This is the prediction. The observation is the B-factor column beside "
       "it, which is the other half of the comparison."),
    _e("fluctuations", "colour_bfactor",
       "Colour by the deposited B-factor",
       "The observed column — the thing the network is being checked against, "
       "on the same structure.", "measured",
       "A cryo-EM B-factor absorbs local resolution, sharpening and the "
       "refinement's own restraints; a grouped or predicted column cannot "
       "answer the question at all, which is why the result gates it first."),

    # ---------------------------------------------------------- paralogue --
    _e("paralogue", "overlay_core_piezo2",
       "Superpose PIEZO2 on the pore module",
       "6KG7 fitted on the outer helix, cap, inner helix and CTD alone, with "
       "the blades then measured rather than fitted.", "measured",
       "Not a fit of the whole protein: by residue number the same pair gives "
       "a confident 47.9 A, because the two proteins do not share a numbering "
       "— this one corresponds them through a real alignment."),
    _e("paralogue", "companion_piezo2",
       "Draw PIEZO2 beside it",
       "6KG7 in its own colour in the same frame, unsuperposed.", "measured",
       "Two structures in one frame is not a comparison — it is two "
       "structures in one frame. The analyses run on the primary only."),

    # ------------------------------------------------------- whole family --
    _e("homology", "overlay_core_worm",
       "Superpose PEZO-1 on the pore module",
       "The nematode channel, 9ZIS, fitted on the pore module through a global "
       "alignment.", "measured",
       "Sequence identity here is inside Rost's twilight zone, so the "
       "alignment behind this correspondence is itself the weakest link — "
       "which is why the family result reports mode overlap and not "
       "identity."),
    _e("homology", "component_pore",
       "Show the part being compared",
       "The pore module of the loaded entry — the region every cross-species "
       "fit in this result is anchored on.", "curated",
       "Curated boundaries for PIEZO1. The other family members have no "
       "curated domains at all, which is why the comparison is anchored on "
       "this entry's."),

    # ------------------------------------------------------- interactions --
    _e("interactions", "contacts",
       "Draw the contacts between the atoms they were found between",
       "One cylinder per contact, coloured by kind, from this same analysis "
       "rather than a second implementation.", "measured",
       "Contacts of this entry in this state. A closed-state entry cannot "
       "show the open-state salt bridges."),
    _e("interactions", "component_pore",
       "Show the pore module, where the gating contacts are",
       "The same inventory drawn on a quarter of the atoms, so the pore's own "
       "contacts are not buried in the blades'.", "curated",
       "Hiding atoms changes nothing about what was detected: the counts in "
       "the table are for the whole assembly either way."),

    # --------------------------------------------------------- the tag ----
    _e("labelling", "halotag",
       "Draw the tags that are being labelled",
       "A HaloTag at each of the three cytosolic C-termini, which is what the "
       "Binomial(3, p) is over.", "modelled",
       "There is no structure of the fusion. The three sites are a model, and "
       "so is the distance any dye would sit at."),
    _e("fusion", "halotag",
       "Draw the fusion on the structure",
       "The tag, its linker seams and the envelope of positions the linker "
       "admits.", "modelled",
       "There is no structure of the PIEZO1-HaloTag fusion. Everything drawn "
       "here is construction."),
    _e("fusion", "component_vestibule",
       "Show where it is attached",
       "The cytoplasmic vestibule and the C-terminal domain the tag hangs "
       "off.", "curated",
       "The anchor residue is real; where the tag goes from there is not."),

    # ------------------------------------------------- full-length model ---
    _e("hybrid", "full_length",
       "Draw the full-length model",
       "The graft on screen with its seam marked and the predicted part banded "
       "by pLDDT.", "modelled",
       "A complete-looking trimer is the hazard here. The amber banner is what "
       "tells you part of it is prediction."),
    _e("hybrid", "colour_plddt",
       "Colour by AlphaFold's own confidence",
       "The prediction's pLDDT on the model, on the fixed band scale.",
       "modelled",
       "pLDDT is the prediction's confidence in a LOCAL structure. It says "
       "nothing about whether the blade is in the right place relative to the "
       "core, which is what the PAE analysis found it does not settle."),

    # ------------------------------------------------ calcium nanodomain ---
    _e("nanodomain", "nanodomain",
       "Draw the field around the pore exit",
       "Shells at decade concentrations from the measured cytosolic mouth. A "
       "shut structure draws nothing and says why.", "modelled",
       "The two surfaces that carry the result — 90% occupancy and the K_D — "
       "are reported rather than drawn, because the whole channel sits inside "
       "both."),
    _e("nanodomain", "component_vestibule",
       "Show the mouth the calcium comes out of",
       "The cytoplasmic vestibule, with the constrictions the exit is measured "
       "from.", "curated",
       "The source is a point on the axis, not a residue: the model is "
       "spherically symmetric and does not know the protein is there."),

    # ---------------------------------------------------------- ligands ----
    _e("ligands", "pockets",
       "Draw the cavities a ligand could occupy",
       "The top-ranked pockets as the alpha spheres they were detected with.",
       "measured",
       "A cavity is geometry, not a binding site — and ligands were excluded "
       "before detection, so a drawn pocket may sit on a resolved lipid."),
    _e("ligands", "highlight_yoda1",
       "Mark the proposed Yoda1 pocket",
       "The three residues whose substitution changes the Yoda1 response, in "
       "this entry's numbering.", "curated",
       "INFERRED from docking and mutagenesis. No PIEZO structure with a bound "
       "modulator has ever been deposited, so nothing here was observed."),
    _e("ligands", "highlight_pip2",
       "Mark the PIP2-binding cluster",
       "The polybasic cluster the lipid is proposed to sit on.", "curated",
       "A cluster of lysines is not a bound lipid. The evidence level is on "
       "the result beside it."),
    # ------------------------------------------- variants on the structure --
    _e("variant_structures", "load_variant",
       "Load the one entry that resolves its own mutation (8YFG)",
       "R2456H, the only deposited variant entry whose mutated residue is "
       "modelled.", "measured",
       "One entry. Three of the other four share a single set of coordinates, "
       "and all of them are closed."),
    _e("variant_structures", "highlight_r2456",
       "Mark R2456 on the structure",
       "The position, marked on all three protomers in this entry's own "
       "numbering.", "curated",
       "Marking a position is not showing a mutation: unless 8YFG is loaded, "
       "what is under the marker is the wild-type residue."),
    _e("paired_variant", "overlay_variant",
       "Superpose the variant on the pore module",
       "8YFG fitted on the pore module of the entry on screen, so the "
       "comparison is made where the difference is claimed to be.", "measured",
       "n = 1. The wild-type entries differ among themselves by more than the "
       "variant differs from them, which is the result."),
    _e("paired_variant", "highlight_r2456",
       "Mark the substituted position",
       "R2456, in this entry's numbering.", "curated",
       "The variant falls inside the wild-type range on both metrics — a "
       "marker says where to look, not that anything is different there."),
    _e("paired_variant", "load_variant",
       "Load the variant entry itself (8YFG)",
       "The deposited structure of R2456H.", "measured",
       "Its own numbering is repaired to 0.999 rather than 1.000, and that is "
       "correct: R2456H is a real residue change a numbering fix must not "
       "absorb."),

    # ------------------------------------------------- the imported census --
    _e("family", "colour_constraint",
       "Colour the model by what the census measured",
       "Per-residue evolutionary constraint over 174 orthologues, painted on "
       "the structure — the census's central claim as a picture.", "imported",
       "NOT MEASURED HERE. An unscored residue is grey, and grey is where "
       "coverage ran out, not where constraint is low."),
    _e("constraint", "colour_constraint",
       "Colour the model by evolutionary constraint",
       "The per-residue track on a scale fixed at 0-1, so two entries stay "
       "comparable.", "imported",
       "Whose numbers these are: the census's, not this project's. An "
       "unscored residue must not be read as an unconstrained one."),
    _e("constraint", "component_anchor",
       "Show the most constrained domain",
       "The anchor, which comes out top of this project's own partition at "
       "0.83 against a whole-protein 0.65.", "curated",
       "The ordering is the result; the boundary is ours and sits up to 141 "
       "residues from the census's own."),
    _e("constraint", "component_blade",
       "Show the least constrained part",
       "THU1-THU9, the blades, where the census measures constraint lowest.",
       "curated",
       "Coverage is also worst here, and low coverage and low constraint are "
       "different things — which is why the unscored residues are grey rather "
       "than dark."),

    # ------------------------------------------------------ where disease is
    _e("disease", "component_pore",
       "Show the region the enrichment is about",
       "The pore module on this project's boundaries — the region pathogenic "
       "missense is tested against.", "curated",
       "On the census's boundaries the same region is 120 residues larger, "
       "and those residues carry six pathogenic positions. Which boundary is "
       "drawn decides the answer."),
    _e("disease", "highlight_pathogenic",
       "Mark the pathogenic pore positions",
       "The nine PIEZO1 positions of the census's fourteen — the other five "
       "are PIEZO2's and have no place on this structure — converted into "
       "this entry's numbering.", "imported",
       "Variants are found where people look. This is a map of reported "
       "disease, not of where disease can occur."),
    _e("disease", "colour_constraint",
       "Colour by the score that classifies them",
       "The constraint track the AUC 0.82 classifier is built on, on the "
       "structure.", "imported",
       "The negative set is population variation, not a clinical judgement of "
       "benignity — the classifier separates constrained from unconstrained, "
       "not pathogenic from benign."),

    # ------------------------------------------------- core and periphery --
    _e("coreperiphery", "overlay_core_piezo2",
       "Superpose PIEZO2 on the pore module",
       "The comparison this result is made of, drawn: cores on top of each "
       "other, blades left where they fall.", "measured",
       "A core-only fit is DIRECTIONAL. It answers where the blades land "
       "given that the pores are superposed, and a pair whose cores do not "
       "fit gets no ratio at all."),
    _e("coreperiphery", "overlay_core_flat",
       "Superpose PIEZO1's own flattened state",
       "7WLT against 7WLU by the same route — the gating transition measured "
       "as core-conserved and periphery-free, at 19x.", "measured",
       "The two are the same protein in two states, so this is what the "
       "channel's own motion looks like by this measure — not a divergence "
       "between proteins."),
    _e("coreperiphery", "highlight_equivalent",
       "Mark the two equivalent positions",
       "PIEZO1 R2456 and R2488 — the positions the census pairs with PIEZO2 "
       "disease residues, in this entry's numbering.", "imported",
       "The evidence is the REGISTER, not the proximity: after a pore-module "
       "fit the median aligned core pair is 2.5 A apart, so proximity is what "
       "every pair gives."),
    _e("coreperiphery", "component_pore",
       "Show what the fit was made on",
       "The pore module — the only part of the two structures that entered "
       "the superposition.", "curated",
       "Everything outside this selection is measurement rather than fit, "
       "which is the whole construction."),

    # ------------------------------------------------------------ piezo3 ---
    _e("piezo3", "load_piezo3_model",
       "Load the only piezo3 coordinates that exist",
       "The AlphaFold model of the zebrafish protein — human piezo3 has been "
       "a pseudogene since before the primate radiation.", "modelled",
       "A predicted MONOMER. The trimer in the result was assembled on "
       "somebody else's template, and 96% of its departure from planarity is "
       "that template's."),
    _e("piezo3", "highlight_pathogenic",
       "Mark the positions it keeps",
       "The nine PIEZO1 pathogenic pore positions on the entry you have "
       "loaded, for comparison with the piezo3 residues listed in the "
       "result.", "imported",
       "Identity at a position is not function. No current has ever been "
       "recorded from any piezo3."),

    # -------------------------------------------- Guo & MacKinnon 2017 -----
    _e("guo2017", "load_guo_entry",
       "Load the paper's own entry (6B3R)",
       "Every number in this replication was measured on it.", "measured",
       "Mouse numbering, like most of the catalogue — a residue number copied "
       "from here into a human-numbered paper is wrong by up to 26."),
    _e("guo2017", "dome_surface",
       "Draw the dome the mechanism rests on",
       "The fitted sphere cap out to the footprint radius and its own flat "
       "projection; the gap between them is the excess area.", "measured",
       "The far-field footprint is deliberately not drawn: linear theory "
       "overestimates it 3.65x at this contact slope."),
    _e("guo2017", "planar_membrane",
       "Draw Figure 4a's two planes",
       "The planar membrane fitted to one protomer, and the same construction "
       "on the trimer — the contrast is the paper's point.", "measured",
       "Every point set has a best-fit plane. Read the residual, not the "
       "lines."),
    _e("guo2017", "micelle",
       "Draw Figure 4b's envelope, as a model",
       "The surface a fixed distance outside the hydrophobic transmembrane "
       "belt.", "modelled",
       "Not the observed micelle density. The thickness is a parameter and "
       "carries no information; only the curvature is a measurement."),
    _e("guo2017", "morph",
       "Build the flattening the figure idealises",
       "The curved-to-flat interpolation between two deposited endpoints, "
       "which is what Figure 7c draws as a smooth sweep.", "modelled",
       "An INTERPOLATION between two structures, not a trajectory: nothing "
       "here says how the channel gets from one to the other."),

    # -------------------------------------------------- Liu et al. 2025 ----
    _e("liu2025", "load_open",
       "Load the intermediate-open state (8IXO)",
       "The structure this paper contributes, and the one an axial model "
       "refuses on a neck their Figure 5 says is bypassed.", "measured",
       "Intermediate-open is their reading of the density; an axial pore "
       "profile of it does not conduct in this model, which is the "
       "disagreement rather than a correction."),
    _e("liu2025", "load_flat",
       "Load the flattened state (7WLU)",
       "The fourth of their four states, all of which are deposited and "
       "catalogued here.", "measured",
       "Flattened is a state of the dome, not an open channel — the wetting "
       "verdict is computed separately and reported beside it."),
    _e("liu2025", "component_cap",
       "Show the cap and its gates",
       "The lateral cap gates their Figure 2 measures opening, on the entry "
       "you have loaded.", "curated",
       "The three cap gates are where their model has ions entering. This "
       "project's axial profile does not contain that route at all."),
    _e("liu2025", "component_md",
       "Show the construct they simulated",
       "Pore module, beam and lateral plug gate — the system their molecular "
       "dynamics ran on.", "curated",
       "This project prepares that construct and does not run it. No number "
       "here comes from a trajectory."),
    _e("liu2025", "morph",
       "Build the transition between two of their states",
       "The interpolation between the curved and flattened endpoints.",
       "modelled",
       "Their paper argues for an intermediate on the way; an interpolation "
       "passes through one by construction and is not evidence of it."),
    _e("liu2025", "pore_surface",
       "Draw the pore of the state you have loaded",
       "The probe spheres, so the constriction their Figure 2D plots acquires "
       "a location on the structure.", "measured",
       "The axial route only. Their Figure 5 says the current does not go "
       "down the axis at all."),
    # ------------------------------------------- the rest of the displays --
    _e("permeation", "component_gate",
       "Show the transmembrane gate up close",
       "The three curated gate residues and the helices they sit on.",
       "curated",
       "A gate is a chemical constriction here, not a physical one: this one "
       "is at or above the water radius in every entry in the catalogue."),
    _e("permeation", "colour_hydrophobicity",
       "Colour by hydrophobicity",
       "The Kyte-Doolittle scale on a FIXED range, which is what the wetting "
       "verdict is computed from.", "measured",
       "Fixed at +/-4.5 deliberately: an auto-ranged map paints a uniformly "
       "polar loop in full orange, and two structures coloured that way "
       "cannot be compared with each other."),
    _e("permeation", "colour_electrostatics",
       "Colour by surface electrostatics",
       "Screened Coulomb from formal charges — the field an ion crossing this "
       "pore would feel.", "modelled",
       "NOT APBS: no dielectric boundary, no ion-exclusion layer, no partial "
       "charges. All three omissions push the same way, so the magnitude is a "
       "lower bound. Read the sign and the pattern."),
    _e("permeation", "companion_open",
       "Draw the conducting entry beside this one",
       "11ZC alongside, in its own colour and the same frame.", "measured",
       "Two structures in one frame is not a superposition and not a "
       "comparison — it is two structures in one frame."),
    _e("permeation", "component_whole",
       "Put the whole assembly back",
       "Undo a component selection and draw the trimer again.", "curated",
       "Nothing was ever removed from the analyses: hiding is a display "
       "choice, and the numbers were the same throughout."),
    _e("nanodomain", "highlight_ctd",
       "Mark the constrictions the ions leave through",
       "The curated cytoplasmic vestibule constrictions, in this entry's "
       "numbering.", "curated",
       "The calcium source is a point on the axis, not these residues: the "
       "Green's function is spherically symmetric and knows nothing about "
       "them."),
    _e("paralogue", "load_piezo2",
       "Load PIEZO2 itself (6KG7)",
       "The paralogue entry, which resolves all 38 transmembrane helices "
       "where a PIEZO1 entry resolves 22-26.", "measured",
       "Its numbering is PIEZO2's. This project has no curated annotation for "
       "it, so the Annotation panel comes up empty and says why."),
    _e("constraint", "allosteric_path",
       "Draw the route the mechanical result is about",
       "The cheapest blade-to-gate path through the correlation graph — the "
       "coupling that predicts constraint at partial rho = 0.29.", "measured",
       "A drawn line reads as unique and is not: the best route avoiding it "
       "costs 1.001x on 8YEZ, and the status line says so."),
    _e("disease", "colour_domain",
       "Colour by the partition itself",
       "This project's domain boundaries on the structure — the partition the "
       "enrichment is computed under.", "curated",
       "The census's boundaries are not these. Where they disagree is 120 "
       "residues carrying six pathogenic positions, which is what moves the "
       "odds ratio."),
    _e("coreperiphery", "overlay_flat",
       "Superpose the same pair on everything shared",
       "7WLU fitted on all shared C-alphas rather than on the pore module — "
       "the contrast that shows what a core-only fit is doing.", "measured",
       "Neither fit is the right one. They answer different questions, and "
       "the ratio between the two RMSDs is what the result reports."),
    _e("liu2025", "load_curved",
       "Load the closed state (7WLT)",
       "The first of their four states.", "measured",
       "Closed is a gating state, not a dome shape: this entry is also the "
       "curved one, and the two claims are measured separately here."),
    _e("liu2025", "load_intermediate",
       "Load the intermediate state (8IXN)",
       "The state between closed and intermediate-open in their scheme.",
       "measured",
       "A structure is one point, and a scheme is an ordering of points. "
       "Nothing here shows the channel passing through it."),
    _e("liu2025", "companion_flat",
       "Draw the flattened state beside this one",
       "7WLU alongside, so the dome difference is visible without a fit.",
       "measured",
       "Side by side is not superposed: what looks like a shift may be how "
       "each entry happens to sit in its own frame, which is why the "
       "canonical framing is applied to both."),
    _e("guo2017", "component_beam",
       "Show the beam, the lever the mechanism needs",
       "The beam and the lateral plug — the long helix Figure 7b colours "
       "separately from the cross-helices.", "curated",
       "The paper gives no residue range for the cross-helices anywhere; "
       "this project finds them from coordinates, and the beam is excluded "
       "from that set by design."),
)
