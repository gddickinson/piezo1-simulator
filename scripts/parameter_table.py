"""The parameter table itself — every number a calculation depends on.

Split out of ``build_parameters.py`` to keep both files under the project's
500-line limit. This module is **data**: the list and nothing else, so the whole
parameter set can be read and diffed without the validation machinery in the
way. The provenance gate that decides whether it may be written lives in
``build_parameters.py``.

Each entry needs a unit, bounds, a ``kind``, a ``description`` and a
``citation`` — the last resolving in ``references.json`` or naming one of the
sentinels, which then obliges a ``source_note`` saying why.
"""

from __future__ import annotations

__all__ = ["P"]

from parameter_table_figures import FIGURE_PARAMETERS
from parameter_table_tags import TAG_PARAMETERS

#: kind: physical | empirical | method | convention
P = [
    # ---------------------------------------------------------- membrane ---
    dict(key="membrane.kappa", name="Bending modulus", value=20.0, unit="k_BT",
         minimum=5.0, maximum=60.0, kind="physical", category="Membrane mechanics",
         citation="haselwandter2018",
         source_note="20-25 k_BT for a typical bilayer; 20 used for PIEZO1",
         description="Helfrich bending rigidity of the bilayer."),
    dict(key="membrane.tension", name="Resting membrane tension", value=0.1,
         unit="k_BT/nm^2", minimum=0.0, maximum=2.0, kind="physical",
         category="Membrane mechanics", citation="haselwandter2018",
         source_note="0.1 k_BT/nm^2 = 0.42 mN/m, a mildly tensioned membrane",
         description="Lateral tension entering the Helfrich functional."),
    dict(key="membrane.small_slope_limit", name="Small-slope validity limit",
         value=0.5, unit="", minimum=0.05, maximum=2.0, kind="method",
         category="Membrane mechanics", citation="method_choice",
         source_note="the Monge expansion drops terms of order |grad h|^2, so "
                     "at slope 0.5 the neglected terms are already ~25%",
         description="Slope above which the linearised solver is flagged."),
    dict(key="membrane.piezo1_contact_slope", name="PIEZO1 dome contact slope",
         value=1.992, unit="", minimum=0.1, maximum=4.0, kind="empirical",
         category="Membrane mechanics", citation="measured_here",
         source_note="a/sqrt(R^2-a^2) from the fitted 7WLT cap; 63 degrees",
         description="Slope where the dome meets the flat bilayer."),

    # -------------------------------------------------------------- dome ---
    dict(key="dome.delta_g0", name="Intrinsic gating free energy", value=9.7,
         unit="k_BT", minimum=0.0, maximum=40.0, kind="empirical",
         category="Dome energetics", citation="cox2016",
         source_note="fitted with dA = 8 nm^2 to give T50 = 5.1 mN/m",
         description="Closed-to-open free energy at zero tension."),
    dict(key="dome.delta_area", name="Gating area change", value=8.0,
         unit="nm^2", minimum=1.0, maximum=120.0, kind="empirical",
         category="Dome energetics", citation="cox2016",
         source_note="8 +- 1 nm^2 from Boltzmann fits to tension-response",
         description="In-plane area change on opening."),
    dict(key="dome.published_radius_closed", name="Published closed dome radius",
         value=10.2, unit="nm", minimum=5.0, maximum=20.0, kind="empirical",
         category="Dome energetics", citation="haselwandter2018",
         source_note="the regression target for our geometry pipeline",
         description="Radius of curvature of the closed-state dome."),

    # ---------------------------------------------------------- kinetics ---
    dict(key="kinetics.sigma_50", name="Tension midpoint", value=1.4,
         unit="mN/m", minimum=0.1, maximum=10.0, kind="empirical",
         category="Gating kinetics", citation="young2023",
         source_note="four-state TENSION model", description="Rate midpoint."),
    dict(key="kinetics.b", name="Tension sensitivity", value=0.8, unit="mN/m",
         minimum=0.05, maximum=5.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="exponential tension scale",
         description="e-fold tension for the rate constants."),
    dict(key="kinetics.k1_0", name="C->O rate prefactor", value=5.1, unit="1/s",
         minimum=0.0, maximum=1000.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="k1 = 5.1 exp(sigma/b)",
         description="Opening rate at zero tension."),
    dict(key="kinetics.km1_0", name="O->C rate prefactor", value=5.0, unit="1/s",
         minimum=0.0, maximum=1000.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="", description="Closing rate."),
    dict(key="kinetics.k2", name="O->I1 rate", value=8.0, unit="1/s",
         minimum=0.0, maximum=1000.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="", description="Inactivation entry."),
    dict(key="kinetics.km2", name="I1->O rate", value=0.4, unit="1/s",
         minimum=0.0, maximum=1000.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="", description="Inactivation exit."),
    dict(key="kinetics.k3_0", name="I1->I2 rate prefactor", value=34.6, unit="1/s",
         minimum=0.0, maximum=1000.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="k3 = 34.6 exp(-sigma/b)",
         description="Deep-inactivation entry."),
    dict(key="kinetics.k4", name="I2->C rate", value=4.0, unit="1/s",
         minimum=0.0, maximum=1000.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="", description="Recovery."),
    dict(key="kinetics.km4", name="C->I2 rate", value=0.6, unit="1/s",
         minimum=0.0, maximum=1000.0, kind="empirical", category="Gating kinetics",
         citation="young2023", source_note="k-3 follows by microscopic "
         "reversibility rather than being fitted", description="Reverse recovery."),
    dict(key="kinetics.conductance_pS", name="Unitary conductance", value=30.0,
         unit="pS", minimum=1.0, maximum=200.0, kind="empirical",
         category="Gating kinetics", citation="shi2020",
         source_note="25-30 pS; also Vaisey & MacKinnon 2026",
         description="Single-channel conductance."),
    dict(key="kinetics.wt_tau_ms", name="Wild-type inactivation tau", value=8.6,
         unit="ms", minimum=0.1, maximum=500.0, kind="empirical",
         category="Gating kinetics", citation="bae2013",
         source_note="whole-cell; mutants are calibrated by FOLD CHANGE against "
                     "this, never by absolute tau across preparations",
         description="Reference inactivation time constant."),
    dict(key="kinetics.t50_measured", name="Measured half-activation tension",
         value=2.7, unit="mN/m", minimum=0.1, maximum=20.0, kind="empirical",
         category="Gating kinetics", citation="lewis2015",
         source_note="2.7 +- 0.1 mN/m cell-attached",
         description="Validation target for the emergent T50."),

    # ------------------------------------------------------ elastic network --
    dict(key="anm.cutoff", name="Spring cutoff", value=15.0, unit="A",
         minimum=6.0, maximum=30.0, kind="method", category="Elastic network",
         citation="atilgan2001",
         source_note="the usual choice for large membrane proteins",
         description="C-alpha pairs closer than this are connected."),
    dict(key="anm.d0", name="Inverse-power spring scale", value=7.5, unit="A",
         minimum=1.0, maximum=20.0, kind="method", category="Elastic network",
         citation="method_choice",
         source_note="distance scale in gamma ~ (d0/d)^2; reproduces measured "
                     "fluctuations better than a uniform cutoff",
         description="Reference distance for distance-weighted springs."),
    dict(key="anm.gamma", name="Spring constant", value=1.0, unit="arb",
         minimum=0.01, maximum=100.0, kind="method", category="Elastic network",
         citation="method_choice",
         source_note="sets the energy scale only; overlaps and mode shapes are "
                     "invariant to it", description="Uniform spring stiffness."),
    dict(key="anm.n_modes", name="Modes to solve", value=20, unit="",
         minimum=6, maximum=200, kind="method", category="Elastic network",
         citation="method_choice",
         source_note="the first 6 per connected component are rigid-body and "
                     "are discarded", description="Number of low modes."),
    dict(key="anm.symmetry_tolerance", name="C3 character tolerance", value=0.25,
         unit="", minimum=0.01, maximum=1.0, kind="method",
         category="Elastic network", citation="method_choice",
         source_note="how close a rotation character must sit to +1 to call a "
                     "mode A-symmetric", description="Symmetry labelling tolerance."),
    dict(key="anm.n_protomers", name="Protomers", value=3, unit="",
         minimum=1, maximum=6, kind="physical", category="Elastic network",
         citation="ge2015", source_note="PIEZO1 is a C3 homotrimer",
         description="Subunits in the assembly."),

    # -------------------------------------------------------------- pore ---
    dict(key="pore.step", name="Slice spacing", value=1.0, unit="A",
         minimum=0.25, maximum=5.0, kind="method", category="Pore geometry",
         citation="method_choice", source_note="sampling at 1 A resolves the CTD constrictions, which are ~1 A wide",
         description="Spacing of slices along the conduction axis."),
    dict(key="pore.leash", name="Probe leash", value=8.0, unit="A",
         minimum=1.0, maximum=30.0, kind="method", category="Pore geometry",
         citation="method_choice",
         source_note="MANDATORY: without it the clearance function has no "
                     "interior maximum and the answer escapes to ~6000 A",
         description="Maximum probe offset from the axis."),
    dict(key="pore.search", name="Search radius", value=18.0, unit="A",
         minimum=5.0, maximum=50.0, kind="method", category="Pore geometry",
         citation="method_choice", source_note="wide enough to include every atom that could bound the probe at the leash limit",
         description="Neighbour radius when maximising clearance."),
    dict(key="pore.constriction_threshold", name="Constriction threshold",
         value=3.0, unit="A", minimum=0.5, maximum=10.0, kind="method",
         category="Pore geometry", citation="method_choice",
         source_note="a local minimum narrower than this is a candidate gate",
         description="Radius below which a minimum is reported."),
    dict(key="pore.ion_radius", name="Conducting ion radius", value=1.6, unit="A",
         minimum=0.5, maximum=5.0, kind="physical", category="Pore geometry",
         citation="smart1996hole",
         source_note="HOLE's convention for a permeant cation",
         description="Radius a pore must exceed to pass an ion sterically."),

    # ------------------------------------------------------- hydration ----
    dict(key="hydration.energy_threshold", name="Dewetting energy threshold",
         value=2.6, unit="kJ/mol", minimum=0.5, maximum=20.0, kind="empirical",
         category="Hydrophobic gating", citation="rao2019heuristic",
         source_note="1 RT; the contour CHAP draws",
         description="Water free energy above which a residue is flagged."),
    dict(key="hydration.closed_cutoff", name="Closed-gate score cutoff",
         value=0.55, unit="", minimum=0.05, maximum=3.0, kind="empirical",
         category="Hydrophobic gating", citation="rao2019heuristic",
         source_note="optimal cutoff from their ROC analysis, AUROC 0.91",
         description="Sum-of-distances above which a channel is called closed."),
    dict(key="hydration.kernel_bandwidth", name="Hydrophobicity bandwidth",
         value=0.35, unit="nm", minimum=0.05, maximum=2.0, kind="method",
         category="Hydrophobic gating", citation="klesse2019chap",
         source_note="CHAP's default; the landscape was built with it, so "
                     "changing it indexes the grid with a different coordinate",
         description="Gaussian kernel width along the pore axis."),
    dict(key="hydration.water_radius", name="Water molecule radius", value=0.15,
         unit="nm", minimum=0.05, maximum=0.5, kind="physical",
         category="Hydrophobic gating", citation="rao2019heuristic",
         source_note="below this a pore is shut sterically and the wetting "
                     "question does not arise", description="Water radius."),
    dict(key="hydration.max_radius", name="Scored pore radius limit", value=7.0,
         unit="A", minimum=1.0, maximum=20.0, kind="method",
         category="Hydrophobic gating", citation="klesse2019chap",
         source_note="CHAP restricts the score to residues lining the narrow "
                     "part (0.7 nm)", description="Widest slice that is scored."),

    # ---------------------------------------------------------- hybrid ----
    dict(key="hybrid.plddt_confident", name="AlphaFold confidence threshold",
         value=70.0, unit="pLDDT", minimum=0.0, maximum=100.0, kind="convention",
         category="Hybrid model", citation="convention",
         source_note="AlphaFold's own banding calls 70-90 confident and below "
                     "70 low; the project adopts that boundary rather than "
                     "choosing its own, so the reported confident fraction "
                     "means what a reader of AlphaFold DB expects",
         description="pLDDT at or above which a predicted residue is treated "
                     "as placed."),
    dict(key="hybrid.anchor_window", name="Hybrid graft anchor window",
         value=200, unit="residues", minimum=20, maximum=2000, kind="method",
         category="Hybrid model", citation="method_choice",
         source_note="the graft is anchored on residues NEAR the seam, not on "
                     "the whole overlap: fitting all 1279 shared residues gives "
                     "19.0 A RMSD because the AlphaFold and cryo-EM blades are "
                     "different conformations, and spreading that error into "
                     "the join misplaces the graft. A 200-residue window uses "
                     "110 residues and fits to 2.4 A",
         description="Residues either side of the seam used to place the "
                     "predicted distal blade."),

    # -------------------------------------------------------- geometry ----
    dict(key="geometry.sphere_trim", name="Sphere-fit trim fraction", value=0.15,
         unit="fraction", minimum=0.0, maximum=0.5, kind="method",
         category="Dome geometry", citation="method_choice",
         source_note="worst-fitting points ignored, so a few outlying helices "
                     "do not drag the fitted curvature",
         description="Fraction of points trimmed when fitting the sphere."),
    dict(key="geometry.radial_bins", name="Radial profile bins", value=40,
         unit="", minimum=5, maximum=200, kind="method",
         category="Dome geometry", citation="method_choice", source_note="40 bins over a ~10 nm dome gives ~2.5 A resolution, finer than the C-alpha spacing",
         description="Bins in the radial height profile."),
    dict(key="geometry.min_ca_per_protomer", name="Minimum C-alphas per protomer",
         value=300, unit="", minimum=20, maximum=2000, kind="method",
         category="Dome geometry", citation="method_choice",
         source_note="deposited entries often carry a short peptide that would "
                     "otherwise be mistaken for a fourth subunit",
         description="Chain length below which a chain is not a protomer."),

    # ---------------------------------------------------- interactions ----
    dict(key="sasa.probe_radius", name="SASA probe radius", value=1.4, unit="A",
         minimum=0.5, maximum=5.0, kind="physical", category="Interactions",
         citation="shrake1973", source_note="the standard water probe",
         description="Rolling-probe radius for solvent-accessible surface."),
    dict(key="sasa.n_points", name="SASA sphere points", value=256, unit="",
         minimum=32, maximum=4096, kind="method", category="Interactions",
         citation="shrake1973", source_note="numerical quadrature density",
         description="Points per atom in the Shrake-Rupley quadrature."),
    dict(key="interactions.min_sequence_separation",
         name="Minimum sequence separation", value=2, unit="residues",
         minimum=0, maximum=20, kind="method", category="Interactions",
         citation="method_choice",
         source_note="excludes trivial i,i+1 backbone contacts",
         description="Residues closer than this in sequence are ignored."),

    # -------------------------------------------------------- pockets ----
    dict(key="pockets.r_min", name="Alpha-sphere minimum radius", value=3.0,
         unit="A", minimum=1.0, maximum=10.0, kind="method", category="Pockets",
         citation="leguilloux2009", source_note="fpocket's default",
         description="Smallest alpha sphere kept."),
    dict(key="pockets.r_max", name="Alpha-sphere maximum radius", value=5.5,
         unit="A", minimum=2.0, maximum=20.0, kind="method", category="Pockets",
         citation="leguilloux2009", source_note="fpocket's default",
         description="Largest alpha sphere kept."),
    dict(key="pockets.min_neighbours", name="Burial neighbour count", value=30,
         unit="", minimum=0, maximum=200, kind="method", category="Pockets",
         citation="method_choice",
         source_note="without this filter the top pocket percolates over the "
                     "whole surface at 408,000 A^3 across 601 residues",
         description="Protein atoms required near a sphere for it to count."),
    dict(key="pockets.neighbour_radius", name="Burial neighbour radius",
         value=8.0, unit="A", minimum=2.0, maximum=20.0, kind="method",
         category="Pockets", citation="method_choice", source_note="one shell beyond the largest alpha sphere, so burial is judged on contacting protein",
         description="Radius within which burial neighbours are counted."),
    dict(key="pockets.ligand_cutoff", name="Ligand contact cutoff", value=4.5,
         unit="A", minimum=2.0, maximum=10.0, kind="method", category="Pockets",
         citation="method_choice", source_note="a heavy-atom contact distance; slightly beyond a hydrogen bond",
         description="Distance defining a residue-ligand contact."),

    # ------------------------------------------------------ allostery ----
    dict(key="allostery.contact_cutoff", name="Network contact cutoff",
         value=10.0, unit="A", minimum=4.0, maximum=25.0, kind="method",
         category="Allostery", citation="method_choice", source_note="the standard C-alpha contact distance for residue networks; shorter than the ANM cutoff because a network edge is a contact, not a spring",
         description="C-alpha separation defining a network edge."),
    dict(key="allostery.min_correlation", name="Minimum edge correlation",
         value=0.001, unit="", minimum=0.0, maximum=1.0, kind="method",
         category="Allostery", citation="method_choice",
         source_note="guards against log(0) in the edge weight",
         description="Correlation below which an edge is dropped."),

    # --------------------------------------------------- conservation ----
    dict(key="conservation.taxon", name="Ortholog clade", value=7742, unit="NCBI taxid",
         minimum=1, maximum=3000000, kind="method", category="Conservation",
         citation="method_choice", source_note="7742 = Vertebrata",
         description="Clade searched for orthologs."),
    dict(key="conservation.min_coverage", name="Minimum ortholog coverage",
         value=0.7, unit="fraction", minimum=0.0, maximum=1.0, kind="method",
         category="Conservation", citation="method_choice",
         source_note="below this a value measures the alignment rather than "
                     "selection pressure on the residue",
         description="Fraction of orthologs that must align at a position."),

    dict(key="conservation.constrained_threshold",
         name="Constrained-position threshold", value=0.9, unit="fraction",
         minimum=0.5, maximum=1.0, kind="method", category="Conservation",
         citation="method_choice",
         source_note="a position must be this conserved to be called "
                     "constrained when ranking candidates with no variant",
         description="Conservation above which a position counts as constrained."),
    dict(key="sasa.n_points_fast", name="SASA points, feature table", value=64,
         unit="", minimum=16, maximum=1024, kind="method",
         category="Interactions", citation="shrake1973",
         source_note="coarser quadrature for the 1279-residue feature table, "
                     "where relative SASA is used as a rank rather than a value",
         description="Quadrature density when building the feature table."),
    dict(key="measure.hydrophobicity_radius",
         name="Pore hydrophobicity contact radius", value=8.0, unit="A",
         minimum=2.0, maximum=20.0, kind="method", category="Interactions",
         citation="method_choice",
         source_note="Kyte-Doolittle averaging radius for the legacy pore "
                     "hydrophobicity profile; the CHAP path uses "
                     "hydration.kernel_bandwidth instead",
         description="Radius within which lining residues are averaged."),
    dict(key="pockets.buriedness_reach", name="Buriedness ray length", value=14.0,
         unit="A", minimum=3.0, maximum=40.0, kind="method", category="Pockets",
         citation="method_choice",
         source_note="how far a cast ray looks for protein before calling that "
                     "direction open", description="Ray length for buriedness."),
    dict(key="pockets.buriedness_clearance", name="Buriedness ray clearance",
         value=2.6, unit="A", minimum=0.5, maximum=10.0, kind="method",
         category="Pockets", citation="method_choice",
         source_note="a ray passing within this of an atom centre counts as "
                     "blocked; roughly a heavy-atom radius",
         description="Miss distance below which a ray is blocked."),
    dict(key="pockets.cluster_distance", name="Alpha-sphere cluster distance",
         value=2.0, unit="A", minimum=0.5, maximum=10.0, kind="method",
         category="Pockets", citation="leguilloux2009",
         source_note="single-linkage distance joining alpha spheres into one "
                     "pocket", description="Clustering distance for pockets."),
    dict(key="pockets.lining_cutoff", name="Pocket lining cutoff", value=5.0,
         unit="A", minimum=2.0, maximum=15.0, kind="method", category="Pockets",
         citation="method_choice",
         source_note="distance from an alpha sphere within which a residue is "
                     "called lining", description="Pocket lining distance."),

    # ------------------------------------------ substitution perturbation --
    dict(key="substitution.contact_length", name="Side-chain reach", value=6.0,
         unit="A", minimum=2.0, maximum=15.0, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="a C-alpha network places contacts at C-alpha separation, "
                     "but a side-chain change is felt by what the side chain "
                     "touches; roughly an extended side-chain length",
         description="Distance over which a contact feels a side-chain change."),
    dict(key="substitution.weight_volume", name="Packing weight", value=1.0,
         unit="", minimum=0.0, maximum=5.0, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="the original volume term, kept at unit weight so the new "
                     "terms are additions rather than a reweighting",
         description="How strongly a volume change scales a contact."),
    dict(key="substitution.weight_charge", name="Charge weight", value=0.6,
         unit="", minimum=0.0, maximum=5.0, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="an ENM spring is an effective stiffness standing in for "
                     "packing, hydrogen bonds and ion pairs together; losing a "
                     "salt bridge is a large fraction of a contact's stiffness",
         description="How strongly a charge change scales a contact to a "
                     "charged partner."),
    dict(key="substitution.weight_hbond", name="Hydrogen-bond weight",
         value=0.3, unit="", minimum=0.0, maximum=5.0, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="weaker than the charge term because a hydrogen bond "
                     "contributes less to an effective stiffness than an ion "
                     "pair at the same separation",
         description="How strongly a change in donor/acceptor capacity scales "
                     "a contact."),
    dict(key="substitution.weight_proline", name="Proline backbone weight",
         value=0.5, unit="", minimum=-2.0, maximum=5.0, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="proline restrains phi and removes a backbone donor, "
                     "stiffening sequence-local contacts specifically",
         description="Stiffening applied to sequence-local contacts when "
                     "proline is introduced."),
    dict(key="substitution.proline_span", name="Proline effect span", value=4,
         unit="residues", minimum=1, maximum=20, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="one turn of helix either side, the range over which a "
                     "backbone restraint is felt",
         description="Sequence separation within which proline stiffens."),
    dict(key="substitution.weight_glycine", name="Glycine weight", value=-0.4,
         unit="", minimum=-5.0, maximum=2.0, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="negative: removing the side chain leaves nothing to "
                     "mediate the contact, so it softens",
         description="Softening applied when glycine is introduced."),
    dict(key="substitution.min_scale", name="Minimum spring scale", value=0.05,
         unit="", minimum=0.001, maximum=1.0, kind="method",
         category="Substitution model", citation="method_choice",
         source_note="a spring may weaken but never invert, or the quadratic "
                     "form stops being an energy and the Hessian stops being "
                     "positive semi-definite",
         description="Floor on the per-contact spring scale."),

    dict(key="crosscheck.min_pair_separation",
         name="Minimum pair separation, distance-space overlap", value=8.0,
         unit="A", minimum=2.0, maximum=40.0, kind="method",
         category="Cross-checks", citation="method_choice",
         source_note="pairs closer than this barely change length during the "
                     "transition, so they add noise rather than signal to the "
                     "distance-space correlation",
         description="Shortest C-alpha pair used in the superposition-free "
                     "overlap."),

    # -------------------------------------------------------- statistics --
    dict(key="stats.alpha", name="Significance level", value=0.05, unit="",
         minimum=0.001, maximum=0.2, kind="convention", category="Statistics",
         citation="convention", source_note="pre-registered for both variant tests",
         description="Type-I error rate."),
    dict(key="stats.target_power", name="Target power", value=0.8, unit="",
         minimum=0.5, maximum=0.99, kind="convention", category="Statistics",
         citation="convention", source_note="the conventional 80%; both power statements in the project are quoted at it",
         description="Power at which the minimum detectable effect is quoted."),
    dict(key="stats.n_permutations", name="Permutation shuffles", value=10000,
         unit="", minimum=99, maximum=1000000, kind="method",
         category="Statistics", citation="method_choice",
         source_note="p-values use the (r+1)/(n+1) convention, so they can "
                     "never be exactly zero", description="Label shuffles."),
    dict(key="stats.n_bootstrap", name="Bootstrap resamples", value=10000,
         unit="", minimum=99, maximum=1000000, kind="method",
         category="Statistics", citation="method_choice", source_note="10k resamples give a percentile interval stable to about the third decimal",
         description="Resamples for the effect-size confidence interval."),

]

P += TAG_PARAMETERS
P += FIGURE_PARAMETERS
