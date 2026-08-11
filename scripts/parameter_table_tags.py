"""Parameters for the HaloTag work — the fusion geometry and the labelling.

Split out of ``parameter_table.py`` to keep both under the project's 500-line
limit, and along a real seam: everything here describes the *tag* rather than
the channel, and was added in Rounds 31 and 32.
"""

from __future__ import annotations

__all__ = ["TAG_PARAMETERS"]

TAG_PARAMETERS = [
    # ---------------------------------------------------- HaloTag fusion ----
    dict(key="fusion.linker_residues", name="HaloTag linker length", value=10,
         unit="residues", minimum=0, maximum=40, kind="convention",
         category="HaloTag fusion", citation="unverified",
         source_note="The construct this project is meant to describe does not "
                     "state its linker, and neither does the halotag_binding_sim "
                     "project it comes from. 10 residues is a mid-range flexible "
                     "linker for a C-terminal fusion. It sets how far the tag can "
                     "reach, so it is the single assumption the fusion geometry "
                     "rests on and must be varied, not trusted.",
         description="Residues between PIEZO1's C-terminus and the HaloTag N-terminus."),
    dict(key="fusion.residue_extension", name="Extension per linker residue",
         value=3.5, unit="A", minimum=1.0, maximum=4.0, kind="convention",
         category="HaloTag fusion", citation="convention",
         source_note="Ca-Ca rise of a fully extended polypeptide, the standard "
                     "value. It is an upper bound: a real flexible linker is "
                     "coiled and reaches less, so the envelope errs generous.",
         description="Maximum reach contributed by each linker residue."),
    dict(key="fusion.grid_spacing", name="Accessible-volume grid spacing",
         value=2.0, unit="A", minimum=0.5, maximum=10.0, kind="method",
         category="HaloTag fusion", citation="method_choice",
         source_note="resolution of the envelope; finer costs cubically and "
                     "does not change where the envelope sits",
         description="Grid step for sampling accessible tag positions."),
    dict(key="fusion.pore_mouth_radius", name="Pore-mouth radius",
         value=15.0, unit="A", minimum=5.0, maximum=40.0, kind="method",
         category="HaloTag fusion", citation="method_choice",
         source_note="how close to the conduction axis an atom must be to count "
                     "as lining the cytosolic mouth rather than sitting on a "
                     "blade. The CTD bundle is well inside it and the blades "
                     "far outside, so the answer is flat in between; taking the "
                     "lowest atom anywhere instead finds a blade tip and moves "
                     "the pore exit by over 10 nm",
         description="Radius about the axis defining the cytosolic pore mouth."),
    dict(key="fusion.clash_clearance", name="Tag-channel clearance", value=2.0,
         unit="A", minimum=0.0, maximum=10.0, kind="method",
         category="HaloTag fusion", citation="method_choice",
         source_note="added to the tag's own radius when rejecting positions, "
                     "so contact rather than overlap is what is excluded",
         description="Extra clearance required between tag body and channel."),
    dict(key="fusion.pose_contact_distance", name="Tag-atom contact distance",
         value=3.4, unit="A", minimum=2.0, maximum=6.0, kind="method",
         category="HaloTag fusion", citation="method_choice",
         source_note="twice Bondi's carbon van der Waals radius, so a pair at "
                     "this separation is touching. It applies only to the "
                     "drawn fold: it sets how many atoms are reported in "
                     "contact with the channel and which of the undetermined "
                     "spins is drawn, and changes no modelled position. The "
                     "finding is insensitive to it — no spin is contact-free "
                     "on 8YEZ, 7WLT or 11ZC anywhere from 2.0 to 4.0 A, "
                     "except 7WLT below 3.4 A.",
         description="Separation below which a drawn tag atom counts as "
                     "touching the channel."),

    # ------------------------------------------------- HaloTag labelling ----
    # Imported from the halotag_binding_sim project in Round 32. The model is
    # three equations: exposure E(t), per-site p(t) = a*(1-exp(-k_on*E)), and
    # Binomial(3, p) over the trimer. Everything below feeds one of them.
    dict(key="labelling.k_on", name="HaloTag covalent on-rate", value=2.7e6,
         unit="1/(M s)", minimum=1e4, maximum=1e8, kind="physical",
         category="HaloTag labelling", citation="los2008halotag",
         source_note="apparent second-order rate for the HaloTag-chloroalkane "
                     "reaction, which is covalent and does not reverse",
         description="Second-order rate constant for ligand capture by one tag."),
    dict(key="labelling.n_sites", name="HaloTags per channel", value=3,
         unit="count", minimum=1, maximum=6, kind="physical",
         category="HaloTag labelling", citation="bertaccini2025piezo1",
         source_note="PIEZO1 is a homotrimer, so a C-terminal fusion puts one "
                     "tag on each protomer; this is the exponent in p^3",
         description="Independent labelling sites on one channel."),
    dict(key="labelling.k_perm_live", name="Ligand access rate, living cell",
         value=0.008333333333333333, unit="1/s", minimum=1e-5, maximum=1e7,
         kind="method", category="HaloTag labelling", citation="unverified",
         source_note="tau ~ 2 min for a cell-permeable JF dye to equilibrate "
                     "across the plasma membrane. A model estimate in the "
                     "source project too, not a measurement; it sets how much "
                     "of the labelling lag is transport rather than chemistry.",
         description="Membrane permeation rate for the ligand in a live cell."),
    dict(key="labelling.partition_live", name="Ligand partition, living cell",
         value=1.0, unit="ratio", minimum=0.0, maximum=10.0, kind="empirical",
         category="HaloTag labelling", citation="grimm2015jf",
         source_note="JF dyes are cell-permeable, so intracellular concentration "
                     "reaches the bath concentration at steady state",
         description="Steady-state intracellular / bath ligand concentration."),
    dict(key="labelling.active_fraction", name="Reactive fraction of tags",
         value=1.0, unit="fraction", minimum=0.0, maximum=1.0,
         kind="method", category="HaloTag labelling", citation="unverified",
         source_note="assumed fully reactive in a living cell. It caps p, so "
                     "the fully-labelled asymptote is active_fraction^3 and no "
                     "incubation time can beat it; fixation is believed to "
                     "lower it but this project has not measured that.",
         description="Fraction of tags chemically able to react."),
    dict(key="labelling.concentration", name="Bath ligand concentration",
         value=2e-7, unit="M", minimum=1e-12, maximum=1e-3, kind="convention",
         category="HaloTag labelling", citation="bertaccini2025piezo1",
         source_note="200 nM, the standard live-cell JF646 labelling condition",
         description="Applied ligand concentration in the bath."),
    dict(key="labelling.incubation_time", name="Incubation time", value=1800.0,
         unit="s", minimum=0.0, maximum=86400.0, kind="convention",
         category="HaloTag labelling", citation="bertaccini2025piezo1",
         source_note="30 min, the standard live-cell labelling protocol",
         description="Ligand incubation time before wash."),
    dict(key="labelling.brightness_noise_cv", name="Per-dye brightness CV",
         value=0.15, unit="fraction", minimum=0.0, maximum=2.0, kind="method",
         category="HaloTag labelling", citation="method_choice",
         source_note="spread applied to each dye's contribution when turning "
                     "integer dye counts into a predicted amplitude histogram; "
                     "it decides whether the 1/2/3-dye levels stay resolvable",
         description="Coefficient of variation of one dye's brightness."),

    # ------------------------------------------------------- permeation ----
    # A 1-D PNP model over the measured pore profile. The transport constants
    # are bulk-water measurements; the two confinement corrections beneath them
    # are not measured, and they are what the answer is most sensitive to.
    dict(key="permeation.temperature", name="Temperature", value=310.0,
         unit="K", minimum=270.0, maximum=330.0, kind="convention",
         category="Permeation", citation="convention",
         source_note="37 C, the condition PIEZO1 conductance is usually recorded at",
         description="Absolute temperature entering the thermal voltage."),
    dict(key="permeation.bath_concentration", name="Bath salt concentration",
         value=0.15, unit="M", minimum=0.001, maximum=3.0, kind="convention",
         category="Permeation", citation="gnanasambandam2015",
         source_note="symmetric ~150 mM monovalent, the standard recording "
                     "condition for the 25-30 pS figure",
         description="Symmetric bath concentration of the permeant salt."),
    dict(key="permeation.test_voltage", name="Test voltage", value=0.06,
         unit="V", minimum=0.001, maximum=0.2, kind="method",
         category="Permeation", citation="method_choice",
         source_note="60 mV; small enough that the I-V is still linear, so the "
                     "chord conductance equals the slope conductance",
         description="Voltage at which the conductance is evaluated."),
    dict(key="permeation.diffusion_cation", name="Cation bulk diffusivity",
         value=1.96e-9, unit="m^2/s", minimum=1e-11, maximum=1e-8,
         kind="convention", category="Permeation", citation="convention",
         source_note="K+ in water at 25 C; the standard tabulated value",
         description="Bulk diffusion coefficient of the permeant cation."),
    dict(key="permeation.diffusion_anion", name="Anion bulk diffusivity",
         value=2.03e-9, unit="m^2/s", minimum=1e-11, maximum=1e-8,
         kind="convention", category="Permeation", citation="convention",
         source_note="Cl- in water at 25 C; the standard tabulated value",
         description="Bulk diffusion coefficient of the counter-anion."),
    dict(key="permeation.diffusion_calcium", name="Calcium bulk diffusivity",
         value=0.79e-9, unit="m^2/s", minimum=1e-11, maximum=1e-8,
         kind="convention", category="Permeation", citation="convention",
         source_note="Ca2+ in water at 25 C; slower than K+ because it is more "
                     "strongly hydrated",
         description="Bulk diffusion coefficient of calcium."),
    dict(key="permeation.radius_cation", name="Cation crystal radius",
         value=1.38, unit="A", minimum=0.3, maximum=4.0, kind="convention",
         category="Permeation", citation="convention",
         source_note="K+ Shannon radius. The crystal rather than hydrated "
                     "radius, because an ion can shed water to cross a "
                     "constriction - which is itself a modelling choice",
         description="Radius subtracted from the pore when computing access."),
    dict(key="permeation.radius_anion", name="Anion crystal radius", value=1.81,
         unit="A", minimum=0.3, maximum=4.0, kind="convention",
         category="Permeation", citation="convention",
         source_note="Cl- Shannon radius", description="Anion crystal radius."),
    dict(key="permeation.radius_calcium", name="Calcium crystal radius",
         value=1.00, unit="A", minimum=0.3, maximum=4.0, kind="convention",
         category="Permeation", citation="convention",
         source_note="Ca2+ Shannon radius", description="Calcium crystal radius."),
    dict(key="permeation.diffusion_scale", name="In-pore diffusivity scaling",
         value=0.5, unit="fraction", minimum=0.01, maximum=1.0, kind="method",
         category="Permeation", citation="unverified",
         source_note="Ions move more slowly in a pore than in bulk water, and "
                     "by how much is not measured for PIEZO1. Together with "
                     "the ion radius this is what the computed conductance is "
                     "most sensitive to, so it must be swept, not trusted.",
         description="In-pore diffusivity as a fraction of the bulk value."),
    dict(key="permeation.permittivity_pore", name="In-pore relative permittivity",
         value=40.0, unit="", minimum=2.0, maximum=80.0, kind="method",
         category="Permeation", citation="unverified",
         source_note="water in a nanopore is less polarisable than in bulk "
                     "(80); 40 is a common compromise. It only matters where "
                     "there is fixed charge to screen.",
         description="Relative permittivity used in the Poisson equation."),
    # The sodium set exists because the published permeability ratio was
    # measured in NaCl. The default "cation" above is potassium — smaller
    # diffusivity and larger radius are both wrong for Na+, and comparing a
    # potassium model against a sodium measurement would let it disagree for a
    # reason that has nothing to do with the pore.
    dict(key="permeation.diffusion_sodium", name="Sodium bulk diffusivity",
         value=1.33e-9, unit="m^2/s", minimum=1e-11, maximum=1e-8,
         kind="convention", category="Permeation", citation="convention",
         source_note="Na+ in water at 25 C; the standard tabulated value. "
                     "Slower than K+ because it holds its hydration shell "
                     "more tightly.",
         description="Bulk diffusion coefficient of sodium."),
    dict(key="permeation.radius_sodium", name="Sodium crystal radius",
         value=1.02, unit="A", minimum=0.3, maximum=4.0, kind="convention",
         category="Permeation", citation="convention",
         source_note="Na+ Shannon radius, six-coordinate",
         description="Sodium crystal radius."),
    dict(key="permeation.dilution_high", name="Dilution potential, high side",
         value=0.15, unit="M", minimum=0.001, maximum=3.0, kind="convention",
         category="Permeation", citation="coste2015",
         source_note="150 mM NaCl intracellular, the concentrated side of the "
                     "protocol P_Cl/P_Na was measured with",
         description="Bath concentration on the concentrated side."),
    dict(key="permeation.dilution_low", name="Dilution potential, dilute side",
         value=0.03, unit="M", minimum=0.001, maximum=3.0, kind="convention",
         category="Permeation", citation="coste2015",
         source_note="30 mM NaCl extracellular (with sucrose for osmolarity), "
                     "the dilute side of the same protocol",
         description="Bath concentration on the dilute side."),
    dict(key="permeation.published_pcl_pna", name="Published P_Cl/P_Na",
         value=0.14, unit="", minimum=0.0, maximum=10.0, kind="empirical",
         category="Permeation", citation="coste2015",
         source_note="mPiezo1 permeates chloride with P_Cl/P_Na = 0.14, from "
                     "the dilution potential inverted through GHK. Below one, "
                     "so the channel prefers cations — but far from zero, so "
                     "it is a weakly selective cation channel and not an "
                     "ion-exchange membrane.",
         description="Measured chloride-over-sodium permeability ratio."),
    dict(key="permeation.published_conductance", name="Published unitary conductance",
         value=27.5, unit="pS", minimum=1.0, maximum=200.0, kind="empirical",
         category="Permeation", citation="coste2010piezo",
         source_note="25-30 pS for PIEZO1 with monovalent cations; the "
                     "midpoint is the comparison target for the PNP model",
         description="Measured single-channel conductance."),

    # --------------------------------------------------- pore fixed charge ----
    # Where a charged side chain can put its charge, measured from its own
    # C-alpha with the chain fully extended. These decide which residues count
    # as pore-lining when only C-alpha coordinates exist, which is the case for
    # the one open structure. Extended values on purpose: the test is meant to
    # err towards inclusion, so a residue it rejects cannot line the pore in
    # any rotamer.
    dict(key="pore_charge.reach_asp", name="Aspartate side-chain reach",
         value=3.9, unit="A", minimum=1.0, maximum=12.0, kind="convention",
         category="Pore fixed charge", citation="convention",
         source_note="Ca to carboxylate-carbon distance in an extended Asp, "
                     "from standard side-chain geometry",
         description="How far Asp can put its charge from its own C-alpha."),
    dict(key="pore_charge.reach_glu", name="Glutamate side-chain reach",
         value=5.1, unit="A", minimum=1.0, maximum=12.0, kind="convention",
         category="Pore fixed charge", citation="convention",
         source_note="one methylene longer than Asp",
         description="How far Glu can put its charge from its own C-alpha."),
    dict(key="pore_charge.reach_lys", name="Lysine side-chain reach",
         value=6.4, unit="A", minimum=1.0, maximum=12.0, kind="convention",
         category="Pore fixed charge", citation="convention",
         source_note="Ca to NZ in an all-trans lysine",
         description="How far Lys can put its charge from its own C-alpha."),
    dict(key="pore_charge.reach_arg", name="Arginine side-chain reach",
         value=7.3, unit="A", minimum=1.0, maximum=12.0, kind="convention",
         category="Pore fixed charge", citation="convention",
         source_note="Ca to the guanidinium centre in an extended arginine; "
                     "the longest ionisable side chain there is",
         description="How far Arg can put its charge from its own C-alpha."),
    dict(key="pore_charge.smoothing", name="Axial charge smoothing",
         value=3.0, unit="A", minimum=0.5, maximum=15.0, kind="method",
         category="Pore fixed charge", citation="method_choice",
         source_note="A charge's height along the axis is known to about a "
                     "side-chain length, and on the open structure only from "
                     "C-alpha. Placing it on a single slice would claim a "
                     "precision the coordinates do not have and would make the "
                     "density depend on the slice spacing. Total charge is "
                     "conserved whatever this is set to.",
         description="Gaussian width each charge is spread over along the axis."),
    dict(key="pore_charge.max_concentration",
         name="Counterion packing ceiling", value=10.0, unit="M",
         minimum=0.1, maximum=100.0, kind="method",
         category="Pore fixed charge", citation="method_choice",
         source_note="Close-packed hydrated K+ at a 3.3 A radius is about "
                     "8 M, so 10 M is where a continuum treatment of point "
                     "ions has certainly stopped describing a solution. It "
                     "only ever flags a result — no calculation reads it as "
                     "an input, and nothing is clipped to it.",
         description="In-pore concentration above which the model is flagged."),

    # ------------------------------------------------------- full-length model ----
    dict(key="full_length.gap_anchor_window", name="Gap anchor window",
         value=25.0, unit="residues", minimum=3.0, maximum=200.0,
         kind="method", category="Full-length model", citation="method_choice",
         source_note="How many resolved residues either side of an "
                     "unresolved stretch are used to place the prediction "
                     "into it. Both flanks are required, so the fill is "
                     "interpolated rather than cantilevered like the blade; "
                     "the flank RMSD is reported so a bad fit is visible "
                     "rather than absorbed.",
         description="Resolved residues each side of a gap used to anchor its fill."),
    dict(key="full_length.max_gap", name="Longest gap filled",
         value=120.0, unit="residues", minimum=1.0, maximum=2000.0,
         kind="method", category="Full-length model", citation="method_choice",
         source_note="A local fit at two flanks says nothing about what a very "
                     "long insert does in between, and a 400-residue loop "
                     "placed that way would look like structure. Gaps longer "
                     "than this are left empty and counted, so the model "
                     "reports what it did not fill.",
         description="Gaps longer than this are reported rather than filled."),

    # -------------------------------------------------- paralogue control ----
    dict(key="paralogue.min_sequence_identity",
         name="Numbering identification floor", value=0.9, unit="fraction",
         minimum=0.5, maximum=1.0, kind="method",
         category="Paralogue control", citation="method_choice",
         source_note="A deposited file's residue names must agree with the "
                     "sequence its numbering belongs to at essentially every "
                     "position. Measured across all 21 entries the winner "
                     "scores 1.000 and the runner-up below 0.25, so anything "
                     "under 0.9 means the file is not in the numbering it was "
                     "matched to and the comparison must not proceed.",
         description="Identity below which a numbering identification is refused."),
    dict(key="paralogue.min_identity_margin",
         name="Numbering identification margin", value=0.5, unit="fraction",
         minimum=0.0, maximum=1.0, kind="method",
         category="Paralogue control", citation="method_choice",
         source_note="A confident answer and a coin flip look identical once "
                     "only the winner is reported. The measured margin is at "
                     "least 0.75 on every entry, so half is a wide floor that "
                     "would still catch a genuinely ambiguous file.",
         description="How far the best reference must beat the second."),
    dict(key="paralogue.tm_boundary_slack",
         name="Transmembrane boundary slack", value=5.0, unit="residues",
         minimum=0.0, maximum=30.0, kind="method",
         category="Paralogue control", citation="method_choice",
         source_note="How far outside a UniProt helix an aligned midpoint may "
                     "land and still count as that helix. UniProt boundaries "
                     "are predictions with a residue or two of play and the "
                     "alignment adds its own. Measured: 37 of 38 helices "
                     "agree, unchanged from 0 to 8 residues of slack; above "
                     "that the count falls, and it falls because adjacent "
                     "PIEZO2 helices sit as little as 2 residues apart, so "
                     "wide windows overlap and the match stops being unique. "
                     "That is a limit of the criterion, not a disagreement, "
                     "and 5 sits in the middle of the flat range.",
         description="Tolerance when asking which helix an aligned residue is in."),

    dict(key="paralogue.max_splice_shift", name="Splice shift search range",
         value=200.0, unit="residues", minimum=1.0, maximum=2000.0,
         kind="method", category="Paralogue control", citation="method_choice",
         source_note="How far a splice-isoform renumbering is searched for "
                     "when the straight identification fails. The one case "
                     "here is 6LQI at 24 residues; the range is wide because "
                     "the search is only ever entered after a failure and a "
                     "wrong answer is refused by requiring both halves to "
                     "agree essentially perfectly.",
         description="Largest constant renumbering searched for after a splice site."),

    dict(key="paralogue.min_mismatch_block", name="Mismatch block floor",
         value=10.0, unit="residues", minimum=2.0, maximum=200.0,
         kind="method", category="Paralogue control", citation="method_choice",
         source_note="How long a run of disagreeing residues must be before "
                     "it is reported as a numbering fault rather than as real "
                     "sequence difference. An engineered variant changes one "
                     "residue and must not be flagged; the measured case is a "
                     "block of 91 in four human entries, so ten is far from "
                     "either.",
         description="Shortest run of disagreement reported as a block."),

    # ---------------------------------------------- fluctuation comparison ----
    dict(key="fluctuation.min_distinct_fraction",
         name="B-factor resolution floor", value=0.10, unit="fraction",
         minimum=0.0, maximum=1.0, kind="method",
         category="Fluctuation comparison", citation="method_choice",
         source_note="Below one distinct B-factor value per ten residues the "
                     "column was refined in groups and cannot resolve "
                     "per-residue mobility, so a correlation against it would "
                     "measure the grouping rather than the network. Measured "
                     "here: 3JAC and 6BPZ carry ~212 distinct values over "
                     "~2700 C-alphas and are the only entries this excludes.",
         description="Distinct B-factor values per residue below which an "
                     "entry is refused rather than compared."),

    # ------------------------------------------------- calcium nanodomain ----
    dict(key="nanodomain.d_calcium", name="Cytosolic calcium diffusivity",
         value=2.2e-10, unit="m^2/s", minimum=1e-12, maximum=1e-8,
         kind="empirical", category="Calcium nanodomain",
         citation="allbritton1992",
         source_note="220 um^2/s in cytosol - roughly a third of the value in "
                     "water, because cytosolic buffers and obstacles retard it",
         description="Effective diffusion coefficient of free Ca2+ in cytosol."),
    dict(key="nanodomain.buffer_kon", name="Buffer association rate", value=1e8,
         unit="1/(M s)", minimum=1e5, maximum=1e10, kind="empirical",
         category="Calcium nanodomain", citation="naraghi1997",
         source_note="a fast mobile endogenous buffer; with the concentration "
                     "below this sets the screening length",
         description="On-rate of the dominant cytosolic calcium buffer."),
    dict(key="nanodomain.buffer_concentration", name="Free buffer concentration",
         value=1e-4, unit="M", minimum=1e-7, maximum=1e-1, kind="empirical",
         category="Calcium nanodomain", citation="naraghi1997",
         source_note="~100 uM free mobile buffer, a typical cytosolic value. "
                     "Only the product k_on*[B] matters, through lambda.",
         description="Concentration of free calcium buffer in the cytosol."),
    dict(key="nanodomain.sensor_kd", name="Calcium sensor dissociation constant",
         value=2e-7, unit="M", minimum=1e-9, maximum=1e-3, kind="empirical",
         category="Calcium nanodomain", citation="tsien1980bapta",
         source_note="~0.2 uM for a BAPTA-based indicator; the threshold the "
                     "nanodomain is compared against",
         description="Kd of the BAPTA-based calcium sensor on the tag."),
    dict(key="nanodomain.calcium_current_fraction",
         name="Calcium share of the unitary current", value=0.05, unit="fraction",
         minimum=0.0, maximum=1.0, kind="method",
         category="Calcium nanodomain", citation="unverified",
         source_note="PIEZO1 is weakly cation-selective; the PNP model gives "
                     "under 5% at 2 mM Ca2+. Not measured for this channel, and "
                     "it scales the nanodomain linearly, so it is swept.",
         description="Fraction of the single-channel current carried by Ca2+."),
    dict(key="nanodomain.resting_calcium", name="Resting cytosolic calcium",
         value=1e-7, unit="M", minimum=1e-9, maximum=1e-5, kind="convention",
         category="Calcium nanodomain", citation="convention",
         source_note="~100 nM, the standard resting cytosolic level; the "
                     "baseline the nanodomain sits on top of",
         description="Bulk resting free calcium concentration."),
]
