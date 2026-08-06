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
]
