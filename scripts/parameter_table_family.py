"""Parameters for the imported PIEZO-family census, split off at a real seam.

Every parameter in ``parameter_table.py`` describes this project's own model —
a bending modulus, a probe radius, a rate constant. Every parameter here
describes **how an external result is joined to it**: how many draws the null
uses, how far a correlation must sit above that null before it is reported as
anything, what counts as surviving the burial confound.

That is a different kind of number, and it needs its own justification. None of
these is a property of PIEZO1. Each is a decision about how much evidence is
enough, and putting them in the registry is what stops those decisions being
made silently inside a function and then quoted as a finding.
"""

from __future__ import annotations

__all__ = ["FAMILY_PARAMETERS"]

FAMILY_PARAMETERS = [
    dict(key="family.null_draws", name="Circular-shift null draws", value=200,
         unit="", minimum=20, maximum=10000, kind="method",
         category="Family census", citation="method_choice",
         source_note="each draw is one Spearman over ~1,300 residues, about a "
                     "millisecond; 200 gives a null standard deviation stable "
                     "to the third decimal, which is finer than the z "
                     "thresholds it is compared against",
         description="Random circular shifts of the constraint track used as "
                     "the null for a constraint-versus-mechanics correlation. "
                     "A shift is used rather than a permutation because both "
                     "series are strongly autocorrelated along the chain and a "
                     "permutation null is far too easy to beat."),
    dict(key="family.min_null_z", name="Minimum z over the shift null",
         value=2.0, unit="sigma", minimum=1.0, maximum=6.0, kind="method",
         category="Family census", citation="method_choice",
         source_note="two sigma, stated in advance and applied to every "
                     "feature alike; the comparison tests eight features at "
                     "once, so a single one clearing this is not a result on "
                     "its own and the report says how many cleared",
         description="How far a correlation must sit from the circular-shift "
                     "null before it is reported as surviving it."),
    dict(key="family.burial_retention", name="Correlation kept past burial",
         value=0.5, unit="fraction", minimum=0.1, maximum=1.0, kind="method",
         category="Family census", citation="method_choice",
         source_note="a mechanical feature must keep half its rank "
                     "correlation once burial is partialled out. Chosen before "
                     "the numbers were looked at, because burial explains most "
                     "conservation in most proteins and a threshold picked "
                     "afterwards would be a threshold picked to pass",
         description="Fraction of its raw Spearman a feature must retain with "
                     "relative SASA and contact number held fixed."),
    dict(key="family.population_af_floor",
         name="Population-variant allele-frequency floor", value=1e-4,
         unit="fraction", minimum=0.0, maximum=0.1, kind="method",
         category="Family census", citation="method_choice",
         source_note="the comparator has to be variation that exists in "
                     "people rather than variation somebody classified as "
                     "harmless. A floor of 1e-4 keeps roughly the positions "
                     "seen more than a handful of times across gnomAD's "
                     "cohort and drops the singletons, which are dominated by "
                     "sequencing depth rather than by tolerance",
         description="Minimum gnomAD allele frequency for a missense position "
                     "to count as population variation in the disease-geography "
                     "test."),
    dict(key="family.distal_last_thu", name="Last THU counted as distal blade",
         value=6, unit="", minimum=1, maximum=8, kind="convention",
         category="Family census", citation="convention",
         source_note="PIEZO1's blade is nine four-TM units; the census split "
                     "distal from proximal with a single chain cut, and this "
                     "project splits it by unit. Six is where the literature "
                     "puts the boundary, and it is registered rather than "
                     "written into the function because the census's own "
                     "distal-versus-proximal finding depends on where the "
                     "split falls and a reader should be able to move it",
         description="THU1 up to this index are the distal blade; the rest are "
                     "proximal."),
    dict(key="family.motif_window", name="Deep-conservation window width",
         value=9, unit="residues", minimum=3, maximum=51, kind="method",
         category="Family census", citation="method_choice",
         source_note="the census reports 'three short windows around the "
                     "pore' as what is conserved to whole-family depth. Nine "
                     "residues is between two and three turns of an alpha "
                     "helix, wide enough that a single well-conserved position "
                     "cannot carry a window and narrow enough to resolve one "
                     "helix face from another",
         description="Smoothing width for finding the stretches of human "
                     "PIEZO1 most conserved across the whole family."),
    dict(key="family.core_rmsd_ceiling", name="Core superposition ceiling",
         value=6.0, unit="angstrom", minimum=1.0, maximum=20.0, kind="method",
         category="Family census", citation="method_choice",
         source_note="the census measured 3.86 A over the pore module between "
                     "a piezo3 prediction and mouse Piezo1; this is the value "
                     "above which a core superposition is reported as not "
                     "having converged on a common core rather than as a "
                     "measurement of one",
         description="C-alpha RMSD above which a pore-module superposition is "
                     "not treated as evidence of a shared core."),
    dict(key="family.splay_ratio", name="Periphery-to-core splay ratio",
         value=2.0, unit="", minimum=1.0, maximum=20.0, kind="method",
         category="Family census", citation="method_choice",
         source_note="the finding being tested is that the cores agree while "
                     "the blades splay; the ratio has to exceed one for the "
                     "statement to mean anything and two is the margin at "
                     "which it is reported rather than noted",
         description="Blade RMSD divided by core RMSD, above which a pair is "
                     "reported as core-conserved and periphery-free."),
    dict(key="family.equivalent_ca_cutoff",
         name="Equivalent-position agreement cutoff", value=5.0,
         unit="angstrom", minimum=1.0, maximum=20.0, kind="method",
         category="Family census", citation="method_choice",
         source_note="two residues an alignment calls equivalent are reported "
                     "as occupying the same place only if their C-alpha land "
                     "within this distance after a core superposition. Five "
                     "angstrom is roughly one residue's step along a helix, so "
                     "a pair inside it cannot be a register error",
         description="C-alpha separation below which two aligned positions in "
                     "different paralogues are reported as the same place."),
]
