"""The family-comparison parameters, split off at a real seam.

Every other parameter in this project describes a *channel* — a bending
modulus, a rate constant, a probe radius. These describe an **alignment**, and
the difference matters: an alignment score is a property of a scoring scheme
and a null model, not of a protein, so each of these has to name the convention
it inherits or the comparison stops being reproducible.

The three that carry weight are ``twilight_identity``, ``null_replicates`` and
``min_z``. Together they are the whole reason this project does not ship a
homology search: measured against a composition-matched shuffle, the percent
identity between human PIEZO1 and the plant or amoebal PIEZOs sits **three and a
half to six standard deviations** above its own null, while the local alignment
score on the same pair sits **eighty-five to a hundred and two**. The homology is
overwhelming and the percentage is nearly all noise. Anything reporting one
without the other is inviting a confident wrong reading.
"""

from __future__ import annotations

__all__ = ["HOMOLOGY_PARAMETERS"]

HOMOLOGY_PARAMETERS = [
    dict(key="homology.gap_open", name="Alignment gap-open score", value=-11.0,
         unit="", minimum=-30.0, maximum=-1.0, kind="convention",
         category="Family homology", citation="henikoff1992",
         source_note="the BLOSUM62 default pairing; changing it changes every "
                     "identity in the family matrix, so it is stated rather "
                     "than left to an aligner default that may move",
         description="Gap-opening penalty for every alignment in the family "
                     "comparison."),
    dict(key="homology.gap_extend", name="Alignment gap-extend score",
         value=-1.0, unit="", minimum=-10.0, maximum=-0.1, kind="convention",
         category="Family homology", citation="henikoff1992",
         source_note="the BLOSUM62 default pairing, as above",
         description="Gap-extension penalty."),
    dict(key="homology.twilight_identity", name="Twilight-zone identity",
         value=0.30, unit="fraction", minimum=0.1, maximum=0.6,
         kind="empirical", category="Family homology", citation="rost1999",
         source_note="Rost's boundary below which pairwise identity stops "
                     "separating homologues from unrelated pairs. Six of the "
                     "nine PIEZOs fall below it against human PIEZO1, which "
                     "is why this number is in the registry and on screen "
                     "rather than in a footnote",
         description="Identity below which a percentage is not evidence of "
                     "homology on its own."),
    dict(key="homology.null_replicates", name="Shuffled-null replicates",
         value=20, unit="", minimum=5, maximum=1000, kind="method",
         category="Family homology", citation="method_choice",
         source_note="each replicate is a full pairwise alignment of ~2,500 "
                     "residues, about 0.1 s; 20 gives a null standard "
                     "deviation stable to the third decimal, which is finer "
                     "than the effects being judged",
         description="Composition-matched shuffles used to place each "
                     "alignment statistic against its own null."),
    dict(key="homology.min_z", name="Minimum z over the null", value=3.0,
         unit="sigma", minimum=1.0, maximum=10.0, kind="method",
         category="Family homology", citation="method_choice",
         source_note="a statistic within 3 sigma of what a shuffled sequence "
                     "of the same composition gives is reported as not "
                     "distinguishable from chance. Deliberately not a p-value: "
                     "the null here is a sampling distribution of a score, not "
                     "a hypothesis test on data",
         description="How far above its shuffled null a statistic must sit "
                     "before it is called signal."),
    dict(key="homology.site_window", name="Site reliability window", value=101,
         unit="residues", minimum=5, maximum=401, kind="method",
         category="Family homology", citation="method_choice",
         source_note="chosen by a power scan, not by taste. The reliability "
                     "test needs enough columns to separate a real local "
                     "alignment from a chance one, and at the family's actual "
                     "divergence a narrow window has none: measured on the "
                     "transmembrane gate, z runs 1.3-2.6 at width 31 for "
                     "every non-mammalian member — refusing a mapping that is "
                     "visibly right — and reaches 3.3-4.5 at 101 for the "
                     "worm, fly and plant while Dictyostelium stays at 2.7-2.9 "
                     "and is correctly still refused. Wider trades locality "
                     "for power: at 101 the test answers 'is this region in "
                     "register', not 'is this residue', and the docstring "
                     "says so",
         description="Aligned columns centred on a mapped position whose "
                     "substitution score decides whether the mapping can be "
                     "trusted."),
]

#: Building a trimer from one protomer. Kept beside the family parameters
#: because the reason the feature exists is a family one: the only structural
#: representation of a non-animal PIEZO is a monomer, and everything measuring
#: a dome needs three protomers.
HOMOLOGY_PARAMETERS += [
    dict(key="assembly.min_corresponding", name="Residues to place a protomer",
         value=200, unit="residues", minimum=20, maximum=2000, kind="method",
         category="Family homology", citation="method_choice",
         source_note="a Kabsch fit will happily superpose twenty atoms and "
                     "report a small RMSD while placing the fold anywhere; "
                     "200 is roughly a 4-TM unit, below which the placement "
                     "is not constrained by the architecture. A monomer that "
                     "cannot reach it against a template is refused rather "
                     "than assembled loosely",
         description="Corresponding C-alphas required before a protomer may "
                     "be placed on a template chain."),
    dict(key="assembly.clash_distance", name="Inter-protomer clash cutoff",
         value=2.5, unit="A", minimum=1.0, maximum=5.0, kind="convention",
         category="Family homology", citation="method_choice",
         source_note="below the ~3.0-3.4 A of a real heavy-atom contact and "
                     "above a hydrogen bond's 2.7-3.2 A donor-acceptor "
                     "distance, so it counts atoms that are interpenetrating "
                     "rather than touching. Nothing models the interface, so "
                     "this count is the honest symptom of a template that "
                     "does not fit the protein being assembled",
         description="Heavy-atom separation below which two atoms in "
                     "different assembled protomers count as clashing."),
]

HOMOLOGY_PARAMETERS += [
    dict(key="assembly.core_fraction", name="Core descent fraction",
         value=0.6, unit="fraction", minimum=0.1, maximum=0.95, kind="method",
         category="Family homology", citation="method_choice",
         source_note="how fast the search shrinks while no residue is yet "
                     "within the cutoff. Only a descent rate — the core "
                     "itself is defined by `assembly.core_cutoff`, and this "
                     "exists because starting from a 19 A global fit nothing "
                     "is within 3 A and a distance criterion alone never "
                     "starts",
         description="Fraction of the current core kept per descent cycle."),
    dict(key="assembly.refit_cycles", name="Core-fit rejection cycles",
         value=20, unit="", minimum=1, maximum=100, kind="method",
         category="Family homology", citation="method_choice",
         source_note="superpose, drop the worst deviations, refit. The loop "
                     "exits on its own when the kept set stops changing, so "
                     "this is a runaway bound rather than a tuned depth",
         description="Outlier-rejection rounds when finding the rigid core."),
    dict(key="assembly.core_cutoff", name="Core-fit outlier cutoff",
         value=3.0, unit="A", minimum=0.5, maximum=15.0, kind="method",
         category="Family homology", citation="method_choice",
         source_note="a residue deviating more than this after superposition "
                     "is not following the template and is dropped. A "
                     "distance rather than a fraction, because a fraction "
                     "drove the core to its floor on every catalogue entry — "
                     "200 of 2,500 residues fitted to 1.2 A, which is not a "
                     "core but the 200 that agree best, and always exists. "
                     "With a cutoff the surviving count is a measurement of "
                     "how much of the protomer the template accounts for",
         description="Post-superposition deviation above which a residue is "
                     "excluded from the rigid core."),
]
