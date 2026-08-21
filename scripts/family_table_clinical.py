"""The census findings that rest on a patient record, a coordinate file or a genome.

Split from ``family_table.py`` at its length limit and along a real seam: every
statement there is a comparison of amino acids between species, and every
statement here needed something else — ClinVar and UniProt disease annotation,
a superposition, or a 194-genome sweep.

Same shape as ``family_table.SEQUENCE_FINDINGS``; see that module for the field
meanings and for how ``check`` entries are verified on every build.
"""

from __future__ import annotations

__all__ = ["CLINICAL_FINDINGS", "EQUIVALENT_POSITIONS", "PORE_MODULE_PATHOGENIC",
           "CENSUS"]


def _row(**kw) -> dict:
    return kw


#: The fourteen pathogenic positions inside the pore module, in the numbering of
#: the gene each belongs to, with the zebrafish piezo3 residue aligned to each.
#: Every one is identical in piezo3 — which is the finding, and also the reason
#: the list is carried rather than a count: a count cannot be drawn on a
#: structure, and these can.
#:
#: ``msa_col`` is the column of the census project's 117-sequence family
#: alignment. Two columns carry a position from *both* disease genes; those are
#: the equivalent pairs below.
PORE_MODULE_PATHOGENIC = (
    # gene,    resi, aa, element,      msa_col, piezo3_resi, piezo3_aa
    ("PIEZO1", 2070, "K", "outer_helix", 8802, 2192, "K"),
    ("PIEZO1", 2088, "R", "pore_linker", 8827, 2210, "R"),
    ("PIEZO1", 2117, "E", "pore_linker", 8881, 2239, "E"),
    ("PIEZO1", 2127, "T", "pore_linker", 8891, 2249, "T"),
    ("PIEZO1", 2225, "M", "CED", 9214, 2347, "M"),
    ("PIEZO1", 2430, "P", "CED", 9808, 2560, "P"),
    ("PIEZO1", 2456, "R", "inner_helix", 9902, 2586, "R"),
    ("PIEZO1", 2488, "R", "CTD", 9935, 2618, "R"),
    ("PIEZO1", 2496, "E", "CTD", 10001, 2626, "E"),
    ("PIEZO2", 2295, "W", "outer_helix", 8798, 2188, "W"),
    ("PIEZO2", 2406, "G", "pore_linker", 9122, 2299, "G"),
    ("PIEZO2", 2686, "R", "inner_helix", 9902, 2586, "R"),
    ("PIEZO2", 2718, "R", "CTD", 9935, 2618, "R"),
    ("PIEZO2", 2739, "S", "CTD", 10063, 2639, "S"),
)

#: The two alignment columns that carry a pathogenic position in both disease
#: genes: the same residue of the same machine, mutated to cause two different
#: human diseases. Kept as an explicit pair because it is the one census result
#: that can be *seen* on a structure, and this project holds coordinates for
#: both proteins.
EQUIVALENT_POSITIONS = (
    {
        "msa_col": 9902,
        "element": "inner_helix",
        "piezo1": 2456, "piezo1_aa": "R", "piezo1_disease": "DHS1 (hereditary xerocytosis)",
        "piezo2": 2686, "piezo2_aa": "R", "piezo2_disease": "DA3 Gordon / MWKS",
        "piezo3": 2586, "piezo3_aa": "R",
        "note": "Both gain-of-function; PIEZO1 R2456H is the best-characterised "
                "slow-inactivating variant in the family and 8YFG resolves it.",
    },
    {
        "msa_col": 9935,
        "element": "CTD",
        "piezo1": 2488, "piezo1_aa": "R", "piezo1_disease": "DHS1 (hereditary xerocytosis)",
        "piezo2": 2718, "piezo2_aa": "R", "piezo2_disease": "DA5 distal arthrogryposis",
        "piezo3": 2618, "piezo3_aa": "R",
        "note": "Both in the C-terminal intracellular domain, the most "
                "identical element between the paralogues (91%).",
    },
)

#: The eukaryote-wide census, as the four numbers a structural project can use:
#: what the family's range actually is, so that this project's nine reviewed
#: members can be described as a sample of something with a measured size.
CENSUS = {
    "n_proteins": 8329,
    "by_kingdom": {"Metazoa": 6592, "Viridiplantae": 1397,
                   "protists": 318, "Fungi": 7},
    #: The four kingdom counts sum to 8,314 against a stated total of 8,329.
    #: Recorded rather than reconciled: it is the census's arithmetic, the gap
    #: is 0.2%, and quietly adjusting either number here would make this copy
    #: disagree with the project it is a copy of. A test asserts the gap is
    #: small and stays where it is.
    "unassigned_to_kingdom": 15,
    "archaeal_proteomes_searched": 628,
    "archaeal_hits": 0,
    "fungal_proteomes_searched": 1501,
    "fungal_genuine": 1,
    "fungal_genuine_species": "Rozella allomycis",
    "vertebrate_genomes": 194,
    "note": (
        "PIEZO is a eukaryotic invention, inherited almost universally, and its "
        "absences split into the real (fungi after the Rozella split, "
        "intracellular parasites) and the bibliographic (arthropod and mollusc "
        "reference proteomes whose genes are in the genome but not the "
        "proteome). The bacterial statement is a 300-genus sample, not a census."
    ),
}

CLINICAL_FINDINGS = (
    {
        "key": "disease_in_the_pore_module",
        "session": "S24",
        "kind": "clinical",
        "title": "Human disease concentrates in the part of the family that never changed",
        "statement": (
            "Of all labelled human missense positions in PIEZO1 and PIEZO2, the "
            "pore module - 17-18% of the protein - holds 14 of 43 pathogenic "
            "positions (33%) and 16 of 146 benign ones (11%): odds ratio 3.9, "
            "one-sided Fisher P = 0.0014. Pathogenic variation piles into the "
            "region the family conserved for half a billion years; benign "
            "variation, if anything, avoids it."
        ),
        "numbers": {
            "n_pathogenic_in_module": 14, "n_pathogenic": 43,
            "n_benign_in_module": 16, "n_benign": 146,
            "odds_ratio": 3.922, "p_fisher": 0.00135,
        },
        "source": "alignments/structure_stats.tsv",
        "check": [
            ("table", _row(quantity="pore_module_pathogenic"), "value", 14),
            ("table", _row(quantity="pore_module_benign"), "value", 16),
            ("table", _row(quantity="pore_module_odds_ratio"), "value", 3.922),
            ("table", _row(quantity="pore_module_p_fisher"), "value", 0.00135),
        ],
        "here": "analysis.disease_geography.pore_module_enrichment",
        "caveat": (
            "Variants are found where people look. The two genes are sequenced "
            "for different clinical reasons, and the benign set is dominated by "
            "population data whose ascertainment differs from the pathogenic "
            "set's. The test says where the labelled positions sit, not where "
            "disease-causing change is possible."
        ),
    },
    {
        "key": "same_residue_two_genes",
        "session": "S24",
        "kind": "clinical",
        "title": "Two disease genes mutate the same residue",
        "statement": (
            "The fourteen pathogenic pore-module positions occupy only twelve "
            "alignment columns. PIEZO1 R2456, behind a human red-cell disorder, "
            "and PIEZO2 R2686, behind Gordon syndrome, are the same residue of "
            "the same machine in two different genes - as are PIEZO1 R2488 and "
            "PIEZO2 R2718. The paralogues have been separate for half a billion "
            "years and disease still finds the same place in both."
        ),
        "numbers": {"n_positions": 14, "n_columns": 12, "n_shared": 2},
        "source": "alignments/alignment_stats.tsv",
        "check": [
            ("table", _row(alignment="pore_variants", scope="pathogenic_positions"),
             "n_sequences", 14),
            ("table", _row(alignment="pore_variants", scope="pathogenic_positions"),
             "n_columns", 12),
        ],
        "here": "analysis.equivalent_positions.locate",
        "caveat": (
            "An alignment column is a claim about correspondence, not a "
            "measurement of one. Whether the two arginines occupy the same "
            "place in space is a separate question - which this project can "
            "answer, and does, by superposing the two channels."
        ),
    },
    {
        "key": "piezo3_kept_disease_residues",
        "session": "S17",
        "kind": "clinical",
        "title": "The residues that make people ill are the residues piezo3 kept",
        "statement": (
            "Carrying labelled human positions across to piezo3: at positions "
            "where a human mutation causes disease piezo3 has the same amino "
            "acid 75% of the time (PIEZO1) and 70% (PIEZO2); at benign "
            "positions 31% and 26%; across aligned columns as a whole 45%. "
            "Disease positions sit above background, benign ones below, and "
            "positions of unknown significance land exactly on it. Of the 14 "
            "pathogenic positions inside the pore module, piezo3 has the "
            "identical residue at all 14."
        ),
        "numbers": {
            "piezo1_pathogenic_frac": 0.75, "piezo1_benign_frac": 0.3068,
            "piezo1_background_frac": 0.4529, "piezo1_odds_ratio": 3.624,
            "piezo1_p_fisher": 0.00612, "pore_module_identical": 14,
        },
        "source": "constraint/piezo3_audit_summary.tsv",
        "check": [
            ("table", _row(gene="PIEZO1", class_bucket="P/LP"),
             "frac_identical", 0.75),
            ("table", _row(gene="PIEZO1", class_bucket="B/LB"),
             "frac_identical", 0.3068),
            ("table", _row(gene="PIEZO1", class_bucket="VUS"),
             "frac_identical", 0.4356),
        ],
        "here": "analysis.piezo3.kept_positions",
        "caveat": (
            "No functional test has been done. The gene is transcribed, "
            "spliced and tissue-patterned, but nobody has measured whether it "
            "conducts. These are predictions about a channel nobody has "
            "recorded from."
        ),
    },
    {
        "key": "constraint_predicts_pathogenicity",
        "session": "S17",
        "kind": "clinical",
        "title": "The per-position constraint separates pathogenic from benign at AUC 0.91",
        "statement": (
            "Because the clinical labels are independent of the evolutionary "
            "scores, they also test them. Used as a classifier on human PIEZO1, "
            "the deep constraint track separates pathogenic from benign "
            "missense at AUC 0.914 (21 vs 107 positions, Mann-Whitney "
            "P = 1.1e-9)."
        ),
        "numbers": {"auc": 0.9143, "n_pathogenic": 21, "n_benign": 107,
                    "p_mannwhitney": 1.07e-09},
        "source": "constraint/variant_constraint_test.tsv",
        "check": [
            ("table", _row(gene="PIEZO1", layer="deep"), "auc", 0.9143),
            ("table", _row(gene="PIEZO1", layer="deep"), "n_pathogenic", 21),
        ],
        "here": "analysis.disease_geography.constraint_classifier",
        "caveat": (
            "The score recovers judgements clinicians already made; it has not "
            "been shown to make new ones correctly. The pathogenic set is "
            "small and non-randomly located."
        ),
    },
    {
        "key": "pore_module_superposition",
        "session": "S24",
        "kind": "structure",
        "title": "The pore module is where the channels agree in space, not just in sequence",
        "statement": (
            "Superposing the predicted zebrafish piezo3 model on the "
            "experimental cryo-EM mouse Piezo1 structure 6B3R *by the pore "
            "module alone* puts the cores at 3.86 A RMSD over 448 matched "
            "C-alpha while the blades splay visibly apart. The core-conserved / "
            "periphery-free result stops being a profile and becomes a shape."
        ),
        "numbers": {"rmsd": 3.86, "n_matched": 448},
        "source": "alignments/structure_stats.tsv",
        "check": [
            ("table", _row(quantity="piezo3_vs_6B3R_pore_rmsd"), "value", 3.86),
            ("table", _row(quantity="piezo3_vs_6B3R_pore_rmsd"), "total", 448),
        ],
        "here": "analysis.core_periphery.compare",
        "caveat": (
            "One of the two structures is a prediction, and a prediction of a "
            "single protomer at that. A low core RMSD between a model and an "
            "experiment is partly a statement about what the predictor was "
            "trained on."
        ),
    },
    {
        "key": "family_range",
        "session": "S20",
        "kind": "census",
        "title": "PIEZO is a eukaryotic invention, and its absences split in two",
        "statement": (
            "8,329 PIEZO-like proteins across Metazoa, Viridiplantae, protists "
            "and Fungi; zero hits in 628 archaeal reference proteomes and a "
            "300-genus bacterial sample under every protocol tried. Among 1,501 "
            "fungal reference proteomes the only genuine PIEZO is Rozella "
            "allomycis, the deepest fungal branch - everything after that "
            "split, Dikarya included, lacks the family."
        ),
        "numbers": {"n_proteins": 8329, "archaeal_hits": 0,
                    "fungal_proteomes": 1501},
        "source": "FINDINGS.md",
        "check": [("document", "8,329"), ("document", "1,501"),
                  ("document", "628")],
        "here": "recorded only - this project holds ten family references, of "
                "which nine are reviewed, and does not sample genomes. The "
                "census is what says how large a sample that is",
        "caveat": (
            "The bacterial statement is a sample, not a census. Iterative "
            "profile search was measured to break at scale on the two largest "
            "databases, so the metazoan and protist totals are profile-sweep "
            "figures with 69 qualified candidates held out."
        ),
    },
    {
        "key": "piezo3_is_ancient_and_lost",
        "session": "S13/S15",
        "kind": "census",
        "title": "piezo3 is as old as its siblings, and humans lost it before the primates",
        "statement": (
            "Both duplications creating PIEZO1, PIEZO2 and piezo3 sit at or "
            "before the origin of jawed vertebrates, roughly 460-560 Ma. piezo3 "
            "has since been lost eleven times; the human remnant is the "
            "catalogued pseudogene PIEZO1P2, still holding 24% of the protein "
            "in recognisable pieces, and the loss reconstructs to the primate "
            "stem 88-74 Ma."
        ),
        "numbers": {"duplication_ma_low": 460, "duplication_ma_high": 560,
                    "n_losses": 11, "human_remnant_frac": 0.24},
        "source": "FINDINGS.md",
        "check": [("document", "PIEZO1P2"), ("document", "88 and 74 million")],
        "here": "recorded only - no piezo3 coordinates exist for any mammal",
        "caveat": (
            "The loss count depends on the question: 5 if asked whether the "
            "locus is still there, 16 if asked whether a working gene is. "
            "Eleven sits between them, and the sensitivity matrix runs 4-22."
        ),
    },
    {
        "key": "pfew_not_a_signature",
        "session": "review 2026-08-19",
        "kind": "sequence",
        "title": "A motif quoted as a family signature is not one",
        "statement": (
            "The PFEW motif reported as absolutely conserved across protozoan "
            "PIEZO homologues in the 2013 survey does not occur in human "
            "PIEZO1, nor in any of the 117 representative sequences spanning "
            "protists to primates. What is conserved that deeply are three "
            "short windows around the pore."
        ),
        "numbers": {"n_sequences_searched": 117, "n_occurrences": 0},
        "source": "FINDINGS.md",
        "check": [("document", "PFEW")],
        "here": "analysis.family_motifs.motif_scan",
        "caveat": (
            "A negative about one literal four-residue string. It does not "
            "show the 2013 survey was wrong about its own sequences, only that "
            "the motif does not generalise to the family as now sampled."
        ),
    },
    {
        "key": "codon_selection_agrees",
        "session": "S17",
        "kind": "sequence",
        "title": "Counting codons rather than amino acids gives the same answer",
        "statement": (
            "Site-wise dN/dS (FEL) over the trimmed codon alignment puts the "
            "pore region as the most protected part of the gene and the "
            "pore-lining helix as the most protected element in it - the same "
            "ordering the amino-acid comparison gives, from a statistic that "
            "shares nothing with it. Across all three genes there is exactly "
            "one site anywhere with a hint of positive selection, against more "
            "than three thousand under purifying selection."
        ),
        "numbers": {"piezo1_anchor_frac_purifying": 0.6306,
                    "piezo1_whole_frac_purifying": 0.4949},
        "source": "constraint/selection_by_domain.tsv",
        "check": [
            ("table", _row(paralog="PIEZO1", domain="anchor"),
             "frac_purifying_q05", 0.6306),
            ("table", _row(paralog="PIEZO1", domain="WHOLE_PROTEIN"),
             "frac_purifying_q05", 0.4949),
        ],
        "here": "analysis.family_constraint.selection_track",
        "caveat": (
            "FEL was run on 57 codon sequences against 121-192 protein "
            "sequences, so its resolution is coarser and its coverage lower "
            "(70-74% of positions carry a usable rate)."
        ),
    },
)
