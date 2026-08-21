"""What is imported from the ``piezo_genes`` census project, and what it claims.

Authored content, in the idiom of ``parameter_table.py`` and ``ligand_table.py``:
the table is data so the whole import can be read and diffed without running the
validation machinery in :mod:`build_family_findings`.

**Why a table of claims rather than a copy of the files.** The census project
produced tens of megabytes of intermediate results; almost none of it is a
finding. What travels here is the small set of statements a structural project
can *use* or *test*, each one carrying the number it rests on and a selector
into the source file that produced it. The build re-reads every number from the
source and refuses to write if one has moved — so this file cannot quietly
become a stale paraphrase of a project that has since been corrected.

Each finding declares:

``key``           stable identifier, used by the loader and the tests
``session``       the census session that established it (S6 … S24)
``kind``          ``sequence`` | ``clinical`` | ``structure`` | ``census``
``statement``     the finding in one sentence, in the census project's own terms
``numbers``       the values claimed, verified against ``source`` on every build
``source``        file under the census project's ``results/``
``check``         how to re-read each number: (selector, column)
``here``          what this project does with it — a module path, or the reason
                  it can only be recorded
``caveat``        what the finding does not establish
"""

from __future__ import annotations

__all__ = ["FINDINGS", "SEQUENCE_FINDINGS", "CONSTRAINT_GENES",
           "TABLE_IMPORTS", "SOURCE_PROJECT"]

#: The census project this is imported from. Its results directory is the only
#: thing read; nothing under its data root (bulk genomes, proteomes, structures)
#: is required, and the build says so when it is missing.
SOURCE_PROJECT = {
    "name": "piezo_genes",
    "title": "PIEZO family gene census",
    "default_path": "../piezo_genes",
    "results_subdir": "results",
    "note": (
        "A 194-genome, eukaryote-wide census of the PIEZO family: what the "
        "family's true range is, that vertebrates have a third paralogue "
        "(piezo3) the databases largely missed, and which parts of the protein "
        "half a billion years of evolution has refused to change."
    ),
}

#: Per-residue constraint tracks imported in full, one per gene. ``reference``
#: is the sequence the residue numbers belong to — this is the join to every
#: annotation in this project, so it must be an accession we already hold or
#: can state the numbering of.
CONSTRAINT_GENES = (
    {
        "gene": "PIEZO1",
        "accession": "Q92508",
        "numbering": "human",
        "length": 2521,
        "source": "constraint/constraint_PIEZO1_Q92508.tsv",
        "n_orthologues": 174,
    },
    {
        "gene": "PIEZO2",
        "accession": "Q9H5I5",
        "numbering": "human_piezo2",
        "length": 2752,
        "source": "constraint/constraint_PIEZO2_Q9H5I5.tsv",
        "n_orthologues": 192,
    },
    {
        "gene": "piezo3",
        "accession": "A0AC58GFC9",
        "numbering": "zebrafish_piezo3",
        "length": 2652,
        "source": "constraint/constraint_piezo3_A0AC58GFC9.tsv",
        "n_orthologues": 121,
    },
)

#: Whole tables copied across, because every row is used rather than one number.
#: ``key`` names them in the resource; ``columns`` is the subset kept.
TABLE_IMPORTS = (
    {
        "key": "domain_map",
        "source": "constraint/domain_map.tsv",
        "columns": ("paralog", "domain", "mouse_start", "mouse_end",
                    "ref_start", "ref_end", "n_residues", "aligned_frac", "note"),
        "note": "The census project's own domain partition, in each gene's own "
                "numbering. Kept so its per-domain numbers can be read on the "
                "boundaries that produced them, never to replace ours.",
    },
    {
        "key": "constraint_by_domain",
        "source": "constraint/constraint_by_domain.tsv",
        "columns": ("paralog", "domain", "n_sites", "mean_jsd", "median_jsd",
                    "vs_whole_protein", "p_greater_than_distal_blades"),
        "note": "Within-paralogue constraint, domain by domain.",
    },
    {
        "key": "paralogue_identity",
        "source": "constraint/paralog_identity_by_domain.tsv",
        "columns": ("pair", "frame", "domain", "n_columns", "identity",
                    "whole_protein_identity", "delta_vs_whole"),
        "note": "What the three paralogues kept of each other, domain by domain.",
    },
    {
        "key": "variant_constraint_auc",
        "source": "constraint/variant_constraint_test.tsv",
        "columns": ("gene", "layer", "n_pathogenic", "n_benign",
                    "mean_jsd_pathogenic", "mean_jsd_benign", "auc",
                    "p_mannwhitney"),
        "note": "The constraint score used as a variant classifier.",
    },
    {
        "key": "piezo3_audit",
        "source": "constraint/piezo3_audit_summary.tsv",
        "columns": ("gene", "contrast", "class_bucket", "n_positions",
                    "n_identical_in_piezo3", "frac_identical",
                    "comparator_frac_identical", "comparator_n",
                    "odds_ratio", "p_fisher"),
        "note": "Whether piezo3 kept the human residue, split by clinical class.",
    },
    {
        "key": "selection_by_domain",
        "source": "constraint/selection_by_domain.tsv",
        "columns": ("paralog", "domain", "n_sites", "median_beta", "mean_beta",
                    "frac_purifying_q05", "median_omega_usable"),
        "note": "The same question asked in codons rather than amino acids "
                "(FEL site-wise dN/dS), which is why it is worth carrying: it "
                "shares no statistic with the JSD track.",
    },
)


def _row(**kw) -> dict:
    """A selector: the row of the source table whose columns equal ``kw``."""
    return kw


#: The sequence-level statements. The clinical, structural and census ones
#: live in ``family_table_clinical.py``, split at this file's length limit and
#: along a real seam: everything here is a comparison of amino acids, and
#: everything there rests on a patient record, a coordinate file or a genome.
#:
#: Every statement carries the number it rests on.
#:
#: ``check`` entries are verified on every build. ``table`` checks re-read a
#: cell from a source TSV; ``document`` checks require the number to appear in
#: the census project's own ``FINDINGS.md``, which is how the claims that were
#: never reduced to a table (the census totals, the loss counts, the dates) are
#: held to the same standard as the ones that were.
SEQUENCE_FINDINGS = (
    {
        "key": "core_is_pore",
        "session": "S17",
        "kind": "sequence",
        "title": "The conserved core is the pore, in every paralogue independently",
        "statement": (
            "Compared against 121-192 orthologues of itself, each of the three "
            "vertebrate PIEZO genes is least constrained across the blades and "
            "most constrained across the pore machinery - the anchor, the two "
            "pore helices and the C-terminal domain. The ordering is the same "
            "computed separately in PIEZO1, PIEZO2 and piezo3."
        ),
        "numbers": {
            "piezo1_distal_blade_mean_jsd": 0.6561,
            "piezo1_pore_linker_mean_jsd": 0.8219,
            "piezo1_inner_helix_mean_jsd": 0.8074,
            "piezo1_pore_module_p": 8.13e-25,
        },
        "source": "constraint/constraint_by_domain.tsv",
        "check": [
            ("table", _row(paralog="PIEZO1", domain="distal_blades"),
             "mean_jsd", 0.6561),
            ("table", _row(paralog="PIEZO1", domain="pore_linker"),
             "mean_jsd", 0.8219),
            ("table", _row(paralog="PIEZO1", domain="inner_helix"),
             "mean_jsd", 0.8074),
        ],
        "here": "analysis.family_constraint.domain_constraint",
        "caveat": (
            "Measured on the census project's own domain partition, which is "
            "not this project's. Recomputed here on domains.json, which is a "
            "different boundary set and therefore a real test rather than a "
            "restatement."
        ),
    },
    {
        "key": "paralogues_kept_the_machine",
        "session": "S17",
        "kind": "sequence",
        "title": "Half a billion years rewrote the blades and left the pore",
        "statement": (
            "PIEZO1, PIEZO2 and piezo3 now agree at 45-48% of positions overall, "
            "but at 73-91% across the anchor, the pore helices, the pore linker "
            "and the CTD, and *below* the whole-protein figure across both blade "
            "regions - in all three pairwise comparisons."
        ),
        "numbers": {
            "piezo1_piezo2_whole": 0.4743,
            "piezo1_piezo2_inner_helix": 0.8519,
            "piezo1_piezo2_ctd": 0.9123,
            "piezo1_piezo3_whole": 0.4529,
            "piezo1_piezo3_outer_helix": 0.8519,
            "piezo1_piezo2_distal_blades": 0.4410,
        },
        "source": "constraint/paralog_identity_by_domain.tsv",
        "check": [
            ("table", _row(pair="PIEZO1_vs_PIEZO2", domain="inner_helix"),
             "identity", 0.8519),
            ("table", _row(pair="PIEZO1_vs_PIEZO2", domain="CTD"),
             "identity", 0.9123),
            ("table", _row(pair="PIEZO1_vs_piezo3", domain="outer_helix"),
             "identity", 0.8519),
            ("table", _row(pair="PIEZO1_vs_PIEZO2", domain="distal_blades"),
             "identity", 0.4410),
        ],
        "here": "analysis.family_constraint.paralogue_asymmetry",
        "caveat": (
            "An identity is a statement about sequence. Whether the two "
            "channels agree in *space* over the same region is a separate "
            "question, answered here by analysis.core_periphery."
        ),
    },
    {
        "key": "ced_is_the_exception",
        "session": "S17",
        "kind": "sequence",
        "title": "The lid is the one piece of the pore the paralogues did not keep",
        "statement": (
            "The CED - the extracellular cap over the pore - is the only part of "
            "the pore machinery that is *less* alike between paralogues than the "
            "protein as a whole (35-40%), and inside PIEZO2 it is less "
            "constrained than even the blades. The three genes share a pore and "
            "differ in its lid."
        ),
        "numbers": {
            "piezo1_piezo2_ced_identity": 0.4017,
            "piezo1_piezo3_ced_identity": 0.3777,
            "piezo2_piezo3_ced_identity": 0.3534,
            "piezo2_ced_mean_jsd": 0.5685,
            "piezo2_distal_blade_mean_jsd": 0.6601,
        },
        "source": "constraint/paralog_identity_by_domain.tsv",
        "check": [
            ("table", _row(pair="PIEZO1_vs_PIEZO2", domain="CED"),
             "identity", 0.4017),
            ("table", _row(pair="PIEZO2_vs_piezo3", domain="CED"),
             "identity", 0.3534),
        ],
        "here": "analysis.family_constraint.paralogue_asymmetry",
        "caveat": (
            "A warning against reading 'piezo3 has PIEZO's pore' as 'piezo3 "
            "behaves like PIEZO1'. In PIEZO1 and PIEZO2 the cap influences ion "
            "selectivity and inactivation, so this is where a functional "
            "difference is most likely to live."
        ),
    },
    {
        "key": "distal_beats_proximal",
        "session": "review 2026-08-19",
        "kind": "sequence",
        "title": "The distal blade is more conserved than the proximal blade",
        "statement": (
            "Within every paralogue the *distal* blade - the part furthest from "
            "the pore - carries higher constraint than the proximal blade. That "
            "is the opposite of what 'peripheral means dispensable' predicts, "
            "and it sits directly against the 2026 nematode structures reporting "
            "the pore-distal blade dispensable for mechanoactivation."
        ),
        "numbers": {
            "piezo1_distal": 0.6561,
            "piezo1_proximal": 0.5583,
            "piezo2_distal": 0.6601,
            "piezo2_proximal": 0.6330,
            "piezo3_distal": 0.6082,
            "piezo3_proximal": 0.5443,
        },
        "source": "constraint/constraint_by_domain.tsv",
        "check": [
            ("table", _row(paralog="PIEZO1", domain="proximal_blades"),
             "mean_jsd", 0.5583),
            ("table", _row(paralog="piezo3", domain="proximal_blades"),
             "mean_jsd", 0.5443),
        ],
        "here": "analysis.family_constraint.blade_gradient",
        "caveat": (
            "The census project's proximal_blades band includes the beam, and "
            "its distal band is much the larger of the two. The comparison is "
            "repeated here on this project's own THU boundaries."
        ),
    },
)


try:                                            # running as a script
    from family_table_clinical import CLINICAL_FINDINGS
except ImportError:                             # imported as scripts.family_table
    from .family_table_clinical import CLINICAL_FINDINGS

#: Everything imported, sequence findings first. The order is the order the
#: report and the CLI print them in.
FINDINGS = SEQUENCE_FINDINGS + CLINICAL_FINDINGS
