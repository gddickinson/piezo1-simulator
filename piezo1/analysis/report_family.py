"""Report entries for the imported PIEZO-family census.

Split from ``report_validation.py`` at that file's length limit and along a real
seam. Those entries validate this project's model against something it did not
choose. These do something different again: they take an **external project's
result** and ask what this project's coordinates and physics say about it.

That difference decides what every entry here returns. Each carries three things
a normal analysis does not need — the source it came from, what the census
concluded, and what changes when the question is asked here — because the
failure mode is not a wrong number, it is a reader taking an imported statement
for one this project measured.
"""

from __future__ import annotations

from ..core.structure import Structure

__all__ = ["analysis_family", "analysis_constraint", "analysis_disease",
           "analysis_core_periphery", "analysis_piezo3"]


def analysis_family(st: Structure, species: str, **kw) -> dict:
    """Every imported finding, with what this project does with each."""
    from ..core.family import load_family_findings

    findings = load_family_findings()
    return {
        "source": findings.source,
        "source_note": findings.provenance["source_note"],
        "verified": findings.provenance["verified"],
        "n_findings": len(findings.findings),
        "findings": [
            {"key": f.key, "session": f.session, "kind": f.kind,
             "title": f.title, "statement": f.statement,
             "explored_here": f.here, "caveat": f.caveat,
             "numbers_verified": f.n_checks}
            for f in findings.findings],
        "census": findings.census,
        "caveat": ("every statement here is the census project's, re-verified "
                   "against its own source files at import time and not "
                   "measured by this project. The entries below measure."),
    }


def analysis_constraint(st: Structure, species: str, **kw) -> dict:
    """Evolutionary constraint on this entry, and on our own domain partition."""
    from .family_constraint import (blade_gradient, census_domain_constraint,
                                    compare_with_own_conservation,
                                    constraint_on_structure, domain_constraint)

    placed = constraint_on_structure(st)
    if not placed:
        return {"error": getattr(placed, "reason", "could not be placed"),
                "caveat": "the constraint track is human PIEZO1's and cannot "
                          "be read on an entry in another numbering"}

    domains = domain_constraint()
    ordered = sorted((d for d in domains if d.mean is not None),
                     key=lambda d: -d.mean)
    blades = blade_gradient()
    cross = compare_with_own_conservation()
    return {
        "numbering": placed.numbering,
        "residues_scored": placed.n_residues_scored,
        "coverage": round(placed.coverage, 4),
        "converted_through_alignment": placed.converted,
        "most_constrained_domains": [
            {"domain": d.domain, "category": d.category, "mean": round(d.mean, 4),
             "vs_whole_protein": round(d.vs_whole, 4), "n": d.n_residues}
            for d in ordered[:5]],
        "least_constrained_domains": [
            {"domain": d.domain, "category": d.category, "mean": round(d.mean, 4),
             "vs_whole_protein": round(d.vs_whole, 4), "n": d.n_residues}
            for d in ordered[-3:]],
        "census_partition": census_domain_constraint(),
        "blade_gradient": {k: v for k, v in blades.items() if k != "per_thu"},
        "per_thu": blades["per_thu"],
        "own_conservation_cross_check": (
            None if cross is None else
            {"n": cross.n, "spearman": round(cross.spearman, 4),
             "pearson": round(cross.pearson, 4),
             "n_orthologues_census": cross.n_orthologues_census,
             "n_orthologues_ours": cross.n_orthologues_ours,
             "note": cross.note}),
        "caveat": ("the per-residue values are the census's; the domain "
                   "partition, the THU split and the cross-check are this "
                   "project's. Where the two partitions disagree, the "
                   "disagreement is the result and both are printed."),
    }


def analysis_disease(st: Structure, species: str, **kw) -> dict:
    """Does pathogenic variation concentrate in the pore module? Re-tested here."""
    from .disease_geography import both_partitions, constraint_classifier

    both = both_partitions()
    classifier = constraint_classifier()
    return {
        "results": {name: {
            "region": r.region, "region_fraction": round(r.region_fraction, 4),
            "pathogenic_in": r.pathogenic_in, "pathogenic_out": r.pathogenic_out,
            "comparator_in": r.comparator_in, "comparator_out": r.comparator_out,
            "odds_ratio": r.odds_ratio, "p_fisher": r.p_value,
            "significant": r.significant, "summary": r.summary()}
            for name, r in both["results"].items()},
        "verdict": both["verdict"],
        "boundary_disagreement": both["boundaries"],
        "disputed_band": both["disputed"],
        "census": both["census"],
        "constraint_classifier": (
            None if classifier is None else
            {"auc_here": round(classifier.auc, 4),
             "auc_census": classifier.census_auc,
             "n_positive": classifier.n_positive,
             "n_negative": classifier.n_negative,
             "agrees_with_census": classifier.agrees_with_census,
             "note": classifier.note}),
        "caveat": ("PIEZO1 only against a population comparator, where the "
                   "census pooled PIEZO1 and PIEZO2 against ClinVar benign "
                   "labels. An independent re-test, not a reproduction."),
    }


def analysis_core_periphery(st: Structure, species: str,
                            partner: str = "6KG7", **kw) -> dict:
    """Fit this entry on a partner by the pore module, then measure the blades."""
    from ..config import STRUCTURE_DIR
    from .core_periphery import compare
    from .equivalent_positions import locate

    path = STRUCTURE_DIR / f"{partner}.cif"
    if not path.exists():
        return {"error": f"{partner} is not downloaded"}
    other = Structure.from_file(path)
    result = compare(st, other, "loaded", partner)
    if not result:
        return {"error": result.reason}

    out = {
        "partner": partner,
        "n_core": result.n_core, "core_rmsd_A": round(result.core_rmsd, 3),
        "n_periphery": result.n_periphery,
        "periphery_rmsd_A": (None if result.periphery_rmsd is None
                             else round(result.periphery_rmsd, 3)),
        "splay_ratio": (None if result.splay_ratio is None
                        else round(result.splay_ratio, 3)),
        "whole_rmsd_A": round(result.whole_rmsd, 3),
        "core_converged": result.core_converged,
        "cross_paralogue": result.cross_paralogue,
        "correspondence": result.note,
        "summary": result.summary(),
    }
    if result.cross_paralogue:
        equivalence = locate(st, other, "loaded", partner)
        out["equivalent_positions"] = {
            "verdict": equivalence.verdict,
            "control_median_A": (None if equivalence.control_median is None
                                 else round(equivalence.control_median, 3)),
            "pairs": [{"label": p.label, "element": p.element,
                       "distance_A": (None if p.distance is None
                                      else round(p.distance, 3)),
                       "register_offset": p.register_offset,
                       "alignment_agrees": p.alignment_agrees,
                       "same_place": p.same_place}
                      for p in equivalence.pairs],
        }
    out["caveat"] = (
        "a core-only fit is directional: it asks where the blades land given "
        "that the pores are superposed, and it can fail. A prediction against "
        "an experiment splays far more than two experiments do, even for the "
        "same protein, so a large splay against a model says more about the "
        "model than about the two proteins.")
    return out


def analysis_piezo3(st: Structure, species: str, **kw) -> dict:
    """The third vertebrate PIEZO, run through this project's own pipeline."""
    from .piezo3 import (best_paralogue_template, fold_comparison,
                         kept_positions, template_survey)
    from .piezo3_channel import build_channel

    try:
        kept = kept_positions()
        fold = fold_comparison()
        fits = template_survey()
        channel = build_channel(best_paralogue_template(fits))
    except FileNotFoundError as exc:
        return {"error": str(exc)}

    return {
        "kept_pathogenic_positions": {
            "n": len(kept),
            "n_kept": sum(k.kept for k in kept),
            "n_agreeing_with_census_record": sum(k.agrees_with_census for k in kept),
            "positions": [
                {"gene": k.gene, "human": f"{k.human_aa}{k.human_resi}",
                 "element": k.element, "piezo3_model_resi": k.model_resi,
                 "piezo3_aa": k.model_aa, "kept": k.kept}
                for k in kept],
        },
        "fold_comparison": (fold.summary() if fold else fold.reason),
        "templates": [f.summary() for f in fits],
        "channel": {
            "template": channel.template,
            "template_protein": channel.template_protein,
            "identity": round(channel.identity, 3),
            "clashes": channel.clashes,
            "borrowed_fraction": channel.borrowed,
            "piezo3": channel.piezo3.summary(),
            "comparison": (None if channel.comparison is None
                           else channel.comparison.summary()),
            "conducts": channel.piezo3.conducts,
            "wetting_score": channel.piezo3.wetting_score,
        },
        "verdict": channel.verdict,
        "caveats": list(channel.caveats),
    }
