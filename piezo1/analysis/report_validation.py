"""Report entries that **validate** the model rather than measure with it.

Split from ``report_tags.py`` at the length limit and along a real seam. Every
other analysis in the registry asks the model a question — how wide is the
pore, where would a tag sit, what current flows. These two ask whether the
model is right about something it did not choose: the B-factors the entry was
deposited with, and a paralogue it was never fitted to.

That difference shows in what they return. Both carry a control that needs none
of the machinery being validated — contact number for the fluctuation, a
shuffled site correspondence for the paralogue — because a validation without
one measures whether two things are both about proteins.
"""

from __future__ import annotations

from ..core.structure import Structure

__all__ = ["analysis_fluctuations", "analysis_paralogue",
           "analysis_homology"]


def analysis_fluctuations(st: Structure, species: str,
                          n_modes: int | None = None, **kw) -> dict:
    """The elastic network against this entry's own B-factors.

    Reported with the control beside it, because a residue with many
    neighbours moves less whether or not there is a normal mode anywhere near
    it. An entry whose column cannot be interpreted returns the reason instead
    of a number.
    """
    from .fluctuations import assess_b_factors, compare_fluctuations

    quality = assess_b_factors(st)
    out = {"n_residues_in_column": quality.n_residues,
           "n_distinct_values": quality.n_distinct,
           "distinct_fraction": quality.distinct_fraction,
           "b_range": [quality.minimum, quality.maximum],
           "column_is_plddt": quality.is_confidence,
           "usable": quality.usable, "reason": quality.reason}
    if not quality.usable:
        out["note"] = ("no comparison was made; the reason above is the "
                       "result, not a failure to compute")
        return out

    result = compare_fluctuations(st, n_modes=n_modes)
    if not result.available:
        out.update({"usable": False, "reason": result.quality.reason})
        return out
    out.update({
        "n_residues_compared": len(result.residues),
        "n_modes": result.n_modes,
        "pearson": result.pearson_r,
        "spearman": result.spearman_r,
        "control_contact_number_pearson": result.control_pearson,
        "control_contact_number_spearman": result.control_spearman,
        "beats_control": result.beats_control,
        "control_inverted": result.control_inverted,
        "by_mode_count": result.by_mode_count,
        "note": ("Spearman is the number to read: the relationship is "
                 "monotone but not linear, so Pearson is dominated by a few "
                 "very mobile residues. A negative control means this entry's "
                 "B-factor rises with burial, which no mobility does — there "
                 "the column is the problem, not the network.")})
    return out


def analysis_paralogue(st: Structure, species: str, **kw) -> dict:
    """PIEZO1 against PIEZO2 — the only control on whether this is the fold.

    Runs on the **loaded** structure when it is a PIEZO1 entry, against 6KG7.
    The dome is reported twice, naively and coverage-matched, because the
    difference between those two rows is the result.
    """
    from .paralogue import compare, identify_numbering

    identity = identify_numbering(st)
    if identity.is_piezo2:
        return {"error": "this is a PIEZO2 entry; load a PIEZO1 structure and "
                         "the comparison runs against the best PIEZO2 partner "
                         "available — the same species where there is one"}
    result = compare(piezo1_pdb=st.name)
    if "error" in result:
        return result

    dome, modes = result["dome"], result["modes"]
    return {
        "piezo1": result["piezo1"], "piezo2": result["piezo2"],
        "tm_helices_agreeing_by_index":
            f"{result['tm_correspondence']['n_agree']} of "
            f"{result['tm_correspondence']['n_helices']}",
        "sequence_identity": result["tm_correspondence"]["identity"],
        "resolved_tm_helices": {"piezo1": dome["n_helices_piezo1"],
                                "piezo2": dome["n_helices_piezo2"],
                                "shared": len(dome["shared_helices"])},
        "dome_naive": [row.summary() for row in dome["naive"]],
        "dome_coverage_matched": [row.summary()
                                  for row in dome["coverage_matched"]],
        "gating_mode_overlap": modes.best_overlap,
        "gating_mode_in_piezo2_symmetric_subspace":
            modes.symmetric_subspace_overlap,
        "shuffled_control": modes.shuffled_control,
        "beats_control": modes.beats_control,
        "n_matched_sites": modes.n_sites,
        "protomer_order": list(modes.protomer_order),
        "superposition_rmsd_A": modes.superposition_rmsd,
        "note": ("the naive dome rows differ mostly in how much blade each "
                 "entry resolves; on the coverage-matched rows the two "
                 "proteins are indistinguishable, and PIEZO1's candidate "
                 "gating mode is present in PIEZO2 — this mechanism is a "
                 "property of the fold rather than of PIEZO1"),
    }


def analysis_homology(st: Structure, species: str, **kw) -> dict:
    """The loaded entry against every PIEZO the catalogue can reach.

    The widest form of the question ``analysis_paralogue`` asks of PIEZO2
    alone. It belongs in this file for the same reason that one does: it tests
    the model against something it was never fitted to — and it is the only
    entry in the registry whose answer can *fall*, since a mode overlap that
    survives one duplication need not survive 800 Myr.

    **It reports a range per protein, not one comparison, and that is the
    result.** Run against a single partner it produced 0.980 for dPIEZO from
    7WLT and **0.189** for the same two proteins from 8YEZ. Which entry is
    picked decides the number, so a representative pair is a cherry-pick and
    the range is the measurement. PIEZO2 comes back 0.80-0.98 across six pairs
    with every one beating its control, which is the positive control that
    makes the invertebrate spread an instability rather than a broken
    instrument.

    Runs no sequence matrix: that costs two minutes for the nulls and is a
    property of the family rather than of the loaded structure, so it lives at
    ``piezo1.cli homology`` where it is asked for deliberately.
    """
    from ..core.numbering_check import identify_numbering
    from .homology_structure import comparable_entries, mode_overlap_spread

    identity = identify_numbering(st)
    partners = {protein: entries
                for protein, entries in comparable_entries().items()
                if protein != identity.protein and entries}
    if not partners:
        return {"error": f"no catalogued homologue of {identity.protein} with "
                         f"three protomers is downloaded"}

    rows = []
    for protein, entries in partners.items():
        spread = mode_overlap_spread([st.name], entries, protein)
        if not spread.n_pairs:
            continue
        rows.append({
            "protein": protein, "n_entries": spread.n_pairs,
            "entries": [b for _a, b in spread.pairs],
            "gating_mode_overlap_low": round(spread.low, 3),
            "gating_mode_overlap_high": round(spread.high, 3),
            "shuffled_control_max": round(max(spread.controls), 3),
            "n_beating_control": spread.n_beating_control,
            "stable": spread.stable,
            "verdict": spread.summary(),
        })
    return {
        "loaded": {"pdb": st.name, "protein": identity.protein,
                   "numbering": identity.reference},
        "homologues": rows,
        "caveat": ("Read the RANGE, not a number. The overlap depends on which "
                   "deposited entry each side uses — across the catalogue the "
                   "same two proteins give anywhere from 0.18 to 0.98 for the "
                   "invertebrates — so a single pair is not a property of the "
                   "proteins. Sequence identity is not the column to read "
                   "either: below 0.30 it is inside Rost's twilight zone and "
                   "barely above what a shuffled sequence gives (see "
                   "docs/HOMOLOGY_SEARCH.md)."),
        "note": ("PIEZO2 is the positive control and it is stable, so the "
                 "instrument can say yes. Where the invertebrate range is "
                 "wide, that is a real limit on what these structures "
                 "support — a result, not a failure."),
    }
