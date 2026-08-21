"""Charts for results this project imported, recorded or reproduced.

Split from :mod:`piezo1.ui.exhibit_plots` at the length limit and along the
catalogue's own seam. The rules are that module's: nothing is recomputed, and
no reference value is typed in here — it comes from the result or from the
parameter registry.

That second rule bites hardest on this half. The control ranges that make the
splay result readable, and the census's own published odds ratio, are recorded
in ``docs/FAMILY.md`` and pinned by tests. Where the analysis puts one in its
result, it is drawn; where it does not, the note points at the window's caveat
rather than repeating a number that would then have two homes.
"""

from __future__ import annotations

from .exhibit_plots import PALETTE, bar_series
from .exhibits import ChartData, Reference, empty_chart

__all__ = ["FAMILY_BUILDERS"]


# --------------------------------------------------------------------------
# The census itself
# --------------------------------------------------------------------------

def census_kingdoms(data: dict) -> ChartData:
    census = data.get("census") or {}
    by_kingdom = census.get("by_kingdom") or {}
    if not by_kingdom:
        return empty_chart("the imported census carries no kingdom breakdown")
    names = list(by_kingdom)
    total = census.get("n_proteins")
    counted = sum(by_kingdom.values()) + census.get("unassigned_to_kingdom", 0)
    gap = "" if total is None else (
        f" The four kingdoms plus the unassigned come to {counted:,} against a "
        f"stated total of {total:,}; the difference is recorded rather than "
        f"reconciled, because it is the census's arithmetic and not ours.")
    return ChartData(
        title="PIEZO proteins found, by kingdom", y_label="proteins (log)",
        categories=names,
        series=[bar_series(names, [by_kingdom[k] for k in names], "proteins",
                      PALETTE[2])],
        log_y=True,
        note=(f"{census.get('archaeal_proteomes_searched', 0):,} archaeal "
              f"proteomes searched, {census.get('archaeal_hits', 0)} hits."
              + gap))


def constraint_domains(data: dict) -> ChartData:
    most = data.get("most_constrained_domains") or []
    least = data.get("least_constrained_domains") or []
    if not most:
        return empty_chart("no per-domain constraint in the result")
    rows = list(most) + list(least)
    names = [f"{r['domain']} ({r['category']})" for r in rows]
    whole = rows[0]["mean"] - rows[0].get("vs_whole_protein", 0.0)
    return ChartData(
        title="Constraint by domain, on this project's boundaries",
        y_label="mean JSD", categories=names,
        series=[bar_series(names, [r["mean"] for r in rows], "mean constraint",
                      PALETTE[0])],
        references=[Reference(whole, "whole protein")],
        note=(f"{data.get('residues_scored', 0)} residues scored, coverage "
              f"{data.get('coverage', 0.0):.1%}"
              + (", carried through the alignment map because this entry is "
                 "not in the track's own numbering."
                 if data.get("converted_through_alignment") else ".")
              + " The values are the census's; the partition is ours."))


# --------------------------------------------------------------------------
# Where disease sits
# --------------------------------------------------------------------------

def disease_odds(data: dict) -> ChartData:
    results = data.get("results") or {}
    if not results:
        return empty_chart("no enrichment in the result")
    names, values, notes = [], [], []
    for key in ("census", "ours"):
        row = results.get(key)
        if not row:
            continue
        names.append(row["region"])
        values.append(row["odds_ratio"])
        notes.append(f"{row['region']}: P = {row['p_fisher']:.3g}"
                     + ("" if row["significant"] else " (not significant)"))
    census = (data.get("census") or {}).get("numbers") or {}
    if "odds_ratio" in census:
        names.append("the census's own result")
        values.append(census["odds_ratio"])
        notes.append(f"the census: P = {census.get('p_fisher', float('nan')):.3g}, "
                     f"two genes pooled against ClinVar benign labels")
    return ChartData(
        title="Pathogenic missense in the pore module", y_label="odds ratio",
        categories=names,
        series=[bar_series(names, values, "odds ratio", PALETTE[4])],
        references=[Reference(1.0, "no enrichment")],
        note="; ".join(notes) + ". The first two bars are one re-test under "
             "two partitions, not two findings.")


def disease_boundaries(data: dict) -> ChartData:
    rows = data.get("boundary_disagreement") or []
    if not rows:
        return empty_chart("the two partitions agree everywhere")
    names = [r["element"] for r in rows]
    band = data.get("disputed_band") or {}
    span = band.get("span") or []
    return ChartData(
        title="How far each boundary moves between the two partitions",
        y_label="residues", categories=names,
        series=[bar_series(names, [r["start_offset"] for r in rows], "start",
                      PALETTE[0]),
                bar_series(names, [r["end_offset"] for r in rows], "end",
                      PALETTE[1])],
        references=[Reference(0.0, "the partitions agree")],
        note=(f"The disputed band is {band.get('n_residues', 0)} residues"
              + (f" ({span[0]}-{span[1]})" if len(span) == 2 else "")
              + f" and carries {len(band.get('pathogenic') or [])} pathogenic "
                f"positions, which is what moves the odds ratio."))


# --------------------------------------------------------------------------
# Core and periphery
# --------------------------------------------------------------------------

def splay_ratio(data: dict) -> ChartData:
    if data.get("splay_ratio") is None:
        return empty_chart(
            "no splay ratio: the cores did not superpose, and a ratio against "
            "a core that did not fit is arithmetic rather than a measurement.")
    names = ["core", "periphery", "whole"]
    values = [data.get("core_rmsd_A"), data.get("periphery_rmsd_A"),
              data.get("whole_rmsd_A")]
    names = [n for n, v in zip(names, values) if v is not None]
    values = [v for v in values if v is not None]
    return ChartData(
        title=f"Fitted on the pore module alone: {data.get('partner', '')}",
        y_label="RMSD (A)", categories=names,
        series=[bar_series(names, values, "after a core-only fit", PALETTE[3])],
        note=(f"Splay ratio {data['splay_ratio']:.2f} over "
              f"{data.get('n_core', 0)} core and {data.get('n_periphery', 0)} "
              f"peripheral C-alpha. What that ratio means is in this window's "
              f"caveat: experimental cross-paralogue pairs and an AlphaFold "
              f"monomer of the same protein occupy different bands, and the "
              f"control is the finding."))


# --------------------------------------------------------------------------
# piezo3
# --------------------------------------------------------------------------

def piezo3_positions(data: dict) -> ChartData:
    kept = data.get("kept_pathogenic_positions") or {}
    if not kept:
        return empty_chart("no pathogenic-position table in the result")
    names = ["pathogenic positions", "identical in piezo3",
             "agreeing with the census"]
    values = [kept.get("n", 0), kept.get("n_kept", 0),
              kept.get("n_agreeing_with_census_record", 0)]
    elements = sorted({p["element"] for p in kept.get("positions", [])})
    return ChartData(
        title="The pore positions human disease strikes", y_label="positions",
        categories=names,
        series=[bar_series(names, values, "count", PALETTE[2])],
        note=(f"Elements involved: {', '.join(elements) or 'none recorded'}. "
              f"The third bar is a check against the imported record, not a "
              f"second measurement. Identity at a position is not function."))


def piezo3_template(data: dict) -> ChartData:
    channel = data.get("channel") or {}
    if "borrowed_fraction" not in channel:
        return empty_chart("no assembled trimer in the result")
    names = ["borrowed from the template", "template sequence identity"]
    values = [channel["borrowed_fraction"], channel.get("identity", 0.0)]
    return ChartData(
        title=f"Assembled on {channel.get('template', 'a template')} "
              f"({channel.get('template_protein', 'unknown')})",
        y_label="fraction", categories=names,
        series=[bar_series(names, values, "fraction", PALETTE[1])],
        references=[Reference(1.0, "all of it the template's")],
        note=(f"{channel.get('clashes', 0):,} inter-protomer clashes, against "
              f"3-8 in a real deposited trimer. At this borrowed fraction a "
              f"dome radius measured on the assembly is mostly a measurement "
              f"of the template."))


# --------------------------------------------------------------------------
# The prediction record, and the two papers
# --------------------------------------------------------------------------

def record_effects(data: dict) -> ChartData:
    tests = [t for t in (data.get("tests") or [])
             if t.get("cliffs_delta") is not None]
    if not tests:
        return empty_chart("no pre-registered test in the record carries an "
                           "effect size")
    names = [f"Round {t['round']}" for t in tests]
    rows = []
    for t in tests:
        p = t.get("p_value")
        rows.append(f"Round {t['round']}: "
                    + (f"p = {p:.3g}" if p is not None else
                       "no primary p recorded")
                    + f", n = {t.get('n_gof', 0)} vs {t.get('n_lof', 0)}")
    return ChartData(
        title=data.get("central_claim", "the central claim"),
        y_label="Cliff's delta", categories=names,
        series=[bar_series(names, [t["cliffs_delta"] for t in tests],
                      "measured effect", PALETTE[4])],
        references=[Reference(0.0, "no effect")],
        note="; ".join(rows)
             + ". Each was fixed in writing before it ran, and each failed to "
               "reject. A bar's direction is not a trend, and Round 22 "
               "records no primary p because its interval already spanned "
               "zero.")


def paper_coverage(data: dict) -> ChartData:
    """One builder for both papers: the two registries report the same census
    of panels under different key names, and duplicating the builder would let
    the two pictures drift apart."""
    counts = data.get("by_status") or data.get("coverage") or {}
    counts = {k: v for k, v in counts.items()
              if k != "total" and isinstance(v, (int, float))}
    if not counts:
        return empty_chart("no panel census in the result")
    paper = (data.get("paper") or {}).get("journal", "")
    names = [key.replace("_", " ") for key in counts]
    refused = data.get("cannot_replicate") or {}
    return ChartData(
        title=f"Panels, by what this project can do with them  {paper}".strip(),
        y_label="panels", categories=names,
        series=[bar_series(names, list(counts.values()), "panels", PALETTE[0])],
        note=(f"{len(refused)} panels carry a specific reason for being "
              f"refused — electrophysiology, a cryo-EM map, micrographs or a "
              f"trajectory. An analogue is a different quantity, not a weaker "
              f"reproduction."))


FAMILY_BUILDERS = {
    "census_kingdoms": census_kingdoms,
    "constraint_domains": constraint_domains,
    "disease_odds": disease_odds,
    "disease_boundaries": disease_boundaries,
    "splay_ratio": splay_ratio,
    "piezo3_positions": piezo3_positions,
    "piezo3_template": piezo3_template,
    "record_effects": record_effects,
    "paper_coverage": paper_coverage,
}
