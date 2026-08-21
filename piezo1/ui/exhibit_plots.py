"""Charts drawn from the result already in the window.

**Nothing here recomputes anything.** Each builder takes the dict the analysis
returned — the same object the table beside it is formatted from — and turns
part of it into a picture. That is the one rule this module has, and it is the
rule the drawn overlays already follow: ``pore_controller`` reads the analysis
object rather than re-running the profiler, because a picture and a table of
two different runs is the kind of disagreement nobody notices.

The second rule follows from the first: **a reference line comes from the
parameter registry or from the result, never from a literal typed here.** A
number written into a chart would be a fourth home for it, beside the code, the
documentation and the claims registry — which is the reason ``tour.py`` computes
every number it states instead of quoting one.

A builder that cannot draw returns an empty chart carrying the reason. Missing
keys are ordinary: a shut pore has no current, an entry with no partner has no
comparison, and saying so is better than a flat line at zero.
"""

from __future__ import annotations

from .exhibits import ChartData, Reference, Series, empty_chart

__all__ = ["BUILDERS", "build_chart", "PALETTE", "bar_series"]

#: Kept in step with `profile_plot`, which is the other chart in the
#: application; two plotting palettes would read as two kinds of thing.
PALETTE = ("#6fb1ff", "#f2a65a", "#7fd18a", "#c678dd", "#e06c75", "#8a919e")


def _p(key: str) -> float:
    from ..parameters import PARAMETERS

    return float(PARAMETERS.value(key))


def bar_series(categories, values, name, color=PALETTE[0], axis=0) -> Series:
    """One bar series over the categories a chart already declares."""
    return Series(name=name, x=list(range(len(values))),
                  y=[float(v) for v in values], kind="bar", color=color,
                  axis=axis)


# --------------------------------------------------------------------------
# Ion permeation
# --------------------------------------------------------------------------

def permeation_conductance(data: dict) -> ChartData:
    if not data.get("conducting"):
        return empty_chart(
            "This pore does not conduct in the model, so there is no "
            f"conductance to draw. {data.get('blocked_by', 'no reason given')}. "
            "Every mechanism found is listed in the table; the simulation "
            "beside this shows what the radius would have to be.")
    values, names = [data["conductance_pS"]], ["drift-diffusion"]
    if "independent_check_pS" in data:
        values.append(data["independent_check_pS"])
        names.append("series-resistor check")
    if "calcium_2mM_pS" in data:
        values.append(data["calcium_2mM_pS"])
        names.append("with 2 mM calcium")
    chart = ChartData(
        title="Unitary conductance", y_label="pS", categories=names,
        series=[bar_series(names, values, "this entry")],
        references=[Reference(_p("permeation.published_conductance"),
                              "published single-channel conductance")],
        note="The two inputs that dominate this number — in-pore diffusivity "
             "and ion radius — are unmeasured, and span 16-94 pS between "
             "them. The series-resistor bar is a closed-form check that shares "
             "no code with the solver, not a second measurement.")
    return chart


# --------------------------------------------------------------------------
# Elastic network against the deposited B-factors
# --------------------------------------------------------------------------

def fluctuation_controls(data: dict) -> ChartData:
    if not data.get("usable", True):
        return empty_chart(f"This entry's B-factor column cannot answer the "
                           f"question: {data.get('reason', 'no reason given')}")
    if "pearson" not in data:
        return empty_chart("no correlation in the result to draw")
    names = ["Pearson", "Spearman"]
    return ChartData(
        title="Elastic network against the deposited B-factors",
        y_label="correlation", categories=names,
        series=[
            bar_series(names, [data["pearson"], data["spearman"]],
                  "elastic network", PALETTE[0]),
            bar_series(names, [data.get("control_contact_number_pearson", 0.0),
                          data.get("control_contact_number_spearman", 0.0)],
                  "contact-number control", PALETTE[5]),
        ],
        references=[Reference(0.0, "no relationship")],
        note="The control counts neighbours and uses no network at all. Where "
             "it wins, the agreement was burial rather than mechanism — which "
             "is why both statistics are shown and not the better one.")


# --------------------------------------------------------------------------
# PIEZO2 and the wider family
# --------------------------------------------------------------------------

def paralogue_overlap(data: dict) -> ChartData:
    if "gating_mode_overlap" not in data:
        return empty_chart("no mode comparison in the result")
    names = ["gating mode overlap", "in symmetric subspace", "shuffled control"]
    values = [data["gating_mode_overlap"],
              data.get("gating_mode_in_piezo2_symmetric_subspace", 0.0),
              data.get("shuffled_control", 0.0)]
    return ChartData(
        title="PIEZO1's gating mode in PIEZO2", y_label="overlap",
        categories=names,
        series=[bar_series(names, values, "overlap", PALETTE[0])],
        note="The control shuffles the correspondence and recomputes. Without "
             "it an overlap of 0.8 has nothing to be large compared with.")


def paralogue_coverage(data: dict) -> ChartData:
    resolved = data.get("resolved_tm_helices") or {}
    if not resolved:
        return empty_chart("no helix counts in the result")
    names = list(resolved)
    return ChartData(
        title="Transmembrane helices each entry resolves",
        y_label="helices", categories=names,
        series=[bar_series(names, [resolved[k] for k in names], "resolved",
                      PALETTE[1])],
        note="This is the whole of the naive dome difference. The "
             "coverage-matched comparison uses the shared set, which is the "
             "only one both entries model.")


def homology_range(data: dict) -> ChartData:
    partners = data.get("homologues") or []
    if not partners:
        return empty_chart("no partner entries were comparable")
    names = [p["protein"] for p in partners]
    return ChartData(
        title="Gating-mode overlap across the family", y_label="overlap",
        categories=names,
        series=[
            bar_series(names, [p.get("gating_mode_overlap_low", 0.0) for p in partners],
                  "lowest pair", PALETTE[5]),
            bar_series(names, [p.get("gating_mode_overlap_high", 0.0) for p in partners],
                  "highest pair", PALETTE[0]),
            bar_series(names, [p.get("shuffled_control_max", 0.0) for p in partners],
                  "shuffled control", PALETTE[4]),
        ],
        note="Each partner is a different deposited entry with its own "
             "coverage, so the bars are not comparable across proteins. The "
             "gap between the lowest and highest pair is why this is reported "
             "as a range.")


# --------------------------------------------------------------------------
# Contacts
# --------------------------------------------------------------------------

def interaction_counts(data: dict) -> ChartData:
    counts = data.get("counts") or {}
    if not counts:
        return empty_chart("no contacts in the result")
    names = list(counts)
    return ChartData(
        title="Contact inventory", y_label="count (log)", categories=names,
        series=[bar_series(names, [counts[k] for k in names], "contacts",
                      PALETTE[2])],
        log_y=True,
        note="Log axis: hydrogen bonds outnumber disulfides by three orders of "
             "magnitude. Most hydrogen bonds are backbone i to i+4 — the helix "
             "the cartoon is already drawing.")


# --------------------------------------------------------------------------
# HaloTag labelling and geometry
# --------------------------------------------------------------------------

def labelling_histogram(data: dict) -> ChartData:
    protocol = data.get("dye_histogram_at_protocol") or []
    if not protocol:
        return empty_chart("no occupancy distribution in the result")
    names = [f"{i} dyes" for i in range(len(protocol))]
    series = [bar_series(names, protocol, "at the protocol", PALETTE[0])]
    reduced = (data.get("if_90_percent_reactive") or {}).get("dye_histogram")
    if reduced:
        series.append(bar_series(names, reduced, "if 90% of tags are reactive",
                            PALETTE[1]))
    return ChartData(
        title="Dyes per channel", y_label="fraction of channels",
        categories=names, series=series,
        note="Binomial(3, p) over the three sites. The second series is the "
             "whole argument for unreactive tags: it is a mixture, not a "
             "shifted mean.")


def fusion_geometry(data: dict) -> ChartData:
    if "tag_centre_to_pore_exit_nm" not in data:
        return empty_chart("no fusion geometry in the result")
    low, high = (data.get("envelope_range_nm") or [None, None])[:2]
    names, values = [], []
    for label, value in (("tag to pore exit", data["tag_centre_to_pore_exit_nm"]),
                         ("envelope nearest", low),
                         ("envelope median", data.get("envelope_median_nm")),
                         ("envelope furthest", high),
                         ("tag radius of gyration",
                          (data.get("tag_radius_A") or 0.0) / 10.0)):
        if value is not None:
            names.append(label)
            values.append(value)
    return ChartData(
        title="Where the tag can sit", y_label="nm", categories=names,
        series=[bar_series(names, values, "modelled", PALETTE[1])],
        note="A region, not a pose. The envelope is every position the linker "
             "admits without the tag clashing into the channel; nothing "
             "chooses between them.")


def hybrid_split(data: dict) -> ChartData:
    if "seam_rmsd_A" not in data:
        return empty_chart("no graft in the result")
    names = ["seam fit", "global disagreement"]
    return ChartData(
        title="The graft, in Angstrom", y_label="RMSD (A)", categories=names,
        series=[bar_series(names, [data["seam_rmsd_A"], data.get("global_rmsd_A", 0.0)],
                      "experiment against prediction", PALETTE[3])],
        note=(f"{data.get('predicted_residues', 0)} residues are prediction, "
              f"{data.get('confident_fraction', 0.0):.0%} of them above pLDDT "
              f"70. The seam is anchored locally, which is why it fits; the "
              f"second bar is what a good local fit hides."))


# --------------------------------------------------------------------------
# Modulators
# --------------------------------------------------------------------------

def ligand_potency(_data: dict) -> ChartData:
    """Read from the curated resource rather than the formatted result.

    The result formats potency as prose ("EC50 26.6 uM (...)"), and parsing a
    sentence back into a number is how a plot ends up disagreeing with the
    table it sits beside. The resource is where both come from.
    """
    from ..core.ligands import load_ligands

    ligands = load_ligands()
    names, values, labels = [], [], []
    for item in ligands.ligands:
        if not item.potency:
            continue
        micromolar = float(item.potency["value"])
        if item.potency.get("unit") == "nM":
            micromolar /= 1000.0
        names.append(item.name)
        values.append(micromolar)
        labels.append(f"{item.name}: {item.potency['measure']}")
    if not values:
        return empty_chart("no modulator carries a measured potency")
    missing = [x.name for x in ligands.ligands if not x.potency]
    return ChartData(
        title="Measured potency", y_label="uM (log)", categories=names,
        series=[bar_series(names, values, "potency", PALETTE[4])], log_y=True,
        note=("Not one quantity: " + "; ".join(labels) + ". "
              + (f"No potency is curated for {', '.join(missing)}. "
                 if missing else "")
              + "A smaller bar is a more potent compound."))


# --------------------------------------------------------------------------
# What the deposited variant structures support
# --------------------------------------------------------------------------

def variant_survey(data: dict) -> ChartData:
    entries = data.get("entries") or []
    usable = [e for e in entries if e.get("bottleneck_A") is not None]
    if not usable:
        return empty_chart("no variant entry yielded a pore profile")
    names = [f"{e['pdb']} {e.get('variant') or ''}".strip() for e in usable]
    series = [bar_series(names, [e["bottleneck_A"] for e in usable],
                    "bottleneck radius (A)", PALETTE[0])]
    scores = [e.get("wetting_score") for e in usable]
    if all(s is not None for s in scores):
        series.append(bar_series(names, scores, "wetting score", PALETTE[1], axis=1))
    duplicates = [f"{e['pdb']} = {', '.join(e['duplicates'])}"
                  for e in usable if e.get("duplicates")]
    return ChartData(
        title="Every deposited variant entry", y_label="bottleneck (A)",
        categories=names, series=series,
        references=[Reference(_p("hydration.water_radius") * 10.0,
                              "radius of a water molecule")],
        note=("Entries sharing identical coordinates: "
              + ("; ".join(duplicates) if duplicates else "none")
              + ". A bottleneck below the water radius is shut on geometry "
                "alone, before any hydrophobicity is considered."))


def _paired_metric(data: dict, metric: str, unit: str) -> ChartData:
    block = (data.get("metrics") or {}).get(metric)
    if not block:
        return empty_chart(f"no {metric} in the result")
    low, high = block["wild_type_range"][:2]
    return ChartData(
        title=f"{data.get('variant', 'variant')} against wild type",
        y_label=unit, categories=[data.get("variant", "variant")],
        series=[bar_series([metric], [block["variant"]], metric, PALETTE[4])],
        references=[Reference(low, "wild-type range", high=high)],
        note=("inside the wild-type range"
              if block.get("within_wild_type_range") else
              "outside the wild-type range")
             + f"; three independent wild-type entries span {high - low:.3g} "
               f"{unit}, and the largest difference to any of them is "
               f"{block.get('largest_variant_difference', float('nan')):.3g}.")


def paired_bottleneck(data: dict) -> ChartData:
    return _paired_metric(data, "bottleneck_A", "A")


def paired_wetting(data: dict) -> ChartData:
    return _paired_metric(data, "wetting_score", "score")


# --------------------------------------------------------------------------

BUILDERS = {
    "permeation_conductance": permeation_conductance,
    "fluctuation_controls": fluctuation_controls,
    "paralogue_overlap": paralogue_overlap,
    "paralogue_coverage": paralogue_coverage,
    "homology_range": homology_range,
    "interaction_counts": interaction_counts,
    "labelling_histogram": labelling_histogram,
    "fusion_geometry": fusion_geometry,
    "hybrid_split": hybrid_split,
    "ligand_potency": ligand_potency,
    "variant_survey": variant_survey,
    "paired_bottleneck": paired_bottleneck,
    "paired_wetting": paired_wetting,
}


def build_chart(name: str, data: dict) -> ChartData:
    """Dispatch, and never raise into the window.

    A chart that cannot be drawn is a panel saying why. A traceback out of a
    paint path takes the application with it, which is the reason
    ``TourStep.report`` catches too.
    """
    from .exhibit_plots_family import FAMILY_BUILDERS

    builder = BUILDERS.get(name) or FAMILY_BUILDERS.get(name)
    if builder is None:
        return empty_chart(f"no chart builder named {name!r}")
    try:
        return builder(data or {})
    except Exception as exc:                      # noqa: BLE001 — see docstring
        return empty_chart(f"could not draw this chart: "
                           f"{type(exc).__name__}: {exc}")
