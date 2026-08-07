#!/usr/bin/env python
"""Execute the Round 48 test exactly as pre-registered.

Runs only what ``docs/PREREGISTRATION_ROUND48.md`` specifies, in the order it
specifies, and writes the result whatever it shows.

Usage::

    python scripts/run_validation_round48.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.design import benjamini_hochberg  # noqa: E402
from piezo1.analysis.features import build_feature_table  # noqa: E402
from piezo1.analysis.validation import (auroc, bootstrap_cliffs_delta,  # noqa: E402
                                        cliffs_delta, permutation_test)
from piezo1.analysis.variant_sets import build_analysis_set  # noqa: E402
from piezo1.config import DERIVED_DIR, STRUCTURE_DIR  # noqa: E402
from piezo1.core.sequence import human_to_mouse  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402

#: Pre-registered §5.3 — the only position carrying both directions, and the
#: only one carrying more than one variant at all.
EXCLUDED_POSITIONS = {2456}

#: Pre-registered §4 — the reference structure, frozen.
REFERENCE = "7WLT"

#: Pre-registered §7. ``alternative`` is stated on the LoF sample: "less" means
#: H1 puts LoF *below* GoF. The negative control has no proposed direction and
#: is therefore two-sided.
PRIMARY = ("relative_sasa", "less", "LoF positions more buried")
SECONDARY = [
    ("conservation", "greater", "LoF positions more conserved"),
    ("prs_gate_response", "less", "GoF positions more coupled to the gate"),
    ("gating_amplitude", "less", "GoF positions move more along the gating mode"),
    ("distance_to_gate", "greater", "GoF positions closer to the gate"),
    ("distance_to_axis", "two-sided", "NEGATIVE CONTROL — no direction proposed"),
]

#: Pre-registered §6 — below this the round is exploratory, whatever p says.
CONFIRMATORY_THRESHOLD = 0.495


def positions(levels) -> dict:
    """One direction per position, per §5.2, with §5.3 applied."""
    entries = build_analysis_set(levels=levels).missense().entries
    by_position: dict[int, set] = {}
    for entry in entries:
        by_position.setdefault(entry.residue, set()).add(entry.classification)
    return {p: next(iter(c)) for p, c in by_position.items()
            if len(c) == 1 and p not in EXCLUDED_POSITIONS}


def feature_values(table, column: str, labels: dict) -> dict:
    """Column value at each labelled position, converted human → mouse.

    Never a constant offset: the curated variants are human-numbered and 7WLT is
    mouse, so every position goes through the alignment map.
    """
    by_residue = table.as_dict(column)
    out = {}
    for human in labels:
        mouse = human_to_mouse(human)
        if mouse is not None and mouse in by_residue:
            out[human] = float(by_residue[mouse])
    return out


def test(values: dict, labels: dict, name: str, alternative: str,
         hypothesis: str) -> dict | None:
    """Cliff's δ for LoF vs GoF, in the pre-registered direction."""
    lof = np.array([values[p] for p, d in labels.items()
                    if d == "LoF" and np.isfinite(values.get(p, np.nan))])
    gof = np.array([values[p] for p, d in labels.items()
                    if d == "GoF" and np.isfinite(values.get(p, np.nan))])
    if len(lof) < 3 or len(gof) < 3:
        return None
    permutation = permutation_test(lof, gof, alternative=alternative)
    effect = bootstrap_cliffs_delta(lof, gof)
    scores = np.concatenate([lof, gof])
    positive = np.concatenate([np.zeros(len(lof), bool), np.ones(len(gof), bool)])
    return {"endpoint": name, "hypothesis": hypothesis,
            "alternative": alternative,
            "n_lof": int(len(lof)), "n_gof": int(len(gof)),
            "median_lof": float(np.median(lof)), "median_gof": float(np.median(gof)),
            "cliffs_delta": float(cliffs_delta(lof, gof)),
            "ci_low": float(effect.ci_low), "ci_high": float(effect.ci_high),
            "p_value": float(permutation.p_value),
            "auroc": float(auroc(scores, positive))}


def variance_share(table, column: str) -> dict:
    """The §2 ceiling, measured rather than asserted.

    The roadmap required the between-position variance share to be reported
    alongside any result. For a feature computed on the wild-type structure it
    is 1.0 exactly, and the demonstration is R2456: four variants, three GoF and
    one LoF, all receiving the identical value.
    """
    import numpy as np

    from piezo1.analysis.substitution import variance_decomposition

    values = table.as_dict(column)
    entries = build_analysis_set(
        levels=("measured", "disease_mechanism")).missense().entries
    positions_, values_, shared = [], [], {}
    for entry in entries:
        mouse = human_to_mouse(entry.residue)
        if mouse in values:
            positions_.append(entry.residue)
            values_.append(values[mouse])
            shared.setdefault(entry.residue, []).append(
                (entry.label, entry.classification, values[mouse]))
    split = variance_decomposition(np.array(positions_), np.array(values_))
    multi = {p: v for p, v in shared.items() if len(v) > 1}
    return {"column": column, "n_variants": len(positions_),
            "n_positions": len(set(positions_)),
            "between_fraction": float(split.between_fraction),
            "within_fraction": float(split.within_fraction),
            "multi_variant_positions": {str(p): v for p, v in multi.items()}}


def main() -> int:
    path = STRUCTURE_DIR / f"{REFERENCE}.cif"
    if not path.exists():
        print(f"{REFERENCE} not downloaded; run python -m piezo1.io.fetch")
        return 1

    combined = positions(("measured", "disease_mechanism"))
    measured_only = positions(("measured",))

    print(f"Building the feature table on {REFERENCE} (this takes a few minutes)")
    table = build_feature_table(Structure.from_file(path))

    # §5.4 — positions not modelled in the reference are excluded, and counted.
    primary_values = feature_values(table, PRIMARY[0], combined)
    modelled = {p: d for p, d in combined.items() if p in primary_values}
    unmodelled = sorted(set(combined) - set(modelled))
    n_gof = sum(1 for d in modelled.values() if d == "GoF")
    print(f"\nPositions: {len(combined)} directional, "
          f"{len(modelled)} modelled ({n_gof} GoF, {len(modelled) - n_gof} LoF)")
    print(f"Excluded, not modelled in {REFERENCE}: {len(unmodelled)} {unmodelled}")
    print(f"Excluded, both directions: {sorted(EXCLUDED_POSITIONS)}")

    # ------------------------------------------------------------- primary
    column, alternative, hypothesis = PRIMARY
    primary = test(primary_values, modelled, column, alternative, hypothesis)
    print(f"\nPRIMARY — {column}: {hypothesis}")
    print(f"  median LoF {primary['median_lof']:.4f} vs "
          f"GoF {primary['median_gof']:.4f}")
    print(f"  Cliff's delta {primary['cliffs_delta']:+.3f}  "
          f"95% CI [{primary['ci_low']:+.3f}, {primary['ci_high']:+.3f}]")
    print(f"  p = {primary['p_value']:.4f}   AUROC {primary['auroc']:.3f}   "
          f"n = {primary['n_lof']} LoF vs {primary['n_gof']} GoF")

    # §8 — all three clauses, fixed in advance.
    clauses = {"p < 0.05": primary["p_value"] < 0.05,
               "delta negative": primary["cliffs_delta"] < 0,
               "CI excludes zero": not (primary["ci_low"] <= 0 <= primary["ci_high"])}
    reject = all(clauses.values())
    print("\n  decision clauses (all three required):")
    for name, passed in clauses.items():
        print(f"    {name:20s} {'pass' if passed else 'FAIL'}")
    print(f"  >>> {'REJECT H0' if reject else 'FAIL TO REJECT H0'}")

    confirmatory = abs(primary["cliffs_delta"]) >= CONFIRMATORY_THRESHOLD
    print(f"  >>> {'CONFIRMATORY' if confirmatory else 'EXPLORATORY'} "
          f"(|delta| {abs(primary['cliffs_delta']):.3f} vs the "
          f"{CONFIRMATORY_THRESHOLD} threshold fixed in section 6)")

    # ----------------------------------------------------------- secondary
    family = []
    for column, alternative, hypothesis in SECONDARY:
        values = feature_values(table, column, modelled)
        entry = test(values, modelled, column, alternative, hypothesis)
        if entry is not None:
            family.append(entry)
    values = feature_values(table, PRIMARY[0], measured_only)
    entry = test(values, {p: d for p, d in measured_only.items() if p in values},
                 f"{PRIMARY[0]} (measured only)", PRIMARY[1],
                 "the primary on measured labels only")
    if entry is not None:
        family.append(entry)

    # ---- the ceiling from §2, measured (roadmap: report the variance share)
    share = variance_share(table, PRIMARY[0])
    print(f"\nWITHIN-POSITION VARIANCE SHARE of {PRIMARY[0]}: "
          f"{share['within_fraction']:.6f} "
          f"({share['n_variants']} variants at {share['n_positions']} positions)")
    for position, group in share["multi_variant_positions"].items():
        print(f"  position {position} carries {len(group)} variants, all valued "
              f"{group[0][2]:.6f}:")
        for label, direction, _ in group:
            print(f"    {label} ({direction})")

    correction = benjamini_hochberg([f["p_value"] for f in family], alpha=0.05)
    print(f"\nSECONDARY FAMILY ({len(family)} tests, BH q=0.05)")
    print(f"  {'endpoint':34s} {'n':>9s} {'delta':>7s} {'p':>7s} {'q':>7s}  sig")
    for item, q, rejected in zip(family, correction.adjusted, correction.rejected):
        item["q_value"], item["significant"] = float(q), bool(rejected)
        print(f"  {item['endpoint']:34s} {item['n_lof']:4d}/{item['n_gof']:<4d} "
              f"{item['cliffs_delta']:+7.3f} {item['p_value']:7.4f} "
              f"{q:7.4f}  {'YES' if rejected else 'no'}")

    out = {"preregistration": "docs/PREREGISTRATION_ROUND48.md",
           "reference_structure": REFERENCE,
           "n_positions_directional": len(combined),
           "n_positions_modelled": len(modelled),
           "unmodelled_positions": unmodelled,
           "excluded_positions": sorted(EXCLUDED_POSITIONS),
           "variance_share": share,
           "primary": primary,
           "decision": "reject" if reject else "fail_to_reject",
           "decision_clauses": clauses,
           "status": "confirmatory" if confirmatory else "exploratory",
           "confirmatory_threshold": CONFIRMATORY_THRESHOLD,
           "secondary": family,
           "any_secondary_significant": bool(any(f["significant"] for f in family))}
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    path = DERIVED_DIR / "validation_round48.json"
    path.write_text(json.dumps(out, indent=1, default=float))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
