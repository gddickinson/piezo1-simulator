#!/usr/bin/env python
"""Execute the Round 41 test exactly as pre-registered.

Runs only what ``docs/PREREGISTRATION_ROUND41.md`` specifies, in the order it
specifies, and writes the result whatever it shows.

Usage::

    python scripts/run_validation_round41.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.design import benjamini_hochberg  # noqa: E402
from piezo1.analysis.gnomad import (GnomadClient,  # noqa: E402
                                    missense_density)
from piezo1.analysis.validation import (auroc, bootstrap_cliffs_delta,  # noqa: E402
                                        cliffs_delta, permutation_test)
from piezo1.analysis.variant_sets import build_analysis_set  # noqa: E402
from piezo1.config import DERIVED_DIR  # noqa: E402

#: Pre-registered §5.3 — the only position carrying both directions.
EXCLUDED_POSITIONS = {2456}


def positions(levels) -> dict:
    """One direction per position, per §5.2, with §5.3 applied."""
    entries = build_analysis_set(levels=levels).missense().entries
    by_position: dict[int, set] = {}
    for entry in entries:
        by_position.setdefault(entry.residue, set()).add(entry.classification)
    return {p: next(iter(c)) for p, c in by_position.items()
            if len(c) == 1 and p not in EXCLUDED_POSITIONS}


def test(values: dict, labels: dict, name: str) -> dict | None:
    """Cliff's δ for LoF vs GoF, one-sided in the hypothesised direction.

    H1 says LoF positions are MORE constrained, i.e. a LOWER missense rate, so
    the test asks whether the LoF sample sits below the GoF one.
    """
    lof = np.array([values[p] for p, d in labels.items()
                    if d == "LoF" and np.isfinite(values.get(p, np.nan))])
    gof = np.array([values[p] for p, d in labels.items()
                    if d == "GoF" and np.isfinite(values.get(p, np.nan))])
    if len(lof) < 3 or len(gof) < 3:
        return None
    permutation = permutation_test(lof, gof, alternative="less")
    effect = bootstrap_cliffs_delta(lof, gof)
    scores = np.concatenate([lof, gof])
    positive = np.concatenate([np.zeros(len(lof), bool), np.ones(len(gof), bool)])
    return {"endpoint": name, "n_lof": int(len(lof)), "n_gof": int(len(gof)),
            "cliffs_delta": float(cliffs_delta(lof, gof)),
            "ci_low": float(effect.ci_low), "ci_high": float(effect.ci_high),
            "p_value": float(permutation.p_value),
            "auroc": float(auroc(scores, positive)),
            "median_lof": float(np.median(lof)),
            "median_gof": float(np.median(gof))}


def main() -> int:
    client = GnomadClient()
    constraint = client.constraint()
    variants = client.variants()
    if constraint is None or variants is None:
        raise SystemExit("gnomAD unavailable and not cached")

    print("GENE-LEVEL (recorded in the pre-registration, §2)")
    print(f"  {constraint.summary()}")

    combined = positions(("measured", "disease_mechanism"))
    measured = positions(("measured",))
    print(f"\nanalysis set: {len(combined)} positions "
          f"({sum(v == 'GoF' for v in combined.values())} GoF, "
          f"{sum(v == 'LoF' for v in combined.values())} LoF); "
          f"position 2456 excluded as pre-registered")

    density = missense_density(variants, window=25)
    rate = density.local_rate()
    values = {p: float(rate[p - 1]) for p in combined}

    primary = test(values, combined, "local missense rate, ±25, combined")
    if primary is None:
        raise SystemExit("primary endpoint has too few positions")
    print("\nPRIMARY — local missense rate (±25), LoF vs GoF positions")
    print(f"  n = {primary['n_lof']} LoF vs {primary['n_gof']} GoF")
    print(f"  median  LoF {primary['median_lof']:.3f} | "
          f"GoF {primary['median_gof']:.3f} missense/residue")
    print(f"  Cliff's δ {primary['cliffs_delta']:+.3f} "
          f"[{primary['ci_low']:+.3f}, {primary['ci_high']:+.3f}]")
    print(f"  p = {primary['p_value']:.4f}   AUROC {primary['auroc']:.3f}")

    reject = (primary["p_value"] < 0.05 and primary["cliffs_delta"] < 0
              and primary["ci_high"] < 0)
    print(f"\n  DECISION: {'REJECT H0' if reject else 'FAIL TO REJECT H0'}")

    family = []
    measured_values = {p: float(rate[p - 1]) for p in measured}
    entry = test(measured_values, measured, "same, measured labels only")
    if entry:
        family.append(entry)

    allele = missense_density(variants, window=25).allele_count
    smoothed = np.convolve(allele, np.ones(51), mode="same") / 51.0
    entry = test({p: float(smoothed[p - 1]) for p in combined}, combined,
                 "summed allele count, ±25")
    if entry:
        family.append(entry)

    for window in (10, 50):
        other = missense_density(variants, window=window).local_rate()
        entry = test({p: float(other[p - 1]) for p in combined}, combined,
                     f"local missense rate, ±{window}")
        if entry:
            family.append(entry)

    entry = test({p: float(density.observed[p - 1]) for p in combined}, combined,
                 "raw per-residue count (negative control)")
    if entry:
        family.append(entry)

    correction = benjamini_hochberg([f["p_value"] for f in family], alpha=0.05)
    print(f"\nSECONDARY FAMILY ({len(family)} tests, BH q=0.05)")
    print(f"  {'endpoint':40s} {'n':>9s} {'delta':>7s} {'p':>7s} {'q':>7s}  sig")
    for item, q, rejected in zip(family, correction.adjusted, correction.rejected):
        item["q_value"], item["significant"] = float(q), bool(rejected)
        print(f"  {item['endpoint']:40s} {item['n_lof']:4d}/{item['n_gof']:<4d} "
              f"{item['cliffs_delta']:+7.3f} {item['p_value']:7.4f} "
              f"{q:7.4f}  {'YES' if rejected else 'no'}")

    out = {"preregistration": "docs/PREREGISTRATION_ROUND41.md",
           "gene_constraint": {"loeuf": constraint.loeuf, "pli": constraint.pli,
                               "mis_z": constraint.mis_z,
                               "oe_mis": constraint.oe_mis,
                               "lof_intolerant": constraint.lof_intolerant,
                               "missense_depleted": constraint.missense_depleted},
           "n_positions": len(combined), "excluded_positions": sorted(EXCLUDED_POSITIONS),
           "primary": primary, "decision": "reject" if reject else "fail_to_reject",
           "secondary": family,
           "any_secondary_significant": bool(any(f["significant"] for f in family))}
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    path = DERIVED_DIR / "validation_round41.json"
    path.write_text(json.dumps(out, indent=1, default=float))
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
