#!/usr/bin/env python
"""Run the blind variant validation exactly as pre-registered.

Executes the protocol in ``docs/PREREGISTRATION.md`` and writes
``docs/VALIDATION.md``. Nothing here may be tuned in response to the result:
the predictor, the inclusion criteria, the statistic and the decision rule were
all fixed before this script was first run.

    python scripts/run_validation.py
"""

from __future__ import annotations

import json
import sys
from dataclasses import asdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.ensemble import build_ensemble  # noqa: E402
from piezo1.analysis.validation import (auroc, bootstrap_cliffs_delta,  # noqa: E402
                                        permutation_test)
from piezo1.analysis.variant_impact import VariantImpactModel  # noqa: E402
from piezo1.core.annotations import load_annotations  # noqa: E402

N_PERM = 10000
N_BOOT = 10000
SEED = 0


def build_predictions():
    ens = build_ensemble(species="mouse", min_common=900)
    pca = ens.pca()
    model = VariantImpactModel(coords=ens.members[0].coords,
                               residues=np.tile(ens.residues, 3),
                               gating_vector=pca.components[0])
    ann = load_annotations("human")
    preds = model.predict_all(ann.variants, annotations=ann)
    by_id = {v.id: v for v in ann.variants}
    rows = []
    for p, v in zip(preds, [v for v in ann.variants if v.residue is not None]):
        rows.append({"pred": p, "variant": v})
    return ens, pca, model, ann, rows, by_id


def apply_inclusion(rows):
    """Pre-registered inclusion criteria, with the counts recorded."""
    counts = {"total": len(rows)}
    kept, dropped = [], {}

    def drop(row, reason):
        dropped.setdefault(reason, []).append(row["variant"].label)

    for row in rows:
        p, v = row["pred"], row["variant"]
        if v.classification not in ("GoF", "LoF"):
            drop(row, f"class {v.classification}")
            continue
        if not p.modelled:
            drop(row, "residue not in the resolved network")
            continue
        if not (p.wt_aa and p.mut_aa and len(p.wt_aa) == 1 and len(p.mut_aa) == 1):
            drop(row, "not a single-residue substitution")
            continue
        kept.append(row)
    counts["included"] = len(kept)
    counts["dropped"] = {k: len(v) for k, v in dropped.items()}
    counts["dropped_labels"] = dropped
    return kept, counts


def primary_test(kept, field: str = "ddg_gating") -> dict:
    gof = np.array([r["pred"].__dict__[field] for r in kept
                    if r["variant"].classification == "GoF"])
    lof = np.array([r["pred"].__dict__[field] for r in kept
                    if r["variant"].classification == "LoF"])
    perm = permutation_test(gof, lof, N_PERM, alternative="less", seed=SEED)
    effect = bootstrap_cliffs_delta(gof, lof, N_BOOT, seed=SEED)
    scores = np.concatenate([gof, lof])
    labels = np.concatenate([np.ones(len(gof), bool), np.zeros(len(lof), bool)])
    # AUROC for ranking GoF *above* LoF on -ddG (softer = more GoF-like).
    return {"field": field, "n_gof": len(gof), "n_lof": len(lof),
            "mean_gof": float(gof.mean()), "mean_lof": float(lof.mean()),
            "median_gof": float(np.median(gof)), "median_lof": float(np.median(lof)),
            "permutation": asdict(perm), "effect": asdict(effect),
            "auroc": auroc(-scores, labels)}


def decision(primary: dict) -> str:
    p = primary["permutation"]["p_value"]
    lo, hi = primary["effect"]["ci_low"], primary["effect"]["ci_high"]
    excludes_zero = (lo > 0) or (hi < 0)
    if p < 0.05 and excludes_zero:
        return "H1 supported"
    if p < 0.05 and not excludes_zero:
        return "H0 not rejected (p < 0.05 but the effect interval spans zero)"
    return "H0 not rejected"


def main() -> int:
    print("Running the pre-registered blind validation.\n")
    ens, pca, model, ann, rows, _ = build_predictions()
    kept, counts = apply_inclusion(rows)

    print(f"gating coordinate: ensemble PC1, {pca.variance_explained[0]:.1%} variance")
    print(f"variants considered: {counts['total']}")
    for reason, n in sorted(counts["dropped"].items(), key=lambda kv: -kv[1]):
        print(f"   dropped {n:3d}  {reason}")
    print(f"   INCLUDED {counts['included']}")

    if counts["included"] < 6:
        print("\nToo few variants survive inclusion for the test to mean anything.")
        return 1

    primary = primary_test(kept, "ddg_gating")
    secondary_norm = primary_test(kept, "ddg_normalised")

    print(f"\nGoF n={primary['n_gof']}  mean ddG {primary['mean_gof']:+.3e}")
    print(f"LoF n={primary['n_lof']}  mean ddG {primary['mean_lof']:+.3e}")
    print(f"\nPRIMARY  difference {primary['permutation']['observed']:+.3e}")
    print(f"   one-sided permutation p = {primary['permutation']['p_value']:.4f}")
    print(f"   Cliff's delta {primary['effect']['delta']:+.3f} "
          f"[{primary['effect']['ci_low']:+.3f}, {primary['effect']['ci_high']:+.3f}] "
          f"({primary['effect']['interpretation']})")
    print(f"   AUROC {primary['auroc']:.3f}")
    print(f"\n   DECISION: {decision(primary)}")

    print(f"\nSECONDARY (ddg_normalised): p = "
          f"{secondary_norm['permutation']['p_value']:.4f}, "
          f"delta {secondary_norm['effect']['delta']:+.3f}, "
          f"AUROC {secondary_norm['auroc']:.3f}")

    # Domain stratification, pre-specified to test the provenance caveat.
    from collections import Counter
    dom_gof = Counter(r["pred"].domain for r in kept
                      if r["variant"].classification == "GoF")
    dom_lof = Counter(r["pred"].domain for r in kept
                      if r["variant"].classification == "LoF")
    print(f"\nDomain distribution  GoF: {dict(dom_gof)}")
    print(f"                     LoF: {dict(dom_lof)}")

    # R2456H, recorded as a single case that cannot support a conclusion.
    r2456 = [r for r in kept if r["variant"].residue == 2456]
    for r in r2456:
        p = r["pred"]
        print(f"\nR2456 case: {r['variant'].label} ({r['variant'].classification}) "
              f"ddG {p.ddg_gating:+.3e} -> predicted {p.direction}")

    out = {"counts": counts, "primary": primary, "secondary_normalised": secondary_norm,
           "decision": decision(primary),
           "domain_gof": dict(dom_gof), "domain_lof": dict(dom_lof),
           "gating_variance": float(pca.variance_explained[0]),
           "n_structures": len(ens), "n_residues": int(len(ens.residues)),
           "r2456": [{"label": r["variant"].label,
                      "classification": r["variant"].classification,
                      "ddg": r["pred"].ddg_gating,
                      "direction": r["pred"].direction} for r in r2456],
           "per_variant": [{"label": r["variant"].label,
                            "classification": r["variant"].classification,
                            "domain": r["pred"].domain,
                            "ddg": r["pred"].ddg_gating,
                            "ddg_normalised": r["pred"].ddg_normalised}
                           for r in kept]}
    dest = Path("data/derived/validation_round7.json")
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
