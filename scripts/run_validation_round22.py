#!/usr/bin/env python
"""Execute the Round 22 exploratory test exactly as pre-registered.

Runs only what `docs/PREREGISTRATION_ROUND22.md` specifies, in the order it
specifies, and writes the result whatever it shows. Nothing here chooses an
analysis after seeing a number.

Usage::

    python scripts/run_validation_round22.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.design import (benjamini_hochberg,  # noqa: E402
                                    leave_one_out, power_curve)
from piezo1.analysis.external import ProtVarClient  # noqa: E402
from piezo1.analysis.validation import (auroc, bootstrap_cliffs_delta,  # noqa: E402
                                        cliffs_delta, interpret_delta,
                                        permutation_test)
from piezo1.config import DERIVED_DIR, HUMAN_ACC  # noqa: E402
from piezo1.core.annotations import load_annotations  # noqa: E402

#: Sign applied so that every feature is expected *positive for LoF* under the
#: pre-registered mechanism, before any of them is looked at.
ORIENTATION = {"foldx_ddg": +1.0, "alphamissense": +1.0, "eve": +1.0,
               "esm1b": +1.0, "conservation": +1.0, "mechanical": -1.0}


def collect() -> tuple[list[dict], dict]:
    """Assemble the analysis set under the pre-registered inclusion criteria."""
    ann = load_annotations("human")
    directional = [v for v in ann.variants if v.classification in ("GoF", "LoF")]
    dropped = []
    included = []
    for v in directional:
        single = (v.residue and v.wt_aa and v.mut_aa
                  and len(v.wt_aa) == 1 and len(v.mut_aa) == 1
                  and v.mut_aa.isalpha())
        (included if single else dropped).append(v)

    client = ProtVarClient(offline=True)
    rows = []
    for v in included:
        scores = client.scores(HUMAN_ACC, v.residue, v.mut_aa)
        entry = client.foldx(HUMAN_ACC, v.residue).get(v.mut_aa) or {}
        rows.append({
            "label": v.label, "residue": int(v.residue),
            "wt": v.wt_aa, "mut": v.mut_aa,
            "classification": v.classification,
            "is_lof": v.classification == "LoF",
            "domain": v.domain or "",
            "foldx_ddg": entry.get("ddg"),
            "alphamissense": getattr(scores, "alphamissense", None),
            "eve": getattr(scores, "eve", None),
            "esm1b": getattr(scores, "esm1b", None),
            "conservation": getattr(scores, "conservation", None),
        })
    counts = {
        "curated": len(ann.variants),
        "directional": len(directional),
        "dropped_not_single_substitution": len(dropped),
        "dropped_labels": [v.label for v in dropped],
        "dropped_lof": sum(1 for v in dropped if v.classification == "LoF"),
        "included": len(rows),
        "included_gof": sum(1 for r in rows if not r["is_lof"]),
        "included_lof": sum(1 for r in rows if r["is_lof"]),
    }
    return rows, counts


def attach_mechanical(rows: list[dict]) -> int:
    """Join the Round 7 mechanical ΔΔG where it exists, for the secondary family."""
    path = DERIVED_DIR / "validation_round7.json"
    if not path.exists():
        return 0
    lookup = {e["label"]: e.get("ddg")
              for e in json.loads(path.read_text())["per_variant"]}
    n = 0
    for row in rows:
        row["mechanical"] = lookup.get(row["label"])
        n += row["mechanical"] is not None
    return n


def test_feature(rows: list[dict], key: str, alternative: str) -> dict:
    """One pre-registered comparison. LoF is the positive class throughout."""
    usable = [r for r in rows if r.get(key) is not None]
    lof = np.array([r[key] for r in usable if r["is_lof"]], dtype=float)
    gof = np.array([r[key] for r in usable if not r["is_lof"]], dtype=float)
    if len(lof) < 2 or len(gof) < 2:
        return {"feature": key, "error": "too few in one class",
                "n_lof": len(lof), "n_gof": len(gof)}

    # `permutation_test(a, b, alternative="less")` asks whether mean(a) < mean(b).
    # The hypothesis is LoF > GoF, so pass (gof, lof) with "less".
    perm = permutation_test(gof, lof, n_permutations=10000,
                            alternative=("less" if alternative == "lof_higher"
                                         else "two-sided"), seed=22)
    delta = cliffs_delta(lof, gof)
    effect = bootstrap_cliffs_delta(lof, gof, n_bootstrap=10000, seed=22)
    scores = np.array([r[key] for r in usable], dtype=float)
    labels = np.array([r["is_lof"] for r in usable], dtype=bool)
    return {
        "feature": key, "alternative": alternative,
        "n_lof": int(len(lof)), "n_gof": int(len(gof)),
        "mean_lof": float(lof.mean()), "mean_gof": float(gof.mean()),
        "median_lof": float(np.median(lof)), "median_gof": float(np.median(gof)),
        "p_value": float(perm.p_value),
        "cliffs_delta": float(delta),
        "ci_low": float(effect.ci_low), "ci_high": float(effect.ci_high),
        "interpretation": interpret_delta(delta),
        "excludes_zero": bool(effect.excludes_zero),
        "auroc": float(auroc(scores, labels)),
    }


def combined(rows: list[dict], keys: list[str]) -> dict:
    """Equal-weight standardised sum, scored out of sample. No fitting."""
    usable = [r for r in rows if all(r.get(k) is not None for k in keys)]
    if len(usable) < 8:
        return {"error": "too few complete rows", "n": len(usable)}
    matrix = np.array([[ORIENTATION[k] * float(r[k]) for k in keys]
                       for r in usable])
    labels = np.array([r["is_lof"] for r in usable], dtype=bool)
    result = leave_one_out(matrix, labels)
    return {"features": keys, "n": len(usable),
            "n_lof": int(labels.sum()), "n_gof": int((~labels).sum()),
            "auroc_in_sample": result.auroc_in,
            "auroc_leave_one_out": result.auroc_out,
            "optimism": result.optimism}


def main() -> int:
    rows, counts = collect()
    n_mech = attach_mechanical(rows)
    print(f"included {counts['included']} variants "
          f"({counts['included_gof']} GoF, {counts['included_lof']} LoF); "
          f"{counts['dropped_not_single_substitution']} dropped as not single "
          f"substitutions, {counts['dropped_lof']} of them LoF")
    print(f"mechanical ddG joined for {n_mech}")

    print("\n--- PRIMARY (pre-registered): FoldX ddG, LoF > GoF ---")
    primary = test_feature(rows, "foldx_ddg", "lof_higher")
    for k in ("n_lof", "n_gof", "mean_lof", "mean_gof", "p_value",
              "cliffs_delta", "ci_low", "ci_high", "interpretation", "auroc"):
        print(f"  {k:16s} {primary.get(k)}")

    print("\n--- SECONDARY family (BH corrected together) ---")
    secondary = []
    for key, alternative in (("alphamissense", "two_sided"),
                             ("eve", "two_sided"),
                             ("esm1b", "two_sided"),
                             ("conservation", "two_sided"),
                             ("mechanical", "gof_lower")):
        if key == "mechanical" and n_mech == 0:
            continue
        result = test_feature(rows, key, "lof_higher"
                              if alternative == "gof_lower" else alternative)
        result["alternative"] = alternative
        secondary.append(result)
        if "error" in result:
            print(f"  {key:16s} {result['error']}")
        else:
            print(f"  {key:16s} n={result['n_lof']}+{result['n_gof']:2d} "
                  f"p={result['p_value']:.3f} delta={result['cliffs_delta']:+.3f} "
                  f"[{result['ci_low']:+.3f},{result['ci_high']:+.3f}] "
                  f"AUROC={result['auroc']:.3f}")

    usable = [r for r in secondary if "error" not in r]
    family = benjamini_hochberg([r["p_value"] for r in usable],
                                [r["feature"] for r in usable], alpha=0.05)
    print("\n  BH-adjusted:")
    for row in family.table():
        print(f"    {row['name']:16s} p={row['p']:.3f}  q={row['q']:.3f}  "
              f"{'significant' if row['significant'] else ''}")

    print("\n--- COMBINED (equal weights, leave-one-out) ---")
    keys = ["foldx_ddg", "alphamissense", "eve", "esm1b", "conservation"]
    comb = combined(rows, keys)
    for k, v in comb.items():
        print(f"  {k}: {v}")

    print("\n--- POWER at the achieved design ---")
    power = power_curve(counts["included_gof"], counts["included_lof"],
                        deltas=[-0.2, -0.3, -0.43, -0.61],
                        n_simulations=3000, n_permutations=999, seed=22)
    for d, p in zip(power.deltas, power.power):
        print(f"  |delta| {abs(d):.2f} -> power {p:.3f}")

    out = {"counts": counts, "primary": primary, "secondary": secondary,
           "bh": family.table(), "combined": comb,
           "power": {"deltas": [float(d) for d in power.deltas],
                     "power": [float(p) for p in power.power],
                     "mde_80": power.detectable(0.8)},
           "preregistration": "docs/PREREGISTRATION_ROUND22.md",
           "design": "exploratory — effect sizes, not decisions"}
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    path = DERIVED_DIR / "validation_round22.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"\nwritten to {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
