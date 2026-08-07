#!/usr/bin/env python
"""Execute the Round 36 test exactly as pre-registered.

Runs only what ``docs/PREREGISTRATION_ROUND36.md`` specifies, in the order it
specifies, and writes the result whatever it shows. Nothing here chooses an
analysis after seeing a number: the primary endpoint, the secondary family, the
inclusion criteria and the decision rule were all committed first, in their own
commit.

Usage::

    python scripts/run_validation_round36.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.design import benjamini_hochberg  # noqa: E402
from piezo1.analysis.external import ProtVarClient  # noqa: E402
from piezo1.analysis.validation import (auroc, bootstrap_cliffs_delta,  # noqa: E402
                                        cliffs_delta, permutation_test)
from piezo1.analysis.variant_sets import build_analysis_set  # noqa: E402
from piezo1.config import DERIVED_DIR, HUMAN_ACC  # noqa: E402
from piezo1.core.sequence import human_sequence, human_to_mouse, mouse_to_human  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.io.registry import load_registry  # noqa: E402
from piezo1.structure.protomers import protomer_blocks  # noqa: E402
from piezo1.structure.superpose import kabsch, match_protomers  # noqa: E402

#: Pre-registered §4.4 — the one variant whose curated and inferred directions
#: disagree. Excluded here rather than after seeing its score.
EXCLUDED = {"V598M"}

#: Sign making every feature *positive for LoF* under the pre-registered
#: mechanism, fixed before any of them was computed. The mechanical predictor is
#: hypothesised NEGATIVE for GoF (softening), hence −1.
ORIENTATION = {"mechanical": -1.0, "mechanical_volume_only": -1.0,
               "foldx_ddg": +1.0, "alphamissense": +1.0, "eve": +1.0,
               "esm1b": +1.0, "conservation": +1.0}


def _structure(pdb: str) -> Structure:
    record = load_registry().get(pdb)
    if record is None or not record.available:
        raise SystemExit(f"{pdb} not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(record.path)


def build_model(substitution_aware: bool):
    """The frozen predictor of §3, on the curved→flat gating vector."""
    from piezo1.analysis.variant_impact import VariantImpactModel

    curved, flat = _structure("7WLT"), _structure("7WLU")
    _cb, cr = protomer_blocks(curved)
    _fb, fr = protomer_blocks(flat)
    common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))

    def resample(st):
        out = []
        for chain in st.chains:
            mask = st.mask_ca() & (st.chain == chain)
            if mask.sum() < 300:
                continue
            index = {int(r): i for i, r in enumerate(st.res_seq[mask])}
            xyz = st.xyz[mask]
            if all(r in index for r in common):
                out.append(np.array([xyz[index[r]] for r in common], float))
        return out[:3]

    cb, fb = resample(curved), resample(flat)
    fb = [fb[i] for i in match_protomers(cb, fb).order]
    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    displacement = (((np.vstack(fb) - centroid) @ rotation.T + translation)
                    - np.vstack(cb))

    reference = human_sequence()
    sequence = {}
    for residue in common:
        human = mouse_to_human(int(residue))
        if human and 1 <= human <= len(reference):
            sequence[int(residue)] = reference[human - 1]

    return VariantImpactModel(coords=np.vstack(cb), residues=np.tile(common, 3),
                              gating_vector=displacement, sequence=sequence,
                              substitution_aware=substitution_aware)


def collect() -> tuple[list[dict], dict]:
    """Assemble the analysis set under the pre-registered criteria (§4, §5)."""
    combined = build_analysis_set(
        levels=("measured", "disease_mechanism")).missense()
    measured = {e.label for e in build_analysis_set(
        levels=("measured",)).missense().entries}

    aware = build_model(substitution_aware=True)
    plain = build_model(substitution_aware=False)
    client = ProtVarClient(offline=True)

    rows, dropped = [], {"excluded_by_name": 0, "no_mouse_equivalent": 0,
                         "not_modelled": 0}
    for entry in combined.entries:
        if entry.label in EXCLUDED:
            dropped["excluded_by_name"] += 1
            continue
        mouse = human_to_mouse(entry.residue)
        if mouse is None:
            dropped["no_mouse_equivalent"] += 1
            continue
        prediction = aware.predict(mouse, entry.wt_aa, entry.mut_aa)
        if not prediction.modelled or not np.isfinite(prediction.gating_cost_change):
            dropped["not_modelled"] += 1
            continue
        control = plain.predict(mouse, entry.wt_aa, entry.mut_aa)

        # A mutant MUST be passed: a position-only query returns nineteen
        # unlabelled entries per predictor, so only conservation is meaningful
        # without it.
        external = client.scores(HUMAN_ACC, entry.residue, entry.mut_aa)
        rows.append({
            "label": entry.label,
            "classification": entry.classification,
            "evidence": entry.evidence,
            "measured": entry.label in measured,
            "mechanical": float(prediction.gating_cost_change),
            "mechanical_volume_only": float(control.gating_cost_change),
            "foldx_ddg": getattr(external, "foldx_ddg", None),
            "alphamissense": getattr(external, "alphamissense", None),
            "eve": getattr(external, "eve", None),
            "esm1b": getattr(external, "esm1b", None),
            "conservation": getattr(external, "conservation", None),
        })
    return rows, dropped


def test_feature(rows: list[dict], feature: str, subset=None) -> dict | None:
    """One endpoint: Cliff's delta, permutation p, bootstrap CI, AUROC."""
    use = [r for r in rows if subset is None or subset(r)]
    gof = np.array([r[feature] for r in use
                    if r["classification"] == "GoF" and r[feature] is not None],
                   dtype=float)
    lof = np.array([r[feature] for r in use
                    if r["classification"] == "LoF" and r[feature] is not None],
                   dtype=float)
    if len(gof) < 3 or len(lof) < 3:
        return None

    sign = ORIENTATION[feature]
    a, b = sign * gof, sign * lof            # oriented so H1 predicts a < b
    permutation = permutation_test(a, b, alternative="less")
    effect = bootstrap_cliffs_delta(a, b)
    # auroc() takes scores and a boolean mask, not two groups: passing two
    # arrays cast the second to all-True, left no negatives, and returned nan.
    scores = np.concatenate([a, b])
    positive = np.concatenate([np.zeros(len(a), bool), np.ones(len(b), bool)])
    return {
        "feature": feature, "n_gof": int(len(gof)), "n_lof": int(len(lof)),
        "cliffs_delta": float(cliffs_delta(a, b)),
        "ci_low": float(effect.ci_low), "ci_high": float(effect.ci_high),
        "p_value": float(permutation.p_value),
        "auroc": float(auroc(scores, positive)),
        "median_gof": float(np.median(gof)), "median_lof": float(np.median(lof)),
    }


def main() -> int:
    rows, dropped = collect()
    print(f"analysis set: {len(rows)} variants "
          f"({sum(r['classification'] == 'GoF' for r in rows)} GoF, "
          f"{sum(r['classification'] == 'LoF' for r in rows)} LoF); "
          f"dropped {dropped}")

    # ---- primary endpoint (§7), tested uncorrected -------------------------
    primary = test_feature(rows, "mechanical")
    if primary is None:
        raise SystemExit("primary endpoint has too few variants to test")
    print("\nPRIMARY — substitution-aware mechanical ddG, combined set")
    print(f"  n = {primary['n_gof']} GoF vs {primary['n_lof']} LoF")
    print(f"  Cliff's delta {primary['cliffs_delta']:+.3f} "
          f"[{primary['ci_low']:+.3f}, {primary['ci_high']:+.3f}]")
    print(f"  p = {primary['p_value']:.4f}   AUROC {primary['auroc']:.3f}")

    # ---- decision rule (§8), applied exactly as written --------------------
    reject = (primary["p_value"] < 0.05 and primary["cliffs_delta"] < 0
              and primary["ci_high"] < 0)
    print(f"\n  DECISION: {'REJECT H0' if reject else 'FAIL TO REJECT H0'}")

    # ---- secondary family (§7), all corrected together ---------------------
    family = []
    measured_only = test_feature(rows, "mechanical", subset=lambda r: r["measured"])
    if measured_only is not None:
        measured_only["feature"] = "mechanical (measured labels only)"
        family.append(measured_only)
    untestable = []
    for feature in ("mechanical_volume_only", "foldx_ddg", "alphamissense",
                    "eve", "esm1b", "conservation"):
        result = test_feature(rows, feature)
        if result is not None:
            family.append(result)
        else:
            # Pre-registered but not runnable on the available data. Recorded
            # rather than silently dropped: §7 forbids removing a test from the
            # family after the fact, and "could not be run" is a different
            # statement from "was not significant".
            have = sum(1 for r in rows if r[feature] is not None)
            untestable.append({"feature": feature, "values_available": have,
                               "reason": "too few variants carry this score"})

    correction = benjamini_hochberg([f["p_value"] for f in family], alpha=0.05)
    print(f"\nSECONDARY FAMILY ({len(family)} tests, Benjamini-Hochberg q=0.05)")
    print(f"  {'endpoint':38s} {'n':>9s} {'delta':>7s} {'p':>7s} {'q':>7s}  sig")
    for entry, q, rejected in zip(family, correction.adjusted, correction.rejected):
        entry["q_value"], entry["significant"] = float(q), bool(rejected)
        print(f"  {entry['feature']:38s} {entry['n_gof']:4d}/{entry['n_lof']:<4d} "
              f"{entry['cliffs_delta']:+7.3f} {entry['p_value']:7.4f} "
              f"{q:7.4f}  {'YES' if rejected else 'no'}")

    out = {
        "preregistration": "docs/PREREGISTRATION_ROUND36.md",
        "n_variants": len(rows), "dropped": dropped,
        "primary": primary, "decision": "reject" if reject else "fail_to_reject",
        "secondary": family,
        "any_secondary_significant": bool(any(f["significant"] for f in family)),
        "untestable_endpoints": untestable,
        "variants": rows,
    }
    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    path = DERIVED_DIR / "validation_round36.json"
    path.write_text(json.dumps(out, indent=1, default=float))
    if untestable:
        print("\nPRE-REGISTERED BUT NOT RUNNABLE (reported, not dropped):")
        for entry in untestable:
            print(f"  {entry['feature']:38s} {entry['values_available']}/"
                  f"{len(rows)} variants carry a value")
    print(f"\nwrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
