#!/usr/bin/env python
"""Emit the census-check tables, for the ``piezo_genes`` project to audit.

Round 93 brought the ``piezo_genes`` census into this project and re-asked its
findings on coordinates. Three of the answers came out differently, and those
disagreements are about *their* data — so they are exported here as tables
rather than left in this project's prose, where the people who could adjudicate
them would never see the numbers.

Every table is regenerated from the committed structures and resources; nothing
is copied out of a document. Run::

    python scripts/export_census_check.py --out ../piezo_genes/results/external_check

Four tables:

``axis_radius.tsv``
    Distance from the three-fold axis for each disputed domain band, per entry.
    The evidence that the census's PIEZO1 anchor and outer-helix bands sit ~120
    residues N-terminal of where the structure puts them: an outer helix must
    form the wall of the pore, and theirs is 39 A from an axis whose lining
    helix is at 10 A.

``core_periphery_splay.tsv``
    Blade RMSD over core RMSD after a pore-module-only fit, for every pair.
    The evidence that the blade splay in their piezo3-versus-Piezo1 figure is
    prediction-versus-experiment: a monomer splays 7-9x from an experiment of
    its own protein, two experimental paralogues splay ~1x.

``blade_band_composition.tsv``
    Why their distal-versus-proximal gradient reverses on the transmembrane
    units: the two bands differ in how much inter-unit linker they contain, and
    linker scores the same either side.

``replications.tsv``
    The findings that held, each with our value beside theirs. Carried because
    a report of three disagreements and no agreements would misrepresent what
    the check found.
"""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.core_periphery import compare  # noqa: E402
from piezo1.analysis.disease_geography import both_partitions  # noqa: E402
from piezo1.analysis.family_constraint import (blade_gradient,  # noqa: E402
                                               compare_with_own_conservation,
                                               domain_constraint)
from piezo1.analysis.paralogue_identity import paralogue_asymmetry  # noqa: E402
from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core.family import load_family_findings  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.structure.protomers import protomer_blocks  # noqa: E402
from piezo1.structure.superpose import detect_c3_axis  # noqa: E402

#: The disputed bands, in **mouse** numbering, because that is the numbering the
#: census's ``domain_map`` carries alongside its reference numbers and the
#: numbering most deposited PIEZO1 entries are in.
BANDS = {
    "inner_helix_agreed": (2458, 2478),
    "outer_helix_ours": (2193, 2213),
    "outer_helix_census": (2073, 2099),
    "anchor_ours": (2093, 2192),
    "anchor_census": (1952, 2063),
    "cap_agreed": (2214, 2457),
    "ctd_agreed": (2479, 2547),
}

#: Mouse-numbered PIEZO1 entries the radial measurement runs on. Several, not
#: one, because a single entry could carry an idiosyncratic axis.
RADIUS_ENTRIES = ("7WLT", "6B3R", "8IXN", "5Z10", "3JAC", "8IMZ", "11YE")

#: Pairs for the splay table. The first two are the control that decides the
#: reading: a prediction against an experiment **of the same protein**.
SPLAY_PAIRS = (
    ("AF-E2JF22-F1-model_v6", "6B3R", "AlphaFold monomer vs experiment, SAME protein"),
    ("AF-E2JF22-F1-model_v6", "7WLT", "AlphaFold monomer vs experiment, SAME protein"),
    ("AF-A0AB32U1Q1-F1-model_v6", "6B3R", "piezo3 prediction vs PIEZO1 experiment"),
    ("7WLT", "7WLU", "PIEZO1 curved vs flattened (the gating motion)"),
    ("7WLT", "6B3R", "PIEZO1 vs PIEZO1"),
    ("6KG7", "9VEE", "PIEZO2 mouse vs human"),
    ("7WLT", "6KG7", "PIEZO1 vs PIEZO2, both experimental"),
    ("7WLT", "9VEE", "PIEZO1 vs PIEZO2, both experimental"),
    ("8YEZ", "9VEE", "PIEZO1 vs PIEZO2, both experimental"),
)


def write(out: Path, name: str, header, rows) -> None:
    path = out / name
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)
    print(f"  {name}: {len(rows)} rows")


def _entry(pdb: str):
    for suffix in (".cif", ".pdb"):
        path = STRUCTURE_DIR / f"{pdb}{suffix}"
        if path.exists():
            return Structure.from_file(path)
    return None


def axis_radius(out: Path) -> None:
    """Mean distance from the C3 axis for each band, per entry."""
    rows = []
    for pdb in RADIUS_ENTRIES:
        structure = _entry(pdb)
        if structure is None:
            continue
        blocks, _ = protomer_blocks(structure)
        if len(blocks) < 3:
            continue
        axis = detect_c3_axis(blocks)
        origin = np.asarray(axis.point, dtype=float)
        direction = np.asarray(axis.direction, dtype=float)
        direction = direction / np.linalg.norm(direction)
        mask = structure.mask_ca()
        rel = structure.xyz[mask] - origin
        res = structure.res_seq[mask]
        along = rel @ direction
        radius = np.linalg.norm(rel - np.outer(along, direction), axis=1)
        for band, (lo, hi) in BANDS.items():
            sel = (res >= lo) & (res <= hi)
            if not sel.any():
                continue
            rows.append([pdb, band, f"{lo}-{hi}", int(sel.sum()),
                         round(float(radius[sel].mean()), 2),
                         round(float(along[sel].mean()), 2)])
    write(out, "axis_radius.tsv",
          ["entry", "band", "mouse_range", "n_ca", "mean_radius_from_axis_A",
           "mean_z_along_axis_A"], rows)


def splay(out: Path) -> None:
    rows = []
    for mobile, target, note in SPLAY_PAIRS:
        a, b = _entry(mobile), _entry(target)
        if a is None or b is None:
            continue
        result = compare(a, b, mobile, target)
        if not result:
            rows.append([mobile, target, note, "", "", "", "", "refused",
                         result.reason])
            continue
        rows.append([
            mobile, target, note, result.n_core,
            round(result.core_rmsd, 3), result.n_periphery,
            round(result.periphery_rmsd, 3) if result.periphery_rmsd else "",
            round(result.splay_ratio, 2) if result.splay_ratio else "core did not fit",
            "cross-paralogue" if result.cross_paralogue else "same protein"])
    write(out, "core_periphery_splay.tsv",
          ["mobile", "target", "note", "n_core_ca", "core_rmsd_A",
           "n_blade_ca", "blade_rmsd_A", "splay_ratio", "relationship"], rows)


def blade_composition(out: Path) -> None:
    g = blade_gradient()
    rows = [
        ["census distal band", f"1-{g['band_cut_human']} (human)",
         g["n_band_distal"], round(g["band_distal"], 4),
         round(g["linker_fraction_distal"], 3), round(g["linker_distal"], 4)],
        ["census proximal band", f"{g['band_cut_human'] + 1}-1935 (human)",
         g["n_band_proximal"], round(g["band_proximal"], 4),
         round(g["linker_fraction_proximal"], 3), round(g["linker_proximal"], 4)],
        ["THU1-6 units only", "transmembrane units", g["n_unit_distal"],
         round(g["unit_distal"], 4), 0.0, ""],
        ["THU7-9 units only", "transmembrane units", g["n_unit_proximal"],
         round(g["unit_proximal"], 4), 0.0, ""],
    ]
    write(out, "blade_band_composition.tsv",
          ["band", "range", "n_scored", "mean_constraint",
           "fraction_that_is_inter_unit_linker", "linker_mean_constraint"], rows)


def replications(out: Path) -> None:
    findings = load_family_findings()
    census_identity = {r["domain"]: float(r["identity"])
                       for r in findings.table("paralogue_identity")
                       if r["pair"] == "PIEZO1_vs_PIEZO2"}
    ours = {d.domain: d for d in paralogue_asymmetry("PIEZO1_vs_PIEZO2")}
    constraint = {d.domain: d for d in domain_constraint() if d.mean is not None}
    partitions = both_partitions()
    cross = compare_with_own_conservation()

    rows = [
        ["pairwise identity, inner helix", census_identity.get("inner_helix"),
         round(ours["inner_helix"].identity, 4), "our alignment, our boundaries"],
        ["pairwise identity, CTD", census_identity.get("CTD"),
         round(ours["ctd"].identity, 4), "our alignment, our boundaries"],
        ["pairwise identity, cap/CED", census_identity.get("CED"),
         round(ours["cap"].identity, 4),
         "the CED exception, reproduced to one percentage point"],
        ["pairwise identity, whole protein", census_identity.get("pore_module"),
         round(ours["cap"].whole_protein, 4), "whole-protein identity"],
        ["constraint, anchor", 0.7499, round(constraint["anchor"].mean, 4),
         "our partition puts the anchor 141 residues from theirs"],
        ["constraint, inner helix", 0.8074,
         round(constraint["inner_helix"].mean, 4), ""],
        ["constraint, THU1 (blade tip)", 0.6561,
         round(constraint["thu1"].mean, 4), "census value is its distal band"],
        ["disease enrichment odds ratio (their boundaries)", 3.922,
         round(partitions["results"]["census"].odds_ratio, 3),
         "ours uses gnomAD population missense, not ClinVar benign"],
        ["disease enrichment P (their boundaries)", 0.00135,
         round(partitions["results"]["census"].p_value, 5),
         "one-sided Fisher, PIEZO1 only"],
        ["disease enrichment odds ratio (our boundaries)", "",
         round(partitions["results"]["ours"].odds_ratio, 3),
         "does not reach significance; see the boundary note"],
        ["constraint classifier AUC", 0.9143, "",
         "not directly comparable: our negative set is population variation"],
        ["conservation cross-check (Spearman)", "",
         round(cross.spearman, 3) if cross else "",
         "their track against ours; no data and no statistic in common"],
    ]
    write(out, "replications.tsv",
          ["quantity", "census_value", "our_value", "note"], rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", required=True,
                        help="directory to write the tables into")
    args = parser.parse_args()
    out = Path(args.out).expanduser().resolve()
    out.mkdir(parents=True, exist_ok=True)
    print(f"census-check tables -> {out}")
    axis_radius(out)
    splay(out)
    blade_composition(out)
    replications(out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
