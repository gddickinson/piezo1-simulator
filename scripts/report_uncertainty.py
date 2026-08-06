#!/usr/bin/env python
"""Attach an interval to each headline number, and say which kind it is.

Rounds 18–28 repeatedly found a recorded number stated with more confidence than
it had earned. This reports the spread on each, distinguishing a bootstrap
confidence interval (resampling data) from a sensitivity range (varying a method
choice) — a cutoff has no sampling distribution, and calling its spread a
confidence interval would be a second kind of overconfidence.

Usage::

    python scripts/report_uncertainty.py
    python scripts/report_uncertainty.py --resamples 100      # quicker
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "tests"))

from piezo1.analysis.uncertainty import (bootstrap,  # noqa: E402
                                         parameter_range, sensitivity)
from piezo1.config import DERIVED_DIR, STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.structure.protomers import protomer_blocks  # noqa: E402


def _shared_blocks(structure, residues):
    out = []
    for chain in structure.chains:
        mask = structure.mask_ca() & (structure.chain == chain)
        if mask.sum() < 300:
            continue
        index = {int(r): i for i, r in enumerate(structure.res_seq[mask])}
        xyz = structure.xyz[mask]
        if all(r in index for r in residues):
            out.append(np.array([xyz[index[r]] for r in residues], float))
    return out[:3]


def dome_radius(resamples: int):
    from piezo1.structure.geometry import fit_sphere
    from test_geometry import _tm_surface

    structure = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
    surface = _tm_surface(structure, "mouse")
    return bootstrap(
        lambda idx: fit_sphere(surface[idx], iterations=4, trim=0.15).radius / 10.0,
        surface, n_resamples=resamples, seed=0,
        what="dome radius of curvature (7WLT), nm",
        note="resamples the transmembrane surface points; says nothing about "
             "whether a sphere is the right shape")


def gating_overlap():
    from piezo1.physics.anm import ANM
    from piezo1.structure.superpose import kabsch, match_protomers

    curved = Structure.from_file(STRUCTURE_DIR / "7WLT.cif")
    flat = Structure.from_file(STRUCTURE_DIR / "7WLU.cif")
    _c, cr = protomer_blocks(curved)
    _f, fr = protomer_blocks(flat)
    common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))
    cb, fb = _shared_blocks(curved, common), _shared_blocks(flat, common)
    fb = [fb[i] for i in match_protomers(cb, fb).order]
    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    displacement = ((((np.vstack(fb) - centroid) @ rotation.T + translation)
                     - np.vstack(cb))).ravel()

    def overlap_at(cutoff):
        anm = ANM.from_trimer(cb, cutoff=cutoff, spring="inverse_square").build()
        modes = anm.calc_modes(n_modes=30)
        anm.label_symmetry(modes)
        values = np.abs(np.asarray(modes.overlap(displacement), float))
        symmetric = np.array([s == "A" for s in modes.symmetry])
        return float(values[symmetric].max()) if symmetric.any() else 0.0

    return sensitivity(overlap_at, [10.0, 12.0, 13.5, 15.0, 16.5, 18.0, 20.0],
                       reference=15.0, knob="anm.cutoff",
                       what="lowest A-mode overlap with the gating transition")


def ensemble_pc1(resamples: int):
    from piezo1.analysis.ensemble import build_ensemble

    ensemble = build_ensemble(species="mouse", min_common=900)
    coords = np.array([m.coords for m in ensemble.members])

    def variance_explained(index):
        matrix = coords[index].reshape(len(index), -1)
        matrix = matrix - matrix.mean(axis=0)
        spectrum = np.linalg.svd(matrix, full_matrices=False,
                                 compute_uv=False) ** 2
        return float(spectrum[0] / spectrum.sum())

    return bootstrap(variance_explained, coords, n_resamples=resamples, seed=1,
                     what="ensemble PC1 variance explained",
                     note="resamples structures; ten is a small sample and the "
                          "interval says so")


def footprint_energy():
    from piezo1.physics.elastica import solve_elastica
    from piezo1.physics.membrane import MembraneParameters

    def energy(_kappa):
        return solve_elastica(8.691, 1.992, MembraneParameters()).energy

    return parameter_range(energy, "membrane.kappa", [20.0, 22.5, 25.0],
                           what="nonlinear footprint energy at the 7WLT "
                                "geometry, k_BT")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--resamples", type=int, default=400)
    args = parser.parse_args()

    results = [
        ("dome.radius_7wlt", dome_radius(args.resamples), "nm", 3),
        ("anm.gating_overlap", gating_overlap(), "", 3),
        ("ensemble.pc1_variance", ensemble_pc1(args.resamples), "", 4),
        ("elastica.energy", footprint_energy(), "k_BT", 2),
    ]

    print("Headline numbers with their spread\n")
    out = {}
    for key, spread, unit, digits in results:
        print(f"{key}")
        print(f"   {spread.summary(unit, digits)}")
        if getattr(spread, "note", ""):
            print(f"   caveat: {spread.note}")
        print()
        out[key] = {"kind": spread.kind, "estimate": spread.estimate,
                    "low": spread.low, "high": spread.high, "unit": unit,
                    "what": spread.what}

    print("None of these captures MODEL error — whether a sphere is the right")
    print("shape for the dome, or springs the right physics. A bootstrap tells")
    print("you how well a fit is determined, not whether it was the right fit.")

    DERIVED_DIR.mkdir(parents=True, exist_ok=True)
    path = DERIVED_DIR / "uncertainty.json"
    path.write_text(json.dumps(out, indent=1))
    print(f"\nwritten to {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
