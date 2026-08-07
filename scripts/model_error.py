#!/usr/bin/env python
"""Estimate systematic (model-form) error and compare it with sampling error.

`uncertainty.py` reports three kinds of spread and says on each that model error
is not among them. This estimates it where a second defensible model exists, and
answers the question that matters: which of the two dominates?

Usage::

    python scripts/model_error.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.claims import _tm_surface  # noqa: E402
from piezo1.analysis.model_error import (compare_with_sampling,  # noqa: E402
                                         dome_model_error, fit_spheroid,
                                         pore_convention_error,
                                         spring_model_error)
from piezo1.analysis.uncertainty import bootstrap  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.io.registry import load_registry  # noqa: E402
from piezo1.structure.geometry import fit_sphere  # noqa: E402
from piezo1.structure.protomers import protomer_blocks  # noqa: E402
from piezo1.structure.superpose import (detect_c3_axis, kabsch,  # noqa: E402
                                        match_protomers)


def _load(pdb: str) -> Structure:
    record = load_registry().get(pdb)
    if record is None or not record.available:
        raise SystemExit(f"{pdb} not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(record.path)


def main() -> int:
    curved, flat = _load("7WLT"), _load("7WLU")
    blocks, curved_res = protomer_blocks(curved)
    axis = detect_c3_axis(blocks)
    surface = _tm_surface(curved, "mouse")

    print("=" * 72)
    print("DOME CURVATURE — sphere vs oblate spheroid")
    sphere = fit_sphere(surface)
    spheroid = fit_spheroid(surface, axis.direction)
    residual = np.abs(np.linalg.norm(surface - sphere.center, axis=1)
                      - sphere.radius)
    print(f"  sphere   R = {sphere.radius / 10:.3f} nm, "
          f"geometric rmse {np.sqrt((residual ** 2).mean()):.3f} A")
    print(f"  spheroid a = {spheroid.equatorial / 10:.3f} nm, "
          f"c = {spheroid.polar / 10:.3f} nm, "
          f"flattening {spheroid.flattening:+.3f}, "
          f"apex {spheroid.apex_curvature / 10:.3f} nm, "
          f"geometric rmse {spheroid.rmse:.3f} A")
    dome = dome_model_error(surface, axis, sphere.radius / 10)
    print(f"  {dome.summary()}")

    # Sampling interval for the same quantity: resample the surface points.
    def radius_of(index):
        return fit_sphere(surface[index]).radius / 10.0

    sampling = bootstrap(radius_of, surface, n_resamples=400,
                         what="dome radius of curvature")
    print(f"  sampling  {sampling.summary(unit='nm')}")
    verdict = compare_with_sampling(dome, sampling)
    print(f"  --> {verdict['verdict']}")

    print()
    print("=" * 72)
    print("GATING OVERLAP — three published spring models")
    _fb, flat_res = protomer_blocks(flat)
    common = np.array(sorted(set(curved_res.tolist()) & set(flat_res.tolist())))
    curved_blocks = [b[np.isin(curved_res, common)] for b in blocks]
    flat_blocks = [b[np.isin(flat_res, common)] for b in _fb]
    flat_blocks = [flat_blocks[i]
                   for i in match_protomers(curved_blocks, flat_blocks).order]
    rotation, translation, centroid = kabsch(np.vstack(flat_blocks),
                                             np.vstack(curved_blocks))
    displacement = (((np.vstack(flat_blocks) - centroid) @ rotation.T
                     + translation) - np.vstack(curved_blocks))
    springs = spring_model_error(curved_blocks, displacement)
    print(f"  {springs.summary()}")

    print()
    print("=" * 72)
    print("PORE BOTTLENECK — Apollonius vs uniform probe")
    for radius in (1.4, 1.7, 2.0):
        pore = pore_convention_error(curved, axis, uniform_radius=radius)
        values = list(pore.values.values())
        print(f"  uniform {radius} A: Apollonius {values[0]:.3f} vs "
              f"uniform {values[1]:.3f} A  (spread {pore.spread:.3f})")
    print(f"  {pore.note}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
