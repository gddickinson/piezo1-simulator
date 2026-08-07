#!/usr/bin/env python
"""Run the four method-level cross-checks and report where the routes land.

Companion to ``scripts/crosscheck_chain.py``, which does the same for the three
headline physics results. Each check re-derives a quantity by a route sharing no
machinery with the pipeline's.

Usage::

    python scripts/crosscheck_methods.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.analysis.conservation import (OrthologSet,  # noqa: E402
                                          conservation_profile, load_orthologs)
from piezo1.analysis.crosscheck_methods import (check_conservation,  # noqa: E402
                                                check_pca, check_pore_radius,
                                                check_sasa,
                                                conservation_by_kmer_anchoring)
from piezo1.analysis.ensemble import build_ensemble  # noqa: E402
from piezo1.analysis.measure import sasa  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.io.registry import load_registry  # noqa: E402
from piezo1.structure.frame import apply_frame, canonical_transform  # noqa: E402
from piezo1.structure.pore import pore_profile  # noqa: E402
from piezo1.structure.protomers import protomer_blocks  # noqa: E402
from piezo1.structure.superpose import detect_c3_axis  # noqa: E402

#: Ortholog length window. A match-maximising anchor cannot bridge a large
#: length difference, so the comparison is restricted to orthologs where both
#: routes are being asked the same question.
LENGTH_WINDOW = 60
N_ORTHOLOGS = 8


def _load(pdb: str) -> Structure:
    record = load_registry().get(pdb)
    if record is None or not record.available:
        raise SystemExit(f"{pdb} not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(record.path)


def main() -> int:
    checks = []

    print("pore radius — pattern search vs brute force")
    structure = _load("8YEZ")
    framed = apply_frame(structure, canonical_transform(structure))
    blocks, _ = protomer_blocks(framed)
    axis = detect_c3_axis(blocks)
    profile = pore_profile(framed, axis, step=1.0)
    checks.append(check_pore_radius(framed, axis, profile))

    print("SASA — Shrake-Rupley vs Monte-Carlo")
    cap = _load("4RAX")                       # the cap domain: small enough to be exact
    mask = cap.mask_protein() & ~cap.hetero
    checks.append(check_sasa(cap, sasa(cap, mask=mask).total, mask=mask))

    print("conservation — Needleman-Wunsch vs k-mer seeds")
    orthologs = load_orthologs()
    reference = orthologs.members[0]
    close = [m for m in orthologs.members[1:]
             if abs(m.length - len(reference.sequence)) <= LENGTH_WINDOW][:N_ORTHOLOGS]
    subset = OrthologSet(members=[reference] + close, meta=dict(orthologs.meta))
    primary = conservation_profile(subset, use_cache=False)
    alternative = conservation_by_kmer_anchoring(
        reference.sequence, [m.sequence for m in close], window=120)
    checks.append(check_conservation(primary.identity, alternative))

    print("PCA — SVD vs power iteration")
    ensemble = build_ensemble(species="mouse")
    pca = ensemble.pca()
    matrix = np.array([m.coords.ravel() for m in ensemble.members])
    checks.append(check_pca(pca.eigenvalues[0], pca.components[0], matrix))

    print()
    for check in checks:
        print(check.summary())
        print(f"         {check.primary_route}")
        print(f"      vs {check.alternative_route}")
        if check.note:
            print(f"         {check.note}")
    disagreeing = [c for c in checks if not c.agrees]
    print(f"\n{len(checks) - len(disagreeing)}/{len(checks)} agree within tolerance")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
