"""Counting atoms buried in each other — the honest symptom of a bad model.

Split from :mod:`piezo1.structure.assembly` at the 500-line limit and along a
real seam: this is a geometric utility about a pair of atoms and knows nothing
about trimers, templates or PIEZO.

It exists because :mod:`~piezo1.structure.assembly` models no inter-protomer
contact at all. If the monomer it places is shaped differently from the template
protomer it is placed onto, the copies interpenetrate, and that count is the
measurement that says so — reported rather than relaxed away.

**Calibrated, because the number means nothing on its own.** An assembly scores
in the thousands, which is only interpretable against what a real trimer scores
on the same counter: 6B3R gives **8**, 7WLT **3** and 9ZIS **6**. Pinned in
``tests/test_assembly.py``, without which "10,162 clashes" would be equally
consistent with a counter that overcounts everything.
"""

from __future__ import annotations

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["count_clashes"]


def count_clashes(structure: Structure) -> int:
    """Heavy-atom pairs closer than the cutoff between different protomers.

    A grid rather than a full pairwise matrix: three copies of a 20,000-atom
    monomer is 3.6 x 10^9 pairs, which is neither necessary nor affordable.
    """
    cutoff = float(_P.value("assembly.clash_distance"))
    keep = structure.element != "H"
    xyz = structure.xyz[keep]
    chain = structure.chain[keep]
    if xyz.size == 0:
        return 0

    cell = np.floor(xyz / cutoff).astype(np.int64)
    buckets: dict = {}
    for i, key in enumerate(map(tuple, cell)):
        buckets.setdefault(key, []).append(i)

    total = 0
    offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1)
               for dz in (-1, 0, 1)]
    for key, members in buckets.items():
        near = []
        for offset in offsets:
            near.extend(buckets.get((key[0] + offset[0], key[1] + offset[1],
                                     key[2] + offset[2]), ()))
        if not near:
            continue
        block = np.asarray(near)
        for i in members:
            other = block[block > i]
            if other.size == 0:
                continue
            different = chain[other] != chain[i]
            if not different.any():
                continue
            candidate = other[different]
            distance = np.linalg.norm(xyz[candidate] - xyz[i], axis=1)
            total += int((distance < cutoff).sum())
    return total
