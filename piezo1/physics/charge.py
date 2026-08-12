"""A current and a count of ions are the same statement, in two units.

Small on purpose. These two lived in :mod:`piezo1.render.flux`, where the
animation needed them — and that was fine until :mod:`piezo1.physics.martini`
and :mod:`piezo1.analysis.liu2025_permeation` needed them too, at which point
``physics`` was importing from ``render`` and the project's one-way dependency
arrow was pointing backwards. ``tests/test_architecture.py`` caught it, which
is what it is for.

They belong here rather than in :mod:`piezo1.physics.permeation`, which is
already at the length limit, and rather than in ``_pnp_kernels``, which is
about turning equations into a matrix. ``render.flux`` re-exports both, so
every existing import keeps working.
"""

from __future__ import annotations

__all__ = ["ELEMENTARY_CHARGE", "ion_rate"]

#: Coulombs. Exact by the 2019 SI definition, like ``F_FARADAY`` in
#: :mod:`piezo1.physics._pnp_kernels`, which is why it is a module constant
#: rather than a registered parameter — it is a definition, not a measurement.
ELEMENTARY_CHARGE = 1.602176634e-19


def ion_rate(current_pA: float, valence: int = 1) -> float:
    """Ions per second carried by ``current_pA`` at the given valence.

    A divalent ion carries twice the charge, so the same current is half as
    many ions — which matters because PIEZO1 is calcium-permeable and an
    animation would otherwise overstate the particle count for Ca2+.
    """
    if valence <= 0:
        raise ValueError("valence must be positive")
    return abs(current_pA) * 1e-12 / (ELEMENTARY_CHARGE * valence)
