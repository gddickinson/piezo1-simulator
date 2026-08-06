"""Backwards-compatible re-export of the protomer helpers.

These moved to :mod:`piezo1.structure.protomers` when ``structure/frame.py``
needed them: importing them from ``ui`` would have pointed the dependency arrow
backwards. Kept as a shim so existing imports do not break.
"""

from __future__ import annotations

from ..structure.protomers import (MIN_CA_PER_PROTOMER, modelled_residues,
                                   protomer_blocks, well_resolved_chains)

__all__ = ["well_resolved_chains", "modelled_residues", "protomer_blocks",
           "MIN_CA_PER_PROTOMER"]
