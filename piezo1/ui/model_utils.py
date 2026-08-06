"""Structural helpers the window needs, kept out of it and testable.

Both of these ask the same question — which residues are resolved in *all three*
protomers — and both are load-bearing. A block set built from chains that do not
share a residue basis produces an elastic network of the wrong size, and every
result downstream is then a plausible wrong number rather than an error.
"""

from __future__ import annotations

import numpy as np

from ..core.structure import Structure

__all__ = ["well_resolved_chains", "modelled_residues", "protomer_blocks"]

#: A protomer needs at least this many C-alphas to count. Deposited entries
#: often carry a short peptide or a partial chain that would otherwise be
#: mistaken for a fourth subunit.
MIN_CA_PER_PROTOMER = 300


def well_resolved_chains(st: Structure) -> list[str]:
    return [ch for ch in st.chains
            if (st.mask_ca() & (st.chain == ch)).sum() > MIN_CA_PER_PROTOMER]


def modelled_residues(st: Structure) -> set[int]:
    """Residue numbers resolved in every well-resolved chain."""
    per = [set(st.res_seq[st.mask_ca() & (st.chain == ch)].tolist())
           for ch in well_resolved_chains(st)]
    return set.intersection(*per) if per else set()


def protomer_blocks(st: Structure) -> tuple[list[np.ndarray], np.ndarray]:
    """Equal-length C-alpha blocks per protomer, plus the residues they span.

    The residue array is returned rather than left implicit because anything
    mapping a per-site result back onto the model needs it, and rebuilding it
    separately is how the two fall out of step.
    """
    chains = [(st.xyz[st.mask_ca() & (st.chain == ch)],
               st.res_seq[st.mask_ca() & (st.chain == ch)])
              for ch in well_resolved_chains(st)]
    if len(chains) < 3:
        return [], np.array([], dtype=np.int64)

    common = set(chains[0][1].tolist())
    for _, seq in chains[1:]:
        common &= set(seq.tolist())
    residues = np.array(sorted(common))
    blocks = [xyz[np.searchsorted(seq, residues)].astype(np.float64)
              for xyz, seq in chains[:3]]
    return blocks, residues
