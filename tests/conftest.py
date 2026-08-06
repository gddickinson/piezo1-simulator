"""Shared fixtures.

Tests that need real structures are skipped rather than failed when the data
has not been fetched, so the suite still runs on a fresh clone.
"""

import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402


def _require(pdb: str) -> Path:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return path


@pytest.fixture(scope="session")
def human_structure():
    return Structure.from_file(_require("8YEZ"))


@pytest.fixture(scope="session")
def curved_structure():
    return Structure.from_file(_require("7WLT"))


@pytest.fixture(scope="session")
def flat_structure():
    return Structure.from_file(_require("7WLU"))


def protomer_blocks(st, n=3):
    """Equal-length, identically ordered C-alpha blocks per protomer."""
    import numpy as np
    chains = []
    for ch in st.chains:
        m = st.mask_ca() & (st.chain == ch)
        if m.sum() > 300:
            chains.append((st.xyz[m], st.res_seq[m]))
    if len(chains) < n:
        return None, None
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:n]:
        common &= set(seq.tolist())
    arr = np.array(sorted(common))
    return [xyz[np.searchsorted(seq, arr)].astype(float)
            for xyz, seq in chains[:n]], arr
