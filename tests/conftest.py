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


@pytest.fixture(scope="session")
def open_profile():
    """Pore profile of 11ZC, the one downloaded entry with an open pore."""
    path = STRUCTURE_DIR / "11ZC.cif"
    if not path.exists():
        pytest.skip("11ZC.cif not downloaded — run python -m piezo1.io.fetch")
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    st = Structure.from_file(path)
    st = apply_frame(st, canonical_transform(st))
    blocks, _ = protomer_blocks(st)
    return pore_profile(st, detect_c3_axis(blocks), step=1.0)


@pytest.fixture(scope="session")
def qt_app():
    """A QApplication with the GL surface format already configured.

    The format has to be set before the first QOpenGLWidget exists, so it
    cannot be left to the window's own constructor in a test process that may
    already have made one.
    """
    pytest.importorskip("PyQt6.QtWidgets")
    from PyQt6.QtWidgets import QApplication

    from piezo1.ui.gl_widget import configure_surface_format

    app = QApplication.instance()
    if app is None:
        configure_surface_format()
        try:
            app = QApplication([])
        except Exception as exc:                       # pragma: no cover
            pytest.skip(f"no Qt platform available: {exc}")
    return app


@pytest.fixture(scope="session")
def structure_by_id():
    """Load any catalogued entry by id, or ``None`` when it is not downloaded.

    Returns ``None`` rather than skipping, so a test needing two structures can
    decide for itself whether a missing one is fatal.
    """
    cache: dict[str, object] = {}

    def load(pdb: str):
        if pdb not in cache:
            path = STRUCTURE_DIR / f"{pdb}.cif"
            cache[pdb] = Structure.from_file(path) if path.exists() else None
        return cache[pdb]

    return load


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
