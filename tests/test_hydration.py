"""Pore hydration and the hydrophobic-gating heuristic.

The point of this round is that **radius alone does not decide conduction**, so
the tests that matter are the ones that could not be passed by a radius
threshold in disguise. The controls below hold the geometry fixed and vary only
the chemistry.

Everything needs the CHAP grid, which is downloaded rather than committed;
tests skip rather than fail without it.
"""

import numpy as np
import pytest

from piezo1.analysis import hydration as H
from piezo1.core.annotations import load_annotations
from piezo1.structure.pore import pore_profile
from piezo1.structure.superpose import detect_c3_axis
from conftest import protomer_blocks


@pytest.fixture(scope="module")
def grid():
    g = H.load_grid()
    if not g.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")
    return g


def _profile(st):
    return pore_profile(st, detect_c3_axis(protomer_blocks(st)[0]), step=1.0)


@pytest.fixture(scope="module")
def closed_profile(human_structure):
    return _profile(human_structure)


@pytest.fixture(scope="module")
def closed(human_structure, closed_profile, grid):
    return H.predict_wetting(human_structure, closed_profile, grid)


# --------------------------------------------------------------------------
# The published landscape, read correctly
# --------------------------------------------------------------------------

def test_grid_axes_match_the_publication(grid):
    assert grid.energy.shape == (100, 100)
    assert grid.radius.min() == pytest.approx(0.10, abs=1e-6)
    assert grid.radius.max() == pytest.approx(0.60, abs=1e-6)
    assert np.isfinite(grid.energy).all(), "grid did not reshape cleanly"


def test_contour_reproduces_the_published_critical_radii(grid):
    """Rao et al.: hydrophilic pores wet below 0.2 nm, hydrophobic ones can
    hold a barrier out to ~0.4 nm. Recovering that from the grid is the check
    that our contour extraction indexes the landscape the right way round."""
    contour = grid.contour()
    assert len(contour) > 20
    hydrophilic = contour[contour[:, 0] < -0.1]
    hydrophobic = contour[contour[:, 0] > 0.2]
    assert hydrophilic[:, 1].max() < 0.2
    assert 0.35 < hydrophobic[:, 1].max() < 0.45
    # Critical radius must rise with hydrophobicity, which is the whole model.
    assert np.corrcoef(contour[:, 0], contour[:, 1])[0, 1] > 0.9


def test_energy_falls_with_radius_and_rises_with_hydrophobicity(grid):
    mid = grid.energy.shape[0] // 2
    assert grid.energy[mid, 0] > grid.energy[mid, -1]
    narrow = grid.energy[:, 5]
    assert narrow[-1] > narrow[0]


def test_units_are_converted_at_the_boundary(grid):
    """The grid is nm; this project is Angstrom. A missed conversion is a
    factor of ten and would put every pore off the top of the landscape."""
    at_3A = float(np.atleast_1d(grid.energy_at(0.0, 3.0))[0])
    at_3nm = float(np.atleast_1d(grid.energy_at(0.0, 30.0))[0])
    assert at_3A > at_3nm
    assert at_3nm == pytest.approx(grid.energy[np.abs(grid.hydrophobicity).argmin(), -1])


# --------------------------------------------------------------------------
# Hydrophobicity profile
# --------------------------------------------------------------------------

def test_hydrophobicity_spans_the_grid_range(human_structure, closed_profile):
    """Smoothing in 3-D around the probe centre instead of along the pore axis
    collapses this to a narrow band near zero (-0.12 to +0.02 when it was
    wrong), which then indexes the landscape with the wrong coordinate. A
    profile that does not use most of the grid's range is the symptom."""
    h = H.hydrophobicity_profile_chap(human_structure, closed_profile)
    ok = np.isfinite(h)
    assert ok.sum() > 20
    assert h[ok].min() < -0.3
    assert h[ok].max() > 0.15


def test_hydrophobicity_uses_the_wimley_white_scale():
    assert len(H.WIMLEY_WHITE_NORMALISED) == 20
    assert H.WIMLEY_WHITE_NORMALISED["ASP"] == -1.0
    assert H.WIMLEY_WHITE_NORMALISED["PHE"] > H.WIMLEY_WHITE_NORMALISED["SER"]
    assert H.WIMLEY_WHITE_NORMALISED["ILE"] == H.WIMLEY_WHITE_NORMALISED["LEU"]
    assert max(H.WIMLEY_WHITE_NORMALISED.values()) <= 1.0
    assert min(H.WIMLEY_WHITE_NORMALISED.values()) >= -1.0


def test_wider_bandwidth_smooths_more(human_structure, closed_profile):
    tight = H.hydrophobicity_profile_chap(human_structure, closed_profile,
                                          bandwidth_nm=0.15)
    broad = H.hydrophobicity_profile_chap(human_structure, closed_profile,
                                          bandwidth_nm=1.0)
    assert np.nanstd(tight) > np.nanstd(broad)


# --------------------------------------------------------------------------
# The result the round asked for
# --------------------------------------------------------------------------

def test_closed_structure_is_called_non_conductive(closed):
    assert closed.score > H.CLOSED_SCORE_CUTOFF
    assert closed.hydrophobic_gate
    assert not closed.conductive


def test_flat_structure_is_called_conductive(grid):
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure
    path = STRUCTURE_DIR / "11ZC.cif"
    if not path.exists():
        pytest.skip("11ZC not downloaded")
    st = Structure.from_file(path)
    pred = H.predict_wetting(st, _profile(st), grid)
    assert pred.score <= H.CLOSED_SCORE_CUTOFF
    assert not pred.hydrophobic_gate
    assert not pred.sterically_occluded
    assert pred.conductive


def test_the_verdict_is_hydrophobicity_not_narrowness(human_structure,
                                                      closed_profile, grid,
                                                      monkeypatch):
    """The control that makes this round worth doing.

    Hold every radius fixed and replace the hydrophobicity scale with a uniform
    hydrophilic value. If the closed call were really a radius threshold in
    disguise, the score would not move. It collapses to zero.
    """
    monkeypatch.setattr(H, "WIMLEY_WHITE_NORMALISED",
                        {k: -1.0 for k in H.WIMLEY_WHITE_NORMALISED})
    relabelled = H.predict_wetting(human_structure, closed_profile, grid)
    assert relabelled.score == pytest.approx(0.0, abs=1e-9)
    assert not relabelled.hydrophobic_gate


def test_dewetting_is_flagged_where_water_would_still_fit(closed):
    """Hydrophobic gating means blocked *without* steric occlusion.

    If every flagged residue were narrower than a water molecule the heuristic
    would be reporting sterics under another name.
    """
    passable = [p for p in closed.dewetted
                if p.radius / 10.0 > H.WATER_RADIUS_NM]
    assert passable, "no dewetted residue is sterically passable"
    assert max(p.radius for p in passable) / 10.0 > 0.3


def test_flagged_residues_recover_the_curated_gate(closed):
    """Independent check: the heuristic sees only coordinates and a
    hydrophobicity scale, never the annotation. It should still land on the
    residues curated as the gate and the cytoplasmic constrictions."""
    ann = load_annotations("human")
    curated: set[int] = set()
    for group in ann.residue_groups:
        label = group.label.lower()
        if "gate" in label or "pore-lining" in label or "constriction" in label:
            curated |= set(group.residues)
    flagged = {p.residue for p in closed.dewetted}
    assert flagged & curated, f"flagged {sorted(flagged)} misses {sorted(curated)}"


def test_matched_radius_opposite_verdict(closed, grid):
    """The cleanest statement of hydrophobic gating available in these data.

    The closed structure carries dewetted residues at ~0.32 nm, essentially the
    same radius as the flat structure's *bottleneck* (0.330 nm), which is
    called wet. Same radius, opposite verdict — decided by hydrophobicity.
    """
    wide = [p for p in closed.dewetted if p.radius / 10.0 > 0.30]
    assert wide, "expected dewetted residues at ~0.32 nm in the closed state"
    assert all(p.hydrophobicity > 0.1 for p in wide)


# --------------------------------------------------------------------------
# Honest limits
# --------------------------------------------------------------------------

def test_steric_and_hydrophobic_closure_are_reported_separately(grid):
    """A limitation found by testing more than the two structures asked for.

    7WLU and 8IXO have 0.098 nm bottlenecks — far too narrow for water — but
    hydrophilic linings, so the Rao score alone calls them open. The heuristic
    answers "would water dewet here?", not "does water fit here?". Merging the
    two into one verdict would hide that, so both are exposed.
    """
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure
    path = STRUCTURE_DIR / "7WLU.cif"
    if not path.exists():
        pytest.skip("7WLU not downloaded")
    st = Structure.from_file(path)
    pred = H.predict_wetting(st, _profile(st), grid)
    assert not pred.hydrophobic_gate, "7WLU has no hydrophobic gate"
    assert pred.sterically_occluded, "7WLU is far too narrow for water"
    assert not pred.conductive
    assert "sterically occluded" in pred.verdict


def test_missing_grid_degrades_rather_than_raises(human_structure,
                                                  closed_profile):
    empty = H.HydrationGrid(np.array([]), np.array([]), np.array([[]]),
                            source="not downloaded")
    pred = H.predict_wetting(human_structure, closed_profile, empty)
    assert not pred.available
    assert "unavailable" in pred.verdict
    assert np.isnan(pred.score)


def test_provenance_is_recorded(closed):
    assert "31235590" in H.RAO_CITATION
    assert "MIT" in H.CHAP_CITATION
    assert closed.meta["threshold_kJ"] == 2.6
    assert closed.meta["cutoff"] == 0.55
    assert "Rao" in closed.meta["citation"]
