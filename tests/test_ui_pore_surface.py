"""The drawn pore, and the three ways a picture of it could be confidently wrong.

Nothing here computes a pore. The controller reads the profile the Analysis
panel produced and draws its spheres, so these tests are about the join — which
is where the equivalent feature went wrong last time: `interaction_controller`
keyed its colours on a string the analysis does not emit and silently lost five
sixths of what it drew.

Three failure modes are checked rather than assumed:

* the bands could be applied in the wrong order, which is exactly what
  `plddt_band_colors` did until Round 76 — every atom painted the first band's
  colour, and the picture looked deliberate;
* the drawn spheres could be a smoothed or resampled version of the measured
  ones, so identity is asserted, not closeness;
* the caption could describe a different run from the picture, so the counts
  come through the same call the colours do.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.parameters import PARAMETERS as _P  # noqa: E402
from piezo1.ui.pore_controller import (BAND_COLORS, BAND_NAMES,  # noqa: E402
                                       LINING_COLOR, NARROW_COLOR,
                                       OPEN_COLOR, PoreSurfaceController,
                                       TIGHT_COLOR, band_index, radius_colors)


# --------------------------------------------------------------------------
# The bands, calibrated on values whose band is known by inspection
# --------------------------------------------------------------------------

def test_each_radius_lands_in_the_band_it_belongs_to():
    ion = _P.value("pore.ion_radius")
    threshold = _P.value("pore.constriction_threshold")
    assert ion < threshold, "the bands assume the ion radius is the tighter cut"

    below = ion / 2.0
    between = 0.5 * (ion + threshold)
    above = threshold * 2.0
    assert list(band_index([below, between, above])) == [0, 1, 2]


def test_the_band_edges_belong_to_the_band_they_open():
    """A radius exactly at a threshold clears it. Stated, because the
    alternative is defensible and the two differ by one slice at the gate."""
    ion = _P.value("pore.ion_radius")
    threshold = _P.value("pore.constriction_threshold")
    assert list(band_index([ion, threshold])) == [1, 2]


def test_bands_are_applied_ascending_so_each_takes_the_highest_it_clears():
    """The Round 76 defect, in the one place it could recur.

    `plddt_band_colors` applied its bands in declaration order, so the
    ``>= 0.0`` pass repainted every atom and the whole model came out the
    lowest band's colour. Applied that way here, a wide-open pore would be
    entirely red — a perfectly plausible picture of a shut channel.
    """
    wide = np.full(50, _P.value("pore.constriction_threshold") * 3.0)
    colours = radius_colors(wide)
    assert np.allclose(colours, np.asarray(OPEN_COLOR, np.float32)), \
        "a wide pore came out in a narrower band's colour — bands are being " \
        "applied in the wrong order"
    narrow = np.full(50, 0.1)
    assert np.allclose(radius_colors(narrow),
                       np.asarray(NARROW_COLOR, np.float32))


def test_the_three_bands_are_three_distinct_colours():
    """Two bands sharing a colour is a band that cannot be seen."""
    assert len({NARROW_COLOR, TIGHT_COLOR, OPEN_COLOR}) == 3
    assert len(BAND_COLORS) == len(BAND_NAMES) == 3
    assert LINING_COLOR not in BAND_COLORS, \
        "the lining marker must not be mistakable for a probe sphere"


def test_the_bands_follow_the_registry_rather_than_a_literal():
    """An override in the parameters dialog must move the picture."""
    original = _P.value("pore.ion_radius")
    radius = [original * 1.5]
    assert int(band_index(radius)[0]) == 1
    try:
        _P.set_value("pore.ion_radius", original * 2.0)
        assert int(band_index(radius)[0]) == 0, \
            "raising the ion radius did not reclassify a slice that is now " \
            "too narrow for it"
    finally:
        _P.reset()
    assert int(band_index(radius)[0]) == 1


# --------------------------------------------------------------------------
# What is drawn is what was measured
# --------------------------------------------------------------------------

class _FakeWindow:
    """Just enough window for the controller's read-only paths."""

    def __init__(self, structure, profile, hydration=None):
        self.structure = structure
        self.analysis = type("A", (), {"pore": profile,
                                       "hydration": hydration})()
        self.status = ""

    def _set_status(self, text):
        self.status = text


@pytest.fixture(scope="module")
def controller(open_profile):
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure
    from piezo1.structure.frame import apply_frame, canonical_transform

    path = STRUCTURE_DIR / "11ZC.cif"
    if not path.exists():
        pytest.skip("11ZC.cif not downloaded — run python -m piezo1.io.fetch")
    st = Structure.from_file(path)
    st = apply_frame(st, canonical_transform(st))
    win = _FakeWindow(st, open_profile)
    return PoreSurfaceController(win)


def test_the_controller_draws_the_analysis_profile_not_a_copy(controller,
                                                              open_profile):
    """Identity, not equality. Two pore profiles of the same structure that
    differ by a parameter override would both look reasonable on screen."""
    assert controller.profile is open_profile


def test_every_drawn_sphere_actually_fits_where_it_is_drawn(controller):
    """The calibration: a sphere is only the pore radius if it is empty.

    A picture drawn from the wrong frame, a stale profile or a resampled
    centre gives spheres that sit inside the protein — and a radius plot
    cannot show that, which is the whole reason this feature exists. Checked
    against the van der Waals radii the profiler itself measures clearance
    with, so a sphere overlapping an atom is a real failure and not a
    convention mismatch.
    """
    from scipy.spatial import cKDTree

    structure = controller.win.structure
    profile = controller.profile
    tree = cKDTree(structure.xyz)
    vdw = structure.vdw_radii()

    def worst_overlap(centers, radii):
        worst, examined = 0.0, 0
        for centre, radius in zip(centers, radii):
            near = tree.query_ball_point(centre, float(radius) + vdw.max())
            if not near:
                continue
            examined += 1
            clearance = (np.linalg.norm(structure.xyz[near] - centre, axis=1)
                         - vdw[near]).min()
            worst = max(worst, float(radius) - clearance)
        return worst, examined

    overlap, examined = worst_overlap(profile.centers, profile.radius)
    assert examined > 20, \
        f"only {examined} spheres had any atom near them; this proves nothing"
    assert overlap < 1e-6, (
        f"a drawn probe sphere overlaps an atom by {overlap:.3f} A — the "
        f"spheres are not the ones the profile was measured with")

    # The same check on centres moved 5 A off the axis, so a pass above is a
    # measurement rather than a check that cannot fail.
    shifted, _ = worst_overlap(profile.centers + np.array([5.0, 0.0, 0.0]),
                               profile.radius)
    assert shifted > 1.0, \
        "moving every probe 5 A sideways did not put one inside an atom, so " \
        "this test would pass on a picture drawn in the wrong place"


def test_the_lining_marker_is_the_bottleneck_residues_in_every_protomer(controller):
    """Three copies, because the constriction is on the three-fold axis.

    Marking one protomer's residues would put the gate visually off-axis,
    which is a claim about asymmetry that nothing measured.
    """
    lining = controller.profile.bottleneck_lining()
    if not lining:
        pytest.skip("this profile records no lining residues")
    coords = controller.lining_coords()
    structure = controller.win.structure
    expected = int((structure.mask_ca()
                    & np.isin(structure.res_seq, np.asarray(lining))).sum())
    assert len(coords) == expected
    assert len(coords) >= 3 * len(lining) - 2, \
        "the bottleneck lining was marked in fewer copies than there are " \
        "protomers resolving it"


# --------------------------------------------------------------------------
# The caption
# --------------------------------------------------------------------------

def test_the_counts_on_the_status_line_come_from_the_drawn_colours(controller):
    counts = controller.band_counts()
    colours = radius_colors(controller.profile.radius)
    for i, name in enumerate(BAND_NAMES):
        drawn = int(np.all(np.isclose(colours, np.asarray(BAND_COLORS[i],
                                                          np.float32)),
                           axis=1).sum())
        assert counts[name] == drawn, \
            f"the caption says {counts[name]} {name} slices and {drawn} are " \
            f"drawn that colour"
    assert sum(counts.values()) == len(controller.profile.radius)


def test_the_status_line_refuses_to_let_a_radius_stand_for_a_verdict(controller):
    """Both caveats are load-bearing and neither is optional.

    A probe sphere is the space left over, not the wall; and Round 19's whole
    point is that a lumen wide enough for an ion can still be shut. Without
    the second, a blue pore reads as an open channel.
    """
    line = controller.status_line().lower()
    assert "not the pore wall" in line
    assert "wetting" in line or "dewet" in line
    assert "bottleneck" in line


def test_an_empty_controller_says_so_rather_than_raising():
    win = _FakeWindow(None, None)
    controller = PoreSurfaceController(win)
    assert controller.band_counts() == {}
    assert "no pore profile" in controller.status_line()
    assert len(controller.lining_coords()) == 0
