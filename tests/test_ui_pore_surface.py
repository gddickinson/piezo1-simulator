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
                                       TIGHT_COLOR, band_index,
                                       drawn_slice_mask, radius_colors)


# --------------------------------------------------------------------------
# The trim, calibrated on profiles whose right answer is known by construction
# --------------------------------------------------------------------------

def _leash() -> float:
    return _P.value("pore.leash")


def test_the_escaped_ends_are_trimmed_and_only_the_ends():
    """The shape the real profiles have: a lumen with the probe ballooning
    out of both mouths once nothing but the tether is bounding it."""
    leash = _leash()
    radius = np.array([3 * leash, 2 * leash] + [2.0] * 10 + [2 * leash])
    mask = drawn_slice_mask(radius)
    assert list(mask) == [False, False] + [True] * 10 + [False]


def test_a_wide_interior_is_kept():
    """The case that killed the first version of this rule.

    Keeping the contiguous in-leash run around the bottleneck cut **71**
    slices off 11ZC's upper vestibule to remove 5 at the bottom. A wide slice
    with protein on both sides of it is a vestibule, not an escape, and the
    profile alone cannot tell them apart — but their *position* can.
    """
    leash = _leash()
    radius = np.array([3 * leash] + [2.0] * 4 + [2 * leash] * 6 + [2.0] * 4)
    mask = drawn_slice_mask(radius)
    assert not mask[0]
    assert mask[1:].all(), "an over-leash vestibule in the middle was trimmed"


def test_a_profile_that_never_escapes_is_drawn_whole():
    """Every closed entry is like this, so the trim must be able to do nothing."""
    radius = np.full(40, 2.0)
    assert drawn_slice_mask(radius).all()


def test_a_profile_that_is_entirely_escaped_is_drawn_whole_rather_than_emptied():
    """There is no end to trim towards, and an empty picture with no
    explanation is worse than a wide one with a caption."""
    assert drawn_slice_mask(np.full(20, 5 * _leash())).all()


def test_the_bottleneck_can_never_be_trimmed_away():
    """The invariant that makes this display-only.

    A trimmed slice is wider than the leash, so it cannot be the minimum
    unless every slice is — and then nothing is trimmed. Driven over random
    profiles rather than argued, because the argument is the kind that is
    right until the rule changes.
    """
    rng = np.random.default_rng(11)
    leash = _leash()
    for _ in range(200):
        radius = rng.uniform(0.5, 3 * leash, size=rng.integers(5, 60))
        mask = drawn_slice_mask(radius)
        assert mask[int(np.argmin(radius))], \
            "the narrowest slice — the whole answer — was trimmed off"


def test_the_trim_follows_the_registry():
    radius = np.array([_leash() * 1.5] + [1.0] * 5)
    assert not drawn_slice_mask(radius)[0]
    original = _P.value("pore.leash")
    try:
        _P.set_value("pore.leash", original * 2.0)
        assert drawn_slice_mask(radius).all(), \
            "doubling the leash did not bring the escaped end back inside it"
    finally:
        _P.reset()


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
    keep = controller.drawn_mask()
    # The colours of the spheres that reach the screen, which is what the
    # caption is describing — the trimmed ends are neither drawn nor counted.
    colours = radius_colors(np.asarray(controller.profile.radius)[keep])
    for i, name in enumerate(BAND_NAMES):
        drawn = int(np.all(np.isclose(colours, np.asarray(BAND_COLORS[i],
                                                          np.float32)),
                           axis=1).sum())
        assert counts[name] == drawn, \
            f"the caption says {counts[name]} {name} slices and {drawn} are " \
            f"drawn that colour"
    assert sum(counts.values()) == int(keep.sum())


def test_the_open_entry_loses_its_bulb_and_keeps_its_vestibule(controller):
    """The measurement on the entry the trim was written for.

    11ZC's probe balloons to 12.2 A below the channel — a bulb hanging under
    the protein — while the wide part *above* the gate is a real vestibule
    with protein around it. Both facts are pinned, because a rule that
    removed the second to remove the first would be worse than no rule.
    """
    profile = controller.profile
    keep = controller.drawn_mask()
    radius = np.asarray(profile.radius)
    dropped = ~keep
    if not dropped.any():
        pytest.skip("this profile never escapes; nothing to say about trimming")

    assert radius[dropped].min() > _P.value("pore.leash")
    assert dropped.sum() < 0.25 * len(radius), (
        f"{int(dropped.sum())} of {len(radius)} slices were trimmed; that is "
        f"a rule removing the pore rather than its escaped ends")
    # Only the ends: the kept slices must be one contiguous block.
    kept = np.flatnonzero(keep)
    assert np.array_equal(kept, np.arange(kept[0], kept[-1] + 1))
    # And the wide interior survives.
    interior = radius[kept] > _P.value("pore.leash")
    assert interior.sum() > 0, (
        "no over-leash slice survived, so this test cannot tell the "
        "end-trimming rule from the contiguous-run rule it replaced")


def test_trimming_does_not_touch_the_measurement(controller):
    """Display only. The profile is the Analysis panel's object and the plot
    is drawn from it, so a trim that mutated it would silently change the
    number the picture is meant to illustrate."""
    profile = controller.profile
    before = (profile.bottleneck_radius, profile.bottleneck_z,
              len(profile.z), float(np.asarray(profile.radius).sum()))
    controller.drawn_mask()
    controller.status_line()
    after = (profile.bottleneck_radius, profile.bottleneck_z,
             len(profile.z), float(np.asarray(profile.radius).sum()))
    assert before == after


def test_the_status_line_says_what_was_trimmed_and_that_it_changed_nothing(
        controller):
    line = controller.status_line()
    if (~controller.drawn_mask()).any():
        assert "NOT drawn" in line
        assert "bulk solvent" in line
        assert "bottleneck are unchanged" in line
        assert f"of {len(controller.profile.z)} probe spheres" in line


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
