"""The helix detector and the cross-helices of Figure 7b.

The detector is calibrated on analytically generated helices of known rise,
radius and turn *before* it is pointed at a structure. The turn criterion is
what earns its place: rise and radius alone passed 41% of the windows of a
synthetic random walk, because a walk with a fixed step length looks locally
like a helix on both.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.structure.architecture import (IDEAL_HELIX_RADIUS,
                                           IDEAL_HELIX_RISE, IDEAL_HELIX_TURN,
                                           cross_helices, cross_helix_scan,
                                           helical_segments, helical_windows,
                                           ideal_helix)
from piezo1.structure.architecture import (_RADIUS_TOLERANCE, _RISE_TOLERANCE,
                                           _TURN_TOLERANCE)


def _passes(ca: np.ndarray) -> float:
    rise, radius, _, turn = helical_windows(ca)
    if not len(rise):
        return 0.0
    good = ((np.abs(rise - IDEAL_HELIX_RISE) <= _RISE_TOLERANCE)
            & (np.abs(radius - IDEAL_HELIX_RADIUS) <= _RADIUS_TOLERANCE)
            & (turn <= _TURN_TOLERANCE))
    return float(good.mean())


def test_the_estimator_is_unbiased_on_an_ideal_helix():
    """The textbook values, recovered — not absorbed into a tolerance.

    Taking the window's own principal axis instead of the sub-window centroids
    returns 1.63 A and 2.07 A here, and the tolerances would then have had to
    carry a systematic error as well as a real spread.
    """
    rise, radius, _, turn = helical_windows(ideal_helix(30))
    assert rise.mean() == pytest.approx(IDEAL_HELIX_RISE, abs=0.06)
    assert radius.mean() == pytest.approx(IDEAL_HELIX_RADIUS, abs=0.1)
    assert turn.max() < _TURN_TOLERANCE


def test_an_ideal_alpha_helix_passes_everywhere():
    assert _passes(ideal_helix(30)) == 1.0


@pytest.mark.parametrize("label,rise,radius,turn", [
    ("3-10 helix", 2.0, 1.9, 120.0),
    ("pi helix", 1.1, 2.8, 87.0),
    ("beta strand", 3.3, 0.9, 180.0),
    ("left-handed alpha", IDEAL_HELIX_RISE, IDEAL_HELIX_RADIUS, -100.0),
])
def test_other_regular_geometries_are_rejected(label, rise, radius, turn):
    """Each fails at least one criterion — a filter that admitted them all
    would not be measuring 'is this an alpha helix'."""
    assert _passes(ideal_helix(30, rise=rise, radius=radius, turn=turn)) == 0.0


@pytest.mark.parametrize("seed", [0, 1, 2, 3])
def test_a_random_coil_is_rejected(seed):
    """The case the turn criterion was added for."""
    rng = np.random.default_rng(seed)
    step = rng.normal(size=(60, 3))
    step /= np.linalg.norm(step, axis=1, keepdims=True)
    coil = np.cumsum(step * 3.8, axis=0)
    assert _passes(coil) < 0.05


def test_the_turn_criterion_is_doing_the_work():
    """Without it, the coil would pass a substantial fraction of windows."""
    rng = np.random.default_rng(0)
    step = rng.normal(size=(200, 3))
    step /= np.linalg.norm(step, axis=1, keepdims=True)
    coil = np.cumsum(step * 3.8, axis=0)
    rise, radius, _, turn = helical_windows(coil)
    without = ((np.abs(rise - IDEAL_HELIX_RISE) <= _RISE_TOLERANCE)
               & (np.abs(radius - IDEAL_HELIX_RADIUS) <= _RADIUS_TOLERANCE))
    with_turn = without & (turn <= _TURN_TOLERANCE)
    assert without.mean() > 0.15, (
        "if rise and radius alone already rejected a coil, the turn criterion "
        "would be untested decoration")
    assert with_turn.mean() < 0.05


def test_the_ideal_helix_generator_is_what_it_claims():
    """Calibrating against a generator nobody checked would prove nothing."""
    helix = ideal_helix(20)
    rise = np.diff(helix[:, 2])
    assert np.allclose(rise, IDEAL_HELIX_RISE)
    assert np.allclose(np.linalg.norm(helix[:, :2], axis=1),
                       IDEAL_HELIX_RADIUS)
    angle = np.degrees(np.arctan2(helix[1, 1], helix[1, 0]))
    assert angle == pytest.approx(IDEAL_HELIX_TURN, abs=1e-6)


# --------------------------------------------------------------------------
# On the real structure
# --------------------------------------------------------------------------

def test_the_beam_is_found_as_one_helix_not_two(structure_6b3r):
    """The beam bends, fails the turn test for a window, and came back split.

    A bent helix arriving as two overlapping segments would double-count it in
    every figure and put a spurious short helix beside it.
    """
    segments = helical_segments(structure_6b3r, "A", 1290, 1370)
    assert len(segments) == 1, [(s.start, s.end) for s in segments]
    beam = segments[0]
    assert beam.n_residues > 60, "Guo & MacKinnon call it 66 amino acids"
    assert beam.rise == pytest.approx(IDEAL_HELIX_RISE, abs=0.2)
    assert beam.radius == pytest.approx(IDEAL_HELIX_RADIUS, abs=0.3)


def test_the_beam_is_not_reported_as_a_cross_helix(structure_6b3r):
    """Figure 7b colours it red and the cross-helices yellow."""
    found = cross_helices(structure_6b3r, "mouse")
    assert found, "some cross-helices should be found"
    for segment in found:
        assert not (segment.start <= 1365 and segment.end >= 1300), (
            f"the beam ({segment.start}-{segment.end}) must be excluded")


def test_every_protomer_finds_the_same_cross_helices(structure_6b3r):
    """A C3-symmetric structure that did not would mean the detector is
    picking up noise rather than architecture."""
    found = cross_helices(structure_6b3r, "mouse")
    per_chain: dict[str, set] = {}
    for segment in found:
        per_chain.setdefault(segment.chain, set()).add(
            (segment.start, segment.end))
    assert len(per_chain) == 3
    assert len(set(map(frozenset, per_chain.values()))) == 1


def test_the_paper_says_at_least_one_per_linker_and_that_holds(structure_6b3r):
    found = cross_helices(structure_6b3r, "mouse")
    assert len(found) % 3 == 0
    assert len(found) // 3 >= 1


def test_the_threshold_is_reported_against_a_scan_not_asserted(structure_6b3r):
    """The count must vary smoothly with the cut, so a reader can see how
    sharp it is rather than trusting 55 degrees."""
    scan = cross_helix_scan(structure_6b3r, "mouse")
    assert scan.counts[0] >= scan.counts[-1]
    assert scan.counts == tuple(sorted(scan.counts, reverse=True))
    assert scan.linker_tilts and scan.transmembrane_tilts
    assert (np.median(scan.linker_tilts)
            > np.median(scan.transmembrane_tilts) + 15.0), (
        "'perpendicular to the TM helices' must be a real distinction")
