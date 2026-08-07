"""The particle animation, and the two ways it could lie.

Round 33 built the permeation physics and deferred the animation. The physics
is the easy part: a channel passes ~10^7 ions per second, so any watchable
stream is about a millionfold slow, and an unlabelled stream would misstate the
rate by six orders of magnitude.

The second way it could lie is worse. Every deposited human structure is closed,
and the permeation result is gated by the wetting verdict — so a trickle of ions
through a shut gate would contradict the project's own structural result while
looking like a demonstration of it.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.render.flux import (DISPLAY_PARTICLES_PER_SECOND,
                                ELEMENTARY_CHARGE, FluxTimebase, ion_rate,
                                timebase)


# ------------------------------------------------------------ the rate

def test_the_ion_rate_matches_the_definition_of_current():
    """Calibrated against the arithmetic, not against itself.

    One ampere is one coulomb per second, so 1 pA is 1e-12 / e ions per second
    for a monovalent ion. Checked at a round number before the real current is
    trusted through it.
    """
    assert ion_rate(1.0) == pytest.approx(1e-12 / ELEMENTARY_CHARGE, rel=1e-12)
    assert ion_rate(1.0) == pytest.approx(6.2415e6, rel=1e-4)


def test_a_divalent_ion_carries_the_same_current_with_half_the_ions():
    """Matters because PIEZO1 is calcium-permeable."""
    assert ion_rate(2.0, valence=2) == pytest.approx(ion_rate(2.0, valence=1) / 2)


def test_a_nonsensical_valence_is_refused():
    with pytest.raises(ValueError, match="valence"):
        ion_rate(1.0, valence=0)


def test_the_current_a_channel_actually_passes_needs_a_millionfold_slowdown():
    """The measured fact that makes the label non-optional."""
    base = timebase(1.65)
    assert base.ions_per_second == pytest.approx(1.03e7, rel=0.02)
    assert 1e5 < base.slowdown < 1e7, base.slowdown
    assert base.particles_per_second == DISPLAY_PARTICLES_PER_SECOND


# ----------------------------------------------------- the honest label

def test_the_statement_names_the_slowdown_and_the_particle_meaning():
    text = timebase(1.65).statement()
    assert "ions/s" in text
    assert "slower than real time" in text
    assert "1 particle = 1 ion" in text


def test_one_particle_is_one_ion_by_construction():
    """The alternative convention would make a viewer count the wrong thing."""
    assert timebase(1.65).ions_per_particle == 1.0


def test_the_slowdown_scales_with_the_current():
    """Twice the current, twice as slow — the label must track the physics."""
    slow_1 = timebase(1.0).slowdown
    slow_2 = timebase(2.0).slowdown
    assert slow_2 == pytest.approx(2 * slow_1, rel=1e-9)


# --------------------------------------------- a shut pore animates nothing

def test_a_non_conducting_pore_produces_no_particles():
    base = timebase(3.0, conducting=False, reason="hydrophobic gate")
    assert not base.conducting
    assert base.particles_per_second == 0.0
    assert base.ions_per_second == 0.0
    assert base.current_pA == 0.0, (
        "a non-conducting time base must not carry a current a caller could draw")


def test_a_shut_pore_states_why_instead_of_showing_a_rate():
    base = timebase(3.0, conducting=False, reason="sum d = 0.90 > 0.55")
    text = base.statement()
    assert "no conduction" in text
    assert "0.90" in text
    assert "slower than real time" not in text, (
        "a shut pore must not advertise a time base it is not using")


def test_zero_current_is_treated_as_not_conducting_even_if_claimed_open():
    """Defensive: a solver returning zero must not divide into a slowdown."""
    base = timebase(0.0, conducting=True)
    assert not base.conducting
    assert base.slowdown == 1.0


def test_a_default_reason_is_supplied_rather_than_left_blank():
    assert timebase(0.0, conducting=False).reason


# ------------------------------------------------------------ the controller

def test_the_controller_uses_the_shared_gated_computation():
    """The GUI must not reach its own conclusion about whether this conducts.

    This began as a test that read the controller's source for
    ``predict_wetting``. That passed while the units were wrong by 10^12, and
    then broke the moment the gating moved somewhere testable — which is the
    argument against source-inspection tests. It now asserts that the
    controller calls the shared function, and the *behaviour* of that function
    is checked below on real coordinates.
    """
    import inspect

    from piezo1.ui import ion_flux_controller

    source = inspect.getsource(ion_flux_controller.IonFluxController._prepare)
    assert "timebase_for_structure" in source, (
        "the controller must use the shared gated computation, not its own")


def test_the_controller_is_reachable_from_the_view_menu():
    from pathlib import Path

    menus = Path(__file__).resolve().parents[1] / "piezo1" / "ui" / "menus.py"
    text = menus.read_text()
    assert "ion_flux.show" in text
    assert "MILLIONFOLD" in text.upper(), (
        "the menu tooltip must warn about the time base before it is switched on")


def test_the_particles_are_not_coloured_like_any_element():
    """They represent a rate, not resolved atoms."""
    from piezo1.core.structure import ELEMENT_COLORS
    from piezo1.ui.ion_flux_controller import ION_COLOR

    for colour in ELEMENT_COLORS.values():
        assert not np.allclose(np.asarray(colour, float)[:3],
                               np.asarray(ION_COLOR, float), atol=0.02)


# ------------------------------- the real structures, on real coordinates

def _structure(pdb):
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure

    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb} not downloaded")
    return Structure.from_file(path)


def test_the_closed_human_structure_animates_nothing():
    """8YEZ is shut on both criteria; the animation must agree with the pipeline."""
    from piezo1.render.flux import timebase_for_structure

    base = timebase_for_structure(_structure("8YEZ"))
    assert not base.conducting
    assert base.particles_per_second == 0.0
    assert "non-conductive" in base.reason or "occluded" in base.reason


def test_the_open_like_structure_gives_a_watchable_rate():
    """11ZC conducts, and the units are the ones that bit on the first attempt.

    `PermeationResult.current` is in AMPERES. The controller's first version
    read a `current_pA` field that does not exist; had it existed and been
    misread, the rate would have been wrong by 10^12 and the animation would
    have looked plausible either way.
    """
    from piezo1.render.flux import timebase_for_structure

    base = timebase_for_structure(_structure("11ZC"))
    assert base.conducting
    assert 1.0 < base.current_pA < 10.0, (
        f"{base.current_pA} pA is not a single-channel current — check units")
    assert 1e6 < base.ions_per_second < 1e8
    assert 1e4 < base.slowdown < 1e8


def test_the_rate_declares_the_conductance_disagreement_it_inherits():
    """The solver records 41 pS against a published 25-30, reported not tuned.

    An animation calibrated to a number 1.5x too large, shown without saying
    so, is the confident-wrong-picture failure Round 50 audited for.
    """
    from piezo1.parameters import PARAMETERS
    from piezo1.render.flux import timebase_for_structure

    base = timebase_for_structure(_structure("11ZC"))
    assert base.caveat, "the inherited overestimate must be stated"
    assert "pS" in base.caveat
    published = PARAMETERS.value("permeation.published_conductance")
    assert f"{published:.0f}" in base.caveat
    assert base.caveat in base.statement()


def test_a_rate_that_matched_the_measurement_would_carry_no_caveat():
    """Calibration: the caveat must be able to be absent, or it says nothing."""
    from piezo1.render.flux import timebase

    assert timebase(1.65).caveat == ""
