"""Turning a computed current into a particle animation that does not lie.

A single PIEZO1 channel at its measured conductance passes on the order of
**10⁷ ions per second**. Nothing can be drawn at that rate: at 60 frames per
second it is 150,000 ions per frame, and any animation a viewer can follow is
therefore slowed by roughly a millionfold.

That is fine, and it is only fine if the display says so. An unlabelled stream
of particles reads as "this is what it looks like", which is wrong by six
orders of magnitude — the same class of error as a confidence interval quoted
where a model spread dominates.

So this module computes the time base rather than a speed: given the current
the permeation solver produced, how many ions per second that is, and what
slowdown makes it watchable. :class:`FluxTimebase` carries the statement the
HUD must show alongside the animation, in the spirit of the morph clock, which
reports a fraction along the path because a seconds axis would imply kinetics
the morph does not contain. Here the opposite holds — the current *is* a rate,
so a time base is meaningful and the honest thing is to name the factor.

**A closed pore animates nothing.** The permeation result is gated by the
wetting verdict, and a channel predicted not to conduct must show no particles
at all rather than a slow trickle. Drawing ions through a closed gate would be
a confidently wrong picture of the project's own central structural result.

**And it must say which constriction refused it.** 17 of the 19 deposited
PIEZO1 entries are refused, and in none of them is the narrowest point at the
transmembrane gate — it is beyond it, at the cytoplasmic constriction or above
it in the cap. Liu et al. 2025, whose intermediate-open entry is in this
catalogue, report that neck as bypassed by lateral portals an axial model does
not contain. "Sterically occluded" alone reads as a shut gate;
:mod:`piezo1.analysis.pore_regions` supplies the location and the gate's own
radius beside it.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["FluxTimebase", "ELEMENTARY_CHARGE", "ion_rate", "timebase",
           "timebase_for_structure", "DISPLAY_PARTICLES_PER_SECOND"]

#: Re-exported from :mod:`piezo1.physics.charge`, where they moved in Round 84d
#: so that `physics` and `analysis` could use them without importing `render`
#: and reversing the project's dependency arrow. Every existing import of
#: ``piezo1.render.flux.ion_rate`` keeps working.
from ..physics.charge import ELEMENTARY_CHARGE, ion_rate  # noqa: F401

#: Ions crossing per displayed second. Chosen so individual particles are
#: separable to the eye; it sets the slowdown rather than being derived from it.
DISPLAY_PARTICLES_PER_SECOND = 12.0


@dataclass(frozen=True)
class FluxTimebase:
    """How a real ionic current maps onto an animation a person can watch."""

    current_pA: float
    ions_per_second: float
    particles_per_second: float
    slowdown: float                 # real seconds represented per shown second
    conducting: bool = True
    reason: str = ""                # why not, when not conducting
    caveat: str = ""                # what the rate itself inherits

    @property
    def ions_per_particle(self) -> float:
        """One, by construction — each particle is one ion, shown slowly.

        The alternative convention, one particle standing for many ions at real
        speed, is the same arithmetic and a worse picture: a viewer counting
        particles would be counting the wrong thing.
        """
        return 1.0

    def statement(self) -> str:
        """The line the HUD shows. Never omitted, never abbreviated away."""
        if not self.conducting:
            return f"no conduction — {self.reason}"
        text = (f"{self.current_pA:.2f} pA = {self.ions_per_second:.2e} ions/s · "
                f"shown {self.slowdown:.0e}× slower than real time · "
                f"1 particle = 1 ion")
        return f"{text} · {self.caveat}" if self.caveat else text

    def frames_for(self, seconds: float, fps: float) -> int:
        return max(1, int(round(seconds * fps)))


def timebase(current_pA: float, valence: int = 1,
             particles_per_second: float | None = None,
             conducting: bool = True, reason: str = "",
             caveat: str = "") -> FluxTimebase:
    """The animation time base for a computed current.

    ``conducting=False`` produces a time base that animates nothing and carries
    the reason, so a caller cannot accidentally draw a stream through a pore
    the wetting model says is shut.
    """
    if particles_per_second is None:
        particles_per_second = DISPLAY_PARTICLES_PER_SECOND

    if not conducting or current_pA == 0.0:
        return FluxTimebase(
            current_pA=0.0, ions_per_second=0.0, particles_per_second=0.0,
            slowdown=1.0, conducting=False,
            reason=reason or "the pore is predicted not to conduct")

    ions = ion_rate(current_pA, valence)
    return FluxTimebase(
        current_pA=float(current_pA), ions_per_second=ions,
        particles_per_second=float(particles_per_second),
        slowdown=ions / particles_per_second, conducting=True, caveat=caveat)


def timebase_for_structure(structure, profile=None, grid=None,
                           pathway: str = "axial",
                           voltage: float | None = None) -> FluxTimebase:
    """The time base for a loaded structure, gated by the wetting verdict.

    Kept here rather than in the controller so it can be exercised on real
    coordinates without Qt. The first version lived in the controller and read
    a field that does not exist (``current_pA``); had it existed, the units
    would have been wrong by 10¹², because
    :class:`~piezo1.physics.permeation.PermeationResult` reports **amperes**.
    A test that only inspected the controller's source would not have caught
    either, and did not.
    """
    from ..analysis.hydration import load_grid, predict_wetting
    from ..analysis.pore_regions import describe_bottleneck
    from ..physics.permeation import default_species, solve_pnp
    from ..structure.pore import pore_profile
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis

    from ..physics.conduction_path import conduction_path

    if profile is None:
        blocks, _ = protomer_blocks(structure)
        profile = pore_profile(structure, detect_c3_axis(blocks))

    # `axial` returns the same object, so the default path is untouched.
    path = conduction_path(structure, profile, pathway)
    profile = path.profile

    try:
        verdict = predict_wetting(structure, profile,
                                  grid=grid if grid is not None else load_grid())
    except Exception as exc:
        return timebase(0.0, conducting=False,
                        reason=f"wetting not evaluated: {exc}")

    if verdict.hydrophobic_gate or verdict.sterically_occluded:
        # Say where. Every refusal in the catalogue is on a constriction BEYOND
        # the gate, which the axial model must pass and the real channel leaves
        # sideways before reaching — so a bare "sterically occluded" reads as a
        # shut gate and is wrong about which constriction it is.
        reason = verdict.summary()
        try:
            reason += " · " + describe_bottleneck(structure, profile).sentence()
        except Exception:                      # never lose the verdict itself
            pass
        if not path.is_axial or path.refused:
            reason += " · " + path.caveat()
        return timebase(0.0, conducting=False, reason=reason)

    try:
        # `wetting` is the second argument. This passed `default_species()`
        # there until Round 84d, so the solver's own blocking check silently
        # saw a list instead of a verdict and skipped — harmless only because
        # the caller had already gated on the same verdict two lines up.
        result = solve_pnp(profile, verdict, voltage=voltage,
                           species=default_species())
    except Exception as exc:
        return timebase(0.0, conducting=False, reason=f"permeation failed: {exc}")

    # The solver's conductance is a recorded DISAGREEMENT with the measurement
    # (about 41 pS against a published 25-30), reported rather than tuned. The
    # animation therefore runs faster than a real channel by that factor, and
    # says so — a stream calibrated to a number 1.5x too large is exactly the
    # kind of confident wrong picture Round 50 audited for.
    from ..parameters import PARAMETERS

    published = PARAMETERS.value("permeation.published_conductance")
    modelled = float(result.conductance) * 1e12
    caveat = ""
    if published > 0 and abs(modelled / published - 1.0) > 0.15:
        caveat = (f"rate from the model's {modelled:.0f} pS, "
                  f"{modelled / published:.1f}x the measured {published:.0f} pS")
    if not path.is_axial:
        caveat = (caveat + " · " if caveat else "") + path.caveat()

    # `.current` is in amperes; the time base is stated in picoamperes.
    return timebase(float(result.current) * 1e12, conducting=True, caveat=caveat)
