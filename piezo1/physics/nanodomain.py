"""The calcium nanodomain a channel makes, and what a tag sitting in it sees.

An open channel is a point source of calcium. Diffusion carries it away and
cytosolic buffers absorb it, and the balance is a standing gradient that reaches
steady state in microseconds — far faster than a channel stays open, so the
static solution is the right one.

.. math::

    [\\mathrm{Ca}](r) = \\frac{i_\\mathrm{Ca}}{4\\pi F D r}\\,e^{-r/\\lambda},
    \\qquad \\lambda = \\sqrt{\\frac{D}{k_\\mathrm{on}^{B}[B]}}

**Where the 4π comes from**, since it is the one place a factor of two hides:
the calcium *ion* flux is :math:`i_\\mathrm{Ca}/(zF)` with :math:`z = 2`, and a
channel in a membrane releases into a half-space, so the flux spreads over
:math:`2\\pi r^2` rather than :math:`4\\pi r^2`. The two twos cancel into
:math:`4\\pi F` rather than the :math:`8\\pi F` a full-space singly-charged
source would give.

This is the same screened Green's function as the membrane footprint already in
:mod:`piezo1.physics.membrane` — a 1/r source with an exponential cutoff — and
it needs exactly two numbers this project already produces: the tag distance
from :mod:`piezo1.structure.fusion` and the unitary current from
:mod:`piezo1.physics.permeation`.

**What the model is for.** A JF646-BAPTA HaloTag reports calcium through a
sensor with a sub-micromolar Kd. If the nanodomain at the tag is far above that
Kd, the sensor is saturated whenever its own channel opens, and puncta
brightness then reports *how many tags are labelled and how often the channel
opens* — not local calcium amplitude. That is a falsifiable claim, and
:meth:`Nanodomain.falsifiers` states what would break it.

**The linearised buffer is an approximation** that fails close in, where the
source depletes free buffer and the true concentration exceeds this estimate.
It therefore errs *low* near the channel, which for a saturation argument is the
safe direction.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from .permeation import F_FARADAY

__all__ = ["Nanodomain", "calcium_at", "screening_length", "saturation",
           "distance_for_occupancy", "sweep"]


def screening_length(diffusivity: float | None = None,
                     buffer_kon: float | None = None,
                     buffer_concentration: float | None = None) -> float:
    """λ = sqrt(D / (k_on·[B])), in metres.

    The distance over which buffer absorbs the gradient. Only the *product*
    ``k_on·[B]`` enters, which is why the two are not separately identifiable
    from a nanodomain measurement.
    """
    diffusivity = _P.value("nanodomain.d_calcium") if diffusivity is None else diffusivity
    buffer_kon = _P.value("nanodomain.buffer_kon") if buffer_kon is None else buffer_kon
    buffer_concentration = (_P.value("nanodomain.buffer_concentration")
                            if buffer_concentration is None else buffer_concentration)
    rate = buffer_kon * buffer_concentration
    if rate <= 0:
        return np.inf
    return float(np.sqrt(diffusivity / rate))


def calcium_at(r_m, current_A: float, diffusivity: float | None = None,
               length: float | None = None, resting: float | None = None):
    """Free calcium at distance ``r`` from a source carrying ``current_A``.

    Returns molar concentration, with the resting level added — the nanodomain
    sits on top of the bulk rather than replacing it, which matters only when
    the two are comparable.
    """
    diffusivity = _P.value("nanodomain.d_calcium") if diffusivity is None else diffusivity
    resting = _P.value("nanodomain.resting_calcium") if resting is None else resting
    length = screening_length() if length is None else length

    r = np.asarray(r_m, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        # mol/m^3 from A / (C/mol * m^2/s * m); /1000 to molar.
        molar = (current_A / (4.0 * np.pi * F_FARADAY * diffusivity * r)
                 * np.exp(-r / length)) / 1000.0
    return np.where(r > 0, molar + resting, np.inf)


def saturation(concentration_M, kd: float | None = None):
    """Fractional occupancy of a 1:1 sensor: ``[Ca] / ([Ca] + Kd)``."""
    kd = _P.value("nanodomain.sensor_kd") if kd is None else kd
    c = np.asarray(concentration_M, dtype=float)
    return c / (c + kd)


def distance_for_occupancy(target: float, current_A: float,
                           kd: float | None = None) -> float:
    """How far from the source the sensor drops to ``target`` occupancy.

    Solved by bisection on a monotonically falling profile. The answer is the
    headline of the falsification argument: if it is far larger than any
    plausible tag distance, the prediction cannot be escaped by moving the tag.
    """
    kd = _P.value("nanodomain.sensor_kd") if kd is None else kd
    if not 0.0 < target < 1.0:
        raise ValueError("target occupancy must be strictly between 0 and 1")

    # Resting calcium alone holds the sensor part-occupied, so occupancy has a
    # floor it can never fall below however far away the channel is. Asking for
    # a target under that floor has no answer, and saying so is better than
    # returning the largest distance searched.
    floor = float(saturation(_P.value("nanodomain.resting_calcium"), kd))
    if target <= floor:
        return float("inf")

    lo, hi = 1e-10, 1e-3                       # 0.1 nm to 1 mm
    if float(saturation(calcium_at(lo, current_A), kd)) < target:
        return float("nan")                    # never reached, even at the mouth
    if float(saturation(calcium_at(hi, current_A), kd)) > target:
        return float("inf")                    # still saturated a millimetre out
    for _ in range(200):
        mid = np.sqrt(lo * hi)                 # geometric: the profile is ~1/r
        if float(saturation(calcium_at(mid, current_A), kd)) > target:
            lo = mid
        else:
            hi = mid
    return float(np.sqrt(lo * hi))


@dataclass
class Nanodomain:
    """The gradient one open channel makes, evaluated where the tag sits."""

    current_A: float                    # total unitary current
    calcium_fraction: float
    distance_m: float                   # representative tag distance
    envelope_m: tuple = ()              # (min, max) reachable tag distances
    meta: dict = field(default_factory=dict)

    @property
    def calcium_current_A(self) -> float:
        return self.current_A * self.calcium_fraction

    @property
    def screening_length_m(self) -> float:
        return screening_length()

    @property
    def concentration_M(self) -> float:
        return float(calcium_at(self.distance_m, self.calcium_current_A))

    @property
    def occupancy(self) -> float:
        return float(saturation(self.concentration_M))

    @property
    def resting_occupancy(self) -> float:
        """Occupancy from bulk resting calcium alone, with no channel open.

        Not a detail: at 100 nM resting against a 0.2 uM Kd this is already 33%,
        so the sensor's dynamic range is 33% to ~100% rather than 0 to 100%, and
        "saturated" has to mean well above that floor.
        """
        return float(saturation(_P.value("nanodomain.resting_calcium")))

    @property
    def saturated(self) -> bool:
        """Occupancy above 90%: the sensor reports opening, not amplitude."""
        return self.occupancy > 0.9

    def envelope_range(self) -> tuple:
        """Concentration and occupancy across the whole accessible envelope.

        Round 31 established that the tag samples a region rather than sitting
        at a point, and that the centroid is not the mean distance. Quoting a
        single number here would inherit a precision the placement does not
        have, so the range is reported alongside it.
        """
        if not self.envelope_m:
            return ()
        near, far = min(self.envelope_m), max(self.envelope_m)
        c_far = float(calcium_at(far, self.calcium_current_A))
        c_near = float(calcium_at(near, self.calcium_current_A))
        return (c_far, c_near, float(saturation(c_far)),
                float(saturation(c_near)))

    def falsifiers(self) -> dict:
        """What would have to be true for the saturation claim to fail.

        A prediction that cannot be broken by any parameter in its own model is
        not obviously a prediction, so each of these is a number someone could
        go and measure.
        """
        half = distance_for_occupancy(0.5, self.calcium_current_A)
        # Ca fraction that would leave the sensor half-occupied where it sits.
        kd = _P.value("nanodomain.sensor_kd")
        resting = _P.value("nanodomain.resting_calcium")
        needed = (kd - resting) * 1000.0 * (
            4.0 * np.pi * F_FARADAY * _P.value("nanodomain.d_calcium")
            * self.distance_m) / np.exp(-self.distance_m / self.screening_length_m)
        fraction = needed / self.current_A if self.current_A else np.inf
        # Buffer strength that would pull lambda down to the tag distance.
        rate = _P.value("nanodomain.d_calcium") / self.distance_m ** 2
        return {
            "distance_for_half_occupancy_m": half,
            "calcium_fraction_for_half_occupancy": float(fraction),
            "buffer_rate_for_lambda_equal_distance_per_s": float(rate),
            "buffer_concentration_needed_M": float(
                rate / _P.value("nanodomain.buffer_kon")),
        }

    def summary(self) -> str:
        return (f"{self.concentration_M * 1e6:.1f} uM at "
                f"{self.distance_m * 1e9:.1f} nm "
                f"(lambda {self.screening_length_m * 1e9:.0f} nm); sensor "
                f"{self.occupancy:.1%} occupied — "
                f"{'saturated' if self.saturated else 'graded'}")


def sweep(current_A: float, distances_nm=None, fractions=None,
          buffer_concentrations=None) -> list:
    """Occupancy across tag distance, calcium share and buffering.

    The three quantities the answer could plausibly depend on, varied together
    rather than one at a time, because the claim is about whether *any*
    combination escapes saturation.
    """
    distances_nm = np.asarray(
        [2, 4, 6, 10, 20] if distances_nm is None else distances_nm, float)
    fractions = np.asarray(
        [0.005, 0.02, 0.05, 0.2] if fractions is None else fractions, float)
    buffer_concentrations = np.asarray(
        [1e-5, 1e-4, 1e-3, 1e-2] if buffer_concentrations is None
        else buffer_concentrations, float)

    out = []
    for buffer in buffer_concentrations:
        length = screening_length(buffer_concentration=buffer)
        for fraction in fractions:
            for distance in distances_nm:
                concentration = float(calcium_at(
                    distance * 1e-9, current_A * fraction, length=length))
                out.append({
                    "distance_nm": float(distance),
                    "calcium_fraction": float(fraction),
                    "buffer_M": float(buffer),
                    "lambda_nm": length * 1e9,
                    "calcium_uM": concentration * 1e6,
                    "occupancy": float(saturation(concentration)),
                })
    return out
