"""Ion selectivity: the reversal potential a charged pore has, and what it means.

A conductance says how *much* current flows. Selectivity says *what carries it*,
and it is the one pore property that a fixed charge — and only a fixed charge —
can produce. Between identical baths every channel reverses at zero volts
whatever it is made of, so selectivity has to be measured across a gradient.

**The protocol is the published one.** Coste et al. 2015 measured mPiezo1's
chloride permeability as a dilution potential: 150 mM NaCl inside, 30 mM NaCl
outside, reversal potential read off the I-V and inverted through the
Goldman-Hodgkin-Katz voltage equation, giving **P_Cl/P_Na = 0.14**. This runs
the same protocol on the model — same salt, same two concentrations, same
inversion — so the two numbers are comparable rather than merely both about
selectivity.

**Sign conventions, stated once.** ``V`` here is the membrane potential in the
electrophysiological sense, cytosolic minus extracellular. Which end of the pore
profile is cytosolic is *measured* by :func:`piezo1.physics.pore_charge
.cytosolic_end` rather than assumed, because getting it backwards turns a
cation-selective channel into an anion-selective one and the number still looks
perfectly reasonable.

**Access resistance is deliberately absent.** At the reversal potential no
current flows, so there is no IR drop anywhere in series and the answer does not
depend on the mouths. That is why the root is found on the pore current itself
rather than on the series-corrected one, which is scrubbed of its sign.

**What it is calibrated against.** An uncharged pore is not selectivity-free: it
still develops a liquid-junction potential from the two ions' different
mobilities, and the ratio it returns must be the mobility ratio, near one and on
the anion's side because Cl- is the faster ion. A negatively charged pore must
move it below one, a positively charged pore above. Those three are known
answers, they are checked before any measured charge is believed, and the middle
one is what would fail if a sign were wrong anywhere in the chain.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from ._pnp_kernels import F_FARADAY, R_GAS
from .permeation import IonSpecies, solve_pnp, sodium_species

__all__ = ["SelectivityResult", "ghk_permeability_ratio", "reversal_potential",
           "measure_selectivity", "REVERSAL_BRACKET_V"]

#: Voltage range searched for the reversal potential. +/-200 mV comfortably
#: contains the Nernst limit of the protocol (-43 mV for a perfectly
#: cation-selective pore at 5:1), and a root outside it would mean the model
#: had left physiology rather than that the bracket was too tight.
REVERSAL_BRACKET_V = 0.2


@dataclass
class SelectivityResult:
    """A permeability ratio, the potential it came from, and its caveats."""

    reversal_V: float
    p_anion_over_cation: float
    cation: str
    anion: str
    high_M: float
    low_M: float
    cytosolic_index: int
    published: float | None = None
    converged: bool = True
    meta: dict = field(default_factory=dict)

    @property
    def reversal_mV(self) -> float:
        return self.reversal_V * 1e3

    @property
    def cation_selective(self) -> bool:
        """Whether the model prefers cations at all — the qualitative claim."""
        return self.p_anion_over_cation < 1.0

    @property
    def ratio_to_published(self) -> float | None:
        if not self.published:
            return None
        return self.p_anion_over_cation / self.published

    def summary(self) -> str:
        published = ("" if self.published is None
                     else f" against a published {self.published:.2f}")
        return (f"P_{self.anion}/P_{self.cation} = "
                f"{self.p_anion_over_cation:.3f}{published}; reversal "
                f"{self.reversal_mV:+.1f} mV at "
                f"{self.high_M * 1e3:.0f}/{self.low_M * 1e3:.0f} mM")


def ghk_permeability_ratio(reversal_V: float, inside_M: float, outside_M: float,
                           temperature: float | None = None) -> float:
    """Invert the GHK voltage equation for ``P_cation / P_anion``.

    For a single 1:1 salt at two concentrations,

    .. math:: V = \\frac{RT}{F}\\ln\\frac{P_+ c_o + P_- c_i}
                                          {P_+ c_i + P_- c_o}

    which rearranges to a closed form. Two limits fix the conventions and are
    worth stating because they are the ones the tests check: a perfectly
    cation-selective pore gives the cation's Nernst potential
    ``(RT/F) ln(c_o/c_i)``, and an entirely unselective one gives zero.
    """
    temperature = temperature or _P.value("permeation.temperature")
    thermal = R_GAS * temperature / F_FARADAY
    x = float(np.exp(reversal_V / thermal))
    numerator = inside_M - x * outside_M
    denominator = x * inside_M - outside_M
    if denominator == 0.0:
        return np.inf
    return numerator / denominator


def reversal_potential(profile, fixed_charge=None,
                       species: list[IonSpecies] | None = None,
                       cytosolic_index: int = 0, wetting=None,
                       tol: float = 1e-9) -> tuple[float, bool]:
    """Membrane potential at which the net pore current is zero, in volts.

    Returns ``(V, converged)``. The solver's applied potential is zero at the
    first slice and ``voltage`` at the last, so the membrane potential —
    cytosolic minus extracellular — is ``-voltage`` when the first slice is the
    cytosolic one and ``+voltage`` when it is not. That single sign is the whole
    orientation dependence, and it is why ``cytosolic_index`` is required rather
    than defaulted silently.

    Bisection rather than Newton: the current is monotone in voltage for a
    passive pore, each evaluation is a full Gummel solve, and a bracketed method
    cannot be thrown by the discontinuity at a blocked slice.
    """
    species = species or sodium_species()
    sign = -1.0 if cytosolic_index == 0 else 1.0

    def current(membrane_V: float) -> float:
        result = solve_pnp(profile, wetting, voltage=sign * membrane_V,
                           species=species, fixed_charge=fixed_charge)
        if result.blocked_by:
            raise ValueError(f"no current to reverse: {result.blocked_by}")
        return result.pore_current

    lo, hi = -REVERSAL_BRACKET_V, REVERSAL_BRACKET_V
    f_lo, f_hi = current(lo), current(hi)
    if f_lo == 0.0:
        return lo, True
    if f_hi == 0.0:
        return hi, True
    if np.sign(f_lo) == np.sign(f_hi):
        # Not a failure of the search: a current of one sign across +/-200 mV
        # means the model has no reversal potential in any physiological range.
        return float("nan"), False

    for _ in range(80):
        mid = 0.5 * (lo + hi)
        f_mid = current(mid)
        if np.sign(f_mid) == np.sign(f_lo):
            lo, f_lo = mid, f_mid
        else:
            hi = mid
        if hi - lo < tol:
            break
    return 0.5 * (lo + hi), True


def measure_selectivity(profile, fixed_charge=None, cytosolic_index: int = 0,
                        wetting=None, high: float | None = None,
                        low: float | None = None,
                        species: list[IonSpecies] | None = None
                        ) -> SelectivityResult:
    """Run the published dilution protocol on this pore and invert it.

    The concentrated bath goes on the cytosolic side and the dilute one on the
    extracellular side, which is how Coste et al. arranged it. Everything the
    answer depends on — which ions, which concentrations, which end is which —
    is carried on the result, because the ratio alone is not interpretable
    without them.
    """
    high = _P.value("permeation.dilution_high") if high is None else high
    low = _P.value("permeation.dilution_low") if low is None else low
    if species is None:
        left, right = ((high, low) if cytosolic_index == 0 else (low, high))
        species = sodium_species(left, right)

    cation = next(s for s in species if s.valence > 0)
    anion = next(s for s in species if s.valence < 0)

    reversal, converged = reversal_potential(
        profile, fixed_charge=fixed_charge, species=species,
        cytosolic_index=cytosolic_index, wetting=wetting)

    if not converged or not np.isfinite(reversal):
        ratio = float("nan")
    else:
        forward = ghk_permeability_ratio(reversal, high, low)
        ratio = float("inf") if forward == 0.0 else 1.0 / forward

    return SelectivityResult(
        reversal_V=float(reversal), p_anion_over_cation=float(ratio),
        cation=cation.name, anion=anion.name, high_M=high, low_M=low,
        cytosolic_index=cytosolic_index,
        published=_P.value("permeation.published_pcl_pna"),
        converged=converged,
        meta={"protocol": f"{high * 1e3:.0f} mM {cation.name}{anion.name} "
                          f"cytosolic / {low * 1e3:.0f} mM extracellular, "
                          f"reversal potential inverted through GHK",
              "temperature_K": _P.value("permeation.temperature"),
              "note": "GHK assumes a constant field and independent ions, "
                      "which the model itself does not — the inversion is used "
                      "because it is what the published number was derived "
                      "with, not because it is the better theory"})
