"""Ion permeation: 1-D Poisson-Nernst-Planck over the measured pore profile.

The pore radius profile from :mod:`piezo1.structure.pore` says how wide the
conduction pathway is at every height. This turns that into a current.

**The model.** Steady-state drift-diffusion for each ion species down a channel
of varying cross-section, with the electrostatic potential solved
self-consistently:

.. math::

    J_i &= -D_i A_i(z)\\left[\\frac{dc_i}{dz}
            + \\frac{z_i F}{RT} c_i \\frac{d\\phi}{dz}\\right], \\quad
            \\frac{dJ_i}{dz} = 0 \\\\
    \\frac{d}{dz}\\left(\\varepsilon A \\frac{d\\phi}{dz}\\right)
        &= -A\\left(F\\sum_i z_i c_i + \\rho_\\text{fixed}\\right)

Discretised with the Scharfetter-Gummel flux, which stays stable when drift
dominates diffusion — a centred difference oscillates there and produces
negative concentrations.

**The Poisson half does not converge here, and that is a result rather than a
bug.** Feeding the solved potential back through Poisson diverges: measured on
11ZC the potential went -0.37 V, then -171 V, then -2e16 V. Adding the screening
derivative that a proper Newton step needs (:func:`_poisson_newton_step`) makes
the operator negative-definite and still does not converge — the update plateaus
and the potential swings +/-1.5 V. The reason is physical. In 150 mM the Debye
length is **5.7-8.1 A** while PIEZO1's open bottleneck radius is **3.3 A**, so
the double layers from opposite walls overlap completely and the pore has no
electroneutral core for a Gummel map to relax onto.

So the potential is solved in the **electroneutral limit** instead
(:func:`_ohmic_potential`): current continuity,
:math:`\\nabla\\cdot(\\sigma A\\nabla\\phi) = 0`, with the local conductivity set
by the local concentrations. That converges, and it agrees with the independent
closed-form :func:`series_conductance` to 1.5% (41.0 vs 40.4 pS on 11ZC). The
Poisson machinery is kept, and :func:`debye_length` is reported on every result,
because the honest statement is that a continuum treatment of a 3-angstrom pore
is at the edge of its validity — not that it was never attempted.

**The pore is not the only resistor.** A short wide pore is limited as much by
the spreading resistance of its own mouths, so Hall's access resistance is added
in series at each end.

**Gating comes from Round 19, not from here.** A pore wide enough to pass an ion
still carries no current if it has dewetted, because the ion's hydration shell is
not there to be shed into. :func:`blocking_mechanisms` returns *every* reason
rather than the first, because 8YEZ is shut both sterically and by dewetting
while 7WLU is shut only sterically, and collapsing that would hide the
comparison this round exists to make.

**What it gives, against what was measured.** 41.0 pS for the open 11ZC against
a published 25-30 pS — high by about half. Sweeping the two *unmeasured*
confinement parameters over plausible ranges (in-pore diffusivity 0.25-1.0 of
bulk, ion radius 1.0-2.0 A) moves it across **16-94 pS**, which straddles the
measurement. The model can therefore be made to agree, but only by choosing
values for two things nobody has measured, so agreement would be tuning rather
than prediction. All-atom MD remains a non-goal.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P
from ._pnp_kernels import _bernoulli, _solve_with_dirichlet

__all__ = ["IonSpecies", "PermeationResult", "solve_pnp", "conductance",
           "series_conductance", "access_resistance", "default_species",
           "blocking_mechanisms", "debye_length", "F_FARADAY", "R_GAS"]

#: Faraday constant, C/mol. Definitional since the 2019 SI redefinition.
F_FARADAY = 96485.33212
#: Molar gas constant, J/(mol K). Definitional since 2019.
R_GAS = 8.314462618


@dataclass(frozen=True)
class IonSpecies:
    """One permeant species: how it drifts, diffuses and how big it is."""

    name: str
    valence: int
    diffusivity: float                 # m^2/s, bulk
    radius: float                      # A, crystal
    concentration: float               # M, in both baths


def default_species(calcium: float = 0.0) -> list[IonSpecies]:
    """A symmetric monovalent bath, optionally with calcium added.

    Monovalent by default because the 25-30 pS the model is checked against was
    measured that way. Calcium is added as an extra species rather than
    replacing the salt, which is what a real recording does.
    """
    salt = _P.value("permeation.bath_concentration")
    scale = _P.value("permeation.diffusion_scale")
    species = [
        IonSpecies("cation", +1, _P.value("permeation.diffusion_cation") * scale,
                   _P.value("permeation.radius_cation"), salt),
        IonSpecies("anion", -1, _P.value("permeation.diffusion_anion") * scale,
                   _P.value("permeation.radius_anion"), salt),
    ]
    if calcium > 0:
        species.append(IonSpecies(
            "calcium", +2, _P.value("permeation.diffusion_calcium") * scale,
            _P.value("permeation.radius_calcium"), calcium))
        # Balance the extra charge with anion, so the bath is electroneutral.
        species[1] = IonSpecies(species[1].name, -1, species[1].diffusivity,
                                species[1].radius, salt + 2 * calcium)
    return species


@dataclass
class PermeationResult:
    """A current, and everything needed to disbelieve it."""

    current: float                     # A, at `voltage`
    conductance: float                 # S
    voltage: float
    z: np.ndarray                      # m, along the axis
    radius: np.ndarray                 # m
    potential: np.ndarray              # V
    concentrations: dict = field(default_factory=dict)
    fluxes: dict = field(default_factory=dict)
    blocked_by: str | None = None
    access_ohm: float = 0.0
    pore_ohm: float = 0.0
    converged: bool = True
    iterations: int = 0
    meta: dict = field(default_factory=dict)

    @property
    def conductance_pS(self) -> float:
        return self.conductance * 1e12

    @property
    def is_conducting(self) -> bool:
        return self.blocked_by is None and self.conductance > 0.0

    def summary(self) -> str:
        if self.blocked_by:
            return f"no current: {self.blocked_by}"
        return (f"{self.conductance_pS:.1f} pS "
                f"({self.current * 1e12:.1f} pA at {self.voltage * 1e3:.0f} mV); "
                f"pore {self.pore_ohm / 1e9:.1f} GOhm + access "
                f"{self.access_ohm / 1e9:.1f} GOhm")


def access_resistance(mouth_radius: float, conductivity: float) -> float:
    """Hall's spreading resistance for one circular pore mouth, in ohm.

    ``R = 1 / (4 * sigma * a)``. Applied at both ends. For a pore only a few
    nanometres long this is not a correction — it can be most of the resistance.
    """
    if mouth_radius <= 0 or conductivity <= 0:
        return np.inf
    return 1.0 / (4.0 * conductivity * mouth_radius)


def _bulk_conductivity(species: list[IonSpecies], temperature: float) -> float:
    """Ohmic conductivity of the bath, S/m, from the Nernst-Einstein relation."""
    total = 0.0
    for s in species:
        # mol/m^3 = 1000 * M
        total += (s.valence ** 2 * F_FARADAY ** 2 * s.diffusivity
                  * s.concentration * 1000.0) / (R_GAS * temperature)
    return total


def _accessible_area(radius_m: np.ndarray, ion_radius_A: float) -> np.ndarray:
    """Cross-section an ion of the given radius can occupy, m^2.

    A hard-sphere exclusion: the ion's centre cannot come closer to the wall
    than its own radius. Crude at atomic scale — it is the same approximation
    that makes a continuum model of a 3 angstrom pore arguable at all — and it
    is one of the two knobs the answer is most sensitive to.
    """
    usable = np.maximum(radius_m - ion_radius_A * 1e-10, 0.0)
    return np.pi * usable ** 2


def series_conductance(z: np.ndarray, radius: np.ndarray,
                       species: list[IonSpecies] | None = None,
                       temperature: float | None = None) -> dict:
    """Conductance from the series-resistor formula, independent of the solver.

    .. math:: R = \\int \\frac{dz}{\\sigma A(z)} + 2 R_\\text{access}

    This is what PNP must reduce to at low voltage with no fixed charge, and it
    is derived here without touching the Scharfetter-Gummel discretisation, the
    Gummel loop, or the Poisson equation. If the two disagree, one of them is
    wrong — which is the point of having it.
    """
    species = species or default_species()
    temperature = temperature or _P.value("permeation.temperature")
    sigma = _bulk_conductivity(species, temperature)

    # Each species carries current through its own accessible area, so the local
    # conductance per unit length is the sum over species.
    per_length = np.zeros_like(z)
    for s in species:
        area = _accessible_area(radius, s.radius)
        per_length += ((s.valence ** 2 * F_FARADAY ** 2 * s.diffusivity
                        * s.concentration * 1000.0) / (R_GAS * temperature)) * area

    if np.any(per_length <= 0):
        return {"conductance": 0.0, "pore_ohm": np.inf, "access_ohm": 0.0,
                "blocked": "a slice is too narrow for any ion to enter"}

    pore_ohm = float(np.trapezoid(1.0 / per_length, z))
    mouth = float(max(radius[0], radius[-1]))
    access_ohm = 2.0 * access_resistance(mouth, sigma)
    total = pore_ohm + access_ohm
    return {"conductance": 1.0 / total, "pore_ohm": pore_ohm,
            "access_ohm": access_ohm, "conductivity": sigma}


def _solve_nernst_planck(z, area, potential, species, temperature, bath):
    """Steady concentration profile and the constant flux it carries.

    Scharfetter-Gummel: the flux between two nodes is exact for a linear
    potential and constant coefficients across the interval, so the scheme does
    not go unstable when drift dominates. A centred difference does, and shows
    it as negative concentrations rather than as an error.
    """
    n = len(z)
    h = np.diff(z)
    thermal = R_GAS * temperature / F_FARADAY
    delta = species.valence * np.diff(potential) / thermal      # (n-1,)

    face_area = 0.5 * (area[:-1] + area[1:])
    a = species.diffusivity * face_area / h                     # (n-1,)

    lower = np.zeros(n)
    diag = np.zeros(n)
    upper = np.zeros(n)
    rhs = np.zeros(n)

    b_plus = _bernoulli(delta)          # multiplies c_{k+1}
    b_minus = _bernoulli(-delta)        # multiplies c_k

    for k in range(1, n - 1):
        left, right = a[k - 1], a[k]
        lower[k] = left * b_minus[k - 1]
        diag[k] = -(left * b_plus[k - 1] + right * b_minus[k])
        upper[k] = right * b_plus[k]

    diag[0] = diag[-1] = 1.0
    rhs[0] = rhs[-1] = bath

    matrix = np.zeros((n, n))
    matrix[np.arange(n), np.arange(n)] = diag
    matrix[np.arange(1, n), np.arange(n - 1)] = lower[1:]
    matrix[np.arange(n - 1), np.arange(1, n)] = upper[:-1]

    concentration = _solve_with_dirichlet(matrix, rhs)
    flux = a[0] * (concentration[0] * b_minus[0] - concentration[1] * b_plus[0])
    return concentration, float(flux)


def debye_length(species: list[IonSpecies], temperature: float,
                 permittivity_relative: float) -> float:
    """Screening length of the bath, in metres.

    Reported with every result because it decides whether a continuum treatment
    of this pore means anything. In 150 mM it is 5.7-8.1 A depending on the
    permittivity assumed, against a PIEZO1 bottleneck radius of 3.3 A — so the
    double layers from opposite walls overlap completely and the pore has no
    electroneutral core.
    """
    ionic_strength = 0.5 * sum(s.valence ** 2 * s.concentration * 1000.0
                               for s in species)
    if ionic_strength <= 0:
        return np.inf
    return float(np.sqrt(permittivity_relative * 8.8541878128e-12 * R_GAS
                         * temperature
                         / (2.0 * F_FARADAY ** 2 * ionic_strength)))


def _ohmic_potential(z, areas, concentrations, species, temperature,
                     v_left, v_right):
    """Potential from current continuity — the electroneutral limit of PNP.

    Where the double layers do not overlap, the pore interior is neutral and
    Poisson reduces to :math:`\\nabla\\cdot(\\sigma A \\nabla\\phi) = 0`, with
    the local conductivity set by the local concentrations. That is what is
    solved here, and it is well-posed: a Laplace problem with positive
    coefficients, so the Gummel loop around it converges.

    The full Poisson coupling is implemented in :func:`_poisson_newton_step` but
    is **not** used by default, and the reason is physics rather than
    convenience — see :func:`debye_length` and the module docstring.
    """
    n = len(z)
    h = np.diff(z)
    conductance_face = np.zeros(n - 1)
    for s in species:
        c = concentrations[s.name]
        local = ((s.valence ** 2 * F_FARADAY ** 2 * s.diffusivity)
                 / (R_GAS * temperature)) * c * areas[s.name]
        conductance_face += 0.5 * (local[:-1] + local[1:]) / h

    matrix = np.zeros((n, n))
    rhs = np.zeros(n)
    for k in range(1, n - 1):
        left, right = conductance_face[k - 1], conductance_face[k]
        matrix[k, k - 1] = left
        matrix[k, k] = -(left + right)
        matrix[k, k + 1] = right
    matrix[0, 0] = matrix[-1, -1] = 1.0
    rhs[0], rhs[-1] = v_left, v_right
    return _solve_with_dirichlet(matrix, rhs)


def solve_pnp(profile, wetting=None, voltage: float | None = None,
              species: list[IonSpecies] | None = None,
              fixed_charge: np.ndarray | None = None,
              max_iterations: int = 200, tol: float = 1e-10,
              relaxation: float = 0.4) -> PermeationResult:
    """Solve 1-D PNP over a measured pore profile.

    ``wetting`` is a :class:`~piezo1.analysis.hydration.WettingPrediction`. When
    it says the pore is shut, the result is zero current with the reason
    recorded — the continuum equations would happily push ions through a
    dewetted pore, because nothing in them knows about the hydration shell.
    """
    voltage = _P.value("permeation.test_voltage") if voltage is None else voltage
    species = species or default_species()
    temperature = _P.value("permeation.temperature")
    permittivity = _P.value("permeation.permittivity_pore") * 8.8541878128e-12

    z = np.asarray(profile.z, dtype=float) * 1e-10          # A -> m
    radius = np.asarray(profile.radius, dtype=float) * 1e-10
    order = np.argsort(z)
    z, radius = z[order], radius[order]

    reasons = blocking_mechanisms(wetting, radius, species)
    if reasons:
        return PermeationResult(
            current=0.0, conductance=0.0, voltage=voltage, z=z, radius=radius,
            potential=np.linspace(0.0, voltage, len(z)),
            blocked_by=" AND ".join(reasons),
            meta={"n_slices": len(z), "mechanisms": reasons,
                  "n_mechanisms": len(reasons),
                  "min_radius_A": float(radius.min()) * 1e10})

    geometric_area = np.pi * radius ** 2
    areas = {s.name: _accessible_area(radius, s.radius) for s in species}
    charge_fixed = (np.zeros_like(z) if fixed_charge is None
                    else np.asarray(fixed_charge, dtype=float))

    potential = np.linspace(0.0, voltage, len(z))
    concentrations, fluxes = {}, {}
    converged, used = False, 0

    for used in range(1, max_iterations + 1):
        for s in species:
            concentration, flux = _solve_nernst_planck(
                z, areas[s.name], potential, s, temperature,
                s.concentration * 1000.0)
            concentrations[s.name] = concentration
            fluxes[s.name] = flux

        updated = _ohmic_potential(z, areas, concentrations, species,
                                   temperature, 0.0, voltage)
        change = float(np.max(np.abs(updated - potential)))
        potential = (1.0 - relaxation) * potential + relaxation * updated
        if change < tol * max(abs(voltage), 1e-12):
            converged = True
            break

    current = sum(s.valence * F_FARADAY * fluxes[s.name] for s in species)

    # Access resistance sits outside the solved region, so it is applied to the
    # solved pore current as a series correction rather than as a boundary term.
    sigma = _bulk_conductivity(species, temperature)
    lam = debye_length(species, temperature,
                       _P.value("permeation.permittivity_pore"))
    access_ohm = 2.0 * access_resistance(float(max(radius[0], radius[-1])), sigma)
    pore_ohm = abs(voltage / current) if current != 0 else np.inf
    total_ohm = pore_ohm + access_ohm
    conductance_value = 1.0 / total_ohm if np.isfinite(total_ohm) else 0.0

    return PermeationResult(
        current=voltage / total_ohm if np.isfinite(total_ohm) else 0.0,
        conductance=conductance_value, voltage=voltage, z=z, radius=radius,
        potential=potential, concentrations=concentrations, fluxes=fluxes,
        access_ohm=access_ohm, pore_ohm=pore_ohm, converged=converged,
        iterations=used,
        meta={"n_slices": len(z), "conductivity_S_per_m": sigma,
              "species": [s.name for s in species],
              "diffusion_scale": _P.value("permeation.diffusion_scale"),
              "debye_length_A": lam * 1e10,
              "min_radius_A": float(radius.min()) * 1e10,
              "double_layers_overlap": bool(lam > float(radius.min())),
              "note": "continuum model of an atomic-scale pore; the two "
                      "confinement parameters are unmeasured"})


def blocking_mechanisms(wetting, radius: np.ndarray,
                        species: list[IonSpecies]) -> list[str]:
    """**Every** reason no current flows, not just the first one found.

    Round 19 established that "does an ion fit?" and "would water stay?" are
    different questions, and the structures answer them differently: 8YEZ is
    shut *both* sterically and by dewetting, 7WLU only sterically. Returning on
    the first match would have collapsed that distinction and made the two look
    identical — which is precisely the comparison this round exists to report.
    """
    reasons = []
    smallest = min(s.radius for s in species if s.valence > 0)
    if float(radius.min()) <= smallest * 1e-10:
        reasons.append(
            f"sterically occluded: the narrowest slice is "
            f"{radius.min() * 1e10:.2f} A, below the {smallest:.2f} A radius "
            f"of the smallest permeant cation")
    if wetting is not None and getattr(wetting, "available", False):
        if wetting.hydrophobic_gate:
            reasons.append(
                f"hydrophobic gate: the lining would dewet "
                f"(wetting score {wetting.score:.2f} > cutoff)")
        if wetting.sterically_occluded and not reasons:
            reasons.append("sterically occluded (wetting analysis)")
    return reasons


def conductance(profile, wetting=None, **kwargs) -> float:
    """Unitary conductance in picosiemens — the number to compare with 25-30 pS."""
    return solve_pnp(profile, wetting, **kwargs).conductance_pS
