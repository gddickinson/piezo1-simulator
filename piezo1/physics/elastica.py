"""Nonlinear membrane mechanics: the axisymmetric Euler–elastica footprint.

:mod:`piezo1.physics.membrane` solves the *linearised* Helfrich problem in the
Monge gauge, which assumes |∇h| ≪ 1. PIEZO1 violates that assumption badly: the
measured dome meets the bilayer at a contact slope near 2.0, about 63°, so the
terms the expansion drops are larger than the ones it keeps. This module drops
the expansion instead.

**Formulation.** Parametrise the meridian by arc length ``s``, with ``r(s)`` the
distance from the symmetry axis, ``z(s)`` the height and ``ψ(s)`` the angle of
the tangent to the horizontal. The two principal curvatures are then exact:

.. math::   c_1 = \\dot\\psi, \\qquad c_2 = \\frac{\\sin\\psi}{r}

and the Helfrich energy, with no expansion anywhere, is

.. math::   E = \\int \\left[\\frac{\\kappa}{2}(c_1+c_2)^2 + \\gamma\\right]
                 2\\pi r\\, ds \\; - \\; \\gamma A_{\\rm proj}

Writing ``M = c₁ + c₂`` and imposing ``ṙ = cos ψ`` with a Lagrange multiplier
``η``, the Euler–Lagrange equations reduce to a first-order system:

.. math::

    \\dot r = \\cos\\psi, \\quad
    \\dot z = \\sin\\psi, \\quad
    \\dot\\psi = M - \\frac{\\sin\\psi}{r}, \\\\
    \\dot M = \\frac{\\eta \\sin\\psi}{\\kappa r}, \\quad
    \\dot\\eta = \\frac{\\kappa}{2}M^2 + \\gamma
                - \\frac{\\kappa M \\sin\\psi}{r}

**The conserved quantity.** The Lagrangian has no explicit ``s`` dependence, so
its Hamiltonian is a constant of the motion:

.. math::   H = \\frac{\\kappa}{2}M^2 r - \\kappa M \\sin\\psi
                + \\eta\\cos\\psi - \\gamma r

``H`` is (up to 2π) the **axial force** transmitted through the membrane at
radius ``r``. For an inclusion nobody is pulling on, ``H = 0``. That is imposed
as a boundary condition, and how far ``H`` drifts from zero along the solved
profile is then a free accuracy diagnostic — see
:attr:`ElasticaSolution.force_residual`, which runs at ~1e-11 in practice.

**Units** follow :mod:`~piezo1.physics.membrane`: nm, k_BT, k_BT/nm².
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_bvp
from scipy.special import k0, k1

from .membrane import (MembraneParameters, analytic_energy, decay_length,
                       solve_footprint)
from ..parameters import PARAMETERS as _P

__all__ = ["ElasticaSolution", "solve_elastica", "LinearComparison",
           "compare_with_linear", "shape_equations", "axial_force",
           "PIEZO1_CONTACT_SLOPE"]

#: Contact slope of the measured 7WLT dome where it meets the bilayer, from
#: ``DomeModel.footprint``'s spherical-cap geometry. 63 degrees. Recorded here
#: because it is the number that makes the linear theory inapplicable.
PIEZO1_CONTACT_SLOPE = _P.value("membrane.piezo1_contact_slope")


def shape_equations(s, y, kappa: float, tension: float) -> np.ndarray:
    """Right-hand side of the axisymmetric Helfrich shape equations.

    ``y`` is ``[r, z, psi, M, eta]``; ``s`` is unused (the system is autonomous)
    but kept in the signature for :func:`scipy.integrate.solve_bvp`.
    """
    r, _z, psi, M, eta = y
    r_safe = np.maximum(r, 1e-9)
    sin_p, cos_p = np.sin(psi), np.cos(psi)
    return np.vstack([
        cos_p,
        sin_p,
        M - sin_p / r_safe,
        eta * sin_p / (kappa * r_safe),
        0.5 * kappa * M ** 2 + tension - kappa * M * sin_p / r_safe,
    ])


def axial_force(y: np.ndarray, kappa: float, tension: float) -> np.ndarray:
    """The conserved Hamiltonian H, i.e. the axial force per 2π.

    Zero for a free inclusion. Deviation from zero along a computed profile
    measures integration error, not physics.
    """
    r, _z, psi, M, eta = y
    return (0.5 * kappa * M ** 2 * r - kappa * M * np.sin(psi)
            + eta * np.cos(psi) - tension * r)


@dataclass
class ElasticaSolution:
    """A fully nonlinear axisymmetric footprint."""

    s: np.ndarray                  # arc length along the meridian, nm
    r: np.ndarray                  # nm
    z: np.ndarray                  # nm
    psi: np.ndarray                # rad, negative where the membrane falls away
    curvature: np.ndarray          # M = c1 + c2, 1/nm
    eta: np.ndarray                # Lagrange multiplier, k_BT/nm
    params: MembraneParameters
    r0: float = 0.0
    slope: float = 0.0
    converged: bool = True
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------- energetics

    @property
    def bending_energy(self) -> float:
        """∫ (κ/2)(c₁+c₂)² dA in k_BT, with the exact area element."""
        return float(np.trapezoid(
            0.5 * self.params.kappa * self.curvature ** 2
            * 2.0 * np.pi * self.r, self.s))

    @property
    def area(self) -> float:
        """True area of the deformed annulus, nm²."""
        return float(np.trapezoid(2.0 * np.pi * self.r, self.s))

    @property
    def projected_area(self) -> float:
        return float(np.pi * (self.r[-1] ** 2 - self.r0 ** 2))

    @property
    def excess_area(self) -> float:
        """Area stored in the deformation, nm².

        Exact here: true area minus projected area. The linear theory's
        ``∫½|∇h|²dA`` is the leading term of this, and at PIEZO1's slope it is
        not a good approximation to it.
        """
        return self.area - self.projected_area

    @property
    def energy(self) -> float:
        """Total excess Helfrich energy, k_BT."""
        return self.bending_energy + self.params.tension * self.excess_area

    # ------------------------------------------------------------ diagnostics

    @property
    def force_residual(self) -> float:
        """max |H| along the profile. Should be ~0; it is an error estimate."""
        y = np.vstack([self.r, self.z, self.psi, self.curvature, self.eta])
        return float(np.abs(axial_force(y, self.params.kappa,
                                        self.params.tension)).max())

    @property
    def max_slope(self) -> float:
        return float(np.abs(np.tan(self.psi)).max())

    @property
    def contact_angle_deg(self) -> float:
        return float(np.degrees(abs(self.psi[0])))

    @property
    def contact_height(self) -> float:
        return float(self.z[0])

    @property
    def outer_radius(self) -> float:
        return float(self.r[-1])

    def height_profile(self, n: int = 400) -> tuple[np.ndarray, np.ndarray]:
        """Resample as h(r) on a uniform radial grid, for plotting or rendering.

        Valid only because ψ never reaches ±90° for the slopes of interest; a
        genuinely overhanging profile has no single-valued h(r) and this would
        silently fold it.
        """
        if np.abs(self.psi).max() >= np.pi / 2:
            raise ValueError("profile overhangs; h(r) is not single valued")
        r = np.linspace(self.r0, self.r[-1], n)
        return r, np.interp(r, self.r, self.z)

    def surface(self, n_angular: int = 96):
        """Revolve into a Cartesian surface, matching FootprintSolution."""
        r, h = self.height_profile()
        theta = np.linspace(0.0, 2 * np.pi, n_angular, endpoint=False)
        R, T = np.meshgrid(r, theta, indexing="ij")
        Z = np.repeat(h[:, None], n_angular, axis=1)
        return R * np.cos(T), R * np.sin(T), Z


def _initial_guess(s: np.ndarray, r0: float, slope: float, kappa: float,
                   tension: float) -> np.ndarray:
    """Seed the BVP with the linear solution.

    ``h = A K₀(r/λ)`` satisfies ∇²h = h/λ², so the curvature guess is free.
    ``η ≈ γr`` is the far-field form implied by H = 0 with ψ → 0.
    """
    lam = decay_length(kappa, tension)
    amp = slope * lam / k1(r0 / lam)
    r = r0 + s
    h = amp * k0(r / lam)
    return np.vstack([r, h, np.arctan(-amp * k1(r / lam) / lam),
                      h / lam ** 2, tension * r])


def solve_elastica(r0: float, slope: float,
                   params: MembraneParameters | None = None,
                   n: int = 600, arc_span: float | None = None,
                   tol: float = 1e-8, max_nodes: int = 200_000,
                   continuation: bool = True) -> ElasticaSolution:
    """Solve the nonlinear footprint for a fixed contact slope at ``r0``.

    Parameters match :func:`~piezo1.physics.membrane.solve_footprint` so the
    two can be compared directly: ``slope`` is ``-dh/dr`` at the inclusion,
    positive meaning the membrane falls away from it.

    ``continuation`` walks the slope up in steps, each solution seeding the
    next. The linear guess alone is usually enough — it reaches slope 2.0
    unaided — but it is not enough at every ``r₀``/λ, and a BVP that fails to
    converge from a bad guess returns a *plausible* wrong shape rather than an
    error, so the safety net is on by default.
    """
    p = params or MembraneParameters()
    kappa, tension = p.kappa, p.tension
    lam = p.decay_length
    if not np.isfinite(lam):
        raise ValueError("zero tension gives an infinite footprint")
    if slope < 0:
        raise ValueError("slope is a magnitude; pass a positive value")
    if arc_span is None:
        arc_span = 12.0 * lam

    s = np.linspace(0.0, arc_span, n)
    guess = _initial_guess(s, r0, min(slope, 0.1), kappa, tension)

    ladder = [slope]
    if continuation and slope > 0.25:
        ladder = list(np.linspace(0.1, slope, max(4, int(slope / 0.4) + 3)))

    def bc_for(target: float):
        psi0 = -np.arctan(target)

        def bc(ya, yb):
            return np.array([ya[0] - r0,
                             ya[2] - psi0,
                             float(axial_force(ya, kappa, tension)),
                             yb[2],
                             yb[1]])
        return bc

    sol = None
    for target in ladder:
        sol = solve_bvp(lambda ss, yy: shape_equations(ss, yy, kappa, tension),
                        bc_for(target), s, guess, tol=tol, max_nodes=max_nodes)
        s, guess = sol.x, sol.y

    out = ElasticaSolution(
        s=sol.x, r=sol.y[0], z=sol.y[1], psi=sol.y[2], curvature=sol.y[3],
        eta=sol.y[4], params=p, r0=r0, slope=slope,
        converged=bool(sol.status == 0),
        meta={"status": int(sol.status), "message": sol.message,
              "n_nodes": int(sol.x.size), "lambda": lam,
              "arc_span": float(arc_span), "tol": tol,
              "continuation_steps": len(ladder),
              "formulation": "arc-length Euler-elastica, H=0"})
    out.meta["force_residual"] = out.force_residual
    return out


# --------------------------------------------------------------------------
# Linear vs nonlinear
# --------------------------------------------------------------------------

@dataclass
class LinearComparison:
    """Side-by-side of the linearised and nonlinear footprint at one slope."""

    r0: float
    slope: float
    params: MembraneParameters
    linear_energy: float
    nonlinear_energy: float
    linear_excess_area: float
    nonlinear_excess_area: float
    contact_angle_deg: float
    force_residual: float

    @property
    def energy_ratio(self) -> float:
        """Nonlinear / linear. Below 1 means the linear theory overestimates."""
        return self.nonlinear_energy / self.linear_energy

    @property
    def area_ratio(self) -> float:
        return self.nonlinear_excess_area / self.linear_excess_area

    @property
    def linear_error(self) -> float:
        """Relative error of the linear energy, |E_lin − E_nl| / E_nl."""
        return abs(self.linear_energy - self.nonlinear_energy) / self.nonlinear_energy

    def summary(self) -> str:
        return (f"contact slope {self.slope:.2f} ({self.contact_angle_deg:.0f} deg): "
                f"linear {self.linear_energy:.2f} kT vs nonlinear "
                f"{self.nonlinear_energy:.2f} kT "
                f"(linear is {self.linear_energy / self.nonlinear_energy:.2f}x); "
                f"excess area {self.linear_excess_area:.0f} vs "
                f"{self.nonlinear_excess_area:.0f} nm^2")


def compare_with_linear(r0: float, slope: float,
                        params: MembraneParameters | None = None,
                        n: int = 600) -> LinearComparison:
    """Solve both theories on matched geometry and report the difference.

    The linear excess area is integrated out to the *same* outer radius the
    nonlinear solve reached, since the two use different independent variables
    and an unmatched cut-off would confound the comparison with a truncation.
    """
    p = params or MembraneParameters()
    nl = solve_elastica(r0, slope, p, n=n)
    lin = solve_footprint(r0, slope, p, r_max=nl.outer_radius, n=max(n, 800))
    return LinearComparison(
        r0=r0, slope=slope, params=p,
        linear_energy=analytic_energy(r0, slope, p.kappa, p.tension),
        nonlinear_energy=nl.energy,
        linear_excess_area=lin.excess_area(),
        nonlinear_excess_area=nl.excess_area,
        contact_angle_deg=nl.contact_angle_deg,
        force_residual=nl.force_residual)
