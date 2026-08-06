"""Membrane mechanics: the Helfrich footprint around a curved inclusion.

PIEZO1's dome does not stop at the protein boundary. The bilayer around it is
bent too, and that **footprint** decays over a length λ = √(κ/γ) — about 14 nm
for a typical membrane. Haselwandter & MacKinnon's central result is that this
surrounding deformation, not the dome itself, dominates the channel's tension
sensitivity, because it is where most of the tension-dependent area lives.

The model is the linearised Helfrich (Canham–Helfrich) energy in the Monge
gauge, where the membrane mid-surface is described by a height field h(r):

.. math::

    E = \\int \\left[ \\frac{\\kappa}{2}(\\nabla^2 h)^2
                    + \\frac{\\gamma}{2}|\\nabla h|^2 \\right] dA

Minimising it gives the Euler–Lagrange equation

.. math::   \\kappa \\nabla^4 h - \\gamma \\nabla^2 h = 0

whose axisymmetric solution decaying at infinity is a modified Bessel function,
``h(r) = A K₀(r/λ)`` with ``λ = √(κ/γ)``. That exact solution is what the
numerical solver here is checked against.

**Units.** Lengths in nanometres, energies in k_BT, tension in k_BT/nm². The
conversion to the experimentalist's unit is ``1 k_BT/nm² = 4.114 mN/m`` at
298 K, provided by :func:`mnm_to_kt_per_nm2` and its inverse — mixing the two
is an easy and silent way to be wrong by a factor of four.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.special import k0, k1

from ..config import KT_ROOM

__all__ = ["MembraneParameters", "FootprintSolution", "decay_length",
           "analytic_profile", "analytic_energy", "solve_footprint",
           "footprint_energy", "kt_per_nm2_to_mnm", "mnm_to_kt_per_nm2",
           "KT_PER_NM2_IN_MNM"]

#: 1 k_BT/nm² expressed in mN/m at 298.15 K.
KT_PER_NM2_IN_MNM = KT_ROOM / 1e-18 * 1e3      # J/m^2 -> mN/m


def kt_per_nm2_to_mnm(value):
    """k_BT/nm² to mN/m. Accepts a scalar or an array."""
    out = np.asarray(value, dtype=float) * KT_PER_NM2_IN_MNM
    return float(out) if out.ndim == 0 else out


def mnm_to_kt_per_nm2(value):
    """mN/m to k_BT/nm². Accepts a scalar or an array."""
    out = np.asarray(value, dtype=float) / KT_PER_NM2_IN_MNM
    return float(out) if out.ndim == 0 else out


@dataclass
class MembraneParameters:
    """Elastic properties of the bilayer.

    Defaults are the values Haselwandter & MacKinnon 2018 used for PIEZO1:
    κ = 20 k_BT and γ = 0.1 k_BT/nm², which is 0.41 mN/m — a mildly tensioned
    membrane, not a resting one.
    """

    kappa: float = 20.0                  # bending modulus, k_BT
    tension: float = 0.1                 # k_BT/nm^2
    label: str = "Haselwandter & MacKinnon 2018"

    @classmethod
    def from_mnm(cls, kappa: float, tension_mnm: float,
                 label: str = "") -> "MembraneParameters":
        return cls(kappa=kappa, tension=mnm_to_kt_per_nm2(tension_mnm),
                   label=label or f"gamma = {tension_mnm} mN/m")

    @property
    def tension_mnm(self) -> float:
        return kt_per_nm2_to_mnm(self.tension)

    @property
    def decay_length(self) -> float:
        """λ = √(κ/γ), the range of the membrane footprint, in nm."""
        return decay_length(self.kappa, self.tension)


def decay_length(kappa: float, tension: float) -> float:
    """λ = √(κ/γ) in nm, for κ in k_BT and γ in k_BT/nm²."""
    if tension <= 0:
        return float("inf")
    return float(np.sqrt(kappa / tension))


# --------------------------------------------------------------------------
# Exact solution
# --------------------------------------------------------------------------

def analytic_profile(r: np.ndarray, r0: float, slope: float,
                     lam: float) -> np.ndarray:
    """Exact footprint height for a fixed contact slope at ``r0``.

    ``h(r) = A·K₀(r/λ)`` with ``A`` set by ``-h'(r0) = slope``. Note the sign
    convention: a positive ``slope`` means the membrane falls away from the
    inclusion, which is what a dome bulging towards the cytoplasm produces.
    """
    r = np.asarray(r, dtype=float)
    amp = slope * lam / k1(r0 / lam)
    return amp * k0(r / lam)


def analytic_energy(r0: float, slope: float, kappa: float,
                    tension: float) -> float:
    """Footprint energy in k_BT for a fixed contact slope.

    Because the Euler–Lagrange equation is satisfied everywhere, the Helfrich
    integral collapses to a boundary term at the inclusion. Writing the exact
    profile as ``h = A K₀(r/λ)`` with ``A = sλ/K₁(x)`` and ``x = r₀/λ``, the
    energy is ``π κ A²/λ² · x K₀(x) K₁(x)``, which simplifies to

    .. math::

        E = \\pi \\kappa\\, s^2 \\frac{r_0}{\\lambda}
            \\frac{K_0(r_0/\\lambda)}{K_1(r_0/\\lambda)}

    Note the direction of the Bessel ratio. Writing it the other way up —
    ``K₁/K₀`` — is an easy slip and gives an answer 2.5x too large at
    ``x = 0.7``, which is exactly the regime PIEZO1 sits in. It was caught here
    by integrating the functional over the *exact* profile and finding the two
    disagreed; that cross-check is now a test.
    """
    lam = decay_length(kappa, tension)
    x = r0 / lam
    return float(np.pi * kappa * slope ** 2 * x * k0(x) / k1(x))


# --------------------------------------------------------------------------
# Numerical solver
# --------------------------------------------------------------------------

@dataclass
class FootprintSolution:
    """Height field of the deformed membrane around an inclusion."""

    r: np.ndarray                 # nm, from the inclusion edge outwards
    h: np.ndarray                 # nm
    params: MembraneParameters
    r0: float = 0.0
    slope: float = 0.0
    energy: float = 0.0           # k_BT
    meta: dict = field(default_factory=dict)

    #: Slope above which the linearised Monge gauge stops being trustworthy.
    #: The expansion drops terms of order |grad h|^2 relative to 1, so at a
    #: slope of 0.5 the neglected terms are already ~25%.
    SMALL_SLOPE_LIMIT: float = 0.5

    @property
    def max_slope(self) -> float:
        return float(np.abs(np.gradient(self.h, self.r, edge_order=2)).max())

    @property
    def within_linear_regime(self) -> bool:
        return self.max_slope <= self.SMALL_SLOPE_LIMIT

    def validity_note(self) -> str:
        """Plain statement of whether the linear theory applies here.

        PIEZO1's measured dome meets the bilayer at a slope near 2, i.e. about
        63 degrees, which is far outside the small-slope expansion. The
        equations still solve and return numbers; those numbers are not
        quantitatively reliable, and saying so is more useful than reporting
        them as though they were.
        """
        s = self.max_slope
        if s <= self.SMALL_SLOPE_LIMIT:
            return (f"max slope {s:.2f} - within the small-slope regime, "
                    f"linear Helfrich theory applies")
        return (f"max slope {s:.2f} ({np.degrees(np.arctan(s)):.0f} degrees) "
                f"EXCEEDS the small-slope limit of {self.SMALL_SLOPE_LIMIT}. "
                f"The linearised Monge-gauge energy neglects terms of order "
                f"|grad h|^2, so energies and areas here are indicative of "
                f"scale and trend only, not quantitative. A full nonlinear "
                f"Helfrich or Euler-elastica treatment is needed for numbers.")

    @property
    def decay_length(self) -> float:
        return self.params.decay_length

    def excess_area(self) -> float:
        """Area stored in the deformation, relative to a flat annulus (nm²).

        To first order the excess area of a gently sloped height field is
        ``∫ ½|∇h|² dA``. This is the area released to the reservoir when the
        footprint flattens, and therefore the part of the tension response the
        footprint contributes.
        """
        dh = np.gradient(self.h, self.r, edge_order=2)
        return float(np.trapezoid(0.5 * dh ** 2 * 2.0 * np.pi * self.r, self.r))

    def fitted_decay_length(self) -> float:
        """Recover λ from the profile itself, by fitting ln(h√r) vs r.

        For large ``r``, ``K₀(r/λ) ~ e^{-r/λ}√(πλ/2r)``, so ``ln(h√r)`` is
        linear in ``r`` with slope ``-1/λ``. Fitting it back is an independent
        check that the solver reproduced the right physics rather than merely
        something smooth and decaying.
        """
        from scipy.optimize import curve_fit
        lam = self.decay_length
        # Least-squares fit of the actual K0 form over the region where the
        # numerical solution is well resolved. The asymptotic
        # ln(h*sqrt(r)) ~ -r/lam shortcut is tempting but fits the far tail,
        # which has fallen to ~1e-5 of the peak and is dominated by
        # discretisation noise; it returned decay lengths scattered between 9
        # and 60 nm for a true value of 14.
        window = (self.r >= self.r0) & (self.r < self.r0 + 4.0 * lam)
        r, h = self.r[window], self.h[window]
        if len(r) < 8 or np.abs(h).max() < 1e-12:
            return float("nan")

        def model(rr, amp, length):
            return amp * k0(rr / max(length, 1e-6))

        try:
            popt, _ = curve_fit(model, r, h, p0=[float(h[0] / k0(self.r0 / lam)), lam],
                                maxfev=20000)
        except Exception:
            return float("nan")
        return float(abs(popt[1]))

    def surface(self, n_angular: int = 96) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Revolve the radial profile into a Cartesian surface for rendering."""
        theta = np.linspace(0.0, 2 * np.pi, n_angular, endpoint=False)
        R, T = np.meshgrid(self.r, theta, indexing="ij")
        Z = np.repeat(self.h[:, None], n_angular, axis=1)
        return R * np.cos(T), R * np.sin(T), Z


def _radial_operators(r: np.ndarray) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    """First and second derivative matrices on a uniform radial grid."""
    n = len(r)
    dr = r[1] - r[0]
    d1 = sp.diags([-1.0, 1.0], [-1, 1], shape=(n, n), format="lil") / (2 * dr)
    d2 = sp.diags([1.0, -2.0, 1.0], [-1, 0, 1], shape=(n, n), format="lil") / dr ** 2
    # One-sided differences at the ends; the boundary rows are overwritten by
    # the boundary conditions anyway, but leaving them singular breaks the
    # operator product used to build the biharmonic term.
    # Second-order one-sided stencils. A first-order (h1-h0)/dr boundary
    # condition drags the whole solution down to first order even though the
    # interior scheme is second order.
    d1[0, 0], d1[0, 1], d1[0, 2] = -3.0 / (2 * dr), 4.0 / (2 * dr), -1.0 / (2 * dr)
    d1[-1, -3], d1[-1, -2], d1[-1, -1] = 1.0 / (2 * dr), -4.0 / (2 * dr), 3.0 / (2 * dr)
    d2[0, 0], d2[0, 1], d2[0, 2] = 1.0 / dr ** 2, -2.0 / dr ** 2, 1.0 / dr ** 2
    d2[-1, -3], d2[-1, -2], d2[-1, -1] = 1.0 / dr ** 2, -2.0 / dr ** 2, 1.0 / dr ** 2
    return d1.tocsr(), d2.tocsr()


def solve_footprint(r0: float, slope: float,
                    params: MembraneParameters | None = None,
                    r_max: float | None = None, n: int = 1200,
                    height: float | None = None) -> FootprintSolution:
    """Solve κ∇⁴h − γ∇²h = 0 on an annulus by finite differences.

    The fourth-order equation is solved as a **coupled pair of second-order
    equations** in ``h`` and the mean curvature ``u = ∇²h``:

    .. math::   \\nabla^2 h - u = 0, \\qquad \\kappa \\nabla^2 u - \\gamma u = 0

    Forming the biharmonic operator directly as ``L @ L`` is the obvious
    approach and it does not work: squaring the discrete Laplacian squares its
    condition number, and the resulting solution converged to a profile with a
    47 nm decay length where the exact answer is 14 nm, with a 59% energy error
    that did not improve on grid refinement. The coupled form keeps every
    operator second-order and converges cleanly.

    Parameters
    ----------
    r0:
        Radius of the inclusion, nm. For PIEZO1 the footprint starts at the
        edge of the dome, so this is the dome's projected radius.
    slope:
        Contact slope at ``r0``, dimensionless (``-dh/dr``). Positive means the
        membrane falls away from the inclusion.
    height:
        Contact height at ``r0``, nm. Defaults to the value the exact
        decaying solution takes for this slope, which is the physically
        consistent choice for an isolated inclusion.
    r_max:
        Outer radius. Defaults to ``r0 + 12λ``, far enough that the K₀ tail has
        decayed to ~1e-5 of its boundary value.
    """
    p = params or MembraneParameters()
    lam = p.decay_length
    if not np.isfinite(lam):
        raise ValueError("zero tension gives an infinite footprint; "
                         "the linearised model does not apply")
    if r_max is None:
        r_max = r0 + 12.0 * lam

    r = np.linspace(r0, r_max, n)
    d1, d2 = _radial_operators(r)
    lap = (d2 + sp.diags(1.0 / r) @ d1).tolil()

    amp = slope * lam / k1(r0 / lam)
    if height is None:
        height = float(amp * k0(r0 / lam))

    eye = sp.identity(n, format="lil")
    zero = sp.lil_matrix((n, n))

    # Block rows: [ lap  -I ] [h]   [0]
    #             [ 0    K  ] [u] = [0]   with K = kappa*lap - gamma*I
    a11, a12 = lap.copy(), (-eye).tolil()
    a21, a22 = zero.copy(), (p.kappa * lap - p.tension * eye).tolil()
    rhs = np.zeros(2 * n)

    def set_row(block_row, other, i: int, row_h, row_u, value: float) -> None:
        block_row[i, :] = 0.0
        other[i, :] = 0.0
        if row_h is not None:
            block_row[i, :] = row_h
        if row_u is not None:
            other[i, :] = row_u
        rhs[i if block_row is a11 else n + i] = value

    # Two conditions at the inclusion, two in the far field.
    set_row(a11, a12, 0, np.eye(1, n, 0).ravel(), None, height)
    set_row(a11, a12, n - 1, np.eye(1, n, n - 1).ravel(), None, 0.0)
    set_row(a22, a21, 0, None, None, -slope)
    a21[0, :] = d1[0, :].toarray().ravel()
    set_row(a22, a21, n - 1, None, None, 0.0)
    a21[n - 1, :] = d1[n - 1, :].toarray().ravel()

    operator = sp.bmat([[a11, a12], [a21, a22]], format="csr")
    solution = spla.spsolve(operator, rhs)
    h = solution[:n]

    sol = FootprintSolution(r=r, h=h, params=p, r0=r0, slope=slope,
                            meta={"n": n, "r_max": r_max, "lambda": lam,
                                  "contact_height": height,
                                  "formulation": "coupled second order"})
    sol.energy = footprint_energy(sol)
    return sol


def footprint_energy(solution: FootprintSolution) -> float:
    """Helfrich energy of a footprint solution, in k_BT, by direct integration.

    Integrates the functional itself rather than using the boundary-term
    shortcut, so it remains valid for profiles that do not satisfy the
    Euler–Lagrange equation — a numerically imperfect solution, or a profile
    imported from simulation.
    """
    r, h = solution.r, solution.h
    p = solution.params
    dh = np.gradient(h, r, edge_order=2)
    d2h = np.gradient(dh, r, edge_order=2)
    laplacian = d2h + dh / np.maximum(r, 1e-12)
    density = 0.5 * p.kappa * laplacian ** 2 + 0.5 * p.tension * dh ** 2
    return float(np.trapezoid(density * 2.0 * np.pi * r, r))
