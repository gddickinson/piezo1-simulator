"""Re-derive each headline result by a route that shares no machinery with it.

Round 18 is the reason this exists. The nonlinear footprint solver integrated
Euler–Lagrange equations derived by hand, and a boundary-value solver converges
happily onto wrong equations. What caught the error was evaluating the *exact
functional in a different gauge* — a check that reused none of the derivation
being checked.

That generalises. A test written from the same understanding as the code shares
its blind spots: if the author misread the physics, the test encodes the
misreading. An independent re-derivation does not, and where two routes agree
the agreement means something.

Each function here answers a question the main pipeline already answers, and
answers it a different way:

* **dome curvature** — the pipeline fits a four-parameter sphere by algebraic
  least squares. Here the curvature comes from the coefficient of a parabola
  through the radial height profile about the C3 axis. Different parametrisation,
  different fit, different failure modes.
* **mode overlap** — the pipeline superposes the two structures with Kabsch and
  projects a Cartesian displacement onto the modes. Here the comparison is made
  in **pairwise distances**, which are invariant to rotation and translation, so
  no superposition happens at all. A superposition bug cannot survive this.
* **T₅₀** — the pipeline builds a four-state generator, exponentiates it and
  bisects on the peak of the response. Here the same rate constants are fed to
  an analytic steady state obtained by solving the linear system directly.

Disagreement is the point. Agreement is worth recording precisely because these
routes could have disagreed.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P
from ..physics.kinetics import STATE_NAMES as _STATE_NAMES

#: Index of the open state, derived from the scheme's own state names rather
#: than written as a literal — if the state order ever changes, a hardcoded
#: index would silently start reading the wrong occupancy.
OPEN_STATE = _STATE_NAMES.index("O")

__all__ = ["CrossCheck", "dome_curvature_by_parabola",
           "dome_curvature_by_cap_geometry",
           "gating_overlap_by_distances", "t50_by_steady_state",
           "t50_by_ode_integration",
           "compare"]


@dataclass
class CrossCheck:
    """Two routes to the same quantity, and how far apart they land."""

    quantity: str
    primary: float
    alternative: float
    unit: str = ""
    tolerance: float = 0.10          # relative, for the verdict only
    primary_route: str = ""
    alternative_route: str = ""
    note: str = ""

    @property
    def difference(self) -> float:
        return self.alternative - self.primary

    @property
    def relative(self) -> float:
        scale = max(abs(self.primary), 1e-12)
        return abs(self.difference) / scale

    @property
    def agrees(self) -> bool:
        return self.relative <= self.tolerance

    def summary(self) -> str:
        mark = "agree" if self.agrees else "DISAGREE"
        return (f"{mark:8s} {self.quantity:34s} "
                f"{self.primary:9.4f} vs {self.alternative:9.4f} "
                f"{self.unit:5s} ({self.relative:.1%})")


# --------------------------------------------------------------------------
# Dome curvature without fitting a sphere
# --------------------------------------------------------------------------

def dome_curvature_by_parabola(surface: np.ndarray, axis) -> float:
    """Radius of curvature from the radial height profile, in nm.

    For a spherical cap of radius ``R`` the height near the apex obeys
    ``h(r) = h₀ − r²/(2R)``, so a straight-line fit of ``h`` against ``r²``
    gives ``R = −1/(2·slope)``. No sphere is fitted, no centre is solved for,
    and the four-parameter algebraic linearisation the main pipeline uses plays
    no part.

    The trade is that this uses the C3 axis, which the sphere fit does not — so
    the two routes fail in different ways rather than one being strictly better.
    """
    surface = np.asarray(surface, dtype=np.float64)
    height = axis.project(surface)
    radial = axis.radial(surface)

    # Orient so the cap opens the same way regardless of the axis direction.
    if np.polyfit(radial ** 2, height, 1)[0] > 0:
        height = -height

    slope = np.polyfit(radial ** 2, height, 1)[0]
    if slope >= 0:
        return float("nan")
    return float(-1.0 / (2.0 * slope) / 10.0)


def dome_curvature_by_cap_geometry(surface: np.ndarray, axis) -> float:
    """Radius from the exact spherical-cap relation, in nm.

    A point on a cap of radius ``R`` with apex at the origin satisfies
    ``h² + 2hR + r² = 0``, so ``R = −(h² + r²)/(2h)`` exactly, per point. No
    expansion, no four-parameter fit — just the algebra inverted and the
    per-point estimates combined by median, which is robust to the few surface
    points that sit off the fitted cap.

    This exists because :func:`dome_curvature_by_parabola` turned out to be a
    *small-slope* approximation, and PIEZO1's dome is not shallow. On synthetic
    caps of known radius the parabola is 0.6% low at an 8.6° contact angle and
    **25.8% low at 63.4°**, which is where PIEZO1 sits. A cross-check that
    disagrees because the checking route is the invalid one is still useful —
    but only once you know that is why.
    """
    surface = np.asarray(surface, dtype=np.float64)
    height = axis.project(surface)
    radial = axis.radial(surface)
    if np.polyfit(radial ** 2, height, 1)[0] > 0:
        height = -height
    height = height - height.max()          # put the apex at zero

    usable = height < -1e-6
    if usable.sum() < 4:
        return float("nan")
    estimates = -(height[usable] ** 2 + radial[usable] ** 2) / (2.0 * height[usable])
    return float(np.median(estimates) / 10.0)


# --------------------------------------------------------------------------
# Mode overlap without superposing anything
# --------------------------------------------------------------------------

def _pair_indices(n_sites: int, n_pairs: int, cutoff_lo: float,
                  coords: np.ndarray, seed: int = 0):
    """Site pairs far enough apart to carry information about the motion."""
    rng = np.random.default_rng(seed)
    i = rng.integers(0, n_sites, n_pairs * 3)
    j = rng.integers(0, n_sites, n_pairs * 3)
    keep = i != j
    i, j = i[keep], j[keep]
    separation = np.linalg.norm(coords[i] - coords[j], axis=1)
    good = separation > cutoff_lo
    return i[good][:n_pairs], j[good][:n_pairs]


def gating_overlap_by_distances(closed: np.ndarray, open_: np.ndarray,
                                modes, mode_index: int | None = None,
                                n_pairs: int = 4000, seed: int = 0) -> float:
    """Overlap computed in pairwise distances, so no superposition is needed.

    Distances between sites are invariant to rotation and translation. The
    observed transition changes them by ``Δd_ij``; a normal mode changes them,
    to first order, by ``(u_i − u_j)·ê_ij``. The correlation between those two
    vectors is an overlap that never superposes anything — so a protomer
    mismatch or a Kabsch sign error cannot pass through it.

    Returns the best absolute correlation over the A-symmetric modes, matching
    what the main pipeline reports.
    """
    closed = np.asarray(closed, dtype=np.float64)
    open_ = np.asarray(open_, dtype=np.float64)
    if closed.shape != open_.shape:
        raise ValueError(f"shape mismatch: {closed.shape} vs {open_.shape}")

    i, j = _pair_indices(len(closed), n_pairs,
                         _P.value("crosscheck.min_pair_separation"),
                         closed, seed)
    delta = closed[i] - closed[j]
    distance = np.linalg.norm(delta, axis=1)
    unit = delta / distance[:, None]
    observed = np.linalg.norm(open_[i] - open_[j], axis=1) - distance

    # `symmetry` is a numpy array; `x or []` on one raises "truth value of an
    # array is ambiguous". The project has hit this before.
    raw = getattr(modes, "symmetry", None)
    symmetry = [] if raw is None else list(raw)
    if mode_index is not None:
        candidates = [mode_index]
    else:
        candidates = [k for k, s in enumerate(symmetry) if s == "A"]
        if not candidates:
            candidates = list(range(modes.n_modes))

    best = 0.0
    for k in candidates:
        vectors = modes.vectors[k]
        predicted = np.einsum("ij,ij->i", vectors[i] - vectors[j], unit)
        if predicted.std() < 1e-15 or observed.std() < 1e-15:
            continue
        best = max(best, abs(float(np.corrcoef(predicted, observed)[0, 1])))
    return best


# --------------------------------------------------------------------------
# T50 without exponentiating a generator
# --------------------------------------------------------------------------

def t50_by_steady_state(model, lo: float = 0.0, hi: float = 20.0,
                        n: int = 4000) -> float:
    """Half-activation from the analytic steady state, in mN/m.

    The main pipeline integrates the master equation with a matrix exponential
    and bisects on the *peak* of the transient. Here the steady-state
    occupancy is obtained by solving ``Qᵀp = 0`` with a normalisation row —
    a linear solve, no exponential, no time stepping — and the half-maximum is
    read off that curve.

    These are **not the same quantity**: a peak response and an equilibrium
    occupancy differ whenever inactivation is fast. The check is therefore not
    that they match to three digits but that they land in the same place, which
    would fail immediately if the rate matrix were assembled wrongly.
    """
    tension = np.linspace(lo, hi, n)
    occupancy = np.empty(n)
    for index, value in enumerate(tension):
        generator = np.asarray(model.rate_matrix(float(value)), dtype=float)
        size = generator.shape[0]
        # Steady state: the left null vector of Q, pinned by sum(p) = 1.
        system = np.vstack([generator.T[:-1], np.ones(size)])
        target = np.zeros(size)
        target[-1] = 1.0
        try:
            solution = np.linalg.lstsq(system, target, rcond=None)[0]
        except np.linalg.LinAlgError:
            occupancy[index] = np.nan
            continue
        occupancy[index] = solution[OPEN_STATE]

    good = np.isfinite(occupancy)
    if not good.any():
        return float("nan")
    peak = occupancy[good].max()
    if peak <= 0:
        return float("nan")
    reached = np.flatnonzero(good & (occupancy >= 0.5 * peak))
    if len(reached) == 0:
        return float("nan")
    first = int(reached[0])
    if first == 0:
        return float(tension[0])
    x0, x1 = tension[first - 1], tension[first]
    y0, y1 = occupancy[first - 1], occupancy[first]
    return float(x0 + (0.5 * peak - y0) * (x1 - x0) / max(y1 - y0, 1e-15))


def t50_by_ode_integration(model, lo: float = 0.5, hi: float = 12.0,
                           n: int = 40, duration: float = 0.5) -> float:
    """Half-activation from direct ODE integration, in mN/m.

    The main pipeline forms ``expm(Qt)`` and reads the peak occupancy off the
    resulting trajectory. Here the same master equation is integrated by an
    adaptive Runge–Kutta solver instead. Same quantity, different numerics, no
    matrix exponential — so an error in the exponential or in how the
    trajectory is sampled cannot survive both.

    This replaced an earlier attempt that used the analytic steady state, which
    turned out to compute a **different quantity**: at equilibrium this channel
    sits ~96% inactivated at every tension, so steady-state open occupancy runs
    only 0.030 to 0.036 and has no half-maximum. T₅₀ is necessarily a property
    of the peak transient — which is also what a patch-clamp measures.
    """
    from scipy.integrate import solve_ivp

    tensions = np.linspace(lo, hi, n)
    peaks = np.empty(n)
    for index, tension in enumerate(tensions):
        generator = np.asarray(model.rate_matrix(float(tension)), dtype=float)
        # Start from the resting distribution, as the step protocol does.
        resting = np.asarray(model.rate_matrix(0.0), dtype=float)
        size = resting.shape[0]
        system = np.vstack([resting.T[:-1], np.ones(size)])
        target = np.zeros(size)
        target[-1] = 1.0
        start = np.linalg.lstsq(system, target, rcond=None)[0]

        solution = solve_ivp(lambda _t, p: generator.T @ p, (0.0, duration),
                             start, method="LSODA", rtol=1e-9, atol=1e-12,
                             dense_output=True)
        times = np.linspace(0.0, duration, 2000)
        peaks[index] = solution.sol(times)[OPEN_STATE].max()

    target_level = 0.5 * peaks.max()
    reached = np.flatnonzero(peaks >= target_level)
    if len(reached) == 0:
        return float("nan")
    first = int(reached[0])
    if first == 0:
        return float(tensions[0])
    x0, x1 = tensions[first - 1], tensions[first]
    y0, y1 = peaks[first - 1], peaks[first]
    return float(x0 + (target_level - y0) * (x1 - x0) / max(y1 - y0, 1e-15))


def compare(checks: list[CrossCheck], verbose: bool = True) -> list[CrossCheck]:
    if verbose:
        for check in checks:
            print("  " + check.summary())
            if check.note:
                print(f"           {check.note}")
    return checks
