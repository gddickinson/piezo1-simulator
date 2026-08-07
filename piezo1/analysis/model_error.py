"""Systematic error from the choice of *model*, not from the data or a knob.

:mod:`piezo1.analysis.uncertainty` reports three kinds of spread and says on
every one of them that **model error is not captured**. Bootstrapping a sphere
fit tells you how well a sphere is determined; it cannot tell you whether a
sphere was the right shape to fit. This module estimates that missing term where
the project has a second defensible model to compare against.

Three places it can be done honestly:

* **dome geometry** — a sphere has one radius. PIEZO1's dome is not obviously
  spherical, so an oblate **spheroid** is fitted as well and the two curvatures
  compared. If the spheroid's equatorial and polar radii differ materially, the
  single number the project reports is an average over a shape it does not have.
* **elastic network** — the spring model is a choice, not a measurement.
  ``uniform``, ``inverse_square`` and ``inverse_sixth`` all appear in the
  literature and all are defensible. The gating-mode overlap is recomputed under
  each.
* **pore radius** — the pipeline uses the **Apollonius** convention, where each
  atom carries its own van der Waals radius. The commoner convention in the
  channel literature is a **uniform probe**: every atom the same size. Both are
  used in published work and they do not give the same number.

**The point is the comparison, not the spread.** A model spread far larger than
the sampling interval means the published confidence interval is measuring the
wrong thing — and this project reports several intervals of exactly that kind.
:func:`compare_with_sampling` puts the two side by side and names which one
dominates, which is the sentence a reader actually needs.

Deliberately *not* claimed: this is a lower bound on model error. Two models
disagreeing bounds it from below; two models agreeing does not bound it from
above, because both may be wrong in the same direction. That limitation is
stated on every result rather than left to be remembered.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["ModelError", "SpheroidFit", "fit_spheroid", "dome_model_error",
           "spring_model_error", "pore_convention_error",
           "compare_with_sampling"]


@dataclass
class ModelError:
    """Spread of a result over a choice of *model form*.

    Deliberately a separate type from :class:`~piezo1.analysis.uncertainty.
    Sensitivity`, which varies a knob within one model. Changing the spring
    exponent is a knob; changing a sphere into a spheroid is a different
    physical claim about the object, and conflating the two would let a
    "sensitivity range" quietly stand in for "we do not know the shape".
    """

    quantity: str
    values: dict = field(default_factory=dict)   # model name -> value
    unit: str = ""
    reference: str = ""                          # which model the project uses
    note: str = ""

    kind = "model error"

    @property
    def estimate(self) -> float:
        if self.reference and self.reference in self.values:
            return float(self.values[self.reference])
        return float(np.median(list(self.values.values())))

    @property
    def low(self) -> float:
        return float(min(self.values.values()))

    @property
    def high(self) -> float:
        return float(max(self.values.values()))

    @property
    def spread(self) -> float:
        return self.high - self.low

    @property
    def relative_spread(self) -> float:
        return self.spread / max(abs(self.estimate), 1e-12)

    def summary(self, digits: int = 3) -> str:
        parts = ", ".join(f"{k} {v:.{digits}f}" for k, v in self.values.items())
        return (f"{self.quantity}: {parts}{' ' + self.unit if self.unit else ''}"
                f" — model spread {self.spread:.{digits}f} "
                f"({self.relative_spread:.1%}). "
                f"A LOWER BOUND: agreement between two models does not bound "
                f"the error from above.")


# --------------------------------------------------------------------------
# Dome: a sphere, or a spheroid?
# --------------------------------------------------------------------------

@dataclass
class SpheroidFit:
    """An axis-aligned oblate spheroid fitted to a point cloud."""

    center: np.ndarray
    equatorial: float          # semi-axis perpendicular to the symmetry axis
    polar: float               # semi-axis along it
    rmse: float
    n_points: int

    @property
    def flattening(self) -> float:
        """0 for a sphere; positive when the polar axis is the shorter one."""
        return 1.0 - self.polar / self.equatorial

    @property
    def apex_curvature(self) -> float:
        """Radius of curvature at the pole — what a sphere fit approximates.

        For a spheroid this is ``a^2 / c``, which is *not* either semi-axis.
        Quoting the polar semi-axis instead is the obvious mistake and gives a
        number that is too small by the flattening.
        """
        return float(self.equatorial ** 2 / self.polar)


def fit_spheroid(points: np.ndarray, axis_direction=None) -> SpheroidFit:
    """Least-squares oblate spheroid with its symmetry axis along ``axis_direction``.

    Solved **exactly and linearly**, with no iteration. In the axis frame the
    implicit equation

        A(x1^2 + x2^2) + B x3^2 + C x1 + D x2 + E x3 + F = 0

    is linear in its six coefficients, and completing the square recovers the
    centre and both semi-axes. The null vector of the design matrix gives the
    coefficients up to scale, which is all that is needed.

    The first version of this alternated between the centre and the semi-axes
    with a hand-rolled gradient step, and it did not work: on a *known* spheroid
    with a full surface it returned a = 163 for a true 100, c = 98 for a true
    60, and put the centre 89 A away. Both axes were inflated by the same
    factor, which is the signature of a drifting centre rather than a shape
    error. Fitting a model with a broken fitter and reporting the disagreement
    as "model error" would have been a confident wrong answer.
    """
    points = np.asarray(points, dtype=np.float64)
    direction = (np.array([0.0, 0.0, 1.0]) if axis_direction is None
                 else np.asarray(axis_direction, dtype=np.float64))
    direction = direction / np.linalg.norm(direction)

    # Work in a frame whose third axis is the symmetry axis.
    helper = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(helper, direction))) > 0.9:
        helper = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(direction, helper)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(direction, e1)
    origin = points.mean(axis=0)
    local = np.column_stack([(points - origin) @ e1, (points - origin) @ e2,
                             (points - origin) @ direction])

    x1, x2, x3 = local[:, 0], local[:, 1], local[:, 2]
    design = np.column_stack([x1 ** 2 + x2 ** 2, x3 ** 2, x1, x2, x3,
                              np.ones(len(local))])
    # Smallest singular vector: the coefficients are defined only up to scale.
    _, _, vt = np.linalg.svd(design, full_matrices=False)
    a_, b_, c_, d_, e_, f_ = vt[-1]

    if a_ == 0 or b_ == 0:
        raise ValueError("degenerate spheroid fit: no quadratic term")
    p1, p2, p3 = -c_ / (2 * a_), -d_ / (2 * a_), -e_ / (2 * b_)
    rhs = a_ * (p1 ** 2 + p2 ** 2) + b_ * p3 ** 2 - f_
    if rhs / a_ <= 0 or rhs / b_ <= 0:
        raise ValueError("degenerate spheroid fit: not an ellipsoid")
    equatorial = float(np.sqrt(rhs / a_))
    polar = float(np.sqrt(rhs / b_))

    centre = origin + p1 * e1 + p2 * e2 + p3 * direction
    delta = points - centre
    along = delta @ direction
    radial = np.linalg.norm(delta - np.outer(along, direction), axis=1)
    # Geometric residual, not the implicit one: comparable with a sphere's RMSE.
    implicit = (radial / equatorial) ** 2 + (along / polar) ** 2 - 1.0
    gradient = np.sqrt((2 * radial / equatorial ** 2) ** 2
                       + (2 * along / polar ** 2) ** 2)
    residual = np.abs(implicit) / np.maximum(gradient, 1e-12)
    return SpheroidFit(center=centre, equatorial=equatorial, polar=polar,
                       rmse=float(np.sqrt(np.mean(residual ** 2))),
                       n_points=len(points))


def dome_model_error(surface: np.ndarray, axis, sphere_radius: float
                     ) -> ModelError:
    """Radius of curvature under a sphere and under an oblate spheroid.

    Both in nm, matching how the project reports dome geometry. The spheroid's
    contribution is its *apex* curvature, which is the quantity a sphere fit is
    approximating — not either semi-axis.
    """
    fit = fit_spheroid(np.asarray(surface, dtype=np.float64),
                       axis_direction=axis.direction)
    return ModelError(
        quantity="dome radius of curvature",
        values={"sphere": float(sphere_radius),
                "spheroid (apex)": fit.apex_curvature / 10.0},
        unit="nm", reference="sphere",
        note=(f"spheroid flattening {fit.flattening:+.3f}, "
              f"semi-axes {fit.equatorial / 10:.2f} / {fit.polar / 10:.2f} nm, "
              f"implicit RMSE {fit.rmse:.4f}"))


# --------------------------------------------------------------------------
# Elastic network: which spring model?
# --------------------------------------------------------------------------

def spring_model_error(blocks: list, displacement: np.ndarray,
                       n_modes: int = 20, cutoff: float = 15.0) -> ModelError:
    """Gating-mode overlap under each spring model in the literature.

    ``uniform``, ``inverse_square`` and ``inverse_sixth`` are all published and
    all defensible. The project reports the ``inverse_square`` number, so that
    is the reference; the spread across the three is the systematic term the
    bootstrap over atoms cannot see.
    """
    from ..physics.anm import SPRING_MODELS, ANM

    values = {}
    for spring in SPRING_MODELS:
        anm = ANM.from_trimer(blocks, cutoff=cutoff, spring=spring).build()
        modes = anm.calc_modes(n_modes=n_modes)
        values[spring] = float(modes.cumulative_overlap(displacement)[-1])
    return ModelError(
        quantity="cumulative gating overlap", values=values,
        reference="inverse_square",
        note=f"{n_modes} modes, {cutoff:.0f} A cutoff; all three springs appear "
             f"in the elastic-network literature")


# --------------------------------------------------------------------------
# Pore: whose radius does the probe use?
# --------------------------------------------------------------------------

def pore_convention_error(structure, axis, step: float | None = None,
                          uniform_radius: float = 1.7) -> ModelError:
    """Bottleneck radius under the Apollonius and uniform-probe conventions.

    The pipeline gives every atom its own van der Waals radius (Apollonius). The
    commoner convention in the channel literature — and what HOLE does — treats
    every atom as the same size. Both are published; they measure subtly
    different things, and the difference is a systematic error in whichever
    number is quoted.
    """
    if step is None:
        step = _P.value("pore.step")
    from ..structure.pore import pore_profile

    apollonius = pore_profile(structure, axis, step=step)

    # Same pipeline, every atom given one radius. Done by substituting the
    # radii rather than by writing a second profiler, so the only thing that
    # differs between the two numbers is the convention.
    original = type(structure).vdw_radii
    try:
        type(structure).vdw_radii = lambda self: np.full(
            self.n_atoms, uniform_radius, dtype=np.float32)
        uniform = pore_profile(structure, axis, step=step)
    finally:
        type(structure).vdw_radii = original

    return ModelError(
        quantity="pore bottleneck radius",
        values={"Apollonius (per-atom vdW)": float(apollonius.bottleneck_radius),
                f"uniform probe ({uniform_radius} A)":
                    float(uniform.bottleneck_radius)},
        unit="A", reference="Apollonius (per-atom vdW)",
        note=(f"probe {uniform_radius} A vs carbon's 1.70 A. On 7WLT the two "
              f"conventions agree EXACTLY at 1.70 because the bottleneck "
              f"lining is carbon; away from it the gap is just the offset "
              f"(0.30 A at both 1.40 and 2.00). So this is not a fixed "
              f"systematic error but a restatement of the probe radius, and "
              f"the Apollonius refinement buys nothing at a carbon-lined "
              f"constriction"))


# --------------------------------------------------------------------------
# The comparison that is the point
# --------------------------------------------------------------------------

def compare_with_sampling(model: ModelError, sampling) -> dict:
    """Put the model spread beside the sampling interval and say which wins.

    ``sampling`` is anything from :mod:`piezo1.analysis.uncertainty` carrying
    ``low`` and ``high``. When the model spread dominates, the confidence
    interval is answering a question nobody asked — it says how well the wrong
    shape is determined.
    """
    sampling_width = float(sampling.high - sampling.low)
    model_width = model.spread
    if sampling_width <= 0 and model_width <= 0:
        ratio = float("nan")
    elif sampling_width <= 0:
        ratio = float("inf")
    else:
        ratio = model_width / sampling_width

    dominant = ("model" if ratio > 1.0 else "sampling"
                if np.isfinite(ratio) else "model")
    return {
        "quantity": model.quantity,
        "model_spread": model_width,
        "sampling_width": sampling_width,
        "ratio": ratio,
        "dominant": dominant,
        "verdict": (
            f"{dominant} error dominates ({model_width:.4g} vs "
            f"{sampling_width:.4g}, {ratio:.1f}x)"
            + ("; the confidence interval measures how well one model is "
               "determined, not whether it is the right model"
               if dominant == "model" else "")),
    }
