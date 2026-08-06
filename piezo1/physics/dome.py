"""Dome-model energetics: turning an area change into a gating probability.

PIEZO1 is gated by an **area change**. Its curved blades bend the bilayer into
a dome whose projected in-plane area is smaller than its surface area;
flattening the dome under tension increases that projected area, and membrane
tension acting through the change does work on the channel:

.. math::   \\Delta E = -T \\cdot \\Delta A

Combined with an intrinsic free-energy difference ΔG₀ favouring the closed
state, the open probability follows a two-state Boltzmann:

.. math::   P_{open}(T) = \\frac{1}{1 + \\exp[(\\Delta G_0 - T \\Delta A)/k_BT]}

so the half-activation tension is simply ``T₅₀ = ΔG₀ / ΔA``.

**The discrepancy this module makes explicit.** ΔA measured two different ways
disagrees by more than an order of magnitude:

* **Functionally**, from Boltzmann fits to tension–response curves: 6–20 nm²
  (Cox et al. 2016 give 8 ± 1 nm² with ΔG₀ = 9.7 ± 1.5 k_BT).
* **Structurally**, from the dome geometry: ~40 nm² of nanodome excess area
  (Dixit 2025), ~120 nm² projected (Guo & MacKinnon 2017), up to ~300 nm² of
  in-plane expansion (Yang 2022).

These are *different quantities* and must not be quoted interchangeably. The
functional ΔA is the area change along the reaction coordinate between the two
gating states; the structural numbers measure the whole protein-plus-footprint
deformation, only part of which is coupled to the gate. :meth:`DomeModel.report`
prints both and the resulting T₅₀ for each, so the gap is visible rather than
buried.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from .membrane import (MembraneParameters, kt_per_nm2_to_mnm,
                       mnm_to_kt_per_nm2, solve_footprint)
from ..parameters import PARAMETERS as _P

__all__ = ["DomeModel", "DomeGeometrySummary", "gating_energy",
           "open_probability", "half_activation_tension",
           "PUBLISHED_AREA_ESTIMATES"]


#: Published ΔA estimates, in nm², with what each actually measures.
PUBLISHED_AREA_ESTIMATES = {
    "cox2016_functional": {
        "value": 8.0, "uncertainty": 1.0, "kind": "functional",
        "note": "Boltzmann fit to a tension-response curve; the area change "
                "along the gating coordinate.",
        "source": "Cox et al. Nat Commun 2016 (PMID 26785635)"},
    "dixit2025_nanodome": {
        "value": 40.0, "uncertainty": None, "kind": "structural",
        "note": "Simulated nanodome excess area at zero tension.",
        "source": "Dixit, Noe & Weikl eLife 2025"},
    "guo2017_projected": {
        "value": 120.0, "uncertainty": None, "kind": "structural",
        "note": "Difference between mid-plane and projected dome area.",
        "source": "Guo & MacKinnon eLife 2017 (PMID 29231809)"},
    "yang2022_inplane": {
        "value": 300.0, "uncertainty": None, "kind": "structural",
        "note": "In-plane areal expansion on flattening.",
        "source": "Yang et al. Nature 2022 (PMID 35388220)"},
}


@dataclass
class DomeGeometrySummary:
    """The geometric numbers a :class:`DomeModel` needs, in nm."""

    radius_of_curvature: float
    dome_area: float                 # nm^2, curved mid-membrane surface
    projected_area: float            # nm^2, its shadow on the membrane plane
    dome_depth: float = 0.0
    footprint_radius: float = 0.0

    @property
    def excess_area(self) -> float:
        """Area released if the dome flattened completely, nm²."""
        return self.dome_area - self.projected_area

    @classmethod
    def from_measurement(cls, dome) -> "DomeGeometrySummary":
        """Build from a :class:`piezo1.structure.geometry.DomeGeometry`.

        That object reports Angstrom; this one works in nm, matching the
        membrane literature. The conversion happens here, once.
        """
        return cls(radius_of_curvature=dome.radius_of_curvature / 10.0,
                   dome_area=dome.dome_area / 100.0,
                   projected_area=dome.projected_area / 100.0,
                   dome_depth=dome.dome_depth / 10.0,
                   footprint_radius=dome.footprint_radius / 10.0)


# --------------------------------------------------------------------------
# Two-state energetics
# --------------------------------------------------------------------------

def gating_energy(tension: float | np.ndarray, delta_area: float,
                  delta_g0: float) -> np.ndarray:
    """ΔG(T) = ΔG₀ − T·ΔA, in k_BT. Tension in k_BT/nm², area in nm²."""
    return np.asarray(delta_g0 - np.asarray(tension, float) * delta_area)


def open_probability(tension: float | np.ndarray, delta_area: float,
                     delta_g0: float) -> np.ndarray:
    """Two-state Boltzmann open probability."""
    dg = gating_energy(tension, delta_area, delta_g0)
    return 1.0 / (1.0 + np.exp(np.clip(dg, -500, 500)))


def half_activation_tension(delta_area: float, delta_g0: float) -> float:
    """T₅₀ = ΔG₀/ΔA, in k_BT/nm². Where ΔG = 0 and P_open = ½."""
    if delta_area <= 0:
        return float("inf")
    return float(delta_g0 / delta_area)


# --------------------------------------------------------------------------
# The model
# --------------------------------------------------------------------------

@dataclass
class DomeModel:
    """Tension-driven gating of a dome-shaped mechanosensor.

    Defaults are Cox et al. 2016's functionally determined parameters, which
    are the ones that reproduce measured tension-response curves.
    """

    delta_area: float = field(
        default_factory=lambda: _P.value("dome.delta_area"))              # nm^2, area change on gating
    delta_g0: float = field(
        default_factory=lambda: _P.value("dome.delta_g0"))                # k_BT, intrinsic closed-state bias
    geometry: DomeGeometrySummary | None = None
    membrane: MembraneParameters = field(default_factory=MembraneParameters)
    label: str = "Cox et al. 2016 functional parameters"
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------ responses

    def open_probability(self, tension_mnm: float | np.ndarray) -> np.ndarray:
        t = mnm_to_kt_per_nm2(np.asarray(tension_mnm, float))
        return open_probability(t, self.delta_area, self.delta_g0)

    def gating_energy(self, tension_mnm: float | np.ndarray) -> np.ndarray:
        t = mnm_to_kt_per_nm2(np.asarray(tension_mnm, float))
        return gating_energy(t, self.delta_area, self.delta_g0)

    @property
    def half_activation_mnm(self) -> float:
        """T₅₀ in mN/m — the number a patch-clamp experiment reports."""
        return kt_per_nm2_to_mnm(
            half_activation_tension(self.delta_area, self.delta_g0))

    def sensitivity(self, tension_mnm: float) -> float:
        """dP_open/dT at a given tension, per mN/m. Steepest at T₅₀."""
        p = float(self.open_probability(tension_mnm))
        return float(p * (1.0 - p) * self.delta_area / kt_per_nm2_to_mnm(1.0))

    # -------------------------------------------------------- the footprint

    def footprint(self, contact_slope: float | None = None,
                  n: int = 1200):
        """Solve the surrounding membrane deformation for this dome.

        The contact slope is derived from the dome geometry if not given: for a
        spherical cap of radius ``R`` meeting the plane at in-plane radius
        ``a``, the boundary slope is ``a/√(R² − a²)``.
        """
        if self.geometry is None:
            raise ValueError("no geometry attached; measure a dome first")
        g = self.geometry
        a = np.sqrt(max(g.projected_area, 0.0) / np.pi)
        if contact_slope is None:
            r = g.radius_of_curvature
            contact_slope = float(a / np.sqrt(max(r * r - a * a, 1e-9)))
        return solve_footprint(a, contact_slope, self.membrane, n=n)

    def footprint_area(self, **kw) -> float:
        """Excess area stored in the surrounding membrane, nm² (linear theory).

        **Not quantitative for PIEZO1.** The measured dome meets the bilayer at
        about 63°, where the small-slope expansion behind this number breaks
        down; it overestimates by ~3.5×. Use :meth:`footprint_area_nonlinear`.
        Kept because it is what the closed-form Bessel solution gives and the
        comparison between the two is itself informative.
        """
        return self.footprint(**kw).excess_area()

    def contact_slope(self) -> float:
        """Slope where the dome meets the flat bilayer, from the cap geometry."""
        if self.geometry is None:
            raise ValueError("no geometry attached; measure a dome first")
        g = self.geometry
        a = np.sqrt(max(g.projected_area, 0.0) / np.pi)
        r = g.radius_of_curvature
        return float(a / np.sqrt(max(r * r - a * a, 1e-9)))

    def footprint_nonlinear(self, contact_slope: float | None = None,
                            n: int = 600):
        """Solve the footprint without the small-slope approximation.

        This is the quantitative version: an arc-length Euler–elastica solve
        with exact principal curvatures. See :mod:`piezo1.physics.elastica`.
        """
        from .elastica import solve_elastica
        if self.geometry is None:
            raise ValueError("no geometry attached; measure a dome first")
        a = np.sqrt(max(self.geometry.projected_area, 0.0) / np.pi)
        slope = self.contact_slope() if contact_slope is None else contact_slope
        return solve_elastica(a, slope, self.membrane, n=n)

    def footprint_area_nonlinear(self, **kw) -> float:
        """Excess area stored in the surrounding membrane, nm², exactly."""
        return self.footprint_nonlinear(**kw).excess_area

    def gating_area_change(self, open_geometry: "DomeGeometrySummary",
                           nonlinear: bool = True) -> dict:
        """ΔA between this (closed) dome and an open one, with the footprint.

        The two-state model is driven by the **change** in projected in-plane
        area, not by any absolute area. Two things contribute:

        * the dome's own projected area grows as it flattens;
        * the surrounding footprint flattens too, releasing the excess area it
          was storing into projected area.

        The second term is where Round 18's correction bites. At the closed
        state's 63° contact slope the linearised footprint overstates the
        stored area 3.5×, and because the open state sits at 40° — where the
        linear theory is much less wrong — the *difference* is overstated far
        more than either endpoint. Passing ``nonlinear=False`` reproduces the
        linear route for comparison.
        """
        if self.geometry is None:
            raise ValueError("no geometry attached; measure a dome first")
        from .elastica import compare_with_linear

        closed = self.geometry
        dome_term = open_geometry.projected_area - closed.projected_area

        def footprint_area(geometry) -> tuple[float, float]:
            radius = float(np.sqrt(max(geometry.projected_area, 0.0) / np.pi))
            slope = float(radius / np.sqrt(
                max(geometry.radius_of_curvature ** 2 - radius ** 2, 1e-9)))
            comparison = compare_with_linear(radius, slope, self.membrane)
            return (comparison.nonlinear_excess_area,
                    comparison.linear_excess_area)

        closed_nl, closed_lin = footprint_area(closed)
        open_nl, open_lin = footprint_area(open_geometry)
        footprint_term = ((closed_nl - open_nl) if nonlinear
                          else (closed_lin - open_lin))

        total = dome_term + footprint_term
        return {
            "delta_area_nm2": total,
            "dome_term_nm2": dome_term,
            "footprint_term_nm2": footprint_term,
            "footprint_term_linear_nm2": closed_lin - open_lin,
            "footprint_term_nonlinear_nm2": closed_nl - open_nl,
            "closed_footprint_nm2": closed_nl if nonlinear else closed_lin,
            "open_footprint_nm2": open_nl if nonlinear else open_lin,
            "nonlinear": nonlinear,
            "t50_mnm": kt_per_nm2_to_mnm(
                half_activation_tension(total, self.delta_g0)),
        }

    def compare_gating_area_routes(self, open_geometry) -> list[dict]:
        """Every route to ΔA side by side, with the T₅₀ each implies.

        The point of the table is that improving the membrane physics does
        **not** close the structural-versus-functional gap. It is a real
        correction in the right direction and it is nowhere near enough, which
        says the gap is about *what is being measured* rather than about the
        accuracy of the membrane model.
        """
        nonlinear = self.gating_area_change(open_geometry, nonlinear=True)
        linear = self.gating_area_change(open_geometry, nonlinear=False)
        rows = [
            {"route": "functional (Cox 2016 fit)",
             "delta_area_nm2": self.delta_area,
             "t50_mnm": self.half_activation_mnm,
             "note": "the number that reproduces measured tension-response"},
            {"route": "structural, dome projected area only",
             "delta_area_nm2": nonlinear["dome_term_nm2"],
             "t50_mnm": kt_per_nm2_to_mnm(half_activation_tension(
                 nonlinear["dome_term_nm2"], self.delta_g0)),
             "note": "ignores the surrounding membrane entirely"},
            {"route": "structural, dome + linear footprint",
             "delta_area_nm2": linear["delta_area_nm2"],
             "t50_mnm": linear["t50_mnm"],
             "note": "the pre-Round-18 route; the footprint term is 6.5x too "
                     "large"},
            {"route": "structural, dome + nonlinear footprint",
             "delta_area_nm2": nonlinear["delta_area_nm2"],
             "t50_mnm": nonlinear["t50_mnm"],
             "note": "correct membrane physics; still far from the functional "
                     "number"},
        ]
        return rows

    def compare_footprint_theories(self) -> dict:
        """Linear vs nonlinear footprint for this dome, side by side.

        Reports the correction factor, which for PIEZO1's geometry is large
        enough to change a qualitative conclusion, not merely a digit.
        """
        from .elastica import compare_with_linear
        if self.geometry is None:
            raise ValueError("no geometry attached; measure a dome first")
        a = np.sqrt(max(self.geometry.projected_area, 0.0) / np.pi)
        c = compare_with_linear(a, self.contact_slope(), self.membrane)
        return {
            "inclusion_radius_nm": a,
            "contact_slope": c.slope,
            "contact_angle_deg": c.contact_angle_deg,
            "linear_energy_kT": c.linear_energy,
            "nonlinear_energy_kT": c.nonlinear_energy,
            "linear_excess_area_nm2": c.linear_excess_area,
            "nonlinear_excess_area_nm2": c.nonlinear_excess_area,
            "linear_overestimate_energy": 1.0 / c.energy_ratio,
            "linear_overestimate_area": 1.0 / c.area_ratio,
            "dome_excess_area_nm2": self.geometry.excess_area,
            "axial_force_residual": c.force_residual,
        }

    # -------------------------------------------------------------- report

    def compare_area_estimates(self) -> list[dict]:
        """T₅₀ implied by each published ΔA, holding ΔG₀ fixed.

        The point of this table is that the structural areas, taken at face
        value as gating areas, predict half-activation tensions far below what
        is measured. They are measuring something else.
        """
        rows = []
        for key, rec in PUBLISHED_AREA_ESTIMATES.items():
            t50 = kt_per_nm2_to_mnm(
                half_activation_tension(rec["value"], self.delta_g0))
            rows.append({"key": key, "delta_area_nm2": rec["value"],
                         "kind": rec["kind"], "t50_mnm": t50,
                         "source": rec["source"], "note": rec["note"]})
        return rows

    def report(self) -> str:
        lines = [
            f"Dome model: {self.label}",
            f"  dA = {self.delta_area:g} nm^2, dG0 = {self.delta_g0:g} kT",
            f"  -> T50 = {self.half_activation_mnm:.2f} mN/m",
        ]
        if self.geometry is not None:
            g = self.geometry
            lines += [
                f"  geometry: R_c = {g.radius_of_curvature:.1f} nm, "
                f"dome {g.dome_area:.0f} nm^2, projected {g.projected_area:.0f} nm^2, "
                f"excess {g.excess_area:.0f} nm^2",
            ]
        lines.append("  T50 implied by each published area estimate:")
        for row in self.compare_area_estimates():
            lines.append(f"    {row['key']:24s} dA={row['delta_area_nm2']:6.1f} nm^2 "
                         f"({row['kind']:10s}) -> T50 = {row['t50_mnm']:6.2f} mN/m")
        lines.append("  Measured T50: 2.7 +/- 0.1 (cell-attached), "
                     "4.7 +/- 0.3 (inside-out), 5.1 +/- 0.2 mN/m")
        return "\n".join(lines)
