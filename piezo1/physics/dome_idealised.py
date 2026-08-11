"""Guo & MacKinnon's idealised dome — the geometry Figure 7 is drawn from.

Figure 7 of Guo & MacKinnon (eLife 2017;6:e33660) does not measure a dome. It
*idealises* one: a spherical mid-plane surface of radius 10.2 nm, centred
4.0 nm above the plane the membrane returns to, carrying a bilayer 3.6 nm
thick. Every number the paper's argument rests on then follows from those two
lengths by closed-form spherical-cap geometry:

=========================  ===============  ==============================
quantity                   paper            this module, from R and c
=========================  ===============  ==============================
dome opening diameter      "about 18 nm"    18.77 nm      (2a)
dome depth                 "about 6 nm"     6.20 nm       (R - c)
mid-plane surface area     400 nm^2         397.3 nm^2    (2 pi R h)
projected area             280 nm^2         276.6 nm^2    (pi a^2)
area released on           120 nm^2         120.7 nm^2    (difference)
  complete flattening
bending energy             "~150 k_BT"      152.7 k_BT    (2 kappa A / R^2)
=========================  ===============  ==============================

That table is the reason this module exists and is also its calibration: the
answers are known independently of any code here, so an implementation that
reproduces them is doing spherical-cap geometry correctly, and one that does
not is broken in a way a real structure would have hidden.

**What is idealised and what is measured.** The idealised dome is *not* a
measurement of 6B3R and must never be reported as one. It is a shape the
authors chose because it makes the energetics tractable, and they say so:
"the shape matches closely but not perfectly the hydrophobic boundaries of
Piezo, but it is sufficient for the following discussion". Our own measured
dome (:func:`piezo1.structure.geometry.measure_dome`) gives different numbers,
and :func:`compare_with_measured` reports both side by side — the gap between
an idealisation and a measurement is a result, not an error to be tuned away.

**Figure 7c made quantitative.** The paper's panel (c) is a schematic: two
planes and a caption saying the projected area grows as the dome flattens.
:func:`flattening_series` computes it. Flattening at constant membrane area is
a one-parameter family in the polar half-angle ``theta``, so the projected
area, the bending energy and the tension work are all functions of a single
coordinate, and the free-energy balance the paper writes as equation 3 can be
plotted rather than asserted.

Units: nm for lengths and areas throughout, k_BT for energies, k_BT/nm^2 for
tension — matching the paper. Angstrom appears only at the boundary with
:mod:`piezo1.structure.geometry`, which reports Angstrom, and is converted in
:func:`compare_with_measured`.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["IdealisedDome", "FlatteningPoint", "guo2017_dome",
           "flattening_series", "compare_with_measured",
           "spherical_cap_from_measurement", "PUBLISHED_FIGURE7"]


#: What the paper states, panel by panel, with the tolerance each is quoted to.
#: The paper rounds — "about 18 nm", "400 nm^2" — so the tolerance is the
#: rounding, not an error bar. Read by the tests and by the panel registry.
PUBLISHED_FIGURE7 = {
    "radius_nm": (10.2, 0.0, "Fig 7a: mid-plane semi-sphere radius"),
    "center_height_nm": (4.0, 0.0, "Fig 7a: centred 4.0 nm above the plane"),
    "thickness_nm": (3.6, 0.0, "Fig 7a: membrane 3.6 nm thick"),
    "opening_diameter_nm": (18.0, 1.0, "Results: 'diameter ... about 18 nm'"),
    "depth_nm": (6.0, 0.5, "Results: 'a depth of about 6 nm'"),
    "dome_area_nm2": (400.0, 5.0, "Fig 7-S1a: total mid-plane surface area"),
    "projected_area_nm2": (280.0, 5.0, "Fig 7-S1a: projected area A_proj"),
    "delta_area_nm2": (120.0, 5.0, "Discussion: 400 - 280 = 120 nm^2"),
    "bending_energy_kT": (150.0, 10.0, "Fig 7-S1b: 'approximately 150 k_BT'"),
    "stabilisation_kT": (42.0, 2.0,
                         "Discussion: 42 k_BT at one tenth lytic tension"),
    "lytic_tension_kT_per_nm2": (3.5, 0.1, "Discussion, after Rawicz 2000"),
}


# --------------------------------------------------------------------------
# The cap
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class IdealisedDome:
    """A spherical cap hanging below a plane, plus the bilayer it carries.

    The sphere has radius ``radius`` and its centre sits ``center_height``
    above the projection plane; the dome is the part of that sphere lying
    below the plane. With the paper's numbers the apex is 6.2 nm below the
    plane, which is the depth the micelle density shows.

    ``center_height`` must be smaller than ``radius`` or the sphere does not
    reach the plane and there is no cap; that is raised rather than clamped,
    because a silently truncated cap would report a plausible area.
    """

    radius: float                    # nm, mid-plane sphere radius R
    center_height: float             # nm, sphere centre above the plane, c
    thickness: float = field(        # nm, bilayer thickness
        default_factory=lambda: _P.value("membrane.thickness"))
    kappa: float = field(            # k_BT, bending modulus
        default_factory=lambda: _P.value("membrane.kappa"))

    def __post_init__(self) -> None:
        if self.radius <= 0:
            raise ValueError("dome radius must be positive")
        if not (0.0 <= self.center_height < self.radius):
            raise ValueError(
                f"centre height {self.center_height} nm must lie in "
                f"[0, {self.radius}) or the sphere never meets the plane")

    # ------------------------------------------------------------ geometry

    @property
    def rim_radius(self) -> float:
        """In-plane radius ``a`` where the cap meets the plane, nm."""
        return float(np.sqrt(self.radius ** 2 - self.center_height ** 2))

    @property
    def opening_diameter(self) -> float:
        """Diameter of the dome opening, nm. The paper's "about 18 nm"."""
        return 2.0 * self.rim_radius

    @property
    def depth(self) -> float:
        """Apex-to-plane depth ``R - c``, nm. The paper's "about 6 nm"."""
        return float(self.radius - self.center_height)

    @property
    def polar_angle(self) -> float:
        """Half-angle subtended at the sphere centre, radians.

        The single coordinate the flattening family is parameterised by:
        ``cos(theta) = c / R``, so ``theta`` runs from its closed-state value
        down to zero as the dome flattens.
        """
        return float(np.arccos(np.clip(self.center_height / self.radius,
                                       -1.0, 1.0)))

    @property
    def area(self) -> float:
        """Curved mid-plane area of the cap, nm^2. ``2 pi R h``."""
        return float(2.0 * np.pi * self.radius * self.depth)

    @property
    def projected_area(self) -> float:
        """Shadow of the cap on the plane, nm^2. ``pi a^2``."""
        return float(np.pi * self.rim_radius ** 2)

    @property
    def excess_area(self) -> float:
        """Area released if the cap flattened completely, nm^2.

        The paper's ``dA_proj`` of 120 nm^2 — an *upper bound* on the gating
        area change, reached only if PIEZO1 becomes fully co-planar with the
        membrane when it opens, which nothing establishes that it does.
        """
        return self.area - self.projected_area

    @property
    def contact_slope(self) -> float:
        """``dz/dr`` where the cap meets the plane, dimensionless ``a/c``.

        Compare ``membrane.piezo1_contact_slope`` (1.992, from the fitted 7WLT
        cap): the idealised dome is the steeper of the two, so the small-slope
        membrane theory is even less applicable to it.
        """
        if self.center_height <= 0:
            return float("inf")
        return float(self.rim_radius / self.center_height)

    @property
    def contact_angle_deg(self) -> float:
        return float(np.degrees(np.arctan(self.contact_slope)))

    def leaflet_radii(self) -> tuple[float, float]:
        """Radii of the two leaflet surfaces, nm — inner then outer.

        Concentric with the mid-plane sphere at half the bilayer thickness on
        each side. This is what Figure 7a draws as two translucent shells.
        """
        half = 0.5 * self.thickness
        return (self.radius - half, self.radius + half)

    # ------------------------------------------------------------- energy

    @property
    def bending_energy(self) -> float:
        """Helfrich bending energy of the cap, k_BT.

        For a sphere the mean curvature is uniform at ``H = 1/R``, so
        ``E = (kappa/2) * (2H)^2 * A = 2 * kappa * A / R^2``.

        This is the cap **only**. The paper is explicit that its ~150 k_BT
        excludes "the curvature back into the membrane plane, which will
        substantially increase the bending energy further" — that skirt is the
        Helfrich footprint, and :mod:`piezo1.physics.elastica` is what computes
        it. Spontaneous curvature is taken as zero, as the paper does, and the
        Gaussian term is omitted: by Gauss-Bonnet it contributes
        ``2 pi kappa_G (1 - cos theta)`` and depends on a modulus this project
        does not have a value for.
        """
        return float(2.0 * self.kappa * self.area / self.radius ** 2)

    def stabilisation(self, tension: float) -> float:
        """Open-state stabilisation ``gamma * dA_proj`` at a tension, k_BT.

        ``tension`` in k_BT/nm^2. At one tenth of the lytic tension the paper
        quotes 42 k_BT, which is this product at its 120 nm^2.
        """
        return float(tension * self.excess_area)

    # -------------------------------------------------------------- report

    def as_dict(self) -> dict:
        return {
            "radius_nm": self.radius,
            "center_height_nm": self.center_height,
            "thickness_nm": self.thickness,
            "rim_radius_nm": self.rim_radius,
            "opening_diameter_nm": self.opening_diameter,
            "depth_nm": self.depth,
            "polar_angle_deg": float(np.degrees(self.polar_angle)),
            "dome_area_nm2": self.area,
            "projected_area_nm2": self.projected_area,
            "delta_area_nm2": self.excess_area,
            "contact_slope": self.contact_slope,
            "contact_angle_deg": self.contact_angle_deg,
            "bending_energy_kT": self.bending_energy,
            "kappa_kT": self.kappa,
        }

    def summary(self) -> str:
        return (f"R = {self.radius:.1f} nm, opening {self.opening_diameter:.1f} nm, "
                f"depth {self.depth:.1f} nm | area {self.area:.0f} nm^2, "
                f"projected {self.projected_area:.0f} nm^2, "
                f"released {self.excess_area:.0f} nm^2 | "
                f"bending {self.bending_energy:.0f} kT")


def guo2017_dome() -> IdealisedDome:
    """The dome exactly as Figure 7a specifies it.

    Both lengths come from the parameter registry so that an override moves
    every derived number with them, and so the source of each is recorded.
    """
    return IdealisedDome(radius=_P.value("dome.published_radius_closed"),
                         center_height=_P.value("dome.idealised_center_height"))


# --------------------------------------------------------------------------
# Figure 7c: flattening, as arithmetic rather than a schematic
# --------------------------------------------------------------------------

@dataclass
class FlatteningPoint:
    """One point along the flattening coordinate."""

    polar_angle_deg: float
    radius_nm: float
    rim_radius_nm: float
    depth_nm: float
    projected_area_nm2: float
    #: Projected-area gain relative to the closed dome — the paper's dA_proj.
    delta_projected_nm2: float
    bending_energy_kT: float
    #: Bending-energy change relative to the closed dome. Negative: flattening
    #: a curved membrane always releases bending energy, which is the paper's
    #: statement that dG_bend is favourable for closed -> open.
    delta_bending_kT: float

    def free_energy(self, tension: float, delta_g_prot: float) -> float:
        """Equation 3 at this degree of flattening, k_BT.

        ``dG = (dG_prot + dG_bend) - gamma * dA_proj``, with both differences
        taken from the closed state. ``tension`` in k_BT/nm^2.
        """
        return float(delta_g_prot + self.delta_bending_kT
                     - tension * self.delta_projected_nm2)


def flattening_series(dome: IdealisedDome | None = None,
                      n: int = 60) -> list[FlatteningPoint]:
    """Flatten the dome at constant membrane area and report what changes.

    The constraint is the one the paper's Figure 7c draws: the dome does not
    gain or lose bilayer as it flattens, it *transfers* out-of-plane area into
    the plane. Holding the mid-plane area ``A`` fixed, a cap of polar half-angle
    ``theta`` has ``R = sqrt(A / (2 pi (1 - cos theta)))``, so the whole family
    is swept by ``theta`` alone, from the closed dome's value down towards
    zero.

    The last point is not exactly flat: at ``theta = 0`` the radius diverges,
    so the sweep stops just short and the residual is reported through
    ``delta_projected_nm2`` approaching the closed dome's ``excess_area``.
    """
    dome = dome or guo2017_dome()
    area = dome.area
    closed_theta = dome.polar_angle
    # Stop short of zero: R -> infinity there. A degree of arc leaves the
    # projected area within 0.01% of the fully flat value, so nothing visible
    # is lost and no division blows up.
    thetas = np.linspace(closed_theta, np.radians(1.0), n)

    out: list[FlatteningPoint] = []
    for theta in thetas:
        one_minus_cos = max(1.0 - np.cos(theta), 1e-12)
        radius = float(np.sqrt(area / (2.0 * np.pi * one_minus_cos)))
        rim = float(radius * np.sin(theta))
        depth = float(radius * one_minus_cos)
        projected = float(np.pi * rim ** 2)
        bending = float(2.0 * dome.kappa * area / radius ** 2)
        out.append(FlatteningPoint(
            polar_angle_deg=float(np.degrees(theta)),
            radius_nm=radius, rim_radius_nm=rim, depth_nm=depth,
            projected_area_nm2=projected,
            delta_projected_nm2=projected - dome.projected_area,
            bending_energy_kT=bending,
            delta_bending_kT=bending - dome.bending_energy))
    return out


# --------------------------------------------------------------------------
# Idealisation against measurement
# --------------------------------------------------------------------------

def spherical_cap_from_measurement(radius_nm: float, rim_radius_nm: float
                                   ) -> IdealisedDome:
    """The cap with this sphere radius that meets the plane at this rim.

    Lets a *measured* dome be expressed in the paper's own parameterisation, so
    the two can be compared in the same terms rather than through two different
    definitions of area. Raises if the rim is wider than the sphere, which
    means the fit and the footprint disagree about the same structure and the
    numbers must not be silently reconciled.
    """
    if rim_radius_nm >= radius_nm:
        raise ValueError(
            f"rim radius {rim_radius_nm:.2f} nm is not smaller than the sphere "
            f"radius {radius_nm:.2f} nm — no spherical cap has both")
    height = float(np.sqrt(radius_nm ** 2 - rim_radius_nm ** 2))
    return IdealisedDome(radius=radius_nm, center_height=height)


def compare_with_measured(dome_geometry) -> dict:
    """Guo & MacKinnon's idealisation beside our measurement of a structure.

    ``dome_geometry`` is a :class:`piezo1.structure.geometry.DomeGeometry`,
    which reports Angstrom; everything here is nm.

    The comparison is deliberately blunt. The idealised cap is a shape chosen
    for tractability and the measured one follows the transmembrane helices
    wherever they go, so the two areas are not measurements of the same thing:
    ours integrates the real radial profile out to the outermost helix, which
    on a full-length model reaches well past the point where the paper's cap
    has already met the plane. Both are reported with the ratio, and nothing
    here adjusts either to agree.
    """
    ideal = guo2017_dome()
    measured_radius = dome_geometry.radius_of_curvature / 10.0
    measured_rim = dome_geometry.footprint_radius / 10.0
    measured = {
        "radius_nm": measured_radius,
        "rim_radius_nm": measured_rim,
        "depth_nm": dome_geometry.dome_depth / 10.0,
        "dome_area_nm2": dome_geometry.dome_area / 100.0,
        "projected_area_nm2": dome_geometry.projected_area / 100.0,
        "delta_area_nm2": dome_geometry.excess_area / 100.0,
    }
    # Express the measurement as a cap where that is possible, so the bending
    # energy is computed the same way for both. A profile whose footprint
    # exceeds its own radius of curvature has no cap form; say so.
    try:
        as_cap = spherical_cap_from_measurement(measured_radius, measured_rim)
        measured["bending_energy_kT"] = as_cap.bending_energy
        measured["cap_form"] = as_cap.as_dict()
        cap_note = ""
    except ValueError as exc:
        measured["bending_energy_kT"] = None
        measured["cap_form"] = None
        cap_note = str(exc)

    ideal_d = ideal.as_dict()
    ratios = {}
    for key in ("radius_nm", "depth_nm", "dome_area_nm2",
                "projected_area_nm2", "delta_area_nm2"):
        denominator = ideal_d.get(key)
        if denominator:
            ratios[key] = measured[key] / denominator

    return {
        "idealised": ideal_d,
        "measured": measured,
        "measured_over_idealised": ratios,
        "cap_form_note": cap_note,
        "caveat": (
            "The idealised dome is a shape Guo & MacKinnon chose to make the "
            "energetics tractable, not a fit to 6B3R. Its areas and ours "
            "integrate different surfaces and are not expected to agree; the "
            "difference is a statement about idealisation, not an error."),
    }
