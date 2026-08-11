"""Membrane-dome geometry for PIEZO channels.

PIEZO1 does not sit flat in the bilayer: its three curved blades bend the
membrane into a dome that bulges towards the cytoplasm.  Flattening that dome
under tension increases the projected area of the protein-plus-membrane unit,
and it is that area change, multiplied by tension, that supplies the gating
energy (Guo & MacKinnon 2017; Haselwandter & MacKinnon 2018).

Everything needed to measure the dome from coordinates lives here:

* :func:`fit_sphere`          - algebraic sphere fit with an iterative refit.
* :func:`radial_profile`      - axial height as a function of distance from the
                                three-fold axis, the raw dome cross-section.
* :func:`measure_dome`        - the full :class:`DomeGeometry` summary.

Reference values to sanity-check against (closed/curved state):
    radius of curvature   ~ 100-120 Angstrom   (10.2 nm, Haselwandter 2018;
                                                11.8 nm, Vaisey 2026)
    dome depth            ~ 60-70 Angstrom     (Chong 2021; Shi 2020)
    mid-bilayer dome area ~ 390 nm^2           (Haselwandter 2018)
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

from .superpose import SymmetryAxis, detect_c3_axis

#: A chain shorter than this is a bound peptide rather than a protomer.
#: Matches the floor the report used inline before this was extracted.
MIN_CA_FOR_SURFACE = 300

__all__ = ["tm_surface_points", "tm_surface_by_chain", "fit_sphere", "SphereFit",
           "radial_profile", "RadialProfile", "DomeGeometry",
           "measure_dome", "spherical_cap_area", "sphere_points"]


def sphere_points(n: int) -> np.ndarray:
    """Golden-spiral points on a unit sphere — even, deterministic coverage.

    Lived in ``analysis.measure`` until Round 84, where it was fine while only
    the SASA used it. The screened-Coulomb surface in ``physics.electrostatics``
    needs the same point set — it is the same Shrake-Rupley construction, and
    two independent copies would be two surfaces that could silently differ —
    and ``physics`` importing ``analysis`` points the dependency arrow
    backwards. ``analysis.measure`` re-exports it, so the SASA is unchanged
    down to the bit.
    """
    i = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    return np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi),
                     np.cos(phi)], axis=1)


# --------------------------------------------------------------------------
# Sphere fitting
# --------------------------------------------------------------------------

@dataclass
class SphereFit:
    center: np.ndarray
    radius: float
    rmse: float
    n_points: int

    def signed_deviation(self, xyz: np.ndarray) -> np.ndarray:
        return np.linalg.norm(np.asarray(xyz) - self.center, axis=1) - self.radius


def fit_sphere(points: np.ndarray, iterations: int = 3,
               trim: float = 0.0) -> SphereFit:
    """Least-squares sphere through ``points``.

    Uses the standard linearisation ``|x|^2 = 2 x.c + (r^2 - |c|^2)``, which is
    exact for the algebraic residual and an excellent starting point for the
    geometric one.  ``trim`` optionally discards that outer fraction of the
    worst-fitting points on each refit, which stops flexible loops from
    dominating a fit that is meant to describe the transmembrane surface.
    """
    pts = np.asarray(points, dtype=np.float64)
    if len(pts) < 4:
        raise ValueError("need at least four points to fit a sphere")

    keep = np.ones(len(pts), dtype=bool)
    center = pts.mean(axis=0)
    radius = float(np.linalg.norm(pts - center, axis=1).mean())

    for _ in range(max(1, iterations)):
        p = pts[keep]
        a = np.hstack([2.0 * p, np.ones((len(p), 1))])
        b = (p * p).sum(axis=1)
        sol, *_ = np.linalg.lstsq(a, b, rcond=None)
        center = sol[:3]
        radius = float(np.sqrt(max(sol[3] + center @ center, 0.0)))
        resid = np.abs(np.linalg.norm(pts - center, axis=1) - radius)
        if trim > 0.0:
            cutoff = np.quantile(resid, 1.0 - trim)
            keep = resid <= cutoff
        else:
            keep = np.ones(len(pts), dtype=bool)

    resid = np.linalg.norm(pts[keep] - center, axis=1) - radius
    return SphereFit(center=center, radius=radius,
                     rmse=float(np.sqrt((resid ** 2).mean())), n_points=int(keep.sum()))


# --------------------------------------------------------------------------
# Radial height profile
# --------------------------------------------------------------------------

@dataclass
class RadialProfile:
    """Mean axial height in concentric annuli about the symmetry axis."""

    r: np.ndarray             # annulus centre, Angstrom
    z: np.ndarray             # mean height, Angstrom
    z_std: np.ndarray
    count: np.ndarray

    @property
    def valid(self) -> np.ndarray:
        return self.count > 0

    def depth(self) -> float:
        """Height difference between the outermost and innermost valid bins."""
        v = self.valid
        if v.sum() < 2:
            return 0.0
        return float(self.z[v][-1] - self.z[v][0])


def radial_profile(xyz: np.ndarray, axis: SymmetryAxis, n_bins: int = 40,
                   r_max: float | None = None) -> RadialProfile:
    """Height above the axis origin, binned by distance from the axis."""
    xyz = np.asarray(xyz, dtype=np.float64)
    r = axis.radial(xyz)
    z = axis.project(xyz)
    if r_max is None:
        r_max = float(np.quantile(r, 0.995))
    edges = np.linspace(0.0, r_max, n_bins + 1)
    idx = np.clip(np.digitize(r, edges) - 1, 0, n_bins - 1)
    inside = r <= r_max

    count = np.bincount(idx[inside], minlength=n_bins)
    zsum = np.bincount(idx[inside], weights=z[inside], minlength=n_bins)
    z2sum = np.bincount(idx[inside], weights=z[inside] ** 2, minlength=n_bins)
    with np.errstate(invalid="ignore", divide="ignore"):
        zmean = np.where(count > 0, zsum / np.maximum(count, 1), np.nan)
        zvar = np.where(count > 1, z2sum / np.maximum(count, 1) - zmean ** 2, 0.0)
    return RadialProfile(r=0.5 * (edges[:-1] + edges[1:]), z=zmean,
                         z_std=np.sqrt(np.maximum(zvar, 0.0)), count=count)


# --------------------------------------------------------------------------
# Dome summary
# --------------------------------------------------------------------------

def spherical_cap_area(radius: float, height: float) -> float:
    """Curved surface area of a spherical cap of given sphere radius and rise."""
    return float(2.0 * np.pi * radius * max(height, 0.0))


@dataclass
class DomeGeometry:
    """Quantitative description of the membrane dome formed by one channel."""

    axis: SymmetryAxis
    sphere: SphereFit
    profile: RadialProfile
    #: Radius of curvature of the mid-membrane surface, Angstrom.
    radius_of_curvature: float = 0.0
    #: Peak-to-trough height of the dome, Angstrom.
    dome_depth: float = 0.0
    #: In-plane radius out to which the protein deforms the bilayer, Angstrom.
    footprint_radius: float = 0.0
    #: Area of the curved mid-membrane surface, Angstrom^2.
    dome_area: float = 0.0
    #: Area of that surface projected onto the membrane plane, Angstrom^2.
    projected_area: float = 0.0
    n_atoms_used: int = 0
    notes: dict = field(default_factory=dict)

    @property
    def excess_area(self) -> float:
        """Area released to the bilayer if the dome flattened completely."""
        return self.dome_area - self.projected_area

    def summary(self) -> str:
        return (
            f"radius of curvature {self.radius_of_curvature / 10:.1f} nm | "
            f"dome depth {self.dome_depth / 10:.1f} nm | "
            f"footprint radius {self.footprint_radius / 10:.1f} nm | "
            f"dome area {self.dome_area / 100:.0f} nm^2 | "
            f"projected {self.projected_area / 100:.0f} nm^2 | "
            f"excess {self.excess_area / 100:.0f} nm^2"
        )


def measure_dome(xyz_by_chain: list[np.ndarray], surface_xyz: np.ndarray,
                 n_bins: int | None = None, trim: float | None = None) -> DomeGeometry:
    """Measure dome geometry from transmembrane-surface coordinates.

    Parameters
    ----------
    xyz_by_chain:
        Three equivalent coordinate sets, one per protomer, used only to
        recover the three-fold axis.
    surface_xyz:
        The points that trace the mid-membrane surface — in practice the
        C-alpha atoms of the transmembrane helices of all three protomers.
    trim:
        Fraction of worst-fitting points ignored when fitting the sphere.
    """
    if n_bins is None:
        n_bins = int(_P.value("geometry.radial_bins"))
    if trim is None:
        trim = _P.value("geometry.sphere_trim")
    axis = detect_c3_axis(xyz_by_chain)
    surface = np.asarray(surface_xyz, dtype=np.float64)

    sphere = fit_sphere(surface, iterations=4, trim=trim)

    # Point the axis from the sphere centre towards the surface centroid so
    # that "up" is out of the dome regardless of chain ordering.
    outward = surface.mean(axis=0) - sphere.center
    if np.dot(outward, axis.direction) < 0:
        axis = SymmetryAxis(point=axis.point, direction=-axis.direction,
                            order=axis.order, rmsd=axis.rmsd, angle_deg=axis.angle_deg)

    # Re-seat the axis origin at the dome apex so heights are measured from it.
    apex = sphere.center + axis.direction * sphere.radius
    axis = SymmetryAxis(point=apex, direction=axis.direction, order=axis.order,
                        rmsd=axis.rmsd, angle_deg=axis.angle_deg)

    prof = radial_profile(surface, axis, n_bins=n_bins)
    valid = prof.valid & np.isfinite(prof.z)
    r_valid = prof.r[valid]
    z_valid = prof.z[valid]

    footprint_r = float(r_valid.max()) if len(r_valid) else 0.0
    depth = float(np.nanmax(z_valid) - np.nanmin(z_valid)) if len(z_valid) else 0.0

    # Curved area by revolving the measured profile: dA = 2*pi*r*ds.
    if len(r_valid) > 1:
        dr = np.diff(r_valid)
        dz = np.diff(z_valid)
        ds = np.sqrt(dr ** 2 + dz ** 2)
        r_mid = 0.5 * (r_valid[1:] + r_valid[:-1])
        dome_area = float(2.0 * np.pi * np.sum(r_mid * ds))
    else:
        dome_area = 0.0
    projected_area = float(np.pi * footprint_r ** 2)

    return DomeGeometry(
        axis=axis, sphere=sphere, profile=prof,
        radius_of_curvature=sphere.radius, dome_depth=depth,
        footprint_radius=footprint_r, dome_area=dome_area,
        projected_area=projected_area, n_atoms_used=len(surface),
        notes={"sphere_rmse": sphere.rmse, "c3_rmsd": axis.rmsd,
               "c3_angle_deg": axis.angle_deg},
    )


def tm_surface_points(structure, reference: str, keep=None
                      ) -> tuple[np.ndarray, set]:
    """Mid-membrane surface points, one per transmembrane helix per protomer.

    The single definition of what :func:`measure_dome` is given. It lived in
    two places until Round 83 — inline in the report and again in the claims
    registry — which was survivable while only PIEZO1 was measured and stopped
    being so the moment a second protein had to be compared against the same
    surface: two definitions is exactly how a comparison ends up measuring how
    each side was defined.

    ``reference`` names the committed UniProt resource whose transmembrane
    features are used, so a human structure must not be read with mouse
    numbering and a PIEZO2 structure must not be read with either.
    ``keep`` restricts to a set of 1-based helix indices, which is how two
    entries resolving different amounts of blade are coverage-matched.

    Returns the points and the set of helix indices that contributed from all
    three protomers, so a caller can see what it actually measured.
    """
    by_chain, resolved = tm_surface_by_chain(structure, reference, keep)
    points = [p for chain_points in by_chain.values()
              for p in chain_points.values()]
    return np.array(points), resolved


def tm_surface_by_chain(structure, reference: str, keep=None
                        ) -> tuple[dict[str, dict[int, np.ndarray]], set]:
    """The same surface points, kept separate per protomer and per helix.

    :func:`tm_surface_points` pools them, which is right for fitting one dome
    to a trimer and wrong for any question about a *single* protomer — and
    Guo & MacKinnon 2017 asks exactly that one in Figure 4a: a subunit taken
    out of the trimer sits in a plane, the trimer does not. Answering it needs
    the same surface definition applied one chain at a time, so both callers
    share this loop rather than each writing down what "the membrane surface"
    means.

    Returns ``{chain: {helix_index: point}}`` and the set of helix indices
    present in all three protomers.
    """
    import json

    from ..config import RESOURCE_DIR

    helices = json.loads(
        (RESOURCE_DIR / f"uniprot_{reference}.json").read_text())["transmembrane"]
    out: dict[str, dict[int, np.ndarray]] = {}
    contributed: dict[int, int] = {}
    for chain in structure.chains:
        mask = structure.mask_ca() & (structure.chain == chain)
        if mask.sum() < MIN_CA_FOR_SURFACE:
            continue
        xyz, seq = structure.xyz[mask], structure.res_seq[mask]
        found: dict[int, np.ndarray] = {}
        for index, helix in enumerate(helices, start=1):
            if keep is not None and index not in keep:
                continue
            middle = 0.5 * (helix["start"] + helix["end"])
            half = max(2.0, (helix["end"] - helix["start"]) / 6.0)
            selected = (seq >= middle - half) & (seq <= middle + half)
            if selected.sum() >= 3:
                found[index] = xyz[selected].mean(axis=0)
                contributed[index] = contributed.get(index, 0) + 1
        if found:
            out[chain] = found
    resolved = {k for k, v in contributed.items() if v >= 3}
    return out, resolved
