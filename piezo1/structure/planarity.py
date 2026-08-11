"""How badly PIEZO1 fits a flat membrane — Guo & MacKinnon 2017, Figure 4a.

Figure 4a draws one protomer with two grey lines across it: "approximate
locations of planar membrane interfaces". The trimer in Figure 4b gets no such
lines, and the text says why:

    "A single subunit of Piezo removed from the trimer can be positioned
    reasonably well into the plane of a lipid membrane. However, the detergent
    micelle containing a trimer is curved into a dome shape."

and later, as the pivot of the whole discussion:

    "We think it is significant that a single subunit is compatible with a
    planar membrane, whereas the trimer is not."

That is a claim about *fit residuals*, drawn as a picture. This module measures
it. Fit a plane to one protomer's transmembrane band, then to all three, and
compare the residuals — a real comparison, with the same surface definition on
both sides, taken from :func:`piezo1.structure.geometry.tm_surface_by_chain`
so that "the membrane surface" is defined once for the dome and for this.

Two more numbers the paper states in words and this recovers:

* **The beam sits at about 60 degrees to the pore axis**, "instead of 90 as we
  would expect if the trimer were located in a fully flattened membrane".
* **The arms project about 30 degrees out of the plane defined by the pore** —
  the same statement, complemented.

**What makes the comparison mean something.** "The trimer fits worse" is not
by itself a finding, so the residual is decomposed into the two things that
could produce it:

* **within-protomer** — how far one protomer's own band departs from its own
  best plane. This is the number Figure 4a's grey lines are drawn against.
* **arrangement** — what is left when each protomer is *made* exactly planar
  by projecting it onto its own plane and the three are left where the C3
  symmetry puts them. If the arms did not tilt out of the pore's plane this
  term would be zero by construction, so it is a measurement of the tilt and
  not of least squares having more points to disappoint.

The two add in quadrature to the trimer residual, and the residual of that
identity is reported: a decomposition that does not close has missed a third
contribution and should not be believed.

Point count is deliberately *not* controlled for, because it needs no control:
replicating a point set leaves a least-squares plane and its RMSD exactly
unchanged. An earlier version of this module reported that replication as a
"control", which was a tautology dressed as evidence — it agreed with the
protomer RMSD to every digit on every structure, because it could not do
anything else.

**Coverage is the thing that will fool you here, and it did.** Measured on
whatever each entry resolves, 6B3R's arrangement term is 17.2 A and Saotome's
6BPZ is 4.7 A — which reads as two structures of the same protein disagreeing
about whether it is curved. They do not. 6BPZ resolves 14 transmembrane
helices and 6B3R resolves 26, and restricting both to the 14 they share brings
them to 4.5 and 4.7 A. The non-planarity is carried almost entirely by the
distal blade, and an entry that does not resolve the blade cannot see it.

So :func:`blade_dependence` is the function to reach for when comparing
entries, and :func:`planarity` records which helices it used in
``meta["helices_used"]`` so that two results are never compared without
checking. This is the same trap ``analysis/paralogue.py`` was written after,
in a different place.

Angstrom throughout, matching the coordinates. Angles in degrees.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

from .geometry import tm_surface_by_chain
from .superpose import SymmetryAxis, detect_c3_axis

__all__ = ["PlaneFit", "PlanarityComparison", "fit_plane", "planarity",
           "beam_angle", "BeamGeometry", "blade_dependence", "BladeDependence",
           "PROXIMAL_FIRST_HELIX"]

#: First transmembrane helix of the pore-proximal module. TM25 opens THU7, the
#: most distal unit every deposited PIEZO1 entry in this project's catalogue
#: resolves in all three protomers; TM13-24 (THU4-6) are resolved by some
#: entries and not others. Splitting here is what makes two entries comparable.
PROXIMAL_FIRST_HELIX = int(_P.value("architecture.proximal_first_helix"))


# --------------------------------------------------------------------------
# Plane fitting
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class PlaneFit:
    """A least-squares plane through a point set."""

    point: np.ndarray            # a point on the plane (the centroid)
    normal: np.ndarray           # unit normal
    rmsd: float                  # root-mean-square out-of-plane deviation, A
    max_deviation: float         # worst single point, A
    n_points: int

    def deviation(self, xyz: np.ndarray) -> np.ndarray:
        """Signed distance of each point from the plane, Angstrom."""
        return (np.asarray(xyz, dtype=np.float64) - self.point) @ self.normal

    def thickness_needed(self) -> float:
        """Slab thickness that would contain every point, Angstrom.

        The honest version of "can this sit in a membrane": not the RMS, but
        the full span, since a bilayer has to contain all of it at once.
        """
        return float(self.max_deviation * 2.0)


def fit_plane(points: np.ndarray) -> PlaneFit:
    """Total-least-squares plane: the normal is the smallest principal axis.

    Uses the SVD of the centred coordinates rather than a regression of z on
    (x, y), which would have no answer for a plane containing the z axis and a
    badly conditioned one for a plane near it. PIEZO1's blades approach both.
    """
    pts = np.asarray(points, dtype=np.float64)
    if len(pts) < 3:
        raise ValueError("need at least three points to fit a plane")
    centroid = pts.mean(axis=0)
    _, _, vh = np.linalg.svd(pts - centroid, full_matrices=False)
    normal = vh[-1]
    normal = normal / np.linalg.norm(normal)
    dev = (pts - centroid) @ normal
    return PlaneFit(point=centroid, normal=normal,
                    rmsd=float(np.sqrt((dev ** 2).mean())),
                    max_deviation=float(np.abs(dev).max()), n_points=len(pts))


# --------------------------------------------------------------------------
# One protomer against three
# --------------------------------------------------------------------------

@dataclass
class PlanarityComparison:
    """Figure 4a's claim, as residuals."""

    #: Plane fitted to each protomer's transmembrane band, by chain.
    per_protomer: dict[str, PlaneFit]
    #: Plane fitted to all three at once.
    trimer: PlaneFit
    #: Mean of the per-protomer RMSDs — the within-protomer term.
    protomer_rmsd: float
    #: Residual left when each protomer is flattened onto its own plane and
    #: the three are left where the symmetry puts them. Zero if the arms lie
    #: in the pore's plane; this is the term that carries the dome.
    arrangement_rmsd: float
    #: ``sqrt(within^2 + arrangement^2) - trimer``, Angstrom. Near zero means
    #: the two terms account for the whole residual.
    decomposition_residual: float
    #: How much worse the trimer is than one protomer.
    ratio: float
    #: Angle between each protomer's plane normal and the C3 axis, degrees —
    #: the tilt the arms carry out of the pore's plane.
    protomer_tilt_deg: dict[str, float]
    axis: SymmetryAxis
    n_helices: int
    meta: dict = field(default_factory=dict)

    @property
    def mean_tilt_deg(self) -> float:
        values = list(self.protomer_tilt_deg.values())
        return float(np.mean(values)) if values else float("nan")

    @property
    def supports_paper(self) -> bool:
        """True when the trimer's excess non-planarity is the *arrangement*.

        Guo & MacKinnon's claim is specifically that a subunit is compatible
        with a planar membrane and the assembly is not — so what has to
        dominate is the term that appears only on assembly. Requiring the
        arrangement term to exceed the within-protomer one makes the test
        capable of failing, and on the flattened structure (7WLU) it does.
        """
        return self.arrangement_rmsd > self.protomer_rmsd

    def summary(self) -> str:
        return (f"protomer plane RMSD {self.protomer_rmsd:.1f} A | "
                f"arrangement {self.arrangement_rmsd:.1f} A | "
                f"trimer {self.trimer.rmsd:.1f} A, {self.ratio:.1f}x worse | "
                f"arms tilt {self.mean_tilt_deg:.0f} deg out of the pore plane")


def planarity(structure, reference: str, keep=None) -> PlanarityComparison:
    """Measure Figure 4a: a protomer in a plane, the trimer not.

    ``reference`` names the committed UniProt resource whose transmembrane
    features define the band — the same argument
    :func:`piezo1.structure.geometry.tm_surface_points` takes, and for the same
    reason: reading a human structure with mouse numbering would move the
    band, not fail.
    """
    by_chain, resolved = tm_surface_by_chain(structure, reference, keep)
    if len(by_chain) < 3:
        raise ValueError(
            f"need three well-resolved protomers, found {len(by_chain)}")
    if not resolved:
        raise ValueError("no transmembrane helix is present in all three "
                         "protomers; nothing is comparable")

    # Only helices present in every protomer, so the three fits are of the
    # same thing and the trimer fit is not dominated by whichever chain
    # resolves the most blade.
    per_protomer, chain_points = {}, {}
    for chain, helices in by_chain.items():
        pts = np.array([helices[i] for i in sorted(resolved) if i in helices])
        if len(pts) < 3:
            continue
        chain_points[chain] = pts
        per_protomer[chain] = fit_plane(pts)
    if len(per_protomer) < 3:
        raise ValueError("fewer than three protomers share enough helices")

    pooled = np.vstack(list(chain_points.values()))
    trimer = fit_plane(pooled)

    axis = detect_c3_axis([np.asarray(p) for p in chain_points.values()])

    # Arrangement term: make each protomer exactly planar by projecting it
    # onto its own best-fit plane, leave the three where they are, refit. What
    # survives is non-planarity that only exists because three tilted copies
    # are arranged around an axis — which is the thing Figure 4a is about.
    flattened = []
    for chain, pts in chain_points.items():
        fit = per_protomer[chain]
        flattened.append(pts - np.outer(fit.deviation(pts), fit.normal))
    arrangement = fit_plane(np.vstack(flattened))

    protomer_rmsd = float(np.mean([f.rmsd for f in per_protomer.values()]))
    predicted = float(np.hypot(protomer_rmsd, arrangement.rmsd))
    tilts = {}
    for chain, fit in per_protomer.items():
        cosine = abs(float(np.dot(fit.normal, axis.direction)))
        tilts[chain] = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    return PlanarityComparison(
        per_protomer=per_protomer, trimer=trimer,
        protomer_rmsd=protomer_rmsd, arrangement_rmsd=arrangement.rmsd,
        decomposition_residual=predicted - trimer.rmsd,
        ratio=trimer.rmsd / protomer_rmsd if protomer_rmsd > 0 else float("inf"),
        protomer_tilt_deg=tilts, axis=axis, n_helices=len(resolved),
        meta={"helices_used": sorted(resolved),
              "reference": reference,
              "bilayer_thickness_A": _P.value("membrane.thickness") * 10.0,
              "caveat": ("A plane through the transmembrane band is the "
                         "loosest possible test of planarity: it asks whether "
                         "any flat slab contains the helix mid-points, not "
                         "whether a bilayer could adopt that shape.")})


# --------------------------------------------------------------------------
# Where the non-planarity actually lives
# --------------------------------------------------------------------------

@dataclass
class BladeDependence:
    """How much of the trimer's non-planarity the distal blade carries.

    Computed on **one** structure, which is the point: it needs no second
    entry and so cannot be confounded by what a second entry resolves.
    """

    full: PlanarityComparison
    #: Same measurement restricted to helices from ``split_at`` onwards.
    proximal: PlanarityComparison
    split_at: int
    #: Helices dropped by the restriction, i.e. the distal blade this entry
    #: resolves. Empty means the entry resolves nothing distal and the two
    #: measurements are the same one reported twice.
    distal_helices: tuple[int, ...]

    @property
    def blade_share(self) -> float:
        """Fraction of the arrangement term lost when the blade is dropped.

        NaN when the entry resolves no distal helices — the honest answer to
        "how much does the blade contribute" for a structure that does not
        have one is *unmeasurable*, not zero.
        """
        if not self.distal_helices or self.full.arrangement_rmsd <= 0:
            return float("nan")
        return float(1.0 - self.proximal.arrangement_rmsd
                     / self.full.arrangement_rmsd)

    def summary(self) -> str:
        if not self.distal_helices:
            return (f"no helix distal to TM{self.split_at} resolved in all "
                    f"three protomers — the blade's contribution cannot be "
                    f"measured on this entry")
        return (f"arrangement {self.full.arrangement_rmsd:.1f} A over "
                f"{self.full.n_helices} helices falls to "
                f"{self.proximal.arrangement_rmsd:.1f} A over "
                f"{self.proximal.n_helices} once the {len(self.distal_helices)} "
                f"distal ones are dropped — the blade carries "
                f"{100 * self.blade_share:.0f}% of it")


def blade_dependence(structure, reference: str,
                     split_at: int | None = None) -> BladeDependence:
    """Split the non-planarity into a blade term and a pore-module term.

    Answers the question two entries with different coverage cannot be asked
    directly: is the trimer non-planar because of its arms, or because of the
    part every structure resolves?

    For 6B3R the answer is the arms, by a wide margin — which is what makes
    Guo & MacKinnon's Figure 4a a statement about the blade, and what makes
    comparing a blade-resolving entry with one that stops at TM25 meaningless.
    """
    if split_at is None:
        split_at = int(_P.value("architecture.proximal_first_helix"))
    _, resolved = tm_surface_by_chain(structure, reference)
    proximal_set = {i for i in resolved if i >= split_at}
    distal = tuple(sorted(i for i in resolved if i < split_at))
    if len(proximal_set) < 3:
        raise ValueError(
            f"only {len(proximal_set)} helices at or beyond TM{split_at} are "
            f"resolved in all three protomers; nothing to compare against")
    return BladeDependence(
        full=planarity(structure, reference),
        proximal=planarity(structure, reference, keep=proximal_set),
        split_at=split_at, distal_helices=distal)


# --------------------------------------------------------------------------
# The beam angle the paper states
# --------------------------------------------------------------------------

def _long_helix_axis(ca: np.ndarray) -> np.ndarray:
    """Axis of a helix from its C-alpha trace, N to C.

    The first principal component. That is the same estimate
    ``analysis.measure.helix_axis`` makes, written out here rather than
    imported because ``structure`` importing ``analysis`` points the project's
    dependency arrow backwards — the guard in ``test_architecture.py`` caught
    exactly that when this function was first written as an import.

    Duplicating an estimate is normally the wrong trade. It is acceptable here
    because the beam is 66 residues, roughly eighteen turns, where the
    principal axis and the true helix axis agree to well under a degree; the
    short-window case, where they do not, has its own unbiased estimator in
    ``structure.architecture.helical_windows``.
    """
    x = np.asarray(ca, dtype=np.float64)
    if len(x) < 4:
        raise ValueError("need at least four C-alpha positions")
    _, _, vh = np.linalg.svd(x - x.mean(axis=0), full_matrices=False)
    axis = vh[0]
    if np.dot(x[-1] - x[0], axis) < 0:
        axis = -axis
    return axis / np.linalg.norm(axis)


@dataclass
class BeamGeometry:
    """Angle between the beam helix and the central pore axis."""

    #: One entry per protomer, degrees from the C3 axis.
    angle_deg: dict[str, float]
    #: Residue range used, in the structure's own numbering.
    residue_range: tuple[int, int]
    numbering: str

    @property
    def mean_deg(self) -> float:
        values = list(self.angle_deg.values())
        return float(np.mean(values)) if values else float("nan")

    @property
    def out_of_plane_deg(self) -> float:
        """90 minus the beam angle — the arms' projection out of the plane.

        The paper gives both forms: the beam "about 60 degrees instead of 90",
        and the arms projecting "approximately 30 degrees out of the plane
        defined by the pore".
        """
        return 90.0 - self.mean_deg


def beam_angle(structure, axis: SymmetryAxis | None = None,
               residue_range: tuple[int, int] = (1300, 1365),
               numbering: str = "mouse") -> BeamGeometry:
    """Angle between each protomer's beam helix and the pore axis.

    The default range is the beam as Guo & MacKinnon define it, in **mouse**
    numbering (residues 1300-1365; the same helix is human 1305-1370, and
    ``domains.json`` carries both). Passing a range in the wrong numbering
    would select a real but different set of residues and return a plausible
    angle, so the numbering is named in the result rather than assumed.

    The axis is the helix's own principal axis, not the vector between its
    ends: the beam is long and slightly bent, and end-to-end would report the
    chord.
    """
    lo, hi = residue_range
    ca = structure.mask_ca()
    angles: dict[str, float] = {}
    blocks = []
    for chain in structure.chains:
        mask = ca & (structure.chain == chain) & \
            (structure.res_seq >= lo) & (structure.res_seq <= hi)
        if mask.sum() >= 8:
            blocks.append(structure.xyz[mask])
    if axis is None:
        if len(blocks) < 3:
            raise ValueError(
                f"only {len(blocks)} of the chains have 8 or more C-alpha "
                f"atoms in {lo}-{hi}, so the three-fold axis cannot be "
                f"recovered from the beams. The usual cause is a residue "
                f"range in the wrong numbering system: these numbers are "
                f"{numbering}, and a structure deposited in another one has "
                f"real residues there that are not the beam.")
        axis = detect_c3_axis(blocks)

    for chain in structure.chains:
        mask = ca & (structure.chain == chain) & \
            (structure.res_seq >= lo) & (structure.res_seq <= hi)
        if mask.sum() < 8:
            continue
        direction = _long_helix_axis(structure.xyz[mask])
        cosine = abs(float(np.dot(direction, axis.direction)))
        angles[chain] = float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))

    if not angles:
        raise ValueError(
            f"no chain has at least 8 C-alpha atoms in {lo}-{hi} "
            f"({numbering} numbering) — is this the right numbering system?")
    return BeamGeometry(angle_deg=angles, residue_range=(lo, hi),
                        numbering=numbering)
