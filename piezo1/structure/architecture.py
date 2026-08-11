"""Helical elements measured from coordinates — the cross-helices of Figure 7b.

Guo & MacKinnon 2017 colour three things in Figure 7b: the protein grey, the
beam red, and the **cross-helices** yellow. The beam has a residue range in the
paper (mouse 1300-1365, and ``domains.json`` carries it). The cross-helices do
not — anywhere. All the paper says is:

    "linkers between 4-TM units contain at least one helix that runs
    perpendicular to the TM helices and to the extended arms. These 'cross'
    helices are mostly hydrophobic and located inside the micelle density, near
    the intracellular interface."

So reproducing that panel means finding them by the property they are named
for: helical, in a linker between two 4-TM units, and running across rather
than through the membrane. That is what :func:`cross_helices` does.

**Why the helix test is geometric.** The obvious route — ask the renderer's
``assign_secondary_structure`` — would have ``structure`` importing from
``render``, against the one-way dependency this project keeps. The test here is
self-contained and, more usefully, calibratable: an alpha helix has a C-alpha
rise of 1.5 A per residue along its own axis, its C-alphas sit 2.3 A off that
axis, and it turns 100 degrees per residue right-handedly. All three are
checked against analytically generated helices before anything real is
measured, and each of a 3-10 helix, a pi helix, a beta strand, a left-handed
helix and a random coil fails at least one — which is what makes it a test
rather than a filter.

**What the threshold is doing.** "Perpendicular" is a word; 55 degrees is a
number, and it is a choice. :func:`cross_helix_scan` reports the count against
the threshold so the separation can be seen rather than trusted — on 6B3R the
transmembrane helices sit well below it and the linker helices well above, and
if that ever stopped being true the count would drift smoothly rather than the
picture quietly changing.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

from .planarity import fit_plane
from .superpose import SymmetryAxis, detect_c3_axis

__all__ = ["HelicalSegment", "helical_windows", "helical_segments",
           "cross_helices", "cross_helix_scan", "IDEAL_HELIX_RISE",
           "IDEAL_HELIX_RADIUS", "IDEAL_HELIX_TURN", "ideal_helix"]

#: The ideal alpha helix, read from the registry so that the geometry the
#: detector is calibrated against and the geometry it tests for cannot drift
#: apart. Module-level names are kept because the tests and the calibration
#: refer to them; ``_criteria()`` re-reads the registry at call time, so an
#: override takes effect on the next call rather than at import.
IDEAL_HELIX_RISE = _P.value("helix.rise")
IDEAL_HELIX_RADIUS = _P.value("helix.radius")
IDEAL_HELIX_TURN = _P.value("helix.turn")
_RISE_TOLERANCE = _P.value("helix.rise_tolerance")
_RADIUS_TOLERANCE = _P.value("helix.radius_tolerance")
_TURN_TOLERANCE = _P.value("helix.turn_tolerance")


def _criteria() -> tuple[float, float, float, float, float, float]:
    """Rise, radius and turn with their tolerances, resolved now."""
    return (_P.value("helix.rise"), _P.value("helix.rise_tolerance"),
            _P.value("helix.radius"), _P.value("helix.radius_tolerance"),
            _P.value("helix.turn"), _P.value("helix.turn_tolerance"))

#: The turn criterion is what excludes a random coil, and it is doing real
#: work: rise and radius alone passed 41% of the windows of a synthetic
#: 3.8-Angstrom-step random walk, because a walk with a fixed step length looks
#: locally like a helix on both. A coil has no consistent sense of rotation
#: about any axis, so requiring the turn to be right-handed and near 100
#: degrees at *every* step of the window removes it without touching the real
#: helices.


def ideal_helix(n: int, rise: float | None = None,
                radius: float | None = None,
                turn: float | None = None) -> np.ndarray:
    """C-alpha trace of an ideal alpha helix along +z, for calibration.

    100 degrees per residue and 1.5 A rise are the textbook values; the
    resulting 3.6 residues per turn and 5.4 A pitch are what the detector must
    recognise.
    """
    rise = _P.value("helix.rise") if rise is None else rise
    radius = _P.value("helix.radius") if radius is None else radius
    turn = _P.value("helix.turn") if turn is None else turn
    angle = np.radians(turn) * np.arange(n)
    return np.column_stack([radius * np.cos(angle), radius * np.sin(angle),
                            rise * np.arange(n)])


@dataclass
class HelicalSegment:
    """A run of C-alpha atoms behaving like a helix."""

    chain: str
    start: int                 # first residue number
    end: int                   # last residue number
    axis: np.ndarray           # unit vector, N-to-C
    center: np.ndarray
    rise: float                # A per residue along the axis
    radius: float              # mean C-alpha offset from the axis
    length: float              # end-to-end along the axis, A

    @property
    def n_residues(self) -> int:
        return self.end - self.start + 1

    def tilt_to(self, normal: np.ndarray) -> float:
        """Angle to a reference direction, degrees, folded into 0-90.

        Folded because a helix has no preferred N-to-C sense for this purpose:
        one running "down" across the membrane is as much a cross-helix as one
        running "up".
        """
        cosine = abs(float(np.dot(self.axis, normal)
                           / (np.linalg.norm(normal) or 1.0)))
        return float(np.degrees(np.arccos(np.clip(cosine, -1.0, 1.0))))


def helical_windows(ca: np.ndarray, window: int | None = None
                    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Per-window rise, radius, axis and minimum turn for a C-alpha trace.

    The window slides one residue at a time and each is judged on its own, so a
    helix that bends is still recognised throughout — which matters for the
    beam, 66 residues long and visibly curved.

    The axis comes from the centroids of consecutive 4-residue sub-windows,
    not from the principal axis of the window itself. Four consecutive
    C-alphas of an alpha helix span 400 degrees of turn, so their centroid sits
    essentially *on* the helix axis, and a line through those centroids is the
    axis. Taking the window's own principal axis instead tilts it towards the
    chord and biases the estimator: on an ideal helix it returns a rise of 1.63
    A and a radius of 2.07 A rather than the textbook 1.50 and 2.30, which
    would have forced the tolerances to absorb a systematic error and left them
    unable to say what they were really excluding.
    """
    window = int(_P.value("helix.window")) if window is None else int(window)
    ca = np.asarray(ca, dtype=np.float64)
    n = len(ca)
    if n < window:
        return np.zeros(0), np.zeros(0), np.zeros((0, 3)), np.zeros(0)
    sub = 4
    rises, radii, axes, turns = [], [], [], []
    for i in range(n - window + 1):
        block = ca[i:i + window]
        if window >= sub + 1:
            centroids = np.array([block[k:k + sub].mean(axis=0)
                                  for k in range(window - sub + 1)])
            base = centroids - centroids.mean(axis=0)
            _, _, vh = np.linalg.svd(base, full_matrices=False)
            axis = vh[0]
            if np.dot(centroids[-1] - centroids[0], axis) < 0:
                axis = -axis
        else:
            axis = block[-1] - block[0]
            axis = axis / (np.linalg.norm(axis) or 1.0)
        centred = block - block.mean(axis=0)
        along = centred @ axis
        perpendicular = centred - np.outer(along, axis)
        rises.append(float((along[-1] - along[0]) / (window - 1)))
        radii.append(float(np.linalg.norm(perpendicular, axis=1).mean()))
        axes.append(axis)

        # Signed rotation about the axis between consecutive residues. The
        # worst step in the window is kept, not the mean: a helix must turn
        # the same way every time, and averaging lets one reversal hide.
        step = np.degrees(np.arctan2(
            np.einsum("ij,j->i", np.cross(perpendicular[:-1], perpendicular[1:]),
                      axis),
            np.einsum("ij,ij->i", perpendicular[:-1], perpendicular[1:])))
        turns.append(float(np.abs(step - _P.value("helix.turn")).max())
                     if len(step) else 180.0)
    return np.array(rises), np.array(radii), np.array(axes), np.array(turns)


def helical_segments(structure, chain: str, lo: int, hi: int,
                     window: int | None = None, min_length: int | None = None
                     ) -> list[HelicalSegment]:
    """Helical runs within a residue range of one chain.

    Contiguity is enforced on **residue numbers**, not on array order: a
    cryo-EM model has unmodelled loops, and a window spanning one would join
    two helices across a gap and report a straight rod with a plausible rise.
    """
    if min_length is None:
        min_length = int(_P.value("cross_helix.min_length"))
    window = int(_P.value("helix.window")) if window is None else int(window)
    rise_target, rise_tol, radius_target, radius_tol, _, turn_tol = _criteria()
    mask = (structure.mask_ca() & (structure.chain == chain)
            & (structure.res_seq >= lo) & (structure.res_seq <= hi))
    if mask.sum() < window:
        return []
    order = np.argsort(structure.res_seq[mask])
    seq = structure.res_seq[mask][order]
    xyz = structure.xyz[mask][order].astype(np.float64)

    # Split at any break in the residue numbering.
    breaks = np.flatnonzero(np.diff(seq) != 1) + 1
    out: list[HelicalSegment] = []
    runs: list = []
    for block_seq, block_xyz in zip(np.split(seq, breaks), np.split(xyz, breaks)):
        if len(block_seq) < window:
            continue
        rise, radius, axes, turn = helical_windows(block_xyz, window)
        good = ((np.abs(rise - rise_target) <= rise_tol)
                & (np.abs(radius - radius_target) <= radius_tol)
                & (turn <= turn_tol))
        start = None
        for i, flag in enumerate(list(good) + [False]):
            if flag and start is None:
                start = i
            elif not flag and start is not None:
                first, last = start, i - 1 + window - 1
                if last - first + 1 >= min_length:
                    runs.append((block_seq, block_xyz, first, min(last, len(block_seq) - 1)))
                start = None

    # Merge runs that overlap or abut. A helix with a bend in the middle fails
    # the turn test for a window or two and comes back as two runs whose
    # residue spans overlap — the beam did exactly that, arriving as 1299-1317
    # and 1315-1365. Reporting both would double-count it and put a spurious
    # short helix in every figure.
    for block_seq, block_xyz, first, last in _merge_runs(runs):
        span = block_xyz[first:last + 1]
        centroids = np.array([span[k:k + 4].mean(axis=0)
                              for k in range(max(1, len(span) - 3))])
        base = centroids - centroids.mean(axis=0)
        if len(centroids) >= 2:
            _, _, vh = np.linalg.svd(base, full_matrices=False)
            axis = vh[0]
            if np.dot(centroids[-1] - centroids[0], axis) < 0:
                axis = -axis
        else:
            axis = span[-1] - span[0]
            axis = axis / (np.linalg.norm(axis) or 1.0)
        centred = span - span.mean(axis=0)
        along = centred @ axis
        # Rise and radius are the mean of the *local* windows, not of the span
        # against one global axis. The beam is 66 residues and visibly bent, so
        # a single axis puts its C-alphas 4.6 A off it and would report the
        # project's best-known helix as not helical. The axis stays global
        # because that is what the tilt is measured against.
        local_rise, local_radius, _, _ = helical_windows(span, window)
        out.append(HelicalSegment(
            chain=chain, start=int(block_seq[first]), end=int(block_seq[last]),
            axis=axis, center=span.mean(axis=0),
            rise=float(local_rise.mean()) if len(local_rise) else float("nan"),
            radius=float(local_radius.mean()) if len(local_radius) else float("nan"),
            length=float(along.max() - along.min())))
    return sorted(out, key=lambda s: s.start)


def _merge_runs(runs: list) -> list:
    """Coalesce overlapping or abutting index runs within the same block."""
    merged: list = []
    for block_seq, block_xyz, first, last in sorted(
            runs, key=lambda r: (id(r[0]), r[2])):
        if merged and merged[-1][0] is block_seq and first <= merged[-1][3] + 1:
            previous = merged[-1]
            merged[-1] = (previous[0], previous[1], previous[2],
                          max(previous[3], last))
        else:
            merged.append((block_seq, block_xyz, first, last))
    return merged


# --------------------------------------------------------------------------
# The cross-helices
# --------------------------------------------------------------------------

def _linker_ranges(reference: str, n_helices: int | None = None,
                   period: int | None = None) -> list[tuple[int, int]]:
    """Residue ranges between consecutive 4-TM units, in ``reference`` numbering.

    A linker runs from the end of a unit's fourth helix to the start of the
    next unit's first. Only the inter-unit linkers are returned — the paper is
    specific that the cross-helices sit "between 4-TM units", and the loops
    inside a unit are 5 and 11 residues long on average, too short to hold one.
    """
    import json

    from ..config import RESOURCE_DIR

    # Derived, not chosen: the period is the repeat `hydropathy` measures and
    # `domains.json` is built on, and the helix count is however many complete
    # units of it the reference has.
    period = 4 if period is None else int(period)
    all_helices = sorted(json.loads(
        (RESOURCE_DIR / f"uniprot_{reference}.json").read_text())["transmembrane"],
        key=lambda t: t["start"])
    if n_helices is None:
        n_helices = ((len(all_helices) - 2) // period) * period
    helices = all_helices[:n_helices]
    return [(helices[i]["end"] + 1, helices[i + 1]["start"] - 1)
            for i in range(period - 1, len(helices) - 1, period)]


def _named_ranges(reference: str, ids: tuple[str, ...]) -> list[tuple[int, int]]:
    """Residue ranges of named domains, in ``reference`` numbering."""
    import json

    from ..config import RESOURCE_DIR

    key = "mouse" if reference.startswith("mouse") else "human"
    domains = json.loads(
        (RESOURCE_DIR / "domains.json").read_text())["domains"]
    return [(d[key]["start"], d[key]["end"]) for d in domains if d["id"] in ids]


#: Elements that are helical, sit in an inter-unit linker, and are *not*
#: cross-helices — because Figure 7b colours them separately. The beam is red
#: there and yellow would be wrong; it also happens to be the longest helix in
#: the structure, so leaving it in would dominate every count and every mean.
EXCLUDED_FROM_CROSS = ("beam", "coiled_coil")


def cross_helices(structure, reference: str, axis: SymmetryAxis | None = None,
                  min_tilt: float | None = None, local_normal: bool = True,
                  exclude: tuple[str, ...] = EXCLUDED_FROM_CROSS
                  ) -> list[HelicalSegment]:
    """Helices in the inter-unit linkers that run across the membrane.

    ``local_normal`` measures the tilt against each protomer's own best-fit
    membrane plane rather than against the three-fold axis. That is the
    faithful reading of "perpendicular to the TM helices": the arms are
    themselves tilted about 30 degrees out of the pore's plane, so a helix
    genuinely perpendicular to its local membrane would come out at 60 degrees
    to the global axis and could be missed. Passing ``False`` uses the global
    axis and reports a different, also-defensible set; the two are compared in
    :func:`cross_helix_scan`.
    """
    from .geometry import tm_surface_by_chain

    if min_tilt is None:
        min_tilt = _P.value("cross_helix.min_tilt_deg")
    by_chain, resolved = tm_surface_by_chain(structure, reference)
    if axis is None:
        blocks = [np.array([pts[i] for i in sorted(resolved) if i in pts])
                  for pts in by_chain.values()]
        blocks = [b for b in blocks if len(b) >= 3]
        if len(blocks) < 3:
            raise ValueError("need three protomers to recover the axis")
        axis = detect_c3_axis(blocks)

    banned = _named_ranges(reference, exclude)

    def is_excluded(segment: HelicalSegment) -> bool:
        return any(segment.start <= hi and segment.end >= lo
                   for lo, hi in banned)

    found: list[HelicalSegment] = []
    for chain, helix_points in by_chain.items():
        points = np.array([helix_points[i] for i in sorted(resolved)
                           if i in helix_points])
        if local_normal and len(points) >= 3:
            normal = fit_plane(points).normal
        else:
            normal = axis.direction
        for lo, hi in _linker_ranges(reference):
            for segment in helical_segments(structure, chain, lo, hi):
                if segment.tilt_to(normal) >= min_tilt and not is_excluded(segment):
                    found.append(segment)
    return sorted(found, key=lambda s: (s.chain, s.start))


@dataclass
class CrossHelixScan:
    """Count against threshold, plus what the transmembrane helices score."""

    thresholds: tuple[float, ...]
    counts: tuple[int, ...]
    #: Tilt of every helical segment found in a linker, degrees.
    linker_tilts: tuple[float, ...]
    #: Tilt of the annotated transmembrane helices, for contrast.
    transmembrane_tilts: tuple[float, ...]
    default_threshold: float
    meta: dict = field(default_factory=dict)

    @property
    def separated(self) -> bool:
        """True when no transmembrane helix reaches the linker helices' median.

        The claim the threshold rests on. If a transmembrane helix tilts as far
        as the median cross-helix then "runs across the membrane" is not
        picking out a distinct population and the yellow helices in Figure 7b
        are a choice of cut rather than a feature.
        """
        if not self.linker_tilts or not self.transmembrane_tilts:
            return False
        return max(self.transmembrane_tilts) < float(np.median(self.linker_tilts))


def cross_helix_scan(structure, reference: str) -> CrossHelixScan:
    """How sharp the cross-helix threshold is on this structure."""
    from .geometry import tm_surface_by_chain
    import json

    from ..config import RESOURCE_DIR

    by_chain, resolved = tm_surface_by_chain(structure, reference)
    annotated = sorted(json.loads(
        (RESOURCE_DIR / f"uniprot_{reference}.json").read_text())["transmembrane"],
        key=lambda t: t["start"])

    linker_tilts: list[float] = []
    tm_tilts: list[float] = []
    for chain, helix_points in by_chain.items():
        points = np.array([helix_points[i] for i in sorted(resolved)
                           if i in helix_points])
        if len(points) < 3:
            continue
        normal = fit_plane(points).normal
        for lo, hi in _linker_ranges(reference):
            for segment in helical_segments(structure, chain, lo, hi):
                linker_tilts.append(segment.tilt_to(normal))
        for index in sorted(resolved):
            if index > len(annotated):
                continue
            helix = annotated[index - 1]
            for segment in helical_segments(structure, chain, helix["start"],
                                            helix["end"], min_length=5):
                tm_tilts.append(segment.tilt_to(normal))

    grid = tuple(float(t) for t in range(0, 91, 5))
    counts = tuple(int(sum(1 for t in linker_tilts if t >= threshold))
                   for threshold in grid)
    return CrossHelixScan(
        thresholds=grid, counts=counts,
        linker_tilts=tuple(linker_tilts), transmembrane_tilts=tuple(tm_tilts),
        default_threshold=float(_P.value("cross_helix.min_tilt_deg")),
        meta={"n_linkers": len(_linker_ranges(reference)),
              "reference": reference,
              "normal": "each protomer's own best-fit membrane plane"})
