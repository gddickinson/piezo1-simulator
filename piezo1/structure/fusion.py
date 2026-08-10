"""PIEZO1–HaloTag fusion geometry, modelled and labelled as such.

Imaging work on PIEZO1 fuses **HaloTag** to its cytosolic C-terminus — three
copies per trimer, one per protomer. To reason about what a tagged channel looks
like, or what a calcium dye attached to it would see, the tags have to be
somewhere in space.

**There is no structure of the fusion.** The tag's own fold is experimental
(PDB 6U32, 1.8 A, with a TMR-HaloTag ligand covalently bound) and PIEZO1's is
experimental, but their *relative* placement is not measured and the linker
between them is flexible. So this module does not produce a pose and call it a
structure. It produces an **accessible volume** — the region the tag centre can
occupy without clashing — in the manner used for FRET position modelling, and
treats any single placement inside it as one draw from an ensemble.

Measured inputs, from the two experimental structures:

* PIEZO1's C-terminus (human residue 2521) sits at **z = -60.4 A** in the
  canonical frame, 23.4 A off the conduction axis, with the pore exit at
  z = -62.7 A.
* HaloTag has a radius of gyration of **17.6 A**, its N-terminus — the end the
  fusion attaches — lies **19.9 A** from its centre, and the bound ligand sits
  18.6 A from that centre.

The one number that is *assumed* rather than measured is the linker length,
which no source for this construct states. It is registered as
``fusion.linker_residues`` with the ``unverified`` sentinel, because it sets how
far the tag can reach and therefore drives every distance reported here. Vary
it; do not trust it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from ..config import STRUCTURE_DIR
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .protomers import protomer_blocks, well_resolved_chains
from .superpose import SymmetryAxis, detect_c3_axis, rotation_matrix

__all__ = ["HaloTag", "AccessibleVolume", "FusionModel", "load_halotag",
           "cterm_anchors", "accessible_volume", "build_fusion",
           "HALOTAG_PDB", "SOLVENT_NAMES"]

#: HaloTag with a TMR-HaloTag ligand covalently bound, 1.8 A. Chosen over an
#: apo structure because the conjugate is what a labelling experiment produces.
HALOTAG_PDB = "6U32"

#: Het residues that are crystallisation solvent rather than the bound dye.
#: Shared with ``fusion_pose``, which draws the same atom set it measures.
SOLVENT_NAMES = ("HOH", "CL", "NA", "K", "MG", "CA", "ZN", "SO4")


@dataclass
class HaloTag:
    """The tag's own experimental structure, with its attachment point."""

    structure: Structure
    anchor: np.ndarray            # N-terminal CA — the end fused to PIEZO1
    centre: np.ndarray            # centroid of the protein atoms
    ligand: np.ndarray | None     # centroid of the bound dye, if resolved
    radius_of_gyration: float
    max_extent: float

    @property
    def anchor_to_centre(self) -> float:
        return float(np.linalg.norm(self.centre - self.anchor))

    @property
    def anchor_to_ligand(self) -> float:
        if self.ligand is None:
            return float("nan")
        return float(np.linalg.norm(self.ligand - self.anchor))

    def summary(self) -> str:
        return (f"HaloTag {HALOTAG_PDB}: Rg {self.radius_of_gyration:.1f} A, "
                f"N-terminus {self.anchor_to_centre:.1f} A from centre, "
                f"ligand {self.anchor_to_ligand:.1f} A from the N-terminus")


@dataclass
class AccessibleVolume:
    """Where the tag centre can sit without clashing, given a flexible linker."""

    points: np.ndarray            # (n, 3) grid points that survive
    anchor: np.ndarray
    spacing: float
    reach: float
    tag_radius: float
    n_before_clash: int = 0
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.points)

    @property
    def volume(self) -> float:
        """nm^3 — grid count times cell volume."""
        return float(len(self.points) * self.spacing ** 3 / 1000.0)

    @property
    def occluded_fraction(self) -> float:
        """How much of the tether's reach the channel body takes away."""
        if not self.n_before_clash:
            return float("nan")
        return 1.0 - len(self.points) / self.n_before_clash

    @property
    def centroid(self) -> np.ndarray:
        """Mean accessible position — the representative the model reports.

        An extremum of the envelope would be an equally valid single draw, but
        it would also be the most favourable one, and quoting it as *the* tag
        position would overstate how far the tag reaches.
        """
        return (self.points.mean(axis=0) if len(self.points)
                else self.anchor.copy())

    def distances_from(self, point: np.ndarray) -> np.ndarray:
        return np.linalg.norm(self.points - np.asarray(point, dtype=float),
                              axis=1)

    def summary(self) -> str:
        return (f"{len(self.points)} positions, {self.volume:.1f} nm^3 within "
                f"{self.reach:.0f} A of the anchor; "
                f"{self.occluded_fraction:.0%} of the reach blocked by the channel")


@dataclass
class FusionModel:
    """PIEZO1 with a HaloTag modelled at each protomer's C-terminus."""

    tag_centres: np.ndarray       # (n, 3)
    anchors: np.ndarray           # (n, 3)
    volume: AccessibleVolume | None = None
    axis: SymmetryAxis | None = None
    pore_exit: np.ndarray | None = None
    anchor_residues: tuple[int, ...] = ()
    meta: dict = field(default_factory=dict)

    @property
    def n_tags(self) -> int:
        return len(self.tag_centres)

    def c3_deviation(self) -> float:
        """Largest departure from perfect three-fold symmetry, in A.

        The tags are placed by rotating one solution about the recovered axis,
        so this should be numerically zero. It is measured rather than assumed
        because a placement that quietly broke the symmetry would look right in
        a picture and be wrong in every calculation that followed.
        """
        if self.axis is None or self.n_tags < 2:
            return float("nan")
        worst = 0.0
        for step in range(1, self.n_tags):
            matrix = rotation_matrix(self.axis.direction,
                                     2.0 * np.pi * step / self.n_tags)
            expected = ((self.tag_centres[0] - self.axis.point) @ matrix.T
                        + self.axis.point)
            worst = max(worst, float(np.linalg.norm(
                expected - self.tag_centres[step])))
        return worst

    def pore_exit_distances(self) -> np.ndarray:
        """Tag centre to pore exit, in nm — the number the calcium work needs."""
        if self.pore_exit is None:
            return np.full(self.n_tags, np.nan)
        return np.linalg.norm(self.tag_centres - self.pore_exit, axis=1) / 10.0

    def seams(self) -> np.ndarray:
        """Anchor-to-tag-centre segments, ``(n, 2, 3)``.

        The seam is the part with no experimental support, so it is returned
        explicitly for a renderer to draw differently from the rest.
        """
        return np.stack([self.anchors, self.tag_centres], axis=1)

    def summary(self) -> str:
        d = self.pore_exit_distances()
        return (f"{self.n_tags} HaloTags modelled; C3 deviation "
                f"{self.c3_deviation():.3f} A; tag centre "
                f"{np.nanmean(d):.1f} nm from the pore exit; "
                f"clearance {self.meta.get('min_clearance', float('nan')):.1f} A "
                f"— PLACEMENT IS A MODEL, not a measurement")


def load_halotag(pdb: str = HALOTAG_PDB) -> HaloTag:
    """Read the tag's experimental structure and locate its attachment point."""
    path = STRUCTURE_DIR / f"{pdb.upper()}.cif"
    if not path.exists():
        raise FileNotFoundError(
            f"{pdb} not downloaded — run python -m piezo1.io.fetch")
    structure = Structure.from_file(path)

    protein = structure.mask_protein() & ~structure.hetero
    coords = structure.xyz[protein].astype(np.float64)
    centre = coords.mean(axis=0)

    ca = structure.mask_ca() & ~structure.hetero
    numbers = structure.res_seq[ca]
    # The fusion joins PIEZO1's C-terminus to the tag's N-terminus, so it is
    # the lowest-numbered resolved residue that anchors it.
    anchor = structure.xyz[ca & (structure.res_seq == numbers.min())][0]

    ligand_mask = structure.hetero & ~np.isin(structure.res_name,
                                              list(SOLVENT_NAMES))
    ligand = (structure.xyz[ligand_mask].astype(np.float64).mean(axis=0)
              if ligand_mask.any() else None)

    return HaloTag(
        structure=structure, anchor=anchor.astype(np.float64), centre=centre,
        ligand=ligand,
        radius_of_gyration=float(np.sqrt(((coords - centre) ** 2).sum(1).mean())),
        max_extent=float(np.linalg.norm(coords - centre, axis=1).max()))


def cterm_anchors(structure: Structure) -> tuple[np.ndarray, tuple[int, ...]]:
    """The last resolved C-alpha of each protomer — where a fusion attaches.

    Returns the coordinates and the residue numbers they came from. The residue
    numbers are returned because they differ between entries: a construct
    truncated before 2521 anchors the tag somewhere else entirely, and a
    distance quoted without saying which residue it started from cannot be
    compared with anything.
    """
    coords, residues = [], []
    for chain in well_resolved_chains(structure)[:3]:
        mask = (structure.mask_ca() & ~structure.hetero
                & (structure.chain == chain))
        if not mask.any():
            continue
        seq, xyz = structure.res_seq[mask], structure.xyz[mask]
        last = int(np.argmax(seq))
        coords.append(xyz[last].astype(np.float64))
        residues.append(int(seq[last]))
    return np.array(coords), tuple(residues)


def accessible_volume(host: Structure, anchor: np.ndarray, tag: HaloTag,
                      linker_residues: int | None = None,
                      spacing: float | None = None) -> AccessibleVolume:
    """Positions the tag centre can occupy without clashing with the host.

    The linker is treated as a flexible tether: its maximum reach is the number
    of residues times a per-residue extension, plus the tag's own
    anchor-to-centre offset. Within that reach, a position is kept if a sphere
    of the tag's radius of gyration placed there does not overlap the channel.

    This is the accessible-volume construction used in FRET position modelling.
    It deliberately gives a **region rather than a pose**: the linker is
    flexible, the tag samples an ensemble, and drawing one orientation would
    imply a precision the model does not have.

    Two approximations worth stating. The tag is a sphere of its radius of
    gyration, not its full 30 A extent, so a real tag in an unfavourable
    orientation clashes where this says it does not. And the reach uses a fully
    extended linker, which no real flexible linker adopts. Both err towards a
    *larger* envelope, so a distance this model calls reachable is an upper
    bound.
    """
    linker_residues = int(_P.value("fusion.linker_residues")
                          if linker_residues is None else linker_residues)
    spacing = float(_P.value("fusion.grid_spacing")
                    if spacing is None else spacing)
    per_residue = _P.value("fusion.residue_extension")
    clearance = _P.value("fusion.clash_clearance")

    anchor = np.asarray(anchor, dtype=np.float64)
    reach = linker_residues * per_residue + tag.anchor_to_centre
    tag_radius = tag.radius_of_gyration

    steps = np.arange(-reach, reach + spacing, spacing)
    grid = np.stack(np.meshgrid(steps, steps, steps, indexing="ij"), axis=-1)
    points = grid.reshape(-1, 3) + anchor
    points = points[np.linalg.norm(points - anchor, axis=1) <= reach]
    n_before = len(points)

    protein = host.mask_protein() & ~host.hetero
    tree = cKDTree(host.xyz[protein].astype(np.float64))
    # A tag centre closer than (its own radius + clearance) to any channel atom
    # would have the tag body overlapping the channel.
    hits = tree.query_ball_point(points, tag_radius + clearance,
                                 return_length=True)
    points = points[hits == 0]

    return AccessibleVolume(
        points=points, anchor=anchor, spacing=spacing, reach=reach,
        tag_radius=tag_radius, n_before_clash=n_before,
        meta={"linker_residues": linker_residues,
              "per_residue_extension": per_residue,
              "anchor_to_centre": tag.anchor_to_centre,
              "clearance": clearance,
              "note": "accessible volume, not a pose; the linker is flexible "
                      "and the tag samples an ensemble"})


def _pore_exit(host: Structure, axis: SymmetryAxis,
               anchors: np.ndarray) -> np.ndarray:
    """Where the conduction pathway opens into the cytosol.

    Two traps here, both of which produce a confident wrong number.

    Taking the lowest protein atom *anywhere* is wrong: in a curved trimer that
    is a distal blade tip, far off-axis. Restricting to atoms within
    ``PORE_MOUTH_RADIUS`` of the axis finds the CTD bundle instead, which is
    what a permeating ion actually leaves through.

    And the sign of ``SymmetryAxis.direction`` is not fixed. It comes from the
    rotation operator relating two protomers, and re-detecting the axis on an
    already-canonically-framed structure returns −z as readily as +z — it does
    for 7WLT and 8YFG but not for 8YEZ. Trusting it put those two structures'
    pore exit at the *extracellular* end, and the tag then appeared to sit
    15-16 nm from it against 3.9 nm for the same construct on 8YEZ. So the
    cytosolic direction is taken from the C-terminal anchors, which are
    intracellular by topology, rather than from the axis.
    """
    protein = host.mask_protein() & ~host.hetero
    xyz = host.xyz[protein].astype(np.float64)

    direction = np.asarray(axis.direction, dtype=np.float64)
    direction = direction / np.linalg.norm(direction)
    if float(((np.atleast_2d(anchors) - axis.point) @ direction).mean()) < 0.0:
        direction = -direction

    near = axis.radial(xyz) <= _P.value("fusion.pore_mouth_radius")
    if not near.any():                       # nothing on-axis; fall back
        near = np.ones(len(xyz), dtype=bool)
    reach = float(((xyz[near] - axis.point) @ direction).max())
    return axis.point + reach * direction


def build_fusion(host: Structure, tag: HaloTag | None = None,
                 axis: SymmetryAxis | None = None,
                 linker_residues: int | None = None,
                 spacing: float | None = None) -> FusionModel:
    """Model one tag per protomer, placed C3-symmetrically.

    The accessible volume is computed once, for the first anchor, and its
    centroid is then **rotated** onto the other protomers rather than searched
    for independently. A homotrimer carries three identical copies of the
    construct, so an independent search would let grid noise break a symmetry
    the molecule actually has.
    """
    if tag is None:
        tag = load_halotag()
    anchors, residues = cterm_anchors(host)
    if len(anchors) < 1:
        raise ValueError("no C-terminal anchor found — is this a protomer?")

    blocks, _ = protomer_blocks(host)
    if axis is None:
        if len(blocks) < 3:
            raise ValueError("need a trimer to place tags symmetrically")
        axis = detect_c3_axis(blocks)

    volume = accessible_volume(host, anchors[0], tag,
                               linker_residues=linker_residues,
                               spacing=spacing)
    if len(volume) == 0:
        raise RuntimeError(
            "no clash-free position for the tag; the linker may be too short "
            "for it to clear the channel")

    centre = volume.centroid
    centres = [centre]
    for step in range(1, len(anchors)):
        matrix = rotation_matrix(axis.direction,
                                 2.0 * np.pi * step / len(anchors))
        centres.append((centre - axis.point) @ matrix.T + axis.point)
    centres = np.array(centres)

    protein = host.mask_protein() & ~host.hetero
    pore_exit = _pore_exit(host, axis, anchors)

    tree = cKDTree(host.xyz[protein].astype(np.float64))
    clearance = float(min(tree.query(c)[0] for c in centres))

    return FusionModel(
        tag_centres=centres, anchors=anchors, volume=volume, axis=axis,
        pore_exit=pore_exit, anchor_residues=residues,
        meta={"tag_pdb": HALOTAG_PDB,
              "accessible_volume_nm3": volume.volume,
              "accessible_positions": len(volume),
              "linker_residues": volume.meta["linker_residues"],
              "min_clearance": clearance,
              "tag_radius": tag.radius_of_gyration,
              "clashes": bool(clearance < tag.radius_of_gyration),
              "note": "PLACEMENT IS A MODEL — there is no structure of the "
                      "fusion, and the linker is flexible"})
