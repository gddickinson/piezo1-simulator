"""The real HaloTag fold placed on the channel, as the arbitrary draw it is.

`fusion.py` deliberately refuses to produce a pose. It gives an **accessible
volume** — the region the tag centre can occupy — because the linker is flexible
and no experiment has determined where the tag sits. Drawing a single sphere of
the tag's radius of gyration is the visual form of that refusal.

But a sphere is also unreadable: it says nothing about the fold, the dye, or how
much room the tag really needs. So this module places the *experimental* HaloTag
structure (6U32) at the modelled centre and states exactly which degrees of
freedom that placement does and does not determine:

* **Position — from the model.** The tag centre goes to the centroid of the
  accessible volume, the same point the sphere occupies, so every distance the
  fusion model reports still describes the thing on screen.
* **Two rotations — from the model.** The tag is turned so its own N-terminus,
  the residue the linker joins, points back at PIEZO1's C-terminus. That is
  forced by the geometry: any other orientation stretches the linker further
  than it needs to go.
* **The third rotation — undetermined.** Spinning the tag about the seam moves
  no endpoint and costs the linker nothing, so nothing in the model picks it.
  It is a free angle, and the interface lets a user turn it precisely so that
  its arbitrariness is visible rather than asserted in a caption.

**The measurement this makes possible.** `fusion.accessible_volume` treats the
tag as a sphere of its radius of gyration (17.6 A) and says so, noting that the
real fold reaches 30.0 A and so "clashes where this says it does not". That was
an acknowledged approximation; with the real coordinates in place it becomes a
number. Over 36 spins at the modelled centre:

=========  ===================  ========================  ==================
Structure  sphere model         orientations clearing     worst body contact
                                the channel
=========  ===================  ========================  ==================
7WLT       fits (24.0 A)        27 / 36                   9 atoms
8YFG       fits (21.7 A)        7 / 36                    12 atoms
8YEZ       fits (21.5 A)        1 / 36                    23 atoms
11ZC       clashes (15.7 A)     0 / 36                    37 atoms
=========  ===================  ========================  ==================

The two models agree on the question that matters — the one structure whose
sphere clearance falls below the radius of gyration is the one where no
orientation of the real fold clears — while the sphere is generous about *how
much* room there is, by the margin its own docstring predicted. Neither is
quietly corrected: shrinking the envelope until the fold fits everywhere would
move the tag position every other number in the fusion work is computed from.

**The contact that had to be excluded first.** Counting every touching atom
said 0/36 orientations clear on all four structures, which looked like the fold
contradicting the sphere. It did not. The persistent contact was the tag's own
N-terminal residue against PIEZO1's C-terminus — the two ends of the linker,
which this module's placement rule deliberately points at each other. A rule
cannot then report its own construction as a finding, so the anchor residue is
excluded and only the fold body is counted. That single correction is what
turned a manufactured disagreement into the agreement above.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .fusion import SOLVENT_NAMES, FusionModel, HaloTag, load_halotag
from .superpose import rotation_matrix

__all__ = ["TagPose", "align_rotation", "drawable_mask", "place_tag",
           "spin_scan", "best_spin", "pose_for_display", "SPIN_SAMPLES"]

#: How finely the undetermined spin is sampled when reporting what fraction of
#: orientations touch the channel, and when choosing the one to draw. A
#: reporting resolution: finer sampling sharpens the fraction and the chosen
#: angle, it cannot change which orientations are admissible.
SPIN_SAMPLES = 36


@dataclass
class TagPose:
    """One drawn orientation of the tag, with what it cost in contacts."""

    coords: np.ndarray            # (n_tags, n_atoms, 3)
    radii: np.ndarray             # (n_atoms,) van der Waals
    ligand: np.ndarray            # (n_atoms,) bool — the covalently bound dye
    body: np.ndarray              # (n_atoms,) bool — not the anchor residue
    anchors: np.ndarray           # (n_tags, 3) placed tag N-termini
    seams: np.ndarray             # (n_tags, 2, 3) channel C-term -> tag N-term
    spin: float                   # radians about the seam; undetermined
    linker_gap: float             # A, channel C-terminus to tag N-terminus
    touching: np.ndarray          # (n_atoms,) bool — which atoms touch
    contact_distance: float
    meta: dict = field(default_factory=dict)

    @property
    def contacts(self) -> int:
        """Every tag atom within the contact distance of the channel."""
        return int(self.touching.sum())

    @property
    def body_contacts(self) -> int:
        """Contacts excluding the anchor residue — the ones that mean anything.

        The placement rule points the tag's N-terminal residue at PIEZO1's
        C-terminus, so a contact there is the rule's own construction rather
        than a property of the fold. Counting it made all four structures look
        as though no orientation clears; see the module docstring.
        """
        return int((self.touching & self.body).sum())

    @property
    def attachment_contacts(self) -> int:
        return self.contacts - self.body_contacts

    @property
    def n_tags(self) -> int:
        return len(self.coords)

    @property
    def n_atoms(self) -> int:
        return self.coords.shape[1]

    @property
    def clears(self) -> bool:
        """Whether the fold body avoids the channel in this orientation."""
        return self.body_contacts == 0

    def summary(self) -> str:
        return (f"HaloTag fold ({self.meta.get('tag_pdb', '?')}) at spin "
                f"{np.degrees(self.spin):.0f} deg: {self.body_contacts} of "
                f"{self.n_atoms} atoms within {self.contact_distance:.1f} A of "
                f"the channel, excluding the attachment residue; linker spans "
                f"{self.linker_gap:.1f} A. "
                f"ORIENTATION ABOUT THE SEAM IS UNDETERMINED")


def align_rotation(source: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Rotation taking ``source`` onto ``target`` about their common normal.

    Rodrigues' formula on the shortest arc. The antiparallel case has no
    shortest arc — every perpendicular axis works — so it is handled explicitly
    rather than left to divide by a vanishing sine.
    """
    a = np.asarray(source, dtype=np.float64)
    b = np.asarray(target, dtype=np.float64)
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)

    cross = np.cross(a, b)
    sine = float(np.linalg.norm(cross))
    cosine = float(a @ b)
    if sine < 1e-12:
        if cosine > 0.0:
            return np.eye(3)
        perpendicular = np.array([1.0, 0.0, 0.0])
        if abs(a @ perpendicular) > 0.9:
            perpendicular = np.array([0.0, 1.0, 0.0])
        axis = np.cross(a, perpendicular)
        return rotation_matrix(axis / np.linalg.norm(axis), np.pi)

    skew = np.array([[0.0, -cross[2], cross[1]],
                     [cross[2], 0.0, -cross[0]],
                     [-cross[1], cross[0], 0.0]])
    return np.eye(3) + skew + skew @ skew * ((1.0 - cosine) / sine ** 2)


def drawable_mask(structure: Structure) -> np.ndarray:
    """The tag atoms worth drawing: the protein and its bound dye, not solvent.

    Crystallisation additives are dropped rather than drawn, because a chloride
    from 6U32's buffer sitting beside PIEZO1 would read as a permeating ion.
    """
    return ~np.isin(structure.res_name, list(SOLVENT_NAMES))


def _channel_tree(host: Structure) -> cKDTree:
    """Channel heavy atoms, cached for the spin scan.

    The scan places the tag 36 times against the same channel, and rebuilding a
    30,000-point tree each time costs more than the queries do. One entry only,
    holding a reference to the host so its identity cannot be recycled, and
    keyed on the coordinate array as well so a reframed structure misses.
    """
    key = (id(host), id(host.xyz))
    if _TREE_CACHE.get("key") != key:
        protein = host.mask_protein() & ~host.hetero
        _TREE_CACHE.clear()
        _TREE_CACHE.update(
            key=key, host=host,
            tree=cKDTree(host.xyz[protein].astype(np.float64)))
    return _TREE_CACHE["tree"]


_TREE_CACHE: dict = {}


def _tag_atoms(tag: HaloTag) -> tuple[np.ndarray, ...]:
    """Coordinates, radii, and the two masks the pose needs: dye and body.

    The body mask drops the anchor residue — the lowest-numbered resolved one,
    the same residue `load_halotag` takes its anchor from — because the
    placement aims it at the channel and so cannot report touching it.
    """
    keep = drawable_mask(tag.structure)
    residues = tag.structure.res_seq[keep]
    protein = ~(tag.structure.hetero[keep])
    anchor_residue = residues[protein].min()
    return (tag.structure.xyz[keep].astype(np.float64),
            tag.structure.vdw_radii()[keep].astype(np.float64),
            tag.structure.hetero[keep],
            residues != anchor_residue)


def place_tag(host: Structure, model: FusionModel, tag: HaloTag | None = None,
              spin: float = 0.0, contact: float | None = None) -> TagPose:
    """Put the tag's real coordinates at each modelled centre.

    The tag is moved rigidly — rotated about its own centre, then translated —
    so its internal geometry is untouched and it remains the deposited
    structure rather than a model of one. The copies on the other two protomers
    are C3 images of the first, matching how `build_fusion` placed the centres.
    """
    tag = tag or load_halotag()
    contact = float(_P.value("fusion.pose_contact_distance")
                    if contact is None else contact)

    coords, radii, ligand, body = _tag_atoms(tag)
    channel_anchor = np.asarray(model.anchors[0], dtype=np.float64)
    centre = np.asarray(model.tag_centres[0], dtype=np.float64)

    # The seam runs from the channel's C-terminus to the tag's N-terminus, so
    # the tag's own centre->N-terminus vector must point back down it.
    seam = channel_anchor - centre
    rotation = align_rotation(tag.anchor - tag.centre, seam)
    axis = seam / np.linalg.norm(seam)
    rotation = rotation_matrix(axis, float(spin)) @ rotation

    placed = (coords - tag.centre) @ rotation.T + centre
    tag_anchor = (tag.anchor - tag.centre) @ rotation.T + centre

    frames, anchors = [placed], [tag_anchor]
    if model.axis is not None:
        for step in range(1, len(model.tag_centres)):
            turn = rotation_matrix(model.axis.direction,
                                   2.0 * np.pi * step / len(model.tag_centres))
            about = model.axis.point
            frames.append((placed - about) @ turn.T + about)
            anchors.append((tag_anchor - about) @ turn.T + about)

    frames = np.asarray(frames)
    anchors = np.asarray(anchors)
    tree = _channel_tree(host)
    touching = tree.query_ball_point(frames[0], contact, return_length=True) > 0

    return TagPose(
        coords=frames, radii=radii, ligand=ligand, body=body, anchors=anchors,
        seams=np.stack([np.asarray(model.anchors)[:len(frames)], anchors],
                       axis=1),
        spin=float(spin),
        linker_gap=float(np.linalg.norm(seam) - tag.anchor_to_centre),
        touching=touching, contact_distance=contact,
        meta={"tag_pdb": model.meta.get("tag_pdb"),
              "max_extent": tag.max_extent,
              "radius_of_gyration": tag.radius_of_gyration,
              "note": "the deposited tag structure, rigidly placed; the spin "
                      "about the seam is undetermined and this is one draw"})


def spin_scan(host: Structure, model: FusionModel, tag: HaloTag | None = None,
              n_angles: int | None = None,
              contact: float | None = None) -> np.ndarray:
    """Body contacts at each sampled spin — the free angle, measured.

    A count per angle rather than a verdict, because the interesting statement
    is the shape of the whole curve: how many orientations clear the channel
    says how much room the tag really has, where a single yes/no would not.
    Counts exclude the anchor residue, for the reason `body_contacts` gives.
    """
    n_angles = int(SPIN_SAMPLES if n_angles is None else n_angles)
    angles = np.arange(n_angles) * 2.0 * np.pi / n_angles
    tag = tag or load_halotag()
    return np.array([
        place_tag(host, model, tag, spin=a, contact=contact).body_contacts
        for a in angles])


def best_spin(host: Structure, model: FusionModel, tag: HaloTag | None = None,
              n_angles: int | None = None,
              contact: float | None = None) -> tuple[float, np.ndarray]:
    """The least-contacting spin, and the whole scan it came from.

    Minimising contacts is a **drawing convention**, not a determination: the
    model does not prefer this angle, it is simply the draw that puts least of
    the tag inside the channel. Where several orientations clear it, the first
    is taken and there is nothing to choose between them. The scan is returned
    alongside so a caller can report how much choice there was.
    """
    counts = spin_scan(host, model, tag, n_angles=n_angles, contact=contact)
    n_angles = len(counts)
    return float(np.argmin(counts) * 2.0 * np.pi / n_angles), counts


def pose_for_display(host: Structure, model: FusionModel,
                     tag: HaloTag | None = None, spin: float | None = None,
                     contact: float | None = None) -> TagPose:
    """The pose an interface should draw, with the scan recorded on it.

    Pass ``spin`` to see a particular orientation; leave it out and the
    least-contacting one is chosen. Either way ``meta['clear_spins']`` says how
    many of the sampled orientations let the fold body avoid the channel
    entirely, so a caption can state whether the tag fits rather than implying
    it by drawing one that does.
    """
    tag = tag or load_halotag()
    if spin is None:
        spin, counts = best_spin(host, model, tag, contact=contact)
    else:
        counts = spin_scan(host, model, tag, contact=contact)
    pose = place_tag(host, model, tag, spin=spin, contact=contact)
    pose.meta.update(clear_spins=int((counts == 0).sum()),
                     spins_sampled=int(len(counts)),
                     fewest_contacts=int(counts.min()),
                     most_contacts=int(counts.max()))
    return pose
