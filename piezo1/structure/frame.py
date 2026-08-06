"""Putting every structure into the same frame.

Deposited PIEZO entries arrive in whatever orientation the depositor's
refinement happened to end in. Nothing about a PDB frame is canonical: 8YEZ,
7WLT and 6B3R sit at unrelated angles in unrelated places. Loading one after
another therefore looks like the molecule jumping around, and overlaying two of
them shows nothing until they are superposed.

Two ways to fix that, offered separately because they answer different
questions:

* :func:`canonical_transform` — put *this* structure into a frame defined by its
  own symmetry: three-fold axis on **+z**, cytosolic side at **−z**, axis
  through the origin. It needs no other structure, works for PIEZO2 and for
  mouse entries, and is the honest choice when the two structures are not
  comparable residue-by-residue. Different conformational states still differ,
  because they genuinely do — this removes the framing, not the biology.

* :func:`reference_transform` — least-squares superpose onto a nominated
  reference over the C-alphas they share. This is what maximises overlap, and it
  is only meaningful when the residue numbering means the same thing in both,
  so it refuses a cross-species pair rather than silently aligning by number.

**The three-fold ambiguity is not a defect.** Fixing the roll about z can only
ever be defined up to 120°, because the molecule is C3-symmetric and a 120°
rotation maps it onto itself. Any of the three choices gives the same picture.

Convention throughout, matching :func:`~piezo1.structure.superpose.align_axis_to_z`::

    xyz_new = (xyz + translation) @ rotation.T
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from .protomers import protomer_blocks
from .superpose import (align_axis_to_z, detect_c3_axis, kabsch,
                        rmsd, rotation_matrix)

__all__ = ["Frame", "canonical_transform", "reference_transform",
           "apply_frame", "standardise", "ALIGNMENT_MODES"]

#: What the UI offers. ``deposited`` is the identity — the file as it came.
ALIGNMENT_MODES = ("deposited", "canonical", "reference")

#: Fraction of the highest-numbered resolved residues counted as "the cytosolic
#: end" when fixing the z sign. A fraction rather than a fixed count, so it
#: behaves the same on a 1,300-residue and a 2,500-residue protomer.
#:
#: **It has to be small.** The obvious choice of 10% is wrong: PIEZO's
#: extracellular cap runs to roughly residue 2457 of 2547, so the top tenth
#: straddles the cap *and* the intracellular CTD and its mean z depends on
#: which of the two is better resolved. On 7WLU and 11ZC that flips the sign and
#: the structure loads upside down while still reporting a perfect C3 fit.
#: Checked against the last 15 residues across all 20 downloaded entries, 0.02
#: and 0.05 agree everywhere, 0.10 fails on two and 0.20 fails on nineteen.
CTERM_FRACTION = 0.02


@dataclass
class Frame:
    """A rigid transform, plus what it was derived from and how well it fitted."""

    rotation: np.ndarray                  # (3, 3)
    translation: np.ndarray               # (3,) applied BEFORE the rotation
    mode: str = "canonical"
    rmsd: float | None = None             # reference mode only
    n_atoms_fitted: int = 0
    axis_rmsd: float | None = None        # canonical mode: C3 fit quality
    reordered: bool = False               # reference mode: protomers permuted
    note: str = ""
    meta: dict = field(default_factory=dict)

    @property
    def is_identity(self) -> bool:
        return (np.allclose(self.rotation, np.eye(3), atol=1e-9)
                and np.allclose(self.translation, 0.0, atol=1e-9))

    def apply(self, xyz: np.ndarray) -> np.ndarray:
        return (np.asarray(xyz, dtype=np.float64)
                + self.translation) @ self.rotation.T

    def summary(self) -> str:
        if self.mode == "deposited":
            return "deposited frame (unchanged)"
        if self.mode == "reference":
            extra = ", protomers reordered" if self.reordered else ""
            return (f"superposed on reference: {self.rmsd:.2f} Å RMSD over "
                    f"{self.n_atoms_fitted} C-alphas{extra}")
        return (f"canonical frame: C3 axis on +z "
                f"(axis fit {self.axis_rmsd:.2f} Å), cytosolic side at −z")


def _identity(note: str = "") -> Frame:
    return Frame(rotation=np.eye(3), translation=np.zeros(3), mode="deposited",
                 note=note)


def canonical_transform(structure: Structure,
                        reference: Structure | None = None) -> Frame:
    """Frame defined by the structure's own three-fold symmetry.

    Three-fold axis on +z through the origin, and the z sign chosen so the
    cytosolic end — the highest-numbered resolved residues — is at negative z.
    Without that sign rule the axis direction recovered from the symmetry
    operator is arbitrary, and half the structures would load upside down.

    That leaves the roll about z. Given a ``reference`` already in this frame
    and in the same numbering system, it is solved for in closed form by
    :func:`_choose_roll`; without one it is set by putting protomer 0's centroid
    on +x, which is reproducible but depends on which residues are resolved.

    Measured over the twenty downloaded entries, this takes structures that sit
    29–147 Å apart in their deposited frames to within 0.9–25 Å, and in every
    case to within about a angstrom of what an unconstrained least-squares
    superposition achieves — 19.77 Å against 19.70 Å for 7WLU on 7WLT, 0.91 Å
    against 0.91 Å for 8ZU8 on 8YEZ. What is left over is real conformational
    difference, not framing. The one large residual, 6KG7 at 57.8 Å, is PIEZO2
    rather than PIEZO1, and its free superposition is no better at 55.4 Å.
    """
    blocks, residues = protomer_blocks(structure)
    if len(blocks) < 3:
        return _identity("not a trimer — no three-fold axis to align to")

    axis = detect_c3_axis(blocks)

    # Sign the axis by topology, not by whatever the eigenvector solver
    # returned: PIEZO's C terminus is intracellular, so it must finish at
    # negative z. The test applies the candidate transform and measures where
    # the C-terminal residues actually land, rather than predicting it from the
    # axis projection. Predicting it read the C terminus off the protomer
    # blocks, whose per-chain `searchsorted` indexing assumes each chain's
    # residues are already in ascending order — on 7WLU they are not, and the
    # structure loaded upside down while reporting a perfect C3 fit.
    mask = structure.mask_ca() & ~structure.hetero
    seq, ca = structure.res_seq[mask], structure.xyz[mask].astype(np.float64)
    cterm = ca[seq >= np.quantile(seq, 1.0 - CTERM_FRACTION)]

    rot, trans = align_axis_to_z(axis, flip=False)
    flip = bool(len(cterm) and ((cterm + trans) @ rot.T)[:, 2].mean() > 0.0)
    if flip:
        rot, trans = align_axis_to_z(axis, flip=True)

    # Roll: put protomer 0's centroid on +x, giving a reproducible starting
    # point. Which of the three it is gets settled below.
    centroid = ((blocks[0].mean(axis=0) + trans) @ rot.T)
    angle = float(np.arctan2(centroid[1], centroid[0]))
    base = rotation_matrix(np.array([0.0, 0.0, 1.0]), -angle) @ rot

    extra, chosen_rmsd = 0.0, None
    if reference is not None:
        extra, chosen_rmsd = _choose_roll(structure, base, trans, reference)
    rotation = (rotation_matrix(np.array([0.0, 0.0, 1.0]), extra) @ base
                if extra else base)

    return Frame(rotation=rotation, translation=trans, mode="canonical",
                 axis_rmsd=float(axis.rmsd), n_atoms_fitted=len(blocks[0]),
                 rmsd=chosen_rmsd,
                 meta={"flipped": flip, "roll_deg": float(np.degrees(-angle)),
                       "roll_fitted_deg": float(np.degrees(extra)),
                       "roll_from_reference": reference is not None,
                       "c3_angle_deg": float(axis.angle_deg)})


def _choose_roll(structure: Structure, base: np.ndarray, trans: np.ndarray,
                 reference: Structure) -> tuple[float, float | None]:
    """Pick the roll about z that best matches a reference, in closed form.

    The starting rule — protomer 0's centroid on +x — is **coverage dependent**,
    and that turns out to matter. 7WLU resolves from residue 576 and 7WLT only
    from 784, so the extra blade density swings the centroid's azimuth by 47.7°
    between two entries of the same protein. That is not a multiple of 120°, so
    no choice among the three symmetry-equivalent rolls can undo it, and the two
    structures sit 71.6 Å apart while a free superposition reaches 19.7 Å.

    So when a reference is available the roll is solved for rather than picked
    from three. Rotation about a fixed axis is a one-parameter Procrustes
    problem with a closed-form optimum::

        θ* = atan2( Σ (xₘ y_t − yₘ x_t),  Σ (xₘ x_t + yₘ y_t) )

    This is a constrained fit, and worth being clear about: the axis and the
    z sign still come from the molecule's own symmetry, so the frame is not a
    least-squares superposition. Only the roll — the one degree of freedom the
    symmetry genuinely leaves free — is set by the reference.

    Returns the angle in radians and the RMSD it achieved, or ``(0.0, None)``
    when the two share too few residues for the comparison to mean anything.
    """
    best_angle, best_rmsd = 0.0, np.inf
    for perm in PERMUTATIONS:
        moving, target, n_common = _corresponding_ca(structure, reference, perm)
        if n_common < 3:
            return 0.0, None
        placed = (moving + trans) @ base.T
        x_m, y_m = placed[:, 0], placed[:, 1]
        x_t, y_t = target[:, 0], target[:, 1]
        angle = float(np.arctan2((x_m * y_t - y_m * x_t).sum(),
                                 (x_m * x_t + y_m * y_t).sum()))
        rolled = placed @ rotation_matrix(np.array([0.0, 0.0, 1.0]), angle).T
        value = float(np.sqrt(((rolled - target) ** 2).sum(1).mean()))
        if value < best_rmsd:
            best_angle, best_rmsd = angle, value
    return best_angle, (None if not np.isfinite(best_rmsd) else best_rmsd)


#: The two protomer correspondences worth testing. All six permutations of three
#: protomers fall into just two classes — the cyclic ones and the reversed ones
#: — and within a class every member scores identically, because a cyclic
#: relabelling is the same thing as rotating the trimer by 120°. Measured on
#: 8YFG against 8YEZ: cyclic 72.88 Å in all three, reversed 12.45 Å in all
#: three. So one representative of each class is a complete search.
#:
#: The reversed class is not hypothetical. 8YFG and 8ZU3 both present chains
#: A, B, D, but the depositors numbered them round the ring in opposite senses,
#: and taking the labels at face value costs 60 Å of apparent RMSD.
PERMUTATIONS = ((0, 1, 2), (0, 2, 1))


def _corresponding_ca(mobile: Structure, reference: Structure,
                      perm: tuple[int, int, int] = (0, 1, 2)):
    """Matched C-alpha arrays for the two structures, protomer by protomer.

    Keying on residue number alone is not enough, and getting that wrong is
    quiet rather than loud: a residue number occurs once **per protomer**, so a
    ``dict`` keyed on it keeps whichever chain happened to be read last and
    throws away two thirds of the atoms. Worse, if the two structures' chains
    are ordered differently, the survivors are not equivalent to each other and
    the fit is between the wrong atoms.

    So the correspondence is built on the protomer blocks, restricted to the
    residues both structures resolve, and concatenated in block order. Returns
    ``(moving, target)``, both ``(3 * n_common, 3)``, or empty arrays when there
    is no usable overlap.
    """
    empty = (np.zeros((0, 3)), np.zeros((0, 3)), 0)
    mob_blocks, mob_res = protomer_blocks(mobile)
    ref_blocks, ref_res = protomer_blocks(reference)
    if len(mob_blocks) < 3 or len(ref_blocks) < 3:
        return empty

    common = np.array(sorted(set(mob_res.tolist()) & set(ref_res.tolist())))
    if len(common) < 3:
        return empty

    mob_keep = np.isin(mob_res, common)
    ref_keep = np.isin(ref_res, common)
    moving = np.concatenate([mob_blocks[i][mob_keep] for i in perm])
    target = np.concatenate([b[ref_keep] for b in ref_blocks])
    return moving, target, len(common)


def reference_transform(mobile: Structure, reference: Structure,
                        same_numbering: bool = True) -> Frame:
    """Least-squares superpose ``mobile`` onto ``reference`` over shared C-alphas.

    ``same_numbering`` must be true for the result to mean anything: residue 1000
    is a different residue in human PIEZO1 and mouse Piezo1, and matching by
    number across species would superpose non-equivalent atoms and report a
    confident RMSD for it.

    Deposited chain labels are not a reliable guide to rotational order, so the
    protomer correspondence is recovered rather than assumed — a trimer fitted
    with its subunits one third out of register gives a large RMSD that looks
    like a conformational difference.
    """
    if not same_numbering:
        return _identity("reference alignment needs a shared numbering system")

    # Chain labels are not a reliable guide to rotational order, so both
    # correspondence classes are fitted and the better kept. A trimer fitted
    # with its subunits out of register gives a large RMSD that reads as a
    # conformational difference rather than a bookkeeping error.
    best = None
    for perm in PERMUTATIONS:
        moving, fixed, n_common = _corresponding_ca(mobile, reference, perm)
        if n_common < 3:
            return _identity(f"only {n_common} shared residues — too few to fit")
        rotation, translation, centroid = kabsch(moving, fixed)
        fitted = (moving - centroid) @ rotation.T + translation
        error = float(rmsd(fitted, fixed))
        if best is None or error < best[0]:
            best = (error, perm, rotation, translation, centroid, len(moving),
                    n_common)

    error, perm, rotation, translation, centroid, n_fitted, n_common = best
    # kabsch centres the mobile set, rotates, then adds the target centroid:
    #   fitted = (x - centroid) @ R.T + translation
    # This module's convention is (x + t) @ R.T, so fold the post-rotation
    # shift back through the rotation:  t = -centroid + translation @ R.
    offset = -centroid + translation @ rotation
    return Frame(rotation=rotation, translation=offset, mode="reference",
                 rmsd=error, n_atoms_fitted=n_fitted,
                 reordered=(perm != (0, 1, 2)),
                 meta={"shared_residues": n_common, "protomer_order": perm})


def apply_frame(structure: Structure, frame: Frame,
                name: str | None = None) -> Structure:
    """Return a copy of ``structure`` with the transform applied to every atom.

    Every atom, not just the fitted C-alphas: ligands, lipids and ions have to
    travel with the protein they are bound to.
    """
    xyz = frame.apply(structure.xyz)
    return structure.copy_with_coords(xyz.astype(structure.xyz.dtype),
                                      name=name or structure.name)


def standardise(structure: Structure, mode: str = "canonical",
                reference: Structure | None = None,
                same_numbering: bool = True) -> tuple[Structure, Frame]:
    """Put ``structure`` into the requested frame, returning it and the transform.

    Falls back to the canonical frame — with the reason recorded on the returned
    :class:`Frame` — when a reference alignment was asked for but cannot be
    computed. Silently returning the deposited frame would look identical to
    success in the viewport.
    """
    if mode not in ALIGNMENT_MODES:
        raise ValueError(f"unknown alignment mode {mode!r}; "
                         f"expected one of {ALIGNMENT_MODES}")
    if mode == "deposited":
        return structure, _identity()

    # A reference only disambiguates the roll if it is in the same frame and the
    # same numbering system; otherwise the comparison is between residues that
    # are not equivalent, and it would pick a roll for a reason that is not real.
    roll_reference = reference if (reference is not None and same_numbering) else None

    if mode == "reference":
        if reference is None:
            frame = canonical_transform(structure)
            frame.note = "no reference loaded — fell back to the canonical frame"
            return apply_frame(structure, frame), frame
        frame = reference_transform(structure, reference,
                                    same_numbering=same_numbering)
        if frame.mode == "deposited":       # refused; say so and still align
            reason = frame.note
            frame = canonical_transform(structure, reference=roll_reference)
            frame.note = f"{reason} — fell back to the canonical frame"
        return apply_frame(structure, frame), frame

    frame = canonical_transform(structure, reference=roll_reference)
    return apply_frame(structure, frame), frame
