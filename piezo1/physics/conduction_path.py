"""Which part of the axis the ions are allowed to use.

Every conduction number this project has produced treats the pore as a
one-dimensional channel running from bulk solvent, down the three-fold axis,
to bulk solvent on the other side. That is what HOLE-style profiling measures
and what :func:`piezo1.physics.permeation.solve_pnp` integrates over, and for
most channels it is right.

**For PIEZO1 it is wrong at both ends, and the literature says so.** Liu et al.
2025 (Neuron 113:590-604, PMID 39719701) report that Na+ reaches the cap
vestibule through *three lateral cap gates* rather than through the top of the
cap, which "remains closed above the residue R2295 position among all the
structures"; and that having crossed the transmembrane gate into the inner
vestibule, the ions leave through *intracellular lateral portals* rather than
through the vertical constriction neck, which "remains closed" even in their
intermediate-open structure. Their own 10 us simulations put 37 Na+ through a
single lateral portal and none through the neck.

Measured here, that is not a subtlety. On every deposited PIEZO1 entry the
axial profile is pinched below the water radius at **R2295 and its immediate
neighbours** at the top, and at the curated CTD constrictions at the bottom -
so an axial model refuses every structure, including the intermediate-open
8IXO whose transmembrane gate has demonstrably dilated (V2476 side-chain
diagonal 14.2 A against 7.7 A on the curved 7WLT, reproducing the paper's
7 -> 14 A).

So the pathway is a **choice**, and this module makes it one instead of an
assumption. ``axial`` is the default and returns the profile object unchanged,
so every recorded number is reproduced bit for bit and no claim moves. The
lateral options truncate the path at the curated endpoints, which is the
smallest change that lets the model represent the route the paper describes.

**What a lateral option does not do.** It does not model the portal. The
truncated end slice becomes the mouth, and the Hall access resistance is then
computed from a radius that is the *pore's* at that height, not the portal's -
which this project does not measure. A lateral result is therefore an upper
bound on how much the portal restricts the current, and
:meth:`ConductionPath.caveat` says so wherever it is shown.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass

import numpy as np

__all__ = ["PATHWAYS", "ConductionPath", "conduction_path", "PATHWAY_LABELS"]

#: The four routes. ``axial`` is today's model and the default everywhere.
PATHWAYS = ("axial", "lateral_entry", "lateral_exit", "lateral")

PATHWAY_LABELS = {
    "axial": "Axial (bulk to bulk down the three-fold axis)",
    "lateral_entry": "Lateral entry (in through the cap gates)",
    "lateral_exit": "Lateral exit (out through the intracellular portals)",
    "lateral": "Lateral entry and exit (Liu et al. 2025)",
}

#: Curated groups bounding the axial segment, both added in Round 84d with the
#: paper as their source. Named rather than numbered here so the residue has
#: one definition in the project and it carries a citation.
ENTRY_GROUP = "cap_constriction"        # R2295 (mouse) / R2279 (human)
EXIT_GROUP = "ctd_constriction"         # M2493/F2494, P2536/E2537 (mouse)


@dataclass(frozen=True)
class ConductionPath:
    """The stretch of the profile the conduction model may use."""

    profile: object                     # a PoreProfile, possibly truncated
    pathway: str = "axial"
    entry: str = "axial"                # how ions get in
    exit: str = "axial"                 # how they get out
    dropped_entry: int = 0              # slices removed at the extracellular end
    dropped_exit: int = 0               # slices removed at the cytosolic end
    refused: str = ""                   # why the truncation could not be made

    @property
    def is_axial(self) -> bool:
        return self.dropped_entry == 0 and self.dropped_exit == 0

    def caveat(self) -> str:
        """What the chosen pathway does and does not claim. Never omitted."""
        if self.refused:
            return (f"pathway {self.pathway!r} refused ({self.refused}) — "
                    f"the whole axis is being used, as in the axial model")
        if self.is_axial:
            return ("axial pathway: ions must traverse the whole axis, "
                    "including the closed cap top and the cytoplasmic "
                    "constriction neck that Liu et al. 2025 report are "
                    "bypassed laterally")
        parts = []
        if self.dropped_entry:
            parts.append(f"{self.dropped_entry} slices above the cap "
                         f"constriction")
        if self.dropped_exit:
            parts.append(f"{self.dropped_exit} slices below the cytoplasmic "
                         f"constriction")
        return (f"lateral pathway: {' and '.join(parts)} excluded — the "
                f"portal itself is NOT modelled, so the mouth radius is the "
                f"pore's at that height and the current is an upper bound")

    def summary(self) -> str:
        r = np.asarray(self.profile.radius, dtype=float)
        return (f"{PATHWAY_LABELS.get(self.pathway, self.pathway)}: "
                f"{len(r)} slices, narrowest {r.min() / 10:.3f} nm")


def conduction_path(structure, profile, pathway: str = "axial"
                    ) -> ConductionPath:
    """Restrict ``profile`` to the stretch ``pathway`` says ions use.

    ``axial`` returns the *same profile object*, so a caller that has not asked
    for anything else cannot get a different number than it did before this
    module existed. A test asserts identity, not equality.

    A truncation that cannot be made — numbering unreadable, endpoint residue
    unresolved — is **refused with the reason recorded** and the full axis
    returned, rather than silently truncating at whatever slice happened to be
    nearest. A confident lateral number on an entry whose cap is not resolved
    would be exactly the kind of quiet wrong answer this project audits for.
    """
    if pathway not in PATHWAYS:
        raise ValueError(f"unknown pathway {pathway!r}; expected one of "
                         f"{', '.join(PATHWAYS)}")
    if pathway == "axial":
        return ConductionPath(profile=profile, pathway=pathway)

    from ..core.annotations import load_annotations
    from ..core.numbering_check import piezo1_numbering
    from .pore_charge import cytosolic_end

    numbering = piezo1_numbering(structure)
    if numbering is None:
        return ConductionPath(profile=profile, pathway=pathway,
                              refused="numbering not readable from the "
                                      "coordinates; not a PIEZO1 entry")
    if profile.axis is None:
        return ConductionPath(profile=profile, pathway=pathway,
                              refused="profile carries no axis, so the two "
                                      "ends cannot be told apart")

    annotations = load_annotations(numbering)
    entry_res = _group_residues(annotations, ENTRY_GROUP)
    exit_res = _group_residues(annotations, EXIT_GROUP)

    # Extracellular first, so "before" and "after" mean what they say.
    order = (np.argsort(-np.asarray(profile.z))
             if cytosolic_end(structure, profile.axis) == 0
             else np.argsort(np.asarray(profile.z)))

    first, last = 0, len(order)
    entry_kind, exit_kind = "axial", "axial"

    if pathway in ("lateral_entry", "lateral"):
        hits = [k for k, i in enumerate(order)
                if entry_res.intersection(profile.slices[i].lining)]
        if not hits:
            return ConductionPath(profile=profile, pathway=pathway,
                                  refused="the cap constriction residue is not "
                                          "resolved, so where the cap gates "
                                          "open cannot be located")
        first = hits[-1] + 1
        entry_kind = "lateral cap gates"

    if pathway in ("lateral_exit", "lateral"):
        hits = [k for k, i in enumerate(order)
                if k >= first and exit_res.intersection(profile.slices[i].lining)]
        if not hits:
            return ConductionPath(profile=profile, pathway=pathway,
                                  refused="no cytoplasmic constriction residue "
                                          "lines the pore below the entry, so "
                                          "where the portals open cannot be "
                                          "located")
        last = hits[0]
        exit_kind = "intracellular lateral portals"

    if last - first < _MIN_SLICES:
        return ConductionPath(profile=profile, pathway=pathway,
                              refused=f"only {max(last - first, 0)} slices "
                                      f"survive the truncation, which is not a "
                                      f"pore")

    keep = np.sort(order[first:last])
    return ConductionPath(
        profile=_subset(profile, keep), pathway=pathway,
        entry=entry_kind, exit=exit_kind,
        dropped_entry=first, dropped_exit=len(order) - last)


#: A truncated path shorter than this is not a channel, it is a gap between two
#: annotations. Not a physical quantity — a sanity floor on the arithmetic.
_MIN_SLICES = 5


def _group_residues(annotations, group_id: str) -> set:
    group = annotations.group(group_id)
    return set(group.residues) if group is not None else set()


def _subset(profile, keep: np.ndarray):
    """A PoreProfile over ``keep`` only, with its provenance carried."""
    meta = dict(getattr(profile, "meta", {}) or {})
    meta["truncated_to"] = int(len(keep))
    return dataclasses.replace(
        profile,
        z=np.asarray(profile.z)[keep],
        radius=np.asarray(profile.radius)[keep],
        centers=np.asarray(profile.centers)[keep],
        slices=[profile.slices[i] for i in keep],
        meta=meta)
