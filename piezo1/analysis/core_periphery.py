"""The core agrees in space; the periphery does not. Measured, on any pair.

The census's last structural result is that superposing a predicted zebrafish
piezo3 model on the cryo-EM mouse Piezo1 structure **by the pore module alone**
puts the cores at 3.86 A over 448 C-alpha while the blades splay visibly apart —
turning a per-residue conservation profile into a shape.

This module generalises that from one pair to any pair this project holds
coordinates for, including across paralogues, and reports the thing the census
could only show as a picture: the ratio of blade RMSD to core RMSD after a
core-only fit.

**Why the fit is core-only and why that is not cheating.** Fitting on everything
shared spreads the error: a pair with a rigid common core and mobile blades
comes out looking uniformly mediocre, and a pair with no common core at all
comes out looking the same. Fitting on the core and *then measuring the
periphery* asks a directional question — given that these two channels' pores
are superposed, where are their blades? — and it can fail. If the core does not
superpose, :attr:`Comparison.core_converged` is False and the splay ratio is not
reported, because a ratio taken against a core that did not fit is arithmetic
rather than a measurement.

**Correspondence is never assumed.** Within one protein, residue numbers are
used after :func:`piezo1.core.numbering_check.identify_numbering` says which
numbering the entry is in. Across PIEZO1 and PIEZO2 it goes through
:func:`piezo1.analysis.paralogue.paralogue_map`, which is a real global
alignment — this project has a standing rule that no constant offset relates two
PIEZO sequences, and there is no exception here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.annotations import is_annotated, load_annotations
from ..core.numbering_check import identify_numbering
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.protomers import well_resolved_chains
from ..structure.superpose import kabsch, rmsd as _rmsd, superpose

__all__ = ["Comparison", "CoreFit", "Refusal", "compare", "core_fit",
           "core_residues", "periphery_residues", "correspondence",
           "CORE_DOMAINS", "PERIPHERY_DOMAINS"]

#: The core, in this project's domain vocabulary — the pore module.
CORE_DOMAINS = ("outer_helix", "cap", "inner_helix", "ctd")

#: The periphery. THU1–7 are the distal and proximal blade units; the anchor and
#: beam are excluded because they are the coupling between the two and belong to
#: neither side of a core-versus-periphery statement.
PERIPHERY_DOMAINS = ("thu1", "thu2", "thu3", "thu4", "thu5", "thu6", "thu7")


@dataclass(frozen=True)
class Refusal:
    """Why a pair could not be compared. Never a number."""

    reason: str

    def __bool__(self) -> bool:
        return False


@dataclass(frozen=True)
class Comparison:
    """One pair, fitted on the core and measured everywhere."""

    mobile_id: str
    target_id: str
    mobile_numbering: str
    target_numbering: str
    n_core: int
    core_rmsd: float
    n_periphery: int
    periphery_rmsd: float | None
    n_shared: int
    whole_rmsd: float
    cross_paralogue: bool
    note: str = ""

    def __bool__(self) -> bool:
        return True

    @property
    def core_converged(self) -> bool:
        return self.core_rmsd <= _P.value("family.core_rmsd_ceiling")

    @property
    def splay_ratio(self) -> float | None:
        """Blade RMSD over core RMSD, or None if the core did not fit."""
        if not self.core_converged or not self.periphery_rmsd or not self.core_rmsd:
            return None
        return self.periphery_rmsd / self.core_rmsd

    @property
    def core_conserved(self) -> bool:
        ratio = self.splay_ratio
        return ratio is not None and ratio >= _P.value("family.splay_ratio")

    def summary(self) -> str:
        if not self.core_converged:
            return (f"{self.mobile_id} on {self.target_id}: the pore modules do "
                    f"not superpose ({self.core_rmsd:.2f} A over {self.n_core} "
                    f"C-alpha, ceiling "
                    f"{_P.value('family.core_rmsd_ceiling'):.1f}); no splay "
                    f"ratio is reported because there is no fitted core to take "
                    f"one against")
        ratio = self.splay_ratio
        return (f"{self.mobile_id} on {self.target_id}: core {self.core_rmsd:.2f} A "
                f"over {self.n_core} C-alpha, blades {self.periphery_rmsd:.2f} A "
                f"over {self.n_periphery} after the same fit — splay "
                f"{ratio:.1f}x. Fitting on all {self.n_shared} shared C-alpha "
                f"instead gives {self.whole_rmsd:.2f} A")


def _domain_residues(names, species: str) -> set[int]:
    ann = load_annotations(species)
    out: set[int] = set()
    for dom in ann.domains:
        if dom.id in names and dom.start and dom.end:
            out.update(range(dom.start, dom.end + 1))
    return out


def core_residues(species: str = "human") -> set[int]:
    return _domain_residues(CORE_DOMAINS, species)


def periphery_residues(species: str = "human") -> set[int]:
    return _domain_residues(PERIPHERY_DOMAINS, species)


def _identity(structure: Structure) -> tuple[str | None, str | None]:
    """(protein, numbering) for an entry, or (None, None) if unrecognised."""
    result = identify_numbering(structure)
    if result is None:
        return None, None
    numbering = getattr(result, "numbering", None) or getattr(result, "reference", None)
    protein = getattr(result, "protein", None)
    return protein, numbering


def correspondence(mobile_numbering: str, target_numbering: str) -> dict | None:
    """A residue-number map from the mobile entry's numbering to the target's.

    ``None`` means "the same numbers mean the same residue". A dict is returned
    when the two entries are in different numberings, and it is always built
    from a real alignment — this project holds that no constant offset relates
    any two PIEZO sequences, including human to mouse within one paralogue.
    """
    if mobile_numbering == target_numbering:
        return None
    from ..core.sequence import NumberingMap
    from .paralogue import _reference

    left, right = _reference(mobile_numbering), _reference(target_numbering)
    numbering = NumberingMap.from_sequences(left["sequence"], right["sequence"],
                                            mobile_numbering, target_numbering)
    return dict(numbering.a_to_b)


def _protomer_ca(structure: Structure) -> tuple[np.ndarray, np.ndarray]:
    """C-alpha coordinates and residue numbers of the best-resolved protomer."""
    chains = well_resolved_chains(structure)
    if not chains:
        mask = structure.mask_ca()
        return structure.xyz[mask], structure.res_seq[mask]
    best, size = chains[0], -1
    for chain in chains:
        mask = structure.mask_ca() & (structure.chain == chain)
        if int(mask.sum()) > size:
            best, size = chain, int(mask.sum())
    mask = structure.mask_ca() & (structure.chain == best)
    return structure.xyz[mask], structure.res_seq[mask]


@dataclass(frozen=True)
class CoreFit:
    """The measurement **and the transform that produced it**.

    Added so the picture and the number cannot be of two different fits. The
    GUI superposes a partner on the pore module alone and needs the rotation;
    :func:`compare` needs the RMSDs. Deriving them separately would be two
    implementations of one fit, which is the disagreement ``pore_controller``
    and ``interaction_controller`` are both written to avoid.
    """

    comparison: Comparison
    rotation: np.ndarray
    translation: np.ndarray
    centroid: np.ndarray
    #: Residue numbers, in the **mobile** entry's numbering, that the fit used.
    fitted_residues: tuple[int, ...]
    #: Blade residues measured after it, same numbering.
    measured_residues: tuple[int, ...]
    #: Mobile residue number -> distance to its partner once the cores are on
    #: top of each other. This is the splay, per residue.
    deviation: dict
    #: The same distances keyed by the **target's** residue numbers. Both are
    #: kept because a caller colouring one structure by them must use that
    #: structure's own numbers, and across paralogues the two sets are
    #: different numbers for the same positions.
    deviation_target: dict


def core_fit(mobile: Structure, target: Structure, mobile_id: str = "mobile",
             target_id: str = "target") -> CoreFit | Refusal:
    """Fit one protomer on the other by the pore module and keep the transform.

    Everything :func:`compare` reports comes out of this, so a superposition
    drawn on screen is the one the table describes.
    """
    mobile_protein, mobile_numbering = _identity(mobile)
    target_protein, target_numbering = _identity(target)
    if mobile_numbering is None or target_numbering is None:
        return Refusal("one of the entries could not be identified against any "
                       "PIEZO reference; a superposition of unknown "
                       "correspondence is a picture, not a measurement")

    mapping = correspondence(mobile_numbering, target_numbering)
    mob_xyz, mob_res = _protomer_ca(mobile)
    tgt_xyz, tgt_res = _protomer_ca(target)
    target_index = {int(r): i for i, r in enumerate(tgt_res)}

    pairs = []
    for i, resi in enumerate(mob_res):
        mapped = int(resi) if mapping is None else mapping.get(int(resi))
        if mapped is None:
            continue
        j = target_index.get(int(mapped))
        if j is not None:
            pairs.append((i, j, int(resi)))
    if len(pairs) < 50:
        return Refusal(f"only {len(pairs)} residues correspond between "
                       f"{mobile_id} and {target_id}; too few to fit anything")

    # WHICH ENTRY'S NUMBERS THE CURATED RANGES ARE IN. The domains exist in
    # human and mouse PIEZO1 only, so the core has to be selected in whichever
    # of the two entries is one of those and mapped to the other through the
    # correspondence already established above. Selecting it by the mobile's
    # numbers whatever the mobile is looked right and was not: on a PIEZO2
    # mobile, "mouse" is in "mouse_piezo2", so PIEZO1's mouse ranges were
    # indexed straight into PIEZO2's numbering — the same mis-framing
    # `paralogue_identity` refuses outright, where it moves the cap from 0.35
    # to 0.85 and announces nothing.
    frame_is_mobile = is_annotated(str(mobile_numbering))
    frame = str(mobile_numbering) if frame_is_mobile else str(target_numbering)
    carried_note = ""
    if frame_is_mobile or is_annotated(frame):
        species = "mouse" if frame == "mouse" else "human"
        core = core_residues(species)
        periphery = periphery_residues(species)
    else:
        # Neither entry is PIEZO1, so the curated ranges have to be *carried*
        # into the frame rather than indexed into it — two PIEZO2 entries have
        # a pore module, and it is not at PIEZO1's numbers. Mapped through the
        # same alignment machinery the correspondence above uses.
        carried = correspondence("human", frame)
        if carried is None:
            return Refusal(
                f"neither {mobile_id} ({mobile_numbering}) nor {target_id} "
                f"({target_numbering}) is in a numbering with curated domains, "
                f"and the pore module could not be carried into either")
        core = {carried[r] for r in core_residues("human") if r in carried}
        periphery = {carried[r] for r in periphery_residues("human")
                     if r in carried}
        # Said out loud on the result: neither entry has curated domains, so
        # where the pore module *is* in these two rests on the alignment as
        # well as the fit. Below Rost's line that is the weaker of the two.
        carried_note = (f"; pore module carried into {frame} by alignment, "
                        f"not curated there")
    # `pairs` carries (mobile index, target index, mobile residue); the frame
    # decides which of the two numbers the membership test is made on.
    framed = [(i, j, (r if frame_is_mobile else int(tgt_res[j])))
              for i, j, r in pairs]
    core_idx = [(i, j) for i, j, r in framed if r in core]
    peri_idx = [(i, j) for i, j, r in framed if r in periphery]
    if len(core_idx) < 30:
        return Refusal(f"only {len(core_idx)} pore-module residues are shared; "
                       f"a core-only fit needs a core")

    mob_core = mob_xyz[[i for i, _ in core_idx]]
    tgt_core = tgt_xyz[[j for _, j in core_idx]]
    # kabsch rather than superpose so the transform itself survives; the
    # fitted coordinates and the RMSD below are what `superpose` computes from
    # the same three arrays, and a test asserts the two routes agree exactly.
    rotation, translation, centroid = kabsch(mob_core, tgt_core)
    core_rmsd = _rmsd((mob_core - centroid) @ rotation.T + translation, tgt_core)

    # The same transform, applied to the blades. Not a second fit: the whole
    # question is where the blades land once the pores are on top of each other.
    fitted_all = (mob_xyz - centroid) @ rotation.T + translation
    peri_rmsd = None
    if len(peri_idx) >= 20:
        d = fitted_all[[i for i, _ in peri_idx]] - tgt_xyz[[j for _, j in peri_idx]]
        peri_rmsd = float(np.sqrt((d * d).sum(axis=1).mean()))

    _, whole_rmsd = superpose(mob_xyz[[i for i, _, _ in pairs]],
                              tgt_xyz[[j for _, j, _ in pairs]])
    comparison = Comparison(
        mobile_id=mobile_id, target_id=target_id,
        mobile_numbering=str(mobile_numbering), target_numbering=str(target_numbering),
        n_core=len(core_idx), core_rmsd=float(core_rmsd),
        n_periphery=len(peri_idx), periphery_rmsd=peri_rmsd,
        n_shared=len(pairs), whole_rmsd=float(whole_rmsd),
        cross_paralogue=mobile_protein != target_protein,
        note=("correspondence by residue number" if mapping is None else
              f"correspondence by global alignment, {mobile_numbering} to "
              f"{target_numbering}") + carried_note)

    core_set = {int(mob_res[i]) for i, _ in core_idx}
    peri_set = {int(mob_res[i]) for i, _ in peri_idx}
    distances = {(i, j): float(np.linalg.norm(fitted_all[i] - tgt_xyz[j]))
                 for i, j, _ in pairs}
    deviation = {int(mob_res[i]): d for (i, _j), d in distances.items()}
    deviation_target = {int(tgt_res[j]): d for (_i, j), d in distances.items()}
    return CoreFit(comparison=comparison, rotation=rotation,
                   translation=translation, centroid=centroid,
                   fitted_residues=tuple(sorted(core_set)),
                   measured_residues=tuple(sorted(peri_set)),
                   deviation=deviation, deviation_target=deviation_target)


def compare(mobile: Structure, target: Structure, mobile_id: str = "mobile",
            target_id: str = "target") -> Comparison | Refusal:
    """Fit one protomer on the other by the pore module, then measure the blades.

    The numbers only; :func:`core_fit` is the same computation with the
    transform kept.
    """
    fit = core_fit(mobile, target, mobile_id, target_id)
    return fit if isinstance(fit, Refusal) else fit.comparison
