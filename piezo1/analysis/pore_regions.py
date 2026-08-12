"""*Where* along the pore the narrow points are — not whether it conducts.

:mod:`piezo1.analysis.hydration` answers whether a pore conducts, and answers
it with a single number: the minimum radius over the whole axial profile, taken
with the hydrophobicity score. That is the right criterion for a
one-dimensional conduction model, and it is unreadable on its own, because
"sterically occluded" invites exactly one reading — *the gate is shut* — and
across this project's catalogue that reading is wrong in every case.

Measured over the 19 deposited PIEZO1 entries, with the transmembrane gate
located from the curated ``hydrophobic_gate`` residues rather than assumed —
18 locate it, and 6LQI is refused because it is deposited in a splice
isoform's numbering:

- the global bottleneck lies **beyond the gate** — at the cytoplasmic
  constriction in 16, above it in the cap in 2;
- it lies **at the transmembrane gate in none** of them;
- the narrowest point *within* the gate is 2.4–4.7 Å, at or above the 1.5 Å
  water radius the steric criterion uses, in every entry.

That is not a defect of the wetting heuristic — it is a property of the
channel, and it is the published one. Liu et al. 2025 (Neuron 113:590–604)
report that in the intermediate-**open** S2472E structure (8IXO) the
transmembrane gate dilates while "the vertical constriction neck remains
closed", because "the lateral portals rather than the constriction neck
represent the major ion-permeation routes". An axial profile contains no
lateral portals, so it must pass through a constriction the real channel goes
around.

The coordinates carry the dilation, on the paper's own measure: the V2476
side-chain diagonal is **7.7 Å on 7WLT and 14.2 Å on 8IXO**, against the
7 Å → 14 Å the paper states. What our *pore radius* makes of it is weaker —
8IXO's gate is 3.52 Å against 2.4–3.2 Å for the well-resolved curved entries,
but the two entries with wider gates still (7WLU 4.67 Å, 3JAC 4.34 Å) are the
two worst-resolved, so radius at the gate is confounded with resolution and
does not separate the states on its own. What is unambiguous is the refusal:
8IXO's lining clears the Rao cutoff at 0.31, its gate is 3.52 Å, and it is
turned down on a 0.98 Å neck at E2537 — the constriction the paper says is
bypassed.

**Three numbers, not one.** The two constrictions flanking the gate are within
0.02 Å of each other on 8IXO, so which of them is the global minimum flips with
the frame the structure is loaded in. Reporting the narrowest point in each
region instead is parameter-free, does not flip, and says the whole thing.

This module changes no verdict; a test enforces that.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.annotations import Annotations

__all__ = ["Constriction", "Bottleneck", "GATE_GROUP", "gate_numbering",
           "gate_mask", "describe_bottleneck"]

#: The curated group naming the transmembrane gate: human I2447/V2450/F2454,
#: mouse I2473/V2476/F2480, sourced to Yang et al. eLife 2025. Used rather than
#: a residue range written down here, so the gate has one definition in this
#: project and it is the one with a citation attached.
GATE_GROUP = "hydrophobic_gate"


@dataclass(frozen=True)
class Constriction:
    """The narrowest slice of one region of the profile."""

    radius: float                       # Angstrom
    z: float                            # Angstrom along the conduction axis
    lining: tuple[int, ...] = ()
    lining_names: tuple[str, ...] = ()

    def label(self) -> str:
        """The lining residues as ``VAL2476`` labels.

        Safe to zip only because :class:`piezo1.structure.pore.PoreSlice` keeps
        the two tuples parallel — they were separately sorted sets until
        Round 84c, which labelled 8YEZ's worst dewetted residue GLU2510 when
        it is PRO2510.
        """
        return "/".join(f"{n}{r}" for n, r
                        in zip(self.lining_names, self.lining)) or "?"

    def text(self) -> str:
        return f"{self.radius / 10:.3f} nm at {self.label()}"


@dataclass(frozen=True)
class Bottleneck:
    """The narrow points of an axial pore profile, by region."""

    gate: Constriction | None = None      # within the curated gate
    above: Constriction | None = None     # extracellular of the gate
    below: Constriction | None = None     # cytoplasmic of the gate
    numbering: str | None = None
    reason: str = ""                      # why a region is missing, when it is

    @property
    def narrowest_region(self) -> str:
        """Which region holds the global minimum: gate, above, below, unknown.

        The three regions partition the profile, so this really is the global
        minimum's region and not merely the smallest of three samples.
        """
        found = {k: c for k, c in (("gate", self.gate), ("above", self.above),
                                   ("below", self.below)) if c is not None}
        if not found:
            return "unknown"
        return min(found, key=lambda k: found[k].radius)

    @property
    def blocked_beyond_the_gate(self) -> bool:
        """True when the narrowest point is not in the gate — the usual case."""
        return self.narrowest_region in ("above", "below")

    def sentence(self) -> str:
        """One line placing every narrow point, gate first.

        The gate goes first deliberately. It is the constriction a reader means
        by "is the channel open", and quoting only the global minimum — which
        is somewhere else in every entry measured — is what made the refusal
        read as a statement about the gate.
        """
        if self.gate is None:
            return self.reason or "gate not located"
        parts = [f"transmembrane gate {self.gate.text()}"]
        if self.above is not None:
            parts.append(f"above it {self.above.text()}")
        if self.below is not None:
            parts.append(f"below it, towards the cytosol, {self.below.text()}")
        text = "; ".join(parts)
        if self.numbering:
            text += f" ({self.numbering} numbering)"
        return text


def gate_numbering(structure) -> str | None:
    """Which numbering the file is in, or ``None`` if it is not PIEZO1's.

    Read off the coordinates by :func:`piezo1.core.numbering_check.
    identify_numbering`, which scores an entry's own residue names against all
    six reference sequences — not taken from the registry, and not decided
    here. The first version of this function did decide it here, by checking
    the three gate residues' names at the human and mouse numbers, and read
    **mouse PIEZO2 (6KG7) as human PIEZO1**: three positions is not enough
    evidence to separate two proteins when the expected residues are Ile, Val
    and Phe, which any transmembrane helix is full of.

    A ``None`` here is a refusal, not a default to human, and it is also
    returned for a PIEZO1 entry the reference cannot read by number: 6LQI is
    deposited in the Piezo1.1 isoform's own numbering, so its residue 2476 is
    not the gate — measured, its "V2476" side chains sit 31 Å apart against
    7.7 Å on 7WLT.
    """
    from ..core.numbering_check import piezo1_numbering

    return piezo1_numbering(structure)


def gate_mask(structure, profile, numbering: str | None = None) -> np.ndarray:
    """Which slices of ``profile`` are lined by a curated gate residue.

    ``numbering`` is passed in by :func:`describe_bottleneck`, which has
    already measured it; left out, it is measured here.
    """
    if numbering is None:
        numbering = gate_numbering(structure)
    mask = np.zeros(len(profile.slices), dtype=bool)
    if numbering is None:
        return mask
    numbers = {e[numbering] for e in _gate_detail() if e.get(numbering)}
    for i, sl in enumerate(profile.slices):
        mask[i] = bool(numbers.intersection(sl.lining))
    return mask


def describe_bottleneck(structure, profile) -> Bottleneck:
    """Locate the narrowest point of each region relative to the gate.

    Which side is cytoplasmic is **measured**, not read off a residue number:
    the gate's own slices give a span along the conduction axis and
    :func:`piezo1.physics.pore_charge.cytosolic_end` gives the direction. A
    residue-number rule would have called the beam residue 1412 — which lines
    the narrowest slice on four entries — extracellular, because it is numbered
    below the pore.
    """
    from ..physics.pore_charge import cytosolic_end

    numbering = gate_numbering(structure)
    mask = gate_mask(structure, profile, numbering)
    if not mask.any():
        return Bottleneck(numbering=numbering,
                          reason="no slice is lined by a curated gate residue "
                                 "— not a PIEZO1 entry, or the gate is not "
                                 "resolved")
    if profile.axis is None:
        return Bottleneck(gate=_narrowest(profile, mask), numbering=numbering,
                          reason="profile carries no axis, so the two sides of "
                                 "the gate cannot be told apart")

    # The gate is a *stretch* of the pore, not just the slices a gate residue
    # happens to touch: a narrow slice between two gate residues belongs to the
    # gate, and leaving it out of all three regions would mean the region
    # holding the smallest radius was not always the global minimum's.
    low, high = float(profile.z[mask].min()), float(profile.z[mask].max())
    span = (profile.z >= low) & (profile.z <= high)
    try:
        cytosolic_low = cytosolic_end(structure, profile.axis) == 0
    except ValueError as exc:
        return Bottleneck(gate=_narrowest(profile, span), numbering=numbering,
                          reason=f"cytosolic end not measurable: {exc}")

    lower, upper = profile.z < low, profile.z > high
    cytosolic, extracellular = (lower, upper) if cytosolic_low else (upper, lower)
    return Bottleneck(gate=_narrowest(profile, span),
                      above=_narrowest(profile, extracellular),
                      below=_narrowest(profile, cytosolic),
                      numbering=numbering)


def _narrowest(profile, mask: np.ndarray) -> Constriction | None:
    """The narrowest slice among ``mask``, or None when it selects nothing."""
    if not mask.any():
        return None
    index = int(np.flatnonzero(mask)[np.argmin(profile.radius[mask])])
    sl = profile.slices[index]
    return Constriction(radius=float(profile.radius[index]),
                        z=float(profile.z[index]), lining=tuple(sl.lining),
                        lining_names=tuple(sl.lining_names))


def _gate_detail() -> tuple[dict, ...]:
    group = Annotations().group(GATE_GROUP)
    return group.detail if group is not None else ()
