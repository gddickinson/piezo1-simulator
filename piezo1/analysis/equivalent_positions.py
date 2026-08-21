"""Are PIEZO1 R2456 and PIEZO2 R2686 the same place, or only the same column?

The census found that the fourteen pathogenic positions inside the pore module
occupy only twelve alignment columns: PIEZO1 R2456 (hereditary xerocytosis) and
PIEZO2 R2686 (Gordon syndrome) fall in one, PIEZO1 R2488 and PIEZO2 R2718 in
another. Half a billion years apart, and disease finds the same place in both.

**That is a claim about an alignment.** An alignment always produces an answer;
whether two aligned residues occupy the same position in space is a different
question, and it is the one this project can settle. So the test here has three
parts, in order, and the order matters:

1. **Does an independent alignment agree?** The census's column comes from its
   117-sequence family alignment. This project builds its own pairwise global
   alignment between the two paralogues for other purposes entirely
   (:func:`piezo1.analysis.paralogue.paralogue_map`). If that map does not pair
   the same two residues, the structural test should not be run at all.
2. **Do the C-alpha land on top of each other** once the two pore modules are
   superposed?
3. **Is that distance small compared with what any aligned pair gives?** Two
   superposed pore modules put *everything* near everything; a 3 A separation is
   only evidence if the median aligned core pair is much further apart. The
   control is the whole distribution, and it is what makes the answer capable of
   being "no".

Step 3 is the calibration. Without it this module would report a small number
for every pair it was handed and would have no way to say a claim was wrong.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.family import EquivalentPosition, load_family_findings
from ..core.numbering_check import identify_numbering
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from ..structure.superpose import superpose
from .core_periphery import _protomer_ca, core_residues, correspondence

__all__ = ["PositionPair", "EquivalenceReport", "locate", "alignment_agrees",
           "map_position"]


@dataclass(frozen=True)
class PositionPair:
    """One claimed-equivalent pair, measured on coordinates."""

    label: str
    element: str
    piezo1_human: int
    piezo2_human: int
    piezo1_resi: int | None
    piezo2_resi: int | None
    alignment_agrees: bool
    distance: float | None
    percentile: float | None
    nearest_resi: int | None = None
    nearest_distance: float | None = None
    piezo1_disease: str = ""
    piezo2_disease: str = ""
    note: str = ""

    @property
    def is_nearest(self) -> bool:
        """Is the claimed partner the closest residue of the other paralogue?

        The decisive form of the test. The percentile control below establishes
        that the whole pore module superposes, so *any* aligned core pair comes
        out close; what distinguishes a correct correspondence from a register
        error one residue along is whether the claimed partner is the nearest
        thing there.
        """
        return (self.nearest_resi is not None
                and self.nearest_resi == self.piezo2_resi)

    @property
    def register_offset(self) -> int | None:
        """How many residues away the nearest partner is from the claimed one.

        Zero is exact. A cross-paralogue superposition of two ~3.5 A structures
        cannot resolve better than about one residue along a helix, so an
        offset of 1 confirms the correspondence at the resolution available and
        an offset of 5 does not.
        """
        if self.nearest_resi is None or self.piezo2_resi is None:
            return None
        return self.nearest_resi - self.piezo2_resi

    @property
    def same_place(self) -> bool:
        """Within the cutoff, alignment agreeing, and nearest to within one."""
        if self.distance is None or not self.alignment_agrees:
            return False
        offset = self.register_offset
        return (self.distance <= _P.value("family.equivalent_ca_cutoff")
                and offset is not None and abs(offset) <= 1)

    def summary(self) -> str:
        if self.distance is None:
            return f"{self.label}: not measurable — {self.note}"
        agree = "agrees" if self.alignment_agrees else "DISAGREES"
        pct = ("" if self.percentile is None
               else f", {self.percentile:.0%} of aligned core pairs are closer")
        if self.nearest_resi is None:
            near = ""
        elif self.is_nearest:
            near = "; it is the nearest residue of the other paralogue"
        else:
            near = (f"; the nearest residue of the other paralogue is "
                    f"{self.nearest_resi} at {self.nearest_distance:.2f} A "
                    f"({self.register_offset:+d} in register)")
        return (f"{self.label} ({self.element}): C-alpha {self.distance:.2f} A "
                f"apart after a pore-module superposition{pct}{near}; this "
                f"project's own alignment {agree} with the census column")


@dataclass(frozen=True)
class EquivalenceReport:
    """Every claimed pair on one structure pair, with the control distribution."""

    piezo1_id: str
    piezo2_id: str
    pairs: tuple
    n_control: int
    control_median: float | None
    control_p10: float | None
    core_rmsd: float | None
    note: str = ""

    @property
    def all_agree(self) -> bool:
        return bool(self.pairs) and all(p.alignment_agrees for p in self.pairs)

    @property
    def all_same_place(self) -> bool:
        return bool(self.pairs) and all(p.same_place for p in self.pairs)

    @property
    def verdict(self) -> str:
        if not self.pairs:
            return self.note
        measured = [p for p in self.pairs if p.distance is not None]
        if not measured:
            return "no claimed pair could be measured on these coordinates"
        exact = sum(p.is_nearest for p in measured)
        within_one = sum(p.same_place for p in measured)
        agree = sum(p.alignment_agrees for p in measured)
        worst = max(p.distance for p in measured)
        return (f"{agree} of {len(measured)} pairs are confirmed by this "
                f"project's own alignment; {exact} put the claimed partner "
                f"nearest in space exactly and {within_one} to within one "
                f"residue, the furthest C-alpha at {worst:.2f} A. The whole "
                f"pore module superposes here — the median aligned core pair "
                f"is {self.control_median:.2f} A apart — so a small distance "
                f"alone is not the evidence; the register is. One residue is "
                f"also the resolution a cross-paralogue fit at "
                f"{self.core_rmsd:.1f} A can claim")

    def summary(self) -> str:
        if not self.pairs:
            return f"{self.piezo1_id} vs {self.piezo2_id}: {self.note}"
        head = (f"{self.piezo1_id} vs {self.piezo2_id}: pore modules superposed "
                f"at {self.core_rmsd:.2f} A; the median aligned core pair is "
                f"{self.control_median:.2f} A apart and the closest tenth are "
                f"within {self.control_p10:.2f} A")
        return head + "\n" + "\n".join("  " + p.summary() for p in self.pairs)


def map_position(residue: int, from_numbering: str, to_numbering: str) -> int | None:
    """One residue number carried between two PIEZO references by alignment."""
    if from_numbering == to_numbering:
        return residue
    mapping = correspondence(from_numbering, to_numbering)
    return None if mapping is None else mapping.get(residue)


def alignment_agrees(pair: EquivalentPosition) -> bool:
    """Does this project's own PIEZO1-PIEZO2 alignment pair the same residues?

    The census's column and this project's pairwise map are built from different
    sequence sets by different algorithms. Agreement is the precondition for the
    structural test; disagreement would mean the two projects are talking about
    different residues and no superposition could settle it.
    """
    mapped = map_position(pair.piezo1, "human", "human_piezo2")
    return mapped == pair.piezo2


def locate(piezo1_structure: Structure, piezo2_structure: Structure,
           piezo1_id: str = "PIEZO1", piezo2_id: str = "PIEZO2",
           seed: int = 0) -> EquivalenceReport:
    """Measure every claimed-equivalent pair on one PIEZO1/PIEZO2 structure pair."""
    findings = load_family_findings()
    claimed = findings.equivalent

    ident1 = identify_numbering(piezo1_structure)
    ident2 = identify_numbering(piezo2_structure)
    if ident1 is None or ident2 is None:
        return EquivalenceReport(piezo1_id, piezo2_id, (), 0, None, None, None,
                                 note="one entry could not be identified")
    if ident1.protein != "PIEZO1" or ident2.protein != "PIEZO2":
        return EquivalenceReport(
            piezo1_id, piezo2_id, (), 0, None, None, None,
            note=(f"this test needs one PIEZO1 and one PIEZO2 entry; got "
                  f"{ident1.protein} and {ident2.protein}"))

    n1, n2 = ident1.reference, ident2.reference
    mapping = correspondence(n1, n2)
    xyz1, res1 = _protomer_ca(piezo1_structure)
    xyz2, res2 = _protomer_ca(piezo2_structure)
    index2 = {int(r): i for i, r in enumerate(res2)}

    pairs_idx = []
    for i, resi in enumerate(res1):
        mapped = int(resi) if mapping is None else mapping.get(int(resi))
        if mapped is None:
            continue
        j = index2.get(int(mapped))
        if j is not None:
            pairs_idx.append((i, j, int(resi)))

    species = "mouse" if "mouse" in str(n1) else "human"
    core = core_residues(species)
    core_idx = [(i, j) for i, j, r in pairs_idx if r in core]
    if len(core_idx) < 30:
        return EquivalenceReport(piezo1_id, piezo2_id, (), 0, None, None, None,
                                 note="too few shared pore-module residues to fit")

    mob = xyz1[[i for i, _ in core_idx]]
    tgt = xyz2[[j for _, j in core_idx]]
    _, core_rmsd = superpose(mob, tgt)
    fitted, _ = superpose(mob, tgt, apply_to=xyz1)

    # The control: every aligned core pair's C-alpha separation after the same
    # fit. A claimed pair has to beat this distribution, not merely be small.
    d = fitted[[i for i, _ in core_idx]] - xyz2[[j for _, j in core_idx]]
    control = np.sqrt((d * d).sum(axis=1))

    index1 = {int(r): i for i, r in enumerate(res1)}
    measured = []
    for claim in claimed:
        r1 = map_position(claim.piezo1, "human", n1)
        r2 = map_position(claim.piezo2, "human_piezo2", n2)
        agrees = alignment_agrees(claim)
        i = index1.get(r1) if r1 else None
        j = index2.get(r2) if r2 else None
        if i is None or j is None:
            measured.append(PositionPair(
                label=claim.label, element=claim.element,
                piezo1_human=claim.piezo1, piezo2_human=claim.piezo2,
                piezo1_resi=r1, piezo2_resi=r2, alignment_agrees=agrees,
                distance=None, percentile=None,
                piezo1_disease=claim.piezo1_disease,
                piezo2_disease=claim.piezo2_disease,
                note=("one of the two residues is not resolved in its entry"
                      if r1 and r2 else
                      "a residue number could not be carried between the "
                      "two references")))
            continue
        vec = fitted[i] - xyz2[j]
        dist = float(np.sqrt((vec * vec).sum()))
        # Which residue of the other paralogue is actually closest. A register
        # error of one residue along a helix moves a C-alpha ~1.5 A, which the
        # 5 A cutoff would not catch and this does.
        deltas = xyz2 - fitted[i]
        all_d = np.sqrt((deltas * deltas).sum(axis=1))
        k = int(np.argmin(all_d))
        measured.append(PositionPair(
            label=claim.label, element=claim.element,
            piezo1_human=claim.piezo1, piezo2_human=claim.piezo2,
            piezo1_resi=r1, piezo2_resi=r2, alignment_agrees=agrees,
            distance=dist,
            percentile=float((control < dist).mean()),
            nearest_resi=int(res2[k]), nearest_distance=float(all_d[k]),
            piezo1_disease=claim.piezo1_disease,
            piezo2_disease=claim.piezo2_disease,
            note=claim.note))

    return EquivalenceReport(
        piezo1_id=piezo1_id, piezo2_id=piezo2_id, pairs=tuple(measured),
        n_control=int(control.size), control_median=float(np.median(control)),
        control_p10=float(np.percentile(control, 10)),
        core_rmsd=float(core_rmsd),
        note=("distances are C-alpha to C-alpha after a superposition using "
              "shared pore-module residues only"))
