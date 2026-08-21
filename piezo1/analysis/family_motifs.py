"""A motif quoted as a family signature, and what is actually conserved that deeply.

The 2013 protozoan survey reported ``PFEW`` as absolutely conserved across the
PIEZO homologues it examined, and the phrase has been repeated since as a family
signature. The census could not find it in human PIEZO1 nor in any of its 117
representative sequences.

That is a negative about a four-letter string, and negatives about strings are
cheap to get wrong — a search that looks in the wrong place, or in one sequence,
proves nothing. So this module does two things:

1. **Searches every reference this project holds**, which is now ten proteins
   spanning human to *Dictyostelium*, and reports the count per protein rather
   than a total. A total of zero is indistinguishable from a broken search; ten
   explicit zeros beside a positive control are not.
2. **Carries a positive control.** :func:`motif_scan` is run on a motif taken
   from human PIEZO1's own sequence before the absent one is believed. A search
   that cannot find a string known to be there is not evidence that another
   string is absent.

What *is* conserved to that depth is then measured rather than asserted, from
the census's whole-family alignment layer: the windows where the family track
stays high across the deepest comparison available.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..core.family import load_constraint
from ..core.numbering_check import PROTEIN_NAMES, REFERENCES, reference_entry

__all__ = ["MotifHit", "MotifScan", "ConservedWindow", "QUOTED_MOTIFS",
           "motif_scan", "control_motif", "deep_windows"]

#: Motifs the literature has called family signatures. Each is checked, and the
#: result is a count per protein rather than a verdict.
QUOTED_MOTIFS = {
    "PFEW": ("reported absolutely conserved across protozoan PIEZO homologues "
             "in the 2013 survey and repeated since as a family signature"),
}


@dataclass(frozen=True)
class MotifHit:
    reference: str
    protein: str
    length: int
    positions: tuple

    @property
    def count(self) -> int:
        return len(self.positions)


@dataclass(frozen=True)
class MotifScan:
    """One motif against every reference sequence this project holds."""

    motif: str
    hits: tuple
    note: str = ""

    @property
    def total(self) -> int:
        return sum(h.count for h in self.hits)

    @property
    def n_proteins(self) -> int:
        return len(self.hits)

    @property
    def present_in(self) -> tuple:
        return tuple(h.protein for h in self.hits if h.count)

    def summary(self) -> str:
        if not self.total:
            return (f"{self.motif} does not occur in any of the "
                    f"{self.n_proteins} PIEZO reference sequences this project "
                    f"holds ({sum(h.length for h in self.hits):,} residues in "
                    f"total)")
        # Three references are all called PIEZO1 (human, mouse, rat), so the
        # reference name is what distinguishes them, not the protein name.
        where = ", ".join(f"{h.reference} x{h.count}"
                          for h in self.hits if h.count)
        return (f"{self.motif} occurs {self.total} times, in "
                f"{len(self.present_in)} of {self.n_proteins} references: {where}")


@dataclass(frozen=True)
class ConservedWindow:
    """A stretch of human PIEZO1 that stays conserved at family depth."""

    start: int
    end: int
    mean: float
    sequence: str
    domain: str

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def summary(self) -> str:
        return (f"{self.start}-{self.end} ({self.domain}): {self.sequence} "
                f"[mean family constraint {self.mean:.2f}]")


def motif_scan(motif: str = "PFEW") -> MotifScan:
    """Every occurrence of a literal motif in every PIEZO reference we hold."""
    motif = motif.upper()
    hits = []
    for name in REFERENCES:
        sequence = reference_entry(name)["sequence"].upper()
        positions, start = [], sequence.find(motif)
        while start != -1:
            positions.append(start + 1)               # 1-based residue number
            start = sequence.find(motif, start + 1)
        hits.append(MotifHit(reference=name, protein=PROTEIN_NAMES.get(name, name),
                             length=len(sequence), positions=tuple(positions)))
    return MotifScan(motif=motif, hits=tuple(hits),
                     note=QUOTED_MOTIFS.get(motif, "not a quoted motif"))


def control_motif(residue: int = 2456, length: int = 4) -> MotifScan:
    """A motif taken from human PIEZO1's own sequence, so the search can succeed.

    Without this the absence of ``PFEW`` is unfalsifiable: a scan that returns
    zero everywhere looks identical whether the motif is absent or the reader is
    broken. Defaults to the four residues around R2456 because that is a
    position the rest of this subsystem is about.
    """
    sequence = reference_entry("human")["sequence"]
    start = max(0, residue - 1)
    return motif_scan(sequence[start:start + length])


def deep_windows(n: int = 5, width: int | None = None,
                 track: str = "family_jsd") -> list[ConservedWindow]:
    """The stretches of human PIEZO1 most conserved at whole-family depth.

    ``family_jsd`` is the census's deepest layer — columns where human PIEZO1 is
    being compared with a plant and an amoeba — so it is the right track for
    "what has survived since the root of the eukaryotes", and the wrong one for
    anything about vertebrates, where it is far too coarse.

    Returned windows do not overlap: the top window's neighbourhood is excluded
    before the next is taken, or the answer would be five views of one peak.
    """
    from ..core.annotations import load_annotations
    from ..parameters import PARAMETERS as _P

    width = int(_P.value("family.motif_window")) if width is None else int(width)
    constraint = load_constraint("PIEZO1", track)
    values = constraint.values.copy()
    kernel = np.ones(width) / width
    smooth = np.convolve(np.nan_to_num(values, nan=0.0), kernel, mode="same")
    smooth[np.isnan(values)] = -np.inf
    annotations = load_annotations("human")

    chosen = []
    for _ in range(n):
        centre = int(np.argmax(smooth))
        if not np.isfinite(smooth[centre]):
            break
        start = max(1, centre + 1 - width // 2)
        end = min(constraint.length, start + width - 1)
        window = values[start - 1:end]
        domain = annotations.domain_at(start + width // 2)
        chosen.append(ConservedWindow(
            start=start, end=end, mean=float(np.nanmean(window)),
            sequence=constraint.sequence[start - 1:end],
            domain=(domain.name if domain else "unassigned")))
        lo = max(0, centre - width)
        hi = min(smooth.size, centre + width + 1)
        smooth[lo:hi] = -np.inf
    return chosen
