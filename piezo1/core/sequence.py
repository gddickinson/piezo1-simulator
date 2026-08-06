"""Sequences and cross-species residue numbering.

This module exists because of one recurring, silent source of error in the
PIEZO literature: **mechanistic papers number residues by mouse Piezo1 (2547 aa)
while clinical papers number by human PIEZO1 (2521 aa)**, and the offset between
them is not constant — insertions and deletions shift it along the chain.

Quoting a "Yoda1 pocket residue A2094" onto a human structure without going
through an alignment is how a viewer ends up highlighting the wrong helix. So
every conversion in this project goes through :class:`NumberingMap`, which is
built from an actual global alignment and cached.
"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from ..config import HUMAN_ACC, MOUSE_ACC, RESOURCE_DIR

__all__ = ["align_global", "NumberingMap", "load_numbering_map",
           "human_sequence", "mouse_sequence", "SPECIES_ACCESSION"]

SPECIES_ACCESSION = {"human": HUMAN_ACC, "mouse": MOUSE_ACC}

_BLOSUM62_FALLBACK_MATCH = 2.0
_BLOSUM62_FALLBACK_MISMATCH = -1.0


def _uniprot_resource(species: str) -> dict:
    path = RESOURCE_DIR / f"uniprot_{species}.json"
    if not path.exists():
        raise FileNotFoundError(
            f"{path} missing — run scripts/build_uniprot_annotations.py"
        )
    return json.loads(path.read_text())


def human_sequence() -> str:
    return _uniprot_resource("human")["sequence"]


def mouse_sequence() -> str:
    return _uniprot_resource("mouse")["sequence"]


# --------------------------------------------------------------------------
# Alignment
# --------------------------------------------------------------------------

def align_global(seq_a: str, seq_b: str) -> tuple[str, str]:
    """Global alignment of two protein sequences, returned as gapped strings.

    Uses Biopython's ``PairwiseAligner`` with BLOSUM62 when available; falls
    back to a simple match/mismatch scheme otherwise.  For two ~2500-residue
    PIEZO1 orthologs that are ~80% identical this is fast and unambiguous.
    """
    try:
        from Bio import Align
        from Bio.Align import substitution_matrices

        aligner = Align.PairwiseAligner()
        aligner.mode = "global"
        try:
            aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
            aligner.open_gap_score = -11.0
            aligner.extend_gap_score = -1.0
        except Exception:  # pragma: no cover - matrix data missing
            aligner.match_score = _BLOSUM62_FALLBACK_MATCH
            aligner.mismatch_score = _BLOSUM62_FALLBACK_MISMATCH
            aligner.open_gap_score = -10.0
            aligner.extend_gap_score = -0.5
        aln = aligner.align(seq_a, seq_b)[0]
        return str(aln[0]), str(aln[1])
    except ImportError as exc:  # pragma: no cover
        raise ImportError("Biopython is required for sequence alignment") from exc


# --------------------------------------------------------------------------
# Numbering map
# --------------------------------------------------------------------------

@dataclass
class NumberingMap:
    """Bidirectional residue-number mapping between two orthologs.

    ``a_to_b[i]`` is the residue number in sequence B aligned to residue ``i``
    of sequence A, or ``None`` where A has an insertion relative to B.
    """

    species_a: str
    species_b: str
    a_to_b: dict[int, int | None]
    b_to_a: dict[int, int | None]
    identity: float
    aligned_length: int
    gaps: int

    # ------------------------------------------------------------ conversion

    def to_b(self, residue: int) -> int | None:
        return self.a_to_b.get(int(residue))

    def to_a(self, residue: int) -> int | None:
        return self.b_to_a.get(int(residue))

    def convert(self, residue: int, source: str) -> int | None:
        """Convert ``residue`` from ``source`` species to the other species."""
        if source == self.species_a:
            return self.to_b(residue)
        if source == self.species_b:
            return self.to_a(residue)
        raise ValueError(f"unknown species {source!r}; expected "
                         f"{self.species_a!r} or {self.species_b!r}")

    def convert_range(self, start: int, end: int, source: str) -> tuple[int | None, int | None]:
        """Convert a residue range, walking inwards past unaligned ends."""
        fwd = self.a_to_b if source == self.species_a else self.b_to_a
        lo = next((fwd[i] for i in range(start, end + 1) if fwd.get(i) is not None), None)
        hi = next((fwd[i] for i in range(end, start - 1, -1) if fwd.get(i) is not None), None)
        return lo, hi

    def offset_at(self, residue: int, source: str) -> int | None:
        """Signed numbering offset at one position, for diagnostics."""
        other = self.convert(residue, source)
        return None if other is None else other - residue

    def offset_profile(self, source: str, step: int = 1) -> list[tuple[int, int | None]]:
        """``[(residue, offset), ...]`` — shows where the offset changes."""
        fwd = self.a_to_b if source == self.species_a else self.b_to_a
        return [(i, (v - i) if v is not None else None)
                for i, v in sorted(fwd.items())[::step]]

    # ---------------------------------------------------------------- build

    @classmethod
    def from_sequences(cls, seq_a: str, seq_b: str,
                       species_a: str, species_b: str) -> "NumberingMap":
        ga, gb = align_global(seq_a, seq_b)
        a_to_b: dict[int, int | None] = {}
        b_to_a: dict[int, int | None] = {}
        ia = ib = 0
        matches = 0
        gaps = 0
        for ca, cb in zip(ga, gb):
            if ca != "-" and cb != "-":
                ia += 1
                ib += 1
                a_to_b[ia] = ib
                b_to_a[ib] = ia
                matches += ca == cb
            elif ca != "-":
                ia += 1
                a_to_b[ia] = None
                gaps += 1
            elif cb != "-":
                ib += 1
                b_to_a[ib] = None
                gaps += 1
        aligned = min(len(seq_a), len(seq_b))
        return cls(species_a=species_a, species_b=species_b,
                   a_to_b=a_to_b, b_to_a=b_to_a,
                   identity=matches / max(aligned, 1),
                   aligned_length=len(ga), gaps=gaps)

    # ----------------------------------------------------------- persistence

    def to_json(self) -> dict:
        return {
            "species_a": self.species_a,
            "species_b": self.species_b,
            "identity": self.identity,
            "aligned_length": self.aligned_length,
            "gaps": self.gaps,
            # Stored as arrays indexed from residue 1; null marks an insertion.
            "a_to_b": [self.a_to_b.get(i) for i in range(1, max(self.a_to_b) + 1)],
            "b_to_a": [self.b_to_a.get(i) for i in range(1, max(self.b_to_a) + 1)],
        }

    @classmethod
    def from_json(cls, data: dict) -> "NumberingMap":
        return cls(
            species_a=data["species_a"], species_b=data["species_b"],
            a_to_b={i + 1: v for i, v in enumerate(data["a_to_b"])},
            b_to_a={i + 1: v for i, v in enumerate(data["b_to_a"])},
            identity=data["identity"], aligned_length=data["aligned_length"],
            gaps=data["gaps"],
        )


_CACHE_NAME = "numbering_human_mouse.json"
_MAP: NumberingMap | None = None


def load_numbering_map(rebuild: bool = False) -> NumberingMap:
    """Load (or build and cache) the human↔mouse PIEZO1 numbering map."""
    global _MAP
    if _MAP is not None and not rebuild:
        return _MAP
    path: Path = RESOURCE_DIR / _CACHE_NAME
    if path.exists() and not rebuild:
        _MAP = NumberingMap.from_json(json.loads(path.read_text()))
        return _MAP
    _MAP = NumberingMap.from_sequences(
        human_sequence(), mouse_sequence(), "human", "mouse"
    )
    path.write_text(json.dumps(_MAP.to_json()))
    return _MAP


def human_to_mouse(residue: int) -> int | None:
    return load_numbering_map().to_b(residue)


def mouse_to_human(residue: int) -> int | None:
    return load_numbering_map().to_a(residue)
