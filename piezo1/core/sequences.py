"""Named sequences and their comparison — the data behind the sequence viewer.

Qt-free on purpose, like the measurement logic, so it can be tested without a
display and driven from a notebook.

Three kinds of sequence appear here and they are **not** interchangeable:

* **UniProt** — the reference protein, numbered 1..N. This is what variant
  numbering means.
* **Structure** — what a deposited model actually resolves, which is a subset
  with gaps. Treating it as the reference silently renumbers everything.
* **Translated CDS** — the protein implied by the downloaded coding sequence.
  Real DNA, so codons and silent variants are representable.

Keeping the distinction explicit is the point: a sequence viewer that shows
"the sequence" without saying which one is a numbering bug waiting to happen.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..config import SEQUENCE_DIR
from .sequence import align_global

__all__ = ["NamedSequence", "SequenceComparison", "translate", "CODON_TABLE",
           "load_named_sequences", "compare_sequences", "read_fasta"]

_BASES = "TCAG"
_AAS = ("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG")
#: The standard genetic code, built rather than typed to avoid transcription
#: slips in a 64-entry table.
CODON_TABLE: dict[str, str] = {
    b1 + b2 + b3: _AAS[i]
    for i, (b1, b2, b3) in enumerate(
        (a, b, c) for a in _BASES for b in _BASES for c in _BASES)
}


def translate(dna: str) -> str:
    """Translate a coding sequence. Unknown codons become ``X``."""
    dna = dna.upper().replace("U", "T")
    return "".join(CODON_TABLE.get(dna[i:i + 3], "X")
                   for i in range(0, len(dna) - 2, 3))


def read_fasta(path) -> str:
    """Sequence body of a single-record FASTA, header discarded."""
    return "".join(line.strip() for line in path.read_text().splitlines()
                   if not line.startswith(">"))


@dataclass
class NamedSequence:
    """One sequence with the numbering it is expressed in."""

    key: str
    label: str
    letters: str
    kind: str = "protein"            # protein | dna
    numbering: str = "human"         # which residue-numbering system
    #: Residue number of each letter. Explicit because a structure sequence is
    #: not 1..N — it has gaps, and assuming otherwise renumbers every variant.
    positions: list[int] = field(default_factory=list)
    dna: str = ""                    # coding sequence, when one is known
    source: str = ""
    note: str = ""

    def __post_init__(self) -> None:
        if not self.positions:
            self.positions = list(range(1, len(self.letters) + 1))

    def __len__(self) -> int:
        return len(self.letters)

    @property
    def has_gaps(self) -> bool:
        return any(b - a != 1 for a, b in zip(self.positions, self.positions[1:]))

    def at(self, residue: int) -> str | None:
        try:
            return self.letters[self.positions.index(residue)]
        except ValueError:
            return None

    def index_of(self, residue: int) -> int | None:
        try:
            return self.positions.index(residue)
        except ValueError:
            return None

    def codon(self, residue: int) -> str:
        """The three bases coding this residue, or empty if no DNA is loaded.

        Uses the residue's **ordinal position in this sequence**, not its
        residue number: a structure sequence starting at residue 576 does not
        start at codon 576 of the transcript.
        """
        index = self.index_of(residue)
        if index is None or not self.dna:
            return ""
        return self.dna[index * 3:index * 3 + 3]

    def segment(self, start: int, end: int) -> str:
        lo, hi = sorted((start, end))
        return "".join(letter for letter, pos in zip(self.letters, self.positions)
                       if lo <= pos <= hi)


@dataclass
class SequenceComparison:
    """A pairwise alignment and what it found."""

    a: NamedSequence
    b: NamedSequence
    aligned_a: str
    aligned_b: str
    identity: float
    n_identical: int
    n_mismatch: int
    n_gap: int
    #: (column, residue in a or None, residue in b or None, letter a, letter b)
    differences: list[tuple[int, int | None, int | None, str, str]] = \
        field(default_factory=list)

    @property
    def length(self) -> int:
        return len(self.aligned_a)

    def summary(self) -> str:
        return (f"{self.a.label} vs {self.b.label}: {self.identity * 100:.1f}% "
                f"identity over {self.length} columns — {self.n_identical} "
                f"identical, {self.n_mismatch} substitutions, {self.n_gap} gaps")

    def column_positions(self, column: int) -> tuple[int | None, int | None]:
        """Residue numbers in each sequence at an alignment column."""
        ia = len(self.aligned_a[:column].replace("-", ""))
        ib = len(self.aligned_b[:column].replace("-", ""))
        pa = (self.a.positions[ia] if self.aligned_a[column] != "-"
              and ia < len(self.a.positions) else None)
        pb = (self.b.positions[ib] if self.aligned_b[column] != "-"
              and ib < len(self.b.positions) else None)
        return pa, pb


def compare_sequences(a: NamedSequence, b: NamedSequence,
                      method: str = "global") -> SequenceComparison:
    """Align two sequences and enumerate every difference.

    ``method`` is ``"global"`` (Needleman–Wunsch through
    :func:`piezo1.core.sequence.align_global`) or ``"positional"``, which pairs
    residues by **residue number** with no gaps inserted. Positional is the
    honest choice when both sequences are already in the same numbering — an
    alignment can slide a run of residues to buy score, and then a "difference"
    is an alignment artefact rather than a substitution.
    """
    if method == "positional":
        shared = sorted(set(a.positions) & set(b.positions))
        aligned_a = "".join(a.at(p) or "-" for p in shared)
        aligned_b = "".join(b.at(p) or "-" for p in shared)
    elif method == "global":
        aligned_a, aligned_b = align_global(a.letters, b.letters)
    else:
        raise ValueError(f"unknown method {method!r}")

    identical = mismatch = gap = 0
    differences = []
    for column, (x, y) in enumerate(zip(aligned_a, aligned_b)):
        if x == "-" or y == "-":
            gap += 1
            continue
        if x == y:
            identical += 1
            continue
        mismatch += 1
        if method == "positional":
            shared_list = sorted(set(a.positions) & set(b.positions))
            pa = pb = shared_list[column]
        else:
            pa, pb = None, None
        differences.append((column, pa, pb, x, y))

    scored = identical + mismatch
    comparison = SequenceComparison(
        a=a, b=b, aligned_a=aligned_a, aligned_b=aligned_b,
        identity=(identical / scored if scored else 0.0),
        n_identical=identical, n_mismatch=mismatch, n_gap=gap,
        differences=differences)
    if method == "global":
        comparison.differences = [
            (c, *comparison.column_positions(c), x, y)
            for c, _pa, _pb, x, y in differences]
    return comparison


def load_named_sequences(structure=None) -> list[NamedSequence]:
    """Everything available to compare: references, translations, structure.

    Missing downloads are skipped rather than raising, so the viewer opens with
    whatever is present and says what is not.
    """
    from .sequence import human_sequence, mouse_sequence

    out: list[NamedSequence] = []
    for key, label, getter, numbering, accession in (
            ("uniprot_human", "UniProt human (Q92508)", human_sequence,
             "human", "Q92508"),
            ("uniprot_mouse", "UniProt mouse (E2JF22)", mouse_sequence,
             "mouse", "E2JF22")):
        try:
            letters = getter()
        except Exception:
            continue
        out.append(NamedSequence(key=key, label=label, letters=letters,
                                 numbering=numbering,
                                 source=f"UniProt {accession}"))

    for species, transcript in (("human", "ENST00000301015"),
                                ("mouse", "ENSMUST00000156333")):
        path = SEQUENCE_DIR / f"{transcript}_{species}_PIEZO1_cds.fasta"
        if not path.exists():
            continue
        dna = read_fasta(path)
        protein = translate(dna).rstrip("*")
        out.append(NamedSequence(
            key=f"cds_{species}", label=f"Translated CDS {species} ({transcript})",
            letters=protein, numbering=species, dna=dna,
            source=f"Ensembl {transcript}",
            note="translation of the reference-genome transcript"))

    if structure is not None:
        for chain in structure.chains:
            # Residue-level, not atom-level: one letter and one number per
            # residue, taken from the same mask so they cannot fall out of
            # step with each other.
            keep = (~structure.residue_hetero) & (structure.residue_chain == chain)
            if keep.sum() < 100:
                continue
            letters = structure.one_letter_sequence(chain)
            positions = [int(p) for p in structure.residue_seq[keep]]
            if len(letters) != len(positions):
                continue
            out.append(NamedSequence(
                key=f"structure_{chain}",
                label=f"{structure.name} chain {chain} (resolved)",
                letters=letters, positions=positions,
                source=f"PDB {structure.name}",
                note="resolved residues only — has gaps"))
    return out
