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


def _family_sequences():
    """``(key, label, letters, numbering, accession)`` for every family member.

    Reads the committed UniProt resources directly rather than importing
    ``analysis.homology``: this is ``core``, and the dependency arrow points
    the other way. The cost is that the ordering is stated here — human and
    mouse first, because they are what a PIEZO1 viewer is usually comparing —
    and a test checks the two lists hold the same members.
    """
    import json

    from ..config import RESOURCE_DIR

    order = ["human", "mouse", "rat", "human_piezo2", "mouse_piezo2",
             "zebrafish_piezo3",
             "worm_piezo", "fly_piezo", "plant_piezo", "dicty_piezo"]
    for key in order:
        path = RESOURCE_DIR / f"uniprot_{key}.json"
        if not path.exists():
            continue
        try:
            data = json.loads(path.read_text())
        except (ValueError, OSError):
            continue
        from .numbering_check import PROTEIN_NAMES

        # UniProt gives Arabidopsis PIEZO no gene name at all, so the label
        # falls back to what this project calls the protein rather than to the
        # resource key, which would put "plant_piezo" in a combo box.
        gene = data.get("gene") or PROTEIN_NAMES.get(key, key)
        organism = data.get("organism") or "?"
        yield (f"uniprot_{key}",
               f"{gene} — {organism} ({data['accession']}, "
               f"{data['length']} aa)",
               data["sequence"], key, data["accession"])


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

    And it is nonsense when they are not, which is why positional **raises**
    across different numbering systems. Since the viewer began offering all
    nine reviewed PIEZOs this stopped being hypothetical: no two of the nine
    share a length, so pairing human residue 2447 with plant residue 2447
    would compare two unrelated positions and report about 2,000 confident
    substitutions. A global alignment is the only defensible route there, and
    below Rost's line even that needs the reliability gate in
    ``analysis.alignment_windows``.
    """
    if method == "positional" and a.numbering != b.numbering:
        raise ValueError(
            f"cannot compare {a.label!r} and {b.label!r} position by position: "
            f"they are numbered in different systems ({a.numbering!r} and "
            f"{b.numbering!r}), so equal residue numbers are not corresponding "
            f"residues. Use method='global'.")
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
    """Everything available to compare: the family, translations, structure.

    Missing downloads are skipped rather than raising, so the viewer opens with
    whatever is present and says what is not.

    **All ten PIEZO references are offered, not just human and mouse.** Until
    Round 89 this returned two, so the comparison tool could align human
    against mouse and nothing else — while the project held annotation for six
    proteins and the family has ten. Each carries its own ``numbering``, which
    is the whole reason they can be listed together safely: **no two of the ten
    share a length**, so a residue number means nothing without the sequence it
    belongs to, and ``compare_sequences`` refuses to pair them positionally
    across different numbering.

    Nine are reviewed; the tenth, zebrafish piezo3, is TrEMBL — the third
    vertebrate paralogue, which human lost before the primate radiation.
    """
    out: list[NamedSequence] = []
    for key, label, letters, numbering, accession in _family_sequences():
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
        # Which numbering the file is actually in, measured from its own
        # residue names. It used to take the default "human", which was
        # harmless while the only references offered were human and mouse and
        # became a live hazard the moment nine were: a worm entry labelled
        # "human" could be compared position-by-position against human PIEZO1
        # and would report two thousand confident substitutions.
        try:
            from .numbering_check import identify_numbering

            chain_numbering = identify_numbering(structure).reference
        except Exception:
            chain_numbering = "unknown"
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
                numbering=chain_numbering,
                source=f"PDB {structure.name}",
                note=f"resolved residues only — has gaps; numbered as "
                     f"{chain_numbering}"))
    return out
