"""The PIEZO family, and what a sequence can and cannot say about it.

Nine proteins. That is the whole reviewed family — ``reviewed:true AND
protein_name:piezo`` against UniProt returns exactly these — and the fact that
it is *enumerable* is the most consequential thing about it for this project.
A homology search is a tool for when the answer set is unknown. Here it is
known, pinned, and small enough to hold in one table:

===============  ==========  =======  ====  ==============================
key              accession   length   TM    what it adds
===============  ==========  =======  ====  ==============================
human            Q92508        2521    38   the reference
mouse            E2JF22        2547    38   most of the mechanism literature
rat              Q0KL00        2535    38   most of the electrophysiology
human_piezo2     Q9H5I5        2752    38   the paralogue, same species
mouse_piezo2     Q8CD54        2822    38   the paralogue with a structure
worm_piezo       A0A061ACU2    2442    36   ~800 Myr, and two structures
fly_piezo        M9MSG8        2551    40   ~800 Myr, the other branch
plant_piezo      F4IN58        2462    35   not an animal
dicty_piezo      Q54S52        3080    35   not an animal, not a plant
===============  ==========  =======  ====  ==============================

**Every number here is reported beside its own null, and that is the point of
the module rather than a decoration.**

A global BLOSUM62 alignment of two unrelated 2,500-residue proteins returns
about 22% identity for nothing but composition and the alignment's own freedom
to slide. **15 of the 36 pairs in this family fall below Rost's 30% line**, and
for those the percentage is mostly that.

The extreme case is PEZO-1 against Arabidopsis PIEZO. Identity **0.238**;
composition-matched shuffles of the same partner give **0.225 ± 0.009**, so
**z = 1.5** — a reader shown "24% identical" would call it weakly related, and
a scrambled sequence gives 22.5%. Run the *local* alignment on the identical
pair and the score is **391 against a null of 51, z = 64**. The homology is not
marginal at all. Two pairs are like this outright (the other is rat against
pzoA, z = 2.6) and every non-animal pair is within a few sigma.

**Percent identity is simply the wrong statistic below that line**, which is
what the line means, and why this module refuses to report one without the
other.

That measurement is also the answer to whether this application needs a BLAST
client, and the answer is no — see ``docs/HOMOLOGY_SEARCH.md``. What it needed
was BLAST's *statistic*, which is here.

**What the family does for the science.** ``analysis.paralogue`` asks whether
the gating mechanism is PIEZO1's or the fold's, and can only ask it of PIEZO2 —
one duplication, entirely within vertebrates. The invertebrates widen that to
~800 Myr. A plant and an amoeba widen it to the root of the eukaryotes, and
they arrive with a specific, checkable consequence: **the 38-helix architecture
is not universal.** Worm has 36, fly 40, plant and Dictyostelium 35. Anything
that transfers a helix index across the family by number is wrong, and
``homology_sites`` is how a position crosses instead.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

from ..config import RESOURCE_DIR
from ..core.numbering_check import (NON_ANIMAL_REFERENCES, PIEZO1_REFERENCES,
                                    PIEZO2_REFERENCES, PROTEIN_NAMES,
                                    REFERENCES, reference_entry)
from ..parameters import PARAMETERS as _P

__all__ = ["FamilyMember", "family", "member", "Relationship", "relationship",
           "FamilyMatrix", "family_matrix", "align_pair", "shuffled_null",
           "GROUPS", "group_of"]


#: How the family divides. Used for labelling and for the ordering of the
#: matrix; never for deciding an answer, which is what the alignment is for.
GROUPS = {
    "PIEZO1": PIEZO1_REFERENCES,
    "PIEZO2": PIEZO2_REFERENCES,
    "invertebrate": ("worm_piezo", "fly_piezo"),
    "non-animal": NON_ANIMAL_REFERENCES,
}


def group_of(key: str) -> str:
    for name, members in GROUPS.items():
        if key in members:
            return name
    return "unclassified"


# --------------------------------------------------------------------------
# The family
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class FamilyMember:
    """One reviewed PIEZO, read from its committed UniProt resource.

    ``n_transmembrane`` is carried because it is the field that refuses the
    obvious shortcut: five of the nine have 38 helices and four do not, so
    "TM12 of PIEZO1" has no counterpart in the plant protein that can be found
    by counting.
    """

    key: str
    accession: str
    protein: str
    organism: str
    gene: str
    length: int
    n_transmembrane: int
    group: str

    @property
    def label(self) -> str:
        return f"{self.protein} ({self.organism})"

    def summary(self) -> str:
        return (f"{self.key}: {self.protein}, {self.organism}, "
                f"{self.length} aa, {self.n_transmembrane} TM "
                f"[{self.accession}]")


_FAMILY: tuple[FamilyMember, ...] | None = None


def family() -> tuple[FamilyMember, ...]:
    """Every reviewed PIEZO, in :data:`GROUPS` order.

    Built from the committed resources rather than from a list written here,
    so the family cannot disagree with the sequences the rest of the project
    reads. Missing a resource raises rather than silently shortening the
    family: a comparison quietly run over eight of nine would report a smaller
    range of divergence and look like a tighter result.
    """
    global _FAMILY
    if _FAMILY is not None:
        return _FAMILY

    order = [k for group in GROUPS.values() for k in group]
    order += [k for k in REFERENCES if k not in order]
    members = []
    for key in order:
        path = RESOURCE_DIR / f"uniprot_{key}.json"
        if not path.exists():
            raise FileNotFoundError(
                f"{path.name} is missing — run "
                f"scripts/build_uniprot_annotations.py")
        data = json.loads(path.read_text())
        members.append(FamilyMember(
            key=key, accession=data["accession"],
            protein=PROTEIN_NAMES.get(key, "unknown"),
            organism=data.get("organism", "?"), gene=data.get("gene") or "?",
            length=int(data["length"]),
            n_transmembrane=int(data["n_transmembrane"]),
            group=group_of(key)))
    _FAMILY = tuple(members)
    return _FAMILY


def member(key: str) -> FamilyMember:
    for m in family():
        if m.key == key:
            return m
    raise KeyError(f"{key!r} is not a family member; have "
                   f"{[m.key for m in family()]}")


def _sequence(key: str) -> str:
    return reference_entry(key)["sequence"]


# --------------------------------------------------------------------------
# Alignment, and the null it is read against
# --------------------------------------------------------------------------

def _aligner(mode: str):
    from Bio import Align
    from Bio.Align import substitution_matrices

    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = _P.value("homology.gap_open")
    aligner.extend_gap_score = _P.value("homology.gap_extend")
    return aligner


def align_pair(a: str, b: str) -> tuple[float, int, float]:
    """``(identity, aligned_columns, local_score)`` for two sequences.

    The two statistics come from two different alignments on purpose. Identity
    is read off the **global** alignment, because that is the number a reader
    means by "percent identity" and the number a family matrix is expected to
    contain. The score is the **local** one, because that is the statistic that
    still separates these proteins from noise once identity has stopped.
    """
    glob = _aligner("global").align(a, b)[0]
    top, bottom = str(glob[0]), str(glob[1])
    pairs = [(x, y) for x, y in zip(top, bottom) if x != "-" and y != "-"]
    matched = sum(x == y for x, y in pairs)
    local = float(_aligner("local").score(a, b))
    return (matched / max(len(pairs), 1), len(pairs), local)


def shuffled_null(a: str, b: str, replicates: int | None = None,
                  seed: int = 0) -> tuple[np.ndarray, np.ndarray]:
    """Both statistics for ``a`` against shuffles of ``b``.

    The shuffle preserves ``b``'s **amino-acid composition** exactly, which is
    what makes it the right null: a composition-matched sequence is what the
    aligner would score if there were no homology at all, and the difference
    between the two is therefore the part of the number that is about these
    proteins. A null drawn from a uniform amino-acid distribution would be
    easier to beat and would flatter every pair in the matrix.

    ``b`` is shuffled rather than ``a`` so that a whole row of the matrix is
    scored against nulls of the same query, which is how the row is read.
    """
    replicates = int(_P.value("homology.null_replicates")) \
        if replicates is None else replicates
    rng = np.random.default_rng(seed)
    letters = np.frombuffer(b.encode(), dtype="S1")
    identities, scores = [], []
    for _ in range(replicates):
        shuffled = rng.permutation(letters).tobytes().decode()
        identity, _columns, local = align_pair(a, shuffled)
        identities.append(identity)
        scores.append(local)
    return np.asarray(identities), np.asarray(scores)


def _z(value: float, null: np.ndarray) -> float:
    sd = float(null.std(ddof=1))
    return float("inf") if sd == 0 else (value - float(null.mean())) / sd


@dataclass
class Relationship:
    """How two family members are related, each statistic beside its null.

    Both are carried because on the distant pairs they disagree about the
    strength of the evidence by more than an order of magnitude, and the
    disagreement is the finding rather than a nuisance.
    """

    a: str
    b: str
    identity: float
    aligned_columns: int
    local_score: float
    null_identity: float
    null_identity_sd: float
    null_local: float
    null_local_sd: float
    replicates: int
    meta: dict = field(default_factory=dict)

    @property
    def identity_z(self) -> float:
        sd = self.null_identity_sd
        return float("inf") if sd == 0 else (self.identity - self.null_identity) / sd

    @property
    def local_z(self) -> float:
        sd = self.null_local_sd
        return float("inf") if sd == 0 else (self.local_score - self.null_local) / sd

    @property
    def in_twilight_zone(self) -> bool:
        """Below Rost's line, where a percentage stops being evidence."""
        return self.identity < _P.value("homology.twilight_identity")

    @property
    def identity_beats_null(self) -> bool:
        return self.identity_z >= _P.value("homology.min_z")

    @property
    def local_beats_null(self) -> bool:
        return self.local_z >= _P.value("homology.min_z")

    @property
    def verdict(self) -> str:
        """What this pair supports, in words that do not overstate it."""
        if not self.local_beats_null:
            return "no detectable homology"
        if not self.in_twilight_zone:
            return "homologous; identity is itself informative"
        if self.identity_beats_null:
            return ("homologous by local score; identity is in the twilight "
                    "zone and barely above its own null")
        return ("homologous by local score; identity is indistinguishable "
                "from chance and must not be quoted alone")

    def summary(self) -> str:
        return (f"{self.a} vs {self.b}: identity {self.identity:.3f} "
                f"(null {self.null_identity:.3f}, z {self.identity_z:.1f}) "
                f"over {self.aligned_columns} columns; local score "
                f"{self.local_score:.0f} (null {self.null_local:.0f}, "
                f"z {self.local_z:.1f}) — {self.verdict}")


def relationship(a: str, b: str, replicates: int | None = None,
                 seed: int = 0) -> Relationship:
    """Align two family members and place both statistics against their null."""
    seq_a, seq_b = _sequence(a), _sequence(b)
    identity, columns, local = align_pair(seq_a, seq_b)
    null_identity, null_local = shuffled_null(seq_a, seq_b, replicates, seed)
    return Relationship(
        a=a, b=b, identity=identity, aligned_columns=columns,
        local_score=local,
        null_identity=float(null_identity.mean()),
        null_identity_sd=float(null_identity.std(ddof=1)),
        null_local=float(null_local.mean()),
        null_local_sd=float(null_local.std(ddof=1)),
        replicates=len(null_identity),
        meta={"gap_open": _P.value("homology.gap_open"),
              "gap_extend": _P.value("homology.gap_extend"),
              "matrix": "BLOSUM62",
              "null": "composition-matched shuffle of the second sequence"})


# --------------------------------------------------------------------------
# The whole family at once
# --------------------------------------------------------------------------

@dataclass
class FamilyMatrix:
    """Every pair, with the identity and the local score kept apart."""

    keys: tuple[str, ...]
    identity: np.ndarray            # symmetric, 1.0 on the diagonal
    identity_z: np.ndarray
    local_z: np.ndarray
    relationships: dict = field(default_factory=dict)

    def get(self, a: str, b: str) -> Relationship | None:
        return self.relationships.get((a, b)) or self.relationships.get((b, a))

    @property
    def n_in_twilight(self) -> int:
        return sum(1 for r in self.relationships.values() if r.in_twilight_zone)

    @property
    def n_identity_indistinguishable(self) -> int:
        """Pairs whose identity is within the null, and whose score is not.

        The count this module exists to be able to state. Every one of these
        is a pair a reader would call weakly related from the percentage and
        strongly related from the alignment.
        """
        return sum(1 for r in self.relationships.values()
                   if r.local_beats_null and not r.identity_beats_null)

    def distances(self) -> np.ndarray:
        """``1 - identity``, the matrix a tree is built from."""
        return 1.0 - self.identity

    def row(self, key: str) -> list[Relationship]:
        return [r for r in (self.get(key, other) for other in self.keys)
                if r is not None]

    def summary(self) -> str:
        return (f"{len(self.keys)} PIEZOs, {len(self.relationships)} pairs; "
                f"{self.n_in_twilight} below the "
                f"{_P.value('homology.twilight_identity'):.0%} twilight line, "
                f"of which {self.n_identity_indistinguishable} have an "
                f"identity indistinguishable from chance while their local "
                f"alignment score is not")


_MATRIX_CACHE: dict = {}


def family_matrix(keys=None, replicates: int | None = None,
                  seed: int = 0) -> FamilyMatrix:
    """The whole comparison: 36 pairs, each with its null. **About two minutes.**

    The cost is the null, not the answer — 36 alignments take three seconds and
    the 720 shuffled ones take the rest. That is the price of not reporting a
    bare percentage, and it is why the GUI runs this on a worker thread and the
    CLI says what it is doing.

    Memoised per process rather than cached to disk. A stored identity matrix
    would be a second place for the family to be wrong; this one is derived
    entirely from committed sequences and a registered scoring scheme, so
    recomputing it is the cheaper kind of correctness.
    """
    members = list(keys) if keys is not None else [m.key for m in family()]
    signature = (tuple(members), replicates, seed,
                 _P.value("homology.gap_open"),
                 _P.value("homology.gap_extend"),
                 _P.value("homology.null_replicates"))
    if signature in _MATRIX_CACHE:
        return _MATRIX_CACHE[signature]
    n = len(members)
    identity = np.eye(n)
    identity_z = np.full((n, n), np.inf)
    local_z = np.full((n, n), np.inf)
    relationships: dict = {}
    for i in range(n):
        for j in range(i + 1, n):
            rel = relationship(members[i], members[j], replicates, seed)
            relationships[(members[i], members[j])] = rel
            identity[i, j] = identity[j, i] = rel.identity
            identity_z[i, j] = identity_z[j, i] = rel.identity_z
            local_z[i, j] = local_z[j, i] = rel.local_z
    matrix = FamilyMatrix(keys=tuple(members), identity=identity,
                          identity_z=identity_z, local_z=local_z,
                          relationships=relationships)
    _MATRIX_CACHE[signature] = matrix
    return matrix
