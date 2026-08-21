"""The PIEZO family census, as this project can read it.

Loads ``resources/family_findings.json`` and ``resources/family_constraint.json``
— thirteen statements and three per-residue tracks imported from the
``piezo_genes`` census project (see ``scripts/build_family_findings.py`` for the
gate they came through).

**Why this is its own module rather than part of ``annotations``.** Everything
in :mod:`piezo1.core.annotations` is *our* curated annotation: authored here,
corrected here, and answerable here. Everything here is somebody else's result.
That difference has to survive contact with the code, because the two fail in
different ways — a wrong domain boundary is a bug we can fix, and a superseded
census number is a quotation that has to be re-imported. So every record
carries its provenance, :class:`Finding` states what this project does with it
and what it does not establish, and nothing here is silently merged into an
annotation record.

**Numbering.** The constraint tracks are indexed by residue number *in their own
reference sequence* — human Q92508, human Q9H5I5, zebrafish A0AC58GFC9. Asking
for a value at a residue number is only meaningful once you know which of those
you are in, so :meth:`ConstraintTrack.value` takes a residue number and the
track knows its own numbering; converting between numbering systems goes through
:mod:`piezo1.core.sequence` as everywhere else, never by arithmetic.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from ..config import RESOURCE_DIR

__all__ = ["Finding", "ConstraintTrack", "EquivalentPosition", "PathogenicPore",
           "FamilyFindings", "load_family_findings", "load_constraint",
           "CONSTRAINT_GENES", "TRACK_NAMES", "constraint_numbering"]

#: The genes a per-residue track exists for, and the numbering each is in.
#: A gene absent from this map has no track — which is most of the family, and
#: is the honest state: the census scored the three vertebrate paralogues.
CONSTRAINT_GENES = {
    "PIEZO1": "human",
    "PIEZO2": "human_piezo2",
    "piezo3": "zebrafish_piezo3",
}

#: The three depths the constraint was computed at, weakest evidence last.
#: ``deep`` is one locus per genome across the paralogue's own orthologues
#: (121–192 sequences); ``vert`` is the vertebrate subset; ``family`` is the
#: whole-family alignment, where a column is comparing PIEZO1 with a protist.
TRACK_NAMES = ("deep_jsd", "vert_jsd", "family_jsd")


@dataclass(frozen=True)
class Finding:
    """One imported statement, with what it rests on and what it does not show."""

    key: str
    session: str
    kind: str                 # sequence | clinical | structure | census
    title: str
    statement: str
    numbers: dict
    source: str
    here: str
    caveat: str
    n_checks: int = 0

    @property
    def recorded_only(self) -> bool:
        """True when this project can hold the finding but not re-run it."""
        return self.here.startswith("recorded only")

    @property
    def module(self) -> str | None:
        """The module that explores it here, or None if it is recorded only."""
        return None if self.recorded_only else self.here.rsplit(".", 1)[0]


@dataclass(frozen=True)
class EquivalentPosition:
    """One alignment column carrying a pathogenic position in both disease genes."""

    msa_col: int
    element: str
    piezo1: int
    piezo1_aa: str
    piezo1_disease: str
    piezo2: int
    piezo2_aa: str
    piezo2_disease: str
    piezo3: int
    piezo3_aa: str
    note: str = ""

    @property
    def label(self) -> str:
        return (f"PIEZO1 {self.piezo1_aa}{self.piezo1} = "
                f"PIEZO2 {self.piezo2_aa}{self.piezo2}")


@dataclass(frozen=True)
class PathogenicPore:
    """A pathogenic position inside the pore module, and piezo3's residue at it."""

    gene: str
    resi: int
    aa: str
    element: str
    msa_col: int
    piezo3_resi: int
    piezo3_aa: str

    @property
    def kept_by_piezo3(self) -> bool:
        return self.aa == self.piezo3_aa

    @property
    def label(self) -> str:
        return f"{self.gene} {self.aa}{self.resi}"


@dataclass(frozen=True)
class ConstraintTrack:
    """Per-residue evolutionary constraint for one gene, in its own numbering.

    ``values`` is 1-indexed by residue through :meth:`value`; the raw array is
    0-indexed and the same length as the protein. A position with no score is
    ``nan`` rather than zero — an unscored residue and a completely
    unconstrained one are different statements, and colouring them the same
    would paint "we did not measure this" as "nothing here matters".
    """

    gene: str
    accession: str
    numbering: str
    length: int
    n_orthologues: int
    sequence: str
    track: str
    values: np.ndarray                      # (length,), nan where unscored
    reliable: np.ndarray                    # (length,) bool
    in_pore_module: np.ndarray              # (length,) bool
    census_domain: tuple = field(default=())

    def value(self, resi: int) -> float | None:
        """The score at a residue number, or None if unscored or out of range."""
        if not 1 <= resi <= self.length:
            return None
        v = float(self.values[resi - 1])
        return None if np.isnan(v) else v

    def residue(self, resi: int) -> str | None:
        if not 1 <= resi <= self.length:
            return None
        return self.sequence[resi - 1] or None

    @property
    def n_scored(self) -> int:
        return int(np.count_nonzero(~np.isnan(self.values)))

    def mask_for(self, residues) -> np.ndarray:
        """Boolean mask over the protein for an iterable of residue numbers."""
        mask = np.zeros(self.length, dtype=bool)
        for r in residues:
            if 1 <= int(r) <= self.length:
                mask[int(r) - 1] = True
        return mask

    def mean_over(self, residues) -> float | None:
        """Mean score over residue numbers, ignoring unscored ones."""
        mask = self.mask_for(residues) & ~np.isnan(self.values)
        if not mask.any():
            return None
        return float(self.values[mask].mean())


@dataclass(frozen=True)
class FamilyFindings:
    """Everything imported from the census, with its provenance attached."""

    provenance: dict
    findings: tuple
    tables: dict
    pathogenic_pore: tuple
    equivalent: tuple
    census: dict

    def by_key(self, key: str) -> Finding | None:
        for f in self.findings:
            if f.key == key:
                return f
        return None

    def by_kind(self, kind: str) -> tuple:
        return tuple(f for f in self.findings if f.kind == kind)

    @property
    def keys(self) -> tuple:
        return tuple(f.key for f in self.findings)

    @property
    def source(self) -> str:
        return (f"{self.provenance['source_project']} "
                f"@ {self.provenance['source_commit']}")

    def table(self, key: str) -> list:
        """Rows of one imported table, or an empty list if it was not imported."""
        return self.tables.get(key, {}).get("rows", [])

    def pathogenic_in(self, gene: str) -> tuple:
        return tuple(p for p in self.pathogenic_pore if p.gene == gene)


def constraint_numbering(gene: str) -> str | None:
    """Which numbering a gene's track is in, or None if there is no track."""
    return CONSTRAINT_GENES.get(gene)


@lru_cache(maxsize=1)
def load_family_findings() -> FamilyFindings:
    """The imported census, cached. Raises if the resource is missing."""
    path = RESOURCE_DIR / "family_findings.json"
    data = json.loads(path.read_text())
    findings = tuple(
        Finding(key=f["key"], session=f["session"], kind=f["kind"],
                title=f["title"], statement=f["statement"], numbers=f["numbers"],
                source=f["source"], here=f["here"], caveat=f["caveat"],
                n_checks=f.get("n_checks", 0))
        for f in data["findings"])
    pathogenic = tuple(PathogenicPore(**p) for p in data["pore_module_pathogenic"])
    equivalent = tuple(EquivalentPosition(**e) for e in data["equivalent_positions"])
    return FamilyFindings(provenance=data["provenance"], findings=findings,
                          tables=data["tables"], pathogenic_pore=pathogenic,
                          equivalent=equivalent, census=data["census"])


@lru_cache(maxsize=8)
def load_constraint(gene: str = "PIEZO1",
                    track: str = "deep_jsd") -> ConstraintTrack:
    """One gene's per-residue constraint track.

    Raises :class:`KeyError` for a gene with no track rather than returning an
    empty one: "this protein was not scored" and "this protein scored zero
    everywhere" must not be reachable through the same value.
    """
    if gene not in CONSTRAINT_GENES:
        raise KeyError(
            f"no constraint track for {gene!r}; the census scored "
            f"{', '.join(sorted(CONSTRAINT_GENES))}")
    if track not in TRACK_NAMES:
        raise KeyError(f"unknown track {track!r}; have {', '.join(TRACK_NAMES)}")
    data = json.loads((RESOURCE_DIR / "family_constraint.json").read_text())
    g = data["genes"][gene]
    raw = g[track]
    values = np.array([np.nan if v is None else float(v) for v in raw], dtype=float)
    return ConstraintTrack(
        gene=gene, accession=g["accession"], numbering=g["numbering"],
        length=g["length"], n_orthologues=g["n_orthologues"],
        sequence=g["sequence"], track=track, values=values,
        reliable=np.array(g["reliable"], dtype=bool),
        in_pore_module=np.array(g["in_pore_module"], dtype=bool),
        census_domain=tuple(g["census_domain"]))
