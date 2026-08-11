"""The catalogue of available PIEZO structures."""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

from ..config import RESOURCE_DIR, STRUCTURE_DIR

__all__ = ["StructureRecord", "Registry", "load_registry"]


@dataclass(frozen=True)
class StructureRecord:
    pdb: str
    file: str
    species: str
    state: str
    gating: str
    resolution: float | None
    released: str
    title: str
    journal: str | None
    year: int | None
    pmid: str | None
    doi: str | None
    emdb: tuple[str, ...]
    n_atoms: int
    n_protomers: int
    protomer_chains: tuple[dict, ...]
    ligands: tuple[str, ...]
    note: str
    recommended_for: tuple[str, ...] = ()
    #: "PIEZO1" or "PIEZO2", **measured** at build time by scoring the file's
    #: own residue names against every reference sequence rather than curated.
    #: Defaults to "unknown" so a stale resource is visibly stale instead of
    #: silently claiming PIEZO1 — which is what the species field did until
    #: Round 83, when PIEZO2 filed as "mouse" slipped past the overlay guard.
    protein: str = "unknown"

    @property
    def path(self) -> Path:
        return STRUCTURE_DIR / self.file

    @property
    def available(self) -> bool:
        return self.path.exists()

    @property
    def residue_range(self) -> tuple[int, int] | None:
        if not self.protomer_chains:
            return None
        c = self.protomer_chains[0]
        return int(c["first"]), int(c["last"])

    @property
    def is_piezo2(self) -> bool:
        return self.protein.upper() == "PIEZO2"

    @property
    def numbering_species(self) -> str:
        """Which numbering the deposited residue numbers use."""
        return "human" if self.species == "human" else "mouse"

    def label(self) -> str:
        res = f"{self.resolution:.2f} A" if self.resolution else "n/a"
        return f"{self.pdb}  ({self.species}, {self.state}, {res})"

    def citation(self) -> str:
        bits = [b for b in (self.journal, str(self.year) if self.year else None)
                if b]
        cite = " ".join(bits)
        if self.pmid:
            cite += f"  PMID {self.pmid}"
        return cite


@dataclass
class Registry:
    entries: list[StructureRecord] = field(default_factory=list)

    def __len__(self) -> int:
        return len(self.entries)

    def __iter__(self):
        return iter(self.entries)

    def get(self, pdb: str) -> StructureRecord | None:
        pdb = pdb.upper()
        for e in self.entries:
            if e.pdb == pdb:
                return e
        return None

    def available(self) -> list[StructureRecord]:
        return [e for e in self.entries if e.available]

    def by_species(self, species: str) -> list[StructureRecord]:
        return [e for e in self.entries if e.species == species]

    def by_state(self, state: str) -> list[StructureRecord]:
        return [e for e in self.entries if e.state == state]

    def recommended(self, purpose: str) -> list[StructureRecord]:
        return [e for e in self.entries if purpose in e.recommended_for]

    def default(self) -> StructureRecord | None:
        for e in self.recommended("default"):
            if e.available:
                return e
        avail = self.available()
        return avail[0] if avail else None

    def morph_pairs(self) -> list[tuple[StructureRecord, StructureRecord]]:
        """Curved/flat endpoint pairs that share a species and a paper."""
        starts = self.recommended("morph_start")
        ends = self.recommended("morph_end")
        pairs = []
        for a in starts:
            for b in ends:
                if a.species != b.species or not (a.available and b.available):
                    continue
                # Same study, or at least the same construct family.
                if a.pmid and a.pmid == b.pmid:
                    pairs.append((a, b))
        # Fall back to any same-species pairing if no shared-paper pair exists.
        if not pairs:
            pairs = [(a, b) for a in starts for b in ends
                     if a.species == b.species and a.available and b.available]
        return pairs


@lru_cache(maxsize=1)
def load_registry() -> Registry:
    path = RESOURCE_DIR / "structures.json"
    if not path.exists():
        return Registry()
    data = json.loads(path.read_text())
    entries = []
    for e in data["entries"]:
        entries.append(StructureRecord(
            pdb=e["pdb"], file=e["file"], species=e.get("species", "unknown"),
            state=e.get("state", "unclassified"), gating=e.get("gating", "unknown"),
            resolution=e.get("resolution"), released=e.get("released", ""),
            title=e.get("title", ""), journal=e.get("journal"),
            year=e.get("year"), pmid=str(e["pmid"]) if e.get("pmid") else None,
            doi=e.get("doi"), emdb=tuple(e.get("emdb", [])),
            n_atoms=e.get("n_atoms", 0), n_protomers=e.get("n_protomers", 0),
            protomer_chains=tuple(e.get("protomer_chains", [])),
            ligands=tuple(e.get("ligands", [])), note=e.get("note", ""),
            recommended_for=tuple(e.get("recommended_for", [])),
            protein=e.get("protein", "unknown"),
        ))
    return Registry(entries=entries)
