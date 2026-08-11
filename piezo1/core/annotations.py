"""Loading and querying the curated annotation resources.

Everything in ``piezo1/resources`` is reached through this module, so the rest
of the application never parses JSON or worries about which numbering system a
file uses. Each record keeps its provenance and confidence, and the API is
built so that "we do not know" is representable — a domain with a null range or
a variant that no structure resolves comes back as such rather than silently
becoming a wrong answer.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from functools import lru_cache

import numpy as np

from ..config import RESOURCE_DIR

__all__ = ["Domain", "ResidueGroup", "Variant", "Annotations", "load_annotations"]


@dataclass(frozen=True)
class Domain:
    id: str
    name: str
    category: str
    start: int | None
    end: int | None
    color: str
    description: str
    source: str
    confidence: str
    mouse_start: int | None = None
    mouse_end: int | None = None
    #: True for an element that lies *inside* another rather than partitioning
    #: the chain beside it. Guo & MacKinnon's cuff — the elbow and base within
    #: the anchor, the hairpin and PE helix within the CTD — are named features
    #: of a region this project already had, not a re-partition of it.
    #: :meth:`Annotations.domain_at` skips them, so "which domain is this
    #: residue in" keeps its old answer; ask :meth:`sub_elements_at` for the
    #: finer one. Adding them without this flag silently moved the anchor from
    #: first to thirteenth in the allosteric-betweenness ranking, because
    #: `domain_at` returns the *smallest* containing domain.
    sub_element: bool = False
    extra: dict = field(default_factory=dict)

    def contains(self, residue: int) -> bool:
        return (self.start is not None and self.end is not None
                and self.start <= residue <= self.end)

    @property
    def length(self) -> int:
        if self.start is None or self.end is None:
            return 0
        return self.end - self.start + 1


@dataclass(frozen=True)
class ResidueGroup:
    """A named set of individual functional residues."""

    id: str
    label: str
    category: str
    description: str
    source: str
    color: str
    evidence: str
    residues: tuple[int, ...]
    detail: tuple[dict, ...] = ()

    def non_conserved(self) -> list[int]:
        return [d["human"] for d in self.detail if d.get("conserved") is False]


@dataclass(frozen=True)
class Variant:
    id: str
    residue: int | None
    wt_aa: str | None
    mut_aa: str | None
    hgvs: str
    classification: str | None
    phenotype: str | None
    functional_effect: str | None
    effect_magnitude: str | None
    domain: str | None
    pmid: str | None
    source_url: str | None
    confidence: str
    modelled_in: tuple[str, ...]
    mouse_residue: int | None
    conserved: bool | None
    sequence_verified: bool

    @property
    def label(self) -> str:
        if self.wt_aa and self.residue and self.mut_aa:
            return f"{self.wt_aa}{self.residue}{self.mut_aa}"
        return self.hgvs or self.id

    def is_modelled_in(self, pdb: str) -> bool:
        return pdb.upper() in self.modelled_in


class Annotations:
    """Aggregated access to domains, functional residues and variants."""

    def __init__(self, species: str = "human") -> None:
        self.species = species
        self.domains: list[Domain] = []
        self.residue_groups: list[ResidueGroup] = []
        self.variants: list[Variant] = []
        self.meta: dict = {}
        self._load()

    # ------------------------------------------------------------- loading

    def _read(self, name: str) -> dict | None:
        path = RESOURCE_DIR / name
        if not path.exists():
            return None
        return json.loads(path.read_text())

    def _load(self) -> None:
        data = self._read("domains.json")
        if data:
            self.meta["domains"] = data.get("numbering", {})
            for d in data["domains"]:
                span = d.get(self.species, {})
                mouse = d.get("mouse", {})
                self.domains.append(Domain(
                    id=d["id"], name=d["name"], category=d["category"],
                    start=span.get("start"), end=span.get("end"),
                    color=d["color"], description=d["description"],
                    source=d["source"], confidence=d["confidence"],
                    mouse_start=mouse.get("start"), mouse_end=mouse.get("end"),
                    sub_element=bool(d.get("sub_element", False)),
                    extra={k: v for k, v in d.items()
                           if k in ("rule", "thu_index", "distal", "transmembrane")},
                ))

        data = self._read("functional_residues.json")
        if data:
            for g in data["groups"]:
                self.residue_groups.append(ResidueGroup(
                    id=g["id"], label=g["label"], category=g["category"],
                    description=g["description"], source=g["source"],
                    color=g["color"], evidence=g.get("evidence", "experimental"),
                    residues=tuple(r["human"] if self.species == "human"
                                   else r["mouse"]
                                   for r in g["residues"]
                                   if r.get("human") is not None),
                    detail=tuple(g["residues"]),
                ))

        data = self._read("variants.json")
        if data:
            self.meta["variant_coverage"] = data.get("coverage_summary", {})
            for v in data["variants"]:
                self.variants.append(Variant(
                    id=v["id"], residue=v.get("residue"),
                    wt_aa=v.get("wt_aa"), mut_aa=v.get("mut_aa"),
                    hgvs=v.get("hgvs_protein", ""),
                    classification=v.get("classification"),
                    phenotype=v.get("phenotype"),
                    functional_effect=v.get("functional_effect"),
                    effect_magnitude=v.get("effect_magnitude"),
                    domain=v.get("domain"), pmid=str(v.get("pmid") or ""),
                    source_url=v.get("source_url"),
                    confidence=v.get("confidence", "unknown"),
                    modelled_in=tuple(v.get("modelled_in", [])),
                    mouse_residue=v.get("mouse_residue"),
                    conserved=v.get("conserved"),
                    sequence_verified=bool(v.get("sequence_verified", False)),
                ))

    # ------------------------------------------------------------- queries

    def domain_at(self, residue: int) -> Domain | None:
        """Most specific *partitioning* domain containing ``residue``.

        Sub-elements are skipped — see :attr:`Domain.sub_element`. The chain's
        architecture is a partition and callers rely on it being one; the cuff
        elements are named features inside two of those parts.
        """
        hits = [d for d in self.domains
                if d.contains(residue) and not d.sub_element]
        return min(hits, key=lambda d: d.length) if hits else None

    def sub_elements_at(self, residue: int) -> list[Domain]:
        """Named features inside the domain this residue belongs to."""
        return [d for d in self.domains
                if d.sub_element and d.contains(residue)]

    def domains_by_category(self, category: str) -> list[Domain]:
        return [d for d in self.domains if d.category == category]

    def variants_at(self, residue: int) -> list[Variant]:
        return [v for v in self.variants if v.residue == residue]

    def variants_by_class(self, classification: str) -> list[Variant]:
        return [v for v in self.variants if v.classification == classification]

    def variant_classes(self) -> list[str]:
        seen: list[str] = []
        for v in self.variants:
            if v.classification and v.classification not in seen:
                seen.append(v.classification)
        return sorted(seen)

    def group(self, group_id: str) -> ResidueGroup | None:
        for g in self.residue_groups:
            if g.id == group_id:
                return g
        return None

    def annotate_residue(self, residue: int) -> dict:
        """Everything known about one residue, for a tooltip or info panel."""
        domain = self.domain_at(residue)
        return {
            "residue": residue,
            "domain": domain.name if domain else None,
            "domain_id": domain.id if domain else None,
            "domain_confidence": domain.confidence if domain else None,
            "groups": [g.label for g in self.residue_groups if residue in g.residues],
            "variants": [
                {"label": v.label, "classification": v.classification,
                 "phenotype": v.phenotype, "effect": v.functional_effect,
                 "pmid": v.pmid, "modelled_in": list(v.modelled_in)}
                for v in self.variants_at(residue)
            ],
        }

    # --------------------------------------------------------- vectorised

    def residue_mask(self, res_seq: np.ndarray, residues) -> np.ndarray:
        return np.isin(res_seq, np.asarray(list(residues), dtype=np.int32))

    def variant_residues(self, classification: str | None = None) -> list[int]:
        return sorted({v.residue for v in self.variants
                       if v.residue is not None
                       and (classification is None
                            or v.classification == classification)})


@lru_cache(maxsize=4)
def load_annotations(species: str = "human") -> Annotations:
    """Load (and cache) the annotation set for one species' numbering."""
    return Annotations(species=species)
