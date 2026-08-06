"""Assembling an analysis set from sources of different evidential weight.

Round 20 established that the binding constraint on this project's central
claim is data. Round 27 went looking for more, and found that the extra
variants are **not the same kind of evidence** as the original ones.

Three levels, kept apart rather than pooled:

* ``measured`` — the direction was determined by electrophysiology. This is what
  the original 68 curated variants carry where they carry anything.
* ``disease_mechanism`` — the direction is *inferred* from which disease the
  variant causes, because dehydrated hereditary stomatocytosis is dominant
  gain-of-function and generalised lymphatic dysplasia is recessive loss of
  function. Weaker: it assumes the variant acts by the usual mechanism for its
  disease.
* ``ambiguous`` — reported under both diseases, so no direction is claimed.

**Why the split matters more than the count.** Pooling them would let 20
mechanism-inferred variants outvote 26 measured ones while looking like a
larger, better study. Any analysis using this module chooses its level
explicitly and reports what it chose.

The original curated resource is **not modified**. Round 7's and Round 22's
recorded results reference it, and silently growing the set underneath a frozen
result would invalidate it without anything appearing to change.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

from ..config import RESOURCE_DIR
from ..core.annotations import load_annotations

__all__ = ["VariantEntry", "VariantSet", "load_clinvar_variants",
           "build_analysis_set", "EVIDENCE_LEVELS"]

#: Ordered weakest-last, so "at least this strong" is a prefix of the list.
EVIDENCE_LEVELS = ("measured", "disease_mechanism")


@dataclass
class VariantEntry:
    """One variant with a direction and the strength of evidence for it."""

    label: str
    residue: int
    wt_aa: str
    mut_aa: str
    classification: str          # GoF | LoF
    evidence: str                # measured | disease_mechanism
    kind: str                    # missense | nonsense | frameshift
    source: str
    citation: str = ""
    conflict: str = ""

    @property
    def is_missense(self) -> bool:
        return self.kind == "missense"


@dataclass
class VariantSet:
    """A directional set assembled at a stated evidence level."""

    entries: list[VariantEntry] = field(default_factory=list)
    levels: tuple[str, ...] = ()
    excluded: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.entries)

    def missense(self) -> "VariantSet":
        return VariantSet([e for e in self.entries if e.is_missense],
                          self.levels, dict(self.excluded))

    def counts(self) -> dict:
        out: dict[str, int] = {}
        for entry in self.entries:
            out[entry.classification] = out.get(entry.classification, 0) + 1
        return out

    def by_evidence(self) -> dict:
        out: dict[str, int] = {}
        for entry in self.entries:
            out[entry.evidence] = out.get(entry.evidence, 0) + 1
        return out

    def summary(self) -> str:
        counts = self.counts()
        return (f"{len(self.entries)} variants "
                f"({counts.get('GoF', 0)} GoF, {counts.get('LoF', 0)} LoF) "
                f"at evidence level(s) {', '.join(self.levels)}; "
                f"by evidence {self.by_evidence()}"
                + (f"; excluded {self.excluded}" if self.excluded else ""))


def load_clinvar_variants() -> list[dict]:
    """The ClinVar-derived resource, or an empty list if not built."""
    path = RESOURCE_DIR / "variants_clinvar.json"
    if not path.exists():
        return []
    return json.loads(path.read_text())["variants"]


def _from_curated() -> list[VariantEntry]:
    entries = []
    for variant in load_annotations("human").variants:
        if variant.classification not in ("GoF", "LoF"):
            continue
        if not (variant.residue and variant.wt_aa and variant.mut_aa):
            continue
        single = (len(variant.wt_aa) == 1 and len(variant.mut_aa) == 1
                  and variant.mut_aa.isalpha())
        kind = "missense" if single else "other"
        entries.append(VariantEntry(
            label=variant.label, residue=int(variant.residue),
            wt_aa=variant.wt_aa, mut_aa=variant.mut_aa,
            classification=variant.classification, evidence="measured",
            kind=kind, source="curated", citation=str(variant.pmid or "")))
    return entries


def build_analysis_set(levels=("measured",),
                       include_ambiguous: bool = False) -> VariantSet:
    """Assemble a directional variant set at the requested evidence levels.

    ``levels`` defaults to ``("measured",)`` — the original curated set — so a
    caller who does not think about evidence strength gets the conservative
    answer rather than the largest one.
    """
    levels = tuple(levels)
    unknown = set(levels) - set(EVIDENCE_LEVELS)
    if unknown:
        raise ValueError(f"unknown evidence level(s): {sorted(unknown)}")

    entries: list[VariantEntry] = []
    excluded: dict[str, int] = {}

    if "measured" in levels:
        entries.extend(_from_curated())
    seen = {e.label for e in entries}

    if "disease_mechanism" in levels:
        for record in load_clinvar_variants():
            if record["ambiguous_direction"] and not include_ambiguous:
                excluded["ambiguous direction"] = (
                    excluded.get("ambiguous direction", 0) + 1)
                continue
            if not record["classification"]:
                excluded["no direction implied"] = (
                    excluded.get("no direction implied", 0) + 1)
                continue
            if record["label"] in seen:
                excluded["already curated"] = (
                    excluded.get("already curated", 0) + 1)
                continue
            seen.add(record["label"])
            entries.append(VariantEntry(
                label=record["label"], residue=int(record["residue"]),
                wt_aa=record["wt_aa"], mut_aa=record["mut_aa"],
                classification=record["classification"],
                evidence="disease_mechanism", kind=record["kind"],
                source="clinvar", citation=record.get("mechanism_citation", "")))

    return VariantSet(entries=entries, levels=levels, excluded=excluded)


def disagreements() -> list[dict]:
    """Variants where the curated direction and the inferred one differ.

    Reported rather than resolved. A disagreement here means the literature is
    mixed, and picking a side by fiat would hide that from any test built on
    the set.
    """
    curated = {e.label: e for e in _from_curated()}
    out = []
    for record in load_clinvar_variants():
        entry = curated.get(record["label"])
        if entry is None or not record["classification"]:
            continue
        if record["ambiguous_direction"]:
            continue
        if entry.classification != record["classification"]:
            out.append({"label": record["label"],
                        "curated": entry.classification,
                        "inferred": record["classification"],
                        "conditions": record.get("conditions", [])})
    return out
