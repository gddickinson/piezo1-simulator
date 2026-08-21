"""Where human disease sits on the channel, tested against this project's own data.

The census's sharpest clinical result is that pathogenic missense variation
piles into the pore module — 17–18% of the protein holding 33% of the pathogenic
positions and 11% of the benign ones, odds ratio 3.9, one-sided Fisher
P = 0.0014 — and that this is the same region half a billion years of evolution
refused to change.

**This module re-runs it here, and the re-run is not a formality.** Three things
differ, and each is a way the original could have been an artefact of its own
inputs:

1. **The boundaries are ours.** The pore module is defined from ``domains.json``
   — outer helix, cap, inner helix, CTD — which is built from UniProt topology
   and Guo & MacKinnon's named elements, not from the census's mouse-literature
   bands. The two disagree about where the anchor ends.
2. **The pathogenic set is ours.** 68 hand-curated variants plus 232 ClinVar
   pathogenic/likely-pathogenic entries, each already gated against Q92508.
3. **The benign comparator is different in kind.** The census used ClinVar's
   benign labels. This uses gnomAD population missense — variation that exists
   in people rather than variation somebody classified as harmless. That is a
   *better* control for the ascertainment problem the census names in its own
   caveat, because a population frequency is not a clinical judgement and does
   not inherit the clinician's idea of where to look.

It is also PIEZO1 only, where the census pooled PIEZO1 and PIEZO2. So this is a
partial, independent re-test rather than a reproduction, and
:attr:`Enrichment.note` says so wherever the number is printed.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from math import comb

import numpy as np

from ..config import RESOURCE_DIR
from ..core.annotations import Annotations, load_annotations
from ..core.sequence import human_sequence
from ..core.family import ConstraintTrack, load_constraint, load_family_findings
from ..parameters import PARAMETERS as _P
from .validation import auroc

__all__ = ["Enrichment", "ClassifierResult", "PORE_MODULE_DOMAINS",
           "PARTITIONS", "pore_module_residues", "pathogenic_positions",
           "population_positions", "fisher_exact_greater",
           "pore_module_enrichment", "both_partitions", "boundary_disagreement",
           "constraint_classifier", "census_comparison"]

#: The two answers to "which residues are the pore module", and the reason both
#: are offered. Measured on human PIEZO1, the two projects' boundaries agree to
#: within four residues on the cap and within twelve on the inner helix and CTD
#: — and place the **anchor 141 residues apart and the outer helix 120 apart**.
#: A test whose answer depends on that choice has to be run both ways, the way
#: :mod:`piezo1.analysis.paralogue` reports its dome comparison naively and
#: coverage-matched: the gap between the two answers is the result.
PARTITIONS = ("ours", "census")

#: The pore module, in this project's own domain vocabulary. Matches the
#: census's ``OH + CED + IH + CTD`` in intent; the residue ranges are ours.
#: Sub-elements (elbow, base, hairpin, PE helix) are not listed because they lie
#: inside these and would double-count.
PORE_MODULE_DOMAINS = ("outer_helix", "cap", "inner_helix", "ctd")


@dataclass(frozen=True)
class Enrichment:
    """A 2x2 test of whether pathogenic positions concentrate in a region."""

    region: str
    n_region: int
    n_protein: int
    pathogenic_in: int
    pathogenic_out: int
    comparator_in: int
    comparator_out: int
    odds_ratio: float | None
    p_value: float
    comparator: str
    note: str = ""

    @property
    def region_fraction(self) -> float:
        return self.n_region / self.n_protein if self.n_protein else 0.0

    @property
    def pathogenic_fraction(self) -> float:
        total = self.pathogenic_in + self.pathogenic_out
        return self.pathogenic_in / total if total else 0.0

    @property
    def comparator_fraction(self) -> float:
        total = self.comparator_in + self.comparator_out
        return self.comparator_in / total if total else 0.0

    @property
    def significant(self) -> bool:
        return self.p_value <= _P.value("stats.alpha")

    def summary(self) -> str:
        odds = "undefined" if self.odds_ratio is None else f"{self.odds_ratio:.2f}"
        return (f"{self.region} is {self.region_fraction:.0%} of the protein and "
                f"holds {self.pathogenic_in} of "
                f"{self.pathogenic_in + self.pathogenic_out} pathogenic "
                f"positions ({self.pathogenic_fraction:.0%}) against "
                f"{self.comparator_fraction:.0%} of {self.comparator}: "
                f"odds ratio {odds}, one-sided Fisher P = {self.p_value:.4g}")


@dataclass(frozen=True)
class ClassifierResult:
    """The census constraint used as a classifier on *our* variant positions."""

    auc: float
    n_positive: int
    n_negative: int
    mean_positive: float
    mean_negative: float
    track: str
    census_auc: float | None = None
    note: str = ""

    @property
    def agrees_with_census(self) -> bool:
        if self.census_auc is None:
            return False
        return abs(self.auc - self.census_auc) <= 0.1


def fisher_exact_greater(a: int, b: int, c: int, d: int) -> float:
    """One-sided Fisher exact P for the table [[a, b], [c, d]].

    Written out rather than imported for the same reason every other statistic
    in this project is: the convention decides the result. The tail summed is
    "at least ``a`` of the first row's total in the first column", which is the
    direction the census tested — pathogenic variation *concentrating* in the
    module — and testing the other tail on the same data would be a different
    claim entirely.
    """
    n = a + b + c + d
    row1, col1 = a + b, a + c
    total = comb(n, col1)
    if total == 0:
        return 1.0
    upper = min(row1, col1)
    p = sum(comb(row1, k) * comb(n - row1, col1 - k)
            for k in range(a, upper + 1))
    return min(1.0, p / total)


def _odds_ratio(a: int, b: int, c: int, d: int) -> float | None:
    if b == 0 or c == 0:
        return None if (a == 0 or d == 0) else float("inf")
    return (a * d) / (b * c)


def pore_module_residues(annotations: Annotations | None = None,
                         partition: str = "ours") -> set[int]:
    """Every residue number in the pore module, on the chosen boundaries.

    ``ours`` is ``domains.json`` — outer helix, cap, inner helix, CTD.
    ``census`` is the imported ``domain_map`` row the census computed its own
    numbers on, which additionally sweeps in the 117 residues it calls
    ``pore_linker`` and starts the module 120 residues earlier.
    """
    if partition not in PARTITIONS:
        raise ValueError(f"partition must be one of {PARTITIONS}, not {partition!r}")
    if partition == "census":
        row = next((r for r in load_family_findings().table("domain_map")
                    if r["paralog"] == "PIEZO1" and r["domain"] == "pore_module"),
                   None)
        if row is None:
            raise KeyError("the census domain_map holds no PIEZO1 pore_module row")
        return set(range(int(row["ref_start"]), int(row["ref_end"]) + 1))
    ann = annotations or load_annotations("human")
    residues: set[int] = set()
    for dom in ann.domains:
        if dom.id in PORE_MODULE_DOMAINS and dom.start and dom.end:
            residues.update(range(dom.start, dom.end + 1))
    return residues


def boundary_disagreement() -> list[dict]:
    """Element by element, how far apart the two projects put each boundary.

    Reported because the pore-module test is sensitive to it and because the
    *pattern* is informative: the cap, inner helix and CTD agree to a few
    residues while the anchor and outer helix are shifted by more than a
    hundred, which is the signature of two ranges taken from different papers
    rather than of a systematic numbering error.
    """
    ours = {d.id: d for d in load_annotations("human").domains}
    pairs = (("anchor", "anchor"), ("outer_helix", "outer_helix"),
             ("CED", "cap"), ("inner_helix", "inner_helix"), ("CTD", "ctd"))
    rows = []
    for census_name, our_name in pairs:
        row = next((r for r in load_family_findings().table("domain_map")
                    if r["paralog"] == "PIEZO1" and r["domain"] == census_name),
                   None)
        dom = ours.get(our_name)
        if row is None or dom is None or dom.start is None:
            continue
        rows.append({
            "element": our_name,
            "census": (int(row["ref_start"]), int(row["ref_end"])),
            "ours": (dom.start, dom.end),
            "start_offset": dom.start - int(row["ref_start"]),
            "end_offset": dom.end - int(row["ref_end"]),
        })
    return rows


def _missense(record: dict) -> bool:
    """Only substitutions. A nonsense or frameshift position says nothing about
    where *this* residue matters — it removes everything downstream — and
    including them would put every truncating variant's position into whichever
    domain the stop codon happens to fall in."""
    wt, mut = record.get("wt_aa"), record.get("mut_aa")
    return bool(wt and mut and mut not in ("*", "-") and wt != mut
                and len(mut) == 1 and mut.isalpha())


def pathogenic_positions() -> set[int]:
    """Distinct residue numbers carrying a pathogenic missense variant.

    Both curated sets, pooled: ``variants.json`` (hand-curated, wild-type
    verified against Q92508) and ``variants_clinvar.json`` (ClinVar
    pathogenic/likely-pathogenic past the same gate). Directions are not used —
    the question is where disease *is*, not which way it pushes the channel.
    """
    positions: set[int] = set()
    curated = json.loads((RESOURCE_DIR / "variants.json").read_text())["variants"]
    for v in curated:
        if _missense(v) and v.get("classification") in ("GoF", "LoF"):
            positions.add(int(v["residue"]))
    clinvar = json.loads(
        (RESOURCE_DIR / "variants_clinvar.json").read_text())["variants"]
    for v in clinvar:
        if _missense(v) and str(v.get("significance", "")).lower().startswith(
                ("pathogenic", "likely pathogenic")):
            positions.add(int(v["residue"]))
    return positions


def population_positions(min_af: float | None = None,
                         offline: bool = True) -> tuple[set[int], str]:
    """Residue positions carrying population missense variation in gnomAD.

    Returns ``(positions, note)``. The note is not decoration: if the gnomAD
    cache is absent this returns an empty set and says so, rather than falling
    back to a length-based background — a background computed from the
    protein's own composition would make the test a statement about domain
    sizes and would still print a P value.
    """
    from .gnomad import GnomadClient, _protein_position

    floor = _P.value("family.population_af_floor") if min_af is None else min_af
    client = GnomadClient(offline=offline)
    variants = client.variants("PIEZO1")
    if not variants:
        return set(), ("gnomAD is not cached and this run is offline; no "
                       "population comparator available")
    positions: set[int] = set()
    for v in variants:
        if v.get("consequence") != "missense_variant":
            continue
        pos = _protein_position(v.get("hgvsp"))
        if pos is None:
            continue
        ac = an = 0
        for source in ("exome", "genome"):
            block = v.get(source) or {}
            ac += block.get("ac") or 0
            an += block.get("an") or 0
        if an and (ac / an) >= floor:
            positions.add(int(pos))
    return positions, (f"gnomAD r4 missense at allele frequency >= {floor:g}, "
                       f"{len(positions)} positions")


def pore_module_enrichment(annotations: Annotations | None = None,
                           offline: bool = True,
                           partition: str = "ours") -> Enrichment:
    """Does pathogenic missense concentrate in the pore module? Re-tested here."""
    ann = annotations or load_annotations("human")
    module = pore_module_residues(ann, partition=partition)
    # Measured, not written down. `Annotations` carries no sequence, and the
    # first draft here fell back to a literal 2521 — a protein length hard-coded
    # into an analysis module, which is the shape of error the parameter
    # registry exists to prevent and which the audit did not catch because it
    # sat inside a conditional expression.
    length = len(human_sequence())
    pathogenic = pathogenic_positions()
    comparator, note = population_positions(offline=offline)
    # A position carrying both a pathogenic variant and common population
    # variation belongs to the pathogenic side: the whole table is about where
    # disease is, and counting it twice would dilute the very contrast tested.
    comparator = comparator - pathogenic

    a = len(pathogenic & module)
    b = len(pathogenic - module)
    c = len(comparator & module)
    d = len(comparator - module)
    return Enrichment(
        region=f"pore module ({partition} boundaries)",
        n_region=len(module), n_protein=length,
        pathogenic_in=a, pathogenic_out=b, comparator_in=c, comparator_out=d,
        odds_ratio=_odds_ratio(a, b, c, d),
        p_value=fisher_exact_greater(a, b, c, d) if (c + d) else 1.0,
        comparator="gnomAD population missense",
        note=(f"PIEZO1 only, on the {partition} domain boundaries, against a "
              "population comparator rather than ClinVar benign labels. The "
              f"census pooled PIEZO1 and PIEZO2. {note}"))


def both_partitions(offline: bool = True) -> dict:
    """The enrichment under each partition, and whether the answer survives.

    The honest form of the result. If the two partitions disagree about
    significance then what has been measured is a boundary choice, and saying
    so is worth more than picking the one that agrees with the census.
    """
    results = {p: pore_module_enrichment(offline=offline, partition=p)
               for p in PARTITIONS}
    significant = {p: r.significant for p, r in results.items()}
    directions = {p: r.pathogenic_fraction > r.comparator_fraction
                  for p, r in results.items()}
    if all(significant.values()):
        verdict = "the enrichment holds on both partitions"
    elif not any(significant.values()):
        verdict = ("the enrichment does not reach significance on either "
                   "partition in PIEZO1 alone; the direction is "
                   + ("the same as the census's on both"
                      if all(directions.values()) else "not consistent"))
    else:
        held = [p for p, s in significant.items() if s]
        verdict = (f"the enrichment reaches significance only on the "
                   f"{', '.join(held)} boundaries, so what is being measured "
                   f"here is partly the boundary choice")
    # Where the two partitions differ, and what is in the gap. The disputed
    # band is not empty of disease — it is the reason the two answers differ —
    # so naming its pathogenic residues turns "boundary-dependent" from an
    # excuse into something a reader can check on the structure.
    ours_set = pore_module_residues(partition="ours")
    census_set = pore_module_residues(partition="census")
    disputed = census_set - ours_set
    pathogenic = pathogenic_positions()
    return {"results": results, "verdict": verdict,
            "boundaries": boundary_disagreement(),
            "disputed": {
                "n_residues": len(disputed),
                "span": (min(disputed), max(disputed)) if disputed else None,
                "pathogenic": sorted(pathogenic & disputed),
                "note": ("residues the census counts as pore module and this "
                         "project does not; on our boundaries they are the "
                         "anchor and the THU9 end of the blade"),
            },
            "census": census_comparison()}


def constraint_classifier(track: ConstraintTrack | None = None,
                          offline: bool = True) -> ClassifierResult | None:
    """The census constraint as a classifier on *our* pathogenic positions.

    The census reported AUC 0.914 separating its ClinVar pathogenic from its
    ClinVar benign positions. Here the positives are this project's curated
    pathogenic set and the negatives are gnomAD population positions, so a
    similar AUC is evidence the score generalises past the labels it was
    checked on, and a much lower one would say it had been reading the label
    set rather than the protein.
    """
    track = track or load_constraint("PIEZO1")
    positives = pathogenic_positions()
    negatives, note = population_positions(offline=offline)
    negatives = negatives - positives
    if not negatives or not positives:
        return None
    pos = np.array([v for v in (track.value(r) for r in sorted(positives))
                    if v is not None])
    neg = np.array([v for v in (track.value(r) for r in sorted(negatives))
                    if v is not None])
    if pos.size < 5 or neg.size < 5:
        return None
    census = next((r for r in load_family_findings().table("variant_constraint_auc")
                   if r["gene"] == "PIEZO1" and r["layer"] == "deep"), None)
    # auroc takes one score array and a boolean mask, not two arrays.
    scores = np.concatenate([pos, neg])
    is_positive = np.concatenate([np.ones(pos.size, bool), np.zeros(neg.size, bool)])
    return ClassifierResult(
        auc=float(auroc(scores, is_positive)), n_positive=int(pos.size),
        n_negative=int(neg.size), mean_positive=float(pos.mean()),
        mean_negative=float(neg.mean()), track=track.track,
        census_auc=float(census["auc"]) if census else None,
        note=("positives are this project's curated pathogenic missense; "
              f"negatives are {note}. The census's own figure used ClinVar "
              "benign labels as the negative set."))


def census_comparison() -> dict:
    """The census's own 2x2, read back so the two can be printed side by side."""
    finding = load_family_findings().by_key("disease_in_the_pore_module")
    return {"statement": finding.statement, "numbers": dict(finding.numbers),
            "caveat": finding.caveat, "source": finding.source}
