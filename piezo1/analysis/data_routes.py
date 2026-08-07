"""What each remaining route to more phenotyped variants would actually yield.

Every test of the central claim has ended in "not enough phenotyped variants".
That sentence is not actionable. This module turns it into counts: for each
route that needs no new electrophysiology, how many usable variants it adds,
what it would cost, and — where the route is worth taking — **which specific
variants to curate**.

The two routes the roadmap named have both already been run. Round 41 tested
population constraint and returned null; Round 45 harvested the open-access
corpus and found 35 candidates, none carrying a direction. So the costing here
is of what is left.

**The correction this module exists to record.** The Round 50 review counted 40
positions carrying more than one variant and called that "a real design". That
count was wrong for the purpose: it included nonsense variants, insertions,
positions where the second variant has no direction, and one position where two
sources disagree. Filtered to what a within-position comparison actually needs —
**two or more missense variants at one position, each with a direction, with no
source conflict** — the count is **one**. Round 48's original figure was right.

That matters because a within-position design was the one route Round 47 left
open, and its optimism came from a number that does not survive its own
filters.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["Route", "CurationTarget", "usable_positions", "one_label_away",
           "assess_routes", "DIRECTIONS"]

DIRECTIONS = frozenset({"GoF", "LoF"})


@dataclass(frozen=True)
class CurationTarget:
    """One variant whose direction, if resolved, would unlock a position."""

    position: int
    needs: str               # the variant lacking a direction
    unlocks_with: str        # the variant at that position that already has one
    existing_direction: str

    def summary(self) -> str:
        return (f"{self.needs} — position {self.position} already has "
                f"{self.unlocks_with} ({self.existing_direction})")


@dataclass
class Route:
    """One way to get more usable variants, with its measured yield."""

    name: str
    yields: int
    cost: str
    status: str              # done | open | blocked
    note: str = ""
    targets: list = field(default_factory=list)

    def summary(self) -> str:
        return f"[{self.status}] {self.name}: +{self.yields} — {self.cost}"


def _missense(label, wt, mut, kind=None) -> bool:
    """A single-residue substitution, excluding nonsense and indels.

    A nonsense variant truncates the protein and an insertion changes its
    length; neither is a substitution a structural predictor can score against
    another substitution at the same position. Both were being counted before.
    """
    if kind is not None and kind != "missense":
        return False
    if not label or not wt or not mut:
        return False
    return (len(str(wt)) == 1 and len(str(mut)) == 1
            and "*" not in str(label) and "-" not in str(label))


def _variants_by_position() -> tuple[dict, set]:
    """``{position: {label: direction}}`` over both sources, plus conflicts.

    A label appearing in both sources with different directions is recorded as
    a conflict and its position excluded — an unresolved disagreement is not
    evidence, and `variant_sets.disagreements()` already surfaces it for the
    reader.
    """
    import json

    from ..config import RESOURCE_DIR
    from ..core.annotations import load_annotations

    by_position: dict[int, dict] = {}
    conflicts: set = set()

    for variant in load_annotations("human").variants:
        if _missense(variant.label, variant.wt_aa, variant.mut_aa):
            by_position.setdefault(variant.residue, {})[variant.label] = \
                variant.classification

    payload = json.loads((RESOURCE_DIR / "variants_clinvar.json").read_text())
    entries = payload["variants"] if isinstance(payload, dict) else payload
    for entry in entries:
        if not _missense(entry.get("label"), entry.get("wt_aa"),
                         entry.get("mut_aa"), entry.get("kind")):
            continue
        position, label = entry.get("residue"), entry.get("label")
        direction = entry.get("classification")
        existing = by_position.setdefault(position, {}).get(label)
        if (existing in DIRECTIONS and direction in DIRECTIONS
                and existing != direction):
            conflicts.add((position, label))
        elif existing is None:
            by_position[position][label] = direction

    return by_position, conflicts


def usable_positions() -> dict:
    """Positions a within-position comparison could actually use.

    Two or more missense variants, each carrying a direction, at a position
    with no source conflict. Measured: **one**.
    """
    by_position, conflicts = _variants_by_position()
    conflicted = {p for p, _ in conflicts}
    out = {}
    for position, labels in by_position.items():
        if position in conflicted:
            continue
        directional = {l: d for l, d in labels.items() if d in DIRECTIONS}
        if len(directional) >= 2:
            out[position] = directional
    return out


def discriminating_positions() -> dict:
    """Usable positions carrying **both** directions — the informative ones."""
    return {p: d for p, d in usable_positions().items()
            if len(set(d.values())) > 1}


def one_label_away() -> list:
    """Variants whose direction would unlock a new within-position pair.

    This is the actionable output: not "we need more variants" but a named,
    countable list of which ones.
    """
    by_position, conflicts = _variants_by_position()
    conflicted = {p for p, _ in conflicts}
    targets = []
    for position, labels in sorted(by_position.items()):
        if position in conflicted or len(labels) < 2:
            continue
        have = [(l, d) for l, d in labels.items() if d in DIRECTIONS]
        lack = [l for l, d in labels.items() if d not in DIRECTIONS]
        if len(have) == 1 and lack:
            label, direction = have[0]
            for missing in sorted(lack):
                targets.append(CurationTarget(
                    position=position, needs=missing, unlocks_with=label,
                    existing_direction=direction))
    return targets


def assess_routes() -> list:
    """Every route that needs no new electrophysiology, with its yield."""
    from ..core.annotations import load_annotations

    curated = load_annotations("human").variants
    engineered = [v for v in curated if v.classification == "engineered"]
    targets = one_label_away()
    usable = usable_positions()

    return [
        Route(name="Population constraint (gnomAD)", yields=0,
              cost="already spent", status="done",
              note="Round 41 ran it and returned null; the pre-registered "
                   "negative control was indistinguishable from the predictor. "
                   "Constraint is not a direction and cannot become one."),
        Route(name="Published supplementary tables (literature harvest)",
              yields=0, cost="already spent", status="done",
              note="Round 45 found 35 candidate substitutions not in the "
                   "curated set. None carries a direction, and only two have "
                   "any electrophysiological measurement behind them."),
        Route(name="Within-position pairs, as they stand",
              yields=len(usable), cost="no curation — already available",
              status="open",
              note="Measured: one position, R2456. The Round 50 review said 40, "
                   "which counted nonsense variants, insertions, undirected "
                   "second variants and one source conflict."),
        Route(name="Curate the variants one label from a pair",
              yields=len(targets), cost=f"{len(targets)} literature reads",
              status="open", targets=targets,
              note="The cheapest route with any yield, and the only one that "
                   "is a named finite list rather than a search. The yield is "
                   "an UPPER BOUND: two of the three are curated as VUS "
                   "precisely because the evidence to direct them was not "
                   "found, and Round 45 measured how thin that literature is. "
                   "Even at full yield this reaches four positions, which is "
                   "not a design."),
        Route(name="Admit the engineered variants", yields=len(engineered),
              cost="no curation — a decision about admissibility",
              status="blocked",
              note="All 15 carry a measured functional effect, but the effect "
                   "is a change in conductance or selectivity, not gain or "
                   "loss of mechanosensitive function. Whether one may stand "
                   "for the other is a scientific question, not a data one."),
    ]


# --------------------------------------------------------------------------
# The same question, split by how the direction was established
# --------------------------------------------------------------------------
#
# Round 62. `variant_sets` refuses to pool `measured` (electrophysiology) with
# `disease_mechanism` (inferred from which disease the variant causes), so a
# within-position design has a different ceiling at each level. Counting only
# the pooled total would overstate what a confirmatory test could use.

def positions_by_evidence() -> dict:
    """Usable and reachable shared positions, per evidence level.

    Returns ``{level: {"usable": [...], "reachable": [...]}}`` where *usable*
    means the position already carries both directions at that level, and
    *reachable* adds the positions one curated label away whose existing
    partner is at that level — a pair cannot be stronger than its weaker half.
    """
    from .variant_sets import build_analysis_set

    out = {}
    for levels in (("measured",), ("measured", "disease_mechanism")):
        entries = build_analysis_set(levels=levels).missense().entries
        by_position: dict[int, dict] = {}
        for entry in entries:
            by_position.setdefault(entry.residue, {})[entry.label] = \
                entry.classification
        usable = sorted(p for p, labels in by_position.items()
                        if len({d for d in labels.values()
                                if d in DIRECTIONS}) > 1)

        # A target only reaches this level if the partner it joins is at it.
        at_level = {e.label for e in entries}
        extra = sorted({t.position for t in one_label_away()
                        if t.unlocks_with in at_level})

        out["+".join(levels)] = {
            "usable": usable,
            "reachable": sorted(set(usable) | set(extra)),
        }
    return out


def evidence_summary() -> dict:
    """The Round 62 answer as counts, against the Round 61 requirement."""
    from .feasibility import paired_positions_required

    levels = positions_by_evidence()
    # The most optimistic requirement Round 61 computed: delta 0.8, meaning the
    # predictor orders nine pairs in ten correctly.
    best_case = paired_positions_required(0.8, n_simulations=1500).positions
    return {
        "by_level": {name: {"usable": len(v["usable"]),
                            "reachable": len(v["reachable"])}
                     for name, v in levels.items()},
        "best_case_requirement": best_case,
        "reachable_anywhere": max(len(v["reachable"]) for v in levels.values()),
    }
