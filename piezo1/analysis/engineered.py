"""The fifteen engineered variants: what each measures, and whether it counts.

Round 54 marked these `blocked` on a scientific question rather than a data
one: **may a change in conductance or selectivity stand for gain or loss of
mechanosensitive function?** They are the only measured functional effects the
project holds that no analysis set uses, so leaving the question unanswered was
leaving evidence on the floor without saying why.

**The answer is no**, and it is not an argument from caution — this project's
own curated annotations demonstrate that the axes dissociate:

- **A2078W**: "Yoda1 sensitivity severely reduced *while mechanosensitivity to
  stretch is retained*". Chemical agonist response and mechanical activation
  come apart at a single residue.
- **KKKK2166-**: "selectively removes inactivation *without changing mechanical
  sensitivity*". Inactivation and mechanosensitivity come apart too.

So a variant can lose most of what an assay measures while its mechanosensation
is untouched. The disease phenotypes are about mechanotransduction —
gain-of-function hereditary xerocytosis variants show slowed inactivation and
larger mechanically-evoked currents; loss-of-function lymphatic dysplasia
variants show reduced or absent ones. A residue in the selectivity filter that
halves unitary conductance changes how much current flows *once the channel is
open*, which is a different question from how readily force opens it.

**But five of the fifteen do carry a mechanosensitivity direction**, and those
are admissible in principle: S1335A, S1335V, A1718W and P2113A all raise the
threshold or desensitise the channel mechanically, and S2446E stabilises an open
intermediate.

**And admitting them changes nothing**, which is the measurement that settles
the round rather than the reasoning. None sits at a position carrying a
directional curated variant. Position 1335 does hold a pair — S1335A and
S1335V — but *both raise the mechanical threshold*, so it is a same-direction
pair and a within-position test needs a discriminating one. Admitting all five
adds **zero** discriminating positions to the one Round 62 counted.

If they were ever admitted they would enter at their own evidence level,
never pooled with `measured`; `variant_sets` already refuses that.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["EngineeredVerdict", "ENGINEERED", "AXES", "by_axis",
           "mechanosensitivity_variants", "admitting_would_add", "decision"]

#: What an assay can measure, and whether it is the axis the disease phenotypes
#: are defined on. Only one of these is.
AXES = {
    "mechanosensitivity": ("force threshold or stretch response — the axis "
                           "gain/loss of function is defined on", True),
    "agonist_only": ("response to Yoda1 or Jedi, a chemical agonist", False),
    "permeation_only": ("unitary conductance, selectivity or pore block", False),
    "inactivation_only": ("inactivation kinetics with mechanics unchanged", False),
    "structural_probe": ("cysteine accessibility — a topology tool, not a "
                         "phenotype", False),
}


@dataclass(frozen=True)
class EngineeredVerdict:
    """One engineered variant, and what its measured effect is *about*."""

    label: str
    axis: str
    basis: str                 # the curated functional_effect it rests on

    def __post_init__(self) -> None:
        if self.axis not in AXES:
            raise ValueError(f"unknown axis {self.axis!r}")

    @property
    def is_a_direction(self) -> bool:
        return AXES[self.axis][1]


ENGINEERED: tuple = (
    # --- the axis the disease phenotypes are defined on
    EngineeredVerdict("S1335A", "mechanosensitivity",
                      "more force required to open"),
    EngineeredVerdict("S1335V", "mechanosensitivity",
                      "raised mechanical activation threshold"),
    EngineeredVerdict("A1718W", "mechanosensitivity",
                      "stretch-activated currents lost"),
    EngineeredVerdict("P2113A", "mechanosensitivity",
                      "desensitises the channel mechanically"),
    EngineeredVerdict("S2446E", "mechanosensitivity",
                      "stabilises an intermediate open state"),

    # --- chemical agonist only, and one of them proves the dissociation
    EngineeredVerdict("A2075W", "agonist_only",
                      "Yoda1 sensitivity partially diminished"),
    EngineeredVerdict("A2078W", "agonist_only",
                      "Yoda1 sensitivity severely reduced WHILE "
                      "mechanosensitivity to stretch is retained"),

    # --- what happens once the channel is already open
    EngineeredVerdict("E2117A", "permeation_only", "halves unitary conductance"),
    EngineeredVerdict("E2117D", "permeation_only",
                      "increases unitary conductance and ruthenium-red affinity"),
    EngineeredVerdict("E2117K", "permeation_only",
                      "abolishes ruthenium-red block, lowers conductance, "
                      "increases anion permeability"),
    EngineeredVerdict("E2117Q", "permeation_only", "reduces Ca2+ selectivity"),
    EngineeredVerdict("E2470X", "permeation_only",
                      "alters Ca2+ selectivity, unitary conductance and RR "
                      "blockade"),

    # --- the second dissociation
    EngineeredVerdict("KKKK2166-", "inactivation_only",
                      "selectively removes inactivation WITHOUT changing "
                      "mechanical sensitivity"),

    # --- tools, not phenotypes
    EngineeredVerdict("Y2438C", "structural_probe",
                      "cysteine-accessible, establishing the last TM as the "
                      "pore-lining inner helix"),
    EngineeredVerdict("M2441C", "structural_probe",
                      "cysteine-accessible pore-facing position"),
)


def by_axis(axis: str) -> list:
    return [v for v in ENGINEERED if v.axis == axis]


def mechanosensitivity_variants() -> list:
    """The five whose measured effect is on the axis that defines direction."""
    return [v for v in ENGINEERED if v.is_a_direction]


def admitting_would_add() -> dict:
    """Discriminating positions gained by admitting the admissible five.

    The measurement that settles the round. Reasoning about whether an assay
    *may* stand in is worth writing down; whether it would *help* is worth
    measuring, and it does not.
    """
    from ..core.annotations import load_annotations
    from .data_routes import DIRECTIONS, _variants_by_position

    curated = {v.label: v for v in load_annotations("human").variants}
    by_position, _ = _variants_by_position()

    gained, same_direction = [], []
    for verdict in mechanosensitivity_variants():
        variant = curated.get(verdict.label)
        if variant is None:
            continue
        neighbours = {l: d for l, d in by_position.get(variant.residue, {}).items()
                      if d in DIRECTIONS}
        siblings = [v for v in mechanosensitivity_variants()
                    if v.label != verdict.label
                    and curated.get(v.label)
                    and curated[v.label].residue == variant.residue]
        if neighbours:
            gained.append(variant.residue)
        elif siblings:
            same_direction.append(variant.residue)

    return {"discriminating_positions_gained": sorted(set(gained)),
            "same_direction_pairs": sorted(set(same_direction)),
            "n_admissible": len(mechanosensitivity_variants())}


def decision() -> dict:
    """The round's answer, as data."""
    added = admitting_would_add()
    return {
        "may_permeation_stand_for_direction": False,
        "why": ("A2078W loses Yoda1 sensitivity while retaining stretch "
                "response, and KKKK2166- removes inactivation without changing "
                "mechanical sensitivity — the project's own annotations show "
                "the axes dissociate, so an assay on one does not report the "
                "other."),
        "admissible_in_principle": [v.label for v in mechanosensitivity_variants()],
        "would_add_discriminating_positions":
            len(added["discriminating_positions_gained"]),
        "note": ("Admitting the five changes nothing measurable: none sits at a "
                 "position with a directional curated variant, and position "
                 "1335's pair is same-direction. If they were ever admitted "
                 "they would enter at their own evidence level."),
    }
