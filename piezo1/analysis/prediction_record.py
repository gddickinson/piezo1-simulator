"""What a variant prediction from this project is actually worth.

The mechanical ΔΔG is the project's central claim, and it has been tested three
times under pre-registration and failed three times. Until now that lived in
`docs/` and the CLI, so a user could select a variant in the application, see a
number, and have no way of knowing what the number is entitled to claim.

This module is the honest answer, as data rather than prose, so the same record
drives the GUI, the CLI and the tests and they cannot drift apart. It is
deliberately Qt-free.

**The record is frozen.** Each entry is a result that was pre-registered and
then recorded; per `NEGATIVE_RESULT_PROTOCOL.md` a recorded result is superseded
by a new entry, never edited. :func:`verify_record` re-reads the stored Round 36
result and fails if the numbers here have drifted from it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

__all__ = ["ValidationEntry", "VALIDATION_RECORD", "headline", "what_it_means",
           "evidence_levels", "verify_record", "variant_evidence"]


@dataclass(frozen=True)
class ValidationEntry:
    """One pre-registered test of the central claim, and what it concluded."""

    round: int
    predictor: str
    n_gof: int
    n_lof: int
    cliffs_delta: float
    p_value: float | None
    power_at_large: float | None
    conclusion: str
    document: str

    @property
    def n(self) -> int:
        return self.n_gof + self.n_lof

    @property
    def rejected(self) -> bool:
        return False              # all three failed to reject; see `conclusion`

    def summary(self) -> str:
        p = "—" if self.p_value is None else f"{self.p_value:.3f}"
        return (f"Round {self.round}: {self.predictor}, "
                f"{self.n_gof} GoF vs {self.n_lof} LoF, "
                f"Cliff's δ {self.cliffs_delta:+.3f}, p = {p} — "
                f"{self.conclusion}")


#: The three tests, in order. Numbers taken from the recorded documents, not
#: recomputed here — `verify_record` is what keeps them honest.
VALIDATION_RECORD: tuple[ValidationEntry, ...] = (
    ValidationEntry(
        round=7, predictor="elastic-network ΔΔG (volume-based)",
        n_gof=16, n_lof=9, cliffs_delta=-0.083, p_value=0.234,
        power_at_large=0.60, conclusion="failed to reject",
        document="docs/VALIDATION.md"),
    ValidationEntry(
        round=22, predictor="FoldX ΔΔG (external)",
        n_gof=20, n_lof=6, cliffs_delta=-0.211, p_value=None,
        power_at_large=None, conclusion="failed to reject; nothing survived BH",
        document="docs/VALIDATION_ROUND22.md"),
    ValidationEntry(
        round=36, predictor="substitution-aware mechanical ΔΔG",
        n_gof=19, n_lof=15, cliffs_delta=-0.249, p_value=0.405,
        power_at_large=0.84, conclusion="failed to reject; nothing survived BH",
        document="docs/VALIDATION_ROUND36.md"),
)

#: What the two evidence levels mean, in the words a user needs.
EVIDENCE_LEVELS = {
    "measured": ("Direction determined by electrophysiology — the strongest "
                 "evidence this project holds."),
    "disease_mechanism": ("Direction *inferred* from which disease the variant "
                          "causes, not measured. Weaker: it assumes the variant "
                          "acts by the usual mechanism for that disease."),
    "none": "No direction is claimed for this variant.",
}


def evidence_levels() -> dict:
    return dict(EVIDENCE_LEVELS)


def headline() -> str:
    """One sentence, for putting next to a number in the interface."""
    latest = VALIDATION_RECORD[-1]
    return (f"Three pre-registered tests, three nulls. Latest (Round "
            f"{latest.round}): Cliff's δ {latest.cliffs_delta:+.3f}, "
            f"p = {latest.p_value:.3f} on {latest.n} variants. This score does "
            f"not predict gain- versus loss-of-function.")


def what_it_means() -> list[str]:
    """The standing caveats, in the order a reader needs them.

    Kept as a list rather than a paragraph so the interface can show the first
    two beside a number and the rest on request, without anyone rewriting them
    into something softer.
    """
    return [
        "The mechanical ΔΔG has been tested three times under pre-registration "
        "and has failed to separate gain- from loss-of-function every time.",
        "The Round 36 design had 84% power for a LARGE effect and 50% for a "
        "medium one, so the null excludes a large effect and not a medium one.",
        "The point estimate has grown across the three tests (−0.083, −0.211, "
        "−0.249) which is suggestive and is NOT evidence: at that effect size "
        "134 variants would be needed and 34 were available.",
        "Round 47 costed the ceiling: even if every variant this project could "
        "ever curate were included, 59 would be reachable, where the smallest "
        "detectable effect is 0.356 against the 0.249 observed. More data of "
        "the kind that could exist would not settle this.",
        "Two other predictor families were pre-registered and also failed: "
        "population constraint (Round 41) and wild-type structural context — "
        "burial, conservation, gate coupling (Round 48). Five tests, five "
        "nulls, five predictors.",
        "The binding constraint is data, not method. Only one deposited variant "
        "structure resolves its own mutation, and all four are gain-of-function.",
        "Use the score to ask which residues sit in mechanically coupled "
        "positions. Do not use it to assign a direction to a variant.",
    ]


def variant_evidence(label: str) -> dict:
    """Everything known about one variant's *evidence*, not its structure.

    Returns the recorded direction, which evidence level supports it, and
    whether the curated and ClinVar-inferred directions disagree — the last
    because a disagreement is the strongest signal that a label is unreliable,
    and it is otherwise buried in a CLI report.
    """
    from .variant_sets import build_analysis_set, disagreements

    combined = build_analysis_set(levels=("measured", "disease_mechanism"))
    entry = next((e for e in combined.entries if e.label == label), None)
    conflict = next((d for d in disagreements() if d["label"] == label), None)

    if entry is None:
        return {"label": label, "direction": None, "evidence": "none",
                "evidence_note": EVIDENCE_LEVELS["none"], "conflict": conflict,
                "in_analysis_set": False}
    return {
        "label": label,
        "direction": entry.classification,
        "evidence": entry.evidence,
        "evidence_note": EVIDENCE_LEVELS.get(entry.evidence, ""),
        "citation": entry.citation,
        "conflict": conflict,
        "in_analysis_set": True,
    }


def verify_record() -> dict:
    """Check the frozen numbers against the stored Round 36 result.

    A record that drifts from the run that produced it is worse than no record,
    because it looks like provenance. This is the same discipline
    `analysis.claims` applies to the documented numbers.
    """
    import json

    from ..config import DERIVED_DIR

    path = DERIVED_DIR / "validation_round36.json"
    if not path.exists():
        return {"checked": False,
                "reason": "run scripts/run_validation_round36.py first"}

    stored = json.loads(path.read_text())["primary"]
    entry = next(e for e in VALIDATION_RECORD if e.round == 36)
    drift = {
        "cliffs_delta": abs(stored["cliffs_delta"] - entry.cliffs_delta),
        "p_value": abs(stored["p_value"] - (entry.p_value or 0.0)),
        "n_gof": abs(stored["n_gof"] - entry.n_gof),
        "n_lof": abs(stored["n_lof"] - entry.n_lof),
    }
    ok = (drift["cliffs_delta"] < 0.005 and drift["p_value"] < 0.02
          and drift["n_gof"] == 0 and drift["n_lof"] == 0)
    return {"checked": True, "agrees": ok, "drift": drift}


@dataclass
class PredictionContext:
    """A variant's mechanical score, wrapped in what it is worth."""

    label: str
    ddg_gating: float | None = None
    modelled: bool = False
    evidence: dict = field(default_factory=dict)
    caveats: list = field(default_factory=list)
    note: str = ""

    def summary(self) -> str:
        if not self.modelled:
            return f"{self.label}: residue not modelled — no prediction."
        return (f"{self.label}: ΔΔG(gating) = {self.ddg_gating:+.4g}. "
                f"{headline()}")
