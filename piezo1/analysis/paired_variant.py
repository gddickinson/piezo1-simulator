"""The one variant-versus-wild-type structural comparison this project can make.

Round 34 established that of four deposited variant entries only **8YFG
(R2456H)** resolves its own mutation and is coordinate-distinct. That leaves
exactly one pair worth measuring, and n = 1 supports no inference on its own.

**What makes a single pair interpretable is a control.** The question is not
"does R2456H differ from wild type?" — two structures always differ — but
"does it differ by more than wild-type entries differ among *themselves*?"
This project has three independent human entries that carry arginine at 2456,
so the wild-type-to-wild-type spread is measurable.

**The duplicate trap.** 8YFC and 9VMX have byte-identical coordinates to 8ZU3.
Including them would add two zero-difference pairs to the control, narrowing the
wild-type spread and making the variant look more distinct than it is. They are
excluded by coordinate fingerprint rather than by name, so a future duplicate is
caught automatically.

**Measured answer: R2456H is not distinguishable.** Its bottleneck radius and
wetting score both fall *inside* the wild-type range, and its largest difference
from any wild-type entry is smaller than the largest difference between two
wild-type entries — on both metrics. Which is not surprising once stated: every
deposited human structure is closed, and R2456H is a gating variant whose
phenotype is slowed inactivation. A closed structure need not show it.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["StructuralMetrics", "PairedComparison", "measure", "compare",
           "WILD_TYPE_CANDIDATES", "VARIANT_ENTRY"]

#: Human entries carrying arginine at 2456 — the wild-type reference set. Two
#: of these are coordinate duplicates and are removed at run time.
WILD_TYPE_CANDIDATES = ("8YEZ", "8ZU3", "8ZU8", "8YFC", "9VMX")

#: The only deposited variant entry that resolves its own mutation.
VARIANT_ENTRY = "8YFG"


@dataclass
class StructuralMetrics:
    """What can be measured on a closed structure and compared across entries."""

    pdb: str
    bottleneck_A: float
    wetting_score: float
    hydrophobic_gate: bool
    sterically_occluded: bool
    fingerprint: str = ""

    def as_vector(self) -> np.ndarray:
        return np.array([self.bottleneck_A, self.wetting_score], dtype=float)


@dataclass
class PairedComparison:
    """One variant against a wild-type set, with the within-set spread."""

    variant: StructuralMetrics
    wild_type: list = field(default_factory=list)
    excluded_duplicates: list = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def _spread(self, attribute: str) -> float:
        values = [getattr(m, attribute) for m in self.wild_type]
        return float(max(values) - min(values)) if len(values) > 1 else float("nan")

    def _largest_difference(self, attribute: str) -> float:
        value = getattr(self.variant, attribute)
        return float(max(abs(value - getattr(m, attribute))
                         for m in self.wild_type)) if self.wild_type else float("nan")

    def within_range(self, attribute: str) -> bool:
        """Does the variant fall inside the wild-type range for this measure?"""
        values = [getattr(m, attribute) for m in self.wild_type]
        value = getattr(self.variant, attribute)
        return bool(min(values) <= value <= max(values))

    def report(self) -> dict:
        out = {}
        for attribute in ("bottleneck_A", "wetting_score"):
            spread = self._spread(attribute)
            largest = self._largest_difference(attribute)
            out[attribute] = {
                "variant": float(getattr(self.variant, attribute)),
                "wild_type_range": [
                    float(min(getattr(m, attribute) for m in self.wild_type)),
                    float(max(getattr(m, attribute) for m in self.wild_type))],
                "wild_type_spread": spread,
                "largest_variant_difference": largest,
                "within_wild_type_range": self.within_range(attribute),
                "exceeds_wild_type_spread": bool(largest > spread),
            }
        return out

    @property
    def distinguishable(self) -> bool:
        """True only if the variant exceeds the wild-type spread on some measure.

        Deliberately generous: *any* measure would do. It is still false.
        """
        return any(v["exceeds_wild_type_spread"] for v in self.report().values())

    def summary(self) -> str:
        report = self.report()
        parts = []
        for name, values in report.items():
            parts.append(
                f"{name}: variant {values['variant']:.3f}, wild type "
                f"{values['wild_type_range'][0]:.3f}-"
                f"{values['wild_type_range'][1]:.3f} "
                f"(spread {values['wild_type_spread']:.3f}, largest variant "
                f"difference {values['largest_variant_difference']:.3f})")
        return (f"{self.variant.pdb} against {len(self.wild_type)} independent "
                f"wild-type entries; " + "; ".join(parts)
                + f". Distinguishable: {self.distinguishable}. "
                f"n = 1 — this supports no inference about R2456H, only a "
                f"statement about what the structures show.")


def measure(pdb: str) -> StructuralMetrics | None:
    """Pore and wetting metrics for one entry, in the canonical frame."""
    from ..core.structure import Structure
    from ..io.registry import load_registry
    from ..structure.frame import apply_frame, canonical_transform
    from ..structure.pore import pore_profile
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis
    from .hydration import load_grid, predict_wetting
    from .variant_structures import coordinate_fingerprint

    record = load_registry().get(pdb)
    if record is None or not record.available:
        return None
    raw = Structure.from_file(record.path)
    framed = apply_frame(raw, canonical_transform(raw))
    blocks, _ = protomer_blocks(framed)
    profile = pore_profile(framed, detect_c3_axis(blocks), step=1.0)

    grid = load_grid()
    if not grid.available:
        return None
    wetting = predict_wetting(framed, profile, grid)
    return StructuralMetrics(
        pdb=pdb, bottleneck_A=float(profile.bottleneck_radius),
        wetting_score=float(wetting.score),
        hydrophobic_gate=bool(wetting.hydrophobic_gate),
        sterically_occluded=bool(wetting.sterically_occluded),
        fingerprint=coordinate_fingerprint(raw))


def compare(variant: str = VARIANT_ENTRY,
            wild_type=WILD_TYPE_CANDIDATES) -> PairedComparison | None:
    """Measure the variant against the independent wild-type entries.

    Duplicates are removed by coordinate fingerprint, not by name: two entries
    with identical atoms contribute one comparison, and adding the second would
    narrow the control spread with a difference of exactly zero.
    """
    variant_metrics = measure(variant)
    if variant_metrics is None:
        return None

    seen, kept, duplicates = {}, [], []
    for pdb in wild_type:
        metrics = measure(pdb)
        if metrics is None:
            continue
        if metrics.fingerprint in seen:
            duplicates.append((pdb, seen[metrics.fingerprint]))
            continue
        seen[metrics.fingerprint] = pdb
        kept.append(metrics)

    return PairedComparison(
        variant=variant_metrics, wild_type=kept,
        excluded_duplicates=duplicates,
        meta={"note": "n = 1 variant structure; the wild-type spread is the "
                      "only thing that makes a single pair interpretable"})
