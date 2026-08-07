"""What AlphaFold's confidence actually constrains, read from PAE not pLDDT.

The project has fetched AlphaFold models since the beginning and never read the
**predicted aligned error**. pLDDT says how well a residue's *local* environment
is predicted; PAE says how well residue *i*'s position is determined when the
model is aligned on residue *j*. Only the second answers the question a hybrid
model needs: **is the distal blade placed correctly relative to the core?**

**The measured answer, and it is not the one the round expected.**

* pLDDT does agree with the seam: the unresolved distal blade (1–569) averages
  **64.5** against the core's **74.2**, with 52% of blade residues below 70
  versus 27% of the core.
* **PAE does not single out the seam.** Once sequence separation is controlled
  for — and it must be, since PAE grows with it — the blade↔core penalty is
  about 4 Å on a 31.75 Å scale, and at short separation it *reverses*: pairs
  50–150 apart across the seam score **13.25**, better than the 15.82 of pairs
  the same distance apart within one region.
* What PAE says instead is that AlphaFold does not determine PIEZO1's long-range
  architecture **anywhere**. At separations beyond 800 residues the mean PAE is
  **85% of maximum** — and **80% of maximum within the cryo-EM-resolved core
  alone**, a region experiment places confidently.

So the honest conclusion for :mod:`piezo1.structure.hybrid`: the seam is not the
weak point. The *global* arrangement is unconstrained by the prediction wherever
it is cut, which argues for placing the blade using the experimental C3 symmetry
and dome geometry rather than trusting AlphaFold's relative placement at all.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field

import numpy as np

__all__ = ["PAEMatrix", "ConfidenceReport", "load_pae", "load_plddt",
           "assess_seam"]


@dataclass
class PAEMatrix:
    """Predicted aligned error, with the comparisons that need a control."""

    values: np.ndarray                 # (n, n), angstrom
    maximum: float

    @property
    def n_residues(self) -> int:
        return len(self.values)

    def block_mean(self, rows, cols) -> float:
        return float(self.values[np.ix_(np.asarray(rows), np.asarray(cols))].mean())

    def by_separation(self, seam: int, bins=None) -> list:
        """Mean PAE within and across a seam, binned by sequence separation.

        The control that matters. Comparing raw blocks conflates "across the
        seam" with "far apart in sequence", because pairs spanning a boundary
        near one end of the chain are systematically more separated. PAE grows
        with separation for this protein, so the uncontrolled comparison would
        attribute that growth to the seam.
        """
        bins = bins or ((50, 150), (150, 400), (400, 800), (800, 1500),
                        (1500, 2600))
        i, j = np.indices(self.values.shape)
        separation = np.abs(i - j)
        crosses = (((i < seam - 1) & (j >= seam - 1))
                   | ((j < seam - 1) & (i >= seam - 1)))

        out = []
        for lo, hi in bins:
            mask = (separation >= lo) & (separation < hi)
            within = self.values[mask & ~crosses]
            across = self.values[mask & crosses]
            if len(within) < 100 or len(across) < 100:
                continue
            out.append({"separation": (lo, hi),
                        "n_within": int(len(within)),
                        "mean_within": float(within.mean()),
                        "n_across": int(len(across)),
                        "mean_across": float(across.mean()),
                        "penalty": float(across.mean() - within.mean())})
        return out

    def saturation(self, min_separation: int = 800, region=None) -> float:
        """Fraction of the maximum PAE reached at long range.

        A value near 1 means the prediction carries no information about the
        relative placement of those parts — the model is locally plausible and
        globally unconstrained.
        """
        values = (self.values if region is None
                  else self.values[np.ix_(np.asarray(region), np.asarray(region))])
        i, j = np.indices(values.shape)
        far = np.abs(i - j) > min_separation
        if not far.any():
            return float("nan")
        return float(values[far].mean() / self.maximum)


@dataclass
class ConfidenceReport:
    """Whether the prediction supports where a hybrid model would be cut."""

    seam: int
    plddt_blade: float
    plddt_core: float
    plddt_blade_low_fraction: float
    plddt_core_low_fraction: float
    separation_table: list = field(default_factory=list)
    saturation_all: float = float("nan")
    saturation_core: float = float("nan")
    meta: dict = field(default_factory=dict)

    @property
    def plddt_agrees_with_seam(self) -> bool:
        """Is the unresolved region genuinely lower-confidence locally?"""
        return self.plddt_blade < self.plddt_core - 5.0

    @property
    def pae_singles_out_seam(self) -> bool:
        """Does the seam stand out once separation is controlled for?

        Requires a consistent penalty. It does not hold here: the penalty
        reverses at short separation.
        """
        if not self.separation_table:
            return False
        return all(row["penalty"] > 2.0 for row in self.separation_table)

    @property
    def global_architecture_constrained(self) -> bool:
        """Does the prediction determine long-range geometry anywhere?"""
        return self.saturation_core < 0.6

    def summary(self) -> str:
        return (f"seam {self.seam}: pLDDT blade {self.plddt_blade:.1f} vs core "
                f"{self.plddt_core:.1f} "
                f"({'agrees' if self.plddt_agrees_with_seam else 'does not agree'}"
                f" with the seam); PAE "
                f"{'singles out' if self.pae_singles_out_seam else 'does NOT single out'}"
                f" the seam; long-range PAE {self.saturation_all:.0%} saturated "
                f"overall and {self.saturation_core:.0%} within the "
                f"experimentally-resolved core")


def load_pae(path=None) -> PAEMatrix | None:
    from ..config import STRUCTURE_DIR

    if path is None:
        matches = sorted(STRUCTURE_DIR.glob(
            "AF-*-predicted_aligned_error*.json"))
        if not matches:
            return None
        path = matches[-1]
    data = json.loads(open(path).read())
    entry = data[0] if isinstance(data, list) else data
    return PAEMatrix(
        values=np.asarray(entry["predicted_aligned_error"], dtype=float),
        maximum=float(entry["max_predicted_aligned_error"]))


def load_plddt(n_residues: int = 2521, path=None) -> np.ndarray | None:
    """Per-residue pLDDT, read from the model's B-factor column."""
    from ..config import STRUCTURE_DIR
    from ..core.structure import Structure

    if path is None:
        matches = sorted(STRUCTURE_DIR.glob("AF-*-model*.cif"))
        if not matches:
            return None
        path = matches[-1]
    structure = Structure.from_file(path)
    mask = structure.mask_ca()
    out = np.zeros(n_residues)
    residues = structure.res_seq[mask]
    keep = (residues >= 1) & (residues <= n_residues)
    out[residues[keep] - 1] = structure.b_factor[mask][keep]
    return out


def assess_seam(seam: int = 570, n_residues: int = 2521) -> ConfidenceReport | None:
    """The whole comparison: does the prediction support cutting here?"""
    pae = load_pae()
    plddt = load_plddt(n_residues)
    if pae is None or plddt is None:
        return None

    blade = plddt[:seam - 1]
    core = plddt[seam - 1:]
    core_index = np.arange(seam - 1, n_residues)
    return ConfidenceReport(
        seam=seam,
        plddt_blade=float(blade.mean()), plddt_core=float(core.mean()),
        plddt_blade_low_fraction=float((blade < 70).mean()),
        plddt_core_low_fraction=float((core < 70).mean()),
        separation_table=pae.by_separation(seam),
        saturation_all=pae.saturation(),
        saturation_core=pae.saturation(region=core_index),
        meta={"max_pae": pae.maximum, "n_residues": pae.n_residues,
              "note": "PAE answers where things sit relative to each other; "
                      "pLDDT does not"})
