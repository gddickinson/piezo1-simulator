"""Pore hydration: does a sterically open pore actually conduct?

Radius alone is a poor predictor of conduction. A pore can be wide enough for a
hydrated ion and still block, because a hydrophobic neck expels liquid water —
**hydrophobic gating**. Rao et al. 2019 (PNAS 116:13989, PMID 31235590)
quantified this over ~200 channel structures and ~600 MD simulations, and found
the critical radius for wetting depends almost linearly on local hydrophobicity:
hydrophilic pores hydrate below 0.2 nm, while strongly hydrophobic ones can
present a barrier at radii up to ~0.4 nm.

Their result is a free-energy landscape over (hydrophobicity, radius), and the
heuristic reads:

1. compute the pore radius and a kernel-smoothed hydrophobicity profile;
2. for each pore-lining residue, look up the water free energy at its
   (hydrophobicity, radius);
3. flag residues above **1 RT = 2.6 kJ/mol**;
4. score = **sum of shortest distances** from the flagged points to the
   2.6 kJ/mol contour;
5. **Σd > 0.55 ⟹ closed.**

**We use the published landscape, not a redrawing of it.** The grid ships in the
CHAP repository under the **MIT licence**, so it can be used directly; guessing
a boundary off a figure would have been a silent correctness problem of exactly
the kind Round 17 dealt with. It is downloaded, not committed
(``python -m piezo1.io.fetch``), and every analysis degrades to "unavailable"
without it.

Reported AUROC for this heuristic is 0.91, against **0.59 for minimum radius
alone** — which is the whole reason the round exists.

**Units.** Coordinates are Å throughout this project; the CHAP grid is in nm.
The conversion happens at the boundary, in :meth:`HydrationGrid.energy`.
"""

from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np


from ..config import DERIVED_DIR
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["WIMLEY_WHITE_NORMALISED", "HydrationGrid", "load_grid",
           "hydrophobicity_profile_chap", "predict_wetting",
           "WettingPrediction", "LiningPoint", "CHAP_CITATION", "RAO_CITATION",
           "ENERGY_THRESHOLD_KJ", "CLOSED_SCORE_CUTOFF", "KERNEL_BANDWIDTH_NM",
           "WATER_RADIUS_NM", "pore_facing_residues"]

CHAP_CITATION = ("Klesse G, Rao S, Sansom MSP, Tucker SJ. CHAP: a versatile "
                 "tool for the structural and functional annotation of ion "
                 "channel pores. J Mol Biol 2019;431:3353-3365. PMID 31220459. "
                 "Software MIT licensed.")
RAO_CITATION = ("Rao S, Klesse G, Stansfeld PJ, Tucker SJ, Sansom MSP. A "
                "heuristic derived from analysis of the ion channel structural "
                "proteome permits the rapid identification of hydrophobic "
                "gates. PNAS 2019;116:13989-13995. PMID 31235590.")

#: 1 RT at room temperature, in kJ/mol. The contour CHAP draws.
ENERGY_THRESHOLD_KJ = _P.value("hydration.energy_threshold")
#: Sum-of-distances above which Rao et al. call a structure closed.
CLOSED_SCORE_CUTOFF = _P.value("hydration.closed_cutoff")
#: CHAP's default hydrophobicity kernel bandwidth, nm.
KERNEL_BANDWIDTH_NM = _P.value("hydration.kernel_bandwidth")
#: Radius of a water molecule, nm (Rao et al. 2019). Below this a pore is shut
#: for steric reasons and the wetting question does not arise.
WATER_RADIUS_NM = _P.value("hydration.water_radius")

#: Wimley–White whole-residue hydrophobicity, normalised to [-1, 1].
#: Transcribed from CHAP's ``share/data/hydrophobicity/wimley_white_1996.json``
#: (MIT) so that our hydrophobicity axis is the same axis the published grid was
#: built on. Using a different scale — Kyte–Doolittle, say — would index the
#: grid with the wrong coordinate and return confident nonsense.
WIMLEY_WHITE_NORMALISED: dict[str, float] = {
    "ALA": -0.13692946, "ARG": -0.41493776, "ASN": -0.17842324,
    "ASP": -1.00000000, "CYS": -0.09128631, "GLN": -0.07883817,
    "GLU": -0.66804979, "GLY": -0.47302905, "HIS": -0.56846473,
    "ILE": 0.33609959, "LEU": 0.33609959, "LYS": -0.75103734,
    "MET": 0.18257261, "PHE": 0.24066390, "PRO": 0.12863071,
    "SER": -0.13692946, "THR": -0.04564315, "TRP": 0.09958506,
    "TYR": -0.09543568, "VAL": 0.21991701,
}
#: CHAP's fallback for anything not in the scale (e.g. modified residues).
HYDROPHOBICITY_FALLBACK = 0.0


# --------------------------------------------------------------------------
# The published landscape
# --------------------------------------------------------------------------

@dataclass
class HydrationGrid:
    """Water free energy over (hydrophobicity, pore radius), from Rao 2019.

    ``hydrophobicity`` and ``radius`` are the unique axis values; ``energy`` is
    the ``(n_h, n_r)`` array in kJ/mol.
    """

    hydrophobicity: np.ndarray
    radius: np.ndarray             # nm
    energy: np.ndarray             # kJ/mol
    source: str = "CHAP heuristic_grid.json (MIT)"

    @property
    def available(self) -> bool:
        return self.energy.size > 0

    def energy_at(self, hydrophobicity, radius_angstrom) -> np.ndarray:
        """Nearest-grid-point lookup, matching CHAP's own scorer.

        ``radius_angstrom`` is in Å because that is this project's unit; the
        grid is in nm. Points outside the grid are clamped to its edge, which
        is deliberate: a radius of 1.2 nm is off the top of the landscape and
        unambiguously wet, and a nearest-edge lookup says so.
        """
        h = np.atleast_1d(np.asarray(hydrophobicity, dtype=float))
        r = np.atleast_1d(np.asarray(radius_angstrom, dtype=float)) / 10.0
        ih = np.abs(self.hydrophobicity[None, :] - h[:, None]).argmin(axis=1)
        ir = np.abs(self.radius[None, :] - r[:, None]).argmin(axis=1)
        out = self.energy[ih, ir]
        return out if out.size > 1 else out

    def contour(self, level: float = ENERGY_THRESHOLD_KJ) -> np.ndarray:
        """The (hydrophobicity, radius) curve at a given energy, as ``(n, 2)``.

        Extracted by linear interpolation along each hydrophobicity column,
        which is single-valued here because energy falls monotonically with
        radius at fixed hydrophobicity. Returned in nm to match the grid.
        """
        pts = []
        for i, h in enumerate(self.hydrophobicity):
            column = self.energy[i]
            crossings = np.flatnonzero((column[:-1] - level)
                                       * (column[1:] - level) <= 0)
            if crossings.size == 0:
                continue
            j = int(crossings[0])
            e0, e1 = column[j], column[j + 1]
            if e1 == e0:
                r = self.radius[j]
            else:
                frac = (level - e0) / (e1 - e0)
                r = self.radius[j] + frac * (self.radius[j + 1] - self.radius[j])
            pts.append((h, r))
        return np.array(pts, dtype=float)


def load_grid(path: Path | None = None) -> HydrationGrid:
    """Load the downloaded CHAP grid, or an empty one if it is not there."""
    path = Path(path or (DERIVED_DIR / "chap_heuristic_grid.json"))
    if not path.exists():
        return HydrationGrid(np.array([]), np.array([]), np.array([[]]),
                             source="not downloaded")
    raw = json.loads(path.read_text())
    h = np.array([e["hydrophobicity"] for e in raw], dtype=float)
    r = np.array([e["radius"] for e in raw], dtype=float)
    e = np.array([e["energy"] for e in raw], dtype=float)
    hs, rs = np.unique(h), np.unique(r)
    energy = np.full((hs.size, rs.size), np.nan)
    energy[np.searchsorted(hs, h), np.searchsorted(rs, r)] = e
    return HydrationGrid(hs, rs, energy)


# --------------------------------------------------------------------------
# Hydrophobicity along the pore
# --------------------------------------------------------------------------

def pore_facing_residues(structure: Structure, profile) -> list[tuple[int, str, float]]:
    """The pore-lining residues, each with its position along the pore axis.

    Returns ``(residue number, residue name, z)`` with ``z`` in Å, taken from
    the **side-chain** centroid where there is one — it is the side chain that
    faces the lumen and sets the local hydrophobicity, and using CA instead
    shifts a bulky residue like Phe several Å off its actual contribution.
    """
    axis = profile.axis
    prot = structure.mask_protein() & ~structure.hetero
    xyz = structure.xyz[prot].astype(np.float64)
    names = structure.res_name[prot]
    seqs = structure.res_seq[prot]
    chains = structure.chain[prot]
    atoms = structure.atom_name[prot]
    backbone = np.isin(atoms, ["N", "CA", "C", "O", "OXT"])

    lining: set[int] = set()
    for sl in profile.slices:
        lining.update(int(r) for r in sl.lining)

    out: list[tuple[int, str, float]] = []
    for chain in np.unique(chains):
        in_chain = chains == chain
        for residue in sorted(lining):
            m = in_chain & (seqs == residue)
            if not m.any():
                continue
            side = m & ~backbone
            pts = xyz[side] if side.any() else xyz[m]
            z = float(axis.project(pts.mean(axis=0)[None, :])[0])
            out.append((residue, str(names[m][0]), z))
    return out


def hydrophobicity_profile_chap(structure: Structure, profile,
                                bandwidth_nm: float = KERNEL_BANDWIDTH_NM
                                ) -> np.ndarray:
    """Kernel-smoothed pore hydrophobicity on the normalised Wimley–White scale.

    CHAP builds this as a Nadaraya–Watson average of the pore-facing residues'
    hydrophobicity, smoothed **along the pore coordinate** with a Gaussian
    kernel of bandwidth 0.35 nm.

    Smoothing in 3-D around the probe centre instead is the obvious shortcut
    and it destroys the measurement: a 1.85 nm neighbourhood pulls in the whole
    shell of residues surrounding the lumen, and the profile collapses to a
    narrow band near zero (−0.12 to +0.02 here) where the published grid spans
    −0.45 to +0.30. The landscape is then indexed by a coordinate that is not
    the one it was built on, and every energy read out of it is wrong while
    looking entirely reasonable.

    Returns one value per slice of ``profile``, NaN where no residue is near.
    """
    residues = pore_facing_residues(structure, profile)
    if not residues:
        return np.full(len(profile.z), np.nan)

    z_res = np.array([r[2] for r in residues])
    values = np.array([WIMLEY_WHITE_NORMALISED.get(r[1],
                                                   HYDROPHOBICITY_FALLBACK)
                       for r in residues])
    bandwidth = bandwidth_nm * 10.0            # nm -> Angstrom

    out = np.full(len(profile.z), np.nan)
    for i, z in enumerate(profile.z):
        weights = np.exp(-0.5 * ((z - z_res) / bandwidth) ** 2)
        total = weights.sum()
        if total > 1e-9:
            out[i] = float((weights * values).sum() / total)
    return out


# --------------------------------------------------------------------------
# The prediction
# --------------------------------------------------------------------------

@dataclass
class LiningPoint:
    """One pore-lining residue placed on the (hydrophobicity, radius) plane."""

    residue: int
    name: str
    z: float                  # position along the pore axis, Angstrom
    radius: float             # Angstrom
    hydrophobicity: float
    energy: float             # kJ/mol
    distance: float = 0.0     # to the 1 RT contour, 0 unless above threshold

    @property
    def above_threshold(self) -> bool:
        return self.energy > ENERGY_THRESHOLD_KJ


@dataclass
class WettingPrediction:
    """Verdict for one structure, with the evidence that produced it.

    **Two independent ways to be shut, kept apart on purpose.** The Rao
    heuristic answers "would water dewet here?", not "does water fit here?".
    Those are different questions, and PIEZO1 has structures that separate
    them: 7WLU and 8IXO have 0.098 nm bottlenecks — far too narrow for a water
    molecule — but hydrophilic linings, so their Σd is small and the heuristic
    alone would call them conductive. Reporting a single merged verdict would
    hide that; :attr:`hydrophobic_gate` and :attr:`sterically_occluded` are
    therefore both exposed and :attr:`conductive` requires neither.
    """

    score: float
    points: list[LiningPoint] = field(default_factory=list)
    min_radius: float = float("nan")      # Angstrom
    available: bool = True
    meta: dict = field(default_factory=dict)

    @property
    def hydrophobic_gate(self) -> bool:
        """Rao et al.'s rule: Σd above 0.55 means a closed hydrophobic gate."""
        return self.score > CLOSED_SCORE_CUTOFF

    @property
    def sterically_occluded(self) -> bool:
        """Narrower than a water molecule, so shut regardless of chemistry.

        Rao et al. quote the radius of a water molecule as ~0.15 nm, which is
        also the bottom of the range their landscape covers.
        """
        return bool(np.isfinite(self.min_radius)
                    and self.min_radius / 10.0 < WATER_RADIUS_NM)

    @property
    def conductive(self) -> bool:
        return not (self.hydrophobic_gate or self.sterically_occluded)

    @property
    def verdict(self) -> str:
        if not self.available:
            return "unavailable (CHAP grid not downloaded)"
        if self.conductive:
            return "conductive"
        reasons = []
        if self.sterically_occluded:
            reasons.append("sterically occluded")
        if self.hydrophobic_gate:
            reasons.append("hydrophobic gate")
        return "non-conductive (" + " + ".join(reasons) + ")"

    @property
    def dewetted(self) -> list[LiningPoint]:
        """The residues carrying the score, worst first — i.e. the gate."""
        return sorted((p for p in self.points if p.above_threshold),
                      key=lambda p: -p.distance)

    def summary(self) -> str:
        if not self.available:
            return self.verdict
        n = len(self.dewetted)
        worst = self.dewetted[0] if n else None
        text = (f"score {self.score:.2f} (cutoff {CLOSED_SCORE_CUTOFF}), "
                f"bottleneck {self.min_radius / 10:.3f} nm -> "
                f"{self.verdict}; {n} residue(s) above 1 RT")
        if worst is not None:
            text += (f", worst {worst.name}{worst.residue} "
                     f"r={worst.radius / 10:.2f} nm h={worst.hydrophobicity:+.2f}")
        return text


def predict_wetting(structure: Structure, profile,
                    grid: HydrationGrid | None = None,
                    max_radius: float = 7.0,
                    bandwidth_nm: float = KERNEL_BANDWIDTH_NM) -> WettingPrediction:
    """Apply the Rao et al. 2019 heuristic to a pore profile.

    ``max_radius`` (Å) mirrors CHAP's own restriction of the score to residues
    lining the narrow part of the pore — 0.7 nm there. Wide vestibules are
    always wet and would only dilute the score.
    """
    grid = grid if grid is not None else load_grid()
    if not grid.available:
        return WettingPrediction(score=float("nan"), available=False,
                                 meta={"reason": "CHAP grid not downloaded; "
                                                 "run python -m piezo1.io.fetch"})

    hydro = hydrophobicity_profile_chap(structure, profile,
                                        bandwidth_nm=bandwidth_nm)
    contour = grid.contour()          # (n, 2) in (hydrophobicity, nm)

    points: list[LiningPoint] = []
    seen: set[int] = set()
    for i, sl in enumerate(profile.slices):
        if not np.isfinite(hydro[i]) or sl.radius > max_radius:
            continue
        energy = float(np.atleast_1d(grid.energy_at(hydro[i], sl.radius))[0])
        for residue, name in zip(sl.lining, sl.lining_names):
            if residue in seen:
                continue
            seen.add(int(residue))
            pt = LiningPoint(residue=int(residue), name=str(name),
                             z=float(sl.z), radius=float(sl.radius),
                             hydrophobicity=float(hydro[i]), energy=energy)
            if pt.above_threshold:
                pt.distance = _distance_to_contour(
                    hydro[i], sl.radius / 10.0, contour)
            points.append(pt)

    score = float(sum(p.distance for p in points))
    return WettingPrediction(
        score=score, points=points,
        min_radius=float(np.min(profile.radius)) if len(profile.radius) else float("nan"),
        meta={"n_lining": len(points), "grid": grid.source,
              "bandwidth_nm": bandwidth_nm, "max_radius_A": max_radius,
              "threshold_kJ": ENERGY_THRESHOLD_KJ,
              "cutoff": CLOSED_SCORE_CUTOFF, "citation": RAO_CITATION})


def _distance_to_contour(hydrophobicity: float, radius_nm: float,
                         contour: np.ndarray) -> float:
    """Shortest distance from a point to the 1 RT contour, in grid units.

    CHAP projects orthogonally onto the nearest contour segment. Distance is
    measured in the plotted (hydrophobicity, radius-in-nm) plane, which is what
    makes the 0.55 cutoff dimensionally meaningful — the two axes are treated
    as comparable, and rescaling either would silently change the verdict.
    """
    if contour.size == 0:
        return 0.0
    point = np.array([hydrophobicity, radius_nm])
    best = float(np.hypot(*(contour - point).T).min())
    # Refine against the segments, not just the vertices.
    a, b = contour[:-1], contour[1:]
    seg = b - a
    length2 = (seg ** 2).sum(axis=1)
    ok = length2 > 1e-15
    if ok.any():
        t = np.clip(((point - a[ok]) * seg[ok]).sum(axis=1) / length2[ok], 0, 1)
        proj = a[ok] + t[:, None] * seg[ok]
        best = min(best, float(np.hypot(*(proj - point).T).min()))
    return best
