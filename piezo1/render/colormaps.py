"""Colour schemes for molecular representations.

Every scheme maps a :class:`~piezo1.core.structure.Structure` (plus optional
annotation) onto an ``(n_atoms, 3)`` float array in 0-1 RGB. Schemes are
deliberately dark-background-first: the viewport is near-black, so mid-to-light
saturated hues read best and pure white is reserved for highlights.
"""

from __future__ import annotations

import colorsys
import json
from dataclasses import dataclass

import numpy as np

from ..config import RESOURCE_DIR
from ..core.structure import Structure

__all__ = ["hex_to_rgb", "chain_colors", "domain_colors", "bfactor_colors",
           "plddt_colors", "value_colors", "uniform_color", "SEQUENCE_COLORS",
           "DomainPalette", "load_domain_palette", "PLDDT_BANDS"]


def hex_to_rgb(value: str) -> tuple[float, float, float]:
    v = value.lstrip("#")
    return tuple(int(v[i:i + 2], 16) / 255.0 for i in (0, 2, 4))  # type: ignore[return-value]


#: Chain palette — distinguishable, and consistent between light and dark.
SEQUENCE_COLORS = [
    (0.31, 0.60, 0.94), (0.95, 0.55, 0.28), (0.44, 0.80, 0.45),
    (0.90, 0.40, 0.45), (0.66, 0.52, 0.92), (0.35, 0.78, 0.80),
    (0.93, 0.75, 0.30), (0.80, 0.50, 0.72),
]

#: AlphaFold's own pLDDT bands, so our colouring matches the AFDB viewer.
PLDDT_BANDS = [
    (90.0, (0.05, 0.34, 0.75)),   # very high
    (70.0, (0.40, 0.79, 0.94)),   # confident
    (50.0, (0.99, 0.86, 0.36)),   # low
    (0.0, (0.99, 0.49, 0.25)),    # very low
]


def uniform_color(structure: Structure, color=(0.55, 0.62, 0.75)) -> np.ndarray:
    return np.tile(np.asarray(color, np.float32), (structure.n_atoms, 1))


def chain_colors(structure: Structure) -> np.ndarray:
    out = np.zeros((structure.n_atoms, 3), np.float32)
    for i, ch in enumerate(structure.chains):
        out[structure.chain == ch] = SEQUENCE_COLORS[i % len(SEQUENCE_COLORS)]
    return out


# --------------------------------------------------------------------------
# Domains
# --------------------------------------------------------------------------

@dataclass
class DomainPalette:
    """Domain definitions with colours, in one species' numbering."""

    species: str
    domains: list[dict]

    def color_of(self, residue: int) -> tuple[float, float, float] | None:
        for d in self.domains:
            span = d[self.species]
            if span["start"] is None or span["end"] is None:
                continue
            if span["start"] <= residue <= span["end"]:
                return hex_to_rgb(d["color"])
        return None

    def lookup_table(self, length: int) -> np.ndarray:
        """Per-residue colour table indexed by residue number (1-based)."""
        table = np.tile(np.array([0.42, 0.45, 0.52], np.float32), (length + 1, 1))
        # Later, more specific domains win: sort by span width, widest first.
        ordered = sorted(
            (d for d in self.domains
             if d[self.species]["start"] and d[self.species]["end"]),
            key=lambda d: d[self.species]["end"] - d[self.species]["start"],
            reverse=True,
        )
        for d in ordered:
            s, e = d[self.species]["start"], d[self.species]["end"]
            lo, hi = max(1, int(s)), min(length, int(e))
            if lo <= hi:
                table[lo:hi + 1] = hex_to_rgb(d["color"])
        return table

    def domain_of(self, residue: int) -> dict | None:
        best = None
        for d in self.domains:
            span = d[self.species]
            if span["start"] is None or span["end"] is None:
                continue
            if span["start"] <= residue <= span["end"]:
                width = span["end"] - span["start"]
                if best is None or width < best[0]:
                    best = (width, d)
        return best[1] if best else None


def load_domain_palette(species: str = "human") -> DomainPalette:
    data = json.loads((RESOURCE_DIR / "domains.json").read_text())
    return DomainPalette(species=species, domains=data["domains"])


def domain_colors(structure: Structure, palette: DomainPalette,
                  max_residue: int = 3000) -> np.ndarray:
    table = palette.lookup_table(max_residue)
    idx = np.clip(structure.res_seq, 0, max_residue)
    out = table[idx]
    # Ligands and other heteroatoms keep element colouring.
    het = structure.hetero
    if het.any():
        out[het] = structure.element_colors()[het]
    return out.astype(np.float32)


# --------------------------------------------------------------------------
# Scalar mappings
# --------------------------------------------------------------------------

def _viridis(t: np.ndarray) -> np.ndarray:
    """Compact viridis approximation — perceptually uniform, colour-blind safe."""
    t = np.clip(t, 0.0, 1.0)[:, None]
    c0 = np.array([0.2777, 0.0055, 0.3341])
    c1 = np.array([0.1050, 1.4046, 1.3845])
    c2 = np.array([-0.3308, 0.2148, 0.0952])
    c3 = np.array([-4.6342, -5.7991, -19.3324])
    c4 = np.array([6.2282, 14.1799, 56.6905])
    c5 = np.array([4.7763, -13.7451, -65.3532])
    c6 = np.array([-5.4354, 4.6456, 26.3125])
    return np.clip(c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6))))),
                   0.0, 1.0)


def value_colors(values: np.ndarray, vmin: float | None = None,
                 vmax: float | None = None) -> np.ndarray:
    v = np.asarray(values, dtype=np.float64)
    lo = np.nanpercentile(v, 2) if vmin is None else vmin
    hi = np.nanpercentile(v, 98) if vmax is None else vmax
    if hi <= lo:
        hi = lo + 1.0
    return _viridis((v - lo) / (hi - lo)).astype(np.float32)


def bfactor_colors(structure: Structure) -> np.ndarray:
    return value_colors(structure.b_factor)


def plddt_colors(structure: Structure) -> np.ndarray:
    """AlphaFold confidence colouring; ``b_factor`` holds pLDDT in AFDB files."""
    out = np.zeros((structure.n_atoms, 3), np.float32)
    b = structure.b_factor
    for threshold, color in PLDDT_BANDS:
        out[b >= threshold] = color
    return out


def distinct_colors(n: int, saturation: float = 0.62,
                    value: float = 0.92) -> np.ndarray:
    """``n`` evenly spaced hues — for labelling arbitrary groups."""
    return np.array([colorsys.hsv_to_rgb(i / max(n, 1), saturation, value)
                     for i in range(n)], dtype=np.float32)
