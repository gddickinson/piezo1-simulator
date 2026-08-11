"""The questions Figure 4—figure supplement 1 asks of the surface potential.

Split from :mod:`piezo1.physics.electrostatics` at the project's length limit
and along a real seam: that module is *how to compute a screened potential
anywhere*, and this one is a small number of specific questions about PIEZO1's
surface.

The important one is :func:`patch_interaction`. The paper's claim is not that
two patches are coloured differently — it is that they attract, and do so
between protomers: "these electrostatic interactions appear to stabilize the
trimeric assembly in its curved conformation". A surface colour cannot say
that, because the potential at a patch includes the patch's own charges and
every other charge in the protein.
"""

from __future__ import annotations

import numpy as np

from .electrostatics import (bjerrum_length, debye_length, formal_charges,
                             surface_potential)

__all__ = ["patch_potential", "patch_interaction", "compare_conventions"]


def patch_potential(structure, residues, chains=None,
                    **kw) -> dict:
    """Mean surface potential over a named set of residues.

    Built for Figure 4-figure supplement 1, which names two charged patches —
    the cap's E2257/E2258/D2264 and the loops' R1761/R1762/R1269 — and asks
    whether they are oppositely charged where they meet. Residue numbers are
    in the structure's own numbering and nothing here converts them, so a
    caller reading a human entry with mouse numbers gets a real but different
    patch; state which numbering the numbers are in.
    """
    mask = structure.mask_residues(residues, chains)
    mask &= structure.mask_protein() & (~structure.hetero)
    if not mask.any():
        raise ValueError(f"no atoms match residues {sorted(residues)}"
                         + (f" in chains {chains}" if chains else "")
                         + " — wrong numbering system?")
    result = surface_potential(structure, mask=mask, **kw)
    if not len(result.points):
        return {"residues": sorted(residues), "n_surface_points": 0,
                "mean_potential_kT_per_e": float("nan"),
                "note": "the patch has no solvent-accessible surface"}
    return {
        "residues": sorted(residues),
        "n_surface_points": int(len(result.points)),
        "mean_potential_kT_per_e": float(result.potential.mean()),
        "median_potential_kT_per_e": float(np.median(result.potential)),
        "min_potential_kT_per_e": float(result.potential.min()),
        "max_potential_kT_per_e": float(result.potential.max()),
        "debye_length_A": result.debye_length,
        "method": result.meta["method"],
    }


def patch_interaction(structure, residues_a, residues_b,
                      ionic_strength: float | None = None) -> dict:
    """Screened electrostatic interaction between two named residue sets.

    Figure 4-figure supplement 1's claim is not that two patches are coloured
    differently — it is that they *attract*: "these electrostatic interactions
    appear to stabilize the trimeric assembly in its curved conformation". A
    surface colour cannot say that, because the potential at a patch includes
    the patch's own charges and every other charge in the protein. This
    computes the pairwise term alone:

    .. math::  U = \\ell_B \\sum_{i \\in A} \\sum_{j \\in B}
                   z_i z_j \\frac{e^{-r_{ij}/\\lambda_D}}{r_{ij}}

    in k_BT. Negative is attractive. Same-chain and cross-chain contributions
    are reported separately, because the paper's specific claim is that the
    contacts are **domain-swapped** — so a favourable term that turned out to
    be within one protomer would not support it.

    Charges are the formal ones on the named residues only. This is an
    interaction energy in a uniform screened medium, not a binding free energy:
    it has no desolvation term, so its magnitude is an overestimate of what the
    contact is worth even though its sign is reliable.
    """
    all_charges = formal_charges(structure)
    debye = debye_length(ionic_strength)
    bjerrum = bjerrum_length()

    set_a, set_b = set(int(r) for r in residues_a), set(int(r) for r in residues_b)

    def pick(wanted: set) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        keep, chains = [], []
        for k, label in enumerate(all_charges.label):
            chain, rest = label.split("/", 1)
            number = int("".join(ch for ch in rest if ch.isdigit()) or -1)
            if number in wanted:
                keep.append(k)
                chains.append(chain)
        keep = np.array(keep, dtype=int)
        return keep, np.array(chains), all_charges.charge[keep] if len(keep) else np.zeros(0)

    idx_a, chain_a, q_a = pick(set_a)
    idx_b, chain_b, q_b = pick(set_b)
    if not len(idx_a) or not len(idx_b):
        raise ValueError(
            f"one of the patches has no ionisable residue in this structure "
            f"({len(idx_a)} and {len(idx_b)} found) — wrong numbering system?")

    xyz_a, xyz_b = all_charges.xyz[idx_a], all_charges.xyz[idx_b]
    r = np.maximum(np.linalg.norm(xyz_a[:, None, :] - xyz_b[None, :, :], axis=-1),
                   1e-6)
    screen = np.exp(-r / debye) if np.isfinite(debye) else np.ones_like(r)
    pair = bjerrum * np.outer(q_a, q_b) * screen / r
    # A charge paired with itself has r = 0 and is not an interaction. This
    # only arises when the two patches overlap, which is exactly the case a
    # self-repulsion control uses — and without this the control returned
    # +6.4e7 k_BT from the clamped 1e-6 separation rather than a number that
    # could be read.
    overlap = idx_a[:, None] == idx_b[None, :]
    pair = np.where(overlap, 0.0, pair)
    same = chain_a[:, None] == chain_b[None, :]

    return {
        "residues_a": sorted(set_a), "residues_b": sorted(set_b),
        "n_charges_a": int(len(idx_a)), "n_charges_b": int(len(idx_b)),
        "energy_kT": float(pair.sum()),
        "same_chain_kT": float(pair[same].sum()),
        "cross_chain_kT": float(pair[~same].sum()),
        "attractive": bool(pair.sum() < 0),
        "domain_swapped": bool(abs(pair[~same].sum()) > abs(pair[same].sum())),
        "closest_pair_A": float(r.min()),
        "debye_length_A": debye,
        "caveat": ("a screened-Coulomb interaction energy in a uniform "
                   "dielectric: the sign is reliable, the magnitude is an "
                   "overestimate because nothing here pays desolvation"),
    }


def compare_conventions(structure, mask: np.ndarray | None = None) -> list[dict]:
    """The same surface under the choices that are open to argument.

    Ionic strength, histidine protonation and whether the modelled chain ends
    carry terminal charges. None of them is settled by the paper, and the point
    of the table is to show which ones could change a reading of the picture
    and which cannot.
    """
    rows = []
    for label, kw in (
            ("Figure 4c: 150 mM, neutral His", {}),
            ("no salt (unscreened)", {"ionic_strength": 0.0}),
            ("500 mM", {"ionic_strength": 0.5}),
            ("His at +0.1 (pH 7.4)", {"histidine": 0.1}),
            ("His fully protonated", {"histidine": 1.0}),
            ("modelled chain ends charged", {"include_termini": True})):
        result = surface_potential(structure, mask=mask, **kw)
        rows.append({
            "convention": label,
            "n_charges": result.meta["n_charges"],
            "net_charge": result.meta["net_charge"],
            "debye_length_A": result.debye_length,
            "mean_kT_per_e": float(result.potential.mean()),
            "p5_kT_per_e": float(np.percentile(result.potential, 5)),
            "p95_kT_per_e": float(np.percentile(result.potential, 95)),
            "fraction_saturated": result.fraction_saturated(),
        })
    return rows
