#!/usr/bin/env python
"""Regenerate the replicable panels of Liu et al. 2025 from our own structures.

Each figure is stamped with the panel it replicates and, for an analogue, with
the caveat **burnt into the image** rather than left in a caption that can be
cropped off — the same rule ``make_guo2017_figures.py`` follows, and for the
same reason: these are pictures that will end up beside the published ones.

Usage::

    python scripts/make_liu2025_figures.py
    python scripts/make_liu2025_figures.py --only iv --outdir /tmp
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib                                              # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                # noqa: E402
import numpy as np                                             # noqa: E402

from piezo1.analysis.liu2025 import PAPER                      # noqa: E402
from piezo1.analysis.liu2025_panels import (STATES, cap_gate_loop_span,  # noqa: E402
                                            cavity_volumes, curvature_radius,
                                            load_state, pore_radius_profile,
                                            axis_length, v2476_diagonal)
from piezo1.analysis.liu2025_permeation import (PUBLISHED, VOLTAGES,  # noqa: E402
                                                sweep_voltages)

COLORS = {"PIEZO1-Curved": "#8a8a8a", "S2472E-Curved": "#5fbf6a",
          "S2472E-Intermediate": "#c0308a", "PIEZO1-Flattened": "#3b4d9e"}
CAVEAT_STYLE = dict(fontsize=7, color="#a33", style="italic", wrap=True)


def _stamp(fig, panel: str, caveat: str = "") -> None:
    fig.text(0.99, 0.005, f"Replicates {panel} — {PAPER['journal']}",
             fontsize=6, va="bottom", ha="right", color="#666")
    if caveat:
        fig.text(0.01, 0.015, "ANALOGUE: " + caveat, va="bottom",
                 **CAVEAT_STYLE)


def figure_iv(outdir: Path) -> Path:
    """Figure 5E: current-voltage, against their ~20 pS slope."""
    fig, ax = plt.subplots(figsize=(5.2, 4.0))
    for state in STATES:
        sweep = sweep_voltages(load_state(state), pathway="lateral")
        if not sweep.conducts:
            continue
        ax.plot(sweep.voltages, sweep.currents_pA, "o-", label=
                f"{state} ({sweep.slope_pS():.0f} pS)", color=COLORS[state])
    published = np.array(VOLTAGES) * PUBLISHED["slope_conductance_pS"] * 1e12 * 1e-12
    ax.plot(VOLTAGES, published, "k--", lw=1.2,
            label=f"Liu et al. MD ({PUBLISHED['slope_conductance_pS']:.0f} pS)")
    ax.axhline(0, color="#bbb", lw=0.6)
    ax.axvline(0, color="#bbb", lw=0.6)
    ax.set_xlabel("Transmembrane potential (V)")
    ax.set_ylabel("Current (pA)")
    ax.set_title("Figure 5E — I–V through the pore, lateral pathway")
    ax.legend(fontsize=7, frameon=False)
    _stamp(fig, "Figure 5E",
           "1-D continuum drift-diffusion, not molecular dynamics. No explicit "
           "ions, no water, and the lateral portal is excluded rather than "
           "modelled — so the current is an upper bound.")
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    return _save(fig, outdir / "liu2025_5e_iv.png")


def figure_pore(outdir: Path) -> Path:
    """Figure 2D: pore radius along the axis, four states."""
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    for state in STATES:
        data = pore_radius_profile(load_state(state))
        z = data["z"] - np.median(data["z"])
        ax.plot(data["radius"], z, lw=1.4, color=COLORS[state], label=state)
    ax.axvline(2.0, color="#999", ls=":", lw=1.0)
    ax.text(2.05, ax.get_ylim()[0] + 4, "their x-axis starts here",
            fontsize=7, color="#666")
    ax.set_xlabel("Pore radius (Å)")
    ax.set_ylabel("Position along the axis (Å)")
    ax.set_title("Figure 2D — pore radius profile")
    ax.legend(fontsize=7, frameon=False)
    _stamp(fig, "Figure 2D")
    fig.text(0.01, 0.015, "Every entry is pinched to about 1 Å at R2295, "
                          "below the published plot's own axis.",
             fontsize=7, color="#444", va="bottom")
    fig.tight_layout(rect=(0, 0.06, 1, 0.96))
    return _save(fig, outdir / "liu2025_2d_pore.png")


def figure_distances(outdir: Path) -> Path:
    """Figures 2B, 2E and 3F/H: measured against published, side by side."""
    rows = []
    for state in STATES:
        structure = load_state(state)
        loops = cap_gate_loop_span(structure)
        rows.append((state, axis_length(structure), v2476_diagonal(structure),
                     loops["A2328-P2382"], loops["D2326-E2383"]))

    labels = ["R2295–E2537\naxis (Å)", "V2476\ndiagonal (Å)",
              "A2328–P2382\n(Å)", "D2326–E2383\n(Å)"]
    published = [(110.0, 100.0), (7.0, 14.0), (4.3, 16.2), (4.8, 12.8)]

    fig, axes = plt.subplots(1, 4, figsize=(10.5, 3.6))
    for j, (ax, label, (pub_c, pub_i)) in enumerate(zip(axes, labels, published)):
        values = [r[j + 1] for r in rows]
        ax.bar(range(len(rows)), [v if v is not None else 0 for v in values],
               color=[COLORS[r[0]] for r in rows])
        ax.axhline(pub_c, color="#333", ls="--", lw=1.0)
        ax.axhline(pub_i, color="#c0308a", ls="--", lw=1.0)
        ax.set_xticks(range(len(rows)))
        ax.set_xticklabels([r[0].replace("-", "\n") for r in rows],
                           fontsize=5.5, rotation=30, ha="right")
        ax.set_title(label, fontsize=8)
        ax.tick_params(labelsize=7)
    axes[0].set_ylabel("Å")
    fig.suptitle("Figures 2B, 2E, 3F/3H — measured (bars) vs published "
                 "(dashed: curved black, intermediate magenta)", fontsize=9)
    _stamp(fig, "Figures 2B, 2E, 3F/3H")
    fig.tight_layout(rect=(0, 0.06, 1, 0.90))
    return _save(fig, outdir / "liu2025_distances.png")


def figure_volumes(outdir: Path) -> Path:
    """Figure 2G: cavity volumes."""
    names = ["CV", "EV", "MV", "IV"]
    fig, ax = plt.subplots(figsize=(5.6, 3.8))
    width = 0.2
    for i, state in enumerate(STATES):
        volumes = cavity_volumes(load_state(state))
        ax.bar(np.arange(len(names)) + i * width,
               [volumes.get(n, 0.0) for n in names], width,
               label=state, color=COLORS[state])
    ax.set_xticks(np.arange(len(names)) + 1.5 * width)
    ax.set_xticklabels(names)
    ax.set_ylabel("Volume (Å$^3$)")
    ax.set_title("Figure 2G — cavity volumes")
    ax.legend(fontsize=7, frameon=False)
    _stamp(fig, "Figure 2G",
           "a solid of revolution about the measured pore path, so circular by "
           "construction and an over-estimate wherever the lumen is not. The "
           "direction of each change reproduces; the values do not.")
    fig.tight_layout(rect=(0, 0.09, 1, 0.96))
    return _save(fig, outdir / "liu2025_2g_volumes.png")


def figure_curvature(outdir: Path) -> Path:
    """Figure 6: curvature radius — the panel that disagrees."""
    published = {"PIEZO1-Curved": 11.0, "S2472E-Curved": 14.0,
                 "S2472E-Intermediate": 32.0, "PIEZO1-Flattened": 117.0}
    measured = {s: curvature_radius(load_state(s)).get("radius_nm", np.nan)
                for s in STATES}

    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    x = np.arange(len(STATES))
    ax.bar(x - 0.2, [published[s] for s in STATES], 0.4, label="published",
           color="#bbb")
    ax.bar(x + 0.2, [measured[s] for s in STATES], 0.4, label="our sphere fit",
           color=[COLORS[s] for s in STATES])
    ax.set_yscale("log")
    ax.set_xticks(x)
    ax.set_xticklabels([s.replace("-", "\n") for s in STATES], fontsize=7)
    ax.set_ylabel("Mid-plane curvature radius (nm, log scale)")
    ax.set_title("Figure 6 — curvature radius: where our fitter saturates")
    ax.legend(fontsize=7, frameon=False)
    _stamp(fig, "Figure 6 (R)",
           "our sphere fit agrees at 9.7 nm on the curved state, where it was "
           "calibrated against Guo & MacKinnon's 10.2, and saturates on the "
           "flat ones: 18.4 nm against their 117. Fitting a sphere to a nearly "
           "flat surface is ill-conditioned. Reported, not adjusted.")
    fig.tight_layout(rect=(0, 0.1, 1, 0.96))
    return _save(fig, outdir / "liu2025_6_curvature.png")


FIGURES = {"iv": figure_iv, "pore": figure_pore, "distances": figure_distances,
           "volumes": figure_volumes, "curvature": figure_curvature}


def _save(fig, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=160)
    plt.close(fig)
    print(f"  wrote {path}")
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--outdir", type=Path, default=Path("docs/img"))
    ap.add_argument("--only", choices=sorted(FIGURES), action="append")
    args = ap.parse_args()

    wanted = args.only or sorted(FIGURES)
    failures = 0
    for name in wanted:
        try:
            FIGURES[name](args.outdir)
        except FileNotFoundError as exc:
            print(f"  skipped {name}: {exc}")
        except Exception as exc:                              # noqa: BLE001
            print(f"  FAILED {name}: {type(exc).__name__}: {exc}")
            failures += 1
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
