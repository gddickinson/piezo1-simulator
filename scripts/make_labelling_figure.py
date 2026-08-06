#!/usr/bin/env python
"""The labelling figure: the p^3 curve and the 1:2:3-dye histogram.

Kept out of ``make_figures.py`` because that script renders the molecule through
OpenGL and this is a plot; mixing the two would drag matplotlib into the render
path. Writes ``docs/img/labelling.png``.

Usage::

    python scripts/make_labelling_figure.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402

from piezo1.analysis.labelling import (LabellingConditions,               # noqa: E402
                                       occupancy_distribution,
                                       population_summary,
                                       site_labelled_fraction,
                                       time_to_fraction)
from piezo1.config import PROJECT_ROOT   # noqa: E402

OUT = PROJECT_ROOT / "docs" / "img"
INK, MUTED = "#e8edf5", "#9aa3b2"
DYE_COLORS = ["#4a5162", "#c26a4a", "#c9a227", "#4a9e78"]


def _style(ax) -> None:
    ax.set_facecolor("#161a21")
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)
    ax.title.set_color(INK)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    standard = LabellingConditions()

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), facecolor="#0f1216")
    for ax in axes:
        _style(ax)

    # (a) per-site vs whole-channel, at the standard protocol.
    t = np.linspace(0.0, 600.0, 2000)
    result = population_summary(t, standard)
    ax = axes[0]
    ax.plot(t / 60, result.p_site, color="#5b8dd6", lw=2, label="per site  $p$")
    ax.plot(t / 60, result.fully_labelled, color="#4a9e78", lw=2,
            label="all three  $p^3$")
    ax.plot(t / 60, result.detectable, color=MUTED, lw=1.2, ls="--",
            label=r"$\geq$1 dye")
    t99 = time_to_fraction(0.99, standard)
    ax.axvline(t99 / 60, color="#c26a4a", lw=1, ls=":")
    ax.annotate(f"99% at {t99 / 60:.1f} min", (t99 / 60, 0.35),
                xytext=(6, 0), textcoords="offset points",
                color="#c26a4a", fontsize=8)
    ax.set_xlabel("time (min)")
    ax.set_ylabel("labelled fraction")
    ax.set_title("a  200 nM JF646, live cell", loc="left", fontsize=10)
    ax.set_ylim(0, 1.04)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower right")

    # (b) the amplification: a per-site shortfall is cubed.
    ax = axes[1]
    p = np.linspace(0, 1, 400)
    ax.plot(p, p, color=MUTED, lw=1, ls="--", label="per site")
    ax.plot(p, p ** 3, color="#4a9e78", lw=2, label="$p^3$")
    for value in (0.8, 0.9):
        ax.plot([value, value, 0], [0, value ** 3, value ** 3], color="#c26a4a",
                lw=0.9, ls=":")
        ax.annotate(f"{value:.1f} → {value ** 3:.2f}", (0.02, value ** 3 + 0.02),
                    color="#c26a4a", fontsize=8)
    ax.set_xlabel("per-site labelled fraction $p$")
    ax.set_ylabel("channels fully labelled")
    ax.set_title("b  every site must bind", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left")

    # (c) the mixture, and what it takes to produce one.
    ax = axes[2]
    cases = [
        ("200 nM\n30 min", occupancy_distribution(
            float(site_labelled_fraction(1800.0, standard)), 3)),
        ("200 nM\n1 min", occupancy_distribution(
            float(site_labelled_fraction(60.0, standard)), 3)),
        ("90% of tags\nreactive", occupancy_distribution(
            float(site_labelled_fraction(
                6 * 3600.0, LabellingConditions(active_fraction=0.9))), 3)),
        ("60% of tags\nreactive", occupancy_distribution(
            float(site_labelled_fraction(
                6 * 3600.0, LabellingConditions(active_fraction=0.6))), 3)),
    ]
    x = np.arange(len(cases))
    bottom = np.zeros(len(cases))
    for k in range(4):
        heights = np.array([c[1][k] for c in cases])
        ax.bar(x, heights, 0.62, bottom=bottom, color=DYE_COLORS[k],
               label=f"{k} dye" + ("s" if k != 1 else ""),
               edgecolor="#0f1216", lw=0.5)
        bottom += heights
    ax.set_xticks(x)
    ax.set_xticklabels([c[0] for c in cases], fontsize=8)
    ax.set_ylabel("fraction of channels")
    ax.set_title("c  a dye mixture needs unreactive tags", loc="left",
                 fontsize=10)
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK,
              loc="lower left", ncol=4, bbox_to_anchor=(0.0, -0.32))

    fig.text(0.005, 0.005,
             "kinetics imported unchanged from halotag_binding_sim; "
             "reproduced to machine precision. Site positions are a model.",
             color=MUTED, fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    path = OUT / "labelling.png"
    fig.savefig(path, dpi=160, facecolor="#0f1216")
    print(f"  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
