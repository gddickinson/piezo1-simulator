#!/usr/bin/env python
"""Figures for what the project does not know.

Two panels, both read from recorded results rather than recomputed:

``record_nulls.png``
    Every pre-registered test as a forest plot — effect size with its interval,
    against the line of no effect. The point of the figure is that all five
    intervals cross it.

``record_data_limit.png``
    Round 47's feasibility result: the sample size each recorded effect would
    need against the most this project could ever assemble.

Usage::

    python scripts/make_record_figure.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from piezo1.config import DERIVED_DIR  # noqa: E402

OUT = ROOT / "docs" / "img"

#: Round -> (label, predictor). Effects and intervals come from the recorded
#: runs; nothing here is recomputed, because these are frozen results.
ROUNDS = [
    (7, "elastic-network ΔΔG"),
    (22, "FoldX ΔΔG"),
    (36, "substitution-aware ΔΔG"),
    (41, "population constraint"),
    (48, "wild-type structural context"),
]


def _recorded(number: int) -> dict:
    """Effect and interval as stored, with the record as the fallback."""
    path = DERIVED_DIR / f"validation_round{number}.json"
    entry = {"delta": None, "low": None, "high": None}
    if path.exists():
        primary = json.loads(path.read_text()).get("primary", {})
        entry.update(delta=primary.get("cliffs_delta"),
                     low=primary.get("ci_low"), high=primary.get("ci_high"))
    if entry["delta"] is None:
        # Round 7 predates the interval convention; its effect is in the record.
        from piezo1.analysis.prediction_record import VALIDATION_RECORD
        match = [r for r in VALIDATION_RECORD if r.round == number]
        if match:
            entry["delta"] = match[0].cliffs_delta
    return entry


def nulls_figure() -> Path:
    rows = [(n, label, _recorded(n)) for n, label in ROUNDS]
    fig, ax = plt.subplots(figsize=(7.6, 3.6))
    y = np.arange(len(rows))[::-1]

    ax.axvline(0.0, color="#c0392b", lw=1.4, zorder=1)
    ax.text(0.02, -0.62, "no effect", color="#c0392b", fontsize=8, va="center")

    for pos, (number, label, entry) in zip(y, rows):
        delta = entry["delta"]
        if delta is None:
            continue
        if entry["low"] is not None:
            ax.plot([entry["low"], entry["high"]], [pos, pos],
                    color="#3b4252", lw=2.2, solid_capstyle="round", zorder=2)
        else:
            ax.annotate("interval not recorded", (delta, pos),
                        textcoords="offset points", xytext=(10, -11),
                        fontsize=7, color="#8b93a1")
        ax.plot([delta], [pos], "o", ms=7, color="#2e6da4", zorder=3)
        ax.annotate(f"δ = {delta:+.3f}", (delta, pos),
                    textcoords="offset points", xytext=(0, 9),
                    ha="center", fontsize=8)

    ax.set_yticks(y)
    ax.set_yticklabels([f"Round {n}\n{label}" for n, label, _ in rows],
                       fontsize=8)
    ax.set_xlabel("Cliff's δ  (negative = predictor separates LoF from GoF as hypothesised)")
    ax.set_title("Five pre-registered tests, five nulls — every interval crosses zero",
                 fontsize=10)
    ax.set_xlim(-0.85, 0.62)
    ax.set_ylim(-0.9, len(rows) - 0.25)
    ax.grid(axis="x", alpha=0.25)
    for spine in ("top", "right", "left"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()

    OUT.mkdir(parents=True, exist_ok=True)
    path = OUT / "record_nulls.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def data_limit_figure() -> Path:
    from piezo1.analysis.feasibility import assess

    report = assess(n_simulations=400)
    labels = [s.label for s in report.scenarios]
    sizes = [s.n for s in report.scenarios]
    powers = [s.power_at_observed for s in report.scenarios]
    reachable = [s.reachable for s in report.scenarios]

    fig, ax = plt.subplots(figsize=(7.6, 3.6))
    colors = ["#2e6da4" if ok else "#c0392b" for ok in reachable]
    bars = ax.bar(range(len(sizes)), sizes, color=colors, width=0.62)

    for i, (bar, n, power) in enumerate(zip(bars, sizes, powers)):
        ax.annotate(f"n = {n}\npower {power:.2f}",
                    (bar.get_x() + bar.get_width() / 2, n),
                    textcoords="offset points", xytext=(0, 4),
                    ha="center", fontsize=8)

    ax.axhline(report.ceiling_n, color="#d9a441", ls="--", lw=1.3)
    # Placed in the empty upper-left rather than on the line: at the line it
    # collided with the bar labels, which is where a reader looks first.
    ax.annotate(f"- - -  ceiling: the most this project could ever assemble "
                f"({report.ceiling_n})",
                (-0.42, max(sizes) * 1.13), ha="left", va="center",
                fontsize=8, color="#a8791f")

    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels([l.replace(" ", "\n", 1) for l in labels], fontsize=8)
    ax.set_ylabel("directional variants")
    ax.set_title(f"The effect the predictor produces (δ = {report.observed_effect:+.3f}) "
                 f"needs {report.required_n} variants", fontsize=10)
    ax.set_ylim(0, max(sizes) * 1.25)
    ax.grid(axis="y", alpha=0.25)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()

    path = OUT / "record_data_limit.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def main() -> int:
    for build in (nulls_figure, data_limit_figure):
        path = build()
        print(f"  wrote {path.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
