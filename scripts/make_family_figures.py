#!/usr/bin/env python
"""Figures for the imported PIEZO-family census and what this project says back.

Four panels, each carrying the thing its number would otherwise be read without.

``family_constraint_by_domain.png``
    Constraint domain by domain on *this project's* boundaries, with the cap
    marked as the exception. The point is the ordering: pore machinery high,
    blade units low.

``family_blade_gradient.png``
    Why the census's distal-versus-proximal result is band composition. The two
    bands, what each is made of, and the reversal on the units.

``family_splay.png``
    The control that reinterprets the census's structural finding: an AlphaFold
    monomer splays further from an experimental structure of its *own* protein
    than two experimental paralogues do from each other.

``family_mechanics.png``
    Every mechanical feature against constraint, raw and with burial held fixed,
    with burial's own correlation drawn as the line each has to be read against.

Usage::

    python scripts/make_family_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import PROJECT_ROOT, STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402

OUT = PROJECT_ROOT / "docs" / "img"
CORE = "#2c6fb5"
BLADE = "#c26a3a"
GREY = "#7a7a7a"


def _save(fig, name: str) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT / name, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote {name}")


def constraint_by_domain() -> None:
    from piezo1.analysis.family_constraint import domain_constraint

    rows = sorted((d for d in domain_constraint() if d.mean is not None),
                  key=lambda d: d.mean)
    names = [d.domain for d in rows]
    values = [d.mean for d in rows]
    colours = [CORE if d.category in ("pore", "gate", "lever") else
               ("#8e7cc3" if d.category == "cap" else BLADE) for d in rows]
    whole = rows[0].mean - rows[0].vs_whole

    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    ax.barh(names, values, color=colours)
    ax.axvline(whole, color="k", ls="--", lw=1,
               label=f"whole protein ({whole:.2f})")
    ax.set_xlabel("mean constraint (Jensen-Shannon divergence, 174 orthologues)")
    ax.set_xlim(0, 0.95)
    ax.set_title("What half a billion years protected, on this project's domains\n"
                 "pore machinery blue · cap purple · blade units orange",
                 fontsize=10)
    ax.legend(loc="lower right", fontsize=8)
    ax.text(0.01, -0.16,
            "Values are the piezo_genes census's; the partition is ours. The two put the "
            "anchor 141 residues apart.",
            transform=ax.transAxes, fontsize=7, color=GREY)
    _save(fig, "family_constraint_by_domain.png")


def blade_gradient() -> None:
    from piezo1.analysis.family_constraint import blade_gradient as gradient

    g = gradient()
    fig, (left, right) = plt.subplots(1, 2, figsize=(9.2, 4.2))

    left.bar(["distal\nband", "proximal\nband"],
             [g["band_distal"], g["band_proximal"]], color=[BLADE, "#d9a066"])
    left.bar(["distal\nband", "proximal\nband"],
             [g["linker_distal"], g["linker_proximal"]], color=GREY, alpha=0.55,
             label="linker within the band")
    left.set_ylabel("mean constraint")
    left.set_ylim(0, 0.85)
    left.set_title("the census's chain-cut bands", fontsize=10)
    left.legend(fontsize=7, loc="upper right")
    for i, frac in enumerate([g["linker_fraction_distal"],
                              g["linker_fraction_proximal"]]):
        left.text(i, 0.03, f"{frac:.0%} linker", ha="center", fontsize=8,
                  color="white", fontweight="bold")

    right.bar(["THU1-6\n(distal)", "THU7-9\n(proximal)"],
              [g["unit_distal"], g["unit_proximal"]], color=[BLADE, "#d9a066"])
    right.set_ylim(0, 0.85)
    right.set_title("the transmembrane units alone", fontsize=10)
    right.annotate("the ordering reverses", xy=(1, g["unit_proximal"]),
                   xytext=(0.15, 0.80), fontsize=8,
                   arrowprops=dict(arrowstyle="->", lw=1))

    fig.suptitle("The distal-versus-proximal blade gradient is band composition",
                 fontsize=11)
    fig.text(0.5, -0.03,
             "Linker scores 0.517 and 0.515 either side, so it is not more conserved at one end; "
             "the band with more of it scores lower.",
             ha="center", fontsize=7, color=GREY)
    _save(fig, "family_blade_gradient.png")


def splay() -> None:
    from piezo1.analysis.core_periphery import compare

    def entry(name):
        path = STRUCTURE_DIR / f"{name}.cif"
        return Structure.from_file(path) if path.exists() else None

    # Labels must be unique: matplotlib's categorical axis merges duplicates,
    # which silently drops a bar and would hide half the control.
    pairs = [
        ("7WLT", "7WLU", "7WLT -> 7WLU\nPIEZO1 curved to flattened", CORE),
        ("AF-E2JF22-F1-model_v6", "6B3R",
         "AF monomer -> 6B3R\nSAME PROTEIN", "#b23a48"),
        ("AF-E2JF22-F1-model_v6", "7WLT",
         "AF monomer -> 7WLT\nSAME PROTEIN", "#b23a48"),
        ("7WLT", "6B3R", "7WLT -> 6B3R\nPIEZO1 to PIEZO1", GREY),
        ("6KG7", "9VEE", "6KG7 -> 9VEE\nPIEZO2 mouse to human", GREY),
        ("7WLT", "6KG7", "7WLT -> 6KG7\nDIFFERENT PARALOGUES", "#3d8c5f"),
        ("7WLT", "9VEE", "7WLT -> 9VEE\nDIFFERENT PARALOGUES", "#3d8c5f"),
    ]
    labels, ratios, colours = [], [], []
    for a, b, label, colour in pairs:
        left, right = entry(a), entry(b)
        if left is None or right is None:
            continue
        result = compare(left, right, a, b)
        if not result or result.splay_ratio is None:
            continue
        labels.append(label)
        ratios.append(result.splay_ratio)
        colours.append(colour)

    order = np.argsort(ratios)
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    ax.barh([labels[i] for i in order], [ratios[i] for i in order],
            color=[colours[i] for i in order])
    ax.axvline(1.0, color="k", lw=1, ls=":")
    ax.set_xlabel("blade RMSD / core RMSD, after fitting on the pore module alone")
    ax.set_title("The blades 'splaying' is the prediction, not the paralogue",
                 fontsize=11)
    ax.tick_params(axis="y", labelsize=8)
    ax.text(0.01, -0.20,
            "A prediction splays 7-9x from an experiment of its own protein; two experimental "
            "paralogues splay about 1x.",
            transform=ax.transAxes, fontsize=7, color=GREY)
    _save(fig, "family_splay.png")


def mechanics() -> None:
    from piezo1.analysis.constraint_mechanics import couple

    path = STRUCTURE_DIR / "7WLT.cif"
    if not path.exists():
        print("  skipping mechanics: 7WLT is not downloaded")
        return
    result = couple(Structure.from_file(path), structure_id="7WLT")
    features = sorted(result.features, key=lambda f: abs(f.partial_spearman))
    names = [f.feature for f in features]
    raw = [abs(f.spearman) for f in features]
    partial = [abs(f.partial_spearman) for f in features]
    burial = max(abs(v) for v in result.burial_alone.values())

    y = np.arange(len(names))
    fig, ax = plt.subplots(figsize=(7.6, 4.6))
    ax.barh(y + 0.19, raw, height=0.36, color="#b8cbe0", label="raw |rho|")
    ax.barh(y - 0.19, partial, height=0.36, color=CORE,
            label="with burial held fixed")
    ax.axvline(burial, color="#b23a48", lw=1.4, ls="--",
               label=f"burial alone ({burial:.2f})")
    ax.set_yticks(y, names, fontsize=8)
    ax.set_xlabel("|Spearman| against evolutionary constraint")
    ax.set_title("Does mechanics explain what evolution protected?", fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    for i, f in enumerate(features):
        if f.survives_correction and f.survives_burial:
            ax.text(partial[i] + 0.012, i - 0.19, "*", fontsize=11,
                    va="center", color=CORE)
    ax.text(0.01, -0.20,
            "* survives the circular-shift null, the multiple-comparison correction and the "
            "burial control. Null is a shift, not a permutation.",
            transform=ax.transAxes, fontsize=7, color=GREY)
    _save(fig, "family_mechanics.png")


def main() -> int:
    print("family figures ->", OUT)
    constraint_by_domain()
    blade_gradient()
    splay()
    mechanics()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
