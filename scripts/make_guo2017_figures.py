#!/usr/bin/env python
"""Regenerate the replicable panels of Guo & MacKinnon 2017 as figures.

Writes into ``docs/img/guo2017/``. Every figure is stamped with the panel it
replicates, and every panel that is an **analogue** rather than a replication
carries its caveat burnt into the figure rather than left in a caption that can
be cropped off — a projection of an atomic model is not a 2D class average, and
a screened-Coulomb surface is not APBS.

Panels needing OpenGL (the ribbon views, the drawn dome, the pore surface) are
not here: ``scripts/make_figures.py`` and ``scripts/make_model_figures.py``
render those, and this script is matplotlib only so it runs headless without a
GL context.

Usage::

    python scripts/make_guo2017_figures.py
    python scripts/make_guo2017_figures.py --structure 7WLT
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib                                          # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                            # noqa: E402
import numpy as np                                         # noqa: E402

from piezo1.analysis.guo2017 import PANELS, coverage  # noqa: E402
from piezo1.analysis.hydropathy import (hydropathy_profile,      # noqa: E402
                                        load_reference)
from piezo1.analysis.projection import project_views              # noqa: E402
from piezo1.analysis.topology import build_topology, unit_extent  # noqa: E402
from piezo1.config import PROJECT_ROOT, STRUCTURE_DIR             # noqa: E402
from piezo1.core.structure import Structure                       # noqa: E402
from piezo1.physics.dome_idealised import (flattening_series,     # noqa: E402
                                           guo2017_dome)
from piezo1.physics.dome import open_probability                  # noqa: E402

OUT = PROJECT_ROOT / "docs" / "img" / "guo2017"
BACKGROUND = "#14181d"
FOREGROUND = "#c8ccd4"


def _style() -> None:
    plt.rcParams.update({
        "figure.facecolor": BACKGROUND, "axes.facecolor": BACKGROUND,
        "savefig.facecolor": BACKGROUND, "text.color": FOREGROUND,
        "axes.labelcolor": FOREGROUND, "xtick.color": FOREGROUND,
        "ytick.color": FOREGROUND, "axes.edgecolor": "#4a525e",
        "figure.dpi": 130, "font.size": 9,
    })


def _stamp(fig, panel_key: str, caveat: str = "") -> None:
    panel = next((p for p in PANELS if p.key == panel_key), None)
    label = panel.label if panel else panel_key
    tag = "" if panel is None or panel.status == "replicated" else \
        f"  —  {panel.status.upper()}"
    fig.text(0.005, 0.985, f"replicates {label}{tag}", va="top", ha="left",
             fontsize=7, color="#8a919e")
    if caveat:
        fig.text(0.005, 0.005, caveat, va="bottom", ha="left", fontsize=6.5,
                 color="#e8a33d", wrap=True)


# --------------------------------------------------------------------------

def figure_7d() -> Path:
    """The four theoretical activation curves."""
    fig, ax = plt.subplots(figsize=(4.4, 3.2))
    tension = np.linspace(0.0, 2.5, 400)
    for delta_g, colour in ((20.0, "#e05252"), (40.0, "#5b8def")):
        for delta_area in (60.0, 20.0):
            ax.plot(tension, open_probability(tension, delta_area, delta_g),
                    color=colour, lw=1.6)
            t50 = delta_g / delta_area
            ax.annotate(f"{delta_area:.0f} nm²", xy=(t50, 0.5),
                        xytext=(t50 + 0.05, 0.62), fontsize=7, color=colour)
    ax.set_xlabel("Tension (k$_B$T/nm²)")
    ax.set_ylabel("P$_o$")
    ax.set_xlim(0, 2.5)
    ax.set_ylim(0, 1.02)
    ax.set_title("red: ΔG = 20 k$_B$T   blue: 40 k$_B$T", fontsize=8)
    _stamp(fig, "7d")
    fig.tight_layout()
    path = OUT / "figure_7d_activation.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def figure_7c() -> Path:
    """Projected area and bending energy against flattening — Figure 7c, computed."""
    dome = guo2017_dome()
    series = flattening_series(dome, n=120)
    angle = [p.polar_angle_deg for p in series]

    fig, (left, right) = plt.subplots(1, 2, figsize=(7.4, 3.1))
    left.plot(angle, [p.projected_area_nm2 for p in series], color="#5b8def")
    left.axhline(dome.area, ls="--", lw=1.0, color="#8a919e")
    left.text(angle[0], dome.area + 4, f"total membrane area {dome.area:.0f} nm²",
              fontsize=7, color="#8a919e")
    left.set_xlabel("polar half-angle (degrees) — flattening →")
    left.set_ylabel("projected area (nm²)")
    left.invert_xaxis()

    right.plot(angle, [p.delta_bending_kT for p in series], color="#e05252",
               label="ΔG$_{bend}$")
    right.plot(angle, [-0.35 * p.delta_projected_nm2 for p in series],
               color="#5fbf7f", label="−γΔA$_{proj}$ at 0.1× lytic")
    right.axhline(0.0, lw=0.8, color="#4a525e")
    right.set_xlabel("polar half-angle (degrees) — flattening →")
    right.set_ylabel("energy change (k$_B$T)")
    right.invert_xaxis()
    right.legend(fontsize=7, facecolor=BACKGROUND, edgecolor="#4a525e",
                 labelcolor=FOREGROUND)

    _stamp(fig, "7c",
           "Flattening at constant membrane area. Complete flattening releases "
           f"the whole {dome.bending_energy:.0f} k_BT of bending energy — with the "
           "paper's own 20-40 k_BT for (dG_prot + dG_bend), that puts dG_prot "
           "near +180 k_BT.")
    fig.tight_layout(rect=(0, 0.07, 1, 0.97))
    path = OUT / "figure_7c_flattening.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def figure_3_supplements(reference: str = "mouse") -> Path:
    """Hydropathy across the chain, in the paper's three blocks."""
    profile = hydropathy_profile(reference=reference)
    helices = sorted(load_reference(reference)["transmembrane"],
                     key=lambda t: t["start"])
    blocks = [(1, 900), (901, 1800), (1801, profile.sequence_length)]

    fig, axes = plt.subplots(len(blocks), 1, figsize=(7.6, 5.4))
    for ax, (lo, hi) in zip(axes, blocks):
        mask = (profile.position >= lo) & (profile.position <= hi)
        ax.axhline(0.0, lw=0.7, color="#4a525e")
        ax.fill_between(profile.position[mask], 0, profile.value[mask],
                        where=profile.value[mask] > 0, color="#5b8def",
                        alpha=0.65, lw=0)
        ax.fill_between(profile.position[mask], 0, profile.value[mask],
                        where=profile.value[mask] <= 0, color="#e05252",
                        alpha=0.5, lw=0)
        for helix in helices:
            if helix["end"] < lo or helix["start"] > hi:
                continue
            ax.axvspan(helix["start"], helix["end"], color="#f2c14e",
                       alpha=0.16, lw=0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(-3.4, 3.4)
        ax.set_ylabel("hydropathy")
    axes[-1].set_xlabel(f"residue ({reference} numbering)")
    axes[0].set_title(f"Kyte–Doolittle, window {profile.window}; "
                      f"shaded bands are the annotated transmembrane helices",
                      fontsize=8)
    _stamp(fig, "3-S1")
    fig.tight_layout()
    path = OUT / "figure_3_supplements_hydropathy.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def figure_6b(structure) -> Path:
    """The pore-radius profile, against the three published radii."""
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis
    from piezo1.analysis.guo2017_mechanism import PORE_CONSTRICTIONS

    blocks, _ = protomer_blocks(structure)
    profile = pore_profile(structure, detect_c3_axis(blocks))
    fig, ax = plt.subplots(figsize=(3.4, 5.0))
    ax.plot(profile.radius, profile.z, color="#5b8def", lw=1.4)
    ax.axvline(4.0, ls="--", lw=0.9, color="#8a919e")
    # In axes coordinates, not data: the y axis is inverted, so anchoring the
    # label to z.min() put it above the top spine and z.max() below the bottom.
    ax.text(4.15 / 8.0, 0.97, "4 Å: TEA-permeable", fontsize=6.5,
            color="#8a919e", rotation=90, va="top", ha="left",
            transform=ax.transAxes)
    for offset, (residue, published) in enumerate(
            sorted(PORE_CONSTRICTIONS.items())):
        touching = [s for s in profile.slices if residue in s.lining]
        if not touching:
            continue
        best = min(touching, key=lambda s: s.radius)
        ax.plot([published], [best.z], "o", ms=4, color="#e05252")
        ax.plot([best.radius], [best.z], "o", ms=4, color="#5fbf7f")
        # Stagger the labels: the three constrictions are within 8 A of each
        # other along the axis and their text overlapped.
        ax.annotate(str(residue), xy=(best.radius, best.z),
                    xytext=(best.radius + 0.6, best.z + (offset - 1) * 5.0),
                    fontsize=7, color=FOREGROUND, va="center",
                    arrowprops=dict(arrowstyle="-", lw=0.5, color="#4a525e"))
    ax.set_xlabel("pore radius (Å)")
    ax.set_ylabel("position along the pore axis (Å)")
    ax.set_xlim(0, 8)
    ax.invert_yaxis()
    ax.set_title(f"{structure.name}: green ours, red published (HOLE)",
                 fontsize=8)
    _stamp(fig, "6b",
           "Our profiler is an independent Apollonius implementation; the "
           "published radii came from HOLE. The ~0.6 A systematic offset is "
           "reported rather than absorbed.")
    fig.tight_layout(rect=(0, 0.05, 1, 0.97))
    path = OUT / "figure_6b_pore.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def figure_2ab(structure) -> Path:
    """Simulated projections, top and side — the class-average analogue."""
    views = project_views(structure)
    fig, axes = plt.subplots(1, 2, figsize=(6.6, 3.5))
    for ax, (name, projection) in zip(axes, views.items()):
        ax.imshow(projection.image.T, cmap="Greys_r", origin="lower")
        bar = projection.scale_bar_pixels(10.0)
        y = projection.image.shape[1] * 0.08
        x0 = projection.image.shape[0] * 0.08
        ax.plot([x0, x0 + bar], [y, y], color="white", lw=2.5)
        ax.text(x0, y * 1.6, "10 nm", color="white", fontsize=7)
        ax.set_title(f"{name} view", fontsize=8)
        ax.axis("off")
    _stamp(fig, "2ab",
           "SIMULATED PROJECTION OF AN ATOMIC MODEL, not a 2D class average. "
           "No CTF, no defocus, no solvent and no detergent micelle - and the "
           "published side view's envelope is substantially micelle.")
    # Leave room at the top: the panel titles and the "replicates" stamp
    # collided at 0.97, and the stamp is the part a reader must not miss.
    fig.tight_layout(rect=(0, 0.06, 1, 0.93))
    path = OUT / "figure_2ab_projection.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def figure_3a(structure, reference: str = "mouse") -> Path:
    """The topology diagram, with the 4-TM units boxed as in Figure 3b."""
    topology = build_topology(reference, structure)
    fig, ax = plt.subplots(figsize=(11.0, 2.9))
    half = topology.meta["membrane_half"]
    x0, x1 = topology.meta["x_range"]
    ax.axhspan(-half, half, color="#2a3038", zorder=0)

    for element in topology.elements:
        if element.kind == "loop":
            continue
        alpha = 0.28 if element.resolved is False else 1.0
        ax.add_patch(plt.Rectangle(
            (element.x0, element.y0), element.x1 - element.x0,
            element.y1 - element.y0, facecolor=element.color, alpha=alpha,
            edgecolor="none", zorder=2))
        if element.kind == "tm_helix":
            ax.text(0.5 * (element.x0 + element.x1), 0.0, str(element.helix),
                    ha="center", va="center", fontsize=5.5,
                    color="#0d1014" if element.resolved is not False else "#8a919e",
                    zorder=3)
        else:
            ax.text(0.5 * (element.x0 + element.x1),
                    0.5 * (element.y0 + element.y1), element.label,
                    ha="center", va="center", fontsize=6, color=FOREGROUND,
                    zorder=3)

    for unit, (lo, hi) in unit_extent(topology).items():
        ax.add_patch(plt.Rectangle(
            (lo - 0.12, -half - 0.35), hi - lo + 0.24, 2 * half + 0.7,
            fill=False, edgecolor="#ff5f56", lw=0.9, zorder=4))
        ax.text(lo, half + 0.45, f"THU{unit}", fontsize=6.5, color="#ff5f56")

    ax.set_xlim(x0 - 0.6, x1 + 0.6)
    ax.set_ylim(-3.0, half + 2.4)
    ax.axis("off")
    ax.text(x0, half + 0.05, "extracellular", fontsize=6.5, color="#8a919e",
            va="bottom")
    ax.text(x0, -half - 0.05, "cytoplasmic", fontsize=6.5, color="#8a919e",
            va="top")
    ax.set_title(topology.summary() + f"  |  numbering: {topology.numbering}",
                 fontsize=8)
    _stamp(fig, "3a",
           "Helices this entry does not model are drawn faint rather than "
           "dropped: dropping one would put TM13 where TM1 belongs and "
           "silently renumber every helix after it.")
    fig.tight_layout(rect=(0, 0.07, 1, 0.96))
    path = OUT / "figure_3a_topology.png"
    fig.savefig(path)
    plt.close(fig)
    return path


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--structure", default="6B3R",
                    help="entry to measure on (default 6B3R, the paper's own)")
    args = ap.parse_args()

    _style()
    OUT.mkdir(parents=True, exist_ok=True)
    report = coverage()
    print(report["summary"])

    written = [figure_7d(), figure_7c(), figure_3_supplements()]

    path = STRUCTURE_DIR / f"{args.structure}.cif"
    if path.exists():
        structure = Structure.from_file(path)
        from piezo1.io.registry import load_registry
        record = load_registry().get(args.structure)
        reference = record.numbering_species if record else "mouse"
        written += [figure_3a(structure, reference), figure_6b(structure),
                    figure_2ab(structure)]
    else:
        print(f"  {args.structure} not downloaded — skipping the panels that "
              f"need coordinates. Run: python -m piezo1.io.fetch")

    for item in written:
        print(f"  wrote {item.relative_to(PROJECT_ROOT)}")
    print(f"\n{len(written)} figures in {OUT.relative_to(PROJECT_ROOT)}")
    print("Panels needing OpenGL (ribbons, the drawn dome, the pore surface) "
          "come from scripts/make_figures.py and make_model_figures.py.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
