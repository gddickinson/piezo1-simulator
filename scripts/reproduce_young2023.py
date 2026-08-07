#!/usr/bin/env python
"""Reproduce the Young et al. 2023 four-state tension response end to end.

Round 40's integration test. The rate constants are Young's, taken from the
registry; the solver, the tension protocol and the time-constant extraction are
this project's. The output is checked against **two other papers'** measurements,
which is what makes it a test rather than a restatement:

* half-activation tension against Lewis et al. 2015 (2.7 +/- 0.1 mN/m)
* inactivation time constant against Bae et al. 2013 (8.6 +/- 0.4 ms)

Agreement on one and disagreement on the other is the result, and both are
reported. Writes ``docs/img/young2023_response.png``.

Usage::

    python scripts/reproduce_young2023.py
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt        # noqa: E402
import numpy as np                     # noqa: E402

from piezo1.config import PROJECT_ROOT  # noqa: E402
from piezo1.parameters import PARAMETERS  # noqa: E402
from piezo1.physics.kinetics import GatingModel  # noqa: E402

OUT = PROJECT_ROOT / "docs" / "img"
INK, MUTED = "#e8edf5", "#9aa3b2"
STEPS = (1.0, 2.0, 3.0, 5.0, 8.0)


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
    model = GatingModel()
    measured_t50 = PARAMETERS.value("kinetics.t50_measured")
    measured_tau = PARAMETERS.value("kinetics.wt_tau_ms")

    print("Young et al. 2023 rates -> this project's solver")
    print(f"  provenance: {model.provenance}")
    print()
    print(f"  {'tension':>8s} {'peak Po':>9s} {'tau (ms)':>10s}")
    responses = {}
    for tension in STEPS:
        result = model.step(tension, duration=1.0, n_points=8000)
        responses[tension] = result
        print(f"  {tension:7.1f}  {float(np.max(result.occupancy[:, 1])):9.4f}"
              f" {result.inactivation_tau() * 1e3:10.2f}")

    t50 = model.half_activation()
    tau = responses[5.0].inactivation_tau() * 1e3
    print()
    print(f"  T50 (peak)      {t50:6.3f} mN/m   vs Lewis 2015 "
          f"{measured_t50} +/- 0.1     -> {abs(t50 - measured_t50) / measured_t50:.1%}")
    print(f"  tau_inact @5    {tau:6.2f} ms     vs Bae 2013   "
          f"{measured_tau} +/- 0.4     -> {tau / measured_tau:.1f}x")
    scale = model.calibrate_k2_for_tau(measured_tau / 1000.0, hi=40.0)
    print(f"  k2 scale needed to reach {measured_tau} ms: {scale:.2f}x "
          f"(k2 {PARAMETERS.value('kinetics.k2'):.0f} -> "
          f"{PARAMETERS.value('kinetics.k2') * scale:.0f} /s)")

    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.0), facecolor="#0f1216")
    for ax in axes:
        _style(ax)

    colours = plt.cm.viridis(np.linspace(0.15, 0.9, len(STEPS)))
    ax = axes[0]
    for colour, tension in zip(colours, STEPS):
        result = responses[tension]
        ax.plot(result.time * 1e3, result.occupancy[:, 1], color=colour, lw=1.6,
                label=f"{tension:.0f} mN/m")
    ax.set_xlim(0, 300)
    ax.set_xlabel("time (ms)")
    ax.set_ylabel("open probability")
    ax.set_title("a  tension steps, Young 2023 rates", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)

    ax = axes[1]
    tensions = np.linspace(0, 12, 120)
    peak = np.array([float(np.ravel(model.peak_open_probability(t))[0])
                     for t in tensions])
    ax.plot(tensions, peak / peak.max(), color="#5b8dd6", lw=2, label="model")
    ax.axvline(t50, color="#4a9e78", lw=1.2, ls="--",
               label=f"model T50 {t50:.2f}")
    ax.axvspan(measured_t50 - 0.1, measured_t50 + 0.1, color="#c26a4a",
               alpha=0.35, label=f"Lewis 2015 {measured_t50}±0.1")
    ax.set_xlabel("membrane tension (mN/m)")
    ax.set_ylabel("normalised peak $P_o$")
    ax.set_title("b  half-activation — AGREES", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower right")

    ax = axes[2]
    taus = [responses[t].inactivation_tau() * 1e3 for t in STEPS]
    ax.plot(STEPS, taus, "o-", color="#c26a4a", lw=2, label="model")
    ax.axhspan(measured_tau - 0.4, measured_tau + 0.4, color="#4a9e78",
               alpha=0.4, label=f"Bae 2013 {measured_tau}±0.4 ms")
    ax.set_yscale("log")
    ax.set_xlabel("membrane tension (mN/m)")
    ax.set_ylabel("inactivation $\\tau$ (ms)")
    ax.set_title("c  inactivation — DISAGREES ~8.5x", loc="left", fontsize=10)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK)

    fig.text(0.005, 0.005,
             "Rate constants from Young et al. 2023; solver, protocol and "
             "time-constant extraction from this project. Checked against two "
             "independent papers.", color=MUTED, fontsize=7)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    path = OUT / "young2023_response.png"
    fig.savefig(path, dpi=160, facecolor="#0f1216")
    print(f"\n  wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
