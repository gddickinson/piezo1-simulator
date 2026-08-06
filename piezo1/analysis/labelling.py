"""HaloTag labelling of the PIEZO1 trimer, imported from ``halotag_binding_sim``.

The kinetics are **not new here**. They come from a companion project which
worked out what concentration and incubation time actually label all three tags
on a channel, and this module exists to put those numbers on the real structure:
the three sites are the three C-terminal anchors from
:mod:`piezo1.structure.fusion`, not three abstract slots.

The model is three equations, in the source project's own notation:

.. math::

    E(t)      &= \\text{partition} \\cdot [L] \\cdot
                 \\left(t - \\frac{1 - e^{-k_\\text{perm} t}}{k_\\text{perm}}\\right) \\\\
    p(t)      &= a \\cdot \\left(1 - e^{-k_\\text{on} E(t)}\\right) \\\\
    P_k(t)    &= \\binom{3}{k} p^k (1-p)^{3-k}

``E`` is cumulative ligand exposure (M·s), ``p`` the per-site labelled fraction,
and the last line the occupancy of a channel carrying three independent sites.
Living, fixed and in-vitro conditions differ **only** through ``k_perm``,
``partition`` and the active fraction ``a``.

**Why p³ is the whole point.** Every site must bind for a channel to be fully
labelled, so a modest per-site shortfall is cubed: p = 0.9 leaves only 0.73 of
channels fully labelled, and the rest appear as a mixture of one-, two- and
three-dye puncta. That mixture is directly comparable to the multi-level
amplitude histograms reported for JF646-BAPTA puncta.

**This module must not disagree with the source.** Any divergence is an import
error, not a discovery, so :func:`compare_with_source` re-runs the original
functions when that project is on the path and reports the largest difference.
The stochastic sampler reproduces the source's draw order exactly — two
``rng.random((n, 3))`` calls, reactivity first — because a different order gives
a statistically identical but numerically different population, and "exactly"
was the criterion.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from math import comb

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["LabellingConditions", "LabellingResult", "site_labelled_fraction",
           "cumulative_exposure", "occupancy_distribution",
           "fully_labelled_fraction", "detectable_fraction",
           "population_summary", "simulate_population", "predicted_brightness",
           "time_to_fraction", "label_sites", "compare_with_source",
           "SOURCE_PROJECT"]

#: Where the kinetics came from. Recorded so a result can name its origin.
SOURCE_PROJECT = "halotag_binding_sim"


@dataclass
class LabellingConditions:
    """Everything that distinguishes one labelling protocol from another.

    Defaults resolve from the registry at construction, so an override in the
    parameters dialog is picked up by the next call rather than at import.
    """

    concentration: float = field(
        default_factory=lambda: _P.value("labelling.concentration"))
    k_on: float = field(default_factory=lambda: _P.value("labelling.k_on"))
    k_perm: float = field(
        default_factory=lambda: _P.value("labelling.k_perm_live"))
    partition: float = field(
        default_factory=lambda: _P.value("labelling.partition_live"))
    active_fraction: float = field(
        default_factory=lambda: _P.value("labelling.active_fraction"))
    n_sites: int = field(
        default_factory=lambda: int(_P.value("labelling.n_sites")))
    name: str = "live cell"

    @property
    def asymptote(self) -> float:
        """The fully-labelled fraction at infinite time.

        ``active_fraction ** n_sites``. Worth having as a property because no
        incubation time can beat it: if fixation leaves 90% of tags reactive,
        73% of channels is the ceiling and waiting longer does nothing.
        """
        return float(self.active_fraction ** self.n_sites)

    def summary(self) -> str:
        return (f"{self.name}: {self.concentration * 1e9:.0f} nM, "
                f"k_on {self.k_on:.2e} 1/(M s), "
                f"tau_perm {1.0 / self.k_perm:.0f} s, "
                f"ceiling {self.asymptote:.3f}")


def cumulative_exposure(t, conditions: LabellingConditions | None = None):
    """Cumulative intracellular ligand exposure ``E(t)``, in M·s.

    The single quantity the chemistry depends on, which is why live, fixed and
    in-vitro all collapse onto one curve once expressed through it.
    """
    conditions = conditions or LabellingConditions()
    t = np.asarray(t, dtype=float)
    if conditions.k_perm <= 0:
        raise ValueError("k_perm must be > 0")
    steady = conditions.partition * conditions.concentration
    # Written as (1 - exp(-kt))/k rather than expanded, so that large k*t
    # underflows cleanly to E = steady*t instead of losing precision.
    relaxation = (1.0 - np.exp(-conditions.k_perm * t)) / conditions.k_perm
    return steady * (t - relaxation)


def site_labelled_fraction(t, conditions: LabellingConditions | None = None):
    """Per-site labelled fraction ``p(t)``, in ``[0, active_fraction]``."""
    conditions = conditions or LabellingConditions()
    exposure = cumulative_exposure(t, conditions)
    return conditions.active_fraction * (1.0 - np.exp(-conditions.k_on * exposure))


def occupancy_distribution(p_site, n_sites: int | None = None):
    """Binomial occupancy over ``k = 0..n_sites``, last axis.

    The 1:2:3-dye mixture: what fraction of channels carry exactly ``k`` dyes.
    """
    n_sites = int(_P.value("labelling.n_sites")) if n_sites is None else n_sites
    p = np.asarray(p_site, dtype=float)
    ks = np.arange(n_sites + 1)
    coefficients = np.array([comb(n_sites, int(k)) for k in ks])
    p_expanded = p[..., None]
    return (coefficients * (p_expanded ** ks)
            * ((1.0 - p_expanded) ** (n_sites - ks)))


def fully_labelled_fraction(p_site, n_sites: int | None = None):
    """``p ** n_sites`` — every site bound."""
    n_sites = int(_P.value("labelling.n_sites")) if n_sites is None else n_sites
    return np.asarray(p_site, dtype=float) ** n_sites


def detectable_fraction(p_site, n_sites: int | None = None):
    """``1 - (1 - p) ** n_sites`` — at least one dye, i.e. visible at all."""
    n_sites = int(_P.value("labelling.n_sites")) if n_sites is None else n_sites
    return 1.0 - (1.0 - np.asarray(p_site, dtype=float)) ** n_sites


@dataclass
class LabellingResult:
    """Time-resolved labelling of a population of channels."""

    t: np.ndarray
    p_site: np.ndarray
    occupancy: np.ndarray            # (len(t), n_sites + 1)
    fully_labelled: np.ndarray
    detectable: np.ndarray
    mean_dyes: np.ndarray
    conditions: LabellingConditions
    meta: dict = field(default_factory=dict)

    def at(self, time: float) -> dict:
        """The state at one time, by nearest sample."""
        i = int(np.argmin(np.abs(self.t - time)))
        return {"t": float(self.t[i]), "p_site": float(self.p_site[i]),
                "occupancy": self.occupancy[i].tolist(),
                "fully_labelled": float(self.fully_labelled[i]),
                "detectable": float(self.detectable[i]),
                "mean_dyes": float(self.mean_dyes[i])}

    def summary(self) -> str:
        final = self.at(float(self.t[-1]))
        return (f"after {final['t'] / 60:.0f} min: per-site "
                f"{final['p_site']:.3f}, fully labelled "
                f"{final['fully_labelled']:.3f}, detectable "
                f"{final['detectable']:.3f}, "
                f"{final['mean_dyes']:.2f} dyes per channel")


def population_summary(t, conditions: LabellingConditions | None = None
                       ) -> LabellingResult:
    """Everything the analytical model says about a labelling time course."""
    conditions = conditions or LabellingConditions()
    t = np.asarray(t, dtype=float)
    p = site_labelled_fraction(t, conditions)
    n = conditions.n_sites
    return LabellingResult(
        t=t, p_site=p, occupancy=occupancy_distribution(p, n),
        fully_labelled=fully_labelled_fraction(p, n),
        detectable=detectable_fraction(p, n),
        mean_dyes=n * p, conditions=conditions,
        meta={"source": SOURCE_PROJECT, "asymptote": conditions.asymptote})


def time_to_fraction(target: float, conditions: LabellingConditions | None = None,
                     t_max: float = 6 * 3600.0) -> float:
    """Time until the fully-labelled fraction reaches ``target``.

    Returns ``inf`` when the target is above the asymptote — which is a real
    answer, not a failure: a condition leaving 10% of tags unreactive can never
    fully label 80% of channels however long it runs.
    """
    conditions = conditions or LabellingConditions()
    if conditions.asymptote < target:
        return float("inf")
    lo, hi = 0.0, float(t_max)
    if fully_labelled_fraction(site_labelled_fraction(hi, conditions),
                               conditions.n_sites) < target:
        return float("inf")
    for _ in range(200):                       # bisection; p(t) is monotonic
        mid = 0.5 * (lo + hi)
        value = fully_labelled_fraction(
            site_labelled_fraction(mid, conditions), conditions.n_sites)
        if value < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def simulate_population(snapshot_times, conditions: LabellingConditions | None = None,
                        n_channels: int = 5000, seed: int = 0) -> dict:
    """Sample a finite population of channels, dye by dye.

    The sites are independent and never unbind, so a Gillespie walk would be
    wasted work: a site is labelled by time ``t`` with probability ``p(t)``, and
    drawing one uniform per site and comparing it against ``p(t)`` is the exact
    equivalent. Drawing the uniform *once* and reusing it across snapshots is
    what makes a single channel's history monotonic — resampling each snapshot
    would let a labelled site become unlabelled again.

    The two draws and their order are the source project's, so the populations
    match number for number rather than merely in distribution.
    """
    conditions = conditions or LabellingConditions()
    n_sites = conditions.n_sites
    snapshot_times = np.asarray(snapshot_times, dtype=float)
    rng = np.random.default_rng(seed)

    reactive = rng.random((n_channels, n_sites)) < conditions.active_fraction
    thresholds = rng.random((n_channels, n_sites))

    dye_counts = np.zeros((len(snapshot_times), n_channels), dtype=int)
    hist = np.zeros((len(snapshot_times), n_sites + 1))
    fully = np.zeros(len(snapshot_times))
    mean_dyes = np.zeros(len(snapshot_times))

    # p already includes active_fraction, but reactivity is drawn separately,
    # so divide it back out to avoid applying the same cap twice.
    ceiling = conditions.active_fraction if conditions.active_fraction > 0 else 1.0
    for i, t in enumerate(snapshot_times):
        p_reactive = float(np.clip(
            site_labelled_fraction(t, conditions) / ceiling, 0.0, 1.0))
        counts = (reactive & (thresholds < p_reactive)).sum(axis=1)
        dye_counts[i] = counts
        for k in range(n_sites + 1):
            hist[i, k] = float(np.mean(counts == k))
        fully[i] = float(np.mean(counts == n_sites))
        mean_dyes[i] = float(counts.mean())

    return {"t": snapshot_times, "dye_counts": dye_counts, "hist": hist,
            "fully_labeled": fully, "mean_dyes": mean_dyes,
            "n_sites": n_sites, "n_channels": n_channels}


def predicted_brightness(dye_counts, per_dye_intensity: float = 1.0,
                         background: float = 0.0, noise_cv: float | None = None,
                         seed: int = 1):
    """Turn integer dye counts into a predicted puncta-amplitude histogram.

    ``brightness = background + k * I * (1 + N(0, cv))``. The point of the
    exercise is whether the 1-, 2- and 3-dye levels stay separable once
    photophysical spread is included, so the CV is a registered parameter
    rather than a literal.
    """
    noise_cv = (_P.value("labelling.brightness_noise_cv")
                if noise_cv is None else noise_cv)
    rng = np.random.default_rng(seed)
    counts = np.asarray(dye_counts, dtype=float)
    signal = counts * per_dye_intensity
    return background + signal + rng.normal(0.0, noise_cv, size=counts.shape) * signal


def label_sites(fusion_model, t: float | None = None,
                conditions: LabellingConditions | None = None,
                n_channels: int = 5000, seed: int = 0) -> dict:
    """Put the labelling statistics onto the real tag positions.

    This is what the structure adds over the source project: the three sites are
    the three modelled tag centres from :mod:`piezo1.structure.fusion`, so an
    occupancy draw can be shown per site in space rather than as a bare count.

    The sites are treated as **equivalent**, because they are: the trimer is
    C3-symmetric, all three C-termini sit at the same height and radius, and the
    ligand reaches them from the same cytosolic pool. Position therefore changes
    where a dye is drawn, not how likely it is to be there.
    """
    conditions = conditions or LabellingConditions()
    t = _P.value("labelling.incubation_time") if t is None else float(t)

    p = float(site_labelled_fraction(t, conditions))
    population = simulate_population([t], conditions, n_channels=n_channels,
                                     seed=seed)
    centres = np.asarray(fusion_model.tag_centres, dtype=float)
    n_sites = min(conditions.n_sites, len(centres))

    rng = np.random.default_rng(seed)
    occupied = rng.random(n_sites) < p
    return {"t": t, "p_site": p,
            "tag_centres": centres[:n_sites],
            "anchor_residues": list(fusion_model.anchor_residues[:n_sites]),
            "occupied": occupied,
            "n_dyes": int(occupied.sum()),
            "occupancy": occupancy_distribution(p, conditions.n_sites).tolist(),
            "fully_labelled": float(fully_labelled_fraction(p, conditions.n_sites)),
            "hist": population["hist"][0].tolist(),
            "conditions": conditions.summary(),
            "note": ("site positions are a model (see structure.fusion); the "
                     "three sites are equivalent by C3 symmetry, so geometry "
                     "sets where a dye is drawn, not whether it binds")}


def compare_with_source(times=None, conditions: LabellingConditions | None = None
                        ) -> dict:
    """Re-run the original ``halotag_sim`` functions and report the difference.

    Returns ``{"available": False}`` when that project is not on the path, which
    is the normal case for a fresh clone — the equations are vendored here so
    this project stays self-contained. When it *is* present, every quantity must
    agree to machine precision: this is an import, and a divergence would mean
    the import is wrong, not that anything was discovered.
    """
    import importlib.util
    import sys
    from pathlib import Path

    root = Path(__file__).resolve().parents[3] / SOURCE_PROJECT
    if importlib.util.find_spec("halotag_sim") is None:
        if not (root / "halotag_sim" / "__init__.py").exists():
            return {"available": False, "reason": f"{SOURCE_PROJECT} not found"}
        sys.path.insert(0, str(root))

    try:
        from halotag_sim.kinetics import site_labeled_fraction as source_p
        from halotag_sim.parameters import CellCondition, KineticParams
        from halotag_sim.stochastic import simulate_population as source_pop
        from halotag_sim.trimer import (fully_labeled_fraction as source_full,
                                        occupancy_distribution as source_occ)
    except Exception as exc:                      # pragma: no cover
        return {"available": False, "reason": f"import failed: {exc}"}

    conditions = conditions or LabellingConditions()
    times = np.asarray(
        np.linspace(0.0, 2 * 3600.0, 241) if times is None else times, float)

    kp = KineticParams(k_on=conditions.k_on, n_sites=conditions.n_sites)
    cell = CellCondition(k_perm=conditions.k_perm, partition=conditions.partition,
                         active_fraction=conditions.active_fraction)

    theirs_p = np.asarray(source_p(times, conditions.concentration, kp, cell))
    ours_p = np.asarray(site_labelled_fraction(times, conditions))

    theirs_full = np.asarray(source_full(theirs_p, conditions.n_sites))
    ours_full = np.asarray(fully_labelled_fraction(ours_p, conditions.n_sites))

    theirs_occ = np.asarray(source_occ(theirs_p, conditions.n_sites))
    ours_occ = np.asarray(occupancy_distribution(ours_p, conditions.n_sites))

    snapshots = [0.0, 300.0, 900.0, 1800.0, 3600.0]
    theirs_pop = source_pop(conditions.concentration, kp, cell, snapshots,
                            n_channels=2000, seed=0)
    ours_pop = simulate_population(snapshots, conditions, n_channels=2000, seed=0)

    return {"available": True,
            "n_times": int(times.size),
            "max_abs_diff_p_site": float(np.max(np.abs(theirs_p - ours_p))),
            "max_abs_diff_fully_labelled": float(np.max(np.abs(theirs_full - ours_full))),
            "max_abs_diff_occupancy": float(np.max(np.abs(theirs_occ - ours_occ))),
            "dye_counts_identical": bool(np.array_equal(
                theirs_pop["dye_counts"], ours_pop["dye_counts"])),
            "max_abs_diff_hist": float(np.max(np.abs(
                theirs_pop["hist"] - ours_pop["hist"]))),
            "source": SOURCE_PROJECT}
