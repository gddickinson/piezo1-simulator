"""Study design: could the test have worked, and did we test too many things?

:mod:`piezo1.analysis.validation` answers "did the predictor separate the
groups?". This module answers the two questions that have to be asked *around*
that one for the answer to mean anything:

* **Power.** A null result from an underpowered design rules out nothing. With
  16 gain-of-function and 9 loss-of-function variants — the Round 7 design —
  only a large effect is detectable at all, so "no separation" has to be read
  as "no *large* separation". Stating that is the difference between an honest
  null and an overclaimed one.
* **Multiplicity.** Round 22 will have several candidate predictors
  (AlphaMissense, EVE, ESM-1b, FoldX, mechanical coupling, conservation, and
  combinations). Testing all of them and reporting the best is how a null
  becomes a false positive. The primary endpoint must be named in advance and
  everything else corrected.

Plus **leave-one-out cross-validation**, so that a score which happens to fit
the 25 variants we have is distinguished from one that would generalise.

Nothing here touches the variant labels. It operates on group sizes and on
supplied arrays, so it can be — and is — run *before* an analysis rather than
after a disappointing one.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.special import ndtr, ndtri

from .validation import auroc, cliffs_delta, permutation_test

__all__ = ["PowerResult", "power_curve", "minimum_detectable_effect",
           "sample_size_for", "shift_for_delta", "delta_for_shift",
           "benjamini_hochberg", "MultipleComparisons", "leave_one_out",
           "LeaveOneOutResult"]


# --------------------------------------------------------------------------
# Effect size conversions
# --------------------------------------------------------------------------

def delta_for_shift(shift: float) -> float:
    """Cliff's delta produced by a location shift between two normals.

    For ``a ~ N(0,1)`` and ``b ~ N(shift,1)``, ``P(a > b) = Φ(−shift/√2)``, so
    ``δ = 2Φ(−shift/√2) − 1``. Used to quote power against the same effect-size
    scale the pre-registration commits to, rather than an abstract shift.
    """
    return float(2.0 * ndtr(-shift / np.sqrt(2.0)) - 1.0)


def shift_for_delta(delta: float) -> float:
    """Inverse of :func:`delta_for_shift`."""
    if not -1.0 < delta < 1.0:
        raise ValueError("Cliff's delta must lie strictly inside (-1, 1)")
    return float(-np.sqrt(2.0) * ndtri((delta + 1.0) / 2.0))


# --------------------------------------------------------------------------
# Power
# --------------------------------------------------------------------------

@dataclass
class PowerResult:
    """Simulated power of the pre-registered test at a range of effect sizes."""

    deltas: np.ndarray            # Cliff's delta (negative = predicted direction)
    power: np.ndarray             # probability of rejecting H0
    n_a: int
    n_b: int
    alpha: float
    n_simulations: int
    meta: dict = field(default_factory=dict)

    def power_at(self, delta: float) -> float:
        return float(np.interp(delta, self.deltas, self.power))

    def detectable(self, target_power: float = 0.8) -> float:
        """Smallest |δ| reaching ``target_power``; NaN if never reached."""
        order = np.argsort(np.abs(self.deltas))
        d, p = np.abs(self.deltas)[order], self.power[order]
        reached = np.flatnonzero(p >= target_power)
        if reached.size == 0:
            return float("nan")
        i = int(reached[0])
        if i == 0:
            return float(d[0])
        return float(np.interp(target_power, [p[i - 1], p[i]], [d[i - 1], d[i]]))

    def summary(self) -> str:
        mde = self.detectable()
        text = (f"n = {self.n_a} vs {self.n_b}, alpha = {self.alpha}: ")
        if np.isnan(mde):
            return text + "80% power is not reached at any effect size tested"
        return text + f"80% power requires |Cliff's delta| >= {mde:.2f}"


def _permutation_p_values(data: np.ndarray, n_a: int, n_permutations: int,
                          rng: np.random.Generator,
                          alternative: str = "less") -> np.ndarray:
    """One-sided permutation p-values for many datasets at once.

    The statistic is ``mean(a) − mean(b)``. Writing ``S`` for the total and
    ``S_a`` for the first group's sum, it equals
    ``S_a (1/n_a + 1/n_b) − S/n_b``, which is **monotonically increasing in
    S_a alone**. So the permutation distribution of the statistic is the
    permutation distribution of a subset sum, and only subset sums need be
    computed. That identity is what makes a proper simulation-based power
    analysis cheap enough to run; there is a test asserting this fast path
    agrees with :func:`~piezo1.analysis.validation.permutation_test`.
    """
    n_sim, n = data.shape
    observed = data[:, :n_a].sum(axis=1)

    counts = np.zeros(n_sim, dtype=int)
    # Permutations are independent of the data, so one set can be reused across
    # simulations; each simulated p-value remains individually valid.
    for start in range(0, n_permutations, 256):
        block = min(256, n_permutations - start)
        idx = np.argsort(rng.random((block, n)), axis=1)[:, :n_a]
        sums = data[:, idx].sum(axis=2)              # (n_sim, block)
        if alternative == "less":
            counts += (sums <= observed[:, None]).sum(axis=1)
        else:
            counts += (sums >= observed[:, None]).sum(axis=1)
    return (counts + 1) / (n_permutations + 1)


def power_curve(n_a: int, n_b: int, deltas=None, n_simulations: int = 2000,
                n_permutations: int = 999, alpha: float = 0.05,
                alternative: str = "less", pool: np.ndarray | None = None,
                seed: int = 0) -> PowerResult:
    """Simulated power of the pre-registered permutation test.

    Parameters
    ----------
    n_a, n_b:
        Group sizes. For Round 7 these were 16 gain-of-function and 9
        loss-of-function.
    deltas:
        Cliff's delta values to evaluate. Negative is the predicted direction
        (group *a* lower than group *b*).
    pool:
        If given, simulate by **resampling this observed distribution** and
        shifting one group, instead of drawing from a normal. Real ΔΔG values
        are heavy tailed, and a difference-in-means test loses power badly on
        heavy tails, so the normal model would flatter the design.
    """
    deltas = np.asarray(deltas if deltas is not None
                        else np.linspace(-0.95, -0.05, 19), dtype=float)
    rng = np.random.default_rng(seed)
    n = n_a + n_b

    if pool is not None:
        pool = np.asarray(pool, dtype=float)
        scale = float(np.std(pool))
        draw = lambda size: rng.choice(pool, size=size, replace=True)  # noqa: E731
    else:
        scale = 1.0
        draw = lambda size: rng.standard_normal(size)                  # noqa: E731

    power = np.empty(len(deltas))
    achieved = np.empty(len(deltas))
    for i, delta in enumerate(deltas):
        shift = shift_for_delta(float(delta)) * scale
        data = draw((n_simulations, n))
        # `shift_for_delta` is defined with a ~ N(0,1) and b ~ N(shift,1), so
        # the displacement goes on group *b*. Adding it to group a instead
        # inverts the effect and reports the power of detecting the opposite
        # sign — caught by the achieved-delta diagnostic below, which is why
        # it is computed rather than assumed.
        data[:, n_a:] += shift
        p = _permutation_p_values(data, n_a, n_permutations, rng, alternative)
        power[i] = float((p < alpha).mean())
        achieved[i] = float(np.mean([
            cliffs_delta(row[:n_a], row[n_a:]) for row in data[:200]]))

    return PowerResult(
        deltas=deltas, power=power, n_a=n_a, n_b=n_b, alpha=alpha,
        n_simulations=n_simulations,
        meta={"n_permutations": n_permutations, "alternative": alternative,
              "model": "resampled from pool" if pool is not None else "normal",
              "achieved_delta": achieved,
              "statistic": "difference in means (as pre-registered)"})


def minimum_detectable_effect(n_a: int, n_b: int, target_power: float = 0.8,
                              **kw) -> float:
    """Smallest |Cliff's delta| the design can detect at ``target_power``."""
    return power_curve(n_a, n_b, **kw).detectable(target_power)


def sample_size_for(delta: float, target_power: float = 0.8,
                    ratio: float = 1.0, alpha: float = 0.05,
                    n_simulations: int = 1000, n_permutations: int = 499,
                    max_n: int = 400, seed: int = 0,
                    alternative: str = "less") -> int:
    """Group size ``n_a`` needed to detect an effect of size ``|delta|``.

    ``n_b = ratio·n_a``. Returns ``max_n`` if the target is not reached within
    the search range, which is a statement about the design rather than a
    failure.

    **Only the magnitude of ``delta`` is used.** A one-sided test is powered
    against an effect in the direction it tests, so the sign is set by
    ``alternative`` rather than by the caller. Taking the sign literally was a
    real bug: ``power_curve`` defaults to ``alternative="less"``, so a caller
    asking the natural question — "how many variants for a large effect?" with
    ``delta=0.43`` — injected the effect *against* the alternative, got roughly
    zero power at every size, and received ``max_n`` back. The docstring
    presented that as "not reachable", so the wrong answer looked like a
    finding about the design.
    """
    magnitude = abs(float(delta))
    signed = -magnitude if alternative == "less" else magnitude

    lo, hi = 4, max_n
    while lo < hi:
        mid = (lo + hi) // 2
        result = power_curve(mid, max(2, int(round(mid * ratio))),
                             deltas=[signed], n_simulations=n_simulations,
                             n_permutations=n_permutations, alpha=alpha,
                             alternative=alternative, seed=seed)
        if result.power[0] >= target_power:
            hi = mid
        else:
            lo = mid + 1
    return int(lo)


# --------------------------------------------------------------------------
# Multiplicity
# --------------------------------------------------------------------------

@dataclass
class MultipleComparisons:
    """Family of tests with false-discovery-rate control."""

    names: list[str]
    p_values: np.ndarray
    adjusted: np.ndarray
    rejected: np.ndarray
    alpha: float
    primary: str | None = None

    @property
    def n_significant(self) -> int:
        return int(self.rejected.sum())

    def table(self) -> list[dict]:
        return [{"name": n, "p": float(p), "q": float(q),
                 "significant": bool(r), "primary": n == self.primary}
                for n, p, q, r in zip(self.names, self.p_values,
                                      self.adjusted, self.rejected)]

    def summary(self) -> str:
        head = (f"{self.n_significant}/{len(self.names)} significant at "
                f"FDR {self.alpha}")
        if self.primary:
            i = self.names.index(self.primary)
            head += (f"; primary endpoint {self.primary} p={self.p_values[i]:.3f}"
                     f" (uncorrected, as pre-specified)")
        return head


def benjamini_hochberg(p_values, names=None, alpha: float = 0.05,
                       primary: str | None = None) -> MultipleComparisons:
    """Benjamini–Hochberg FDR control over a family of tests.

    Chosen over Bonferroni because the candidate predictors are correlated —
    AlphaMissense, EVE and ESM-1b all read evolutionary signal — and
    controlling the family-wise error rate across correlated tests is so
    conservative it would guarantee another null.

    A **primary** endpoint may be named. Its uncorrected p-value is what the
    pre-registered decision rule uses; the correction governs the secondary
    family. That is the standard arrangement and it must be fixed in advance,
    because choosing the primary after seeing the p-values is precisely the
    error the correction exists to prevent.
    """
    p = np.asarray(p_values, dtype=float)
    if np.any((p < 0) | (p > 1)):
        raise ValueError("p-values must lie in [0, 1]")
    names = list(names) if names is not None else [f"test_{i}" for i in range(len(p))]
    if primary is not None and primary not in names:
        raise ValueError(f"primary endpoint {primary!r} is not in the family")

    order = np.argsort(p)
    ranked = p[order]
    n = len(p)
    # Step-up: q_(i) = min over j >= i of  n/j * p_(j), enforced monotone.
    scaled = ranked * n / np.arange(1, n + 1)
    adjusted_sorted = np.minimum.accumulate(scaled[::-1])[::-1]
    adjusted = np.empty(n)
    adjusted[order] = np.minimum(adjusted_sorted, 1.0)

    return MultipleComparisons(
        names=names, p_values=p, adjusted=adjusted,
        rejected=adjusted < alpha, alpha=alpha, primary=primary)


# --------------------------------------------------------------------------
# Cross-validation
# --------------------------------------------------------------------------

@dataclass
class LeaveOneOutResult:
    """Out-of-sample performance of a score, one held-out variant at a time."""

    predictions: np.ndarray
    labels: np.ndarray
    auroc_out: float
    auroc_in: float
    n: int
    meta: dict = field(default_factory=dict)

    @property
    def optimism(self) -> float:
        """In-sample minus out-of-sample AUROC — how much fitting flattered it."""
        return self.auroc_in - self.auroc_out

    def summary(self) -> str:
        return (f"AUROC in-sample {self.auroc_in:.3f}, leave-one-out "
                f"{self.auroc_out:.3f} (optimism {self.optimism:+.3f}, n={self.n})")


def leave_one_out(scores, labels, combine=None) -> LeaveOneOutResult:
    """Leave-one-out cross-validation of a scoring rule.

    ``scores`` is ``(n_samples, n_features)`` or ``(n_samples,)``; ``labels`` is
    boolean, True for the positive class. ``combine`` maps a training set to a
    scoring function; the default standardises each feature on the training
    fold and sums, which is the simplest defensible combination and has no
    fitted parameters beyond location and scale.

    The point is honesty about a combined predictor. Fitting weights on 25
    variants and reporting the resulting AUROC would measure how well 25 points
    can be fitted, not whether the rule generalises. **Every** step that looks
    at the labels — including standardisation — must happen inside the fold.
    """
    scores = np.asarray(scores, dtype=float)
    if scores.ndim == 1:
        scores = scores[:, None]
    labels = np.asarray(labels, dtype=bool)
    if len(scores) != len(labels):
        raise ValueError("scores and labels differ in length")
    if labels.all() or not labels.any():
        raise ValueError("need both classes present")

    def default_combine(train_x, train_y):
        centre = np.nanmean(train_x, axis=0)
        spread = np.nanstd(train_x, axis=0)
        spread[spread < 1e-12] = 1.0
        del train_y                       # deliberately unsupervised
        return lambda x: np.nansum((x - centre) / spread, axis=-1)

    combine = combine or default_combine
    n = len(scores)
    out = np.empty(n)
    for i in range(n):
        keep = np.ones(n, dtype=bool)
        keep[i] = False
        score_fn = combine(scores[keep], labels[keep])
        out[i] = float(np.atleast_1d(score_fn(scores[i]))[0])

    # In-sample scores go through the same one-row-at-a-time path, so a
    # ``combine`` that is not vectorised cannot make the two arms disagree.
    full = combine(scores, labels)
    in_sample = np.array([float(np.atleast_1d(full(row))[0]) for row in scores])
    return LeaveOneOutResult(
        predictions=out, labels=labels,
        auroc_out=auroc(out, labels), auroc_in=auroc(in_sample, labels),
        n=n, meta={"n_features": scores.shape[1]})

