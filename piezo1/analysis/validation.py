"""Non-parametric statistics for the blind variant test.

The tests here are the ones named in ``docs/PREREGISTRATION.md`` and nothing
else. They are deliberately assumption-light: group sizes are small and
unequal, and a mechanical model produces heavy-tailed scores, so a t-test would
be the wrong instrument even before the multiple-comparisons problem.

Everything is implemented directly rather than pulled from a stats package so
that the exact convention — one-sided direction, tie handling, whether the
observed statistic is counted in the permutation null — is visible and testable
rather than inherited.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["PermutationResult", "EffectSize", "permutation_test",
           "cliffs_delta", "bootstrap_cliffs_delta", "auroc",
           "interpret_delta"]


@dataclass
class PermutationResult:
    observed: float
    p_value: float
    null_mean: float
    null_std: float
    n_permutations: int
    alternative: str
    n_a: int
    n_b: int
    meta: dict = field(default_factory=dict)

    @property
    def significant(self) -> bool:
        return self.p_value < 0.05


@dataclass
class EffectSize:
    delta: float
    ci_low: float
    ci_high: float
    interpretation: str
    n_bootstrap: int = 0

    @property
    def excludes_zero(self) -> bool:
        return (self.ci_low > 0) or (self.ci_high < 0)


def permutation_test(a: np.ndarray, b: np.ndarray, n_permutations: int | None = None,
                     alternative: str = "less", seed: int = 0
                     ) -> PermutationResult:
    """Difference in means of ``a`` minus ``b``, tested by label shuffling.

    ``alternative="less"`` tests whether ``mean(a) < mean(b)``, which is the
    pre-registered direction: gain-of-function variants should *soften* the
    gating coordinate and so sit lower in ΔΔG than loss-of-function ones.

    The observed statistic is included in the null count — the ``(r + 1) /
    (n + 1)`` convention — so a p-value can never be exactly zero, which it
    should not be from a finite number of shuffles.
    """
    if n_permutations is None:
        n_permutations = int(_P.value("stats.n_permutations"))
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or len(b) < 2:
        raise ValueError(f"need at least two values per group, got {len(a)}, {len(b)}")

    observed = float(a.mean() - b.mean())
    pooled = np.concatenate([a, b])
    n_a = len(a)
    rng = np.random.default_rng(seed)

    null = np.empty(n_permutations)
    for i in range(n_permutations):
        rng.shuffle(pooled)
        null[i] = pooled[:n_a].mean() - pooled[n_a:].mean()

    if alternative == "less":
        count = int((null <= observed).sum())
    elif alternative == "greater":
        count = int((null >= observed).sum())
    elif alternative == "two-sided":
        count = int((np.abs(null) >= abs(observed)).sum())
    else:
        raise ValueError(f"unknown alternative {alternative!r}")

    return PermutationResult(
        observed=observed, p_value=(count + 1) / (n_permutations + 1),
        null_mean=float(null.mean()), null_std=float(null.std()),
        n_permutations=n_permutations, alternative=alternative,
        n_a=n_a, n_b=len(b))


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    """Cliff's delta: P(a > b) − P(a < b), in [−1, 1].

    Non-parametric and robust to the outliers a mechanical model will produce,
    which is why it was pre-registered instead of Cohen's d. Ties contribute
    zero to both terms.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) == 0 or len(b) == 0:
        return float("nan")
    comparison = a[:, None] - b[None, :]
    greater = int((comparison > 0).sum())
    less = int((comparison < 0).sum())
    return (greater - less) / float(len(a) * len(b))


def bootstrap_cliffs_delta(a: np.ndarray, b: np.ndarray,
                           n_bootstrap: int | None = None, seed: int = 0,
                           alpha: float = 0.05) -> EffectSize:
    """Cliff's delta with a percentile bootstrap confidence interval."""
    if n_bootstrap is None:
        n_bootstrap = int(_P.value("stats.n_bootstrap"))
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    point = cliffs_delta(a, b)
    rng = np.random.default_rng(seed)
    draws = np.empty(n_bootstrap)
    for i in range(n_bootstrap):
        draws[i] = cliffs_delta(rng.choice(a, size=len(a), replace=True),
                                rng.choice(b, size=len(b), replace=True))
    lo, hi = np.percentile(draws, [100 * alpha / 2, 100 * (1 - alpha / 2)])
    return EffectSize(delta=point, ci_low=float(lo), ci_high=float(hi),
                      interpretation=interpret_delta(point),
                      n_bootstrap=n_bootstrap)


def interpret_delta(delta: float) -> str:
    """Romano et al.'s conventional thresholds for Cliff's delta."""
    d = abs(delta)
    if d < 0.147:
        return "negligible"
    if d < 0.33:
        return "small"
    if d < 0.474:
        return "medium"
    return "large"


def auroc(scores: np.ndarray, positive: np.ndarray) -> float:
    """Area under the ROC curve, by the rank (Mann–Whitney) identity.

    ``positive`` is a boolean mask. Returns the probability that a randomly
    chosen positive scores above a randomly chosen negative, with ties counted
    as one half.
    """
    scores = np.asarray(scores, dtype=float)
    positive = np.asarray(positive, dtype=bool)
    n_pos = int(positive.sum())
    n_neg = int((~positive).sum())
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    order = np.argsort(scores, kind="mergesort")
    ranks = np.empty(len(scores), dtype=float)
    ranks[order] = np.arange(1, len(scores) + 1)
    # Average ranks within ties so a constant predictor scores exactly 0.5.
    _, inverse, counts = np.unique(scores, return_inverse=True, return_counts=True)
    sums = np.bincount(inverse, weights=ranks)
    ranks = (sums / counts)[inverse]
    return float((ranks[positive].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))
