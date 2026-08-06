"""Intervals for the numbers this project reports.

Rounds 18 to 28 repeatedly found that a recorded number was stated with more
confidence than it had earned — a footprint area wrong by 3.5×, a null result
that could only ever have excluded a large effect, a T₅₀ eighteen times off. A
point estimate invites exactly that. This module attaches a spread.

**Three kinds of spread, kept apart because they mean different things.**

* :class:`Bootstrap` — resample the *data* (atoms in a fit, structures in an
  ensemble). This is a genuine confidence interval for sampling variability, and
  the only one of the three that deserves the name.
* :class:`Sensitivity` — vary a *method* choice, such as the elastic-network
  cutoff. There is no sampling distribution here; a cutoff is not a random
  variable. Reporting the spread as a confidence interval would be a second kind
  of overconfidence, so it is named differently and says so.
* :class:`ParameterRange` — vary a registered parameter across its published
  range, e.g. κ from 20 to 25 k_BT. Propagated uncertainty from an input, not a
  statement about this dataset.

**What none of them captures** is model error: whether a sphere is the right
shape for the dome, whether springs are the right physics. Bootstrapping a
sphere fit tells you how well the sphere is determined, not whether a sphere was
the right question. That limitation is stated on every result rather than left
for the reader to remember.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Callable, Sequence

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["Bootstrap", "Sensitivity", "ParameterRange", "bootstrap",
           "sensitivity", "parameter_range", "format_with_interval"]


@dataclass
class Bootstrap:
    """A confidence interval from resampling the data."""

    estimate: float
    low: float
    high: float
    n_resamples: int
    #: 1 - stats.alpha. Derived rather than repeated: two copies of a
    #: significance level drift, and this project already registers one.
    level: float = field(default_factory=lambda: 1.0 - _P.value("stats.alpha"))
    what: str = ""
    note: str = ""
    samples: np.ndarray = field(default_factory=lambda: np.zeros(0))

    kind = "confidence interval"

    @property
    def width(self) -> float:
        return self.high - self.low

    @property
    def relative_width(self) -> float:
        scale = max(abs(self.estimate), 1e-12)
        return self.width / scale

    def contains(self, value: float) -> bool:
        return self.low <= value <= self.high

    def summary(self, unit: str = "", digits: int = 3) -> str:
        return (f"{self.estimate:.{digits}f} "
                f"[{self.low:.{digits}f}, {self.high:.{digits}f}]"
                f"{' ' + unit if unit else ''} "
                f"({int(self.level * 100)}% CI, {self.n_resamples} resamples)")


@dataclass
class Sensitivity:
    """The spread of a result across a *method* choice.

    Deliberately not called an interval. A network cutoff has no sampling
    distribution, so quoting a percentile of it as a confidence interval would
    claim a statistical meaning it does not have.
    """

    estimate: float
    values: np.ndarray
    settings: Sequence
    knob: str = ""
    what: str = ""

    kind = "sensitivity range"

    @property
    def low(self) -> float:
        return float(np.min(self.values))

    @property
    def high(self) -> float:
        return float(np.max(self.values))

    @property
    def spread(self) -> float:
        return self.high - self.low

    @property
    def relative_spread(self) -> float:
        return self.spread / max(abs(self.estimate), 1e-12)

    def summary(self, unit: str = "", digits: int = 3) -> str:
        return (f"{self.estimate:.{digits}f} "
                f"(range {self.low:.{digits}f}–{self.high:.{digits}f}"
                f"{' ' + unit if unit else ''} over {self.knob} "
                f"{list(self.settings)}) — sensitivity, not a confidence "
                f"interval")


@dataclass
class ParameterRange(Sensitivity):
    """Spread from varying a registered parameter over its published range."""

    kind = "parameter range"

    def summary(self, unit: str = "", digits: int = 3) -> str:
        return (f"{self.estimate:.{digits}f} "
                f"(range {self.low:.{digits}f}–{self.high:.{digits}f}"
                f"{' ' + unit if unit else ''} over {self.knob} "
                f"{list(self.settings)}) — propagated from an input")


def bootstrap(statistic: Callable[[np.ndarray], float], data,
              n_resamples: int = 400, level: float | None = None,
              seed: int = 0, what: str = "", note: str = "") -> Bootstrap:
    """Percentile bootstrap of ``statistic`` over rows of ``data``.

    Resamples row indices with replacement, so ``data`` may be a coordinate
    array, a list of structures, or anything the statistic understands. A
    resample that fails is dropped and the count reported, rather than being
    silently replaced by the point estimate — which would narrow the interval.
    """
    level = (1.0 - _P.value("stats.alpha")) if level is None else level
    rng = np.random.default_rng(seed)
    n = len(data)
    if n < 3:
        raise ValueError(f"need at least three observations, got {n}")

    estimate = float(statistic(np.arange(n)))
    values = []
    for _ in range(n_resamples):
        index = rng.integers(0, n, n)
        try:
            value = float(statistic(index))
        except Exception:
            continue
        if np.isfinite(value):
            values.append(value)

    samples = np.asarray(values)
    if len(samples) < max(20, n_resamples // 10):
        raise RuntimeError(
            f"only {len(samples)} of {n_resamples} resamples succeeded; the "
            f"statistic is too fragile for a bootstrap to mean anything")
    alpha = (1.0 - level) / 2.0
    return Bootstrap(estimate=estimate,
                     low=float(np.quantile(samples, alpha)),
                     high=float(np.quantile(samples, 1.0 - alpha)),
                     n_resamples=len(samples), level=level, what=what,
                     note=note, samples=samples)


def sensitivity(statistic: Callable, settings: Sequence, reference=None,
                knob: str = "", what: str = "") -> Sensitivity:
    """Evaluate ``statistic`` at each method setting and report the spread."""
    values = []
    for setting in settings:
        try:
            values.append(float(statistic(setting)))
        except Exception:
            values.append(np.nan)
    array = np.asarray(values, dtype=float)
    good = np.isfinite(array)
    if not good.any():
        raise RuntimeError("no setting produced a usable value")
    estimate = (float(statistic(reference)) if reference is not None
                else float(np.median(array[good])))
    return Sensitivity(estimate=estimate, values=array[good],
                       settings=[s for s, ok in zip(settings, good) if ok],
                       knob=knob or "setting", what=what)


def parameter_range(statistic: Callable[[float], float], key: str,
                    values: Sequence[float], what: str = "") -> ParameterRange:
    """Vary a registered parameter over a stated range and report the spread.

    The parameter is restored afterwards, including if the statistic raises —
    leaving the registry modified would make every later number in the session
    incomparable with the documentation, which is exactly what the override
    tracking exists to prevent.
    """
    from ..parameters import PARAMETERS

    was_default = PARAMETERS.is_default(key)
    original = PARAMETERS.value(key)
    results = []
    try:
        estimate = float(statistic(original))
        for value in values:
            PARAMETERS.set_value(key, value)
            try:
                results.append(float(statistic(value)))
            except Exception:
                results.append(np.nan)
    finally:
        if was_default:
            PARAMETERS.reset(key)
        else:
            PARAMETERS.set_value(key, original)

    array = np.asarray(results, dtype=float)
    good = np.isfinite(array)
    return ParameterRange(estimate=estimate, values=array[good],
                          settings=[v for v, ok in zip(values, good) if ok],
                          knob=key, what=what)


def format_with_interval(value: float, spread, unit: str = "",
                         digits: int = 2) -> str:
    """One-line rendering that names which kind of spread it is."""
    if spread is None:
        return f"{value:.{digits}f}{' ' + unit if unit else ''} (no interval)"
    return spread.summary(unit=unit, digits=digits)
