"""Does changing a registered parameter change the answer?

:mod:`provenance_chain` establishes that some code *reads* a parameter. That is
a static fact and it is not the same as the parameter mattering. A key can be
read into a variable that is then shadowed, passed to a function that ignores
it, or resolved on a branch nothing takes — and the static check would call it
wired.

This module settles it the only way that is not an inference: **override the
parameter, recompute, and see whether the number moved.** Round 49 needed
exactly this to prove that ``pore.step`` was inert, because reading the source
had already been misleading once.

The measurement has two halves and both matter:

- ``moved`` — the value changed under the override, so the parameter reaches
  the computation.
- ``restored`` — resetting the registry returns the original value **exactly**.
  A parameter that changes the answer but does not restore it is worse than one
  that does nothing, because it leaves the process in a state where the
  documented number can no longer be reproduced.

A parameter that fails ``restored`` is a bug in the caching, not in the wiring,
and this is where such a fault would surface first.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["ParameterEffect", "measure_effect", "probe_effects"]


@dataclass
class ParameterEffect:
    """One parameter, one probe, and whether the override reached the number."""

    key: str
    baseline: float | None = None
    overridden: float | None = None
    restored_value: float | None = None
    override_used: float | None = None
    error: str = ""

    @property
    def moved(self) -> bool:
        if self.baseline is None or self.overridden is None:
            return False
        return abs(self.overridden - self.baseline) > 1e-9

    @property
    def restored(self) -> bool:
        if self.baseline is None or self.restored_value is None:
            return False
        return abs(self.restored_value - self.baseline) <= 1e-12

    @property
    def ok(self) -> bool:
        return not self.error and self.moved and self.restored

    def summary(self) -> str:
        if self.error:
            return f"{self.key}: could not probe — {self.error}"
        return (f"{self.key}: {self.baseline:.6g} -> {self.overridden:.6g} "
                f"at {self.override_used:g} "
                f"({'moves' if self.moved else 'NO EFFECT'}, "
                f"{'restores' if self.restored else 'DOES NOT RESTORE'})")


def measure_effect(key: str, probe, override: float) -> ParameterEffect:
    """Override ``key``, run ``probe``, and report whether the answer moved.

    ``probe`` must return a float and must not itself set parameters. The
    registry is restored even when the probe raises, because leaving an
    override in place would silently corrupt every later measurement in the
    process — including the ones this module makes about other parameters.
    """
    from ..parameters import PARAMETERS, reset, set_value

    result = ParameterEffect(key=key, override_used=override)
    had = PARAMETERS.overrides()
    try:
        result.baseline = float(probe())
        set_value(key, override)
        result.overridden = float(probe())
    except Exception as exc:
        result.error = f"{type(exc).__name__}: {exc}"
    finally:
        reset()
        for existing, value in had.items():
            set_value(existing, value)
    if not result.error:
        try:
            result.restored_value = float(probe())
        except Exception as exc:  # pragma: no cover — probe was fine above
            result.error = f"{type(exc).__name__}: {exc}"
    return result


def probe_effects(probes: dict) -> list:
    """Run every ``{key: (probe, override)}`` pair and return the results."""
    return [measure_effect(key, probe, override)
            for key, (probe, override) in sorted(probes.items())]
