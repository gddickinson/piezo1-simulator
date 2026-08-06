"""The single source of truth for every number a calculation depends on.

Loaded from ``resources/parameters.json``, which is authored and validated by
``scripts/build_parameters.py`` behind a provenance gate: a parameter must cite
a key that exists in ``references.json``, or declare itself a method choice and
say why.

**Why a registry rather than literals in the modules.** A constant written into
a function default is invisible. You cannot list them, you cannot show them to a
user, and you cannot tell which ones came from a paper and which somebody
guessed once. Every value here carries its unit, its bounds, its provenance and
its default, so the whole parameter set can be inspected, edited and — crucially
— *reported alongside any result it produced*.

**Overrides are tracked, not silent.** Changing a value is allowed; changing one
without the change being visible in every number that follows is not. Anything
that records a result — :mod:`piezo1.analysis.report`, and the documentation
claims in :mod:`piezo1.analysis.claims` — asks :func:`overrides` and refuses to
present a number as reproducing the literature when it was computed with
something else. That is the whole reason this module exists rather than a bag of
module-level globals.

Reading a value::

    from piezo1.parameters import PARAMETERS as P
    kappa = P.value("membrane.kappa")

Modules take ``None`` defaults and resolve through here at call time, so a
change takes effect on the next call rather than at import.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass, field

from .config import RESOURCE_DIR

__all__ = ["Parameter", "ParameterRegistry", "PARAMETERS", "value",
           "set_value", "reset", "overrides", "resolve"]


@dataclass(frozen=True)
class Parameter:
    """One documented number."""

    key: str
    name: str
    default: float
    unit: str
    kind: str                 # physical | empirical | method | convention
    category: str
    citation: str
    description: str
    source_note: str = ""
    minimum: float | None = None
    maximum: float | None = None

    @property
    def cited(self) -> bool:
        """Whether this cites the literature rather than a method choice."""
        return self.citation not in {"derived", "measured_here", "method_choice",
                                     "convention", "unverified"}

    def clamp(self, candidate: float) -> float:
        if self.minimum is not None:
            candidate = max(candidate, self.minimum)
        if self.maximum is not None:
            candidate = min(candidate, self.maximum)
        return candidate

    def format(self, current: float) -> str:
        text = f"{self.name}: {current:g}"
        if self.unit:
            text += f" {self.unit}"
        if current != self.default:
            text += f" (default {self.default:g})"
        return text


@dataclass
class ParameterRegistry:
    """Defaults from the resource, plus whatever the session has overridden."""

    parameters: dict[str, Parameter] = field(default_factory=dict)
    _overrides: dict[str, float] = field(default_factory=dict)
    sentinels: dict[str, str] = field(default_factory=dict)

    # ------------------------------------------------------------- reading

    def __contains__(self, key: str) -> bool:
        return key in self.parameters

    def __len__(self) -> int:
        return len(self.parameters)

    def get(self, key: str) -> Parameter:
        try:
            return self.parameters[key]
        except KeyError:
            raise KeyError(
                f"unknown parameter {key!r}. Every number a calculation depends "
                f"on must be registered — add it to scripts/build_parameters.py "
                f"with a citation and rebuild.") from None

    def value(self, key: str) -> float:
        """Current value: the override if there is one, else the default."""
        parameter = self.get(key)
        return self._overrides.get(key, parameter.default)

    def default(self, key: str) -> float:
        return self.get(key).default

    def is_default(self, key: str) -> bool:
        return key not in self._overrides

    def categories(self) -> list[str]:
        return sorted({p.category for p in self.parameters.values()})

    def in_category(self, category: str) -> list[Parameter]:
        return sorted((p for p in self.parameters.values()
                       if p.category == category), key=lambda p: p.key)

    # ------------------------------------------------------------- writing

    def set_value(self, key: str, candidate: float) -> float:
        """Override a parameter. Values outside the declared bounds are clamped.

        Clamping rather than raising is deliberate: these are exposed in a UI
        with free-text entry, and a typo should not take the application down.
        The clamped value is returned so the caller can show what was applied.
        """
        parameter = self.get(key)
        applied = parameter.clamp(float(candidate))
        if applied == parameter.default:
            self._overrides.pop(key, None)
        else:
            self._overrides[key] = applied
        return applied

    def reset(self, key: str | None = None) -> None:
        if key is None:
            self._overrides.clear()
        else:
            self._overrides.pop(key, None)

    def overrides(self) -> dict[str, float]:
        """Every parameter currently differing from its documented default."""
        return dict(self._overrides)

    @property
    def modified(self) -> bool:
        return bool(self._overrides)

    def override_summary(self) -> str:
        if not self._overrides:
            return "all parameters at their documented defaults"
        parts = [f"{key}={v:g} (default {self.default(key):g})"
                 for key, v in sorted(self._overrides.items())]
        return f"{len(parts)} parameter(s) overridden: " + "; ".join(parts)

    # --------------------------------------------------------- persistence

    def as_dict(self, only_overrides: bool = True) -> dict:
        if only_overrides:
            return dict(self._overrides)
        return {key: self.value(key) for key in self.parameters}

    def apply(self, values: dict) -> list[str]:
        """Apply a set of overrides, returning any keys that were not known."""
        unknown = []
        for key, candidate in (values or {}).items():
            if key in self.parameters:
                self.set_value(key, candidate)
            else:
                unknown.append(key)
        return unknown

    def provenance_rows(self) -> list[dict]:
        """Every parameter with its current value and where it came from."""
        return [{"key": p.key, "name": p.name, "value": self.value(p.key),
                 "default": p.default, "unit": p.unit, "kind": p.kind,
                 "category": p.category, "citation": p.citation,
                 "source_note": p.source_note,
                 "overridden": not self.is_default(p.key)}
                for p in sorted(self.parameters.values(), key=lambda x: x.key)]


def _load() -> ParameterRegistry:
    path = RESOURCE_DIR / "parameters.json"
    if not path.exists():
        return ParameterRegistry()
    raw = json.loads(path.read_text())
    registry = ParameterRegistry(sentinels=raw.get("sentinels", {}))
    for entry in raw["parameters"]:
        registry.parameters[entry["key"]] = Parameter(
            key=entry["key"], name=entry["name"], default=entry["value"],
            unit=entry["unit"], kind=entry["kind"], category=entry["category"],
            citation=entry["citation"], description=entry["description"],
            source_note=entry.get("source_note", ""),
            minimum=entry.get("minimum"), maximum=entry.get("maximum"))

    # A headless override file, for reproducing someone else's parameter set.
    # Deliberately not the GUI's persisted settings: a value someone once typed
    # into a dialog must not silently change what a script computes.
    override_path = os.environ.get("PIEZO1_PARAMETERS")
    if override_path and os.path.exists(override_path):
        unknown = registry.apply(json.loads(open(override_path).read()))
        if unknown:
            print(f"warning: PIEZO1_PARAMETERS lists unknown keys: {unknown}")
    return registry


#: The process-wide registry.
PARAMETERS = _load()


def value(key: str) -> float:
    return PARAMETERS.value(key)


def set_value(key: str, candidate: float) -> float:
    return PARAMETERS.set_value(key, candidate)


def reset(key: str | None = None) -> None:
    PARAMETERS.reset(key)


def overrides() -> dict[str, float]:
    return PARAMETERS.overrides()


def resolve(candidate: float | None, key: str) -> float:
    """``candidate`` if the caller supplied one, else the registered value.

    The pattern every module uses. Functions take ``None`` rather than a
    literal default so the registry is consulted **at call time** — an override
    then takes effect on the next call instead of requiring a reimport.
    """
    return PARAMETERS.value(key) if candidate is None else candidate
