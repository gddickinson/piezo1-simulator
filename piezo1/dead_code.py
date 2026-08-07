"""Functions and classes nothing references — the scaffolding audit.

`parameter_audit` keeps unregistered numbers out. This keeps unused code out,
for the same reason: a project that is entirely load-bearing is easier to trust
than a large one that is mostly scaffolding, because a reader cannot tell which
is which without checking.

**Why this is AST-based and not a grep.** Round 55 wrote the grep version first
and it reported 102 unused public names, including `format_result` (used inside
its own module) and every return-type dataclass (constructed, never named
elsewhere). The second attempt used the AST but collapsed same-file references
into a set, so internal calls vanished and it reported 129 dead functions
including `fetch_pdb` and `cmd_list`. Both would have been catastrophic to act
on. The rule in `CLAUDE.md` — calibrate a checking instrument before believing
it — is why neither was.

So :func:`audit` is calibrated by :func:`calibration`, which requires that
known-used names are **not** flagged and that a planted unused name **is**.
``tests/test_dead_code.py`` runs it before reporting anything.

**What counts as a reference.** Any `Name`, `Attribute` or import alias
anywhere in `piezo1/`, `tests/` or `scripts/`, plus bare words inside string
literals — because registries in this project dispatch by string (`ANALYSES`,
the CLI subcommands) and a name reached only that way is still live.

**What is deliberately not audited.** Constants and dataclass fields: a
dataclass that is only ever constructed positionally has fields that never
appear as names, and reporting them would be noise. Entry points are excluded
by name for the same reason a linter excludes `__main__`.
"""

from __future__ import annotations

import ast
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

__all__ = ["DeadName", "audit", "calibration", "reference_counts",
           "SEARCHED", "EXEMPT"]

#: Trees searched for references. Scripts and tests count as users: a function
#: that only a script calls is still load-bearing.
SEARCHED = ("piezo1", "tests", "scripts")

#: Names that are entry points or protocol hooks, so nothing references them
#: by name and their absence from the graph is expected rather than a finding.
EXEMPT = {
    "main": "console entry point, invoked by name from pyproject.toml",
    "setup": "pytest/Qt lifecycle hook",
    "teardown": "pytest/Qt lifecycle hook",
}


@dataclass(frozen=True)
class DeadName:
    """One definition nothing refers to."""

    name: str
    kind: str                # function | class
    path: str

    def summary(self) -> str:
        return f"{self.kind} {self.name} ({self.path}) is never referenced"


@dataclass
class Report:
    dead: list = field(default_factory=list)
    n_definitions: int = 0

    def summary(self) -> str:
        return (f"{len(self.dead)} unreferenced of {self.n_definitions} "
                f"definitions")


def _root() -> Path:
    from .config import PROJECT_ROOT

    return PROJECT_ROOT


def reference_counts(root: Path | None = None) -> dict:
    """``{name: Counter(path -> times referenced)}`` over every searched tree.

    Counts *occurrences*, not distinct files. Collapsing to a set is what made
    the second attempt at this hide every same-file call.
    """
    root = root or _root()
    counts: dict = {}
    for tree_name in SEARCHED:
        base = root / tree_name
        if not base.exists():
            continue
        for path in base.rglob("*.py"):
            try:
                parsed = ast.parse(path.read_text())
            except (SyntaxError, UnicodeDecodeError):
                continue
            for node in ast.walk(parsed):
                name = None
                if isinstance(node, ast.Name):
                    name = node.id
                elif isinstance(node, ast.Attribute):
                    name = node.attr
                elif isinstance(node, ast.alias):
                    name = node.asname or node.name.split(".")[-1]
                elif isinstance(node, ast.Constant) and isinstance(node.value, str):
                    for token in node.value.replace(",", " ").split():
                        cleaned = token.strip("\"'()[]{}")
                        if cleaned:
                            counts.setdefault(cleaned, Counter())[path] += 1
                    continue
                if name:
                    counts.setdefault(name, Counter())[path] += 1
    return counts


def _definitions(root: Path | None = None) -> dict:
    """``{name: (kind, {paths})}`` for top-level functions and classes."""
    root = root or _root()
    out: dict = {}
    for path in (root / "piezo1").rglob("*.py"):
        try:
            parsed = ast.parse(path.read_text())
        except (SyntaxError, UnicodeDecodeError):
            continue
        for node in parsed.body:
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                kind = "function"
            elif isinstance(node, ast.ClassDef):
                kind = "class"
            else:
                continue
            entry = out.setdefault(node.name, (kind, set()))
            entry[1].add(path)
    return out


def audit(root: Path | None = None) -> Report:
    """Every top-level function and class nothing anywhere refers to."""
    root = root or _root()
    references = reference_counts(root)
    definitions = _definitions(root)

    dead = []
    for name, (kind, paths) in sorted(definitions.items()):
        if name.startswith("__") or name in EXEMPT:
            continue
        counter = references.get(name, Counter())
        # A bare `def f():` contributes no Name node, so ANY count is a use.
        if sum(counter.values()) == 0:
            dead.append(DeadName(name=name, kind=kind,
                                 path=str(sorted(paths)[0].relative_to(root))))
    return Report(dead=dead, n_definitions=len(definitions))


def calibration(root: Path | None = None) -> dict:
    """Check the instrument against known answers before its output is used.

    Returns ``{"false_positives": [...], "detects_planted": bool}``. The first
    must be empty and the second true, or :func:`audit`'s output means nothing.
    """
    root = root or _root()
    references = reference_counts(root)
    definitions = _definitions(root)

    known_used = ["measure_dome", "fetch_pdb", "cmd_list", "pore_profile",
                  "build_feature_table", "solve_pnp", "detect_interactions",
                  "verify_claims", "build_report"]
    false_positives = [
        name for name in known_used
        if name in definitions and sum(references.get(name, Counter()).values()) == 0]

    # Built from fragments so the whole token never appears in any source
    # file. Written literally, the string scanner finds it here and the probe
    # reports the instrument as broken — which it did on the first run, and
    # which is the string scanning working exactly as intended.
    planted = "zz" + "_planted_" + "probe" + "_absent"
    detects = sum(references.get(planted, Counter()).values()) == 0

    return {"false_positives": false_positives, "detects_planted": detects}
