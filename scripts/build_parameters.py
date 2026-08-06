#!/usr/bin/env python
"""Author ``piezo1/resources/parameters.json`` behind a provenance gate.

Every number a calculation depends on is declared here with its unit, its
bounds, and **where it came from**. The gate is that ``citation`` must either
resolve to a key in ``references.json`` or be one of the explicit sentinels
below — and a sentinel obliges the entry to explain itself in ``source_note``.

This is deliberately the same shape as ``build_variants.py``: authored content,
validated on the way out, committed as a resource rather than buried in code.
Keeping it as data means the whole parameter set can be read, reviewed and
diffed without opening a single module.

Usage::

    python scripts/build_parameters.py
    python scripts/build_parameters.py --check     # validate without writing
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from piezo1.config import RESOURCE_DIR  # noqa: E402

#: Citations that are not literature references. Each obliges a ``source_note``.
SENTINELS = {
    "derived": "computed from other quantities in this project",
    "measured_here": "measured by this project from deposited coordinates",
    "method_choice": "an algorithmic choice, not a measured quantity",
    "convention": "a community convention rather than a measurement",
    "unverified": "value in use but its source could not be confirmed",
}

from parameter_table import P  # noqa: E402


def validate(entries: list[dict]) -> list[str]:
    """Every entry must be complete and every citation must resolve."""
    refs = {e["key"] for e in json.loads(
        (RESOURCE_DIR / "references.json").read_text())["references"]}
    problems = []
    seen = set()
    for entry in entries:
        key = entry.get("key", "<missing>")
        if key in seen:
            problems.append(f"{key}: duplicate")
        seen.add(key)
        for field in ("key", "name", "value", "unit", "kind", "category",
                      "citation", "description"):
            if field not in entry:
                problems.append(f"{key}: missing {field}")
        if entry.get("kind") not in ("physical", "empirical", "method",
                                     "convention"):
            problems.append(f"{key}: bad kind {entry.get('kind')!r}")
        citation = entry.get("citation", "")
        if citation in SENTINELS:
            if not entry.get("source_note"):
                problems.append(
                    f"{key}: citation '{citation}' requires a source_note "
                    f"({SENTINELS[citation]})")
        elif citation not in refs:
            problems.append(f"{key}: citation {citation!r} is not in "
                            f"references.json — add the reference first")
        lo, hi = entry.get("minimum"), entry.get("maximum")
        if lo is not None and hi is not None:
            if not (lo <= entry["value"] <= hi):
                problems.append(f"{key}: value {entry['value']} outside "
                                f"[{lo}, {hi}]")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--check", action="store_true",
                    help="validate without writing")
    args = ap.parse_args()

    problems = validate(P)
    if problems:
        print(f"{len(problems)} problem(s):")
        for problem in problems:
            print(f"  {problem}")
        return 1

    by_category: dict[str, int] = {}
    for entry in P:
        by_category[entry["category"]] = by_category.get(entry["category"], 0) + 1
    print(f"{len(P)} parameters validated across {len(by_category)} categories")
    for category, count in sorted(by_category.items()):
        print(f"  {category:24s} {count}")

    cited = sum(1 for e in P if e["citation"] not in SENTINELS)
    print(f"\n{cited} cite a published reference; {len(P) - cited} are "
          f"method choices or conventions, each with a stated reason")

    if args.check:
        return 0
    dest = RESOURCE_DIR / "parameters.json"
    dest.write_text(json.dumps(
        {"schema": 1, "sentinels": SENTINELS, "parameters": P}, indent=1))
    print(f"\nwrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
