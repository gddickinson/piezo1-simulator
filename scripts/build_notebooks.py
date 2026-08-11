#!/usr/bin/env python
"""Author the example notebooks, behind a gate that runs every code cell.

The `notebooks/` directory was created on the first day of the project and
never filled. Git does not record empty directories, so it existed only on the
machine that made it — a fresh clone did not have it at all. What went in its
place was `docs/NOTEBOOK.md`, which documents the headless API in prose.

Prose is not a substitute, because prose cannot be executed. These notebooks
are, and this script is what keeps them honest:

* the cell content lives in `notebook_content*.py` as ordinary Python, so it
  can be reviewed and diffed instead of being buried in JSON;
* **every code cell is executed, in order, in one namespace, before anything is
  written.** A notebook that raises is not published;
* the notebooks carry `assert`s on the numbers they quote, so running one
  checks the science rather than just the syntax.

They ship **without stored outputs**. A committed output is a number nobody
verifies, which is the exact failure this project spends most of its machinery
avoiding — and it goes stale silently, because nothing recomputes it. Run them.

    python scripts/build_notebooks.py           # verify and write
    python scripts/build_notebooks.py --check   # verify only, write nothing
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from notebook_content import NOTEBOOKS as CORE  # noqa: E402
from notebook_content_analysis import NOTEBOOKS as ANALYSIS  # noqa: E402

OUT = ROOT / "notebooks"
ALL = {**CORE, **ANALYSIS}

#: nbformat 4.4, written directly. `nbformat` is a convenience, not a
#: requirement, and adding a dependency so a build script can emit JSON it
#: already knows how to emit would be a poor trade.
NBFORMAT = (4, 4)


def _cell(kind: str, source: str) -> dict:
    lines = source.strip("\n").splitlines(keepends=True)
    if kind == "markdown":
        return {"cell_type": "markdown", "metadata": {}, "source": lines}
    return {"cell_type": "code", "execution_count": None, "metadata": {},
            "outputs": [], "source": lines}


def build(notebook: dict) -> dict:
    return {
        "cells": [_cell(kind, source) for kind, source in notebook["cells"]],
        "metadata": {
            "kernelspec": {"display_name": "Python 3 (piezo1)",
                           "language": "python", "name": "python3"},
            "language_info": {"name": "python", "version": "3.11"},
            "piezo1": {"title": notebook["title"],
                       "needs_data": notebook.get("needs_data", True)},
        },
        "nbformat": NBFORMAT[0],
        "nbformat_minor": NBFORMAT[1],
    }


def run_cells(notebook: dict) -> tuple[bool, str, float]:
    """Execute every code cell in order, sharing one namespace.

    This is what a reader does when they open the notebook and press run, so it
    is what the gate has to reproduce. Running the cells in a single namespace
    also catches the commonest notebook defect: a cell that only works because
    of a name defined in a cell the author later moved or deleted.
    """
    namespace: dict = {"__name__": "__notebook__"}
    started = time.time()
    for index, (kind, source) in enumerate(notebook["cells"], start=1):
        if kind != "code":
            continue
        try:
            exec(compile(source, f"<cell {index}>", "exec"), namespace)
        except Exception as exc:
            return False, f"cell {index}: {type(exc).__name__}: {exc}", \
                time.time() - started
    return True, "", time.time() - started


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--check", action="store_true",
                        help="run the cells but write nothing")
    parser.add_argument("--only", choices=sorted(ALL))
    args = parser.parse_args()

    OUT.mkdir(parents=True, exist_ok=True)
    names = [args.only] if args.only else list(ALL)
    failures = []

    for name in names:
        notebook = ALL[name]
        ok, error, seconds = run_cells(notebook)
        mark = "ok  " if ok else "FAIL"
        print(f"  {mark} {name:34s} {seconds:6.1f}s  {error}")
        if not ok:
            failures.append(name)
            continue
        if not args.check:
            path = OUT / f"{name}.ipynb"
            path.write_text(json.dumps(build(notebook), indent=1) + "\n")

    if failures:
        print(f"\n{len(failures)} notebook(s) did not run: "
              f"{', '.join(failures)}\nNothing was written for them.")
        return 1
    print(f"\nall {len(names)} notebooks ran"
          + ("" if args.check else f" and were written to {OUT}"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
