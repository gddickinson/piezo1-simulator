#!/usr/bin/env python
"""Run the project from a genuinely empty clone, and report what breaks.

Round 60 did this by hand and found three defects that are invisible on a
developer machine, because the failing path is never taken once the files are on
disk: eight tests that failed instead of skipping, an Ensembl download that had
been broken for months, and a documented ceiling that silently fell from 59 to
34 when the corpus was absent.

Doing it by hand once finds the bugs that exist today. This makes it one
command, so the next cache-shaped defect is found by running something rather
than by remembering to.

**The specific thing it checks.** On an empty clone every data-dependent test
must **skip**, not fail. `conftest.py` states that rule and Round 60 found eight
tests breaking it. A skip is fine and expected; a failure is a reproducibility
bug. So this exits non-zero on failures and says nothing about skips beyond
counting them.

Usage::

    python scripts/cold_clone_check.py              # empty clone only, no network
    python scripts/cold_clone_check.py --fetch      # also fetch and re-run
    python scripts/cold_clone_check.py --keep       # leave the clone for inspection
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import time
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

@dataclass
class Step:
    """One stage of the chain, with what it cost and what it produced."""

    name: str
    seconds: float = 0.0
    passed: int = 0
    failed: int = 0
    skipped: int = 0
    errors: int = 0
    detail: str = ""
    ok: bool = True

    def line(self) -> str:
        mark = "ok  " if self.ok else "FAIL"
        counts = ""
        if self.passed or self.failed or self.skipped:
            counts = (f"{self.passed} passed, {self.failed} failed, "
                      f"{self.skipped} skipped")
        return (f"  {mark} {self.name:34s} {self.seconds:6.1f}s  "
                f"{counts or self.detail}")


@dataclass
class ColdCloneReport:
    steps: list = field(default_factory=list)
    clone: str = ""

    @property
    def failures(self) -> list:
        return [s for s in self.steps if not s.ok]

    @property
    def total_seconds(self) -> float:
        return sum(s.seconds for s in self.steps)

    def summary(self) -> str:
        if self.failures:
            return (f"{len(self.failures)} step(s) failed in "
                    f"{self.total_seconds:.0f}s — "
                    + "; ".join(s.name for s in self.failures))
        return (f"all {len(self.steps)} steps clean in "
                f"{self.total_seconds:.0f}s")


def _counts(output: str) -> tuple[int, int, int, int]:
    """(passed, failed, skipped, errors) from pytest's summary line.

    Searched for anywhere in the output rather than taken from the last line:
    with ``-q`` and no failures the final line is a progress bar, so reading
    the tail found zero of everything and reported a clean run as broken.

    And pytest is invoked without a verbosity flag, because ``pytest.ini``
    already sets ``-q`` — passing another made it ``-qq``, which suppresses the
    summary line entirely and left the counts at zero however they were parsed.
    """
    summary = ""
    for line in reversed(output.strip().splitlines()):
        if re.search(r"\d+ (passed|failed|skipped|error)", line):
            summary = line
            break

    def find(word):
        match = re.search(rf"(\d+) {word}", summary)
        return int(match.group(1)) if match else 0

    return find("passed"), find("failed"), find("skipped"), find("error")


def _run(name: str, command: list, cwd: Path, expect_tests: bool = False) -> Step:
    started = time.time()
    result = subprocess.run(command, cwd=cwd, capture_output=True, text=True)
    step = Step(name=name, seconds=time.time() - started)
    if expect_tests:
        step.passed, step.failed, step.skipped, step.errors = _counts(result.stdout)
        # pytest's exit code is the signal; the counts are for the report. A
        # cold clone may skip freely, and 0 means nothing failed or errored —
        # including collection errors, which produce no FAILED lines at all.
        step.ok = result.returncode == 0
        if not step.ok:
            step.detail = "\n".join(
                l for l in result.stdout.splitlines()
                if l.startswith("FAILED") or l.startswith("ERROR"))[:2000]
    else:
        step.ok = result.returncode == 0
        step.detail = (result.stdout.strip().splitlines() or [""])[-1][:160]
        if not step.ok:
            step.detail = (result.stderr.strip() or result.stdout.strip())[-400:]
    return step


def check(fetch: bool = False, keep: bool = False,
          source: Path | None = None) -> ColdCloneReport:
    """Clone to a temporary directory, run with no data, and report."""
    source = source or ROOT
    workdir = Path(tempfile.mkdtemp(prefix="piezo1-cold-"))
    clone = workdir / "clone"
    report = ColdCloneReport(clone=str(clone))

    try:
        report.steps.append(_run(
            "git clone", ["git", "clone", "--quiet", str(source), str(clone)],
            cwd=workdir))
        if not report.steps[-1].ok:
            return report

        for name in ("ref", "data"):
            shutil.rmtree(clone / name, ignore_errors=True)

        report.steps.append(_run(
            "suite on the empty clone", [sys.executable, "-m", "pytest"],
            cwd=clone, expect_tests=True))

        if fetch:
            report.steps.append(_run(
                "python -m piezo1.io.fetch",
                [sys.executable, "-m", "piezo1.io.fetch"], cwd=clone))
            report.steps.append(_run(
                "suite with data", [sys.executable, "-m", "pytest"],
                cwd=clone, expect_tests=True))
        return report
    finally:
        if not keep:
            shutil.rmtree(workdir, ignore_errors=True)
        else:
            print(f"\nclone kept at {clone}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--fetch", action="store_true",
                        help="also download the data and re-run the suite")
    parser.add_argument("--keep", action="store_true",
                        help="leave the temporary clone in place")
    args = parser.parse_args()

    print("Running the project from an empty clone.\n")
    report = check(fetch=args.fetch, keep=args.keep)
    for step in report.steps:
        print(step.line())
        if not step.ok and step.detail:
            for line in step.detail.splitlines()[:12]:
                print(f"         {line}")

    print(f"\n{report.summary()}")
    if report.failures:
        print("\nOn an empty clone a data-dependent test must SKIP, not fail.\n"
              "A failure here is a reproducibility bug: it means something "
              "works only because\nthe data happens to be on disk. See "
              "docs/NEGATIVE_RESULT_PROTOCOL.md and Round 60.")
    return 1 if report.failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
