"""The cold-clone check, and the calibration that makes its "clean" mean something.

Round 60 ran this by hand and found three defects invisible on a developer
machine. Round 74 makes it one command — and the command is worth nothing unless
it can fail, so these tests plant failures and require it to notice.

It earned its place on its first run: it caught two tests written in Rounds 65
and 67 that failed on an empty clone instead of skipping, which is the exact bug
Round 60 found eight times.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "cold_clone_check.py"

sys.path.insert(0, str(ROOT / "scripts"))
import cold_clone_check as ccc  # noqa: E402


# ------------------------------------------- calibrate the count parsing

def test_it_reads_a_normal_summary_line():
    passed, failed, skipped, errors = ccc._counts(
        "collected 900 items\n\n900 passed, 12 skipped in 41.02s\n")
    assert (passed, failed, skipped, errors) == (900, 0, 12, 0)


def test_it_reads_a_summary_with_failures_and_errors():
    passed, failed, skipped, errors = ccc._counts(
        "3 failed, 500 passed, 20 skipped, 2 errors in 30.1s")
    assert (passed, failed, skipped, errors) == (500, 3, 20, 2)


def test_it_finds_the_summary_when_a_progress_bar_comes_last():
    """The bug that made a clean run report as broken.

    With `-q` and no failures the last line is a progress bar, so reading the
    tail found zero of everything.
    """
    output = "800 passed, 40 skipped in 38.36s\n....ss..   [100%]\n"
    assert ccc._counts(output)[0] == 800


def test_missing_counts_do_not_crash_it():
    assert ccc._counts("") == (0, 0, 0, 0)
    assert ccc._counts("no summary here") == (0, 0, 0, 0)


# ---------------------------------------- the pass/fail signal is the exit code

def test_a_failing_suite_is_reported_as_a_failed_step(tmp_path):
    """Planted: the check must notice a test that fails on an empty clone."""
    step = ccc._run("planted failure",
                    [sys.executable, "-c", "import sys; sys.exit(1)"],
                    cwd=tmp_path, expect_tests=True)
    assert not step.ok


def test_a_clean_suite_is_reported_as_ok(tmp_path):
    step = ccc._run("planted pass",
                    [sys.executable, "-c",
                     "print('600 passed, 30 skipped in 1.0s')"],
                    cwd=tmp_path, expect_tests=True)
    assert step.ok
    assert step.passed == 600 and step.skipped == 30


def test_skips_alone_never_count_as_failure(tmp_path):
    """The whole point: a cold clone may skip freely."""
    step = ccc._run("all skipped",
                    [sys.executable, "-c",
                     "print('0 passed, 900 skipped in 1.0s')"],
                    cwd=tmp_path, expect_tests=True)
    assert step.ok, "skipping is expected on an empty clone, not a failure"


def test_the_report_names_which_step_failed():
    report = ccc.ColdCloneReport(steps=[
        ccc.Step(name="git clone", ok=True, seconds=2.0),
        ccc.Step(name="suite on the empty clone", ok=False, seconds=50.0),
    ])
    assert report.failures
    assert "suite on the empty clone" in report.summary()
    assert report.total_seconds == pytest.approx(52.0)


# --------------------------------------------------- the script as shipped

def test_it_does_not_add_a_second_quiet_flag():
    """`pytest.ini` already sets -q; a second makes it -qq and hides the summary."""
    text = SCRIPT.read_text()
    assert '"pytest", "-q"' not in text
    assert "pytest.ini" in text, "the reason should be recorded where it bit"


def test_it_removes_both_data_directories_before_running():
    text = SCRIPT.read_text()
    assert '"ref", "data"' in text or "('ref', 'data')" in text


def test_it_explains_what_a_failure_means():
    text = SCRIPT.read_text()
    assert "must SKIP, not fail" in text
    assert "reproducibility bug" in text


def test_the_script_runs_and_reports(tmp_path):
    """End to end, on the real repository. Slow but this is the whole point."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT)], cwd=ROOT,
        capture_output=True, text=True, timeout=1800)
    assert "git clone" in result.stdout
    assert "suite on the empty clone" in result.stdout
    assert result.returncode == 0, (
        f"the project no longer runs clean from an empty clone:\n"
        f"{result.stdout[-2000:]}")
