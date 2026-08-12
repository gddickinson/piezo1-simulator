"""The tier map is a partition, and the selector is calibrated.

A shorter suite is a selection, and a selection can silently drop things: a
test file no situation runs would decay without anyone deciding that. So the
load-bearing test here is the partition — every ``test_*.py`` in exactly one
tier — which makes forgetting to assign a new test file a failure rather than
a quiet gap. The selector itself is calibrated the way this project calibrates
every instrument: shown to say *no* (a science file deselected from a quick
run, an unknown tier refused loudly) before its yes is believed.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

TESTS_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = TESTS_DIR.parent

from tiers import TIERS, all_files, files_for  # noqa: E402


def test_every_test_file_is_in_exactly_one_tier():
    """The partition. A file in no tier is a test no situation runs; a file
    in two would double-run and make the tier sizes lie."""
    on_disk = {p.name for p in TESTS_DIR.glob("test_*.py")}
    mapped = all_files()

    unassigned = on_disk - mapped
    assert not unassigned, (
        f"test files in no tier — assign them in tests/tiers.py: "
        f"{sorted(unassigned)}")

    vanished = mapped - on_disk
    assert not vanished, (
        f"tier entries with no file on disk — remove or rename them: "
        f"{sorted(vanished)}")

    seen: dict[str, str] = {}
    for tier, (_desc, files) in TIERS.items():
        for name in files:
            assert name not in seen, (
                f"{name} is in both {seen[name]!r} and {tier!r}")
            seen[name] = tier


def test_the_quick_tier_stays_out_of_qt_and_gl():
    """`quick` promises to run anywhere in under a minute; a Qt or GL suite
    misfiled into it breaks both halves of that promise at once."""
    _desc, quick = TIERS["quick"]
    offenders = [n for n in quick
                 if n.startswith(("test_ui_", "test_render_"))
                 and n != "test_tiers.py"]
    assert not offenders, f"Qt/GL suites filed under quick: {offenders}"


def test_an_unknown_tier_is_refused_rather_than_selecting_nothing():
    """A typo that selected nothing would report a green, empty run."""
    try:
        files_for(["quick", "scince"])
    except KeyError as exc:
        assert "scince" in str(exc)
    else:
        raise AssertionError("an unknown tier name was accepted")


def test_the_selector_actually_deselects():
    """The calibration: run the real conftest hook over one quick file and
    one science file, and the science file must not be collected.

    Without this, every tier target could quietly run the whole suite —
    which is the one failure mode that would never be noticed, because the
    only symptom is the time the user was trying to save.
    """
    def collect(suite):
        result = subprocess.run(
            [sys.executable, "-m", "pytest", "--collect-only", "-q",
             "--suite", suite,
             str(TESTS_DIR / "test_tiers.py"), str(TESTS_DIR / "test_anm.py")],
            capture_output=True, text=True, cwd=PROJECT_ROOT, timeout=300)
        assert result.returncode == 0, result.stdout + result.stderr
        return result.stdout

    # Two-sided, or a hook that selects everything would pass one half.
    quick = collect("quick")
    assert "test_tiers.py" in quick and "test_anm" not in quick, (
        "a science file survived a --suite quick run:\n" + quick)
    science = collect("science")
    assert "test_anm.py" in science and "test_tiers.py" not in science, (
        "a quick file survived a --suite science run:\n" + science)
