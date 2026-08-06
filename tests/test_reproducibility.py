"""Packaging, and the guard against documentation drift.

The project states a lot of specific numbers in prose — a dome radius, a mode
overlap, two null results. Prose does not fail a test suite, so a solver rewrite
can leave a document confidently asserting something the code stopped producing.
`piezo1.analysis.claims` is the registry that notices; these tests check the
registry itself is sound, and that the packaging actually describes the project.
"""

import tomllib
from pathlib import Path

import pytest

from piezo1.analysis.claims import (CLAIMS, Claim, claims_by_cost,
                                    verify_claims)

ROOT = Path(__file__).resolve().parent.parent


# --------------------------------------------------------------------------
# The claims registry
# --------------------------------------------------------------------------

def test_every_claim_is_well_formed():
    keys = [c.key for c in CLAIMS]
    assert len(set(keys)) == len(keys), "duplicate claim keys"
    for claim in CLAIMS:
        assert claim.tolerance > 0, f"{claim.key} has no tolerance"
        assert claim.cost in ("fast", "medium", "slow")
        assert claim.description and claim.document
        assert callable(claim.compute)


def test_claimed_documents_exist():
    """A claim pointing at a document that is not there cannot be checked."""
    for claim in CLAIMS:
        assert (ROOT / claim.document).exists(), \
            f"{claim.key} cites missing {claim.document}"


def test_tolerances_are_tight_enough_to_be_meaningful():
    """A tolerance wide enough to admit anything is not a check.

    Every claim must pin its value to better than 10% — otherwise the registry
    would report success while the number moved substantially.
    """
    for claim in CLAIMS:
        if claim.expected == 0.0:
            continue
        assert claim.relative(claim.expected + claim.tolerance) <= 0.10, \
            f"{claim.key} tolerance is {claim.tolerance} on {claim.expected}"


def test_recorded_results_are_frozen():
    """Null results must be marked so drift cannot be fixed by editing prose."""
    frozen = {c.key for c in CLAIMS if c.frozen}
    assert {"round7.p_value", "round7.auroc", "round22.cliffs_delta",
            "round22.auroc"} <= frozen


def test_drift_is_detected():
    """The detector has to detect. A claim that always passes is decoration."""
    wrong = Claim("t.wrong", "deliberately wrong", 1.0, 0.01, "", "README.md",
                  lambda: 2.0)
    right = Claim("t.right", "correct", 1.0, 0.01, "", "README.md",
                  lambda: 1.0)
    results = verify_claims([wrong, right], verbose=False)
    assert not results[0].ok
    assert results[1].ok


def test_missing_data_is_reported_as_skipped_not_failed():
    """A claim that cannot run because data is absent is not documentation
    drift, and conflating the two would make a fresh clone look broken."""
    def explode():
        raise FileNotFoundError("not downloaded")

    claim = Claim("t.missing", "needs data", 1.0, 0.01, "", "README.md",
                  explode)
    result = verify_claims([claim], verbose=False)[0]
    assert not result.ok
    assert "FileNotFoundError" in result.error
    assert result.value is None


def test_cost_filter_is_monotone():
    # Claim is a dataclass holding a callable, so compare by key.
    fast = {c.key for c in claims_by_cost("fast")}
    medium = {c.key for c in claims_by_cost("medium")}
    slow = {c.key for c in claims_by_cost("slow")}
    assert fast <= medium <= slow
    assert len(slow) == len(CLAIMS)
    assert all(c.cost == "fast" for c in claims_by_cost("fast"))


def test_the_fast_claims_still_reproduce():
    """The cheap half of the registry, run for real on every suite run."""
    results = verify_claims(claims_by_cost("fast"), verbose=False)
    drifted = [r for r in results if not r.ok and not r.error]
    assert not drifted, [
        f"{r.claim.key}: documented {r.claim.expected}, computed {r.value}"
        for r in drifted]


# --------------------------------------------------------------------------
# Packaging
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def pyproject():
    path = ROOT / "pyproject.toml"
    if not path.exists():
        pytest.skip("pyproject.toml missing")
    return tomllib.loads(path.read_text())


def test_project_metadata_is_complete(pyproject):
    project = pyproject["project"]
    assert project["name"] == "piezo1-simulator"
    assert project["requires-python"].startswith(">=3.1")
    assert project["dependencies"], "the engine must declare its dependencies"


def test_gui_dependencies_are_optional(pyproject):
    """Everything below `render` runs headless, which is what makes the science
    testable without a display. PyQt must not be a hard requirement."""
    core = " ".join(pyproject["project"]["dependencies"]).lower()
    assert "pyqt" not in core and "moderngl" not in core
    assert "pyqt6" in " ".join(
        pyproject["project"]["optional-dependencies"]["gui"]).lower()


def test_entry_points_resolve(pyproject):
    import importlib
    for table in ("scripts", "gui-scripts"):
        for _name, target in pyproject["project"].get(table, {}).items():
            module, function = target.split(":")
            assert hasattr(importlib.import_module(module), function), target


def test_curated_resources_are_packaged_but_downloads_are_not(pyproject):
    """`piezo1/resources/` is authored content with provenance and ships;
    `ref/` and `data/` are downloads and must never be distributed."""
    package_data = pyproject["tool"]["setuptools"]["package-data"]
    assert "*.json" in package_data["piezo1.resources"]
    assert not any(key.startswith(("ref", "data"))
                   for key in package_data)


def test_environment_locks_exist_and_pin_versions():
    conda_lock = ROOT / "environment.lock.yml"
    pip_lock = ROOT / "requirements.lock.txt"
    if not conda_lock.exists():
        pytest.skip("locks not generated")
    text = conda_lock.read_text()
    assert "name: piezo1" in text
    assert text.count("=") > 100, "lock should pin many versions"
    for line in pip_lock.read_text().splitlines():
        if line.strip() and not line.startswith("#"):
            assert "==" in line or "@" in line, f"unpinned: {line}"


def test_makefile_exposes_the_reproduction_targets():
    makefile = ROOT / "Makefile"
    if not makefile.exists():
        pytest.skip("Makefile missing")
    text = makefile.read_text()
    for target in ("reproduce:", "verify:", "test:", "fetch:", "env:",
                   "sizes:"):
        assert target in text, f"Makefile has no {target} target"
