"""Every module must import on its own.

Round 78 found `piezo1.analysis.report_tags` unimportable in a fresh
interpreter: it imported a helper from `report`, and `report` imports
`report_tags` at the bottom of the file, so the second import found a
half-initialised module. The whole suite passed regardless, because something
always imported `report` first — and that is exactly why a cycle survives. It
only bites the person who reaches for one module directly, which is a notebook
user or anyone following `docs/NOTEBOOK.md`.

The helper turned out to be a duplicate of `structure.protomers.protomer_blocks`
— identical output on every trimer, differing only in its sentinel for a
non-trimer and in hardcoding the 300-C-alpha floor rather than using
`well_resolved_chains`. Deleting it fixed the duplication and the cycle at once.
"""

from __future__ import annotations

import importlib
import pkgutil
import sys

import pytest

import piezo1

#: Modules that legitimately need something the headless engine does not ship.
#: Each states why, so the exemption is a decision rather than a blanket skip.
NEEDS_EXTRA = {
    "piezo1.ui": "PyQt6",
    "piezo1.render": "moderngl",
}


def package_modules() -> list[str]:
    found = []
    for info in pkgutil.walk_packages(piezo1.__path__, prefix="piezo1."):
        found.append(info.name)
    return sorted(found)


def _optional(name: str) -> str | None:
    for prefix, extra in NEEDS_EXTRA.items():
        if name == prefix or name.startswith(prefix + "."):
            return extra
    return None


def _import_alone(name: str) -> None:
    """Import ``name`` with no piezo1 module already loaded, then put it back.

    Clearing `sys.modules` is what reproduces a fresh interpreter, and it is
    the whole point: with the package already imported, a cycle is invisible
    because the module it needs is sitting in the cache fully built.

    The **restore** is not tidiness, it is correctness, and getting it wrong
    broke seven unrelated tests. `piezo1.parameters.PARAMETERS` is a singleton.
    Other test modules bind it at import time; if this check leaves freshly
    created module objects in the cache, the code under test resolves a
    *different* registry from the one those tests hold, so `set_value` appears
    to do nothing and overrides seem not to take effect. The original objects
    are therefore saved by identity and put back exactly, and the fresh ones
    discarded.
    """
    saved = {m: sys.modules[m] for m in list(sys.modules)
             if m == "piezo1" or m.startswith("piezo1.")}
    for module in saved:
        del sys.modules[module]
    try:
        importlib.import_module(name)
    finally:
        for module in [m for m in list(sys.modules)
                       if m == "piezo1" or m.startswith("piezo1.")]:
            del sys.modules[module]
        sys.modules.update(saved)


@pytest.mark.parametrize("name", package_modules())
def test_the_module_imports_by_itself(name):
    extra = _optional(name)
    if extra is not None:
        pytest.importorskip(extra, reason=f"{name} needs {extra}")
    try:
        _import_alone(name)
    except ImportError as exc:
        if "partially initialized" in str(exc) or "circular" in str(exc):
            pytest.fail(f"{name} has an import cycle: {exc}")
        raise


def test_the_check_leaves_the_singletons_it_borrowed(tmp_path):
    """The damage this file did before the restore was written.

    Seven unrelated tests failed — overrides not taking effect, recorded
    parameters coming back empty — because clearing `sys.modules` left the rest
    of the session resolving a *different* `PARAMETERS` object from the one
    those tests had bound at import time. Identity, not equality, is what
    matters here, so that is what is asserted.
    """
    from piezo1.parameters import PARAMETERS
    from piezo1.analysis import claims

    before_registry = PARAMETERS
    before_module = sys.modules["piezo1.parameters"]
    before_claims = sys.modules["piezo1.analysis.claims"]

    _import_alone("piezo1.analysis.report_tags")

    assert sys.modules["piezo1.parameters"] is before_module
    assert sys.modules["piezo1.analysis.claims"] is before_claims
    assert sys.modules["piezo1.parameters"].PARAMETERS is before_registry
    assert claims.CLAIMS is sys.modules["piezo1.analysis.claims"].CLAIMS


def test_an_override_still_works_after_the_check_has_run():
    """The symptom, reproduced directly rather than left to another file."""
    from piezo1.parameters import PARAMETERS, reset, set_value

    _import_alone("piezo1.analysis.report")
    try:
        set_value("anm.cutoff", 17.0)
        assert PARAMETERS.value("anm.cutoff") == 17.0, (
            "an override stopped taking effect after the import check ran")
    finally:
        reset()


def test_the_check_would_catch_a_real_cycle(tmp_path, monkeypatch):
    """Calibration. A check that cannot fail asserts nothing.

    Two modules that import each other at module scope, built here rather than
    described, so this test fails if `_import_alone` ever stops reproducing a
    fresh interpreter.
    """
    package = tmp_path / "cyclic_probe"
    package.mkdir()
    (package / "__init__.py").write_text("")
    (package / "first.py").write_text("from .second import THING\nVALUE = 1\n")
    (package / "second.py").write_text("from .first import VALUE\nTHING = 2\n")
    monkeypatch.syspath_prepend(str(tmp_path))

    for loaded in [m for m in sys.modules if m.startswith("cyclic_probe")]:
        del sys.modules[loaded]
    with pytest.raises(ImportError) as caught:
        importlib.import_module("cyclic_probe.first")
    assert "partially initialized" in str(caught.value) or \
           "circular" in str(caught.value)


def test_the_analysis_package_is_importable_module_by_module():
    """The one that actually broke, named so the regression is legible."""
    for name in ("piezo1.analysis.report_tags", "piezo1.analysis.report",
                 "piezo1.analysis.claims"):
        _import_alone(name)


def test_no_duplicate_of_protomer_blocks_comes_back():
    """The helper that caused it. One implementation, in `structure`."""
    import inspect

    from piezo1.analysis import report, report_tags
    from piezo1.structure.protomers import protomer_blocks

    for module in (report, report_tags):
        source = inspect.getsource(module)
        assert "def _protomer_blocks" not in source, (
            f"{module.__name__} has its own copy of protomer_blocks again")
        assert "protomer_blocks" in source, f"{module.__name__} must use the shared one"
    assert protomer_blocks.__module__ == "piezo1.structure.protomers"
