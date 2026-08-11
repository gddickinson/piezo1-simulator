"""The example notebooks: that they exist, that they run, that they are current.

`notebooks/` was created on the first day of the project and stayed empty for
its whole life. Git does not track empty directories, so it was never even in a
clone — it existed only on the machine that made it, which is why nothing ever
noticed.

These tests execute the **committed `.ipynb` files**, not the content modules
that generate them. That is deliberate: the builder already runs its own source
before publishing, so re-running it here would prove nothing about what a
reader downloads. What a reader downloads is the JSON.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
NOTEBOOKS = ROOT / "notebooks"


def notebook_files() -> list[Path]:
    return sorted(NOTEBOOKS.glob("*.ipynb"))


def code_cells(path: Path) -> list[str]:
    document = json.loads(path.read_text())
    return ["".join(cell["source"]) for cell in document["cells"]
            if cell["cell_type"] == "code"]


@pytest.fixture(scope="module")
def data_available() -> bool:
    from piezo1.io.registry import load_registry

    return bool(load_registry().available())


# ------------------------------------------------------------ they exist

def test_the_directory_is_not_empty():
    """The original complaint, made permanent.

    An empty directory is invisible to git, so this is the only thing that can
    notice if the notebooks are ever removed without the folder going with them.
    """
    assert NOTEBOOKS.is_dir(), "the notebooks directory is gone"
    assert notebook_files(), (
        "notebooks/ is empty again; either add notebooks or delete the folder, "
        "because an empty directory is not in the repository at all")
    assert (NOTEBOOKS / "README.md").exists(), "no index for the notebooks"


def test_there_is_one_for_each_thing_a_reader_needs():
    names = {p.stem for p in notebook_files()}
    assert {"01_first_look", "02_gating_motion", "03_pore_to_current",
            "04_variants_and_the_null"} <= names


# --------------------------------------------------------- they are valid

@pytest.mark.parametrize("path", notebook_files(), ids=lambda p: p.stem)
def test_the_file_is_a_valid_notebook(path):
    document = json.loads(path.read_text())
    assert document["nbformat"] == 4
    assert document["cells"], "a notebook with no cells"
    for cell in document["cells"]:
        assert cell["cell_type"] in ("code", "markdown")
        assert isinstance(cell["source"], list)
        if cell["cell_type"] == "code":
            assert "outputs" in cell and "execution_count" in cell


@pytest.mark.parametrize("path", notebook_files(), ids=lambda p: p.stem)
def test_no_outputs_are_committed(path):
    """A stored output is a number nobody recomputes.

    It would go stale silently and read as authoritative while doing it, which
    is the failure this project spends most of its machinery avoiding. The
    notebooks assert their numbers instead.
    """
    document = json.loads(path.read_text())
    for index, cell in enumerate(document["cells"], start=1):
        if cell["cell_type"] == "code":
            assert cell["outputs"] == [], (
                f"{path.name} cell {index} carries a stored output")
            assert cell["execution_count"] is None


@pytest.mark.parametrize("path", notebook_files(), ids=lambda p: p.stem)
def test_every_notebook_explains_itself_before_it_computes(path):
    """A wall of code with no prose teaches nothing."""
    document = json.loads(path.read_text())
    assert document["cells"][0]["cell_type"] == "markdown"
    kinds = [c["cell_type"] for c in document["cells"]]
    assert kinds.count("markdown") >= kinds.count("code") - 1, (
        f"{path.name} is mostly code; these are meant to be read")


# ------------------------------------------------------------- they run

@pytest.mark.parametrize("path", notebook_files(), ids=lambda p: p.stem)
def test_the_committed_notebook_runs(path, data_available):
    """Execute the cells in order in one namespace, as a reader would.

    One namespace rather than one per cell, because the commonest notebook
    defect is a cell that only works thanks to a name defined in a cell the
    author later moved or deleted.

    The notebooks assert their own numbers, so a passing run here is also a
    check that the dome radius, the mode symmetry and the null result have not
    drifted.
    """
    if not data_available:
        pytest.skip("no structures downloaded; run python -m piezo1.io.fetch")

    namespace: dict = {"__name__": "__notebook__"}
    for index, source in enumerate(code_cells(path), start=1):
        try:
            exec(compile(source, f"{path.name}:cell {index}", "exec"), namespace)
        except Exception as exc:                       # pragma: no cover
            pytest.fail(f"{path.name} cell {index} raised "
                        f"{type(exc).__name__}: {exc}")


def test_the_registry_is_left_as_it_was_found(data_available):
    """`03` overrides `pore.step` to show what a parameter does.

    A notebook that left the registry modified would poison every later cell
    and, in this suite, every later test — and `verify_claims` refuses to run
    against a modified registry, so the damage would surface far from here.
    """
    if not data_available:
        pytest.skip("no structures downloaded; run python -m piezo1.io.fetch")
    from piezo1.parameters import PARAMETERS, overrides

    before = PARAMETERS.value("pore.step")
    namespace: dict = {"__name__": "__notebook__"}
    for source in code_cells(NOTEBOOKS / "03_pore_to_current.ipynb"):
        exec(compile(source, "<cell>", "exec"), namespace)

    assert not overrides(), f"the notebook left overrides in place: {overrides()}"
    assert PARAMETERS.value("pore.step") == before


# ------------------------------------------------- they have not drifted

@pytest.mark.parametrize("path", notebook_files(), ids=lambda p: p.stem)
def test_the_committed_file_matches_its_content_module(path):
    """Editing the .ipynb directly is a trap: the next build overwrites it.

    This catches that before the edit is lost, and catches the reverse — a
    content module changed without rebuilding.
    """
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        from build_notebooks import ALL, build
    finally:
        sys.path.pop(0)

    assert path.stem in ALL, (
        f"{path.name} is committed but no content module produces it")
    expected = json.dumps(build(ALL[path.stem]), indent=1) + "\n"
    assert path.read_text() == expected, (
        f"{path.name} differs from its source. Run "
        f"`python scripts/build_notebooks.py` — and if you edited the .ipynb "
        f"by hand, move the change into scripts/notebook_content*.py first.")


def test_every_content_module_notebook_is_committed():
    sys.path.insert(0, str(ROOT / "scripts"))
    try:
        from build_notebooks import ALL
    finally:
        sys.path.pop(0)

    missing = [name for name in ALL if not (NOTEBOOKS / f"{name}.ipynb").exists()]
    assert not missing, f"authored but never built: {missing}"


def test_the_notebooks_are_reachable_from_the_documentation():
    """A folder nobody is told about is the state this started in."""
    readme = (ROOT / "README.md").read_text()
    assert "notebooks/" in readme, "the README does not mention the notebooks"
    index = (NOTEBOOKS / "README.md").read_text()
    for path in notebook_files():
        assert path.name in index, f"{path.name} is missing from the index"
