"""The MD scaffold, and the one property that makes it safe to ship unrun.

A module that prepares a simulation nobody has run is a standing invitation to
one specific mistake: quoting its manifest as though it were a result. Every
test here is about that boundary. The preparation itself is checked too — the
construct has to be their construct — but the load-bearing assertions are that
no code path produces a :class:`MartiniRun` without a file behind it, and that
asking for results when there are none raises instead of estimating.
"""

from __future__ import annotations

import json

import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.physics.martini import (SEGMENTS, SYSTEM, TOOLCHAIN, MartiniRun,
                                    load_results, prepare, results_available,
                                    write_inputs)


def _load(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


@pytest.fixture(scope="module")
def prepared():
    return prepare(_load("8IXO"), entry="8IXO")


# ------------------------------------------- the boundary: input is not result

def test_a_prepared_system_declares_itself_input_only(prepared):
    manifest = prepared.manifest()
    assert manifest["is_input_only"] is True
    assert manifest["results"] is None
    assert "no simulation has been run" in manifest["note"].lower()


def test_the_summary_cannot_be_mistaken_for_a_result(prepared):
    assert "INPUT ONLY" in prepared.summary()
    assert "nothing has been simulated" in prepared.summary()


def test_a_prepared_directory_has_no_results(tmp_path, prepared):
    directory = write_inputs(prepared, tmp_path / "run")
    assert (directory / "manifest.json").exists()
    assert not results_available(directory), (
        "a prepared directory must never look like a finished one")


def test_asking_for_results_that_do_not_exist_raises(tmp_path, prepared):
    """No fallback, deliberately.

    A function that returned an estimate here would put a continuum number
    where a reader expects a simulated one, which is the single thing this
    module exists to prevent.
    """
    directory = write_inputs(prepared, tmp_path / "run")
    with pytest.raises(FileNotFoundError, match="prepared, not"):
        load_results(directory)


def test_nothing_constructs_a_run_except_load_results():
    """Checked against the source, because this is a discipline not a type.

    `MartiniRun` is a plain dataclass and Python will not stop anyone building
    one. What can be enforced is that *this module* does not.
    """
    import inspect

    from piezo1.physics import martini

    source = inspect.getsource(martini)
    constructions = source.count("MartiniRun(")
    assert constructions == 1, (
        f"MartiniRun is constructed {constructions} times; only load_results "
        f"may build one")
    assert "return MartiniRun(" in inspect.getsource(martini.load_results)


# --------------------------------------------- calibration: results do load

def test_a_real_results_file_loads_and_converts_to_a_current(tmp_path):
    """The other half — the reader must actually work, or the refusals above
    are the only behaviour and the module could never be used."""
    directory = tmp_path / "finished"
    directory.mkdir()
    (directory / "results.json").write_text(json.dumps({
        "voltages_V": [-0.1, -0.5],
        "permeated": [4, 20],
        "duration_ns": 1000.0,
        "source": "test fixture, not a real trajectory",
        "structure": "8IXO"}))

    assert results_available(directory)
    run = load_results(directory)
    assert isinstance(run, MartiniRun)
    currents = run.currents_pA()
    # 20 ions in 1 us is 20 * e / 1e-6 s = 3.2 pA.
    assert currents[1] == pytest.approx(3.2, rel=0.02)
    assert currents[1] > currents[0]


@pytest.mark.parametrize("missing", ["voltages_V", "permeated", "duration_ns",
                                     "source"])
def test_an_incomplete_results_file_is_refused(tmp_path, missing):
    """A results file must say what was run and for how long."""
    payload = {"voltages_V": [-0.5], "permeated": [20], "duration_ns": 1000.0,
               "source": "fixture"}
    payload.pop(missing)
    directory = tmp_path / "bad"
    directory.mkdir()
    (directory / "results.json").write_text(json.dumps(payload))
    with pytest.raises((ValueError, KeyError)):
        load_results(directory)


def test_a_count_per_voltage_is_required(tmp_path):
    directory = tmp_path / "mismatched"
    directory.mkdir()
    (directory / "results.json").write_text(json.dumps({
        "voltages_V": [-0.1, -0.25, -0.5], "permeated": [4, 20],
        "duration_ns": 1000.0, "source": "fixture"}))
    with pytest.raises(ValueError, match="per voltage"):
        load_results(directory)


# ----------------------------------------------------- it is *their* construct

def test_the_construct_is_the_one_they_simulated(prepared):
    """Pore module, beam and lateral plug gate — not the deposited trimer."""
    assert set(SEGMENTS) == {"pore_module", "beam", "lateral_plug_gate"}
    assert SEGMENTS["pore_module"] == (1956, 2547)
    assert SEGMENTS["lateral_plug_gate"] == (1401, 1421)
    for segment in SEGMENTS:
        assert prepared.kept[segment], f"{segment} kept nothing"


def test_the_truncation_keeps_their_construct_and_drops_the_blades(prepared):
    """Their construct is 664 residues per protomer, not the 1,400 8IXO models.

    Measured: 650 of the 664 requested resolve, and the result is 15,699 atoms
    against the trimer's 32,112 — 49%. Half, not a tenth, because the pore
    module is most of what an entry resolves once the distal blade is gone;
    the saving is in the blades, which is exactly what they truncated away.
    """
    whole = _load("8IXO")
    assert prepared.n_atoms < 0.6 * whole.n_atoms
    assert 600 < prepared.n_residues < 700, (
        f"{prepared.n_residues} residues is not their ~664-residue construct")
    assert prepared.n_residues < len(set(whole.res_seq.tolist())) / 2


def test_the_lateral_plug_gate_is_kept(prepared):
    """Their Figure 5G deletes it, so it must be there to delete."""
    kept = prepared.kept["lateral_plug_gate"]
    assert min(kept) >= 1401 and max(kept) <= 1421


def test_a_structure_whose_numbering_cannot_be_read_is_refused():
    """The segments are residue numbers; a wrong reading keeps the wrong 600."""
    with pytest.raises(ValueError, match="numbering"):
        prepare(_load("6KG7"), entry="6KG7")


def test_the_manifest_names_the_tools_rather_than_pretending_to_be_them(prepared):
    manifest = prepared.manifest()
    assert manifest["toolchain"] == list(TOOLCHAIN)
    assert any("martinize" in t for t in TOOLCHAIN)
    assert manifest["conditions"]["voltages_V"] == SYSTEM["voltages_V"]
    assert "STAR Methods" in manifest["conditions"]["source"]


def test_the_voltages_are_the_ones_the_paper_swept(prepared):
    from piezo1.analysis.liu2025_permeation import VOLTAGES

    assert tuple(SYSTEM["voltages_V"]) == tuple(VOLTAGES), (
        "the scaffold and the continuum analogue must sweep the same voltages, "
        "or the two cannot be put on one axis when the run exists")
