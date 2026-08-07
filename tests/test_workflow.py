"""Session persistence, provenance-stamped reports, and the headless CLI."""

import json
from pathlib import Path

import pytest

from piezo1 import __version__
from piezo1.analysis.report import (ANALYSES, build_report, collect_provenance)
from piezo1.cli import build_parser, main
from piezo1.io.session import (SESSION_FORMAT, Session, load_session,
                               save_session)


# --------------------------------------------------------------------------
# Sessions
# --------------------------------------------------------------------------

def test_session_round_trip(tmp_path):
    s = Session(structure="8YEZ", style="tube", color_by="chain",
                selected_residues=[2456, 2117], selection_label="R2456H",
                analyses={"dome": {"species": "human"}})
    path = save_session(s, tmp_path / "session.json")
    back = load_session(path)
    assert back.structure == "8YEZ"
    assert back.selected_residues == [2456, 2117]
    assert back.analyses == {"dome": {"species": "human"}}
    assert back.software_version == __version__
    assert back.saved_at


def test_session_stores_no_coordinates_or_results(tmp_path):
    """A session records what you were looking at, not the data.

    Storing results would let a saved file drift silently out of step with the
    code that produced it.
    """
    path = save_session(Session(structure="8YEZ"), tmp_path / "s.json")
    raw = json.loads(path.read_text())
    assert "coords" not in raw and "xyz" not in raw
    assert path.stat().st_size < 4000


def test_session_rejects_a_newer_format():
    with pytest.raises(ValueError, match="newer than this build"):
        Session.from_dict({"format_version": SESSION_FORMAT + 1})


def test_session_ignores_unknown_keys():
    """A file from a newer minor version should still open with what it can."""
    s = Session.from_dict({"structure": "7WLT", "some_future_field": 42,
                           "format_version": SESSION_FORMAT})
    assert s.structure == "7WLT"


def test_missing_session_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        load_session(tmp_path / "absent.json")


def test_session_describe_is_informative():
    s = Session(structure="8YEZ", selected_residues=[1, 2],
                selection_label="gate", analyses={"pore": {}})
    text = s.describe()
    assert "8YEZ" in text and "gate" in text and "pore" in text


# --------------------------------------------------------------------------
# Provenance and reports
# --------------------------------------------------------------------------

def test_provenance_records_the_environment(human_structure):
    p = collect_provenance(human_structure, {"cutoff": 15.0})
    assert p.software_version == __version__
    assert p.timestamp.endswith("+00:00")
    assert p.python and p.platform
    assert "numpy" in p.libraries
    assert p.structure["name"] == "8YEZ"
    assert p.structure["n_atoms"] > 30000
    assert p.parameters["cutoff"] == 15.0
    # Registry metadata should be attached where the structure is known.
    assert p.structure.get("citation")


def test_report_runs_the_named_analyses(human_structure):
    report = build_report(human_structure, analyses=["dome", "pore"])
    assert set(report.results) == {"dome", "pore"}
    assert report.results["dome"]["radius_of_curvature_nm"] > 5
    assert report.results["pore"]["bottleneck_radius_A"] > 0
    assert not report.provenance.warnings


def test_report_records_an_unknown_analysis_rather_than_failing(human_structure):
    report = build_report(human_structure, analyses=["dome", "nonsense"])
    assert "dome" in report.results
    assert any("nonsense" in w for w in report.provenance.warnings)


def test_report_captures_analysis_failure(human_structure, monkeypatch):
    """A failing analysis must be recorded, not swallow the whole report."""
    def boom(*a, **k):
        raise RuntimeError("deliberate")
    monkeypatch.setitem(ANALYSES, "dome", boom)
    report = build_report(human_structure, analyses=["dome", "pore"])
    assert "deliberate" in report.results["dome"]["error"]
    assert report.results["pore"]["bottleneck_radius_A"] > 0
    assert any("dome failed" in w for w in report.provenance.warnings)


def test_report_markdown_and_json_agree(human_structure, tmp_path):
    report = build_report(human_structure, analyses=["dome"])
    md = report.to_markdown(tmp_path / "r.md")
    data = json.loads((report.to_json(tmp_path / "r.json")).read_text())
    assert "Provenance" in md and "8YEZ" in md
    assert data["results"]["dome"]["radius_of_curvature_nm"] == \
        report.results["dome"]["radius_of_curvature_nm"]
    assert data["provenance"]["software_version"] == __version__


def test_report_json_survives_numpy_types(human_structure, tmp_path):
    report = build_report(human_structure, analyses=["modes"])
    path = report.to_json(tmp_path / "modes.json")
    json.loads(path.read_text())          # must not raise


# --------------------------------------------------------------------------
# CLI
# --------------------------------------------------------------------------

def test_json_flag_works_after_the_subcommand():
    """Users type `cli dome 8YEZ --json`; a top-level-only flag rejects that."""
    for argv in (["dome", "8YEZ", "--json"], ["--json", "dome", "8YEZ"]):
        args = build_parser().parse_args(argv)
        assert getattr(args, "json", False) is True, argv
    # Absent unless given; main() supplies the default.
    args = build_parser().parse_args(["dome", "8YEZ"])
    assert getattr(args, "json", False) is False


def test_cli_list_runs(capsys):
    assert main(["list", "--species", "human"]) == 0
    out = capsys.readouterr().out
    assert "8YEZ" in out


def test_cli_dome_emits_json(capsys):
    from piezo1.config import STRUCTURE_DIR

    if not (STRUCTURE_DIR / "8YEZ.cif").exists():
        pytest.skip("8YEZ not downloaded; run python -m piezo1.io.fetch")
    assert main(["dome", "8YEZ", "--json"]) == 0
    data = json.loads(capsys.readouterr().out)
    assert 5.0 < data["radius_of_curvature_nm"] < 20.0
    assert "Haselwandter" in data["reference"]


def test_cli_variants_reports_coverage(capsys):
    assert main(["variants", "--classification", "GoF"]) == 0
    out = capsys.readouterr().out
    assert "resolved in no human structure" in out


def test_cli_unknown_structure_exits_with_advice():
    with pytest.raises(SystemExit) as exc:
        main(["dome", "ZZZZ"])
    assert "fetch" in str(exc.value)


def test_every_analysis_is_reachable_from_the_cli():
    """The CLI and the report share one registry, so they cannot diverge."""
    parser = build_parser()
    for name in ANALYSES:
        args = parser.parse_args([name, "8YEZ"])
        assert args.command == name
