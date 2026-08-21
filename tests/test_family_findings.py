"""The imported PIEZO-family census: the resource, the loader, and the gate.

Nothing here touches coordinates. What it guards is the thing an import can get
wrong that a computation cannot: **quoting somebody else's number after they
have corrected it**, and **joining their per-residue track to the wrong
sequence**. Both fail silently and both produce a plausible answer.

The load-bearing tests are therefore the two negatives — that the build refuses
a source it cannot verify against, and that the constraint track's own amino
acids match this project's copy of Q92508 residue by residue.
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

from piezo1.config import RESOURCE_DIR
from piezo1.core.family import (CONSTRAINT_GENES, ConstraintTrack,
                                load_constraint, load_family_findings)
from piezo1.core.sequence import human_sequence

ROOT = Path(__file__).resolve().parent.parent


@pytest.fixture(scope="module")
def findings():
    return load_family_findings()


@pytest.fixture(scope="module")
def track():
    return load_constraint("PIEZO1")


# --------------------------------------------------------------- the resource

def test_every_finding_says_what_is_done_with_it_here(findings):
    """An imported statement with no destination is a quotation, not a tool."""
    for finding in findings.findings:
        assert finding.here, f"{finding.key} does not say what is done with it"
        assert finding.caveat, f"{finding.key} states no caveat"
        assert len(finding.statement) > 80, f"{finding.key}: statement too thin"
        if not finding.recorded_only:
            module = finding.module
            assert module and module.startswith("analysis."), (
                f"{finding.key} claims to be explored at {finding.here!r}, "
                f"which is not a module path")


def test_every_module_a_finding_names_actually_exists(findings):
    """A finding pointing at a module nobody wrote is worse than no pointer."""
    import importlib

    for finding in findings.findings:
        if finding.recorded_only:
            continue
        module, callable_name = finding.here.rsplit(".", 1)
        imported = importlib.import_module(f"piezo1.{module}")
        assert hasattr(imported, callable_name), (
            f"{finding.key} names {finding.here}, which does not exist")


def test_recorded_only_findings_say_why_rather_than_nothing(findings):
    recorded = [f for f in findings.findings if f.recorded_only]
    assert recorded, "at least the census totals cannot be re-run here"
    for finding in recorded:
        assert len(finding.here) > len("recorded only") + 10, (
            f"{finding.key} is recorded only and does not say why")


def test_every_finding_had_its_numbers_re_read_at_build_time(findings):
    for finding in findings.findings:
        assert finding.n_checks >= 1, (
            f"{finding.key} carries no check, so nothing would notice if the "
            f"census corrected it")


def test_provenance_names_the_source_commit(findings):
    assert findings.provenance["source_project"] == "piezo_genes"
    assert findings.provenance["source_commit"]
    assert findings.source.count("@") == 1


# --------------------------------------------------- the numbering join (gate)

def test_the_constraint_track_is_in_the_numbering_it_claims(track):
    """The one check that would catch an off-by-an-indel join.

    Every residue of the imported track carries the amino acid the census read
    it as. If those disagree with this project's own copy of Q92508 at even one
    position, the track is not in human PIEZO1 numbering and every per-residue
    colouring built on it is wrong somewhere.
    """
    human = human_sequence()
    assert track.length == len(human) == 2521
    mismatches = [i for i, (a, b) in enumerate(zip(track.sequence, human), 1)
                  if a and a != b]
    assert not mismatches, f"{len(mismatches)} residues disagree, first at {mismatches[:5]}"


def test_a_gene_with_no_track_raises_rather_than_returning_an_empty_one():
    """'Not scored' and 'scored zero everywhere' must not share a value."""
    with pytest.raises(KeyError):
        load_constraint("PIEZO4")
    for gene in CONSTRAINT_GENES:
        assert isinstance(load_constraint(gene), ConstraintTrack)


def test_unscored_residues_are_nan_and_not_zero(track):
    """Zero means 'free to change', which is the opposite of 'not measured'."""
    piezo2 = load_constraint("PIEZO2")
    assert np.isnan(piezo2.values).sum() + piezo2.n_scored == piezo2.length
    assert track.value(0) is None and track.value(track.length + 1) is None


def test_the_track_agrees_with_the_findings_it_is_quoted_in(track, findings):
    """R2456 is quoted throughout this subsystem; it must be a real residue."""
    assert track.residue(2456) == "R"
    for pair in findings.equivalent:
        assert track.residue(pair.piezo1) == pair.piezo1_aa


# ------------------------------------------------------------- the build gate

def test_the_build_refuses_a_source_it_cannot_verify_against(tmp_path):
    """The gate that stops a rebuild re-stamping the resource with nothing behind it.

    Pointed at a directory with no census in it, the importer must exit
    non-zero and write nothing — not fall back to the committed copy and report
    success, which would make the provenance line a fiction.
    """
    before = (RESOURCE_DIR / "family_findings.json").read_bytes()
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_family_findings.py"),
         "--source", str(tmp_path)],
        capture_output=True, text=True, timeout=300)
    assert result.returncode != 0
    assert "REFUSED" in result.stderr
    assert (RESOURCE_DIR / "family_findings.json").read_bytes() == before


def test_the_build_verifies_against_the_census_when_it_is_present():
    """The positive half. Skips when the census project is not on disk."""
    source = ROOT.parent / "piezo_genes"
    if not (source / "results").is_dir():
        pytest.skip("the piezo_genes census is not beside this project")
    result = subprocess.run(
        [sys.executable, str(ROOT / "scripts" / "build_family_findings.py"),
         "--check", "--source", str(source)],
        capture_output=True, text=True, timeout=600)
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout


# --------------------------------------------------------- the claimed content

def test_the_fourteen_pathogenic_pore_positions_are_all_kept_by_piezo3(findings):
    positions = findings.pathogenic_pore
    assert len(positions) == 14
    assert all(p.kept_by_piezo3 for p in positions)
    assert {p.gene for p in positions} == {"PIEZO1", "PIEZO2"}


def test_exactly_two_alignment_columns_carry_both_disease_genes(findings):
    """The census's twelve-columns-for-fourteen-positions result, as data."""
    columns = {p.msa_col for p in findings.pathogenic_pore}
    assert len(columns) == 12
    shared = [c for c in columns
              if len({p.gene for p in findings.pathogenic_pore
                      if p.msa_col == c}) == 2]
    assert len(shared) == 2
    assert {p.label for p in findings.equivalent} == {
        "PIEZO1 R2456 = PIEZO2 R2686", "PIEZO1 R2488 = PIEZO2 R2718"}


def test_the_census_scope_is_carried_as_numbers_not_adjectives(findings):
    census = findings.census
    assert census["n_proteins"] > 8000
    assert census["archaeal_hits"] == 0
    # The census's four kingdom counts sum to 15 short of its own stated
    # total. Carried as a recorded gap rather than reconciled: adjusting a
    # number here would make this copy disagree with the project it copies.
    counted = sum(census["by_kingdom"].values())
    assert counted + census["unassigned_to_kingdom"] == census["n_proteins"]
    assert census["unassigned_to_kingdom"] / census["n_proteins"] < 0.01


# ------------------------------------------------- the motif, and its control

def test_the_absent_motif_search_has_a_control_that_finds_something():
    """An absence is only evidence if the search could have found something.

    A scan returning zero looks identical whether the motif is absent or the
    reader is broken, so the positive control runs first: a motif taken from
    human PIEZO1's own sequence must be found, in the three PIEZO1 references
    and nowhere else. Only then is the reported absence of ``PFEW`` worth
    anything.
    """
    from piezo1.analysis.family_motifs import control_motif, motif_scan

    control = control_motif()
    assert control.total >= 3, "the search cannot find a motif known to be there"
    assert set(control.present_in) == {"PIEZO1"}

    absent = motif_scan("PFEW")
    assert absent.total == 0
    assert absent.n_proteins == 10
    assert sum(h.length for h in absent.hits) > 25000
    assert "does not occur" in absent.summary()

    # And it is not a search that returns zero for everything: a single residue
    # that every protein has must be found in every one of them.
    assert len(motif_scan("G").present_in) == absent.n_proteins


def test_what_is_conserved_to_family_depth_is_the_pore_machinery():
    """The positive half of the motif finding.

    The census's claim is that what is conserved across the whole family is not
    a quoted four-letter motif but three short windows around the pore. Checked
    here on the family-depth track: the most conserved windows must land in the
    pore machinery rather than the blades, and they must not all be one peak
    seen five times.
    """
    from piezo1.analysis.family_motifs import deep_windows

    windows = deep_windows(6)
    assert len(windows) == 6
    pore_words = ("anchor", "helix", "c-terminal", "cap")
    hits = [w for w in windows
            if any(word in w.domain.lower() for word in pore_words)]
    assert len(hits) >= 5, [w.domain for w in windows]
    # Non-overlapping: five views of one peak would not be five windows.
    spans = sorted((w.start, w.end) for w in windows)
    for (_, end), (start, _) in zip(spans, spans[1:]):
        assert start > end, spans
    # All of them past residue 2000, which is where the pore module begins.
    assert min(w.start for w in windows) > 2000
