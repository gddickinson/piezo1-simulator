"""What deposited MD can contribute — a null about data availability.

Written to fail if the situation improves. If MemProtMD ever ingests a modern
PIEZO structure, or the one it has gains the blade basic clusters, these tests
say so rather than letting a stale limitation stand.

Runs from the cache, so the suite needs no network.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.external_md import (CONTROLS, MD_SOURCES, assess,
                                         lipid_site_coverage,
                                         memprotmd_coverage)


@pytest.fixture(scope="module")
def coverage():
    result = memprotmd_coverage(offline=True)
    if not result.checked:
        pytest.skip(result.note)
    return result


def test_memprotmd_holds_almost_no_piezo_structures(coverage):
    """1 of 21. If this rises, Round 42 is worth revisiting."""
    assert coverage.present == ["3JAC"]
    assert len(coverage.absent) == 20
    assert coverage.fraction < 0.10

    # The structures this project actually uses are the absent ones.
    for pdb in ("7WLT", "7WLU", "6B3R", "8YEZ"):
        assert pdb in coverage.absent


def test_the_probe_had_a_working_control(coverage):
    """Without it, "everything absent" and "the request is wrong" look identical.

    This project has already been caught by that class of error — a PMID that
    resolved to an unrelated paper, a frame mismatch that reported a closed
    channel as conducting. An absence is only evidence if the instrument is
    known to detect presence.
    """
    assert coverage.controls_ok
    assert len(CONTROLS) >= 2


def test_the_one_available_entry_cannot_answer_the_question():
    """3JAC resolves the PIP2 cluster and none of the three blade clusters.

    A simulation of a model that omits the lipid-binding residues cannot report
    their lipid contacts, however good the simulation is.
    """
    sites = lipid_site_coverage("3JAC")
    if "error" in sites:
        pytest.skip(sites["error"])

    assert sites["lipid_residues_resolved"] == 4
    assert sites["lipid_residues_total"] == 15
    assert not sites["can_address_lipid_contacts"]

    groups = sites["groups"]
    pip2 = next(v for k, v in groups.items() if "PIP2" in k)
    assert pip2["complete"], "the PIP2 cluster is the one group it does resolve"
    blades = [v for k, v in groups.items() if "Basic cluster" in k]
    assert len(blades) == 3
    assert all(v["resolved"] == 0 for v in blades)


def test_every_named_source_is_recorded_as_unusable_with_a_reason():
    """The roadmap named three; each is recorded with why, not just that."""
    assert set(MD_SOURCES) == {"MemProtMD", "Zenodo", "GPCRmd"}
    for name, entry in MD_SOURCES.items():
        assert entry["usable"] is False, name
        assert entry["why"], f"{name} must say why, not merely that"
        assert entry["piezo_coverage"], name


def test_the_assessment_refuses_to_claim_the_comparison_is_possible():
    result = assess(offline=True)
    if not result["memprotmd"]["checked"]:
        pytest.skip("coverage not cached")
    assert result["comparison_possible"] is False
    assert "1 of 21" in result["note"]


def test_a_network_failure_is_not_recorded_as_an_absence():
    """The failure mode that would have manufactured this round's conclusion.

    An offline run must return "not checked", never an empty coverage that
    looks like a measured absence.
    """
    result = memprotmd_coverage(offline=True, cache_name="does_not_exist")
    assert not result.checked
    assert result.present == [] and result.absent == []
    assert "offline" in result.note
