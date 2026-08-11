"""The panel registry: what reproduces, what does not, and whether it says so.

The registry's value is as much in its refusals as in its replications, so most
of this file is about them: that every panel which is not fully replicated
carries a reason, that the reasons are specific rather than boilerplate, and
that an analogue can never be mistaken for the thing it stands in for.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.guo2017 import (PANELS, PAPER, STATUSES, coverage,
                                     not_replicable, panel_by_key, replicate,
                                     replicate_all)


# --------------------------------------------------------------------------
# The registry itself
# --------------------------------------------------------------------------

def test_the_paper_is_identified_correctly():
    """6B3R is Guo & MacKinnon's; 6BPZ is Saotome's, and they were swapped
    in this project's bibliography seed until Round 84."""
    assert PAPER["pdb"] == "6B3R"
    assert PAPER["pmid"] == "29231809"
    assert "mouse" in PAPER["numbering"]

    import json

    from piezo1.config import RESOURCE_DIR
    references = {r["key"]: r for r in json.loads(
        (RESOURCE_DIR / "references.json").read_text())["references"]}
    assert "6B3R" in references["guo2017"]["topic"]
    assert "6BPZ" in references["saotome2018"]["topic"]


def test_every_panel_key_is_unique():
    keys = [p.key for p in PANELS]
    assert len(keys) == len(set(keys))


def test_every_panel_that_is_not_replicated_gives_a_reason():
    """An unexplained refusal is indistinguishable from an oversight."""
    for panel in PANELS:
        if panel.status == "replicated":
            continue
        assert panel.reason, f"{panel.key} has no reason"
        assert len(panel.reason) > 60, (
            f"{panel.key}'s reason is too short to be specific: "
            f"{panel.reason!r}")


def test_every_unreplicable_panel_names_what_it_would_need():
    for panel in not_replicable():
        assert panel.needs, f"{panel.key} does not say what is missing"
        assert panel.compute is None


def test_every_replicated_panel_names_the_module_that_does_it():
    """So INTERFACE.md and the registry cannot drift apart."""
    import importlib

    for panel in PANELS:
        if panel.status == "not_replicable":
            continue
        assert panel.module, f"{panel.key} names no module"
        importlib.import_module(panel.module)


def test_a_panel_cannot_be_registered_with_a_bad_status():
    from piezo1.analysis.guo2017 import Panel

    with pytest.raises(ValueError, match="bad status"):
        Panel("x", "1", "", "t", "s", "invented")
    with pytest.raises(ValueError, match="needs a reason"):
        Panel("x", "1", "", "t", "s", "analogue")
    with pytest.raises(ValueError, match="not_replicable but has a callable"):
        Panel("x", "1", "", "t", "s", "not_replicable", reason="a" * 70,
              compute=lambda *a, **k: {})


def test_the_coverage_summary_adds_up():
    report = coverage()
    assert sum(report["by_status"].values()) == len(PANELS)
    assert set(report["by_status"]) == set(STATUSES)
    assert 0.0 < report["replicable_fraction"] < 1.0, (
        "a coverage of 0 or 1 would mean the registry is not discriminating")
    assert report["by_status"]["not_replicable"] >= 6, (
        "the cryo-EM and liposome panels are genuinely out of reach; if this "
        "drops, either a map has been added or a refusal has been softened")


def test_an_unknown_key_lists_the_known_ones():
    with pytest.raises(KeyError, match="7-S1"):
        panel_by_key("nonesuch")


# --------------------------------------------------------------------------
# The analogues must not read as replications
# --------------------------------------------------------------------------

def test_the_projection_never_claims_to_be_a_class_average(structure_6b3r):
    result = replicate("2ab", structure=structure_6b3r)
    assert result["status"] == "analogue"
    assert result["result"]["is_experimental"] is False
    assert "not a 2D class average" in result["result"]["caveat"]


def test_the_electrostatic_surface_never_claims_to_be_apbs(structure_6b3r):
    result = replicate("4c", structure=structure_6b3r)
    assert result["status"] == "analogue"
    assert "APBS" in result["result"]["caveat"]
    assert "Debye" in result["result"]["method"]


def test_a_panel_that_cannot_be_run_returns_its_reason_not_an_error():
    result = replicate("2e")
    assert result["status"] == "not_replicable"
    assert "half map" in result["reason"]
    assert "result" not in result


# --------------------------------------------------------------------------
# The replications, against what the paper states
# --------------------------------------------------------------------------

def test_figure_7_supplement_reproduces_every_stated_number():
    result = replicate("7-S1")["result"]
    assert result["n_checked"] >= 10
    assert result["n_agreeing"] == result["n_checked"], [
        c for c in result["checks"] if not c["agrees"]]


def test_figure_7d_reproduces_the_four_activation_curves():
    """The curves are fixed by dG/dA, so this checks the two-state model."""
    result = replicate("7d")["result"]
    assert result["t50_values"] == {
        "dG20_dA20": pytest.approx(1.0), "dG20_dA60": pytest.approx(1 / 3),
        "dG40_dA20": pytest.approx(2.0), "dG40_dA60": pytest.approx(2 / 3)}
    for curve in result["curves"]:
        assert curve["p_open_at_t50"] == pytest.approx(0.5, abs=0.01)


def test_figure_6b_agrees_on_which_residues_constrict_but_not_on_the_radii(
        structure_6b3r):
    """A real methodological difference, reported rather than absorbed.

    Guo & MacKinnon used HOLE; this project's profiler is an independent
    implementation. The three named residues come out systematically wider.
    That is worth pinning: if it ever became zero the profiler would have been
    fitted to the paper, and if it grew a lot something has drifted.
    """
    result = replicate("6b", structure=structure_6b3r)["result"]
    assert result["n_compared"] == 3
    assert not result["conductive"], "6B3R is a closed structure"
    assert result["bottleneck_radius_A"] < 1.5
    assert 0.2 < result["mean_offset_A"] < 1.2, (
        f"offset from HOLE is {result['mean_offset_A']:.2f} A")
    for row in result["constrictions"]:
        assert row["measured_radius_A"] > row["published_radius_A"]


def test_figure_6b_agrees_with_the_published_spacing_between_constrictions(
        structure_6b3r):
    """The absolute axis origin differs; the spacing is comparable."""
    result = replicate("6b", structure=structure_6b3r)["result"]
    by_residue = {r["residue"]: r for r in result["constrictions"]}
    measured = (by_residue[2536]["displacement_relative_A"]
                - by_residue[2493]["displacement_relative_A"])
    assert measured == pytest.approx(
        result["published"]["spacing_A"]["M2493_to_P2536"], abs=4.0)


def test_figure_4_supplement_reproduces_the_domain_swapped_interface(
        structure_6b3r):
    result = replicate("4-S1", structure=structure_6b3r)["result"]
    assert result["n_contacts"] >= 3
    assert result["all_domain_swapped"]
    assert result["attractive"]
    kinds = {c["kind"] for c in result["contacts"]}
    assert {"salt_bridge", "hydrogen_bond"} <= kinds, (
        "the paper says both hydrogen bonds and salt bridges")
    # E2257-R1762 is found by both conventions; D2264-R1761 by neither, at
    # 6.43 A charge-centroid separation against a 5.5 A cutoff, though its
    # closest atoms are 4.58 A apart. The two conventions disagree about it
    # and the paper does not say which it used.
    pairs = {tuple(sorted((c["a"][2:], c["b"][2:]))) for c in result["contacts"]}
    assert any("GLU2257" in a and "ARG1762" in b or
               "ARG1762" in a and "GLU2257" in b for a, b in pairs)


def test_figure_4a_reproduces_the_beam_angle(structure_6b3r):
    result = replicate("4a", structure=structure_6b3r)["result"]
    assert result["beam_angle_deg"] == pytest.approx(60.0, abs=8.0)
    assert result["arms_out_of_plane_deg"] == pytest.approx(30.0, abs=12.0)
    assert result["supports_paper"]


def test_figure_3_supplement_supports_the_four_tm_repeat():
    result = replicate("3-S1")["result"]
    assert result["repeat"]["supported"]
    assert result["repeat"]["n_units"] == result["published"]["n_units"]
    assert result["repeat"]["z"] > 3.0


def test_figure_3a_marks_what_this_entry_does_not_resolve(structure_6b3r):
    """Figure 3a greys out TM1-12, and so must ours — from the coordinates."""
    result = replicate("3a", structure=structure_6b3r)["result"]
    assert result["n_helices"] == 38 and result["n_units"] == 9
    unresolved = set(result["unresolved_helices"])
    assert set(range(1, 13)) <= unresolved, (
        "6B3R models nothing before residue 577; TM1-12 must be marked")
    assert 38 not in unresolved and 37 not in unresolved


def test_figure_6_supplement_places_the_cuff_around_the_pore(structure_6b3r):
    """Each named element must actually be near the axis, in this entry."""
    result = replicate("6-S1cd", structure=structure_6b3r)["result"]
    by_id = {row["id"]: row for row in result["elements"]}
    assert {"elbow", "base", "hairpin", "pore_extension"} <= set(by_id)
    for name in ("elbow", "base", "hairpin", "pore_extension"):
        assert by_id[name]["n_ca_resolved"] > 10, name
        assert by_id[name]["radial_min_A"] < 40.0, (
            f"{name} should surround the pore, not sit on a blade")
    # The derived PE helix is the one that has to hug the axis: it is defined
    # as the pore-lining segment, so a range that had picked up the hairpin
    # instead would show here.
    assert by_id["pore_extension"]["radial_max_A"] < 20.0
    assert by_id["hairpin"]["radial_max_A"] > 20.0
    assert by_id["pore_extension"]["confidence"] == "medium"


# --------------------------------------------------------------------------
# The whole run
# --------------------------------------------------------------------------

def test_every_replicable_panel_runs(structure_6b3r):
    report = replicate_all(structure=structure_6b3r)
    failures = {k: v["error"] for k, v in report["panels"].items()
                if "error" in v}
    assert not failures, failures
    for panel in PANELS:
        if panel.compute is None:
            continue
        assert "result" in report["panels"][panel.key], panel.key


def test_a_failing_panel_does_not_abort_the_others():
    """Twelve panels need no structure; one missing entry must not hide them."""
    report = replicate_all(structure=None, keys=["7-S1", "6b", "7d"])
    assert "result" in report["panels"]["7-S1"]
    assert "result" in report["panels"]["7d"]


def test_every_surface_quoting_the_coverage_agrees_with_the_registry():
    """The count appears in five places and a sixth panel must move them all.

    Round 59 added the same guard for the null-result count across the tour,
    the README and the conclusion, after they disagreed. This is the same
    ratchet for the replication coverage: the module docstring, the menu
    tooltip, the README, SCIENCE.md and the in-application help all state it.
    """
    from pathlib import Path

    import piezo1.analysis.guo2017 as registry

    report = coverage()
    total = report["n_panels"]
    replicated = report["by_status"]["replicated"]
    root = Path(__file__).resolve().parents[1]

    spelled = {28: "twenty-eight", 31: "thirty-one", 32: "thirty-two"}
    sources = {
        "module docstring": registry.__doc__,
        "menu tooltip": (root / "piezo1/ui/menus.py").read_text(),
        "help topic": (root / "piezo1/ui/help_topics_paper.py").read_text(),
        "README": (root / "README.md").read_text(),
        "SCIENCE.md": (root / "docs/SCIENCE.md").read_text(),
    }
    for name, text in sources.items():
        assert (str(total) in text or spelled.get(total, "\0") in text), (
            f"{name} does not state the panel count of {total}")
        assert (str(replicated) in text
                or spelled.get(replicated, "\0") in text
                or "Sixteen" in text), (
            f"{name} does not state that {replicated} panels reproduce")


def test_the_registry_is_reachable_from_the_shared_analysis_registry():
    from piezo1.analysis.report import ANALYSES

    assert "guo2017" in ANALYSES
