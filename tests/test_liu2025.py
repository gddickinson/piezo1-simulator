"""The Liu et al. 2025 registry, and the panels that put numbers beside theirs.

Structured like ``test_guo2017.py``, because the registry is: the refusals are
half the deliverable, so most of what is checked here is that each one states a
*specific* reason and cannot quietly acquire a callable.

The replications are checked against the paper's own stated numbers rather than
against a stored copy of our answers. That is what makes them a test of the
measurement instead of a regression pin — if our pore profiler or dome fitter
changes, these should move, and the ones that move away from the published
value should fail.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.liu2025 import (PANELS, PAPER, STATUSES, coverage,
                                     not_replicable, panel_by_key, replicate)
from piezo1.analysis.liu2025_panels import (STATES, axis_length,
                                            cap_gate_loop_span,
                                            cavity_volumes, curvature_radius,
                                            load_state, spring_linker_span,
                                            v2476_diagonal)
from piezo1.config import STRUCTURE_DIR


def _state(name: str):
    if not (STRUCTURE_DIR / f"{STATES[name]}.cif").exists():
        pytest.skip(f"{STATES[name]} not downloaded — run python -m piezo1.io.fetch")
    return load_state(name)


# --------------------------------------------------------------- the registry

def test_every_panel_has_a_valid_status():
    for panel in PANELS:
        assert panel.status in STATUSES, panel.key


def test_every_refusal_states_a_specific_reason():
    """A length floor, because "not possible" is not a reason.

    The same guard `test_guo2017` uses: a short refusal reads as an oversight
    and is indistinguishable from one.
    """
    for panel in PANELS:
        if panel.status == "replicated":
            continue
        assert len(panel.reason) > 60, (
            f"{panel.key}: reason is too short to be one — {panel.reason!r}")


def test_a_refusal_can_never_carry_a_callable():
    for panel in not_replicable():
        assert panel.compute is None, panel.key


def test_an_analogue_can_never_read_as_a_replication():
    """The dangerous status. Every analogue must say what it is *not*."""
    for panel in PANELS:
        if panel.status != "analogue":
            continue
        assert panel.reason, panel.key
        result = replicate(panel.key)
        assert result["status"] == "analogue"
        assert result["reason"] == panel.reason


def test_the_refusals_name_what_they_would_need():
    for panel in not_replicable():
        assert panel.needs, f"{panel.key}: refused without saying what is missing"


def test_the_coverage_counts_add_up():
    c = coverage()
    assert c["total"] == len(PANELS)
    assert c["replicated"] + c["analogue"] + c["not_replicable"] == c["total"]
    assert c["not_replicable"] >= 10, (
        "this paper is half electrophysiology and cryo-EM; a low refusal count "
        "means something has been filed as reproducible that is not")


def test_the_paper_is_identified_well_enough_to_find():
    for key in ("pmid", "doi", "journal", "pdb", "numbering"):
        assert PAPER[key], key
    assert "8IXO" in PAPER["pdb"]


def test_the_electrophysiology_is_refused_rather_than_modelled():
    """Figure 1 and 3O are patch clamp, and the gating model must not stand in.

    This project *has* a Markov gating model, fitted to published rates. Using
    it to draw their inactivation panel would be reporting its own input back
    as a replication.
    """
    for key in ("1a", "1c", "3o"):
        panel = panel_by_key(key)
        assert panel.status == "not_replicable"
        assert any("patch" in n for n in panel.needs), panel.needs


def test_the_md_panels_are_refused_and_name_the_simulation():
    """5F and 5G are what the Martini scaffold exists to prepare for."""
    for key in ("5f", "5g"):
        panel = panel_by_key(key)
        assert panel.status == "not_replicable"
        assert any("MD" in n for n in panel.needs), panel.needs


# ------------------------------------------------- the measurements, vs theirs

def test_the_pore_axis_shortens_by_the_published_amount():
    """Figure 2B: 110 A curved, 100 A intermediate."""
    curved = axis_length(_state("PIEZO1-Curved"))
    intermediate = axis_length(_state("S2472E-Intermediate"))
    assert curved == pytest.approx(110.0, abs=3.0)
    assert intermediate == pytest.approx(100.0, abs=5.0)
    assert curved - intermediate > 8.0, (
        "the shortening is the measurement; its size is the claim")


def test_the_tm_gate_diagonal_opens_from_7_to_14_angstrom():
    """Figure 2E, on side chains — the panel's whole point.

    They read this as clearing the 9-12 A a hydrated Na+ needs, so the
    threshold is checked too: the curved state must fall short of it and the
    intermediate must clear it.
    """
    curved = v2476_diagonal(_state("PIEZO1-Curved"))
    intermediate = v2476_diagonal(_state("S2472E-Intermediate"))
    assert curved == pytest.approx(7.0, abs=1.5)
    assert intermediate == pytest.approx(14.0, abs=1.5)
    assert curved < 9.0 < 12.0 < intermediate


def test_the_cap_gate_loops_separate_by_the_published_amounts():
    """Figures 3F and 3H: 4.3 -> 16.2 A and 4.8 -> 12.8 A, across subunits."""
    curved = cap_gate_loop_span(_state("PIEZO1-Curved"))
    intermediate = cap_gate_loop_span(_state("S2472E-Intermediate"))
    assert curved["A2328-P2382"] == pytest.approx(4.3, abs=1.5)
    assert intermediate["A2328-P2382"] == pytest.approx(16.2, abs=2.0)
    assert curved["D2326-E2383"] == pytest.approx(4.8, abs=1.5)
    assert intermediate["D2326-E2383"] == pytest.approx(12.8, abs=2.5)


def test_the_compressed_spring_matches_the_published_separation():
    """Figure 2F: Y2464 at 17 A in the intermediate state."""
    span = spring_linker_span(_state("S2472E-Intermediate"))
    assert span["Y2464"] == pytest.approx(17.0, abs=2.0)


def test_the_cavity_volumes_reproduce_the_direction_and_not_the_value():
    """Figure 2G, filed as an analogue for exactly this reason.

    CV, EV and MV grow into the intermediate state and IV does not. The
    absolute volumes are a solid of revolution and over-estimate, so the test
    asserts the ordering and explicitly refuses to assert the values.
    """
    curved = cavity_volumes(_state("PIEZO1-Curved"))
    intermediate = cavity_volumes(_state("S2472E-Intermediate"))
    for cavity in ("CV", "EV", "MV"):
        assert intermediate[cavity] > curved[cavity], cavity
    assert intermediate["IV"] == pytest.approx(curved["IV"], rel=0.5), (
        "the inner vestibule is the one they report as comparable")


def test_the_curvature_radius_agrees_where_it_was_calibrated_and_not_beyond():
    """Figure 6, and the disagreement recorded rather than adjusted.

    Our sphere fit was calibrated on Guo & MacKinnon's 10.2 nm for a curved
    dome and reproduces it. Asked for a nearly flat surface it saturates:
    they report 117 nm for the flattened state and we give under 25. Fitting a
    sphere to a flat surface is ill-conditioned, and under-estimating a large
    radius is how that fails.
    """
    curved = curvature_radius(_state("PIEZO1-Curved"))["radius_nm"]
    flat = curvature_radius(_state("PIEZO1-Flattened"))["radius_nm"]

    assert curved == pytest.approx(11.0, abs=2.0), (
        "the curved state is where this fitter was calibrated")
    assert flat > curved, "the flattened state must at least come out flatter"
    assert flat < 40.0, (
        f"the fit now gives {flat:.0f} nm against their 117 — if it has stopped "
        f"saturating, the Figure 6 panel can be promoted from analogue")


def test_the_replication_is_ordered_the_way_the_paper_is():
    """Across all four states, not just the two the numbers are quoted for.

    The S2472E-Curved control matters: their claim is that the mutation alone
    does not open the gate, and it should therefore sit with the curved state
    rather than with the intermediate.
    """
    diagonals = {s: v2476_diagonal(_state(s)) for s in STATES}
    assert diagonals["S2472E-Curved"] < diagonals["S2472E-Intermediate"]
    assert abs(diagonals["S2472E-Curved"]
               - diagonals["PIEZO1-Curved"]) < 2.0, (
        "the mutation alone should not open the gate")


def test_running_a_panel_stamps_the_paper_on_the_result():
    result = replicate("2b")
    assert result["paper"]["pmid"] == "39719701"
    assert set(result["result"]) == set(STATES)
    assert all(np.isfinite(v["axis_A"]) for v in result["result"].values())
