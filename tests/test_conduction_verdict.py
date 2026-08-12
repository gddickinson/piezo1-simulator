"""The composed conduction verdict, and the property that makes it safe.

Round 84d evaluated both halves of the wetting verdict on whichever profile the
pathway produced. That refused every entry on the axial route and almost none
on a lateral one, because the Rao score is a sum over lining residues and
truncating the path truncates the sum. Round 84f reads each criterion off the
profile it is calibrated on.

The load-bearing test is not that the new rule separates the states — it is
that **the axial pathway is unchanged on every downloaded entry**, because
every conduction number this project has recorded was computed that way. The
separation is a *check* on the rule, and it is checked against somebody else's
data (Liu et al.'s Figure 5D ordering) rather than against a stored copy of
ours.
"""

from __future__ import annotations

import json

import numpy as np
import pytest

from piezo1.analysis.conduction import RULE, conduction_verdict
from piezo1.analysis.hydration import load_grid, predict_wetting
from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure


def _grid():
    grid = load_grid()
    if not grid.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")
    return grid


def _load(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


def _profile(structure):
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    blocks, _ = protomer_blocks(structure)
    return pore_profile(structure, detect_c3_axis(blocks))


def _piezo1_entries():
    path = STRUCTURE_DIR.parent / "structures.json"
    del path
    from piezo1.config import RESOURCE_DIR

    entries = json.loads((RESOURCE_DIR / "structures.json").read_text())["entries"]
    return [e["pdb"] for e in entries
            if e.get("protein") == "PIEZO1"
            and e.get("state") not in ("predicted", "fragment")]


# ---------------------------------------------------- the rule comes first

def test_the_rule_is_stated_and_names_both_profiles():
    """A model-form change has to say what it is before what it produces."""
    assert "axial" in RULE and "route" in RULE
    assert "cutoff" in RULE, "the reason the chemistry stays on the full axis"
    assert len(RULE) > 120


def test_the_module_records_how_the_rule_was_reached():
    """Three compositions were tried and the third worked. That is the wrong
    order, and the module says so rather than presenting the result as the
    derivation."""
    from piezo1.analysis import conduction

    text = conduction.__doc__
    assert "third" in text and "wrong order" in text
    assert "superseded" in text, (
        "the Round 84d lateral numbers must be marked superseded, not deleted")


# ------------------------------------------- the property that keeps it safe

@pytest.mark.parametrize("pdb", _piezo1_entries())
def test_the_axial_pathway_is_unchanged(pdb):
    """Every recorded conduction number was computed this way.

    Checked on every entry rather than one, because the two profiles being the
    same object is what makes it true and a single entry would not show a case
    where they were not.
    """
    grid = _grid()
    structure = _load(pdb)
    profile = _profile(structure)
    before = predict_wetting(structure, profile, grid=grid)
    after = conduction_verdict(structure, profile, "axial", grid=grid)

    assert after.conductive == before.conductive
    assert after.hydrophobic_gate == before.hydrophobic_gate
    assert after.sterically_occluded == before.sterically_occluded
    assert after.steric_radius == before.min_radius


def test_the_axial_pathway_uses_one_profile_for_both_halves():
    grid = _grid()
    structure = _load("8IXO")
    profile = _profile(structure)
    verdict = conduction_verdict(structure, profile, "axial", grid=grid)
    assert verdict.path.profile is profile


# ------------------------------- calibration: the composition can say either

def test_a_pore_can_be_refused_on_chemistry_with_a_passable_route():
    """The case the old composition could not produce.

    7WLT's lining dewets — score 1.35 over the whole axis — while its lateral
    route is 2.10 A, comfortably past a water molecule. The verdict must be
    'no' and must be 'no' for the *chemistry*.
    """
    grid = _grid()
    structure = _load("7WLT")
    verdict = conduction_verdict(structure, _profile(structure), "lateral",
                                 grid=grid)
    assert not verdict.conductive
    assert verdict.hydrophobic_gate
    assert not verdict.sterically_occluded, (
        "the route is passable; the refusal is about the lining")
    assert "dewet" in " ".join(verdict.reasons)


def test_a_pore_can_be_refused_on_the_route_with_passable_chemistry():
    """The other half, or the composition asserts nothing.

    8IXO's lining clears the cutoff at 0.31 and its *axial* route is 0.98 A —
    so on the axial pathway it is refused for sterics alone.
    """
    grid = _grid()
    structure = _load("8IXO")
    verdict = conduction_verdict(structure, _profile(structure), "axial",
                                 grid=grid)
    assert not verdict.conductive
    assert not verdict.hydrophobic_gate
    assert verdict.sterically_occluded
    assert "narrows" in " ".join(verdict.reasons)


def test_the_summary_says_which_profile_decided_each_half():
    grid = _grid()
    structure = _load("8IXO")
    text = conduction_verdict(structure, _profile(structure), "lateral",
                              grid=grid).summary()
    assert "full axis" in text
    assert "lateral route" in text


# ------------------------------------------ the check, against their Figure 5D

def test_the_states_separate_the_way_their_figure_5d_orders_them():
    """Liu et al.'s Figure 5D at -0.5 V over a microsecond: PIEZO1-Curved
    passes ~0 Na+, PIEZO1-Flattened ~10, S2472E-Intermediate ~20.

    Reported as a check on the rule, not as its derivation. The part that would
    be easy to get wrong is that the *flattened* state conducts at all — it
    does, in their data and here, because its transmembrane gate is dilated
    while its cap gate is shut.
    """
    from piezo1.analysis.liu2025_permeation import sweep_voltages

    grid = _grid()
    slopes = {}
    for name, pdb in (("curved", "7WLT"), ("flattened", "7WLU"),
                      ("intermediate", "8IXO")):
        structure = _load(pdb)
        verdict = conduction_verdict(structure, _profile(structure), "lateral",
                                     grid=grid)
        if not verdict.conductive:
            slopes[name] = 0.0
            continue
        slopes[name] = sweep_voltages(structure, pathway="lateral",
                                      grid=grid).slope_pS()

    assert slopes["curved"] == 0.0, "the curved state must not conduct"
    assert slopes["flattened"] > 0.0, (
        "the flattened state conducts in their Figure 5D and must here")
    assert slopes["intermediate"] > slopes["flattened"], (
        f"the ordering is wrong: {slopes}")


def test_the_curved_state_is_refused_almost_everywhere():
    """1 of 16, and the exception is a coverage artefact rather than a state.

    3JAC is 4.8 A with 346 unnamed residues; nothing in it is narrow, so its
    score is 0.06 and it passes on chemistry. Named here so a future change
    that fixes it is visible rather than looking like an improvement in the
    model.
    """
    grid = _grid()
    conducting = []
    from piezo1.config import RESOURCE_DIR

    entries = json.loads((RESOURCE_DIR / "structures.json").read_text())["entries"]
    curved = [e["pdb"] for e in entries
              if e.get("protein") == "PIEZO1" and e.get("state") == "curved"]
    for pdb in curved:
        structure = _load(pdb)
        if conduction_verdict(structure, _profile(structure), "lateral",
                              grid=grid).conductive:
            conducting.append(pdb)
    assert conducting == ["3JAC"], (
        f"expected only the 4.8 A entry to slip through; got {conducting}")


# ------------------------------------------------ the mistake, pinned

def test_the_solver_is_given_the_composed_verdict_not_the_axial_one():
    """The bug this had on its first attempt, worth a test of its own.

    Passing the raw full-axis prediction into `solve_pnp` re-imposes the axial
    steric block on a route chosen to avoid it, and the sweep silently returned
    **0.0 pS** for both 8IXO and 7WLU while the verdict said they conduct.
    """
    from piezo1.analysis.liu2025_permeation import sweep_voltages

    grid = _grid()
    for pdb in ("8IXO", "7WLU"):
        sweep = sweep_voltages(_load(pdb), pathway="lateral", grid=grid)
        assert sweep.conducts, f"{pdb} should conduct on the lateral route"
        assert sweep.slope_pS() > 1.0, (
            f"{pdb} came back at {sweep.slope_pS():.1f} pS — the solver is "
            f"being handed the axial verdict again")


def test_an_entry_whose_pathway_is_refused_does_not_silently_conduct():
    """6LQI's splice numbering makes the truncation impossible."""
    grid = _grid()
    structure = _load("6LQI")
    verdict = conduction_verdict(structure, _profile(structure), "lateral",
                                 grid=grid)
    assert verdict.path.refused
    assert not verdict.conductive, (
        "a refused pathway falls back to the whole axis, which is occluded")
