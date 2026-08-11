"""The drawn blade-to-gate route, and the instrument that stops it lying.

A single line through a protein is the most persuasive picture this
application can draw and the easiest one to over-read. Two things make it
honest, and both are tested before the real structure is touched:

* **the endpoints.** A route from *any* blade residue to the gate is a
  five-step hop from whichever blade residue happens to sit nearest the pore,
  and it never goes near the beam. The lever claim is about force arriving from
  far out, so the source is the most distal blade unit the entry resolves —
  which depends on the entry, so it is measured and reported rather than
  written down;
* **the degeneracy measurement.** `alternative_cost` re-runs the search with
  the drawn route's own edges removed. It is a checking instrument, so it is
  calibrated here on two graphs whose answer is known by construction: one
  where a single bridge is the only way across, and one where the routes are
  interchangeable. A check that cannot return "unique" would make every route
  look degenerate, and a check that cannot return "degenerate" would make
  every route look unique.
"""

from __future__ import annotations

import os

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from piezo1.analysis.allostery import allosteric_path  # noqa: E402
from piezo1.core.annotations import load_annotations  # noqa: E402
from piezo1.ui.path_controller import (MIN_SOURCE_SITES,  # noqa: E402
                                       AllostericPathController,
                                       alternative_cost, path_endpoints)


# --------------------------------------------------------------------------
# Calibrating the degeneracy instrument on known answers
# --------------------------------------------------------------------------

def _chain_with_bridge():
    """Two clusters joined by one pair of residues. One way across, by design."""
    left = np.array([[0.0, y, 0.0] for y in range(4)])          # sites 0-3
    bridge = np.array([[10.0, 0.0, 0.0], [20.0, 0.0, 0.0]])     # sites 4-5
    right = np.array([[30.0, y, 0.0] for y in range(4)])        # sites 6-9
    coords = np.vstack([left, bridge, right])
    n = len(coords)
    dcc = np.full((n, n), 0.9)
    np.fill_diagonal(dcc, 1.0)
    return coords, dcc, np.arange(n)


#: The two clusters of :func:`_chain_with_bridge`, by site index.
LEFT, RIGHT = [0, 1, 2, 3], [6, 7, 8, 9]


def test_a_route_with_one_bridge_is_reported_as_unique():
    """The instrument must be able to say "there is no alternative"."""
    coords, dcc, residues = _chain_with_bridge()
    # A cutoff that joins each cluster internally and each bridge step, but
    # nothing that skips the bridge.
    path = allosteric_path(coords, dcc, LEFT, RIGHT, residues,
                           contact_cutoff=11.0)
    cost = alternative_cost(coords, dcc, LEFT, RIGHT, residues, path.sites,
                            contact_cutoff=11.0)
    assert not np.isfinite(cost), (
        "removing the only bridge left a finite route, so this check cannot "
        "distinguish a unique path from a degenerate one")


def test_a_lattice_of_equivalent_routes_is_reported_as_degenerate():
    """And it must be able to say "this is one of many"."""
    xs, ys = np.meshgrid(np.arange(8) * 5.0, np.arange(5) * 5.0)
    coords = np.stack([xs.ravel(), ys.ravel(), np.zeros(xs.size)], axis=1)
    n = len(coords)
    dcc = np.full((n, n), 0.9)
    np.fill_diagonal(dcc, 1.0)
    residues = np.arange(n)
    source = np.flatnonzero(coords[:, 0] == 0.0).tolist()
    target = np.flatnonzero(coords[:, 0] == 35.0).tolist()
    path = allosteric_path(coords, dcc, source, target, residues,
                           contact_cutoff=6.0)
    cost = alternative_cost(coords, dcc, source, target, residues, path.sites,
                            contact_cutoff=6.0)
    assert np.isfinite(cost)
    assert cost / path.cost < 1.5, (
        f"the best alternative cost {cost / path.cost:.2f}x on a lattice "
        f"where the routes are interchangeable")


def test_the_correlation_matrix_is_left_exactly_as_it_was_found():
    """The alternative search must not touch the correlations at all.

    An earlier version suppressed the drawn route's entries in place and
    restored them afterwards. Restoring is easy to get almost right, and a
    matrix left with a few suppressed entries would bias every later call
    without raising — a slightly different shortest path is still a shortest
    path. Deleting graph edges instead means there is nothing to restore.
    """
    coords, dcc, residues = _chain_with_bridge()
    before = dcc.copy()
    path = allosteric_path(coords, dcc, LEFT, RIGHT, residues,
                           contact_cutoff=11.0)
    alternative_cost(coords, dcc, LEFT, RIGHT, residues, path.sites,
                     contact_cutoff=11.0)
    assert np.array_equal(dcc, before)


def test_an_empty_route_returns_no_answer_rather_than_a_number():
    coords, dcc, residues = _chain_with_bridge()
    assert np.isnan(alternative_cost(coords, dcc, [0], [1], residues, [0],
                                     contact_cutoff=11.0))


# --------------------------------------------------------------------------
# The endpoints
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def annotations():
    return load_annotations("human")


def test_the_source_is_the_most_distal_blade_unit_that_is_resolved(annotations):
    """Which unit that is depends on the entry, so it is derived, not fixed.

    A deposited structure resolves nothing before ~570, so THU1-THU3 are
    absent and the search has to start at THU4. A full-length model resolves
    them, and then it must start further out — otherwise the extra coordinates
    would change nothing about the claim they were added to support.
    """
    deposited = np.tile(np.arange(570, 2521), 3)
    _, _, name = path_endpoints(deposited, annotations)
    assert name.startswith("THU4"), name

    full_length = np.tile(np.arange(1, 2522), 3)
    _, _, name = path_endpoints(full_length, annotations)
    assert name.startswith("THU1"), (
        f"with the distal blade present the search still starts at {name}")


def test_a_few_stray_residues_do_not_count_as_a_blade_unit(annotations):
    """The threshold has to do something, or the "most distal" unit is
    whichever one has a single modelled residue at the edge of the map."""
    residues = np.concatenate([np.arange(1, MIN_SOURCE_SITES),
                               np.tile(np.arange(570, 2521), 3)])
    source, _, name = path_endpoints(residues, annotations)
    assert name.startswith("THU4"), (
        f"{MIN_SOURCE_SITES - 1} sites of THU1 were accepted as the source")
    assert len(source) >= MIN_SOURCE_SITES


def test_the_target_is_the_curated_gate_in_every_protomer(annotations):
    residues = np.tile(np.arange(570, 2521), 3)
    _, target, _ = path_endpoints(residues, annotations)
    gate = set(annotations.group("hydrophobic_gate").residues)
    assert target, "no gate residues found at all"
    assert set(residues[np.asarray(target)]) <= gate
    assert len(target) == 3 * len(gate & set(range(570, 2521))), \
        "the gate was not found in all three protomers"


def test_the_endpoints_follow_the_numbering_the_file_is_in():
    """Mouse and human blade ranges differ, and not by a constant.

    Reading a mouse entry with human ranges is the Round 86 hazard in the one
    place it would produce a confident picture rather than an error.
    """
    human = load_annotations("human")
    mouse = load_annotations("mouse")
    thu4_human = next(d for d in human.domains if d.name.startswith("THU4"))
    thu4_mouse = next(d for d in mouse.domains if d.name.startswith("THU4"))
    assert thu4_human.start != thu4_mouse.start, (
        "the two numbering systems agree on THU4, so this test proves nothing "
        "about which one is being used")


# --------------------------------------------------------------------------
# The refusals
# --------------------------------------------------------------------------

class _FakeWindow:
    def __init__(self, protein="PIEZO1", modes=object(), structure=object()):
        self.structure = structure
        self.record = type("R", (), {"protein": protein,
                                     "numbering_species": "human"})()
        self.modes = modes
        self.viewport = type("V", (), {"scene": object(),
                                       "update": lambda self: None})()
        self.status = ""

    def _set_status(self, text):
        self.status = text


def test_a_paralogue_entry_is_refused_rather_than_approximated():
    """PIEZO1 annotation on a PIEZO2 entry would give a route between two
    arbitrary places, and it would look exactly like a real one."""
    controller = AllostericPathController(_FakeWindow(protein="PIEZO2"))
    refusal = controller.refusal()
    assert "PIEZO2" in refusal
    controller.show(True)
    assert not controller.visible
    assert "PIEZO2" in controller.win.status


def test_without_modes_it_says_where_to_get_them():
    controller = AllostericPathController(_FakeWindow(modes=None))
    assert "Physics" in controller.refusal()


def test_a_piezo1_entry_with_modes_is_not_refused():
    """The refusals must not be the only outcome, or they prove nothing."""
    assert AllostericPathController(_FakeWindow()).refusal() == ""


# --------------------------------------------------------------------------
# What the picture says
# --------------------------------------------------------------------------

def _result(cost=0.22, alternative=0.2183):
    return {"residues": [617, 856, 1046, 2038, 2126, 2454],
            "sites": [0, 1, 2, 3, 4, 5], "cost": cost,
            "correlations": [0.98, 0.99, 0.985, 0.995, 0.99],
            "coords": np.arange(18, dtype=float).reshape(6, 3),
            "alternative_cost": alternative,
            "source_name": "THU4 (TM13-TM16)"}


def test_the_weakest_step_is_the_one_that_stands_out():
    """Averaging the tube's colour would hide the link the route barely made,
    which is the one place a shortest path is worth disbelieving."""
    controller = AllostericPathController(_FakeWindow())
    controller.result = _result()
    colours = controller.step_colors()
    correlations = np.asarray(controller.result["correlations"])
    assert len(colours) == len(correlations)
    weakest = int(np.argmin(correlations))
    strongest = int(np.argmax(correlations))
    from piezo1.ui.path_controller import PATH_COLOR, WEAK_COLOR

    assert np.allclose(colours[weakest], WEAK_COLOR, atol=1e-6)
    assert np.allclose(colours[strongest], PATH_COLOR, atol=1e-6)


def test_a_route_of_uniform_correlation_does_not_paint_a_false_weak_link():
    controller = AllostericPathController(_FakeWindow())
    controller.result = _result()
    controller.result["correlations"] = [0.99] * 5
    colours = controller.step_colors()
    assert len(set(map(tuple, np.round(colours, 6)))) == 1


def test_the_status_line_reports_how_far_from_unique_the_route_is():
    controller = AllostericPathController(_FakeWindow())
    controller.result = _result()
    line = controller.status_line()
    assert "REPRESENTATIVE" in line, \
        "a near-identical alternative was found and the caption did not say so"
    assert "THU4" in line
    assert "not a measured signal" in line


def test_a_genuinely_unique_route_is_described_as_one():
    controller = AllostericPathController(_FakeWindow())
    controller.result = _result(alternative=float("inf"))
    line = controller.status_line()
    assert "only way" in line
    assert "REPRESENTATIVE" not in line
