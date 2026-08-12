"""Where the pore is pinched, and whether the instrument that says so can say no.

This module exists because a refusal read wrongly. "Sterically occluded" is
what the ion-flux animation showed for 17 of 19 deposited PIEZO1 entries, and
it reads as *the gate is shut* — which is not what was measured in any of them.
So the question "which constriction refused it" is now answered, and the
answer is a measuring instrument like any other: calibrated on planted cases
whose answer is known before it is believed against real coordinates.

Three ways it could lie, each checked first:

1. It could report "beyond the gate" whatever it is given. A narrow point
   planted *inside* the gate must come back as ``gate``.
2. It could have the two sides the wrong way round, which is invisible —
   "above" and "below" are equally plausible words. The cytosolic side is
   checked against curated CTD-constriction residues, which the geometry never
   sees.
3. It could locate the gate by trusting the registry's numbering. It does not,
   and a PIEZO2 entry must be refused rather than read in PIEZO1's numbering.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from piezo1.analysis.hydration import WATER_RADIUS_NM
from piezo1.analysis.pore_regions import (Bottleneck, describe_bottleneck,
                                          gate_mask, gate_numbering)
from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.core.annotations import load_annotations

MOUSE_ENTRIES = ("7WLT", "7WLU", "8IMZ", "8IXO", "11ZC")
HUMAN_ENTRIES = ("8YEZ", "8YFC", "8ZU8")


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


@pytest.fixture(scope="module")
def curved():
    st = _load("7WLT")
    return st, _profile(st)


def _planted(profile, index: int, radius: float = 0.4):
    """The same profile with its minimum moved to a chosen slice."""
    wide = np.full_like(profile.radius, 5.0)
    wide[index] = radius
    return dataclasses.replace(profile, radius=wide)


# ------------------------------------------------------ calibration: it can say

def test_a_narrow_point_planted_in_the_gate_is_reported_as_the_gate(curved):
    """The case the real data never produces, so it has to be planted.

    Without this, "never at the gate" across 19 entries would be equally
    consistent with a classifier that cannot return ``gate`` at all.
    """
    structure, profile = curved
    mask = gate_mask(structure, profile)
    assert mask.any(), "7WLT resolves the gate; the mask is the premise here"

    planted = _planted(profile, int(np.flatnonzero(mask)[0]))
    found = describe_bottleneck(structure, planted)
    assert found.narrowest_region == "gate"
    assert not found.blocked_beyond_the_gate


def test_a_narrow_point_planted_outside_the_gate_is_not(curved):
    """The other half: both answers must be reachable from the same input."""
    structure, profile = curved
    mask = gate_mask(structure, profile)
    below_gate = np.flatnonzero(profile.z > profile.z[mask].max())
    if not len(below_gate):
        pytest.skip("no slices beyond the gate on this entry")

    found = describe_bottleneck(structure, _planted(profile, int(below_gate[-1])))
    assert found.narrowest_region in ("above", "below")
    assert found.blocked_beyond_the_gate


# ---------------------------------------------- calibration: the sides are right

@pytest.mark.parametrize("pdb", ["7WLU", "8IMZ"])
def test_the_cytosolic_side_is_the_one_with_the_curated_ctd_constriction(pdb):
    """The direction is measured; this checks it against annotation instead.

    ``cytosolic_end`` orients the profile from where the C-terminal residues
    sit. The curated ``ctd_constriction`` group — mouse M2493/F2494 and
    P2536/E2537 — is independent of that geometry, so if the two sides were
    swapped the constriction named below the gate would be a cap residue and
    this fails. Both entries are chosen because their cytosolic constriction
    *is* one of the curated residues, which is the only case where the check
    can be made without a second judgement call.
    """
    structure = _load(pdb)
    found = describe_bottleneck(structure, _profile(structure))
    curated = set(load_annotations("mouse").group("ctd_constriction").residues)

    assert found.below is not None and found.above is not None
    assert set(found.below.lining) & curated, (
        f"{pdb}: the cytosolic constriction is {found.below.label()}, "
        f"which is not a curated CTD constriction — the sides may be swapped")
    assert not (set(found.above.lining) & curated), (
        f"{pdb}: {found.above.label()} is a CTD residue on the extracellular "
        f"side, which is the signature of an inverted axis")


# --------------------------------------------------- calibration: numbering

@pytest.mark.parametrize("pdb", MOUSE_ENTRIES)
def test_mouse_entries_are_read_in_mouse_numbering(pdb):
    assert gate_numbering(_load(pdb)) == "mouse"


@pytest.mark.parametrize("pdb", HUMAN_ENTRIES)
def test_human_entries_are_read_in_human_numbering(pdb):
    assert gate_numbering(_load(pdb)) == "human"


@pytest.mark.parametrize("pdb", ["6KG7", "9VEE", "9ZIT"])
def test_a_paralogue_is_refused_rather_than_read_in_piezo1s_numbering(pdb):
    """The instrument must be able to return nothing.

    PIEZO2's gate is not at PIEZO1's residue numbers, and reading a PIEZO2
    entry through PIEZO1's annotation is the hazard `numbering_check` exists
    for. A confident location would be worse than none.

    This is the test that found the first implementation wrong. It decided the
    numbering here, by checking the three gate residues' names at the human and
    mouse positions, and called **mouse PIEZO2 human PIEZO1** — Ile, Val and
    Phe are not three distinguishing observations. The job now goes to
    `identify_numbering`, which scores every residue against all six
    references.
    """
    structure = _load(pdb)
    assert gate_numbering(structure) is None

    found = describe_bottleneck(structure, _profile(structure))
    assert found.gate is None
    assert found.reason, "a refusal must say why"
    assert found.narrowest_region == "unknown"
    assert found.sentence() == found.reason


def test_a_splice_isoforms_numbering_is_refused_too():
    """6LQI is PIEZO1 and still cannot be read by residue number.

    It is deposited in the Piezo1.1 isoform's own numbering, +24 after the
    splice site, so its "V2476" is not the gate. Measured rather than trusted:
    the three copies of that residue's side chain sit 31 A apart, against
    7.7 A on 7WLT — so a located gate here would have been a confident,
    entirely reasonable-looking number about the wrong residue.
    """
    structure = _load("6LQI")
    assert gate_numbering(structure) is None
    assert describe_bottleneck(structure, _profile(structure)).gate is None


# ------------------------------------------------------------ the measurement

@pytest.mark.parametrize("pdb", MOUSE_ENTRIES + HUMAN_ENTRIES)
def test_the_block_is_never_at_the_transmembrane_gate(pdb):
    """The finding, pinned so it cannot decay quietly.

    The axial conduction model refuses these entries on a constriction the
    channel is not thought to conduct through: Liu et al. 2025 report that the
    vertical neck stays closed even in the intermediate-open state, because
    the lateral portals carry the current. If a future change ever puts the
    global minimum *at* the gate, that is a real result and this test should
    be the thing that surfaces it.
    """
    structure = _load(pdb)
    found = describe_bottleneck(structure, _profile(structure))
    assert found.gate is not None, f"{pdb}: the gate should be locatable"
    assert found.narrowest_region != "gate", (
        f"{pdb}: the narrowest point is at the gate for the first time")
    assert found.gate.radius / 10.0 >= WATER_RADIUS_NM, (
        f"{pdb}: the gate itself is {found.gate.radius:.2f} A, below the "
        f"water radius — the refusal really would be about the gate")


def test_the_open_prone_entry_reproduces_the_published_gate_dilation():
    """Liu et al. 2025 Figure 2E, on their own measure rather than ours.

    They report the V2476 side-chain diagonal widening from ~7 A in
    PIEZO1-Curved (7WLT) to ~14 A in the intermediate-open S2472E (8IXO),
    "exceeding the size barrier of ~9-12 A required for wetting and ion
    permeation of hydrophobic pore". The coordinates carry it.

    This matters because it is the independent evidence that 8IXO's refusal is
    not a statement about the gate: the gate really has opened, on the metric
    the people who solved it used.
    """
    import itertools

    def diagonal(pdb: str) -> float:
        st = _load(pdb)
        mask = (st.res_seq == 2476) & np.isin(st.atom_name, ["CG1", "CG2"])
        per = {ch: st.xyz[mask & (st.chain == ch)].astype(float).mean(axis=0)
               for ch in sorted(set(st.chain[mask].tolist()))}
        assert len(per) == 3, f"{pdb}: V2476 not resolved in all three protomers"
        return float(np.mean([np.linalg.norm(per[a] - per[b])
                              for a, b in itertools.combinations(per, 2)]))

    curved_diagonal, intermediate_diagonal = diagonal("7WLT"), diagonal("8IXO")
    assert curved_diagonal == pytest.approx(7.7, abs=0.5)
    assert intermediate_diagonal == pytest.approx(14.2, abs=0.7)


def test_the_open_prone_entry_is_refused_beyond_its_open_gate():
    """8IXO clears hydrophobicity, has a wide gate, and is still turned down.

    The refusal is a 0.98 A neck at E2537 — the vertical constriction the paper
    reports as *remaining closed* in this very structure, because the lateral
    portals carry the current. The axial model has no lateral portals in it, so
    it must pass a constriction the channel goes around.
    """
    from piezo1.analysis.hydration import load_grid, predict_wetting

    grid = load_grid()
    if not grid.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")

    structure = _load("8IXO")
    profile = _profile(structure)
    found = describe_bottleneck(structure, profile)
    verdict = predict_wetting(structure, profile, grid=grid)

    assert not verdict.hydrophobic_gate, "8IXO's lining clears the Rao cutoff"
    assert verdict.sterically_occluded, "and it is refused on radius alone"
    assert found.blocked_beyond_the_gate
    assert found.gate.radius / 10.0 > WATER_RADIUS_NM * 2, (
        f"the gate is {found.gate.radius:.2f} A — comfortably passable")
    assert 2537 in (found.below.lining if found.below else ()), (
        f"expected the E2537 neck; got {found.below.label() if found.below else None}")


def test_the_gate_radius_alone_does_not_separate_the_states():
    """The honest limit of our own measure, recorded rather than glossed.

    A first attempt at the test above asserted 8IXO had the widest gate in the
    catalogue. It does not: 7WLU (4.67 A) and 3JAC (4.34 A) are wider and are
    the two worst-resolved entries in the set, at 6.81 and 4.8 A. Gate radius
    is confounded with resolution, which is why the dilation is stated on the
    paper's side-chain measure above and not on this one.
    """
    radii = {}
    for pdb in ("7WLU", "3JAC", "8IXO", "7WLT"):
        structure = _load(pdb)
        found = describe_bottleneck(structure, _profile(structure))
        radii[pdb] = found.gate.radius
    assert radii["7WLU"] > radii["8IXO"] and radii["3JAC"] > radii["8IXO"], (
        f"the confound has gone away: {radii} — re-examine whether gate radius "
        f"can now carry the claim on its own")


def test_locating_the_block_changes_no_verdict(curved):
    """This module reports; it must not decide.

    The wetting prediction is recomputed either side of a `describe_bottleneck`
    call and must be bit-identical, so the diagnostic cannot become a second
    route to a conduction answer.
    """
    from piezo1.analysis.hydration import load_grid, predict_wetting

    grid = load_grid()
    if not grid.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")

    structure, profile = curved
    before = predict_wetting(structure, profile, grid=grid)
    describe_bottleneck(structure, profile)
    after = predict_wetting(structure, profile, grid=grid)
    assert before.score == after.score
    assert before.min_radius == after.min_radius
    assert before.conductive == after.conductive


def test_an_empty_bottleneck_still_produces_a_sentence():
    """Defensive: the caller puts this on the status line, so it must not raise."""
    assert Bottleneck().sentence() == "gate not located"
