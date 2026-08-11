"""Putting the pore's own charge on the axis, and the two things that decides.

Which residues count, and which way round the structure is. Both are checked
against known answers before the real thing is measured, because both fail
quietly: a reach test that admits everything reports a pore made of whatever is
nearby, and an orientation read the wrong way turns a cation-selective channel
into an anion-selective one without anything looking odd.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.parameters import PARAMETERS
from piezo1.physics.pore_charge import (AVOGADRO, CHARGE, REACH_KEYS,
                                        ChargedGroup, charged_groups,
                                        cytosolic_end, map_charge, pore_charge)


@pytest.fixture(scope="module")
def framed_open():
    """11ZC in the canonical frame, with its axis and profile."""
    path = STRUCTURE_DIR / "11ZC.cif"
    if not path.exists():
        pytest.skip("11ZC.cif not downloaded — run python -m piezo1.io.fetch")
    from piezo1.structure.frame import apply_frame, canonical_transform
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    st = Structure.from_file(path)
    st = apply_frame(st, canonical_transform(st))
    blocks, _ = protomer_blocks(st)
    axis = detect_c3_axis(blocks)
    return st, axis, pore_profile(st, axis, step=1.0)


class _Profile:
    def __init__(self, z, radius):
        self.z = np.asarray(z, dtype=float)
        self.radius = np.asarray(radius, dtype=float)


# ------------------------------------------- why C-alpha, stated as a measurement

def test_the_one_open_structure_has_no_side_chains(framed_open):
    """The constraint the whole module is shaped around.

    11ZC is the only downloaded entry whose pore is open and the only one
    deposited without side chains — N, CA, C, O and nothing else. So the
    charged group's own coordinates do not exist exactly where they are needed,
    which is why every entry is measured from C-alpha instead of only this one.
    """
    st, _, _ = framed_open
    names = set(np.unique(st.atom_name).tolist())
    assert names == {"N", "CA", "C", "O"}
    assert not (names & {"OE1", "OE2", "NZ", "NH1", "CB"})


# ------------------------------------------------ calibration: the reach test

def test_the_reach_test_uses_the_residues_own_side_chain_length():
    """Arginine reaches further than aspartate, and the test knows it."""
    for name in CHARGE:
        near = ChargedGroup(1, name, "A", 0.0, 5.0, 4.0, CHARGE[name])
        assert near.lines_pore, f"{name} should reach a wall 1 A away"
    close = ChargedGroup(1, "ASP", "A", 0.0, 8.0, 4.0, -1.0)
    far = ChargedGroup(1, "ARG", "A", 0.0, 8.0, 4.0, +1.0)
    assert not close.lines_pore, "Asp cannot reach 4 A"
    assert far.lines_pore, "Arg can"
    assert far.reach > close.reach


def test_a_residue_on_the_outside_of_the_protein_is_excluded(framed_open):
    """The instrument has to be able to say no.

    The blade's basic clusters sit 60 to 120 A from the conduction axis. If any
    of them appeared in the pore's charge the criterion would be measuring
    something other than the pore.
    """
    st, axis, profile = framed_open
    groups = charged_groups(st, profile, axis, mode="lining", species="mouse")
    assert groups, "nothing at all was found, so the test asserts nothing"
    assert max(g.radial for g in groups) < 25.0
    # Mouse 629/630/633 and 1724/1727/1728 are curated basic clusters on the
    # blade; every one of them is far outside.
    numbers = {g.res_seq for g in groups}
    assert not (numbers & {629, 630, 633, 1724, 1727, 1728})


def test_curated_is_a_subset_of_lining(framed_open):
    st, axis, profile = framed_open
    curated = charged_groups(st, profile, axis, mode="curated", species="mouse")
    lining = charged_groups(st, profile, axis, mode="lining", species="mouse")
    keys = {(g.res_seq, g.chain, round(g.z, 3)) for g in lining}
    assert curated
    assert len(lining) > len(curated)
    for group in curated:
        assert (group.res_seq, group.chain, round(group.z, 3)) in keys


def test_three_of_the_four_selectivity_glutamates_are_not_in_the_lumen(framed_open):
    """A measurement, and one the paper that named them agrees with.

    The annotation calls four glutamates selectivity determinants. On the open
    structure only E2461 (mouse E2487) is within side-chain reach of the lumen;
    E2117, E2469 and E2470 are 6.5 to 13 A too far. Coste et al. reached the
    same conclusion about E2117 from function alone — "may not lie in the
    selectivity filter but could be located close enough to the pore to
    allosterically modulate its properties" — so the geometry and the
    electrophysiology agree without either having been fitted to the other.
    """
    st, axis, profile = framed_open
    admitted = {g.res_seq for g in charged_groups(st, profile, axis,
                                                  mode="curated",
                                                  species="mouse")}
    assert 2487 in admitted                    # human E2461
    for excluded in (2133, 2495, 2496):        # human E2117, E2469, E2470
        assert excluded not in admitted
    # And one the annotation never called a selectivity residue is in: mouse
    # E2537 (human E2511), curated as a CTD constriction.
    assert 2537 in admitted


def test_the_numbering_goes_through_the_alignment_not_a_constant(framed_open):
    """Reading a mouse structure with human numbers must find something else.

    The human-to-mouse offset is 16 at E2117 and 26 at E2461, so no constant
    could convert both. Asking for human numbering on this mouse entry is
    therefore a different question, and it has to give a different answer —
    if it gave the same one, the numbering would not be doing anything.
    """
    st, axis, profile = framed_open
    mouse = {g.res_seq for g in charged_groups(st, profile, axis,
                                               mode="curated", species="mouse")}
    human = {g.res_seq for g in charged_groups(st, profile, axis,
                                               mode="curated", species="human")}
    assert mouse != human


# -------------------------------------------- calibration: which end is which

def test_the_cytosolic_end_is_measured_and_can_come_out_either_way(framed_open):
    """Flip the structure and the answer must flip.

    A check that always returns the same index is not a check. The canonical
    frame puts PIEZO1's cytosolic C-terminal domain at negative z, so the first
    slice is the cytosolic one; reflecting the coordinates through the origin
    has to move it to the last.
    """
    import dataclasses

    st, axis, _ = framed_open
    assert cytosolic_end(st, axis) == 0

    flipped = dataclasses.replace(st, xyz=-st.xyz)
    assert cytosolic_end(flipped, axis) == -1


# ------------------------------------------------- the density it turns into

def test_the_smoothing_conserves_total_charge():
    """Total charge is a fact; how it is spread along z is a modelling choice.

    Integrating the density over the lumen must return the elementary charges
    put in, whatever the kernel width and whatever the slice spacing — so no
    conclusion can rest on either.
    """
    for step in (0.5, 1.0, 2.0):
        z = np.arange(-30.0, 30.0 + step, step)
        profile = _Profile(z, np.full_like(z, 5.0))
        groups = [ChargedGroup(1, "GLU", chain, 0.0, 6.0, 5.0, -1.0)
                  for chain in "ABC"]
        for smoothing in (1.5, 3.0, 6.0):
            charge = map_charge(groups, profile, smoothing=smoothing)
            area = np.pi * 5.0 ** 2                       # A^2
            total = np.trapezoid(charge.density * area, z) * AVOGADRO / 1e30
            assert total == pytest.approx(-3.0, rel=1e-3), (step, smoothing)


def test_the_density_is_zero_where_there_is_no_charge():
    z = np.arange(-40.0, 41.0, 1.0)
    profile = _Profile(z, np.full_like(z, 5.0))
    charge = map_charge([ChargedGroup(1, "GLU", "A", 0.0, 6.0, 5.0, -1.0)],
                        profile, smoothing=2.0)
    assert charge.density[0] == pytest.approx(0.0, abs=1e-30)
    assert charge.density[len(z) // 2] < 0.0
    assert charge.net_charge == -1.0


def test_an_empty_set_of_groups_gives_a_pore_the_solver_treats_as_neutral():
    z = np.arange(0.0, 41.0, 1.0)
    charge = map_charge([], _Profile(z, np.full_like(z, 5.0)))
    assert np.all(charge.density == 0.0)
    assert charge.n_groups == 0 and charge.net_charge == 0.0


def test_a_vanishing_radius_cannot_manufacture_an_unbounded_density():
    """One carboxylate divided by a closing lumen would otherwise diverge."""
    z = np.arange(-10.0, 11.0, 1.0)
    radius = np.abs(z) * 0.3                    # closes to zero at the middle
    charge = map_charge([ChargedGroup(1, "GLU", "A", 0.0, 2.0, 0.0, -1.0)],
                        _Profile(z, radius), smoothing=2.0)
    assert np.all(np.isfinite(charge.density))
    floor = PARAMETERS.value("pore.ion_radius")
    peak = charge.peak_density * AVOGADRO / 1e30 * np.pi * floor ** 2
    assert peak < 1.0, "the peak line density cannot exceed the charge itself"


# ------------------------------------------------------ the measured profile

def test_the_measured_charge_on_the_open_structure(framed_open):
    """Round 81's numbers, pinned.

    Six charges on the curated route, all of them carboxylate, from two residue
    numbers; forty-odd on the geometric route with the sign of the net reversed
    by the extracellular cap. The two are not variants of one answer, and the
    pin exists so that stays visible.
    """
    st, axis, profile = framed_open
    curated = pore_charge(st, profile, axis, mode="curated", species="mouse")
    assert curated.n_groups == 6
    assert curated.net_charge == -6.0
    assert {row["res_seq"] for row in curated.residue_summary()} == {2487, 2537}
    assert all(row["copies"] == 3 for row in curated.residue_summary())

    lining = pore_charge(st, profile, axis, mode="lining", species="mouse")
    assert lining.n_groups > 40
    assert lining.net_charge > 0.0, (
        "the geometric route is net positive, which the curated one is not — "
        "if that has changed the two routes no longer disagree in kind")
    assert curated.peak_density > lining.peak_density


def test_the_charge_parameters_are_registered_with_a_stated_reason():
    for key in ("pore_charge.reach_asp", "pore_charge.reach_glu",
                "pore_charge.reach_lys", "pore_charge.reach_arg",
                "pore_charge.smoothing", "pore_charge.max_concentration"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, f"{key} is not registered"
        assert parameter.source_note, key
    reaches = [PARAMETERS.value(f"pore_charge.reach_{n}")
               for n in ("asp", "glu", "lys", "arg")]
    assert reaches == sorted(reaches), "reach should grow with side-chain length"


def test_the_reach_lookup_and_the_key_table_cannot_drift():
    """The four keys are written out twice, on purpose, so this checks them.

    A parameter reached only through a dictionary is invisible to
    `provenance_chain.resolved_keys`, which scans for the call rather than
    running it, and all four of these were reported as read by nothing when
    they were first added. Writing the calls out fixes that and creates a
    second source of truth, which this closes.
    """
    for name, key in REACH_KEYS.items():
        group = ChargedGroup(1, name, "A", 0.0, 0.0, 0.0, CHARGE[name])
        assert group.reach == PARAMETERS.value(key), name
    with pytest.raises(KeyError):
        ChargedGroup(1, "ALA", "A", 0.0, 0.0, 0.0, 0.0).reach
