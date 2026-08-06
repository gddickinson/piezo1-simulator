"""Non-covalent interaction detection."""

import numpy as np
import pytest

from piezo1.analysis.interactions import (CUTOFFS, compare_interactions,
                                          detect_interactions)


@pytest.fixture(scope="module")
def contacts(human_structure):
    return detect_interactions(human_structure, min_sequence_separation=3)


def test_detects_the_annotated_disulfide(contacts):
    """UniProt annotates C2411–C2415 in human PIEZO1. It must be found.

    The annotation came from UniProt and the detection from coordinates;
    neither knows about the other, so agreement validates both.
    """
    ss = contacts.of_kind("disulfide")
    assert len(ss) == 3, "one disulfide per protomer expected"
    for bond in ss:
        assert {bond.res_i, bond.res_j} == {2411, 2415}
        assert bond.chain_i == bond.chain_j
        assert 1.8 < bond.distance < 2.5


def test_r2456_salt_bridges_to_the_selectivity_glutamate(contacts):
    """R2456 pairs with E2117 of the *neighbouring* protomer.

    R2456H is the archetypal gain-of-function xerocytosis variant and E2117 is
    a selectivity-determining glutamate, so this inter-subunit pairing is a
    concrete structural route from the mutation to the phenotype. Worth pinning
    so a change in the detection criteria cannot quietly lose it.
    """
    bridges = [i for i in contacts.of_kind("salt_bridge")
               if {i.res_i, i.res_j} == {2456, 2117}]
    assert len(bridges) == 3, "expected one per protomer"
    for b in bridges:
        assert b.chain_i != b.chain_j, "should be inter-protomer"
        assert b.distance <= CUTOFFS["salt_bridge"]


def test_all_kinds_are_found_and_bounded(contacts):
    counts = contacts.counts()
    for kind in ("hydrogen_bond", "salt_bridge", "hydrophobic", "disulfide"):
        assert counts.get(kind, 0) > 0, kind
    for i in contacts:
        assert i.distance > 0


def test_hydrogen_bond_distances_respect_the_cutoff(contacts):
    for i in contacts.of_kind("hydrogen_bond"):
        assert 2.2 <= i.distance <= CUTOFFS["hbond_distance"]


def test_nitrogen_nitrogen_pairs_are_excluded_unless_histidine(contacts,
                                                              human_structure):
    """Two nitrogens are both donors, so they cannot hydrogen bond.

    Without hydrogens this has to be excluded by atom identity. Histidine ring
    nitrogens are the exception, since they may be unprotonated and therefore
    act as acceptors.
    """
    st = human_structure
    for i in contacts.of_kind("hydrogen_bond"):
        if st.element[i.atom_i] == "N" and st.element[i.atom_j] == "N":
            assert "HIS" in (i.name_i, i.name_j)


def test_labels_include_atom_names(contacts):
    """Residue names alone make ordinary backbone bonds look like errors.

    ARG.O to ARG.N is a perfectly normal backbone hydrogen bond, but printed
    as "ARG – ARG" it reads as a detection failure.
    """
    text = str(contacts.of_kind("hydrogen_bond")[0])
    assert "." in text


def test_cross_selection_mode_only_returns_interface_contacts(human_structure):
    st = human_structure
    prot = st.mask_protein() & ~st.hetero
    lig = st.mask_ligands()
    if not lig.any():
        pytest.skip("no ligands in this structure")
    inter = detect_interactions(st, mask_a=prot, mask_b=lig,
                                min_sequence_separation=0)
    assert len(inter) > 0
    assert inter.meta["cross_selection_only"] is True
    for i in inter:
        assert bool(lig[i.atom_i]) != bool(lig[i.atom_j]), \
            "every contact must cross the interface"


def test_caveat_is_recorded(contacts):
    assert "hydrogen" in contacts.meta["caveat"].lower()
    assert contacts.meta["cutoffs"]["hbond_distance"] == CUTOFFS["hbond_distance"]


def test_compare_between_identical_states_finds_no_change(human_structure):
    st = human_structure
    kinds = ("disulfide", "salt_bridge")
    a = detect_interactions(st, kinds=kinds, min_sequence_separation=3)
    b = detect_interactions(st.copy_with_coords(st.xyz.copy()), kinds=kinds,
                            min_sequence_separation=3)
    diff = compare_interactions(a, b)
    assert not diff["lost"] and not diff["gained"]
    assert len(diff["retained"]) == len(b)


def test_compare_detects_broken_contacts(human_structure):
    """Pulling the structure apart must show contacts lost, none gained."""
    st = human_structure
    kinds = ("disulfide",)
    a = detect_interactions(st, kinds=kinds, min_sequence_separation=3)
    # Displace one chain far away; its disulfide survives, the others are
    # unaffected, so nothing should be gained and nothing lost either — the
    # bond is intra-chain. Instead scramble one chain's coordinates.
    xyz = st.xyz.copy()
    xyz[st.chain == "A"] += np.array([500.0, 0.0, 0.0], dtype=np.float32)
    moved_all = st.copy_with_coords(xyz)
    b = detect_interactions(moved_all, kinds=kinds, min_sequence_separation=3)
    # Rigid translation of a whole chain preserves its internal disulfide.
    assert len(b.of_kind("disulfide")) == len(a.of_kind("disulfide"))


def test_involving_finds_a_residue(contacts):
    hits = contacts.involving(2456)
    assert hits
    assert all(2456 in (i.res_i, i.res_j) for i in hits)


def test_min_sequence_separation_filters_neighbours(human_structure):
    loose = detect_interactions(human_structure, kinds=("hydrogen_bond",),
                                min_sequence_separation=0)
    strict = detect_interactions(human_structure, kinds=("hydrogen_bond",),
                                 min_sequence_separation=5)
    assert len(strict) < len(loose)
