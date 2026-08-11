"""The drawn contacts, and the way this feature silently loses most of them.

Nothing here computes anything: the controller calls `detect_interactions` and
draws what comes back. So the tests are about the join between the two, which
is where it went wrong — the first version keyed its colours on "hbond" where
the analysis says "hydrogen_bond", and **7,984 of 9,863 contacts** failed to
draw while the status line confidently reported the rest. Nothing raised.
"""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

pytest.importorskip("PyQt6.QtWidgets")

from piezo1.config import STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.ui.interaction_controller import (DEFAULT_KINDS,  # noqa: E402
                                              KIND_COLORS, KIND_RADIUS,
                                              InteractionController, _family)


@pytest.fixture(scope="module")
def contacts():
    path = STRUCTURE_DIR / "4RAX.cif"
    if not path.exists():
        pytest.skip("4RAX.cif not downloaded — run python -m piezo1.io.fetch")
    from piezo1.analysis.interactions import detect_interactions

    return detect_interactions(Structure.from_file(path))


def test_every_kind_the_analysis_emits_has_a_colour(contacts):
    """The guard for the bug this feature shipped with.

    A kind with no colour is not drawn and nothing says so, which is the worst
    available failure: a picture that looks complete and is missing five
    sixths of its content.
    """
    emitted = {_family(c.kind) for c in contacts.interactions}
    assert emitted, "no contacts found at all; the fixture proves nothing"
    missing = emitted - set(KIND_COLORS)
    assert not missing, f"kinds with no colour, so silently undrawn: {missing}"
    assert set(KIND_COLORS) == set(KIND_RADIUS)


def test_the_hydrogen_bonds_are_the_ones_that_were_lost(contacts):
    """Named, because the failure was specific and would recur the same way."""
    kinds = {c.kind for c in contacts.interactions}
    assert "hydrogen_bond" in kinds, (
        "the analysis calls them 'hydrogen_bond'; if that has been renamed, "
        "KIND_COLORS must be renamed with it or they stop being drawn")
    assert "hbond" not in kinds


def test_the_numerous_kinds_are_off_by_default_and_counted_anyway(contacts):
    """Off because they repeat the ribbon — but never silently dropped.

    Most hydrogen bonds are backbone i to i+4, which is what the cartoon
    already draws; turning them on adds thousands of lines that say nothing
    new. Whatever is hidden must still be counted, or the picture reads as the
    whole inventory.
    """
    controller = InteractionController.__new__(InteractionController)
    controller.result = contacts
    controller.kinds = set(DEFAULT_KINDS)

    assert "hydrogen_bond" not in DEFAULT_KINDS
    assert "hydrophobic" not in DEFAULT_KINDS
    assert set(DEFAULT_KINDS) < set(KIND_COLORS)

    counted = controller.counts()
    assert counted.get("hydrogen_bond", 0) > sum(
        counted.get(k, 0) for k in DEFAULT_KINDS), (
        "if the specific contacts ever outnumber the hydrogen bonds, the "
        "reason for this default has gone")
    line = controller.status_line({k: [] for k in controller.kinds})
    assert "hidden" in line, "what is not drawn has to be reported as not drawn"


def test_the_status_line_carries_both_caveats(contacts):
    controller = InteractionController.__new__(InteractionController)
    controller.result = contacts
    controller.kinds = set(KIND_COLORS)
    line = controller.status_line()
    assert "hydrogens" in line, "no deposited entry has them"
    assert "THIS state" in line, "contacts belong to one conformation"


def test_the_drawing_uses_the_analysis_atom_indices(contacts):
    """A contact must be drawn between the atoms it was found between.

    Re-deriving the endpoints from residue numbers would be a second
    implementation and would put a bond on the wrong atom of a side chain.
    """
    path = STRUCTURE_DIR / "4RAX.cif"
    structure = Structure.from_file(path)
    for contact in contacts.interactions[:200]:
        assert 0 <= contact.atom_i < structure.n_atoms
        assert 0 <= contact.atom_j < structure.n_atoms
        assert int(structure.res_seq[contact.atom_i]) == contact.res_i
        assert int(structure.res_seq[contact.atom_j]) == contact.res_j


def test_a_disulfide_is_drawn_thicker_than_a_hydrophobic_contact():
    """The radii encode how strong a claim each kind is."""
    assert KIND_RADIUS["disulfide"] > KIND_RADIUS["salt_bridge"]
    assert KIND_RADIUS["salt_bridge"] > KIND_RADIUS["hydrophobic"]
