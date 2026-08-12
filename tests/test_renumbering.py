"""Applying the numbering corrections, and the case that must stay untouched.

Round 83's identification instrument found two register errors and reported
them. Nothing applied them, so five entries were read at residue numbers that
point at the wrong residue — and everything this project reads by number
(transmembrane helices, domain boundaries, variants) was wrong inside the
affected range.

A renumberer is a rewriter, which makes it the most dangerous kind of
instrument here: applied to a file that is already right it would silently
corrupt it. So the load-bearing test is the **null** — 8YEZ resolves the same
767-857 region as the four entries that carry the error and must come back
untouched — and the positive cases are checked by the identity they reach, not
by a shift copied out of the roadmap.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.core.numbering_check import (apply_renumbering,
                                         canonical_renumbering, renumber)

#: The four human entries carrying the 767-857 register error, the splice case,
#: and the controls that must not move.
CARRIES_ERROR = ("8ZU3", "8YFC", "9VMX", "8YFG")
UNAFFECTED = ("8YEZ", "7WLT", "3JAC", "6B3R")


def _load(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


# ------------------------------------------- the null: it must be able to say no

@pytest.mark.parametrize("pdb", UNAFFECTED)
def test_an_entry_that_is_already_right_is_left_alone(pdb):
    """The test that makes the rest mean anything.

    8YEZ resolves the same region as the four entries that carry the register
    error and is numbered correctly throughout. A renumberer that touched it
    would be inventing a correction, and would corrupt every entry it was
    pointed at.
    """
    structure = _load(pdb)
    correction = canonical_renumbering(structure)
    assert not correction.needed, correction.summary()
    assert correction.n_corrected == 0

    fixed, report = renumber(structure)
    assert fixed is structure, "an untouched file must be returned unchanged"
    assert not report.needed
    assert "no correction needed" in report.summary()


# --------------------------------------------- the register error, by identity

@pytest.mark.parametrize("pdb", CARRIES_ERROR)
def test_the_register_error_is_repaired_to_full_agreement(pdb):
    """Checked by the identity reached, not by a remembered shift.

    8YFG is the exception that proves it is measuring rather than pattern
    matching: it carries the R2456H mutation, a genuine residue change, so it
    reaches 0.999 and **must not** reach 1.000.
    """
    correction = canonical_renumbering(_load(pdb))
    assert correction.needed
    assert correction.identity_before == pytest.approx(0.932, abs=0.005)
    floor = 0.998 if pdb == "8YFG" else 0.9999
    assert correction.identity_after >= floor, correction.summary()
    if pdb == "8YFG":
        assert correction.identity_after < 1.0, (
            "8YFG carries a real substitution; a perfect score would mean the "
            "renumbering had absorbed it")


@pytest.mark.parametrize("pdb", CARRIES_ERROR)
def test_exactly_the_91_affected_residues_move(pdb):
    """The roadmap records 91 residues, 767-857. Checked against the file."""
    structure = _load(pdb)
    correction = canonical_renumbering(structure)
    assert correction.n_residues == 91, correction.summary()

    moved = apply_renumbering(structure.res_seq, correction.shifts)
    changed = np.asarray(structure.res_seq)[moved != structure.res_seq]
    assert changed.min() == 767 and changed.max() == 857


def test_the_span_may_cover_unresolved_gaps_and_says_so():
    """The extension runs past the data, which is safe and must be stated.

    It grows while the corrected identity does not fall, and 8ZU3 resolves
    nothing between 713-766 or 858-914 — so the span reaches wider than the
    error while moving only the 91 residues that exist.
    """
    structure = _load("8ZU3")
    correction = canonical_renumbering(structure)
    low, high, shift = correction.shifts[0]
    assert shift == 22
    assert low < 767 and high > 857, "the span did not extend"
    assert correction.n_residues == 91, "but only the resolved residues moved"
    assert "unresolved gaps" in correction.summary()


# ------------------------------------------------------------- the splice case

def test_the_splice_isoform_is_repaired_to_full_agreement():
    """6LQI is deposited in the Piezo1.1 isoform's own continuous numbering.

    0.447 over the whole file, because everything after the splice site is 24
    low. This is the entry `pore_regions` and `conduction_path` refuse outright,
    and the correction is what would let them read it.
    """
    correction = canonical_renumbering(_load("6LQI"))
    assert correction.needed
    assert correction.identity_before == pytest.approx(0.447, abs=0.01)
    assert correction.identity_after >= 0.9999, correction.summary()
    assert correction.shifts[0][2] == 24
    assert correction.n_residues > 700


# --------------------------------------------- what the numbering reaches

def test_correcting_the_numbering_moves_the_dome_and_not_the_pore():
    """The roadmap's own validation ask: say which quantities it reaches.

    The dome is fitted to the **annotated** transmembrane helices, so a
    register error inside the blade puts the wrong atoms in the fit. The pore
    profile is purely geometric and reads no residue number at all, so it must
    not move by a single Angstrom.
    """
    from piezo1.structure.geometry import measure_dome, tm_surface_points
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    structure = _load("8ZU3")
    fixed, correction = renumber(structure)
    assert correction.needed

    def dome(st):
        points, resolved = tm_surface_points(st, "human")
        blocks, _ = protomer_blocks(st)
        return measure_dome(blocks, points).radius_of_curvature, len(resolved)

    def bottleneck(st):
        blocks, _ = protomer_blocks(st)
        return float(pore_profile(st, detect_c3_axis(blocks)).radius.min())

    before, n_before = dome(structure)
    after, n_after = dome(fixed)
    assert n_after > n_before, (
        "the correction should recover helices that were being missed")
    assert abs(after - before) > 1.0, (
        f"the dome radius should move: {before / 10:.2f} -> {after / 10:.2f} nm")

    assert bottleneck(fixed) == bottleneck(structure), (
        "the pore profile reads no residue number and must not move")


def test_no_frozen_claim_depends_on_an_affected_entry():
    """So the correction supersedes nothing, and that is worth recording.

    If a claim ever starts using one of these entries, this fails and the
    supersession has to be handled deliberately.
    """
    from piezo1.analysis.claims import CLAIMS

    affected = set(CARRIES_ERROR) | {"6LQI"}
    for claim in CLAIMS:
        text = f"{claim.key} {claim.description}"
        assert not (affected & set(text.split())), (
            f"claim {claim.key} names an entry whose numbering just changed")


def test_applying_a_shift_is_pure_arithmetic_on_the_range():
    """The primitive, checked without any structure at all."""
    numbers = np.array([10, 20, 30, 40])
    out = apply_renumbering(numbers, [(20, 30, 5)])
    assert list(out) == [10, 25, 35, 40]
    assert list(numbers) == [10, 20, 30, 40], "the input was mutated"
