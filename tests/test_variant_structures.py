"""What the deposited variant structures can support — pinned, because it is a null.

A negative result decays quietly: someone adds a structure, the denominator
changes, and the recorded limitation stops being true without anything failing.
These tests hold the three facts Round 34 rests on, each measured rather than
asserted, and each stated so that it *fails* if the situation improves — which
is the outcome to want.
"""

from __future__ import annotations

import pytest

from piezo1.analysis.variant_structures import (VARIANT_ENTRIES,
                                                coordinate_fingerprint,
                                                survey_variant_structures)
from piezo1.config import STRUCTURE_DIR
from piezo1.core.structure import Structure


@pytest.fixture(scope="module")
def survey():
    missing = [p for p in VARIANT_ENTRIES
               if not (STRUCTURE_DIR / f"{p}.cif").exists()]
    if missing:
        pytest.skip(f"not downloaded: {missing}")
    return survey_variant_structures()


# ---------------------------------------------------- fact 1: all are closed

def test_no_deposited_human_structure_conducts(survey):
    """So a *difference* in conductance cannot be measured between them.

    This is what stopped Round 34's intended comparison. If it ever fails, an
    open human structure has appeared and the comparison becomes possible —
    which is a reason to revisit the round, not to relax the test.
    """
    assert survey.coverage()["all_closed"]
    for entry in survey.entries:
        assert entry.conductance_pS == 0.0, entry.pdb
        assert entry.bottleneck_A < 1.4, (
            f"{entry.pdb} bottleneck {entry.bottleneck_A:.2f} A is wide enough "
            f"for a cation; it may now conduct")


# ------------------------------------- fact 2: the variants are mostly absent

def test_only_one_entry_actually_shows_its_mutation(survey):
    """A1988 and E756 are unmodelled in the entries named for them.

    A structure that does not resolve its own mutation cannot show what the
    mutation does, however good it is elsewhere.
    """
    coverage = survey.coverage()
    assert coverage["deposited_variant_entries"] == 4
    assert coverage["resolve_their_own_mutation"] == 1

    shown = [e for e in survey.variants if e.shows_mutation]
    assert [e.pdb for e in shown] == ["8YFG"]
    assert shown[0].variant == "R2456H"
    assert shown[0].observed == "HIS"


def test_the_wild_type_residue_is_present_where_the_mutant_is_not(survey):
    """R2456 is arginine in every other entry — the control that makes the
    histidine in 8YFG meaningful rather than a numbering accident."""
    for pdb in ("8YEZ", "8ZU3", "8ZU8"):
        structure = Structure.from_file(STRUCTURE_DIR / f"{pdb}.cif")
        mask = (structure.mask_protein() & ~structure.hetero
                & (structure.res_seq == 2456))
        assert mask.any(), pdb
        assert str(structure.res_name[mask][0]) == "ARG", pdb


# ------------------------------------ fact 3: three entries share one model

def test_three_entries_have_identical_coordinates(survey):
    """8ZU3, 8YFC and 9VMX deposit the same protein model.

    Separate depositions, separate titles, separate files with different
    checksums — so this is a fact about the depositions, not about the
    download. It means they cannot distinguish anything from one another.
    """
    groups = survey.duplicate_groups()
    assert groups == [["8YFC", "8ZU3", "9VMX"]]

    fingerprints = {p: coordinate_fingerprint(
        Structure.from_file(STRUCTURE_DIR / f"{p}.cif"))
        for p in ("8ZU3", "8YFC", "9VMX", "8YEZ", "8ZU8", "8YFG")}
    assert len({fingerprints[p] for p in ("8ZU3", "8YFC", "9VMX")}) == 1
    # ...and they are genuinely different from the others.
    assert len(set(fingerprints.values())) == 4


def test_the_identical_entries_are_different_files():
    """Guards against the far more likely explanation: our own download cache.

    If the fetcher had written one file under three names the coordinates would
    also be identical, and the conclusion above would be about this project
    rather than about the PDB.
    """
    import hashlib

    from piezo1.io.registry import load_registry

    if not load_registry().available():
        pytest.skip("no structures downloaded; run python -m piezo1.io.fetch")

    digests = {}
    for pdb in ("8ZU3", "8YFC", "9VMX"):
        data = (STRUCTURE_DIR / f"{pdb}.cif").read_bytes()
        digests[pdb] = hashlib.md5(data).hexdigest()
        assert f"data_{pdb}".encode() in data[:4096], (
            f"{pdb}.cif does not identify itself as {pdb}")
    assert len(set(digests.values())) == 3, "the three files are not distinct"


# ------------------------------------------------------- what follows from it

def test_only_one_direction_is_represented(survey):
    """Every deposited variant structure is gain-of-function.

    So this route cannot discriminate direction even in principle: there is no
    loss-of-function structure to contrast against. That is the finding, and it
    is why Round 34 reports coverage instead of a comparison.
    """
    coverage = survey.coverage()
    assert coverage["directions_available"] == ["GoF"]
    assert coverage["distinct_variants_shown"] == 1


def test_coverage_is_reported_against_the_curated_set(survey):
    coverage = survey.coverage()
    # The curated set carries 68 variants; these are the ones with a direction.
    assert coverage["curated_variants"] > 30
    assert coverage["informative"] == 1
    assert coverage["informative"] < coverage["curated_variants"] / 30


def test_every_structure_gets_the_same_treatment(survey):
    """A difference between entries must not come from a difference in method."""
    assert len(survey.entries) == len(VARIANT_ENTRIES)
    for entry in survey.entries:
        assert entry.fingerprint
        assert entry.n_protein_atoms > 30000
        assert entry.mechanisms, f"{entry.pdb} recorded no blocking mechanism"


def test_summary_states_the_limitation_rather_than_hiding_it(survey):
    text = survey.summary()
    assert "4 deposited" in text
    assert "GoF" in text and "LoF" in text
    assert "All closed: True" in text
