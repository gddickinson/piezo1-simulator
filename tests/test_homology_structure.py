"""The structural side of the family comparison, and the cherry-pick it caught.

Split from ``test_homology.py`` at the 500-line limit and along a real seam:
everything there is about *sequences* and needs no coordinates, and everything
here loads structures and builds elastic networks.

The test that matters most is ``test_a_single_entry_pair_is_not_a_measurement``.
The overlap between PIEZO1's gating mode and dPIEZO's was 0.980 from one pair of
deposited entries and **0.189** from another, and the first had already been
written into ``docs/HOMOLOGY_SEARCH.md`` before the second was run. Every claim
here is therefore made about a *range* over entry pairs, with PIEZO2 as the
positive control that keeps "not stable" from being a verdict the instrument
gives to everything.
"""

from __future__ import annotations

import numpy as np
import pytest


def _registry_path(pdb: str):
    from piezo1.io.registry import load_registry

    record = load_registry().get(pdb)
    if record is None or not record.available:
        pytest.skip(f"{pdb} not downloaded — run python -m piezo1.io.fetch")
    return record.path


# --------------------------------------------------------------------------
# Structural comparison across the family
# --------------------------------------------------------------------------

def test_index_pairing_is_refused_where_the_helix_counts_differ():
    from piezo1.analysis.homology_structure import index_pairing_valid

    valid, reason = index_pairing_valid("human", "mouse_piezo2")
    assert valid and "38" in reason

    for partner in ("worm_piezo", "fly_piezo", "plant_piezo", "dicty_piezo"):
        valid, reason = index_pairing_valid("human", partner)
        assert not valid, partner
        assert "arithmetic rather than structure" in reason


def test_a_monomer_is_refused_with_the_reason():
    from piezo1.analysis.homology_structure import compare_structures

    result = compare_structures("7WLT", "AF-F4IN58-F1-MODEL_V6")
    assert not result.ok
    assert "protomer" in result.error and "monomer" in result.error.lower()


def test_two_entries_of_the_same_protein_are_refused_as_a_state_comparison():
    from piezo1.analysis.homology_structure import compare_structures

    result = compare_structures("7WLT", "7WLU")
    assert not result.ok and "same protein" in result.error


def test_a_single_entry_pair_is_not_a_measurement_of_the_two_proteins():
    """The cherry-pick this nearly published, pinned in the units it appeared in.

    7WLT against 9W7X gives a gating-mode overlap of 0.98 with dPIEZO. 8YEZ
    against the *same* 9W7X gives 0.19. Same two proteins, same method. The
    first number had already been written into docs/HOMOLOGY_SEARCH.md before
    the second was measured.
    """
    from piezo1.analysis.homology_structure import compare_structures

    high = compare_structures("7WLT", "9W7X", with_dome=False)
    low = compare_structures("8YEZ", "9W7X", with_dome=False)
    assert high.ok and low.ok
    assert high.b.renumbered, "9W7X must be renumbered before it is compared"
    assert high.modes.best_overlap > 0.9
    assert low.modes.best_overlap < 0.5
    assert not low.modes.beats_control, (
        "the low pair does not even clear its own shuffled control, which is "
        "what makes quoting the high one a cherry-pick rather than a range")


def test_the_spread_reports_instability_and_piezo2_is_the_positive_control():
    """A 'not stable' verdict means nothing unless something comes back stable.

    PIEZO2 is that something: six entry pairs, every one clearing its control,
    all within 0.2 of each other. Without it, 'not stable' everywhere would be
    equally consistent with an instrument that cannot say yes — the failure
    Round 85 added the positive-control rule for.
    """
    from piezo1.analysis.homology_structure import mode_overlap_spread

    piezo1 = ["7WLT", "8YEZ", "6B3R"]
    paralogue = mode_overlap_spread(piezo1, ["6KG7", "9VEE"], "PIEZO2")
    assert paralogue.n_pairs == 6
    assert paralogue.n_beating_control == paralogue.n_pairs
    assert paralogue.stable, paralogue.summary()
    assert paralogue.low > 0.7

    worm = mode_overlap_spread(piezo1, ["9UOY", "9ZIS"], "PEZO-1")
    assert not worm.stable, worm.summary()
    assert worm.high - worm.low > 0.5
    assert "NOT STABLE" in worm.summary()


def test_the_mode_summary_names_the_partner_it_actually_compared():
    """The label was hard-coded to PIEZO2 and printed for dPIEZO and PEZO-1.

    The numbers were right and the sentence was wrong, which nothing raises on.
    """
    from piezo1.analysis.homology_structure import compare_structures

    fly = compare_structures("7WLT", "9W7X", with_dome=False)
    assert "dPIEZO" in fly.modes.summary()
    assert "PIEZO2" not in fly.modes.summary()


def test_comparable_entries_are_best_resolved_first():
    from piezo1.analysis.homology_structure import comparable_entries
    from piezo1.io.registry import load_registry

    registry = load_registry()
    for _protein, entries in comparable_entries().items():
        sizes = [max(c["n_ca"] for c in registry.get(pdb).protomer_chains)
                 for pdb in entries]
        assert sizes == sorted(sizes, reverse=True), entries


# --------------------------------------------------------------------------
# The sequence viewer
# --------------------------------------------------------------------------

def test_the_viewer_offers_the_whole_family():
    from piezo1.core.numbering_check import REFERENCES
    from piezo1.core.sequences import load_named_sequences

    keys = {s.key for s in load_named_sequences()}
    assert {f"uniprot_{k}" for k in REFERENCES} <= keys


def test_positional_comparison_is_refused_across_numbering_systems():
    """Pairing human 2447 with plant 2447 would report ~2,000 substitutions."""
    from piezo1.core.sequences import compare_sequences, load_named_sequences

    by = {s.key: s for s in load_named_sequences()}
    with pytest.raises(ValueError, match="different systems"):
        compare_sequences(by["uniprot_human"], by["uniprot_plant_piezo"],
                          method="positional")
    # And the global route still works, so the refusal is a redirection.
    comparison = compare_sequences(by["uniprot_human"],
                                   by["uniprot_plant_piezo"])
    assert 0.2 < comparison.identity < 0.35


def test_a_structure_chain_carries_the_numbering_it_is_actually_in():
    from piezo1.core.sequences import load_named_sequences
    from piezo1.core.structure import Structure

    sequences = load_named_sequences(
        Structure.from_file(_registry_path("9UOY")))
    chains = [s for s in sequences if s.key.startswith("structure_")]
    assert chains and all(s.numbering == "worm_piezo" for s in chains)


# --------------------------------------------------------------------------
# The AlphaFold isoform guard
# --------------------------------------------------------------------------

def test_the_alphafold_fetch_refuses_an_isoform_model_offline():
    """Calibrated without the network, on the payload that caused the bug.

    Q9H5I5 returns isoform 3 (709 aa) first and has no canonical model at all.
    The old code took ``entries[0]`` and wrote a quarter-length model of the
    wrong isoform, as a well-formed mmCIF nothing downstream could question.
    """
    import piezo1.io.fetch as F

    payload = [{"uniprotAccession": "Q9H5I5-3", "uniprotSequence": "M" * 709,
                "cifUrl": "https://example.invalid/iso3.cif"},
               {"uniprotAccession": "Q9H5I5-2", "uniprotSequence": "M" * 2689,
                "cifUrl": "https://example.invalid/iso2.cif"}]
    canonical = [e for e in payload
                 if e["uniprotAccession"] == "Q9H5I5"]
    assert not canonical, "the guard's premise: no canonical entry is offered"

    # And the constant records the absence rather than leaving a gap in a list.
    assert F.HUMAN_PIEZO2_ACC in F.ALPHAFOLD_UNAVAILABLE
    assert F.HUMAN_PIEZO2_ACC not in F.ALPHAFOLD_ACCESSIONS
    assert "709" in F.ALPHAFOLD_UNAVAILABLE[F.HUMAN_PIEZO2_ACC]


def test_no_isoform_alphafold_model_reached_the_structure_directory():
    """A model named ``AF-<acc>-<n>-F1`` is an isoform and must not be present."""
    import re

    from piezo1.config import STRUCTURE_DIR

    for path in STRUCTURE_DIR.glob("AF-*.cif"):
        assert not re.match(r"AF-[A-Z0-9]+-\d+-F1", path.name), path.name


def test_the_default_prediction_model_is_named_not_globbed():
    """``sorted(glob)[-1]`` was right only by luck and PIEZO2 sorts after it."""
    from piezo1.analysis import prediction_confidence as PC

    assert PC.DEFAULT_MODEL_ACCESSION == "Q92508"
    assert "AF-Q92508" < "AF-Q9H5I5", (
        "the ordering that made the old glob dangerous has changed; the "
        "reasoning in the comment needs revisiting")


def test_predicted_msf_and_plddt_still_load_the_human_model():
    from piezo1.analysis.prediction_confidence import load_plddt

    values = load_plddt()
    if values is None:
        pytest.skip("AlphaFold model not downloaded")
    assert values.shape == (2521,)
    assert np.count_nonzero(values) > 2400
