"""Named sequences, translation and comparison.

Qt-free, like the sequence logic itself. The load-bearing checks are that the
three kinds of sequence are kept distinct and that a downloaded coding sequence
really does translate to the reference protein — a viewer that showed a
back-translation would look identical and mean nothing.
"""

import numpy as np
import pytest

from piezo1.core.sequences import (CODON_TABLE, NamedSequence,
                                   compare_sequences, load_named_sequences,
                                   translate)


def test_genetic_code_is_complete_and_correct():
    assert len(CODON_TABLE) == 64
    assert CODON_TABLE["ATG"] == "M"
    assert CODON_TABLE["TGG"] == "W"
    assert [CODON_TABLE[c] for c in ("TAA", "TAG", "TGA")] == ["*", "*", "*"]
    # Six leucine and six arginine codons, six serine — a standard-code check
    # that catches a table built with the bases in the wrong order.
    from collections import Counter
    counts = Counter(CODON_TABLE.values())
    assert counts["L"] == 6 and counts["R"] == 6 and counts["S"] == 6
    assert counts["M"] == 1 and counts["W"] == 1


def test_translation_handles_partial_and_unknown_codons():
    assert translate("ATGGCC") == "MA"
    assert translate("ATGGC") == "M"          # trailing partial codon dropped
    assert translate("ATGNNN") == "MX"
    assert translate("aug".replace("u", "t").upper()) == "M"


def test_positions_default_to_one_based():
    seq = NamedSequence("k", "label", "MADS")
    assert seq.positions == [1, 2, 3, 4]
    assert seq.at(1) == "M" and seq.at(4) == "S"
    assert seq.at(99) is None
    assert not seq.has_gaps


def test_gapped_sequence_keeps_its_own_numbering():
    """A structure sequence starts at 570 and has holes.

    Treating its offsets as residue numbers renumbers every variant, which is
    the bug this class exists to prevent.
    """
    seq = NamedSequence("s", "structure", "MAD", positions=[570, 571, 600])
    assert seq.has_gaps
    assert seq.at(600) == "D"
    assert seq.index_of(600) == 2
    assert seq.at(572) is None
    assert seq.segment(570, 571) == "MA"


def test_codon_uses_ordinal_position_not_residue_number():
    """A sequence starting at residue 570 does not start at codon 570."""
    seq = NamedSequence("s", "s", "MA", positions=[570, 571], dna="ATGGCC")
    assert seq.codon(570) == "ATG"
    assert seq.codon(571) == "GCC"
    assert seq.codon(999) == ""


def test_positional_comparison_needs_no_alignment():
    a = NamedSequence("a", "a", "MADS")
    b = NamedSequence("b", "b", "MAGS")
    result = compare_sequences(a, b, "positional")
    assert result.n_mismatch == 1
    assert result.n_gap == 0
    assert result.identity == pytest.approx(0.75)
    assert result.differences[0][1] == 3          # residue number
    assert result.differences[0][3:] == ("D", "G")


def test_unknown_method_is_rejected():
    a = NamedSequence("a", "a", "MA")
    with pytest.raises(ValueError, match="unknown method"):
        compare_sequences(a, a, "smith-waterman")


# --------------------------------------------------------------------------
# Against the real downloaded data
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def sequences():
    seqs = load_named_sequences()
    if not seqs:
        pytest.skip("no sequences downloaded")
    return {s.key: s for s in seqs}


def test_human_cds_translates_exactly_to_uniprot(sequences):
    """The strongest check available that the right transcript was fetched.

    If the CDS were the wrong isoform, or a back-translation, or off by a
    codon, this would not be 100%.
    """
    if "cds_human" not in sequences:
        pytest.skip("human CDS not downloaded")
    result = compare_sequences(sequences["uniprot_human"],
                               sequences["cds_human"], "positional")
    assert result.identity == 1.0, result.summary()
    assert result.n_mismatch == 0
    assert len(sequences["cds_human"]) == 2521


def test_mouse_cds_differs_from_trembl_at_three_positions(sequences):
    """Recorded rather than hidden: the Ensembl reference-genome transcript is
    not byte-identical to the TrEMBL entry. Three residues in 2547."""
    if "cds_mouse" not in sequences:
        pytest.skip("mouse CDS not downloaded")
    result = compare_sequences(sequences["uniprot_mouse"],
                               sequences["cds_mouse"], "positional")
    assert result.n_mismatch == 3
    assert result.identity > 0.998
    assert [d[1] for d in result.differences] == [147, 229, 1572]


def test_human_and_mouse_are_about_eighty_percent_identical(sequences):
    result = compare_sequences(sequences["uniprot_human"],
                               sequences["uniprot_mouse"], "global")
    assert 0.75 < result.identity < 0.90, result.summary()
    assert result.n_gap > 0, "an interspecies alignment should need gaps"


def test_positional_comparison_refuses_mixed_numbering_silently_wrong(sequences):
    """Human and mouse share no numbering, so pairing by residue number is
    meaningless — and since Round 89 it **raises** rather than merely looking
    implausible.

    The old behaviour was to run and let the identity collapse to something a
    reader would notice. That was tolerable while the viewer offered two
    sequences; it is not now that it offers nine, because human against plant
    PIEZO would report about two thousand confident substitutions between
    positions that have nothing to do with each other, and 'the number looks
    low' is not a guard.
    """
    with pytest.raises(ValueError, match="different systems"):
        compare_sequences(sequences["uniprot_human"],
                          sequences["uniprot_mouse"], "positional")

    # The redirection has to work, or the refusal is just a removed feature.
    result = compare_sequences(sequences["uniprot_human"],
                               sequences["uniprot_mouse"], "global")
    assert result.identity > 0.8


def test_structure_sequence_is_gapped_and_starts_where_the_model_does(
        human_structure):
    seqs = {s.key: s for s in load_named_sequences(human_structure)}
    chain = seqs.get("structure_A")
    if chain is None:
        pytest.skip("structure chain A not available")
    assert chain.has_gaps, "a deposited model does not resolve everything"
    assert chain.positions[0] >= 500, "8YEZ starts around residue 570"
    assert len(chain.letters) == len(chain.positions)
