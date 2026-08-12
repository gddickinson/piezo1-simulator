"""PIEZO2 as the control, and the three ways this comparison could lie.

It could read a file in the wrong numbering — mouse Piezo2 is 2,822 aa against
human PIEZO2's 2,752 and mouse Piezo1's 2,547, so no constant offset relates
any pair. It could pair transmembrane helices by index when the alignment says
otherwise. And it could compare two domes traced over different amounts of
blade and report the coverage difference as a shape difference, which is what
the naive measurement does.

Each is checked before the result is believed, and the last one is not
hypothetical: it changes PIEZO2's dome depth from 8.5 nm to 5.6 nm.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import RESOURCE_DIR, STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.parameters import PARAMETERS
from piezo1.core.numbering_check import (PIEZO2_REFERENCES, REFERENCES,
                                             identify_numbering,
                                             mismatch_blocks)
from piezo1.analysis.paralogue import (compare, dome_comparison,
                                       mode_comparison, paralogue_map,
                                       tm_index_correspondence)


def _require(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


# ------------------------------------------- calibration: which protein is it

def test_every_entry_identifies_its_own_numbering_unambiguously():
    """The known-answer case: an entry must match one reference and no other.

    This is a measurement rather than a label, and it is the one that decides
    whether the transmembrane helices are read in the right place. A file whose
    numbering belongs to a sequence agrees with it at every position; against
    any other sequence it falls to the background rate for unrelated protein.
    """
    expected = {"7WLT": "mouse", "6B3R": "mouse", "11ZC": "mouse",
                "8YEZ": "human", "8ZU3": "human", "6KG7": "mouse_piezo2"}
    floor = PARAMETERS.value("paralogue.min_sequence_identity")
    margin = PARAMETERS.value("paralogue.min_identity_margin")

    for pdb, reference in expected.items():
        identity = identify_numbering(_require(pdb))
        assert identity.reference == reference, identity.summary()
        assert identity.identity > floor, identity.summary()
        assert identity.margin > margin, identity.summary()
        assert identity.confident


def test_the_identification_can_say_the_wrong_thing_if_given_the_wrong_file():
    """An instrument that always agrees with the label asserts nothing.

    Renumbering a PIEZO1 entry by a constant is exactly the mistake the
    identification exists to catch, so it has to fail on one.
    """
    import dataclasses

    st = _require("7WLT")
    shifted = dataclasses.replace(st, res_seq=(st.res_seq + 300).astype(st.res_seq.dtype))
    identity = identify_numbering(shifted)
    assert not identity.confident, identity.summary()
    assert identity.identity < PARAMETERS.value("paralogue.min_sequence_identity")


def test_piezo2_is_identified_as_a_paralogue_and_piezo1_is_not():
    assert identify_numbering(_require("6KG7")).is_piezo2
    assert not identify_numbering(_require("7WLT")).is_piezo2
    assert set(PIEZO2_REFERENCES) < set(REFERENCES)


# ----------------------------------- calibration: the helix pairing by index

def test_the_alignment_confirms_pairing_helices_by_index():
    """37 of 38, with the exception named.

    Coverage matching pairs TM *k* with TM *k*, which is an architectural claim
    about the two proteins rather than a measurement. Checking it against a
    real global alignment is cheap and it is the difference between a
    comparison and an assumption.
    """
    result = tm_index_correspondence()
    assert result["n_helices"] == 38
    assert result["n_agree"] >= 36, result["disagree"]
    assert 0.35 < result["identity"] < 0.65, (
        "paralogue identity should be around a half; a very different number "
        "means the wrong sequences are being aligned")
    for row in result["disagree"]:
        assert row["tm"] in (29,), f"a new helix now disagrees: {row}"


def test_the_helix_pairing_does_not_depend_on_the_boundary_slack():
    """A conclusion resting on a tolerance is a conclusion about the tolerance.

    37 of 38 from zero slack to eight residues of it, so the default of five
    is not doing any work. Above that the count falls — and the reason is
    worth recording because it is not disagreement: adjacent PIEZO2 helices
    sit as little as two residues apart, so a wide window matches two helices
    at once and the *unique* match the criterion asks for stops existing.
    """
    counts, ambiguous = [], []
    try:
        for slack in (0.0, 2.0, 5.0, 8.0):
            PARAMETERS.set_value("paralogue.tm_boundary_slack", slack)
            counts.append(tm_index_correspondence()["n_agree"])
        PARAMETERS.set_value("paralogue.tm_boundary_slack", 15.0)
        wide = tm_index_correspondence()
        ambiguous = [r for r in wide["disagree"] if len(r["landed_in"]) > 1]
    finally:
        PARAMETERS.reset("paralogue.tm_boundary_slack")

    assert len(set(counts)) == 1, f"the answer moved with the slack: {counts}"
    assert wide["n_agree"] < counts[0]
    assert len(ambiguous) > 5, (
        "a wide window should fail by matching two helices at once; if it "
        "fails some other way the criterion is not what this claims")


def test_no_constant_offset_relates_the_two_proteins():
    """The roadmap's standing requirement, measured rather than asserted."""
    numbering = paralogue_map()
    offsets = {numbering.offset_at(r, "mouse") for r in range(600, 2400, 50)}
    offsets.discard(None)
    assert len(offsets) > 5, (
        f"only {len(offsets)} distinct offsets across the chain; if a constant "
        f"worked the alignment would not be needed")


# ------------------------------------- the finding: coverage, not shape

def test_the_naive_dome_difference_is_a_coverage_artefact():
    """Round 83's result, and the reason the naive comparison is kept.

    Measured directly the two domes look very different. Restricted to the
    transmembrane helices both entries resolve, PIEZO2 falls into the PIEZO1
    range. Both rows are reported because the gap between them is the finding.
    """
    first, second = _require("7WLT"), _require("6KG7")
    result = dome_comparison(first, second, "7WLT", "6KG7")

    assert result["n_helices_piezo2"] == 38, (
        "6KG7 resolves every transmembrane helix, which is what makes the "
        "naive comparison unfair and the matched one possible")
    assert result["n_helices_piezo1"] < result["n_helices_piezo2"]
    assert len(result["shared_helices"]) == result["n_helices_piezo1"]

    naive_1, naive_2 = result["naive"]
    matched_1, matched_2 = result["coverage_matched"]

    # Naive: PIEZO2's dome looks half again as deep.
    assert naive_2.depth_nm > 1.5 * naive_1.depth_nm
    # Matched: the difference mostly goes away, and the excess area reverses.
    assert matched_2.depth_nm < 1.3 * matched_1.depth_nm
    assert matched_2.excess_nm2 < matched_1.excess_nm2
    assert abs(matched_2.radius_nm - matched_1.radius_nm) < 1.0

    # The PIEZO1 row must not move when it is restricted to helices it already
    # supplied all of; if it does, the matching is doing something else.
    assert matched_1.depth_nm == pytest.approx(naive_1.depth_nm, rel=1e-9)


def test_the_matched_radius_of_curvature_is_inside_the_piezo1_range():
    """A result about generality: on the dome the two are indistinguishable."""
    piezo2 = _require("6KG7")
    radii = []
    for pdb in ("7WLT", "6B3R"):
        result = dome_comparison(_require(pdb), piezo2, pdb, "6KG7")
        one, two = result["coverage_matched"]
        radii.append((one.radius_nm, two.radius_nm))
    piezo1_radii = [a for a, _ in radii]
    for _, matched in radii:
        assert min(piezo1_radii) - 1.0 < matched < max(piezo1_radii) + 1.0


# ------------------------------------- the finding: the mode is the fold's

def test_the_gating_mode_is_present_in_piezo2():
    """The headline, with the control that makes it mean something.

    PIEZO1's lowest symmetric mode overlaps a single PIEZO2 symmetric mode
    strongly, and almost all of it lies inside PIEZO2's symmetric subspace. The
    shuffled-correspondence control says that is not what any two elastic
    networks of this size give each other.
    """
    result = mode_comparison(_require("7WLT"), _require("6KG7"))

    assert result.n_sites > 1000
    assert result.superposition_rmsd < 8.0, (
        "the two folds should superpose; a large RMSD means the site "
        "correspondence is wrong rather than that the proteins differ")
    assert result.best_overlap > 0.6
    assert result.best_symmetry == "A", (
        "the overlap must land on a symmetric mode — an E mode cannot couple "
        "to isotropic tension, so it could not be a gating coordinate")
    assert result.symmetric_subspace_overlap > 0.85
    assert result.beats_control
    assert result.shuffled_control < 0.4


def test_the_protomer_order_is_searched_rather_than_read_off_the_labels():
    """It is not the identity here, which is exactly why it is searched.

    This project has twice been given a wrong answer by trusting deposited
    chain labels for rotational order. Across a paralogue there is no reason
    at all for them to agree, and they do not.
    """
    result = mode_comparison(_require("7WLT"), _require("6KG7"))
    assert result.protomer_order != (0, 1, 2), (
        "if the order has become the identity, this test no longer "
        "demonstrates anything and the claim needs re-checking")
    assert sorted(result.protomer_order) == [0, 1, 2]


def test_both_networks_put_a_symmetric_mode_low_in_the_spectrum():
    """The weaker, coverage-free version of the same statement."""
    result = mode_comparison(_require("7WLT"), _require("6KG7"))
    assert "A" in result.piezo1_symmetry and "A" in result.piezo2_symmetry
    assert result.piezo1_symmetry.index("A") <= 4
    assert result.piezo2_symmetry.index("A") <= 4
    assert len(result.piezo1_symmetry) == len(result.piezo2_symmetry)


# ------------------------------------------------- the registry was wrong

def test_the_registry_note_matches_what_the_file_contains():
    """6KG7's note said "resolves residues 8-823". It resolves 8-2822.

    Written as a test rather than only fixed, because the note was wrong for
    long enough that the entry was treated as a blade fragment and excluded
    from a comparison it was the best available structure for.
    """
    from piezo1.io.registry import load_registry

    record = load_registry().get("6KG7")
    if record is None or not record.available:
        pytest.skip("6KG7 not downloaded")
    structure = Structure.from_file(record.path)
    mask = structure.mask_ca() & (structure.chain == structure.chains[0])
    residues = np.sort(structure.res_seq[mask])

    assert len(residues) == 1817
    assert int(residues.min()) == 8 and int(residues.max()) == 2822
    assert "8-823" not in record.note
    assert "2,822" in record.note or "2822" in record.note
    # ...and it really is more than any PIEZO1 entry resolves.
    # ...and it resolves more than any PIEZO1 entry, which is the claim the
    # note makes. The invertebrates are excluded: 9UOY resolves more still,
    # and it is not a PIEZO1.
    for other in load_registry():
        # Experimental PIEZO1 entries only. The AlphaFold monomers cover the
        # whole chain by construction — they are a prediction, not a thing
        # anyone resolved — and the invertebrates are not PIEZO1.
        if (other.pdb == "6KG7" or not other.available
                or other.state in ("fragment", "predicted")
                or other.protein != "PIEZO1"):
            continue
        st = Structure.from_file(other.path)
        m = st.mask_ca() & (st.chain == st.chains[0])
        assert int(m.sum()) <= len(residues), other.pdb


def test_the_registry_contains_only_piezo_structures():
    """The catalogue must not sweep up whatever else has been downloaded.

    `build_structure_registry` globs every `.cif` in the structure directory,
    which was harmless until Round 31 downloaded the HaloTag crystal structure
    for the fusion geometry. The next rebuild — Round 83's — swept 6U32 in as
    an unclassified entry of "unknown" species, and it took two unrelated tests
    failing to notice: the ligand audit found the tag's dye and reported it as
    a possible undocumented modulator.

    Every entry must therefore be identifiable as a PIEZO, which is a
    measurement the paralogue module already makes.
    """
    from piezo1.io.registry import load_registry

    from piezo1.core.numbering_check import PROTEIN_NAMES

    for record in load_registry():
        # Read from the naming table rather than a list copied into the test:
        # the family grew from four proteins to five when Arabidopsis PIEZO
        # was added, and a copied list makes that a test failure instead of a
        # fact about the catalogue.
        assert record.protein in set(PROTEIN_NAMES.values()), record.pdb
        if not record.available:
            continue
        identity = identify_numbering(Structure.from_file(record.path))
        assert identity.explained, f"{record.pdb}: {identity.summary()}"
        assert record.protein == identity.protein, (
            f"{record.pdb}: registry says {record.protein}, coordinates say "
            f"{identity.protein}")


def test_unassigned_residues_are_not_counted_as_disagreement():
    """3JAC carries 346 UNK in 918 C-alphas, and matches at 1.000 on the rest.

    ``AA3TO1`` maps UNK to X, so membership is not the test. Counting them as
    mismatches put 3JAC at 0.623 — indistinguishable from a numbering error,
    and it was read as one until the mismatches were looked at and every single
    one turned out to be an unknown.
    """
    identity = identify_numbering(_require("3JAC"))
    assert identity.n_unassigned > 300
    assert identity.identity == pytest.approx(1.0)
    assert identity.confident
    assert "unassigned" in identity.summary()


def test_four_human_entries_carry_a_block_numbered_22_low():
    """A localised register error that the whole-file identity hides.

    8ZU3, 8YFC, 9VMX and 8YFG score 0.932 — past any sensible floor — and the
    7% they are missing is not spread out. It is residues 767-857, every one
    disagreeing, every one agreeing again if read 22 higher. 8YEZ resolves the
    same region and does not have it, so this is a property of those
    depositions rather than of the sequence.
    """
    affected = []
    for pdb in ("8ZU3", "8YFC", "9VMX", "8YFG"):
        identity = identify_numbering(_require(pdb))
        assert identity.confident, "it is still confidently the right protein"
        assert not identity.clean, f"{pdb} should report its block"
        assert identity.blocks
        for block in identity.blocks:
            assert block.repaired_by == 22
            assert block.repaired_identity == pytest.approx(1.0)
            assert 760 <= block.start and block.end <= 860
        affected.append(sum(b.n_residues for b in identity.blocks))
    assert all(n > 70 for n in affected), affected

    # The entry that resolves the same region without the fault.
    clean = identify_numbering(_require("8YEZ"))
    assert clean.clean and clean.identity == pytest.approx(1.0)


def test_a_point_difference_is_not_reported_as_a_numbering_fault():
    """An engineered variant changes one residue and must not be flagged."""
    import json

    st = _require("7WLT")
    mask = st.mask_ca() & (st.chain == st.chains[0])
    sequence = json.loads(
        (RESOURCE_DIR / "uniprot_mouse.json").read_text())["sequence"]
    names = st.res_name[mask].copy()
    names[5] = "TRP"                      # a single substitution, wherever it lands
    assert mismatch_blocks(st.res_seq[mask], names, sequence) == []


def test_6lqi_is_numbered_in_the_splice_isoforms_own_coordinates():
    """A defect the identification found, and a live one.

    6LQI is the Piezo1.1 isoform, which lacks residues 1382-1405. It is
    deposited in the isoform's **own continuous numbering**, so every residue
    past the splice site is 24 lower than its canonical counterpart. Agreement
    with canonical mouse Piezo1 is 1.000 before the site and 0.058 after it —
    and 1.000 again once shifted by +24.

    That is more than half the chain, and it means any annotation this project
    applies to 6LQI by residue number — a transmembrane helix, a domain
    boundary, a variant — is in the wrong place there. The comparison refuses
    the entry rather than measuring it; fixing the rest is its own round.
    """
    identity = identify_numbering(_require("6LQI"))
    assert not identity.confident
    assert identity.explained and identity.splice is not None
    shift = identity.splice
    assert shift.offset == 24
    assert shift.identity_before > 0.99 and shift.identity_after > 0.99
    assert shift.n_after > 700
    assert 1300 < shift.breakpoint < 1420

    result = compare(piezo1_pdb="6LQI", piezo2_pdb="6KG7")
    assert "error" in result and "canonical numbering" in result["error"]


def test_the_splice_detector_refuses_a_file_that_is_merely_wrong():
    """A partial rescue would let any badly numbered file be explained away."""
    import dataclasses

    from piezo1.core.numbering_check import detect_splice

    st = _require("7WLT")
    mask = st.mask_ca() & (st.chain == st.chains[0])
    numbers, names = st.res_seq[mask], st.res_name[mask]
    sequence = json_sequence("mouse")

    # A genuine splice: canonical, then a constant shift.
    shifted = np.where(numbers > 1500, numbers - 24, numbers)
    found = detect_splice(shifted, names, sequence)
    assert found is not None and found.offset == 24

    # Not a splice: the second half scrambled rather than shifted.
    rng = np.random.default_rng(0)
    scrambled = numbers.copy().astype(int)
    tail = scrambled > 1500
    scrambled[tail] = rng.permutation(scrambled[tail])
    assert detect_splice(scrambled, names, sequence) is None

    # Not a splice: an entry that already matches has nothing to detect.
    assert detect_splice(numbers, names, sequence) is None


def json_sequence(name: str) -> str:
    import json

    return json.loads((RESOURCE_DIR / f"uniprot_{name}.json").read_text())["sequence"]


def test_the_piezo2_annotation_resources_are_committed_and_parallel():
    """Both proteins' transmembrane lists must come from the same source.

    Otherwise the two dome measurements differ by how their membrane surface
    was defined rather than by the shape of the membrane.
    """
    import json

    helices, lengths = {}, {}
    for name in REFERENCES:
        path = RESOURCE_DIR / f"uniprot_{name}.json"
        assert path.exists(), f"{path.name} is not committed"
        data = json.loads(path.read_text())
        helices[name] = data["n_transmembrane"]
        lengths[name] = data["length"]

    # The five vertebrate PIEZOs share the 38-helix architecture this project's
    # domain table is built on. Nothing else does — 36 for the worm, 40 for the
    # fly, 35 for the plant and the amoeba — which is exactly why they are the
    # interesting control and why nothing here may transfer a helix index to
    # them by number.
    assert all(helices[n] == 38 for n in
               ("human", "mouse", "rat", "human_piezo2", "mouse_piezo2"))
    assert all(helices[n] != 38 for n in
               ("worm_piezo", "fly_piezo", "plant_piezo", "dicty_piezo"))
    assert lengths == {"human": 2521, "mouse": 2547, "rat": 2535,
                       "human_piezo2": 2752, "mouse_piezo2": 2822,
                       "worm_piezo": 2442, "fly_piezo": 2551,
                       "plant_piezo": 2462, "dicty_piezo": 3080}
    # Nine lengths, no two the same and no constant offset between any pair.
    # The point of stating it: every one of these is a numbering system, and a
    # residue number quoted without one is not a residue.
    assert len(set(lengths.values())) == len(lengths)


def test_compare_refuses_the_arguments_the_wrong_way_round():
    result = compare(piezo1_pdb="6KG7", piezo2_pdb="7WLT")
    assert "error" in result and "expected" in result["error"]


def test_the_paralogue_parameters_are_registered_with_their_reasons():
    for key in ("paralogue.min_sequence_identity", "paralogue.min_identity_margin",
                "paralogue.tm_boundary_slack"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, key
        assert parameter.citation == "method_choice" and parameter.source_note
