"""The family comparison, with every instrument calibrated before it is believed.

Three things here are checking instruments in this project's sense — they exist
to judge whether another number may be trusted — and each is driven against a
case whose answer is known without any of the code being tested:

* the **shuffled null** behind every identity and local score;
* the **alignment window**, which decides whether a mapped position may be read;
* the **numbering identification**, which grew from six references to nine and
  had to be shown not to have lost its margin.

The last of those is the one that would have been easiest to skip. Rat Piezo1 is
94.2% identical to mouse Piezo1 as a sequence, so adding it to ``REFERENCES``
looks certain to collapse the margin ``NumberingIdentity.confident`` requires
and make eleven catalogue entries unreadable. It does not, and the reason is
worth pinning: the identification reads residue names at their own *numbers*
rather than aligning, and the twelve-residue length difference puts every
position past the first indel out of register. Measured, not assumed.
"""

from __future__ import annotations

import random

import numpy as np
import pytest

from piezo1.analysis import homology as H
from piezo1.analysis import homology_sites as HS
from piezo1.analysis.alignment_windows import (alignment_windows,
                                               window_identity, window_score)
from piezo1.core.numbering_check import (PIEZO1_REFERENCES, PROTEIN_NAMES,
                                         REFERENCES, identify_numbering,
                                         reference_entry)
from piezo1.parameters import PARAMETERS

RNG = random.Random(20260812)


def _random_protein(n: int) -> str:
    return "".join(RNG.choice("ACDEFGHIKLMNPQRSTVWY") for _ in range(n))


def _shuffle(sequence: str) -> str:
    letters = list(sequence)
    RNG.shuffle(letters)
    return "".join(letters)


# --------------------------------------------------------------------------
# The family itself
# --------------------------------------------------------------------------

def test_the_family_is_nine_and_every_member_is_committed():
    members = H.family()
    assert len(members) == 9
    assert {m.key for m in members} == set(REFERENCES)
    for m in members:
        assert m.length > 2000 and m.accession and m.organism != "?"
        assert m.protein == PROTEIN_NAMES[m.key]


def test_no_two_members_share_a_length_or_a_helix_count_architecture():
    """The reason nothing may transfer a residue number or a helix index.

    Nine lengths, all different — so a residue number without its numbering
    system is not a residue. And four distinct helix counts, so 'TM12' is not
    a portable name either: only the five vertebrate members have 38.
    """
    members = H.family()
    assert len({m.length for m in members}) == 9

    vertebrate = {m.key for m in members if m.group in ("PIEZO1", "PIEZO2")}
    assert all(m.n_transmembrane == 38 for m in members if m.key in vertebrate)
    assert all(m.n_transmembrane != 38 for m in members
               if m.key not in vertebrate)
    assert {m.n_transmembrane for m in members} == {35, 36, 38, 40}


def test_a_missing_resource_raises_rather_than_shortening_the_family(monkeypatch):
    """A family quietly missing a member reports a narrower range of divergence.

    That would look like a tighter result rather than a broken one, which is
    why the loader raises.
    """
    monkeypatch.setattr(H, "_FAMILY", None)
    monkeypatch.setattr(H, "REFERENCES", REFERENCES + ("nonexistent_piezo",))
    monkeypatch.setattr(H, "GROUPS", dict(H.GROUPS,
                                          invented=("nonexistent_piezo",)))
    with pytest.raises(FileNotFoundError):
        H.family()
    H._FAMILY = None


# --------------------------------------------------------------------------
# Calibration: the shuffled null
# --------------------------------------------------------------------------

def test_the_null_is_calibrated_on_two_known_answers():
    """A sequence against itself, and against something unrelated.

    Both answers are known without running any of the code under test: a
    sequence is 100% identical to itself, and two independently drawn random
    proteins are related only by the amino-acid alphabet.
    """
    sequence = _random_protein(600)
    identity, columns, local = H.align_pair(sequence, sequence)
    assert identity == pytest.approx(1.0)
    assert columns == len(sequence)

    other = _random_protein(600)
    unrelated_identity, _columns, unrelated_local = H.align_pair(sequence, other)
    assert unrelated_identity < 0.35
    assert unrelated_local < local / 5


def test_the_null_can_say_no_when_there_is_no_homology():
    """The case that makes the null a measurement rather than a formality.

    Two unrelated random proteins must come back *not* distinguishable from
    the shuffled null. A null that only ever said yes would pass every real
    pair and assert nothing.
    """
    a, b = _random_protein(700), _random_protein(700)
    null_identity, null_local = H.shuffled_null(a, b, replicates=6, seed=1)
    _identity, _columns, local = H.align_pair(a, b)
    z = (local - null_local.mean()) / null_local.std(ddof=1)
    assert z < PARAMETERS.value("homology.min_z"), (
        "two unrelated random sequences were called homologous")


def test_the_null_preserves_composition():
    """A uniform-alphabet null would flatter every pair in the matrix."""
    sequence = reference_entry("human")["sequence"]
    shuffled = _shuffle(sequence)
    assert sorted(shuffled) == sorted(sequence)
    assert shuffled != sequence


# --------------------------------------------------------------------------
# The measured result
# --------------------------------------------------------------------------

@pytest.mark.parametrize("partner,floor", [("mouse", 0.80), ("rat", 0.80),
                                           ("human_piezo2", 0.45)])
def test_the_close_members_are_identifiable_from_identity_alone(partner, floor):
    relationship = H.relationship("human", partner, replicates=5)
    assert relationship.identity > floor
    assert not relationship.in_twilight_zone
    assert relationship.identity_beats_null and relationship.local_beats_null


@pytest.mark.parametrize("partner", ["plant_piezo", "dicty_piezo"])
def test_the_distant_members_are_homologous_and_their_identity_is_not_evidence(
        partner):
    """The finding the whole module exists to be able to state.

    The local score must be overwhelming and the identity must be inside the
    twilight zone — checked in both directions, because either half alone
    would be a different and weaker claim.
    """
    # Registered default replicate count, because the documented numbers were
    # produced at it and the null's spread depends on it.
    relationship = H.relationship("human", partner)
    assert relationship.in_twilight_zone
    assert relationship.local_beats_null
    assert relationship.local_z > 8 * relationship.identity_z, (
        "the two statistics no longer disagree by a large factor; the "
        "documented argument in docs/HOMOLOGY_SEARCH.md rests on that gap")


def test_at_least_one_pair_has_an_identity_indistinguishable_from_chance():
    """PEZO-1 against Arabidopsis PIEZO: 23.8% real, ~22.5% shuffled.

    Pinned because it is the single most persuasive row in the decision
    document. If it ever stops being true the document needs rewriting, not
    the test relaxing.
    """
    relationship = H.relationship("worm_piezo", "plant_piezo", replicates=10)
    assert relationship.identity < PARAMETERS.value("homology.twilight_identity")
    assert not relationship.identity_beats_null
    assert relationship.local_beats_null
    assert "must not be quoted alone" in relationship.verdict


# --------------------------------------------------------------------------
# Calibration: the alignment window
# --------------------------------------------------------------------------

def test_the_window_is_calibrated_on_a_sequence_against_itself():
    sequence = _random_protein(500)
    assert window_identity(sequence, sequence, 250) == pytest.approx(1.0)
    assert window_score(sequence, sequence, 250) > 3.0


def test_the_window_finds_a_planted_block_that_whole_sequence_identity_misses():
    """The calibration that makes the window worth having.

    A single conserved block inside otherwise scrambled sequence. The window
    centred on the block must find it; the window far from it must not; and
    the whole-sequence identity must be unable to tell the two apart. Without
    this the window is just a smaller version of the number it is supposed to
    improve on.
    """
    source = _random_protein(1200)
    block = slice(500, 620)
    partner = list(_random_protein(1200))
    partner[block] = list(source[block])
    partner = "".join(partner)

    inside = window_score(source, partner, 560)
    outside = window_score(source, partner, 150)
    assert inside > 2.0 and inside > 2 * outside, (inside, outside)
    # Identity sees the block too, but far less sharply — which is the point.
    # The gap between the two statistics is what earns the score its place.
    inside_identity = window_identity(source, partner, 560)
    outside_identity = window_identity(source, partner, 150)
    assert inside_identity > 1.5 * outside_identity, (inside_identity,
                                                      outside_identity)
    assert (inside / max(outside, 1e-9)) > (inside_identity / outside_identity)

    # And the statistic it replaces cannot see it: the block is 10% of the
    # sequence, so whole-sequence identity moves by a few points at most.
    whole, _columns, _local = H.align_pair(source, partner)
    assert whole < 0.4


def test_the_window_refuses_a_composition_matched_shuffle():
    """The instrument must be able to say no, on real PIEZO sequences."""
    human = reference_entry("human")["sequence"]
    windows = alignment_windows("human")
    # A window in the middle of a true self-alignment is far above its null.
    _mapped, score, _identity = windows.at(2447)
    assert (score - windows.null_mean) / windows.null_sd > 10

    # The null itself is built from a shuffle and must sit near zero z.
    assert windows.null_sd > 0
    assert abs(windows.null_mean) < abs(score)
    assert len(human) == 2521


# --------------------------------------------------------------------------
# Curated sites across the family
# --------------------------------------------------------------------------

def test_human_against_itself_is_the_positive_control():
    """Every curated position must come back identical, or nothing else counts."""
    result = HS.report(targets=["human"])
    for row in result.rows:
        assert row.n_unreliable == 0, row.summary()
        assert row.identity == 1.0, row.summary()


def test_an_unreadable_group_reports_none_rather_than_zero_conserved():
    """A group whose every position is unreliable must not read as 0% conserved.

    Perfect non-conservation is the most confident possible way of saying
    nothing, and it is what a naive implementation reports.
    """
    result = HS.report(targets=["dicty_piezo"], groups=["cap_gate_loops"])
    row = result.rows[0]
    assert row.n_readable == 0
    assert row.identity is None
    assert "can be trusted in" in row.summary()
    assert "none of 4 positions" in row.summary()


def test_the_gate_is_readable_out_to_the_plant_and_refused_in_the_amoeba():
    """The measured boundary, and the reason the window width is 101.

    Pinned in both directions: if the worm/fly/plant stop being readable the
    width has lost its power again, and if Dictyostelium starts being readable
    the threshold has become one that cannot refuse.
    """
    result = HS.report(groups=["hydrophobic_gate"])
    readable = {r.target: r.n_readable for r in result.rows}
    assert all(readable[k] == 3 for k in
               ("human", "mouse", "rat", "worm_piezo", "fly_piezo",
                "plant_piezo"))
    assert readable["dicty_piezo"] == 0


def test_the_anchor_brake_is_the_only_universally_identical_group():
    """Human P2113/F2114, identical in every member where it can be read.

    Including Dictyostelium, which refuses the pore gate. Recorded because it
    is a result rather than a formality: the one curated site this project
    holds that has survived the whole eukaryotic tree is not in the pore.
    """
    result = HS.report()
    assert result.universal() == ["anchor_brake"]


def test_the_cap_does_not_travel_outside_the_vertebrates():
    """A limit on this project's own annotation, not a fact about the cap."""
    result = HS.report(groups=["cap_gate", "cap_gate_loops",
                               "cap_constriction"])
    for key in ("worm_piezo", "fly_piezo", "plant_piezo", "dicty_piezo"):
        rows = [r for r in result.rows if r.target == key]
        assert all(r.n_readable <= 1 for r in rows), \
            [r.summary() for r in rows]


# --------------------------------------------------------------------------
# The numbering references, widened from six to nine
# --------------------------------------------------------------------------

def test_adding_rat_does_not_collapse_the_margin_on_a_mouse_entry(structure_7wlt):
    """94.2% identical as a sequence, 0.066 as a numbering.

    The measurement that made it safe to widen ``REFERENCES``. If this ever
    starts scoring high, ``confident`` fails on every mouse entry and
    ``paralogue.compare`` and ``pore_regions`` refuse them all.
    """
    identity = identify_numbering(structure_7wlt)
    assert identity.reference == "mouse"
    assert identity.identity == pytest.approx(1.0, abs=1e-3)
    assert identity.scores["rat"] < 0.25
    assert identity.margin > PARAMETERS.value("paralogue.min_identity_margin")

    # And the sequences really are nearly identical, so the test is about the
    # instrument rather than about the two proteins being distant.
    relationship = H.relationship("mouse", "rat", replicates=3)
    assert relationship.identity > 0.90


def test_every_downloaded_entry_still_identifies_against_nine_references():
    from piezo1.core.structure import Structure
    from piezo1.io.registry import load_registry

    for record in load_registry():
        if not record.available:
            continue
        identity = identify_numbering(Structure.from_file(record.path))
        assert identity.explained, f"{record.pdb}: {identity.summary()}"
        assert identity.reference in REFERENCES


def test_the_piezo1_references_are_the_three_mammals():
    assert set(PIEZO1_REFERENCES) == {"human", "mouse", "rat"}
    assert all(PROTEIN_NAMES[k] == "PIEZO1" for k in PIEZO1_REFERENCES)


# --------------------------------------------------------------------------
# The parameters this all rests on
# --------------------------------------------------------------------------

def test_every_homology_parameter_is_registered_with_a_reason():
    keys = [k for k in PARAMETERS.parameters if k.startswith("homology.")]
    assert len(keys) == 6, keys
    for key in keys:
        parameter = PARAMETERS.get(key)
        assert parameter is not None and parameter.description
        assert parameter.source_note, key


def test_the_twilight_threshold_is_the_published_one():
    assert PARAMETERS.value("homology.twilight_identity") == pytest.approx(0.30)
    assert PARAMETERS.get("homology.twilight_identity").citation == "rost1999"


def test_the_matrix_is_memoised_but_not_cached_to_disk():
    """A stored identity matrix would be a second place for the family to be wrong."""
    from piezo1.config import RESOURCE_DIR

    assert not list(RESOURCE_DIR.glob("*homology*"))
    assert not list(RESOURCE_DIR.glob("*identity_matrix*"))
    first = H.family_matrix(keys=["human", "mouse"], replicates=3)
    second = H.family_matrix(keys=["human", "mouse"], replicates=3)
    assert first is second


def test_window_width_was_chosen_by_a_power_scan_not_by_taste():
    """The narrow window must genuinely have failed, or the width is unmotivated.

    Drives the instrument at the width the first version used and asserts it
    refuses the fly gate — the mapping that is visibly right at offset 0.
    """
    original = PARAMETERS.value("homology.site_window")
    try:
        PARAMETERS.set_value("homology.site_window", 31)
        HS_windows = alignment_windows.__wrapped__ if hasattr(
            alignment_windows, "__wrapped__") else alignment_windows
        import piezo1.analysis.alignment_windows as AW

        AW._WINDOWS.clear()
        windows = HS_windows("fly_piezo")
        _mapped, score, _identity = windows.at(2454)
        narrow_z = (score - windows.null_mean) / windows.null_sd
        assert narrow_z < PARAMETERS.value("homology.min_z")

        PARAMETERS.set_value("homology.site_window", original)
        AW._WINDOWS.clear()
        windows = HS_windows("fly_piezo")
        _mapped, score, _identity = windows.at(2454)
        wide_z = (score - windows.null_mean) / windows.null_sd
        assert wide_z >= PARAMETERS.value("homology.min_z")
        assert wide_z > narrow_z
    finally:
        PARAMETERS.reset()
        import piezo1.analysis.alignment_windows as AW

        AW._WINDOWS.clear()



def _registry_path(pdb: str):
    from piezo1.io.registry import load_registry

    record = load_registry().get(pdb)
    if record is None or not record.available:
        pytest.skip(f"{pdb} not downloaded — run python -m piezo1.io.fetch")
    return record.path


@pytest.fixture
def structure_7wlt():
    from piezo1.core.structure import Structure

    return Structure.from_file(_registry_path("7WLT"))
