"""Getting a number out, and the three ways the file could lie downstream.

An exported file is read somewhere none of this project's guards can reach, so
the checks are about what the *file* says rather than what the writer meant.

1. **The round trip is checked numerically**, element by element against the
   source array — a file that merely opens proves nothing.
2. **Unscored must not read as zero.** A residue the analysis could not score
   arriving in PyMOL indistinguishable from one that scored zero is a confident
   wrong number in somebody else's session.
3. **The column quantises and has a range.** A value that does not fit must be
   refused rather than truncated, and the loss that *is* accepted has to be
   stated rather than assumed away.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.core.export import (COLUMN_RANGE, COLUMN_STEP, COLUMN_WIDTH,
                                SENTINEL, fits_column, read_scalar_pdb,
                                write_scalar_pdb)


@pytest.fixture(scope="module")
def structure():
    path = STRUCTURE_DIR / "8YEZ.cif"
    if not path.exists():
        pytest.skip("8YEZ.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


@pytest.fixture(scope="module")
def scores(structure):
    """A reproducible per-residue scalar over half the residues."""
    residues = sorted(set(int(v) for v in structure.res_seq))
    rng = np.random.default_rng(11)
    chosen = residues[::2]
    return {r: float(v) for r, v in zip(chosen, rng.random(len(chosen)))}


# ------------------------------------------------- the round trip, numerically

def test_the_written_column_matches_the_source_array_atom_by_atom(
        tmp_path, structure, scores):
    """Not "it opens". Every atom, against the value it was given."""
    report = write_scalar_pdb(structure, scores, tmp_path / "s.pdb", scale=100.0)
    occupancy, b = read_scalar_pdb(report.path)

    assert len(b) == structure.n_atoms
    expected = np.full(structure.n_atoms, np.nan)
    for residue, value in scores.items():
        expected[structure.res_seq == residue] = value * 100.0

    measured = np.isfinite(expected)
    assert np.allclose(b[measured], expected[measured], atol=COLUMN_STEP / 2), (
        "the B-factor column does not carry the values it was given")
    assert measured.sum() == report.n_measured


def test_the_occupancy_column_flags_exactly_the_scored_atoms(
        tmp_path, structure, scores):
    report = write_scalar_pdb(structure, scores, tmp_path / "s.pdb")
    occupancy, b = read_scalar_pdb(report.path)

    scored = np.isin(structure.res_seq, sorted(scores))
    assert np.array_equal(occupancy > 0, scored)
    assert (b[~scored] == SENTINEL).all()


# ------------------------------------- unscored is distinguishable from zero

def test_a_residue_scored_zero_is_not_confused_with_an_unscored_one(
        tmp_path, structure):
    """The trap this module exists to avoid, constructed explicitly.

    One residue scores exactly 0.0 and its neighbour is absent from the
    mapping. In the written file they must differ in **both** columns, or a
    reader in another program cannot tell them apart.
    """
    residues = sorted(set(int(v) for v in structure.res_seq))
    zero, absent = residues[0], residues[1]
    report = write_scalar_pdb(structure, {zero: 0.0, residues[2]: 0.5},
                              tmp_path / "s.pdb")
    occupancy, b = read_scalar_pdb(report.path)

    at_zero = structure.res_seq == zero
    at_absent = structure.res_seq == absent
    assert (b[at_zero] == 0.0).all() and (occupancy[at_zero] == 1.0).all()
    assert (b[at_absent] == SENTINEL).all() and (occupancy[at_absent] == 0.0).all()
    assert b[at_zero][0] != b[at_absent][0], (
        "a genuine zero and an unscored residue write the same B-factor")


def test_the_summary_states_the_unscored_count_and_how_to_select_them(
        tmp_path, structure, scores):
    report = write_scalar_pdb(structure, scores, tmp_path / "s.pdb")
    text = report.summary()
    assert "unscored" in text
    assert "occupancy 0.00" in text
    assert "cannot be read as a score of zero" in text


def test_the_header_says_the_column_is_not_a_b_factor(tmp_path, structure,
                                                      scores):
    """Someone opening this in six months has only the file to go on."""
    report = write_scalar_pdb(structure, scores, tmp_path / "s.pdb",
                              scale=100.0, name="conservation")
    header = report.path.read_text().splitlines()[0]
    assert "NOT A B-FACTOR" in header
    assert "CONSERVATION" in header
    assert "scale x100" in report.path.read_text()


# --------------------------------------------- the format's limits are stated

def test_the_quantisation_is_reported_rather_than_assumed_away(
        tmp_path, structure, scores):
    """A 0-1 score written raw has about a hundred levels; scaled, ten thousand.

    Both are the truth about the file, and the caller has to be able to see
    which one they got.
    """
    raw = write_scalar_pdb(structure, scores, tmp_path / "raw.pdb")
    scaled = write_scalar_pdb(structure, scores, tmp_path / "scaled.pdb",
                              scale=100.0)
    assert raw.levels < 200
    assert scaled.levels > 5000
    assert "distinguishable levels" in raw.summary()


def test_a_value_the_column_cannot_hold_is_refused_not_truncated(
        tmp_path, structure, scores):
    with pytest.raises(ValueError, match="cannot hold"):
        write_scalar_pdb(structure, scores, tmp_path / "s.pdb", scale=1e6)


def test_a_value_colliding_with_the_sentinel_is_refused(tmp_path, structure):
    """Otherwise the flag and the data become the same number."""
    residues = sorted(set(int(v) for v in structure.res_seq))
    with pytest.raises(ValueError, match="sentinel"):
        write_scalar_pdb(structure, {residues[0]: -99.99, residues[1]: 0.5},
                         tmp_path / "s.pdb")


def test_an_empty_mapping_is_refused(tmp_path, structure):
    with pytest.raises(ValueError, match="no atom carries a value"):
        write_scalar_pdb(structure, {}, tmp_path / "s.pdb")


def test_a_per_atom_array_of_the_wrong_length_is_refused(tmp_path, structure):
    with pytest.raises(ValueError, match="one value per atom"):
        write_scalar_pdb(structure, np.zeros(7), tmp_path / "s.pdb")


# ----------------------------------------------- it works on a real analysis

def test_a_real_computed_scalar_survives_the_trip(tmp_path, structure):
    """End to end on something the project actually computes.

    The wetting prediction's per-residue energies are a scalar a user would
    plausibly want to colour by in ChimeraX, and they are exactly the kind that
    covers only part of the chain — which is why the unscored flag matters.
    """
    from piezo1.analysis.hydration import load_grid, predict_wetting
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    grid = load_grid()
    if not grid.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")

    blocks, _ = protomer_blocks(structure)
    profile = pore_profile(structure, detect_c3_axis(blocks))
    prediction = predict_wetting(structure, profile, grid=grid)
    scores = {p.residue: p.energy for p in prediction.points}
    assert scores, "the wetting prediction scored nothing to export"

    report = write_scalar_pdb(structure, scores, tmp_path / "wetting.pdb",
                              scale=10.0, name="wetting energy")
    occupancy, b = read_scalar_pdb(report.path)

    assert 0 < report.n_measured < report.n_atoms, (
        "a partial scalar is the interesting case; this one covered everything")
    assert report.n_unmeasured == int((occupancy == 0).sum())
    for residue, value in list(scores.items())[:20]:
        here = structure.res_seq == residue
        assert np.allclose(b[here], value * 10.0, atol=COLUMN_STEP / 2)


def test_the_file_is_readable_by_our_own_parser(tmp_path, structure, scores):
    """It has to be a PDB file, not merely a file with PDB-shaped lines."""
    from piezo1.io.cif_reader import read_pdb_atoms

    report = write_scalar_pdb(structure, scores, tmp_path / "s.pdb", scale=100.0)
    parsed = read_pdb_atoms(report.path)
    assert len(parsed["xyz"]) == structure.n_atoms
    assert np.allclose(parsed["xyz"], structure.xyz, atol=5e-4)


# ------------------------------------------- the field width, derived not recalled

def test_the_column_limits_are_what_actually_fits_six_characters():
    """The constant was wrong on the first attempt, so it is now measured.

    ``-999.99`` and ``9999.99`` are both seven characters. Writing either
    overflows the B-factor field and shifts every column after it, producing a
    file that still parses and is wrong.
    """
    assert fits_column(COLUMN_RANGE[0]) and fits_column(COLUMN_RANGE[1])
    assert not fits_column(-999.99), "seven characters; it does not fit"
    assert not fits_column(9999.99), "seven characters; it does not fit"
    assert len(f"{COLUMN_RANGE[1]:{COLUMN_WIDTH}.2f}") == COLUMN_WIDTH


def test_every_written_line_keeps_the_pdb_column_layout(tmp_path, structure,
                                                        scores):
    """The check that would have caught the overflow: read the element back.

    An overflowing B-factor pushes the element symbol out of columns 77-78, so
    a file with a shifted column has a blank or wrong element there.
    """
    report = write_scalar_pdb(structure, scores, tmp_path / "s.pdb", scale=100.0)
    elements = set()
    for line in report.path.read_text().splitlines():
        if line.startswith(("ATOM  ", "HETATM")):
            assert len(line) >= 78, f"line too short: {line!r}"
            elements.add(line[76:78].strip())
    assert elements and "" not in elements, (
        f"a column has shifted; element field reads {elements}")
