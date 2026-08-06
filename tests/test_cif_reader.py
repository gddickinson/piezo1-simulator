"""Parser correctness.

The regression that motivates most of this file: the original tokenizer treated
only space and tab as whitespace, so the trailing newline of each _atom_site row
became an extra token and silently shifted every subsequent row by one column.
It surfaced far downstream as an int() failure. Column alignment is therefore
asserted explicitly.
"""

import numpy as np
import pytest

from piezo1.core import Structure
from piezo1.io.cif_reader import _tokenize, read_cif_atoms


def test_tokenize_ignores_trailing_newline():
    row = "ATOM   1     N N   . GLU A 1 570  ? 283.282 378.167 269.483 1.00 77.66 ?  570  GLU A N   1 \n"
    assert len(_tokenize(row)) == 21
    assert _tokenize(row)[0] == "ATOM"
    assert _tokenize(row)[-1] == "1"


def test_tokenize_handles_quotes():
    assert _tokenize("A 'two words' B") == ["A", "two words", "B"]
    assert _tokenize('X "O5\'" Y') == ["X", "O5'", "Y"]


def test_tokenize_unquoted_apostrophe_inside_token():
    # An atom name like C1' is legal unquoted because it does not *start* with
    # a quote character.
    assert _tokenize("C1' N") == ["C1'", "N"]


def test_columns_stay_aligned(human_structure):
    st = human_structure
    # Every element must be a plausible chemical symbol, which it would not be
    # if columns had shifted.
    assert set(np.unique(st.element)) <= {
        "C", "N", "O", "S", "P", "SE", "H", "F", "CL", "NA", "K", "MG", "CA",
        "ZN", "FE", "MN", "CU", "GD", "BR", "I"}
    # Residue numbers must be sane.
    prot = st.res_seq[~st.hetero]
    assert prot.min() >= 1 and prot.max() <= 3000


def test_reader_finds_expected_content(human_structure):
    st = human_structure
    assert st.n_atoms > 30000
    assert len(st.chains) >= 3
    assert st.mask_ca().sum() > 3000
    assert st.b_factor.min() >= 0.0


def test_residue_index_consistency(human_structure):
    st = human_structure
    assert len(st.res_first) == st.n_residues
    assert st.res_atom_index.max() == st.n_residues - 1
    # res_first must point at genuine residue starts
    assert (st.res_seq[st.res_first] == st.residue_seq).all()


def test_subset_and_copy_preserve_shape(human_structure):
    st = human_structure
    ca = st.subset(st.mask_ca())
    assert ca.n_atoms == st.mask_ca().sum()
    assert ca.n_residues == ca.n_atoms
    moved = st.copy_with_coords(st.xyz + 1.0)
    # Coordinates are float32, so a +1 A shift on a ~300 A coordinate has
    # about 3e-5 A of representation error. Compare with an absolute tolerance.
    assert np.allclose(moved.xyz - st.xyz, 1.0, atol=1e-3)
    assert moved.n_atoms == st.n_atoms


def test_one_letter_sequence(human_structure):
    seq = human_structure.one_letter_sequence(chain="A")
    assert len(seq) > 1000
    assert set(seq) <= set("ACDEFGHIKLMNPQRSTVWYX")
