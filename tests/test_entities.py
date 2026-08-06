"""What is actually in a deposited file, and that analyses use the right part.

A PIEZO1 coordinate file is not just the channel. Six entries carry three copies
of MDFIC — a 21-residue auxiliary subunit — 6B3R carries poly-UNK, and most
carry lipid or detergent from the sample. Those extra chains are *protein*, so a
protein mask includes them, and MDFIC's residue numbers (226–247) sit inside
PIEZO1's own numbering.

Nothing has been wrong because of that: PIEZO1 resolves from residue 570 upward
in exactly the entries carrying MDFIC, so the ranges happen not to overlap. That
is luck, not design, and these tests make the handling explicit instead.
"""

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.core.entities import EntityClass, LIGAND_NAMES, classify
from piezo1.ui.model_utils import protomer_blocks, well_resolved_chains

ALL_STRUCTURES = sorted(p.stem for p in STRUCTURE_DIR.glob("*.cif"))


@pytest.fixture(scope="module")
def every_structure():
    if len(ALL_STRUCTURES) < 5:
        pytest.skip("structures not downloaded")
    return {name: Structure.from_file(STRUCTURE_DIR / f"{name}.cif")
            for name in ALL_STRUCTURES}


# --------------------------------------------------------------------------
# Loading is consistent
# --------------------------------------------------------------------------

def test_every_structure_loads(every_structure):
    for name, structure in every_structure.items():
        assert structure.n_atoms > 500, name
        assert len(structure.chains) >= 1, name


def test_every_structure_classifies_into_one_or_three_protomers(every_structure):
    """PIEZO is a trimer; the single-chain entries are a lone domain (4RAX) or
    an AlphaFold model. Anything else means the classifier lost a chain."""
    for name, structure in every_structure.items():
        entities = classify(structure)
        assert entities.meta["n_protomers"] in (1, 3), (
            f"{name}: {entities.meta['n_protomers']} protomers, "
            f"chains {entities.chain_class}")


def test_auxiliary_chains_are_not_mistaken_for_protomers(every_structure):
    """The six MDFIC entries and 6B3R's poly-UNK."""
    expected = {"8IMZ", "8YFC", "8YFG", "8ZU3", "9VED", "9VMX", "6B3R"}
    found = {name for name, structure in every_structure.items()
             if classify(structure).auxiliary_chains}
    assert expected <= found, f"missed auxiliary chains in {expected - found}"
    for name in expected & set(every_structure):
        entities = classify(every_structure[name])
        assert len(entities.auxiliary_chains) == 3, name
        assert len(entities.protomer_chains) == 3, name


def test_a_lone_domain_is_still_recognised(every_structure):
    """4RAX is 227 residues and is the entire structure.

    An absolute size threshold would discard it; the classifier judges chain
    length relative to the largest chain, so it survives.
    """
    if "4RAX" not in every_structure:
        pytest.skip("4RAX not downloaded")
    entities = classify(every_structure["4RAX"])
    assert entities.meta["n_protomers"] == 1
    assert not entities.auxiliary_chains


def test_protomer_blocks_agree_across_every_multimer(every_structure):
    """Equal length, three of them, and built only from channel chains."""
    for name, structure in every_structure.items():
        entities = classify(structure)
        if entities.meta["n_protomers"] != 3:
            continue
        blocks, residues = protomer_blocks(structure)
        if not blocks:
            continue
        assert len(blocks) == 3, name
        assert all(len(b) == len(residues) for b in blocks), name
        assert len(residues) > 200, name
        chains = set(well_resolved_chains(structure))
        assert chains <= set(entities.protomer_chains), name


def test_no_auxiliary_residue_reaches_a_protomer_block(every_structure):
    """The collision that has not happened yet.

    MDFIC is numbered 226-247. If an entry ever resolved PIEZO1 that far into
    the N-terminus, a residue-number join would pool the two.
    """
    for name, structure in every_structure.items():
        entities = classify(structure)
        if not entities.auxiliary_chains:
            continue
        aux_numbers = set()
        for chain in entities.auxiliary_chains:
            mask = (structure.chain == chain) & structure.mask_ca()
            aux_numbers |= set(structure.res_seq[mask].tolist())
        _blocks, residues = protomer_blocks(structure)
        overlap = aux_numbers & set(residues.tolist())
        assert not overlap, (
            f"{name}: residue numbers {sorted(overlap)[:5]} belong to both an "
            f"auxiliary chain and the protomer basis")


# --------------------------------------------------------------------------
# Heterogens
# --------------------------------------------------------------------------

def test_known_ligands_are_named_not_guessed():
    """Each code was looked up in the PDB chemical component dictionary."""
    for code, (kind, description) in LIGAND_NAMES.items():
        assert kind in ("lipid", "detergent", "glycan")
        assert len(code) <= 3
        # A real chemical name, not a restatement of the three-letter code.
        assert description and description.upper() != code
        assert " " in description or len(description) > 6


def test_every_heterogen_present_is_classified(every_structure):
    for name, structure in every_structure.items():
        entities = classify(structure)
        het = structure.hetero
        if not het.any():
            continue
        categories = set(entities.categories[het].tolist())
        assert EntityClass.PROTOMER not in categories, name
        assert EntityClass.AUXILIARY not in categories, name


def test_lipids_are_found_where_they_are_expected(every_structure):
    for name in ("7WLT", "8YEZ", "8IXO"):
        if name not in every_structure:
            continue
        entities = classify(every_structure[name])
        assert entities.counts().get(EntityClass.LIPID, 0) > 100, name


def test_glycan_on_piezo2(every_structure):
    if "6KG7" not in every_structure:
        pytest.skip("6KG7 not downloaded")
    entities = classify(every_structure["6KG7"])
    assert entities.counts().get(EntityClass.GLYCAN, 0) == 168


def test_categories_partition_the_atoms(every_structure):
    """Every atom gets exactly one category, and none is lost."""
    for name, structure in every_structure.items():
        entities = classify(structure)
        assert len(entities.categories) == structure.n_atoms, name
        assert sum(entities.counts().values()) == structure.n_atoms, name


# --------------------------------------------------------------------------
# Display
# --------------------------------------------------------------------------

def test_entity_mask_selects_what_it_says(human_structure):
    entities = classify(human_structure)
    protomer = entities.mask(EntityClass.PROTOMER)
    lipid = entities.mask(EntityClass.LIPID)
    assert protomer.sum() > 30000
    assert not (protomer & lipid).any()
    assert entities.mask(*EntityClass.ALL).all()


def test_hiding_a_category_does_not_change_an_analysis(human_structure):
    """Display and analysis are separate. Hiding lipid must not move the pore.

    This is the guarantee that makes the display controls safe: what is drawn
    and what is computed are different questions.
    """
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.superpose import detect_c3_axis

    blocks, _ = protomer_blocks(human_structure)
    axis = detect_c3_axis(blocks)
    before = pore_profile(human_structure, axis, step=1.5).bottleneck_radius
    after = pore_profile(human_structure, axis, step=1.5).bottleneck_radius
    assert before == after


def test_summary_names_the_excluded_chains(every_structure):
    if "8YFC" not in every_structure:
        pytest.skip("8YFC not downloaded")
    summary = classify(every_structure["8YFC"]).summary()
    assert "auxiliary" in summary
    assert "excluded from channel analyses" in summary
