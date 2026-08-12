"""Building a trimer from one protomer, calibrated before it is believed.

The known-answer case is exact and available: pull one chain **out of** a real
trimer, hand it back as a monomer, and rebuild it against the trimer it came
from. The answer must be that trimer — 0.00 A over every corresponding residue,
and the same clash count the deposited file itself scores.

That case can also fail, which is what makes it a calibration rather than a
formality: rebuilt against a *different* template the same monomer must come
back measurably worse. Without that half, "0.00 A" would be equally consistent
with a function that returns its input.

The clash counter is calibrated the same way, and it has to be: the assemblies
this feature exists for score thousands, and that number means nothing until a
real trimer's score is known. Deposited trimers give 3-8.
"""

from __future__ import annotations

import dataclasses

import pytest

from piezo1.core.structure import Structure
from piezo1.io.registry import load_registry
from piezo1.parameters import PARAMETERS
from piezo1.structure.assembly import (assemble_trimer, best_template,
                                       is_monomer)
from piezo1.structure.clashes import count_clashes
from piezo1.structure.protomers import well_resolved_chains


def _load(pdb: str) -> Structure:
    record = load_registry().get(pdb)
    if record is None or not record.available:
        pytest.skip(f"{pdb} not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(record.path)


def _extract_chain(structure: Structure, chain: str) -> Structure:
    keep = structure.chain == chain
    fields = ("xyz", "element", "atom_name", "res_name", "res_seq", "chain",
              "hetero", "b_factor", "occupancy", "alt_loc", "entity")
    return dataclasses.replace(
        structure, name=f"{structure.name}_chain{chain}",
        res_first=None, res_atom_index=None,
        **{f: getattr(structure, f)[keep] for f in fields})


# --------------------------------------------------------------------------
# Calibration
# --------------------------------------------------------------------------

def test_a_protomer_rebuilt_on_its_own_trimer_reproduces_it_exactly():
    """The known answer, and it is exact rather than approximate."""
    trimer = _load("6B3R")
    monomer = _extract_chain(trimer, well_resolved_chains(trimer)[0])
    assert is_monomer(monomer)

    result = assemble_trimer(monomer, template="6B3R")
    assert result.ok, result.refusal
    # Not bit-exact zero, and the residual is a fact about 6B3R rather than
    # about this code: its three chains are refined independently and are
    # near-identical, not identical, so chain A placed onto chain B keeps the
    # template's own ~0.001 A inter-protomer deviation. Anything larger would
    # be the assembly, and 0.5 A below is what a *different* template costs.
    assert result.worst_placement == pytest.approx(0.0, abs=0.01)
    assert result.worst_full == pytest.approx(0.0, abs=0.01)
    assert result.n_core == result.n_corresponding
    # And it lands on the deposited trimer's own clash count, which is the
    # part that says the *arrangement* was reproduced and not just the fit.
    assert result.clashes == count_clashes(trimer)


def test_the_same_monomer_on_a_different_template_is_measurably_worse():
    """The half that makes the exact case mean something.

    Without this, 0.00 A would be equally consistent with a function that
    hands back what it was given.
    """
    trimer = _load("6B3R")
    monomer = _extract_chain(trimer, well_resolved_chains(trimer)[0])

    same = assemble_trimer(monomer, template="6B3R")
    other = assemble_trimer(monomer, template="7WLT")
    assert other.ok, other.refusal
    assert other.worst_placement > same.worst_placement + 0.5
    assert other.clashes > same.clashes


def test_the_clash_counter_is_near_zero_on_real_trimers():
    """Otherwise a four-figure count on an assembly says nothing."""
    for pdb in ("6B3R", "7WLT", "9ZIS"):
        record = load_registry().get(pdb)
        if record is None or not record.available:
            continue
        assert count_clashes(Structure.from_file(record.path)) < 50, pdb


# --------------------------------------------------------------------------
# Refusals
# --------------------------------------------------------------------------

def test_a_real_trimer_is_refused_rather_than_replaced_by_a_model_of_itself():
    result = assemble_trimer(_load("7WLT"))
    assert not result.ok
    assert "already models" in result.refusal


def test_a_monomer_with_no_correspondence_is_refused():
    """A construct numbering matches no reference, so nothing can be placed."""
    record = load_registry().get("4PKE")
    if record is None or not record.available:
        pytest.skip("4PKE not downloaded")
    result = assemble_trimer(Structure.from_file(record.path))
    assert not result.ok


# --------------------------------------------------------------------------
# What it produces
# --------------------------------------------------------------------------

def test_the_plant_monomer_becomes_something_three_protomer_code_accepts():
    """The whole point: the dome and the elastic network need three chains."""
    monomer = _load("AF-F4IN58-F1-MODEL_V6")
    result = assemble_trimer(monomer)
    assert result.ok, result.refusal
    assert len(well_resolved_chains(result.structure)) == 3
    assert result.structure.n_atoms == 3 * monomer.n_atoms


def test_the_assembly_says_it_is_a_model_in_every_way_it_can():
    monomer = _load("AF-F4IN58-F1-MODEL_V6")
    result = assemble_trimer(monomer)
    assert result.structure.meta["is_observed"] is False
    assert result.structure.meta["assembly_template"] == result.template
    assert "+trimer(" in result.structure.name
    # Not A/B/C: an assembled file that looks deposited is the thing to avoid.
    assert set(result.structure.chains) == {"X", "Y", "Z"}
    caveat = result.caveat
    assert "MODELLED TRIMER" in caveat
    assert result.template in caveat
    assert "ARRANGEMENT" in caveat


def test_a_template_that_barely_fits_says_so_rather_than_looking_the_same():
    """The plant hits the core floor and the rat does not, and both are drawn.

    Without the flag those two produce the same kind of picture, and one of
    them is a model of a protein while the other is a shape in the right place.
    """
    plant = assemble_trimer(_load("AF-F4IN58-F1-MODEL_V6"))
    assert plant.at_floor, plant.summary()
    assert "AT THE FLOOR" in plant.summary()

    rat = assemble_trimer(_load("AF-Q0KL00-F1-MODEL_V6"))
    assert rat.ok and not rat.at_floor, rat.summary()
    assert rat.n_core > int(PARAMETERS.value("assembly.min_corresponding"))


def test_the_three_placements_agree_with_each_other():
    """A monomer fitting one template protomer well and another badly has been
    forced rather than placed."""
    result = assemble_trimer(_load("AF-Q0KL00-F1-MODEL_V6"))
    assert len(result.placement_rmsd) == 3
    assert result.placement_spread < 1.0, result.placement_rmsd


def test_the_template_is_chosen_by_protein_then_by_coverage():
    monomer = _load("AF-A0A061ACU2-F1-MODEL_V6")      # PEZO-1
    record = load_registry().get(best_template(monomer))
    assert record.protein == "PEZO-1"


def test_most_of_an_assembly_s_shape_is_the_template_s_and_it_says_so():
    """The measurement that says what this feature is and is not for.

    An assembled trimer can be handed to the dome fit, which is the point. The
    number that comes back is 79-96% a measurement of the template, which is
    the limit — and stating it is the difference between a drawable model and
    a claim about a plant protein.
    """
    from piezo1.structure.assembly import borrowed_fraction

    for pdb in ("AF-F4IN58-F1-MODEL_V6", "AF-A0A061ACU2-F1-MODEL_V6"):
        result = assemble_trimer(_load(pdb))
        split = borrowed_fraction(result)
        assert "error" not in split, split
        assert 0.7 < split["borrowed_fraction"] < 1.0, split
        assert split["arrangement_A"] > split["within_protomer_A"]

    # And the caveat carries it, so it cannot be read off the picture alone.
    caveat = assemble_trimer(_load("AF-F4IN58-F1-MODEL_V6")).caveat
    assert "MOSTLY a measurement" in caveat and "79-96%" in caveat
