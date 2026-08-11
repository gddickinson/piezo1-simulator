"""Full-length models: built for every state, and consumable by every analysis.

The graft existed as coordinates for one protomer, which is the right unit for
asking how well it fits and the wrong unit for everything else — no C3 axis, no
protomer blocks, no lumen. So the distal blade, half the protein by residue
count, was in no measurement this project makes.

These check the two things that decide whether the new model is usable: that it
is a :class:`~piezo1.core.Structure` the existing pipeline cannot tell apart
from a deposited one, and that it can never be *mistaken* for one.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.parameters import PARAMETERS
from piezo1.structure.full_length import (FILL_MODES, build_full_length,
                                          is_full_length, predicted_mask,
                                          resolved_gaps)


def _require(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


# ------------------------------------------------------- it builds, in every state

@pytest.mark.parametrize("pdb", ["8YEZ", "7WLT", "11ZC", "8IXO"])
def test_a_full_length_model_builds_for_every_state(pdb):
    """Curved, flat, flattened and intermediate — the request was all of them."""
    experimental = _require(pdb)
    model = build_full_length(experimental, "full")

    assert model.structure.n_atoms > experimental.n_atoms
    assert len(model.seams) == 3, "one graft per protomer"
    assert model.n_predicted_residues > 1000
    assert is_full_length(model.structure)
    assert model.structure.name.endswith("+AF-full")


def test_the_four_modes_are_nested():
    """Each mode adds what its label says and nothing else."""
    experimental = _require("7WLT")
    counts = {}
    for key, _label, _tip in FILL_MODES:
        model = build_full_length(experimental, key)
        counts[key] = model.n_predicted_residues

    assert counts["none"] == 0
    assert counts["gaps"] > 0 and counts["blade"] > 0
    assert counts["full"] == counts["gaps"] + counts["blade"]


def test_an_unknown_mode_is_refused():
    with pytest.raises(ValueError, match="mode must be"):
        build_full_length(_require("7WLT"), "everything")


# ------------------------------- the pipeline cannot tell it from a deposited file

def test_the_analyses_run_on_it_unchanged():
    """The whole point: no analysis needs to know this model exists.

    `protomer_blocks` locates residues with `searchsorted`, which on unsorted
    input returns the wrong atoms rather than failing — so the built model is
    ordered by chain and residue like a deposited file. This is the test that
    caught that: the first version raised an IndexError, and a slightly
    different residue set would have produced a silently scrambled protomer.
    """
    from piezo1.structure.geometry import measure_dome, tm_surface_points
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    model = build_full_length(_require("7WLT"), "full")
    built = model.structure

    blocks, residues = protomer_blocks(built)
    assert len(blocks) == 3
    assert len(residues) > 2000, "the full-length model has more shared residues"
    for block in blocks:
        assert len(block) == len(residues)

    axis = detect_c3_axis(blocks)
    assert axis.angle_deg == pytest.approx(120.0, abs=1.0)
    assert axis.rmsd < 2.0

    points, resolved = tm_surface_points(built, "mouse")
    assert len(resolved) == 38, (
        "the point of a full-length model is that every transmembrane helix "
        "is present; the deposited entry resolves 22")
    dome = measure_dome(blocks, points)
    assert 5.0 < dome.radius_of_curvature / 10.0 < 20.0


def test_residues_come_out_sorted_within_each_chain():
    """The ordering `searchsorted` depends on, asserted rather than assumed."""
    built = build_full_length(_require("7WLT"), "full").structure
    for chain in built.chains:
        mask = built.mask_ca() & (built.chain == chain)
        numbers = built.res_seq[mask]
        assert np.all(np.diff(numbers) >= 0), f"chain {chain} is not sorted"


def test_the_elastic_network_runs_on_it():
    """Slow but not impossible — and the symmetry labelling still works."""
    from piezo1.physics.anm import ANM
    from piezo1.structure.protomers import protomer_blocks

    blocks, _ = protomer_blocks(build_full_length(_require("7WLT"), "full").structure)
    anm = ANM.from_trimer(blocks).build()
    modes = anm.calc_modes(n_modes=6)
    anm.label_symmetry(modes)
    assert modes.n_modes == 6
    assert set(modes.symmetry) <= {"A", "E"}


# --------------------------------------------------- it cannot be mistaken for data

def test_provenance_is_derived_so_it_cannot_go_stale():
    """A stored mask would survive a subset and describe the wrong atoms."""
    model = build_full_length(_require("7WLT"), "full")
    built = model.structure
    predicted = predicted_mask(built)
    assert predicted.sum() > 0
    assert predicted.sum() < built.n_atoms

    half = built.subset(np.arange(0, built.n_atoms, 2))
    again = predicted_mask(half)
    assert len(again) == half.n_atoms
    assert again.sum() == predicted[::2].sum(), (
        "the mask must be recomputed from the subset, not carried")


def test_predicted_atoms_carry_plddt_and_experimental_ones_do_not():
    """The B-factor column tells the two apart without any extra field."""
    model = build_full_length(_require("7WLT"), "full")
    built = model.structure
    predicted = predicted_mask(built)
    values = built.b_factor[predicted]
    assert values.min() >= 0.0 and values.max() <= 100.0, "pLDDT is a percentage"
    assert 0.0 < model.confident_fraction < 1.0


def test_the_name_and_the_warnings_say_what_it_is():
    model = build_full_length(_require("7WLT"), "full")
    assert "+AF" in model.structure.name
    warnings = model.warnings()
    assert warnings and "PREDICTION" in warnings[0]
    assert any("pLDDT" in w for w in warnings)
    assert "predicted" in model.summary()


def test_deposited_only_adds_nothing_and_says_nothing():
    """The default must be indistinguishable from not using this at all."""
    experimental = _require("7WLT")
    model = build_full_length(experimental, "none")
    assert model.structure.n_atoms == experimental.n_atoms
    assert model.n_predicted_residues == 0
    assert not is_full_length(model.structure)
    assert not predicted_mask(model.structure).any()


# ---------------------------------------------------------------- the gap filling

def test_gaps_are_anchored_on_both_flanks_and_long_ones_are_refused():
    """An internal fill is interpolated; the blade is not. That is the difference.

    A gap too long to be placed by a local fit at its two ends is left empty
    and *counted*, because a 400-residue insert positioned that way would look
    exactly like structure.
    """
    experimental = _require("7WLT")
    longest = PARAMETERS.value("full_length.max_gap")
    gaps = resolved_gaps(experimental, experimental.chains[0])
    assert gaps, "7WLT has unresolved stretches inside its range"

    model = build_full_length(experimental, "gaps")
    assert model.meta["gaps_filled"] > 0
    assert model.meta["gaps_skipped"] > 0, (
        "some gap should exceed the ceiling; if none does, the ceiling is "
        "not doing anything and the test asserts nothing")
    # Every skip is recorded with where it was, so the model reports what it
    # did not do rather than looking complete.
    assert all(len(entry) == 3 for entry in model.meta["skipped"])
    assert model.meta["gap_flank_rmsd_max"] < 8.0, (
        "a flank fit this bad means the prediction and the experiment are "
        "different conformations there, not that the gap was filled")
    reasons = [last - first + 1 > longest
               for _c, first, last in model.meta["skipped"]]
    assert any(reasons), (
        "at least one skip should be the length ceiling rather than a missing "
        "flank, or the ceiling is untested")


def test_the_blade_is_reported_as_extrapolated_and_the_gaps_as_interpolated():
    """The two halves are different claims and the module says which is which."""
    labels = dict((key, tip) for key, _label, tip in FILL_MODES)
    assert "interpolated" in labels["gaps"]
    assert "extrapolated" in labels["blade"]


def test_the_blade_symmetry_is_an_output_not_an_input():
    """Each protomer is grafted independently, so the C3 is a measurement.

    Replicating one graft by the measured C3 would have made this number zero
    by construction and meaningless. It is not zero on 11ZC.
    """
    tight = build_full_length(_require("7WLT"), "blade")
    assert tight.blade_c3_deviation < 1.0

    loose = build_full_length(_require("11ZC"), "blade")
    assert loose.blade_c3_deviation > 1.0, (
        "11ZC is the low-resolution flat entry; if its independently placed "
        "blades agree perfectly, symmetry is being imposed somewhere")


def test_the_full_length_parameters_are_registered_with_their_reasons():
    for key in ("full_length.gap_anchor_window", "full_length.max_gap"):
        parameter = PARAMETERS.get(key)
        assert parameter is not None, key
        assert parameter.citation == "method_choice" and parameter.source_note
