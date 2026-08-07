"""The full-length model, and the two numbers that keep it honest.

Aim A1 asks for cryo-EM plus prediction for the unresolved blade, with the seam
visible rather than hidden. The model builds — and the interesting part is what
it reports about itself.

The graft is anchored on residues *near* the seam. Fitting the whole 1279-residue
overlap instead gives 19.0 Å RMSD, because the AlphaFold and cryo-EM blades are
different conformations of a long flexible arm, and spreading that error into
the join misplaces the graft. Anchoring locally fits the seam to 2.4 Å — and
then the far end of the blade sits 75 Å from the experiment, which is the honest
caveat a good local fit would otherwise hide.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.structure.hybrid import (DEFAULT_PREDICTED_MODEL,
                                     build_hybrid_model)


@pytest.fixture(scope="module")
def hybrid():
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure

    for name in ("8YEZ.cif", DEFAULT_PREDICTED_MODEL):
        if not (STRUCTURE_DIR / name).exists():
            pytest.skip(f"{name} not downloaded; run python -m piezo1.io.fetch")
    return build_hybrid_model(Structure.from_file(STRUCTURE_DIR / "8YEZ.cif"))


# ------------------------------------------------------ what it joins

def test_it_grafts_exactly_what_the_experiment_does_not_resolve(hybrid):
    assert hybrid.seam_residue == 570
    assert hybrid.meta["grafted_residues"] == 569
    assert hybrid.res_seq[hybrid.predicted].max() < hybrid.seam_residue
    assert hybrid.res_seq[hybrid.experimental_only].min() >= hybrid.seam_residue


def test_every_atom_carries_its_provenance(hybrid):
    """The point of the module: no measurement may silently average across."""
    assert set(np.unique(hybrid.source)) == {"experimental", "predicted"}
    assert hybrid.predicted.sum() + hybrid.experimental_only.sum() == len(hybrid.xyz)
    assert np.isnan(hybrid.plddt[hybrid.experimental_only]).all()
    assert np.isfinite(hybrid.plddt[hybrid.predicted]).all()


def test_one_protomer_not_the_trimer(hybrid):
    """The prediction is a monomer; grafting onto three would assume placement."""
    assert hybrid.meta["chain"]
    assert hybrid.meta["shared_residues"] == 1279


# --------------------------------------------- the local anchor is the point

def test_the_seam_is_anchored_locally_and_fits_well(hybrid):
    assert hybrid.meta["anchor_residues"] < hybrid.meta["shared_residues"]
    assert hybrid.overlap_rmsd < 4.0, (
        f"the seam fit has degraded to {hybrid.overlap_rmsd:.2f} A")


def test_a_global_anchor_is_much_worse_and_that_is_why_it_is_not_used():
    """The measurement that justifies the anchor window, not an assertion."""
    from piezo1.config import STRUCTURE_DIR
    from piezo1.core import Structure

    if not (STRUCTURE_DIR / "8YEZ.cif").exists():
        pytest.skip("8YEZ not downloaded")
    experimental = Structure.from_file(STRUCTURE_DIR / "8YEZ.cif")
    local = build_hybrid_model(experimental)
    globalfit = build_hybrid_model(experimental, anchor_window=10_000)

    assert globalfit.overlap_rmsd > 4 * local.overlap_rmsd, (
        f"global {globalfit.overlap_rmsd:.1f} A vs local "
        f"{local.overlap_rmsd:.1f} A — the anchor window has stopped mattering")
    assert globalfit.overlap_rmsd == pytest.approx(19.0, abs=1.5)


def test_the_confidence_threshold_is_registered_not_invented():
    from piezo1.parameters import PARAMETERS

    parameter = PARAMETERS.get("hybrid.plddt_confident")
    assert parameter.default == 70.0
    assert parameter.kind == "convention"
    assert "AlphaFold" in parameter.source_note


def test_the_anchor_window_is_a_registered_parameter():
    from piezo1.parameters import PARAMETERS

    parameter = PARAMETERS.get("hybrid.anchor_window")
    assert parameter.default == 200
    assert parameter.kind == "method"
    assert "19.0" in parameter.source_note, (
        "the reason for the window is the measurement that motivated it")


# ------------------------------------------------- what it says about itself

def test_it_reports_the_global_disagreement_a_good_local_fit_would_hide(hybrid):
    """2.4 Å at the seam, 75 Å overall. Both numbers or neither."""
    assert hybrid.global_rmsd > 10 * hybrid.overlap_rmsd
    text = hybrid.summary()
    assert f"{hybrid.overlap_rmsd:.2f}" in text
    assert f"{hybrid.global_rmsd:.1f}" in text


def test_the_prediction_is_weakest_where_it_is_relied_on(hybrid):
    """Mean pLDDT 64.5 on the graft against 74.2 where the experiment covers."""
    confident = hybrid.confident_prediction.sum() / hybrid.predicted.sum()
    assert confident < 0.7, f"{confident:.0%} above the pLDDT threshold"


def test_it_warns_before_it_is_used(hybrid):
    warnings = hybrid.warnings()
    joined = " ".join(warnings).lower()
    assert "pldd" in joined
    assert "not the same conformation" in joined
    assert "model" in joined and "experimental_only" in joined


def test_the_default_selection_for_an_analysis_excludes_the_graft(hybrid):
    assert hybrid.experimental_only.sum() > hybrid.predicted.sum()
    assert not hybrid.predicted[hybrid.experimental_only].any()


def test_it_refuses_a_structure_it_cannot_graft_onto():
    from piezo1.core.structure import Structure

    n = 10
    tiny = Structure(
        name="tiny", xyz=np.zeros((n, 3)), element=np.array(["C"] * n),
        atom_name=np.array(["CA"] * n), res_name=np.array(["ALA"] * n),
        res_seq=np.arange(1, n + 1), chain=np.array(["A"] * n),
        hetero=np.zeros(n, bool), b_factor=np.zeros(n),
        occupancy=np.ones(n), alt_loc=np.array([""] * n),
        entity=np.zeros(n, int))
    with pytest.raises(ValueError, match="graft onto|C-alphas"):
        build_hybrid_model(tiny)
