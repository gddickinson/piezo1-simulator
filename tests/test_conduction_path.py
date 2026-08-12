"""The pathway option, and the one thing it must never do.

Every conduction number this project has recorded was computed on the axial
route. This module lets a caller pick a different one, so the load-bearing test
is not that the lateral options work — it is that the default is **the same
object**, not merely the same number, so nothing recorded can drift because a
new option exists.

The rest is the usual discipline: a truncation is only believed after it has
been shown to cut where it says it cuts, and to refuse when it cannot.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.hydration import WATER_RADIUS_NM, load_grid, predict_wetting
from piezo1.config import STRUCTURE_DIR
from piezo1.core import Structure
from piezo1.core.annotations import load_annotations
from piezo1.physics.conduction_path import PATHWAYS, conduction_path


def _load(pdb: str) -> Structure:
    path = STRUCTURE_DIR / f"{pdb}.cif"
    if not path.exists():
        pytest.skip(f"{pdb}.cif not downloaded — run python -m piezo1.io.fetch")
    return Structure.from_file(path)


def _profile(structure):
    from piezo1.structure.pore import pore_profile
    from piezo1.structure.protomers import protomer_blocks
    from piezo1.structure.superpose import detect_c3_axis

    blocks, _ = protomer_blocks(structure)
    return pore_profile(structure, detect_c3_axis(blocks))


@pytest.fixture(scope="module")
def intermediate():
    """8IXO — Liu et al.'s intermediate-open S2472E structure."""
    structure = _load("8IXO")
    return structure, _profile(structure)


# ----------------------------------------------- the default cannot move

def test_the_axial_pathway_returns_the_same_profile_object(intermediate):
    """Identity, not equality.

    Every recorded conduction result was computed on the full axis. If the
    default returned a *copy*, a future edit to the copying could move those
    numbers without anything in this suite noticing.
    """
    structure, profile = intermediate
    assert conduction_path(structure, profile, "axial").profile is profile


def test_the_axial_pathway_drops_nothing_and_says_so(intermediate):
    structure, profile = intermediate
    path = conduction_path(structure, profile, "axial")
    assert path.is_axial
    assert path.dropped_entry == path.dropped_exit == 0
    assert "axial pathway" in path.caveat()
    assert "closed cap top" in path.caveat(), (
        "the axial caveat must name what it is making the ions traverse")


def test_an_unknown_pathway_is_refused_rather_than_defaulted(intermediate):
    structure, profile = intermediate
    with pytest.raises(ValueError, match="unknown pathway"):
        conduction_path(structure, profile, "sideways")


# ------------------------------------------- calibration: it cuts where it says

def test_lateral_entry_cuts_above_the_cap_constriction(intermediate):
    """The cut is checked against the curated residue, not against a count.

    R2295 is the residue Liu et al. name as the closed top of the cap, and the
    truncated path must start below the last slice that residue lines — not
    near it, below it.
    """
    structure, profile = intermediate
    cap = set(load_annotations("mouse").group("cap_constriction").residues)
    path = conduction_path(structure, profile, "lateral_entry")
    assert not path.refused, path.refused
    assert path.dropped_entry > 0
    assert not any(cap.intersection(sl.lining) for sl in path.profile.slices), (
        "the cap constriction is still on the truncated path")


def test_lateral_exit_cuts_below_the_cytoplasmic_constriction(intermediate):
    structure, profile = intermediate
    neck = set(load_annotations("mouse").group("ctd_constriction").residues)
    path = conduction_path(structure, profile, "lateral_exit")
    assert not path.refused, path.refused
    assert path.dropped_exit > 0
    assert not any(neck.intersection(sl.lining) for sl in path.profile.slices)


def test_the_full_lateral_route_is_the_two_cuts_together(intermediate):
    structure, profile = intermediate
    entry = conduction_path(structure, profile, "lateral_entry")
    exit_ = conduction_path(structure, profile, "lateral_exit")
    both = conduction_path(structure, profile, "lateral")
    assert both.dropped_entry == entry.dropped_entry
    assert both.dropped_exit == exit_.dropped_exit
    assert len(both.profile.z) < min(len(entry.profile.z), len(exit_.profile.z))


@pytest.mark.parametrize("pathway", PATHWAYS)
def test_every_pathway_keeps_the_slices_in_order(intermediate, pathway):
    """`solve_pnp` sorts by z, but a scrambled profile would break the
    lining-based checks above and every plot drawn from it."""
    structure, profile = intermediate
    z = np.asarray(conduction_path(structure, profile, pathway).profile.z)
    assert np.all(np.diff(z) > 0)


# --------------------------------------------------- calibration: it can say no

def test_a_paralogue_is_refused_rather_than_truncated_somewhere():
    """The ranges are residue numbers, so a wrong reading cuts the wrong place.

    A confident lateral current on a PIEZO2 entry — cut at whatever happens to
    sit at PIEZO1's cap residue — is exactly the failure mode this project
    audits for, and it is worse than no answer.
    """
    structure = _load("6KG7")
    path = conduction_path(structure, _profile(structure), "lateral")
    assert path.refused
    assert "numbering" in path.refused
    assert path.profile is not None, "a refusal still returns the full axis"
    assert "refused" in path.caveat()


def test_the_splice_isoform_entry_is_refused_too():
    """6LQI is PIEZO1 and still cannot be read by residue number."""
    structure = _load("6LQI")
    assert conduction_path(structure, _profile(structure), "lateral").refused


# ------------------------------------------------------------- the measurement

def test_the_lateral_route_opens_the_intermediate_structure(intermediate):
    """The finding, and the reason the option exists.

    8IXO's transmembrane gate has demonstrably dilated — the V2476 side-chain
    diagonal is 14.2 A against 7.7 A on the curved 7WLT — and the axial model
    still refuses it, on constrictions Liu et al. report are bypassed. On their
    route it conducts.
    """
    grid = load_grid()
    if not grid.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")

    structure, profile = intermediate
    axial = conduction_path(structure, profile, "axial")
    lateral = conduction_path(structure, profile, "lateral")

    assert not predict_wetting(structure, axial.profile,
                               grid=grid).conductive
    assert predict_wetting(structure, lateral.profile, grid=grid).conductive
    assert lateral.profile.radius.min() / 10.0 >= WATER_RADIUS_NM


def test_composing_the_verdict_on_the_truncated_profile_fails_to_separate():
    """Why `analysis.conduction` exists, demonstrated rather than asserted.

    Round 84d evaluated **both** halves of the wetting verdict on the truncated
    profile and reported "the lateral route does not separate open from closed"
    as an honest negative, pinned by a test. It was not a property of the
    channel. The Rao score is a sum over lining residues carried by narrow
    slices, and the truncation removes them, so no entry anywhere reaches the
    cutoff and the verdict falls back on a residual radius.

    This keeps the broken composition as a *demonstration*: curved entries
    still pass it, which is the defect. The composed verdict in
    `analysis.conduction` refuses them, and `test_conduction_verdict` checks
    that against Liu et al.'s Figure 5D ordering.
    """
    grid = load_grid()
    if not grid.available:
        pytest.skip("CHAP grid not downloaded — run python -m piezo1.io.fetch")

    from piezo1.analysis.conduction import conduction_verdict

    naive, composed = [], []
    for pdb in ("8IXO", "7WLT", "6B3R", "8IMZ"):
        structure = _load(pdb)
        profile = _profile(structure)
        path = conduction_path(structure, profile, "lateral")
        if predict_wetting(structure, path.profile, grid=grid).conductive:
            naive.append(pdb)
        if conduction_verdict(structure, profile, "lateral",
                              grid=grid).conductive:
            composed.append(pdb)

    assert "8IXO" in naive and "8IXO" in composed
    assert len(naive) > 1, (
        "the broken composition no longer lets curved entries through, so this "
        "demonstration has stopped demonstrating anything")
    assert composed == ["8IXO"], (
        f"the composed verdict should refuse every curved entry; got {composed}")


def test_the_caveat_says_the_portal_is_not_modelled(intermediate):
    """The one thing a lateral current must never be read as.

    The truncated end slice becomes the mouth and its radius is the pore's,
    not the portal's — which this project does not measure. So the current is
    an upper bound, and the sentence saying so cannot be optional.
    """
    structure, profile = intermediate
    caveat = conduction_path(structure, profile, "lateral").caveat()
    assert "NOT modelled" in caveat
    assert "upper bound" in caveat
