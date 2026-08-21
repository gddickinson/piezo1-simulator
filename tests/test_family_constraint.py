"""The census's constraint read on this project's structures, domains and variants.

Three instruments live here and each is calibrated before it is believed:

* :func:`domain_constraint` and :func:`paralogue_asymmetry` **replicate** census
  results on boundaries the census did not choose. A replication is only worth
  something if it could have failed, so the tests check the agreement *and* the
  one place it does not agree.
* :func:`blade_gradient` **contradicts** a census finding. A contradiction needs
  more support than an agreement, so its mechanism is checked, not just its
  direction.
* :func:`pore_module_enrichment` runs a Fisher test, which is checked against a
  hand-computable table before it is run on real counts.
"""

from __future__ import annotations

import numpy as np
import pytest

from piezo1.analysis.disease_geography import (both_partitions,
                                               boundary_disagreement,
                                               constraint_classifier,
                                               fisher_exact_greater,
                                               pathogenic_positions,
                                               pore_module_residues)
from piezo1.analysis.family_constraint import (StructureRefusal, blade_gradient,
                                               compare_with_own_conservation,
                                               constraint_on_structure,
                                               domain_constraint,
                                               paralogue_asymmetry)
from piezo1.core.family import load_constraint, load_family_findings


# --------------------------------------------------- replication on our domains

def test_the_pore_machinery_is_the_most_constrained_on_our_boundaries():
    """The census's headline, re-measured on a partition it did not choose.

    Its bands and ours put the anchor 141 residues apart, so this is a test of
    the protein rather than of a boundary choice.
    """
    rows = {d.domain: d for d in domain_constraint() if d.mean is not None}
    pore_side = ("anchor", "ctd", "inner_helix", "outer_helix")
    blade_side = ("thu1", "thu2")
    for name in pore_side:
        assert rows[name].vs_whole > 0.1, f"{name} is not above the average"
    for name in blade_side:
        assert rows[name].vs_whole < 0.01, f"{name} is not at or below it"
    assert min(rows[n].mean for n in pore_side) > max(rows[n].mean for n in blade_side)


def test_the_cap_is_the_exception_the_census_says_it_is():
    """The CED is pore machinery and is *less* alike between paralogues than
    the protein as a whole. Checked in both pairs this project can frame."""
    for pair in ("PIEZO1_vs_PIEZO2", "PIEZO1_vs_piezo3"):
        rows = {d.domain: d for d in paralogue_asymmetry(pair)}
        assert not rows["cap"].above_average, f"{pair}: cap is not the exception"
        for name in ("inner_helix", "ctd", "anchor"):
            assert rows[name].above_average, f"{pair}: {name} should be above"


def test_our_identities_reproduce_the_census_where_it_matters():
    """Different alignment, different boundaries, same answer — or a real gap."""
    census = {r["domain"]: float(r["identity"])
              for r in load_family_findings().table("paralogue_identity")
              if r["pair"] == "PIEZO1_vs_PIEZO2"}
    ours = {d.domain: d.identity for d in paralogue_asymmetry("PIEZO1_vs_PIEZO2")}
    # The cap is the finding, and it lands within a percentage point.
    assert abs(ours["cap"] - census["CED"]) < 0.02
    for name, census_name in (("inner_helix", "inner_helix"), ("ctd", "CTD")):
        assert abs(ours[name] - census[census_name]) < 0.06, name


def test_a_pair_we_hold_no_boundaries_for_is_refused_rather_than_answered():
    """Framing in PIEZO2 and indexing with PIEZO1's ranges gave the cap 0.85
    where the census measures 0.35, and announced nothing."""
    with pytest.raises(KeyError):
        paralogue_asymmetry("PIEZO2_vs_piezo3")


# ------------------------------------------- the contradiction, and its mechanism

def test_the_distal_blade_finding_is_a_property_of_the_bands_not_the_blades():
    """The census's distal-beats-proximal result, contradicted with a mechanism.

    Three things have to hold together for the contradiction to be worth
    anything: its bands must reproduce here (or the import is wrong, not the
    finding), the ordering must reverse on the units, and the linker must score
    the same either side — which is what makes it composition rather than
    biology.
    """
    result = blade_gradient()
    assert abs(result["band_distal"] - result["census_distal"]) < 0.02
    assert abs(result["band_proximal"] - result["census_proximal"]) < 0.02
    assert result["holds_on_census_bands"] is True
    assert result["holds_on_thu_units"] is False
    assert abs(result["linker_distal"] - result["linker_proximal"]) < 0.01
    assert result["linker_fraction_proximal"] > 2 * result["linker_fraction_distal"]
    assert "boundary-dependent" in result["verdict"]


def test_moving_the_distal_proximal_split_moves_the_answer():
    """The split is a registered parameter because the finding turns on it."""
    from piezo1.parameters import PARAMETERS as _P

    base = blade_gradient()["unit_distal"]
    _P.set_value("family.distal_last_thu", 3)
    try:
        moved = blade_gradient()["unit_distal"]
    finally:
        _P.reset("family.distal_last_thu")
    assert moved != pytest.approx(base)
    assert blade_gradient()["unit_distal"] == pytest.approx(base)


# --------------------------------------------------------- the numbering gate

def test_the_track_is_refused_on_an_entry_it_cannot_be_read_on(structure_by_id):
    """A PIEZO2 entry coloured by PIEZO1's constraint would look entirely
    plausible and be wrong at every residue."""
    piezo2 = structure_by_id("6KG7")
    if piezo2 is None:
        pytest.skip("6KG7 is not downloaded")
    result = constraint_on_structure(piezo2)
    assert isinstance(result, StructureRefusal)
    assert not result
    assert "numbering" in result.reason


def test_a_mouse_entry_is_carried_through_the_alignment_not_an_offset(curved_structure):
    result = constraint_on_structure(curved_structure)
    assert result and result.numbering == "mouse" and result.converted
    assert result.coverage > 0.9
    # An offset would be uniform. The map is not: mouse 2473 is human 2447,
    # a 26-residue difference, while mouse 600 is human 594.
    from piezo1.core.sequence import mouse_to_human
    assert mouse_to_human(2473) - 2473 != mouse_to_human(600) - 600


def test_unscored_atoms_are_nan_on_the_structure_not_zero(curved_structure):
    result = constraint_on_structure(curved_structure)
    values = np.asarray(result.per_atom)
    assert np.isnan(values).any(), "some atoms must be unscored on a real entry"
    assert np.nanmin(values) > 0.0, "a real JSD is never exactly zero here"


# ------------------------------------------------------------- the cross-check

def test_our_own_conservation_and_the_census_agree():
    """Two evolutionary measurements sharing no data and no statistic.

    Ours is Shannon entropy over 61 fetched vertebrate orthologues; theirs is
    Jensen-Shannon divergence over 174 genome-backed loci. Strong agreement is
    mutual support; this test exists so that if either route changes, the
    disagreement surfaces here rather than inside a ranking.
    """
    result = compare_with_own_conservation()
    if result is None:
        pytest.skip("no fetched orthologues; run the conservation analysis first")
    assert result.n > 1000
    assert result.spearman > 0.7, result.summary()
    assert result.agree


# --------------------------------------------------- the Fisher test, calibrated

def test_fisher_exact_matches_a_table_computed_by_hand():
    """Before it is run on real counts. The tail summed decides the answer, and
    testing the other one would be a different claim entirely."""
    # 2x2 [[a=3, b=1], [c=1, d=3]]: P(X >= 3) = (C(4,3)C(4,1) + C(4,4)C(4,0))/C(8,4)
    assert fisher_exact_greater(3, 1, 1, 3) == pytest.approx((16 + 1) / 70)
    # A table with no association gives a large P; a perfect one a small P.
    assert fisher_exact_greater(5, 5, 5, 5) > 0.5
    assert fisher_exact_greater(10, 0, 0, 10) < 0.001
    # It is a probability.
    for table in ((3, 1, 1, 3), (0, 4, 4, 0), (1, 1, 1, 1)):
        assert 0.0 <= fisher_exact_greater(*table) <= 1.0


def test_the_two_partitions_disagree_and_the_disputed_band_holds_disease():
    """The honest form of the re-test: which answer you get depends on where
    the pore module is said to start, and the 120 residues in dispute are not
    empty of disease."""
    result = both_partitions()
    ours, census = result["results"]["ours"], result["results"]["census"]
    assert census.n_region > ours.n_region
    assert census.odds_ratio > ours.odds_ratio
    assert census.significant and not ours.significant
    assert "boundary" in result["verdict"]
    disputed = result["disputed"]
    assert disputed["n_residues"] == 120
    assert len(disputed["pathogenic"]) >= 5


def test_the_boundary_disagreement_is_confined_to_the_anchor_end():
    """The pattern is the evidence that these are two papers' ranges rather
    than a systematic numbering error: three elements agree to a few residues
    and two are shifted by more than a hundred."""
    rows = {r["element"]: r for r in boundary_disagreement()}
    assert abs(rows["cap"]["start_offset"]) < 10
    assert abs(rows["ctd"]["start_offset"]) < 20
    assert rows["anchor"]["start_offset"] > 100
    assert rows["outer_helix"]["start_offset"] > 100


def test_the_constraint_score_still_classifies_against_a_different_negative_set():
    """The census scored 0.914 against ClinVar benign labels. Here the negatives
    are gnomAD population variation — not a clinical judgement at all — so
    agreement is evidence the score generalises past the labels it was checked
    against."""
    result = constraint_classifier()
    if result is None:
        pytest.skip("no cached gnomAD data")
    assert result.n_positive >= 20 and result.n_negative >= 50
    assert result.auc > 0.7
    assert result.agrees_with_census
    assert result.mean_positive > result.mean_negative


def test_pathogenic_positions_exclude_truncating_variants():
    """A nonsense variant's position says nothing about that residue — it
    removes everything downstream — and counting it would file every truncation
    under whichever domain the stop codon lands in."""
    positions = pathogenic_positions()
    assert 40 > len(positions) >= 20
    module = pore_module_residues()
    assert positions & module, "some pathogenic positions are in the pore module"
    track = load_constraint("PIEZO1")
    assert all(track.value(p) is not None for p in positions)
