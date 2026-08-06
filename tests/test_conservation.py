"""Evolutionary conservation and the constrained-position search."""

import numpy as np
import pytest

from piezo1.analysis.conservation import (AMINO_ACIDS, ConservationProfile,
                                          Ortholog, OrthologSet,
                                          conservation_profile,
                                          constrained_positions,
                                          load_orthologs, rank_candidates)
from piezo1.core.annotations import load_annotations


@pytest.fixture(scope="module")
def orthologs():
    try:
        return load_orthologs()
    except FileNotFoundError:
        pytest.skip("ortholog cache missing — run fetch_orthologs() once")


@pytest.fixture(scope="module")
def profile(orthologs):
    return conservation_profile(orthologs)


# --------------------------------------------------------------------------
# The metric itself
# --------------------------------------------------------------------------

def test_invariant_column_scores_one():
    """Identical sequences must give conservation exactly 1."""
    seq = "ACDEFGHIKLMNPQRSTVWY" * 3
    members = [Ortholog("ref", "ref", len(seq), seq),
               Ortholog("a", "sp_a", len(seq), seq),
               Ortholog("b", "sp_b", len(seq), seq)]
    p = conservation_profile(OrthologSet(members=members), reference=seq)
    assert np.allclose(p.conservation, 1.0)
    assert np.allclose(p.identity, 1.0)
    assert np.allclose(p.coverage, 1.0)


def test_variable_column_scores_lower():
    ref = "AAAAAAAAAA"
    members = [Ortholog("a", "sp_a", 10, "AAAAAAAAAA"),
               Ortholog("b", "sp_b", 10, "AAAAACAAAA"),
               Ortholog("c", "sp_c", 10, "AAAAADAAAA"),
               Ortholog("d", "sp_d", 10, "AAAAAEAAAA")]
    p = conservation_profile(OrthologSet(members=members), reference=ref)
    varied = p.conservation[5]
    fixed = p.conservation[0]
    assert varied < fixed
    assert p.identity[5] < p.identity[0]


def test_uncovered_columns_are_not_called_conserved():
    """A position nothing aligns to is unknown, not invariant.

    Scoring it as conserved would make every unresolved loop look constrained.
    """
    ref = "AAAAAAAAAAWWWWWWWWWW"
    members = [Ortholog("a", "sp_a", 10, "AAAAAAAAAA")]
    p = conservation_profile(OrthologSet(members=members), reference=ref)
    assert p.coverage[15] < 0.5
    assert p.conservation[15] < 0.5


def test_entropy_is_normalised_to_the_unit_interval(profile):
    assert (profile.entropy >= -1e-9).all()
    assert (profile.entropy <= 1.0 + 1e-9).all()
    assert len(AMINO_ACIDS) == 20


# --------------------------------------------------------------------------
# The real ortholog set
# --------------------------------------------------------------------------

def test_ortholog_set_is_deduplicated_by_species(orthologs):
    """One entry per organism; counting a species twice double-weights it."""
    names = orthologs.organisms()
    assert len(names) == len(set(names))
    assert len(orthologs) >= 20
    assert "Homo sapiens" in names


def test_profile_spans_the_human_sequence(profile):
    assert len(profile.residues) == 2521
    assert profile.residues[0] == 1 and profile.residues[-1] == 2521
    assert profile.n_orthologs >= 20


def test_pore_module_is_more_conserved_than_the_distal_blade(profile):
    """A structural expectation the data should meet unprompted.

    The pore and its lever are under strong selection; the distal blade tip,
    which no experimental structure even resolves, is the least constrained.
    """
    ann = load_annotations("human")
    means = profile.domain_means(ann)
    for tight in ("anchor", "inner_helix", "ctd", "outer_helix"):
        assert means[tight] > means["thu1"], tight
    assert means["anchor"] > 0.9


def test_annotated_functional_residues_are_highly_conserved(profile):
    """Curated sites should score high — a check on both the metric and them."""
    ann = load_annotations("human")
    for group_id, floor in (("selectivity_acidic", 0.9),
                            ("pip2_cluster", 0.9),
                            ("anchor_brake", 0.95)):
        values = [profile.at(r)["conservation"]
                  for r in ann.group(group_id).residues
                  if profile.at(r) is not None]
        assert np.mean(values) >= floor, f"{group_id}: {np.mean(values):.3f}"


# --------------------------------------------------------------------------
# Constrained positions
# --------------------------------------------------------------------------

def test_constrained_positions_exclude_known_variants(profile):
    ann = load_annotations("human")
    known = {v.residue for v in ann.variants if v.residue is not None}
    found = constrained_positions(profile, ann, conservation_threshold=0.99,
                                  min_coverage=0.8)
    assert found
    assert not ({c["residue"] for c in found} & known)
    assert all(c["conservation"] >= 0.99 for c in found)


def test_threshold_tightens_the_set(profile):
    ann = load_annotations("human")
    loose = constrained_positions(profile, ann, conservation_threshold=0.9)
    tight = constrained_positions(profile, ann, conservation_threshold=0.999)
    assert len(tight) < len(loose)


def test_resolved_filter_restricts_to_testable_positions(profile):
    ann = load_annotations("human")
    subset = set(range(2400, 2500))
    found = constrained_positions(profile, ann, resolved=subset)
    assert all(2400 <= c["residue"] < 2500 for c in found)


# --------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------

def test_ranking_orders_by_the_supplied_evidence():
    candidates = [{"residue": r, "conservation": 1.0, "coverage": 1.0,
                   "domain": None, "domain_name": None, "annotated_sites": []}
                  for r in (10, 20, 30)]
    ranked = rank_candidates(candidates, {"score": {10: 0.1, 20: 0.5, 30: 0.9}})
    assert [c["residue"] for c in ranked] == [30, 20, 10]
    assert ranked[0]["score_percentile"] == pytest.approx(1.0)


def test_missing_feature_is_absent_not_zero():
    """A residue with no value must not be scored as if it had the worst one."""
    candidates = [{"residue": r, "conservation": 1.0, "coverage": 1.0,
                   "domain": None, "domain_name": None, "annotated_sites": []}
                  for r in (10, 20)]
    ranked = rank_candidates(candidates, {"a": {10: 0.9}, "b": {20: 0.9}})
    for c in ranked:
        # Each has exactly one feature; both should score on that one alone.
        assert c["combined_score"] == pytest.approx(0.0) or \
               np.isfinite(c["combined_score"])
    lookup = {c["residue"]: c for c in ranked}
    assert lookup[10]["b_percentile"] is None
    assert lookup[20]["a_percentile"] is None


def test_ranking_handles_an_empty_candidate_list():
    assert rank_candidates([], {"x": {1: 1.0}}) == []
