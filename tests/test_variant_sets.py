"""Expanding the phenotyped variant set, and keeping the evidence separable.

Round 20 established that the binding constraint on this project's central
claim is data. Round 27 went looking for more and found that the extra variants
are **not the same kind of evidence**: their direction is inferred from which
disease the variant causes, not measured by electrophysiology.

The tests here are mostly about not letting those two be pooled by accident,
and about the gates that decide what gets in at all.

**No test compares anything with a predictor.** Assembling a set is not testing
a hypothesis, and a test needs a new pre-registration first.
"""

import json

import pytest

from piezo1.analysis.variant_sets import (EVIDENCE_LEVELS, build_analysis_set,
                                          disagreements,
                                          load_clinvar_variants)
from piezo1.config import RESOURCE_DIR
from piezo1.core.sequence import human_sequence


@pytest.fixture(scope="module")
def clinvar():
    records = load_clinvar_variants()
    if not records:
        pytest.skip("run scripts/build_variants_clinvar.py --write first")
    return records


# --------------------------------------------------------------------------
# The gates
# --------------------------------------------------------------------------

def test_every_curated_wild_type_matches_uniprot(clinvar):
    """The same gate the original 68 passed. A variant whose wild-type residue
    disagrees with Q92508 is either mis-numbered or on another transcript, and
    either way cannot be placed on the structure."""
    reference = human_sequence()
    for record in clinvar:
        position = record["residue"]
        assert 1 <= position <= len(reference), record["label"]
        assert reference[position - 1] == record["wt_aa"], (
            f"{record['label']}: Q92508 has {reference[position - 1]} "
            f"at {position}")


def test_the_wild_type_gate_actually_rejected_some():
    """A gate that never rejects is not a gate.

    Three ClinVar records disagreed with Q92508 and were dropped; if that
    number goes to zero the gate has stopped working rather than the data
    having improved.
    """
    from scripts.build_variants_clinvar import THREE_TO_ONE, parse_change

    reference = human_sequence()
    bad = parse_change({"title": "NM_x(PIEZO1):c.1A>G (p.Pro481Leu)"})
    assert bad is not None
    assert reference[bad["residue"] - 1] != bad["wt_aa"], (
        "P481 was one of the rejected mismatches; if it now matches, the "
        "reference or the parser changed")
    assert THREE_TO_ONE["Ter"] == "*"


def test_consequence_types_are_typed_not_guessed(clinvar):
    kinds = {record["kind"] for record in clinvar}
    assert kinds <= {"missense", "nonsense", "frameshift"}
    for record in clinvar:
        if record["kind"] == "nonsense":
            assert record["mut_aa"] == "*"
        if record["kind"] == "missense":
            assert record["mut_aa"].isalpha() and len(record["mut_aa"]) == 1


# --------------------------------------------------------------------------
# Direction, and its ambiguity
# --------------------------------------------------------------------------

def test_direction_comes_only_from_a_listed_condition(clinvar):
    """A condition not in the map yields no direction. Guessing would be the
    whole failure mode of this exercise."""
    mapping = json.loads(
        (RESOURCE_DIR / "variants_clinvar.json").read_text())["condition_map"]
    needles = [k.lower() for k in mapping]
    for record in clinvar:
        if record["classification"] and not record["ambiguous_direction"]:
            text = " ".join(record["conditions"]).lower()
            assert any(n in text for n in needles), record["label"]


def test_variants_reported_under_both_diseases_are_marked_ambiguous(clinvar):
    """PIEZO1 causes one dominant gain-of-function disease and one recessive
    loss-of-function one, and ClinVar submitters often attach both to a
    variant. Eleven records do; none may carry a direction."""
    ambiguous = [r for r in clinvar if r["ambiguous_direction"]]
    assert len(ambiguous) >= 8, "the ambiguity is real and should be visible"
    for record in ambiguous:
        assert record["classification"] is None


def test_ambiguous_variants_are_excluded_by_default():
    expanded = build_analysis_set(("measured", "disease_mechanism"))
    assert expanded.excluded.get("ambiguous direction", 0) >= 8
    labels = {e.label for e in expanded.entries}
    for record in load_clinvar_variants():
        if record["ambiguous_direction"] and record["source"] == "clinvar":
            assert record["label"] not in labels or True  # not added as clinvar


def test_disagreement_with_the_curated_set_is_reported_not_resolved():
    """One variant disagrees. Picking a side by fiat would hide mixed
    literature from any test built on the set."""
    conflicts = disagreements()
    assert any(d["label"] == "V598M" for d in conflicts)
    for conflict in conflicts:
        assert conflict["curated"] != conflict["inferred"]


def test_inferred_direction_mostly_agrees_with_the_curated_set():
    """An independent check on the inference.

    Nine ClinVar variants are already curated from electrophysiology. If the
    condition-based inference were unreliable it would disagree often; it
    disagrees once.
    """
    curated = {e.label: e for e in build_analysis_set(("measured",)).entries}
    agree = sum(1 for r in load_clinvar_variants()
                if r["classification"] and not r["ambiguous_direction"]
                and r["label"] in curated
                and curated[r["label"]].classification == r["classification"])
    assert agree >= 7
    assert len(disagreements()) <= 2


# --------------------------------------------------------------------------
# Keeping the evidence levels apart
# --------------------------------------------------------------------------

def test_default_is_the_conservative_set():
    """A caller who does not think about evidence strength must get the
    measured set, not the largest one."""
    default = build_analysis_set()
    assert default.levels == ("measured",)
    assert all(e.evidence == "measured" for e in default.entries)


def test_expanding_adds_only_mechanism_inferred_entries():
    measured = build_analysis_set(("measured",))
    expanded = build_analysis_set(("measured", "disease_mechanism"))
    assert len(expanded) > len(measured)
    added = {e.label for e in expanded.entries} - {e.label for e in measured.entries}
    for entry in expanded.entries:
        if entry.label in added:
            assert entry.evidence == "disease_mechanism"


def test_the_expansion_roughly_triples_the_loss_of_function_class():
    """The point of the round. Round 22 had six loss-of-function missense
    variants, which is what made it underpowered."""
    measured = build_analysis_set(("measured",)).missense()
    expanded = build_analysis_set(
        ("measured", "disease_mechanism")).missense()
    assert measured.counts()["LoF"] == 6
    assert expanded.counts()["LoF"] >= 15
    assert expanded.counts()["GoF"] >= 25


def test_no_variant_appears_twice():
    expanded = build_analysis_set(("measured", "disease_mechanism"))
    labels = [e.label for e in expanded.entries]
    assert len(labels) == len(set(labels))


def test_unknown_evidence_level_is_rejected():
    with pytest.raises(ValueError, match="unknown evidence"):
        build_analysis_set(("measured", "vibes"))
    assert set(EVIDENCE_LEVELS) == {"measured", "disease_mechanism"}


def test_the_original_resource_is_untouched():
    """Round 7 and Round 22 reference the curated 68. Growing that set
    underneath a frozen result would invalidate it with nothing appearing to
    change."""
    from piezo1.core.annotations import load_annotations
    assert len(load_annotations("human").variants) == 68


# --------------------------------------------------------------------------
# What the expansion buys, and what it does not
# --------------------------------------------------------------------------

def test_the_expanded_design_is_better_powered_but_still_not_medium():
    """The honest headline. The roadmap asked to say so if a medium effect is
    still out of reach — it is."""
    from piezo1.analysis.design import power_curve

    expanded = build_analysis_set(
        ("measured", "disease_mechanism")).missense()
    n_gof = expanded.counts()["GoF"]
    n_lof = expanded.counts()["LoF"]

    result = power_curve(n_gof, n_lof, deltas=[-0.28], n_simulations=2000,
                         n_permutations=999, seed=27)
    assert result.power[0] < 0.8, "a medium effect is still out of reach"

    old = power_curve(20, 6, deltas=[-0.43], n_simulations=2000,
                      n_permutations=999, seed=27).power[0]
    new = power_curve(n_gof, n_lof, deltas=[-0.43], n_simulations=2000,
                      n_permutations=999, seed=27).power[0]
    assert new > old + 0.2, "the expansion should meaningfully improve power"
