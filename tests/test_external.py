"""External predictors via the ProtVar API.

These run **offline against the on-disk cache**, so the suite neither needs a
network nor hammers a public service. Anything not cached is skipped rather
than fetched.

As in Round 16, nothing here compares a predictor against the variant
phenotype labels. That comparison needs the new pre-registration in Round 22.
"""

import json

import pytest

from piezo1.analysis.external import (PROTVAR_BASE, PROTVAR_CITATION,
                                      PROTVAR_LICENCE, ExternalScores,
                                      ProtVarClient)
from piezo1.core.annotations import load_annotations


@pytest.fixture(scope="module")
def offline():
    """A client that will only ever read the cache."""
    client = ProtVarClient(offline=True)
    if not any(client.cache_dir.glob("*.json")):
        pytest.skip("ProtVar cache empty — run an online annotation once")
    return client


# --------------------------------------------------------------------------
# Licence and attribution
# --------------------------------------------------------------------------

def test_licence_and_attribution_are_recorded():
    """The whole reason for using this API is its licence; record it."""
    assert "CC BY 4.0" in PROTVAR_LICENCE
    assert "creativecommons.org" in PROTVAR_LICENCE
    assert "ProtVar" in PROTVAR_CITATION
    assert PROTVAR_BASE.startswith("https://www.ebi.ac.uk/ProtVar")
    assert ExternalScores("Q92508", 1).licence == PROTVAR_LICENCE


# --------------------------------------------------------------------------
# The disambiguation trap
# --------------------------------------------------------------------------

def test_position_query_returns_only_conservation(offline):
    """Without ``mt`` the payload cannot be attributed to a substitution.

    The endpoint returns nineteen entries per predictor, one per possible
    mutation, in an order the response never states. Guessing an alphabetical
    ordering would be a silent correctness bug, so a position-only query keeps
    just the conservation score, which is genuinely position-level.
    """
    result = offline.scores("Q92508", 2456)
    if result is None:
        pytest.skip("position 2456 not cached")
    assert result.conservation is not None
    assert result.alphamissense is None
    assert result.eve is None
    assert result.esm1b is None


def test_mutant_query_returns_one_score_each(offline):
    result = offline.scores("Q92508", 2456, "H")
    if result is None:
        pytest.skip("R2456H not cached")
    assert result.mutant == "H"
    assert result.alphamissense is not None
    assert result.eve is not None
    assert result.esm1b is not None
    assert 0.0 <= result.alphamissense <= 1.0
    assert 0.0 <= result.eve <= 1.0


def test_different_mutations_give_different_scores(offline):
    """The point of the ``mt`` parameter: substitutions must differ."""
    scores = {}
    for mutant in ("H", "C", "K", "P"):
        r = offline.scores("Q92508", 2456, mutant)
        if r and r.alphamissense is not None:
            scores[mutant] = r.alphamissense
    if len(scores) < 3:
        pytest.skip("not enough R2456 substitutions cached")
    assert len(set(scores.values())) == len(scores), scores


def test_foldx_is_keyed_by_mutation(offline):
    """Unlike /score, this endpoint labels every entry with its mutatedType."""
    table = offline.foldx("Q92508", 2456)
    if not table:
        pytest.skip("FoldX for 2456 not cached")
    assert set(table) <= set("ACDEFGHIKLMNPQRSTVWY")
    for mutant, entry in table.items():
        assert entry["wild_type"] == "R", f"{mutant}: wrong wild type"
        assert isinstance(entry["ddg"], float)
    # Proline into a helix should be among the most destabilising.
    if "P" in table and "K" in table:
        assert table["P"]["ddg"] > table["K"]["ddg"]


# --------------------------------------------------------------------------
# Behaviour
# --------------------------------------------------------------------------

def test_offline_client_never_reaches_the_network(offline, monkeypatch):
    def explode(*args, **kwargs):
        raise AssertionError("offline client attempted a network call")
    monkeypatch.setattr("urllib.request.urlopen", explode)
    offline.scores("Q92508", 2456, "H")
    offline.scores("Q92508", 999999, "A")     # certainly not cached


def test_missing_data_degrades_rather_than_raises(offline):
    """A missing external score should weaken an analysis, not abort it."""
    assert offline.scores("Q92508", 999999, "A") is None
    assert offline.foldx("Q92508", 999999) == {}
    assert offline.pockets("Q92508", 999999) == []


def test_cache_is_json_and_reusable(offline):
    files = list(offline.cache_dir.glob("*.json"))
    assert files
    for path in files[:5]:
        json.loads(path.read_text())      # must not raise


def test_stats_are_tracked(offline):
    before = dict(offline.stats)
    offline.scores("Q92508", 2456, "H")
    after = offline.stats
    assert sum(after.values()) > sum(before.values())


# --------------------------------------------------------------------------
# Cross-validation of our own numbering
# --------------------------------------------------------------------------

def test_protvar_confirms_our_wild_type_residues(offline):
    """An independent check on the residue numbering, from outside the project.

    ProtVar reports the wild-type residue it holds for each position. If our
    variant table had drifted against Q92508 — an off-by-one anywhere in the
    numbering work — this is where it would show.
    """
    ann = load_annotations("human")
    checked = mismatches = 0
    for variant in ann.variants:
        if not (variant.residue and variant.wt_aa and variant.mut_aa):
            continue
        if len(variant.wt_aa) != 1 or len(variant.mut_aa) != 1:
            continue
        table = offline.foldx("Q92508", variant.residue)
        entry = table.get(variant.mut_aa)
        if not entry or not entry.get("wild_type"):
            continue
        checked += 1
        if entry["wild_type"] != variant.wt_aa:
            mismatches += 1
    if checked < 10:
        pytest.skip(f"only {checked} variants cached")
    assert mismatches == 0, f"{mismatches} of {checked} wild types disagree"


def test_nonsense_variants_get_only_conservation(offline):
    """Missense predictors have nothing to say about a stop codon.

    Their absence for such variants is correct behaviour, not missing data.
    """
    ann = load_annotations("human")
    stops = [v for v in ann.variants
             if v.residue and v.mut_aa in ("*", "X", None) and v.wt_aa]
    if not stops:
        pytest.skip("no nonsense variants in the table")
    for variant in stops[:5]:
        result = offline.scores("Q92508", variant.residue)
        if result is None:
            continue
        assert result.alphamissense is None
