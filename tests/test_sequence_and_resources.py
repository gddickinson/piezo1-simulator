"""Cross-species numbering and the integrity of the shipped resources."""

import json

import pytest

from piezo1.config import RESOURCE_DIR
from piezo1.core.annotations import load_annotations
from piezo1.core.sequence import (human_sequence, load_numbering_map,
                                  mouse_sequence)
from piezo1.io.registry import load_registry


def test_reference_sequence_lengths():
    assert len(human_sequence()) == 2521
    assert len(mouse_sequence()) == 2547


def test_numbering_map_identity_and_coverage():
    nm = load_numbering_map()
    assert 0.80 < nm.identity < 0.86
    assert len(nm.a_to_b) == 2521
    assert len(nm.b_to_a) == 2547


@pytest.mark.parametrize("human,mouse", [
    (1718, 1718),   # Yoda1 pocket - same number by coincidence
    (2075, 2091),   # Yoda1 pocket
    (2078, 2094),   # Yoda1 pocket
    (2117, 2133),   # selectivity glutamate
    (2470, 2496),   # CTD glutamate - human E2496 is a DIFFERENT residue
    (2456, 2482),   # R2456H / mouse R2482
    (2166, 2182),   # PIP2 lysine cluster
    (2169, 2185),
    (1335, 1330),   # S1330 in mouse numbering
    (2446, 2472),   # S2472E phosphomimetic
])
def test_known_cross_species_equivalences(human, mouse):
    nm = load_numbering_map()
    assert nm.to_b(human) == mouse
    assert nm.to_a(mouse) == human
    # And the amino acid must agree, which is the real check.
    assert human_sequence()[human - 1] == mouse_sequence()[mouse - 1]


def test_offset_is_not_constant():
    """If this ever becomes constant, something has gone badly wrong."""
    nm = load_numbering_map()
    offsets = {nm.to_b(r) - r for r in (100, 600, 1000, 1700, 2000, 2500)
               if nm.to_b(r) is not None}
    assert len(offsets) > 3


def test_domains_resource_is_self_consistent():
    ann = load_annotations("human")
    assert len(ann.domains) >= 15
    seq_len = len(human_sequence())
    for d in ann.domains:
        if d.start is None:
            continue
        assert 1 <= d.start <= d.end <= seq_len, d.id
        assert d.confidence in ("high", "medium", "low")
        assert d.source
    # The pore module must come out of the topology in the right order.
    ids = {d.id: d for d in ann.domains}
    assert ids["outer_helix"].end < ids["cap"].start
    assert ids["cap"].end < ids["inner_helix"].start
    assert ids["inner_helix"].end < ids["ctd"].start
    assert ids["cap"].length > 200


def test_functional_residues_match_the_sequence():
    data = json.loads((RESOURCE_DIR / "functional_residues.json").read_text())
    seq = human_sequence()
    for group in data["groups"]:
        for r in group["residues"]:
            assert r["verified"], f"{group['id']} {r['human']}"
            assert seq[r["human"] - 1] == r["human_aa"]


def test_variants_are_sequence_verified():
    ann = load_annotations("human")
    assert len(ann.variants) >= 50
    seq = human_sequence()
    unverified = [v for v in ann.variants
                  if v.residue and v.wt_aa and len(v.wt_aa) == 1
                  and seq[v.residue - 1] != v.wt_aa]
    assert not unverified, [v.id for v in unverified]


def test_registry_entries_resolve_to_files():
    reg = load_registry()
    assert len(reg) >= 15
    for rec in reg.available():
        assert rec.path.exists()
        assert rec.species in ("human", "mouse", "unknown")
    assert reg.default() is not None
