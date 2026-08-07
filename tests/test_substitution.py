"""Substitution-aware spring perturbation, and its pre-registered criterion.

Round 7's blind test failed with a precise diagnostic: 99.8% of the mechanical
ΔΔG's variance was *between position*. The cause was algebraic, not statistical
— the old model scaled every contact of a residue by one number, so

    ΔΔG = (s − 1) · Q(position)

is a rank-one product in which the substitution is a scalar. Four substitutions
at R2456 could then only differ by a factor.

**No test here compares anything with a phenotype.** Round 22's protocol
requires a new pre-registration before the labels are touched again, and this is
method development, not a hypothesis test.
"""

import collections

import numpy as np
import pytest

from piezo1.analysis.substitution import (PROPERTIES, contact_scales,
                                          substitution_summary,
                                          variance_decomposition)
from piezo1.parameters import PARAMETERS

#: The pre-registered threshold, fixed in ROADMAP.md before the work was done.
WITHIN_POSITION_TARGET = 0.20


# --------------------------------------------------------------------------
# Residue properties
# --------------------------------------------------------------------------

def test_property_table_is_complete_and_sane():
    assert len(PROPERTIES) == 20
    assert PROPERTIES["G"].volume < PROPERTIES["A"].volume < PROPERTIES["W"].volume
    assert PROPERTIES["D"].charge == -1.0 and PROPERTIES["E"].charge == -1.0
    assert PROPERTIES["K"].charge == 1.0 and PROPERTIES["R"].charge == 1.0
    assert 0.0 < PROPERTIES["H"].charge < 1.0, "histidine is partly protonated"
    assert PROPERTIES["P"].special == "P" and PROPERTIES["G"].special == "G"
    assert PROPERTIES["L"].hydrophobic > PROPERTIES["D"].hydrophobic


def test_summary_reports_what_changed():
    summary = substitution_summary("R", "P")
    assert summary["known"]
    assert summary["charge_change"] == pytest.approx(-1.0)
    assert summary["introduces_proline"]
    assert not substitution_summary("R", "K")["introduces_proline"]
    assert not substitution_summary("X", "K")["known"]


# --------------------------------------------------------------------------
# The degeneracy that had to be broken
# --------------------------------------------------------------------------

def test_scales_are_per_contact_not_one_number():
    partners = list("EDRKASGLFV")
    distances = np.linspace(5.0, 14.0, 10)
    scales = contact_scales("R", "H", partners, distances,
                            np.arange(10) * 12.0)
    assert scales.shape == (10,)
    assert scales.std() > 0.05, "a per-contact model must vary across contacts"


def test_different_substitutions_perturb_different_contacts():
    """The core repair. Under the old model these patterns were identical up to
    a constant factor, so every substitution ranked positions the same way."""
    partners = list("EDRKASGLFV")
    distances = np.linspace(5.0, 14.0, 10)
    separation = np.array([1, 2, 3, 20, 40, 60, 80, 100, 120, 140], float)

    patterns = {m: contact_scales("R", m, partners, distances, separation)
                for m in "HKPC"}
    for a, b in (("H", "K"), ("H", "P"), ("K", "C"), ("P", "C")):
        correlation = np.corrcoef(patterns[a], patterns[b])[0, 1]
        assert correlation < 0.99, (
            f"R->{a} and R->{b} still have proportional patterns "
            f"(r = {correlation:.4f}); the separability is not broken")


def test_charge_change_is_felt_only_at_charged_partners():
    """R->C removes a positive charge. That should matter beside an aspartate
    and not beside a leucine."""
    distances = np.full(2, 6.0)
    scales = contact_scales("R", "C", ["D", "L"], distances, np.array([50, 50.]))
    assert scales[0] < scales[1], (
        "losing a salt bridge should soften the contact to D more than to L")


def test_proline_stiffens_only_sequence_local_contacts():
    distances = np.full(4, 7.0)
    separation = np.array([1.0, 2.0, 50.0, 80.0])
    scales = contact_scales("A", "P", ["A"] * 4, distances, separation)
    assert scales[0] > scales[2], "proline should stiffen i±1 more than i±50"
    assert scales[1] > scales[3]


def test_glycine_softens_what_the_side_chain_was_touching():
    distances = np.full(3, 6.0)
    to_glycine = contact_scales("L", "G", ["A"] * 3, distances)
    to_similar = contact_scales("L", "I", ["A"] * 3, distances)
    assert to_glycine.mean() < to_similar.mean()


def test_springs_may_weaken_but_never_invert():
    """A negative spring makes the Hessian indefinite and the quadratic form
    stops being an energy."""
    floor = PARAMETERS.value("substitution.min_scale")
    for wt in "RKDEWLG":
        for mut in "GPCDEKR":
            scales = contact_scales(wt, mut, list("DEKRAG"),
                                    np.full(6, 5.0), np.arange(6, dtype=float))
            assert (scales >= floor - 1e-12).all(), f"{wt}->{mut}: {scales}"


def test_unknown_residues_fall_back_to_no_change():
    scales = contact_scales("X", "A", ["A"], np.array([6.0]))
    assert np.allclose(scales, 1.0)


def test_distant_contacts_barely_feel_a_side_chain_change():
    near = contact_scales("A", "W", ["A"], np.array([5.0]))[0]
    far = contact_scales("A", "W", ["A"], np.array([14.0]))[0]
    assert abs(near - 1.0) > abs(far - 1.0)


# --------------------------------------------------------------------------
# Variance decomposition
# --------------------------------------------------------------------------

def test_decomposition_identity_holds():
    rng = np.random.default_rng(0)
    positions = rng.integers(0, 5, 40)
    values = rng.normal(size=40) + positions
    d = variance_decomposition(positions, values)
    assert d.between + d.within == pytest.approx(d.total, rel=1e-12)
    assert d.within_fraction + d.between_fraction == pytest.approx(1.0)


def test_decomposition_recognises_a_purely_positional_score():
    positions = np.repeat(np.arange(6), 3)
    values = positions.astype(float)          # identical within each position
    d = variance_decomposition(positions, values)
    assert d.within_fraction == pytest.approx(0.0)
    assert d.n_multi == 6


def test_decomposition_recognises_a_purely_substitutional_score():
    positions = np.repeat(np.arange(6), 3)
    values = np.tile([-1.0, 0.0, 1.0], 6)     # same spread at every position
    d = variance_decomposition(positions, values)
    assert d.within_fraction == pytest.approx(1.0)


# --------------------------------------------------------------------------
# The pre-registered criterion, on the real structures
# --------------------------------------------------------------------------

@pytest.fixture(scope="module")
def gating_model(curved_structure, flat_structure):
    from piezo1.analysis.variant_impact import VariantImpactModel
    from piezo1.core.sequence import human_sequence, mouse_to_human
    from piezo1.structure.superpose import kabsch, match_protomers
    from piezo1.ui.model_utils import protomer_blocks

    _cb, cr = protomer_blocks(curved_structure)
    _fb, fr = protomer_blocks(flat_structure)
    common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))

    def resample(st):
        out = []
        for chain in st.chains:
            mask = st.mask_ca() & (st.chain == chain)
            if mask.sum() < 300:
                continue
            index = {int(r): i for i, r in enumerate(st.res_seq[mask])}
            xyz = st.xyz[mask]
            if all(r in index for r in common):
                out.append(np.array([xyz[index[r]] for r in common], float))
        return out[:3]

    cb, fb = resample(curved_structure), resample(flat_structure)
    if len(cb) < 3 or len(fb) < 3:
        pytest.skip("need three protomers in both endpoints")
    fb = [fb[i] for i in match_protomers(cb, fb).order]
    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    displacement = (((np.vstack(fb) - centroid) @ rotation.T + translation)
                    - np.vstack(cb))

    reference = human_sequence()
    sequence = {}
    for residue in common:
        human = mouse_to_human(int(residue))
        if human and 1 <= human <= len(reference):
            sequence[int(residue)] = reference[human - 1]

    def make(aware: bool):
        return VariantImpactModel(
            coords=np.vstack(cb), residues=np.tile(common, 3),
            gating_vector=displacement, sequence=sequence,
            substitution_aware=aware)
    return make


def _score(make, aware, subset):
    from piezo1.core.sequence import human_to_mouse
    model = make(aware)
    positions, values = [], []
    for human_residue, variants in subset.items():
        mouse_residue = human_to_mouse(human_residue)
        if mouse_residue is None:
            continue
        for variant in variants:
            prediction = model.predict(mouse_residue, variant.wt_aa,
                                       variant.mut_aa)
            if prediction.modelled and np.isfinite(prediction.gating_cost_change):
                positions.append(human_residue)
                values.append(prediction.gating_cost_change)
    return np.array(positions), np.array(values)


@pytest.fixture(scope="module")
def multiply_substituted():
    from piezo1.core.annotations import load_annotations
    by_position = collections.defaultdict(list)
    for variant in load_annotations("human").variants:
        if not (variant.residue and variant.wt_aa and variant.mut_aa):
            continue
        if len(variant.wt_aa) != 1 or len(variant.mut_aa) != 1:
            continue
        if not variant.mut_aa.isalpha():
            continue
        if variant.label not in {v.label for v in by_position[variant.residue]}:
            by_position[variant.residue].append(variant)
    return {k: v for k, v in by_position.items() if len(v) > 1}


def test_multiply_substituted_positions_exist(multiply_substituted):
    assert 2456 in multiply_substituted
    assert len(multiply_substituted[2456]) == 4
    assert len(multiply_substituted) >= 5


def test_within_position_variance_meets_the_criterion(gating_model,
                                                      multiply_substituted):
    """The pre-registered success criterion for this round.

    Measured on the multiply-substituted positions, as ROADMAP.md specified.
    Singly-substituted positions contribute exactly zero within-variance by
    construction, so including them would drive the statistic down for reasons
    that have nothing to do with the model.
    """
    positions, values = _score(gating_model, True, multiply_substituted)
    decomposition = variance_decomposition(positions, values)
    assert decomposition.within_fraction > WITHIN_POSITION_TARGET, (
        f"criterion not met: {decomposition.summary()}")


def test_the_new_model_beats_the_old_one_on_its_own_diagnostic(
        gating_model, multiply_substituted):
    old = variance_decomposition(
        *_score(gating_model, False, multiply_substituted))
    new = variance_decomposition(
        *_score(gating_model, True, multiply_substituted))
    assert old.within_fraction < 0.10, "the old model was positional"
    assert new.within_fraction > 5 * old.within_fraction


def test_falling_back_without_a_sequence_reproduces_the_old_behaviour(
        gating_model, multiply_substituted):
    """Round 7's recorded result must not move.

    Its script passes no sequence, so the model keeps the uniform scale. If
    that changed, a frozen validation number would silently move.
    """
    from piezo1.analysis.variant_impact import VariantImpactModel
    from piezo1.core.sequence import human_to_mouse

    aware = gating_model(True)
    without = VariantImpactModel(
        coords=aware.coords, residues=aware.residues,
        gating_vector=aware.gating_vector)          # no sequence
    residue = human_to_mouse(2456)
    a = without.predict(residue, "R", "H")
    b = gating_model(False).predict(residue, "R", "H")
    assert a.gating_cost_change == pytest.approx(b.gating_cost_change, rel=1e-12)
    assert "uniform" in a.note
