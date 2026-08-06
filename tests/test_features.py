"""The per-residue feature table.

Every check here is structural. **No phenotype comparison appears in this
file**: the Round 7 blind test returned a null result that stands as recorded,
and re-testing these features against the variant labels requires a new
pre-registration written first.
"""

import numpy as np
import pytest

from piezo1.analysis.allostery import perturbation_response
from piezo1.analysis.features import (FEATURE_NOTES, MAX_ASA,
                                      build_feature_table)
from piezo1.physics.anm import ANM


@pytest.fixture(scope="module")
def features(human_structure):
    return build_feature_table(human_structure)


def test_table_shape_and_completeness(features):
    assert len(features) > 1000
    assert len(features.names) >= 10
    for name in features.names:
        assert len(features.column(name)) == len(features)
        assert name in FEATURE_NOTES, f"{name} has no documented meaning"


def test_max_asa_table_is_complete():
    assert len(MAX_ASA) == 20
    assert MAX_ASA["GLY"] < MAX_ASA["ALA"] < MAX_ASA["TRP"]


def test_features_are_finite_where_they_should_be(features):
    for name in ("prs_gate_response", "prs_coupling", "gating_amplitude",
                 "msf", "distance_to_gate", "n_contacts"):
        values = features.column(name)
        assert np.isfinite(values).all(), name


def test_relative_sasa_is_a_fraction(features):
    values = features.column("relative_sasa")
    ok = np.isfinite(values)
    assert ok.sum() > 100
    assert values[ok].min() >= 0.0
    assert values[ok].max() <= 1.5     # a little over 1 is possible at termini


def test_gate_response_falls_off_with_distance(features):
    """A structural expectation the table must satisfy.

    Pushing a residue near the gate moves the gate more than pushing one far
    away. If this correlation ever went positive the response matrix would be
    wrong.
    """
    d = features.column("distance_to_gate")
    p = features.column("prs_gate_response")
    ok = np.isfinite(d) & np.isfinite(p)
    assert np.corrcoef(d[ok], p[ok])[0, 1] < -0.3


def test_gating_mode_used_is_symmetric(features):
    """The gating coordinate must be an A mode, not whichever came first.

    Isotropic tension is C3-symmetric, so only A modes can couple to it.
    """
    assert features.meta["gating_mode_symmetry"] == "A"


def test_conservation_reproduces_the_domain_ranking(features):
    """Joining conservation onto the table must not scramble it.

    The anchor is the most conserved domain (Round 9); if the residue join were
    off by even one position that ordering would break.
    """
    values = features.column("conservation")
    if not np.isfinite(values).any():
        pytest.skip("ortholog cache missing")
    by_domain: dict[str, list[float]] = {}
    for residue, value in zip(features.residues, values):
        if np.isfinite(value):
            by_domain.setdefault(features.domains[int(residue)], []).append(value)
    means = {k: float(np.mean(v)) for k, v in by_domain.items() if len(v) > 5}
    assert means["anchor"] > 0.95
    assert means["anchor"] >= max(means.values()) - 1e-9


def test_prs_matrix_is_symmetric(human_structure):
    """Effectiveness and sensitivity coincide on the raw matrix.

    ‖C_ij‖_F equals ‖C_ji‖_F because the Frobenius norm is transpose
    invariant, so reporting both as independent features would double-count.
    The table therefore carries one honestly named `prs_coupling` column.
    """
    st = human_structure
    chains = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() > 300:
            chains.append((st.xyz[m], st.res_seq[m]))
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:3]:
        common &= set(seq.tolist())
    arr = np.array(sorted(common))
    blocks = [x[np.searchsorted(s, arr)].astype(float) for x, s in chains[:3]]
    anm = ANM.from_trimer(blocks, cutoff=15.0).build()
    modes = anm.calc_modes(n_modes=12)

    raw = perturbation_response(modes, np.tile(arr, 3), normalise=False)
    assert raw.is_symmetric
    assert np.allclose(raw.effectiveness, raw.sensitivity, rtol=1e-9)


def test_normalisation_makes_effectiveness_degenerate(human_structure):
    """Which is why the table takes coupling from the unnormalised matrix."""
    st = human_structure
    chains = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() > 300:
            chains.append((st.xyz[m], st.res_seq[m]))
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:3]:
        common &= set(seq.tolist())
    arr = np.array(sorted(common))
    blocks = [x[np.searchsorted(s, arr)].astype(float) for x, s in chains[:3]]
    anm = ANM.from_trimer(blocks, cutoff=15.0).build()
    modes = anm.calc_modes(n_modes=12)

    normalised = perturbation_response(modes, np.tile(arr, 3), normalise=True)
    e = normalised.effectiveness
    assert e.std() / e.mean() < 0.02, "effectiveness should be near-constant"
    s = normalised.sensitivity
    assert s.std() / s.mean() > 0.2, "sensitivity should still vary"


def test_no_two_columns_are_the_same_quantity(features):
    """A near-perfect correlation means a redundant column."""
    matrix, names = features.correlations()
    for i, a in enumerate(names):
        for j, b in enumerate(names):
            if i < j and np.isfinite(matrix[i, j]):
                assert abs(matrix[i, j]) < 0.99, f"{a} and {b} are the same"


def test_lookup_and_percentile(features):
    residue = int(features.residues[len(features) // 2])
    row = features.get(residue)
    assert row["residue"] == residue
    assert set(features.names) <= set(row)
    assert features.get(999999) is None

    pct = features.percentile("prs_gate_response")
    ok = np.isfinite(pct)
    assert pct[ok].min() == pytest.approx(0.0)
    assert pct[ok].max() == pytest.approx(1.0)


def test_as_dict_matches_the_column(features):
    table = features.as_dict("conservation")
    residue = int(features.residues[10])
    assert table[residue] == pytest.approx(
        float(features.column("conservation")[10]), nan_ok=True)


def test_csv_export_has_a_row_per_residue(features):
    text = features.to_csv()
    lines = text.strip().splitlines()
    assert len(lines) == len(features) + 1
    assert lines[0].startswith("residue,domain,")


def test_metadata_carries_the_caveat(features):
    assert "predictor" in features.meta["caveat"]
    assert features.meta["n_residues"] == len(features)
