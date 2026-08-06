"""Perturbation response scanning, cross-correlation and allosteric pathways."""

import numpy as np
import pytest

from piezo1.analysis.allostery import (allosteric_path, build_network,
                                       covariance_blocks_norm,
                                       cross_correlation, detour_cost,
                                       path_betweenness,
                                       perturbation_response)
from piezo1.core.annotations import load_annotations
from piezo1.physics.anm import ANM


@pytest.fixture(scope="module")
def network(human_structure):
    """An elastic network over 8YEZ with its residue and coordinate bookkeeping."""
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
    modes = anm.calc_modes(n_modes=30)
    anm.label_symmetry(modes)
    return {"coords": np.vstack(blocks), "residues": np.tile(arr, 3),
            "modes": modes, "per": len(arr)}


@pytest.fixture(scope="module")
def dcc(network):
    return cross_correlation(network["modes"])


@pytest.fixture(scope="module")
def prs(network):
    return perturbation_response(network["modes"], network["residues"])


@pytest.fixture(scope="module")
def graph(network, dcc):
    return build_network(network["coords"], dcc, contact_cutoff=10.0)


def _sites(residues, lo, hi):
    return [i for i, r in enumerate(residues) if lo <= r <= hi]


# --------------------------------------------------------------------------
# Covariance-derived matrices
# --------------------------------------------------------------------------

def test_cross_correlation_is_a_proper_correlation(dcc):
    assert np.allclose(np.diag(dcc), 1.0, atol=1e-9)
    assert np.allclose(dcc, dcc.T, atol=1e-9)
    assert dcc.min() >= -1.0 - 1e-9 and dcc.max() <= 1.0 + 1e-9


def test_anticorrelated_motion_exists(dcc):
    """A lever must have parts moving oppositely, or it is not a lever."""
    assert dcc.min() < -0.2


def test_prs_matrix_is_square_and_positive(prs, network):
    n = len(network["residues"])
    assert prs.matrix.shape == (n, n)
    assert (prs.matrix >= 0).all()
    assert prs.meta["n_sites"] == n


def test_chunking_does_not_change_the_answer(network):
    """The chunked covariance must equal an unchunked one exactly."""
    modes = network["modes"]
    small = type(modes)(eigenvalues=modes.eigenvalues[:5],
                        vectors=modes.vectors[:5, :200],
                        n_sites=200)
    a = covariance_blocks_norm(small, chunk=32)
    b = covariance_blocks_norm(small, chunk=200)
    assert np.allclose(a, b, rtol=1e-12, atol=1e-14)


def test_effectiveness_and_sensitivity_have_the_right_shape(prs, network):
    n = len(network["residues"])
    assert prs.effectiveness.shape == (n,)
    assert prs.sensitivity.shape == (n,)
    assert len(prs.top_effectors(5)) == 5
    assert len(prs.top_sensors(5)) == 5


def test_gate_response_is_dominated_by_the_pore_module(prs, network):
    """Pushing the inner helix should move the gate more than pushing the cap.

    A weak result on its own — the inner helix *is* next to the gate — but if
    it ever failed, the response matrix would be meaningless.
    """
    ann = load_annotations("human")
    res = network["residues"]
    gate = [i for i, r in enumerate(res)
            if r in set(ann.group("hydrophobic_gate").residues)]
    response = prs.per_residue(prs.response_at(gate))

    ih = np.mean([response[r] for r in response if 2432 <= r <= 2452])
    cap = np.mean([response[r] for r in response if 2198 <= r <= 2431])
    assert ih > cap


# --------------------------------------------------------------------------
# Pathways
# --------------------------------------------------------------------------

def test_path_connects_source_to_target(network, dcc):
    res = network["residues"]
    ann = load_annotations("human")
    gate = [i for i, r in enumerate(res)
            if r in set(ann.group("hydrophobic_gate").residues)]
    blade = _sites(res, 570, 1302)
    path = allosteric_path(network["coords"], dcc, blade, gate, res)
    assert len(path) > 2
    assert 570 <= path.residues[0] <= 1302
    assert path.residues[-1] in set(ann.group("hydrophobic_gate").residues)
    assert path.cost > 0
    assert len(path.correlations) == len(path) - 1


def test_anchor_lies_on_the_optimal_route(network, dcc, graph):
    """The lever model predicts the anchor transmits blade motion to the gate.

    Measured as a *detour cost*: forcing the path through the anchor should
    cost essentially nothing, meaning it is already on the best route.
    """
    res = network["residues"]
    ann = load_annotations("human")
    gate = [i for i, r in enumerate(res)
            if r in set(ann.group("hydrophobic_gate").residues)]
    blade = _sites(res, 570, 1302)

    anchor = detour_cost(graph, blade, gate, _sites(res, 2077, 2176))
    cap = detour_cost(graph, blade, gate, _sites(res, 2198, 2431))

    assert anchor["penalty"] < 0.005, "anchor should be on the optimal path"
    assert cap["penalty"] > anchor["penalty"], "the cap is not a force-transmission route"


def test_beam_is_a_near_degenerate_alternative(network, dcc, graph):
    """The beam should be a viable route, even if not the single cheapest.

    This is the honest form of the lever prediction. The beam does not appear
    on the single shortest path, but forcing the route through it costs almost
    nothing, so it is a parallel channel rather than an excluded one.
    """
    res = network["residues"]
    ann = load_annotations("human")
    gate = [i for i, r in enumerate(res)
            if r in set(ann.group("hydrophobic_gate").residues)]
    blade = _sites(res, 570, 1302)
    beam = _sites(res, 1305, 1370)
    assert beam, "beam residues must be resolved for this test to mean anything"

    d = detour_cost(graph, blade, gate, beam)
    assert d["penalty"] < 0.05, f"beam detour penalty {d['penalty']:.3f}"
    assert d["via"] >= d["direct"] - 1e-9, (
        "a constrained path can never be cheaper than the unconstrained one; "
        "if it is, the two legs are not joining at a shared via-point")


def test_detour_never_beats_the_direct_path(network, dcc, graph):
    """The invariant that caught a real error.

    Computing source-to-X and X-to-target independently lets each leg pick its
    best endpoints, which on a C3 trimer can be in different protomers. Done
    that way the detour came out *cheaper* than the shortest path, which is
    impossible.
    """
    res = network["residues"]
    blade = _sites(res, 570, 1302)
    gate = _sites(res, 2440, 2460)
    for lo, hi in ((1305, 1370), (2077, 2176), (2198, 2431), (2453, 2521)):
        via = _sites(res, lo, hi)
        if not via:
            continue
        d = detour_cost(graph, blade, gate, via)
        assert d["via"] >= d["direct"] - 1e-9, f"region {lo}-{hi}"


def test_betweenness_ranks_the_pore_module_highest(network, graph):
    """Aggregated over many paths, the anchor and CTD should dominate.

    Betweenness is used rather than a single shortest path because one
    marginally better edge reroutes a single path entirely.
    """
    res = network["residues"]
    ann = load_annotations("human")
    gate = [i for i, r in enumerate(res)
            if r in set(ann.group("hydrophobic_gate").residues)]
    blade = _sites(res, 570, 1302)
    bet = path_betweenness(graph, blade, gate, res, max_pairs=200)
    assert bet

    by_domain: dict[str, float] = {}
    for r, v in bet.items():
        d = ann.domain_at(r)
        key = d.id if d else "none"
        by_domain[key] = by_domain.get(key, 0.0) + v

    assert "anchor" in by_domain
    ranked = sorted(by_domain, key=lambda k: -by_domain[k])
    assert "anchor" in ranked[:3], f"anchor ranked {ranked.index('anchor') + 1}"
    assert by_domain.get("cap", 0.0) < by_domain["anchor"]


def test_disconnected_target_raises(network, dcc):
    res = network["residues"]
    with pytest.raises(ValueError):
        allosteric_path(network["coords"], dcc, [0], [1], res,
                        contact_cutoff=0.1)
