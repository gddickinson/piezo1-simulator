"""Every headline number in the documentation, with the code that regenerates it.

The project records a lot of specific numbers — a dome radius, a mode overlap, a
half-activation tension, two null results. Those live in prose, and prose does
not fail a test suite. A parameter change, a solver rewrite or a re-fetched
structure can quietly leave the documentation asserting something the code no
longer produces, and nothing would notice.

This is the registry that notices. Each :class:`Claim` names a quantity, the
value the docs state, the tolerance, where it is written, and a callable that
recomputes it from scratch. ``python scripts/reproduce.py --verify`` runs them
and reports every drift.

**A claim is not a test.** Tests assert that the code behaves; these assert that
the *documentation is still true of the code*. The two fail for different
reasons and both are worth having.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from typing import Callable

import numpy as np

__all__ = ["Claim", "ClaimResult", "CLAIMS", "verify_claims", "claims_by_cost"]


@dataclass
class Claim:
    """One documented number and how to regenerate it."""

    key: str
    description: str
    expected: float
    tolerance: float
    unit: str
    document: str
    compute: Callable[[], float]
    #: Rough runtime: "fast" (< 2 s), "medium" (< 30 s), "slow" (minutes).
    cost: str = "fast"
    #: The published value this is compared against, where one exists.
    published: str = ""
    #: Set when the claim is a *recorded result* that must never be silently
    #: revised — a null result, or a superseded number kept for the record.
    frozen: bool = False

    def relative(self, value: float) -> float:
        scale = max(abs(self.expected), 1e-12)
        return abs(value - self.expected) / scale


@dataclass
class ClaimResult:
    claim: Claim
    value: float | None
    seconds: float
    error: str = ""
    #: Parameter overrides in force when this was computed. Non-empty means the
    #: number did not come from the documented parameter set and must not be
    #: compared with the documentation as though it had.
    parameters: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        if self.error or self.value is None:
            return False
        return abs(self.value - self.claim.expected) <= self.claim.tolerance

    @property
    def comparable(self) -> bool:
        """False when overrides mean this cannot be read as reproduction."""
        return not self.parameters

    def line(self) -> str:
        if self.error:
            return f"  SKIP {self.claim.key:28s} {self.error}"
        mark = "ok  " if self.ok else "DRIFT"
        return (f"  {mark} {self.claim.key:28s} {self.value:12.4f} "
                f"vs {self.claim.expected:10.4f} {self.claim.unit:8s} "
                f"({self.seconds:5.1f} s)")


# --------------------------------------------------------------------------
# The computations
# --------------------------------------------------------------------------

def _structure(pdb: str):
    from ..config import STRUCTURE_DIR
    from ..core.structure import Structure
    path = STRUCTURE_DIR / f"{pdb.upper()}.cif"
    if not path.exists():
        raise FileNotFoundError(f"{pdb} not downloaded")
    return Structure.from_file(path)


def _blocks(structure):
    from ..ui.model_utils import protomer_blocks
    return protomer_blocks(structure)


def _tm_surface(structure, species: str) -> np.ndarray:
    import json

    from ..config import RESOURCE_DIR
    tms = json.loads((RESOURCE_DIR / f"uniprot_{species}.json").read_text())
    points = []
    for chain in structure.chains:
        mask = structure.mask_ca() & (structure.chain == chain)
        if mask.sum() < 300:
            continue
        xyz, seq = structure.xyz[mask], structure.res_seq[mask]
        for tm in tms["transmembrane"]:
            mid = 0.5 * (tm["start"] + tm["end"])
            half = max(2.0, (tm["end"] - tm["start"]) / 6.0)
            sel = (seq >= mid - half) & (seq <= mid + half)
            if sel.sum() >= 3:
                points.append(xyz[sel].mean(axis=0))
    return np.array(points)


def _dome_radius() -> float:
    from ..structure.geometry import measure_dome
    st = _structure("7WLT")
    blocks, _ = _blocks(st)
    return measure_dome(blocks, _tm_surface(st, "mouse")).radius_of_curvature / 10.0


def _gating_overlap() -> float:
    """Overlap of the lowest symmetric mode with the observed transition."""
    from ..physics.anm import ANM
    from ..structure.superpose import match_protomers

    curved, flat = _structure("7WLT"), _structure("7WLU")
    cb, cr = _blocks(curved)
    fb, fr = _blocks(flat)
    common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))

    def resample(st, residues):
        out = []
        for chain in st.chains:
            mask = st.mask_ca() & (st.chain == chain)
            if mask.sum() < 300:
                continue
            index = {int(r): i for i, r in enumerate(st.res_seq[mask])}
            xyz = st.xyz[mask]
            out.append(np.array([xyz[index[r]] for r in residues], dtype=float))
        return out[:3]

    cb, fb = resample(curved, common), resample(flat, common)
    match = match_protomers(cb, fb)
    fb = [fb[i] for i in match.order]

    anm = ANM.from_trimer(cb, cutoff=15.0, spring="inverse_square").build()
    modes = anm.calc_modes(n_modes=30)
    anm.label_symmetry(modes)

    from ..structure.superpose import kabsch
    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    fitted = (np.vstack(fb) - centroid) @ rotation.T + translation
    displacement = (fitted - np.vstack(cb)).ravel()

    # `overlap` returns the whole spectrum at once; take the best A mode.
    overlaps = np.abs(np.asarray(modes.overlap(displacement), dtype=float))
    symmetric = np.array([s == "A" for s in modes.symmetry])
    return float(overlaps[symmetric].max()) if symmetric.any() else 0.0


def _t50_kinetics() -> float:
    from ..physics.kinetics import GatingModel
    return float(GatingModel().half_activation())


def _membrane_tension() -> float:
    from ..physics.membrane import MembraneParameters, kt_per_nm2_to_mnm
    return kt_per_nm2_to_mnm(MembraneParameters(kappa=20.0,
                                                tension=20.0 / 14.0 ** 2).tension)


def _decay_length_recovered() -> float:
    from ..physics.membrane import MembraneParameters, solve_footprint
    params = MembraneParameters(kappa=20.0, tension=20.0 / 14.0 ** 2)
    return float(solve_footprint(10.0, 0.2, params).fitted_decay_length())


def _elastica_energy() -> float:
    from ..physics.elastica import solve_elastica
    from ..physics.membrane import MembraneParameters
    return float(solve_elastica(8.691, 1.992, MembraneParameters()).energy)


def _elastica_area() -> float:
    from ..physics.elastica import solve_elastica
    from ..physics.membrane import MembraneParameters
    return float(solve_elastica(8.691, 1.992, MembraneParameters()).excess_area)


def _linear_overestimate() -> float:
    from ..physics.elastica import compare_with_linear
    from ..physics.membrane import MembraneParameters
    return 1.0 / compare_with_linear(8.691, 1.992, MembraneParameters()).energy_ratio


def _pore_bottleneck(pdb: str) -> float:
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    st = _structure(pdb)
    blocks, _ = _blocks(st)
    return float(pore_profile(st, detect_c3_axis(blocks), step=1.0)
                 .bottleneck_radius)


def _wetting_score(pdb: str) -> float:
    from .hydration import load_grid, predict_wetting
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    grid = load_grid()
    if not grid.available:
        raise FileNotFoundError("CHAP grid not downloaded")
    st = _structure(pdb)
    blocks, _ = _blocks(st)
    profile = pore_profile(st, detect_c3_axis(blocks), step=1.0)
    return float(predict_wetting(st, profile, grid).score)


def _cds_identity() -> float:
    from ..core.sequences import compare_sequences, load_named_sequences
    seqs = {s.key: s for s in load_named_sequences()}
    if "cds_human" not in seqs:
        raise FileNotFoundError("human CDS not downloaded")
    return compare_sequences(seqs["uniprot_human"], seqs["cds_human"],
                             "positional").identity


def _recorded(path: str, *keys) -> float:
    import json

    from ..config import DERIVED_DIR
    file = DERIVED_DIR / path
    if not file.exists():
        raise FileNotFoundError(f"{path} not generated")
    node = json.loads(file.read_text())
    for key in keys:
        node = node[key]
    return float(node)


# --------------------------------------------------------------------------
# The registry
# --------------------------------------------------------------------------

CLAIMS: list[Claim] = [
    Claim("dome.radius_7wlt", "Dome radius of curvature, 7WLT",
          9.72, 0.15, "nm", "docs/SCIENCE.md", _dome_radius, "medium",
          published="10.2 nm (Haselwandter & MacKinnon 2018)"),
    Claim("anm.gating_overlap", "Lowest A-mode overlap with the transition",
          0.705, 0.05, "", "docs/SCIENCE.md", _gating_overlap, "slow"),
    Claim("kinetics.t50", "Emergent half-activation tension",
          2.71, 0.05, "mN/m", "docs/SCIENCE.md", _t50_kinetics, "fast",
          published="2.7 +- 0.1 mN/m (Lewis & Grandl 2015)"),
    Claim("membrane.tension", "Tension implied by kappa=20 kT, lambda=14 nm",
          0.420, 0.005, "mN/m", "docs/SCIENCE.md", _membrane_tension, "fast",
          published="0.42 mN/m"),
    Claim("membrane.lambda_recovered", "Decay length recovered from the solver",
          13.998, 0.05, "nm", "docs/SCIENCE.md", _decay_length_recovered,
          "fast", published="14.0 nm by construction"),
    Claim("elastica.energy", "Nonlinear footprint energy at the 7WLT geometry",
          25.27, 0.3, "kT", "docs/SCIENCE.md", _elastica_energy, "medium"),
    Claim("elastica.excess_area", "Nonlinear footprint excess area",
          178.9, 2.0, "nm^2", "docs/SCIENCE.md", _elastica_area, "medium"),
    Claim("elastica.linear_overestimate",
          "Factor by which linear theory overestimates the footprint",
          3.65, 0.06, "x", "docs/SCIENCE.md", _linear_overestimate, "medium"),
    Claim("pore.bottleneck_8yez", "Closed-state pore bottleneck",
          0.95, 0.05, "A", "docs/SCIENCE.md",
          lambda: _pore_bottleneck("8YEZ"), "medium"),
    Claim("pore.bottleneck_11zc", "Flat-state pore bottleneck",
          3.30, 0.10, "A", "docs/SCIENCE.md",
          lambda: _pore_bottleneck("11ZC"), "medium"),
    Claim("hydration.score_8yez", "Rao 2019 wetting score, closed 8YEZ",
          0.82, 0.05, "", "docs/SCIENCE.md",
          lambda: _wetting_score("8YEZ"), "medium"),
    Claim("hydration.score_11zc", "Rao 2019 wetting score, flat 11ZC",
          0.00, 0.01, "", "docs/SCIENCE.md",
          lambda: _wetting_score("11ZC"), "medium"),
    Claim("sequence.cds_identity",
          "Human CDS translation against UniProt Q92508",
          1.0, 1e-9, "fraction", "SESSION_LOG.md", _cds_identity, "fast"),

    # --- recorded results, never to be silently revised --------------------
    Claim("round7.p_value", "Round 7 blind test p-value",
          0.234, 0.02, "", "docs/VALIDATION.md",
          lambda: _recorded("validation_round7.json", "primary",
                            "permutation", "p_value"),
          "fast", frozen=True),
    Claim("round7.auroc", "Round 7 AUROC",
          0.5417, 0.001, "", "docs/VALIDATION.md",
          lambda: _recorded("validation_round7.json", "primary", "auroc"),
          "fast", frozen=True),
    Claim("round22.cliffs_delta", "Round 22 primary effect size",
          -0.2105, 0.001, "", "docs/VALIDATION_ROUND22.md",
          lambda: _recorded("validation_round22.json", "primary",
                            "cliffs_delta"),
          "fast", frozen=True),
    Claim("round22.auroc", "Round 22 primary AUROC",
          0.3947, 0.001, "", "docs/VALIDATION_ROUND22.md",
          lambda: _recorded("validation_round22.json", "primary", "auroc"),
          "fast", frozen=True),
]


def claims_by_cost(max_cost: str = "slow") -> list[Claim]:
    """Claims no more expensive than ``max_cost``."""
    order = {"fast": 0, "medium": 1, "slow": 2}
    ceiling = order[max_cost]
    return [c for c in CLAIMS if order[c.cost] <= ceiling]


def verify_claims(claims: list[Claim] | None = None,
                  verbose: bool = True,
                  allow_overrides: bool = False) -> list[ClaimResult]:
    """Recompute each claim and compare with what the documentation says.

    **Refuses to run against overridden parameters.** Every documented number
    was produced with the registry at its defaults, so recomputing with a
    changed value would report drift that the user caused — and the obvious
    reading of that report would be that the code is broken. Pass
    ``allow_overrides=True`` to compare deliberately, in which case every
    result is annotated with the parameter set that produced it.
    """
    from ..parameters import PARAMETERS

    changed = PARAMETERS.overrides()
    if changed and not allow_overrides:
        raise RuntimeError(
            "cannot verify documented numbers against modified parameters — "
            + PARAMETERS.override_summary()
            + ". Reset them (piezo1.parameters.reset()) or pass "
              "allow_overrides=True to compare on purpose.")

    results = []
    for claim in (claims if claims is not None else CLAIMS):
        start = time.time()
        try:
            value = float(claim.compute())
            error = ""
        except Exception as exc:
            value, error = None, f"{type(exc).__name__}: {exc}"
        result = ClaimResult(claim, value, time.time() - start, error,
                             parameters=dict(changed))
        results.append(result)
        if verbose:
            print(result.line(), flush=True)
    return results
