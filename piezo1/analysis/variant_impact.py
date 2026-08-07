"""Predicting a variant's mechanical effect on gating.

The question: given a mutation, does it make the gating motion **harder** or
**easier**? An elastic network gives a direct answer, because the energy needed
to traverse a known displacement is a quadratic form in the Hessian:

.. math::   E = \\tfrac{1}{2} d^{T} H d

Mutating a residue changes the springs at that residue, so the change in the
cost of the *same* observed gating motion is

.. math::   \\Delta\\Delta G_{gating} = \\tfrac{1}{2} d^{T} (H_{mut} - H_{wt}) d

A positive value means the mutation stiffens the gating coordinate — the
channel should be harder to open, i.e. loss of function. A negative value means
it softens it — easier to open, i.e. gain of function.

Two things make this cheap and clean. ``H_mut − H_wt`` is non-zero only for
contacts involving the mutated residue, so the quadratic form costs
``O(contacts)`` rather than a re-diagonalisation. And ``d`` is not invented: it
is the gating coordinate established in Round 4, where PC1 of the experimental
ensemble accounted for 90% of the variance and matched an A-symmetric mode at
0.804 overlap.

**Pre-registration.** The perturbation model below — spring constants scaled by
the fractional side-chain volume change — was chosen on physical grounds and
fixed *before* any variant's phenotype label was consulted. That matters: the
blind test in Round 7 is only blind if the predictor was not tuned against the
answer. If a different model is tried later, it must be reported as a second
hypothesis rather than substituted silently.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

__all__ = ["RESIDUE_VOLUME", "CouplingScore", "VariantImpactModel",
           "spring_scale_from_volume"]

#: Side-chain volumes in Å3 (Zamyatnin 1972), by one-letter code.
RESIDUE_VOLUME = {
    "A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8,
    "E": 138.4, "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6,
    "M": 162.9, "F": 189.9, "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8,
    "Y": 193.6, "V": 140.0,
}


def spring_scale_from_volume(wt_aa: str, mut_aa: str,
                             sensitivity: float = 1.0) -> float:
    """Multiplicative change in local spring constants for a substitution.

    A larger side chain packs harder against its neighbours and stiffens them;
    a smaller one leaves a cavity and loosens them. The scale is taken as
    ``1 + sensitivity·(V_mut − V_wt)/V_wt``, floored so a spring can never go
    negative.

    This is a deliberately simple physical model with one parameter, fixed in
    advance. Volume is used rather than a substitution matrix because the
    elastic network is a *mechanical* model: it knows about packing, not about
    evolutionary exchangeability.
    """
    v_wt = RESIDUE_VOLUME.get(wt_aa.upper())
    v_mut = RESIDUE_VOLUME.get(mut_aa.upper())
    if v_wt is None or v_mut is None:
        return 1.0
    return float(max(0.05, 1.0 + sensitivity * (v_mut - v_wt) / v_wt))


@dataclass
class CouplingScore:
    """How much a substitution changes the elastic cost of the gating motion.

    **This is not a prediction of gain or loss of function**, and the name says
    so because for four rounds it did not. The class was called
    ``CouplingScore`` and carried a ``direction`` property returning
    "stiffening (LoF-like)" or "softening (GoF-like)" — a reading that five
    pre-registered tests have failed to support, and that Rounds 47 and 54
    showed cannot be supported by any data this project could obtain.

    What the number legitimately says: this position is mechanically coupled to
    the gating motion by this much, and this substitution perturbs that
    coupling by this much. That is a statement about *where a residue sits in
    the mechanics*, which is useful and is a different claim.
    """

    residue: int
    wt_aa: str | None
    mut_aa: str | None
    gating_cost_change: float = 0.0    # arbitrary ENM energy units
    cost_change_normalised: float = 0.0   # per unit of local strain, see below
    spring_scale: float = 1.0
    n_contacts: int = 0
    local_strain: float = 0.0          # |d| at the mutated residue
    prs_gate_response: float = float("nan")
    betweenness: float = float("nan")
    dcc_to_gate: float = float("nan")
    domain: str | None = None
    modelled: bool = True
    note: str = ""

    @property
    def sign(self) -> str:
        """Whether the substitution stiffens or softens the gating motion.

        Deliberately **not** called ``direction``, and deliberately not mapped
        onto gain or loss of function. Stiffening the elastic network is not
        the same as impairing the channel: R2456H, R2456K and R2456P are
        gain-of-function and R2456C is loss, and this model gives all four
        nearly the same number.
        """
        if not np.isfinite(self.gating_cost_change) or self.gating_cost_change == 0.0:
            return "neutral"
        return "stiffening" if self.gating_cost_change > 0 else "softening"

    def as_dict(self) -> dict:
        return {k: getattr(self, k) for k in
                ("residue", "wt_aa", "mut_aa", "gating_cost_change",
                 "cost_change_normalised",
                 "spring_scale", "n_contacts", "local_strain",
                 "prs_gate_response", "betweenness", "dcc_to_gate", "domain",
                 "modelled", "note")}


@dataclass
class VariantImpactModel:
    """Scores variants by how much they change the cost of the gating motion.

    Parameters
    ----------
    coords:
        ``(n_sites, 3)`` C-alpha coordinates, protomer blocks contiguous.
    residues:
        Residue number per site.
    gating_vector:
        ``(n_sites, 3)`` displacement defining the gating coordinate. Use the
        observed curved→flat difference, or PC1 of the experimental ensemble.
    cutoff, spring, d0, gamma:
        Elastic network parameters; must match the network the gating vector
        was validated against.
    """

    coords: np.ndarray
    residues: np.ndarray
    gating_vector: np.ndarray
    cutoff: float = 15.0
    spring: str = "inverse_square"
    d0: float = 7.5
    gamma: float = 1.0
    sensitivity: float = 1.0
    #: Residue number -> one-letter code. Required for the substitution-aware
    #: perturbation: a per-contact scale needs to know what each partner IS,
    #: and without it the model falls back to the old scalar behaviour.
    sequence: dict = field(default_factory=dict)
    substitution_aware: bool = True
    meta: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.coords = np.ascontiguousarray(self.coords, dtype=np.float64)
        self.gating_vector = np.ascontiguousarray(self.gating_vector,
                                                  dtype=np.float64)
        if self.gating_vector.shape != self.coords.shape:
            raise ValueError(
                f"gating vector is {self.gating_vector.shape} but there are "
                f"{len(self.coords)} sites")
        # Normalise the gating direction so ddG is per unit of motion and does
        # not depend on how the difference vector happened to be scaled.
        norm = np.linalg.norm(self.gating_vector)
        if norm > 0:
            self.gating_vector = self.gating_vector / norm
        self._tree = cKDTree(self.coords)
        self._neighbours: dict[int, np.ndarray] = {}

    # ------------------------------------------------------------ mechanics

    def _contacts_of(self, site: int) -> np.ndarray:
        if site not in self._neighbours:
            idx = np.asarray(self._tree.query_ball_point(self.coords[site],
                                                         self.cutoff), dtype=int)
            self._neighbours[site] = idx[idx != site]
        return self._neighbours[site]

    def _spring(self, dist: np.ndarray) -> np.ndarray:
        from ..physics.anm import SPRING_MODELS
        return SPRING_MODELS[self.spring](dist, self.gamma, d0=self.d0)

    def quadratic_form_at(self, site: int, scale) -> float:
        """``½ dᵀ(H_mut − H_wt) d`` for scaling one residue's springs.

        For an anisotropic network, each contact ``(i, j)`` contributes
        ``½ k_ij [(d_i − d_j)·û_ij]²`` to the deformation energy. Scaling that
        contact's stiffness by ``s`` therefore changes the energy by
        ``½ (s − 1) k_ij [(d_i − d_j)·û_ij]²`` — only the projection of the
        relative displacement *along the contact* matters, which is exactly the
        anisotropy the model exists to capture.
        """
        neighbours = self._contacts_of(site)
        if len(neighbours) == 0:
            return 0.0
        delta = self.coords[neighbours] - self.coords[site]
        dist = np.linalg.norm(delta, axis=1)
        good = dist > 1e-9
        if not good.any():
            return 0.0
        unit = delta[good] / dist[good, None]
        k = self._spring(dist[good])
        rel = self.gating_vector[neighbours][good] - self.gating_vector[site]
        projection = np.einsum("ij,ij->i", rel, unit)
        # `scale` may be one number or one per contact. The array form is what
        # breaks the rank-one separability that made the old score blind to
        # which substitution occurred.
        factors = np.asarray(scale, dtype=float)
        if factors.ndim:
            factors = factors[good]
        return float(0.5 * np.sum((factors - 1.0) * k * projection ** 2))

    # ------------------------------------------------------------- scoring

    def _partner_context(self, site: int):
        """Partner residue letters, distances and sequence separations."""
        neighbours = self._contacts_of(site)
        if len(neighbours) == 0:
            return [], np.zeros(0), np.zeros(0)
        delta = self.coords[neighbours] - self.coords[site]
        distances = np.linalg.norm(delta, axis=1)
        here = int(self.residues[site])
        partner_numbers = self.residues[neighbours]
        letters = [self.sequence.get(int(r), "X") for r in partner_numbers]
        separation = partner_numbers.astype(float) - here
        return letters, distances, separation

    def contact_scales_at(self, site: int, wt_aa: str, mut_aa: str):
        """Per-contact spring scales for a substitution at one site."""
        from .substitution import contact_scales
        letters, distances, separation = self._partner_context(site)
        if len(distances) == 0:
            return np.ones(0)
        return contact_scales(wt_aa, mut_aa, letters, distances, separation)

    def sites_for(self, residue: int) -> np.ndarray:
        return np.flatnonzero(self.residues == residue)

    def predict(self, residue: int, wt_aa: str | None = None,
                mut_aa: str | None = None,
                spring_scale: float | None = None) -> CouplingScore:
        """Score one variant.

        The mutation is applied at *every* protomer, because a homotrimer
        carries three copies of it and the physiological channel is mutated
        throughout. Scoring one protomer would underestimate the effect
        threefold and, worse, break the C3 symmetry the gating coordinate
        depends on.
        """
        sites = self.sites_for(residue)
        if len(sites) == 0:
            return CouplingScore(
                residue=residue, wt_aa=wt_aa, mut_aa=mut_aa, modelled=False,
                note="residue not resolved in this structure")

        per_contact = (self.substitution_aware and self.sequence
                       and wt_aa and mut_aa and spring_scale is None)
        if spring_scale is None:
            spring_scale = (spring_scale_from_volume(wt_aa, mut_aa,
                                                     self.sensitivity)
                            if wt_aa and mut_aa else 1.0)

        if per_contact:
            ddg = sum(self.quadratic_form_at(
                int(s), self.contact_scales_at(int(s), wt_aa, mut_aa))
                for s in sites)
        else:
            ddg = sum(self.quadratic_form_at(int(s), spring_scale)
                      for s in sites)
        n_contacts = int(sum(len(self._contacts_of(int(s))) for s in sites))
        strain = float(np.linalg.norm(self.gating_vector[sites], axis=1).mean())

        # A residue in a rigid region can only ever score near zero, however
        # drastic the substitution. Normalising by the local strain separates
        # "this mutation is mechanically mild" from "this position barely moves".
        denom = max(strain ** 2 * max(n_contacts, 1), 1e-30)
        return CouplingScore(
            residue=residue, wt_aa=wt_aa, mut_aa=mut_aa,
            gating_cost_change=float(ddg), cost_change_normalised=float(ddg / denom),
            spring_scale=float(spring_scale), n_contacts=n_contacts,
            local_strain=strain,
            note=("per-contact perturbation" if per_contact
                  else "uniform spring scale"))

    def predict_all(self, variants, annotations=None) -> list[CouplingScore]:
        """Score an iterable of :class:`piezo1.core.annotations.Variant`."""
        out = []
        for v in variants:
            if v.residue is None:
                continue
            wt = v.wt_aa if v.wt_aa and len(v.wt_aa) == 1 else None
            mut = v.mut_aa if v.mut_aa and len(v.mut_aa) == 1 else None
            p = self.predict(v.residue, wt, mut)
            if annotations is not None:
                d = annotations.domain_at(v.residue)
                p.domain = d.id if d else None
            if wt is None or mut is None:
                p.note = ("not a single-residue substitution; spring scale "
                          "left at 1.0 so ddG is zero by construction")
            out.append(p)
        return out

    # -------------------------------------------------------- extra features

    def attach_network_features(self, predictions: list[CouplingScore],
                                prs=None, betweenness: dict | None = None,
                                dcc: np.ndarray | None = None,
                                gate_sites=None) -> list[CouplingScore]:
        """Add the allostery-derived features from Round 5, where available."""
        gate_response = None
        if prs is not None and gate_sites is not None:
            gate_response = prs.per_residue(prs.response_at(gate_sites))
        for p in predictions:
            if gate_response is not None:
                p.prs_gate_response = float(gate_response.get(p.residue, np.nan))
            if betweenness is not None:
                p.betweenness = float(betweenness.get(p.residue, 0.0))
            if dcc is not None and gate_sites is not None:
                sites = self.sites_for(p.residue)
                if len(sites):
                    p.dcc_to_gate = float(
                        np.mean(np.abs(dcc[np.ix_(sites, list(gate_sites))])))
        return predictions
