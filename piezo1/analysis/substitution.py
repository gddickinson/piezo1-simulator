"""Substitution-aware spring perturbation for the elastic network.

**The problem this exists to fix.** Round 7's blind test failed with a precise
diagnostic: 99.8% of the mechanical ΔΔG's variance was *between position*, not
between substitutions. The cause is visible in the algebra rather than the data.
The old model scaled every contact of the mutated residue by one number, so

.. math::   \\Delta\\Delta G = (s - 1) \\; Q(\\text{position})

— a rank-one product in which the substitution enters only as a multiplicative
scalar. Four substitutions at R2456 therefore differ only by a factor, and rank
them all the same way. No amount of refining ``s`` can fix that; the separability
has to go.

**The fix.** Scale each contact *individually*, by properties of both the new
residue and the partner it touches. A charge change matters at contacts with
oppositely charged partners and nowhere else; proline stiffens the backbone at
sequence-local contacts; glycine removes a side chain and so weakens whatever
the side chain was touching. Different substitutions then perturb *different
subsets* of contacts, and the energy is no longer a scalar times a positional
quantity.

**What this is and is not.** An elastic network is a mechanical model, and its
spring constants are an effective stiffness standing in for packing, hydrogen
bonds, salt bridges and hydrophobic contact together. Letting that effective
stiffness depend on charge and hydrogen-bonding capacity is a statement about
what the springs represent, not a claim to have added electrostatics. It is
stated that way rather than dressed up.

Every weight is a registered parameter with bounds and a stated basis.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["ResidueProperties", "PROPERTIES", "contact_scales",
           "substitution_summary", "variance_decomposition",
           "VarianceDecomposition"]


@dataclass(frozen=True)
class ResidueProperties:
    """Side-chain properties that change what a residue's contacts are like."""

    volume: float          # A^3, Zamyatnin 1972 side-chain volumes
    charge: float          # formal charge at pH 7
    donors: int            # side-chain hydrogen-bond donors
    acceptors: int         # side-chain hydrogen-bond acceptors
    hydrophobic: float     # normalised Wimley-White, -1..1
    special: str = ""      # "P" restrains the backbone, "G" removes the side chain


#: Volumes are side-chain volumes; charges are formal at pH 7 (histidine is
#: given +0.1 because it is only partly protonated there); donor and acceptor
#: counts are side-chain only, since the backbone is common to every residue.
PROPERTIES: dict[str, ResidueProperties] = {
    "A": ResidueProperties(88.6, 0.0, 0, 0, -0.137),
    "R": ResidueProperties(173.4, +1.0, 5, 0, -0.415),
    "N": ResidueProperties(114.1, 0.0, 2, 2, -0.178),
    "D": ResidueProperties(111.1, -1.0, 0, 4, -1.000),
    "C": ResidueProperties(108.5, 0.0, 1, 0, -0.091),
    "Q": ResidueProperties(143.8, 0.0, 2, 2, -0.079),
    "E": ResidueProperties(138.4, -1.0, 0, 4, -0.668),
    "G": ResidueProperties(60.1, 0.0, 0, 0, -0.473, special="G"),
    "H": ResidueProperties(153.2, +0.1, 1, 1, -0.568),
    "I": ResidueProperties(166.7, 0.0, 0, 0, 0.336),
    "L": ResidueProperties(166.7, 0.0, 0, 0, 0.336),
    "K": ResidueProperties(168.6, +1.0, 3, 0, -0.751),
    "M": ResidueProperties(162.9, 0.0, 0, 0, 0.183),
    "F": ResidueProperties(189.9, 0.0, 0, 0, 0.241),
    "P": ResidueProperties(112.7, 0.0, 0, 0, 0.129, special="P"),
    "S": ResidueProperties(89.0, 0.0, 1, 2, -0.137),
    "T": ResidueProperties(116.1, 0.0, 1, 2, -0.046),
    "W": ResidueProperties(227.8, 0.0, 1, 0, 0.100),
    "Y": ResidueProperties(193.6, 0.0, 1, 1, -0.095),
    "V": ResidueProperties(140.0, 0.0, 0, 0, 0.220),
}


def _closeness(distance: np.ndarray) -> np.ndarray:
    """How strongly a contact feels a side-chain change.

    A C-alpha network places every contact at its C-alpha separation, but a
    side-chain substitution is felt most by whatever the side chain actually
    touches. This falls off over a side-chain length so that a 6 Å contact is
    perturbed far more than a 14 Å one — without it every contact inside the
    cutoff would respond equally, which is part of why the old model was so
    positional.
    """
    length = _P.value("substitution.contact_length")
    return np.exp(-np.maximum(distance - length, 0.0) / length)


def contact_scales(wt_aa: str, mut_aa: str, partner_aa, distances,
                   sequence_separation=None) -> np.ndarray:
    """Per-contact multiplicative spring change for one substitution.

    Returns one factor per contact rather than a single number. That is the
    whole point: a scalar makes ΔΔG separable into (substitution) × (position)
    and therefore blind to which substitution occurred.
    """
    distances = np.asarray(distances, dtype=float)
    partners = [str(p).upper()[:1] for p in partner_aa]
    n = len(distances)
    if n == 0:
        return np.zeros(0)

    wt = PROPERTIES.get(wt_aa.upper()[:1])
    mut = PROPERTIES.get(mut_aa.upper()[:1])
    if wt is None or mut is None:
        return np.ones(n)

    close = _closeness(distances)
    scale = np.ones(n)

    # --- packing: a bigger side chain presses on its neighbours -------------
    packing = _P.value("substitution.weight_volume")
    scale += packing * ((mut.volume - wt.volume) / wt.volume) * close

    # --- charge: only felt at contacts with charged partners ----------------
    charge_weight = _P.value("substitution.weight_charge")
    partner_charge = np.array([PROPERTIES[p].charge if p in PROPERTIES else 0.0
                               for p in partners])
    # A favourable pair (opposite signs) contributes stiffness; losing it
    # softens that contact, gaining one stiffens it. The product q_i·q_j is
    # negative when favourable, hence the sign.
    before = -wt.charge * partner_charge
    after = -mut.charge * partner_charge
    scale += charge_weight * (after - before) * close

    # --- hydrogen bonding: donor meets acceptor -----------------------------
    hbond_weight = _P.value("substitution.weight_hbond")
    partner_donors = np.array([PROPERTIES[p].donors if p in PROPERTIES else 0
                               for p in partners], dtype=float)
    partner_acceptors = np.array([PROPERTIES[p].acceptors if p in PROPERTIES
                                  else 0 for p in partners], dtype=float)
    pairing_before = (wt.donors * partner_acceptors
                      + wt.acceptors * partner_donors)
    pairing_after = (mut.donors * partner_acceptors
                     + mut.acceptors * partner_donors)
    denominator = np.maximum(pairing_before, 1.0)
    scale += hbond_weight * ((pairing_after - pairing_before)
                             / denominator) * close

    # --- proline: restrains the backbone, felt by sequence neighbours --------
    if mut.special == "P" and sequence_separation is not None:
        separation = np.abs(np.asarray(sequence_separation, dtype=float))
        local = separation <= _P.value("substitution.proline_span")
        scale += _P.value("substitution.weight_proline") * local

    # --- glycine: no side chain left to mediate anything ---------------------
    if mut.special == "G" and wt.special != "G":
        scale += _P.value("substitution.weight_glycine") * close

    # A spring may weaken but never invert; the network must stay positive
    # semi-definite or the quadratic form stops being an energy.
    return np.maximum(scale, _P.value("substitution.min_scale"))


def substitution_summary(wt_aa: str, mut_aa: str) -> dict:
    """What changed, for reporting alongside a prediction."""
    wt = PROPERTIES.get(wt_aa.upper()[:1])
    mut = PROPERTIES.get(mut_aa.upper()[:1])
    if wt is None or mut is None:
        return {"known": False}
    return {
        "known": True,
        "volume_change": mut.volume - wt.volume,
        "charge_change": mut.charge - wt.charge,
        "donor_change": mut.donors - wt.donors,
        "acceptor_change": mut.acceptors - wt.acceptors,
        "hydrophobicity_change": mut.hydrophobic - wt.hydrophobic,
        "introduces_proline": mut.special == "P" and wt.special != "P",
        "introduces_glycine": mut.special == "G" and wt.special != "G",
    }


# --------------------------------------------------------------------------
# The success criterion, decided in advance
# --------------------------------------------------------------------------

@dataclass
class VarianceDecomposition:
    """One-way decomposition of a score across positions and substitutions.

    The measure Round 7 diagnosed the old model with, and the criterion this
    round was set in advance: **within-position variance must exceed 20% of
    total**, or the approach has failed and that is the result.
    """

    total: float
    between: float
    within: float
    n_positions: int
    n_values: int
    n_multi: int

    @property
    def within_fraction(self) -> float:
        return self.within / self.total if self.total > 0 else 0.0

    @property
    def between_fraction(self) -> float:
        return self.between / self.total if self.total > 0 else 0.0

    def summary(self) -> str:
        return (f"{self.n_values} scores at {self.n_positions} positions "
                f"({self.n_multi} multiply substituted): "
                f"within-position {self.within_fraction:.1%}, "
                f"between-position {self.between_fraction:.1%}")


def variance_decomposition(positions, values) -> VarianceDecomposition:
    """Split the variance of ``values`` into between- and within-position parts.

    Uses the standard one-way sum-of-squares identity: total = between +
    within, where between is the spread of the position means about the grand
    mean and within is the spread inside each position.
    """
    positions = np.asarray(positions)
    values = np.asarray(values, dtype=float)
    finite = np.isfinite(values)
    positions, values = positions[finite], values[finite]
    if len(values) < 2:
        return VarianceDecomposition(0.0, 0.0, 0.0, 0, len(values), 0)

    grand = values.mean()
    total = float(((values - grand) ** 2).sum())
    between = 0.0
    within = 0.0
    unique = np.unique(positions)
    multi = 0
    for position in unique:
        group = values[positions == position]
        between += len(group) * (group.mean() - grand) ** 2
        within += float(((group - group.mean()) ** 2).sum())
        multi += len(group) > 1
    return VarianceDecomposition(total=total, between=float(between),
                                 within=float(within), n_positions=len(unique),
                                 n_values=len(values), n_multi=multi)
