"""Predicted fluctuation against the deposited B-factor — the standard ANM check.

Every structure this project loads carries a ``b_factor`` for every atom, and
until Round 82 no analysis had read one. That is a conspicuous omission: the
first thing anyone does with an elastic network is ask whether its predicted
mean-square fluctuation tracks the crystallographic or cryo-EM B-factor, and
this project's central mechanism claim rests on such a network.

**The prediction already existed.** :meth:`piezo1.physics.anm.ModeSet.msf` sums
``|v_k|^2 / lambda_k`` over modes and is consumed by the feature table and by
the fluctuation colouring. What was missing was the *comparison*, not the
quantity — a distinction worth stating because it changes what this round is:
not new physics, a missing validation.

**Why the answer is more likely to be about the data than about the model.**
An observed B-factor is a temperature factor only in the sense that
refinement put it there. In a cryo-EM map it absorbs local resolution,
sharpening, per-particle scaling and the refinement's own restraints, and
several entries here were refined with **grouped** B-factors — 3JAC carries
212 distinct values over 2754 C-alphas, which is one value per thirteen
residues and cannot resolve per-residue mobility whatever the network says. So
the comparison is gated on the column first:

- a **uniform** column says nothing;
- a **grouped** column resolves fewer residues than it has;
- an **AlphaFold** model carries **pLDDT** in that field, which is a
  *confidence* and runs the other way — high where the model is certain, which
  is where a real B-factor would be low. Comparing against it would produce a
  confident negative correlation and mean nothing.

Each is refused with the reason rather than averaged into a number.

**What the correlation is worth even when the column is good.** The prediction
here is truncated: the sparse solver returns the lowest few tens of modes, and
a full B-factor prediction sums over all ``3N - 6``. Truncation over-weights
the slowest, most collective motions, so the comparison is reported against
mode count rather than at one value of it.

**And a control, because the obvious explanation is not the network.** A
residue with many neighbours moves less, and that is true of any packed solid
with no normal modes in it at all. So every correlation here is reported
beside the same correlation for **contact number** — a predictor that needs no
Hessian, no eigenvalues and no gating coordinate. If the network does not beat
it, the agreement is burial wearing a mechanism's clothes, and that is a thing
worth knowing before quoting the number.

**Measured, over the whole downloaded catalogue.** Eighteen of twenty-one
entries can answer; three cannot and say why. The network's median rank
correlation is **0.74** against the control's **0.32**, and it beats the
control on **13 of 15** entries whose column behaves like a mobility at all.
On *Pearson* the same comparison is 0.48 against 0.39 and the network wins only
**9 of 15** — so the honest summary is that the elastic network orders PIEZO1's
residues by mobility considerably better than burial does, and predicts the
*size* of the mobility barely better. Both are reported; neither is quoted
alone.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["BFactorQuality", "FluctuationComparison", "assess_b_factors",
           "observed_b_factors", "predicted_msf", "contact_number",
           "compare_fluctuations", "survey_fluctuations", "pearson",
           "spearman", "B_TO_MSF"]

#: B = (8 pi^2 / 3) <dr^2>, the definition of an isotropic temperature factor.
#: A unit conversion, not a fitted constant — and it cancels out of every
#: correlation here, so it is carried for readability rather than for the
#: answer.
B_TO_MSF = 3.0 / (8.0 * np.pi ** 2)


def pearson(a: np.ndarray, b: np.ndarray) -> float:
    """Product-moment correlation, written out rather than imported.

    Two lines, and the project's convention is that a statistic whose
    definition decides a result is visible where it is used.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    if len(a) < 2 or np.std(a) == 0.0 or np.std(b) == 0.0:
        return float("nan")
    a = a - a.mean()
    b = b - b.mean()
    return float(a @ b / np.sqrt((a @ a) * (b @ b)))


def _ranks(values: np.ndarray) -> np.ndarray:
    """Ranks with ties averaged — the convention Spearman's rho assumes.

    Ties matter here rather than being a formality: a grouped B-factor column
    is mostly ties, and ranking them arbitrarily would invent an ordering the
    refinement never asserted.
    """
    values = np.asarray(values, dtype=float)
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(len(values), dtype=float)
    ranks[order] = np.arange(len(values), dtype=float)
    sorted_values = values[order]
    start = 0
    for i in range(1, len(values) + 1):
        if i == len(values) or sorted_values[i] != sorted_values[start]:
            if i - start > 1:
                ranks[order[start:i]] = ranks[order[start:i]].mean()
            start = i
    return ranks


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Rank correlation. Reported beside Pearson because they disagree here."""
    return pearson(_ranks(a), _ranks(b))


@dataclass
class BFactorQuality:
    """Whether an entry's B-factor column can be compared with anything.

    The gate, not the result. Every field is a measurement of the column
    itself, made before any network is built, so a refusal costs nothing and
    cannot be influenced by the answer it would have given.
    """

    n_residues: int
    n_distinct: int
    minimum: float
    maximum: float
    is_confidence: bool = False
    usable: bool = True
    reason: str = ""

    @property
    def span(self) -> float:
        return self.maximum - self.minimum

    @property
    def distinct_fraction(self) -> float:
        return self.n_distinct / max(self.n_residues, 1)

    def summary(self) -> str:
        verdict = "usable" if self.usable else f"NOT usable: {self.reason}"
        return (f"{self.n_distinct} distinct over {self.n_residues} residues "
                f"({self.distinct_fraction:.0%}), {self.minimum:.1f}-"
                f"{self.maximum:.1f}; {verdict}")


@dataclass
class FluctuationComparison:
    """Predicted against observed, with everything needed to disbelieve it."""

    residues: np.ndarray
    predicted: np.ndarray            # mean-square fluctuation, arbitrary units
    observed: np.ndarray             # B-factor, A^2
    pearson_r: float
    spearman_r: float
    n_modes: int
    quality: BFactorQuality
    #: The same correlations for contact number, which needs no network at all.
    control_pearson: float = float("nan")
    control_spearman: float = float("nan")
    by_mode_count: dict = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @property
    def available(self) -> bool:
        return bool(len(self.residues))

    @property
    def beats_control(self) -> bool:
        """Whether the network explains more than burial does, on ranks."""
        if not self.available or not np.isfinite(self.control_spearman):
            return False
        return self.spearman_r > self.control_spearman

    @property
    def control_inverted(self) -> bool:
        """The control came out *negative*, which is a verdict on the column.

        More neighbours means less mobile in any packed solid. An entry whose
        B-factor rises with burial has a column that is not reporting mobility,
        whatever it is reporting, and beating the control there is not evidence
        of anything. Measured on 8YEZ and 8ZU8, where the control reaches
        -0.60 — so those two are excluded from the headline rather than
        counted as wins.
        """
        return bool(np.isfinite(self.control_spearman)
                    and self.control_spearman < 0.0)

    def summary(self) -> str:
        if not self.available:
            return f"no comparison: {self.quality.reason}"
        return (f"r = {self.pearson_r:+.3f} (Spearman {self.spearman_r:+.3f}) "
                f"over {len(self.residues)} residues, {self.n_modes} modes; "
                f"contact-number control {self.control_spearman:+.3f}")


def _is_predicted_model(structure: Structure) -> bool:
    """Whether this entry's B column holds pLDDT rather than a B-factor.

    Decided by **provenance**, not by the values: an AlphaFold model is
    identified by being one, and no threshold on a number can tell a confidence
    from a temperature factor without already assuming the answer. The
    independent confirmation that this gate points the right way is that the
    predicted models anti-correlate where the deposited ones do not, which
    ``test_fluctuations`` measures rather than assumes.
    """
    name = str(getattr(structure, "name", "") or "")
    return name.upper().startswith("AF-") or "model_v" in name


def observed_b_factors(structure: Structure) -> tuple[np.ndarray, np.ndarray]:
    """C-alpha B-factors and their residue numbers, in file order."""
    mask = structure.mask_ca()
    return structure.res_seq[mask].astype(int), structure.b_factor[mask].astype(float)


def assess_b_factors(structure: Structure) -> BFactorQuality:
    """Can this column be compared at all? Measured before anything is built."""
    residues, values = observed_b_factors(structure)
    quality = BFactorQuality(
        n_residues=len(values),
        n_distinct=int(len(np.unique(values))) if len(values) else 0,
        minimum=float(values.min()) if len(values) else 0.0,
        maximum=float(values.max()) if len(values) else 0.0,
        is_confidence=_is_predicted_model(structure))

    floor = _P.value("fluctuation.min_distinct_fraction")
    if quality.is_confidence:
        quality.usable, quality.reason = False, (
            "the B column of a predicted model holds pLDDT, a confidence that "
            "runs opposite to a temperature factor")
    elif quality.n_residues < 2:
        quality.usable, quality.reason = False, "too few residues"
    elif quality.span <= 0.0:
        quality.usable, quality.reason = False, (
            "the column is uniform, so it carries no per-residue information")
    elif quality.distinct_fraction < floor:
        quality.usable, quality.reason = False, (
            f"grouped refinement: {quality.n_distinct} distinct values over "
            f"{quality.n_residues} residues is below the "
            f"{floor:.0%} floor, so the column cannot resolve per-residue "
            f"mobility")
    return quality


def predicted_msf(modes, n_protomers: int = 3,
                  n_modes: int | None = None) -> np.ndarray:
    """Mean-square fluctuation per residue, averaged over the protomer copies.

    The three protomers are identical by construction and the observed
    B-factors are not, so averaging the prediction and comparing it against
    each copy's own B-factor would count the same prediction three times. The
    observation is averaged the same way; see :func:`compare_fluctuations`.
    """
    msf = modes.msf(n_modes)
    per = len(msf) // n_protomers
    return msf[:per * n_protomers].reshape(n_protomers, per).mean(axis=0)


def contact_number(coords: np.ndarray, cutoff: float | None = None,
                   n_protomers: int = 3) -> np.ndarray:
    """Neighbours within the network cutoff, averaged over protomers.

    The control predictor. It is the crudest possible statement of why a
    B-factor varies — a buried residue has more neighbours and moves less — and
    it uses none of the elastic network: no Hessian, no eigenvalues, no modes.
    Negated on return so that, like a fluctuation, larger means more mobile.
    """
    from scipy.spatial import cKDTree

    cutoff = _P.value("anm.cutoff") if cutoff is None else cutoff
    coords = np.asarray(coords, dtype=float)
    tree = cKDTree(coords)
    # Minus one for the residue's own entry in its own neighbourhood.
    counts = np.array([len(tree.query_ball_point(p, cutoff)) - 1
                       for p in coords], dtype=float)
    per = len(counts) // n_protomers
    return -counts[:per * n_protomers].reshape(n_protomers, per).mean(axis=0)


def compare_fluctuations(structure: Structure, blocks=None, residues=None,
                         n_modes: int | None = None,
                         cutoff: float | None = None,
                         mode_counts: tuple = (10, 20, 40)
                         ) -> FluctuationComparison:
    """Correlate the network's predicted fluctuation with the deposited B-factor.

    ``blocks`` and ``residues`` are what
    :func:`piezo1.structure.protomers.protomer_blocks` returns; they are
    accepted so a caller that already has them does not pay for them twice.

    Returns an empty comparison — not an exception — when the column cannot be
    used, because "this entry cannot answer" is a result and the survey needs
    to be able to print it.
    """
    from ..physics.anm import ANM
    from ..structure.protomers import protomer_blocks

    quality = assess_b_factors(structure)
    if not quality.usable:
        return FluctuationComparison(
            residues=np.zeros(0, dtype=int), predicted=np.zeros(0),
            observed=np.zeros(0), pearson_r=float("nan"),
            spearman_r=float("nan"), n_modes=0, quality=quality,
            meta={"structure": structure.name})

    if blocks is None or residues is None:
        blocks, residues = protomer_blocks(structure)
    if not blocks:
        quality.usable = False
        quality.reason = "needs three well-resolved protomers"
        return FluctuationComparison(
            residues=np.zeros(0, dtype=int), predicted=np.zeros(0),
            observed=np.zeros(0), pearson_r=float("nan"),
            spearman_r=float("nan"), n_modes=0, quality=quality,
            meta={"structure": structure.name})

    n_modes = int(_P.value("anm.n_modes")) if n_modes is None else n_modes
    cutoff = _P.value("anm.cutoff") if cutoff is None else cutoff
    anm = ANM.from_trimer(blocks, cutoff=cutoff).build()
    modes = anm.calc_modes(n_modes=n_modes)

    observed = _observed_per_residue(structure, residues, len(blocks))
    predicted = predicted_msf(modes, n_protomers=len(blocks))
    control = contact_number(anm.coords, cutoff, n_protomers=len(blocks))
    keep = np.isfinite(observed) & np.isfinite(predicted)

    by_count = {}
    for count in mode_counts:
        if count <= modes.n_modes:
            trial = predicted_msf(modes, len(blocks), n_modes=count)
            by_count[count] = pearson(trial[keep], observed[keep])

    return FluctuationComparison(
        residues=np.asarray(residues)[keep], predicted=predicted[keep],
        observed=observed[keep],
        pearson_r=pearson(predicted[keep], observed[keep]),
        spearman_r=spearman(predicted[keep], observed[keep]),
        n_modes=modes.n_modes, quality=quality,
        control_pearson=pearson(control[keep], observed[keep]),
        control_spearman=spearman(control[keep], observed[keep]),
        by_mode_count=by_count,
        meta={"structure": structure.name, "cutoff": cutoff,
              "n_protomers": len(blocks),
              "note": "predicted fluctuation is truncated to the modes solved "
                      "for, which over-weights the slowest; see by_mode_count"})


def _observed_per_residue(structure: Structure, residues, n_protomers: int
                          ) -> np.ndarray:
    """Mean B-factor of each shared residue number across the protomers.

    Averaged rather than taken from one chain because the prediction is
    C3-symmetric by construction: comparing a symmetric prediction against one
    protomer's own column would charge the model for an asymmetry it was never
    given the freedom to have.
    """
    numbers, values = observed_b_factors(structure)
    chains = structure.chain[structure.mask_ca()]
    out = np.full(len(residues), np.nan)
    index = {}
    for number, value, chain in zip(numbers, values, chains):
        index.setdefault(int(number), []).append(float(value))
    for i, number in enumerate(residues):
        got = index.get(int(number))
        if got:
            out[i] = float(np.mean(got[:max(n_protomers, 1)]))
    return out


def survey_fluctuations(entries=None, n_modes: int | None = None) -> list[dict]:
    """Run the comparison over the downloaded catalogue and report every row.

    Entries whose column cannot be used appear with their reason rather than
    being dropped, because how many entries could not answer is part of the
    answer.
    """
    from ..io.registry import load_registry
    from ..structure.protomers import protomer_blocks

    records = [e for e in load_registry() if e.available]
    if entries is not None:
        wanted = {e.upper() for e in entries}
        records = [e for e in records if e.pdb.upper() in wanted]

    rows = []
    for record in records:
        structure = Structure.from_file(record.path)
        quality = assess_b_factors(structure)
        row = {"pdb": record.pdb, "species": record.species,
               "state": record.state, "usable": quality.usable,
               "reason": quality.reason,
               "distinct_fraction": quality.distinct_fraction}
        if quality.usable:
            blocks, residues = protomer_blocks(structure)
            comparison = compare_fluctuations(structure, blocks, residues,
                                              n_modes=n_modes)
            row.update({"usable": comparison.available,
                        "reason": comparison.quality.reason,
                        "n_residues": len(comparison.residues),
                        "pearson": comparison.pearson_r,
                        "spearman": comparison.spearman_r,
                        "control_pearson": comparison.control_pearson,
                        "control_spearman": comparison.control_spearman,
                        "beats_control": comparison.beats_control,
                        "control_inverted": comparison.control_inverted,
                        "by_mode_count": comparison.by_mode_count})
        rows.append(row)
    return rows
