"""Does the elastic network explain what half a billion years refused to change?

The census established *that* PIEZO's pore machinery is conserved and its blades
are not. It could not ask *why*, because a sequence alignment contains no
mechanics. This project has the mechanics and no evolutionary depth of its own
worth the name. Joining them is the one question neither could ask alone:

    **Is a residue's evolutionary constraint predicted by how mechanically
    coupled it is — or only by how buried it is?**

The second half of that sentence is the whole difficulty. Buried residues are
conserved in every protein ever studied, for reasons that have nothing to do
with mechanotransduction, and burial correlates with almost every mechanical
quantity the elastic network produces. A raw correlation between constraint and
gating-mode amplitude is therefore worth very little on its own. So every
feature is reported three ways: its own correlation, its correlation with burial
partialled out, and the correlation burial alone achieves.

**The null is a circular shift, not a shuffle.** Constraint is strongly
autocorrelated along the chain — neighbouring residues are in the same helix,
under the same pressure — and so is every structural feature. A permutation null
destroys that autocorrelation and is far too easy to beat: it makes almost any
comparison look significant. Shifting the whole track by a random offset and
wrapping preserves the autocorrelation exactly while destroying the
correspondence, which is the null this comparison actually needs.

Reported, never asserted: :func:`couple` returns the numbers and
:attr:`MechanicalCoupling.verdict` states what they support, including when the
answer is "nothing beyond burial".
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.family import ConstraintTrack, load_constraint
from ..core.numbering_check import piezo1_numbering
from ..core.sequence import mouse_to_human
from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .features import ResidueFeatures, build_feature_table
from .fluctuations import pearson, spearman

__all__ = ["FeatureCoupling", "MechanicalCoupling", "couple",
           "align_track_to_features", "circular_shift_null", "partial_spearman",
           "MECHANICAL_FEATURES", "BURIAL_FEATURES"]

#: The mechanical quantities tested against constraint. Each is a property of
#: the elastic network or the geometry, and none of them saw a single sequence
#: from another species — which is what makes the comparison independent.
MECHANICAL_FEATURES = ("gating_amplitude", "prs_coupling", "prs_gate_response",
                       "betweenness", "dcc_to_gate", "msf",
                       "distance_to_gate", "distance_to_axis")

#: The confound. Burial is why most conserved residues are conserved, in every
#: protein; a mechanical claim has to survive it being held fixed.
BURIAL_FEATURES = ("relative_sasa", "n_contacts")


@dataclass(frozen=True)
class FeatureCoupling:
    """One mechanical feature against constraint, with burial held fixed."""

    feature: str
    n: int
    spearman: float
    pearson: float
    partial_spearman: float
    null_mean: float
    null_sd: float
    z: float
    p_empirical: float = 1.0
    q_value: float = 1.0

    @property
    def survives_null(self) -> bool:
        """Beyond the registered margin on the circular-shift null.

        The margin is a z, but eight features are tested at once, so
        :attr:`q_value` is what the report leads with. A feature clearing the
        z and failing the correction is exactly the situation the correction
        exists for, and both are carried so it can be seen.
        """
        return abs(self.z) >= _P.value("family.min_null_z")

    @property
    def survives_correction(self) -> bool:
        return self.q_value <= _P.value("stats.alpha")

    @property
    def survives_burial(self) -> bool:
        """Keeps the registered fraction of its rank correlation past burial."""
        if self.spearman == 0:
            return False
        return (abs(self.partial_spearman)
                >= _P.value("family.burial_retention") * abs(self.spearman))


@dataclass(frozen=True)
class MechanicalCoupling:
    """Every feature's answer, plus the two controls that make them readable."""

    structure_id: str
    gene: str
    track: str
    numbering: str
    n_residues: int
    features: tuple = ()
    burial_alone: dict = field(default_factory=dict)
    own_conservation: float | None = None
    note: str = ""

    def by_name(self, name: str) -> FeatureCoupling | None:
        for f in self.features:
            if f.feature == name:
                return f
        return None

    @property
    def best(self) -> FeatureCoupling | None:
        if not self.features:
            return None
        return max(self.features, key=lambda f: abs(f.partial_spearman))

    @property
    def verdict(self) -> str:
        if not self.features:
            return "no residue could be scored; nothing to say"
        burial = max((abs(v) for v in self.burial_alone.values()), default=0.0)
        survivors = [f for f in self.features
                     if f.survives_correction and f.survives_burial]
        if not survivors:
            return (f"no mechanical quantity predicts constraint beyond burial "
                    f"(burial alone reaches rho = {burial:.2f}); on this "
                    f"structure the conserved core is explained by being "
                    f"buried, not by being coupled")
        best = max(survivors, key=lambda f: abs(f.partial_spearman))
        return (f"{best.feature} keeps rho = {best.partial_spearman:.2f} with "
                f"burial held fixed (raw {best.spearman:.2f}, shift-null "
                f"q = {best.q_value:.3f}); burial alone reaches {burial:.2f}. "
                f"{len(survivors)} of {len(self.features)} mechanical features "
                f"survive the null, the correction and the burial control")


def _ranks(values: np.ndarray) -> np.ndarray:
    order = np.argsort(values, kind="mergesort")
    ranks = np.empty(values.size, dtype=float)
    ranks[order] = np.arange(values.size, dtype=float)
    # tie-average, so a flat feature cannot manufacture a correlation
    unique, inverse, counts = np.unique(values, return_inverse=True,
                                        return_counts=True)
    if counts.max() > 1:
        sums = np.zeros(unique.size)
        np.add.at(sums, inverse, ranks)
        ranks = (sums / counts)[inverse]
    return ranks


def partial_spearman(x: np.ndarray, y: np.ndarray,
                     controls: np.ndarray) -> float:
    """Rank correlation of x and y with ``controls`` (n, k) regressed out.

    Ranks first, then least squares — so it is a partial *Spearman*, not a
    partial Pearson on raw values. The distinction matters here: relative SASA
    is bounded at zero and heavily skewed, and a linear control on the raw
    values leaves most of the burial signal in the residuals.
    """
    rx, ry = _ranks(np.asarray(x, float)), _ranks(np.asarray(y, float))
    c = np.column_stack([_ranks(np.asarray(col, float))
                         for col in np.asarray(controls, float).T])
    design = np.column_stack([np.ones(c.shape[0]), c])
    def residual(v):
        coef, *_ = np.linalg.lstsq(design, v, rcond=None)
        return v - design @ coef
    return pearson(residual(rx), residual(ry))


def align_track_to_features(features: ResidueFeatures, track: ConstraintTrack,
                            numbering: str) -> tuple[np.ndarray, np.ndarray]:
    """Constraint values lined up with the feature table's residue list.

    Returns ``(values, mask)`` — the score for each row of the table, and which
    rows got one. A mouse entry is converted through the alignment map; nothing
    is inferred from an offset.
    """
    values = np.full(features.residues.size, np.nan)
    for i, resi in enumerate(features.residues):
        human = int(resi) if numbering == "human" else mouse_to_human(int(resi))
        if human is None:
            continue
        v = track.value(human)
        if v is not None:
            values[i] = v
    return values, ~np.isnan(values)


def circular_shift_null(values: np.ndarray, feature: np.ndarray,
                        n: int = 200, seed: int = 0) -> np.ndarray:
    """Spearman under ``n`` random circular shifts of the constraint track.

    Wrapping the track along the residue axis keeps its autocorrelation exactly
    — the reason a plain permutation is the wrong null here — while destroying
    which residue each value belongs to. Returns the whole distribution, so the
    caller can take both a z and an empirical p from it; summarising to a mean
    and a standard deviation here would throw away the tail the p needs.
    """
    rng = np.random.default_rng(seed)
    size = values.size
    if size < 10:
        return np.zeros(0)
    offsets = rng.integers(1, size, size=n)
    return np.array([spearman(np.roll(values, int(off)), feature)
                     for off in offsets])


def couple(structure: Structure, track: ConstraintTrack | None = None,
           features: ResidueFeatures | None = None,
           n_null: int | None = None, seed: int = 0,
           structure_id: str = "") -> MechanicalCoupling:
    """Correlate every mechanical feature with the census constraint.

    ``features`` may be supplied to reuse a table that has already been built —
    it costs a full elastic-network solve — and is otherwise computed here with
    this project's own conservation column included, so the two evolutionary
    routes appear side by side.
    """
    track = track or load_constraint("PIEZO1")
    numbering = piezo1_numbering(structure)
    if numbering not in ("human", "mouse"):
        return MechanicalCoupling(
            structure_id=structure_id, gene=track.gene, track=track.track,
            numbering=str(numbering), n_residues=0,
            note=("refused: this entry is not in a numbering the human PIEZO1 "
                  "constraint track can be read at"))
    if features is None:
        features = build_feature_table(structure)
    n_null = int(n_null if n_null is not None else _P.value("family.null_draws"))

    values, mask = align_track_to_features(features, track, numbering)
    if mask.sum() < 20:
        return MechanicalCoupling(
            structure_id=structure_id, gene=track.gene, track=track.track,
            numbering=numbering, n_residues=int(mask.sum()),
            note="fewer than 20 residues could be scored; nothing to correlate")

    y = values[mask]
    controls = np.column_stack([
        features.columns[name][mask] for name in BURIAL_FEATURES
        if name in features.columns])
    results = []
    for name in MECHANICAL_FEATURES:
        column = features.columns.get(name)
        if column is None:
            continue
        x = column[mask]
        if np.allclose(x, x[0]):
            continue
        null = circular_shift_null(y, x, n=n_null, seed=seed)
        rho = spearman(y, x)
        null_mean = float(null.mean()) if null.size else 0.0
        null_sd = float(null.std()) if null.size else 0.0
        z = (rho - null_mean) / null_sd if null_sd > 0 else 0.0
        # (r + 1) / (n + 1), the project's convention everywhere a permutation
        # p is taken: it can never be exactly zero, which a ratio of counts can.
        extreme = int(np.count_nonzero(np.abs(null - null_mean) >= abs(rho - null_mean)))
        p_emp = (extreme + 1) / (null.size + 1) if null.size else 1.0
        results.append(FeatureCoupling(
            feature=name, n=int(mask.sum()), spearman=rho,
            pearson=pearson(y, x),
            partial_spearman=(partial_spearman(y, x, controls)
                              if controls.size else rho),
            null_mean=null_mean, null_sd=null_sd, z=float(z),
            p_empirical=float(p_emp)))

    # Eight features are tested against one track. Without a correction the
    # question "does mechanics predict constraint" gets eight chances to say
    # yes, which is the standing reason this project applies BH to any family
    # it reports together.
    if results:
        from .design import benjamini_hochberg
        bh = benjamini_hochberg([r.p_empirical for r in results],
                                names=[r.feature for r in results])
        q = dict(zip(bh.names, bh.adjusted))
        results = [FeatureCoupling(**{**r.__dict__, "q_value": float(q[r.feature])})
                   for r in results]

    burial = {name: spearman(y, features.columns[name][mask])
              for name in BURIAL_FEATURES if name in features.columns}
    own = (spearman(y, features.columns["conservation"][mask])
           if "conservation" in features.columns else None)
    return MechanicalCoupling(
        structure_id=structure_id, gene=track.gene, track=track.track,
        numbering=numbering, n_residues=int(mask.sum()),
        features=tuple(results), burial_alone=burial, own_conservation=own,
        note=("constraint from the census's deep alignment; every mechanical "
              "column from this project's elastic network, which saw no "
              "sequence from any other species"))
