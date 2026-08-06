"""Principal component analysis of the experimental structure ensemble.

The deposited PIEZO structures are not independent snapshots of a simulation —
they are the conformations that experiment has actually managed to trap. The
principal components of that set are therefore an empirical map of the
conformational space the protein is known to explore, and comparing them with
the elastic-network modes is the strongest available test of whether the
network model describes real motion or merely plausible motion.

Four traps have to be handled before any of that is meaningful, and each of
them silently produces a number rather than an error:

1. **Species.** Human structures are numbered by Q92508, mouse by E2JF22, and
   the offset between them is not constant. Everything is converted to human
   numbering through :mod:`piezo1.core.sequence` before intersecting.
2. **Coverage.** Entries resolve different residues. Intersecting all 20 usable
   PIEZO1 entries leaves only 325 residues, because two poorly-ordered early
   structures drag the whole set down. :func:`build_ensemble` reports what each
   entry costs and can drop the worst offenders.
3. **Protomer correspondence.** Deposited chain labels do not reliably indicate
   rotational order around the three-fold axis, so correspondence is
   established by superposition, not by label.
4. **Paralogues.** PDB 6KG7 is *PIEZO2*, not PIEZO1. It is excluded by default;
   including it would put a 40%-identity paralogue into an ensemble meant to
   describe one protein's motion.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..core.sequence import load_numbering_map
from ..core.structure import Structure
from ..io.registry import StructureRecord, load_registry
from ..structure.superpose import match_protomers, superpose

__all__ = ["EnsembleMember", "StructureEnsemble", "PCAResult",
           "build_ensemble", "subspace_overlap", "rwsip"]

#: Entries excluded from ensemble analysis by default, with the reason. Pass
#: an explicit ``exclude`` mapping to :func:`build_ensemble` to override.
DEFAULT_EXCLUSIONS = {
    "6KG7": "PIEZO2, a paralogue - not the same protein",
    "3JAC": "large poly-UNK regions with arbitrary residue numbering",
    "6LQI": ("Piezo1.1 splice isoform, missing residues 1382-1405. Its "
             "difference from the rest is a sequence difference, not a "
             "conformational one, and it dominates a whole principal "
             "component on its own: including it splits the gating "
             "coordinate across PC1 (58%) and PC2 (36%), where excluding it "
             "gives a single clean PC1 at 90%"),
}


@dataclass
class EnsembleMember:
    pdb: str
    species: str
    state: str
    coords: np.ndarray            # (n_sites, 3), superposed, human numbering
    protomer_order: tuple[int, ...]
    handedness_flipped: bool
    rmsd_to_reference: float
    resolution: float | None = None


@dataclass
class StructureEnsemble:
    """A set of structures on a common residue basis, mutually superposed."""

    members: list[EnsembleMember]
    residues: np.ndarray          # human residue numbers, per protomer
    n_protomers: int = 3
    reference: str = ""
    excluded: dict[str, str] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.members)

    @property
    def n_sites(self) -> int:
        return len(self.residues) * self.n_protomers

    @property
    def coords(self) -> np.ndarray:
        """``(n_structures, n_sites, 3)``."""
        return np.stack([m.coords for m in self.members])

    @property
    def labels(self) -> list[str]:
        return [m.pdb for m in self.members]

    def states(self) -> list[str]:
        return [m.state for m in self.members]

    def pairwise_rmsd(self) -> np.ndarray:
        x = self.coords
        n = len(x)
        out = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                d = x[i] - x[j]
                out[i, j] = out[j, i] = np.sqrt((d * d).sum(axis=1).mean())
        return out

    # ------------------------------------------------------------------ PCA

    def pca(self, n_components: int | None = None) -> "PCAResult":
        """Principal components of the ensemble, by SVD of the mean-centred set.

        SVD of the ``(n_structures, 3N)`` matrix is used rather than forming the
        ``3N x 3N`` covariance: with 3N in the tens of thousands and only ~20
        structures, the covariance is enormous and rank-deficient, while the
        SVD costs almost nothing and is numerically better behaved.
        """
        x = self.coords.reshape(len(self), -1)
        mean = x.mean(axis=0)
        centred = x - mean
        # Full covariance eigenvalues are (singular values^2)/(n-1).
        u, s, vt = np.linalg.svd(centred, full_matrices=False)
        var = s ** 2 / max(len(self) - 1, 1)
        k = len(s) if n_components is None else min(n_components, len(s))
        return PCAResult(
            eigenvalues=var[:k],
            components=vt[:k].reshape(k, -1, 3),
            mean=mean.reshape(-1, 3),
            projections=(u[:, :k] * s[:k]),
            total_variance=float(var.sum()),
            labels=list(self.labels),
            states=self.states(),
        )


@dataclass
class PCAResult:
    eigenvalues: np.ndarray       # variance along each component
    components: np.ndarray        # (n_comp, n_sites, 3)
    mean: np.ndarray              # (n_sites, 3)
    projections: np.ndarray       # (n_structures, n_comp)
    total_variance: float
    labels: list[str] = field(default_factory=list)
    states: list[str] = field(default_factory=list)

    @property
    def n_components(self) -> int:
        return len(self.eigenvalues)

    @property
    def variance_explained(self) -> np.ndarray:
        return self.eigenvalues / max(self.total_variance, 1e-30)

    def cumulative_variance(self) -> np.ndarray:
        return np.cumsum(self.variance_explained)

    def amplitude(self, index: int) -> float:
        """Norm of the displacement along a component, in Angstrom.

        This is a length in the full 3N-dimensional space, so for a
        thousands-of-sites model it is a large number and does *not* mean any
        atom moves that far. Use :meth:`rms_displacement` for the per-site
        figure, which is the one to quote.
        """
        return float(np.sqrt(self.eigenvalues[index]))

    def rms_displacement(self, index: int) -> float:
        """Per-site RMS displacement along a component, in Angstrom.

        The interpretable version of :meth:`amplitude`: divide out the number
        of sites so the number means what a reader assumes it means.
        """
        return float(np.sqrt(self.eigenvalues[index] / len(self.components[index])))

    def collectivity(self, index: int) -> float:
        a = (self.components[index] ** 2).sum(axis=1)
        a = a / a.sum()
        nz = a[a > 0]
        return float(np.exp(-np.sum(nz * np.log(nz))) / len(a))

    # ------------------------------------------------- comparison with modes

    def overlap_with_modes(self, modes, n_pcs: int | None = None,
                           n_modes: int | None = None) -> np.ndarray:
        """``|cos angle|`` between each principal component and each ANM mode."""
        npc = self.n_components if n_pcs is None else min(n_pcs, self.n_components)
        nm = modes.n_modes if n_modes is None else min(n_modes, modes.n_modes)
        a = self.components[:npc].reshape(npc, -1)
        b = modes.vectors[:nm].reshape(nm, -1)
        if a.shape[1] != b.shape[1]:
            raise ValueError(
                f"principal components span {a.shape[1] // 3} sites but the "
                f"mode set spans {b.shape[1] // 3}. Build the elastic network "
                f"on the ensemble's own residue basis.")
        a = a / np.linalg.norm(a, axis=1, keepdims=True)
        b = b / np.linalg.norm(b, axis=1, keepdims=True)
        return np.abs(a @ b.T)

    def best_mode_for(self, modes, pc: int = 0) -> tuple[int, float]:
        row = self.overlap_with_modes(modes, n_pcs=pc + 1)[pc]
        j = int(np.argmax(row))
        return j, float(row[j])

    def cumulative_overlap_with_modes(self, modes, pc: int = 0,
                                      n_modes: int | None = None) -> np.ndarray:
        row = self.overlap_with_modes(modes, n_pcs=pc + 1, n_modes=n_modes)[pc]
        return np.sqrt(np.cumsum(row ** 2))


# --------------------------------------------------------------------------
# Subspace comparison
# --------------------------------------------------------------------------

def subspace_overlap(a_vectors: np.ndarray, b_vectors: np.ndarray) -> float:
    """Unweighted subspace overlap of two mode sets, in [0, 1].

    ``(1/D) Σ_ij (u_i·v_j)²`` for D vectors each — 1 when the two subspaces
    coincide, ~D/3N when they are unrelated.
    """
    a = a_vectors.reshape(len(a_vectors), -1)
    b = b_vectors.reshape(len(b_vectors), -1)
    a = a / np.linalg.norm(a, axis=1, keepdims=True)
    b = b / np.linalg.norm(b, axis=1, keepdims=True)
    return float(((a @ b.T) ** 2).sum() / len(a))


def rwsip(a_vectors: np.ndarray, a_weights: np.ndarray,
          b_vectors: np.ndarray, b_weights: np.ndarray,
          n: int | None = None) -> float:
    """Root weighted square inner product (Carnevale et al. 2007).

    Like :func:`subspace_overlap` but weighting each pair by the amplitude the
    two models assign to those directions, so a mode that neither model
    actually populates cannot inflate the score:

    .. math::

        RWSIP = \\left[\\frac{\\sum_{ij}\\sqrt{w^A_i w^B_j}(u_i\\cdot v_j)^2}
                            {\\sum_i \\sqrt{w^A_i w^B_i}}\\right]^{1/2}

    For PCA the weight is the variance; for an elastic network it is ``1/λ``,
    since mean-square fluctuation along a mode scales that way. Passing raw ANM
    eigenvalues instead would weight the *stiffest* modes most, which is
    backwards.
    """
    k = min(len(a_vectors), len(b_vectors)) if n is None else n
    a = a_vectors[:k].reshape(k, -1)
    b = b_vectors[:k].reshape(k, -1)
    a = a / np.linalg.norm(a, axis=1, keepdims=True)
    b = b / np.linalg.norm(b, axis=1, keepdims=True)
    wa = np.asarray(a_weights, float)[:k]
    wb = np.asarray(b_weights, float)[:k]
    dots = (a @ b.T) ** 2
    weight = np.sqrt(np.outer(wa, wb))
    denom = np.sqrt(wa * wb).sum()
    return float(np.sqrt((weight * dots).sum() / max(denom, 1e-30)))


# --------------------------------------------------------------------------
# Construction
# --------------------------------------------------------------------------

def _protomer_sets(st: Structure, min_ca: int = 300):
    out = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() > min_ca:
            out.append((st.xyz[m].astype(np.float64), st.res_seq[m]))
    return out[:3]


def _to_human(residues: np.ndarray, species: str, nm) -> dict[int, int]:
    """Map this structure's residue numbers to human numbering."""
    if species == "human":
        return {int(r): int(r) for r in residues}
    out = {}
    for r in residues:
        h = nm.to_a(int(r))
        if h is not None:
            out[int(r)] = h
    return out


def build_ensemble(records: list[StructureRecord] | None = None,
                   species: str | None = None,
                   exclude: dict[str, str] | None = None,
                   min_common: int = 0,
                   reference: str | None = None,
                   verbose: bool = False) -> StructureEnsemble:
    """Assemble structures onto a common residue basis in human numbering.

    Parameters
    ----------
    species:
        Restrict to ``"human"`` or ``"mouse"``. Mixing them is allowed and the
        numbering is handled, but the sequences genuinely differ at ~18% of
        positions, so some of the resulting variance is evolutionary rather
        than conformational.
    min_common:
        Drop any entry that would reduce the shared residue count below this.
        Entries are considered worst-coverage-first, and every drop is recorded
        in :attr:`StructureEnsemble.excluded`.
    """
    nm = load_numbering_map()
    reg = load_registry()
    recs = records if records is not None else [
        r for r in reg.available() if r.state != "fragment"]
    if species:
        recs = [r for r in recs if r.species == species]

    excluded = dict(DEFAULT_EXCLUSIONS if exclude is None else exclude)
    recs = [r for r in recs if r.pdb not in excluded]

    loaded: list[tuple[StructureRecord, list, dict[int, int]]] = []
    for rec in recs:
        st = Structure.from_file(rec.path)
        blocks = _protomer_sets(st)
        if len(blocks) < 3:
            excluded[rec.pdb] = f"only {len(blocks)} well-resolved protomers"
            continue
        shared = set(blocks[0][1].tolist())
        for _, seq in blocks[1:]:
            shared &= set(seq.tolist())
        mapping = _to_human(np.array(sorted(shared)), rec.numbering_species, nm)
        loaded.append((rec, blocks, mapping))

    if len(loaded) < 3:
        raise ValueError(f"only {len(loaded)} usable structures")

    def common_of(items) -> set[int]:
        return set.intersection(*[set(m.values()) for _, _, m in items])

    if min_common:
        # Greedily drop whichever entry is costing the most coverage.
        while len(loaded) > 3 and len(common_of(loaded)) < min_common:
            costs = []
            for i in range(len(loaded)):
                rest = loaded[:i] + loaded[i + 1:]
                costs.append((len(common_of(rest)), i))
            gain, idx = max(costs)
            if gain <= len(common_of(loaded)):
                break
            rec = loaded[idx][0]
            excluded[rec.pdb] = (f"restricted the shared basis to "
                                 f"{len(common_of(loaded))} residues")
            loaded.pop(idx)

    common = np.array(sorted(common_of(loaded)))
    if len(common) < 30:
        raise ValueError(f"only {len(common)} residues shared across the set")

    # Reference: prefer the requested entry, else the best-resolution one.
    order = sorted(range(len(loaded)),
                   key=lambda i: (loaded[i][0].pdb != (reference or ""),
                                  loaded[i][0].resolution or 99))
    loaded = [loaded[i] for i in order]
    ref_rec = loaded[0][0]

    members: list[EnsembleMember] = []
    ref_blocks = None
    for rec, blocks, mapping in loaded:
        inverse = {h: own for own, h in mapping.items()}
        own_numbers = np.array([inverse[h] for h in common])
        picked = []
        for xyz, seq in blocks:
            idx = np.searchsorted(seq, own_numbers)
            picked.append(xyz[idx])
        if ref_blocks is None:
            ref_blocks = picked
            coords = np.vstack(picked)
            match_order, flipped, err = (0, 1, 2), False, 0.0
        else:
            match = match_protomers(picked, ref_blocks)
            match_order, flipped = match.order, match.handedness_flipped
            stacked = np.vstack([picked[i] for i in match.order])
            coords, err = superpose(stacked, np.vstack(ref_blocks))
        members.append(EnsembleMember(
            pdb=rec.pdb, species=rec.species, state=rec.state, coords=coords,
            protomer_order=tuple(match_order), handedness_flipped=flipped,
            rmsd_to_reference=float(err), resolution=rec.resolution))
        if verbose:
            print(f"  {rec.pdb:6s} {rec.state:13s} rmsd {err:6.2f} A"
                  f"{'  (protomers reversed)' if flipped else ''}")

    return StructureEnsemble(
        members=members, residues=common, reference=ref_rec.pdb,
        excluded=excluded,
        meta={"n_structures": len(members), "n_residues": len(common),
              "species_filter": species,
              "numbering": "human Q92508"})
