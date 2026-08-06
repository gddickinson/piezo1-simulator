"""Anisotropic network model (ANM) for large membrane proteins.

An ANM treats the protein as a network of beads (one per C-alpha) joined by
harmonic springs between neighbours. Diagonalising the resulting Hessian gives
the normal modes: the low-frequency ones describe the large, collective,
functionally relevant motions. For PIEZO1 this is exactly the right tool —
the gating motion is a whole-blade flattening that no all-atom simulation can
reach in accessible time, but that falls out of an ANM in seconds.

Two things here are specific to PIEZO1 and worth noting:

**Sparse everything.** A PIEZO1 trimer is ~3900 resolved C-alphas (7500 for a
full-length model), so the Hessian is 11700-22500 square. Dense
diagonalisation is wasteful when only the lowest ~40 modes matter; we build a
sparse Hessian and use shift-invert Lanczos.

**Symmetry labelling.** The channel is a C3 homotrimer, so every mode carries an
irreducible-representation label: ``A`` (three-fold symmetric) or ``E``
(the degenerate pair). This is not decoration. Isotropic membrane tension is
itself C3-symmetric, so *only* ``A`` modes can couple to it at first order.
Identifying the lowest ``A`` mode identifies the candidate gating coordinate.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.spatial import cKDTree
from ..parameters import PARAMETERS as _P

__all__ = ["ANM", "ModeSet", "build_hessian", "SPRING_MODELS"]


# --------------------------------------------------------------------------
# Spring models
# --------------------------------------------------------------------------

def _gamma_uniform(dist: np.ndarray, gamma: float, **_: float) -> np.ndarray:
    """Classical Hinsen/Atilgan ANM: every contact has the same stiffness."""
    return np.full_like(dist, gamma)


def _gamma_inverse_square(dist: np.ndarray, gamma: float,
                          d0: float = field(default_factory=lambda: _P.value("anm.d0")), **_: float) -> np.ndarray:
    """Parameter-free ANM weighting, gamma * (d0/d)^2 (Yang et al. 2009).

    Distant contacts are softened, which suppresses the spurious stiffening
    that a hard cutoff produces at the edge of a large, hollow assembly like
    the PIEZO1 blade.
    """
    return gamma * (d0 / np.maximum(dist, 1e-6)) ** 2


def _gamma_inverse_sixth(dist: np.ndarray, gamma: float,
                         d0: float = 7.5, **_: float) -> np.ndarray:
    """Steeper 1/d^6 weighting, closer to a van der Waals falloff."""
    return gamma * (d0 / np.maximum(dist, 1e-6)) ** 6


SPRING_MODELS = {
    "uniform": _gamma_uniform,
    "inverse_square": _gamma_inverse_square,
    "inverse_sixth": _gamma_inverse_sixth,
}


# --------------------------------------------------------------------------
# Hessian construction
# --------------------------------------------------------------------------

def build_hessian(coords: np.ndarray, cutoff: float = 15.0, gamma: float = 1.0,
                  spring: str = "inverse_square", d0: float = 7.5) -> sp.csr_matrix:
    """Assemble the sparse 3N x 3N ANM Hessian.

    For a pair ``(i, j)`` separated by ``r`` with unit vector ``u``, the 3x3
    off-diagonal super-element is ``-gamma_ij * outer(u, u)``; diagonal
    super-elements are minus the sum of their row, which enforces the six
    zero-frequency rigid-body modes exactly.
    """
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    n = len(coords)
    if n < 4:
        raise ValueError("need at least four sites for an ANM")
    if spring not in SPRING_MODELS:
        raise ValueError(f"unknown spring model {spring!r}; "
                         f"choose from {sorted(SPRING_MODELS)}")

    tree = cKDTree(coords)
    pairs = np.asarray(list(tree.query_pairs(cutoff)), dtype=np.int64)
    if len(pairs) == 0:
        raise ValueError(f"no contacts within {cutoff} A — cutoff too small")

    i, j = pairs[:, 0], pairs[:, 1]
    diff = coords[j] - coords[i]
    dist = np.linalg.norm(diff, axis=1)
    unit = diff / dist[:, None]
    k = SPRING_MODELS[spring](dist, gamma, d0=d0)

    # Outer products for every pair: (n_pairs, 3, 3)
    blocks = -k[:, None, None] * unit[:, :, None] * unit[:, None, :]

    rows = np.empty(len(pairs) * 9 * 2, dtype=np.int64)
    cols = np.empty_like(rows)
    vals = np.empty(len(rows), dtype=np.float64)

    a, b = np.meshgrid(np.arange(3), np.arange(3), indexing="ij")
    a, b = a.ravel(), b.ravel()
    m = len(pairs)
    # Upper block (i, j)
    rows[:m * 9] = (3 * i[:, None] + a[None, :]).ravel()
    cols[:m * 9] = (3 * j[:, None] + b[None, :]).ravel()
    vals[:m * 9] = blocks.reshape(m, 9).ravel()
    # Lower block (j, i) — the Hessian is symmetric.
    rows[m * 9:] = (3 * j[:, None] + a[None, :]).ravel()
    cols[m * 9:] = (3 * i[:, None] + b[None, :]).ravel()
    vals[m * 9:] = blocks.transpose(0, 2, 1).reshape(m, 9).ravel()

    hessian = sp.coo_matrix((vals, (rows, cols)), shape=(3 * n, 3 * n)).tocsr()

    # Diagonal super-elements: minus the sum of each site's off-diagonal 3x3
    # blocks. Accumulating per (a, b) component rather than by row sum is what
    # keeps the six rigid-body modes exactly at zero.
    diag_rows, diag_cols, diag_vals = [], [], []
    for aa in range(3):
        for bb in range(3):
            contrib = np.zeros(n)
            np.add.at(contrib, i, blocks[:, aa, bb])
            np.add.at(contrib, j, blocks[:, aa, bb])
            diag_rows.append(3 * np.arange(n) + aa)
            diag_cols.append(3 * np.arange(n) + bb)
            diag_vals.append(-contrib)
    diag = sp.coo_matrix(
        (np.concatenate(diag_vals),
         (np.concatenate(diag_rows), np.concatenate(diag_cols))),
        shape=(3 * n, 3 * n),
    ).tocsr()
    return (hessian + diag).tocsr()


# --------------------------------------------------------------------------
# Mode set
# --------------------------------------------------------------------------

@dataclass
class ModeSet:
    """Normal modes of an elastic network.

    ``vectors`` is ``(n_modes, n_sites, 3)`` so a mode can be added straight to
    a coordinate array. ``eigenvalues`` are in the arbitrary units set by
    ``gamma``; only their ratios are meaningful.
    """

    eigenvalues: np.ndarray
    vectors: np.ndarray
    n_sites: int
    symmetry: np.ndarray | None = None      # 'A' or 'E' per mode
    character: np.ndarray | None = None     # raw C3 character, for diagnostics
    meta: dict = field(default_factory=dict)

    @property
    def n_modes(self) -> int:
        return len(self.eigenvalues)

    @property
    def frequencies(self) -> np.ndarray:
        """Mode frequencies in arbitrary units (sqrt of the eigenvalue)."""
        return np.sqrt(np.maximum(self.eigenvalues, 0.0))

    def mode(self, index: int, amplitude: float = 1.0) -> np.ndarray:
        """Displacement field of one mode, scaled to a peak of ``amplitude`` A."""
        v = self.vectors[index]
        peak = np.linalg.norm(v, axis=1).max()
        return v * (amplitude / peak) if peak > 0 else v

    def msf(self, n_modes: int | None = None) -> np.ndarray:
        """Mean-square fluctuation per site, summed over modes as 1/lambda."""
        k = self.n_modes if n_modes is None else min(n_modes, self.n_modes)
        w = 1.0 / np.maximum(self.eigenvalues[:k], 1e-12)
        return np.einsum("m,mij,mij->i", w, self.vectors[:k], self.vectors[:k])

    def overlap(self, displacement: np.ndarray) -> np.ndarray:
        """Cosine overlap of each mode with an observed displacement field.

        The classic test of whether an elastic network "knows about" a real
        conformational change: overlap of the lowest few modes with the
        experimentally observed difference vector.
        """
        d = np.asarray(displacement, dtype=np.float64).reshape(-1)
        norm = np.linalg.norm(d)
        if norm == 0:
            return np.zeros(self.n_modes)
        d = d / norm
        flat = self.vectors.reshape(self.n_modes, -1)
        flat = flat / np.linalg.norm(flat, axis=1, keepdims=True)
        return np.abs(flat @ d)

    def cumulative_overlap(self, displacement: np.ndarray) -> np.ndarray:
        """Root-sum-square overlap accumulated over increasing mode count."""
        o = self.overlap(displacement)
        return np.sqrt(np.cumsum(o ** 2))

    def collectivity(self, index: int) -> float:
        """Bruschweiler's kappa: the fraction of sites a mode actually moves."""
        a = (self.vectors[index] ** 2).sum(axis=1)
        a = a / a.sum()
        nz = a[a > 0]
        return float(np.exp(-np.sum(nz * np.log(nz))) / self.n_sites)


# --------------------------------------------------------------------------
# The model
# --------------------------------------------------------------------------

@dataclass
class ANM:
    """Anisotropic network model over a set of C-alpha coordinates."""

    coords: np.ndarray
    cutoff: float = field(default_factory=lambda: _P.value("anm.cutoff"))
    gamma: float = field(default_factory=lambda: _P.value("anm.gamma"))
    spring: str = "inverse_square"
    d0: float = 7.5
    hessian: sp.csr_matrix | None = field(default=None, repr=False)

    def __post_init__(self) -> None:
        self.coords = np.ascontiguousarray(self.coords, dtype=np.float64)

    @property
    def n_sites(self) -> int:
        return len(self.coords)

    def build(self) -> "ANM":
        self.hessian = build_hessian(self.coords, cutoff=self.cutoff,
                                     gamma=self.gamma, spring=self.spring,
                                     d0=self.d0)
        return self

    def n_components(self) -> int:
        """Number of disconnected pieces in the contact network.

        This is not a curiosity. Each disconnected piece contributes its own
        six rigid-body modes, so a model containing a detached fragment — an
        isolated peptide, a stray ligand chain, a blade tip separated by a long
        unmodelled loop — yields 6 x N near-zero eigenvalues. Discarding only
        six then leaves rigid-body motions masquerading as the lowest
        functional modes, which is a silent and thoroughly misleading failure.
        """
        from scipy.sparse.csgraph import connected_components
        tree = cKDTree(self.coords)
        adj = tree.sparse_distance_matrix(tree, self.cutoff, output_type="coo_matrix")
        n, _ = connected_components(adj, directed=False)
        return int(n)

    def calc_modes(self, n_modes: int = 20, zero_modes: int | None = None,
                   shift: float = -1e-6, tol: float = 0.0) -> ModeSet:
        """Lowest non-trivial modes by shift-invert Lanczos.

        The rigid-body eigenvalues are discarded. ``zero_modes`` defaults to
        ``6 * n_components``, which is the correct count for a network that may
        not be fully connected; pass an explicit value to override.
        ``shift`` sits just below zero so shift-invert converges onto the low
        end of the spectrum.
        """
        if self.hessian is None:
            self.build()
        assert self.hessian is not None
        n_comp = self.n_components()
        if zero_modes is None:
            zero_modes = 6 * n_comp
        k = n_modes + zero_modes
        n = 3 * self.n_sites
        if k >= n - 1:
            raise ValueError(f"requested {k} modes from a {n}-dimensional problem")

        vals, vecs = spla.eigsh(self.hessian.astype(np.float64), k=k,
                                sigma=shift, which="LM", tol=tol)
        order = np.argsort(vals)
        vals, vecs = vals[order], vecs[:, order]

        # Discard the rigid-body modes, identified by near-zero eigenvalue.
        keep = slice(zero_modes, k)
        vals = vals[keep]
        vecs = vecs[:, keep]

        vectors = vecs.T.reshape(len(vals), self.n_sites, 3)
        return ModeSet(eigenvalues=vals, vectors=vectors, n_sites=self.n_sites,
                       meta={"cutoff": self.cutoff, "spring": self.spring,
                             "gamma": self.gamma, "d0": self.d0,
                             "n_components": n_comp,
                             "zero_modes_dropped": zero_modes,
                             "rigid_body_eigenvalues": np.sort(
                                 spla.eigsh(self.hessian, k=zero_modes,
                                            sigma=shift, which="LM",
                                            return_eigenvectors=False))
                             .tolist()})

    # ------------------------------------------------------------- symmetry

    def label_symmetry(self, modes: ModeSet, n_protomers: int = 3,
                       tolerance: float = 0.25) -> ModeSet:
        """Classify modes by their C3 character.

        The sites must be ordered protomer-by-protomer with identical ordering
        within each protomer — which is what
        :meth:`piezo1.physics.anm.ANM.from_trimer` guarantees.

        The character is ``<u, S u>`` where ``S`` cyclically permutes protomers
        and rotates the vectors by 120 degrees about the symmetry axis. It is
        ``+1`` for a totally symmetric ``A`` mode and ``-1/2`` for the ``E``
        pair.
        """
        n = self.n_sites
        per = n // n_protomers
        if per * n_protomers != n:
            raise ValueError("site count is not divisible by the protomer count")

        axis, origin = self._symmetry_axis(per, n_protomers)
        from ..structure.superpose import rotation_matrix
        rot = rotation_matrix(axis, 2.0 * np.pi / n_protomers)

        chars = np.empty(modes.n_modes)
        for m in range(modes.n_modes):
            u = modes.vectors[m]
            su = np.empty_like(u)
            for p in range(n_protomers):
                src = slice(p * per, (p + 1) * per)
                dst = slice(((p + 1) % n_protomers) * per,
                            ((p + 1) % n_protomers + 1) * per)
                su[dst] = u[src] @ rot.T
            chars[m] = float(np.sum(u * su) / np.sum(u * u))

        labels = np.where(np.abs(chars - 1.0) < tolerance, "A", "E")
        modes.symmetry = labels
        modes.character = chars
        return modes

    def _symmetry_axis(self, per: int, n_protomers: int) -> tuple[np.ndarray, np.ndarray]:
        from ..structure.superpose import detect_c3_axis
        chains = [self.coords[p * per:(p + 1) * per] for p in range(n_protomers)]
        ax = detect_c3_axis(chains)
        return ax.direction, ax.point

    # ------------------------------------------------------------ factories

    @classmethod
    def from_trimer(cls, chain_coords: list[np.ndarray], **kwargs) -> "ANM":
        """Build from per-protomer coordinate blocks, preserving their order.

        Keeping the protomer blocks contiguous and identically ordered is what
        makes :meth:`label_symmetry` valid, so trimeric models should always be
        constructed this way rather than by concatenating arbitrary selections.
        """
        counts = {len(c) for c in chain_coords}
        if len(counts) != 1:
            raise ValueError(f"protomers have differing site counts: {counts}")
        return cls(np.vstack(chain_coords), **kwargs)
