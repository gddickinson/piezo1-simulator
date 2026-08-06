"""Allostery: how force applied at one place moves another.

PIEZO1's blades are up to 100 Å from the gate they open. The question this
module answers is mechanical rather than descriptive — push here, what moves
there, and along which route does the signal travel?

Three standard analyses, all derived from the elastic network's covariance
matrix ``C = H⁺ = Σ_k (1/λ_k) v_k v_kᵀ``:

**Perturbation Response Scanning** (Atilgan & Atilgan, *PLoS Comput Biol* 2009).
Apply a unit force at residue *i* and measure the displacement everywhere else.
Averaged over random force directions, the mean-square response at *j* is
``‖C_ij‖²_F``, the squared Frobenius norm of the 3×3 covariance block. Averaging
that matrix along its rows gives each residue's **effectiveness** — how well
pushing there moves the rest of the protein — and along its columns its
**sensitivity** — how much it moves when pushed elsewhere. Effectors and
sensors respectively.

**Dynamic cross-correlation.** ``DCC_ij = tr(C_ij)/√(tr(C_ii)·tr(C_jj))``, in
[−1, 1]. Positive means the two residues move together, negative means they
move oppositely — which for a lever is exactly what you expect either side of
the fulcrum.

**Allosteric pathway.** Following Sethi et al. (*PNAS* 2009), build a graph over
residues in spatial contact, weight each edge ``−log|DCC_ij|``, and take the
shortest path. Strongly correlated neighbours are "cheap" to traverse, so the
shortest path is the route along which motion is most reliably transmitted.

The covariance is never formed in full: for a PIEZO1 trimer it would be an
11466² matrix, about a gigabyte. Everything here works from the modes directly,
in chunks.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

__all__ = ["PRSResult", "AllostericPath", "covariance_blocks_norm",
           "perturbation_response", "cross_correlation", "allosteric_path",
           "effectiveness_profile", "build_network", "path_betweenness",
           "detour_cost"]


# --------------------------------------------------------------------------
# Covariance-derived matrices
# --------------------------------------------------------------------------

def _weighted_vectors(modes, n_modes: int | None = None) -> np.ndarray:
    """``(n_sites, 3, n_modes)`` of modes scaled by ``1/√λ``.

    The covariance is ``Σ_k v_k v_kᵀ / λ_k``, so folding ``1/√λ`` into the
    vectors makes every downstream product a plain contraction. Softer modes
    dominate, which is the physics: a mode with a small eigenvalue is one the
    protein actually moves along.
    """
    k = modes.n_modes if n_modes is None else min(n_modes, modes.n_modes)
    w = 1.0 / np.sqrt(np.maximum(modes.eigenvalues[:k], 1e-12))
    return np.transpose(modes.vectors[:k], (1, 2, 0)) * w[None, None, :]


def covariance_blocks_norm(modes, n_modes: int | None = None,
                           chunk: int = 128) -> np.ndarray:
    """``‖C_ij‖²_F`` for every residue pair — the raw PRS response matrix.

    Computed in row chunks so peak memory stays at ``chunk × n_sites × 9``
    rather than the ``n_sites² × 9`` a direct einsum would need.
    """
    a = _weighted_vectors(modes, n_modes)          # (n, 3, k)
    n = len(a)
    out = np.empty((n, n), dtype=np.float64)
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        # block[i, j, alpha, beta] = sum_k a[i, alpha, k] a[j, beta, k]
        block = np.einsum("iak,jbk->ijab", a[start:stop], a, optimize=True)
        out[start:stop] = (block ** 2).sum(axis=(2, 3))
    return out


def cross_correlation(modes, n_modes: int | None = None,
                      chunk: int = 256) -> np.ndarray:
    """Normalised dynamic cross-correlation, in [−1, 1]."""
    a = _weighted_vectors(modes, n_modes)
    n = len(a)
    trace = np.empty((n, n), dtype=np.float64)
    for start in range(0, n, chunk):
        stop = min(start + chunk, n)
        # tr(C_ij) = sum_alpha sum_k a[i, alpha, k] a[j, alpha, k]
        trace[start:stop] = np.einsum("iak,jak->ij", a[start:stop], a,
                                      optimize=True)
    diag = np.sqrt(np.maximum(np.diag(trace), 1e-30))
    return np.clip(trace / np.outer(diag, diag), -1.0, 1.0)


# --------------------------------------------------------------------------
# Perturbation response scanning
# --------------------------------------------------------------------------

@dataclass
class PRSResult:
    """Response matrix and the effector/sensor profiles derived from it."""

    matrix: np.ndarray            # (n, n), response at j to a force at i
    residues: np.ndarray          # residue numbers, one per site
    chains: np.ndarray | None = None
    meta: dict = field(default_factory=dict)

    @property
    def effectiveness(self) -> np.ndarray:
        """Row average: how well pushing here moves everything else."""
        m = self.matrix.copy()
        np.fill_diagonal(m, 0.0)
        return m.mean(axis=1)

    @property
    def sensitivity(self) -> np.ndarray:
        """Column average: how much this moves when pushed from elsewhere."""
        m = self.matrix.copy()
        np.fill_diagonal(m, 0.0)
        return m.mean(axis=0)

    def response_at(self, target_sites) -> np.ndarray:
        """Response at a chosen set of sites to a force at each site in turn.

        This is the question the project actually cares about: which residues,
        when pushed, move the *gate*.
        """
        idx = np.asarray(list(target_sites), dtype=int)
        m = self.matrix.copy()
        np.fill_diagonal(m, 0.0)
        return m[:, idx].mean(axis=1)

    def top_effectors(self, n: int = 10) -> list[tuple[int, float]]:
        order = np.argsort(self.effectiveness)[::-1][:n]
        return [(int(self.residues[i]), float(self.effectiveness[i])) for i in order]

    def top_sensors(self, n: int = 10) -> list[tuple[int, float]]:
        order = np.argsort(self.sensitivity)[::-1][:n]
        return [(int(self.residues[i]), float(self.sensitivity[i])) for i in order]

    def per_residue(self, values: np.ndarray) -> dict[int, float]:
        """Collapse per-site values onto residue numbers, averaging protomers."""
        out: dict[int, list[float]] = {}
        for res, v in zip(self.residues, values):
            out.setdefault(int(res), []).append(float(v))
        return {k: float(np.mean(v)) for k, v in out.items()}


def perturbation_response(modes, residues: np.ndarray,
                          chains: np.ndarray | None = None,
                          n_modes: int | None = None,
                          normalise: bool = True) -> PRSResult:
    """Run PRS over an elastic network.

    ``normalise`` divides each row by its own mean, which removes the trivial
    effect that surface residues respond more than buried ones simply because
    they are less restrained, and leaves the *pattern* of propagation.
    """
    matrix = covariance_blocks_norm(modes, n_modes)
    if normalise:
        row = matrix.mean(axis=1, keepdims=True)
        matrix = matrix / np.maximum(row, 1e-30)
    return PRSResult(matrix=matrix, residues=np.asarray(residues),
                     chains=chains,
                     meta={"n_modes": n_modes or modes.n_modes,
                           "normalised": normalise,
                           "n_sites": len(matrix)})


def effectiveness_profile(prs: PRSResult) -> dict[int, float]:
    return prs.per_residue(prs.effectiveness)


# --------------------------------------------------------------------------
# Allosteric pathways
# --------------------------------------------------------------------------

@dataclass
class AllostericPath:
    """A route through the residue network from a source to a target."""

    sites: list[int]                  # site indices along the path
    residues: list[int]               # residue numbers along the path
    cost: float                       # summed -log|DCC| edge weights
    correlations: list[float] = field(default_factory=list)
    meta: dict = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.sites)

    def passes_through(self, lo: int, hi: int) -> bool:
        """Whether the path visits any residue in ``[lo, hi]``."""
        return any(lo <= r <= hi for r in self.residues)

    def residues_in(self, lo: int, hi: int) -> list[int]:
        return [r for r in self.residues if lo <= r <= hi]

    def summary(self) -> str:
        return (f"{len(self)} steps, cost {self.cost:.2f}, "
                f"residues {self.residues[0]} -> {self.residues[-1]}")


def allosteric_path(coords: np.ndarray, dcc: np.ndarray,
                    source_sites, target_sites,
                    residues: np.ndarray,
                    contact_cutoff: float = 10.0,
                    min_correlation: float = 1e-3) -> AllostericPath:
    """Shortest correlation-weighted path from a source to a target set.

    Edges join residues within ``contact_cutoff`` and cost ``−log|DCC_ij|``, so
    a pair that moves together is cheap to cross and an uncorrelated pair is
    expensive. The shortest path is therefore the route along which motion is
    transmitted most reliably, not merely the shortest way through space.
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import dijkstra

    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    if dcc.shape != (n, n):
        raise ValueError(f"correlation matrix is {dcc.shape}, expected ({n}, {n})")

    pairs = np.asarray(list(cKDTree(coords).query_pairs(contact_cutoff)),
                       dtype=int)
    if len(pairs) == 0:
        raise ValueError(f"no residue pairs within {contact_cutoff} A")

    corr = np.abs(dcc[pairs[:, 0], pairs[:, 1]])
    weight = -np.log(np.maximum(corr, min_correlation))
    graph = coo_matrix(
        (np.concatenate([weight, weight]),
         (np.concatenate([pairs[:, 0], pairs[:, 1]]),
          np.concatenate([pairs[:, 1], pairs[:, 0]]))),
        shape=(n, n)).tocsr()

    sources = np.asarray(list(source_sites), dtype=int)
    targets = np.asarray(list(target_sites), dtype=int)
    dist, predecessors = dijkstra(graph, indices=sources,
                                  return_predecessors=True, directed=False)

    sub = dist[:, targets]
    if not np.isfinite(sub).any():
        raise ValueError("no connected path between the source and target sets")
    si, ti = np.unravel_index(int(np.nanargmin(np.where(np.isfinite(sub), sub, np.inf))),
                              sub.shape)
    start, end = int(sources[si]), int(targets[ti])

    walk = [end]
    while walk[-1] != start:
        prev = predecessors[si, walk[-1]]
        if prev < 0:
            raise ValueError("path reconstruction failed")
        walk.append(int(prev))
    walk.reverse()

    corrs = [float(abs(dcc[a, b])) for a, b in zip(walk, walk[1:])]
    return AllostericPath(
        sites=walk, residues=[int(residues[i]) for i in walk],
        cost=float(sub[si, ti]), correlations=corrs,
        meta={"contact_cutoff": contact_cutoff,
              "n_edges": int(len(pairs)),
              "source_site": start, "target_site": end,
              "weight": "-log|DCC|"})


def build_network(coords: np.ndarray, dcc: np.ndarray,
                  contact_cutoff: float = 10.0,
                  min_correlation: float = 1e-3):
    """Correlation-weighted contact graph, as a sparse matrix."""
    from scipy.sparse import coo_matrix
    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    pairs = np.asarray(list(cKDTree(coords).query_pairs(contact_cutoff)), dtype=int)
    if len(pairs) == 0:
        raise ValueError(f"no residue pairs within {contact_cutoff} A")
    w = -np.log(np.maximum(np.abs(dcc[pairs[:, 0], pairs[:, 1]]), min_correlation))
    return coo_matrix(
        (np.concatenate([w, w]),
         (np.concatenate([pairs[:, 0], pairs[:, 1]]),
          np.concatenate([pairs[:, 1], pairs[:, 0]]))), shape=(n, n)).tocsr()


def detour_cost(graph, source_sites, target_sites, via_sites) -> dict:
    """Cheapest source→target path forced through ``via_sites``, versus free.

    The obvious way to ask "does the signal go through X" is to compute
    source→X and X→target separately and add them. That is wrong, and wrong in
    a way that produces an impossible answer: each leg independently picks its
    best endpoints, which for a C3 trimer may be in *different protomers*, so
    the two legs do not join up. Done that way the "detour" came out cheaper
    than the unconstrained shortest path, which cannot happen.

    The correct comparison shares a single via-point: minimise
    ``d(source → v) + d(v → target)`` over ``v`` in ``via_sites``.
    """
    from scipy.sparse.csgraph import dijkstra
    src = np.asarray(list(source_sites), dtype=int)
    tgt = np.asarray(list(target_sites), dtype=int)
    via = np.asarray(list(via_sites), dtype=int)

    d_from_src = dijkstra(graph, indices=src, directed=False).min(axis=0)
    d_to_tgt = dijkstra(graph, indices=tgt, directed=False).min(axis=0)

    direct = float(np.nanmin(d_from_src[tgt]))
    through = d_from_src[via] + d_to_tgt[via]
    if not np.isfinite(through).any():
        return {"direct": direct, "via": float("inf"), "penalty": float("inf"),
                "best_via_site": None}
    best = int(np.nanargmin(np.where(np.isfinite(through), through, np.inf)))
    return {"direct": direct, "via": float(through[best]),
            "penalty": float(through[best] - direct),
            "best_via_site": int(via[best])}


def path_betweenness(graph, source_sites, target_sites, residues: np.ndarray,
                     max_pairs: int = 400, seed: int = 0) -> dict[int, float]:
    """How often each residue lies on a shortest source→target path.

    A single shortest path is a fragile thing to draw conclusions from — one
    slightly-better edge reroutes it entirely. Counting residue occurrences
    over many source/target pairs is the standard dynamical-network-analysis
    answer and is far more stable.
    """
    from scipy.sparse.csgraph import dijkstra
    rng = np.random.default_rng(seed)
    src = np.asarray(list(source_sites), dtype=int)
    tgt = np.asarray(list(target_sites), dtype=int)
    if len(src) * len(tgt) > max_pairs:
        n_src = max(1, int(np.sqrt(max_pairs * len(src) / max(len(tgt), 1))))
        src = rng.choice(src, size=min(n_src, len(src)), replace=False)
        n_tgt = max(1, max_pairs // max(len(src), 1))
        tgt = rng.choice(tgt, size=min(n_tgt, len(tgt)), replace=False)

    dist, pred = dijkstra(graph, indices=src, return_predecessors=True,
                          directed=False)
    counts: dict[int, int] = {}
    n_paths = 0
    for i in range(len(src)):
        for t in tgt:
            if not np.isfinite(dist[i, t]):
                continue
            node, guard = int(t), 0
            while node != src[i] and guard < 10000:
                counts[int(residues[node])] = counts.get(int(residues[node]), 0) + 1
                node = int(pred[i, node])
                if node < 0:
                    break
                guard += 1
            n_paths += 1
    if n_paths == 0:
        return {}
    return {k: v / n_paths for k, v in counts.items()}
