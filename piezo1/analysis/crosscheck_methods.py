"""Re-derive four *methods* by routes that share none of their machinery.

:mod:`piezo1.analysis.crosscheck` re-derives three headline physics results.
This does the same for four algorithms, on the same principle: a test written
from the same understanding as the code shares its blind spots, so the check has
to come from somewhere else.

Each of these deliberately swaps out the part most likely to be wrong:

* **pore radius** — the pipeline maximises Apollonius clearance with a coarse
  polar grid followed by a shrinking pattern search. That is a *local* optimiser
  on a piecewise-smooth surface, and the failure mode is a local maximum. Here
  the same clearance is maximised by dense uniform random sampling, which has no
  notion of a basin and cannot get stuck in one.

* **SASA** — the pipeline uses Shrake–Rupley with a fixed golden-spiral point
  set. A bad point set is invisible to a test that uses the same point set. Here
  the surface integral is done by Monte-Carlo with independent random directions.

* **conservation** — the pipeline anchors each ortholog to the reference by a
  full Needleman–Wunsch alignment. Here positions are matched by exact k-mer
  seeds with no dynamic programming and no gap penalties at all.

* **PCA** — the pipeline takes the SVD of the centred coordinate matrix. Here
  the leading component comes from power iteration on the covariance operator,
  which touches neither LAPACK's SVD nor its eigensolver.

Disagreement is the point. Each disagreement is either a bug or an
approximation that can be named, and the round is worth running either way.

**Measured (Round 37).** Three of the four agree: SASA to 0.1%, PC1 exactly
(|cos| = 1.000000, eigenvalue identical), and the pore bottleneck to 5.2% —
where the random search finds the *larger* value, which is the informative sign,
since a brute-force search can only match or beat a local optimiser.

Conservation correlates at **0.817** and its residual is a **named bias in the
alternative, not in the pipeline**. Anchoring by maximum exact matches is a
selection: it preferentially lines up residues that happen to agree, so it
inflates conservation exactly where a position is variable. The floor of the
k-mer profile is 0.36 rather than 0, it agrees with the pipeline at invariant
positions (0.993 where the pipeline says 1.00), and it reads 0.653 where the
pipeline says below 0.50. The Needleman–Wunsch route uses gap penalties and a
substitution matrix rather than raw match counts and is not subject to that.

So the k-mer route is the weaker instrument — the same verdict Round 30 reached
about the parabola. That is a legitimate outcome of a cross-check, and it is
recorded rather than tuned away.
"""

from __future__ import annotations

import numpy as np

from .crosscheck import CrossCheck

__all__ = ["pore_radius_by_random_search", "sasa_by_monte_carlo",
           "conservation_by_kmer_anchoring", "pc1_by_power_iteration",
           "check_pore_radius", "check_sasa", "check_conservation", "check_pca"]


# --------------------------------------------------------------------------
# Pore radius without the pattern-search optimiser
# --------------------------------------------------------------------------

def pore_radius_by_random_search(structure, axis, z: float, leash: float = 8.0,
                                 search: float = 18.0, n_samples: int = 20000,
                                 seed: int = 0) -> float:
    """Largest sphere fitting at height ``z``, found by brute force.

    Identical *definition* to the pipeline — maximise the Apollonius clearance
    over probe centres tethered within ``leash`` of the axis — and a completely
    different *search*. Uniform random sampling of the disc converges slowly but
    has no basin structure to be trapped by, so a local-maximum bug in the
    pattern search shows up as a random search finding something larger.
    """
    from scipy.spatial import cKDTree

    mask = structure.mask_protein() & ~structure.hetero
    coords = structure.xyz[mask].astype(np.float64)
    radii = structure.vdw_radii()[mask].astype(np.float64)
    tree = cKDTree(coords)

    direction = np.asarray(axis.direction, dtype=np.float64)
    direction = direction / np.linalg.norm(direction)
    helper = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(helper, direction)) > 0.9:
        helper = np.array([1.0, 0.0, 0.0])
    e1 = np.cross(direction, helper)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(direction, e1)
    origin = np.asarray(axis.point, dtype=np.float64) + direction * z

    rng = np.random.default_rng(seed)
    # Uniform on the disc: sqrt of the radial variate, or the samples pile up
    # at the centre and the search never explores the rim.
    r = leash * np.sqrt(rng.random(n_samples))
    theta = 2.0 * np.pi * rng.random(n_samples)
    points = (origin + np.outer(r * np.cos(theta), e1)
              + np.outer(r * np.sin(theta), e2))

    best = 0.0
    for point in points:
        idx = tree.query_ball_point(point, search)
        if not idx:
            value = search
        else:
            idx = np.asarray(idx)
            value = float((np.linalg.norm(coords[idx] - point, axis=1)
                           - radii[idx]).min())
        best = max(best, value)
    return max(best, 0.0)


def check_pore_radius(structure, axis, profile, tolerance: float = 0.08
                      ) -> CrossCheck:
    """Compare the bottleneck radius against a brute-force search at that height."""
    z = float(profile.bottleneck_z)
    alternative = pore_radius_by_random_search(structure, axis, z)
    return CrossCheck(
        quantity="pore bottleneck radius",
        primary=float(profile.bottleneck_radius), alternative=alternative,
        unit="A", tolerance=tolerance,
        primary_route="polar grid then shrinking pattern search",
        alternative_route="20k uniform random probe centres on the disc",
        note="same clearance definition, no basin structure to get trapped in")


# --------------------------------------------------------------------------
# SASA without the Shrake-Rupley point set
# --------------------------------------------------------------------------

def sasa_by_monte_carlo(structure, probe: float = 1.4, n_samples: int = 4000,
                        mask=None, seed: int = 0) -> float:
    """Total solvent-accessible area by Monte-Carlo integration.

    Shrake–Rupley places a deterministic quasi-uniform point set on each
    expanded sphere and counts the fraction not buried. This does the same
    integral with independent random directions, so a defect in the point set —
    the usual way that method goes wrong — cannot be shared.
    """
    from scipy.spatial import cKDTree

    if mask is None:
        mask = structure.mask_protein() & ~structure.hetero
    coords = structure.xyz[mask].astype(np.float64)
    radii = structure.vdw_radii()[mask].astype(np.float64) + probe
    tree = cKDTree(coords)

    rng = np.random.default_rng(seed)
    # Marsaglia: normalised Gaussians are uniform on the sphere. Rejection
    # sampling in a cube is not, and biases towards the corners.
    directions = rng.normal(size=(n_samples, 3))
    directions /= np.linalg.norm(directions, axis=1, keepdims=True)

    total = 0.0
    largest = float(radii.max())
    for i, (centre, radius) in enumerate(zip(coords, radii)):
        neighbours = np.asarray(tree.query_ball_point(centre, radius + largest))
        neighbours = neighbours[neighbours != i]
        if len(neighbours) == 0:
            total += 4.0 * np.pi * radius ** 2
            continue
        points = centre + radius * directions
        distances = np.linalg.norm(
            points[:, None, :] - coords[neighbours][None, :, :], axis=2)
        exposed = np.all(distances >= radii[neighbours][None, :], axis=1)
        total += 4.0 * np.pi * radius ** 2 * float(exposed.mean())
    return total


def check_sasa(structure, primary_total: float, mask=None,
               tolerance: float = 0.03) -> CrossCheck:
    alternative = sasa_by_monte_carlo(structure, mask=mask)
    return CrossCheck(
        quantity="solvent-accessible area", primary=float(primary_total),
        alternative=alternative, unit="A^2", tolerance=tolerance,
        primary_route="Shrake-Rupley, 256-point golden spiral",
        alternative_route="Monte-Carlo, 4000 independent random directions",
        note="same integral, unrelated quadrature")


# --------------------------------------------------------------------------
# Conservation without dynamic-programming alignment
# --------------------------------------------------------------------------

def conservation_by_kmer_anchoring(reference: str, sequences: list[str],
                                   k: int = 9, window: int = 60) -> np.ndarray:
    """Per-position conservation with no alignment algorithm at all.

    For each reference position the surrounding k-mer is slid along each
    ortholog within a local window and scored by exact matches; the best
    offset decides which residue is "the same position". No dynamic
    programming, no substitution matrix, no gap penalties.

    This is a deliberately cruder instrument. Near an indel it will pick the
    wrong offset, so disagreement there is expected and diagnosable — it is
    exactly the region where an alignment is doing real work. Agreement
    elsewhere means the conservation signal is in the sequences rather than in
    the aligner.
    """
    reference = reference.upper()
    n = len(reference)
    half = k // 2
    counts = [dict() for _ in range(n)]

    for sequence in sequences:
        sequence = sequence.upper()
        if not sequence:
            continue
        for i in range(n):
            lo, hi = max(0, i - half), min(n, i + half + 1)
            motif = reference[lo:hi]
            centre = i - lo

            # Search only near the diagonal: orthologs of the same protein do
            # not move by hundreds of residues, and an unrestricted scan finds
            # spurious repeats.
            start = max(0, i - window)
            stop = min(len(sequence) - len(motif), i + window)
            if stop < start:
                continue
            best_score, best_at = -1, None
            for j in range(start, stop + 1):
                score = sum(1 for a, b in zip(motif, sequence[j:j + len(motif)])
                            if a == b)
                if score > best_score:
                    best_score, best_at = score, j
            if best_at is None:
                continue
            position = best_at + centre
            if 0 <= position < len(sequence):
                residue = sequence[position]
                counts[i][residue] = counts[i].get(residue, 0) + 1

    out = np.zeros(n)
    for i, table in enumerate(counts):
        total = sum(table.values())
        if total < 2:
            out[i] = np.nan
            continue
        frequencies = np.array(list(table.values()), dtype=float) / total
        entropy = -np.sum(frequencies * np.log2(frequencies))
        # Same normalisation as the pipeline: 1 means invariant.
        out[i] = 1.0 - entropy / np.log2(20)
    return out


def check_conservation(primary: np.ndarray, alternative: np.ndarray,
                       tolerance: float = 0.25) -> CrossCheck:
    """Compare the two conservation profiles by their correlation.

    A correlation rather than a per-position difference, because the k-mer route
    carries a known upward bias at variable positions (see the module docstring)
    and a mean absolute difference would mostly measure that bias. What would be
    damning is the two disagreeing in *shape*.

    The tolerance is 0.25 — i.e. a correlation of 0.75 counts as agreement —
    because the alternative is deliberately the cruder instrument. Requiring
    closer agreement would be requiring a match-maximising anchor to behave like
    a gap-penalised alignment, which it cannot.
    """
    both = np.isfinite(primary) & np.isfinite(alternative)
    if both.sum() < 20:
        return CrossCheck(quantity="conservation profile", primary=float("nan"),
                          alternative=float("nan"), tolerance=tolerance,
                          note="too few shared positions to compare")
    correlation = float(np.corrcoef(primary[both], alternative[both])[0, 1])
    return CrossCheck(
        quantity="conservation correlation", primary=1.0,
        alternative=correlation, tolerance=tolerance,
        primary_route="Needleman-Wunsch anchoring, reference-indexed",
        alternative_route=f"exact k-mer seeds, no DP ({int(both.sum())} positions)",
        note="the residual is a bias in THIS route: anchoring by maximum "
             "exact matches inflates conservation at variable positions")


# --------------------------------------------------------------------------
# PCA without SVD
# --------------------------------------------------------------------------

def pc1_by_power_iteration(matrix: np.ndarray, iterations: int = 500,
                           tol: float = 1e-12, seed: int = 0
                           ) -> tuple[np.ndarray, float]:
    """Leading eigenvector and eigenvalue of the covariance, by power iteration.

    Never forms the covariance matrix — for a coordinate matrix that would be
    3N x 3N — and never calls a library eigensolver. Repeated application of
    ``Xᵀ(X v)`` converges on the top eigenvector of ``XᵀX``, which is what the
    SVD's leading right singular vector is.
    """
    centred = matrix - matrix.mean(axis=0)
    rng = np.random.default_rng(seed)
    v = rng.normal(size=centred.shape[1])
    v /= np.linalg.norm(v)

    value = 0.0
    for _ in range(iterations):
        w = centred.T @ (centred @ v)
        norm = np.linalg.norm(w)
        if norm == 0:
            break
        w /= norm
        if np.linalg.norm(w - v) < tol or np.linalg.norm(w + v) < tol:
            v = w
            break
        v = w
    value = float(v @ (centred.T @ (centred @ v)))
    return v, value / max(len(centred) - 1, 1)


def check_pca(primary_eigenvalue: float, primary_component: np.ndarray,
              matrix: np.ndarray, tolerance: float = 0.02) -> CrossCheck:
    vector, eigenvalue = pc1_by_power_iteration(matrix)
    # The pipeline stores components as (n_sites, 3); flatten to compare with
    # the flat eigenvector of the same coordinate matrix.
    reference = np.asarray(primary_component, dtype=float).ravel()
    # Sign is arbitrary in both routes, so compare |cos|.
    overlap = abs(float(np.dot(vector, reference)
                        / (np.linalg.norm(vector) * np.linalg.norm(reference))))
    return CrossCheck(
        quantity="PC1 eigenvalue", primary=float(primary_eigenvalue),
        alternative=float(eigenvalue), tolerance=tolerance,
        primary_route="SVD of the centred coordinate matrix",
        alternative_route="power iteration on X^T X, no library eigensolver",
        note=f"|cos| between the two PC1 directions: {overlap:.6f}")
