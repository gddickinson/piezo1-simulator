"""Pocket detection by Delaunay alpha spheres.

The fpocket construction (Le Guilloux, Schmidtke & Tuffery, *BMC Bioinformatics*
2009), reimplemented in numpy so nothing has to be installed or licensed.

The idea is geometric and rather elegant. Take the Delaunay tetrahedralisation
of the atom centres; every tetrahedron has a circumsphere that touches its four
atoms and, by the Delaunay property, contains no other atom. That empty sphere
is an **alpha sphere**, and its radius says what kind of space it occupies:

* **too small** (< ~3 Å) — the sphere sits inside the packed protein interior;
* **too large** (> ~6.5 Å) — it is out in bulk solvent;
* **in between** — it is in a cavity big enough for a small molecule.

Clustering the surviving spheres gives the pockets. No grid, no probe rolling,
no parameters beyond the two radii and a clustering distance.

**A caution about what a pocket means.** A geometric cavity is not a binding
site. This module finds voids; whether a ligand binds one is a separate
question that geometry alone cannot answer, and the ranking below is by size
and burial, not by affinity.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import Delaunay, cKDTree

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P

__all__ = ["AlphaSpheres", "Pocket", "find_pockets", "alpha_spheres",
           "ligand_contact_residues"]


@dataclass
class AlphaSpheres:
    """The filtered alpha spheres of a structure."""

    centers: np.ndarray           # (n, 3)
    radii: np.ndarray             # (n,)
    vertices: np.ndarray          # (n, 4) atom indices each sphere touches
    r_min: float
    r_max: float
    n_total: int = 0              # before filtering

    def __len__(self) -> int:
        return len(self.radii)


@dataclass
class Pocket:
    """A cluster of alpha spheres, i.e. a candidate binding cavity."""

    index: int
    centers: np.ndarray
    radii: np.ndarray
    residues: tuple[int, ...] = ()
    chains: tuple[str, ...] = ()
    atom_indices: np.ndarray = field(default_factory=lambda: np.zeros(0, int))
    buriedness: float = 0.0
    meta: dict = field(default_factory=dict)

    @property
    def n_spheres(self) -> int:
        return len(self.radii)

    @property
    def center(self) -> np.ndarray:
        return self.centers.mean(axis=0)

    @property
    def volume(self) -> float:
        """Union volume of the alpha spheres, by Monte-Carlo integration.

        The spheres overlap heavily, so summing 4/3πr³ would overcount badly —
        typically by several fold.
        """
        return float(self.meta.get("volume", 0.0))

    def contains_residues(self, residues) -> set[int]:
        return set(residues) & set(self.residues)

    def summary(self) -> str:
        return (f"pocket {self.index}: {self.n_spheres} spheres, "
                f"{self.volume:.0f} A^3, {len(self.residues)} lining residues, "
                f"buriedness {self.buriedness:.2f}")


# --------------------------------------------------------------------------
# Alpha spheres
# --------------------------------------------------------------------------

def _circumspheres(points: np.ndarray, simplices: np.ndarray
                   ) -> tuple[np.ndarray, np.ndarray]:
    """Circumcentre and circumradius of every tetrahedron, vectorised.

    For vertices ``p0..p3`` the centre satisfies ``|c − p0|² = |c − pi|²``,
    which linearises to ``2(pi − p0)·c = |pi|² − |p0|²`` — a 3×3 solve per
    tetrahedron, done in one batched ``np.linalg.solve``.
    """
    p = points[simplices]                       # (n, 4, 3)
    p0 = p[:, 0, :]
    a = 2.0 * (p[:, 1:, :] - p0[:, None, :])    # (n, 3, 3)
    sq = (p ** 2).sum(axis=2)                   # (n, 4)
    b = sq[:, 1:] - sq[:, [0]]                  # (n, 3)

    # Degenerate (near-flat) tetrahedra give singular systems; mark and skip.
    det = np.linalg.det(a)
    ok = np.abs(det) > 1e-8
    centers = np.full((len(simplices), 3), np.nan)
    if ok.any():
        # numpy 2 treats a 2-D right-hand side as a single matrix rather than a
        # batch of vectors, so the trailing axis has to be explicit.
        centers[ok] = np.linalg.solve(a[ok], b[ok][..., None])[..., 0]
    radii = np.linalg.norm(centers - p0, axis=1)
    return centers, radii


def alpha_spheres(coords: np.ndarray, r_min: float = 3.0,
                  r_max: float = 5.5, min_neighbours: int = 30,
                  neighbour_radius: float = 8.0) -> AlphaSpheres:
    """Delaunay alpha spheres of a point set, filtered by radius and burial.

    The radius filter alone is not enough on a large, open protein. PIEZO1 is a
    curved propeller with enormous solvent grooves between its blades, and with
    a radius filter only, single-linkage clustering percolates the whole
    surface into one object — the top "pocket" came out at 408 000 A^3 with 601
    lining residues, which is the protein's exterior, not a cavity.

    ``min_neighbours`` therefore requires each sphere to have that many atoms
    within ``neighbour_radius``, which discards spheres sitting on the open
    surface and leaves the enclosed ones. Together with a tighter ``r_max``
    this stops the percolation.
    """
    coords = np.ascontiguousarray(coords, dtype=np.float64)
    if len(coords) < 5:
        raise ValueError("need at least five atoms for a tetrahedralisation")
    tri = Delaunay(coords)
    centers, radii = _circumspheres(coords, tri.simplices)
    keep = np.isfinite(radii) & (radii >= r_min) & (radii <= r_max)
    if min_neighbours > 0 and keep.any():
        tree = cKDTree(coords)
        counts = np.array([len(x) for x in
                           tree.query_ball_point(centers[keep], neighbour_radius)])
        idx = np.flatnonzero(keep)
        keep[idx[counts < min_neighbours]] = False
    return AlphaSpheres(centers=centers[keep], radii=radii[keep],
                        vertices=tri.simplices[keep], r_min=r_min, r_max=r_max,
                        n_total=len(tri.simplices))


def _cluster(centers: np.ndarray, distance: float) -> np.ndarray:
    """Single-linkage clustering by connected components of a contact graph."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components
    n = len(centers)
    if n == 0:
        return np.zeros(0, dtype=int)
    pairs = np.asarray(list(cKDTree(centers).query_pairs(distance)), dtype=int)
    if len(pairs) == 0:
        return np.arange(n)
    graph = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])),
                       shape=(n, n))
    _, labels = connected_components(graph, directed=False)
    return labels


def _monte_carlo_volume(centers: np.ndarray, radii: np.ndarray,
                        n_samples: int = 6000, seed: int = 0) -> float:
    """Union volume of overlapping spheres."""
    if len(centers) == 0:
        return 0.0
    lo = (centers - radii[:, None]).min(axis=0)
    hi = (centers + radii[:, None]).max(axis=0)
    box = float(np.prod(hi - lo))
    if box <= 0:
        return 0.0
    rng = np.random.default_rng(seed)
    pts = rng.uniform(lo, hi, size=(n_samples, 3))
    inside = np.zeros(n_samples, dtype=bool)
    r2 = radii * radii
    # Chunked over spheres so the (n_samples x n_spheres) distance array never
    # has to exist all at once for a large pocket. Two refinements, both exact:
    # squared distances (the square root was never needed, since d <= r and
    # d^2 <= r^2 decide the same way), and only points still outside are
    # carried into the next block. Together ~1.9x, bit-identical.
    for start in range(0, len(centers), 400):
        outside = ~inside
        if not outside.any():
            break
        block = centers[start:start + 400]
        block_r2 = r2[start:start + 400]
        delta = pts[outside][:, None, :] - block[None, :, :]
        d2 = np.einsum("ijk,ijk->ij", delta, delta)
        inside[outside] = (d2 <= block_r2[None, :]).any(axis=1)
    return box * float(inside.mean())


def _buriedness(center: np.ndarray, coords: np.ndarray, tree: cKDTree,
                n_rays: int = 42, reach: float | None = None,
                clearance: float | None = None) -> float:
    """Fraction of directions from a point that run into protein.

    A cavity buried in the core scores near 1; a dimple on the surface scores
    low. Rays are sampled on a golden spiral so the result is deterministic.
    """
    clearance = _P.value("pockets.buriedness_clearance") if clearance is None else clearance
    reach = _P.value("pockets.buriedness_reach") if reach is None else reach
    i = np.arange(n_rays) + 0.5
    phi = np.arccos(1.0 - 2.0 * i / n_rays)
    theta = np.pi * (1.0 + 5.0 ** 0.5) * i
    dirs = np.stack([np.cos(theta) * np.sin(phi),
                     np.sin(theta) * np.sin(phi), np.cos(phi)], axis=1)
    steps = np.arange(2.0, reach, 1.5)
    blocked = 0
    for d in dirs:
        probe = center[None, :] + steps[:, None] * d[None, :]
        if (tree.query(probe)[0] < clearance).any():
            blocked += 1
    return blocked / n_rays


def find_pockets(structure: Structure, r_min: float = 3.0, r_max: float = 5.5,
                 cluster_distance: float | None = None, min_spheres: int = 20,
                 min_neighbours: int = 30,
                 protein_only: bool = True, max_pockets: int = 30,
                 lining_cutoff: float | None = None) -> list[Pocket]:
    """Detect and rank cavities.

    Ligands are excluded by default: a bound lipid fills the very cavity we are
    trying to find, so leaving it in makes the site disappear. That matters
    here — PIEZO1's structures contain resolved lipids sitting in exactly the
    pockets of interest.
    """
    lining_cutoff = _P.value("pockets.lining_cutoff") if lining_cutoff is None else lining_cutoff
    cluster_distance = _P.value("pockets.cluster_distance") if cluster_distance is None else cluster_distance
    mask = np.ones(structure.n_atoms, dtype=bool)
    if protein_only:
        mask &= structure.mask_protein() & ~structure.hetero
    coords = structure.xyz[mask].astype(np.float64)
    res_seq = structure.res_seq[mask]
    chains = structure.chain[mask]

    spheres = alpha_spheres(coords, r_min, r_max,
                            min_neighbours=min_neighbours)
    if len(spheres) == 0:
        return []

    labels = _cluster(spheres.centers, cluster_distance)
    tree = cKDTree(coords)
    pockets: list[Pocket] = []
    for label in np.unique(labels):
        sel = labels == label
        if int(sel.sum()) < min_spheres:
            continue
        centers, radii = spheres.centers[sel], spheres.radii[sel]
        neighbours = tree.query_ball_point(centers, lining_cutoff)
        atoms = np.unique(np.concatenate([np.asarray(n, dtype=int)
                                          for n in neighbours if len(n)])) \
            if any(len(n) for n in neighbours) else np.zeros(0, int)
        pockets.append(Pocket(
            index=0, centers=centers, radii=radii,
            residues=tuple(sorted(set(int(r) for r in res_seq[atoms]))),
            chains=tuple(sorted(set(str(c) for c in chains[atoms]))),
            atom_indices=atoms,
            buriedness=_buriedness(centers.mean(axis=0), coords, tree),
            meta={"volume": _monte_carlo_volume(centers, radii),
                  "n_spheres": int(sel.sum()),
                  "mean_radius": float(radii.mean())}))

    # Rank by volume weighted by burial: a large but open groove is less
    # interesting than a smaller enclosed cavity.
    pockets.sort(key=lambda p: -(p.volume * (0.5 + p.buriedness)))
    for i, p in enumerate(pockets[:max_pockets], start=1):
        p.index = i
    return pockets[:max_pockets]


# --------------------------------------------------------------------------
# Resolved ligands
# --------------------------------------------------------------------------

def ligand_contact_residues(structure: Structure, cutoff: float = 4.5
                            ) -> dict[str, dict]:
    """Residues contacting each resolved ligand or lipid, by chemical name."""
    lig = structure.mask_ligands()
    if not lig.any():
        return {}
    prot = structure.mask_protein() & ~structure.hetero
    tree = cKDTree(structure.xyz[prot].astype(np.float64))
    res_seq = structure.res_seq[prot]
    chains = structure.chain[prot]

    out: dict[str, dict] = {}
    lig_idx = np.flatnonzero(lig)
    for name in sorted(set(structure.res_name[lig].tolist())):
        sel = lig_idx[structure.res_name[lig_idx] == name]
        hits = tree.query_ball_point(structure.xyz[sel].astype(np.float64), cutoff)
        atoms = np.unique(np.concatenate([np.asarray(h, dtype=int)
                                          for h in hits if len(h)])) \
            if any(len(h) for h in hits) else np.zeros(0, int)
        out[str(name)] = {
            "n_copies": int(len(set(structure.res_seq[sel].tolist()))),
            "n_atoms": int(len(sel)),
            "residues": sorted(set(int(r) for r in res_seq[atoms])),
            "chains": sorted(set(str(c) for c in chains[atoms])),
            "cutoff": cutoff,
        }
    return out
