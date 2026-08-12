"""Pore-radius profiling along the ion conduction pathway.

Answers the question a viewer cannot: *is this structure open?* At each height
along the three-fold axis we find the largest sphere that fits in the pore
without overlapping any atom — the same quantity HOLE and CHAP report — and the
narrowest such sphere is the gate.

HOLE itself is not usable here: it has no Apple-Silicon build, and
MDAnalysis's ``hole2`` wrapper is an empty stub as of 2.10. This module is a
self-contained numpy implementation.

**The leash matters.** The clearance function has no interior maximum: a probe
free to wander leaves the pore sideways and grows without bound (unconstrained,
on real PIEZO1 coordinates, it escapes to R ≈ 6188 Å — it has found the bulk
solvent outside the protein, which is a true maximum and a useless answer). The
probe centre is therefore tethered within ``leash`` Å of the axis, which is
what makes the profile mean "radius of *the pore*".
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.spatial import cKDTree

from ..core.structure import Structure
from ..parameters import PARAMETERS as _P
from .superpose import SymmetryAxis

__all__ = ["PoreProfile", "pore_profile", "PoreSlice"]


@dataclass
class PoreSlice:
    z: float                  # height along the axis, Angstrom
    radius: float             # largest fitting sphere radius, Angstrom
    center: np.ndarray        # probe centre in world coordinates
    lining: tuple[int, ...] = ()      # residue numbers touching the probe
    lining_names: tuple[str, ...] = ()   # parallel to `lining`, same order


@dataclass
class PoreProfile:
    """Radius as a function of position along the conduction axis."""

    z: np.ndarray
    radius: np.ndarray
    centers: np.ndarray
    slices: list[PoreSlice] = field(default_factory=list)
    axis: SymmetryAxis | None = None
    meta: dict = field(default_factory=dict)

    # ------------------------------------------------------------- summary

    @property
    def bottleneck_index(self) -> int:
        return int(np.argmin(self.radius))

    @property
    def bottleneck_radius(self) -> float:
        return float(self.radius.min())

    @property
    def bottleneck_z(self) -> float:
        return float(self.z[self.bottleneck_index])

    def bottleneck_lining(self) -> tuple[int, ...]:
        return self.slices[self.bottleneck_index].lining if self.slices else ()

    def constrictions(self, threshold: float | None = None,
                      min_separation: float = 5.0) -> list[PoreSlice]:
        """Local minima narrower than ``threshold`` Å, i.e. candidate gates.

        ``threshold`` resolves from the registry at call time when not given,
        so an override takes effect on the next call rather than at import.

        3.0 Å is the conventional cut: a sphere of that radius is roughly a
        hydrated ion, so anything narrower is a barrier to permeation.
        """
        if threshold is None:
            threshold = _P.value("pore.constriction_threshold")
        out: list[PoreSlice] = []
        r = self.radius
        for i in range(1, len(r) - 1):
            if r[i] <= r[i - 1] and r[i] <= r[i + 1] and r[i] < threshold:
                if out and abs(self.z[i] - out[-1].z) < min_separation:
                    if r[i] < out[-1].radius:
                        out[-1] = self.slices[i]
                    continue
                out.append(self.slices[i])
        return out

    def is_conductive(self, ion_radius: float | None = None) -> bool:
        """Whether a dehydrated cation could pass the narrowest point.

        1.6 Å is roughly the Pauling radius of a hydrated-then-partly-stripped
        Na+/K+; this is a coarse indicator, not a permeation calculation.
        Resolved at call time, so an override applies to the next call.
        """
        if ion_radius is None:
            ion_radius = _P.value("pore.ion_radius")
        return self.bottleneck_radius >= ion_radius

    def summary(self) -> str:
        lining = ", ".join(str(r) for r in self.bottleneck_lining()[:6])
        return (f"bottleneck {self.bottleneck_radius:.2f} A at z = "
                f"{self.bottleneck_z:.1f} A, lined by {lining or 'n/a'}")


# --------------------------------------------------------------------------
# Core computation
# --------------------------------------------------------------------------

def _clearance(point: np.ndarray, tree: cKDTree, radii: np.ndarray,
               search: float) -> tuple[float, np.ndarray]:
    """Distance from ``point`` to the nearest atom *surface*.

    Returns ``(clearance, indices_of_touching_atoms)``. The Apollonius
    (additively weighted) distance is ``|p - c_i| - r_i``, so atoms of
    different size are handled correctly rather than by a single probe radius.
    """
    idx = tree.query_ball_point(point, search)
    if not idx:
        return search, np.zeros(0, dtype=int)
    idx = np.asarray(idx)
    d = np.linalg.norm(tree.data[idx] - point, axis=1) - radii[idx]
    best = float(d.min())
    touching = idx[d < best + 0.6]
    return best, touching


def _optimise_slice(z: float, axis_point: np.ndarray, e1: np.ndarray,
                    e2: np.ndarray, direction: np.ndarray, tree: cKDTree,
                    radii: np.ndarray, leash: float, search: float,
                    coarse: int = 7, refine_steps: int = 40
                    ) -> tuple[float, np.ndarray, np.ndarray]:
    """Largest fitting sphere in one plane, with the probe tethered to the axis.

    A coarse polar grid inside the leash finds the basin, then a shrinking
    local search refines it. Gradient methods are avoided because the clearance
    function is piecewise-smooth — its gradient jumps whenever the nearest atom
    changes — and a derivative-free refinement is both simpler and more robust.
    """
    origin = axis_point + direction * z

    # Coarse polar sampling, including the axis itself.
    offsets = [np.zeros(2)]
    for ring in range(1, coarse + 1):
        rr = leash * ring / coarse
        n_theta = max(6, int(2 * np.pi * rr / 1.5))
        for t in np.linspace(0, 2 * np.pi, n_theta, endpoint=False):
            offsets.append(np.array([rr * np.cos(t), rr * np.sin(t)]))
    offsets = np.array(offsets)

    best_val, best_off = -np.inf, np.zeros(2)
    for off in offsets:
        p = origin + off[0] * e1 + off[1] * e2
        val, _ = _clearance(p, tree, radii, search)
        if val > best_val:
            best_val, best_off = val, off

    # Local refinement by shrinking pattern search.
    step = leash / coarse
    for _ in range(refine_steps):
        improved = False
        for d in ((step, 0), (-step, 0), (0, step), (0, -step)):
            cand = best_off + np.array(d)
            if np.linalg.norm(cand) > leash:
                continue
            p = origin + cand[0] * e1 + cand[1] * e2
            val, _ = _clearance(p, tree, radii, search)
            if val > best_val:
                best_val, best_off, improved = val, cand, True
        if not improved:
            step *= 0.5
            if step < 0.02:
                break

    centre = origin + best_off[0] * e1 + best_off[1] * e2
    _, touching = _clearance(centre, tree, radii, search)
    return max(best_val, 0.0), centre, touching


def pore_profile(structure: Structure, axis: SymmetryAxis,
                 z_min: float | None = None, z_max: float | None = None,
                 step: float | None = None, leash: float | None = None,
                 search: float | None = None, protein_only: bool = True,
                 ) -> PoreProfile:
    """Compute the pore-radius profile of a channel about its symmetry axis.

    Parameters
    ----------
    structure:
        The model. Ligands and lipids are excluded by default: a lipid acyl
        chain lying in the vestibule would otherwise be reported as a
        constriction of the protein.
    axis:
        The conduction axis, normally from
        :func:`piezo1.structure.superpose.detect_c3_axis`.
    z_min, z_max:
        Range along the axis, relative to ``axis.point``. Defaults span the
        atoms lying within ``leash`` of the axis.
    step:
        Slice spacing in Angstrom.
    leash:
        Maximum distance of the probe centre from the axis. See the module
        docstring — without this the answer is meaningless.
    """
    # Resolved here rather than in the signature, so an override takes effect
    # on the next call instead of being frozen at import.
    if step is None:
        step = _P.value("pore.step")
    if leash is None:
        leash = _P.value("pore.leash")
    if search is None:
        search = _P.value("pore.search")

    mask = structure.mask_protein() if protein_only else np.ones(
        structure.n_atoms, dtype=bool)
    if protein_only:
        mask &= ~structure.hetero
    xyz = structure.xyz[mask].astype(np.float64)
    radii = structure.vdw_radii()[mask].astype(np.float64)
    res_seq = structure.res_seq[mask]
    res_name = structure.res_name[mask]
    if len(xyz) < 100:
        raise ValueError("too few atoms to profile a pore")

    tree = cKDTree(xyz)
    direction = axis.direction / np.linalg.norm(axis.direction)

    # An orthonormal frame perpendicular to the axis.
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(tmp, direction)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    e1 = np.cross(direction, tmp)
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(direction, e1)

    proj = axis.project(xyz)
    radial = axis.radial(xyz)
    near = radial < leash + 12.0
    if z_min is None:
        z_min = float(np.percentile(proj[near], 1)) if near.any() else float(proj.min())
    if z_max is None:
        z_max = float(np.percentile(proj[near], 99)) if near.any() else float(proj.max())

    zs = np.arange(z_min, z_max + step * 0.5, step)
    radius = np.zeros(len(zs))
    centers = np.zeros((len(zs), 3))
    slices: list[PoreSlice] = []

    for i, z in enumerate(zs):
        r, c, touching = _optimise_slice(z, axis.point, e1, e2, direction,
                                         tree, radii, leash, search)
        radius[i] = r
        centers[i] = c
        # Numbers and names must come out of the SAME de-duplication, or they
        # are not parallel and every consumer that zips them mislabels the
        # pore. They were two independently sorted sets until Round 84c: on
        # 8YEZ that reported the worst dewetted residue as GLU2510 when it is
        # PRO2510, and on 11YE as LEU2427 when it is ILE2296 — a different
        # residue, not just a different name. `zip` also truncates to the
        # shorter tuple, so a slice touching three residues that share a name
        # contributed one lining point instead of three.
        if len(touching):
            unique, first = np.unique(res_seq[touching], return_index=True)
            lining_res = tuple(int(v) for v in unique)
            lining_nm = tuple(str(v) for v in res_name[touching][first])
        else:
            lining_res, lining_nm = (), ()
        slices.append(PoreSlice(z=float(z), radius=float(r), center=c,
                                lining=lining_res, lining_names=lining_nm))

    return PoreProfile(z=zs, radius=radius, centers=centers, slices=slices,
                       axis=axis,
                       meta={"leash": leash, "step": step,
                             "n_atoms": int(len(xyz)),
                             "protein_only": protein_only,
                             "note": ("Radius is the largest sphere fitting "
                                      "without overlapping any van der Waals "
                                      "surface, with the probe tethered within "
                                      f"{leash} A of the axis.")})
