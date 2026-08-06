"""Rigid-body superposition and C3 symmetry handling.

PIEZO1 is a homotrimer, and almost every analysis in this project needs either
(a) two structures put in a common frame, or (b) the molecular three-fold axis
that defines "up" through the membrane.  Both live here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "kabsch", "superpose", "rmsd", "SymmetryAxis", "detect_c3_axis",
    "rotation_matrix", "align_axis_to_z", "match_protomers", "ProtomerMatch",
]


def kabsch(mobile: np.ndarray, target: np.ndarray,
           weights: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Optimal rigid transform mapping ``mobile`` onto ``target``.

    Returns ``(rotation, translation, mobile_centroid)`` such that
    ``(mobile - mobile_centroid) @ rotation.T + translation`` is the fitted
    coordinate set.  Reflections are excluded.
    """
    mobile = np.asarray(mobile, dtype=np.float64)
    target = np.asarray(target, dtype=np.float64)
    if mobile.shape != target.shape:
        raise ValueError(f"shape mismatch: {mobile.shape} vs {target.shape}")
    if weights is None:
        w = np.ones(len(mobile))
    else:
        w = np.asarray(weights, dtype=np.float64)
    w = w / w.sum()

    cm = (mobile * w[:, None]).sum(axis=0)
    ct = (target * w[:, None]).sum(axis=0)
    p = mobile - cm
    q = target - ct

    cov = (p * w[:, None]).T @ q
    v, _s, wt = np.linalg.svd(cov)
    d = np.sign(np.linalg.det(v @ wt))
    correction = np.diag([1.0, 1.0, d])
    rot = (v @ correction @ wt).T          # applies as x @ rot.T
    return rot, ct, cm


def superpose(mobile: np.ndarray, target: np.ndarray,
              apply_to: np.ndarray | None = None,
              weights: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Superpose and return ``(transformed_coordinates, rmsd)``."""
    rot, ct, cm = kabsch(mobile, target, weights)
    fitted = (mobile - cm) @ rot.T + ct
    err = rmsd(fitted, target)
    if apply_to is not None:
        return (np.asarray(apply_to, dtype=np.float64) - cm) @ rot.T + ct, err
    return fitted, err


def rmsd(a: np.ndarray, b: np.ndarray) -> float:
    d = np.asarray(a, dtype=np.float64) - np.asarray(b, dtype=np.float64)
    return float(np.sqrt((d * d).sum(axis=1).mean()))


def rotation_matrix(axis: np.ndarray, angle: float) -> np.ndarray:
    """Rodrigues rotation matrix for ``angle`` radians about a unit ``axis``."""
    axis = np.asarray(axis, dtype=np.float64)
    axis = axis / np.linalg.norm(axis)
    x, y, z = axis
    c, s = np.cos(angle), np.sin(angle)
    C = 1.0 - c
    return np.array([
        [x * x * C + c,     x * y * C - z * s, x * z * C + y * s],
        [y * x * C + z * s, y * y * C + c,     y * z * C - x * s],
        [z * x * C - y * s, z * y * C + x * s, z * z * C + c],
    ])


@dataclass
class SymmetryAxis:
    """A molecular rotation axis: a point, a direction and its fit quality."""

    point: np.ndarray            # a point on the axis (3,)
    direction: np.ndarray        # unit vector (3,)
    order: int = 3
    rmsd: float = 0.0            # rmsd of the symmetry-mate superposition
    angle_deg: float = 0.0       # recovered rotation angle

    def project(self, xyz: np.ndarray) -> np.ndarray:
        """Signed distance of each point along the axis, relative to ``point``."""
        return (np.asarray(xyz) - self.point) @ self.direction

    def radial(self, xyz: np.ndarray) -> np.ndarray:
        """Perpendicular distance of each point from the axis."""
        d = np.asarray(xyz) - self.point
        along = np.outer(d @ self.direction, self.direction)
        return np.linalg.norm(d - along, axis=1)


def _axis_from_rotation(rot: np.ndarray) -> tuple[np.ndarray, float]:
    """Extract the rotation axis and angle from a 3x3 rotation matrix."""
    angle = np.arccos(np.clip((np.trace(rot) - 1.0) / 2.0, -1.0, 1.0))
    # Eigenvector with eigenvalue +1 is the axis; use SVD of (R - I) for stability.
    w, v = np.linalg.eig(rot)
    idx = int(np.argmin(np.abs(w - 1.0)))
    axis = np.real(v[:, idx])
    axis = axis / np.linalg.norm(axis)
    # Fix the sign so the rotation is right-handed about ``axis``.
    skew = np.array([rot[2, 1] - rot[1, 2], rot[0, 2] - rot[2, 0], rot[1, 0] - rot[0, 1]])
    if np.dot(skew, axis) < 0:
        axis = -axis
    return axis, float(np.degrees(angle))


def detect_c3_axis(chain_coords: list[np.ndarray]) -> SymmetryAxis:
    """Recover the three-fold axis from three equivalent chains.

    ``chain_coords`` must be a list of three ``(n, 3)`` arrays holding the same
    atoms, in the same order, for each protomer.  The axis is obtained from the
    rotation that maps protomer 0 onto protomer 1, and its position is fixed by
    requiring that the three chain centroids lie on a circle around it.
    """
    if len(chain_coords) != 3:
        raise ValueError("detect_c3_axis expects exactly three chains")
    n = min(len(c) for c in chain_coords)
    a, b = chain_coords[0][:n], chain_coords[1][:n]
    rot, ct, cm = kabsch(a, b)
    fitted = (a - cm) @ rot.T + ct
    err = rmsd(fitted, b)
    axis, angle = _axis_from_rotation(rot)

    # The axis passes through the centroid of the three chain centroids.
    centroids = np.array([c[:n].mean(axis=0) for c in chain_coords])
    point = centroids.mean(axis=0)

    # Orient +axis away from the centroid of the intracellular-most atoms so
    # that the sign convention is stable; callers may flip it explicitly.
    return SymmetryAxis(point=point, direction=axis, order=3, rmsd=err, angle_deg=angle)


@dataclass
class ProtomerMatch:
    """How the protomers of one trimer correspond to those of another."""

    order: tuple[int, ...]        # index into the mobile set for each target slot
    rmsd: float                   # trimer RMSD under this correspondence
    all_rmsd: dict                # every correspondence tried, for diagnostics
    handedness_flipped: bool      # True if the cyclic order had to be reversed


def match_protomers(mobile_blocks: list[np.ndarray],
                    target_blocks: list[np.ndarray]) -> ProtomerMatch:
    """Find which protomer of ``mobile`` corresponds to which of ``target``.

    Deposited chain labels are **not** a reliable guide: two structures of the
    same C3 trimer can label their chains in opposite rotational order around
    the symmetry axis. Superposing them by label then produces a nonsensical
    RMSD (we saw 71 A between PDB 7WLT and 7WLU, versus 19.7 A once the
    correspondence was fixed), and any difference vector derived from it is
    meaningless.

    Only the three cyclic rotations and the three reversed ones are chemically
    sensible, so all six are tried and the best is returned.
    """
    if len(mobile_blocks) != 3 or len(target_blocks) != 3:
        raise ValueError("match_protomers expects two sets of three protomers")
    if {len(b) for b in mobile_blocks + target_blocks} != {len(target_blocks[0])}:
        raise ValueError("all protomer blocks must have the same site count")

    target = np.vstack(target_blocks)
    candidates = [(0, 1, 2), (1, 2, 0), (2, 0, 1),      # same handedness
                  (0, 2, 1), (2, 1, 0), (1, 0, 2)]      # reversed
    scores: dict[tuple[int, ...], float] = {}
    for perm in candidates:
        mobile = np.vstack([mobile_blocks[i] for i in perm])
        _, err = superpose(mobile, target)
        scores[perm] = err

    best = min(scores, key=scores.get)  # type: ignore[arg-type]
    return ProtomerMatch(order=best, rmsd=scores[best],
                         all_rmsd={k: round(v, 3) for k, v in scores.items()},
                         handedness_flipped=best in candidates[3:])


def align_axis_to_z(axis: SymmetryAxis, flip: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Build ``(rotation, translation)`` putting ``axis`` on +z through the origin.

    Apply as ``xyz_new = (xyz + translation) @ rotation.T``.
    """
    d = np.asarray(axis.direction, dtype=np.float64)
    d = d / np.linalg.norm(d)
    if flip:
        d = -d
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(d, z)
    s = np.linalg.norm(v)
    if s < 1e-9:
        rot = np.eye(3) if d[2] > 0 else np.diag([1.0, -1.0, -1.0])
    else:
        c = float(np.dot(d, z))
        vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
        rot = np.eye(3) + vx + vx @ vx * ((1 - c) / (s * s))
    return rot, -np.asarray(axis.point, dtype=np.float64)
