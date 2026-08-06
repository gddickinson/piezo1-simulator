"""Backbone splines, rotation-minimising frames and secondary structure.

Cartoon rendering needs three things from a chain of C-alpha positions: a
smooth curve through them, a coordinate frame that slides along that curve
without spinning, and a per-residue secondary-structure label to set the ribbon
profile.

The frame is built by **parallel transport** rather than from carbonyl vectors.
The traditional approach — orienting the ribbon with the peptide C=O — suffers
the well-known "flip" problem, where successive carbonyls alternate direction
in strands and the ribbon shears 180 degrees between residues; viewers then
patch it with an ad-hoc sign correction. Parallel transport sidesteps the whole
issue, needs only C-alpha coordinates (so it works on the many poly-alanine and
backbone-only regions of PIEZO1 cryo-EM models), and produces a demonstrably
untwisted ribbon.
"""

from __future__ import annotations

import numpy as np

__all__ = ["catmull_rom", "resample_uniform", "parallel_transport_frames",
           "assign_secondary_structure", "smooth_path", "SS_COIL", "SS_HELIX",
           "SS_STRAND"]

SS_COIL, SS_HELIX, SS_STRAND = 0, 1, 2


# --------------------------------------------------------------------------
# Splines
# --------------------------------------------------------------------------

def catmull_rom(points: np.ndarray, subdivisions: int = 6,
                tension: float = 0.5) -> np.ndarray:
    """Catmull-Rom spline through ``points``.

    Interpolating (the curve passes through every control point) and local
    (moving one atom only disturbs the ribbon nearby), which is what a
    molecular cartoon wants. Endpoints are handled by reflecting the second
    and penultimate points.
    """
    p = np.asarray(points, dtype=np.float64)
    n = len(p)
    if n < 2:
        return p.copy()
    if n == 2:
        t = np.linspace(0.0, 1.0, subdivisions + 1)[:, None]
        return p[0] + t * (p[1] - p[0])

    ext = np.vstack([2 * p[0] - p[1], p, 2 * p[-1] - p[-2]])
    t = np.linspace(0.0, 1.0, subdivisions, endpoint=False)[:, None]
    t2, t3 = t * t, t * t * t

    p0, p1, p2, p3 = ext[:-3], ext[1:-2], ext[2:-1], ext[3:]
    m1 = tension * (p2 - p0)
    m2 = tension * (p3 - p1)

    # Hermite basis, evaluated for every segment at once.
    h00 = (2 * t3 - 3 * t2 + 1)[None, :, :]
    h10 = (t3 - 2 * t2 + t)[None, :, :]
    h01 = (-2 * t3 + 3 * t2)[None, :, :]
    h11 = (t3 - t2)[None, :, :]

    seg = (h00 * p1[:, None, :] + h10 * m1[:, None, :]
           + h01 * p2[:, None, :] + h11 * m2[:, None, :])
    out = seg.reshape(-1, 3)
    return np.vstack([out, p[-1]])


def smooth_path(points: np.ndarray, passes: int = 2, weight: float = 0.5) -> np.ndarray:
    """Laplacian smoothing, used to tame the zig-zag of a raw C-alpha trace."""
    p = np.array(points, dtype=np.float64, copy=True)
    for _ in range(passes):
        if len(p) < 3:
            break
        mid = 0.5 * (p[:-2] + p[2:])
        p[1:-1] = (1.0 - weight) * p[1:-1] + weight * mid
    return p


def resample_uniform(path: np.ndarray, spacing: float) -> np.ndarray:
    """Resample a polyline to approximately constant arc-length spacing."""
    p = np.asarray(path, dtype=np.float64)
    if len(p) < 2:
        return p.copy()
    seg = np.linalg.norm(np.diff(p, axis=0), axis=1)
    arc = np.concatenate([[0.0], np.cumsum(seg)])
    total = arc[-1]
    if total <= 0:
        return p.copy()
    n = max(int(np.ceil(total / max(spacing, 1e-6))), 1)
    target = np.linspace(0.0, total, n + 1)
    out = np.empty((len(target), 3))
    for k in range(3):
        out[:, k] = np.interp(target, arc, p[:, k])
    return out


# --------------------------------------------------------------------------
# Frames
# --------------------------------------------------------------------------

def parallel_transport_frames(path: np.ndarray,
                              initial_normal: np.ndarray | None = None
                              ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Rotation-minimising frames along a polyline.

    Returns ``(tangent, normal, binormal)``, each ``(n, 3)``. The normal is
    carried forward by the minimal rotation that maps one tangent onto the
    next, so the frame never spins about the curve.
    """
    p = np.asarray(path, dtype=np.float64)
    n = len(p)
    if n < 2:
        z = np.zeros((n, 3))
        return z, z, z

    tangent = np.empty_like(p)
    tangent[1:-1] = p[2:] - p[:-2]
    tangent[0] = p[1] - p[0]
    tangent[-1] = p[-1] - p[-2]
    lengths = np.linalg.norm(tangent, axis=1, keepdims=True)
    tangent = tangent / np.maximum(lengths, 1e-9)

    normal = np.empty_like(p)
    if initial_normal is None:
        seed = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(seed, tangent[0])) > 0.9:
            seed = np.array([0.0, 1.0, 0.0])
    else:
        seed = np.asarray(initial_normal, dtype=np.float64)
    v = seed - tangent[0] * np.dot(seed, tangent[0])
    nv = np.linalg.norm(v)
    normal[0] = v / nv if nv > 1e-9 else np.array([0.0, 0.0, 1.0])

    for i in range(1, n):
        t0, t1 = tangent[i - 1], tangent[i]
        axis = np.cross(t0, t1)
        s = np.linalg.norm(axis)
        if s < 1e-9:
            normal[i] = normal[i - 1]
        else:
            axis = axis / s
            angle = np.arctan2(s, float(np.dot(t0, t1)))
            c, sn = np.cos(angle), np.sin(angle)
            v = normal[i - 1]
            # Rodrigues rotation of the previous normal onto the new tangent.
            normal[i] = (v * c + np.cross(axis, v) * sn
                         + axis * np.dot(axis, v) * (1.0 - c))
        normal[i] -= tangent[i] * np.dot(normal[i], tangent[i])
        ln = np.linalg.norm(normal[i])
        normal[i] = normal[i] / ln if ln > 1e-9 else normal[i - 1]

    binormal = np.cross(tangent, normal)
    return tangent, normal, binormal


# --------------------------------------------------------------------------
# Secondary structure from C-alpha geometry alone
# --------------------------------------------------------------------------

def assign_secondary_structure(ca: np.ndarray, smooth: int = 2) -> np.ndarray:
    """Label each residue coil/helix/strand from C-alpha geometry.

    A P-SEA-style rule set: alpha helices and beta strands have characteristic
    i to i+3 and i to i+4 C-alpha distances and a characteristic C-alpha
    pseudo-torsion, and those three numbers separate the two cleanly without
    needing backbone hydrogen bonds. This matters here because large parts of
    the PIEZO1 cryo-EM models are backbone-only, where a hydrogen-bond method
    such as DSSP cannot run at all.

    Reference geometry (Labesse et al. 1997):
        helix   d13 ~ 5.3 A, d14 ~ 6.2 A, torsion ~ +50 deg
        strand  d13 ~ 10.4 A, d14 ~ 13.0 A, torsion ~ -170 deg
    """
    ca = np.asarray(ca, dtype=np.float64)
    n = len(ca)
    ss = np.full(n, SS_COIL, dtype=np.int8)
    if n < 5:
        return ss

    def dist(k: int) -> np.ndarray:
        d = np.full(n, np.nan)
        d[: n - k] = np.linalg.norm(ca[k:] - ca[:-k], axis=1)
        return d

    d2, d3, d4 = dist(2), dist(3), dist(4)

    # C-alpha pseudo-torsion over four consecutive residues.
    tor = np.full(n, np.nan)
    if n >= 4:
        b1 = ca[1:-2] - ca[:-3]
        b2 = ca[2:-1] - ca[1:-2]
        b3 = ca[3:] - ca[2:-1]
        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)
        m1 = np.cross(n1, b2 / np.linalg.norm(b2, axis=1, keepdims=True))
        x = np.einsum("ij,ij->i", n1, n2)
        y = np.einsum("ij,ij->i", m1, n2)
        tor[: n - 3] = np.degrees(np.arctan2(y, x))

    helix = (
        (np.abs(d3 - 5.3) < 0.6) & (np.abs(d4 - 6.4) < 1.0)
        & (np.abs(d2 - 5.5) < 0.8)
        & (np.abs(tor - 50.0) < 30.0)
    )
    strand = (
        (np.abs(d3 - 10.4) < 1.3) & (np.abs(d4 - 13.0) < 1.8)
        & (np.abs(d2 - 6.7) < 0.7)
        & ((np.abs(tor - 195.0) < 45.0) | (np.abs(tor + 170.0) < 45.0))
    )

    # A residue's label applies to the window it starts, so paint it forward.
    for i in np.flatnonzero(np.nan_to_num(helix)):
        ss[i:min(i + 4, n)] = SS_HELIX
    for i in np.flatnonzero(np.nan_to_num(strand)):
        ss[i:min(i + 3, n)] = np.where(ss[i:min(i + 3, n)] == SS_HELIX,
                                       SS_HELIX, SS_STRAND)

    for _ in range(smooth):
        ss = _despeckle(ss)
    return ss


def _despeckle(ss: np.ndarray, min_run: int = 3) -> np.ndarray:
    """Remove runs shorter than ``min_run``; isolated labels are noise."""
    out = ss.copy()
    n = len(ss)
    start = 0
    for i in range(1, n + 1):
        if i == n or ss[i] != ss[start]:
            if i - start < min_run and ss[start] != SS_COIL:
                out[start:i] = SS_COIL
            start = i
    return out
