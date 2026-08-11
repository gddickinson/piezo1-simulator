"""Mesh construction: swept tubes, cartoon ribbons and the membrane surface.

Every mesh is returned as ``(positions, normals, colors, indices)`` ready to be
handed to :class:`piezo1.render.primitives.MeshBatch`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .spline import (SS_COIL, SS_HELIX, SS_STRAND, catmull_rom,
                     parallel_transport_frames, smooth_path)

__all__ = ["Mesh", "build_tube", "build_cartoon", "build_membrane_mesh",
           "build_disc", "build_sphere", "CARTOON_PROFILE"]


@dataclass
class Mesh:
    positions: np.ndarray
    normals: np.ndarray
    colors: np.ndarray
    indices: np.ndarray

    @property
    def n_vertices(self) -> int:
        return len(self.positions)

    @property
    def n_triangles(self) -> int:
        return len(self.indices) // 3

    @classmethod
    def empty(cls) -> "Mesh":
        return cls(np.zeros((0, 3), np.float32), np.zeros((0, 3), np.float32),
                   np.zeros((0, 3), np.float32), np.zeros(0, np.int32))

    def concat(self, other: "Mesh") -> "Mesh":
        if other.n_vertices == 0:
            return self
        if self.n_vertices == 0:
            return other
        return Mesh(
            np.vstack([self.positions, other.positions]),
            np.vstack([self.normals, other.normals]),
            np.vstack([self.colors, other.colors]),
            np.concatenate([self.indices, other.indices + self.n_vertices]),
        )


#: Cross-section half-width and half-height, in Angstrom, per secondary
#: structure class. Helices and strands get a flat ribbon; coil gets a thin
#: round tube.
CARTOON_PROFILE = {
    SS_COIL: (0.35, 0.35),
    SS_HELIX: (1.30, 0.28),
    SS_STRAND: (1.25, 0.25),
}


def _cross_section(sides: int) -> np.ndarray:
    """Unit ellipse sample points in the (normal, binormal) plane."""
    theta = np.linspace(0.0, 2.0 * np.pi, sides, endpoint=False)
    return np.stack([np.cos(theta), np.sin(theta)], axis=1)


def _sweep(path: np.ndarray, half_width: np.ndarray, half_height: np.ndarray,
           colors: np.ndarray, sides: int = 10,
           cap: bool = True) -> Mesh:
    """Sweep an elliptical cross-section along a path."""
    n = len(path)
    if n < 2:
        return Mesh.empty()
    tangent, normal, binormal = parallel_transport_frames(path)
    unit = _cross_section(sides)

    # (n, sides, 3)
    offs = (unit[None, :, 0, None] * half_width[:, None, None] * normal[:, None, :]
            + unit[None, :, 1, None] * half_height[:, None, None] * binormal[:, None, :])
    positions = path[:, None, :] + offs
    norms = offs / np.maximum(np.linalg.norm(offs, axis=2, keepdims=True), 1e-9)

    positions = positions.reshape(-1, 3)
    norms = norms.reshape(-1, 3)
    vcolors = np.repeat(colors, sides, axis=0)

    # Quads between consecutive rings, split into triangles.
    i = np.arange(n - 1)[:, None] * sides
    j = np.arange(sides)[None, :]
    jn = (j + 1) % sides
    a = i + j
    b = i + jn
    c = i + sides + jn
    d = i + sides + j
    tris = np.concatenate([
        np.stack([a, b, c], axis=-1).reshape(-1, 3),
        np.stack([a, c, d], axis=-1).reshape(-1, 3),
    ]).astype(np.int32)

    mesh = Mesh(positions.astype(np.float32), norms.astype(np.float32),
                vcolors.astype(np.float32), tris.ravel())

    if cap:
        mesh = mesh.concat(_cap(path[0], -tangent[0], positions[:sides],
                                colors[0], flip=True))
        mesh = mesh.concat(_cap(path[-1], tangent[-1],
                                positions[(n - 1) * sides: n * sides],
                                colors[-1], flip=False))
    return mesh


def _cap(center: np.ndarray, normal: np.ndarray, ring: np.ndarray,
         color: np.ndarray, flip: bool) -> Mesh:
    sides = len(ring)
    positions = np.vstack([center[None, :], ring])
    normals = np.tile(normal, (sides + 1, 1))
    colors = np.tile(color, (sides + 1, 1))
    j = np.arange(sides)
    jn = (j + 1) % sides
    tris = (np.stack([np.zeros_like(j), j + 1, jn + 1], axis=-1)
            if not flip else
            np.stack([np.zeros_like(j), jn + 1, j + 1], axis=-1))
    return Mesh(positions.astype(np.float32), normals.astype(np.float32),
                colors.astype(np.float32), tris.ravel().astype(np.int32))


def build_tube(ca: np.ndarray, colors: np.ndarray, radius: float = 0.6,
               sides: int = 10, subdivisions: int = 6,
               smoothing: int = 1) -> Mesh:
    """A constant-radius tube through a C-alpha trace."""
    ca = np.asarray(ca, dtype=np.float64)
    if len(ca) < 2:
        return Mesh.empty()
    path = catmull_rom(smooth_path(ca, passes=smoothing), subdivisions)
    col = _interpolate_colors(colors, len(path))
    hw = np.full(len(path), radius)
    return _sweep(path, hw, hw, col, sides=sides)


def build_cartoon(ca: np.ndarray, ss: np.ndarray, colors: np.ndarray,
                  sides: int = 12, subdivisions: int = 6,
                  scale: float = 1.0, arrows: bool = True) -> Mesh:
    """Classic cartoon: flat ribbons for helices and strands, tube for coil.

    The ribbon profile is interpolated along the spline rather than switched
    abruptly at residue boundaries, so helix ends taper into coil instead of
    stepping.
    """
    ca = np.asarray(ca, dtype=np.float64)
    ss = np.asarray(ss)
    n = len(ca)
    if n < 3:
        return Mesh.empty()

    # Helices are smoothed harder: the raw C-alpha trace of a helix is a coil
    # of radius ~2.3 A, and a ribbon following it exactly looks like a spring.
    base = smooth_path(ca, passes=1)
    helix_smoothed = smooth_path(ca, passes=4, weight=0.65)
    is_helix = (ss == SS_HELIX)[:, None]
    control = np.where(is_helix, helix_smoothed, base)

    path = catmull_rom(control, subdivisions)
    m = len(path)
    idx = np.clip((np.arange(m) / max(m - 1, 1) * (n - 1)).astype(int), 0, n - 1)

    hw = np.array([CARTOON_PROFILE[int(s)][0] for s in ss]) * scale
    hh = np.array([CARTOON_PROFILE[int(s)][1] for s in ss]) * scale

    if arrows:
        hw = hw.copy()
        # Widen the last two residues of each strand into an arrowhead.
        strand = ss == SS_STRAND
        ends = np.flatnonzero(strand & ~np.roll(strand, -1))
        for e in ends:
            lo = max(e - 1, 0)
            hw[lo:e + 1] = np.linspace(2.4, 0.4, e + 1 - lo) * scale

    hw_path = np.interp(np.arange(m), np.linspace(0, m - 1, n), hw)
    hh_path = np.interp(np.arange(m), np.linspace(0, m - 1, n), hh)
    col = np.asarray(colors, dtype=np.float64)[idx]
    return _sweep(path, hw_path, hh_path, col, sides=sides)


def _interpolate_colors(colors: np.ndarray, m: int) -> np.ndarray:
    colors = np.asarray(colors, dtype=np.float64)
    n = len(colors)
    if n == m:
        return colors
    idx = np.clip((np.arange(m) / max(m - 1, 1) * (n - 1)).astype(int), 0, n - 1)
    return colors[idx]


# --------------------------------------------------------------------------
# Membrane
# --------------------------------------------------------------------------

def build_membrane_mesh(height_fn, r_max: float, n_radial: int = 96,
                        n_angular: int = 128, thickness: float = 0.0,
                        color: tuple[float, float, float] = (0.32, 0.45, 0.72),
                        axis: np.ndarray | None = None,
                        origin: np.ndarray | None = None) -> Mesh:
    """Surface of revolution ``z = height_fn(r)`` about a symmetry axis.

    Used both for the mid-bilayer surface of the PIEZO1 dome and for the two
    leaflet surfaces, by passing a non-zero ``thickness`` offset.
    """
    origin = np.zeros(3) if origin is None else np.asarray(origin, float)
    axis = np.array([0.0, 0.0, 1.0]) if axis is None else np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    e1 = np.cross(axis, [0.0, 0.0, 1.0])
    if np.linalg.norm(e1) < 1e-6:
        e1 = np.cross(axis, [0.0, 1.0, 0.0])
    e1 /= np.linalg.norm(e1)
    e2 = np.cross(axis, e1)

    r = np.linspace(0.0, r_max, n_radial)
    theta = np.linspace(0.0, 2.0 * np.pi, n_angular, endpoint=False)
    R, T = np.meshgrid(r, theta, indexing="ij")
    Z = height_fn(R) + thickness

    pts = (origin[None, None, :]
           + R[..., None] * (np.cos(T)[..., None] * e1 + np.sin(T)[..., None] * e2)
           + Z[..., None] * axis[None, None, :])

    # Analytic normal of a surface of revolution: dz/dr sets the tilt.
    dz = np.gradient(height_fn(r), r, edge_order=2)
    DZ = np.broadcast_to(dz[:, None], R.shape)
    radial = np.cos(T)[..., None] * e1 + np.sin(T)[..., None] * e2
    nrm = axis[None, None, :] - DZ[..., None] * radial
    nrm /= np.maximum(np.linalg.norm(nrm, axis=2, keepdims=True), 1e-9)

    positions = pts.reshape(-1, 3)
    normals = nrm.reshape(-1, 3)
    colors = np.tile(np.asarray(color, dtype=np.float64), (len(positions), 1))

    i = np.arange(n_radial - 1)[:, None] * n_angular
    j = np.arange(n_angular)[None, :]
    jn = (j + 1) % n_angular
    a, b = i + j, i + jn
    c, d = i + n_angular + jn, i + n_angular + j
    tris = np.concatenate([
        np.stack([a, b, c], axis=-1).reshape(-1, 3),
        np.stack([a, c, d], axis=-1).reshape(-1, 3),
    ]).astype(np.int32)

    return Mesh(positions.astype(np.float32), normals.astype(np.float32),
                colors.astype(np.float32), tris.ravel())


def build_sphere(center, radius: float, n_lat: int = 32, n_lon: int = 64,
                 color: tuple[float, float, float] = (0.4, 0.6, 0.9)) -> Mesh:
    """A UV sphere as a *mesh*, which the impostor sphere batch cannot be.

    The renderer already draws spheres far better than this — ray-cast
    impostors, pixel-exact at any zoom. What they cannot be is transparent:
    :class:`SphereBatch` has no alpha, because an atom never needed one. A
    field iso-surface does: it encloses the protein and would otherwise hide
    everything it is drawn around.
    """
    center = np.asarray(center, dtype=float)
    lat = np.linspace(0.0, np.pi, n_lat)
    lon = np.linspace(0.0, 2.0 * np.pi, n_lon, endpoint=False)
    LAT, LON = np.meshgrid(lat, lon, indexing="ij")
    normals = np.stack([np.sin(LAT) * np.cos(LON),
                        np.sin(LAT) * np.sin(LON),
                        np.cos(LAT)], axis=-1)
    positions = center[None, None, :] + radius * normals

    i = np.arange(n_lat - 1)[:, None] * n_lon
    j = np.arange(n_lon)[None, :]
    jn = (j + 1) % n_lon
    a, b = i + j, i + jn
    c, d = i + n_lon + jn, i + n_lon + j
    tris = np.concatenate([
        np.stack([a, b, c], axis=-1).reshape(-1, 3),
        np.stack([a, c, d], axis=-1).reshape(-1, 3),
    ]).astype(np.int32)

    flat = positions.reshape(-1, 3)
    return Mesh(flat.astype(np.float32),
                normals.reshape(-1, 3).astype(np.float32),
                np.tile(np.asarray(color, np.float32), (len(flat), 1)),
                tris.ravel())


def build_disc(radius: float, n_angular: int = 128,
               color: tuple[float, float, float] = (0.3, 0.4, 0.6),
               axis: np.ndarray | None = None,
               origin: np.ndarray | None = None) -> Mesh:
    """A flat disc, used as the undeformed reference membrane plane."""
    return build_membrane_mesh(lambda r: np.zeros_like(r), radius,
                               n_radial=2, n_angular=n_angular,
                               color=color, axis=axis, origin=origin)
