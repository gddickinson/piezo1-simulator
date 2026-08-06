"""Trackball camera for the molecular viewport.

Orientation is held as a quaternion rather than Euler angles so that repeated
dragging never gimbal-locks and never accumulates roll. The camera orbits a
pivot point (normally the molecular centre) at a settable distance, which is
the interaction model every molecular viewer uses.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

__all__ = ["Camera", "quat_multiply", "quat_to_matrix", "look_at", "perspective"]


# --------------------------------------------------------------------------
# Quaternion helpers
# --------------------------------------------------------------------------

def quat_multiply(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Hamilton product of two ``(w, x, y, z)`` quaternions."""
    aw, ax, ay, az = a
    bw, bx, by, bz = b
    return np.array([
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ], dtype=np.float64)


def quat_to_matrix(q: np.ndarray) -> np.ndarray:
    """Rotation matrix from a unit ``(w, x, y, z)`` quaternion."""
    q = q / np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
        [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
        [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
    ], dtype=np.float64)


def quat_from_axis_angle(axis: np.ndarray, angle: float) -> np.ndarray:
    axis = np.asarray(axis, dtype=np.float64)
    n = np.linalg.norm(axis)
    if n < 1e-12:
        return np.array([1.0, 0.0, 0.0, 0.0])
    axis = axis / n
    s = np.sin(angle / 2.0)
    return np.array([np.cos(angle / 2.0), axis[0] * s, axis[1] * s, axis[2] * s])


# --------------------------------------------------------------------------
# Matrices
# --------------------------------------------------------------------------

def perspective(fov_deg: float, aspect: float, near: float, far: float) -> np.ndarray:
    f = 1.0 / np.tan(np.radians(fov_deg) / 2.0)
    m = np.zeros((4, 4), dtype=np.float64)
    m[0, 0] = f / max(aspect, 1e-6)
    m[1, 1] = f
    m[2, 2] = (far + near) / (near - far)
    m[2, 3] = (2.0 * far * near) / (near - far)
    m[3, 2] = -1.0
    return m


def orthographic(half_height: float, aspect: float, near: float, far: float) -> np.ndarray:
    h = max(half_height, 1e-6)
    w = h * max(aspect, 1e-6)
    m = np.eye(4, dtype=np.float64)
    m[0, 0] = 1.0 / w
    m[1, 1] = 1.0 / h
    m[2, 2] = -2.0 / (far - near)
    m[2, 3] = -(far + near) / (far - near)
    return m


def look_at(eye: np.ndarray, target: np.ndarray, up: np.ndarray) -> np.ndarray:
    f = target - eye
    f = f / np.linalg.norm(f)
    s = np.cross(f, up)
    ns = np.linalg.norm(s)
    if ns < 1e-9:
        up = np.array([0.0, 0.0, 1.0])
        s = np.cross(f, up)
        ns = np.linalg.norm(s)
    s = s / ns
    u = np.cross(s, f)
    m = np.eye(4, dtype=np.float64)
    m[0, :3], m[1, :3], m[2, :3] = s, u, -f
    m[:3, 3] = -m[:3, :3] @ eye
    return m


# --------------------------------------------------------------------------
# Camera
# --------------------------------------------------------------------------

@dataclass
class Camera:
    """Orbiting trackball camera."""

    pivot: np.ndarray = field(default_factory=lambda: np.zeros(3))
    distance: float = 300.0
    rotation: np.ndarray = field(default_factory=lambda: np.array([1.0, 0.0, 0.0, 0.0]))
    fov: float = 35.0
    aspect: float = 1.0
    orthographic: bool = False
    #: Near/far are recomputed from the scene radius each frame.
    scene_radius: float = 200.0
    pan: np.ndarray = field(default_factory=lambda: np.zeros(3))

    # ---------------------------------------------------------------- setup

    def frame(self, coords: np.ndarray, margin: float = 1.06) -> "Camera":
        """Point the camera at a coordinate set and pull back to include it.

        The pull-back distance uses the *true* bounding-sphere radius, not half
        the bounding-box diagonal. PIEZO1 is a wide, flat propeller, and the
        diagonal overestimates its radius badly enough to leave the molecule
        floating in the middle of an empty viewport.
        """
        coords = np.asarray(coords, dtype=np.float64)
        if len(coords) == 0:
            return self
        lo, hi = coords.min(axis=0), coords.max(axis=0)
        self.pivot = 0.5 * (lo + hi)
        self.scene_radius = float(np.linalg.norm(coords - self.pivot, axis=1).max()) or 100.0

        # Project into the *current* camera orientation and solve for the
        # distance that just contains everything. For a point at camera-frame
        # offset (x, y, z) from the pivot, staying inside the frustum requires
        #     d >= z + |x| / tan(half_horizontal)   and the same for y,
        # so the tight distance is the maximum of that over all points. Framing
        # by bounding sphere instead leaves a flat molecule like PIEZO1 filling
        # barely half the viewport.
        local = (coords - self.pivot) @ quat_to_matrix(self.rotation).T
        half_v = np.radians(self.fov) / 2.0
        half_h = np.arctan(np.tan(half_v) * max(self.aspect, 1e-3))
        need_x = local[:, 2] + np.abs(local[:, 0]) / max(np.tan(half_h), 1e-6)
        need_y = local[:, 2] + np.abs(local[:, 1]) / max(np.tan(half_v), 1e-6)
        self.distance = float(max(need_x.max(), need_y.max()) * margin)
        self.distance = max(self.distance, self.scene_radius * 0.2)
        self.pan = np.zeros(3)
        return self

    # ----------------------------------------------------------- navigation

    def orbit(self, dx: float, dy: float, speed: float = 4.0) -> None:
        """Rotate by a normalised drag in screen space."""
        if dx == 0.0 and dy == 0.0:
            return
        # Rotate about camera-space axes so dragging always feels the same.
        q_yaw = quat_from_axis_angle(np.array([0.0, 1.0, 0.0]), dx * speed)
        q_pitch = quat_from_axis_angle(np.array([1.0, 0.0, 0.0]), dy * speed)
        self.rotation = quat_multiply(quat_multiply(q_yaw, q_pitch), self.rotation)
        self.rotation /= np.linalg.norm(self.rotation)

    def roll(self, angle: float) -> None:
        q = quat_from_axis_angle(np.array([0.0, 0.0, 1.0]), angle)
        self.rotation = quat_multiply(q, self.rotation)
        self.rotation /= np.linalg.norm(self.rotation)

    def zoom(self, factor: float) -> None:
        self.distance = float(np.clip(self.distance * factor,
                                      self.scene_radius * 0.05,
                                      self.scene_radius * 40.0))

    def translate(self, dx: float, dy: float) -> None:
        """Pan in the camera plane; ``dx``/``dy`` are fractions of the viewport."""
        scale = self.distance * np.tan(np.radians(self.fov) / 2.0) * 2.0
        rot = quat_to_matrix(self.rotation)
        right, up = rot[0], rot[1]
        self.pan = self.pan - right * dx * scale * self.aspect - up * dy * scale

    def spin(self, degrees: float) -> None:
        """Rotate about the vertical screen axis — used for the auto-turntable."""
        self.orbit(np.radians(degrees) / 4.0, 0.0)

    # -------------------------------------------------------------- matrices

    @property
    def eye(self) -> np.ndarray:
        rot = quat_to_matrix(self.rotation)
        return self.pivot + self.pan + rot[2] * self.distance

    def view_matrix(self) -> np.ndarray:
        rot = quat_to_matrix(self.rotation)
        m = np.eye(4, dtype=np.float64)
        m[:3, :3] = rot
        m[:3, 3] = -rot @ (self.pivot + self.pan)
        m[2, 3] -= self.distance
        return m

    def clip_planes(self) -> tuple[float, float]:
        near = max(self.distance - self.scene_radius * 2.5, self.scene_radius * 0.01)
        far = self.distance + self.scene_radius * 2.5
        return near, far

    def projection_matrix(self) -> np.ndarray:
        near, far = self.clip_planes()
        if self.orthographic:
            half = self.distance * np.tan(np.radians(self.fov) / 2.0)
            return orthographic(half, self.aspect, near, far)
        return perspective(self.fov, self.aspect, near, far)

    def matrices(self) -> tuple[np.ndarray, np.ndarray]:
        return self.view_matrix(), self.projection_matrix()

    # ---------------------------------------------------------------- picking

    def ray_through_pixel(self, x_ndc: float, y_ndc: float) -> tuple[np.ndarray, np.ndarray]:
        """World-space ``(origin, direction)`` for a point in NDC [-1, 1]."""
        proj = self.projection_matrix()
        view = self.view_matrix()
        inv = np.linalg.inv(proj @ view)
        near = inv @ np.array([x_ndc, y_ndc, -1.0, 1.0])
        far = inv @ np.array([x_ndc, y_ndc, 1.0, 1.0])
        near = near[:3] / near[3]
        far = far[:3] / far[3]
        d = far - near
        return near, d / np.linalg.norm(d)
