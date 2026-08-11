"""Scene graph and frame rendering.

A :class:`Scene` owns the GL context wrapper, the shader library, the camera
and an ordered set of named batches. The application adds and removes batches;
the scene handles uniforms, render order and blending.

Render order matters: opaque impostors and meshes first with depth writes on,
then transparent geometry (the membrane) back-to-front with depth writes off,
so the bilayer can be seen through without z-fighting against the protein.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from ..config import RenderSettings
from .camera import Camera
from .primitives import Batch, CylinderBatch, MeshBatch, ShaderLibrary, SphereBatch

try:
    import moderngl
except ImportError:  # pragma: no cover
    moderngl = None  # type: ignore[assignment]

__all__ = ["Scene", "Light"]


@dataclass
class Light:
    """A single directional key light, specified in camera space."""

    direction: tuple[float, float, float] = (0.35, 0.45, 0.82)
    ambient: float = 0.30
    shininess: float = 28.0


class Scene:
    """Holds every drawable and renders one frame."""

    def __init__(self, ctx, settings: RenderSettings | None = None) -> None:
        self.ctx = ctx
        self.settings = settings or RenderSettings()
        self.library = ShaderLibrary(ctx)
        self.camera = Camera(fov=self.settings.fov_degrees)
        self.light = Light()
        self.batches: dict[str, Batch] = {}
        self.order: list[str] = []
        self.radius_scale = 1.0
        self.viewport = (0, 0, 100, 100)
        self._configure_state()

    # ------------------------------------------------------------ GL state

    def _configure_state(self) -> None:
        self.ctx.enable(moderngl.DEPTH_TEST)
        self.ctx.enable(moderngl.CULL_FACE)
        self.ctx.cull_face = "back"
        self.ctx.blend_func = moderngl.SRC_ALPHA, moderngl.ONE_MINUS_SRC_ALPHA

    # ------------------------------------------------------------- batches

    def add(self, batch: Batch, name: str | None = None) -> Batch:
        key = name or batch.name
        if key in self.batches:
            self.remove(key)
        self.batches[key] = batch
        self.order.append(key)
        return batch

    def remove(self, name: str) -> None:
        batch = self.batches.pop(name, None)
        if batch is not None:
            batch.release()
        if name in self.order:
            self.order.remove(name)

    def clear(self) -> None:
        for name in list(self.order):
            self.remove(name)

    def get(self, name: str) -> Batch | None:
        return self.batches.get(name)

    def set_visible(self, name: str, visible: bool) -> None:
        batch = self.batches.get(name)
        if batch is not None:
            batch.visible = visible

    # ------------------------------------------------------------- factories

    def spheres(self, name: str) -> SphereBatch:
        return self.add(SphereBatch(self.ctx, self.library, name=name), name)

    def cylinders(self, name: str) -> CylinderBatch:
        return self.add(CylinderBatch(self.ctx, self.library, name=name), name)

    def mesh(self, name: str, two_sided: bool = True,
             transparent: bool = False) -> MeshBatch:
        return self.add(MeshBatch(self.ctx, self.library, name=name,
                                  two_sided=two_sided, transparent=transparent), name)

    # ------------------------------------------------------------- uniforms

    def uniforms(self) -> dict:
        view, proj = self.camera.matrices()
        near, far = self.camera.clip_planes()
        bg = self.settings.background
        # Fog spans the depth range actually occupied by the molecule, so the
        # cue strength does not change as the user zooms.
        d = self.camera.distance
        r = self.camera.scene_radius
        return {
            "u_view": tuple(view.T.astype("f4").ravel()),
            "u_proj": tuple(proj.T.astype("f4").ravel()),
            "u_normal_matrix": tuple(view[:3, :3].T.astype("f4").ravel()),
            "u_light_dir": self.light.direction,
            "u_ambient": self.light.ambient,
            "u_shininess": self.light.shininess,
            "u_fog_color": (bg[0], bg[1], bg[2]),
            "u_fog_near": float(max(d - r * 0.9, near)),
            "u_fog_far": float(min(d + r * 1.5, far)),
            "u_fog_strength": 0.55 if self.settings.depth_cue else 0.0,
            "u_radius_scale": float(self.radius_scale),
            "u_outline": 0,
        }

    # --------------------------------------------------------------- render

    def resize(self, width: int, height: int) -> None:
        width, height = max(width, 1), max(height, 1)
        self.viewport = (0, 0, width, height)
        self.camera.aspect = width / height

    def render(self) -> None:
        self.ctx.viewport = self.viewport
        self.ctx.clear(*self.settings.background)
        self._configure_state()
        uniforms = self.uniforms()

        opaque = [n for n in self.order
                  if not getattr(self.batches[n], "transparent", False)]
        transparent = [n for n in self.order
                       if getattr(self.batches[n], "transparent", False)]

        self.ctx.disable(moderngl.BLEND)
        for name in opaque:
            batch = self.batches[name]
            if batch.visible and batch.count:
                # Impostor batches declare `cull = False`: their bounding quad
                # is oriented by the geometry rather than by the camera, so
                # culling discards it outright. See `CylinderBatch`.
                if getattr(batch, "cull", True):
                    self.ctx.enable(moderngl.CULL_FACE)
                else:
                    self.ctx.disable(moderngl.CULL_FACE)
                batch.set_uniforms(uniforms)
                batch.render()
        self.ctx.enable(moderngl.CULL_FACE)

        if transparent:
            self.ctx.enable(moderngl.BLEND)
            self.ctx.depth_mask = False
            self.ctx.disable(moderngl.CULL_FACE)
            for name in transparent:
                batch = self.batches[name]
                if batch.visible and batch.count:
                    batch.set_uniforms(uniforms)
                    batch.render()
            self.ctx.depth_mask = True
            self.ctx.enable(moderngl.CULL_FACE)
            self.ctx.disable(moderngl.BLEND)

    # ---------------------------------------------------------------- utils

    def frame_all(self) -> None:
        """Point the camera at every sphere batch currently in the scene."""
        pts: list[np.ndarray] = []
        for batch in self.batches.values():
            if isinstance(batch, SphereBatch) and batch.count and batch.vbo:
                data = np.frombuffer(batch.vbo.read(), dtype=SphereBatch.DTYPE)
                pts.append(np.asarray(data["center"][:batch.count]))
        if pts:
            self.camera.frame(np.vstack(pts))

    def release(self) -> None:
        self.clear()
        self.library.release()
