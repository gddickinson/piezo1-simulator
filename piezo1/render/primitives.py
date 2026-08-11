"""GPU draw batches: impostor spheres, impostor cylinders and triangle meshes.

Each batch owns its OpenGL buffers and knows how to upload numpy arrays into
them. Uploading is separated from drawing so that an animation can push new
coordinates every frame without rebuilding vertex arrays.

All three batch types share the uniform block applied by
:class:`piezo1.render.scene.Scene`, so lighting and depth cueing stay
consistent between representations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

try:
    import moderngl
except ImportError:  # pragma: no cover - the GUI cannot run without it
    moderngl = None  # type: ignore[assignment]

__all__ = ["ShaderLibrary", "SphereBatch", "CylinderBatch", "MeshBatch", "Batch"]

SHADER_DIR = Path(__file__).resolve().parent / "shaders"


class ShaderLibrary:
    """Loads and caches GLSL programs from ``render/shaders``."""

    def __init__(self, ctx) -> None:
        self.ctx = ctx
        self._cache: dict[str, "moderngl.Program"] = {}

    def get(self, name: str):
        if name not in self._cache:
            vert = (SHADER_DIR / f"{name}.vert").read_text()
            frag = (SHADER_DIR / f"{name}.frag").read_text()
            self._cache[name] = self.ctx.program(vertex_shader=vert,
                                                 fragment_shader=frag)
        return self._cache[name]

    def release(self) -> None:
        for prog in self._cache.values():
            prog.release()
        self._cache.clear()


@dataclass
class Batch:
    """Common state for a drawable batch."""

    name: str = "batch"
    visible: bool = True
    count: int = 0
    ctx: object = field(default=None, repr=False)
    program: object = field(default=None, repr=False)
    vbo: object = field(default=None, repr=False)
    vao: object = field(default=None, repr=False)

    def set_uniforms(self, uniforms: dict) -> None:
        """Apply the shared uniform dictionary, skipping absent names."""
        prog = self.program
        if prog is None:
            return
        for key, value in uniforms.items():
            member = prog.get(key, None)
            if member is None:
                continue
            try:
                member.value = value
            except Exception:
                # A uniform optimised out of one program but present in another
                # is normal; never let it break the frame.
                pass

    def release(self) -> None:
        for obj in (self.vao, self.vbo):
            if obj is not None:
                try:
                    obj.release()
                except Exception:
                    pass
        self.vao = self.vbo = None
        self.count = 0


# --------------------------------------------------------------------------
# Spheres
# --------------------------------------------------------------------------

class SphereBatch(Batch):
    """Instanced sphere impostors: 4 vertices per atom, ray-cast per pixel."""

    #: interleaved per-instance record: centre, radius, colour, flags
    DTYPE = np.dtype([("center", "f4", 3), ("radius", "f4"),
                      ("color", "f4", 3), ("flags", "f4")])

    def __init__(self, ctx, library: ShaderLibrary, name: str = "spheres") -> None:
        super().__init__(name=name, ctx=ctx)
        self.program = library.get("sphere")
        self._capacity = 0

    def upload(self, centers: np.ndarray, radii: np.ndarray,
               colors: np.ndarray, flags: np.ndarray | None = None) -> None:
        n = len(centers)
        data = np.empty(n, dtype=self.DTYPE)
        data["center"] = centers
        data["radius"] = radii
        data["color"] = colors
        data["flags"] = 0.0 if flags is None else flags

        if self.vbo is None or n > self._capacity:
            self.release_buffers()
            self.vbo = self.ctx.buffer(data.tobytes(), dynamic=True)
            self._capacity = n
            self.vao = self.ctx.vertex_array(
                self.program,
                [(self.vbo, "3f 1f 3f 1f/i",
                  "in_center", "in_radius", "in_color", "in_flags")],
            )
        else:
            self.vbo.write(data.tobytes())
        self.count = n

    def update_centers(self, centers: np.ndarray) -> None:
        """Fast path for animation: rewrite only the position field."""
        if self.vbo is None or self.count != len(centers):
            raise RuntimeError("upload() must precede update_centers()")
        stride = self.DTYPE.itemsize
        buf = np.frombuffer(self.vbo.read(), dtype=self.DTYPE).copy()
        buf["center"] = centers
        self.vbo.write(buf.tobytes())
        del stride

    def release_buffers(self) -> None:
        super().release()
        self._capacity = 0

    def render(self) -> None:
        if self.visible and self.count and self.vao is not None:
            self.vao.render(moderngl.TRIANGLE_STRIP, vertices=4, instances=self.count)


# --------------------------------------------------------------------------
# Cylinders
# --------------------------------------------------------------------------

class CylinderBatch(Batch):
    """Instanced cylinder impostors with per-end colours.

    **Never culled.** The bounding quad is built around the cylinder's own axis,
    so which way it faces depends on the direction of the bond rather than on
    the camera — and with back-face culling on, the whole quad is thrown away
    before the fragment stage ever runs. That is why ball-and-stick drew balls
    only and the HaloTag seam drew nothing: the impostors were being culled,
    not mis-shaded. Culling is meaningless for an impostor in any case, because
    the fragment shader ray-casts the real surface and writes its own depth.
    """

    #: Read by :meth:`piezo1.render.scene.Scene.render`.
    cull = False

    DTYPE = np.dtype([("start", "f4", 3), ("end", "f4", 3), ("radius", "f4"),
                      ("color_a", "f4", 3), ("color_b", "f4", 3)])

    def __init__(self, ctx, library: ShaderLibrary, name: str = "cylinders") -> None:
        super().__init__(name=name, ctx=ctx)
        self.program = library.get("cylinder")
        self._capacity = 0

    def upload(self, starts: np.ndarray, ends: np.ndarray, radii: np.ndarray,
               colors_a: np.ndarray, colors_b: np.ndarray | None = None) -> None:
        n = len(starts)
        data = np.empty(n, dtype=self.DTYPE)
        data["start"] = starts
        data["end"] = ends
        data["radius"] = radii
        data["color_a"] = colors_a
        data["color_b"] = colors_a if colors_b is None else colors_b

        if self.vbo is None or n > self._capacity:
            self.release_buffers()
            self.vbo = self.ctx.buffer(data.tobytes(), dynamic=True)
            self._capacity = n
            self.vao = self.ctx.vertex_array(
                self.program,
                [(self.vbo, "3f 3f 1f 3f 3f/i",
                  "in_start", "in_end", "in_radius", "in_color_a", "in_color_b")],
            )
        else:
            self.vbo.write(data.tobytes())
        self.count = n

    def release_buffers(self) -> None:
        super().release()
        self._capacity = 0

    def render(self) -> None:
        if self.visible and self.count and self.vao is not None:
            self.vao.render(moderngl.TRIANGLE_STRIP, vertices=4, instances=self.count)


# --------------------------------------------------------------------------
# Meshes
# --------------------------------------------------------------------------

class MeshBatch(Batch):
    """Indexed triangle mesh with per-vertex colour and alpha."""

    DTYPE = np.dtype([("position", "f4", 3), ("normal", "f4", 3),
                      ("color", "f4", 3), ("alpha", "f4")])

    def __init__(self, ctx, library: ShaderLibrary, name: str = "mesh",
                 two_sided: bool = True, transparent: bool = False) -> None:
        super().__init__(name=name, ctx=ctx)
        self.program = library.get("mesh")
        self.two_sided = two_sided
        self.transparent = transparent
        self.ibo = None
        self._n_index = 0

    def upload(self, positions: np.ndarray, normals: np.ndarray,
               colors: np.ndarray, indices: np.ndarray,
               alpha: np.ndarray | float = 1.0) -> None:
        n = len(positions)
        data = np.empty(n, dtype=self.DTYPE)
        data["position"] = positions
        data["normal"] = normals
        data["color"] = colors
        data["alpha"] = alpha

        self.release()
        self.vbo = self.ctx.buffer(data.tobytes())
        idx = np.ascontiguousarray(indices, dtype="i4").ravel()
        self.ibo = self.ctx.buffer(idx.tobytes())
        self._n_index = len(idx)
        self.vao = self.ctx.vertex_array(
            self.program,
            [(self.vbo, "3f 3f 3f 1f", "in_position", "in_normal",
              "in_color", "in_alpha")],
            index_buffer=self.ibo,
            index_element_size=4,
        )
        self.count = n

    def release(self) -> None:
        if self.ibo is not None:
            try:
                self.ibo.release()
            except Exception:
                pass
            self.ibo = None
        self._n_index = 0
        super().release()

    def render(self) -> None:
        if self.visible and self._n_index and self.vao is not None:
            self.program["u_two_sided"].value = 1 if self.two_sided else 0
            self.vao.render(moderngl.TRIANGLES)
