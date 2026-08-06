"""The OpenGL viewport widget.

Hosts a moderngl context inside a ``QOpenGLWidget`` and forwards mouse and
wheel input to the scene camera.

One Qt-specific detail matters and is easy to get wrong: ``QOpenGLWidget`` does
**not** render to the window's default framebuffer. It renders into an FBO that
Qt owns and later composites. moderngl must therefore be pointed at
``defaultFramebufferObject()`` at the start of every frame, or the widget draws
into framebuffer 0 and the user sees nothing.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QPoint, Qt, QTimer, pyqtSignal
from PyQt6.QtGui import QSurfaceFormat
from PyQt6.QtOpenGLWidgets import QOpenGLWidget

from ..config import RenderSettings
from ..render.scene import Scene

__all__ = ["ViewportWidget", "configure_surface_format"]


def configure_surface_format(settings: RenderSettings | None = None) -> None:
    """Request an OpenGL core profile before any widget is created.

    Must be called before ``QApplication`` is constructed. macOS gives 4.1 core
    as its maximum, which is exactly what the impostor shaders need.
    """
    s = settings or RenderSettings()
    fmt = QSurfaceFormat()
    fmt.setVersion(s.gl_major, s.gl_minor)
    fmt.setProfile(QSurfaceFormat.OpenGLContextProfile.CoreProfile)
    fmt.setDepthBufferSize(24)
    fmt.setStencilBufferSize(8)
    fmt.setSamples(s.samples)
    fmt.setSwapBehavior(QSurfaceFormat.SwapBehavior.DoubleBuffer)
    QSurfaceFormat.setDefaultFormat(fmt)


class ViewportWidget(QOpenGLWidget):
    """Interactive 3D viewport."""

    #: Emitted once the GL context exists and a :class:`Scene` has been created.
    scene_ready = pyqtSignal(object)
    #: Emitted with an atom index (or -1) when the user clicks.
    atom_picked = pyqtSignal(int)
    #: Emitted with a human-readable status string.
    status = pyqtSignal(str)

    def __init__(self, settings: RenderSettings | None = None, parent=None) -> None:
        super().__init__(parent)
        self.settings = settings or RenderSettings()
        self.ctx = None
        self.scene: Scene | None = None
        self._last_pos: QPoint | None = None
        self._buttons = Qt.MouseButton.NoButton
        self._pick_source: np.ndarray | None = None
        self._spin_speed = 0.0
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._on_tick)
        self._timer.start(int(1000 / max(self.settings.target_fps, 1)))
        self._animations: list = []

    # ------------------------------------------------------------------ GL

    def initializeGL(self) -> None:
        import moderngl
        self.ctx = moderngl.create_context(require=410)
        self.scene = Scene(self.ctx, self.settings)
        info = self.ctx.info
        self.status.emit(f"OpenGL {info['GL_VERSION']} — {info['GL_RENDERER']}")
        self.scene_ready.emit(self.scene)

    def resizeGL(self, w: int, h: int) -> None:
        if self.scene is not None:
            ratio = self.devicePixelRatioF()
            self.scene.resize(int(w * ratio), int(h * ratio))

    def paintGL(self) -> None:
        if self.scene is None or self.ctx is None:
            return
        # See the module docstring: Qt owns the target framebuffer.
        self.ctx.detect_framebuffer(self.defaultFramebufferObject()).use()
        self.scene.render()

    # ------------------------------------------------------------ animation

    def add_animation(self, callback) -> None:
        """Register ``callback(dt_seconds) -> bool``; return False to stop."""
        self._animations.append(callback)

    def clear_animations(self) -> None:
        self._animations.clear()

    def set_spin(self, degrees_per_second: float) -> None:
        self._spin_speed = degrees_per_second

    def _on_tick(self) -> None:
        dt = 1.0 / max(self.settings.target_fps, 1)
        dirty = False
        if self._spin_speed and self.scene is not None:
            self.scene.camera.spin(self._spin_speed * dt)
            dirty = True
        if self._animations:
            still = []
            for cb in self._animations:
                try:
                    if cb(dt) is not False:
                        still.append(cb)
                except Exception as exc:  # keep the UI alive
                    self.status.emit(f"animation stopped: {exc}")
            self._animations = still
            dirty = True
        if dirty:
            self.update()

    # ---------------------------------------------------------------- input

    def mousePressEvent(self, event) -> None:
        self._last_pos = event.position().toPoint()
        self._buttons = event.buttons()

    def mouseReleaseEvent(self, event) -> None:
        pos = event.position().toPoint()
        if (self._last_pos is not None
                and (pos - self._last_pos).manhattanLength() < 3
                and event.button() == Qt.MouseButton.LeftButton):
            self._pick_at(pos)
        self._buttons = Qt.MouseButton.NoButton

    def mouseMoveEvent(self, event) -> None:
        if self._last_pos is None or self.scene is None:
            return
        pos = event.position().toPoint()
        dx = (pos.x() - self._last_pos.x()) / max(self.width(), 1)
        dy = (pos.y() - self._last_pos.y()) / max(self.height(), 1)
        self._last_pos = pos
        buttons = event.buttons()
        mods = event.modifiers()

        if buttons & Qt.MouseButton.LeftButton:
            if mods & Qt.KeyboardModifier.ShiftModifier:
                self.scene.camera.translate(dx, -dy)
            else:
                self.scene.camera.orbit(dx, dy)
            self.update()
        elif buttons & Qt.MouseButton.MiddleButton:
            self.scene.camera.translate(dx, -dy)
            self.update()
        elif buttons & Qt.MouseButton.RightButton:
            self.scene.camera.zoom(1.0 + dy * 2.0)
            self.update()

    def wheelEvent(self, event) -> None:
        if self.scene is None:
            return
        delta = event.angleDelta().y() / 120.0
        self.scene.camera.zoom(0.9 ** delta)
        self.update()

    def keyPressEvent(self, event) -> None:
        if self.scene is None:
            return
        key = event.key()
        cam = self.scene.camera
        if key == Qt.Key.Key_R:
            self.scene.frame_all()
        elif key == Qt.Key.Key_O:
            cam.orthographic = not cam.orthographic
            self.status.emit("orthographic" if cam.orthographic else "perspective")
        elif key == Qt.Key.Key_Space:
            self.set_spin(0.0 if self._spin_speed else 30.0)
        elif key in (Qt.Key.Key_Plus, Qt.Key.Key_Equal):
            self.scene.radius_scale *= 1.15
        elif key == Qt.Key.Key_Minus:
            self.scene.radius_scale /= 1.15
        else:
            super().keyPressEvent(event)
            return
        self.update()

    # -------------------------------------------------------------- picking

    def set_pick_source(self, coords: np.ndarray | None) -> None:
        """Supply the coordinate array that clicks should be tested against."""
        self._pick_source = None if coords is None else np.asarray(coords, np.float64)

    def _pick_at(self, pos: QPoint) -> None:
        if self.scene is None or self._pick_source is None or not len(self._pick_source):
            return
        x_ndc = 2.0 * pos.x() / max(self.width(), 1) - 1.0
        y_ndc = 1.0 - 2.0 * pos.y() / max(self.height(), 1)
        origin, direction = self.scene.camera.ray_through_pixel(x_ndc, y_ndc)

        rel = self._pick_source - origin
        along = rel @ direction
        # Only consider atoms in front of the near plane.
        valid = along > 0
        if not valid.any():
            self.atom_picked.emit(-1)
            return
        perp = np.linalg.norm(rel - np.outer(along, direction), axis=1)
        perp[~valid] = np.inf
        # Among atoms the ray passes close to, take the nearest to the camera.
        hits = np.flatnonzero(perp < 2.2)
        if len(hits) == 0:
            self.atom_picked.emit(-1)
            return
        best = hits[np.argmin(along[hits])]
        self.atom_picked.emit(int(best))
