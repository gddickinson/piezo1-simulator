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
from PyQt6.QtGui import QColor, QFont, QPainter, QPen, QSurfaceFormat
from PyQt6.QtOpenGLWidgets import QOpenGLWidget
from PyQt6.QtWidgets import QWidget

from ..config import RenderSettings
from ..render.scene import Scene
from .hud import HudOverlay

__all__ = ["ViewportWidget", "configure_surface_format", "CLICK_SLOP",
           "nearest_hit", "PRIMARY_SOURCE"]

#: How far the mouse may move between press and release and still count as a
#: click rather than a drag, in pixels. A few pixels of tremor is normal on a
#: trackpad; a rotation is tens to hundreds.
CLICK_SLOP = 3

#: How close, in Angstrom, the pick ray must pass to an atom centre to hit it.
#: Generous relative to a carbon so that a cartoon ribbon — which draws no atom
#: where the ray lands — still picks the residue the user aimed at.
PICK_RADIUS = 2.2

#: The name the primary structure's atoms are picked under. Every other pick
#: source is a drawn feature — a modelled tag, an extra structure — registered
#: by whichever controller drew it, so that clicking anything on screen
#: identifies it rather than silently identifying whatever primary atom
#: happens to lie behind it.
PRIMARY_SOURCE = "structure"


def nearest_hit(sources: dict, origin: np.ndarray, direction: np.ndarray,
                radius: float = PICK_RADIUS, masks: dict | None = None):
    """The atom nearest the camera along a pick ray, across named sources.

    ``sources`` maps a name to an ``(n, 3)`` coordinate array. ``masks``, when
    given, maps a source name to a boolean per-atom array; atoms marked False
    cannot be hit — they are the ones the user has hidden, and a click must
    answer for what is *on screen*, not for what the arrays hold. Returns
    ``(name, index)`` for the winning atom, or ``None`` when the ray passes
    nothing. Pure geometry, deliberately free of Qt and GL, so the rule that
    decides *which* thing a click hits — nearest visible wins, whatever drew
    it — can be tested on arrays whose answer is known by construction.
    """
    best = None                                   # (along, name, index)
    for name, coords in sources.items():
        if coords is None or len(coords) == 0:
            continue
        rel = np.asarray(coords, np.float64) - origin
        along = rel @ direction
        valid = along > 0                          # in front of the camera
        mask = None if masks is None else masks.get(name)
        if mask is not None:
            valid = valid & np.asarray(mask, dtype=bool)
        if not valid.any():
            continue
        perp = np.linalg.norm(rel - np.outer(along, direction), axis=1)
        perp[~valid] = np.inf
        hits = np.flatnonzero(perp < radius)
        if len(hits) == 0:
            continue
        i = int(hits[np.argmin(along[hits])])
        if best is None or along[i] < best[0]:
            best = (float(along[i]), name, i)
    return None if best is None else (best[1], best[2])


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
    #: Emitted with a feature name and atom index when the click lands on a
    #: registered feature — a modelled tag, an extra structure — instead of
    #: the primary model.
    feature_picked = pyqtSignal(str, int)
    #: Emitted with a human-readable status string.
    status = pyqtSignal(str)
    #: Right-click: the position to pop up at, and the atom under it (-1 none).
    context_requested = pyqtSignal(QPoint, int)

    def __init__(self, settings: RenderSettings | None = None, parent=None) -> None:
        super().__init__(parent)
        self.settings = settings or RenderSettings()
        self.ctx = None
        self.scene: Scene | None = None
        self._last_pos: QPoint | None = None
        self._press_pos: QPoint | None = None
        self._buttons = Qt.MouseButton.NoButton
        self._pick_source: np.ndarray | None = None
        #: Which primary atoms are actually drawn. The pick source is the full
        #: atom array; hiding a category or choosing a component must also
        #: stop the hidden atoms answering clicks, or a click identifies an
        #: invisible lipid in front of the visible helix the user aimed at.
        self._pick_mask: np.ndarray | None = None
        #: Extra pickable coordinate sets, keyed by feature name. Registered
        #: by the controller that draws each feature and dropped when it
        #: clears, so a click can never identify something not on screen.
        self._feature_sources: dict[str, np.ndarray] = {}
        self._spin_speed = 0.0
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        #: World-space labels drawn as a 2D overlay: (xyz, text, colour).
        self.overlay_labels: list = []
        self.measure_mode = False
        # Text is drawn by a transparent child widget rather than by QPainter
        # inside paintGL. Mixing QPainter with a moderngl context in the same
        # paint call is fragile - the two fight over GL state, and the text
        # simply failed to appear. A sibling widget has its own paint event and
        # cannot interfere with the scene at all.
        self.overlay = _LabelOverlay(self)
        #: Scale bar, animation clock and measured readouts. Same reasoning as
        #: the label overlay: a sibling widget, not QPainter inside paintGL.
        self.hud = HudOverlay(self)

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
        # The scale bar depends on camera distance, so it would go stale on
        # every zoom if it were not repainted alongside the scene.
        self._refresh_overlays()

    def set_overlay_labels(self, labels) -> None:
        self.overlay_labels = list(labels)
        self.overlay.update()

    def resizeEvent(self, event) -> None:
        super().resizeEvent(event)
        self.overlay.setGeometry(self.rect())
        self.hud.setGeometry(self.rect())

    # ------------------------------------------------------------ animation

    def add_animation(self, callback) -> None:
        """Register ``callback(dt_seconds) -> bool``; return False to stop."""
        self._animations.append(callback)

    def clear_animations(self) -> None:
        self._animations.clear()

    def set_spin(self, degrees_per_second: float) -> None:
        self._spin_speed = degrees_per_second

    def _refresh_overlays(self) -> None:
        """Repaint the 2-D overlays. The scale bar depends on camera distance,
        so it goes stale on every zoom unless it is repainted with the scene."""
        self.overlay.update()
        self.hud.update()

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
        self._press_pos = self._last_pos
        self._buttons = event.buttons()

    def mouseReleaseEvent(self, event) -> None:
        """A click picks; a drag rotates and must not.

        The two are told apart by how far the mouse moved *since the button
        went down*, which is why the press position is kept separately.
        Comparing against ``_last_pos`` cannot work: ``mouseMoveEvent`` updates
        it on every step to compute the drag delta, so by release time it is
        the release position and the distance is always zero — every rotation
        ended in a pick. That was merely noisy while a pick only rewrote the
        status bar; once a pick highlights, it repaints the selection each time
        the user turns the structure.
        """
        pos = event.position().toPoint()
        still = (self._press_pos is not None
                 and (pos - self._press_pos).manhattanLength() < CLICK_SLOP)
        if still and event.button() == Qt.MouseButton.LeftButton:
            self._pick_at(pos)
        elif still and event.button() == Qt.MouseButton.RightButton:
            # Right-drag zooms, so the menu is on the *click* — the same
            # distinction the left button already makes between picking and
            # rotating. A user who drags to zoom never gets a menu they did
            # not ask for.
            #
            # Resolved across every pick source, not only the primary: the
            # menu's residue entries are annotation read by primary atom
            # index, so when a drawn feature is nearest the click the menu
            # opens generic (-1) rather than naming the occluded residue
            # behind the tag the user actually aimed at.
            hit = self.hit_at(pos)
            index = (hit[1] if hit is not None and hit[0] == PRIMARY_SOURCE
                     else -1)
            self.context_requested.emit(pos, index)
        self._press_pos = None
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
        self._pick_mask = None                     # a new array, a new mask

    def set_pick_mask(self, mask: np.ndarray | None) -> None:
        """Restrict picking to the atoms currently drawn. None means all."""
        if (mask is not None and self._pick_source is not None
                and len(mask) != len(self._pick_source)):
            # A mask for a different array would silently shift every pick;
            # dropping it fails towards "everything pickable", which is the
            # pre-mask behaviour rather than a new wrong one.
            mask = None
        self._pick_mask = None if mask is None else np.asarray(mask, dtype=bool)

    def set_feature_pick_source(self, name: str,
                                coords: np.ndarray | None) -> None:
        """Register (or, with ``None``, drop) a feature's pickable atoms."""
        if coords is None:
            self._feature_sources.pop(name, None)
        else:
            self._feature_sources[name] = np.asarray(coords, np.float64)

    def clear_feature_pick_sources(self) -> None:
        self._feature_sources.clear()

    def hit_at(self, pos: QPoint):
        """What is under ``pos``, across everything pickable on screen.

        Returns ``(source_name, index)`` — ``PRIMARY_SOURCE`` for the loaded
        structure — or ``None`` when nothing is hit. Nearest to the camera
        wins whatever drew it, because "what did I click" has one honest
        answer: the thing in front.
        """
        if self.scene is None:
            return None
        sources: dict[str, np.ndarray] = {}
        if self._pick_source is not None and len(self._pick_source):
            sources[PRIMARY_SOURCE] = self._pick_source
        sources.update(self._feature_sources)
        if not sources:
            return None
        x_ndc = 2.0 * pos.x() / max(self.width(), 1) - 1.0
        y_ndc = 1.0 - 2.0 * pos.y() / max(self.height(), 1)
        origin, direction = self.scene.camera.ray_through_pixel(x_ndc, y_ndc)
        masks = ({PRIMARY_SOURCE: self._pick_mask}
                 if self._pick_mask is not None else None)
        return nearest_hit(sources, origin, direction, masks=masks)

    def _pick_at(self, pos: QPoint) -> None:
        """Identify the thing under ``pos`` and announce it as a selection."""
        hit = self.hit_at(pos)
        if hit is not None:
            name, index = hit
            if name == PRIMARY_SOURCE:
                self.atom_picked.emit(index)
            else:
                self.feature_picked.emit(name, index)
        elif ((self._pick_source is not None and len(self._pick_source))
              or self._feature_sources):
            # Something is pickable and the click missed it all; the status
            # line should say so rather than leave the click unanswered.
            self.atom_picked.emit(-1)

    def atom_at(self, pos: QPoint) -> int | None:
        """The atom under ``pos``, or -1 for none — without announcing it.

        Separate from ``_pick_at`` because the right-click menu needs to know
        what it was opened on without that counting as a selection: opening a
        menu and then dismissing it must leave the model exactly as it was.
        ``None`` means there is nothing to pick against at all.
        """
        if self.scene is None or self._pick_source is None or not len(self._pick_source):
            return None
        x_ndc = 2.0 * pos.x() / max(self.width(), 1) - 1.0
        y_ndc = 1.0 - 2.0 * pos.y() / max(self.height(), 1)
        origin, direction = self.scene.camera.ray_through_pixel(x_ndc, y_ndc)
        masks = ({PRIMARY_SOURCE: self._pick_mask}
                 if self._pick_mask is not None else None)
        hit = nearest_hit({PRIMARY_SOURCE: self._pick_source},
                          origin, direction, masks=masks)
        return -1 if hit is None else hit[1]


class _LabelOverlay(QWidget):
    """Transparent child widget that draws world-anchored text labels."""

    def __init__(self, viewport: "ViewportWidget") -> None:
        super().__init__(viewport)
        self.viewport = viewport
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self.setGeometry(viewport.rect())

    def project(self, xyz) -> tuple[float, float] | None:
        """World point to widget pixels, or None if it is not on screen."""
        scene = self.viewport.scene
        if scene is None:
            return None
        view, proj = scene.camera.matrices()
        clip = (proj @ view) @ np.append(np.asarray(xyz, dtype=float), 1.0)
        if clip[3] <= 0:
            return None                     # behind the camera
        ndc = clip[:3] / clip[3]
        if not (-1.15 < ndc[0] < 1.15 and -1.15 < ndc[1] < 1.15):
            return None
        return ((ndc[0] * 0.5 + 0.5) * self.width(),
                (1.0 - (ndc[1] * 0.5 + 0.5)) * self.height())

    def paintEvent(self, event) -> None:
        labels = self.viewport.overlay_labels
        if not labels:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setFont(QFont("Helvetica", 11))
        for xyz, text, colour in labels:
            pos = self.project(xyz)
            if pos is None:
                continue
            x, y = int(pos[0]), int(pos[1])
            # A dark halo, so labels stay readable over pale cartoon ribbons.
            painter.setPen(QPen(QColor(8, 10, 16, 220), 3))
            for dx, dy in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                painter.drawText(x + dx, y + dy, text)
            painter.setPen(QPen(QColor(*colour)))
            painter.drawText(x, y, text)
        painter.end()
