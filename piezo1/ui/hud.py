"""Heads-up display drawn over the viewport: scale bar, clock, readouts.

Three things a figure or a recording needs and a bare 3-D view does not give:

* a **scale bar**, so a rendered image carries its own units. A structure
  screenshot without one states nothing quantitative, and readers cannot
  reconstruct the scale from a perspective projection by eye;
* a **time counter** for animations, since a morph or a mode sweep is
  meaningless without knowing where in the cycle a frame sits — and, for a
  morph, whether the axis is time at all (it is not; see below);
* **readouts** of whatever has been measured, so the number and the thing it
  describes appear in the same frame.

Drawn with QPainter on a transparent child widget rather than in GL. The
renderer owns depth and the ray-cast impostors write ``gl_FragDepth``; mixing
2-D overlay drawing into that pipeline fought the depth buffer when it was tried
for text labels in an earlier round, and lost silently.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPainter, QPen
from PyQt6.QtWidgets import QWidget

__all__ = ["HudOverlay", "HudSettings", "nice_scale_length"]

#: Lengths a scale bar is allowed to take, in Angstrom. Chosen so the label is
#: always a round number a reader can hold in their head.
NICE_LENGTHS = [1, 2, 5, 10, 20, 50, 100, 200, 500, 1000]


def nice_scale_length(target: float) -> float:
    """Largest tidy length not exceeding ``target``, in Angstrom."""
    usable = [v for v in NICE_LENGTHS if v <= target]
    return float(usable[-1]) if usable else float(NICE_LENGTHS[0])


@dataclass
class HudSettings:
    """What the overlay shows. Every element is independently switchable."""

    scale_bar: bool = True
    clock: bool = False
    readouts: bool = True
    structure_name: bool = True
    orientation_axes: bool = False
    font_scale: float = 1.0
    corner: str = "bottom-left"          # where the scale bar sits

    #: Named readouts the user has chosen to display, in order.
    fields: list[str] = field(default_factory=list)

    def as_dict(self) -> dict:
        return {"scale_bar": self.scale_bar, "clock": self.clock,
                "readouts": self.readouts, "structure_name": self.structure_name,
                "orientation_axes": self.orientation_axes,
                "font_scale": self.font_scale, "corner": self.corner,
                "fields": list(self.fields)}

    @classmethod
    def from_dict(cls, data: dict) -> "HudSettings":
        known = {k: v for k, v in (data or {}).items()
                 if k in cls.__dataclass_fields__}
        return cls(**known)


class HudOverlay(QWidget):
    """Transparent overlay: scale bar, animation clock and measured values."""

    MARGIN = 16

    def __init__(self, viewport) -> None:
        super().__init__(viewport)
        self.viewport = viewport
        self.settings = HudSettings()
        self.readouts: dict[str, str] = {}
        self.structure_name = ""
        #: Set when part of what is drawn is prediction rather than experiment.
        #: Deliberately not a member of `HudSettings`: it is not a preference.
        self.provenance = ""
        self.clock_text = ""
        self.clock_note = ""
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents, True)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground, True)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self.setGeometry(viewport.rect())

    # ------------------------------------------------------------------ data

    def set_readout(self, key: str, text: str) -> None:
        """Record a measured value. Empty text removes it."""
        if text:
            self.readouts[key] = text
        else:
            self.readouts.pop(key, None)
        self.update()

    def clear_readouts(self) -> None:
        self.readouts.clear()
        self.update()

    def set_clock(self, text: str, note: str = "") -> None:
        self.clock_text, self.clock_note = text, note
        self.update()

    # ---------------------------------------------------------------- scale

    def world_per_pixel(self) -> float:
        """Angstrom covered by one logical pixel at the camera's pivot depth.

        Exact only in the plane through the pivot, which is where the molecule
        is. A perspective scale bar is always a statement about one depth, and
        this one says which.
        """
        scene = getattr(self.viewport, "scene", None)
        if scene is None or self.height() <= 0:
            return 0.0
        camera = scene.camera
        if getattr(camera, "orthographic", False):
            visible_height = 2.0 * camera.distance * np.tan(
                np.radians(camera.fov) / 2.0)
        else:
            visible_height = 2.0 * camera.distance * np.tan(
                np.radians(camera.fov) / 2.0)
        return float(visible_height / self.height())

    # ---------------------------------------------------------------- paint

    def paintEvent(self, event) -> None:            # noqa: N802 (Qt naming)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        base = max(8, int(11 * self.settings.font_scale))
        font = QFont()
        font.setPointSize(base)
        painter.setFont(font)

        if self.settings.structure_name and self.structure_name:
            self._draw_title(painter, base)
        if self.settings.readouts and self.readouts:
            self._draw_readouts(painter, base)
        if self.settings.clock and self.clock_text:
            self._draw_clock(painter, base)
        if self.settings.scale_bar:
            self._draw_scale_bar(painter, base)
        if self.settings.orientation_axes:
            self._draw_axes(painter, base)
        if self.provenance:
            self._draw_provenance(painter, base)
        painter.end()

    def _shadowed(self, painter: QPainter, rect: QRectF, flags, text: str,
                  colour: str = "#f0f3f8") -> None:
        """Text with a dark offset copy, so it stays legible on any background."""
        painter.setPen(QPen(QColor(0, 0, 0, 190)))
        painter.drawText(rect.translated(1.0, 1.0), flags, text)
        painter.setPen(QPen(QColor(colour)))
        painter.drawText(rect, flags, text)

    def _draw_provenance(self, painter: QPainter, base: int) -> None:
        """An amber banner whenever part of what is drawn is a prediction.

        Not switchable, and not in the display-options dialog with the other
        readouts. Every other element of this overlay is a convenience; this
        one is the difference between a measurement and a model, and a user who
        has turned the title off and come back to the window an hour later has
        no other way to tell. It is the same reasoning as the HaloTag fold's
        status line, which also cannot be suppressed.
        """
        font = painter.font()
        font.setBold(True)
        painter.setFont(font)
        rect = QRectF(self.MARGIN, self.height() - self.MARGIN - base * 2.4,
                      self.width() - 2 * self.MARGIN, base * 2.0)
        self._shadowed(painter, rect,
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom,
                       self.provenance, colour="#ffb454")
        font.setBold(False)
        painter.setFont(font)

    def _draw_title(self, painter: QPainter, base: int) -> None:
        font = painter.font()
        font.setPointSize(int(base * 1.35))
        font.setBold(True)
        painter.setFont(font)
        self._shadowed(painter,
                       QRectF(self.MARGIN, self.MARGIN, self.width() - 2 * self.MARGIN,
                              base * 2.2),
                       Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
                       self.structure_name)
        font.setPointSize(base)
        font.setBold(False)
        painter.setFont(font)

    def _draw_readouts(self, painter: QPainter, base: int) -> None:
        line = base * 1.7
        top = self.MARGIN + (base * 2.6 if self.settings.structure_name
                             and self.structure_name else 0)
        keys = self.settings.fields or list(self.readouts)
        for i, key in enumerate([k for k in keys if k in self.readouts]):
            self._shadowed(
                painter,
                QRectF(self.MARGIN, top + i * line, self.width() - 2 * self.MARGIN,
                       line),
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
                self.readouts[key], "#c8ccd4")

    def _draw_clock(self, painter: QPainter, base: int) -> None:
        font = painter.font()
        font.setPointSize(int(base * 1.25))
        painter.setFont(font)
        box = QRectF(self.width() - 260 - self.MARGIN, self.MARGIN, 260,
                     base * 2.0)
        self._shadowed(painter, box,
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTop,
                       self.clock_text, "#f2a65a")
        if self.clock_note:
            font.setPointSize(max(7, int(base * 0.85)))
            painter.setFont(font)
            self._shadowed(painter, box.translated(0, base * 2.0),
                           Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTop,
                           self.clock_note, "#8a919e")
        font.setPointSize(base)
        painter.setFont(font)

    def _draw_scale_bar(self, painter: QPainter, base: int) -> None:
        scale = self.world_per_pixel()
        if scale <= 0:
            return
        # Aim for about a fifth of the viewport, then round down to a tidy
        # number so the label is readable rather than exact-but-ugly.
        target = 0.2 * self.width() * scale
        length_a = nice_scale_length(target)
        pixels = length_a / scale
        if pixels < 24 or pixels > self.width() * 0.8:
            return

        y = self.height() - self.MARGIN - base * 2.2
        x = (self.MARGIN if "left" in self.settings.corner
             else self.width() - self.MARGIN - pixels)

        painter.setPen(QPen(QColor(0, 0, 0, 190), 5.0, Qt.PenStyle.SolidLine,
                            Qt.PenCapStyle.FlatCap))
        painter.drawLine(QPointF(x, y), QPointF(x + pixels, y))
        painter.setPen(QPen(QColor("#f0f3f8"), 3.0, Qt.PenStyle.SolidLine,
                            Qt.PenCapStyle.FlatCap))
        painter.drawLine(QPointF(x, y), QPointF(x + pixels, y))
        for end in (x, x + pixels):
            painter.drawLine(QPointF(end, y - 5), QPointF(end, y + 5))

        label = (f"{length_a:.0f} Å" if length_a < 10
                 else f"{length_a / 10:.0f} nm ({length_a:.0f} Å)")
        self._shadowed(painter, QRectF(x, y + 4, max(pixels, 120), base * 1.9),
                       Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
                       label)

    def _draw_axes(self, painter: QPainter, base: int) -> None:
        """A small orientation gnomon, rotated with the camera."""
        scene = getattr(self.viewport, "scene", None)
        if scene is None:
            return
        from ..render.camera import quat_to_matrix
        rotation = quat_to_matrix(scene.camera.rotation)
        size = 30.0
        origin = QPointF(self.width() - self.MARGIN - size - 14,
                         self.height() - self.MARGIN - size - 14)
        for axis, colour, name in ((0, "#f26d6d", "x"), (1, "#7ed67e", "y"),
                                   (2, "#6fb1ff", "z")):
            direction = np.zeros(3)
            direction[axis] = 1.0
            local = rotation @ direction
            tip = QPointF(origin.x() + float(local[0]) * size,
                          origin.y() - float(local[1]) * size)
            painter.setPen(QPen(QColor(colour), 2.0))
            painter.drawLine(origin, tip)
            self._shadowed(painter, QRectF(tip.x() - 7, tip.y() - 9, 16, 14),
                           Qt.AlignmentFlag.AlignCenter, name, colour)
