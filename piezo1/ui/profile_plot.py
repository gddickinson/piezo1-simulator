"""A small XY profile plot, drawn with QPainter.

The pore profile has to be read *against* the hydrophobicity profile — that is
the whole point of Round 19 — so the widget supports two independent y-axes and
threshold lines, and reports where the user clicked so the 3-D view can follow.

Deliberately not matplotlib or pyqtgraph: this needs one plot type, has to
repaint inside a dock at interactive rates, and must match the application's
dark theme. A charting dependency would be more code to configure than to
write, and both drag in their own event-loop integration.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QPainter, QPainterPath, QPen
from PyQt6.QtWidgets import QSizePolicy, QWidget

__all__ = ["ProfilePlot", "Trace", "Marker"]


@dataclass
class Trace:
    """One curve. ``axis`` selects the left (0) or right (1) y-scale."""

    name: str
    x: np.ndarray
    y: np.ndarray
    color: str = "#6fb1ff"
    axis: int = 0
    fill: bool = False
    width: float = 1.8


@dataclass
class Marker:
    """A labelled horizontal or vertical reference line."""

    value: float
    label: str = ""
    color: str = "#8a919e"
    vertical: bool = False
    axis: int = 0
    dashed: bool = True


@dataclass
class _Axis:
    lo: float = 0.0
    hi: float = 1.0
    label: str = ""
    color: str = "#c8ccd4"

    def span(self) -> float:
        return max(self.hi - self.lo, 1e-9)


class ProfilePlot(QWidget):
    """Two-axis line plot with click-to-locate."""

    position_clicked = pyqtSignal(float)

    MARGIN_L = 52
    MARGIN_R = 52
    MARGIN_T = 14
    MARGIN_B = 34

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.traces: list[Trace] = []
        self.markers: list[Marker] = []
        self.x_label = ""
        self.left = _Axis(label="", color="#6fb1ff")
        self.right = _Axis(label="", color="#f2a65a")
        self._x = _Axis()
        self._cursor: float | None = None
        self.setMinimumHeight(190)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)

    # ------------------------------------------------------------------ data

    def clear(self) -> None:
        self.traces, self.markers = [], []
        self._cursor = None
        self.update()

    def set_data(self, traces: list[Trace], markers: list[Marker] | None = None,
                 x_label: str = "", left_label: str = "",
                 right_label: str = "") -> None:
        self.traces = [t for t in traces if len(t.x) and len(t.y)]
        self.markers = list(markers or [])
        self.x_label = x_label
        self.left.label = left_label
        self.right.label = right_label
        self._autoscale()
        self.update()

    def set_cursor(self, x: float | None) -> None:
        self._cursor = x
        self.update()

    def _autoscale(self) -> None:
        if not self.traces:
            return
        xs = np.concatenate([t.x for t in self.traces])
        self._x.lo, self._x.hi = float(np.min(xs)), float(np.max(xs))
        for index, axis in ((0, self.left), (1, self.right)):
            ys = [t.y[np.isfinite(t.y)] for t in self.traces if t.axis == index]
            ys = [y for y in ys if len(y)]
            if not ys:
                axis.lo, axis.hi = 0.0, 1.0
                continue
            joined = np.concatenate(ys)
            lo, hi = float(joined.min()), float(joined.max())
            pad = 0.08 * max(hi - lo, 1e-6)
            # Only pull a left axis down to zero — a radius of 0 is meaningful,
            # whereas hydrophobicity is signed and anchoring it at 0 would
            # squash the part of the range that carries the signal.
            axis.lo = min(0.0, lo - pad) if index == 0 else lo - pad
            axis.hi = hi + pad

    # ----------------------------------------------------------------- paint

    def _plot_rect(self) -> QRectF:
        return QRectF(self.MARGIN_L, self.MARGIN_T,
                      max(self.width() - self.MARGIN_L - self.MARGIN_R, 10),
                      max(self.height() - self.MARGIN_T - self.MARGIN_B, 10))

    def _to_px(self, rect: QRectF, x: float, y: float, axis: int) -> QPointF:
        ax = self.left if axis == 0 else self.right
        fx = (x - self._x.lo) / self._x.span()
        fy = (y - ax.lo) / ax.span()
        return QPointF(rect.left() + fx * rect.width(),
                       rect.bottom() - fy * rect.height())

    def paintEvent(self, event) -> None:      # noqa: N802 (Qt naming)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#181b21"))
        rect = self._plot_rect()

        if not self.traces:
            painter.setPen(QPen(QColor("#6f7684")))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "no profile computed")
            painter.end()
            return

        self._draw_grid(painter, rect)
        for marker in self.markers:
            self._draw_marker(painter, rect, marker)
        for trace in self.traces:
            self._draw_trace(painter, rect, trace)
        self._draw_cursor(painter, rect)
        self._draw_legend(painter, rect)
        painter.end()

    def _draw_grid(self, painter: QPainter, rect: QRectF) -> None:
        font = QFont(); font.setPointSize(8); painter.setFont(font)
        painter.setPen(QPen(QColor("#2a2f38"), 1))
        painter.drawRect(rect)
        for i in range(1, 5):
            y = rect.top() + i * rect.height() / 5.0
            painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
        for i in range(1, 6):
            x = rect.left() + i * rect.width() / 6.0
            painter.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))

        painter.setPen(QPen(QColor("#8a919e")))
        for i in range(6):
            x = self._x.lo + i * self._x.span() / 5.0
            px = rect.left() + i * rect.width() / 5.0
            painter.drawText(QRectF(px - 26, rect.bottom() + 3, 52, 14),
                             Qt.AlignmentFlag.AlignCenter, f"{x:.0f}")
        if self.x_label:
            painter.drawText(QRectF(rect.left(), rect.bottom() + 16,
                                    rect.width(), 14),
                             Qt.AlignmentFlag.AlignCenter, self.x_label)

        for index, axis in ((0, self.left), (1, self.right)):
            if not any(t.axis == index for t in self.traces):
                continue
            painter.setPen(QPen(QColor(axis.color)))
            for i in range(6):
                value = axis.lo + i * axis.span() / 5.0
                py = rect.bottom() - i * rect.height() / 5.0
                box = (QRectF(2, py - 7, self.MARGIN_L - 6, 14) if index == 0
                       else QRectF(rect.right() + 4, py - 7, self.MARGIN_R - 6, 14))
                align = (Qt.AlignmentFlag.AlignRight if index == 0
                         else Qt.AlignmentFlag.AlignLeft)
                painter.drawText(box, align | Qt.AlignmentFlag.AlignVCenter,
                                 f"{value:.2f}")
            if axis.label:
                painter.save()
                painter.translate(12 if index == 0 else self.width() - 8,
                                  rect.center().y())
                painter.rotate(-90)
                painter.drawText(QRectF(-rect.height() / 2, -10,
                                        rect.height(), 14),
                                 Qt.AlignmentFlag.AlignCenter, axis.label)
                painter.restore()

    def _draw_trace(self, painter: QPainter, rect: QRectF, trace: Trace) -> None:
        ok = np.isfinite(trace.y)
        if not ok.any():
            return
        path = QPainterPath()
        started = False
        for xv, yv, good in zip(trace.x, trace.y, ok):
            if not good:
                started = False           # break the line rather than bridge a gap
                continue
            point = self._to_px(rect, float(xv), float(yv), trace.axis)
            if started:
                path.lineTo(point)
            else:
                path.moveTo(point)
                started = True
        painter.setPen(QPen(QColor(trace.color), trace.width))
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawPath(path)

        if trace.fill:
            filled = QPainterPath(path)
            axis = self.left if trace.axis == 0 else self.right
            base = self._to_px(rect, float(trace.x[ok][-1]), axis.lo, trace.axis)
            filled.lineTo(base)
            filled.lineTo(self._to_px(rect, float(trace.x[ok][0]), axis.lo,
                                      trace.axis))
            filled.closeSubpath()
            colour = QColor(trace.color); colour.setAlpha(38)
            painter.fillPath(filled, colour)

    def _draw_marker(self, painter: QPainter, rect: QRectF,
                     marker: Marker) -> None:
        pen = QPen(QColor(marker.color), 1.2)
        if marker.dashed:
            pen.setStyle(Qt.PenStyle.DashLine)
        painter.setPen(pen)
        if marker.vertical:
            x = self._to_px(rect, marker.value, self.left.lo, 0).x()
            painter.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))
            box = QRectF(x + 3, rect.top() + 2, 90, 13)
        else:
            y = self._to_px(rect, self._x.lo, marker.value, marker.axis).y()
            painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
            box = QRectF(rect.left() + 4, y - 14, 150, 13)
        if marker.label:
            font = QFont(); font.setPointSize(7); painter.setFont(font)
            painter.drawText(box, Qt.AlignmentFlag.AlignLeft
                             | Qt.AlignmentFlag.AlignVCenter, marker.label)

    def _draw_cursor(self, painter: QPainter, rect: QRectF) -> None:
        if self._cursor is None:
            return
        x = self._to_px(rect, self._cursor, self.left.lo, 0).x()
        if rect.left() <= x <= rect.right():
            painter.setPen(QPen(QColor("#f0f3f8"), 1.0))
            painter.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))

    def _draw_legend(self, painter: QPainter, rect: QRectF) -> None:
        font = QFont(); font.setPointSize(8); painter.setFont(font)
        y = rect.top() + 4
        for trace in self.traces:
            painter.setPen(QPen(QColor(trace.color), 2.4))
            painter.drawLine(QPointF(rect.right() - 116, y + 6),
                             QPointF(rect.right() - 100, y + 6))
            painter.setPen(QPen(QColor("#c8ccd4")))
            painter.drawText(QRectF(rect.right() - 95, y, 92, 13),
                             Qt.AlignmentFlag.AlignLeft, trace.name)
            y += 14

    # ----------------------------------------------------------------- input

    def mousePressEvent(self, event) -> None:      # noqa: N802
        rect = self._plot_rect()
        if not self.traces or not rect.contains(event.position()):
            return
        frac = (event.position().x() - rect.left()) / max(rect.width(), 1e-9)
        value = self._x.lo + frac * self._x.span()
        self.set_cursor(value)
        self.position_clicked.emit(float(value))
