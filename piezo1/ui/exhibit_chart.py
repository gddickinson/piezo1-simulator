"""The chart the explore window paints, in the idiom of ``profile_plot``.

Deliberately not matplotlib, for the reasons that module gives: this has to
repaint inside a window at interactive rates while a slider moves, and it has
to match the application's own theme. What it adds over ``ProfilePlot`` is what
the exhibits need and the pore profile never did — **bars**, **points**,
**log axes** and **reference bands** — so a categorical result (contacts by
kind, panels by status) and a decade-spanning one (potency, calcium) are both
drawable without inventing a second convention for either.

Two axes, because several exhibits are one quantity read against another in
different units: calcium against sensor occupancy, bottleneck radius against
wetting score. ``log_y`` applies to the **left** axis only; the right one is
always linear, and saying so here is cheaper than a second flag nobody sets.
"""

from __future__ import annotations

import math

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import QColor, QFont, QPainter, QPalette, QPen
from PyQt6.QtWidgets import QSizePolicy, QWidget

from .exhibits import ChartData

__all__ = ["ChartView"]

#: Fallbacks only. Ink is taken from the widget's own palette at paint time,
#: so the chart follows the interface theme the way every other panel does —
#: Round 90 added a light theme, and a chart with dark-theme text baked in is
#: grey on white there.
GRID = "#2a2f38"
TEXT = "#c8ccd4"
MUTED = "#8a919e"


def _finite(values):
    return [v for v in values if v is not None and math.isfinite(v)]


class _Scale:
    """One axis: data range in, pixels out, with an optional log transform."""

    def __init__(self, log: bool = False) -> None:
        self.lo, self.hi, self.log = 0.0, 1.0, log

    def fit(self, values, pad: float = 0.08, include_zero: bool = False) -> None:
        values = _finite(values)
        if self.log:
            values = [v for v in values if v > 0]
        if not values:
            self.lo, self.hi = (0.0, 1.0)
            return
        lo, hi = min(values), max(values)
        if include_zero and not self.log:
            lo, hi = min(lo, 0.0), max(hi, 0.0)
        if self.log:
            lo, hi = math.log10(lo), math.log10(hi)
        if hi - lo < 1e-12:
            lo, hi = lo - 0.5, hi + 0.5
        span = hi - lo
        self.lo, self.hi = lo - pad * span, hi + pad * span

    def to_unit(self, value: float) -> float:
        if value is None or not math.isfinite(value):
            return float("nan")
        if self.log:
            if value <= 0:
                return float("nan")
            value = math.log10(value)
        return (value - self.lo) / max(self.hi - self.lo, 1e-12)

    def ticks(self, count: int = 5) -> list[tuple[float, str]]:
        out = []
        for index in range(count + 1):
            unit = index / count
            raw = self.lo + unit * (self.hi - self.lo)
            value = 10.0 ** raw if self.log else raw
            out.append((unit, _label(value)))
        return out


def _label(value: float) -> str:
    if value == 0:
        return "0"
    magnitude = abs(value)
    if magnitude >= 1e4 or magnitude < 1e-2:
        return f"{value:.1e}"
    if magnitude >= 100:
        return f"{value:.0f}"
    return f"{value:.3g}"


class ChartView(QWidget):
    """Paints one :class:`~piezo1.ui.exhibits.ChartData`."""

    MARGIN_L = 62
    MARGIN_R = 62
    MARGIN_T = 26
    MARGIN_B = 52

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.chart = ChartData()
        self._x = _Scale()
        self._left = _Scale()
        self._right = _Scale()
        self.setMinimumHeight(240)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self._ink_text = QColor(TEXT)
        self._ink_muted = QColor(MUTED)
        self._ink_grid = QColor(GRID)

    def _ink(self) -> None:
        """Text, grid and muted ink from the palette, once per paint."""
        base = self.palette().color(QPalette.ColorRole.WindowText)
        self._ink_text = QColor(base)
        self._ink_muted = QColor(base)
        self._ink_muted.setAlpha(160)
        self._ink_grid = QColor(base)
        self._ink_grid.setAlpha(52)

    # ------------------------------------------------------------------ data

    def set_chart(self, chart: ChartData) -> None:
        self.chart = chart or ChartData()
        self._rescale()
        self.update()

    def _rescale(self) -> None:
        chart = self.chart
        self._x = _Scale(chart.log_x)
        self._left = _Scale(chart.log_y)
        self._right = _Scale(False)
        bars = [s for s in chart.series if s.kind == "bar"]
        xs = [v for s in chart.series for v in s.x]
        if bars and chart.categories:
            # Categorical: one slot per category, with room either side.
            self._x.lo, self._x.hi = -0.6, len(chart.categories) - 0.4
        else:
            vertical = [r.value for r in chart.references if r.vertical]
            self._x.fit(xs + vertical)
        for index, scale in ((0, self._left), (1, self._right)):
            values = [v for s in chart.series if s.axis == index for v in s.y]
            if index == 0:
                for ref in chart.references:
                    if not ref.vertical:
                        values += [ref.value] + ([ref.high] if ref.high else [])
            if not values:
                continue
            # A bar is read from zero, so its axis has to contain zero or the
            # picture exaggerates every difference on it.
            scale.fit(values, include_zero=any(s.kind == "bar" and
                                               s.axis == index
                                               for s in chart.series))

    # --------------------------------------------------------------- paint

    def _rect(self) -> QRectF:
        return QRectF(self.MARGIN_L, self.MARGIN_T,
                      max(self.width() - self.MARGIN_L - self.MARGIN_R, 10.0),
                      max(self.height() - self.MARGIN_T - self.MARGIN_B, 10.0))

    def _point(self, rect: QRectF, x: float, y: float, axis: int) -> QPointF:
        scale = self._right if axis == 1 else self._left
        return QPointF(rect.left() + self._x.to_unit(x) * rect.width(),
                       rect.bottom() - scale.to_unit(y) * rect.height())

    def paintEvent(self, event) -> None:      # noqa: N802 (Qt naming)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        font = QFont()
        font.setPointSize(9)
        painter.setFont(font)
        self._ink()
        rect = self._rect()

        if self.chart.empty:
            painter.setPen(QPen(self._ink_muted))
            painter.drawText(self.rect().adjusted(24, 24, -24, -24),
                             int(Qt.AlignmentFlag.AlignCenter
                                 | Qt.TextFlag.TextWordWrap),
                             self.chart.note or "nothing to draw")
            painter.end()
            return

        self._grid(painter, rect)
        for ref in self.chart.references:
            self._band(painter, rect, ref)
        for series in self.chart.series:
            if series.kind == "bar":
                self._bars(painter, rect, series)
            elif series.kind == "point":
                self._points(painter, rect, series)
            else:
                self._line(painter, rect, series)
        # Lines and labels go on *top* of the data: a reference band behind a
        # single wide bar is invisible, and the whole point of the band is
        # that the bar has to be read against it.
        for ref in self.chart.references:
            self._reference(painter, rect, ref)
        self._legend(painter, rect)
        painter.end()

    def _grid(self, painter: QPainter, rect: QRectF) -> None:
        painter.setPen(QPen(self._ink_grid, 1))
        painter.drawRect(rect)
        painter.setPen(QPen(self._ink_text))

        for unit, text in self._left.ticks():
            y = rect.bottom() - unit * rect.height()
            painter.setPen(QPen(self._ink_grid, 1, Qt.PenStyle.DotLine))
            painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
            painter.setPen(QPen(self._ink_text))
            painter.drawText(QRectF(0, y - 8, self.MARGIN_L - 6, 16),
                             int(Qt.AlignmentFlag.AlignRight
                                 | Qt.AlignmentFlag.AlignVCenter), text)

        if any(s.axis == 1 for s in self.chart.series):
            painter.setPen(QPen(self._ink_muted))
            for unit, text in self._right.ticks():
                y = rect.bottom() - unit * rect.height()
                painter.drawText(QRectF(rect.right() + 6, y - 8,
                                        self.MARGIN_R - 8, 16),
                                 int(Qt.AlignmentFlag.AlignLeft
                                     | Qt.AlignmentFlag.AlignVCenter), text)

        painter.setPen(QPen(self._ink_text))
        if self.chart.categories:
            # Each label is confined to its own slot and wraps inside it.
            # Given a fixed width they overlapped, and two overlapping domain
            # names read as one domain nobody has heard of.
            slot = rect.width() / max(len(self.chart.categories), 1)
            for index, name in enumerate(self.chart.categories):
                centre = rect.left() + self._x.to_unit(index) * rect.width()
                painter.drawText(QRectF(centre - slot / 2.0, rect.bottom() + 4,
                                        slot, self.MARGIN_B - 6),
                                 int(Qt.AlignmentFlag.AlignHCenter
                                     | Qt.AlignmentFlag.AlignTop
                                     | Qt.TextFlag.TextWordWrap), name)
        else:
            for unit, text in self._x.ticks(4):
                x = rect.left() + unit * rect.width()
                painter.drawText(QRectF(x - 45, rect.bottom() + 4, 90, 16),
                                 int(Qt.AlignmentFlag.AlignHCenter), text)
            painter.setPen(QPen(self._ink_muted))
            painter.drawText(QRectF(rect.left(), rect.bottom() + 22,
                                    rect.width(), 18),
                             int(Qt.AlignmentFlag.AlignHCenter),
                             self.chart.x_label)
        painter.setPen(QPen(self._ink_muted))
        painter.drawText(QRectF(rect.left(), 4, rect.width(), 18),
                         int(Qt.AlignmentFlag.AlignLeft), self.chart.y_label)

    def _ref_label(self, painter: QPainter, box: QRectF, text: str,
                   align) -> None:
        """A reference label on a plate of the widget's own background.

        Bars and curves move with the data, so there is no corner that is
        reliably free: "no enrichment" landed on empty space under one
        partition and on top of a bar under the other. The plate keeps the
        label legible over whatever is behind it without hiding the data.
        """
        if not text:
            return
        metrics = painter.fontMetrics()
        width = min(metrics.horizontalAdvance(text) + 10, box.width())
        left = (box.right() - width if align == Qt.AlignmentFlag.AlignRight
                else box.left())
        plate = QColor(self.palette().color(QPalette.ColorRole.Window))
        plate.setAlpha(205)
        painter.fillRect(QRectF(left, box.top(), width, box.height()), plate)
        pen = painter.pen()
        painter.drawText(QRectF(left + 4, box.top(), width - 8, box.height()),
                         int(align | Qt.AlignmentFlag.AlignVCenter), text)
        painter.setPen(pen)

    def _band(self, painter: QPainter, rect: QRectF, ref) -> None:
        """The shaded part of a reference range, drawn before the data."""
        if ref.vertical or ref.high is None:
            return
        top = self._point(rect, self._x.lo, max(ref.value, ref.high), 0).y()
        bottom = self._point(rect, self._x.lo, min(ref.value, ref.high), 0).y()
        if not (math.isfinite(top) and math.isfinite(bottom)):
            return
        colour = QColor(ref.color)
        colour.setAlpha(46)
        painter.fillRect(QRectF(rect.left(), top, rect.width(),
                                max(bottom - top, 1.0)), colour)

    def _reference(self, painter: QPainter, rect: QRectF, ref) -> None:
        pen = QPen(QColor(ref.color), 1.2, Qt.PenStyle.DashLine)
        painter.setPen(pen)
        if ref.vertical:
            unit = self._x.to_unit(ref.value)
            if not math.isfinite(unit):
                return
            x = rect.left() + unit * rect.width()
            painter.drawLine(QPointF(x, rect.top()), QPointF(x, rect.bottom()))
            # Labelled along the bottom rather than the top: the legend is in
            # the top-right corner, and two vertical markers near it wrote
            # over each other and over the legend.
            width = max(rect.right() - x - 8, 60.0)
            self._ref_label(painter, QRectF(x + 4, rect.bottom() - 19, width,
                                            17), ref.label,
                            Qt.AlignmentFlag.AlignLeft)
            return
        if ref.high is not None:
            edge = self._left.to_unit(ref.high)
            if math.isfinite(edge):
                y = rect.bottom() - edge * rect.height()
                painter.drawLine(QPointF(rect.left(), y),
                                 QPointF(rect.right(), y))
        unit = self._left.to_unit(ref.value)
        if not math.isfinite(unit):
            return
        y = rect.bottom() - unit * rect.height()
        painter.drawLine(QPointF(rect.left(), y), QPointF(rect.right(), y))
        # Right-aligned: bars grow from the left, so a label there lands on
        # top of the first one and is read against the wrong colour.
        self._ref_label(painter, QRectF(rect.right() - 306, y - 17, 300, 16),
                        ref.label, Qt.AlignmentFlag.AlignRight)

    # ------------------------------------------------------------- series

    def _bar_slots(self):
        bars = [s for s in self.chart.series if s.kind == "bar"]
        return bars, max(len(bars), 1)

    def _bars(self, painter: QPainter, rect: QRectF, series) -> None:
        bars, count = self._bar_slots()
        # By identity: `Series` is a dataclass, so two series with the same
        # numbers compare equal and `index` would put both in one slot.
        index = next(i for i, s in enumerate(bars) if s is series)
        slot = rect.width() / max(len(self.chart.categories) or len(series.x), 1)
        # Capped: with one category the slot is the whole panel, and a bar
        # that wide reads as a filled background rather than a value.
        width = min(slot * 0.72 / count, 96.0)
        scale = self._right if series.axis == 1 else self._left
        if scale.log:
            # There is no zero on a log axis, so a bar is read from the floor
            # of the drawn range. Computing it from zero gives log(0) and a
            # bar with no geometry at all — which paints nothing and looks
            # exactly like a category with no data.
            zero = rect.bottom()
        else:
            zero = (rect.bottom()
                    - max(min(scale.to_unit(0.0), 1.0), 0.0) * rect.height())
        colour = QColor(series.color or "#6fb1ff")
        for position, value in zip(series.x, series.y):
            centre = rect.left() + self._x.to_unit(position) * rect.width()
            left = centre - width * count / 2.0 + index * width
            unit = scale.to_unit(value)
            if not math.isfinite(unit):
                continue
            top = rect.bottom() - unit * rect.height()
            painter.fillRect(QRectF(left, min(top, zero), width,
                                    abs(zero - top)), colour)

    def _line(self, painter: QPainter, rect: QRectF, series) -> None:
        painter.setPen(QPen(QColor(series.color or "#6fb1ff"), 1.8))
        previous = None
        for x, y in zip(series.x, series.y):
            point = self._point(rect, x, y, series.axis)
            if math.isfinite(point.x()) and math.isfinite(point.y()):
                if previous is not None:
                    painter.drawLine(previous, point)
                previous = point
            else:
                previous = None

    def _points(self, painter: QPainter, rect: QRectF, series) -> None:
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor(series.color or "#6fb1ff"))
        for x, y in zip(series.x, series.y):
            point = self._point(rect, x, y, series.axis)
            if math.isfinite(point.x()) and math.isfinite(point.y()):
                painter.drawEllipse(point, 3.2, 3.2)
        painter.setBrush(Qt.BrushStyle.NoBrush)

    def _legend(self, painter: QPainter, rect: QRectF) -> None:
        if len(self.chart.series) < 2:
            return
        names = [s.name + (" (right)" if s.axis == 1 else "")
                 for s in self.chart.series]
        metrics = painter.fontMetrics()
        # Sized to the longest entry: a fixed box clipped "free energy of the
        # transition" to "...transiti", which is worse than no legend.
        width = min(max(metrics.horizontalAdvance(n) for n in names) + 24,
                    max(rect.width() - 40, 80.0))
        y = rect.top() + 6
        for series, name in zip(self.chart.series, names):
            painter.fillRect(QRectF(rect.right() - width - 8, y + 3, 10, 10),
                             QColor(series.color or "#6fb1ff"))
            painter.setPen(QPen(self._ink_text))
            painter.drawText(QRectF(rect.right() - width + 8, y, width, 16),
                             int(Qt.AlignmentFlag.AlignLeft), name)
            y += 16
