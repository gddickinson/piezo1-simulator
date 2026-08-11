"""The topology diagram as a widget — Guo & MacKinnon 2017, Figure 3a and 3b.

Draws one protomer's membrane topology: the transmembrane helices as cylinders
crossing a membrane band, the loops arcing above and below, the cap as a box
and the beam and cuff elements below. Groups of four helices — the 4-TM units
the whole architecture is built from — can be boxed and labelled, which is what
Figure 3b's red boxes do.

Three things it does that the published figure cannot:

* **It follows the loaded structure.** Every element is marked resolved or not
  from the entry's own coordinates, and an unresolved helix is drawn hollow.
  Figure 3a greys out TM1-12 for 6B3R; load 7WLT and a different set greys out.
* **The boxes are a selection.** Boxing a unit selects its residues on the
  3-D model, so "which part of the blade is this" is answered by looking
  rather than by counting helices in a picture.
* **It is one hover away from the numbers.** Every element carries its residue
  range in the structure's own numbering, and the numbering is named, because
  the paper is in mouse and most of this project's entries are not.

QPainter on a plain widget, in the idiom of :mod:`piezo1.ui.profile_plot` and
:mod:`piezo1.ui.sequence_view`: one plot type, has to repaint inside a window
at interactive rates, and must match the application's dark theme. The layout
itself is not decided here — it comes from :func:`piezo1.analysis.topology
.build_topology`, so the diagram and any analysis of the same architecture
cannot disagree.
"""

from __future__ import annotations

from PyQt6.QtCore import QPointF, QRectF, Qt, pyqtSignal
from PyQt6.QtGui import (QBrush, QColor, QFont, QPainter, QPainterPath,
                         QPen)
from PyQt6.QtWidgets import QSizePolicy, QToolTip, QWidget

from ..analysis.topology import Topology, build_topology, unit_extent

__all__ = ["TopologyView", "MEMBRANE_COLOR", "UNRESOLVED_ALPHA"]

#: The membrane band. Two lines and a wash, matching Figure 4a's grey rules
#: rather than a solid slab, so the helices stay readable through it.
MEMBRANE_COLOR = "#2a3038"
MEMBRANE_EDGE = "#4a525e"
#: An unresolved helix is drawn hollow at this alpha rather than omitted —
#: omitting it would silently renumber the diagram, which is the one thing a
#: topology figure must never do.
UNRESOLVED_ALPHA = 70
SELECTION_COLOR = "#ff5f56"
LABEL_COLOR = "#c8ccd4"
MUTED_COLOR = "#8a919e"


class TopologyView(QWidget):
    """Figure 3's diagram, for the loaded structure."""

    #: (residue_lo, residue_hi, numbering) for anything clicked or boxed.
    residues_selected = pyqtSignal(int, int, str)
    #: A short line for the status bar.
    status = pyqtSignal(str)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(280)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.topology: Topology | None = None
        #: Unit indices to box, as Figure 3b does. Empty draws no boxes.
        self.boxed: set[int] = set()
        #: Element kinds to draw. The cuff and beam are off by default: Figure
        #: 3a shows them, and they crowd a narrow window.
        self.shown_kinds: set[str] = {"tm_helix", "loop", "box", "bar", "cuff"}
        self._hover: int | None = None
        self._scale = 1.0
        self._origin = QPointF(0.0, 0.0)

    # ------------------------------------------------------------- contents

    def set_topology(self, topology: Topology | None) -> None:
        self.topology = topology
        self._hover = None
        self.update()

    def set_structure(self, structure, reference: str = "mouse") -> None:
        """Rebuild the diagram for a structure, or for none.

        ``reference`` is the committed UniProt resource whose transmembrane
        table defines the helices. It must match the structure's numbering: a
        human entry read with the mouse table would place every helix a few
        residues out and the diagram would look entirely plausible.
        """
        self.set_topology(build_topology(reference, structure))

    def set_boxed(self, units) -> None:
        self.boxed = {int(u) for u in units}
        self.update()
        if self.topology is not None and self.boxed:
            self.residues_selected.emit(*self.boxed_range())

    def boxed_range(self) -> tuple[int, int, str]:
        """Residue span covered by the boxed units, in the diagram's numbering.

        The span from the first boxed unit's start to the last one's end —
        contiguous even when the selection is not, which is deliberate: the
        3-D view highlights a range, and a gapped selection drawn as one range
        would be a lie about what is highlighted. The window says so.
        """
        if self.topology is None or not self.boxed:
            return (0, 0, "")
        helices = [e for e in self.topology.of_kind("tm_helix")
                   if e.unit in self.boxed]
        if not helices:
            return (0, 0, self.topology.numbering)
        return (min(e.start for e in helices), max(e.end for e in helices),
                self.topology.numbering)

    def available_units(self) -> list[int]:
        if self.topology is None:
            return []
        return sorted({e.unit for e in self.topology.of_kind("tm_helix")
                       if e.unit is not None})

    # -------------------------------------------------------------- layout

    def _layout(self) -> tuple[float, QPointF]:
        """Scale and origin mapping layout units to pixels."""
        if self.topology is None:
            return 1.0, QPointF(0.0, 0.0)
        x0, x1 = self.topology.meta["x_range"]
        margin = 18.0
        ys = [v for e in self.topology.elements for v in (e.y0, e.y1)]
        y0, y1 = (min(ys), max(ys)) if ys else (-1.0, 1.0)
        span_x = max(x1 - x0, 1e-6)
        span_y = max(y1 - y0, 1e-6)
        scale = min((self.width() - 2 * margin) / span_x,
                    (self.height() - 2 * margin - 16.0) / span_y)
        used_x = span_x * scale
        used_y = span_y * scale
        origin = QPointF((self.width() - used_x) / 2 - x0 * scale,
                         (self.height() - used_y) / 2 + y1 * scale)
        return scale, origin

    def _rect(self, element) -> QRectF:
        scale, origin = self._scale, self._origin
        left = origin.x() + element.x0 * scale
        right = origin.x() + element.x1 * scale
        # y is inverted: positive is extracellular, which is up on screen.
        top = origin.y() - element.y1 * scale
        bottom = origin.y() - element.y0 * scale
        return QRectF(left, min(top, bottom), max(right - left, 1.0),
                      max(abs(bottom - top), 1.0))

    # ------------------------------------------------------------- painting

    def paintEvent(self, event) -> None:            # noqa: N802 (Qt naming)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.fillRect(self.rect(), QColor("#14181d"))
        if self.topology is None:
            painter.setPen(QColor(MUTED_COLOR))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "No structure loaded")
            painter.end()
            return

        self._scale, self._origin = self._layout()
        self._draw_membrane(painter)
        self._draw_boxes(painter)
        for element in self.topology.elements:
            if element.kind not in self.shown_kinds:
                continue
            if element.kind == "loop":
                self._draw_loop(painter, element)
        for index, element in enumerate(self.topology.elements):
            if element.kind not in self.shown_kinds or element.kind == "loop":
                continue
            self._draw_block(painter, element, hovered=index == self._hover)
        self._draw_legend(painter)
        painter.end()

    def _draw_membrane(self, painter: QPainter) -> None:
        half = self.topology.meta["membrane_half"]
        x0, x1 = self.topology.meta["x_range"]
        top = self._origin.y() - half * self._scale
        bottom = self._origin.y() + half * self._scale
        band = QRectF(self._origin.x() + x0 * self._scale - 10.0, top,
                      (x1 - x0) * self._scale + 20.0, bottom - top)
        painter.fillRect(band, QColor(MEMBRANE_COLOR))
        painter.setPen(QPen(QColor(MEMBRANE_EDGE), 1.0))
        painter.drawLine(QPointF(band.left(), top), QPointF(band.right(), top))
        painter.drawLine(QPointF(band.left(), bottom),
                         QPointF(band.right(), bottom))
        painter.setPen(QColor(MUTED_COLOR))
        painter.setFont(QFont(self.font().family(), 8))
        painter.drawText(QPointF(band.left() + 2, top - 4), "extracellular")
        painter.drawText(QPointF(band.left() + 2, bottom + 12), "cytoplasmic")

    def _draw_boxes(self, painter: QPainter) -> None:
        """Figure 3b's red boxes around the selected 4-TM units."""
        if not self.boxed:
            return
        extents = unit_extent(self.topology)
        half = self.topology.meta["membrane_half"]
        pen = QPen(QColor(SELECTION_COLOR), 1.6)
        painter.setPen(pen)
        painter.setBrush(QBrush(QColor(255, 95, 86, 26)))
        painter.setFont(QFont(self.font().family(), 8, QFont.Weight.Bold))
        for unit in sorted(self.boxed):
            if unit not in extents:
                continue
            lo, hi = extents[unit]
            left = self._origin.x() + lo * self._scale - 4.0
            right = self._origin.x() + hi * self._scale + 4.0
            top = self._origin.y() - (half + 0.35) * self._scale
            bottom = self._origin.y() + (half + 0.35) * self._scale
            painter.drawRoundedRect(QRectF(left, top, right - left,
                                           bottom - top), 4.0, 4.0)
            painter.drawText(QPointF(left + 2, top - 3), f"THU{unit}")

    def _draw_loop(self, painter: QPainter, element) -> None:
        rect = self._rect(element)
        colour = QColor(element.color)
        if element.resolved is False:
            colour.setAlpha(UNRESOLVED_ALPHA)
            painter.setPen(QPen(colour, 1.2, Qt.PenStyle.DotLine))
        else:
            painter.setPen(QPen(colour, 1.2))
        path = QPainterPath()
        start = QPointF(rect.left(), self._origin.y()
                        - element.y0 * self._scale)
        end = QPointF(rect.right(), self._origin.y()
                      - element.y0 * self._scale)
        apex = self._origin.y() - element.y1 * self._scale
        path.moveTo(start)
        path.cubicTo(QPointF(start.x(), apex), QPointF(end.x(), apex), end)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawPath(path)

    def _draw_block(self, painter: QPainter, element,
                    hovered: bool = False) -> None:
        rect = self._rect(element)
        colour = QColor(element.color)
        if element.resolved is False:
            colour.setAlpha(UNRESOLVED_ALPHA)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.setPen(QPen(colour, 1.2, Qt.PenStyle.DashLine))
        else:
            painter.setBrush(QBrush(colour))
            painter.setPen(QPen(colour.darker(160), 1.0))
        if hovered:
            painter.setPen(QPen(QColor("#ffffff"), 1.6))
        radius = 3.0 if element.kind == "tm_helix" else 2.0
        painter.drawRoundedRect(rect, radius, radius)

        if element.kind == "tm_helix" and rect.width() > 9:
            painter.setPen(QColor("#0d1014") if element.resolved is not False
                           else QColor(MUTED_COLOR))
            painter.setFont(QFont(self.font().family(), 7))
            painter.drawText(rect, Qt.AlignmentFlag.AlignCenter,
                             str(element.helix))
        elif element.kind in ("box", "bar", "cuff"):
            painter.setPen(QColor(LABEL_COLOR))
            painter.setFont(QFont(self.font().family(), 8))
            painter.drawText(rect, Qt.AlignmentFlag.AlignCenter, element.label)

    def _draw_legend(self, painter: QPainter) -> None:
        painter.setPen(QColor(MUTED_COLOR))
        painter.setFont(QFont(self.font().family(), 8))
        unresolved = self.topology.unresolved_helices
        parts = [self.topology.summary()]
        if unresolved:
            parts.append(f"dashed = not modelled ({len(unresolved)} helices)")
        parts.append(f"numbering: {self.topology.numbering}")
        painter.drawText(QRectF(8, self.height() - 16, self.width() - 16, 14),
                         Qt.AlignmentFlag.AlignLeft, "  |  ".join(parts))

    # -------------------------------------------------------------- picking

    def _element_at(self, position) -> int | None:
        if self.topology is None:
            return None
        for index, element in enumerate(self.topology.elements):
            if element.kind not in self.shown_kinds or element.kind == "loop":
                continue
            if self._rect(element).adjusted(-1, -1, 1, 1).contains(position):
                return index
        return None

    def mouseMoveEvent(self, event) -> None:        # noqa: N802 (Qt naming)
        index = self._element_at(event.position())
        if index != self._hover:
            self._hover = index
            self.update()
        if index is None:
            QToolTip.hideText()
            return
        element = self.topology.elements[index]
        state = ("not modelled in this entry" if element.resolved is False
                 else f"{element.n_modelled} C-alpha modelled")
        unit = f", THU{element.unit}" if element.unit else ""
        QToolTip.showText(event.globalPosition().toPoint(),
                          f"{element.label}{unit}\n"
                          f"{element.start}-{element.end} "
                          f"({self.topology.numbering})\n{state}", self)

    def mousePressEvent(self, event) -> None:       # noqa: N802 (Qt naming)
        index = self._element_at(event.position())
        if index is None:
            return
        element = self.topology.elements[index]
        if (event.button() == Qt.MouseButton.LeftButton
                and element.unit is not None
                and event.modifiers() & Qt.KeyboardModifier.ShiftModifier):
            # Shift-click toggles the whole unit's box, which is the quick way
            # to reproduce Figure 3b without going to the checkboxes.
            self.boxed ^= {element.unit}
            self.update()
            if self.boxed:
                self.residues_selected.emit(*self.boxed_range())
            return
        self.residues_selected.emit(element.start, element.end,
                                    self.topology.numbering)
        self.status.emit(
            f"{element.label}: {element.start}-{element.end} "
            f"({self.topology.numbering} numbering)"
            + ("" if element.resolved is not False
               else " — not modelled in this entry"))

    def leaveEvent(self, event) -> None:            # noqa: N802 (Qt naming)
        self._hover = None
        self.update()
