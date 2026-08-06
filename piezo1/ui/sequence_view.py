"""A monospace sequence grid with click-and-drag selection.

Drawn with QPainter rather than assembled from per-residue widgets: PIEZO1 is
2521 residues and a widget per residue is 2521 widgets to lay out on every
resize. Painting is one pass over the visible rows.

Selection is reported in **residue numbers**, never in offsets into the string.
A structure sequence starts at residue 570 and has gaps in it, so an offset is
only meaningful alongside the sequence it came from — and passing offsets to
something that expects residue numbers is the numbering bug this project keeps
finding in other guises.
"""

from __future__ import annotations

from dataclasses import dataclass

from PyQt6.QtCore import QRectF, Qt, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPen
from PyQt6.QtWidgets import QSizePolicy, QWidget

__all__ = ["SequenceView", "ResidueStyle"]

#: Fallback residue colouring, by broad chemistry.
CHEMISTRY = {
    **{a: "#6fb1ff" for a in "AVLIMFWPC"},      # hydrophobic
    **{a: "#7ed67e" for a in "STNQGY"},          # polar
    **{a: "#f26d6d" for a in "DE"},              # acidic
    **{a: "#f2a65a" for a in "KRH"},             # basic
}


@dataclass
class ResidueStyle:
    """Per-residue decoration supplied by the window."""

    background: str = ""
    foreground: str = ""
    underline: str = ""
    tooltip: str = ""


class SequenceView(QWidget):
    """Wrapped monospace sequence with numbering, selection and a DNA track."""

    selection_changed = pyqtSignal(int, int)      # first, last residue
    residue_hovered = pyqtSignal(int)

    GUTTER = 78
    PAD = 8

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.sequence = None
        self.styles: dict[int, ResidueStyle] = {}
        self.show_dna = False
        self.selection: tuple[int, int] | None = None
        self._anchor: int | None = None
        self._columns = 60
        self._metrics: QFontMetricsF | None = None

        font = QFont("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self.setFont(font)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    # ------------------------------------------------------------------ data

    def set_sequence(self, sequence, styles: dict[int, ResidueStyle] | None = None,
                     show_dna: bool = False) -> None:
        self.sequence = sequence
        self.styles = styles or {}
        self.show_dna = bool(show_dna and sequence is not None and sequence.dna)
        self.selection = None
        self._anchor = None
        self._relayout()
        self.update()

    def set_styles(self, styles: dict[int, ResidueStyle]) -> None:
        self.styles = styles or {}
        self.update()

    def select(self, first: int, last: int | None = None) -> None:
        if self.sequence is None:
            return
        last = first if last is None else last
        self.selection = (min(first, last), max(first, last))
        self.update()
        self.selection_changed.emit(*self.selection)

    def selected_residues(self) -> list[int]:
        if self.sequence is None or self.selection is None:
            return []
        lo, hi = self.selection
        return [p for p in self.sequence.positions if lo <= p <= hi]

    # ---------------------------------------------------------------- layout

    def _cell(self) -> float:
        metrics = QFontMetricsF(self.font())
        self._metrics = metrics
        return max(metrics.horizontalAdvance("W"), 7.0) + 2.0

    def _row_height(self) -> float:
        metrics = self._metrics or QFontMetricsF(self.font())
        rows = 3.0 if self.show_dna else 1.0
        return metrics.height() * rows + 14.0

    def _relayout(self) -> None:
        cell = self._cell()
        usable = max(self.width() - self.GUTTER - 2 * self.PAD, cell * 10)
        # Round to a multiple of ten so the ruler above the sequence lands on
        # tens, which is how anyone actually counts along a sequence.
        self._columns = max(10, int(usable / cell) // 10 * 10)
        if self.sequence is not None:
            rows = -(-len(self.sequence) // self._columns)
            self.setMinimumHeight(int(rows * self._row_height() + 2 * self.PAD + 18))

    def resizeEvent(self, event) -> None:            # noqa: N802
        super().resizeEvent(event)
        self._relayout()

    def _residue_at(self, x: float, y: float) -> int | None:
        if self.sequence is None:
            return None
        cell = self._cell()
        row = int((y - self.PAD - 18) // self._row_height())
        col = int((x - self.GUTTER - self.PAD) // cell)
        if row < 0 or col < 0 or col >= self._columns:
            return None
        index = row * self._columns + col
        if 0 <= index < len(self.sequence):
            return self.sequence.positions[index]
        return None

    # ----------------------------------------------------------------- paint

    def paintEvent(self, event) -> None:             # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing)
        painter.fillRect(self.rect(), QColor("#14171c"))
        if self.sequence is None:
            painter.setPen(QPen(QColor("#6f7684")))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                             "no sequence loaded")
            painter.end()
            return

        painter.setFont(self.font())
        cell = self._cell()
        metrics = self._metrics
        line = metrics.height()
        row_h = self._row_height()
        selection = self.selection

        for index, (letter, residue) in enumerate(
                zip(self.sequence.letters, self.sequence.positions)):
            row, col = divmod(index, self._columns)
            x = self.GUTTER + self.PAD + col * cell
            y = self.PAD + 18 + row * row_h

            if col == 0:
                painter.setPen(QPen(QColor("#6f7684")))
                painter.drawText(QRectF(self.PAD, y, self.GUTTER - 12, line),
                                 Qt.AlignmentFlag.AlignRight
                                 | Qt.AlignmentFlag.AlignVCenter, str(residue))

            style = self.styles.get(residue)
            box = QRectF(x, y, cell, line)
            selected = selection is not None and selection[0] <= residue <= selection[1]
            if selected:
                painter.fillRect(box, QColor(122, 167, 255, 90))
            elif style and style.background:
                colour = QColor(style.background)
                colour.setAlpha(120)
                painter.fillRect(box, colour)

            if style and style.foreground:
                painter.setPen(QPen(QColor(style.foreground)))
            else:
                painter.setPen(QPen(QColor(CHEMISTRY.get(letter, "#c8ccd4"))))
            painter.drawText(box, Qt.AlignmentFlag.AlignCenter, letter)

            if style and style.underline:
                painter.setPen(QPen(QColor(style.underline), 2.0))
                painter.drawLine(QRectF(x + 1, y + line - 1, cell - 2, 1).topLeft(),
                                 QRectF(x + 1, y + line - 1, cell - 2, 1).topRight())

            if self.show_dna:
                codon = self.sequence.codon(residue)
                if codon:
                    small = QFont(self.font())
                    small.setPointSizeF(max(6.0, self.font().pointSizeF() * 0.62))
                    painter.setFont(small)
                    painter.setPen(QPen(QColor("#7ed67e" if selected else "#5d6470")))
                    painter.drawText(QRectF(x - cell * 0.6, y + line, cell * 2.2,
                                            line * 1.8),
                                     Qt.AlignmentFlag.AlignCenter, codon)
                    painter.setFont(self.font())

            # Ruler tick every ten residues, above the first row of each block.
            if residue % 10 == 0:
                painter.setPen(QPen(QColor("#3a4048")))
                painter.drawText(QRectF(x - cell, y - 14, cell * 3, 12),
                                 Qt.AlignmentFlag.AlignCenter, str(residue))
        painter.end()

    # ----------------------------------------------------------------- input

    def mousePressEvent(self, event) -> None:        # noqa: N802
        residue = self._residue_at(event.position().x(), event.position().y())
        if residue is None:
            return
        self._anchor = residue
        self.select(residue, residue)

    def mouseMoveEvent(self, event) -> None:         # noqa: N802
        residue = self._residue_at(event.position().x(), event.position().y())
        if residue is None:
            return
        self.residue_hovered.emit(residue)
        if self._anchor is not None and event.buttons():
            self.select(self._anchor, residue)

    def mouseReleaseEvent(self, event) -> None:      # noqa: N802
        self._anchor = None
