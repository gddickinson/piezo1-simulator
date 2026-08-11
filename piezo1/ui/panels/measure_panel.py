"""Click-to-measure: distances, angles and dihedrals in the 3D view.

Picking is **armed**, not always-on, because a click in the viewport already
means "tell me about this residue" and a measurement tool that quietly consumed
those clicks would break inspection. The cost of that choice is that a user who
has not found the Measure button sees clicks do nothing here — which is exactly
what happened — so the arm state is now visible in three places at once: the
button's own label, the hint under it, and the table, which shows a pick the
moment it is made rather than only when the measurement completes.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (QAbstractItemView, QComboBox, QFileDialog,
                             QHBoxLayout, QHeaderView, QLabel, QPushButton,
                             QTableWidget, QTableWidgetItem, QVBoxLayout,
                             QWidget)

from ...analysis.measure import MEASUREMENT_KINDS, MeasurementSet

__all__ = ["MeasurePanel", "PENDING_COLOR"]

#: Pending picks are dimmed so an incomplete selection is never read as a
#: result. The same blue the viewport draws the pick markers in.
PENDING_COLOR = QColor(120, 200, 255)


def _a(kind: str) -> str:
    """"a distance", "an angle" — the measurement kinds start with both."""
    return f"{'an' if kind[:1] in 'aeiou' else 'a'} {kind}"


class MeasurePanel(QWidget):
    """Drives a :class:`MeasurementSet` from clicks in the viewport."""

    mode_changed = pyqtSignal(bool)
    measurements_changed = pyqtSignal()
    status = pyqtSignal(str)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.set = MeasurementSet()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)

        row = QHBoxLayout()
        self.kind_combo = QComboBox()
        self.kind_combo.setToolTip('Distance needs two atoms, angle three, dihedral four.')
        self.kind_combo.addItems(sorted(MEASUREMENT_KINDS))
        self.kind_combo.setCurrentText("distance")
        self.kind_combo.currentTextChanged.connect(self._set_kind)
        row.addWidget(self.kind_combo, 1)

        self.arm_button = QPushButton("Start picking")

        self.arm_button.setToolTip(
            "Clicks in the 3-D view normally identify a residue. Press this to\n"
            "send them here instead, then click atoms in the viewport: each\n"
            "appears in the table below and is marked in blue on the model.\n"
            "The regression case is the C2411–C2415 disulfide at 2.04 Å.")
        self.arm_button.setCheckable(True)
        self.arm_button.toggled.connect(self._arm)
        row.addWidget(self.arm_button)
        layout.addLayout(row)

        self.hint = QLabel()
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color:#9aa3b2; font-size:11px;")
        layout.addWidget(self.hint)
        self._update_hint()

        self.table = QTableWidget(0, 3)
        self.table.setHorizontalHeaderLabels(["Selection", "Value", "Units"])
        self.table.verticalHeader().setVisible(False)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.horizontalHeader().setSectionResizeMode(
            0, QHeaderView.ResizeMode.Stretch)
        layout.addWidget(self.table, 1)

        row = QHBoxLayout()
        for label, slot in (("Delete", self._delete), ("Clear all", self._clear),
                            ("Copy", self._copy), ("Export CSV…", self._export)):
            b = QPushButton(label)
            b.clicked.connect(slot)
            row.addWidget(b)
        layout.addLayout(row)

    # ------------------------------------------------------------- control

    @property
    def armed(self) -> bool:
        return self.arm_button.isChecked()

    def _set_kind(self, kind: str) -> None:
        self.set.set_kind(kind)
        self._update_hint()
        self._refill()
        self.measurements_changed.emit()

    def _arm(self, on: bool) -> None:
        self.arm_button.setText("Picking — click atoms" if on else "Start picking")
        if not on:
            self.set.clear_pending()
        self._update_hint()
        self._refill()
        self.mode_changed.emit(on)
        self.measurements_changed.emit()

    def _update_hint(self) -> None:
        """What to do next, which depends entirely on whether picking is armed.

        The old text — "Pick 2 atoms for a distance" — described the goal and
        not the step, so it read as an instruction that clicking would work.
        """
        kind = self.set.kind
        needed = self.set.required
        if not self.armed:
            self.hint.setText(
                f"Press <b>Start picking</b>, then click {needed} atoms in the "
                f"3-D view for {_a(kind)}. Until then a click identifies the "
                f"residue in the status bar instead.")
            return
        picked = len(self.set.pending)
        self.hint.setText(
            f"<b>Click atoms in the 3-D view.</b> {picked} of {needed} picked "
            f"for {_a(kind)}; picks are blue on the model.")

    # ------------------------------------------------------------- picking

    def add_pick(self, index: int, position, label: str) -> None:
        """Route a viewport click into the measurement set."""
        if not self.armed:
            return
        before = len(self.set.pending)
        result = self.set.add_atom(index, position, label)
        if len(self.set.pending) == before and result is None:
            self.status.emit("same atom picked twice — ignored")
            return
        if result is None:
            remaining = self.set.required - len(self.set.pending)
            self.status.emit(f"{label} selected; {remaining} more to go")
        else:
            self.status.emit(str(result))
        self._update_hint()
        self._refill()
        self.measurements_changed.emit()

    # -------------------------------------------------------------- table

    def _refill(self) -> None:
        """Completed measurements, then whatever is picked but not yet used.

        The pending rows are the point: without them the panel stayed empty
        until the last atom of a measurement was clicked, so a user who had
        selected one atom of two had no confirmation anywhere but the status
        bar — and no way to see *which* atom the tool thought they meant.
        """
        rows = self.set.measurements
        pending = list(self.set.pending_labels)
        self.table.setRowCount(len(rows) + (1 if pending else 0))

        for r, m in enumerate(rows):
            for c, text in enumerate((" – ".join(m.labels),
                                      f"{m.value:.2f}", m.units)):
                item = QTableWidgetItem(text)
                if c == 1:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(r, c, item)

        if not pending:
            return
        remaining = self.set.required - len(pending)
        cells = (" – ".join(pending), "…",
                 f"{remaining} more for {_a(self.set.kind)}")
        for c, text in enumerate(cells):
            item = QTableWidgetItem(text)
            item.setForeground(QBrush(PENDING_COLOR))
            if c == 1:
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                      | Qt.AlignmentFlag.AlignVCenter)
            self.table.setItem(len(rows), c, item)

    def _delete(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()},
                      reverse=True)
        # The pending row sits past the end of `measurements`; deleting it
        # means abandoning the half-made selection, not removing a result.
        for r in rows:
            if r < len(self.set.measurements):
                self.set.remove(r)
            else:
                self.set.clear_pending()
        self._update_hint()
        self._refill()
        self.measurements_changed.emit()

    def _clear(self) -> None:
        self.set.clear()
        self._update_hint()
        self._refill()
        self.measurements_changed.emit()

    def _copy(self) -> None:
        from PyQt6.QtWidgets import QApplication
        text = self.set.to_text()
        QApplication.clipboard().setText(text)
        self.status.emit(f"copied {len(self.set.measurements)} measurements")

    def _export(self) -> None:
        if not self.set.measurements:
            self.status.emit("nothing to export")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export measurements",
                                              "measurements.csv", "CSV (*.csv)")
        if path:
            with open(path, "w") as fh:
                fh.write(self.set.to_csv())
            self.status.emit(f"wrote {path}")

    # ------------------------------------------------------------- overlay

    def overlay_labels(self) -> list:
        """World-anchored labels for the viewport, plus any pending picks."""
        out = [(m.anchor, f"{m.value:.2f} {m.units}", (255, 214, 61))
               for m in self.set.measurements]
        for label, position in zip(self.set.pending_labels,
                                   self.set.pending_positions):
            out.append((position, label, (120, 200, 255)))
        return out

    def highlighted_atoms(self) -> list[int]:
        picked = list(self.set.pending)
        for m in self.set.measurements:
            picked.extend(m.atoms)
        return picked
