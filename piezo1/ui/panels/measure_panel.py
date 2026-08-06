"""Click-to-measure: distances, angles and dihedrals in the 3D view."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QAbstractItemView, QComboBox, QFileDialog,
                             QHBoxLayout, QHeaderView, QLabel, QPushButton,
                             QTableWidget, QTableWidgetItem, QVBoxLayout,
                             QWidget)

from ...analysis.measure import MEASUREMENT_KINDS, MeasurementSet

__all__ = ["MeasurePanel"]


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

        self.arm_button = QPushButton("Measure")

        self.arm_button.setToolTip('Arm picking, then click atoms in the viewport.\nThe regression case is the C2411-C2415 disulfide at 2.04 A.')
        self.arm_button.setCheckable(True)
        self.arm_button.toggled.connect(self._arm)
        row.addWidget(self.arm_button)
        layout.addLayout(row)

        self.hint = QLabel("Pick 2 atoms for a distance.")
        self.hint.setWordWrap(True)
        self.hint.setStyleSheet("color:#9aa3b2; font-size:11px;")
        layout.addWidget(self.hint)

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
        n = self.set.required
        self.hint.setText(f"Pick {n} atoms for a{'n' if kind[0] in 'aeiou' else ''} {kind}.")
        self.measurements_changed.emit()

    def _arm(self, on: bool) -> None:
        self.arm_button.setText("Measuring — click atoms" if on else "Measure")
        if not on:
            self.set.clear_pending()
        self.mode_changed.emit(on)
        self.measurements_changed.emit()

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
            self._refill()
        self.measurements_changed.emit()

    # -------------------------------------------------------------- table

    def _refill(self) -> None:
        rows = self.set.measurements
        self.table.setRowCount(len(rows))
        for r, m in enumerate(rows):
            for c, text in enumerate((" – ".join(m.labels),
                                      f"{m.value:.2f}", m.units)):
                item = QTableWidgetItem(text)
                if c == 1:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                self.table.setItem(r, c, item)

    def _delete(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()},
                      reverse=True)
        for r in rows:
            self.set.remove(r)
        self._refill()
        self.measurements_changed.emit()

    def _clear(self) -> None:
        self.set.clear()
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
