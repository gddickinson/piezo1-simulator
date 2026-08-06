"""The parameter editor: every number a calculation uses, with its provenance.

One row per registered parameter, grouped by category, showing the current
value, the documented default, the unit, and the citation it came from. The
full reference is on the tooltip, so "where does 2.6 kJ/mol come from" is one
hover away rather than a grep through the source.

**Overridden rows are marked, and the marking is not decorative.** A number
computed with a changed parameter is not comparable with the values in
`docs/SCIENCE.md`, and the rest of the application enforces that: reports carry
a banner, and the documentation-claims verifier refuses to run at all. The
colouring here is the user-facing half of the same rule.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (QAbstractItemView, QCheckBox, QDialog,
                             QDialogButtonBox, QHBoxLayout, QLabel, QLineEdit,
                             QPushButton, QTreeWidget, QTreeWidgetItem,
                             QVBoxLayout)

from ..parameters import PARAMETERS

__all__ = ["ParametersDialog"]

COLUMNS = ["Parameter", "Value", "Default", "Unit", "Source"]
OVERRIDDEN = QColor("#f2a65a")
CITED = QColor("#7ed67e")
UNCITED = QColor("#8a919e")


class ParametersDialog(QDialog):
    """Browse and edit the registry."""

    changed = pyqtSignal()

    def __init__(self, references: dict | None = None, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Parameters — values, defaults and sources")
        self.resize(1040, 680)
        self.references = references or {}

        layout = QVBoxLayout(self)
        layout.addWidget(QLabel(
            "Every number the calculations use. Edit the <b>Value</b> column to "
            "change one; values outside the declared range are clamped.<br>"
            "<span style='color:#f2a65a'>Amber</span> marks a parameter that "
            "differs from its documented default — results computed with one "
            "are not comparable with the numbers in <tt>docs/SCIENCE.md</tt>."))

        row = QHBoxLayout()
        self.filter_edit = QLineEdit()
        self.filter_edit.setPlaceholderText(
            "Filter by name, key, unit or citation…")
        self.filter_edit.textChanged.connect(self._apply_filter)
        row.addWidget(self.filter_edit, 1)
        self.modified_only = QCheckBox("Show only modified")
        self.modified_only.toggled.connect(self._apply_filter)
        row.addWidget(self.modified_only)
        layout.addLayout(row)

        self.tree = QTreeWidget()
        self.tree.setColumnCount(len(COLUMNS))
        self.tree.setHeaderLabels(COLUMNS)
        self.tree.setAlternatingRowColors(True)
        self.tree.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked
                                  | QAbstractItemView.EditTrigger.SelectedClicked)
        self.tree.itemChanged.connect(self._on_edited)
        layout.addWidget(self.tree, 1)

        self.status = QLabel("")
        self.status.setWordWrap(True)
        layout.addWidget(self.status)

        buttons = QHBoxLayout()
        reset_one = QPushButton("Reset selected")
        reset_one.clicked.connect(self._reset_selected)
        buttons.addWidget(reset_one)
        reset_all = QPushButton("Reset all to documented defaults")
        reset_all.clicked.connect(self._reset_all)
        buttons.addWidget(reset_all)
        buttons.addStretch(1)
        box = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        box.rejected.connect(self.accept)
        box.accepted.connect(self.accept)
        buttons.addWidget(box)
        layout.addLayout(buttons)

        self._populate()

    # ------------------------------------------------------------- building

    def _populate(self) -> None:
        self.tree.blockSignals(True)
        self.tree.clear()
        for category in PARAMETERS.categories():
            parent = QTreeWidgetItem([category])
            parent.setFirstColumnSpanned(True)
            parent.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self.tree.addTopLevelItem(parent)
            for parameter in PARAMETERS.in_category(category):
                parent.addChild(self._row(parameter))
            parent.setExpanded(True)
        for column in range(len(COLUMNS)):
            self.tree.resizeColumnToContents(column)
        self.tree.blockSignals(False)
        self._refresh_status()

    def _row(self, parameter) -> QTreeWidgetItem:
        current = PARAMETERS.value(parameter.key)
        source = parameter.citation
        item = QTreeWidgetItem([parameter.name, f"{current:g}",
                                f"{parameter.default:g}", parameter.unit,
                                source])
        item.setData(0, Qt.ItemDataRole.UserRole, parameter.key)
        item.setFlags(item.flags() | Qt.ItemFlag.ItemIsEditable)

        bounds = ""
        if parameter.minimum is not None and parameter.maximum is not None:
            bounds = f"\nallowed range: {parameter.minimum:g} – {parameter.maximum:g}"
        reference = self.references.get(parameter.citation)
        provenance = (f"\n\nSource: {reference}" if reference
                      else f"\n\nSource: {parameter.citation}")
        if parameter.source_note:
            provenance += f"\n{parameter.source_note}"
        tip = (f"{parameter.key}\n{parameter.description}"
               f"\nkind: {parameter.kind}{bounds}{provenance}")
        for column in range(len(COLUMNS)):
            item.setToolTip(column, tip)
        self._paint(item, parameter)
        return item

    def _paint(self, item: QTreeWidgetItem, parameter) -> None:
        overridden = not PARAMETERS.is_default(parameter.key)
        for column in range(len(COLUMNS)):
            item.setForeground(column, QBrush(
                OVERRIDDEN if overridden
                else (CITED if parameter.cited else UNCITED)))
        font = item.font(0)
        font.setBold(overridden)
        for column in range(len(COLUMNS)):
            item.setFont(column, font)

    # -------------------------------------------------------------- editing

    def _on_edited(self, item: QTreeWidgetItem, column: int) -> None:
        if column != 1:
            return
        key = item.data(0, Qt.ItemDataRole.UserRole)
        if not key:
            return
        parameter = PARAMETERS.get(key)
        try:
            requested = float(item.text(1))
        except ValueError:
            self.status.setText(
                f"<span style='color:#f26d6d'>{item.text(1)!r} is not a "
                f"number — {parameter.name} left at "
                f"{PARAMETERS.value(key):g}</span>")
            self._rewrite(item, key, parameter)
            return

        applied = PARAMETERS.set_value(key, requested)
        if applied != requested:
            self.status.setText(
                f"<span style='color:#f2a65a'>{requested:g} is outside the "
                f"allowed range for {parameter.name}; clamped to "
                f"{applied:g} {parameter.unit}</span>")
        else:
            self.status.setText("")
        self._rewrite(item, key, parameter)
        self._refresh_status(keep=bool(self.status.text()))
        self.changed.emit()

    def _rewrite(self, item: QTreeWidgetItem, key: str, parameter) -> None:
        self.tree.blockSignals(True)
        item.setText(1, f"{PARAMETERS.value(key):g}")
        self._paint(item, parameter)
        self.tree.blockSignals(False)

    def _reset_selected(self) -> None:
        for item in self.tree.selectedItems():
            key = item.data(0, Qt.ItemDataRole.UserRole)
            if key:
                PARAMETERS.reset(key)
                self._rewrite(item, key, PARAMETERS.get(key))
        self._refresh_status()
        self.changed.emit()

    def _reset_all(self) -> None:
        PARAMETERS.reset()
        self._populate()
        self.changed.emit()

    # ------------------------------------------------------------- filtering

    def _apply_filter(self) -> None:
        needle = self.filter_edit.text().strip().lower()
        only_modified = self.modified_only.isChecked()
        for i in range(self.tree.topLevelItemCount()):
            parent = self.tree.topLevelItem(i)
            visible_children = 0
            for j in range(parent.childCount()):
                child = parent.child(j)
                key = child.data(0, Qt.ItemDataRole.UserRole)
                parameter = PARAMETERS.get(key)
                haystack = " ".join([parameter.key, parameter.name,
                                     parameter.unit, parameter.citation,
                                     parameter.category]).lower()
                show = (not needle or needle in haystack)
                if only_modified and PARAMETERS.is_default(key):
                    show = False
                child.setHidden(not show)
                visible_children += show
            parent.setHidden(visible_children == 0)

    def _refresh_status(self, keep: bool = False) -> None:
        if keep:
            return
        if PARAMETERS.modified:
            self.status.setText(
                f"<b style='color:#f2a65a'>{PARAMETERS.override_summary()}</b>"
                "<br>Reports built now will carry a warning banner, and the "
                "documentation verifier will refuse to run.")
        else:
            self.status.setText(
                f"{len(PARAMETERS)} parameters, all at their documented "
                f"defaults.")
