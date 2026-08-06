"""The Analysis dock: pore, pockets, conservation and allostery.

Everything here is presentation and signals only. The heavy work lives in
:class:`piezo1.ui.analysis_controller.AnalysisController`, which runs it off the
GUI thread — a pocket search or a PRS scan takes seconds, and a frozen window
during a scientific calculation reads as a crash.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox,
                             QFormLayout, QGroupBox, QHBoxLayout, QLabel,
                             QPushButton, QSpinBox, QTableWidget,
                             QTableWidgetItem, QTabWidget, QVBoxLayout,
                             QWidget)

from ..profile_plot import Marker, ProfilePlot, Trace

__all__ = ["AnalysisPanel"]


def _table(headers: list[str]) -> QTableWidget:
    table = QTableWidget(0, len(headers))
    table.setHorizontalHeaderLabels(headers)
    table.verticalHeader().setVisible(False)
    table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    table.horizontalHeader().setStretchLastSection(True)
    return table


class AnalysisPanel(QWidget):
    """Tabs for the analyses that were previously CLI-only."""

    pore_requested = pyqtSignal()
    pockets_requested = pyqtSignal(int)
    conservation_requested = pyqtSignal()
    allostery_requested = pyqtSignal()
    residues_selected = pyqtSignal(object, str)
    focus_requested = pyqtSignal(object)
    color_requested = pyqtSignal(str, bool)      # which scalar, on/off
    pore_position_picked = pyqtSignal(float)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.tabs = QTabWidget()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        layout.addWidget(self.tabs)

        self.tabs.addTab(self._build_pore(), "Pore")
        self.tabs.addTab(self._build_pockets(), "Pockets")
        self.tabs.addTab(self._build_residue_scalars(), "Residue maps")

        self._pore = None
        self._pockets: list = []

    # ------------------------------------------------------------------ pore

    def _build_pore(self) -> QWidget:
        page = QWidget()
        box = QVBoxLayout(page)

        row = QHBoxLayout()
        self.pore_button = QPushButton("Compute pore profile")
        self.pore_button.setToolTip(
            "Radius of the largest sphere that fits at each height along the\n"
            "conduction axis, with the probe tethered near the axis.\n"
            "Without that leash the answer escapes to ~6000 A.")
        self.pore_button.clicked.connect(self.pore_requested.emit)
        row.addWidget(self.pore_button)
        self.hydro_check = QCheckBox("hydrophobicity")
        self.hydro_check.setChecked(True)
        self.hydro_check.setToolTip(
            "Overlay the CHAP hydrophobicity profile. Radius alone predicts "
            "conduction at AUROC 0.59; with hydrophobicity, 0.91.")
        self.hydro_check.toggled.connect(self._replot)
        row.addWidget(self.hydro_check)
        box.addLayout(row)

        self.pore_plot = ProfilePlot()
        self.pore_plot.position_clicked.connect(self.pore_position_picked.emit)
        box.addWidget(self.pore_plot, 1)

        self.pore_label = QLabel("no profile computed")
        self.pore_label.setWordWrap(True)
        self.pore_label.setTextFormat(Qt.TextFormat.RichText)
        box.addWidget(self.pore_label)

        self.constriction_table = _table(["z (Å)", "radius (Å)", "lining"])
        self.constriction_table.setToolTip(
            "Local minima in the radius profile — candidate gates.\n"
            "Select one to highlight the residues lining it.")
        self.constriction_table.itemSelectionChanged.connect(
            self._on_constriction)
        box.addWidget(self.constriction_table, 1)
        return page

    def set_pore(self, profile, hydration=None, hydrophobicity=None) -> None:
        self._pore = (profile, hydration, hydrophobicity)
        self._replot()

        if profile is None:
            self.pore_label.setText("no profile computed")
            self.constriction_table.setRowCount(0)
            return

        text = (f"<b>bottleneck {profile.bottleneck_radius:.2f} Å</b> "
                f"({profile.bottleneck_radius / 10:.3f} nm) at "
                f"z = {profile.bottleneck_z:.1f} Å")
        if hydration is not None and hydration.available:
            text += (f"<br>Rao 2019 score <b>{hydration.score:.2f}</b> "
                     f"(closed above 0.55) → <b>{hydration.verdict}</b>")
        self.pore_label.setText(text)

        rows = profile.constrictions()
        self.constriction_table.setRowCount(len(rows))
        for i, sl in enumerate(rows):
            lining = ", ".join(str(r) for r in sl.lining[:6])
            for col, value in enumerate((f"{sl.z:.1f}", f"{sl.radius:.2f}",
                                         lining)):
                self.constriction_table.setItem(i, col,
                                                QTableWidgetItem(value))
        self.constriction_table.resizeColumnsToContents()

    def _replot(self) -> None:
        if self._pore is None or self._pore[0] is None:
            self.pore_plot.clear()
            return
        profile, hydration, hydro = self._pore
        traces = [Trace("radius (Å)", profile.z, profile.radius,
                        "#6fb1ff", axis=0, fill=True)]
        markers = [Marker(profile.bottleneck_z, "bottleneck", "#f26d6d",
                          vertical=True),
                   # A water molecule is ~1.5 Å in radius; below this the pore
                   # is shut on sterics alone, whatever its chemistry.
                   Marker(1.5, "water radius 1.5 Å", "#8a919e", axis=0)]
        if self.hydro_check.isChecked() and hydro is not None:
            traces.append(Trace("hydrophobicity", profile.z, hydro,
                                "#f2a65a", axis=1))
        self.pore_plot.set_data(
            traces, markers, x_label="z along conduction axis (Å)",
            left_label="pore radius (Å)",
            right_label="Wimley–White (normalised)")

    def _on_constriction(self) -> None:
        rows = self.constriction_table.selectionModel().selectedRows()
        if not rows or self._pore is None:
            return
        sl = self._pore[0].constrictions()[rows[0].row()]
        if sl.lining:
            self.residues_selected.emit(list(sl.lining),
                                        f"constriction at z = {sl.z:.1f} Å")
            self.pore_plot.set_cursor(sl.z)

    # --------------------------------------------------------------- pockets

    def _build_pockets(self) -> QWidget:
        page = QWidget()
        box = QVBoxLayout(page)
        form = QFormLayout()
        self.pocket_count = QSpinBox()
        self.pocket_count.setRange(1, 40)
        self.pocket_count.setValue(10)
        self.pocket_count.setToolTip("How many of the largest pockets to keep")
        form.addRow("keep top", self.pocket_count)
        box.addLayout(form)

        self.pockets_button = QPushButton("Find pockets")
        self.pockets_button.setToolTip(
            "Delaunay alpha-sphere detection with a burial filter. Volumes are\n"
            "Monte-Carlo unions, not sums of overlapping spheres.")
        self.pockets_button.clicked.connect(
            lambda: self.pockets_requested.emit(self.pocket_count.value()))
        box.addWidget(self.pockets_button)

        self.pocket_table = _table(["#", "volume (Å³)", "buried", "residues"])
        self.pocket_table.setToolTip("Select a pocket to highlight and frame it")
        self.pocket_table.itemSelectionChanged.connect(self._on_pocket)
        box.addWidget(self.pocket_table, 1)

        self.pocket_label = QLabel("no pockets computed")
        self.pocket_label.setWordWrap(True)
        box.addWidget(self.pocket_label)
        return page

    def set_pockets(self, pockets: list) -> None:
        self._pockets = list(pockets)
        self.pocket_table.setRowCount(len(self._pockets))
        for i, pocket in enumerate(self._pockets):
            residues = ", ".join(str(r) for r in sorted(pocket.residues)[:6])
            values = (str(i + 1), f"{pocket.volume:.0f}",
                      f"{pocket.buriedness:.2f}", residues)
            for col, value in enumerate(values):
                self.pocket_table.setItem(i, col, QTableWidgetItem(value))
        self.pocket_table.resizeColumnsToContents()
        self.pocket_label.setText(
            f"{len(self._pockets)} pockets — select one to highlight it"
            if self._pockets else "no pockets found")

    def _on_pocket(self) -> None:
        rows = self.pocket_table.selectionModel().selectedRows()
        if not rows or not self._pockets:
            return
        pocket = self._pockets[rows[0].row()]
        self.residues_selected.emit(sorted(pocket.residues),
                                    f"pocket {rows[0].row() + 1}")
        self.focus_requested.emit(sorted(pocket.residues))

    # -------------------------------------------------- per-residue scalars

    def _build_residue_scalars(self) -> QWidget:
        page = QWidget()
        box = QVBoxLayout(page)

        conservation = QGroupBox("Conservation")
        inner = QVBoxLayout(conservation)
        self.conservation_button = QPushButton("Compute conservation")
        self.conservation_button.setToolTip(
            "Per-residue Shannon entropy over vertebrate PIEZO1 orthologs,\n"
            "one sequence per species. Fetches on first use.\n"
            "Positions covered by under 70% of orthologs are dropped.")
        self.conservation_button.clicked.connect(
            self.conservation_requested.emit)
        inner.addWidget(self.conservation_button)
        self.conservation_label = QLabel("not computed")
        self.conservation_label.setWordWrap(True)
        inner.addWidget(self.conservation_label)
        box.addWidget(conservation)

        allostery = QGroupBox("Allostery (perturbation response)")
        inner = QVBoxLayout(allostery)
        self.allostery_button = QPushButton("Scan coupling to the gate")
        self.allostery_button.setToolTip(
            "Needs normal modes — compute them in the Physics panel first.")
        self.allostery_button.clicked.connect(self.allostery_requested.emit)
        inner.addWidget(self.allostery_button)
        self.allostery_label = QLabel("not computed")
        self.allostery_label.setWordWrap(True)
        inner.addWidget(self.allostery_label)
        box.addWidget(allostery)

        colour = QGroupBox("Colour the model by")
        inner = QVBoxLayout(colour)
        self.scalar_combo = QComboBox()
        self.scalar_combo.addItem("— off —", "")
        inner.addWidget(self.scalar_combo)
        self.scalar_check = QCheckBox("apply to structure")
        self.scalar_check.setToolTip(
            "Colour the 3-D model by the selected per-residue value.\n"
            "Residues with no value take the map's minimum, not zero.")
        self.scalar_check.toggled.connect(self._emit_colour)
        self.scalar_combo.currentIndexChanged.connect(
            lambda _: self._emit_colour(self.scalar_check.isChecked()))
        inner.addWidget(self.scalar_check)
        self.scalar_label = QLabel("")
        self.scalar_label.setWordWrap(True)
        inner.addWidget(self.scalar_label)
        box.addWidget(colour)

        self.top_table = _table(["residue", "value"])
        self.top_table.itemSelectionChanged.connect(self._on_top_residue)
        box.addWidget(self.top_table, 1)
        return page

    def add_scalar(self, key: str, label: str) -> None:
        if self.scalar_combo.findData(key) < 0:
            self.scalar_combo.addItem(label, key)
        self.scalar_combo.setCurrentIndex(self.scalar_combo.findData(key))

    def current_scalar(self) -> str:
        return str(self.scalar_combo.currentData() or "")

    def _emit_colour(self, on: bool) -> None:
        self.color_requested.emit(self.current_scalar(), bool(on))

    def set_top_residues(self, pairs: list[tuple[int, float]],
                         note: str = "") -> None:
        self.top_table.setRowCount(len(pairs))
        for i, (residue, value) in enumerate(pairs):
            self.top_table.setItem(i, 0, QTableWidgetItem(str(residue)))
            self.top_table.setItem(i, 1, QTableWidgetItem(f"{value:.3f}"))
        self.top_table.resizeColumnsToContents()
        self.scalar_label.setText(note)

    def _on_top_residue(self) -> None:
        rows = self.top_table.selectionModel().selectedRows()
        if not rows:
            return
        item = self.top_table.item(rows[0].row(), 0)
        if item is not None:
            self.residues_selected.emit([int(item.text())], "selected residue")

    # ---------------------------------------------------------------- status

    def set_busy(self, busy: bool, what: str = "") -> None:
        for button in (self.pore_button, self.pockets_button,
                       self.conservation_button, self.allostery_button):
            button.setEnabled(not busy)
        if busy and what:
            self.set_message(what, f"computing {what}…")

    def set_message(self, which: str, text: str) -> None:
        target = {"conservation": self.conservation_label,
                  "allostery": self.allostery_label,
                  "pockets": self.pocket_label,
                  "pore": self.pore_label}.get(which)
        if target is not None:
            target.setText(text)
