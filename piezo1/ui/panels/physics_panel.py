"""Physics controls: dome geometry, normal modes and mode animation.

This panel is where the application stops being a viewer and starts being an
instrument. It measures the membrane dome from the coordinates on screen,
builds an elastic network model, and animates the low-frequency modes — with
each mode labelled by its C3 symmetry, because only the symmetric ones can
couple to isotropic membrane tension.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QComboBox, QDoubleSpinBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QProgressBar, QPushButton,
                             QSlider, QSpinBox, QVBoxLayout, QWidget)

__all__ = ["PhysicsPanel"]


class PhysicsPanel(QWidget):
    """Dome measurement and elastic-network mode analysis."""

    measure_dome_requested = pyqtSignal()
    compute_modes_requested = pyqtSignal(dict)
    mode_selected = pyqtSignal(int)
    animate_toggled = pyqtSignal(bool)
    amplitude_changed = pyqtSignal(float)
    color_by_mode_requested = pyqtSignal(bool)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(10)

        # ------------------------------------------------------------- dome
        box = QGroupBox("Membrane dome")
        v = QVBoxLayout(box)
        self.measure_button = QPushButton("Measure dome geometry")
        self.measure_button.clicked.connect(self.measure_dome_requested.emit)
        v.addWidget(self.measure_button)
        self.dome_label = QLabel("Not measured.")
        self.dome_label.setWordWrap(True)
        self.dome_label.setTextFormat(Qt.TextFormat.RichText)
        self.dome_label.setStyleSheet(
            "color:#c3cad8; background:#161a24; border:1px solid #262c3a;"
            "border-radius:5px; padding:7px; font-size:11px;")
        v.addWidget(self.dome_label)
        layout.addWidget(box)

        # ---------------------------------------------------- elastic model
        box = QGroupBox("Elastic network model")
        form = QFormLayout(box)

        self.cutoff_spin = QDoubleSpinBox()
        self.cutoff_spin.setRange(8.0, 25.0)
        self.cutoff_spin.setValue(15.0)
        self.cutoff_spin.setSuffix(" Å")
        self.cutoff_spin.setToolTip("C-alpha contact cutoff for the spring network")
        form.addRow("Cutoff", self.cutoff_spin)

        self.spring_combo = QComboBox()
        self.spring_combo.addItems(["inverse_square", "uniform", "inverse_sixth"])
        self.spring_combo.setToolTip(
            "Distance weighting of the springs. inverse_square (parameter-free "
            "ANM) avoids the spurious stiffening a hard cutoff causes at the "
            "edge of a large hollow assembly such as the PIEZO1 blade.")
        form.addRow("Spring model", self.spring_combo)

        self.nmodes_spin = QSpinBox()
        self.nmodes_spin.setRange(6, 60)
        self.nmodes_spin.setValue(24)
        form.addRow("Modes", self.nmodes_spin)

        self.compute_button = QPushButton("Compute normal modes")
        self.compute_button.clicked.connect(self._emit_compute)
        form.addRow(self.compute_button)

        self.progress = QProgressBar()
        self.progress.setVisible(False)
        self.progress.setRange(0, 0)
        form.addRow(self.progress)
        layout.addWidget(box)

        # ------------------------------------------------------------ modes
        box = QGroupBox("Mode animation")
        v = QVBoxLayout(box)
        self.mode_combo = QComboBox()
        self.mode_combo.setEnabled(False)
        self.mode_combo.currentIndexChanged.connect(self._on_mode)
        v.addWidget(self.mode_combo)

        self.mode_label = QLabel("Compute modes to begin.")
        self.mode_label.setWordWrap(True)
        self.mode_label.setTextFormat(Qt.TextFormat.RichText)
        self.mode_label.setStyleSheet("color:#9aa3b2; font-size:11px;")
        v.addWidget(self.mode_label)

        row = QHBoxLayout()
        row.addWidget(QLabel("Amplitude"))
        self.amp_slider = QSlider(Qt.Orientation.Horizontal)
        self.amp_slider.setRange(1, 60)
        self.amp_slider.setValue(18)
        self.amp_slider.valueChanged.connect(
            lambda x: self.amplitude_changed.emit(float(x)))
        row.addWidget(self.amp_slider, 1)
        self.amp_value = QLabel("18 Å")
        self.amp_slider.valueChanged.connect(
            lambda x: self.amp_value.setText(f"{x} Å"))
        row.addWidget(self.amp_value)
        v.addLayout(row)

        row = QHBoxLayout()
        self.animate_button = QPushButton("Animate")
        self.animate_button.setCheckable(True)
        self.animate_button.setEnabled(False)
        self.animate_button.toggled.connect(self._on_animate)
        row.addWidget(self.animate_button)

        self.color_button = QPushButton("Colour by displacement")
        self.color_button.setCheckable(True)
        self.color_button.setEnabled(False)
        self.color_button.toggled.connect(self.color_by_mode_requested.emit)
        row.addWidget(self.color_button)
        v.addLayout(row)
        layout.addWidget(box)

        layout.addStretch(1)
        self._modes = None

    # ------------------------------------------------------------- emitters

    def _emit_compute(self) -> None:
        self.compute_modes_requested.emit({
            "cutoff": self.cutoff_spin.value(),
            "spring": self.spring_combo.currentText(),
            "n_modes": self.nmodes_spin.value(),
        })

    def _on_mode(self, index: int) -> None:
        if self._modes is None or index < 0:
            return
        self.mode_selected.emit(index)
        self._describe(index)

    def _on_animate(self, on: bool) -> None:
        self.animate_button.setText("Stop" if on else "Animate")
        self.animate_toggled.emit(on)

    # -------------------------------------------------------------- updates

    def set_busy(self, busy: bool) -> None:
        self.progress.setVisible(busy)
        self.compute_button.setEnabled(not busy)
        self.measure_button.setEnabled(not busy)

    def set_dome(self, dome, reference: str = "") -> None:
        if dome is None:
            self.dome_label.setText("Could not measure — need three protomers.")
            return
        rc = dome.radius_of_curvature / 10.0
        self.dome_label.setText(
            f"<b>radius of curvature {rc:.1f} nm</b><br>"
            f"dome depth {dome.dome_depth / 10:.1f} nm · "
            f"footprint radius {dome.footprint_radius / 10:.1f} nm<br>"
            f"dome area {dome.dome_area / 100:.0f} nm² · "
            f"projected {dome.projected_area / 100:.0f} nm² · "
            f"excess {dome.excess_area / 100:.0f} nm²<br>"
            f"<span style='color:#7f8798'>C3 axis recovered at "
            f"{dome.notes.get('c3_angle_deg', 0):.2f}° "
            f"(RMSD {dome.notes.get('c3_rmsd', 0):.2f} Å); "
            f"sphere fit RMSE {dome.notes.get('sphere_rmse', 0):.1f} Å"
            f"{reference}</span>")

    def set_modes(self, modes) -> None:
        self._modes = modes
        self.mode_combo.blockSignals(True)
        self.mode_combo.clear()
        if modes is None:
            self.mode_combo.setEnabled(False)
            self.animate_button.setEnabled(False)
            self.color_button.setEnabled(False)
            self.mode_combo.blockSignals(False)
            return
        for i in range(modes.n_modes):
            sym = modes.symmetry[i] if modes.symmetry is not None else "?"
            self.mode_combo.addItem(
                f"Mode {i + 1}   {sym}   λ={modes.eigenvalues[i]:.5f}")
        self.mode_combo.blockSignals(False)
        self.mode_combo.setEnabled(True)
        self.animate_button.setEnabled(True)
        self.color_button.setEnabled(True)
        self.mode_combo.setCurrentIndex(0)
        self._describe(0)

    def _describe(self, index: int) -> None:
        m = self._modes
        if m is None:
            return
        sym = m.symmetry[index] if m.symmetry is not None else "?"
        char = m.character[index] if m.character is not None else float("nan")
        coll = m.collectivity(index)
        if sym == "A":
            meaning = ("three-fold symmetric — <b>can</b> couple to isotropic "
                       "membrane tension, so this is a candidate gating "
                       "coordinate")
        else:
            meaning = ("degenerate E pair — symmetry forbids first-order "
                       "coupling to isotropic tension")
        self.mode_label.setText(
            f"symmetry <b>{sym}</b> (character {char:+.3f}) · "
            f"collectivity {coll:.3f}<br>{meaning}")
