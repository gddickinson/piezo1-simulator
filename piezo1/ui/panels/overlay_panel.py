"""Controls for overlaying a second structure on the loaded one."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QPushButton, QVBoxLayout,
                             QWidget)

from ...render.representations import Style

__all__ = ["OverlayPanel"]

#: Representations worth using for a second structure. Cartoon and backbone
#: read clearly against a reference; spheres would bury it.
OVERLAY_STYLES = [("Backbone trace", Style.BACKBONE),
                  ("Tube", Style.TUBE),
                  ("Cartoon", Style.CARTOON)]

SUPERPOSITION_MODES = [
    ("Match protomers, then fit", "protomer"),
    ("Fit on chain A only", "chain"),
]


class OverlayPanel(QWidget):
    """Choose a second structure, superpose it, and control how it is drawn."""

    overlay_requested = pyqtSignal(str, str)      # pdb, mode
    clear_requested = pyqtSignal()
    style_changed = pyqtSignal(object)
    visibility_changed = pyqtSignal(bool)
    deviation_colour_requested = pyqtSignal(bool)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        box = QVBoxLayout(self)

        chooser = QGroupBox("Second structure")
        form = QFormLayout(chooser)
        self.structure_combo = QComboBox()
        self.structure_combo.setToolTip(
            "Superposed onto the loaded structure by shared residue numbers.\n"
            "Cross-species pairs are refused: human and mouse numbering do not\n"
            "correspond, so the fit would be confidently wrong.")
        form.addRow("Overlay", self.structure_combo)

        self.mode_combo = QComboBox()
        for label, key in SUPERPOSITION_MODES:
            self.mode_combo.addItem(label, key)
        self.mode_combo.setToolTip(
            "Deposited entries do not label protomers in a consistent\n"
            "rotational order. Matching them first puts 7WLU onto 7WLT at\n"
            "12.3 Å instead of 90.7 Å.")
        form.addRow("Superposition", self.mode_combo)
        box.addWidget(chooser)

        row = QHBoxLayout()
        self.load_button = QPushButton("Overlay")
        self.load_button.clicked.connect(
            lambda: self.overlay_requested.emit(
                self.structure_combo.currentData() or "",
                self.mode_combo.currentData()))
        row.addWidget(self.load_button)
        self.clear_button = QPushButton("Remove")
        self.clear_button.clicked.connect(self.clear_requested.emit)
        row.addWidget(self.clear_button)
        box.addLayout(row)

        display = QGroupBox("Display")
        inner = QVBoxLayout(display)
        self.style_combo = QComboBox()
        for label, style in OVERLAY_STYLES:
            self.style_combo.addItem(label, style)
        self.style_combo.currentIndexChanged.connect(
            lambda _: self.style_changed.emit(self.style_combo.currentData()))
        inner.addWidget(self.style_combo)

        self.visible_check = QCheckBox("Show overlay")
        self.visible_check.setChecked(True)
        self.visible_check.toggled.connect(self.visibility_changed.emit)
        inner.addWidget(self.visible_check)

        self.deviation_check = QCheckBox("Colour reference by deviation")
        self.deviation_check.setToolTip(
            "Colour the loaded structure by how far the overlay sits from it\n"
            "at each residue — where the two actually differ.")
        self.deviation_check.toggled.connect(
            self.deviation_colour_requested.emit)
        inner.addWidget(self.deviation_check)
        box.addWidget(display)

        self.result_label = QLabel("no overlay")
        self.result_label.setWordWrap(True)
        self.result_label.setTextFormat(Qt.TextFormat.RichText)
        box.addWidget(self.result_label)
        box.addStretch(1)

    def set_choices(self, records, exclude: str = "") -> None:
        self.structure_combo.clear()
        for record in records:
            if record.pdb == exclude:
                continue
            self.structure_combo.addItem(
                f"{record.pdb} — {record.state}", record.pdb)

    def set_result(self, result, error: str = "") -> None:
        if error:
            self.result_label.setText(f"<span style='color:#f26d6d'>{error}</span>")
            return
        if result is None:
            self.result_label.setText("no overlay")
            return
        text = (f"<b>RMSD {result.rmsd:.2f} Å</b> over {result.n_common} "
                f"common C-alphas<br>max deviation "
                f"{result.meta.get('max_deviation', 0.0):.1f} Å")
        if result.reordered:
            text += (f"<br><span style='color:#f2a65a'>protomers rematched "
                     f"{result.protomer_order}</span> — by chain label this "
                     f"would have been {result.rmsd_by_label:.1f} Å")
        self.result_label.setText(text)
