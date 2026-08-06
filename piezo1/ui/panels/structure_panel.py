"""Structure chooser and appearance controls."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QGroupBox,
                             QLabel, QSlider, QVBoxLayout, QWidget)

from ...io.registry import StructureRecord, load_registry
from ...render.representations import ColorBy, Style

__all__ = ["StructurePanel"]

STYLE_LABELS = [
    ("Cartoon", Style.CARTOON),
    ("Tube", Style.TUBE),
    ("Backbone trace", Style.BACKBONE),
    ("Spheres (van der Waals)", Style.SPHERES),
    ("Ball and stick", Style.BALL_AND_STICK),
]

COLOR_LABELS = [
    ("Domain", ColorBy.DOMAIN),
    ("Chain / protomer", ColorBy.CHAIN),
    ("Secondary structure", ColorBy.SECONDARY),
    ("B-factor", ColorBy.BFACTOR),
    ("AlphaFold pLDDT", ColorBy.PLDDT),
    ("Element", ColorBy.ELEMENT),
    ("Computed value", ColorBy.VALUE),
    ("Uniform", ColorBy.UNIFORM),
]


class StructurePanel(QWidget):
    """Pick a structure and control how it is drawn."""

    structure_requested = pyqtSignal(str)
    style_changed = pyqtSignal(object)
    color_changed = pyqtSignal(object)
    ligands_toggled = pyqtSignal(bool)
    radius_changed = pyqtSignal(float)
    spin_toggled = pyqtSignal(bool)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.registry = load_registry()
        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.setSpacing(10)

        # ---------------------------------------------------------- structure
        box = QGroupBox("Structure")
        form = QFormLayout(box)
        self.species_combo = QComboBox()
        self.species_combo.addItems(["all", "human", "mouse"])
        self.species_combo.currentTextChanged.connect(self._repopulate)
        form.addRow("Species", self.species_combo)

        self.structure_combo = QComboBox()
        self.structure_combo.setMinimumWidth(260)
        self.structure_combo.currentIndexChanged.connect(self._on_structure)
        form.addRow("Entry", self.structure_combo)

        self.detail = QLabel("—")
        self.detail.setWordWrap(True)
        self.detail.setTextFormat(Qt.TextFormat.RichText)
        self.detail.setOpenExternalLinks(True)
        self.detail.setStyleSheet("color: #9aa3b2; font-size: 11px;")
        form.addRow(self.detail)
        layout.addWidget(box)

        # -------------------------------------------------------- appearance
        box = QGroupBox("Appearance")
        form = QFormLayout(box)
        self.style_combo = QComboBox()
        for label, _ in STYLE_LABELS:
            self.style_combo.addItem(label)
        self.style_combo.currentIndexChanged.connect(
            lambda i: self.style_changed.emit(STYLE_LABELS[i][1]))
        form.addRow("Representation", self.style_combo)

        self.color_combo = QComboBox()
        for label, _ in COLOR_LABELS:
            self.color_combo.addItem(label)
        self.color_combo.currentIndexChanged.connect(
            lambda i: self.color_changed.emit(COLOR_LABELS[i][1]))
        form.addRow("Colour by", self.color_combo)

        self.radius_slider = QSlider(Qt.Orientation.Horizontal)
        self.radius_slider.setRange(20, 250)
        self.radius_slider.setValue(100)
        self.radius_slider.valueChanged.connect(
            lambda v: self.radius_changed.emit(v / 100.0))
        form.addRow("Atom size", self.radius_slider)

        self.ligand_check = QCheckBox("Show lipids and ligands")
        self.ligand_check.setChecked(True)
        self.ligand_check.toggled.connect(self.ligands_toggled.emit)
        form.addRow(self.ligand_check)

        self.spin_check = QCheckBox("Auto-rotate")
        self.spin_check.toggled.connect(self.spin_toggled.emit)
        form.addRow(self.spin_check)
        layout.addWidget(box)

        layout.addStretch(1)
        self._repopulate()

    # ------------------------------------------------------------- helpers

    def _visible_records(self) -> list[StructureRecord]:
        species = self.species_combo.currentText()
        records = self.registry.available()
        if species != "all":
            records = [r for r in records if r.species == species]
        return records

    def _repopulate(self) -> None:
        self.structure_combo.blockSignals(True)
        self.structure_combo.clear()
        self._records = self._visible_records()
        for r in self._records:
            self.structure_combo.addItem(r.label(), r.pdb)
        self.structure_combo.blockSignals(False)
        if self._records:
            self._on_structure(0)

    def _on_structure(self, index: int) -> None:
        if not (0 <= index < len(self._records)):
            return
        rec = self._records[index]
        rng = rec.residue_range
        span = f"{rng[0]}–{rng[1]}" if rng else "n/a"
        doi = (f' · <a style="color:#7aa7ff" href="https://doi.org/{rec.doi}">DOI</a>'
               if rec.doi else "")
        lig = ", ".join(rec.ligands) if rec.ligands else "none"
        self.detail.setText(
            f"<b>{rec.state}</b> · {rec.gating} · {rec.n_protomers} protomers<br>"
            f"residues {span} ({rec.numbering_species} numbering)<br>"
            f"ligands: {lig}<br>"
            f"<i>{rec.note}</i><br>{rec.citation()}{doi}"
        )
        self.structure_requested.emit(rec.pdb)

    def select(self, pdb: str) -> None:
        for i, r in enumerate(self._records):
            if r.pdb == pdb.upper():
                self.structure_combo.setCurrentIndex(i)
                return

    def current_record(self) -> StructureRecord | None:
        i = self.structure_combo.currentIndex()
        return self._records[i] if 0 <= i < len(self._records) else None

    def set_state(self, style: str | None = None, color_by: str | None = None,
                  ligands: bool | None = None,
                  radius_scale: float | None = None) -> None:
        """Restore appearance controls from a saved session.

        Signals are blocked while the widgets are set and re-emitted once at
        the end, so restoring four settings triggers one rebuild of the scene
        rather than four — the intermediate states are not views anybody asked
        for, and on a 120k-atom trimer they are visible as flicker.
        """
        for widget in (self.style_combo, self.color_combo, self.ligand_check,
                       self.radius_slider):
            widget.blockSignals(True)
        try:
            if style is not None:
                for i, (_, value) in enumerate(STYLE_LABELS):
                    if value.value == style:
                        self.style_combo.setCurrentIndex(i)
                        break
            if color_by is not None:
                for i, (_, value) in enumerate(COLOR_LABELS):
                    if value.value == color_by:
                        self.color_combo.setCurrentIndex(i)
                        break
            if ligands is not None:
                self.ligand_check.setChecked(bool(ligands))
            if radius_scale is not None:
                self.radius_slider.setValue(int(round(radius_scale * 100)))
        finally:
            for widget in (self.style_combo, self.color_combo,
                           self.ligand_check, self.radius_slider):
                widget.blockSignals(False)

        self.style_changed.emit(self.current_style())
        self.color_changed.emit(self.current_color())
        self.ligands_toggled.emit(self.ligand_check.isChecked())
        self.radius_changed.emit(self.radius_slider.value() / 100.0)

    def current_style(self):
        return STYLE_LABELS[self.style_combo.currentIndex()][1]

    def current_color(self):
        return COLOR_LABELS[self.color_combo.currentIndex()][1]

    @property
    def radius_spin(self):
        """Alias so callers can read the atom-size control by one name."""
        return self.radius_slider

