"""Structure chooser and appearance controls."""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QFormLayout, QGroupBox,
                             QHBoxLayout, QLabel, QSlider, QVBoxLayout,
                             QWidget)

from ...io.registry import StructureRecord, load_registry
from ...structure.full_length import FILL_MODES
from ...render.representations import ColorBy, Style

__all__ = ["StructurePanel", "FILTERS"]

#: The catalogue is 21 entries and one combo of 21 is a list to scroll, not a
#: choice to make. Each filter is a *field of the record*, and its options are
#: read from the data rather than written here — a new state or a second PIEZO2
#: entry appears in the box without anyone remembering to add it.
#:
#: ``protein`` is first because it is the one that changes what you are looking
#: at rather than which view of it: PIEZO2 is a different molecule, and it was
#: reachable only by knowing that 6KG7 is the one mouse entry that is not
#: PIEZO1.
FILTERS = (
    ("protein", "Protein", "PIEZO1 or the PIEZO2 paralogue (6KG7). Measured\n"
                           "from each file's own residue names, not curated."),
    ("species", "Species", "Which numbering the deposited residue numbers use.\n"
                           "The human↔mouse offset is not constant."),
    ("state", "State", "Curved, flattened, intermediate — or a fragment of\n"
                       "the molecule rather than the whole trimer."),
    ("gating", "Gating", "What the pore is doing. Every deposited human entry\n"
                         "is closed; 11ZC is the only open-like one."),
)

STYLE_LABELS = [
    ("Cartoon", Style.CARTOON),
    ("Tube", Style.TUBE),
    ("Backbone trace", Style.BACKBONE),
    ("Spheres (van der Waals)", Style.SPHERES),
    ("Balls", Style.BALLS),
    ("Sticks", Style.STICKS),
    ("Ball and stick", Style.BALL_AND_STICK),
]

COLOR_LABELS = [
    ("Domain", ColorBy.DOMAIN),
    ("Chain / protomer", ColorBy.CHAIN),
    ("Secondary structure", ColorBy.SECONDARY),
    ("B-factor", ColorBy.BFACTOR),
    ("AlphaFold pLDDT", ColorBy.PLDDT),
    ("Hydrophobicity (Kyte-Doolittle)", ColorBy.HYDROPHOBICITY),
    ("Element", ColorBy.ELEMENT),
    ("Computed value", ColorBy.VALUE),
    ("Uniform", ColorBy.UNIFORM),
]


class StructurePanel(QWidget):
    """Pick a structure and control how it is drawn."""

    structure_requested = pyqtSignal(str)
    fill_changed = pyqtSignal(str)              # which prediction to splice in
    style_changed = pyqtSignal(object)
    color_changed = pyqtSignal(object)
    ligands_toggled = pyqtSignal(bool)
    radius_changed = pyqtSignal(float)
    entities_changed = pyqtSignal(object)      # frozenset of visible classes
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
        # Two filters per row: four rows of one combo each would push the
        # entry chooser off the bottom of a normally-sized dock.
        self.filter_combos: dict[str, QComboBox] = {}
        for left, right in zip(FILTERS[::2], FILTERS[1::2]):
            row = QWidget()
            line = QHBoxLayout(row)
            line.setContentsMargins(0, 0, 0, 0)
            line.setSpacing(6)
            labels = []
            for field, title, tip in (left, right):
                combo = QComboBox()
                combo.addItems(self._options(field))
                combo.setToolTip(tip)
                combo.currentTextChanged.connect(self._repopulate)
                self.filter_combos[field] = combo
                line.addWidget(combo, 1)
                labels.append(title)
            form.addRow(" / ".join(labels), row)

        #: Backwards-compatible name: the species filter was the only one.
        self.species_combo = self.filter_combos["species"]

        self.structure_combo = QComboBox()
        self.structure_combo.setMinimumWidth(260)
        self.structure_combo.setToolTip(
            "The curated PIEZO catalogue, narrowed by the filters above.\n"
            "Each entry states its gating state, resolved residue range,\n"
            "numbering species and citation.")
        self.structure_combo.currentIndexChanged.connect(self._on_structure)
        form.addRow("Entry", self.structure_combo)

        self.count = QLabel("")
        self.count.setStyleSheet("color: #7f8798; font-size: 11px;")
        form.addRow("", self.count)

        # Which model to build from the chosen entry. It sits with the entry
        # rather than under View because it decides *what is loaded*, not how
        # it is drawn — every analysis, animation and measurement then runs on
        # whatever this says, and none of them has to know about it.
        self.fill_combo = QComboBox()
        for _key, label, _tip in FILL_MODES:
            self.fill_combo.addItem(label)
        self.fill_combo.setToolTip(
            "\n\n".join(f"{label}: {tip}" for _k, label, tip in FILL_MODES))
        self.fill_combo.currentIndexChanged.connect(
            lambda _i: self.fill_changed.emit(self.current_fill()))
        form.addRow("Completeness", self.fill_combo)

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
        self.style_combo.setToolTip(
            "Cartoon uses secondary structure assigned from C-alpha geometry\n"
            "(P-SEA), because most of these entries have no backbone atoms.\n"
            "Spheres and sticks are ray-cast impostors, so they are exact.")
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
        self.radius_slider.setToolTip("Scale every atom radius; visual only")
        form.addRow("Atom size", self.radius_slider)

        self.entity_box = QGroupBox("What is in this file")
        self.entity_box.setToolTip(
            "A deposited PIEZO1 file contains more than the channel: lipids and\n"
            "detergent from the sample, glycan, and in six entries the MDFIC\n"
            "auxiliary subunit. Each can be shown or hidden independently.")
        self.entity_layout = QVBoxLayout(self.entity_box)
        self.entity_summary = QLabel("no structure loaded")
        self.entity_summary.setWordWrap(True)
        self.entity_summary.setStyleSheet("color:#8a919e;")
        self.entity_layout.addWidget(self.entity_summary)
        self.entity_checks: dict = {}
        form.addRow(self.entity_box)

        # What is on screen, which is not the same question as what is in the
        # file. Hidden while only one structure is drawn, so the normal case
        # carries no extra furniture.
        self.displayed_label = QLabel()
        self.displayed_label.setWordWrap(True)
        self.displayed_label.setToolTip(
            "Structures currently drawn. The first is the primary one; every\n"
            "analysis runs on it, whatever else is displayed. Turn extra\n"
            "structures on and off under View.")
        self.displayed_label.setVisible(False)
        form.addRow(self.displayed_label)

        self.ligand_check = QCheckBox("Show lipids and ligands")
        self.ligand_check.setChecked(True)
        self.ligand_check.setToolTip(
            "Show resolved lipids and ligands as spheres. These are modelled\n"
            "densities, not docked poses.")
        self.ligand_check.toggled.connect(self.ligands_toggled.emit)
        form.addRow(self.ligand_check)

        self.spin_check = QCheckBox("Auto-rotate")
        self.spin_check.setToolTip(
            "Rotate continuously. Speed is set in Options.")
        self.spin_check.toggled.connect(self.spin_toggled.emit)
        form.addRow(self.spin_check)
        layout.addWidget(box)

        layout.addStretch(1)
        self._repopulate()

    # ------------------------------------------------------------- helpers

    def _options(self, field: str) -> list[str]:
        """The values this field actually takes, so the box cannot go stale."""
        values = {str(getattr(r, field, "")) for r in self.registry.available()}
        return ["all"] + sorted(v for v in values if v)

    def _visible_records(self) -> list[StructureRecord]:
        records = self.registry.available()
        for field, combo in self.filter_combos.items():
            wanted = combo.currentText()
            if wanted != "all":
                records = [r for r in records
                           if str(getattr(r, field, "")) == wanted]
        return records

    def clear_filters(self) -> None:
        """Put every filter back to "all", without rebuilding four times."""
        for combo in self.filter_combos.values():
            combo.blockSignals(True)
            combo.setCurrentIndex(0)
            combo.blockSignals(False)
        self._repopulate()

    def _repopulate(self) -> None:
        self.structure_combo.blockSignals(True)
        self.structure_combo.clear()
        self._records = self._visible_records()
        for r in self._records:
            self.structure_combo.addItem(r.label(), r.pdb)
        self.structure_combo.blockSignals(False)

        total = len(self.registry.available())
        self.count.setText(f"{len(self._records)} of {total} downloaded entries")
        if self._records:
            self._on_structure(0)
        else:
            # Leaving the previous structure's details on screen under an empty
            # chooser reads as "this is what you are looking at", and it is not.
            self.detail.setText(
                "<b>No entry matches these filters.</b><br>"
                "The catalogue has one PIEZO2 entry and one fragment, so some "
                "combinations are empty by construction.")

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
        """Show ``pdb``, clearing the filters if they are hiding it.

        A silent no-op here would be a trap: ``--structure 6KG7``, a restored
        session and the morph controller all arrive through this, and with a
        filter set they would leave whatever was already loaded on screen
        while appearing to have honoured the request.
        """
        pdb = pdb.upper()
        if not any(r.pdb == pdb for r in self._records):
            if any(r.pdb == pdb for r in self.registry.available()):
                self.clear_filters()
        for i, r in enumerate(self._records):
            if r.pdb == pdb:
                self.structure_combo.setCurrentIndex(i)
                return

    def current_fill(self) -> str:
        return FILL_MODES[max(self.fill_combo.currentIndex(), 0)][0]

    def set_fill(self, mode: str) -> None:
        """Set the completeness without asking for a reload."""
        for index, (key, _label, _tip) in enumerate(FILL_MODES):
            if key == mode:
                self.fill_combo.blockSignals(True)
                self.fill_combo.setCurrentIndex(index)
                self.fill_combo.blockSignals(False)
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

    def set_displayed(self, names: list[str], colors: dict | None = None) -> None:
        """Show which structures are on screen, and in which colour.

        Only appears when there is more than one, since with a single structure
        the answer is already in the chooser above.
        """
        colors = colors or {}
        if len(names) < 2:
            self.displayed_label.setVisible(False)
            return

        parts = [f"<b>{names[0]}</b> (primary)"]
        for name in names[1:]:
            rgb = colors.get(name)
            swatch = ""
            if rgb is not None:
                hexcode = "#%02x%02x%02x" % tuple(int(255 * c) for c in rgb)
                swatch = (f"<span style='color:{hexcode};'>&#9632;</span> ")
            parts.append(f"{swatch}{name}")
        self.displayed_label.setText("Displayed: " + ", ".join(parts))
        self.displayed_label.setVisible(True)

    def set_entities(self, entity_map) -> None:
        """Rebuild the per-category checkboxes for the loaded structure.

        Only categories actually present get a control. A permanent list with
        most entries greyed out would say nothing about what you are looking
        at, and what is in the file is exactly the thing worth surfacing.
        """
        from PyQt6.QtWidgets import QCheckBox

        from ...core.entities import EntityClass

        for widget in self.entity_checks.values():
            self.entity_layout.removeWidget(widget)
            widget.deleteLater()
        self.entity_checks = {}

        if entity_map is None:
            self.entity_summary.setText("no structure loaded")
            return

        counts = entity_map.counts()
        self.entity_summary.setText(entity_map.summary())
        for key in entity_map.present():
            label = f"{EntityClass.LABELS[key]} ({counts[key]:,} atoms)"
            check = QCheckBox(label)
            check.setChecked(True)
            check.toggled.connect(lambda _on: self._emit_entities())
            self.entity_layout.addWidget(check)
            self.entity_checks[key] = check
        self._emit_entities()

    def _emit_entities(self) -> None:
        self.entities_changed.emit(
            frozenset(k for k, w in self.entity_checks.items() if w.isChecked()))

    def visible_entities(self) -> frozenset:
        return frozenset(k for k, w in self.entity_checks.items()
                         if w.isChecked())
