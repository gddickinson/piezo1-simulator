"""The topology window: Figure 3 for whichever structure is loaded.

A separate top-level window rather than a dock, for the same reason the
sequence window is one — 38 helices in a row wants the full width, and it is a
thing you switch to rather than glance at.

What it adds over the drawing itself:

* **Unit selection.** A checkbox per 4-TM unit, plus the pore module. Ticking
  one boxes it in the diagram exactly as Figure 3b's red boxes do, and selects
  its residues on the 3-D model, so the box is a question you can follow up
  rather than an annotation.
* **A resolved/complete switch.** The diagram can be drawn against what the
  entry models, or against the whole 38-helix architecture with the unmodelled
  parts dashed. The second is Figure 3a; the first is what you have.
* **PNG export**, because a topology diagram is something people put in talks.

The window states which numbering the residue ranges are in, on the diagram and
beside the export. It is the mouse table for a mouse entry and the human one
for a human entry, and a range copied out of here into a paper written in the
other convention would be wrong by up to 26 residues.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QImage, QPainter
from PyQt6.QtWidgets import (QCheckBox, QFileDialog, QGroupBox, QHBoxLayout,
                             QLabel, QMainWindow, QPushButton, QScrollArea,
                             QVBoxLayout, QWidget)

from .topology_view import TopologyView

__all__ = ["TopologyWindow"]

#: Element kinds the user can turn on and off, with the label each gets.
KIND_TOGGLES = (("tm_helix", "Transmembrane helices"),
                ("loop", "Loops"),
                ("box", "Cap (CED)"),
                ("bar", "Beam"),
                ("cuff", "Cuff (elbow, base, hairpin, PE)"))


class TopologyWindow(QMainWindow):
    """Monomer topology in the membrane, with selectable 4-TM units."""

    residues_selected = pyqtSignal(int, int, str)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("PIEZO1 — monomer topology")
        self.resize(1180, 520)
        self._unit_boxes: dict[int, QCheckBox] = {}

        self.view = TopologyView()
        self.view.residues_selected.connect(self.residues_selected)
        self.view.status.connect(self._set_status)

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addWidget(self._controls())
        layout.addWidget(self.view, 1)
        self.caption = QLabel()
        self.caption.setWordWrap(True)
        self.caption.setStyleSheet("color: #8a919e;")
        layout.addWidget(self.caption)
        self.setCentralWidget(central)
        self.statusBar().showMessage("Load a structure to draw its topology")
        self._update_caption()

    # ------------------------------------------------------------- controls

    def _controls(self) -> QWidget:
        bar = QWidget()
        row = QHBoxLayout(bar)
        row.setContentsMargins(0, 0, 0, 0)

        self.units_box = QGroupBox("4-TM units to box (Figure 3b)")
        self.units_layout = QHBoxLayout(self.units_box)
        self.units_layout.setContentsMargins(6, 2, 6, 2)
        scroll = QScrollArea()
        scroll.setWidget(self.units_box)
        scroll.setWidgetResizable(True)
        scroll.setMaximumHeight(64)
        scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        row.addWidget(scroll, 3)

        kinds = QGroupBox("Show")
        kinds_layout = QVBoxLayout(kinds)
        kinds_layout.setContentsMargins(6, 2, 6, 2)
        kinds_layout.setSpacing(1)
        self._kind_boxes: dict[str, QCheckBox] = {}
        for kind, label in KIND_TOGGLES:
            box = QCheckBox(label)
            box.setChecked(True)
            box.toggled.connect(self._kinds_changed)
            kinds_layout.addWidget(box)
            self._kind_boxes[kind] = box
        row.addWidget(kinds, 1)

        buttons = QWidget()
        column = QVBoxLayout(buttons)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(3)
        for text, slot in (("Select all units", self.select_all_units),
                           ("Clear", self.clear_units),
                           ("Export PNG…", self.export_png)):
            button = QPushButton(text)
            button.clicked.connect(slot)
            column.addWidget(button)
        column.addStretch(1)
        row.addWidget(buttons)
        return bar

    def _kinds_changed(self) -> None:
        self.view.shown_kinds = {kind for kind, box in self._kind_boxes.items()
                                 if box.isChecked()}
        self.view.update()

    # -------------------------------------------------------------- content

    def set_structure(self, structure, reference: str = "mouse") -> None:
        """Draw the topology of a structure, in its own numbering."""
        self.view.set_structure(structure, reference)
        self._rebuild_units()
        self._update_caption()
        topology = self.view.topology
        if topology is not None:
            self._set_status(topology.summary())

    def _rebuild_units(self) -> None:
        while self.units_layout.count():
            item = self.units_layout.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self._unit_boxes.clear()
        for unit in self.view.available_units():
            box = QCheckBox(f"THU{unit}")
            box.toggled.connect(self._units_changed)
            self.units_layout.addWidget(box)
            self._unit_boxes[unit] = box
        self.units_layout.addStretch(1)

    def _units_changed(self) -> None:
        chosen = {unit for unit, box in self._unit_boxes.items()
                  if box.isChecked()}
        self.view.set_boxed(chosen)
        self._update_caption()
        if chosen:
            lo, hi, numbering = self.view.boxed_range()
            gapped = sorted(chosen) != list(range(min(chosen), max(chosen) + 1))
            self._set_status(
                f"THU{', '.join(str(u) for u in sorted(chosen))}: "
                f"{lo}-{hi} ({numbering} numbering)"
                + (" — selection is not contiguous, so the model highlights "
                   "the whole span between the first and last" if gapped
                   else ""))

    def select_all_units(self) -> None:
        for box in self._unit_boxes.values():
            box.setChecked(True)

    def clear_units(self) -> None:
        for box in self._unit_boxes.values():
            box.setChecked(False)

    def _update_caption(self) -> None:
        topology = self.view.topology
        if topology is None:
            self.caption.setText(
                "Load a structure from the Model panel; the diagram follows it.")
            return
        unresolved = topology.unresolved_helices
        self.caption.setText(
            f"Nine 4-TM units plus the pore module (TM37-38), after "
            f"Guo & MacKinnon 2017 Figure 3. Residue ranges are in "
            f"**{topology.numbering}** numbering. "
            + (f"{len(unresolved)} helices are drawn dashed because "
               f"{topology.structure} does not model them — they are kept in "
               f"the diagram so the helix numbering is not silently shifted. "
               if unresolved else "")
            + "Shift-click a helix to box its unit; click any element to "
              "select it on the model.")

    def _set_status(self, message: str) -> None:
        self.statusBar().showMessage(message)

    # --------------------------------------------------------------- export

    def export_png(self, checked: bool = False,
                   path: str | None = None) -> str | None:
        """Render the diagram to a PNG at twice the on-screen size."""
        if path is None:
            name = getattr(self.view.topology, "structure", None) or "topology"
            path, _ = QFileDialog.getSaveFileName(
                self, "Export topology", f"{name}_topology.png",
                "PNG image (*.png)")
        if not path:
            return None
        scale = 2
        image = QImage(self.view.width() * scale, self.view.height() * scale,
                       QImage.Format.Format_ARGB32)
        image.fill(0)
        painter = QPainter(image)
        painter.scale(scale, scale)
        self.view.render(painter)
        painter.end()
        image.save(path)
        self._set_status(f"wrote {path}")
        return path
