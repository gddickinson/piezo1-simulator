"""The sequence window: browse, annotate, select onto the model, and compare.

A separate top-level window rather than a dock, because sequence work wants the
full width — 2521 residues at 60 columns is 42 rows — and because it is a task
you switch to rather than glance at.

Selection is two-way: dragging across the sequence highlights those residues in
the 3-D model, and anything selected elsewhere in the application shows up here.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QAbstractItemView, QCheckBox, QComboBox,
                             QHBoxLayout, QLabel, QLineEdit, QMainWindow,
                             QPushButton, QScrollArea, QSplitter, QTableWidget,
                             QTableWidgetItem, QTabWidget, QVBoxLayout, QWidget)

from ..core.annotations import load_annotations
from ..core.sequences import compare_sequences, load_named_sequences
from .sequence_view import ResidueStyle, SequenceView

__all__ = ["SequenceWindow"]

#: Colouring schemes for the sequence, beyond the default chemistry palette.
COLOUR_MODES = [("Chemistry", "chemistry"), ("Domain", "domain"),
                ("Variants", "variants"), ("Functional sites", "sites"),
                ("Plain", "plain")]


class SequenceWindow(QMainWindow):
    """Browse and compare the sequences behind the model."""

    residues_selected = pyqtSignal(object, str)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("PIEZO1 — sequences")
        self.resize(1100, 720)
        self.sequences: list = []
        self.by_key: dict = {}
        self.annotations = load_annotations("human")
        self._comparison = None

        tabs = QTabWidget()
        tabs.addTab(self._build_browser(), "Browse")
        tabs.addTab(self._build_compare(), "Compare")
        self.setCentralWidget(tabs)
        self.statusBar().showMessage("no structure loaded")

    # --------------------------------------------------------------- browser

    def _build_browser(self) -> QWidget:
        page = QWidget()
        box = QVBoxLayout(page)

        row = QHBoxLayout()
        self.sequence_combo = QComboBox()
        self.sequence_combo.setToolTip(
            "UniProt is the reference numbering variants are quoted in.\n"
            "A structure sequence contains only resolved residues and has gaps.\n"
            "A translated CDS carries the real codons.")
        self.sequence_combo.currentIndexChanged.connect(self._refresh)
        row.addWidget(QLabel("Sequence"), 0)
        row.addWidget(self.sequence_combo, 2)

        self.colour_combo = QComboBox()
        for label, key in COLOUR_MODES:
            self.colour_combo.addItem(label, key)
        self.colour_combo.setCurrentIndex(1)
        self.colour_combo.currentIndexChanged.connect(self._refresh)
        row.addWidget(QLabel("Colour by"), 0)
        row.addWidget(self.colour_combo, 1)

        self.dna_check = QCheckBox("Show codons")
        self.dna_check.setToolTip(
            "Real coding sequence from the Ensembl canonical transcript,\n"
            "not a back-translation. Only available for sequences that have one.")
        self.dna_check.toggled.connect(self._refresh)
        row.addWidget(self.dna_check)
        box.addLayout(row)

        find = QHBoxLayout()
        self.find_edit = QLineEdit()
        self.find_edit.setPlaceholderText(
            "Go to residue number, or find a motif such as KKKK…")
        self.find_edit.returnPressed.connect(self._find)
        find.addWidget(self.find_edit, 1)
        button = QPushButton("Go")
        button.clicked.connect(self._find)
        find.addWidget(button)
        self.select_button = QPushButton("Show selection on the model")
        self.select_button.setToolTip(
            "Highlight the selected residues in the 3-D view")
        self.select_button.clicked.connect(self._emit_selection)
        find.addWidget(self.select_button)
        box.addLayout(find)

        self.view = SequenceView()
        self.view.selection_changed.connect(self._on_selection)
        self.view.residue_hovered.connect(self._on_hover)
        area = QScrollArea()
        area.setWidget(self.view)
        area.setWidgetResizable(True)
        box.addWidget(area, 1)

        self.detail = QLabel("")
        self.detail.setWordWrap(True)
        self.detail.setTextFormat(Qt.TextFormat.RichText)
        self.detail.setMinimumHeight(56)
        box.addWidget(self.detail)
        return page

    # --------------------------------------------------------------- compare

    def _build_compare(self) -> QWidget:
        page = QWidget()
        box = QVBoxLayout(page)

        row = QHBoxLayout()
        self.compare_a = QComboBox()
        self.compare_b = QComboBox()
        self.method_combo = QComboBox()
        self.method_combo.addItem("Global alignment (Needleman–Wunsch)", "global")
        self.method_combo.addItem("By residue number (no gaps)", "positional")
        self.method_combo.setToolTip(
            "Global alignment inserts gaps to maximise score, which is right\n"
            "across species. By residue number pairs positions directly, which\n"
            "is the honest choice when both sequences already share a numbering:\n"
            "an aligner can slide residues to buy score and invent differences.")
        for widget, label in ((self.compare_a, "A"), (self.compare_b, "B")):
            row.addWidget(QLabel(label))
            row.addWidget(widget, 2)
        row.addWidget(self.method_combo, 2)
        run = QPushButton("Compare")
        run.clicked.connect(self._run_comparison)
        row.addWidget(run)
        box.addLayout(row)

        self.compare_summary = QLabel("no comparison run")
        self.compare_summary.setWordWrap(True)
        self.compare_summary.setTextFormat(Qt.TextFormat.RichText)
        box.addWidget(self.compare_summary)

        splitter = QSplitter(Qt.Orientation.Vertical)
        self.diff_table = QTableWidget(0, 5)
        self.diff_table.setHorizontalHeaderLabels(
            ["column", "residue A", "A", "B", "residue B"])
        self.diff_table.verticalHeader().setVisible(False)
        self.diff_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.diff_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self.diff_table.itemSelectionChanged.connect(self._on_difference)
        splitter.addWidget(self.diff_table)

        self.alignment_view = SequenceView()
        alignment_area = QScrollArea()
        alignment_area.setWidget(self.alignment_view)
        alignment_area.setWidgetResizable(True)
        splitter.addWidget(alignment_area)
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 1)
        box.addWidget(splitter, 1)
        return page

    # ------------------------------------------------------------------ data

    def set_structure(self, structure, species: str = "human") -> None:
        """Reload every available sequence for the current structure."""
        self.annotations = load_annotations(species)
        self.sequences = load_named_sequences(structure)
        self.by_key = {s.key: s for s in self.sequences}

        for combo in (self.sequence_combo, self.compare_a, self.compare_b):
            combo.blockSignals(True)
            combo.clear()
            for s in self.sequences:
                combo.addItem(s.label, s.key)
            combo.blockSignals(False)
        if len(self.sequences) > 1:
            self.compare_b.setCurrentIndex(1)
        self._refresh()
        missing = [] if any(s.dna for s in self.sequences) else ["coding sequence"]
        self.statusBar().showMessage(
            f"{len(self.sequences)} sequences"
            + (f" — {', '.join(missing)} not downloaded; "
               "run python -m piezo1.io.fetch" if missing else ""))

    def _current(self):
        return self.by_key.get(self.sequence_combo.currentData())

    def _styles(self, sequence) -> dict[int, ResidueStyle]:
        mode = self.colour_combo.currentData()
        styles: dict[int, ResidueStyle] = {}
        if sequence is None or mode in ("plain", "chemistry"):
            return {} if mode != "plain" else {
                p: ResidueStyle(foreground="#c8ccd4") for p in sequence.positions
            } if sequence else {}

        if mode == "domain":
            from ..render.colormaps import load_domain_palette
            palette = load_domain_palette(sequence.numbering)
            for entry in palette.domains:
                span = entry[palette.species]
                if span["start"] is None or span["end"] is None:
                    continue
                for residue in range(int(span["start"]), int(span["end"]) + 1):
                    styles[residue] = ResidueStyle(background=entry["color"],
                                                   tooltip=entry.get("name", ""))
        elif mode == "variants":
            for variant in self.annotations.variants:
                if variant.residue is None:
                    continue
                colour = {"GoF": "#f26d6d", "LoF": "#6fb1ff",
                          "engineered": "#8a919e"}.get(
                              variant.classification, "#f2a65a")
                styles[variant.residue] = ResidueStyle(
                    background=colour, tooltip=variant.label)
        elif mode == "sites":
            for group in self.annotations.residue_groups:
                for residue in group.residues:
                    styles[residue] = ResidueStyle(background="#7ed67e",
                                                  tooltip=group.label)
        return styles

    def _refresh(self) -> None:
        sequence = self._current()
        self.dna_check.setEnabled(bool(sequence and sequence.dna))
        self.view.set_sequence(sequence, self._styles(sequence),
                               show_dna=self.dna_check.isChecked())
        if sequence is not None:
            note = f" · {sequence.note}" if sequence.note else ""
            gaps = " · has gaps" if sequence.has_gaps else ""
            self.detail.setText(
                f"<b>{sequence.label}</b> — {len(sequence)} residues, "
                f"{sequence.numbering} numbering, source {sequence.source}"
                f"{note}{gaps}")

    # ------------------------------------------------------------ selection

    def _on_selection(self, first: int, last: int) -> None:
        sequence = self._current()
        if sequence is None:
            return
        count = len(self.view.selected_residues())
        text = (f"<b>{first}–{last}</b> ({count} residues): "
                f"{sequence.segment(first, last)[:60]}")
        if first == last:
            text = self._describe_residue(sequence, first)
        self.detail.setText(text)

    def _describe_residue(self, sequence, residue: int) -> str:
        letter = sequence.at(residue) or "?"
        parts = [f"<b>{letter}{residue}</b>"]
        codon = sequence.codon(residue)
        if codon:
            parts.append(f"codon <tt>{codon}</tt>")
        domain = self.annotations.domain_at(residue)
        if domain is not None:
            parts.append(f"domain <b>{domain.name}</b>")
        for variant in self.annotations.variants:
            if variant.residue == residue:
                parts.append(f"variant <b>{variant.label}</b> "
                             f"({variant.classification})")
        for group in self.annotations.residue_groups:
            if residue in group.residues:
                parts.append(f"site: {group.label}")
        return " · ".join(parts)

    def _on_hover(self, residue: int) -> None:
        self.statusBar().showMessage(f"residue {residue}")

    def _emit_selection(self) -> None:
        residues = self.view.selected_residues()
        if not residues:
            self.statusBar().showMessage("nothing selected")
            return
        sequence = self._current()
        label = f"sequence {residues[0]}–{residues[-1]}"
        if sequence is not None and sequence.numbering != "human":
            # Converting here rather than in the model keeps the one sanctioned
            # conversion path in core.sequence and out of the UI.
            from ..core.sequence import mouse_to_human
            converted = [mouse_to_human(r) for r in residues]
            residues = [r for r in converted if r is not None]
            label += " (converted from mouse numbering)"
        self.residues_selected.emit(residues, label)
        self.statusBar().showMessage(f"{len(residues)} residues sent to the model")

    def show_residues(self, residues) -> None:
        """Reflect a selection made elsewhere in the application."""
        if not residues:
            return
        self.view.select(min(residues), max(residues))

    def _find(self) -> None:
        sequence = self._current()
        text = self.find_edit.text().strip()
        if sequence is None or not text:
            return
        if text.isdigit():
            residue = int(text)
            if sequence.index_of(residue) is None:
                self.statusBar().showMessage(
                    f"residue {residue} is not in this sequence")
                return
            self.view.select(residue, residue)
            self.statusBar().showMessage(f"residue {residue}")
            return
        index = sequence.letters.upper().find(text.upper())
        if index < 0:
            self.statusBar().showMessage(f"motif {text!r} not found")
            return
        first = sequence.positions[index]
        last = sequence.positions[min(index + len(text) - 1,
                                      len(sequence.positions) - 1)]
        self.view.select(first, last)
        self.statusBar().showMessage(f"motif {text!r} at {first}–{last}")

    # ------------------------------------------------------------ comparison

    def _run_comparison(self) -> None:
        a = self.by_key.get(self.compare_a.currentData())
        b = self.by_key.get(self.compare_b.currentData())
        if a is None or b is None or a is b:
            self.compare_summary.setText("choose two different sequences")
            return
        method = self.method_combo.currentData()
        if method == "positional" and a.numbering != b.numbering:
            self.compare_summary.setText(
                "<b>Pairing by residue number needs a shared numbering.</b> "
                f"{a.label} is {a.numbering}, {b.label} is {b.numbering} — "
                "use the global alignment instead.")
            return

        self.statusBar().showMessage("aligning…")
        comparison = compare_sequences(a, b, method)
        self._comparison = comparison
        self.compare_summary.setText(
            f"<b>{comparison.identity * 100:.1f}% identity</b> over "
            f"{comparison.length} columns — {comparison.n_identical} identical, "
            f"{comparison.n_mismatch} substitutions, {comparison.n_gap} gaps")

        rows = comparison.differences[:2000]
        self.diff_table.setRowCount(len(rows))
        for i, (column, pa, pb, x, y) in enumerate(rows):
            for col, value in enumerate((str(column), str(pa or "-"), x, y,
                                         str(pb or "-"))):
                self.diff_table.setItem(i, col, QTableWidgetItem(value))
        self.diff_table.resizeColumnsToContents()
        self.alignment_view.set_sequence(
            a, {p: ResidueStyle(background="#f26d6d")
                for _c, p, _pb, _x, _y in comparison.differences if p})
        self.statusBar().showMessage(comparison.summary())

    def _on_difference(self) -> None:
        rows = self.diff_table.selectionModel().selectedRows()
        if not rows or self._comparison is None:
            return
        item = self.diff_table.item(rows[0].row(), 1)
        if item is None or item.text() == "-":
            return
        residue = int(item.text())
        self.alignment_view.select(residue, residue)
        self.residues_selected.emit([residue], f"difference at {residue}")
