"""Domains, functional-residue groups and disease variants.

The variant table is the scientific centre of the application, so it is
deliberately honest about coverage: a variant whose residue is not resolved in
the structure on screen is shown greyed out with the reason, rather than
selected to no visible effect.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtGui import QBrush, QColor
from PyQt6.QtWidgets import (QAbstractItemView, QComboBox, QHBoxLayout,
                             QHeaderView, QLabel, QLineEdit, QListWidget,
                             QListWidgetItem, QTableWidget, QTableWidgetItem,
                             QTabWidget, QVBoxLayout, QWidget)

from ...core.annotations import (Annotations, annotation_gap,
                                 load_annotations)

__all__ = ["AnnotationPanel"]

CLASS_COLORS = {
    "GoF": "#ff6b6b",
    "LoF": "#5b9dff",
    "VUS": "#b0b6c4",
    "blood-group": "#8ad35e",
    "engineered": "#ffd93d",
    "benign": "#6f7684",
}


class AnnotationPanel(QWidget):
    """Browse domains, functional sites and variants; select them on screen."""

    residues_selected = pyqtSignal(object, str)   # (list[int], label)
    selection_cleared = pyqtSignal()
    focus_requested = pyqtSignal(object)          # list[int]

    def __init__(self, species: str = "human", parent=None) -> None:
        super().__init__(parent)
        self._species = species
        self.annotations: Annotations = load_annotations(species)
        self._modelled: set[int] = set()
        self._current_pdb = ""

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)
        self.tabs = QTabWidget()
        layout.addWidget(self.tabs)

        self.tabs.addTab(self._build_domains(), "Domains")
        self.tabs.addTab(self._build_sites(), "Sites")
        self.tabs.addTab(self._build_variants(), "Variants")

        self.info = QLabel("Select an item to highlight it on the structure.")
        self.info.setWordWrap(True)
        self.info.setTextFormat(Qt.TextFormat.RichText)
        self.info.setOpenExternalLinks(True)
        self.info.setStyleSheet(
            "color:#c3cad8; background:#161a24; border:1px solid #262c3a;"
            "border-radius:5px; padding:7px; font-size:11px;")
        self.info.setMinimumHeight(64)
        layout.addWidget(self.info)

    # ------------------------------------------------------------- domains

    def _build_domains(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(4, 4, 4, 4)
        self.domain_list = QListWidget()
        self._refill_domains()
        self.domain_list.currentItemChanged.connect(self._on_domain)
        v.addWidget(self.domain_list)
        return w

    def _refill_domains(self) -> None:
        self.domain_list.clear()
        for d in self.annotations.domains:
            if d.start is None:
                continue
            item = QListWidgetItem(f"{d.name}   {d.start}–{d.end}")
            item.setData(Qt.ItemDataRole.UserRole, d.id)
            item.setForeground(QBrush(QColor(d.color)))
            if d.confidence != "high":
                item.setToolTip(f"confidence: {d.confidence}")
            self.domain_list.addItem(item)

    def _on_domain(self, item: QListWidgetItem | None) -> None:
        if item is None:
            return
        did = item.data(Qt.ItemDataRole.UserRole)
        d = next((x for x in self.annotations.domains if x.id == did), None)
        if d is None or d.start is None:
            return
        residues = list(range(d.start, d.end + 1))
        self.residues_selected.emit(residues, d.name)
        self.focus_requested.emit(residues)
        self.info.setText(
            f"<b>{d.name}</b> &nbsp; residues {d.start}–{d.end} "
            f"(mouse {d.mouse_start}–{d.mouse_end})<br>{d.description}<br>"
            f"<span style='color:#7f8798'>source: {d.source} · "
            f"confidence: {d.confidence}</span>")

    # --------------------------------------------------------------- sites

    def _build_sites(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(4, 4, 4, 4)
        self.site_list = QListWidget()
        self._refill_sites()
        self.site_list.currentItemChanged.connect(self._on_site)
        v.addWidget(self.site_list)
        return w

    def _refill_sites(self) -> None:
        self.site_list.clear()
        for g in self.annotations.residue_groups:
            res = ", ".join(str(r) for r in g.residues[:6])
            more = "…" if len(g.residues) > 6 else ""
            item = QListWidgetItem(f"{g.label}  ({res}{more})")
            item.setData(Qt.ItemDataRole.UserRole, g.id)
            item.setForeground(QBrush(QColor(g.color)))
            self.site_list.addItem(item)

    def _on_site(self, item: QListWidgetItem | None) -> None:
        if item is None:
            return
        g = self.annotations.group(item.data(Qt.ItemDataRole.UserRole))
        if g is None:
            return
        self.residues_selected.emit(list(g.residues), g.label)
        self.focus_requested.emit(list(g.residues))
        detail = " ".join(
            f"{d['human_aa']}{d['human']}"
            + ("" if d.get("conserved") else
               f"<span style='color:#ffb454'>(mouse {d['mouse_aa']}{d['mouse']})</span>")
            for d in g.detail)
        self.info.setText(
            f"<b>{g.label}</b><br>{detail}<br>{g.description}<br>"
            f"<span style='color:#7f8798'>evidence: {g.evidence} · "
            f"source: {g.source}</span>")

    # ------------------------------------------------------------ variants

    def _build_variants(self) -> QWidget:
        w = QWidget()
        v = QVBoxLayout(w)
        v.setContentsMargins(4, 4, 4, 4)

        row = QHBoxLayout()
        self.class_filter = QComboBox()
        self.class_filter.setToolTip('Filter variants by classification: gain-of-function,\nloss-of-function, engineered, uncertain, or blood-group.')
        self.class_filter.addItem("all classes")
        self.class_filter.addItems(self.annotations.variant_classes())
        self.class_filter.currentTextChanged.connect(self._refill_variants)
        row.addWidget(self.class_filter)
        self.search = QLineEdit()
        self.search.setPlaceholderText("filter by residue, phenotype…")
        self.search.textChanged.connect(self._refill_variants)
        row.addWidget(self.search, 1)
        v.addLayout(row)

        self.variant_table = QTableWidget(0, 4)
        self.variant_table.setHorizontalHeaderLabels(
            ["Variant", "Class", "Domain", "Phenotype"])
        self.variant_table.verticalHeader().setVisible(False)
        self.variant_table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self.variant_table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        hh = self.variant_table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(3, QHeaderView.ResizeMode.Stretch)
        self.variant_table.itemSelectionChanged.connect(self._on_variant)
        v.addWidget(self.variant_table)

        self.coverage_note = QLabel("")
        self.coverage_note.setStyleSheet("color:#7f8798; font-size:10px;")
        v.addWidget(self.coverage_note)

        self._refill_variants()
        return w

    def _refill_variants(self) -> None:
        want = self.class_filter.currentText()
        needle = self.search.text().strip().lower()
        rows = []
        for var in self.annotations.variants:
            if want != "all classes" and var.classification != want:
                continue
            if needle and needle not in " ".join(
                    str(x).lower() for x in
                    (var.label, var.phenotype, var.domain, var.functional_effect)):
                continue
            rows.append(var)
        self._variant_rows = rows

        self.variant_table.setRowCount(len(rows))
        for r, var in enumerate(rows):
            resolved = (not self._modelled) or (var.position in self._modelled)
            for c, text in enumerate((var.label, var.classification or "",
                                      var.domain or "", var.phenotype or "")):
                item = QTableWidgetItem(str(text))
                if c == 1 and var.classification in CLASS_COLORS:
                    item.setForeground(QBrush(QColor(CLASS_COLORS[var.classification])))
                if not resolved:
                    item.setForeground(QBrush(QColor("#565d6b")))
                    item.setToolTip(
                        f"residue {var.position} is not modelled in "
                        f"{self._current_pdb}")
                self.variant_table.setItem(r, c, item)

        if self._modelled:
            missing = sum(1 for v in rows
                          if v.position not in self._modelled)
            self.coverage_note.setText(
                f"{len(rows) - missing} of {len(rows)} shown variants are "
                f"resolved in {self._current_pdb}; {missing} greyed out are not "
                f"modelled there.")

    def _on_variant(self) -> None:
        rows = self.variant_table.selectionModel().selectedRows()
        if not rows:
            return
        var = self._variant_rows[rows[0].row()]
        if var.position is None:
            self.info.setText(f"<b>{var.label}</b> — no mapped residue number.")
            return
        # The variant keeps its published name and is marked at its position
        # in THIS entry's numbering — R2456H is R2456H everywhere, and on a
        # mouse entry it sits at 2482.
        self.residues_selected.emit([var.position], var.label)
        self.focus_requested.emit([var.position])

        modelled = ", ".join(var.modelled_in) if var.modelled_in else None
        warn = ("" if modelled else
                "<br><span style='color:#ff9f43'>⚠ not resolved in any human "
                "PIEZO1 structure — nothing will be highlighted</span>")
        cons = ("" if var.conserved is not False else
                f"<br><span style='color:#ffb454'>not conserved: mouse residue "
                f"{var.mouse_residue}</span>")
        pmid = (f" · <a style='color:#7aa7ff' "
                f"href='https://pubmed.ncbi.nlm.nih.gov/{var.pmid}/'>PMID "
                f"{var.pmid}</a>" if var.pmid else "")
        self.info.setText(
            f"<b>{var.label}</b> &nbsp; {var.classification} &nbsp; "
            f"{var.domain or ''}<br>"
            f"{var.phenotype or ''}<br>{var.functional_effect or ''}"
            f"{' — ' + var.effect_magnitude if var.effect_magnitude else ''}"
            f"{cons}{warn}{self._evidence_html(var.label)}<br>"
            f"<span style='color:#7f8798'>mouse equivalent: "
            f"{var.mouse_residue or 'n/a'} · resolved in: {modelled or 'none'}"
            f"{pmid}</span>")

    def _evidence_html(self, label: str) -> str:
        """How the direction was established, and where the two sources differ.

        The classification alone reads as a fact. It is not: for 20 of the 46
        directional variants the direction is *inferred* from which disease the
        variant causes rather than measured, and for one the two sources
        disagree outright. Both belong next to the label rather than in a report
        the user will not open.
        """
        from ...analysis.prediction_record import variant_evidence

        try:
            evidence = variant_evidence(label)
        except Exception:
            return ""
        if not evidence["in_analysis_set"]:
            return ""

        colour = {"measured": "#7fd18a",
                  "disease_mechanism": "#ffb454"}.get(evidence["evidence"], "#7f8798")
        text = (f"<br><span style='color:{colour}'>direction "
                f"<b>{evidence['evidence'].replace('_', ' ')}</b></span>"
                f"<span style='color:#7f8798'> — {evidence['evidence_note']}"
                "</span>")
        if evidence["conflict"]:
            conflict = evidence["conflict"]
            text += (f"<br><span style='color:#ff6b6b'>⚠ sources disagree: "
                     f"curated says {conflict['curated']}, the disease "
                     f"mechanism implies {conflict['inferred']} — this project "
                     f"reports the disagreement rather than resolving it"
                     f"</span>")
        return text

    # -------------------------------------------------------------- context

    def set_species(self, species: str) -> None:
        """Read the annotation in **this entry's** numbering.

        The panel used to load human annotation once, at construction, and
        never ask again — so on a mouse entry, which is most of the catalogue,
        every domain range, site and variant in it was a human residue number
        applied to mouse coordinates. The offset is not constant and reaches 26
        residues: picking "Transmembrane hydrophobic gate" highlighted human
        2447/2450/2454 where the gate of the structure on screen is mouse
        2473/2476/2480. It is the same defect Round 93 found in
        ``analysis/features.py`` and Round 89b found in the functional-residue
        loader, in the one panel a user reads first.

        A numbering with no curated annotation — PIEZO2, PEZO-1, dPIEZO — now
        yields **empty lists and the reason**, rather than PIEZO1's annotation
        drawn on another protein.
        """
        if species == self._species:
            return
        self._species = species
        self.annotations = load_annotations(species)
        self._refill_domains()
        self._refill_sites()
        self.class_filter.blockSignals(True)
        self.class_filter.clear()
        self.class_filter.addItem("all classes")
        self.class_filter.addItems(self.annotations.variant_classes())
        self.class_filter.blockSignals(False)
        self._refill_variants()
        gap = annotation_gap(species)
        self.info.setText(
            f"<span style='color:#ffb454'>No curated annotation in {species} "
            f"numbering.</span><br>{gap}" if gap else
            f"Annotation shown in <b>{species}</b> numbering — the numbering "
            f"of the structure on screen.")

    def species(self) -> str:
        return self._species

    def set_structure_context(self, pdb: str, modelled_residues: set[int]) -> None:
        """Tell the panel which residues the displayed structure resolves."""
        self._current_pdb = pdb
        self._modelled = set(modelled_residues)
        self._refill_variants()
