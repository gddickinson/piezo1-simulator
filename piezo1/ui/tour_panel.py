"""The guided tour, driving the real application rather than describing it.

Each step sets up the view and then triggers the *same* analysis the Analysis
and Physics panels use, so the number it reports is the number the application
computes. Nothing is narrated from a literal — see :mod:`piezo1.tour`.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt, pyqtSignal
from PyQt6.QtWidgets import (QHBoxLayout, QLabel, QListWidget, QProgressBar,
                             QPushButton, QSplitter, QTextBrowser, QVBoxLayout,
                             QWidget)

from ..parameters import PARAMETERS
from ..tour import TOUR

__all__ = ["TourPanel"]


class TourPanel(QWidget):
    """Step list, prose, and the live measurement for the current step."""

    step_requested = pyqtSignal(int)

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 6, 6, 6)

        header = QLabel(
            "<b>Guided tour of the mechanism.</b> Every number below is "
            "computed when the step runs — none is written into the text.")
        header.setWordWrap(True)
        layout.addWidget(header)

        splitter = QSplitter(Qt.Orientation.Vertical)

        self.steps = QListWidget()
        for step in TOUR:
            self.steps.addItem(step.title)
        self.steps.currentRowChanged.connect(self._on_row)
        self.steps.setMaximumHeight(150)
        splitter.addWidget(self.steps)

        self.body = QTextBrowser()
        self.body.setOpenExternalLinks(True)
        splitter.addWidget(self.body)
        splitter.setStretchFactor(1, 1)
        layout.addWidget(splitter, 1)

        self.measurement = QLabel("")
        self.measurement.setWordWrap(True)
        self.measurement.setTextFormat(Qt.TextFormat.RichText)
        self.measurement.setStyleSheet(
            "background:#1c2027; border:1px solid #2a2f38; padding:8px;")
        self.measurement.setMinimumHeight(72)
        layout.addWidget(self.measurement)

        self.sources = QLabel("")
        self.sources.setWordWrap(True)
        self.sources.setStyleSheet("color:#8a919e;")
        layout.addWidget(self.sources)

        self.progress = QProgressBar()
        self.progress.setRange(0, len(TOUR))
        self.progress.setTextVisible(True)
        layout.addWidget(self.progress)

        row = QHBoxLayout()
        self.back_button = QPushButton("← Back")
        self.back_button.clicked.connect(lambda: self._move(-1))
        row.addWidget(self.back_button)
        self.next_button = QPushButton("Next →")
        self.next_button.clicked.connect(lambda: self._move(1))
        row.addWidget(self.next_button)
        layout.addLayout(row)

        self.setMinimumWidth(300)

    # ------------------------------------------------------------- navigation

    def start(self) -> None:
        self.steps.setCurrentRow(0)

    def _move(self, delta: int) -> None:
        row = self.steps.currentRow() + delta
        if 0 <= row < len(TOUR):
            self.steps.setCurrentRow(row)

    def _on_row(self, row: int) -> None:
        if not (0 <= row < len(TOUR)):
            return
        step = TOUR[row]
        self.body.setHtml(_style(step.body_html()))
        self.measurement.setText(
            "<i>running…</i>" if step.run or step.measure else "")
        self.progress.setValue(row + 1)
        self.progress.setFormat(f"step {row + 1} of {len(TOUR)}")
        self.back_button.setEnabled(row > 0)
        self.next_button.setEnabled(row < len(TOUR) - 1)
        self._show_sources(step)
        self.step_requested.emit(row)

    def _show_sources(self, step) -> None:
        if not step.cites:
            self.sources.setText("")
            return
        parts = []
        for key in step.cites:
            try:
                parameter = PARAMETERS.get(key)
            except KeyError:
                continue
            value = f"{PARAMETERS.value(key):g}"
            unit = f" {parameter.unit}" if parameter.unit else ""
            parts.append(f"{parameter.name} {value}{unit} [{parameter.citation}]")
        self.sources.setText("Values used: " + "; ".join(parts) if parts else "")

    # --------------------------------------------------------------- results

    def set_measurement(self, text: str) -> None:
        self.measurement.setText(text or "")

    def current_step(self):
        row = self.steps.currentRow()
        return TOUR[row] if 0 <= row < len(TOUR) else None


def _style(body: str) -> str:
    return f"""<html><head><style>
      body {{ color:#c8ccd4; font-size:13px; line-height:1.55; }}
      b {{ color:#f0f3f8; }}
      i {{ color:#8a919e; }}
    </style></head><body>{body}</body></html>"""
