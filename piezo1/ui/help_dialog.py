"""The Help window: a topic list, the text, and links to the shipped documents."""

from __future__ import annotations

import subprocess
import sys

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QDesktopServices
from PyQt6.QtCore import QUrl
from PyQt6.QtWidgets import (QDialog, QHBoxLayout, QLabel, QListWidget,
                             QPushButton, QSplitter, QTableWidget,
                             QTableWidgetItem, QTabWidget, QTextBrowser,
                             QVBoxLayout, QWidget)

from ..config import PROJECT_ROOT
from .help_content import DOC_LINKS, SHORTCUTS, TOPICS, topic_html

__all__ = ["HelpDialog", "open_document"]


def open_document(relative: str) -> bool:
    """Open a shipped document in whatever the system uses for it.

    Returns False when the file is not there — several documents are generated
    rather than committed, and a dialog that silently does nothing is worse
    than one that says the file has not been built yet.
    """
    path = PROJECT_ROOT / relative
    if not path.exists():
        return False
    if sys.platform == "darwin":
        subprocess.Popen(["open", str(path)])
    elif sys.platform.startswith("win"):
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))
    else:
        subprocess.Popen(["xdg-open", str(path)])
    return True


class HelpDialog(QDialog):
    """Feature guide, keyboard shortcuts and the document index."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("PIEZO1 Simulator — Help")
        self.resize(920, 620)
        # Not modal: the point of the guide is to be read while driving the
        # application, and a modal dialog would block exactly that.
        self.setModal(False)

        tabs = QTabWidget()
        tabs.addTab(self._build_guide(), "Feature guide")
        tabs.addTab(self._build_shortcuts(), "Keyboard and mouse")
        tabs.addTab(self._build_documents(), "Documents")

        layout = QVBoxLayout(self)
        layout.addWidget(tabs)
        close = QPushButton("Close")
        close.clicked.connect(self.accept)
        row = QHBoxLayout()
        row.addStretch(1)
        row.addWidget(close)
        layout.addLayout(row)

    # ----------------------------------------------------------------- guide

    def _build_guide(self) -> QWidget:
        self.topics = QListWidget()
        self.topics.setMaximumWidth(210)
        for title, _ in TOPICS:
            self.topics.addItem(title)

        self.body = QTextBrowser()
        self.body.setOpenExternalLinks(True)
        self.topics.currentRowChanged.connect(self._show_topic)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self.topics)
        splitter.addWidget(self.body)
        splitter.setStretchFactor(1, 1)

        page = QWidget()
        box = QVBoxLayout(page)
        box.addWidget(splitter)
        self.topics.setCurrentRow(0)
        return page

    def show_topic_named(self, name: str) -> None:
        """Open the guide at a named topic, for context-sensitive help."""
        for i, (title, _) in enumerate(TOPICS):
            if title.lower().startswith(name.lower()):
                self.topics.setCurrentRow(i)
                return

    def _show_topic(self, row: int) -> None:
        if 0 <= row < len(TOPICS):
            title, body = TOPICS[row]
            self.body.setHtml(topic_html(title, body))

    # ------------------------------------------------------------- shortcuts

    def _build_shortcuts(self) -> QWidget:
        table = QTableWidget(len(SHORTCUTS), 2)
        table.setHorizontalHeaderLabels(["Input", "Action"])
        table.verticalHeader().setVisible(False)
        table.horizontalHeader().setStretchLastSection(True)
        for i, (key, action) in enumerate(SHORTCUTS):
            table.setItem(i, 0, QTableWidgetItem(key))
            table.setItem(i, 1, QTableWidgetItem(action))
        table.resizeColumnsToContents()
        table.setEditTriggers(table.EditTrigger.NoEditTriggers)

        page = QWidget()
        box = QVBoxLayout(page)
        box.addWidget(table)
        return page

    # ------------------------------------------------------------- documents

    def _build_documents(self) -> QWidget:
        page = QWidget()
        box = QVBoxLayout(page)
        box.addWidget(QLabel(
            "The project's written record. These open in your system viewer."))

        self.doc_table = QTableWidget(len(DOC_LINKS), 2)
        self.doc_table.setHorizontalHeaderLabels(["Document", "Contains"])
        self.doc_table.verticalHeader().setVisible(False)
        self.doc_table.horizontalHeader().setStretchLastSection(True)
        for i, (title, path, description) in enumerate(DOC_LINKS):
            item = QTableWidgetItem(title)
            item.setData(Qt.ItemDataRole.UserRole, path)
            if not (PROJECT_ROOT / path).exists():
                item.setText(f"{title}  (not generated)")
                item.setForeground(Qt.GlobalColor.gray)
            self.doc_table.setItem(i, 0, item)
            self.doc_table.setItem(i, 1, QTableWidgetItem(description))
        self.doc_table.resizeColumnsToContents()
        self.doc_table.setEditTriggers(self.doc_table.EditTrigger.NoEditTriggers)
        self.doc_table.cellDoubleClicked.connect(self._open_row)
        box.addWidget(self.doc_table)

        self.doc_status = QLabel("")
        self.doc_status.setStyleSheet("color:#8a919e;")
        box.addWidget(self.doc_status)

        button = QPushButton("Open selected document")
        button.clicked.connect(
            lambda: self._open_row(self.doc_table.currentRow(), 0))
        box.addWidget(button)
        return page

    def _open_row(self, row: int, _column: int = 0) -> None:
        item = self.doc_table.item(row, 0) if row >= 0 else None
        if item is None:
            return
        path = item.data(Qt.ItemDataRole.UserRole)
        if open_document(path):
            self.doc_status.setText(f"opened {path}")
        else:
            self.doc_status.setText(
                f"{path} has not been generated — see the roadmap for which "
                f"script builds it")
