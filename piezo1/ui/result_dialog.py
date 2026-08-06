"""A readable window for an analysis that returns numbers rather than a picture.

Several analyses produce a table, not something to draw on the model — the ion
current, the interaction inventory, the variant-structure survey. They were
reachable only from the command line, which meant a GUI user could not get at
them at all. This renders one of the shared ``ANALYSES`` results as formatted
text, so nothing computable is command-line-only.

Deliberately generic: it takes whatever dict the analysis returns. A bespoke
panel per analysis would look better and would be one more thing to fall out of
step with the function it displays.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (QDialog, QDialogButtonBox, QLabel, QTextEdit,
                             QVBoxLayout)

__all__ = ["ResultDialog", "format_result"]


def format_result(data, indent: int = 0) -> list[str]:
    """Render a nested result as aligned key/value lines."""
    pad = "  " * indent
    lines: list[str] = []
    if isinstance(data, dict):
        width = max((len(str(k)) for k in data), default=0)
        for key, value in data.items():
            if isinstance(value, dict):
                lines.append(f"{pad}{key}:")
                lines += format_result(value, indent + 1)
            elif isinstance(value, (list, tuple)):
                if value and isinstance(value[0], (dict, list)):
                    lines.append(f"{pad}{key}:")
                    for item in value:
                        lines += format_result(item, indent + 1)
                else:
                    lines.append(f"{pad}{str(key):{width}} : "
                                 + ", ".join(_number(v) for v in value))
            else:
                lines.append(f"{pad}{str(key):{width}} : {_number(value)}")
    elif isinstance(data, (list, tuple)):
        for item in data:
            lines += format_result(item, indent)
    else:
        lines.append(f"{pad}{_number(data)}")
    return lines


def _number(value) -> str:
    if isinstance(value, bool) or value is None:
        return str(value)
    if isinstance(value, float):
        # Enough digits to be useful, not so many as to imply precision.
        return f"{value:.4g}"
    return str(value)


class ResultDialog(QDialog):
    """Non-modal window showing one analysis result."""

    def __init__(self, title: str, data, caveat: str = "", parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(False)
        self.resize(720, 560)

        layout = QVBoxLayout(self)
        if caveat:
            note = QLabel(caveat)
            note.setWordWrap(True)
            note.setStyleSheet("color:#d9a441;")
            layout.addWidget(note)

        view = QTextEdit()
        view.setReadOnly(True)
        font = QFont("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        view.setFont(font)
        view.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        view.setPlainText("\n".join(format_result(data)))
        layout.addWidget(view, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        buttons.accepted.connect(self.close)
        layout.addWidget(buttons)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
