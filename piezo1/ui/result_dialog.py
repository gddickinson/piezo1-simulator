"""A readable window for an analysis that returns numbers rather than a picture.

Several analyses produce a table, not something to draw on the model — the ion
current, the interaction inventory, the variant-structure survey. They were
reachable only from the command line, which meant a GUI user could not get at
them at all. This renders one of the shared ``ANALYSES`` results as formatted
text, so nothing computable is command-line-only.

Deliberately generic: it takes whatever dict the analysis returns. A bespoke
panel per analysis would look better and would be one more thing to fall out of
step with the function it displays.

Every one of these windows carries an **Explore** button. The table is where a
piece of reasoning ends, and the figure it came from, the model behind the
number and the same result drawn on the structure were all reachable only from
``docs/img`` and the scripts — which is the gap Round 34 closed for the
analyses themselves. What each analysis has to show is declared in
:mod:`piezo1.ui.exhibits`; the button opens
:class:`piezo1.ui.explore_window.ExploreWindow` and hands it the result already
displayed here, so nothing is recomputed and the two windows cannot be of
different runs.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont
from PyQt6.QtWidgets import (QDialog, QDialogButtonBox, QLabel, QPushButton,
                             QTextEdit, QVBoxLayout)

__all__ = ["ResultDialog", "format_result", "provenance_line"]


def provenance_line(structure_name: str = "", species: str = "") -> str:
    """Which structure, and which parameter set, a number was computed under.

    Stamped onto every result window. Two of the three hazards Round 50 set out
    to audit are exactly this: a number read against the wrong structure when
    several are displayed, and a number read as the documented one when the
    registry has been overridden. Neither is visible in a table of values, and
    this window is non-modal, so it can outlive the state that produced it.
    """
    from ..parameters import PARAMETERS

    where = structure_name or "unknown structure"
    if species:
        where += f" ({species} numbering)"
    if PARAMETERS.modified:
        return (f"{where} · ⚠ NON-DEFAULT PARAMETERS: "
                f"{PARAMETERS.override_summary()} — not comparable with "
                f"docs/SCIENCE.md")
    return f"{where} · parameters at documented defaults"


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

    def __init__(self, title: str, data, caveat: str = "", parent=None,
                 structure_name: str = "", species: str = "",
                 explore=None, n_exhibits: int = 0) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(False)
        self.resize(720, 560)

        layout = QVBoxLayout(self)

        # Recorded at construction, which is when the numbers were computed —
        # the window outlives the state, so a later override must not rewrite
        # what this says.
        self.provenance = provenance_line(structure_name, species)
        stamp = QLabel(self.provenance)
        stamp.setWordWrap(True)
        stamp.setStyleSheet(
            "color:#f2a65a;font-weight:bold;" if "NON-DEFAULT" in self.provenance
            else "color:#6f7684;")
        layout.addWidget(stamp)

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

        # A table is the end of a piece of reasoning; the button is the rest of
        # it. Added on the left of the box rather than as a link in the text,
        # because a result that has figures and a live model behind it should
        # not depend on the reader knowing they exist.
        self.explore_button = QPushButton("Explore these findings…")
        self.explore_button.setToolTip(
            "Figures, charts built from these very numbers, models you can "
            "drive, and the same result drawn on the structure.")
        if explore is None or n_exhibits == 0:
            self.explore_button.setEnabled(False)
            self.explore_button.setToolTip(
                "Nothing is registered to illustrate this result yet.")
        else:
            self.explore_button.clicked.connect(lambda: explore())
        buttons.addButton(self.explore_button,
                          QDialogButtonBox.ButtonRole.ActionRole)
        layout.addWidget(buttons)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, False)
