"""A status line that cannot resize the window.

The status bar carries this application's caveats, and several of them are
deliberately long — the HaloTag fold's "orientation is undetermined", the
micelle's "not the observed density", the pore surface's two sentences about
what a probe sphere is not. They are long because they have to be: a picture
that reads as a measurement needs a whole sentence to say what it is instead,
and a test guards each one against being shortened away.

Put in a plain ``QLabel``, that had a consequence nobody had looked for. A
non-wrapping ``QLabel`` reports its **full text width** as its size hint, a
``QStatusBar`` passes that up as a minimum, and ``QMainWindow`` honours it — so
running an analysis with a 600-character caveat *widened the window*, in the
worst cases past the edge of the monitor, where the title bar and the panel
controls become unreachable. The caveats were making the application unusable
in proportion to how careful they were.

So the display and the text are separated:

* :meth:`StatusMessage.text` returns the **whole** message. Every guard that
  checks a caveat is present reads this, so eliding cannot weaken one.
* the widget *draws* an elided version, and reports a size hint that ignores
  the message entirely, so nothing it is given can move the window;
* the full text is reachable three ways — the tooltip, a click, and the
  **Messages** window, which keeps a scrollable history so a caveat that
  flashed past during a long analysis can still be read.

The history matters more than it looks. Several controllers set a status while
they work and another when they finish, and before this the first was simply
lost.
"""

from __future__ import annotations

from PyQt6.QtCore import QSize, Qt, pyqtSignal
from PyQt6.QtGui import QFontMetrics
from PyQt6.QtWidgets import (QDialog, QHBoxLayout, QLabel, QPlainTextEdit,
                             QPushButton, QSizePolicy, QVBoxLayout)

__all__ = ["StatusMessage", "StatusLog", "HISTORY_LIMIT", "MINIMUM_WIDTH"]

#: How many past messages the log keeps. Enough to cover a session's worth of
#: analyses without the window becoming its own scrolling problem.
HISTORY_LIMIT = 200

#: The width the widget claims it needs, in pixels — deliberately small and
#: independent of the message. This is the number that stops a long caveat
#: from widening the window.
MINIMUM_WIDTH = 60


class StatusMessage(QLabel):
    """Status text: full in :meth:`text`, elided on screen, never resizing."""

    clicked = pyqtSignal()

    def __init__(self, text: str = "", parent=None) -> None:
        super().__init__(parent)
        self._full = ""
        self._rendered_at = -1
        self.history: list[str] = []
        # Ignored horizontally: the layout may give it whatever is left over,
        # and it will never ask for more.
        self.setSizePolicy(QSizePolicy.Policy.Ignored,
                           QSizePolicy.Policy.Preferred)
        self.setTextInteractionFlags(Qt.TextInteractionFlag.NoTextInteraction)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setText(text)

    # --------------------------------------------------------------- content

    def setText(self, text: str) -> None:          # noqa: N802 (Qt naming)
        """Record the whole message and display as much of it as fits."""
        text = "" if text is None else str(text)
        if text and (not self.history or self.history[-1] != text):
            self.history.append(text)
            del self.history[:-HISTORY_LIMIT]
        self._full = text
        # The tooltip is the cheapest way to the rest of a long message, and
        # the only one that needs no click.
        self.setToolTip(text if self._is_elided(text) else "")
        self._render()

    def text(self) -> str:                          # noqa: N802 (Qt naming)
        """The **whole** message, not what is on screen.

        Every guard asserting that a caveat cannot be omitted reads this, so
        the elision is display-only and cannot make one of them pass vacuously.
        """
        return self._full

    @property
    def is_elided(self) -> bool:
        return self._is_elided(self._full)

    # --------------------------------------------------------------- display

    def _is_elided(self, text: str) -> bool:
        if not text:
            return False
        return QFontMetrics(self.font()).horizontalAdvance(text) > self._room()

    def _room(self) -> int:
        return max(self.width() - 8, 1)

    def _render(self) -> None:
        self._rendered_at = self.width()
        elided = QFontMetrics(self.font()).elidedText(
            self._full, Qt.TextElideMode.ElideRight, self._room())
        super().setText(elided)

    def _render_if_resized(self) -> None:
        """Re-elide when the widget has changed width since the last render.

        Belt and braces beside ``resizeEvent``, which is not delivered to a
        widget that has never been shown — so a resize before the first show
        would otherwise leave a message elided to the old width, and stay that
        way until something else happened to repaint it.
        """
        if self.width() != self._rendered_at:
            self._render()

    def displayed_text(self) -> str:
        """What is actually painted — elided when the message does not fit."""
        self._render_if_resized()
        return QLabel.text(self)

    def resizeEvent(self, event) -> None:           # noqa: N802 (Qt naming)
        super().resizeEvent(event)
        self._render()

    def paintEvent(self, event) -> None:            # noqa: N802 (Qt naming)
        self._render_if_resized()
        super().paintEvent(event)

    # ------------------------------------------------------------- geometry

    def minimumSizeHint(self) -> QSize:             # noqa: N802 (Qt naming)
        """Independent of the message — the whole point of this class."""
        return QSize(MINIMUM_WIDTH, super().minimumSizeHint().height())

    def sizeHint(self) -> QSize:                    # noqa: N802 (Qt naming)
        return QSize(MINIMUM_WIDTH, super().sizeHint().height())

    # ---------------------------------------------------------------- events

    def mousePressEvent(self, event) -> None:       # noqa: N802 (Qt naming)
        self.clicked.emit()
        super().mousePressEvent(event)


class StatusLog(QDialog):
    """The full current message, and the scrollable history behind it."""

    def __init__(self, source: StatusMessage, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("PIEZO1 — messages")
        self.resize(760, 420)
        self.source = source

        layout = QVBoxLayout(self)
        self.view = QPlainTextEdit()
        self.view.setReadOnly(True)
        self.view.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        layout.addWidget(self.view, 1)

        row = QHBoxLayout()
        note = QLabel("Newest first. The status bar shortens long messages to "
                      "fit; nothing is lost here.")
        note.setStyleSheet("color:#8a919e;")
        note.setWordWrap(True)
        row.addWidget(note, 1)
        for label, slot in (("Copy", self._copy), ("Refresh", self.refresh),
                            ("Close", self.close)):
            button = QPushButton(label)
            button.clicked.connect(slot)
            row.addWidget(button)
        layout.addLayout(row)
        self.refresh()

    def refresh(self) -> None:
        entries = list(reversed(self.source.history)) or ["(no messages yet)"]
        self.view.setPlainText("\n\n".join(entries))
        self.view.moveCursor(self.view.textCursor().MoveOperation.Start)

    def _copy(self) -> None:
        from PyQt6.QtWidgets import QApplication

        clipboard = QApplication.clipboard()
        if clipboard is not None:
            clipboard.setText(self.view.toPlainText())


def install(window, bar) -> StatusMessage:
    """Put a :class:`StatusMessage` and its log button on a status bar.

    Returns the message widget, which the window keeps as ``status_label`` —
    the name every caller already uses, and it still answers ``setText`` and
    ``text``.
    """
    message = StatusMessage("Starting…")
    bar.addWidget(message, 1)

    button = QPushButton("⋯")
    button.setFixedWidth(26)
    button.setToolTip("Show the full message and the history (long messages "
                      "are shortened to fit the status bar)")
    button.setFlat(True)
    bar.addPermanentWidget(button)

    def open_log() -> None:
        log = getattr(window, "_status_log", None)
        if log is None:
            log = StatusLog(message, window)
            window._status_log = log
        log.refresh()
        log.show()
        log.raise_()
        log.activateWindow()

    button.clicked.connect(open_log)
    message.clicked.connect(open_log)
    return message
