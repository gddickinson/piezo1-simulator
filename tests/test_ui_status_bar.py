"""The status line must not be able to resize the window.

This project's status bar carries its caveats, and several are deliberately
long — a picture that reads as a measurement needs a whole sentence to say what
it is instead. In a plain ``QLabel`` that had a consequence nobody had looked
for: a non-wrapping label reports its full text width as a size hint, a
``QStatusBar`` passes it up as a minimum, and ``QMainWindow`` honours it, so
running an analysis *widened the window* — past the edge of the monitor in the
worst cases, where the title bar and the panel controls become unreachable.

The measured contrast, from :func:`test_a_plain_label_is_the_bug`: at 1500
characters a plain label demands a window over 12,000 pixels wide, and the
replacement demands the same width at 0 characters as at 4000.

The other half of the file is the part that keeps the fix honest. Eliding is a
*display* change; every guard in this suite that asserts a caveat cannot be
omitted reads ``status_label.text()``, and if that returned the shortened
string those guards would start passing vacuously. So ``text()`` returns the
whole message and ``displayed_text()`` returns what is painted.
"""

from __future__ import annotations

import pytest

from piezo1.ui.status_bar import HISTORY_LIMIT, StatusLog, StatusMessage

LONG = ("MODELLED MICELLE, NOT THE OBSERVED DENSITY. Figure 4b is the "
        "unsharpened cryo-EM map at 6 sigma; this is the surface 9 A outside "
        "the hydrophobic belt. ") * 4


@pytest.fixture
def message(qt_app):
    widget = StatusMessage()
    widget.resize(300, 20)
    return widget


# --------------------------------------------------------------------------
# The bug, and that it is fixed
# --------------------------------------------------------------------------

def test_a_plain_label_is_the_bug(qt_app):
    """The behaviour being replaced, measured rather than described.

    Without this, the fix below is a change with no demonstrated reason.
    """
    from PyQt6.QtWidgets import QLabel, QMainWindow, QStatusBar

    window = QMainWindow()
    bar = QStatusBar()
    window.setStatusBar(bar)
    label = QLabel()
    bar.addWidget(label, 1)
    window.resize(400, 300)
    window.show()
    qt_app.processEvents()

    label.setText("")
    qt_app.processEvents()
    narrow = window.minimumSizeHint().width()
    label.setText("X" * 1500)
    qt_app.processEvents()
    wide = window.minimumSizeHint().width()
    assert wide > narrow + 5000, (
        f"a plain QLabel should force the window to {wide} px; if it no "
        f"longer does, Qt has changed and this fix may be unnecessary")


def test_the_window_minimum_ignores_the_message_length(qt_app):
    """The fix: the same minimum at 0 characters as at 4000."""
    from PyQt6.QtWidgets import QMainWindow, QStatusBar

    from piezo1.ui import status_bar

    window = QMainWindow()
    bar = QStatusBar()
    window.setStatusBar(bar)
    label = status_bar.install(window, bar)
    window.resize(400, 300)
    window.show()
    qt_app.processEvents()

    widths = []
    for n in (0, 200, 1500, 4000):
        label.setText("X" * n)
        qt_app.processEvents()
        widths.append(window.minimumSizeHint().width())
    assert len(set(widths)) == 1, (
        f"the window minimum moved with the message: {widths}")


def test_the_widget_never_asks_for_more_than_its_floor(message):
    from piezo1.ui.status_bar import MINIMUM_WIDTH

    message.setText(LONG)
    assert message.minimumSizeHint().width() == MINIMUM_WIDTH
    assert message.sizeHint().width() == MINIMUM_WIDTH


# --------------------------------------------------------------------------
# Nothing is lost
# --------------------------------------------------------------------------

def test_the_full_text_survives_and_the_display_is_shortened(message):
    """The load-bearing distinction. Every caveat guard reads ``text()``."""
    message.setText(LONG)
    assert message.text() == LONG
    assert message.displayed_text() != LONG
    assert message.is_elided
    assert len(message.displayed_text()) < len(LONG)


def test_a_short_message_is_not_shortened_and_has_no_tooltip(message):
    message.setText("dome: 9.7 nm")
    assert message.displayed_text() == "dome: 9.7 nm"
    assert not message.is_elided
    assert message.toolTip() == ""


def test_a_long_message_puts_the_whole_thing_in_the_tooltip(message):
    message.setText(LONG)
    assert message.toolTip() == LONG


def test_widening_the_widget_shows_more_of_the_message(qt_app, message):
    message.setText(LONG)
    narrow = len(message.displayed_text())
    message.resize(1400, 20)
    qt_app.processEvents()
    assert len(message.displayed_text()) > narrow, (
        "the elision must follow the widget's width, or the message would "
        "stay truncated on a wide monitor for no reason")


def test_the_history_keeps_messages_that_flashed_past(message):
    """Controllers set one status while working and another when done, and
    before this the first was simply lost."""
    message.setText("building the modelled micelle envelope…")
    message.setText(LONG)
    assert len(message.history) == 2
    assert message.history[0].startswith("building")
    assert message.history[-1] == LONG


def test_a_repeated_message_is_not_recorded_twice(message):
    message.setText("same")
    message.setText("same")
    assert message.history == ["same"]


def test_the_history_is_bounded(message):
    for i in range(HISTORY_LIMIT + 50):
        message.setText(f"message {i}")
    assert len(message.history) == HISTORY_LIMIT
    assert message.history[-1] == f"message {HISTORY_LIMIT + 49}"


def test_the_log_shows_the_history_newest_first(qt_app, message):
    message.setText("first")
    message.setText("second")
    log = StatusLog(message)
    text = log.view.toPlainText()
    assert text.index("second") < text.index("first")
    assert "first" in text and "second" in text


def test_the_log_says_so_when_there_is_nothing(qt_app, message):
    log = StatusLog(message)
    assert "no messages yet" in log.view.toPlainText()


# --------------------------------------------------------------------------
# In the real window
# --------------------------------------------------------------------------

def test_the_window_keeps_the_whole_caveat(qt_app):
    """`status_label.text()` is what the caveat guards read."""
    from piezo1.ui.main_window import MainWindow

    window = MainWindow()
    window.resize(1000, 700)
    window._set_status(LONG)
    assert window.status_label.text() == LONG
    assert window.status_label.displayed_text() != LONG


def test_running_an_analysis_does_not_widen_the_window(qt_app):
    """The reported symptom, end to end."""
    from piezo1.ui.main_window import MainWindow

    window = MainWindow()
    window.resize(1100, 700)
    window.show()
    before = window.width()
    for _ in range(5):
        window._set_status(LONG * 3)
    assert window.width() == before
