"""The four content panes an exhibit can be shown in.

Split from :mod:`piezo1.ui.explore_window` at the seam ``fold_view`` uses: the
window owns the lifecycle — which exhibit is selected, what the header says,
which structure and result it was opened for — and this owns the geometry of
showing one. Nothing here reaches back into the application.

Each pane is written around the way its kind of content misleads:

* a **figure** may not be on disk, and a broken image would look like a broken
  application rather than an ungenerated file, so it degrades to the command
  that builds it;
* a **chart** can be empty for a good reason, so the note is part of the pane
  and not decoration;
* a **simulation** invites the reading that the curve was measured, so the
  registry default is marked on every control and the caveat sits under the
  plot where it cannot be scrolled away;
* a **model** button changes the 3-D view, which is the most persuasive
  surface in the application, so it reports what it did.
"""

from __future__ import annotations

import math

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont, QMovie, QPixmap
from PyQt6.QtWidgets import (QGridLayout, QHBoxLayout, QLabel, QPushButton,
                             QScrollArea, QSizePolicy, QSlider, QVBoxLayout,
                             QWidget)

from .exhibit_chart import ChartView

__all__ = ["FigurePane", "ChartPane", "SimulationPane", "ModelPane",
           "note_label"]

NOTE_STYLE = "color:#8a919e;"
CAVEAT_STYLE = "color:#d9a441;"


def note_label(text: str, style: str = NOTE_STYLE) -> QLabel:
    label = QLabel(text)
    label.setWordWrap(True)
    label.setStyleSheet(style)
    label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return label


def _mono(text: str) -> QLabel:
    label = QLabel(text)
    font = QFont("Menlo")
    font.setStyleHint(QFont.StyleHint.Monospace)
    label.setFont(font)
    label.setWordWrap(True)
    label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
    return label


class FigurePane(QWidget):
    """A generated figure, scaled to the pane, or the command that makes it."""

    def __init__(self, exhibit, parent=None) -> None:
        super().__init__(parent)
        self.exhibit = exhibit
        self._pixmap: QPixmap | None = None
        self._movie: QMovie | None = None

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.image = QLabel()
        self.image.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self.image.setSizePolicy(QSizePolicy.Policy.Ignored,
                                 QSizePolicy.Policy.Ignored)
        area = QScrollArea()
        area.setWidgetResizable(True)
        area.setWidget(self.image)
        layout.addWidget(area, 1)
        layout.addWidget(_mono(f"docs/img/{exhibit.figure}  ·  rebuild with:  "
                               f"{exhibit.rebuild}"))
        self._load()

    def _load(self) -> None:
        path = self.exhibit.figure_file()
        if path is None:
            self.image.setText(
                f"This figure has not been generated in this clone.\n\n"
                f"{self.exhibit.rebuild}\n\n"
                f"Figures under docs/img are regenerable outputs; a missing "
                f"one is not an error.")
            self.image.setStyleSheet(NOTE_STYLE)
            self.image.setWordWrap(True)
            return
        if path.suffix.lower() == ".gif":
            # An animation is the one figure kind that shows a *process*, so it
            # is played rather than shown as its first frame.
            self._movie = QMovie(str(path))
            self.image.setMovie(self._movie)
            self._movie.start()
            return
        self._pixmap = QPixmap(str(path))
        self._rescale()

    def _rescale(self) -> None:
        if self._pixmap is None or self._pixmap.isNull():
            return
        width = max(self.image.width() - 8, 120)
        scaled = self._pixmap.scaledToWidth(
            min(width, self._pixmap.width()),
            Qt.TransformationMode.SmoothTransformation)
        self.image.setPixmap(scaled)

    def resizeEvent(self, event) -> None:      # noqa: N802 (Qt naming)
        super().resizeEvent(event)
        self._rescale()

    def close_pane(self) -> None:
        if self._movie is not None:
            self._movie.stop()


class ChartPane(QWidget):
    """A chart built from the result already in the window."""

    def __init__(self, chart, parent=None) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        self.view = ChartView()
        self.view.set_chart(chart)
        layout.addWidget(self.view, 1)
        if chart.note and not chart.empty:
            layout.addWidget(note_label(chart.note))

    def close_pane(self) -> None:
        pass


class SimulationPane(QWidget):
    """Sliders driving a model, and the curve they produce.

    The controls are built from :class:`~piezo1.ui.exhibit_models.Control`, so
    a control that *is* a registered parameter starts at the registry's own
    value and says so. Recomputation is debounced rather than run on every
    pixel of slider travel: one of these reads coordinates, and a redraw per
    mouse-move event would make the slider feel broken.
    """

    STEPS = 400
    DEBOUNCE_MS = 40

    def __init__(self, simulation, context, parent=None) -> None:
        super().__init__(parent)
        self.simulation = simulation
        self.context = context
        self._sliders: dict[str, QSlider] = {}
        self._readouts: dict[str, QLabel] = {}
        self._pristine: set[str] = {c.key for c in simulation.controls}

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._controls())
        self.view = ChartView()
        layout.addWidget(self.view, 1)
        self.note = note_label("")
        layout.addWidget(self.note)
        layout.addWidget(note_label(
            f"SENSITIVITY, NOT A MEASUREMENT. {simulation.caveat} Moving a "
            f"control changes nothing outside this window — the parameter "
            f"registry is not written to, so reports and the claims verifier "
            f"are unaffected.", CAVEAT_STYLE))

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.timeout.connect(self._recompute)
        self._recompute()

    def _controls(self) -> QWidget:
        box = QWidget()
        grid = QGridLayout(box)
        grid.setContentsMargins(0, 0, 0, 6)
        for row, control in enumerate(self.simulation.controls):
            label = QLabel(control.label + (f" ({control.unit})"
                                            if control.unit else ""))
            slider = QSlider(Qt.Orientation.Horizontal)
            slider.setRange(0, self.STEPS)
            slider.setValue(self._to_slider(control, control.start()))
            readout = QLabel()
            readout.setMinimumWidth(120)
            readout.setStyleSheet(NOTE_STYLE)
            self._sliders[control.key] = slider
            self._readouts[control.key] = readout
            grid.addWidget(label, row, 0)
            grid.addWidget(slider, row, 1)
            grid.addWidget(readout, row, 2)
            # Connected *after* the initial position is set, so an untouched
            # control stays "pristine" and reports the registry's own value
            # rather than the nearest of 400 slider steps to it. The opening
            # curve is then the model at its documented defaults, exactly.
            slider.valueChanged.connect(self._changed)
            if control.parameter:
                grid.addWidget(note_label(f"registered as {control.parameter}"),
                               row, 3)
        return box

    # ---------------------------------------------------- slider arithmetic

    def _to_slider(self, control, value: float) -> int:
        low, high = control.low, control.high
        value = min(max(value, low), high)
        if control.log and low > 0:
            fraction = ((math.log10(value) - math.log10(low))
                        / (math.log10(high) - math.log10(low)))
        else:
            fraction = (value - low) / max(high - low, 1e-12)
        return int(round(fraction * self.STEPS))

    def _from_slider(self, control, position: int) -> float:
        fraction = position / self.STEPS
        if control.log and control.low > 0:
            return float(10.0 ** (math.log10(control.low) + fraction *
                                  (math.log10(control.high)
                                   - math.log10(control.low))))
        return float(control.low + fraction * (control.high - control.low))

    def _half_step(self, control) -> float:
        """Half the value change one slider step makes, at the current
        position — which on a log control depends on where it is."""
        position = self._sliders[control.key].value()
        low = self._from_slider(control, max(position - 1, 0))
        high = self._from_slider(control, min(position + 1, self.STEPS))
        return abs(high - low) / 2.0

    def values(self) -> dict:
        return {c.key: (c.start() if c.key in self._pristine
                        else self._from_slider(c, self._sliders[c.key].value()))
                for c in self.simulation.controls}

    def _changed(self) -> None:
        sender = self.sender()
        for key, slider in self._sliders.items():
            if slider is sender:
                self._pristine.discard(key)
        self._timer.start(self.DEBOUNCE_MS)

    def _recompute(self) -> None:
        from .exhibit_models import run_simulation

        values = self.values()
        for control in self.simulation.controls:
            value = values[control.key]
            default = control.start()
            # Within half a slider step *is* the default: the slider quantises,
            # so an exact comparison would say "moved" about a control nobody
            # had touched, and the registry's own value would never be marked.
            at_default = (control.key in self._pristine
                          or abs(value - default) <= self._half_step(control))
            self._readouts[control.key].setText(
                f"{value:.4g}" + ("  (default)" if at_default
                                  else f"  (default {default:.4g})"))
        chart = run_simulation(self.simulation.key, values, self.context)
        self.view.set_chart(chart)
        self.note.setText(chart.note)
        self.note.setVisible(bool(chart.note) and not chart.empty)

    def close_pane(self) -> None:
        self._timer.stop()


class ModelPane(QWidget):
    """A button that turns on the overlay, and what it will draw."""

    def __init__(self, exhibit, action, parent=None, reason: str = "") -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 8, 0, 0)
        self.status = note_label("")
        row = QHBoxLayout()
        self.button = QPushButton(action.label if action else "unavailable")
        self.button.setEnabled(action is not None)
        self.button.clicked.connect(self._run)
        row.addWidget(self.button)
        row.addStretch(1)
        layout.addLayout(row)
        layout.addWidget(self.status)
        layout.addStretch(1)
        self.exhibit, self.action = exhibit, action
        if action is None:
            self.status.setText(reason or "no control in this window draws that")

    def _run(self) -> None:
        try:
            self.status.setText(self.action.run())
        except Exception as exc:               # noqa: BLE001 — never crash
            self.status.setText(f"could not draw it: "
                                f"{type(exc).__name__}: {exc}")

    def close_pane(self) -> None:
        pass
