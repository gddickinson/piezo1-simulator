"""The window behind every result window's **Explore** button.

An analysis window is a table. This is what the table came out of: the figure,
the curve the model traces when its inputs move, and the same thing drawn on
the structure. One window per analysis, listing that analysis's exhibits from
:mod:`piezo1.ui.exhibits`, with the content panes in
:mod:`piezo1.ui.exhibit_panes`.

Two things it must carry and cannot be allowed to lose:

**The provenance stamp.** It is passed in from the result window rather than
recomputed, for the reason ``ResultDialog`` records: both windows are
non-modal, and a stamp read later could describe a registry that has since
moved. A chart with no stamp is worse than a table with none — nobody reads a
picture sceptically.

**The one-control rule for anything drawn on the model.** An exhibit that puts
something in the 3-D view drives the *same* menu entry or panel button the user
would click, never the controller behind it. Otherwise the overlay would be on
while its menu entry said it was off, and there would be two answers to what is
being drawn. The displays themselves — load, superpose, show one component,
mark a residue set, recolour, morph — and the control each presses are in
:mod:`piezo1.ui.model_actions`, re-exported here because that is where the
window's callers look for them.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QDialogButtonBox, QLabel, QListWidget,
                             QListWidgetItem, QMainWindow, QSplitter,
                             QVBoxLayout, QWidget)

from .exhibit_panes import (ChartPane, FigurePane, ModelPane, SimulationPane,
                            note_label)
from .exhibits import BASES, exhibits_for
from .model_actions import (MODEL_ACTIONS, BoundAction, ModelActionSpec,
                            find_button, find_menu_action)

__all__ = ["ExploreWindow", "MODEL_ACTIONS", "ModelActionSpec",
           "BoundAction", "find_menu_action", "find_button"]


def _empty(layout) -> None:
    """Remove and destroy everything in a layout, nested layouts included."""
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget is not None:
            widget.setParent(None)
            widget.deleteLater()
        child = item.layout()
        if child is not None:
            _empty(child)
            child.deleteLater()


class ExploreWindow(QMainWindow):
    """Figures, charts, simulations and overlays for one analysis."""

    def __init__(self, analysis: str, title: str, data: dict,
                 window=None, provenance: str = "", structure=None,
                 species: str = "human") -> None:
        super().__init__(window)
        self.setWindowTitle(f"Explore — {title}")
        self.resize(1080, 720)
        self.analysis = analysis
        self.data = data or {}
        self.win = window
        self.structure = structure
        self.species = species
        self.exhibits = exhibits_for(analysis)
        self._pane = None

        central = QWidget()
        outer = QVBoxLayout(central)
        outer.setContentsMargins(10, 10, 10, 10)
        if provenance:
            stamp = QLabel(provenance)
            stamp.setWordWrap(True)
            stamp.setStyleSheet(
                "color:#f2a65a;font-weight:bold;" if "NON-DEFAULT" in provenance
                else "color:#6f7684;")
            outer.addWidget(stamp)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        self.list = QListWidget()
        self.list.setMaximumWidth(360)
        # Elision defeats the wrapping: a title long enough to need two lines
        # is exactly the one that gets a "..." instead of them.
        self.list.setTextElideMode(Qt.TextElideMode.ElideNone)
        # Wrapped, with the kind on its own line: padded to a column the
        # titles were clipped mid-word, and "How far would the constricti"
        # tells a reader nothing about what the exhibit is.
        self.list.setWordWrap(True)
        for item in self.exhibits:
            entry = QListWidgetItem(f"{item.kind.upper()}\n{item.title}")
            entry.setToolTip(item.what)
            self.list.addItem(entry)
        self.list.currentRowChanged.connect(self._select)
        splitter.addWidget(self.list)

        self.content = QWidget()
        self.content_layout = QVBoxLayout(self.content)
        self.content_layout.setContentsMargins(12, 0, 0, 0)
        splitter.addWidget(self.content)
        splitter.setStretchFactor(1, 1)
        # Given only a stretch factor the splitter sizes the list by its hint
        # and the two-line titles come back elided; this gives it the width
        # the wrapping needs.
        splitter.setSizes([320, 760])
        outer.addWidget(splitter, 1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.close)
        outer.addWidget(buttons)
        self.setCentralWidget(central)

        if self.exhibits:
            self.list.setCurrentRow(0)
        else:
            self.content_layout.addWidget(note_label(
                "Nothing is registered for this analysis yet."))

    # ------------------------------------------------------------- content

    def _clear(self) -> None:
        """Empty the content side completely before rebuilding it.

        Recursive, and it has to be: taking a *nested layout* out of a layout
        leaves its widgets parented to the panel, so they keep their old
        geometry and go on painting. The symptom was the previous exhibit's
        one-line chip left hanging over an otherwise blank panel — which no
        test caught, because ``findChildren`` finds a widget whether or not it
        is where it belongs.
        """
        if self._pane is not None and hasattr(self._pane, "close_pane"):
            self._pane.close_pane()
        self._pane = None
        _empty(self.content_layout)

    def _select(self, row: int) -> None:
        if not 0 <= row < len(self.exhibits):
            return
        self._clear()
        exhibit = self.exhibits[row]

        heading = QLabel(exhibit.title)
        heading.setWordWrap(True)
        heading.setStyleSheet("font-weight:bold;font-size:14px;")
        self.content_layout.addWidget(heading)

        chip = QLabel(f"{exhibit.kind} · {BASES[exhibit.basis]}")
        chip.setStyleSheet("color:#6f7684;")
        chip.setWordWrap(True)
        self.content_layout.addWidget(chip)
        self.content_layout.addWidget(note_label(exhibit.what, "color:#c8ccd4;"))

        self._pane = self._build_pane(exhibit)
        self.content_layout.addWidget(self._pane, 1)
        self.content_layout.addWidget(
            note_label(f"NOT THIS: {exhibit.not_this}", "color:#d9a441;"))
        # Shown and laid out explicitly. Rebuilding the panel of an already
        # visible window left every new widget unshown at its default 640x480
        # geometry — the panel came up blank in the running application while
        # the very first exhibit, built before the window was shown, was fine.
        for index in range(self.content_layout.count()):
            item = self.content_layout.itemAt(index).widget()
            if item is not None:
                item.show()
        self.content_layout.activate()

    def _build_pane(self, exhibit) -> QWidget:
        if exhibit.kind == "figure":
            return FigurePane(exhibit)
        if exhibit.kind == "chart":
            from .exhibit_plots import build_chart

            return ChartPane(build_chart(exhibit.plot, self.data))
        if exhibit.kind == "simulation":
            from .exhibit_models import Context, SIMULATIONS

            simulation = SIMULATIONS.get(exhibit.simulation)
            if simulation is None:
                return note_label(f"no simulation named "
                                  f"{exhibit.simulation!r}")
            context = Context(structure=self.structure, result=self.data,
                              species=self.species)
            return SimulationPane(simulation, context)
        spec = MODEL_ACTIONS.get(exhibit.action)
        bound = (BoundAction(spec, self.win, result=self.data) if spec
                 else None)
        reason = ("this window was opened without the application window, so "
                  "there is nothing to draw on" if self.win is None else "")
        return ModelPane(exhibit,
                         bound if bound and bound.resolved else None,
                         reason=reason)

    def closeEvent(self, event) -> None:      # noqa: N802 (Qt naming)
        self._clear()
        super().closeEvent(event)
