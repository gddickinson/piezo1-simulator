"""Dock management: creation, layout persistence and reset.

Every panel is a full dock — movable to any edge, floatable into its own
window, closable, and tabbable with its neighbours. Qt gives all of that for
free once the features and allowed areas are set, but *only* if they are set;
the defaults leave a dock stuck in the two areas it was created for.

The part Qt does not give free is getting back. A user who floats four panels
and closes two has no route to the original arrangement, so the default layout
is captured once at startup — before the window is shown, while it is still
pristine — and **Reset layout** restores exactly that.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from PyQt6.QtCore import QByteArray, QSettings, Qt
from PyQt6.QtWidgets import QDockWidget, QScrollArea, QWidget

__all__ = ["DockManager", "DockSpec"]

ALL_AREAS = (Qt.DockWidgetArea.LeftDockWidgetArea
             | Qt.DockWidgetArea.RightDockWidgetArea
             | Qt.DockWidgetArea.TopDockWidgetArea
             | Qt.DockWidgetArea.BottomDockWidgetArea)

FULL_FEATURES = (QDockWidget.DockWidgetFeature.DockWidgetMovable
                 | QDockWidget.DockWidgetFeature.DockWidgetFloatable
                 | QDockWidget.DockWidgetFeature.DockWidgetClosable)


@dataclass
class DockSpec:
    """How one panel is presented."""

    key: str
    title: str
    widget: QWidget
    area: Qt.DockWidgetArea
    tooltip: str = ""
    tabify_with: str = ""
    scrollable: bool = True


@dataclass
class DockManager:
    """Owns the window's docks and its layout state."""

    window: object
    docks: dict[str, QDockWidget] = field(default_factory=dict)
    _default_state: QByteArray | None = None
    _settings_key: str = "layout/state"

    # ----------------------------------------------------------------- build

    def add(self, spec: DockSpec) -> QDockWidget:
        dock = QDockWidget(spec.title, self.window)   # type: ignore[arg-type]
        # objectName is what saveState/restoreState key on. Without it Qt
        # silently declines to restore the dock and the layout comes back
        # half-applied, which looks like a corrupt settings file.
        dock.setObjectName(f"dock_{spec.key}")
        dock.setAllowedAreas(ALL_AREAS)
        dock.setFeatures(FULL_FEATURES)
        dock.setWidget(self._wrap(spec.widget) if spec.scrollable
                       else spec.widget)
        if spec.tooltip:
            dock.setToolTip(spec.tooltip)
        self.window.addDockWidget(spec.area, dock)    # type: ignore[attr-defined]

        if spec.tabify_with and spec.tabify_with in self.docks:
            self.window.tabifyDockWidget(                # type: ignore[attr-defined]
                self.docks[spec.tabify_with], dock)
        self.docks[spec.key] = dock
        return dock

    @staticmethod
    def _wrap(widget: QWidget) -> QScrollArea:
        """Let a dock be shorter than the panel inside it.

        Qt will not shrink a dock below its content's minimum, so without this
        the tallest panel sets a floor on the whole window height and a user on
        a short display cannot reach the bottom of the window at all.
        """
        area = QScrollArea()
        area.setWidget(widget)
        area.setWidgetResizable(True)
        area.setFrameShape(QScrollArea.Shape.NoFrame)
        area.setMinimumWidth(240)
        return area

    # ---------------------------------------------------------------- layout

    def capture_default(self) -> None:
        """Remember the arrangement the application ships with."""
        if self._default_state is None:
            self._default_state = self.window.saveState()  # type: ignore[attr-defined]

    def reset(self) -> None:
        """Return every panel to its shipped place, visible and docked."""
        if self._default_state is None:
            return
        for dock in self.docks.values():
            dock.setFloating(False)
            dock.show()
        self.window.restoreState(self._default_state)      # type: ignore[attr-defined]
        for dock in self.docks.values():
            dock.show()

    def show_all(self) -> None:
        for dock in self.docks.values():
            dock.show()

    def hide_all(self) -> None:
        for dock in self.docks.values():
            dock.hide()

    def float_all(self, floating: bool = True) -> None:
        for dock in self.docks.values():
            dock.setFloating(floating)

    # ------------------------------------------------------------ persistence

    def save(self, settings: QSettings) -> None:
        settings.setValue(self._settings_key, self.window.saveState())  # type: ignore[attr-defined]
        settings.setValue("layout/geometry",
                          self.window.saveGeometry())                   # type: ignore[attr-defined]

    def restore(self, settings: QSettings) -> bool:
        """Reapply a remembered layout. Returns whether anything was applied.

        Geometry is restored first and then clamped, because a layout saved on
        a large external monitor and reopened on a laptop would otherwise put
        the window mostly off-screen — the same failure the fixed startup size
        used to cause, arriving by a different route.
        """
        geometry = settings.value("layout/geometry")
        state = settings.value(self._settings_key)
        applied = False
        if isinstance(geometry, QByteArray) and not geometry.isEmpty():
            applied = bool(self.window.restoreGeometry(geometry))        # type: ignore[attr-defined]
            self._clamp_to_screen()
        if isinstance(state, QByteArray) and not state.isEmpty():
            applied = bool(self.window.restoreState(state)) or applied   # type: ignore[attr-defined]
        return applied

    def _clamp_to_screen(self) -> None:
        window = self.window
        screen = window.screen()                                         # type: ignore[attr-defined]
        if screen is None:
            return
        available = screen.availableGeometry()
        frame = window.frameGeometry()                                   # type: ignore[attr-defined]
        if available.contains(frame):
            return
        width = min(frame.width(), available.width())
        height = min(frame.height(), available.height())
        window.resize(width, height)                                     # type: ignore[attr-defined]
        window.move(                                                     # type: ignore[attr-defined]
            max(available.left(), min(frame.left(), available.right() - width)),
            max(available.top(), min(frame.top(), available.bottom() - height)))

    # ------------------------------------------------------------- menu glue

    def view_actions(self) -> list:
        """One checkable action per dock, wired by Qt to its visibility."""
        actions = []
        for key, dock in self.docks.items():
            action = dock.toggleViewAction()
            action.setToolTip(f"Show or hide the {dock.windowTitle()} panel")
            actions.append((key, action))
        return actions
