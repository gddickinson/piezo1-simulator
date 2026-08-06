"""Remembered settings and the menu actions that change them.

A mixin rather than a separate object because every one of these is a small
piece of window behaviour that the menus call directly; routing them through a
controller would add a layer without removing a dependency. Split out of
`main_window.py` to keep it under the project's 500-line limit.

Everything persists through ``QSettings``, so the application reopens the way it
was left — including, deliberately, the panel layout, which is the setting a
user is most likely to have spent time arranging.
"""

from __future__ import annotations

import json

from PyQt6.QtWidgets import QTabWidget

from .help_dialog import HelpDialog, open_document

__all__ = ["PreferencesMixin"]


class PreferencesMixin:
    """Menu handlers for layout, options and help. Mixed into MainWindow."""

    #: What selecting an annotation does to the camera. Default is to leave it
    #: alone: on a 2500-residue trimer the view jumping on every click in a
    #: list is disorienting, and the highlight is usually visible anyway.
    FOCUS_MODES = ("none", "centre", "frame")

    def focus_mode(self) -> str:
        return str(self.settings.value("options/focus_mode", "none", type=str))

    def _set_focus_mode(self, mode: str) -> None:
        self.settings.setValue("options/focus_mode", mode)
        self._set_status({
            "none": "selections will not move the view",
            "centre": "the view will centre on each selection",
            "frame": "the view will centre and zoom on each selection",
        }[mode])


    # ------------------------------------------------------- menu handlers

    def _reset_layout(self) -> None:
        """Put every panel back where it started, including closed ones."""
        self.docks.reset()
        self._set_status("panel layout reset")

    def _show_help(self, topic: str = "") -> None:
        """Open the guide, reusing the window if it is already up."""
        if self._help is None:
            self._help = HelpDialog(self)
        if topic == "shortcuts":
            tabs = self._help.findChild(QTabWidget)
            if tabs is not None:
                tabs.setCurrentIndex(1)       # the shortcuts table
        elif topic:
            self._help.show_topic_named(topic)
        self._help.show()
        self._help.raise_()
        self._help.activateWindow()

    def _open_document(self, path: str) -> None:
        if open_document(path):
            self._set_status(f"opened {path}")
        else:
            self._set_status(f"{path} has not been generated yet")

    def _set_remember_layout(self, on: bool) -> None:
        self.settings.setValue("options/remember_layout", bool(on))
        self._set_status("layout will be remembered on exit" if on
                         else "layout will not be remembered")

    def _set_show_hints(self, on: bool) -> None:
        self.hint_label.setVisible(bool(on))
        self.settings.setValue("options/show_hints", bool(on))

    def _spin_speed(self) -> float:
        return float(self.settings.value("options/spin_speed", 28.0, type=float))

    def _set_spin_speed(self, speed: float) -> None:
        self.settings.setValue("options/spin_speed", float(speed))
        self.viewport.set_spin(speed if self.structure_panel.spin_check.isChecked()
                               else 0.0)
        self._set_status(f"spin speed {speed:.0f} deg/s")

    def _reset_options(self) -> None:
        self.settings.clear()
        self.docks.reset()
        self._size_to_screen()
        self._set_status("options and layout restored to defaults")

    # ------------------------------------------------------- presentation

    def _toggle_fullscreen(self) -> None:
        self.presentation.toggle()
        self._set_status("presentation mode — press F11 or Escape to leave"
                         if self.presentation.active else "windowed")

    def _show_display_options(self) -> None:
        from .presentation import DisplayOptionsDialog
        dialog = DisplayOptionsDialog(self.viewport.hud, self)
        if dialog.exec():
            dialog.apply()
            self.settings.setValue("hud/settings",
                                   json.dumps(self.viewport.hud.settings.as_dict()))
            self.presentation.refresh()
            self._set_status("display options updated")

    def _restore_hud_settings(self) -> None:
        from .hud import HudSettings
        raw = self.settings.value("hud/settings", "", type=str)
        if not raw:
            return
        try:
            self.viewport.hud.settings = HudSettings.from_dict(json.loads(raw))
        except (ValueError, TypeError):
            # A settings file written by an older build should not stop the
            # application starting; fall back to defaults silently.
            pass

    def _show_sequences(self) -> None:
        """Open the sequence window, creating it on first use."""
        from .sequence_window import SequenceWindow
        if self._sequence_window is None:
            self._sequence_window = SequenceWindow(self)
            self._sequence_window.residues_selected.connect(self._highlight)
        if self.structure is not None:
            species = (self.record.numbering_species if self.record else "human")
            self._sequence_window.set_structure(self.structure, species)
        self._sequence_window.show()
        self._sequence_window.raise_()
        self._sequence_window.activateWindow()
