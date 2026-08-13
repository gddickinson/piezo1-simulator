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

__all__ = ["PreferencesMixin", "BACKGROUNDS"]

#: Viewport background choices, as (key, label, RGBA). The default is the
#: RenderSettings default, asserted below rather than assumed: two copies of
#: one colour drift, and a "default" that silently differed from a fresh
#: start would make the option lie. Coarse steps on purpose — the useful
#: question is dark room versus manuscript figure, not a colour picker.
BACKGROUNDS = (
    ("midnight", "Midnight (default)", (0.055, 0.063, 0.086, 1.0)),
    ("black", "Black", (0.0, 0.0, 0.0, 1.0)),
    ("slate", "Dark grey", (0.16, 0.17, 0.20, 1.0)),
    ("pearl", "Light grey", (0.82, 0.84, 0.87, 1.0)),
    ("white", "White", (1.0, 1.0, 1.0, 1.0)),
)


def _assert_default_matches() -> None:
    from ..config import RenderSettings
    assert BACKGROUNDS[0][2] == RenderSettings().background, (
        "the 'default' background option differs from RenderSettings — "
        "choosing it would change the picture")


_assert_default_matches()


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


    # ----------------------------------------------------------- appearance

    def background_key(self) -> str:
        key = str(self.settings.value("options/background", "midnight",
                                      type=str))
        return key if any(k == key for k, _l, _c in BACKGROUNDS) else "midnight"

    def _set_background(self, key: str) -> None:
        if not any(k == key for k, _l, _c in BACKGROUNDS):
            return
        self.settings.setValue("options/background", key)
        self._apply_background()
        self._set_status(f"viewport background: {key}. The scale bar and "
                         f"readouts keep their dark halo, so they stay "
                         f"legible on any of these.")

    def _apply_background(self) -> None:
        """Push the stored choice into the render settings.

        The scene reads ``settings.background`` every frame for both the
        clear and the depth-cue fog — the fog *is* the background colour, so
        changing one without the other would haze everything toward the old
        colour. Sharing the settings object is what keeps them one value.
        """
        colour = dict((k, c) for k, _l, c in BACKGROUNDS)[self.background_key()]
        self.viewport.settings.background = tuple(colour)
        self.viewport.update()

    def ui_theme(self) -> str:
        from .theme import THEMES
        key = str(self.settings.value("options/ui_theme", "dark", type=str))
        return key if any(k == key for k, _l in THEMES) else "dark"

    def _set_ui_theme(self, key: str) -> None:
        from .theme import THEMES
        if not any(k == key for k, _l in THEMES):
            return
        self.settings.setValue("options/ui_theme", key)
        self._apply_ui_theme()
        self._set_status(f"interface theme: {key}. The viewport background "
                         f"is a separate option, beside this one.")

    def _apply_ui_theme(self) -> None:
        from PyQt6.QtWidgets import QApplication

        from .theme import apply_theme
        app = QApplication.instance()
        if app is not None:
            apply_theme(app, self.ui_theme())

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
        # The two appearance options act immediately when chosen, so they
        # must also act immediately when forgotten — a reset that left a
        # white viewport behind would not be a reset.
        self._apply_background()
        self._apply_ui_theme()
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

    def _show_topology(self) -> None:
        """Open the monomer topology diagram, creating it on first use.

        The reference passed is the structure's **own** numbering, not mouse:
        the diagram places helices by residue number, so reading a human entry
        with the mouse transmembrane table would shift every one of them by up
        to 26 residues and the picture would look entirely reasonable.
        """
        from .topology_window import TopologyWindow
        if self._topology_window is None:
            self._topology_window = TopologyWindow(self)
            self._topology_window.residues_selected.connect(
                self._highlight_topology_range)
        if self.structure is not None:
            species = (self.record.numbering_species if self.record else "human")
            self._topology_window.set_structure(self.structure, species)
        self._topology_window.show()
        self._topology_window.raise_()
        self._topology_window.activateWindow()

    def _highlight_topology_range(self, lo: int, hi: int,
                                  numbering: str) -> None:
        """Select a residue range from the topology diagram on the model.

        All three protomers, because a topology diagram is of *a* protomer and
        the number it gives is a residue number rather than one copy of it —
        the same rule an annotation click follows.
        """
        self._highlight(list(range(int(lo), int(hi) + 1)),
                        f"topology {lo}-{hi} ({numbering})")

    def _show_parameters(self) -> None:
        """Open the parameter registry editor."""
        from ..parameters import PARAMETERS
        from .parameters_dialog import ParametersDialog

        dialog = ParametersDialog(self._reference_lookup(), self)
        dialog.changed.connect(self._on_parameters_changed)
        dialog.exec()
        self._on_parameters_changed()
        del PARAMETERS

    @staticmethod
    def _reference_lookup() -> dict:
        """Citation key -> formatted reference, for the parameter tooltips."""
        import json

        from ..config import RESOURCE_DIR
        path = RESOURCE_DIR / "references.json"
        if not path.exists():
            return {}
        out = {}
        for entry in json.loads(path.read_text())["references"]:
            bits = [entry.get("authors", ""), entry.get("title", ""),
                    entry.get("journal", ""), str(entry.get("year", ""))]
            out[entry["key"]] = " ".join(b for b in bits if b)
            if entry.get("pmid"):
                out[entry["key"]] += f" PMID {entry['pmid']}"
        return out

    def _on_parameters_changed(self) -> None:
        """Make a non-default parameter set visible wherever a number appears.

        Silent overrides are the failure this whole mechanism exists to
        prevent, so the status bar says so persistently rather than once.
        """
        from ..parameters import PARAMETERS
        if PARAMETERS.modified:
            self._set_status("⚠ " + PARAMETERS.override_summary())
            self.hint_label.setText(
                f"⚠ {len(PARAMETERS.overrides())} non-default parameter(s)")
            self.hint_label.setStyleSheet("color:#f2a65a;")
            self.hint_label.setVisible(True)
        else:
            self._set_status("parameters at their documented defaults")
            self.hint_label.setText(
                "drag rotate · shift+drag pan · wheel zoom · R reset · "
                "space spin · click to identify")
            self.hint_label.setStyleSheet("color:#6f7684;")

    def _start_tour(self) -> None:
        """Show the tour dock and begin at the first step."""
        dock = self.docks.docks.get("tour")
        if dock is not None:
            dock.show()
            dock.raise_()
        self.tour_panel.start()
        self._set_status("guided tour — every number is measured, not quoted")
