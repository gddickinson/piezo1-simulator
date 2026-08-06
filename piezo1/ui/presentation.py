"""Full-screen presentation mode and the readouts shown over the viewport.

Two related jobs. **Presentation mode** hides the panels and menu so the 3-D
view fills the screen, for demonstrating or recording. **Readouts** put the
numbers that have been measured into the same frame as the thing they describe,
because a screenshot of a structure with the measurement in a side panel loses
the measurement the moment it is cropped.

Which readouts appear is a choice, so it is a dialog rather than a fixed list —
what belongs on screen when showing a pore differs from a mode animation.
"""

from __future__ import annotations

from PyQt6.QtCore import Qt
from PyQt6.QtWidgets import (QCheckBox, QComboBox, QDialog, QDialogButtonBox,
                             QDoubleSpinBox, QFormLayout, QGroupBox, QLabel,
                             QListWidget, QListWidgetItem, QVBoxLayout)

__all__ = ["PresentationController", "DisplayOptionsDialog", "READOUTS"]

#: Every value the HUD can show, with the label used in the chooser. The key is
#: what :meth:`PresentationController.refresh` writes into the overlay.
READOUTS: list[tuple[str, str]] = [
    ("structure", "Structure and gating state"),
    ("selection", "Current selection"),
    ("dome", "Dome radius of curvature"),
    ("pore", "Pore bottleneck and wetting verdict"),
    ("mode", "Normal mode, symmetry and eigenvalue"),
    ("pockets", "Pocket count and largest volume"),
    ("conservation", "Conservation coverage"),
    ("camera", "Camera distance and scale"),
]


class PresentationController:
    """Owns full-screen state and keeps the HUD readouts current."""

    def __init__(self, window) -> None:
        self.win = window
        self._restore: dict = {}

    # ------------------------------------------------------------ full screen

    @property
    def active(self) -> bool:
        return bool(self._restore)

    def toggle(self) -> None:
        self.leave() if self.active else self.enter()

    def enter(self) -> None:
        """Hide the furniture and fill the screen with the model."""
        if self.active:
            return
        win = self.win
        self._restore = {
            "docks": {k: d.isVisible() for k, d in win.docks.docks.items()},
            "menu": win.menuBar().isVisible(),
            "status": win.statusBar().isVisible(),
        }
        for dock in win.docks.docks.values():
            dock.hide()
        win.menuBar().setVisible(False)
        win.statusBar().setVisible(False)
        win.showFullScreen()
        win.viewport.setFocus()
        self.refresh()

    def leave(self) -> None:
        """Put everything back exactly as it was, including hidden panels."""
        if not self.active:
            return
        win = self.win
        win.showNormal()
        win.menuBar().setVisible(self._restore.get("menu", True))
        win.statusBar().setVisible(self._restore.get("status", True))
        # Restore each panel's *previous* visibility rather than showing them
        # all: a user who had closed a panel before presenting should not find
        # it reopened afterwards.
        for key, was_visible in self._restore.get("docks", {}).items():
            dock = win.docks.docks.get(key)
            if dock is not None:
                dock.setVisible(was_visible)
        self._restore = {}
        self.refresh()

    # --------------------------------------------------------------- readouts

    def refresh(self) -> None:
        """Recompute every enabled readout from current application state."""
        win = self.win
        hud = win.viewport.hud
        hud.clear_readouts()
        hud.structure_name = ""
        if win.structure is None:
            hud.update()
            return

        enabled = set(hud.settings.fields) or {"structure", "selection"}
        record = win.record

        if hud.settings.structure_name:
            hud.structure_name = win.structure.name

        if "structure" in enabled and record is not None:
            hud.set_readout("structure",
                            f"{record.pdb} · {record.state} · {record.gating}")
        if "selection" in enabled and win.selected_residues:
            shown = ", ".join(str(r) for r in win.selected_residues[:6])
            more = ("…" if len(win.selected_residues) > 6 else "")
            hud.set_readout("selection",
                            f"{win.selection_label or 'selection'}: {shown}{more}")

        if "dome" in enabled:
            dome = getattr(win.physics_panel, "dome_result", None)
            if dome is not None:
                hud.set_readout(
                    "dome",
                    f"dome R_c {dome.radius_of_curvature / 10:.1f} nm "
                    f"(published closed 10.2 nm)")

        if "pore" in enabled and win.analysis.pore is not None:
            text = f"pore bottleneck {win.analysis.pore.bottleneck_radius:.2f} Å"
            wetting = win.analysis.hydration
            if wetting is not None and wetting.available:
                text += f" · {wetting.verdict}"
            hud.set_readout("pore", text)

        if "mode" in enabled and win.modes is not None:
            index = getattr(win.physics, "_mode_index", 0)
            symmetry = ""
            if win.modes.symmetry is not None and index < len(win.modes.symmetry):
                symmetry = f" ({win.modes.symmetry[index]})"
            hud.set_readout(
                "mode",
                f"mode {index + 1}{symmetry} · "
                f"eigenvalue {win.modes.eigenvalues[index]:.4f}")

        if "pockets" in enabled and win.analysis.pockets:
            largest = max(p.volume for p in win.analysis.pockets)
            hud.set_readout("pockets",
                            f"{len(win.analysis.pockets)} pockets · "
                            f"largest {largest:.0f} Å³")

        if "conservation" in enabled and "conservation" in win.analysis.scalars:
            values = win.analysis.scalars["conservation"]
            hud.set_readout("conservation",
                            f"conservation over {len(values)} residues")

        if "camera" in enabled and win.viewport.scene is not None:
            per_pixel = hud.world_per_pixel()
            hud.set_readout(
                "camera",
                f"camera {win.viewport.scene.camera.distance:.0f} Å · "
                f"{per_pixel:.2f} Å/px")
        hud.update()


class DisplayOptionsDialog(QDialog):
    """Choose what the overlay shows."""

    def __init__(self, hud, parent=None) -> None:
        super().__init__(parent)
        self.hud = hud
        self.setWindowTitle("Display options")
        self.resize(430, 560)
        layout = QVBoxLayout(self)

        elements = QGroupBox("Overlay elements")
        form = QFormLayout(elements)
        self.scale_check = QCheckBox("Scale bar")
        self.scale_check.setToolTip(
            "A bar of known length, so a screenshot carries its own units.\n"
            "Exact in the plane through the camera pivot, where the molecule is.")
        self.clock_check = QCheckBox("Animation time counter")
        self.clock_check.setToolTip(
            "Elapsed time and frame for animations. A morph shows its position\n"
            "as a fraction, not a time, because a morph is an interpolation\n"
            "rather than a trajectory.")
        self.name_check = QCheckBox("Structure name")
        self.readout_check = QCheckBox("Measured values")
        self.axes_check = QCheckBox("Orientation axes")
        for widget, attribute in ((self.scale_check, "scale_bar"),
                                  (self.clock_check, "clock"),
                                  (self.name_check, "structure_name"),
                                  (self.readout_check, "readouts"),
                                  (self.axes_check, "orientation_axes")):
            widget.setChecked(getattr(hud.settings, attribute))
            form.addRow(widget)
        layout.addWidget(elements)

        placement = QGroupBox("Placement and size")
        pform = QFormLayout(placement)
        self.corner_combo = QComboBox()
        self.corner_combo.addItems(["bottom-left", "bottom-right"])
        self.corner_combo.setCurrentText(hud.settings.corner)
        pform.addRow("Scale bar corner", self.corner_combo)
        self.font_spin = QDoubleSpinBox()
        self.font_spin.setRange(0.6, 3.0)
        self.font_spin.setSingleStep(0.1)
        self.font_spin.setValue(hud.settings.font_scale)
        self.font_spin.setToolTip(
            "Larger text for a projector or a recording")
        pform.addRow("Text size", self.font_spin)
        layout.addWidget(placement)

        layout.addWidget(QLabel("Values to show (tick to display):"))
        self.field_list = QListWidget()
        chosen = set(hud.settings.fields)
        for key, label in READOUTS:
            item = QListWidgetItem(label)
            item.setData(Qt.ItemDataRole.UserRole, key)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if key in chosen
                               else Qt.CheckState.Unchecked)
            self.field_list.addItem(item)
        layout.addWidget(self.field_list, 1)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def apply(self) -> None:
        settings = self.hud.settings
        settings.scale_bar = self.scale_check.isChecked()
        settings.clock = self.clock_check.isChecked()
        settings.structure_name = self.name_check.isChecked()
        settings.readouts = self.readout_check.isChecked()
        settings.orientation_axes = self.axes_check.isChecked()
        settings.corner = self.corner_combo.currentText()
        settings.font_scale = float(self.font_spin.value())
        settings.fields = [
            self.field_list.item(i).data(Qt.ItemDataRole.UserRole)
            for i in range(self.field_list.count())
            if self.field_list.item(i).checkState() == Qt.CheckState.Checked]
        self.hud.update()
