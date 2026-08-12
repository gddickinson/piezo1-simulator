"""Session save/load and report export, driven from the File menu.

Kept out of :mod:`piezo1.ui.main_window` because that file is the application
shell and was already close to the project's 500-line limit; this is a distinct
concern with its own failure modes.

A session records **what was being looked at** — structure, style, colouring,
camera, selection, which analyses had been run — and never coordinates or
results. That was decided in `io/session.py` and is worth restating here: a
session that embedded results would go stale the moment any parameter changed,
and would quietly present an old number against a new structure.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from PyQt6.QtWidgets import QFileDialog, QMessageBox

from ..config import PROJECT_ROOT
from ..io.session import Session, load_session, save_session

__all__ = ["SessionController"]


class SessionController:
    """File-menu actions: save, load, and export an analysis report."""

    def __init__(self, window) -> None:
        self.win = window
        self.path: Path | None = None

    def export_scalar(self) -> None:
        """Write the scalar currently colouring the model into a PDB B-factor.

        Exports the **raw residue map**, not `view.values`: the displayed array
        has unmeasured residues filled to the map floor, and a reader in PyMOL
        could not tell those from a genuinely low score. Unscored atoms go out
        with occupancy 0.00 instead.
        """
        from PyQt6.QtWidgets import QFileDialog, QMessageBox

        from ..core.export import write_scalar_pdb

        analysis = getattr(self.win, "analysis", None)
        key = getattr(analysis, "coloured_key", "")
        structure = getattr(self.win, "structure", None)
        if structure is None or not key or key not in analysis.scalars:
            QMessageBox.information(
                self.win, "Nothing to export",
                "Colour the model by a computed value first — Analysis panel, "
                "then pore, pockets, conservation or the response scan. The "
                "export writes that scalar into the B-factor column.")
            return

        path, _ = QFileDialog.getSaveFileName(
            self.win, "Export coloured structure",
            f"{structure.name}_{key}.pdb", "PDB (*.pdb)")
        if not path:
            return
        values = analysis.scalars[key]
        span = max(abs(v) for v in values.values()) or 1.0
        # Use the column properly: two decimals over a 0-1 score is a hundred
        # levels, over 0-100 it is ten thousand. The factor goes in the header.
        scale = 1.0 if span > 50 else 100.0 / span
        try:
            report = write_scalar_pdb(structure, values, path, scale=scale,
                                      name=key)
        except ValueError as exc:
            QMessageBox.warning(self.win, "Export failed", str(exc))
            return
        self.win._set_status(report.summary())


    # -------------------------------------------------------------- sessions

    def _capture(self) -> Session:
        win = self.win
        camera = win.viewport.scene.camera

        # Which analyses had been run, with the parameters they used — enough
        # to re-run them, never the results themselves.
        analyses: dict = {}
        if win.analysis.pore is not None:
            analyses["pore"] = {
                "bottleneck_radius_A": float(win.analysis.pore.bottleneck_radius)}
        if win.analysis.pockets:
            analyses["pockets"] = {"n": len(win.analysis.pockets)}
        for key in win.analysis.scalars:
            analyses[key] = {"n_residues": len(win.analysis.scalars[key])}
        if win.modes is not None:
            analyses["modes"] = {"n_modes": int(win.modes.n_modes)}

        return Session(
            structure=win.structure.name if win.structure else "",
            species=(win.record.numbering_species if win.record else "human"),
            style=win._current_style().value,
            color_by=win._current_color().value,
            show_ligands=win.structure_panel.ligand_check.isChecked(),
            radius_scale=float(win.structure_panel.radius_slider.value()) / 100.0,
            camera_rotation=[float(v) for v in camera.rotation],
            camera_pivot=[float(v) for v in camera.pivot],
            camera_distance=float(camera.distance),
            selected_residues=sorted(win.selected_residues),
            selection_label=win.selection_label,
            analyses=analyses)

    def save(self) -> None:
        if self.win.structure is None:
            self.win._set_status("nothing to save — load a structure first")
            return
        start = str(self.path or (PROJECT_ROOT / "session.json"))
        name, _ = QFileDialog.getSaveFileName(
            self.win, "Save session", start, "Session (*.json)")
        if not name:
            return
        try:
            save_session(self._capture(), Path(name))
        except Exception as exc:
            self._error("Could not save the session", exc)
            return
        self.path = Path(name)
        self.win._set_status(f"session saved to {name}")

    def load(self) -> None:
        start = str(self.path or PROJECT_ROOT)
        name, _ = QFileDialog.getOpenFileName(
            self.win, "Load session", start, "Session (*.json)")
        if not name:
            return
        try:
            session = load_session(Path(name))
        except Exception as exc:
            self._error("Could not read that session", exc)
            return
        self.path = Path(name)
        self.apply(session)

    def apply(self, session: Session) -> None:
        """Restore a session onto the window.

        The structure loads first and synchronously, because style, colouring
        and selection all address it; applying them to whatever happened to be
        loaded before would silently produce a valid-looking wrong view.
        """
        win = self.win
        if session.structure and (win.structure is None
                                  or win.structure.name != session.structure):
            win.load_structure(session.structure)
        if win.structure is None:
            win._set_status(f"session referenced {session.structure}, "
                            "which could not be loaded")
            return

        win.structure_panel.set_state(style=session.style,
                                      color_by=session.color_by,
                                      ligands=session.show_ligands,
                                      radius_scale=session.radius_scale)
        camera = win.viewport.scene.camera
        camera.rotation = np.array(session.camera_rotation, dtype=float)
        camera.pivot = np.array(session.camera_pivot, dtype=float)
        camera.distance = float(session.camera_distance)
        win.viewport.update()

        if session.selected_residues:
            win._highlight(list(session.selected_residues),
                           session.selection_label or "restored selection")
        win._set_status(
            f"session restored — {session.structure}, "
            f"{len(session.selected_residues)} residues selected"
            + (f", analyses recorded: {', '.join(sorted(session.analyses))}"
               if session.analyses else ""))

    # --------------------------------------------------------------- reports

    def export_report(self) -> None:
        win = self.win
        if win.structure is None:
            win._set_status("nothing to report — load a structure first")
            return
        name, selected = QFileDialog.getSaveFileName(
            win, "Export analysis report",
            str(PROJECT_ROOT / f"{win.structure.name}_report.md"),
            "Markdown (*.md);;JSON (*.json)")
        if not name:
            return

        from ..analysis.report import build_report
        analyses = ["dome", "pore"]
        if win.modes is not None:
            analyses.append("modes")
        if win.analysis.hydration is not None:
            analyses.append("hydration")

        win._set_status(f"building report ({', '.join(analyses)})…")
        try:
            report = build_report(win.structure,
                                  species=(win.record.numbering_species
                                           if win.record else "human"),
                                  analyses=analyses)
            path = Path(name)
            if path.suffix == ".json" or "JSON" in (selected or ""):
                report.to_json(path.with_suffix(".json"))
            else:
                report.to_markdown(path.with_suffix(".md"))
        except Exception as exc:
            self._error("Could not build the report", exc)
            return
        win._set_status(f"report written to {name}")

    # ----------------------------------------------------------------- error

    def _error(self, title: str, exc: Exception) -> None:
        QMessageBox.warning(self.win, title, f"{type(exc).__name__}: {exc}")
        self.win._set_status(f"{title.lower()} — {exc}")
