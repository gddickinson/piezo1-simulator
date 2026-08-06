"""The application shell.

Wires the viewport to the control panels and owns the loaded model, the
annotation set and the physics results. Long computations run on a worker
thread so the viewport never stalls.
"""

from __future__ import annotations

import time

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QKeySequence
from PyQt6.QtWidgets import (QApplication, QDockWidget, QFileDialog, QLabel,
                             QMainWindow, QMessageBox, QStatusBar)

from ..config import SETTINGS, STRUCTURE_DIR
from ..core.annotations import load_annotations
from ..core.structure import Structure
from ..io.registry import load_registry
from ..render.representations import ColorBy, MolecularView, Style
from .gl_widget import ViewportWidget
from .morph_controller import MorphController
from .physics_controller import PhysicsController
from .panels.annotation_panel import AnnotationPanel
from .panels.physics_panel import PhysicsPanel
from .panels.structure_panel import StructurePanel
from .theme import apply_dark_theme

__all__ = ["MainWindow"]


class MainWindow(QMainWindow):
    """Top-level window."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("PIEZO1 Dynamic Structural Simulator")
        self.resize(1680, 1000)

        self.registry = load_registry()
        self.annotations = load_annotations("human")
        self.structure: Structure | None = None
        self.view: MolecularView | None = None
        self.record = None
        self.modes = None
        self._mode_blocks: list[np.ndarray] = []

        self.viewport = ViewportWidget(SETTINGS.render)
        self.setCentralWidget(self.viewport)
        self.viewport.scene_ready.connect(self._on_scene_ready)
        self.viewport.status.connect(self._set_status)
        self.viewport.atom_picked.connect(self._on_pick)

        self._build_docks()
        self._build_menu()

        bar = QStatusBar()
        self.setStatusBar(bar)
        self.status_label = QLabel("Starting…")
        bar.addWidget(self.status_label, 1)
        self.hint_label = QLabel("drag rotate · shift+drag pan · wheel zoom · "
                                 "R reset · space spin · click to identify")
        self.hint_label.setStyleSheet("color:#6f7684;")
        bar.addPermanentWidget(self.hint_label)

    # ----------------------------------------------------------------- setup

    def _build_docks(self) -> None:
        self.structure_panel = StructurePanel()
        self.structure_panel.structure_requested.connect(self.load_structure)
        self.structure_panel.style_changed.connect(self._set_style)
        self.structure_panel.color_changed.connect(self._set_color)
        self.structure_panel.ligands_toggled.connect(self._set_ligands)
        self.structure_panel.radius_changed.connect(self._set_radius)
        self.structure_panel.spin_toggled.connect(
            lambda on: self.viewport.set_spin(28.0 if on else 0.0))

        self.annotation_panel = AnnotationPanel("human")
        self.annotation_panel.residues_selected.connect(self._highlight)
        self.annotation_panel.focus_requested.connect(self._focus_residues)

        self.physics_panel = PhysicsPanel()
        self.physics = PhysicsController(self)
        self.physics_panel.measure_dome_requested.connect(self.physics.measure_dome)
        self.physics_panel.compute_modes_requested.connect(self.physics.compute_modes)
        self.physics_panel.mode_selected.connect(self.physics.select_mode)
        self.physics_panel.animate_toggled.connect(self.physics.animate_mode)
        self.physics_panel.amplitude_changed.connect(self.physics.set_amplitude)
        self.physics_panel.color_by_mode_requested.connect(self.physics.color_by_mode)
        self.morph_controller = MorphController(self)
        self.physics_panel.morph_requested.connect(self.morph_controller.build)
        self.physics_panel.morph_position_changed.connect(
            self.morph_controller.show_frame)
        self.physics_panel.morph_play_toggled.connect(self.morph_controller.play)

        left = QDockWidget("Model", self)
        left.setWidget(self.structure_panel)
        left.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea
                             | Qt.DockWidgetArea.RightDockWidgetArea)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, left)

        right = QDockWidget("Annotation", self)
        right.setWidget(self.annotation_panel)
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, right)

        physics = QDockWidget("Physics", self)
        physics.setWidget(self.physics_panel)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, physics)

        self.resizeDocks([left, physics], [420, 520], Qt.Orientation.Vertical)
        self.resizeDocks([right], [400], Qt.Orientation.Horizontal)

    def _build_menu(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        act = QAction("&Open structure…", self)
        act.setShortcut(QKeySequence.StandardKey.Open)
        act.triggered.connect(self._open_file)
        file_menu.addAction(act)
        file_menu.addSeparator()
        act = QAction("&Quit", self)
        act.setShortcut(QKeySequence.StandardKey.Quit)
        act.triggered.connect(self.close)
        file_menu.addAction(act)

        view_menu = self.menuBar().addMenu("&View")
        for label, cb in (("Reset camera", lambda: self._reset_camera()),
                          ("Clear highlight", lambda: self._highlight([], ""))):
            a = QAction(label, self)
            a.triggered.connect(cb)
            view_menu.addAction(a)

        help_menu = self.menuBar().addMenu("&Help")
        a = QAction("About", self)
        a.triggered.connect(self._about)
        help_menu.addAction(a)

    # ------------------------------------------------------------- lifecycle

    def _on_scene_ready(self, scene) -> None:
        pairs = [(a.pdb, b.pdb) for a, b in self.registry.morph_pairs()]
        self.physics_panel.set_morph_pairs(pairs)
        rec = self.registry.default()
        if rec is not None:
            self.structure_panel.select(rec.pdb)
        else:
            self._set_status("No structures found — run scripts/fetch_data.py")

    def _set_status(self, text: str) -> None:
        self.status_label.setText(text)

    # ---------------------------------------------------------------- loading

    def load_structure(self, pdb: str) -> None:
        if self.viewport.scene is None:
            return
        rec = self.registry.get(pdb)
        if rec is None or not rec.available:
            self._set_status(f"{pdb} not available locally")
            return
        t0 = time.time()
        try:
            st = Structure.from_file(rec.path)
        except Exception as exc:
            QMessageBox.warning(self, "Load failed", f"{rec.path}\n\n{exc}")
            return

        if self.view is not None:
            self.view.clear()
        self.record = rec
        self.structure = st
        self.modes = None
        self.physics_panel.set_modes(None)
        self.physics.reset()

        self.view = MolecularView(self.viewport.scene, st, name=rec.pdb)
        self.view.set_species(rec.numbering_species)
        self.view.style = self._current_style()
        self.view.color_by = self._current_color()
        self.view.rebuild()

        self.viewport.set_pick_source(st.xyz)
        self._reset_camera()

        modelled = self._modelled_residues(st)
        self.annotation_panel.set_structure_context(rec.pdb, modelled)
        self._mode_blocks = self._protomer_blocks(st)

        self._set_status(
            f"{rec.pdb}: {st.n_atoms:,} atoms, {st.n_residues:,} residues, "
            f"{len(self._mode_blocks)} protomers · loaded in "
            f"{time.time() - t0:.2f} s · {rec.numbering_species} numbering")
        self.viewport.update()

    def _modelled_residues(self, st: Structure) -> set[int]:
        per = []
        for ch in st.chains:
            m = st.mask_ca() & (st.chain == ch)
            if m.sum() > 300:
                per.append(set(st.res_seq[m].tolist()))
        return set.intersection(*per) if per else set()

    def _protomer_blocks(self, st: Structure) -> list[np.ndarray]:
        """Equal-length C-alpha blocks, one per protomer, identically ordered."""
        chains = []
        for ch in st.chains:
            m = st.mask_ca() & (st.chain == ch)
            if m.sum() > 300:
                chains.append((st.xyz[m], st.res_seq[m]))
        if len(chains) < 3:
            return []
        common = set(chains[0][1].tolist())
        for _, seq in chains[1:]:
            common &= set(seq.tolist())
        common_arr = np.array(sorted(common))
        blocks = []
        for xyz, seq in chains[:3]:
            idx = np.searchsorted(seq, common_arr)
            blocks.append(xyz[idx].astype(np.float64))
        return blocks

    def _open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open structure", str(STRUCTURE_DIR),
            "Structures (*.cif *.pdb *.cif.gz *.pdb.gz)")
        if not path:
            return
        st = Structure.from_file(path)
        if self.view is not None:
            self.view.clear()
        self.structure = st
        self.view = MolecularView(self.viewport.scene, st, name=st.name)
        self.view.rebuild()
        self.viewport.set_pick_source(st.xyz)
        self._mode_blocks = self._protomer_blocks(st)
        self._reset_camera()
        self._set_status(f"{st.name}: {st.n_atoms:,} atoms")

    # -------------------------------------------------------------- styling

    def _current_style(self) -> Style:
        from .panels.structure_panel import STYLE_LABELS
        return STYLE_LABELS[self.structure_panel.style_combo.currentIndex()][1]

    def _current_color(self) -> ColorBy:
        from .panels.structure_panel import COLOR_LABELS
        return COLOR_LABELS[self.structure_panel.color_combo.currentIndex()][1]

    def _set_style(self, style: Style) -> None:
        if self.view is None:
            return
        self.view.style = style
        self.view.rebuild()
        self.viewport.update()

    def _set_color(self, color: ColorBy) -> None:
        if self.view is None:
            return
        self.view.color_by = color
        self.view.rebuild()
        self.viewport.update()

    def _set_ligands(self, on: bool) -> None:
        if self.view is None:
            return
        self.view.ligands_as_spheres = on
        self.view.rebuild()
        self.viewport.update()

    def _set_radius(self, scale: float) -> None:
        if self.viewport.scene is not None:
            self.viewport.scene.radius_scale = scale
            self.viewport.update()

    def _reset_camera(self) -> None:
        if self.viewport.scene is None or self.structure is None:
            return
        cam = self.viewport.scene.camera
        cam.orbit(0.0, -0.42)          # look down the three-fold axis a little
        cam.frame(self.structure.xyz)
        self.viewport.update()

    # ------------------------------------------------------------ highlight

    def _highlight(self, residues, label: str) -> None:
        if self.view is None or self.structure is None:
            return
        st = self.structure
        if not residues:
            self.view.highlight = None
        else:
            self.view.highlight = np.isin(
                st.res_seq, np.asarray(list(residues), dtype=np.int32))
        # Highlight is drawn by the sphere shader, so show the selected atoms
        # as spheres on top of whatever representation is active.
        self.viewport.scene.remove(f"{self.view.name}:selection")
        if residues:
            mask = self.view.highlight
            n = int(mask.sum())
            if n:
                batch = self.viewport.scene.spheres(f"{self.view.name}:selection")
                batch.upload(st.xyz[mask], st.vdw_radii()[mask] * 1.05,
                             np.tile(np.array([1.0, 0.83, 0.2], np.float32), (n, 1)),
                             np.ones(n, np.float32))
                self._set_status(f"{label}: {n} atoms highlighted")
            else:
                self._set_status(f"{label}: no atoms — residues not modelled "
                                 f"in {self.record.pdb if self.record else '?'}")
        self.viewport.update()

    def _focus_residues(self, residues) -> None:
        if self.structure is None or self.viewport.scene is None or not residues:
            return
        mask = np.isin(self.structure.res_seq,
                       np.asarray(list(residues), dtype=np.int32))
        if mask.any():
            self.viewport.scene.camera.pivot = self.structure.xyz[mask].mean(axis=0)
            self.viewport.update()

    def _on_pick(self, index: int) -> None:
        if index < 0 or self.structure is None:
            self._set_status("nothing under the cursor")
            return
        st = self.structure
        res = int(st.res_seq[index])
        info = self.annotations.annotate_residue(res)
        bits = [f"{st.res_name[index]}{res} chain {st.chain[index]} "
                f"atom {st.atom_name[index]}"]
        if info["domain"]:
            bits.append(f"domain: {info['domain']}")
        if info["groups"]:
            bits.append("sites: " + ", ".join(info["groups"]))
        for v in info["variants"]:
            bits.append(f"variant {v['label']} ({v['classification']})")
        self._set_status("   ·   ".join(bits))

    def _about(self) -> None:
        QMessageBox.about(
            self, "PIEZO1 Dynamic Structural Simulator",
            "<h3>PIEZO1 Dynamic Structural Simulator</h3>"
            "<p>An interactive, physics-driven model of the PIEZO1 "
            "mechanosensitive ion channel.</p>"
            "<p>Structures from the RCSB PDB; sequence and features from "
            "UniProt; predicted structure from AlphaFold DB. Every annotation "
            "carries its provenance — see the panels for sources.</p>")

    def closeEvent(self, event) -> None:
        self.physics.cleanup()
        if self.viewport.scene is not None:
            self.viewport.makeCurrent()
            self.viewport.scene.release()
        super().closeEvent(event)


def main() -> int:
    import sys
    from .gl_widget import configure_surface_format
    configure_surface_format(SETTINGS.render)
    app = QApplication(sys.argv)
    apply_dark_theme(app)
    win = MainWindow()
    win.show()
    return app.exec()
