"""The application shell.

Wires the viewport to the control panels and owns the loaded model, the
annotation set and the physics results. Long computations run on a worker
thread so the viewport never stalls. Startup and argument parsing live in
:mod:`piezo1.ui.app`.
"""

from __future__ import annotations

import time

import numpy as np
from PyQt6.QtCore import Qt
from PyQt6.QtGui import QAction, QKeySequence
from PyQt6.QtWidgets import (QApplication, QLabel, QMainWindow, QMessageBox,
                             QStatusBar)

from ..config import SETTINGS
from ..core.annotations import annotation_gap, load_annotations
from ..core.structure import Structure
from ..io.registry import load_registry
from ..render.representations import ColorBy, MolecularView, Style
from ..structure.frame import ALIGNMENT_MODES
from .analysis_controller import AnalysisController
from .docks import DockManager, DockSpec
from .gl_widget import ViewportWidget
from .menus import build_menus, make_settings
from .alignment import AlignmentMixin
from .companions import CompanionMixin
from .completeness import CompletenessMixin
from .tabular_analyses import TabularAnalysisMixin
from .fusion_controller import FusionController
from .dome_controller import DomeController
from .hybrid_controller import HybridController
from .interaction_controller import InteractionController
from .component_controller import ComponentController
from .ion_flux_controller import IonFluxController
from .micelle_controller import MicelleController
from .nanodomain_controller import NanodomainController
from .path_controller import AllostericPathController
from .pocket_controller import PocketController
from .planar_membrane_controller import PlanarMembraneController
from .potential_controller import ElectrostaticColourController
from .pore_controller import PoreSurfaceController
from .appearance import AppearanceMixin
from .preferences import PreferencesMixin
from .presentation import PresentationController
from .context_menu import ContextMenuMixin
from .selection import PICK_HINTS, SelectionMixin
from ..structure.protomers import modelled_residues, protomer_blocks
from .morph_controller import MorphController
from .physics_controller import PhysicsController
from . import status_bar
from .session_controller import SessionController
from .overlay_controller import OverlayController
from .panels.analysis_panel import AnalysisPanel
from .panels.overlay_panel import OverlayPanel
from .tour_controller import TourController
from .tour_panel import TourPanel
from .panels.annotation_panel import AnnotationPanel
from .panels.measure_panel import MeasurePanel
from .panels.physics_panel import PhysicsPanel
from .panels.structure_panel import StructurePanel

__all__ = ["MainWindow", "PICK_HINTS"]



class MainWindow(AlignmentMixin, CompanionMixin, CompletenessMixin,
                 TabularAnalysisMixin,
                 AppearanceMixin, PreferencesMixin, SelectionMixin,
                 ContextMenuMixin, QMainWindow):
    """Top-level window."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle("PIEZO1 Dynamic Structural Simulator")
        self._size_to_screen()

        self.settings = make_settings()
        self._help = None
        self._sequence_window = None
        self._topology_window = None
        self.registry = load_registry()
        self.annotations = load_annotations("human")
        self.structure: Structure | None = None
        self.view: MolecularView | None = None
        self.record = None
        self.modes = None
        #: The spliced model on screen, or None when the deposited
        #: coordinates are being shown unaltered. Read by the status
        #: line and by every result window, so the two can never
        #: disagree about which is loaded.
        self.full_length = None
        #: The modelled trimer on screen, or None. Read by the provenance
        #: banner and the status prefix, so neither can disagree with what is
        #: drawn — the same rule `full_length` follows.
        self.assembly = None
        self._fill_refusal = ""
        # How freshly loaded structures are framed, and what they are framed
        # against. The reference is the first structure aligned this session, so
        # everything after it lands in the same place.
        self.alignment_mode = self.settings.value(
            "options/alignment_mode", "canonical", type=str)
        if self.alignment_mode not in ALIGNMENT_MODES:
            self.alignment_mode = "canonical"
        self._alignment_reference: Structure | None = None
        self._alignment_species: str | None = None
        self._mode_blocks: list[np.ndarray] = []
        self._mode_residues: np.ndarray = np.array([], dtype=np.int64)
        self.selected_residues: list[int] = []
        self.selection_label: str = ""
        self._pick_hints = 0

        self.viewport = ViewportWidget(SETTINGS.render)
        self.setCentralWidget(self.viewport)
        self.viewport.scene_ready.connect(self._on_scene_ready)
        self.viewport.status.connect(self._set_status)
        self.viewport.atom_picked.connect(self._on_pick)
        self.viewport.feature_picked.connect(self._on_feature_pick)
        self.viewport.context_requested.connect(self._show_context_menu)

        self._build_docks()
        self.fusion = FusionController(self)
        self.hybrid = HybridController(self)
        self.dome_surface = DomeController(self)
        self.contacts = InteractionController(self)
        self.ion_flux = IonFluxController(self)
        self.components = ComponentController(self)
        # Each of these draws something an analysis already measured, and each
        # reads that analysis's own result object rather than recomputing —
        # so a picture and the panel beside it can never be of different runs.
        self.pore_surface = PoreSurfaceController(self)
        self.pocket_view = PocketController(self)
        self.path = AllostericPathController(self)
        self.nanodomain = NanodomainController(self)
        self.micelle = MicelleController(self)
        self.planar_membrane = PlanarMembraneController(self)
        self.electrostatics = ElectrostaticColourController(self)
        self.session = SessionController(self)
        self.presentation = PresentationController(self)
        self._build_menu()

        # Captured before any saved layout is applied, so Reset always has the
        # arrangement the application ships with rather than whatever the user
        # last left.
        self._restore_hud_settings()
        self.docks.capture_default()
        if self.settings.value("options/remember_layout", True, type=bool):
            self.docks.restore(self.settings)

        bar = QStatusBar()
        self.setStatusBar(bar)
        # A plain QLabel here reports its full text width as a size hint,
        # which QStatusBar passes up as a minimum and QMainWindow honours — so
        # a long caveat widened the window, in the worst cases past the edge of
        # the screen. `StatusMessage` keeps the whole text and displays an
        # elided one; the "..." button opens the scrollable history.
        self.status_label = status_bar.install(self, bar)
        self.hint_label = QLabel("drag rotate · shift+drag pan · wheel zoom · "
                                 "R reset · space spin · click to identify")
        self.hint_label.setStyleSheet("color:#6f7684;")
        self.hint_label.setVisible(
            self.settings.value("options/show_hints", True, type=bool))
        bar.addPermanentWidget(self.hint_label)

    # ----------------------------------------------------------------- setup

    #: The layout was designed at this size; it is a preference, not a demand.
    PREFERRED_SIZE = (1680, 1000)

    def _size_to_screen(self) -> None:
        """Open at the preferred size, or the screen's, whichever is smaller.

        A hard ``resize(1680, 1000)`` puts the title bar off the top of a
        1080p or a scaled laptop display, and on some window managers the
        window then cannot be moved or resized back. Clamping to the available
        geometry — which already excludes the menu bar and dock — keeps every
        edge reachable.
        """
        screen = self.screen() or QApplication.primaryScreen()
        if screen is None:                      # headless; nothing to clamp to
            self.resize(*self.PREFERRED_SIZE)
            return
        available = screen.availableGeometry()
        width = min(self.PREFERRED_SIZE[0], int(available.width() * 0.95))
        height = min(self.PREFERRED_SIZE[1], int(available.height() * 0.95))
        self.resize(width, height)
        self.move(available.left() + (available.width() - width) // 2,
                  available.top() + (available.height() - height) // 2)

    def _build_docks(self) -> None:
        self.structure_panel = StructurePanel()
        self.structure_panel.structure_requested.connect(self.load_structure)
        self.structure_panel.fill_changed.connect(self._reload_for_fill)
        self.structure_panel.style_changed.connect(self._set_style)
        self.structure_panel.color_changed.connect(self._set_color)
        self.structure_panel.ligands_toggled.connect(self._set_ligands)
        self.structure_panel.ligand_style_changed.connect(self._set_ligand_style)
        self.structure_panel.radius_changed.connect(self._set_radius)
        self.structure_panel.entities_changed.connect(self._set_entities)
        self.structure_panel.spin_toggled.connect(
            lambda on: self.viewport.set_spin(self._spin_speed() if on else 0.0))

        self.annotation_panel = AnnotationPanel("human")
        self.annotation_panel.residues_selected.connect(self._highlight)
        self.annotation_panel.focus_requested.connect(self._focus_residues)

        self.measure_panel = MeasurePanel()
        self.measure_panel.measurements_changed.connect(self._refresh_measurements)
        self.measure_panel.status.connect(self._set_status)
        self.measure_panel.mode_changed.connect(self._set_measure_mode)

        self.physics_panel = PhysicsPanel()
        self.physics = PhysicsController(self)
        self.physics_panel.measure_dome_requested.connect(self.physics.measure_dome)
        self.physics_panel.compute_modes_requested.connect(self.physics.compute_modes)
        self.physics_panel.mode_selected.connect(self.physics.select_mode)
        self.physics_panel.animate_toggled.connect(self.physics.animate_mode)
        self.physics_panel.amplitude_changed.connect(self.physics.set_amplitude)
        self.physics_panel.color_by_mode_requested.connect(self.physics.color_by_mode)
        self.physics_panel.color_by_fluctuation_requested.connect(
            self.physics.color_by_fluctuation)
        self.morph_controller = MorphController(self)
        self.physics_panel.morph_requested.connect(self.morph_controller.build)
        self.physics_panel.morph_position_changed.connect(
            self.morph_controller.show_frame)
        self.physics_panel.morph_play_toggled.connect(self.morph_controller.play)

        self.analysis_panel = AnalysisPanel()
        self.analysis = AnalysisController(self)
        self.analysis_panel.pore_requested.connect(self.analysis.compute_pore)
        self.analysis_panel.pockets_requested.connect(
            self.analysis.compute_pockets)
        self.analysis_panel.conservation_requested.connect(
            self.analysis.compute_conservation)
        self.analysis_panel.allostery_requested.connect(
            self.analysis.compute_allostery)
        self.analysis_panel.residues_selected.connect(self._highlight)
        self.analysis_panel.focus_requested.connect(self._focus_residues)
        self.analysis_panel.color_requested.connect(self.analysis.apply_color)
        self.analysis_panel.pore_position_picked.connect(
            self.analysis.focus_pore_position)

        self.overlay_panel = OverlayPanel()
        self.overlay = OverlayController(self)
        self.overlay_panel.overlay_requested.connect(self.overlay.load)
        self.overlay_panel.clear_requested.connect(self.overlay.clear)
        self.overlay_panel.style_changed.connect(self.overlay.set_style)
        self.overlay_panel.visibility_changed.connect(self.overlay.set_visible)
        self.overlay_panel.deviation_colour_requested.connect(
            self.overlay.color_reference_by_deviation)

        self.tour_panel = TourPanel()
        self.tour = TourController(self)
        self.tour_panel.step_requested.connect(self.tour.run_step)

        self.docks = DockManager(self)
        for spec in (
            DockSpec("model", "Model", self.structure_panel,
                     Qt.DockWidgetArea.LeftDockWidgetArea,
                     "Choose a structure and how it is drawn"),
            DockSpec("physics", "Physics", self.physics_panel,
                     Qt.DockWidgetArea.LeftDockWidgetArea,
                     "Dome geometry, normal modes and morphing"),
            DockSpec("annotation", "Annotation", self.annotation_panel,
                     Qt.DockWidgetArea.RightDockWidgetArea,
                     "Domains, functional sites and curated variants"),
            DockSpec("analysis", "Analysis", self.analysis_panel,
                     Qt.DockWidgetArea.RightDockWidgetArea,
                     "Pore, pockets, conservation and allostery",
                     tabify_with="annotation"),
            DockSpec("measure", "Measure", self.measure_panel,
                     Qt.DockWidgetArea.RightDockWidgetArea,
                     "Click-to-measure distances, angles and dihedrals",
                     tabify_with="annotation"),
            DockSpec("overlay", "Overlay", self.overlay_panel,
                     Qt.DockWidgetArea.RightDockWidgetArea,
                     "Superpose a second structure and compare them",
                     tabify_with="annotation"),
            DockSpec("tour", "Guided tour", self.tour_panel,
                     Qt.DockWidgetArea.LeftDockWidgetArea,
                     "Walk the mechanism, with every number measured live",
                     tabify_with="physics"),
        ):
            self.docks.add(spec)
        self.docks.docks["annotation"].raise_()

        # Proportional, not absolute: a fixed 420+520 split needs a 940 px
        # column, which a laptop does not have once menu and status bars are
        # taken.
        usable = max(self.height() - 120, 320)
        self.resizeDocks([self.docks.docks["model"], self.docks.docks["physics"]],
                         [int(usable * 0.45), int(usable * 0.55)],
                         Qt.Orientation.Vertical)
        self.resizeDocks([self.docks.docks["annotation"]],
                         [min(400, max(self.width() // 4, 280))],
                         Qt.Orientation.Horizontal)

    def _build_menu(self) -> None:
        build_menus(self)


    # ------------------------------------------------------------- lifecycle

    def _on_scene_ready(self, scene) -> None:
        pairs = [(a.pdb, b.pdb) for a, b in self.registry.morph_pairs()]
        self.physics_panel.set_morph_pairs(pairs)
        rec = self.registry.default()
        if rec is not None:
            self.structure_panel.select(rec.pdb)
        else:
            self._set_status("No structures found — run scripts/fetch_data.py")

    def _reload_for_fill(self, _mode: str) -> None:
        """Rebuild the current entry at the newly chosen completeness."""
        if self.record is not None:
            self.load_structure(self.record.pdb)

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

        st, self.full_length = self._apply_fill(st)

        st, frame = self._standardise(st, rec)

        # Order matters. `overlay.clear()` finishes by rebuilding the primary
        # view to undo any deviation colouring, so it has to run while that view
        # is still the old one — and *before* the old batches are removed.
        # Clearing first meant the rebuild put them straight back, leaving the
        # previous structure on screen for good. That was survivable only while
        # deposited frames sat 100 Å apart; once everything is framed
        # canonically the two superimpose and the ghost is unmistakable.
        self.overlay.clear()
        if self.multi_structure:
            # Keep what was there, in its own colour, instead of discarding it.
            self.demote_to_companion(incoming=rec.pdb)
        else:
            if self.view is not None:
                self.view.clear()
                self.view = None
            self.clear_companions()

        self.record = rec
        self.structure = st
        self.modes = None
        self.physics_panel.set_modes(None)
        self.physics.reset()
        self.analysis.reset()
        self.fusion.clear()
        # The full-length overlay is a model of the entry being replaced, and
        # was missing from this list: it survived a structure change, drawn
        # (and, now that features answer clicks, pickable) over the new entry.
        self.hybrid.clear()
        # The surface was fitted to the structure being replaced, so
        # leaving it up would draw one model's dome over another's.
        self.dome_surface.clear()
        self.contacts.clear()
        # Same rule for everything else drawn from a measurement: the pore, the
        # pockets, the route and the calcium field all belong to the structure
        # being replaced, and one entry's bottleneck drawn inside another's
        # lumen is a picture nothing on screen would contradict.
        self.pore_surface.clear()
        self.pocket_view.clear()
        self.path.reset()
        self.nanodomain.clear()
        # And the morph, missing from this list until Round 87: it is a stored
        # displacement field plus a stored base array, both of the old entry.
        self.morph_controller.reset()

        self.view = MolecularView(self.viewport.scene, st, name=st.name)
        self.view.set_species(rec.numbering_species)
        self.view.style = self._current_style()
        self.view.color_by = self._current_color()
        self.view.ligand_style = self.structure_panel.current_ligand_style()
        self.view.rebuild()

        self.viewport.set_pick_source(st.xyz)
        self.structure_panel.set_entities(self.view.entity_map())
        self._show_provenance()
        # A component chosen on the previous entry is re-applied to this one,
        # because the view is a *new* MolecularView and would otherwise show
        # everything while the menu still said one part — the selector and the
        # picture disagreeing is worse than losing the selection.
        self.components.refresh()
        self.presentation.refresh()
        self.overlay_panel.set_choices(self.registry.entries, exclude=rec.pdb)
        if self._sequence_window is not None:
            self._sequence_window.set_structure(st, rec.numbering_species)
        if self._topology_window is not None:
            self._topology_window.set_structure(st, rec.numbering_species)
        # Atom indices are per-structure, so measurements taken on the previous
        # one would silently point at different atoms.
        self.measure_panel.set.clear()
        self._refresh_measurements()
        self._reset_camera()

        modelled = modelled_residues(st)
        self.annotation_panel.set_structure_context(rec.pdb, modelled)
        self._mode_blocks, self._mode_residues = protomer_blocks(st)

        # Provenance FIRST and in capitals — a status line is read left to
        # right and truncated on the right. The annotation gap rides with it:
        # grey reads as "no domains here", not "none for this protein".
        gap = annotation_gap(rec.numbering_species)
        self._set_status(
            f"{self._provenance_prefix()}"
            + (f"NO ANNOTATION FOR THIS PROTEIN — {gap}. " if gap else "")
            + f"{st.name}: {st.n_atoms:,} atoms, "
            f"{st.n_residues:,} residues, {len(self._mode_blocks)} protomers "
            f"· loaded in {time.time() - t0:.2f} s · {rec.numbering_species} "
            f"numbering"
            + (f" · {frame.summary()}" if not frame.is_identity else ""))
        self._refresh_displayed()
        self.viewport.update()

    def _about(self) -> None:
        QMessageBox.about(
            self, "PIEZO1 Dynamic Structural Simulator",
            "<h3>PIEZO1 Dynamic Structural Simulator</h3>"
            "<p>An interactive, physics-driven model of the PIEZO1 "
            "mechanosensitive ion channel.</p>"
            "<p>Structures from the RCSB PDB; sequence and features from "
            "UniProt; predicted structure from AlphaFold DB. Every annotation "
            "carries its provenance — see the panels for sources.</p>")

    def keyPressEvent(self, event) -> None:      # noqa: N802 (Qt naming)
        """Escape leaves presentation mode.

        Without it a user who enters full screen with the menu bar hidden has
        no visible way out, which on some window managers means force-quitting.
        """
        if (event.key() == Qt.Key.Key_Escape and self.presentation.active):
            self._toggle_fullscreen()
            return
        super().keyPressEvent(event)

    def closeEvent(self, event) -> None:
        if self.settings.value("options/remember_layout", True, type=bool):
            self.docks.save(self.settings)
        self.physics.cleanup()
        self.analysis.cleanup()
        self.overlay.cleanup()
        if self.viewport.scene is not None:
            self.viewport.makeCurrent()
            self.viewport.scene.release()
        super().closeEvent(event)
