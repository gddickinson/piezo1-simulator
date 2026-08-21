"""Standardised framing, file opening and camera reset.

Split out of :mod:`piezo1.ui.main_window` to keep it under the project's
500-line limit. These belong together: all three decide *where the model sits*
rather than what it is, and all three were involved in the same two faults —
the camera composing a relative orbit on every load, so one structure came back
at a different angle each visit, and deposited frames being unrelated, so two
structures never overlapped.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtWidgets import QFileDialog, QMessageBox

from ..config import STRUCTURE_DIR
from ..core.annotations import load_annotations
from ..core.structure import Structure
from ..render.representations import MolecularView
from ..structure.frame import ALIGNMENT_MODES, Frame
from ..structure.protomers import modelled_residues, protomer_blocks

__all__ = ["AlignmentMixin"]


class _OpenedFile:
    """Stand-in registry record for a file opened outside the catalogue.

    ``_standardise`` needs a numbering species to decide whether a reference
    alignment is meaningful; a file opened from disk has no catalogue entry to
    supply one.
    """

    __slots__ = ("pdb", "numbering_species")

    def __init__(self, pdb: str, numbering_species: str) -> None:
        self.pdb = pdb
        self.numbering_species = numbering_species



class AlignmentMixin:
    """Frame-standardising, file opening and camera reset for the window."""

    def _standardise(self, st: Structure, rec):
        """Put a freshly loaded structure into the chosen frame.

        Deposited entries sit at unrelated orientations, so without this two
        structures never overlap and switching between them looks like the
        molecule jumping. The reference is whatever is currently loaded, which
        makes the *first* structure of a session define the frame and every
        later one line up with it.
        """
        from ..structure.frame import standardise

        mode = getattr(self, "alignment_mode", "canonical")
        if mode == "deposited":
            return st, Frame(rotation=np.eye(3), translation=np.zeros(3),
                             mode="deposited")

        reference, same_numbering = self._alignment_reference, True
        if reference is not None and self._alignment_species is not None:
            same_numbering = (self._alignment_species == rec.numbering_species)
        try:
            aligned, frame = standardise(st, mode=mode, reference=reference,
                                         same_numbering=same_numbering)
        except Exception as exc:                       # never block a load
            self._set_status(f"{rec.pdb}: alignment failed ({exc}) — "
                             f"showing the deposited frame")
            return st, Frame(rotation=np.eye(3), translation=np.zeros(3),
                             mode="deposited")

        # The first successfully aligned structure becomes the reference, so a
        # session settles on one frame instead of drifting with each load.
        if self._alignment_reference is None and not frame.is_identity:
            self._alignment_reference = aligned
            self._alignment_species = rec.numbering_species
        return aligned, frame

    def set_alignment_mode(self, mode: str) -> None:
        """Change the frame and reload, so the change is visible immediately."""
        if mode not in ALIGNMENT_MODES:
            raise ValueError(f"unknown alignment mode {mode!r}")
        self.alignment_mode = mode
        self.settings.setValue("options/alignment_mode", mode)
        self._alignment_reference = None
        self._alignment_species = None
        if self.record is not None:
            self.load_structure(self.record.pdb)

    def _open_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Open structure", str(STRUCTURE_DIR),
            "Structures (*.cif *.pdb *.cif.gz *.pdb.gz)")
        if not path:
            return
        try:
            st = Structure.from_file(path)
        except Exception as exc:
            QMessageBox.warning(self, "Load failed", f"{path}\n\n{exc}")
            return
        self.show_opened_structure(st)

    def show_opened_structure(self, st) -> None:
        """Display a structure opened from outside the catalogue.

        Split from the dialog so the whole path can be driven in a test —
        which is how Round 91 found that this path cleared almost nothing:
        it had a two-entry copy of `load_structure`'s clear list, so every
        overlay it lacked survived opening a file, and its ordering ran
        `overlay.clear()` after the old view was cleared, which is the
        documented resurrection bug alive in a second place.
        """
        # An opened file has no registry record, so it has no declared numbering
        # species. Assume human — every reference alignment then goes through
        # the same check as a catalogued entry rather than skipping it.
        species = self.record.numbering_species if self.record else "human"
        st, frame = self._standardise(st, _OpenedFile(st.name, species))

        # Overlay first, while the old view still exists — its clear ends by
        # rebuilding the primary view to undo deviation colouring, and run
        # after `view.clear()` it would put the old batches straight back.
        # That is the exact resurrection `load_structure` documents; this
        # path had it live until Round 91.
        self.overlay.clear()
        if self.multi_structure:
            self.demote_to_companion(incoming=None)
        # Demotion needs a registry record to name the companion by; a
        # previous *opened* file has none, so the old view can still be live
        # here and must not be left drawn under the new one.
        if self.view is not None:
            self.view.clear()
            self.view = None
        if not self.multi_structure:
            self.clear_companions()
        self.record = None
        self.structure = st
        self.modes = None
        # An opened file is deposited coordinates: any spliced-model or
        # assembly state left from the previous entry would put the amber
        # PART PREDICTED banner over a file it does not describe.
        self.full_length = None
        self.assembly = None
        self._fill_refusal = ""
        self.physics_panel.set_modes(None)
        # The same list `load_structure` clears, not a shorter copy: this
        # path's own list once held two entries, and every overlay it lacked
        # survived opening a file — the dome, the contacts, the morph, all
        # drawn for a structure no longer on screen.
        self._clear_structure_overlays()

        self.view = MolecularView(self.viewport.scene, st, name=st.name)
        self.view.set_species(species)
        self.view.style = self._current_style()
        self.view.color_by = self._current_color()
        self.view.rebuild()

        self.viewport.set_pick_source(st.xyz)
        self._refresh_pick_mask()
        self.structure_panel.set_entities(self.view.entity_map())
        # A component chosen on the previous entry is re-applied, exactly as
        # `load_structure` does: the new view shows everything, and a menu
        # saying one part while the picture shows the whole is worse than
        # losing the selection.
        self.components.refresh()
        self.presentation.refresh()
        self.overlay_panel.set_choices(self.registry.entries)
        if self._sequence_window is not None:
            self._sequence_window.set_structure(st, species)
        # Atom indices are per-structure; measurements from the previous one
        # would silently point at different atoms.
        self.measure_panel.set.clear()
        self._refresh_measurements()
        # The same numbering the view was told, so the panel cannot be reading
        # one protein's annotation onto another's coordinates.
        self.annotations = load_annotations(species)
        self.annotation_panel.set_species(species)
        self.annotation_panel.set_structure_context(st.name, modelled_residues(st))
        self._mode_blocks, self._mode_residues = protomer_blocks(st)
        self._reset_camera()
        self._set_status(
            f"{st.name}: {st.n_atoms:,} atoms, {st.n_residues:,} residues"
            + (f" · {frame.summary()}" if not frame.is_identity else ""))
        self.viewport.update()

    # -------------------------------------------------------------- styling

    def _reset_camera(self) -> None:
        if self.viewport.scene is None or self.structure is None:
            return
        cam = self.viewport.scene.camera
        # Absolute, not relative. `orbit` composes with the current rotation, so
        # calling it here made every structure load add another 0.42 rad of
        # pitch — which is why reloading the same entry came back tilted
        # differently each time.
        cam.set_orientation(0.0, -0.42)   # look down the three-fold axis a little
        cam.frame(self.structure.xyz)
        self.viewport.update()

