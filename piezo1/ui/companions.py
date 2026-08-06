"""Showing more than one structure at once.

Loading a structure normally **replaces** what is on screen, which is what you
want almost always: two PIEZO1 entries in the same frame sit on top of each
other, and a stale one left behind reads as extra density rather than as a
second model. So multi-structure display is opt-in, and while it is on the
window says which structures are drawn.

This is deliberately not the same thing as :mod:`piezo1.ui.overlay_controller`.
That superposes *one* nominated structure onto the loaded one, reports the RMSD
and can colour the reference by per-residue deviation — it is a measurement.
This is a display setting: several models drawn together in the shared canonical
frame, each in its own colour, with nothing computed between them.

**Companions are display only.** Every analysis runs on the primary structure,
whatever else happens to be drawn, for the same reason hiding a lipid cannot
change a pore profile: what is on screen must never decide what is computed.
"""

from __future__ import annotations

from dataclasses import dataclass

from ..core.structure import Structure
from ..render.representations import ColorBy, MolecularView, Style

__all__ = ["CompanionMixin", "COMPANION_COLORS", "Companion"]

#: Distinguishable colours for companions, in order. The primary structure keeps
#: whatever colouring the user chose; companions are deliberately flat and muted
#: so the primary stays the thing being looked at.
COMPANION_COLORS = [
    (0.90, 0.55, 0.30),      # orange
    (0.45, 0.75, 0.55),      # green
    (0.80, 0.45, 0.70),      # magenta
    (0.50, 0.65, 0.90),      # blue
    (0.85, 0.80, 0.40),      # yellow
]


@dataclass
class Companion:
    """One extra structure drawn alongside the primary one."""

    pdb: str
    structure: Structure
    view: MolecularView
    color: tuple
    species: str


class CompanionMixin:
    """Additional structures drawn alongside the loaded one."""

    def _companions(self) -> dict:
        # Created lazily so the mixin needs nothing from __init__.
        if not hasattr(self, "_companion_map"):
            self._companion_map: dict[str, Companion] = {}
        return self._companion_map

    @property
    def multi_structure(self) -> bool:
        return bool(self.settings.value("options/multi_structure", False,
                                        type=bool))

    def set_multi_structure(self, on: bool) -> None:
        """Turn multi-structure display on or off.

        Turning it *off* drops the companions rather than leaving them drawn
        with no way to name them — a setting that says "one structure" while
        three are on screen would be worse than either state.
        """
        self.settings.setValue("options/multi_structure", bool(on))
        if not on:
            self.clear_companions()
        self._refresh_displayed()

    # ------------------------------------------------------------ management

    def displayed_structures(self) -> list[str]:
        """Everything currently drawn, primary first."""
        primary = [self.record.pdb] if self.record is not None else []
        return primary + list(self._companions())

    def add_companion(self, pdb: str, structure: Structure | None = None,
                      species: str | None = None) -> None:
        """Draw ``pdb`` alongside the primary structure, in the same frame.

        ``structure`` may be passed when the caller already has it framed — as
        :meth:`demote_to_companion` does when the previous primary stays on
        screen — so it is not read and re-aligned a second time.
        """
        if self.viewport.scene is None:
            return
        if self.record is not None and pdb == self.record.pdb:
            self._set_status(f"{pdb} is already the primary structure")
            return
        if pdb in self._companions():
            self._set_status(f"{pdb} is already displayed")
            return

        record = self.registry.get(pdb)
        if structure is None:
            if record is None or not record.available:
                self._set_status(f"{pdb} not available locally")
                return
            try:
                structure = Structure.from_file(record.path)
            except Exception as exc:
                self._set_status(f"{pdb} failed to load: {exc}")
                return
            # Same framing as the primary, or they would not sit on each other.
            structure, _frame = self._standardise(structure, record)
        species = species or (record.numbering_species if record else "human")

        used = {c.color for c in self._companions().values()}
        color = next((c for c in COMPANION_COLORS if c not in used),
                     COMPANION_COLORS[len(self._companions())
                                      % len(COMPANION_COLORS)])

        view = MolecularView(self.viewport.scene, structure, name=f"extra:{pdb}")
        view.set_species(species)
        view.style = Style.BACKBONE
        view.color_by = ColorBy.UNIFORM
        view.uniform_color = color
        view.ligands_as_spheres = False
        view.rebuild()

        self._companions()[pdb] = Companion(
            pdb=pdb, structure=structure, view=view, color=color,
            species=species)
        self._refresh_displayed()
        self.viewport.update()

    def demote_to_companion(self, incoming: str | None = None) -> None:
        """Keep the current primary on screen when another one is loaded.

        This is what makes multi-structure display reachable without a second
        way of opening a file: with the option on, loading simply stops
        discarding what was there. The structure is reused as it stands — it is
        already in the shared frame — rather than read from disk again.

        ``incoming`` is the structure about to become primary. It is needed
        because nothing else prevents a structure being drawn twice: reloading
        the entry already on screen would demote it and then re-add it as the
        primary, and promoting a companion would leave the old copy behind. Both
        put two identical models in the same place, which is invisible on screen
        and wrong in the list of what is displayed.
        """
        if incoming is not None:
            # Promoting something already drawn as a companion: drop that copy.
            self.remove_companion(incoming)

        if self.record is None or self.structure is None:
            return
        pdb, structure = self.record.pdb, self.structure
        species = self.record.numbering_species

        if self.view is not None:
            self.view.clear()
            self.view = None
        self.record = None                   # so add_companion does not refuse

        if pdb == incoming:                  # reloading the same entry
            return
        self.add_companion(pdb, structure=structure, species=species)

    def remove_companion(self, pdb: str) -> None:
        companion = self._companions().pop(pdb, None)
        if companion is None:
            return
        companion.view.clear()
        self._refresh_displayed()
        self.viewport.update()

    def clear_companions(self) -> None:
        for companion in list(self._companions().values()):
            companion.view.clear()
        self._companions().clear()
        self._refresh_displayed()
        if self.viewport.scene is not None:
            self.viewport.update()

    # -------------------------------------------------------------- reporting

    def _refresh_displayed(self) -> None:
        """Tell the user what is on screen, and let them take it off again."""
        shown = self.displayed_structures()
        panel = getattr(self, "structure_panel", None)
        if panel is not None and hasattr(panel, "set_displayed"):
            panel.set_displayed(
                shown,
                {pdb: c.color for pdb, c in self._companions().items()})
        if len(shown) > 1:
            self._set_status(f"displaying {len(shown)} structures: "
                             + ", ".join(shown))
