"""Showing one part of the assembly, with the residues that matter on it.

A PIEZO1 trimer fills the viewport with blades. Almost every question is about
one part of it, and this draws that part: the backbone of the chosen component
in whatever representation is selected, and the curated residues on it as
ball-and-stick on top.

**It hides, it does not subset.** The structure object is untouched and every
analysis still runs on the whole assembly — the pore profile, the dome, the
modes. That separation is deliberate: this project's rule is that what is drawn
never decides what is computed, and a component selector that quietly changed
the input to the dome fit would be the most confusing possible way to break it.
The status line says which component is shown for the same reason.

The residues drawn on top come from the same curated groups the analyses read,
so a highlighted residue is one the annotation names and not one this module
chose. Where a component's own domains do not contain a highlighted residue —
the PIP2 lysines sit just outside THU9 — the residue is added to the backbone
rather than left floating.
"""

from __future__ import annotations

import numpy as np

__all__ = ["ComponentController", "HIGHLIGHT_COLOR", "HIGHLIGHT_RADIUS",
           "BOND_RADIUS", "HIGHLIGHT_STYLES"]

NAME = "component"

#: Warm gold, deliberately unlike the domain palette and unlike the ion stream:
#: these are annotation, not structure and not a rate.
HIGHLIGHT_COLOR = (1.0, 0.78, 0.28)

#: Angstrom. Ball-and-stick sized rather than van der Waals, so a side chain
#: reads as a side chain and does not swallow the backbone behind it.
HIGHLIGHT_RADIUS = 0.55
BOND_RADIUS = 0.28

#: Two heavy atoms closer than this in the same residue are bonded, near enough
#: for drawing a side chain. Not a chemistry claim — `analysis.interactions` is
#: where bond geometry is decided.
_BOND_CUTOFF = 1.95

#: How the highlighted residues are drawn. All three stay gold — the colour is
#: what says "annotation, not structure" — and none changes which residues are
#: highlighted, which comes from the curated groups alone.
HIGHLIGHT_STYLES = (
    ("ball_and_stick", "Ball and stick"),
    ("sticks", "Sticks"),
    ("spheres", "Spheres (van der Waals)"),
)


class ComponentController:
    """Draws one named component and the curated residues on it."""

    def __init__(self, window) -> None:
        self.win = window
        self.key = "whole"
        self.selection = None
        #: A key from HIGHLIGHT_STYLES.
        self.highlight_style = "ball_and_stick"

    def set_style(self, key: str) -> None:
        """Restyle the highlighted residues. Changes nothing about which."""
        if key not in {k for k, _label in HIGHLIGHT_STYLES}:
            return
        self.highlight_style = key
        if self.selection is not None:
            self._draw()

    # ------------------------------------------------------------ lifecycle

    def show(self, key: str) -> None:
        """Switch to a component. ``whole`` puts everything back."""
        from ..structure.components import component_masks

        structure = getattr(self.win, "structure", None)
        if structure is None:
            self._set_status("load a structure first")
            return

        self.key = key
        try:
            self.selection = component_masks(structure, key)
        except KeyError as exc:
            self._set_status(str(exc))
            return

        self._apply()
        self._draw()
        self._announce()

    def clear(self) -> None:
        """Remove the drawn residues, leaving the model as it was."""
        scene = getattr(self.win.viewport, "scene", None)
        if scene is not None:
            for name in list(getattr(scene, "batches", {})):
                if name.startswith(NAME):
                    scene.remove(name)

    def refresh(self) -> None:
        """Re-apply after a reload, so a component survives changing entry."""
        if self.key != "whole":
            self.show(self.key)

    # --------------------------------------------------------------- drawing

    def _apply(self) -> None:
        """Restrict the molecular view to the component's residues."""
        view = getattr(self.win, "view", None)
        if view is None or self.selection is None:
            return
        selection = self.selection
        view.set_visible_residues(None if selection.component.is_whole
                                  else selection.residues)

    def _draw(self) -> None:
        self.clear()
        scene = getattr(self.win.viewport, "scene", None)
        structure = getattr(self.win, "structure", None)
        if scene is None or structure is None or self.selection is None:
            return
        mask = self.selection.highlight
        if not mask.any():
            return

        xyz = structure.xyz[mask].astype(np.float32)
        if self.highlight_style == "spheres":
            radii = structure.vdw_radii()[mask].astype(np.float32)
        elif self.highlight_style == "sticks":
            # Joint spheres at the bond radius, so the sticks meet cleanly
            # without reading as balls — the same convention the main view's
            # sticks style follows.
            radii = np.full(len(xyz), BOND_RADIUS, np.float32)
        else:
            radii = np.full(len(xyz), HIGHLIGHT_RADIUS, np.float32)
        scene.spheres(f"{NAME}:atoms").upload(
            xyz, radii,
            np.tile(np.array(HIGHLIGHT_COLOR, np.float32), (len(xyz), 1))
              .reshape(len(xyz), 3))

        if self.highlight_style != "spheres":
            starts, ends = self._bonds(structure, mask)
            colour = np.tile(np.array(HIGHLIGHT_COLOR, np.float32),
                             (len(starts), 1)).reshape(len(starts), 3)
            scene.cylinders(f"{NAME}:bonds").upload(
                starts, ends, np.full(len(starts), BOND_RADIUS, np.float32),
                colour)
        self.win.viewport.update()

    def _bonds(self, structure, mask: np.ndarray):
        """Within-residue heavy-atom bonds, so the side chains read as sticks.

        Within a residue only. A bond drawn between two residues would suggest
        a contact, which is what `interaction_controller` draws with published
        criteria — two overlays claiming the same thing by different rules is
        exactly the duplication this project bans.
        """
        index = np.flatnonzero(mask)
        if len(index) < 2:
            return (np.zeros((0, 3), np.float32),) * 2
        xyz = structure.xyz[index].astype(np.float64)
        residue = structure.res_seq[index]
        chain = structure.chain[index]

        starts, ends = [], []
        for key in set(zip(residue.tolist(), chain.tolist())):
            here = np.flatnonzero((residue == key[0]) & (chain == key[1]))
            if len(here) < 2:
                continue
            points = xyz[here]
            gap = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
            a, b = np.where(np.triu(gap < _BOND_CUTOFF, k=1))
            starts.extend(points[a])
            ends.extend(points[b])
        if not starts:
            return (np.zeros((0, 3), np.float32),) * 2
        return (np.asarray(starts, np.float32), np.asarray(ends, np.float32))

    # ------------------------------------------------------------ reporting

    def _announce(self) -> None:
        if self.selection is None:
            return
        text = self.selection.summary()
        if self.selection.note:
            text += " · " + self.selection.note
        if not self.selection.component.is_whole:
            text += (" · hidden, not removed: every analysis still runs on the "
                     "whole assembly")
        self._set_status(text)

    def _set_status(self, text: str) -> None:
        setter = getattr(self.win, "_set_status", None)
        if callable(setter):
            setter(text)
