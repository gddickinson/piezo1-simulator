"""How the model is drawn: style, colouring, atom size and what is shown.

Split out of `main_window.py` at the 500-line limit. These are the handlers the
Model panel drives, and they share one property worth stating: **none of them
changes what is computed.** Hiding the lipids or the auxiliary subunit changes
the picture and nothing else — the analyses always run on the channel
protomers, whatever is currently visible.
"""

from __future__ import annotations

from ..render.representations import ColorBy, Style
from .panels.structure_panel import COLOR_LABELS, STYLE_LABELS

__all__ = ["AppearanceMixin"]


class AppearanceMixin:
    """Style, colour, atom size and entity visibility. Mixed into MainWindow."""

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

    def _set_entities(self, visible) -> None:
        """Show or hide whole categories of content — lipid, detergent, the
        MDFIC auxiliary subunit — without touching the analyses, which always
        use the channel protomers regardless of what is drawn."""
        if self.view is None:
            return
        self.view.visible_entities = frozenset(visible)
        self.view.rebuild()
        self.viewport.update()
        hidden = [k for k in self.view.entity_map().present() if k not in visible]
        self._set_status("showing everything in the file" if not hidden
                         else f"hidden: {', '.join(hidden)}")

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
