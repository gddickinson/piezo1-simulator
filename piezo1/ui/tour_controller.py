"""Runs a tour step: sets the view, triggers the analysis, reports the result.

The controller deliberately calls the *same* controllers the panels use rather
than computing anything itself. A tour that recomputed would be a second
implementation, and a teaching tool that quietly disagreed with the application
it is teaching would be worse than none.
"""

from __future__ import annotations

from PyQt6.QtCore import QTimer

from ..tour import TOUR

__all__ = ["TourController"]


class TourController:
    """Applies a tour step to the window and collects its measurement."""

    def __init__(self, window) -> None:
        self.win = window
        self.results: dict = {}
        self._pending: str = ""

    def run_step(self, index: int) -> None:
        if not (0 <= index < len(TOUR)):
            return
        step = TOUR[index]
        win = self.win

        if step.structure and (win.structure is None
                               or win.structure.name != step.structure):
            win.structure_panel.select(step.structure)
            # Loading is synchronous, but the scene rebuild is not; give Qt a
            # turn before measuring so the view the user sees matches the
            # number they are about to read.
            QTimer.singleShot(120, lambda: self._after_load(step))
            return
        self._after_load(step)

    def _after_load(self, step) -> None:
        win = self.win
        if step.style or step.color_by:
            win.structure_panel.set_state(style=step.style or None,
                                          color_by=step.color_by or None)
        if step.highlight:
            win._highlight(list(step.highlight), f"tour: {step.key}")
        elif step.highlight_group:
            self._highlight_group(step.highlight_group)

        self._pending = step.key
        if step.run:
            self._trigger(step.run)
        self._report(step)

    def _highlight_group(self, needle: str) -> None:
        residues: list[int] = []
        for group in self.win.annotations.residue_groups:
            if needle.lower() in group.label.lower():
                residues.extend(group.residues)
        for domain in self.win.annotations.domains:
            if needle.lower() in domain.name.lower():
                residues.extend(range(domain.start, domain.end + 1))
        if residues:
            self.win._highlight(sorted(set(residues)), f"tour: {needle}")

    # ---------------------------------------------------------- the analyses

    def _trigger(self, what: str) -> None:
        win = self.win
        if what == "dome":
            win.physics.measure_dome()
            dome = getattr(win.physics_panel, "dome_result", None)
            if dome is not None:
                self.results["dome"] = dome
        elif what == "pore":
            win.analysis.compute_pore()
        elif what == "modes":
            if win.modes is None:
                win.physics.compute_modes(
                    {"cutoff": 15.0, "spring": "inverse_square", "n_modes": 20})
            else:
                self.results["modes"] = win.modes
        elif what == "footprint":
            self.results["footprint"] = self._footprint()

    def _footprint(self):
        """The linear-versus-nonlinear comparison at the measured geometry."""
        import numpy as np

        from ..physics.dome import DomeGeometrySummary, DomeModel
        dome = self.results.get("dome")
        if dome is None:
            self.win.physics.measure_dome()
            dome = getattr(self.win.physics_panel, "dome_result", None)
            if dome is None:
                return None
            self.results["dome"] = dome
        model = DomeModel(geometry=DomeGeometrySummary.from_measurement(dome))
        from ..physics.elastica import compare_with_linear
        radius = float(np.sqrt(max(model.geometry.projected_area, 0.0) / np.pi))
        return compare_with_linear(radius, model.contact_slope(), model.membrane)

    # ------------------------------------------------------------- reporting

    def collect(self) -> dict:
        """Pull in anything the panels have computed since the step started."""
        win = self.win
        if win.modes is not None:
            self.results["modes"] = win.modes
        if win.analysis.pore is not None:
            self.results["pore"] = win.analysis.pore
        if win.analysis.hydration is not None:
            self.results["hydration"] = win.analysis.hydration
        dome = getattr(win.physics_panel, "dome_result", None)
        if dome is not None:
            self.results["dome"] = dome
        return self.results

    def _report(self, step) -> None:
        self.win.tour_panel.set_measurement(step.report(self.collect()))
        # Threaded analyses finish after the step returns, so look again.
        if step.run in ("pore", "modes"):
            QTimer.singleShot(2500, lambda: self._late_report(step))

    def _late_report(self, step) -> None:
        if self._pending == step.key:
            self.win.tour_panel.set_measurement(step.report(self.collect()))
