"""Colour the model by electrostatic potential — Guo & MacKinnon 2017, Fig 4c.

The published panel is APBS: a Poisson-Boltzmann solve on a grid with a
dielectric boundary between a low-dielectric protein interior and
high-dielectric solvent. This is not that, and the difference has a visible
consequence rather than a theoretical one — see
:mod:`piezo1.physics.electrostatics`. What is drawn here is the
linear-superposition Debye-Hückel potential from formal charges through a
uniform solvent dielectric, evaluated on the solvent-accessible surface and
averaged onto each atom.

**Three things this had to get right to be worth drawing at all.**

*The scale is fixed.* Red at −5 k_BT/e, white at zero, blue at +5, the
convention and the saturation points of Figure 4c. That is why it uses
:data:`ColorBy.POTENTIAL` rather than the existing ``VALUE`` colouring, which
auto-ranges to the 2nd and 98th percentiles: an auto-ranged potential map
paints an almost-neutral protein in full red and blue and cannot be compared
with a published surface or with the same protein in another state.

*Buried atoms are not coloured as neutral by accident.* An atom with no
accessible surface has no surface potential, and the honest colour for "not
measured" is the same white as "measured and zero" only if the status line says
how many there are. It does.

*It says it is not APBS.* On 6B3R nothing on the surface reaches ±5 k_BT/e,
where the published panel visibly saturates — the expected direction of the
approximation, and the number the status line leads with.
"""

from __future__ import annotations

import numpy as np

from ..render.representations import ColorBy

__all__ = ["ElectrostaticColourController"]


class ElectrostaticColourController:
    """Owns the potential colouring under View -> Colour by electrostatics."""

    def __init__(self, window) -> None:
        self.win = window
        self.result = None
        self._on = False

    @property
    def visible(self) -> bool:
        return self._on

    def reset(self) -> None:
        """Forget the colouring without repainting anything.

        For the structure-replacement path: the view this would repaint is
        about to be discarded, but ``result`` and ``_on`` describe the old
        entry, and a stale result is a status line quoting one structure's
        potential over another.
        """
        self._on = False
        self.result = None

    def show(self, on: bool) -> None:
        window = self.win
        if window.view is None or window.structure is None:
            if on:
                window._set_status("load a structure first")
            return
        if not on:
            self._on = False
            self.result = None
            window.view.color_by = window._current_color()
            window.view.values = None
            window.view.rebuild()
            window.viewport.update()
            return

        # The two mode colourings drive `view.values` through the same slot;
        # leaving one lit while this paints would show a control describing a
        # colour that is not on screen.
        physics = getattr(window, "physics", None)
        if physics is not None:
            for button in ("color_button", "fluctuation_button"):
                physics._untick(button)

        window._set_status("computing the surface potential…")
        try:
            self._compute()
        except Exception as exc:                        # noqa: BLE001
            window._set_status(
                f"electrostatics failed: {type(exc).__name__}: {exc}")
            return

        window.view.values = self.result.atom_potential
        window.view.color_by = ColorBy.POTENTIAL
        window.view.rebuild()
        window.viewport.update()
        self._on = True
        window._set_status(self.status_line())

    def _compute(self) -> None:
        from ..physics.electrostatics import surface_potential

        self.result = surface_potential(self.win.structure)

    # ------------------------------------------------------------- reporting

    def status_line(self) -> str:
        """What must be said whenever the surface is coloured.

        Leads with the fact that it is not APBS, because the picture is the
        most recognisable thing in the paper and nothing about a red-and-blue
        protein distinguishes the two methods.
        """
        result = self.result
        if result is None:
            return "no potential colouring"
        potential = np.asarray(result.potential, dtype=float)
        buried = int(np.isnan(result.atom_potential).sum())
        total = len(result.atom_potential)
        return (
            f"NOT APBS: screened Coulomb (Debye-Huckel) from formal charges "
            f"through a UNIFORM dielectric — no dielectric boundary, no "
            f"ion-exclusion layer, no partial charges, all three "
            f"under-estimating |phi|. Scale fixed at +-{result.scale:g} "
            f"k_BT/e as in Figure 4c, and "
            f"{100 * result.fraction_saturated():.1f}% of the surface reaches "
            f"it where the published panel visibly saturates · "
            f"5-95% span {np.percentile(potential, 5):+.2f} to "
            f"{np.percentile(potential, 95):+.2f}, net charge "
            f"{result.meta['net_charge']:+.0f} e, Debye length "
            f"{result.debye_length:.1f} A · {buried} of {total} atoms have no "
            f"accessible surface and are painted neutral, which is the same "
            f"colour as a measured zero")
