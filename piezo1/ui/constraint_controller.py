"""Colour the model by how much evolution has refused to change each residue.

Under View -> Colour by constraint. The values are the ``piezo_genes`` census's
per-residue Jensen-Shannon divergence over 174 genome-backed PIEZO1 orthologues;
the picture is what a per-residue table cannot be, which is the reason to draw
it at all: the census's central claim is that the conserved part *is the pore
machinery*, and one look at a coloured trimer either shows that or it does not.

**Three things this had to get right, and they are the same three the
electrostatic colouring had to get right for different reasons.**

*The scale is fixed at 0 to 1.* Jensen-Shannon divergence is already bounded
and its absolute value is meaningful — 0.8 is high constraint whatever else is
on screen. :data:`ColorBy.VALUE` auto-ranges to the 2nd and 98th percentiles,
which would repaint the same protein differently depending on how much blade
happened to be resolved and would make two entries incomparable. So this uses
the fixed value scale.

*A residue with no score is not painted as unconstrained.* The blade tips are
where alignment coverage is worst and where the census's own claim is that
constraint is *low*; painting "not scored" in the same colour as "scored low"
would turn a coverage hole into a finding. Unscored atoms are held out and the
status line counts them.

*It refuses an entry it cannot read.* The track is human PIEZO1's. A mouse
entry goes through the alignment map; a PIEZO2 entry, PEZO-1, dPIEZO or 6LQI's
splice numbering is refused with the reason rather than coloured by whatever
sits at those numbers in PIEZO1.
"""

from __future__ import annotations

import numpy as np

from ..render.representations import ColorBy

__all__ = ["ConstraintColourController", "SCALE"]

#: The fixed colour range. JSD is bounded in [0, 1] and the census's own
#: domain means run 0.54-0.83, so the whole scale is used and no entry can
#: rescale it.
SCALE = (0.0, 1.0)


class ConstraintColourController:
    """Owns the constraint colouring under View -> Colour by constraint."""

    def __init__(self, window) -> None:
        self.win = window
        self.result = None
        self._on = False

    @property
    def visible(self) -> bool:
        return self._on

    def reset(self) -> None:
        """Forget the colouring without repainting anything.

        For the structure-replacement path: ``result`` describes the entry
        being discarded, and a stale one is a status line quoting one
        structure's coverage over another's picture.
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

        # Every value colouring drives `view.values` through one slot. Leaving
        # another lit would show a control describing a colour not on screen.
        physics = getattr(window, "physics", None)
        if physics is not None:
            for button in ("color_button", "fluctuation_button"):
                physics._untick(button)
        electrostatics = getattr(window, "electrostatics", None)
        if electrostatics is not None and getattr(electrostatics, "visible", False):
            electrostatics.show(False)

        window._set_status("placing the family constraint track…")
        try:
            self._compute()
        except Exception as exc:                            # noqa: BLE001
            window._set_status(
                f"constraint colouring failed: {type(exc).__name__}: {exc}")
            return
        if self.result is None or not self.result:
            reason = getattr(self.result, "reason", "no residue could be scored")
            window._set_status(f"REFUSED: {reason}")
            return

        window.view.values = self.result.per_atom
        window.view.color_by = ColorBy.CONSTRAINT
        window.view.rebuild()
        window.viewport.update()
        self._on = True
        window._set_status(self.status_line())

    def _compute(self) -> None:
        from ..analysis.family_constraint import constraint_on_structure

        self.result = constraint_on_structure(self.win.structure)

    # ------------------------------------------------------------- reporting

    def status_line(self) -> str:
        """What must be said whenever the model is coloured this way.

        Leads with whose numbers these are. A coloured trimer looks exactly
        like something this application measured, and it is not — the picture
        is the join, and the values came from a 194-genome census that ran
        somewhere else.
        """
        result = self.result
        if result is None:
            return "no constraint colouring"
        values = np.asarray(result.per_atom, dtype=float)
        finite = values[~np.isnan(values)]
        unscored = int(np.isnan(values).sum())
        route = ("read directly at human numbering" if not result.converted
                 else "residue numbers carried mouse->human through the "
                      "alignment map, never by an offset")
        return (
            f"NOT MEASURED HERE: per-residue Jensen-Shannon divergence over "
            f"174 genome-backed PIEZO1 orthologues, from the piezo_genes "
            f"census · scale FIXED at {SCALE[0]:g}-{SCALE[1]:g} so two entries "
            f"stay comparable · {result.n_residues_scored} of "
            f"{result.n_residues} residues scored "
            f"({result.coverage:.0%}), {route} · median "
            f"{np.median(finite):.2f}, 5-95% span {np.percentile(finite, 5):.2f} "
            f"to {np.percentile(finite, 95):.2f} · {unscored} atoms are "
            f"UNSCORED and held out of the colouring rather than painted at "
            f"zero, because the blade tips are where coverage is worst and "
            f"where low constraint is the claim being made")
