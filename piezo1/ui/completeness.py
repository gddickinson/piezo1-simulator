"""What is drawn, and what it is — the Completeness selector's behaviour.

Split from `main_window.py` at the 500-line limit and along a real seam: these
four methods are one concern. Three of them decide what gets spliced into the
displayed structure, and the fourth is the sentence that must accompany it.

The rule they share is the reason the selector is worth having. Whatever is
chosen here becomes `self.structure`, and every analysis, animation and
measurement then runs on it **without knowing** — which is exactly why the
provenance cannot be optional. A dome radius measured on a part-predicted chain
or on a modelled trimer looks identical to one measured on deposited
coordinates, so the amber banner and the status prefix are not decoration; they
are the only thing distinguishing the two.

Both refusals put the selector **back** to deposited-only rather than loading
something the user did not ask for under a name that says they did. And they
are refusals of different kinds: a full-length graft is refused for want of a
prediction, while a trimer assembly is refused because the entry already has
three protomers — not a failure at all, and a message that says so.
"""

from __future__ import annotations

__all__ = ["CompletenessMixin"]


class CompletenessMixin:
    """Splicing, assembly, and the provenance that has to travel with them."""

    def _show_provenance(self) -> None:
        """Put the amber banner up, or take it down, to match what is loaded."""
        hud = getattr(self.viewport, "hud", None)
        if hud is None:
            return
        model = getattr(self, "full_length", None)
        assembly = getattr(self, "assembly", None)
        if assembly is not None:
            hud.provenance = (f"⚠ MODELLED TRIMER — arrangement taken from "
                              f"{assembly.template}, not measured")
        else:
            hud.provenance = ("" if model is None else
                              f"⚠ {model.n_predicted_residues} residues are "
                              f"ALPHAFOLD PREDICTION, not experiment")
        hud.update()

    def _provenance_prefix(self) -> str:
        """What every status line about the loaded model begins with."""
        assembly = getattr(self, "assembly", None)
        if assembly is not None:
            return (f"MODELLED TRIMER on {assembly.template}"
                    + (" AT THE CORE FLOOR" if assembly.at_floor else "")
                    + f", {assembly.clashes:,} clashes · ")
        model = getattr(self, "full_length", None)
        if model is None:
            refusal = getattr(self, "_fill_refusal", "")
            return f"DEPOSITED ({refusal}) · " if refusal else "DEPOSITED · "
        return (f"PART PREDICTED ({model.n_predicted_residues} residues, "
                f"{model.confident_fraction:.0%} above pLDDT 70) · ")

    def _apply_fill(self, st):
        """Splice in whatever the Completeness selector asks for.

        Returns ``(structure, model_or_None)``. A refusal — a PIEZO2 entry, a
        fragment — puts the selector back to deposited-only and says why,
        rather than loading something the user did not ask for under a name
        that says they did.
        """
        from .panels.structure_panel import FILL_MODES

        self._fill_refusal = ""
        self.assembly = None
        mode = self.structure_panel.current_fill()
        if mode == "none":
            return st, None
        if mode == "trimer":
            return self._assemble_trimer(st)
        from ..structure.full_length import build_full_length
        try:
            model = build_full_length(st, mode)
        except (ValueError, FileNotFoundError) as exc:
            # Kept on the window rather than only set as the status, because
            # the load finishes by writing its own status over this one — and
            # a refusal the user never sees is a selector that silently
            # disagrees with what is on screen.
            self.structure_panel.set_fill("none")
            self._fill_refusal = f"cannot build a full-length model — {exc}"
            return st, None
        label = dict((k, v) for k, v, _ in FILL_MODES)[mode]
        self._set_status(f"{model.structure.name}: {label} — {model.summary()}")
        return model.structure, model

    def _assemble_trimer(self, st):
        """Build the modelled trimer, or put the selector back and say why.

        Kept beside `_apply_fill` rather than inside it because the refusals
        are different in kind: a full-length graft is refused for want of a
        prediction, and this is refused because the entry already *has* three
        protomers, which is not a failure at all.
        """
        from ..structure.assembly import assemble_trimer

        result = assemble_trimer(st)
        if not result.ok:
            self.structure_panel.set_fill("none")
            self._fill_refusal = f"cannot assemble a trimer — {result.refusal}"
            return st, None
        self.assembly = result
        self._set_status(result.summary())
        return result.structure, None
