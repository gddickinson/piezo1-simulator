"""Draw the contacts the interaction analysis finds, between the atoms it found them between.

The inventory has been available since Round 21 and only ever as a table. A
table is the wrong medium for a contact: what a reader wants to know about
R2456-E2117 is not that it exists but *where* — which protomers it joins, and
whether it sits at the gate, the anchor or the blade. Six hundred rows of
residue numbers cannot answer that and one picture can.

**Each kind gets its own colour, because they are not the same evidence.** A
disulfide is a covalent bond that the map resolved. A salt bridge is two
charged groups within a published cutoff. A hydrophobic contact is two carbons
near each other, which is a much weaker statement.

**Two kinds are off by default and the reason is not that there are many.**
8YEZ yields 9,863 contacts, of which 7,984 are hydrogen bonds — and most of
those are backbone i to i+4, which *are* the secondary structure the cartoon is
already drawing. Turning them on repeats the ribbon in eight thousand green
lines. Hydrophobic contacts are off because they are the weakest claim in the
set. What remains is 270 cylinders, each of which says something the ribbon
does not. Whatever is hidden is counted on the status line rather than dropped.

**Nothing here is a new measurement.** The controller calls
:func:`piezo1.analysis.interactions.detect_interactions` and draws exactly what
it returns: same cutoffs, same geometry, same donor-acceptor rules. If the
picture and the table ever disagree, the picture is wrong.

**Two things the drawing must not imply.** Deposited entries carry no
hydrogens, so every criterion here is heavy-atom based and a drawn hydrogen
bond is an inference from geometry rather than an observed proton. And a
contact is a property of *this* structure in *this* state — a closed entry does
not show the open state's salt bridges, and a residue whose side chain is
unresolved cannot contribute one at all. Both are on the status line.
"""

from __future__ import annotations

import numpy as np

__all__ = ["InteractionController", "KIND_COLORS", "KIND_RADIUS",
           "DEFAULT_KINDS", "NAME"]

NAME = "interactions"

#: One colour per kind, dimmest for the weakest claim. Keyed by the kind strings `detect_interactions` actually emits, not by
#: names invented here. The first version used "hbond" where the analysis says
#: "hydrogen_bond", so 7,984 of them — five sixths of everything found —
#: silently failed to draw while the status line cheerfully reported the rest.
#: `test_ui_interactions` now fails if any emitted kind has no colour.
KIND_COLORS = {
    "disulfide": (0.98, 0.85, 0.20),        # covalent, and resolved
    "salt_bridge": (0.35, 0.75, 1.00),
    "hydrogen_bond": (0.55, 0.95, 0.60),
    "cation_pi": (1.00, 0.55, 0.85),
    "pi_stack": (0.85, 0.60, 1.00),
    "hydrophobic": (0.55, 0.55, 0.60),      # many, and weak
}

#: Cylinder radius per kind, Angstrom. A disulfide is drawn thickest because it
#: is the only one that is a bond rather than a proximity.
KIND_RADIUS = {
    "disulfide": 0.22,
    "salt_bridge": 0.16,
    "hydrogen_bond": 0.12,
    "cation_pi": 0.14,
    "pi_stack": 0.14,
    "hydrophobic": 0.07,
}


#: What is drawn unless the caller says otherwise: the specific contacts, not
#: the numerous ones. On 8YEZ this is 270 cylinders against 9,863.
#:
#: **Hydrogen bonds are off, and not because there are a lot of them.** Most of
#: the 7,984 are backbone i to i+4 — they *are* the secondary structure, which
#: the cartoon representation is already drawing. Adding them puts eight
#: thousand lines on screen that repeat what the ribbon said. Hydrophobic
#: contacts are off because they are the weakest claim in the set: two carbons
#: near each other.
DEFAULT_KINDS = ("disulfide", "salt_bridge", "cation_pi", "pi_stack")


def _family(kind: str) -> str:
    """Collapse any sub-kind onto the colours above.

    The analysis may distinguish parallel from T-shaped pi stacking, which
    matters in the table and would be two near-identical purples on screen.
    Anything that is not recognised is returned unchanged so it shows up as a
    missing colour rather than being quietly folded into another family.
    """
    if kind.startswith("pi_stack"):
        return "pi_stack"
    return kind


class InteractionController:
    """Owns the drawn contacts under View -> Interactions."""

    def __init__(self, window) -> None:
        self.win = window
        self.result = None
        self.kinds: set[str] = set(DEFAULT_KINDS)
        self._names: list[str] = []

    @property
    def visible(self) -> bool:
        return bool(self._names)

    def show(self, on: bool) -> None:
        if not on:
            self.clear()
            return
        if self.win.structure is None or self.win.viewport.scene is None:
            self.win._set_status("load a structure first")
            return
        self._build()

    def clear(self) -> None:
        scene = self.win.viewport.scene
        if scene is not None:
            for name in self._names:
                scene.remove(name)
        self._names = []
        self.win.viewport.update()

    def set_kinds(self, kinds) -> None:
        """Choose which families are drawn, and redraw if anything is on."""
        self.kinds = {k for k in kinds if k in KIND_COLORS}
        if self.visible:
            self._build()

    # ------------------------------------------------------------- building

    def _build(self) -> None:
        from ..analysis.interactions import detect_interactions

        self.clear()
        self.win._set_status("finding contacts…")
        try:
            self.result = detect_interactions(self.win.structure)
        except Exception as exc:
            self.win._set_status(
                f"interactions failed: {type(exc).__name__}: {exc}")
            return
        self._draw()

    def _draw(self) -> None:
        structure = self.win.structure
        scene = self.win.viewport.scene
        by_family: dict[str, list] = {}
        for contact in self.result.interactions:
            family = _family(contact.kind)
            if family in self.kinds:
                by_family.setdefault(family, []).append(contact)

        for family, contacts in by_family.items():
            starts = np.array([structure.xyz[c.atom_i] for c in contacts],
                              dtype=np.float32)
            ends = np.array([structure.xyz[c.atom_j] for c in contacts],
                            dtype=np.float32)
            colour = np.tile(np.array(KIND_COLORS[family], np.float32),
                             (len(contacts), 1))
            radii = np.full(len(contacts), KIND_RADIUS[family], np.float32)
            name = f"{NAME}:{family}"
            batch = scene.cylinders(name)
            batch.upload(starts, ends, radii, colour, colour)
            self._names.append(name)

        self.win.viewport.update()
        self.win._set_status(self.status_line(by_family))

    # ------------------------------------------------------------ reporting

    def counts(self) -> dict:
        """How many of each family the analysis found, drawn or not."""
        out: dict[str, int] = {}
        if self.result is None:
            return out
        for contact in self.result.interactions:
            family = _family(contact.kind)
            out[family] = out.get(family, 0) + 1
        return out

    def status_line(self, drawn: dict | None = None) -> str:
        """What is on screen, and the two things it must not be read as.

        The caveats are not optional and are the same two the tabular view
        carries: no deposited entry has hydrogens, so a hydrogen bond here is
        geometry rather than an observed proton; and these are the contacts of
        this structure in this state, so a closed entry cannot show the open
        state's.
        """
        if self.result is None:
            return "no interactions computed"
        drawn = self.counts() if drawn is None else {
            k: len(v) for k, v in drawn.items()}
        shown = ", ".join(f"{n} {k.replace('_', ' ')}"
                          for k, n in sorted(drawn.items(), key=lambda kv: -kv[1]))
        hidden = sum(v for k, v in self.counts().items() if k not in self.kinds)
        extra = f" · {hidden} more hidden" if hidden else ""
        return (f"contacts drawn: {shown or 'none'}{extra} · heavy-atom "
                f"criteria (no deposited entry has hydrogens, so a hydrogen "
                f"bond here is geometry, not an observed proton) · these are "
                f"the contacts of THIS structure in THIS state")
