"""Draw the pockets the detector finds, as the alpha spheres it found them with.

The pocket list has been a table of residue numbers since Round 13. A residue
list is the wrong description of a cavity: it says what lines the pocket and
nothing about its *shape*, and shape is the whole reason anyone looks for one.
The Yoda1 result — a groove rather than an enclosed cavity — is a statement
about shape that the table could not make.

**Alpha spheres are the detection, so alpha spheres are what is drawn.**
:func:`piezo1.analysis.pockets.find_pockets` clusters spheres that fit between
four atoms and touch none; a pocket *is* its cluster. Drawing a smoothed
surface instead would show a boundary the method never computed.

**They overlap, heavily, and that is not a drawing artefact.** Summing the
sphere volumes overcounts several-fold, which is why `Pocket.volume` is a
Monte-Carlo union rather than a sum. A viewer counting spheres is counting the
sampling, not the cavity.

**The cavity on screen may already be occupied.** Detection runs on protein
atoms with ligands excluded, deliberately: a resolved lipid fills the very
pocket being looked for, and leaving it in makes the site vanish. So a drawn
pocket can sit exactly on top of a lipid that is also drawn, and that is the
detector working as designed rather than a clash.

**A cavity is not a binding site.** `ligands.json` carries a graded site
evidence level for exactly this reason, and only one of six modulators has a
residue-level site at all — at `docking_md`, which is inferred rather than
observed. No deposited PIEZO entry contains a bound modulator. A pocket here is
a geometric statement.
"""

from __future__ import annotations

import numpy as np

__all__ = ["PocketController", "POCKET_COLORS", "NAME", "DEFAULT_TOP"]

NAME = "pockets"

#: One colour per drawn pocket, in rank order. Deliberately unlike the chain
#: palette and unlike the interaction kinds — a pocket is neither.
POCKET_COLORS = (
    (0.95, 0.45, 0.55),
    (0.45, 0.85, 0.65),
    (0.60, 0.55, 0.95),
    (0.95, 0.75, 0.35),
    (0.40, 0.80, 0.90),
    (0.85, 0.55, 0.85),
    (0.70, 0.85, 0.40),
    (0.90, 0.60, 0.40),
)

#: How many of the ranked pockets are drawn unless the caller says otherwise.
#: The detector returns up to 30; drawing all of them fills the protein with
#: spheres and the ranking stops meaning anything.
DEFAULT_TOP = 5


class PocketController:
    """Owns the drawn pockets under View -> Pockets."""

    def __init__(self, window) -> None:
        self.win = window
        self.top = DEFAULT_TOP
        self._names: list[str] = []
        #: Set while waiting for the Analysis panel's pocket run to land.
        self.pending = False

    # ----------------------------------------------------------------- state

    @property
    def visible(self) -> bool:
        return bool(self._names)

    @property
    def pockets(self) -> list:
        """The pockets being drawn — the Analysis controller's own list.

        Read rather than recomputed, so the picture and the table are the same
        objects. A second `find_pockets` call would also be a second ranking,
        and two rankings that disagree is a defect with no symptom.
        """
        return self.win.analysis.pockets

    def show(self, on: bool) -> None:
        self.pending = False
        if not on:
            self.clear()
            return
        if self.win.structure is None or self.win.viewport.scene is None:
            self.win._set_status("load a structure first")
            return
        if not self.pockets:
            self.pending = True
            self.win._set_status("detecting pockets to draw them…")
            self.win.analysis.compute_pockets()
            return
        self._draw()

    def clear(self) -> None:
        scene = self.win.viewport.scene
        if scene is not None:
            for name in self._names:
                scene.remove(name)
        self._names = []
        self.win.viewport.update()

    def refresh(self) -> None:
        """Redraw from the current list. Called when a pocket run finishes."""
        if not (self.pending or self.visible):
            return
        self.pending = False
        if not self.pockets or self.win.viewport.scene is None:
            return
        self._draw()

    # -------------------------------------------------------------- building

    def drawn(self) -> list:
        return list(self.pockets)[:self.top]

    def _draw(self) -> None:
        self.clear()
        scene = self.win.viewport.scene
        for i, pocket in enumerate(self.drawn()):
            centers = np.asarray(pocket.centers, dtype=np.float32)
            radii = np.asarray(pocket.radii, dtype=np.float32)
            colour = np.tile(
                np.array(POCKET_COLORS[i % len(POCKET_COLORS)], np.float32),
                (len(centers), 1))
            name = f"{NAME}:{pocket.index}"
            batch = scene.spheres(name)
            batch.upload(centers, radii, colour)
            self._names.append(name)
        self.win.viewport.update()
        self.win._set_status(self.status_line())

    # ------------------------------------------------------------- reporting

    def status_line(self) -> str:
        """What is drawn, and the two readings it must not invite.

        Neither is optional. Spheres overlap, so counting them is not a
        volume — that is why `Pocket.volume` integrates the union. And a
        geometric cavity is not a binding site: no deposited PIEZO entry
        contains a bound modulator, and the one residue-level site this
        project holds is inferred from docking.
        """
        drawn = self.drawn()
        if not drawn:
            return "no pockets detected"
        total = len(self.pockets)
        spheres = sum(p.n_spheres for p in drawn)
        listed = " · ".join(
            f"#{p.index} {p.volume:.0f} A^3 (buriedness {p.buriedness:.2f})"
            for p in drawn[:3])
        hidden = f" · {total - len(drawn)} lower-ranked not drawn" if total > len(drawn) else ""
        return (f"pockets: {len(drawn)} of {total} drawn as {spheres} alpha "
                f"spheres · {listed}{hidden} · the spheres OVERLAP, so their "
                f"count is the sampling and the volume is their union · a "
                f"cavity is geometry, not a binding site: no deposited PIEZO "
                f"entry holds a bound modulator · ligands were excluded before "
                f"detection, so a drawn pocket may sit on a resolved lipid")
