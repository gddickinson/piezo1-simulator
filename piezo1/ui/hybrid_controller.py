"""Drawing the full-length model, with the graft impossible to mistake.

`structure/hybrid.py` has existed since Round 65 and nothing could build one
outside a notebook — the same exposure gap Round 58 found for the coupling
score. This puts it under **View → Full-length model**.

**The danger is specific and it is why the drawing is the way it is.** Cryo-EM
resolves roughly residues 570–2521; the distal blade, 569 residues, is
AlphaFold. Drawn in one colour the two are indistinguishable, and a picture of
a complete PIEZO1 trimer is exactly the confident-wrong-picture failure the
Round 50 hazard audit was written for. So:

* the experimental core is drawn in one **flat grey** — deliberately dull, and
  deliberately not any of the chain or domain colourings, so it reads as
  background rather than as a result;
* the graft is coloured by **pLDDT in AlphaFold's own bands**, which is the
  signal a structural biologist already reads as "this is a prediction". Only
  48% of it clears the confident threshold, and at these colours that is
  visible rather than stated;
* the **seam** is marked, because the join is where the model is least
  trustworthy and it is otherwise invisible;
* the status line cannot be omitted. It carries the residue range, the
  confident fraction and the 75 Å by which the two models disagree away from
  the seam — the number a good local fit hides.
"""

from __future__ import annotations

import numpy as np

__all__ = ["HybridController", "EXPERIMENTAL_COLOR", "SEAM_COLOR",
           "SEAM_RADIUS", "PLDDT_BANDS", "HYBRID_STYLES"]

#: How the model may be drawn. Presentation only — every style keeps the
#: grey-versus-pLDDT colouring, the seam marker and the status line, because
#: those are what stop a complete-looking trimer from reading as measured.
#: The ribbon styles run through the C-alphas the model records, so they show
#: the chain's path where the sphere cloud shows its bulk.
HYBRID_STYLES = (
    ("spheres", "Atom spheres"),
    ("tube", "Backbone tube"),
    ("backbone", "Backbone trace"),
)

#: Flat and dull on purpose: the experimental part is context here, and any of
#: the real colourings would compete with the confidence bands that matter.
EXPERIMENTAL_COLOR = (0.45, 0.47, 0.52)
#: The join, marked because it is where the model is weakest.
SEAM_COLOR = (1.0, 0.25, 0.55)
SEAM_RADIUS = 4.5

#: AlphaFold's own bands, reused rather than reinvented so the picture matches
#: what the AFDB viewer shows for the same file — and so a fault in the banding
#: is fixed in one place. There was one: see `plddt_band_colors`.
from ..render.colormaps import PLDDT_BANDS, plddt_band_colors  # noqa: E402

NAME = "hybrid"


class HybridController:
    """Builds and draws the experimental-plus-predicted model."""

    def __init__(self, window) -> None:
        self.win = window
        self.model = None
        #: A key from HYBRID_STYLES.
        self.style = "spheres"

    @property
    def visible(self) -> bool:
        return self.model is not None

    def set_style(self, key: str) -> None:
        """Redraw the model in a chosen representation. Display only."""
        if key not in {k for k, _label in HYBRID_STYLES}:
            return
        self.style = key
        if self.visible:
            self._draw()

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
            for key in list(scene.batches):
                if key.startswith(f"{NAME}:"):
                    scene.remove(key)
            self.win.viewport.update()
        self.model = None

    # --------------------------------------------------------------- building

    def _build(self) -> None:
        from ..structure.hybrid import build_hybrid_model

        self.win._set_status("building the full-length model…")
        try:
            self.model = build_hybrid_model(self.win.structure)
        except FileNotFoundError:
            self.win._set_status(
                "the AlphaFold model is not downloaded — "
                "run python -m piezo1.io.fetch")
            self.model = None
            return
        except (ValueError, RuntimeError) as exc:
            self.win._set_status(f"full-length model failed: {exc}")
            self.model = None
            return
        self._draw()

    def _colors(self) -> np.ndarray:
        """Grey for what was measured, confidence bands for what was not."""
        model = self.model
        colors = np.tile(np.float32(EXPERIMENTAL_COLOR), (len(model.xyz), 1))
        predicted = model.predicted
        colors[predicted] = plddt_band_colors(model.plddt[predicted])
        return colors

    def _draw(self) -> None:
        scene = self.win.viewport.scene
        if scene is None or self.model is None:
            return
        for key in list(scene.batches):
            if key.startswith(f"{NAME}:"):
                scene.remove(key)

        if self.style == "spheres" or not self._draw_ribbon(scene):
            model = self.model
            xyz = np.asarray(model.xyz, dtype=np.float32)
            batch = scene.spheres(f"{NAME}:atoms")
            batch.upload(xyz, np.full(len(xyz), 1.6, np.float32),
                         self._colors(), np.zeros(len(xyz), np.float32))

        seam = self._seam_point()
        if seam is not None:
            marker = scene.spheres(f"{NAME}:seam")
            marker.upload(seam[None, :].astype(np.float32),
                          np.float32([SEAM_RADIUS]),
                          np.float32([SEAM_COLOR]), np.float32([1.0]))

        self.win._set_status(self.status_line())
        self.win.viewport.update()

    def _draw_ribbon(self, scene) -> bool:
        """The model as a backbone ribbon through its recorded C-alphas.

        The colours are the same per-atom array the sphere cloud uses, read at
        the C-alphas — grey stays grey and each pLDDT band keeps its colour,
        so restyling cannot soften the prediction's own uncertainty signal.
        Returns False when the model carries no usable C-alpha record, and the
        caller falls back to the sphere cloud rather than drawing nothing.
        """
        from ..render.geometry_builders import Mesh, build_tube

        model = self.model
        ca = getattr(model, "ca", None)
        if ca is None or int(np.asarray(ca).sum()) < 4:
            return False
        ca = np.asarray(ca, dtype=bool)
        res = model.res_seq[ca]
        order = np.argsort(res, kind="stable")
        xyz = np.asarray(model.xyz, np.float64)[ca][order]
        col = self._colors()[ca][order].astype(np.float64)
        res = res[order]

        # Split at unmodelled gaps, exactly as the main view does, so the
        # ribbon does not draw a straight bar across a missing loop.
        breaks = np.flatnonzero(np.diff(res) > 1) + 1
        mesh = Mesh.empty()
        for seg in np.split(np.arange(len(res)), breaks):
            if len(seg) < 4:
                continue
            if self.style == "tube":
                part = build_tube(xyz[seg], col[seg], radius=0.85, sides=10)
            else:
                part = build_tube(xyz[seg], col[seg], radius=0.35, sides=6,
                                  subdivisions=3)
            mesh = mesh.concat(part)
        if not mesh.n_vertices:
            return False
        scene.mesh(f"{NAME}:ribbon", two_sided=True).upload(
            mesh.positions, mesh.normals, mesh.colors, mesh.indices)
        return True

    def _seam_point(self) -> np.ndarray | None:
        model = self.model
        if model.seam_residue is None:
            return None
        at_seam = model.res_seq == model.seam_residue
        if not at_seam.any():
            return None
        return np.asarray(model.xyz[at_seam].mean(axis=0), dtype=float)

    def status_line(self) -> str:
        """What must be on screen whenever the model is.

        A test enforces that this names the grafted region as predicted: a
        full-length trimer whose 569 predicted residues look experimental is
        the failure this whole feature has to avoid.
        """
        model = self.model
        predicted = model.predicted
        n_pred = int(predicted.sum())
        residues = model.res_seq[predicted]
        confident = (float(model.confident_prediction.sum()) / n_pred
                     if n_pred else float("nan"))
        return (f"Full-length model: residues "
                f"{int(residues.min())}–{int(residues.max())} are PREDICTED "
                f"(AlphaFold, coloured by pLDDT — {confident:.0%} confident), "
                f"the rest experimental (grey). Seam at "
                f"{model.seam_residue} fits to {model.overlap_rmsd:.1f} Å, but "
                f"the two models differ by {model.global_rmsd:.0f} Å overall — "
                f"the graft is placed, NOT validated.")
