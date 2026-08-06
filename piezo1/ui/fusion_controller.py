"""Drawing the modelled HaloTag fusion on the loaded structure.

Rounds 31 and 32 computed where a C-terminal HaloTag would sit and how likely
each of the three is to carry a dye, but nothing put it on screen — the numbers
were reachable only through ``python -m piezo1.cli fusion``. This draws them.

**Everything here is a model, and it is drawn so as to look like one.** There is
no structure of the fusion. The tag body is a sphere of the tag's radius of
gyration, not its fold; the linker is a straight seam, not a conformation; and
the accessible-volume cloud is shown precisely so that a single sphere is not
mistaken for a determined position. The channel itself stays in its own
colouring, so what is measured and what is modelled never share a style.
"""

from __future__ import annotations

import numpy as np

__all__ = ["FusionController", "TAG_COLOR", "SEAM_COLOR", "ENVELOPE_COLOR",
           "DYE_COLOR"]

#: Deliberately unlike any colouring the channel uses, so a tag is never read as
#: part of the experimental structure.
TAG_COLOR = (0.85, 0.55, 0.25)
SEAM_COLOR = (0.95, 0.80, 0.35)
ENVELOPE_COLOR = (0.35, 0.55, 0.70)
DYE_COLOR = (0.90, 0.30, 0.45)

NAME = "halotag"


class FusionController:
    """Builds and draws the fusion model for the currently loaded structure."""

    def __init__(self, window) -> None:
        self.win = window
        self.model = None
        self.show_envelope = False
        self.show_dyes = False

    # ------------------------------------------------------------- lifecycle

    @property
    def visible(self) -> bool:
        return self.model is not None

    def toggle(self, on: bool) -> None:
        self.show(on)

    def show(self, on: bool) -> None:
        if not on:
            self.clear()
            return
        if self.win.structure is None or self.win.viewport.scene is None:
            self.win._set_status("load a structure first")
            return
        self._build()

    def set_envelope(self, on: bool) -> None:
        self.show_envelope = bool(on)
        if self.visible:
            self._draw()

    def set_dyes(self, on: bool) -> None:
        self.show_dyes = bool(on)
        if self.visible:
            self._draw()

    def clear(self) -> None:
        scene = self.win.viewport.scene
        if scene is not None:
            for key in list(scene.batches):
                if key.startswith(f"{NAME}:"):
                    scene.remove(key)
        self.model = None
        if self.win.viewport.scene is not None:
            self.win.viewport.update()

    # --------------------------------------------------------------- building

    def _build(self) -> None:
        from ..structure.fusion import build_fusion, load_halotag

        try:
            tag = load_halotag()
        except FileNotFoundError:
            self.win._set_status(
                "6U32 not downloaded — run python -m piezo1.io.fetch")
            return

        self.win._set_status("modelling the HaloTag fusion…")
        try:
            # The channel is already in the canonical frame if the alignment
            # option is on; if it is not, the fusion is still built on whatever
            # frame is displayed, so the tags land on the model the user sees.
            self.model = build_fusion(self.win.structure, tag)
        except (ValueError, RuntimeError) as exc:
            self.win._set_status(f"HaloTag model failed: {exc}")
            self.model = None
            return

        self._draw()
        distances = self.model.pore_exit_distances()
        self.win._set_status(
            f"HaloTag ×{self.model.n_tags} modelled at residue "
            f"{self.model.anchor_residues[0]} — tag centre "
            f"{distances[0]:.1f} nm from the pore exit, "
            f"{self.model.volume.volume:.0f} nm³ accessible. "
            f"MODEL, not a structure.")

    def _draw(self) -> None:
        scene = self.win.viewport.scene
        if scene is None or self.model is None:
            return
        for key in list(scene.batches):
            if key.startswith(f"{NAME}:"):
                scene.remove(key)

        centres = np.asarray(self.model.tag_centres, dtype=np.float32)
        radius = float(self.model.meta["tag_radius"])

        # The tag body. One sphere of the radius of gyration: the fold is known
        # but its orientation on the channel is not, so drawing the real fold
        # would imply a pose that has not been determined.
        bodies = scene.spheres(f"{NAME}:tags")
        bodies.upload(centres,
                      np.full(len(centres), radius, np.float32),
                      np.tile(np.float32(TAG_COLOR), (len(centres), 1)),
                      np.zeros(len(centres), np.float32))

        # The seam — the part with no experimental support at all.
        seams = np.asarray(self.model.seams(), dtype=np.float32)
        starts, ends = seams[:, 0], seams[:, 1]
        seam_batch = scene.cylinders(f"{NAME}:seam")
        colour = np.tile(np.float32(SEAM_COLOR), (len(starts), 1))
        seam_batch.upload(starts, ends, np.full(len(starts), 1.2, np.float32),
                          colour, colour)

        if self.show_envelope:
            self._draw_envelope(scene)
        if self.show_dyes:
            self._draw_dyes(scene, centres)
        self.win.viewport.update()

    def _draw_envelope(self, scene) -> None:
        """The accessible volume, as a thinned point cloud.

        Thinned because the envelope holds around 30,000 grid points and drawing
        every one buries the channel. The stride is uniform, so the cloud stays
        an honest picture of the region's shape.
        """
        points = np.asarray(self.model.volume.points, dtype=np.float32)
        if len(points) == 0:
            return
        stride = max(1, len(points) // 4000)
        points = points[::stride]

        # One cloud per protomer: the envelope was solved for the first anchor,
        # so the other two are its C3 images.
        from ..structure.superpose import rotation_matrix
        axis = self.model.axis
        clouds = [points]
        for step in (1, 2):
            matrix = rotation_matrix(axis.direction, 2.0 * np.pi * step / 3.0)
            clouds.append(((points - axis.point) @ matrix.T
                           + axis.point).astype(np.float32))
        cloud = np.concatenate(clouds)

        batch = scene.spheres(f"{NAME}:envelope")
        batch.upload(cloud, np.full(len(cloud), 0.6, np.float32),
                     np.tile(np.float32(ENVELOPE_COLOR), (len(cloud), 1)),
                     np.zeros(len(cloud), np.float32))

    def _draw_dyes(self, scene, centres: np.ndarray) -> None:
        """Which tags carry a dye, drawn from the labelling model.

        A single draw from the occupancy distribution, not a prediction that
        these particular tags are the labelled ones — at a saturating protocol
        all three are, and the interest is in the sub-saturating case.
        """
        from ..analysis.labelling import label_sites

        try:
            result = label_sites(self.model)
        except Exception as exc:
            self.win._set_status(f"labelling unavailable: {exc}")
            return
        occupied = np.asarray(result["occupied"], dtype=bool)
        if not occupied.any():
            return
        spots = centres[:len(occupied)][occupied]
        batch = scene.spheres(f"{NAME}:dyes")
        batch.upload(spots, np.full(len(spots), 3.0, np.float32),
                     np.tile(np.float32(DYE_COLOR), (len(spots), 1)),
                     np.zeros(len(spots), np.float32))
        self.win._set_status(
            f"{result['n_dyes']} of {len(occupied)} tags labelled "
            f"(per-site p = {result['p_site']:.3f})")
