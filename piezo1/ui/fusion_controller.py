"""Drawing the modelled HaloTag fusion on the loaded structure.

Rounds 31 and 32 computed where a C-terminal HaloTag would sit and how likely
each of the three is to carry a dye, but nothing put it on screen — the numbers
were reachable only through ``python -m piezo1.cli fusion``. This draws them.

**Everything here is a model, and it is drawn so as to look like one.** There is
no structure of the fusion. The linker is a straight seam, not a conformation,
and the accessible-volume cloud is shown precisely so that a single position is
not mistaken for a determined one. The channel itself stays in its own
colouring, so what is measured and what is modelled never share a style.

The tag body can be drawn two ways, and the difference is the point:

* **as a sphere** of the tag's radius of gyration — the honest picture of a
  position that is modelled and an orientation that is unknown;
* **as its real fold**, the deposited 6U32 coordinates rigidly placed at the
  same centre. More informative and more dangerous, because a drawn fold looks
  like a determined pose. So the spin about the seam stays a free angle the
  user can turn, the atoms that touch the channel are coloured as such, and the
  status line says how many of the sampled orientations clear it at all.
"""

from __future__ import annotations

import numpy as np

__all__ = ["FusionController", "TAG_COLOR", "SEAM_COLOR", "ENVELOPE_COLOR",
           "DYE_COLOR", "CONTACT_COLOR", "TAG_ATOM_SCALE", "FOLD_STYLES"]

#: Chosen to sit apart from the channel — but only *mostly*, and the difference
#: matters now the fold is drawn. Measured against the chain palette, TAG_COLOR
#: is 0.10 from its orange and DYE_COLOR 0.10 from its red, which is not far;
#: nor is there anywhere to move to, since the eight chain hues plus a dark
#: background leave no free colour that is still visible. So colour is *not*
#: what keeps a modelled tag from being read as experimental structure. The
#: status line is, which is why the fold cannot be drawn without it.
TAG_COLOR = (0.85, 0.55, 0.25)
SEAM_COLOR = (0.95, 0.80, 0.35)
ENVELOPE_COLOR = (0.35, 0.55, 0.70)
DYE_COLOR = (0.90, 0.30, 0.45)
#: Tag atoms inside the channel, so the reported contact count is also visible.
CONTACT_COLOR = (0.95, 0.15, 0.15)

#: Fraction of the van der Waals radius the fold is drawn at. Full radii make a
#: featureless blob and the project's ball-and-stick 0.42 A makes a dot cloud;
#: half fills the fold enough to read its shape. Presentation only.
TAG_ATOM_SCALE = 0.5

NAME = "halotag"

#: How the tag's real fold may be drawn once it is shown at all. Keys other
#: than the default name a :class:`~piezo1.render.representations.Style`, so
#: the fold is restyled by the same machinery as the channel and cannot drift
#: from it. Presentation only: every style is the same rigidly placed 6U32,
#: at the same undetermined spin, with the same status-line caveat — a cartoon
#: of the fold is no more a determined pose than a sphere cloud of it.
FOLD_STYLES = (
    ("atoms", "Atom spheres"),
    ("cartoon", "Cartoon"),
    ("tube", "Ribbon tube"),
    ("backbone", "Backbone trace"),
    ("sticks", "Sticks"),
    ("ball_and_stick", "Ball and stick"),
)


class FusionController:
    """Builds and draws the fusion model for the currently loaded structure."""

    def __init__(self, window) -> None:
        self.win = window
        self.model = None
        self.pose = None
        self.show_envelope = False
        self.show_dyes = False
        self.show_atoms = False
        #: A key from FOLD_STYLES. Applies only while the fold is drawn; the
        #: radius-of-gyration sphere is not a style but a statement, and it
        #: stays a sphere.
        self.fold_style = "atoms"
        #: None means "whichever orientation clears the channel best"; once the
        #: user turns it, an explicit angle they chose.
        self.spin = None
        #: Coordinates to model and measure against, when they are not the
        #: loaded ones. A morph moves the channel without replacing the
        #: structure, and the anchor is a C-alpha of the channel — so a tag
        #: placed against the loaded coordinates would hang in space beside a
        #: flattened dome. ``None`` means "the loaded structure".
        self.host = None

    # ------------------------------------------------------------- lifecycle

    @property
    def visible(self) -> bool:
        return self.model is not None

    @property
    def _host(self):
        return self.win.structure if self.host is None else self.host

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
        # If a morph is loaded, the tag has to follow it. Built here it would
        # sit at the start position while the channel is part-way along the
        # path — the tag and the anchor it is attached to visibly apart.
        morph = getattr(self.win, "morph_controller", None)
        if self.model is not None and morph is not None:
            morph.refresh_fusion()

    def model_for(self, host):
        """The fusion model for a given set of coordinates, without drawing.

        Separate from :meth:`_build` because the morph needs one of these per
        frame and must not touch the status line, the scene or the cached
        ``self.model`` while doing it.
        """
        from ..structure.fusion import build_fusion, load_halotag
        try:
            return build_fusion(host, load_halotag())
        except (FileNotFoundError, ValueError, RuntimeError):
            return None

    def set_frame(self, model, host) -> None:
        """Draw a supplied model against supplied coordinates.

        The morph's entry point. Everything else about the drawing — sphere or
        fold, envelope, dyes, the spin — stays exactly as the user set it, so
        scrubbing the slider cannot silently change what is being shown.
        """
        if model is None:
            return
        self.model = model
        self.host = host
        self._draw()

    def set_envelope(self, on: bool) -> None:
        self.show_envelope = bool(on)
        if self.visible:
            self._draw()

    def set_dyes(self, on: bool) -> None:
        self.show_dyes = bool(on)
        if self.visible:
            self._draw()

    def set_atoms(self, on: bool) -> None:
        """Swap between the radius-of-gyration sphere and the real fold."""
        self.show_atoms = bool(on)
        self.spin = None
        if self.visible:
            self._draw()

    def set_fold_style(self, key: str) -> None:
        """Choose how the fold is drawn — cartoon, sticks, or the atom cloud.

        Stored even when nothing is on screen, so the choice survives toggling
        the fold off and on. It never applies to the radius-of-gyration
        sphere, which is drawn that way to claim exactly what the model
        determined, not as a preference.
        """
        if key not in {k for k, _label in FOLD_STYLES}:
            return
        self.fold_style = key
        if self.visible and self.show_atoms:
            self._draw()
        elif self.visible:
            self.win._set_status(
                "fold style stored — it applies once the tag structure is "
                "shown (View → HaloTag fusion → Show tag structure)")

    def rotate_tags(self) -> None:
        """Turn the tag about the seam by one sampling step.

        This exists to make the undetermined degree of freedom *visible*. A
        caption saying the orientation is arbitrary is easy to skip; watching
        the fold spin while the channel, the anchor and every reported distance
        stay put is not.
        """
        if not self.visible or not self.show_atoms:
            self.win._set_status("show the tag structure first")
            return
        from ..structure.fusion_pose import SPIN_SAMPLES

        current = self.pose.spin if self.pose is not None else 0.0
        self.spin = (current + 2.0 * np.pi / SPIN_SAMPLES) % (2.0 * np.pi)
        self._draw()

    def clear(self) -> None:
        unregister = getattr(self.win, "unregister_pick_feature", None)
        if callable(unregister):
            unregister(NAME)
        scene = self.win.viewport.scene
        if scene is not None:
            for key in list(scene.batches):
                if key.startswith(f"{NAME}:"):
                    scene.remove(key)
        self.model = None
        self.pose = None
        self.host = None
        morph = getattr(self.win, "morph_controller", None)
        if morph is not None:
            morph.refresh_fusion()          # drops the per-frame models
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
        self.host = None
        try:
            # The channel is already in the canonical frame if the alignment
            # option is on; if it is not, the fusion is still built on whatever
            # frame is displayed, so the tags land on the model the user sees.
            self.model = build_fusion(self.win.structure, tag)
        except (ValueError, RuntimeError) as exc:
            self.win._set_status(f"HaloTag model failed: {exc}")
            self.model = None
            return

        distances = self.model.pore_exit_distances()
        self._draw(default_status=(
            f"HaloTag ×{self.model.n_tags} modelled at residue "
            f"{self.model.anchor_residues[0]} — tag centre "
            f"{distances[0]:.1f} nm from the pore exit, "
            f"{self.model.volume.volume:.0f} nm³ accessible. "
            f"MODEL, not a structure."))

    def _draw(self, default_status: str = "") -> None:
        scene = self.win.viewport.scene
        if scene is None or self.model is None:
            return
        for key in list(scene.batches):
            if key.startswith(f"{NAME}:"):
                scene.remove(key)

        centres = np.asarray(self.model.tag_centres, dtype=np.float32)
        radius = float(self.model.meta["tag_radius"])

        status = default_status
        if self.show_atoms:
            seams, status = self._draw_fold(scene)
        else:
            # One sphere of the radius of gyration. The fold is known but its
            # orientation on the channel is not, so a sphere is the shape that
            # claims exactly what the model determined.
            self.pose = None
            bodies = scene.spheres(f"{NAME}:tags")
            bodies.upload(centres,
                          np.full(len(centres), radius, np.float32),
                          np.tile(np.float32(TAG_COLOR), (len(centres), 1)),
                          np.zeros(len(centres), np.float32))
            seams = np.asarray(self.model.seams(), dtype=np.float32)

        # The seam — the part with no experimental support at all.
        starts, ends = seams[:, 0], seams[:, 1]
        seam_batch = scene.cylinders(f"{NAME}:seam")
        colour = np.tile(np.float32(SEAM_COLOR), (len(starts), 1))
        seam_batch.upload(starts, ends, np.full(len(starts), 1.2, np.float32),
                          colour, colour)

        self._register_picking(centres)

        if self.show_envelope:
            self._draw_envelope(scene)
        if self.show_dyes:
            # The dye count is informative, the fold's line carries the caveat.
            # So the count prefixes it rather than replacing it: a drawn fold
            # must never be on screen without the statement that its
            # orientation is undetermined.
            dyes = self._draw_dyes(scene, centres)
            if dyes:
                status = f"{dyes}. {status}" if self.show_atoms else dyes
        if status:
            self.win._set_status(status)
        self.win.viewport.update()

    def _draw_fold(self, scene) -> tuple[np.ndarray, str]:
        """The deposited tag structure at the modelled centres.

        Returns the seams to draw and the status line. The seams now run to
        each tag's own N-terminus rather than to its centre: with the fold
        shown, the linker's far end is a real atom rather than a notional
        middle.
        """
        from ..structure.fusion_pose import pose_for_display

        # `_host`, not the loaded structure: the contact count is what colours
        # tag atoms red, and counting them against a channel that is somewhere
        # else along a morph would report contacts with a shape not on screen.
        pose = pose_for_display(self._host, self.model, spin=self.spin)
        self.pose = pose

        colours = np.tile(np.float32(TAG_COLOR), (pose.n_atoms, 1))
        colours[pose.ligand] = DYE_COLOR
        colours[pose.touching & pose.body] = CONTACT_COLOR

        if self.fold_style == "atoms" or not self._draw_fold_view(scene, pose,
                                                                  colours):
            batch = scene.spheres(f"{NAME}:fold")
            batch.upload(
                pose.coords.reshape(-1, 3).astype(np.float32),
                np.tile(pose.radii * TAG_ATOM_SCALE,
                        pose.n_tags).astype(np.float32),
                np.tile(colours, (pose.n_tags, 1)),
                np.zeros(pose.n_atoms * pose.n_tags, np.float32))

        clear = pose.meta["clear_spins"]
        return np.asarray(pose.seams, dtype=np.float32), (
            f"HaloTag {pose.meta['tag_pdb']} fold at spin "
            f"{np.degrees(pose.spin):.0f}°, {pose.body_contacts} atoms inside "
            f"the channel (red). {clear} of {pose.meta['spins_sampled']} "
            f"orientations clear it. THE SPIN IS UNDETERMINED — this is one "
            f"draw; turn it with View → HaloTag fusion → Turn tag orientation.")

    def _register_picking(self, centres: np.ndarray) -> None:
        """Let clicks identify the drawn tag, saying what it is.

        The describe text leads with MODELLED and, for the fold, repeats that
        the spin is undetermined — a click is an identification, and a tag
        atom identified like a deposited one would be the confident wrong
        answer the status line exists to prevent.
        """
        register = getattr(self.win, "register_pick_feature", None)
        if not callable(register):
            return
        if self.show_atoms and self.pose is not None:
            pose = self.pose
            labels = self._fold_labels(pose.n_atoms)

            def describe(i, labels=labels, n=pose.n_atoms):
                return (f"{labels[i % n]} (tag {i // n + 1}) — MODELLED "
                        f"position of the deposited 6U32 fold; the spin about "
                        f"the linker is UNDETERMINED")
            register(NAME, pose.coords.reshape(-1, 3), describe)
        else:
            def describe(i):
                return (f"HaloTag tag {i + 1} centre — MODELLED position "
                        f"(radius-of-gyration sphere); no structure of the "
                        f"fusion exists")
            register(NAME, centres, describe)

    def _fold_labels(self, n_atoms: int) -> list:
        from ..structure.fusion import load_halotag
        from ..structure.fusion_pose import drawable_mask

        try:
            tag = load_halotag().structure
            base = tag.subset(drawable_mask(tag))
        except (FileNotFoundError, ValueError):
            base = None
        if base is None or base.n_atoms != n_atoms:
            return [f"HaloTag atom {i}" for i in range(n_atoms)]
        return [f"HaloTag {rn}{int(rs)} atom {an}"
                for rn, rs, an in zip(base.res_name, base.res_seq,
                                      base.atom_name, strict=True)]

    def _draw_fold_view(self, scene, pose, colours: np.ndarray) -> bool:
        """The fold in a chosen representation, through the same machinery
        that styles the channel.

        Builds a real :class:`Structure` of the three placed tags — one chain
        per copy, so bonds and cartoon traces stay within a tag — and hands it
        to a `MolecularView` whose `color_override` carries the same per-atom
        colours the sphere cloud uses. The contact atoms stay red and the dye
        stays its own colour in every style, because those colours are the
        visible half of the reported numbers, not decoration.

        Returns False when the placed atoms cannot be matched back to the tag
        file, in which case the caller draws the sphere cloud instead: a fold
        silently missing from the screen is worse than one in the wrong style.
        """
        from ..core.structure import Structure
        from ..render.representations import MolecularView, Style
        from ..structure.fusion import load_halotag
        from ..structure.fusion_pose import drawable_mask

        try:
            style = Style(self.fold_style)
            tag = load_halotag().structure
            base = tag.subset(drawable_mask(tag))
        except (FileNotFoundError, ValueError):
            return False
        if base.n_atoms != pose.n_atoms:
            return False

        fields = {f: np.concatenate([getattr(base, f)] * pose.n_tags)
                  for f in Structure._ARRAY_FIELDS if f != "xyz"}
        # One chain label per copy: cross-chain bonds are skipped and cartoon
        # traces are per chain, so this is what keeps the three tags separate.
        fields["chain"] = np.concatenate(
            [np.full(base.n_atoms, str(i + 1)) for i in range(pose.n_tags)])
        placed = Structure(
            xyz=pose.coords.reshape(-1, 3).astype(np.float32),
            name="halotag-fold", **fields)
        placed._build_residue_index()

        view = MolecularView(
            scene, placed, name=f"{NAME}:fold", style=style,
            color_override=np.tile(colours, (pose.n_tags, 1)))
        # In ribbon styles the dye would vanish with the side chains; the
        # ligand pass keeps it, in its own colour. In atom styles it is
        # already in the atoms batch, and drawing it twice adds nothing.
        view.ligands_as_spheres = style in (Style.CARTOON, Style.TUBE,
                                            Style.BACKBONE)
        view.rebuild()
        return True

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

    def _draw_dyes(self, scene, centres: np.ndarray) -> str:
        """Which tags carry a dye, drawn from the labelling model.

        A single draw from the occupancy distribution, not a prediction that
        these particular tags are the labelled ones — at a saturating protocol
        all three are, and the interest is in the sub-saturating case.

        With the fold drawn the marker goes on the dye 6U32 actually resolves,
        which is where a fluorophore sits relative to the tag; with only the
        sphere it goes at the centre, because that is all the sphere knows.
        """
        from ..analysis.labelling import label_sites

        try:
            result = label_sites(self.model)
        except Exception as exc:
            return f"labelling unavailable: {exc}"
        occupied = np.asarray(result["occupied"], dtype=bool)
        if not occupied.any():
            return ""
        if self.pose is not None and self.pose.ligand.any():
            sites = np.array([self.pose.coords[i][self.pose.ligand].mean(axis=0)
                              for i in range(self.pose.n_tags)],
                             dtype=np.float32)
        else:
            sites = centres
        spots = sites[:len(occupied)][occupied[:len(sites)]]
        batch = scene.spheres(f"{NAME}:dyes")
        batch.upload(spots, np.full(len(spots), 3.0, np.float32),
                     np.tile(np.float32(DYE_COLOR), (len(spots), 1)),
                     np.zeros(len(spots), np.float32))
        return (f"{result['n_dyes']} of {len(occupied)} tags labelled "
                f"(per-site p = {result['p_site']:.3f})")
