"""Drives the curved-to-flattened conformational morph.

Split out of :mod:`piezo1.ui.main_window` to keep both files readable. The
controller owns the trajectory, the atom-to-site mapping and the playback
animation; the window owns everything else.

A morph is drawn as a **displacement field**: :meth:`MorphController.show_frame`
adds ``frame - frames[0]`` to the coordinates already on screen, rather than
replacing them. That is what carries side chains, ligands and everything else
the C-alpha path does not describe. It also means the path and the picture must
live in the same frame and share the same site indexing, and both of those were
wrong:

* the endpoints were read from disk, in the deposited frame, while the viewport
  shows the canonical one — 180 degrees apart on 7WLT, so every atom was pushed
  the wrong way and the flattened endpoint landed **36 A** from 7WLU, further
  than the 19.7 A change itself;
* atoms at a residue outside the shared basis were tied to a site by *residue
  number*, which put each bound lipid on a C-alpha a median of **64.8 A** away
  (the nearest site is 6.2 A) and flung it across the model with the CTD.

Both are now closed by construction: the path is built from the coordinates
that are displayed, and anything without a site of its own takes the nearest
one in space.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtWidgets import QApplication

from ..core.structure import Structure
from ..physics.anm import ANM

__all__ = ["MorphController"]


class MorphController:
    """Builds and plays back a morph between two experimental endpoints."""

    def __init__(self, window) -> None:
        self.win = window
        self.trajectory = None
        self.residues = None
        self._map = None
        self._base = None
        self._start = None
        self._phase = 0.0
        self._fraction = 0.0
        #: One solved HaloTag model per frame, or None when the tag is not
        #: shown. Solved rather than interpolated — see `refresh_fusion`.
        self._fusion_frames = None
        self._fusion_clashing = 0
        #: Tag-to-pore-exit on the deposited *end* entry, so the far end of the
        #: path can be quoted beside the thing it is standing in for.
        self._fusion_endpoint = None

    def reset(self) -> None:
        """Drop a path built on a structure that is no longer on screen.

        :meth:`show_frame` adds a stored displacement field to a stored base
        array, both captured when the morph was built. Loading another entry
        replaces the view and its atoms while the slider stays enabled, so
        dragging it would push one structure's motion onto another's
        coordinates — or, when the two differ in atom count, fail a long way
        from the cause.
        """
        panel = self.win.physics_panel
        if panel.morph_play.isChecked():
            panel.morph_play.blockSignals(True)
            panel.morph_play.setChecked(False)
            panel.morph_play.blockSignals(False)
            self.play(False)
        self.trajectory = None
        self.residues = None
        self._map = None
        self._base = None
        self._start = None
        self._phase = 0.0
        self._fraction = 0.0
        self._fusion_frames = None
        self._fusion_clashing = 0
        self._fusion_endpoint = None
        panel.morph_slider.blockSignals(True)
        panel.morph_slider.setValue(0)
        panel.morph_slider.blockSignals(False)
        panel.set_morph(None)

    def build(self, params: dict) -> None:
        """Interpolate between two experimental endpoints and load the path."""
        from ..structure.morph import morph, prepare_endpoints

        start_rec = self.win.registry.get(params["start"])
        end_rec = self.win.registry.get(params["end"])
        if start_rec is None or end_rec is None:
            self.win._set_status("morph endpoints not available")
            return
        self.win._set_status(f"building {params['method']} morph "
                         f"{start_rec.pdb} → {end_rec.pdb}…")
        QApplication.processEvents()
        # Captured before anything is cleared: loading the start entry runs
        # `load_structure`, which clears the fusion, so by the time the path
        # exists there is no way left to ask whether the tag was being shown.
        wanted_tag = self.win.fusion.visible
        self.reset()

        # The start structure goes on screen *first*, and the path is then
        # built from the coordinates that are on screen. Reading it from disk
        # here instead gave a displacement field in the deposited frame while
        # the viewport shows the canonical one, and `show_frame` applied it
        # regardless — see the module docstring. Going through the structure
        # panel rather than calling load_structure directly keeps the chooser
        # in sync with what is actually displayed.
        self.win.structure_panel.select(start_rec.pdb)
        if self.win.record is None or self.win.record.pdb != start_rec.pdb:
            self.win.load_structure(start_rec.pdb)
        st_a = self.win.structure
        if st_a is None or self.win.record is None \
                or self.win.record.pdb != start_rec.pdb:
            self.win._set_status(
                f"{start_rec.pdb} could not be loaded — morph not built")
            return

        try:
            st_b = Structure.from_file(end_rec.path)
            a_blocks, a_res = self._blocks_and_residues(st_a)
            b_blocks, b_res = self._blocks_and_residues(st_b)
            start, end, common, info = prepare_endpoints(a_blocks, a_res,
                                                         b_blocks, b_res)
            modes = None
            if params["method"] == "modal":
                per = len(common)
                anm = ANM.from_trimer([start[p * per:(p + 1) * per]
                                       for p in range(3)], cutoff=15.0).build()
                modes = anm.calc_modes(n_modes=30)
            traj = morph(start, end, n_frames=41, method=params["method"],
                         modes=modes)
        except Exception as exc:
            self.win._set_status(f"morph failed — {type(exc).__name__}: {exc}")
            return

        self.trajectory = traj
        self.residues = common
        self._start = st_a
        self._phase = 0.0
        self._fraction = 0.0
        self._build_map(st_a, common)
        self.win.physics_panel.set_morph(traj)
        # Frame on the widest point of the path, so nothing swings out of view
        # as the user drags the slider.
        self._frame_whole_path()
        # The tag, put back and made to follow. The deposited end entry's own
        # model is measured first — a distance is invariant under the rigid
        # placement, so it can be read off the file as deposited — so that the
        # far end of the path can be quoted beside the thing it stands in for.
        if wanted_tag:
            end_model = self.win.fusion.model_for(st_b)
            self._fusion_endpoint = (
                None if end_model is None
                else float(end_model.pore_exit_distances()[0]))
            self.win.fusion.show(True)     # rebuilds, then calls refresh_fusion
        # What the far end of the slider actually is. `restrained` and `linear`
        # land on the target's C-alpha positions exactly; `modal` is confined
        # to the elastic-network subspace and deliberately stops short, so the
        # two cannot be described by one sentence.
        captured = traj.meta.get("fraction_captured_by_modes")
        ends = (f"the path ends on {end_rec.pdb}'s C-alpha positions"
                if captured is None else
                f"the elastic-network subspace captures {captured:.0%} of the "
                f"change, so the path stops short of {end_rec.pdb}")
        self.win._set_status(
            f"morph {start_rec.pdb} → {end_rec.pdb}: {len(traj)} frames, "
            f"{info['n_common_residues']} common residues, endpoint RMSD "
            f"{info['endpoint_rmsd']:.1f} Å, worst bond error "
            f"{traj.bond_error.max():.2f} Å. Only C-alphas are interpolated — "
            f"{ends}, and every other atom is carried with its residue, so the "
            f"last frame is not the deposited {end_rec.pdb}"
            + (" (protomer order was reversed)"
               if info["handedness_flipped"] else "")
            + self._fusion_note())

    def _frame_whole_path(self) -> None:
        """Fit the camera to the union of the first, middle and last frames.

        Framing on the start structure alone would let the blades swing outside
        the viewport as the dome flattens, since the flat state is
        substantially wider.
        """
        if self.trajectory is None or self.win.viewport.scene is None:
            return
        tr = self.trajectory
        sample = np.vstack([tr.frames[0], tr.frames[len(tr) // 2], tr.frames[-1]])
        self.win.viewport.scene.camera.frame(sample)
        self.win.viewport.update()

    def _blocks_and_residues(self, st: Structure):
        chains = []
        for ch in st.chains:
            m = st.mask_ca() & (st.chain == ch)
            if m.sum() > 300:
                chains.append((st.xyz[m], st.res_seq[m]))
        common = set(chains[0][1].tolist())
        for _, seq in chains[1:3]:
            common &= set(seq.tolist())
        arr = np.array(sorted(common))
        return ([xyz[np.searchsorted(seq, arr)].astype(np.float64)
                 for xyz, seq in chains[:3]], arr)

    def _build_map(self, st: Structure, common: np.ndarray) -> None:
        """Map every atom of the displayed structure onto a morph site.

        ``st`` is the structure the path was built from, so the chain list here
        and the one in :meth:`_blocks_and_residues` are derived from the same
        object and cannot disagree about which block a protomer is.

        An atom at a residue *in* the shared basis takes its own residue's
        site, which keeps whole residues together. Everything else — bound
        lipids, an ion, a residue one endpoint resolves and the other does not
        — takes the **nearest site in space**. Rounding by residue number
        instead tied 7WLT's 1,407 lipid atoms to a C-alpha a median of 64.8 A
        away, purely because a lipid is numbered 2601 and the last shared
        residue is 2546, and they then travelled with the CTD tip: one of the
        largest displacements in the whole motion.
        """
        from scipy.spatial import cKDTree

        per = len(common)
        sites = self.trajectory.frames[0]
        chains = [ch for ch in st.chains
                  if (st.mask_ca() & (st.chain == ch)).sum() > 300][:3]
        mapping = np.full(st.n_atoms, -1, dtype=np.int64)
        for p, ch in enumerate(chains):
            idx = np.where(st.chain == ch)[0]
            res = st.res_seq[idx]
            pos = np.searchsorted(common, res)
            exact = (pos < per) & (common[np.clip(pos, 0, per - 1)] == res)
            mapping[idx[exact]] = p * per + pos[exact]
        loose = mapping < 0
        if loose.any():
            mapping[loose] = cKDTree(sites).query(st.xyz[loose])[1]
        self._map = mapping
        self._base = st.xyz.copy()

    def _coords_at(self, frame: np.ndarray) -> np.ndarray:
        """The drawn coordinates for one point on the path.

        The one expression both the picture and anything measuring it go
        through. The cast comes before the addition rather than after, and that
        is not a detail: `float32 + float32` and `float64 then round` differ by
        a few thousandths of an Angstrom, and a tag solved against one set of
        coordinates while drawn on the other is two structures again.
        """
        delta = (frame - self.trajectory.frames[0])[self._map]
        return self._base + delta.astype(np.float32)

    def _frame_structure(self, index: int) -> Structure:
        """One stored frame as a structure, for anything that measures it.

        The atoms are the start entry's; only their positions differ.
        """
        return self._start.copy_with_coords(
            self._coords_at(self.trajectory.frames[index]))

    def refresh_fusion(self) -> None:
        """Solve the HaloTag placement once per frame, or drop it.

        The tag is anchored to a C-alpha of the channel, so it has to travel
        with the morph or it hangs in space beside a flattened dome. Carrying it
        *rigidly* with its anchor is not enough: the placement is the centroid
        of the region the tag can occupy without clashing, and flattening takes
        a third of that region away (242 → 177 nm³ on 7WLT → 7WLU). A rigid
        carry misses that and lands about 7 Å out at the far end.

        So the model is re-solved on every frame. That costs ~140 ms a frame —
        the whole path is about six seconds — which is why it is done once, here,
        rather than in :meth:`show_frame`. It is skipped entirely unless the tag
        is actually being shown, and it holds each frame's envelope rather than
        recomputing it when the cloud is switched on: 20-30 MB for a path, which
        buys a drawing cost of nothing per frame.
        """
        if self.trajectory is None or not self.win.fusion.visible:
            self._fusion_frames = None
            return
        n = len(self.trajectory)
        frames, clashing = [], 0
        for i in range(n):
            if i % 8 == 0:
                self.win._set_status(
                    f"modelling the HaloTag along the path — frame {i + 1} "
                    f"of {n}…")
                QApplication.processEvents()
            model = self.win.fusion.model_for(self._frame_structure(i))
            frames.append(model)
            if model is not None and model.meta.get("clashes"):
                clashing += 1
        self._fusion_frames = frames
        self._fusion_clashing = clashing
        self.show_frame(self._fraction)

    def _fusion_note(self) -> str:
        """What the carried tag is and is not, measured rather than asserted."""
        frames = [f for f in (self._fusion_frames or []) if f is not None]
        if not frames:
            return ""
        first, last = frames[0], frames[-1]
        note = (f" HaloTag carried on its anchor at residue "
                f"{first.anchor_residues[0]}: accessible volume "
                f"{first.volume.volume:.0f} → {last.volume.volume:.0f} nm³ "
                f"along the path, re-solved at each of {len(frames)} frames")
        if self._fusion_clashing:
            note += (f"; the representative position is closer to the channel "
                     f"than the tag's own radius of gyration at "
                     f"{self._fusion_clashing} of them, which is a result about "
                     f"the flattened state and not a drawing error")
        # The far end carries the *start* entry's side chains and unshared
        # residues, and the pore exit is whichever atom reaches furthest down
        # the axis — so the distance there is not the deposited entry's own.
        if self._fusion_endpoint is not None:
            note += (f". Tag-to-pore-exit at the far end "
                     f"{last.pore_exit_distances()[0]:.2f} nm against "
                     f"{self._fusion_endpoint:.2f} nm on the deposited entry "
                     f"itself — the gap is the carried side chains, so quote "
                     f"the deposited value, not this one")
        # This line replaces whatever the fusion controller last wrote, and a
        # drawn fold may never be on screen without the statement that its spin
        # is undetermined. So the statement comes along rather than being lost.
        if self.win.fusion.show_atoms:
            note += (". THE SPIN IS UNDETERMINED — the fold is one draw; turn "
                     "it with View → HaloTag fusion → Turn tag orientation")
        return note

    def show_frame(self, fraction: float) -> None:
        if self.trajectory is None or self.win.view is None:
            return
        self._fraction = float(np.clip(fraction, 0.0, 1.0))
        # Whole residues are carried with their C-alpha site — see `_coords_at`.
        self.win.view.update_coords(
            self._coords_at(self.trajectory.at(fraction)))
        if self._fusion_frames is not None:
            # The nearest solved frame, not an interpolation between two of
            # them: a tag centre interpolated between two solved positions can
            # pass through the channel, and every stored one is a position the
            # envelope actually admitted. At 41 frames the step is 2.5% of the
            # path.
            i = int(round(self._fraction * (len(self._fusion_frames) - 1)))
            self.win.fusion.set_frame(self._fusion_frames[i],
                                      self.win.view.structure)
        self.win.viewport.update()

    def play(self, on: bool) -> None:
        self.win.viewport.clear_animations()
        self.win.physics_panel.morph_play.setText("Stop" if on else "Play")
        hud = self.win.viewport.hud
        if not on or self.trajectory is None:
            hud.set_clock("")
            return

        self._elapsed = 0.0
        self._frames = 0

        def step(dt: float) -> bool:
            self._phase += dt * 0.35
            self._elapsed += dt
            self._frames += 1
            # Ping-pong so the user sees both directions of the transition.
            f = abs((self._phase % 2.0) - 1.0)
            self.win.physics_panel.morph_slider.blockSignals(True)
            self.win.physics_panel.morph_slider.setValue(int(f * 100))
            self.win.physics_panel.morph_slider.blockSignals(False)
            self.show_frame(f)
            # Reported as a fraction along the path, never as a time: a morph
            # is an interpolation between two endpoints, not a trajectory, and
            # a seconds axis would imply kinetics the model does not contain.
            hud.set_clock(
                f"{f * 100:5.1f}% along the path   frame {self._frames}",
                f"interpolation, not a trajectory · {self._elapsed:.1f} s "
                f"of playback")
            return True

        self.win.viewport.add_animation(step)

