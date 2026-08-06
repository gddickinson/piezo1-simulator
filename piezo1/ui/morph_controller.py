"""Drives the curved-to-flattened conformational morph.

Split out of :mod:`piezo1.ui.main_window` to keep both files readable. The
controller owns the trajectory, the atom-to-site mapping and the playback
animation; the window owns everything else.
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
        self._phase = 0.0

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
        try:
            st_a = Structure.from_file(start_rec.path)
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

        # Show the morph on the start structure. Going through the structure
        # panel rather than calling load_structure directly keeps the chooser
        # in sync with what is actually on screen.
        self.win.structure_panel.select(start_rec.pdb)
        if self.win.record is None or self.win.record.pdb != start_rec.pdb:
            self.win.load_structure(start_rec.pdb)
        self.trajectory = traj
        self.residues = common
        self._phase = 0.0
        self._build_map(common)
        self.win.physics_panel.set_morph(traj)
        # Frame on the widest point of the path, so nothing swings out of view
        # as the user drags the slider.
        self._frame_whole_path()
        self.win._set_status(
            f"morph {start_rec.pdb} → {end_rec.pdb}: {len(traj)} frames, "
            f"{info['n_common_residues']} common residues, endpoint RMSD "
            f"{info['endpoint_rmsd']:.1f} Å, worst bond error "
            f"{traj.bond_error.max():.2f} Å"
            + (" (protomer order was reversed)"
               if info["handedness_flipped"] else ""))

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

    def _build_map(self, common: np.ndarray) -> None:
        """Map every atom of the displayed structure onto a morph site."""
        st = self.win.structure
        assert st is not None
        per = len(common)
        chains = [ch for ch in st.chains
                  if (st.mask_ca() & (st.chain == ch)).sum() > 300][:3]
        mapping = np.zeros(st.n_atoms, dtype=np.int64)
        for p, ch in enumerate(chains):
            sel = st.chain == ch
            pos = np.clip(np.searchsorted(common, st.res_seq[sel]), 0, per - 1)
            mapping[sel] = p * per + pos
        self._map = np.clip(mapping, 0, 3 * per - 1)
        self._base = self.win.structure.xyz.copy()

    def show_frame(self, fraction: float) -> None:
        if self.trajectory is None or self.win.view is None:
            return
        frame = self.trajectory.at(fraction)
        start = self.trajectory.frames[0]
        # Carry whole residues with their C-alpha site.
        delta = (frame - start)[self._map]
        self.win.view.update_coords(self._base + delta.astype(np.float32))
        self.win.viewport.update()

    def play(self, on: bool) -> None:
        self.win.viewport.clear_animations()
        self.win.physics_panel.morph_play.setText("Stop" if on else "Play")
        if not on or self.trajectory is None:
            return

        def step(dt: float) -> bool:
            self._phase += dt * 0.35
            # Ping-pong so the user sees both directions of the transition.
            f = abs((self._phase % 2.0) - 1.0)
            self.win.physics_panel.morph_slider.blockSignals(True)
            self.win.physics_panel.morph_slider.setValue(int(f * 100))
            self.win.physics_panel.morph_slider.blockSignals(False)
            self.show_frame(f)
            return True

        self.win.viewport.add_animation(step)

