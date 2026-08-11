"""Dome measurement and elastic-network mode analysis for the GUI.

Split out of :mod:`piezo1.ui.main_window`. Owns the mode set, the worker
thread that computes it, and the animation that displays it.
"""

from __future__ import annotations

import time

import numpy as np
from PyQt6.QtCore import QObject, QThread, pyqtSignal

from ..config import RESOURCE_DIR
from ..physics.anm import ANM
from ..render.representations import ColorBy
from ..structure.geometry import measure_dome

__all__ = ["PhysicsController", "ModeWorker"]


class ModeWorker(QObject):
    """Runs the elastic-network calculation off the GUI thread."""

    finished = pyqtSignal(object, object, float)   # modes, anm, seconds
    failed = pyqtSignal(str)

    def __init__(self, blocks: list[np.ndarray], params: dict) -> None:
        super().__init__()
        self.blocks = blocks
        self.params = params

    def run(self) -> None:
        try:
            t0 = time.time()
            anm = ANM.from_trimer(self.blocks,
                                  cutoff=self.params["cutoff"],
                                  spring=self.params["spring"]).build()
            modes = anm.calc_modes(n_modes=self.params["n_modes"])
            anm.label_symmetry(modes)
            self.finished.emit(modes, anm, time.time() - t0)
        except Exception as exc:  # surfaced in the status bar
            self.failed.emit(f"{type(exc).__name__}: {exc}")



class PhysicsController:
    """Dome measurement and normal-mode analysis, driven by the physics panel."""

    def __init__(self, window) -> None:
        self.win = window
        self.anm = None
        self._mode_index = 0
        self._amplitude = 18.0
        self._phase = 0.0
        self._base_coords = None
        self._thread = None
        self._ca_map = None
        self._ca_map_for = None

    def reset(self) -> None:
        self.anm = None
        self._base_coords = None
        self._ca_map = None
        self._ca_map_for = None

    def measure_dome(self) -> None:
        if not self.win._mode_blocks or self.win.structure is None:
            self.win.physics_panel.set_dome(None)
            return
        species = self.win.record.numbering_species if self.win.record else "human"
        surface = self._tm_surface(species)
        if surface is None or len(surface) < 12:
            self.win.physics_panel.set_dome(None)
            return
        try:
            dome = measure_dome(self.win._mode_blocks, surface)
        except Exception as exc:
            self.win._set_status(f"dome measurement failed: {exc}")
            return
        ref = ("<br>published closed-state value: 10.2 nm "
               "(Haselwandter &amp; MacKinnon 2018)")
        self.win.physics_panel.set_dome(dome, ref)
        self.win._set_status(f"dome: {dome.summary()}")

    def _tm_surface(self, species: str) -> np.ndarray | None:
        """Mid-membrane surface points: the centre of every TM helix."""
        import json
        path = RESOURCE_DIR / f"uniprot_{species}.json"
        if not path.exists() or self.win.structure is None:
            return None
        tms = json.loads(path.read_text())["transmembrane"]
        st = self.win.structure
        pts = []
        for ch in st.chains:
            m = st.mask_ca() & (st.chain == ch)
            if m.sum() < 300:
                continue
            xyz, seq = st.xyz[m], st.res_seq[m]
            for tm in tms:
                mid = 0.5 * (tm["start"] + tm["end"])
                half = max(2.0, (tm["end"] - tm["start"]) / 6.0)
                sel = (seq >= mid - half) & (seq <= mid + half)
                if sel.sum() >= 3:
                    pts.append(xyz[sel].mean(axis=0))
        return np.array(pts) if pts else None

    def compute_modes(self, params: dict) -> None:
        if not self.win._mode_blocks:
            self.win._set_status("need three equivalent protomers for an ANM")
            return
        if self._thread is not None:
            return
        n = len(self.win._mode_blocks[0]) * 3
        self.win._set_status(f"building elastic network over {n:,} C-alpha sites…")
        self.win.physics_panel.set_busy(True)

        self._thread = QThread()
        self._worker = ModeWorker(self.win._mode_blocks, params)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_modes)
        self._worker.failed.connect(self._on_modes_failed)
        self._thread.start()

    def cleanup(self) -> None:
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait()
            self._thread = None
        self.win.physics_panel.set_busy(False)

    def _on_modes(self, modes, anm, seconds: float) -> None:
        self.cleanup()
        self.win.modes = modes
        self.anm = anm
        self._base_coords = np.vstack(self.win._mode_blocks)
        self.win.physics_panel.set_modes(modes)
        n_a = int((modes.symmetry == "A").sum()) if modes.symmetry is not None else 0
        self.win._set_status(
            f"{modes.n_modes} modes in {seconds:.1f} s · {n_a} symmetric (A), "
            f"{modes.n_modes - n_a} degenerate (E) · "
            f"lowest eigenvalue {modes.eigenvalues[0]:.5f}")

    def _on_modes_failed(self, message: str) -> None:
        self.cleanup()
        self.win._set_status(f"mode calculation failed — {message}")

    def select_mode(self, index: int) -> None:
        self._mode_index = index
        self._phase = 0.0

    def set_amplitude(self, value: float) -> None:
        self._amplitude = value

    #: Angular frequency of the displayed mode sweep, rad/s. A display choice:
    #: the true frequency of an elastic-network mode is not physical time, and
    #: the clock says so rather than implying a real period.
    SWEEP_RATE = 1.6

    def animate_mode(self, on: bool) -> None:
        self.win.viewport.clear_animations()
        hud = self.win.viewport.hud
        if not on or self.win.modes is None or self.win.structure is None:
            if self.win.view is not None and self._base_coords is not None:
                self._apply_mode_displacement(0.0)
            hud.set_clock("")
            return

        self._elapsed = 0.0
        self._frames = 0
        period = 2.0 * np.pi / self.SWEEP_RATE

        def step(dt: float) -> bool:
            self._phase += dt * self.SWEEP_RATE
            self._elapsed += dt
            self._frames += 1
            self._apply_mode_displacement(np.sin(self._phase))
            cycle = (self._elapsed % period) / period
            hud.set_clock(
                f"{self._elapsed:6.2f} s   frame {self._frames}",
                f"mode {self._mode_index + 1} sweep, "
                f"{cycle * 100:.0f}% of a {period:.1f} s display cycle "
                f"(not a physical period)")
            return True

        self.win.viewport.add_animation(step)

    def _apply_mode_displacement(self, scale: float) -> None:
        """Displace the C-alpha network and carry whole residues with it."""
        if self.win.modes is None or self.win.structure is None or self.win.view is None:
            return
        disp = self.win.modes.mode(self._mode_index, self._amplitude) * scale
        st = self.win.structure
        # Map the per-protomer C-alpha displacement back onto every atom by
        # nearest C-alpha within the same chain.
        if self._ca_map is None or self._ca_map_for != st.name:
            self._build_ca_map()
        moved = st.xyz.copy()
        flat = disp.reshape(-1, 3)
        moved += flat[self._ca_map]
        self.win.view.update_coords(moved)
        self.win.viewport.update()

    def _build_ca_map(self) -> None:
        """For each atom, the index of its protomer-block C-alpha site."""
        st = self.win.structure
        assert st is not None
        blocks_per = len(self.win._mode_blocks[0])
        chains = [ch for ch in st.chains
                  if (st.mask_ca() & (st.chain == ch)).sum() > 300][:3]
        mapping = np.zeros(st.n_atoms, dtype=np.int64)
        for p, ch in enumerate(chains):
            m = st.mask_ca() & (st.chain == ch)
            ca_seq = st.res_seq[m]
            common = np.array(sorted(set.intersection(
                *[set(st.res_seq[st.mask_ca() & (st.chain == c)].tolist())
                  for c in chains])))
            atom_mask = st.chain == ch
            # nearest common residue for every atom of this chain
            pos = np.searchsorted(common, st.res_seq[atom_mask])
            pos = np.clip(pos, 0, len(common) - 1)
            mapping[atom_mask] = p * blocks_per + pos
            del ca_seq
        self._ca_map = np.clip(mapping, 0, 3 * blocks_per - 1)
        self._ca_map_for = st.name

    def color_by_mode(self, on: bool) -> None:
        if self.win.view is None or self.win.modes is None or self.win.structure is None:
            return
        if on:
            self._untick("fluctuation_button")
        if not on:
            self.win.view.color_by = self.win._current_color()
            self.win.view.values = None
        else:
            if self._ca_map is None:
                self._build_ca_map()
            mag = np.linalg.norm(self.win.modes.vectors[self._mode_index], axis=1)
            self.win.view.values = mag[self._ca_map]
            self.win.view.color_by = ColorBy.VALUE
        self.win.view.rebuild()
        self.win.viewport.update()

    def color_by_fluctuation(self, on: bool) -> None:
        """Colour by the network's predicted mean-square fluctuation.

        The quantity is `ModeSet.msf` — Σ|v|²/λ over every computed mode — read
        rather than recomputed, so this is the same array
        :func:`piezo1.analysis.fluctuations.predicted_msf` correlates against
        the deposited B-factors. If the two ever disagreed, the picture would
        be showing something the validation never tested.

        Deliberately **not** averaged over the protomers, unlike
        `predicted_msf`. That average exists because the *observation* is one
        B-factor per residue per chain and counting one prediction three times
        would be wrong; here the three copies are on screen separately, and
        folding them would paint a difference that the model does not have as
        if it were agreement.
        """
        if self.win.view is None or self.win.modes is None or self.win.structure is None:
            return
        if on:
            self._untick("color_button")
        if not on:
            self.win.view.color_by = self.win._current_color()
            self.win.view.values = None
        else:
            if self._ca_map is None:
                self._build_ca_map()
            msf = np.asarray(self.win.modes.msf(), dtype=float)
            self.win.view.values = msf[self._ca_map]
            self.win.view.color_by = ColorBy.VALUE
            self.win._set_status(self.fluctuation_line())
        self.win.view.rebuild()
        self.win.viewport.update()

    def _untick(self, button: str) -> None:
        """Turn the other value-colouring off without re-entering its handler.

        Both buttons drive `ColorBy.VALUE` through the same `view.values`
        slot, so leaving one checked while the other paints would show a lit
        button describing a colour that is not on screen.
        """
        widget = getattr(self.win.physics_panel, button, None)
        if widget is not None and widget.isChecked():
            widget.blockSignals(True)
            widget.setChecked(False)
            widget.blockSignals(False)

    def fluctuation_line(self) -> str:
        """What the colours mean, and what says whether to believe them.

        The scale is arbitrary — a mean-square fluctuation in a network with an
        unfitted spring constant has no units anyone can compare with an
        Angstrom. What makes it a prediction rather than a picture is the
        measured agreement with the deposited B-factors, and that lives one
        menu away, so this points at it rather than quoting a number the
        colouring did not compute.
        """
        modes = self.win.modes
        n = 0 if modes is None else modes.n_modes
        return (f"colouring by predicted mean-square fluctuation over {n} "
                f"modes (sum |v|^2/lambda) · the scale is ARBITRARY: this "
                f"network has no fitted spring constant, so only the ordering "
                f"means anything · whether the ordering is right is measured "
                f"in Analysis -> Fluctuation vs B-factor, where the network's "
                f"median Spearman is 0.74 against a burial-only control's "
                f"0.32 — but on Pearson 0.48 against 0.39")

