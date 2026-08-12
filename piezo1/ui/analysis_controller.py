"""Drives the Analysis dock: pore, pockets, conservation, allostery.

Every one of these takes seconds to tens of seconds on a 2500-residue trimer, so
they all run on a worker thread. A window that stops repainting during a
scientific calculation is indistinguishable from one that has crashed, and the
user's only recourse is to kill it.

The results are per-residue scalars, which are pushed into the renderer through
the existing ``ColorBy.VALUE`` path rather than a new one — the same mechanism
mode-displacement colouring already uses.
"""

from __future__ import annotations

import numpy as np
from PyQt6.QtCore import QObject, QThread, pyqtSignal

from ..render.representations import ColorBy

__all__ = ["AnalysisController", "AnalysisWorker"]


class AnalysisWorker(QObject):
    """Runs one named analysis off the GUI thread."""

    finished = pyqtSignal(str, object)
    failed = pyqtSignal(str, str)

    def __init__(self, kind: str, payload: dict) -> None:
        super().__init__()
        self.kind = kind
        self.payload = payload

    def run(self) -> None:
        try:
            self.finished.emit(self.kind, getattr(self, f"_{self.kind}")())
        except Exception as exc:                      # shown in the status bar
            self.failed.emit(self.kind, f"{type(exc).__name__}: {exc}")

    # ------------------------------------------------------------- analyses

    def _pore(self) -> dict:
        from ..analysis.hydration import (hydrophobicity_profile_chap,
                                          load_grid, predict_wetting)
        from ..structure.pore import pore_profile
        from ..structure.superpose import detect_c3_axis

        st = self.payload["structure"]
        axis = detect_c3_axis(self.payload["blocks"])
        profile = pore_profile(st, axis, step=self.payload.get("step", 1.0))

        grid = load_grid()
        hydro = hydration = None
        if grid.available:
            hydro = hydrophobicity_profile_chap(st, profile)
            hydration = predict_wetting(st, profile, grid)
        return {"profile": profile, "hydrophobicity": hydro,
                "hydration": hydration}

    def _pockets(self) -> dict:
        from ..analysis.pockets import find_pockets
        pockets = find_pockets(self.payload["structure"])
        return {"pockets": pockets[:self.payload.get("top", 10)]}

    def _conservation(self) -> dict:
        from ..analysis.conservation import conservation_profile, fetch_orthologs
        orthologs = fetch_orthologs()
        profile = conservation_profile(orthologs)
        # Positions where too few orthologs align carry a conservation value
        # that is really a statement about the alignment, not about selection
        # pressure. Dropping them is why the map has gaps rather than a
        # confident flat band across the unresolved blade.
        keep = profile.coverage >= 0.7
        return {"profile": profile,
                "values": {int(r): float(v) for r, v, ok
                           in zip(profile.residues, profile.conservation, keep)
                           if ok and np.isfinite(v)}}

    def _path(self) -> dict:
        """The blade-to-gate route, and how far it is from being unique.

        The degeneracy number is the point of doing this here rather than in
        the drawing code. A single line on screen reads as *the* pathway, so
        the same search is re-run with this route's own steps deleted from the
        graph and the best remaining route is costed. Both numbers reach the
        status line together.

        The degeneracy measurement itself lives in
        :func:`piezo1.ui.path_controller.alternative_cost`, so it can be
        calibrated on a graph whose answer is known — a route with one bridge
        must come back as unique — rather than only on the trimer, where
        nothing independently says what the right answer is.
        """
        from ..analysis.allostery import allosteric_path, cross_correlation
        from .path_controller import alternative_cost, path_endpoints

        coords = np.vstack(self.payload["blocks"])
        residues = np.tile(np.asarray(self.payload["residues"]), 3)
        dcc = cross_correlation(self.payload["modes"])
        source, target, source_name = path_endpoints(
            residues, self.payload["annotations"])
        if not source or not target:
            raise ValueError(
                "no blade or gate residues are resolved in this structure")
        path = allosteric_path(coords, dcc, source, target, residues)
        alternative = alternative_cost(coords, dcc, source, target, residues,
                                       path.sites)

        return {"path": path, "residues": list(path.residues),
                "sites": list(path.sites), "cost": float(path.cost),
                "correlations": list(path.correlations),
                "coords": coords[np.asarray(path.sites, dtype=int)],
                "alternative_cost": alternative, "source_name": source_name}

    def _allostery(self) -> dict:
        from ..analysis.allostery import perturbation_response
        modes = self.payload["modes"]
        residues = np.asarray(self.payload["residues"])
        prs = perturbation_response(modes, np.tile(residues, 3),
                                    normalise=False)
        coupling = np.asarray(prs.sensitivity, dtype=float)
        n = len(residues)
        # Average the three protomers: the trimer is C3 symmetric, so a
        # per-protomer difference is numerical noise rather than biology, and
        # colouring one protomer differently from its mates looks like a bug.
        folded = coupling[:n * 3].reshape(3, n).mean(axis=0)
        return {"prs": prs,
                "values": {int(r): float(v) for r, v in zip(residues, folded)}}


class AnalysisController:
    """Owns the worker thread and pushes results into the panel and renderer."""

    #: Human-readable names for the scalars that can colour the model.
    SCALAR_LABELS = {"conservation": "Conservation (Shannon, vertebrates)",
                     "allostery": "Mechanical coupling to the gate (PRS)"}

    def __init__(self, window) -> None:
        self.win = window
        self.scalars: dict[str, dict[int, float]] = {}
        self.pore = None
        self.hydration = None
        self.pockets: list = []
        self._thread: QThread | None = None
        self._worker: AnalysisWorker | None = None
        self._active = ""

    def reset(self) -> None:
        self.scalars.clear()
        self.pore = self.hydration = None
        self.pockets = []
        panel = self.win.analysis_panel
        panel.set_pore(None)
        panel.set_pockets([])
        panel.set_top_residues([])

    # ------------------------------------------------------------- launching

    def _start(self, kind: str, payload: dict) -> None:
        if self._thread is not None:
            self.win._set_status(f"{self._active} already running")
            return
        if self.win.structure is None:
            return
        self._active = kind
        self.win.analysis_panel.set_busy(True, kind)
        self.win._set_status(f"computing {kind}…")

        self._thread = QThread()
        self._worker = AnalysisWorker(kind, payload)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_finished)
        self._worker.failed.connect(self._on_failed)
        self._thread.start()

    def compute_pore(self) -> None:
        if not self.win._mode_blocks:
            self.win._set_status("pore profile needs three resolved protomers")
            return
        self._start("pore", {"structure": self.win.structure,
                             "blocks": self.win._mode_blocks})

    def compute_pockets(self, top: int = 10) -> None:
        self._start("pockets", {"structure": self.win.structure, "top": top})

    def compute_conservation(self) -> None:
        self._start("conservation", {"structure": self.win.structure})

    def compute_path(self) -> None:
        if self.win.modes is None or not self.win._mode_blocks:
            self.win._set_status(
                "the allosteric path needs normal modes — compute them in "
                "Physics first")
            return
        self._start("path", {"modes": self.win.modes,
                             "blocks": self.win._mode_blocks,
                             "residues": self.win._mode_residues,
                             "annotations": self.win.path.annotations()})

    def compute_allostery(self) -> None:
        if self.win.modes is None:
            self.win._set_status(
                "allostery needs normal modes — compute them in Physics first")
            return
        self._start("allostery", {"modes": self.win.modes,
                                  "residues": self.win._mode_residues})

    def cleanup(self) -> None:
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait()
            self._thread = None
            self._worker = None
        self.win.analysis_panel.set_busy(False)

    # -------------------------------------------------------------- results

    def _on_failed(self, kind: str, message: str) -> None:
        self.cleanup()
        self.win.analysis_panel.set_message(kind, f"failed — {message}")
        self.win._set_status(f"{kind} failed — {message}")

    def _on_finished(self, kind: str, result: dict) -> None:
        self.cleanup()
        panel = self.win.analysis_panel
        if kind == "pore":
            self.pore = result["profile"]
            self.hydration = result["hydration"]
            panel.set_pore(self.pore, self.hydration, result["hydrophobicity"])
            note = f"bottleneck {self.pore.bottleneck_radius:.2f} Å"
            if self.hydration is not None and self.hydration.available:
                note += f" · {self.hydration.verdict}"
            self.win._set_status(note)
            # The drawn pore reads this same profile object, so the picture
            # and the plot can never be of different runs.
            self.win.pore_surface.refresh()
            # The nanodomain's source point is the cytosolic end of this same
            # profile, so it waits on the same run.
            self.win.nanodomain.refresh()

        elif kind == "pockets":
            self.pockets = result["pockets"]
            panel.set_pockets(self.pockets)
            self.win._set_status(f"{len(self.pockets)} pockets")
            self.win.pocket_view.refresh()

        elif kind == "path":
            self.win.path.refresh(result)

        elif kind in ("conservation", "allostery"):
            self.scalars[kind] = result["values"]
            panel.add_scalar(kind, self.SCALAR_LABELS[kind])
            values = result["values"]
            top = sorted(values.items(), key=lambda kv: -kv[1])[:15]
            panel.set_top_residues(top, self._describe(kind, values))
            panel.set_message(kind, self._describe(kind, values))
            self.win._set_status(f"{kind}: {len(values)} residues scored")

    def _describe(self, kind: str, values: dict[int, float]) -> str:
        array = np.array(list(values.values()))
        if kind == "conservation":
            return (f"{len(values)} residues · mean {array.mean():.3f} · "
                    f"{int((array > 0.95).sum())} above 0.95")
        return (f"{len(values)} residues · coupling "
                f"{array.min():.3g}–{array.max():.3g}")

    # -------------------------------------------------------------- colouring

    def apply_color(self, key: str, on: bool) -> None:
        """Colour the model by a per-residue scalar via ``ColorBy.VALUE``."""
        view = self.win.view
        if view is None or self.win.structure is None:
            return
        if not on or not key or key not in self.scalars:
            view.color_by = self.win._current_color()
            view.values = None
            self.coloured_key = ""
        else:
            view.values = self.residue_values_to_atoms(self.scalars[key])
            view.color_by = ColorBy.VALUE
            # Remembered so the export can write the RAW residue map rather
            # than `view.values`, which has unmeasured residues filled to the
            # map floor for display. Exporting the filled array would make an
            # unscored residue indistinguishable from a low-scoring one in
            # somebody else's viewer.
            self.coloured_key = key
        view.rebuild()
        self.win.viewport.update()

    #: Which scalar is currently painted on the model, or "".
    coloured_key = ""

    def residue_values_to_atoms(self, values: dict[int, float]) -> np.ndarray:
        """Spread a per-residue scalar over atoms.

        Residues with no value get the minimum rather than zero: the colour
        scale is normalised over whatever it is handed, so injecting a zero
        into a conservation map that runs 0.6–1.0 would rescale the entire
        legend around positions that were never measured.
        """
        st = self.win.structure
        assert st is not None
        floor = min(values.values()) if values else 0.0
        lookup = np.full(int(st.res_seq.max()) + 2, floor, dtype=np.float32)
        for residue, value in values.items():
            if 0 <= residue < len(lookup):
                lookup[residue] = value
        return lookup[np.clip(st.res_seq, 0, len(lookup) - 1)]

    def focus_pore_position(self, z: float) -> None:
        """Select the residues lining the pore at a clicked height."""
        if self.pore is None:
            return
        index = int(np.argmin(np.abs(self.pore.z - z)))
        sl = self.pore.slices[index]
        if sl.lining:
            self.win._highlight(list(sl.lining),
                                f"pore lining at z = {sl.z:.1f} Å "
                                f"(radius {sl.radius:.2f} Å)")
