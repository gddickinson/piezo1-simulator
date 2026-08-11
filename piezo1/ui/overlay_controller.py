"""Superpose and display a second structure over the loaded one.

The comparison this project keeps needing — curved against flattened, wild type
against a variant, human against mouse — is two structures in one frame. Doing
it by eye across two windows is useless; the difference between the closed and
open PIEZO1 dome is a few Angstrom over a 280 Å propeller.

**The trap this handles.** Deposited entries do not label protomers in a
consistent rotational order — Round 4 found four entries labelled the reverse
way round. Overlaying 7WLT with 7WLU here, correspondence search rematches the
protomers to (2, 1, 0); taken at chain-label face value the two would sit 90.7 Å
apart against 12.3 Å once matched. A viewer that silently trusted labels would
present a wrong overlay that looks like an enormous conformational change,
which is exactly the conclusion someone would want to draw from it. So the
correspondence is searched and both numbers are reported.
"""

from __future__ import annotations

import dataclasses
from dataclasses import dataclass, field

import numpy as np
from PyQt6.QtCore import QObject, QThread, pyqtSignal

from ..core.structure import Structure
from ..render.representations import ColorBy, MolecularView, Style
from .model_utils import protomer_blocks

__all__ = ["OverlayController", "OverlayWorker", "OverlayResult"]

OVERLAY_NAME = "overlay"


@dataclass
class OverlayResult:
    """A superposed second structure and how well it fitted."""

    structure: Structure
    rotation: np.ndarray
    translation: np.ndarray
    centroid: np.ndarray
    rmsd: float
    rmsd_by_label: float
    n_common: int
    protomer_order: tuple[int, ...] = (0, 1, 2)
    per_residue: dict[int, float] = field(default_factory=dict)
    meta: dict = field(default_factory=dict)

    @property
    def reordered(self) -> bool:
        return tuple(self.protomer_order) != (0, 1, 2)

    def summary(self) -> str:
        text = (f"{self.structure.name}: RMSD {self.rmsd:.2f} Å over "
                f"{self.n_common} common C-alphas")
        if self.reordered:
            text += (f" · protomers rematched {self.protomer_order} "
                     f"(by chain label it would be {self.rmsd_by_label:.1f} Å)")
        return text


class OverlayWorker(QObject):
    """Superposition off the GUI thread — a trimer fit is not instant."""

    finished = pyqtSignal(object)
    failed = pyqtSignal(str)

    def __init__(self, reference: Structure, mobile: Structure,
                 mode: str = "protomer") -> None:
        super().__init__()
        self.reference = reference
        self.mobile = mobile
        self.mode = mode

    def run(self) -> None:
        try:
            self.finished.emit(self._superpose())
        except Exception as exc:
            self.failed.emit(f"{type(exc).__name__}: {exc}")

    # ------------------------------------------------------------ the fit

    def _common_ca(self, a: Structure, b: Structure):
        """C-alpha coordinates for residues resolved in both, chain A of each.

        Residue numbers are the join key, which is only valid within one
        numbering system. Pairs that are not — a different species, the PIEZO2
        paralogue, or a splice isoform numbered in its own coordinates — are
        refused by :meth:`OverlayController._numbering_refusal` before reaching
        here, rather than producing a confident wrong fit.
        """
        def table(st):
            mask = st.mask_ca() & (st.chain == st.chains[0])
            return dict(zip(st.res_seq[mask].tolist(), st.xyz[mask]))

        ta, tb = table(a), table(b)
        shared = sorted(set(ta) & set(tb))
        if len(shared) < 20:
            raise ValueError(f"only {len(shared)} residues in common")
        return (np.array([ta[r] for r in shared], dtype=float),
                np.array([tb[r] for r in shared], dtype=float), shared)

    def _shared_blocks(self, st, residues):
        """C-alpha blocks for one structure on a *given* residue basis."""
        blocks = []
        for chain in st.chains:
            mask = st.mask_ca() & (st.chain == chain)
            if mask.sum() < 300:
                continue
            seq, xyz = st.res_seq[mask], st.xyz[mask]
            index = {int(r): i for i, r in enumerate(seq)}
            if not all(r in index for r in residues):
                continue
            blocks.append(np.array([xyz[index[r]] for r in residues],
                                   dtype=float))
        return blocks[:3]

    def _superpose(self) -> OverlayResult:
        from ..structure.superpose import kabsch, match_protomers, rmsd

        order = (0, 1, 2)
        by_label = float("nan")
        if self.mode == "protomer":
            # The basis must be shared across BOTH structures. Building each
            # one's blocks independently gives different site counts whenever
            # the two entries resolve different residues, which is almost
            # always — and match_protomers cannot compare them.
            _ref, ref_res = protomer_blocks(self.reference)
            _mob, mob_res = protomer_blocks(self.mobile)
            common = sorted(set(ref_res.tolist()) & set(mob_res.tolist()))
            ref_blocks = self._shared_blocks(self.reference, common)
            mob_blocks = self._shared_blocks(self.mobile, common)
            if len(ref_blocks) == 3 and len(mob_blocks) == 3:
                match = match_protomers(ref_blocks, mob_blocks)
                order = tuple(match.order)
                by_label = float(rmsd(np.vstack(ref_blocks),
                                      np.vstack(mob_blocks)))

        fixed, moving, shared = self._common_ca(self.reference, self.mobile)
        # kabsch returns (rotation, translation, mobile_centroid); the fitted
        # coordinates are (x - centroid) @ R.T + t, so the centroid has to be
        # carried through to the full-structure transform below.
        rotation, translation, centroid = kabsch(moving, fixed)
        fitted = (moving - centroid) @ rotation.T + translation
        value = float(np.sqrt(((fitted - fixed) ** 2).sum(axis=1).mean()))
        deviation = {int(r): float(np.linalg.norm(f - x))
                     for r, f, x in zip(shared, fitted, fixed)}

        return OverlayResult(
            structure=self.mobile, rotation=rotation, translation=translation,
            centroid=centroid,
            rmsd=value,
            rmsd_by_label=by_label if np.isfinite(by_label) else value,
            n_common=len(shared), protomer_order=order,
            per_residue=deviation,
            meta={"mode": self.mode,
                  "reference": self.reference.name,
                  "max_deviation": max(deviation.values()) if deviation else 0.0})


class OverlayController:
    """Owns the second view, its styling and the superposition thread."""

    #: Distinguishable defaults. The reference keeps whatever the user chose;
    #: the overlay is deliberately a single flat colour so the two are told
    #: apart at a glance rather than by squinting at two domain palettes.
    DEFAULT_COLOUR = (0.98, 0.62, 0.30)

    def __init__(self, window) -> None:
        self.win = window
        self.view: MolecularView | None = None
        self.result: OverlayResult | None = None
        self._thread: QThread | None = None
        self._worker: OverlayWorker | None = None

    # -------------------------------------------------------------- lifecycle

    def load(self, pdb: str, mode: str = "protomer") -> None:
        if self.win.structure is None:
            self.win._set_status("load a reference structure first")
            return
        if self._thread is not None:
            return
        from ..config import STRUCTURE_DIR
        path = STRUCTURE_DIR / f"{pdb.upper()}.cif"
        if not path.exists():
            self.win._set_status(f"{pdb} is not downloaded")
            return

        mobile = Structure.from_file(path)
        refusal = self._numbering_refusal(pdb, mobile)
        if refusal:
            self.win._set_status(refusal)
            return

        self.win._set_status(f"superposing {pdb}…")
        self._thread = QThread()
        self._worker = OverlayWorker(self.win.structure, mobile, mode)
        self._worker.moveToThread(self._thread)
        self._thread.started.connect(self._worker.run)
        self._worker.finished.connect(self._on_ready)
        self._worker.failed.connect(self._on_failed)
        self._thread.start()

    def _numbering_refusal(self, pdb: str, mobile: Structure) -> str:
        """Why these two must not be joined on residue number, or "".

        Residue numbers are the join key, so the two entries have to be in the
        *same* numbering. That was checked against the registry's species field
        until it was noticed that 6KG7 — PIEZO2 — is filed as "mouse" exactly
        like every mouse Piezo1 entry, so the guard passed it: overlaying 7WLT
        on 6KG7 returned a confident 47.9 A over 920 "matched" C-alphas of
        which 6% were even the same amino acid, where the alignment-based
        comparison gives 4.36 A over 3,708. A number, in the right units,
        wrong by more than tenfold.

        So the check is now the measurement rather than the label, and it
        covers three cases with one question: a different species, a different
        protein, and a file numbered in a splice isoform's own coordinates.
        """
        from ..analysis.numbering_check import identify_numbering

        reference = identify_numbering(self.win.structure)
        moving = identify_numbering(mobile)
        name = self.win.structure.name or "the reference"

        for entry, identity in ((pdb, moving), (name, reference)):
            if not identity.confident:
                detail = (identity.splice.summary() if identity.splice
                          else identity.summary())
                return (f"{entry} is not in canonical numbering — {detail} — "
                        f"so residue numbers cannot be the join key and this "
                        f"overlay is refused")

        if moving.reference == reference.reference:
            return ""

        paralogue = moving.is_piezo2 != reference.is_piezo2
        why = ("different proteins" if paralogue else "different species")
        hint = (" — use Analysis → PIEZO2 comparison, which matches the two "
                "through a real alignment" if paralogue else "")
        return (f"{pdb} is {moving.reference} and {name} is "
                f"{reference.reference}: {why}, so residue numbers do not "
                f"correspond and this overlay is refused{hint}")

    def cleanup(self) -> None:
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait()
            self._thread = None
            self._worker = None

    def clear(self) -> None:
        """Remove the overlay and any per-residue deviation colouring."""
        self.cleanup()
        scene = self.win.viewport.scene
        if scene is not None:
            for key in list(scene.batches):
                if key.startswith(f"{OVERLAY_NAME}:"):
                    scene.remove(key)
        self.view = None
        self.result = None
        if self.win.view is not None:
            self.win.view.values = None
            self.win.view.color_by = self.win._current_color()
            self.win.view.rebuild()
        self.win.viewport.update()

    # ---------------------------------------------------------------- results

    def _on_failed(self, message: str) -> None:
        self.cleanup()
        self.win._set_status(f"overlay failed — {message}")
        if hasattr(self.win, "overlay_panel"):
            self.win.overlay_panel.set_result(None, message)

    def _on_ready(self, result: OverlayResult) -> None:
        self.cleanup()
        self.result = result

        # Structure is a dataclass held immutable by convention, so the moved
        # copy is a replace() rather than an in-place edit: the caller's model
        # must keep its own coordinates.
        moved = dataclasses.replace(
            result.structure,
            xyz=((result.structure.xyz - result.centroid) @ result.rotation.T
                 + result.translation).astype(np.float32))
        scene = self.win.viewport.scene
        if scene is None:
            return
        for key in list(scene.batches):
            if key.startswith(f"{OVERLAY_NAME}:"):
                scene.remove(key)

        self.view = MolecularView(scene, moved, name=OVERLAY_NAME)
        self.view.set_species(
            self.win.record.numbering_species if self.win.record else "human")
        self.view.style = Style.BACKBONE
        self.view.color_by = ColorBy.UNIFORM
        self.view.ligands_as_spheres = False
        self.view.rebuild()
        self.win.viewport.update()

        self.win._set_status(result.summary())
        if hasattr(self.win, "overlay_panel"):
            self.win.overlay_panel.set_result(result)

    # --------------------------------------------------------------- styling

    def set_style(self, style: Style) -> None:
        if self.view is not None:
            self.view.style = style
            self.view.rebuild()
            self.win.viewport.update()

    def set_visible(self, on: bool) -> None:
        scene = self.win.viewport.scene
        if scene is None:
            return
        for key, batch in scene.batches.items():
            if key.startswith(f"{OVERLAY_NAME}:"):
                batch.visible = bool(on)
        self.win.viewport.update()

    def color_reference_by_deviation(self, on: bool) -> None:
        """Colour the *reference* by how far the overlay moved at each residue.

        Put on the reference rather than the overlay deliberately: the reference
        is the structure the user is already looking at in full detail, and the
        question being asked is "where did this one change", which is a property
        of a position, not of the second model.
        """
        view = self.win.view
        if view is None or self.win.structure is None:
            return
        if not on or self.result is None:
            view.values = None
            view.color_by = self.win._current_color()
        else:
            view.values = self.win.analysis.residue_values_to_atoms(
                self.result.per_residue)
            view.color_by = ColorBy.VALUE
        view.rebuild()
        self.win.viewport.update()
