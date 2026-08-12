"""Ions moving through the pore at a rate the computed current sets.

Round 33 built the permeation physics and left the animation undone, with the
note that the morph clock's discipline applies. It does, and in an instructive
way: the morph refuses a seconds axis because an interpolation between two
endpoints is not a trajectory, while a current genuinely *is* a rate — so here
a time base is meaningful and the honest thing is to name the factor rather
than avoid it.

The factor is large. A single channel passes ~10⁷ ions per second, so anything
watchable runs about a millionfold slow, and :mod:`piezo1.render.flux` computes
that from the solver's own output. The HUD carries the statement whenever
particles are on screen.

**A pore predicted not to conduct shows nothing.** The permeation result is
gated by the wetting verdict, and every deposited human structure is closed.
Drawing a trickle of ions through a shut gate would contradict the project's
own structural result while looking like a demonstration of it — so the closed
case is a first-class outcome here, with the reason shown in place of the
stream.

That outcome is *common*: of the 19 deposited PIEZO1 entries, 17 are refused,
so "nothing happened" is the answer a user meets almost every time — and the
two that do animate are the two worst-resolved models in the catalogue, 11ZC
at 6.0 Å with backbone atoms only and 3JAC at 4.8 Å with 346 unnamed residues.
Two consequences follow, and both were wrong until Round 84c.

First, the cases that do animate have to actually animate. They did not: the
stream is empty on its first frame, the empty upload raised, and
``_on_tick`` responded by unregistering the animation — so a conducting
structure was indistinguishable from a shut one. See
:mod:`piezo1.render.primitives`.

Second, a refusal has to say **where**. Of the 18 entries whose gate can be
located, the narrowest point is beyond it at the cytoplasmic constriction in
16 and above it in the cap in 2; it is in the transmembrane gate in *none* of
them. A message that says only "sterically occluded" invites the reading that
the gate is shut, which is not what was measured.
:func:`piezo1.render.flux.timebase_for_structure` carries the location and the
transmembrane gate's own radius beside it.
"""

from __future__ import annotations

import numpy as np

__all__ = ["IonFluxController", "NAME", "ION_COLOR", "ION_RADIUS"]

NAME = "ion_flux"

#: Deliberately unlike any element colour in `core.structure.ELEMENT_COLORS`:
#: these are not atoms of the model, they are a rate made visible.
ION_COLOR = (0.35, 0.85, 1.0)

#: Angstrom. Larger than a real ion so the particles read as markers rather
#: than as resolved species — the animation shows a rate, not a structure.
ION_RADIUS = 1.9


class IonFluxController:
    """Drives the particle stream and the time-base readout."""

    def __init__(self, window) -> None:
        self.win = window
        self._timebase = None
        self._positions = np.zeros((0,), dtype=float)   # arc length travelled, A
        self._path = np.zeros((0, 3), dtype=float)      # measured pore centres
        self._arc = np.zeros((0,), dtype=float)         # arc length of each centre
        self._running = False
        self._carry = 0.0
        #: Which route the ions may take. `axial` is what every number this
        #: project has recorded was computed on; the lateral options are Liu
        #: et al. 2025's route and are opt-in from the View menu.
        self.pathway = "axial"
        #: Transmembrane potential, volts. None uses the registered default.
        #: Their Figure 5 sweeps 0, -0.1, -0.25 and -0.5 V.
        self.voltage = None

    def set_pathway(self, pathway: str) -> None:
        self.pathway = pathway
        self._restart()

    def set_voltage(self, voltage: float | None) -> None:
        self.voltage = voltage
        self._restart()

    def _restart(self) -> None:
        """Re-run with the new setting, but only if ions are already drawn."""
        if self._running:
            self.show(False)
            self.show(True)

    @property
    def length(self) -> float:
        """Arc length of the conduction path, Angstrom."""
        return float(self._arc[-1]) if len(self._arc) else 0.0

    # ------------------------------------------------------------ lifecycle

    def show(self, on: bool) -> None:
        if not on:
            self.clear()
            return
        if not self._prepare():
            return
        self._running = True
        self.win.viewport.add_animation(self._step)

    def clear(self) -> None:
        self._running = False
        scene = getattr(self.win.viewport, "scene", None)
        if scene is not None:
            for key in list(getattr(scene, "batches", {})):
                if key.startswith(NAME):
                    scene.remove(key)
        hud = getattr(self.win.viewport, "hud", None)
        if hud is not None:
            hud.set_readout("ion_flux", "")
            hud.set_clock("")
        self._positions = np.zeros((0,), dtype=float)
        self._path = np.zeros((0, 3), dtype=float)
        self._arc = np.zeros((0,), dtype=float)

    # -------------------------------------------------------------- setup

    def _prepare(self) -> bool:
        """Compute the pore path and the time base. False if nothing to show."""
        from ..physics.pore_charge import cytosolic_end
        from ..render.flux import timebase_for_structure
        from ..structure.pore import pore_profile
        from ..structure.protomers import protomer_blocks
        from ..structure.superpose import detect_c3_axis

        structure = getattr(self.win, "structure", None)
        if structure is None:
            self._set_status("load a structure first")
            return False

        blocks, _ = protomer_blocks(structure)
        axis = detect_c3_axis(blocks)
        profile = pore_profile(structure, axis)

        # The same gated computation the headless pipeline runs, so the GUI
        # cannot reach a different conclusion about whether this conducts.
        self._timebase = timebase_for_structure(
            structure, profile=profile, pathway=self.pathway,
            voltage=self.voltage)
        self._announce()
        if not self._timebase.conducting:
            return False

        # The stream is drawn over the part of the pore the chosen pathway
        # actually uses, so a lateral route does not animate ions through the
        # closed ends it just excluded.
        from ..physics.conduction_path import conduction_path

        self._set_path(conduction_path(structure, profile, self.pathway).profile,
                       cytosolic_end(structure, axis))
        self._positions = np.array([], dtype=float)
        self._carry = 0.0
        return len(self._arc) > 1

    def _set_path(self, profile, cytosolic: int) -> None:
        """Store the route the ions take: the probe centres, not the axis.

        The pore is measured as the largest sphere that fits at each height,
        and its centre is leashed to within `pore.leash` (8 A) of the symmetry
        axis rather than pinned to it. On 11ZC the fitted centre sits a median
        0.56 A off the axis and up to the full 8 A, and at 11 of 125 heights
        the axis line falls *outside* the sphere that was fitted there — so
        ions drawn on the straight axis would cross the wall of the pore that
        was measured. Interpolating the centres keeps the stream inside it.

        Direction is measured too. `cytosolic_end` reports which end of the
        profile the C-terminal domain sits at, and the stream is ordered to run
        towards it: `detect_c3_axis` fixes a line, not a sign, so half the
        entries would otherwise show inward current flowing out of the cell.
        """
        centers = np.asarray(profile.centers, dtype=float)
        if cytosolic == 0:
            centers = centers[::-1]
        step = np.linalg.norm(np.diff(centers, axis=0), axis=1)
        self._path = centers
        self._arc = np.concatenate([[0.0], np.cumsum(step)])

    def _announce(self) -> None:
        hud = getattr(self.win.viewport, "hud", None)
        if hud is None:
            return
        hud.set_readout("ion_flux", self._timebase.statement())
        if not self._timebase.conducting:
            hud.set_clock("no ions drawn", self._timebase.reason)
        self._set_status(self._timebase.statement())

    def _set_status(self, text: str) -> None:
        setter = getattr(self.win, "_set_status", None)
        if callable(setter):
            setter(text)

    # ------------------------------------------------------------ animation

    def _step(self, dt: float) -> bool:
        """One frame. Returns False to unregister."""
        if not self._running or self._timebase is None:
            return False

        length = self.length
        if length <= 0:
            return False

        # Particles enter at a rate the current sets, and cross in a fixed
        # display time so the stream reads as continuous.
        crossing_seconds = 2.0
        speed = length / crossing_seconds
        self._positions = self._positions + speed * dt
        self._positions = self._positions[self._positions <= length]

        self._carry += self._timebase.particles_per_second * dt
        emitted = int(self._carry)
        if emitted:
            self._carry -= emitted
            self._positions = np.concatenate(
                [self._positions, np.zeros(emitted, dtype=float)])

        self._draw()
        hud = getattr(self.win.viewport, "hud", None)
        if hud is not None:
            hud.set_clock(f"{len(self._positions)} ions in the pore",
                          self._timebase.statement())
        return True

    def points_at(self, arc: np.ndarray) -> np.ndarray:
        """World coordinates for arc-length positions along the pore path."""
        arc = np.atleast_1d(np.asarray(arc, dtype=float))
        if len(self._arc) < 2 or len(arc) == 0:
            return np.zeros((0, 3), dtype=float)
        return np.column_stack([np.interp(arc, self._arc, self._path[:, k])
                                for k in range(3)])

    def _draw(self) -> None:
        scene = getattr(self.win.viewport, "scene", None)
        if scene is None:
            return
        batch = scene.spheres(f"{NAME}:ions")
        points = self.points_at(self._positions)
        # An empty frame is normal — the stream starts with no ions in the pore
        # — and `upload` treats zero rows as "draw nothing" rather than raising.
        batch.upload(points.astype(np.float32),
                     np.full(len(points), ION_RADIUS, np.float32),
                     np.tile(np.array(ION_COLOR, np.float32),
                             (len(points), 1)).reshape(len(points), 3))
