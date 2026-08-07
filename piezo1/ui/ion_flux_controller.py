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
        self._positions = np.zeros((0,), dtype=float)
        self._axis = None
        self._span = (0.0, 0.0)
        self._running = False
        self._carry = 0.0

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

    # -------------------------------------------------------------- setup

    def _prepare(self) -> bool:
        """Compute the pore axis and the time base. False if nothing to show."""
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
        self._axis = axis
        self._span = (float(profile.z.min()), float(profile.z.max()))

        # The same gated computation the headless pipeline runs, so the GUI
        # cannot reach a different conclusion about whether this conducts.
        self._timebase = timebase_for_structure(structure, profile=profile)
        self._announce()
        if not self._timebase.conducting:
            return False

        self._positions = np.array([], dtype=float)
        self._carry = 0.0
        return True

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

        low, high = self._span
        length = high - low
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

        self._draw(low)
        hud = getattr(self.win.viewport, "hud", None)
        if hud is not None:
            hud.set_clock(f"{len(self._positions)} ions in the pore",
                          self._timebase.statement())
        return True

    def _draw(self, low: float) -> None:
        scene = getattr(self.win.viewport, "scene", None)
        if scene is None or self._axis is None:
            return
        batch = scene.spheres(f"{NAME}:ions")
        if len(self._positions) == 0:
            batch.upload(np.zeros((0, 3), np.float32), np.zeros(0, np.float32),
                         np.zeros((0, 3), np.float32))
            return

        direction = self._axis.direction / np.linalg.norm(self._axis.direction)
        points = (self._axis.point[None, :]
                  + (low + self._positions)[:, None] * direction[None, :])
        batch.upload(points.astype(np.float32),
                     np.full(len(points), ION_RADIUS, np.float32),
                     np.tile(np.array(ION_COLOR, np.float32),
                             (len(points), 1)))
