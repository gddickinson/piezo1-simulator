"""Draw the calcium an open channel puts around its own pore exit.

Round 32 predicted 113.8 µM at the modelled tag against a 0.2 µM sensor Kd, and
concluded the sensor is saturated whenever its channel opens — so puncta
brightness reports labelling stoichiometry and open probability, not calcium
amplitude. That conclusion is about a *distance*, and it has only ever been
stated as one.

**The iso-surfaces are exactly spheres.** The screened buffered-diffusion
Green's function is spherically symmetric about a point source, so nothing is
being idealised into a sphere for drawing; the model has no other shape.

**And the ones that carry the conclusion do not fit in the picture.** At
11ZC's 2.43 pA the sensor is still 90% occupied at **119 nm** and half-occupied
at **372 nm**, against a channel about 30 nm across. Drawn, they are two
enormous shells with a speck inside, and the speck is the protein. So they are
reported as numbers and not drawn — the same rule
:mod:`piezo1.ui.dome_controller` applies to the far-field membrane footprint,
and for the same reason: a surface that swamps the thing it is drawn around
stops being a picture of anything.

What *is* drawn is the near field at the protein's own scale: shells at decade
concentrations, out to a stated multiple of the structure's own size. Those
show the gradient where a tag could actually sit, which is the part a picture
can say better than a number.

**The model does not know the protein is there.** A point source in free
solution, screened by a uniform buffer. The shells therefore pass straight
through the channel and the membrane, because the equation they come from does
too. Clipping them where they meet an atom would imply a boundary condition
that was never imposed.

**Nothing is borrowed.** `report_tags` substitutes 11ZC's current when the
loaded entry is shut, labelled. A picture cannot carry that label convincingly:
a cloud drawn around 8YEZ reads as 8YEZ's, whatever the caption says. So a shut
structure draws **nothing**, and the status line reports the reason — which is
itself Round 34's result, that no deposited human PIEZO1 entry conducts.
"""

from __future__ import annotations

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["NanodomainController", "NAME", "SHELL_CONCENTRATIONS",
           "SCENE_MULTIPLE", "MIN_SHELL_A", "SOURCE_COLOR", "NEAR_COLOR",
           "FAR_COLOR", "SATURATION_OCCUPANCY", "KD_OCCUPANCY"]

NAME = "nanodomain"

#: Shell concentrations in molar, highest first. Decades and half-decades, so
#: the spacing on screen reads as the ~1/r falloff it is.
SHELL_CONCENTRATIONS = (1e-2, 3e-3, 1e-3, 3e-4, 1e-4, 3e-5, 1e-5)

#: A shell is drawn only out to this multiple of the structure's own extent.
#: Beyond it the surface fills the viewport and the protein becomes a speck,
#: which is what the two occupancy surfaces would do at every plausible
#: current — hence they are reported rather than drawn.
SCENE_MULTIPLE = 1.6

#: Below this the shell is inside the pore mouth and draws over the source.
MIN_SHELL_A = 6.0

#: Hottest drawn shell — nearest the mouth, highest concentration.
NEAR_COLOR = (1.00, 0.42, 0.28)
#: Coolest drawn shell.
FAR_COLOR = (0.38, 0.66, 1.00)
#: The point source: the measured cytosolic mouth of the pore.
SOURCE_COLOR = (1.00, 0.95, 0.55)

#: `Nanodomain.saturated` is occupancy > 0.9, so the reported saturation radius
#: is that same threshold and the picture cannot disagree with the boolean.
SATURATION_OCCUPANCY = 0.9
#: Half occupancy is the Kd by definition.
KD_OCCUPANCY = 0.5


class NanodomainController:
    """Owns the drawn calcium field under View -> Calcium nanodomain."""

    def __init__(self, window) -> None:
        self.win = window
        self.model = None
        #: Concentration (M) -> radius (Angstrom), every shell, drawn or not.
        self.shells: dict[float, float] = {}
        #: Occupancy -> radius (Angstrom) for the two surfaces that carry the
        #: conclusion. Kept whether or not they fit, because the status line
        #: reports them either way.
        self.occupancy_radii: dict[float, float] = {}
        self._names: list[str] = []
        self.pending = False

    # ----------------------------------------------------------------- state

    @property
    def visible(self) -> bool:
        return bool(self._names)

    def show(self, on: bool) -> None:
        self.pending = False
        if not on:
            self.clear()
            return
        if self.win.structure is None or self.win.viewport.scene is None:
            self.win._set_status("load a structure first")
            return
        if self.win.analysis.pore is None:
            self.pending = True
            self.win._set_status("measuring the pore to find its exit…")
            self.win.analysis.compute_pore()
            return
        self._build()

    def clear(self) -> None:
        scene = self.win.viewport.scene
        if scene is not None:
            for name in self._names:
                scene.remove(name)
        self._names = []
        self.model = None
        self.shells = {}
        self.occupancy_radii = {}
        self.win.viewport.update()

    def refresh(self) -> None:
        if not (self.pending or self.visible):
            return
        self.pending = False
        if self.win.analysis.pore is None or self.win.viewport.scene is None:
            return
        self._build()

    # ------------------------------------------------------------ the physics

    def source_point(self) -> np.ndarray:
        """Where the calcium enters the cytosol: the pore's cytosolic mouth.

        Which end that is gets *measured*, by the same call
        :mod:`piezo1.physics.selectivity` uses — the sign of that answer is the
        whole content of a permeability ratio, so neither module assumes it.
        """
        from ..physics.pore_charge import cytosolic_end

        profile = self.win.analysis.pore
        end = cytosolic_end(self.win.structure, profile.axis)
        return np.asarray(profile.centers[end], dtype=float)

    def unitary_current(self) -> float:
        """The loaded structure's own current, in amps. Never a borrowed one."""
        from ..physics.permeation import solve_pnp

        result = solve_pnp(self.win.analysis.pore, self.win.analysis.hydration)
        return abs(float(result.current))

    def scene_radius(self, origin) -> float:
        """How far the structure itself reaches from the source, in Angstrom."""
        xyz = self.win.structure.xyz
        return float(np.linalg.norm(xyz - np.asarray(origin), axis=1).max())

    def _build(self) -> None:
        from ..physics.nanodomain import (Nanodomain, distance_for_occupancy,
                                          saturation)

        self.clear()
        try:
            current = self.unitary_current()
        except Exception as exc:
            self.win._set_status(
                f"nanodomain failed: {type(exc).__name__}: {exc}")
            return
        if current <= 0.0:
            self.win._set_status(self.closed_line())
            return

        fraction = _P.value("nanodomain.calcium_current_fraction")
        self.model = Nanodomain(current_A=current, calcium_fraction=fraction,
                                distance_m=1e-9)
        calcium = self.model.calcium_current_A

        # A concentration and an occupancy are the same statement through the
        # 1:1 binding curve, so both inversions go through the one solver
        # rather than a second bisection written here.
        for concentration in SHELL_CONCENTRATIONS:
            metres = distance_for_occupancy(
                float(saturation(concentration)), calcium)
            if np.isfinite(metres):
                self.shells[concentration] = metres * 1e10
        for occupancy in (SATURATION_OCCUPANCY, KD_OCCUPANCY):
            self.occupancy_radii[occupancy] = \
                distance_for_occupancy(occupancy, calcium) * 1e10
        self._draw()

    # -------------------------------------------------------------- building

    def drawable(self, origin=None) -> list[tuple[float, float]]:
        """(concentration, radius) for the shells that fit around the protein.

        The budget is the point of this method. A shell wider than
        `SCENE_MULTIPLE` times the structure's own reach fills the viewport and
        turns the protein into a speck — which is precisely what the two
        occupancy surfaces do, at every current this project can produce. They
        are reported instead, and the status line says they were.
        """
        if not self.shells or self.win.structure is None:
            return []
        origin = self.source_point() if origin is None else origin
        limit = SCENE_MULTIPLE * self.scene_radius(origin)
        return [(c, r) for c, r in sorted(self.shells.items(), reverse=True)
                if MIN_SHELL_A <= r <= limit]

    def _shell_color(self, index: int, total: int) -> tuple:
        if total <= 1:
            return NEAR_COLOR
        t = index / (total - 1)
        near = np.asarray(NEAR_COLOR)
        far = np.asarray(FAR_COLOR)
        return tuple(near + t * (far - near))

    def _draw(self) -> None:
        from ..render.geometry_builders import build_sphere

        scene = self.win.viewport.scene
        origin = self.source_point()
        shells = self.drawable(origin)

        # Outermost first: the transparent pass draws in submission order, so
        # a near shell submitted first would be overdrawn by the far one.
        for i, (concentration, radius) in enumerate(reversed(shells)):
            colour = self._shell_color(len(shells) - 1 - i, len(shells))
            mesh = build_sphere(origin, radius, color=colour)
            name = f"{NAME}:{concentration:g}"
            batch = scene.mesh(name, two_sided=True, transparent=True)
            batch.upload(mesh.positions, mesh.normals, mesh.colors,
                         mesh.indices, 0.16)
            self._names.append(name)

        batch = scene.spheres(f"{NAME}:source")
        batch.upload(origin[None, :].astype(np.float32),
                     np.array([2.2], np.float32),
                     np.array([SOURCE_COLOR], np.float32))
        self._names.append(f"{NAME}:source")

        self.win.viewport.update()
        self.win._set_status(self.status_line())

    # ------------------------------------------------------------- reporting

    def closed_line(self) -> str:
        """Why a shut structure draws nothing, said as the result it is."""
        return ("calcium nanodomain: this structure carries no current, so the "
                "nanodomain is exactly zero and nothing is drawn — that is "
                "Round 34's result, not a failure: no deposited human PIEZO1 "
                "entry conducts. Load 11ZC, the one open-like entry, to see "
                "one. A current is deliberately NOT borrowed here, as the "
                "report does: a cloud drawn around a shut channel reads as "
                "that channel's whatever the caption says.")

    def status_line(self) -> str:
        """The near field drawn, the far field reported, and why it is not drawn."""
        if self.model is None:
            return "no calcium nanodomain"
        shells = self.drawable()
        if shells:
            span = (f"{len(shells)} shells drawn, "
                    f"{shells[0][0] * 1e6:.0f} uM at {shells[0][1] / 10:.1f} nm "
                    f"in to {shells[-1][0] * 1e6:.0f} uM at "
                    f"{shells[-1][1] / 10:.1f} nm out")
        else:
            span = "no shell falls within the structure's own extent"
        saturated = self.occupancy_radii.get(SATURATION_OCCUPANCY, float("nan"))
        half = self.occupancy_radii.get(KD_OCCUPANCY, float("nan"))
        reach = self.scene_radius(self.source_point()) / 10.0
        kd = _P.value("nanodomain.sensor_kd")
        tag = ""
        model = getattr(self.win.fusion, "model", None)
        if model is not None and getattr(model, "pore_exit", None) is not None:
            from ..physics.nanodomain import calcium_at, saturation

            distance = float(model.pore_exit_distances()[0]) * 1e-9
            concentration = float(calcium_at(distance, self.model.calcium_current_A))
            tag = (f" · at the modelled tag, {distance * 1e9:.1f} nm out, "
                   f"{concentration * 1e6:.1f} uM and the sensor is "
                   f"{float(saturation(concentration)):.0%} occupied")
        return (f"calcium nanodomain: {self.model.current_A * 1e12:.2f} pA, "
                f"{self.model.calcium_fraction:.0%} of it calcium · {span} · "
                f"the sensor is 90% occupied out to {saturated / 10:.0f} nm "
                f"and half-occupied at its {kd * 1e6:.1f} uM Kd out to "
                f"{half / 10:.0f} nm — against a channel reaching {reach:.0f} "
                f"nm, so NEITHER is drawn: they would swamp the protein, and "
                f"the whole channel sits inside both{tag} · point source in "
                f"free solution — the shells pass through the protein because "
                f"the equation does too")
