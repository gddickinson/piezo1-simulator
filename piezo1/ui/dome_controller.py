"""Draw the dome that was measured, and the membrane footprint it implies.

The dome is this project's central geometric claim — 9.7 nm of curvature on
7WLT against a published 10.2 — and until now it was four numbers in a status
bar. Nothing showed *where* the fitted surface actually sits relative to the
protein, which is the one thing that would let a viewer disbelieve it: a sphere
fitted to the wrong atoms produces a perfectly reasonable radius and a surface
that visibly misses the transmembrane helices.

**Two surfaces, and both are defined by the measurement.**

The **fitted cap** is a sphere through the mid-points of the transmembrane
helices, drawn out to the footprint radius. The **flat disc** at the rim height
is that cap's own projection onto the membrane plane — the surface the dome
would relax to if it flattened completely. The gap between them *is* the excess
area the gating model is built on, which is otherwise a number in a panel with
no picture attached.

**What is deliberately not drawn: the far-field footprint.** The obvious third
surface is the bilayer relaxing back to flat outside the rim, from the
linearised Helfrich solution. The first version of this drew it, and the
picture is why it does not now. PIEZO1's cap meets the membrane at a slope of
about **1.9** — a 63 degree contact angle — and the linear theory is a
small-slope expansion. Continued from that rim it plunges 158 Å over a 526 Å
skirt that swamps the protein, and Round 18 already measured that it
overestimates this footprint by **3.65x**. A confident, detailed, wrong surface
is the exact failure this project exists to avoid, so the status line states
the footprint radius as a number and draws nothing beyond it.

**Why a surface of revolution is honest.** The dome is C3-symmetric, not
axisymmetric, so a surface of revolution is already an idealisation. What makes
it defensible is that the measurement it is drawing made exactly the same
assumption: :func:`piezo1.structure.geometry.measure_dome` fits a sphere and a
radial profile. Drawing anything more detailed would show structure the number
does not have.
"""

from __future__ import annotations

import numpy as np

__all__ = ["DomeController", "CAP_COLOR", "PROJECTION_COLOR"]

#: The fitted cap — the part that is a measurement.
CAP_COLOR = (0.30, 0.52, 0.80)
#: The flat projection the excess area is measured against.
PROJECTION_COLOR = (0.62, 0.64, 0.70)


class DomeController:
    """Owns the drawn dome surface under View -> Dome surface."""

    def __init__(self, window) -> None:
        self.win = window
        self.geometry = None
        self._names: list[str] = []

    # ----------------------------------------------------------------- state

    @property
    def visible(self) -> bool:
        return bool(self._names)

    def show(self, on: bool) -> None:
        if not on:
            self.clear()
            return
        if self.win.structure is None or self.win.viewport.scene is None:
            self.win._set_status("load a structure first")
            return
        self._build()

    def clear(self) -> None:
        scene = self.win.viewport.scene
        if scene is not None:
            for name in self._names:
                scene.remove(name)
        self._names = []
        self.geometry = None
        self.win.viewport.update()

    # ------------------------------------------------------------- building

    def _build(self) -> None:
        from ..structure.geometry import measure_dome, tm_surface_points

        self.clear()
        record = self.win.record
        reference = record.numbering_species if record else "human"
        blocks = self.win._mode_blocks
        if not blocks:
            self.win._set_status("need three protomers to find the dome axis")
            return
        try:
            points, resolved = tm_surface_points(self.win.structure, reference)
            if len(points) < 4:
                raise ValueError("too few transmembrane helices resolved")
            self.geometry = measure_dome(blocks, points)
        except Exception as exc:
            self.win._set_status(f"dome surface failed: {type(exc).__name__}: {exc}")
            return
        self._draw(len(resolved))

    def _height(self, geometry):
        """``z(r)`` of the fitted cap, in the frame the sphere was fitted in.

        The cap is the lower or upper half of the fitted sphere depending on
        which side of its centre the protein sits, and that is decided by the
        data rather than assumed: the sign is taken from where the measured
        radial profile actually is.
        """
        axis = geometry.axis
        centre = np.asarray(geometry.sphere.center, dtype=float)
        radius = float(geometry.sphere.radius)
        direction = axis.direction / np.linalg.norm(axis.direction)
        centre_height = float((centre - axis.point) @ direction)
        # Which half of the fitted sphere the protein is on, taken from the
        # measured radial profile rather than assumed from the frame: the
        # canonical frame fixes which way is cytosolic, not which way the dome
        # bulges, and a flattened entry can sit either side of its own centre.
        profile = geometry.profile
        heights = np.asarray(profile.z, float)[np.asarray(profile.valid)]
        apex = float(np.mean(heights)) if heights.size else centre_height
        sign = 1.0 if apex > centre_height else -1.0

        def height(r):
            inside = np.clip(radius ** 2 - np.asarray(r, float) ** 2, 0.0, None)
            return centre_height + sign * np.sqrt(inside)

        return height, sign, centre_height

    def _draw(self, n_helices: int) -> None:
        from ..render.geometry_builders import build_membrane_mesh

        geometry = self.geometry
        scene = self.win.viewport.scene
        axis = geometry.axis
        direction = axis.direction / np.linalg.norm(axis.direction)
        height, sign, centre_height = self._height(geometry)
        cap_radius = min(float(geometry.sphere.radius),
                         float(geometry.footprint_radius) or
                         float(geometry.sphere.radius))

        cap = build_membrane_mesh(height, cap_radius, color=CAP_COLOR,
                                  axis=direction, origin=axis.point)
        self._upload("dome_cap", scene, cap, 0.45)

        rim = float(height(cap_radius))
        flat = build_membrane_mesh(lambda r: np.full_like(np.asarray(r, float), rim),
                                   cap_radius, color=PROJECTION_COLOR,
                                   axis=direction, origin=axis.point)
        self._upload("dome_projection", scene, flat, 0.22)

        self.win.viewport.update()
        self.win._set_status(self.status_line(n_helices))

    def _upload(self, name: str, scene, mesh, alpha: float) -> None:
        batch = scene.mesh(name, transparent=True)
        batch.upload(mesh.positions, mesh.normals, mesh.colors,
                     mesh.indices, alpha)
        self._names.append(name)

    # ------------------------------------------------------------- reporting

    def status_line(self, n_helices: int | None = None) -> str:
        """What must be said whenever the surface is on screen.

        Half of what is drawn is a fit to the coordinates and half is the
        solution of an equation with two registered parameters in it. A viewer
        looking at one continuous surface has no way to tell, so this is not
        optional — the same rule the HaloTag fold's status line follows.
        """
        geometry = self.geometry
        if geometry is None:
            return "no dome surface"
        helices = "" if n_helices is None else f"{n_helices} TM helices · "
        return (f"dome surface: {helices}fitted cap R_c "
                f"{geometry.radius_of_curvature / 10:.2f} nm, depth "
                f"{geometry.dome_depth / 10:.2f} nm (sphere RMSE "
                f"{geometry.notes.get('sphere_rmse', float('nan')):.2f} A) · "
                f"grey disc is its flat projection, and the gap between them "
                f"is the {geometry.excess_area / 100:.0f} nm^2 of excess area "
                f"· the membrane beyond {geometry.footprint_radius / 10:.1f} nm "
                f"is NOT drawn: linear Helfrich theory overestimates it 3.65x "
                f"at this contact slope")



