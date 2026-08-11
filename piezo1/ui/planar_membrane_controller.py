"""One protomer in a flat bilayer — Guo & MacKinnon 2017, Figure 4a.

Figure 4a draws a single subunit from the side with two grey lines across it:
"approximate locations of planar membrane interfaces". The trimer in 4b gets no
such lines, and the text says why — "a single subunit of Piezo removed from the
trimer can be positioned reasonably well into the plane of a lipid membrane.
However, the detergent micelle containing a trimer is curved into a dome". That
contrast is the pivot of the paper's whole discussion.

This draws it, and — the part a picture cannot do on its own — it draws the
*same* planes for the trimer on request, so the claim can be seen failing
rather than merely asserted. The planes are the least-squares fit to the
transmembrane band from :mod:`piezo1.structure.planarity`, which is also what
the numbers in the status line come from, so the picture and the measurement
cannot be of different things.

**The honest part is the residual, not the lines.** Any point set has a
best-fit plane; drawing it proves nothing. What the status line carries is how
far the band actually departs from it, and the slab thickness that would be
needed to contain the whole thing — against a real bilayer's 36 A. On 6B3R one
protomer needs about 42 A and the trimer about 60 A, which is the difference
Figure 4a is drawing.
"""

from __future__ import annotations

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["PlanarMembraneController", "LEAFLET_COLOR", "MIDPLANE_COLOR"]

#: The two leaflet interfaces — Figure 4a's grey rules.
LEAFLET_COLOR = (0.62, 0.64, 0.70)
#: The mid-plane, drawn fainter. Not in the published panel; it is what the
#: fit actually returns, and the two interfaces are it plus and minus half a
#: bilayer.
MIDPLANE_COLOR = (0.45, 0.55, 0.72)


class PlanarMembraneController:
    """Owns the flat-membrane planes under View -> Planar membrane."""

    def __init__(self, window) -> None:
        self.win = window
        self.comparison = None
        self.chain: str | None = None
        self._names: list[str] = []

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
        self.comparison = None
        self.win.viewport.update()

    def use_trimer(self, on: bool) -> None:
        """Fit the planes to all three protomers instead of one.

        The control that makes the panel a claim rather than a drawing: the
        same construction on the assembly, where it visibly does not work.
        """
        self._trimer = bool(on)
        if self.visible:
            self._build()

    # ------------------------------------------------------------- building

    def _build(self) -> None:
        from ..structure.planarity import planarity

        was_trimer = getattr(self, "_trimer", False)
        self.clear()
        self._trimer = was_trimer
        record = self.win.record
        reference = record.numbering_species if record else "human"
        try:
            self.comparison = planarity(self.win.structure, reference)
        except Exception as exc:
            self.win._set_status(
                f"planar membrane failed: {type(exc).__name__}: {exc}")
            return
        self._draw()

    def _fit_and_points(self):
        """The plane to draw and the points it was fitted to."""
        from ..structure.geometry import tm_surface_by_chain

        comparison = self.comparison
        record = self.win.record
        reference = record.numbering_species if record else "human"
        by_chain, resolved = tm_surface_by_chain(self.win.structure, reference)
        wanted = sorted(resolved)

        if getattr(self, "_trimer", False):
            points = np.vstack([
                np.array([helices[i] for i in wanted if i in helices])
                for helices in by_chain.values()])
            self.chain = None
            return comparison.trimer, points

        # One protomer: the first, so repeated runs draw the same one.
        chain = sorted(comparison.per_protomer)[0]
        self.chain = chain
        points = np.array([by_chain[chain][i] for i in wanted
                           if i in by_chain[chain]])
        return comparison.per_protomer[chain], points

    def _draw(self) -> None:
        from ..render.geometry_builders import build_disc

        fit, points = self._fit_and_points()
        scene = self.win.viewport.scene
        normal = np.asarray(fit.normal, dtype=float)
        centre = np.asarray(fit.point, dtype=float)

        # Size the discs to the band they are drawn across, plus a margin, so
        # they frame the protomer rather than filling the viewport.
        in_plane = points - centre
        in_plane = in_plane - np.outer(in_plane @ normal, normal)
        radius = float(np.linalg.norm(in_plane, axis=1).max()) * 1.15

        half = 0.5 * _P.value("membrane.thickness") * 10.0
        for name, offset, colour, alpha in (
                ("membrane_upper", +half, LEAFLET_COLOR, 0.30),
                ("membrane_lower", -half, LEAFLET_COLOR, 0.30),
                ("membrane_mid", 0.0, MIDPLANE_COLOR, 0.12)):
            mesh = build_disc(radius, color=colour, axis=normal,
                              origin=centre + normal * offset)
            batch = scene.mesh(name, transparent=True)
            batch.upload(mesh.positions, mesh.normals, mesh.colors,
                         mesh.indices, alpha)
            # A plane seen from below is as informative as from above.
            batch.cull = False
            self._names.append(name)

        self.win.viewport.update()
        self.win._set_status(self.status_line(fit))

    # ------------------------------------------------------------- reporting

    def status_line(self, fit=None) -> str:
        comparison = self.comparison
        if comparison is None:
            return "no membrane plane drawn"
        if fit is None:
            fit = (comparison.trimer if getattr(self, "_trimer", False)
                   else comparison.per_protomer[sorted(comparison.per_protomer)[0]])
        bilayer = _P.value("membrane.thickness") * 10.0
        subject = ("all three protomers" if getattr(self, "_trimer", False)
                   else f"protomer {self.chain}")
        return (
            f"planar membrane fitted to {subject}: RMS departure "
            f"{fit.rmsd:.1f} A, worst {fit.max_deviation:.1f} A, so a slab of "
            f"{fit.thickness_needed():.0f} A would be needed to contain the "
            f"band against a real bilayer's {bilayer:.0f} A · "
            f"protomer {comparison.protomer_rmsd:.1f} A vs trimer "
            f"{comparison.trimer.rmsd:.1f} A "
            f"({comparison.ratio:.1f}x) — Figure 4a's claim, and it holds here "
            f"only because the arms tilt "
            f"{comparison.mean_tilt_deg:.0f} deg out of the pore's plane · "
            f"the planes are a least-squares fit, and every point set has one: "
            f"read the residual, not the lines")
