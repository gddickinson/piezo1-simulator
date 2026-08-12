"""Draw the pore the profile measures, as the probe spheres it was measured with.

The pore radius has existed since Round 12 and only ever as a two-axis plot:
radius against height. A plot answers *how narrow* and cannot answer *where* —
which is the question a viewer actually has when they are told a structure is
shut. The gate is at z = 4 Å tells you nothing unless you already know where
the origin is and which way the axis points.

**What is drawn is exactly what was measured.** :func:`pore_profile` finds, at
each height, the largest sphere that fits without overlapping any atom. Those
spheres are the measurement, so those spheres are what appear — same centres,
same radii, no smoothing and no surface reconstruction. The controller reads
the profile the Analysis panel already computed rather than computing its own,
so the picture and the plot cannot disagree; if they ever do, one of them is
reading a stale object and that is a bug rather than a finding.

**A probe sphere is free volume, not a wall.** The obvious misreading is that
these spheres are the pore's surface. They are its *complement*: the space
left over. A sphere that looks large is a wide pore, and the protein surface
is somewhere outside it at an unstated distance.

**Three bands, and the bands are two registered numbers.** Red is narrower
than `pore.ion_radius` — a bare ion of that radius does not fit. Amber clears
that but not `pore.constriction_threshold`, the conventional hydrated-ion cut.
Blue clears both. The bands are applied **ascending**, so each slice takes the
highest band it clears: applied in declaration order the `>= 0` pass paints
everything the first colour, which is exactly what Round 76 found the pLDDT
colouring had been doing since it was written.

**Geometry is not a verdict.** A pore wide enough for an ion can still be
non-conducting, because a hydrophobic lumen dewets before it occludes — that
is the whole content of Round 19, and this project's own answer for whether a
structure conducts comes from the wetting model rather than from a radius. So
the status line carries the wetting verdict when the grid is available, and
says the radius alone does not settle it when it is not.
"""

from __future__ import annotations

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["PoreSurfaceController", "radius_colors", "band_index",
           "drawn_slice_mask", "NAME", "NARROW_COLOR", "TIGHT_COLOR",
           "OPEN_COLOR", "BAND_COLORS", "BAND_NAMES", "LINING_COLOR",
           "LINING_RADIUS"]

NAME = "pore_surface"

#: Narrower than a bare ion of the registered radius.
NARROW_COLOR = (0.92, 0.28, 0.22)
#: Passes a bare ion, not a hydrated one.
TIGHT_COLOR = (0.96, 0.72, 0.24)
#: Clears the conventional constriction threshold.
OPEN_COLOR = (0.30, 0.62, 0.95)
#: C-alphas of the residues lining the narrowest slice.
LINING_COLOR = (0.98, 0.92, 0.45)
#: Drawn small, because a C-alpha is not the side chain that does the lining.
LINING_RADIUS = 0.9


#: The three bands, narrowest first. Index into both of these together.
BAND_COLORS = (NARROW_COLOR, TIGHT_COLOR, OPEN_COLOR)
BAND_NAMES = ("narrow", "tight", "open")


def band_index(radius, ion_radius: float | None = None,
               threshold: float | None = None) -> np.ndarray:
    """Which band each slice falls in: 0 narrow, 1 tight, 2 open.

    Bands are applied **ascending** so that each slice ends on the highest
    band it clears. Order matters and has bitten this project once already:
    `plddt_band_colors` applied its bands in declaration order, so the
    ``>= 0.0`` pass repainted every atom and the whole model came out the
    "very low confidence" orange. Written this way the same mistake produces a
    visibly wrong picture in a test rather than a plausible one on screen.

    Both thresholds resolve from the registry at call time, so an override in
    the parameters dialog takes effect on the next redraw. This is the single
    place the bands are decided — the colours and the counts on the status
    line both come through here, so the picture cannot say one thing and the
    caption another.
    """
    if ion_radius is None:
        ion_radius = _P.value("pore.ion_radius")
    if threshold is None:
        threshold = _P.value("pore.constriction_threshold")
    radius = np.asarray(radius, dtype=float)
    index = np.zeros(radius.shape, dtype=np.int64)
    index[radius >= ion_radius] = 1
    index[radius >= threshold] = 2
    return index


def radius_colors(radius, ion_radius: float | None = None,
                  threshold: float | None = None) -> np.ndarray:
    """One colour per slice, from the two registered pore thresholds."""
    bands = band_index(radius, ion_radius, threshold)
    return np.asarray(BAND_COLORS, dtype=np.float32)[bands]


def drawn_slice_mask(radius, leash: float | None = None) -> np.ndarray:
    """Which slices are inside the protein, and so worth drawing.

    **Why anything is dropped at all.** The clearance function has no interior
    maximum: a free probe leaves the pore sideways and grows without bound —
    6188 Å on real coordinates, which is the bulk solvent outside the protein.
    `pore_profile` stops that by tethering the probe centre within `pore.leash`
    of the axis, and that is what makes the number mean "radius of *the pore*".
    Past the last atom the tether is the only thing left bounding it, so the
    profile keeps returning numbers that are really "how much empty space is
    near the axis here" — on 11ZC, five spheres up to 12.2 Å across, drawn as a
    bulb hanging below the channel.

    **The criterion is the method's own parameter.** Once the probe's *radius*
    exceeds the leash, the tether has stopped localising anything: a sphere
    that large centred anywhere in the tethered region swallows the axis and
    its whole neighbourhood. That is the point past which the profile is not
    reporting a lumen, and it needs no new measurement to detect.

    **Only the ends are trimmed, and that distinction is the whole design.**
    An over-leash slice in the *middle* is a genuine wide vestibule with
    protein on both sides of it — 11ZC has 43 of them, and the first version of
    this rule (keep the contiguous run around the bottleneck) cut 71 slices off
    the upper vestibule to get rid of 5 at the bottom. So the leading and
    trailing runs are dropped and nothing else is.

    **This is display only.** The profile object, the bottleneck, the plot and
    every reported number are untouched — a trimmed slice is wider than the
    leash and so can never have been the minimum, unless the whole profile is,
    in which case nothing is trimmed at all.
    """
    if leash is None:
        leash = _P.value("pore.leash")
    radius = np.asarray(radius, dtype=float)
    inside = radius <= leash
    mask = np.ones(radius.shape, dtype=bool)
    if not inside.any():
        # Every slice is over the leash, so there is no end to trim towards:
        # drawing all of it and saying so beats drawing nothing.
        return mask
    first, last = int(np.argmax(inside)), int(len(inside) - 1 - np.argmax(inside[::-1]))
    mask[:first] = False
    mask[last + 1:] = False
    return mask


class PoreSurfaceController:
    """Owns the drawn pore under View -> Pore surface."""

    def __init__(self, window) -> None:
        #: How solid the drawn probe spheres are, 0-1. A probe sphere is the
        #: space left over rather than the wall, and at full opacity it hides
        #: the lining residues that define it — which are what anyone looking
        #: at a pore is looking for.
        self.opacity = 1.0
        self.win = window
        self._names: list[str] = []
        #: Set while waiting for the Analysis panel's pore run to land, so the
        #: toggle can be switched on before the profile exists.
        self.pending = False

    # ----------------------------------------------------------------- state

    @property
    def visible(self) -> bool:
        return bool(self._names)

    @property
    def profile(self):
        """The profile being drawn — the Analysis controller's, never a copy."""
        return self.win.analysis.pore

    def show(self, on: bool) -> None:
        self.pending = False
        if not on:
            self.clear()
            return
        if self.win.structure is None or self.win.viewport.scene is None:
            self.win._set_status("load a structure first")
            return
        if self.profile is None:
            # Computed once, by the Analysis panel's own worker, so there is
            # exactly one pore profile in the application at a time.
            self.pending = True
            self.win._set_status("computing the pore profile to draw it…")
            self.win.analysis.compute_pore()
            return
        self._draw()

    def set_opacity(self, value: float) -> None:
        """Change how solid the drawn pore is, redrawing if it is on screen."""
        self.opacity = max(0.0, min(1.0, float(value)))
        if self.visible:
            self._draw()
            self.win.viewport.update()

    def clear(self) -> None:
        scene = self.win.viewport.scene
        if scene is not None:
            for name in self._names:
                scene.remove(name)
        self._names = []
        self.win.viewport.update()

    def refresh(self) -> None:
        """Redraw from the current profile. Called when a pore run finishes."""
        if not (self.pending or self.visible):
            return
        self.pending = False
        if self.profile is None or self.win.viewport.scene is None:
            return
        self._draw()

    # -------------------------------------------------------------- building

    def drawn_mask(self) -> np.ndarray:
        """The slices this controller draws: the profile minus its escaped ends."""
        profile = self.profile
        if profile is None:
            return np.zeros(0, dtype=bool)
        return drawn_slice_mask(profile.radius)

    def _draw(self) -> None:
        self.clear()
        profile = self.profile
        scene = self.win.viewport.scene
        keep = self.drawn_mask()
        centers = np.asarray(profile.centers, dtype=np.float32)[keep]
        radii = np.asarray(profile.radius, dtype=np.float32)[keep]

        batch = scene.spheres(f"{NAME}:probe")
        batch.alpha = float(self.opacity)
        batch.upload(centers, radii, radius_colors(radii).astype(np.float32))
        self._names.append(f"{NAME}:probe")

        lining = self.lining_coords()
        if len(lining):
            batch = scene.spheres(f"{NAME}:lining")
            batch.upload(lining.astype(np.float32),
                         np.full(len(lining), LINING_RADIUS, np.float32),
                         np.tile(np.array(LINING_COLOR, np.float32),
                                 (len(lining), 1)))
            self._names.append(f"{NAME}:lining")

        self.win.viewport.update()
        self.win._set_status(self.status_line())

    def lining_coords(self) -> np.ndarray:
        """C-alphas of the residues lining the narrowest slice, all protomers.

        The bottleneck is the one place in the profile a reader needs to find,
        and marking it by drawing its probe sphere differently would be a lie
        about a radius. Marking the residues around it instead says the same
        thing without touching the measurement.
        """
        profile = self.profile
        structure = self.win.structure
        if profile is None or structure is None:
            return np.zeros((0, 3), dtype=float)
        lining = profile.bottleneck_lining()
        if not lining:
            return np.zeros((0, 3), dtype=float)
        mask = structure.mask_ca() & np.isin(structure.res_seq,
                                             np.asarray(lining, dtype=np.int64))
        return structure.xyz[mask]

    # ------------------------------------------------------------- reporting

    def band_counts(self) -> dict:
        """How many slices fall in each band, for the status line and tests.

        Counted through `band_index`, the same call that picks the colours, so
        the caption cannot drift from the picture it describes.
        """
        profile = self.profile
        if profile is None:
            return {}
        bands = band_index(profile.radius)[self.drawn_mask()]
        return {name: int((bands == i).sum())
                for i, name in enumerate(BAND_NAMES)}

    def status_line(self) -> str:
        """What is on screen, and the two things a radius does not establish.

        Neither caveat is optional. A probe sphere is the space left over
        rather than the wall, and a lumen wide enough for an ion can still be
        shut — this project's own verdict on that comes from the wetting model,
        not from the radius being drawn.
        """
        profile = self.profile
        if profile is None:
            return "no pore profile"
        counts = self.band_counts()
        ion = _P.value("pore.ion_radius")
        threshold = _P.value("pore.constriction_threshold")
        lining = profile.bottleneck_lining()
        where = (f" lined by {', '.join(str(r) for r in lining[:6])}"
                 if lining else "")
        hydration = self.win.analysis.hydration
        if hydration is not None and getattr(hydration, "available", False):
            verdict = f" · wetting model says {hydration.verdict}"
        else:
            verdict = (" · radius alone does not say conducting: a wide "
                       "hydrophobic lumen dewets, and the wetting model is "
                       "what answers that here")
        keep = self.drawn_mask()
        dropped = int((~keep).sum())
        trimmed = (f" · {dropped} slices past each end NOT drawn: wider than "
                   f"the {_P.value('pore.leash'):.0f} A leash, so the tether "
                   f"rather than the protein is bounding them and they are "
                   f"bulk solvent — the profile and the bottleneck are "
                   f"unchanged" if dropped else "")
        return (f"pore surface: {int(keep.sum())} of {len(profile.z)} probe "
                f"spheres, bottleneck {profile.bottleneck_radius:.2f} A at "
                f"z = {profile.bottleneck_z:.1f} A{where}{trimmed} · "
                f"{counts.get('narrow', 0)} red below the {ion:.1f} A ion, "
                f"{counts.get('tight', 0)} amber below the {threshold:.1f} A "
                f"hydrated cut, {counts.get('open', 0)} blue above it · these "
                f"are the probe spheres that FIT, not the pore wall{verdict}")
