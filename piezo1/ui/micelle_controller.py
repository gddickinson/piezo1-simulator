"""Draws the modelled detergent micelle — Guo & MacKinnon 2017, Figure 4b.

Figure 4b is the one panel of that paper whose *whole content* is something
this project does not hold: the unsharpened cryo-EM map at 6 sigma, showing a
digitonin micelle curved into a dome around the trimer. The picture is the
paper's direct evidence that PIEZO1 bends its surroundings, and it is evidence
because it is a **measurement of the detergent**.

What this draws is a construction from the protein's own coordinates — the
surface at a fixed offset outside the hydrophobic belt. So the status line says
so and cannot be omitted, exactly as the HaloTag fold's does, because a
translucent envelope hugging a channel is precisely the kind of picture that
gets screenshotted and captioned "the micelle".

**What is a measurement here and what is not.** The *thickness* of the shell is
a registered parameter and carries no information. The *curvature* is a sphere
fitted to the belt atoms themselves, which is a property of the protein: on
6B3R it comes out at 9.8 nm against the paper's 10.2 nm idealisation and this
project's own 10.8 nm dome fit, and the status line reports all three so the
one number worth reading is next to what it should be compared with.

The envelope is drawn translucent and back-face culled off, so the protein
stays visible inside it — the published panel's ribbon-inside-mesh reading.
"""

from __future__ import annotations

__all__ = ["MicelleController", "MICELLE_COLOR", "MICELLE_ALPHA"]

#: A warm grey with a hint of yellow — detergent, and deliberately not any of
#: the chain palette's colours or the dome cap's blue, both of which are drawn
#: at the same place.
MICELLE_COLOR = (0.82, 0.78, 0.62)
#: Translucent enough to read the ribbon through. The published panel is a
#: wireframe mesh; a low-alpha surface reads the same way and needs no second
#: rendering path.
MICELLE_ALPHA = 0.20


class MicelleController:
    """Owns the drawn micelle envelope under View -> Micelle density."""

    NAME = "micelle_envelope"

    def __init__(self, window) -> None:
        self.win = window
        self.envelope = None
        self._shown = False

    @property
    def visible(self) -> bool:
        return self._shown

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
        if scene is not None and self._shown:
            scene.remove(self.NAME)
        self._shown = False
        self.envelope = None
        self.win.viewport.update()

    # ------------------------------------------------------------- building

    def _build(self) -> None:
        from ..structure.micelle import build_micelle

        self.clear()
        record = self.win.record
        reference = record.numbering_species if record else "human"
        self.win._set_status("building the modelled micelle envelope…")
        try:
            self.envelope = build_micelle(self.win.structure, reference)
        except Exception as exc:
            self.win._set_status(
                f"micelle failed: {type(exc).__name__}: {exc}")
            return
        self._draw()

    def _draw(self) -> None:
        import numpy as np

        envelope = self.envelope
        scene = self.win.viewport.scene
        colours = np.tile(np.asarray(MICELLE_COLOR, dtype=np.float32),
                          (len(envelope.vertices), 1))
        batch = scene.mesh(self.NAME, transparent=True)
        batch.upload(envelope.vertices, envelope.normals, colours,
                     envelope.faces.ravel(), MICELLE_ALPHA)
        # The envelope encloses the protein, so culling either face of it hides
        # half the surface from inside — the same trap the impostor quads hit.
        batch.cull = False
        self._shown = True
        self.win.viewport.update()
        self.win._set_status(self.status_line())

    # ------------------------------------------------------------- reporting

    def status_line(self) -> str:
        """What must be said whenever the envelope is on screen.

        Not optional. The published panel is a density map and this is a
        construction, and nothing about the picture distinguishes them.
        """
        envelope = self.envelope
        if envelope is None:
            return "no micelle drawn"
        curvature = ""
        if envelope.sphere is not None:
            curvature = (f"belt curvature R = {envelope.sphere.radius / 10:.1f} "
                         f"nm (published idealisation 10.2 nm) · ")
        dome = getattr(self.win, "dome_surface", None)
        measured = getattr(dome, "geometry", None) if dome else None
        if measured is not None:
            curvature += (f"our dome fit "
                          f"{measured.radius_of_curvature / 10:.1f} nm · ")
        return (f"MODELLED MICELLE, NOT THE OBSERVED DENSITY. Figure 4b is the "
                f"unsharpened cryo-EM map at 6 sigma; this is the surface "
                f"{envelope.offset:.0f} A outside the hydrophobic belt "
                f"({envelope.n_belt_atoms} apolar transmembrane side-chain "
                f"atoms). {curvature}The shell thickness is a parameter and "
                f"carries no information; only the curvature is a measurement "
                f"of the protein.")
