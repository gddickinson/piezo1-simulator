#!/usr/bin/env python
"""Generate every figure used in the README and documentation.

Figures are rendered on a **common scale and a common orientation**, which
matters for the curved-versus-flat comparison: two panels drawn at whatever
zoom each happened to need, one of them upside down, would exaggerate or hide
the very difference the figure exists to show.

Orientation is fixed by the molecular three-fold axis, with its sign chosen so
that the extracellular cap is always up.

Usage::

    python scripts/make_figures.py
    python scripts/make_figures.py --only dome
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import RESOURCE_DIR, RenderSettings, STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.core.annotations import load_annotations  # noqa: E402
from piezo1.render.camera import quat_from_axis_angle  # noqa: E402
from piezo1.render.representations import ColorBy, MolecularView, Style  # noqa: E402
from piezo1.render.scene import Scene  # noqa: E402
from piezo1.structure.superpose import detect_c3_axis  # noqa: E402

OUT = Path("docs/img")
BG = (14, 17, 24)


def protomer_blocks(st: Structure):
    chains = []
    for ch in st.chains:
        m = st.mask_ca() & (st.chain == ch)
        if m.sum() > 300:
            chains.append((st.xyz[m], st.res_seq[m]))
    if len(chains) < 3:
        return None, None
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:]:
        common &= set(seq.tolist())
    arr = np.array(sorted(common))
    return [xyz[np.searchsorted(seq, arr)].astype(np.float64)
            for xyz, seq in chains[:3]], arr


def cap_up_axis(st: Structure, species: str) -> np.ndarray | None:
    """Three-fold axis, signed so the extracellular cap points along +axis."""
    blocks, _ = protomer_blocks(st)
    if blocks is None:
        return None
    axis = detect_c3_axis(blocks).direction
    ann = load_annotations(species)
    cap = next((d for d in ann.domains if d.id == "cap"), None)
    if cap and cap.start:
        m = st.mask_ca() & (st.res_seq >= cap.start) & (st.res_seq <= cap.end)
        if m.sum() > 30:
            centre = np.vstack(blocks).mean(axis=0)
            if np.dot(st.xyz[m].mean(axis=0) - centre, axis) < 0:
                axis = -axis
    return axis


def orient_profile(camera, axis: np.ndarray) -> None:
    """Put ``axis`` along screen-up with the camera looking across it."""
    target = np.array([0.0, 1.0, 0.0])
    v = np.cross(axis, target)
    s = np.linalg.norm(v)
    if s < 1e-8:
        camera.rotation = np.array([1.0, 0.0, 0.0, 0.0])
        return
    angle = np.arctan2(s, float(np.dot(axis, target)))
    camera.rotation = quat_from_axis_angle(v / s, angle)


def render_panel(pdb: str, species: str, size, view: str = "profile",
                 style: Style = Style.CARTOON, color: ColorBy = ColorBy.DOMAIN,
                 distance: float | None = None, ligands: bool = True):
    """Render one panel; returns ``(PIL image, camera)``."""
    import moderngl
    from PIL import Image

    st = Structure.from_file(STRUCTURE_DIR / f"{pdb}.cif")
    ctx = moderngl.create_standalone_context(require=410)
    scene = Scene(ctx, RenderSettings(samples=4))
    mv = MolecularView(scene, st, name=pdb, style=style, color_by=color)
    mv.set_species(species)
    mv.ligands_as_spheres = ligands
    mv.rebuild()

    w, h = size
    scene.resize(w, h)
    axis = cap_up_axis(st, species)
    if view == "profile" and axis is not None:
        orient_profile(scene.camera, axis)
    elif view == "axial" and axis is not None:
        orient_profile(scene.camera, axis)
        scene.camera.orbit(0.0, -0.5)
    scene.camera.frame(st.xyz)
    if distance is not None:
        scene.camera.distance = distance

    colour = ctx.texture((w, h), 4, samples=4)
    depth = ctx.depth_renderbuffer((w, h), samples=4)
    fbo = ctx.framebuffer(color_attachments=[colour], depth_attachment=depth)
    fbo.use()
    scene.render()
    ctx.finish()
    tex = ctx.texture((w, h), 4)
    resolve = ctx.framebuffer(color_attachments=[tex])
    ctx.copy_framebuffer(resolve, fbo)
    img = Image.frombytes("RGB", (w, h), resolve.read(components=3, alignment=1))
    img = img.transpose(Image.FLIP_TOP_BOTTOM)
    cam_distance = scene.camera.distance
    scene.release()
    ctx.release()
    return img, cam_distance


def label(img, lines, x=18, y=14):
    from PIL import ImageDraw
    d = ImageDraw.Draw(img)
    for i, (text, colour) in enumerate(lines):
        d.text((x, y + i * 15), text, fill=colour)
    return img


def figure_dome_comparison() -> None:
    """Curved versus flat, same scale, same orientation."""
    from PIL import Image
    size = (1000, 660)
    # Render once to learn each panel's natural distance, then use the larger
    # for both so the two states are directly comparable.
    _, d1 = render_panel("7WLT", "mouse", size, ligands=False)
    _, d2 = render_panel("11ZC", "mouse", size, ligands=False)
    shared = max(d1, d2)
    a, _ = render_panel("7WLT", "mouse", size, distance=shared, ligands=False)
    b, _ = render_panel("11ZC", "mouse", size, distance=shared, ligands=False)

    label(a, [("CURVED  7WLT", (232, 240, 255)),
              ("radius of curvature 9.7 nm", (154, 163, 178)),
              ("dome depth 4.9 nm", (154, 163, 178))])
    label(b, [("FLATTENED  11ZC", (232, 240, 255)),
              ("radius of curvature 21.6 nm", (154, 163, 178)),
              ("dome depth 3.5 nm", (154, 163, 178))])

    out = Image.new("RGB", (size[0] * 2 + 6, size[1]), BG)
    out.paste(a, (0, 0))
    out.paste(b, (size[0] + 6, 0))
    out.save(OUT / "curved_vs_flat.png")
    print(f"  wrote {OUT / 'curved_vs_flat.png'}  (shared camera distance "
          f"{shared:.0f} A)")


def figure_hero() -> None:
    img, _ = render_panel("8YEZ", "human", (1800, 1050), view="axial")
    img.save(OUT / "hero_human_piezo1.png")
    print(f"  wrote {OUT / 'hero_human_piezo1.png'}")


def figure_profile() -> None:
    img, _ = render_panel("8YEZ", "human", (1400, 900), view="profile")
    label(img, [("Human PIEZO1  8YEZ  3.3 A", (232, 240, 255)),
                ("dome profile, viewed across the three-fold axis",
                 (154, 163, 178))])
    img.save(OUT / "human_profile.png")
    print(f"  wrote {OUT / 'human_profile.png'}")


def figure_domain_key() -> None:
    """A legend strip explaining the domain colours."""
    from PIL import Image, ImageDraw
    ann = load_annotations("human")
    shown = [d for d in ann.domains
             if d.id in ("thu1", "thu5", "thu9", "beam", "anchor",
                         "outer_helix", "cap", "inner_helix", "ctd")]
    w, rowh = 1400, 30
    img = Image.new("RGB", (w, rowh * len(shown) + 16), BG)
    d = ImageDraw.Draw(img)
    for i, dom in enumerate(shown):
        y = 8 + i * rowh
        d.rectangle([14, y + 6, 40, y + 22], fill=dom.color)
        d.text((54, y + 9), f"{dom.name}", fill=(226, 232, 244))
        d.text((520, y + 9), f"human {dom.start}-{dom.end}", fill=(154, 163, 178))
        d.text((700, y + 9), f"mouse {dom.mouse_start}-{dom.mouse_end}",
               fill=(154, 163, 178))
        d.text((900, y + 9), dom.category, fill=(122, 167, 255))
    img.save(OUT / "domain_key.png")
    print(f"  wrote {OUT / 'domain_key.png'}")


def figure_secondary() -> None:
    img, _ = render_panel("8YEZ", "human", (1200, 800), view="profile",
                          color=ColorBy.SECONDARY, ligands=False)
    label(img, [("Secondary structure assigned from C-alpha geometry",
                 (232, 240, 255)),
                ("helix 77%  ·  strand 10%  ·  coil 13%", (154, 163, 178)),
                ("red helix · yellow strand · grey coil", (154, 163, 178))])
    img.save(OUT / "secondary_structure.png")
    print(f"  wrote {OUT / 'secondary_structure.png'}")


FIGURES = {
    "hero": figure_hero,
    "profile": figure_profile,
    "dome": figure_dome_comparison,
    "key": figure_domain_key,
    "secondary": figure_secondary,
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", choices=sorted(FIGURES), default=None)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    todo = [args.only] if args.only else list(FIGURES)
    for name in todo:
        print(f"{name}:")
        FIGURES[name]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
