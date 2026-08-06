#!/usr/bin/env python
"""Render a structure to a PNG without opening a window.

Useful for three things: checking the shaders on a headless machine, producing
figures for the documentation, and regression-testing the renderer.

Usage::

    python scripts/render_offscreen.py 8YEZ --style cartoon --color domain
    python scripts/render_offscreen.py 7WLT --style spheres --size 1600 1200
"""

from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import RenderSettings, STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.render.representations import ColorBy, MolecularView, Style  # noqa: E402
from piezo1.render.scene import Scene  # noqa: E402


def render(structure_id: str, out: Path, style: str, color: str,
           size: tuple[int, int], view: str, species: str,
           samples: int = 4) -> Path:
    import moderngl
    from PIL import Image

    path = STRUCTURE_DIR / f"{structure_id}.cif"
    if not path.exists():
        raise SystemExit(f"{path} not found — run scripts/fetch_data.py")

    t0 = time.time()
    st = Structure.from_file(path)
    t_load = time.time() - t0

    ctx = moderngl.create_standalone_context(require=410)
    settings = RenderSettings(samples=samples)
    scene = Scene(ctx, settings)

    mv = MolecularView(scene, st, name=structure_id,
                       style=Style(style), color_by=ColorBy(color))
    mv.set_species(species)
    t0 = time.time()
    mv.rebuild()
    t_build = time.time() - t0

    w, h = size
    scene.resize(w, h)
    # Orient first, then frame: the tight-fit calculation depends on which way
    # the molecule is facing.
    if view == "side":
        scene.camera.orbit(0.0, 0.25)
    elif view == "top":
        scene.camera.orbit(0.0, -0.42)
    scene.camera.frame(st.xyz)

    colour = ctx.texture((w, h), 4, samples=samples)
    depth = ctx.depth_renderbuffer((w, h), samples=samples)
    fbo = ctx.framebuffer(color_attachments=[colour], depth_attachment=depth)
    fbo.use()

    t0 = time.time()
    scene.render()
    ctx.finish()
    t_render = time.time() - t0

    # Resolve multisampling into a plain texture before reading back.
    resolve_tex = ctx.texture((w, h), 4)
    resolve = ctx.framebuffer(color_attachments=[resolve_tex])
    ctx.copy_framebuffer(resolve, fbo)
    data = resolve.read(components=3, alignment=1)

    img = Image.frombytes("RGB", (w, h), data).transpose(Image.FLIP_TOP_BOTTOM)
    out.parent.mkdir(parents=True, exist_ok=True)
    img.save(out)

    ribbon = scene.get(f"{structure_id}:ribbon")
    atoms = scene.get(f"{structure_id}:atoms")
    print(f"{structure_id}: {st.n_atoms} atoms, {st.n_residues} residues")
    print(f"  load {t_load:.2f}s | geometry {t_build:.2f}s | render {t_render*1000:.0f} ms")
    if ribbon:
        print(f"  ribbon vertices: {ribbon.count:,}")
    if atoms:
        print(f"  sphere impostors: {atoms.count:,}")
    print(f"  wrote {out}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("structure", nargs="?", default="8YEZ")
    ap.add_argument("--style", default="cartoon",
                    choices=[s.value for s in Style])
    ap.add_argument("--color", default="domain",
                    choices=[c.value for c in ColorBy])
    ap.add_argument("--species", default="human", choices=["human", "mouse"])
    ap.add_argument("--view", default="side", choices=["front", "side", "top"])
    ap.add_argument("--size", nargs=2, type=int, default=[1400, 1000])
    ap.add_argument("-o", "--out", type=Path, default=None)
    args = ap.parse_args()

    out = args.out or Path("docs/img") / f"{args.structure}_{args.style}_{args.color}.png"
    render(args.structure, out, args.style, args.color,
           (args.size[0], args.size[1]), args.view, args.species)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
