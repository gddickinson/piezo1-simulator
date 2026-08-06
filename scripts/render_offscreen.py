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


def _orient_to_symmetry_axis(camera, st, edge_on: bool = True) -> None:
    """Point the camera relative to the molecular three-fold axis.

    ``edge_on`` gives the dome *profile* — the axis vertical on screen and the
    camera looking across it, which is the view that actually shows how curved
    the membrane region is. Otherwise the camera looks straight down the axis.
    """
    from piezo1.render.camera import quat_from_axis_angle, quat_multiply
    from piezo1.structure.superpose import detect_c3_axis

    chains = []
    for ch in st.chains:
        m = st.mask_ca() & (st.chain == ch)
        if m.sum() > 300:
            chains.append((st.xyz[m], st.res_seq[m]))
    if len(chains) < 3:
        return
    common = set(chains[0][1].tolist())
    for _, seq in chains[1:]:
        common &= set(seq.tolist())
    common_arr = np.array(sorted(common))
    blocks = [xyz[np.searchsorted(seq, common_arr)].astype(np.float64)
              for xyz, seq in chains[:3]]
    axis = detect_c3_axis(blocks).direction

    # Build the rotation that carries the molecular axis onto screen-up (+y),
    # then express it as the camera quaternion.
    target = np.array([0.0, 1.0, 0.0]) if edge_on else np.array([0.0, 0.0, 1.0])
    v = np.cross(axis, target)
    s = np.linalg.norm(v)
    if s < 1e-8:
        return
    angle = np.arctan2(s, float(np.dot(axis, target)))
    camera.rotation = quat_multiply(
        quat_from_axis_angle(v / s, angle),
        np.array([1.0, 0.0, 0.0, 0.0]),
    )


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
    if view in ("profile", "axial"):
        _orient_to_symmetry_axis(scene.camera, st, edge_on=(view == "profile"))
    elif view == "side":
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
    ap.add_argument("--view", default="side", choices=["front", "side", "top", "profile", "axial"])
    ap.add_argument("--size", nargs=2, type=int, default=[1400, 1000])
    ap.add_argument("-o", "--out", type=Path, default=None)
    args = ap.parse_args()

    out = args.out or Path("docs/img") / f"{args.structure}_{args.style}_{args.color}.png"
    render(args.structure, out, args.style, args.color,
           (args.size[0], args.size[1]), args.view, args.species)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
