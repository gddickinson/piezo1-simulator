#!/usr/bin/env python
"""Figures for the two models the app draws on top of a structure.

Separate from `make_figures.py` because those render one deposited entry and
these add content to the scene — a grafted prediction, a placed tag — which is
the whole point of the picture. Same OpenGL path the application uses, so what
appears here is what a user sees.

    python scripts/make_model_figures.py            # both
    python scripts/make_model_figures.py --only hybrid
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import RenderSettings, STRUCTURE_DIR  # noqa: E402
from piezo1.core.structure import Structure  # noqa: E402
from piezo1.render.colormaps import plddt_band_colors  # noqa: E402
from piezo1.render.representations import (ColorBy, MolecularView,  # noqa: E402
                                           Style)
from piezo1.render.scene import Scene  # noqa: E402

OUT = Path(__file__).resolve().parent.parent / "docs" / "img"
SIZE = (1400, 950)


def _oriented_camera(scene, structure, coords):
    """Axis vertical, camera looking across it — the view that shows the dome."""
    from piezo1.render.camera import quat_from_axis_angle
    from piezo1.structure.superpose import detect_c3_axis

    blocks = []
    for chain in structure.chains:
        mask = structure.mask_ca() & (structure.chain == chain)
        if mask.sum() > 300:
            blocks.append((structure.xyz[mask], structure.res_seq[mask]))
    if len(blocks) >= 3:
        common = set(blocks[0][1].tolist())
        for _, seq in blocks[1:]:
            common &= set(seq.tolist())
        shared = np.array(sorted(common))
        trio = [xyz[np.searchsorted(seq, shared)].astype(np.float64)
                for xyz, seq in blocks[:3]]
        axis = detect_c3_axis(trio).direction
        target = np.array([0.0, 1.0, 0.0])
        v = np.cross(axis, target)
        s = np.linalg.norm(v)
        if s > 1e-8:
            scene.camera.rotation = quat_from_axis_angle(
                v / s, np.arctan2(s, float(np.dot(axis, target))))
    scene.camera.frame(coords)


def _render(build, coords_for_framing, name: str, captions):
    """Run ``build(scene, structure)``, render, caption and save."""
    import moderngl
    from PIL import Image, ImageDraw

    ctx = moderngl.create_standalone_context(require=410)
    scene = Scene(ctx, RenderSettings(samples=4))
    w, h = SIZE
    scene.resize(w, h)

    structure = build(scene)
    _oriented_camera(scene, structure, coords_for_framing())

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

    draw = ImageDraw.Draw(img)
    for i, (text, rgb) in enumerate(captions):
        draw.text((18, 14 + i * 16), text, fill=rgb)

    path = OUT / f"{name}.png"
    img.save(path)
    print(f"  wrote {path}")
    scene.release()
    ctx.release()


# --------------------------------------------------------------- the figures

def figure_hybrid() -> None:
    """The full-length model: what was measured, and what was predicted.

    The point of the picture is that the two are not the same colour. The
    experimental core is grey; the grafted blade carries AlphaFold's own
    confidence bands, so a viewer sees *where* the prediction is weak rather
    than being told a single number about it.
    """
    from piezo1.structure.hybrid import build_hybrid_model
    from piezo1.ui.hybrid_controller import (EXPERIMENTAL_COLOR, SEAM_COLOR,
                                             SEAM_RADIUS)

    host = Structure.from_file(STRUCTURE_DIR / "8YEZ.cif")
    model = build_hybrid_model(host)

    def build(scene):
        xyz = np.asarray(model.xyz, np.float32)
        colors = np.tile(np.float32(EXPERIMENTAL_COLOR), (len(xyz), 1))
        colors[model.predicted] = plddt_band_colors(model.plddt[model.predicted])
        scene.spheres("hybrid").upload(
            xyz, np.full(len(xyz), 1.7, np.float32), colors,
            np.zeros(len(xyz), np.float32))
        at_seam = model.res_seq == model.seam_residue
        seam = np.asarray(model.xyz[at_seam].mean(axis=0), np.float32)
        scene.spheres("seam").upload(
            seam[None, :], np.float32([SEAM_RADIUS * 2.0]),
            np.float32([SEAM_COLOR]), np.float32([1.0]))
        return host

    n_pred = int(model.predicted.sum())
    confident = float(model.confident_prediction.sum()) / n_pred
    _render(build, lambda: np.asarray(model.xyz), "hybrid_model", [
        ("Full-length PIEZO1, one protomer: measured core, predicted blade",
         (232, 240, 255)),
        (f"grey = cryo-EM (8YEZ), residues {model.seam_residue}-2521",
         (168, 176, 190)),
        (f"coloured = AlphaFold, residues 1-{model.seam_residue - 1}, in its own",
         (168, 176, 190)),
        (f"pLDDT bands: blue >90, cyan 70-90, yellow 50-70, orange <50",
         (168, 176, 190)),
        (f"Only {confident:.0%} of the graft clears 70.", (168, 176, 190)),
        (f"pink = the seam: the models agree to {model.overlap_rmsd:.1f} A there,",
         (255, 150, 170)),
        (f"and differ by {model.global_rmsd:.0f} A over the region they share.",
         (255, 150, 170)),
        ("The grafted blade is a PREDICTION, not a measurement.",
         (255, 190, 120)),
    ])


def figure_halotag_fold() -> None:
    """The HaloTag fusion drawn as the real fold, in one arbitrary orientation."""
    from piezo1.structure.fusion import build_fusion, load_halotag
    from piezo1.structure.fusion_pose import pose_for_display
    from piezo1.ui.fusion_controller import (CONTACT_COLOR, DYE_COLOR,
                                             SEAM_COLOR, TAG_ATOM_SCALE,
                                             TAG_COLOR)

    host = Structure.from_file(STRUCTURE_DIR / "8YEZ.cif")
    tag = load_halotag()
    model = build_fusion(host, tag)
    pose = pose_for_display(host, model, tag)

    def build(scene):
        # Uniform, not chain colouring: the tag's orange sits 0.10 from the
        # chain palette's orange, so a per-protomer picture makes the modelled
        # tag look like part of the experimental trimer — the one thing this
        # figure must not do.
        view = MolecularView(scene, host, name="8YEZ", style=Style.CARTOON,
                             color_by=ColorBy.UNIFORM)
        view.rebuild()
        colours = np.tile(np.float32(TAG_COLOR), (pose.n_atoms, 1))
        colours[pose.ligand] = DYE_COLOR
        colours[pose.touching & pose.body] = CONTACT_COLOR
        scene.spheres("tag").upload(
            pose.coords.reshape(-1, 3).astype(np.float32),
            np.tile(pose.radii * TAG_ATOM_SCALE, pose.n_tags).astype(np.float32),
            np.tile(colours, (pose.n_tags, 1)),
            np.zeros(pose.n_atoms * pose.n_tags, np.float32))
        seams = np.asarray(pose.seams, np.float32)
        colour = np.tile(np.float32(SEAM_COLOR), (len(seams), 1))
        scene.cylinders("seam").upload(
            seams[:, 0], seams[:, 1], np.full(len(seams), 1.4, np.float32),
            colour, colour)
        return host

    frame = np.vstack([host.xyz, pose.coords.reshape(-1, 3)])
    _render(build, lambda: frame, "halotag_fold", [
        ("HaloTag on each cytosolic C-terminus, drawn as its real fold (6U32)",
         (232, 240, 255)),
        ("grey = PIEZO1 (8YEZ)   orange = the tag   pink = its bound TMR dye",
         (168, 176, 190)),
        (f"yellow = the linker, {pose.linker_gap:.0f} A of it, unresolved",
         (168, 176, 190)),
        (f"{pose.meta['clear_spins']} of {pose.meta['spins_sampled']} tag "
         f"orientations clear the channel on this entry", (168, 176, 190)),
        ("Position is modelled; the spin about the linker is UNDETERMINED.",
         (255, 190, 120)),
        ("This is one draw of many.", (255, 190, 120)),
    ])


FIGURES = {"hybrid": figure_hybrid, "halotag": figure_halotag_fold}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    parser.add_argument("--only", choices=sorted(FIGURES))
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    for name in ([args.only] if args.only else list(FIGURES)):
        print(f"{name}:")
        FIGURES[name]()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
