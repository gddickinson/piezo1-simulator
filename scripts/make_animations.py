#!/usr/bin/env python
"""Render the library of PIEZO1 animations.

Each animation shows one mechanism, captioned with what it is and what is
measured, so the file stands on its own when pasted into a talk.

    python scripts/make_animations.py                 # all of them
    python scripts/make_animations.py --only gating
    python scripts/make_animations.py --list
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from piezo1.config import RESOURCE_DIR, STRUCTURE_DIR  # noqa: E402
from piezo1.core import Structure  # noqa: E402
from piezo1.core.annotations import load_annotations  # noqa: E402
from piezo1.physics.anm import ANM  # noqa: E402
from piezo1.render.animation import (AnimationSpec, Animator,  # noqa: E402
                                     ease_in_out, ping_pong)
from piezo1.render.representations import ColorBy, MolecularView, Style  # noqa: E402
from piezo1.structure.morph import morph, prepare_endpoints  # noqa: E402
from piezo1.structure.superpose import detect_c3_axis  # noqa: E402

OUT = Path("docs/anim")

#: GIF size controls, applied to every animation. See Animator.save.
SAVE_OPTS = {"scale": 0.72, "colors": 128}


# --------------------------------------------------------------------------
# Shared helpers
# --------------------------------------------------------------------------

def load(pdb: str) -> Structure:
    return Structure.from_file(STRUCTURE_DIR / f"{pdb}.cif")


def protomers(st: Structure):
    chains = []
    for c in st.chains:
        m = st.mask_ca() & (st.chain == c)
        if m.sum() > 300:
            chains.append((st.xyz[m], st.res_seq[m]))
    common = set(chains[0][1].tolist())
    for _, s in chains[1:3]:
        common &= set(s.tolist())
    arr = np.array(sorted(common))
    return ([x[np.searchsorted(s, arr)].astype(float) for x, s in chains[:3]],
            arr)


def cap_up_axis(st: Structure, species: str):
    blocks, _ = protomers(st)
    axis = detect_c3_axis(blocks)
    ann = load_annotations(species)
    cap = next((d for d in ann.domains if d.id == "cap"), None)
    if cap and cap.start:
        m = st.mask_ca() & (st.res_seq >= cap.start) & (st.res_seq <= cap.end)
        if m.sum() > 30:
            centre = np.vstack(blocks).mean(axis=0)
            if np.dot(st.xyz[m].mean(axis=0) - centre, axis.direction) < 0:
                axis.direction = -axis.direction
    return axis


def orient_profile(camera, axis) -> None:
    from piezo1.render.camera import quat_from_axis_angle
    target = np.array([0.0, 1.0, 0.0])
    v = np.cross(axis.direction, target)
    s = np.linalg.norm(v)
    if s < 1e-8:
        return
    camera.rotation = quat_from_axis_angle(
        v / s, np.arctan2(s, float(np.dot(axis.direction, target))))


def atom_site_map(st: Structure, common: np.ndarray) -> np.ndarray:
    """Map every atom onto its protomer-block C-alpha site."""
    per = len(common)
    chains = [c for c in st.chains
              if (st.mask_ca() & (st.chain == c)).sum() > 300][:3]
    out = np.zeros(st.n_atoms, dtype=np.int64)
    for p, c in enumerate(chains):
        sel = st.chain == c
        pos = np.clip(np.searchsorted(common, st.res_seq[sel]), 0, per - 1)
        out[sel] = p * per + pos
    return np.clip(out, 0, 3 * per - 1)


# --------------------------------------------------------------------------
# Animations
# --------------------------------------------------------------------------

def anim_gating_morph(fmt: str = "gif", size=(960, 720),
                      out: Path | None = None) -> Path:
    """The gating transition: curved dome flattening under tension.

    ``size`` and ``out`` exist so the README's smaller copy comes from here
    rather than from someone downscaling a file by hand. It is committed —
    `docs/anim/` is git-ignored, so a GitHub reader would otherwise see a
    broken image — and a committed figure no script can rebuild is a fossil.
    """
    a, b = load("7WLT"), load("7WLU")
    ab, ar = protomers(a)
    bb, br = protomers(b)
    start, end, common, info = prepare_endpoints(ab, ar, bb, br)
    traj = morph(start, end, n_frames=61, method="restrained")

    spec = AnimationSpec(
        name="gating_morph", n_frames=72, size=size, fps=18,
        caption="PIEZO1 gating: the dome flattens under membrane tension",
        subcaption=f"7WLT -> 7WLU  ·  {info['n_common_residues']} residues  ·  "
                   f"{info['endpoint_rmsd']:.0f} A RMSD  ·  interpolation, not a trajectory")
    with Animator(spec) as an:
        view = MolecularView(an.scene, a, name="m", style=Style.CARTOON,
                             color_by=ColorBy.DOMAIN)
        view.set_species("mouse")
        view.ligands_as_spheres = False
        view.rebuild()
        orient_profile(an.scene.camera, cap_up_axis(a, "mouse"))
        an.scene.camera.frame(np.vstack([traj.frames[0], traj.frames[-1]]))
        base = a.xyz.copy()
        amap = atom_site_map(a, common)

        def step(i: int, t: float) -> None:
            delta = (traj.at(t) - traj.frames[0])[amap]
            view.update_coords(base + delta.astype(np.float32))

        an.run(step, schedule=ping_pong(spec.n_frames, hold=5))
        return an.save(out or OUT / f"gating_morph.{fmt}", **SAVE_OPTS)


def anim_normal_mode(mode_index: int = 2, fmt: str = "gif") -> Path:
    """The lowest symmetric elastic-network mode — the gating coordinate."""
    st = load("7WLT")
    blocks, common = protomers(st)
    anm = ANM.from_trimer(blocks, cutoff=15.0).build()
    modes = anm.calc_modes(n_modes=12)
    anm.label_symmetry(modes)
    sym = modes.symmetry[mode_index]
    meaning = ("three-fold symmetric: CAN couple to isotropic tension"
               if sym == "A" else
               "degenerate E pair: symmetry FORBIDS coupling to tension")

    spec = AnimationSpec(
        name=f"mode_{mode_index + 1}", n_frames=48, size=(900, 700), fps=20,
        caption=f"Elastic-network mode {mode_index + 1}  ·  symmetry {sym}",
        subcaption=f"{meaning}  ·  collectivity {modes.collectivity(mode_index):.2f}")
    with Animator(spec) as an:
        view = MolecularView(an.scene, st, name="m", style=Style.TUBE,
                             color_by=ColorBy.VALUE)
        view.set_species("mouse")
        view.ligands_as_spheres = False
        amap = atom_site_map(st, common)
        view.values = np.linalg.norm(modes.vectors[mode_index], axis=1)[amap]
        view.rebuild()
        orient_profile(an.scene.camera, cap_up_axis(st, "mouse"))
        base = st.xyz.copy()
        disp = modes.mode(mode_index, amplitude=22.0)
        an.scene.camera.frame(np.vstack([base + disp[amap], base - disp[amap]]))

        def step(i: int, t: float) -> None:
            view.update_coords(base + (disp[amap] * np.sin(2 * np.pi * t)).astype(np.float32))

        an.run(step, schedule=np.linspace(0.0, 1.0, spec.n_frames))
        return an.save(OUT / f"mode_{mode_index + 1}.{fmt}", **SAVE_OPTS)


def anim_ligand_site(site: str = "yoda1_pocket", fmt: str = "gif") -> Path:
    """Orbit a ligand or lipid site with its residues highlighted."""
    st = load("8YEZ")
    ann = load_annotations("human")
    group = ann.group(site)
    if group is None:
        raise SystemExit(f"unknown site {site!r}")
    mask = np.isin(st.res_seq, np.asarray(group.residues, dtype=np.int32))

    spec = AnimationSpec(
        name=site, n_frames=60, size=(900, 700), fps=20, spin_degrees=360.0,
        caption=group.label,
        subcaption=(" ".join(f"{d['human_aa']}{d['human']}" for d in group.detail)
                    + f"  ·  evidence: {group.evidence}"))
    with Animator(spec) as an:
        view = MolecularView(an.scene, st, name="m", style=Style.CARTOON,
                             color_by=ColorBy.DOMAIN)
        view.rebuild()
        n = int(mask.sum())
        if n:
            batch = an.scene.spheres("site")
            batch.upload(st.xyz[mask], st.vdw_radii()[mask] * 1.1,
                         np.tile(np.array([1.0, 0.84, 0.2], np.float32), (n, 1)),
                         np.ones(n, np.float32))
        # Frame on the site, but keep enough context to see where it sits.
        centre = st.xyz[mask].mean(axis=0) if n else st.center
        an.scene.camera.frame(st.xyz)
        an.scene.camera.pivot = centre
        an.scene.camera.distance *= 0.55
        an.run(lambda i, t: None, schedule=np.linspace(0, 1, spec.n_frames))
        return an.save(OUT / f"site_{site}.{fmt}", **SAVE_OPTS)


def anim_lipid_contacts(fmt: str = "gif") -> Path:
    """The resolved pore lipid and the residues it touches."""
    st = load("8YEZ")
    lig = st.mask_ligands()
    if not lig.any():
        raise SystemExit("no ligands in 8YEZ")
    from piezo1.analysis.interactions import detect_interactions
    prot = st.mask_protein() & ~st.hetero
    inter = detect_interactions(st, mask_a=prot, mask_b=lig,
                                min_sequence_separation=0)
    contact_res = sorted({i.res_i if not lig[i.atom_i] else i.res_j
                          for i in inter})
    mask = np.isin(st.res_seq, np.asarray(contact_res, dtype=np.int32)) & prot

    spec = AnimationSpec(
        name="lipid_contacts", n_frames=60, size=(900, 700), fps=20,
        spin_degrees=360.0,
        caption="Resolved lipid (L9Q) and its protein contacts",
        subcaption=f"{len(inter)} contacts to {len(contact_res)} residues  ·  "
                   f"detected from coordinates, heavy atoms only")
    with Animator(spec) as an:
        view = MolecularView(an.scene, st, name="m", style=Style.CARTOON,
                             color_by=ColorBy.DOMAIN)
        view.ligands_as_spheres = True
        view.rebuild()
        n = int(mask.sum())
        if n:
            b = an.scene.spheres("contacts")
            b.upload(st.xyz[mask], st.vdw_radii()[mask] * 0.9,
                     np.tile(np.array([0.36, 0.92, 0.66], np.float32), (n, 1)),
                     np.zeros(n, np.float32))
        an.scene.camera.frame(st.xyz)
        an.scene.camera.pivot = st.xyz[lig].mean(axis=0)
        an.scene.camera.distance *= 0.5
        an.run(lambda i, t: None, schedule=np.linspace(0, 1, spec.n_frames))
        return an.save(OUT / f"lipid_contacts.{fmt}", **SAVE_OPTS)


def anim_variant(residue: int = 2456, fmt: str = "gif") -> Path:
    """A disease variant in its structural context."""
    st = load("8YEZ")
    ann = load_annotations("human")
    variants = ann.variants_at(residue)
    label = variants[0].label if variants else f"residue {residue}"
    pheno = variants[0].phenotype if variants else ""
    mask = (st.res_seq == residue) & ~st.hetero

    from piezo1.analysis.interactions import detect_interactions
    contacts = detect_interactions(st, min_sequence_separation=3).involving(residue)
    partners = sorted({i.res_j if i.res_i == residue else i.res_i
                       for i in contacts} - {residue})
    pmask = np.isin(st.res_seq, np.asarray(partners, dtype=np.int32)) & ~st.hetero

    spec = AnimationSpec(
        name=f"variant_{residue}", n_frames=60, size=(900, 700), fps=20,
        spin_degrees=360.0,
        caption=f"{label}  ·  {pheno[:60]}",
        subcaption=f"{len(contacts)} contacts to {len(partners)} residues; "
                   f"salt bridge to E2117 of the neighbouring protomer")
    with Animator(spec) as an:
        view = MolecularView(an.scene, st, name="m", style=Style.CARTOON,
                             color_by=ColorBy.DOMAIN)
        view.ligands_as_spheres = False
        view.rebuild()
        for m, colour, scale, name in ((pmask, (0.45, 0.72, 1.0), 0.8, "partners"),
                                       (mask, (1.0, 0.30, 0.35), 1.1, "variant")):
            n = int(m.sum())
            if n:
                b = an.scene.spheres(name)
                b.upload(st.xyz[m], st.vdw_radii()[m] * scale,
                         np.tile(np.array(colour, np.float32), (n, 1)),
                         np.zeros(n, np.float32))
        an.scene.camera.frame(st.xyz)
        an.scene.camera.pivot = st.xyz[mask].mean(axis=0)
        an.scene.camera.distance *= 0.42
        an.run(lambda i, t: None, schedule=np.linspace(0, 1, spec.n_frames))
        return an.save(OUT / f"variant_{residue}.{fmt}", **SAVE_OPTS)


ANIMATIONS = {
    "gating": (anim_gating_morph, "Curved dome flattening (7WLT -> 7WLU)"),
    # Committed to docs/img because the README shows it and docs/anim is
    # git-ignored. Same code path, smaller frame.
    "readme": (lambda fmt="gif": anim_gating_morph(
        fmt, size=(640, 480),
        out=Path(__file__).resolve().parent.parent / "docs" / "img"
        / "gating_morph_small.gif"),
        "The gating morph at README size, into docs/img"),
    "mode": (anim_normal_mode, "Lowest symmetric elastic-network mode"),
    "yoda1": (lambda fmt="gif": anim_ligand_site("yoda1_pocket", fmt),
              "Yoda1 binding pocket"),
    "pip2": (lambda fmt="gif": anim_ligand_site("pip2_cluster", fmt),
             "PIP2-binding lysine cluster"),
    "gate": (lambda fmt="gif": anim_ligand_site("hydrophobic_gate", fmt),
             "Transmembrane hydrophobic gate"),
    "lipid": (anim_lipid_contacts, "Resolved pore lipid and its contacts"),
    "variant": (anim_variant, "R2456H in structural context"),
}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", choices=sorted(ANIMATIONS))
    ap.add_argument("--format", default="gif", choices=["gif", "mp4"])
    ap.add_argument("--list", action="store_true")
    ap.add_argument("--scale", type=float, default=0.72,
                    help="downscale factor for GIF output")
    ap.add_argument("--colors", type=int, default=128,
                    help="GIF palette size")
    args = ap.parse_args()
    SAVE_OPTS["scale"] = args.scale
    SAVE_OPTS["colors"] = args.colors

    if args.list:
        for k, (_, desc) in sorted(ANIMATIONS.items()):
            print(f"  {k:9s} {desc}")
        return 0

    OUT.mkdir(parents=True, exist_ok=True)
    todo = [args.only] if args.only else list(ANIMATIONS)
    for key in todo:
        fn, desc = ANIMATIONS[key]
        print(f"{key}: {desc}")
        try:
            path = fn(fmt=args.format)
            size = path.stat().st_size / 1e6
            print(f"   wrote {path}  ({size:.1f} MB)")
        except Exception as exc:
            print(f"   FAILED {type(exc).__name__}: {exc}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
