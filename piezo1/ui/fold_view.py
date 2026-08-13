"""Drawing and labelling the placed HaloTag fold, in any representation.

Split from :mod:`piezo1.ui.fusion_controller` at the project's 500-line limit
and along a real seam: the controller owns the fusion's *lifecycle* — when the
model is built, what the status line must say, which degrees of freedom the
user may turn — and this module owns the geometry of putting the placed fold
on screen. Nothing here touches the status line or the model; everything is a
pure function of the pose it is handed.
"""

from __future__ import annotations

import numpy as np

__all__ = ["draw_fold_view", "fold_labels"]


def fold_labels(n_atoms: int) -> list:
    """One human-readable label per drawable tag atom, in pose order.

    Falls back to bare indices when the tag file cannot be matched — a label
    is identification, and a wrong residue name is worse than none.
    """
    from ..structure.fusion import load_halotag
    from ..structure.fusion_pose import drawable_mask

    try:
        tag = load_halotag().structure
        base = tag.subset(drawable_mask(tag))
    except (FileNotFoundError, ValueError):
        base = None
    if base is None or base.n_atoms != n_atoms:
        return [f"HaloTag atom {i}" for i in range(n_atoms)]
    return [f"HaloTag {rn}{int(rs)} atom {an}"
            for rn, rs, an in zip(base.res_name, base.res_seq,
                                  base.atom_name, strict=True)]


def draw_fold_view(scene, pose, colours: np.ndarray, style_key: str,
                   name: str) -> bool:
    """The fold in a chosen representation, through the same machinery that
    styles the channel.

    Builds a real :class:`Structure` of the placed tags — one chain per copy,
    so bonds and cartoon traces stay within a tag — and hands it to a
    `MolecularView` whose `color_override` carries the same per-atom colours
    the sphere cloud uses. The contact atoms stay red and the dye stays its
    own colour in every style, because those colours are the visible half of
    the reported numbers, not decoration.

    Returns False when the placed atoms cannot be matched back to the tag
    file, in which case the caller draws the sphere cloud instead: a fold
    silently missing from the screen is worse than one in the wrong style.
    """
    from ..core.structure import Structure
    from ..render.representations import MolecularView, Style
    from ..structure.fusion import load_halotag
    from ..structure.fusion_pose import drawable_mask

    try:
        style = Style(style_key)
        tag = load_halotag().structure
        base = tag.subset(drawable_mask(tag))
    except (FileNotFoundError, ValueError):
        return False
    if base.n_atoms != pose.n_atoms:
        return False

    fields = {f: np.concatenate([getattr(base, f)] * pose.n_tags)
              for f in Structure._ARRAY_FIELDS if f != "xyz"}
    # One chain label per copy: cross-chain bonds are skipped and cartoon
    # traces are per chain, so this is what keeps the three tags separate.
    fields["chain"] = np.concatenate(
        [np.full(base.n_atoms, str(i + 1)) for i in range(pose.n_tags)])
    placed = Structure(
        xyz=pose.coords.reshape(-1, 3).astype(np.float32),
        name="halotag-fold", **fields)
    placed._build_residue_index()

    view = MolecularView(
        scene, placed, name=name, style=style,
        color_override=np.tile(colours, (pose.n_tags, 1)))
    # In ribbon styles the dye would vanish with the side chains; the ligand
    # pass keeps it, in its own colour. In atom styles it is already in the
    # atoms batch, and drawing it twice adds nothing.
    view.ligands_as_spheres = style in (Style.CARTOON, Style.TUBE,
                                        Style.BACKBONE)
    view.rebuild()
    return True
