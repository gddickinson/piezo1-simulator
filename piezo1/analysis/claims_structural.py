"""Claim recomputations that need downloaded coordinates.

Split out of ``claims.py`` at the project's length limit, along a real seam:
everything here reloads a deposited structure and re-runs a pipeline over it,
which is why every claim it serves costs ``medium`` or ``slow`` and skips on a
fresh clone. What stays in ``claims.py`` recomputes from registered parameters
and recorded results alone.

The dependency runs one way — ``claims`` imports this, never the reverse. Round
78 lost a session to a cycle of exactly the shape a shared loader invites.
"""

from __future__ import annotations

import numpy as np

__all__ = ["_structure", "_blocks", "_tm_surface", "_dome_radius",
           "_gating_overlap", "_pore_bottleneck", "_selectivity",
           "_wetting_score", "_within_position_variance"]


def _structure(pdb: str):
    from ..config import STRUCTURE_DIR
    from ..core.structure import Structure
    path = STRUCTURE_DIR / f"{pdb.upper()}.cif"
    if not path.exists():
        raise FileNotFoundError(f"{pdb} not downloaded")
    return Structure.from_file(path)


def _blocks(structure):
    from ..structure.protomers import protomer_blocks
    return protomer_blocks(structure)


def _tm_surface(structure, species: str) -> np.ndarray:
    """Kept as a name because claims and tests call it; the definition moved.

    Round 83 needed one definition of the dome surface rather than two, so
    this delegates to :func:`piezo1.structure.geometry.tm_surface_points`.
    """
    from ..structure.geometry import tm_surface_points

    return tm_surface_points(structure, species)[0]


def _dome_radius() -> float:
    from ..structure.geometry import measure_dome
    st = _structure("7WLT")
    blocks, _ = _blocks(st)
    return measure_dome(blocks, _tm_surface(st, "mouse")).radius_of_curvature / 10.0


def _gating_overlap() -> float:
    """Overlap of the lowest symmetric mode with the observed transition."""
    from ..physics.anm import ANM
    from ..structure.superpose import match_protomers

    curved, flat = _structure("7WLT"), _structure("7WLU")
    cb, cr = _blocks(curved)
    fb, fr = _blocks(flat)
    common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))

    def resample(st, residues):
        out = []
        for chain in st.chains:
            mask = st.mask_ca() & (st.chain == chain)
            if mask.sum() < 300:
                continue
            index = {int(r): i for i, r in enumerate(st.res_seq[mask])}
            xyz = st.xyz[mask]
            out.append(np.array([xyz[index[r]] for r in residues], dtype=float))
        return out[:3]

    cb, fb = resample(curved, common), resample(flat, common)
    match = match_protomers(cb, fb)
    fb = [fb[i] for i in match.order]

    anm = ANM.from_trimer(cb, cutoff=15.0, spring="inverse_square").build()
    modes = anm.calc_modes(n_modes=30)
    anm.label_symmetry(modes)

    from ..structure.superpose import kabsch
    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    fitted = (np.vstack(fb) - centroid) @ rotation.T + translation
    displacement = (fitted - np.vstack(cb)).ravel()

    # `overlap` returns the whole spectrum at once; take the best A mode.
    overlaps = np.abs(np.asarray(modes.overlap(displacement), dtype=float))
    symmetric = np.array([s == "A" for s in modes.symmetry])
    return float(overlaps[symmetric].max()) if symmetric.any() else 0.0


def _pore_bottleneck(pdb: str) -> float:
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    st = _structure(pdb)
    blocks, _ = _blocks(st)
    return float(pore_profile(st, detect_c3_axis(blocks))
                 .bottleneck_radius)


def _selectivity(mode: str | None) -> float:
    """P_Cl/P_Na from the dilution potential, with and without the pore's charge.

    ``mode=None`` is the uncharged baseline, which is not one: the pore is
    already weakly cation-selective from size exclusion, and quoting a charged
    result without it would credit the charge with what geometry had done.
    """
    from ..physics.pore_charge import cytosolic_end, pore_charge
    from ..physics.selectivity import measure_selectivity
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis

    st = _structure("11ZC")
    blocks, _ = _blocks(st)
    axis = detect_c3_axis(blocks)
    profile = pore_profile(st, axis)
    charge = (None if mode is None else
              pore_charge(st, profile, axis, mode=mode, species="mouse").density)
    result = measure_selectivity(profile, fixed_charge=charge,
                                 cytosolic_index=cytosolic_end(st, axis))
    return float(result.p_anion_over_cation)


def _wetting_score(pdb: str) -> float:
    from .hydration import load_grid, predict_wetting
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    grid = load_grid()
    if not grid.available:
        raise FileNotFoundError("CHAP grid not downloaded")
    st = _structure(pdb)
    blocks, _ = _blocks(st)
    profile = pore_profile(st, detect_c3_axis(blocks))
    return float(predict_wetting(st, profile, grid).score)


def _within_position_variance() -> float:
    """Fraction of the mechanical ddG variance that is *within* position.

    The Round 26 criterion, measured on the multiply-substituted positions.
    Round 7's model gave 0.049 here; the per-contact model must exceed 0.20.
    """
    import collections

    from .substitution import variance_decomposition
    from .variant_impact import VariantImpactModel
    from ..core.annotations import load_annotations
    from ..core.sequence import human_sequence, human_to_mouse, mouse_to_human
    from ..structure.superpose import kabsch, match_protomers
    from ..structure.protomers import protomer_blocks

    curved, flat = _structure("7WLT"), _structure("7WLU")
    _cb, cr = protomer_blocks(curved)
    _fb, fr = protomer_blocks(flat)
    common = np.array(sorted(set(cr.tolist()) & set(fr.tolist())))

    def resample(st):
        out = []
        for chain in st.chains:
            mask = st.mask_ca() & (st.chain == chain)
            if mask.sum() < 300:
                continue
            index = {int(r): i for i, r in enumerate(st.res_seq[mask])}
            xyz = st.xyz[mask]
            if all(r in index for r in common):
                out.append(np.array([xyz[index[r]] for r in common], float))
        return out[:3]

    cb, fb = resample(curved), resample(flat)
    fb = [fb[i] for i in match_protomers(cb, fb).order]
    rotation, translation, centroid = kabsch(np.vstack(fb), np.vstack(cb))
    displacement = (((np.vstack(fb) - centroid) @ rotation.T + translation)
                    - np.vstack(cb))

    reference = human_sequence()
    sequence = {}
    for residue in common:
        human = mouse_to_human(int(residue))
        if human and 1 <= human <= len(reference):
            sequence[int(residue)] = reference[human - 1]

    by_position = collections.defaultdict(list)
    for variant in load_annotations("human").variants:
        if not (variant.residue and variant.wt_aa and variant.mut_aa):
            continue
        if len(variant.wt_aa) != 1 or len(variant.mut_aa) != 1:
            continue
        if not variant.mut_aa.isalpha():
            continue
        if variant.label not in {v.label for v in by_position[variant.residue]}:
            by_position[variant.residue].append(variant)
    multi = {k: v for k, v in by_position.items() if len(v) > 1}

    model = VariantImpactModel(coords=np.vstack(cb),
                               residues=np.tile(common, 3),
                               gating_vector=displacement, sequence=sequence)
    positions, values = [], []
    for human_residue, variants in multi.items():
        mouse_residue = human_to_mouse(human_residue)
        if mouse_residue is None:
            continue
        for variant in variants:
            prediction = model.predict(mouse_residue, variant.wt_aa,
                                       variant.mut_aa)
            if prediction.modelled and np.isfinite(prediction.gating_cost_change):
                positions.append(human_residue)
                values.append(prediction.gating_cost_change)
    return variance_decomposition(positions, values).within_fraction


# --------------------------------------------------------------------------
# Guo & MacKinnon 2017 — the paper the dome model comes from
# --------------------------------------------------------------------------

def _guo_panel(key: str, *path):
    """One number out of a replicated panel, measured on the paper's own entry."""
    from .guo2017 import replicate

    node = replicate(key, structure=_structure("6B3R"))["result"]
    for step in path:
        node = node[step]
    return float(node)


def _guo_beam_angle() -> float:
    """Beam-to-pore-axis angle on 6B3R. The paper states 'about 60 degrees'."""
    return _guo_panel("4a", "beam_angle_deg")


def _guo_blade_share() -> float:
    """Fraction of the trimer's non-planarity carried by the distal blade."""
    return _guo_panel("4a", "blade", "share_of_arrangement")


def _guo_pore_offset() -> float:
    """How much wider our pore radii are than the published HOLE ones, A.

    A methodological difference reported rather than absorbed. Zero would mean
    the profiler had been fitted to the paper.
    """
    return _guo_panel("6b", "mean_offset_A")


def _guo_patch_interaction() -> float:
    """Screened interaction between the cap and loop charged patches, k_BT."""
    return _guo_panel("4-S1", "interaction_kT")
