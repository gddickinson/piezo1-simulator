"""The computation behind each replicable panel of Guo & MacKinnon 2017.

Split from :mod:`piezo1.analysis.guo2017`, which is the *registry* — the list
of panels, what each shows, and whether it can be reproduced at all. This
module is the other half: one callable per replicable panel, each returning a
plain dict.

Every function takes the structure it is measured on rather than loading one,
except where the panel is about a sequence and no structure is involved. That
keeps the registry able to run a panel against 7WLT or 11ZC and see what
changes, which is often more informative than reproducing 6B3R.

Nothing here decides whether a panel *can* be replicated; that judgement is
data in the registry, so it can be read without running anything.

Figures 1 to 4 are here — the reconstruction, the topology and the curved
micelle. Figures 6 and 7, the pore and the dome mechanism, are in
:mod:`piezo1.analysis.guo2017_mechanism`: split at the project's length limit
and along the paper's own seam, since everything there is the argument the
paper is *for* rather than the structure it is built on.
"""

from __future__ import annotations

import numpy as np

__all__ = ["panel_2cd", "panel_2ab", "panel_3a", "panel_3_supplements",
           "panel_4a", "panel_4b", "panel_4c", "panel_4_supplement",
           "CED_ACIDIC_PATCH", "LOOP_BASIC_PATCH"]


#: Figure 4-figure supplement 1's two patches, mouse numbering.
CED_ACIDIC_PATCH = (2257, 2258, 2264)
LOOP_BASIC_PATCH = (1269, 1761, 1762)


# --------------------------------------------------------------------------
# Figure 2 — the reconstruction
# --------------------------------------------------------------------------

def panel_2cd(structure, reference: str = "mouse", **kw) -> dict:
    """Ribbon views of the trimer, top and side, chains coloured separately.

    Nothing to compute — this is a rendering. What is returned is what the
    renderer needs and what a reader should be able to check: that there are
    three protomers, that they are related by a real three-fold, and how much
    of the chain is actually drawn.
    """
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis

    blocks, chains = protomer_blocks(structure)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}
    axis = detect_c3_axis(blocks)
    ca = structure.mask_ca()
    return {
        "n_protomers": len(blocks),
        "protomer_chains": list(chains),
        "c3_angle_deg": axis.angle_deg,
        "c3_rmsd_A": axis.rmsd,
        "residues_modelled_per_protomer": int(blocks[0].shape[0]),
        "residue_range": [int(structure.res_seq[ca].min()),
                          int(structure.res_seq[ca].max())],
        "published": {"residues_modelled": 1518, "of_total": 2547,
                      "note": "Guo & MacKinnon state 1518 of 2547 residues "
                              "modelled with the N-terminal 576 missing"},
        "render": {"colour_by": "chain", "views": ["top", "side"],
                   "style": "cartoon"},
    }


def panel_2ab(structure, reference: str = "mouse", **kw) -> dict:
    """Simulated projections standing in for the 2D class averages."""
    from .projection import project_views

    views = project_views(structure)
    return {
        "views": {name: {"shape": list(p.image.shape),
                         "pixel_size_A": p.pixel_size,
                         "resolution_A": p.resolution,
                         "extent_nm": p.extent_A / 10.0,
                         "scale_bar_10nm_px": p.scale_bar_pixels()}
                  for name, p in views.items()},
        "is_experimental": False,
        "caveat": next(iter(views.values())).caveat,
    }


# --------------------------------------------------------------------------
# Figure 3 — topology
# --------------------------------------------------------------------------

def panel_3a(structure=None, reference: str = "mouse", **kw) -> dict:
    """The topology cartoon: 38 helices in nine 4-TM units plus the pore."""
    from .topology import build_topology

    topology = build_topology(reference, structure)
    return topology.as_dict()


def panel_3_supplements(structure=None, reference: str = "mouse", **kw) -> dict:
    """Hydropathy, and the 4-TM repeat the paper infers from it."""
    from .hydropathy import (annotated_hydropathy, hydropathy_profile,
                             repeat_periodicity, threshold_scan)

    profile = hydropathy_profile(reference=reference)
    repeat = repeat_periodicity(reference=reference)
    return {
        "window": profile.window,
        "sequence_length": profile.sequence_length,
        "hydropathy_range": [float(profile.value.min()),
                             float(profile.value.max())],
        "annotated_helices": annotated_hydropathy(profile, reference),
        "threshold_scan": threshold_scan(profile, reference),
        "repeat": {
            "period": repeat.period, "n_units": repeat.n_units,
            "n_helices": repeat.n_helices, "phase": repeat.phase,
            "loop_by_phase": list(repeat.loop_by_phase),
            "contrast_residues": repeat.contrast,
            "control_mean": repeat.control_mean,
            "control_sd": repeat.control_sd, "z": repeat.z,
            "supported": repeat.supported, "summary": repeat.summary()},
        "published": {"n_units": 9, "n_helices": 38,
                      "note": "'a Piezo subunit contains nine 4-TM units "
                              "altogether ... giving a total of 36 TM helices "
                              "... plus two C-terminal TM helices'"},
    }


# --------------------------------------------------------------------------
# Figure 4 — the curved micelle and the surface
# --------------------------------------------------------------------------

def panel_4a(structure, reference: str = "mouse", **kw) -> dict:
    """A protomer fits a plane; the trimer does not."""
    from ..structure.planarity import beam_angle, blade_dependence

    split = blade_dependence(structure, reference)
    result = split.full
    out = {
        "protomer_plane_rmsd_A": result.protomer_rmsd,
        "arrangement_rmsd_A": result.arrangement_rmsd,
        "trimer_plane_rmsd_A": result.trimer.rmsd,
        "decomposition_residual_A": result.decomposition_residual,
        "trimer_over_protomer": result.ratio,
        "arms_out_of_plane_deg": result.mean_tilt_deg,
        "n_helices_used": result.n_helices,
        "helices_used": result.meta["helices_used"],
        "supports_paper": result.supports_paper,
        "blade": {"share_of_arrangement": split.blade_share,
                  "proximal_only_rmsd_A": split.proximal.arrangement_rmsd,
                  "n_distal_helices": len(split.distal_helices),
                  "summary": split.summary()},
        "summary": result.summary(),
    }
    try:
        beam = beam_angle(structure)
        out["beam_angle_deg"] = beam.mean_deg
        out["beam_out_of_plane_deg"] = beam.out_of_plane_deg
        out["beam_per_chain_deg"] = beam.angle_deg
    except ValueError as exc:
        out["beam_angle_deg"] = None
        out["beam_error"] = str(exc)
    out["published"] = {
        "beam_angle_deg": 60.0,
        "arms_out_of_plane_deg": 30.0,
        "note": "'The angle between the central pore axis and the beam in each "
                "subunit is about 60 degrees instead of 90'; the arms 'project "
                "approximately 30 degrees out of the plane defined by the pore'",
    }
    return out


def panel_4b(structure, reference: str = "mouse", **kw) -> dict:
    """The modelled micelle envelope, and the curvature it measures."""
    from ..structure.micelle import build_micelle

    envelope = build_micelle(structure, reference)
    out = {
        "n_belt_atoms": envelope.n_belt_atoms,
        "offset_A": envelope.offset,
        "n_vertices": envelope.n_vertices,
        "area_nm2": envelope.area() / 100.0,
        "volume_nm3": envelope.enclosed_volume() / 1000.0,
        "is_observed": envelope.is_observed,
        "definition": envelope.meta["definition"],
        "caveat": envelope.caveat,
        "summary": envelope.summary(),
    }
    if envelope.sphere is not None:
        out["belt_curvature_nm"] = envelope.sphere.radius / 10.0
        out["belt_fit_rmse_A"] = envelope.sphere.rmse
    out["published"] = {
        "idealised_radius_nm": 10.2,
        "contour": "6 sigma, unsharpened map",
        "note": "Figure 4b is the observed micelle density. This is a "
                "construction from the protein's own coordinates: the shell "
                "thickness is a registered parameter and only the curvature "
                "is a measurement."}
    return out


def panel_4c(structure, reference: str = "mouse", **kw) -> dict:
    """Surface electrostatic potential at 150 mM NaCl."""
    from ..physics.electrostatics import surface_potential

    result = surface_potential(structure)
    phi = result.potential
    return {
        "n_surface_points": int(len(phi)),
        "n_charges": result.meta["n_charges"],
        "net_charge_e": result.meta["net_charge"],
        "debye_length_A": result.debye_length,
        "bjerrum_length_A": result.bjerrum_length,
        "mean_kT_per_e": float(phi.mean()),
        "p5_kT_per_e": float(np.percentile(phi, 5)),
        "p95_kT_per_e": float(np.percentile(phi, 95)),
        "colour_scale_kT_per_e": result.scale,
        "fraction_saturated": result.fraction_saturated(),
        "positive_fraction": float(np.mean(phi > 0)),
        "method": result.meta["method"],
        "caveat": result.meta["not_apbs"],
        "published": {"scale_kT_per_e": 5.0, "ionic_strength_mM": 150,
                      "method": "APBS"},
    }


def panel_4_supplement(structure, reference: str = "mouse", **kw) -> dict:
    """The cap-to-loop charged interface, and whether it attracts."""
    from ..physics.electrostatic_patches import (patch_interaction,
                                                  patch_potential)
    from .interactions import detect_interactions

    acidic, basic = set(CED_ACIDIC_PATCH), set(LOOP_BASIC_PATCH)
    contacts = detect_interactions(structure, min_sequence_separation=3)
    across = [i for i in contacts.interactions
              if ({i.res_i, i.res_j} & acidic) and ({i.res_i, i.res_j} & basic)]

    interaction = patch_interaction(structure, CED_ACIDIC_PATCH, LOOP_BASIC_PATCH)
    return {
        "acidic_patch": sorted(acidic), "basic_patch": sorted(basic),
        "n_contacts": len(across),
        "contacts": [{"kind": i.kind, "distance_A": i.distance,
                      "a": f"{i.chain_i}/{i.name_i}{i.res_i}",
                      "b": f"{i.chain_j}/{i.name_j}{i.res_j}",
                      "domain_swapped": i.chain_i != i.chain_j}
                     for i in sorted(across, key=lambda x: x.distance)],
        "all_domain_swapped": all(i.chain_i != i.chain_j for i in across)
                              if across else False,
        "interaction_kT": interaction["energy_kT"],
        "cross_chain_kT": interaction["cross_chain_kT"],
        "same_chain_kT": interaction["same_chain_kT"],
        "attractive": interaction["attractive"],
        "acidic_surface_kT_per_e":
            patch_potential(structure, CED_ACIDIC_PATCH)["mean_potential_kT_per_e"],
        "basic_surface_kT_per_e":
            patch_potential(structure, LOOP_BASIC_PATCH)["mean_potential_kT_per_e"],
        "published": {
            "pairs": ["E2257-R1762", "D2264-R1761"],
            "note": "'Hydrogen bonds and salt bridges connect E2257 to R1762 "
                    "and D2264 in R1761 in a domain-swapped manner.'"},
    }
