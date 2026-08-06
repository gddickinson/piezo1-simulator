"""Report entries for the tag and permeation work (Rounds 31-33).

Split out of ``report.py`` to keep it under the project's 500-line limit. These
three share a seam: each depends on modules the older analyses do not, and each
carries a caveat about what is modelled rather than measured.
"""

from __future__ import annotations

import numpy as np

from ..core.structure import Structure
from ..parameters import PARAMETERS
from .report import _protomer_blocks

__all__ = ["analysis_fusion", "analysis_labelling", "analysis_permeation"]


def analysis_fusion(st: Structure, species: str, **kw) -> dict:
    """Where a C-terminal HaloTag would sit — a model, and labelled as one.

    Every distance here depends on ``fusion.linker_residues``, which no source
    for the construct states, so the linker is reported alongside the result
    rather than left implicit.
    """
    from ..structure.frame import apply_frame, canonical_transform
    from ..structure.fusion import build_fusion, load_halotag

    try:
        tag = load_halotag()
    except FileNotFoundError as exc:
        return {"error": str(exc)}

    # The fusion geometry is quoted relative to the conduction axis, so the
    # structure has to be in the frame that defines it.
    framed = apply_frame(st, canonical_transform(st))
    try:
        model = build_fusion(framed, tag)
    except (ValueError, RuntimeError) as exc:
        return {"error": str(exc)}

    reachable = model.volume.distances_from(model.pore_exit) / 10.0
    return {"tag_pdb": model.meta["tag_pdb"],
            "anchor_residues": list(model.anchor_residues),
            "linker_residues": model.meta["linker_residues"],
            "n_tags": model.n_tags,
            "c3_deviation_A": model.c3_deviation(),
            "tag_centre_to_pore_exit_nm": float(model.pore_exit_distances()[0]),
            "envelope_median_nm": float(np.median(reachable)),
            "envelope_range_nm": [float(reachable.min()), float(reachable.max())],
            "fraction_in_4_to_6_nm": float(((reachable >= 4.0)
                                            & (reachable <= 6.0)).mean()),
            "accessible_volume_nm3": model.meta["accessible_volume_nm3"],
            "occluded_fraction": model.volume.occluded_fraction,
            "min_clearance_A": model.meta["min_clearance"],
            "tag_radius_A": model.meta["tag_radius"],
            "clashes": model.meta["clashes"],
            "note": model.meta["note"]}


def analysis_labelling(st: Structure, species: str, **kw) -> dict:
    """HaloTag labelling of the three C-terminal sites.

    The kinetics are imported from ``halotag_binding_sim``; this reports them at
    the registered protocol and, where the tag geometry is available, on the
    modelled site positions.
    """
    from .labelling import (LabellingConditions, occupancy_distribution,
                            population_summary, site_labelled_fraction,
                            time_to_fraction)

    conditions = LabellingConditions()
    protocol_t = PARAMETERS.value("labelling.incubation_time")
    times = np.linspace(0.0, max(protocol_t, 3600.0), 601)
    result = population_summary(times, conditions)
    at_protocol = result.at(protocol_t)

    out = {"source": result.meta["source"],
           "conditions": conditions.summary(),
           "concentration_M": conditions.concentration,
           "incubation_time_s": protocol_t,
           "asymptote": conditions.asymptote,
           "p_site_at_protocol": at_protocol["p_site"],
           "fully_labelled_at_protocol": at_protocol["fully_labelled"],
           "detectable_at_protocol": at_protocol["detectable"],
           "mean_dyes_at_protocol": at_protocol["mean_dyes"],
           "dye_histogram_at_protocol": at_protocol["occupancy"],
           "time_to_99_percent_s": time_to_fraction(0.99, conditions)}

    # What an incomplete-reactivity population would look like instead. Reported
    # because the two routes to a dye mixture are easy to conflate and only one
    # of them is available at a saturating concentration.
    reduced = LabellingConditions(active_fraction=0.9)
    ceiling_p = float(site_labelled_fraction(6 * 3600.0, reduced))
    out["if_90_percent_reactive"] = {
        "asymptote": reduced.asymptote,
        "dye_histogram": occupancy_distribution(ceiling_p,
                                                reduced.n_sites).tolist()}

    try:
        from ..structure.frame import apply_frame, canonical_transform
        from ..structure.fusion import build_fusion, load_halotag
        from .labelling import label_sites
        framed = apply_frame(st, canonical_transform(st))
        sites = label_sites(build_fusion(framed, load_halotag()), t=protocol_t)
        out["sites"] = {"anchor_residues": sites["anchor_residues"],
                        "n_dyes_drawn": sites["n_dyes"],
                        "note": sites["note"]}
    except Exception as exc:
        out["sites"] = {"error": str(exc)}
    return out


def analysis_permeation(st: Structure, species: str, step: float = 1.0,
                        **kw) -> dict:
    """Ion current through the pore, gated by the wetting verdict."""
    from ..physics.permeation import (default_species, series_conductance,
                                      solve_pnp)
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    from .hydration import load_grid, predict_wetting

    blocks, _ = _protomer_blocks(st)
    if blocks is None:
        return {"error": "needs three well-resolved protomers"}

    profile = pore_profile(st, detect_c3_axis(blocks), step=step)
    grid = load_grid()
    wetting = predict_wetting(st, profile, grid) if grid.available else None

    result = solve_pnp(profile, wetting)
    series = series_conductance(np.asarray(profile.z) * 1e-10,
                                np.asarray(profile.radius) * 1e-10)
    out = {"bottleneck_radius_A": profile.bottleneck_radius,
           "conducting": result.is_conducting,
           "conductance_pS": result.conductance_pS,
           "published_pS": PARAMETERS.value("permeation.published_conductance"),
           "blocked_by": result.blocked_by,
           "mechanisms": result.meta.get("mechanisms", []),
           "n_mechanisms": len(result.meta.get("mechanisms", []))}
    if result.is_conducting:
        out.update({
            "current_pA": result.current * 1e12,
            "voltage_mV": result.voltage * 1e3,
            "pore_Gohm": result.pore_ohm / 1e9,
            "access_Gohm": result.access_ohm / 1e9,
            "independent_check_pS": series["conductance"] * 1e12,
            "converged": result.converged,
            "debye_length_A": result.meta["debye_length_A"],
            "double_layers_overlap": result.meta["double_layers_overlap"],
            "calcium_2mM_pS": solve_pnp(
                profile, wetting,
                species=default_species(calcium=0.002)).conductance_pS})
    out["note"] = ("continuum model of an atomic-scale pore; the in-pore "
                   "diffusivity and ion radius are unmeasured and the result "
                   "spans 16-94 pS across their plausible ranges")
    return out
