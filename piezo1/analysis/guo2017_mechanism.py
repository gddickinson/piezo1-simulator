"""Figures 6 and 7 — the pore, and the membrane dome mechanism.

Split from :mod:`piezo1.analysis.guo2017_panels` at the project's length limit
and along the paper's own seam. Everything here is the argument Guo & MacKinnon
2017 is *for*: a narrow pore that cannot account for mechanosensitivity by
in-plane expansion, and a dome whose flattening can.

It is also where this project's own numbers come from. The 10.2 nm radius, the
120 nm² and the two-state Boltzmann in ``physics/dome.py`` are all Figure 7.
"""

from __future__ import annotations

import numpy as np

from ..parameters import PARAMETERS as _P

__all__ = ["panel_6ab", "panel_6_supplement", "panel_7a", "panel_7b",
           "panel_7c", "panel_7d", "panel_7_supplement", "PORE_CONSTRICTIONS"]


#: The pore radii Guo & MacKinnon state in the Results, in **mouse**
#: numbering: "At the level of the membrane inner leaflet, the pore radius at
#: positions E2537, P2536 and M2493 is 0.1 A, 0.4 A and 0.3 A, respectively."
#: E2487 and H2490 are labelled in Figure 6a,b without a stated radius.
PORE_CONSTRICTIONS = {2537: 0.1, 2536: 0.4, 2493: 0.3}


# --------------------------------------------------------------------------
# Figure 6 — the pore
# --------------------------------------------------------------------------

def panel_6ab(structure, reference: str = "mouse", **kw) -> dict:
    """The pore-radius profile, against the three radii the paper states."""
    from ..structure.pore import pore_profile
    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis

    blocks, _ = protomer_blocks(structure)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}
    axis = detect_c3_axis(blocks)
    profile = pore_profile(structure, axis)

    # Which end is extracellular sets the sign of the paper's "displacement
    # along pore from the top", so it is measured rather than assumed.
    from ..physics.pore_charge import cytosolic_end
    cytosolic = cytosolic_end(structure, axis)
    top_z = profile.z[-1] if cytosolic == 0 else profile.z[0]
    sign = -1.0 if cytosolic == 0 else 1.0

    rows = []
    for residue, published in sorted(PORE_CONSTRICTIONS.items()):
        touching = [s for s in profile.slices if residue in s.lining]
        if not touching:
            rows.append({"residue": residue, "published_radius_A": published,
                         "measured_radius_A": None,
                         "note": "residue lines no slice of the profile"})
            continue
        best = min(touching, key=lambda s: s.radius)
        rows.append({"residue": residue, "published_radius_A": published,
                     "measured_radius_A": best.radius,
                     "difference_A": best.radius - published,
                     "displacement_from_top_A": float(sign * (best.z - top_z))})

    measured = [r for r in rows if r["measured_radius_A"] is not None]
    # The absolute displacement is not comparable with Figure 6b's axis: the
    # paper measures from the top of the cap, while our profile begins wherever
    # the probe stops being bounded by protein, which on 6B3R is some 60 A
    # higher. The *spacing* between the constrictions is comparable, so it is
    # reported separately and the absolute number is labelled as ours.
    if measured:
        origin = min(r["displacement_from_top_A"] for r in measured)
        for row in measured:
            row["displacement_relative_A"] = \
                row["displacement_from_top_A"] - origin
    return {
        "bottleneck_radius_A": profile.bottleneck_radius,
        "bottleneck_z_A": profile.bottleneck_z,
        "bottleneck_lining": list(profile.bottleneck_lining()),
        "conductive": profile.is_conductive(),
        "cytosolic_end_index": cytosolic,
        "constrictions": rows,
        "mean_offset_A": (float(np.mean([r["difference_A"] for r in measured]))
                          if measured else None),
        "n_compared": len(measured),
        "published": {
            "radii_A": PORE_CONSTRICTIONS,
            "opening_radius_needed_A": 4.0,
            # Read off Figure 6b's axis, where M2493, P2536 and E2537 sit at
            # roughly 47, 55 and 58 A. Spacings, not positions, because only
            # the spacings survive the difference in axis origin.
            "spacing_A": {"M2493_to_P2536": 8.0, "P2536_to_E2537": 3.0},
            "axis_origin": "top of the cap; ours begins where the probe stops "
                           "being bounded, about 60 A higher, so only the "
                           "spacing between constrictions is comparable",
            "note": "'the pore radius at positions E2537, P2536 and M2493 is "
                    "0.1 A, 0.4 A and 0.3 A'; TEA permeation implies the open "
                    "pore reaches at least 4 A. Guo & MacKinnon used HOLE; "
                    "this project's profiler is an independent "
                    "implementation, so a systematic offset between them is "
                    "expected and is reported rather than absorbed."},
    }


def panel_6_supplement(structure, reference: str = "mouse", **kw) -> dict:
    """The cuff: elbow, base, hairpin and the pore-extension helix."""
    from ..core.annotations import load_annotations

    annotations = load_annotations()
    key = "mouse" if reference.startswith("mouse") else "human"
    wanted = ("elbow", "base", "hairpin", "pore_extension", "inner_helix",
              "outer_helix")
    ca = structure.mask_ca()

    from ..structure.protomers import protomer_blocks
    from ..structure.superpose import detect_c3_axis
    blocks, _ = protomer_blocks(structure)
    axis = detect_c3_axis(blocks) if blocks else None

    rows = []
    for domain in annotations.domains:
        if domain.id not in wanted:
            continue
        lo, hi = ((domain.mouse_start, domain.mouse_end) if key == "mouse"
                  else (domain.start, domain.end))
        mask = ca & (structure.res_seq >= lo) & (structure.res_seq <= hi)
        row = {"id": domain.id, "name": domain.name,
               "range": [lo, hi], "numbering": key,
               "n_ca_resolved": int(mask.sum()),
               "confidence": domain.confidence, "source": domain.source}
        if axis is not None and mask.any():
            radial = axis.radial(structure.xyz[mask])
            row["radial_min_A"] = float(radial.min())
            row["radial_mean_A"] = float(radial.mean())
            row["radial_max_A"] = float(radial.max())
        rows.append(row)
    return {
        "elements": rows,
        "published": {
            "elbow": [2116, 2142], "base": [2149, 2175],
            "hairpin": [2501, 2534], "numbering": "mouse",
            "note": "The PE helix is named and drawn but given no residue "
                    "range anywhere in the paper; ours is derived and marked "
                    "medium confidence."},
    }


# --------------------------------------------------------------------------
# Figure 7 — the dome model
# --------------------------------------------------------------------------

def panel_7a(structure=None, reference: str = "mouse", **kw) -> dict:
    """The idealised dome, and our measurement of a real one beside it."""
    from ..physics.dome_idealised import compare_with_measured, guo2017_dome

    dome = guo2017_dome()
    out = {"idealised": dome.as_dict(),
           "leaflet_radii_nm": list(dome.leaflet_radii()),
           "summary": dome.summary()}
    if structure is not None:
        from ..structure.geometry import measure_dome, tm_surface_points
        from ..structure.protomers import protomer_blocks
        blocks, _ = protomer_blocks(structure)
        if blocks:
            points, _ = tm_surface_points(structure, reference)
            out["comparison"] = compare_with_measured(
                measure_dome(blocks, points))
    out["published"] = {
        "radius_nm": 10.2, "center_height_nm": 4.0, "thickness_nm": 3.6,
        "note": "'a semi-sphere-shaped membrane 3.6 nm thick ... radius of "
                "10.2 nm ... centered 4.0 nm above the projection plane'"}
    return out


def panel_7b(structure, reference: str = "mouse", **kw) -> dict:
    """Beam in red, cross-helices in yellow, everything else grey."""
    from ..structure.architecture import cross_helices, cross_helix_scan

    found = cross_helices(structure, reference)
    scan = cross_helix_scan(structure, reference)
    per_chain: dict[str, list] = {}
    for segment in found:
        per_chain.setdefault(segment.chain, []).append(
            {"start": segment.start, "end": segment.end,
             "n_residues": segment.n_residues, "length_A": segment.length,
             "rise_A": segment.rise, "radius_A": segment.radius})
    return {
        "n_cross_helices": len(found),
        "per_protomer": {chain: rows for chain, rows in per_chain.items()},
        "n_per_protomer": (len(found) // max(len(per_chain), 1)),
        "threshold_deg": scan.default_threshold,
        "linker_tilt_median_deg": (float(np.median(scan.linker_tilts))
                                   if scan.linker_tilts else None),
        "transmembrane_tilt_median_deg":
            (float(np.median(scan.transmembrane_tilts))
             if scan.transmembrane_tilts else None),
        "populations_separated": scan.separated,
        "count_vs_threshold": list(zip(scan.thresholds, scan.counts)),
        "beam_excluded": True,
        "published": {
            "beam": [1300, 1365], "numbering": "mouse",
            "note": "'linkers between 4-TM units contain at least one helix "
                    "that runs perpendicular to the TM helices'. No residue "
                    "range for the cross-helices is given anywhere in the "
                    "paper, so this set is measured by that property and is "
                    "this project's, not the authors'."},
    }


def panel_7c(structure=None, reference: str = "mouse", n: int = 40,
             **kw) -> dict:
    """Projected area against flattening — the schematic, computed."""
    from ..physics.dome_idealised import flattening_series, guo2017_dome

    dome = guo2017_dome()
    series = flattening_series(dome, n=n)
    tension = _P.value("dome.lytic_tension") / 10.0
    return {
        "constraint": "mid-plane membrane area held constant while the cap "
                      "flattens; the one-parameter family in the polar "
                      "half-angle",
        "closed_polar_angle_deg": series[0].polar_angle_deg,
        "series": [{"polar_angle_deg": p.polar_angle_deg,
                    "radius_nm": p.radius_nm,
                    "projected_area_nm2": p.projected_area_nm2,
                    "delta_projected_nm2": p.delta_projected_nm2,
                    "bending_energy_kT": p.bending_energy_kT,
                    "delta_bending_kT": p.delta_bending_kT}
                   for p in series],
        "fully_flat": {
            "delta_projected_nm2": series[-1].delta_projected_nm2,
            "delta_bending_kT": series[-1].delta_bending_kT,
            "tension_work_kT": tension * series[-1].delta_projected_nm2},
        "tension_kT_per_nm2": tension,
        "note": ("Flattening completely releases the whole bending energy as "
                 "well as the whole projected-area gain. Taken with the "
                 "paper's own (dG_prot + dG_bend) of 20-40 k_BT, that puts "
                 "dG_prot at roughly +170 to +190 k_BT — a consequence of the "
                 "paper's numbers it does not spell out, and a large intrinsic "
                 "cost for the protein to pay."),
    }


def panel_7d(structure=None, reference: str = "mouse", **kw) -> dict:
    """The four theoretical activation curves."""
    from ..physics.dome import open_probability

    tension = np.linspace(0.0, 2.5, 251)
    curves = []
    for delta_g in (20.0, 40.0):
        for delta_area in (20.0, 60.0):
            probability = open_probability(tension, delta_area, delta_g)
            t50 = delta_g / delta_area
            curves.append({
                "delta_g_kT": delta_g, "delta_area_nm2": delta_area,
                "t50_kT_per_nm2": t50,
                "colour": "red" if delta_g == 20.0 else "blue",
                "p_open_at_t50": float(np.interp(t50, tension, probability)),
                "tension": tension.tolist(),
                "p_open": probability.tolist()})
    return {
        "equation": "P_o = (1 + exp[(dG_prot + dG_bend) - gamma * dA_proj])^-1",
        "curves": curves,
        "t50_values": {f"dG{int(c['delta_g_kT'])}_dA{int(c['delta_area_nm2'])}":
                       c["t50_kT_per_nm2"] for c in curves},
        "published": {
            "delta_g_kT": [20.0, 40.0], "delta_area_nm2": [20.0, 60.0],
            "note": "Figure 7d. The curves are fully determined by the ratio "
                    "dG/dA, so reproducing them is a check that this "
                    "project's two-state model is the paper's."},
    }


def panel_7_supplement(structure=None, reference: str = "mouse", **kw) -> dict:
    """Figure 7-figure supplement 1: the area and energy arithmetic."""
    from ..physics.dome_idealised import PUBLISHED_FIGURE7, guo2017_dome

    dome = guo2017_dome()
    values = dome.as_dict()
    values["stabilisation_kT"] = dome.stabilisation(
        _P.value("dome.lytic_tension") / 10.0)
    values["lytic_tension_kT_per_nm2"] = _P.value("dome.lytic_tension")

    rows = []
    for key, (published, tolerance, source) in PUBLISHED_FIGURE7.items():
        if key not in values:
            continue
        got = values[key]
        rows.append({"quantity": key, "published": published, "computed": got,
                     "difference": got - published, "tolerance": tolerance,
                     "agrees": abs(got - published) <= tolerance,
                     "source": source})
    return {
        "checks": rows,
        "n_agreeing": sum(1 for r in rows if r["agrees"]),
        "n_checked": len(rows),
        "note": ("Every number in Figure 7 and its supplement follows from two "
                 "lengths — a 10.2 nm sphere centred 4.0 nm above the plane — "
                 "by closed-form spherical-cap geometry. The agreement is "
                 "therefore a check on the arithmetic, not a measurement of "
                 "PIEZO1."),
    }
