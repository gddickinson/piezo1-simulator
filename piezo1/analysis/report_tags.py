"""Report entries for the tag and permeation work (Rounds 31-33).

Split out of ``report.py`` to keep it under the project's 500-line limit. These
three share a seam: each depends on modules the older analyses do not, and each
carries a caveat about what is modelled rather than measured.
"""

from __future__ import annotations

import numpy as np

from ..parameters import PARAMETERS as _P

from ..core.structure import Structure
from ..parameters import PARAMETERS
from ..structure.protomers import protomer_blocks

__all__ = ["analysis_fusion", "analysis_labelling", "analysis_permeation",
           "analysis_nanodomain", "analysis_prediction_record",
           "analysis_ligands", "analysis_paired_variant",
           "analysis_fluctuations"]


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


def analysis_permeation(st: Structure, species: str, step: float | None = None,
                        **kw) -> dict:
    """Ion current through the pore, gated by the wetting verdict."""
    if step is None:
        step = _P.value("pore.step")
    from ..physics.permeation import (default_species, series_conductance,
                                      solve_pnp)
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    from .hydration import load_grid, predict_wetting

    blocks, _ = protomer_blocks(st)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}

    axis = detect_c3_axis(blocks)
    profile = pore_profile(st, axis, step=step)
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
        out["fixed_charge"] = _charged_pore(st, profile, axis, species, wetting)
    out["note"] = ("continuum model of an atomic-scale pore; the in-pore "
                   "diffusivity and ion radius are unmeasured and the result "
                   "spans 16-94 pS across their plausible ranges")
    return out


def _charged_pore(st: Structure, profile, axis, species: str, wetting) -> dict:
    """What the pore's own charge does to the current and to selectivity.

    Both routes are reported because they disagree, and the disagreement is
    informative: ``curated`` is the literature's answer to which residues line
    the pore, ``lining`` is the coordinates' answer, and neither is evidence
    for the other. The neutral pore is reported alongside as the baseline — it
    is *already* weakly cation-selective from size exclusion alone, so a
    selective answer is not by itself evidence that the charge did anything.
    """
    from ..physics.permeation import solve_pnp
    from ..physics.pore_charge import MODES, cytosolic_end, pore_charge
    from ..physics.selectivity import measure_selectivity

    end = cytosolic_end(st, axis)
    out: dict = {"cytosolic_end_index": end,
                 "published_pcl_pna": PARAMETERS.value(
                     "permeation.published_pcl_pna")}
    try:
        neutral = measure_selectivity(profile, cytosolic_index=end,
                                      wetting=wetting)
        out["neutral"] = {"p_anion_over_cation": neutral.p_anion_over_cation,
                          "reversal_mV": neutral.reversal_mV,
                          "protocol": neutral.meta["protocol"]}
    except ValueError as exc:
        return {"error": str(exc)}

    for mode in MODES:
        charge = pore_charge(st, profile, axis, mode=mode, species=species)
        result = solve_pnp(profile, wetting, fixed_charge=charge.density)
        rows = charge.residue_summary()
        row = {"n_groups": charge.n_groups,
               "net_charge_e": charge.net_charge,
               # The geometric route names seventeen residues, which is a wall
               # of text in a result window; the curated one names two and the
               # detail is the point there.
               "residues": rows if mode == "curated" else
               [f"{r['res_name']}{r['res_seq']}x{r['copies']}"
                f"{'+' if r['charge'] > 0 else '-'}" for r in rows],
               "peak_density_M": charge.peak_density / 1000.0,
               "conductance_pS": result.conductance_pS,
               "peak_in_pore_M": result.meta.get("peak_in_pore_M"),
               "exceeds_packing_limit": result.meta.get(
                   "exceeds_packing_limit"),
               "electroneutrality_residual": result.meta.get(
                   "electroneutrality_residual"),
               "converged": result.converged}
        selectivity = measure_selectivity(profile, fixed_charge=charge.density,
                                          cytosolic_index=end, wetting=wetting)
        row["p_anion_over_cation"] = selectivity.p_anion_over_cation
        row["reversal_mV"] = selectivity.reversal_mV
        row["cation_selective"] = selectivity.cation_selective
        out[mode] = row

    out["note"] = ("the fixed charge makes the pore cation-selective on both "
                   "routes, which is the direction the measurement has; the "
                   "two routes bracket the published ratio rather than "
                   "reproducing it, and the curated route reaches an in-pore "
                   "concentration a continuum model cannot be believed at")
    return out


def analysis_nanodomain(st: Structure, species: str, step: float | None = None,
                        **kw) -> dict:
    """Calcium at the tag when this channel is open, and whether a sensor saturates.

    Needs the two numbers earlier rounds produce: the unitary current from the
    pore, and the tag distance from the fusion model. A closed structure carries
    no current, so the nanodomain is reported for the **open** reference entry
    with the tag geometry taken from the loaded structure — and says so.
    """
    if step is None:
        step = _P.value("pore.step")
    from ..physics.nanodomain import Nanodomain, screening_length, sweep
    from ..physics.permeation import solve_pnp
    from ..structure.frame import apply_frame, canonical_transform
    from ..structure.fusion import build_fusion, load_halotag
    from ..structure.pore import pore_profile
    from ..structure.superpose import detect_c3_axis
    from .hydration import load_grid, predict_wetting

    grid = load_grid()
    framed = apply_frame(st, canonical_transform(st))
    # The axis must come from the SAME frame as the coordinates it is used with.
    # Detecting it on the unframed structure and applying it to the framed one
    # measured the pore along a line that misses the pore, and reported the
    # closed 8YEZ as carrying 32 pA.
    blocks, _ = protomer_blocks(framed)
    if not blocks:
        return {"error": "needs three well-resolved protomers"}
    profile = pore_profile(framed, detect_c3_axis(blocks), step=step)
    wetting = predict_wetting(framed, profile, grid) if grid.available else None
    result = solve_pnp(profile, wetting)

    current, source = abs(result.current), st.name
    if current == 0.0:
        # Every deposited human structure is closed (Round 34), so a nanodomain
        # for the loaded entry would be exactly zero and uninformative. Borrow
        # the current from the one open-like entry and label it.
        from ..io.registry import load_registry
        record = load_registry().get("11ZC")
        if record is None or not record.available:
            return {"error": "this structure is closed and 11ZC is not "
                             "downloaded to supply an open-state current"}
        open_st = Structure.from_file(record.path)
        open_framed = apply_frame(open_st, canonical_transform(open_st))
        open_blocks, _ = protomer_blocks(open_framed)
        open_profile = pore_profile(open_framed, detect_c3_axis(open_blocks),
                                    step=step)
        open_wetting = (predict_wetting(open_framed, open_profile, grid)
                        if grid.available else None)
        current = abs(solve_pnp(open_profile, open_wetting).current)
        source = "11ZC (open-like)"

    try:
        model = build_fusion(framed, load_halotag())
    except (FileNotFoundError, ValueError, RuntimeError) as exc:
        return {"error": f"tag geometry unavailable: {exc}"}

    reach = model.volume.distances_from(model.pore_exit) / 10.0
    fraction = PARAMETERS.value("nanodomain.calcium_current_fraction")
    nano = Nanodomain(current_A=current, calcium_fraction=fraction,
                      distance_m=float(model.pore_exit_distances()[0]) * 1e-9,
                      envelope_m=(float(reach.min()) * 1e-9,
                                  float(reach.max()) * 1e-9))
    far_c, near_c, far_occ, near_occ = nano.envelope_range()
    rows = sweep(current)
    return {
        "current_source": source,
        "unitary_current_pA": current * 1e12,
        "calcium_fraction": fraction,
        "tag_distance_nm": nano.distance_m * 1e9,
        "screening_length_nm": screening_length() * 1e9,
        "calcium_uM": nano.concentration_M * 1e6,
        "sensor_occupancy": nano.occupancy,
        "resting_occupancy": nano.resting_occupancy,
        "saturated": nano.saturated,
        "across_envelope": {
            "calcium_uM": [far_c * 1e6, near_c * 1e6],
            "occupancy": [far_occ, near_occ]},
        "falsifiers": nano.falsifiers(),
        "sweep_combinations": len(rows),
        "sweep_not_saturated": sum(1 for r in rows if r["occupancy"] < 0.9),
        "note": ("The sensor is saturated whenever its own channel opens, so "
                 "puncta brightness reports labelling stoichiometry and open "
                 "probability rather than calcium amplitude. Tag distance is "
                 "modelled, not measured."),
    }


def analysis_prediction_record(st: Structure, species: str, **kw) -> dict:
    """What a variant prediction from this project is entitled to claim.

    The project's central claim, its three pre-registered tests, and the
    standing caveats — surfaced rather than left in `docs/`. Independent of the
    loaded structure, because the record is about the predictor, not about any
    one entry.
    """
    from .prediction_record import (VALIDATION_RECORD, evidence_levels,
                                    headline, verify_record, what_it_means)

    return {
        "central_claim": ("predicting gain- versus loss-of-function from "
                          "structure"),
        "verdict": headline(),
        "tests": [{"round": e.round, "predictor": e.predictor,
                   "n_gof": e.n_gof, "n_lof": e.n_lof,
                   "cliffs_delta": e.cliffs_delta, "p_value": e.p_value,
                   "power_at_large_effect": e.power_at_large,
                   "conclusion": e.conclusion, "document": e.document}
                  for e in VALIDATION_RECORD],
        "what_this_means": what_it_means(),
        "evidence_levels": evidence_levels(),
        "record_matches_the_stored_run": verify_record(),
    }


def analysis_ligands(st: Structure, species: str, **kw) -> dict:
    """The curated modulators, and how much is known about where they bind.

    Reports the site evidence level with every site, because no PIEZO structure
    with a bound small-molecule modulator has been deposited and a pocket drawn
    on a structure looks exactly like one observed in it.
    """
    from ..core.ligands import load_ligands

    ligands = load_ligands()
    if not len(ligands):
        return {"error": ligands.note}
    return {
        "summary": ligands.summary(),
        "no_bound_structure_exists": not ligands.any_observed_site,
        "site_evidence_levels": ligands.evidence_levels,
        "ligands": [
            {"name": l.name, "role": l.role, "kind": l.kind,
             "formula": l.formula, "molecular_weight": l.molecular_weight,
             "pubchem_cid": l.pubchem_cid, "uniprot": l.uniprot,
             "potency": l.potency_text(), "site": l.site_text(),
             "site_evidence": l.site_evidence,
             "description": l.description}
            for l in ligands.ligands],
        "note": ligands.note,
    }


def analysis_paired_variant(st: Structure, species: str, **kw) -> dict:
    """The one variant-versus-wild-type structural comparison available.

    Independent of the loaded structure: the comparison is between specific
    deposited entries, not whatever happens to be on screen.
    """
    from .paired_variant import compare

    result = compare()
    if result is None or len(result.wild_type) < 2:
        return {"error": "structures or the CHAP grid are not downloaded"}
    return {
        "variant": result.variant.pdb,
        "wild_type_set": [m.pdb for m in result.wild_type],
        "excluded_as_duplicates": [
            {"entry": a, "identical_to": b}
            for a, b in result.excluded_duplicates],
        "metrics": result.report(),
        "distinguishable": result.distinguishable,
        "summary": result.summary(),
        "note": ("n = 1. This says what the deposited structures show, not what "
                 "R2456H does. Every entry is closed, and R2456H's phenotype is "
                 "slowed inactivation, which a closed structure need not show."),
    }


def analysis_hybrid(st: Structure, species: str, **kw) -> dict:
    """The full-length model: experimental core plus the predicted distal blade.

    Reported as two populations that must never be averaged. Cryo-EM resolves
    roughly 570-2521; the rest is AlphaFold, and the numbers that say how much
    to trust it — the fraction of the graft clearing pLDDT 70, and the 75 A by
    which the two models disagree away from the seam — are returned beside the
    atom counts rather than left for a caller to look up.
    """
    from ..structure.hybrid import build_hybrid_model

    try:
        model = build_hybrid_model(st)
    except (FileNotFoundError, ValueError) as exc:
        return {"error": str(exc)}

    predicted = model.predicted
    n_pred = int(predicted.sum())
    confident = int(model.confident_prediction.sum())
    residues = model.res_seq[predicted]
    return {
        "seam_residue": model.seam_residue,
        "experimental_atoms": int(model.experimental_only.sum()),
        "predicted_atoms": n_pred,
        "predicted_residues": int(np.unique(residues).size) if n_pred else 0,
        "predicted_range": ([int(residues.min()), int(residues.max())]
                            if n_pred else None),
        "confident_fraction": confident / n_pred if n_pred else float("nan"),
        "mean_plddt": (float(np.nanmean(model.plddt[predicted]))
                       if n_pred else float("nan")),
        "seam_rmsd_A": model.overlap_rmsd,
        "global_rmsd_A": model.global_rmsd,
        "warnings": model.warnings(),
        "summary": model.summary(),
        "note": ("Two populations, not one structure. Every atom carries its "
                 "source so no analysis can average across the join; the "
                 "grafted blade is a PREDICTION and the seam-local fit says "
                 "nothing about the rest of it."),
    }


def analysis_fluctuations(st: Structure, species: str,
                          n_modes: int | None = None, **kw) -> dict:
    """The elastic network against this entry's own B-factors.

    Reported with the control beside it, because a residue with many
    neighbours moves less whether or not there is a normal mode anywhere near
    it. An entry whose column cannot be interpreted returns the reason instead
    of a number.
    """
    from .fluctuations import assess_b_factors, compare_fluctuations

    quality = assess_b_factors(st)
    out = {"n_residues_in_column": quality.n_residues,
           "n_distinct_values": quality.n_distinct,
           "distinct_fraction": quality.distinct_fraction,
           "b_range": [quality.minimum, quality.maximum],
           "column_is_plddt": quality.is_confidence,
           "usable": quality.usable, "reason": quality.reason}
    if not quality.usable:
        out["note"] = ("no comparison was made; the reason above is the "
                       "result, not a failure to compute")
        return out

    result = compare_fluctuations(st, n_modes=n_modes)
    if not result.available:
        out.update({"usable": False, "reason": result.quality.reason})
        return out
    out.update({
        "n_residues_compared": len(result.residues),
        "n_modes": result.n_modes,
        "pearson": result.pearson_r,
        "spearman": result.spearman_r,
        "control_contact_number_pearson": result.control_pearson,
        "control_contact_number_spearman": result.control_spearman,
        "beats_control": result.beats_control,
        "control_inverted": result.control_inverted,
        "by_mode_count": result.by_mode_count,
        "note": ("Spearman is the number to read: the relationship is "
                 "monotone but not linear, so Pearson is dominated by a few "
                 "very mobile residues. A negative control means this entry's "
                 "B-factor rises with burial, which no mobility does — there "
                 "the column is the problem, not the network.")})
    return out
