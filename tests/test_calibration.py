"""Every checking instrument, calibrated against a case with a known answer.

This project's most expensive mistakes have all been the same shape: an
*alternative* route, built to check the main one, was itself wrong and produced
a plausible number rather than an error. A spheroid fitter that would have
reported 89% model error. A document checker that could not read its own
documents' minus sign. A probe whose "no effect" came from badly chosen
coordinates. In each case the checker disagreed with the pipeline and the
disagreement looked like a finding.

So a checking instrument is a measuring instrument, and an uncalibrated one is
worse than none — it manufactures findings. This file holds the register of
instruments and the calibration each one is subject to, and a guard that no
instrument can be added without one.

The register was assembled by auditing the modules rather than by trusting
their docstrings. Four instruments turned out to have no known-answer case and
are calibrated here for the first time.
"""

from __future__ import annotations

import importlib
import itertools
import re
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent

#: Modules whose public callables exist to *check* something. Every lowercase
#: name in their ``__all__`` must appear in :data:`CALIBRATED` below.
CHECKING_MODULES = (
    "analysis.crosscheck",
    "analysis.hydropathy",
    "structure.planarity",
    "structure.architecture",
    "structure.micelle",
    "analysis.crosscheck_methods",
    "analysis.model_error",
    "analysis.uncertainty",
    "analysis.validation",
    "analysis.design",
    "analysis.provenance_chain",
    "analysis.parameter_effect",
    "analysis.fluctuations",
    # Round 89. Both are instruments in the strict sense: neither measures a
    # property of PIEZO1, both exist to decide whether another number may be
    # believed. ``alignment_windows`` earns its place twice over — it was
    # wrong twice before it was calibrated, once in its statistic and once in
    # its width.
    "analysis.homology",
    "analysis.alignment_windows",
    "analysis.homology_structure",
    "structure.assembly",
    "structure.clashes",
    "dead_code",
)

#: instrument -> the test that drives it against a known answer.
CALIBRATED = {
    # crosscheck: re-derivations of physics results
    "crosscheck.dome_curvature_by_cap_geometry":
        "test_crosscheck.py::test_the_exact_cap_route_is_also_exact_at_every_slope",
    "crosscheck.dome_curvature_by_parabola":
        "test_crosscheck.py::test_the_parabola_route_degrades_with_slope",
    "crosscheck.gating_overlap_by_distances":
        "test_crosscheck.py::test_distance_route_recovers_a_planted_mode",
    "crosscheck.t50_by_ode_integration":
        "test_crosscheck.py::test_ode_integration_reproduces_the_matrix_exponential",
    "crosscheck.t50_by_steady_state":
        "test_crosscheck.py::test_the_steady_state_is_a_different_quantity",
    "crosscheck.compare": "test_crosscheck.py",

    # crosscheck_methods: re-derivations of algorithms
    "crosscheck_methods.pore_radius_by_random_search":
        "test_crosscheck_methods.py (slab with a known gap)",
    "crosscheck_methods.sasa_by_monte_carlo":
        "test_crosscheck_methods.py (isolated sphere, 4 pi r^2)",
    "crosscheck_methods.conservation_by_kmer_anchoring":
        "test_crosscheck_methods.py (random sequences isolate the bias)",
    "crosscheck_methods.pc1_by_power_iteration":
        "test_crosscheck_methods.py (planted principal component)",
    "crosscheck_methods.check_pore_radius": "test_crosscheck_methods.py",
    "crosscheck_methods.check_sasa": "test_crosscheck_methods.py",
    "crosscheck_methods.check_conservation": "test_crosscheck_methods.py",
    "crosscheck_methods.check_pca": "test_crosscheck_methods.py",

    # model_error
    "model_error.fit_spheroid":
        "test_model_error.py (known semi-axes, incl. the original failure)",
    "model_error.spring_model_error":
        "test_calibration.py::test_spring_model_error_recovers_a_planted_mode",
    "model_error.pore_convention_error": "test_model_error.py (the null's mechanism)",
    "model_error.dome_model_error": "test_model_error.py",
    "model_error.compare_with_sampling": "test_model_error.py",

    # uncertainty
    "uncertainty.bootstrap":
        "test_uncertainty.py::test_bootstrap_recovers_a_known_interval",
    "uncertainty.sensitivity":
        "test_calibration.py::test_sensitivity_range_is_exactly_the_known_extremes",
    "uncertainty.parameter_range":
        "test_uncertainty.py::test_parameter_range_restores_even_when_the_statistic_raises",
    "uncertainty.format_with_interval": "test_uncertainty.py",

    # validation
    "validation.permutation_test":
        "test_calibration.py::test_permutation_p_matches_exhaustive_enumeration",
    "validation.cliffs_delta": "test_validation.py::test_cliffs_delta_extremes",
    "validation.bootstrap_cliffs_delta": "test_validation.py",
    "validation.auroc": "test_validation.py::test_auroc_known_cases",
    "validation.interpret_delta": "test_validation.py::test_delta_interpretation_thresholds",

    # homology: the family, and whether a percentage may be believed
    "homology.align_pair":
        "test_homology.py::test_the_null_is_calibrated_on_two_known_answers",
    "homology.shuffled_null":
        "test_homology.py::test_the_null_can_say_no_when_there_is_no_homology",
    "homology.relationship":
        "test_homology.py::test_at_least_one_pair_has_an_identity_indistinguishable_from_chance",
    "homology.family_matrix": "test_homology.py (36 pairs, each against its null)",
    "homology.family": "test_homology.py::test_the_family_is_nine_and_every_member_is_committed",
    "homology.member": "test_homology.py (family lookup)",
    "homology.group_of": "test_homology.py::test_no_two_members_share_a_length_or_a_helix_count_architecture",

    # assembly: is a trimer built from one protomer a model of anything?
    "assembly.assemble_trimer":
        "test_assembly.py::test_a_protomer_rebuilt_on_its_own_trimer_reproduces_it_exactly",
    "assembly.borrowed_fraction":
        "test_assembly.py::test_most_of_an_assembly_s_shape_is_the_template_s_and_it_says_so",
    "assembly.best_template":
        "test_assembly.py::test_the_template_is_chosen_by_protein_then_by_coverage",
    "clashes.count_clashes":
        "test_assembly.py::test_the_clash_counter_is_near_zero_on_real_trimers",
    "assembly.is_monomer":
        "test_assembly.py::test_a_real_trimer_is_refused_rather_than_replaced_by_a_model_of_itself",

    # homology_structure: is one entry pair a measurement of two proteins?
    "homology_structure.mode_overlap_spread":
        "test_homology_structure.py::test_the_spread_reports_instability_and_piezo2_is_the_positive_control",
    "homology_structure.compare_structures":
        "test_homology_structure.py::test_a_single_entry_pair_is_not_a_measurement_of_the_two_proteins",
    "homology_structure.index_pairing_valid":
        "test_homology_structure.py::test_index_pairing_is_refused_where_the_helix_counts_differ",
    "homology_structure.comparable_entries":
        "test_homology_structure.py::test_comparable_entries_are_best_resolved_first",
    "homology_structure.helix_counts":
        "test_homology.py::test_no_two_members_share_a_length_or_a_helix_count_architecture",

    # alignment_windows: is the alignment in register here?
    "alignment_windows.window_score":
        "test_homology.py::test_the_window_finds_a_planted_block_that_whole_sequence_identity_misses",
    "alignment_windows.window_identity":
        "test_homology.py::test_the_window_is_calibrated_on_a_sequence_against_itself",
    "alignment_windows.alignment_windows":
        "test_homology.py::test_the_window_refuses_a_composition_matched_shuffle",
    "alignment_windows.column_scores":
        "test_homology.py (self-alignment scores above its own null)",
    "alignment_windows.window_null_distribution":
        "test_homology.py::test_window_width_was_chosen_by_a_power_scan_not_by_taste",

    # design
    "design.power_curve":
        "test_design.py::test_null_effect_gives_the_nominal_false_positive_rate",
    "design.minimum_detectable_effect":
        "test_calibration.py::test_minimum_detectable_effect_delivers_its_target_power",
    "design.sample_size_for": "test_design.py",
    "design.shift_for_delta": "test_design.py::test_simulated_effect_matches_the_requested_effect",
    "design.delta_for_shift": "test_design.py",
    "design.benjamini_hochberg": "test_design.py (known p-value sets)",
    "design.leave_one_out": "test_design.py::test_leave_one_out_recovers_a_real_signal",

    # provenance_chain
    "provenance_chain.number_in_document":
        "test_provenance_chain.py::test_it_finds_a_number_written_in_any_reasonable_form",
    "provenance_chain.resolved_keys":
        "test_provenance_chain.py::test_a_registered_parameter_that_nothing_reads_is_detected",
    "provenance_chain.unwired_parameters": "test_provenance_chain.py",
    "provenance_chain.record_sources":
        "test_provenance_chain.py::test_it_records_the_parameters_a_computation_actually_reads",
    "provenance_chain.trace": "test_provenance_chain.py",
    "provenance_chain.walk": "test_provenance_chain.py",
    "provenance_chain.git_state": "test_provenance_chain.py",

    # dead_code
    "dead_code.audit":
        "test_dead_code.py::test_the_detector_is_calibrated",
    "dead_code.calibration":
        "test_dead_code.py::test_the_detector_is_calibrated",
    "dead_code.reference_counts":
        "test_dead_code.py::test_same_file_references_count",

    # fluctuations: the elastic network against the deposited B-factors.
    # A checking instrument in the strict sense — it is the standard test of
    # whether the network describes the molecule — so it is registered here
    # rather than treated as an analysis that happens to produce a number.
    "fluctuations.pearson":
        "test_fluctuations.py::test_the_correlations_return_their_analytic_limits",
    "fluctuations.spearman":
        "test_fluctuations.py::test_ranks_average_ties_rather_than_inventing_an_order",
    "fluctuations.predicted_msf":
        "test_fluctuations.py::"
        "test_the_comparison_recovers_a_planted_fluctuation_and_a_bad_network_does_not",
    "fluctuations.contact_number":
        "test_fluctuations.py::"
        "test_contact_number_is_a_control_and_not_a_copy_of_the_prediction",
    "fluctuations.assess_b_factors":
        "test_fluctuations.py::"
        "test_a_predicted_model_is_refused_and_the_gate_points_the_right_way",
    "fluctuations.observed_b_factors":
        "test_fluctuations.py::"
        "test_observed_b_factors_are_read_per_residue_not_per_atom",
    "fluctuations.compare_fluctuations":
        "test_fluctuations.py::test_a_shuffled_observation_correlates_with_nothing",
    "fluctuations.survey_fluctuations":
        "test_fluctuations.py::test_the_survey_of_every_downloaded_entry",

    # hydropathy: the 4-TM repeat test is a checking instrument — it is the
    # evidence for the nine-unit architecture the whole project is built on,
    # including the distal blade the full-length model grafts.
    "hydropathy.repeat_periodicity":
        "test_hydropathy.py::test_the_repeat_test_finds_a_planted_period "
        "and ::test_the_repeat_test_says_no_to_a_true_null",
    "hydropathy.hydropathy_profile":
        "test_hydropathy.py::test_a_uniform_sequence_gives_its_own_value_everywhere",
    "hydropathy.predict_segments":
        "test_hydropathy.py::"
        "test_a_planted_hydrophobic_block_is_found_where_it_was_planted",
    "hydropathy.compare_with_reference":
        "test_hydropathy.py::"
        "test_the_default_threshold_recovers_few_helices_and_that_is_reported",
    "hydropathy.annotated_hydropathy":
        "test_hydropathy.py::test_piezo1_helices_sit_below_the_conventional_membrane_cut",
    "hydropathy.threshold_scan":
        "test_hydropathy.py::test_the_threshold_scan_trades_recall_for_specificity_monotonically",
    "hydropathy.load_reference":
        "test_hydropathy.py::test_a_missing_reference_names_what_it_looked_for",

    # planarity: Figure 4a's claim as residuals. A checking instrument because
    # it decides whether a published statement about shape holds, and because
    # its first version reported a tautology (points replicated three times)
    # as a control.
    "planarity.fit_plane":
        "test_planarity.py::test_a_known_out_of_plane_displacement_is_recovered",
    "planarity.planarity":
        "test_planarity.py::"
        "test_the_flattened_structure_is_the_control_that_makes_this_mean_something",
    "planarity.blade_dependence":
        "test_planarity.py::test_coverage_decides_the_answer_and_the_module_says_so",
    "planarity.beam_angle":
        "test_planarity.py::"
        "test_the_beam_angle_opens_towards_90_when_the_channel_flattens",

    # architecture: the helix detector, calibrated on analytic helices of
    # known rise, radius and turn before it is used on coordinates.
    "architecture.helical_windows":
        "test_cross_helices.py::test_the_estimator_is_unbiased_on_an_ideal_helix",
    "architecture.helical_segments":
        "test_cross_helices.py::test_the_beam_is_found_as_one_helix_not_two",
    "architecture.cross_helices":
        "test_cross_helices.py::test_the_beam_is_not_reported_as_a_cross_helix",
    "architecture.cross_helix_scan":
        "test_cross_helices.py::"
        "test_the_threshold_is_reported_against_a_scan_not_asserted",
    "architecture.ideal_helix":
        "test_cross_helices.py::test_the_ideal_helix_generator_is_what_it_claims",

    # micelle: the modelled Figure 4b envelope. A measuring instrument rather
    # than a checking one, but it produces a curvature that gets compared with
    # a published number, so it is registered on the same terms as the rest.
    "micelle.distance_field":
        "test_ui_membrane_views.py::"
        "test_the_envelope_of_one_point_is_a_sphere_of_the_offset_radius",
    "micelle.build_micelle":
        "test_ui_membrane_views.py::"
        "test_the_micelle_encloses_the_belt_and_says_it_is_a_model",
    "micelle.belt_atoms":
        "test_ui_membrane_views.py::"
        "test_the_belt_is_apolar_transmembrane_side_chains",

    # parameter_effect
    "parameter_effect.measure_effect":
        "test_parameter_effect.py::test_a_parameter_with_no_effect_is_reported_as_such",
    "parameter_effect.probe_effects": "test_parameter_effect.py",
}


# ---------------------------------------------------- the register guard

def test_every_checking_instrument_has_a_calibration():
    """No instrument may be added without a known-answer case.

    This is the mechanism that keeps the rule from decaying into an
    aspiration, in the same way `parameter_audit` does for registered numbers.
    """
    missing = []
    for module in CHECKING_MODULES:
        mod = importlib.import_module(f"piezo1.{module}")
        short = module.split(".")[-1]
        for name in getattr(mod, "__all__", []):
            if name[:1].isupper() or name.startswith("_"):
                continue                      # dataclasses carry no answer
            if f"{short}.{name}" not in CALIBRATED:
                missing.append(f"{short}.{name}")
    assert not missing, (
        "these checking instruments have no calibration registered: "
        + ", ".join(sorted(missing)))


def test_named_calibrating_tests_exist():
    """A register pointing at a test that does not exist proves nothing."""
    missing = []
    for instrument, where in CALIBRATED.items():
        filename = where.split("::")[0].split(" ")[0]
        path = ROOT / filename
        if not path.exists():
            missing.append(f"{instrument} -> {filename} (no such file)")
            continue
        if "::" in where:
            testname = where.split("::")[1].split(" ")[0]
            if not re.search(rf"def {re.escape(testname)}\b", path.read_text()):
                missing.append(f"{instrument} -> {testname} (not in {filename})")
    assert not missing, "\n".join(missing)


# ------------------------------------- the four calibrations added here

def test_spring_model_error_recovers_a_planted_mode():
    """Plant a mode from one spring model; that model must score ~1.

    Previously this instrument was only checked by asserting its spread on real
    structures was "modest" — which is a statement about the answer, not about
    the instrument. If it had been computing overlaps wrongly, a modest spread
    would have looked exactly the same.
    """
    from piezo1.analysis.model_error import spring_model_error
    from piezo1.physics.anm import ANM

    rng = np.random.default_rng(3)
    blocks = [rng.normal(scale=9.0, size=(60, 3)) for _ in range(3)]
    anm = ANM.from_trimer(blocks, cutoff=15.0, spring="inverse_square").build()
    planted = anm.calc_modes(n_modes=12).vectors[0]

    error = spring_model_error(blocks, planted, n_modes=12, cutoff=15.0)
    recovered = error.values["inverse_square"]
    others = [v for k, v in error.values.items() if k != "inverse_square"]

    # The calibration is the SEPARATION, not an absolute threshold. On these
    # random coordinates the low modes are near-degenerate, and ARPACK starts
    # from an unseeded random vector, so the recovered overlap varies run to
    # run between about 0.95 and 1.00 with identical inputs. That is a property
    # of degenerate geometry, not of the instrument: the real structures have
    # well-separated low modes and `anm.gating_overlap` reproduces to every
    # digit. Asserting a tight absolute bound here would make this test flaky
    # for a reason that has nothing to do with what it is checking.
    assert recovered > 0.9, (
        f"the model that generated the displacement must recover it, got "
        f"{recovered:.4f}")
    assert recovered > 1.5 * max(others), (
        f"the generating model must stand clear of the others: "
        f"{recovered:.4f} vs {max(others):.4f}")
    assert error.spread > 0.3, "the instrument must be able to say they differ"


def test_minimum_detectable_effect_delivers_its_target_power():
    """The effect it returns must actually be detectable at the stated power.

    Until now this was only checked against the project's own recorded results,
    which is circular: the same function supplies the power statements those
    rounds are judged by. This closes the loop against the simulator instead.
    """
    from piezo1.analysis.design import minimum_detectable_effect, power_curve

    for n_a, n_b in ((20, 20), (16, 9)):
        delta = minimum_detectable_effect(n_a, n_b, target_power=0.8)
        achieved = power_curve(n_a, n_b, deltas=[-delta], n_simulations=3000,
                               n_permutations=999, seed=7).power[0]
        assert achieved == pytest.approx(0.80, abs=0.06), (
            f"MDE {delta:.3f} at n={n_a}/{n_b} gives power {achieved:.3f}")


def test_sensitivity_range_is_exactly_the_known_extremes():
    """Over known outputs the range must be their min and max, exactly.

    The existing tests check that it refuses to call itself a confidence
    interval — which is about the wording, not the arithmetic.
    """
    from piezo1.analysis.uncertainty import sensitivity

    values = [2.0, 5.0, 3.0]
    result = sensitivity(lambda i: values[i], settings=[0, 1, 2],
                         reference=1, knob="i", what="calibration")
    assert result.low == pytest.approx(2.0)
    assert result.high == pytest.approx(5.0)
    assert result.estimate == pytest.approx(5.0), "estimate is the reference"

    # A single setting has no spread — the degenerate case must not widen it.
    flat = sensitivity(lambda i: 4.0, settings=[0], reference=0)
    assert flat.low == flat.high == pytest.approx(4.0)


def test_permutation_p_matches_exhaustive_enumeration():
    """For samples small enough to enumerate, the sampler must agree.

    Eight values give C(8,4) = 70 partitions, so the exact one-sided p is
    computable directly. The sampled value is very slightly larger because of
    the (r+1)/(n+1) convention, which is deliberately conservative — that
    direction is asserted rather than merely tolerated.
    """
    from piezo1.analysis.validation import permutation_test

    a = np.array([1.0, 2, 3, 4])
    b = np.array([5.0, 6, 7, 8])
    pool = np.concatenate([a, b])
    observed = a.mean() - b.mean()

    hits = total = 0
    for index in itertools.combinations(range(8), 4):
        mask = np.zeros(8, bool)
        mask[list(index)] = True
        if pool[mask].mean() - pool[~mask].mean() <= observed:
            hits += 1
        total += 1
    exact = hits / total

    sampled = permutation_test(a, b, n_permutations=20000, seed=1).p_value
    assert sampled == pytest.approx(exact, abs=0.005)
    assert sampled >= exact, "the (r+1)/(n+1) convention must not anti-conserve"


def test_the_calibrations_added_here_would_fail_on_a_broken_instrument():
    """Each new calibration must be able to reject, or it asserts nothing."""
    from piezo1.analysis.uncertainty import sensitivity

    # A 'sensitivity' that returned the reference for both ends would pass any
    # wording test and fail this one.
    values = [2.0, 5.0, 3.0]
    result = sensitivity(lambda i: values[i], settings=[0, 1, 2], reference=1)
    assert result.low != result.high

    from piezo1.analysis.validation import permutation_test
    identical = np.arange(4, dtype=float)
    assert permutation_test(identical, identical.copy()).p_value > 0.3, (
        "identical samples must not look significant")
